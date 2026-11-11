//! EM fit of the MAGICMAP mixture, ported from the R code in
//! R/magicmap.R (magicmap_1k, estimator "MLE") and R/magicmap_fitting_funcs.R.
//!
//! This is a like-for-like port: the same updates, the same convergence rules
//! and tolerances, and the same SQUAREM acceleration scheme (squarem1 of the
//! SQUAREM R package with the settings magicmap_1k uses), with arithmetic
//! written in the same order as the R expressions. The R code remains the
//! reference; tests/testthat/test-rust-em.R compares the two.
//!
//! Each starting point is fitted independently, so starts can be spread over
//! threads without changing any result.
//!
//! Within one start, every pass over the variants is done in blocks of BLOCK
//! variants: each block is summed on its own, in variant order, and the block
//! sums are added in block order. Threads left over after one per start share
//! the blocks. The blocks do not depend on the number of threads, so neither
//! does any result. With BLOCK variants or fewer there is one block, and the
//! arithmetic is in the same order as the R expressions.

// Index loops and negated float comparisons are kept as they are: they mirror
// the R expressions, including how those treat NaN.
#![allow(clippy::needless_range_loop, clippy::neg_cmp_op_on_partial_ord)]

use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;
use std::thread;

const LN_2PI: f64 = 1.837_877_066_409_345_3;

/// Variants per block (see the module comment).
const BLOCK: usize = 8192;

/// Pointer to an output array that worker threads fill at disjoint positions.
struct SharedOut(*mut f64);
// Safety: every use writes each position from exactly one block.
unsafe impl Send for SharedOut {}
unsafe impl Sync for SharedOut {}
impl SharedOut {
    /// # Safety
    /// `index` must be in bounds and written by one thread only.
    #[inline]
    unsafe fn write(&self, index: usize, value: f64) {
        *self.0.add(index) = value;
    }
}

/// A step of the fit that R would have stopped on (an error or a NaN).
struct FitFailed;

struct Problem<'a> {
    n: usize,
    k: usize,
    x: &'a [f64],
    y: &'a [f64],
    inner_tol: f64,
    max_iter: usize,
    /// threads that share the blocks of one start
    threads: usize,
    // per-variant quantities that do not change during the fit
    sx2: Vec<f64>,
    sy2: Vec<f64>,
    /// re*sx*sy
    c12: Vec<f64>,
    /// sx^2*sy^2 - c12^2
    det: Vec<f64>,
    /// -log(2*pi) - 0.5*log(det)
    dens_const: Vec<f64>,
    /// sx^2*sy^2*(1 - re^2), as written in slope_est_func
    slope_w: Vec<f64>,
}

impl<'a> Problem<'a> {
    #[allow(clippy::too_many_arguments)]
    fn new(
        x: &'a [f64],
        sx: &'a [f64],
        y: &'a [f64],
        sy: &'a [f64],
        re: f64,
        k: usize,
        inner_tol: f64,
        max_iter: usize,
        threads: usize,
    ) -> Self {
        let n = x.len();
        let mut sx2 = Vec::with_capacity(n);
        let mut sy2 = Vec::with_capacity(n);
        let mut c12 = Vec::with_capacity(n);
        let mut det = Vec::with_capacity(n);
        let mut dens_const = Vec::with_capacity(n);
        let mut slope_w = Vec::with_capacity(n);
        for i in 0..n {
            let a = sx[i] * sx[i];
            let b = sy[i] * sy[i];
            let c = re * sx[i] * sy[i];
            let d = a * b - c * c;
            sx2.push(a);
            sy2.push(b);
            c12.push(c);
            det.push(d);
            dens_const.push(-LN_2PI - 0.5 * d.ln());
            slope_w.push(a * b * (1.0 - re * re));
        }
        Problem {
            n,
            k,
            x,
            y,
            inner_tol,
            max_iter,
            threads,
            sx2,
            sy2,
            c12,
            det,
            dens_const,
            slope_w,
        }
    }

    /// component_logdens(): log density of variant i under a line with slope b
    #[inline]
    fn logdens(&self, i: usize, b: f64, xpred: f64) -> f64 {
        let r1 = self.x[i] - xpred;
        let r2 = self.y[i] - b * xpred;
        self.dens_const[i]
            - 0.5 * (self.sy2[i] * r1 * r1 - 2.0 * self.c12[i] * r1 * r2 + self.sx2[i] * r2 * r2)
                / self.det[i]
    }

    /// Runs `f(from, to)` on every block of variants and returns the results
    /// in block order, sharing the blocks among this start's threads.
    fn blocks<R: Send>(&self, f: impl Fn(usize, usize) -> R + Sync) -> Vec<R> {
        let n = self.n;
        let n_blocks = n.div_ceil(BLOCK).max(1);
        let run = |block: usize| f((block * BLOCK).min(n), ((block + 1) * BLOCK).min(n));
        let workers = self.threads.min(n_blocks);
        if workers <= 1 {
            return (0..n_blocks).map(run).collect();
        }
        let per_worker = n_blocks.div_ceil(workers);
        let run = &run;
        thread::scope(|scope| {
            let handles: Vec<_> = (0..workers)
                .map(|w| {
                    scope.spawn(move || {
                        let from = (w * per_worker).min(n_blocks);
                        let to = ((w + 1) * per_worker).min(n_blocks);
                        (from..to).map(run).collect::<Vec<R>>()
                    })
                })
                .collect();
            let mut out = Vec::with_capacity(n_blocks);
            for handle in handles {
                out.extend(handle.join().expect("EM block worker panicked"));
            }
            out
        })
    }

    /// ll_comp_func(): marginal log-likelihood
    fn loglik(&self, b: &[f64], p: &[f64], xpred: &[f64]) -> f64 {
        let mut total = 0.0;
        for part in self.blocks(|from, to| self.loglik_block(b, p, xpred, from, to)) {
            total += part;
        }
        total
    }

    fn loglik_block(&self, b: &[f64], p: &[f64], xpred: &[f64], from: usize, to: usize) -> f64 {
        let k = self.k;
        let mut work = vec![0.0; k];
        let mut total = 0.0;
        for i in from..to {
            let mut mx = f64::NEG_INFINITY;
            let mut any_nan = false;
            for c in 0..k {
                let v = self.logdens(i, b[c], xpred[i]) + p[c].ln();
                work[c] = v;
                if v.is_nan() {
                    any_nan = true;
                } else if v > mx {
                    mx = v;
                }
            }
            if any_nan {
                mx = f64::NAN;
            }
            if mx.is_finite() {
                let mut s = 0.0;
                for value in work.iter().take(k) {
                    s += (value - mx).exp();
                }
                total += mx + s.ln();
            } else {
                total += mx;
            }
        }
        total
    }

    /// estep_posteriors(): post is n x k, stored by component (column-major)
    fn estep(&self, b: &[f64], p: &[f64], xpred: &[f64], post: &mut [f64]) {
        assert_eq!(post.len(), self.n * self.k);
        let out = SharedOut(post.as_mut_ptr());
        self.blocks(|from, to| self.estep_block(b, p, xpred, &out, from, to));
    }

    fn estep_block(&self, b: &[f64], p: &[f64], xpred: &[f64], post: &SharedOut, from: usize, to: usize) {
        let (n, k) = (self.n, self.k);
        let mut work = vec![0.0; k];
        for i in from..to {
            let mut mx = f64::NEG_INFINITY;
            let mut any_nan = false;
            for c in 0..k {
                let v = self.logdens(i, b[c], xpred[i]) + p[c].ln();
                work[c] = v;
                if v.is_nan() {
                    any_nan = true;
                } else if v > mx {
                    mx = v;
                }
            }
            if any_nan {
                mx = f64::NAN;
            }
            let mut s = 0.0;
            for c in 0..k {
                work[c] -= mx;
                s += work[c].exp();
            }
            let log_s = s.ln();
            let mut row = 0.0;
            for c in 0..k {
                let v = (work[c] - log_s).exp();
                work[c] = v;
                row += v;
            }
            for c in 0..k {
                // Safety: c * n + i < n * k, and variant i belongs to this block only
                unsafe { post.write(c * n + i, work[c] / row) };
            }
        }
    }

    /// colMeans(post): the mixing proportions
    fn proportions(&self, post: &[f64], p: &mut [f64]) {
        let (n, k) = (self.n, self.k);
        let parts = self.blocks(|from, to| {
            let mut sums = vec![0.0; k];
            for c in 0..k {
                let mut s = 0.0;
                for i in from..to {
                    s += post[c * n + i];
                }
                sums[c] = s;
            }
            sums
        });
        for c in 0..k {
            let mut s = 0.0;
            for part in &parts {
                s += part[c];
            }
            p[c] = s / n as f64;
        }
    }

    /// slope_est_func()
    fn slopes(&self, xpred: &[f64], post: &[f64], b: &mut [f64]) {
        let (n, k) = (self.n, self.k);
        // per block: numerators of the k slopes, then denominators
        let parts = self.blocks(|from, to| {
            let mut sums = vec![0.0; 2 * k];
            for c in 0..k {
                let mut numer = 0.0;
                let mut denom = 0.0;
                for i in from..to {
                    let w = self.slope_w[i];
                    let pc = post[c * n + i];
                    let xp = xpred[i];
                    numer += (pc * xp * (self.sx2[i] * self.y[i] - self.c12[i] * (self.x[i] - xp))) / w;
                    denom += (pc * self.sx2[i] * (xp * xp)) / w;
                }
                sums[c] = numer;
                sums[k + c] = denom;
            }
            sums
        });
        for (c, slope) in b.iter_mut().enumerate().take(k) {
            let mut numer = 0.0;
            let mut denom = 0.0;
            for part in &parts {
                numer += part[c];
                denom += part[k + c];
            }
            *slope = numer / denom;
        }
    }

    /// xpred_est_func(): update of the latent target effects
    fn xpred(&self, b: &[f64], post: &[f64], out: &mut [f64]) {
        let (n, k) = (self.n, self.k);
        assert_eq!(out.len(), n);
        let out = SharedOut(out.as_mut_ptr());
        let out = &out;
        self.blocks(|from, to| {
            for i in from..to {
                let mut numer = 0.0;
                let mut denom = 0.0;
                for c in 0..k {
                    let pc = post[c * n + i];
                    let bc = b[c];
                    numer += pc
                        * (self.x[i] * (self.sy2[i] - bc * self.c12[i])
                            + self.y[i] * (bc * self.sx2[i] - self.c12[i]));
                    denom += pc * (self.sy2[i] + bc * bc * self.sx2[i] - 2.0 * bc * self.c12[i]);
                }
                // Safety: i < n, and variant i belongs to this block only
                unsafe { out.write(i, numer / denom) };
            }
        });
    }

    /// sum(Q_loglik_matrix)
    fn q_loglik(&self, b: &[f64], p: &[f64], xpred: &[f64], post: &[f64]) -> f64 {
        let (n, k) = (self.n, self.k);
        let parts = self.blocks(|from, to| {
            let mut q = 0.0;
            for c in 0..k {
                let log_p = p[c].ln();
                for i in from..to {
                    q += post[c * n + i] * (self.logdens(i, b[c], xpred[i]) + log_p);
                }
            }
            q
        });
        let mut total = 0.0;
        for part in parts {
            total += part;
        }
        total
    }

    /// The sums over variants that magicmap_vcov() (R/magicmap_variance.R)
    /// needs. Each variant's score and Hessian are formed as in the R code,
    /// with the same order of arithmetic; only the sums over variants are
    /// added in a different order (by block).
    ///
    /// With na = 2k - 1 parameters (k slopes, k - 1 proportions), fills three
    /// na x na matrices: `h_aa` = sum of the per-variant Hessian blocks,
    /// `cross` = crossprod(h_ax, h_ax / h_xx) and `b_mat` = crossprod of the
    /// profiled scores. Returns the number of variants with h_xx >= 0 and the
    /// number with h_xx NaN.
    fn info_sums(
        &self,
        beta: &[f64],
        priorp: &[f64],
        xpred: &[f64],
        h_aa: &mut [f64],
        cross: &mut [f64],
        b_mat: &mut [f64],
    ) -> (usize, usize) {
        let k = self.k;
        let na = 2 * k - 1;
        let npar = na + 1;
        let parts = self.blocks(|from, to| {
            // h_aa, cross, b_mat one after another
            let mut sums = vec![0.0; 3 * na * na];
            let mut nonneg = 0usize;
            let mut nan = 0usize;
            let mut la = vec![0.0; k];
            let mut post = vec![0.0; k];
            let mut u = vec![0.0; k];
            let mut v = vec![0.0; k];
            let mut uu = vec![0.0; k];
            let mut ux = vec![0.0; k];
            let mut vv = vec![0.0; k];
            let mut g = vec![0.0; npar];
            let mut s = vec![0.0; npar];
            let mut h = vec![0.0; npar * npar];
            let mut h_ax = vec![0.0; na];
            let mut s_prof = vec![0.0; na];
            for i in from..to {
                let s1sq = self.sx2[i];
                let s2sq = self.sy2[i];
                let c12 = self.c12[i];
                let d = self.det[i];
                let xp = xpred[i];
                let r1 = self.x[i] - xp;
                for c in 0..k {
                    let b = beta[c];
                    let r2 = self.y[i] - b * xp;
                    let q = (s2sq * (r1 * r1) - 2.0 * c12 * r1 * r2 + s1sq * (r2 * r2)) / d;
                    la[c] = priorp[c].ln() - LN_2PI - 0.5 * d.ln() - 0.5 * q;
                    u[c] = -(xp / d) * (c12 * r1 - s1sq * r2);
                    v[c] = (s2sq * r1 - c12 * (r2 + b * r1) + b * s1sq * r2) / d;
                    uu[c] = -s1sq * (xp * xp) / d;
                    ux[c] = -(c12 * r1 - s1sq * r2 - xp * c12 + b * xp * s1sq) / d;
                    vv[c] = -(s2sq - 2.0 * b * c12 + (b * b) * s1sq) / d;
                }

                // posteriors
                let mut mx = f64::NEG_INFINITY;
                let mut any_nan = false;
                for c in 0..k {
                    if la[c].is_nan() {
                        any_nan = true;
                    } else if la[c] > mx {
                        mx = la[c];
                    }
                }
                if any_nan {
                    mx = f64::NAN;
                }
                let mut total = 0.0;
                for c in 0..k {
                    post[c] = (la[c] - mx).exp();
                    total += post[c];
                }
                for c in 0..k {
                    post[c] /= total;
                }

                // score s and Hessian h (upper triangle) of this variant's
                // log-likelihood in (slopes, proportions 1..k-1, its X)
                s.fill(0.0);
                h.fill(0.0);
                for c in 0..k {
                    let pc = post[c];
                    g.fill(0.0);
                    g[c] = u[c];
                    if k > 1 {
                        if c < k - 1 {
                            g[k + c] = 1.0 / priorp[c];
                        } else {
                            for j in k..na {
                                g[j] = -1.0 / priorp[k - 1];
                            }
                        }
                    }
                    g[npar - 1] = v[c];
                    for a in 0..npar {
                        s[a] += pc * g[a];
                        if g[a] == 0.0 {
                            continue;
                        }
                        for l in a..npar {
                            h[a * npar + l] += pc * g[a] * g[l];
                        }
                    }
                    h[c * npar + c] += pc * uu[c];
                    h[c * npar + npar - 1] += pc * ux[c];
                    h[(npar - 1) * npar + npar - 1] += pc * vv[c];
                    if k > 1 {
                        if c < k - 1 {
                            h[(k + c) * npar + k + c] -= pc / (priorp[c] * priorp[c]);
                        } else {
                            let last = pc / (priorp[k - 1] * priorp[k - 1]);
                            for a in k..na {
                                for l in a..na {
                                    h[a * npar + l] -= last;
                                }
                            }
                        }
                    }
                }
                for a in 0..npar {
                    for l in a..npar {
                        h[a * npar + l] -= s[a] * s[l];
                    }
                }

                // profile out this variant's X
                let h_xx = h[(npar - 1) * npar + npar - 1];
                if h_xx.is_nan() {
                    nan += 1;
                } else if h_xx >= 0.0 {
                    nonneg += 1;
                }
                let s_x = s[npar - 1] / h_xx;
                for a in 0..na {
                    h_ax[a] = h[a * npar + npar - 1];
                    s_prof[a] = s[a] - h_ax[a] * s_x;
                }
                for a in 0..na {
                    for l in a..na {
                        let at = a * na + l;
                        sums[at] += h[a * npar + l];
                        sums[na * na + at] += h_ax[a] * (h_ax[l] / h_xx);
                        sums[2 * na * na + at] += s_prof[a] * s_prof[l];
                    }
                }
            }
            (sums, nonneg, nan)
        });

        let mut sums = vec![0.0; 3 * na * na];
        let mut nonneg = 0;
        let mut nan = 0;
        for (part, part_nonneg, part_nan) in &parts {
            for (total, value) in sums.iter_mut().zip(part) {
                *total += value;
            }
            nonneg += part_nonneg;
            nan += part_nan;
        }
        // fill both triangles (column-major or row-major is then the same)
        for a in 0..na {
            for l in a..na {
                for (m, out) in [&mut *h_aa, &mut *cross, &mut *b_mat].into_iter().enumerate() {
                    let value = sums[m * na * na + a * na + l];
                    out[a * na + l] = value;
                    out[l * na + a] = value;
                }
            }
        }
        (nonneg, nan)
    }

    /// One EM cycle (em_fixed_point, estimator "MLE"): theta = (b, p, X).
    fn fixed_point(
        &self,
        theta: &[f64],
        out: &mut [f64],
        scratch: &mut Scratch,
    ) -> Result<(), FitFailed> {
        let k = self.k;
        let (b_curr, rest) = theta.split_at(k);
        let (p_curr, x_curr) = rest.split_at(k);
        let mut post = std::mem::take(&mut scratch.post);

        // E step
        self.estep(b_curr, p_curr, x_curr, &mut post);

        // M step
        let (b_new, rest) = out.split_at_mut(k);
        let (p_new, x_new) = rest.split_at_mut(k);
        self.proportions(&post, p_new);
        let any_zero = p_new.contains(&0.0);
        if any_zero {
            let mut s = 0.0;
            for value in p_new.iter_mut() {
                if !(*value >= 1e-10) && !value.is_nan() {
                    *value = 1e-10;
                }
                s += *value;
            }
            for value in p_new.iter_mut() {
                *value /= s;
            }
        }
        b_new.copy_from_slice(b_curr);
        x_new.copy_from_slice(x_curr);

        let mut q_old = f64::NEG_INFINITY;
        let mut result = Ok(());
        let mut y_iter = 0;
        while y_iter < self.max_iter {
            y_iter += 1;
            self.slopes(x_new, &post, b_new);
            let mut x_next = std::mem::take(&mut scratch.x_next);
            self.xpred(b_new, &post, &mut x_next);
            x_new.copy_from_slice(&x_next);
            scratch.x_next = x_next;

            let q_new = self.q_loglik(b_new, p_new, x_new, &post);

            let change = ((q_new - q_old) / q_new).abs();
            if change.is_nan() {
                result = Err(FitFailed);
                break;
            }
            if change < self.inner_tol {
                break;
            }
            q_old = q_new;
        }

        scratch.post = post;
        result?;
        if out.iter().any(|v| v.is_nan()) {
            return Err(FitFailed);
        }
        Ok(())
    }
}

struct Scratch {
    post: Vec<f64>,
    x_next: Vec<f64>,
}

impl Scratch {
    fn new(n: usize, k: usize) -> Self {
        Scratch {
            post: vec![0.0; n * k],
            x_next: vec![0.0; n],
        }
    }
}

struct StartFit {
    theta: Vec<f64>,
    loglik: f64,
    fpevals: usize,
    converged: bool,
}

/// SQUAREM::squarem1 as called by magicmap_1k: K = 1, method 3,
/// step.min0 = step.max0 = 1, mstep = 4, objfn.inc = 1, maximising.
fn squarem_fit(
    problem: &Problem,
    start: &[f64],
    em_tol: f64,
    max_iter: usize,
) -> Result<StartFit, FitFailed> {
    let (n, k) = (problem.n, problem.k);
    let len = 2 * k + n;
    let mut scratch = Scratch::new(n, k);
    let objective = |theta: &[f64]| problem.loglik(&theta[..k], &theta[k..2 * k], &theta[2 * k..]);

    let step_min = 1.0_f64;
    let step_max0 = 1.0_f64;
    let mut step_max = 1.0_f64;
    let mstep = 4.0_f64;
    let objfn_inc = 1.0_f64;

    let mut p = start.to_vec();
    let mut p1 = vec![0.0; len];
    let mut p2 = vec![0.0; len];
    let mut p_new = vec![0.0; len];
    let mut p_tmp = vec![0.0; len];
    let mut lold = objective(&p);
    let mut feval = 0usize;

    while feval < max_iter {
        problem.fixed_point(&p, &mut p1, &mut scratch)?;
        feval += 1;
        let mut sr2 = 0.0;
        for j in 0..len {
            let q1 = p1[j] - p[j];
            sr2 += q1 * q1;
        }
        if sr2.sqrt() < em_tol {
            p.copy_from_slice(&p1);
            break;
        }

        problem.fixed_point(&p1, &mut p2, &mut scratch)?;
        feval += 1;
        let mut sq2 = 0.0;
        for j in 0..len {
            let q2 = p2[j] - p1[j];
            sq2 += q2 * q2;
        }
        if sq2.sqrt() < em_tol {
            p.copy_from_slice(&p2);
            break;
        }

        let mut sv2 = 0.0;
        for j in 0..len {
            let q1 = p1[j] - p[j];
            let q2 = p2[j] - p1[j];
            sv2 += (q2 - q1) * (q2 - q1);
        }
        let mut alpha = (sr2 / sv2).sqrt();
        if alpha.is_nan() {
            // R: if (abs(alpha - 1) > 0.01) fails on NA
            return Err(FitFailed);
        }
        alpha = step_min.max(step_max.min(alpha));
        for j in 0..len {
            let q1 = p1[j] - p[j];
            let q2 = p2[j] - p1[j];
            p_new[j] = p[j] + 2.0 * alpha * q1 + alpha * alpha * (q2 - q1);
        }
        let mut extrapolation_failed = false;
        if (alpha - 1.0).abs() > 0.01 {
            let ok = problem.fixed_point(&p_new, &mut p_tmp, &mut scratch).is_ok();
            feval += 1;
            if ok {
                p_new.copy_from_slice(&p_tmp);
            } else {
                extrapolation_failed = true;
            }
        }

        let lnew;
        if extrapolation_failed {
            p_new.copy_from_slice(&p2);
            lnew = objective(&p2);
            if alpha == step_max {
                step_max = step_max0.max(step_max / mstep);
            }
            alpha = 1.0;
        } else {
            let candidate = objective(&p_new);
            if !candidate.is_nan() && lold.is_nan() {
                // R: the monotonicity comparison is NA and `if` fails
                return Err(FitFailed);
            }
            if candidate.is_nan() || candidate < lold - objfn_inc {
                p_new.copy_from_slice(&p2);
                lnew = objective(&p2);
                if alpha == step_max {
                    step_max = step_max0.max(step_max / mstep);
                }
                alpha = 1.0;
            } else {
                lnew = candidate;
            }
        }
        if alpha == step_max {
            step_max *= mstep;
        }
        p.copy_from_slice(&p_new);
        if !lnew.is_nan() {
            lold = lnew;
        }
    }

    let converged = feval < max_iter;
    // magicmap_1k recalculates the log-likelihood at the returned parameters
    let loglik = objective(&p);
    Ok(StartFit {
        theta: p,
        loglik,
        fpevals: feval,
        converged,
    })
}

fn worker_count(n_starts: usize, requested_threads: usize) -> usize {
    let available = thread::available_parallelism()
        .map(|value| value.get())
        .unwrap_or(1);
    requested_threads.max(1).min(available).min(n_starts.max(1))
}

/// Fit from every start. `starts` is an n_starts x 2k column-major matrix
/// (slopes then proportions). Outputs are column-major with one row per start.
/// Returns 0, or 1 if the routine panicked.
///
/// # Safety
/// All pointers must be valid for the lengths implied by n, k and n_starts.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn magicmap_em_starts(
    x: *const f64,
    sx: *const f64,
    y: *const f64,
    sy: *const f64,
    n: usize,
    re: f64,
    k: usize,
    starts: *const f64,
    n_starts: usize,
    em_tol: f64,
    inner_tol: f64,
    max_iter: usize,
    requested_threads: usize,
    b_out: *mut f64,
    p_out: *mut f64,
    x_out: *mut f64,
    ll_out: *mut f64,
    iter_out: *mut i32,
    converged_out: *mut i32,
    status_out: *mut i32,
) -> i32 {
    let result = catch_unwind(AssertUnwindSafe(|| {
        let x = slice::from_raw_parts(x, n);
        let sx = slice::from_raw_parts(sx, n);
        let y = slice::from_raw_parts(y, n);
        let sy = slice::from_raw_parts(sy, n);
        let starts = slice::from_raw_parts(starts, n_starts * 2 * k);
        // one thread per start first; any left over share the blocks within a start
        let workers = worker_count(n_starts, requested_threads);
        let available = thread::available_parallelism()
            .map(|value| value.get())
            .unwrap_or(1);
        let block_threads = (requested_threads.max(1).min(available) / workers.max(1)).max(1);
        let problem = Problem::new(x, sx, y, sy, re, k, inner_tol, max_iter, block_threads);

        let fit_one = |s: usize| -> Option<StartFit> {
            let mut theta = Vec::with_capacity(2 * k + n);
            for j in 0..2 * k {
                theta.push(starts[j * n_starts + s]);
            }
            theta.extend_from_slice(x); // xpred_init <- dat$x
            squarem_fit(&problem, &theta, em_tol, max_iter).ok()
        };

        let mut fits: Vec<Option<StartFit>> = Vec::with_capacity(n_starts);
        if workers <= 1 {
            for s in 0..n_starts {
                fits.push(fit_one(s));
            }
        } else {
            let chunk = n_starts.div_ceil(workers);
            let fit_one = &fit_one;
            let pieces: Vec<Vec<Option<StartFit>>> = thread::scope(|scope| {
                let handles: Vec<_> = (0..workers)
                    .map(|w| {
                        scope.spawn(move || {
                            let from = (w * chunk).min(n_starts);
                            let to = ((w + 1) * chunk).min(n_starts);
                            (from..to).map(fit_one).collect::<Vec<_>>()
                        })
                    })
                    .collect();
                handles
                    .into_iter()
                    .map(|handle| handle.join().expect("EM worker panicked"))
                    .collect()
            });
            for piece in pieces {
                fits.extend(piece);
            }
        }

        let b_out = slice::from_raw_parts_mut(b_out, n_starts * k);
        let p_out = slice::from_raw_parts_mut(p_out, n_starts * k);
        let x_out = slice::from_raw_parts_mut(x_out, n_starts * n);
        let ll_out = slice::from_raw_parts_mut(ll_out, n_starts);
        let iter_out = slice::from_raw_parts_mut(iter_out, n_starts);
        let converged_out = slice::from_raw_parts_mut(converged_out, n_starts);
        let status_out = slice::from_raw_parts_mut(status_out, n_starts);
        for (s, fit) in fits.iter().enumerate() {
            match fit {
                Some(fit) => {
                    for c in 0..k {
                        b_out[c * n_starts + s] = fit.theta[c];
                        p_out[c * n_starts + s] = fit.theta[k + c];
                    }
                    for i in 0..n {
                        x_out[i * n_starts + s] = fit.theta[2 * k + i];
                    }
                    ll_out[s] = fit.loglik;
                    iter_out[s] = fit.fpevals as i32;
                    converged_out[s] = i32::from(fit.converged);
                    status_out[s] = 0;
                }
                None => {
                    for c in 0..k {
                        b_out[c * n_starts + s] = f64::NAN;
                        p_out[c * n_starts + s] = f64::NAN;
                    }
                    for i in 0..n {
                        x_out[i * n_starts + s] = f64::NAN;
                    }
                    ll_out[s] = f64::NAN;
                    iter_out[s] = 0;
                    converged_out[s] = 0;
                    status_out[s] = 1;
                }
            }
        }
    }));
    match result {
        Ok(()) => 0,
        Err(_) => 1,
    }
}

/// Sums over variants for the standard errors (see Problem::info_sums).
/// `h_aa_out`, `cross_out` and `b_out` are (2k-1) x (2k-1) matrices;
/// `counts_out` receives the number of variants with h_xx >= 0 and with h_xx
/// NaN. Returns 0, or 1 if the routine panicked.
///
/// # Safety
/// All pointers must be valid for the lengths implied by n and k.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn magicmap_info_sums(
    x: *const f64,
    sx: *const f64,
    y: *const f64,
    sy: *const f64,
    n: usize,
    re: f64,
    k: usize,
    beta: *const f64,
    priorp: *const f64,
    xpred: *const f64,
    requested_threads: usize,
    h_aa_out: *mut f64,
    cross_out: *mut f64,
    b_out: *mut f64,
    counts_out: *mut i32,
) -> i32 {
    let result = catch_unwind(AssertUnwindSafe(|| {
        let x = slice::from_raw_parts(x, n);
        let sx = slice::from_raw_parts(sx, n);
        let y = slice::from_raw_parts(y, n);
        let sy = slice::from_raw_parts(sy, n);
        let beta = slice::from_raw_parts(beta, k);
        let priorp = slice::from_raw_parts(priorp, k);
        let xpred = slice::from_raw_parts(xpred, n);
        let na = 2 * k - 1;
        let h_aa = slice::from_raw_parts_mut(h_aa_out, na * na);
        let cross = slice::from_raw_parts_mut(cross_out, na * na);
        let b_mat = slice::from_raw_parts_mut(b_out, na * na);
        let counts = slice::from_raw_parts_mut(counts_out, 2);
        let threads = worker_count(usize::MAX, requested_threads);
        // the tolerances and iteration limit are not used by info_sums
        let problem = Problem::new(x, sx, y, sy, re, k, 0.0, 0, threads);
        let (nonneg, nan) = problem.info_sums(beta, priorp, xpred, h_aa, cross, b_mat);
        counts[0] = nonneg.min(i32::MAX as usize) as i32;
        counts[1] = nan.min(i32::MAX as usize) as i32;
    }));
    match result {
        Ok(()) => 0,
        Err(_) => 1,
    }
}
