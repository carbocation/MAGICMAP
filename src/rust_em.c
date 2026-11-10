#include <R.h>
#include <Rinternals.h>
#include <stddef.h>

extern int magicmap_em_starts(
    const double *x,
    const double *sx,
    const double *y,
    const double *sy,
    size_t n,
    double re,
    size_t k,
    const double *starts,
    size_t n_starts,
    double em_tol,
    double inner_tol,
    size_t max_iter,
    size_t requested_threads,
    double *b_out,
    double *p_out,
    double *x_out,
    double *ll_out,
    int *iter_out,
    int *converged_out,
    int *status_out);

static void check_real_vector(SEXP value, R_xlen_t length, const char *name) {
    if (!isReal(value) || XLENGTH(value) != length) {
        error("'%s' must be a double vector of length %lld", name, (long long) length);
    }
}

/* Fit the mixture by EM from each row of `starts` (slopes then proportions).
   Returns list(b, p, x, logLik, iter, converged, status); one row / element
   per start. status 0 = fitted, 1 = the fit failed (as an R error would). */
SEXP C_magicmap_em_starts(
    SEXP x, SEXP sx, SEXP y, SEXP sy, SEXP re, SEXP k, SEXP starts,
    SEXP em_tol, SEXP inner_tol, SEXP max_iter, SEXP threads) {

    if (!isReal(x) || XLENGTH(x) < 1) {
        error("'x' must be a non-empty double vector");
    }
    R_xlen_t n = XLENGTH(x);
    check_real_vector(sx, n, "sx");
    check_real_vector(y, n, "y");
    check_real_vector(sy, n, "sy");

    int components = asInteger(k);
    if (components == NA_INTEGER || components < 1) {
        error("'k' must be a positive whole number");
    }
    if (!isReal(starts) || !isMatrix(starts)) {
        error("'starts' must be a double matrix");
    }
    SEXP dimensions = getAttrib(starts, R_DimSymbol);
    int n_starts = INTEGER(dimensions)[0];
    if (n_starts < 1 || INTEGER(dimensions)[1] != 2 * components) {
        error("'starts' must have at least one row and 2*k columns");
    }

    double re_value = asReal(re);
    double em_tol_value = asReal(em_tol);
    double inner_tol_value = asReal(inner_tol);
    int max_iter_value = asInteger(max_iter);
    int threads_value = asInteger(threads);
    if (!R_FINITE(re_value) || !R_FINITE(em_tol_value) || !R_FINITE(inner_tol_value)) {
        error("'re', 'em_tol' and 'inner_tol' must be finite");
    }
    if (max_iter_value == NA_INTEGER || max_iter_value < 1) {
        error("'max_iter' must be a positive whole number");
    }
    if (threads_value == NA_INTEGER || threads_value < 1) {
        error("'threads' must be a positive whole number");
    }

    SEXP b_out = PROTECT(allocMatrix(REALSXP, n_starts, components));
    SEXP p_out = PROTECT(allocMatrix(REALSXP, n_starts, components));
    SEXP x_out = PROTECT(allocMatrix(REALSXP, n_starts, (int) n));
    SEXP ll_out = PROTECT(allocVector(REALSXP, n_starts));
    SEXP iter_out = PROTECT(allocVector(INTSXP, n_starts));
    SEXP converged_out = PROTECT(allocVector(LGLSXP, n_starts));
    SEXP status_out = PROTECT(allocVector(INTSXP, n_starts));

    int failed = magicmap_em_starts(
        REAL(x), REAL(sx), REAL(y), REAL(sy), (size_t) n, re_value,
        (size_t) components, REAL(starts), (size_t) n_starts,
        em_tol_value, inner_tol_value, (size_t) max_iter_value,
        (size_t) threads_value,
        REAL(b_out), REAL(p_out), REAL(x_out), REAL(ll_out),
        INTEGER(iter_out), LOGICAL(converged_out), INTEGER(status_out));
    if (failed != 0) {
        UNPROTECT(7);
        error("the Rust EM routine failed unexpectedly");
    }

    const char *names[] = {"b", "p", "x", "logLik", "iter", "converged", "status", ""};
    SEXP output = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(output, 0, b_out);
    SET_VECTOR_ELT(output, 1, p_out);
    SET_VECTOR_ELT(output, 2, x_out);
    SET_VECTOR_ELT(output, 3, ll_out);
    SET_VECTOR_ELT(output, 4, iter_out);
    SET_VECTOR_ELT(output, 5, converged_out);
    SET_VECTOR_ELT(output, 6, status_out);
    UNPROTECT(8);
    return output;
}
