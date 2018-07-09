// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// make_grid
void make_grid(double x[], double y[], double *grid_size, int *pixels, int grid[]);
RcppExport SEXP _rSPACE_make_grid(SEXP x[]SEXP, SEXP y[]SEXP, SEXP *grid_sizeSEXP, SEXP *pixelsSEXP, SEXP grid[]SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x[](x[]SEXP);
    Rcpp::traits::input_parameter< double >::type y[](y[]SEXP);
    Rcpp::traits::input_parameter< double >::type *grid_size(*grid_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type *pixels(*pixelsSEXP);
    Rcpp::traits::input_parameter< int >::type grid[](grid[]SEXP);
    make_grid(x[], y[], *grid_size, *pixels, grid[]);
    return R_NilValue;
END_RCPP
}
// sample_ind
void sample_ind(double x[], double y[], int *N, double *buffer, int use[], int *maxid, int new_order[]);
RcppExport SEXP _rSPACE_sample_ind(SEXP x[]SEXP, SEXP y[]SEXP, SEXP *NSEXP, SEXP *bufferSEXP, SEXP use[]SEXP, SEXP *maxidSEXP, SEXP new_order[]SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x[](x[]SEXP);
    Rcpp::traits::input_parameter< double >::type y[](y[]SEXP);
    Rcpp::traits::input_parameter< int >::type *N(*NSEXP);
    Rcpp::traits::input_parameter< double >::type *buffer(*bufferSEXP);
    Rcpp::traits::input_parameter< int >::type use[](use[]SEXP);
    Rcpp::traits::input_parameter< int >::type *maxid(*maxidSEXP);
    Rcpp::traits::input_parameter< int >::type new_order[](new_order[]SEXP);
    sample_ind(x[], y[], *N, *buffer, use[], *maxid, new_order[]);
    return R_NilValue;
END_RCPP
}
// use_surface
void use_surface(double x_wolv[], double y_wolv[], int *N_wolv, double x[], double y[], double snow[], int *pixels, double *sd_x, double *sd_y, double *trunc_cutoff);
RcppExport SEXP _rSPACE_use_surface(SEXP x_wolv[]SEXP, SEXP y_wolv[]SEXP, SEXP *N_wolvSEXP, SEXP x[]SEXP, SEXP y[]SEXP, SEXP snow[]SEXP, SEXP *pixelsSEXP, SEXP *sd_xSEXP, SEXP *sd_ySEXP, SEXP *trunc_cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x_wolv[](x_wolv[]SEXP);
    Rcpp::traits::input_parameter< double >::type y_wolv[](y_wolv[]SEXP);
    Rcpp::traits::input_parameter< int >::type *N_wolv(*N_wolvSEXP);
    Rcpp::traits::input_parameter< double >::type x[](x[]SEXP);
    Rcpp::traits::input_parameter< double >::type y[](y[]SEXP);
    Rcpp::traits::input_parameter< double >::type snow[](snow[]SEXP);
    Rcpp::traits::input_parameter< int >::type *pixels(*pixelsSEXP);
    Rcpp::traits::input_parameter< double >::type *sd_x(*sd_xSEXP);
    Rcpp::traits::input_parameter< double >::type *sd_y(*sd_ySEXP);
    Rcpp::traits::input_parameter< double >::type *trunc_cutoff(*trunc_cutoffSEXP);
    use_surface(x_wolv[], y_wolv[], *N_wolv, x[], y[], snow[], *pixels, *sd_x, *sd_y, *trunc_cutoff);
    return R_NilValue;
END_RCPP
}
// calc_prob
void calc_prob(double use[], int grid[], int detection[], double *detectionP, int *pixels, int *max_grid, double test[]);
RcppExport SEXP _rSPACE_calc_prob(SEXP use[]SEXP, SEXP grid[]SEXP, SEXP detection[]SEXP, SEXP *detectionPSEXP, SEXP *pixelsSEXP, SEXP *max_gridSEXP, SEXP test[]SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type use[](use[]SEXP);
    Rcpp::traits::input_parameter< int >::type grid[](grid[]SEXP);
    Rcpp::traits::input_parameter< int >::type detection[](detection[]SEXP);
    Rcpp::traits::input_parameter< double >::type *detectionP(*detectionPSEXP);
    Rcpp::traits::input_parameter< int >::type *pixels(*pixelsSEXP);
    Rcpp::traits::input_parameter< int >::type *max_grid(*max_gridSEXP);
    Rcpp::traits::input_parameter< double >::type test[](test[]SEXP);
    calc_prob(use[], grid[], detection[], *detectionP, *pixels, *max_grid, test[]);
    return R_NilValue;
END_RCPP
}

RcppExport void calc_prob(void *, void *, void *, void *, void *, void *, void *, void *);
RcppExport void make_grid(void *, void *, void *, void *, void *, void *);
RcppExport void sample_ind(void *, void *, void *, void *, void *, void *, void *, void *);
RcppExport void use_surface(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"_rSPACE_make_grid", (DL_FUNC) &_rSPACE_make_grid, 5},
    {"_rSPACE_sample_ind", (DL_FUNC) &_rSPACE_sample_ind, 7},
    {"_rSPACE_use_surface", (DL_FUNC) &_rSPACE_use_surface, 10},
    {"_rSPACE_calc_prob", (DL_FUNC) &_rSPACE_calc_prob, 7},
    {"calc_prob",   (DL_FUNC) &calc_prob,    8},
    {"make_grid",   (DL_FUNC) &make_grid,    6},
    {"sample_ind",  (DL_FUNC) &sample_ind,   8},
    {"use_surface", (DL_FUNC) &use_surface, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_rSPACE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}