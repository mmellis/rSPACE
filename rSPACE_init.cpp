#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void make_grid(double *, double *, double *, int *, int *);
extern void sample_ind(double *, double *, int *, double *, int *, int *, int *);
extern void use_surface(double *, double *, int *, double *, double *, double *, int *, double *, double *, double *);
extern void calc_prob(double *, int *, int *, double *, int *, int *, double *);


static R_NativePrimitiveArgType make_grid_t[] = {
    REALSXP, REALSXP, REALSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType sample_ind_t[] = {
    REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType use_surface_t[] = {
    REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType calc_prob_t[] = {
    REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP
};

static const R_CMethodDef CEntries[] = {
	{"calc_prob",   (DL_FUNC) &calc_prob,    7, calc_prob_t},
    {"make_grid",   (DL_FUNC) &make_grid,    5, make_grid_t},
    {"sample_ind",  (DL_FUNC) &sample_ind,   7, sample_ind_t},
    {"use_surface", (DL_FUNC) &use_surface, 10, use_surface_t},
    {NULL, NULL, 0}
};

void R_init_rSPACE(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}
