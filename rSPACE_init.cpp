#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void calc_prob(void *, void *, void *, void *, void *, void *, void *);
extern void make_grid(void *, void *, void *, void *, void *);
extern void sample_ind(void *, void *, void *, void *, void *, void *, void *);
extern void use_surface(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"calc_prob",   (DL_FUNC) &calc_prob,    7},
    {"make_grid",   (DL_FUNC) &make_grid,    5},
    {"sample_ind",  (DL_FUNC) &sample_ind,   7},
    {"use_surface", (DL_FUNC) &use_surface, 10},
    {NULL, NULL, 0}
};

void R_init_rSPACE(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}