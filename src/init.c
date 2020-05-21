#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void Romega( void *, void *, void *, void *, void *, void *,
		void *, void *, void *, void *, void *, void *);
extern void RmixQ(void *, void *, void *, void *, void *, void *,
		void *, void *, void *, void *, void *, void *);
extern void RcalcQ(void *, void *, void *, void *, void *, void *,
		void *, void *, void *);
extern void RcalcTau(void *, void *, void *, void *, void *, void *,
		void *, void * ); 
extern void RnegHW( void *, void *, void *, void *, void *, void *,
		void *, void *, void *, void *, void *);
extern void Rlogit( void *, void *, void *);
extern void Rexpit( void *, void *, void *);
extern void RtoNEF( void *, void *, void *, void *, void *, void *);
extern void RfromNEF( void *, void *, void *, void *, void *, void *);
extern void Rzhat( void *, void *, void *, void *, void *, void *,
		void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"Romega",         (DL_FUNC) &Romega,        12},
    {"RmixQ",          (DL_FUNC) &RmixQ,         12},
    {"RcalcQ",         (DL_FUNC) &RcalcQ,       9},
    {"RcalcTau",       (DL_FUNC) &RcalcTau,     8},
    {"RnegHW",       (DL_FUNC) &RnegHW,     11},
    {"Rlogit",       (DL_FUNC) &Rlogit,     3},
    {"Rexpit",       (DL_FUNC) &Rexpit,     3},
    {"RtoNEF",       (DL_FUNC) &RtoNEF,     6},
    {"RfromNEF",       (DL_FUNC) &RfromNEF,     6},
    {"Rzhat",       (DL_FUNC) &Rzhat,     9},
    {NULL, NULL, 0}
};

void R_init_maptpx(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
