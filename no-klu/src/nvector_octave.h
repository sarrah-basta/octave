#include <stdio.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_config.h>
#include <dColVector.h>
#include <builtin-defun-decls.h>

#include "oct-locbuf.h"
#include "quit.h"

#include <oct.h>
#include <ov.h>
#include <ovl.h>



#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Macros NV_CONTENT_S, NV_DATA_S, NV_OWN_DATA_S,
 *        NV_LENGTH_S, and NV_Ith_S
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_C(v)  ( (ColumnVector *)(v->content) )
#define NV_STCONTENT_C(v) ( (ColumnVector)(v->))

#define NV_LENGTH_C(v)   ( NV_CONTENT_C(v)->numel() )

#define NV_DATA_C(v)     ( NV_CONTENT_C(v)->fortran_vec() )
#define NV_Ith_C(v,i)     ( NV_DATA_C(v)[i])
typedef struct _N_VectorContent_Octave *N_VectorContent_Octave;

SUNDIALS_EXPORT N_Vector N_VNew_Octave(int vec_length, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewEmpty_Octave (SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMake_Octave(ColumnVector v, sunindextype length, SUNContext sunctx);

SUNDIALS_EXPORT sunindextype N_VGetLength_Octave(N_Vector v);
SUNDIALS_EXPORT void N_VPrint_Octave(N_Vector v);
SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_Octave(N_Vector v);
SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Octave(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Octave(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Octave(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_Octave(N_Vector v, sunindextype *lrw, sunindextype *liw);
SUNDIALS_EXPORT double *N_VGetArrayPointer_Octave(N_Vector v);
SUNDIALS_EXPORT void N_VSetArrayPointer_Octave(realtype *v_data, N_Vector v);

// /* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_Octave(double a, N_Vector x, double b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Octave(double c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Octave(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Octave(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Octave(double c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Octave(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Octave(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Octave(N_Vector x, double b, N_Vector z);
SUNDIALS_EXPORT double N_VDotProd_Octave(N_Vector x, N_Vector y);
SUNDIALS_EXPORT double N_VMaxNorm_Octave(N_Vector x);
SUNDIALS_EXPORT double N_VWrmsNorm_Octave(N_Vector x, N_Vector w);
SUNDIALS_EXPORT double N_VWrmsNormMask_Octave(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT double N_VMin_Octave(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Octave(double c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_Octave(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_Octave(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT double N_VMinQuotient_Octave(N_Vector num, N_Vector denom);

// /* fused vector operations */
SUNDIALS_EXPORT int N_VScaleAddMulti_Octave(int nvec, double* a, N_Vector x,
                                            N_Vector* Y, N_Vector* Z);

/* vector array operations */
SUNDIALS_EXPORT int N_VLinearSumVectorArray_Octave(int nvec, 
                                                   double a, N_Vector* X,
                                                   double b, N_Vector* Y,
                                                   N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleVectorArray_Octave(int nvec, double* c,
                                               N_Vector* X, N_Vector* Z);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT double N_VWSqrSumLocal_Octave(N_Vector x, N_Vector w);
SUNDIALS_EXPORT double N_VWSqrSumMaskLocal_Octave(N_Vector x, N_Vector w, N_Vector id);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int N_VEnableFusedOps_Octave(N_Vector v, booleantype tf);


#ifdef __cplusplus
}
#endif
