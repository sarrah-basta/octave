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

// bool temper;
// /* Defining the content field */
// struct _N_VectorContent_Octave{
//     int length;  /* vector length       */
//     int own_data;           /* data ownership flag */
//     ColumnVector data_Octave;  /* our column vector */
//     double *data;            /* start of data array */
// };

/*
 * -----------------------------------------------------------------
 * Macros NV_CONTENT_S, NV_DATA_S, NV_OWN_DATA_S,
 *        NV_LENGTH_S, and NV_Ith_S
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_C(v)  ( (ColumnVector *)(v->content) )
#define NV_STCONTENT_C(v) ( (ColumnVector)(v->))

#define NV_LENGTH_C(v)   ( NV_CONTENT_C(v)->numel() )

// #define NV_OWN_DATA_C(v)  (temper)

#define NV_DATA_C(v)     ( NV_CONTENT_C(v)->fortran_vec() )
#define NV_Ith_C(v,i)     ( NV_DATA_C(v)[i])
typedef struct _N_VectorContent_Octave *N_VectorContent_Octave;

SUNDIALS_EXPORT N_Vector N_VNew_Octave(int vec_length, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewEmpty_Octave (SUNContext sunctx);
// SUNDIALS_EXPORT N_Vector N_VMake_Octave(sunindextype length, double *v_data, SUNContext sunctx);

SUNDIALS_EXPORT N_Vector* N_VCloneVectorArray_Octave(int count, N_Vector w);

SUNDIALS_EXPORT N_Vector* N_VCloneVectorArrayEmpty_Octave(int count, N_Vector w);

SUNDIALS_EXPORT void N_VDestroyVectorArray_Octave(N_Vector* vs, int count);

SUNDIALS_EXPORT sunindextype N_VGetLength_Octave(N_Vector v);

SUNDIALS_EXPORT void N_VPrint_Octave(N_Vector v);

SUNDIALS_EXPORT void N_VPrintFile_Octave(N_Vector v, FILE *outfile);

SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_Octave(N_Vector v);
SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Octave(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Octave(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Octave(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_Octave(N_Vector v, sunindextype *lrw, sunindextype *liw);
SUNDIALS_EXPORT double *N_VGetArrayPointer_Octave(N_Vector v);
SUNDIALS_EXPORT void N_VSetArrayPointer_Octave(double *v_data, N_Vector v);

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
SUNDIALS_EXPORT double N_VWL2Norm_Octave(N_Vector x, N_Vector w);
SUNDIALS_EXPORT double N_VL1Norm_Octave(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Octave(double c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_Octave(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_Octave(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT double N_VMinQuotient_Octave(N_Vector num, N_Vector denom);

// /* fused vector operations */
SUNDIALS_EXPORT int N_VLinearCombination_Octave(int nvec, double* c, N_Vector* V,
                                                N_Vector z);
SUNDIALS_EXPORT int N_VScaleAddMulti_Octave(int nvec, double* a, N_Vector x,
                                            N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT int N_VDotProdMulti_Octave(int nvec, N_Vector x,
                                           N_Vector* Y, double* dotprods);

/* vector array operations */
SUNDIALS_EXPORT int N_VLinearSumVectorArray_Octave(int nvec, 
                                                   double a, N_Vector* X,
                                                   double b, N_Vector* Y,
                                                   N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleVectorArray_Octave(int nvec, double* c,
                                               N_Vector* X, N_Vector* Z);
SUNDIALS_EXPORT int N_VConstVectorArray_Octave(int nvecs, double c,
                                               N_Vector* Z);
SUNDIALS_EXPORT int N_VWrmsNormVectorArray_Octave(int nvecs, N_Vector* X,
                                                  N_Vector* W, double* nrm);
SUNDIALS_EXPORT int N_VWrmsNormMaskVectorArray_Octave(int nvecs, N_Vector* X,
                                                      N_Vector* W, N_Vector id,
                                                      double* nrm);
SUNDIALS_EXPORT int N_VScaleAddMultiVectorArray_Octave(int nvec, int nsum,
                                                       double* a,
                                                       N_Vector* X,
                                                       N_Vector** Y,
                                                       N_Vector** Z);
SUNDIALS_EXPORT int N_VLinearCombinationVectorArray_Octave(int nvec, int nsum,
                                                           double* c,
                                                           N_Vector** X,
                                                           N_Vector* Z);

// /* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT double N_VWSqrSumLocal_Octave(N_Vector x, N_Vector w);
SUNDIALS_EXPORT double N_VWSqrSumMaskLocal_Octave(N_Vector x, N_Vector w, N_Vector id);

// /*
//  * -----------------------------------------------------------------
//  * Enable / disable fused vector operations
//  * -----------------------------------------------------------------
//  */

SUNDIALS_EXPORT int N_VEnableFusedOps_Octave(N_Vector v, booleantype tf);

// SUNDIALS_EXPORT int N_VEnableLinearCombination_Octave(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableScaleAddMulti_Octave(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableDotProdMulti_Octave(N_Vector v, booleantype tf);

// SUNDIALS_EXPORT int N_VEnableLinearSumVectorArray_Octave(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableScaleVectorArray_Octave(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableConstVectorArray_Octave(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableWrmsNormVectorArray_Octave(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableWrmsNormMaskVectorArray_Octave(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableScaleAddMultiVectorArray_Octave(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableLinearCombinationVectorArray_Octave(N_Vector v, booleantype tf);


#ifdef __cplusplus
}
#endif
