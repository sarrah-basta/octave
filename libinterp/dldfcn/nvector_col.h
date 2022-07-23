#include <stdio.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_config.h>
#include <dColVector.h>

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

#define NV_LENGTH_C(v)   ( NV_CONTENT_C(v)->numel() )

// #define NV_OWN_DATA_C(v)  (temper)

#define NV_DATA_C(v)     ( NV_CONTENT_C(v)->fortran_vec() )

#define NV_Ith_C(v,i)     ( NV_DATA_C(v)[i])

typedef struct _N_VectorContent_Col *N_VectorContent_Col;

SUNDIALS_EXPORT N_Vector N_VNew_Col(int vec_length, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewEmpty_Col (int length, SUNContext sunctx);
// SUNDIALS_EXPORT N_Vector N_VMake_Col(sunindextype length, double *v_data, SUNContext sunctx);

SUNDIALS_EXPORT N_Vector* N_VCloneVectorArray_Col(int count, N_Vector w);

// SUNDIALS_EXPORT N_Vector* N_VCloneVectorArrayEmpty_Col(int count, N_Vector w);

SUNDIALS_EXPORT void N_VDestroyVectorArray_Col(N_Vector* vs, int count);

SUNDIALS_EXPORT sunindextype N_VGetLength_Col(N_Vector v);

SUNDIALS_EXPORT void N_VPrint_Col(N_Vector v);

SUNDIALS_EXPORT void N_VPrintFile_Col(N_Vector v, FILE *outfile);

SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_Col(N_Vector v);
// SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Col(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Col(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Col(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_Col(N_Vector v, sunindextype *lrw, sunindextype *liw);
SUNDIALS_EXPORT double *N_VGetArrayPointer_Col(N_Vector v);
SUNDIALS_EXPORT void N_VSetArrayPointer_Col(double *v_data, N_Vector v);

// /* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_Col(double a, N_Vector x, double b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Col(double c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Col(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Col(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Col(double c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Col(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Col(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Col(N_Vector x, double b, N_Vector z);
SUNDIALS_EXPORT double N_VDotProd_Col(N_Vector x, N_Vector y);
SUNDIALS_EXPORT double N_VMaxNorm_Col(N_Vector x);
SUNDIALS_EXPORT double N_VWrmsNorm_Col(N_Vector x, N_Vector w);
SUNDIALS_EXPORT double N_VWrmsNormMask_Col(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT double N_VMin_Col(N_Vector x);
SUNDIALS_EXPORT double N_VWL2Norm_Col(N_Vector x, N_Vector w);
SUNDIALS_EXPORT double N_VL1Norm_Col(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Col(double c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_Col(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_Col(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT double N_VMinQuotient_Col(N_Vector num, N_Vector denom);

// /* fused vector operations */
SUNDIALS_EXPORT int N_VLinearCombination_Col(int nvec, double* c, N_Vector* V,
                                                N_Vector z);
SUNDIALS_EXPORT int N_VScaleAddMulti_Col(int nvec, double* a, N_Vector x,
                                            N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT int N_VDotProdMulti_Col(int nvec, N_Vector x,
                                           N_Vector* Y, double* dotprods);

/* vector array operations */
SUNDIALS_EXPORT int N_VLinearSumVectorArray_Col(int nvec, 
                                                   double a, N_Vector* X,
                                                   double b, N_Vector* Y,
                                                   N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleVectorArray_Col(int nvec, double* c,
                                               N_Vector* X, N_Vector* Z);
SUNDIALS_EXPORT int N_VConstVectorArray_Col(int nvecs, double c,
                                               N_Vector* Z);
SUNDIALS_EXPORT int N_VWrmsNormVectorArray_Col(int nvecs, N_Vector* X,
                                                  N_Vector* W, double* nrm);
SUNDIALS_EXPORT int N_VWrmsNormMaskVectorArray_Col(int nvecs, N_Vector* X,
                                                      N_Vector* W, N_Vector id,
                                                      double* nrm);
SUNDIALS_EXPORT int N_VScaleAddMultiVectorArray_Col(int nvec, int nsum,
                                                       double* a,
                                                       N_Vector* X,
                                                       N_Vector** Y,
                                                       N_Vector** Z);
SUNDIALS_EXPORT int N_VLinearCombinationVectorArray_Col(int nvec, int nsum,
                                                           double* c,
                                                           N_Vector** X,
                                                           N_Vector* Z);

// /* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT double N_VWSqrSumLocal_Col(N_Vector x, N_Vector w);
SUNDIALS_EXPORT double N_VWSqrSumMaskLocal_Col(N_Vector x, N_Vector w, N_Vector id);

// /*
//  * -----------------------------------------------------------------
//  * Enable / disable fused vector operations
//  * -----------------------------------------------------------------
//  */

SUNDIALS_EXPORT int N_VEnableFusedOps_Col(N_Vector v, booleantype tf);

// SUNDIALS_EXPORT int N_VEnableLinearCombination_Col(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableScaleAddMulti_Col(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableDotProdMulti_Col(N_Vector v, booleantype tf);

// SUNDIALS_EXPORT int N_VEnableLinearSumVectorArray_Col(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableScaleVectorArray_Col(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableConstVectorArray_Col(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableWrmsNormVectorArray_Col(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableWrmsNormMaskVectorArray_Col(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableScaleAddMultiVectorArray_Col(N_Vector v, booleantype tf);
// SUNDIALS_EXPORT int N_VEnableLinearCombinationVectorArray_Col(N_Vector v, booleantype tf);


#ifdef __cplusplus
}
#endif
