/* generating a custom N_Vector implementation */
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "nvector_col.h"
#include <sundials/sundials_math.h>

extern "C"
{

using namespace octave;

/* Defining required MACROS */
#define ZERO RCONST(0.0)
#define HALF RCONST(0.5)
#define ONE RCONST(1.0)
#define ONEPT5 RCONST(1.5)


  /*
   * -----------------------------------------------------------------
   * Functions exported by nvector_Col
   * -----------------------------------------------------------------
   */

  /* Private functions for special cases of vector operations */
  static void VCopy_Col(N_Vector x, N_Vector z);                              /* z=x       */
  static void VSum_Col(N_Vector x, N_Vector y, N_Vector z);                   /* z=x+y     */
  static void VDiff_Col(N_Vector x, N_Vector y, N_Vector z);                  /* z=x-y     */
  static void VNeg_Col(N_Vector x, N_Vector z);                               /* z=-x      */
  static void VScaleSum_Col(realtype c, N_Vector x, N_Vector y, N_Vector z);  /* z=c(x+y)  */
  static void VScaleDiff_Col(realtype c, N_Vector x, N_Vector y, N_Vector z); /* z=c(x-y)  */
  static void VLin1_Col(realtype a, N_Vector x, N_Vector y, N_Vector z);      /* z=ax+y    */
  static void VLin2_Col(realtype a, N_Vector x, N_Vector y, N_Vector z);      /* z=ax-y    */
  static void Vaxpy_Col(realtype a, N_Vector x, N_Vector y);                  /* y <- ax+y */
  static void VScaleBy_Col(realtype a, N_Vector x);                           /* x <- ax   */

  /* Private functions for special cases of vector array operations */
  static int VSumVectorArray_Col(int nvec, N_Vector *X, N_Vector *Y, N_Vector *Z);                   /* Z=X+Y     */
  static int VDiffVectorArray_Col(int nvec, N_Vector *X, N_Vector *Y, N_Vector *Z);                  /* Z=X-Y     */
  static int VScaleSumVectorArray_Col(int nvec, realtype c, N_Vector *X, N_Vector *Y, N_Vector *Z);  /* Z=c(X+Y)  */
  static int VScaleDiffVectorArray_Col(int nvec, realtype c, N_Vector *X, N_Vector *Y, N_Vector *Z); /* Z=c(X-Y)  */
  static int VLin1VectorArray_Col(int nvec, realtype a, N_Vector *X, N_Vector *Y, N_Vector *Z);      /* Z=aX+Y    */
  static int VLin2VectorArray_Col(int nvec, realtype a, N_Vector *X, N_Vector *Y, N_Vector *Z);      /* Z=aX-Y    */
  static int VaxpyVectorArray_Col(int nvec, realtype a, N_Vector *X, N_Vector *Y);                   /* Y <- aX+Y */

  /*
   * -----------------------------------------------------------------
   * exported functions
   * -----------------------------------------------------------------
   */

  /* ----------------------------------------------------------------
   * Returns vector type ID. Used to identify vector implementation
   * from abstract N_Vector interface.
   */
  N_Vector_ID N_VGetVectorID_Col(N_Vector v)
  {
    return SUNDIALS_NVEC_SERIAL;
  }

  /* ----------------------------------------------------------------------------
   * Function to create a new empty Col vector
   */

  N_Vector N_VNew_Col(int length, SUNContext sunctx)
  {
    N_Vector v;
    void *content;
    // realtype *data;

    if (sunctx == NULL)
      return (NULL);

    /* Create an empty vector object */
    v = NULL;
    v = N_VNewEmpty(sunctx);
    if (v == NULL)
      return (NULL);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = N_VGetVectorID_Col;
  v->ops->nvclone           = N_VClone_Col;
  v->ops->nvcloneempty      = NULL;
  v->ops->nvdestroy         = N_VDestroy_Col;
  v->ops->nvspace           = N_VSpace_Col;
  v->ops->nvgetarraypointer = N_VGetArrayPointer_Col;
  v->ops->nvsetarraypointer = N_VSetArrayPointer_Col;
  v->ops->nvgetlength       = N_VGetLength_Col;

  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_Col;
  v->ops->nvconst        = N_VConst_Col;
  v->ops->nvprod         = N_VProd_Col;
  v->ops->nvdiv          = N_VDiv_Col;
  v->ops->nvscale        = N_VScale_Col;
  v->ops->nvabs          = N_VAbs_Col;
  v->ops->nvinv          = N_VInv_Col;
  v->ops->nvaddconst     = N_VAddConst_Col;
  v->ops->nvdotprod      = N_VDotProd_Col;
  v->ops->nvmaxnorm      = N_VMaxNorm_Col;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_Col;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_Col;
  v->ops->nvmin          = N_VMin_Col;
  v->ops->nvwl2norm      = N_VWL2Norm_Col;
  v->ops->nvl1norm       = N_VL1Norm_Col;
  v->ops->nvcompare      = N_VCompare_Col;
  v->ops->nvinvtest      = N_VInvTest_Col;
  v->ops->nvconstrmask   = N_VConstrMask_Col;
  v->ops->nvminquotient  = N_VMinQuotient_Col;

  /* fused and vector array operations are disabled (NULL) by default */

  /*
    * OPTIONAL operations.
    *
    * These operations provide default implementations that may be overriden.
    */

  /* fused vector operations (optional) */
  v->ops->nvlinearcombination = N_VLinearCombination_Col;
  v->ops->nvscaleaddmulti = N_VScaleAddMulti_Col;
  v->ops->nvdotprodmulti = N_VDotProdMulti_Col;

  // /* vector array operations (optional) */
  v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Col;
  v->ops->nvscalevectorarray = N_VScaleVectorArray_Col;
  v->ops->nvconstvectorarray = N_VConstVectorArray_Col;
  v->ops->nvwrmsnormvectorarray = N_VWrmsNormVectorArray_Col;
  v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_Col;
  v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Col;
  v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Col;

  // /*
  //  * OPTIONAL operations with no default implementation.
  //  */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = NULL;
  v->ops->nvmaxnormlocal     = NULL;
  v->ops->nvminlocal         = NULL;
  v->ops->nvl1normlocal      = NULL;
  v->ops->nvinvtestlocal     = NULL;
  v->ops->nvconstrmasklocal  = N_VConstrMask_Col;
  v->ops->nvminquotientlocal = NULL;
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Col;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Col;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal = NULL;
  v->ops->nvdotprodmultiallreduce = NULL;

  /* XBraid interface operations */
  v->ops->nvbufsize = NULL;
  v->ops->nvbufpack = NULL;
  v->ops->nvbufunpack = NULL;

  // /* debugging functions */
  v->ops->nvprint = N_VPrint_Col;
  v->ops->nvprintfile = N_VPrintFile_Col;

  /* Create data */
  if (length > 0) {

    /* Allocate memory */
    /* Create content */
    content = NULL;
    content = new ColumnVector(length);
    if (content == NULL)
    {
      N_VDestroy(v);
      return (NULL);
    }
    /* Attach content */
    v->content = content;
    if(NV_DATA_C(v) == NULL) { N_VDestroy_Col(v); return(NULL); }

    /* Attach data */
    // NV_OWN_DATA_C(v) = SUNTRUE;

  }
  return (v);
}


/* ----------------------------------------------------------------------------
 * Function to create an array of new serial vectors.
 */

N_Vector* N_VCloneVectorArray_Col(int count, N_Vector w)
{
  N_Vector* vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector*) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VClone_Col(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Col(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_Col
 */

void N_VDestroyVectorArray_Col(N_Vector* vs, int count)
{
  int j;

  for (j = 0; j < count; j++) N_VDestroy_Col(vs[j]);

  free(vs); vs = NULL;

  return;
}

/* ----------------------------------------------------------------------------
 * Function to return number of vector elements
 */
sunindextype N_VGetLength_Col(N_Vector v)
{
  return NV_LENGTH_C(v);
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to stdout
 */

void N_VPrint_Col(N_Vector x)
{
  N_VPrintFile_Col(x, stdout);
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to outfile
 */

void N_VPrintFile_Col(N_Vector x, FILE* outfile)
{
  sunindextype i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);

  for (i = 0; i < N; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%35.32Lg\n", xd[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%19.16g\n", xd[i]);
#else
    fprintf(outfile, "%11.8g\n", xd[i]);
#endif
  }
  fprintf(outfile, "\n");

  return;
}

// /*
//  * -----------------------------------------------------------------
//  * implementation of vector operations
//  * -----------------------------------------------------------------
//  */

N_Vector N_VClone_Col(N_Vector w)
{
  N_Vector v;
  void* content;
  sunindextype length;

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(w->sunctx);
  if (v == NULL) return(NULL);

  /* Attach operations */
  if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

  /* Create content */
  content = NULL;
  content = new ColumnVector(NV_LENGTH_C(w));
  if (content == NULL) { N_VDestroy(v); return(NULL); }

  /* Attach content */
  v->content = content;
  /* Attach data */
  // NV_OWN_DATA_C(v) = SUNTRUE;

  return(v);
}

void N_VDestroy_Col(N_Vector v)
{
  if (v == NULL) return;

  /* free content */
  if (v->content != NULL) {
    // /* free data array if it's owned by the vector */
    // if (NV_OWN_DATA_C(v) && NV_DATA_C(v) != NULL) {
    //   free(NV_DATA_C(v));
    //   NV_DATA_C(v) = NULL;
    // }
    free(v->content);
    v->content = NULL;
  }

  /* free ops and vector */
  if (v->ops != NULL) { free(v->ops); v->ops = NULL; }
  free(v); v = NULL;

  return;
}

void N_VSpace_Col(N_Vector v, sunindextype *lrw, sunindextype *liw)
{
  *lrw = NV_LENGTH_C(v);
  *liw = 1;

  return;
}

realtype *N_VGetArrayPointer_Col(N_Vector v)
{
  return((realtype *) NV_DATA_C(v));
}

void N_VSetArrayPointer_Col(realtype* v_data, N_Vector v)
{
  // if (NV_LENGTH_C(v) > 0) NV_DATA_C(v) = v_data;
  

  return;
}

void N_VdrrayPointer_Col(realtype *v_data, N_Vector v)
{
  if (NV_LENGTH_C(v) > 0) v_data = NV_CONTENT_C(v)->fortran_vec();

  return;
}

void N_VLinearSum_Col(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype c, *xd, *yd, *zd;
  N_Vector v1, v2;
  booleantype test;

  xd = yd = zd = NULL;

  if ((b == ONE) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
    Vaxpy_Col(a,x,y);
    return;
  }

  if ((a == ONE) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
    Vaxpy_Col(b,y,x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE)) {
    VSum_Col(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_Col(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE)) {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_Col(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_Col(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b) {
    VScaleSum_Col(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b) {
    VScaleDiff_Col(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  yd = NV_DATA_C(y);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])+(b*yd[i]);

  return;
}

void N_VConst_Col(realtype c, N_Vector z)
{
  sunindextype i, N;
  realtype *zd;

  zd = NULL;

  N  = NV_LENGTH_C(z);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++) zd[i] = c;

  return;
}

void N_VProd_Col(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  yd = NV_DATA_C(y);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]*yd[i];

  return;
}

void N_VDiv_Col(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  yd = NV_DATA_C(y);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]/yd[i];

  return;
}

void N_VScale_Col(realtype c, N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  if (z == x) {  /* BLAS usage: scale x <- cx */
    VScaleBy_Col(c, x);
    return;
  }

  if (c == ONE) {
    VCopy_Col(x, z);
  } else if (c == -ONE) {
    VNeg_Col(x, z);
  } else {
    N  = NV_LENGTH_C(x);
    xd = NV_DATA_C(x);
    zd = NV_DATA_C(z);
    for (i = 0; i < N; i++)
      zd[i] = c*xd[i];
  }

  return;
}

void N_VAbs_Col(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++)
    zd[i] = SUNRabs(xd[i]);

  return;
}

void N_VInv_Col(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++)
    zd[i] = ONE/xd[i];

  return;
}

void N_VAddConst_Col(N_Vector x, realtype b, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]+b;

  return;
}

realtype N_VDotProd_Col(N_Vector x, N_Vector y)
{
  sunindextype i, N;
  realtype sum, *xd, *yd;

  sum = ZERO;
  xd = yd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  yd = NV_DATA_C(y);

  for (i = 0; i < N; i++)
    sum += xd[i]*yd[i];

  return(sum);
}

realtype N_VMaxNorm_Col(N_Vector x)
{
  sunindextype i, N;
  realtype max, *xd;

  max = ZERO;
  xd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);

  for (i = 0; i < N; i++) {
    if (SUNRabs(xd[i]) > max) max = SUNRabs(xd[i]);
  }

  return(max);
}

realtype N_VWrmsNorm_Col(N_Vector x, N_Vector w)
{
  return(SUNRsqrt(N_VWSqrSumLocal_Col(x, w)/(NV_LENGTH_C(x))));
}

realtype N_VWSqrSumLocal_Col(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  realtype sum, prodi, *xd, *wd;

  sum = ZERO;
  xd = wd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  wd = NV_DATA_C(w);

  for (i = 0; i < N; i++) {
    prodi = xd[i]*wd[i];
    sum += SUNSQR(prodi);
  }

  return(sum);
}

realtype N_VWrmsNormMask_Col(N_Vector x, N_Vector w, N_Vector id)
{
  return(SUNRsqrt(N_VWSqrSumMaskLocal_Col(x, w, id) / (NV_LENGTH_C(x))));
}

realtype N_VWSqrSumMaskLocal_Col(N_Vector x, N_Vector w, N_Vector id)
{
  sunindextype i, N;
  realtype sum, prodi, *xd, *wd, *idd;

  sum = ZERO;
  xd = wd = idd = NULL;

  N  = NV_LENGTH_C(x);
  xd  = NV_DATA_C(x);
  wd  = NV_DATA_C(w);
  idd = NV_DATA_C(id);

  for (i = 0; i < N; i++) {
    if (idd[i] > ZERO) {
      prodi = xd[i]*wd[i];
      sum += SUNSQR(prodi);
    }
  }

  return(sum);
}

realtype N_VMin_Col(N_Vector x)
{
  sunindextype i, N;
  realtype min, *xd;

  xd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);

  min = xd[0];

  for (i = 1; i < N; i++) {
    if (xd[i] < min) min = xd[i];
  }

  return(min);
}

realtype N_VWL2Norm_Col(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  realtype sum, prodi, *xd, *wd;

  sum = ZERO;
  xd = wd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  wd = NV_DATA_C(w);

  for (i = 0; i < N; i++) {
    prodi = xd[i]*wd[i];
    sum += SUNSQR(prodi);
  }

  return(SUNRsqrt(sum));
}

realtype N_VL1Norm_Col(N_Vector x)
{
  sunindextype i, N;
  realtype sum, *xd;

  sum = ZERO;
  xd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);

  for (i = 0; i<N; i++)
    sum += SUNRabs(xd[i]);

  return(sum);
}

void N_VCompare_Col(realtype c, N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++) {
    zd[i] = (SUNRabs(xd[i]) >= c) ? ONE : ZERO;
  }

  return;
}

booleantype N_VInvTest_Col(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;
  booleantype no_zero_found;

  xd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  zd = NV_DATA_C(z);

  no_zero_found = SUNTRUE;
  for (i = 0; i < N; i++) {
    if (xd[i] == ZERO)
      no_zero_found = SUNFALSE;
    else
      zd[i] = ONE/xd[i];
  }

  return no_zero_found;
}

booleantype N_VConstrMask_Col(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype i, N;
  realtype temp;
  realtype *cd, *xd, *md;
  booleantype test;

  cd = xd = md = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  cd = NV_DATA_C(c);
  md = NV_DATA_C(m);

  temp = ZERO;

  for (i = 0; i < N; i++) {
    md[i] = ZERO;

    /* Continue if no constraints were set for the variable */
    if (cd[i] == ZERO)
      continue;

    /* Check if a set constraint has been violated */
    test = (SUNRabs(cd[i]) > ONEPT5 && xd[i]*cd[i] <= ZERO) ||
           (SUNRabs(cd[i]) > HALF   && xd[i]*cd[i] <  ZERO);
    if (test) {
      temp = md[i] = ONE;
    }
  }

  /* Return false if any constraint was violated */
  return (temp == ONE) ? SUNFALSE : SUNTRUE;
}

realtype N_VMinQuotient_Col(N_Vector num, N_Vector denom)
{
  booleantype notEvenOnce;
  sunindextype i, N;
  realtype *nd, *dd, min;

  nd = dd = NULL;

  N  = NV_LENGTH_C(num);
  nd = NV_DATA_C(num);
  dd = NV_DATA_C(denom);

  notEvenOnce = SUNTRUE;
  min = BIG_REAL;

  for (i = 0; i < N; i++) {
    if (dd[i] == ZERO) continue;
    else {
      if (!notEvenOnce) min = SUNMIN(min, nd[i]/dd[i]);
      else {
	min = nd[i]/dd[i];
        notEvenOnce = SUNFALSE;
      }
    }
  }

  return(min);
}

// /*
//  * -----------------------------------------------------------------
//  * fused vector operations
//  * -----------------------------------------------------------------
//  */

int N_VLinearCombination_Col(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  int          i;
  sunindextype j, N;
  realtype*    zd=NULL;
  realtype*    xd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_Col(c[0], X[0], z);
    return(0);
  }

  /* should have called N_VLinearSum */
  if (nvec == 2) {
    N_VLinearSum_Col(c[0], X[0], c[1], X[1], z);
    return(0);
  }

  /* get vector length and data array */
  N  = NV_LENGTH_C(z);
  zd = NV_DATA_C(z);

  /*
   * X[0] += c[i]*X[i], i = 1,...,nvec-1
   */
  if ((X[0] == z) && (c[0] == ONE)) {
    for (i=1; i<nvec; i++) {
      xd = NV_DATA_C(X[i]);
      for (j=0; j<N; j++) {
        zd[j] += c[i] * xd[j];
      }
    }
    return(0);
  }

  /*
   * X[0] = c[0] * X[0] + sum{ c[i] * X[i] }, i = 1,...,nvec-1
   */
  if (X[0] == z) {
    for (j=0; j<N; j++) {
      zd[j] *= c[0];
    }
    for (i=1; i<nvec; i++) {
      xd = NV_DATA_C(X[i]);
      for (j=0; j<N; j++) {
        zd[j] += c[i] * xd[j];
      }
    }
    return(0);
  }

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */
  xd = NV_DATA_C(X[0]);
  for (j=0; j<N; j++) {
    zd[j] = c[0] * xd[j];
  }
  for (i=1; i<nvec; i++) {
    xd = NV_DATA_C(X[i]);
    for (j=0; j<N; j++) {
      zd[j] += c[i] * xd[j];
    }
  }
  return(0);
}


int N_VScaleAddMulti_Col(int nvec, realtype* a, N_Vector x, N_Vector* Y, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VLinearSum */
  if (nvec == 1) {
    N_VLinearSum_Col(a[0], x, ONE, Y[0], Z[0]);
    return(0);
  }

  /* get vector length and data array */
  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z) {
    for (i=0; i<nvec; i++) {
      yd = NV_DATA_C(Y[i]);
      for (j=0; j<N; j++) {
        yd[j] += a[i] * xd[j];
      }
    }
    return(0);
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i=0; i<nvec; i++) {
    yd = NV_DATA_C(Y[i]);
    zd = NV_DATA_C(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = a[i] * xd[j] + yd[j];
    }
  }
  return(0);
}


int N_VDotProdMulti_Col(int nvec, N_Vector x, N_Vector* Y, realtype* dotprods)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VDotProd */
  if (nvec == 1) {
    dotprods[0] = N_VDotProd_Col(x, Y[0]);
    return(0);
  }

  /* get vector length and data array */
  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);

  /* compute multiple dot products */
  for (i=0; i<nvec; i++) {
    yd = NV_DATA_C(Y[i]);
    dotprods[i] = ZERO;
    for (j=0; j<N; j++) {
      dotprods[i] += xd[j] * yd[j];
    }
  }

  return(0);
}


/*
 * -----------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------
 */

int N_VLinearSumVectorArray_Col(int nvec,
                                   realtype a, N_Vector* X,
                                   realtype b, N_Vector* Y,
                                   N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;
  realtype     c;
  N_Vector*   V1;
  N_Vector*   V2;
  booleantype  test;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VLinearSum */
  if (nvec == 1) {
    N_VLinearSum_Col(a, X[0], b, Y[0], Z[0]);
    return(0);
  }

  /* BLAS usage: axpy y <- ax+y */
  if ((b == ONE) && (Z == Y))
    return(VaxpyVectorArray_Col(nvec, a, X, Y));

  /* BLAS usage: axpy x <- by+x */
  if ((a == ONE) && (Z == X))
    return(VaxpyVectorArray_Col(nvec, b, Y, X));

  /* Case: a == b == 1.0 */
  if ((a == ONE) && (b == ONE))
    return(VSumVectorArray_Col(nvec, X, Y, Z));

  /* Cases:                    */
  /*   (1) a == 1.0, b = -1.0, */
  /*   (2) a == -1.0, b == 1.0 */
  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    return(VDiffVectorArray_Col(nvec, V2, V1, Z));
  }

  /* Cases:                                                  */
  /*   (1) a == 1.0, b == other or 0.0,                      */
  /*   (2) a == other or 0.0, b == 1.0                       */
  /* if a or b is 0.0, then user should have called N_VScale */
  if ((test = (a == ONE)) || (b == ONE)) {
    c  = test ? b : a;
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    return(VLin1VectorArray_Col(nvec, c, V1, V2, Z));
  }

  /* Cases:                     */
  /*   (1) a == -1.0, b != 1.0, */
  /*   (2) a != 1.0, b == -1.0  */
  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    return(VLin2VectorArray_Col(nvec, c, V1, V2, Z));
  }

  /* Case: a == b                                                         */
  /* catches case both a and b are 0.0 - user should have called N_VConst */
  if (a == b)
    return(VScaleSumVectorArray_Col(nvec, a, X, Y, Z));

  /* Case: a == -b */
  if (a == -b)
    return(VScaleDiffVectorArray_Col(nvec, a, X, Y, Z));

  /* Do all cases not handled above:                               */
  /*   (1) a == other, b == 0.0 - user should have called N_VScale */
  /*   (2) a == 0.0, b == other - user should have called N_VScale */
  /*   (3) a,b == other, a !=b, a != -b                            */

  /* get vector length */
  N = NV_LENGTH_C(Z[0]);

  /* compute linear sum for each vector pair in vector arrays */
  for (i=0; i<nvec; i++) {
    xd = NV_DATA_C(X[i]);
    yd = NV_DATA_C(Y[i]);
    zd = NV_DATA_C(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = a * xd[j] + b * yd[j];
    }
  }

  return(0);
}


int N_VScaleVectorArray_Col(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    zd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_Col(c[0], X[0], Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_LENGTH_C(Z[0]);

  /*
   * X[i] *= c[i]
   */
  if (X == Z) {
    for (i=0; i<nvec; i++) {
      xd = NV_DATA_C(X[i]);
      for (j=0; j<N; j++) {
        xd[j] *= c[i];
      }
    }
    return(0);
  }

  /*
   * Z[i] = c[i] * X[i]
   */
  for (i=0; i<nvec; i++) {
    xd = NV_DATA_C(X[i]);
    zd = NV_DATA_C(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = c[i] * xd[j];
    }
  }
  return(0);
}


int N_VConstVectorArray_Col(int nvec, realtype c, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    zd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VConst */
  if (nvec == 1) {
    N_VConst_Col(c, Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_LENGTH_C(Z[0]);

  /* set each vector in the vector array to a constant */
  for (i=0; i<nvec; i++) {
    zd = NV_DATA_C(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = c;
    }
  }

  return(0);
}


int N_VWrmsNormVectorArray_Col(int nvec, N_Vector* X, N_Vector* W,
                                  realtype* nrm)
{
  int          i;
  sunindextype j, N;
  realtype*    wd=NULL;
  realtype*    xd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VWrmsNorm */
  if (nvec == 1) {
    nrm[0] = N_VWrmsNorm_Col(X[0], W[0]);
    return(0);
  }

  /* get vector length */
  N = NV_LENGTH_C(X[0]);

  /* compute the WRMS norm for each vector in the vector array */
  for (i=0; i<nvec; i++) {
    xd = NV_DATA_C(X[i]);
    wd = NV_DATA_C(W[i]);
    nrm[i] = ZERO;
    for (j=0; j<N; j++) {
      nrm[i] += SUNSQR(xd[j] * wd[j]);
    }
    nrm[i] = SUNRsqrt(nrm[i]/N);
  }

  return(0);
}


int N_VWrmsNormMaskVectorArray_Col(int nvec, N_Vector* X, N_Vector* W,
                                      N_Vector id, realtype* nrm)
{
  int          i;
  sunindextype j, N;
  realtype*    wd=NULL;
  realtype*    xd=NULL;
  realtype*    idd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VWrmsNorm */
  if (nvec == 1) {
    nrm[0] = N_VWrmsNormMask_Col(X[0], W[0], id);
    return(0);
  }

  /* get vector length and mask data array */
  N   = NV_LENGTH_C(X[0]);
  idd = NV_DATA_C(id);

  /* compute the WRMS norm for each vector in the vector array */
  for (i=0; i<nvec; i++) {
    xd = NV_DATA_C(X[i]);
    wd = NV_DATA_C(W[i]);
    nrm[i] = ZERO;
    for (j=0; j<N; j++) {
      if (idd[j] > ZERO)
        nrm[i] += SUNSQR(xd[j] * wd[j]);
    }
    nrm[i] = SUNRsqrt(nrm[i]/N);
  }

  return(0);
}


int N_VScaleAddMultiVectorArray_Col(int nvec, int nsum, realtype* a,
                                        N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  int          i, j;
  sunindextype k, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  int          retval;
  N_Vector*   YY;
  N_Vector*   ZZ;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);
  if (nsum < 1) return(-1);

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1) {

    /* should have called N_VLinearSum */
    if (nsum == 1) {
      N_VLinearSum_Col(a[0], X[0], ONE, Y[0][0], Z[0][0]);
      return(0);
    }

    /* should have called N_VScaleAddMulti */
    YY = (N_Vector*) malloc(nsum * sizeof(N_Vector));
    ZZ = (N_Vector*) malloc(nsum * sizeof(N_Vector));

    for (j=0; j<nsum; j++) {
      YY[j] = Y[j][0];
      ZZ[j] = Z[j][0];
    }

    retval = N_VScaleAddMulti_Col(nsum, a, X[0], YY, ZZ);

    free(YY);
    free(ZZ);
    return(retval);
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 1) {
    retval = N_VLinearSumVectorArray_Col(nvec, a[0], X, ONE, Y[0], Z[0]);
    return(retval);
  }

  /* ----------------------------
   * Compute multiple linear sums
   * ---------------------------- */

  /* get vector length */
  N = NV_LENGTH_C(X[0]);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z) {
    for (i=0; i<nvec; i++) {
      xd = NV_DATA_C(X[i]);
      for (j=0; j<nsum; j++){
        yd = NV_DATA_C(Y[j][i]);
        for (k=0; k<N; k++) {
          yd[k] += a[j] * xd[k];
        }
      }
    }
    return(0);
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i=0; i<nvec; i++) {
    xd = NV_DATA_C(X[i]);
    for (j=0; j<nsum; j++){
      yd = NV_DATA_C(Y[j][i]);
      zd = NV_DATA_C(Z[j][i]);
      for (k=0; k<N; k++) {
        zd[k] = a[j] * xd[k] + yd[k];
      }
    }
  }
  return(0);
}


int N_VLinearCombinationVectorArray_Col(int nvec, int nsum, realtype* c,
                                           N_Vector** X, N_Vector* Z)
{
  int          i; /* vector arrays index in summation [0,nsum) */
  int          j; /* vector index in vector array     [0,nvec) */
  sunindextype k; /* element index in vector          [0,N)    */
  sunindextype N;
  realtype*    zd=NULL;
  realtype*    xd=NULL;

  int          retval;
  realtype*    ctmp;
  N_Vector*   Y;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);
  if (nsum < 1) return(-1);

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1) {

    /* should have called N_VScale */
    if (nsum == 1) {
      N_VScale_Col(c[0], X[0][0], Z[0]);
      return(0);
    }

    /* should have called N_VLinearSum */
    if (nsum == 2) {
      N_VLinearSum_Col(c[0], X[0][0], c[1], X[1][0], Z[0]);
      return(0);
    }

    /* should have called N_VLinearCombination */
    Y = (N_Vector*) malloc(nsum * sizeof(N_Vector));

    for (i=0; i<nsum; i++) {
      Y[i] = X[i][0];
    }

    retval = N_VLinearCombination_Col(nsum, c, Y, Z[0]);

    free(Y);
    return(retval);
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VScaleVectorArray */
  if (nsum == 1) {

    ctmp = (realtype*) malloc(nvec * sizeof(realtype));

    for (j=0; j<nvec; j++) {
      ctmp[j] = c[0];
    }

    retval = N_VScaleVectorArray_Col(nvec, ctmp, X[0], Z);

    free(ctmp);
    return(retval);
  }

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 2) {
    retval = N_VLinearSumVectorArray_Col(nvec, c[0], X[0], c[1], X[1], Z);
    return(retval);
  }

  /* --------------------------
   * Compute linear combination
   * -------------------------- */

  /* get vector length */
  N = NV_LENGTH_C(Z[0]);

  /*
   * X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
   */
  if ((X[0] == Z) && (c[0] == ONE)) {
    for (j=0; j<nvec; j++) {
      zd = NV_DATA_C(Z[j]);
      for (i=1; i<nsum; i++) {
        xd = NV_DATA_C(X[i][j]);
        for (k=0; k<N; k++) {
          zd[k] += c[i] * xd[k];
        }
      }
    }
    return(0);
  }

  /*
   * X[0][j] = c[0] * X[0][j] + sum{ c[i] * X[i][j] }, i = 1,...,nvec-1
   */
  if (X[0] == Z) {
    for (j=0; j<nvec; j++) {
      zd = NV_DATA_C(Z[j]);
      for (k=0; k<N; k++) {
        zd[k] *= c[0];
      }
      for (i=1; i<nsum; i++) {
        xd = NV_DATA_C(X[i][j]);
        for (k=0; k<N; k++) {
          zd[k] += c[i] * xd[k];
        }
      }
    }
    return(0);
  }

  /*
   * Z[j] = sum{ c[i] * X[i][j] }, i = 0,...,nvec-1
   */
  for (j=0; j<nvec; j++) {
    xd = NV_DATA_C(X[0][j]);
    zd = NV_DATA_C(Z[j]);
    for (k=0; k<N; k++) {
      zd[k] = c[0] * xd[k];
    }
    for (i=1; i<nsum; i++) {
      xd = NV_DATA_C(X[i][j]);
      for (k=0; k<N; k++) {
        zd[k] += c[i] * xd[k];
      }
    }
  }
  return(0);
}


// /*
//  * -----------------------------------------------------------------
//  * private functions for special cases of vector operations
//  * -----------------------------------------------------------------
//  */

static void VCopy_Col(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i];

  return;
}

static void VSum_Col(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  yd = NV_DATA_C(y);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]+yd[i];

  return;
}

static void VDiff_Col(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  yd = NV_DATA_C(y);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]-yd[i];

  return;
}

static void VNeg_Col(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++)
    zd[i] = -xd[i];

  return;
}

static void VScaleSum_Col(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  yd = NV_DATA_C(y);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++)
    zd[i] = c*(xd[i]+yd[i]);

  return;
}

static void VScaleDiff_Col(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  yd = NV_DATA_C(y);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++)
    zd[i] = c*(xd[i]-yd[i]);

  return;
}

static void VLin1_Col(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  yd = NV_DATA_C(y);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])+yd[i];

  return;
}

static void VLin2_Col(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  yd = NV_DATA_C(y);
  zd = NV_DATA_C(z);

  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])-yd[i];

  return;
}

static void Vaxpy_Col(realtype a, N_Vector x, N_Vector y)
{
  sunindextype i, N;
  realtype *xd, *yd;

  xd = yd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  yd = NV_DATA_C(y);

  if (a == ONE) {
    for (i = 0; i < N; i++)
      yd[i] += xd[i];
    return;
  }

  if (a == -ONE) {
    for (i = 0; i < N; i++)
      yd[i] -= xd[i];
    return;
  }

  for (i = 0; i < N; i++)
    yd[i] += a*xd[i];

  return;
}

static void VScaleBy_Col(realtype a, N_Vector x)
{
  sunindextype i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);

  for (i = 0; i < N; i++)
    xd[i] *= a;

  return;
}


/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector array operations
 * -----------------------------------------------------------------
 */

static int VSumVectorArray_Col(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  N = NV_LENGTH_C(X[0]);

  for (i=0; i<nvec; i++) {
    xd = NV_DATA_C(X[i]);
    yd = NV_DATA_C(Y[i]);
    zd = NV_DATA_C(Z[i]);
    for (j=0; j<N; j++)
      zd[j] = xd[j] + yd[j];
  }

  return(0);
}

static int VDiffVectorArray_Col(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  N = NV_LENGTH_C(X[0]);

  for (i=0; i<nvec; i++) {
    xd = NV_DATA_C(X[i]);
    yd = NV_DATA_C(Y[i]);
    zd = NV_DATA_C(Z[i]);
    for (j=0; j<N; j++)
      zd[j] = xd[j] - yd[j];
  }

  return(0);
}

static int VScaleSumVectorArray_Col(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  N = NV_LENGTH_C(X[0]);

  for (i=0; i<nvec; i++) {
    xd = NV_DATA_C(X[i]);
    yd = NV_DATA_C(Y[i]);
    zd = NV_DATA_C(Z[i]);
    for (j=0; j<N; j++)
      zd[j] = c * (xd[j] + yd[j]);
  }

  return(0);
}

static int VScaleDiffVectorArray_Col(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  N = NV_LENGTH_C(X[0]);

  for (i=0; i<nvec; i++) {
    xd = NV_DATA_C(X[i]);
    yd = NV_DATA_C(Y[i]);
    zd = NV_DATA_C(Z[i]);
    for (j=0; j<N; j++)
      zd[j] = c * (xd[j] - yd[j]);
  }

  return(0);
}

static int VLin1VectorArray_Col(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  N = NV_LENGTH_C(X[0]);

  for (i=0; i<nvec; i++) {
    xd = NV_DATA_C(X[i]);
    yd = NV_DATA_C(Y[i]);
    zd = NV_DATA_C(Z[i]);
    for (j=0; j<N; j++)
      zd[j] = (a * xd[j]) + yd[j];
  }

  return(0);
}

static int VLin2VectorArray_Col(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  N = NV_LENGTH_C(X[0]);

  for (i=0; i<nvec; i++) {
    xd = NV_DATA_C(X[i]);
    yd = NV_DATA_C(Y[i]);
    zd = NV_DATA_C(Z[i]);
    for (j=0; j<N; j++)
      zd[j] = (a * xd[j]) - yd[j];
  }

  return(0);
}

static int VaxpyVectorArray_Col(int nvec, realtype a, N_Vector* X, N_Vector* Y)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;

  N = NV_LENGTH_C(X[0]);

  if (a == ONE) {
    for (i=0; i<nvec; i++) {
      xd = NV_DATA_C(X[i]);
      yd = NV_DATA_C(Y[i]);
      for (j=0; j<N; j++)
        yd[j] += xd[j];
    }

    return(0);
  }

  if (a == -ONE) {
    for (i=0; i<nvec; i++) {
      xd = NV_DATA_C(X[i]);
      yd = NV_DATA_C(Y[i]);
      for (j=0; j<N; j++)
        yd[j] -= xd[j];
    }

    return(0);
  }    

  for (i=0; i<nvec; i++) {
    xd = NV_DATA_C(X[i]);
    yd = NV_DATA_C(Y[i]);
    for (j=0; j<N; j++)
      yd[j] += a * xd[j];
  }

  return(0);
}


// /*
//  * -----------------------------------------------------------------
//  * Enable / Disable fused and vector array operations
//  * -----------------------------------------------------------------
//  */

int N_VEnableFusedOps_Col(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_Col;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Col;
    v->ops->nvdotprodmulti      = N_VDotProdMulti_Col;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Col;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_Col;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_Col;
    v->ops->nvwrmsnormvectorarray          = N_VWrmsNormVectorArray_Col;
    v->ops->nvwrmsnormmaskvectorarray      = N_VWrmsNormMaskVectorArray_Col;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Col;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Col;
  } else {
    /* disable all fused vector operations */
    v->ops->nvlinearcombination = NULL;
    v->ops->nvscaleaddmulti     = NULL;
    v->ops->nvdotprodmulti      = NULL;
    /* disable all vector array operations */
    v->ops->nvlinearsumvectorarray         = NULL;
    v->ops->nvscalevectorarray             = NULL;
    v->ops->nvconstvectorarray             = NULL;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = NULL;
    v->ops->nvlinearcombinationvectorarray = NULL;
  }

  /* return success */
  return(0);
}

} /* extern "C" */