// /* generating a custom N_Vector implementation */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "nvector_octave.h"
#include <sundials/sundials_math.h>
#include <builtin-defun-decls.h>
#include <data.h>

#include <oct.h>
#include <ov.h>
#include <ovl.h>


// /* A vector is a structure with an implementation-dependent
//    'content' field, and a pointer to a structure of vector
//    operations corresponding to that implementation. */

// /* Forward reference for pointer to N_Vector object */
// typedef _SUNDIALS_STRUCT_ _N_VectorCol *N_VectorCol;

// /* Forward reference for pointer to N_Vector_v->Ops object */
// typedef _SUNDIALS_STRUCT_ _N_VectorCol_v->Ops *N_VectorCol_v->Ops;

// struct _N_VectorCol_v->Ops{
//   N_Vector_ID  (*nvgetvectorid)(N_VectorCol);
// };

// struct _N_VectorCol{
//     // _N_VectorContent_Octave *N_VectorContent_Octave;
//     void *content;
//     _N_VectorCol_v->Ops v->ops;
//     SUNContext sunctx;
// };

// in header file

/* Defining the content field */
// struct _N_VectorContent_Octave{
//     int length;  /* vector length       */
//     bool own_data;           /* data ownership flag */
//     double *data;            /* data array, will be our column vector */
// };

// typedef struct _N_VectorContent_Octave *N_VectorContent_Octave;

extern "C"
{

using namespace octave;

/* Defining required MACROS */
#define ZERO RCONST(0.0)
#define HALF RCONST(0.5)
#define ONE RCONST(1.0)
#define ONEPT5 RCONST(1.5)

// need comment
// ColumnVector y;
// ColumnVector& y_data = *y;

  // need comment

  /*
   * -----------------------------------------------------------------
   * Functions exported by nvector_Octave
   * -----------------------------------------------------------------
   */

  /* Private functions for special cases of vector operations */
  static void VCopy_Octave(N_Vector x, N_Vector z);                              /* z=x       */
  static void VSum_Octave(N_Vector x, N_Vector y, N_Vector z);                   /* z=x+y     */
  static void VDiff_Octave(N_Vector x, N_Vector y, N_Vector z);                  /* z=x-y     */
  static void VNeg_Octave(N_Vector x, N_Vector z);                               /* z=-x      */
  static void VScaleSum_Octave(realtype c, N_Vector x, N_Vector y, N_Vector z);  /* z=c(x+y)  */
  static void VScaleDiff_Octave(realtype c, N_Vector x, N_Vector y, N_Vector z); /* z=c(x-y)  */
  static void VLin1_Octave(realtype a, N_Vector x, N_Vector y, N_Vector z);      /* z=ax+y    */
  static void VLin2_Octave(realtype a, N_Vector x, N_Vector y, N_Vector z);      /* z=ax-y    */
  static void Vaxpy_Octave(realtype a, N_Vector x, N_Vector y);                  /* y <- ax+y */
  static void VScaleBy_Octave(realtype a, N_Vector x);                           /* x <- ax   */

  /* Private functions for special cases of vector array operations */
  static int VSumVectorArray_Octave(int nvec, N_Vector *X, N_Vector *Y, N_Vector *Z);                   /* Z=X+Y     */
  static int VDiffVectorArray_Octave(int nvec, N_Vector *X, N_Vector *Y, N_Vector *Z);                  /* Z=X-Y     */
  static int VScaleSumVectorArray_Octave(int nvec, realtype c, N_Vector *X, N_Vector *Y, N_Vector *Z);  /* Z=c(X+Y)  */
  static int VScaleDiffVectorArray_Octave(int nvec, realtype c, N_Vector *X, N_Vector *Y, N_Vector *Z); /* Z=c(X-Y)  */
  static int VLin1VectorArray_Octave(int nvec, realtype a, N_Vector *X, N_Vector *Y, N_Vector *Z);      /* Z=aX+Y    */
  static int VLin2VectorArray_Octave(int nvec, realtype a, N_Vector *X, N_Vector *Y, N_Vector *Z);      /* Z=aX-Y    */
  static int VaxpyVectorArray_Octave(int nvec, realtype a, N_Vector *X, N_Vector *Y);                   /* Y <- aX+Y */

  /*
   * -----------------------------------------------------------------
   * exported functions
   * -----------------------------------------------------------------
   */

  /* ----------------------------------------------------------------
   * Returns vector type ID. Used to identify vector implementation
   * from abstract N_Vector interface.
   */
  N_Vector_ID N_VGetVectorID_Octave(N_Vector v)
  {
    return SUNDIALS_NVEC_SERIAL;
  }

  /* ----------------------------------------------------------------------------
   * Function to create a new empty Col vector
   */

    N_Vector N_VNewEmpty_Octave(SUNContext sunctx)
  {
    N_Vector v;

    /* Create an empty vector object */
    v = NULL;
    v = N_VNewEmpty(sunctx);
    if (v == NULL)
      return (NULL);

    /* Attach operations */

    /* constructors, destructors, and utility operations */
    v->ops->nvgetvectorid     = N_VGetVectorID_Octave;
    v->ops->nvclone           = N_VClone_Octave;
    v->ops->nvcloneempty      = N_VCloneEmpty_Octave;
    v->ops->nvdestroy         = N_VDestroy_Octave;
    v->ops->nvspace           = N_VSpace_Octave;
    v->ops->nvgetarraypointer = N_VGetArrayPointer_Octave;
    v->ops->nvsetarraypointer = NULL;
    v->ops->nvgetlength       = N_VGetLength_Octave;

    /* standard vector operations */
    v->ops->nvlinearsum    = N_VLinearSum_Octave;
    v->ops->nvconst        = N_VConst_Octave;
    v->ops->nvprod         = N_VProd_Octave;
    v->ops->nvdiv          = N_VDiv_Octave;
    v->ops->nvscale        = N_VScale_Octave;
    v->ops->nvabs          = N_VAbs_Octave;
    v->ops->nvinv          = N_VInv_Octave;
    v->ops->nvaddconst     = N_VAddConst_Octave;
    v->ops->nvdotprod      = N_VDotProd_Octave;
    v->ops->nvmaxnorm      = N_VMaxNorm_Octave;
    v->ops->nvwrmsnormmask = N_VWrmsNormMask_Octave;
    v->ops->nvwrmsnorm     = N_VWrmsNorm_Octave;
    v->ops->nvmin          = N_VMin_Octave;
    v->ops->nvwl2norm      = NULL;
    v->ops->nvl1norm       = NULL;
    v->ops->nvcompare      = N_VCompare_Octave;
    v->ops->nvinvtest      = N_VInvTest_Octave;
    v->ops->nvconstrmask   = N_VConstrMask_Octave;
    v->ops->nvminquotient  = N_VMinQuotient_Octave;

    /* fused vector operations (optional) */
    v->ops->nvlinearcombination = NULL;
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_Octave;
    v->ops->nvdotprodmulti = NULL;

    // /* vector array operations (optional) */
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Octave;
    v->ops->nvscalevectorarray = N_VScaleVectorArray_Octave;
    v->ops->nvconstvectorarray = NULL;
    v->ops->nvwrmsnormvectorarray = NULL;
    v->ops->nvwrmsnormmaskvectorarray = NULL;
    v->ops->nvscaleaddmultivectorarray = NULL;
    v->ops->nvlinearcombinationvectorarray = NULL;

    // /*
    //  * OPTIONAL operations with no default implementation.
    //  */

    /* local reduction operations */
    v->ops->nvdotprodlocal     = NULL;
    v->ops->nvmaxnormlocal     = NULL;
    v->ops->nvminlocal         = NULL;
    v->ops->nvl1normlocal      = NULL;
    v->ops->nvinvtestlocal     = NULL;
    v->ops->nvconstrmasklocal  = N_VConstrMask_Octave;
    v->ops->nvminquotientlocal = NULL;
    v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Octave;
    v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Octave;

    // /* single buffer reduction operations */
    // v->ops->nvdotprodmultilocal = NULL;
    // v->ops->nvdotprodmultiallreduce = NULL;

    // /* XBraid interface operations */
    // v->ops->nvbufsize = NULL;
    // v->ops->nvbufpack = NULL;
    // v->ops->nvbufunpack = NULL;

    // /* debugging functions */
    v->ops->nvprint = N_VPrint_Octave;
    v->ops->nvprintfile = NULL;

    return(v);
  }

  N_Vector N_VNew_Octave(int length, SUNContext sunctx)
  {
    N_Vector v;
    void *content;
    // realtype *data;

    if (sunctx == NULL)
      return (NULL);

    /* Create an empty vector object */
    v = NULL;
    v = N_VNewEmpty_Octave(sunctx);
    if (v == NULL)
      return (NULL);

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
    if(NV_DATA_C(v) == NULL) { N_VDestroy_Octave(v); return(NULL); }

    /* Attach data */
    // NV_OWN_DATA_C(v) = SUNTRUE;

  }
  return (v);
}
  /* ----------------------------------------------------------------------------
   * Function to create a serial N_Vector with existing ColumnVector component
   */

  N_Vector N_VMake_Octave(ColumnVector cv, SUNContext sunctx)
  {
    N_Vector v;
    void *content;

    v = NULL;
    v = N_VNewEmpty_Octave(sunctx);
    if (v == NULL)
      return (NULL);

    sunindextype length = cv.numel();
    content = NULL;
    content = new ColumnVector(length);
    content = &cv;
    if (content == NULL)
    {
      N_VDestroy(v);
      return (NULL);
    }
    /* Attach content */
    v->content = content;

    return (v);
  }

// /* ----------------------------------------------------------------------------
//  * Function to create an array of new serial vectors.
//  */

// N_Vector* N_VCloneVectorArray_Octave(int count, N_Vector w)
// {
//   N_Vector* vs;
//   int j;

//   if (count <= 0) return(NULL);

//   vs = NULL;
//   vs = (N_Vector*) malloc(count * sizeof(N_Vector));
//   if(vs == NULL) return(NULL);

//   for (j = 0; j < count; j++) {
//     vs[j] = NULL;
//     vs[j] = N_VClone_Octave(w);
//     if (vs[j] == NULL) {
//       N_VDestroyVectorArray_Octave(vs, j-1);
//       return(NULL);
//     }
//   }

//   return(vs);
// }

// // /* ----------------------------------------------------------------------------
// //  * Function to create an array of new serial vectors with NULL data array.
// //  */

// // N_Vector* N_VCloneVectorArrayEmpty_Octave(int count, N_Vector w)
// // {
// //   N_Vector* vs;
// //   int j;

// //   if (count <= 0) return(NULL);

// //   vs = NULL;
// //   vs = (N_Vector*) malloc(count * sizeof(N_Vector));
// //   if(vs == NULL) return(NULL);

// //   for (j = 0; j < count; j++) {
// //     vs[j] = NULL;
// //     vs[j] = N_VCloneEmpty_Octave(w);
// //     if (vs[j] == NULL) {
// //       N_VDestroyVectorArray_Octave(vs, j-1);
// //       return(NULL);
// //     }
// //   }

// //   return(vs);
// // }

// /* ----------------------------------------------------------------------------
//  * Function to free an array created with N_VCloneVectorArray_Octave
//  */

// void N_VDestroyVectorArray_Octave(N_Vector* vs, int count)
// {
//   int j;

//   for (j = 0; j < count; j++) N_VDestroy_Octave(vs[j]);

//   free(vs); vs = NULL;

//   return;
// }

/* ----------------------------------------------------------------------------
 * Function to return number of vector elements
 */
sunindextype N_VGetLength_Octave(N_Vector v)
{
  return NV_LENGTH_C(v);
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to stdout
 */

// void N_VPrint_Octave(N_Vector x)
// {
//   N_VPrintFile_Octave(x, stdout);
// }

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to outfile
 */

void N_VPrint_Octave(N_Vector x)
{
  sunindextype i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);

  ColumnVector *xv;
  xv = static_cast<ColumnVector *> NV_CONTENT_C(x);

  std::cout<<(*xv);
// #if defined(SUNDIALS_EXTENDED_PRECISION)
//     fprintf(outfile, "%35.32Lg\n", xd[i]);
// #elif defined(SUNDIALS_DOUBLE_PRECISION)
//     std::cout<<(*xv);
// #else
//     fprintf(outfile, "%11.8g\n", xd[i]);
// #endif
  // fprintf(outfile, "\n");

  return;
}

// /*
//  * -----------------------------------------------------------------
//  * implementation of vector operations
//  * -----------------------------------------------------------------
//  */

N_Vector N_VCloneEmpty_Octave(N_Vector w)
{
  N_Vector v;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty_Octave(w->sunctx);
  if (v == NULL) return(NULL);

  /* Attach operations */
  if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

  return(v);
}

N_Vector N_VClone_Octave(N_Vector w)
{
  void* content;
  sunindextype length;

  /* Create vector */
  N_Vector v = N_VCloneEmpty_Octave(w);
  if (v == NULL) return(NULL);

  /* Create content */
  content = NULL;
  content = new ColumnVector(NV_LENGTH_C(w));
  if (content == NULL) { N_VDestroy(v); return(NULL); }

  /* Attach content */
  v->content = content;
  /* Claim data */
  // NV_OWN_DATA_C(v) = SUNTRUE;

  return(v);
}

void N_VDestroy_Octave(N_Vector v)
{
  if (v == NULL) return;

  /* free content */
  if (v->content != NULL) {
    // /* free data array if it's owned by the vector */
    // if (NV_OWN_DATA_C(v) && NV_DATA_C(v) != NULL) {
    //   free(NV_DATA_C(v));
    //   NV_DATA_C(v) = NULL;
    // }
    // free(v->content);
    v->content = NULL;
  }

  /* free ops and vector */
  if (v->ops != NULL) { free(v->ops); v->ops = NULL; }
  free(v); v = NULL;

  return;
}

void N_VSpace_Octave(N_Vector v, sunindextype *lrw, sunindextype *liw)
{
  *lrw = NV_LENGTH_C(v);
  *liw = 1;

  return;
}

realtype *N_VGetArrayPointer_Octave(N_Vector v)
{
  return((realtype *) NV_DATA_C(v));
}

void N_VSetArrayPointer_Octave(realtype* v_data, N_Vector v)
{
  // if (NV_LENGTH_C(v) > 0) NV_DATA_C(v) = v_data;
  

  return;
}

void N_VdrrayPointer_Octave(realtype *v_data, N_Vector v)
{
  if (NV_LENGTH_C(v) > 0) v_data = NV_CONTENT_C(v)->fortran_vec();

  return;
}

/*
 * Linear combination of two vectors: z = a*x + b*y
 */
void N_VLinearSum_Octave(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{

    ColumnVector *xv, *yv, *zv;
    xv = static_cast <ColumnVector *> NV_CONTENT_C(x);
    yv = static_cast <ColumnVector *> NV_CONTENT_C(y);
    zv = static_cast <ColumnVector *> NV_CONTENT_C(z);
  
    ColumnVector *v1, *v2;
    booleantype test;

  // /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE)) {
    (*zv) = (*xv) + (*yv);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    if(test == 0)
    (*zv) = (*yv) - (*xv);
    else
    (*zv) = (*xv) - (*yv);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE)) {
    if(test == 1)
    (*zv) = (*xv) + b * (*yv);
    else
    (*zv) = a * (*xv) + (*yv);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE)) {
    if(test == 1)
    (*zv) = b * (*yv) - (*xv);
    else
    (*zv) = a * (*xv) - (*yv);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b) {
    (*zv) = a * ((*yv) + (*xv));
    return;
  }

  /* Case: a == -b */

  if (a == -b) {
    (*zv) = b * ((*yv) - (*xv));
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  (*zv) = a * (*xv) + b * (*yv);

  return;
}

/*
 * Set all vector elements to a constant: z[i] = c
 */
void N_VConst_Octave(realtype c, N_Vector z)
{
  ColumnVector *zv;
  zv = static_cast<ColumnVector *> NV_CONTENT_C(z);

  zv->fill(c);

  return;
}

/*
 * Elementwise multiply vectors: z[i] = x[i]*y[i]
 */
void N_VProd_Octave(N_Vector x, N_Vector y, N_Vector z)
{

  ColumnVector *xv, *zv, *yv;
  xv = static_cast<ColumnVector *> NV_CONTENT_C(x);
  zv = static_cast<ColumnVector *> NV_CONTENT_C(z);
  yv = static_cast<ColumnVector *> NV_CONTENT_C(y);
  *zv = product((*xv),(*yv));

  return;
}

void N_VDiv_Octave(N_Vector x, N_Vector y, N_Vector z)
{

  ColumnVector *xv, *zv, *yv;
  xv = static_cast<ColumnVector *> NV_CONTENT_C(x);
  zv = static_cast<ColumnVector *> NV_CONTENT_C(z);
  yv = static_cast<ColumnVector *> NV_CONTENT_C(y);
  *zv = quotient((*xv),(*yv));

  return;
}

void N_VScale_Octave(realtype c, N_Vector x, N_Vector z)
{
  ColumnVector *xv = new ColumnVector(NV_LENGTH_C(x));
  ColumnVector *zv = new ColumnVector(NV_LENGTH_C(z));
  xv = static_cast <ColumnVector *> NV_CONTENT_C(x);
  zv = static_cast <ColumnVector *> NV_CONTENT_C(z);

  // checks if both NVectors z and x
  // are stored in the same place
  if (zv == xv) {  /* BLAS usage: scale x <- cx */
    (*xv) *= c;
    // Not using this results in poor iterative
    // algorthm performance
    return;
  }
  ColumnVector temp(NV_LENGTH_C(x));
  temp.fill(1);
  if (c == ONE) {
    // (*zv) += (*xv);
    VCopy_Octave(x, z);  // substituting with (*zv) = (*xv) leads to wrong results.
  } else if (c == -ONE) {
    (*zv) -= (*xv);
  } else {
    (*zv) = c * (*xv);
  }

  return;

  // sunindextype i, N;
  // realtype *xd, *zd;

  // xd = zd = NULL;

  // if (z == x) {  /* BLAS usage: scale x <- cx */
  //   VScaleBy_Octave(c, x);
  //   return;
  // }

  // if (c == ONE) {
  //   VCopy_Octave(x, z);
  // } else if (c == -ONE) {
  //   VNeg_Octave(x, z);
  // } else {
  //   N  = NV_LENGTH_C(x);
  //   xd = NV_DATA_C(x);
  //   zd = NV_DATA_C(z);
  //   for (i = 0; i < N; i++)
  //     zd[i] = c*xd[i];
  // }

  // return;
}

/*
 * Elementwise absolute value: z[i] = |x[i]|
 */
void N_VAbs_Octave(N_Vector x, N_Vector z)
{
  ColumnVector *xv,*zv;
  xv = static_cast <ColumnVector *> NV_CONTENT_C(x);
  zv = static_cast <ColumnVector *> NV_CONTENT_C(z);

  (*zv) = xv->abs();
  return;
}

void N_VInv_Octave(N_Vector x, N_Vector z)
{
  ColumnVector *xv,*zv;
  xv = static_cast <ColumnVector *> NV_CONTENT_C(x);
  zv = static_cast <ColumnVector *> NV_CONTENT_C(z);

  (*zv) = ONE / (*xv);

  return;
}

void N_VAddConst_Octave(N_Vector x, realtype b, N_Vector z)
{
  ColumnVector *xv,*zv;
  xv = static_cast <ColumnVector *> NV_CONTENT_C(x);
  zv = static_cast <ColumnVector *> NV_CONTENT_C(z);

  (*zv) = (*xv) +b;
  return;
}
/*
 * Scalar product of vectors x and y
 */
realtype N_VDotProd_Octave(N_Vector x, N_Vector y)
{
  ColumnVector *xv,*yv;
  realtype sum;
  xv = static_cast <ColumnVector *> NV_CONTENT_C(x);
  yv = static_cast <ColumnVector *> NV_CONTENT_C(y);
  sum = (*xv).transpose() * (*yv);
  return(sum);

}

/*
 * Max norm (L infinity) of vector x
 */
realtype N_VMaxNorm_Octave(N_Vector x)
{
  ColumnVector *xv, abret;
  realtype ret;
  xv = static_cast <ColumnVector *> NV_CONTENT_C(x);
  abret = static_cast <ColumnVector> (xv->abs()); 
  ColumnVector *abp = &(abret);
  ret = abp->max();
  
  return(ret);
}

realtype N_VWrmsNorm_Octave(N_Vector x, N_Vector w)
{
  octave_value_list ov = ovl((N_VWSqrSumLocal_Octave(x, w)/(NV_LENGTH_C(x))));
  return((octave::Fsqrt(ov,1)(0)).double_value());
}

realtype N_VWSqrSumLocal_Octave(N_Vector x, N_Vector w)
{
  N_Vector prod = N_VClone_Octave(x);
  ColumnVector *xv,*yv, *pv;
  realtype sum;
  xv = static_cast <ColumnVector *> NV_CONTENT_C(x);
  yv = static_cast <ColumnVector *> NV_CONTENT_C(w);
  pv = static_cast <ColumnVector *> NV_CONTENT_C(prod);

  (*pv) = product((*xv),(*yv));
  octave_value_list ov = ovl((*pv));
  sum = (octave::Fsumsq(ov,1)(0)).double_value(); 

  return(sum);
}

realtype N_VWrmsNormMask_Octave(N_Vector x, N_Vector w, N_Vector id)
{
  octave_value_list ov = ovl((N_VWSqrSumMaskLocal_Octave(x, w, id)/(NV_LENGTH_C(x))));
  return((octave::Fsqrt(ov,1)(0)).double_value());
}

realtype N_VWSqrSumMaskLocal_Octave(N_Vector x, N_Vector w, N_Vector id)
{
  N_Vector prod = N_VClone_Octave(x);
  N_Vector prod2 = N_VClone_Octave(x);
  ColumnVector *xv,*yv, *mv, *pv, *pv2;
  realtype sum;
  xv = static_cast <ColumnVector *> NV_CONTENT_C(x);
  yv = static_cast <ColumnVector *> NV_CONTENT_C(w);
  mv = static_cast <ColumnVector *> NV_CONTENT_C(id);    //mask vector
  pv = static_cast <ColumnVector *> NV_CONTENT_C(prod);
  pv2 = static_cast <ColumnVector *> NV_CONTENT_C(prod2);

  const octave_value_list ov = octave_value_list({(*mv),(*xv),(*mv)});
  const octave_value_list ov2 = octave_value_list({(*mv),(*yv),(*mv)});
  (*pv) = (octave::Fmerge(ov2,1)(0)).column_vector_value();
  (*pv2) = (octave::Fmerge(ov2,1)(0)).column_vector_value();
  sum = N_VWSqrSumLocal_Octave(prod, prod2);

  return(sum);
}

/* Minimum value in Vector */
realtype N_VMin_Octave(N_Vector x)
{
  ColumnVector *xv;
  realtype min;
  xv = static_cast <ColumnVector *> NV_CONTENT_C(x);

  min = xv->min();
  return(min);
}

// realtype N_VWL2Norm_Octave(N_Vector x, N_Vector w)
// {
//   sunindextype i, N;
//   realtype sum, prodi, *xd, *wd;

//   sum = ZERO;
//   xd = wd = NULL;

//   N  = NV_LENGTH_C(x);
//   xd = NV_DATA_C(x);
//   wd = NV_DATA_C(w);

//   for (i = 0; i < N; i++) {
//     prodi = xd[i]*wd[i];
//     sum += SUNSQR(prodi);
//   }

//   return(SUNRsqrt(sum));
// }

// realtype N_VL1Norm_Octave(N_Vector x)
// {
//   sunindextype i, N;
//   realtype sum, *xd;

//   sum = ZERO;
//   xd = NULL;

//   N  = NV_LENGTH_C(x);
//   xd = NV_DATA_C(x);

//   for (i = 0; i<N; i++)
//     sum += SUNRabs(xd[i]);

//   return(sum);
// }

void N_VCompare_Octave(realtype c, N_Vector x, N_Vector z)
{
//   ColumnVector *xv,*zv;
//   xv = static_cast <ColumnVector *> NV_CONTENT_C(x);
//   zv = static_cast <ColumnVector *> NV_CONTENT_C(z);
// std::cout << "This is a matrix:" << std::endl
// << (*xv) << std::endl;
// octave_value_list in;
// in(0) = *xv;
// in(1) = c;
// octave_value_list out = octave::Fnth_element(in, 1);
// ColumnVector norm_of_the_matrix = out(0).column_vector_value();
// std::cout << "This is the norm of the matrix:" << std::endl
// << norm_of_the_matrix
// << std::endl;
// return;





  // const octave_value_list ov = ovl((*xv),(c));
  // ColumnVector a;
  // // printf("start error");
  // // a = (a >= c).column_vector_value();
  // octave::Fge(ov,1);
  // // (*zv) = a;
  // // ((*xv) >= (*zv));

  // // (*zv) = retval(0).column_vector_value();

  // // return;

  // octave_value_list retval;
  // ColumnVector v ;
  // ColumnVector w;
  // retval(0) = (v >= w);
  // return;
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

booleantype N_VInvTest_Octave(N_Vector x, N_Vector z)
{
  SUNContext ctx;
  sunindextype i, N;
  realtype *xd, *zd;
  booleantype no_zero_found = SUNTRUE ;

  xd = zd = NULL;

  N  = NV_LENGTH_C(x);
  xd = NV_DATA_C(x);
  zd = NV_DATA_C(z);

  N_Vector y = N_VClone(x);

  ColumnVector *xv, *zv, *yv;
  xv = static_cast <ColumnVector *> NV_CONTENT_C(x);
  zv = static_cast <ColumnVector *> NV_CONTENT_C(z);
  yv = static_cast <ColumnVector *> NV_CONTENT_C(y);

  printf("%ld",xv->nnz());
  
  (*zv) = ONE /(*xv);
  if(xv->nnz() != xv->numel())
    no_zero_found = SUNFALSE;
  if(no_zero_found==SUNTRUE)
    (*zv) = ONE /(*xv);
  else {
    (*yv) = (*zv);
    const octave_value_list ov = ovl((*zv));
    (*yv) = (octave::Fisfinite(ov,1)(0)).column_vector_value();
    const octave_value_list ov2 = octave_value_list({(*yv),(*zv),(*yv)});
    (*zv) = (octave::Fmerge(ov2,1)(0)).column_vector_value();
    std::cout<<(*zv)<<" "<<(*yv);
  }

  return no_zero_found;
}

//Should I try writing a lambda function like trilinois
booleantype N_VConstrMask_Octave(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype i, N;
  realtype temp;
  realtype *cd, *xd, *md;
  booleantype test;

  ColumnVector *cv, *xv, *mv;
  xv = static_cast <ColumnVector *> NV_CONTENT_C(x);
  cv = static_cast <ColumnVector *> NV_CONTENT_C(c);
  mv = static_cast <ColumnVector *> NV_CONTENT_C(m);

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

/*
 * Find minimum quotient: minq  = min ( num[i]/denom[i]), denom[i] != 0.
 */
realtype N_VMinQuotient_Octave(N_Vector num, N_Vector denom)
{
    N_Vector z = N_VClone_Octave(num);
    ColumnVector *xv, *zv, *yv;
    realtype min;
    min = BIG_REAL;

    xv = static_cast<ColumnVector *> NV_CONTENT_C(num);
    yv = static_cast<ColumnVector *> NV_CONTENT_C(denom);
    zv = static_cast<ColumnVector *> NV_CONTENT_C(z);
    *zv = quotient((*xv),(*yv));
    min = (zv)->min();

  return(min);
}

// /*
//  * -----------------------------------------------------------------
//  * fused vector operations
//  * -----------------------------------------------------------------
//  */

int N_VLinearCombination_Octave(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  int          i;
  sunindextype j, N;
  realtype*    zd=NULL;
  realtype*    xd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_Octave(c[0], X[0], z);
    return(0);
  }

  /* should have called N_VLinearSum */
  if (nvec == 2) {
    N_VLinearSum_Octave(c[0], X[0], c[1], X[1], z);
    return(0);
  }

  ColumnVector *xv, *zv, *yv;

  xv = static_cast<ColumnVector *> NV_CONTENT_C(X[0]);
  zv = static_cast<ColumnVector *> NV_CONTENT_C(z);

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */

  (*zv) = (*c)*(*xv);
  return(0);
}


int N_VScaleAddMulti_Octave(int nvec, realtype* a, N_Vector x, N_Vector* Y, N_Vector* Z)
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
    N_VLinearSum_Octave(a[0], x, ONE, Y[0], Z[0]);
    return(0);
  }

  ColumnVector *xv, *zv, *yv;

  yv = static_cast<ColumnVector *> NV_CONTENT_C(Y[0]);
  zv = static_cast<ColumnVector *> NV_CONTENT_C(Z[0]);
  xv = static_cast<ColumnVector *> NV_CONTENT_C(x);

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  (*zv) = (*yv) + (*a)*(*xv);
  return(0);
}


// int N_VDotProdMulti_Octave(int nvec, N_Vector x, N_Vector* Y, realtype* dotprods)
// {
//   int          i;
//   sunindextype j, N;
//   realtype*    xd=NULL;
//   realtype*    yd=NULL;

//   /* invalid number of vectors */
//   if (nvec < 1) return(-1);

//   /* should have called N_VDotProd */
//   if (nvec == 1) {
//     dotprods[0] = N_VDotProd_Octave(x, Y[0]);
//     return(0);
//   }

//   /* get vector length and data array */
//   N  = NV_LENGTH_C(x);
//   xd = NV_DATA_C(x);

//   /* compute multiple dot products */
//   for (i=0; i<nvec; i++) {
//     yd = NV_DATA_C(Y[i]);
//     dotprods[i] = ZERO;
//     for (j=0; j<N; j++) {
//       dotprods[i] += xd[j] * yd[j];
//     }
//   }

//   return(0);
// }


/*
 * -----------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------
 */

int N_VLinearSumVectorArray_Octave(int nvec,
                                   realtype a, N_Vector* X,
                                   realtype b, N_Vector* Y,
                                   N_Vector* Z)
{

  ColumnVector *xv, *zv, *yv;
  int          i;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VLinearSum */
  if (nvec == 1) {
    N_VLinearSum_Octave(a, X[0], b, Y[0], Z[0]);
    return(0);
  }                   

  for (i=0; i<nvec; i++) { 
    yv = static_cast<ColumnVector *> NV_CONTENT_C(Y[i]);
    zv = static_cast<ColumnVector *> NV_CONTENT_C(Z[i]);
    xv = static_cast<ColumnVector *> NV_CONTENT_C(X[i]);
    (*zv) = a * (*xv) + b * (*yv);
  }

  return(0);
}


int N_VScaleVectorArray_Octave(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    zd=NULL;
  ColumnVector *xv, *zv;

  xv = static_cast<ColumnVector *> NV_CONTENT_C(X[0]);
  zv = static_cast<ColumnVector *> NV_CONTENT_C(Z[0]);

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_Octave(c[0], X[0], Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_LENGTH_C(Z[0]);

  /*
   * X[i] *= c[i]
   */
  // only checking if NVectors X and Z 
  // are stored in the same address
  if (xv == zv) { 
    // (*xv) *= (*c); //this should work but doesn't, in both N_VScale
    //                // instead results in convergence failures and wrong answer
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
  (*zv) = (*c) * (*xv); 
  return(0);
}


// int N_VConstVectorArray_Octave(int nvec, realtype c, N_Vector* Z)
// {
//   int          i;
//   sunindextype j, N;
//   realtype*    zd=NULL;

//   /* invalid number of vectors */
//   if (nvec < 1) return(-1);

//   /* should have called N_VConst */
//   if (nvec == 1) {
//     N_VConst_Octave(c, Z[0]);
//     return(0);
//   }

//   /* get vector length */
//   N = NV_LENGTH_C(Z[0]);

//   /* set each vector in the vector array to a constant */
//   for (i=0; i<nvec; i++) {
//     zd = NV_DATA_C(Z[i]);
//     for (j=0; j<N; j++) {
//       zd[j] = c;
//     }
//   }

//   return(0);
// }


// int N_VWrmsNormVectorArray_Octave(int nvec, N_Vector* X, N_Vector* W,
//                                   realtype* nrm)
// {
//   int          i;
//   sunindextype j, N;
//   realtype*    wd=NULL;
//   realtype*    xd=NULL;

//   /* invalid number of vectors */
//   if (nvec < 1) return(-1);

//   /* should have called N_VWrmsNorm */
//   if (nvec == 1) {
//     nrm[0] = N_VWrmsNorm_Octave(X[0], W[0]);
//     return(0);
//   }

//   /* get vector length */
//   N = NV_LENGTH_C(X[0]);

//   /* compute the WRMS norm for each vector in the vector array */
//   for (i=0; i<nvec; i++) {
//     xd = NV_DATA_C(X[i]);
//     wd = NV_DATA_C(W[i]);
//     nrm[i] = ZERO;
//     for (j=0; j<N; j++) {
//       nrm[i] += SUNSQR(xd[j] * wd[j]);
//     }
//     nrm[i] = SUNRsqrt(nrm[i]/N);
//   }

//   return(0);
// }


// int N_VWrmsNormMaskVectorArray_Octave(int nvec, N_Vector* X, N_Vector* W,
//                                       N_Vector id, realtype* nrm)
// {
//   int          i;
//   sunindextype j, N;
//   realtype*    wd=NULL;
//   realtype*    xd=NULL;
//   realtype*    idd=NULL;

//   /* invalid number of vectors */
//   if (nvec < 1) return(-1);

//   /* should have called N_VWrmsNorm */
//   if (nvec == 1) {
//     nrm[0] = N_VWrmsNormMask_Octave(X[0], W[0], id);
//     return(0);
//   }

//   /* get vector length and mask data array */
//   N   = NV_LENGTH_C(X[0]);
//   idd = NV_DATA_C(id);

//   /* compute the WRMS norm for each vector in the vector array */
//   for (i=0; i<nvec; i++) {
//     xd = NV_DATA_C(X[i]);
//     wd = NV_DATA_C(W[i]);
//     nrm[i] = ZERO;
//     for (j=0; j<N; j++) {
//       if (idd[j] > ZERO)
//         nrm[i] += SUNSQR(xd[j] * wd[j]);
//     }
//     nrm[i] = SUNRsqrt(nrm[i]/N);
//   }

//   return(0);
// }


// int N_VScaleAddMultiVectorArray_Octave(int nvec, int nsum, realtype* a,
//                                         N_Vector* X, N_Vector** Y, N_Vector** Z)
// {
//   int          i, j;
//   sunindextype k, N;
//   realtype*    xd=NULL;
//   realtype*    yd=NULL;
//   realtype*    zd=NULL;

//   int          retval;
//   N_Vector*   YY;
//   N_Vector*   ZZ;

//   /* invalid number of vectors */
//   if (nvec < 1) return(-1);
//   if (nsum < 1) return(-1);

//   /* ---------------------------
//    * Special cases for nvec == 1
//    * --------------------------- */

//   if (nvec == 1) {

//     /* should have called N_VLinearSum */
//     if (nsum == 1) {
//       N_VLinearSum_Octave(a[0], X[0], ONE, Y[0][0], Z[0][0]);
//       return(0);
//     }

//     /* should have called N_VScaleAddMulti */
//     YY = (N_Vector*) malloc(nsum * sizeof(N_Vector));
//     ZZ = (N_Vector*) malloc(nsum * sizeof(N_Vector));

//     for (j=0; j<nsum; j++) {
//       YY[j] = Y[j][0];
//       ZZ[j] = Z[j][0];
//     }

//     retval = N_VScaleAddMulti_Octave(nsum, a, X[0], YY, ZZ);

//     free(YY);
//     free(ZZ);
//     return(retval);
//   }

//   /* --------------------------
//    * Special cases for nvec > 1
//    * -------------------------- */

//   /* should have called N_VLinearSumVectorArray */
//   if (nsum == 1) {
//     retval = N_VLinearSumVectorArray_Octave(nvec, a[0], X, ONE, Y[0], Z[0]);
//     return(retval);
//   }

//   /* ----------------------------
//    * Compute multiple linear sums
//    * ---------------------------- */

//   /* get vector length */
//   N = NV_LENGTH_C(X[0]);

//   /*
//    * Y[i][j] += a[i] * x[j]
//    */
//   if (Y == Z) {
//     for (i=0; i<nvec; i++) {
//       xd = NV_DATA_C(X[i]);
//       for (j=0; j<nsum; j++){
//         yd = NV_DATA_C(Y[j][i]);
//         for (k=0; k<N; k++) {
//           yd[k] += a[j] * xd[k];
//         }
//       }
//     }
//     return(0);
//   }

//   /*
//    * Z[i][j] = Y[i][j] + a[i] * x[j]
//    */
//   for (i=0; i<nvec; i++) {
//     xd = NV_DATA_C(X[i]);
//     for (j=0; j<nsum; j++){
//       yd = NV_DATA_C(Y[j][i]);
//       zd = NV_DATA_C(Z[j][i]);
//       for (k=0; k<N; k++) {
//         zd[k] = a[j] * xd[k] + yd[k];
//       }
//     }
//   }
//   return(0);
// }


// int N_VLinearCombinationVectorArray_Octave(int nvec, int nsum, realtype* c,
//                                            N_Vector** X, N_Vector* Z)
// {
//   int          i; /* vector arrays index in summation [0,nsum) */
//   int          j; /* vector index in vector array     [0,nvec) */
//   sunindextype k; /* element index in vector          [0,N)    */
//   sunindextype N;
//   realtype*    zd=NULL;
//   realtype*    xd=NULL;

//   int          retval;
//   realtype*    ctmp;
//   N_Vector*   Y;

//   /* invalid number of vectors */
//   if (nvec < 1) return(-1);
//   if (nsum < 1) return(-1);

//   /* ---------------------------
//    * Special cases for nvec == 1
//    * --------------------------- */

//   if (nvec == 1) {

//     /* should have called N_VScale */
//     if (nsum == 1) {
//       N_VScale_Octave(c[0], X[0][0], Z[0]);
//       return(0);
//     }

//     /* should have called N_VLinearSum */
//     if (nsum == 2) {
//       N_VLinearSum_Octave(c[0], X[0][0], c[1], X[1][0], Z[0]);
//       return(0);
//     }

//     /* should have called N_VLinearCombination */
//     Y = (N_Vector*) malloc(nsum * sizeof(N_Vector));

//     for (i=0; i<nsum; i++) {
//       Y[i] = X[i][0];
//     }

//     retval = N_VLinearCombination_Octave(nsum, c, Y, Z[0]);

//     free(Y);
//     return(retval);
//   }

//   /* --------------------------
//    * Special cases for nvec > 1
//    * -------------------------- */

//   /* should have called N_VScaleVectorArray */
//   if (nsum == 1) {

//     ctmp = (realtype*) malloc(nvec * sizeof(realtype));

//     for (j=0; j<nvec; j++) {
//       ctmp[j] = c[0];
//     }

//     retval = N_VScaleVectorArray_Octave(nvec, ctmp, X[0], Z);

//     free(ctmp);
//     return(retval);
//   }

//   /* should have called N_VLinearSumVectorArray */
//   if (nsum == 2) {
//     retval = N_VLinearSumVectorArray_Octave(nvec, c[0], X[0], c[1], X[1], Z);
//     return(retval);
//   }

//   /* --------------------------
//    * Compute linear combination
//    * -------------------------- */

//   /* get vector length */
//   N = NV_LENGTH_C(Z[0]);

//   /*
//    * X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
//    */
//   if ((X[0] == Z) && (c[0] == ONE)) {
//     for (j=0; j<nvec; j++) {
//       zd = NV_DATA_C(Z[j]);
//       for (i=1; i<nsum; i++) {
//         xd = NV_DATA_C(X[i][j]);
//         for (k=0; k<N; k++) {
//           zd[k] += c[i] * xd[k];
//         }
//       }
//     }
//     return(0);
//   }

//   /*
//    * X[0][j] = c[0] * X[0][j] + sum{ c[i] * X[i][j] }, i = 1,...,nvec-1
//    */
//   if (X[0] == Z) {
//     for (j=0; j<nvec; j++) {
//       zd = NV_DATA_C(Z[j]);
//       for (k=0; k<N; k++) {
//         zd[k] *= c[0];
//       }
//       for (i=1; i<nsum; i++) {
//         xd = NV_DATA_C(X[i][j]);
//         for (k=0; k<N; k++) {
//           zd[k] += c[i] * xd[k];
//         }
//       }
//     }
//     return(0);
//   }

//   /*
//    * Z[j] = sum{ c[i] * X[i][j] }, i = 0,...,nvec-1
//    */
//   for (j=0; j<nvec; j++) {
//     xd = NV_DATA_C(X[0][j]);
//     zd = NV_DATA_C(Z[j]);
//     for (k=0; k<N; k++) {
//       zd[k] = c[0] * xd[k];
//     }
//     for (i=1; i<nsum; i++) {
//       xd = NV_DATA_C(X[i][j]);
//       for (k=0; k<N; k++) {
//         zd[k] += c[i] * xd[k];
//       }
//     }
//   }
//   return(0);
// }


// /*
//  * -----------------------------------------------------------------
//  * private functions for special cases of vector operations
//  * -----------------------------------------------------------------
//  */

static void VScaleBy_Octave(realtype a, N_Vector x)
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

static void VCopy_Octave(N_Vector x, N_Vector z) 
// used in N_VScale and if just do (*zv)=(*xv) results in wrong answer.
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

static void VSum_Octave(N_Vector x, N_Vector y, N_Vector z)
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

static void VDiff_Octave(N_Vector x, N_Vector y, N_Vector z)
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

static void VNeg_Octave(N_Vector x, N_Vector z)
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

static void VScaleSum_Octave(realtype c, N_Vector x, N_Vector y, N_Vector z)
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

static void VScaleDiff_Octave(realtype c, N_Vector x, N_Vector y, N_Vector z)
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

static void VLin1_Octave(realtype a, N_Vector x, N_Vector y, N_Vector z)
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

static void VLin2_Octave(realtype a, N_Vector x, N_Vector y, N_Vector z)
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

// static void Vaxpy_Octave(realtype a, N_Vector x, N_Vector y)
// {
//   sunindextype i, N;
//   realtype *xd, *yd;

//   xd = yd = NULL;

//   N  = NV_LENGTH_C(x);
//   xd = NV_DATA_C(x);
//   yd = NV_DATA_C(y);

//   if (a == ONE) {
//     for (i = 0; i < N; i++)
//       yd[i] += xd[i];
//     return;
//   }

//   if (a == -ONE) {
//     for (i = 0; i < N; i++)
//       yd[i] -= xd[i];
//     return;
//   }

//   for (i = 0; i < N; i++)
//     yd[i] += a*xd[i];

//   return;
// }

/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector array operations
 * -----------------------------------------------------------------
 */

static int VSumVectorArray_Octave(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z)
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

static int VDiffVectorArray_Octave(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z)
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

static int VScaleSumVectorArray_Octave(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z)
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

static int VScaleDiffVectorArray_Octave(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z)
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

static int VLin1VectorArray_Octave(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z)
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

static int VLin2VectorArray_Octave(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z)
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

static int VaxpyVectorArray_Octave(int nvec, realtype a, N_Vector* X, N_Vector* Y)
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

int N_VEnableFusedOps_Octave(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

   if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = NULL;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Octave;
    v->ops->nvdotprodmulti      = NULL;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Octave;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_Octave;
    v->ops->nvconstvectorarray             = NULL;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = NULL;
    v->ops->nvlinearcombinationvectorarray = NULL;
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


// int N_VEnableLinearCombination_Octave(N_Vector v, booleantype tf)
// {
//   /* check that vector is non-NULL */
//   if (v == NULL) return(-1);

//   /* check that ops structure is non-NULL */
//   if (v->ops == NULL) return(-1);

//   /* enable/disable operation */
//   if (tf)
//     v->ops->nvlinearcombination = N_VLinearCombination_Octave;
//   else
//     v->ops->nvlinearcombination = NULL;

//   /* return success */
//   return(0);
// }

// int N_VEnableScaleAddMulti_Octave(N_Vector v, booleantype tf)
// {
//   /* check that vector is non-NULL */
//   if (v == NULL) return(-1);

//   /* check that ops structure is non-NULL */
//   if (v->ops == NULL) return(-1);

//   /* enable/disable operation */
//   if (tf)
//     v->ops->nvscaleaddmulti = N_VScaleAddMulti_Octave;
//   else
//     v->ops->nvscaleaddmulti = NULL;

//   /* return success */
//   return(0);
// }

// int N_VEnableDotProdMulti_Octave(N_Vector v, booleantype tf)
// {
//   /* check that vector is non-NULL */
//   if (v == NULL) return(-1);

//   /* check that ops structure is non-NULL */
//   if (v->ops == NULL) return(-1);

//   /* enable/disable operation */
//   if (tf)
//     v->ops->nvdotprodmulti = N_VDotProdMulti_Octave;
//   else
//     v->ops->nvdotprodmulti = NULL;

//   /* return success */
//   return(0);
// }

// int N_VEnableLinearSumVectorArray_Octave(N_Vector v, booleantype tf)
// {
//   /* check that vector is non-NULL */
//   if (v == NULL) return(-1);

//   /* check that ops structure is non-NULL */
//   if (v->ops == NULL) return(-1);

//   /* enable/disable operation */
//   if (tf)
//     v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Octave;
//   else
//     v->ops->nvlinearsumvectorarray = NULL;

//   /* return success */
//   return(0);
// }

// int N_VEnableScaleVectorArray_Octave(N_Vector v, booleantype tf)
// {
//   /* check that vector is non-NULL */
//   if (v == NULL) return(-1);

//   /* check that ops structure is non-NULL */
//   if (v->ops == NULL) return(-1);

//   /* enable/disable operation */
//   if (tf)
//     v->ops->nvscalevectorarray = N_VScaleVectorArray_Octave;
//   else
//     v->ops->nvscalevectorarray = NULL;

//   /* return success */
//   return(0);
// }

// int N_VEnableConstVectorArray_Octave(N_Vector v, booleantype tf)
// {
//   /* check that vector is non-NULL */
//   if (v == NULL) return(-1);

//   /* check that ops structure is non-NULL */
//   if (v->ops == NULL) return(-1);

//   /* enable/disable operation */
//   if (tf)
//     v->ops->nvconstvectorarray = N_VConstVectorArray_Octave;
//   else
//     v->ops->nvconstvectorarray = NULL;

//   /* return success */
//   return(0);
// }

// int N_VEnableWrmsNormVectorArray_Octave(N_Vector v, booleantype tf)
// {
//   /* check that vector is non-NULL */
//   if (v == NULL) return(-1);

//   /* check that ops structure is non-NULL */
//   if (v->ops == NULL) return(-1);

//   /* enable/disable operation */
//   if (tf)
//     v->ops->nvwrmsnormvectorarray = N_VWrmsNormVectorArray_Octave;
//   else
//     v->ops->nvwrmsnormvectorarray = NULL;

//   /* return success */
//   return(0);
// }

// int N_VEnableWrmsNormMaskVectorArray_Octave(N_Vector v, booleantype tf)
// {
//   /* check that vector is non-NULL */
//   if (v == NULL) return(-1);

//   /* check that ops structure is non-NULL */
//   if (v->ops == NULL) return(-1);

//   /* enable/disable operation */
//   if (tf)
//     v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_Octave;
//   else
//     v->ops->nvwrmsnormmaskvectorarray = NULL;

//   /* return success */
//   return(0);
// }

// int N_VEnableScaleAddMultiVectorArray_Octave(N_Vector v, booleantype tf)
// {
//   /* check that vector is non-NULL */
//   if (v == NULL) return(-1);

//   /* check that ops structure is non-NULL */
//   if (v->ops == NULL) return(-1);

//   /* enable/disable operation */
//   if (tf)
//     v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Octave;
//   else
//     v->ops->nvscaleaddmultivectorarray = NULL;

//   /* return success */
//   return(0);
// }

// int N_VEnableLinearCombinationVectorArray_Octave(N_Vector v, booleantype tf)
// {
//   /* check that vector is non-NULL */
//   if (v == NULL) return(-1);

//   /* check that ops structure is non-NULL */
//   if (v->ops == NULL) return(-1);

//   /* enable/disable operation */
//   if (tf)
//     v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Octave;
//   else
//     v->ops->nvlinearcombinationvectorarray = NULL;

//   /* return success */
//   return(0);
// }

} /* extern "C" */