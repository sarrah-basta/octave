/* generating a custom N_Vector implementation */
#include <stdio.h>
#include <stdlib.h>
// #include <stddlib.h>
#include <vector>
#include <nvector/nvector_col.h>

// #include "dColVector.h"
// #include "dMatrix.h"
// #include "dSparse.h"
// #include "f77-fcn.h"
// #include "lo-utils.h"

// #include "Cell.h"
// #include "defun-dld.h"
// #include "error.h"
// #include "errwarn.h"
// #include "oct-map.h"
// #include "ov.h"
// #include "ovl.h"
// #include "pager.h"
// #include "parse.h"

/* for custom implementation of N_Vector using dColVector */
// typedef struct N_Vector *NColVector;

// #ifdef __cplusplus  /* wrapper to enable C++ usage */
// extern "C" {
// #endif

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
//     // _N_VectorContent_Col *N_VectorContent_Col;
//     void *content;
//     _N_VectorCol_v->Ops v->ops;
//     SUNContext sunctx;
// };

//in header file 

/* Defining the content field */
// struct _N_VectorContent_Col{
//     int length;  /* vector length       */
//     bool own_data;           /* data ownership flag */
//     double *data;            /* data array, will be our column vector */
// };

// typedef struct _N_VectorContent_Col *N_VectorContent_Col;

extern "C" {

using namespace sundials;

/*
 * -----------------------------------------------------------------
 * Macros NV_CONTENT_S, NV_DATA_S, NV_OWN_DATA_S,
 *        NV_LENGTH_S, and NV_Ith_S
 * -----------------------------------------------------------------
 */
#define NV_CONTENT_C(v)  ( (N_VectorContent_Col)(v->content) )

#define NV_LENGTH_C(v)   ( NV_CONTENT_C(v)->length )

#define NV_OWN_DATA_C(v) ( NV_CONTENT_C(v)->own_data )

#define NV_DATA_C(v)     ( NV_CONTENT_C(v)->data )

#define NV_Ith_C(v,i)    ( NV_DATA_C(v)[i] )
/*
 * -----------------------------------------------------------------
 * Functions exported by nvector_Col
 * -----------------------------------------------------------------
 */

/* Defining required MACROS */
#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
// N_Vector_ID N_VGetVectorID_Col(N_Vector v)
// {
//   return SUNDIALS_NVEC_CUSTOM;
// }

/* ----------------------------------------------------------------------------
 * Function to create a new empty Col vector
 */

N_Vector N_VNewEmpty_Col(int length, SUNContext sunctx)
{
    N_Vector v;
    N_VectorContent_Col content;
    // N_VectorContent_Col content;

    if (sunctx == NULL) return(NULL);

    /* Create an empty vector object */
    v = NULL;
    v = N_VNewEmpty(sunctx);
    if (v == NULL) return(NULL);

    /* Attach operations */

    /* initialize operations to NULL */

  /*
   * REQUIRED operations.
   *
   * These must be implemented by derivations of the generic N_Vector.
   */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid           = N_VGetVectorID_Col;
  v->ops->nvclone                 = NULL;
  v->ops->nvcloneempty            = NULL;
  v->ops->nvdestroy               = NULL;
  v->ops->nvspace                 = NULL;
  v->ops->nvgetarraypointer       = NULL;
  v->ops->nvgetdevicearraypointer = NULL;
  v->ops->nvsetarraypointer       = NULL;
  v->ops->nvgetcommunicator       = NULL;
  v->ops->nvgetlength             = NULL;

  /* standard vector operations */
  v->ops->nvlinearsum    = NULL;
  v->ops->nvconst        = NULL;
  v->ops->nvprod         = NULL;
  v->ops->nvdiv          = NULL;
  v->ops->nvscale        = NULL;
  v->ops->nvabs          = NULL;
  v->ops->nvinv          = NULL;
  v->ops->nvaddconst     = NULL;
  v->ops->nvdotprod      = NULL;
  v->ops->nvmaxnorm      = NULL;
  v->ops->nvwrmsnorm     = NULL;
  v->ops->nvwrmsnormmask = NULL;
  v->ops->nvmin          = NULL;
  v->ops->nvwl2norm      = NULL;
  v->ops->nvl1norm       = NULL;
  v->ops->nvcompare      = NULL;
  v->ops->nvinvtest      = NULL;
  v->ops->nvconstrmask   = NULL;
  v->ops->nvminquotient  = NULL;

  /*
   * OPTIONAL operations.
   *
   * These operations provide default implementations that may be overriden.
   */

  /* fused vector operations (optional) */
  v->ops->nvlinearcombination = NULL;
  v->ops->nvscaleaddmulti     = NULL;
  v->ops->nvdotprodmulti      = NULL;

  /* vector array operations (optional) */
  v->ops->nvlinearsumvectorarray         = NULL;
  v->ops->nvscalevectorarray             = NULL;
  v->ops->nvconstvectorarray             = NULL;
  v->ops->nvwrmsnormvectorarray          = NULL;
  v->ops->nvwrmsnormmaskvectorarray      = NULL;
  v->ops->nvscaleaddmultivectorarray     = NULL;
  v->ops->nvlinearcombinationvectorarray = NULL;

  /*
   * OPTIONAL operations with no default implementation.
   */

  /* local reduction operations (optional) */
  v->ops->nvdotprodlocal      = NULL;
  v->ops->nvmaxnormlocal      = NULL;
  v->ops->nvminlocal          = NULL;
  v->ops->nvl1normlocal       = NULL;
  v->ops->nvinvtestlocal      = NULL;
  v->ops->nvconstrmasklocal   = NULL;
  v->ops->nvminquotientlocal  = NULL;
  v->ops->nvwsqrsumlocal      = NULL;
  v->ops->nvwsqrsummasklocal  = NULL;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal     = NULL;
  v->ops->nvdotprodmultiallreduce = NULL;

  /* XBraid interface operations */
  v->ops->nvbufsize   = NULL;
  v->ops->nvbufpack   = NULL;
  v->ops->nvbufunpack = NULL;

  /* debugging functions */
  v->ops->nvprint     = NULL;
  v->ops->nvprintfile = NULL;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Col) malloc(sizeof *content);
  if (content == NULL) { N_VDestroy(v); return(NULL); }

  /* Attach content */
  v->content = content;

  /* initialize content and sunctx to NULL */
  content->length = length;
  content->own_data = SUNFALSE;
  content->data     = NULL;

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a new custom implementation column vector
 */
// N_Vector N_VNew_Col(int length)
// {
//   N_Vector v;
//   double *data;

//   v = NULL;
//   v = N_VNewEmpty_Col(length);
//   if (v == NULL) return(NULL);

//   /* Create data */
//   if (length > 0) {

//     /* Allocate memory */
//     data = NULL;
//     data = (realtype *) malloc(length * sizeof(realtype));
//     if(data == NULL) { N_VDestroy_Serial(v); return(NULL); }

//     /* Attach data */
//     NV_OWN_DATA_S(v) = SUNTRUE;
//     NV_DATA_S(v)     = data;

//   }

//   return(v);
// }

/* ----------------------------------------------------------------------------
 * Function to create a serial N_Vector with user data component
 */

N_Vector N_VMake_Col(sunindextype length, double *v_data, SUNContext sunctx)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Col(length, sunctx);
  if (v == NULL) return(NULL);

  if (length > 0) {
    /* Attach data */
    NV_OWN_DATA_C(v) = SUNFALSE;
    NV_DATA_C(v)     = v_data;
  }

  return(v);
}


} /* extern "C" */