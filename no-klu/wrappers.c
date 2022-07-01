/* generating a custom N_Vector implementation */
#include <stdio.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_math.h>

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

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Defining the content field */
struct _N_VectorContent_Col{
    int length;  /* vector length       */
    bool own_data;           /* data ownership flag */
    double *data;            /* data array          */
};

typedef struct _N_VectorContent_Col *N_VectorContent_Col;

/*
 * -----------------------------------------------------------------
 * Macros NV_CONTENT_S, NV_DATA_S, NV_OWN_DATA_S,
 *        NV_LENGTH_S, and NV_Ith_S
 * -----------------------------------------------------------------
 */

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
N_Vector_ID N_VGetVectorID_Col(N_Vector v)
{
  return SUNDIALS_NVEC_CUSTOM;
}

/* ----------------------------------------------------------------------------
 * Function to create a new empty Col vector
 */

N_Vector N_VNewEmpty_Col(int length)
{
    N_Vector v;
    N_VectorContent_Col content;

    /* Create an empty vector object */
    v = NULL;
    v = N_VNewEmpty();
    if (v == NULL) return(NULL);

    /* Attach operations */

    /* constructors, destructors, and utility operations */
    v->ops->nvgetvectorid     = N_VGetVectorID_Col;
    v->ops->nvclone           = N_VClone_Col;
    v->ops->nvcloneempty      = N_VCloneEmpty_Col;
    v->ops->nvdestroy         = N_VDestroy_Col;
    v->ops->nvspace           = N_VSpace_Col;
    v->ops->nvgetarraypointer = N_VGetArrayPointer_Col;
    // v->ops->nvsetarraypointer = N_VSetArrayPointer_Col;
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
    // v->ops->nvdotprod      = N_VDotProd_Col;
    v->ops->nvmaxnorm      = N_VMaxNorm_Col;
    v->ops->nvwrmsnormmask = N_VWrmsNormMask_Col;
    v->ops->nvwrmsnorm     = N_VWrmsNorm_Col;
    v->ops->nvmin          = N_VMin_Col;
    // v->ops->nvwl2norm      = N_VWL2Norm_Col;
    // v->ops->nvl1norm       = N_VL1Norm_Col;
    v->ops->nvcompare      = N_VCompare_Col;
    // v->ops->nvinvtest      = N_VInvTest_Col;
    v->ops->nvconstrmask   = N_VConstrMask_Col;
    v->ops->nvminquotient  = N_VMinQuotient_Col;

    /* Create content */
    content = NULL;
    content = (N_VectorContent_Col) malloc(sizeof *content);
    if (content == NULL) { N_VDestroy(v); return(NULL); }

    /* Attach content */
    v->content = content;

    /* Initialize content */
    content->length   = length;
    content->own_data = SUNFALSE;
    content->data     = NULL;

    return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a new custom implementation column vector
 */
N_Vector N_VNew_Col(int length)
{
  N_Vector v;
  double *data;

  v = NULL;
  v = N_VNewEmpty_Col(length);
  if (v == NULL) return(NULL);

  /* Create data */
  if (length > 0) {

    /* Allocate memory */
    data = NULL;
    data = (realtype *) malloc(length * sizeof(realtype));
    if(data == NULL) { N_VDestroy_Serial(v); return(NULL); }

    /* Attach data */
    NV_OWN_DATA_S(v) = SUNTRUE;
    NV_DATA_S(v)     = data;

  }

  return(v);
}


}