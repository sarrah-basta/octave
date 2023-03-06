#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include "oct-sundials.h"

OCTAVE_BEGIN_NAMESPACE(octave)

/* --------------------------------------------------------------------------
 * Utility functions
 * -------------------------------------------------------------------------- */
/* Private functions for special cases of vector operations */
static void VCopy_Local (N_Vector x, N_Vector z);

/* --------------------------------------------------------------------------
 * Constructors
 * -------------------------------------------------------------------------- */

/* Returns vector type ID. Used to identify vector implementation
  from abstract N_Vector interface. */
N_Vector_ID N_VGetVectorID ([[maybe_unused]] N_Vector v)
{
  return SUNDIALS_NVEC_CUSTOM;
}

/* Function to create a new empty octave vector */

N_Vector N_VNewEmpty (SP_ARG_SUNCONTEXT)
{
  N_Vector v;

  /* Create an empty vector object */
  v = nullptr;
  v = ::N_VNewEmpty (SP_OCTAVE_SUNCONTEXT);
  if (v == nullptr)
    return nullptr;

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = octave::N_VGetVectorID;
  v->ops->nvclone           = octave::N_VClone;
  v->ops->nvcloneempty      = octave::N_VCloneEmpty;
  v->ops->nvdestroy         = octave::N_VDestroy;
  v->ops->nvspace           = octave::N_VSpace;
  v->ops->nvgetarraypointer = octave::N_VGetArrayPointer;
  v->ops->nvsetarraypointer = nullptr;
  v->ops->nvgetlength       = octave::N_VGetLength;

  /* standard vector operations */
  v->ops->nvlinearsum    = octave::N_VLinearSum;
  v->ops->nvconst        = octave::N_VConst;
  v->ops->nvprod         = octave::N_VProd;
  v->ops->nvdiv          = octave::N_VDiv;
  v->ops->nvscale        = octave::N_VScale;
  v->ops->nvabs          = octave::N_VAbs;
  v->ops->nvinv          = octave::N_VInv;
  v->ops->nvaddconst     = octave::N_VAddConst;
  v->ops->nvdotprod      = octave::N_VDotProd;
  v->ops->nvmaxnorm      = octave::N_VMaxNorm;
  v->ops->nvwrmsnormmask = octave::N_VWrmsNormMask;
  v->ops->nvwrmsnorm     = octave::N_VWrmsNorm;
  v->ops->nvmin          = octave::N_VMin;
  v->ops->nvwl2norm      = nullptr;
  v->ops->nvl1norm       = nullptr;
  v->ops->nvcompare      = octave::N_VCompare;
  v->ops->nvinvtest      = octave::N_VInvTest;
  v->ops->nvconstrmask   = octave::N_VConstrMask;
  v->ops->nvminquotient  = octave::N_VMinQuotient;

  /* fused vector operations (optional) */
  v->ops->nvlinearcombination = nullptr;
  v->ops->nvscaleaddmulti = octave::N_VScaleAddMulti;
  v->ops->nvdotprodmulti = nullptr;

  /* vector array operations (optional) */
  v->ops->nvlinearsumvectorarray = octave::N_VLinearSumVectorArray;
  v->ops->nvscalevectorarray = octave::N_VScaleVectorArray;
  v->ops->nvconstvectorarray = nullptr;
  v->ops->nvwrmsnormvectorarray = nullptr;
  v->ops->nvwrmsnormmaskvectorarray = nullptr;
  v->ops->nvscaleaddmultivectorarray = nullptr;
  v->ops->nvlinearcombinationvectorarray = nullptr;

  /*
   * OPTIONAL operations with no default implementation.
   */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = nullptr;
  v->ops->nvmaxnormlocal     = nullptr;
  v->ops->nvminlocal         = nullptr;
  v->ops->nvl1normlocal      = nullptr;
  v->ops->nvinvtestlocal     = nullptr;
  v->ops->nvconstrmasklocal  = nullptr;
  v->ops->nvminquotientlocal = nullptr;
  v->ops->nvwsqrsumlocal     = octave::N_VWSqrSumLocal;
  v->ops->nvwsqrsummasklocal = octave::N_VWSqrSumMaskLocal;

  /* debugging functions */
  v->ops->nvprint = octave::N_VPrint;
  v->ops->nvprintfile = nullptr;

  return v;
}

/* constructor to create a new nvector implementation
  with a required length */
N_Vector N_VNew (int length ARG_SUNCONTEXT)
{
  N_Vector v;
  void *content;
  #  if defined (HAVE_SUNDIALS_SUNCONTEXT)
  if (m_sunContext == nullptr)
    return nullptr;
  #endif

  /* Create an empty vector object */
  v = nullptr;
  v = octave::N_VNewEmpty (SP_OCTAVE_SUNCONTEXT);
  if (v == nullptr)
    return nullptr;

  /* Create data */
  if (length > 0) {

  /* Allocate memory & Create content */
  content = nullptr;
  content = new ColumnVector (length);
  if (content == nullptr)
    {
      ::N_VDestroy(v);
      return nullptr;
    }
  /* Attach content */
  v->content = content;
  if(NV_DATA_C (v) == nullptr)
    {
      octave::N_VDestroy (v);
      return nullptr;
    }

}
return v;
}
/* Constructor to create a serial N_Vector
  from existing ColumnVector component */

N_Vector N_VMake (const ColumnVector& cv ARG_SUNCONTEXT)
{
  N_Vector v;
  void *content;

  /* Create an empty vector object */
  v = nullptr;
  v = octave::N_VNewEmpty (SP_OCTAVE_SUNCONTEXT);
  if (v == nullptr)
    return nullptr;

  /* Create content to point to allocated memory */
  content = nullptr;
  content = const_cast<ColumnVector *> (&cv);

  if (content == nullptr)
    {
      ::N_VDestroy (v);
      return nullptr;
    }

  /* Attach content */
  v->content = content;
  return v;
}

/* --------------------------------------------------------------------------
 * Vector get, set, and utility functions
 * -------------------------------------------------------------------------- */

/* Function to return the global length of the vector.*/
sunindextype N_VGetLength (N_Vector v)
{
  return NV_LENGTH_C (v);
}

/* Function to print the vector vector to stdout */
void N_VPrint (N_Vector x)
{
  ColumnVector *xv;
  xv = static_cast<ColumnVector *> NV_CONTENT_C (x);

  std::cout<<(*xv);

  return;
}

/* --------------------------------------------------------------------------
 * Vector operations implementations
 * -------------------------------------------------------------------------- */


N_Vector N_VCloneEmpty (N_Vector w)
{
  N_Vector v;

  /* Check input */
  if (w == nullptr)
    return nullptr;

  /* Create vector */
  v = nullptr;
  #  if defined (HAVE_SUNDIALS_SUNCONTEXT)
    v = octave::N_VNewEmpty(w->sunctx);
  #else
    v = octave::N_VNewEmpty();
  #endif
  if (v == nullptr)
    return nullptr;

  /* Attach operations */
  if (N_VCopyOps (w, v))
    {
      ::N_VDestroy (v);
      return nullptr;
    }

  return v;
}

N_Vector N_VClone (N_Vector w)
{
  void* content;

  /* check inputs */
  if (w == nullptr)
    return nullptr;

  /* Create an empty clone vector */
  N_Vector v = octave::N_VCloneEmpty (w);
  if (v == nullptr)
   return nullptr;
  /* Create content */
  content = nullptr;
  content = new ColumnVector (NV_LENGTH_C (w));

  if (content == nullptr)
   {
      ::N_VDestroy (v);
      return nullptr;
    }

  /* Attach content */
  v->content = content;

  return v;
}

void N_VDestroy (N_Vector v)
{
  if (v == nullptr)
   return;

  /* free content */
  if (v->content != nullptr)
   {
      // our ColumnVector was created with 'new', so use 'delete' to destroy it.
      delete NV_CONTENT_C (v);
      v->content = nullptr;
    }

  /* free ops and vector */
  if (v->ops != nullptr)
   {
    free (v->ops);
    v->ops = nullptr;
    }
  free (v);
   v = nullptr;

  return;
}

void N_VSpace (N_Vector v, sunindextype *lrw, sunindextype *liw)
{
  *lrw = NV_LENGTH_C (v);
  *liw = 1;

  return;
}

realtype *N_VGetArrayPointer (N_Vector v)
{
  return reinterpret_cast<realtype *> NV_DATA_C (v);
}

/*
 * Linear combination of two vectors: z = a*x + b*y
 */
void N_VLinearSum (realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{

  ColumnVector *xv, *yv, *zv;
  xv = static_cast <ColumnVector *> NV_CONTENT_C (x);
  yv = static_cast <ColumnVector *> NV_CONTENT_C (y);
  zv = static_cast <ColumnVector *> NV_CONTENT_C (z);
  booleantype test;

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE))
    {
      (*zv) = (*xv) + (*yv);
      return;
    }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE)))
      || ((a == -ONE) && (b == ONE)))
    {
      if(test == 0)
        (*zv) = (*yv) - (*xv);
      else
        (*zv) = (*xv) - (*yv);
      return;
    }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE))
   {
      if(test == 1)
        (*zv) = (*xv) + b * (*yv);
      else
        (*zv) = a * (*xv) + (*yv);
      return;
    }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE))
   {
      if(test == 1)
        (*zv) = b * (*yv) - (*xv);
      else
        (*zv) = a * (*xv) - (*yv);
      return;
    }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b)
    {
      (*zv) = a * ((*yv) + (*xv));
      return;
    }

  /* Case: a == -b */

  if (a == -b)
   {
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
void N_VConst (realtype c, N_Vector z)
{
  ColumnVector *zv;
  zv = static_cast<ColumnVector *> NV_CONTENT_C (z);

  zv->fill (c);

  return;
}

/*
 * Elementwise multiply vectors: z[i] = x[i]*y[i]
 */
void N_VProd (N_Vector x, N_Vector y, N_Vector z)
{

  ColumnVector *xv, *zv, *yv;
  xv = static_cast<ColumnVector *> NV_CONTENT_C (x);
  zv = static_cast<ColumnVector *> NV_CONTENT_C (z);
  yv = static_cast<ColumnVector *> NV_CONTENT_C (y);
  *zv = product((*xv),(*yv));

  return;
}

/*
 * Elementwise divide vectors: z[i] = x[i]/y[i]
 */
void N_VDiv (N_Vector x, N_Vector y, N_Vector z) 
{

  ColumnVector *xv, *zv, *yv;
  xv = static_cast<ColumnVector *> NV_CONTENT_C (x);
  zv = static_cast<ColumnVector *> NV_CONTENT_C (z);
  yv = static_cast<ColumnVector *> NV_CONTENT_C (y);
  *zv = quotient((*xv),(*yv));

  return;
}

/*
 * Scale vector: z = c*x
 */
void N_VScale (realtype c, N_Vector x, N_Vector z)
{
  ColumnVector *xv = const_cast <ColumnVector *> NV_CONTENT_C (x);
  ColumnVector *zv = const_cast <ColumnVector *> NV_CONTENT_C (z);

  // checks if both NVectors z and x
  // are stored in the same place
  if (zv == xv)
   {  /* BLAS usage: scale x <- cx */
      (*xv) *= c;
      // Not using this results in poor iterative
      // algorthm performance
      return;
    }
  if (c == ONE)
    VCopy_Local (x, z);  // substituting with (*zv) = (*xv) leads to wrong results.
  else if (c == -ONE)
    (*zv) -= (*xv);
  else
    (*zv) = c * (*xv);

  return;

}

/*
 * Elementwise absolute value: z[i] = |x[i]|
 */
void N_VAbs (N_Vector x, N_Vector z)
{
  ColumnVector *xv,*zv;
  xv = static_cast <ColumnVector *> NV_CONTENT_C (x);
  zv = static_cast <ColumnVector *> NV_CONTENT_C (z);

  (*zv) = xv->abs();
  return;
}

/*
 * Elementwise inverse: z[i] = 1/x[i]
 */
void N_VInv (N_Vector x, N_Vector z)
{
  ColumnVector *xv,*zv;
  xv = static_cast <ColumnVector *> NV_CONTENT_C (x);
  zv = static_cast <ColumnVector *> NV_CONTENT_C (z);

  (*zv) = ONE / (*xv);

  return;
}

/*
 * Add constant: z = x + b
 */
void N_VAddConst (N_Vector x, realtype b, N_Vector z)
{
  ColumnVector *xv,*zv;
  xv = static_cast <ColumnVector *> NV_CONTENT_C (x);
  zv = static_cast <ColumnVector *> NV_CONTENT_C (z);

  (*zv) = (*xv) + b;
  return;
}

/*
 * Scalar product of vectors x and y
 */
realtype N_VDotProd (N_Vector x, N_Vector y)
{
  ColumnVector *xv,*yv;
  realtype sum;
  xv = static_cast <ColumnVector *> NV_CONTENT_C (x);
  yv = static_cast <ColumnVector *> NV_CONTENT_C (y);
  sum = (*xv).transpose() * (*yv);
  return sum;
}

/*
 * Max norm (L infinity) of vector x
 */
realtype N_VMaxNorm (N_Vector x)
{
  ColumnVector *xv, abret;
  realtype ret;
  xv = static_cast <ColumnVector *> NV_CONTENT_C (x);
  abret = static_cast <ColumnVector> (xv->abs());
  ColumnVector *abp = &(abret);
  ret = abp->max();

  return ret;
}

/*
 * Weighted RMS norm
 */
realtype N_VWrmsNorm (N_Vector x, N_Vector w)
{
  octave_value_list ov = ovl ((octave::N_VWSqrSumLocal (x, w)/(NV_LENGTH_C (x))));
  return (octave::Fsqrt (ov,1)(0)).double_value();
}

/*
 * local weighted squared sum
 */
realtype N_VWSqrSumLocal (N_Vector x, N_Vector w)
{
  N_Vector prod = octave::N_VClone (x);
  ColumnVector *xv,*yv, *pv;
  realtype sum;
  xv = static_cast <ColumnVector *> NV_CONTENT_C (x);
  yv = static_cast <ColumnVector *> NV_CONTENT_C (w);
  pv = static_cast <ColumnVector *> NV_CONTENT_C (prod);

  (*pv) = product ((*xv),(*yv));
  octave_value_list ov = ovl ((*pv));
  sum = (octave::Fsumsq (ov,1)(0)).double_value();
  return sum;
}

/*
 * Masked weighted RMS norm
 */
realtype N_VWrmsNormMask (N_Vector x, N_Vector w, N_Vector id)
{
  octave_value_list ov = ovl((octave::N_VWSqrSumMaskLocal (x, w, id)/(NV_LENGTH_C (x))));
  return (octave::Fsqrt (ov,1)(0)).double_value();
}

/*
 * local weighted masked squared sum
 */
realtype N_VWSqrSumMaskLocal (N_Vector x, N_Vector w, N_Vector id)
{
  N_Vector prod = octave::N_VClone (x);
  N_Vector prod2 = octave::N_VClone (x);
  ColumnVector *xv,*yv, *mv, *pv, *pv2;
  realtype sum;
  xv = static_cast <ColumnVector *> NV_CONTENT_C (x);
  yv = static_cast <ColumnVector *> NV_CONTENT_C (w);
  mv = static_cast <ColumnVector *> NV_CONTENT_C (id);    //mask vector
  pv = static_cast <ColumnVector *> NV_CONTENT_C (prod);
  pv2 = static_cast <ColumnVector *> NV_CONTENT_C (prod2);

  const octave_value_list ov = octave_value_list ({(*mv),(*xv),(*mv)});
  const octave_value_list ov2 = octave_value_list ({(*mv),(*yv),(*mv)});
  (*pv) = (octave::Fmerge(ov2,1)(0)).column_vector_value();
  (*pv2) = (octave::Fmerge(ov2,1)(0)).column_vector_value();
  sum = octave::N_VWSqrSumLocal(prod, prod2);

  return sum;
}

/* Minimum value in Vector */
realtype N_VMin (N_Vector x)
{
  ColumnVector *xv;
  realtype min;
  xv = static_cast <ColumnVector *> NV_CONTENT_C (x);

  min = xv->min();
  return min;
}

/*
 * Elementwise z[i] = |x[i]| >= c ? 1 : 0
 */
void N_VCompare (realtype c, N_Vector x, N_Vector z)
{
  ColumnVector *xv,*zv;
  xv = static_cast <ColumnVector *> NV_CONTENT_C (x);
  zv = static_cast <ColumnVector *> NV_CONTENT_C (z);
  /*
  * FIXME : The comparisons here shouldn't be low-level,
  * but octave::Fge which should be used gives an interpreter error
  */
  for (int i = 0; i < (xv)->numel(); i++) 
      (*zv)(i) = ((xv->abs())(i) >= c) ? ONE : ZERO;
    

  return;
}

/*
 * Elementwise inverse with zero checking: z[i] = 1/x[i], x[i] != 0
 */
booleantype N_VInvTest (N_Vector x, N_Vector z)
{
  N_Vector y = octave::N_VClone (x);

  ColumnVector *xv, *zv, *yv;
  xv = static_cast <ColumnVector *> NV_CONTENT_C (x);
  zv = static_cast <ColumnVector *> NV_CONTENT_C (z);
  yv = static_cast <ColumnVector *> NV_CONTENT_C (y);
  booleantype no_zero_found = SUNTRUE ;
  (*zv) = ONE / (*xv);
  if(xv->nnz() != xv->numel())
    no_zero_found = SUNFALSE;
  if(no_zero_found==SUNTRUE)
    (*zv) = ONE / (*xv);
  else
   {
      (*yv) = (*zv);
      const octave_value_list ov = ovl ((*zv));
      (*yv) = (octave::Fisfinite (ov,1)(0)).column_vector_value();
      const octave_value_list ov2 = octave_value_list ({(*yv),(*zv),(*yv)});
      (*zv) = (octave::Fmerge (ov2,1)(0)).column_vector_value();
    }

  return no_zero_found;
}

/* FIXME : Same problem as Compare */
booleantype N_VConstrMask (N_Vector c, N_Vector x, N_Vector m)
{
  realtype temp;
  booleantype test;

  ColumnVector *cv, *xv, *mv; 
  xv = static_cast <ColumnVector *> NV_CONTENT_C (x);
  cv = static_cast <ColumnVector *> NV_CONTENT_C (c);
  mv = static_cast <ColumnVector *> NV_CONTENT_C (m);

  for (sunindextype i = 0; i < NV_LENGTH_C (x); i++) {
    (*mv)(i) = ZERO;

    /* Continue if no constraints were set for the variable */
    if ((*cv)(i) == ZERO)
      continue;

    /* Check if a set constraint has been violated */
    test = ((cv->abs())(i) > ONEPT5 && (*xv)(i)*(*cv)(i) <= ZERO)
             || ((cv->abs())(i)  > HALF && (*xv)(i)*(*cv)(i) <  ZERO);
    if (test)
     {
        temp = (*mv)(i) = ONE;
      }
  }

  /* Return false if any constraint was violated */
  return (temp == ONE) ? SUNFALSE : SUNTRUE;
}

/*
 * Find minimum quotient: minq  = min ( num[i]/denom[i]), denom[i] != 0.
 */
realtype N_VMinQuotient (N_Vector num, N_Vector denom)
{
    N_Vector z = octave::N_VClone (num);
    ColumnVector *xv, *zv, *yv;
    realtype min;
    min = BIG_REAL;

    xv = static_cast<ColumnVector *> NV_CONTENT_C (num);
    yv = static_cast<ColumnVector *> NV_CONTENT_C (denom);
    zv = static_cast<ColumnVector *> NV_CONTENT_C (z);
    *zv = quotient ((*xv),(*yv));
    min = (zv)->min();

  return min;
}

/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */

int N_VLinearCombination (int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  /* invalid number of vectors */
  if (nvec < 1)
    return -1;

  /* should have called N_VScale */
  if (nvec == 1)
   {
      octave::N_VScale (c[0], X[0], z);
      return 0;
    }

  /* should have called N_VLinearSum */
  if (nvec == 2)
   {
      octave::N_VLinearSum (c[0], X[0], c[1], X[1], z);
      return 0;
    }

  ColumnVector *xv, *zv;

  xv = static_cast<ColumnVector *> NV_CONTENT_C (X[0]);
  zv = static_cast<ColumnVector *> NV_CONTENT_C (z);

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */

  (*zv) = (*c)*(*xv);
  return 0;
}


int N_VScaleAddMulti (int nvec, realtype* a, N_Vector x, N_Vector* Y, N_Vector* Z)
{
  /* invalid number of vectors */
  if (nvec < 1) return -1;

  /* should have called N_VLinearSum */
  if (nvec == 1)
   {
      octave::N_VLinearSum (a[0], x, ONE, Y[0], Z[0]);
      return 0;
    }

  ColumnVector *xv, *zv, *yv;

  yv = static_cast<ColumnVector *> NV_CONTENT_C (Y[0]);
  zv = static_cast<ColumnVector *> NV_CONTENT_C (Z[0]);
  xv = static_cast<ColumnVector *> NV_CONTENT_C (x);

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  (*zv) = (*yv) + (*a)*(*xv);
  return 0;
}

/*
 * -----------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------
 */

int N_VLinearSumVectorArray (int nvec,
                                   realtype a, N_Vector* X,
                                   realtype b, N_Vector* Y,
                                   N_Vector* Z)
{

  ColumnVector *xv, *zv, *yv;
  int          i;

  /* invalid number of vectors */
  if (nvec < 1)
   return -1;

  /* should have called N_VLinearSum */
  if (nvec == 1)
   {
      octave::N_VLinearSum (a, X[0], b, Y[0], Z[0]);
      return 0;
    }

  for (i=0; i<nvec; i++)
   {
      yv = static_cast<ColumnVector *> NV_CONTENT_C (Y[i]);
      zv = static_cast<ColumnVector *> NV_CONTENT_C (Z[i]);
      xv = static_cast<ColumnVector *> NV_CONTENT_C (X[i]);
      (*zv) = a * (*xv) + b * (*yv);
    }

  return 0;
}


int N_VScaleVectorArray (int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  sunindextype i,j;
  realtype    *xd=nullptr;
  ColumnVector *xv, *zv;

  xv = static_cast<ColumnVector *> NV_CONTENT_C (X[0]);
  zv = static_cast<ColumnVector *> NV_CONTENT_C (Z[0]);

  /* invalid number of vectors */
  if (nvec < 1) return -1;

  /* should have called N_VScale */
  if (nvec == 1)
   {
      octave::N_VScale(c[0], X[0], Z[0]);
      return 0;
    }

  /*
   * X[i] *= c[i]
   */
  /* only checking if NVectors X and Z
   are stored in the same address */
  if (xv == zv)
   {
      // (*xv) *= (*c); //this should work but doesn't, in both N_VScale
      //                // instead results in convergence failures and wrong answer
      for (i=0; i<nvec; i++)
       {
          xd = NV_DATA_C (X[i]);
          for (j=0; j<(zv)->numel(); j++)
            xd[j] *= c[i];
        }
      return 0;
    }

  /*
   * Z[i] = c[i] * X[i]
   */
  (*zv) = (*c) * (*xv);
  return 0;
}

/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector operations
 * -----------------------------------------------------------------
 */


static void VCopy_Local (N_Vector x, N_Vector z)
// used in N_VScale and if just do (*zv)=(*xv) results in wrong answer.
{
  ColumnVector *xv, *zv;

  xv = static_cast<ColumnVector *> NV_CONTENT_C (x);
  zv = static_cast<ColumnVector *> NV_CONTENT_C (z);

  for (sunindextype i = 0; i < xv->numel(); i++)
    (*zv)(i) = (*xv)(i);

  return;
}

/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

int N_VEnableFusedOps (N_Vector v, booleantype tf)
{
  /* check that vector is non-nullptr */
  if (v == nullptr) return -1;

  /* check that ops structure is non-nullptr */
  if (v->ops == nullptr) return -1;

   if (tf)
    {
      /* enable all fused vector operations */
      v->ops->nvlinearcombination = nullptr;
      v->ops->nvscaleaddmulti     = octave::N_VScaleAddMulti;
      v->ops->nvdotprodmulti      = nullptr;
      /* enable all vector array operations */
      v->ops->nvlinearsumvectorarray         = octave::N_VLinearSumVectorArray;
      v->ops->nvscalevectorarray             = octave::N_VScaleVectorArray;
      v->ops->nvconstvectorarray             = nullptr;
      v->ops->nvwrmsnormvectorarray          = nullptr;
      v->ops->nvwrmsnormmaskvectorarray      = nullptr;
      v->ops->nvscaleaddmultivectorarray     = nullptr;
      v->ops->nvlinearcombinationvectorarray = nullptr;
    }
  else
   {
      /* disable all fused vector operations */
      v->ops->nvlinearcombination = nullptr;
      v->ops->nvscaleaddmulti     = nullptr;
      v->ops->nvdotprodmulti      = nullptr;
      /* disable all vector array operations */
      v->ops->nvlinearsumvectorarray         = nullptr;
      v->ops->nvscalevectorarray             = nullptr;
      v->ops->nvconstvectorarray             = nullptr;
      v->ops->nvwrmsnormvectorarray          = nullptr;
      v->ops->nvwrmsnormmaskvectorarray      = nullptr;
      v->ops->nvscaleaddmultivectorarray     = nullptr;
      v->ops->nvlinearcombinationvectorarray = nullptr;
    }

  /* return success */
  return 0;
}


/*
* -----------------------------------------------------------------
* exported functions for octave implementation of SUNMATRIX_SPARSE
* -----------------------------------------------------------------
*/

/* Constructor function to create a new sparse matrix
  based on dimensions */

SUNMatrix SUNSparseMatrix (sunindextype M, sunindextype N,
                          sunindextype NNZ ARG_SUNCONTEXT)
{
  SUNMatrix A;
  void *content;

  /* return with nullptr matrix on illegal input */
  if ((M <= 0) || (N <= 0) || (NNZ < 0))
    return nullptr;

  /* Create an empty matrix object */
  A = nullptr;
  A = SUNMatNewEmpty (SP_OCTAVE_SUNCONTEXT);
  if (A == nullptr)
    return nullptr;

  // /* Attach required operations */
  A->ops->getid = octave::SUNMatGetID_Sparse;
  A->ops->clone = nullptr;
  A->ops->destroy = octave::SUNMatDestroy_Sparse;
  A->ops->zero = octave::SUNMatZero_Sparse;
  A->ops->copy = nullptr;
  A->ops->scaleadd = nullptr;
  A->ops->scaleaddi = nullptr;
  A->ops->matvec = nullptr;
  A->ops->space = nullptr;

  /* Create content */
  content = nullptr;
  // creating a new SparseMatrix Octave data structure with rows and columns as specified
  content = new SparseMatrix (M, N, NNZ);
  if (content == nullptr)
    {
      SUNMatDestroy (A);
      return nullptr;
    }

  /* Attach content */
  A->content = content;
  /* Fill content */
  return A;
}

/* Functions to access the contents
  of the sparse matrix structure */

sunindextype SUNSparseMatrix_Rows (SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
    return SM_ROWS_O(A);
  else
    return SUNMAT_ILL_INPUT;
}

sunindextype SUNSparseMatrix_Columns (SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
    return SM_COLUMNS_O(A);
  else
    return SUNMAT_ILL_INPUT;
}

sunindextype SUNSparseMatrix_NNZ (SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
    return SM_NNZ_O(A);
  else
    return SUNMAT_ILL_INPUT;
}

sunindextype SUNSparseMatrix_NP (SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
    return SM_NP_O(A);
  else
    return SUNMAT_ILL_INPUT;
}

int SUNSparseMatrix_SparseType (SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
    return SM_SPARSETYPE_O(A);
  else
    return SUNMAT_ILL_INPUT;
}

realtype *SUNSparseMatrix_Data (SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
    return SM_DATA_O(A);
  else
    return nullptr;
}

sunindextype *SUNSparseMatrix_IndexValues (SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
    return SM_INDEXVALS_O(A);
  else
    return nullptr;
}

sunindextype *SUNSparseMatrix_IndexPointers (SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
    return SM_INDEXPTRS_O(A);
  else
    return nullptr;
}

/*
 * -----------------------------------------------------------------
 * implementation of matrix operations
 * -----------------------------------------------------------------
 */

SUNMatrix_ID SUNMatGetID_Sparse ([[maybe_unused]] SUNMatrix A)
{
  return SUNMATRIX_CUSTOM;
}

void SUNMatDestroy_Sparse (SUNMatrix A)
{
  if (A == nullptr)
    return;

  //   /* free content */
  if (A->content != nullptr)
    {
      /* deleting object created using new */
      delete SM_CONTENT_O(A);
    }
  A->content = nullptr;

  /* free ops and matrix */
  if (A->ops)
    {
      free(A->ops);
      A->ops = nullptr;
    }
  free(A);
  A = nullptr;

  return;
}

int SUNMatZero_Sparse (SUNMatrix A)
{
  SparseMatrix *ptr = new SparseMatrix (SM_ROWS_O(A),SM_COLS_O(A));
  SparseMatrix *am = static_cast <SparseMatrix *> SM_CONTENT_O(A);
  *am = *ptr;

  return SUNMAT_SUCCESS;
}

/* Function to print the sparse matrix */

void SUNSparseMatrix_Print (SUNMatrix A)
{
  SparseMatrix *am = static_cast <SparseMatrix *>SM_CONTENT_O(A);
  std::cout<<(*am);
  return;
}


/* ----------------------------------------------------------------------------
 * Function to create a new General Octave linear solver
 */

SUNLinearSolver SUNLinSol_Gen (N_Vector y, SUNMatrix A ARG_SUNCONTEXT)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_GEN content;

  if (octave::SUNSparseMatrix_Rows(A) != octave::SUNSparseMatrix_Columns(A))
   return nullptr;

  if ( ::N_VGetVectorID(y) != SUNDIALS_NVEC_CUSTOM )
    return nullptr;

  if (octave::SUNSparseMatrix_Rows(A) != ::N_VGetLength(y))
   return nullptr;


  /* Create an empty linear solver */
  S = nullptr;
  S = SUNLinSolNewEmpty (SP_OCTAVE_SUNCONTEXT);
  if (S == nullptr)
    return nullptr;

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_Gen;
  S->ops->getid      = SUNLinSolGetID_Gen;
  S->ops->initialize = nullptr;
  S->ops->setup      = nullptr;
  S->ops->solve      = SUNLinSolSolve_Gen;
  S->ops->lastflag   = nullptr;
  S->ops->space      = nullptr;
  S->ops->free       = SUNLinSolFree_Gen;

    /* Create content */
  content = nullptr;
  content = reinterpret_cast<SUNLinearSolverContent_GEN> (malloc(sizeof *content));
  if (content == nullptr) { SUNLinSolFree(S); return nullptr; }

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->last_flag       = 0;
  content->first_factorize = 1;

  return S;
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_Gen ([[maybe_unused]] SUNLinearSolver S)
{
  return SUNLINEARSOLVER_DIRECT;
}


SUNLinearSolver_ID SUNLinSolGetID_Gen ([[maybe_unused]] SUNLinearSolver S)
{
  return SUNLINEARSOLVER_CUSTOM;
}


int SUNLinSolSolve_Gen (SUNLinearSolver S, SUNMatrix A, N_Vector x,
                       N_Vector b, [[maybe_unused]] realtype tol)
{
  /* check for valid inputs */
  if ( (A == nullptr) || (S == nullptr) || (x == nullptr) || (b == nullptr) )
    return SUNLS_MEM_NULL;
  ColumnVector *xv = static_cast <ColumnVector *> NV_CONTENT_C (x);
  ColumnVector *zv = static_cast <ColumnVector *> NV_CONTENT_C (b);

  SparseMatrix *am = static_cast<SparseMatrix *> SM_CONTENT_O (A);

  (*xv) = am->solve((*zv));

  return 1;
}

int SUNLinSolFree_Gen (SUNLinearSolver S)
{
  /* return with success if already freed */
  if (S == nullptr)
   return SUNLS_SUCCESS;

  /* delete generic structures */
  if (S->ops)
    {
      free(S->ops);
      S->ops = nullptr;
    }
  free(S); S = nullptr;
  return SUNLS_SUCCESS;
}

OCTAVE_END_NAMESPACE(octave)
