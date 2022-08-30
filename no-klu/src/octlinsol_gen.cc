/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * Based on codes <solver>_Gen.c, written by Carol Woodward @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the KLU implementation of
 * the OCTLinSol package.
 * -----------------------------------------------------------------*/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "octlinsol_gen.h"
#include <sundials/sundials_math.h>

#define ZERO      RCONST(0.0)
#define ONE       RCONST(1.0)
#define TWO       RCONST(2.0)
#define TWOTHIRDS RCONST(0.666666666666666666666666666666667)


/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new General Octave linear solver
 */

SUNLinearSolver OCTLinSol_Gen(N_Vector y, SUNMatrix A, SUNContext sunctx)
{
  SUNLinearSolver S;
  OCTLinearSolverContent_GEN content;

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  // if (SUNMatGetID(A) != SUNMATRIX_SPARSE) return(NULL);
  if (SUNMatGetID(A) == SUNMATRIX_SPARSE)
  {
    if (OCTSparseMatrix_Rows(A) != OCTSparseMatrix_Columns(A)) return(NULL);

    // if ( (N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL) &&
    //      (N_VGetVectorID(y) != SUNDIALS_NVEC_OPENMP) &&
    //      (N_VGetVectorID(y) != SUNDIALS_NVEC_PTHREADS) )
    //   return(NULL);

    if (OCTSparseMatrix_Rows(A) != N_VGetLength(y)) return(NULL);
  }
  else if (SUNMatGetID(A) == SUNMATRIX_DENSE)
  {
    sunindextype MatrixRows;
    if (OCTDenseMatrix_Rows(A) != OCTDenseMatrix_Columns(A)) return(NULL);

    if ( (N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL) &&
        (N_VGetVectorID(y) != SUNDIALS_NVEC_OPENMP) &&
        (N_VGetVectorID(y) != SUNDIALS_NVEC_PTHREADS) )
      return(NULL);

    MatrixRows = OCTDenseMatrix_Rows(A);
    if (MatrixRows != N_VGetLength(y)) return(NULL);
  }

    /* Create an empty linear solver */
    S = NULL;
    S = SUNLinSolNewEmpty(sunctx);
    if (S == NULL) return(NULL);

    /* Attach operations */
    S->ops->gettype    = OCTLinSolGetType_Gen;
    S->ops->getid      = OCTLinSolGetID_Gen;
    S->ops->initialize = NULL;
    S->ops->setup      = NULL;
    S->ops->solve      = OCTLinSolSolve_Gen;
    S->ops->lastflag   = NULL;
    S->ops->space      = NULL;
    S->ops->free       = OCTLinSolFree_Gen;

      /* Create content */
    content = NULL;
    content = (OCTLinearSolverContent_GEN) malloc(sizeof *content);
    if (content == NULL) { SUNLinSolFree(S); return(NULL); }

    /* Attach content */
    S->content = content;

    /* Fill content */
    content->last_flag       = 0;
    content->first_factorize = 1;

  // #if defined(SUNDIALS_INT64_T)
  //   if (OCTSparseMatrix_SparseType(A) == CSC_MAT) {
  //     content->klu_solver = (KLUSolveFn) &klu_l_solve;
  //   } else {
  //     content->klu_solver = (KLUSolveFn) &klu_l_tsolve;
  //   }
  // #elif defined(SUNDIALS_INT32_T)
  //   if (OCTSparseMatrix_SparseType(A) == CSC_MAT) {
  //     content->klu_solver = &klu_solve;
  //   } else {
  //     content->klu_solver = &klu_tsolve;
  //   }
  // #else  /* incompatible sunindextype for KLU */
  // #error  Incompatible sunindextype for KLU
  // #endif

    return(S);
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type OCTLinSolGetType_Gen(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DIRECT);
}


SUNLinearSolver_ID OCTLinSolGetID_Gen(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_CUSTOM);
}


int OCTLinSolSolve_Gen(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                       N_Vector b, realtype tol)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
  {
    /* check for valid inputs */
    if ( (A == NULL) || (S == NULL) || (x == NULL) || (b == NULL) )
      return(SUNLS_MEM_NULL);
    
    std::cout<<"id "<<N_VGetVectorID(b)<<" len ";
    sunindextype nx = N_VGetLength_Octave(x);
    sunindextype nb = N_VGetLength_Octave(b);
    /* copy b into x */
    // N_VScale_Octave(ONE, b, x);
    std::cout<<"151 id "<<N_VGetVectorID(b)<<" len "<<nx<<"\n";
    ColumnVector *xv = new ColumnVector(nx);
    ColumnVector *zv = new ColumnVector(nb);

    xv = static_cast <ColumnVector *> NV_CONTENT_C(x);
    zv = static_cast <ColumnVector *> NV_CONTENT_C(b);
    
    // realtype *puntx = NV_DATA_C(x);
    // realtype *puntb = NV_DATA_C(b);

    //   for (octave_f77_int_type i = 0; i < nx; i++){
    //     std::cout<<puntx[i]<<" "<<puntb[i]<<"\n";
    //     (*xv)(i) = puntx[i];
    //     std::cout<<i<<" here";
    //   }
    //   for (octave_f77_int_type i = 0; i < nb; i++){
    //     std::cout<<puntx[i]<<" "<<puntb[i]<<"\n";
    //     (*zv)(i) = puntb[i];
    //     std::cout<<i<<" here";
    //   }

      // std::cout<<"and now \n"<<(*xv)<<"done \n";
    // ColumnVector xv = NVecToCol(x,NV_LENGTH_C(x));
    // ColumnVector zv = NVecToCol(x,NV_LENGTH_C(b));
    // std::cout<<"id "<<N_VGetVectorID(x)<<" len "<<nx<<"\n";
    // N_VPrint(x);
    // N_VPrint(b);
    SparseMatrix *am;
    am = static_cast<SparseMatrix *> SM_CONTENT_O(A);
    std::cout<<"183 id "<<N_VGetVectorID(b)<<"\n";

    printf("\n\n solving once over \n\n");

    (*xv) = am->solve((*zv));
    // N_VPrint(x);
    // N_VPrint(b);
    std::cout<<(*xv);

    return 1;
  }
  // else if ( SUNMatGetID(A) == SUNMATRIX_DENSE )
  // {
  //   /* check for valid inputs */
  //   if ( (A == NULL) || (S == NULL) || (x == NULL) || (b == NULL) )
  //     return(SUNLS_MEM_NULL);

  //   /* copy b into x */
  //   N_VScale_Octave(ONE, b, x);

  //   ColumnVector *xv, *zv;
  //   xv = static_cast<ColumnVector *> NV_CONTENT_C(x);
  //   zv = static_cast<ColumnVector *> NV_CONTENT_C(b);

  //   Matrix *am;
  //   am = static_cast<Matrix *> SM_CONTENT_D(A);

  //   (*xv) = am->solve((*zv));
  // }
  
}

  ColumnVector NVecToCol (N_Vector& v, octave_f77_int_type n)
    {
      ColumnVector data (n);
      realtype *punt = NV_DATA_C(v);

      for (octave_f77_int_type i = 0; i < n; i++)
        data(i) = punt[i];

      // &data = static_cast <ColumnVector *> NV_CONTENT_C(v);
      return data;
    }

int OCTLinSolFree_Gen(SUNLinearSolver S)
{
  /* return with success if already freed */
  if (S == NULL) return(SUNLS_SUCCESS);

  /* delete generic structures */
  if (S->ops) {
    free(S->ops);
    S->ops = NULL;
  }
  free(S); S = NULL;
  return(SUNLS_SUCCESS);
}



