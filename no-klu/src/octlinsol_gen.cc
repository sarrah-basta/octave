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

  if (OCTSparseMatrix_Rows(A) != OCTSparseMatrix_Columns(A)) return(NULL);

  if ( N_VGetVectorID(y) != SUNDIALS_NVEC_CUSTOM ) 
    return(NULL);

  if (OCTSparseMatrix_Rows(A) != N_VGetLength(y)) return(NULL);
 

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
  /* check for valid inputs */
  if ( (A == NULL) || (S == NULL) || (x == NULL) || (b == NULL) )
    return(SUNLS_MEM_NULL);
  ColumnVector *xv = static_cast <ColumnVector *> NV_CONTENT_C(x);
  ColumnVector *zv = static_cast <ColumnVector *> NV_CONTENT_C(b);

  SparseMatrix *am = static_cast<SparseMatrix *> SM_CONTENT_O(A);

  (*xv) = am->solve((*zv));

  return 1;
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



