/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * Based on sundials_GenOCTLinSol_Gen_impl.h and arkode_GenOCTLinSol_Gen.h/cvode_GenOCTLinSol_Gen.h/...
 *     code, written by Carol S. Woodward @ LLNL
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
 * This is the header file for the KLU implementation of the
 * OCTLinSol module, OCTLinSol_Gen.
 *
 * Note:
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 * -----------------------------------------------------------------
 */


#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
// #include <sunmatrix/sunmatrix_sparse.h>
#include "octmatrix_sparse.h"
#include "octmatrix_dense.h"
#include "nvector_octave.h"
#include <dSparse.h>
#include <dMatrix.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

struct _OCTLinearSolverContent_GEN {
  int              last_flag;
  int              first_factorize;
};

typedef struct _OCTLinearSolverContent_GEN *OCTLinearSolverContent_GEN;
/* -------------------------------------
 * Exported Functions for OCTLINSOL_GEN
 ------------ */

SUNDIALS_EXPORT SUNLinearSolver OCTLinSol_Gen(N_Vector y, SUNMatrix A, SUNContext sunctx);
// SUNDIALS_EXPORT int OCTLinSol_GenReInit(SUNLinearSolver S, SUNMatrix A,
//                                         sunindextype nnz, int reinit_type);
// SUNDIALS_EXPORT int OCTLinSol_GenSetOrdering(SUNLinearSolver S,
//                                              int ordering_choice);


/* -----------------------------------------------
 *  Implementations of SUNLinearSolver operations
 * ----------------------------------------------- */

SUNDIALS_EXPORT SUNLinearSolver_Type OCTLinSolGetType_Gen(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID OCTLinSolGetID_Gen(SUNLinearSolver S);
// SUNDIALS_EXPORT int OCTLinSolInitialize_GenOCTLinSol_Gen(SUNLinearSolver S);
// SUNDIALS_EXPORT int OCTLinSolSetup_GenOCTLinSol_Gen(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int OCTLinSolSolve_Gen(SUNLinearSolver S, SUNMatrix A,
                                       N_Vector x, N_Vector b, realtype tol);
// SUNDIALS_EXPORT sunindextype OCTLinSolLastFlag_GenOCTLinSol_Gen(SUNLinearSolver S);
// SUNDIALS_EXPORT int OCTLinSolSpace_GenOCTLinSol_Gen(SUNLinearSolver S,
                                      //  long int *lenrwLS,
                                      //  long int *leniwLS);
SUNDIALS_EXPORT int OCTLinSolFree_Gen(SUNLinearSolver S);


#ifdef __cplusplus
}
#endif
