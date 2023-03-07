// Copyright (C) 2016-2023 The Octave Project Developers
//
// See the file COPYRIGHT.md in the top-level directory of this
// distribution or <https://octave.org/copyright/>.
//
// This file is part of Octave.
//
// Octave is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Octave is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Octave; see the file COPYING.  If not, see
// <https://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////

/* -----------------------------------------------------------------
 * generating custom SUNDIALS implementations
 * to be able to allow Octave's native classes and
 * solvers to be utilized by the SUNDIALS IDA Class
 * ----------------------------------------------------------------- */

#ifndef OCTAVE_OCT_SUNDIALS_H
#define OCTAVE_OCT_SUNDIALS_H 1
#include "octave-config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <dColVector.h>
#include <dSparse.h>
#include <dMatrix.h>
#include "oct-locbuf.h"
#include "quit.h"

#include <oct.h>
#include <builtin-defun-decls.h>
#include <data.h>
#include <ov.h>
#include <ovl.h>

/* FIXME : This flag currently checks for a
  lot of Sundials features, such as nvecserial,
  etc, which might not be needed for just these files.
  Maybe a simpler configure flag should be made for
  this. */

#if defined(HAVE_SUNDIALS)
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#endif

/* -----------------------------------------------------------------
 * Useful mathematical macros for computation
 * ----------------------------------------------------------------- */
#define ZERO RCONST(0.0)
#define HALF RCONST(0.5)
#define ONE RCONST(1.0)
#define ONEPT5 RCONST(1.5)
#define TWO RCONST(2.0)
#define TWOTHIRDS RCONST(0.666666666666666666666666666666667)

/* -----------------------------------------------------------------
 * SUNContext : Since SUNDIALS 6.0.0, all SUNDIALS objects
 * (vectors, linear and nonlinear solvers, matrices, etc.)
 * hold a reference to a common simulation context object
 * defined by the SUNContext class. The macro ensures these
 * declarations maintain backward compatibility with versions
 * before 6.0.0 as well.
 * ----------------------------------------------------------------- */
#if (HAVE_SUNDIALS_SUNCONTEXT)
#define SP_ARG_SUNCONTEXT SUNContext m_sunContext
#define ARG_SUNCONTEXT , SUNContext m_sunContext
#define OCTAVE_SUNCONTEXT , m_sunContext
#define SP_OCTAVE_SUNCONTEXT m_sunContext
#else
#define ARG_SUNCONTEXT
#define OCTAVE_SUNCONTEXT
#define SP_ARG_SUNCONTEXT
#define SP_OCTAVE_SUNCONTEXT
#endif

/*
 * -----------------------------------------------------------------
 * Macros for access to octave implementation of N_Vector
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_C(v) (reinterpret_cast<ColumnVector *>(v->content))
#define NV_LENGTH_C(v) (NV_CONTENT_C(v)->numel())
#define NV_DATA_C(v) (NV_CONTENT_C(v)->fortran_vec())
#define NV_Ith_C(v, i) (NV_DATA_C(v)[i])

OCTAVE_BEGIN_NAMESPACE(octave)

/* -----------------------------------------------------------------
 * NVECTOR API functions exported for Octave's custom implementation
 * ----------------------------------------------------------------- */

N_Vector N_VNew (int vec_length ARG_SUNCONTEXT);
N_Vector N_VNewEmpty (SP_ARG_SUNCONTEXT);
N_Vector N_VMake (const ColumnVector &v ARG_SUNCONTEXT);

sunindextype N_VGetLength (N_Vector v);
void N_VPrint (N_Vector v);
N_Vector_ID N_VGetVectorID (N_Vector v);
N_Vector N_VCloneEmpty (N_Vector w);
N_Vector N_VClone (N_Vector w);
void N_VDestroy (N_Vector v);
void N_VSpace (N_Vector v, sunindextype *lrw, sunindextype *liw);
double *N_VGetArrayPointer (N_Vector v);

/* standard vector operations */
void N_VLinearSum (double a, N_Vector x, double b, N_Vector y, N_Vector z);
void N_VConst (double c, N_Vector z);
void N_VProd (N_Vector x, N_Vector y, N_Vector z);
void N_VDiv (N_Vector x, N_Vector y, N_Vector z);
void N_VScale (double c, N_Vector x, N_Vector z);
void N_VAbs (N_Vector x, N_Vector z);
void N_VInv (N_Vector x, N_Vector z);
void N_VAddConst (N_Vector x, double b, N_Vector z);
double N_VDotProd (N_Vector x, N_Vector y);
double N_VMaxNorm (N_Vector x);
double N_VWrmsNorm (N_Vector x, N_Vector w);
double N_VWrmsNormMask (N_Vector x, N_Vector w, N_Vector id);
double N_VMin (N_Vector x);
void N_VCompare (double c, N_Vector x, N_Vector z);
booleantype N_VInvTest (N_Vector x, N_Vector z);
booleantype N_VConstrMask (N_Vector c, N_Vector x, N_Vector m);
double N_VMinQuotient (N_Vector num, N_Vector denom);

/* fused vector operations */
int N_VScaleAddMulti (int nvec, double *a, N_Vector x,
                     N_Vector *Y, N_Vector *Z);

/* vector array operations */
int N_VLinearSumVectorArray (int nvec,
                            double a, N_Vector *X,
                            double b, N_Vector *Y,
                            N_Vector *Z);
int N_VScaleVectorArray (int nvec, double *c,
                        N_Vector *X, N_Vector *Z);

/* OPTIONAL local reduction kernels (no parallel communication) */
// double
double N_VWSqrSumLocal (N_Vector x, N_Vector w);
// double
double N_VWSqrSumMaskLocal (N_Vector x, N_Vector w,
                           N_Vector id);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */
int N_VEnableFusedOps (N_Vector v, booleantype tf);

/* ------------------------------------------------------
 * Macros for access to Octave implementation of SUNMATRIX
 * ----------------------------------------------------- */
#define CSC_MAT 0 /* Matrix Type Definition */
#define SM_CONTENT_O(A) (reinterpret_cast<SparseMatrix *>(A->content))
#define SM_ROWS_O(A) (SM_CONTENT_O(A)->rows())
#define SM_COLS_O(A) (SM_CONTENT_O(A)->cols())
#define SM_COLUMNS_O(A) (SM_CONTENT_O(A)->columns())
#define SM_NNZ_O(A) (SM_CONTENT_O(A)->nnz())
/* directly setting =N as Octave only supports CSC Representation */
#define SM_NP_O(A) (SM_CONTENT_O(A)->columns())
#define SM_SPARSETYPE_O(A) (CSC_MAT)
#define SM_DATA_O(A) (SM_CONTENT_O(A)->data())
#define SM_INDEXVALS_O(A) (SM_CONTENT_O(A)->ridx())
#define SM_INDEXPTRS_O(A) (SM_CONTENT_O(A)->cidx())

/* ----------------------------------------
 * SUNMATRIX API functions exported for Octave's
 * custom implementation
 * ---------------------------------------- */

SUNMatrix SUNSparseMatrix (sunindextype M, sunindextype N,
                          sunindextype NNZ ARG_SUNCONTEXT);
void SUNSparseMatrix_Print (SUNMatrix A);

sunindextype SUNSparseMatrix_Rows (SUNMatrix A);
sunindextype SUNSparseMatrix_Columns (SUNMatrix A);
sunindextype SUNSparseMatrix_NNZ (SUNMatrix A);
sunindextype SUNSparseMatrix_NP (SUNMatrix A);
int SUNSparseMatrix_SparseType (SUNMatrix A);
realtype *SUNSparseMatrix_Data (SUNMatrix A);
sunindextype *SUNSparseMatrix_IndexValues (SUNMatrix A);
sunindextype *SUNSparseMatrix_IndexPointers (SUNMatrix A);

SUNMatrix_ID SUNMatGetID_Sparse (SUNMatrix A);
void SUNMatDestroy_Sparse (SUNMatrix A);
int SUNMatZero_Sparse (SUNMatrix A);

/*
 * ----------------------------------------------------------------------------
 * Octave implementation of SUNLinearSolver
 * ----------------------------------------------------------------------------
 */
struct _SUNLinearSolverContent_GEN
{
  int last_flag;
  int first_factorize;
};

/* --------------------------------------------------
 * Declarations of Octave's custom implementation of
 * a Matrix Linear Solver
 * ------------------------------------------------- */
typedef struct _SUNLinearSolverContent_GEN *SUNLinearSolverContent_GEN;

/* -------------------------------------
 * SUNLINEARSOLVER API functions exported for Octave's
 * custom implementation
 * -----------------------------------*/

SUNLinearSolver SUNLinSol_Gen (N_Vector y, SUNMatrix A ARG_SUNCONTEXT);
SUNLinearSolver_Type SUNLinSolGetType_Gen (SUNLinearSolver S);
SUNLinearSolver_ID SUNLinSolGetID_Gen (SUNLinearSolver S);
int SUNLinSolSolve_Gen (SUNLinearSolver S, SUNMatrix A,
                       N_Vector x, N_Vector b, realtype tol);
int SUNLinSolFree_Gen (SUNLinearSolver S);

OCTAVE_END_NAMESPACE(octave)
