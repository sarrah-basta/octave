/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner @ LLNL
 * Based on code sundials_Dense.h by: Carol Woodward and
 *     Slaven Peles @ LLNL, and Daniel R. Reynolds @ SMU
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
 * This is the header file for the Dense implementation of the
 * SUNMATRIX module, SUNMATRIX_Dense.
 *
 * Notes:
 *   - The definition of the generic SUNMatrix structure can be found
 *     in the header file sundials_matrix.h.
 *   - The definition of the type 'realtype' can be found in the
 *     header file sundials_types.h, and it may be changed (at the
 *     configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type 'booleantype' and 'indextype'.
 * -----------------------------------------------------------------
 */


#include <stdio.h>
#include <dMatrix.h>
#include <sundials/sundials_matrix.h>
#include <dColVector.h>
#include <builtin-defun-decls.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* ---------------------------------------
 * Macros for access to SUNMATRIX_Dense
 * --------------------------------------- */

#define SM_CONTENT_D(A)     ( (Matrix *)(A->content) )

#define SM_ROWS_D(A)        ( SM_CONTENT_D(A)->rows() )

#define SM_COLUMNS_D(A)        ( SM_CONTENT_D(A)->cols() )

#define SM_LDATA_D(A)       ( int(SM_ROWS_D(A)*SM_COLUMNS_D(A)) )

#define SM_DATA_D(A)        ( SM_CONTENT_D(A)->fortran_vec() )

#define SM_COLUMN_D(A,j)    ( (SM_CONTENT_D(A)->column(j)).fortran_vec() )

#define SM_ELEMENT_D(A,i,j)  ( SM_CONTENT_D(A)->elem(i,j) )

/* ----------------------------------------
 * Exported Functions for SUNMATRIX_Dense
 * ---------------------------------------- */

SUNDIALS_EXPORT SUNMatrix OCTDenseMatrix(sunindextype M, sunindextype N,
                                        SUNContext sunctx);

SUNDIALS_EXPORT void OCTDenseMatrix_Print(SUNMatrix A, FILE* outfile);

SUNDIALS_EXPORT sunindextype OCTDenseMatrix_Rows(SUNMatrix A);
SUNDIALS_EXPORT sunindextype OCTDenseMatrix_Columns(SUNMatrix A);
SUNDIALS_EXPORT realtype *OCTDenseMatrix_Column(SUNMatrix A, sunindextype j);
SUNDIALS_EXPORT sunindextype OCTDenseMatrix_NNZ(SUNMatrix A);
SUNDIALS_EXPORT sunindextype OCTDenseMatrix_NP(SUNMatrix A);
SUNDIALS_EXPORT realtype* OCTDenseMatrix_Data(SUNMatrix A);
SUNDIALS_EXPORT sunindextype OCTDenseMatrix_LData(SUNMatrix A);

SUNDIALS_EXPORT SUNMatrix_ID OCTMatGetID_Dense(SUNMatrix A);
SUNDIALS_EXPORT SUNMatrix OCTMatClone_Dense(SUNMatrix A);
SUNDIALS_EXPORT void OCTMatDestroy_Dense (SUNMatrix A);
SUNDIALS_EXPORT int OCTMatZero_Dense (SUNMatrix A);
SUNDIALS_EXPORT int OCTMatCopy_Dense (SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int OCTMatScaleAddI_Dense(realtype c, SUNMatrix A);
SUNDIALS_EXPORT int OCTMatMatvec_Dense(SUNMatrix A, N_Vector x, N_Vector y);
SUNDIALS_EXPORT int OCTMatSpace_Dense(SUNMatrix A, long int *lenrw, long int *leniw);


#ifdef __cplusplus
}
#endif


