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
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_band.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* ---------------------------------------
 * Macros for access to SUNMATRIX_Dense
 * --------------------------------------- */

#define SM_CONTENT_S(A)     ( (Matrix *)(A->content) )

#define SM_ROWS_S(A)        ( SM_CONTENT_S(A)->rows() )

#define SM_COLS_S(A)        ( SM_CONTENT_S(A)->cols() )

#define SM_LDATA_S(A)       ( int(SM_ROWS_S(A)*SM_COLS_S(A)) )

#define SM_DATA_S(A)        ( SM_CONTENT_S(A)->fortran_vec() )

#define SM_ELEMENT_S(A,i,j)  ( SM_CONTENT_S(A)->elem(i,j) )

/* ----------------------------------------
 * Exported Functions for SUNMATRIX_Dense
 * ---------------------------------------- */

SUNDIALS_EXPORT SUNMatrix OCTDenseMatrix(sunindextype M, sunindextype N,
                                        SUNContext sunctx);


// SUNDIALS_EXPORT int OCTDenseMatrix_ToCSR(const SUNMatrix A, SUNMatrix* Bout);
// SUNDIALS_EXPORT int OCTDenseMatrix_ToCSC(const SUNMatrix A, SUNMatrix* Bout);

// SUNDIALS_EXPORT int OCTDenseMatrix_Realloc(SUNMatrix A);

SUNDIALS_EXPORT int OCTDenseMatrix_Reallocate(SUNMatrix A, sunindextype NNZ);

SUNDIALS_EXPORT void OCTDenseMatrix_Print(SUNMatrix A, FILE* outfile);

SUNDIALS_EXPORT sunindextype OCTDenseMatrix_Rows(SUNMatrix A);
SUNDIALS_EXPORT sunindextype OCTDenseMatrix_Columns(SUNMatrix A);
SUNDIALS_EXPORT sunindextype OCTDenseMatrix_NNZ(SUNMatrix A);
SUNDIALS_EXPORT sunindextype OCTDenseMatrix_NP(SUNMatrix A);
SUNDIALS_EXPORT realtype* OCTDenseMatrix_Data(SUNMatrix A);
SUNDIALS_EXPORT sunindextype OCTDenseMatrix_LData(SUNMatrix A);

SUNDIALS_EXPORT SUNMatrix_ID OCTMatGetID_Dense(SUNMatrix A);
SUNDIALS_EXPORT SUNMatrix OCTMatClone_Dense(SUNMatrix A);
SUNDIALS_EXPORT void OCTMatDestroy_Dense (SUNMatrix A);
SUNDIALS_EXPORT int OCTMatZero_Dense (SUNMatrix A);
SUNDIALS_EXPORT int OCTMatCopy_Dense (SUNMatrix A, SUNMatrix B);
// SUNDIALS_EXPORT int SUNMatScaleAdd_Dense(realtype c, SUNMatrix A, SUNMatrix B);
// SUNDIALS_EXPORT int SUNMatScaleAddI_Dense(realtype c, SUNMatrix A);
SUNDIALS_EXPORT int OCTMatMatvec_Dense(SUNMatrix A, N_Vector x, N_Vector y);
SUNDIALS_EXPORT int OCTMatSpace_Dense(SUNMatrix A, long int *lenrw, long int *leniw);


#ifdef __cplusplus
}
#endif


