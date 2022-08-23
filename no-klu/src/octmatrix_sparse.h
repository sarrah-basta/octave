/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner @ LLNL
 * Based on code sundials_sparse.h by: Carol Woodward and
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
 * This is the header file for the sparse implementation of the
 * SUNMATRIX module, SUNMATRIX_SPARSE.
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

#include <dSparse.h>
#include <sundials/sundials_matrix.h>
#include "octmatrix_dense.h"
// #include <sunmatrix/sunmatrix_dense.h>
// #include <sunmatrix/sunmatrix_band.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ------------------------
 * Matrix Type Definitions
 * ------------------------ */

#define CSC_MAT 0
#define CSR_MAT 1



/* ---------------------------------------
 * Macros for access to SUNMATRIX_SPARSE
 * --------------------------------------- */

#define SM_CONTENT_S(A)     ( (SparseMatrix *)(A->content) )

#define SM_ROWS_S(A)        ( SM_CONTENT_S(A)->rows() )

#define SM_COLUMNS_S(A)     ( SM_CONTENT_S(A)->columns() )

#define SM_NNZ_S(A)         ( SM_CONTENT_S(A)->nnz() )

#define SM_NP_S(A)          ( SM_CONTENT_S(A)->columns() ) //directly setting =N as we are only doing CSC Representation

#define SM_SPARSETYPE_S(A)  ( CSC_MAT ) 

#define SM_DATA_S(A)        ( SM_CONTENT_S(A)->data() )

#define SM_INDEXVALS_S(A)   ( SM_CONTENT_S(A)->ridx() )

#define SM_INDEXPTRS_S(A)   ( SM_CONTENT_S(A)->cidx() )

/* ----------------------------------------
 * Exported Functions for SUNMATRIX_SPARSE
 * ---------------------------------------- */

SUNDIALS_EXPORT SUNMatrix OCTSparseMatrix(sunindextype M, sunindextype N,
                                          sunindextype NNZ, int sparsetype,
                                          SUNContext sunctx);

SUNDIALS_EXPORT SUNMatrix OCTSparseMatrix_Make(SparseMatrix sm, SUNContext sunctx);

SUNDIALS_EXPORT SUNMatrix OCTSparseFromDenseMatrix(SUNMatrix A,
                                                   realtype droptol,
                                                   int sparsetype);

SUNDIALS_EXPORT SUNMatrix OCTSparseFromBandMatrix(SUNMatrix A,
                                                  realtype droptol,
                                                  int sparsetype);

// SUNDIALS_EXPORT int OCTSparseMatrix_ToCSR(const SUNMatrix A, SUNMatrix* Bout);
// SUNDIALS_EXPORT int OCTSparseMatrix_ToCSC(const SUNMatrix A, SUNMatrix* Bout);

// SUNDIALS_EXPORT int OCTSparseMatrix_Realloc(SUNMatrix A);

SUNDIALS_EXPORT int OCTSparseMatrix_Reallocate(SUNMatrix A, sunindextype NNZ);

SUNDIALS_EXPORT void OCTSparseMatrix_Print(SUNMatrix A, FILE* outfile);

SUNDIALS_EXPORT sunindextype OCTSparseMatrix_Rows(SUNMatrix A);
SUNDIALS_EXPORT sunindextype OCTSparseMatrix_Columns(SUNMatrix A);
SUNDIALS_EXPORT sunindextype OCTSparseMatrix_NNZ(SUNMatrix A);
SUNDIALS_EXPORT sunindextype OCTSparseMatrix_NP(SUNMatrix A);
SUNDIALS_EXPORT int OCTSparseMatrix_SparseType(SUNMatrix A);
SUNDIALS_EXPORT realtype* OCTSparseMatrix_Data(SUNMatrix A);
SUNDIALS_EXPORT sunindextype* OCTSparseMatrix_IndexValues(SUNMatrix A);
SUNDIALS_EXPORT sunindextype* OCTSparseMatrix_IndexPointers(SUNMatrix A);

SUNDIALS_EXPORT SUNMatrix_ID OCTMatGetID_Sparse(SUNMatrix A);
SUNDIALS_EXPORT SUNMatrix OCTMatClone_Sparse(SUNMatrix A);
SUNDIALS_EXPORT void OCTMatDestroy_Sparse (SUNMatrix A);
SUNDIALS_EXPORT int OCTMatZero_Sparse (SUNMatrix A);
SUNDIALS_EXPORT int OCTMatCopy_Sparse (SUNMatrix A, SUNMatrix B);
// SUNDIALS_EXPORT int SUNMatScaleAdd_Sparse(realtype c, SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int OCTMatScaleAddI_Sparse(realtype c, SUNMatrix A);
SUNDIALS_EXPORT int OCTMatMatvec_Sparse(SUNMatrix A, N_Vector x, N_Vector y);
SUNDIALS_EXPORT int OCTMatSpace_Sparse(SUNMatrix A, long int *lenrw, long int *leniw);


#ifdef __cplusplus
}
#endif


