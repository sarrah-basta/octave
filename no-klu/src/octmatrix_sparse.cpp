/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner @ LLNL
 * Based on code sundials_sparse.c by: Carol Woodward and
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
 * This is the implementation file for the sparse implementation of
 * the SUNMATRIX package.
 * -----------------------------------------------------------------
 */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

// #include <sunmatrix/SUNMATRIX_CUSTOM.h>
#include "octmatrix_sparse.h"
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_math.h>

extern "C"
{

  using namespace octave;

  #define ZERO RCONST(0.0)
  #define ONE RCONST(1.0)

  /*
   * -----------------------------------------------------------------
   * exported functions
   * -----------------------------------------------------------------
   */

  /*
   * ==================================================================
   * Private function prototypes (functions working on SlsMat)
   * ==================================================================
   */

  /* ----------------------------------------------------------------------------
   * Function to create a new sparse matrix
   */

  SUNMatrix OCTSparseMatrix(sunindextype M, sunindextype N,
                            sunindextype NNZ, int sparsetype,
                            SUNContext sunctx)
  {
    SUNMatrix A;
    void *content;

    /* return with NULL matrix on illegal input */
    if ((M <= 0) || (N <= 0) || (NNZ < 0))
      return (NULL);
    if ((sparsetype != CSC_MAT) && (sparsetype != CSR_MAT))
      return (NULL);

    /* Create an empty matrix object */
    A = NULL;
    A = SUNMatNewEmpty(sunctx);
    if (A == NULL)
      return (NULL);

    // /* Attach operations */
    A->ops->getid = OCTMatGetID_Sparse;
    A->ops->clone = NULL;
    A->ops->destroy = OCTMatDestroy_Sparse;
    A->ops->zero = OCTMatZero_Sparse;
    A->ops->copy = NULL;
    A->ops->scaleadd = NULL;
    A->ops->scaleaddi = NULL;
    A->ops->matvec = NULL;
    A->ops->space = NULL;

    /* Create content */
    content = NULL;
    // creating a new SparseMatrix Octave data structure with rows and columns as specified
    content = new SparseMatrix(M, N, NNZ);
    if (content == NULL)
    {
      SUNMatDestroy(A);
      return (NULL);
    }

    /* Attach content */
    A->content = content;
        //   /* Fill content */
        return (A);
  }

  /* ----------------------------------------------------------------------------
   * Function to create a serial N_Vector with existing ColumnVector component
   */

  SUNMatrix OCTSparseMatrix_Make(SparseMatrix sm, SUNContext sunctx)
  {
    SUNMatrix A;
    void *content;

    /* Create an empty matrix object */
    A = NULL;
    A = SUNMatNewEmpty(sunctx);
    if (A == NULL)
      return (NULL);

    // /* Attach operations */
    A->ops->getid = OCTMatGetID_Sparse;
    A->ops->clone = NULL;
    A->ops->destroy = OCTMatDestroy_Sparse;
    A->ops->zero = OCTMatZero_Sparse;
    A->ops->copy = NULL;
    A->ops->scaleadd = NULL;
    A->ops->scaleaddi = NULL;
    A->ops->matvec = NULL;
    A->ops->space = NULL;

    /* Create content */
    content = NULL;
    // creating a new SparseMatrix Octave data structure with rows and columns as specified
    content = &sm;
    if (content == NULL)
    {
      SUNMatDestroy(A);
      return (NULL);
    }

    /* Attach content */
    A->content = content;

    //   /* Fill content */
    return (A);
  }
  /* ----------------------------------------------------------------------------
   * Functions to access the contents of the sparse matrix structure
   */

  sunindextype OCTSparseMatrix_Rows(SUNMatrix A)
  {
    if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
      return SM_ROWS_O(A);
    else
      return SUNMAT_ILL_INPUT;
  }

  sunindextype OCTSparseMatrix_Columns(SUNMatrix A)
  {
    if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
      return SM_COLUMNS_O(A);
    else
      return SUNMAT_ILL_INPUT;
  }

  sunindextype OCTSparseMatrix_NNZ(SUNMatrix A)
  {
    if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
      return SM_NNZ_O(A);
    else
      return SUNMAT_ILL_INPUT;
  }

  sunindextype OCTSparseMatrix_NP(SUNMatrix A)
  {
    if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
      return SM_NP_O(A);
    else
      return SUNMAT_ILL_INPUT;
  }

  int OCTSparseMatrix_SparseType(SUNMatrix A)
  {
    if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
      return SM_SPARSETYPE_O(A);
    else
      return SUNMAT_ILL_INPUT;
  }

  realtype *OCTSparseMatrix_Data(SUNMatrix A)
  {
    if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
      return SM_DATA_O(A);
    else
      return NULL;
  }

  sunindextype *OCTSparseMatrix_IndexValues(SUNMatrix A)
  {
    if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
      return SM_INDEXVALS_O(A);
    else
      return NULL;
  }

  sunindextype *OCTSparseMatrix_IndexPointers(SUNMatrix A)
  {
    if (SUNMatGetID(A) == SUNMATRIX_CUSTOM)
      return SM_INDEXPTRS_O(A);
    else
      return NULL;
  }

  // /*
  //  * -----------------------------------------------------------------
  //  * implementation of matrix operations
  //  * -----------------------------------------------------------------
  //  */

  SUNMatrix_ID OCTMatGetID_Sparse(SUNMatrix A)
  {
    return SUNMATRIX_CUSTOM;
  }

  void OCTMatDestroy_Sparse(SUNMatrix A)
  {
    if (A == NULL)
      return;

    //   /* free content */
    if (A->content != NULL)
    {
      /* deleting object created using new */
      delete(A->content);
    }

    /* free ops and matrix */
    if (A->ops)
    {
      free(A->ops);
      A->ops = NULL;
    }
    free(A);
    A = NULL;

    return;
  }

  int OCTMatZero_Sparse(SUNMatrix A)
  {
    // sunindextype i;
    SparseMatrix *am;
    SparseMatrix *ptr = new SparseMatrix(SM_ROWS_O(A),SM_COLS_O(A));
    // SparseMatrix *ptr = &nm;
    am = static_cast <SparseMatrix *> SM_CONTENT_O(A);
    am = ptr;

    return SUNMAT_SUCCESS;
  }

  // Debugging Functions
  /* ----------------------------------------------------------------------------
   * Function to print the sparse matrix
   */

  void OCTSparseMatrix_Print(SUNMatrix A)
  {
    SparseMatrix *am = static_cast <SparseMatrix *>SM_CONTENT_O(A);
    std::cout<<(*am);
    return;
  }

} // extern "C"