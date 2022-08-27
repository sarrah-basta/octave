/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner @ LLNL
 * Based on code sundials_Dense.c by: Carol Woodward and
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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

// #include <sunmatrix/SUNMATRIX_DENSE.h>
#include "octmatrix_dense.h"
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_math.h>

extern "C"
{

  using namespace octave;

  #define ZERO RCONST(0.0)
  #define ONE RCONST(1.0)

  /* Private function prototypes */
  static booleantype SMCompatible_Dense(SUNMatrix A, SUNMatrix B);
  static booleantype SMCompatible2_Dense(SUNMatrix A, N_Vector x, N_Vector y);
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
   * Function to create a new dense matrix
   */

  SUNMatrix OCTDenseMatrix(sunindextype M, sunindextype N,
                         SUNContext sunctx)
  {
    SUNMatrix A;
    void *content;

    /* return with NULL matrix on illegal input */
    if ((M <= 0) || (N <= 0))
      return (NULL);

    /* Create an empty matrix object */
    A = NULL;
    A = SUNMatNewEmpty(sunctx);
    if (A == NULL)
      return (NULL);

    // /* Attach operations */
    A->ops->getid = OCTMatGetID_Dense;
    A->ops->clone = OCTMatClone_Dense;
    A->ops->destroy = OCTMatDestroy_Dense;
    A->ops->zero = OCTMatZero_Dense;
    A->ops->copy = OCTMatCopy_Dense;
    A->ops->scaleadd = NULL;
    A->ops->scaleaddi = OCTMatScaleAddI_Dense;
    A->ops->matvec = OCTMatMatvec_Dense;
    A->ops->space = OCTMatSpace_Dense;

    /* Create content */
    content = NULL;
    // creating a new DenseMatrix Octave data structure with rows and columns as specified
    content = new Matrix(M, N);
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

  sunindextype OCTDenseMatrix_Rows(SUNMatrix A)
  {
    if (SUNMatGetID(A) == SUNMATRIX_DENSE)
      return SM_ROWS_D(A);
    else
      return SUNMAT_ILL_INPUT;
  }

  sunindextype OCTDenseMatrix_Columns(SUNMatrix A)
  {
    if (SUNMatGetID(A) == SUNMATRIX_DENSE)
      return SM_COLUMNS_D(A);
    else
      return SUNMAT_ILL_INPUT;
  }

  realtype *OCTDenseMatrix_Column(SUNMatrix A, sunindextype j)
  {
    if (SUNMatGetID(A) == SUNMATRIX_DENSE)
      return ((realtype *) SM_COLUMN_D(A,j));
    else
      return NULL;
  }

  realtype *OCTDenseMatrix_Data(SUNMatrix A)
  {
    if (SUNMatGetID(A) == SUNMATRIX_DENSE)
      return SM_DATA_D(A);
    else
      return NULL;
  }

  sunindextype OCTDenseMatrix_LData(SUNMatrix A)
  {
    if (SUNMatGetID(A) == SUNMATRIX_DENSE)
      return SM_LDATA_D(A);
    else
      return NULL;
  }

  // Function required to test

  /* ----------------------------------------------------------------------------
   * Function to create a new sparse matrix from an existing dense matrix
   * by copying all nonzero values into the sparse matrix structure.  Returns NULL
   * if the request for matrix storage cannot be satisfied.
   */

  /*
   * -----------------------------------------------------------------
   * implementation of matrix operations
   * -----------------------------------------------------------------
   */

  SUNMatrix_ID OCTMatGetID_Dense(SUNMatrix A)
  {
    return SUNMATRIX_DENSE;
  }

  SUNMatrix OCTMatClone_Dense(SUNMatrix A)
  {
    SUNMatrix B = OCTDenseMatrix(SM_ROWS_D(A), SM_COLUMNS_D(A),
                                A->sunctx);
    return (B);
  }

  void OCTMatDestroy_Dense(SUNMatrix A)
  {
    if (A == NULL)
      return;

    //   /* free content */
    if (A->content != NULL)
    {
      /* free content struct */
      delete A->content;
      A->content = NULL;
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

  int OCTMatZero_Dense(SUNMatrix A)
  {
    // sunindextype i;
    Matrix *am;
    am = static_cast <Matrix *> SM_CONTENT_D(A);
    sunindextype i,j;
    //Used change_capacity to directly remove memory
    //that was occupied by zeros.
    // But as documented this is an expensive approach,
    // and also changes size of data arrays
    // am->change_capacity(0);

    // So should naive approach to replace
    // with zeros be used instead
    /* Perform operation */
    for(i=0; i < SM_ROWS_D(A); i++ ){
      for(j=0; j<SM_COLUMNS_D(A); j++){
        (*am)((j-1)*i + 1) = 0;
      }
    }
    return SUNMAT_SUCCESS;
  }


  int OCTMatSpace_Dense(SUNMatrix A, long int *lenrw, long int *leniw)
  {
    *lenrw = SM_LDATA_D(A);
    *leniw = 3 + SM_COLUMNS_D(A);
    return SUNMAT_SUCCESS;
  }

  int OCTMatCopy_Dense(SUNMatrix A, SUNMatrix B)
  {
    sunindextype i,j;

    /* Verify that A and B are compatible */
    if (!SMCompatible_Dense(A, B))
      return SUNMAT_ILL_INPUT;

    /* Perform operation */
    for (j=0; j<SM_COLUMNS_D(A); j++)
      for (i=0; i<SM_ROWS_D(A); i++)
        SM_ELEMENT_D(B,i,j) = SM_ELEMENT_D(A,i,j);

    return SUNMAT_SUCCESS;
  }

  // Debugging Functions
  /* ----------------------------------------------------------------------------
   * Function to print the sparse matrix
   */

  void OCTDenseMatrix_Print(SUNMatrix A, FILE *outfile)
  {
    sunindextype i, j;
    char *matrixtype;
    char *indexname;

    /* should not be called unless A is a sparse matrix;
       otherwise return immediately */
    if (SUNMatGetID(A) != SUNMATRIX_DENSE)
      return;

    std::cout<<A;
    return;
  }

  // optional matrix functions

  int OCTMatScaleAddI_Dense(realtype c, SUNMatrix A)
  {
    sunindextype i, j;

    /* Perform operation */
    for (j=0; j<SM_COLUMNS_D(A); j++)
      for (i=0; i<SM_ROWS_D(A); i++) {
        SM_ELEMENT_D(A,i,j) *= c;
        if (i == j)
          SM_ELEMENT_D(A,i,j) += ONE;
      }
    return SUNMAT_SUCCESS;

  }

  //only needed temporarily for test function
  int OCTMatMatvec_Dense(SUNMatrix A, N_Vector x, N_Vector y)
  {
    sunindextype i, j;
    realtype *col_j, *xd, *yd;

    /* Verify that A, x and y are compatible */
    if (!SMCompatible2_Dense(A, x, y))
      return SUNMAT_ILL_INPUT;

    /* access vector data (return if failure) */
    xd = N_VGetArrayPointer(x);
    yd = N_VGetArrayPointer(y);
    if ((xd == NULL) || (yd == NULL) || (xd == yd))
      return SUNMAT_MEM_FAIL;

    /* Perform operation */
    for (i=0; i<SM_ROWS_D(A); i++)
      yd[i] = ZERO;
    for(j=0; j<SM_COLUMNS_D(A); j++) {
      // col_j = (*SM_COLUMNS_D(A))(j);
      for (i=0; i<SM_ROWS_D(A); i++)
        yd[i] += SM_ELEMENT_D(A,i,j)*xd[j];
    }
    return SUNMAT_SUCCESS;
  }

  /*
   * =================================================================
   * private functions
   * =================================================================
   */
  static booleantype SMCompatible_Dense(SUNMatrix A, SUNMatrix B)
{
  /* both matrices must be SUNMATRIX_DENSE */
  if (SUNMatGetID(A) != SUNMATRIX_DENSE)
    return SUNFALSE;
  if (SUNMatGetID(B) != SUNMATRIX_DENSE)
    return SUNFALSE;

  /* both matrices must have the same shape */
  if (SM_ROWS_D(A) != SM_ROWS_D(B))
    return SUNFALSE;
  if (SM_COLUMNS_D(A) != SM_COLUMNS_D(B))
    return SUNFALSE;

  return SUNTRUE;
}

  /* -----------------------------------------------------------------
   * Function to check compatibility of a SUNMatrix object with two
   * N_Vectors (A*x = b)
   */

  static booleantype SMCompatible2_Dense(SUNMatrix A, N_Vector x, N_Vector y)
  {
    /* vectors must implement N_VGetArrayPointer */
    if ((x->ops->nvgetarraypointer == NULL) ||
        (y->ops->nvgetarraypointer == NULL))
      return SUNFALSE;

    /* Verify that the dimensions of A, x, and y agree */
    if ((OCTDenseMatrix_Columns(A) != N_VGetLength(x)) ||
        (OCTDenseMatrix_Rows(A) != N_VGetLength(y)))
      return SUNFALSE;

    return SUNTRUE;
  }

} // extern "C"