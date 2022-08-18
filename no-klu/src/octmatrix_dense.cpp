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
  static int Matvec_DenseCSC(SUNMatrix A, N_Vector x, N_Vector y);
  static int Matvec_DenseCSR(SUNMatrix A, N_Vector x, N_Vector y);
  static int format_convert(const SUNMatrix A, SUNMatrix B);

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

    /* ----------------------------------------------------------------------------
 * Function to reallocate internal sparse matrix storage arrays so that the
 * resulting sparse matrix has storage for a specified number of nonzeros.
 * Returns 0 on success and 1 on failure (e.g. if A does not have sparse type,
 * or if nnz is negative)
 */

// int OCTDenseMatrix_Reallocate(SUNMatrix A, sunindextype NNZ)
// {
//   /* check for valid matrix type */
//   if (SUNMatGetID(A) != SUNMATRIX_DENSE)  return SUNMAT_ILL_INPUT;

//   /* check for valid nnz */
//   if (NNZ < 0)  return SUNMAT_ILL_INPUT;

//   /* perform reallocation */
//   // SM_INDEXVALS_S(A) = (sunindextype *) realloc(SM_INDEXVALS_S(A), NNZ*sizeof(sunindextype));
//   // SM_DATA_D(A) = (realtype *) realloc(SM_DATA_D(A), NNZ*sizeof(realtype));
//   // SM_NNZ_S(A) = NNZ;

//   Matrix *am;
//   am = static_cast <Matrix *> SM_CONTENT_D(A);

//   //Used change_capacity to directly remove memory
//   //that was occupied by zeros.
//   // But as documented this is an expensive approach,
//   // and also changes size of data arrays
//   am->change_capacity(NNZ);
//   return SUNMAT_SUCCESS;
// }

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

  // int SUNMatScaleAdd_Dense(realtype c, SUNMatrix A, SUNMatrix B)
  // {
  //   sunindextype j, i, p, nz, newvals, M, N, cend;
  //   booleantype newmat;
  //   sunindextype *w, *Ap, *Ai, *Bp, *Bi, *Cp, *Ci;
  //   realtype *x, *Ax, *Bx, *Cx;
  //   SUNMatrix C;

  //   /* Verify that A and B are compatible */
  //   if (!SMCompatible_Dense(A, B))
  //     return SUNMAT_ILL_INPUT;

  //   /* store shortcuts to matrix dimensions (M is inner dimension, N is outer) */
  //   if (SM_DenseTYPE_S(A) == CSC_MAT) {
  //     M = SM_ROWS_D(A);
  //     N = SM_COLUMNS_D(A);
  //   }
  //   else {
  //     M = SM_COLUMNS_D(A);
  //     N = SM_ROWS_D(A);
  //   }

  //   /* access data arrays from A and B (return if failure) */
  //   Ap = Ai = Bp = Bi = NULL;
  //   Ax = Bx = NULL;
  //   if (SM_INDEXPTRS_S(A))  Ap = SM_INDEXPTRS_S(A);
  //   else  return(SUNMAT_MEM_FAIL);
  //   if (SM_INDEXVALS_S(A))  Ai = SM_INDEXVALS_S(A);
  //   else  return(SUNMAT_MEM_FAIL);
  //   if (SM_DATA_D(A))       Ax = SM_DATA_D(A);
  //   else  return(SUNMAT_MEM_FAIL);
  //   if (SM_INDEXPTRS_S(B))  Bp = SM_INDEXPTRS_S(B);
  //   else  return(SUNMAT_MEM_FAIL);
  //   if (SM_INDEXVALS_S(B))  Bi = SM_INDEXVALS_S(B);
  //   else  return(SUNMAT_MEM_FAIL);
  //   if (SM_DATA_D(B))       Bx = SM_DATA_D(B);
  //   else  return(SUNMAT_MEM_FAIL);

  //   /* create work arrays for row indices and nonzero column values */
  //   w = (sunindextype *) malloc(M * sizeof(sunindextype));
  //   x = (realtype *) malloc(M * sizeof(realtype));

  //   /* determine if A already contains the sparsity pattern of B */
  //   newvals = 0;
  //   for (j=0; j<N; j++) {

  //     /* clear work array */
  //     for (i=0; i<M; i++)  w[i] = 0;

  //     /* scan column of A, incrementing w by one */
  //     for (i=Ap[j]; i<Ap[j+1]; i++)
  //       w[Ai[i]] += 1;

  //     /* scan column of B, decrementing w by one */
  //     for (i=Bp[j]; i<Bp[j+1]; i++)
  //       w[Bi[i]] -= 1;

  //     /* if any entry of w is negative, A doesn't contain B's sparsity,
  //        so increment necessary storage counter */
  //     for (i=0; i<M; i++)
  //       if (w[i] < 0)  newvals += 1;
  //   }

  //   /* If extra nonzeros required, check whether A has sufficient storage space
  //      for new nonzero entries (so B can be inserted into existing storage) */
  //   newmat = SUNFALSE;   /* no reallocation needed */
  //   if (newvals > (SM_NNZ_S(A) - Ap[N]))
  //     newmat = SUNTRUE;

  //   /* perform operation based on existing/necessary structure */

  //   /*   case 1: A already contains sparsity pattern of B */
  //   if (newvals == 0) {

  //     /* iterate through columns, adding matrices */
  //     for (j=0; j<N; j++) {

  //       /* clear work array */
  //       for (i=0; i<M; i++)
  //         x[i] = ZERO;

  //       /* scan column of B, updating work array */
  //       for (i = Bp[j]; i < Bp[j+1]; i++)
  //         x[Bi[i]] = Bx[i];

  //       /* scan column of A, updating array entries appropriately */
  //       for (i = Ap[j]; i < Ap[j+1]; i++)
  //         Ax[i] = c*Ax[i] + x[Ai[i]];

  //     }

  //   /*   case 2: A has sufficient storage, but does not already contain B's sparsity */
  //   } else if (!newmat) {

  //     /* determine storage location where last column (row) should end */
  //     nz = Ap[N] + newvals;

  //     /* store pointer past last column (row) from original A,
  //        and store updated value in revised A */
  //     cend = Ap[N];
  //     Ap[N] = nz;

  //     /* iterate through columns (rows) backwards */
  //     for (j=N-1; j>=0; j--) {

  //       /* clear out temporary arrays for this column (row) */
  //       for (i=0; i<M; i++) {
  //         w[i] = 0;
  //         x[i] = RCONST(0.0);
  //       }

  //       /* iterate down column (row) of A, collecting nonzeros */
  //       for (p=Ap[j]; p<cend; p++) {
  //         w[Ai[p]] += 1;         /* indicate that row (column) is filled */
  //         x[Ai[p]] = c*Ax[p];    /* collect/scale value */
  //       }

  //       /* iterate down column of B, collecting nonzeros */
  //       for (p=Bp[j]; p<Bp[j+1]; p++) {
  //         w[Bi[p]] += 1;       /* indicate that row is filled */
  //         x[Bi[p]] += Bx[p];   /* collect value */
  //       }

  //       /* fill entries of A with this column's (row's) data */
  //       for (i=M-1; i>=0; i--) {
  //         if ( w[i] > 0 ) {
  //           Ai[--nz] = i;
  //           Ax[nz] = x[i];
  //         }
  //       }

  //       /* store ptr past this col (row) from orig A, update value for new A */
  //       cend = Ap[j];
  //       Ap[j] = nz;

  //     }

  //   /*   case 3: A must be reallocated with sufficient storage */
  //   } else {

  //     /* create new matrix for sum */
  //     C = OCTDenseMatrix(SM_ROWS_D(A), SM_COLUMNS_D(A),
  //                         Ap[N] + newvals, SM_DenseTYPE_S(A), A->sunctx);

  //     /* access data from CSR structures (return if failure) */
  //     Cp = Ci = NULL;
  //     Cx = NULL;
  //     if (SM_INDEXPTRS_S(C))  Cp = SM_INDEXPTRS_S(C);
  //     else  return(SUNMAT_MEM_FAIL);
  //     if (SM_INDEXVALS_S(C))  Ci = SM_INDEXVALS_S(C);
  //     else  return(SUNMAT_MEM_FAIL);
  //     if (SM_DATA_D(C))       Cx = SM_DATA_D(C);
  //     else  return(SUNMAT_MEM_FAIL);

  //     /* initialize total nonzero count */
  //     nz = 0;

  //     /* iterate through columns (rows) */
  //     for (j=0; j<N; j++) {

  //       /* set current column (row) pointer to current # nonzeros */
  //       Cp[j] = nz;

  //       /* clear out temporary arrays for this column (row) */
  //       for (i=0; i<M; i++) {
  //         w[i] = 0;
  //         x[i] = RCONST(0.0);
  //       }

  //       /* iterate down column of A, collecting nonzeros */
  //       for (p=Ap[j]; p<Ap[j+1]; p++) {
  //         w[Ai[p]] += 1;         /* indicate that row is filled */
  //         x[Ai[p]] = c*Ax[p];    /* collect/scale value */
  //       }

  //       /* iterate down column of B, collecting nonzeros */
  //       for (p=Bp[j]; p<Bp[j+1]; p++) {
  //         w[Bi[p]] += 1;       /* indicate that row is filled */
  //         x[Bi[p]] += Bx[p];   /* collect value */
  //       }

  //       /* fill entries of C with this column's data */
  //       for (i=0; i<M; i++) {
  //         if ( w[i] > 0 ) {
  //           Ci[nz] = i;
  //           Cx[nz++] = x[i];
  //         }
  //       }
  //     }

  //     /* indicate end of data */
  //     Cp[N] = nz;

  //     /* update A's structure with C's values; nullify C's pointers */
  //     // SM_NNZ_S(A) = SM_NNZ_S(C);

  //     free(SM_DATA_D(A));
  //     // SM_DATA_D(A) = SM_DATA_D(C);
  //     // SM_DATA_D(C) = NULL;

  //     free(SM_INDEXVALS_S(A));
  //     // SM_INDEXVALS_S(A) = SM_INDEXVALS_S(C);
  //     // SM_INDEXVALS_S(C) = NULL;

  //     free(SM_INDEXPTRS_S(A));
  //     // SM_INDEXPTRS_S(A) = SM_INDEXPTRS_S(C);
  //     // SM_INDEXPTRS_S(C) = NULL;

  //     /* clean up */
  //     OCTMatDestroy_Dense(C);

  //   }

  //   /* clean up */
  //   free(w);
  //   free(x);

  //   /* return success */
  //   return(0);

  // }

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

  // /* -----------------------------------------------------------------
  //  * Computes y=A*x, where A is a CSC SUNMATRIX_DENSE of dimension MxN, x is a
  //  * compatible N_Vector object of length N, and y is a compatible
  //  * N_Vector object of length M.
  //  *
  //  * Returns 0 if successful, 1 if unsuccessful (failed memory access, or both
  //  * x and y are the same vector).
  //  */
  // int Matvec_DenseCSC(SUNMatrix A, N_Vector x, N_Vector y)
  // {
  //   sunindextype i, j;
  //   sunindextype *Ap, *Ai;
  //   realtype *Ax, *xd, *yd;

  //   /* access data from CSC structure (return if failure) */
  //   Ap = SM_INDEXPTRS_S(A);
  //   Ai = SM_INDEXVALS_S(A);
  //   Ax = SM_DATA_D(A);
  //   if ((Ap == NULL) || (Ai == NULL) || (Ax == NULL))
  //     return SUNMAT_MEM_FAIL;

  //   /* access vector data (return if failure) */
  //   xd = N_VGetArrayPointer(x);
  //   yd = N_VGetArrayPointer(y);
  //   if ((xd == NULL) || (yd == NULL) || (xd == yd) )
  //     return SUNMAT_MEM_FAIL;

  //   /* initialize result */
  //   for (i=0; i<SM_ROWS_D(A); i++)
  //     yd[i] = 0.0;

  //   /* iterate through matrix columns */
  //   for (j=0; j<SM_COLUMNS_D(A); j++) {

  //     /* iterate down column of A, performing product */
  //     for (i=Ap[j]; i<Ap[j+1]; i++)
  //       yd[Ai[i]] += Ax[i]*xd[j];

  //   }

  //   return SUNMAT_SUCCESS;
  // }

  // /* -----------------------------------------------------------------
  //  * Computes y=A*x, where A is a CSR SUNMATRIX_DENSE of dimension MxN, x is a
  //  * compatible N_Vector object of length N, and y is a compatible
  //  * N_Vector object of length M.
  //  *
  //  * Returns 0 if successful, -1 if unsuccessful (failed memory access).
  //  */
  // // int Matvec_DenseCSR(SUNMatrix A, N_Vector x, N_Vector y)
  // // {
  // //   sunindextype i, j;
  // //   sunindextype *Ap, *Aj;
  // //   realtype *Ax, *xd, *yd;

  // //   /* access data from CSR structure (return if failure) */
  // //   Ap = SM_INDEXPTRS_S(A);
  // //   Aj = SM_INDEXVALS_S(A);
  // //   Ax = SM_DATA_D(A);
  // //   if ((Ap == NULL) || (Aj == NULL) || (Ax == NULL))
  // //     return SUNMAT_MEM_FAIL;

  // //   /* access vector data (return if failure) */
  // //   xd = N_VGetArrayPointer(x);
  // //   yd = N_VGetArrayPointer(y);
  // //   if ((xd == NULL) || (yd == NULL) || (xd == yd))
  // //     return SUNMAT_MEM_FAIL;

  // //   /* initialize result */
  // //   for (i=0; i<SM_ROWS_D(A); i++)
  // //     yd[i] = 0.0;

  // //   /* iterate through matrix rows */
  // //   for (i=0; i<SM_ROWS_D(A); i++) {

  // //     /* iterate along row of A, performing product */
  // //     for (j=Ap[i]; j<Ap[i+1]; j++)
  // //       yd[i] += Ax[j]*xd[Aj[j]];

  // //   }

  // //   return(0);
  // // }

  // /* -----------------------------------------------------------------
  //  * Copies A into a matrix B in the opposite format of A.
  //  * Returns 0 if successful, nonzero if unsuccessful.
  //  */
  // // int format_convert(const SUNMatrix A, SUNMatrix B)
  // // {
  // //     realtype *Ax, *Bx;
  // //     sunindextype *Ap, *Aj;
  // //     sunindextype *Bp, *Bi;
  // //     sunindextype n_row, n_col, nnz;
  // //     sunindextype n, col, csum, row, last;

  // //     if (SM_DenseTYPE_S(A) == SM_DenseTYPE_S(B))
  // //       return OCTMatCopy_Dense(A, B);

  // //     Ap = SM_INDEXPTRS_S(A);
  // //     Aj = SM_INDEXVALS_S(A);
  // //     Ax = SM_DATA_D(A);

  // //     n_row = (SM_DenseTYPE_S(A) == CSR_MAT) ? SM_ROWS_D(A) : SM_COLUMNS_D(A);
  // //     n_col = (SM_DenseTYPE_S(A) == CSR_MAT) ? SM_COLUMNS_D(A) : SM_ROWS_D(A);

  // //     Bp = SM_INDEXPTRS_S(B);
  // //     Bi = SM_INDEXVALS_S(B);
  // //     Bx = SM_DATA_D(B);

  // //     nnz = Ap[n_row];

  // //     OCTMatZero_Dense(B);

  // //     /* compute number of non-zero entries per column (if CSR) or per row (if CSC) of A */
  // //     for (n = 0; n < nnz; n++)
  // //     {
  // //         Bp[Aj[n]]++;
  // //     }

  // //     /* cumualtive sum the nnz per column to get Bp[] */
  // //     for (col = 0, csum = 0; col < n_col; col++)
  // //     {
  // //         sunindextype temp  = Bp[col];
  // //         Bp[col] = csum;
  // //         csum += temp;
  // //     }
  // //     Bp[n_col] = nnz;

  // //     for (row = 0; row < n_row; row++)
  // //     {
  // //         sunindextype jj;
  // //         for (jj = Ap[row]; jj < Ap[row+1]; jj++)
  // //         {
  // //             sunindextype col  = Aj[jj];
  // //             sunindextype dest = Bp[col];

  // //             Bi[dest] = row;
  // //             Bx[dest] = Ax[jj];

  // //             Bp[col]++;
  // //         }
  // //     }

  // //     for (col = 0, last = 0; col <= n_col; col++)
  // //     {
  // //         sunindextype temp  = Bp[col];
  // //         Bp[col] = last;
  // //         last    = temp;
  // //     }

  // //     return 0;
  // // }

} // extern "C"