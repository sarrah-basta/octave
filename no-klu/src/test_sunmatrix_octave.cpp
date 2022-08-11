/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner @ LLNL
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
 * This is the testing routine to check the SUNMatrix Sparse module
 * implementation.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
// #include <sunmatrix/sunmatrix_sparse.h>
#include "octmatrix_sparse.h"
// #include <nvector/nvector_Octave.h>
#include "nvector_octave.h"
#include <sundials/sundials_math.h>
#include "test_sunmatrix.h"


/* ----------------------------------------------------------------------
 * Main SUNMatrix Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int          fails=0;                    /* counter for test failures  */
  sunindextype matrows, matcols;           /* matrix dims                */
  int          mattype;                    /* matrix storage type        */
  N_Vector     x, y, z;                    /* test vectors               */
  realtype*    vecdata;                    /* pointers to vector data    */
  SUNMatrix    A, B, C, D, I;              /* test matrices              */
  realtype*    matdata;                    /* pointer to matrix data     */
  sunindextype i, j, k, kstart, kend, N, uband, lband;
  sunindextype *colptrs, *rowindices;
  sunindextype *rowptrs, *colindices;
  int          print_timing, square;
  SUNContext   sunctx;

  if (SUNContext_Create(NULL, &sunctx)) {
    printf("ERROR: SUNContext_Create failed\n");
    return(-1);
  }

  /* check input and set vector length */
  // if (argc < 5){
  //   printf("ERROR: FOUR (4) Input required: matrix rows, matrix cols, matrix type (0/1), print timing \n");
  //   return(-1);
  // }

  // matrows = (sunindextype) atol(argv[1]);
  // if (matrows < 1) {
  //   printf("ERROR: number of rows must be a positive integer\n");
  //   return(-1);
  // }
  matrows = 100;

  // matcols = (sunindextype) atol(argv[2]);
  // if (matcols < 1) {
  //   printf("ERROR: number of cols must be a positive integer\n");
  //   return(-1);
  // }
  matcols = 100;

  // k = (sunindextype) atol(argv[3]);
  // if ((k != 0) && (k != 1)) {
  //   printf("ERROR: matrix type must be 0 or 1\n");
  //   return(-1);
  // }
  // mattype = (k == 0) ? CSC_MAT : CSR_MAT;
  mattype = CSC_MAT;

  // print_timing = atoi(argv[4]);
  print_timing = 1;
  SetTiming(print_timing);

  square = (matrows == matcols) ? 1 : 0;
  printf("\nSparse matrix test: size %ld by %ld, type = %idone, %ld\n",
         (long int) matrows, (long int) matcols, mattype, 6758080645);


  /* Initialize vectors and matrices to NULL */
  x = NULL;
  y = NULL;
  z = NULL;
  A = NULL;
  B = NULL;
  C = NULL;
  D = NULL;
  I = NULL;

  /* check creating sparse matrix from dense matrix */
  B = SUNDenseMatrix(5,6, sunctx);

  matdata = SUNDenseMatrix_Data(B);
  matdata[2]  = RCONST(1.0);    /* [ 0 2 0 0 7 0 ] */
  matdata[5]  = RCONST(2.0);    /* [ 0 0 4 0 8 0 ] */
  matdata[9]  = RCONST(3.0);    /* [ 1 0 0 0 0 0 ] */
  matdata[11] = RCONST(4.0);    /* [ 0 0 5 6 0 0 ] */
  matdata[13] = RCONST(5.0);    /* [ 0 3 0 0 0 9 ] */
  matdata[18] = RCONST(6.0);
  matdata[20] = RCONST(7.0);
  matdata[21] = RCONST(8.0);
  matdata[29] = RCONST(9.0);

  if (mattype == CSR_MAT) {

    /* Check CSR */
    C = OCTSparseMatrix(5, 6, 9, CSR_MAT, sunctx);
    rowptrs = OCTSparseMatrix_IndexPointers(C);
    colindices = OCTSparseMatrix_IndexValues(C);
    matdata = OCTSparseMatrix_Data(C);
    rowptrs[0] = 0;
    matdata[0] = RCONST(2.0);   colindices[0] = 1;
    matdata[1] = RCONST(7.0);   colindices[1] = 4;
    rowptrs[1] = 2;
    matdata[2] = RCONST(4.0);   colindices[2] = 2;
    matdata[3] = RCONST(8.0);   colindices[3] = 4;
    rowptrs[2] = 4;
    matdata[4] = RCONST(1.0);   colindices[4] = 0;
    rowptrs[3] = 5;
    matdata[5] = RCONST(5.0);   colindices[5] = 2;
    matdata[6] = RCONST(6.0);   colindices[6] = 3;
    rowptrs[4] = 7;
    matdata[7] = RCONST(3.0);   colindices[7] = 1;
    matdata[8] = RCONST(9.0);   colindices[8] = 5;
    rowptrs[5] = 9;

    A = OCTSparseFromDenseMatrix(B, ZERO, CSR_MAT);
    fails += check_matrix(A, C, 1e-15);

    if (fails) {
      printf("FAIL: SUNMatrix SparseFromDense CSR conversion failed\n");
      return(1);
    }

    SUNMatDestroy(A);
    SUNMatDestroy(C);

  } else {

    /* Check CSC */
    D = OCTSparseMatrix(5, 6, 9, CSC_MAT, sunctx);
    colptrs = OCTSparseMatrix_IndexPointers(D);
    rowindices = OCTSparseMatrix_IndexValues(D);
    matdata = OCTSparseMatrix_Data(D);
    colptrs[0] = 0;
    matdata[0] = RCONST(1.0);   rowindices[0] = 2;
    colptrs[1] = 1;
    matdata[1] = RCONST(2.0);   rowindices[1] = 0;
    matdata[2] = RCONST(3.0);   rowindices[2] = 4;
    colptrs[2] = 3;
    matdata[3] = RCONST(4.0);   rowindices[3] = 1;
    matdata[4] = RCONST(5.0);   rowindices[4] = 3;
    colptrs[3] = 5;
    matdata[5] = RCONST(6.0);   rowindices[5] = 3;
    colptrs[4] = 6;
    matdata[6] = RCONST(7.0);   rowindices[6] = 0;
    matdata[7] = RCONST(8.0);   rowindices[7] = 1;
    colptrs[5] = 8;
    matdata[8] = RCONST(9.0);   rowindices[8] = 4;
    colptrs[6] = 9;

    A = OCTSparseFromDenseMatrix(B, 1e-15, CSC_MAT);
    fails += check_matrix(A, D, 1e-15);

    if (fails) {
      printf("FAIL: SUNMatrix SparseFromDense CSC conversion failed\n");
      return(1);
    }

    SUNMatDestroy(A);
    SUNMatDestroy(D);

  }
  SUNMatDestroy(B);


  /* Create/fill I matrix */
  I = NULL;
  if (square) {
    I = OCTSparseMatrix(matrows, matcols, matcols, mattype, sunctx);
    matdata    = OCTSparseMatrix_Data(I);
    colindices = OCTSparseMatrix_IndexValues(I);
    rowptrs    = OCTSparseMatrix_IndexPointers(I);
    for(i=0; i<matrows; i++) {
      matdata[i] = ONE;
      colindices[i] = i;
      rowptrs[i] = i;
    }
    rowptrs[matrows] = matrows;
  }

  /* Create/fill random dense matrices, create sparse from them */
  C = SUNDenseMatrix(matrows, matcols, sunctx);
  D = SUNDenseMatrix(matrows, matcols, sunctx);
  for (k=0; k<3*matrows; k++) {
    i = rand() % matrows;
    j = rand() % matcols;
    matdata = SUNDenseMatrix_Column(D,j);
    matdata[i] = (realtype) rand() / (realtype) RAND_MAX;
  }
  for (k=0; k<matrows; k++) {
    i = rand() % matrows;
    j = rand() % matcols;
    matdata = SUNDenseMatrix_Column(C,j);
    matdata[i] = (realtype) rand() / (realtype) RAND_MAX;
  }
  A = OCTSparseFromDenseMatrix(C, ZERO, mattype);
  B = OCTSparseFromDenseMatrix(D, ZERO, mattype);

  /* Create vectors and fill */
  x = N_VNew_Octave(matcols, sunctx);
  y = N_VNew_Octave(matrows, sunctx);
  z = N_VNew_Octave(matrows, sunctx);
  vecdata = N_VGetArrayPointer(x);
  for(i=0; i<matcols; i++)
    vecdata[i] = (realtype) rand() / (realtype) RAND_MAX;
  if (SUNMatMatvec(C, x, y) != 0) {
    printf("FAIL: SUNMatrix module Dense matvec failure \n \n");
    SUNMatDestroy(A);  SUNMatDestroy(B);
    SUNMatDestroy(C);  SUNMatDestroy(D);
    N_VDestroy(x);  N_VDestroy(y);  N_VDestroy(z);
    if (square)
      SUNMatDestroy(I);
    return(1);
  }
  if (SUNMatMatvec(D, x, z) != 0) {
    printf("FAIL: SUNMatrix module Dense matvec failure \n \n");
    SUNMatDestroy(A);  SUNMatDestroy(B);
    SUNMatDestroy(C);  SUNMatDestroy(D);
    N_VDestroy(x);  N_VDestroy(y);  N_VDestroy(z);
    if (square)
      SUNMatDestroy(I);
    return(1);
  }

  /* SUNMatrix Tests */
  fails += Test_SUNMatGetID(A, SUNMATRIX_SPARSE, 0);
  fails += Test_SUNMatClone(A, 0);
  fails += Test_SUNMatCopy(A, 0);
  fails += Test_SUNMatZero(A, 0);
  // fails += Test_SUNMatScaleAdd(A, I, 0);
  // fails += Test_SUNMatScaleAdd2(A, B, x, y, z);
  // if (square) {
  //   fails += Test_SUNMatScaleAddI(A, I, 0);
  //   fails += Test_SUNMatScaleAddI2(A, x, y);
  // }
  // fails += Test_SUNMatMatvec(A, x, y, 0);
  fails += Test_SUNMatSpace(A, 0);
  // if (mattype == CSR_MAT) {
  //   fails += Test_OCTSparseMatrixToCSC(A);
  // } else {
  //   fails += Test_OCTSparseMatrixToCSR(A);
  // }

  /* Print result */
  if (fails) {
    printf("FAIL: SUNMatrix module failed %i tests \n \n", fails);
    printf("\nA =\n");
    OCTSparseMatrix_Print(A,stdout);
    printf("\nB =\n");
    OCTSparseMatrix_Print(B,stdout);
    if (square) {
      printf("\nI =\n");
      OCTSparseMatrix_Print(I,stdout);
    }
    // printf("\nx =\n");
    // N_VPrint_Octave(x);
    // printf("\ny =\n");
    // N_VPrint_Octave(y);
    // printf("\nz =\n");
    // N_VPrint_Octave(z);
  } else {
    printf("SUCCESS: SUNMatrix module passed all tests \n \n");
  }

  /* Free vectors and matrices */
  N_VDestroy(x);
  N_VDestroy(y);
  N_VDestroy(z);
  SUNMatDestroy(A);
  SUNMatDestroy(B);
  SUNMatDestroy(C);
  SUNMatDestroy(D);
  if (square)
    SUNMatDestroy(I);

  SUNContext_Free(&sunctx);

  return(fails);

  return(0);
}


/* ----------------------------------------------------------------------
 * Check matrix
 * --------------------------------------------------------------------*/
int check_matrix(SUNMatrix A, SUNMatrix B, realtype tol)
{
  int failure = 0;
  realtype *Adata, *Bdata;
  sunindextype *Aindexptrs, *Bindexptrs;
  sunindextype *Aindexvals, *Bindexvals;
  sunindextype i, ANP, BNP, Annz, Bnnz;

  /* get matrix pointers */
  Adata = OCTSparseMatrix_Data(A);
  Aindexptrs = OCTSparseMatrix_IndexPointers(A);
  Aindexvals = OCTSparseMatrix_IndexValues(A);
  ANP = OCTSparseMatrix_NP(A);
  Annz = Aindexptrs[ANP];

  Bdata = OCTSparseMatrix_Data(B);
  Bindexptrs = OCTSparseMatrix_IndexPointers(B);
  Bindexvals = OCTSparseMatrix_IndexValues(B);
  BNP = OCTSparseMatrix_NP(B);
  Bnnz = Bindexptrs[BNP];

  /* matrices must have same sparsetype, shape and actual data lengths */
  if (SUNMatGetID(A) != SUNMatGetID(B)) {
    printf(">>> ERROR: check_matrix: Different storage types (%d vs %d)\n",
           SUNMatGetID(A), SUNMatGetID(B));
    return(1);
  }
  if (OCTSparseMatrix_SparseType(A) != OCTSparseMatrix_SparseType(B)) {
    printf(">>> ERROR: check_matrix: Different storage types (%d vs %d)\n",
           OCTSparseMatrix_SparseType(A), OCTSparseMatrix_SparseType(B));
    return(1);
  }
  if (OCTSparseMatrix_Rows(A) != OCTSparseMatrix_Rows(B)) {
    printf(">>> ERROR: check_matrix: Different numbers of rows (%ld vs %ld)\n",
           (long int) OCTSparseMatrix_Rows(A), (long int) OCTSparseMatrix_Rows(B));
    return(1);
  }
  if (OCTSparseMatrix_Columns(A) != OCTSparseMatrix_Columns(B)) {
    printf(">>> ERROR: check_matrix: Different numbers of columns (%ld vs %ld)\n",
           (long int) OCTSparseMatrix_Columns(A),
           (long int) OCTSparseMatrix_Columns(B));
    return(1);
  }
  if (Annz != Bnnz) {
    printf(">>> ERROR: check_matrix: Different numbers of nonzeos (%ld vs %ld)\n",
           (long int) Annz, (long int) Bnnz);
    return(1);
  }

  /* compare sparsity patterns */
  for (i=0; i<ANP; i++)
    failure += (Aindexptrs[i] != Bindexptrs[i]);
  if (failure > ZERO) {
    printf(">>> ERROR: check_matrix: Different indexptrs \n");
    return(1);
  }
  for (i=0; i<Annz; i++)
    failure += (Aindexvals[i] != Bindexvals[i]);
  if (failure > ZERO) {
    printf(">>> ERROR: check_matrix: Different indexvals \n");
    return(1);
  }

  /* compare matrix values */
  for(i=0; i<Annz; i++)
    failure += SUNRCompareTol(Adata[i], Bdata[i], tol);
  if (failure > ZERO) {
    printf(">>> ERROR: check_matrix: Different entries \n");
    return(1);
  }

  return(0);
}

int check_matrix_entry(SUNMatrix A, realtype val, realtype tol)
{
  int failure = 0;
  realtype *Adata;
  sunindextype *indexptrs;
  sunindextype i, NP;

  /* get data pointer */
  Adata = OCTSparseMatrix_Data(A);

  /* compare data */
  indexptrs = OCTSparseMatrix_IndexPointers(A);
  NP = OCTSparseMatrix_NP(A);
  for(i=0; i < indexptrs[NP]; i++){
    failure += SUNRCompareTol(Adata[i], val, tol);
  }

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

int check_vector(N_Vector x, N_Vector y, realtype tol)
{
  int failure = 0;
  realtype *xdata, *ydata;
  sunindextype xldata, yldata;
  sunindextype i;

  /* get vector data */
  xdata = N_VGetArrayPointer(x);
  ydata = N_VGetArrayPointer(y);

  /* check data lengths */
  xldata = N_VGetLength_Octave(x);
  yldata = N_VGetLength_Octave(y);

  if (xldata != yldata) {
    printf(">>> ERROR: check_vector: Different data array lengths \n");
    return(1);
  }

  /* check vector data */
  for(i=0; i < xldata; i++){
    failure += SUNRCompareTol(xdata[i], ydata[i], tol);
  }

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

booleantype has_data(SUNMatrix A)
{
  realtype *Adata = OCTSparseMatrix_Data(A);
  if (Adata == NULL)
    return SUNFALSE;
  else
    return SUNTRUE;
}

booleantype is_square(SUNMatrix A)
{
  if (OCTSparseMatrix_Rows(A) == OCTSparseMatrix_Columns(A))
    return SUNTRUE;
  else
    return SUNFALSE;
}

void sync_device(SUNMatrix A)
{
  /* not running on GPU, just return */
  return;
}
