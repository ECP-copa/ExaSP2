/// \file
/// Sparse data structure

#ifndef __SPARSEMATRIX_H
#define __SPARSEMATRIX_H

#include <stdio.h>

#include "decomposition.h"
#include "mytype.h"

/// Matrix type
typedef enum {
  DENSE,
  SPARSE} MatrixType;

/// Sparse Matrix struct
///
///
typedef struct DataMatrixSt
{
  MatrixType mtype; //!< matrix type (dense/sparse)
  int hsize;        //!< denotes size of dense matrix would be hsize x hsize
  int msize;        //!< max number of non-zeroes per row
  int bandwidth;    //!< bandwidth of matrix
  int* iia;         //!< ELLPACK-R vector containing number of non-zero elements per row
#ifdef CONTIG_MATRIX
  int*  jjcontig;      //!< 1-D array mapping to jja
#endif
  int** jja;           //!< ELLPACK-R 2-D array containing indeces of values per row

  real_t maxEval;      //!< max Gershgorin limit
  real_t minEval;      //!< min Gershgorin limit
  real_t maxMinusMin;  //!< max-min difference
#ifdef CONTIG_MATRIX
  real_t*  valcontig;  //!< 1-D array mapping to val
#endif
  real_t** val;        //!< ELLPACK-R 2-D array containing non-zero elements by row
} DataMatrix;

// Allocate/deallocate matrices
DataMatrix* initDataMatrix(const MatrixType matType, int hsize, int msize);
void destroyDataMatrix(struct DataMatrixSt* dataMatrix);

// Other functions
int nnzStart(int hsize, int msize);
void sparsity(struct DataMatrixSt* dataMatrix);
void gershgorin(struct DataMatrixSt* dataMatrix, struct DomainSt* domain);
void normalize(struct DataMatrixSt* dataMatrix);
void trace(struct DataMatrixSt* dataMatrix, struct DomainSt* domain, real_t* tr, real_t* tr2);

#endif
