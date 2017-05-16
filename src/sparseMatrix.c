/// \File
/// Routines for creating and manipulating sparse matrices.

#ifdef MMATH_SPARSE

#include "sparseMatrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "performance.h"
#include "parallel.h"
#include "constants.h"

/// \details
/// Adjust number of non-zeroes
int nnzStart(int hsize, int msize)
{
  int M = msize;
  if (M == 0) M = hsize;
  if ((M % 32) > 0) M += (32 - (M % 32));
  if (M > hsize) M = hsize;
  if (printRank()) printf("Adjusted M = %d\n", M);

  return M;
}

/// \details
/// Allocate space for sparse matrix
DataMatrix* initDataMatrix(const MatrixType mtype, int hsize, int msize)
{
  if (mtype == DENSE)
  {
    printf("This DataMatrix implementation does not support dense\n");
    exit;
  }
  else if (mtype == SPARSE)
  {
    DataMatrix* dataMatrix = (DataMatrix*)malloc(sizeof(DataMatrix));

    // hsize is number of rows and msize is the max number of non-zeroes per row
    dataMatrix->mtype = mtype;
    dataMatrix->hsize = hsize;
    dataMatrix->msize = msize;
  
    // iia holds the number of non-zeroes in each row
    dataMatrix->iia = (int*)malloc(hsize*sizeof(int));

#ifdef CONTIG_MATRIX
    dataMatrix->jjcontig = (int*)malloc(hsize*msize*sizeof(int));
    dataMatrix->jja = (int**)malloc(hsize*sizeof(int*));
#pragma omp parallel for
    for (int i = 0; i < hsize; i++)
    {
      dataMatrix->jja[i] = &(dataMatrix->jjcontig[i*msize]);
    }

    // Zero counts of non-zeroes per row and indices
    memset(dataMatrix->jjcontig, 0, hsize*msize*sizeof(int));
  
    dataMatrix->valcontig = (real_t*)malloc(hsize*msize*sizeof(real_t));
    dataMatrix->val = (real_t**)malloc(hsize*sizeof(real_t*));

#pragma omp parallel for
    for (int i = 0; i < hsize; i++)
    {
      dataMatrix->val[i] = &(dataMatrix->valcontig[i*msize]);
    }
  
    // Zero non-zero values
    memset(dataMatrix->valcontig, ZERO, hsize*msize*sizeof(real_t));
 
#else

    // jja contains the column index for each non-zero value
    dataMatrix->jja = (int**)malloc(hsize*sizeof(int*));
    for (int i = 0; i < hsize; i++)
    {
      dataMatrix->jja[i] = (int*)malloc(msize*sizeof(int)); 
    }

    // val contains the non-zeroes
    dataMatrix->val = (real_t**)malloc(hsize*sizeof(real_t*));
    for (int i = 0; i < hsize; i++)
    {
      dataMatrix->val[i] = (real_t*)malloc(msize*sizeof(real_t));
    }
#endif

    // Zero counts of non-zeroes per row
    memset(dataMatrix->iia, 0, hsize*sizeof(int));

    // Used for normalization
    dataMatrix->maxEval = ZERO;
    dataMatrix->minEval = ZERO;
    dataMatrix->maxMinusMin = ZERO;

    // Matrix bandwidth
    dataMatrix->bandwidth = 0;

    return dataMatrix;
  }
}

/// \details
/// Deallocate space for sparse matrix
void destroyDataMatrix(struct DataMatrixSt* dataMatrix)
{
  int hsize = dataMatrix->hsize;

  free(dataMatrix->iia);

#ifdef CONTIG_MATRIX
  free(dataMatrix->jjcontig);
  free(dataMatrix->jja);

  free(dataMatrix->valcontig);
  free(dataMatrix->val);
#else
  for (int i = 0; i < hsize; i++)
  {
    //free(spmatrix->jja[i]);
  }
  free(dataMatrix->jja);

  for (int i = 0; i < hsize; i++)
  {
    free(dataMatrix->val[i]);
  }
  free(dataMatrix->val);
#endif

  dataMatrix->hsize = 0;
  dataMatrix->msize = 0;
  dataMatrix->bandwidth = 0;

  dataMatrix->minEval = ZERO;
  dataMatrix->maxEval = ZERO;
  dataMatrix->maxMinusMin = ZERO;
}

/// \details
/// Calculate sparcity statistics for a sparse matrix
void sparsity(struct DataMatrixSt* dataMatrix)
{
  int hsize = dataMatrix->hsize;
  int hValCount=0;
  int hDist[hsize];

  memset(hDist, 0, hsize*sizeof(int));

  for (int i = 0; i < hsize; i++)
  {
    hValCount += dataMatrix->iia[i];
    if (dataMatrix->iia[i] > 0)
      hDist[dataMatrix->iia[i]] += 1;
  }

  if (printRank())
  {
    printf("\nSparsity:\nInitial sparsity = %d, fraction = %e, Avg per row = %f\n", hValCount,
       (real_t)hValCount/(real_t)(hsize*hsize), (real_t)hValCount/(real_t)hsize);
    int maxRowCount = 0;
    for (int i = 0; i < hsize; i++)
    {
       maxRowCount = MAX(maxRowCount, dataMatrix->iia[i]);
    }
    printf("Max per row = %d\n", maxRowCount);
    for (int i = 0; i < hsize; i++)
    {
      if (hDist[i] > 0)
         printf("I = %d, count = %d, fraction = %f\n", 
           i, hDist[i], (real_t)hDist[i]/(real_t)hsize);
    }
  }
}

/// \details
/// Calculate gershgorin bounds for sparse matrix
void gershgorin(struct DataMatrixSt* dataMatrix, struct DomainSt* domain)
{
  int hsize = dataMatrix->hsize;
  real_t eMin = 10000;
  real_t eMax = -10000;

  real_t sumP, sumM, maxMinusMin;

#pragma omp parallel for private(sumM,sumP) reduction(max:eMax) reduction(min:eMin)
  for(int i = 0; i < hsize; i++)
  {
    sumM = 0.0;

    for(int j = 0; j < dataMatrix->iia[i]; j++) {
      real_t hx = ABS(dataMatrix->val[i][j]);
      sumM += hx;
      if (dataMatrix->jja[i][j] == i)
      {
        sumP = dataMatrix->val[i][j];
        sumM -= hx;
      }
    }
    eMax = ((eMax < (sumP + sumM)) ? sumP + sumM : eMax);
    eMin = ((eMin > (sumP - sumM)) ? sumP - sumM : eMin);

  }

  // Determine eMax and eMin across ranks
#ifdef DO_MPI
  if (getNRanks() > 1)
  {
    startTimer(reduceCommTimer);
    minRealReduce(&eMin);
    stopTimer(reduceCommTimer);
    collectCounter(reduceCounter, sizeof(real_t));
   
    startTimer(reduceCommTimer);
    maxRealReduce(&eMax);
    stopTimer(reduceCommTimer);
    collectCounter(reduceCounter, sizeof(real_t));
  }
#endif
    
  maxMinusMin = eMax-eMin;

  if (printRank()) 
    printf("\nGershgorin:\nNew  eMax, eMin = %e, %e\n", eMax, eMin); // GERSGORIN BOUNDS;

  dataMatrix->maxEval = eMax;
  dataMatrix->minEval = eMin;
  dataMatrix->maxMinusMin = maxMinusMin;
}

/// \details
/// Normalize a matrix in sparse format using the gershgorin estimates
void normalize(struct DataMatrixSt* dataMatrix)
{
  int hsize = dataMatrix->hsize;
  int sumIia = 0;
  int maxIia = 0;

#pragma omp parallel for reduction(+:sumIia) reduction(max:maxIia)
  for(int i = 0; i < hsize; i++)
  {
    for(int j = 0; j < dataMatrix->iia[i]; j++)
    {
      if (dataMatrix->jja[i][j] == i) {
        dataMatrix->val[i][j] = (dataMatrix->maxEval - dataMatrix->val[i][j])/dataMatrix->maxMinusMin;
      } else {
        dataMatrix->val[i][j] = -dataMatrix->val[i][j]/dataMatrix->maxMinusMin;
      }
    }
    sumIia += dataMatrix->iia[i];
    maxIia = MAX(maxIia, dataMatrix->iia[i]);
 }

 // WE NOW HAVE X = (eMax*I-H)/(eMax-eMin)
 if (printRank() && debug == 1)
   printf("Initial sparsity normalized = %d, fraction = %e,  avg = %g, max = %d\n",
       sumIia, (real_t)sumIia/(real_t)(hsize*hsize), (real_t)sumIia/(real_t)hsize, maxIia);
}

/// \details
/// Calculate trace and trace^2 for a sparse matrix.
void trace(struct DataMatrixSt* dataMatrix, struct DomainSt* domain, real_t* tr, real_t* tr2)
{
  int hsize = dataMatrix->hsize;
  real_t trace = ZERO;
  real_t trace2 = ZERO;

#pragma omp parallel for reduction(+:trace, trace2)
  for(int i = domain->localRowMin; i < domain->localRowMax; i++)
  {
#ifdef POS1
    // Diagonal values are in first position
    trace += dataMatrix->val[i][0];
    trace2 += dataMatrix->val[i][0] * dataMatrix->val[i][0];
#else
    for(int j = 0; j < dataMatrix->iia[i]; j++)
    {
      if (i == dataMatrix->jja[i][j])
      {
        trace += dataMatrix->val[i][j];
        trace2 += dataMatrix->val[i][j] * dataMatrix->val[i][j];
      }
    }
#endif
  }
  
  *tr = trace;
  *tr2 = trace2;
}

#endif
