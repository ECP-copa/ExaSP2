/// \file
/// Sparse matrix functions.

#ifdef MMATH_SPARSE

#include "sparseMath.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "sparseMatrix.h"
#include "parallel.h"
#include "constants.h"

// \details
// Sparse matrix multiply X^2
void multiplyX2(real_t* trX, real_t* trX2, struct DataMatrixSt* xmatrix, struct DataMatrixSt* x2matrix, struct DomainSt* domain)
{
  int hsize = xmatrix->hsize;
  int ix[hsize];
  real_t x[hsize];
  real_t traceX = ZERO;
  real_t traceX2 = ZERO;

  memset(ix, 0, hsize*sizeof(int));
  memset(x, ZERO, hsize*sizeof(real_t));

#pragma omp parallel for firstprivate(ix,x) reduction(+:traceX,traceX2)
  for(int i = domain->localRowMin; i < domain->localRowMax; i++)  
  // CALCULATES THRESHOLDED X^2
  {
    int l = 0;
    for(int jp = 0; jp < xmatrix->iia[i]; jp++)
    {
      real_t a = xmatrix->val[i][jp];
      int j = xmatrix->jja[i][jp];
      if (j == i)
      {
        traceX += a;
      }
      for(int kp = 0; kp < xmatrix->iia[j]; kp++)
      {
        int k = xmatrix->jja[j][kp];
        if (ix[k] == 0)
        {
          x[k] = ZERO;
          x2matrix->jja[i][l] = k;
          ix[k] = i+1;
          l++;
        }
        x[k] = x[k] + a * xmatrix->val[j][kp]; // TEMPORARY STORAGE VECTOR LENGTH FULL N
      }
    }

    int ll = 0;
    for(int j = 0; j < l; j++)
    {
      int jp = x2matrix->jja[i][j];
      real_t xtmp = x[jp];
      if (jp == i)
      {
        traceX2 += xtmp;
        x2matrix->val[i][ll] = xtmp;
        x2matrix->jja[i][ll] = jp;
        ll++;
      }
      else if(ABS(xtmp) > eps)
      {
        x2matrix->val[i][ll] = xtmp;
        x2matrix->jja[i][ll] = jp;
        ll++;
      }
      ix[jp] = 0;
      x[jp] = ZERO;
    }
    x2matrix->iia[i] = ll;
  }

  *trX = traceX;
  *trX2 = traceX2;

}

// \details
// Sparse matrix add - 2*x - x^2
void add(struct DataMatrixSt* xmatrix, struct DataMatrixSt* x2matrix, struct DomainSt* domain)
{
  int hsize = xmatrix->hsize;
  int ix[hsize];
  real_t x[hsize];

  memset(ix, 0, hsize*sizeof(int));
  memset(x, ZERO, hsize*sizeof(real_t));

#pragma omp parallel for firstprivate(x,ix)
  for (int i = domain->localRowMin; i < domain->localRowMax; i++) // X = 2X-X^2
  {
    int l = 0;
    for(int jp = 0; jp < xmatrix->iia[i]; jp++)
    {
      int k = xmatrix->jja[i][jp];
      if (ix[k] == 0)
      {
        x[k] = ZERO;
        ix[k] = i+1;
        xmatrix->jja[i][l] = k;
        l++;
      }
      x[k] = x[k] + TWO*xmatrix->val[i][jp];
    }

    for(int jp = 0; jp < x2matrix->iia[i]; jp++)
    {
      int k = x2matrix->jja[i][jp];
      if (ix[k] == 0)
      {
        x[k] = ZERO;
        ix[k] = i+1;
        xmatrix->jja[i][l] = k;
        l++;
      }
      x[k] = x[k] - x2matrix->val[i][jp];
    }
    xmatrix->iia[i] = l;

    int ll = 0;
    for(int jp = 0; jp < l; jp++)
    {
      real_t xTmp = x[xmatrix->jja[i][jp]];
      if (ABS(xTmp) > eps) // THIS THRESHOLDING COULD BE IGNORED!?
      {
        xmatrix->val[i][ll] = xTmp;
        xmatrix->jja[i][ll] = xmatrix->jja[i][jp];
        ll++;
      }
      x[xmatrix->jja[i][jp]] = ZERO;
      ix[xmatrix->jja[i][jp]] = 0;
    }
    xmatrix->iia[i] = ll;
  }

}

// \details
// Sparse matrix set X^2 
void setX2(struct DataMatrixSt* xmatrix, struct DataMatrixSt* x2matrix, struct DomainSt* domain)
{
  int hsize = xmatrix->hsize;
  
#pragma omp parallel for
  for(int i = domain->localRowMin; i < domain->localRowMax; i++) // X = X^2
  {
    real_t xtmp;
    int ll = 0;
    for(int jp = 0; jp < x2matrix->iia[i]; jp++)
    {
      xtmp = x2matrix->val[i][jp];
      if (ABS(xtmp) > eps)
      {
        xmatrix->val[i][ll] = xtmp;
        xmatrix->jja[i][ll] = x2matrix->jja[i][jp];
        ll++;
      }
    }
    xmatrix->iia[i] = ll;
  }
}

// \details
// Multiple a Sparse matrix by a scalar
void multiplyScalar(struct DataMatrixSt* xmatrix, struct DomainSt* domain, real_t scalar)
{
#pragma omp parallel for
  for (int i = domain->localRowMin; i < domain->localRowMax; i++) 
  {
    for (int j = 0; j < xmatrix->iia[i]; j++)
    {
      xmatrix->val[i][j] *= scalar;
    }
  }
}

#endif
