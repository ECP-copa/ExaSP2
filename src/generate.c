/// \File
/// Routines for generating synthetic Hamiltonian matrices.

#ifdef MMATH_SPARSE

#include "generate.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrixMath.h"
#include "random.h"
#include "constants.h"

/// \details
/// Generate Hamiltonian sparse matrix
DataMatrix* generateHMatrix(const MatrixType mtype, int N, int M, real_t a, real_t alpha)
{
  DataMatrix* dataMatrix = initDataMatrix(mtype, N, M);

  //
  // hx = amplitude * random[0,1] * exp(-alpha * (i-j)^2)
  //
  // This will produce a hamiltonian similar to polyethylene.
  //
  // A more general version based on atom positions:
  //
  // hx = amplitude * random[0,1] * exp(-alpha * (Ri-Rj)^2)
  //      where (Ri - Rj) is the distance between atoms
  //

  uint64_t seed = mkSeed(N, M);

  int l, mstart, mend;
  real_t hx;

  int nnz = 0;
  for (int i = 0; i < N; i++)
  {
    l = 0;
    mstart = MAX(0, i-M+1);
    mend = MIN(N, i+M);
    for (int j = mstart; j < mend; j++)
    {
      hx = a * lcg61(&seed) * exp(-alpha * (i - j) * (i - j));
      if (j == i || ABS(hx) > eps)
      {
        dataMatrix->jja[i][l] = j;
        dataMatrix->val[i][l] = hx;
        l++;
        nnz++;
      }
    }
    dataMatrix->iia[i] = l;
  }

  printf("Generated H Matrix nnz = %d avg nnz/row = %d\n", nnz, nnz/N);

  return dataMatrix;
}

#endif
