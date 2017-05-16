/// \file
/// Generate Hamiltonian matrix.

#ifndef __GENERATE_H
#define __GENERATE_H

#include "matrixMath.h"
#include "mytype.h"

/// Generate Hamiltonian given N and M
DataMatrix* generateHMatrix(const MatrixType mtype, int N, int M, real_t a, real_t alpha);

#endif
