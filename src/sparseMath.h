/// \file
/// Sparse matrix math routines.

#ifndef __SPARSEMATH_H
#define __SPARSEMATH_H

#include "sparseMatrix.h"

#include <stdio.h>

#include "decomposition.h"
#include "mytype.h"

void multiplyX2(real_t* trX, real_t* trX2, struct DataMatrixSt* xmatrix, struct DataMatrixSt* x2matrix, struct DomainSt* domain);
void add(struct DataMatrixSt* xmatrix, struct DataMatrixSt* x2matrix, struct DomainSt* domain);
void setX2(struct DataMatrixSt* xmatrix, struct DataMatrixSt* x2matrix, struct DomainSt* domain);
void multiplyScalar(struct DataMatrixSt* xmatrix, struct DomainSt* domain, real_t scalar);

#endif
