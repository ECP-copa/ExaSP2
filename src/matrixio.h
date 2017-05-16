/// \file
/// Matrix I/O.

#ifndef __MATRIXIO_H
#define __MATRIXIO_H

#include "sparseMatrix.h"
#include "mytype.h"

/// Read/write hamiltonian and density matrices.
void writeDataPattern(char* fname, struct DataMatrixSt* dataMatrix, real_t hthresh);
void readMTX(char* fname, struct DataMatrixSt* dataMatrix);
void writeMTX(char* fname, struct DataMatrixSt* dataMatrix);

#endif
