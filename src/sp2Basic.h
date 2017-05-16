/// \file
/// SP2 loop functions.

#ifndef __SP2BASIC_H
#define __SP2BASIC_H

#include <stdio.h>

#include "decomposition.h"
#include "matrixMath.h"
#include "mytype.h"

void sp2Loop(struct DataMatrixSt* xmatrix, struct DomainSt* domain);
void reportResults(int iter, struct DataMatrixSt* xmatrix, struct DataMatrixSt* x2matrix, struct DomainSt* domain);

#endif
