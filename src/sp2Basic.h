/// \file
/// SP2 loop functions.

#ifndef __SP2BASIC_H
#define __SP2BASIC_H

#include "bml.h"

#include <stdio.h>

#include "mytype.h"

void normalize(bml_matrix_t* h_bml);
void sp2Loop(bml_matrix_t* h_bml, bml_matrix_t* rho_bml, real_t threshold, real_t bndfil, int minsp2iter, int maxsp2iter, real_t idemTol);
void reportResults(int iter, bml_matrix_t* rho_bml, bml_matrix_t* x2_bml);

#endif
