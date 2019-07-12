/// \file
/// SP2 loop functions.

#ifndef __SP2BASIC_H
#define __SP2BASIC_H

#include "bml.h"

#include <stdio.h>

#include "mytype.h"

void normalize(bml_matrix_t* h_bml);

void sp2Loop(const bml_matrix_t* h_bml, 
             bml_matrix_t* rho_bml, 
             const real_t nocc, 
             const int minsp2iter, 
             const int maxsp2iter, 
             const real_t idemTol,
             const real_t threshold);

void reportResults(const int iter,
                   const bml_matrix_t* rho_bml, 
                   const bml_matrix_t* x2_bml);

#endif
