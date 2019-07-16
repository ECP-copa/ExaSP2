/// \file
/// SP2 loop functions.

#ifndef __SP2IMP_H
#define __SP2IMP_H

#include "bml.h"

#include <stdio.h>

#include "mytype.h"

void normalize(bml_matrix_t* h_bml, const real_t cnst, const real_t mu);

void implicit_recursiveLoops(const bml_matrix_t* h_bml, 
             bml_matrix_t* p_bml,
	     const real_t beta,
             const real_t mu, 
             const int rec_steps, 
             const real_t threshold);

void conjugateGradient(const bml_matrix_t* A_bml,
                       const bml_matrix_t* b_bml,
                       bml_matrix_t* x_bml,
                       const real_t cg_tol,
                       const real_t threshold); 

void reportResults(const int iter,
                   const bml_matrix_t* rho_bml, 
                   const bml_matrix_t* x2_bml);

#endif
