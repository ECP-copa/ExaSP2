/// \file
/// SP2 Fermi functions.

#ifndef __SP2FERMI_H
#define __SP2FERMI_H

#include "bml.h"

#include <stdio.h>

#include "mytype.h"

void normalize(bml_matrix_t* h_bml, 
               const real_t h1, 
               const real_t hN, 
               const real_t mu);

void sp2Init(const bml_matrix_t* h_bml, 
             bml_matrix_t* rho_bml, 
             const int nsteps, 
             const real_t nocc, 
             real_t* mu, 
             real_t* beta, 
             int* sgnlist, 
             real_t* h1,
             real_t* hN,
             const real_t tscale,
             const real_t occErrLimit, 
             const real_t traceLimit, 
             const real_t threshold);

void sp2Loop(const bml_matrix_t* h_bml,
             bml_matrix_t* rho_bml,
             const int nsteps, 
             const real_t nocc,
             real_t* mu, 
             const real_t beta, 
             const int* sgnlist,
             const real_t h1, 
             const real_t hN,
             const int osteps,
             const real_t eps, 
             const real_t traceLimit,
             const real_t threshold);

void reportResults(const int iter, 
                   const bml_matrix_t* rho_bml, 
                   const bml_matrix_t* x2_bml);

#endif
