/// \file
/// SP2 loop.

#ifdef SP2_IMP

#include "bml.h"

#include "sp2Imp.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "performance.h"
#include "parallel.h"
#include "constants.h"

/// \details
/// Normalize a Hamiltonian matrix prior to running the implicit recursive algorithm.
/// 
/// X0 = cnst(mu*I - F) + 0.5I
/// 
void normalize(bml_matrix_t* h_bml, const real_t cnst, const real_t mu)
{
  real_t alpha, beta; 
  alpha = -cnst;
  beta = cnst*mu + HALF;
  bml_scale_add_identity(h_bml, alpha, beta, ZERO);
}

/// \details
/// The implicit recursive expansion algorithm.
void implicit_recursiveLoops(const bml_matrix_t* h_bml, 
             bml_matrix_t* p_bml, 
             const int rec_steps, 
             const real_t linsysTol, 
             const real_t threshold)
{

  startTimer(sp2LoopTimer);

  int N = bml_get_N(h_bml);
  int M = bml_get_M(h_bml);
  bml_matrix_type_t bml_type = bml_get_type(h_bml);
  bml_matrix_precision_t precision = bml_get_precision(h_bml);
  bml_distribution_mode_t dmode = bml_get_distribution_mode(h_bml);    

  bml_matrix_t* I_bml = bml_identity_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* xtmp_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* p2_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* ai_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* x_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* a_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);

  const real_t cnst = beta/pow(2,2+rec_steps);
  const int ns_iter = 4;
  int i,j;

  // Do normalization
  startTimer(normTimer);
  bml_copy(h_bml, p_bml);
  normalize(p_bml,cnst,mu);
  stopTimer(normTimer);
  
  for (i=1; i <= rec_steps; i++) {

    bml_multiply_x2(p_bml, p2_bml, threshold);
    bml_copy(p2_bml, a_bml);
    bml_add(a_bml, p_bml, ONE, MINUS_ONE, threshold);
    bml_add(a_bml, I_bml, TWO, ONE);

    if (i == 1) {
    	ai_bml = bml_inverse(a_bml);
    }
    else {
	for (j=1; j <= ns_iter; j++) {
		bml_copy(ai_bml, xtmp_bml);
		bml_multiply(ai_bml, a_bml, x_bml, MINUS_ONE, ZERO, threshold);
		bml_multiply(x_bml, xtmp_bml, ai_bml, ONE, TWO, threshold);
	}
    }
    bml_multiply(ai_bml, p2_bml, p_bml, ONE, ZERO, threshold);
  }  
 
  // Report results
 /* reportResults(iter, rho_bml, x2_bml);*/


void reportResults(const int iter, 
                   const bml_matrix_t* p_bml, 
                   const bml_matrix_t* x2_bml)
{
    int sumIIA = 0;
    int sumIIC = 0;
//  int sumIIA = bml_get_sparsity(rho_bml);
//  int sumIIC = bml_get_sparsity(x2_bml);
  int maxIIA = bml_get_bandwidth(p_bml);
  int maxIIC = bml_get_bandwidth(x2_bml);


  if (bml_printRank())
  {
    printf("\nResults:\n");
    printf("X2 Sparsity CCN = %d, fraction = %e avg = %g, max = %d\n", sumIIC, 
      (real_t)sumIIC/(real_t)(N_i*N_i), (real_t)sumIIC/(real_t)N_i, maxIIC);

    printf("RHO Sparsity AAN = %d, fraction = %e avg = %g, max = %d\n", sumIIA, 
      (real_t)sumIIA/(real_t)(N_i*N_i), (real_t)sumIIA/(real_t)N_i, maxIIA);

    printf("Number of iterations = %d\n", iter);
  }
}

#endif
