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
	     const real_t beta,
	     const real_t mu, 
             const int rec_steps,    
             const real_t threshold)
{
  startTimer(sp2LoopTimer);

  const int method = 0;

  int N = bml_get_N(h_bml);
  int M = bml_get_M(h_bml);
  bml_matrix_type_t bml_type = bml_get_type(h_bml);
  bml_matrix_precision_t precision = bml_get_precision(h_bml);
  bml_distribution_mode_t dmode = bml_get_distribution_mode(h_bml);    

  bml_matrix_t* xtmp_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* p2_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* ai_bml = bml_identity_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* x_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* a_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* w_bml;
  bml_matrix_t* I_bml;
  if (method == 1) {  
     I_bml = bml_identity_matrix(bml_type, precision, N, M, dmode);
     w_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  }

  const real_t cnst = beta/pow(2,2+rec_steps);
  const real_t cg_tol = 0.000001;
  int i,j;
  real_t norm;

  // Normalize hamiltonian 
  startTimer(normTimer);
  bml_copy(h_bml, p_bml);
  normalize(p_bml, cnst, mu);
  stopTimer(normTimer);
  
  for (i=1; i <= rec_steps; i++) {
    
    // Set up linear system 
    startTimer(linsyssetupTimer);
    bml_multiply_x2(p_bml, p2_bml, threshold);
    bml_copy(p2_bml, a_bml);
    bml_add(a_bml, p_bml, ONE, MINUS_ONE, threshold);
    bml_scale_add_identity(a_bml, TWO, ONE, threshold);
  //  bml_add(a_bml, I_bml, TWO, ONE, threshold);
    stopTimer(linsyssetupTimer);
   
    // Newton-Schulz-method
    if (method == 1) {
	   
		startTimer(nsiterTimer);
		j = 1;
		norm = 1.0;
		// Find inverse with Newton-Schulz iteration
		while (norm > 0.01)   {

			// First time use Conjugate Gradient for inverse starting guess 
			if (i == 1 && j == 1) { 
		     		startTimer(inverseTimer);
				conjugateGradient(a_bml, I_bml, ai_bml, x_bml, xtmp_bml, w_bml, cg_tol, threshold);
				stopTimer(inverseTimer);
			}
			bml_copy(ai_bml, xtmp_bml);
			bml_multiply(ai_bml, a_bml, x_bml, MINUS_ONE, ZERO, threshold);
			bml_multiply(x_bml, xtmp_bml, ai_bml, ONE, TWO, threshold);
			startTimer(normTimer);
			bml_add(xtmp_bml, ai_bml, ONE, MINUS_ONE, threshold);
			norm = bml_fnorm(xtmp_bml);
			stopTimer(normTimer);
			j++;

		}
	    // Get next density matrix  
	    bml_multiply(ai_bml, p2_bml, p_bml, ONE, ZERO, threshold);
	    stopTimer(nsiterTimer);
     }
     // Conjugate gradient method
     else {
	    conjugateGradient(a_bml, p2_bml, p_bml, x_bml, xtmp_bml, ai_bml, cg_tol, threshold);
     }
  }
  
  bml_print_bml_matrix(p_bml, 0, 10, 0, 10);
 
  bml_deallocate(&ai_bml);
  bml_deallocate(&a_bml);
  bml_deallocate(&p2_bml);
  bml_deallocate(&xtmp_bml);
  bml_deallocate(&x_bml);
  if (method == 1) { 
     bml_deallocate(&w_bml);
     bml_deallocate(&I_bml);
  }
  stopTimer(sp2LoopTimer);	 
}
 
  // Report results
 /* reportResults(iter, rho_bml, x2_bml);*/

void conjugateGradient(const bml_matrix_t* A_bml, 
	               bml_matrix_t* p2_bml, 
	               bml_matrix_t* p_bml, 
                       bml_matrix_t* d_bml,
                       bml_matrix_t* wtmp_bml,
                       bml_matrix_t* w_bml,
	               const real_t cg_tol,
	               const real_t threshold) {

  bml_matrix_type_t bml_type = bml_get_type(p_bml);
  int N = bml_get_N(p_bml);
  int M = bml_get_M(p_bml);
  bml_matrix_precision_t precision = bml_get_precision(p_bml);
  bml_distribution_mode_t dmode = bml_get_distribution_mode(p_bml);

  real_t alpha, beta, r_norm_old, r_norm_new;
  int k = 0;

  bml_multiply(A_bml, p_bml, p2_bml, MINUS_ONE, ONE, threshold);
  r_norm_old = bml_sum_squares(p2_bml);
  r_norm_new = r_norm_old;
  
  while (cg_tol < r_norm_new) {
    
    k++;
    if (k == 1) {
      bml_copy(p2_bml, d_bml);
    }
    else {
      beta = r_norm_new/r_norm_old;
      bml_add(d_bml, p2_bml, beta, ONE, threshold);
    }
    bml_multiply_AB(A_bml, d_bml, wtmp_bml, threshold);
    bml_multiply_AB(d_bml, wtmp_bml, w_bml, threshold);
    alpha = r_norm_new/bml_trace(w_bml);
    bml_add(p_bml, d_bml, ONE, alpha, threshold);
    bml_add(p2_bml, wtmp_bml, ONE, -alpha, threshold);
    r_norm_old = r_norm_new;
    r_norm_new = bml_sum_squares(p2_bml);
  }
}


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