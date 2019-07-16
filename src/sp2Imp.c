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

  bml_matrix_t* I_bml = bml_identity_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* xtmp_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* p2_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* ai_bml = bml_identity_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* x_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* a_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);  


  const real_t cnst = beta/pow(2,2+rec_steps);
  const int ns_iter = 4;
  const real_t cg_tol = 0.1;
  int i,j;
  real_t norm;

  // Do normalization
  startTimer(normTimer);
  bml_copy(h_bml, p_bml);
  normalize(p_bml,cnst,mu);
  stopTimer(normTimer);
  
  for (i=1; i <= rec_steps; i++) {
   
    printf("i = %d\n",i );
    //bml_print_bml_matrix(p_bml, 0, 10, 0, 10);
    startTimer(x2Timer);
    bml_multiply_x2(p_bml, p2_bml, threshold);
    stopTimer(x2Timer);
    bml_copy(p2_bml, a_bml);
    bml_add(a_bml, p_bml, ONE, MINUS_ONE, threshold);
    bml_add(a_bml, I_bml, TWO, ONE, threshold);
   
    if (method == 0) {
	    if (i == 1) {
		printf("Starting cg\n");
		startTimer(inverseTimer);
		//ai_bml = bml_inverse(a_bml);
		bml_print_bml_matrix(ai_bml, 0, 10, 0, 10);
		conjugateGradient(a_bml, I_bml, ai_bml, cg_tol, threshold);
		bml_print_bml_matrix(ai_bml, 0, 10, 0, 10);
		stopTimer(inverseTimer);
		printf("Done with cg\n");
	    }
	    else {
		startTimer(nsiterTimer);
		for (j=1; j <= ns_iter; j++) {
		
			bml_copy(ai_bml, xtmp_bml);
			bml_multiply(ai_bml, a_bml, x_bml, MINUS_ONE, ZERO, threshold);
			bml_multiply(x_bml, xtmp_bml, ai_bml, ONE, TWO, threshold);
			startTimer(normTimer);
			bml_add(xtmp_bml, ai_bml, ONE, MINUS_ONE, threshold);
			norm = bml_fnorm(xtmp_bml);
			printf("j = %d, norm = %lf\n", j, norm);
			stopTimer(normTimer); 
		}
		stopTimer(nsiterTimer);
	    }
	    startTimer(mmTimer);
	    bml_multiply(ai_bml, p2_bml, p_bml, ONE, ZERO, threshold);
	    stopTimer(mmTimer);
     }
     else {
	    conjugateGradient(a_bml, p2_bml, p_bml, cg_tol, threshold);
     }
  }
  
  bml_print_bml_matrix(p_bml, 0, 10, 0, 10);
 
  bml_deallocate(&ai_bml);
  bml_deallocate(&a_bml);
  bml_deallocate(&I_bml);
  bml_deallocate(&p2_bml);
  bml_deallocate(&xtmp_bml);
  bml_deallocate(&x_bml);
  stopTimer(sp2LoopTimer);	 
}
 
  // Report results
 /* reportResults(iter, rho_bml, x2_bml);*/

void conjugateGradient(const bml_matrix_t* A_bml, 
	               const bml_matrix_t* b_bml, 
	               bml_matrix_t* x_bml, 
	               const real_t cg_tol,
	               const real_t threshold) {

  bml_matrix_type_t bml_type = bml_get_type(x_bml);
  int N = bml_get_N(x_bml);
  int M = bml_get_M(x_bml);
  bml_matrix_precision_t precision = bml_get_precision(x_bml);
  bml_distribution_mode_t dmode = bml_get_distribution_mode(x_bml);

  bml_matrix_t* r_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* p_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* wtmp_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* w_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  real_t alpha, beta, r_norm_old, r_norm_new;
  int k = 0;

  bml_multiply_AB(A_bml, x_bml, r_bml, threshold);
  bml_add(r_bml, b_bml, MINUS_ONE, ONE, threshold);
  r_norm_old = bml_sum_squares(r_bml);
  r_norm_new = r_norm_old;
  
  while (cg_tol < r_norm_new) {
    printf("%lf\n", r_norm_new);
    k++;
    if (k == 1) {
      bml_copy(r_bml, p_bml);
    }
    else {
      beta = r_norm_new/r_norm_old;
      bml_add(p_bml, r_bml, ONE, beta, threshold);
    }
    bml_multiply_AB(A_bml, p_bml, wtmp_bml, threshold);
    bml_multiply_AB(p_bml, wtmp_bml, w_bml, threshold);
    alpha = r_norm_new/bml_trace(w_bml);
    bml_add(x_bml, p_bml, ONE, alpha, threshold);
    bml_add(r_bml, wtmp_bml, ONE, -alpha, threshold);
    r_norm_old = r_norm_new;
    r_norm_new = bml_sum_squares(r_bml);
  }
  printf("%lf\n", r_norm_new);
  bml_deallocate(&r_bml);
  bml_deallocate(&p_bml);
  bml_deallocate(&w_bml);
  bml_deallocate(&wtmp_bml);
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
