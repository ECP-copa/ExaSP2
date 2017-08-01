/// \file
/// SP2 loop.

#ifdef SP2_FERMI

#include "bml.h"

#include "sp2Fermi.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "performance.h"
#include "parallel.h"
#include "constants.h"

/// \details
/// Normalize a Hamiltonian matrix prior to running the SP2 Fermi algorithm.
/// 
/// X0 = ((hN-mu) * I - H) / (hN - h1)
///  or X0 = (hN*I-H0-mu*I)/(hN-h1)
/// 
/// where h1 and hN are scaled Gershgorin bounds.
///
void normalize(bml_matrix_t* h_bml, 
               const real_t h1, 
               const real_t hN, 
               const real_t mu)
{
  real_t alpha, beta, maxMinusMinInverse;

  maxMinusMinInverse = ONE / (hN - h1);
  alpha = MINUS_ONE;
  beta = hN - mu;

  bml_scale_add_identity(h_bml, alpha, beta, ZERO);
  bml_scale_inplace(maxMinusMinInverse, h_bml);

}

/// \details
/// Finite temperature truncated SP2 Fermi init.
/// The second order spectral projection algorithm.
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
             const real_t threshold)
{

  int N = bml_get_N(h_bml);
  int M = bml_get_M(h_bml);
  bml_matrix_type_t bml_type = bml_get_type(h_bml);
  bml_matrix_precision_t precision = bml_get_precision(h_bml);
  bml_distribution_mode_t dmode = bml_get_distribution_mode(h_bml);

  // Calculate gershgorin bounds and rescale
  real_t* gbnd = bml_gershgorin(h_bml);
  *mu = 0.5 * (gbnd[1] + gbnd[0]);
  *h1 = tscale * gbnd[0];
  *hN = tscale * gbnd[1];
  bml_free_memory(gbnd);

  real_t traceX, traceX0, traceX1, traceX2;
  real_t lambda = ZERO;
  real_t occErr = ONE;
  int firstTime = 1;

  real_t* trace = bml_allocate_memory(2*sizeof(real_t));

  bml_matrix_t* i_bml = bml_identity_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* x1_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* x2_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);
  bml_matrix_t* tmp_bml = bml_zero_matrix(bml_type, precision, N, M, dmode);

  while (occErr > occErrLimit)
  {
    startTimer(normTimer);
    bml_copy(h_bml, rho_bml);
    normalize(rho_bml, *h1, *hN, *mu);
    stopTimer(normTimer);

    bml_copy(i_bml, x1_bml);
    real_t sfactor = MINUS_ONE / (hN - h1);
    bml_scale_inplace(sfactor, x1_bml);

    for (int i = 0; i < nsteps; i++)
    {
      startTimer(x2InitTimer);
      trace = bml_multiply_x2(rho_bml, x2_bml, threshold);
      stopTimer(x2InitTimer);
      traceX0 = trace[0];
      traceX2 = trace[1];

      // First time through determine sequence branching
      if (firstTime == 1)
      {
        if (ABS(traceX2-nocc) < ABS(TWO * traceX0 - traceX2 - nocc))
          sgnlist[i] = -1;
        else
          sgnlist[i] = 1;
      }

      // X1 = X1 + sgnlist(i)*(X1 - X0*X1 - X1*X0)
      // if sgnlist == 1, X1 = 2 * X1 - (X0*X1 + X1*X0)
      // if sgnlist == -1, X1 = X0*X1 + X1*X0
      // 
      // tmp = X0*X1 + X1*X0
      startTimer(mmInitTimer);
      bml_multiply(rho_bml, x1_bml, tmp_bml, ONE, ZERO, threshold);
      stopTimer(mmInitTimer);
      startTimer(mmInitTimer);
      bml_multiply(x1_bml, rho_bml, tmp_bml, ONE, ONE, threshold);
      stopTimer(mmInitTimer);

      if (sgnlist[i] == 1)
      {
        // X1 = 2 * X1 - tmp
        startTimer(xaddInitTimer);
        bml_add(x1_bml, tmp_bml, TWO, MINUS_ONE, threshold);
        stopTimer(xaddInitTimer);
      }
      else
        bml_copy(tmp_bml, x1_bml);

      // X0 = X0 + sgnlist(i)*(X0 - X0_2)
      // if sgnlist == 1, X0 = 2.0*X0 - X0_2
      // if sgnlist == -1, X0 = X0_2
      // 
      if (sgnlist[i] == 1)
      {
        // X0 = 2 * X0 - X2
        startTimer(xaddInitTimer);
        bml_add(rho_bml, x2_bml, TWO, MINUS_ONE, threshold);
        stopTimer(xaddInitTimer);
      }
      else
        bml_copy(x2_bml, rho_bml);

    }

    firstTime = 0;
    traceX0 = bml_trace(rho_bml);
    traceX1 = bml_trace(x1_bml);
    occErr = ABS(nocc - traceX0);

    // Newton=Rhapson step to correct for occupation
    if (ABS(traceX1) > traceLimit)
      lambda = (nocc - traceX0) / traceX1;
    else
      lambda = ZERO;

    *mu += lambda;
    printf("mu = %lg occErr = %lg occErrLimit = %lg\n", *mu, occErr, occErrLimit);
  }

  bml_free_memory(trace);
  bml_deallocate(&x2_bml);

  // X0*(I-X0)
  // I = I - X0
  startTimer(xaddInitTimer);
  bml_add(i_bml, rho_bml, ONE, MINUS_ONE, threshold);
  stopTimer(xaddInitTimer);

  // tmp = X0*I
  startTimer(mmInitTimer);
  bml_multiply(rho_bml, i_bml, tmp_bml, ONE, ZERO, threshold);
  stopTimer(mmInitTimer);

  traceX = bml_trace(tmp_bml);
  traceX1 = bml_trace(x1_bml);

  if (ABS(traceX) > traceLimit)
    *beta = -traceX1 / traceX;
  else
    *beta = MINUS_THOUSAND;

  // X = 2 * X
  bml_scale_inplace(TWO, rho_bml);

  bml_deallocate(&tmp_bml);
  bml_deallocate(&i_bml);
  bml_deallocate(&x1_bml);

} 

/// \details
/// Finite temperature truncated SP2 Fermi loop.
/// The second order spectral projection algorithm.
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
             const real_t occLimit,
             const real_t traceLimit,
             const real_t threshold)
{
  bml_matrix_type_t matrix_type = bml_get_type(h_bml);
  int N = bml_get_N(h_bml);
  int M = bml_get_M(h_bml);
  bml_matrix_precision_t precision = bml_get_precision(h_bml);
  bml_distribution_mode_t dmode = bml_get_distribution_mode(h_bml);

  real_t* trace = bml_allocate_memory(2*sizeof(real_t));

  bml_matrix_t* i_bml = bml_identity_matrix(matrix_type, precision, N, M, dmode);
  bml_matrix_t* dx_bml = bml_zero_matrix(matrix_type, precision, N, M, dmode);
  bml_matrix_t* x2_bml = bml_zero_matrix(matrix_type, precision, N, M, dmode);

  real_t traceX0, traceX2, traceDX, lambda;
  real_t occErr = ONE + occLimit;
  int iter = 0;

  while ((osteps == 0 && occErr > occLimit) ||
         (osteps > 0 && iter < osteps))
  {
    iter += 1;
    startTimer(normTimer);
    bml_copy(h_bml, rho_bml);
    normalize(rho_bml, h1, hN, *mu);
    stopTimer(normTimer);

    for (int i = 0; i < nsteps; i++)
    {
      startTimer(x2Timer);
      trace = bml_multiply_x2(rho_bml, x2_bml, threshold);
      stopTimer(x2Timer);
      traceX0 = trace[0];
      traceX2 = trace[1];

      // X0 = X0 + sgnlist(i)*(X0 - X0_2)
      if (sgnlist[i] == 1)
      {
        startTimer(xaddTimer);
        bml_add(rho_bml, x2_bml, TWO, MINUS_ONE, threshold);
        stopTimer(xaddTimer);
      }
      else
        bml_copy(x2_bml, rho_bml);
    }

    traceX0 = bml_trace(rho_bml);
    occErr = ABS(nocc - traceX0);

    // DX = -beta*X0*(I-X0)
    bml_copy(i_bml, x2_bml);
    startTimer(xaddTimer);
    bml_add(x2_bml, rho_bml, ONE, MINUS_ONE, threshold);
    stopTimer(xaddTimer);
    startTimer(mmTimer);
    bml_multiply(rho_bml, x2_bml, dx_bml, -beta, ZERO, threshold);
    stopTimer(mmTimer);
    traceDX = bml_trace(dx_bml);

    // Newton-Rhapson step to correct for occupation
    if (ABS(traceDX) > traceLimit)
      lambda = (nocc - traceX0) / traceDX;
    else
      lambda = ZERO;

    *mu += lambda;

  }

  // Correction for occupation
  startTimer(xaddTimer);
  bml_add(rho_bml, dx_bml, ONE, lambda, threshold);
  stopTimer(xaddTimer);

  // X = 2*X
  bml_scale_inplace(TWO, rho_bml);

  bml_free_memory(trace);
  bml_deallocate(&i_bml);
  bml_deallocate(&x2_bml);
  bml_deallocate(&dx_bml);
}

/// \details
/// Report density matrix results
void reportResults(const int iter, 
                   const bml_matrix_t* rho_bml, 
                   const bml_matrix_t* x2_bml)
{
    int sumIIA = 0;
    int sumIIC = 0;
//  int sumIIA = bml_get_sparsity(rho_bml);
//  int sumIIC = bml_get_sparsity(x2_bml);
  int maxIIA = bml_get_bandwidth(rho_bml);
  int maxIIC = bml_get_bandwidth(x2_bml);

#ifdef DO_MPI
  // Collect number of non-zeroes and max non-zeroes per row across ranks
/*
  if (bml_getNRanks() > 1)
  {
    startTimer(reduceCommTimer);
    addIntReduce2(&sumIIA, &sumIIC);
    stopTimer(reduceCommTimer);
    collectCounter(reduceCounter, 2 * sizeof(int));

    startTimer(reduceCommTimer);
    maxIntReduce2(&maxIIA, &maxIIC);
    stopTimer(reduceCommTimer);
    collectCounter(reduceCounter, 2 * sizeof(int));
  }
*/
#endif

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
