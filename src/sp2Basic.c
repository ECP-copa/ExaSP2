/// \file
/// SP2 loop.

#ifdef SP2_BASIC

#include "bml.h"

#include "sp2Basic.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "performance.h"
#include "parallel.h"
#include "constants.h"

/// \details
/// Normalize a Hamiltonian matrix prior to running the SP2 algorithm.
/// 
/// X0 = (e_max * I - H) / (e_max - e_min)
/// 
/// where e_max and e_min are obtained using the Gershgorin circle theorem.
///
void normalize(bml_matrix_t* h_bml)
{
  real_t alpha, beta, maxMinusMin;

  real_t* gbnd = bml_gershgorin(h_bml);
  maxMinusMin = gbnd[1] - gbnd[0];
  alpha = MINUS_ONE / maxMinusMin;
  beta = gbnd[1] / maxMinusMin;
  bml_scale_add_identity(h_bml, alpha, beta, ZERO);

  bml_free_memory(gbnd);
}

/// \details
/// The second order spectral projection algorithm.
void sp2Loop(const bml_matrix_t* h_bml, 
             bml_matrix_t* rho_bml, 
             const real_t nocc,
             const int minsp2iter, 
             const int maxsp2iter, 
             const real_t idemTol, 
             const real_t threshold)
{
  //DataExchange* dataExchange;

  startTimer(sp2LoopTimer);

#ifdef DO_MPI
//  if (getNRanks() > 1)
//    dataExchange = initDataExchange(domain);
#endif

  // Do gershgorin normalization
  startTimer(normTimer);
  bml_copy(h_bml, rho_bml);
  normalize(rho_bml);
  stopTimer(normTimer);
  
  // Basic SP2 algorithm
 
  real_t idempErr = ZERO;
  real_t idempErr1 = ZERO;
  real_t idempErr2 = ZERO;

  real_t trX = ZERO;
  real_t trX2 = ZERO;

  real_t tr2XX2, trXOLD, limDiff;

  real_t* trace;

  int iter = 0;
  int breakLoop = 0;

  if (bml_printRank() && debug_i == 1)
    printf("\nSP2Loop:\n");

  // X2 <- X
  bml_matrix_t* x2_bml = bml_copy_new(rho_bml);

  while ( breakLoop == 0 && iter < maxsp2iter )
  {
    trX = bml_trace(rho_bml);

#ifdef DO_MPI
    if (bml_getNRanks() > 1)
    {
      startTimer(exchangeTimer);
//      exchangeSetup(dataExchange, xmatrix, domain);
      stopTimer(exchangeTimer);
    }
#endif

    // Matrix multiply X^2
    startTimer(x2Timer);
    trace = bml_multiply_x2(rho_bml, x2_bml, threshold);
    trX2 = trace[1];
    bml_free_memory(trace);
    stopTimer(x2Timer);

#ifdef DO_MPI
    // Reduce trace of X and X^2 across all processors
    if (bml_getNRanks() > 1)
    {
      startTimer(reduceCommTimer);
      addRealReduce2(&trX, &trX2);
      stopTimer(reduceCommTimer);
      collectCounter(reduceCounter, 2 * sizeof(real_t));
    }
#endif

    if (bml_printRank() && debug_i == 1) 
      printf("iter = %d  trX = %e  trX2 = %e\n", iter, trX, trX2);
 
    tr2XX2 = TWO*trX - trX2;
    trXOLD = trX;
    limDiff = ABS(trX2 - nocc) - ABS(tr2XX2 - nocc);

    if (limDiff > idemTol) 
    {
      // X = 2 * X - X^2
      trX = TWO * trX - trX2;

      startTimer(xaddTimer);
      bml_add(rho_bml, x2_bml, TWO, MINUS_ONE, threshold);
      stopTimer(xaddTimer);
    }
    else if (limDiff < -idemTol)
    {
      // X = X^2
      trX = trX2;

      startTimer(xsetTimer);
      bml_copy(x2_bml, rho_bml);
      stopTimer(xsetTimer);
    }
    else 
    {
      trX = trXOLD;
      breakLoop = 1;
    }
         
    idempErr2 = idempErr1;
    idempErr1 = idempErr;
    idempErr = ABS(trX - trXOLD);    

    iter++;

    if (iter >= minsp2iter && (idempErr >= idempErr2)) breakLoop = 1;

    // Exchange matrix pieces across processors
#ifdef DO_MPI
    if (bml_getNRanks() > 1)
    {
      startTimer(exchangeTimer);
      //exchangeData(dataExchange, xmatrix, domain);
      stopTimer(exchangeTimer);
    }
#endif
  }

  stopTimer(sp2LoopTimer);

  // Multiply by 2
  bml_scale_inplace(TWO, rho_bml);
  
  // Report results
  reportResults(iter, rho_bml, x2_bml);

#ifdef DO_MPI
  // Gather matrix to processor 0
  if (bml_getNRanks() > 1)
    //allGatherData(dataExchange, xmatrix, domain);
#endif

#ifdef DO_MPI
  if (bml_getNRanks() > 1)
    //destroyDataExchange(dataExchange);
#endif
  bml_deallocate(&x2_bml);
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
