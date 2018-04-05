/// \File
/// Contains frequently used constants

#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <math.h>

#include "mytype.h"

#define ABS fabs

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) < (B) ? (A) : (B))

/// Numerical values used
static const real_t MINUS_ONE = -1.0;
static const real_t ZERO = 0.0;
static const real_t ONE = 1.0;
static const real_t TWO = 2.0;
static const real_t THREE = 3.0;
static const real_t HALF = 0.5;
static const real_t MINUS_THOUSAND = -1000.0;

static const int HUNDRED = 100;

#ifdef MAIN_FILE
int msparse_i;
int N_i;
int M_i;
int mtype_i;
int minsp2iter_i;
int maxsp2iter_i;
int nsteps_i;
int osteps_i;
int debug_i;

real_t nocc_i; 
real_t eps_i; 
real_t idemTol_i; 
real_t bndfil_i;
real_t beta_i;
real_t tscale_i;
real_t occLimit_i;
real_t traceLimit_i;
#else
extern int msparse_i;
extern int N_i;
extern int M_i;
extern int mtype_i;
extern int minsp2iter_i;
extern int maxsp2iter_i;
extern int nsteps_i;
extern int osteps_i;
extern int debug_i;

extern real_t nocc_i;          
extern real_t eps_i;          
extern real_t idemTol_i;
extern real_t bndfil_i;
extern real_t beta_i;
extern real_t tscale_i;
extern real_t occLimit_i;
extern real_t traceLimit_i;
#endif

#endif
