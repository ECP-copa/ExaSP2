/// \file
/// SP2 main program.
///
/// \mainpage ExaSp2: A Linear Algebra Electronic Structure Mini-app
///
/// ExaSp2 is a reference implementation of typical sparse linear algebra
/// algorithms and workloads for a quantum molecular dynamics (QMD)
/// electronic structure code, such as the Los Alamos Transferable
/// Tight-binding for Energetics (LATTE) code. The algorithm is based
/// on a recursive second-order Fermi-Operator expansion method (SP2)
/// and is tailored for density functional based tight-binding calculations
/// of non-metallic material systems. It is composed of a series of
/// generalized matrix-matrix multiplications. It is created and maintained by
/// The Co-Design Center for Particle Application
/// (CoPA).  https://github.com/ECP-copa.  The
/// code is intended to serve as a vehicle for co-design by allowing
/// others to extend and/or reimplement it as needed to test performance of
/// new architectures, programming models, etc.
///
/// The current version of ExaSp2 is available from:
/// https://github.com/ECP-copa/ExaSP2
///
/// Table of Contents
/// =================
///
/// Click on the links below to browse the ExaSp2 documentation.
///
/// \subpage pg_qmd_basics
///
/// \subpage pg_building_cosp2
///
/// \subpage pg_running_cosp2
///
/// \subpage pg_measuring_performance
///
/// \subpage pg_problem_selection_and_scaling
///
/// \subpage pg_cosp2_architecture
///
/// \subpage pg_optimization_targets
///

#define MAIN_FILE

#include "bml.h"

#include <stdio.h>
#include <omp.h>

#include "sp2Solver.h"
#include "parallel.h"
#include "performance.h"
#include "mycommand.h"
#include "constants.h"

/// \details
/// Adjust number of non-zeroes
int nnzStart(const int hsize,
             const int msize)
{
  int M = msize;
  if (M == 0) M = hsize;
  if ((M % 32) > 0) M += (32 - (M % 32));
  if (M > hsize) M = hsize;
  if (bml_printRank()) printf("Adjusted M = %d\n", M);

  return M;
}

/// \details
/// Initialize h matrix
bml_matrix_t* initSimulation(const Command cmd)
{
  bml_matrix_t* h_bml;

  bml_matrix_type_t matrix_type = cmd.mtype;
  bml_matrix_precision_t precision = double_real;
  bml_distribution_mode_t dmode = sequential;

  // Read in size of hamiltonian matrix and max number of non-zeroes
  if (bml_printRank()) printf("N = %d M = %d\n", N_i, msparse_i);

  // Calculate M - max number of non-zeroes per row
  M_i = nnzStart(N_i, msparse_i);  

  if (cmd.gen == 0)
  {
    // Allocate  input hamiltonian matrix
    h_bml = bml_zero_matrix(matrix_type, precision, N_i, M_i, dmode);
    startTimer(readhTimer);
    bml_read_bml_matrix(h_bml, cmd.hmatName);
    stopTimer(readhTimer);
  }

  else
  {
    // Banded Hamiltonian is generated
    h_bml = bml_banded_matrix(matrix_type, precision, N_i, M_i, dmode);
  }

  return h_bml;
}


int main(int argc,
         char** argv)
{
  // Start
  bml_init(&argc, &argv);
  profileStart(totalTimer);
  profileStart(loopTimer);
  if (bml_printRank()) printf("ExaSp2: SP2 Loop\n");

  // Read in command line parameters
  Command cmd = parseCommandLine(argc, argv);
  msparse_i = cmd.M;
  N_i = cmd.N;
  mtype_i = cmd.mtype;
  minsp2iter_i = cmd.minsp2iter;
  maxsp2iter_i = cmd.maxsp2iter;
  nsteps_i = cmd.nsteps;
  debug_i = cmd.debug;
  dout_i = cmd.dout;

  nocc_i = cmd.nocc;
  eps_i = cmd.eps;
  idemTol_i = cmd.idemTol;
  bndfil_i = cmd.bndfil;
  mu_i = cmd.mu;
  beta_i = cmd.beta; 
  tscale_i = cmd.tscale;
  occLimit_i = cmd.occLimit;
  traceLimit_i = cmd.traceLimit;

  if (bml_printRank())
  {
    printf("\nParameters:\n");
    printf("msparse = %d  N = %d\n", msparse_i, N_i);
    printf("minsp2iter = %d  maxsp2iter = %d\n", minsp2iter_i, maxsp2iter_i);
    printf("nsteps = %d  osteps = %d\n", nsteps_i, osteps_i);
    printf("mtype = %d  hmatName = %s\n", cmd.mtype, cmd.hmatName);
    printf("eps = %lg  tscale = %lg\n", eps_i, tscale_i);
    printf("idemTol = %lg  bndfil = %lg  beta = %lg  mu = %lg\n", idemTol_i, bndfil_i, beta_i, mu_i);
    printf("occLimit = %lg  traceLimit = %lg\n", occLimit_i, traceLimit_i);
    printf("debug = %d  dout = %d\n\n", debug_i, dout_i);
  }

  // Initialize
  startTimer(preTimer);
  bml_matrix_t* h_bml = initSimulation(cmd);
  bml_matrix_type_t matrix_type = bml_get_type(h_bml);
  bml_matrix_precision_t precision = bml_get_precision(h_bml);
  bml_distribution_mode_t dmode = bml_get_distribution_mode(h_bml);
  bml_matrix_t* rho_bml = bml_zero_matrix(matrix_type, precision, N_i, M_i, dmode);
  stopTimer(preTimer);

  // Run SP2 variant
  if (nocc_i <= 0.0)
  {
    nocc_i = bndfil_i * N_i;
  }
  printf("nocc = %lg\n", nocc_i);

#ifdef SP2_IMP
  printf("Calling Implicit Fermi\n"); 
  implicit_recursiveLoops(h_bml, rho_bml, beta_i, mu_i, nsteps_i, eps_i);
#endif

#ifdef SP2_BASIC
  printf("Calling Basic\n");
  // Perform SP2 loop
  sp2Loop(h_bml, rho_bml, nocc_i, minsp2iter_i, maxsp2iter_i, idemTol_i, eps_i);
#endif

#ifdef SP2_FERMI
  printf("Calling Fermi\n");
  real_t mu = ZERO;
  real_t beta = beta_i;
  int* sgnlist = bml_allocate_memory(nsteps_i*sizeof(int));
  real_t h1 = ZERO;
  real_t hN = ZERO;
  real_t kbt = ZERO;

  // Perform truncated SP2 Fermi initialization followed by Fermi
  printf("sp2Init start: mu = %lg beta = %lg \n", mu, beta);
  startTimer(sp2InitTimer);
  sp2Init(h_bml, rho_bml, nsteps_i, nocc_i, &mu, &beta, sgnlist, &h1, &hN,
    tscale_i, occLimit_i, traceLimit_i, eps_i); 
  stopTimer(sp2InitTimer);

  kbt = ABS(ONE / beta);
  printf("sp2Init complete: mu = %lg beta = %lg kbt = %lg\n", mu, beta, kbt);

  startTimer(sp2LoopTimer);
  sp2Loop(h_bml, rho_bml, nsteps_i, nocc_i, &mu, beta, sgnlist, h1, hN,
    osteps_i, eps_i, traceLimit_i, eps_i);
  stopTimer(sp2LoopTimer);

  printf("sp2Loop complete: mu = %lg beta = %lg kbt = %lg\n", mu, beta, kbt);

  bml_free_memory(sgnlist);

#endif
  // Done
  profileStop(totalTimer);
  profileStop(loopTimer);

  /// Show timing results
  printPerformanceResults(N_i, 0);

  /// Write out density matrix
  if (bml_printRank() && dout_i == 1)
  {
    bml_write_bml_matrix(rho_bml, "dmatrix.out.mtx");
  }

  /// Deallocate matrices, etc.
  bml_deallocate(&h_bml);
  bml_deallocate(&rho_bml);

  bml_shutdown();

  return 0;
}

// --------------------------------------------------------------
//
//
/// \page pg_building_exasp2 Building ExaSP2
///
/// ExaSP2 is written with portability in mind and should compile using
/// practically any compiler that implements the C99 standard.  You will
/// need to create a Makefile by copying the sample provided with the
/// distribution (Makefile.vanilla).
///
///     $ cp Makefile.vanilla Makefile
///
/// and use the make command to build the code
///
///    $ make
///
/// The sample Makefile will compile the code on many platforms.  See
/// comments in Makefile.vanilla for information about specifying the
/// name of the C compiler, and/or additional compiler switches that
/// might be necessary for your platform.
///
/// The main options available in the Makefile are toggling single/double
/// precision and enabling/disabling MPI. In the event MPI is not
/// available, setting the DO_MPI flag to OFF will create a purely
/// serial build (you will likely also need to change the setting of
/// CC).
///
/// The makefile should handle all the dependency checking needed, via
/// makedepend.
///
/// 'make clean' removes the object and dependency files.
///
/// 'make distclean' additionally removes the executable file and the
/// documentation files.

//--------------------------------------------------------------


// --------------------------------------------------------------


/// \page pg_measuring_performance Measuring Performance
///
/// ExaSP2 implements a simple and extensible system of internal timers
/// and counters to measure the performance profile of the code.
/// As explained in performance.c, it is easy to create additional
/// timers and counters and associate them with code regions of specific
/// interest.  In addition, the getTime() and getTick() functions can be
/// easily reimplemented to take advantage of platform specific timing
/// resources.
///
/// A performance report is printed at the end of each simulation.
///
///
/// Counters for Rank 0
///        Counter          Calls    Avg/Call(MB)         Total(MB)
/// _________________________________________________________________
/// reduce                    35           0.0000            0.0005
/// send                      75           1.0479           78.5938
/// recv                      75           1.0485           78.6367
///
/// Counter Statistics Across 16 Ranks:
///         Counter      Rank: Min(MB)        Rank: Max(MB)       Avg(MB)      Stdev(MB)
/// _______________________________________________________________________________________
/// reduce               0:      0.0005       0:      0.0005        0.0005        0.0000
/// send                 0:     78.5938      14:     78.6784       78.6164        0.0313
/// recv                 4:     78.5955      14:     78.6768       78.6164        0.0251
///
///
/// Timings for Rank 0
///         Timer        # Calls    Avg/Call (s)   Total (s)    % Loop
/// ___________________________________________________________________
/// total                      1       0.7124        0.7124      100.00
/// loop                       1       0.7124        0.7124      100.00
///   pre                      1       0.1119        0.1119       15.71
///     readh                  1       0.0945        0.0945       13.27
///   sp2Loop                  1       0.5505        0.5505       77.28
///     norm                   1       0.0028        0.0028        0.40
///     x2                    31       0.0114        0.3527       49.51
///     xadd                  14       0.0011        0.0150        2.10
///     xset                  17       0.0010        0.0162        2.27
///     exchange              62       0.0022        0.1389       19.50
///     reduceComm            35       0.0005        0.0158        2.22
///
/// Timing Statistics Across 16 Ranks:
///         Timer        Rank: Min(s)       Rank: Max(s)      Avg(s)    Stdev(s)
/// _____________________________________________________________________________
/// total                3:    0.7068       8:    0.7239      0.7157      0.0042
/// loop                 3:    0.7068       8:    0.7239      0.7157      0.0042
///   pre               12:    0.1105       2:    0.1127      0.1115      0.0006
///     readh            4:    0.0944       8:    0.0951      0.0947      0.0002
///   sp2Loop            6:    0.5501       1:    0.5508      0.5504      0.0002
///     norm             2:    0.0026      14:    0.0030      0.0028      0.0001
///     x2               7:    0.3511       2:    0.3551      0.3531      0.0010
///     xadd            10:    0.0145      14:    0.0154      0.0149      0.0002
///     xset             9:    0.0155       2:    0.0166      0.0159      0.0003
///     exchange        10:    0.1338       1:    0.1394      0.1367      0.0017
///     reduceComm      14:    0.0081      10:    0.0219      0.0173      0.0039
///
///
/// This report consists of two blocks each for counters and timers.  The
/// first block for timers lists the absolute wall clock time spent in each
/// timer on rank 0 of the job. The first block for counters lists the sizes
/// of communication messages on rank 0 of the job. The second block for
/// timers and counters reports minimum, maximum, average, and standard
/// deviation of times or sizes across all tasks. The ranks where the
/// minimum and maximum values occured are also reported to aid in identifying
/// hotspots or load imbalances.
///
/// *** Architecture/Configuration for above timing numbers:
/// A cluster with dual-socket 8-core Intel Xeon SandyBridge E5-2670 processors.
/// Each node has a total of 16 cores and shares 64 GB of RAM.
///

// --------------------------------------------------------------


/// \page pg_problem_selection_and_scaling Problem Selection and Scaling
///
/// ExaSP2 is a reference linear algebra electronic structure simulation code as used in
/// materials science.
///
/// Problem Specification
/// ======================
///
/// ExaSP2 requires an input Hamiltonian matrix. A generator is included that produces
/// a synthetic banded matrix similar to polyethylene. N (number of orbitals or rows)
/// and M (number of non-zeroes per row) must be specified. The extent of the band is
/// determined by N and M. SP2 matrix multiplication adds non-zero elements per row.
/// You may have to experiment to get a reasonable M.
///
/// Given a user-specified amplitude and alpha (--amp and --alpha command line options):
///
/// value = amplitude * random[0,1] * exp(-alpha * (i-j)^2)
///
/// A more general version based on atom positions:
///
/// value = amplitude * random[0,1] * exp(-alpha * (Ri-Rj)^2)
///         where (Ri - Rj) is the distance between atoms
///
/// Two Hamiltonian matrix examples are given in data/ that represent polyethylene
/// chains of 512 and 1024 molecules. A good M value (or number of non-zeroes per
/// row) is 256 for each.
///
/// Scaling Studies in ExaSP2
/// =======================
///
/// ExaSP2 implements a simple domain decomposition dividing the Hamiltonian matrix
/// of N rows into sub-matrices of chunks of rows, each owned by an MPI rank.
/// Each domain is a single-program multiple data (SPMD) partition of the larger problem.
///
/// Weak Scaling
/// -----------
///
/// A weak scaling test fixes the amount of work per processor and
/// compares the execution time over number of processors. Weak scaling
/// keeps the ratio of inter-processor communication (surface) to
/// intra-processor work (volume) fixed. The amount of inter-processor
/// work scales with the number of processors in the domain. Run increasingly
/// larger problems over more processors, ex. size N on 1 MPI rank, 2N on
/// 2 MPI ranks, etc.
///
///
/// Strong Scaling
/// ---------------
///
/// A strong scaling test fixes the total problem size and compares the
/// execution time for different numbers of processors. Strong scaling
/// increases the ratio of inter-processor communication (surface) to
/// intra-processor work (volume). Run with increased numbers of processors,
/// ex. 1,2,4,8,16 MPI ranks.
///
///
/// \page pg_exasp2_architecture ExaSP2 Architecture
///
/// Program Flow
/// ============
///
/// We have attempted to make the program flow in ExaSP2 1.0 as simple and
/// transparent as possible.  The main program consists of initialization,
/// the sparse SP2 loop, and end processing.
///
/// The initialization initializes MPI (if requested), reads in command line options,
/// reads or generates the sparse Hamiltonian matrix. The SP2 loop is as follows.
///
/// SP2 Algorithm:
///
///   breakLoop = 0
///   iter = 0
///   idempErr = idempErr1 = idempErr2 = 0
///
///   do while (breakLoop == 0 && iter < 100)
///   {
///     trX = 0
///     trX2 = 0
///
///     exchange data information (if running MPI)
///
///     X2 = X * X
///
///     reduce trX and trX2 (if running MPI)
///
///     tr2XX2 = 2.0 * trX - trX2
///     trXOLD = trX
///     limDiff = |trX2 - occ| - |tr2XX2 - occ|
///
///     if (limDiff > idemTol)
///     {
///       X = 2.0 * X - X2
///       trX = 2.0 * trX - trX2
///     }
///     else
///     {
///       X = X2
///       trX = trX2
///     }
///
///     idempErr2 = idempErr1
///     idempErr1 = idempErr
///     idempErr = |trX - trXOLD|
///
///     iter++
///
///     if (iter >= 25 && (idempErr >= idempErr2)) breakLoop = 1
///
///     exchange relevant sparse matrix chunks of rows (if running MPI)
///
///   }
///
/// End processing consists of printing performance results and deleting
/// structures used.
///
/// Key Data Structures
/// ==================
/// Data structures used are the following.
/// - SparseMatrixSt Sparse matrix storage structure based on ELLPACK-R.
/// - DomainSt Decomposition of chunks of rows across MPI ranks.
/// - DataExchangeSt Structure for communication of data.
///

/// \page pg_optimization_targets Optimization Targets
///
/// Computation
/// ============
///
/// The computational effort of electronic structure codes is usually focused on the
/// density matrix build. Non-metalic systems produce sparse Hamiltonian matrices.
/// The SP2 algorithm takes advantage of efficient sparse matrix-matrix multiplication
/// and addition. Using the ELLPACK-R storage structure and multi-threading using
/// OpenMP, results in sub-second performance.
///
/// Communication
/// =============
///
/// As the number of atoms per MPI rank decreases, the communication
/// routines will start to require a significant fraction of the
/// run time.  The main communication routine in ExaSP2 is dataExchange().
/// The halo data exchange allows communication of matrix row chunks between
/// ranks that require them for computation. All-to-all is not required
/// to do the SP2 algorithm.
///
/// Initially, all ranks have the entire Hamiltonian matrix.
/// Information on which chunks are to be exchanged is communicated at
/// the start of each SP2 iteration. Chunks are exchanged at the end of
/// each iteration. Once the SP2 algorithm is done, the remaining chunks
/// are exchanged, so each rank contains the entire new density matrix.
///

// --------------------------------------------------------------


/// \page pg_qmd_basics QMD Basics
///
/// The molecular dynamics (MD) computer simulation method is a well
/// established and important tool for the study of the dynamical
/// properties of liquids, solids, and other systems of interest in
/// Materials Science and Engineering, Chemistry and Biology. A material
/// is represented in terms of atoms and molecules. The method of MD
/// simulation involves the evaluation of the force acting on each atom
/// due to all other atoms in the system and the numerical integration
/// of the Newtonian equations of motion. Though MD was initially
/// developed to compute the equilibrium thermodynamic behavior of
/// materials (equation of state), most recent applications have used MD
/// to study non-equilibrium processes.
///
/// Wikipeda offers a basic introduction to molecular dynamics with
/// many references:
///
/// http://en.wikipedia.org/wiki/Molecular_dynamics
///
/// Quantum Molecular Dynamics (QMD) simulations capture the making and
/// breaking of covalent bonds, charge transfer between species of
/// differing electronegativities, and long-range electrostatic interactions.
/// Electronic structure of atoms and molecules is modeled explicitly.
/// It provides the most accurate and reliable descriptions of interatomic bonding.
/// The prohibitive computational cost has prevented widespread use. Today,
/// better algorithms, and multi-core and GPU implementations are important
/// paths forward.
///
/// Wikipedia background on computational chemistry:
///
/// http://en.wikipedia.org/wiki/Computational_chemistry
///
/// References:
///
/// E. H. Rubensson, A. M. N. Niklasson, 2014,
/// Interior eigenvalues from density matrix expansions in quantum mechanical
/// molecular dynamics, SIAM J. Sci. Comput. Vol. 36, No. 2, pp. B147–B170.
///
/// P. Souvatzis, A. M. N. Niklasson, 2014,
/// First principles molecular dynamics without self-consistent field
/// optimization, J. Chem. Physics 140, 044117.
///
/// M. J. Cawkwell, E. J. Sanville, S. M. Mniszewski, A. M. N. Niklasson, 2012,
/// Computing the density matrix in electronic structure theory on
/// graphics processing units, J. Chem. Theory Comput., 8, 4094−4101.
///
/// M. J. Cawkwell, A. M. N. Niklasson, 2012,
/// Energy conserving, linear scaling Born-Oppenheimer molecular dynamics,
/// J. Chem. Physics 137, 134105.
///
/// A. M. N. Niklasson, P., Steneteg, N. Bock, 2011,
/// Extended Lagrangian free energy molecular dynamics,
/// J. Chem. Physics 135, 164111.
///
/// A. M. N. Niklasson, 2008,
/// Extended Born-Oppenheimer molecular dynamics,
/// PRL 100, 123004.
///
/// A. M. N. Niklasson, 2002,
/// Expansion algorithm for the density matrix,
/// Phys. Rev. B 66, 155115.
///
