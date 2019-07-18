# ExaSP2 

C17099

ExaSP2 is a reference implementation of typical linear algebra
algorithms and workloads for a quantum molecular dynamics (QMD)
electronic structure code. The algorithm is based on a recursive
second-order Fermi-Operator expansion method (SP2) and is tailored
for density functional based tight-binding calculations of
material systems. The SP2 algorithm variants are part of the Los Alamos 
Transferable Tight-binding for Energetics (LATTE) code, based on 
a matrix expansion of the Fermi operator in a recursive series 
of generalized matrix-matrix multiplications. 
It is created and maintained by Co-Design Center for Particle 
Applications (CoPA). The code is intended
to serve as a vehicle for co-design by allowing others to extend 
and/or reimplement as needed to test performance of new 
architectures, programming models, etc.

ExaSP2 has been written to use the data structures and matrix operations
in the Basic Matrix Library (BML). It is a composable and extensible 
framework for SP2 algorithms. A set of application programming interfaces
(APIs) (under development) are available for SP2 solvers, data decomposition
(BML), data exchange (BML), matrix operations (BML), and parallel 
communications. New implementations can be added easily. New data structures, 
decompositions, and matrix operations are added through the BML.

At build time, one can compose an SP2 test program by choosing the 
implementation for the SP2 solver (ex. BASIC) and parallel communication 
(ex. MPI). At run time choices can be made as input parameters for the matrix
type (ex. dense, sparse), data decomposition (ex. 1-D, 2-D, GRAPH), and data 
exchange (ex. HALO).

# Current Capabilities

## Electronic Structure Solvers:
 * BASIC    - original SP2 algorithm for calculation of the density matrix at zero temperature (default)
 * FERMI    - truncated SP2 for finite temperature and partial occupation
 * IMP      - implicit recursive expansion for finite temperature and partial occupation
 * RESPONSE - quantum perturbation theory for response properties at zero and finite temperature (future)
 * KERNEL   - fast Newton solver for non-linear equations (future))
 * CHEBYSHEV- Chebyshev kernel polynomial method (future)

## Data Decomposition:
 * 1-D   - chunks of rows/columns (default)
 * 2-D   - blocks (future)
 * GRAPH - graph-partitioned sub-matrices (future)

## Data Exchange:
 * HALO - exchange halo data (future)

## Matrix Type: representations and operations
 * SPARSE - using BML ELLPACK format (default)
 * DENSE  - using BML DENSE format

## Parallel Communication:
 * NONE - serial (default)
 * MPI  - parallel

## Hamiltonian Matrix Input:
 * In *.mtx format

# Compilation

## Dependencies
 * Required: C
 * Required: BML (https://github.com/lanl/bml)
 * Required: BLAS/LAPACK (ex. MKL)
 * Optional: MPI, other libraries (per implementation choices)

Build in the src directory.
Copy and modify Makefile for the architecture and environment.
```
cp Makefile.vanilla Makefile
```

Build default version.
```
make clean;make
```

Build custom version.
```
make clean;make SP2SOLVER=BASIC PARALLEL=MPI
```

# Running

Run the default serial version: (generates random sparse Hamiltonian)
```
export OMP_NUM_THREADS=16; ./bin/ExaSP2-serial
```
Run the default parallel version: (generates random sparse Hamiltonian)
```
export OMP_NUM_THREADS=16; mpirun -np 16 --map-by node -x OMP_NUM_THREADS ./bin/ExaSP2-parallel
```

Run the parallel version with input arguments: (Hamiltonian from file, sparse matrix type)
```
export OMP_NUM_THREADS=16; mpirun -np 16 --map-by node -x OMP_NUM_THREADS ./bin/ExaSP2-parallel --hmatName ./data/hmatrix.1024.mtx --mtype 2 --N 12288 --M 256
```
