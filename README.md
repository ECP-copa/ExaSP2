# ExaSP2 

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

ExaSP2 is a composable and extensible framework for SP2 algorithms.
A set of application programming interfaces (APIs) (under development) are 
available for SP2 solvers, data decomposition, data exchange, matrix math, 
and parallel communications. Default implementations are provided. New 
implementations can be added easily.

At build time, one can compose an SP2 test program by choosing the 
implementations for the SP2 solver (ex. BASIC), data decomposition (ex. ROW),
data exchange (ex. HALO), matrix math (ex. SPARSE), and parallel communication
(ex. MPI).

# Current Capabilities

## SP2 Solvers:
 * BASIC    - original SP2 algorithm for calculation of the density matrix at zero temperature (default)
 * FERMI    - truncated SP2 for finute temperature and partial occupation (future)
 * RESPONSE - quantum perturbation theory for response properties at zero and finite temperature (future)
 * KERNEL   - fast Newton solver for non-linear equations (future))
 * CHEBYSHEV- Chebyshev polynomial (future)

## Data Decomposition:
 * ROW   - chunks of rows (default)
 * GRAPH - graph-partitioned sub-matrices (future)

## Data Exchange:
 * HALO - exchange halo data (default)

## Matrix Math: representations and operations
 * SPARSE - only sparse (default)
 * BML    - Basic Matrix Library (future)

## Parallel Communication:
 * NONE - serial (default)
 * MPI  - parallel

## Hamiltonian Matrix Input:
 * In *.mtx format

# Compilation

## Dependencies
 * Required: C
 * Optional: MPI, BLAS library, other libraries (per implementation choices)

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
make clean;make SP2SOLVER=BASIC DECOMPOSITION=ROW, DATAEXCHANGE=HALO MATRIXMATH=SPARSE PARALLEL=MPI
```

# Running

Run the default serial version: (generates random Hamiltonian)
```
export OMP_NUM_THREADS=16; ./bin/ExaSP2-serial
```
Run the default parallel version: (generates random Hamiltonian)
```
export OMP_NUM_THREADS=16; mpirun -np 16 --map-by node -x OMP_NUM_THREADS ./bin/ExaSP2-parallel
```

Run the parallel version with input arguments: (Hamiltonian from file)
```
export OMP_NUM_THREADS=16; mpirun -np 16 --map-by node -x OMP_NUM_THREADS ./bin/ExaSP2-parallel --hmatName ./data/hmatrix.1024.mtx --N 12288 --M 256
```
