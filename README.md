This is ExaSP2 version 1.0
=========================

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

Copyright and license

ExaSP2 is a composable and extensible framework for SP2 algorithms.
A set of application programming interfaces (APIs) (under development) are 
available for SP2 solvers, data decomposition, data exchange, matrix math, 
and parallel communications. Default implementations are provided. New 
implementations can be added easily.

At build time, one can compose an SP2 test program by choosing the 
implementations for the SP2 solver (ex. BASIC), data decomposition (ex. ROW),
data exchange (ex. HALO), matrix math (ex. SPARSE), and parallel communication
(ex. MPI).

Required: C
Optional: MPI, BLAS library, other libraries

To build:

cd src
cp Makefile.vanilla Makefile

Modify the Makefile for the architecture and environment. 

For default:
make clean;make

For custom:
make clean;make SP2SOLVER=BASIC DECOMPOSITION=ROW, DATAEXCHANGE=HALO MATRIXMATH=SPARSE PARALLEL=MPI

To run default serial version: (generates random Hamiltonian)
export OMP_NUM_THREADS=16; ./bin/ExaSP2-serial

To run default parallel version: (generates random Hamiltonian)
export OMP_NUM_THREADS=16; mpirun -np 16 --map-by node -x OMP_NUM_THREADS ./bin/ExaSP2-parallel

To run parallel version with input arguments: (Hamiltonian from file)
export OMP_NUM_THREADS=16; mpirun -np 16 --map-by node -x OMP_NUM_THREADS ./bin/ExaSP2-parallel --hmatName ./data/hmatrix.1024.mtx --N 12288 --M 256
