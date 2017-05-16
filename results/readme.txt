This directory contains example results from running CoSP2 on the 2 polyethylene
datasets given in the data directory. Both runs use 2 MPI ranks and 16 OpenMPI
threads/rank. Note that additional parameters may be required to run on specific
architectures for affinity, etc. 

Polyethylene 512:
Run line - export OMP_NUM_THREADS=16; mpirun -np 2 ./CoSP2Parallel --hmatName data/hmatrix.512.mtx --N 6144 --M 256 --debug 1 --dout 1
Text output - poly_512_sp2.out
Output density matrix in Matrix Market format - dmatrix.512.out.mtx

Polyethylene 1024:
Run line - export OMP_NUM_THREADS=16; mpirun -np 2 ./CoSP2Parallel --hmatName data/hmatrix.1024.mtx --N 12288 --M 256 --debug 1 --dout 1
Text output - poly_1024_sp2.out 
Output density matrix in Matrix Market format - dmatrix.1024.out.mtx 

Note: These results are from the CoSP2 proxy app.
