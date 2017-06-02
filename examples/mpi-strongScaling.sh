# Simple strong scaling study using polyethylene chain of 1024 atoms 
# and 12,288 orbitals using MPI
export OMP_NUM_THREADS=1 
mpirun -np 1 -x OMP_NUM_THREADS ./bin/ExaSP2-parallel --hmatName ./data/hmatrix.1024.mtx --N 12288 --M 260 
mpirun -np 2 -x OMP_NUM_THREADS ./bin/ExaSP2-parallel --hmatName ./data/hmatrix.1024.mtx --N 12288 --M 260
mpirun -np 4 -x OMP_NUM_THREADS ./bin/ExaSP2-parallel --hmatName ./data/hmatrix.1024.mtx --N 12288 --M 260
mpirun -np 8 -x OMP_NUM_THREADS ./bin/ExaSP2-parallel --hmatName ./data/hmatrix.1024.mtx --N 12288 --M 260
mpirun -np 16 -x OMP_NUM_THREADS ./bin/ExaSP2-parallel --hmatName ./data/hmatrix.1024.mtx --N 12288 --M 260
