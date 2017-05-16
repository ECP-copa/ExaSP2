# Simple scaling study using polyethylene chain of 1024 atoms and 
# 12,288 orbitals using MPI + OpenMP.
export OMPI_MCA_opal_paffinity_alone=0
export OMP_NUM_THREADS=16
mpirun -np 1 --map-by node -x OMP_NUM_THREADS ../bin/ExaSP2-parallel --hmatName data/hmatrix.1024.mtx --N 12288 --M 256
mpirun -np 2 --map-by node -x OMP_NUM_THREADS ../bin/ExaSP2-parallel --hmatName data/hmatrix.1024.mtx --N 12288 --M 256
mpirun -np 4 --map-by node -x OMP_NUM_THREADS ../bin/ExaSP2-parallel --hmatName data/hmatrix.1024.mtx --N 12288 --M 256
mpirun -np 8 --map-by node -x OMP_NUM_THREADS ../bin/ExaSP2-parallel --hmatName data/hmatrix.1024.mtx --N 12288 --M 256
mpirun -np 16 --map-by node -x OMP_NUM_THREADS ../bin/ExaSP2-parallel --hmatName data/hmatrix.1024.mtx --N 12288 --M 256
