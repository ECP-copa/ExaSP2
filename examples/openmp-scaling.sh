# Simple scaling study using polyethylene chain of 1024 atoms 
# and 12,288 orbitals using OpenMP
export OMP_NUM_THREADS=1; ./bin/ExaSP2-serial --hmatName data/hmatrix.1024.mtx --N 12288 --M 256
export OMP_NUM_THREADS=2; ./bin/ExaSP2-serial --hmatName data/hmatrix.1024.mtx --N 12288 --M 256
export OMP_NUM_THREADS=4; ./bin/ExaSP2-serial --hmatName data/hmatrix.1024.mtx --N 12288 --M 256
export OMP_NUM_THREADS=8; ./bin/ExaSP2-serial --hmatName data/hmatrix.1024.mtx --N 12288 --M 256
export OMP_NUM_THREADS=16; ./bin/ExaSP2-serial --hmatName data/hmatrix.1024.mtx --N 12288 --M 256
