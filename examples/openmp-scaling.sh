# Simple scaling study using polyethylene chain of 1024 atoms 
# and 12,288 orbitals using OpenMP
echo
echo "Nodes = 1  Threads = 1"
export OMP_NUM_THREADS=1; ./bin/ExaSP2-serial-BASIC --hmatName data/poly_chain.1024.mtx --N 12288 --M 260
echo
echo "Nodes = 1  Threads = 2"
export OMP_NUM_THREADS=2; ./bin/ExaSP2-serial-BASIC --hmatName data/poly_chain.1024.mtx --N 12288 --M 260
echo
echo "Nodes = 1  Threads = 4"
export OMP_NUM_THREADS=4; ./bin/ExaSP2-serial-BASIC --hmatName data/poly_chain.1024.mtx --N 12288 --M 260
echo
echo "Nodes = 1  Threads = 8"
export OMP_NUM_THREADS=8; ./bin/ExaSP2-serial-BASIC --hmatName data/poly_chain.1024.mtx --N 12288 --M 260
echo
echo "Nodes = 1  Threads = 16"
export OMP_NUM_THREADS=16; ./bin/ExaSP2-serial-BASIC --hmatName data/poly_chain.1024.mtx --N 12288 --M 260
