# Simple scaling study using polyethylene chain of 1024 atoms 
# and 12,288 orbitals using OpenMP
echo
echo "Nodes = 1  Threads = 1"
export OMP_NUM_THREADS=1; ~/ExaSP2/bin/ExaSP2-serial-IMP --hmatName ~/ExaSP2/data/poly_chain.1024.mtx --N 12288 --M 3000 --beta 4 --mu 0.2 --dout 1 --nsteps 10 --eps 1e-9
echo
echo "Nodes = 1  Threads = 2"
export OMP_NUM_THREADS=2; ~/ExaSP2/bin/ExaSP2-serial-IMP --hmatName ~/ExaSP2/data/poly_chain.1024.mtx --N 12288 --M 3000 --beta 4 --mu 0.2 --dout 1 --nsteps 10 --eps 1e-9
echo
echo "Nodes = 1  Threads = 4"
export OMP_NUM_THREADS=4; ~/ExaSP2/bin/ExaSP2-serial-IMP --hmatName ~/ExaSP2/data/poly_chain.1024.mtx --N 12288 --M 3000 --beta 4 --mu 0.2 --dout 1 --nsteps 10 --eps 1e-9
echo
echo "Nodes = 1  Threads = 8"
export OMP_NUM_THREADS=8; ~/ExaSP2/bin/ExaSP2-serial-IMP --hmatName ~/ExaSP2/data/poly_chain.1024.mtx --N 12288 --M 3000 --beta 4 --mu 0.2 --dout 1 --nsteps 10 --eps 1e-9
echo
echo "Nodes = 1  Threads = 16"
export OMP_NUM_THREADS=16; ~/ExaSP2/bin/ExaSP2-serial-IMP --hmatName ~/ExaSP2/data/poly_chain.1024.mtx --N 12288 --M 3000 --beta 4 --mu 0.2 --dout 1 --nsteps 10 --eps 1e-9
