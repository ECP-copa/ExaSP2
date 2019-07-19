#export OMP_NUM_THREADS=16;./bin/ExaSP2-serial-FERMI --hmatName data/nm.100.last.mtx --N 1900 --M 1900 --bndfil 0.63157894736842102 --beta 0.92552875303437621  --nsteps 18 --occLimit 1e-5

#export OMP_NUM_THREADS=8;~/ExaSP2/bin/ExaSP2-serial-IMP --hmatName ~/ExaSP2/data/poly_chain.1024.mtx --N 12288 --M 3000 --beta 4 --mu 0.2 --dout 1 --nsteps 10 --eps 1e-9 

#export OMP_NUM_THREADS=1;~/ExaSP2/bin/ExaSP2-serial-IMP --hmatName ~/ExaSP2/data/poly_chain.512.mtx --N 6144 --M 3000 --beta 1 --mu 0.2 --dout 1 --nsteps 10 --eps 1e-9

export OMP_NUM_THREADS=16;~/ExaSP2/bin/ExaSP2-serial-IMP --hmatName ~/ExaSP2/data/trpcage_8K.mtx --N 16863 --M 10000 --nsteps 10 --beta 10 --mu 0.2 --dout 1 --eps 1e-9
