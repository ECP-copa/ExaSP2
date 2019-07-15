#export OMP_NUM_THREADS=16;./bin/ExaSP2-serial-FERMI --hmatName data/nm.100.last.mtx --N 1900 --M 1900 --bndfil 0.63157894736842102 --beta 0.92552875303437621  --nsteps 18 --occLimit 1e-5

#export OMP_NUM_THREADS=16;./bin/ExaSP2-serial-FERMI --hmatName data/nm.200.last.mtx --N 3800 --M 3800 --bndfil 0.63157894736842102 --beta 0.90663935914396565 --nsteps 18 --occLimit 1e-5

#export OMP_NUM_THREADS=16;./bin/ExaSP2-serial-FERMI --hmatName data/poly_chain.1024.mtx --N 12288 --M 1000 --nsteps 20 --occLimit 1e-5

export OMP_NUM_THREADS=1;~/ExaSP2/bin/ExaSP2-serial-FERMI --hmatName ~/ExaSP2/data/F_ort_matrix_0008.mtx --N 1204 --beta 0.5 --nocc 0.2 --dout 1 
