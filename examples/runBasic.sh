echo "Nodes = " $1 
echo "Threads = " $2
#export OMPI_MCA_opal_paffinity_alone=0
export OMP_NUM_THREADS=$2
mpirun -np $1 --map-by node -x OMP_NUM_THREADS ./bin/ExaSP2-parallel-BASIC
