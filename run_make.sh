#!/bin/bash

#-- use module to set mpicc etc, on Mars
#module load intel/2019.5
#module load mpi/mpich/3.3.1_intel_2019.5
#module load netcdf-c/4.4.1

#-- add mpi to PATH if module is not used, on server1
MPI_ROOT=/share/apps/gnu-4.8.5/mpich-3.3/
export PATH=$MPI_ROOT/bin:$PATH
export NETCDFROOT=/share/apps/gnu-4.8.5/disable-netcdf-4.4.1

echo "mpicc and mpicxx will invoke:"
mpicc -show
mpicxx -show

echo
echo "start to make ..."
make -f Makefile 
