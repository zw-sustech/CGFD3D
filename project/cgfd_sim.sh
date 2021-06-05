#!/bin/bash

set -e
printf "\nUse 4 CPUs on following nodes:\n"
cat /export/home/lihl/CGFD3D-elastic/project/hostlist

#-- simulation
printf "\nStart simualtion ...\n";
#time /share/apps/gnu-4.8.5/mpich-3.3/bin/mpiexec -f /export/home/lihl/CGFD3D-elastic/project/hostlist -np 4 /export/home/lihl/CGFD3D-elastic/fd3dpmltopo_ew 
time /share/apps/gnu-4.8.5/mpich-3.3/bin/mpiexec -machinefile /export/home/lihl/CGFD3D-elastic/project/hostlist -np 4 /export/home/lihl/CGFD3D-elastic/cgfdm3d_elastic_mpi /export/home/lihl/CGFD3D-elastic/project/test.json 100
if [ 0 -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi

