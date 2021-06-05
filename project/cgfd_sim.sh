#!/bin/bash

set -e
printf "\nUse 9 CPUs on following nodes:\n"
cat /export/home/wangyh/CGFD3D-elastic/project/hostlist

#-- simulation
printf "\nStart simualtion ...\n";
#time /export/apps/gnu-4.8.5/mpich-3.3/bin/mpiexec -f /export/home/wangyh/CGFD3D-elastic/project/hostlist -np 9 /export/home/wangyh/CGFD3D-elastic/fd3dpmltopo_ew 
time /export/apps/gnu-4.8.5/mpich-3.3/bin/mpiexec -machinefile /export/home/wangyh/CGFD3D-elastic/project/hostlist -np 9 /export/home/wangyh/CGFD3D-elastic/cgfdm3d_elastic_mpi /export/home/wangyh/CGFD3D-elastic/project/test.json 100
if [ 0 -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi

