#!/bin/bash

set -x

#. /share/apps/intel-2019.3/intel/bin/compilervars.sh intel64
#export FI_PROVIDER=tcp
#MPIDIR=/share/apps/intel-2019.3/mpich-3.3

MPIDIR=/export/apps/gnu-4.8.5/mpich-3.3

EXEC_DIR=/home/zhangw/code/zwlab/CGFD3D-elastic
EXEC_WAVE=$EXEC_DIR/main_curv_col_el_3d

PROJDIR=/home/zhangw/work/wpsfd_mpi/01

parfile=${PROJDIR}/test.json

#NUMPROCS=`grep ^num_threads_per_dim ${MAINCONF} | sed 's/=/ /g' | awk '{print $2 * $3 * $4}'`
NUMPROCS=2


#/share/apps/DDT/ddt-19.0.5/bin/ddt  \
#    $EXEC_DIR/seis3d_station $MAINCONF

/export/apps/DDT/ddt-19.0.5/bin/ddt  \
    $MPIDIR/bin/mpiexec -np $NUMPROCS $EXEC_WAVE $parfile 100
