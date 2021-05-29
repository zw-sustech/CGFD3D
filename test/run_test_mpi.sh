#!/bin/bash

#set -x
set -e

date

#-- source intel lib
#source /share/apps/intel-2019.3/intel/bin/compilervars.sh intel64
#export FI_PROVIDER=tcp
#MPIDIR=/share/apps/intel-2019.3/mpich-3.3

#-- sever related dir
MPIDIR=/export/apps/gnu-4.8.5/mpich-3.3
EXEC_MATLAB=/export/apps/Matlab/R2015b/bin/matlab2015

#-- user related dir
EXEC_DIR=/home/zhangw/code/zwlab/CGFD3D-elastic
EXEC_WAVE=$EXEC_DIR/cgfdm3d_elastic_mpi

#-- project related dir
PROJDIR=/export/home/zhangw/work/cgfd_opt/11_mpi
export PROJDIR

par_file=${PROJDIR}/test.json
output_dir=${PROJDIR}/output
grid_dir=${PROJDIR}/output
media_dir=${PROJDIR}/output

#-- create dir
mkdir -p $PROJDIR
mkdir -p $output_dir
mkdir -p $grid_dir
mkdir -p $media_dir

#----------------------------------------------------------------------
#-- create hostlist for mpirun
#----------------------------------------------------------------------
cat << ieof > ${PROJDIR}/hostlist
server1
ieof

#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > $par_file
{
  "grid_name" : "blk_1",
  "number_of_total_grid_points_x" : 200,
  "number_of_total_grid_points_y" : 200,
  "number_of_total_grid_points_z" : 60,

  "number_of_mpiprocs_x" : 2,
  "number_of_mpiprocs_y" : 2,

  "size_of_time_step" : 0.01,
  "number_of_time_steps" : 500,

  "boundary_condition" : {
      "x_left"   : "cfspml",
      "x_right"  : "cfspml",
      "y_front"  : "cfspml",
      "y_back"   : "cfspml",
      "z_bottom" : "cfspml",
      "z_top"    : "free"
  },

  "#boundary_condition" : [
      "none",
      "none",
      "none",
      "none",
      "none",
      "free"
  ],

  "cfspml" : {
      "number_of_layers" : [
        5,5,5,5,5,0
      ],
      "alpha_max" : [
        0.0,0.0,0.0,0.0,0.0,0.0
      ],
      "beta_max" : [
        1.0,1.0,1.0,1.0,1.0,1.0
      ],
      "pml_velocity" : [
        3000.0, 3000.0, 3000.0, 3000.0, 3000.0, 3000.0
      ]
  },

  "#-- input_grid_type" : "vmap cartesian grid import",
  "input_grid_type" : "cartesian",
  "cartesian_grid_origin"   : [ 0.0, 0.0, -6000.0 ],
  "cartesian_grid_stepsize" : [ 100.0, 100.0, 100.0 ],
  "#input_vmap_file" : "test.vmap",

  "#-- input_metric_type" : "calculate grid import",
  "input_metric_type" : "calculate",
  "#input_metric_file" : "test.sta",

  "#-- input_medium_type" : "grid layer import infunc",
  "input_medium_type" : "infunc",
  "medium_grid_Vp_file" : "mod1_Vp.grd",
  "medium_grid_Vs_file" : "mod1_Vs.grd",
  "medium_grid_rho_file" : "mod1_rho.grd",
  "medium_grid_Qp_file" : "mod1_Qp.grd",
  "medium_grid_Qs_file" : "mod1_Qs.grd",

  "input_source_file" : "test.src",

  "input_receiver_file" : "test.rcv",

  "receiver_line" : [
    {
      "name" : "line_x_1",
      "grid_index_start"    : [  0, 49, 59 ],
      "grid_index_incre"    : [  1,  0,  0 ],
      "count"    : 100
    },
    {
      "name" : "line_y_1",
      "grid_index_start"    : [ 19, 49, 59 ],
      "grid_index_incre"    : [  0,  1,  0 ],
      "count"    : 20
    } 
  ],

  "slice" : {
      "x_index" : [ 19, 49, 59 ],
      "y_index" : [ 20, 50, 60 ],
      "z_index" : [ 29, 59 ]
  },

  "snapshot" : [
    {
      "name" : "volume_vel",
      "grid_index_start" : [ 10, 0, 0 ],
      "grid_index_count" : [ 40,50, 50 ],
      "grid_index_incre" : [  2, 2, 1 ],
      "time_index_start" : 0,
      "time_index_incre" : 1,
      "save_velocity" : 1,
      "save_stress"   : 0,
      "save_strain"   : 0
    }
  ],

  "check_nan_every_nummber_of_steps" : 0,
  "output_all" : 0,

  "grid_dir"   : "$grid_dir",
  "media_dir"  : "$media_dir",
  "output_dir" : "$output_dir"
}
ieof

echo "+ created $par_file"

#----------------------------------------------------------------------
#-- create source list
#----------------------------------------------------------------------
#cat << ieof > ${PROJDIR}/source.list
#ieof
#echo "+ created source.list"

#
#-------------------------------------------------------------------------------
#-- Performce simulation
#-------------------------------------------------------------------------------
#

#-- get np
NUMPROCS_X=`grep number_of_mpiprocs_x ${par_file} | sed 's/:/ /g' | sed 's/,/ /g' | awk '{print $2}'`
NUMPROCS_Y=`grep number_of_mpiprocs_y ${par_file} | sed 's/:/ /g' | sed 's/,/ /g' | awk '{print $2}'`
NUMPROCS=$(( NUMPROCS_X*NUMPROCS_Y ))
echo $NUMPROCS_X $NUMPROCS_Y $NUMPROCS

cat << ieof > ${PROJDIR}/cgfd_sim.sh
#!/bin/bash

set -e
printf "\nUse $NUMPROCS CPUs on following nodes:\n"
cat ${PROJDIR}/hostlist

#-- simulation
printf "\nStart simualtion ...\n";
#time $MPIDIR/bin/mpiexec -f ${PROJDIR}/hostlist -np $NUMPROCS $EXEC_DIR/fd3dpmltopo_ew $MAINCONF
time $MPIDIR/bin/mpiexec -machinefile ${PROJDIR}/hostlist -np $NUMPROCS $EXEC_WAVE $par_file 100
if [ $? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi

ieof

#-------------------------------------------------------------------------------
#-- start run
#-------------------------------------------------------------------------------

chmod 755 ${PROJDIR}/cgfd_sim.sh
${PROJDIR}/cgfd_sim.sh
if [ $? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi

date

#
#-------------------------------------------------------------------------------
#-- plot results
#-------------------------------------------------------------------------------
#

#EXEC_Plot_Script=./hill3d_plot.sh

#if [[ -n "${Plot_Result}" && ${Plot_Result} -eq 1 ]]; then
#   printf "\nPlot result using %s\n", ${EXEC_Plot_Script}
#   time ${EXEC_Plot_Script}
#else
#   printf "\nDo not plot result\n"
#fi
#
#date

# vim:ft=conf:ts=4:sw=4:nu:et:ai:
