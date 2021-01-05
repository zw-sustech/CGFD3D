#!/bin/bas

#set -x
set -e

date

#-- source intel lib
#source /share/apps/intel-2019.3/intel/bin/compilervars.sh intel64
#export FI_PROVIDER=tcp
#MPIDIR=/share/apps/intel-2019.3/mpich-3.3

#-- sever related dir
MPIDIR=/export/apps/gnu-4.8.5/mpich-3.3
#EXEC_MATLAB=/export/apps/Matlab/R2015b/bin/matlab2015

#-- user related dir
EXEC_DIR=/export/home/wangyh/CGFD3D-elastic
EXEC_WAVE=$EXEC_DIR/cgfdm3d_elastic_mpi

#-- project related dir
proj_dir=/export/home/wangyh/CGFD3D-elastic/project

export proj_dir

output_dir=$proj_dir/output
mkdir -p $output_dir

#EXEC_Plot_Script=./hill3d_plot.sh

#-- conf file names
par_file=${proj_dir}/test.json

#-- create prj dir
mkdir -p $proj_dir

#----------------------------------------------------------------------
#-- create hostlist for mpirun
#----------------------------------------------------------------------
cat << ieof > ${proj_dir}/hostlist
server1
ieof

#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > $par_file
{
  "number_of_total_grid_points_x" : 100,
  "number_of_total_grid_points_y" : 100,
  "number_of_total_grid_points_z" : 60,

  "number_of_mpiprocs_x" : 3,
  "number_of_mpiprocs_y" : 3,

  "size_of_time_step" : 0.01,
  "number_of_time_steps" : 500,

  "boundary_condition" : [
      "cfspml",
      "cfspml",
      "cfspml",
      "cfspml",
      "cfspml",
      "free"
  ],

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
      "grid_index_start" : [  0, 0,  0 ],
      "grid_index_count" : [ 100,100, 50 ],
      "grid_index_incre" : [  1, 1, 1 ],
      "time_index_start" : 0,
      "time_index_incre" : 1,
      "save_velocity" : 1,
      "save_stress"   : 0,
      "save_strain"   : 0
    }
  ],

  "check_nan_every_nummber_of_steps" : 0,
  "output_all" : 0,

  "output_dir" : "$output_dir"
}
ieof

echo "+ created $par_file"

#----------------------------------------------------------------------
#-- create source list
#----------------------------------------------------------------------
#cat << ieof > ${proj_dir}/source.list
#ieof
#echo "+ created source.list"

#-------------------------------------------------------------------------------
#-- Performce simulation
#-------------------------------------------------------------------------------
#

#-- get np
NUMPROCS_X=`grep number_of_mpiprocs_x ${par_file} | sed 's/:/ /g' | sed 's/,/ /g' | awk '{print $2}'`
NUMPROCS_Y=`grep number_of_mpiprocs_y ${par_file} | sed 's/:/ /g' | sed 's/,/ /g' | awk '{print $2}'`
NUMPROCS=$(( NUMPROCS_X*NUMPROCS_Y ))
echo $NUMPROCS_X $NUMPROCS_Y $NUMPROCS

#-- gen run script
cat << ieof > ${proj_dir}/cgfd_sim.sh
#!/bin/bash

set -e
printf "\nUse $NUMPROCS CPUs on following nodes:\n"
cat ${proj_dir}/hostlist

#-- simulation
printf "\nStart simualtion ...\n";
#time $MPIDIR/bin/mpiexec -f ${proj_dir}/hostlist -np $NUMPROCS $EXEC_DIR/fd3dpmltopo_ew $MAINCONF
time $MPIDIR/bin/mpiexec -machinefile ${proj_dir}/hostlist -np $NUMPROCS $EXEC_WAVE $par_file 100
if [ $? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi

ieof

#-- start run
chmod 755 ${proj_dir}/cgfd_sim.sh
${proj_dir}/cgfd_sim.sh
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

#EXEC_MATLAB=/share/apps/Matlab/R2015b/bin/matlab2015
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

