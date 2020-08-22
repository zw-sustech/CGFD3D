#!/bin/bash

#set -x
set -e

date

#-- source intel lib
#source /share/apps/intel-2019.3/intel/bin/compilervars.sh intel64
#export FI_PROVIDER=tcp
#MPIDIR=/share/apps/intel-2019.3/mpich-3.3

#-- sever related dir
MPIDIR=/share/apps/gnu-4.8.5/mpich-3.3
EXEC_MATLAB=/share/apps/Matlab/R2015b/bin/matlab2015

#-- user related dir
EXEC_DIR=/home/zhangw/code/zwlab/CGFD3D-elastic
EXEC_WAVE=$EXEC_DIR/cgfdm3d_elastic_mpi

#-- project related dir
PROJDIR=/export/home/zhangw/work/cgfd_init_run
export PROJDIR

output_dir=$PROJDIR/output
mkdir -p $output_dir

EXEC_Plot_Script=./hill3d_plot.sh

#-- conf file names
parfile=${PROJDIR}/test.json

#-- create prj dir
mkdir -p $PROJDIR

#--
#Run_Sim_Only=1
Plot_Result=1

if [[ -z "${Run_Sim_Only}" || ${Run_Sim_Only} -eq 0 ]]; then

#----------------------------------------------------------------------
#-- create hostlist for mpirun
#----------------------------------------------------------------------
cat << ieof > ${PROJDIR}/hostlist
server1
ieof

#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > $parfile
{
  "grid_name" : "blk_1",
  "number_of_total_grid_points_x" : 100,
  "number_of_total_grid_points_y" : 100,
  "number_of_total_grid_points_z" : 60,

  "number_of_mpiprocs_x" : 1,
  "number_of_mpiprocs_y" : 1,

  "size_of_time_step" : 0.01,
  "number_of_time_steps" : 1500,

  "coord_by_cartesian" : 1,
  "cartesian_grid_x0" : 0.0,
  "cartesian_grid_y0" : 0.0,
  "cartesian_grid_z0" : -6000.0,
  "cartesian_grid_dx" : 100.0,
  "cartesian_grid_dy" : 100.0,
  "cartesian_grid_dz" : 100.0,

  "metric_by_import" : 0,

  "medium_by_import" : 0,

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
  "#snapshot" : [
    {
      "name" : "snap_sz",
      "grid_index_start"             : [  1, 1, 50 ],
      "grid_index_count"             : [ 100,100, 1 ],
      "grid_index_stride"            : [  1, 1, 1 ],
      "time_step_start_count_stride" : [  0, 2000, 1 ]
    },
    {
      "name" : "snap_sy",
      "grid_index_start"             : [  1, 50, 1 ],
      "grid_index_count"             : [ 100, 1, 60 ],
      "grid_index_stride"            : [  1, 1, 1 ],
      "time_step_start_count_stride" : [  0, 2000, 1 ]
    },
    {
      "name" : "snap_sx",
      "grid_index_start"             : [ 50, 1, 1 ],
      "grid_index_count"             : [  1,100,60 ],
      "grid_index_stride"            : [  1, 1, 1 ],
      "time_step_start_count_stride" : [  0, 2000, 1 ]
    }
  ],

  "output_dir" : "$output_dir"
}
ieof

echo "+ created $parfile"

#----------------------------------------------------------------------
#-- create source list
#----------------------------------------------------------------------
cat << ieof > ${PROJDIR}/source.list
# input format:
# number_of_subevent
#     e.g.:  1
# event_type evt_axis_is_grid sx1 sy1 sz1 if_relocate_evt_to_surface
#   # event_type: 1:force,2:moment_in_tensor,3:moment_in_angle
#   # evt_axis_is_grid: 0:axis, 1:grid
#   # if_relocate_evt_to_surface: 0 false, 1 true
#     e.g.:  1 1 0.0 -800.0 9999.0 0
# if_force:
#  fx1 fy1 fz1 |f|
#     e.g.: 0 0 1 1e16
# if moment_in_tensor
#  mxx myy mzz mxy mxz myz m0 m0_need_multiple_mu
#     e.g.: 1.0 1.0 1.0 0.0 0.0 0.0 1e16 1
# if_moment_in_angle:
#  strike dip rake m0 m0_need_multiple_mu
#     e.g.: 45 90 90 1e16 0
# stf_name/mrf_name(ricker,sampled)
#     e.g.: ricker
# if_ricker_etc:
#   ricker_fc ricker_t0
#     e.g.: 1.5 1.0
#   other stf
#    gauss: width, center_t0
#    bell: width, start_t0
#    triangle: width, start_t0
#    step: start_t0
#    delta: delta_t
# if_sample:
#   nt dt
#   stf_1
#   stf_2
#     e.g.: 1000 0.01
#           0.0
#           0.0
1
# evt1
#1 1 0.0 -800.0 9999.0 1
#0.0 0.0 1.0 1e16
#ricker
#1.5 1.0
##- evt2
2 1 -800 0.0 427.3 0
1.0 1.0 1.0 0.0 0.0 0.0 1.0e+16 0
ricker_deriv
1.5 1.0
#gauss
#0.223 0.75
#- ev for interp
#1 1 -800 0.0 427.3 0
#0.0 0.0 1.0 1e16
#sampled
#151        0.0010000
#       0.0000000
#       0.0012798
#       0.0031477
#       0.0058203
#       0.0095688
#       0.0147228
#       0.0216696
#       0.0308483
#       0.0427371
#       0.0578327
#       0.0766224
#       0.0995493
#       0.1269731
#       0.1591297
#       0.1960926
#       0.2377433
#       0.2837514
#       0.3335714
#       0.3864563
#       0.4414885
#       0.4976267
#       0.5537650
#       0.6087972
#       0.6616821
#       0.7115021
#       0.7575102
#       0.7991609
#       0.8361238
#       0.8682803
#       0.8957042
#       0.9186311
#       0.9374208
#       0.9525164
#       0.9644052
#       0.9735839
#       0.9805307
#       0.9856847
#       0.9894332
#       0.9921058
#       0.9939737
#       0.9952535
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
#       0.9957702
ieof

echo "+ created source.list"

#
#-------------------------------------------------------------------------------
#-- create project directory structure
#-------------------------------------------------------------------------------
#

#CHECKPOINT_ROOT=`grep checkingpoint_dir ${MAINCONF} | awk '{print $3}'`
#GRID_ROOT=`grep grid_dir ${MAINCONF} | awk '{print $3}'`
#METRIC_ROOT=`grep metric_dir ${MAINCONF} | awk '{print $3}'`
#MEDIA_ROOT=`grep medium_dir ${MAINCONF} | awk '{print $3}'`
#SOURCE_ROOT=`grep source_dir ${MAINCONF} | awk '{print $3}'`
#STATION_ROOT=`grep station_dir ${MAINCONF} | awk '{print $3}'`
#OUTPUT_ROOT=`grep output_dir ${MAINCONF} | awk '{print $3}'`
#
#if [[ -z ${CHECKPOINT_ROOT} || \
#     -z ${GRID_ROOT} || \
#     -z ${METRIC_ROOT} || \
#     -z ${MEDIA_ROOT} || \
#     -z ${SOURCE_ROOT} || \
#     -z ${STATION_ROOT} || \
#     -z ${OUTPUT_ROOT} \
#     ]]; then
#    echo "!!!Error: one or more _ROOT empty"
#    echo " MAINCONF=${MAINCONF}"
#    exit
#fi
#
##set -x
#
#mkdir -p ${CHECKPOINT_ROOT}
#mkdir -p ${GRID_ROOT}
#mkdir -p ${METRIC_ROOT}
#mkdir -p ${MEDIA_ROOT}
#mkdir -p ${SOURCE_ROOT}
#mkdir -p ${STATION_ROOT}
#mkdir -p ${OUTPUT_ROOT}
#
#echo "0 0 0 # checkpoint, syncpoint, new nt" > ${PROJDIR}/checkpoint.dat
##set +x


##-------------------------------------------------------------------------------
##-- pre-processing
##-------------------------------------------------------------------------------
#
#cat << ieof > ${PROJDIR}/fd3dtopoew_prep.sh
##!/bin/bash
##
#set -e
##
##-- generate grid
#echo "+ generate grid ..."
#${EXEC_DIR}/seis3d_grid $MAINCONF
#
#if [ $? -ne 0 ]; then
#    printf "\ngrid generation fail! stop!\n"
#    exit 1
#fi
#
##-- generate metric
#echo "+ calculate grid metric..."
#${EXEC_DIR}/seis3d_metric $MAINCONF
#
#if [ $? -ne 0 ]; then
#    printf "\ncoordinate mapping calculation fail! stop!\n"
#    exit 1
#fi
#
##-- generate media
#echo "+ evaluate media ..."
#${EXEC_DIR}/seis3d_media $MAINCONF
#if [ $? -ne 0 ]; then
#    printf "\ndiscretize media fail! stop!\n"
#    exit 1
#fi
#
##-- generate source
#echo "+ discrete source ..."
#${EXEC_DIR}/seis3d_source $MAINCONF
#if [ $? -ne 0 ]; then
#    printf "\nlocate source fail! stop!\n"
#    exit 1
#fi
#
##-- generate station
#echo "+ discrete station ..."
#${EXEC_DIR}/seis3d_station $MAINCONF
#if [ $? -ne 0 ]; then
#    printf "\nlocate station fail! stop!\n"
#    exit 1
#fi
#
#ieof
#
##-- start prep run
#chmod 755 ${PROJDIR}/fd3dtopoew_prep.sh
#${PROJDIR}/fd3dtopoew_prep.sh
#
#if [ $? -ne 0 ]; then
#    printf "\nlocate source fail! stop!\n"
#    exit 1
#fi
#
#date

else
   echo "\nSkip preprocessing, only perform simulation"
fi # Run_Sim_Only

#
#-------------------------------------------------------------------------------
#-- Performce simulation
#-------------------------------------------------------------------------------
#

#-- get np
#NUMPROCS=`grep ^num_threads_per_dim ${MAINCONF} | sed 's/=/ /g' | awk '{print $2 * $3 * $4}'`
NUMPROCS=1

cat << ieof > ${PROJDIR}/cgfd_sim.sh
#!/bin/bash

set -e
printf "\nUse $NUMPROCS CPUs on following nodes:\n"
cat ${PROJDIR}/hostlist

#-- simulation
printf "\nStart simualtion ...\n";
#time $MPIDIR/bin/mpiexec -f ${PROJDIR}/hostlist -np $NUMPROCS $EXEC_DIR/fd3dpmltopo_ew $MAINCONF
time $MPIDIR/bin/mpiexec -machinefile ${PROJDIR}/hostlist -np $NUMPROCS $EXEC_WAVE $parfile 100
if [ $? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi

ieof

#-- start sim run
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

#if [[ -n "${Plot_Result}" && ${Plot_Result} -eq 1 ]]; then
#   printf "\nPlot result using %s\n", ${EXEC_Plot_Script}
#   time ${EXEC_Plot_Script}
#else
#   printf "\nDo not plot result\n"
#fi
#
#date

# vim:ft=conf:ts=4:sw=4:nu:et:ai:
