#!/bin/bash

#set -x
set -e

date

#-- source intel lib
#source /share/apps/intel-2019.3/intel/bin/compilervars.sh intel64
#export FI_PROVIDER=tcp
#MPIDIR=/share/apps/intel-2019.3/mpich-3.3

#-- system related dir
MPIDIR=/share/apps/gnu-4.8.5/mpich-3.3

#-- program related dir
EXEC_DIR=/home/zhangw/code/zwlab/CGFD3D-elastic
EXEC_WAVE=$EXEC_DIR/cgfdm3d_elastic_mpi

#-- conf
PROJDIR=/export/home/zhangw/work/cgfd_opt/14_jsonnew
PAR_FILE=${PROJDIR}/test.json
GRID_DIR=${PROJDIR}/output
MEDIA_DIR=${PROJDIR}/output
SOURCE_DIR=${PROJDIR}/output
OUTPUT_DIR=${PROJDIR}/output

#-- create dir
mkdir -p $PROJDIR
mkdir -p $OUTPUT_DIR
mkdir -p $GRID_DIR
mkdir -p $MEDIA_DIR

#----------------------------------------------------------------------
#-- create hostlist for mpirun
#----------------------------------------------------------------------
cat << ieof > ${PROJDIR}/hostlist
server1
ieof

#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > $PAR_FILE
{
  "number_of_total_grid_points_x" : 100,
  "number_of_total_grid_points_y" : 100,
  "number_of_total_grid_points_z" : 60,

  "number_of_mpiprocs_x" : 1,
  "number_of_mpiprocs_y" : 1,

  "size_of_time_step" : 0.01,
  "number_of_time_steps" : 500,

  "boundary_x_left" : {
      "cfspml" : {
          "number_of_layers" : 5,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 3000.0
          }
      },
  "boundary_x_right" : {
      "cfspml" : {
          "number_of_layers" : 5,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 3000.0
          }
      },
  "boundary_y_front" : {
      "cfspml" : {
          "number_of_layers" : 5,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 3000.0
          }
      },
  "boundary_y_back" : {
      "cfspml" : {
          "number_of_layers" : 5,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 3000.0
          }
      },
  "boundary_z_bottom" : {
      "cfspml" : {
          "number_of_layers" : 5,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 3000.0
          }
      },
  "boundary_z_top" : {
      "free" : "timg"
      },

  "grid_generation_method" : {
      "#import" : "$GRID_DIR",
      "cartesian" : {
        "origin"  : [0.0, 0.0, -6000.0 ],
        "inteval" : [ 100.0, 100.0, 100.0 ]
      },
      "#layer_interp" : {
        "in_grid_layer_file" : "$PROJDIR/test/grid_layer.gdlay",
        "refine_factor" : [ 2, 2, 2 ]
      }
  },
  "is_export_grid" : 1,
  "grid_export_dir"   : "$GRID_DIR",

  "metric_calculation_method" : {
      "#import" : "$GRID_DIR",
      "calculate" : 1
  },
  "is_export_metric" : 1,

  "media_input" : {
      "#import" : "$MEDIA_DIR",
      "code_generate" : 1,
      "#in_3lay_file" : "$PROJDIR/test/hill3d.md3lay",
      "#in_3grd_file" : "$PROJDIR/test/hill3d.md3grd"
  },
  "is_export_media" : 1,
  "media_export_dir"  : "$MEDIA_DIR",

  "source_input" : {
      "code_generate" : 1,
      "in_source_file" : "$PROJDIR/test/forceandmoment.fdsrc"
  },
  "is_export_source" : 1,
  "source_export_dir"  : "$SOURCE_DIR",

  "output_dir" : "$OUTPUT_DIR",

  "receiver_file" : "surface.fdrcv",

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
      "grid_index_start" : [ 10, 1, 0 ],
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
  "output_all" : 0 
}
ieof

echo "+ created $PAR_FILE"

#----------------------------------------------------------------------
#-- create source list
#----------------------------------------------------------------------
#cat << ieof > ${PROJDIR}/source.list
#ieof
#echo "+ created source.list"

#-------------------------------------------------------------------------------
#-- Performce simulation
#-------------------------------------------------------------------------------
#

#-- get np
NUMPROCS_X=`grep number_of_mpiprocs_x ${PAR_FILE} | sed 's/:/ /g' | sed 's/,/ /g' | awk '{print $2}'`
NUMPROCS_Y=`grep number_of_mpiprocs_y ${PAR_FILE} | sed 's/:/ /g' | sed 's/,/ /g' | awk '{print $2}'`
NUMPROCS=$(( NUMPROCS_X*NUMPROCS_Y ))
echo $NUMPROCS_X $NUMPROCS_Y $NUMPROCS

#-- gen run script
cat << ieof > ${PROJDIR}/cgfd_sim.sh
#!/bin/bash

set -e
printf "\nUse $NUMPROCS CPUs on following nodes:\n"
cat ${PROJDIR}/hostlist

#-- simulation
printf "\nStart simualtion ...\n";
#time $MPIDIR/bin/mpiexec -f ${PROJDIR}/hostlist -np $NUMPROCS $EXEC_DIR/fd3dpmltopo_ew $MAINCONF
time $MPIDIR/bin/mpiexec -machinefile ${PROJDIR}/hostlist -np $NUMPROCS $EXEC_WAVE $PAR_FILE 100
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
