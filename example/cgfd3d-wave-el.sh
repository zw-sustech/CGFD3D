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
#EXEC_WAVE=/home/zhangw/code/zwlab/CGFD3D-elastic/main_curv_col_el_3d
EXEC_WAVE=`pwd`/../main_curv_col_el_3d
echo "EXEC_WAVE=$EXEC_WAVE"

#-- input dir
INPUTDIR=`pwd`

#-- output and conf
PROJDIR=~/work/cgfd3d-wave-el/05grd
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
  "number_of_total_grid_points_x" : 120,
  "number_of_total_grid_points_y" : 100,
  "number_of_total_grid_points_z" : 60,

  "number_of_mpiprocs_x" : 1,
  "number_of_mpiprocs_y" : 1,

  "#size_of_time_step" : 0.008,
  "#size_of_time_step" : 0.015,
  "#number_of_time_steps" : 500,
  "time_window_length" : 4,
  "check_stability" : 1,

  "boundary_x_left" : {
      "cfspml" : {
          "number_of_layers" : 5,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 5000.0
          }
      },
  "boundary_x_right" : {
      "cfspml" : {
          "number_of_layers" : 5,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 5000.0
          }
      },
  "boundary_y_front" : {
      "cfspml" : {
          "number_of_layers" : 5,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 5000.0
          }
      },
  "boundary_y_back" : {
      "cfspml" : {
          "number_of_layers" : 5,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 5000.0
          }
      },
  "boundary_z_bottom" : {
      "cfspml" : {
          "number_of_layers" : 5,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 5000.0
          }
      },
  "boundary_z_top" : {
      "free" : "timg"
      },

  "grid_generation_method" : {
      "#import" : "$GRID_DIR",
      "cartesian" : {
        "origin"  : [0.0, 0.0, -5900.0 ],
        "inteval" : [ 100.0, 100.0, 100.0 ]
      },
      "#layer_interp" : {
        "in_grid_layer_file" : "$INPUTDIR/prep_grid/random_topo.gdlay",
        "refine_factor" : [ 1, 1, 1 ],
        "horizontal_start_index" : [ 3, 3 ],
        "vertical_last_to_top" : 0
      }
  },
  "is_export_grid" : 1,
  "grid_export_dir"   : "$GRID_DIR",

  "metric_calculation_method" : {
      "#import" : "$GRID_DIR",
      "calculate" : 1
  },
  "is_export_metric" : 1,

  "medium" : {
      "type" : "elastic_iso",
      "#type" : "elastic_vti",
      "#code" : "func_name_here",
      "#import" : "$MEDIA_DIR",
      "#infile_layer" : "$INPUTDIR/prep_medium/basin_el_iso.md3lay",
      "infile_grid" : "$INPUTDIR/prep_medium/topolay_el_iso.md3grd",
      "equivalent_medium_method" : "har"
  },
  "is_export_media" : 1,
  "media_export_dir"  : "$MEDIA_DIR",

  "#visco_config" : {
      "type" : "graves_Qs",
      "Qs_freq" : 1.0
  },

  "source_input" : {
      "in_par" : {
         "name" : "evt_by_par",
         "source" : [
            {
                "index" : [ 80, 49, 50 ],
                "#coord" : [ 4000, 4000, -1000 ],
                "wavelet_name" : "ricker",
                "ricker_center_frequency" : 2.0,
                "ricker_peak_time" : 0.5,
                "#wavelet_name" : "gaussian",
                "#gaussian_rms_width" : 2.0,
                "#gaussian_peak_time" : 0.5,
                "start_time" : 0.0,
                "end_time"   : 1.0,
                "#force_vector" : [ 0.0, 0.0, 1e16],
                "moment_tensor" : [ 1e16, 1e16, 1e16, 0.0, 0.0, 0.0]
            }
         ]
      },
      "#in_source_file" : "$INPUTDIR/source.anasrc"
  },
  "is_export_source" : 1,
  "source_export_dir"  : "$SOURCE_DIR",

  "output_dir" : "$OUTPUT_DIR",

  "in_station_file" : "$INPUTDIR/station.list",

  "receiver_line" : [
    {
      "name" : "line_x_1",
      "grid_index_start"    : [  0, 49, 59 ],
      "grid_index_incre"    : [  1,  0,  0 ],
      "grid_index_count"    : 20
    },
    {
      "name" : "line_y_1",
      "grid_index_start"    : [ 19, 49, 59 ],
      "grid_index_incre"    : [  0,  1,  0 ],
      "grid_index_count"    : 20
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
      "grid_index_start" : [ 0, 0, 0 ],
      "grid_index_count" : [ 120,100, 60 ],
      "grid_index_incre" : [  1, 1, 1 ],
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

printf "\nStart simualtion ...\n";
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

# vim:ft=conf:ts=4:sw=4:nu:et:ai: