#!/bin/bash

#set -x
set -e

date

#-- system related dir
MPIDIR=/data3/lihl/software/openmpi-gnu-4.1.2

#-- program related dir
EXEC_WAVE=`pwd`/../main_curv_col_el_3d
echo "EXEC_WAVE=$EXEC_WAVE"

#-- input dir
INPUTDIR=`pwd`

#-- output and conf
PROJDIR=~/work/cgfd3d-wave-el/10src
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
#-- create in soruce file
#----------------------------------------------------------------------
cat << ieof > $INPUTDIR/test.src
# name of this input source
event_1
# number of source
1
# flag for stf
#  1st value : 0 analytic stf or 1 discrete values
#    for analytical, 2nd value is time length
#    for discrete,  2nd value is dt and 3rd is nt
#    e.g.,
# 0 4.0 
# 1 0.05 20
0 1.0
# flag for source component and mechanism format
#  1st value: source components, 1(force), 2(momoment), 3(force+moment)
#  2nd value: mechanism format for moment source: 0 moment, 1 angle + mu + D + A
2 0
# flag for location
#   1st value: 0 computational coordiate, 1 physical coordinate
#   2nd value: third coordinate is 0 axis or 1 depth 
0 0
# location of each source
#   sx sy sz
80.2 49.3 50.5
#80 49 50
#8020 4930 -950
# stf and cmp
0.0 ricker 2.0 0.5   # t0  stf_name  ricker_fc ricker_t0
1e16  1e16  1e16 0 0 0 
ieof

#cp -f $INPUTDIR/test_by_values.src $INPUTDIR/test.src

#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > $PAR_FILE
{
  "number_of_total_grid_points_x" : 300,
  "number_of_total_grid_points_y" : 300,
  "number_of_total_grid_points_z" : 60,

  "number_of_mpiprocs_x" : 2,
  "number_of_mpiprocs_y" : 2,

  "#size_of_time_step" : 0.008,
  "#size_of_time_step" : 0.020,
  "#number_of_time_steps" : 500,
  "time_window_length" : 4,
  "check_stability" : 1,

  "boundary_x_left" : {
      "cfspml" : {
          "number_of_layers" : 10,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 7000.0
          }
      },
  "boundary_x_right" : {
      "cfspml" : {
          "number_of_layers" : 10,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 7000.0
          }
      },
  "boundary_y_front" : {
      "cfspml" : {
          "number_of_layers" : 10,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 7000.0
          }
      },
  "boundary_y_back" : {
      "cfspml" : {
          "number_of_layers" : 10,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 7000.0
          }
      },
  "boundary_z_bottom" : {
      "cfspml" : {
          "number_of_layers" : 10,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 7000.0
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
        "in_grid_layer_file" : "$INPUTDIR/prep_grid/tangshan_area_topo.gdlay",
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
      "code" : "func_name_here",
      "#import" : "$MEDIA_DIR",
      "#infile_layer" : "$INPUTDIR/prep_medium/basin_el_iso.md3lay",
      "#infile_grid" : "$INPUTDIR/prep_medium/topolay_el_iso.md3grd",
      "equivalent_medium_method" : "loc",
      "#equivalent_medium_method" : "har"
  },
  "is_export_media" : 1,
  "media_export_dir"  : "$MEDIA_DIR",

  "#visco_config" : {
      "type" : "graves_Qs",
      "Qs_freq" : 1.0
  },

  "in_source_file" : "$INPUTDIR/test.src",
  "is_export_source" : 1,
  "source_export_dir"  : "$SOURCE_DIR",

  "output_dir" : "$OUTPUT_DIR",

  "in_station_file" : "$INPUTDIR/station.list",

  "#receiver_line" : [
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

  "#slice" : {
      "x_index" : [ 190 ],
      "y_index" : [ 90 ],
      "z_index" : [ 59 ]
  },

  "snapshot" : [
    {
      "name" : "volume_vel",
      "grid_index_start" : [ 0, 0, 59 ],
      "grid_index_count" : [ 300,300, 1 ],
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

printf "\nStart simualtion ...\n";
time $MPIDIR/bin/mpiexec -machinefile ${PROJDIR}/hostlist -np $NUMPROCS $EXEC_WAVE $PAR_FILE 11
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
