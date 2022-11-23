#!/bin/bash

#set -x

date

#-- system related dir, from module env or manually set
MPIDIR=/share/apps/gnu-4.8.5/mpich-3.3
#MPIDIR=$MPI_ROOT

#-- program related dir
CUR_DIR=`pwd`

EXE_DIR=`dirname ${CUR_DIR}`
EXEC_WAVE=${EXE_DIR}/main_curv_col_el_3d
echo "EXEC_WAVE=$EXEC_WAVE"

#-- output and conf
PROJDIR=~/work/cgfd3d-wave-el/cfspml
EVTNM=codetest
echo "PROJDIR=${PROJDIR}"
echo "EVTNM=${EVTNM}"

PAR_FILE=${PROJDIR}/test.json
GRID_DIR=${PROJDIR}/output
MEDIA_DIR=${PROJDIR}/output
SOURCE_DIR=${PROJDIR}/output
OUTPUT_DIR=${PROJDIR}/output

#RUN_SCRIPT_FILE=${PROJDIR}/runscript.sh
RUN_SCRIPT_FILE=${PROJDIR}/runscript.lsf
echo "RUN_SCRIPT_FILE=${RUN_SCRIPT_FILE}"

#-- create dir
mkdir -p $PROJDIR
mkdir -p $OUTPUT_DIR
mkdir -p $GRID_DIR
mkdir -p $MEDIA_DIR

#-- tmp dir should be in local, here is for test
TMP_DIR=${PROJDIR}/scratch
mkdir -p ${TMP_DIR}

#----------------------------------------------------------------------
#-- grid and mpi configurations
#----------------------------------------------------------------------

#-- total x grid points
NX=300
#-- total y grid points
NY=300
#-- total z grid points
NZ=60
#-- total x mpi procs
NPROCS_X=2
#-- total y mpi procs
NPROCS_Y=2
#-- total mpi procs
NPROCS=$(( NPROCS_X*NPROCS_Y ))

#----------------------------------------------------------------------
#-- create in soruce file
#----------------------------------------------------------------------

create_source_file()
{
cat << ieof > ${PROJDIR}/test.src
# name of this input source
${EVTNM}
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
#  2nd value: mechanism format for moment source:
#       0 : 6 moment components, 
#       1 : 3 angles, mu, slip rate or D, A (mu <0 means to use internal mu value)
2 0
# flag for location
#   1st value: 0 computational coordiate, 1 physical coordinate
#   2nd value: third coordinate is 0 axis or 1 depth 
0 1
# location of each source
#   sx sy sz
80.2 49.3 5.5
#80 49 50
#8020 4930 -950
# stf and cmp
0.0 ricker 2.0 0.5   # t0  stf_name  ricker_fc ricker_t0
0  0  0 0 0 0 
ieof

echo "+ created $PROJDIR/test.src"
}

#----------------------------------------------------------------------
#-- create station list file
#----------------------------------------------------------------------

create_station_file()
{
cat << ieof > $PROJDIR/test.station
# number of station
1
# name is_grid_indx is_3dim_depth  x y z
r1  0  1  1000 1000 0
ieof

echo "+ created $PROJDIR/test.station"
}

#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------

create_json_file()
{
cat << ieof > $PAR_FILE
{
  "number_of_total_grid_points_x" : $NX,
  "number_of_total_grid_points_y" : $NY,
  "number_of_total_grid_points_z" : $NZ,

  "number_of_mpiprocs_x" : $NPROCS_X,
  "number_of_mpiprocs_y" : $NPROCS_Y,

  "#size_of_time_step" : 0.008,
  "#size_of_time_step" : 0.020,
  "#number_of_time_steps" : 500,
  "time_window_length" : 6,
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
        "in_grid_layer_file" : "${CUR_DIR}/prep_grid/tangshan_area_topo.gdlay",
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
      "input_way" : "code",
      "#input_way" : "binfile",
      "binfile" : {
        "size"    : [1101, 1447, 1252],
        "spacing" : [-10, 10, 10],
        "origin"  : [0.0,0.0,0.0],
        "dim1" : "z",
        "dim2" : "x",
        "dim3" : "y",
        "Vp" : "${CUR_DIR}/prep_medium/seam_Vp.bin",
        "Vs" : "${CUR_DIR}/prep_medium/seam_Vs.bin",
        "rho" : "${CUR_DIR}/prep_medium/seam_rho.bin"
      },
      "code" : "func_name_here",
      "#import" : "$MEDIA_DIR",
      "#infile_layer" : "${CUR_DIR}/prep_medium/basin_el_iso.md3lay",
      "#infile_grid" : "${CUR_DIR}/prep_medium/topolay_el_iso.md3grd",
      "equivalent_medium_method" : "loc",
      "#equivalent_medium_method" : "har"
  },
  "is_export_media" : 1,
  "media_export_dir"  : "$MEDIA_DIR",

  "#visco_config" : {
      "type" : "graves_Qs",
      "Qs_freq" : 1.0
  },

  "in_source_file" : "$PROJDIR/test.src",
  "is_export_source" : 1,
  "source_export_dir"  : "$SOURCE_DIR",

  "in_ddsource_file" : "${CUR_DIR}/prep_source/event_3moment_srcdd.nc",

  "output_dir" : "$OUTPUT_DIR",
  "tmp_dir"    : "$TMP_DIR",

  "in_station_file" : "$PROJDIR/test.station",

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

  "slice" : {
      "x_index" : [ 90, 149 ],
      "y_index" : [ 90, 149 ],
      "z_index" : [ 30, 59 ]
  },

  "snapshot" : [
    {
      "name" : "volume_vel",
      "grid_index_start" : [ 0, 0, $(( NZ - 1 )) ],
      "grid_index_count" : [ $NX, $NY, 1 ],
      "grid_index_incre" : [  1, 1, 1 ],
      "time_index_start" : 0,
      "time_index_incre" : 1,
      "save_velocity" : 1,
      "save_stress"   : 1,
      "save_strain"   : 0
    }
  ],

  "check_nan_every_nummber_of_steps" : 0,
  "output_all" : 0 
}
ieof

echo "+ created $PAR_FILE"
}

#-------------------------------------------------------------------------------
#-- generate run script
#-------------------------------------------------------------------------------

#-- for lsf
create_script_lsf()
{
cat << ieof > ${RUN_SCRIPT_FILE}
#!/bin/bash

#-- Job Name
#BSUB -J ${EVTNM} 
#-- Queue name
#--  mars: normal, large;  earth: short, large, long
#BSUB -q normal
#-- requires number of cores (default: 1)
#BSUB -n ${NPROCS}
#-- Other specification
##BSUB -R "hname!=c013"
#-- Merge stderr with stdout, %J is the job-id
#BSUB -o stdout.%J

printf "\nUse $NPROCS CPUs on following nodes:\n"

MPI_CMD="$MPIDIR/bin/mpiexec -np ${NPROCS} ${EXEC_WAVE} ${PAR_FILE} 110"

printf "\nStart simualtion ...\n";
printf "%s\n\n" "\${MPI_CMD}";

time \${MPI_CMD};
if [ \$? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi
ieof
}

#-- for pbs
create_script_pbs()
{
cat << ieof > ${RUN_SCRIPT_FILE}
#!/bin/bash

#-- Job Name
#PBS -N ${EVTNM} 
#-- Queue name, should check the naming of the system
#PBS -q normal
#-- requires number of cores
#--  be aware that pbs requires to specify number of nodes, ncpus per node
#PBS -l select=1:ncpus=${NPROCS}
#-- Merge stderr with stdout
#PBS -j oe

printf "\nUse $NPROCS CPUs on following nodes:\n"
printf "%s " \`cat ${PBS_NODEFILE} | sort\`;

MPI_CMD="$MPIDIR/bin/mpiexec -np ${NPROCS} ${EXEC_WAVE} ${PAR_FILE} 110"

printf "\nStart simualtion ...\n";
printf "%s\n\n" "\${MPI_CMD}";

time \${MPI_CMD};
if [ \$? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi
ieof
}

#-- for directly run
create_script_run()
{
#-- create hostlist
cat << ieof > ${PROJDIR}/hostlist
server1
ieof

#-- create runscript
cat << ieof > ${RUN_SCRIPT_FILE}
#!/bin/bash

printf "\nUse $NPROCS CPUs on following nodes:\n"
printf "%s " \`cat ${PROJDIR}/hostlist | sort\`;

MPI_CMD="$MPIDIR/bin/mpiexec -machinefile ${PROJDIR}/hostlist -np $NPROCS $EXEC_WAVE $PAR_FILE 1000"

printf "\nStart simualtion ...\n";
printf "%s\n\n" "'\${MPI_CMD}'";

time \${MPI_CMD};
if [ \$? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi
ieof

chmod 755 ${RUN_SCRIPT_FILE}
}

#-------------------------------------------------------------------------------
# submit or run
#-------------------------------------------------------------------------------

#-- common steps
create_json_file;

create_source_file;

create_station_file;

#-- run with lsf
#echo "sumbit to lsf ..."
#create_script_lsf;
#bsub < ${RUN_SCRIPT_FILE}

#-- run with pbs
#echo "sumbit to pbs ..."
#create_script_pbs;
#qsub ${RUN_SCRIPT_FILE}

#-- directly run
echo "start run script ..."
create_script_run;
${RUN_SCRIPT_FILE}

date

# vim:ts=4:sw=4:nu:et:ai:
