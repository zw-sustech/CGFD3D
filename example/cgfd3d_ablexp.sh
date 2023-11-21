#!/bin/bash

#set -x

date

#-- system related dir, from module env or manually set
#MPIDIR=/share/apps/gnu-4.8.5/mpich-3.3
MPIDIR=$MPI_ROOT

#-- program related dir
CUR_DIR=`pwd`

EXE_DIR=`dirname ${CUR_DIR}`
EXEC_WAVE=${EXE_DIR}/main_curv_col_el_3d
echo "EXEC_WAVE=$EXEC_WAVE"

#-- output and conf
PROJDIR=~/work/cgfd3d/sac_nc
EVTNM=codetest
echo "PROJDIR=${PROJDIR}"
echo "EVTNM=${EVTNM}"

PAR_FILE=${PROJDIR}/test.json
GRID_DIR=${PROJDIR}/fd_grid
METRIC_DIR=${PROJDIR}/fd_grid
MEDIA_DIR=${PROJDIR}/fd_medium
SOURCE_DIR=${PROJDIR}/fd_src
OUTPUT_DIR=${PROJDIR}/output

#RUN_SCRIPT_FILE=${PROJDIR}/runscript.sh
RUN_SCRIPT_FILE=${PROJDIR}/runscript.lsf
echo "RUN_SCRIPT_FILE=${RUN_SCRIPT_FILE}"

#-- create dir
mkdir -p $PROJDIR
mkdir -p $OUTPUT_DIR
mkdir -p $GRID_DIR
mkdir -p $METRIC_DIR
mkdir -p $MEDIA_DIR
mkdir -p $SOURCE_DIR

#-- tmp dir should be in local, here is for test
TMP_DIR=${PROJDIR}/scratch
mkdir -p ${TMP_DIR}

#----------------------------------------------------------------------
#-- grid and mpi configurations
#----------------------------------------------------------------------

#-- total x grid points
NX=200
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
1.0e16  1.0e16  1.0e16 0 0 0 
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
5
# name is_physical_coord is_3dim_depth  x y z
r1  1  1  1000 1000 0
r2  1  1  2000 2000 0
r3  1  1  3000 3000 0
r4  1  1  4000 4000 0
r5  1  1  5000 5000 0
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
  "time_window_length" : 5,
  "check_stability" : 1,

  "boundary_x_left" : {
      "ablexp" : {
          "number_of_layers" : 20,
          "ref_vel"  : 7000.0
          }
      },
  "boundary_x_right" : {
      "ablexp" : {
          "number_of_layers" : 20,
          "ref_vel"  : 7000.0
          }
      },
  "boundary_y_front" : {
      "ablexp" : {
          "number_of_layers" : 20,
          "ref_vel"  : 7000.0
          }
      },
  "boundary_y_back" : {
      "ablexp" : {
          "number_of_layers" : 20,
          "ref_vel"  : 7000.0
          }
      },
  "boundary_z_bottom" : {
      "ablexp" : {
          "number_of_layers" : 20,
          "ref_vel"  : 7000.0
          }
      },
  "boundary_z_top" : {
      "free" : "timg"
      },


  "#==": "set one and only one grid input method from following three ways",
  "#grid_import_from_previous_run" : "$GRID_DIR",
  "grid_generate_by_cartesian" : {
      "origin"  : [0.0, 0.0, -5900.0 ],
      "inteval" : [ 100.0, 100.0, 100.0 ]
      },
  "#grid_generate_by_gdlay" : {
        "file_name" : "${CUR_DIR}/prep_grid/tangshan_area_topo.gdlay",
        "refine_factor" : [ 1, 1, 1 ],
        "horizontal_start_index" : [ 3, 3 ],
        "vertical_last_to_top" : 0
  },
  "grid_export_to_dir"   : "$GRID_DIR",


  "#metric_import_from_previous_run" : "$METRIC_DIR",
  "metric_export_to_dir"   : "$METRIC_DIR",


  "medium_type" : "elastic_iso",

  "#==": "set one and only one medium input method from following five ways",
  "#medium_import_from_previous_run" : "$MEDIA_DIR",
  "medium_create_in_code" : "func_name_here",
  "#medium_read_from_binfile" : {
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
  "#medium_read_from_mdlay" : "${CUR_DIR}/prep_medium/basin_el_iso.md3lay",
  "#medium_read_from_mdgrd" : "${CUR_DIR}/prep_medium/basin_el_iso.md3grd",

  "medium_equivalent_method" : "loc",
  "#medium_equivalent_method" : "har",
  "media_export_to_dir"  : "$MEDIA_DIR",

  "#visco_config" : {
      "type" : "gmb",
      "Qs_freq" : 1.0,
      "number_of_maxwell" : 2,
      "max_freq" : 10.0,
      "min_freq" : 0.1,
      "refer_freq" : 1.0
  },

  "in_source_file" : "$PROJDIR/test.src",
  "source_export_to_dir"  : "$SOURCE_DIR",

  "#in_ddsource_file" : "${CUR_DIR}/prep_source/event_3moment_srcdd.nc",

  "!-- station or line output" : "--!",
  "receiver" : {
      "save_velocity" : 1,
      "save_stress"   : 1,
      "save_strain"   : 1,
      "time_step_per_save" : 10,

      "station_file"        : "$PROJDIR/test.station",
      "station_save_by_sac" : 0,

      "receiver_line" : [
        {
          "name" : "line_x_1",
          "grid_index_start"    : [  0, 49, 59 ],
          "grid_index_incre"    : [  10,  0,  0 ],
          "grid_index_count"    : 20
        },
        {
          "name" : "line_y_1",
          "grid_index_start"    : [ 19, 0, 59 ],
          "grid_index_incre"    : [  0, 10,  0 ],
          "grid_index_count"    : 30
        } 
      ]
  },

  "snapshot" : [
    {
      "name" : "snap_surf",
      "grid_index_start" : [ 0, 0, $(( NZ - 1 )) ],
      "grid_index_count" : [ $NX, $NY, 1 ],
      "grid_index_incre" : [  1, 1, 1 ],
      "time_index_start" : 0,
      "time_index_incre" : 1,
      "save_velocity" : 1,
      "save_stress"   : 1,
      "save_strain"   : 0,
      "save_coord"    : 1
    },
    {
      "name" : "snap_vol",
      "grid_index_start" : [ 0, 0, 0 ],
      "grid_index_count" : [ $NX, $NY, $NZ ],
      "grid_index_incre" : [  1, 1, 1 ],
      "time_index_start" : 0,
      "time_index_incre" : 1,
      "save_velocity" : 1,
      "save_stress"   : 1,
      "save_strain"   : 0,
      "save_coord"    : 1
    }
  ],

  "output_dir" : "$OUTPUT_DIR",
  "tmp_dir"    : "$TMP_DIR",

  "parallel_netcdf" : 1,
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
echo "sumbit to lsf ..."
create_script_lsf;
bsub < ${RUN_SCRIPT_FILE}

#-- run with pbs
#echo "sumbit to pbs ..."
#create_script_pbs;
#qsub ${RUN_SCRIPT_FILE}

#-- directly run
#echo "start run script ..."
#create_script_run;
#${RUN_SCRIPT_FILE}

date

# vim:ts=4:sw=4:nu:et:ai:
