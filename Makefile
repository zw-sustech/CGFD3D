################################################################################
#  Makefile for CGFDM3D package
#
#  Author: ZHANG Wei <zhangwei@sustech.edu.cn>
#  Copyright (C) Wei ZHANG, 2020. All Rights Reserved.
#
################################################################################

# known issues:
#  .h dependency is not included, use make cleanall

#-------------------------------------------------------------------------------
# compiler
#-------------------------------------------------------------------------------

CC     :=  /share/apps/gnu-4.8.5/mpich-3.3/bin/mpicc
CXX    :=  /share/apps/gnu-4.8.5/mpich-3.3/bin/mpicxx
NETCDF :=  /share/apps/gnu-4.8.5/disable-netcdf-4.4.1

# CC     :=  /home/zl/software/mpich3.2/bin/mpicc
# CXX    :=  /home/zl/software/mpich3.2/bin/mpicxx
# NETCDF :=  /home/zl/software/netcdf-4.4.1

#-- 
CFLAGS := -I$(NETCDF)/include -I./lib/ -I./forward/ -I./media/ $(CFLAGS)

#- debug
#CFLAGS   := -g $(CFLAGS)
#- O3
CFLAGS   := -O3 $(CFLAGS)
CPPFLAGS := -O2 -std=c++11 $(CPPFLAGS)

#- static
#LDFLAGS := $(NETCDF)/lib/libnetcdf.a -lm -static $(LDFLAGS)
LDFLAGS := -lm  $(LDFLAGS) $(NETCDF)/lib/libnetcdf.a
#- dynamic
#LDFLAGS := -L$(NETCDF)/lib -lnetcdf -lm $(LDFLAGS)

#-------------------------------------------------------------------------------
# target
#-------------------------------------------------------------------------------

# special vars:
# 	$@ The file name of the target
# 	$< The names of the first prerequisite
#   $^ The names of all the prerequisites 

cgfdm3d_elastic_mpi: \
		cJSON.o sacLib.o fdlib_mem.o fdlib_math.o  \
		fd_t.o par_t.o interp.o\
		media_layer2model.o \
		media_geometry3d.o \
	  media_interpolation.o \
		media_read_interface_file.o \
		gd_curv.o md_el_iso.o wf_el_1st.o \
		abs_funcs.o src_funcs.o io_funcs.o \
		sv_eliso1st_curv_macdrp.o \
		cgfdm3d_elastic_main.o
	$(CXX) -o $@ $^ $(LDFLAGS)


media_geometry3d.o: media/media_geometry3d.cpp 
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_layer2model.o: media/media_layer2model.cpp
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_interpolation.o: media/media_interpolation.cpp
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_read_interface_file.o: media/media_read_interface_file.cpp
	${CXX} -c -o $@ $(CPPFLAGS) $<
cJSON.o: lib/cJSON.c
	${CC} -c -o $@ $(CFLAGS) $<
sacLib.o: lib/sacLib.c
	${CC} -c -o $@ $(CFLAGS) $<
fdlib_mem.o: lib/fdlib_mem.c
	${CC} -c -o $@ $(CFLAGS) $<
fdlib_math.o: lib/fdlib_math.c
	${CC} -c -o $@ $(CFLAGS) $<
fd_t.o: forward/fd_t.c
	${CC} -c -o $@ $(CFLAGS) $<
par_t.o: forward/par_t.c
	${CC} -c -o $@ $(CFLAGS) $<
interp.o: forward/interp.c
	${CC} -c -o $@ $(CFLAGS) $<
gd_curv.o: forward/gd_curv.c
	${CC} -c -o $@ $(CFLAGS) $<
md_el_iso.o: forward/md_el_iso.c
	${CC} -c -o $@ $(CFLAGS) $<
wf_el_1st.o: forward/wf_el_1st.c
	${CC} -c -o $@ $(CFLAGS) $<
abs_funcs.o: forward/abs_funcs.c
	${CC} -c -o $@ $(CFLAGS) $<
src_funcs.o: forward/src_funcs.c
	${CC} -c -o $@ $(CFLAGS) $<
io_funcs.o: forward/io_funcs.c
	${CC} -c -o $@ $(CFLAGS) $<
sv_eliso1st_curv_macdrp.o: forward/sv_eliso1st_curv_macdrp.c
	${CC} -c -o $@ $(CFLAGS) $<
cgfdm3d_elastic_main.o: forward/cgfdm3d_elastic_main.c
	${CC} -c -o $@ $(CFLAGS) $<

cleanexe:
	rm -f cgfdm3d_elastic_mpi
cleanobj:
	rm -f *.o
cleanall: cleanexe cleanobj
	echo "clean all"
distclean: cleanexe cleanobj
	echo "clean all"
