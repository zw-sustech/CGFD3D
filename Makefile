################################################################################
#  Makefile for CGFDM3D package
#
#  Author: ZHANG Wei <zhangwei@sustech.edu.cn>
#  Copyright (C) Wei ZHANG, 2020. All Rights Reserved.
#
################################################################################

#-------------------------------------------------------------------------------
# compiler
#-------------------------------------------------------------------------------

CC     :=  /share/apps/gnu-4.8.5/mpich-3.3/bin/mpicc
NETCDF :=  /share/apps/gnu-4.8.5/disable-netcdf-4.4.1

#-- 
CFLAGS := -I$(NETCDF)/include -I./lib/ -I./forward/ $(CFLAGS)

#- debug
CFLAGS   := -g $(CFLAGS)
#- O3
#CFLAGS   := -O3 $(CFLAGS)

#- static
#LDFLAGS := $(NETCDF)/lib/libnetcdf.a -lm -static $(LDFLAGS)
LDFLAGS := -lm $(LDFLAGS) $(NETCDF)/lib/libnetcdf.a
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
		cJSON.o fdlib_mem.o fdlib_math.o  \
		fd_t.o par_t.o \
		gd_curv.o md_el_iso.o wf_el_1st.o \
		abs_funcs.o src_funcs.o \
		sv_eliso1st_curv_macdrp.o \
		cgfdm3d_elastic_main.o
	$(CC) -o $@ $^ $(LDFLAGS)

cJSON.o: lib/cJSON.c
	${CC} -c -o $@ $(CFLAGS) $<
fdlib_mem.o: lib/fdlib_mem.c
	${CC} -c -o $@ $(CFLAGS) $<
fdlib_math.o: lib/fdlib_math.c
	${CC} -c -o $@ $(CFLAGS) $<
fd_t.o: forward/fd_t.c
	${CC} -c -o $@ $(CFLAGS) $<
par_t.o: forward/par_t.c
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
