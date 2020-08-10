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
CFLAGS := -c -I$(NETCDF)/include $(CFLAGS)

#- debug
CFLAGS   := $(if $(WithOMP),-O0,-g -fcheck=all -finit-integer=0 -finit-real=zero) $(CFLAGS)
#- O3
#CFLAGS   := -O3 $(CFLAGS)

#- static
#LDFLAGS := $(NETCDF)/lib/libnetcdf.a -static $(LDFLAGS)
#- dynamic
LDFLAGS := -L$(NETCDF)/lib -lnetcdf $(LDFLAGS)

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
	$(CC) -o $@ $(LDFLAGS) $^

cJSON.o: lib/cJSON.c
	${CC} -c $@ $(CFLAGS) $<
fdlib_mem.o: lib/fdlib_mem.c
	${CC} -c $@ $(CFLAGS) $<
fdlib_math.o: lib/fdlib_math.c
	${CC} -c $@ $(CFLAGS) $<
fd_t.o: forward/fd_t.c
	${CC} -c $@ $(CFLAGS) $<
par_t.o: forward/par_t.c
	${CC} -c $@ $(CFLAGS) $<
gd_curv.o: forward/gd_curv.c
	${CC} -c $@ $(CFLAGS) $<
md_el_iso.o: forward/md_el_iso.c
	${CC} -c $@ $(CFLAGS) $<
wf_el_1st.o: forward/wf_el_1st.c
	${CC} -c $@ $(CFLAGS) $<
abs_funcs.o: forward/abs_funcs.c
	${CC} -c $@ $(CFLAGS) $<
src_funcs.o: forward/src_funcs.c
	${CC} -c $@ $(CFLAGS) $<
sv_eliso1st_curv_macdrp.o: forward/sv_eliso1st_curv_macdrp.c
	${CC} -c $@ $(CFLAGS) $<
cgfdm3d_elastic_main.o: forward/cgfdm3d_elastic_main.c
	${CC} -c $@ $(CFLAGS) $<

cleanexe:
	rm -f cgfdm3d_elastic_mpi
cleanobj:
	rm -f *.o
cleanall: cleanexe cleanobj
	echo "clean all"
distclean: cleanexe cleanobj
	echo "clean all"
