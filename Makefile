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

#-- 
CFLAGS := -I$(NETCDF)/include -I./lib/ -I./forward/ -I./media/  $(CFLAGS)

#- debug
#CFLAGS   := -g $(CFLAGS)
#CPPFLAGS := -g -std=c++11 $(CPPFLAGS)
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

default: main_curv_col_el_3d
#default: main_curv_col_ac_3d
#default: main_cart_col_el_3d
#default: main_cart_stg_el_3d
#default: main_cart_stg_ac_3d

all: main_curv_col_el_3d \
     main_curv_col_ac_3d \
     main_cart_col_el_3d \
     main_cart_stg_el_3d \
     main_cart_stg_ac_3d

main_curv_col_el_3d: \
		cJSON.o sacLib.o fdlib_mem.o fdlib_math.o  \
		fd_t.o par_t.o interp.o mympi_t.o \
		media_utility.o \
		media_layer2model.o \
		media_grid2model.o \
		media_geometry3d.o \
		media_read_file.o \
		gd_info.o gd_t.o md_t.o wav_t.o \
		bdry_free.o bdry_pml.o src_t.o io_funcs.o \
		blk_t.o \
		sv_eq1st_curv_col.o \
		sv_eq1st_curv_col_ac_iso.o \
		sv_eq1st_curv_col_el_aniso.o sv_eq1st_curv_col_el_iso.o \
		main_curv_col_el_3d.o
	$(CXX) -o $@ $^ $(LDFLAGS)

main_curv_col_ac_3d: \
		cJSON.o sacLib.o fdlib_mem.o fdlib_math.o  \
		fd_t.o par_t.o interp.o mympi_t.o \
		media_utility.o \
		media_layer2model.o \
		media_grid2model.o \
		media_geometry3d.o \
		media_read_file.o \
		gd_info.o gd_t.o md_t.o wav_t.o \
		bdry_free.o bdry_pml.o src_t.o io_funcs.o \
		blk_t.o \
		sv_eq1st_curv_col.o \
		sv_eq1st_curv_col_ac_iso.o \
		sv_eq1st_curv_col_el_aniso.o sv_eq1st_curv_col_el_iso.o \
		main_curv_col_ac_3d.o
	$(CXX) -o $@ $^ $(LDFLAGS)

main_cart_col_el_3d: \
		cJSON.o sacLib.o fdlib_mem.o fdlib_math.o  \
		fd_t.o par_t.o interp.o mympi_t.o \
		media_utility.o \
		media_layer2model.o \
		media_grid2model.o \
		media_geometry3d.o \
		media_read_file.o \
		gd_info.o gd_t.o md_t.o wav_t.o \
		bdry_free.o bdry_pml.o src_t.o io_funcs.o \
		blk_t.o \
		sv_eq1st_cart_col.o \
		sv_eq1st_cart_col_el_iso.o \
		main_cart_col_el_3d.o
	$(CXX) -o $@ $^ $(LDFLAGS)

main_cart_stg_el_3d: \
		cJSON.o sacLib.o fdlib_mem.o fdlib_math.o  \
		fd_t.o par_t.o interp.o mympi_t.o \
		media_utility.o \
		media_layer2model.o \
		media_grid2model.o \
		media_geometry3d.o \
		media_read_file.o \
		gd_info.o gd_t.o md_t.o wav_t.o \
		bdry_free.o bdry_pml.o src_t.o io_funcs.o \
		blk_t.o \
		sv_eq1st_cart_stg_el_iso.o \
		main_cart_stg_el_3d.o
	$(CXX) -o $@ $^ $(LDFLAGS)

main_cart_stg_ac_3d: \
		cJSON.o sacLib.o fdlib_mem.o fdlib_math.o  \
		fd_t.o par_t.o interp.o mympi_t.o \
		media_utility.o \
		media_layer2model.o \
		media_grid2model.o \
		media_geometry3d.o \
		media_read_file.o \
		gd_info.o gd_t.o md_t.o wav_t.o \
		bdry_free.o bdry_pml.o src_t.o io_funcs.o \
		blk_t.o \
		sv_eq1st_cart_stg_ac_iso.o \
		main_cart_stg_ac_3d.o
	$(CXX) -o $@ $^ $(LDFLAGS)

media_geometry3d.o: media/media_geometry3d.cpp 
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_utility.o: media/media_utility.cpp 
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_layer2model.o: media/media_layer2model.cpp
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_grid2model.o: media/media_grid2model.cpp
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_read_file.o: media/media_read_file.cpp
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
mympi_t.o: forward/mympi_t.c
	${CC} -c -o $@ $(CFLAGS) $<
gd_info.o: forward/gd_info.c
	${CC} -c -o $@ $(CFLAGS) $<
gd_t.o: forward/gd_t.c
	${CC} -c -o $@ $(CFLAGS) $<
md_t.o: forward/md_t.c
	${CC} -c -o $@ $(CFLAGS) $<
wav_t.o: forward/wav_t.c
	${CC} -c -o $@ $(CFLAGS) $<
bdry_pml.o: forward/bdry_pml.c
	${CC} -c -o $@ $(CFLAGS) $<
bdry_free.o: forward/bdry_free.c
	${CC} -c -o $@ $(CFLAGS) $<
src_t.o: forward/src_t.c
	${CC} -c -o $@ $(CFLAGS) $<
io_funcs.o: forward/io_funcs.c
	${CC} -c -o $@ $(CFLAGS) $<
blk_t.o: forward/blk_t.c
	${CC} -c -o $@ $(CFLAGS) $<
sv_eq1st_curv_col.o:          forward/sv_eq1st_curv_col.c
	${CC} -c -o $@ $(CFLAGS) $<
sv_eq1st_curv_col_el_iso.o:   forward/sv_eq1st_curv_col_el_iso.c
	${CC} -c -o $@ $(CFLAGS) $<
sv_eq1st_curv_col_el_aniso.o: forward/sv_eq1st_curv_col_el_aniso.c
	${CC} -c -o $@ $(CFLAGS) $<
sv_eq1st_curv_col_ac_iso.o:   forward/sv_eq1st_curv_col_ac_iso.c
	${CC} -c -o $@ $(CFLAGS) $<
sv_eq1st_cart_col.o:          forward/sv_eq1st_cart_col.c
	${CC} -c -o $@ $(CFLAGS) $<
sv_eq1st_cart_col_el_iso.o:   forward/sv_eq1st_cart_col_el_iso.c
	${CC} -c -o $@ $(CFLAGS) $<
sv_eq1st_cart_stg_el_iso.o:   forward/sv_eq1st_cart_stg_el_iso.c
	${CC} -c -o $@ $(CFLAGS) $<
sv_eq1st_cart_stg_ac_iso.o:   forward/sv_eq1st_cart_stg_ac_iso.c
	${CC} -c -o $@ $(CFLAGS) $<

main_curv_col_el_3d.o: forward/main_curv_col_el_3d.c
	${CC} -c -o $@ $(CFLAGS) $<
main_curv_col_ac_3d.o: forward/main_curv_col_ac_3d.c
	${CC} -c -o $@ $(CFLAGS) $<
main_cart_col_el_3d.o: forward/main_cart_col_el_3d.c
	${CC} -c -o $@ $(CFLAGS) $<
main_cart_stg_el_3d.o: forward/main_cart_stg_el_3d.c
	${CC} -c -o $@ $(CFLAGS) $<
main_cart_stg_ac_3d.o: forward/main_cart_stg_ac_3d.c
	${CC} -c -o $@ $(CFLAGS) $<

cleanexe:
	rm -f main_curv_col_el_3d
	rm -f main_curv_col_ac_3d
	rm -f main_cart_col_el_3d
	rm -f main_cart_stg_el_3d
	rm -f main_cart_stg_ac_3d
cleanobj:
	rm -f *.o
cleanall: cleanexe cleanobj
	echo "clean all"
distclean: cleanexe cleanobj
	echo "clean all"
