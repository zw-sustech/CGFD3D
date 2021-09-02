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

#CC     :=  /share/apps/gnu-4.8.5/mpich-3.3/bin/mpicc
#CXX    :=  /share/apps/gnu-4.8.5/mpich-3.3/bin/mpicxx
#NETCDF :=  /share/apps/gnu-4.8.5/disable-netcdf-4.4.1
CXX := g++
CC := gcc

#-- 
#CFLAGS := -I$(NETCDF)/include -I./lib/ -I./forward/ -I./media/ $(CFLAGS)

#- debug
#CFLAGS   := -g $(CFLAGS)
#- O3
CFLAGS   := -O3 $(CFLAGS)
CPPFLAGS := -g -O3 -std=c++11 -I../ $(CPPFLAGS)

#- static
#LDFLAGS := $(NETCDF)/lib/libnetcdf.a -lm -static $(LDFLAGS)
#LDFLAGS := -lm  $(LDFLAGS) $(NETCDF)/lib/libnetcdf.a
#- dynamic
#LDFLAGS := -L$(NETCDF)/lib -lnetcdf -lm $(LDFLAGS)

#-------------------------------------------------------------------------------
# target
#-------------------------------------------------------------------------------

# special vars:
# 	$@ The file name of the target
# 	$< The names of the first prerequisite
#   $^ The names of all the prerequisites 

test: \
	test.o \
	media_utility.o \
	media_layer2model.o \
	media_grid2model.o \
	media_geometry3d.o \
	media_read_file.o 
	$(CXX) -o $@ $^ $(LDFLAGS)


test.o: test.cpp 
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_geometry3d.o: ../media_geometry3d.cpp 
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_utility.o: ../media_utility.cpp 
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_layer2model.o: ../media_layer2model.cpp
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_grid2model.o: ../media_grid2model.cpp
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_read_file.o: ../media_read_file.cpp
	${CXX} -c -o $@ $(CPPFLAGS) $<

cleanexe:
	rm -f test
cleanobj:
	rm -f *.o
cleanall: cleanexe cleanobj
	echo "clean all"
distclean: cleanexe cleanobj
	echo "clean all"
