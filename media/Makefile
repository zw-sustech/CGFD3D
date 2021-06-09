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

CC     :=  mpicc
CXX    :=  mpicxx

#-- 
#CFLAGS := -I$(NETCDF)/include -I./lib/ -I./forward/ -I./media/ $(CFLAGS)

#- debug
CFLAGS   := -g $(CFLAGS)
#- O3
CFLAGS   := -O3 $(CFLAGS)
CPPFLAGS := -O2 -std=c++11 $(CPPFLAGS)

#- static
#LDFLAGS := $(NETCDF)/lib/libnetcdf.a -lm -static $(LDFLAGS)
#LDFLAGS := -lm  $(LDFLAGS) $(NETCDF)/lib/libnetcdf.a
LDFLAGS := -lm  
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
	test_isPointInHexahedron.o \
	media_geometry3d.o
	$(CXX) -o $@ $^ $(LDFLAGS)


media_geometry3d.o: media_geometry3d.cpp 
	${CXX} -c -o $@ $(CPPFLAGS) $<
test_isPointInHexahedron.o: test_isPointInHexahedron.c 
	${CC} -c -o $@ $(CFLAGS) $<

cleanexe:
	rm -f test
cleanobj:
	rm -f *.o
cleanall: cleanexe cleanobj
	echo "clean all"
distclean: cleanexe cleanobj
	echo "clean all"
