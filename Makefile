# Makefile for SPECFEM3D_GEOTECH
# REVISION
#   HNG, Jan 17,2014

default: semgeotech

all: createdir partmesh semgeotech psemgeotech

createdir:
	(mkdir -p bin; mkdir -p input; mkdir -p output; mkdir -p partition; mkdir -p tmp)

clean:
	(cd src; make clean)

partmesh: 
	(cd src; make $@)

semgeotech: 
	(cd src; make $@)

psemgeotech: 
	(cd src; make $@)

