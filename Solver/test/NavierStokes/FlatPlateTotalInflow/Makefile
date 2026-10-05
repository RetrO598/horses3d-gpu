#####################################################################
## HORSES - Laminar flat plate with a total-condition inlet
#####################################################################
HOME_DIR?=$(abspath ../../..)
SOLVER?=$(HOME_DIR)/bin/horses3d.ns
PYTHON?=python3

.DEFAULT_GOAL := setup
.PHONY: setup mesh run smoke clean

setup:
	$(MAKE) -C SETUP ns HOME_DIR=$(HOME_DIR)
	mkdir -p RESULTS

mesh: MESH/flatplate_65x65x1.msh

MESH/flatplate_65x65x1.msh: MESH/mesh_flatplate_65x65.su2 SETUP/extrude_su2.py
	$(PYTHON) SETUP/extrude_su2.py $< $@ --span 0.003048

run: setup mesh
	$(SOLVER) FlatPlate.control

smoke: setup mesh
	sed -e 's/Number of time steps             = 100000/Number of time steps             = 5/' \
	    -e 's/Output Interval                  = 1000/Output Interval                  = 5/' \
	    -e 's@RESULTS/FlatPlate.hsol@RESULTS/Smoke.hsol@g' FlatPlate.control > Smoke.control
	$(SOLVER) Smoke.control

clean:
	$(MAKE) -C SETUP allclean
	rm -f Smoke.control
