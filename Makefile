#
# Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D
# Copyright (C) 2012  Jaroslaw Piwonski, CAU, jpi@informatik.uni-kiel.de
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# object files
M3DOBJSC = \
	src/metos3d_debug.o \
	src/metos3d_util.o \
	src/metos3d_geometry.o \
	src/metos3d_load.o \
	src/metos3d_bgc.o \
	src/metos3d_transport.o \
	src/metos3d_timestep.o \
	src/metos3d_solver.o \
	src/metos3d_init.o \
	src/metos3d.o

# BGC model name
BGCMODELNAME = $(notdir $(BGC:%/=%))
BGCMODELFILE = model.o
M3DOBJSBGC = $(addprefix $(BGC)/, $(BGCMODELFILE))

# executable name
PROGRAMBASE = metos3d-simpack-
PROGRAM = ${PROGRAMBASE}${BGCMODELNAME}.exe

ALL: $(PROGRAM)
CFLAGS = -DBGC=metos3dbgc_ -DBGCINIT=metos3dbgcinit_ -DBGCFINAL=metos3dbgcfinal_
FFLAGS = 
CLEANFILES = $(M3DOBJSBGC) $(M3DOBJSC) $(PROGRAM)

ifdef PETSC_DIR
ifdef BGC
include $(PETSC_DIR)/conf/variables
include $(PETSC_DIR)/conf/rules
$(PROGRAM): $(M3DOBJSBGC) $(M3DOBJSC) chkopts
	-$(CLINKER) -o $@ $(M3DOBJSBGC) $(M3DOBJSC) $(PETSC_LIB)
else
$(PROGRAM):
	@echo '###'
	@echo '### Please choose a BGC model.'
	@echo '###'
	@echo '###   Recommended:'
	@echo '###'
	@echo '###     1. Link the parent model directory into the current directory.'
	@echo '###     2. Tell the Makefile, where to find your model.'
	@echo '###     3. Run the executable with an option file.'
	@echo '###'
	@echo '###   Example:'
	@echo '###'
	@echo '###     $$> ln -s ../model/ model'
	@echo '###     $$> make BGC=model/I-Cs/'
	@echo '###     $$> ./metos3d-simpack-I-Cs.exe model/I-Cs/option/local.jserver.option.I-Cs.simpack.txt'
	@echo '###'
endif
else
$(PROGRAM):
	@echo '###'
	@echo '### Please set the PETSc environment (PETSC_DIR, PETSC_ARCH) variables.'
	@echo '###'
	@echo '###   Best practice:'
	@echo '###'
	@echo '###     1. Prepare a file with export commands and store it in the "petsc" directory.'
	@echo '###     2. Use the shell source operator "." to set the variables.'
	@echo '###'
	@echo '###   Example:'
	@echo '###'
	@echo '###     $$> . petsc/local.jserver.petsc.opt.txt'
	@echo '###'
endif
