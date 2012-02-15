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
#   Makefile
#

# object files
OBJSC = \
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
BGCWORK = $(BGC:%/=%)
BGCMODELNAME = $(BGCWORK:model/%=%)
BGCMODELFILE = model.o
OBJSBGC = $(addprefix $(BGCWORK)/, $(BGCMODELFILE))

# executable name
PROGRAMBASE = metos3d-simpack-
PROGRAM = ${PROGRAMBASE}${BGCMODELNAME}.exe

ALL: $(PROGRAM)
CFLAGS = -DBGC=metos3dbgc_ -DBGCINIT=metos3dbgcinit_ -DBGCFINAL=metos3dbgcfinal_
FFLAGS = 
CLEANFILES = $(OBJSBGC) $(OBJSC) $(PROGRAM)

ifdef PETSC_DIR
ifdef BGC
include $(PETSC_DIR)/conf/variables
include $(PETSC_DIR)/conf/rules
$(PROGRAM): $(OBJSBGC) $(OBJSC) chkopts
	-$(CLINKER) -o $@ $(OBJSBGC) $(OBJSC) $(PETSC_LIB)
else
$(PROGRAM):
	@echo '###'
	@echo '### Please choose a BGC model.'
	@echo '###'
endif
else
$(PROGRAM):
	@echo '###'
	@echo '### Please set the PETSc environment variables.'
	@echo '###'
endif
