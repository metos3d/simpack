/*
 * Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D
 * Copyright (C) 2012  Jaroslaw Piwonski, CAU, jpi@informatik.uni-kiel.de
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  metos3d_init.c
 *
 */

#include "metos3d_init.h"

#pragma mark -
#undef  kDebugLevel
#define kDebugLevel kDebugLevel1

#undef  __FUNCT__
#define __FUNCT__ "Metos3DInitWithFilePath"
PetscErrorCode
Metos3DInitWithFilePath(Metos3D *metos3d, const char *filePath)
{
    MPI_Comm    comm = PETSC_COMM_WORLD;
    PetscInt    debug;
    PetscFunctionBegin;
    // store communicator
    metos3d->comm = comm;
    // read in options file
    PetscOptionsInsertFile(comm, filePath, PETSC_TRUE);
    // read in debug option
    metos3d->debugLevel = kDebugLevel0;
    Metos3DUtilOptionsGetInt(metos3d, "-Metos3DDebugLevel", &debug);
    // store user debug level
    metos3d->debugLevel = debug;
    Metos3DDebug(debug, kDebugLevel, comm, F3S, "Metos3DInitWithFilePath", "filePath:", filePath);
    // init geometry, load, bgc, trasport, timestep, solver, ...
    Metos3DGeometryInit(metos3d);
    Metos3DLoadInit(metos3d);
    Metos3DBGCInit(metos3d);
    Metos3DTransportInit(metos3d);
    Metos3DTimeStepInit(metos3d);
    Metos3DSolverInit(metos3d);
    // wait for all processors  
    PetscBarrier(PETSC_NULL);
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DFinal"
PetscErrorCode
Metos3DFinal(Metos3D *metos3d)
{
    MPI_Comm    comm    = metos3d->comm;
    PetscInt    debug   = metos3d->debugLevel;
    PetscFunctionBegin;
    Metos3DDebug(debug, kDebugLevel, comm, "Metos3DFinal\n");
    // wait for all processors  
    PetscBarrier(PETSC_NULL);
    // final ..., solver, timestep, transport, bgc, load, geometry
    Metos3DSolverFinal(metos3d);
    Metos3DTimeStepFinal(metos3d);
    Metos3DTransportFinal(metos3d);
    Metos3DBGCFinal(metos3d);    
    Metos3DLoadFinal(metos3d);
    Metos3DGeometryFinal(metos3d);
    PetscFunctionReturn(0);
}
