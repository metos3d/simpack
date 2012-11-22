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

#undef  kDebugLevel
#define kDebugLevel kDebugLevel1

#undef  __FUNCT__
#define __FUNCT__ "Metos3DInitWithFilePath"
PetscErrorCode
Metos3DInitWithFilePath(Metos3D *metos3d, const char *filePath)
{
    PetscFunctionBegin;
    // debug start, kDebugLevel1
    PetscGetTime(&metos3d->startTime[kDebugLevel1]);
    // store communicator
    // read in options file
    metos3d->comm = PETSC_COMM_WORLD;
    PetscOptionsInsertFile(metos3d->comm, filePath, PETSC_TRUE);
    // read in debug option
    metos3d->debugLevel = kDebugLevel0;
    Metos3DUtilOptionsGetInt(metos3d, "-Metos3DDebugLevel", &metos3d->debugLevel);
    // set level 0 as lowest possible
    // set level 4 as highest possible
    metos3d->debugLevel = (metos3d->debugLevel < kDebugLevel0 ? kDebugLevel0 : metos3d->debugLevel);
    metos3d->debugLevel = (metos3d->debugLevel > kDebugLevel4 ? kDebugLevel4 : metos3d->debugLevel);
    // init geometry, load, bgc, trasport, timestep, solver, ...
    Metos3DGeometryInit(metos3d);
    Metos3DLoadInit(metos3d);
    Metos3DBGCInit(metos3d);
    Metos3DTransportInit(metos3d);
    Metos3DTimeStepInit(metos3d);
    Metos3DSolverInit(metos3d);
    // wait for all processors
    PetscBarrier(PETSC_NULL);
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, F3S, "Metos3DInitWithFilePath", "filePath:", filePath);
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DFinal"
PetscErrorCode
Metos3DFinal(Metos3D *metos3d)
{
    PetscFunctionBegin;
//    // debug start
//    PetscGetTime(&metos3d->startTime[kDebugLevel]);
    // wait for all processors
    PetscBarrier(PETSC_NULL);
    // final ..., solver, timestep, transport, bgc, load, geometry
    Metos3DSolverFinal(metos3d);
    Metos3DTimeStepFinal(metos3d);
    Metos3DTransportFinal(metos3d);
    Metos3DBGCFinal(metos3d);    
    Metos3DLoadFinal(metos3d);
    Metos3DGeometryFinal(metos3d);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DFinal\n");
    PetscFunctionReturn(0);
}
