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
 *  metos3d_type.h
 *
 */

#ifndef METOS3D_TYPE_H
#define METOS3D_TYPE_H 1

#include "petscsnes.h"

#define max(a, b) (((a) > (b)) ? (a) : (b))
//#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
//#define sign( x ) ( ((x) > 0.0) ? 1.0 : ( ((x) < 0.0) ? -1.0 : 0.0 ) )

/*
 *  Metos3D context data type.
 */
typedef struct {
    // communication
    MPI_Comm        comm;
    // debug
    PetscInt        debugLevel;
    PetscLogDouble  startTime[4];
    PetscInt        fileFormatPrefixCount;
    char            *fileFormatPrefix[PETSC_MAX_PATH_LEN];
    char            *filePrefix;
    PetscInt        moduloStepCount;
    PetscInt        *moduloStep;
    // geometry
    PetscInt        profileCount;
    PetscInt        *profileStart;
    PetscInt        *profileEnd;
    PetscInt        profileLengthMax;
    PetscInt        vectorLength;
    // load
    PetscInt        processCount;
    PetscInt        processMine;
    PetscInt        vectorLengthLocal;
    PetscInt        vectorLengthPrevious;
    PetscInt        profileCountLocal;
    PetscInt        profileCountPrevious;
    PetscInt        *profileStartLocal;
    PetscInt        *profileEndLocal;
    // bgc
    PetscInt        tracerCount;
    Vec             *y0;
    Vec             *ybgcinBD;
    Vec             *ybgcoutBD;
    PetscInt        parameterCount;
    PetscReal       *u0;
    PetscInt        boundaryConditionCount;
    PetscInt        *boundaryConditionCountArray;
    Vec             **boundaryConditionArray;
    Vec             *bgcbc;
    Vec             *bgcbcBD;
    PetscInt        domainConditionCount;
    PetscInt        *domainConditionCountArray;
    Vec             **domainConditionArray;
    Vec             *bgcdc;
    Vec             *bgcdcBD;
    // transport
    PetscInt        matrixCount;
    Mat             *matrixExplicitArray;
    Mat             *matrixImplicitArray;
    Mat             Aework;
    Mat             Aiwork;
    // time step
    PetscScalar     timeStep;
    PetscScalar     timeStepStart;
    PetscInt        timeStepCount;
    // solver
    SNES            snes;
    Vec             *fnBD;
    Mat             JShell;
} Metos3D;

#endif /* !METOS3D_TYPE_H */

