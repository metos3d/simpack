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
 *  metos3d_system.c
 *
 */

#include "metos3d_geometry.h"

#undef  kDebugLevel
#define kDebugLevel kDebugLevel2

#undef  __FUNCT__
#define __FUNCT__ "Metos3DGeometryInit"
PetscErrorCode
Metos3DGeometryInit(Metos3D *metos3d)
{
    MPI_Comm    comm = metos3d->comm;
    // work vars
    
//    char        geometryInputDirectory  [PETSC_MAX_PATH_LEN];
//    char        profileFile             [PETSC_MAX_PATH_LEN];
//    char        volumeFile              [PETSC_MAX_PATH_LEN];
//    char        filePath                [PETSC_MAX_PATH_LEN];
//
//    PetscFunctionBegin;
//    Metos3DDebug(metos3d, kDebugLevel, "Metos3DGeometryInit\n");
//
//    // geometry
//    Metos3DUtilOptionsGetString(metos3d, "-Metos3DGeometryInputDirectory", geometryInputDirectory);
//    Metos3DUtilOptionsGetString(metos3d, "-Metos3DProfileFile", profileFile);
//    Metos3DUtilOptionsGetString(metos3d, "-Metos3DVolumeFile", volumeFile);
//    
//    Metos3DFlag(PETSC_FALSE, geometryInputDirectory);

    char        geometryType[PETSC_MAX_PATH_LEN];    
    PetscBool   flag;
    PetscFunctionBegin;
//    // debug start
//    PetscGetTime(&metos3d->startTime[kDebugLevel]);
    // geometry type
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DGeometryType", geometryType);
    PetscStrcmp("Profile", geometryType, &flag);
    if (flag == PETSC_TRUE) {
        // work vars
        char        geometryInputDirectory  [PETSC_MAX_PATH_LEN];    
        char        profileStartFile        [PETSC_MAX_PATH_LEN];
        char        profileEndFile          [PETSC_MAX_PATH_LEN];
        char        filePath                [PETSC_MAX_PATH_LEN];
        PetscViewer viewer;
        PetscInt    profileCount;
        PetscInt    vectorLength;
        PetscInt    nlayermax, iprof;
        // profiles
        Metos3DUtilOptionsGetString(metos3d, "-Metos3DProfileInputDirectory", geometryInputDirectory);
        Metos3DUtilOptionsGetString(metos3d, "-Metos3DProfileIndexStartFile", profileStartFile);
        Metos3DUtilOptionsGetString(metos3d, "-Metos3DProfileIndexEndFile", profileEndFile);
        // get start indices (1 indexed)
        // concat file path
        sprintf(filePath, "%s%s", geometryInputDirectory, profileStartFile);
//        Metos3DDebug(metos3d, kDebugLevel, F3S, "Metos3DGeometryInit", "filePath:", filePath);
        PetscViewerBinaryOpen(comm, filePath, FILE_MODE_READ, &viewer);
        // read profile count
        PetscViewerBinaryRead(viewer, (void*)&profileCount, 1, PETSC_INT);
        // store
        metos3d->profileCount = profileCount;
//        Metos3DDebug(metos3d, kDebugLevel, F2SD, "Metos3DGeometryInit", "profileCount:", profileCount);
        PetscMalloc(profileCount*sizeof(PetscInt), &metos3d->profileStart);
        // read all indices
        PetscViewerBinaryRead(viewer, (void*)metos3d->profileStart, profileCount, PETSC_INT);
        PetscViewerDestroy(&viewer);
        // get end indices (1 indexed)
        // concat file path
        sprintf(filePath, "%s%s", geometryInputDirectory, profileEndFile);
//        Metos3DDebug(metos3d, kDebugLevel, F3S, "Metos3DGeometryInit", "filePath:", filePath);
        PetscViewerBinaryOpen(comm, filePath, FILE_MODE_READ, &viewer);
        // read profile count
        PetscViewerBinaryRead(viewer, (void*)&profileCount, 1, PETSC_INT);
        // store
        metos3d->profileCount = profileCount;
//        Metos3DDebug(metos3d, kDebugLevel, F2SD, "Metos3DGeometryInit", "profileCount:", profileCount);
        PetscMalloc(profileCount*sizeof(PetscInt), &metos3d->profileEnd);
        // read all indices
        PetscViewerBinaryRead(viewer, (void*)metos3d->profileEnd, profileCount, PETSC_INT);
        PetscViewerDestroy(&viewer);
        // determine vector length 
        vectorLength = metos3d->profileEnd[profileCount-1];
        // store
        metos3d->vectorLength = vectorLength;
//        Metos3DDebug(metos3d, kDebugLevel, F2SD, "Metos3DGeometryInit", "vectorLength:", vectorLength);
        // determine max length of profile
        nlayermax = 0;
        for (iprof = 0; iprof < profileCount; iprof++) {
            nlayermax = max(nlayermax, metos3d->profileEnd[iprof]-metos3d->profileStart[iprof]+1);
        }
        // store
        metos3d->profileLengthMax = nlayermax;
//        Metos3DDebug(metos3d, kDebugLevel, F2SD, "Metos3DGeometryInit", "profileLengthMax:", nlayermax);
    }
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DGeometryInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DGeometryFinal"
PetscErrorCode
Metos3DGeometryFinal(Metos3D *metos3d)
{
    // work vars
    char        geometryType[PETSC_MAX_PATH_LEN];    
    PetscBool   flag;
    PetscFunctionBegin;
//    // debug start
//    PetscGetTime(&metos3d->startTime[kDebugLevel]);
    // geometry
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DGeometryType", geometryType);
    PetscStrcmp("Profile", geometryType, &flag);
    if (flag == PETSC_TRUE) {
        PetscFree(metos3d->profileStart);
        PetscFree(metos3d->profileEnd);
    }
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DGeometryFinal\n");    
    PetscFunctionReturn(0);
}
