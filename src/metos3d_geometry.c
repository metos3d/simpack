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
    char        geometryInputDirectory  [PETSC_MAX_PATH_LEN];
    char        profileFile             [PETSC_MAX_PATH_LEN];
    char        volumeFile              [PETSC_MAX_PATH_LEN];
    char        filePath                [PETSC_MAX_PATH_LEN];
    MatInfo     info;
    Vec         rowSum, rowMax;
    PetscScalar sum, max;
    PetscInt    profileCount;
    PetscInt    vectorLength;
    PetscInt    nlayermax;
    PetscFunctionBegin;
    // read file path and names
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DGeometryInputDirectory", geometryInputDirectory);
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DProfileFile", profileFile);
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DVolumeFile", volumeFile);
    // create and read profile matrix
    sprintf(filePath, "%s%s", geometryInputDirectory, profileFile);
    MatCreate(comm, &metos3d->Profiles);
    Metos3DUtilMatrixLoad(metos3d, filePath, &metos3d->Profiles);
    // create and read relative volumes vector
    sprintf(filePath, "%s%s", geometryInputDirectory, volumeFile);
    VecCreate(comm, &metos3d->rVolumes);
    Metos3DUtilVectorLoad(metos3d, filePath, &metos3d->rVolumes);
    // number of profiles
    MatGetInfo(metos3d->Profiles, MAT_GLOBAL_SUM, &info);
    profileCount = info.nz_allocated;
    metos3d->profileCount = profileCount;
    Metos3DDebug(metos3d, kDebugLevel+1, F2SD, "Metos3DGeometryInit", "profileCount:", profileCount);
    // vector length
    MatGetVecs(metos3d->Profiles, PETSC_NULL, &rowSum);
    MatGetRowSum(metos3d->Profiles, rowSum);
    VecSum(rowSum, &sum);
    vectorLength = (int)sum;
    metos3d->vectorLength = vectorLength;
    Metos3DDebug(metos3d, kDebugLevel+1, F2SD, "Metos3DGeometryInit", "vectorLength:", vectorLength);
    VecDestroy(&rowSum);
    // max profile length
    MatGetVecs(metos3d->Profiles, PETSC_NULL, &rowMax);
    MatGetRowMaxAbs(metos3d->Profiles, rowMax, PETSC_NULL);
    VecMax(rowMax, PETSC_NULL, &max);
    nlayermax = (int)max;
    metos3d->profileLengthMax = nlayermax;
    Metos3DDebug(metos3d, kDebugLevel+1, F2SD, "Metos3DGeometryInit", "profileLengthMax:", nlayermax);
    VecDestroy(&rowMax);    
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DGeometryInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DGeometryFinal"
PetscErrorCode
Metos3DGeometryFinal(Metos3D *metos3d)
{
    PetscFunctionBegin;
    VecDestroy(&metos3d->rVolumes);
    MatDestroy(&metos3d->Profiles);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DGeometryFinal\n");    
    PetscFunctionReturn(0);
}
