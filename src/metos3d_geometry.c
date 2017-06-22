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
 */

#include "metos3d_geometry.h"

#undef  kDebugLevel
#define kDebugLevel kDebugLevel2

#undef  __FUNCT__
#define __FUNCT__ "Metos3DGeometryInit"
PetscErrorCode
Metos3DGeometryInit(Metos3D *metos3d)
{
    // work vars
    char        geometryType[PETSC_MAX_PATH_LEN];
    PetscBool   flag;
    PetscFunctionBegin;
    // geometry type
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DGeometryType", geometryType);
    PetscStrcmp("Profile", geometryType, &flag);
    if (flag == PETSC_TRUE) {
        // work vars
        char                geometryInputDirectory  [PETSC_MAX_PATH_LEN];
        char                maskFile                [PETSC_MAX_PATH_LEN];
        char                volumeFile              [PETSC_MAX_PATH_LEN];
        char                filePath                [PETSC_MAX_PATH_LEN];
        Mat                 landSeaMask;
        PetscInt            nrow, ncols, irow, icol, ival, idxoff, valoff;
        const PetscInt      *cols;
        const PetscScalar   *vals;
        PetscViewer         viewer;
        PetscInt            profileCount, nlayermax, vectorLength;
        
        // options
        Metos3DUtilOptionsGetString(metos3d, "-Metos3DProfileInputDirectory", geometryInputDirectory);
        Metos3DUtilOptionsGetString(metos3d, "-Metos3DProfileMaskFile", maskFile);
        
        // profiles from land-sea mask
        sprintf(filePath, "%s%s", geometryInputDirectory, maskFile);
        Metos3DDebug(metos3d, kDebugLevel, F3S, "Metos3DGeometryInit", "filePath:", filePath);
        // create mat on every proc
        MatCreate(PETSC_COMM_SELF, &landSeaMask);
        PetscViewerBinaryOpen(PETSC_COMM_SELF, filePath, FILE_MODE_READ, &viewer);
        MatLoad(landSeaMask, viewer);
        PetscViewerDestroy(&viewer);
        // get land-sea mask dimensions
        MatGetSize(landSeaMask, &nrow, NULL);
        // go through rows and collect profileCount
        profileCount = 0;
        for (irow = 0; irow < nrow; irow++) {
            MatGetRow(landSeaMask, irow, &ncols, NULL, NULL);
            profileCount += ncols;
            MatRestoreRow(landSeaMask, irow, &ncols, NULL, NULL);
        }
        // allocate memory for indices
        PetscMalloc(profileCount * sizeof(PetscInt), &metos3d->profileStart);
        PetscMalloc(profileCount * sizeof(PetscInt), &metos3d->profileEnd);
        // go through rows, collect incides, determine max
        nlayermax = 0;
        idxoff = 0;
        valoff = 0;
        for (irow = 0; irow < nrow; irow++) {
            // get row
            MatGetRow(landSeaMask, irow, &ncols, &cols, &vals);
            // go through columns
            for (icol = 0; icol < ncols; icol++) {
                // cast to int
                ival = (PetscInt)vals[icol];
                // check max
                nlayermax = max(nlayermax, ival);
                // compute starts and ends
                metos3d->profileStart[idxoff] = valoff + 1;
                metos3d->profileEnd[idxoff] = valoff + ival;
                // increase offsets
                valoff += ival;
                idxoff++;
            }
            // restore row
            MatRestoreRow(landSeaMask, irow, &ncols, &cols, &vals);
        }
        // free land sea mask
        MatDestroy(&landSeaMask);
        // set vectors length as last end index
        vectorLength = metos3d->profileEnd[profileCount - 1];
        // store
        // profile count
        metos3d->profileCount = profileCount;
        Metos3DDebug(metos3d, kDebugLevel, F2SD, "Metos3DGeometryInit", "profileCount:", profileCount);
        // max length of profile
        metos3d->profileLengthMax = nlayermax;
        Metos3DDebug(metos3d, kDebugLevel, F2SD, "Metos3DGeometryInit", "profileLengthMax:", nlayermax);
        // vectorLength
        metos3d->vectorLength = vectorLength;
        Metos3DDebug(metos3d, kDebugLevel, F2SD, "Metos3DGeometryInit", "vectorLength:", vectorLength);
        // volumes
        // read in option
        Metos3DUtilOptionsGetString(metos3d, "-Metos3DProfileVolumeFile", volumeFile);
        // create file path
        sprintf(filePath, "%s%s", geometryInputDirectory, volumeFile);
        Metos3DDebug(metos3d, kDebugLevel, F3S, "Metos3DGeometryInit", "filePath:", filePath);
        // create seq vector and read in data
        VecCreateSeq(PETSC_COMM_SELF, vectorLength, &metos3d->volumes);
        PetscViewerBinaryOpen(PETSC_COMM_SELF, filePath, FILE_MODE_READ, &viewer);
        VecLoad(metos3d->volumes, viewer);
        PetscViewerDestroy(&viewer);
        // compute square root of volumes
        VecSqrtAbs(metos3d->volumes);
    }
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
    // geometry type
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DGeometryType", geometryType);
    PetscStrcmp("Profile", geometryType, &flag);
    if (flag == PETSC_TRUE) {
        PetscFree(metos3d->profileStart);
        PetscFree(metos3d->profileEnd);
        VecDestroy(&metos3d->volumes);
    }
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DGeometryFinal\n");
    PetscFunctionReturn(0);    
}
