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

#include "metos3d_transport.h"

#undef  kDebugLevel
#define kDebugLevel kDebugLevel2

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTransportInit"
PetscErrorCode
Metos3DTransportInit(Metos3D *metos3d)
{
    // work vars
    char        transportType[PETSC_MAX_PATH_LEN];    
    PetscBool   flag;
    PetscFunctionBegin;
    // transport type
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DTransportType", transportType);
    PetscStrcmp("Matrix", transportType, &flag);
    if (flag == PETSC_TRUE) {
        // matrix
        Metos3DTransportMatrixInit(metos3d);
    }
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTransportInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTransportFinal"
PetscErrorCode
Metos3DTransportFinal(Metos3D *metos3d)
{
    // work vars
    char        transportType[PETSC_MAX_PATH_LEN];    
    PetscBool   flag;
    PetscFunctionBegin;
    // transport type
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DTransportType", transportType);
    PetscStrcmp("Matrix", transportType, &flag);
    if (flag == PETSC_TRUE) {
        // matrix
        Metos3DTransportMatrixFinal(metos3d);
    }
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTransportFinal\n");
    PetscFunctionReturn(0);
}

#undef  kDebugLevel
#define kDebugLevel kDebugLevel3

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTransportMatrixInit"
PetscErrorCode
Metos3DTransportMatrixInit(Metos3D *metos3d)
{
    // work vars
    PetscInt    n_vec     = metos3d->vectorLength;
    PetscInt    n_vec_loc = metos3d->vectorLengthLocal;
    char        inputDirectory  [PETSC_MAX_PATH_LEN];
    char        formatexp       [PETSC_MAX_PATH_LEN];
    char        formatimp       [PETSC_MAX_PATH_LEN];
    PetscInt    nmat, imat;
    Mat         *Ae, *Ai;
    char        fileFormat[PETSC_MAX_PATH_LEN];
    char        filePath  [PETSC_MAX_PATH_LEN];
    PetscFunctionBegin;
    // input directory
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DMatrixInputDirectory", inputDirectory);
    // matrixCount
    Metos3DUtilOptionsGetInt(metos3d, "-Metos3DMatrixCount", &nmat);
    metos3d->matrixCount = nmat;
    // format
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DMatrixExplicitFileFormat", formatexp);
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DMatrixImplicitFileFormat", formatimp);
    // create mat
    PetscMalloc(nmat*sizeof(Mat), &Ae);
    PetscMalloc(nmat*sizeof(Mat), &Ai);
    metos3d->matrixExplicitArray = Ae;
    metos3d->matrixImplicitArray = Ai;
    // parse format
    Metos3DUtilFormatParse(metos3d, formatexp);
    Metos3DUtilFormatParse(metos3d, formatimp);
    // load all explicit mats
    for (imat = 0; imat < nmat; imat++)
    {
        // file path
        sprintf(fileFormat, "%s%s", "%s", formatexp);
        sprintf(filePath, fileFormat, inputDirectory, imat);
        // create matrix, set size and type
        MatCreate(metos3d->comm, &Ae[imat]);
        MatSetSizes(Ae[imat], n_vec_loc, n_vec_loc, n_vec, n_vec);
        MatSetType(Ae[imat], MATAIJ);
        // load matrix
        Metos3DUtilMatrixLoad(metos3d, filePath, &Ae[imat]);
    }
    // work matrix
    MatDuplicate(Ae[0], MAT_DO_NOT_COPY_VALUES, &metos3d->Aework);
    // load all implicit mats
    for (imat = 0; imat < nmat; imat++)
    {
        sprintf(fileFormat, "%s%s", "%s", formatimp);
        sprintf(filePath, fileFormat, inputDirectory, imat);        
        // create matrix, set size and type
        MatCreate(metos3d->comm, &Ai[imat]);
        MatSetSizes(Ai[imat], n_vec_loc, n_vec_loc, n_vec, n_vec);
        MatSetType(Ai[imat], MATAIJ);
        // load matrix
        Metos3DUtilMatrixLoad(metos3d, filePath, &Ai[imat]);
    }
    // work matrix
    MatDuplicate(Ai[0], MAT_DO_NOT_COPY_VALUES, &metos3d->Aiwork);
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTransportMatrixInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTransportMatrixFinal"
PetscErrorCode
Metos3DTransportMatrixFinal(Metos3D *metos3d)
{
    // work vars
    PetscInt    nmat    = metos3d->matrixCount;
    PetscInt    imat;
    PetscFunctionBegin;
    // matrix
    for (imat = 0; imat < nmat; imat++)
    {
        MatDestroy(&metos3d->matrixExplicitArray[imat]);
        MatDestroy(&metos3d->matrixImplicitArray[imat]);
    }
    PetscFree(metos3d->matrixExplicitArray);
    PetscFree(metos3d->matrixImplicitArray);
    // work
    MatDestroy(&metos3d->Aework);
    MatDestroy(&metos3d->Aiwork);
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTransportMatrixFinal\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTransport"
PetscErrorCode
Metos3DTransport(Metos3D *metos3d, PetscScalar t, PetscInt n_mat, Mat *A, PetscInt n_tracer, Vec *y_in, Vec *y_out, Mat *A_work)
{
    // work vars
    PetscInt    i, i_alpha, i_beta;
    PetscScalar alpha, beta;
    PetscFunctionBegin;
    // interpolate
    Metos3DUtilInterpolate(metos3d, t, n_mat, &i_alpha, &alpha, &i_beta, &beta);
    // A_work = A[i_alpha]
    // A_work = alpha * A_work
    // A_work = A_work + beta * A[i_beta]
    MatCopy(A[i_alpha], *A_work, SAME_NONZERO_PATTERN);
    MatScale(*A_work, alpha);
    MatAXPY(*A_work, beta, A[i_beta], SAME_NONZERO_PATTERN);
    // apply matrix
    for (i=0; i<n_tracer; i++) {
        MatMult(*A_work, y_in[i], y_out[i]);
    }
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTransport\n");
    PetscFunctionReturn(0);
}
