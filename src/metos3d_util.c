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

#include "metos3d_util.h"

#undef  kDebugLevel
#define kDebugLevel kDebugLevel3

#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilVectorLoad"
PetscErrorCode
Metos3DUtilVectorLoad(Metos3D *metos3d, char *filePath, Vec *v)
{
    MPI_Comm comm = metos3d->comm;
    PetscViewer viewer;
    PetscFunctionBegin;
    PetscViewerBinaryOpen(comm, filePath, FILE_MODE_READ, &viewer);
    VecLoad(*v, viewer);
    PetscViewerDestroy(&viewer);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, F3S, "Metos3DUtilVectorLoad", "filePath:", filePath);
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilVectorView"
PetscErrorCode
Metos3DUtilVectorView(Metos3D *metos3d, char *filePath, Vec *v)
{
    MPI_Comm comm = metos3d->comm;
    // work vars
    PetscViewer     viewer;
    PetscFunctionBegin;
    // open file, write vector
    PetscViewerBinaryOpen(comm, filePath, FILE_MODE_WRITE, &viewer);
    VecView(*v, viewer);
    PetscViewerDestroy(&viewer);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, F3S, "Metos3DUtilVectorView", "filePath:", filePath);
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilFormatParse"
PetscErrorCode
Metos3DUtilFormatParse(Metos3D *metos3d, char *string)
{
    char *where;
    PetscFunctionBegin;
    // PETSc does not read in '%', we use '$' instead.
    // IMPORTANT LIMITATION: This means you cannot use '$' in your filename.
    where = NULL;
    PetscStrchr(string, '$', &where);
    while (where != NULL) {
        *where = '%';
        PetscStrchr(string, '$', &where);            
    }
    // debug
    Metos3DDebug(metos3d, kDebugLevel, F3S, "Metos3DUtilFormatParse", "format:", string);
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilMatrixLoad"
PetscErrorCode
Metos3DUtilMatrixLoad(Metos3D *metos3d, char *filePath, Mat *A)
{
    MPI_Comm comm = metos3d->comm;
    PetscViewer viewer;
    PetscFunctionBegin;
    PetscViewerBinaryOpen(comm, filePath, FILE_MODE_READ, &viewer);
    MatLoad(*A, viewer);
    PetscViewerDestroy(&viewer);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, F3S, "Metos3DUtilMatrixLoad", "filePath:", filePath);
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilOptionsGetInt"
PetscErrorCode
Metos3DUtilOptionsGetInt(Metos3D *metos3d, const char *optionName, PetscInt *ivalue)
{
    PetscBool   flag = PETSC_FALSE;
    char        message[PETSC_MAX_PATH_LEN];
    PetscFunctionBegin;
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, optionName, ivalue, &flag);
    sprintf(message, "Please provide the '%s' option", optionName);
    Metos3DFlag(flag, message);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, F4SD, "Metos3DUtilOptionsGetInt", "optionName:", optionName, "value:", *ivalue);
    PetscFunctionReturn(0);
}
#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilOptionsGetScalar"
PetscErrorCode
Metos3DUtilOptionsGetScalar(Metos3D *metos3d, const char *optionName, PetscScalar *dvalue)
{
    PetscBool   flag = PETSC_FALSE;
    char        message[PETSC_MAX_PATH_LEN];
    PetscFunctionBegin;
    PetscOptionsGetScalar(PETSC_NULL, PETSC_NULL, optionName, dvalue, &flag);
    sprintf(message, "Please provide the '%s' option", optionName);
    Metos3DFlag(flag, message);    
    // debug
    Metos3DDebug(metos3d, kDebugLevel, F4SE, "Metos3DUtilOptionsGetScalar", "optionName:", optionName, "value:", *dvalue);
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilOptionsGetRealArray"
PetscErrorCode
Metos3DUtilOptionsGetRealArray(Metos3D *metos3d, const char *optionName, PetscInt *nmax, PetscReal *dvalue)
{
    PetscBool   flag = PETSC_FALSE;
    char        message[PETSC_MAX_PATH_LEN];
    PetscInt    i;
    PetscFunctionBegin;
    PetscOptionsGetRealArray(PETSC_NULL, PETSC_NULL, optionName, dvalue, nmax, &flag);
    sprintf(message, "Please provide the '%s' option", optionName);
    Metos3DFlag(flag, message);
    for (i=0; i<(*nmax); i++)
    {
        Metos3DDebug(metos3d, kDebugLevel, F4SE, "Metos3DUtilOptionsGetRealArray", "optionName:", optionName, "value:", dvalue[i]);
    }
    // debug
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilOptionsGetString"
PetscErrorCode
Metos3DUtilOptionsGetString(Metos3D *metos3d, const char *optionName, char *string)
{
    PetscBool   flag = PETSC_FALSE;
    char        message[PETSC_MAX_PATH_LEN];
    PetscFunctionBegin;
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, optionName, string, PETSC_MAX_PATH_LEN, &flag);
    sprintf(message, "Please provide the '%s' option", optionName);
    Metos3DFlag(flag, message);    
    // debug
    Metos3DDebug(metos3d, kDebugLevel, F5S, "Metos3DUtilOptionsGetString", "optionName:", optionName, "value:", string);
    PetscFunctionReturn(0);
}

#undef  kDebugLevel
#define kDebugLevel kDebugLevel4

#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilVecCopySeparateToDiagonal"
PetscErrorCode
Metos3DUtilVecCopySeparateToDiagonal(Metos3D *metos3d, PetscInt n, PetscInt nvecloc, Vec *y, Vec *yBD)
{
    // geometry
    PetscInt        nprofloc    = metos3d->profileCountLocal;
    PetscInt        *istartloc  = metos3d->profileStartLocal;
    PetscInt        *iendloc    = metos3d->profileEndLocal;
    // work vars
    PetscInt        i, iprof, ii, length;
    PetscScalar     **yarray;
    PetscScalar     *yarrayBD;
    PetscFunctionBegin;
    // get arrays
    VecGetArrays(y, n, &yarray);
    VecGetArray(*yBD, &yarrayBD);
    for (i = 0; i < n; i++) {
        for (iprof = 0; iprof < nprofloc; iprof++) {
            // length
            length = iendloc[iprof]-istartloc[iprof]+1;
            // start, end are 1 indexed
            for (ii = 0; ii < length; ii++) {
                yarrayBD[i*length + ii + n*(istartloc[iprof]-1)] = yarray[i][ii + istartloc[iprof]-1];
            }
        }
    }
    // restore arrays
    VecRestoreArrays(y, n, &yarray);
    VecRestoreArray(*yBD, &yarrayBD);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DUtilVecCopySeparateToDiagonal\n");
    PetscFunctionReturn( 0 );
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilVecCopyDiagonalToSeparate"
PetscErrorCode
Metos3DUtilVecCopyDiagonalToSeparate(Metos3D *metos3d, PetscInt n, PetscInt nvecloc, Vec *yBD, Vec *y)
{
    // geometry
    PetscInt        nprofloc    = metos3d->profileCountLocal;
    PetscInt        *istartloc  = metos3d->profileStartLocal;
    PetscInt        *iendloc    = metos3d->profileEndLocal;
    // work vars
    PetscInt        i, iprof, ii, length;
    PetscScalar     **yarray;
    PetscScalar     *yarrayBD;
    PetscFunctionBegin;
    // get arrays
    VecGetArray(*yBD, &yarrayBD);
    VecGetArrays(y, n, &yarray);
    for (i = 0; i < n; i++) {
        for (iprof = 0; iprof < nprofloc; iprof++) {
            // length
            length = iendloc[iprof]-istartloc[iprof]+1;
            // start, end are 1 indexed
            for (ii = 0; ii < length; ii++) {
                yarray[i][ii + istartloc[iprof]-1] = yarrayBD[i*length + ii + n*(istartloc[iprof]-1)];
            }
        }
    }
    // restore arrays
    VecRestoreArray(*yBD, &yarrayBD);
    VecRestoreArrays(y, n, &yarray);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DUtilVecCopyDiagonalToSeparate\n");
    PetscFunctionReturn( 0 );
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilVecCopySeparateToDiagonalProfile"
PetscErrorCode
Metos3DUtilVecCopySeparateToDiagonalProfile(Metos3D *metos3d, PetscInt n, PetscInt nloc, Vec *y, Vec *yBD)
{
    // work vars
    PetscInt        i, iloc;
    PetscScalar     **yarray;
    PetscScalar     *yarrayBD;
    PetscFunctionBegin;
    // get arrays
    VecGetArrays(y, n, &yarray);
    VecGetArray(*yBD, &yarrayBD);
    for (i = 0; i < n; i++) {
        for (iloc = 0; iloc < nloc; iloc++ ) {
            yarrayBD[i + n*iloc] = yarray[i][iloc];
        }
    }
    VecRestoreArrays(y, n, &yarray);
    VecRestoreArray(*yBD, &yarrayBD);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DUtilVecCopySeparateToDiagonalProfile\n");
    PetscFunctionReturn( 0 );
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilVecCreateAndSetValue"
PetscErrorCode
Metos3DUtilVecCreateAndSetValue(Metos3D *metos3d, PetscInt ntracer, PetscInt nvec, PetscInt nvecloc, Vec **v, PetscScalar val)
{
    MPI_Comm    comm = metos3d->comm;
    // work vars
    PetscInt    itracer;
    PetscFunctionBegin;
    // create array
    PetscMalloc(ntracer*sizeof(Vec), v);
    for (itracer = 0; itracer < ntracer; itracer++)
    {
        // create vector
        VecCreateMPI(comm, nvecloc, nvec, &(*v)[itracer]);
        VecSet((*v)[itracer], val);
        VecAssemblyBegin((*v)[itracer]);
        VecAssemblyEnd  ((*v)[itracer]);
    }
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DUtilVecCreateAndSetValue\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilVecCreateAndSetValues"
PetscErrorCode
Metos3DUtilVecCreateAndSetValues(Metos3D *metos3d, PetscInt ntracer, PetscInt nvec, PetscInt nvecloc, Vec **v, PetscScalar *val)
{
    MPI_Comm    comm = metos3d->comm;
    // work vars
    PetscInt    itracer;
    PetscFunctionBegin;
    PetscMalloc(ntracer*sizeof(Vec), v);
    for (itracer = 0; itracer < ntracer; itracer++)
    {
        VecCreateMPI(comm, nvecloc, nvec, &(*v)[itracer]);
        VecSet((*v)[itracer], val[itracer]);
        VecAssemblyBegin((*v)[itracer]);
        VecAssemblyEnd((*v)[itracer]);
    }
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DUtilVecCreateAndSetValues\n");
    PetscFunctionReturn(0);
}    

#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilInterpolate"
PetscErrorCode
Metos3DUtilInterpolate(Metos3D *metos3d, PetscReal t, PetscInt N, PetscInt *i_alpha, PetscScalar *alpha, PetscInt *i_beta, PetscScalar *beta)
{
    PetscFunctionBegin;
    // we assume t is (t mod 1.0)!
    if (N > 1) {
        PetscScalar w;
        w = t * N + 0.5;
        *beta    = fmod(w, 1.0);
        *i_beta  = (PetscInt)fmod(floor(w), (PetscScalar)N);
        *alpha   = (1.0 - *beta);
        *i_alpha = (PetscInt)fmod(floor(w) + N - 1, (PetscScalar)N);
    } else {
        *i_alpha    = 0;
        *i_beta     = 0;
        *beta       = 0.0;
        *alpha      = 1.0;        
    }
    // debug
    Metos3DDebug(metos3d, kDebugLevel, FS2SESD, "Metos3DUtilInterpolate", "alpha:", *alpha, "i_alpha:", *i_alpha, "beta:", *beta, "i_beta:", *i_beta);
    PetscFunctionReturn(0);
}

