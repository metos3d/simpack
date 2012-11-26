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
 *  metos3d_util.c
 *
 * This file incorporates work covered by the following copyright and
 * permission notice:
 *
 *     Copyright (C) 2012, Samar Khatiwala, spk@ldeo.columbia.edu
 *
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or
 *     sell copies of the Software, and to permit persons to whom
 *     the Software is furnished to do so, subject to the following
 *     conditions:
 *
 *     The above copyright notice and this permission notice shall
 *     be included in all copies or substantial portions of the Software.
 *
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
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
//    // debug start
//    PetscGetTime(&metos3d->startTime[kDebugLevel]);
    PetscViewerBinaryOpen(comm, filePath, FILE_MODE_READ, &viewer);
    VecLoad(*v, viewer);
    PetscViewerDestroy(&viewer);
    // debug stop
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
//    // debug start
//    PetscGetTime(&metos3d->startTime[kDebugLevel]);
    // open file, write vector
    PetscViewerBinaryOpen(comm, filePath, FILE_MODE_WRITE, &viewer);
    VecView(*v, viewer);
    PetscViewerDestroy(&viewer);
    // debug stop
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
//    // debug start
//    PetscGetTime(&metos3d->startTime[kDebugLevel]);
    // PETSc does not read in '%', we use '$' instead.
    // IMPORTANT LIMITATION: This means you cannot use '$' in your filename.
    where = NULL;
    PetscStrchr(string, '$', &where);
    while (where != NULL) {
        *where = '%';
        PetscStrchr(string, '$', &where);            
    }
    // debug stop
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
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, F3S, "Metos3DUtilMatrixLoad", "filePath:", filePath);
    PetscFunctionReturn(0);
}

/*
 *  Metos3DUtilMatrixLoadAndCreate
 *  
 *  Creates and loads a PETSc matrix for a given distribution.
 *  This routine is based on 'MatLoadIntoMatrix2' from
 *  Samar Khatiwala's (Columbia University, spk@ldeo.columbia.edu)
 *  Transport Matrix Method package:
 *
 *  http://www.ldeo.columbia.edu/~spk/Research/TMM/MIT_Matrix_Global_2.8deg.tar.gz
 *
 */
#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilMatrixLoadAndCreate"
PetscErrorCode
Metos3DUtilMatrixLoadAndCreate(Metos3D *metos3d, char *filePath, Mat *A)
{
    MPI_Comm        comm    = metos3d->comm;
    // system
    PetscInt        n_vec   = metos3d->vectorLength;
    PetscInt        n_proc  = metos3d->processCount;
    PetscInt        myproc  = metos3d->processMine;
    // system load
    PetscInt        n_vec_loc   = metos3d->vectorLengthLocal;
    PetscInt        n_vec_prev  = metos3d->vectorLengthPrevious;
    // work
    PetscInt        i, ii, fd1, fd2, header[4], nnz, n_col_max, nnz_loc, nnz_prev;
    off_t           off, offset;
    PetscInt        *n_col_per_row_loc, *nnz_all, *jj;
    PetscScalar     *aa;
    PetscFunctionBegin;
//    // debug start
//    PetscGetTime(&metos3d->startTime[kDebugLevel]);
    // fd1 open, read header
    PetscBinaryOpen(filePath, FILE_MODE_READ, &fd1);
    PetscBinaryRead(fd1, (void*)header, 4, PETSC_INT);
    nnz = header[3];
    // fd1 seek to local part and read
    off = PETSC_BINARY_INT_SIZE*(4+n_vec_prev);
    PetscBinarySeek(fd1, off, PETSC_BINARY_SEEK_SET, &offset);
    PetscMalloc(n_vec_loc*sizeof(PetscInt), &n_col_per_row_loc);
    PetscBinaryRead(fd1, (void*)n_col_per_row_loc, n_vec_loc, PETSC_INT);
    // determine maximum memory needed
    n_col_max = 0;
    nnz_loc   = 0;
    for (i=0; i<n_vec_loc; i++)
    {
        n_col_max = max(n_col_max, n_col_per_row_loc[i]);
        nnz_loc   = nnz_loc + n_col_per_row_loc[i];
    }
    // compute seek
    PetscMalloc(n_proc*sizeof(PetscInt), &nnz_all);
    MPI_Allgather(&nnz_loc, 1, MPIU_INT, nnz_all, 1, MPIU_INT, comm);
    nnz_prev = 0;
    for (i=0; i<myproc; i++)
    {
        nnz_prev = nnz_prev + nnz_all[i];
    }
    // fd1 seek
    off = PETSC_BINARY_INT_SIZE * (4 + n_vec + nnz_prev);
    PetscBinarySeek(fd1, off, PETSC_BINARY_SEEK_SET, &offset);
    // fd2 open and seek
    PetscBinaryOpen(filePath, FILE_MODE_READ, &fd2);
    off = PETSC_BINARY_INT_SIZE * (4 + n_vec + nnz) + PETSC_BINARY_SCALAR_SIZE * nnz_prev;
    PetscBinarySeek(fd2, off, PETSC_BINARY_SEEK_SET, &offset);
    // create matrix    
    // set type to MATAIJ, compatible with sequential (MATSEQAIJ) and parallel (MATMPIAIJ)
    MatCreate(comm, A);
    MatSetSizes(*A, n_vec_loc, n_vec_loc, n_vec, n_vec);
    MatSetType(*A, MATAIJ);
    // preallocate memory
    MatSeqAIJSetPreallocation(*A, n_col_max, PETSC_NULL);
    MatMPIAIJSetPreallocation(*A, n_col_max, PETSC_NULL, n_col_max, PETSC_NULL);
    // read and set matrix
    PetscMalloc(n_col_max*sizeof(PetscInt), &jj);
    PetscMalloc(n_col_max*sizeof(PetscScalar), &aa);
    for (i=0; i<n_vec_loc; i++)
    {
        PetscBinaryRead(fd1, jj, n_col_per_row_loc[i], PETSC_INT);
        PetscBinaryRead(fd2, aa, n_col_per_row_loc[i], PETSC_SCALAR);
        ii = n_vec_prev + i;
        MatSetValues(*A, 1, &ii, n_col_per_row_loc[i], jj, aa, INSERT_VALUES);
    }
    // assembly
    MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
    // free memory
    PetscFree(aa);
    PetscFree(jj);
    PetscFree(n_col_per_row_loc);
    PetscFree(nnz_all);
    // close files
    PetscBinaryClose(fd1);
    PetscBinaryClose(fd2);
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, F3S, "Metos3DUtilMatrixLoadAndCreate", "filePath:", filePath);
    PetscFunctionReturn( 0 );
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DUtilOptionsGetInt"
PetscErrorCode
Metos3DUtilOptionsGetInt(Metos3D *metos3d, const char *optionName, PetscInt *ivalue)
{
    PetscBool   flag = PETSC_FALSE;
    char        message[PETSC_MAX_PATH_LEN];
    PetscFunctionBegin;
//    // debug start
//    PetscGetTime(&metos3d->startTime[kDebugLevel]);
    PetscOptionsGetInt(PETSC_NULL, optionName, ivalue, &flag);
    sprintf(message, "Please provide the '%s' option", optionName);
    Metos3DFlag(flag, message);
    // debug stop
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
//    // debug start
//    PetscGetTime(&metos3d->startTime[kDebugLevel]);
    PetscOptionsGetScalar(PETSC_NULL, optionName, dvalue, &flag);
    sprintf(message, "Please provide the '%s' option", optionName);
    Metos3DFlag(flag, message);    
    // debug stop
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
//    // debug start
//    PetscGetTime(&metos3d->startTime[kDebugLevel]);
    PetscOptionsGetRealArray(PETSC_NULL, optionName, dvalue, nmax, &flag);
    sprintf(message, "Please provide the '%s' option", optionName);
    Metos3DFlag(flag, message);
    for (i=0; i<(*nmax); i++)
    {
        Metos3DDebug(metos3d, kDebugLevel, F4SE, "Metos3DUtilOptionsGetRealArray", "optionName:", optionName, "value:", dvalue[i]);
    }
    // debug stop
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
//    // debug start
//    PetscGetTime(&metos3d->startTime[kDebugLevel]);
    PetscOptionsGetString(PETSC_NULL, optionName, string, PETSC_MAX_PATH_LEN, &flag);
    sprintf(message, "Please provide the '%s' option", optionName);
    Metos3DFlag(flag, message);    
    // debug stop
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
    //    // debug start
    //    PetscGetTime(&metos3d->startTime[kDebugLevel]);
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
    // debug stop
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
    //    // debug start
    //    PetscGetTime(&metos3d->startTime[kDebugLevel]);
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
    // debug stop
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
    //    // debug start
    //    PetscGetTime(&metos3d->startTime[kDebugLevel]);
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
    // debug stop
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
    //    // debug start
    //    PetscGetTime(&metos3d->startTime[kDebugLevel]);
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
    // debug stop
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
    //    // debug start
    //    PetscGetTime(&metos3d->startTime[kDebugLevel]);
    //
    PetscMalloc(ntracer*sizeof(Vec), v);
    for (itracer = 0; itracer < ntracer; itracer++)
    {
        VecCreateMPI(comm, nvecloc, nvec, &(*v)[itracer]);
        VecSet((*v)[itracer], val[itracer]);
        VecAssemblyBegin((*v)[itracer]);
        VecAssemblyEnd((*v)[itracer]);
    }
    // debug stop
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
//        PetscScalar w;
//        w = t * N + 0.5;
//        *beta    = fmod(w, 1.0);
//        *i_beta  = (PetscInt)fmod(floor(w), (PetscScalar)N);
//        *alpha   = (1.0 - *beta);
//        *i_alpha = (PetscInt)fmod(floor(w) + N - 1, (PetscScalar)N);
        // work vars
        PetscInt    i_work;
        i_work      = (((int)floor(t*2.0*N))+1)/2;
        *i_alpha    = (i_work-1+N)%N;
        *i_beta     = i_work%N;
        *beta       = fmod(N*t+0.5, 1.0);
        *alpha      = 1.0-(*beta);
    } else {
        *i_alpha    = 0;
        *i_beta     = 0;
        *beta       = 0.0;
        *alpha      = 1.0;        
    }
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, FS2SESD, "Metos3DUtilInterpolate", "alpha:", *alpha, "i_alpha:", *i_alpha, "beta:", *beta, "i_beta:", *i_beta);
    PetscFunctionReturn(0);
}

