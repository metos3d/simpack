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
 *  metos3d_model.c
 *
 */

#include "metos3d_timestep.h"

#pragma mark -
#undef  kDebugLevel
#define kDebugLevel kDebugLevel2

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepInit"
PetscErrorCode
Metos3DTimeStepInit(Metos3D *metos3d)
{
    PetscFunctionBegin;
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTimeStepInit\n");
    // options
    Metos3DUtilOptionsGetScalar(metos3d, "-Metos3DTimeStep", &metos3d->timeStep);
    Metos3DUtilOptionsGetScalar(metos3d, "-Metos3DTimeStepStart", &metos3d->timeStepStart);
    Metos3DUtilOptionsGetInt(metos3d, "-Metos3DTimeStepCount", &metos3d->timeStepCount);
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepFinal"
PetscErrorCode
Metos3DTimeStepFinal(Metos3D *metos3d)
{
    PetscFunctionBegin;
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTimeStepFinal\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepFunction"
PetscErrorCode
Metos3DTimeStepFunction(SNES snes, Vec ynBD, Vec fnBD, void *ctx)
{
    Metos3D     *metos3d    = (Metos3D*)ctx;
    // geometry
    PetscInt    nvec        = metos3d->vectorLength;
    // bgc
    PetscInt    ntracer     = metos3d->tracerCount;
    // load
    PetscInt    nvecloc     = metos3d->vectorLengthLocal;
    // parameter
    PetscInt    nparam      = metos3d->parameterCount;
    PetscReal   *u0         = metos3d->u0;
    // work vars
    PetscInt    itracer;
    Vec         *yin, *yinold, *yout;
    PetscFunctionBegin;
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTimeStepFunction\n");
    // wait for all processors  
    PetscBarrier(PETSC_NULL);
    // create work vectors
    Metos3DUtilVecCreateAndSetValue(metos3d, ntracer, nvec, nvecloc, &yin, 0.0);
    Metos3DUtilVecCreateAndSetValue(metos3d, ntracer, nvec, nvecloc, &yinold, 0.0);
    Metos3DUtilVecCreateAndSetValue(metos3d, ntracer, nvec, nvecloc, &yout, 0.0);
    // yin      = ynBD
    // yinold   = ynBD
    // yout     = Phi(yin)
    // yout     = - yout + yinold;
    // fnBD     = yout
    Metos3DUtilVecCopyDiagonalToSeparate(metos3d, ntracer, nvecloc, &ynBD, yin);
    Metos3DUtilVecCopyDiagonalToSeparate(metos3d, ntracer, nvecloc, &ynBD, yinold);
    Metos3DTimeStepPhi(metos3d, yin, yout, nparam, u0);
    for (itracer = 0; itracer < ntracer; itracer++) VecAYPX(yout[itracer], -1.0, yinold[itracer]);
    Metos3DUtilVecCopySeparateToDiagonal(metos3d, ntracer, nvecloc, yout, &fnBD);
    // free work vectors
    VecDestroyVecs(ntracer, &yin);
    VecDestroyVecs(ntracer, &yinold);
    VecDestroyVecs(ntracer, &yout);
    // wait for all processors  
    PetscBarrier(PETSC_NULL);
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepPhi"
PetscErrorCode
Metos3DTimeStepPhi(Metos3D *metos3d, Vec *yin, Vec *yout, PetscInt nparam, PetscReal *u0)
{
    PetscInt    npref   = metos3d->fileFormatPrefixCount;
    // load
    PetscInt    nvecloc     = metos3d->vectorLengthLocal;
    // bgc
    PetscInt    ntracer     = metos3d->tracerCount;
    Vec         *ybgcinBD   = metos3d->ybgcinBD;
    Vec         *ybgcoutBD  = metos3d->ybgcoutBD;
    // time step
    PetscScalar t0          = metos3d->timeStepStart;
    PetscScalar dt          = metos3d->timeStep;
    PetscInt    stepmax     = metos3d->timeStepCount;    
    // work vars
    PetscScalar t;
    PetscInt    itracer, istep;
    Vec         *ywork;
    PetscFunctionBegin;
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTimeStepPhi\n");
    // prepare work vector
    VecDuplicateVecs(*yin, ntracer, &ywork);
    // init
    istep = 0;
    // time step
    // we project all to a period of 1.0
    t = fmod(t0, 1.0);
    // bgc
    Metos3DUtilVecCopySeparateToDiagonal(metos3d, ntracer, nvecloc, yin, ybgcinBD);
    VecZeroEntries(*ybgcoutBD);
    Metos3DBGCStepInit(metos3d, t, dt, ybgcinBD, ybgcoutBD, nparam, u0);
    // step
    while(istep < stepmax)
    {
        // work vars
        char filePrefixFormat[PETSC_MAX_PATH_LEN];    
        char filePrefix      [PETSC_MAX_PATH_LEN];    
        // yout = Phi(yi)
        // yin = yout
        Metos3DTimeStepPhiStep(metos3d, t, dt, istep, yin, yout, ywork, nparam, u0);
        for(itracer = 0; itracer < ntracer; itracer++) VecCopy(yout[itracer], yin[itracer]);
        
        // file prefix
        if (npref > 1) {
            PetscInt modstep = metos3d->moduloStepCount;
            if (modstep > 0) {
                PetscInt imodstep = metos3d->moduloStep[1];
                if (istep%imodstep == 0) {
                    sprintf(filePrefixFormat, "%s%s", metos3d->filePrefix, metos3d->fileFormatPrefix[1]);
                    sprintf(filePrefix, filePrefixFormat, istep);
                    // output
                    Metos3DBGCOutputPrefix(metos3d, filePrefix, ntracer, yin);
                }
            } else {
                sprintf(filePrefixFormat, "%s%s", metos3d->filePrefix, metos3d->fileFormatPrefix[1]);
                sprintf(filePrefix, filePrefixFormat, istep);
                // output
                Metos3DBGCOutputPrefix(metos3d, filePrefix, ntracer, yin);
            }
        }
        
        // step forward
        istep++;
        t = fmod(t0 + istep*dt, 1.0);
    }
    // yin = yout
    for(itracer = 0; itracer < ntracer; itracer++) VecCopy(yin[itracer], yout[itracer]);
    // final    
    Metos3DUtilVecCopySeparateToDiagonal(metos3d, ntracer, nvecloc, yin, ybgcinBD);
    VecZeroEntries(*ybgcoutBD);
    Metos3DBGCStepFinal(metos3d, t, dt, ybgcinBD, ybgcoutBD, nparam, u0);
    // free work vector
    VecDestroyVecs(ntracer, &ywork);
    PetscFunctionReturn(0);
}

#pragma mark -
#undef  kDebugLevel
#define kDebugLevel kDebugLevel3

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepPhiStep"
PetscErrorCode
Metos3DTimeStepPhiStep(Metos3D *metos3d, PetscScalar t, PetscScalar dt, PetscInt istep, Vec *yin, Vec *yout, Vec *ywork, PetscInt nparam, PetscReal *u0)
{
    // bgc
    PetscInt    ntracer     = metos3d->tracerCount;
    Vec         *ybgcinBD   = metos3d->ybgcinBD;
    Vec         *ybgcoutBD  = metos3d->ybgcoutBD;
    // load
    PetscInt    nvecloc     = metos3d->vectorLengthLocal;
    // transport
    PetscInt    nmat        = metos3d->matrixCount;
    Mat         *Ae         = metos3d->matrixExplicitArray;
    Mat         *Ai         = metos3d->matrixImplicitArray;
    Mat         *Aework     = &metos3d->Aework;
    Mat         *Aiwork     = &metos3d->Aiwork;
    // work vars
    PetscInt    itracer;
    PetscFunctionBegin;
    Metos3DDebug(metos3d, kDebugLevel, FSSDSE, "Metos3DTimeStepPhiStep", "istep:", istep, "t:", t);
    // bgc
    // ybgcinBD  = yin
    // ybgcoutBD = 0
    // ybgcoutBD = bgc(ybgcinBD)
    // yout      = ybgcoutBD
    Metos3DUtilVecCopySeparateToDiagonal(metos3d, ntracer, nvecloc, yin, ybgcinBD);
    VecZeroEntries(*ybgcoutBD);
    Metos3DBGCStep(metos3d, t, dt, ybgcinBD, ybgcoutBD, nparam, u0);
    Metos3DUtilVecCopyDiagonalToSeparate(metos3d, ntracer, nvecloc, ybgcoutBD, yout);
    // transport
    // ywork = Ae*yin
    // ywork = ywork + yout
    // yout  = Ai*ywork
    Metos3DTransport(metos3d, t, nmat, Ae, ntracer, yin, ywork, Aework);
    for (itracer = 0; itracer < ntracer; itracer++) VecAXPY(ywork[itracer], 1.0, yout[itracer]);
    Metos3DTransport(metos3d, t, nmat, Ai, ntracer, ywork, yout, Aiwork);
    PetscFunctionReturn(0);
}
