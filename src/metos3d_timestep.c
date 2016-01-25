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

#include "metos3d_timestep.h"

#undef  kDebugLevel
#define kDebugLevel kDebugLevel2

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepInit"
PetscErrorCode
Metos3DTimeStepInit(Metos3D *metos3d)
{
    PetscFunctionBegin;
    // register event
    PetscLogEventRegister("TimeStepPhi", 0, &metos3d->eventTimeStepPhi);
    // options
    Metos3DUtilOptionsGetScalar(metos3d, "-Metos3DTimeStep", &metos3d->timeStep);
    Metos3DUtilOptionsGetScalar(metos3d, "-Metos3DTimeStepStart", &metos3d->timeStepStart);
    Metos3DUtilOptionsGetInt(metos3d, "-Metos3DTimeStepCount", &metos3d->timeStepCount);
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTimeStepInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepFinal"
PetscErrorCode
Metos3DTimeStepFinal(Metos3D *metos3d)
{
    PetscFunctionBegin;
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTimeStepFinal\n");
    PetscFunctionReturn(0);
}

#undef  kDebugLevel
#define kDebugLevel kDebugLevel3

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
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTimeStepFunction\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepPhi"
PetscErrorCode
Metos3DTimeStepPhi(Metos3D *metos3d, Vec *yin, Vec *yout, PetscInt nparam, PetscReal *u0)
{
    // bgc
    PetscInt    npref       = metos3d->fileFormatPrefixCount;
    PetscInt    ntracer     = metos3d->tracerCount;
    // time step
    PetscScalar t0          = metos3d->timeStepStart;
    PetscScalar dt          = metos3d->timeStep;
    PetscInt    nstep       = metos3d->timeStepCount;
    // work vars
    PetscScalar tj;
    PetscInt    itracer, istep;
    Vec         *ywork;
    PetscFunctionBegin;
    // wait for all processors
    PetscBarrier(PETSC_NULL);
    // start log event
    PetscLogEventBegin(metos3d->eventTimeStepPhi, 0, 0, 0, 0);
    // prepare work vector
    VecDuplicateVecs(*yin, ntracer, &ywork);
    // initial point in time, project time to [0,1[
    tj = fmod(t0, 1.0);
    // init bgc, yout not set (yet)
    Metos3DBGCStepInit(metos3d, tj, dt, yin, yout, nparam, u0);
    // time step loop
    for (istep = 0; istep < nstep; istep++) {
        // point in time
        tj = fmod(t0 + istep*dt, 1.0);
        // work vars
        char filePrefixFormat[PETSC_MAX_PATH_LEN];    
        char filePrefix      [PETSC_MAX_PATH_LEN];    
        // file prefix
        if (npref > 1) {
            if ((metos3d->spinupStep + 1)%metos3d->moduloStep[0] == 0) {
                PetscInt modstep = metos3d->moduloStepCount;
                if (modstep > 0) {
                    PetscInt imodstep = metos3d->moduloStep[1];
                    if ((istep+1)%imodstep == 0) {
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
        }
        // yout = Phi(yi)
        // yin = yout
        Metos3DTimeStepPhiStep(metos3d, tj, dt, istep, yin, yout, ywork, nparam, u0);
        for(itracer = 0; itracer < ntracer; itracer++) VecCopy(yout[itracer], yin[itracer]);
    }
    // final bgc, yout not set (yet)
    Metos3DBGCStepFinal(metos3d, tj, dt, yin, yout, nparam, u0);
    // free work vector
    VecDestroyVecs(ntracer, &ywork);
    // wait for all processors
    PetscBarrier(PETSC_NULL);
    // stop log event
    PetscLogEventEnd(metos3d->eventTimeStepPhi, 0, 0, 0, 0);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DTimeStepPhi\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DTimeStepPhiStep"
PetscErrorCode
Metos3DTimeStepPhiStep(Metos3D *metos3d, PetscScalar t, PetscScalar dt, PetscInt istep, Vec *yin, Vec *yout, Vec *ywork, PetscInt nparam, PetscReal *u0)
{
    // bgc
    PetscInt    ntracer     = metos3d->tracerCount;
    // transport
    PetscInt    nmat        = metos3d->matrixCount;
    Mat         *Ae         = metos3d->matrixExplicitArray;
    Mat         *Ai         = metos3d->matrixImplicitArray;
    Mat         *Aework     = &metos3d->Aework;
    Mat         *Aiwork     = &metos3d->Aiwork;
    // work vars
    PetscInt    itracer;
    PetscFunctionBegin;
    // bgc
    Metos3DBGCStep(metos3d, t, dt, yin, yout, nparam, u0);
    // transport
    // ywork = Ae*yin
    // ywork = ywork + yout
    // yout  = Ai*ywork
    Metos3DTransport(metos3d, t, nmat, Ae, ntracer, yin, ywork, Aework);
    for (itracer = 0; itracer < ntracer; itracer++) VecAXPY(ywork[itracer], 1.0, yout[itracer]);
    Metos3DTransport(metos3d, t, nmat, Ai, ntracer, ywork, yout, Aiwork);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, FSSDSE, "Metos3DTimeStepPhiStep", "istep:", istep, "t:", t);
    PetscFunctionReturn(0);
}
