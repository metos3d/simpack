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

#include "metos3d_solver.h"

#undef  kDebugLevel
#define kDebugLevel kDebugLevel1

#undef  __FUNCT__
#define __FUNCT__ "Metos3DSolver"
PetscErrorCode
Metos3DSolver(Metos3D *metos3d)
{
    // bgc
    PetscInt    ntracer = metos3d->tracerCount;
    PetscInt    nvec    = metos3d->vectorLength;
    PetscInt    nvecloc = metos3d->vectorLengthLocal;
    // work vars
    char        solverType[PETSC_MAX_PATH_LEN];
    char        message   [PETSC_MAX_PATH_LEN];
    PetscBool   flag;
    PetscFunctionBegin;
    // solver type
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DSolverType", solverType);
    PetscStrcmp("Newton", solverType, &flag);
    if (flag == PETSC_TRUE) {
        // work vars
        Vec *ystarBD;
        // init to zero before start
        metos3d->fileFormatPrefixCount = 0;
        metos3d->moduloStepCount = 0;
        // SNES
        // ystarBD = y0
        Metos3DUtilVecCreateAndSetValue(metos3d, 1, ntracer*nvec, ntracer*nvecloc, &ystarBD, 0.0);
        Metos3DUtilVecCopySeparateToDiagonal(metos3d, ntracer, nvecloc, metos3d->y0, ystarBD);
        // solve
        SNESSolve(metos3d->snes, PETSC_NULL, *ystarBD);
        // y0 = ystarBD
        Metos3DUtilVecCopyDiagonalToSeparate(metos3d, ntracer, nvecloc, ystarBD, metos3d->y0);
        VecDestroyVecs(1, &ystarBD);
        // output
        Metos3DBGCOutput(metos3d, ntracer, metos3d->y0);
    } else {
        PetscStrcmp("Spinup", solverType, &flag);
        if (flag == PETSC_TRUE) {
            // monitor file format
            PetscInt    filemax, ifile;
            PetscInt    stepmax;
            // init to zero before start
            metos3d->fileFormatPrefixCount = 0;
            metos3d->moduloStepCount = 0;
            filemax = 2;
            PetscOptionsGetStringArray(PETSC_NULL, PETSC_NULL, "-Metos3DSpinupMonitorFileFormatPrefix", metos3d->fileFormatPrefix, &filemax, &flag);
            metos3d->fileFormatPrefixCount = filemax;
            if (filemax > 0) {
                for (ifile = 0; ifile < filemax; ifile++) {
                    Metos3DUtilFormatParse(metos3d, metos3d->fileFormatPrefix[ifile]);
                }
                PetscMalloc(PETSC_MAX_PATH_LEN * sizeof(char), &metos3d->filePrefix);
            }
            
            // monitor modulo step
            if (filemax > 0) {
                stepmax = 2;
                PetscMalloc(stepmax*sizeof(PetscInt), &metos3d->moduloStep);
                PetscOptionsGetIntArray(PETSC_NULL, PETSC_NULL, "-Metos3DSpinupMonitorModuloStep", metos3d->moduloStep, &stepmax, &flag);
                metos3d->moduloStepCount = stepmax;
            }
            // solver
            Metos3DSolverSpinup(metos3d);
            
            // file format
            if (filemax > 0) {
                PetscFree(metos3d->moduloStep);
                PetscFree(metos3d->filePrefix);
                for (ifile = 0; ifile < filemax; ifile++) {
                    PetscFree(metos3d->fileFormatPrefix[ifile]);
                }
            }
        } else {
            sprintf(message, "Unknown: -Metos3DSolverType %s", solverType);
            Metos3DFlag(flag, message);
        }
    }
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DSolver\n");
    PetscFunctionReturn(0);
}

#undef  kDebugLevel
#define kDebugLevel kDebugLevel2

#undef  __FUNCT__
#define __FUNCT__ "Metos3DSolverInit"
PetscErrorCode
Metos3DSolverInit(Metos3D *metos3d)
{
    MPI_Comm    comm        = metos3d->comm;
    // bgc
    PetscInt    ntracer     = metos3d->tracerCount;
    // geometry
    PetscInt    nvec        = metos3d->vectorLength;
    // load
    PetscInt    nvecloc     = metos3d->vectorLengthLocal;
    // work vars
    char        solverType[PETSC_MAX_PATH_LEN];    
    char        message   [PETSC_MAX_PATH_LEN];
    PetscInt    itracer;
    PetscBool   flag;
    PetscFunctionBegin;
    // solver type
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DSolverType", solverType);
    PetscStrcmp("Newton", solverType, &flag);
    if (flag == PETSC_TRUE) {
        // SNES
        // create vectors
        Metos3DUtilVecCreateAndSetValue(metos3d, 1, ntracer*nvec, ntracer*nvecloc, &metos3d->fnBD, 0.0);
        // snes
        SNESCreate(comm, &metos3d->snes);
        // function
        SNESSetFunction(metos3d->snes, *metos3d->fnBD, Metos3DTimeStepFunction, (void*)metos3d);
        // mat shell
        MatCreateSNESMF(metos3d->snes, &metos3d->JShell);
        SNESSetJacobian(metos3d->snes, metos3d->JShell, metos3d->JShell, Metos3DSolverFormJacobian, (void*)metos3d);
        // option
        SNESSetOptionsPrefix(metos3d->snes, "Metos3DNewton_");
        SNESSetFromOptions(metos3d->snes);
    } else {
        PetscStrcmp("Spinup", solverType, &flag);
        if (flag == PETSC_TRUE) {
            // Spinup
            // prepare norm weights
            // create MPI vectors
            Metos3DUtilVecCreateAndSetValue(metos3d, ntracer, nvec, nvecloc, &metos3d->normWeights, 0.0);
            Metos3DUtilVecCreateAndSetValue(metos3d, 1, ntracer*nvec, ntracer*nvecloc, &metos3d->normWeightsBD, 0.0);
            // create scatter context
            IS is;
            VecScatter normWeightsScatter;
            ISCreateStride(PETSC_COMM_SELF, nvec, 0, 1, &is);
            VecScatterCreate(metos3d->volumes, is, metos3d->normWeights[0], is, &normWeightsScatter);
            // loop over tracers
            // copy from volumes to ntracer times normWeights
            for (itracer = 0; itracer < ntracer; itracer++) {
                VecScatterBegin(normWeightsScatter, metos3d->volumes, metos3d->normWeights[itracer], INSERT_VALUES, SCATTER_FORWARD);
                VecScatterEnd(normWeightsScatter, metos3d->volumes, metos3d->normWeights[itracer], INSERT_VALUES, SCATTER_FORWARD);
            }
            // destroy scatter context
            ISDestroy(&is);
            VecScatterDestroy(&normWeightsScatter);
            // copy from separate to block diagonal
            Metos3DUtilVecCopySeparateToDiagonal(metos3d, ntracer, nvecloc, metos3d->normWeights, metos3d->normWeightsBD);
        } else {
            sprintf(message, "Unknown: -Metos3DSolverType %s", solverType);
            Metos3DFlag(flag, message);
        }
    }
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DSolverInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DSolverFinal"
PetscErrorCode
Metos3DSolverFinal(Metos3D *metos3d)
{
    // work vars
    char        solverType[PETSC_MAX_PATH_LEN];    
    char        message   [PETSC_MAX_PATH_LEN];    
    PetscBool   flag;
    PetscFunctionBegin;
    // solver type
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DSolverType", solverType);
    PetscStrcmp("Newton", solverType, &flag);
    if (flag == PETSC_TRUE) {
        // SNES
        VecDestroyVecs(1, &metos3d->fnBD);
        MatDestroy(&metos3d->JShell);
        SNESDestroy(&metos3d->snes);
    } else {
        PetscStrcmp("Spinup", solverType, &flag);
        if (flag == PETSC_TRUE) {
            PetscInt ntracer = metos3d->tracerCount;
            // Spinup
            // norm weights
            VecDestroyVecs(1, &metos3d->normWeightsBD);
            VecDestroyVecs(ntracer, &metos3d->normWeights);
            
        } else {
            sprintf(message, "Unknown: -Metos3DSolverType %s", solverType);
            Metos3DFlag(flag, message);
        }
    }
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DSolverFinal\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DSolverSpinup"
PetscErrorCode
Metos3DSolverSpinup(Metos3D *metos3d)
{
    PetscInt    npref   = metos3d->fileFormatPrefixCount;
    // bgc
    PetscInt    ntracer = metos3d->tracerCount;
    PetscInt    nvec    = metos3d->vectorLength;
    PetscInt    nvecloc = metos3d->vectorLengthLocal;
    // parameter
    PetscInt    nparam  = metos3d->parameterCount;
    PetscReal   *u0     = metos3d->u0;
    // solver
    Vec         *normWeightsBD = metos3d->normWeightsBD;
    // work vars
    PetscInt    nstep, istep, itracer;
    PetscScalar acc;
    PetscBool   flag        = PETSC_FALSE;
    PetscBool   countFlag   = PETSC_FALSE;
    PetscBool   accFlag     = PETSC_FALSE;
    PetscBool   monFlag     = PETSC_FALSE;
    Vec         *yin, *yout;
    Vec         *ystarBD, *yworkBD;
    PetscReal   ynorm, ynormweight;
    PetscFunctionBegin;
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-Metos3DSpinupCount", &nstep, &flag);
    if (flag == PETSC_TRUE) {
        // count
        countFlag = PETSC_TRUE;
        Metos3DDebug(metos3d, kDebugLevel3, F4SD, "Metos3DSolverSpinup", "optionName:", "-Metos3DSpinupCount", "value:", nstep);
    }
    PetscOptionsGetScalar(PETSC_NULL, PETSC_NULL, "-Metos3DSpinupTolerance", &acc, &flag);
    if (flag == PETSC_TRUE) {
        // acc
        accFlag = PETSC_TRUE;
        Metos3DDebug(metos3d, kDebugLevel3, F4SE, "Metos3DSolverSpinup", "optionName:", "-Metos3DSpinupTolerance", "value:", acc);
    }
    if ((countFlag == PETSC_FALSE) && (accFlag == PETSC_FALSE)) {
        char message[PETSC_MAX_PATH_LEN];    
        sprintf(message, "Please provide the '%s' or '%s' option", "-Metos3DSpinupCount", "-Metos3DSpinupTolerance");
        Metos3DFlag(PETSC_FALSE, message);
    }
    // create vectors
    Metos3DUtilVecCreateAndSetValue(metos3d, ntracer, nvec, nvecloc, &yin, 0.0);
    Metos3DUtilVecCreateAndSetValue(metos3d, ntracer, nvec, nvecloc, &yout, 0.0);    
    Metos3DUtilVecCreateAndSetValue(metos3d, 1, ntracer*nvec, ntracer*nvecloc, &ystarBD, 0.0);    
    Metos3DUtilVecCreateAndSetValue(metos3d, 1, ntracer*nvec, ntracer*nvecloc, &yworkBD, 0.0);    
    // init
    Metos3DUtilVecCopySeparateToDiagonal(metos3d, ntracer, nvecloc, metos3d->y0, ystarBD);
    Metos3DUtilVecCopyDiagonalToSeparate(metos3d, ntracer, nvecloc, ystarBD, yin);
    // monitor
    PetscOptionsGetBool(PETSC_NULL, PETSC_NULL, "-Metos3DSpinupMonitor", &monFlag, &flag);
    
    // bgc init
    PetscScalar t0 = metos3d->timeStepStart;
    PetscScalar dt = metos3d->timeStep;
    PetscScalar tj;
    tj = fmod(t0, 1.0);
    Metos3DBGCStepInit(metos3d, tj, dt, yin, yout, nparam, u0);
    
    // spinup
    istep = 0;
    ynorm = 1.e300;
    while (1) {
        // check conditions
        if ((countFlag == PETSC_FALSE) && (accFlag == PETSC_TRUE)) {
            if (ynorm < acc) break;
        }
        if ((countFlag == PETSC_TRUE) && (accFlag == PETSC_FALSE)) {
            if (istep >= nstep) break;
        }
        if ((countFlag == PETSC_TRUE) && (accFlag == PETSC_TRUE)) {
            if ((ynorm < acc) || (istep >= nstep)) break;
        }
        // file prefix
        if (npref > 0) {
            sprintf(metos3d->filePrefix, metos3d->fileFormatPrefix[0], istep);
        }
        // monitor
        if (monFlag) PetscTime(&metos3d->startTime[kDebugLevel0]);

        // set value for output control
        metos3d->spinupStep = istep;
        // phi step
        // yin = yout
        Metos3DTimeStepPhi(metos3d, yin, yout, nparam, u0);
        for(itracer = 0; itracer < ntracer; itracer++) {
            VecCopy(yout[itracer], yin[itracer]);
        }
        
        // monitor
        if (monFlag) {
            // yworkBD = yin
            // ystarBD = ystarBD - yworkBD
            // ynorm
            // ystarBD = ystarBD .* normWeights
            // ynormweight
            // ystarBD = yworkBD
            Metos3DUtilVecCopySeparateToDiagonal(metos3d, ntracer, nvecloc, yin, yworkBD);
            VecAXPY(*ystarBD, -1.0, *yworkBD);
            VecNorm(*ystarBD, NORM_2, &ynorm);
            VecPointwiseMult(*ystarBD, *normWeightsBD, *ystarBD);
            VecNorm(*ystarBD, NORM_2, &ynormweight);
            VecCopy(*yworkBD, *ystarBD);
            // print out
            Metos3DDebug(metos3d, kDebugLevel0, FDSEE, istep, "Spinup Function norm", ynorm, ynormweight);
        }
        
        // step counter
        istep++;
    }
    
    // bgc final
    Metos3DBGCStepFinal(metos3d, tj, dt, yin, yout, nparam, u0);

    // output
    Metos3DBGCOutput(metos3d, ntracer, yin);
    // free work vectors
    VecDestroyVecs(ntracer, &yin);
    VecDestroyVecs(ntracer, &yout);
    VecDestroyVecs(1, &ystarBD);
    VecDestroyVecs(1, &yworkBD);
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DSolverSpinup\n");

    PetscFunctionReturn(0);
}

#undef  kDebugLevel
#define kDebugLevel kDebugLevel4

#undef  __FUNCT__
#define __FUNCT__ "Metos3DSolverFormJacobian"
PetscErrorCode
Metos3DSolverFormJacobian(SNES snes, Vec ynBD, Mat JShell, Mat PCShell, void *ctx)
{
    Metos3D     *metos3d = (Metos3D*)ctx;
    PetscFunctionBegin;
    // assemble matrix after MatSetBase from SNES
    MatAssemblyBegin(JShell, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(JShell, MAT_FINAL_ASSEMBLY);
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DSolverFormJacobian\n");
    PetscFunctionReturn(0);
}
