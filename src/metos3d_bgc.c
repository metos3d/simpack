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

#include "metos3d_bgc.h"

#undef  kDebugLevel
#define kDebugLevel kDebugLevel2

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCInit"
PetscErrorCode
Metos3DBGCInit(Metos3D *metos3d)
{
    PetscFunctionBegin;
    // register event
    PetscLogEventRegister("BGCStep", 0, &metos3d->eventBGCStep);
    // init tracer, boundary, domain, parameter
    Metos3DBGCTracerInit(metos3d);
    
    Metos3DBGCDiagnosticInit(metos3d);
//    Option left: name:-Metos3DTracerDiagnosticCount value: 1
//    Option left: name:-Metos3DTracerDiagnosticName value: fP
//    Option left: name:-Metos3DTracerDiagnosticOutputDirectory value: work/
    
    Metos3DBGCParameterInit(metos3d);
    Metos3DBGCBoundaryConditionInit(metos3d);
    Metos3DBGCDomainConditionInit(metos3d);
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCFinal"
PetscErrorCode
Metos3DBGCFinal(Metos3D *metos3d)
{
    PetscFunctionBegin;
    // final parameter, domain, boundary, tracer
    Metos3DBGCDomainConditionFinal(metos3d);
    Metos3DBGCBoundaryConditionFinal(metos3d);
    Metos3DBGCParameterFinal(metos3d);
    
    Metos3DBGCDiagnosticFinal(metos3d);

    Metos3DBGCTracerFinal(metos3d);
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCFinal\n");
    PetscFunctionReturn(0);
}

#undef  kDebugLevel
#define kDebugLevel kDebugLevel3



#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCDiagnosticInit"
PetscErrorCode
Metos3DBGCDiagnosticInit(Metos3D *metos3d)
{
    
    // work vars
    PetscInt    ndiag, idiag;
    PetscInt    nmax;
    PetscBool   flag;
    char        *diagFileNames[PETSC_MAX_PATH_LEN];
    
    // diag count
    Metos3DUtilOptionsGetInt(metos3d, "-Metos3DTracerDiagnosticCount", &ndiag);
    metos3d->diagCount = ndiag;
    // diag name

    // file name
    nmax = ndiag;
    PetscOptionsGetStringArray(PETSC_NULL, PETSC_NULL, "-Metos3DTracerDiagnosticName", diagFileNames, &nmax, &flag);
    if (flag == PETSC_TRUE) {
        for(idiag = 0; idiag < ndiag; idiag++) {
            Metos3DDebug(metos3d, kDebugLevel, F3S, "Metos3DTracerDiagnosticName", "filename:", diagFileNames[idiag]);
        }
    }
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCDiagnosticInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCDiagnosticFinal"
PetscErrorCode
Metos3DBGCDiagnosticFinal(Metos3D *metos3d)
{
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCDiagnosticFinal\n");
    PetscFunctionReturn(0);
}



#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCTracerInit"
PetscErrorCode
Metos3DBGCTracerInit(Metos3D *metos3d)
{
    // load
    PetscInt    nvec    = metos3d->vectorLength;
    PetscInt    nvecloc = metos3d->vectorLengthLocal;
    // work vars
    PetscInt    ntracer, itracer;
    PetscInt    nmax;
    PetscBool   flag;
    char        initFileNameFormat  [PETSC_MAX_PATH_LEN];
    PetscFunctionBegin;
    // tracer count
    Metos3DUtilOptionsGetInt(metos3d, "-Metos3DTracerCount", &ntracer);
    metos3d->tracerCount = ntracer;
    // read tracer init options
    // order:
    // 1. file name format
    // 2. file names
    // 3. constant values
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-Metos3DTracerInitFileFormat", initFileNameFormat, PETSC_MAX_PATH_LEN, &flag);
    if (flag == PETSC_TRUE) {
        // work vars
        char    tracerInputDirectory[PETSC_MAX_PATH_LEN];
        char    format              [PETSC_MAX_PATH_LEN];
        char    filePath            [PETSC_MAX_PATH_LEN];
        // input directory, vector, format
        Metos3DUtilOptionsGetString(metos3d, "-Metos3DTracerInputDirectory", tracerInputDirectory);
        Metos3DUtilVecCreateAndSetValue(metos3d, ntracer, nvec, nvecloc, &metos3d->y0, 0.0);
        Metos3DUtilFormatParse(metos3d, initFileNameFormat);
        for (itracer = 0; itracer < ntracer; itracer++)
        {
            sprintf(format, "%s%s", tracerInputDirectory, initFileNameFormat);
            sprintf(filePath, format, itracer);
            Metos3DUtilVectorLoad(metos3d, filePath, &metos3d->y0[itracer]);
        }
    } else {
        // work vars
        char    tracerInputDirectory[PETSC_MAX_PATH_LEN];
        char    *initFileNames      [PETSC_MAX_PATH_LEN];
        char    filePath            [PETSC_MAX_PATH_LEN];
        // file name
        nmax = ntracer;
        PetscOptionsGetStringArray(PETSC_NULL, PETSC_NULL, "-Metos3DTracerInitFile", initFileNames, &nmax, &flag);
        if (flag == PETSC_TRUE)
        {
            // input directory, vector, file name
            Metos3DUtilOptionsGetString(metos3d, "-Metos3DTracerInputDirectory", tracerInputDirectory);
            Metos3DUtilVecCreateAndSetValue(metos3d, ntracer, nvec, nvecloc, &metos3d->y0, 0.0);
            for (itracer = 0; itracer < ntracer; itracer++)
            {
                sprintf(filePath, "%s%s", tracerInputDirectory, initFileNames[itracer]);                
                Metos3DUtilVectorLoad(metos3d, filePath, &metos3d->y0[itracer]);
                PetscFree(initFileNames[itracer]);
            }
        }
        else
        {
            // work vars
            PetscReal   *y0value;
            // value
            nmax = ntracer;
            PetscMalloc(nmax*sizeof(PetscReal),&y0value);
            // array, vector
            Metos3DUtilOptionsGetRealArray(metos3d, "-Metos3DTracerInitValue", &nmax, y0value);
            Metos3DUtilVecCreateAndSetValues(metos3d, ntracer, nvec, nvecloc, &metos3d->y0, y0value);
            PetscFree(y0value);
        }
    }
    // create work vectors
    Metos3DUtilVecCreateAndSetValue(metos3d, 1, ntracer*nvec, ntracer*nvecloc, &metos3d->ybgcinBD, 0.0);
    Metos3DUtilVecCreateAndSetValue(metos3d, 1, ntracer*nvec, ntracer*nvecloc, &metos3d->ybgcoutBD, 0.0);
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCTracerInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCParameterInit"
PetscErrorCode
Metos3DBGCParameterInit(Metos3D *metos3d)
{
    // work vars
    PetscInt    nparam;
    PetscFunctionBegin;
    // parameter count
    Metos3DUtilOptionsGetInt(metos3d, "-Metos3DParameterCount", &nparam);
    metos3d->parameterCount = nparam;
    if (nparam > 0) {
        // work vars
        PetscInt    nmax;
        // real array
        nmax = nparam;
        PetscMalloc(nmax*sizeof(PetscReal),&metos3d->u0);
        Metos3DUtilOptionsGetRealArray(metos3d, "-Metos3DParameterValue", &nmax, metos3d->u0);
    }
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCParameterInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCBoundaryConditionInit"
PetscErrorCode
Metos3DBGCBoundaryConditionInit(Metos3D *metos3d)
{
    // work vars
    PetscInt    nbc;
    PetscFunctionBegin;
    // count
    Metos3DUtilOptionsGetInt(metos3d, "-Metos3DBoundaryConditionCount", &nbc);
    metos3d->boundaryConditionCount = nbc;
    if (nbc > 0) {
        // geometry
        PetscInt    nprof       = metos3d->profileCount;
        // load
        PetscInt    nprofloc    = metos3d->profileCountLocal;
        // work vars
        PetscInt    *nbcs;
        char        inputDirectory  [PETSC_MAX_PATH_LEN];
        char        *conditionName  [PETSC_MAX_PATH_LEN];
        char        message         [PETSC_MAX_PATH_LEN];
        char        optionName      [PETSC_MAX_PATH_LEN];
        char        fileFormat      [PETSC_MAX_PATH_LEN];
        char        filePath        [PETSC_MAX_PATH_LEN];
        char        fileName        [PETSC_MAX_PATH_LEN];
        PetscInt    nmax, ibc, ibcs;
        PetscBool   flag;
        // count array
        PetscMalloc(nbc*sizeof(PetscInt),&nbcs);
        metos3d->boundaryConditionCountArray = nbcs;
        // input directory
        Metos3DUtilOptionsGetString(metos3d, "-Metos3DBoundaryConditionInputDirectory", inputDirectory);
        // name
        nmax = nbc;
        PetscOptionsGetStringArray(PETSC_NULL, PETSC_NULL, "-Metos3DBoundaryConditionName", conditionName, &nmax, &flag);
        sprintf(message, "Please provide the '%s' option", "-Metos3DBoundaryConditionName");
        Metos3DFlag(flag, message);
        // 
        PetscMalloc(nbc*sizeof(Vec*), &metos3d->boundaryConditionArray);
        // condition
        for (ibc = 0; ibc < nbc; ibc++) {
            Metos3DDebug(metos3d, kDebugLevel, F3S, "Metos3DBGCBoundaryConditionInit", "name:", conditionName[ibc]);
            // count
            sprintf(optionName, "-Metos3D%sCount", conditionName[ibc]);
            Metos3DUtilOptionsGetInt(metos3d, optionName, &nbcs[ibc]);
            // name
            sprintf(optionName, "-Metos3D%sFileFormat", conditionName[ibc]);
            Metos3DUtilOptionsGetString(metos3d, optionName, fileFormat);
            Metos3DUtilFormatParse(metos3d, fileFormat);
            // create vector
            Metos3DUtilVecCreateAndSetValue(metos3d, nbcs[ibc], nprof, nprofloc, &metos3d->boundaryConditionArray[ibc], 0.0);
            // files
            for (ibcs = 0; ibcs < nbcs[ibc]; ibcs++) {
                // file path
                sprintf(fileName, fileFormat, ibcs);
                sprintf(filePath, "%s%s", inputDirectory, fileName);
                // load
                Metos3DUtilVectorLoad(metos3d, filePath, &metos3d->boundaryConditionArray[ibc][ibcs]);
            }
            // name
            PetscFree(conditionName[ibc]);
        }
        // bgc api
        Metos3DUtilVecCreateAndSetValue(metos3d, nbc, nprof, nprofloc, &metos3d->bgcbc, 0.0);
        Metos3DUtilVecCreateAndSetValue(metos3d, 1, nbc*nprof, nbc*nprofloc, &metos3d->bgcbcBD, 0.0);
    }
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCBoundaryConditionInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCDomainConditionInit"
PetscErrorCode
Metos3DBGCDomainConditionInit(Metos3D *metos3d)
{
    // work vars
    PetscInt    ndc;
    PetscFunctionBegin;
    // count
    Metos3DUtilOptionsGetInt(metos3d, "-Metos3DDomainConditionCount", &ndc);
    metos3d->domainConditionCount = ndc;
    if (ndc > 0) {
        // geometry
        PetscInt    nvec    = metos3d->vectorLength;
        // load
        PetscInt    nvecloc = metos3d->vectorLengthLocal;
        // work vars
        PetscInt    *ndcs;
        char        inputDirectory  [PETSC_MAX_PATH_LEN];
        char        *conditionName  [PETSC_MAX_PATH_LEN];
        char        message         [PETSC_MAX_PATH_LEN];
        char        optionName      [PETSC_MAX_PATH_LEN];
        char        fileFormat      [PETSC_MAX_PATH_LEN];
        char        filePath        [PETSC_MAX_PATH_LEN];
        char        fileName        [PETSC_MAX_PATH_LEN];
        PetscInt    nmax, idc, idcs;
        PetscBool   flag;
        // count array
        PetscMalloc(ndc*sizeof(PetscInt),&ndcs);
        metos3d->domainConditionCountArray = ndcs;
        // input directory
        Metos3DUtilOptionsGetString(metos3d, "-Metos3DDomainConditionInputDirectory", inputDirectory);
        // name
        nmax = ndc;
        PetscOptionsGetStringArray(PETSC_NULL, PETSC_NULL, "-Metos3DDomainConditionName", conditionName, &nmax, &flag);
        sprintf(message, "Please provide the '%s' option", "-Metos3DDomainConditionName");
        Metos3DFlag(flag, message);
        // 
        PetscMalloc(ndc*sizeof(Vec*), &metos3d->domainConditionArray);
        // condition
        for (idc = 0; idc < ndc; idc++) {
            Metos3DDebug(metos3d, kDebugLevel, F3S, "Metos3DBGCDomainConditionInit", "name:", conditionName[idc]);
            // count
            sprintf(optionName, "-Metos3D%sCount", conditionName[idc]);
            Metos3DUtilOptionsGetInt(metos3d, optionName, &ndcs[idc]);
            // name
            sprintf(optionName, "-Metos3D%sFileFormat", conditionName[idc]);
            Metos3DUtilOptionsGetString(metos3d, optionName, fileFormat);
            Metos3DUtilFormatParse(metos3d, fileFormat);
            // create vector
            Metos3DUtilVecCreateAndSetValue(metos3d, ndcs[idc], nvec, nvecloc, &metos3d->domainConditionArray[idc], 0.0);
            // files
            for (idcs = 0; idcs < ndcs[idc]; idcs++) {
                // file path
                sprintf(fileName, fileFormat, idcs);
                sprintf(filePath, "%s%s", inputDirectory, fileName);
                // load
                Metos3DUtilVectorLoad(metos3d, filePath, &metos3d->domainConditionArray[idc][idcs]);
            }            
            // name
            PetscFree(conditionName[idc]);
        }
        // bgc api
        Metos3DUtilVecCreateAndSetValue(metos3d, ndc, nvec, nvecloc, &metos3d->bgcdc, 0.0);
        Metos3DUtilVecCreateAndSetValue(metos3d, 1, ndc*nvec, ndc*nvecloc, &metos3d->bgcdcBD, 0.0);
    }
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCDomainConditionInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCDomainConditionFinal"
PetscErrorCode
Metos3DBGCDomainConditionFinal(Metos3D *metos3d)
{
    // work vars
    PetscInt    ndc     = metos3d->domainConditionCount;
    PetscFunctionBegin;
    if (ndc > 0) {
        // work vars
        PetscInt    idc;
        PetscInt    *ndcs = metos3d->domainConditionCountArray;
        Vec         **dcs = metos3d->domainConditionArray;
        // loop
        for (idc = 0; idc < ndc; idc++) {
            // vec array
            VecDestroyVecs(ndcs[idc], &dcs[idc]);
        }
        PetscFree(dcs);
        PetscFree(ndcs);
        // bgc api
        VecDestroyVecs(ndc, &metos3d->bgcdc);
        VecDestroyVecs(1, &metos3d->bgcdcBD);
    }
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCDomainConditionFinal\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCBoundaryConditionFinal"
PetscErrorCode
Metos3DBGCBoundaryConditionFinal(Metos3D *metos3d)
{
    // work vars
    PetscInt    nbc     = metos3d->boundaryConditionCount;
    PetscFunctionBegin;
    if (nbc > 0) {
        // work vars
        PetscInt    ibc;
        PetscInt    *nbcs = metos3d->boundaryConditionCountArray;
        Vec         **bcs = metos3d->boundaryConditionArray;
        // loop
        for (ibc = 0; ibc < nbc; ibc++) {
            // vec array
            VecDestroyVecs(nbcs[ibc], &bcs[ibc]);
        }
        // count array
        // value array
        PetscFree(bcs);
        PetscFree(nbcs);
        // bgc api
        VecDestroyVecs(nbc, &metos3d->bgcbc);
        VecDestroyVecs(1, &metos3d->bgcbcBD);
    }
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCBoundaryConditionFinal\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCParameterFinal"
PetscErrorCode
Metos3DBGCParameterFinal(Metos3D *metos3d)
{
    // parameter
    PetscInt    nparam  = metos3d->parameterCount;
    PetscFunctionBegin;
    if (nparam > 0) {
        PetscFree(metos3d->u0);    
    }
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCParameterFinal\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCTracerFinal"
PetscErrorCode
Metos3DBGCTracerFinal(Metos3D *metos3d)
{
    // bgc
    PetscInt    ntracer = metos3d->tracerCount;
    PetscFunctionBegin;
    // initial value vector
    VecDestroyVecs(ntracer, &metos3d->y0);
    // work vectors
    VecDestroyVecs(1, &metos3d->ybgcinBD);
    VecDestroyVecs(1, &metos3d->ybgcoutBD);
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCTracerFinal\n");
    PetscFunctionReturn(0);
}

#undef  kDebugLevel
#define kDebugLevel kDebugLevel3

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCStepInit"
PetscErrorCode
Metos3DBGCStepInit(Metos3D *metos3d, PetscScalar t, PetscScalar dt, Vec *yin, Vec *yout, PetscInt nparam, PetscReal *u0)
{
    // load
    PetscInt    nvecloc     = metos3d->vectorLengthLocal;
    PetscInt    nprofloc    = metos3d->profileCountLocal;
    PetscInt    *istartloc  = metos3d->profileStartLocal;
    PetscInt    *iendloc    = metos3d->profileEndLocal;
    // bgc
    PetscInt    ntracer     = metos3d->tracerCount;
    PetscInt    nbc         = metos3d->boundaryConditionCount;
    PetscInt    ndc         = metos3d->domainConditionCount;
    // work vars
    Vec         *ybgcinBD   = metos3d->ybgcinBD;
    Vec         *ybgcoutBD  = metos3d->ybgcoutBD;
    PetscInt    iprof, nlayer;
    PetscScalar *bcarray    = PETSC_NULL;
    PetscScalar *dcarray    = PETSC_NULL;
    PetscScalar *yinarray   = PETSC_NULL;
    PetscScalar *youtarray  = PETSC_NULL;
    PetscFunctionBegin;
    // copy to block diagonal and zero
    Metos3DUtilVecCopySeparateToDiagonal(metos3d, ntracer, nvecloc, yin, ybgcinBD);
    VecZeroEntries(*ybgcoutBD);
    // boundary condition step
    // domain condition step
    Metos3DBGCStepBoundaryCondition(metos3d, t);    
    Metos3DBGCStepDomainCondition(metos3d, t);
    // bc, dc, yin, yout array
    if (nbc > 0) VecGetArray(*metos3d->bgcbcBD, &bcarray);
    if (ndc > 0) VecGetArray(*metos3d->bgcdcBD, &dcarray);
    VecGetArray(*ybgcinBD, &yinarray);
    VecGetArray(*ybgcoutBD, &youtarray);
    // go through profiles
    for (iprof = 0; iprof < nprofloc; iprof++)
    {
        // layer
        nlayer = iendloc[iprof]-istartloc[iprof]+1;        
        //
        // call BGC C API
        //
#ifdef BGCINIT
        BGCINIT(
                (int*)&ntracer,                                     // tracer count
                (int*)&nlayer,                                      // layer count
                (int*)&nparam,                                      // parameter count
                (int*)&nbc,                                         // boundary condition count
                (int*)&ndc,                                         // domain condition count
                (double*)&dt,                                       // time step
                (double*)&youtarray[ntracer*(istartloc[iprof]-1)],  // source minus sink
                (double*)&t,                                        // point in time
                (double*)&yinarray [ntracer*(istartloc[iprof]-1)],  // tracer
                (double*)u0,                                        // parameter
                (double*)&bcarray  [nbc*iprof],                     // boundary condition
                (double*)&dcarray  [ndc*(istartloc[iprof]-1)]       // domain condition
                );
#endif        
    }
    // bc, dc, yin, yout array
    if (nbc > 0) VecRestoreArray(*metos3d->bgcbcBD, &bcarray);
    if (ndc > 0) VecRestoreArray(*metos3d->bgcdcBD, &dcarray);
    VecRestoreArray(*ybgcinBD, &yinarray);
    VecRestoreArray(*ybgcoutBD, &youtarray);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, F2SE, "Metos3DBGCStepInit", "t:", t);
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCStepFinal"
PetscErrorCode
Metos3DBGCStepFinal(Metos3D *metos3d, PetscScalar t, PetscScalar dt, Vec *yin, Vec *yout, PetscInt nparam, PetscReal *u0)
{
    // load
    PetscInt    nvecloc     = metos3d->vectorLengthLocal;
    PetscInt    nprofloc    = metos3d->profileCountLocal;
    PetscInt    *istartloc  = metos3d->profileStartLocal;
    PetscInt    *iendloc    = metos3d->profileEndLocal;
    // bgc
    PetscInt    ntracer     = metos3d->tracerCount;
    PetscInt    nbc         = metos3d->boundaryConditionCount;
    PetscInt    ndc         = metos3d->domainConditionCount;
    // work vars
    Vec         *ybgcinBD   = metos3d->ybgcinBD;
    Vec         *ybgcoutBD  = metos3d->ybgcoutBD;
    PetscInt    iprof, nlayer;
    PetscScalar *bcarray    = PETSC_NULL;
    PetscScalar *dcarray    = PETSC_NULL;
    PetscScalar *yinarray   = PETSC_NULL;
    PetscScalar *youtarray  = PETSC_NULL;
    PetscFunctionBegin;
    // copy to block diagonal and zero
    Metos3DUtilVecCopySeparateToDiagonal(metos3d, ntracer, nvecloc, yin, ybgcinBD);
    VecZeroEntries(*ybgcoutBD);
    // boundary condition step
    // domain condition step
    Metos3DBGCStepBoundaryCondition(metos3d, t);
    Metos3DBGCStepDomainCondition(metos3d, t);
    // bc, dc, yin, yout array
    if (nbc > 0) VecGetArray(*metos3d->bgcbcBD, &bcarray);
    if (ndc > 0) VecGetArray(*metos3d->bgcdcBD, &dcarray);
    VecGetArray(*ybgcinBD, &yinarray);
    VecGetArray(*ybgcoutBD, &youtarray);
    // go through profiles
    for (iprof = 0; iprof < nprofloc; iprof++)
    {
        // layer
        nlayer = iendloc[iprof]-istartloc[iprof]+1;        
        //
        // call BGC C API
        //
#ifdef BGCFINAL
        BGCFINAL(
                (int*)&ntracer,                                     // tracer count
                (int*)&nlayer,                                      // layer count
                (int*)&nparam,                                      // parameter count
                (int*)&nbc,                                         // boundary condition count
                (int*)&ndc,                                         // domain condition count
                (double*)&dt,                                       // time step
                (double*)&youtarray[ntracer*(istartloc[iprof]-1)],  // source minus sink
                (double*)&t,                                        // point in time
                (double*)&yinarray [ntracer*(istartloc[iprof]-1)],  // tracer
                (double*)u0,                                        // parameter
                (double*)&bcarray  [nbc*iprof],                     // boundary condition
                (double*)&dcarray  [ndc*(istartloc[iprof]-1)]       // domain condition
                );
#endif        
    }
    // bc, dc, yin, yout array
    if (nbc > 0) VecRestoreArray(*metos3d->bgcbcBD, &bcarray);
    if (ndc > 0) VecRestoreArray(*metos3d->bgcdcBD, &dcarray);
    VecRestoreArray(*ybgcinBD, &yinarray);
    VecRestoreArray(*ybgcoutBD, &youtarray);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, F2SE, "Metos3DBGCStepFinal", "t:", t);
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCStep"
PetscErrorCode
Metos3DBGCStep(Metos3D *metos3d, PetscScalar t, PetscScalar dt, Vec *yin, Vec *yout, PetscInt nparam, PetscReal *u0)
{
    // load
    PetscInt    nvecloc     = metos3d->vectorLengthLocal;
    PetscInt    nprofloc    = metos3d->profileCountLocal;
    PetscInt    *istartloc  = metos3d->profileStartLocal;
    PetscInt    *iendloc    = metos3d->profileEndLocal;
    // bgc
    PetscInt    ntracer     = metos3d->tracerCount;
    PetscInt    nbc         = metos3d->boundaryConditionCount;
    PetscInt    ndc         = metos3d->domainConditionCount;
    // work vars
    Vec         *ybgcinBD   = metos3d->ybgcinBD;
    Vec         *ybgcoutBD  = metos3d->ybgcoutBD;
    PetscInt    iprof, nlayer;
    PetscScalar *bcarray    = PETSC_NULL;
    PetscScalar *dcarray    = PETSC_NULL;
    PetscScalar *yinarray   = PETSC_NULL;
    PetscScalar *youtarray  = PETSC_NULL;
    PetscFunctionBegin;
    // start log event
    PetscLogEventBegin(metos3d->eventBGCStep, 0, 0, 0, 0);
    // yin, copy to block diagonal, yout, zero entries
    Metos3DUtilVecCopySeparateToDiagonal(metos3d, ntracer, nvecloc, yin, ybgcinBD);
    VecZeroEntries(*ybgcoutBD);
    // boundary condition step
    // domain condition step
    Metos3DBGCStepBoundaryCondition(metos3d, t);
    Metos3DBGCStepDomainCondition(metos3d, t);
    // bc, dc, yin, yout array
    if (nbc > 0) VecGetArray(*metos3d->bgcbcBD, &bcarray);
    if (ndc > 0) VecGetArray(*metos3d->bgcdcBD, &dcarray);
    VecGetArray(*ybgcinBD, &yinarray);
    VecGetArray(*ybgcoutBD, &youtarray);
    // go through profiles
    for (iprof = 0; iprof < nprofloc; iprof++)
    {
        // layer
        nlayer = iendloc[iprof]-istartloc[iprof]+1;        
        //
        // call BGC C API
        //
#ifdef BGC
        BGC(
            (int*)&ntracer,                                     // tracer count
            (int*)&nlayer,                                      // layer count
            (int*)&nparam,                                      // parameter count
            (int*)&nbc,                                         // boundary condition count
            (int*)&ndc,                                         // domain condition count
            
            (int*)&ndc,                                         // diag count
            
            (double*)&dt,                                       // time step
            (double*)&youtarray[ntracer*(istartloc[iprof]-1)],  // source minus sink
            (double*)&t,                                        // point in time
            (double*)&yinarray [ntracer*(istartloc[iprof]-1)],  // tracer
            (double*)u0,                                        // parameter
            (double*)&bcarray  [nbc*iprof],                     // boundary condition
            (double*)&dcarray  [ndc*(istartloc[iprof]-1)],      // domain condition

            (double*)&dcarray  [ndc*(istartloc[iprof]-1)]       // diag variable

            );
#endif        
    }
    // bc, dc, yin, yout array
    if (nbc > 0) VecRestoreArray(*metos3d->bgcbcBD, &bcarray);
    if (ndc > 0) VecRestoreArray(*metos3d->bgcdcBD, &dcarray);
    VecRestoreArray(*ybgcinBD, &yinarray);
    VecRestoreArray(*ybgcoutBD, &youtarray);
    // copy back to separate vectors
    Metos3DUtilVecCopyDiagonalToSeparate(metos3d, ntracer, nvecloc, ybgcoutBD, yout);
    // stop log event
    PetscLogEventEnd(metos3d->eventBGCStep, 0, 0, 0, 0);
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCStep\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCOutput"
PetscErrorCode
Metos3DBGCOutput(Metos3D *metos3d, PetscInt n, Vec *v)
{
    // work vars
    PetscInt    nmax, i_tracer;
    PetscBool   flag = PETSC_FALSE;
    char        outputDirectory     [PETSC_MAX_PATH_LEN];
    char        *outputFileNames    [PETSC_MAX_PATH_LEN];
    char        outputFileNameFormat[PETSC_MAX_PATH_LEN];
    char        format              [PETSC_MAX_PATH_LEN];
    char        filePath            [PETSC_MAX_PATH_LEN];
    PetscFunctionBegin;
    // output tracer vectors to disk
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DTracerOutputDirectory", outputDirectory);
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-Metos3DTracerOutputFileFormat", outputFileNameFormat, PETSC_MAX_PATH_LEN, &flag);
    if (flag == PETSC_TRUE) {
        Metos3DUtilFormatParse(metos3d, outputFileNameFormat);
        for (i_tracer = 0; i_tracer < n; i_tracer++)
        {
            sprintf(format, "%s%s", outputDirectory, outputFileNameFormat);
            sprintf(filePath, format, i_tracer);
            Metos3DUtilVectorView(metos3d, filePath, &v[i_tracer]);
        }
    } else {
        nmax = n;
        PetscOptionsGetStringArray(PETSC_NULL, PETSC_NULL, "-Metos3DTracerOutputFile", outputFileNames, &nmax, &flag);
        if (flag == PETSC_TRUE) {
            for (i_tracer = 0; i_tracer < n; i_tracer++)
            {
                sprintf(filePath, "%s%s", outputDirectory, outputFileNames[i_tracer]);
                Metos3DUtilVectorView(metos3d, filePath, &v[i_tracer]);
                PetscFree(outputFileNames[i_tracer]);
            }
        }
    }
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCOutput\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCOutputPrefix"
PetscErrorCode
Metos3DBGCOutputPrefix(Metos3D *metos3d, char *prefix, PetscInt n, Vec *v)
{
    // work vars
    PetscInt    nmax, i_tracer;
    PetscBool   flag = PETSC_FALSE;
    char        outputDirectory     [PETSC_MAX_PATH_LEN];
    char        *outputFileNames    [PETSC_MAX_PATH_LEN];
    char        outputFileNameFormat[PETSC_MAX_PATH_LEN];
    char        format              [PETSC_MAX_PATH_LEN];
    char        filePath            [PETSC_MAX_PATH_LEN];
    PetscFunctionBegin;
    // output tracer vectors to disk
    Metos3DUtilOptionsGetString(metos3d, "-Metos3DTracerOutputDirectory", outputDirectory);
    PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-Metos3DTracerOutputFileFormat", outputFileNameFormat, PETSC_MAX_PATH_LEN, &flag);
    if (flag == PETSC_TRUE) {
        Metos3DUtilFormatParse(metos3d, outputFileNameFormat);
        for (i_tracer = 0; i_tracer < n; i_tracer++)
        {
            sprintf(format, "%s%s%s", outputDirectory, prefix, outputFileNameFormat);
            sprintf(filePath, format, i_tracer);
            Metos3DUtilVectorView(metos3d, filePath, &v[i_tracer]);
        }
    } else {
        nmax = n;
        PetscOptionsGetStringArray(PETSC_NULL, PETSC_NULL, "-Metos3DTracerOutputFile", outputFileNames, &nmax, &flag);
        if (flag == PETSC_TRUE) {
            for (i_tracer = 0; i_tracer < n; i_tracer++)
            {
                sprintf(filePath, "%s%s%s", outputDirectory, prefix, outputFileNames[i_tracer]);
                Metos3DUtilVectorView(metos3d, filePath, &v[i_tracer]);
                PetscFree(outputFileNames[i_tracer]);
            }
        }
    }
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCOutputPrefix\n");
    PetscFunctionReturn(0);
}

#undef  kDebugLevel
#define kDebugLevel kDebugLevel4

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCStepBoundaryCondition"
PetscErrorCode
Metos3DBGCStepBoundaryCondition(Metos3D *metos3d, PetscScalar t)
{
    // bgc
    PetscInt    nbc         = metos3d->boundaryConditionCount;
    PetscInt    *nbcs       = metos3d->boundaryConditionCountArray;
    Vec         **bcs       = metos3d->boundaryConditionArray;
    Vec         *bgcbc      = metos3d->bgcbc;
    // work vars
    PetscInt    ibc;
    PetscFunctionBegin;
    // check count
    if (nbc > 0) {
        // prepare boundary condition
        for (ibc = 0; ibc < nbc; ibc++) {
            // work vars
            PetscInt    ialpha, ibeta;
            PetscScalar alpha, beta;
            // interpolate
            Metos3DUtilInterpolate(metos3d, t, nbcs[ibc], &ialpha, &alpha, &ibeta, &beta);
            // set
            // bgcbc[ibc] = bcs[ibc][ialpha]
            // bgcbc[ibc] = alpha*bgcbc[ibc]
            // bgcbc[ibc] = bgcbc[ibc] + beta*bcs[ibc][ibeta]
            VecCopy (bcs[ibc][ialpha], bgcbc[ibc]);
            VecScale(bgcbc[ibc], alpha);
            VecAXPY(bgcbc[ibc], beta, bcs[ibc][ibeta]);
        }
        // copy to one vector
        PetscInt    nprofloc = metos3d->profileCountLocal;
        Metos3DUtilVecCopySeparateToDiagonalProfile(metos3d, nbc, nprofloc, bgcbc, metos3d->bgcbcBD);
    }
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCStepBoundaryCondition\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DBGCStepDomainCondition"
PetscErrorCode
Metos3DBGCStepDomainCondition(Metos3D *metos3d, PetscScalar t)
{
    // bgc
    PetscInt    ndc         = metos3d->domainConditionCount;
    PetscInt    *ndcs       = metos3d->domainConditionCountArray;
    Vec         **dcs       = metos3d->domainConditionArray;
    Vec         *bgcdc      = metos3d->bgcdc;
    // work vars
    PetscInt    idc;
    PetscFunctionBegin;
    // check count
    if (ndc > 0) {
        // prepare boundary condition
        for (idc = 0; idc < ndc; idc++) {
            // work vars
            PetscInt    ialpha, ibeta;
            PetscScalar alpha, beta;
            // interpolate
            Metos3DUtilInterpolate(metos3d, t, ndcs[idc], &ialpha, &alpha, &ibeta, &beta);
            // set
            // bgcdc[idc] = dcs[idc][ialpha]
            // bgcdc[idc] = alpha*bgcdc[idc]
            // bgcdc[idc] = bgcdc[idc] + beta*dcs[idc][ibeta]
            VecCopy (dcs[idc][ialpha], bgcdc[idc]);
            VecScale(bgcdc[idc], alpha);
            VecAXPY(bgcdc[idc], beta, dcs[idc][ibeta]);
        }
        // copy to one vector
        PetscInt    nvecloc = metos3d->vectorLengthLocal;
        Metos3DUtilVecCopySeparateToDiagonal(metos3d, ndc, nvecloc, bgcdc, metos3d->bgcdcBD);
    }
    // debug
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DBGCStepDomainCondition\n");
    PetscFunctionReturn(0);
}
