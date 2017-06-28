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

#ifndef METOS3D_BGC_H
#define METOS3D_BGC_H 1

#include "metos3d_load.h"

//
//  BGC API
//

// define BGCINIT, BGCBEGIN, BGC, BGCEND, BGCFINAL Makefile

#ifdef BGCINIT
extern void BGCINIT (int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, int*, double*);
#endif

#ifdef BGCBEGIN
extern void BGCBEGIN(int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, int*, double*);
#endif

#ifdef BGC
extern void BGC     (int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, int*, double*);
#endif

#ifdef BGCEND
extern void BGCEND  (int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, int*, double*);
#endif

#ifdef BGCFINAL
extern void BGCFINAL(int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, int*, double*);
#endif

// init
extern PetscErrorCode Metos3DBGCInit(Metos3D*);
extern PetscErrorCode Metos3DBGCTracerInit(Metos3D*);
extern PetscErrorCode Metos3DBGCDiagInit(Metos3D*);
extern PetscErrorCode Metos3DBGCParameterInit(Metos3D*);
extern PetscErrorCode Metos3DBGCBoundaryConditionInit(Metos3D*);
extern PetscErrorCode Metos3DBGCDomainConditionInit(Metos3D*);
// final
extern PetscErrorCode Metos3DBGCFinal(Metos3D*);
extern PetscErrorCode Metos3DBGCDomainConditionFinal(Metos3D*);
extern PetscErrorCode Metos3DBGCBoundaryConditionFinal(Metos3D*);
extern PetscErrorCode Metos3DBGCParameterFinal(Metos3D*);
extern PetscErrorCode Metos3DBGCDiagFinal(Metos3D*);
extern PetscErrorCode Metos3DBGCTracerFinal(Metos3D*);
// step
extern PetscErrorCode Metos3DBGCStepBegin(Metos3D*, PetscReal, PetscReal, Vec*, Vec*, PetscInt, PetscReal*);
extern PetscErrorCode Metos3DBGCStepEnd(Metos3D*, PetscReal, PetscReal, Vec*, Vec*, PetscInt, PetscReal*);
extern PetscErrorCode Metos3DBGCStep(Metos3D*, PetscReal, PetscReal, Vec*, Vec*, PetscInt, PetscReal*);
extern PetscErrorCode Metos3DBGCStepBoundaryCondition(Metos3D*, PetscReal);
extern PetscErrorCode Metos3DBGCStepDomainCondition(Metos3D*, PetscReal);
// output
extern PetscErrorCode Metos3DBGCOutput(Metos3D*, PetscInt, Vec*);
extern PetscErrorCode Metos3DBGCOutputPrefix(Metos3D*, char*, PetscInt, Vec*);
extern PetscErrorCode Metos3DBGCTracerMonitor(Metos3D*, Vec*);
extern PetscErrorCode Metos3DBGCDiagMonitor(Metos3D*);

#endif /* !METOS3D_BGC_H */


