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
 *  metos3d_util.h
 *
 */

#ifndef METOS3D_UTIL_H
#define METOS3D_UTIL_H 1

#include "metos3d_debug.h"
#include "math.h"

// load
extern PetscErrorCode Metos3DUtilVectorLoad(Metos3D*, char*, Vec*);
extern PetscErrorCode Metos3DUtilMatrixLoad(Metos3D*, char*, Mat*);
extern PetscErrorCode Metos3DUtilMatrixLoadAndCreate(Metos3D*, char*, Mat*);
// save
extern PetscErrorCode Metos3DUtilVectorView(Metos3D*, char*, Vec*);
// format
extern PetscErrorCode Metos3DUtilFormatParse(Metos3D*, char*);
// copy
// vector
extern PetscErrorCode Metos3DUtilVecCopySeparateToDiagonal(Metos3D*, PetscInt, PetscInt, Vec*, Vec*);
extern PetscErrorCode Metos3DUtilVecCopyDiagonalToSeparate(Metos3D*, PetscInt, PetscInt, Vec*, Vec*);
// profile
extern PetscErrorCode Metos3DUtilVecCopySeparateToDiagonalProfile(Metos3D*, PetscInt, PetscInt, Vec*, Vec*);
// create
extern PetscErrorCode Metos3DUtilVecCreateAndSetValue(Metos3D*, PetscInt, PetscInt, PetscInt, Vec**, PetscScalar);
extern PetscErrorCode Metos3DUtilVecCreateAndSetValues(Metos3D*, PetscInt, PetscInt, PetscInt, Vec**, PetscScalar*);
// options
extern PetscErrorCode Metos3DUtilOptionsGetInt(Metos3D*, const char*, PetscInt*);
extern PetscErrorCode Metos3DUtilOptionsGetScalar(Metos3D*, const char*, PetscScalar*);
extern PetscErrorCode Metos3DUtilOptionsGetRealArray(Metos3D*, const char*, PetscInt*, PetscReal*);
extern PetscErrorCode Metos3DUtilOptionsGetString(Metos3D*, const char*, char*);
// interpolate
extern PetscErrorCode Metos3DUtilInterpolate(Metos3D*, PetscReal, PetscInt, PetscInt*, PetscScalar*, PetscInt*, PetscScalar*);

#endif /* !METOS3D_UTIL_H */
