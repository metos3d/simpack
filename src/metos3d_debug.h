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
 *  metos3d_debug.h
 *
 */

#ifndef METOS3D_DEBUG_H
#define METOS3D_DEBUG_H 1

#include "metos3d_type.h"
#include <stdarg.h>

#define kDebugLevel0 0
#define kDebugLevel1 1
#define kDebugLevel2 2
#define kDebugLevel3 3
#define kDebugLevel4 4

// formats
#define F2S     "%-40s %18s\n"                                                        // 2 x string
#define F2SE    "%-40s %18s %-8e\n"                                                   // 2 x string, 1 x double
#define F3S     "%-40s %18s %-42s\n"                                                // 3 x string
#define F2SD    "%-40s %18s %-8d\n"                                                 // 2 x string, 1 x int
#define FSSDSE  "%-40s %18s %-8d %14s %-8e\n"                                       // 2 x string, 1 x int, 1 x string, 1 x double
#define FS2SESD "%-40s %18s %-8e %14s %-8d %14s %-8e %14s %-8d\n"
#define F4SD    "%-40s %18s %-42s %14s %-8d\n"                                      // 4 x string, 1 x int
#define F4SE    "%-40s %18s %-42s %14s %-8e\n"                                      // 4 x string, 1 x double
#define F5S     "%-40s %18s %-42s %14s %-42s\n"                                     // 5 x string
#define FS5SD   "%-40s %18s %-8d, %14s %-8d, %14s %-8d, %14s %-8d, %14s %-8d\n"     // 1 x string, 5 x (1 x string, 1 x int)
#define FDSE    "%0004d %s %.12e\n"

#define SYNCFS5SD "              ... %-40s     %14s %-8d %14s %-8d %14s %-8d %14s %-8d %14s %-8d\n"  // 1 x string, 5 x (1 x string, 1 x int)

extern PetscErrorCode Metos3DDebug(Metos3D*, PetscInt, const char*, ...);
extern PetscErrorCode Metos3DFlag(PetscBool, char*);

#define Metos3DSynchronizedDebug(du,dl,c,m,a) {if(du>=dl) PetscSynchronizedPrintf(c,m,a);}
#define Metos3DSynchronizedDebug2(du,dl,c,m,a1,a2) {if(du>=dl) PetscSynchronizedPrintf(c,m,a1,a2);}
#define Metos3DSynchronizedDebug3(du,dl,c,m,a1,a2,a3) {if(du>=dl) PetscSynchronizedPrintf(c,m,a1,a2,a3);}
#define Metos3DSynchronizedDebug11(du,dl,c,m,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11) {if(du>=dl) PetscSynchronizedPrintf(c,m,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11);}

#endif /* !METOS3D_DEBUG_H */
