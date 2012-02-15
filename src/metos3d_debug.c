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
 *  metos3d_debug.c
 *
 */

#include "metos3d_debug.h"

#undef  __FUNCT__
#define __FUNCT__ "Metos3DDebug"
PetscErrorCode
Metos3DDebug(PetscInt debugLevel, PetscInt level, MPI_Comm comm, const char *format, ...)
{
    char newformat[PETSC_MAX_PATH_LEN];
    PetscFunctionBegin;
    if (debugLevel >= level) {
        PetscMPIInt rank;
        MPI_Comm_rank(comm, &rank);
        if (!rank) {
            va_list args;
            va_start(args, format);
            // add indent depending on level
            if (level==1) sprintf(newformat, "---> %s", format);
            if (level==2) sprintf(newformat, " --> %s", format);
            if (level>=3) sprintf(newformat, "  -> %s", format);
            PetscVFPrintf(PETSC_STDOUT, newformat, args);
            va_end(args);
        }
    }
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DFlag"
PetscErrorCode
Metos3DFlag(PetscTruth flag, char *message)
{
    PetscFunctionBegin;
    if (flag == PETSC_FALSE) {
        PetscPrintf(PETSC_COMM_WORLD, "\n");
        PetscPrintf(PETSC_COMM_WORLD, "### ERROR:\n");
        PetscPrintf(PETSC_COMM_WORLD, "### ERROR: %s\n", message);
        PetscPrintf(PETSC_COMM_WORLD, "### ERROR:\n");
        PetscPrintf(PETSC_COMM_WORLD, "\n");
        PetscEnd();
    }
    PetscFunctionReturn(0);
}
