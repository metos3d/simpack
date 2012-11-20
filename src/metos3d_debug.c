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
Metos3DDebug(Metos3D *metos3d, PetscInt level, const char *format, ...)
{
    PetscFunctionBegin;
    if (metos3d->debugLevel >= level) {
        PetscMPIInt rank;
        MPI_Comm_rank(metos3d->comm, &rank);
        if (!rank) {
            // work vars
            va_list         args;
            PetscLogDouble  endTime, elapsedTime;//, secs, mins, hours;
            char            newformat[PETSC_MAX_PATH_LEN];
            // get end time
            // compute elapsed time
            // and init new timing
            PetscGetTime(&endTime);
            elapsedTime = endTime - metos3d->startTime[level];
//            secs = fmod(elapsedTime, 60.0);
//            mins = fmod(floor(elapsedTime/60.0), 60.0);
//            hours = floor(elapsedTime/3600.0);
//            PetscGetTime(&metos3d->startTime);
            // set new format
//            if (metos3d->debugLevel > level)
//            {
//                // without time
//                if (level==0) sprintf(newformat, "               > %s", format);
//                if (level==1) sprintf(newformat, "               > %s", format);
//                if (level==2) sprintf(newformat, "             : +%s", format);
//                if (level==3) sprintf(newformat, "             : +-%s", format);
//            }
//            else
//            {
//                // with time
//                if (level==0) sprintf(newformat, "%12.3fs: %s", elapsedTime, format);
//                if (level==1) sprintf(newformat, "%12.3fs: %s", elapsedTime, format);
//                if (level==2) sprintf(newformat, "%12.3fs:   %s", elapsedTime, format);
//                if (level==3) sprintf(newformat, "%12.3fs: %s", elapsedTime, format);
//            }
//            // add indent depending on level
//            if (level==1) sprintf(newformat, "%16.2f sec -> %s", elapsedTime, format);
//            if (level==2) sprintf(newformat, "%16.2f sec: >%s", elapsedTime, format);
//            if (level>=3) sprintf(newformat, "%16.2f sec: ->%s", elapsedTime, format);
            
            
//            sprintf(newformat, "%.0fh%02.0fm%.3fs %s", hours, mins, secs, format);
            sprintf(newformat, "%10.1fs %s", elapsedTime, format);
            
            // get variable arg list
            va_start(args, format);
            PetscVFPrintf(PETSC_STDOUT, newformat, args);
            va_end(args);
        }
    }
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DFlag"
PetscErrorCode
Metos3DFlag(PetscBool flag, char *message)
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
