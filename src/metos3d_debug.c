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
            PetscLogDouble  endTime, elapsedTime;
            char            newformat[PETSC_MAX_PATH_LEN];
            // get end time
            // compute elapsed time
            // and init new timing
            PetscTime(&endTime);
            if (level == 0)
                elapsedTime = endTime - metos3d->startTime[kDebugLevel0];
            else
                elapsedTime = endTime - metos3d->startTime[kDebugLevel1];
            
            if (metos3d->debugLevel == 0) sprintf(newformat, "%12.3fs %s", elapsedTime, format);
            if (metos3d->debugLevel == 1)
            {
                if (level==0) sprintf(newformat, "%12.3fs %s", elapsedTime, format);
                if (level==1) sprintf(newformat, "%12.3fs %s", elapsedTime, format);
            }
            if (metos3d->debugLevel == 2)
            {
                if (level==0) sprintf(newformat, "%12.3fs    %s", elapsedTime, format);
                if (level==1) sprintf(newformat, "%12.3fs .  %s", elapsedTime, format);
                if (level==2) sprintf(newformat, "%12.3fs .. %s", elapsedTime, format);
            }
            if (metos3d->debugLevel == 3)
            {
                if (level==0) sprintf(newformat, "%12.3fs     %s", elapsedTime, format);
                if (level==1) sprintf(newformat, "%12.3fs .   %s", elapsedTime, format);
                if (level==2) sprintf(newformat, "%12.3fs ..  %s", elapsedTime, format);
                if (level==3) sprintf(newformat, "%12.3fs ... %s", elapsedTime, format);
            }
            if (metos3d->debugLevel == 4)
            {
                if (level==0) sprintf(newformat, "%12.3fs      %s", elapsedTime, format);
                if (level==1) sprintf(newformat, "%12.3fs .    %s", elapsedTime, format);
                if (level==2) sprintf(newformat, "%12.3fs ..   %s", elapsedTime, format);
                if (level==3) sprintf(newformat, "%12.3fs ...  %s", elapsedTime, format);
                if (level==4) sprintf(newformat, "%12.3fs .... %s", elapsedTime, format);
            }
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
