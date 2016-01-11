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

#include "metos3d.h"

int
main(int argc, char **args)
{
    PetscInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);
    if (argc >= 2) {
        // data type
        Metos3D metos3d;
        // init solver context
        // solve
        // final solver context
        Metos3DInitWithFilePath(&metos3d, args[1]);
        Metos3DSolver(&metos3d);
        Metos3DFinal(&metos3d);
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "### ERROR:\n");
        PetscPrintf(PETSC_COMM_WORLD, "### ERROR: %s\n", "Please provide an option file!");
        PetscPrintf(PETSC_COMM_WORLD, "### ERROR:\n");
    }
    PetscFinalize();
    return 0;
}
