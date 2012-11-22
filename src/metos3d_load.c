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
 *  metos3d_load.c
 *
 */

#include "metos3d_load.h"

#undef  kDebugLevel
#define kDebugLevel kDebugLevel2

#undef  __FUNCT__
#define __FUNCT__ "Metos3DLoadInit"
PetscErrorCode
Metos3DLoadInit(Metos3D *metos3d)
{
    MPI_Comm    comm    = metos3d->comm;
    PetscInt    debug   = metos3d->debugLevel;
    // geometry
    PetscInt    nprof   = metos3d->profileCount;
    PetscInt    nvec    = metos3d->vectorLength;
    PetscInt    *istart = metos3d->profileStart;
    PetscInt    *iend   = metos3d->profileEnd;
    // work vars
    PetscInt    nproc, myproc, nvecopt, iproc, iprof;
    PetscInt    *sumvecloc, *sumvecprev, *sumprofloc, *sumprofprev;
    PetscInt    nvecloc, nvecprev, nprofloc, nprofprev;
    PetscInt    *istartloc, *iendloc;
    PetscFunctionBegin;
//    // debug start
//    PetscGetTime(&metos3d->startTime[kDebugLevel]);
    // determine process count and my process number
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &myproc);
    // store
    metos3d->processCount = nproc;
    metos3d->processMine = myproc;
//    Metos3DDebug(metos3d, kDebugLevel, F2SD, "Metos3DLoadInit", "processCount:", nproc);
    // compute simple load distribution
    PetscMalloc(nproc*sizeof(PetscInt), &sumvecloc);
    PetscMalloc(nproc*sizeof(PetscInt), &sumvecprev);
    PetscMalloc(nproc*sizeof(PetscInt), &sumprofloc);
    PetscMalloc(nproc*sizeof(PetscInt), &sumprofprev);
    // zero memory
    PetscMemzero(sumvecloc, nproc*sizeof(PetscInt));
    PetscMemzero(sumvecprev, nproc*sizeof(PetscInt));
    PetscMemzero(sumprofloc, nproc*sizeof(PetscInt));
    PetscMemzero(sumprofprev, nproc*sizeof(PetscInt));
    // compute local first
    nvecopt = nvec/nproc;
    iproc   = 0;
    for (iprof = 0; iprof < nprof; iprof++)
    {
        // sum profiles and vector elements per process
        sumvecloc [iproc] = sumvecloc [iproc] + iend[iprof] - istart[iprof] + 1;
        sumprofloc[iproc] = sumprofloc[iproc] + 1;
        if (sumvecloc[iproc] >= nvecopt) iproc++;
    }
    // compute previous then
    sumvecprev [0] = 0;
    sumprofprev[0] = 0;
    for (iproc = 1; iproc < nproc; iproc++)
    {
        sumvecprev [iproc] = sumvecprev [iproc-1] + sumvecloc [iproc-1];
        sumprofprev[iproc] = sumprofprev[iproc-1] + sumprofloc[iproc-1];
    }
    // compute local indices (1 indexed, same as global)
    nvecloc   = sumvecloc[myproc];
    nprofloc  = sumprofloc[myproc];
    nvecprev  = sumvecprev[myproc];
    nprofprev = sumprofprev[myproc];
    // free sums
    PetscFree(sumvecloc);
    PetscFree(sumvecprev);
    PetscFree(sumprofloc);
    PetscFree(sumprofprev);
    // alloc memory
    PetscMalloc(nprofloc*sizeof(PetscInt), &istartloc);
    PetscMalloc(nprofloc*sizeof(PetscInt), &iendloc);
    for (iprof = 0; iprof < nprofloc; iprof++)
    {
        istartloc[iprof] = istart[iprof+nprofprev] - nvecprev;
        iendloc  [iprof] = iend  [iprof+nprofprev] - nvecprev;
    }
    // store
    metos3d->vectorLengthLocal      = nvecloc;
    metos3d->vectorLengthPrevious   = nvecprev;
    metos3d->profileCountLocal      = nprofloc;
    metos3d->profileCountPrevious   = nprofprev;
    metos3d->profileStartLocal      = istartloc;
    metos3d->profileEndLocal        = iendloc;
    // debug
    Metos3DSynchronizedDebug11(debug, kDebugLevel3, comm, SYNCFS5SD, "Metos3DLoadInit", "myproc:", myproc,
        "nvecloc:", nvecloc, "nvecprev:", nvecprev, "nprofloc:", nprofloc, "nprofprev:", nprofprev);
    PetscSynchronizedFlush(comm);
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DLoadInit\n");
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "Metos3DLoadFinal"
PetscErrorCode
Metos3DLoadFinal(Metos3D *metos3d)
{
    PetscFunctionBegin;
//    // debug start
//    PetscGetTime(&metos3d->startTime[kDebugLevel]);
    // local indices
    PetscFree(metos3d->profileStartLocal);
    PetscFree(metos3d->profileEndLocal);
    // debug stop
    Metos3DDebug(metos3d, kDebugLevel, "Metos3DLoadFinal\n");
    PetscFunctionReturn(0);
}
