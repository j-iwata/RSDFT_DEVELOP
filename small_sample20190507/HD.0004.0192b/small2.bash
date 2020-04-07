#!/bin/bash
#PJM --rsc-list "node=4"
#PJM --rsc-list "elapse=2:00:00"
#PJM --name "output"
#PJM --rsc-list "rscgrp=small"
#PJM -s

#PJM --mpi "use-rankdir"
#PJM --stgin-dir "./ ./"
#PJM --stgout-dir "rank=* %r:./ ./ recursive=2"

. /work/system/Env_base
export TIMEX="/usr/bin/time -p"
NTHREADS=8
export PARALLEL=${NTHREADS}
export OMP_NUM_THREADS=${NTHREADS}
export FLIB_FASTOMP=TRUE

# timex
TIMEX="/usr/bin/time -p"
#PROF="fapp"
#PROF_ARGS="-C -Hevent_number=10,7,7,11,49,49,31,0 -d Fprofd_${PJM_JOBID}"
#PROF_ARGS="-C -Ihwm -L 1 -Hevent=Instructions -d Fprofd_${PJM_JOBID}"
#MSH="/home/system/bin/msh"
MPIEXEC="mpiexec"
#MPIEXEC_ARGS="-of-proc part.o"${PJM_JOBID}

# Large Page
#LPG="lpgparm"
#LPG_ARGS="-t 256MB -s 256MB -h 256MB -d 256MB -p 256MB"

# Execution
LD="./rsdft.x"

${TIMEX} ${PROF} ${PROF_ARGS} ${MPIEXEC} ${MPIEXEC_ARGS} ${LPG} ${LPG_ARGS} ${LD}
#${TIMEX} ${PROF} ${PROF_ARGS}_${NTHREADS} ${MPIEXEC} ${LPG} ${LPG_ARGS} ${LD} > fDGEMM_Kt${NTHREADS}.txt
