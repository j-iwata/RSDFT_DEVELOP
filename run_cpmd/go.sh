#!/bin/sh -x
#PJM -L elapse=00:30:00
#PJM -L node=1:noncont
#PJM --mpi proc=4
#PJM -j
#PJM -S

export LANG=C

export NPROCS=4
export PARALLEL=12
export OMP_NUM_THREADS=${PARALLEL}
export OMP_STACKSIZE=12m
export FLIB_FASTOMP=TRUE

LANG_HOME=/opt/local/aa64own
export PATH=${LANG_HOME}/bin:${PATH}
export LD_LIBRARY_PATH=${LANG_HOME}/lib64:${LD_LIBRARY_PATH}
export OPAL_PREFIX=${LANG_HOME}

#
# execution
#

/bin/cp -p fort.1.p${NPROCS} fort.1

/usr/bin/time -p mpiexec -n ${NPROCS} -of-proc stdout ../srcT2/rsdft.x

