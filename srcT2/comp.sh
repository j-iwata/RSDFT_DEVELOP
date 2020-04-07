#!/bin/bash -x

export PATH=${LANG_HOME}/bin:${PATH}
export LD_LIBRARY_PATH=${LANG_HOME}/ib64:${LD_LIBRARY_PATH}
export OPAL_PREFIX=${LANG_HOME}

#make clean
make
