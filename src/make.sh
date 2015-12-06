#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS -j oe
#PBS -o OUT
#PBS -N make

CURRENT=$PBS_O_WORKDIR
#CURRENT=$PWD
cd $CURRENT

time make -j 16 > make.log
