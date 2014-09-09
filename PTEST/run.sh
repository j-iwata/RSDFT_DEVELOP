#!/bin/sh

mpif90 -O1 CallSpeed.f90 > o1.out
#mpif90 -O2 CallSpeed.f90 > o2.out
#mpif90 -O3 CallSpeed.f90 > o3.out

for i in 100 1000 10000 100000
do
  echo $i > fort.1
  mpirun -n 1 a.out > log_$i
  call_time=`grep call log_$i | cut -c 7-28`
  if_time=`grep if log_$i | cut -c 7-28`
  echo $i $call_time >> call_time
  echo $i $if_time >> if_time
done
