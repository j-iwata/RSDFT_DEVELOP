#!/bin/sh

grep 'use ' *.f90 | cut -d'!' -f2 | cut -d',' -f1 > use.log1
grep 'f90' use.log1 > USELOG
flineMax=`wc -l USELOG | tr -d [:alpha:]`
echo $flineMax
read lineMax
echo $lineMax
prevUsingModule='null'
prev2UsingModule='null'
prevUsedModule='null'

echo '#producing explicit make dependency' > Makefile.dep_o
echo '#producing explicit make dependency' > Makefile.dep_mod

for i in `seq 1 ${lineMax}`
do
  usingModule=`head -n $i use.log2 | tail -n 1 | cut -d'.' -f1`
  usedModule=`head -n $i use.log2 | tail -n 1 | sed 's/.*use //'`
  if `test ${usingModule} = ${prevUsingModule}` ; then
    echo "      ${usedModule}.mod \\" >> Makefile.dep_o
    echo "      ${usedModule}.o \\" >> Makefile.dep_mod
  elif `test ${prev2UsingModule} = ${prevUsingModule}` ; then
    sed -ie '$d' Makefile.dep_o
    sed -ie '$d' Makefile.dep_mod
    echo "      ${prevUsedModule}.mod" >> Makefile.dep_o
    echo "      ${prevUsedModule}.o" >> Makefile.dep_mod
    echo >> Makefile.dep_o
    echo >> Makefile.dep_mod
    echo "${usingModule}.o : ${usedModule}.mod \\" >> Makefile.dep_o
    echo "${usingModule}.mod : ${usedModule}.o \\" >> Makefile.dep_mod
  else
    sed -ie '$d' Makefile.dep_o
    sed -ie '$d' Makefile.dep_mod
    echo "${prevUsingModule}.o : ${prevUsedModule}.mod" >> Makefile.dep_o
    echo "${prevUsingModule}.mod : ${prevUsedModule}.o" >> Makefile.dep_mod
    echo >> Makefile.dep_o
    echo >> Makefile.dep_mod
    echo "${usingModule}.o : ${usedModule}.mod \\" >> Makefile.dep_o
    echo "${usingModule}.mod : ${usedModule}.o \\" >> Makefile.dep_mod
  fi
  prev2UsingModule=$prevUsingModule
  prevUsingModule=$usingModule
  prevUsedModule=$usedModule
done

gsed -i 's///g' Makefile.dep_o
gsed -i 's///g' Makefile.dep_mod
exit 0
sed -ie 'd' Makefile.dep_o
sed -ie 'd' Makefile.dep_mod
