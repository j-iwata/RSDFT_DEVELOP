#!/bin/sh

DIRPATH=$1

rm USELOG
rm Makefile.dep_o
rm Makefile.dep_mod

for FILE in ${DIRPATH}*.f90
do
  echo ${FILE} | cut -d',' -f1 >> USELOG
  grep ' use ' ${FILE} | cut -d',' -f1 | sed 's/.*use //' >> USELOG
  echo >> USELOG
done
wc -l USELOG
read lineMax
echo $lineMax

prev_kf90='null'

for i in `seq 1 ${lineMax}`
do
  currentLine=`head -n $i USELOG | tail -n 1`
  kf90=`echo $currentLine | cut -d'.' -f2`
  if `test $kf90 = 'f90'` ; then
    moduleName=`echo $currentLine | cut -d'.' -f1`
    echo "${moduleName}.o : ${moduleName}.f90 \\" >> Makefile.dep_o
    echo "${moduleName}.mod : ${moduleName}.f90 ${moduleName}.o \\" >> Makefile.dep_mod
  elif `test $prev_kf90 = 'f90'` ; then
    echo "${moduleName}.o : ${moduleName}.f90 "
  else
    echo not f90
  fi
  prev_kf90=${kf90}
done
exit 0

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
