#!/bin/sh

head -n 1 *.f90 > module.log
grep '==>' module.log | cut -d' ' -f2 | cut -d'.' -f1 > module_head.log
grep 'MODULE' module.log | sed 's/.*MODULE //' > module_name.log

gsed -i 's///g' module.log
gsed -i 's///g' module_head.log
gsed -i 's///g' module_name.log

diff module_head.log module_name.log
exit 0



grep 'f90' use.log1 > USELOG
flineMax=`wc -l USELOG | tr -d [:alpha:]`
echo $flineMax
read lineMax
echo $lineMax
prevUsingModule='null'
prev2UsingModule='null'
prevUsedModule='null'

echo '#producing explicit make dependency' > Makefile.dep

for i in `seq 1 ${lineMax}`
do
  usingModule=`head -n $i use.log2 | tail -n 1 | cut -d'.' -f1`
  usedModule=`head -n $i use.log2 | tail -n 1 | sed 's/.*use //'`
  if `test ${usingModule} = ${prevUsingModule}` ; then
    echo "      ${usedModule}.o ${usedModule}.mod \\" >> Makefile.dep
  elif `test ${prev2UsingModule} = ${prevUsingModule}` ; then
    sed -ie '$d' Makefile.dep
    echo "      ${prevUsedModule}.o ${prevUsedModule}.mod" >> Makefile.dep
    echo >> Makefile.dep
    echo "${usingModule}.o ${usingModule}.mod : ${usedModule}.o ${usedModule}.mod \\" >> Makefile.dep
  else
    sed -ie '$d' Makefile.dep
    echo "${prevUsingModule} ${prevUsingModule}.mod : ${prevUsedModule}.o ${prevUsedModule}.mod" >> Makefile.dep
    echo >> Makefile.dep
    echo "${usingModule}.o ${usingModule}.mod : ${usedModule}.o ${prevUsedModule}.mod \\" >> Makefile.dep
  fi
  prev2UsingModule=$prevUsingModule
  prevUsingModule=$usingModule
  prevUsedModule=$usedModule
done

gsed -i 's///g' Makefile.dep
