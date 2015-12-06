#!/bin/sh

alias sed='gsed'

DIRPATH=$1

rm NameChecker

for FILE in ${DIRPATH}*.f90
do
  nameFile=`echo ${FILE} | cut -d'.' -f1`
  nameModule=`grep -i 'MODULE' ${FILE} | head -n 1 | sed 's/MODULE //i' | cut -d' ' -f1 | sed 's///g'`
  if [ -z "$nameModule" ] ; then
    echo "'$FILE' is not a module." >> NameChecker
    echo "" >> NameChecker
    continue
  fi
  if `test ${nameFile} = ${nameModule}` ; then
    continue
  else
    echo "'${FILE}' has a different module name." >> NameChecker
    echo "'${nameModule}' is its module name" >> NameChecker
    echo "" >> NameChecker
  fi
done

#sed -ie 's///g' SUBLOG

exit 0

