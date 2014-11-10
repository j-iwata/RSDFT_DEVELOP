#!/bin/sh

DIRPATH=$1

rm SUBLOG

for FILE in ${DIRPATH}*.f90
do
  echo ${FILE}  >> SUBLOG
  grep 'SUBROUTINE ' ${FILE} >> SUBLOG
  grep 'subroutine ' ${FILE} >> SUBLOG
  grep 'Subroutine ' ${FILE} >> SUBLOG
  echo >> SUBLOG
done
sed -ie '/ END SUBROUTINE/d' SUBLOG
sed -ie '/ end subroutine/d' SUBLOG
sed -ie '/ End Subroutine/d' SUBLOG
sed -ie '/ End subroutine/d' SUBLOG
sed -ie '/!/d' SUBLOG
sed -ie 's///g' SUBLOG
sed -ie 's/SUBROUTINE//g' SUBLOG
sed -ie 's/Subroutine//g' SUBLOG
sed -ie 's/subroutine//g' SUBLOG

exit 0

sed -ie '/omp_lib/d' USELOG
wc -l USELOG
read lineMax
echo $lineMax

for i in `seq 1 ${lineMax}`
do
  currentLine=`head -n $i USELOG | tail -n 1`
  dummy=`echo $currentLine`1
  if `test $dummy = 1` ; then
    continue
  fi
  j=`expr $i + 1`
  nextLine=`head -n $j USELOG | tail -n 1`1
  kf90=`echo $currentLine | cut -d'.' -f2`
  kf90=${kf90}1
  moduleName=`echo $currentLine | cut -d'.' -f1`
  if `test $nextLine = 1` ; then
    if `test $kf90 = f901` ; then
      echo >> Makefile.dep_o
      echo >> Makefile.dep_mod
      echo "${moduleName}.o : ${moduleName}.f90" >> Makefile.dep_o
      echo '	$(FC) $(FFLAGS) -c $<' >> Makefile.dep_o
      echo "${moduleName}.mod : ${moduleName}.f90 ${moduleName}.o" >> Makefile.dep_mod
      echo '	$(FC) $(FFLAGS) -c $<' >> Makefile.dep_mod
    else
      echo "      ${moduleName}.mod" >> Makefile.dep_o
      echo '	$(FC) $(FFLAGS) -c $<' >> Makefile.dep_o
    fi
  else
    if `test $kf90 = f901` ; then
      echo >> Makefile.dep_o
      echo >> Makefile.dep_mod
      echo "${moduleName}.o : ${moduleName}.f90 \\" >> Makefile.dep_o
      echo "${moduleName}.mod : ${moduleName}.f90 ${moduleName}.o" >> Makefile.dep_mod
    else
      echo "      ${moduleName}.mod \\" >> Makefile.dep_o
    fi
  fi
done

sed -i 's///g' Makefile.dep_o
sed -i 's///g' Makefile.dep_mod

exit 0

