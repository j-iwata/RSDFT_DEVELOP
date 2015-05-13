mkdir HOGE
for i in *.f90; 
do 
     cat $i | sed s/MPI_COMPLEX16/MPI_DOUBLE_COMPLEX/g | sed s/mpi_complex16/MPI_DOUBLE_COMPLEX/g > HOGE/$i; 
done
mv HOGE/*f90 .
rm -rf HOGE
