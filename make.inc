########################################################################
##### COMPILE & LINK WILL BE:                                      #####
#####   $(FC) -c $(FFLAGS) <SRC>                                   #####
#####   $(FC) $(LFLAGS) <OBJ> $(LAPACK_L)                          #####
##### PLEASE CHECK COMPLETE OPTIONS BY YOURSELF                    #####
########################################################################
FC = mpif90
FFLAGS =
LFLAGS =
LAPACK_L = $(DIR1)/ext_sub_blas.o $(DIR1)/ext_sub_blas_tmp.o \
           $(DIR1)/ext_sub_lapack_c16.o \
           $(DIR1)/ext_sub_lapack_r8.o \
           $(DIR1)/ext_sub_lapack_u.o \
           $(DIR1)/ext_sub_lapack_c16_2.o \
           $(DIR1)/ext_sub_lapack_aux.o
FFLAGS1 = -I../
#FFLAGS1 = -D_DRSDFT_ -I../
#FFLAGS1 = -D_SPLINE_ -I../
#FFLAGS1 = -D_DRSDFT_ -D_SPLINE -I../
########################################################################
########################################################################
##### INTEL                                                        #####
#FFLAGS = -fast -xHost -cpp
#FFLAGS = -fast -xHost -openmp -cpp
#FFLAGS = -fast -xAVX -openmp -cpp
#FFLAGS = -C -openmp -cpp
FFLAGS0 = -C -cpp
LFLAGS = -openmp
#LAPACK_L = -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 \
           -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
           -lguide -lpthread
LAPACK_L = -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 \
           -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
           -lpthread
#LAPACK_L = -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 \
           -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
########################################################################
##### PGF                                                          #####
#FFLAGS = -O3
#FFLAGS = -fastsse
#FFLAGS = -C
#FFLAGS = -fastsse -mp -Mpreprocess -DOMP
#FFLAGS = -C -mp -Mpreprocess -DOMP -DPAPI
#FFLAGS = -C -Mchkfpstk -Mchkptr -Mchkstk -Minform=inform -mp \
         -Mpreprocess -DOMP
#LFLAGS = -mp
#FFLAGS = -C -Mconcur -cpp
#LFLAGS = -Mconcur
#FFLAGS = -fastsse -O3 -tp amd64 -mp -Mpreprocess -DOMP
#FFLAGS = -fastsse -O3 -tp amd64e -mp -Mpreprocess -DOMP
#FFLAGS = -fastsse -O3 -tp barcelona-64 -mp -Mpreprocess -DOMP
#FFLAGS = -fastsse -O3 -tp shanghai-64 -mp -Mpreprocess -DOMP
#LFLAGS = -mp
########################################################################
##### PACS-CS @ Tsukuba                                            #####
#LAPACK_L = -L/opt/intel/cmkl/7.1/lib/em64t -lmkl_scalapack \
           -lmkl_blacs -lmkl_lapack -lmkl_em64t -lguide -lpthread
########################################################################
##### Fujitsu PrimeQuest @ IMS                                     #####
#FC = mpifrt
#FFLAGS = -Am -X9 -Fixed -Kfast -Kparallel -Kpmsg -Qp
#LFLAGS =
#LAPACK_L = -SCALAPACK -KOMP -SSL2BLAMP
########################################################################
##### Hitachi SR16000 @ IMS                                        #####
#FC = mpif90
#FFLAGS = -Oss -nolimit -64 -noparallel
#LFLAGS = -64 -noparallel
#FFLAGS = -Oss -nolimit -64 -omp -parallel -cpp
#LFLAGS = -64 -omp -parallel
#LAPACK_L = -lscalapack -lblacs -lblacsF77init -llapack -lblas
#LAPACK_L = -lf90s -lscalapack -lblacs -lblacsF77init -llapack -lblas
########################################################################
##### Fujitsu FX1 @ Nagoya                                         #####
#FC = mpifrt
#FFLAGS = -KOMP -Cpp -DOMP -DDRSDFT -DRSMOL
#FFLAGS = -KOMP -Cpp -DOMP
#FFLAGS = -KOMP -Cpp -DOMP -Kfast,noeval
#FFLAGS = -Hefosux -KOMP -Cpp -DOMP
#FFLAGS = -KOMP -Cpp -DOMP -Kparallel
#FFLAGS = -Kfast
#FFLAGS = -Kfast,noeval
#FFLAGS = -Ha
#FFLAGS = -Hefosux
#LFLAGS = -KOMP
#LAPACK_L = -SSL2 -SCALAPACK
#LAPACK_L = -SSL2BLAMP -SCALAPACK
########################################################################
##### NEC TX7 @ Tohoku                                             #####
#FFLAGS = -O3 -openmp
#LFLAGS = -O3 -openmp
#LAPACK_L = -L/opt/MathKeisan/lib -lscalapack -lblacsF90init -lblacs \
           -lblacsF90init -llapack -lblas -lmpi
########################################################################
##### NEC SX-9 @ Tohoku                                            #####
#FC = sxmpif90
#FFLAGS = -Ep -Chopt -pi -Wf"-pvctl vwork=stack"
#LFLAGS =
#LAPACK_L = -lscalapack -lblacsF90init -lblacs -lblacsF90init -llapack \
           -lblas -lmpi
########################################################################
##### NEC SX-9 @ Osaka                                             #####
#FC = sxmpif90
#FFLAGS = -Ep -Chopt -pi -Wf"-pvctl vwork=stack" -cfsx9
#LFLAGS = -cfsx9
#LAPACK_L = -lscalapack -lblacsF90init -lblacs -lblacsF90init -llapack \
           -lblas -lmpi
########################################################################
##### NEC SX-9 @ ISSP                                              #####
#FC = sxmpif90
#FFLAGS = -Ep -Chopt -pi -Wf"-pvctl vwork=stack"
#LFLAGS =
#LAPACK_L = -lscalapack -lblacsF90init -lblacs -lblacsF90init -llapack \
           -lblas -lmpi
#LAPACK_L = -lscalapack_64 -lblacsF90init_64 -lblacs_64 \
           -lblacsF90init_64 -llapack_64 -lblas_64 -lmpiw
########################################################################
##### SGI Altix ICE 8400EX @ ISSP                                  #####
#FC = ifort
#FFLAGS = -openmp -cpp -DOMP
#FFLAGS = -fast -xSSE4.2 -openmp -cpp -DOMP
#LFLAGS = -openmp
#LAPACK_L = -lmkl_scalapack_lp64 -lmkl_blacs_sgimpt_lp64 -lmkl_lapack \
           -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmpi
#LAPACK_L = -lmkl_scalapack_lp64 -lmkl_blacs_sgimpt_lp64 -lmkl_lapack \
           -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lmpi
########################################################################
##### T2K-Tsukuba (INTEL)                                          #####
#FC = mpif90
#FFLAGS = -CB
#FFLAGS = -CB -static
#FFLAGS = -O3 -cpp
#FFLAGS = -O4
#FFLAGS = -tpp7
#FFLAGS = -axW
#FFLAGS = -C -openmp -cpp -DOMP -DPAPI
#FFLAGS = -C -openmp -cpp -DOMP
#LFLAGS = -openmp
#FFLAGS = -C -parallel -cpp
#LFLAGS = -parallel
#LAPACK_L = -L/opt/intel/mkl/10.0.1.014/lib/em64t -lmkl_scalapack \
           -lmkl_blacs_lp64 -lmkl_em64t -lguide
#LDIR = /global/usr/mpi/mvapich2-1.2-10-01/intel/lib
#LDIR = /opt/mpi/intel/mvapich2-1.5.1-medium/lib
#LAPACK_L = $(LDIR)/libscalapack.a $(LDIR)/blacsF77init_MPI-LINUX-0.a \
           $(LDIR)/blacs_MPI-LINUX-0.a \
           $(LDIR)/blacsCinit_MPI-LINUX-0.a -lacml
#LAPACK_L = $(LDIR)/libscalapack.a $(LDIR)/blacsF77init_MPI-LINUX-0.a \
           $(LDIR)/blacs_MPI-LINUX-0.a \
           $(LDIR)/blacsCinit_MPI-LINUX-0.a \
           -L/opt/intel/mkl/10.0.1.014/lib/em64t -lmkl_em64t -lguide
########################################################################
##### T2K-Tsukuba (PGF)                                            #####
#FC = mpif90
#FFLAGS = -O3
#FFLAGS = -fastsse
#FFLAGS = -C
#FFLAGS = -O3 -Mpreprocess -DDRSDFT -DRSMOL
#FFLAGS = -fastsse -mp -Mpreprocess -DOMP
#FFLAGS = -C -mp -Mpreprocess -DOMP -DPAPI
#FFLAGS = -C -Mchkfpstk -Mchkptr -Mchkstk -Minform=inform -mp \
         -Mpreprocess -DOMP
#LFLAGS = -mp
#FFLAGS = -C -Mconcur -cpp
#LFLAGS = -Mconcur
#LDIR = /opt/pgi/linux86-64/7.2/mpi/mvapich/lib
#LDIR = /global/usr/mpi/mvapich2-1.2-10-01/pgi/lib
#LDIR = /opt/mpi/pgi/mvapich2-1.4.1/lib
#LDIR = /opt/mpi/pgi/mvapich2-1.5.1-medium/lib
#LAPACK_L = $(LDIR)/libscalapack.a $(LDIR)/blacsF77init_MPI-LINUX-0.a \
           $(LDIR)/blacs_MPI-LINUX-0.a \
           $(LDIR)/blacsCinit_MPI-LINUX-0.a -lacml
#LAPACK_L = $(LDIR)/libscalapack.a $(LDIR)/blacsF77init_MPI-LINUX-0.a \
           $(LDIR)/blacs_MPI-LINUX-0.a \
           $(LDIR)/blacsCinit_MPI-LINUX-0.a \
           /opt/pgi/linux86-64/7.2-2/lib/libacml_mp.a
# /opt/papi/papi-3.6.2/pgi/lib/libpapi.a
########################################################################
##### T2K-Tokyo                                                    #####
#FC = mpif90
#FFLAGS = -Oss -noparallel -autoinline -hf90
#LFLAGS =
#FFLAGS = -O3 -ipo -xP
#LFLAGS = -ipo
#FFLAGS = -fastsse -O3 -tp barcelona-64
#LFLAGS = -fastsse -O3
#LAPACK_L = -llapack_sc -lblas_sc
#LAPACK_L = -L/opt/intel/mkl/10.1.0.015/lib/em64t -lmkl_scalapack_lp64 \
           -lmkl_blacs_lp64 -lmkl_lapack -lmkl_intel_lp64 \
           -lmkl_sequential -lmkl_core -lpthread
#LAPACK_L = -lacml
########################################################################
#####  K computer                                                  #####
#FC = mpifrtpx
#FFLAGS = -Hs -g -Kvisimpact,array_private,auto,ilfunc,ocl,NOFLTLD \
         -Kopenmp -Am -X9 -V -Ec -Qa,d,i,p,t,x -Koptmsg=2 -V \
         -Cpp -Ccpp -DOMP -D_FPCOLL_
#FFLAGS = -Kfast,array_private,auto,ilfunc,ocl,NOFLTLD,region_extension\
         -Kopenmp -Am -X9 -V -Ec -Qa,d,i,p,t,x -Koptmsg=2 -V \
         -Cpp -Ccpp -DOMP -D_FPCOLL_
#FFLAGS = -Kvisimpact,array_private,auto,ilfunc,ocl,NOFLTLD \
         -Kopenmp -Am -X9 -V -Ec -Qa,d,i,p,t,x -Koptmsg=2 -V \
         -Cpp -Ccpp -DOMP -D_FPCOLL_
#FFLAGS = -Kvisimpact,preex,array_private,auto,ocl,NOFLTLD \
        -Ksimd=2 -Kopenmp -Ktemparraystack -Koptmsg=2 \
        -X9 -V -Ec -Qa,d,i,p,t,x -Cpp -Ccpp
#LFLAGS = ${FFLAGS} -SSL2BLAMP -SCALAPACK -Ntl_trt
#LAPACK_L =
########################################################################
########################################################################
#
FFLAGS = $(FFLAGS0) $(FFLAGS1)
#
########################################################################
.SUFFIXES:
.SUFFIXES: .o .f .f90
.f.o :
	$(FC) $(FFLAGS) -c $<
.f90.o :
	$(FC) $(FFLAGS) -c $<
