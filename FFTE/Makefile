include ../make.inc

#FFLAGS_FFTE = -fast -xHOST -openmp

all : pzfft3dv.o fft235.o kernel.o zfft3d.o pzfft3d.o

#kernel.o:
#	$(FC) $(FFLAGS_FFTE) -c kernel.f
#pzfft3dv.o:
#	$(FC) $(FFLAGS_FFTE) -c pzfft3dv.f
#fft235.o:
#	$(FC) $(FFLAGS_FFTE) -c fft235.f
#zfft3d.o:
#	$(FC) $(FFLAGS_FFTE) -c zfft3d.f
#pzfft3d.o:
#	$(FC) $(FFLAGS_FFTE) -c pzfft3d.f

clean :
	rm -r -f *.o mpif.h *.lst
