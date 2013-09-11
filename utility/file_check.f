!--------1---------2---------3---------4---------5---------6---------7--

      PROGRAM file_check
      implicit none
      include 'mpif.h'
      integer :: k,n,i,j,ML,MB,MB1,MB2,MBZ,n1,n2,ML0,irank
      integer :: ML1,ML2,ML3,ierr,i1,i2,i3,j1,j2,j3
      integer :: itmp(7),a1,a2,a3,b1,b2,b3
      integer :: nprocs,myrank
      integer,allocatable :: LL(:,:)
      real(8) :: fs,Va,ax,dV
      real(8),allocatable :: rtmp(:),occ(:,:)
      complex(8),allocatable :: utmp(:)

      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world,nprocs,ierr)
      call mpi_comm_rank(mpi_comm_world,myrank,ierr)

      MBZ=4
      ax=205.25978403d0
      Va=ax*ax*0.05d0*ax

      goto 1

      open(1,file='wf.dat',form='unformatted')
      read(1) ML,ML1,ML2,ML3
      read(1) MB,MB1,MB2

      dV=Va/ML

      write(*,*) "Va,dV=",Va,dV
      write(*,*) "ML =",ML
      write(*,*) "MB,MB1,MB2=",MB,MB1,MB2

      allocate( LL(3,ML) ) ; LL=0
      allocate( occ(MB,MBZ) ) ; occ=0.d0

      read(1) LL(:,:)

      read(1) occ(:,:)

      write(*,*) "sum(occ)=",sum(occ)

      allocate( utmp(ML) )

      do k=1,MBZ
      do n=1,MB

         read(1) utmp

         if ( k==4 ) then
         write(*,*) n,k,sum(abs(utmp)**2)*dV
         end if

      end do
      end do

      close(1)

      deallocate( utmp )

 1    continue


      open(1,file='vrho.dat',form='unformatted')
      read(1) ML,ML1,ML2,ML3

      write(*,*) "ML =",ML
      write(*,*) "ML1,ML2,ML3=",ML1,ML2,ML3

      dV=Va/ML

      if ( .not.allocated(LL) ) then
         allocate( LL(3,ML) )
         LL=0
      end if

      read(1) LL(:,:)

      allocate( rtmp(ML) ) ; rtmp=0.d0

      read(1) rtmp(:)

      write(*,*) "sum(rho)=",sum(rtmp)*dV
      write(*,*) "min,max=",minval(rtmp),maxval(rtmp)

      read(1) rtmp(:)
      write(*,*) "min,max (Vloc)=",minval(rtmp),maxval(rtmp)

      read(1) rtmp(:)
      write(*,*) "min,max (Vh)  =",minval(rtmp),maxval(rtmp)

      read(1) rtmp(:)
      write(*,*) "min,max (Vxc) =",minval(rtmp),maxval(rtmp)


      deallocate( rtmp )

      close(1)


      if ( allocated(occ) ) deallocate( occ )
      if ( allocated(LL)  ) deallocate( LL )

      call mpi_finalize(ierr)

      stop
      END PROGRAM file_check
