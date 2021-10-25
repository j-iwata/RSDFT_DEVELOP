program gen_restart_dat
  implicit none
  integer,parameter :: u=10
  integer :: nb,nk,ns,ib,ik,is,np(7)
  real(8) :: ne,x
  real(8),allocatable :: occ(:,:,:)
  write(*,*) "Input Nband, Nk, and Nspin = ? ? ?"
  read(*,*) nb, nk, ns
  write(*,*) "Input # of electrons = ?"
  read(*,*) ne

  allocate( occ(nb,nk,ns) ); occ=0.0d0

  x = 0.0d0
  do is = 1, ns
    do ik = 1, nk
      do ib = 1, nb
        x = x + 2.0d0/ns
        if ( x < ne+1.d-10 ) occ(ib,ik,is)=2.0d0/ns
      end do
    end do
  end do

  rewind u
  write(u,*) nb, nk, ns
  do is = 1, ns
    do ik = 1, nk
      do ib = 1, nb
        write(u,'(3i6,2x,g22.15)') ib,ik,is,occ(ib,ik,is)
      end do
    end do
  end do

  write(*,*) "Input MPI process partition = ? ? ? ? ? ? ?"
  np=1
  read(*,*) np(1:7)
  write(u,'(7i6)') np

end program gen_restart_dat
