MODULE ps_gth_module

  use atom_module, only: Nelement
  use pseudopot_module

  implicit none

  PRIVATE
  PUBLIC :: read_ps_gth,hnl,no,Rcloc,Rps0

  real(8),allocatable :: Rps0(:,:)
  integer,allocatable :: no(:,:)
  real(8),allocatable :: Rcloc(:),hnl(:,:,:)

CONTAINS

  SUBROUTINE read_ps_gth(rank)
    implicit none
    integer,intent(IN) :: rank
    integer,parameter :: lrefmax=2
    integer,parameter :: unit_pp=34
    integer :: i,j,ielm,ierr
    character(2) :: name
    real(8) :: znuc,rcnl,work(2)
    include 'mpif.h'
    integer :: MMr,iorb,L,n
    real(8) :: pi,rmax,dr,ep,gamma,const1,const2,v,r

    allocate( Rcloc(Nelement) ) ; Rcloc=0.d0
    allocate( parloc(1:4,Nelement) ) ; parloc=0.d0
    allocate( Zps(Nelement)  ) ; Zps=0.d0
    allocate( norb(Nelement) ) ; norb=0
    allocate( lo(4,Nelement) ) ; lo=0
    allocate( no(4,Nelement) ) ; no=0
    allocate( Rps(4,Nelement) ) ; Rps=0.d0
    allocate( Rps0(4,Nelement) ) ; Rps0=0.d0
    allocate( hnl(2,0:1,Nelement) ) ; hnl=0.d0
    allocate( inorm(4,Nelement) ) ; inorm=0

    if ( rank == 0 ) then
       write(*,'(a60," read_ps_gth")') repeat("-",60)
       do ielm=1,Nelement
          open(unit_pp,file=file_ps(ielm),status='old')
          read(unit_pp,*) name,znuc,Zps(ielm)
          read(unit_pp,*) Rcloc(ielm),parloc(1:4,ielm)
          write(*,*) name,znuc,Zps(ielm)
          write(*,'(1x,f14.8,4f14.8)') Rcloc(ielm),parloc(1:4,ielm)
          do i=1,lrefmax
             work(1:2)=0.d0
             read(unit_pp,*,END=10) rcnl,work(1:2)
             do j=1,2
                if ( work(j) /= 0.d0 ) then
                   norb(ielm)=norb(ielm)+1
                   Rps0(norb(ielm),ielm)=rcnl
                   lo(norb(ielm),ielm)=i-1
                   no(norb(ielm),ielm)=j
                   hnl(j,i-1,ielm)=work(j)
                end if
             end do
             write(*,'(1x,f14.8,4f14.8)')Rps0(norb(ielm),ielm),hnl(1:2,i-1,ielm)
          end do
10        continue
          close(unit_pp)
       end do
    end if

    call MPI_BCAST(Rcloc,Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(parloc,4*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Zps,Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(norb,Nelement,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(lo,4*Nelement,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(no,4*Nelement,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Rps0,4*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(hnl,4*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    do ielm=1,Nelement
       do iorb=1,norb(ielm)
          L=lo(iorb,ielm)
          n=no(iorb,ielm)
          if ( hnl(n,L,ielm) < 0.d0 ) then
             inorm(iorb,ielm)=-1
             hnl(n,L,ielm) = abs( hnl(n,L,ielm) )
          else
             inorm(iorb,ielm)=1
          end if
       end do
    end do

    rmax = 30.d0
    MMr  = 3000
    dr   = rmax/MMr
    pi   = acos(-1.d0)
    ep   = 1.d-10

    do ielm=1,Nelement
    do iorb=1,norb(ielm)
       n=no(iorb,ielm)
       L=lo(iorb,ielm)
       gamma=sqrt(Pi)
       do i=1,L+2*n-1
          gamma=gamma*(i-0.5d0)
       end do
       const1=sqrt(2.d0)/(Rps0(iorb,ielm)**(L+2*n-0.5d0)*sqrt(gamma))
       const2=0.5d0/(Rps0(iorb,ielm)*Rps0(iorb,ielm))
       do i=1,MMr
          r=i*dr
          v=const1*r**(L+2*n-2)*exp(-r*r*const2)
          if ( abs(v) < ep ) then
             Rps(iorb,ielm)=max( Rps(iorb,ielm),r )
             exit
          end if
       end do
       if ( rank == 0 ) then
          write(*,*) ielm,iorb,n,L,Rps(iorb,ielm)
       end if
     end do ! iorb
     end do ! ielm

  END SUBROUTINE read_ps_gth


END MODULE ps_gth_module
