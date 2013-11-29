MODULE ps_gth_module

  use atom_module, only: Nelement
  use pseudopot_module

  implicit none

  PRIVATE
  PUBLIC :: read_ps_gth,hnl,hnml,knml,no,Rcloc,Rps0

  real(8),allocatable :: Rps0(:,:)
  integer,allocatable :: no(:,:)
  real(8),allocatable :: Rcloc(:),hnl(:,:,:),knl(:,:,:)
  real(8),allocatable :: hnml(:,:,:,:),knml(:,:,:,:)

CONTAINS

  SUBROUTINE read_ps_gth(rank)
    implicit none
    integer,intent(IN) :: rank
    integer,parameter :: lrefmax=3
    integer,parameter :: unit_pp=34
    integer :: i,j,ielm,ierr
    character(2) :: name
    real(8) :: znuc,rcnl,work(3),work2(3)
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
    allocate( hnl(3,0:2,Nelement) ) ; hnl=0.d0
    allocate( knl(3,1:2,Nelement) ) ; knl=0.d0
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
             work(1:3)=0.d0
             work2(1:3)=0.d0
             read(unit_pp,*,END=10) rcnl,work(1:3)
             if ( pselect == 5 .and. i > 1 ) read(unit_pp,*,END=10) work2(1:3)
             do j=1,3
                if ( work(j) /= 0.d0 ) then
                   norb(ielm)=norb(ielm)+1
                   Rps0(norb(ielm),ielm)=rcnl
                   lo(norb(ielm),ielm)=i-1
                   no(norb(ielm),ielm)=j
                   hnl(j,i-1,ielm)=work(j)
                end if
                if ( work2(j) /= 0.d0 ) then
                   knl(j,i-1,ielm)=work2(j)
                end if
             end do
             write(*,'(1x,f14.8,4f14.8)')Rps0(norb(ielm),ielm),hnl(1:3,i-1,ielm)
             if ( i > 1 ) write(*,'(1x,14x,4f14.8)') knl(1:3,i-1,ielm)
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
    call MPI_BCAST(hnl,9*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(knl,6*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    allocate( hnml(3,3,0:2,Nelement) ) ; hnml=0.0d0
    allocate( knml(3,3,1:2,Nelement) ) ; knml=0.0d0

    do ielm=1,Nelement
       hnml(1,1,0,ielm)= hnl(1,0,ielm)
       hnml(2,2,0,ielm)= hnl(2,0,ielm)
       hnml(3,3,0,ielm)= hnl(3,0,ielm)
       hnml(1,2,0,ielm)=-0.5d0*sqrt(3.d0/5.d0)*hnml(2,2,0,ielm)
       hnml(1,3,0,ielm)= 0.5d0*sqrt(5.d0/21.d0)*hnml(3,3,0,ielm)
       hnml(2,3,0,ielm)=-0.5d0*sqrt(100.d0/63.d0)*hnml(3,3,0,ielm)
       hnml(1,1,1,ielm)= hnl(1,1,ielm)
       hnml(2,2,1,ielm)= hnl(2,1,ielm)
       hnml(3,3,1,ielm)= hnl(3,1,ielm)
       hnml(1,2,1,ielm)=-0.5d0*sqrt(5.d0/7.d0)*hnml(2,2,1,ielm)
       hnml(1,3,1,ielm)=sqrt(35.d0/11.d0)/6.d0*hnml(3,3,1,ielm)
       hnml(2,3,1,ielm)=-14.d0/sqrt(11.d0)/6.d0*hnml(3,3,1,ielm)
       hnml(1,1,2,ielm)= hnl(1,2,ielm)
       hnml(2,2,2,ielm)= hnl(2,2,ielm)
       hnml(3,3,2,ielm)= hnl(3,2,ielm)
       hnml(1,2,2,ielm)=-0.5d0*sqrt(7.d0/9.d0)*hnml(2,2,2,ielm)
       hnml(1,3,2,ielm)= 0.5d0*sqrt(63.d0/143.d0)*hnml(3,3,2,ielm)
       hnml(2,3,2,ielm)=-0.5d0*18.d0/sqrt(143.d0)*hnml(3,3,2,ielm)
!
       knml(1,1,1,ielm)= knl(1,1,ielm)
       knml(2,2,1,ielm)= knl(2,1,ielm)
       knml(3,3,1,ielm)= knl(3,1,ielm)
       knml(1,2,1,ielm)=-0.5d0*sqrt(5.d0/7.d0)*knml(2,2,1,ielm)
       knml(1,3,1,ielm)=sqrt(35.d0/11.d0)/6.d0*knml(3,3,1,ielm)
       knml(2,3,1,ielm)=-14.d0/sqrt(11.d0)/6.d0*knml(3,3,1,ielm)
       knml(1,1,2,ielm)= knl(1,2,ielm)
       knml(2,2,2,ielm)= knl(2,2,ielm)
       knml(3,3,2,ielm)= knl(3,2,ielm)
       knml(1,2,2,ielm)=-0.5d0*sqrt(7.d0/9.d0)*knml(2,2,2,ielm)
       knml(1,3,2,ielm)= 0.5d0*sqrt(63.d0/143.d0)*knml(3,3,2,ielm)
       knml(2,3,2,ielm)=-0.5d0*18.d0/sqrt(143.d0)*knml(3,3,2,ielm)
       do L=0,lrefmax-1
          do j=1,3
          do i=1,j-1
             hnml(j,i,L,ielm)=hnml(i,j,L,ielm)
             if ( L >= 1 ) knml(j,i,L,ielm)=knml(i,j,L,ielm)
          end do
          end do
       end do
    end do ! ielm

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
