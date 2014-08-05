MODULE ps_gth_module

  use atom_module, only: Nelement
  use pseudopot_module

  implicit none

  PRIVATE
  PUBLIC :: read_ps_gth,hnl,Rcloc,Rps0

  real(8),allocatable :: Rps0(:,:)
  real(8),allocatable :: Rcloc(:),hnl(:,:,:),knl(:,:,:)

  integer,parameter :: lrefmax=3
  integer,parameter :: unit_pp=34

CONTAINS


  SUBROUTINE read_ps_gth(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: i,j,ielm,ierr
    integer :: MMr,iorb,L,n
    character(2) :: name
    character(30) :: cbuf
    real(8) :: znuc,rcnl,work1(3),work2(3)
    real(8) :: pi,rmax,dr,ep,gamma,const1,const2,v,r
    include 'mpif.h'

    allocate( Rcloc(Nelement)        ) ; Rcloc=0.0d0
    allocate( parloc(1:4,Nelement)   ) ; parloc=0.0d0
    allocate( Zps(Nelement)          ) ; Zps=0.0d0
    allocate( norb(Nelement)         ) ; norb=0
    allocate( lo(6,Nelement)         ) ; lo=0
    allocate( no(6,Nelement)         ) ; no=0
    allocate( Rps(6,Nelement)        ) ; Rps=0.0d0
    allocate( Rps0(6,Nelement)       ) ; Rps0=0.0d0
    allocate( hnl(3,0:2,Nelement)    ) ; hnl=0.0d0
    allocate( knl(3,1:2,Nelement)    ) ; knl=0.0d0
    allocate( inorm(6,Nelement)      ) ; inorm=0
    allocate( hnml(3,3,0:2,Nelement) ) ; hnml=0.0d0
    allocate( knml(3,3,1:2,Nelement) ) ; knml=0.0d0

    if ( rank == 0 ) then

       write(*,'(a60," read_ps_gth")') repeat("-",60)

       do ielm=1,Nelement

          open(unit_pp,file=file_ps(ielm),status='old')

          read(unit_pp,'(a)') cbuf

          if ( index( cbuf, "KRACK" ) /= 0 ) then

             call Read_KrackFormat( znuc, Zps(ielm), Rcloc(ielm) &
                  , parloc(1,ielm), Rps0(1,ielm), lo(1,ielm), no(1,ielm) &
                  , hnml(1,1,0,ielm), norb(ielm) )

          else

             backspace(unit_pp)

             read(unit_pp,*) name,znuc,Zps(ielm)
             read(unit_pp,*) Rcloc(ielm),parloc(1:4,ielm)
             write(*,*) name,znuc,Zps(ielm)
             write(*,'(1x,f14.8,4f14.8)') Rcloc(ielm),parloc(1:4,ielm)
             do i=1,lrefmax
                work1(1:3)=0.0d0
                work2(1:3)=0.0d0
                read(unit_pp,*,END=10) rcnl,work1(1:3)
                if ( pselect == 5 .and. i > 1 ) then
                   read(unit_pp,*,END=10) work2(1:3)
                end if
                do j=1,3
                   if ( work1(j) /= 0.0d0 ) then
                      norb(ielm) = norb(ielm) + 1
                      Rps0(norb(ielm),ielm) = rcnl
                      lo(norb(ielm),ielm)   = i-1
                      no(norb(ielm),ielm)   = j
                      hnl(j,i-1,ielm)       = work1(j)
                   end if
                   if ( work2(j) /= 0.0d0 ) then
                      knl(j,i-1,ielm)=work2(j)
                   end if
                end do
                write(*,'(1x,f14.8,4f14.8)') &
                     Rps0(norb(ielm),ielm),hnl(1:3,i-1,ielm)
                if ( i > 1 ) write(*,'(1x,14x,4f14.8)') knl(1:3,i-1,ielm)
             end do ! i
10           continue

             call calc_nondiagonal( hnl(1,0,ielm), knl(1,1,ielm) &
                  , hnml(1,1,0,ielm), knml(1,1,1,ielm) )

          end if ! [ Format ]

          close(unit_pp)

       end do ! ielm

    end if ! [ rank == 0 ]

    call MPI_BCAST(Rcloc,Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(parloc,4*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Zps,Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(norb,Nelement,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(lo,6*Nelement,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(no,6*Nelement,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Rps0,6*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(hnl,9*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(knl,6*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(hnml,27*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(knml,18*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    do ielm=1,Nelement
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
    end do ! ielm

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


  SUBROUTINE calc_nondiagonal( hnl, knl, hnml, knml )
    implicit none
    real(8),intent(IN)    :: hnl(3,0:2), knl(3,1:2)
    real(8),intent(INOUT) :: hnml(3,3,0:2), knml(3,3,1:2)
    hnml(1,1,0)= hnl(1,0)
    hnml(2,2,0)= hnl(2,0)
    hnml(3,3,0)= hnl(3,0)
    hnml(1,2,0)=-0.5d0*sqrt(3.d0/5.d0)*hnml(2,2,0)
    hnml(1,3,0)= 0.5d0*sqrt(5.d0/21.d0)*hnml(3,3,0)
    hnml(2,3,0)=-0.5d0*sqrt(100.d0/63.d0)*hnml(3,3,0)
    hnml(1,1,1)= hnl(1,1)
    hnml(2,2,1)= hnl(2,1)
    hnml(3,3,1)= hnl(3,1)
    hnml(1,2,1)=-0.5d0*sqrt(5.d0/7.d0)*hnml(2,2,1)
    hnml(1,3,1)=sqrt(35.d0/11.d0)/6.d0*hnml(3,3,1)
    hnml(2,3,1)=-14.d0/sqrt(11.d0)/6.d0*hnml(3,3,1)
    hnml(1,1,2)= hnl(1,2)
    hnml(2,2,2)= hnl(2,2)
    hnml(3,3,2)= hnl(3,2)
    hnml(1,2,2)=-0.5d0*sqrt(7.d0/9.d0)*hnml(2,2,2)
    hnml(1,3,2)= 0.5d0*sqrt(63.d0/143.d0)*hnml(3,3,2)
    hnml(2,3,2)=-0.5d0*18.d0/sqrt(143.d0)*hnml(3,3,2)
!
    knml(1,1,1)= knl(1,1)
    knml(2,2,1)= knl(2,1)
    knml(3,3,1)= knl(3,1)
    knml(1,2,1)=-0.5d0*sqrt(5.d0/7.d0)*knml(2,2,1)
    knml(1,3,1)=sqrt(35.d0/11.d0)/6.d0*knml(3,3,1)
    knml(2,3,1)=-14.d0/sqrt(11.d0)/6.d0*knml(3,3,1)
    knml(1,1,2)= knl(1,2)
    knml(2,2,2)= knl(2,2)
    knml(3,3,2)= knl(3,2)
    knml(1,2,2)=-0.5d0*sqrt(7.d0/9.d0)*knml(2,2,2)
    knml(1,3,2)= 0.5d0*sqrt(63.d0/143.d0)*knml(3,3,2)
    knml(2,3,2)=-0.5d0*18.d0/sqrt(143.d0)*knml(3,3,2)
  END SUBROUTINE calc_nondiagonal


  SUBROUTINE Read_KrackFormat( znuc, Zps, Rcloc, parloc, Rps0, lo, no &
       , hnml, norb )
    implicit none
    character(2) :: name
    real(8) :: znuc,Zps,hnml(3,3,0:2)
    real(8) :: Rcloc,parloc(4),Rps0(6)
    integer :: lo(6), no(6),norb
    integer :: i, j, k
    real(8) :: work(3), rcnl
    character(30) :: cbuf

    read(unit_pp,*) name,znuc,Zps
    read(unit_pp,*) Rcloc,parloc(1:4)
    write(*,*) name,znuc,Zps
    write(*,'(1x,f14.8,4f14.8)') Rcloc, parloc(1:4)

    norb = 0

    do i=1,lrefmax

       do j=1,3

          work(:) = 0.0d0

          read(unit_pp,*,END=10) cbuf, rcnl, work(1:3)

          if ( cbuf(1:1) == "s" .or. cbuf(1:1) == "S" .or. &
               cbuf(1:1) == "p" .or. cbuf(1:1) == "P" .or. &
               cbuf(1:1) == "d" .or. cbuf(1:1) == "D" ) then

             norb       = norb + 1
             Rps0(norb) = rcnl
             lo(norb)   = i-1
             no(norb)   = j

             do k=1,3
                if ( j+k-1 > 3 ) exit
                if ( work(k) /= 0.0d0 ) hnml(j,j+k-1,i-1) = work(k)
             end do

             write(*,'(1x,a1,2x,2i4,f14.8)') cbuf, norb, lo(norb), Rps0(norb)

             write(*,'(12x,3f14.8)') ( hnml(j,k,i-1), k=1,3 )

          else

             norb       = norb + 1
             Rps0(norb) = rcnl
             lo(norb)   = i-1
             no(norb)   = j

             backspace(unit_pp)
             read(unit_pp,*) work(1:3)

             do k=1,3
                if ( j+k-1 > 3 ) exit
                if ( work(k) /= 0.0d0 ) hnml(j,j+k-1,i-1) = work(k)
             end do

             write(*,'(12x,3f14.8)') ( hnml(j,k,i-1), k=1,3 )

          end if

          read(unit_pp,*,END=10) cbuf
          backspace(unit_pp)
          if ( cbuf(1:1) == "s" .or. cbuf(1:1) == "S" .or. &
               cbuf(1:1) == "p" .or. cbuf(1:1) == "P" .or. &
               cbuf(1:1) == "d" .or. cbuf(1:1) == "D" ) exit

       end do ! j

    end do ! i
10  continue

  END SUBROUTINE Read_KrackFormat


END MODULE ps_gth_module
