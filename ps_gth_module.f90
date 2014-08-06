MODULE ps_gth_module

  implicit none

  PRIVATE
  PUBLIC :: read_ps_gth, ps_gth

  integer :: unit_pp
  integer,parameter :: lrefmax=6

  TYPE gth
     integer :: norb
     integer :: inorm(lrefmax)
     integer :: lo(lrefmax)
     integer :: no(lrefmax)
     real(8) :: znuc
     real(8) :: Rc(lrefmax)
     real(8) :: hnl(3,0:2)
     real(8) :: knl(3,1:2)
     real(8) :: hnml(3,3,0:2)
     real(8) :: knml(3,3,1:2)
     real(8) :: Dij(lrefmax,lrefmax)
     real(8) :: parloc(4)
     real(8) :: Rcloc
  END TYPE gth

  type(gth) :: ps_gth

CONTAINS


  SUBROUTINE read_ps_gth( unit_ps, ippform )
    implicit none
    integer,intent(IN) :: unit_ps, ippform
    integer :: i,j
    integer :: MMr,iorb,L,n
    character(2) :: name
    character(30) :: cbuf
    real(8) :: znuc,rcnl,work1(3),work2(3)
    integer :: norb,lo(lrefmax),no(lrefmax),inorm(lrefmax)
    real(8) :: Zps,parloc(4),Rcloc,Rps0(lrefmax)
    real(8) :: hnml(3,3,0:2), knml(3,3,2), hnl(3,0:2), knl(3,2)
    logical :: iflag_hgh

    unit_pp = unit_ps

    ps_gth%norb     =0
    ps_gth%znuc     =0.0d0
    ps_gth%Rc(:)    =0.0d0
    ps_gth%inorm(:) =0
    ps_gth%lo(:)    =0
    ps_gth%no(:)    =0
    ps_gth%Dij(:,:) =0.0d0
    ps_gth%parloc(:)=0.0d0
    ps_gth%Rcloc    =0.0d0
    ps_gth%hnl(:,:) =0.0d0
    ps_gth%knl(:,:) =0.0d0
    ps_gth%hnml(:,:,:) =0.0d0
    ps_gth%knml(:,:,:) =0.0d0

    write(*,'(a60," read_ps_gth")') repeat("-",60)

    read(unit_pp,'(a)') cbuf

    if ( index( cbuf, "KRACK" ) /= 0 ) then

       call Read_KrackFormat( znuc,Zps,Rcloc,parloc,Rps0,lo,no,hnml,norb )

    else ! [ Original Format ]

       backspace(unit_pp)

       read(unit_pp,*) name,znuc,Zps
       read(unit_pp,*) Rcloc, parloc(1:4)
       write(*,*) name,znuc,Zps
       write(*,'(1x,f14.8,4f14.8)') Rcloc,parloc(1:4)

       j=0
       do i=1,lrefmax
          cbuf=""
          read(unit_pp,'(a)',END=99) cbuf
          if ( cbuf /= "" ) j=j+1
       end do
99     continue
       if ( j <= 2 ) then
          iflag_hgh = .false.
          write(*,*) "GTH format",j
       else
          iflag_hgh=.true.
          write(*,*) "HGH format",j
       end if
       do j=1,i
          backspace(unit_pp)
       end do

       norb=0
       do i=1,lrefmax

          work1(1:3)=0.0d0
          work2(1:3)=0.0d0
          read(unit_pp,*,END=10) rcnl,work1(1:3)

          if ( iflag_hgh .and. i > 1 ) read(unit_pp,*,END=10) work2(1:3)

          do j=1,3
             if ( work1(j) /= 0.0d0 ) then
                norb       = norb + 1
                Rps0(norb) = rcnl
                lo(norb)   = i-1
                no(norb)   = j
                hnl(j,i-1) = work1(j)
             end if
             if ( work2(j) /= 0.0d0 ) knl(j,i-1)=work2(j)
          end do

          write(*,'(1x,f14.8,4f14.8)') Rps0(norb),hnl(1:3,i-1)
          if ( i > 1 ) write(*,'(1x,14x,4f14.8)') knl(1:3,i-1)

       end do ! i
10     continue

       if ( iflag_hgh ) call calc_nondiagonal( hnl, knl, hnml, knml )

    end if ! [ Format ]

    if ( iflag_hgh ) then

       do L=0,maxval(lo)
          do j=1,3
          do i=1,j-1
             hnml(j,i,L)=hnml(i,j,L)
             if ( L >= 1 ) knml(j,i,L)=knml(i,j,L)
          end do
          end do
       end do

    else

       do iorb=1,norb
          L=lo(iorb)
          n=no(iorb)
          if ( hnl(n,L) < 0.0d0 ) then
             inorm(iorb)=-1
             hnl(n,L) = abs( hnl(n,L) )
          else
             inorm(iorb)=1
          end if
       end do

    end if

    ps_gth%znuc        = Zps
    ps_gth%parloc(:)   = parloc(:)
    ps_gth%Rcloc       = Rcloc
    ps_gth%Rc(:)       = Rps0(:)
    ps_gth%lo(:)       = lo(:)
    ps_gth%no(:)       = no(:)
    ps_gth%norb        = norb
    ps_gth%hnl(:,:)    = hnl(:,:)
    ps_gth%knl(:,:)    = knl(:,:)
    ps_gth%hnml(:,:,:) = hnml(:,:,:)
    ps_gth%knml(:,:,:) = knml(:,:,:)
    ps_gth%inorm(:)    = inorm(:)

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


  SUBROUTINE Read_KrackFormat( znuc,Zps,Rcloc,parloc,Rps0,lo,no,hnml,norb )
    implicit none
    character(2) :: name
    real(8) :: znuc,Zps,hnml(3,3,0:2)
    real(8) :: Rcloc,parloc(4),Rps0(lrefmax)
    integer :: lo(lrefmax), no(lrefmax),norb
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
