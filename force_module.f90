MODULE force_module

  use ps_local_module
  use ps_nloc2_module
  use ps_pcc_module
  use ewald_module
  use watch_module
  use parallel_module, only: disp_switch_parallel,myrank

  use atom_module, only: opt_constrain

  implicit none

  PRIVATE
  PUBLIC :: calc_force
  logical :: isConstraint=.true.,firstConstraint=.true.

CONTAINS

  SUBROUTINE calc_force(MI,force)
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force(3,MI)
    real(8),allocatable :: work(:,:)
integer :: a
real(8),allocatable :: work1(:,:),work2(:,:),work3(:,:)
    real(8) :: ctt(0:3),ett(0:3)

    force(:,:) = 0.d0

    ctt(:)=0.d0
    ett(:)=0.d0

    allocate( work(3,MI) )
allocate( work1(3,MI) ) ; work1=0.d0
allocate( work2(3,MI) ) ; work2=0.d0
allocate( work3(3,MI) ) ; work3=0.d0

    call watch(ctt(0),ett(0))

#ifdef _FFTE_
    call calc_force_ps_local_ffte(MI,work)
#else
!    call calc_force_ps_local(MI,work)
call calc_force_ps_local(MI,work1)
#endif
!    force = force + work

    if ( flag_pcc_0 ) then
       call calc_force_ps_pcc(MI,work)
       force = force + work
    end if

    call watch(ctt(1),ett(1))

!    call calc_force_ps_nloc2(MI,work)
call calc_force_ps_nloc2(MI,work2)
!    force = force + work

    call watch(ctt(2),ett(2))

!    call calc_force_ewald(MI,work)
call calc_force_ewald(MI,work3)
!    force = force + work

    call watch(ctt(3),ett(3))

do a=1,MI
if ( myrank==0 ) write(200,'(I5,A9,3g20.7)') a,'local',work1(1:3,a)
if ( myrank==0 ) write(200,'(I5,A9,3g20.7)') a,' nloc',work2(1:3,a)
if ( myrank==0 ) write(200,'(I5,A9,3g20.7)') a,'ewald',work3(1:3,a)
enddo
    force=work1+work2+work3
if ( myrank==0 ) write(200,*) '--------------------'
do a=1,MI
if ( myrank==0 ) write(200,'(I5,A9,3g20.7)') a,'total',force(1:3,a)
enddo
if ( myrank==0 ) write(200,*) '--------------------'

! --- constraint & symmetry ---

!    call_symforce
    if (isConstraint) then
      call constraint(MI,force)
    endif

    deallocate( work )
deallocate( work1 )
deallocate( work2 )
deallocate( work3 )

    if ( disp_switch_parallel ) then
       write(*,*) "time(force1)",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(force2)",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(force3)",ctt(3)-ctt(2),ett(3)-ett(2)
    end if

  END SUBROUTINE calc_force

  SUBROUTINE constraint(MI,force_)
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(INOUT) :: force_(3,MI)
    integer :: m
    integer :: ia
    
    if (firstConstraint) then
      m=maxval(opt_constrain)
      if (m==0) then
        isConstraint=.false.
        return
      endif
      firstConstraint=.false.
    endif

    do ia=1,MI
      select case (opt_constrain(ia))
      case default
        stop
      case (0)
      case (3000)
        force_(1:3,ia)=0.d0
      case (1011)
        force_(1,ia)=0.d0
      case (1101)
        force_(2,ia)=0.d0
      case (1110)
        force_(3,ia)=0.d0
      case (2001)
        force_(1,ia)=0.d0
        force_(2,ia)=0.d0
      case (2010)
        force_(1,ia)=0.d0
        force_(3,ia)=0.d0
      case (2100)
        force_(2,ia)=0.d0
        force_(3,ia)=0.d0
      end select
    enddo
    return
  END SUBROUTINE constraint

END MODULE force_module
