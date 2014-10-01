MODULE force_module

  use ps_local_module
  use ps_nloc2_module
  use ps_pcc_module
  use force_ewald_module
  use watch_module
  use parallel_module, only: disp_switch_parallel
  use vdw_grimme_module, only: calc_F_vdw_grimme
  use atom_module, only: aa_atom, md_atom

  implicit none

  PRIVATE
  PUBLIC :: init_force, calc_force

  integer :: Ntim
  real(8),allocatable :: tim(:,:,:)

  include 'mpif.h'

CONTAINS


  SUBROUTINE init_force( rank, unit )
    implicit none
    integer :: rank, unit
    if ( disp_switch_parallel ) &
         write(*,'(1x,a60," init_force")') repeat("-",60)
    Ntim=0
    Ntim=maxval( md_atom )
    if ( rank == 0 ) write(*,*) "Ntim=",Ntim
    if ( Ntim <= 0 ) return
    allocate( tim(3,3,Ntim) ) ; tim=0.0d0
    tim(1,1,1) = 1.0d0
    tim(2,2,1) = 1.0d0
    tim(3,3,1) = 1.0d0
    call read_force( rank, unit )
  END SUBROUTINE init_force


  SUBROUTINE read_force( rank, unit )
    implicit none
    integer,intent(IN) :: rank, unit
    integer :: i,j,k,kr,ierr,mchk
    character(3) :: cbuf,ckey
    mchk=0
    if ( rank == 0 ) then
       write(*,'(1x,a40," read_force")') repeat("-",40)
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital( cbuf, ckey )
          if ( ckey == "TIM" ) then
             do j=1,Ntim
                do k=1,3
                   read(unit,*) ( tim(kr,k,j), kr=1,3 )
                   mchk=mchk+1
                end do ! k
             end do ! j
             exit
          end if
       end do ! i
999    continue
       if ( mchk /= 3*Ntim ) then
          write(*,*) "WARNING: Insufficient data for tim-matrix !!!"
          write(*,*) "mchk, 3*Ntim =",mchk,3*Ntim
       end if
       do j=1,Ntim
          write(*,'(1x,"tim(",i1,")",2x,3f10.5)') j,(tim(kr,1,j),kr=1,3)
          do k=2,3
             write(*,'(1x,6x,2x,3f10.5)') (tim(kr,k,j),kr=1,3)
          end do
       end do
    end if
    call mpi_bcast( tim, 9*Ntim, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
  END SUBROUTINE read_force


  SUBROUTINE calc_force( MI, force )
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force(3,MI)
    real(8),allocatable :: work(:,:)
    real(8) :: ctt(0:4),ett(0:4)

    force(:,:) = 0.d0

    ctt(:)=0.d0
    ett(:)=0.d0

    allocate( work(3,MI) )

    call watch(ctt(0),ett(0))

#ifdef _FFTE_
    call calc_force_ps_local_ffte(MI,work)
#else
    call calc_force_ps_local(MI,work)
#endif
    force = force + work

    if ( flag_pcc_0 ) then
       call calc_force_ps_pcc(MI,work)
       force = force + work
    end if

    call watch(ctt(1),ett(1))

    call calc_force_ps_nloc2(MI,work)
    force = force + work

    call watch(ctt(2),ett(2))

    call calc_force_ewald(MI,work)
    force = force + work

    call watch(ctt(3),ett(3))

    call calc_F_vdw_grimme( MI, aa_atom, force )

    call watch(ctt(4),ett(4))

! --- constraint & symmetry ---

!    call_symforce
    if ( Ntim > 0 ) call constraint_force( MI, force )

    deallocate( work )

    if ( disp_switch_parallel ) then
       write(*,*) "time(force1)",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(force2)",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(force3)",ctt(3)-ctt(2),ett(3)-ett(2)
       write(*,*) "time(force4)",ctt(4)-ctt(3),ett(4)-ett(3)
    end if

  END SUBROUTINE calc_force


  SUBROUTINE constraint_force( MI, force )
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(INOUT) :: force(3,MI)
    integer :: i,j
    do i=1,MI
       j=md_atom(i)
       if ( j > 0 ) then
          force(:,i) = matmul( tim(:,:,j), force(:,i) )
       end if
    end do
  END SUBROUTINE constraint_force


END MODULE force_module
