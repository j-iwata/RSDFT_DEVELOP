module test_force_module

  use watch_module
  use parallel_module
  use atom_module

  use force_local_sol_module, only: calc_force_local_sol
  use force_ewald_module
  use nonlocal_module
  use ps_pcc_force_module, only: calc_ps_pcc_force
  use ps_pcc_module, only: flag_pcc_0

  use force_local_mol_module
  use force_ion_mol_module
  use force_nloc2_mol_module

  implicit none

  private
  public :: test_force

  logical :: disp_switch
  integer :: MB

contains


  subroutine test_force( systype, disp_switch_in )
    implicit none
    integer,intent(in) :: systype
    logical,intent(in) :: disp_switch_in

    call write_border( 0, " test_force(start)" )

    disp_switch = disp_switch_in

    MB = sum( ir_band )

    select case(systype)
    case default
       if ( disp_switch ) write(*,'(a40,"test_force(solid)")') repeat("-",40)
       call test_force_sol
    case(1)
       if ( disp_switch ) write(*,'(a40,"test_force(mol)")') repeat("-",40)
       call test_force_mol
    end select

    call write_border( 0, " test_force(end)" )

  end subroutine test_force


  subroutine test_force_sol
    implicit none
    real(8),allocatable :: force(:,:),forcet(:,:)
    integer :: i,ii,MBD_org
    real(8) :: ttmp(2),ttt(2,0:9)

    allocate( force(3,Natom)  ) ; force=0.d0
    allocate( forcet(3,Natom) ) ; forcet=0.d0

    call watchb( ttmp, barrier='on' ); ttt=0.0d0

    call calc_force_local_sol( Natom, force )

    call watchb( ttmp, ttt(:,1), barrier='on' )

    forcet(:,:)=forcet(:,:)+force(:,:)

    call watchb( ttmp, ttt(:,5), barrier='on' )

    if ( disp_switch ) then
       write(*,*) "(flocal)"
       do i=1,Natom
          write(*,'(1x,i6,3g20.10,i6)') i,force(1:3,i),myrank
       end do
    end if

    if ( flag_pcc_0 ) then
       call watchb( ttmp, barrier='on' )
       call calc_ps_pcc_force( Natom, force )
       call watchb( ttmp, ttt(:,4), barrier='on' )
       forcet(:,:)=forcet(:,:)+force(:,:)
       if ( disp_switch ) then
          write(*,*) "(fpcc)"
          do i=1,Natom
             write(*,'(1x,i6,3g20.10,i6)') i,force(1:3,i),myrank
          end do
       end if
    end if

    call watchb( ttmp, barrier='on' )

    call calc_force_ewald(Natom,force)

    call watchb( ttmp, ttt(:,2), barrier='on' )

    forcet(:,:)=forcet(:,:)+force(:,:)

    call watchb( ttmp, ttt(:,5), barrier='on' )

    if ( disp_switch ) then
       write(*,*) "(fewald)"
       do i=1,Natom
          write(*,'(1x,i6,3g20.10,i6)') i,force(1:3,i),myrank
       end do
    end if

! --- Nonlocal Part ---

    ttt(:,3)=0.0d0

    call watchb( ttmp, barrier='on' )

    call calc_force_nonlocal(Natom,force)

    call watchb( ttmp, ttt(:,3), barrier='on' )
    if ( myrank == 0 ) then
       write(*,*) "MB_d_nl / MB =",MB_d_nl, " / ",MB
       write(*,*) "time(nloc)=",ttt(:,3)
    end if

    MBD_org = MB_d_nl

    MB_d_nl = 1

    do ii=0,10

    if ( (ii>0.and.MB_d_nl>=MB) .or. MB_d_nl>MBD_org ) then
       exit
    else if ( ii == 0 ) then
       MB_d_nl=1
    else if ( mod(MB,MB_d_nl*2) == 0 ) then
       MB_d_nl=MB_d_nl*2
    else if ( mod(MB,MB_d_nl*3) == 0 ) then
       MB_d_nl=MB_d_nl*3
    else if ( mod(MB,MB_d_nl*5) == 0 ) then
       MB_d_nl=MB_d_nl*5
    else if ( mod(MB,MB_d_nl*7) == 0 ) then
       MB_d_nl=MB_d_nl*7
    else if ( MB_d_nl >= MB ) then
       MB_d_nl=MB
    end if

    ttt(:,3)=0.0d0

    call watchb( ttmp, barrier='on' )

    call calc_force_nonlocal(Natom,force)

    call watchb( ttmp, ttt(:,3), barrier='on' )

    if ( myrank == 0 ) then
       write(*,*) "MB_d_nl / MB =",MB_d_nl, " / ",MB
       write(*,*) "time(nloc)=",ttt(:,3)
    end if

    end do ! ii

    call watchb( ttmp, barrier='on' )

    forcet(:,:)=forcet(:,:)+force(:,:)

    call watchb( ttmp, ttt(:,5), barrier='on' )

    if ( disp_switch ) then
       write(*,*) "(fnloc)"
       do i=1,Natom
          write(*,'(1x,i6,3g20.10,i6)') i,force(1:3,i),myrank
       end do
    end if

    MB_d_nl = MBD_org

! ---

    if ( disp_switch ) then
       write(*,*) "ftot"
       do i=1,Natom
          write(*,'(1x,i6,3g20.10,i6)') i,forcet(1:3,i),myrank
       end do
    end if

    if ( myrank == 0 ) then
       write(*,*) "time(loc) =",ttt(:,1)
       write(*,*) "time(ewld)=",ttt(:,2)
       if ( flag_pcc_0 ) write(*,*) "time(pcc) =",ttt(:,4)
       write(*,*) "time(othr)=",ttt(:,5)
    end if

    deallocate( forcet )
    deallocate( force  )

  END SUBROUTINE test_force_sol


  SUBROUTINE test_force_mol
    implicit none
    real(8),allocatable :: force(:,:),forcet(:,:)
    integer :: i

    allocate( force(3,Natom)  ) ; force=0.d0
    allocate( forcet(3,Natom) ) ; forcet=0.d0

    call watcht(disp_switch,"floc",0)
    call calc_force_local_mol(force)
    call watcht(disp_switch,"floc",1)
    forcet(:,:)=forcet(:,:)+force(:,:)

    if ( disp_switch ) then
       do i=1,Natom
          write(*,'(1x,i6,3g20.10,i6)') i,force(1:3,i),myrank
       end do
    end if

    call watcht(disp_switch,"fewl",0)
    call calc_force_ion_mol(force)
    call watcht(disp_switch,"fewl",1)
    forcet(:,:)=forcet(:,:)+force(:,:)

    if ( disp_switch ) then
       do i=1,Natom
          write(*,'(1x,i6,3g20.10,i6)') i,force(1:3,i),myrank
       end do
    end if

    call watcht(disp_switch,"fnlc",0)
    call calc_force_nloc2_mol(force)
    call watcht(disp_switch,"fnlc",1)
    forcet(:,:)=forcet(:,:)+force(:,:)

    if ( disp_switch ) then
       do i=1,Natom
          write(*,'(1x,i6,3g20.10,i6)') i,force(1:3,i),myrank
       end do
    end if

    if ( disp_switch ) then
       write(*,*) "ftot"
       do i=1,Natom
          write(*,'(1x,i6,3g20.10,i6)') i,forcet(1:3,i),myrank
       end do
    end if

    deallocate( forcet )
    deallocate( force  )

  END SUBROUTINE test_force_mol


END MODULE test_force_module
