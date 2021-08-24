module io_read_wf_simple_module

  use parallel_module
  use wf_module
  use rgrid_module
  use rgrid_mol_module, only: LL
  use rsdft_mpi_module
  use watch_module, only: watchb
  use io_read_wf_general_module, only: read_wf_general

  implicit none

  private
  public :: read_wf_simple

contains

  subroutine read_wf_simple( file_wf2, SYStype, b, d_wf_out, z_wf_out, occ_out, kbb_out )
    use grid_module, only: construct_map_1d_to_3d_grid
    use rgrid_mol_module, only: construct_map_1d_to_3d_rgrid_mol
    implicit none
    character(*),intent(in) :: file_wf2
    integer,intent(in) :: SYStype
    type(wfrange),optional,intent(INOUT) :: b
    real(8),allocatable,optional,intent(INOUT) :: d_wf_out(:,:,:,:)
    complex(8),allocatable,optional,intent(INOUT) :: z_wf_out(:,:,:,:)
    real(8),allocatable,optional,intent(INOUT) :: occ_out(:,:,:)
    real(8),allocatable,optional,intent(INOUT) :: kbb_out(:,:)
    include 'mpif.h'
    integer :: ierr,istatus(MPI_STATUS_SIZE,123),irank,jrank
    integer :: itmp(21),n,k,s,i,i1,i2,i3,j1,j2,j3,mx,my,mz,m1,m2,m3
    integer :: ML_tmp,ML1_tmp,ML2_tmp,ML3_tmp,MB_tmp,MB1_tmp,MB2_tmp
    integer :: ML,ML1,ML2,ML3,MB,n1,n2,ML0,MBZ_tmp,MSP_tmp,MMBZ_tmp
    integer :: MB_0,MB_1,MBZ,MBZ_0,MBZ_1,MSP,MSP_0,MSP_1,nn
    integer :: IO_ctrl0, OC, TYPE_WF
    integer,allocatable :: LL2(:,:),LL_tmp(:,:),ir(:),id(:)
    real(8) :: aa0(3,3),bb0(3,3)
    real(8),allocatable :: kbb0(:,:),occ_tmp(:,:,:)
    logical :: flag_related, flag_newformat, disp_on
    logical :: read_wf_is_complex, present_wf_is_complex
    character(5) :: cmyrank
    character(32) :: file_wf_split
    complex(8),allocatable :: utmp(:),utmp3(:,:,:)
    complex(4),allocatable :: utmpSP(:)
    real(8),allocatable :: dtmp(:),dtmp3(:,:,:)
    real(4),allocatable :: dtmpSP(:)
    type(para) :: pinfo0, pinfo1
    real(8) :: ttmp(2),tt(2)

    call write_border( 0, " read_wf_simple(start)" )

    call watchb( ttmp, barrier='on' ); tt=0.0d0

! ---

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    MB  = MB_WF
    n1  = ML_0_WF
    n2  = ML_1_WF
    ML0 = n2 - n1 + 1

    MB_0  = MB_0_WF
    MB_1  = MB_1_WF
    MBZ   = MK_WF
    MBZ_0 = MK_0_WF
    MBZ_1 = MK_1_WF
    MSP   = MS_WF
    MSP_0 = MS_0_WF
    MSP_1 = MS_1_WF

! ---

    if ( myrank == 0 ) then

      itmp(:)=0

      open(3,file=file_wf2,form='unformatted')
      read(3) itmp(1)

      flag_newformat = .false.
      if ( itmp(1) < 0 ) flag_newformat = .true.
      write(*,*) "flag_newformat=",flag_newformat

      if ( flag_newformat ) then
        read(3) itmp(1:4)    ! ML,ML1,ML2,ML3
        read(3) itmp(5:7)    ! MB,MB1,MB2
        read(3) itmp(8:10)   ! MBZ,MSP,MMBZ
        read(3) itmp(11:13)  ! IO_Ctrl0,OC,TYPE_WF
      else
        rewind 3
        read(3) itmp(1:4)
        read(3) itmp(5:7)
        itmp(8)=MBZ
        itmp(9)=MSP
      end if

    end if

    call MPI_Bcast(flag_newformat,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(itmp,size(itmp),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    ML_tmp  = itmp(1)
    ML1_tmp = itmp(2)
    ML2_tmp = itmp(3)
    ML3_tmp = itmp(4)
    MB_tmp  = itmp(5)
    MB1_tmp = itmp(6)
    MB2_tmp = itmp(7)
    MBZ_tmp = itmp(8)
    MSP_tmp = itmp(9)
    MMBZ_tmp= itmp(10)
    IO_ctrl0= itmp(11)
    OC      = itmp(12)
    TYPE_WF = itmp(13)

    call check_disp_switch( disp_on, 0 )
    if ( disp_on ) then
      write(*,*) "ML ,ML_tmp                =",ML ,ML_tmp
      write(*,*) "ML1    ,ML2    ,ML3       =",ML1,ML2,ML3
      write(*,*) "ML1_tmp,ML2_tmp,ML3_tmp   =",ML1_tmp,ML2_tmp,ML3_tmp
      write(*,*) "MB ,MB_tmp                =",MB ,MB_tmp
      write(*,*) "MB1_tmp, MB2_tmp          =",MB1_tmp,MB2_tmp
      if ( flag_newformat ) then
        write(*,*) "MBZ_tmp, MSP_tmp,MMBZ_tmp =",MBZ_tmp,MSP_tmp,MMBZ_tmp
        write(*,*) "IO_ctrl0, OC, TYPE_WF     =",IO_ctrl0, OC, TYPE_WF
        write(*,*) "(new format)"
      else
        write(*,*) "(old format)"
      end if
    end if
    if ( ML_tmp /= ML .or. ML1_tmp /= ML1 .or. &
         ML2_tmp /= ML2 .or. ML3_tmp /= ML3 ) call stop_program("stop@read_wf_simple")

    if ( myrank == 0 ) then
      allocate( LL_tmp(3,ML) ); LL_tmp=0
      read(3) LL_tmp(1:3,1:ML)
      write(*,*) "minval(LL_tmp),maxval(LL_tmp)"
      write(*,*) minval(LL_tmp(1,1:ML) ),maxval( LL_tmp(1,1:ML))
      write(*,*) minval(LL_tmp(2,1:ML) ),maxval( LL_tmp(2,1:ML))
      write(*,*) minval(LL_tmp(3,1:ML) ),maxval( LL_tmp(3,1:ML))
    end if

    if ( flag_newformat ) then

      if ( present(b) ) then
        if ( b%MB == 0 ) then
          b%MB  = MB
          b%MB0 = 1
          b%MB1 = MB
        end if
        if ( b%MK == 0 ) then
          b%MK  = MBZ_tmp
          b%MK0 = 1
          b%MK1 = MBZ_tmp
          b%MMK = MMBZ_tmp
        end if
        if ( b%MS == 0 ) then
          b%MS  = MSP_tmp
          b%MS0 = 1
          b%MS1 = MSP_tmp
        end if
        if ( present(d_wf_out) ) call allocate_b_dwf( b, d_wf_out )
        if ( present(z_wf_out) ) call allocate_b_zwf( b, z_wf_out )
        if ( present(occ_out) ) call allocate_b_occ( b, occ_out )
        if ( present(kbb_out) ) then
          allocate( kbb_out(3,b%MK0:b%MK1) ); kbb_out=0.0d0
        end if
      end if

    end if

! ---

    allocate( occ_tmp(MB_tmp,MBZ_tmp,MSP_tmp) ); occ_tmp=0.0d0

    if ( myrank == 0 ) read(3) occ_tmp

    call MPI_Bcast(occ_tmp,size(occ_tmp),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    
    if ( present(occ_out) ) then
      m1 = b%MB0
      m2 = min( b%MB1, MB_tmp )
      occ_out = 0.0d0
      occ_out(m1:m2,b%MK0:b%MK1,b%MS0:b%MS1) = occ_tmp(m1:m2,b%MK0:b%MK1,b%MS0:b%MS1)
    else
      m1 = min( MB , MB_tmp  )
      m2 = min( MBZ, MBZ_tmp )
      m3 = min( MSP, MSP_tmp )
      occ = 0.0d0
      occ(1:m1,1:m2,1:m3) = occ_tmp(1:m1,1:m2,1:m3)
      if ( myrank == 0 ) write(*,*) "sum(occ)=",sum(occ)
    end if

    deallocate( occ_tmp )

! ---

    if ( flag_newformat ) then

      if ( present(kbb_out) ) then
        if ( myrank == 0 ) read(3) aa0, bb0, kbb_out
        call MPI_Bcast(kbb_out,size(kbb_out),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      else
        if ( myrank == 0 ) read(3)
      end if

      if ( myrank == 0 ) then
        read(3) itmp(14:21)  ! parallel_info(0:7)
        read(3)
        read(3)
        read(3)
        read(3)
        write(*,'(1x,"np(0:7)=",i5,1x,7i4)') itmp(14:21)
      end if
      call MPI_Bcast(itmp(14),8,MPI_Integer,0,MPI_COMM_WORLD,ierr)

! parallelization of the previous calculation
!
      pinfo0%np(0:7) = itmp(14:21) 
      call construct_para( ML_tmp, MB_tmp, MBZ_tmp, MSP_tmp, pinfo0 )
!
! parallelization of the present calculation
!
      call get_np_parallel( pinfo1%np )
      call construct_para( ML, MB, MBZ, MSP, pinfo1 )
!
! -------------------------------------------

    end if !flag_newformat

! ---

    read_wf_is_complex = ( type_wf==0 .or. type_wf==2 )

    present_wf_is_complex = present(z_wf_out)

    if ( myrank == 0 ) then
      if ( SYStype == 0 ) then
        call construct_map_1d_to_3d_grid( LL2 ) 
      else if ( SYStype == 1 ) then
        call construct_map_1d_to_3d_rgrid_mol( LL2 )
      end if
    end if

    if ( .not.(all(pinfo0%np==1).or.all(pinfo1%np==pinfo0%np)) ) then

       call read_wf_general &
       ( pinfo0,pinfo1,ML1,ML2,ML3,MB_tmp,MBZ_tmp,MSP_tmp,MB,MBZ,MSP,type_wf )

    else

! ---

      if ( IO_ctrl0 == 3 ) then
        if ( all(pinfo0%np==1) ) then
          write(cmyrank,'(i5.5)') 0
          IO_ctrl0 = 0
        else
          write(cmyrank,'(i5.5)') myrank
        end if
        file_wf_split = trim(file_wf2)//"."//trim(adjustl(cmyrank))
        if ( myrank == 0 ) close(3)
        open(3,file=file_wf_split,form="unformatted")
      end if

      if ( IO_ctrl0 == 0 ) then
        if ( SYStype == 0 ) then
          if ( read_wf_is_complex ) then
            allocate( utmp3(0:ML1-1,0:ML2-1,0:ML3-1) ); utmp3=(0.0d0,0.0d0)
          else
            allocate( dtmp3(0:ML1-1,0:ML2-1,0:ML3-1) ); dtmp3=0.0d0
          end if
        else if ( SYStype == 1 ) then
          mx=(Ngrid(1)-1)/2
          my=(Ngrid(2)-1)/2
          mz=(Ngrid(3)-1)/2
          if ( read_wf_is_complex ) then
            allocate( utmp3(-mx:mx,-my:my,-mz:mz) ); utmp3=(0.0d0,0.0d0)
          else
            allocate( dtmp3(-mx:mx,-my:my,-mz:mz) ); dtmp3=0.0d0
          end if
        end if
      end if

      nn = ML; if ( IO_ctrl0 == 3 ) nn = n2-n1+1
      if ( read_wf_is_complex ) then
        allocate( utmp(nn) ); utmp=(0.0d0,0.0d0)
      else
        allocate( dtmp(nn) ); dtmp=0.0d0
      end if

      if ( present_wf_is_complex ) then
        if ( .not.allocated(utmp) ) then
          allocate( utmp(nn) ); utmp=(0.0d0,0.0d0)
        end if
      else ![ present_wf is real ]
        if ( .not.allocated(dtmp) ) then
          allocate( dtmp(nn) ); dtmp=0.0d0
        end if
      end if

      select case( type_wf )
      case( 2 ) !single-complex
        allocate( utmpSP(nn) ); utmpSP=(0.0,0.0)
      case( 3 ) !single
        allocate( dtmpSP(nn) ); dtmpSP=0.0
      end select

      do s = 1, MSP_tmp
      do k = 1, MBZ_tmp
      do n = MB1_tmp, MB2_tmp

        if ( IO_ctrl0 == 0 ) then

          call chk_rank_relevance( n,k,s, irank, flag_related )

          if ( myrank == 0 ) then

            select case( type_wf )
            case( 0 )
              read(3) utmp
              do i=1,ML
                i1=LL_tmp(1,i) ; i2=LL_tmp(2,i) ; i3=LL_tmp(3,i)
                utmp3(i1,i2,i3)=utmp(i)
              end do
            case( 1 )
              read(3) dtmp
              do i=1,ML
                i1=LL_tmp(1,i) ; i2=LL_tmp(2,i) ; i3=LL_tmp(3,i)
                dtmp3(i1,i2,i3)=dtmp(i)
              end do
            case( 2 )
              read(3) utmpSP
              do i=1,ML
                i1=LL_tmp(1,i) ; i2=LL_tmp(2,i) ; i3=LL_tmp(3,i)
                utmp3(i1,i2,i3)=utmpSP(i)
              end do
            case( 3 )
              read(3) dtmpSP
              do i=1,ML
                i1=LL_tmp(1,i) ; i2=LL_tmp(2,i) ; i3=LL_tmp(3,i)
                dtmp3(i1,i2,i3)=dtmpSP(i)
              end do
            end select

            if ( read_wf_is_complex ) then
              do i=1,ML
                i1=LL2(1,i) ; i2=LL2(2,i) ; i3=LL2(3,i)
                utmp(i) = utmp3(i1,i2,i3)
              end do
            else
              do i=1,ML
                i1=LL2(1,i) ; i2=LL2(2,i) ; i3=LL2(3,i)
                dtmp(i) = dtmp3(i1,i2,i3)
              end do
            end if

          end if ![ myrank == 0 ]

          call MPI_Barrier( MPI_COMM_WORLD, ierr )

          if ( irank /= 0 .and. myrank_f == 0 ) then

            if ( read_wf_is_complex ) then
              if ( myrank == 0 ) then
                call MPI_Send(utmp,ML,MPI_COMPLEX16,irank,0,MPI_COMM_WORLD,ierr)
              end if
              if ( myrank == irank ) then
                call MPI_Recv(utmp,ML,MPI_COMPLEX16,0,0,MPI_COMM_WORLD,istatus,ierr)
              end if
            else
              if ( myrank == 0 ) then
                call MPI_Send(dtmp,ML,MPI_REAL8,irank,0,MPI_COMM_WORLD,ierr)
              end if
              if ( myrank == irank ) then
                call MPI_Recv(dtmp,ML,MPI_REAL8,0,0,MPI_COMM_WORLD,istatus,ierr)
              end if
            end if

          end if

          call MPI_Barrier( MPI_COMM_WORLD, ierr )

          if ( flag_related ) then
            if ( present_wf_is_complex ) then
              if ( .not.read_wf_is_complex ) utmp=dtmp
              call MPI_Scatterv(utmp,ir_grid,id_grid,MPI_COMPLEX16, &
              z_wf_out(n1,n,k,s),ML0,MPI_COMPLEX16,0,comm_grid,ierr)
            else
              if ( read_wf_is_complex ) dtmp=real(utmp)
              call MPI_Scatterv(dtmp,ir_grid,id_grid,MPI_REAL8, &
              d_wf_out(n1,n,k,s),ML0,MPI_REAL8,0,comm_grid,ierr)
            end if
          end if

        else if ( IO_ctrl0 == 3 ) then

          call get_corresponding_rank( myrank_g, n,k,s, pinfo1, irank )
          call get_corresponding_rank( myrank_g, n,k,s, pinfo0, jrank )

          if ( irank == myrank .or. jrank == myrank ) then

            select case( type_wf )
            case( 0 ) !double-complex

              if ( jrank == myrank ) then
                read(3) utmp(:)
              end if
              if ( irank == jrank ) then
                if ( present_wf_is_complex ) then
                  z_wf_out(:,n,k,s) = utmp(:)
                else
                  d_wf_out(:,n,k,s) = utmp(:)
                end if
              else if ( irank == myrank ) then
                call MPI_Recv(utmp,nn,MPI_COMPLEX16,jrank,0,MPI_COMM_WORLD,istatus,ierr)
                if ( present_wf_is_complex ) then
                  z_wf_out(:,n,k,s) = utmp(:)
                else
                  d_wf_out(:,n,k,s) = utmp(:)
                end if
              else if ( jrank == myrank ) then
                call MPI_Send(utmp,nn,MPI_COMPLEX16,irank,0,MPI_COMM_WORLD,ierr)
              end if

            case( 1 ) !double

              if ( jrank == myrank ) then
                read(3) dtmp(:)
              end if
              if ( irank == jrank ) then ![ irank == jrank == myrank ]
                if ( present_wf_is_complex ) then
                  z_wf_out(:,n,k,s) = dtmp(:)
                else
                  d_wf_out(:,n,k,s) = dtmp(:)
                end if
              else if ( irank == myrank ) then
                call MPI_Recv(dtmp,nn,MPI_REAL8,jrank,0,MPI_COMM_WORLD,istatus,ierr)
                if ( present_wf_is_complex ) then
                  z_wf_out(:,n,k,s) = dtmp(:)
                else
                  d_wf_out(:,n,k,s) = dtmp(:)
                end if
              else if ( jrank == myrank ) then
                call MPI_Send(dtmp,nn,MPI_REAL8,irank,0,MPI_COMM_WORLD,ierr)
              end if

            case( 2 ) !single-complex

              if ( jrank == myrank ) then
                read(3) utmpSP(:)
              end if
              if ( irank == jrank ) then
                if ( present_wf_is_complex ) then
                  z_wf_out(:,n,k,s) = utmpSP(:)
                else
                  d_wf_out(:,n,k,s) = utmpSP(:)
                end if
              else if ( irank == myrank ) then
                call MPI_Recv(utmpSP,nn,MPI_COMPLEX,jrank,0,MPI_COMM_WORLD,istatus,ierr)
                if ( present_wf_is_complex ) then
                  z_wf_out(:,n,k,s) = utmpSP(:)
                else
                  d_wf_out(:,n,k,s) = utmpSP(:)
                end if
              else if ( jrank == myrank ) then
                call MPI_Send(utmpSP,nn,MPI_COMPLEX,irank,0,MPI_COMM_WORLD,ierr)
              end if

            case( 3 ) !single

              if ( jrank == myrank ) then
                read(3) dtmpSP(:)
              end if
              if ( irank == jrank ) then
                if ( present_wf_is_complex ) then
                  z_wf_out(:,n,k,s) = dtmpSP(:)
                else
                  d_wf_out(:,n,k,s) = dtmpSP(:)
                end if
              else if ( irank == myrank ) then
                call MPI_Recv(dtmpSP,nn,MPI_REAL,jrank,0,MPI_COMM_WORLD,istatus,ierr)
                if ( present_wf_is_complex ) then
                  z_wf_out(:,n,k,s) = dtmpSP(:)
                else
                  d_wf_out(:,n,k,s) = dtmpSP(:)
                end if
              else if ( jrank == myrank ) then
                call MPI_Send(dtmpSP,nn,MPI_REAL,irank,0,MPI_COMM_WORLD,ierr)
              end if

            end select ![ type_wf ]

          end if ![ irank, jrank ]

        end if ![ IO_ctrl0 ]

      end do ! n
      end do ! k
      end do ! s

      if ( IO_ctrl0 == 0 .and. myrank == 0 .or. IO_ctrl0 == 3 ) close(3)

      if ( allocated(utmp) ) deallocate( utmp )
      if ( allocated(dtmp) ) deallocate( dtmp )
      if ( allocated(utmpSP) ) deallocate( utmpSP )
      if ( allocated(dtmpSP) ) deallocate( dtmpSP )
      if ( allocated(utmp3) ) deallocate( utmp3 )
      if ( allocated(dtmp3) ) deallocate( dtmp3 )
      if ( allocated(LL_tmp) ) deallocate( LL_tmp )
      if ( allocated(LL2) ) deallocate( LL2 )

    end if ! read_wf_general/read_wf_simple

    call watchb( ttmp, tt, barrier='on' )

    if ( disp_on ) then
      write(*,*) "read from ",file_wf2
      write(*,*) "time(read_wf_simple)=",tt
    end if

    call write_border( 0, " read_wf_simple(end)" )

    return
  end subroutine read_wf_simple

  subroutine chk_rank_relevance( n,k,s, irank, flag_related )
    implicit none
    integer,intent(in) :: n,k,s
    integer,intent(out) :: irank
    logical,intent(out) :: flag_related
    integer :: i1,i2,i3,j1,j2,j3
    do irank = 0, nprocs-1
      i1=id_band(id_class(irank,4))+1
      j1=id_band(id_class(irank,4))+ir_band(id_class(irank,4))
      i2=id_bzsm(id_class(irank,5))+1
      j2=id_bzsm(id_class(irank,5))+ir_bzsm(id_class(irank,5))
      i3=id_spin(id_class(irank,6))+1
      j3=id_spin(id_class(irank,6))+ir_spin(id_class(irank,6))
      if ( id_grid(id_class(irank,0))==0 .and. &
           i1<=n .and. n<=j1 .and. &
           i2<=k .and. k<=j2 .and. &
           i3<=s .and. s<=j3 ) exit
    end do
    if ( irank >= nprocs ) then
      write(*,*) "ERROR!(stpo@simple_wf_io_read): myrank=",myrank
      call stop_program_f( "stop@chk_rank_relevance" )
    end if
    flag_related = .false.
    if ( MK_0_WF <= k .and. k <= MK_1_WF .and. &
         MB_0_WF <= n .and. n <= MB_1_WF  .and. &
         MS_0_WF <= s .and. s <= MS_1_WF ) flag_related = .true.
  end subroutine chk_rank_relevance

  subroutine get_corresponding_rank( mrank_g, n,k,s, p, irank )
    implicit none
    integer,intent(in) :: mrank_g, n,k,s
    type(para),intent(inout) :: p
    integer,intent(out) :: irank
    integer :: i1,i2,i3,i4,i5,i6,i7
    integer :: n0,n1,k0,k1,s0,s1,irank_g
    irank=-1
    loop: do i7=0,p%np(7)-1
    do i6=0,p%np(6)-1
    do i5=0,p%np(5)-1
    do i4=0,p%np(4)-1
      s0=p%spin%id(i6)+1 ; s1=p%spin%id(i6)+p%spin%ir(i6)
      k0=p%bzsm%id(i5)+1 ; k1=p%bzsm%id(i5)+p%bzsm%ir(i5)
      n0=p%band%id(i4)+1 ; n1=p%band%id(i4)+p%band%ir(i4)
      irank_g=-1
      do i3=0,p%np(3)-1
      do i2=0,p%np(2)-1
      do i1=0,p%np(1)-1
        irank=irank+1
        irank_g=irank_g+1
        if ( s0 <= s .and. s <= s1 .and. &
             k0 <= k .and. k <= k1 .and. &
             n0 <= n .and. n <= n1 .and. irank_g == mrank_g ) exit loop
      end do
      end do
      end do
    end do
    end do
    end do
    end do loop
  end subroutine get_corresponding_rank

end module io_read_wf_simple_module
