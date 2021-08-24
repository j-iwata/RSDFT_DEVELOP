module wf_module

  use parallel_module
  use wf_sub_module
  use rsdft_mpi_module
  use rsdft_allgather_module, only: d_rsdft_allgatherv_div
  use io_tools_module
  use watch_module, only: watchb
  use memory_module, only: check_memory

  implicit none

  PRIVATE
  PUBLIC :: unk,esp,esp0,occ,res,init_wf,test_on_wf,gather_wf,gather_b_wf &
           ,ML_WF, ML_0_WF, ML_1_WF, MB_WF, MB_0_WF, MB_1_WF &
           ,MK_WF, MK_0_WF, MK_1_WF, MS_WF, MS_0_WF, MS_1_WF &
           ,Sunk &
           ,hunk, iflag_hunk, workwf &
           ,allocate_work_wf, deallocate_work_wf
  PUBLIC :: write_esp_wf
  PUBLIC :: write_info_esp_wf
  PUBLIC :: wfrange
  public :: allocate_b_dwf
  public :: allocate_b_zwf
  public :: allocate_b_occ
  PUBLIC :: referred_orbital
  PUBLIC :: set_initial_wf
  public :: control_work_wf

  logical, public :: USE_WORKWF_AT_ESPCAL = .true.
  logical, public :: USE_WORKWF_AT_GS     = .true.
  logical, public :: USE_WORKWF_AT_MATE   = .true.
  logical, public :: USE_WORKWF_AT_ROTV   = .true.
  logical, public :: USE_WORKWF_AT_CG     = .true.
  logical :: backup_workwf_sw(5)=.true.

#ifdef _DRSDFT_
  real(8),parameter :: zero=0.d0
  real(8),allocatable :: unk(:,:,:,:)
  real(8),allocatable :: hunk(:,:,:,:)
  real(8),allocatable :: workwf(:,:)
  integer,parameter :: TYPE_MAIN=MPI_REAL8
#else
  complex(8),parameter :: zero=(0.d0,0.d0)
  complex(8),allocatable :: unk(:,:,:,:)
  complex(8),allocatable :: hunk(:,:,:,:)
  complex(8),allocatable :: workwf(:,:)
  integer,parameter :: TYPE_MAIN=RSDFT_MPI_COMPLEX16
#endif

#ifdef _DRSDFT_
    real(8),allocatable :: Sunk(:,:)
#else
    complex(8),allocatable :: Sunk(:,:)
#endif

  real(8),allocatable :: esp(:,:,:), esp0(:,:,:)
  real(8),allocatable :: occ(:,:,:)
  real(8),allocatable :: res(:,:,:)
  logical,allocatable :: referred_orbital(:,:,:)

  integer :: ML_WF, ML_0_WF, ML_1_WF
  integer :: MB_WF, MB_0_WF, MB_1_WF
  integer :: MK_WF, MK_0_WF, MK_1_WF
  integer :: MS_WF, MS_0_WF, MS_1_WF

  integer :: iwork_wf=0
  integer :: iflag_hunk=0
  integer :: nblock_gather_wf=0
  integer :: ndiv_gather_wf=0

  type wfrange
     integer :: ML,ML0,ML1
     integer :: MB,MB0,MB1
     integer :: MK,MK0,MK1,MMK
     integer :: MS,MS0,MS1
  end type wfrange

contains


  subroutine read_wf
    implicit none
    call IOTools_readIntegerKeyword( "WORKWF", iwork_wf )
    call IOTools_readIntegerKeyword( "NBLOCK_GATHER_WF", nblock_gather_wf )
    call IOTools_readIntegerKeyword( "NDIV_GATHER_WF", ndiv_gather_wf )
  end subroutine read_wf


  subroutine init_wf( SYStype_in )
    implicit none
    integer,optional,intent(IN) :: SYStype_in
    integer :: SYStype

    call write_border( 0, " init_wf(start)" )

    call read_wf

    SYStype=0
    if ( present(SYStype_in) ) SYStype=SYStype_in

    ML_WF   = sum( ir_grid )
    ML_0_WF = id_grid(myrank_g) + 1
    ML_1_WF = id_grid(myrank_g) + ir_grid(myrank_g)

    MB_WF   = sum( ir_band )
    MB_0_WF = id_band(myrank_b) + 1
    MB_1_WF = id_band(myrank_b) + ir_band(myrank_b)

    MK_WF   = sum( ir_bzsm )
    MK_0_WF = id_bzsm(myrank_k) + 1
    MK_1_WF = id_bzsm(myrank_k) + ir_bzsm(myrank_k)

    MS_WF   = sum( ir_spin )
    MS_0_WF = id_spin(myrank_s) + 1
    MS_1_WF = id_spin(myrank_s) + ir_spin(myrank_s)

    if ( allocated(occ) ) deallocate(occ)
    if ( allocated(res) ) deallocate(res)
    if ( allocated(esp) ) deallocate(esp)
    if ( allocated(unk) ) deallocate(unk)
    if ( allocated(referred_orbital) ) deallocate(referred_orbital)

    call check_memory( 'wf',ML_1_WF-ML_0_WF+1, MB_WF &
                       ,MK_1_WF-MK_0_WF+1, MS_1_WF-MS_0_WF+1 )

    allocate( unk(ML_0_WF:ML_1_WF,MB_WF,MK_0_WF:MK_1_WF,MS_0_WF:MS_1_WF) )
    unk=zero
    allocate( esp(MB_WF,MK_WF,MS_WF) )
    esp=0.0d0
    allocate( res(MB_WF,MK_WF,MS_WF) )
    res=0.0d0
    allocate( occ(MB_WF,MK_WF,MS_WF) )
    occ=0.0d0
    allocate( referred_orbital(MB_WF,MK_WF,MS_WF) )
    referred_orbital=.true.

    call random_initial_wf
!    call delta_initial_wf
!    call fft_initial_wf_sub( ML_WF,MB_WF,MK_WF,MS_WF,ML_0_WF,ML_1_WF &
!        ,MB_0_WF,MB_1_WF,MK_0_WF,MK_1_WF,MS_0_WF,MS_1_WF,unk )
!    call random_initial_wf_sub( ML_WF,MB_WF,MK_WF,MS_WF,ML_0_WF,ML_1_WF &
!         ,MB_0_WF,MB_1_WF,MK_0_WF,MK_1_WF,MS_0_WF,MS_1_WF,unk,SYStype )

    if ( iwork_wf == 1 ) then
       call allocate_work_wf( iwork_wf )
       USE_WORKWF_AT_ESPCAL = .false.
       USE_WORKWF_AT_GS     = .false.
       USE_WORKWF_AT_MATE   = .true.
       USE_WORKWF_AT_ROTV   = .false.
       USE_WORKWF_AT_CG     = .false.
       backup_workwf_sw = (/ USE_WORKWF_AT_ESPCAL, &
                             USE_WORKWF_AT_GS,     &
                             USE_WORKWF_AT_MATE,   &
                             USE_WORKWF_AT_ROTV,   &
                             USE_WORKWF_AT_CG      /)
    end if

    call write_border( 0, " init_wf(end)" )

  end subroutine init_wf


  subroutine random_initial_wf
    implicit none
    integer :: s,k,n,i
    integer,allocatable :: ir(:)
    real(8) :: u(2)

    call random_seed( size=n )
    allocate( ir(n) )
    ir(:)=MB_0_WF+ML_0_WF
    call random_seed( put=ir )
    deallocate( ir )
 
    do s=MS_0_WF,MS_1_WF
       do k=MK_0_WF,MK_1_WF
          do n=MB_0_WF,MB_1_WF
             do i=ML_0_WF,ML_1_WF
                call random_number(u)
                unk(i,n,k,s)=dcmplx(u(1),u(2))
             end do
          end do
       end do
    end do

  end subroutine random_initial_wf


  subroutine delta_initial_wf
    implicit none
    integer :: s,k,n,i
    unk=(0.0d0,0.0d0)
    i=0
    do s=1,MS_WF
       do k=1,MK_WF
          do n=1,MB_WF
             i=i+1
             if ( ML_0_WF <= i .and. i <= ML_1_WF .and. &
                  MB_0_WF <= n .and. n <= MB_1_WF .and. &
                  MK_0_WF <= k .and. k <= MK_1_WF .and. &
                  MS_0_WF <= s .and. s <= MS_1_WF ) then
                unk(i,n,k,s)=(1.0d0,0.0d0)
             end if
          end do
       end do
    end do
  end subroutine delta_initial_wf


  SUBROUTINE set_initial_wf( indx, wf )
    implicit none
    character(*),intent(IN) :: indx
#ifdef _DRSDFT_
    real(8),intent(OUT)  :: wf(:,:,:,:)
#else
    complex(8),intent(OUT)  :: wf(:,:,:,:)
#endif
    integer :: s,k,n,i
    real(8) :: u(2)
    select case( indx )
    case default
       do s=1,size(wf,4)
       do k=1,size(wf,3)
       do n=1,size(wf,2)
          do i=1,size(wf,1)
             call random_number(u)
             u=u-0.5d0
             wf(i,n,k,s)=dcmplx(u(1),u(2))
          end do ! i
       end do ! n
       end do ! k
       end do ! s
    end select
  END SUBROUTINE set_initial_wf


  SUBROUTINE test_on_wf(dV,disp_switch)
    implicit none
    real(8),intent(IN) :: dV          ! volume element
    logical,intent(IN) :: disp_switch ! diplay switch
    integer :: s,k,m,n,mm
#ifdef _DRSDFT_
    real(8),allocatable :: uu(:,:)
#else
    complex(8),allocatable :: uu(:,:)
#endif

    allocate( uu(MB_WF,MB_WF) ) ; uu=zero

    mm = size(uu)

    do s=MS_0_WF,MS_1_WF
    do k=MK_0_WF,MK_1_WF
       uu(:,:)=zero
       do n=1,MB_WF
       do m=1,n
#ifdef _DRSDFT_
          uu(m,n)=sum( unk(:,m,k,s)*unk(:,n,k,s) )*dV
#else
          uu(m,n)=sum( conjg(unk(:,m,k,s))*unk(:,n,k,s) )*dV
#endif
       end do ! m
       end do ! n
       call rsdft_allreduce_sum( uu, comm_grid )
       do n=1,MB_WF
       do m=1,n
          if ( disp_switch ) then
             write(320,'(1x,i2,i5,2i7,2g25.16)') s,k,m,n,uu(m,n)
          end if
       end do ! m
       end do ! n
    end do ! k
    end do ! s

    deallocate( uu )

  END SUBROUTINE test_on_wf


  subroutine gather_wf
    implicit none
    integer :: k,s
    call write_border( 1, " gather_wf(start)" )
    do s=MS_0_WF,MS_1_WF
    do k=MK_0_WF,MK_1_WF
       call gather_b_wf( k, s )
    end do
    end do
    call write_border( 1, " gather_wf(end)" )
  end subroutine gather_wf

  subroutine gather_b_wf( k, s )
    implicit none
    integer,intent(in) :: k,s
    integer :: mm,nn,mn
    real(8) :: ttmp(2),tttt(2,3)

    if ( MB_0_WF == 1 .and. MB_1_WF == MB_WF ) return
    call write_border( 1, " gather_b_wf(start)" )
    call watchb( ttmp, barrier='on' ); tttt=0.0d0

    mm=size(unk,1)
    ir_band(:)=ir_band(:)*mm
    id_band(:)=id_band(:)*mm

    nn=size(unk,2)
    mn=mm*nn

    if ( ndiv_gather_wf > 0 ) nblock_gather_wf = max( mn/ndiv_gather_wf, 1 )

    if ( nblock_gather_wf > 0 ) then
#ifdef _DRSDFT_
      call d_rsdft_allgatherv_div( mn, unk(:,:,k,s), ir_band, id_band &
                                 , comm_band, nblock_gather_wf )
      if ( iflag_hunk >= 1 ) then
         call d_rsdft_allgatherv_div( mn, hunk(:,:,k,s), ir_band, id_band &
                                    , comm_band, nblock_gather_wf )
      end if
#else
      call rsdft_allgatherv( unk(:,MB_0_WF:MB_1_WF,k,s) &
           , unk(:,:,k,s), ir_band, id_band, comm_band )
      if ( iflag_hunk >= 1 ) then
         call rsdft_allgatherv( hunk(:,MB_0_WF:MB_1_WF,k,s) &
              , hunk(:,:,k,s), ir_band, id_band, comm_band )
      end if
#endif
    else
      call rsdft_allgatherv( unk(:,MB_0_WF:MB_1_WF,k,s) &
           , unk(:,:,k,s), ir_band, id_band, comm_band )
      if ( iflag_hunk >= 1 ) then
         call rsdft_allgatherv( hunk(:,MB_0_WF:MB_1_WF,k,s) &
              , hunk(:,:,k,s), ir_band, id_band, comm_band )
      end if
    end if

    ir_band(:)=ir_band(:)/mm
    id_band(:)=id_band(:)/mm
    call watchb( ttmp, tttt(:,1), barrier='on' )
    !write(*,'(1x,"gather_b_wf"," (",i1,")",2f8.3)') (mm,tttt(:,mm),mm=1,1)
    call write_border( 1, " gather_b_wf(end)" )

  end subroutine gather_b_wf


  SUBROUTINE allocate_work_wf( iflag )
    implicit none
    integer,intent(IN) :: iflag

    call write_border( 0, " allocate_work_wf(start)" )

    iflag_hunk=iflag
    if ( iwork_wf == 0 ) iflag_hunk=0

    if ( myrank == 0 ) then
       write(*,*) "iflag,iwork_wf,iflag_hunk=",iflag,iwork_wf,iflag_hunk
    end if

    if ( iflag == 1 ) then

       if ( allocated(hunk) ) deallocate(hunk)

       allocate( hunk(ML_0_WF:ML_1_WF,MB_WF,MK_0_WF:MK_1_WF,MS_0_WF:MS_1_WF) )

    else if ( iflag == 2 ) then

       if ( allocated(hunk) ) deallocate(hunk)

       allocate( hunk(ML_0_WF:ML_1_WF,MB_WF,MK_WF,MS_0_WF:MS_1_WF) )

    end if

    hunk(:,:,:,:)=zero

    if ( myrank == 0 ) then
       if ( TYPE_MAIN == RSDFT_MPI_COMPLEX16 ) then
          write(*,*) "size(hunk)(MB)=",size(hunk)*16.d0/1024.d0**2
       else if ( TYPE_MAIN == MPI_REAL8 ) then
          write(*,*) "size(hunk)(MB)=",size(hunk)*8.d0/1024.d0**2
       end if
    end if

    call write_border( 0, " allocate_work_wf(end)" )

  END SUBROUTINE allocate_work_wf


  SUBROUTINE deallocate_work_wf
    implicit none
    if ( allocated(hunk) ) deallocate(hunk)
  END SUBROUTINE deallocate_work_wf


  subroutine write_esp_wf( full_info, index_vbm )
    implicit none
    logical,optional,intent(in) :: full_info
    integer,optional,intent(in) :: index_vbm
    integer :: k,n,s,i,n1,n2,nn,fi(2)
    real(8) :: f(6),evbm,ecbm
    character(57) :: header_string, format_string
    call check_disp_length( i, 0 ) ; if ( i < 1 ) return

    if( present(index_vbm) )then
       evbm=maxval( esp(1:index_vbm,:,:) )
       ecbm=minval( esp(index_vbm+1:,:,:) )
       f(1)=ecbm-evbm
       f(2)=(ecbm-evbm)*27.2116d0
       format_string='(/,1x,"index_vbm, Band gap(ht,eV) :",i6,2f15.5)'
       fi(1)=index_vbm
       call write_int_and_real( format_string,1,fi(1:1),2,f ) 
    end if

    write(header_string,'(a4,a6,a20,2a13,1x)') &
         "k","n","esp(n,k,s)","esp_err  ","occ(n,k,s)  "
    call write_string( "" )
    call write_string( header_string )
    nn=sum(occ)
    n1=max( 1, nn/2-5 )
    n2=min( nn/2+5, size(esp,1) )
    if ( present(full_info) ) then
       if ( full_info ) then
          n1=1
          n2=size(esp,1)
       end if
    end if
    if ( i > 1 ) then
       n1=1
       n2=size(esp,1)
    end if
    format_string='(i4,i6,2(f20.15,2g13.5,1x))'
    do k=1,size(esp,2)
    do n=n1,n2
       if ( all(.not.referred_orbital(n,k,:)) ) cycle
       i=0
       do s=1,size(esp,3)
          i=i+1 ; f(i)=esp(n,k,s)
          i=i+1 ; f(i)=esp(n,k,s)-esp0(n,k,s)
          i=i+1 ; f(i)=occ(n,k,s)
       end do
       fi(1:2)=(/ k, n /)
       call write_int_and_real( format_string, 2, fi, i, f ) 
    end do
    end do
  END SUBROUTINE write_esp_wf


  SUBROUTINE write_info_esp_wf( control )
    implicit none
    integer,intent(IN) :: control
    integer :: s,k,n
    integer,parameter :: u=99
    call write_border( 1, " write_info_esp_wf(start)" )
    if ( control > 0 ) then
       if ( control == 2 ) rewind u
       write(u,*) "Eigenvalues"
       write(u,'(a4,a6,a20,2a13,1x)') &
            "k","n","esp(n,k,s)","esp_err  ","occ(n,k,s)  "
       do k=1,MK_WF
       do n=1,MB_WF
          write(u,'(i4,i6,2(f20.15,2g13.5,1x))') k,n &
               ,(esp(n,k,s),esp(n,k,s)-esp0(n,k,s),occ(n,k,s),s=1,MS_WF)
       end do
       end do
    end if
    call write_border( 1, " write_info_esp_wf(end)" )
  END SUBROUTINE write_info_esp_wf


  subroutine allocate_b_zwf( b, wf )
    implicit none
    type(wfrange),intent(inout) :: b
    complex(8),allocatable,intent(inout) :: wf(:,:,:,:)
    allocate( wf(b%ML0:b%ML1,b%MB0:b%MB1,b%MK0:b%MK1,b%MS0:b%MS1) )
    wf=(0.0d0,0.0d0)
  end subroutine allocate_b_zwf

  subroutine allocate_b_dwf( b, wf )
    implicit none
    type(wfrange),intent(inout) :: b
    real(8),allocatable,intent(inout) :: wf(:,:,:,:)
    allocate( wf(b%ML0:b%ML1,b%MB0:b%MB1,b%MK0:b%MK1,b%MS0:b%MS1) )
    wf=0.0d0
  end subroutine allocate_b_dwf

  subroutine allocate_b_occ( b, occup )
    implicit none
    type(wfrange),intent(inout) :: b
    real(8),allocatable,intent(inout) :: occup(:,:,:)
    allocate( occup(b%MB0:b%MB1,b%MK0:b%MK1,b%MS0:b%MS1) )
    occup=0.0d0
  end subroutine allocate_b_occ


  subroutine control_work_wf( ctrl )
    implicit none
    character(*), intent(in) :: ctrl
    select case( ctrl )
    case( "USE_ALL" )
       USE_WORKWF_AT_ESPCAL = .true.
       USE_WORKWF_AT_GS     = .true.
       USE_WORKWF_AT_MATE   = .true.
       USE_WORKWF_AT_ROTV   = .true.
       USE_WORKWF_AT_CG     = .true.
    case( "DEFAULT_SETTING" )
       USE_WORKWF_AT_ESPCAL = backup_workwf_sw(1)
       USE_WORKWF_AT_GS     = backup_workwf_sw(2)
       USE_WORKWF_AT_MATE   = backup_workwf_sw(3)
       USE_WORKWF_AT_ROTV   = backup_workwf_sw(4)
       USE_WORKWF_AT_CG     = backup_workwf_sw(5)
    case default
       write(*,*) ctrl
       call stop_program("keyword is not defined(stop@cnotrol_work_wf)")
    end select
  end subroutine control_work_wf

end module wf_module
