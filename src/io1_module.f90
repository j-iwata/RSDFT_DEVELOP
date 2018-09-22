module io1_module

  use parallel_module, only: para, RSDFT_MPI_COMPLEX16, myrank &
                           , np_grid, ir_grid, id_grid, comm_grid, myrank_g &
                           , id_bzsm, ir_bzsm, id_spin, ir_spin, nprocs &
                           , id_band, ir_band, id_class, myrank_f, comm_fkmb &
                           , construct_para, get_np_parallel
  use wf_module, only: wfrange, MB_WF,MB_0_WF,MB_1_WF,ML_0_WF,ML_1_WF &
                      ,MS_0_WF, MS_1_WF, MS_WF, MK_0_WF, MK_1_WF, MK_WF &
                      , occ, unk, allocate_b_wf, allocate_b_occ &
                      , gather_wf, test_on_wf
  use rgrid_module, only: Ngrid, Igrid, dV
  use rgrid_mol_module, only: LL
  use rsdft_mpi_module !, only: rsdft_allgatherv

  implicit none

  private
  public :: read_data_io1

contains

  subroutine read_data_io1( file_wf2, SYStype, b, wf_out, occ_out, kbb_out )
    implicit none
    character(*),intent(IN) :: file_wf2
    integer,intent(IN) :: SYStype
    type(wfrange),optional,intent(INOUT) :: b
#ifdef _DRSDFT_
    real(8),allocatable,optional,intent(INOUT) :: wf_out(:,:,:,:)
#else
    complex(8),allocatable,optional,intent(INOUT) :: wf_out(:,:,:,:)
#endif
    real(8),allocatable,optional,intent(INOUT) :: occ_out(:,:,:)
    real(8),allocatable,optional,intent(INOUT) :: kbb_out(:,:)
    include 'mpif.h'
    integer :: ierr,irank,jrank,itag,nreq
    integer :: itmp(21),n,k,s,i,i1,i2,i3,j1,j2,j3,mx,my,mz,m1,m2,m3
    integer :: k1,k2,kk
    integer :: ML_tmp,ML1_tmp,ML2_tmp,ML3_tmp,MB_tmp,MB1_tmp,MB2_tmp
    integer :: ML,ML1,ML2,ML3,MB,n1,n2,ML0,MBZ_tmp,MSP_tmp,MMBZ_tmp
    integer :: MB_0,MB_1,MBZ,MBZ_0,MBZ_1,MSP,MSP_0,MSP_1
    integer :: IO_ctrl0, OC, TYPE_WF
    integer,allocatable :: LL2(:,:),LL_tmp(:,:),ir(:),id(:)
    integer,allocatable :: ireq(:),istatus(:,:)
    real(8) :: aa0(3,3),bb0(3,3)
    real(8),allocatable :: kbb0(:,:),occ_tmp(:,:,:)
    logical :: flag_related, flag_newformat, DISP_SWITCH,flag_open
    character(5) :: cmyrank
    character(32) :: file_wf_split
#ifdef _DRSDFT_
    integer,parameter :: TYPE_MAIN=MPI_REAL8
    real(8),parameter :: zero=0.0d0
    real(8),allocatable :: utmp(:), utmp3(:,:,:)
    real(4),allocatable :: utmpSP(:)
#else
    integer,parameter :: TYPE_MAIN=RSDFT_MPI_COMPLEX16
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: utmp(:), utmp3(:,:,:)
    complex(4),allocatable :: utmpSP(:)
#endif
    real(8),allocatable :: dtmp(:)
    real(4),allocatable :: dtmpSP(:)
    type(para) :: pinfo0, pinfo1
    integer :: jrank_g,jrank_n,jrank_k,jrank_s,s1,s2,l1,l2,ln
    integer :: jrank_mod,irank_g,iirank,check_file_id,current_file_id
    real(8) :: sum0,sum1

    call write_border( 0, "read_data_io1(start)" )
    call check_disp_switch( DISP_SWITCH, 0 )

! ---

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    MB  = MB_WF
    m1  = ML_0_WF
    m2  = ML_1_WF
    ML0 = m2 - m1 + 1

    MB_0  = MB_0_WF
    MB_1  = MB_1_WF
    MBZ   = MK_WF
    MBZ_0 = MK_0_WF
    MBZ_1 = MK_1_WF
    MSP   = MS_WF
    MSP_0 = MS_0_WF
    MSP_1 = MS_1_WF

! ---

    allocate( LL2(3,ML)    ) ; LL2=0
    allocate( LL_tmp(3,ML) ) ; LL_tmp=0

    if ( SYStype == 0 ) then

       i=m1-1
       do i3=Igrid(1,3),Igrid(2,3)
       do i2=Igrid(1,2),Igrid(2,2)
       do i1=Igrid(1,1),Igrid(2,1)
          i=i+1
          LL2(1,i)=i1
          LL2(2,i)=i2
          LL2(3,i)=i3
       end do
       end do
       end do

    else if ( SYStype == 1 ) then

       LL2(1:3,m1:m2) = LL(1:3,m1:m2)

    end if

    allocate( ir(0:np_grid-1), id(0:np_grid-1) )

    ir(0:np_grid-1)=3*ir_grid(0:np_grid-1)
    id(0:np_grid-1)=3*id_grid(0:np_grid-1)
    call rsdft_allgatherv( LL2(:,m1:m2), LL2, ir, id, comm_grid )

    deallocate( id,ir )

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

    call mpi_bcast(flag_newformat,1,mpi_logical,0,mpi_comm_world,ierr)
    call mpi_bcast(itmp,size(itmp),mpi_integer,0,mpi_comm_world,ierr)
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

    if ( DISP_SWITCH ) then
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
         ML2_tmp /= ML2 .or. ML3_tmp /= ML3 ) stop "stop@simple_wf_io_read"

    if ( myrank == 0 ) then
       read(3) LL_tmp(1:3,1:ML)
    end if
    call mpi_bcast(LL_tmp,3*ML,mpi_integer,0,mpi_comm_world,ierr)

    if ( DISP_SWITCH ) then
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
          if ( present(wf_out) ) call allocate_b_wf( b, wf_out )
          if ( present(occ_out) ) call allocate_b_occ( b, occ_out )
          if ( present(kbb_out) ) then
             allocate( kbb_out(3,b%MK0:b%MK1) ) ; kbb_out=0.0d0
          end if
       end if

    end if

! ---

    allocate( occ_tmp(MB_tmp,MBZ_tmp,MSP_tmp) ) ; occ_tmp=0.0d0
    if ( myrank == 0 ) read(3) occ_tmp
    call mpi_bcast(occ_tmp,size(occ_tmp),mpi_real8,0,mpi_comm_world,ierr)
    
    if ( present(occ_out) ) then
       n1=b%MB0
       n2=min(b%MB1,MB_tmp)
       occ_out(n1:n2,b%MK0:b%MK1,b%MS0:b%MS1) &
            = occ_tmp(n1:n2,b%MK0:b%MK1,b%MS0:b%MS1)
    else
       n=min(MB ,MB_tmp )
       k=min(MBZ,MBZ_tmp)
       s=min(MSP,MSP_tmp)
       occ(1:n,1:k,1:s)=occ_tmp(1:n,1:k,1:s)
       if ( myrank == 0 ) then
          write(*,*) "sum(occ)=",sum(occ)
       end if
    end if

    deallocate( occ_tmp )

! ---

    if ( flag_newformat ) then

       if ( present(kbb_out) ) then
          if ( myrank == 0 ) read(3) aa0,bb0,kbb_out
          call mpi_bcast(kbb_out,size(kbb_out),mpi_real8,0,mpi_comm_world,ierr)
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
       call mpi_bcast(itmp(14),8,mpi_integer,0,mpi_comm_world,ierr)

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

    end if

! ---

    if ( SYStype == 0 ) then

       allocate( utmp3(0:ML1-1,0:ML2-1,0:ML3-1) ); utmp3=zero

    else if ( SYStype == 1 ) then

       mx=(Ngrid(1)-1)/2
       my=(Ngrid(2)-1)/2
       mz=(Ngrid(3)-1)/2
       allocate( utmp3(-mx:mx,-my:my,-mz:mz) ); utmp3=0.d0

    end if

! ---

    if ( myrank == 0 ) close(3)

!    current_file_id=0
!    if ( myrank < pinfo0%np(0) ) then
!       write(cmyrank,'(i5.5)') myrank
!       file_wf_split = trim(file_wf2)//"."//trim(adjustl(cmyrank))
!       open(3,file=file_wf_split,form="unformatted")
!    end if

! ---

    allocate( utmp(ML) ); utmp=zero
    if ( OC == 4 .or. OC == 5 ) then
       allocate( utmpSP(ML) ); utmpSP=zero
    end if

    if ( type_wf == 1 ) then ! real-wf
       allocate( dtmp(ML) ); dtmp=0.0d0
       if ( OC == 4 .or. OC == 5 ) then
          allocate( dtmpSP(ML) ); dtmpSP=0.0d0
       end if
    end if

! ---

    allocate( istatus(MPI_STATUS_SIZE,nprocs+1) ); istatus=0
    allocate( ireq(nprocs+1) ); ireq=0

    jrank = -1

    do jrank_s=0,pinfo0%spin%np-1

       s1 = pinfo0%spin%id(jrank_s) + 1
       s2 = s1 + pinfo0%spin%ir(jrank_s) - 1

    do jrank_k=0,pinfo0%bzsm%np-1

       k1 = pinfo0%bzsm%id(jrank_k) + 1
       k2 = k1 + pinfo0%bzsm%ir(jrank_k) - 1

    do jrank_n=0,pinfo0%band%np-1

       n1 = pinfo0%band%id(jrank_n) + 1
       n2 = n1 + pinfo0%band%ir(jrank_n) - 1

    do jrank_g=0,pinfo0%grid%np-1

       j1 = pinfo0%grid%id(jrank_g) + 1
       j2 = j1 + pinfo0%grid%ir(jrank_g) - 1

       jrank = jrank + 1

       jrank_mod = mod( jrank, nprocs )

       if ( jrank_mod == myrank ) then
          write(cmyrank,'(i5.5)') jrank
          file_wf_split = trim(file_wf2)//"."//trim(adjustl(cmyrank))
          open(3,file=file_wf_split,form="unformatted")
       end if

       do s=s1,s2
       do k=k1,k2
       do n=n1,n2

          if ( jrank_mod == myrank ) then

             select case( OC )
             case default ! double precision data

                select case( type_wf )
                case( 0 ) ! complex-wf
                   read(3) utmp(j1:j2)
                case( 1 ) ! real-wf
                   read(3) dtmp(j1:j2)
                   utmp(j1:j2)=dtmp(j1:j2)
                case default
                   write(*,*) "type_wf=",type_wf
                   call stop_program_f("Illegal type_wf(stop@read_data_io1)")
                end select

             case( 4, 5 ) ! single precision data

                select case( type_wf )
                case( 0 ) ! complex-wf
                   read(3) utmpSP(j1:j2)
                   utmp(j1:j2)=utmpSP(j1:j2)
                case( 1 ) ! real-wf
                   read(3) dtmpSP(j1:j2)
                   utmp(j1:j2)=dtmpSP(j1:j2)
                case default
                   write(*,*) "type_wf=",type_wf
                   call stop_program_f("Illegal type_wf(stop@read_data_io1_2)")
                end select

             end select ! OC

          end if !( jrank_mod==myrank )

          nreq=0
          itag=0
          do irank_g=0,pinfo1%grid%np-1
             call get_corresponding_rank( irank_g, n,k,s, pinfo1, irank )
             i1 = pinfo1%grid%id(irank_g) + 1
             i2 = i1 + pinfo1%grid%ir(irank_g) - 1
             if ( j1<=i1.and.i1<=j2 .or. j1<=i2.and.i2<=j2 ) then
                l1 = max( j1, i1 )
                l2 = min( j2, i2 )
                ln = l2 - l1 + 1
                if ( jrank_mod == myrank ) then
                   nreq=nreq+1
                   call MPI_Isend( utmp(l1),ln,TYPE_MAIN &
                        ,irank,itag,MPI_COMM_WORLD,ireq(nreq),ierr )
                end if
                if ( irank == myrank ) then
                   nreq=nreq+1
                   call MPI_Irecv( unk(l1,n,k,s),ln,TYPE_MAIN,jrank_mod &
                        ,itag,MPI_COMM_WORLD,ireq(nreq),ierr )
                end if
             end if
          end do ! irank_g

          if ( nreq > 0 ) call MPI_Waitall( nreq, ireq, istatus, ierr )

       end do ! n
       end do ! k
       end do ! s

       if ( jrank_mod == myrank ) close(3)

    end do ! jrank_g
    end do ! jrank_n
    end do ! jrank_k
    end do ! jrank_s

    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=MB_0 ,MB_1
       call rsdft_allgatherv( unk(m1:m2,n,k,s), utmp &
            ,pinfo1%grid%ir, pinfo1%grid%id, comm_grid )
       do i=1,ML
          utmp3(LL_tmp(1,i),LL_tmp(2,i),LL_tmp(3,i))=utmp(i)
       end do
       do i=m1,m2
          unk(i,n,k,s)=utmp3(LL2(1,i),LL2(2,i),LL2(3,i))
       end do
    end do
    end do
    end do

! ---

    deallocate( ireq )
    deallocate( istatus )
    if ( allocated(dtmpSP) ) deallocate( dtmpSP )
    if ( allocated(dtmp) ) deallocate( dtmp )
    if ( allocated(utmpSP) ) deallocate( utmpSP )
    deallocate( utmp )
    deallocate( utmp3 )
    deallocate( LL_tmp )
    deallocate( LL2 )

#ifdef _DRSDFT_
    call mpi_bcast( unk,size(unk),MPI_REAL8,0,comm_fkmb,ierr )
#else
    call mpi_bcast( unk,size(unk),RSDFT_MPI_COMPLEX16,0,comm_fkmb,ierr )
#endif

    call write_border( 0, "read_data_io1(end)" )

    return
  end subroutine read_data_io1

  subroutine get_corresponding_rank( mrank_g, n,k,s, p, irank )
    implicit none
    integer,intent(IN) :: mrank_g, n,k,s
    type(para),intent(INOUT) :: p
    integer,intent(OUT) :: irank
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

end module io1_module
