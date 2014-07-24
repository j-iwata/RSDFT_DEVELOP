MODULE PSQRijPrep
  use parallel_module, only: myrank,myrank_g,np_grid,node_partition,COMM_GRID
  use aa_module, only: aa
  use atom_module, only: Natom,ki_atom,aa_atom
  use rgrid_module, only: Igrid,Ngrid,Hgrid
  use ps_nloc2_variables
  use VarPSMember
  ! lo,Nelement_,norb,viod
  use VarPSMemberG
  use VarParaPSnonLocG
  use ps_nloc2_module, only: prepMapsTmp
  use pseudopot_module, only: pselect
  use minimal_box_module
  implicit none

!  include 'mpif.h'
  PRIVATE
  PUBLIC :: prepQRijp102
  
  complex(8),allocatable :: qaL(:,:) !qaL(k3max,Lrefmax)

  complex(8),parameter :: z0=(0.d0,0.d0),z1=(1.d0,0.d0),zi=(0.d0,1.d0)
  real(8),parameter :: pi4=16.d0*atan(1.d0)
  real(8),allocatable :: y2a(:,:,:,:)

CONTAINS

  SUBROUTINE allocateQaL
    implicit none
    if (allocated(qaL)) deallocate(qaL)
    allocate( qaL(k2max,max_Lref) ) ; qaL=z0
    return
  END SUBROUTINE allocateQaL


#ifdef _FOR_COMPILE_
  SUBROUTINE prepQRijp102
    implicit none
    return
  END SUBROUTINE prepQRijp102
#else
  SUBROUTINE prepQRijp102
    implicit none
    real(8) :: x,y,z
    integer :: l,m

    integer :: ia,ik,ik1,ik2,ik3,il,ir,ir0
    integer :: i,j,i1,i2,i3
    integer :: iorb1,iorb2
    integer :: NRc,n
    real(8) :: Rps2

    integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab2,ab3
    integer :: ML1,ML2,ML3
    integer :: np1,np2,np3

    real(8) :: r
    integer :: mm1,mm2,mm3
    integer :: MMJJ_0,nrqr

    integer :: nl3vmax
    integer,allocatable :: icheck_tmp5(:,:)
    integer,allocatable :: JJ_tmp(:,:,:,:,:),MJJ_tmp_Q(:,:)
    real(8),allocatable :: uV_tmp(:,:,:)

    integer :: M_irad
    integer,allocatable :: irad(:,:)

    real(8) :: c1,c2,c3
    real(8) :: maxerr
    integer :: c_nzqr_pre

    real(8) :: Rx,Ry,Rz
    integer :: ic1,ic2,ic3,id1,id2,id3
    integer :: k1,k2,k3,i1_0,i2_0,i3_0
    integer :: d1,d2,d3
    real(8) :: r2

    integer :: ll3,mm,m1,m2,iqr
    real(8) :: v,v0,err,err0,ep,QRtmp

    logical,allocatable :: lcheck_tmp1(:,:)

    real(8),allocatable :: QRij_tmp(:,:,:,:)
    integer :: c_nzqr,c_nzqr_0
    integer,allocatable :: itmp(:,:),icheck_tmp1(:),icheck_tmp2(:),icheck_tmp4(:,:,:)
    integer,allocatable :: qr_nsend_tmp(:),nl_rank_map_tmp(:),maps_tmp(:)
    integer,allocatable :: sendmap_tmp(:,:),recvmap_tmp(:,:)
    
    integer :: ierr
    real(8) :: ctt(0:9),ett(0:9)
    integer,allocatable :: ireq(:)
    ,allocatable :: istatus(:,:)

write(400+myrank,*) ">>>>> prepQRijp102"

    INTERFACE
       FUNCTION Ylm(x,y,z,l,m)
         real(8) :: Ylm
         real(8),intent(IN) :: x,y,z
         integer,intent(IN) :: l,m
       END FUNCTION Ylm
    END INTERFACE

    ctt=0.0d0 ; ett=0.0d0
    call watch(ctt(6),ett(6))

  call allocateQaL

    if ( .not.allocated(y2a) ) then
      allocate( y2a(max_psgrd,max_Lref,max_k2,Nelement_) )
      y2a=0.d0
      do ik=1,Nelement_
        do ik2=1,max_k2
          iorb=k2_to_iorb(2,ik2,ik)
          do il=1,max_Lref
            d1=0.d0
            d2=0.d0
            call spline(rad1(1,ik),qrL(1,il,ik2,ik),NRps(iorb,ik),d1,d2,y2a(1,il,ik2,ik))
          end do
        end do
      end do
    end if

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = Igrid(2,1)-Igrid(1,1)+1
    ab2 = Igrid(2,2)-Igrid(1,2)+1
    ab3 = Igrid(2,3)-Igrid(1,3)+1

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    call watch(ctt(7),ett(7))

    r=maxval(Rps)+maxval(Hgrid(1:3))+1.d-8
    call make_minimal_box(r,mm1,mm2,mm3,MMJJ_0)
    
    call watch(ctt(8),ett(8))

    if ( .not.allocated(icheck_tmp5) ) then
      nl3vmax=maxval(nl3v)
      allocate( icheck_tmp5(k1max,Natom)             ) ; icheck_tmp5=0
      allocate( JJ_tmp(6,MMJJ_0,k1max,Natom)         ) ; JJ_tmp=0
      allocate( MJJ_tmp_Q(k1max,Natom)               ) ; MJJ_tmp_Q=0
      allocate( QRij_tmp(MMJJ_0,nl3vmax,k1max,Natom) ) ; QRij_tmp=0.d0
    end if

    call watch(ctt(0),ett(0))
    
#ifndef _SPLINE_
    allocate( irad(0:3000,Nelement_) ) ; irad=0
    M_irad=0
    do ik=1,Nelement_
       NRc=maxval( NRps(:,ik) )
       NRc=min( 3000, NRc )
       m=0
       irad(0,ik)=1
       do ir=1,NRc
          m=int(100.d0*rad1(ir,ik))+1
          irad( m,ik )=ir
       end do
       ir=irad(0,ik)
       do i=1,m
          if ( irad(i,ik)==0 ) then
             irad(i,ik)=ir
             cycle
          end if
          ir=irad(i,ik)
       end do
       irad(m+1:,ik)=ir
       M_irad=max(M_irad,m)
    end do
#endif

    c1          = 1.d0/Ngrid(1)
    c2          = 1.d0/Ngrid(2)
    c3          = 1.d0/Ngrid(3)
    maxerr      = 0.d0
    MMJJ_Q      = 0
    MMJJ_t_Q    = 0
    c_nzqr_pre  = 0

!!$OMP parallel do schedule(dynamic) firstprivate( maxerr ) &
!!$OMP    private( Rx,Ry,Rz,ic1,ic2,ic3,ik,iorb,Rps2,NRc,L,j,i,i1,i2,i3 &
!!$OMP            ,id1,id2,id3,k1,k2,k3,i1_0,i2_0,i3_0,d1,d2,d3,x,y,z,r2,r &
!!$OMP            ,v0,err0,ir0,ir,mm,m1,m2,v,err )
    j=0
    do ia=1,Natom
      Rx = aa(1,1)*aa_atom(1,ia)+aa(1,2)*aa_atom(2,ia)+aa(1,3)*aa_atom(3,ia)
      Ry = aa(2,1)*aa_atom(1,ia)+aa(2,2)*aa_atom(2,ia)+aa(2,3)*aa_atom(3,ia)
      Rz = aa(3,1)*aa_atom(1,ia)+aa(3,2)*aa_atom(2,ia)+aa(3,3)*aa_atom(3,ia)
      ic1 = nint( aa_atom(1,ia)*Ngrid(1) )
      ic2 = nint( aa_atom(2,ia)*Ngrid(2) )
      ic3 = nint( aa_atom(3,ia)*Ngrid(3) )
      ik = ki_atom(ia)
      do ik1=1,N_k1(ik)
        ik2=k1_to_k2(ik1,ik)
        ik3=k1_to_k3(ik1,ik)
        iorb1=k1_to_iorb(1,ik1,ik)
        iorb2=k1_to_iorb(2,ik1,ik)
        Rps2=max(Rps(iorb1,ik)**2,Rps(iorb2,ik)**2)
        NRc =max(NRps(iorb1,ik),NRps(iorb2,ik))
          do i=1,MMJJ_0
            i1 = map_grid_ion(1,i)
            i2 = map_grid_ion(2,i)
            i3 = map_grid_ion(3,i)
            id1 = ic1 + i1
            id2 = ic2 + i2
            id3 = ic3 + i3
            k1=id1/ML1 ; if ( id1<0 ) k1=(id1+1)/ML1-1
            k2=id2/ML2 ; if ( id2<0 ) k2=(id2+1)/ML2-1
            k3=id3/ML3 ; if ( id3<0 ) k3=(id3+1)/ML3-1
            i1_0=id1-k1*ML1
            i2_0=id2-k2*ML2
            i3_0=id3-k3*ML3
            if ( Igrid(1,1) <= i1_0 .and. i1_0 <= Igrid(2,1) .and. &
                Igrid(1,2) <= i2_0 .and. i2_0 <= Igrid(2,2) .and. &
                Igrid(1,3) <= i3_0 .and. i3_0 <= Igrid(2,3) ) then
              j=j+1
              d1 = id1*c1
              d2 = id2*c2
              d3 = id3*c3
              x  = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
              y  = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
              z  = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
              r2 = x*x+y*y+z*z
              if ( r2 > Rps2+1.d-10 ) cycle

              r    = sqrt(r2)

              call get_qaL(r,x,y,z)
              
              do ll3=1,nl3v(ik2,ik)
                l=l3v(ll3,ik2,ik)-1
                v0   = 0.d0
                err0 = 0.d0
                QRtmp= 0.d0

                if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 .or. &
                     abs(z)>1.d-14 .or. L==0 ) then
#ifdef _SPLINE_
                  call splint(rad1(1,ik),qrL(1,l,ik2,ik),y2a(1,l,ik2,ik),NRc,r,v0)
#else
                  ir0=irad( int(100.d0*r),ik )
                  do ir=ir0,NRc
                    if ( r<rad1(ir,ik) ) exit
                  end do
                  if ( ir <= 2 ) then
                    v0=qrL(2,ll3,ik2,ik)
                    if ( ir < 1 ) stop "prep_QRij_p102(0)"
                  else if ( ir <= NRc ) then
                    err0=1.d10
                    do mm=1,20
                      m1=max(1,ir-mm)
                      m2=min(ir+mm,NRc)
                      call polint &
                      (rad1(m1,ik),qrL(m1,ll3,ik2,ik),m2-m1+1,r,v,err)
                      if ( abs(err)<err0 ) then
                        v0=v
                        err0=abs(err)
                        if ( err0<ep ) exit
                      end if
                    end do
                  else
                    write(*,*) "prep_QRij_p102(1)",ir,NRc,qrL(NRc,ll3,ik2,ik)
                    write(*,*) qrL(NRc+1,ll3,ik3,ik),r,rad1(ir,ik)
                    stop
                  end if
                  maxerr=max(maxerr,err0)
#endif
                  QRtmp=v0/pi4* real(qaL(ik3,ll3)/(-zi)**L)
                end if ! x,y,z

                QRij_tmp(j,ik1,ia) = QRij_tmp(j,ik1,ia)+QRtmp
              end do ! ll3

              JJ_tmp(1,j,ik1,ia) = i1_0
              JJ_tmp(2,j,ik1,ia) = i2_0
              JJ_tmp(3,j,ik1,ia) = i3_0
              JJ_tmp(4,j,ik1,ia) = k1
              JJ_tmp(5,j,ik1,ia) = k2
              JJ_tmp(6,j,ik1,ia) = k3
            end if ! Igrid
          end do ! i ( 1 - MMJJ_0 )
          MJJ_tmp_Q(ik1,ia)=j
          c_nzqr_pre=c_nzqr_pre+1
          icheck_tmp5(ik1,ia)=c_nzqr_pre
       end do ! ik1
    end do ! ia
!!$OMP end parallel do
    MMJJ_Q = maxval( MJJ_tmp_Q )

#ifndef _SPLINE_
    deallocate( irad )
#endif
    
    Mqr=0
    do ia=1,Natom
      ik=ki_atom(ia)
      do ik1=1,N_k1(ia)
        Mqr=Mqr+1
      end do
    end do

    allocate( lcheck_tmp1(Mqr,0:np_grid-1) )
    lcheck_tmp1(:,:)=.false.
    iqr=0
    do ia=1,Natom
      ik=ki_atom(ia)
      do ik1=1,N_k1(ia)
        iqr=iqr+1
        if (j>0) then
          lcheck_tmp1(iqr,myrank_g)=.true.
        end if
      end do
    end do
    call mpi_allgather(lcheck_tmp1(1,myrank_g),Mqr,mpi_logical &
                      ,lcheck_tmp1,Mqr,mpi_logical,comm_grid,ierr)
    
    call watch(ctt(1),ett(1))

! for grid-parallel computation

    c_nzqr_0 = N_nzqr

    n=maxval( node_partition(1:3) )
    allocate( itmp(n,3)                   ) ; itmp=0
    allocate( qr_nsend_tmp(0:nprocs_g-1)  ) ; qr_nsend_tmp=0
    allocate( icheck_tmp1(0:nprocs_g-1)   ) ; icheck_tmp1=0
    allocate( icheck_tmp2(0:nprocs_g-1)   ) ; icheck_tmp2=0
    allocate( nl_rank_map_tmp(0:nprocs_g) ) ; nl_rank_map_tmp=-1
    allocate( maps_tmp(c_nzqr_0)          ) ; maps_tmp=0
    allocate( sendmap_tmp(c_nzqr_0,0:nprocs_g-1) ) ; sendmap_tmp=0
    allocate( recvmap_tmp(c_nzqr_0,0:nprocs_g-1) ) ; recvmap_tmp=0

    np1 = node_partition(1)
    np2 = node_partition(2)
    np3 = node_partition(3)

    nrqr=0
    do ia=1,Natom
      ik=ki_atom(ia)
      do ik1=1,N_k1(ik)
        icheck_tmp1(:)=0
        do n=0,np_grid-1
          if ( lcheck_tmp1(lma,n) ) icheck_tmp1(n) = 1
        end do
        icheck_tmp1(myrank_g) = icheck_tmp5(ik1,ia)

        call prepMapsTmp(np1,np2,np3,nprocs_g,itmp,icheck_tmp1,icheck_tmp2)

        if ( icheck_tmp1(myrank_g)/=0 ) then
          if ( icheck_tmp1(myrank_g)>0 ) then
            maps_tmp(icheck_tmp2(myrank_g),1)=icheck_tmp1(myrank_g)
          end if
          maps_tmp(icheck_tmp2(myrank_g),2)=ia
          maps_tmp(icheck_tmp2(myrank_g),3)=ik1
          do n=0,nprocs_g-1
            if ( n==myrank_g .or. icheck_tmp1(n)==0 ) cycle
            qr_nsend_tmp(n)=qr_nsend_tmp(n)+1  
            sendmap_tmp(qr_nsend_tmp(n),n)=icheck_tmp2(myrank_g)
            recvmap_tmp(qr_nsend_tmp(n),n)=icheck_tmp2(n)
            if ( any(nl_rank_map_tmp(0:nrqr)==n) ) cycle
            nrqr=nrqr+1
            nl_rank_map_tmp(nrqr)=n
          end do
        end if
      end do ! ik
    end do ! ia
    
    call watch(ctt(2),ett(2))

    c_nzqr = icheck_tmp2(myrank_g)

    deallocate( itmp )
    deallocate( icheck_tmp1 )
    deallocate( icheck_tmp2 )
    deallocate( icheck_tmp5 )
    deallocate( lcheck_tmp1 )
!===================================================================
    if ( allocated(QRij) ) then
       deallocate( uV )
       deallocate( JJ_MAP )
       deallocate( MJJ_MAP )
       deallocate( MJJ )
       deallocate( iuV )
       deallocate( amap )
       deallocate( lmap )
       deallocate( mmap )
       deallocate( iorbmap )
    end if
    allocate( uV(MMJJ,nzlma)       ) ; uV=0.d0
    allocate( JJ_MAP(6,MMJJ,nzlma) ) ; JJ_MAP=0
    allocate( MJJ_MAP(nzlma)       ) ; MJJ_MAP=0
    allocate( MJJ(nzlma)           ) ; MJJ=0
    allocate( iuV(nzlma)           ) ; iuV=0
    allocate( amap(nzlma)          ) ; amap=0
    allocate( lmap(nzlma)          ) ; lmap=0
    allocate( mmap(nzlma)          ) ; mmap=0
    allocate( iorbmap(nzlma)       ) ; iorbmap=0
!===================================================================
    
    if (allocated(nl_rank_map_Q)) deallocate(nl_rank_map_Q)
    allocate( nl_rank_map_Q(nrqr) )
    do i=1,nrqr
       nl_rank_map_Q(i)=nl_rank_map_tmp(i)
    end do
    deallocate( nl_rank_map_tmp )

!!$OMP parallel do private( a,l,m,iorb,Rx,Ry,Rz,j,i1,i2,i3,k1,k2,k3,d1,d2,d3,x,y,z )
    do iqr=1,c_nzqr
       if ( maps_tmp(lma,1) == 0 ) cycle
       ia    = maps_tmp(iqr,2)
       ik1   = maps_tmp(iqr,3)
       MJJ_MAP_Q(lma) = MJJ_tmp_Q(ik1,ia)
       do j=1,MJJ_MAP(lma)
          QRij(j,lma) = QRij_tmp(j,iorb,a)*Ylm(x,y,z,l,m)
          JJ_MAP(1:6,j,lma) = JJ_tmp(1:6,j,iorb,a)
       end do
    end do
!!$OMP end parallel do

!!$OMP parallel do private( a,l,m,iorb,Rx,Ry,Rz,j,i1,i2,i3,k1,k2,k3,d1,d2,d3,x,y,z )
    do kk1=1,c_nzqr
       l=maps_tmp(kk1)
       if (l==0) cycle
!===================================================================
       ik1=
       ia=
!===================================================================
       MJJ_MAP_Q(kk1) = MJJ_tmp_Q(ik1,ia)
       do j=1,MJJ_MAP_Q(kk1)
!===================================================================
          QRij_tmp(j,kk1) = QRij_tmp(j,iorb,a)
!===================================================================
          JJ_MAP_Q(1:6,j,kk1) = JJ_tmp(1:6,j,ik1,ia)
       end do
    end do
!!$OMP end parallel do

    deallocate( MJJ_tmp_Q )
    deallocate( JJ_tmp    )
    deallocate( QRij_tmp  )
    deallocate( maps_tmp  )

    call watch(ctt(3),ett(3))

    allocate( icheck_tmp4(a1b:b1b,a2b:b2b,a3b:b3b) )
    icheck_tmp4=0
    do kk1=1,c_nzqr
       j=0
       icheck_tmp4=0
       do i=1,MJJ_MAP_Q(kk1)
          i1=JJ_MAP(1,i,kk1)
          i2=JJ_MAP(2,i,kk1)
          i3=JJ_MAP(3,i,kk1)
          if ( icheck_tmp4(i1,i2,i3)==0 ) then
             j=j+1
             icheck_tmp4(i1,i2,i3)=j
          end if
       end do
       MJJ(lma)=j
    end do
    MAXMJJ = maxval( MJJ(1:c_nzqr) )
    deallocate( icheck_tmp4 )

    nl_max_send = maxval( qr_nsend_tmp )

!===================================================================
    if ( allocated(qr_nsend) ) then
       deallocate( qr_nsend )
       deallocate( sendmap_Q )
       deallocate( recvmap_Q )
    end if
    allocate( qr_nsend(0:nprocs_g-1) ) ; qr_nsend=0
    allocate( sendmap_Q(nl_max_send,0:nprocs_g-1) ) ; sendmap_Q=0
    allocate( recvmap_Q(nl_max_send,0:nprocs_g-1) ) ; recvmap_Q=0
!===================================================================

    do n=0,nprocs_g-1
       sendmap_Q(1:nl_max_send,n) = sendmap_tmp(1:nl_max_send,n)
       qr_nsend(n) = qr_nsend_tmp(n)
    end do

!===================================================================
    allocate( ireq(2*nprocs_g) ) ; ireq=0
    allocate( istatus(MPI_STATUS_SIZE,2*nprocs_g) ) ; istatus=
!===================================================================
    nreq=0
    do n=0,nprocs_g-1
       if ( qr_nsend(n)<=0 .or. n==myrank_g ) cycle
       nreq=nreq+1
       call mpi_isend(recvmap_tmp(1,n),qr_nsend(n),mpi_integer,n,1 &
            ,comm_grid,ireq(nreq),ierr)
       nreq=nreq+1
       call mpi_irecv(recvmap_Q(1,n) ,qr_nsend(n),mpi_integer,n,1 &
            ,comm_grid,ireq(nreq),ierr)
    end do
    call mpi_waitall(nreq,ireq,istatus,ierr)
    deallocate( istatus )
    deallocate( ireq )

    deallocate( recvmap_tmp,sendmap_tmp,qr_nsend_tmp )

    call prepThreeWayComm( nrqr,nl_rank_map_Q,nrqr_xyz,num_2_rank_Q )

    call watch(ctt(4),ett(4))

!===================================================================
    call allocate_ps_nloc2(MB_d)
!===================================================================

!===================================================================
    if ( allocated(JJP) ) then
       deallocate( JJP )
       deallocate( uVk )
    end if
    allocate( JJP(MAXMJJ,nzlma) ) ; JJP=0
    allocate( uVk(MAXMJJ,nzlma,MBZ_0:MBZ_1) ) ; uVk=0.d0
!===================================================================

    call prep_uvk_ps_nloc2(MBZ_0,MBZ_1,kbb(1,MBZ_0))

    call watch(ctt(5),ett(5))

    if ( disp_switch_parallel ) then
       write(*,*) "time(prepQRijp102_1)",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(prepQRijp102_2)",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(prepQRijp102_3)",ctt(3)-ctt(2),ett(3)-ett(2)
       write(*,*) "time(prepQRijp102_4)",ctt(4)-ctt(3),ett(4)-ett(3)
       write(*,*) "time(prepQRijp102_5)",ctt(5)-ctt(4),ett(5)-ett(4)
       write(*,*) "time(prepQRijp102_7)",ctt(7)-ctt(6),ett(7)-ett(6)
       write(*,*) "time(prepQRijp102_8)",ctt(8)-ctt(7),ett(8)-ett(7)
       write(*,*) "time(prepQRijp102_9)",ctt(0)-ctt(8),ett(0)-ett(8)
    end if
    
write(400+myrank,*) "<<<<< prepQRijp102"
    return
  END SUBROUTINE prepQRijp102
#endif

!---------------------------------------------------------------------------------------

  SUBROUTINE get_qaL ( r,x,y,z )
    implicit none

    real(8) :: r,x,y,z
    real(8),parameter :: sq3=sqrt(3.d0)
    real(8),parameter :: sq5=sqrt(5.d0)

!----------------------------
    qaL(:,:)   = z0

    qaL(1,1)   = z1
    qaL(3,1)   = z1
    qaL(6,1)   = z1
    qaL(10,1)  = z1
    qaL(15,1)  = z1
    qaL(21,1)  = z1
    qaL(28,1)  = z1
    qaL(36,1)  = z1
    qaL(45,1)  = z1

!--      
    if ( r <= 1.0D-10 ) return
!--      
    x = x/r 
    y = y/r 
    z = z/r 

    qaL(1,1)=1.d0
    qaL(2,1)=sq3*zi*y
    qaL(3,1)=1.d0
    qaL(3,2)=(3.d0*x**2-3.d0*y**2+3.d0*z**2-1.d0)/2.d0
    qaL(4,1)=-sq3*zi*z
    qaL(5,1)=0.d0
    qaL(5,2)=3.d0*y*z
    qaL(6,1)=1.d0
    qaL(6,2)=-(3.d0*z**2-1.d0)
    qaL(7,1)=sq3*zi*x
    qaL(8,1)=0.d0
    qaL(8,2)=-3.d0*x*y
    qaL(9,1)=0.d0
    qaL(9,2)=3.d0*x*z
    qaL(10,1)=1.d0
    qaL(10,2)=-(3.d0*x**2-3.d0*y**2-3.d0*z**2+1.d0)/2.d0
    qaL(11,1)=-sqrt(15.d0)*x*y
    qaL(12,1)=3.d0*sq5*zi*x/5.d0
    qaL(12,2)=(   3.d0*sq5*zi*x * (5.d0*x**2-15.d0*y**2+5.d0*z**2-1.d0) )/20.d0
    qaL(13,1)=0.d0
    qaL(13,2)=3.d0*sq5*zi*x*y*z
    qaL(14,1)=3.d0*sq5*zi*y/5.d0
    qaL(14,2)=-(3.d0*sq5*zi*y*(15.d0*x**2 -5.d0*y**2-5.d0*z**2+1.d0))/(20.d0)
    qaL(15,1)=1.d0
    qaL(15,2)=(5.d0*(3.d0*z**2-1.d0))/7.d0
    qaL(15,3)=-(3.d0*(35.d0*x**4-210.d0*x**2*y**2 +35.d0*y**4-35.d0*z**4+30.d0*z**2-3.d0))/56.d0
    qaL(16,1)=3.d0*sqrt(5.d0/3.0d0)*y*z
    qaL(17,1)=-(3.d0*sq5*zi*z)/5.d0
    qaL(17,2)=-(3.d0*sq5*zi*z *(5.d0*x**2-5.d0*y**2+5.d0*z**2-3.d0))/10.d0
    qaL(18,1)=(3.d0*sq5*zi*y)/5.d0
    qaL(18,2)=-(3.d0*zi*y*(5.d0*z**2-1.d0))/sq5
    qaL(19,1)=0.d0
    qaL(19,2)=3.d0*sq5*zi*x*y*z
    qaL(20,1)=0.d0
    qaL(20,2)=(15.d0*x*z)/7.d0
    qaL(20,3)=(15.d0*x*z *(7.d0*x**2-21.d0*y**2+7.d0*z**2-3.d0))/28.d0
    qaL(21,1)=1.d0
    qaL(21,2)=(5.d0*(3.d0*x**2-3.d0*y**2-3.d0*z**2+1.d0))/14.d0
    qaL(21,3)=-(3.d0*(35.d0*x**2*z**2-5.d0*x**2-35.d0*y**2*z**2 +5.d0*y**2+35.d0*z**4-30.d0*z**2+3.d0))/14.d0
    qaL(22,1)=-(sq5*(3.d0*z**2-1.d0))/2.d0
    qaL(23,1)=-(sqrt(15.d0)*zi*y)/5.d0
    qaL(23,2)=-1.5d0*sqrt(3.d0/5.d0)*zi*y*(5.d0*z**2-1.d0)
    qaL(24,1)=-(2.d0*sqrt(15.d0)*zi*z)/5.d0
    qaL(24,2)=9.d0*zi*z*(5.d0*z**2-3.d0)/(2.d0*sqrt(15.d0))
    qaL(25,1)=-(sqrt(15.d0)*zi*x)/5.d0
    qaL(25,2)=-1.5d0*sqrt(3.d0/5.d0)*zi*x*(5.d0*z**2-1.d0)
    qaL(26,1)=0.d0
    qaL(26,2)=10.d0*sq3*x*y/7.d0
    qaL(26,3)=15.d0*sq3*x*y*(7.d0*z**2-1.d0)/14.d0
    qaL(27,1)=0.d0
    qaL(27,2)=(15.d0*y*z)/(7.d0*sq3)
    qaL(27,3)=-15.d0*sq3*y*z*(7.d0*z**2-3.d0)/14.d0
    qaL(28,1)=1.d0
    qaL(28,2)=-(5.d0*(3.d0*z**2-1.d0))/7.d0
    qaL(28,3)=(9.d0*(35.d0*z**4-30.d0*z**2+3.d0))/28.d0
    qaL(29,1)=(3.d0*sqrt(5.d0/3.d0)*x*z)
    qaL(30,1)=0.d0
    qaL(30,2)=3.d0*sq5*zi*x*y*z
    qaL(31,1)=(3.d0*sq5*zi*x)/5.d0
    qaL(31,2)=-(3.d0*zi*x*(5.d0*z**2-1.d0))/sq5
    qaL(32,1)=-(3.d0*sq5*zi*z)/5.d0
    qaL(32,2)=0.3d0*sq5*zi*z*(5.d0*x**2-5.d0*y**2-5.d0*z**2+3.d0)
    qaL(33,1)=0.d0
    qaL(33,2)=(15.d0*y*z)/7.d0
    qaL(33,3)=-15.d0*y*z *(21.d0*x**2-7.d0*y**2-7.d0*z**2+3.d0)/28.d0
    qaL(34,1)=0.d0
    qaL(34,2)=-(15.d0*x*y)/7.d0
    qaL(34,3)=(15.d0*x*y*(7.d0*z**2-1.d0))/7.d0
    qaL(35,1)=0.d0
    qaL(35,2)=(15.d0*x*z)/(7.d0*sq3)
    qaL(35,3)=-15.d0*sq3*x*z*(7.d0*z**2-3.d0)/14.d0
    qaL(36,1)=1.d0
    qaL(36,2)=-5.d0*(3.d0*x**2-3.d0*y**2+3.d0*z**2-1.d0)/14.d0
    qaL(36,3)=(3.d0*(35.d0*x**2*z**2-5.d0*x**2-35.d0*y**2*z**2 +5.d0*y**2-35.d0*z**4+30.d0*z**2-3.d0))/14.d0
    qaL(37,1)=-(3.d0*sq5*(x**2-y**2))/(2.d0*sq3)
    qaL(38,1)=-(3.d0*sqrt(10.d0)*zi*y)/(5.d0*sqrt(2.d0))
    qaL(38,2)=-0.15d0*sq5*zi*y *(15.d0*x**2-5.d0*y**2+5.d0*z**2-1.d0)
    qaL(39,1)=0.d0
    qaL(39,2)=1.5d0*sq5*zi*z*(x**2-y**2)
    qaL(40,1)=(3.d0*sq5*zi*x)/5.d0
    qaL(40,2)=-0.15d0*sq5*zi*x  *(5.d0*x**2-15.d0*y**2-5.d0*z**2+1.d0)
    qaL(41,1)=0.d0
    qaL(41,2)=0.d0
    qaL(41,3)=7.5d0*x*y*(x**2-y**2)
    qaL(42,1)=0.d0
    qaL(42,2)=-(15.d0*y*z)/7.d0
    qaL(42,3)=-15.d0*y*z *(21.d0*x**2-7.d0*y**2+7.d0*z**2-3.d0)/28.d0
    qaL(43,1)=0.d0
    qaL(43,2)=(15.d0*(x**2-y**2))/(7.d0*sq3)
    qaL(43,3)=15.d0*sq3 *(7.d0*x**2*z**2-x**2-7.d0*y**2*z**2+y**2)/28.d0
    qaL(44,1)=0.d0
    qaL(44,2)=(15.d0*x*z)/7.d0
    qaL(44,3)=-15.d0*x*z *(7.d0*x**2-21.d0*y**2-7.d0*z**2+3.d0)/28.d0
    qaL(45,1)=1.d0
    qaL(45,2)=(5.d0*(3.d0*z**2-1.d0))/7.d0
    qaL(45,3)=3.d0*(35.d0*x**4-210.d0*x**2*y**2 +35.d0*y**4+35.d0*z**4-30.d0*z**2+3.d0)/56.d0

!--    THIS will be done when qaL is used.
!      qaL(:,:)=qaL(:,:)/( (4.d0*Pi)*(-zi)**L )

    return

  END SUBROUTINE get_qaL

END MODULE PSQRijPrep
