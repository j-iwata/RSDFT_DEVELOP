MODULE PSQRijPrep

  use parallel_module, only: myrank,myrank_g,nprocs_g,np_grid,node_partition,COMM_GRID,disp_switch_parallel,MB_d
  use aa_module, only: aa
  use atom_module, only: Natom,ki_atom,aa_atom
  use rgrid_module, only: Igrid,Ngrid,Hgrid
  use VarPSMember
  ! lo,Nelement_,norb,viod
  use VarPSMemberG
  use VarParaPSnonLocG
  use ps_nloc2_module, only: prepMapsTmp
  use pseudopot_module, only: pselect
  use minimal_box_module
  use para_rgrid_comm
  use watch_module
  use array_bound_module, only: ML_0
  use polint_module
  use spline_module

  implicit none

  include 'mpif.h'
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
!    allocate( qaL(k2max,max_Lref) ) ; qaL=z0
    allocate( qaL(45,3) ) ; qaL=z0
    return
  END SUBROUTINE allocateQaL


  SUBROUTINE prepQRijp102
    implicit none
    real(8) :: x,y,z
    integer :: l,m

    integer :: ia,ik,ik1,ik2,ik3,il,ir,ir0
    integer :: i,j,i1,i2,i3
    integer :: iorb,iorb1,iorb2
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
    integer,allocatable :: JJ_tmp(:,:,:,:),MJJ_tmp_Q(:,:)
    real(8),allocatable :: uV_tmp(:,:,:)

    integer :: M_irad
    integer,allocatable :: irad(:,:)

    real(8) :: c1,c2,c3
    real(8) :: maxerr
    integer :: c_nzqr_pre

    real(8) :: Rx,Ry,Rz
    integer :: ic1,ic2,ic3,id1,id2,id3
    integer :: k1,k2,k3,i1_0,i2_0,i3_0
    real(8) :: d1,d2,d3
    real(8) :: r2

    integer :: ll3,mm,m1,m2,iqr
    integer :: j3
    real(8) :: v,v0,err,err0,QRtmp

    logical,allocatable :: lcheck_tmp1(:,:)

    real(8),allocatable :: QRij_tmp(:,:,:)
    real(8),allocatable :: QRij_tmp2(:,:)
    integer :: c_nzqr_0
    integer,allocatable :: itmp(:,:),icheck_tmp1(:),icheck_tmp2(:),icheck_tmp4(:,:,:)
    integer,allocatable :: qr_nsend_tmp(:),nl_rank_map_tmp(:),maps_tmp(:,:)
    integer,allocatable :: sendmap_tmp(:,:),recvmap_tmp(:,:)
    logical,allocatable :: isInThisNode_Q(:,:)
    
    integer :: ierr,nreq
    real(8) :: ctt(0:9),ett(0:9)
    integer,allocatable :: ireq(:)
    integer,allocatable :: istatus(:,:)

    integer :: kk1

    real(8),parameter :: ep=1.d-8

#ifdef _SHOWALL_INIT_
    write(200+myrank,*) ">>>>> prepQRijp102"
#endif

    ctt=0.0d0 ; ett=0.0d0
    call watch(ctt(6),ett(6))

    call allocateQaL

#ifdef _SPLINE_
    if ( .not.allocated(y2a) ) then
      allocate( y2a(max_psgrd,max_Lref,max_k2,Nelement_) )
      y2a=0.d0
      do ik=1,Nelement_
        do ik2=1,N_k2(ik)
          iorb  = k2_to_iorb(2,ik2,ik)
          do il=1,nl3v(ik2,ik)
            d1=0.d0
            d2=0.d0
            call spline(rad1(1,ik),qrL(1,il,ik2,ik),NRps(iorb,ik),d1,d2,y2a(1,il,ik2,ik))
          end do
        end do
      end do
    end if
!    if (myrank==0) write(200,*) 'spline finished'
#endif

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

    nl3vmax = maxval(nl3v)
    if ( .not.allocated(icheck_tmp5) ) allocate( icheck_tmp5(k1max,Natom)     ) ; icheck_tmp5=0
    if ( .not.allocated(JJ_tmp)      ) allocate( JJ_tmp(6,MMJJ_0,k1max,Natom) ) ; JJ_tmp     =0
    if ( .not.allocated(MJJ_tmp_Q)   ) allocate( MJJ_tmp_Q(k1max,Natom)       ) ; MJJ_tmp_Q  =0
    if ( .not.allocated(QRij_tmp)    ) allocate( QRij_tmp(MMJJ_0,k1max,Natom) ) ; QRij_tmp   =0.d0

    call watch(ctt(0),ett(0))
    
#ifndef _SPLINE_
    allocate( irad(0:3000,Nelement_) ) ; irad=0
    M_irad=0
    do ik=1,Nelement_
       NRc  = maxval( NRps(:,ik) )
       NRc  = min( 3000, NRc )
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
    n=maxval(N_k1)
    allocate( isInThisNode_Q(1:Natom,1:n) ) ; isInThisNode_Q=.false.

!call write_qrL(5400,myrank,Nelement_)

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
#ifdef _SHOWALL_MAP_Q_
        write(530+myrank,*) repeat('-',30)
        write(530+myrank,*) "ia,ik,ik1=",ia,ik,ik1
#endif
        j=0
        ik2 = k1_to_k2(ik1,ik)
        ik3 = k1_to_k3(ik1,ik)
        iorb1 = k1_to_iorb(1,ik1,ik)
        iorb2 = k1_to_iorb(2,ik1,ik)
        Rps2  = max(Rps(iorb1,ik)**2,Rps(iorb2,ik)**2)
        NRc   = max(NRps(iorb1,ik),NRps(iorb2,ik))
          do i=1,MMJJ_0
            i1 = map_grid_ion(1,i)
            i2 = map_grid_ion(2,i)
            i3 = map_grid_ion(3,i)
            id1 = ic1 + i1
            id2 = ic2 + i2
            id3 = ic3 + i3
            k1  = id1/ML1 ; if ( id1<0 ) k1=(id1+1)/ML1-1
            k2  = id2/ML2 ; if ( id2<0 ) k2=(id2+1)/ML2-1
            k3  = id3/ML3 ; if ( id3<0 ) k3=(id3+1)/ML3-1
            i1_0  = id1-k1*ML1
            i2_0  = id2-k2*ML2
            i3_0  = id3-k3*ML3
            if ( Igrid(1,1) <= i1_0 .and. i1_0 <= Igrid(2,1) .and. &
                Igrid(1,2) <= i2_0 .and. i2_0 <= Igrid(2,2) .and. &
                Igrid(1,3) <= i3_0 .and. i3_0 <= Igrid(2,3) ) then
              isInThisNode_Q(ia,ik1)=.true.
              d1 = id1*c1
              d2 = id2*c2
              d3 = id3*c3
              x  = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
              y  = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
              z  = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
              r2 = x*x+y*y+z*z
                if ( r2 > Rps2+1.d-10 ) cycle
              j=j+1
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
                  call splint(rad1(1,ik),qrL(1,ll3,ik2,ik),y2a(1,ll3,ik2,ik),NRc,r,v0)
#else
                  ir0=irad( int(100.d0*r),ik )
                  do ir=ir0,NRc
                    if ( r<rad1(ir,ik) ) exit
                  end do
                  if ( ir <= 2 ) then
                    v0  = qrL(2,ll3,ik2,ik)
                    if ( ir < 1 ) stop "prep_QRij_p102(0)"
!!!!!                  else if ( ir <= NRc ) then
                  else if ( ir <= NRc ) then
                    err0=1.d10
!if (ik==2) write(5900+myrank,*) repeat('-',40) 
                    do mm=1,20
                      m1=max(1,ir-mm)
                      m2=min(ir+mm,NRc)
                      call polint(rad1(m1,ik),qrL(m1,ll3,ik2,ik),m2-m1+1,r,v,err)
!if (ik==2) write(5900+myrank,'(3I6,5G20.7)') mm,m1,m2,rad1(m1,ik),qrL(m1,ll3,ik2,ik),v,err,err0
                      if ( abs(err)<err0 ) then
                        v0=v
                        err0=abs(err)
                        if ( err0<ep ) exit
                      end if
                    end do
!if (ik==2) write(5900+myrank,'(I6,G20.7)') ik,v0
                  else
                    write(*,*) "prep_QRij_p102(1)",ir,NRc,qrL(NRc,ll3,ik2,ik)
                    write(*,*) qrL(NRc+1,ll3,ik3,ik),r,rad1(ir,ik)
                    stop
                  end if
                  maxerr=max(maxerr,err0)
#endif
                  QRtmp = v0/pi4* real(qaL(ik3,ll3)/(-zi)**L)
!write(5900+myrank,'(6I6,5G20.7,I6)') ia,ik1,ik2,ll3,ir0,ir,x,y,z,r,v0,NRc
!write(5200+myrank,'(4I6,5G20.7)') ia,ik1,j,ll3,QRtmp,qaL(ik3,ll3),v0
                end if ! x,y,z
!                if (abs(QRtmp)<1.d-10) cycle

                QRij_tmp(j,ik1,ia) = QRij_tmp(j,ik1,ia)+QRtmp
              end do ! ll3
!              if (abs(QRij_tmp(j,ik1,ia))<1.d-10) then
!                if (myrank==0) write(200,*) "some points are deleted"
!                QRij_tmp(j,ik1,ia)=0.d0
!                j=j-1
!                cycle
!              endif
              JJ_tmp(1,j,ik1,ia) = i1_0
              JJ_tmp(2,j,ik1,ia) = i2_0
              JJ_tmp(3,j,ik1,ia) = i3_0
              JJ_tmp(4,j,ik1,ia) = k1
              JJ_tmp(5,j,ik1,ia) = k2
              JJ_tmp(6,j,ik1,ia) = k3
#ifdef _SHOWALL_MAP_Q_
              write(530+myrank,'(2I5,A8,6I4)') i,j," JJ_tmp=",JJ_tmp(1:6,j,ik1,ia)
#endif
            end if ! Igrid
          end do ! i ( 1 - MMJJ_0 )
          MJJ_tmp_Q(ik1,ia)   = j
#ifdef _SHOWALL_MAP_Q_
          write(530+myrank,*) "MJJ_tmp_Q(ik1,ia)=",MJJ_tmp_Q(ik1,ia)
#endif
       end do ! ik1
    end do ! ia
!!$OMP end parallel do
    MMJJ_Q = maxval( MJJ_tmp_Q )

#ifndef _SPLINE_
    deallocate( irad )
#endif
    
    c_nzqr_pre=0
    Mqr=0
    do ia=1,Natom
      ik=ki_atom(ia)
      do ik1=1,N_k1(ik)
        Mqr=Mqr+1
        j=MJJ_tmp_Q(ik1,ia)
        if (j>0) then
!        if ( isInThisNode_Q(ia,ik1) ) then
          c_nzqr_pre = c_nzqr_pre+1
          icheck_tmp5(ik1,ia) = c_nzqr_pre
        endif
      end do
    end do
    deallocate( isInThisNode_Q )

    allocate( lcheck_tmp1(Mqr,0:np_grid-1) )
    lcheck_tmp1(:,:)=.false.
    iqr=0
    do ia=1,Natom
      ik  = ki_atom(ia)
      do ik1=1,N_k1(ik)
        iqr=iqr+1
        j=MJJ_tmp_Q(ik1,ia)
        if (j>0) then
          lcheck_tmp1(iqr,myrank_g) = .true.
        end if
      end do
    end do
    call mpi_allgather(lcheck_tmp1(1,myrank_g),Mqr,mpi_logical,lcheck_tmp1,Mqr,mpi_logical,comm_grid,ierr)
    
    call watch(ctt(1),ett(1))

! for grid-parallel computation

    c_nzqr_0 = N_nzqr
#ifdef _SHOWALL_MAP_Q_
    write(1100+myrank,*) 'c_nzqr=    ',c_nzqr_0
    write(1100+myrank,*) 'N_k1(max)= ',maxval(N_k1)
#endif

    n=maxval( node_partition(1:3) )
    allocate( itmp(n,3)                   ) ; itmp            =0
    allocate( qr_nsend_tmp(0:nprocs_g-1)  ) ; qr_nsend_tmp    =0
    allocate( icheck_tmp1(0:nprocs_g-1)   ) ; icheck_tmp1     =0
    allocate( icheck_tmp2(0:nprocs_g-1)   ) ; icheck_tmp2     =0
    allocate( nl_rank_map_tmp(0:nprocs_g) ) ; nl_rank_map_tmp =-1
    allocate( maps_tmp(c_nzqr_0,5)        ) ; maps_tmp        =0
    allocate( sendmap_tmp(c_nzqr_0,0:nprocs_g-1) ) ; sendmap_tmp =0
    allocate( recvmap_tmp(c_nzqr_0,0:nprocs_g-1) ) ; recvmap_tmp =0

    np1 = node_partition(1)
    np2 = node_partition(2)
    np3 = node_partition(3)

#ifdef _SHOWALL_MAP_Q_
    write(1100+myrank,*) repeat('-',50),'Q start'
    write(1100+myrank,'(2A5,A12)') 'ia','ik1','icheck_tmp5'
    do ia=1,Natom
      do ik1=1,k1max
        write(1100+myrank,'(2I5,I12)') ia,ik1,icheck_tmp5(ik1,ia)
      enddo
    enddo
    write(1100+myrank,*) repeat('=',50),'Q end'
#endif

#ifdef _SHOWALL_MAP_Q_
    write(1100+myrank,*) repeat('-',50),'icheck_tmp2'
#endif
    nrqr=0
    iqr=0
    do ia=1,Natom
      ik  = ki_atom(ia)
      do ik1=1,N_k1(ik)
        iqr=iqr+1
        kk1=kk1map(ik1,ia)
        icheck_tmp1(:)=0
!        call MPI_ALLGATHER(icheck_tmp5(ik1,ia),1,MPI_INTEGER,icheck_tmp1,1,MPI_INTEGER,COMM_GRID,ierr)
        do n=0,np_grid-1
          if ( lcheck_tmp1(iqr,n) ) icheck_tmp1(n) = 1
        end do
        icheck_tmp1(myrank_g) = icheck_tmp5(ik1,ia)

        call prepMapsTmp(np1,np2,np3,nprocs_g,itmp,icheck_tmp1,icheck_tmp2)

#ifdef _SHOWALL_MAP_Q_
        write(1100+myrank,'(2I4," icheck_tmp2(myrank_g)= ",I5)') ia,ik1,icheck_tmp2(myrank_g)
#endif

        if ( icheck_tmp1(myrank_g)/=0 ) then
          if ( icheck_tmp1(myrank_g)>0 ) then
            maps_tmp(icheck_tmp2(myrank_g),1)=icheck_tmp1(myrank_g)
          end if
          maps_tmp(icheck_tmp2(myrank_g),2) = ia
          maps_tmp(icheck_tmp2(myrank_g),3) = ik1
          maps_tmp(icheck_tmp2(myrank_g),4) = nzqr_pair(kk1,1)
          maps_tmp(icheck_tmp2(myrank_g),5) = nzqr_pair(kk1,2)
          do n=0,nprocs_g-1
            if ( n==myrank_g .or. icheck_tmp1(n)==0 ) cycle
            qr_nsend_tmp(n)=qr_nsend_tmp(n)+1  
            sendmap_tmp(qr_nsend_tmp(n),n)  = icheck_tmp2(myrank_g)
            recvmap_tmp(qr_nsend_tmp(n),n)  = icheck_tmp2(n)
            if ( any(nl_rank_map_tmp(0:nrqr)==n) ) cycle
            nrqr=nrqr+1
            nl_rank_map_tmp(nrqr)=n
          end do
        else ! [ icheck_tmp1(myrank_g) == 0 ]
           if ( all(icheck_tmp1 == 0) ) then
              do n=0,nprocs_g-1
                 icheck_tmp2(n) = icheck_tmp2(n) + 1
              end do ! n
           end if
        end if
      end do ! ik
    end do ! ia
    
    call watch(ctt(2),ett(2))

    c_nzqr = icheck_tmp2(myrank_g)

    if ( c_nzqr /= N_nzqr ) then
       write(*,*) "N_nzqr==c_nzqr is necessary"
       write(*,*) "N_nzqr,c_nzqr,myrank=",N_nzqr,c_nzqr,myrank
       stop "stop@prepQRijp102(PSQRijPrep.f90)"
    end if

    deallocate( itmp        )
    deallocate( icheck_tmp1 )
    deallocate( icheck_tmp2 )
    deallocate( icheck_tmp5 )
    deallocate( lcheck_tmp1 )

!===================================================================
    allocate( QRij_tmp2(MMJJ_Q,c_nzqr) ) ; QRij_tmp2  =0.d0

    call allocateJJMAPQ
    ! JJ_MAP_Q,MJJ_MAP_Q,MJJ_Q,nl_rank_map_Q
!===================================================================

    do i=1,nrqr
       nl_rank_map_Q(i)=nl_rank_map_tmp(i)
    end do
    deallocate( nl_rank_map_tmp )

    do iqr=1,c_nzqr
      if (maps_tmp(iqr,1)==0) cycle
      amap_Q(iqr) =maps_tmp(iqr,2)
      k1map_Q(iqr)=maps_tmp(iqr,3)
      lmamap_Q(iqr,1)=maps_tmp(iqr,4)
      lmamap_Q(iqr,2)=maps_tmp(iqr,5)
    end do

!!$OMP parallel do private( a,l,m,iorb,Rx,Ry,Rz,j,i1,i2,i3,k1,k2,k3,d1,d2,d3,x,y,z )
    do iqr=1,c_nzqr
      l   = maps_tmp(iqr,1)
      if (l==0) cycle
      ia  = maps_tmp(iqr,2)
      ik1 = maps_tmp(iqr,3)
      MJJ_MAP_Q(iqr) = MJJ_tmp_Q(ik1,ia)
#ifdef _SHOWALL_MAP_Q_
      write(530+myrank,'(4A5)') 'ia','ik1','l','MJJ_MAP_Q'
      write(530+myrank,'(4I5)') ia,ik1,l,MJJ_MAP_Q(iqr)
#endif
!write(6100+myrank,'(5A6)') 'iqr','ia','ik1','l','MJJ_MAP_Q'
!write(6100+myrank,'(5I6)') iqr,ia,ik1,l,MJJ_MAP_Q(iqr)
      do j=1,MJJ_MAP_Q(iqr)
        QRij_tmp2(j,iqr)    = QRij_tmp(j,ik1,ia)
        JJ_MAP_Q(1:6,j,iqr) = JJ_tmp(1:6,j,ik1,ia)
!write(6000+myrank,'(2I6,G20.7,6I4)') iqr,j,QRij_tmp2(j,iqr),JJ_tmp(1:6,j,ik1,ia)
#ifdef _SHOWALL_MAP_Q_
        write(520+myrank,'(I5,A6,E15.7e2,A8,6I4)') j," QRij=",QRij_tmp2(j,iqr)," JJ_tmp=",JJ_tmp(1:6,j,ik1,ia)
#endif
      end do
    end do
!!$OMP end parallel do

    deallocate( MJJ_tmp_Q )
    deallocate( JJ_tmp    )
    deallocate( QRij_tmp  )
    deallocate( maps_tmp  )

    call watch(ctt(3),ett(3))

    allocate( icheck_tmp4(a1b:b1b,a2b:b2b,a3b:b3b) ) ; icheck_tmp4=0
#ifdef _SHOWALL_MAP_Q_
    write(530+myrank,'(6I3)') a1b,b1b,a2b,b2b,a3b,b3b
#endif
    do iqr=1,c_nzqr
       j=0
       icheck_tmp4=0
       do i=1,MJJ_MAP_Q(iqr)
          i1  = JJ_MAP_Q(1,i,iqr)
          i2  = JJ_MAP_Q(2,i,iqr)
          i3  = JJ_MAP_Q(3,i,iqr)
#ifdef _SHOWALL_MAP_Q_
          write(530+myrank,'(3I3)') i1,i2,i3
#endif
          if ( icheck_tmp4(i1,i2,i3)==0 ) then
             j=j+1
             icheck_tmp4(i1,i2,i3)  = j
          end if
       end do
       MJJ_Q(iqr) = j
    end do
    MAXMJJ_Q = maxval( MJJ_Q(1:c_nzqr) )
!===================================================================
    if ( allocated(QRij) ) then
       deallocate( QRij      )
    end if
    allocate( QRij(MAXMJJ_Q,c_nzqr) ) ; QRij  =0.d0
    if ( allocated(JJP_Q) ) deallocate( JJP_Q )
    allocate( JJP_Q(MAXMJJ_Q,c_nzqr) ) ; JJP_Q=0
!===================================================================
    do iqr=1,c_nzqr
       j=0
       icheck_tmp4=0
       do i=1,MJJ_MAP_Q(iqr)
          i1  = JJ_MAP_Q(1,i,iqr)
          i2  = JJ_MAP_Q(2,i,iqr)
          i3  = JJ_MAP_Q(3,i,iqr)
          QRtmp = QRij_tmp2(i,iqr)
          if ( icheck_tmp4(i1,i2,i3)==0 ) then
            j=j+1
            icheck_tmp4(i1,i2,i3)  = j
            QRij(j,iqr) = QRtmp
            JJP_Q(j,iqr) = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ML_0
!write(6200+myrank,'(4I6,G20.7)') iqr,i,j,JJP_Q(j,iqr),QRij(j,iqr)
          else
            j3  = icheck_tmp4(i1,i2,i3)
            QRij(j3,iqr) = QRij(j3,iqr) + QRtmp
!write(6200+myrank,'(4I6,G20.7)') iqr,i,j3,JJP_Q(j3,iqr),QRij(j3,iqr)
          end if
       end do
    end do

    deallocate( QRij_tmp2   )
    deallocate( icheck_tmp4 )

    nl_max_send_Q = maxval( qr_nsend_tmp )
    call allocateMAPQ(nl_max_send_Q,nprocs_g)
    ! qr_nsend,sendmap_Q,recvmap_Q

    do n=0,nprocs_g-1
       sendmap_Q(1:nl_max_send_Q,n) = sendmap_tmp(1:nl_max_send_Q,n)
       qr_nsend(n)                = qr_nsend_tmp(n)
    end do

    allocate( ireq(2*nprocs_g) ) ; ireq=0
    allocate( istatus(MPI_STATUS_SIZE,2*nprocs_g) ) ; istatus=0

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
    nl_max_send_Q=maxval(qr_nsend)    ! NEED MORE CONSIDERATION!!!!!!!!!
    call allocateBufQ(nl_max_send_Q,nprocs_g,MB_d)
!===================================================================

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
    
#ifdef _SHOWALL_INIT_
    write(400+myrank,*) "<<<<< prepQRijp102"
#endif
    return
  END SUBROUTINE prepQRijp102

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
