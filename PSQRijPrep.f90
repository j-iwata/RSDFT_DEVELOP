MODULE PSQRijPrep
use parallel_module, only: myrank
  use atom_module, only: Natom,ki_atom
  use ps_nloc2_variables
  use VarPSMember
  ! lo,Nelement_,norb,viod
  use VarPSMemberG
  use pseudopot_module, only: pselect
  implicit none
  PRIVATE
  PUBLIC :: prepQRijp102
  
  complex(8),allocatable :: qaL(:,:) !qaL(k3max,Lrefmax)

  complex(8),parameter :: z0=(0.d0,0.d0),z1=(1.d0,0.d0),zi=(0.d0,1.d0)

  real(8),allocatable :: y2a(:,:,:,:),

CONTAINS

  SUBROUTINE prepQRijp102
    implicit none
    real(8) :: x,y,z
    integer :: l,m

    integer :: i,ik,iorb
    integer :: NRc,n
    real(8) :: d1,d2

if (myrank==0) write(400+myrank,*) ">>>>> inside prepQRijp102"

    INTERFACE
       FUNCTION Ylm(x,y,z,l,m)
         real(8) :: Ylm
         real(8),intent(IN) :: x,y,z
         integer,intent(IN) :: l,m
       END FUNCTION Ylm
    END INTERFACE

  call allocateQaL

    Mqr=0
    do i=1,Natom
       ik=ki_atom(i)
       do iorb=1,norb(ik)
          Mlma=Mlma+2*lo(iorb,ik)+1
       end do
    end do

    if ( Mlma <= 0 ) return

    if ( .not.allocated(y2a) .and. pselect /=4 ) then
       NRc=maxval(NRps)
       n=maxval(norb)
       allocate( y2a(NRc,n,Nelement_) )
       y2a=0.d0
       do ik=1,Nelement_
       do 
          d1=0.d0
          d2=0.d0
          call spline(rad1(1,ik),qrL(1,iorb,ik),NRps(iorb,ik),d1,d2,y2a(1,iorb,ik))
       end do
       end do
    end if

    if ( Mlma < nprocs_g ) then
       nzlma_0 = Mlma
    else
       nzlma_0 = min(Mlma*125/nprocs_g,Mlma)
    end if

!    ctt(:)=0.d0
!    ett(:)=0.d0

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
    mm1 = maxval( abs(mcube_grid_ion(:,1)) ) + 1
    mm2 = maxval( abs(mcube_grid_ion(:,2)) ) + 1
    mm3 = maxval( abs(mcube_grid_ion(:,3)) ) + 1

    call watch(ctt(8),ett(8))

    MMJJ_0 = M_grid_ion

    if ( .not.allocated(icheck_tmp3) ) then
    L=maxval(lo)
    n=maxval(norb)
    allocate( icheck_tmp3(Natom,n,2*L+1) ) ; icheck_tmp3=0
    allocate( JJ_tmp(6,MMJJ_0,n,Natom)   ) ; JJ_tmp=0
    allocate( MJJ_tmp(n,Natom)           ) ; MJJ_tmp=0
    allocate( uV_tmp(MMJJ_0,n,Natom)     ) ; uV_tmp=0.d0
    end if

    call watch(ctt(0),ett(0))

    if ( pselect == 4 ) then

       call prep_ps_nloc_gth(Natom,n,L,MMJJ_0,M_grid_ion,map_grid_ion &
                            ,icheck_tmp3,JJ_tmp,MJJ_tmp,uV_tmp,nzlma,MMJJ)

    else

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

    c1                 = 1.d0/Ngrid(1)
    c2                 = 1.d0/Ngrid(2)
    c3                 = 1.d0/Ngrid(3)
    maxerr             = 0
    icheck_tmp3(:,:,:) = 0
    MMJJ               = 0
    nzlma              = 0
    lma                = 0
    lma0               = 0

!$OMP parallel do schedule(dynamic) firstprivate( maxerr ) &
!$OMP    private( Rx,Ry,Rz,ic1,ic2,ic3,ik,iorb,Rps2,NRc,L,j,i,i1,i2,i3 &
!$OMP            ,id1,id2,id3,k1,k2,k3,i1_0,i2_0,i3_0,d1,d2,d3,x,y,z,r2,r &
!$OMP            ,v0,err0,ir0,ir,mm,m1,m2,v,err )
    do a=1,Natom

       Rx = aa(1,1)*aa_atom(1,a)+aa(1,2)*aa_atom(2,a)+aa(1,3)*aa_atom(3,a)
       Ry = aa(2,1)*aa_atom(1,a)+aa(2,2)*aa_atom(2,a)+aa(2,3)*aa_atom(3,a)
       Rz = aa(3,1)*aa_atom(1,a)+aa(3,2)*aa_atom(2,a)+aa(3,3)*aa_atom(3,a)

       ic1 = nint( aa_atom(1,a)*Ngrid(1) )
       ic2 = nint( aa_atom(2,a)*Ngrid(2) )
       ic3 = nint( aa_atom(3,a)*Ngrid(3) )

       ik = ki_atom(a)

       do iorb=1,norb(ik)

          Rps2 = Rps(iorb,ik)**2
          NRc  = NRps(iorb,ik)
          L    = lo(iorb,ik)
          j    = 0

          do i=1,M_grid_ion

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

                d1 = id1*c1
                d2 = id2*c2
                d3 = id3*c3

                x  = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
                y  = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
                z  = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
                r2 = x*x+y*y+z*z

                if ( r2 > Rps2+1.d-10 ) cycle

                r    = sqrt(r2)
                v0   = 0.d0
                err0 = 0.d0

                if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 .or. &
                     abs(z)>1.d-14 .or. L==0 ) then
#ifdef _SPLINE_
                   call splint(rad1(1,ik),viod(1,iorb,ik),y2a,NRc,r,v0)
#else
                   ir0=irad( int(100.d0*r),ik )
                   do ir=ir0,NRc
                      if ( r<rad1(ir,ik) ) exit
                   end do
                   if ( ir <= 2 ) then
                      v0=viod(2,iorb,ik)
                      if ( ir < 1 ) stop "ps_nloc2(0)"
                   else if ( ir <= NRc ) then
                      err0=1.d10
                      do mm=1,20
                         m1=max(1,ir-mm)
                         m2=min(ir+mm,NRc)
                         call polint &
                              (rad1(m1,ik),viod(m1,iorb,ik),m2-m1+1,r,v,err)
                         if ( abs(err)<err0 ) then
                            v0=v
                            err0=abs(err)
                            if ( err0<ep ) exit
                         end if
                      end do
                   else
                      write(*,*) "ps_nloc2(1)",ir,NRc,viod(NRc,iorb,ik)
                      write(*,*) viod(NRc+1,iorb,ik),r,rad1(ir,ik)
                      stop
                   end if
                   maxerr=max(maxerr,err0)
#endif
                end if

                j=j+1
                JJ_tmp(1,j,iorb,a) = i1_0
                JJ_tmp(2,j,iorb,a) = i2_0
                JJ_tmp(3,j,iorb,a) = i3_0
                JJ_tmp(4,j,iorb,a) = k1
                JJ_tmp(5,j,iorb,a) = k2
                JJ_tmp(6,j,iorb,a) = k3
                uV_tmp(j,iorb,a)   = v0
                      
             end if

          end do ! i ( 1 - M_grid_ion )

          MJJ_tmp(iorb,a)=j

       end do ! iorb
    end do ! a
!$OMP end parallel do

#ifndef _SPLINE_
    deallocate( irad )
#endif

    end if ! pselect

    lma=0
    do a=1,Natom
       ik=ki_atom(a)
       do iorb=1,norb(ik)
          j=MJJ_tmp(iorb,a)
          if ( j > 0 ) then
             L=lo(iorb,ik)
             nzlma=nzlma+2*L+1
             do m=1,2*L+1
                lma=lma+1
                icheck_tmp3(a,iorb,m)=lma
             end do
          end if
       end do
    end do
    MMJJ = maxval( MJJ_tmp )

    allocate( lcheck_tmp1(Mlma,0:np_grid-1) )
    lcheck_tmp1(:,:)=.false.
    lma=0
    do a=1,Natom
       ik=ki_atom(a)
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          j=MJJ_tmp(iorb,a)
          do m=1,2*L+1
             lma=lma+1
             if ( j > 0 ) then
                lcheck_tmp1(lma,myrank_g)=.true.
             end if
          end do
       end do
    end do
    call mpi_allgather(lcheck_tmp1(1,myrank_g),Mlma,mpi_logical &
                      ,lcheck_tmp1,Mlma,mpi_logical,comm_grid,ierr)

    call watch(ctt(1),ett(1))

! for grid-parallel computation

    nzlma_0 = min(nzlma_0*2,Mlma)
    nrlma   = 0

    n=maxval( node_partition(1:3) )
    allocate( itmp(n,3) ) ; itmp=0
    allocate( lma_nsend_tmp(0:nprocs_g-1) )
    allocate( icheck_tmp1(0:nprocs_g-1)   )
    allocate( icheck_tmp2(0:nprocs_g-1)   )
    allocate( nl_rank_map_tmp(0:nprocs_g) )
    allocate( maps_tmp(nzlma_0,6) )
    allocate( sendmap_tmp(nzlma_0,0:nprocs_g-1) )
    allocate( recvmap_tmp(nzlma_0,0:nprocs_g-1) )

    maps_tmp(:,:)      = 0
    sendmap_tmp(:,:)   = 0
    recvmap_tmp(:,:)   = 0
    icheck_tmp1(:)     = 0
    icheck_tmp2(:)     = 0
    lma_nsend_tmp(:)   = 0
    nl_rank_map_tmp(:) =-1

    np1 = node_partition(1)
    np2 = node_partition(2)
    np3 = node_partition(3)

    lma=0
    do a=1,Natom
       ik=ki_atom(a)
    do iorb=1,norb(ik)
       L=lo(iorb,ik)
    do m=-L,L
       lma=lma+1

       icheck_tmp1(:)=0
       do n=0,np_grid-1
          if ( lcheck_tmp1(lma,n) ) icheck_tmp1(n) = 1
       end do
       icheck_tmp1(myrank_g) = icheck_tmp3(a,iorb,m+L+1)

       itmp(:,:)=0
       n=-1
       do i3=1,node_partition(3)
       do i2=1,node_partition(2)
       do i1=1,node_partition(1)
          n=n+1
          if ( icheck_tmp1(n) == 0 ) cycle
          itmp(i1,1) = i1
          itmp(i2,2) = i2
          itmp(i3,3) = i3
       end do
       end do
       end do
       k1=count( itmp(:,1)>0 )
       k2=count( itmp(:,2)>0 )
       k3=count( itmp(:,3)>0 )
 
#ifdef _REFACT_
       call #########

#else
       ic1=0
       id1=np1
       do i=1,np1
          if ( ic1==0 .and. itmp(i,1)/=0 ) then
             ic1=i
          else if ( ic1/=0 .and. itmp(i,1)==0 ) then
             id1=i-1
             exit
          end if
       end do
       if ( id1-ic1+1/=k1 ) then
          i1=0
          j1=np1
          do i=id1+1,np1
             if ( i1==0 .and. itmp(i,1)/=0 ) then
                i1=i
             else if ( i1/=0 .and. itmp(i,1)==0 ) then
                j1=i-1
                exit
             end if
          end do
          i1=i1-np1
          j1=j1-np1
          ic1=i1
       end if
       ic2=0
       id2=np2
       do i=1,np2
          if ( ic2==0 .and. itmp(i,2)/=0 ) then
             ic2=i
          else if ( ic2/=0 .and. itmp(i,2)==0 ) then
             id2=i-1
             exit
          end if
       end do
       if ( id2-ic2+1/=k2 ) then
          i2=0
          j2=np2
          do i=id2+1,np2
             if ( i2==0 .and. itmp(i,2)/=0 ) then
                i2=i
             else if ( i2/=0 .and. itmp(i,2)==0 ) then
                j2=i-1
                exit
             end if
          end do
          i2=i2-np2
          j2=j2-np2
          ic2=i2
       end if
       ic3=0
       id3=np3
       do i=1,np3
          if ( ic3==0 .and. itmp(i,3)/=0 ) then
             ic3=i
          else if ( ic3/=0 .and. itmp(i,3)==0 ) then
             id3=i-1
             exit
          end if
       end do
       if ( id3-ic3+1/=k3 ) then
          i3=0
          j3=np3
          do i=id3+1,np3
             if ( i3==0 .and. itmp(i,3)/=0 ) then
                i3=i
             else if ( i3/=0 .and. itmp(i,3)==0 ) then
                j3=i-1
                exit
             end if
          end do
          i3=i3-np3
          j3=j3-np3
          ic3=i3
       end if
       do j3=ic3,id3
       do j2=ic2,id2
       do j1=ic1,id1
          k1=mod(j1+np1-1,np1)+1
          k2=mod(j2+np2-1,np2)+1
          k3=mod(j3+np3-1,np3)+1
          k = k1-1 + (k2-1)*np1 + (k3-1)*np1*np2
          if ( icheck_tmp1(k)==0 ) icheck_tmp1(k)=-1
       end do
       end do
       end do
       do n=0,nprocs_g-1
          if ( icheck_tmp1(n)/=0 ) then
             icheck_tmp2(n)=icheck_tmp2(n)+1
          end if
       end do
#endif

       if ( icheck_tmp1(myrank_g)/=0 ) then
          if ( icheck_tmp1(myrank_g)>0 ) then
             maps_tmp(icheck_tmp2(myrank_g),1)=icheck_tmp1(myrank_g)
          end if
          maps_tmp(icheck_tmp2(myrank_g),2)=inorm(iorb,ik)
          maps_tmp(icheck_tmp2(myrank_g),3)=a
          maps_tmp(icheck_tmp2(myrank_g),4)=L
          maps_tmp(icheck_tmp2(myrank_g),5)=m
          maps_tmp(icheck_tmp2(myrank_g),6)=iorb

          do n=0,nprocs_g-1
             if ( n==myrank_g .or. icheck_tmp1(n)==0 ) cycle
             lma_nsend_tmp(n)=lma_nsend_tmp(n)+1  
             sendmap_tmp(lma_nsend_tmp(n),n)=icheck_tmp2(myrank_g)
             recvmap_tmp(lma_nsend_tmp(n),n)=icheck_tmp2(n)
             if ( any(nl_rank_map_tmp(0:nrlma)==n) ) cycle
             nrlma=nrlma+1
             nl_rank_map_tmp(nrlma)=n
          end do
       end if

    end do ! m
    end do ! iorb
    end do ! a

    call watch(ctt(2),ett(2))

    nzlma = icheck_tmp2(myrank_g)

    deallocate( itmp )
    deallocate( icheck_tmp2 )
    deallocate( icheck_tmp1 )
!    deallocate( icheck_tmp3 )
    deallocate( lcheck_tmp1 )

    if ( allocated(uV) ) then
       deallocate( uV )
       deallocate( JJ_MAP )
       deallocate( MJJ_MAP )
       deallocate( MJJ )
       deallocate( iuV )
       deallocate( amap )
       deallocate( lmap )
       deallocate( mmap )
       deallocate( iorbmap )
       deallocate( nl_rank_map )
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
    allocate( nl_rank_map(nrlma)   ) ; nl_rank_map=-1

    do i=1,nrlma
       nl_rank_map(i)=nl_rank_map_tmp(i)
    end do

    deallocate( nl_rank_map_tmp )

    do lma=1,nzlma
       if ( maps_tmp(lma,1) == 0 ) cycle
       iuV(lma)     = maps_tmp(lma,2)
       amap(lma)    = maps_tmp(lma,3)
       lmap(lma)    = maps_tmp(lma,4)
       mmap(lma)    = maps_tmp(lma,5)
       iorbmap(lma) = maps_tmp(lma,6)
    end do

    c1=1.d0/ML1
    c2=1.d0/ML2
    c3=1.d0/ML3

!$OMP parallel do private( a,l,m,iorb,Rx,Ry,Rz,j,i1,i2,i3,k1,k2,k3,d1,d2,d3,x,y,z )
    do lma=1,nzlma
       if ( maps_tmp(lma,1) == 0 ) cycle
       a    = amap(lma)
       l    = lmap(lma)
       m    = mmap(lma)
       iorb = iorbmap(lma)
       MJJ_MAP(lma) = MJJ_tmp(iorb,a)
       Rx=aa(1,1)*aa_atom(1,a)+aa(1,2)*aa_atom(2,a)+aa(1,3)*aa_atom(3,a)
       Ry=aa(2,1)*aa_atom(1,a)+aa(2,2)*aa_atom(2,a)+aa(2,3)*aa_atom(3,a)
       Rz=aa(3,1)*aa_atom(1,a)+aa(3,2)*aa_atom(2,a)+aa(3,3)*aa_atom(3,a)
       do j=1,MJJ_MAP(lma)
          i1=JJ_tmp(1,j,iorb,a)
          i2=JJ_tmp(2,j,iorb,a)
          i3=JJ_tmp(3,j,iorb,a)
          k1=JJ_tmp(4,j,iorb,a)
          k2=JJ_tmp(5,j,iorb,a)
          k3=JJ_tmp(6,j,iorb,a)
          d1=c1*i1+k1
          d2=c2*i2+k2
          d3=c3*i3+k3
          x = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
          y = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
          z = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
          uV(j,lma) = uV_tmp(j,iorb,a)*Ylm(x,y,z,l,m)
          JJ_MAP(1:6,j,lma) = JJ_tmp(1:6,j,iorb,a)
       end do
    end do
!$OMP end parallel do

!    deallocate( MJJ_tmp )
!    deallocate( JJ_tmp )
!    deallocate( uV_tmp )
    deallocate( maps_tmp )

    call watch(ctt(3),ett(3))

    allocate( icheck_tmp4(a1b:b1b,a2b:b2b,a3b:b3b) )
    icheck_tmp4=0
    do lma=1,nzlma
       j=0
       icheck_tmp4=0
       do i=1,MJJ_MAP(lma)
          i1=JJ_MAP(1,i,lma)
          i2=JJ_MAP(2,i,lma)
          i3=JJ_MAP(3,i,lma)
          if ( icheck_tmp4(i1,i2,i3)==0 ) then
             j=j+1
             icheck_tmp4(i1,i2,i3)=j
          end if
       end do
       MJJ(lma)=j
    end do
    MAXMJJ = maxval( MJJ(1:nzlma) )
    deallocate( icheck_tmp4 )

    nl_max_send = maxval( lma_nsend_tmp )

    if ( allocated(lma_nsend) ) then
       deallocate( lma_nsend )
       deallocate( sendmap )
       deallocate( recvmap )
    end if
    allocate( lma_nsend(0:nprocs_g-1) ) ; lma_nsend=0
    allocate( sendmap(nl_max_send,0:nprocs_g-1) ) ; sendmap=0
    allocate( recvmap(nl_max_send,0:nprocs_g-1) ) ; recvmap=0

    do n=0,nprocs_g-1
       sendmap(1:nl_max_send,n) = sendmap_tmp(1:nl_max_send,n)
       lma_nsend(n) = lma_nsend_tmp(n)
    end do


    allocate( ireq(2*nprocs_g) )
    allocate( istatus(MPI_STATUS_SIZE,2*nprocs_g) )
    nreq=0
    do n=0,nprocs_g-1
       if ( lma_nsend(n)<=0 .or. n==myrank_g ) cycle
       nreq=nreq+1
       call mpi_isend(recvmap_tmp(1,n),lma_nsend(n),mpi_integer,n,1 &
            ,comm_grid,ireq(nreq),ierr)
       nreq=nreq+1
       call mpi_irecv(recvmap(1,n) ,lma_nsend(n),mpi_integer,n,1 &
            ,comm_grid,ireq(nreq),ierr)
    end do
    call mpi_waitall(nreq,ireq,istatus,ierr)
    deallocate( istatus )
    deallocate( ireq )

    deallocate( recvmap_tmp,sendmap_tmp,lma_nsend_tmp )

    call prepThreeWayComm( nrlma,nl_rank_map,nrlma_xyz,num_2_rank )

    call watch(ctt(4),ett(4))

    call allocate_ps_nloc2(MB_d)

    if ( allocated(JJP) ) then
       deallocate( JJP )
       deallocate( uVk )
    end if
    allocate( JJP(MAXMJJ,nzlma) ) ; JJP=0
    allocate( uVk(MAXMJJ,nzlma,MBZ_0:MBZ_1) ) ; uVk=0.d0

    call prep_uvk_ps_nloc2(MBZ_0,MBZ_1,kbb(1,MBZ_0))

    call watch(ctt(5),ett(5))

    if ( disp_switch_parallel ) then
       write(*,*) "time(ps_nloc2_1)",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(ps_nloc2_2)",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(ps_nloc2_3)",ctt(3)-ctt(2),ett(3)-ett(2)
       write(*,*) "time(ps_nloc2_4)",ctt(4)-ctt(3),ett(4)-ett(3)
       write(*,*) "time(ps_nloc2_5)",ctt(5)-ctt(4),ett(5)-ett(4)
       write(*,*) "time(ps_nloc2_7)",ctt(7)-ctt(6),ett(7)-ett(6)
       write(*,*) "time(ps_nloc2_8)",ctt(8)-ctt(7),ett(8)-ett(7)
       write(*,*) "time(ps_nloc2_9)",ctt(0)-ctt(8),ett(0)-ett(8)
    end if
    
if (myrank==0) write(400+myrank,*) "<<<<< end of prepQRijp102"
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
