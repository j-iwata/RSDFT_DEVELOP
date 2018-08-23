module force_ps_nloc2_module

  use array_bound_module, only: MB_0,MB_1,MBZ_0,MBZ_1,MSP_0,MSP_1,ML_0
  use ps_nloc2_variables, only: nzlma, MMJJ, mmap, lmap, amap &
       , iorbmap, JJ_MAP, MJJ_MAP, JJP, MJJ, uVk, nrlma_xyz, iuV, lma_nsend &
       , num_2_rank, recvmap, sendmap, sbufnl, rbufnl, Mlma, TYPE_MAIN, zero
  use var_ps_member, only: ps_type, norb, lo, rad1, dviod, ippform, NRps
  use pseudopot_module, only: pselect
  use parallel_module, only: MB_d, comm_grid
  use atom_module, only: aa_atom, ki_atom, Natom, Nelement
  use aa_module, only: aa
  use force_sub_sub_module, only: gaunt, construct_gaunt_coef_l1
  use wf_module, only: occ, unk
  use ylm_module, only: Ylm
  use rgrid_variables, only: dV, Ngrid, Igrid
  use watch_module, only: watchb
  use spline_module, only: splint, spline
  use ps_nloc2_init_module, only: ps_nloc2_init_derivative
  use ps_nloc_gth_module, only: init_force_ps_nloc_gth
  use bz_module, only: kbb
  use polint_module, only: polint

  implicit none

  PRIVATE
  PUBLIC :: calc_force_ps_nloc2_hp

  type(gaunt) :: yyy
  real(8),allocatable :: y2b(:,:,:)

contains

  subroutine calc_force_ps_nloc2_hp(MI,force2)
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force2(3,MI)
    integer :: i1,i2,i3
    integer :: i,j,k,s,n,ir,iorb,L,L1,L1z,NRc,irank,jrank
    integer :: nreq,max_nreq,ib,ib1,ib2,nnn
    integer :: a,a0,ik,m,lm0,lm1,lma,im,m1,m2
    integer :: ierr,M_irad,ir0,ir1
    integer,allocatable :: ireq(:),istatus(:,:),irad(:,:),ilm1(:,:,:)
    real(8),parameter :: ep=1.d-8
    real(8),save :: Y1(0:3,-3:3,0:4,-4:4)
    real(8),save :: Y2(0:3,-3:3,0:4,-4:4)
    real(8),save :: Y3(0:3,-3:3,0:4,-4:4)
    real(8) :: err,err0,maxerr,Rx,Ry,Rz
    real(8) :: a1,a2,a3,c1,c2,c3,d1,d2,d3
    real(8) :: x,y,z,r,kr,pi2,c
    real(8) :: tmp,tmp0,tmp1
    real(8) :: ttmp(2), tttt(2,12)
    real(8) :: yy1,yy2,yy3
    real(8),allocatable :: work2(:,:),duVdR(:,:,:)
#ifdef _DRSDFT_
    real(8) :: ztmp
    real(8),allocatable :: wtmp5(:,:,:,:,:),vtmp2(:,:,:)
#else
    complex(8) :: ztmp
    complex(8),allocatable :: wtmp5(:,:,:,:,:),vtmp2(:,:,:)
#endif
    logical,save :: flag_Y = .true.
    logical,allocatable :: a_rank(:)
    integer :: ML1,ML2,ML3,i0,iorb0
    integer :: k1,k2,k3,a1b,a2b,a3b,ab1,ab2,ab3
    logical :: disp_sw
    logical, save :: flag_init_force = .true.
    include 'mpif.h'

    call write_border( 1, "calc_force_ps_nloc2_hp(start)" )

    call check_disp_switch( disp_sw, 0 )

    call watchb( ttmp ); tttt=0.0d0

    !flag_init_force = .true. ! MIZUHO-IR for cellopt (comment-out by ji for CPMD performance)
    if ( flag_init_force ) then
       call ps_nloc2_init_derivative
       flag_init_force = .false.
    end if

    call watchb( ttmp, tttt(:,1) )

    if ( Mlma <= 0 ) then
       force2(:,:)=0.0d0
       return
    end if

    call construct_gaunt_coef_L1( 3, yyy )

    pi2 = 2.0d0*acos(-1.0d0)

    maxerr=0.d0

    a1b=Igrid(1,1)
    a2b=Igrid(1,2)
    a3b=Igrid(1,3)
    ab1=Igrid(2,1)-Igrid(1,1)+1
    ab2=Igrid(2,2)-Igrid(1,2)+1
    ab3=Igrid(2,3)-Igrid(1,3)+1

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    c1=1.d0/ML1
    c2=1.d0/ML2
    c3=1.d0/ML3

    call watchb( ttmp, tttt(:,2) )

    if ( .not.allocated(ilm1) ) then
       L1=maxval(lo)+1
       n=maxval(norb)
       allocate( ilm1(0:L1,n,Nelement) ) ; ilm1=0
       do ik=1,Nelement
          lm1=0
          do iorb=1,norb(ik)
             L=lo(iorb,ik)
             do L1=abs(L-1),L+1
                lm1=lm1+1
                ilm1(L1,iorb,ik)=lm1
             end do
          end do
       end do
    end if

    if ( .not.allocated(y2b) .and. all(ippform /= 4) ) then
       lm1=maxval(ilm1)
       NRc=maxval(NRps)
       allocate( y2b(NRc,lm1,Nelement) )
       y2b=0.d0
       do ik=1,Nelement
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          do L1=abs(L-1),L+1
             lm1=ilm1(L1,iorb,ik)
             d1=0.d0
             d2=0.d0
             call spline(rad1(1,ik),dviod(1,lm1,ik),NRps(iorb,ik),d1,d2,y2b(1,lm1,ik))
          end do
       end do
       end do
    end if

    call watchb( ttmp, tttt(:,3) )

    allocate( wtmp5(0:3,nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1) )
    allocate( vtmp2(0:3,nzlma,MB_d) )
    allocate( a_rank(Natom) )
    allocate( duVdR(3,MMJJ,nzlma) )

    call watchb( ttmp, tttt(:,4) )

!$OMP parallel

!$omp master
    call watchb( ttmp )
!$omp end master

!$OMP workshare
    wtmp5=zero
    a_rank(:)=.false.
!$OMP end workshare

!$OMP do private( i1,i2,i3,k1,k2,k3 )
    do a=1,Natom
       i1 = nint( aa_atom(1,a)*ML1 )
       i2 = nint( aa_atom(2,a)*ML2 )
       i3 = nint( aa_atom(3,a)*ML3 )
       k1 = i1/ML1 ; if ( i1<0 ) k1=(i1+1)/ML1-1
       k2 = i2/ML2 ; if ( i2<0 ) k2=(i2+1)/ML2-1
       k3 = i3/ML3 ; if ( i3<0 ) k3=(i3+1)/ML3-1
       i1 = i1 - k1*ML1
       i2 = i2 - k2*ML2
       i3 = i3 - k3*ML3
       if ( Igrid(1,1) <= i1 .and. i1 <= Igrid(2,1) .and. &
            Igrid(1,2) <= i2 .and. i2 <= Igrid(2,2) .and. &
            Igrid(1,3) <= i3 .and. i3 <= Igrid(2,3) ) then
          a_rank(a)=.true.
       end if
    end do
!$OMP end do

!$omp master
    call watchb( ttmp, tttt(:,5) )
!$omp end master

    if ( any( ippform == 4 ) ) then

!$OMP master
       call watchb( ttmp )
!$OMP end master

!$OMP single
       call init_force_ps_nloc_gth &
            (MMJJ,nzlma,amap,lmap,mmap,iorbmap,MJJ_MAP,JJ_MAP,Y1,Y2,Y3,duVdR)
!$OMP end single

!$OMP master
       call watchb( ttmp, tttt(:,6) )
!$OMP end master

    else

!$omp master
    call watchb( ttmp )
!$omp end master

#ifndef _SPLINE_
!$OMP single
    allocate( irad(0:3000,Nelement) )
    irad=0
    M_irad=0
    do ik=1,Nelement
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
!$OMP end single
#endif

!$OMP workshare
    duVdR=0.d0
!$OMP end workshare

!$OMP master
    call watchb( ttmp, tttt(:,7) )
!$OMP end master

!$OMP do schedule(dynamic) firstprivate( maxerr ) &
!$OMP    private( a,L,m,iorb,ik,Rx,Ry,Rz,NRc,d1,d2,d3,x,y,z,r  &
!$OMP            ,ir,ir0,yy1,yy2,yy3,err0,err,tmp0,tmp1,m1,m2  &
!$OMP            ,lma,j,L1,L1z,lm1,im,ir1 )
    do lma=1,nzlma
       a    = amap(lma)
       if ( a <= 0 ) cycle
       L    = lmap(lma)
       m    = mmap(lma)
       iorb = iorbmap(lma)
       ik   = ki_atom(a)
       Rx=aa(1,1)*aa_atom(1,a)+aa(1,2)*aa_atom(2,a)+aa(1,3)*aa_atom(3,a)
       Ry=aa(2,1)*aa_atom(1,a)+aa(2,2)*aa_atom(2,a)+aa(2,3)*aa_atom(3,a)
       Rz=aa(3,1)*aa_atom(1,a)+aa(3,2)*aa_atom(2,a)+aa(3,3)*aa_atom(3,a)
       NRc=NRps(iorb,ik)
!!$OMP parallel do firstprivate( maxerr ) &
!!$OMP             private( d1,d2,d3,x,y,z,r,ir0,yy1,yy2,yy3,lm1,err,err0 &
!!$OMP                     ,tmp0,tmp1,m1,m2,j,L1,im,L1z )
       do j=1,MJJ_MAP(lma)
          d1=c1*JJ_MAP(1,j,lma)+JJ_MAP(4,j,lma)
          d2=c2*JJ_MAP(2,j,lma)+JJ_MAP(5,j,lma)
          d3=c3*JJ_MAP(3,j,lma)+JJ_MAP(6,j,lma)
          x = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
          y = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
          z = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
          r = sqrt(x*x+y*y+z*z)
!#ifndef _SPLINE_
!          ir0=irad( int(100.d0*r),ik )
!          do ir=ir0,NRc
!             if ( r<rad1(ir,ik) ) exit
!          end do
!#endif
          yy1=0.d0
          yy2=0.d0
          yy3=0.d0
          do L1=abs(L-1),L+1
             lm1=ilm1(L1,iorb,ik)
             if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 &
             .or. abs(z)>1.d-14 .or. L1==0 ) then
#ifdef _SPLINE_
                if ( r < rad1(2,ik) ) then
                   tmp0=dviod(2,lm1,ik)/(rad1(2,ik)**2)
                else
                   call splint(rad1(1,ik),dviod(1,lm1,ik),y2b(1,lm1,ik),NRc,r,tmp0)
                   tmp0=tmp0/(r*r)
                end if
#else
                ir0 = 1
                ir1 = NRc
111             if ( ir1 - ir0 > 1 ) then
                   ir = ( ir0 + ir1 )/2
                   if ( rad1(ir,ik) > r ) then
                      ir1 = ir
                   else
                      ir0 = ir
                   end if
                   goto 111
                end if
                ir = ir0
                if ( ir <= 2 ) then
                   err0=0.d0
                   tmp0=dviod(2,lm1,ik)/(rad1(2,ik)**2)
                   if ( ir < 1 ) stop "calc_force_ps_nloc2"
                else if ( ir <= NRc ) then
                   err0=1.d10
                   do im=1,20
                      m1=max(1,ir-im)
                      m2=min(ir+im,NRc)
                      call polint(rad1(m1,ik),dviod(m1,lm1,ik) &
                           ,m2-m1+1,r,tmp1,err)
                      if ( abs(err)<err0 ) then
                         tmp0=tmp1
                         err0=abs(err)
                         if ( err0<ep ) exit
                      end if
                   end do
                   tmp0=tmp0/(r*r)
                else
!                   write(*,*) "force_ps_nloc2",ir,NRc
!                   write(*,*) r,int(100.d0*r),irad( int(100.d0*r),ik )
!                   write(*,*) rad1(NRc,ik)
!                   stop
                   tmp0=0.0d0
                end if
                maxerr=max(maxerr,err0)
#endif
                do L1z=-L1,L1
                   tmp1=tmp0*Ylm(x,y,z,L1,L1z)
                   yy1=yy1+tmp1*yyy%x(L,m,L1,L1z)
                   yy2=yy2+tmp1*yyy%y(L,m,L1,L1z)
                   yy3=yy3-tmp1*yyy%z(L,m,L1,L1z)
                end do
             end if
          end do ! L1

          duVdR(1,j,lma)=yy1
          duVdR(2,j,lma)=yy2
          duVdR(3,j,lma)=yy3

       end do ! j
!!$OMP end parallel do
    end do ! lma
!$OMP end do

!$OMP master
    call watchb( ttmp, tttt(:,8) )
!$OMP end master

#ifndef _SPLINE_
!$OMP single
    deallocate( irad )
!$OMP end single
#endif

    end if

!$omp master
    call watchb( ttmp )
!$omp end master

    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
!$OMP do schedule(dynamic) private( c,i,d1,d2,d3,kr,ztmp,i1,i2,i3 )
    do n=MB_0,MB_1
       if ( occ(n,k,s) == 0.d0 ) cycle
       c=-2.d0*occ(n,k,s)*dV*dV
       do lma=1,nzlma
!          if ( MJJ_MAP(lma) == MJJ(lma) ) then
!#ifdef _DRSDFT_
!          do j=1,MJJ(lma)
!             i=JJP(j,lma)
!             wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s) &
!                  +uVk(j,lma,k)*unk(i,n,k,s)
!             wtmp5(1,lma,n,k,s)=wtmp5(1,lma,n,k,s)+duVdR(1,j,lma)*unk(i,n,k,s)
!             wtmp5(2,lma,n,k,s)=wtmp5(2,lma,n,k,s)+duVdR(2,j,lma)*unk(i,n,k,s)
!             wtmp5(3,lma,n,k,s)=wtmp5(3,lma,n,k,s)+duVdR(3,j,lma)*unk(i,n,k,s)
!          end do
!          wtmp5(0,lma,n,k,s)=iuV(lma)*c*wtmp5(0,lma,n,k,s)
!#else
!          do j=1,MJJ(lma)
!             i=JJP(j,lma)
!             wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s) &
!                  +uVk(j,lma,k)*conjg(unk(i,n,k,s))
!          end do
!          wtmp5(0,lma,n,k,s)=iuV(lma)*c*wtmp5(0,lma,n,k,s)
!          do j=1,MJJ_MAP(lma)
!             i=JJP(j,lma)
!             d1=c1*JJ_MAP(1,j,lma)+JJ_MAP(4,j,lma)
!             d2=c2*JJ_MAP(2,j,lma)+JJ_MAP(5,j,lma)
!             d3=c3*JJ_MAP(3,j,lma)+JJ_MAP(6,j,lma)
!             kr=pi2*(kbb(1,k)*d1+kbb(2,k)*d2+kbb(3,k)*d3)
!             ztmp=dcmplx(cos(kr),sin(kr))*unk(i,n,k,s)
!             wtmp5(1,lma,n,k,s)=wtmp5(1,lma,n,k,s)+duVdR(1,j,lma)*ztmp
!             wtmp5(2,lma,n,k,s)=wtmp5(2,lma,n,k,s)+duVdR(2,j,lma)*ztmp
!             wtmp5(3,lma,n,k,s)=wtmp5(3,lma,n,k,s)+duVdR(3,j,lma)*ztmp
!          end do
!#endif
!          else ! --- MJJ(lma) /= MJJ_MAP(lma) ---
#ifdef _DRSDFT_
          do j=1,MJJ(lma)
             i=JJP(j,lma)
             wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s) &
                  +uVk(j,lma,k)*unk(i,n,k,s)
          end do
          wtmp5(0,lma,n,k,s)=iuV(lma)*c*wtmp5(0,lma,n,k,s)
          do j=1,MJJ_MAP(lma)
             i1=JJ_MAP(1,j,lma)
             i2=JJ_MAP(2,j,lma)
             i3=JJ_MAP(3,j,lma)
             i = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ML_0
             wtmp5(1,lma,n,k,s)=wtmp5(1,lma,n,k,s)+duVdR(1,j,lma)*unk(i,n,k,s)
             wtmp5(2,lma,n,k,s)=wtmp5(2,lma,n,k,s)+duVdR(2,j,lma)*unk(i,n,k,s)
             wtmp5(3,lma,n,k,s)=wtmp5(3,lma,n,k,s)+duVdR(3,j,lma)*unk(i,n,k,s)
          end do
#else
          do j=1,MJJ(lma)
             i=JJP(j,lma)
             wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s) &
                  +uVk(j,lma,k)*conjg(unk(i,n,k,s))
          end do
          wtmp5(0,lma,n,k,s)=iuV(lma)*c*wtmp5(0,lma,n,k,s)
          do j=1,MJJ_MAP(lma)
             i1=JJ_MAP(1,j,lma)
             i2=JJ_MAP(2,j,lma)
             i3=JJ_MAP(3,j,lma)
             i = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ML_0
             d1=c1*i1+JJ_MAP(4,j,lma)
             d2=c2*i2+JJ_MAP(5,j,lma)
             d3=c3*i3+JJ_MAP(6,j,lma)
             kr=pi2*(kbb(1,k)*d1+kbb(2,k)*d2+kbb(3,k)*d3)
             ztmp=dcmplx(cos(kr),sin(kr))*unk(i,n,k,s)
             wtmp5(1,lma,n,k,s)=wtmp5(1,lma,n,k,s)+duVdR(1,j,lma)*ztmp
             wtmp5(2,lma,n,k,s)=wtmp5(2,lma,n,k,s)+duVdR(2,j,lma)*ztmp
             wtmp5(3,lma,n,k,s)=wtmp5(3,lma,n,k,s)+duVdR(3,j,lma)*ztmp
          end do
#endif
!          end if
       end do ! lma
    end do ! n
!$OMP end do
    end do ! k
    end do ! s

!$OMP master
    call watchb( ttmp, tttt(:,9) )
!$OMP end master

!$OMP single
    deallocate( duVdR )
    max_nreq=2*maxval( nrlma_xyz )
    allocate( ireq(max_nreq) )
    allocate( istatus(MPI_STATUS_SIZE,max_nreq) )
!$OMP end single

!$OMP workshare
    force2(:,:)=0.d0
!$OMP end workshare

!$omp master
    call watchb( ttmp, tttt(:,10) )
!$omp end master

!$OMP single
    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=MB_0,MB_1,MB_d

       ib1=n
       ib2=min(ib1+MB_d-1,MB_1)
       nnn=ib2-ib1+1

       if ( occ(n,k,s) == 0.d0 ) cycle

!       call do3StepComm_F(nrlma_xyz,num_2_rank,sendmap,recvmap,lma_nsend &
!                        ,sbufnl,rbufnl,nzlma,ib1,ib2,wtmp5(0,1,ib1,k,s))

       do i=1,6
          select case(i)
          case(1,3,5)
             j=i+1
             vtmp2(:,:,1:nnn)=wtmp5(:,:,ib1:ib2,k,s)
          case(2,4,6)
             j=i-1
          end select
          do m=1,nrlma_xyz(i)
             nreq=0
             irank=num_2_rank(m,i)
             jrank=num_2_rank(m,j)
             if( irank >= 0 )then
                i1=0
                do ib=1,nnn
                do i2=1,lma_nsend(irank)
                do i3=0,3
                   i1=i1+1
                   sbufnl(i1,irank)=vtmp2(i3,sendmap(i2,irank),ib)
                end do
                end do
                end do
                nreq=nreq+1
                call mpi_isend(sbufnl(1,irank),4*lma_nsend(irank)*nnn &
                     ,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
             end if
             if( jrank >= 0 )then
                nreq=nreq+1
                call mpi_irecv(rbufnl(1,jrank),4*lma_nsend(jrank)*nnn &
                     ,TYPE_MAIN,jrank,1,comm_grid,ireq(nreq),ierr)
             end if
             call mpi_waitall(nreq,ireq,istatus,ierr)
             if( jrank >= 0 )then
                i1=0
                do ib=ib1,ib2
                do i2=1,lma_nsend(jrank)
                do i3=0,3
                   i1=i1+1
                   wtmp5(i3,recvmap(i2,jrank),ib,k,s) &
                        =wtmp5(i3,recvmap(i2,jrank),ib,k,s)+rbufnl(i1,jrank)
                end do
                end do
                end do
             end if
          end do ! m
       end do ! i

       do ib=ib1,ib2
       do lma=1,nzlma
          a=amap(lma)
          if ( a <= 0 ) cycle
          if ( a_rank(a) ) then
             force2(1,a)=force2(1,a) &
                  +real(wtmp5(0,lma,ib,k,s)*wtmp5(1,lma,ib,k,s),8)
             force2(2,a)=force2(2,a) &
                  +real(wtmp5(0,lma,ib,k,s)*wtmp5(2,lma,ib,k,s),8)
             force2(3,a)=force2(3,a) &
                  +real(wtmp5(0,lma,ib,k,s)*wtmp5(3,lma,ib,k,s),8)
          end if
       end do
       end do

    end do ! n
    end do ! k
    end do ! s

    deallocate( istatus )
    deallocate( ireq )

    allocate( work2(3,Natom) )
!$OMP end single

!$omp master
    call watchb( ttmp, tttt(:,11) )
!$omp end master

!$OMP workshare
    work2(:,:)=force2(:,:)
!$OMP end workshare

!$OMP single
    call mpi_allreduce(work2,force2,3*Natom &
         ,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    deallocate( work2 )
    deallocate( a_rank )
    deallocate( vtmp2 )
    deallocate( wtmp5 )
!$OMP end single

!$OMP master
    call watchb( ttmp, tttt(:,12) )
!$OMP end master

!$OMP end parallel

!    if ( myrank == 0 ) then
!       do i=1,12
!          write(*,'("time_force_nloc2(",i2,")",2f10.5)') i, tttt(:,i)
!       end do
!    end if

    call write_border( 1, "calc_force_ps_nloc2_hp(end)" )

  end subroutine calc_force_ps_nloc2_hp

end module force_ps_nloc2_module
