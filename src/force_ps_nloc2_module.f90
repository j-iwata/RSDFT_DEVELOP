module force_ps_nloc2_module

  use force_sub_sub_module, only: gaunt, construct_gaunt_coef_l1

  implicit none

  private
  public :: calc_force_ps_nloc2

  type(gaunt) :: yyy

  real(8),allocatable :: y2b(:,:,:)

contains

  subroutine calc_force_ps_nloc2( force2 )
    use parallel_module, only: comm_grid, MB_d_nl
    use rgrid_variables, only: dV, Igrid, Ngrid
    use ps_nloc2_variables, only: rbufnl, sbufnl, TYPE_MAIN, Mlma, nzlma, MMJJ, &
    amap,lmap,mmap,iorbmap,MJJ_MAP,JJ_MAP,JJP,nrlma_xyz,MJJ,lma_nsend, &
    num_2_rank,iuV,sendmap,recvmap,d_uVk,z_uVk,flag_backup_uVunk_ps_nloc2, &
    d_backup_uVunk_ps_nloc2,z_backup_uVunk_ps_nloc2
    use ps_nloc2_init_module, only: ps_nloc2_init_derivative
    use atom_module, only: aa_atom, ki_atom
    use pseudopot_module, only: lo, norb, ippform, NRps
    use var_ps_member, only: rad1, dviod
    use array_bound_module, only: MB_0,MB_1,MBZ_0,MBZ_1,MSP_0,MSP_1,ML_0
    use ylm_module, only: Ylm
    use aa_module, only: aa
    use wf_module, only: occ,unk
    use var_sys_parameter, only: use_real8_wf
    use bz_module, only: kbb
    use ps_nloc2_comm_module, only: d_ps_nloc2_comm, z_ps_nloc2_comm
    use rsdft_allreduce_module, only: rsdft_allreduce
    use spline_module, only: splint, spline
    use ps_nloc_gth_module, only: init_force_ps_nloc_gth
    implicit none
    real(8),intent(inout) :: force2(:,:)
    include 'mpif.h'
    integer :: i1,i2,i3
    integer :: i,j,k,s,n,ir,iorb,L,L1,L1z,NRc,irank,jrank
    integer :: nreq,max_nreq,ib,ib1,ib2,nnn
    integer :: a,ik,m,lm1,lma,im,m1,m2
    integer :: ierr,ir0,ir1
#ifndef _SPLINE_
    integer :: M_irad
    integer,allocatable :: irad(:,:)
#endif
    integer,allocatable :: ireq(:),istatus(:,:),ilm1(:,:,:)
    real(8),parameter :: ep=1.d-8
    real(8),save :: Y1(0:3,-3:3,0:4,-4:4)
    real(8),save :: Y2(0:3,-3:3,0:4,-4:4)
    real(8),save :: Y3(0:3,-3:3,0:4,-4:4)
    real(8) :: err,err0,maxerr,Rx,Ry,Rz
    real(8) :: c1,c2,c3,d1,d2,d3
    real(8) :: x,y,z,r,kr,pi2,c
    real(8) :: tmp0,tmp1
    real(8) :: yy1,yy2,yy3
    real(8),allocatable :: force1(:,:),duVdR(:,:,:)
    real(8),allocatable :: wd(:,:,:,:,:), wd0(:,:,:,:)
    ! real(8),allocatable :: vd(:,:,:)
    complex(8) :: ztmp
    complex(8),allocatable :: wz(:,:,:,:,:), wz0(:,:,:,:)
    ! complex(8),allocatable :: vz(:,:,:)
    logical,save :: flag_Y = .true.
    logical,allocatable :: a_rank(:)
    integer :: ML1,ML2,ML3,Natom,Nelement,MI
    integer :: k1,k2,k3,a1b,a2b,a3b,ab1,ab2,ab3,ofs
    logical :: disp_sw
    logical, save :: flag_init_force = .true.

    call check_disp_switch( disp_sw, 0 )

    !call watchb( ttmp ); tttt=0.0d0

    if ( Mlma <= 0 ) return

    if ( flag_init_force ) then
      call ps_nloc2_init_derivative
      flag_init_force = .false.
    end if

    !call watchb( ttmp, tttt(:,1) )

    call construct_gaunt_coef_L1( 3, yyy )

    Natom = size(aa_atom,2)
    Nelement = maxval(ki_atom)

    pi2 = 2.0d0*acos(-1.0d0)

    maxerr=0.0d0

    ofs = Igrid(1,0)
    a1b = Igrid(1,1)
    a2b = Igrid(1,2)
    a3b = Igrid(1,3)
    ab1 = Igrid(2,1)-Igrid(1,1)+1
    ab2 = Igrid(2,2)-Igrid(1,2)+1
    ab3 = Igrid(2,3)-Igrid(1,3)+1

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    c1 = 1.0d0/ML1
    c2 = 1.0d0/ML2
    c3 = 1.0d0/ML3

    !call watchb( ttmp, tttt(:,2) )

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
      y2b=0.0d0
      do ik=1,Nelement
      do iorb=1,norb(ik)
        L=lo(iorb,ik)
        do L1=abs(L-1),L+1
          lm1=ilm1(L1,iorb,ik)
          d1=0.0d0
          d2=0.0d0
          call spline(rad1(1,ik),dviod(1,lm1,ik),NRps(iorb,ik),d1,d2,y2b(1,lm1,ik))
          end do
      end do
      end do
    end if

    !call watchb( ttmp, tttt(:,3) )

    allocate( a_rank(Natom) ); a_rank=.false.
    allocate( duVdR(3,MMJJ,nzlma) ); duVdR=0.0d0
    allocate( force1(3,Natom) ); force1=0.0d0

    !call watchb( ttmp, tttt(:,4) )

    !$omp parallel

    ! !$omp workshare
    ! a_rank(:)=.false.
    ! !$omp end workshare

    !$omp do private( i1,i2,i3,k1,k2,k3 )
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
    end do !a
    !$omp end do

    !call watchb( ttmp, tttt(:,5), barrier='on' )

    if ( any( ippform == 4 ) ) then

      !$omp master
      !call watchb( ttmp, barrier='on' )
      !$omp end master

      !$omp single
      call init_force_ps_nloc_gth &
            (MMJJ,nzlma,amap,lmap,mmap,iorbmap,MJJ_MAP,JJ_MAP,Y1,Y2,Y3,duVdR)
      !$omp end single

      !$omp master
      !call watchb( ttmp, tttt(:,6), barrier='on' )
      !$omp end master

    else

      !call watchb( ttmp, barrier='on' )

#ifndef _SPLINE_
      !$omp single
      allocate( irad(0:3000,Nelement) )
      irad=0
      M_irad=0
      do ik=1,Nelement
        NRc=maxval( NRps(:,ik) )
        NRc=min( 3000, NRc )
        m=0
        irad(0,ik)=1
        do ir=1,NRc
          m=int(100.0d0*rad1(ir,ik))+1
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
      !$omp end single
#endif

      ! !$omp workshare
      ! duVdR=0.0d0
      ! !$omp end workshare

      !$omp master
      !call watchb( ttmp, tttt(:,7), barrier='on' )
      !$omp end master

      !$omp do schedule(dynamic) firstprivate( maxerr ) &
      !$omp    private( a,L,m,iorb,ik,Rx,Ry,Rz,NRc,d1,d2,d3,x,y,z,r  &
      !$omp            ,ir,ir0,yy1,yy2,yy3,err0,err,tmp0,tmp1,m1,m2  &
      !$omp            ,lma,j,L1,L1z,lm1,im,ir1 )
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
        !!$omp parallel do firstprivate( maxerr ) &
        !!$omp             private( d1,d2,d3,x,y,z,r,ir0,yy1,yy2,yy3,lm1,err,err0 &
        !!$omp                     ,tmp0,tmp1,m1,m2,j,L1,im,L1z )
        do j=1,MJJ_MAP(lma)
          d1=c1*JJ_MAP(1,j,lma)+JJ_MAP(4,j,lma)
          d2=c2*JJ_MAP(2,j,lma)+JJ_MAP(5,j,lma)
          d3=c3*JJ_MAP(3,j,lma)+JJ_MAP(6,j,lma)
          x = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
          y = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
          z = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
          r = sqrt(x*x+y*y+z*z)
          yy1=0.0d0
          yy2=0.0d0
          yy3=0.0d0
          do L1=abs(L-1),L+1
            lm1=ilm1(L1,iorb,ik)
            if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 .or. abs(z)>1.d-14 .or. L1==0 ) then
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
          ! 111 if ( ir1 - ir0 > 1 ) then
          !       ir = ( ir0 + ir1 )/2
          !       if ( rad1(ir,ik) > r ) then
          !         ir1 = ir
          !       else
          !         ir0 = ir
          !       end if
          !       goto 111
          !     end if
              do while ( ir1 - ir0 > 1 )
                ir = ( ir0 + ir1 )/2
                if ( rad1(ir,ik) > r ) then
                  ir1 = ir
                else
                  ir0 = ir
                end if
              end do
              ir = ir0
              if ( ir <= 2 ) then
                err0=0.0d0
                tmp0=dviod(2,lm1,ik)/(rad1(2,ik)**2)
                if ( ir < 1 ) stop "calc_force_ps_nloc2"
              else if ( ir <= NRc ) then
                err0=1.d10
                do im=1,20
                  m1=max(1,ir-im)
                  m2=min(ir+im,NRc)
                  call polint(rad1(m1,ik),dviod(m1,lm1,ik),m2-m1+1,r,tmp1,err)
                  if ( abs(err) < err0 ) then
                    tmp0=tmp1
                    err0=abs(err)
                    if ( err0 < ep ) exit
                  end if
                end do
                tmp0=tmp0/(r*r)
              else
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

          duVdR(1,j,lma) = yy1
          duVdR(2,j,lma) = yy2
          duVdR(3,j,lma) = yy3

        end do !j
        !!$omp end parallel do
      end do !lma
      !$omp end do

      !$omp master
      !call watchb( ttmp, tttt(:,8), barrier='on' )
      !$omp end master

#ifndef _SPLINE_
      !$omp single
      deallocate( irad )
      !$omp end single
#endif

    end if

    !call watchb( ttmp, barrier='on' )

    if ( use_real8_wf() ) then

      !$omp single
      allocate( wd(0:3,nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1) )
      wd=0.0d0
      !$omp end single

      do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
      !$omp do schedule(dynamic) private( c,i,d1,d2,d3,kr,ztmp,i1,i2,i3 )
      do n=MB_0,MB_1
        if ( occ(n,k,s) == 0.0d0 ) cycle
        c=-2.0d0*occ(n,k,s)*dV*dV
        do lma=1,nzlma
          do j=1,MJJ(lma)
            i=JJP(j,lma)
            wd(0,lma,n,k,s) = wd(0,lma,n,k,s) + d_uVk(j,lma,k)*unk(i,n,k,s)
          end do
          wd(0,lma,n,k,s) = iuV(lma)*c*wd(0,lma,n,k,s)
          do j=1,MJJ_MAP(lma)
            i1=JJ_MAP(1,j,lma)
            i2=JJ_MAP(2,j,lma)
            i3=JJ_MAP(3,j,lma)
            i = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ofs
            wd(1,lma,n,k,s) = wd(1,lma,n,k,s) + duVdR(1,j,lma)*unk(i,n,k,s)
            wd(2,lma,n,k,s) = wd(2,lma,n,k,s) + duVdR(2,j,lma)*unk(i,n,k,s)
            wd(3,lma,n,k,s) = wd(3,lma,n,k,s) + duVdR(3,j,lma)*unk(i,n,k,s)
          end do
        end do !lma
      end do !n
      !$omp end do
      end do !k
      end do !s

    else ![ use complex16 wf ]

      !$omp single
      allocate( wz(0:3,nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1) )
      wz=(0.0d0,0.0d0)
      !$omp end single

      do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
      !$omp do schedule(dynamic) private( c,i,d1,d2,d3,kr,ztmp,i1,i2,i3 )
      do n=MB_0,MB_1
        if ( occ(n,k,s) == 0.0d0 ) cycle
        c=-2.0d0*occ(n,k,s)*dV*dV
        do lma=1,nzlma
          do j=1,MJJ(lma)
            i=JJP(j,lma)
            ! wz(0,lma,n,k,s) = wz(0,lma,n,k,s) + z_uVk(j,lma,k)*conjg(unk(i,n,k,s))
            ztmp=unk(i,n,k,s)
            wz(0,lma,n,k,s) = wz(0,lma,n,k,s) + z_uVk(j,lma,k)*conjg(ztmp)
          end do
          wz(0,lma,n,k,s) = iuV(lma)*c*wz(0,lma,n,k,s)
          do j=1,MJJ_MAP(lma)
            i1=JJ_MAP(1,j,lma)
            i2=JJ_MAP(2,j,lma)
            i3=JJ_MAP(3,j,lma)
            i = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ofs
            d1=c1*i1+JJ_MAP(4,j,lma)
            d2=c2*i2+JJ_MAP(5,j,lma)
            d3=c3*i3+JJ_MAP(6,j,lma)
            kr=pi2*(kbb(1,k)*d1+kbb(2,k)*d2+kbb(3,k)*d3)
            ztmp=dcmplx(cos(kr),sin(kr))*unk(i,n,k,s)
            wz(1,lma,n,k,s) = wz(1,lma,n,k,s) + duVdR(1,j,lma)*ztmp
            wz(2,lma,n,k,s) = wz(2,lma,n,k,s) + duVdR(2,j,lma)*ztmp
            wz(3,lma,n,k,s) = wz(3,lma,n,k,s) + duVdR(3,j,lma)*ztmp
          end do
        end do ! lma
      end do ! n
      !$omp end do
      end do ! k
      end do ! s
    
    end if

    !$omp master
    !call watchb( ttmp, tttt(:,9), barrier='on' )
    !$omp end master

    !$omp single
    deallocate( duVdR )
    ! max_nreq=2*maxval( nrlma_xyz )
    ! allocate( ireq(max_nreq) )
    ! allocate( istatus(MPI_STATUS_SIZE,max_nreq) )
    !$omp end single

    !$omp master
    !call watchb( ttmp, tttt(:,10), barrier='on' )
    !$omp end master

    !$omp single
    !call watchb( ttmp, barrier='on' )

    ! if ( use_real8_wf() ) then
    !   allocate( vd(0:3,nzlma,MB_d_nl) ); vd=0.0d0
    ! else
    !   allocate( vz(0:3,nzlma,MB_d_nl) ); vz=(0.0d0,0.0d0)
    ! end if

    do s = MSP_0, MSP_1
    do k = MBZ_0, MBZ_1
    do n = MB_0 , MB_1 , MB_d_nl

      if ( abs(occ(n,k,s)) < 1.0d-8 ) cycle

      ib1 = n
      ib2 = min( ib1 + MB_d_nl - 1, MB_1 )

      if ( use_real8_wf() ) then
        call d_ps_nloc2_comm( wd(:,:,ib1:ib2,k,s) )
      else
        call z_ps_nloc2_comm( wz(:,:,ib1:ib2,k,s) )
      end if

      ! do i=1,6
      !   select case(i)
      !   case(1,3,5)
      !     j=i+1
      !     vd(:,:,1:nnn)=wd(:,:,ib1:ib2,k,s)
      !   case(2,4,6)
      !     j=i-1
      !   end select
      !   do m=1,nrlma_xyz(i)
      !     nreq=0
      !     irank=num_2_rank(m,i)
      !     jrank=num_2_rank(m,j)
      !     if( irank >= 0 )then
      !       i1=0
      !       do ib=1,nnn
      !       do i2=1,lma_nsend(irank)
      !       do i3=0,3
      !         i1=i1+1
      !         sbufnl(i1,irank)=vtmp2(i3,sendmap(i2,irank),ib)
      !       end do
      !       end do
      !       end do
      !       nreq=nreq+1
      !       call mpi_isend(sbufnl(1,irank),4*lma_nsend(irank)*nnn &
      !           ,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
      !     end if
      !     if( jrank >= 0 )then
      !       nreq=nreq+1
      !       call mpi_irecv(rbufnl(1,jrank),4*lma_nsend(jrank)*nnn &
      !           ,TYPE_MAIN,jrank,1,comm_grid,ireq(nreq),ierr)
      !     end if
      !     call mpi_waitall(nreq,ireq,istatus,ierr)
      !     if( jrank >= 0 )then
      !       i1=0
      !       do ib=ib1,ib2
      !       do i2=1,lma_nsend(jrank)
      !       do i3=0,3
      !         i1=i1+1
      !         wtmp5(i3,recvmap(i2,jrank),ib,k,s) &
      !               =wtmp5(i3,recvmap(i2,jrank),ib,k,s)+rbufnl(i1,jrank)
      !       end do
      !       end do
      !       end do
      !     end if
      !   end do ! m
      ! end do ! i

      if ( use_real8_wf() ) then
        do ib = ib1, ib2
        do lma = 1, nzlma
          a = amap(lma); if ( a <= 0 ) cycle
          if ( a_rank(a) ) then
            force1(1,a) = force1(1,a) + real( wd(0,lma,ib,k,s)*wd(1,lma,ib,k,s) )
            force1(2,a) = force1(2,a) + real( wd(0,lma,ib,k,s)*wd(2,lma,ib,k,s) )
            force1(3,a) = force1(3,a) + real( wd(0,lma,ib,k,s)*wd(3,lma,ib,k,s) )
          end if
        end do
        end do
      else ![ use complex16 wf ]
        do ib = ib1, ib2
        do lma = 1, nzlma
          a = amap(lma); if ( a <= 0 ) cycle
          if ( a_rank(a) ) then
            force1(1,a) = force1(1,a) + real( wz(0,lma,ib,k,s)*wz(1,lma,ib,k,s) )
            force1(2,a) = force1(2,a) + real( wz(0,lma,ib,k,s)*wz(2,lma,ib,k,s) )
            force1(3,a) = force1(3,a) + real( wz(0,lma,ib,k,s)*wz(3,lma,ib,k,s) )
          end if
        end do
        end do
      end if

    end do ! n
    end do ! k
    end do ! s

    ! deallocate( istatus )
    ! deallocate( ireq )

    !call watchb( ttmp, tttt(:,11) )
    !$omp end single

    ! ---

    if ( flag_backup_uVunk_ps_nloc2 ) then
      if ( use_real8_wf() ) then
        !$omp single
        allocate( wd0(nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1) )
        wd0=0.0d0
        !$omp end single
        do s = MSP_0, MSP_1
        do k = MBZ_0, MBZ_1
        !$omp do schedule(dynamic) private( n, c )
        do n = MB_0 , MB_1
          if ( occ(n,k,s) == 0.0d0 ) cycle
          c = 1.0d0/( -2.0d0*occ(n,k,s)*dV )
          wd0(:,n,k,s) = c*wd(0,:,n,k,s)
        end do !n
        end do !k
        end do !s
        call d_backup_uVunk_ps_nloc2( wd0 )
        !$omp single
        deallocate( wd0 )
        !$omp end single
      else
        !$omp single
        allocate( wz0(nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1) )
        wz0=(0.0d0,0.0d0)
        !$omp end single
        do s = MSP_0, MSP_1
        do k = MBZ_0, MBZ_1
        !$omp do schedule(dynamic) private( n, c )
        do n = MB_0 , MB_1
          if ( occ(n,k,s) == 0.0d0 ) cycle
          c = 1.0d0/( -2.0d0*occ(n,k,s)*dV )
          wz0(:,n,k,s) = c*wz(0,:,n,k,s)
        end do !n
        end do !k
        end do !s
        call z_backup_uVunk_ps_nloc2( wz0 )
        !$omp single
        deallocate( wz0 )
        !$omp end single
      end if
    end if

    ! ---

    !$omp master
    !call watchb( ttmp )
    !$omp end master


    !$omp single
    if ( allocated(wd) ) deallocate( wd )
    if ( allocated(wz) ) deallocate( wz )
    deallocate( a_rank )
    call rsdft_allreduce( force1 )
    !$omp end single

    !$omp workshare
    ! force2 = force2 + force1
    force2 = force1
    !$omp end workshare

    !$omp master
    !call watchb( ttmp )
    !$omp end master

    !$omp end parallel

    !if ( myrank == 0 ) then
    !   do i=1,12
    !      write(*,'("time_force_nloc2(",i2")') tttt(:,i)
    !   end do
    !end if

  end subroutine calc_force_ps_nloc2

end module force_ps_nloc2_module
