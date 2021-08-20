module ps_nloc2_module

  use aa_module, only: aa
  use atom_module, only: aa_atom, ki_atom
  use array_bound_module, only: ML_0, MB_0, MB_1, MBZ_0, MBZ_1, MSP_0, MSP_1
  use pseudopot_module, only: lo, norb, NRps, ippform, ps_type, pselect
  use var_ps_member, only: viod, rad1, Rps, inorm
  use ps_nloc2_variables, only: amap, lmap, mmap, iorbmap, nzlma, nrlma_xyz &
                               ,num_2_rank, lma_nsend, sendmap, recvmap &
                               ,JJ_MAP, MJJ_MAP, iuV, MJJ, MMJJ, Mlma &
                               ,uV, iamap, nl_rank_map, nl_max_send &
                               ,FLAG_KEEP_uV, FLAG_KEEP_JJ_MAP &
                               ,d_allocate_ps_nloc2, z_allocate_ps_nloc2
  use ps_nloc_gth_module, only: init_ps_nloc_gth, prep_ps_nloc_gth
  use minimal_box_module, only: m_grid_ion, map_grid_ion, mcube_grid_ion, make_minimal_box
  use bz_module, only: kbb
  use watch_module, only: watchb
  use ps_nloc2_op_module, only: reset_init_flag_ps_nloc2_op
  use ylm_module, only: Ylm
  use hsort_module, only: indexx, sort_index_sub
  use spline_module, only: splint, spline
  use memory_module, only: check_memory

  implicit none

  private
  public :: prep_ps_nloc2
  public :: prep_uvk_ps_nloc2
  public :: prep_rvk_ps_nloc2

  real(8),allocatable :: y2a(:,:,:),y2b(:,:,:)
  ! integer,allocatable :: ilm1(:,:,:)

  integer,allocatable :: icheck_tmp3(:,:,:)
  integer,allocatable :: JJ_tmp(:,:,:,:)
  integer,allocatable :: MJJ_tmp(:,:)
  real(8),allocatable :: uV_tmp(:,:,:)

contains


  subroutine prep_ps_nloc2
    use parallel_module, only: MB_d_nl, myrank, node_partition, myrank_g, comm_grid &
                              ,np_grid
    use rgrid_variables, only: Ngrid, Igrid, Hgrid
    use var_sys_parameter, only: use_real8_wf
    use ps_nloc2_variables, only: d_uVk, z_uVk, JJP
    implicit none
    include 'mpif.h'
    integer,allocatable :: icheck_tmp1(:),icheck_tmp2(:),itmp(:,:)
    integer,allocatable :: icheck_tmp4(:,:,:)
    integer,allocatable :: sendmap_tmp(:,:),recvmap_tmp(:,:),ireq(:)
    integer,allocatable :: lma_nsend_tmp(:),maps_tmp(:,:),itmp1(:)
    integer,allocatable :: nl_rank_map_tmp(:),itmp3(:,:)
    integer,allocatable :: itmp2(:),LLp(:,:)
    integer,allocatable :: istatus(:,:)
    integer :: a,i,j,k,L,m,n,mm1,mm2,mm3,m1,m2,k1,k2,k3
    integer :: i1,i2,i3,j1,j2,j3,ik,ir,iorb,mm,ierr,ir0,ir1,irlma
    integer :: ic1,ic2,ic3,id1,id2,id3
    integer :: nzlma_0,NRc,MMJJ_0,lma,lma0,i1_0,i2_0,i3_0
    integer :: nreq, nprocs_g
    real(8),parameter :: ep=1.d-8
    real(8) :: x,y,z,r,Rx,Ry,Rz,Rps2,v,v0,d1,d2,d3,r2
    real(8) :: c1,c2,c3,a1,a2,a3,maxerr,err0,err
    real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
    real(8),allocatable :: work(:)
    !real(8) :: ttmp(2),tttt(2,11)
    integer :: ML1,ML2,ML3,a1b,b1b,a2b,b2b,a3b,b3b
    integer :: ab1,ab2,ab3
    integer :: np1,np2,np3,nrlma,Natom,Nelement
    logical,allocatable :: lcheck_tmp1(:,:), lcheck_tmp2(:,:)
    logical :: disp_sw
    integer :: num_local_atom, max_num_local_atom
    integer,allocatable :: lst_local_atom(:)
    integer :: ia

    !call watchb( ttmp, barrier='on' ); tttt=0.0d0

    call write_border( 0, " prep_ps_nloc2(start)" )
    call check_disp_switch( disp_sw, 0 )

!------------------------------------------ max atom*orb

    nprocs_g = np_grid
    Natom = size(ki_atom)
    Nelement = maxval(ki_atom)

    Mlma=0
    do i=1,Natom
      ik=ki_atom(i)
      do iorb=1,norb(ik)
        Mlma=Mlma+2*lo(iorb,ik)+1
      end do
    end do

    if ( Mlma <= 0 ) return   ! if no orb return
!--------------------------------------------------------

    if ( .not.allocated(y2a) .and. all(ippform /= 4) ) then
      NRc=maxval(NRps)
      n=maxval(norb)
      allocate( y2a(NRc,n,Nelement) )
      y2a=0.0d0
      do ik=1,Nelement
      do iorb=1,norb(ik)
        d1=0.0d0
        d2=0.0d0
        call spline(rad1(1,ik),viod(1,iorb,ik),NRps(iorb,ik) &
             ,d1,d2,y2a(1,iorb,ik))
      end do
      end do
    end if

    if ( all(ippform == 4) ) call init_ps_nloc_gth(disp_sw)

    if ( Mlma < nprocs_g ) then
      nzlma_0 = Mlma
    else
      nzlma_0 = min(Mlma*125/nprocs_g,Mlma)
    end if

!
!      __________
!     /         /|
!    /         / |
!   /         /  |
!  /b2b      /   |
!  -----------   |
!  |    _____|___|
!  |   /b3b  |   /
!  |  /      |  /
!  | /       | /
!  |/        |/
!  ----------
! a1b  ab1  b1b
! a2b
! a3b
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

    !call watchb( ttmp, tttt(:,1), barrier='on' )

!----- make_minimal_box -----
! atom position centered grid
! minimal box for largest atom
! MMJJ_0 : total # of grid points
    r=maxval(Rps)+maxval(Hgrid(1:3))+1.d-8
    call make_minimal_box(r,mm1,mm2,mm3,MMJJ_0)
    mm1 = maxval( abs(mcube_grid_ion(:,1)) ) + 1
    mm2 = maxval( abs(mcube_grid_ion(:,2)) ) + 1
    mm3 = maxval( abs(mcube_grid_ion(:,3)) ) + 1

! IN:   r
! OUT:  m_grid_ion==MMJJ_0
!       map_grid_ion(1:3,1:m_grid_ion)
!       mcube_grid_ion(1_min:2_max,1:3) - max & min value of map_grid_ion
!       mm1,mm2,mm3 : max of abs(mcube_grid_ion)

    !call watchb( ttmp, tttt(:,2), barrier='on' )

    MMJJ_0 = M_grid_ion

! ---

    c1 = 1.0d0/dble(ML1)
    c2 = 1.0d0/dble(ML2)
    c3 = 1.0d0/dble(ML3)

    xmin= 1.0d100
    xmax=-1.0d100
    ymin= 1.0d100
    ymax=-1.0d100
    zmin= 1.0d100
    zmax=-1.0d100
    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=a1b,b1b
      d1=i1*c1
      d2=i2*c2
      d3=i3*c3
      x=aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3
      y=aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3
      z=aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3
      xmin=min(x,xmin)
      xmax=max(x,xmax)
      ymin=min(y,ymin)
      ymax=max(y,ymax)
      zmin=min(z,zmin)
      zmax=max(z,zmax)
    end do
    end do
    end do

    xmin = xmin - r
    xmax = xmax + r
    ymin = ymin - r
    ymax = ymax + r
    zmin = zmin - r
    zmax = zmax + r

    allocate( itmp1(Natom) ); itmp1=0

    n = 0
    loop_a: do a = 1, Natom
      do i3 = -1, 1
      do i2 = -1, 1
      do i1 = -1, 1
        a1 = aa_atom(1,a) + i1
        a2 = aa_atom(2,a) + i2
        a3 = aa_atom(3,a) + i3
        Rx = aa(1,1)*a1 + aa(1,2)*a2 + aa(1,3)*a3
        Ry = aa(2,1)*a1 + aa(2,2)*a2 + aa(2,3)*a3
        Rz = aa(3,1)*a1 + aa(3,2)*a2 + aa(3,3)*a3
        if ( xmin <= Rx .and. Rx <= xmax .and. &
             ymin <= Ry .and. Ry <= ymax .and. &
             zmin <= Rz .and. Rz <= zmax ) then
          n = n + 1
          itmp1(n)=a
          cycle loop_a
        end if
      end do
      end do
      end do
    end do loop_a

    call MPI_Allreduce(n,m1,1,MPI_INTEGER,MPI_SUM,comm_grid,ierr)
    call MPI_Allreduce(n,m2,1,MPI_INTEGER,MPI_MAX,comm_grid,ierr)

    !write(*,'(i4,2x,3(2f12.5,2x),3i6)') myrank,xmin,xmax,ymin,ymax,zmin,zmax,n,m1,m2

    if ( any(ippform==4) ) then
      max_num_local_atom = Natom
      num_local_atom = Natom
    else
      max_num_local_atom = m2
      num_local_atom = n
     !max_num_local_atom = Natom
     !num_local_atom = Natom
    end if

    if ( allocated(lst_local_atom) ) deallocate(lst_local_atom)
    allocate( lst_local_atom(max_num_local_atom) ); lst_local_atom=0

    do i = 1, num_local_atom
      lst_local_atom(i) = itmp1(i)
    end do

    deallocate( itmp1 )

! ---

    !if ( allocated(icheck_tmp3) ) deallocate(icheck_tmp3)
    !if ( allocated(JJ_tmp)      ) deallocate(JJ_tmp)
    !if ( allocated(MJJ_tmp)     ) deallocate(MJJ_tmp)
    !if ( allocated(uV_tmp)      ) deallocate(uV_tmp)
    L=maxval(lo)
    n=maxval(norb)
    allocate( icheck_tmp3(max_num_local_atom,n,2*L+1) ); icheck_tmp3=0
    allocate( JJ_tmp(6,MMJJ_0,n,max_num_local_atom)   ); JJ_tmp=0
    allocate( MJJ_tmp(n,max_num_local_atom)           ); MJJ_tmp=0
    allocate( uV_tmp(MMJJ_0,n,max_num_local_atom)     ); uV_tmp=0.0d0

    !call watchb( ttmp, tttt(:,3), barrier='on' )

    if ( any( ippform == 4 ) ) then

      call prep_ps_nloc_gth(Natom,n,L,MMJJ_0,M_grid_ion,map_grid_ion &
                           ,icheck_tmp3,JJ_tmp,MJJ_tmp,uV_tmp,nzlma,MMJJ)

      if ( .not.all( ippform == 4 ) ) then
        write(*,*) "Mixed use of different pseudopotenial is forbidden"
        stop "stop@prep_ps_nloc2"
      end if

    else

#ifndef _SPLINE_
    allocate( irad(0:3000,Nelement) ) ; irad=0
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
#endif

    c1                 = 1.0d0/Ngrid(1)
    c2                 = 1.0d0/Ngrid(2)
    c3                 = 1.0d0/Ngrid(3)
    maxerr             = 0
    icheck_tmp3(:,:,:) = 0
    MMJJ               = 0
    nzlma              = 0
    lma                = 0
    lma0               = 0

!$OMP parallel do schedule(dynamic) firstprivate( maxerr ) &
!$OMP    private( Rx,Ry,Rz,ic1,ic2,ic3,ik,iorb,Rps2,NRc,L,j,i,i1,i2,i3 &
!$OMP            ,id1,id2,id3,k1,k2,k3,i1_0,i2_0,i3_0,d1,d2,d3,x,y,z,r2,r &
!$OMP            ,v0,err0,ir0,ir,mm,m1,m2,v,err,ir1,a )
    do ia=1,num_local_atom

      a = lst_local_atom(ia)

! Rx,Ry,Rz : atom position in real grid
      Rx = aa(1,1)*aa_atom(1,a)+aa(1,2)*aa_atom(2,a)+aa(1,3)*aa_atom(3,a)
      Ry = aa(2,1)*aa_atom(1,a)+aa(2,2)*aa_atom(2,a)+aa(2,3)*aa_atom(3,a)
      Rz = aa(3,1)*aa_atom(1,a)+aa(3,2)*aa_atom(2,a)+aa(3,3)*aa_atom(3,a)

! ic1,ic2,ic3 : atom position in grid point
      ic1 = nint( aa_atom(1,a)*Ngrid(1) )
      ic2 = nint( aa_atom(2,a)*Ngrid(2) )
      ic3 = nint( aa_atom(3,a)*Ngrid(3) )

      ik = ki_atom(a)

! process for each orbit
      do iorb=1,norb(ik)

        Rps2 = Rps(iorb,ik)**2
        NRc  = NRps(iorb,ik)
        L    = lo(iorb,ik)

!--------------------------------------------------for minimal box
! j : counting the # of grids
! matching minimal box and total grid
        j    = 0
        do i=1,M_grid_ion

! grid point in minimal box
          i1 = map_grid_ion(1,i)
          i2 = map_grid_ion(2,i)
          i3 = map_grid_ion(3,i)

! grid point in total grid
          id1 = ic1 + i1
          id2 = ic2 + i2
          id3 = ic3 + i3

! for parallel computation
! if i1,i2,i3 is negative : for P.B.C.
          k1=id1/ML1 ; if ( id1<0 ) k1=(id1+1)/ML1-1
          k2=id2/ML2 ; if ( id2<0 ) k2=(id2+1)/ML2-1
          k3=id3/ML3 ; if ( id3<0 ) k3=(id3+1)/ML3-1
          i1_0=id1-k1*ML1
          i2_0=id2-k2*ML2
          i3_0=id3-k3*ML3
! if this point is in this node

          if ( Igrid(1,1) <= i1_0 .and. i1_0 <= Igrid(2,1) .and. &
               Igrid(1,2) <= i2_0 .and. i2_0 <= Igrid(2,2) .and. &
               Igrid(1,3) <= i3_0 .and. i3_0 <= Igrid(2,3) ) then

! ratio adjustment
            d1 = id1*c1
            d2 = id2*c2
            d3 = id3*c3

! get the real grid of this point centered from atom position
            x  = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3 - Rx
            y  = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3 - Ry
            z  = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3 - Rz
            r2 = x*x+y*y+z*z

! if this point is out of PS-data then skip
            if ( r2 > Rps2+1.d-10 ) cycle

            r    = sqrt(r2)
            v0   = 0.0d0
            err0 = 0.0d0

! if this point is away from atom position or l==0
! interpolation is needed
            if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 .or. &
                 abs(z)>1.d-14 .or. L==0 ) then
#ifdef _SPLINE_
              call splint(rad1(1,ik),viod(1,iorb,ik),y2a(1,iorb,ik),NRc,r,v0)
#else
!                   ir0=irad( int(100.d0*r),ik )
!                   do ir=ir0,NRc
!                      if ( r < rad1(ir,ik) ) exit
!                   end do
              ir0 = 1
              ir1 = NRc
111           if ( ir1 - ir0 > 1 ) then
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
                v0=viod(2,iorb,ik)
                if ( ir < 1 ) stop "ps_nloc2(0)"
              else if ( ir <= NRc ) then
                err0=1.d10
                do mm=1,20
                  m1=max(1,ir-mm)
                  m2=min(ir+mm,NRc)
                  call polint &
                       (rad1(m1,ik),viod(m1,iorb,ik),m2-m1+1,r,v,err)
                  if ( abs(err) < err0 ) then
                    v0=v
                    err0=abs(err)
                    if ( err0 < ep ) exit
                  end if
                end do
              else
                write(*,*) "ps_nloc2(1)",ir,NRc,viod(NRc,iorb,ik)
                write(*,*) viod(NRc+1,iorb,ik),r,rad1(ir,ik)
                stop ' ERROR : abnormal end'
              end if
              maxerr=max(maxerr,err0)
#endif
            end if

! j : start from 1
! with i1_0,i2_0,i3_0 k1,k2,k3
!                if (abs(v0) < 1.d-10 ) cycle 
            j=j+1
            JJ_tmp(1,j,iorb,ia) = i1_0
            JJ_tmp(2,j,iorb,ia) = i2_0
            JJ_tmp(3,j,iorb,ia) = i3_0
            JJ_tmp(4,j,iorb,ia) = k1
            JJ_tmp(5,j,iorb,ia) = k2
            JJ_tmp(6,j,iorb,ia) = k3
            uV_tmp(j,iorb,ia)   = v0

          end if

        end do ! i ( 1 - M_grid_ion )

!==================================================for minimal box
! if there is not point in this node, j==0

        MJJ_tmp(iorb,ia)=j

      end do ! iorb
    end do ! ia
!$OMP end parallel do

#ifndef _SPLINE_
    deallocate( irad )
#endif

    end if ! ippform

! ---

    allocate( itmp1(Natom) ); itmp1=0
    do ia=1,num_local_atom
      itmp1(lst_local_atom(ia))=ia
    end do

    lma=0
    do a=1,Natom
      ia=itmp1(a)
      if ( ia == 0 ) cycle
      ik=ki_atom(a)
      do iorb=1,norb(ik)
        j=MJJ_tmp(iorb,ia)
        if ( j > 0 ) then
          L=lo(iorb,ik)
          nzlma=nzlma+2*L+1
          do m=1,2*L+1
            lma=lma+1
            icheck_tmp3(ia,iorb,m)=lma
          end do
        end if
      end do
    end do
    MMJJ = maxval( MJJ_tmp )

    allocate( lcheck_tmp1(0:np_grid-1,Mlma) ) ; lcheck_tmp1(:,:)=.false.
    allocate( lcheck_tmp2(0:np_grid-1,Mlma) ) ; lcheck_tmp2(:,:)=.false.
    lma=0
    do a=1,Natom
      ia=itmp1(a)
      ik=ki_atom(a)
      do iorb=1,norb(ik)
        L=lo(iorb,ik)
        j=0
        if ( ia /= 0 ) j=MJJ_tmp(iorb,ia)
        do m=1,2*L+1
          lma=lma+1
          if ( j > 0 ) then
            lcheck_tmp2(myrank_g,lma)=.true.
          end if
        end do
      end do
    end do
    call MPI_ALLREDUCE( lcheck_tmp2, lcheck_tmp1, size(lcheck_tmp1), &
                        MPI_LOGICAL, MPI_LOR, comm_grid, ierr )
    deallocate( lcheck_tmp2 )

    !call watchb( ttmp, tttt(:,4), barrier='on' )

! for grid-parallel computation

    nzlma_0 = min(nzlma_0*2,Mlma)

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

    nrlma=0
    lma=0
    do a=1,Natom
      ia=itmp1(a)
      ik=ki_atom(a)
    do iorb=1,norb(ik)
      L=lo(iorb,ik)
    do m=-L,L
      lma=lma+1

      icheck_tmp1(:)=0
      do n=0,np_grid-1
        if ( lcheck_tmp1(n,lma) ) icheck_tmp1(n) = 1
      end do
      if ( ia /= 0 ) icheck_tmp1(myrank_g) = icheck_tmp3(ia,iorb,m+L+1)

      itmp(:,:)=0
      n=-1
      do i3=1,np3
      do i2=1,np2
      do i1=1,np1
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

      if ( icheck_tmp1(myrank_g) /= 0 ) then
        if ( icheck_tmp1(myrank_g) > 0 ) then
          maps_tmp(icheck_tmp2(myrank_g),1)=icheck_tmp1(myrank_g)
        end if
        maps_tmp(icheck_tmp2(myrank_g),2)=inorm(iorb,ik)
        maps_tmp(icheck_tmp2(myrank_g),3)=a
        maps_tmp(icheck_tmp2(myrank_g),4)=L
        maps_tmp(icheck_tmp2(myrank_g),5)=m
        maps_tmp(icheck_tmp2(myrank_g),6)=iorb

        do n=0,nprocs_g-1
          if ( n == myrank_g .or. icheck_tmp1(n) == 0 ) cycle
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

    !call watchb( ttmp, tttt(:,5), barrier='on' )

    nzlma = icheck_tmp2(myrank_g)

    deallocate( itmp )
    deallocate( icheck_tmp2 )
    deallocate( icheck_tmp1 )
    deallocate( lcheck_tmp1 )
    deallocate( icheck_tmp3 )

    if ( allocated(uV)          ) deallocate( uV )
    if ( allocated(JJ_MAP)      ) deallocate( JJ_MAP )
    if ( allocated(MJJ_MAP)     ) deallocate( MJJ_MAP )
    if ( allocated(MJJ)         ) deallocate( MJJ )
    if ( allocated(iuV)         ) deallocate( iuV )
    if ( allocated(iamap)       ) deallocate( iamap )
    if ( allocated(amap)        ) deallocate( amap )
    if ( allocated(lmap)        ) deallocate( lmap )
    if ( allocated(mmap)        ) deallocate( mmap )
    if ( allocated(iorbmap)     ) deallocate( iorbmap )
    if ( allocated(nl_rank_map) ) deallocate( nl_rank_map )
    !if ( disp_sw ) write(*,*) "(uV)"
    !call check_memory( 8.0d0, MMJJ, nzlma )
    !if ( disp_sw ) write(*,*) "(JJ_MAP)(KEEP_FLAG)", FLAG_KEEP_JJ_MAP
    !call check_memory( 4.0d0, 6, MMJJ, nzlma )
    !if ( disp_sw ) write(*,*) "(MJJ_MAP)"
    !call check_memory( 4.0d0, nzlma )
    allocate( uV(MMJJ,nzlma)       ) ; uV=0.d0
    allocate( JJ_MAP(6,MMJJ,nzlma) ) ; JJ_MAP=0
    allocate( MJJ_MAP(nzlma)       ) ; MJJ_MAP=0
    allocate( MJJ(nzlma)           ) ; MJJ=0
    allocate( iuV(nzlma)           ) ; iuV=0
    allocate( iamap(nzlma)         ) ; iamap=0
    allocate( amap(nzlma)          ) ; amap=0
    allocate( lmap(nzlma)          ) ; lmap=0
    allocate( mmap(nzlma)          ) ; mmap=0
    allocate( iorbmap(nzlma)       ) ; iorbmap=0
    allocate( nl_rank_map(nrlma)   ) ; nl_rank_map=-1

    !call watchb( ttmp, tttt(:,6), barrier='on' )

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
      iamap(lma)   = itmp1(amap(lma))
    end do

    c1=1.0d0/ML1
    c2=1.0d0/ML2
    c3=1.0d0/ML3

!$omp parallel do private( ia,a,l,m,iorb,Rx,Ry,Rz,j,i1,i2,i3,k1,k2,k3,d1,d2,d3,x,y,z )
    do lma=1,nzlma
      if ( maps_tmp(lma,1) == 0 ) cycle
      ia   = iamap(lma)
      a    = amap(lma)
      l    = lmap(lma)
      m    = mmap(lma)
      iorb = iorbmap(lma)
      MJJ_MAP(lma) = MJJ_tmp(iorb,ia)
      Rx=aa(1,1)*aa_atom(1,a)+aa(1,2)*aa_atom(2,a)+aa(1,3)*aa_atom(3,a)
      Ry=aa(2,1)*aa_atom(1,a)+aa(2,2)*aa_atom(2,a)+aa(2,3)*aa_atom(3,a)
      Rz=aa(3,1)*aa_atom(1,a)+aa(3,2)*aa_atom(2,a)+aa(3,3)*aa_atom(3,a)
      do j=1,MJJ_MAP(lma)
        i1=JJ_tmp(1,j,iorb,ia)
        i2=JJ_tmp(2,j,iorb,ia)
        i3=JJ_tmp(3,j,iorb,ia)
        k1=JJ_tmp(4,j,iorb,ia)
        k2=JJ_tmp(5,j,iorb,ia)
        k3=JJ_tmp(6,j,iorb,ia)
        d1=c1*i1+k1
        d2=c2*i2+k2
        d3=c3*i3+k3
        x = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
        y = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
        z = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
        uV(j,lma) = uV_tmp(j,iorb,ia)*Ylm(x,y,z,l,m)
        JJ_MAP(1:6,j,lma) = JJ_tmp(1:6,j,iorb,ia)
      end do ! j
    end do ! lma
!$omp end parallel do

    deallocate( MJJ_tmp )
    deallocate( JJ_tmp )
    deallocate( uV_tmp )
    deallocate( maps_tmp )
    deallocate( itmp1 )

    !call watchb( ttmp, tttt(:,7), barrier='on' )

    allocate( icheck_tmp4(a1b:b1b,a2b:b2b,a3b:b3b) )
    icheck_tmp4=0
    do lma=1,nzlma
      j=0
      icheck_tmp4=0
      do i=1,MJJ_MAP(lma)
        i1=JJ_MAP(1,i,lma)
        i2=JJ_MAP(2,i,lma)
        i3=JJ_MAP(3,i,lma)
        if ( icheck_tmp4(i1,i2,i3) == 0 .and. abs(uV(i,lma)) > 1.0d-14 ) then
          j=j+1
          icheck_tmp4(i1,i2,i3)=j
        end if
      end do
      MJJ(lma)=j
    end do
    deallocate( icheck_tmp4 )

! ---

    nl_max_send = maxval( lma_nsend_tmp )

    if ( allocated(lma_nsend) ) deallocate( lma_nsend )
    if ( allocated(sendmap)   ) deallocate( sendmap )
    if ( allocated(recvmap)   ) deallocate( recvmap )

    allocate( lma_nsend(0:nprocs_g-1) ) ; lma_nsend=0
    allocate( sendmap(nl_max_send,0:nprocs_g-1) ) ; sendmap=0
    allocate( recvmap(nl_max_send,0:nprocs_g-1) ) ; recvmap=0

! ---

    do n=0,nprocs_g-1
      sendmap(1:nl_max_send,n) = sendmap_tmp(1:nl_max_send,n)
      lma_nsend(n) = lma_nsend_tmp(n)
    end do

! ---

    allocate( ireq(2*nprocs_g) )
    allocate( istatus(MPI_STATUS_SIZE,2*nprocs_g) )
    nreq=0
    do n=0,nprocs_g-1
      if ( lma_nsend(n)<=0 .or. n==myrank_g ) cycle
      nreq=nreq+1
      call MPI_Isend(recvmap_tmp(1,n),lma_nsend(n), &
      MPI_INTEGER,n,1,comm_grid,ireq(nreq),ierr)
      nreq=nreq+1
      call MPI_Irecv(recvmap(1,n) ,lma_nsend(n), &
      MPI_INTEGER,n,1,comm_grid,ireq(nreq),ierr)
    end do
    call MPI_Waitall(nreq,ireq,istatus,ierr)
    deallocate( istatus )
    deallocate( ireq )

    deallocate( recvmap_tmp,sendmap_tmp,lma_nsend_tmp )

    allocate( LLp(3,0:nprocs_g-1) )
    n=-1
    do i3=0,node_partition(3)-1
    do i2=0,node_partition(2)-1
    do i1=0,node_partition(1)-1
      n=n+1
      LLp(1,n)=i1
      LLp(2,n)=i2
      LLp(3,n)=i3
    end do
    end do
    end do

    allocate( itmp(3,nrlma) ) ; itmp=0
    allocate( itmp1(nrlma), work(nrlma) )
    allocate( itmp2(nrlma),itmp3(3,nrlma) )

    do irlma=1,nrlma
      n=nl_rank_map(irlma)
      itmp(1,irlma)=LLp(1,n)-LLp(1,myrank_g)
      itmp(2,irlma)=LLp(2,n)-LLp(2,myrank_g)
      itmp(3,irlma)=LLp(3,n)-LLp(3,myrank_g)
    end do

    nrlma_xyz(1:6)=0

    m=0
    n=0
    do i=1,nrlma
      if( itmp(2,i)==0 .and. itmp(3,i)==0 .and. itmp(1,i)>0 )then
        n=n+1
        work(n)=itmp(1,i)
        itmp2(n)=i
      end if
    end do
    if ( n>0 ) then
      call indexx(n,work,itmp1)
      do i=1,n
        j=itmp2( itmp1(i) )
        itmp3(:,m+i)=itmp(:,j)
      end do
    end if
    m=m+n
    nrlma_xyz(1)=nrlma_xyz(1)+n
    n=0
    do i=1,nrlma
      if( itmp(2,i)==0 .and. itmp(3,i)==0 .and. itmp(1,i)<0 )then
        n=n+1
        work(n)=itmp(1,i)
        itmp2(n)=i
      end if
    end do
    if ( n>0 ) then
      call indexx(n,work,itmp1)
      do i=1,n
        j=itmp2(itmp1(i))
        itmp3(:,m+n-i+1)=itmp(:,j)
      end do
    end if
    m=m+n
    nrlma_xyz(2)=nrlma_xyz(2)+n

    n=0
    do i=1,nrlma
      if( itmp(1,i)==0 .and. itmp(3,i)==0 .and. itmp(2,i)>0 )then
        n=n+1
        work(n)=itmp(2,i)
        itmp2(n)=i
      end if
    end do
    if ( n>0 ) then
      call indexx(n,work,itmp1)
      do i=1,n
        j=itmp2( itmp1(i) )
        itmp3(:,m+i)=itmp(:,j)
      end do
    end if
    m=m+n
    nrlma_xyz(3)=nrlma_xyz(3)+n
    n=0
    do i=1,nrlma
      if( itmp(1,i)==0 .and. itmp(3,i)==0 .and. itmp(2,i)<0 )then
        n=n+1
        work(n)=itmp(2,i)
        itmp2(n)=i
      end if
    end do
    if ( n>0 ) then
      call indexx(n,work,itmp1)
      do i=1,n
        j=itmp2(itmp1(i))
        itmp3(:,m+n-i+1)=itmp(:,j)
      end do
    end if
    m=m+n
    nrlma_xyz(4)=nrlma_xyz(4)+n

    n=0
    do i=1,nrlma
      if( itmp(1,i)==0 .and. itmp(2,i)==0 .and. itmp(3,i)>0 )then
        n=n+1
        work(n)=itmp(3,i)
        itmp2(n)=i
      end if
    end do
    if ( n>0 ) then
      call indexx(n,work,itmp1)
      do i=1,n
        j=itmp2( itmp1(i) )
        itmp3(:,m+i)=itmp(:,j)
      end do
    end if
    m=m+n
    nrlma_xyz(5)=nrlma_xyz(5)+n
    n=0
    do i=1,nrlma
      if( itmp(1,i)==0 .and. itmp(2,i)==0 .and. itmp(3,i)<0 )then
        n=n+1
        work(n)=itmp(3,i)
        itmp2(n)=i
      end if
    end do
    if ( n>0 ) then
      call indexx(n,work,itmp1)
      do i=1,n
        j=itmp2(itmp1(i))
        itmp3(:,m+n-i+1)=itmp(:,j)
      end do
    end if
    m=m+n
    nrlma_xyz(6)=nrlma_xyz(6)+n


    n=maxval( nrlma_xyz )
    if ( allocated(num_2_rank) ) deallocate( num_2_rank )
    allocate( num_2_rank(n,6) )
    num_2_rank(:,:)=MPI_PROC_NULL

    m=0
    do i=1,6
      do j=1,nrlma_xyz(i)
        m=m+1
        i1=itmp3(1,m)+LLp(1,myrank_g)
        i2=itmp3(2,m)+LLp(2,myrank_g)
        i3=itmp3(3,m)+LLp(3,myrank_g)
        k = i1 + i2*np1 + i3*np1*np2
        num_2_rank(j,i)=k
      end do
    end do

    deallocate( itmp,itmp1,itmp2,itmp3,work )
    deallocate( LLp )

    do i=1,5,2
      n=max( nrlma_xyz(i),nrlma_xyz(i+1) )
      nrlma_xyz(i)=n
      nrlma_xyz(i+1)=n
    end do

    if ( use_real8_wf() ) then
      call d_allocate_ps_nloc2( MBZ_0, MBZ_1 )
    else
      call z_allocate_ps_nloc2( MBZ_0, MBZ_1 )
    end if
    !call watchb( ttmp, tttt(:,8), barrier='on' )
    ! MAXMJJ=0; if ( nzlma > 0 ) MAXMJJ=maxval( MJJ(1:nzlma) )
    ! call allocate_ps_nloc2( MB_d_nl, itype_nl_sendrecv=2 )
    ! if ( allocated(JJP) ) deallocate( JJP )
    ! ! if ( allocated(uVk) ) deallocate( uVk )
    ! if ( allocated(d_uVk) ) deallocate( d_uVk )
    ! if ( allocated(z_uVk) ) deallocate( z_uVk )
    ! allocate( JJP(MAXMJJ,nzlma) ); JJP=0
    ! !allocate( uVk(MAXMJJ,nzlma,MBZ_0:MBZ_1) ); uVk=(0.0d0,0.0d0)
    ! allocate( d_uVk(MAXMJJ,nzlma,MBZ_0:MBZ_1) ); d_uVk=0.0d0
    ! allocate( z_uVk(MAXMJJ,nzlma,MBZ_0:MBZ_1) ); z_uVk=(0.0d0,0.0d0)
    !if ( disp_sw ) write(*,*) "(JJP)"
    !call check_memory( 4.0d0,MAXMJJ,nzlma )
    !if ( disp_sw ) write(*,*) "(uVk)"
    !call check_memory( 8.0d0,MAXMJJ,MBZ_1-MBZ_0+1 )
    !call watchb( ttmp, tttt(:,9), barrier='on' )

    call prep_uvk_ps_nloc2( MBZ_0, MBZ_1, kbb(1,MBZ_0) )

    if ( .not.FLAG_KEEP_JJ_MAP ) deallocate( JJ_MAP )
    if ( .not.FLAG_KEEP_uV ) deallocate( uV )

    if ( use_real8_wf() ) then
      call sort_index_sub( MJJ, JJP, d_uVk )
    else
      call sort_index_sub( MJJ, JJP, z_uVk )
    end if

    !call watchb( ttmp, tttt(:,10), barrier='on' )

    call reset_init_flag_ps_nloc2_op

    call write_border( 0, " prep_ps_nloc2(end)" )

    !call watchb( ttmp, tttt(:,11), barrier='on' )

! ---

    !if ( myrank == 0 ) then
    !  do i=1,11
    !    write(*,'("time_prep_ps_nloc2(",i2,")",2f10.5)') i, tttt(:,i)
    !  end do
    !end if

  end subroutine prep_ps_nloc2


  subroutine prep_uvk_ps_nloc2( k0, k1, kbb )
    use rgrid_variables, only: Ngrid, Igrid
    use ps_nloc2_variables, only: d_uVk, z_uVk, JJP
    use var_sys_parameter, only: use_real8_wf
    implicit none
    integer,intent(in) :: k0,k1
    real(8),intent(in) :: kbb(3,k0:k1)
    integer :: ofs,a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab2,ab3
    integer :: i,j,k,j3,lma,i1,i2,i3,m1,m2,m3
    integer,allocatable :: icheck_tmp4(:,:,:)
    real(8) :: c1,c2,c3,d1,d2,d3,pi2,kr,u3
    complex(8) :: ztmp0

    ofs = Igrid(1,0)
    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)

    c1=1.0d0/Ngrid(1)
    c2=1.0d0/Ngrid(2)
    c3=1.0d0/Ngrid(3)

    ab1=b1b-a1b+1
    ab2=b2b-a2b+1
    ab3=b3b-a3b+1

    pi2 = 2.0d0*acos(-1.0d0)

    allocate( icheck_tmp4(a1b:b1b,a2b:b2b,a3b:b3b) )
    icheck_tmp4=0

    if ( use_real8_wf() ) then

      do k=k0,k1
        do lma=1,nzlma
          j=0
          icheck_tmp4=0
          do i=1,MJJ_MAP(lma)
            i1=JJ_MAP(1,i,lma)
            i2=JJ_MAP(2,i,lma)
            i3=JJ_MAP(3,i,lma)
            m1=JJ_MAP(4,i,lma)
            m2=JJ_MAP(5,i,lma)
            m3=JJ_MAP(6,i,lma)
            j3=icheck_tmp4(i1,i2,i3)
            u3=uV(i,lma)
            if ( j3 == 0 .and. abs(u3) > 1.0d-14 ) then
              j=j+1
              icheck_tmp4(i1,i2,i3)=j
              d_uVk(j,lma,k) = u3
              JJP(j,lma) = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ofs
            else if ( j3 /= 0 ) then
              d_uVk(j3,lma,k) = d_uVk(j3,lma,k) + u3
            end if
          end do !i
        end do !lma
      end do !k

    else ![ use complex16 wf ]

      do k=k0,k1
        d1=pi2*kbb(1,k)
        d2=pi2*kbb(2,k)
        d3=pi2*kbb(3,k)
        do lma=1,nzlma
          j=0
          icheck_tmp4=0
          do i=1,MJJ_MAP(lma)
            i1=JJ_MAP(1,i,lma)
            i2=JJ_MAP(2,i,lma)
            i3=JJ_MAP(3,i,lma)
            m1=JJ_MAP(4,i,lma)
            m2=JJ_MAP(5,i,lma)
            m3=JJ_MAP(6,i,lma)
            j3=icheck_tmp4(i1,i2,i3)
            u3=uV(i,lma)
            kr=d1*(c1*i1+m1)+d2*(c2*i2+m2)+d3*(c3*i3+m3)
            ztmp0=dcmplx(cos(kr),-sin(kr))*u3
            if ( j3 == 0 .and. abs(u3) > 1.0d-14 ) then
              j=j+1
              icheck_tmp4(i1,i2,i3)=j
              z_uVk(j,lma,k) = ztmp0
              JJP(j,lma) = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ofs
            else if ( j3 /= 0 ) then
              z_uVk(j3,lma,k) = z_uVk(j3,lma,k) + ztmp0
            end if
          end do !i
        end do !lma
      end do !k

    end if

    deallocate( icheck_tmp4 )

  end subroutine prep_uvk_ps_nloc2


  subroutine prep_rvk_ps_nloc2( k0, k1, kbb )
    use rgrid_variables, only: Igrid, Ngrid
    use ps_nloc2_variables, only: z_uVk, xVk, yVk, zVk
    use var_sys_parameter, only: use_real8_wf
    use aa_module, only: aa
    implicit none
    integer,intent(in) :: k0,k1
    real(8),intent(in) :: kbb(3,k0:k1)
    integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab2,ab3
    integer :: i,j,k,j3,lma,i1,i2,i3,i4,m1,m2,m3,m4
    integer,allocatable :: icheck_tmp4(:,:,:)
    real(8) :: a1,a2,a3,c1,c2,c3,d1,d2,d3,pi2,kr,u3,x,y,z
    complex(8) :: ztmp0
    complex(8),parameter :: z0=(0.0d0,0.0d0)

    if ( use_real8_wf() ) return

    m1 =   size( z_uVk, 1 )
    m2 =   size( z_uVk, 2 )
    m3 = lbound( z_uVk, 3 )
    m4 = ubound( z_uVk, 3 )

    if ( allocated(xVk) ) then
      i1 = size(xVk,1)
      i2 = size(xVk,2)
      i3 = lbound(xVk,3)
      i4 = ubound(xVk,3)
      if ( i1/=m1 .or. i2/=m2 .or. i3/=m3 .or. i4/=m4 ) deallocate( zVk, yVk, xVk )
    end if

    if ( .not.allocated(xVk) ) then
      allocate( xVk(m1,m2,m3:m4) ); xVk=z0
      allocate( yVk(m1,m2,m3:m4) ); yVk=z0
      allocate( zVk(m1,m2,m3:m4) ); zVk=z0
    end if

    xVk=z0
    yVk=z0
    zVk=z0

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = b1b-a1b+1
    ab2 = b2b-a2b+1
    ab3 = b3b-a3b+1

    c1 = 1.0d0/Ngrid(1)
    c2 = 1.0d0/Ngrid(2)
    c3 = 1.0d0/Ngrid(3)

    pi2 = 2.0d0*acos(-1.0d0)

    allocate( icheck_tmp4(a1b:b1b,a2b:b2b,a3b:b3b) )
    icheck_tmp4=0

    do k=k0,k1
      d1=pi2*kbb(1,k)
      d2=pi2*kbb(2,k)
      d3=pi2*kbb(3,k)
      do lma=1,nzlma
        j=0
        icheck_tmp4=0
        do i=1,MJJ_MAP(lma)
          i1=JJ_MAP(1,i,lma)
          i2=JJ_MAP(2,i,lma)
          i3=JJ_MAP(3,i,lma)
          m1=JJ_MAP(4,i,lma)
          m2=JJ_MAP(5,i,lma)
          m3=JJ_MAP(6,i,lma)
          a1=c1*i1+m1
          a2=c2*i2+m2
          a3=c3*i3+m3
          j3=icheck_tmp4(i1,i2,i3)
          u3=uV(i,lma)
          kr=d1*a1+d2*a2+d3*a3
          ztmp0=dcmplx(cos(kr),-sin(kr))*u3
          x=a1*aa(1,1)+a2*aa(1,2)+a3*aa(1,3)
          y=a1*aa(2,1)+a2*aa(2,2)+a3*aa(2,3)
          z=a1*aa(3,1)+a2*aa(3,2)+a3*aa(3,3)
          if ( j3 == 0 .and. abs(u3) > 1.0d-14 ) then
            j=j+1
            icheck_tmp4(i1,i2,i3)=j
            xVk(j,lma,k) = x*ztmp0
            yVk(j,lma,k) = y*ztmp0
            zVk(j,lma,k) = z*ztmp0
          else if ( j3 /= 0 ) then
            xVk(j3,lma,k) = xVk(j3,lma,k) + x*ztmp0
            yVk(j3,lma,k) = yVk(j3,lma,k) + y*ztmp0
            zVk(j3,lma,k) = zVk(j3,lma,k) + z*ztmp0
          end if
        end do !i
      end do !lma
    end do !k

    deallocate( icheck_tmp4 )

  end subroutine prep_rvk_ps_nloc2


  subroutine construct_lma_map_ps_nloc2( map )
    implicit none
    integer,intent(inout),allocatable :: map(:)
    integer :: a, ik, iorb, l, m, lma, i, Natom
    Natom=size(ki_atom)
    allocate( map(nzlma) ); map=0
    lma=0
    do a=1,Natom
      ik=ki_atom(a)
    do iorb=1,norb(ik)
      l=lo(iorb,ik)
    do m=-l,l
      lma=lma+1
      do i=1,nzlma
        if ( amap(i)==a .and. iorbmap(i)==iorb .and. mmap(i)==m ) then
          map(i)=lma
          exit
        end if
      end do
    end do ! m
    end do ! iorb
    end do ! a
  end subroutine construct_lma_map_ps_nloc2


end module ps_nloc2_module
