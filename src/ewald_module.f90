MODULE ewald_module

  use ewald_variables, only: eta,mg,mr,LG,LR,ipair,mpair,rrcut,ecut
  use aa_module, only: aa, Va
  use bb_module, only: bb
  use atom_module, only: Natom, aa_atom, ki_atom
  use pseudopot_module, only: Zps
  use parallel_module
  use watch_module
  use hsort_module
  use bberf_module
  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: test_ewald, calc_ewald, cal_ewald, init_ewald

  real(8),allocatable :: zatom(:)
  real(8) :: Qtot,Qtot2,Vcell
  integer,allocatable :: id(:),ir(:)
  real(8) :: ewldg_0=0.0d0, g_0=0.0d0, g_0_0=0.0d0
  real(8) :: ewldr_0=0.0d0, r_0=0.0d0, r_0_0=0.0d0, r_1=0.0d0
  real(8) :: xyz_min(3), xyz_max(3)
  integer :: irange_bak(2,3)=reshape( (/0,0,0,0,0,0/),(/2,3/) )

  real(8) :: factor_eta=10.0d0
  logical :: flag_init=.false.
  logical :: flag_first_time=.true.
  character(12) :: ctrl_ewald=""

!#ifdef _NO_QPRECISION_
  integer,parameter :: QP=kind(0.0d0)
!#else
!  integer,parameter :: QP=kind(0.0q0)
!#endif
  real(QP) :: ewldg16, ewldr16

CONTAINS


  SUBROUTINE init_ewald
    implicit none
    if ( flag_init ) return
    flag_init=.true.
    call IOTools_readReal8Keyword( "FACTEWALD", factor_eta )
    call IOTools_readStringKeyword( "CTRLEWALD", ctrl_ewald )
  END SUBROUTINE init_ewald


  SUBROUTINE cal_ewald(Ewld)
    implicit none
    real(8),intent(OUT) :: Ewld
    real(8) :: ewldg,ewldr,pi
    call write_border( 0, "cal_ewald(start)" )
    pi=acos(-1.d0)
    call calc_ewldg( mg, LG, ewldg, ewldg16 )
    ewldg=ewldg*4.d0*pi/Vcell
    ewldg=ewldg-pi/(eta*Vcell)*Qtot**2
    ewldg=ewldg-2.d0*sqrt(eta/pi)*Qtot2
    ewldr=ewldg
    ewldg=0.0d0
    call calc_ewldr( mr, LR, ewldr )
    Ewld=0.5d0*(ewldg+ewldr)
    call write_border( 0, "cal_ewald(end)" )
  END SUBROUTINE cal_ewald


  SUBROUTINE test_ewald(Ewld)
    implicit none
    real(8),intent(OUT) :: Ewld
    integer,parameter :: maxloop=50,max_loop_r=10,max_loop_g=10
    integer :: loop_eta,loop_g,loop_r,mg0,mr0,i,m1,m2,m3,m,n,j
    integer :: mg_tot_0,mr_tot_0,mr_tot,mg_tot,mr_tmp,mr_dif,ierr
    integer,allocatable :: LR_tmp(:,:)
    real(8),parameter :: epsg=1.d-12, epsr=1.d-12, ep=1.d-10
    real(8),parameter :: factor=1.0d0
    real(8) :: const,const1,ewldr,ewldr0,ewldg,ewldg0,eta_small,eta_big
    real(8) :: r,r0,rr,alpha,g,g0,t,pi,max_bb,max_aa,eta_0,tt0,ttt,gg
    real(8) :: ewldr_tmp,gmax
    real(8) :: ct0,et0,ct1,et1,timg(2),timg0(2),timr(2),timr0(2)
    logical :: flag_conv
    integer,allocatable :: grd_chk(:,:,:),grd_chk0(:,:,:)

    pi=acos(-1.d0)

    Vcell = abs( aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
                +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
                -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2) )

    max_aa=0.d0
    do i=1,3
       t=sqrt( sum(aa(1:3,i)**2) )
       max_aa=max(max_aa,t)
    end do
    max_bb=0.d0
    do i=1,3
       t=sqrt( sum(bb(1:3,i)**2) )
       max_bb=max(max_bb,t)
    end do

    if ( .not.allocated(zatom) ) then
       allocate( zatom(Natom) )
       do i=1,Natom
          zatom(i)=Zps(ki_atom(i))
       end do
    end if

    if ( .not.allocated(id) ) then
       allocate( id(0:nprocs-1) )
       allocate( ir(0:nprocs-1) )
    end if

    if ( .not.allocated(ipair) ) then
       mpair = ( Natom*(Natom+1) )/2
       call prep_parallel(mpair,nprocs,ir,id)
       mpair = ir(myrank)
       allocate( ipair(2,mpair) )
       m=0
       n=0
       do j=1,Natom
       do i=j,Natom
          m=m+1
          if ( m < id(myrank)+1 .or. id(myrank)+ir(myrank) < m ) cycle
          n=n+1
          ipair(1,n)=i
          ipair(2,n)=j
       end do
       end do
    end if

    eta = pi/Vcell**(2.d0/3.d0)

    eta_small = 0.d0
    eta_big   = 0.d0
    eta_0     = 0.d0
    tt0       = 0.d0
    mg_tot_0  = 0
    mr_tot_0  = 0

    Qtot =0.d0
    Qtot2=0.d0
    do i=1,Natom
       Qtot =Qtot +zatom(i)
       Qtot2=Qtot2+zatom(i)**2
    end do

    const1 = 4.d0*Pi/Vcell

    do loop_eta=1,maxloop

!--- G

       const=1.d0/(4.d0*eta)
!       do i=1,10000
!          g=i
!          t=const1*exp(-const*g*g)/(g*g)
!          if ( t <= epsg ) exit
!       end do
       g=12.d0*sqrt(eta)

       ewldg0 = 0.d0
       flag_conv = .false.
       do loop_g=1,max_loop_g

          ecut=g*g

          call search_grid(bb,ecut,m1,m2,m3,mg_tot)
          call prep_parallel(mg_tot,nprocs,ir,id)
          mg=ir(myrank)
          if ( allocated(LG) ) deallocate( LG )
          allocate( LG(3,mg) )
          call construct_grid(bb,ecut,mg_tot,m1,m2,m3,mg,LG)

          call watch(ct0,et0)
          call calc_ewldg(mg,LG,ewldg)
          ewldg=ewldg*4.d0*pi/Vcell
          ewldg=ewldg-pi/(eta*Vcell)*Qtot**2
          ewldg=ewldg-2.d0*sqrt(eta/pi)*Qtot2
          call watch(ct1,et1) ; timg(1)=ct1-ct0 ; timg(2)=et1-et0

          ir(myrank)=mg
          call mpi_allgather &
               (mg,1,mpi_integer,ir,1,mpi_integer,mpi_comm_world,ierr)

          if ( disp_switch_parallel ) then
             write(*,'(1x,2i4,2i8,2g26.15,2g12.5)') &
                  loop_eta,loop_g,mg_tot,mg,ewldg,sqrt(ecut),timg(1:2)
          end if

          if ( abs((ewldg-ewldg0)/ewldg) < 1.d-15 ) then
             flag_conv=.true.
             ewldg=ewldg0
             mg=mg0
             g=g0
             timg(:)=timg0(:)
             exit
          else
             ewldg0=ewldg
             mg0=mg
             g0=g
             g=(g+max_bb)*factor
             timg0(:)=timg(:)
          end if

       end do ! loop_g

       if ( flag_conv ) then
          if ( mg <= 0 ) then
             gg=0.d0
          else
             gg=( bb(1,1)*LG(1,mg)+bb(1,2)*LG(2,mg)+bb(1,3)*LG(3,mg) )**2 &
               +( bb(2,1)*LG(1,mg)+bb(2,2)*LG(2,mg)+bb(2,3)*LG(3,mg) )**2 &
               +( bb(3,1)*LG(1,mg)+bb(3,2)*LG(2,mg)+bb(3,3)*LG(3,mg) )**2
          end if
          call mpi_allreduce(gg,ecut,1,mpi_real8,mpi_max,mpi_comm_world,i)
          call search_grid(bb,ecut,m1,m2,m3,mg_tot)
          call prep_parallel(mg_tot,nprocs,ir,id)
          mg=ir(myrank)
          if ( allocated(LG) ) deallocate( LG )
          allocate( LG(3,mg) )
          call construct_grid(bb,ecut,mg_tot,m1,m2,m3,mg,LG)
       else
          stop "ewldg is not converged"
       end if

       if ( disp_switch_parallel ) then
          write(*,'(1x,"ewald(G)  =",2g24.15,2i8)') ewldg*0.5d0,sqrt(ecut),mg,mg_tot
       end if

!--- R

       alpha = sqrt(eta)
!       do i=1,1000
!          r=i
!          t=bberfc(alpha*r)/r
!          if ( t <= epsr ) exit
!       end do
       r=6.d0/sqrt(eta)

       ewldr0    = 0.d0
!       ewldr     = 0.d0
       ewldr     = ewldg ; ewldg=0.d0
       ewldr_tmp = ewldr
       flag_conv = .false.
       mr_tmp    = 0
       timr(:)   = 0.d0
       do loop_r=1,max_loop_r

          rrcut=r*r

          call watch(ct0,et0)

          call search_grid(aa,rrcut,m1,m2,m3,mr_tot)
!          call prep_parallel(mr_tot,nprocs,ir,id) ; mr=ir(myrank)
          mr=mr_tot ; id(:)=0 ; ir(:)=mr

          allocate( grd_chk(-m1:m1,-m2:m2,-m3:m3) ) ; grd_chk=0

          if ( loop_r == 1 ) then
          else
             allocate( grd_chk0(-m1:m1,-m2:m2,-m3:m3) ) ; grd_chk0=0
             mr_tmp = size(LR,2)
             do i=1,mr_tmp
                grd_chk0( LR(1,i),LR(2,i),LR(3,i) ) = 1
             end do
             m=size(grd_chk)
             call mpi_allreduce(grd_chk0,grd_chk,m &
                               ,mpi_integer,mpi_sum,mpi_comm_world,ierr)
             deallocate( grd_chk0 )
          end if

          if ( allocated(LR) ) deallocate(LR)
          allocate( LR(3,mr) )

          call construct_grid(aa,rrcut,mr_tot,m1,m2,m3,mr,LR)

          if ( allocated(LR_tmp) ) deallocate( LR_tmp )
          allocate( LR_tmp(3,mr) ) ; LR_tmp=0

          if ( loop_r == 1 ) then
             m=mr
             LR_tmp(1:3,1:mr) = LR(1:3,1:mr)
          else
             m=0
             do i=1,mr
                if ( grd_chk(LR(1,i),LR(2,i),LR(3,i)) == 0 ) then
                   m=m+1
                   LR_tmp(1:3,m) = LR(1:3,i)
                end if
             end do
          end if

          deallocate( grd_chk )

          call watch(ct1,et1)
!(1)
!          call watch(ct0,et0)
!          ewldr=ewldr_tmp
!!          ewldr=0.d0
!          call calc_ewldr(mr,LR,ewldr)
!          call watch(ct1,et1)
!          timr(1)=ct1-ct0 ; timr(2)=et1-et0
!(2)
          call watch(ct0,et0)
          call calc_ewldr(m,LR_tmp,ewldr)
          call watch(ct1,et1)
          timr(1)=timr(1)+ct1-ct0 ; timr(2)=timr(2)+et1-et0

          deallocate( LR_tmp )

          if ( disp_switch_parallel ) then
             write(*,'(1x,2i4,3i8,2g26.15,2g12.5)') &
                  loop_eta,loop_r,mr_tot,mr,m,ewldr,sqrt(rrcut),timr(1:2)
          end if

          if ( abs((ewldr-ewldr0)/ewldr) < 1.d-15 ) then
             flag_conv=.true.
             ewldr=ewldr0
             mr=mr0
             r=r0
             timr(:)=timr0(:)
             exit
          else
             ewldr0=ewldr
             mr0=mr
             r0=r
             r=(r+max_aa)*factor
             timr0(:)=timr(:)
          end if

       end do ! loop_r

       if ( flag_conv ) then
          rr=( aa(1,1)*LR(1,mr)+aa(1,2)*LR(2,mr)+aa(1,3)*LR(3,mr) )**2 &
            +( aa(2,1)*LR(1,mr)+aa(2,2)*LR(2,mr)+aa(2,3)*LR(3,mr) )**2 &
            +( aa(3,1)*LR(1,mr)+aa(3,2)*LR(2,mr)+aa(3,3)*LR(3,mr) )**2
          call mpi_allreduce(rr,rrcut,1,mpi_real8,mpi_max,mpi_comm_world,i)
          call search_grid(aa,rrcut,m1,m2,m3,mr_tot)
!          call prep_parallel(mr_tot,nprocs,ir,id) ; mr=ir(myrank)
          mr=mr_tot ; id(:)=0 ; ir(:)=mr
          if ( allocated(LR) ) deallocate( LR )
          allocate( LR(3,mr) )
          call construct_grid(aa,rrcut,mr_tot,m1,m2,m3,mr,LR)
       else
          stop "ewldr is not converged"
       end if

       if ( disp_switch_parallel ) then
          write(*,'(1x,"ewald(G+R)=",2g24.15,2i8)') ewldr*0.5d0,sqrt(rrcut),mr,mr_tot
       end if

       Ewld=0.5d0*(ewldr+ewldg)

       call mpi_allreduce(timg,timg0,2,mpi_real8,mpi_max,mpi_comm_world,i)
       call mpi_allreduce(timr,timr0,2,mpi_real8,mpi_max,mpi_comm_world,i)

       if ( disp_switch_parallel ) then
          write(*,'(1x,i6,4x,"Ewld=",3g22.15,2i8,2x,3g12.5)') &
               loop_eta,Ewld,eta,eta-eta_0,mg_tot,mr_tot &
               ,timg0(2),timr0(2),timg0(2)+timr0(2)
       end if

       if ( abs(eta-eta_0) < ep ) exit
       if ( mg_tot_0 == mg_tot .and. mr_tot_0 == mr_tot ) exit

       mg_tot_0 = mg_tot
       mr_tot_0 = mr_tot

       ttt=timg0(2)+timr0(2)
       if ( loop_eta == 1 ) then
          eta_0 = eta
          tt0 = ttt
          eta = eta*2.d0
       else
          if ( eta_big == 0.d0 ) then
             if ( ttt < tt0 ) then
                eta_small = eta_0
                eta_0 = eta
                tt0 = ttt
                eta = eta*2.d0
             else
                eta_big = eta
                if ( mod(loop_eta,2) == 0 ) then
                   eta=(eta_small+eta_0)/2
                else
                   eta=(eta_0+eta_big)/2
                end if
             end if
          else
             if ( mod(loop_eta,2) == 1 ) then
                if ( ttt > tt0 ) then
                   eta_small = eta
                   eta = (eta_0+eta_big)/2
                else
                   eta_big = eta_0
                   eta_0 = eta
                   tt0 = ttt
                   eta = (eta_0+eta_big)/2
                end if
             else
                if ( ttt > tt0 ) then
                   eta_big = eta
                   eta = (eta_small+eta_0)/2
                else
                   eta_small = eta_0
                   eta_0 = eta
                   tt0 = ttt
                   eta = (eta_small+eta_0)/2
                end if
             end if
          end if
       end if

       if ( loop_eta == maxloop ) eta=eta_0

    end do ! loop_eta

    return
  END SUBROUTINE test_ewald


  SUBROUTINE search_grid(v,vvcut,m1,m2,m3,mv)
    implicit none
    real(8),intent(IN)  :: v(3,3),vvcut
    integer,intent(OUT) :: m1,m2,m3,mv
    integer :: i,m,i1,i2,i3
    real(8) :: vx,vy,vz,vv,vv1,vv2,vv3
    vv1=sum(v(:,1)**2)
    vv2=sum(v(:,2)**2)
    vv3=sum(v(:,3)**2)
    m1=nint( sqrt(vvcut/vv1) )
    m2=nint( sqrt(vvcut/vv2) )
    m3=nint( sqrt(vvcut/vv3) )
    mv=0
    do i=1,100000
       m=0
       do i3=-m3,m3
       do i2=-m2,m2
       do i1=-m1,m1
          vx=v(1,1)*i1+v(1,2)*i2+v(1,3)*i3
          vy=v(2,1)*i1+v(2,2)*i2+v(2,3)*i3
          vz=v(3,1)*i1+v(3,2)*i2+v(3,3)*i3
          vv=vx*vx+vy*vy+vz*vz
          if ( vv < vvcut + 1.d-8 ) then
             m=m+1
             cycle
          end if
          vx=v(1,1)*i1+v(1,2)*i2+v(1,3)*i3-xyz_min(1)
          vy=v(2,1)*i1+v(2,2)*i2+v(2,3)*i3-xyz_min(2)
          vz=v(3,1)*i1+v(3,2)*i2+v(3,3)*i3-xyz_min(3)
          vv=vx*vx+vy*vy+vz*vz
          if ( vv < vvcut + 1.d-8 ) then
             m=m+1
             cycle
          end if
          vx=v(1,1)*i1+v(1,2)*i2+v(1,3)*i3+xyz_max(1)
          vy=v(2,1)*i1+v(2,2)*i2+v(2,3)*i3+xyz_max(2)
          vz=v(3,1)*i1+v(3,2)*i2+v(3,3)*i3+xyz_max(3)
          vv=vx*vx+vy*vy+vz*vz
          if ( vv < vvcut + 1.d-8 ) then
             m=m+1
             cycle
          end if
       end do
       end do
       end do
       if ( m == mv ) exit
       mv=m
       m1=m1+1
       m2=m2+1
       m3=m3+1
    end do
  END SUBROUTINE search_grid


  SUBROUTINE prep_parallel(mv,nprocs,ir,id)
    implicit none
    integer,intent(IN) :: mv,nprocs
    integer,intent(OUT) :: ir(0:nprocs-1),id(0:nprocs-1)
    integer :: i,j
    ir(:)=0
    do i=1,mv
       j=mod(i-1,nprocs)
       ir(j)=ir(j)+1
    end do
    do j=0,nprocs-1
       id(j)=sum(ir(0:j))-ir(j)
    end do
  END SUBROUTINE prep_parallel


  SUBROUTINE construct_grid(v,vvcut,mt,m1,m2,m3,mv,LV)
    implicit none
    real(8),intent(IN)  :: v(3,3),vvcut
    integer,intent(IN)  :: mt,m1,m2,m3,mv
    integer,intent(OUT) :: LV(3,mv)
    integer :: m,n,i1,i2,i3,i,j
    integer,allocatable :: indx(:),LVtmp(:,:)
    real(8) :: vx,vy,vz,vv
    real(8),allocatable :: ll(:)
    allocate( LVtmp(3,mt) ) ; LVtmp=0
    m=0
    n=0
    do i3=-m3,m3
    do i2=-m2,m2
    do i1=-m1,m1
       vx=v(1,1)*i1+v(1,2)*i2+v(1,3)*i3
       vy=v(2,1)*i1+v(2,2)*i2+v(2,3)*i3
       vz=v(3,1)*i1+v(3,2)*i2+v(3,3)*i3
       vv=vx*vx+vy*vy+vz*vz
       if ( vv < vvcut + 1.d-8 ) then
          m=m+1
          LVtmp(1,m)=i1
          LVtmp(2,m)=i2
          LVtmp(3,m)=i3
          if ( m < id(myrank)+1 .or. id(myrank)+ir(myrank) < m ) cycle
          n=n+1
          LV(1,n)=i1
          LV(2,n)=i2
          LV(3,n)=i3
          cycle
       end if
       vx=v(1,1)*i1+v(1,2)*i2+v(1,3)*i3-xyz_min(1)
       vy=v(2,1)*i1+v(2,2)*i2+v(2,3)*i3-xyz_min(2)
       vz=v(3,1)*i1+v(3,2)*i2+v(3,3)*i3-xyz_min(3)
       vv=vx*vx+vy*vy+vz*vz
       if ( vv < vvcut + 1.d-8 ) then
          m=m+1
          LVtmp(1,m)=i1
          LVtmp(2,m)=i2
          LVtmp(3,m)=i3
          if ( m < id(myrank)+1 .or. id(myrank)+ir(myrank) < m ) cycle
          n=n+1
          LV(1,n)=i1
          LV(2,n)=i2
          LV(3,n)=i3
          cycle
       end if
       vx=v(1,1)*i1+v(1,2)*i2+v(1,3)*i3+xyz_max(1)
       vy=v(2,1)*i1+v(2,2)*i2+v(2,3)*i3+xyz_max(2)
       vz=v(3,1)*i1+v(3,2)*i2+v(3,3)*i3+xyz_max(3)
       vv=vx*vx+vy*vy+vz*vz
       if ( vv < vvcut + 1.d-8 ) then
          m=m+1
          LVtmp(1,m)=i1
          LVtmp(2,m)=i2
          LVtmp(3,m)=i3
          if ( m < id(myrank)+1 .or. id(myrank)+ir(myrank) < m ) cycle
          n=n+1
          LV(1,n)=i1
          LV(2,n)=i2
          LV(3,n)=i3
          cycle
       end if
    end do
    end do
    end do
    if ( n /= mv ) then
       write(*,*) "n,m,mv,mt=",n,m,mv,mt,myrank
       stop "stop at construct_grid(ewald1)"
    end if
    if ( m /= mt ) then
       write(*,*) "n,m,mv,mt=",n,m,mv,mt,myrank
       stop "stop at construct_grid(ewald2)"
    end if

    allocate( ll(m),indx(m) )
    do i=1,m
       i1=LVtmp(1,i)
       i2=LVtmp(2,i)
       i3=LVtmp(3,i)
       vx=v(1,1)*i1+v(1,2)*i2+v(1,3)*i3
       vy=v(2,1)*i1+v(2,2)*i2+v(2,3)*i3
       vz=v(3,1)*i1+v(3,2)*i2+v(3,3)*i3
       ll(i)=vx*vx+vy*vy+vz*vz
    end do
    call indexx(m,ll,indx)
    if ( m == n ) then
       do i=1,m
          LV(1:3,i)=LVtmp(1:3,indx(i))
       end do
    else
    n=0
    do i=1,m
       j=mod(i-1,nprocs)
       if ( j == myrank ) then
          n=n+1
          LV(1:3,n)=LVtmp(1:3,indx(i))
       end if
    end do
    end if
    deallocate( indx,ll,LVtmp )

  END SUBROUTINE construct_grid


  SUBROUTINE calc_ewldg( mg, LG, ewldg, qewldg )
    implicit none
    integer,intent(INOUT) :: mg
    integer,intent(IN)    :: LG(3,mg)
    real(8),intent(OUT)   :: ewldg
    real(QP),optional,intent(OUT) :: qewldg
    integer :: i,a,m1
    real(8) :: x,y,z,t,pi2,gx,gy,gz,gg,e0,e1,g1,const
    complex(8) :: sg
    real(QP) :: qe0
    pi2=2.d0*acos(-1.d0)
    const=1.d0/(4.d0*eta)
    qe0=0.0_QP
    e0=0.d0
    e1=0.d0
    g1=0.d0
    m1=0
    do i=1,mg
       x=pi2*LG(1,i)
       y=pi2*LG(2,i)
       z=pi2*LG(3,i)
       if ( all( LG(1:3,i) == 0 ) ) cycle
       sg=(0.0d0,0.0d0)
!$OMP parallel do private( t ) reduction(+:sg)
       do a=1,Natom
          t=x*aa_atom(1,a)+y*aa_atom(2,a)+z*aa_atom(3,a)
          sg=sg+zatom(a)*dcmplx(cos(t),-sin(t))
       end do
!$OMP end parallel do
       gx=bb(1,1)*LG(1,i)+bb(1,2)*LG(2,i)+bb(1,3)*LG(3,i)
       gy=bb(2,1)*LG(1,i)+bb(2,2)*LG(2,i)+bb(2,3)*LG(3,i)
       gz=bb(3,1)*LG(1,i)+bb(3,2)*LG(2,i)+bb(3,3)*LG(3,i)
       gg=gx*gx+gy*gy+gz*gz
       t = exp(-const*gg)/gg
       e0 = e0 + abs(sg)**2*t
       if ( present(qewldg) ) qe0 = qe0 + abs(sg)**2*t
       if ( e0 /= e1 ) then
          e1=e0
          g1=gg
          m1=i
       end if
    end do
    call mpi_allreduce(e0,ewldg,1,mpi_real8,mpi_sum,mpi_comm_world,i)
    if ( present(qewldg) ) then
       e0=qe0
       call mpi_allreduce(e0,e1,1,mpi_real8,mpi_sum,mpi_comm_world,i)
       qewldg=e1
    end if
    mg=m1
  END SUBROUTINE calc_ewldg


  SUBROUTINE calc_ewldr( mr, LR, ewldr, qewldr, icheck )
    implicit none
    integer,intent(IN) :: mr,LR(3,mr)
    real(8),intent(INOUT) :: ewldr
    real(QP),optional,intent(INOUT) :: qewldr
    integer,optional,intent(INOUT) :: icheck(:)
    integer :: i,i1,i2,i1i2,ic1,ic2
    real(8) :: a1_2,a2_2,a3_2,a1_1,a2_1,a3_1,z1,z2
    real(8) :: sum0,sum1,sum2,a1,a2,a3,x,y,z,r,rr,t,alpha,rmax
    real(QP) :: qsum0, qsum1, qsum2
    alpha=sqrt(eta)
    sum2=0.d0
    qsum2=0.0_QP
    rmax=0.0d0
!    ic1=0
!    ic2=0
    do i1i2=1,mpair
       i1=ipair(1,i1i2)
       i2=ipair(2,i1i2)
       a1_2=aa_atom(1,i2)
       a2_2=aa_atom(2,i2)
       a3_2=aa_atom(3,i2)
       z2=zatom(i2)
       a1_1=aa_atom(1,i1)
       a2_1=aa_atom(2,i1)
       a3_1=aa_atom(3,i1)
       z1=zatom(i1)
       sum0=0.d0
       qsum0=0.0_QP
!$OMP parallel do private( a1,a2,a3,x,y,z,rr,r,t ) reduction(+:sum0)
       do i=1,mr
          a1=LR(1,i)+a1_1-a1_2
          a2=LR(2,i)+a2_1-a2_2
          a3=LR(3,i)+a3_1-a3_2
          x=aa(1,1)*a1+aa(1,2)*a2+aa(1,3)*a3
          y=aa(2,1)*a1+aa(2,2)*a2+aa(2,3)*a3
          z=aa(3,1)*a1+aa(3,2)*a2+aa(3,3)*a3
          rr=x*x+y*y+z*z
          if ( rr == 0.0d0 ) cycle
          r=sqrt(rr)
          t=bberfc(alpha*r)/r
          sum0=sum0+t
          !ic1=ic1+1
          !if ( t == 0.0d0 ) then
          !   ic2=ic2+1
          !else
          !   if ( present(icheck) ) icheck(i)=icheck(i)+1
          !   rmax=max(rmax,r)
          !end if
          !if ( present(qewldr) ) qsum0=qsum0+t
          if ( t/=0.0d0 .and. present(icheck) ) icheck(i)=icheck(i)+1
       end do
!$OMP end parallel do
       if ( i1 /= i2 ) sum0=sum0*2.d0
       sum2=sum2+z1*z2*sum0
       if ( present(qewldr) ) then
          if ( i1 /= i2 ) qsum0=qsum0*2.0_QP
          qsum2=qsum2+z1*z2*qsum0
       end if
    end do
    !if ( myrank == 0 ) write(*,'(1x,2i12,2f15.5," E")') ic1,ic2,rmax,sqrt(rrcut)
    call mpi_allreduce(rmax,r_1,1,mpi_real8,mpi_max,mpi_comm_world,i)
    call mpi_allreduce(sum2,sum1,1,mpi_real8,mpi_sum,mpi_comm_world,i)
    ewldr=ewldr+sum1
    if ( present(qewldr) ) then
       sum2=qsum2
       call mpi_allreduce(sum2,sum1,1,mpi_real8,mpi_sum,mpi_comm_world,i)
       qewldr=qewldr+sum1
    end if
  END SUBROUTINE calc_ewldr


  SUBROUTINE calc_ewald(Ewld)
    implicit none
    real(8),intent(OUT) :: Ewld
    integer,parameter :: maxloop=50,max_loop_r=10,max_loop_g=10
    integer :: i,m1,m2,m3,m,n,j,loop_g,loop_r,i1,i2,i3
    integer :: mg_tot_0,mr_tot_0,mr_tot,mg_tot,mr_tmp,mr_dif,ierr
    integer,allocatable :: LR_tmp(:,:)
    real(8),parameter :: epsg=1.d-12, epsr=1.d-12, ep=1.d-10
    real(8),parameter :: factor=1.0d0
    real(8) :: const,const1,ewldr,ewldg
    real(8) :: r,rr,alpha,g,t,pi,max_bb,max_aa,tt0,ttt,gg,min_aa,avg_aa
    real(8) :: ewldr_tmp,gmax,g_0_guess,r_0_guess
    real(8) :: ct0,et0,ct1,et1,timg(2,0:5),timr(2,2),ttmp(2),ttmp0(2)
    integer,allocatable :: grd_chk(:,:,:),grd_chk0(:,:,:)
    integer,allocatable :: LR_bak(:,:), icheck_LR_tmp(:)
    logical :: disp_sw
    integer :: irange(2,3)

    if ( ctrl_ewald == "FIXED_RANGE" ) then
       if ( flag_first_time ) then
       else
          call cal_ewald( Ewld )
          return
       end if
       flag_first_time = .false.
    end if

    call check_disp_switch( disp_sw, 0 )
    call  write_border( 0, "calc_ewald(start)" )

    pi=acos(-1.d0)

    Vcell = abs( aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
                +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
                -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2) )

    min_aa=1.0d100
    max_aa=0.0d0
    avg_aa=0.0d0
    do i=1,3
       t=sqrt( sum(aa(1:3,i)**2) )
       min_aa=min(min_aa,t)
       max_aa=max(max_aa,t)
       avg_aa=avg_aa+t/3.0d0
    end do

    max_bb=0.0d0
    do i=1,3
       t=sqrt( sum(bb(1:3,i)**2) )
       max_bb=max(max_bb,t)
    end do

    if ( .not.allocated(zatom) ) then
       allocate( zatom(Natom) )
       do i=1,Natom
          zatom(i)=Zps(ki_atom(i))
       end do
    end if

    if ( .not.allocated(id) ) then
       allocate( id(0:nprocs-1) )
       allocate( ir(0:nprocs-1) )
    end if

    if ( .not.allocated(ipair) ) then
       mpair = ( Natom*(Natom+1) )/2
       call prep_parallel(mpair,nprocs,ir,id)
       mpair = ir(myrank)
       allocate( ipair(2,mpair) )
       m=0
       n=0
       do j=1,Natom
       do i=j,Natom
          m=m+1
          if ( m < id(myrank)+1 .or. id(myrank)+ir(myrank) < m ) cycle
          n=n+1
          ipair(1,n)=i
          ipair(2,n)=j
       end do
       end do
    end if

    Qtot =0.d0
    Qtot2=0.d0
    do i=1,Natom
       Qtot =Qtot +zatom(i)
       Qtot2=Qtot2+zatom(i)**2
    end do

    const1 = 4.d0*Pi/Vcell

    if ( eta == 0.0d0 ) then
       eta = pi/Vcell**(2.d0/3.d0) * factor_eta
    end if

!--- G

    xyz_min(:)=0.0d0
    xyz_max(:)=0.0d0

    const=1.d0/(4.d0*eta)

    if ( g_0 == 0.0d0 ) then
       g_0=12.d0*sqrt(eta)
    else if ( g_0 > g_0_0*2 ) then
       g_0=12.d0*sqrt(eta)
    end if

    g=g_0

    do loop_g=1,max_loop_g

       ecut=g*g

       call watchb( ttmp0 )
       !call watchb( ttmp ) ; timg=0.0d0

       call search_grid( bb, ecut, m1,m2,m3, mg_tot )

       !call watchb( ttmp, timg(:,1) )

       call prep_parallel( mg_tot, nprocs, ir, id )

       !call watchb( ttmp, timg(:,2) )

       mg=ir(myrank)
       if ( allocated(LG) ) deallocate( LG )
       allocate( LG(3,mg) )

       !call watchb( ttmp, timg(:,3) )

       call construct_grid(bb,ecut,mg_tot,m1,m2,m3,mg,LG)

       !call watchb( ttmp, timg(:,4) )

       call calc_ewldg( mg, LG, ewldg )

       ewldg=ewldg*4.d0*pi/Vcell
       ewldg=ewldg-pi/(eta*Vcell)*Qtot**2
       ewldg=ewldg-2.d0*sqrt(eta/pi)*Qtot2

       !call watchb( ttmp, timg(:,5) )

       call watchb( ttmp0, timg(:,0) )

       !if ( disp_sw ) then
       !   write(*,'(1x,i4,2i8,2g26.15,4g10.3)') &
       !        loop_g,mg_tot,mg,ewldg,ecut,timg(1:2,0)
       !   !call write_watchb( timg(:,1:5), 5, (/"sg","pp","al","cg","eg"/) ) 
       !end if

       if ( abs((ewldg-ewldg_0)/ewldg) < 1.d-15 ) then
          ewldg=ewldg_0
          ecut =g_0*g_0
          exit
       else
          ewldg_0=ewldg
          g_0=g
          g=(g+max_bb)*factor
       end if

    end do ! loop_g

!--- R

    call atom_distance_minmax( xyz_min, xyz_max ) 

    alpha = sqrt(eta)

    if ( r_0 == 0.0d0 ) then
       r_0=6.0d0/sqrt(eta)
    else if ( r_0 > r_0_0*2 ) then
       r_0=6.0d0/sqrt(eta)
    end if

    r=r_0

    ewldr = ewldg
    ewldg = 0.0d0

    do loop_r=1,max_loop_r

       rrcut=r*r

       call watchb( ttmp ) ; timr=0.0d0

       call search_grid( aa, rrcut, m1,m2,m3, mr_tot )
       mr=mr_tot ; id(:)=0 ; ir(:)=mr

       if ( allocated(grd_chk) ) deallocate(grd_chk)
       allocate( grd_chk(-m1:m1,-m2:m2,-m3:m3) ) ; grd_chk=0

       if ( loop_r == 1 ) then
       else

          irange_bak = irange

          do i=1,size(LR,2)
             grd_chk( LR(1,i),LR(2,i),LR(3,i) ) = 1
          end do

       end if

       if ( allocated(LR) ) deallocate(LR)
       allocate( LR(3,mr) ) ; LR=0

       call construct_grid(aa,rrcut,mr_tot,m1,m2,m3,mr,LR)

       if ( allocated(LR_tmp) ) deallocate( LR_tmp )
       allocate( LR_tmp(3,mr) ) ; LR_tmp=0

       if ( loop_r == 1 ) then
          mr_tmp=mr
          LR_tmp(1:3,1:mr) = LR(1:3,1:mr)
       else
          mr_tmp=0
          do i=1,mr
             if ( grd_chk(LR(1,i),LR(2,i),LR(3,i)) == 0 ) then
                mr_tmp=mr_tmp+1
                LR_tmp(1:3,mr_tmp) = LR(1:3,i)
             end if
          end do
       end if

       if ( allocated(icheck_LR_tmp) ) deallocate( icheck_LR_tmp )
       allocate( icheck_LR_tmp(mr_tmp) ) ; icheck_LR_tmp=0

! ---

       call watchb( ttmp, timr(:,1) )

       call calc_ewldr( mr_tmp, LR_tmp, ewldr, icheck=icheck_LR_tmp )

       call watchb( ttmp, timr(:,2) )

! ---

       grd_chk=0
       do i=1,size(LR,2)
          grd_chk( LR(1,i),LR(2,i),LR(3,i) ) = 1
       end do
       do i=1,size(icheck_LR_tmp)
          if ( icheck_LR_tmp(i) == 0 ) grd_chk(LR_tmp(1,i),LR_tmp(2,i),LR_tmp(3,i))=0
       end do

       allocate( grd_chk0(-m1:m1,-m2:m2,-m3:m3) ) ; grd_chk0=0
       grd_chk0=grd_chk
       call mpi_allreduce(grd_chk0,grd_chk,size(grd_chk) &
            ,mpi_integer,mpi_sum,mpi_comm_world,ierr)
       deallocate( grd_chk0 )

! ---

       irange=0
       do i3=-m3,m3
       do i2=-m2,m2
       do i1=-m1,m1
          if ( grd_chk(i1,i2,i3) /= 0 ) then
             irange(1,1) = min( irange(1,1), i1 )
             irange(2,1) = max( irange(2,1), i1 )
             irange(1,2) = min( irange(1,2), i2 )
             irange(2,2) = max( irange(2,2), i2 )
             irange(1,3) = min( irange(1,3), i3 )
             irange(2,3) = max( irange(2,3), i3 )
          end if
       end do
       end do
       end do

! ---

       !if ( disp_sw ) then
       !   write(*,'(1x,i4,4i8,2g26.15,4g10.3)') &
       !        loop_r,mr_tot,mr,mr_tmp,count(icheck_LR_tmp/=0),ewldr,sqrt(rrcut),timr(1:2,1),timr(1:2,2)
       !   write(*,'(1x,6i8)') irange
       !   write(*,'(1x,6i8)') irange_bak
       !end if

       if ( abs((ewldr-ewldr_0)/ewldr) < 1.d-15 ) then
          ewldr=ewldr_0
          rrcut=r_0*r_0
          if ( r_0_0 == 0.0d0 ) r_0_0=r_0
          m1=max( abs(irange_bak(1,1)), abs(irange_bak(2,1)) )
          m2=max( abs(irange_bak(1,2)), abs(irange_bak(2,2)) )
          m3=max( abs(irange_bak(1,3)), abs(irange_bak(2,3)) )
          mr=(2*m1+1)*(2*m2+1)*(2*m3+1)
          if ( disp_sw ) write(*,*) "mr,m1,m2,m3",mr,m1,m2,m3
          deallocate( LR )
          allocate( LR(3,mr) ) ; LR=0
          i=0
          do i3=-m3,m3
          do i2=-m2,m2
          do i1=-m1,m1
             i=i+1
             LR(1,i)=i1
             LR(2,i)=i2
             LR(3,i)=i3
          end do
          end do
          end do
          exit
       else
          ewldr_0=ewldr
          r_0=r_1
          r=(r_1+max_aa*1.5d0)*factor
       end if

    end do ! loop_r

    if ( allocated(icheck_LR_tmp) ) deallocate( icheck_LR_tmp )
    if ( allocated(grd_chk) ) deallocate( grd_chk )

    Ewld=0.5d0*(ewldr+ewldg)

    call write_border( 0, " calc_ewald(end)" )

    return
  END SUBROUTINE calc_ewald


  SUBROUTINE atom_distance_minmax( rmin, rmax )
    implicit none
    real(8),intent(OUT) :: rmin(3),rmax(3)
    real(8) :: x,y,z,a1,a2,a3,tmin(3),tmax(3)
    integer :: i1,i2,i
    tmin=1.d100
    tmax=0.d0
    do i=1,mpair
       i1=ipair(1,i)
       i2=ipair(2,i)
       a1=aa_atom(1,i1)-aa_atom(1,i2)
       a2=aa_atom(2,i1)-aa_atom(2,i2)
       a3=aa_atom(3,i1)-aa_atom(3,i2)
       x=aa(1,1)*a1+aa(1,2)*a2+aa(1,3)*a3
       y=aa(2,1)*a1+aa(2,2)*a2+aa(2,3)*a3
       z=aa(3,1)*a1+aa(3,2)*a2+aa(3,3)*a3
       tmin(1)=min(tmin(1),x)
       tmin(2)=min(tmin(2),y)
       tmin(3)=min(tmin(3),z)
       tmax(1)=max(tmax(1),x)
       tmax(2)=max(tmax(2),y)
       tmax(3)=max(tmax(3),z)
    end do
    call mpi_allreduce(tmin,rmin,3,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,i)
    call mpi_allreduce(tmax,rmax,3,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,i)
  END SUBROUTINE atom_distance_minmax


END MODULE ewald_module
