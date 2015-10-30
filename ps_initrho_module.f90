MODULE ps_initrho_module

  use rgrid_module, only: dV,Ngrid,Igrid,Hgrid
  use ggrid_module
  use atom_module, only: Natom,Nelement, ki_atom,aa_atom
  use strfac_module, only: SGK
  use density_module, only: rho
  use pseudopot_module, only: Mr,rad,rab,cdd,Zps,cdd_coef
  use electron_module, only: Nspin, dspin, NSelectron
  use parallel_module
  use aa_module
  use fft_module

  implicit none

  PRIVATE
  PUBLIC :: init_ps_initrho,construct_ps_initrho,construct_r_ps_initrho

  logical :: flag_initrho_0
  logical,allocatable :: flag_initrho(:)
  real(8),allocatable :: cddg(:,:)

CONTAINS


  SUBROUTINE init_ps_initrho
    implicit none
    integer :: i,ig,ik,MKI
    real(8) :: sum0,G,x,sb,const,Vcell
    real(8),allocatable :: tmp(:)

    allocate( flag_initrho(Nelement) )
    flag_initrho(:)=.false.
    flag_initrho_0 =.false.
    if ( allocated(cdd) ) then
       do ik=1,Nelement
          if ( all(cdd(:,ik)==0.d0) ) cycle
          flag_initrho(ik) = .true.
          flag_initrho_0   = .true.
       end do
    end if

    if ( .not.flag_initrho_0 ) return

    MKI   = Nelement
    Vcell = Ngrid(0)*dV
    const = 1.d0/Vcell
    allocate( cddg(NMGL,MKI) ) ; cddg=0.d0
    do ik=1,MKI
       allocate( tmp(Mr(ik)) )
       do ig=1,NMGL
          G=sqrt(GG(ig))
          if ( G == 0.d0 ) then
             do i=1,Mr(ik)
                tmp(i)=cdd(i,ik)*rab(i,ik)
             end do
          else
             do i=1,Mr(ik)
                x=G*rad(i,ik)
                if ( x<1.d-1 ) then
                   sb=-(1.d0/39916800.d0*x**10-1.d0/362880.d0*x**8 &
                   +1.d0/5040.d0*x**6-1.d0/120.d0*x**4+1.d0/6.d0*x**2-1.d0)
                else
                   sb=sin(x)/x
                end if
                tmp(i)=cdd(i,ik)*sb*rab(i,ik)
             end do
          end if
          call simp(tmp,sum0,2)
          cddg(ig,ik)=sum0*const
       end do
       deallocate( tmp )
    end do ! ik
  END SUBROUTINE init_ps_initrho

  SUBROUTINE simp(f,s,m)
    implicit none
    integer,intent(IN)  :: m
    real(8),intent(IN)  :: f(:)
    real(8),intent(OUT) :: s
    real(8),allocatable :: g(:)
    integer :: i,n,nn,nmax
    n=size(f) ; nmax=int(n/m)*m
    do i=0,m
       nmax=nmax+i ; if ( nmax>=n ) exit
    end do
    allocate( g(nmax) ) ; g(1:n)=f ; if ( nmax>n ) g(n+1:)=0.d0
    select case(m)
    case default
       s = 0.5d0*(f(1)+f(n)) + sum(f(2:n-1))
    case(2)
       s=0.d0
       do i=1,nmax-2,2
          s = s + g(i) + 4.d0*g(i+1) + g(i+2)
       end do
       s=s/3.d0
    case(4)
       s=0.d0
       do i=1,nmax-4,4
          s=s+7*g(i)+32*g(i+1)+12*g(i+2)+32*g(i+3)+7*g(i+4)
       end do
       s=s*2.d0/45.d0
    case(6)
       s=0.d0
       do i=1,nmax-6,6
          s=s+41*g(i)+216*g(i+1)+27*g(i+2)+272*g(i+3) &
               +27*g(i+4)+216*g(i+5)+41*g(i+6)
       end do
       s=s/140.d0
    end select
    deallocate( g )
    return
  END SUBROUTINE simp


  SUBROUTINE construct_ps_initrho
    implicit none
    integer :: a,i,i1,i2,i3,ik,j,MG,ierr
    integer :: ML1,ML2,ML3,ML,ML_0,ML_1
    complex(8),allocatable :: zwork(:,:,:),zwork1(:,:,:),vg(:)
    real(8) :: c,c0
    real(8),allocatable :: rho_tmp(:,:),nelectron_ik(:)

    if ( .not. flag_initrho_0 ) return

   ! call construct_r_ps_initrho
   ! return

    MG   = NGgrid(0)
    ML   = Ngrid(0)
    ML1  = Ngrid(1)
    ML2  = Ngrid(2)
    ML3  = Ngrid(3)
    ML_0 = Igrid(1,0)
    ML_1 = Igrid(2,0)

    if ( .not. allocated(rho) ) then
       allocate( rho(ML_0:ML_1,Nspin) )
    end if
    rho=0.0d0

    allocate( zwork(0:ML1-1,0:ML2-1,0:ML3-1) )
    allocate( vg(MG) ) ; vg=(0.d0,0.d0)
    allocate( rho_tmp(ML_0:ML_1,Nelement) ) ; rho_tmp=0.d0
    allocate( nelectron_ik(Nelement) ) ; nelectron_ik=0.0d0

    do i=1,Natom
       ik=ki_atom(i)
       nelectron_ik(ik)=nelectron_ik(ik)+Zps(ik)
    end do

    call construct_Ggrid(2)

    call init_fft

    do ik=1,Nelement

    do i=MG_0,MG_1
       j=MGL(i)
       vg(i)=cddg(j,ik)*SGK(i,ik)
    end do

    call allgatherv_Ggrid(vg)

    zwork(:,:,:)=(0.d0,0.d0)
    do i=1,NGgrid(0)
       zwork(LLG(1,i),LLG(2,i),LLG(3,i))=vg(i)
    end do

    call backward_fft( zwork, zwork1 )

    i=ML_0-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       rho_tmp(i,ik)=rho_tmp(i,ik)+real( zwork(i1,i2,i3) )
    end do
    end do
    end do

    if ( minval(rho_tmp(:,ik)) < 0.d0 ) then
       write(*,*) "WARNING: rho is negative at some points",minval(rho_tmp(:,ik))
    end if
!    rho=abs(rho)
    where( rho_tmp(:,ik) < 0.d0 )
      rho_tmp(:,ik)=0.d0
    end where

    c0=sum( rho_tmp(:,ik) )*dV
    call mpi_allreduce(c0,c,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
    if ( disp_switch_parallel ) write(*,*) c,nelectron_ik(ik)

    c=nelectron_ik(ik)/c
    rho(:,1) = rho(:,1) + c*rho_tmp(:,ik)

    c0=sum( rho(:,1) )*dV
    call mpi_allreduce(c0,c,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
    if ( disp_switch_parallel ) write(*,*) "sum(rho)*dV=",c

    end do ! ik

    call finalize_fft
    call destruct_Ggrid

    deallocate( nelectron_ik )
    deallocate( vg )
    if ( allocated(zwork1) ) deallocate( zwork1 )
    deallocate( zwork )

    if ( Nspin>1 ) then
       c=1.d0/Nspin
       rho(:,1)=c*rho(:,1)
       do i=2,Nspin
          rho(:,i)=rho(:,1)
       end do
    end if

  END SUBROUTINE construct_ps_initrho


  SUBROUTINE construct_r_ps_initrho
    implicit none
    integer :: iatm,ielm,i,i1,i2,i3,j1,j2,j3,ig,mg,ML_0,ML_1,ierr
    real(8) :: a,b,c,a1,a2,a3,x,y,z,r1,r2,r3,c1,c2,c3,pi4,rr,zi
    real(8),allocatable :: rho_tmp(:)

    if ( Nspin == 2 .and. allocated(dspin) ) then
       call construct_r_spin_ps_initrho
       return
    end if

    ML_0 = Igrid(1,0)
    ML_1 = Igrid(2,0)
    mg   = size(cdd_coef,2)
    c1   = 1.0d0/Ngrid(1)
    c2   = 1.0d0/Ngrid(2)
    c3   = 1.0d0/Ngrid(3)
    pi4  = 4.0d0*acos(-1.0d0)

    if ( .not. allocated(rho) ) then
       allocate( rho(ML_0:ML_1,Nspin) )
    end if
    rho=0.0d0

    allocate( rho_tmp(ML_0:ML_1) )
    rho_tmp=0.0d0

    do ielm=1,Nelement

       rho_tmp(:)=0.0d0
       zi=0.0d0

       do iatm=1,Natom

          if ( ki_atom(iatm) /= ielm ) cycle
          zi=zi+Zps(ielm)
          a1=aa_atom(1,iatm)
          a2=aa_atom(2,iatm)
          a3=aa_atom(3,iatm)

          do j3=-1,1
          do j2=-1,1
          do j1=-1,1
             i=ML_0-1
             do i3=Igrid(1,3),Igrid(2,3)
             do i2=Igrid(1,2),Igrid(2,2)
             do i1=Igrid(1,1),Igrid(2,1)
                r1=i1*c1
                r2=i2*c2
                r3=i3*c3
                x=aa(1,1)*(r1-a1-j1)+aa(1,2)*(r2-a2-j2)+aa(1,3)*(r3-a3-j3)
                y=aa(2,1)*(r1-a1-j1)+aa(2,2)*(r2-a2-j2)+aa(2,3)*(r3-a3-j3)
                z=aa(3,1)*(r1-a1-j1)+aa(3,2)*(r2-a2-j2)+aa(3,3)*(r3-a3-j3)
                rr=x*x+y*y+z*z
                i=i+1
                do ig=1,mg
                   a=cdd_coef(1,ig,ielm)
                   b=cdd_coef(2,ig,ielm)
                   c=cdd_coef(3,ig,ielm)
                   rho_tmp(i)=rho_tmp(i)+(a+b*rr)*exp(-c*rr)
                end do
             end do
             end do
             end do
          end do
          end do
          end do

       end do ! iatm

       a=sum(rho_tmp)*dV
       call mpi_allreduce(a,r1,1,mpi_real8,mpi_sum,comm_grid,ierr)
       a=minval(rho_tmp)
       call mpi_allreduce(a,r2,1,mpi_real8,mpi_min,comm_grid,ierr)
       a=maxval(rho_tmp)
       call mpi_allreduce(a,r3,1,mpi_real8,mpi_max,comm_grid,ierr)
       if ( disp_switch_parallel ) write(*,*) zi,r1,r2,r3
       a=zi/r1
       rho_tmp(:)=a*rho_tmp(:)

       rho(:,1)=rho(:,1)+rho_tmp(:)

       a=sum(rho)*dV
       call mpi_allreduce(a,r1,1,mpi_real8,mpi_sum,comm_grid,ierr)
       a=minval(rho)
       call mpi_allreduce(a,r2,1,mpi_real8,mpi_min,comm_grid,ierr)
       a=maxval(rho)
       call mpi_allreduce(a,r3,1,mpi_real8,mpi_max,comm_grid,ierr)
       if ( disp_switch_parallel ) write(*,*) r1,r2,r3

    end do ! ielm

    deallocate( rho_tmp )

    if ( Nspin > 1 ) then
       c=1.0d0/Nspin
       rho(:,1) = c*rho(:,1)
       do i=2,Nspin
          rho(:,i) = rho(:,1)
       end do
    end if

  END SUBROUTINE construct_r_ps_initrho


  SUBROUTINE construct_r_spin_ps_initrho
    implicit none
    integer :: iatm,ielm,i,i1,i2,i3,j1,j2,j3,ig,mg,ML_0,ML_1,ierr,s
    real(8) :: a,b,c,a1,a2,a3,x,y,z,r1,r2,r3,c1,c2,c3,pi4,rr
    real(8),allocatable :: rho_tmp(:)

    if ( disp_switch_parallel ) then
       write(*,'(a60," const_r_spin_initrho")') repeat("-",60)
    end if

    ML_0 = Igrid(1,0)
    ML_1 = Igrid(2,0)
    mg   = size(cdd_coef,2)
    c1   = 1.0d0/Ngrid(1)
    c2   = 1.0d0/Ngrid(2)
    c3   = 1.0d0/Ngrid(3)
    pi4  = 4.0d0*acos(-1.0d0)

    if ( .not. allocated(rho) ) then
       allocate( rho(ML_0:ML_1,Nspin) )
    end if
    rho=0.0d0

    allocate( rho_tmp(ML_0:ML_1) )
    rho_tmp=0.0d0

    do iatm=1,Natom

       ielm=ki_atom(iatm)

       a1=aa_atom(1,iatm)
       a2=aa_atom(2,iatm)
       a3=aa_atom(3,iatm)

       rho_tmp(:)=0.0d0

       do j3=-1,1
       do j2=-1,1
       do j1=-1,1
          i=ML_0-1
          do i3=Igrid(1,3),Igrid(2,3)
          do i2=Igrid(1,2),Igrid(2,2)
          do i1=Igrid(1,1),Igrid(2,1)
             r1=i1*c1
             r2=i2*c2
             r3=i3*c3
             x=aa(1,1)*(r1-a1-j1)+aa(1,2)*(r2-a2-j2)+aa(1,3)*(r3-a3-j3)
             y=aa(2,1)*(r1-a1-j1)+aa(2,2)*(r2-a2-j2)+aa(2,3)*(r3-a3-j3)
             z=aa(3,1)*(r1-a1-j1)+aa(3,2)*(r2-a2-j2)+aa(3,3)*(r3-a3-j3)
             rr=x*x+y*y+z*z
             i=i+1
             do ig=1,mg
                a=cdd_coef(1,ig,ielm)
                b=cdd_coef(2,ig,ielm)
                c=cdd_coef(3,ig,ielm)
                rho_tmp(i)=rho_tmp(i)+(a+b*rr)*exp(-c*rr)
             end do
          end do
          end do
          end do
       end do
       end do
       end do

       c=( Zps(ielm) + dspin(iatm) )/Nspin /Zps(ielm)
       rho(:,1) = rho(:,1) + c*rho_tmp(:)

       c=( Zps(ielm) - dspin(iatm) )/Nspin /Zps(ielm)
       rho(:,2) = rho(:,2) + c*rho_tmp(:)

    end do ! iatm

    do s=1,Nspin
       a=sum(rho(:,s))*dV
       call mpi_allreduce(a,b,1,mpi_real8,mpi_sum,comm_grid,ierr)
       c=NSelectron(s)/b
       rho(:,s)=c*rho(:,s)
    end do

    do s=1,Nspin
       a=sum( rho(:,s) )*dV
       call mpi_allreduce(a,c1,1,mpi_real8,mpi_sum,comm_grid,ierr)
       b=minval( rho(:,s) )
       call mpi_allreduce(b,c2,1,mpi_real8,mpi_sum,comm_grid,ierr)
       c=maxval( rho(:,s) )
       call mpi_allreduce(c,c3,1,mpi_real8,mpi_sum,comm_grid,ierr)
       if ( disp_switch_parallel ) then
          write(*,'(1x,"sum(init_rho(",i1,"))=",2f15.10)') s,c1,NSelectron(s)
          write(*,'(1x,"min(init_rho(",i1,"))=",f15.10)') s,c2
          write(*,'(1x,"max(init_rho(",i1,"))=",f15.10)') s,c3
       end if
    end do

    deallocate( rho_tmp )

  END SUBROUTINE construct_r_spin_ps_initrho


END MODULE ps_initrho_module
