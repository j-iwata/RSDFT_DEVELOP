MODULE hartree_module

  use hartree_mol_module
  use hartree_esm_module
  use bb_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: E_hartree,Vh,calc_hartree

  real(8) :: E_hartree
  real(8),allocatable :: Vh(:)
  integer :: SYStype=0

CONTAINS


  SUBROUTINE calc_hartree(n1,n2,n3,rho,SYStype_in)
    implicit none
    integer,intent(IN) :: n1,n2,n3
    real(8),intent(IN) :: rho(n1:n2,n3)
    integer,optional,intent(IN) :: SYStype_in

    if ( .not.allocated(Vh) ) then
       allocate( Vh(n1:n2) )
       Vh=0.d0
    end if

    if ( present(SYStype_in) ) then
       SYStype=SYStype_in
    end if

    select case(SYStype)
    case default
#ifdef _FFTE_
       call calc_hartree_sol_ffte(n1,n2,n3,rho)
#else
       call calc_hartree_sol(n1,n2,n3,rho)
#endif
    case(1)
       call calc_hartree_mol(n1,n2,n3,rho,Vh,E_hartree)
    case(3)
       call calc_hartree_esm(n1,n2,n3,rho,Vh,E_hartree)
    end select

  END SUBROUTINE calc_hartree


  SUBROUTINE calc_hartree_sol(n1,n2,n3,rho)
    use rgrid_module, only: Ngrid,Igrid,dV
    use ggrid_module
    use parallel_module
    implicit none
    integer,intent(IN) :: n1,n2,n3
    real(8),intent(IN) :: rho(n1:n2,n3)
    integer :: i,i1,i2,i3,j1,j2,j3,ierr,irank,ispin,a1,a2,a3,b1,b2,b3
    real(8) :: Eh0,pi4,g2,ctt(0:5),ett(0:5)
    real(8),allocatable :: work(:)
    complex(8),parameter :: z0=(0.d0,0.d0)
    complex(8),allocatable :: zwork0(:,:,:),zwork1(:,:,:)
    integer :: ML1,ML2,ML3,ML
    integer :: ifacx(30),ifacy(30),ifacz(30)
    integer,allocatable :: lx1(:),lx2(:),ly1(:),ly2(:),lz1(:),lz2(:)
    complex(8),allocatable :: wsavex(:),wsavey(:),wsavez(:)

    pi4 = 4.d0*acos(-1.d0)

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    ctt(:)=0.d0
    ett(:)=0.d0

    call watch(ctt(0),ett(0))

    allocate( zwork0(0:ML1-1,0:ML2-1,0:ML3-1) )

    allocate( work(ML) )
    work(n1:n2) = rho(n1:n2,1)
    do ispin=2,n3
       work(n1:n2) = work(n1:n2) + rho(n1:n2,ispin)
    end do

    call mpi_allgatherv(work(n1),n2-n1+1,mpi_real8 &
         ,work,ir_grid,id_grid,mpi_real8,comm_grid,ierr)

    i=0
    irank=-1
    do i3=1,node_partition(3)
    do i2=1,node_partition(2)
    do i1=1,node_partition(1)
       irank=irank+1
       do j3=pinfo_grid(5,irank),pinfo_grid(5,irank)+pinfo_grid(6,irank)-1
       do j2=pinfo_grid(3,irank),pinfo_grid(3,irank)+pinfo_grid(4,irank)-1
       do j1=pinfo_grid(1,irank),pinfo_grid(1,irank)+pinfo_grid(2,irank)-1
          i=i+1
          zwork0(j1,j2,j3)=work(i)
       end do
       end do
       end do
    end do
    end do
    end do

    deallocate( work )

    allocate( zwork1(0:ML1-1,0:ML2-1,0:ML3-1) )

    allocate( lx1(ML),lx2(ML),ly1(ML),ly2(ML),lz1(ML),lz2(ML) )
    allocate( wsavex(ML1),wsavey(ML2),wsavez(ML3) )

    call prefft(ML1,ML2,ML3,ML,wsavex,wsavey,wsavez &
         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call watch(ctt(1),ett(1))

    call fft3fx(ML1,ML2,ML3,ML,zwork0,zwork1,wsavex,wsavey,wsavez &
         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call watch(ctt(2),ett(2))

    call construct_Ggrid(2)

    zwork1(:,:,:)=z0
    do i=1,NGgrid(0)
       g2=GG(MGL(i))
       if ( g2 == 0.d0 ) cycle
       i1=LLG(1,i)
       i2=LLG(2,i)
       i3=LLG(3,i)
       zwork1(i1,i2,i3)=zwork0(i1,i2,i3)*pi4/g2
    end do

    call destruct_Ggrid

    call watch(ctt(3),ett(3))

    call fft3bx(ML1,ML2,ML3,ML,zwork1,zwork0,wsavex,wsavey,wsavez &
         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call watch(ctt(4),ett(4))

    deallocate( wsavez,wsavey,wsavex )
    deallocate( lz2,lz1,ly2,ly1,lx2,lx1 )
    deallocate( zwork0 )

    i=n1-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       Vh(i)=real( zwork1(i1,i2,i3) )
    end do
    end do
    end do

    Eh0=0.d0
    do ispin=1,n3
       i=n1-1
       do i3=Igrid(1,3),Igrid(2,3)
       do i2=Igrid(1,2),Igrid(2,2)
       do i1=Igrid(1,1),Igrid(2,1)
          i=i+1
          Eh0 = Eh0 + real( zwork1(i1,i2,i3) )*rho(i,ispin)
       end do
       end do
       end do
    end do
    Eh0=0.5d0*Eh0*dV
    call mpi_allreduce(Eh0,E_hartree,1,mpi_real8,mpi_sum,comm_grid,ierr)

    deallocate( zwork1 )

    call watch(ctt(5),ett(5))

    if ( disp_switch_parallel ) then
       write(*,*) "time(hatree1)=",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(hatree2)=",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(hatree3)=",ctt(3)-ctt(2),ett(3)-ett(2)
       write(*,*) "time(hatree4)=",ctt(4)-ctt(3),ett(4)-ett(3)
       write(*,*) "time(hatree5)=",ctt(5)-ctt(4),ett(5)-ett(4)
    end if

  END SUBROUTINE calc_hartree_sol


  SUBROUTINE calc_hartree_sol_ffte(n1,n2,n3,rho)
    use rgrid_module, only: Ngrid,Igrid,dV
    use ggrid_module
    use parallel_module
    implicit none
    integer,intent(IN) :: n1,n2,n3
    real(8),intent(IN) :: rho(n1:n2,n3)
    integer :: i,i1,i2,i3,j1,j2,j3,ierr,irank,ispin,a1,a2,a3,b1,b2,b3
    real(8) :: Eh0,pi4,g2,ctt(0:5),ett(0:5)
    real(8),allocatable :: work(:)
    complex(8),parameter :: z0=(0.d0,0.d0)
    complex(8),allocatable :: zwork0(:,:,:),zwork1(:,:,:),zwork2(:,:,:)
    integer :: ML1,ML2,ML3,ML
    integer :: ifacx(30),ifacy(30),ifacz(30)
    integer :: comm_fftx, npux, ix
    integer :: comm_ffty, npuy, iy
    integer :: comm_fftz, npuz, iz
    integer :: icolor,j,m
    integer :: MG,ML_0,ML_1,a1b,b1b,a2b,b2b,a3b,b3b,np1,np2,np3
    integer :: MG1,MG2,MG3,NG1,NG2,NG3,MGmax


    pi4 = 4.d0*acos(-1.d0)

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    MG    = NGgrid(0)
    ML_0  = Igrid(1,0)
    ML_1  = Igrid(2,0)
    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    np1 = node_partition(1)
    np2 = node_partition(2)
    np3 = node_partition(3)
    MG1 = -ML1/2    ; MG2=-ML2/2    ; MG3=-ML3/2
    NG1 = (ML1-1)/2 ; NG2=(ML2-1)/2 ; NG3=(ML3-1)/2

    ix=a1b/(ML1/np1); iy=a2b/(ML2/np2); iz=a3b/(ML3/np3);

    icolor=iy+iz*np2
    call mpi_comm_split(comm_grid,icolor, 0, comm_fftx, ierr)
    icolor=iz+ix*nprocs;
    call mpi_comm_split(comm_grid,icolor, 0, comm_ffty, ierr)
    icolor=iy+ix*nprocs;
    call mpi_comm_split(comm_grid,icolor, 0, comm_fftz, ierr)

    call mpi_comm_size(comm_fftx, npux, ierr)
    call mpi_comm_size(comm_ffty, npuy, ierr)
    call mpi_comm_size(comm_fftz, npuz, ierr)
    
    ctt(:)=0.d0
    ett(:)=0.d0

    call watch(ctt(0),ett(0))

    allocate( zwork0(0:ML1-1,0:ML2-1,0:ML3-1) )
    allocate( zwork1(0:ML1-1,a2b:b2b,a3b:b3b),zwork2(0:ML1-1,a2b:b2b,a3b:b3b))
    zwork1(:,:,:)=z0
    zwork2(:,:,:)=z0


    allocate( work(ML) )
    work(n1:n2) = rho(n1:n2,1)
    do ispin=2,n3
       work(n1:n2) = work(n1:n2) + rho(n1:n2,ispin)
    end do

    call mpi_allgatherv(work(n1),n2-n1+1,mpi_real8 &
         ,work,ir_grid,id_grid,mpi_real8,comm_grid,ierr)

    i=0
    irank=-1
    do i3=1,node_partition(3)
    do i2=1,node_partition(2)
    do i1=1,node_partition(1)
       irank=irank+1
       do j3=pinfo_grid(5,irank),pinfo_grid(5,irank)+pinfo_grid(6,irank)-1
       do j2=pinfo_grid(3,irank),pinfo_grid(3,irank)+pinfo_grid(4,irank)-1
       do j1=pinfo_grid(1,irank),pinfo_grid(1,irank)+pinfo_grid(2,irank)-1
          i=i+1
          zwork0(j1,j2,j3)=work(i)
       end do
       end do
       end do
    end do
    end do
    end do

    deallocate( work )

    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=a1b,b1b
       zwork1(i1,i2,i3)=zwork0(i1,i2,i3)
    end do
    end do
    end do

    call mpi_allreduce(zwork1,zwork2,ML1*(b2b-a2b+1)*(b3b-a3b+1),mpi_complex16,mpi_sum,comm_fftx,ierr)


!    allocate( zwork1(0:ML1-1,0:ML2-1,0:ML3-1) )
!    allocate( lx1(ML),lx2(ML),ly1(ML),ly2(ML),lz1(ML),lz2(ML) )
!    allocate( wsavex(ML1),wsavey(ML2),wsavez(ML3) )
!    call prefft(ML1,ML2,ML3,ML,wsavex,wsavey,wsavez &
!         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call watch(ctt(1),ett(1))


    call pzfft3dv(zwork2,zwork1,ML1,ML2,ML3,comm_ffty,comm_fftz,npuy,npuz,0)
    call pzfft3dv(zwork2,zwork1,ML1,ML2,ML3,comm_ffty,comm_fftz,npuy,npuz,-1)


!    call fft3fx(ML1,ML2,ML3,ML,zwork0,zwork1,wsavex,wsavey,wsavez &
!         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call watch(ctt(2),ett(2))

    call construct_Ggrid(2)

    m=0
    do i3=MG3,NG3
    do i2=MG2,NG2
    do i1=MG1,NG1
       if ( i1==0 .and. i2==0 .and. i3==0 ) cycle
         m=m+1
    end do
    end do
    end do

    MGmax = m

    zwork2(:,:,:)=z0
    do i3=MG3,NG3
    do i2=MG2,NG2
    do i1=MG1,NG1
       g2=(bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3)**2 &
  &          +(bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3)**2 &
  &          +(bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3)**2
       if ( g2<=0.d0 .or. (MGmax/=MG .and. Ecut-1.d-14<=g2) ) cycle
       j1=mod(i1+ML1,ML1)
       j2=mod(i2+ML2,ML2)
       j3=mod(i3+ML3,ML3)
       if(a2b <= j2 .and. j2 <= b2b .and. a3b <= j3 .and. j3 <= b3b) then
         zwork2(j1,j2,j3)=zwork1(j1,j2,j3)*pi4/g2
       endif
    end do
    end do
    end do


!    zwork1(:,:,:)=z0
!    do i=1,NGgrid(0)
!       g2=GG(MGL(i))
!       if ( g2 == 0.d0 ) cycle
!       i1=LLG(1,i)
!       i2=LLG(2,i)
!       i3=LLG(3,i)
!       zwork1(i1,i2,i3)=zwork0(i1,i2,i3)*pi4/g2
!    end do

    call destruct_Ggrid

    call watch(ctt(3),ett(3))

!   call fft3bx(ML1,ML2,ML3,ML,zwork1,zwork0,wsavex,wsavey,wsavez &
!        ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call pzfft3dv(zwork2,zwork1,ML1,ML2,ML3,comm_ffty,comm_fftz,npuy,npuz,1)

    call watch(ctt(4),ett(4))

!    deallocate( wsavez,wsavey,wsavex )
!    deallocate( lz2,lz1,ly2,ly1,lx2,lx1 )
    deallocate( zwork0 )

    i=n1-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       Vh(i)=real( zwork1(i1,i2,i3) )
    end do
    end do
    end do

    Eh0=0.d0
    do ispin=1,n3
       i=n1-1
       do i3=Igrid(1,3),Igrid(2,3)
       do i2=Igrid(1,2),Igrid(2,2)
       do i1=Igrid(1,1),Igrid(2,1)
          i=i+1
          Eh0 = Eh0 + real( zwork1(i1,i2,i3) )*rho(i,ispin)
       end do
       end do
       end do
    end do
    Eh0=0.5d0*Eh0*dV
    call mpi_allreduce(Eh0,E_hartree,1,mpi_real8,mpi_sum,comm_grid,ierr)

!    deallocate( zwork1 )
    deallocate( zwork1,zwork2 )
    call mpi_comm_free(comm_fftx,ierr)
    call mpi_comm_free(comm_ffty,ierr)
    call mpi_comm_free(comm_fftz,ierr)

    call watch(ctt(5),ett(5))

    if ( disp_switch_parallel ) then
       write(*,*) "time(hatree1)=",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(hatree2)=",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(hatree3)=",ctt(3)-ctt(2),ett(3)-ett(2)
       write(*,*) "time(hatree4)=",ctt(4)-ctt(3),ett(4)-ett(3)
       write(*,*) "time(hatree5)=",ctt(5)-ctt(4),ett(5)-ett(4)
    end if

  END SUBROUTINE calc_hartree_sol_ffte


END MODULE hartree_module
