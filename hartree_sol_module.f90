MODULE hartree_sol_module

  use hartree_variables, only: E_hartree,Vh
  use rgrid_module, only: Ngrid,Igrid,dV
  use ggrid_module, only: GG,LLG,MGL,NGgrid,construct_ggrid,destruct_ggrid
  use parallel_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: calc_hartree_sol

CONTAINS

  SUBROUTINE calc_hartree_sol(n1,n2,n3,rho)
    implicit none
    integer,intent(IN) :: n1,n2,n3
    real(8),intent(IN) :: rho(n1:n2,n3)
    integer :: i,i1,i2,i3,j1,j2,j3,ierr,irank,ispin
    real(8) :: Eh0,pi4,g2,ctt(0:5),ett(0:5)
    real(8),allocatable :: work(:)
    complex(8),parameter :: z0=(0.d0,0.d0)
    complex(8),allocatable :: zwork0(:,:,:),zwork1(:,:,:)
    integer :: ML1,ML2,ML3,ML
    integer :: ifacx(30),ifacy(30),ifacz(30)
    integer,allocatable :: lx1(:),lx2(:),ly1(:),ly2(:),lz1(:),lz2(:)
    complex(8),allocatable :: wsavex(:),wsavey(:),wsavez(:)
    logical :: disp_sw

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

    call check_disp_switch( disp_sw, 0 )
    if ( disp_sw ) then
       write(*,*) "time(hatree1)=",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(hatree2)=",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(hatree3)=",ctt(3)-ctt(2),ett(3)-ett(2)
       write(*,*) "time(hatree4)=",ctt(4)-ctt(3),ett(4)-ett(3)
       write(*,*) "time(hatree5)=",ctt(5)-ctt(4),ett(5)-ett(4)
    end if

  END SUBROUTINE calc_hartree_sol

END MODULE hartree_sol_module
