MODULE cg_u_module

  use hamiltonian_module
  use cgpc_module
  use localpot2_Smatrix_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: cg_u, init_cg_u

  integer :: ML_0,ML_1
  integer :: MB_0,MB_1
  integer :: MB_d
  real(8) :: dV
  integer :: comm_grid

  include 'mpif.h'

CONTAINS


  SUBROUTINE init_cg_u( n1,n2, m1,m2, dV_in, MB_d_in, comm_in )
    implicit none
    integer,intent(IN) :: n1,n2,m1,m2,MB_d_in,comm_in
    real(8),intent(IN) :: dV_in
    ML_0 = n1
    ML_1 = n2
    MB_0 = m1
    MB_1 = m2
    dV   = dV_in
    MB_d = MB_d_in
    comm_grid = comm_in
  END SUBROUTINE init_cg_u

#ifdef _DRSDFT_

  SUBROUTINE cg_u( k,s,Mcg,igs,unk,esp,res,disp_switch )
    implicit none
    integer,intent(IN) :: k,s,Mcg,igs
    real(8),intent(INOUT) :: unk(ML_0:,:)
    real(8),intent(INOUT) :: esp(:),res(:)
    logical,intent(IN) :: disp_switch
    integer :: ns,ne,nn,n,m,icg,ML0,Nhpsi,Npc,Ncgtot,ierr
    integer :: mm,icmp,i,TYPE_MAIN,n1,n2
    real(8),parameter :: ep0=0.d0
    real(8),parameter :: ep1=1.d-15
    real(8) :: rwork(9),W(2),c,d,r,c1,ct0,ct1,et0,et1,ctt(4),ett(4)
    real(8),allocatable :: sb(:),rb(:),E(:),E1(:),gkgk(:),bk(:)
    complex(8) :: work(9),zphase,ztmp

    real(8),parameter :: zero=0.d0
    real(8),allocatable :: hxk(:,:),hpk(:,:),gk(:,:),Pgk(:,:)
    real(8),allocatable :: pk(:,:),pko(:,:),Sf(:,:)
    real(8),allocatable :: vtmp2(:,:),wtmp2(:,:)
    real(8),allocatable :: utmp2(:,:),btmp2(:,:)
    real(8),allocatable :: uSu(:)
    real(8),allocatable :: unk_tmp(:,:)

    TYPE_MAIN = MPI_REAL8

    ctt(:)=0.d0
    ett(:)=0.d0
    ctt_hamil(:)=0.d0
    ett_hamil(:)=0.d0

    ML0 = ML_1-ML_0+1

    mm  = ML0  ; if (TYPE_MAIN==mpi_complex16) mm=2*ML0
    c1  = 2.d0 ; if (TYPE_MAIN==mpi_complex16) c1=1.d0
    icmp= 1    ; if (TYPE_MAIN==mpi_complex16) icmp=2

    Ncgtot = 0
    Nhpsi  = 0
    Npc    = 0

    n1 = ML_0
    n2 = ML_1

    allocate( hxk(n1:n2,MB_d), hpk(n1:n2,MB_d) )
    allocate( gk(n1:n2,MB_d) , Pgk(n1:n2,MB_d) )
    allocate( pk(n1:n2,MB_d) , pko(n1:n2,MB_d) )
    allocate( sb(MB_d),rb(MB_d) )
    allocate( E(MB_d),E1(MB_d),gkgk(MB_d),bk(MB_d) )
    allocate( vtmp2(6,MB_d),wtmp2(6,MB_d) )
    allocate( utmp2(2,2),btmp2(2,2) )
    allocate( Sf(n1:n2,MB_d) )
    allocate( uSu(MB_d) )
    allocate( unk_tmp(n1:n2,MB_d) )

!$OMP parallel workshare
    res(:)  = 0.d0
    esp(:)  = 0.d0
!$OMP end parallel workshare

    do ns=MB_0,MB_1
       ne=ns
       nn=ne-ns+1

       E1(:)=1.d10

       call watch(ct0,et0)

       call hamiltonian(k,s,unk(:,ns:ne),hxk,n1,n2,ns,ne) ; Nhpsi=Nhpsi+1

       call watch(ct1,et1) ; ctt(1)=ctt(1)+ct1-ct0 ; ett(1)=ett(1)+et1-et0

       do n=1,nn
          call dot_product(unk(n1,n+ns-1),hxk(n1,n),sb(n),dV,mm,1)
       end do

       call watch(ct0,et0) ; ctt(2)=ctt(2)+ct0-ct1 ; ett(2)=ett(2)+et0-et1

       call mpi_allreduce(sb,E,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       call watch(ct1,et1) ; ctt(3)=ctt(3)+ct1-ct0 ; ett(3)=ett(3)+et1-et0

! ---

       do n=1,nn
          call op_localpot2_Smatrix( unk(:,n+ns-1), Sf(:,n) )
       end do

       call watch(ct0,et0) ; ctt(1)=ctt(1)+ct0-ct1 ; ett(1)=ett(1)+et0-et1

       do n=1,nn
          call dot_product(unk(n1,n+ns-1),Sf(n1,n),sb(n),dV,mm,1)
       end do

       call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

       call mpi_allreduce(sb,uSu,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       call watch(ct0,et0) ; ctt(3)=ctt(3)+ct0-ct1 ; ett(3)=ett(3)+et0-et1

! ---

       call watch(ct1,et1)

       do n=1,nn
!$OMP parallel do
          do i=n1,n2
             gk(i,n) = -c1*( hxk(i,n) - E(n)*Sf(i,n) )
          end do
!$OMP end parallel do
       end do

       do n=1,nn
          call op_localpot2_Smatrix( gk(:,n), Sf(:,n) )
          call dot_product(gk(n1,n),Sf(n1,n),sb(n),dV,mm,1)
       end do

       call watch(ct0,et0) ; ctt(2)=ctt(2)+ct0-ct1 ; ett(2)=ett(2)+et0-et1

       call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       call watch(ct1,et1) ; ctt(3)=ctt(3)+ct1-ct0 ; ett(3)=ett(3)+et1-et0

       do icg=1,Mcg+1

          call watch(ct0,et0)

          Ncgtot=Ncgtot+1

          do n=1,nn
!$OMP parallel workshare
             Pgk(n1:n2,n)=gk(n1:n2,n)
!$OMP end parallel workshare
          end do

          res(ns:ne)=rb(1:nn)/c1**2

! --- Convergence check ---

          if ( all(rb(1:nn)<ep0) ) exit
          if ( all(abs(E(1:nn)-E1(1:nn))<ep1) ) exit
          if ( icg==Mcg+1 ) exit

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

! --- Preconditioning ---

          call preconditioning(E,k,s,nn,ML0,unk(:,ns:ne),gk,Pgk)

          call watch(ct0,et0) ; ctt(4)=ctt(4)+ct0-ct1 ; ett(4)=ett(4)+et0-et1

! ---

          do n=1,nn
             call op_localpot2_Smatrix( Pgk(:,n), Sf(:,n) )
             call dot_product(Pgk(n1,n),Pgk(n1,n),sb(n),dV,mm,1)
!             call dot_product(Pgk(n1,n),gk(n1,n),sb(n),dV,mm,1)
          end do

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

          call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

          call watch(ct0,et0) ; ctt(3)=ctt(3)+ct0-ct1 ; ett(3)=ett(3)+et0-et1

          if ( icg==1 ) then
!$OMP parallel workshare
             pk(n1:n2,1:nn) = Pgk(n1:n2,1:nn)
!$OMP end parallel workshare
          else
             do n=1,nn
                bk(n)=rb(n)/gkgk(n)
!$OMP parallel do
                do i=n1,n2
                   pk(i,n)=Pgk(i,n)+bk(n)*pk(i,n)
                end do
!$OMP end parallel do
             end do
          end if
          gkgk(1:nn)=rb(1:nn)

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

          call hamiltonian(k,s,pk,hpk,n1,n2,ns,ne) ; Nhpsi=Nhpsi+1

          call watch(ct0,et0) ; ctt(1)=ctt(1)+ct0-ct1 ; ett(1)=ett(1)+et0-et1

          do n=1,nn
             call op_localpot2_Smatrix( unk(:,n+ns-1), Sf(:,n) )
          end do

          do n=1,nn
             vtmp2(1:6,n)=zero
          end do

          do n=1,nn
             m=n+ns-1
             call dot_product(unk(n1,m),Sf(n1,n),vtmp2(1,n),dV,mm,1)
             call dot_product(pk(n1,n),Sf(n1,n),vtmp2(2,n),dV,mm,icmp)
             call dot_product(unk(n1,m),hxk(n1,n),vtmp2(4,n),dV,mm,1)
             call dot_product(pk(n1,n),hxk(n1,n),vtmp2(5,n),dV,mm,icmp)
             call dot_product(pk(n1,n),hpk(n1,n),vtmp2(6,n),dV,mm,1)
          end do

          do n=1,nn
             call op_localpot2_Smatrix( pk(:,n), Sf(:,n) )
          end do
          do n=1,nn
             call dot_product(pk(n1,n),Sf(n1,n),vtmp2(3,n),dV,mm,1)
          end do

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

          call mpi_allreduce(vtmp2,wtmp2,6*nn,TYPE_MAIN,mpi_sum,comm_grid,ierr)

          call watch(ct0,et0) ; ctt(3)=ctt(3)+ct0-ct1 ; ett(3)=ett(3)+et0-et1

          do n=1,nn
             m=n+ns-1
             btmp2(1,1)=wtmp2(1,n)
             btmp2(2,1)=wtmp2(2,n)
             btmp2(1,2)=wtmp2(2,n)
             btmp2(2,2)=wtmp2(3,n)
             utmp2(1,1)=wtmp2(4,n)
             utmp2(2,1)=wtmp2(5,n)
             utmp2(1,2)=wtmp2(5,n)
             utmp2(2,2)=wtmp2(6,n)
             if (TYPE_MAIN==mpi_complex16) then
                ztmp=btmp2(1,2)
                ztmp=conjg(ztmp)
                btmp2(1,2)=ztmp
                ztmp=utmp2(1,2)
                ztmp=conjg(ztmp)
                utmp2(1,2)=ztmp
                call zhegv(1,'V','U',2,utmp2,2,btmp2,2,W,work,9,rwork,ierr)
                if ( abs(W(1)-E(n))>1.d-1 .and. abs(W(2)-E(n))<=1.d-1 ) then
                   utmp2(1,1)=utmp2(1,2)
                   utmp2(2,1)=utmp2(2,2)
                   W(1)=W(2)
                end if
!- Fix the phase -
                ztmp=utmp2(1,1)
                r=abs(ztmp)
                c=real(ztmp)/r
                d=aimag(ztmp)/r
                zphase=dcmplx(c,-d)
                utmp2(1,1)=utmp2(1,1)*zphase
                utmp2(2,1)=utmp2(2,1)*zphase
             else
                call dsygv(1,'V','U',2,utmp2,2,btmp2,2,W,rwork,9,ierr)
                if ( abs(W(1)-E(n))>1.d-1 .and. abs(W(2)-E(n))<=1.d-1 ) then
                   utmp2(1,1)=utmp2(1,2)
                   utmp2(2,1)=utmp2(2,2)
                   W(1)=W(2)
                end if
!- Fix the phase -
                c=utmp2(1,1)
                if( c<0.d0 ) then
                   utmp2(1,1)=-utmp2(1,1)
                   utmp2(2,1)=-utmp2(2,1)
                end if
             end if

             E1(n)=E(n)
             E(n) =W(1)

             do i=n1,n2
                unk_tmp(i,n) = utmp2(1,1)*unk(i,m) + utmp2(2,1)*pk(i,n)
             end do

             call op_localpot2_Smatrix( unk_tmp(:,n), Sf(:,n) )

!$OMP parallel
!$OMP do
             do i=n1,n2
                hxk(i,n)=utmp2(1,1)*hxk(i,n)+utmp2(2,1)*hpk(i,n)
             end do
!$OMP end do
!$OMP do
             do i=n1,n2
                gk(i,n) = -c1*( hxk(i,n) - W(1)*Sf(i,n) )
             end do
!$OMP end do
!$OMP end parallel

             call op_localpot2_Smatrix( gk(:,n), Sf(:,n) )

             call dot_product(gk(n1,n),Sf(n1,n),sb(n),dV,mm,1)

          end do ! n

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

          call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

          call watch(ct0,et0) ; ctt(3)=ctt(3)+ct0-ct1 ; ett(3)=ett(3)+et0-et1

          do n=1,nn
             m=n+ns-1
             if ( rb(n)/res(m)>1.d8 ) then
                E(n)=E1(n)
                cycle
             end if
!$OMP parallel do
             do i=n1,n2
                unk(i,m) = utmp2(1,1)*unk(i,m) + utmp2(2,1)*pk(i,n)
             end do
!$OMP end parallel do
          end do

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

       end do ! icg

       esp(ns:ne)=E(1:nn)

    end do  ! band-loop

    deallocate( unk_tmp )
    deallocate( uSu )
    deallocate( Sf )
    deallocate( btmp2,utmp2  )
    deallocate( wtmp2,vtmp2  )
    deallocate( bk,gkgk,E1,E )
    deallocate( rb,sb )
    deallocate( pko,pk  )
    deallocate( Pgk,gk  )
    deallocate( hpk,hxk )

    if ( disp_switch ) then
       write(*,*) "time(hamil_kin)",ctt_hamil(1),ett_hamil(1)
       write(*,*) "time(hamil_loc)",ctt_hamil(2),ett_hamil(2)
       write(*,*) "time(hamil_nlc)",ctt_hamil(3),ett_hamil(3)
       write(*,*) "time(hamil_exx)",ctt_hamil(4),ett_hamil(4)
       write(*,*) "time(hamil_cg)",ctt(1),ett(1)
       write(*,*) "time(op_cg_u )",ctt(2),ett(2)
       write(*,*) "time(com_cg_u)",ctt(3),ett(3)
       write(*,*) "time(pc_cg_u )",ctt(4),ett(4)
    end if

  END SUBROUTINE cg_u

#else

  SUBROUTINE cg_u( k,s,Mcg,igs,unk,esp,res,disp_switch )
    implicit none
    integer,intent(IN) :: k,s,Mcg,igs
    complex(8),intent(INOUT) :: unk(ML_0:,:)
    real(8),intent(INOUT) :: esp(:),res(:)
    logical,intent(IN) :: disp_switch
    integer :: ns,ne,nn,n,m,icg,ML0,Nhpsi,Npc,Ncgtot,ierr
    integer :: mm,icmp,i,TYPE_MAIN,n1,n2
    real(8),parameter :: ep0=0.d0
    real(8),parameter :: ep1=1.d-15
    real(8) :: rwork(9),W(2),c,d,r,c1,ct0,ct1,et0,et1,ctt(4),ett(4)
    real(8),allocatable :: sb(:),rb(:),E(:),E1(:),gkgk(:),bk(:)
    complex(8) :: work(9),zphase,ztmp

    complex(8),parameter :: zero=(0.d0,0.d0)
    complex(8),allocatable :: hxk(:,:),hpk(:,:),gk(:,:),Pgk(:,:)
    complex(8),allocatable :: pk(:,:),pko(:,:),Sf(:,:)
    complex(8),allocatable :: vtmp2(:,:),wtmp2(:,:)
    complex(8),allocatable :: utmp2(:,:),btmp2(:,:)
    complex(8),allocatable :: uSu(:)
    complex(8),allocatable :: unk_tmp(:,:)

    TYPE_MAIN = MPI_COMPLEX16

    ctt(:)=0.d0
    ett(:)=0.d0
    ctt_hamil(:)=0.d0
    ett_hamil(:)=0.d0

    ML0 = ML_1-ML_0+1

    mm  = ML0  ; if (TYPE_MAIN==mpi_complex16) mm=2*ML0
    c1  = 2.d0 ; if (TYPE_MAIN==mpi_complex16) c1=1.d0
    icmp= 1    ; if (TYPE_MAIN==mpi_complex16) icmp=2

    Ncgtot = 0
    Nhpsi  = 0
    Npc    = 0

    n1 = ML_0
    n2 = ML_1

    allocate( hxk(n1:n2,MB_d), hpk(n1:n2,MB_d) )
    allocate( gk(n1:n2,MB_d) , Pgk(n1:n2,MB_d) )
    allocate( pk(n1:n2,MB_d) , pko(n1:n2,MB_d) )
    allocate( sb(MB_d),rb(MB_d) )
    allocate( E(MB_d),E1(MB_d),gkgk(MB_d),bk(MB_d) )
    allocate( vtmp2(6,MB_d),wtmp2(6,MB_d) )
    allocate( utmp2(2,2),btmp2(2,2) )
    allocate( Sf(n1:n2,MB_d) )
    allocate( uSu(MB_d) )
    allocate( unk_tmp(n1:n2,MB_d) )

!$OMP parallel workshare
    res(:)  = 0.d0
    esp(:)  = 0.d0
!$OMP end parallel workshare

    do ns=MB_0,MB_1
       ne=ns
       nn=ne-ns+1

       E1(:)=1.d10

! ---

       call watch(ct0,et0)

       call hamiltonian(k,s,unk(:,ns:ne),hxk,n1,n2,ns,ne) ; Nhpsi=Nhpsi+1

       call watch(ct1,et1) ; ctt(1)=ctt(1)+ct1-ct0 ; ett(1)=ett(1)+et1-et0

       do n=1,nn
          call dot_product(unk(n1,n+ns-1),hxk(n1,n),sb(n),dV,mm,1)
       end do

       call watch(ct0,et0) ; ctt(2)=ctt(2)+ct0-ct1 ; ett(2)=ett(2)+et0-et1

       call mpi_allreduce(sb,E,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       call watch(ct1,et1) ; ctt(3)=ctt(3)+ct1-ct0 ; ett(3)=ett(3)+et1-et0

! ---

       do n=1,nn
          call op_localpot2_Smatrix( unk(:,n+ns-1), Sf(:,n) )
       end do

       call watch(ct0,et0) ; ctt(1)=ctt(1)+ct0-ct1 ; ett(1)=ett(1)+et0-et1

       do n=1,nn
          call dot_product(unk(n1,n+ns-1),Sf(n1,n),sb(n),dV,mm,1)
       end do

       call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

       call mpi_allreduce(sb,uSu,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       call watch(ct0,et0) ; ctt(3)=ctt(3)+ct0-ct1 ; ett(3)=ett(3)+et0-et1

! ---

       do n=1,nn
!$OMP parallel do
          do i=n1,n2
             gk(i,n) = -c1*( hxk(i,n) - E(n)*Sf(i,n) )
          end do
!$OMP end parallel do
       end do

       do n=1,nn
          call op_localpot2( gk(:,n), Sf(:,n) )
          call dot_product(gk(n1,n),Sf(n1,n),sb(n),dV,mm,1)
       end do

       call watch(ct0,et0) ; ctt(2)=ctt(2)+ct0-ct1 ; ett(2)=ett(2)+et0-et1

       call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       call watch(ct1,et1) ; ctt(3)=ctt(3)+ct1-ct0 ; ett(3)=ett(3)+et1-et0

       do icg=1,Mcg+1

          Ncgtot=Ncgtot+1

          do n=1,nn
!$OMP parallel workshare
             Pgk(n1:n2,n)=gk(n1:n2,n)
!$OMP end parallel workshare
          end do

          res(ns:ne)=rb(1:nn)/c1**2

! --- Convergence check ---

          if ( all(rb(1:nn)<ep0) ) exit
          if ( all(abs(E(1:nn)-E1(1:nn))<ep1) ) exit
          if ( icg==Mcg+1 ) exit

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

! --- Preconditioning ---

          call preconditioning(E,k,s,nn,ML0,unk(n1,ns),gk,Pgk)

          call watch(ct0,et0) ; ctt(4)=ctt(4)+ct0-ct1 ; ett(4)=ett(4)+et0-et1

! ---

          do n=1,nn
             call op_localpot2_Smatrix( Pgk(:,n), Sf(:,n) )
             call dot_product(Pgk(n1,n),Sf(n1,n),sb(n),dV,mm,1)
!             call dot_product(Pgk(n1,n),gk(n1,n),sb(n),dV,mm,1)
          end do

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

          call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

          call watch(ct0,et0) ; ctt(3)=ctt(3)+ct0-ct1 ; ett(3)=ett(3)+et0-et1

          if ( icg==1 ) then
!$OMP parallel workshare
             pk(n1:n2,1:nn) = Pgk(n1:n2,1:nn)
!$OMP end parallel workshare
          else
             do n=1,nn
                bk(n)=rb(n)/gkgk(n)
!$OMP parallel do
                do i=n1,n2
                   pk(i,n)=Pgk(i,n)+bk(n)*pk(i,n)
                end do
!$OMP end parallel do
             end do
          end if
          gkgk(1:nn)=rb(1:nn)

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

          call hamiltonian(k,s,pk,hpk,n1,n2,ns,ne) ; Nhpsi=Nhpsi+1

          call watch(ct0,et0) ; ctt(1)=ctt(1)+ct0-ct1 ; ett(1)=ett(1)+et0-et1

          do n=1,nn
             call op_localpot2_Smatrix( unk(:,n+ns-1), Sf(:,n) )
          end do

          call watch(ct0,et0) ; ctt(1)=ctt(1)+ct0-ct1 ; ett(1)=ett(1)+et0-et1

          do n=1,nn
             vtmp2(1:6,n)=zero
          end do

          do n=1,nn
             m=n+ns-1
             call dot_product(unk(n1,m),Sf(n1,n),vtmp2(1,n),dV,mm,1)
             call dot_product(pk(n1,n),Sf(n1,n),vtmp2(2,n),dV,mm,icmp)
             call dot_product(unk(n1,m),hxk(n1,n),vtmp2(4,n),dV,mm,1)
             call dot_product(pk(n1,n),hxk(n1,n),vtmp2(5,n),dV,mm,icmp)
             call dot_product(pk(n1,n),hpk(n1,n),vtmp2(6,n),dV,mm,1)
          end do

          do n=1,nn
             call op_localpot2_Smatrix( pk(:,n), Sf(:,n) )
          end do
          do n=1,nn
             call dot_product(pk(n1,n),Sf(n1,n),vtmp2(3,n),dV,mm,1)
          end do

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

          call mpi_allreduce(vtmp2,wtmp2,6*nn,TYPE_MAIN,mpi_sum,comm_grid,ierr)

          call watch(ct0,et0) ; ctt(3)=ctt(3)+ct0-ct1 ; ett(3)=ett(3)+et0-et1

          do n=1,nn
             m=n+ns-1
             btmp2(1,1)=wtmp2(1,n)
             btmp2(2,1)=wtmp2(2,n)
             btmp2(1,2)=wtmp2(2,n)
             btmp2(2,2)=wtmp2(3,n)
             utmp2(1,1)=wtmp2(4,n)
             utmp2(2,1)=wtmp2(5,n)
             utmp2(1,2)=wtmp2(5,n)
             utmp2(2,2)=wtmp2(6,n)
             if (TYPE_MAIN==mpi_complex16) then
                ztmp=btmp2(1,2)
                ztmp=conjg(ztmp)
                btmp2(1,2)=ztmp
                ztmp=utmp2(1,2)
                ztmp=conjg(ztmp)
                utmp2(1,2)=ztmp
                call zhegv(1,'V','U',2,utmp2,2,btmp2,2,W,work,9,rwork,ierr)
                if ( abs(W(1)-E(n))>1.d-1 .and. abs(W(2)-E(n))<=1.d-1 ) then
                   utmp2(1,1)=utmp2(1,2)
                   utmp2(2,1)=utmp2(2,2)
                   W(1)=W(2)
                end if
!- Fix the phase -
                ztmp=utmp2(1,1)
                r=abs(ztmp)
                c=real(ztmp)/r
                d=aimag(ztmp)/r
                zphase=dcmplx(c,-d)
                utmp2(1,1)=utmp2(1,1)*zphase
                utmp2(2,1)=utmp2(2,1)*zphase
             else
                call dsygv(1,'V','U',2,utmp2,2,btmp2,2,W,rwork,9,ierr)
                if ( abs(W(1)-E(n))>1.d-1 .and. abs(W(2)-E(n))<=1.d-1 ) then
                   utmp2(1,1)=utmp2(1,2)
                   utmp2(2,1)=utmp2(2,2)
                   W(1)=W(2)
                end if
!- Fix the phase -
                c=utmp2(1,1)
                if( c<0.d0 ) then
                   utmp2(1,1)=-utmp2(1,1)
                   utmp2(2,1)=-utmp2(2,1)
                end if
             end if

             E1(n)=E(n)
             E(n) =W(1)

             do i=n1,n2
                unk_tmp(i,n) = utmp2(1,1)*unk(i,m) + utmp2(2,1)*pk(i,n)
             end do

             call op_localpot2_Smatrix( unk_tmp(:,n), Sf(:,n) )

!$OMP parallel
!$OMP do
             do i=n1,n2
                hxk(i,n)=utmp2(1,1)*hxk(i,n)+utmp2(2,1)*hpk(i,n)
             end do
!$OMP end do
!$OMP do
             do i=n1,n2
                gk(i,n) = -c1*( hxk(i,n) - W(1)*Sf(i,n) )
             end do
!$OMP end do
!$OMP end parallel

             call op_localpot2_Smatrix( gk(:,n), Sf(:,n) )

             call dot_product(gk(n1,n),Sf(n1,n),sb(n),dV,mm,1)

          end do

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

          call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

          call watch(ct0,et0) ; ctt(3)=ctt(3)+ct0-ct1 ; ett(3)=ett(3)+et0-et1

          do n=1,nn
             m=n+ns-1
             if ( rb(n)/res(m)>1.d8 ) then
                E(n)=E1(n)
                cycle
             end if
!$OMP parallel do
             do i=n1,n2
                unk(i,m) = utmp2(1,1)*unk(i,m) + utmp2(2,1)*pk(i,n)
             end do
!$OMP end parallel do
          end do

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

       end do ! icg

       esp(ns:ne)=E(1:nn)

    end do  ! band-loop

    deallocate( unk_tmp )
    deallocate( uSu )
    deallocate( Sf )
    deallocate( btmp2,utmp2  )
    deallocate( wtmp2,vtmp2  )
    deallocate( bk,gkgk,E1,E )
    deallocate( rb,sb )
    deallocate( pko,pk  )
    deallocate( Pgk,gk  )
    deallocate( hpk,hxk )

    if ( disp_switch ) then
       write(*,*) "time(hamil_kin)",ctt_hamil(1),ett_hamil(1)
       write(*,*) "time(hamil_loc)",ctt_hamil(2),ett_hamil(2)
       write(*,*) "time(hamil_nlc)",ctt_hamil(3),ett_hamil(3)
       write(*,*) "time(hamil_exx)",ctt_hamil(4),ett_hamil(4)
       write(*,*) "time(hamil_cg)",ctt(1),ett(1)
       write(*,*) "time(op_cg_u )",ctt(2),ett(2)
       write(*,*) "time(com_cg_u)",ctt(3),ett(3)
       write(*,*) "time(pc_cg_u )",ctt(4),ett(4)
    end if

  END SUBROUTINE cg_u

#endif

END MODULE cg_u_module
