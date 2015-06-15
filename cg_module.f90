MODULE cg_module

  use rgrid_module, only: zdV,dV
  use hamiltonian_module
  use cgpc_module
  use parallel_module
  use array_bound_module, only: ML_0,ML_1,MB_0,MB_1
  use cg_lobpcg_module, only: init_lobpcg, lobpcg
  use cg_u_module, only: init_cg_u, cg_u
  use wf_module, only: hunk, iflag_hunk
  use kinetic_module, only: SYStype
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: conjugate_gradient,read_cg,Ncg,iswitch_gs,read_oldformat_cg

  integer :: Ncg,iswitch_gs
  integer :: iswitch_cg

CONTAINS


  SUBROUTINE read_cg(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(4) :: cbuf,ckey
    Ncg = 2
    iswitch_gs = 0
    iswitch_cg = 1
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:3) == "NCG" ) then
             backspace(unit)
             read(unit,*) cbuf,Ncg
          else if ( ckey(1:4) == "SWGS" ) then
             backspace(unit)
             read(unit,*) cbuf,iswitch_gs
          else if ( ckey(1:3) == "ICG" ) then
             backspace(unit)
             read(unit,*) cbuf,iswitch_cg
          end if
       end do
999    continue
       write(*,*) "Ncg=",Ncg
       write(*,*) "iswitch_gs=",iswitch_gs
       write(*,*) "iswitch_cg=",iswitch_cg
    end if
    call send_cg(0)
  END SUBROUTINE read_cg


  SUBROUTINE read_oldformat_cg(rank,unit)
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) Ncg,iswitch_gs
       write(*,*) "Ncg=",Ncg
       write(*,*) "iswitch_gs=",iswitch_gs
    end if
    call send_cg(0)
  END SUBROUTINE read_oldformat_cg


  SUBROUTINE send_cg(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    call mpi_bcast(Ncg,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(iswitch_gs,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(iswitch_cg,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_cg


  SUBROUTINE conjugate_gradient(n1,n2,MB,k,s,Mcg,igs,unk,esp,res)
    implicit none
    integer,intent(IN) :: n1,n2,MB,k,s,Mcg,igs
    real(8),intent(INOUT) :: esp(MB),res(MB)
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: unk(n1:n2,MB)
#else
    complex(8),intent(INOUT) :: unk(n1:n2,MB)
#endif
    integer :: ipc

    call init_cgpc( n1, n2, k, s, dV, SYStype, ipc )

    select case( iswitch_cg )
    case default
    case( 1 )
       if ( disp_switch_parallel ) &
            write(*,'("--- CG ( with IPC=",i1," ) ---")') ipc
       call conjugate_gradient_1(n1,n2,MB,k,s,Mcg,igs,unk,esp,res)
    case( 2 )
       if ( disp_switch_parallel ) &
            write(*,'("--- LOBPCG ( with IPC=",i1," ) ---")') ipc
       call init_lobpcg( n1,n2,MB_0,MB_1,dV,MB_d,comm_grid )
       call lobpcg( k,s,Mcg,igs,unk,esp,res )
    case( 3 )
       if ( disp_switch_parallel ) &
            write(*,'("--- CG_U ( with IPC=",i1," ) ---")') ipc
       call init_cg_u( n1,n2,MB_0,MB_1,dV,MB_d,comm_grid )
       call cg_u( k,s,Mcg,igs,unk,esp,res,disp_switch_parallel )
    end select

  END SUBROUTINE conjugate_gradient


#ifdef _DRSDFT_
  SUBROUTINE conjugate_gradient_1(n1,n2,MB,k,s,Mcg,igs,unk,esp,res)
    implicit none
    integer,intent(IN) :: n1,n2,MB,k,s,Mcg,igs
    real(8),intent(INOUT) :: unk(n1:n2,MB)
    real(8),intent(INOUT) :: esp(MB),res(MB)
    integer :: ns,ne,nn,n,m,icg,ML0,Nhpsi,Npc,Ncgtot,ierr
    integer :: mm,icmp,i,TYPE_MAIN,timer_counter
    real(8),parameter :: ep0=0.d0
    real(8),parameter :: ep1=1.d-15
    real(8) :: rwork(9),W(2),c,d,r,c1,ct0,ct1,et0,et1,ctt(4),ett(4)
    real(8) :: timecg(2,16), ttmp(2), ttmp_cg(2)
    real(8),allocatable :: sb(:),rb(:),E(:),E1(:),gkgk(:),bk(:)
    complex(8) :: work(9),zphase,ztmp
    real(8),parameter :: zero=0.d0
    real(8),allocatable :: hxk(:,:),hpk(:,:),gk(:,:),Pgk(:,:)
    real(8),allocatable :: pk(:,:),pko(:,:)
    real(8),allocatable :: vtmp2(:,:),wtmp2(:,:)
    real(8),allocatable :: utmp2(:,:),btmp2(:,:)

    call watchb( ttmp_cg ) ; ttmp(:)=ttmp_cg(:)

    TYPE_MAIN = MPI_REAL8

    timecg(:,:)=0.0d0
    time_hmlt(:,:)=0.0d0
    time_kine(:,:)=0.0d0
    time_nlpp(:,:)=0.0d0
    time_cgpc(:,:)=0.0d0

    ML0 = ML_1-ML_0+1

    mm  = ML0  ; if (TYPE_MAIN==mpi_complex16) mm=2*ML0
    c1  = 2.d0 ; if (TYPE_MAIN==mpi_complex16) c1=1.d0
    icmp= 1    ; if (TYPE_MAIN==mpi_complex16) icmp=2

    Ncgtot = 0
    Nhpsi  = 0
    Npc    = 0

    allocate( hxk(n1:n2,MB_d), hpk(n1:n2,MB_d) )
    allocate( gk(n1:n2,MB_d) , Pgk(n1:n2,MB_d) )
    allocate( pk(n1:n2,MB_d) , pko(n1:n2,MB_d) )
    allocate( sb(MB_d),rb(MB_d) )
    allocate( E(MB_d),E1(MB_d),gkgk(MB_d),bk(MB_d) )
    allocate( vtmp2(6,MB_d),wtmp2(6,MB_d) )
    allocate( utmp2(2,2),btmp2(2,2) )

!$OMP parallel workshare
    res(:)  = 0.d0
    esp(:)  = 0.d0
!$OMP end parallel workshare

    call watchb( ttmp, timecg(:,5) )

    do ns=MB_0,MB_1
       ne=ns
       nn=ne-ns+1

       E1(:)=1.d10

       call watchb( ttmp )

       if ( iflag_hunk >= 1 ) then
!$OMP parallel workshare
          hxk(:,1:nn)=hunk(:,ns:ne,k,s)
!$OMP end parallel workshare
       else
          call hamiltonian(k,s,unk(n1,ns),hxk,n1,n2,ns,ne) ; Nhpsi=Nhpsi+1
       end if

       call watchb( ttmp, timecg(:,1) )

       do n=1,nn
          call dot_product(unk(n1,n+ns-1),hxk(n1,n),sb(n),dV,mm,1)
       end do

       call watchb( ttmp, timecg(:,2) )

       call mpi_allreduce(sb,E,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       call watchb( ttmp, timecg(:,3) )

       do n=1,nn
!$OMP parallel do
          do i=n1,n2
             gk(i,n)=-c1*(hxk(i,n)-E(n)*unk(i,n+ns-1))
          end do
!$OMP end parallel do
          call dot_product(gk(n1,n),gk(n1,n),sb(n),dV,mm,1)
       end do

       call watchb( ttmp, timecg(:,2) )

       call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       call watchb( ttmp, timecg(:,3) )

       do icg=1,Mcg+1

          call watchb( ttmp )

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

          call watchb( ttmp, timecg(:,2) )

! --- Preconditioning ---

          call preconditioning(E,k,s,nn,ML0,unk(n1,ns),gk,Pgk)

          call watchb( ttmp, timecg(:,4) )

! ---

          do n=1,nn
             call dot_product(Pgk(n1,n),gk(n1,n),sb(n),dV,mm,1)
          end do

          call watchb( ttmp, timecg(:,2) )

          call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

          call watchb( ttmp, timecg(:,3) )

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

          call watchb( ttmp, timecg(:,2) )

          call hamiltonian(k,s,pk,hpk,n1,n2,ns,ne) ; Nhpsi=Nhpsi+1

          call watchb( ttmp, timecg(:,1) )

          do n=1,nn
             vtmp2(1:6,n)=zero
             m=n+ns-1
             call dot_product(unk(n1,m),unk(n1,m),vtmp2(1,n),dV,mm,1)
             call dot_product(pk(n1,n),unk(n1,m),vtmp2(2,n),dV,mm,icmp)
             call dot_product(pk(n1,n),pk(n1,n),vtmp2(3,n),dV,mm,1)
             call dot_product(unk(n1,m),hxk(n1,n),vtmp2(4,n),dV,mm,1)
             call dot_product(pk(n1,n),hxk(n1,n),vtmp2(5,n),dV,mm,icmp)
             call dot_product(pk(n1,n),hpk(n1,n),vtmp2(6,n),dV,mm,1)
          end do

          call watchb( ttmp, timecg(:,2) )

          call mpi_allreduce(vtmp2,wtmp2,6*nn,TYPE_MAIN,mpi_sum,comm_grid,ierr)

          call watchb( ttmp, timecg(:,3) )

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

!$OMP parallel
!$OMP do
             do i=n1,n2
                hxk(i,n)=utmp2(1,1)*hxk(i,n)+utmp2(2,1)*hpk(i,n)
             end do
!$OMP end do
!$OMP do
             do i=n1,n2
                gk(i,n) = -c1*( hxk(i,n) &
                     -W(1)*(utmp2(1,1)*unk(i,m)+utmp2(2,1)*pk(i,n)) )
             end do
!$OMP end do
!$OMP end parallel

             call dot_product(gk(n1,n),gk(n1,n),sb(n),dV,mm,1)

          end do ! n

          call watchb( ttmp, timecg(:,2) )

          call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

          call watchb( ttmp, timecg(:,3) )

          do n=1,nn
             m=n+ns-1
             if ( rb(n)/res(m)>1.d8 ) then
                E(n)=E1(n)
                cycle
             end if
!$OMP parallel do
             do i=n1,n2
                unk(i,m)=utmp2(1,1)*unk(i,m)+utmp2(2,1)*pk(i,n)
             end do
!$OMP end parallel do
             if ( iflag_hunk >= 1 ) then
!$OMP parallel do
                do i=n1,n2
                   hunk(i,m,k,s)=hxk(i,n)
                end do
!$OMP end parallel do
             end if
          end do

          call watchb( ttmp, timecg(:,2) )

       end do ! icg

       esp(ns:ne)=E(1:nn)

    end do  ! band-loop

    call watchb( ttmp )

    deallocate( btmp2,utmp2  )
    deallocate( wtmp2,vtmp2  )
    deallocate( bk,gkgk,E1,E )
    deallocate( rb,sb )
    deallocate( pko,pk  )
    deallocate( Pgk,gk  )
    deallocate( hpk,hxk )

    call watchb( ttmp, timecg(:,6) )
    call watchb( ttmp_cg, timecg(:,7) )

    if ( disp_switch_parallel ) then
       write(*,*) "time(hmlt_kin)",( time_hmlt(i,1), i=1,2 )
       write(*,*) "time(hmlt_loc)",( time_hmlt(i,2), i=1,2 )
       write(*,*) "time(hmlt_nlc)",( time_hmlt(i,3), i=1,2 )
       write(*,*) "time(hmlt_exx)",( time_hmlt(i,4), i=1,2 )
       call write_watchb( time_kine(1,1),11, "kine" ) 
       call write_watchb( time_nlpp(1,1), 4, "nlpp" ) 
       call write_watchb( time_cgpc(1,1),13, "cgpc" ) 
       call write_watchb( timecg(1,1), 7, "dcg" ) 
    end if

  END SUBROUTINE conjugate_gradient_1

#else

  SUBROUTINE conjugate_gradient_1(n1,n2,MB,k,s,Mcg,igs,unk,esp,res)
    implicit none
    integer,intent(IN) :: n1,n2,MB,k,s,Mcg,igs
    complex(8),intent(INOUT) :: unk(n1:n2,MB)
    real(8),intent(INOUT) :: esp(MB),res(MB)
    integer :: ns,ne,nn,n,m,icg,ML0,Nhpsi,Npc,Ncgtot,ierr
    integer :: mm,icmp,i,TYPE_MAIN
    real(8),parameter :: ep0=0.d0
    real(8),parameter :: ep1=1.d-15
    real(8) :: rwork(9),W(2),c,d,r,c1,ct0,ct1,et0,et1,ctt(4),ett(4)
    real(8),allocatable :: sb(:),rb(:),E(:),E1(:),gkgk(:),bk(:)
    complex(8) :: work(9),zphase,ztmp

    complex(8),parameter :: zero=(0.d0,0.d0)
    complex(8),allocatable :: hxk(:,:),hpk(:,:),gk(:,:),Pgk(:,:)
    complex(8),allocatable :: pk(:,:),pko(:,:)
    complex(8),allocatable :: vtmp2(:,:),wtmp2(:,:)
    complex(8),allocatable :: utmp2(:,:),btmp2(:,:)

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

    allocate( hxk(n1:n2,MB_d), hpk(n1:n2,MB_d) )
    allocate( gk(n1:n2,MB_d) , Pgk(n1:n2,MB_d) )
    allocate( pk(n1:n2,MB_d) , pko(n1:n2,MB_d) )
    allocate( sb(MB_d),rb(MB_d) )
    allocate( E(MB_d),E1(MB_d),gkgk(MB_d),bk(MB_d) )
    allocate( vtmp2(6,MB_d),wtmp2(6,MB_d) )
    allocate( utmp2(2,2),btmp2(2,2) )

!$OMP parallel workshare
    res(:)  = 0.d0
    esp(:)  = 0.d0
!$OMP end parallel workshare

    do ns=MB_0,MB_1
       ne=ns
       nn=ne-ns+1

       E1(:)=1.d10

       call watch(ct0,et0)

       if ( iflag_hunk >= 1 ) then
!$OMP parallel workshare
          hxk(:,1:nn)=hunk(:,ns:ne,k,s)
!$OMP end parallel workshare
       else
          call hamiltonian(k,s,unk(n1,ns),hxk,n1,n2,ns,ne) ; Nhpsi=Nhpsi+1
       end if

       call watch(ct1,et1) ; ctt(1)=ctt(1)+ct1-ct0 ; ett(1)=ett(1)+et1-et0

       do n=1,nn
          call dot_product(unk(n1,n+ns-1),hxk(n1,n),sb(n),dV,mm,1)
       end do

       call watch(ct0,et0) ; ctt(2)=ctt(2)+ct0-ct1 ; ett(2)=ett(2)+et0-et1

       call mpi_allreduce(sb,E,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       call watch(ct1,et1) ; ctt(3)=ctt(3)+ct1-ct0 ; ett(3)=ett(3)+et1-et0

       do n=1,nn
!$OMP parallel do
          do i=n1,n2
             gk(i,n)=-c1*(hxk(i,n)-E(n)*unk(i,n+ns-1))
          end do
!$OMP end parallel do
          call dot_product(gk(n1,n),gk(n1,n),sb(n),dV,mm,1)
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
             call dot_product(Pgk(n1,n),gk(n1,n),sb(n),dV,mm,1)
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
             vtmp2(1:6,n)=zero
             m=n+ns-1
             call dot_product(unk(n1,m),unk(n1,m),vtmp2(1,n),dV,mm,1)
             call dot_product(pk(n1,n),unk(n1,m),vtmp2(2,n),dV,mm,icmp)
             call dot_product(pk(n1,n),pk(n1,n),vtmp2(3,n),dV,mm,1)
             call dot_product(unk(n1,m),hxk(n1,n),vtmp2(4,n),dV,mm,1)
             call dot_product(pk(n1,n),hxk(n1,n),vtmp2(5,n),dV,mm,icmp)
             call dot_product(pk(n1,n),hpk(n1,n),vtmp2(6,n),dV,mm,1)
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

!$OMP parallel
!$OMP do
             do i=n1,n2
                hxk(i,n)=utmp2(1,1)*hxk(i,n)+utmp2(2,1)*hpk(i,n)
             end do
!$OMP end do
!$OMP do
             do i=n1,n2
                gk(i,n) = -c1*( hxk(i,n) &
                     -W(1)*(utmp2(1,1)*unk(i,m)+utmp2(2,1)*pk(i,n)) )
             end do
!$OMP end do
!$OMP end parallel

             call dot_product(gk(n1,n),gk(n1,n),sb(n),dV,mm,1)

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
                unk(i,m)=utmp2(1,1)*unk(i,m)+utmp2(2,1)*pk(i,n)
             end do
!$OMP end parallel do
             if ( iflag_hunk >= 1 ) then
!$OMP parallel do
                do i=n1,n2
                   hunk(i,m,k,s)=hxk(i,n)
                end do
!$OMP end parallel do
             end if
          end do

          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

       end do ! icg

       esp(ns:ne)=E(1:nn)

    end do  ! band-loop

    deallocate( btmp2,utmp2  )
    deallocate( wtmp2,vtmp2  )
    deallocate( bk,gkgk,E1,E )
    deallocate( rb,sb )
    deallocate( pko,pk  )
    deallocate( Pgk,gk  )
    deallocate( hpk,hxk )

    if ( disp_switch_parallel ) then
!       write(*,*) "time(hamil_kin)",ctt_hamil(1),ett_hamil(1)
!       write(*,*) "time(hamil_loc)",ctt_hamil(2),ett_hamil(2)
!       write(*,*) "time(hamil_nlc)",ctt_hamil(3),ett_hamil(3)
!       write(*,*) "time(hamil_exx)",ctt_hamil(4),ett_hamil(4)
!       write(*,*) "time(hamil_cg)",ctt(1),ett(1)
!       write(*,*) "time(op_cg   )",ctt(2),ett(2)
!       write(*,*) "time(com_cg  )",ctt(3),ett(3)
!       write(*,*) "time(pc_cg   )",ctt(4),ett(4)
    end if

  END SUBROUTINE conjugate_gradient_1

#endif


#ifdef TEST
  SUBROUTINE dot_product(a,b,c,alpha,n,m)
    implicit none
    integer,intent(IN) :: n,m
    integer :: i
    real(8) :: a(*),b(*),c(*),alpha,tmp

    c(1:m)=0.d0

    tmp=0.d0
!$OMP parallel do reduction(+:tmp)
    do i=1,n
       tmp=tmp+a(i)*b(i)
    end do
!$OMP end parallel do
    c(1)=tmp*alpha

    if ( m==2 ) then
       tmp=0.d0
!$OMP parallel do reduction(+:tmp)
       do i=1,n-1,2
          tmp=tmp+a(i)*b(i+1)-a(i+1)*b(i)
       end do
!$OMP end parallel do
       c(m)=tmp*alpha
    end if

    return
  END SUBROUTINE dot_product
#endif

END MODULE cg_module
