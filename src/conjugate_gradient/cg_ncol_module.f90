MODULE cg_ncol_module

  use rgrid_module, only: zdV,dV
  use hamiltonian_module
  use hamiltonian_ncol_module
  use cgpc_module
  use parallel_module, only: RSDFT_MPI_COMPLEX16,comm_grid,id_bzsm,myrank_k
  use watch_module
  use io_tools_module
  use noncollinear_module, only: flag_noncollinear

  implicit none

  PRIVATE
  PUBLIC :: conjugate_gradient_ncol
  PUBLIC :: flag_noncollinear

  integer :: Ncg = 3
  integer :: iswitch_gs = 0
  integer :: iswitch_cg = 1
  logical :: flag_init_read = .true.
  logical :: flag_init_cg = .true.

CONTAINS


  SUBROUTINE init_cg
    implicit none
    logical :: disp
    if ( .not.flag_init_cg ) return
    call check_disp_switch( disp, 0 )
    if ( disp ) then
       write(*,*) "NCG =",Ncg
       write(*,*) "ICG =",iswitch_cg
       write(*,*) "SWGS=",iswitch_gs
    end if
    flag_init_cg=.false.
  END SUBROUTINE init_cg


  SUBROUTINE conjugate_gradient_ncol( n1,n2, MB, k, unk, esp, res )
    implicit none
    integer,intent(IN) :: n1,n2,MB,k
    real(8),intent(INOUT) :: esp(MB),res(MB)
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: unk(:,:,:,:)
#else
    complex(8),intent(INOUT) :: unk(:,:,:,:)
#endif
    type(time) :: tt

    call write_border( 1, " conjugate_gradient_ncol(start)" )

    call start_timer( t_out=tt )

    call init_cg

!    if ( pp_kind == "USPP" ) then
!       call conjugate_gradient_g( n1,n2,MB,k,s,Ncg,unk,esp,res,iswitch_gs )
!    else

       select case( iswitch_cg )
       case default
          call stop_program( "iswitch_cg /= 1 is not available in noncol " )
       case( 1 )
          call conjugate_gradient_ncol_1(n1,n2,MB,k,Ncg,unk,esp,res)
!       case( 2 )
!          call init_lobpcg( n1,n2,MB_0,MB_1,dV,MB_d,comm_grid )
!          call lobpcg( k,s,Ncg,iswitch_gs,unk,esp,res )
       end select

!    end if

    call result_timer( "cg", tt )

    call write_border( 1, " conjugate_gradient_ncol(end)" )

  END SUBROUTINE conjugate_gradient_ncol


  SUBROUTINE conjugate_gradient_ncol_1(n1,n2,MB,k,Mcg,unk,esp,res)

    implicit none
    integer,intent(IN) :: n1,n2,MB,k,Mcg
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: unk(:,:,:,:)
#else
    complex(8),intent(INOUT) :: unk(:,:,:,:)
#endif
    real(8),intent(INOUT) :: esp(MB),res(MB)
    integer :: ns,ne,nn,n,m,icg,ML0,Nhpsi,Npc,Ncgtot,ierr
    integer :: mm,i,s,MS,k0
    real(8),parameter :: ep0=0.d0
    real(8),parameter :: ep1=1.d-15
    real(8) :: rwork(9),W(2),c,d,r,c1
    complex(8) :: sb,rb,E,E1,gkgk,bk
    complex(8) :: work(9),zphase,ztmp
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: xk(:,:,:),hxk(:,:,:),hpk(:,:,:),gk(:,:,:),Pgk(:,:,:)
    complex(8),allocatable :: pk(:,:,:),pko(:,:)
    complex(8) :: vtmp2(2,2),wtmp2(2,2)
    complex(8) :: utmp2(2,2),btmp2(2,2)
    include 'mpif.h'

    ML0 = size(unk,1)
    MS  = size(unk,4)
    k0  = k-id_bzsm(myrank_k)

    mm = 2*ML0
    c1 = 1.0d0

    Ncgtot = 0
    Nhpsi  = 0
    Npc    = 0

    allocate(  xk(ML0,1,MS) ) ;  xk=zero
    allocate( hxk(ML0,1,MS) ) ; hxk=zero
    allocate( hpk(ML0,1,MS) ) ; hpk=zero
    allocate(  gk(ML0,1,MS) ) ;  gk=zero
    allocate( Pgk(ML0,1,MS) ) ; Pgk=zero
    allocate(  pk(ML0,1,MS) ) ;  pk=zero

!$OMP parallel workshare
    res(:)  = 0.d0
    esp(:)  = 0.d0
!$OMP end parallel workshare

    do n=1,MB

       E1=1.d10

       do s=1,MS
          xk(:,1,s) = unk(:,n,k0,s)
       end do

       do s=1,MS
#ifndef _DRSDFT_
          call hamiltonian( xk(:,:,s), hxk(:,:,s), n,k,s ) ; Nhpsi=Nhpsi+1
#endif
       end do
       call hamiltonian_ncol( k, n1,n2, xk, hxk )

       sb = sum( conjg(xk)*hxk )*dV

       call mpi_allreduce(sb,rb,1,MPI_COMPLEX16,mpi_sum,comm_grid,ierr)

       E = sb

       do s=1,MS
!$OMP parallel do
          do i=1,ML0
             gk(i,1,s) = -c1*( hxk(i,1,s) - E*xk(i,1,s) )
          end do
!$OMP end parallel do
       end do

       sb = sum( conjg(gk)*gk )*dV

       call mpi_allreduce(sb,rb,1,MPI_COMPLEX16,mpi_sum,comm_grid,ierr)

       do icg=1,Mcg+1

          Ncgtot=Ncgtot+1

          do s=1,MS
!$OMP parallel workshare
             Pgk(:,1,s)=gk(:,1,s)
!$OMP end parallel workshare
          end do

          res(n)=rb/c1**2

! --- Convergence check ---

          if ( res(n) < ep0 ) exit
          if ( abs(E-E1) < ep1 ) exit
          if ( icg == Mcg+1 ) exit

! --- Preconditioning ---
          W(1)=E
          do s=1,MS
#ifndef _DRSDFT_
             call preconditioning( W(1),k,s,1,ML0,xk(:,:,s),gk(:,:,s),Pgk(:,:,s) )
#endif
          end do
! --- orthogonalization
!          do n=ns,ne
!             call cggs( iswitch_gs, ML0, MB, n, dV, unk(n1,1), Pgk(n1,n-ns+1) )
!          end do
! ---

          sb = sum( conjg(gk)*Pgk )*dV

          call mpi_allreduce(sb,rb,1,mpi_complex16,mpi_sum,comm_grid,ierr)

          if ( icg == 1 ) then

!$OMP parallel workshare
             pk(1:ML0,:,1:MS) = Pgk(1:ML0,:,1:MS)
!$OMP end parallel workshare

          else

             bk = rb/gkgk

             do s=1,MS
!$OMP parallel do
                do i=1,ML0
                   pk(i,1,s) = Pgk(i,1,s) + bk*pk(i,1,s)
                end do
!$OMP end parallel do
             end do

          end if
          gkgk = rb

          do s=1,MS
#ifndef _DRSDFT_
             call hamiltonian( pk(:,:,s), hpk(:,:,s), n,k,s ) ; Nhpsi=Nhpsi+1
#endif
          end do
          call hamiltonian_ncol( k, n1,n2, pk, hpk )

          vtmp2(1,1) = sum( conjg(xk)*xk  )*dV
          vtmp2(1,2) = sum( conjg(xk)*pk  )*dV
          vtmp2(2,1) = sum( conjg(pk)*xk  )*dV
          vtmp2(2,2) = sum( conjg(pk)*pk  )*dV

          wtmp2(1,1) = sum( conjg(xk)*hxk )*dV
          wtmp2(1,2) = sum( conjg(xk)*hpk )*dV
          wtmp2(2,1) = sum( conjg(pk)*hxk )*dV
          wtmp2(2,2) = sum( conjg(pk)*hpk )*dV

          call mpi_allreduce(vtmp2,btmp2,4,MPI_COMPLEX16,mpi_sum,comm_grid,ierr)
          call mpi_allreduce(wtmp2,utmp2,4,MPI_COMPLEX16,mpi_sum,comm_grid,ierr)

          call zhegv(1,'V','U',2,utmp2,2,btmp2,2,W,work,9,rwork,ierr)

          if ( abs(W(1)-E) > 1.d-1 .and. abs(W(2)-E) <= 1.d-1 ) then
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

          E1 = E
          E  = W(1)

          do s=1,MS
!$OMP parallel
!$OMP do
             do i=1,ML0
                hxk(i,1,s)=utmp2(1,1)*hxk(i,1,s)+utmp2(2,1)*hpk(i,1,s)
             end do
!$OMP end do
!$OMP do
             do i=1,ML0
                gk(i,1,s) = -c1*( hxk(i,1,s) &
                               -W(1)*(utmp2(1,1)*xk(i,1,s)+utmp2(2,1)*pk(i,1,s)) )
             end do
!$OMP end do
!$OMP end parallel
          end do ! s

          sb = sum( conjg(gk)*gk )*dV

          call mpi_allreduce(sb,rb,1,MPI_COMPLEX16,mpi_sum,comm_grid,ierr)

          if ( abs(rb)/res(n) > 1.d8 ) then
             E = E1
             cycle
          end if

          do s=1,MS
!$OMP parallel do
             do i=1,ML0
                xk(i,1,s)=utmp2(1,1)*xk(i,1,s)+utmp2(2,1)*pk(i,1,s)
             end do
!$OMP end parallel do
          end do

       end do ! icg

       esp(n)=E
       do s=1,MS
          unk(:,n,k0,s) = xk(:,1,s)
       end do

    end do ! n (band)

    deallocate( pk,Pgk,gk,hpk,hxk )

  END SUBROUTINE conjugate_gradient_ncol_1


END MODULE cg_ncol_module
