MODULE cg_module

  use rgrid_module, only: zdV,dV
  use hamiltonian_module
  use cgpc_module
  use parallel_module
  use array_bound_module, only: ML_0,ML_1,MB_0,MB_1

  implicit none

  PRIVATE
  PUBLIC :: conjugate_gradient,read_cg,send_cg,Ncg,iswitch_gs

  integer :: Ncg,iswitch_gs

CONTAINS


  SUBROUTINE read_cg(unit)
    integer,intent(IN) :: unit
    read(unit,*) Ncg,iswitch_gs
    write(*,*) "Ncg=",Ncg
    write(*,*) "iswitch_gs=",iswitch_gs
  END SUBROUTINE read_cg


  SUBROUTINE send_cg(rank)
    integer,intent(IN) :: rank
    integer :: ierr
    call mpi_bcast(Ncg,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(iswitch_gs,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_cg

#ifdef _DRSDFT_
  SUBROUTINE conjugate_gradient(n1,n2,MB,k,s,Mcg,igs,unk,esp,res)
    integer,intent(IN) :: n1,n2,MB,k,s,Mcg,igs
    real(8),intent(INOUT) :: unk(n1:n2,MB)
    real(8),intent(INOUT) :: esp(MB),res(MB)
    integer :: ns,ne,nn,n,m,icg,ierr
    integer :: mm,i,ld,j,ML0
    real(8),parameter :: ep0=0.d0
    real(8),parameter :: ep1=1.d-15
    real(8) :: rwork(999),W(999),c,d,r
    real(8),allocatable :: sb(:),rb(:),E(:),E1(:)
    real(8) :: work(999),zphase,ztmp,a1,a2,a3

    real(8),allocatable :: v(:,:,:),hv(:,:,:)
    real(8),allocatable :: vt(:,:),ut(:,:),wt(:,:)
    real(8),allocatable :: vv0(:,:,:),vv1(:,:,:)
    real(8),parameter :: zero=0.d0,one=1.d0

!    n1  = ML_0
!    n2  = ML_1
    ML0 = ML_1-ML_0+1
!    MB  = Nband

!    c1  = 1.d0

    ld  = 3*MB_d

    allocate( sb(MB_d),rb(MB_d) )
    allocate( E(MB_d),E1(MB_d) )
    allocate( v(n1:n2,MB_d,4), hv(n1:n2,MB_d,3) )
    allocate( vt(n1:n2,MB_d*3), ut(n1:n2,MB_d*3) )
    allocate( wt(n1:n2,MB_d) )
    allocate( vv0(ld,ld,2),vv1(ld,ld,2) )
    v=zero ; hv=zero

    res(:) = 0.d0
    esp(:) = 0.d0

    do ns=MB_0,MB_1,MB_d

       ne=min(ns+MB_d-1,MB_1)
       nn=ne-ns+1

       E1(:)=1.d10

!$OMP parallel

       do n=1,nn
!$OMP do
          do i=n1,n2
             v(i,n,1) = unk(i,ns+n-1)
          end do
!$OMP end do
       end do

!$OMP end parallel

       call hamiltonian(k,s,v,hv,n1,n2,ns,ne)

!$OMP parallel

       do n=1,nn
!$OMP single
          c=0.d0
!$OMP end single
!$OMP do reduction(+:c)
          do i=n1,n2
             c=c+v(i,n,1)*hv(i,n,1)
          end do
!$OMP end do
!$OMP single
          sb(n)=c*dV
!$OMP end single
       end do

!$OMP single
       call mpi_allreduce(sb,E,nn,mpi_real8,mpi_sum,comm_grid,ierr)
!$OMP end single

       do n=1,nn
!$OMP do
          do i=n1,n2
             v(i,n,4) = -( hv(i,n,1) - E(n)*v(i,n,1) )
          end do
!$OMP end do
!$OMP single
          c=0.d0
!$OMP end single
!$OMP do reduction(+:c)
          do i=n1,n2
             c=c+abs(v(i,n,4))**2
          end do
!$OMP end do
!$OMP single
          sb(n)=c*dV
!$OMP end single
       end do

!$OMP single
       call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)
!$OMP end single
!$OMP end parallel

       do icg=1,Ncg+1

!$OMP parallel

          do n=1,nn
!$OMP do
             do i=n1,n2
                v(i,n,3)=v(i,n,4)
             end do
!$OMP end do
          end do

!$OMP single
          res(ns:ne)=rb(1:nn)
!$OMP end single

!$OMP end parallel

! --- Convergence check ---

          if ( all(rb(1:nn) < ep0) ) exit
          if ( all(abs(E(1:nn)-E1(1:nn)) < ep1) ) exit
          if ( icg == Ncg+1 ) exit

! --- Preconditioning ---

          call preconditioning(E,k,s,nn,ML0,v(n1,1,4),v(n1,1,3))

! ---

!$OMP parallel

          if ( icg == 1 ) then
             do n=1,nn
!$OMP do
                do i=n1,n2
                   v(i,n,2)=v(i,n,3)
                end do
!$OMP end do
             end do
          end if

!$OMP single
          vv0=zero
!$OMP end single

!$OMP end parallel

          if ( icg==1 ) then
             mm=2*nn
             call hamiltonian(k,s,v(n1,1,2),hv(n1,1,2),n1,n2,ns,ne)
          else
             mm=3*nn
             call hamiltonian(k,s,v(n1,1,3),hv(n1,1,3),n1,n2,ns,ne)
          end if

!$OMP parallel

!$OMP single
          m=0
!$OMP end single
          do j=1,3
             do n=1,nn
!$OMP single
                m=m+1
!$OMP end single
!$OMP do
                do i=n1,n2
                   vt(i,m)= v(i,n,j)
                   ut(i,m)=hv(i,n,j)
                end do
!$OMP end do
             end do
          end do

!$OMP end parallel


!$OMP single
          ztmp=0.5d0*dV
!$OMP end single

!          call zherk( 'U','C',mm,ML0,zdV,vt,ML0,zero,vv0(1,1,1),ld)
!          call zher2k('U','C',mm,ML0,ztmp,vt,ML0,ut,ML0,zero,vv0(1,1,2),ld)
          call dsyrk( 'U','T',mm,ML0,zdV,vt,ML0,zero,vv0(1,1,1),ld)
          call dsyr2k('U','T',mm,ML0,ztmp,vt,ML0,ut,ML0,zero,vv0(1,1,2),ld)

!$OMP single
          call mpi_allreduce &
               (vv0,vv1,ld*ld*2,MPI_REAL8,mpi_sum,comm_grid,ierr)
!$OMP end single

!          call zhegv &
!          (1,'V','U',mm,vv1(1,1,2),ld,vv1(1,1,1),ld,W,work,999,rwork,ierr)
          call dsygv &
          (1,'V','U',mm,vv1(1,1,2),ld,vv1(1,1,1),ld,W,work,999,ierr)

!$OMP single
          E1(1:nn)=E(1:nn)
          E(1:nn) =W(1:nn)
!$OMP end single

          call dgemm('N','N',ML0,nn,mm,one,vt,ML0,vv1(1,1,2),ld,zero,wt,ML0)

!$OMP parallel
          do n=1,nn
!$OMP do
             do i=n1,n2
                v(i,n,1) = wt(i,n)
             end do
!$OMP end do
          end do
!$OMP end parallel

          call dgemm('N','N',ML0,nn,mm,one,ut,ML0,vv1(1,1,2),ld,zero,wt,ML0)

!$OMP parallel
          do n=1,nn
!$OMP do
             do i=n1,n2
                hv(i,n,1)=wt(i,n)
             end do
!$OMP end do
          end do
!$OMP end parallel

          call dgemm('N','N',ML0,nn,mm-nn,one,vt(n1,nn+1),ML0 &
                    ,vv1(nn+1,nn+1,2),ld,zero,wt,ML0)

!$OMP parallel
          do n=1,nn
!$OMP do
             do i=n1,n2
                v(i,n,2)=wt(i,n)
             end do
!$OMP end do
          end do
!$OMP end parallel

          call dgemm('N','N',ML0,nn,mm-nn,one,ut(n1,nn+1),ML0 &
                    ,vv1(nn+1,nn+1,2),ld,zero,wt,ML0)
!$OMP parallel
          do n=1,nn
!$OMP do
             do i=n1,n2
                hv(i,n,2)=wt(i,n)
             end do
!$OMP end do
          end do

          do n=1,nn
!$OMP do
             do i=n1,n2
                v(i,n,4) = -( hv(i,n,1) - W(n)*v(i,n,1) )
             end do
!$OMP end do
!$OMP single
             c=0.d0
!$OMP end single
!$OMP do reduction(+:c)
             do i=n1,n2
                c=c+abs(v(i,n,4))**2
             end do
!$OMP end do
!$OMP single
             sb(n)=c*dV
!$OMP end single
          end do

!$OMP single
          call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)
!$OMP end single

!$OMP end parallel

       end do ! icg

!$OMP parallel

!$OMP single
       esp(ns:ne)=E(1:nn)
!$OMP end single

       do n=1,nn
!$OMP do
          do i=n1,n2
             unk(i,ns+n-1) = v(i,n,1)
          end do
!$OMP end do
       end do

!$OMP end parallel

    end do  ! band-loop

    deallocate( vv1, vv0 )
    deallocate( wt )
    deallocate( ut, vt )
    deallocate( v, hv )
    deallocate( E,E1 )
    deallocate( sb,rb )

    return

  END SUBROUTINE conjugate_gradient

#else

  SUBROUTINE conjugate_gradient(n1,n2,MB,k,s,Mcg,igs,unk,esp,res)
    integer,intent(IN) :: n1,n2,MB,k,s,Mcg,igs
    complex(8),intent(INOUT) :: unk(n1:n2,MB)
    real(8),intent(INOUT) :: esp(MB),res(MB)
    integer :: ns,ne,nn,n,m,icg,ierr
    integer :: mm,i,ld,j,ML0
    real(8),parameter :: ep0=0.d0
    real(8),parameter :: ep1=1.d-15
    real(8) :: rwork(999),W(999),c,d,r
    real(8),allocatable :: sb(:),rb(:),E(:),E1(:)
    complex(8) :: work(999),zphase,ztmp,a1,a2,a3

    complex(8),allocatable :: v(:,:,:),hv(:,:,:)
    complex(8),allocatable :: vt(:,:),ut(:,:),wt(:,:)
    complex(8),allocatable :: vv0(:,:,:),vv1(:,:,:)
    complex(8),parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)

!    n1  = ML_0
!    n2  = ML_1
    ML0 = ML_1-ML_0+1
!    MB  = Nband

!    c1  = 1.d0

    ld  = 3*MB_d

    allocate( sb(MB_d),rb(MB_d) )
    allocate( E(MB_d),E1(MB_d) )
    allocate( v(n1:n2,MB_d,4), hv(n1:n2,MB_d,3) )
    allocate( vt(n1:n2,MB_d*3), ut(n1:n2,MB_d*3) )
    allocate( wt(n1:n2,MB_d) )
    allocate( vv0(ld,ld,2),vv1(ld,ld,2) )
    v=zero ; hv=zero

    res(:) = 0.d0
    esp(:) = 0.d0

    do ns=MB_0,MB_1,MB_d

       ne=min(ns+MB_d-1,MB_1)
       nn=ne-ns+1

       E1(:)=1.d10

!$OMP parallel

       do n=1,nn
!$OMP do
          do i=n1,n2
             v(i,n,1) = unk(i,ns+n-1)
          end do
!$OMP end do
       end do

!$OMP end parallel

       call hamiltonian(k,s,v,hv,n1,n2,ns,ne)

!$OMP parallel

       do n=1,nn
!$OMP single
          c=0.d0
!$OMP end single
!$OMP do reduction(+:c)
          do i=n1,n2
             c=c+conjg(v(i,n,1))*hv(i,n,1)
          end do
!$OMP end do
!$OMP single
          sb(n)=c*dV
!$OMP end single
       end do

!$OMP single
       call mpi_allreduce(sb,E,nn,mpi_real8,mpi_sum,comm_grid,ierr)
!$OMP end single

       do n=1,nn
!$OMP do
          do i=n1,n2
             v(i,n,4) = -( hv(i,n,1) - E(n)*v(i,n,1) )
          end do
!$OMP end do
!$OMP single
          c=0.d0
!$OMP end single
!$OMP do reduction(+:c)
          do i=n1,n2
             c=c+abs(v(i,n,4))**2
          end do
!$OMP end do
!$OMP single
          sb(n)=c*dV
!$OMP end single
       end do

!$OMP single
       call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)
!$OMP end single
!$OMP end parallel

       do icg=1,Ncg+1

!$OMP parallel

          do n=1,nn
!$OMP do
             do i=n1,n2
                v(i,n,3)=v(i,n,4)
             end do
!$OMP end do
          end do

!$OMP single
          res(ns:ne)=rb(1:nn)
!$OMP end single

!$OMP end parallel

! --- Convergence check ---

          if ( all(rb(1:nn) < ep0) ) exit
          if ( all(abs(E(1:nn)-E1(1:nn)) < ep1) ) exit
          if ( icg == Ncg+1 ) exit

! --- Preconditioning ---

          call preconditioning(E,k,s,nn,ML0,v(n1,1,4),v(n1,1,3))

! ---

!$OMP parallel

          if ( icg == 1 ) then
             do n=1,nn
!$OMP do
                do i=n1,n2
                   v(i,n,2)=v(i,n,3)
                end do
!$OMP end do
             end do
          end if

!$OMP single
          vv0=zero
!$OMP end single

!$OMP end parallel

          if ( icg==1 ) then
             mm=2*nn
             call hamiltonian(k,s,v(n1,1,2),hv(n1,1,2),n1,n2,ns,ne)
          else
             mm=3*nn
             call hamiltonian(k,s,v(n1,1,3),hv(n1,1,3),n1,n2,ns,ne)
          end if

!$OMP parallel

!$OMP single
          m=0
!$OMP end single
          do j=1,3
             do n=1,nn
!$OMP single
                m=m+1
!$OMP end single
!$OMP do
                do i=n1,n2
                   vt(i,m)= v(i,n,j)
                   ut(i,m)=hv(i,n,j)
                end do
!$OMP end do
             end do
          end do

!$OMP end parallel


!$OMP single
          ztmp=0.5d0*dV
!$OMP end single

          call zherk( 'U','C',mm,ML0,zdV,vt,ML0,zero,vv0(1,1,1),ld)
          call zher2k('U','C',mm,ML0,ztmp,vt,ML0,ut,ML0,zero,vv0(1,1,2),ld)

!$OMP single
          call mpi_allreduce &
               (vv0,vv1,ld*ld*2,MPI_COMPLEX16,mpi_sum,comm_grid,ierr)
!$OMP end single

          call zhegv &
          (1,'V','U',mm,vv1(1,1,2),ld,vv1(1,1,1),ld,W,work,999,rwork,ierr)

!$OMP single
          E1(1:nn)=E(1:nn)
          E(1:nn) =W(1:nn)
!$OMP end single

          call zgemm('N','N',ML0,nn,mm,one,vt,ML0,vv1(1,1,2),ld,zero,wt,ML0)

!$OMP parallel
          do n=1,nn
!$OMP do
             do i=n1,n2
                v(i,n,1) = wt(i,n)
             end do
!$OMP end do
          end do
!$OMP end parallel

          call zgemm('N','N',ML0,nn,mm,one,ut,ML0,vv1(1,1,2),ld,zero,wt,ML0)

!$OMP parallel
          do n=1,nn
!$OMP do
             do i=n1,n2
                hv(i,n,1)=wt(i,n)
             end do
!$OMP end do
          end do
!$OMP end parallel

          call zgemm('N','N',ML0,nn,mm-nn,one,vt(n1,nn+1),ML0 &
                    ,vv1(nn+1,nn+1,2),ld,zero,wt,ML0)

!$OMP parallel
          do n=1,nn
!$OMP do
             do i=n1,n2
                v(i,n,2)=wt(i,n)
             end do
!$OMP end do
          end do
!$OMP end parallel

          call zgemm('N','N',ML0,nn,mm-nn,one,ut(n1,nn+1),ML0 &
                    ,vv1(nn+1,nn+1,2),ld,zero,wt,ML0)
!$OMP parallel
          do n=1,nn
!$OMP do
             do i=n1,n2
                hv(i,n,2)=wt(i,n)
             end do
!$OMP end do
          end do

          do n=1,nn
!$OMP do
             do i=n1,n2
                v(i,n,4) = -( hv(i,n,1) - W(n)*v(i,n,1) )
             end do
!$OMP end do
!$OMP single
             c=0.d0
!$OMP end single
!$OMP do reduction(+:c)
             do i=n1,n2
                c=c+abs(v(i,n,4))**2
             end do
!$OMP end do
!$OMP single
             sb(n)=c*dV
!$OMP end single
          end do

!$OMP single
          call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)
!$OMP end single

!$OMP end parallel

       end do ! icg

!$OMP parallel

!$OMP single
       esp(ns:ne)=E(1:nn)
!$OMP end single

       do n=1,nn
!$OMP do
          do i=n1,n2
             unk(i,ns+n-1) = v(i,n,1)
          end do
!$OMP end do
       end do

!$OMP end parallel

    end do  ! band-loop

    deallocate( vv1, vv0 )
    deallocate( wt )
    deallocate( ut, vt )
    deallocate( v, hv )
    deallocate( E,E1 )
    deallocate( sb,rb )

    return

  END SUBROUTINE conjugate_gradient
#endif

END MODULE cg_module
