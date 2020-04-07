MODULE hamiltonian_matrix_module

  use hsort_module
  use kinetic_variables, only: Md
  use kinetic_sol_module
  use construct_matrix_ps_nloc2_module
  use localpot_module
  use symmetry_module
  use rgrid_module
  use bz_module, only: kbb

  implicit none

  PRIVATE
  PUBLIC :: construct_hamiltonian_matrix

#ifdef _DRSDFT_
  real(8),parameter :: zero=0.0d0
  real(8),allocatable :: Hmat(:,:)
  real(8),allocatable :: w1(:,:), w2(:,:)
#else
  complex(8),parameter :: zero=(0.0d0,0.0d0)
  complex(8),allocatable :: Hmat(:,:)
  complex(8),allocatable :: w1(:,:), w2(:,:)
#endif

  integer :: ML

  real(8),allocatable :: SymMat(:,:)

CONTAINS


  SUBROUTINE construct_hamiltonian_matrix( ML_in )
    implicit none
    integer,intent(IN) :: ML_in
    integer :: i,n,k,s,isym,i1,i2,i3,j,i0,j1,j2,j3
    integer :: k1,k2,k3,nkg
    real(8) :: c,r0(3),r1(3),r,tmp(3),err
    integer,allocatable :: sym_mat(:,:,:)

    write(*,'(a50," construct_hamiltonian_matrix")') repeat("-",50)

    ML = ML_in

    c = 1.0d0/nsym

    allocate( Hmat(ML,ML)   ) ; Hmat=zero
    allocate( SymMat(ML,ML) ) ; SymMat=0.0d0
    allocate( w1(ML,ML)     ) ; w1=zero
    allocate( w2(ML,ML)     ) ; w2=zero

    do s=1,1
    do k=1,1

       Hmat(:,:) = zero

       call construct_matrix_kinetic_sol( k, ML, Hmat )
!       call construct_matrix_localpot( s, ML, Hmat )
!       call construct_matrix_ps_nloc2( k, ML, Hmat )

       n=count(abs(Hmat)>1.d-10)
       write(*,*) n,dble(n)/dble(ML)
!       do i=1,ML
!          write(*,'(1x,2i5,2f20.15)') i,count(abs(Hmat(i,:))>1.d-10) &
!          ,Hmat(i,i)
!       end do
       !call check_hamil
!       i1=Ngrid(1)/2
!       i2=Ngrid(2)/2
!       i3=Ngrid(3)/2
!       i0=1+i1+i2*Ngrid(1)+i3*Ngrid(1)*Ngrid(2)
!       write(*,*) i0,i1,i2,i3
!       j=0
!       do i3=0,Ngrid(3)-1
!       do i2=0,Ngrid(2)-1
!       do i1=0,Ngrid(1)-1
!          i=1+i1+i2*Ngrid(1)+i3*Ngrid(1)*Ngrid(2)
!          if ( abs(Hmat(i0,i)) > 1.d-10 ) then
!             j=j+1
!             w1(1,j)=i1
!             w1(2,j)=i2
!             w1(3,j)=i3
!             w1(4,j)=Hmat(i0,i)
!             if ( i == i0 ) then
!                r0(1)=i1
!                r0(2)=i2
!                r0(3)=i3
!             end if
!          end if
!       end do
!       end do
!       end do
!       do i=1,j
!          r1(1:3)=real(w1(1:3,i))
!          r=sqrt(sum((r1-r0)**2))
!          write(*,'(1x,4i5,2f20.15)') i,nint(real(w1(1:3,i))),real(w1(4,i)),r
!       end do

       if ( isymmetry /= 0 ) then
          call get_mat_symmetry( sym_mat )
          w2(:,:)=zero
          nkg=0
          do isym=1,nsym
             tmp(:)=matmul(kbb(:,k),sym_mat(:,:,isym))
             !write(*,'(1x,i4,3f8.3)') isym,tmp
             err=abs(sum((tmp-kbb(:,k))**2))
             if ( err > 1.d-8 ) cycle
             nkg=nkg+1
             write(*,'(1x,"isym/nsym, nkg=",i3," /",i3,i4)') isym,nsym,nkg
             call construct_matrix_symmetry( isym, ML, SymMat )
             w1(:,:) = matmul( Hmat, transpose(SymMat) )
             w2(:,:) = w2(:,:) + matmul( SymMat, w1 )
          end do
          c=1.0d0/nkg
          Hmat(:,:) = c*w2(:,:)
       end if
       n=count(abs(Hmat)>1.d-10)
       write(*,*) n,dble(n)/dble(ML)
!       do i=1,ML
!          write(*,'(1x,2i5,2f20.15)') i,count(abs(Hmat(i,:))>1.d-10) &
!          ,Hmat(i,i)
!       end do
       !call check_hamil(10)
!       j=0
!       do i3=0,Ngrid(3)-1
!       do i2=0,Ngrid(2)-1
!       do i1=0,Ngrid(1)-1
!          i=1+i1+i2*Ngrid(1)+i3*Ngrid(1)*Ngrid(2)
!          if ( abs(Hmat(i0,i)) > 1.d-10 ) then
!             j=j+1
!             w1(1,j)=i1
!             w1(2,j)=i2
!             w1(3,j)=i3
!             w1(4,j)=Hmat(i0,i)
!          end if
!       end do
!       end do
!       end do
!       do i=1,j
!          r1(1:3)=real(w1(1:3,i))
!          r=sqrt(sum((r1-r0)**2))
!          write(*,'(1x,4i5,2f20.15)') i,nint(real(w1(1:3,i))),real(w1(4,i)),r
!       end do

!       call construct_matrix_kinetic_sol( k, ML, Hmat )
       call construct_matrix_localpot( s, ML, Hmat )
       call construct_matrix_ps_nloc2( k, ML, Hmat )

       call diag_Hmat

    end do ! k
    end do ! s

    deallocate( w2 )
    deallocate( w1 )
    deallocate( SymMat )
    deallocate( Hmat )

  END SUBROUTINE construct_hamiltonian_matrix


  SUBROUTINE diag_Hmat
    implicit none
    integer :: LWORK,LRWORK,ierr,i
    real(8),allocatable :: rwork(:),eigval(:)
    complex(8),allocatable :: work(:)
    allocate( eigval(ML) )
#ifdef _DRSDFT_
    LWORK=3*ML-1
    allocate( rwork(LWORK) )
    call dsyev('N','L',ML,Hmat,ML,eigval,rwork,LWORK,rwork,ierr)
    deallocate( rwork )
#else
    LWORK=2*ML-1
    LRWORK=3*ML-2
    allocate( work(LWORK),rwork(LRWORK) )
    call zheev('N','L',ML,Hmat,ML,eigval,work,LWORK,rwork,ierr)
    deallocate( rwork,work )
#endif
    do i=1,10
       write(*,'(1x,i5,f20.15)') i,eigval(i)
    end do
    do i=1,3
       write(*,'(1x,5x,".")')
    end do
    do i=ML-9,ML
       write(*,'(1x,i5,f20.15)') i,eigval(i)
    end do
    deallocate( eigval )
  END SUBROUTINE diag_Hmat


  SUBROUTINE check_hamil( unit )
    implicit none
    integer,optional,intent(IN) :: unit
    integer :: i1,i2,i3,j1,j2,j3,k1,k2,k3,i,j,i0
    integer,allocatable :: icheck(:,:,:),itmp(:,:),jtmp(:)
    real(8),allocatable :: rcheck(:,:,:),rtmp(:)
    if ( present(unit) ) rewind unit
    allocate( icheck(-2*Md:2*Md,-2*Md:2*Md,-2*Md:2*Md) ) ; icheck=0
    allocate( rcheck(-2*Md:2*Md,-2*Md:2*Md,-2*Md:2*Md) ) ; rcheck=0.0d0
    do i3=0,Ngrid(3)-1
    do i2=0,Ngrid(2)-1
    do i1=0,Ngrid(1)-1
       i=1+i1+i2*Ngrid(1)+i3*Ngrid(1)*Ngrid(2)
       if ( i /= 1 ) cycle
       write(*,*) "grid",i1,i2,i3,i
       i0=0
       do j3=0,Ngrid(3)-1
       do j2=0,Ngrid(2)-1
       do j1=0,Ngrid(1)-1
          j=1+j1+j2*Ngrid(1)+j3*Ngrid(1)*Ngrid(2)
          if ( abs(Hmat(i,j))>1.d-10 ) then
             i0=i0+1
             k1=j1-i1 ; if ( k1 >= Ngrid(1)/2 ) k1=k1-Ngrid(1)
             k2=j2-i2 ; if ( k2 >= Ngrid(2)/2 ) k2=k2-Ngrid(2)
             k3=j3-i3 ; if ( k3 >= Ngrid(3)/2 ) k3=k3-Ngrid(3)
             write(*,'(1x,7i4,2x,2f20.15)') i0,k1,k2,k3,j1-i1,j2-i2,j3-i3 &
                  ,Hmat(i,j),Hmat(i,j)*Hgrid(1)**2
!             if ( present(unit) ) then
!                write(unit,'(1x,3i4,f22.16)') k1,k2,k3,Hmat(i,j)*Hgrid(1)**2
!             end if
             icheck(k1,k2,k3)=icheck(k1,k2,k3)+1
             rcheck(k1,k2,k3)=Hmat(i,j)*Hgrid(1)**2
          end if
       end do
       end do
       end do
    end do
    end do
    end do
    write(*,*) maxval(icheck),sum(icheck),count(icheck/=0)
    i=count(icheck/=0)

    if ( present(unit) ) write(unit,*) i

    allocate( itmp(3,i), jtmp(i), rtmp(i) )
    itmp=0 ; jtmp=0 ; rtmp=0.0d0
    i=0
    do i3=-2*Md,2*Md
    do i2=-2*Md,2*Md
    do i1=-2*Md,2*Md
       if ( icheck(i1,i2,i3) /= 0 ) then
          i=i+1
          itmp(1,i)=i1
          itmp(2,i)=i2
          itmp(3,i)=i3
          rtmp(i)=i1*i1+i2*i2+i3*i3
          if ( present(unit) ) then
             write(unit,'(1x,3i4,f22.16)') i1,i2,i3,rcheck(i1,i2,i3)
          end if
       end if
    end do
    end do
    end do
    call indexx( size(rtmp), rtmp, jtmp )
    do i=1,size(rtmp)
       j=jtmp(i)
       i1=itmp(1,j)
       i2=itmp(2,j)
       i3=itmp(3,j)
       write(*,'(1x,4i4,2f20.10)') i, itmp(:,j), rtmp(j), rcheck(i1,i2,i3)
    end do
    deallocate( itmp, jtmp, rtmp )
    deallocate( rcheck )
    deallocate( icheck )
  END SUBROUTINE check_hamil


END MODULE hamiltonian_matrix_module
