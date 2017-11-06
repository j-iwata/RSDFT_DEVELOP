MODULE kinetic_sym_ini_module

  use symmetry_module, only: basis_conversion_symmetry, calc_mat_bb_symmetry
  use kinetic_variables, only: Md,coef_lap0,coef_lap,ggg,coef_nab,zcoef_kin,const_k2
  use kinetic_sym_module, only: init2
  use bz_module, only: Nbzsm, kbb
  use lattice_module, only: get_inverse_lattice

  implicit none

  PRIVATE
  PUBLIC :: init_kinetic_sym
  PUBLIC :: construct_kinetic_sym
  PUBLIC :: get_mat_kinetic_sym_ini

  integer,allocatable :: sym_mat(:,:,:)
  real(8) :: basis(3,3)

  integer,allocatable :: kgroup_mat(:,:,:,:)
  integer,allocatable :: kgroup_g(:)

CONTAINS


  SUBROUTINE init_kinetic_sym( indx, aa )

    implicit none
    character(*),intent(IN) :: indx
    real(8),intent(IN) :: aa(3,3)
    integer :: i,k

    call write_border( 0, " init_kinetic_sym(start)" )

    select case( indx )
    case( "HEXAGONAL" )
       call set_hexagonal_sym( sym_mat, basis )
    case( "FCC" )
       call set_fcc_sym( sym_mat, basis )
!   case( "BCC" )
!   case( "CUBIC" )
    case default
       write(*,*) "indx= ",indx
       call write_string( "this lattice is undefined(init_kinetic_sym)" )
       return
    end select

    call basis_conversion_symmetry( basis, aa, sym_mat )

    call construct_kgroup
    do k=1,Nbzsm
       call construct_kinetic_sym( kgroup_mat(:,:,1:kgroup_g(k),k), k )
    end do

    call write_border( 0, " init_kinetic_sym(end)" )

  END SUBROUTINE init_kinetic_sym


  SUBROUTINE construct_kgroup

    implicit none
    integer :: k,j,i
    real(8) :: tmp(3),err
    logical :: disp

    call check_disp_switch( disp, 0 )

    allocate( kgroup_mat(3,3,size(sym_mat,3),Nbzsm) ) ; kgroup_mat=0
    allocate( kgroup_g(Nbzsm) ) ; kgroup_g=0

    do k=1,Nbzsm
       j=0
       do i=1,size(sym_mat,3)
          tmp(:) = matmul( kbb(:,k), sym_mat(:,:,i) )
          if ( any(abs(tmp)>0.5d0+1.d-10) ) write(*,*) i,k,tmp
          err = abs(sum((tmp-kbb(:,k))**2))
          if ( err < 1.d-8 ) then
             if ( disp ) write(*,'(1x,2i4,3f10.5,2x,3f10.5)') k,i,kbb(:,k),tmp(:)
             j=j+1
             kgroup_mat(:,:,j,k)=sym_mat(:,:,i)
          end if
       end do ! i
       kgroup_g(k)=j
       if ( disp ) write(*,*) "# of kgroup elements =",j,k
       call chk_grp( kgroup_mat(:,:,1:j,k), basis )
    end do ! k

  END SUBROUTINE construct_kgroup


  SUBROUTINE chk_grp( rga, aa )
    implicit none
    integer,intent(IN) :: rga(:,:,:)
    real(8),intent(IN) :: aa(3,3)
    real(8) :: aa_inv(3,3),tmp0(3,3),tmp1(3,3),err
    real(8),allocatable :: rtmp(:,:,:)
    integer :: n,ns,m
    ns=size(rga,3)
    allocate( rtmp(3,3,ns) ) ; rtmp=0.0d0
    call get_inverse_lattice( aa, aa_inv )
    do n=1,ns
       tmp0(:,:) = rga(:,:,n)
       tmp1(:,:) = matmul( tmp0, aa_inv )
       rtmp(:,:,n) = matmul( aa, tmp1 )
    end do
    do n=1,ns
       tmp0(:,:) = transpose( rtmp(:,:,n) )
       tmp1(:,:) = matmul( rtmp(:,:,n), tmp0 )
       do m=1,ns
          err=sum((tmp0-rtmp(:,:,m))**2)
          if ( err < 1.d-8 ) exit
       end do
       if ( m > ns ) call stop_program("stop@chk_grp in kinetic_sym_ini_module.f90")
!       write(*,*) n,sum(abs(tmp1)),m,ns
    end do
    deallocate( rtmp )
  END SUBROUTINE chk_grp



  SUBROUTINE set_fcc_sym( mat, vec )
    implicit none
    integer,allocatable :: mat(:,:,:)
    real(8),intent(OUT) :: vec(3,3)
    integer :: tmp(3,3),i,j
    allocate( mat(3,3,48) ) ; mat=0
    mat(:,:, 1)=reshape( (/   1,  0,  0,   0,  1,  0,   0,  0,  1 /),(/3,3/) ) 
    mat(:,:, 2)=reshape( (/  -1, -1, -1,   0,  0,  1,   0,  1,  0 /),(/3,3/) ) 
    mat(:,:, 3)=reshape( (/   0,  0,  1,  -1, -1, -1,   1,  0,  0 /),(/3,3/) ) 
    mat(:,:, 4)=reshape( (/   0,  1,  0,   1,  0,  0,  -1, -1, -1 /),(/3,3/) ) 
    mat(:,:, 5)=reshape( (/   0,  1,  0,   0,  0,  1,   1,  0,  0 /),(/3,3/) )
    mat(:,:, 6)=reshape( (/  -1, -1, -1,   1,  0,  0,   0,  0,  1 /),(/3,3/) )
    mat(:,:, 7)=reshape( (/   1,  0,  0,  -1, -1, -1,   0,  1,  0 /),(/3,3/) )
    mat(:,:, 8)=reshape( (/   0,  0,  1,   0,  1,  0,  -1, -1, -1 /),(/3,3/) )
    mat(:,:, 9)=reshape( (/   0,  0,  1,   1,  0,  0,   0,  1,  0 /),(/3,3/) )
    mat(:,:,10)=reshape( (/  -1, -1, -1,   0,  1,  0,   1,  0,  0 /),(/3,3/) )
    mat(:,:,11)=reshape( (/   0,  1,  0,  -1, -1, -1,   0,  0,  1 /),(/3,3/) )
    mat(:,:,12)=reshape( (/   1,  0,  0,   0,  0,  1,  -1, -1, -1 /),(/3,3/) )
    mat(:,:,13)=reshape( (/   0, -1,  0,  -1,  0,  0,   0,  0, -1 /),(/3,3/) )
    mat(:,:,14)=reshape( (/   1,  1,  1,   0,  0, -1,  -1,  0,  0 /),(/3,3/) )
    mat(:,:,15)=reshape( (/   0,  0, -1,   1,  1,  1,   0, -1,  0 /),(/3,3/) )
    mat(:,:,16)=reshape( (/  -1,  0,  0,   0, -1,  0,   1,  1,  1 /),(/3,3/) )
    mat(:,:,17)=reshape( (/  -1,  0,  0,   0,  0, -1,   0, -1,  0 /),(/3,3/) )
    mat(:,:,18)=reshape( (/   1,  1,  1,   0, -1,  0,   0,  0, -1 /),(/3,3/) )
    mat(:,:,19)=reshape( (/   0, -1,  0,   1,  1,  1,  -1,  0,  0 /),(/3,3/) )
    mat(:,:,20)=reshape( (/   0,  0, -1,  -1,  0,  0,   1,  1,  1 /),(/3,3/) )
    mat(:,:,21)=reshape( (/   0,  0, -1,   0, -1,  0,  -1,  0,  0 /),(/3,3/) )
    mat(:,:,22)=reshape( (/   1,  1,  1,  -1,  0,  0,   0, -1,  0 /),(/3,3/) )
    mat(:,:,23)=reshape( (/  -1,  0,  0,   1,  1,  1,   0,  0, -1 /),(/3,3/) )
    mat(:,:,24)=reshape( (/   0, -1,  0,   0,  0, -1,   1,  1,  1 /),(/3,3/) )
    mat(:,:,25)=reshape( (/  -1,  0,  0,   0, -1,  0,   0,  0, -1 /),(/3,3/) )
    mat(:,:,26)=reshape( (/   1,  1,  1,   0,  0, -1,   0, -1,  0 /),(/3,3/) )
    mat(:,:,27)=reshape( (/   0,  0, -1,   1,  1,  1,  -1,  0,  0 /),(/3,3/) )
    mat(:,:,28)=reshape( (/   0, -1,  0,  -1,  0,  0,   1,  1,  1 /),(/3,3/) )
    mat(:,:,29)=reshape( (/   0, -1,  0,   0,  0, -1,  -1,  0,  0 /),(/3,3/) )
    mat(:,:,30)=reshape( (/   1,  1,  1,  -1,  0,  0,   0,  0, -1 /),(/3,3/) )
    mat(:,:,31)=reshape( (/  -1,  0,  0,   1,  1,  1,   0, -1,  0 /),(/3,3/) )
    mat(:,:,32)=reshape( (/   0,  0, -1,   0, -1,  0,   1,  1,  1 /),(/3,3/) )
    mat(:,:,33)=reshape( (/   0,  0, -1,  -1,  0,  0,   0, -1,  0 /),(/3,3/) )
    mat(:,:,34)=reshape( (/   1,  1,  1,   0, -1,  0,  -1,  0,  0 /),(/3,3/) )
    mat(:,:,35)=reshape( (/   0, -1,  0,   1,  1,  1,   0,  0, -1 /),(/3,3/) )
    mat(:,:,36)=reshape( (/  -1,  0,  0,   0,  0, -1,   1,  1,  1 /),(/3,3/) )
    mat(:,:,37)=reshape( (/   0,  1,  0,   1,  0,  0,   0,  0,  1 /),(/3,3/) )
    mat(:,:,38)=reshape( (/  -1, -1, -1,   0,  0,  1,   1,  0,  0 /),(/3,3/) )
    mat(:,:,39)=reshape( (/   0,  0,  1,  -1, -1, -1,   0,  1,  0 /),(/3,3/) )
    mat(:,:,40)=reshape( (/   1,  0,  0,   0,  1,  0,  -1, -1, -1 /),(/3,3/) )
    mat(:,:,41)=reshape( (/   1,  0,  0,   0,  0,  1,   0,  1,  0 /),(/3,3/) )
    mat(:,:,42)=reshape( (/  -1, -1, -1,   0,  1,  0,   0,  0,  1 /),(/3,3/) )
    mat(:,:,43)=reshape( (/   0,  1,  0,  -1, -1, -1,   1,  0,  0 /),(/3,3/) )
    mat(:,:,44)=reshape( (/   0,  0,  1,   1,  0,  0,  -1, -1, -1 /),(/3,3/) )
    mat(:,:,45)=reshape( (/   0,  0,  1,   0,  1,  0,   1,  0,  0 /),(/3,3/) )
    mat(:,:,46)=reshape( (/  -1, -1, -1,   1,  0,  0,   0,  1,  0 /),(/3,3/) )
    mat(:,:,47)=reshape( (/   1,  0,  0,  -1, -1, -1,   0,  0,  1 /),(/3,3/) )
    mat(:,:,48)=reshape( (/   0,  1,  0,   0,  0,  1,  -1, -1, -1 /),(/3,3/) )
    do i=1,size(mat,3)
       tmp(:,:)=mat(:,:,i)
       mat(:,:,i)=transpose( tmp(:,:) )
    end do
    vec(:,1)=(/ 0.0d0,  0.5d0,  0.5d0 /)
    vec(:,2)=(/ 0.5d0,  0.0d0,  0.5d0 /)
    vec(:,3)=(/ 0.5d0,  0.5d0,  0.0d0 /)
!    write(20,*) size(mat,3),1
!    do i=1,size(mat,3)
!       write(20,'(3(1x,3i3),2x,3i3)') (mat(j,:,i),j=1,3),0,0,0
!    end do
  END SUBROUTINE set_fcc_sym

  SUBROUTINE set_hexagonal_sym( mat, vec )
    implicit none
    integer,allocatable :: mat(:,:,:)
    real(8),intent(OUT) :: vec(3,3)
    integer :: tmp(3,3),i
    allocate( mat(3,3,24) ) ; mat=0
    mat(:,:, 1)=reshape( (/ 1, 0, 0,  0, 1, 0,  0, 0, 1 /),(/3,3/) ) ! (+a,+b,+c)
    mat(:,:, 2)=reshape( (/-1,-1, 0,  0, 1, 0,  0, 0, 1 /),(/3,3/) ) ! (-a-b,+b,+c)
    mat(:,:, 3)=reshape( (/-1, 0, 0,  1, 1, 0,  0, 0, 1 /),(/3,3/) ) ! (-a,+a+b,+c)
    mat(:,:, 4)=reshape( (/ 0,-1, 0,  1, 1, 0,  0, 0, 1 /),(/3,3/) ) ! (-b,+a+b,+c)
    mat(:,:, 5)=reshape( (/ 0,-1, 0, -1, 0, 0,  0, 0, 1 /),(/3,3/) ) ! (-b,-a,+c)
    mat(:,:, 6)=reshape( (/ 1, 1, 0, -1, 0, 0,  0, 0, 1 /),(/3,3/) ) ! (+a+b,-a,+c)
    mat(:,:, 7)=reshape( (/-1,-1, 0,  1, 0, 0,  0, 0, 1 /),(/3,3/) ) ! (-a-b,+a,+c)
    mat(:,:, 8)=reshape( (/ 0, 1, 0,  1, 0, 0,  0, 0, 1 /),(/3,3/) ) ! (+b,+a,+c)
    mat(:,:, 9)=reshape( (/-1, 0, 0,  0,-1, 0,  0, 0, 1 /),(/3,3/) ) ! (-a,-b,+c)
    mat(:,:,10)=reshape( (/ 1, 1, 0,  0,-1, 0,  0, 0, 1 /),(/3,3/) ) ! (+a+b,-b,+c)
    mat(:,:,11)=reshape( (/ 1, 0, 0, -1,-1, 0,  0, 0, 1 /),(/3,3/) ) ! (+a,-a-b,+c)
    mat(:,:,12)=reshape( (/ 0, 1, 0, -1,-1, 0,  0, 0, 1 /),(/3,3/) ) ! (+b,-a-b,+c)
    mat(:,:,13)=reshape( (/ 1, 0, 0,  0, 1, 0,  0, 0,-1 /),(/3,3/) ) ! (+a,+b,-c)
    mat(:,:,14)=reshape( (/-1,-1, 0,  0, 1, 0,  0, 0,-1 /),(/3,3/) ) ! (-a-b,+b,-c)
    mat(:,:,15)=reshape( (/-1, 0, 0,  1, 1, 0,  0, 0,-1 /),(/3,3/) ) ! (-a,+a+b,-c)
    mat(:,:,16)=reshape( (/ 0,-1, 0,  1, 1, 0,  0, 0,-1 /),(/3,3/) ) ! (-b,+a+b,-c)
    mat(:,:,17)=reshape( (/ 0,-1, 0, -1, 0, 0,  0, 0,-1 /),(/3,3/) ) ! (-b,-a,-c)
    mat(:,:,18)=reshape( (/ 1, 1, 0, -1, 0, 0,  0, 0,-1 /),(/3,3/) ) ! (+a+b,-a,-c)
    mat(:,:,19)=reshape( (/-1,-1, 0,  1, 0, 0,  0, 0,-1 /),(/3,3/) ) ! (-a-b,+a,-c)
    mat(:,:,20)=reshape( (/ 0, 1, 0,  1, 0, 0,  0, 0,-1 /),(/3,3/) ) ! (+b,+a,-c)
    mat(:,:,21)=reshape( (/-1, 0, 0,  0,-1, 0,  0, 0,-1 /),(/3,3/) ) ! (-a,-b,-c)
    mat(:,:,22)=reshape( (/ 1, 1, 0,  0,-1, 0,  0, 0,-1 /),(/3,3/) ) ! (+a+b,-b,-c)
    mat(:,:,23)=reshape( (/ 1, 0, 0, -1,-1, 0,  0, 0,-1 /),(/3,3/) ) ! (+a,-a-b,-c)
    mat(:,:,24)=reshape( (/ 0, 1, 0, -1,-1, 0,  0, 0,-1 /),(/3,3/) ) ! (+b,-a-b,-c)
    do i=1,size(mat,3)
       tmp(:,:)=mat(:,:,i)
       mat(:,:,i)=transpose( tmp(:,:) )
    end do
    vec(:,1)=(/ sqrt(3.0d0)*0.5d0,  0.5d0, 0.0d0 /)
    vec(:,2)=(/ sqrt(3.0d0)*0.5d0, -0.5d0, 0.0d0 /)
    vec(:,3)=(/ 0.0d0            ,  0.0d0, 1.0d0 /)
  END SUBROUTINE set_hexagonal_sym


  SUBROUTINE construct_kinetic_sym( sym_mat, k )

    implicit none
    integer,intent(IN) :: sym_mat(:,:,:), k
    integer :: i1,i2,i3,j1,j2,j3,j(3),Rj(3),m,n,isym,nsym
#ifdef _DRSDFT_
    real(8),allocatable :: tmp(:,:,:)
#else
    complex(8),allocatable :: tmp(:,:,:)
#endif
    real(8) :: d

    allocate( tmp(-2*Md:2*Md,-2*Md:2*Md,-2*Md:2*Md) ) ; tmp=(0.0d0,0.0d0)

    nsym=size( sym_mat, 3 )

    i1=0
    i2=0
    i3=0

    do isym=1,nsym
       j(1) = i1
       j(2) = i2
       j(3) = i3
       Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
       tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + coef_lap0 + const_k2(k)
    end do

    do m=1,Md

       do isym=1,nsym

#ifdef _DRSDFT_
          j(1) = i1 + m
          j(2) = i2
          j(3) = i3
          Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + coef_lap(1,m)

          j(1) = i1 - m
          j(2) = i2
          j(3) = i3
          Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + coef_lap(1,m)

          j(1) = i1
          j(2) = i2 + m
          j(3) = i3
          Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + coef_lap(2,m)

          j(1) = i1
          j(2) = i2 - m
          j(3) = i3
          Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + coef_lap(2,m)

          j(1) = i1
          j(2) = i2
          j(3) = i3 + m
          Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + coef_lap(3,m)

          j(1) = i1
          j(2) = i2
          j(3) = i3 - m
          Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + coef_lap(3,m)
#else
          j(1) = i1 + m
          j(2) = i2
          j(3) = i3
          Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + zcoef_kin(1,m,k)

          j(1) = i1 - m
          j(2) = i2
          j(3) = i3
          Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + conjg(zcoef_kin(1,m,k))

          j(1) = i1
          j(2) = i2 + m
          j(3) = i3
          Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + zcoef_kin(2,m,k)

          j(1) = i1
          j(2) = i2 - m
          j(3) = i3
          Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + conjg(zcoef_kin(2,m,k))

          j(1) = i1
          j(2) = i2
          j(3) = i3 + m
          Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + zcoef_kin(3,m,k)

          j(1) = i1
          j(2) = i2
          j(3) = i3 - m
          Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + conjg(zcoef_kin(3,m,k))
#endif

       end do ! isym

    end do ! m


    if ( ggg(4) /= 0.0d0 ) then

       do n=1,Md
       do m=1,Md

          d = -ggg(4)*coef_nab(1,m)*coef_nab(2,n)

          do isym=1,nsym

             j(1) = i1 + m
             j(2) = i2 + n
             j(3) = i3
             Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
             tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + d

             j(1) = i1 - m
             j(2) = i2 + n
             j(3) = i3
             Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
             tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) - d

             j(1) = i1 + m
             j(2) = i2 - n
             j(3) = i3
             Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
             tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) - d

             j(1) = i1 - m
             j(2) = i2 - n
             j(3) = i3
             Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
             tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + d

          end do ! isym

       end do ! m
       end do ! n

    end if ![ ggg(4) /= 0.0d0 ]


    if ( ggg(5) /= 0.0d0 ) then

       do n=1,Md
       do m=1,Md

          d = -ggg(5)*coef_nab(2,m)*coef_nab(3,n)

          do isym=1,nsym

             j(1) = i1
             j(2) = i2 + m
             j(3) = i3 + n
             Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
             tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + d

             j(1) = i1
             j(2) = i2 - m
             j(3) = i3 + n
             Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
             tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) - d

             j(1) = i1
             j(2) = i2 + m
             j(3) = i3 - n
             Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
             tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) - d

             j(1) = i1
             j(2) = i2 - m
             j(3) = i3 - n
             Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
             tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + d

          end do ! isym

       end do ! m
       end do ! n

    end if ![ ggg(5) /= 0.0d0 ]


    if ( ggg(6) /= 0.0d0 ) then

       do n=1,Md
       do m=1,Md

          d = -ggg(6)*coef_nab(1,m)*coef_nab(3,n)

          do isym=1,nsym

             j(1) = i1 + m
             j(2) = i2
             j(3) = i3 + n
             Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
             tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + d

             j(1) = i1 - m
             j(2) = i2
             j(3) = i3 + n
             Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
             tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) - d

             j(1) = i1 + m
             j(2) = i2
             j(3) = i3 - n
             Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
             tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) - d

             j(1) = i1 - m
             j(2) = i2
             j(3) = i3 - n
             Rj(:) = matmul( sym_mat(:,:,isym), j(:) )
             tmp(Rj(1),Rj(2),Rj(3)) = tmp(Rj(1),Rj(2),Rj(3)) + d

          end do ! isym

       end do ! m
       end do ! n

    end if ![ ggg(6) /= 0.0d0 ]

!    write(*,*) count(tmp/=0.0d0),count(abs(tmp)>1.d-8),nsym

    tmp=tmp/dble(nsym)

!    call init2( 2*Md, 1, k, tmp )
    call init2( 2*Md, Nbzsm, k, tmp )
 
    deallocate( tmp )

  END SUBROUTINE construct_kinetic_sym


  SUBROUTINE get_mat_kinetic_sym_ini( mata, matb )
    implicit none
    integer,optional,allocatable,intent(INOUT) :: mata(:,:,:)
    real(8),optional,allocatable,intent(INOUT) :: matb(:,:,:)
    integer :: nsym
    if ( .not.allocated(sym_mat) ) return
    nsym=size(sym_mat,3)
    if ( present(mata) ) then
       if ( .not.allocated(mata) ) allocate( mata(3,3,nsym) )
       mata=0
       mata=sym_mat
    end if
    if ( present(matb) ) then
       if ( .not.allocated(matb) ) allocate( matb(3,3,nsym) )
       matb=0.0d0
       call calc_mat_bb_symmetry( mata, matb )
    end if
  END SUBROUTINE get_mat_kinetic_sym_ini


END MODULE kinetic_sym_ini_module
