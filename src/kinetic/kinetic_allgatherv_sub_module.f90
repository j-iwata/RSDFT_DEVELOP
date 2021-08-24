module kinetic_allgatherv_sub_module

  use kinetic_variables, only: coef_lap0, coef_lap, zcoef_kin, coef_nab &
       ,flag_nab, flag_n12, flag_n23, flag_n31, const_k2, ggg, Md
  use grid_module, only: construct_map_3d_to_1d_grid

  implicit none

  private
  public :: kinetic_allgatherv_sub

contains


  subroutine kinetic_allgatherv_sub( Igrid,comm,list_vec,kin_coef )
    implicit none
    integer,intent(in) :: Igrid(2,0:3),comm
    integer,allocatable,intent(out) :: list_vec(:,:)
#ifdef _DRSDFT_
    real(8),allocatable,intent(out) :: kin_coef(:,:,:)
    real(8) :: d
    real(8),allocatable :: Hrow(:)
    real(8),parameter :: zero=0.0d0
#else
    complex(8),allocatable,intent(out) :: kin_coef(:,:,:)
    complex(8) :: d
    complex(8),allocatable :: Hrow(:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
#endif
    integer :: k,i,i0,i1,i2,i3,m,n,j,j1,j2,j3,ierr,Ngrid(0:3)
    integer :: ML1,ML2,ML3,ML,MK,n0,n1,n2,itmp(3,2),jtmp(3,2)
    integer,allocatable :: LLL(:,:,:)
    include 'mpif.h'

    itmp(1:3,1) = Igrid(1,1:3)
    itmp(1:3,2) = Igrid(2,1:3)
    call MPI_Allreduce( itmp(1,1),jtmp(1,1),3,MPI_INTEGER,MPI_MIN,comm,ierr )
    call MPI_Allreduce( itmp(1,2),jtmp(1,2),3,MPI_INTEGER,MPI_MAX,comm,ierr )
    Ngrid(1:3) = jtmp(1:3,2) - jtmp(1:3,1) + 1
    Ngrid(0) = product( Ngrid(1:3) )

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    MK  = size( const_k2 ) - 1

    call construct_map_3d_to_1d_grid( LLL )

    allocate( Hrow(ML) ); Hrow=zero

    n1 = Igrid(2,1) - Igrid(1,1) + 1
    n2 = Igrid(2,2) - Igrid(1,2) + 1
    i0 = Igrid(1,0)
    n0 = Igrid(2,0) - i0 + 1

    do k=1,MK

    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)

       i=1+i1-Igrid(1,1)+(i2-Igrid(1,2))*n1+(i3-Igrid(1,3))*n1*n2

       Hrow(:)=zero

       d = coef_lap0 + const_k2(k)

       Hrow(i+i0-1) = Hrow(i+i0-1) + d

       if ( flag_nab ) then

          do m=1,Md

             j1 = mod(i1+m+ML1,ML1)
             j  = LLL(j1,i2,i3)
             Hrow(j) = Hrow(j) + zcoef_kin(1,m,k)

             j1 = mod(i1-m+ML1,ML1)
             j  = LLL(j1,i2,i3)
             Hrow(j) = Hrow(j) + conjg( zcoef_kin(1,m,k) )

             j2 = mod(i2+m+ML2,ML2)
             j  = LLL(i1,j2,i3)
             Hrow(j) = Hrow(j) + zcoef_kin(2,m,k)

             j2 = mod(i2-m+ML2,ML2)
             j  = LLL(i1,j2,i3)
             Hrow(j) = Hrow(j) + conjg( zcoef_kin(2,m,k) )

             j3 = mod(i3+m+ML3,ML3)
             j  = LLL(i1,i2,j3)
             Hrow(j) = Hrow(j) + zcoef_kin(3,m,k)

             j3 = mod(i3-m+ML3,ML3)
             j  = LLL(i1,i2,j3)
             Hrow(j) = Hrow(j) + conjg( zcoef_kin(3,m,k) )

          end do ! m

       else

          do m=1,Md

             j1 = mod(i1+m+ML1,ML1)
             j  = LLL(j1,i2,i3)
             Hrow(j) = Hrow(j) + coef_lap(1,m)

             j1 = mod(i1-m+ML1,ML1)
             j  = LLL(j1,i2,i3)
             Hrow(j) = Hrow(j) + coef_lap(1,m)

             j2 = mod(i2+m+ML2,ML2)
             j  = LLL(i1,j2,i3)
             Hrow(j) = Hrow(j) + coef_lap(2,m)

             j2 = mod(i2-m+ML2,ML2)
             j  = LLL(i1,j2,i3)
             Hrow(j) = Hrow(j) + coef_lap(2,m)

             j3 = mod(i3+m+ML3,ML3)
             j  = LLL(i1,i2,j3)
             Hrow(j) = Hrow(j) + coef_lap(3,m)

             j3 = mod(i3-m+ML3,ML3)
             j  = LLL(i1,i2,j3)
             Hrow(j) = Hrow(j) + coef_lap(3,m)

          end do ! m

       end if

       if ( flag_n12 .or. flag_n23 .or. flag_n31 ) then

          if ( flag_n12 ) then

             do n=1,Md
             do m=1,Md

                d = -ggg(4)*coef_nab(1,m)*coef_nab(2,n)

                j1 = mod(i1+m+ML1,ML1)
                j2 = mod(i2+n+ML2,ML2)
                j  = LLL(j1,j2,i3)
                Hrow(j) = Hrow(j) + d

                j1 = mod(i1-m+ML1,ML1)
                j2 = mod(i2+n+ML2,ML2)
                j  = LLL(j1,j2,i3)
                Hrow(j) = Hrow(j) - d

                j1 = mod(i1+m+ML1,ML1)
                j2 = mod(i2-n+ML2,ML2)
                j  = LLL(j1,j2,i3)
                Hrow(j) = Hrow(j) - d

                j1 = mod(i1-m+ML1,ML1)
                j2 = mod(i2-n+ML2,ML2)
                j  = LLL(j1,j2,i3)
                Hrow(j) = Hrow(j) + d

             end do ! m
             end do ! n

          end if

          if ( flag_n23 ) then

             do n=1,Md
             do m=1,Md

                d = -ggg(5)*coef_nab(2,m)*coef_nab(3,n)

                j2 = mod(i2+m+ML2,ML2)
                j3 = mod(i3+n+ML3,ML3)
                j  = LLL(i1,j2,j3)
                Hrow(j) = Hrow(j) + d

                j2 = mod(i2-m+ML2,ML2)
                j3 = mod(i3+n+ML3,ML3)
                j  = LLL(i1,j2,j3)
                Hrow(j) = Hrow(j) - d

                j2 = mod(i2+m+ML2,ML2)
                j3 = mod(i3-n+ML3,ML3)
                j  = LLL(i1,j2,j3)
                Hrow(j) = Hrow(j) - d

                j2 = mod(i2-m+ML2,ML2)
                j3 = mod(i3-n+ML3,ML3)
                j  = LLL(i1,j2,j3)
                Hrow(j) = Hrow(j) + d

             end do ! m
             end do ! n

          end if

          if ( flag_n31 ) then

             do n=1,Md
             do m=1,Md

                d = -ggg(6)*coef_nab(1,m)*coef_nab(3,n)

                j1 = mod(i1+m+ML1,ML1)
                j3 = mod(i3+n+ML3,ML3)
                j  = LLL(j1,i2,j3)
                Hrow(j) = Hrow(j) + d

                j1 = mod(i1-m+ML1,ML1)
                j3 = mod(i3+n+ML3,ML3)
                j  = LLL(j1,i2,j3)
                Hrow(j) = Hrow(j) - d

                j1 = mod(i1+m+ML1,ML1)
                j3 = mod(i3-n+ML3,ML3)
                j  = LLL(j1,i2,j3)
                Hrow(j) = Hrow(j) - d

                j1 = mod(i1-m+ML1,ML1)
                j3 = mod(i3-n+ML3,ML3)
                j  = LLL(j1,i2,j3)
                Hrow(j) = Hrow(j) + d

             end do ! m
             end do ! n

          end if

       end if

       if ( .not.allocated(list_vec) ) then
          n=count(Hrow/=zero)
          allocate( list_vec(n,n0) ); list_vec=0
          allocate( kin_coef(n,n0,MK) ); kin_coef=zero
       end if

       n=0
       do j=1,ML
          if ( Hrow(j) /= zero ) then
             n=n+1
             list_vec(n,i)=j
             kin_coef(n,i,k)=Hrow(j)
          end if
       end do

    end do ! i3
    end do ! i2
    end do ! i1

    end do ! k

    deallocate( Hrow )
    deallocate( LLL )

  end subroutine kinetic_allgatherv_sub


end module kinetic_allgatherv_sub_module
