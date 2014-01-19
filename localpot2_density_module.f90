MODULE localpot2_density_module

  use localpot2_variables
  use wf_module
  use array_bound_module
  use rgrid_module
  use electron_module

  implicit none

  PRIVATE
  PUBLIC :: localpot2_density

CONTAINS

  SUBROUTINE localpot2_density(m1,m2,m3,vout)
    implicit none
    integer,intent(IN)  :: m1,m2,m3
    real(8),intent(OUT) :: vout(0:m1-1,0:m2-1,0:m3-1)
    integer :: ml,mb,mk,ms,ML1,ML2,ML3
    integer :: s,k,n,ic1,ic2,ic3,id1,id2,id3,jd1,jd2,jd3
    integer :: itp1,itp2,itp3,i,jc1,jc2,jc3
    integer,allocatable :: LLL(:,:,:)
    complex(8) :: v

    vout(:,:,:)=0.0d0

    ML1=Ngrid(1)
    ML2=Ngrid(2)
    ML3=Ngrid(3)

    allocate( LLL(0:ML1-1,0:ML2-1,0:ML3-1) ) ; LLL=0

    i=0
    do ic3=0,Ngrid(3)-1
    do ic2=0,Ngrid(2)-1
    do ic1=0,Ngrid(1)-1
       i=i+1
       LLL(ic1,ic2,ic3) = i
    end do
    end do
    end do


    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=1,Nband

       do ic3=Igrid(1,3),Igrid(2,3)
       do jd3=0,Ndens_loc-1
          id3=ic3*Ndens_loc+jd3

          do ic2=Igrid(1,2),Igrid(2,2)
          do jd2=0,Ndens_loc-1
             id2=ic2*Ndens_loc+jd2

             do ic1=Igrid(1,1),Igrid(2,1)
             do jd1=0,Ndens_loc-1
                id1=ic1*Ndens_loc+jd1

                v=(0.0d0,0.0d0)
                do itp3=nitp_0,nitp_1
                do itp2=nitp_0,nitp_1
                do itp1=nitp_0,nitp_1

                   jc1 = mod(ic1+itp1+ML1,ML1)
                   jc2 = mod(ic2+itp2+ML2,ML2)
                   jc3 = mod(ic3+itp3+ML3,ML3)
                   i   = LLL(jc1,jc2,jc3)

                   v = v + Clag1(itp1,jd1)*Clag2(itp2,jd2)*Clag3(itp3,jd3) &
                           *unk(i,n,k,s)

                end do
                end do
                end do

                vout(id1,id2,id3) = vout(id1,id2,id3) + occ(n,k,s)*abs(v)**2

             end do
             end do

          end do
          end do

       end do
       end do

    end do
    end do
    end do

    deallocate( LLL )

    rewind 10
    do id1=0,m1-1
       write(10,*) id1,vout(id1,0,0)
    end do

    write(*,*) minval(vout),maxval(vout),count(vout/=0.0d0)
    write(*,*) sum(vout)*dV*Ngrid(0)/size(vout)
    write(*,*) sum(vout(0:m1-2:2,0:m2-2:2,0:m3-2:2))*dV

  END SUBROUTINE localpot2_density

END MODULE localpot2_density_module
