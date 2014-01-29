MODULE localpot2_density_module

  use localpot2_variables
  use wf_module
  use array_bound_module
  use rgrid_module
  use electron_module
  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: localpot2_density

CONTAINS

  SUBROUTINE localpot2_density(m1_0,m2_0,m3_0,vout)
    implicit none
    integer,intent(IN)  :: m1_0,m2_0,m3_0
    real(8),intent(OUT) :: vout(m1_0,m2_0,m3_0)
    integer :: ml,mb,mk,ms,ML1,ML2,ML3,m1,m2,m3
    integer :: s,k,n,ic1,ic2,ic3,id1,id2,id3,jd1,jd2,jd3
    integer :: itp1,itp2,itp3,i,jc1,jc2,jc3,ierr
    integer,allocatable :: LLL(:,:,:),LLLtmp(:,:,:)
#ifdef _DRSDFT_
    integer,parameter :: TYP=MPI_REAL8
    real(8) :: v
    real(8),allocatable :: utmp(:)
#else
    integer,parameter :: TYP=MPI_COMPLEX16
    complex(8) :: v
    complex(8),allocatable :: utmp(:)
#endif

    vout(:,:,:)=0.0d0

    m1=Ngrid_dense(1)
    m2=Ngrid_dense(2)
    m3=Ngrid_dense(3)

    ML1=Ngrid(1)
    ML2=Ngrid(2)
    ML3=Ngrid(3)

    allocate( LLL(0:ML1-1,0:ML2-1,0:ML3-1)    ) ; LLL=0
    allocate( LLLtmp(0:ML1-1,0:ML2-1,0:ML3-1) ) ; LLLtmp=0

    i=ML_0-1
    do ic3=Igrid(1,3),Igrid(2,3)
    do ic2=Igrid(1,2),Igrid(2,2)
    do ic1=Igrid(1,1),Igrid(2,1)
       i=i+1
       LLLtmp(ic1,ic2,ic3) = i
    end do
    end do
    end do
    call MPI_ALLREDUCE(LLLtmp,LLL,Ngrid(0),MPI_INTEGER,MPI_SUM,comm_grid,ierr)

    deallocate( LLLtmp )

    allocate( utmp(Ngrid(0)) ) ; utmp=0.0d0

    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=1,Nband

       call MPI_ALLGATHERV(unk(ML_0,n,k,s),ir_grid(myrank_g),TYP &
            ,utmp,ir_grid,id_grid,TYP,comm_grid,ierr)

       do ic3=Igrid(1,3),Igrid(2,3)
       do jd3=0,Ndens_loc-1
          id3=ic3*Ndens_loc+jd3-Igrid_dense(1,3)+1

          do ic2=Igrid(1,2),Igrid(2,2)
          do jd2=0,Ndens_loc-1
             id2=ic2*Ndens_loc+jd2-Igrid_dense(1,2)+1

             do ic1=Igrid(1,1),Igrid(2,1)
             do jd1=0,Ndens_loc-1
                id1=ic1*Ndens_loc+jd1-Igrid_dense(1,1)+1

                v=(0.0d0,0.0d0)
                do itp3=nitp_0,nitp_1
                do itp2=nitp_0,nitp_1
                do itp1=nitp_0,nitp_1

                   jc1 = mod(ic1+itp1+ML1,ML1)
                   jc2 = mod(ic2+itp2+ML2,ML2)
                   jc3 = mod(ic3+itp3+ML3,ML3)
                   i   = LLL(jc1,jc2,jc3)

                   v = v + Clag1(itp1,jd1)*Clag2(itp2,jd2)*Clag3(itp3,jd3) &
                           *utmp(i)

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

    deallocate( utmp )
    deallocate( LLL )

  END SUBROUTINE localpot2_density

END MODULE localpot2_density_module
