MODULE wf_module

  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: unk,esp,occ,res,init_wf,test_on_wf,gather_wf &
           ,ML_WF, ML_0_WF, ML_1_WF, MB_WF, MB_0_WF, MB_1_WF &
           ,MK_WF, MK_0_WF, MK_1_WF, MS_WF, MS_0_WF, MS_1_WF &
           ,Sunk

#ifdef _DRSDFT_
  real(8),parameter :: zero=0.d0
  real(8),allocatable :: unk(:,:,:,:)
  integer,parameter :: TYPE_MAIN=MPI_REAL8
#else
  complex(8),parameter :: zero=(0.d0,0.d0)
  complex(8),allocatable :: unk(:,:,:,:)
  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
#endif

#ifdef _USPP_
#ifdef _DRSDFT_
    real(8),allocatable :: Sunk(:,:)
#else
    complex(8),allocatable :: Sunk(:,:)
#endif
#endif

  real(8),allocatable :: esp(:,:,:)
  real(8),allocatable :: occ(:,:,:)
  real(8),allocatable :: res(:,:,:)

  integer :: ML_WF, ML_0_WF, ML_1_WF
  integer :: MB_WF, MB_0_WF, MB_1_WF
  integer :: MK_WF, MK_0_WF, MK_1_WF
  integer :: MS_WF, MS_0_WF, MS_1_WF

CONTAINS


  SUBROUTINE init_wf
    implicit none

    if ( myrank == 0 ) write(*,'(a60," init_wf")') repeat("-",60)

    ML_WF   = sum( ir_grid )
    ML_0_WF = id_grid(myrank_g) + 1
    ML_1_WF = id_grid(myrank_g) + ir_grid(myrank_g)

    MB_WF   = sum( ir_band )
    MB_0_WF = id_band(myrank_b) + 1
    MB_1_WF = id_band(myrank_b) + ir_band(myrank_b)

    MK_WF   = sum( ir_bzsm )
    MK_0_WF = id_bzsm(myrank_k) + 1
    MK_1_WF = id_bzsm(myrank_k) + ir_bzsm(myrank_k)

    MS_WF   = sum( ir_spin )
    MS_0_WF = id_spin(myrank_s) + 1
    MS_1_WF = id_spin(myrank_s) + ir_spin(myrank_s)

    if ( allocated(occ) ) deallocate(occ)
    if ( allocated(res) ) deallocate(res)
    if ( allocated(esp) ) deallocate(esp)
    if ( allocated(unk) ) deallocate(unk)

    allocate( unk(ML_0_WF:ML_1_WF,MB_WF,MK_0_WF:MK_1_WF,MS_0_WF:MS_1_WF) )
    unk=zero
    allocate( esp(MB_WF,MK_WF,MS_WF) )
    esp=0.0d0
    allocate( res(MB_WF,MK_WF,MS_WF) )
    res=0.0d0
    allocate( occ(MB_WF,MK_WF,MS_WF) )
    occ=0.0d0

    call random_initial_wf

  END SUBROUTINE init_wf


  SUBROUTINE random_initial_wf
    implicit none
    integer :: s,k,n,i
    integer,allocatable :: ir(:)
    real(8) :: u(2)

    call random_seed( size=n )
    allocate( ir(n) )
    ir(:)=MB_0_WF+ML_0_WF
    call random_seed( put=ir )
    deallocate( ir )
 
    do s=MS_0_WF,MS_1_WF
       do k=MK_0_WF,MK_1_WF
          do n=MB_0_WF,MB_1_WF
             do i=ML_0_WF,ML_1_WF
                call random_number(u)
                unk(i,n,k,s)=dcmplx(u(1),u(2))
             end do
          end do
       end do
    end do
do s=MS_0_WF,MS_1_WF
do k=MK_0_WF,MK_1_WF
do n=MB_0_WF,MB_1_WF
do i=ML_0_WF,ML_1_WF
write(580+myrank,*) s,k,n,i,unk(i,n,k,s)
end do
end do
end do
end do

  END SUBROUTINE random_initial_wf


  SUBROUTINE test_on_wf(dV,disp_switch)
    implicit none
    real(8),intent(IN) :: dV          ! volume element
    logical,intent(IN) :: disp_switch ! diplay switch
    integer :: ierr,s,k,m,n,mm
#ifdef _DRSDFT_
    real(8),allocatable :: uu(:,:)
#else
    complex(8),allocatable :: uu(:,:)
#endif

    allocate( uu(MB_WF,MB_WF) ) ; uu=zero

    mm = size(uu)

    do s=MS_0_WF,MS_1_WF
    do k=MK_0_WF,MK_1_WF
       uu(:,:)=zero
       do n=1,MB_WF
       do m=1,n
#ifdef _DRSDFT_
          uu(m,n)=sum( unk(:,m,k,s)*unk(:,n,k,s) )*dV
#else
          uu(m,n)=sum( conjg(unk(:,m,k,s))*unk(:,n,k,s) )*dV
#endif
       end do ! m
       end do ! n
       call mpi_allreduce(MPI_IN_PLACE,uu,mm,TYPE_MAIN &
            ,MPI_SUM,comm_grid,ierr)
       do n=1,MB_WF
       do m=1,n
          if ( disp_switch ) then
             write(590,'(1x,i2,i5,2i7,2g25.16)') s,k,m,n,uu(m,n)
          end if
       end do ! m
       end do ! n
    end do ! k
    end do ! s

    deallocate( uu )

  END SUBROUTINE test_on_wf


  SUBROUTINE gather_wf
    implicit none
    integer :: k,s,mm,ierr
    mm=ML_1_WF-ML_0_WF+1
    ir_band(:)=ir_band(:)*mm
    id_band(:)=id_band(:)*mm
    do s=MS_0_WF,MS_1_WF
    do k=MK_0_WF,MK_1_WF
       call mpi_allgatherv( unk(ML_0_WF,MB_0_WF,k,s),ir_band(myrank_b) &
            ,TYPE_MAIN,unk(ML_0_WF,1,k,s),ir_band,id_band &
            ,TYPE_MAIN,comm_band,ierr )
    end do
    end do
    ir_band(:)=ir_band(:)/mm
    id_band(:)=id_band(:)/mm
  END SUBROUTINE gather_wf


END MODULE wf_module
