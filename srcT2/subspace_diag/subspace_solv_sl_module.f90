MODULE subspace_solv_sl_module

  use wf_module, only: esp
  use scalapack_module
  use subspace_diag_variables
#ifdef _EIGEN_
  use eigen_libs_mod
#endif
  use parallel_module, only: myrank

  implicit none

  PRIVATE
  PUBLIC :: subspace_solv_sl

CONTAINS


  SUBROUTINE subspace_solv_sl(k,s)
    implicit none
    include 'mpif.h'
    integer,intent(IN) :: k,s
    integer :: itmp(1),LWORK0,LRWORK0,LIWORK0,TRILWMIN,ierr,MB
    integer,save :: LWORK=0,LRWORK=0,LIWORK=0
    integer,allocatable :: iwork(:)
    real(8) :: rtmp(1)
    real(8),allocatable :: rwork(:)
    complex(8) :: ctmp(1)
    complex(8),allocatable :: zwork(:)
!   character(8) :: idiag
    integer :: myrank, i , j, nn, ista
    logical :: dflag
    integer :: i0, j0 
    logical :: idisp=.false.

    call write_border( 1, "subspace_solv_sl(start)" )

    MB = MB_diag
    ierr = 0

    ! get myrank
    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
    if (myrank==0) idisp=.true.

    select case(idiag)
    case('pzheevd')

       if ( LWORK==0 ) then
          call pzheevd('V',UPLO,MB,Hsub,1,1,DESCA,esp(1,k,s),Vsub,1,1 &
                       ,DESCZ,ctmp,-1,rtmp,-1,itmp,-1,ierr)
          LWORK =nint(real(ctmp(1)))
          LRWORK=nint(rtmp(1))
          LIWORK=itmp(1)
       end if
       LWORK =max(LWORK,MB+(NP0+NQ0+MBSIZE)*MBSIZE)
       LRWORK=max(LRWORK,(1+8*MB+2*NPX*NQX)*2)
       LIWORK=max(LIWORK,7*MB+8*NPCOL+2)

       allocate( zwork(LWORK),rwork(LRWORK),iwork(LIWORK) )

       call pzheevd('V',UPLO,MB,Hsub,1,1,DESCA,esp(1,k,s),Vsub,1,1 &
                    ,DESCZ,zwork,LWORK,rwork,LRWORK,iwork,LIWORK,ierr)

       deallocate( iwork,rwork,zwork )

    case('pzheev')

       if ( LWORK==0 ) then
          call pzheev('V',UPLO,MB,Hsub,1,1,DESCA,esp(1,k,s),Vsub,1,1 &
                      ,DESCZ,ctmp,-1,rtmp,-1,ierr)
          LWORK =nint(real(ctmp(1)))
          LRWORK=nint(rtmp(1))
       end if

       LWORK =max(LWORK,(NP0+NQ0+MBSIZE)*MBSIZE+3*MB+MB**2)
       LRWORK=max(LRWORK,4*MB-2)

       allocate( zwork(LWORK),rwork(LRWORK) )

       call pzheev('V',UPLO,MB,Hsub,1,1,DESCA,esp(1,k,s),Vsub,1,1 &
            ,DESCZ,zwork,LWORK,rwork,LRWORK,ierr)

       deallocate( rwork,zwork )

    case('pdsyevd')
       !debug
       write(*,*) ">>>  exec pdsyevd  ",myrank

       TRILWMIN = 3*MB + max( MBSIZE*(NPX+1),3*MBSIZE )
       LRWORK0  = max( 1+6*MB+2*NPX*NQX, TRILWMIN )
       LIWORK0  = 7*MB+8*NPCOL+2

       if ( LRWORK==0 ) then
          call pdsyevd('V',UPLO,MB,Hsub,1,1,DESCA,esp(1,k,s),Vsub,1,1 &
                       ,DESCZ,rtmp,-1,itmp,LIWORK,ierr)
          LRWORK=nint(rtmp(1))
          LRWORK=LRWORK*10
       end if
       LRWORK=max(LRWORK,LRWORK0*10)
       LIWORK=max(LIWORK,LIWORK0)

       allocate( rwork(LRWORK),iwork(LIWORK) )

       call pdsyevd('V',UPLO,MB,Hsub,1,1,DESCA,esp(1,k,s),Vsub,1,1 &
                    ,DESCZ,rwork,LRWORK,iwork,LIWORK,ierr)

       deallocate( iwork,rwork )

    case('pdsyev')

       LRWORK0=0

       if ( LRWORK==0 ) then
          call pdsyev('V',UPLO,MB,Hsub,1,1,DESCA,esp(1,k,s),Vsub,1,1,DESCZ,rtmp,-1,ierr)
          LRWORK=nint(rtmp(1))
       end if
       LRWORK=max(LRWORK,LRWORK0)

       allocate( rwork(LRWORK) )

       call pdsyev('V',UPLO,MB,Hsub,1,1,DESCA,esp(1,k,s),Vsub,1,1,DESCZ,rwork,LRWORK,ierr)

       deallocate( rwork )

    case('check_u')
       !debug (check created symetric matrix)
       TRILWMIN = 3*MB + max( MBSIZE*(NPX+1),3*MBSIZE )
       LRWORK0  = max( 1+6*MB+2*NPX*NQX, TRILWMIN )
       LIWORK0  = 7*MB+8*NPCOL+2

       if ( LRWORK==0 ) then
          call pdsyevd('V','U',MB,Hsub,1,1,DESCA,esp(1,k,s),Vsub,1,1 &
                       ,DESCZ,rtmp,-1,itmp,LIWORK,ierr)
          LRWORK=nint(rtmp(1))
          LRWORK=LRWORK*10
       end if
       LRWORK=max(LRWORK,LRWORK0*10)
       LIWORK=max(LIWORK,LIWORK0)

       allocate( rwork(LRWORK),iwork(LIWORK) )

       call pdsyevd('V','U',MB,Hsub,1,1,DESCA,esp(1,k,s),Vsub,1,1 &
                    ,DESCZ,rwork,LRWORK,iwork,LIWORK,ierr)

       if (myrank==0) then
         write(*,*) "************ Vsub *************"
         do i=1,10
            write(*,'(10(1x,e13.6))') Vsub(i,1:10)
         enddo
         write(10,*) Vsub(1,:)
       endif

       deallocate( iwork,rwork )

    case('eigen_s')
#ifdef _EIGEN_
!      if (idisp) call show_matrix_val(Hsub,lld_r,1,10,1,10,6, "Hsub  ")       !debug

       ! create symetric matrix Hsub
       if (imate==1) call fill_Matrix_scalapack(Hsub ,MB)

!      if (idisp) call show_matrix_val(Hsub,lld_r,1,10,1,10,6, "Hsub2 ")       !debug

       ! re-distribute form block-cyclic to cyclic
       call trans_blkcy2cy(Hsub, Hsub_e, MB)

!      if (idisp) call show_matrix_val(Hsub_e,LLD_R_e,1,10,1,10,6, "Hsub_e")   !debug

       ! solver (eigen_s)
       if (iflag_e) then
!         write(*,*) ">>>  exec eigen_s  ",myrank                              !debug
          call eigen_s(MB,MB,Hsub_e,LLD_R_e,esp(1,k,s),Vsub_e,LLD_R_e)
!         call eigen_sx(MB,MB,Hsub_e,LLD_R_e,esp(1,k,s),Vsub_e,LLD_R_e)
       endif

!      if (idisp) call show_matrix_val(Vsub_e,LLD_R_e,1,10,1,10,6, "Vsub_e")   !debug

       call mpi_barrier(mpi_comm_world,ierr)
       ! form cyclic to block-cyclic 
       call trans_cy2blkcy(Vsub_e, Vsub, MB)

!      if (idisp) call show_matrix_val(Vsub,lld_r,1,10,1,10,6, "Vsub  ")       !debug
#endif
    case default

       write(*,*) "idiag=",idiag
       write(*,*) "This routine is not available!"
       stop

    end select

    !debug eigen-value
    !if (idisp) then
    !   do i=1,MB
    !     write(12,'(i4,3x,d23.16)') i, esp(i,1,1)
    !   enddo
    !endif

    if ( ierr /= 0 ) then
       write(*,*) "ierr,idiag=",ierr,idiag
       stop
    end if

    call write_border( 1, "subspace_solv_sl(end)" )

    return

  END SUBROUTINE subspace_solv_sl

  subroutine show_matrix_val(a, lda, n1,n2, n3,n4, nunit, ntag)
      implicit none
      real(8),intent(in) :: a(lda,*)
      integer,intent(in) :: n1,n2,n3,n4, nunit, lda
      character(6), intent(in) :: ntag
      character(17) :: nfmt      
      integer :: n, i

      n=n4-n3+1
      write(nfmt,'("("i0"e14.6)")') n
 
      write(nunit,*) " ********* ", ntag, " ********* "
      do i=n1,n2
         write(nunit,nfmt) a(i,n3:n4)
      enddo
      return
  end subroutine show_matrix_val


END MODULE subspace_solv_sl_module
