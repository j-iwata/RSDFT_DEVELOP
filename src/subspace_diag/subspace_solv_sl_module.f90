module subspace_solv_sl_module

  use wf_module
  use scalapack_module
  use subspace_diag_variables, only: MB_diag,Hsub,Vsub
  use watch_module

  implicit none

  private
  public :: subspace_solv_sl

contains


  subroutine subspace_solv_sl(k,s)
    implicit none
    integer,intent(IN) :: k,s
    integer :: itmp(1),LWORK0,LRWORK0,LIWORK0,TRILWMIN,ierr,MB
    integer,save :: LWORK=0,LRWORK=0,LIWORK=0
    integer,allocatable :: iwork(:)
    real(8) :: rtmp(1)
    real(8),allocatable :: rwork(:)
    complex(8) :: ctmp(1)
    complex(8),allocatable :: zwork(:)
    character(8) :: idiag
    type(time) :: t
    logical :: disp_sw

#ifndef _LAPACK_

    call write_border( 1, " subspace_solv_sl(start)" )
    call start_timer( t_out=t )

    MB = MB_diag

    ierr = 0

    call check_disp_switch( disp_sw, 0 )

#ifdef _DRSDFT_
    idiag = 'PDSYEVD'
#else
    idiag = 'PZHEEVD'
#endif

    select case( idiag )
    case( 'PZHEEVD' )

      if ( LWORK==0 ) then
        call PZHEEVD('V',UPLO,MB,Hsub,1,1,DESCA,esp(1,k,s),Vsub,1,1 &
                    ,DESCZ,ctmp,-1,rtmp,-1,itmp,-1,ierr)
        LWORK =nint(real(ctmp(1)))
        LRWORK=nint(rtmp(1))
        LIWORK=itmp(1)
        ! if ( disp_sw ) then
        !   write(*,'("(PZHEEVD) Work-array sizes (by query)   : LWORK,LRWORK,LIWORK=",3i8)') LWORK,LRWORK,LIWORK
        !   write(*,'("(PZHEEVD) Work-array sizes (in document): LWORK,LRWORK,LIWORK=",3i8)') &
        !     MB+(NP0+NQ0+MBSIZE)*MBSIZE, 1+9*MB+3*NPX*NQX, 7*MB+8*NPCOL+2
        ! end if
        LWORK =max( LWORK , MB+(NP0+NQ0+MBSIZE)*MBSIZE )
        LRWORK=max( LRWORK, (1+9*MB+3*NPX*NQX) )
        LIWORK=max( LIWORK, 7*MB+8*NPCOL+2 )
        ! if ( disp_sw ) then
        !   write(*,'("(PZHEEVD) Work-array sizes (allocate)   : LWORK,LRWORK,LIWORK=",3i8)') LWORK,LRWORK,LIWORK
        ! end if
      end if

      allocate( zwork(LWORK),rwork(LRWORK),iwork(LIWORK) )

      call PZHEEVD('V',UPLO,MB,Hsub,1,1,DESCA,esp(1,k,s),Vsub,1,1 &
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

    case( 'PDSYEVD' )

      TRILWMIN = 3*MB + max( MBSIZE*(NPX+1),3*MBSIZE )
      LRWORK0  = max( 1+6*MB+2*NPX*NQX, TRILWMIN )
      LIWORK0  = 7*MB+8*NPCOL+2

      if ( LRWORK==0 ) then
        call PDSYEVD('V',UPLO,MB,Hsub,1,1,DESCA,esp(1,k,s),Vsub,1,1 &
                     ,DESCZ,rtmp,-1,itmp,LIWORK,ierr)
        LRWORK=nint(rtmp(1))
        LRWORK=LRWORK*10
      end if
      LRWORK=max(LRWORK,LRWORK0*10)
      LIWORK=max(LIWORK,LIWORK0)

      allocate( rwork(LRWORK),iwork(LIWORK) )

      call PDSYEVD('V',UPLO,MB,Hsub,1,1,DESCA,esp(1,k,s),Vsub,1,1 &
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

    case default

       write(*,*) "idiag=",idiag
       write(*,*) "This routine is not available!"
       stop

    end select

    if ( ierr /= 0 ) then
       write(*,*) "ierr,idiag=",ierr,idiag
       stop
    end if

    call result_timer( "solv_sl", t )
    call write_border( 1, " subspace_solv_sl(end)" )

    return

#endif

  end subroutine subspace_solv_sl

end module subspace_solv_sl_module