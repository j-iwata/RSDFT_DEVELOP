MODULE pseudopot_module

  use var_ps_member
  use ps_read_PSV
  use var_ps_member_g, only: allocatePSG, sendPSG, npq, ddi, qqr, nl3v, l3v, qrL
  use ps_read_TM_module
  use ps_read_YB_module
  use ps_read_UPF_module
  use ps_gth_module
  use io_tools_module
  use virtualH_module
  use simc_module, only: fit_initrho_simc

  implicit none

  PRIVATE
  PUBLIC :: ippform,file_ps,inorm,NRps,norb,Mr,lo,no,vql,cdd,cdc,rad &
           ,anorm,viod,Rps,Zps,parloc,rab,cdd_coef,ps_type,Rcloc &
           ,hnml,knml,hnl,knl

  PUBLIC :: read_pseudopot

  integer,PUBLIC :: pselect = 2
  logical,PUBLIC :: flag_so = .false.

  integer :: Nelement
  integer :: unit_ps,ielm

CONTAINS


  SUBROUTINE read_ppname_pseudopot
    implicit none
    integer :: i
    call IOTools_readIntegerString( "PP", ippform(1), file_ps(1) )
    do i=2,Nelement
       call IOTools_readIntegerString( "PP", ippform(i), file_ps(i), norewind=.true. )
    end do
  END SUBROUTINE read_ppname_pseudopot


  SUBROUTINE read_param_pseudopot
    implicit none
    call IOTools_readIntegerKeyword( "PSELECT", pselect )
    call IOTools_findKeyword( "SPINORBIT", flag_so, flag_bcast=.true. )
  END SUBROUTINE read_param_pseudopot


  SUBROUTINE read_pseudopot( Nelement_in, rank )

    implicit none
    integer,intent(IN) :: Nelement_in, rank
    real(8),allocatable :: psi_(:,:,:),phi_(:,:,:),bet_(:,:,:)
    real(8),allocatable :: ddi_(:,:,:),qqr_(:,:,:)
    integer :: i,j,io,jo,li,lj
    integer :: Lrefmax,Rrefmax,npqmax,nsmpl
    integer :: max_nterms

    call write_border( 0, " read_pseudopot(start)" )

    Nelement = Nelement_in

    Nelement_PP = Nelement
    Nelement_   = Nelement

    if( allocated(ippform) ) deallocate(ippform)  ! MIZUHO-IR for cellopt
    allocate( ippform(Nelement) ) ; ippform=0
    if( allocated(file_ps) ) deallocate(file_ps)  ! MIZUHO-IR for cellopt
    allocate( file_ps(Nelement) ) ; file_ps=""

    call read_ppname_pseudopot

    call read_param_pseudopot

    if ( any(ippform>100) ) pselect=102

    if ( .not.( pselect==2 .or. pselect==3 .or. pselect==102 ) ) then
       stop "invalid pselect(stop@read_param_pseudopot)"
    end if

    if( allocated(ps) ) deallocate(ps)  ! MIZUHO-IR for cellopt
    allocate( ps(Nelement) )

    if ( rank == 0 ) then

       max_psgrd=0
       max_psorb=0

       do ielm=1,Nelement

          unit_ps=33+ielm
          open(unit_ps,FILE=file_ps(ielm),STATUS='old')

          select case( ippform(ielm) )
          case( 1 )

             close(unit_ps)
             open(unit_ps,FILE=file_ps(ielm),form='unformatted',STATUS='old')

             call ps_read_TM( unit_ps, ps(ielm) )

             call ps_allocate( ps(ielm)%Mr, ps(ielm)%norb )
             Mr(ielm)                 = ps(ielm)%Mr
             norb(ielm)               = ps(ielm)%norb
             Zps(ielm)                = ps(ielm)%Zps
             anorm(1:norb(ielm),ielm) = ps(ielm)%anorm(1:norb(ielm))
             inorm(1:norb(ielm),ielm) = ps(ielm)%inorm(1:norb(ielm))
             Rps(1:norb(ielm),ielm)   = ps(ielm)%Rps(1:norb(ielm))
             NRps(1:norb(ielm),ielm)  = ps(ielm)%NRps(1:norb(ielm))
             lo(1:norb(ielm),ielm)    = ps(ielm)%lo(1:norb(ielm))
             vql(1:Mr(ielm),ielm)     = ps(ielm)%vql(1:Mr(ielm))
             cdd(1:Mr(ielm),ielm)     = ps(ielm)%cdd(1:Mr(ielm))
             cdc(1:Mr(ielm),ielm)     = ps(ielm)%cdc(1:Mr(ielm))
             rad(1:Mr(ielm),ielm)     = ps(ielm)%rad(1:Mr(ielm))
             rab(1:Mr(ielm),ielm)     = ps(ielm)%rab(1:Mr(ielm))
             viod(1:Mr(ielm),1:norb(ielm),ielm) &
                  = ps(ielm)%viod(1:Mr(ielm),1:norb(ielm))

          case( 2, 102 )

             call read_PSV( unit_ps, ielm, ps(ielm) )

             call ps_allocate( ps(ielm)%Mr, ps(ielm)%norb )
             Mr(ielm)                 = ps(ielm)%Mr
             norb(ielm)               = ps(ielm)%norb
             Zps(ielm)                = ps(ielm)%Zps
             anorm(1:norb(ielm),ielm) = ps(ielm)%anorm(1:norb(ielm))
             inorm(1:norb(ielm),ielm) = ps(ielm)%inorm(1:norb(ielm))
             Rps(1:norb(ielm),ielm)   = ps(ielm)%Rps(1:norb(ielm))
             NRps(1:norb(ielm),ielm)  = ps(ielm)%NRps(1:norb(ielm))
             lo(1:norb(ielm),ielm)    = ps(ielm)%lo(1:norb(ielm))
             no(1:norb(ielm),ielm)    = ps(ielm)%no(1:norb(ielm))
             vql(1:Mr(ielm),ielm)     = ps(ielm)%vql(1:Mr(ielm))
             cdd(1:Mr(ielm),ielm)     = ps(ielm)%cdd(1:Mr(ielm))
             cdc(1:Mr(ielm),ielm)     = ps(ielm)%cdc(1:Mr(ielm))
             rad(1:Mr(ielm),ielm)     = ps(ielm)%rad(1:Mr(ielm))
             rab(1:Mr(ielm),ielm)     = ps(ielm)%rab(1:Mr(ielm))
             viod(1:Mr(ielm),1:norb(ielm),ielm) &
                  = ps(ielm)%viod(1:Mr(ielm),1:norb(ielm))
             parloc(1:4,ielm)         = ps(ielm)%parloc(1:4)
             nlf(ielm)                = ps(ielm)%nlf
             nrf(1:norb(ielm),ielm)   = ps(ielm)%nrf(1:norb(ielm))

             if ( pselect == 102 ) then !-----> uspp

                Lrefmax = ps(ielm)%nlf
                Rrefmax = maxval( ps(ielm)%nrf )
                npqmax  = (Lrefmax*Rrefmax*(Lrefmax*Rrefmax+1))/2
                nsmpl   = maxval( ps(ielm)%NRps(:) )
                call allocatePSG( Lrefmax,Rrefmax,npqmax,nsmpl,Nelement )

                npq(ielm) = ps(ielm)%npq

                ddi(1:Rrefmax,1:Rrefmax,1:Lrefmax,ielm) &
                     = ps(ielm)%ddi(1:Rrefmax,1:Rrefmax,1:Lrefmax)
                qqr(1:Rrefmax,1:Rrefmax,1:Lrefmax,ielm) &
                     = ps(ielm)%qqr(1:Rrefmax,1:Rrefmax,1:Lrefmax)
                nl3v(1:npqmax,ielm) = ps(ielm)%nl3v(1:npqmax)
                l3v(1:Lrefmax,1:npqmax,ielm) &
                     = ps(ielm)%l3v(1:Lrefmax,1:npqmax)
                qrL(1:nsmpl,1:Lrefmax,1:npqmax,ielm) = &
                     ps(ielm)%qrL(1:nsmpl,1:Lrefmax,1:npqmax)

             else if (  count( ps(ielm)%anorm /= 0.0d0 ) &
                      < count( ps(ielm)%Dij /= 0.0d0 )  ) then !--> MultiRef

                ps_type=1
                anorm(:,ielm)=1.0d0

             end if

          case( 3 )

             call ps_read_YB( unit_ps, ps(ielm) )

             call ps_allocate( ps(ielm)%Mr, ps(ielm)%norb )
             Mr(ielm)                 = ps(ielm)%Mr
             norb(ielm)               = ps(ielm)%norb
             Zps(ielm)                = ps(ielm)%Zps
             anorm(1:norb(ielm),ielm) = ps(ielm)%anorm(1:norb(ielm))
             inorm(1:norb(ielm),ielm) = ps(ielm)%inorm(1:norb(ielm))
             Rps(1:norb(ielm),ielm)   = ps(ielm)%Rps(1:norb(ielm))
             NRps(1:norb(ielm),ielm)  = ps(ielm)%NRps(1:norb(ielm))
             lo(1:norb(ielm),ielm)    = ps(ielm)%lo(1:norb(ielm))
             vql(1:Mr(ielm),ielm)     = ps(ielm)%vql(1:Mr(ielm))
             cdd(1:Mr(ielm),ielm)     = ps(ielm)%cdd(1:Mr(ielm))
             cdc(1:Mr(ielm),ielm)     = ps(ielm)%cdc(1:Mr(ielm))
             rad(1:Mr(ielm),ielm)     = ps(ielm)%rad(1:Mr(ielm))
             rab(1:Mr(ielm),ielm)     = ps(ielm)%rab(1:Mr(ielm))
             viod(1:Mr(ielm),1:norb(ielm),ielm) &
                  = ps(ielm)%viod(1:Mr(ielm),1:norb(ielm))

          case( 4 )

             call read_ps_gth( unit_ps, ps(ielm) )

             call ps_allocate( 1, ps(ielm)%norb )
             norb(ielm)               = ps(ielm)%norb
             Zps(ielm)                = ps(ielm)%Zps
             if ( norb(ielm) /= 0 ) then
             Rps(1:norb(ielm),ielm)   = ps(ielm)%Rps(1:norb(ielm))
             lo(1:norb(ielm),ielm)    = ps(ielm)%lo(1:norb(ielm))
             no(1:norb(ielm),ielm)    = ps(ielm)%no(1:norb(ielm))
             inorm(1:norb(ielm),ielm) = ps(ielm)%inorm(1:norb(ielm))
             end if
             parloc(1:4,ielm)         = ps(ielm)%parloc(1:4)
             Rcloc(ielm)              = ps(ielm)%Rcloc
             hnl(:,:,ielm)            = ps(ielm)%hnl(:,:)
             knl(:,:,ielm)            = ps(ielm)%knl(:,:)
             hnml(:,:,:,ielm)         = ps(ielm)%hnml(:,:,:)
             knml(:,:,:,ielm)         = ps(ielm)%knml(:,:,:)

             if ( any( hnml /= 0.0d0 ) ) ps_type=1

          case( 5 )

             call ps_read_UPF( unit_ps, ps(ielm) )

             call ps_allocate( ps(ielm)%Mr, ps(ielm)%norb )
             Mr(ielm)                 = ps(ielm)%Mr
             norb(ielm)               = ps(ielm)%norb
             Zps(ielm)                = ps(ielm)%Zps
             if ( ps(ielm)%norb > 0 ) then
             anorm(1:norb(ielm),ielm) = ps(ielm)%anorm(1:norb(ielm))
             inorm(1:norb(ielm),ielm) = ps(ielm)%inorm(1:norb(ielm))
             Rps(1:norb(ielm),ielm)   = ps(ielm)%Rps(1:norb(ielm))
             NRps(1:norb(ielm),ielm)  = ps(ielm)%NRps(1:norb(ielm))
             lo(1:norb(ielm),ielm)    = ps(ielm)%lo(1:norb(ielm))
             no(1:norb(ielm),ielm)    = ps(ielm)%no(1:norb(ielm))
             viod(1:Mr(ielm),1:norb(ielm),ielm) &
                                      = ps(ielm)%viod(1:Mr(ielm),1:norb(ielm))
             end if
             vql(1:Mr(ielm),ielm)     = ps(ielm)%vql(1:Mr(ielm))
             cdd(1:Mr(ielm),ielm)     = ps(ielm)%cdd(1:Mr(ielm))
             cdc(1:Mr(ielm),ielm)     = ps(ielm)%cdc(1:Mr(ielm))
             rad(1:Mr(ielm),ielm)     = ps(ielm)%rad(1:Mr(ielm))
             rab(1:Mr(ielm),ielm)     = ps(ielm)%rab(1:Mr(ielm))

             if ( allocated(ps(ielm)%Dij) ) then
                if ( any( ps(ielm)%Dij /= 0.0d0 ) ) then ! Multireference
                   ps_type=1
                   anorm(:,ielm)=1.0d0
                end if
             end if

          case default

             stop "ippform error"

          end select ! ippform

          close(unit_ps)

       end do ! ielm

       write(*,*) "ps_type = ",ps_type
       if ( ps_type == 1 ) then
          write(*,*) "(non-diagonal partrs are in nonlocal pseudopotential)"
          if ( any(ippform==4) ) then
          else
             do ielm=1,Nelement
                do j=1,ps(ielm)%norb
                   jo=ps(ielm)%no(j)
                   lj=ps(ielm)%lo(j)
                do i=1,ps(ielm)%norb
                   io=ps(ielm)%no(i)
                   li=ps(ielm)%lo(i)
                   if ( li /= lj ) cycle
                   hnml(io,jo,li,ielm) = ps(ielm)%Dij(i,j)
                end do
                end do
                do j=1,norb(ielm)
                   viod(:,j,ielm)=viod(:,j,ielm)/sqrt( anorm(j,ielm) )
                end do
             end do ! ielm
          end if
       end if

       write(*,*) "# of Gaussian parameters of initial density"
       max_ngauss=0
       max_nterms=0
       do ielm=1,Nelement
          if ( allocated(ps(ielm)%cdd_coef) ) then
             max_ngauss = max( max_ngauss, ps(ielm)%ngauss )
             max_nterms = max( max_nterms, size(ps(ielm)%cdd_coef,1) )
             write(*,*) ielm, ps(ielm)%ngauss, size(ps(ielm)%cdd_coef,1)
          end if
       end do
       write(*,*) "max_ngauss,max_nterms=",max_ngauss,max_nterms
       if ( max_ngauss > 0 ) then
          allocate( cdd_coef(max_nterms,max_ngauss,Nelement) ) ; cdd_coef=0.0d0
          do ielm=1,Nelement
             do j=1,ps(ielm)%ngauss
                do i=1,size(ps(ielm)%cdd_coef,1)
                   cdd_coef(i,j,ielm)=ps(ielm)%cdd_coef(i,j)
                end do
             end do
          end do
       else
          write(*,*) "Gaussian fit is applied"
          allocate( cdd_coef(3,4,Nelement) ) ; cdd_coef=0.0d0
          do ielm=1,Nelement
             call fit_initrho_simc(rad(:,ielm),cdd(:,ielm),cdd_coef(:,:,ielm))
          end do
       end if

    end if ! [ rank == 0 ]

! --- bcast pseudopotential data

    call send_pseudopot(rank)
    if ( pselect > 100 ) call sendPSG( rank, Nelement )

    do ielm=1,Nelement
       call ps_send_ps1d( ps(ielm) )
       if ( pselect > 100 ) call psg_send_ps1d( ps(ielm) )
    end do

! ---

!    call chk_pot(1,rank)

    call virtualH( Zps, vql )

    call write_border( 0, " read_pseudopot(end)" )

  END SUBROUTINE read_pseudopot


  SUBROUTINE chk_pot(iflag,rank)
    implicit none
    integer,intent(IN) :: iflag,rank
    integer :: u,ielm,i,j
    if ( rank == 0 ) then
       do ielm=1,Nelement
          u=9+ielm
          rewind u
          do i=1,Mr(ielm)
             if ( iflag == 2 ) then
                write(u,'(1x,5f20.10)') rad(i,ielm),cdd(i,ielm),cdc(i,ielm)
             else
                write(u,'(1x,5f20.10)') &
                     rad(i,ielm),vql(i,ielm),(viod(i,j,ielm),j=1,norb(ielm))
             end if
          end do
          write(*,'(1x,"chk_pot(",i1,"): fort.",i2)') iflag,u
       end do
    end if
    stop "stop@chk_pot"
  END SUBROUTINE chk_pot


END MODULE pseudopot_module
