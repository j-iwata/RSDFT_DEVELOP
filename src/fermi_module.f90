MODULE fermi_module

  use bberf_module
  use io_tools_module
  use watch_module, only: start_timer,result_timer

  implicit none

  PRIVATE
  PUBLIC :: calc_fermi
  PUBLIC :: ekbt, efermi, Eentropy

  real(8) :: ekbt, efermi, Eentropy
  integer :: mb1,mb2,kinteg,mp_version

  integer :: nsetocc,isetocc(2)
  real(8) :: setocc(20)

  logical :: first_time = .true.
  real(8),allocatable :: factor(:)
  logical :: flag_read = .true.
  logical,allocatable :: msk(:,:,:)

CONTAINS


  SUBROUTINE read_fermi
    implicit none
    real(8) :: rtmp(22)
    ekbt=1.d-5
    kinteg=5
    nsetocc=0
    isetocc=0
    setocc=0.0d0
    rtmp=0.0d0
    mp_version=0
    call IOTools_readReal8Keyword( "EKBT", ekbt )
    call IOTools_readIntegerKeyword( "KINTEG", kinteg )
    call IOTools_readIntegerKeyword( "MPVER", mp_version )
    call IOTools_readReal8Keywords( "SETOCC", rtmp )
    isetocc(1) = nint( rtmp(1) )
    isetocc(2) = nint( rtmp(2) )
    if ( isetocc(1)>0 .and. isetocc(2)>0 .and. isetocc(1)<=isetocc(2) ) then
       nsetocc = isetocc(2) - isetocc(1) + 1
       setocc(1:nsetocc) = rtmp(3:nsetocc+2)
    end if
    flag_read=.false.
  END SUBROUTINE read_fermi


  SUBROUTINE calc_fermi(iter,nfixed,znel,dspn,esp,wbz,occ)

    implicit none
    integer,intent(IN)  :: iter,nfixed
    real(8),intent(IN)  :: esp(:,:,:),wbz(:),znel,dspn
    real(8),intent(OUT) :: occ(:,:,:)
    real(8) :: oc0
    integer :: n,k,s

    call start_timer

    if ( flag_read ) call read_fermi

    if ( nsetocc > 0 ) then

       occ(:,:,:)=0.0d0
       oc0=2.0d0/size(esp,3)
       do s=1,size(esp,3)
       do k=1,size(esp,2)
          do n=1,isetocc(1)-1
             occ(n,k,s) = wbz(k)*oc0
          end do
          do n=isetocc(1),isetocc(2)
             occ(n,k,s) = setocc(n-isetocc(1)+1)*wbz(k)*oc0
          end do
       end do
       end do

       efermi=maxval( esp(isetocc(2),:,:) )
       Eentropy=0.0d0

    else

       select case( mp_version )
       case default
          call calc_fermi_mp( iter,nfixed,znel,dspn,esp,wbz,occ )
       case( 1 )
          call calc_fermi_mp_1( iter,nfixed,znel,dspn,esp,wbz,occ )
       end select
       call calc_entropy_mp( esp, efermi, ekbt, wbz, Eentropy )

    end if

    call result_timer( "fermi" )

  END SUBROUTINE calc_fermi


  SUBROUTINE calc_fermi_mp(iter,nfixed,znel,dspn,esp,wbz,occ)
    implicit none
    integer,intent(IN)  :: iter,nfixed
    real(8),intent(IN)  :: esp(:,:,:),wbz(:),znel,dspn
    real(8),intent(OUT) :: occ(:,:,:)
    real(8) :: ef1,ef2,ef,ef0,oc0,zne,Ntmp,ff,xx,d,Nele1,Nele2,sqrt_pi
    integer :: n,k,s,ichk_conv,loop,MBN,MBZ,MSP
    real(8),parameter :: eps=1.0d-12, e_mergin=0.0004d0
    integer,parameter :: max_loop=10000
    logical :: disp_switch

    call write_border( 1, " calc_fermi_mp(start)" )

    call check_disp_switch( disp_switch, 0 )

    MBN = size( esp, 1 )
    MBZ = size( esp, 2 )
    MSP = size( esp, 3 )

    if ( first_time ) call prep_ff0

    allocate( msk(size(esp,1),size(esp,2),size(esp,3)) ) ; msk=.true.

! --- bi-section search ---

    if ( MSP == 1 .or. iter > Nfixed ) then

       call get_ef( znel, esp, MSP, ef )
       do loop=1,max_loop
          ef1 = ef - loop*e_mergin
          ef2 = ef + loop*e_mergin
          call calc_num_electron( ef1, esp, MSP, wbz, Nele1 )
          call calc_num_electron( ef2, esp, MSP, wbz, Nele2 )
          d = ( znel - Nele1 )*( znel - Nele2 )
          if ( d <= 0.0d0 ) exit
       end do ! loop
       if ( loop > max_loop ) then
          if ( disp_switch ) then
             write(*,'(1x,5f20.10)') ef,ef1,ef2,Nele1,Nele2
             call stop_program( "stop@calc_fermi(1)" )
          end if
       end if

       ef0 = 1.0d10
       oc0 = 2.0d0/dble(MSP)

       do loop=1,max_loop

          ef = 0.5d0*( ef1 + ef2 )
          if ( ef == ef0 ) goto 100

          Ntmp=0.0d0
          do s=1,size(esp,3)
          do k=1,size(esp,2)
          do n=1,size(esp,1)
             xx = ( esp(n,k,s) - ef )/ekbt
             ff = ff0(kinteg,xx)
             occ(n,k,s) = ff*wbz(k)*oc0
             Ntmp = Ntmp + occ(n,k,s)
          end do
          end do
          end do

          d = Ntmp - znel
          if ( d > eps ) then
             ef2=ef
          else if ( d < -eps ) then
             ef1=ef
          else
             goto 100
          end if

          ef0 = ef

       end do ! loop

    else

       call write_string( "total spin density is fixed!!" )

       ichk_conv=0

       do s=1,MSP

          zne = 0.5d0*znel + (3-2*s)*0.5d0*dspn

          call get_ef( zne, esp(:,:,s:s), MSP, ef )

          do loop=1,max_loop
             ef1 = ef - loop*e_mergin
             ef2 = ef + loop*e_mergin
             call calc_num_electron( ef1, esp(:,:,s:s), MSP, wbz, Nele1 )
             call calc_num_electron( ef2, esp(:,:,s:s), MSP, wbz, Nele2 )
             d = ( zne - Nele1 )*( zne - Nele2 )
             if ( d <= 0.0d0 ) exit
          end do ! loop
          if ( loop > max_loop ) call stop_program( "stop@calc_fermi(2)" )

          do loop=1,max_loop

             ef = 0.5d0*(ef1+ef2)

             Ntmp=0.0d0
             do k=1,size(esp,2)
             do n=1,size(esp,1)
                xx = ( esp(n,k,s) - ef )/ekbt
                ff = ff0(kinteg,xx)
                occ(n,k,s) = ff*wbz(k)
                Ntmp = Ntmp + occ(n,k,s)
             end do
             end do

             d = Ntmp - zne
             if ( d > eps ) then
                ef2=ef
             else if ( d < -eps ) then
                ef1=ef
             else
                ichk_conv = ichk_conv + 1
                exit
             end if

          end do ! loop

       end do ! s

       if ( ichk_conv == 2 ) goto 100

    end if

    if ( abs(Ntmp-zne) > 1.d-10 ) then
       if ( disp_switch ) then
          write(6,*)' EF IS NOT CONVERGED'
          write(6,*)' Check the # of electron, mb1, and mb2'
          write(6,*)' EF1 & EF2=',ef1,ef2
          write(6,*)' Ntmp,zne=',Ntmp,zne,Ntmp-zne
          do s=1,MSP
          do k=1,MBZ
             write(*,*) "s,k =",s,k
             do n=1,MBN
                write(*,*) n,occ(n,k,s),esp(n,k,s)
             end do
          end do
          end do
       end if
       call stop_program( "stop@calc_fermi(3)" )
    end if

100 continue

    efermi = ef

    deallocate( msk )

    call write_border( 1, " calc_fermi_mp(end)" )

  END SUBROUTINE calc_fermi_mp


  SUBROUTINE prep_ff0
    implicit none
    real(8) :: sqrt_pi
    integer :: n
    allocate( factor(0:2*kinteg) ) ; factor=0.0d0
    sqrt_pi=sqrt(acos(-1.0d0))
    factor(0)= 1.0d0/sqrt_pi
    factor(1)=-1.0d0/(4.0d0*sqrt_pi)
    do n=2,2*kinteg
       factor(n)=-factor(n-1)/(4.0d0*n)
    end do
    first_time=.false.
  END SUBROUTINE prep_ff0

  FUNCTION ff0(n,x)
    implicit none
    integer :: n,i
    real(8) :: x,ff,ff0,hp0,hp1,hp2,hp3
    ff0 = 0.5d0*(1.0d0-bberf(x))
    if ( n <= 0 ) return
    hp0 = 1.0d0
    hp1 = 2.0d0*x
    ff  = factor(1)*hp1
    do i=2,n
       hp2 = 2.0d0*x*hp1 - 2.0d0*(2*i-3)*hp0
       hp3 = 2.0d0*x*hp2 - 2.0d0*(2*i-2)*hp1
       ff  = ff + factor(i)*hp3
       hp0 = hp2
       hp1 = hp3
    end do
    ff0 = ff0 + ff*exp(-x*x)
    return
  END FUNCTION ff0


  SUBROUTINE get_ef( Nelectron, esp, msp, ef )
    implicit none
    real(8),intent(IN) :: Nelectron
    real(8),intent(IN) :: esp(:,:,:)
    integer,intent(IN) :: msp
    real(8),intent(OUT) :: ef
    integer :: n,k,s,i
    real(8) :: Ncount,e1,e2
    msk=.true.
    Ncount=0.0d0
    e1=minval(esp)
    loop: do s=1,size(esp,3)
    do k=1,size(esp,2)
    do n=1,size(esp,1)
       e2=minval(esp,mask=msk)
       do i=1,2/msp
          Ncount=Ncount+1.0d0
          if ( abs(Ncount-Nelectron)<1.d-10 ) then
             ef=e2
             exit loop
          else if ( Ncount > Nelectron ) then
             ef=(e1+e2)*0.5d0
             exit loop
          end if
       end do
       msk(n,k,s)=.false.
       e1=e2
    end do
    end do
    end do loop
  END SUBROUTINE get_ef


  SUBROUTINE calc_num_electron( ef, esp, msp, wbz, Nele )
    implicit none
    real(8),intent(IN)  :: ef
    real(8),intent(IN)  :: esp(:,:,:), wbz(:)
    integer,intent(IN)  :: msp
    real(8),intent(OUT) :: Nele
    integer :: n,k,s
    real(8) :: occ0,xx,ff
    occ0=2.0d0/msp
    Nele=0.0d0
    do s=1,size(esp,3)
    do k=1,size(esp,2)
    do n=1,size(esp,1)
       xx = ( esp(n,k,s) - ef )/ekbt
       ff = ff0(kinteg,xx)
       Nele = Nele + ff*wbz(k)*occ0
    end do
    end do
    end do
  END SUBROUTINE calc_num_electron


  SUBROUTINE calc_entropy_mp( esp, ef, kT, wbz, entropy )
    implicit none
    real(8),intent(IN)  :: esp(:,:,:), ef, kT, wbz(:)
    real(8),intent(OUT) :: entropy
    integer :: n,k,s
    real(8) :: x
    entropy=0.0d0
    do s=1,size(esp,3)
    do k=1,size(esp,2)
    do n=1,size(esp,1)
       x=(esp(n,k,s)-ef)/kT
       entropy = entropy + ff0(2*kinteg,x)*wbz(k)
    end do
    end do
    end do
    entropy = entropy*0.5d0*factor(kinteg)*kT
  END SUBROUTINE calc_entropy_mp


  SUBROUTINE calc_fermi_mp_1(iter,nfixed,znel,dspn,esp,wbz,occ)

    implicit none
    integer,intent(IN)  :: iter,nfixed
    real(8),intent(IN)  :: esp(:,:,:),wbz(:),znel,dspn
    real(8),intent(OUT) :: occ(:,:,:)
    real(8) :: ef1,ef2,ef,ef0
    integer :: id,n,k,s,efconv,nn,MB,MBZ,MSP
    real(8) :: zne,octmp,ff,xx
    real(8),parameter :: eps=0.d0
    integer,parameter :: mxcycl=1000
    logical :: disp_switch

    call write_border( 1, " calc_fermi_mp_1(start)" )
    call check_disp_switch( disp_switch, 0 )

    MB  = size( esp, 1 )
    MBZ = size( esp, 2 )
    MSP = size( esp, 3 )

    if ( nsetocc > 0 ) then
       mb1=1
       mb2=MB
       occ(:,:,:)=0.0d0
       nn=nsetocc/MSP
       do s=1,MSP
       do k=1,MBZ
          do n=1,isetocc(1)-1
             occ(n,k,s)=1.0d0
          end do
          occ(isetocc(1):isetocc(2),k,s) = setocc(1+nn*(s-1):nn+nn*(s-1))
       end do
       end do
       ef=maxval( esp(isetocc(2),1:MBZ,1:MSP) )
       goto 100
    end if

    if ( first_time ) call prep_ff0

    mb1 = 1
    mb2 = MB

! Set upper & lower boundarires of Fermi energy

    ef1 = minval( esp(mb1:mb2,1:MBZ,1:MSP) )
    ef2 = maxval( esp(mb1:mb2,1:MBZ,1:MSP) )
    if ( ef1 == ef2 ) then
       ef1 = ef1 - 0.01d0
       ef2 = ef2 + 0.01d0
    end if

!C Safety margin for highly degenerate systems & artificial fault
!C
    ef2 = ef2 + min( ekbt*1.d2, 0.1d0 )

    if ( MSP == 1 .or. iter > Nfixed ) then

       zne = znel - 2.d0*(mb1-1)

       ef0 = 1.d10

       do id=1,mxcycl

          ef = 0.5d0*( ef1 + ef2 )
          if ( ef == ef0 ) goto 100
          octmp = 0.0d0

          do s=1,MSP
          do k=1,MBZ
          do n=mb1,mb2

             xx=(esp(n,k,s)-ef)/ekbt
             ff=ff0(kinteg,xx)

             octmp = octmp + ff*wbz(k)*2.d0/dble(MSP)
             occ(n,k,s)=ff

          end do
          end do
          end do

          if ( octmp-zne > eps ) then
             ef2=ef
          else if ( octmp-zne < -eps ) then
             ef1=ef
          else
             goto 100
          end if

          ef0 = ef

       end do ! id

    else

       efconv=0
       if ( DISP_SWITCH ) then
          write(*,*) "total spin density is fixed!!"
       end if

       do s=1,MSP

          ef1 = minval( esp(mb1:mb2,1:MBZ,s) )
          ef2 = maxval( esp(mb1:mb2,1:MBZ,s) )
          if ( ef1 == ef2 ) then
             ef1 = ef1 - 0.01d0
             ef2 = ef2 + 0.01d0
          end if
          ef2 = ef2 + min( ekbt*1.d2, 0.1d0 )

          zne = 0.5d0*znel + (3-2*s)*0.5d0*dspn

          do id=1,mxcycl

             ef = 0.5d0*(ef1+ef2)
             octmp = 0.0d0

             do n=mb1,mb2

                do k=1,MBZ

                   xx=(esp(n,k,s)-ef)/ekbt
                   ff=ff0(kinteg,xx)

                   octmp=octmp+ff*wbz(k)
                   occ(n,k,s)=ff

                end do

             end do

             if ( octmp-zne > eps ) then
                ef2=ef
             else if ( octmp-zne < -eps ) then
                ef1=ef
             else
                efconv=efconv+1
                exit
             end if

          end do ! id

       end do ! s

       if ( efconv == 2 ) goto 100

    end if

    if ( abs(octmp-zne) > 1.d-10 ) then
       if ( disp_switch ) then
          write(6,*)' EF IS NOT CONVERGED'
          write(6,*)' Check the # of electron, mb1, and mb2'
          write(6,*)' EF1 & EF2=',ef1,ef2
          write(6,*)' octmp,zne=',octmp,zne,octmp-zne
          do s=1,MSP
             do k=1,MBZ
                write(*,*) "s,k =",s,k
                do n=mb1,mb2
                   write(*,*) n,occ(n,k,s),esp(n,k,s)
                end do
             end do
          end do
       end if
       stop'FERMI'
    end if

100 continue

    efermi   = ef
    Eentropy = 0.d0

    do s=1,MSP
    do k=1,MBZ
       do n=1,mb1-1
          occ(n,k,s)=2.d0*wbz(k)/dble(MSP)
       end do
       do n=mb1,mb2
          occ(n,k,s)=2.d0*occ(n,k,s)*wbz(k)/dble(MSP)
       end do
       do n=mb2+1,MB
          occ(n,k,s)=0.d0
       end do
    end do
    end do

    call write_border( 1, " calc_fermi_mp_1(end)" )

  END SUBROUTINE calc_fermi_mp_1


END MODULE fermi_module
