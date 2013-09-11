MODULE fermi_module

  implicit none

  PRIVATE
  PUBLIC :: read_fermi,ekbt,calc_fermi &
           ,efermi, Eentropy

  real(8) :: ekbt, efermi, Eentropy
  integer :: mb1,mb2

CONTAINS


  SUBROUTINE read_fermi(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(4) :: cbuf,ckey
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:4) == "EKBT" ) then
             backspace(unit)
             read(unit,*) cbuf,ekbt
          end if
       end do
999    continue
       write(*,*) "ekbt=",ekbt
    end if
    call send_fermi(0)
  END SUBROUTINE read_fermi


  SUBROUTINE send_fermi(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(ekbt,1,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_fermi


  SUBROUTINE calc_fermi(iter,nfixed,MB,MBZ,MSP,znel,dspn,esp,wbz,occ &
                       ,disp_switch)
!(This routine may assume esp's sorted)
    integer,intent(IN)  :: iter,nfixed,MB,MBZ,MSP
    logical,intent(IN)  :: disp_switch
    real(8),intent(IN)  :: esp(MB,MBZ,MSP),wbz(MBZ),znel,dspn
    real(8),intent(OUT) :: occ(MB,MBZ,MSP)
    real(8) :: ef1,ef2,ef
    integer :: id,m,n,k,s,efconv
    real(8) :: zne,octmp,ff,xx,ff0,c
    real(8),parameter :: eps=0.d0, big=36.9d0
    integer,parameter :: mxcycl=1000

    occ = 0.d0

    mb1 = 1
    mb2 = MB

! Set upper & lower boundarires of Fermi energy

    ef1=minval( esp(mb1:mb2,1:MBZ,1:MSP) )
    ef2=maxval( esp(mb1:mb2,1:MBZ,1:MSP) )

    if ( ef1==ef2 ) then
       ef1=ef1-0.01d0
       ef2=ef2+0.01d0
    end if

! Safety margin for highly degenerate systems & artificial fault

    ef2 = ef2 + min( ekbt*1.d2, 0.1d0 )

    if ( MSP==1 .or. iter>nfixed ) then

       zne = znel - 2.d0*(mb1-1)

       do id=1,mxcycl

          ef=0.5d0*(ef1+ef2)

          octmp=0.0d0
          do s=1,MSP
          do k=1,MBZ
          do n=mb1,mb2
!             if ( ekbt<eps ) then
!                if ( (esp(n,k,s)-ef )>0.0d0 ) then
!                   ff=0.0d0
!                else if( ( esp(n,k,s)-ef)<0.0d0 ) then
!                   ff=1.0d0
!                else
!                   ff=0.5d0
!                end if
!             else
                xx=( esp(n,k,s)-ef )/ekbt
                if ( xx>big ) then
                   ff=0.0d0
!                else if ( xx<-big ) then
!                   ff=1.0d0
                else
                   ff=1.0d0/(exp(xx)+1.0d0)
                end if
!             end if
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
       end do

    else

       efconv=0

       do s=1,MSP

          ef1=minval( esp(mb1:mb2,1:MBZ,s) )
          ef2=maxval( esp(mb1:mb2,1:MBZ,s) )

          if ( ef1==ef2 ) then
             ef1=ef1-0.01d0
             ef2=ef2+0.01d0
          end if
          ef2 = ef2 + min( ekbt*1.d2, 0.1d0 )

          zne=(znel+dble(3-2*s)*dspn)*0.5d0

          do id=1,mxcycl

             ef=0.5d0*(ef1+ef2)

             octmp=0.0d0
             do k=1,MBZ
             do n=mb1,mb2
!                if ( ekbt<eps ) then
!                   if ( (esp(n,k,s)-ef )>0.0d0 ) then
!                      ff=0.0d0
!                   else if( ( esp(n,k,s)-ef)<0.0d0 ) then
!                      ff=1.0d0
!                   else
!                      ff=0.5d0
!                   end if
!                else
                   xx=( esp(n,k,s)-ef )/ekbt
                   if ( xx>big ) then
                      ff=0.0d0
!                   else if ( xx<-big ) then
!                      ff=1.0d0
                   else
                      ff=1.0d0/(exp(xx)+1.0d0)
                   end if
!                end if
                octmp=octmp+ff*wbz(k)
                occ(n,k,s)=ff
             end do
             end do
             if ( octmp-zne>eps ) then
                ef2=ef
             else if ( octmp-zne<-eps ) then
                ef1=ef
             else
                efconv=efconv+1
                exit
             end if
          end do
       end do

       if ( efconv==2 ) goto 100

    end if

    if ( abs(octmp-zne) > 1.d-10 ) then
       if ( disp_switch ) then
          write(6,*)' EF IS NOT CONVERGED'
          write(6,*)' Check the # of electron, mb1, and mb2'
          write(6,*)' EF1 & EF2=',ef1,ef2
          write(6,*)' octmp,zne=',octmp,zne,octmp-zne
!          do s=1,nspin
!          do k=1,MBZ
!             write(*,*) "s,k =",s,k
!             do ib=mb1,mb2
!                write(*,*) ib,occ(ib,k,s),esp(ib,k,s)
!             end do
!          end do
!          end do
       end if
       stop 'FERMI'
    end if

100 continue

    efermi = ef

    Eentropy = 0.d0
    do s=1,MSP
    do k=1,MBZ
       c=wbz(k)*2.d0/MSP
       do n=1,MB
          ff=occ(n,k,s)
          if ( ff == 0.d0 .or. ff == 1.d0 ) cycle
          Eentropy=Eentropy-ekbt*( (1.d0-ff)*log(1.d0-ff)+ff*log(ff) )*c
       end do
    end do
    end do

    if ( disp_switch ) write(*,*) "Eentropy=",Eentropy

    do s=1,MSP
    do k=1,MBZ
       c=wbz(k)*2.d0/MSP
       do n=mb1,mb2
          occ(n,k,s) = occ(n,k,s)*c
       end do
    end do
    end do

!    octmp = sum(occ(1:MB,1:MBZ,1:MSP))
!    octmp = znel-octmp
!    if ( abs(octmp) > eps ) then
!       do s=1,MSP
!  loop1 : do n=mb2,mb1,-1
!          do k=1,MBZ
!             if ( occ(n,k,s)>0.d0 ) then
!                m=n
!                exit loop1
!             end if
!          end do
!          end do loop1
!          do k=1,MBZ
!             occ(m,k,s)=occ(m,k,s)+octmp/dble(MBZ*MSP)
!          end do
!       end do
!    end if

  END SUBROUTINE calc_fermi

END MODULE fermi_module
