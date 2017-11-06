MODULE init_occ_electron_ncol_module

  implicit none

  PRIVATE
  PUBLIC :: init_occ_electron_ncol

CONTAINS

  SUBROUTINE init_occ_electron_ncol(Nelectron,Ndspin,Nbzsm,weight_bz,occ)
    implicit none
    integer,intent(IN) :: Nbzsm
    real(8),intent(IN) :: Nelectron,Ndspin,weight_bz(Nbzsm)
    real(8),intent(OUT) :: occ(:,:,:)
    integer :: n,k,s,Nspin,Nband
    real(8) :: sum0,d,Nel

    if ( all(weight_bz==0.0d0) ) return

    call write_border( 0, " init_occ_electron_ncol(start)" )

    occ=0.0d0

    Nband = size(occ,1)
    Nspin = size(occ,3)

    Nel = Nelectron
    d   = 1.0d0
    if ( Nel < 0.0d0 ) stop "Ndspin is too large !!!"

    do s=1,1
       sum0=0.0d0
       do n=1,Nband
          if ( sum0+d > Nel ) exit
          sum0=sum0+d
          occ(n,1,s)=d
       end do
       if ( sum0 < Nel ) then
          if ( n > Nband ) then
             stop "Nband is small (init_occupation)"
          else
             occ(n,1,s) = Nel-sum0
             write(*,*) Nel,sum0
          end if
       end if
       do k=2,Nbzsm
          do n=1,Nband
             occ(n,k,s)=occ(n,1,s)*weight_bz(k)
          end do
       end do
       do n=1,Nband
          occ(n,1,s)=occ(n,1,s)*weight_bz(1)
       end do
    end do
    sum0=sum(occ)
    if ( abs(sum0-Nelectron)>1.d-10 ) then
       write(*,'(1x,"sum(occ), Nelectron =",2g30.20)') sum(occ),Nelectron
       stop "sum(occ) /= Nelectron !!!"
    end if
    occ(:,:,Nspin)=occ(:,:,1)

    call write_border( 0, " init_occ_electron_ncol(end)" )

  END SUBROUTINE init_occ_electron_ncol

END MODULE init_occ_electron_ncol_module
