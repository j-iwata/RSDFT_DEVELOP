module init_occ_electron_module

  implicit none

  private
  public :: init_occ_electron

contains

  subroutine init_occ_electron(Nelectron,Ndspin,weight_bz,occ)
    implicit none
    real(8),intent(in) :: Nelectron,Ndspin,weight_bz(Nbzsm)
    real(8),intent(out) :: occ(:,:,:)
    integer :: n,k,s,Nspin,Nbzsm,Nband
    real(8) :: sum0,d,Nel
    logical :: disp

    if ( all(weight_bz==0.0d0) ) return

    call write_border( 0, " init_occ_electron(start)" )

    call check_disp_switch( disp, 0 )

    occ=0.0d0

    Nband = size(occ,1)
    Nbzsm = size(occ,1)
    Nspin = size(occ,3)

    d=2.0d0/Nspin
    do s=1,Nspin
      Nel = 0.5d0*d*Nelectron + (3-2*s)*0.5d0*Ndspin
      if ( Nel < 0.0d0 ) stop "Ndspin is too large !!!"
      sum0=0.0d0
      do n=1,Nband
        if ( sum0+d > Nel ) exit
        sum0=sum0+d
        occ(n,1,s)=d
      end do
      if ( sum0 < Nel ) then
        if ( n > Nband ) then
          call stop_program( "Nband is small (init_occupation_electron)" )
        else
          occ(n,1,s) = Nel-sum0
          if ( disp ) then
            write(*,*) "Nel, sum0 =",Nel,sum0
            write(*,*) "occ,ib_vbm,s=",occ(n,1,s),n,s
          end if
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

    call write_border( 0, " init_occ_electron(end)" )

  end subroutine init_occ_electron

end module init_occ_electron_module
