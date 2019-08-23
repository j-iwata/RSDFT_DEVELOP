module allel_module

  use atom_module, only: zn_atom
  use ggrid_module, only: GG
  use ps_local_variables, only: vqlg
  use io_tools_module

  implicit none

  private
  public :: init_allel
  public :: init_ae_local_allel

  logical,public :: flag_allel=.false.

  real(8) :: Vcell

contains

  subroutine init_allel( Zps, Vcell_in )
    implicit none
    real(8),intent(out) :: Zps(:)
    real(8),intent(in) ::Vcell_in
    call IOTools_findKeyword( "ALLEL", flag_allel, flag_bcast=.true. )
    if ( flag_allel ) then
       Zps=zn_atom
       Vcell=Vcell_in
    end if
  end subroutine init_allel

  subroutine init_ae_local_allel
    implicit none
    integer :: Nelm, ig, NMGL, ik
    real(8) :: pi, const, G2
    NMGL=size( GG )
    Nelm=size( zn_atom )
    allocate( vqlg(NMGL,Nelm) ); vqlg=0.0d0
    pi=acos(-1.0d0)
    do ik=1,Nelm
       const=-4.0d0*pi/Vcell*Zn_atom(ik)
       do ig=1,NMGL
          G2=GG(ig)
          if ( G2 == 0.0d0 ) cycle
          vqlg(ig,ik)=vqlg(ig,ik)+const/G2
       end do
    end do
  end subroutine init_ae_local_allel

end module allel_module
