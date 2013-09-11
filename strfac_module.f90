MODULE strfac_module

  use Ggrid_module
  use atom_module

  implicit none

  PRIVATE
  PUBLIC :: SGK, construct_strfac, destruct_strfac

  complex(8),allocatable :: SGK(:,:)

CONTAINS

  SUBROUTINE construct_strfac

    integer :: a,i,i1,i2,i3,ik,j,ierr,MG
    real(8) :: Gr,pi2,a1,a2,a3

    pi2 = 2.d0*acos(-1.d0)
    MG  = NGgrid(0)

    allocate( SGK(MG,Nelement) ) ; SGK=(0.d0,0.d0)

    call construct_Ggrid(1)

    do a=1,Natom
       ik=ki_atom(a)
       a1=pi2*aa_atom(1,a)
       a2=pi2*aa_atom(2,a)
       a3=pi2*aa_atom(3,a)
       do i=MG_0,MG_1
          Gr=LLG(1,i)*a1+LLG(2,i)*a2+LLG(3,i)*a3
          SGK(i,ik)=SGK(i,ik)+dcmplx(cos(Gr),-sin(Gr))
       end do
    end do

    call destruct_Ggrid

  END SUBROUTINE construct_strfac

  SUBROUTINE destruct_strfac
    deallocate( SGK )
  END SUBROUTINE destruct_strfac

END MODULE strfac_module
