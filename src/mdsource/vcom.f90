!-----------------------------------------------------------------------
!     Set center of motion off
!-----------------------------------------------------------------------
subroutine vcom(kine0)
  use atom_module, only: Natom,ki_atom,zn_atom
  use cpmd_variables, only: pmass,AMU,Velocity
  implicit none
  real(8),intent(out)   :: kine0
  real(8) rescale,pm,kine,vw(4)
  integer i,k
  vw(:)=0.d0
  do i=1,Natom
     pm=pmass(zn_atom(ki_atom(i)))*AMU
     do k=1,3
        vw(k)=vw(k)+Velocity(k,i)*pm
     enddo
     vw(4)=vw(4)+pm
  enddo
  call calkin(kine0)
  vw(1:3)=vw(1:3)/vw(4)
  do i=1,Natom
     do k=1,3
        Velocity(k,i)=Velocity(k,i)-vw(k)
     enddo
  end do
  call calkin(kine)
  rescale=sqrt(kine0/kine)
  Velocity(:,:)=Velocity(:,:)*rescale
  return
end subroutine vcom
