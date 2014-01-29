MODULE esm_ion_module

  use aa_module
  use atom_module
  use pseudopot_module, only: Zps
  use modified_bessel_module
  use esm_rgrid_module, only: Rsize2

  implicit none

  PRIVATE
  PUBLIC :: calc_energy_esm_ion

CONTAINS

  SUBROUTINE calc_energy_esm_ion
    implicit none
    integer :: t1,t2,ma,ikz,ima,ikz_max,ima_max
    real(8) :: pi2,kz,mp,x1,x2,y1,y2,z1,z2,r1,r2,p1,p2,dk,k,rl,rr
    real(8) :: qq,const,ka
    complex(8) :: phase_k,phase_m,ze

    pi2=2.d0*acos(-1.0d0)
    dk=pi2/aa(3,3)
    const=-1.d0/aa(3,3)

    ikz_max=30
    ima_max=5

    ze=(0.d0,0.d0)
    do t2=1,Natom
       x2=aa(1,1)*aa_atom(1,t2)+aa(1,2)*aa_atom(2,t2)+aa(1,3)*aa_atom(3,t2)
       y2=aa(2,1)*aa_atom(1,t2)+aa(2,2)*aa_atom(2,t2)+aa(2,3)*aa_atom(3,t2)
       z2=aa(3,1)*aa_atom(1,t2)+aa(3,2)*aa_atom(2,t2)+aa(3,3)*aa_atom(3,t2)
       r2=sqrt(x2*x2+y2*y2)
       if ( r2 == 0.0d0 ) then
          p2 = 0.0d0
       else
          p2 = sign( acos(x2/r2),y2 )
       end if
       do t1=1,Natom
!          if ( t1 == t2 ) cycle
          x1=aa(1,1)*aa_atom(1,t1)+aa(1,2)*aa_atom(2,t1)+aa(1,3)*aa_atom(3,t1)
          y1=aa(2,1)*aa_atom(1,t1)+aa(2,2)*aa_atom(2,t1)+aa(2,3)*aa_atom(3,t1)
          z1=aa(3,1)*aa_atom(1,t1)+aa(3,2)*aa_atom(2,t1)+aa(3,3)*aa_atom(3,t1)
          r1=sqrt(x1*x1+y1*y1)
          qq=Zps(ki_atom(t1))*Zps(ki_atom(t2))*const
          if ( r1 == 0.0d0 ) then
             p1 = 0.0d0
          else
             p1 = sign( acos(x1/r1),y1 )
          end if
          if ( r1 <= r2 ) then
             rl=r1
             rr=r2
          else
             rl=r2
             rr=r1
          end if
          do ikz=-ikz_max,ikz_max
             k=dk*ikz
             kz=k*(z2-z1)
             phase_k = dcmplx(cos(kz),sin(kz))
             ka=abs(k)
          do ima=-ima_max,ima_max
             mp=ima*(p2-p1)
             phase_m = dcmplx(cos(mp),sin(mp))
             ma=abs(ima)
             ze=ze+qq*phase_k*phase_m*gm(ma,ka,rl,rr)
             write(*,*) ikz,ima,gm(ma,ka,rl,rr)
          end do ! ima
          end do ! ikz
       end do ! t1
    end do ! t2

    write(*,*) "ze=",ze

  END SUBROUTINE calc_energy_esm_ion

  FUNCTION gm(m,k,rl,rr)
    implicit none
    real(8) :: gm
    real(8),intent(IN) :: k,rl,rr
    integer,intent(IN) :: m
    real(8),parameter :: ktmp=1.d-14,rtmp=1.d-14
    if ( k <= 1.d-14 ) then
       if ( rr <= 1.d-14 ) then
          gm=bessi(m,ktmp)*( &
               bessk(m,ktmp*Rsize2)/bessi(m,ktmp*Rsize2)*bessi(m,ktmp) ) &
              -bessi(m,ktmp)*bessk(m,ktmp)
       else
          gm=bessi(m,ktmp*rl) &
               *( bessk(m,ktmp*Rsize2)/bessi(m,ktmp*Rsize2)*bessi(m,ktmp*rr) &
                 -bessk(m,ktmp*rr) )
       end if
    else
       if ( rr <= 1.d-14 ) then
          gm=bessi(m,k*rtmp)*( &
               bessi(m,k*Rsize2)/bessi(m,k*Rsize2)*bessk(m,k*rtmp) ) &
               -bessi(m,k*rtmp)*bessk(m,k*rtmp)
       else
          gm=bessi(m,k*rl) &
               *( bessi(m,k*Rsize2)/bessi(m,k*Rsize2)*bessk(m,k*rr) &
                 -bessk(m,k*rr) )
       end if
    end if
  END FUNCTION gm

END MODULE esm_ion_module
