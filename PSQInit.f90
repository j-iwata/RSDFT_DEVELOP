MODULE PSQInit
  use VarPSQG

  PRIVATE
  PUBLIC :: 

  implicit none

CONTAINS

!-----------------------------------------------
  SUBROUTINE initKtoKPSQ()
    implicit none
    integer :: ik,l1,i1,m1
    integer :: k1,nr1,nr2
    
    call allocateKtoK( k1max,Nelement,Rrefmax,Lrefmax )

    do ik=1,Nelement
        k1=0
        nr1=0
        do l1=1,nl(ik)
            do i1=1,nr(l1,ik)
                nr1=nr1+1
                do m1=1,2*l1-1
                    nr2=0
                    do 

  END SUBROUTINE initKtoKPSQ


END MODULE PSQInit
