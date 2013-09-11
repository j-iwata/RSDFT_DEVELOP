!--------1---------2---------3---------4---------5---------6---------7--
!
! Spherical Harmonic Funciton
!
      SUBROUTINE prep_hartree1
      use global_variables
      implicit none
      integer :: i,j,k,lm,L,a,m,n,n1,n2,ML0
      integer :: m1,m2,MK0
      logical :: flag_alloc(2)
      real(8),parameter :: eps=1.d-20
      real(8) :: r,x,y,z,const,mem,memax,H
      real(8) :: ctime0,ctime1,etime0,etime1,ct(2),et(2)

      INTERFACE
         FUNCTION Ylm(x,y,z,l,m)
         real(8) :: Ylm
         real(8),intent(IN) :: x,y,z
         integer,intent(IN) :: l,m
         END FUNCTION Ylm
      END INTERFACE

      call watch(ctime0,etime0) ; ct=0.d0 ; et=0.d0
      if (DISP_SWITCH) write(*,'(a60," prep_shf")') repeat("-",60)

      n1  = idisp(myrank)+1
      n2  = idisp(myrank)+ircnt(myrank)
      ML0 = ircnt(myrank)
      m1  = idisp2(myrank)+1+ML_irreducible
      m2  = idisp2(myrank)+ircnt2(myrank)
      MK0 = ircnt2(myrank)

      lmmax_HT=(Lmax_HT+1)**2

      H = H1

      mem=0.d0
      memax=0.d0

      flag_alloc(1:2)=.false.

      if (DISP_SWITCH) then
         write(*,*) "Lmax_HT  =",Lmax_HT
         write(*,*) "lmmax_HT =",lmmax_HT
      end if

!- allocate ------------------------------
      call gv_alloc("prep_hartree1")
!-----------------------------------------

!- allocate ----------------------------------------------------
      if ( .not.allocated(LL) ) then
         flag_alloc(1)=.true.
         allocate( LL(3,n1:n2) ) ; LL=0
         mem=mem+bsintg*size(LL) ; memax=max(mem,memax)
         call Make_GridMap_1(LL,n1,n2)
      end if
      if ( .not.allocated(KK) ) then
         flag_alloc(2)=.true.
         allocate( KK(3,m1:m2) ) ; KK=0
         mem=mem+bsintg*size(KK) ; memax=max(mem,memax)
         call Make_GridMap_2(KK,m1,m2)
      end if
!---------------------------------------------------------------


      do i=n1,n2
         x=LL(1,i)*H ; y=LL(2,i)*H ; z=LL(3,i)*H
         r=sqrt(x*x+y*y+z*z)
         lm=0
         do L=0,Lmax_HT
            do m=-L,L
               lm=lm+1
               shf1(i,lm)=Ylm(x,y,z,l,m)*r**L
            end do
         end do
      end do

      do i=m1,m2
         x=KK(1,i)*H ; y=KK(2,i)*H ; z=KK(3,i)*H
         r=sqrt(x*x+y*y+z*z)
         lm=0
         do L=0,Lmax_HT
            const=4.d0*Pi/(2.d0*L+1.d0)
            do m=-L,L
               lm=lm+1
               shf2(lm,i)=Ylm(x,y,z,l,m)/r**(L+1)*const
            end do
         end do
      end do

!- deallocate -----------------------------------
      if ( flag_alloc(2) ) then
         mem=mem-bsintg*size(KK)
         deallocate( KK )
      end if
      if ( flag_alloc(1) ) then
         mem=mem-bsintg*size(LL)
         deallocate( LL )
      end if
!------------------------------------------------

      call watch(ctime1,etime1)
      if (DISP_SWITCH) then
         write(*,*) "TIME(PREP_HARTREE1)=",ctime1-ctime0,etime1-etime0
         write(*,*) "MEO=",MEO
         write(*,*) "MEM(MB)",memax*B2MB,mem
      end if

      return
      END SUBROUTINE prep_hartree1

!--------1---------2---------3---------4---------5---------6---------7--
!
! Division of Area for Multipole Expansion (MEO=2)
!
      SUBROUTINE prep_hartree2
      use global_variables
      implicit none
      real(8),allocatable :: ra(:)
      real(8) :: ctime0,ctime1,etime0,etime1,x,y,z,r2
      real(8) :: mem,memax,H
      integer :: i,a,m,n,i1,i2,i3,n1,n2,ierr
      integer,allocatable :: itmp(:),jtmp(:)

      if ( isymmetry==1 ) then
         call prep_hartree2_sym
         return
      end if

      call watch(ctime0,etime0)

      if (DISP_SWITCH) write(*,'(a60," prep_hartree2")') repeat("-",60)

      n1=idisp(myrank)+1
      n2=idisp(myrank)+ircnt(myrank)

      H = H1

      mem=0.d0
      memax=0.d0

!- allocate -------------------------------------------------
      n=maxval(ircnt)
      allocate( itmp(n1:n2) ) ; itmp=0
      allocate( jtmp(MI)    ) ; jtmp=0
      allocate( LL(3,n1:n2) ) ; LL=0
      allocate( ra(n1:n2)   ) ; ra=1.d10
      mem=mem+bsintg*(size(itmp)+size(jtmp))
      mem=mem+bsintg*size(LL)
      mem=mem+bdreal*size(ra)
      memax=max(mem,memax)
!------------------------------------------------------------

      call Make_GridMap_1(LL,n1,n2)

      do a=1,MI
      do i=n1,n2
         x=LL(1,i)*H-asi(1,a)
         y=LL(2,i)*H-asi(2,a)
         z=LL(3,i)*H-asi(3,a)
         r2=x*x+y*y+z*z
         if ( r2<ra(i) ) then
            ra(i)=r2
            itmp(i)=a
         end if
      end do
      end do

!- deallocate -----------------------------------
      mem=mem-bdreal*size(ra) ; deallocate( ra )
!------------------------------------------------

      do a=1,MI
         jtmp(a)=count(itmp==a) ! # of grid points near the atom a
      end do

      NMadv=count( jtmp>0 )     ! # of atoms (or # of regions) in my rank
      maxMdv=maxval(jtmp)       ! max # of grid points around each atom

      if ( DISP_SWITCH ) then
         write(*,*) "NMadv,maxMdv,ML,sum(jtmp)=",NMadv,maxMdv,ML,sum(jtmp)
      end if

!- allocate ---------------------------
! Ixyz,Mdv,adv
      call gv_alloc("prep_hartree2")
!--------------------------------------

      n=0
      do a=1,MI
         if ( jtmp(a)>0 ) then
            n=n+1
            Mdv(n)=jtmp(a)
            adv(n)=a
         end if
      end do

      if ( n/=NMadv ) call stop_program

      do n=1,NMadv
         a=adv(n)
         m=0
         do i=n1,n2
            if ( itmp(i)==a ) then
               m=m+1
               Ixyz(m,n)=i
            end if
         end do
      end do

!- deallocate ---------------------------------------------------------
      mem=mem-bsintg*size(LL) ; deallocate(LL)
      mem=mem-bsintg*(size(jtmp)+size(itmp)) ; deallocate( jtmp,itmp )
!----------------------------------------------------------------------

      call watch(ctime1,etime1)
      if (DISP_SWITCH) then
         write(*,*) "TIME(PREP_HARTREE2)=",ctime1-ctime0,etime1-etime0
         write(*,*) "MEO=",MEO
         write(*,*) "MEM(MB)",memax*B2MB,mem
      end if

      return
      END SUBROUTINE prep_hartree2
