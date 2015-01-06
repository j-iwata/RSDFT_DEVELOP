!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE calc_virial
      use global_variables
      implicit none
      integer :: n1,n2,i,m,ierr,i1,i2,i3
      real(8),allocatable :: wt0(:,:,:),wt(:,:,:),vt(:)
      real(8) :: Tc,sum0,x,y,z

      n1=idisp(myrank)+1
      n2=idisp(myrank)+ircnt(myrank)

      allocate( wt(-Mx:Mx,-My:My,-Mz:Mz) ) ; wt=0.d0
      allocate( wt0(-Mx:Mx,-My:My,-Mz:Mz) ) ; wt0=0.d0

      do i=n1,n2
         wt0(LL(1,i),LL(2,i),LL(3,i))=Vxc(i)
      end do

      m=size(wt)
      call mpi_allreduce(wt0,wt,m,mpi_real8,mpi_sum,mpi_comm_world,ierr)

      deallocate( wt0 )
      allocate( vt(n1:n2) ) ; vt=0.d0

      do i=n1,n2
         i1=LL(1,i) ; i2=LL(2,i) ; i3=LL(3,i)
         x=i1*H
         do m=1,Md
            vt(i)=vt(i)+nab(m)*( x*( wt(i1+m,i2,i3)-wt(i1-m,i2,i3) ) & 
                                +y*( wt(i1,i2+m,i3)-wt(i1,i2-m,i3) ) &
                                +z*( wt(i1,i2,i3+m)-wt(i1,i2,i3-m) ) )
         end do
      end do

      sum0=sum( rho(n1:n2)*vt(n1:n2) )*dV
      call mpi_allreduce(sum0,Tc,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

      Tc=-Exc-Tc

      if ( DISP_SWITCH ) then
         write(*,*) "Tc=",Tc
         write(*,*) "Exc=",Exc
         write(*,*) "Tc+Exc=",Tc+Exc
      end if

      return
      END SUBROUTINE calc_virial

!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE calc_selfene
      use global_variables
      integer :: n1,n2,n,iloop
      real(8),allocatable :: occ_tmp(:)
      real(8) :: ff

      if ( ifac_vloc/=0 ) goto 900

      allocate( occ_tmp(MST) )
      occ_tmp(1:MST)=occ(1:MST)

      if ( myrank==0 ) rewind 100

      do iloop=1,1

!      ff=0.1d0*iloop
!      occ(2)=2.d0-ff/3.d0
!      occ(3)=2.d0-ff/3.d0
!      occ(4)=2.d0-ff/3.d0

      call psi_rho
      esp(:)=0.d0
      call Total_Energy

      if ( myrank==0 ) then
         write(100,'(1x,f8.4,7f15.8)') -ff,Etot,E_hartree,Ex,Ec,Exc,Ekin,Enl
      end if

      end do

      deallocate( occ_tmp )

      return

 900  write(*,*) "ifac_vloc=",ifac_vloc," stop at calc_selfint"
      call stop_program

      END SUBROUTINE calc_selfene

!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE calc_selfint
      use global_variables
      integer :: n1,n2,n,maxiter
      real(8),allocatable :: vr_tmp(:,:),e_tmp(:),e_si(:,:)
      real(8) :: tolconv

      n1      = idisp(myrank)+1
      n2      = idisp(myrank)+ircnt(myrank)
      tolconv = 1.d-25
      maxiter = 2000

      if ( ifac_vloc/=0 ) goto 900

      allocate( vr_tmp(n1:n2,3),e_tmp(4),e_si(5,MST) )
      vr_tmp(n1:n2,1)=rho(n1:n2)
      vr_tmp(n1:n2,2)=Vh(n1:n2)
      vr_tmp(n1:n2,3)=Vxc(n1:n2)
      e_tmp(1)=E_hartree
      e_tmp(2)=Exc
      e_tmp(3)=Ex
      e_tmp(4)=Ec
      e_si(:,:)=0.d0

      if ( myrank==0 ) rewind 100

      do iloop=1,1 !20

!      ff=0.1d0*iloop
!      occ(2)=2.d0-ff/3.d0
!      occ(3)=2.d0-ff/3.d0
!      occ(4)=2.d0-ff/3.d0
!      occ(4)=2.d0-ff

      sum0=0.d0
      do n=1,MST
         rho(n1:n2)=occ(n)*psi(n1:n2,n)*psi(n1:n2,n)
         Vh=0.d0
         call hartree(rho,Vh,tolconv,maxiter,0)
         call Exc_Cor
         e_si(1,n)=E_hartree
         e_si(2,n)=Exc
         e_si(3,n)=Ex
         e_si(4,n)=Ec
         e_si(3,n)=sum(Vion(n1:n2)*rho(n1:n2))*dV
         sum0=sum0+E_hartree+Exc
      end do
      if (DISP_SWITCH) then
!         do n=1,MST
!            write(*,'(1x,i6,2x,5f15.8)') n,e_si(1:5,n)
!         end do
         write(*,*) "Total SI =",sum0
      end if

      if ( myrank==0 ) write(100,'(1x,f8.4,12f15.8)') znext,e_si(1:3,1:4)

      end do

      rho(n1:n2)=vr_tmp(n1:n2,1)
      Vh(n1:n2) =vr_tmp(n1:n2,2)
      Vxc(n1:n2)=vr_tmp(n1:n2,3)
      E_hartree=e_tmp(1)
      Exc=e_tmp(2)
      Ex=e_tmp(3)
      Ec=e_tmp(4)

      deallocate( e_tmp,vr_tmp )

      return

 900  write(*,*) "ifac_vloc=",ifac_vloc," stop at calc_selfint"
      call stop_program

      END SUBROUTINE calc_selfint
