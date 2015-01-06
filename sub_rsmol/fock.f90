!--------1---------2---------3---------4---------5---------6---------7--

SUBROUTINE Fock(n1,n2,ib,psi,tpsi)
  use global_variables
  implicit none
  integer :: n1,n2,ib,k,s,i,n
#ifdef DRSDFT
  real(8) :: psi(n1:n2),tpsi(n1:n2)
#else
  complex(8) :: psi(n1:n2),tpsi(n1:n2)
#endif
  intent(IN) :: n1,n2,ib,psi
  intent(INOUT) :: tpsi
  real(8),allocatable :: trho(:),tvht(:)
  real(8) :: mem,memax,ctime0,ctime1,etime0,etime1,c
  logical :: DISP_SWITCH_TMP

!  call bwatch(ctime0,etime0)

  mem=0.d0
  memax=0.d0

  DISP_SWITCH_TMP=DISP_SWITCH
  DISP_SWITCH=.false.

  k=1
  s=1

  select case(iflag_hf)
  case default

     do i=n1,n2
        tpsi(i)=tpsi(i)+VFunk(i,ib,k,s)
     end do

  case(2)

     allocate( trho(n1:n2) ) ; trho=0.d0 ; mem=mem+bdreal*size(trho) ; memax=max(mem,memax)
     allocate( tvht(n1:n2) ) ; tvht=0.d0 ; mem=mem+bdreal*size(tvht) ; memax=max(mem,memax)
 
     do n=1,MB

        if ( abs(occ(n,k,s))<1.d-10 ) cycle

        do i=n1,n2
!           trho(i)=unk(i,n,k,s)*psi(i)
           trho(i)=unk_hf(i,n,k,s)*psi(i)
        end do

        call Hartree_mol(n1,n2,trho,tvht,1.d-25,2000,0)

        c=occ(n,k,s)*0.5d0
        do i=n1,n2
!           tpsi(i)=tpsi(i)-c*tvht(i)*unk(i,n,k,s)
           tpsi(i)=tpsi(i)-c*tvht(i)*unk_hf(i,n,k,s)
        end do

     end do

     mem=mem-bdreal*size(tvht) ; deallocate( tvht )
     mem=mem-bdreal*size(trho) ; deallocate( trho )

  case(22)

     allocate( trho(n1:n2) ) ; trho=0.d0 ; mem=mem+bdreal*size(trho) ; memax=max(mem,memax)
     allocate( tvht(n1:n2) ) ; tvht=0.d0 ; mem=mem+bdreal*size(tvht) ; memax=max(mem,memax)
 
     do n=1,MB

        if ( abs(occ(n,k,s))<1.d-10 ) cycle

        do i=n1,n2
           trho(i)=unk(i,n,k,s)*psi(i)
!           trho(i)=unk_hf(i,n,k,s)*psi(i)
        end do

        call Hartree_mol(n1,n2,trho,tvht,1.d-25,2000,0)

        c=occ(n,k,s)*0.5d0
        do i=n1,n2
           tpsi(i)=tpsi(i)-c*tvht(i)*unk(i,n,k,s)
!           tpsi(i)=tpsi(i)-c*tvht(i)*unk_hf(i,n,k,s)
        end do

     end do

     mem=mem-bdreal*size(tvht) ; deallocate( tvht )
     mem=mem-bdreal*size(trho) ; deallocate( trho )

  case(3)

     allocate( trho(n1:n2) ) ; trho=0.d0 ; mem=mem+bdreal*size(trho) ; memax=max(mem,memax)
     allocate( tvht(n1:n2) ) ; tvht=0.d0 ; mem=mem+bdreal*size(tvht) ; memax=max(mem,memax)
 
     do n=1,MB

        if ( abs(occ(n,k,s))<1.d-10 ) cycle

        do i=n1,n2
           trho(i)=unk_hf(i,n,k,s)*psi(i)
        end do

        call Hartree_mol(n1,n2,trho,tvht,1.d-25,2000,0)

        c=beta*occ(n,k,s)*0.5d0
        do i=n1,n2
           tpsi(i)=tpsi(i)-c*tvht(i)*unk_hf(i,n,k,s)
        end do

        do i=n1,n2
           trho(i)=unk_hf1(i,n,k,s)*psi(i)
        end do

        call Hartree_mol(n1,n2,trho,tvht,1.d-25,2000,0)

        c=(1.d0-beta)*occ(n,k,s)*0.5d0
        do i=n1,n2
           tpsi(i)=tpsi(i)-c*tvht(i)*unk_hf1(i,n,k,s)
        end do

     end do

     mem=mem-bdreal*size(tvht) ; deallocate( tvht )
     mem=mem-bdreal*size(trho) ; deallocate( trho )

  end select
   
  DISP_SWITCH=DISP_SWITCH_TMP

!  call bwatch(ctime1,etime1)

!  if (DISP_SWITCH) then
!     write(*,*) "TIME(FOCK)=",ctime1-ctime0,etime1-etime0
!     write(*,*) " MEM(MB)=",memax*B2MB,mem
!  end if

  return

END SUBROUTINE Fock
  
