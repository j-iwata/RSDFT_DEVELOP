MODULE fock_module

  use xc_hybrid_module
  use bb_module
  use bz_module
  use electron_module
  use array_bound_module
  use wf_module
  use fock_fft_module

  implicit none

  PRIVATE
  PUBLIC :: Fock, op_fock

CONTAINS

  SUBROUTINE op_fock(k,s,n1,n2,ib1,ib2,tpsi,htpsi)
    implicit none
    integer,intent(IN) :: k,s,n1,n2,ib1,ib2
    complex(8),intent(IN) :: tpsi(n1:n2,ib1:ib2)
    complex(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    integer :: ib
    if ( iflag_hse == 0 ) return
    do ib=ib1,ib2
       call Fock(k,s,n1,n2,ib,tpsi(n1,ib),htpsi(n1,ib))
    end do
  END SUBROUTINE op_fock

  SUBROUTINE Fock(k,s,n1,n2,ib,psi,tpsi)
    implicit none
    integer :: k,s,n1,n2,ib
    integer :: q,tr,i,n
    real(8) :: k_fock(3),q_fock(3)
    complex(8) :: psi(n1:n2),tpsi(n1:n2)
    intent(IN) :: k,s,n1,n2,ib,psi
    intent(INOUT) :: tpsi
    complex(8),allocatable :: trho(:),tvht(:)
    real(8) :: mem,memax,ctime0,ctime1,etime0,etime1,c
    complex(8),parameter :: zero=(0.0d0,0.0d0)

    if ( iflag_pbe0   > 0 ) iflag_hf=iflag_pbe0
    if ( iflag_hse    > 0 ) iflag_hf=iflag_hse
    if ( iflag_lcwpbe > 0 ) iflag_hf=iflag_lcwpbe

    select case(iflag_hf)
    case default

       do i=n1,n2
          tpsi(i)=tpsi(i)+VFunk(i,ib,k,s)
       end do

    case(2) ! Case of SCF iterations
  
       k_fock(:)=bb(:,1)*kbb(1,k)+bb(:,2)*kbb(2,k)+bb(:,3)*kbb(3,k)

       allocate( trho(n1:n2) ) ; trho=zero
       allocate( tvht(n1:n2) ) ; tvht=zero

       do tr=0,tr_switch ! For time-reversal symmetry
 
          do q=1,MBZ

             q_fock(:)=bb(:,1)*kbb(1,q)+bb(:,2)*kbb(2,q)+bb(:,3)*kbb(3,q)
             if ( tr == 1 ) then
                q_fock(:) = -q_fock(:)
             end if

             do n=1,MB

                if ( abs(occ(n,q,s)) < 1.d-10 ) cycle

                c=occ_factor*occ(n,q,s)*alpha_hf

                if ( tr == 0 ) then
#ifndef _DRSDFT_
                   do i=n1,n2
                      trho(i)=conjg(unk_hf(i,n,q,s))*psi(i)
                   end do
#endif
                else
                   do i=n1,n2
                      trho(i)=unk_hf(i,n,q,s)*psi(i)
                   end do
                end if

                call Fock_fft(n1,n2,k_fock,q_fock,trho,tvht,tr)
!                call Fock_fft_parallel(n1,n2,k_fock,q_fock,trho,tvht,tr)

                if ( tr == 0 ) then
                   do i=n1,n2
                      tpsi(i)=tpsi(i)-c*tvht(i)*unk_hf(i,n,q,s)
                   end do
                else
#ifndef _DRSDFT_
                   do i=n1,n2
                      tpsi(i)=tpsi(i)-c*tvht(i)*conjg(unk_hf(i,n,q,s))
                   end do
#endif
                end if

             end do ! n

          end do ! q

       end do ! tr

       deallocate( trho )
       deallocate( tvht )

    case(22) ! Case of 1st SCF iteraction

       k_fock(:)=bb(:,1)*kbb(1,k)+bb(:,2)*kbb(2,k)+bb(:,3)*kbb(3,k)

       allocate( trho(n1:n2) ) ; trho=zero
       allocate( tvht(n1:n2) ) ; tvht=zero

       do tr=0,tr_switch ! For time-reversal symmetry

          do q=1,MBZ

             q_fock(:)=bb(:,1)*kbb(1,q)+bb(:,2)*kbb(2,q)+bb(:,3)*kbb(3,q)
             if ( tr == 1 ) then
                q_fock(:) = -q_fock(:)
             end if

             do n=1,MB

                if ( abs(occ(n,q,s)) < 1.d-10 ) cycle

                c=occ_factor*occ(n,q,s)

                if ( tr == 0 ) then
#ifndef _DRSDFT_
                   do i=n1,n2
                      trho(i)=conjg(unk(i,n,q,s))*psi(i)
                   end do
#endif
                else
                   do i=n1,n2
                      trho(i)=unk(i,n,q,s)*psi(i)
                   end do
                end if

                call Fock_fft(n1,n2,k_fock,q_fock,trho,tvht,tr)
!                call Fock_fft_parallel(n1,n2,k_fock,q_fock,trho,tvht,tr)

                if ( tr == 0 ) then
                   do i=n1,n2
                      tpsi(i)=tpsi(i)-c*tvht(i)*unk(i,n,q,s)
                   end do
                else
#ifndef _DRSDFT_
                   do i=n1,n2
                      tpsi(i)=tpsi(i)-c*tvht(i)*conjg(unk(i,n,q,s))
                   end do
#endif
                end if

             end do ! n

          end do ! q

       end do ! tr

       deallocate( trho )
       deallocate( tvht )

    case(3) ! Case of Band Calculation

       k_fock(:)=bb(:,1)*kbb(1,k)+bb(:,2)*kbb(2,k)+bb(:,3)*kbb(3,k)

       allocate( trho(n1:n2) ) ; trho=zero
       allocate( tvht(n1:n2) ) ; tvht=zero

       do tr=0,tr_switch ! For time-reversal symmetry
           
          do q=1,MBZ_fock

             q_fock(:)=bb(:,1)*kbb_fock(1,q)+bb(:,2)*kbb_fock(2,q)+bb(:,3)*kbb_fock(3,q)
             if ( tr == 1 ) then
                q_fock(:) = -q_fock(:)
             end if

             do n=1,MB_fock

                if ( abs(occ_fock(n,q,s)) < 1.d-10 ) cycle

                c=occ_factor*occ_fock(n,q,s)

                if ( tr == 0 ) then
#ifndef _DRSDFT_
                   do i=n1,n2
                      trho(i)=conjg(unk_hf(i,n,q,s))*psi(i)
                   end do
#endif
                else
                   do i=n1,n2
                      trho(i)=unk_hf(i,n,q,s)*psi(i)
                   end do
                end if
           
                call Fock_fft(n1,n2,k_fock,q_fock,trho,tvht,tr)
!                call Fock_fft_parallel(n1,n2,k_fock,q_fock,trho,tvht,tr)
          
                if ( tr == 0 ) then
                   do i=n1,n2
                      tpsi(i)=tpsi(i)-c*tvht(i)*unk_hf(i,n,q,s)
                   end do
                else
#ifndef _DRSDFT_
                   do i=n1,n2
                      tpsi(i)=tpsi(i)-c*tvht(i)*conjg(unk_hf(i,n,q,s))
                   end do
#endif
                end if

             end do ! n

          end do ! q

       end do ! tr

       deallocate( trho )
       deallocate( tvht ) 

    end select
  
    if ( iflag_pbe0   > 0 ) iflag_hf=0
    if ( iflag_hse    > 0 ) iflag_hf=0
    if ( iflag_lcwpbe > 0 ) iflag_hf=0
 
    return

  END SUBROUTINE Fock

END MODULE fock_module
