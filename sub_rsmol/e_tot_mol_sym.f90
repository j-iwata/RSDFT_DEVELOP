!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE E_tot_mol_sym(iflag)
      use global_variables
      implicit none
      integer,intent(IN) :: iflag
      integer :: n,k,s,n1,n2,ns,ne,ierr,ir,nd,id,ie,ml0
      real(8) :: Enl,Eloc,sum0,EXA,EXA0,e,Eeig,Eloc_ion,Etot2
      real(8) :: ctime0,ctime1,etime0,etime1,mem,memax,fac,nop
      real(8) :: s0(9),s1(9)
      real(8),allocatable :: esp0(:,:,:),work(:),work1(:),sb(:,:)
      logical :: flag_alloc

      call watch(ctime0,etime0)

      n1 = idisp(myrank) + 1
      n2 = idisp(myrank) + ircnt(myrank)
      ML0=n2-n1+1

      Etot     = 0.d0
      Ekin     = 0.d0
      Enl      = 0.d0
      Eloc_ion = 0.d0
      Eloc     = 0.d0
      Eeig     = 0.d0
      s0(:)    = 0.d0
      s1(:)    = 0.d0
      flag_alloc = .false.
      mem=0.d0
      memax=0.d0

!
! --- KS eigenvalues, Kinetic & Nonlocal ---
!
      if ( iflag>0 ) then

         allocate( esp0(MB,MBZ,nspin) ) ; mem=mem+bdreal*MB*MBZ*nspin
         nd=maxval( irdim(1:Nir) )
         allocate( sb(4,nd) )
         mem=mem+bdreal*size(sb) ; memax=max(mem,memax)
         allocate( vtmp2(n1:n2,nd) )
         mem=mem+bdreal*size(vtmp2) ; memax=max(mem,memax)

         esp0(:,:,:) = 0.d0
         esp(:,:,:)  = 0.d0
         tsp         = 0.d0
         elocsp      = 0.d0
         enlsp       = 0.d0

         ctime_hpsi=0.d0 ; ct_hpsi=0.d0
         etime_hpsi=0.d0 ; et_hpsi=0.d0
         nop_hpsi(:)=0.d0

         select case(SYStype)
         case default

            do s=MSP_0,MSP_1
            do k=MBZ_0,MBZ_1
            do ns=MB_0,MB_1,MB_d
               ne=min(ns+MB_d-1,MB_1)
               call hpsi_spe(k,s,ns,ne,esp0(ns,k,s))
            end do
            end do
            end do

         case(1,2)

            if ( .not.allocated(LL) ) then
               allocate( LL(3,n1:n2) ) ; LL=0
               call Make_GridMap_1(LL,n1,n2)
               flag_alloc=.true.
            end if
            do s=MSP_0,MSP_1
            do k=MBZ_0,MBZ_1

               if ( isymmetry==1 ) then
                  do ir=Nir_0,Nir_1
                     N_excited=iedim(ir)
                     nd=irdim(ir)
                     fac=dV*nsym/dble(nd)
                     do ie=1,N_excited
                        do id=1,nd
                           n=irlabel2n(ir,id,ie)
                           vtmp2(n1:n2,id)=unk(n1:n2,n,k,s)
                        end do
                        call hpsi_spe_mol_sym(k,s,vtmp2,n1,n2,nd,ir,sb)
                        do id=1,nd
                           n=irlabel2n(ir,id,ie)
                           tsp=tsp+occ(n,k,s)*sb(1,id)
                           elocsp=elocsp+occ(n,k,s)*sb(2,id)
                           enlsp=enlsp+occ(n,k,s)*sb(3,id)
                           esp0(n,k,s)=sb(4,id)
                        end do
                     end do ! ie
                  end do ! ir
               else
                  do ns=MB_0,MB_1,MB_d
                     ne=min(ns+MB_d-1,MB_1)
                     call hpsi_spe_mol(k,s,ns,ne,esp0(ns,k,s),10)
                  end do
               end if
            end do ! k
            end do ! s
            if ( flag_alloc ) deallocate(LL)

         end select

         call mpi_allreduce(esp0,esp,MB*MBZ*nspin,mpi_real8,mpi_sum,mpi_comm_world,ierr)

!         if(DISP_SWITCH)then
!            do s=1,nspin
!               do k=1,MBZ
!                  do n=1,MB
!                     write(*,*) n,k,s,esp(n,k,s)
!                  end do
!               end do
!            end do
!         end if

         mem=mem-bdmain*size(vtmp2) ; deallocate( vtmp2 )
         mem=mem-bdreal*size(sb) ; deallocate( sb )
         mem=mem-bdreal*size(esp0) ; deallocate( esp0 )

         s0(1)=tsp
         s0(2)=enlsp
         call mpi_allreduce(s0,s1,2,mpi_real8,mpi_sum,mpi_comm_world,ierr)

      end if

      Eeig = sum( occ(1:MB,1:MBZ,1:nspin)*esp(1:MB,1:MBZ,1:nspin) )

!
! --- Ion local, Hartree, Exchange-Correlation ---
!

      if ( SYStype/=0 .and. isymmetry==1 ) then
         do s=MSP_0,MSP_1
            s1(3) = s1(3) + sum(rho(n1:n2,s)*Vion(n1:n2)*weight(n1:n2))*dV
            s1(4) = s1(4) + sum(rho(n1:n2,s)*Vloc(n1:n2,s)*weight(n1:n2))*dV
            s1(5) = s1(5) + 0.5d0*sum(rho(n1:n2,s)*Vh(n1:n2)*weight(n1:n2))*dV
         end do
      else
         do s=MSP_0,MSP_1
            s1(3) = s1(3) + sum(rho(n1:n2,s)*Vion(n1:n2))*dV
            s1(4) = s1(4) + sum(rho(n1:n2,s)*Vloc(n1:n2,s))*dV
            s1(5) = s1(5) + 0.5d0*sum(rho(n1:n2,s)*Vh(n1:n2))*dV
         end do
      end if
      call mpi_allreduce(s1(3),s0(3),7,mpi_real8,mpi_sum,comm_grid,ierr)
      call mpi_allreduce(s0(3),s1(3),3,mpi_real8,mpi_sum,comm_spin,ierr)
      Ekin      = s1(1)
      Enl       = s1(2)
      Eloc_ion  = s1(3)
      Eloc      = s1(4)
      E_hartree = s1(5)

      if ( ifac_vloc>0 ) then
         allocate( work(n1:n2),work1(n1:n2) )
         work(n1:n2)=rho(n1:n2,1)
         do s=2,nspin
            work(n1:n2)=work(n1:n2)+rho(n1:n2,s)
         end do
         work1(n1:n2)=Vh(n1:n2)+Vion0(n1:n2)-Vion(n1:n2)
         call Hartree_mol(n1,n2,work,work1,0.d0,2000,2*ifac_vloc)
         deallocate( work1,work )
      end if

      Etot  = Eeig - Eloc + Eloc_ion + E_hartree + E_exchange + E_correlation + Ewld
      Etot2 = Ekin + E_hartree + E_exchange + E_correlation + Eloc_ion + Enl + Ewld

      call watch(ctime1,etime1)
      if ( DISP_SWITCH ) then
         write(*,*) "TIME(E_TOT_MOL_SYM)=",ctime1-ctime0,etime1-etime0
         write(*,*) '(II) ',Ewld ,    '(CNT)',const_g0
         write(*,*) '(NL) ',Enl,      '(LOC)',Eloc_ion
         write(*,*) '(KI) ',Ekin,     '(SP) ',Eeig
         write(*,*) '(EX) ',E_exchange,       '(CO) ',E_correlation
         write(*,*) '(EXC)',EXC,      '(HT) ',E_hartree
         write(*,*) '(TO) ',Etot,     '(TO2)',Etot2
      end if

      return
      END SUBROUTINE E_tot_mol_sym
