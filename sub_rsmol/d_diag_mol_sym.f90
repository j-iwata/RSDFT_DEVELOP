!--------1---------2---------3---------4---------5---------6---------7---------8---------9

      SUBROUTINE diag_mol_sym
      use global_variables
      implicit none

      integer,allocatable :: ids(:),irc(:)
      integer :: i,k,m,ms,me,mm,ML0,n,ns,ne,nn,nme,n1,n2,mrnk,ierr
      integer :: nns,nne,nnn,mms,mme,ib1,ib2,MBLKH,s,nd,Nexmax
      complex(8) :: ztmp
      complex(8),allocatable :: work(:)
      real(8) :: ctime0,ctime1,etime0,etime1
      real(8) :: ct(7),et(7),ct0,et0,ct00,et00,ct1,et1,ct2,et2,ct3,et3
      real(8) :: nop(7),s0(30),s1(30,3),nop1,memax,mem
      integer :: j,info,loc(1),ii,WORK1,WORK2
      integer,save :: LWORK=0,LIWORK,LRWORK
      integer,allocatable :: iwork(:)
      real(8),allocatable :: rwork(:)
      character(1) :: check

      integer :: ie,ir,id,jd,isym,lt(3),je,kd,matmp
      real(8) :: c,d,sum0,fac
      real(8),allocatable :: w(:,:,:,:),u(:,:,:,:)
      real(8),allocatable :: Atmp(:,:)

      if ( MB==1 ) return

      if ( .not.(SYStype/=0 .or. isymmetry==1) ) then
         write(*,*) "this is only for mol-sym calc."
         call stop_program1("d_diag_mol_sym",1)
      end if

      call watch(ctime0,etime0)
      ct(:)=0.d0 ; et(:)=0.d0 ; ctime_hpsi=0.d0 ; etime_hpsi=0.d0

      nop(:) = 0.d0
      nop1   = 0.d0

      n1    = idisp(myrank)+1
      n2    = idisp(myrank)+ircnt(myrank)
      ML0   = ircnt(myrank)
      nme   = (MB*MB+MB)/2
      memax = 0.d0
      mem   = 0.d0

!- allocate ----------------------------------------------------
!      allocate( utmp2(MB,MB) ) ; utmp2=zero
!      mem=mem+bdmain*MB*MB ; memax=max(memax,mem)
!---------------------------------------------------------------

!- allocate ----------------------------------------------------
!      allocate( ids(0:np_band-1),irc(0:np_band-1) ) ; ids=0 ; irc=0
!      mem=mem+bsintg*np_band*2 ; memax=max(memax,mem)
!---------------------------------------------------------------

!      ids(0:np_band-1) = id_band(0:np_band-1)*ML0
!      irc(0:np_band-1) = ir_band(0:np_band-1)*ML0
!      mrnk            = id_class(myrank,4)

      do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1

!         call mpi_allgatherv(unk(n1,MB_0,k,s),irc(mrnk),TYPE_MAIN &
!                            ,unk(n1,1,k,s),irc,ids,TYPE_MAIN,comm_band,ierr)

!
! --- matrix elements ---
!

!- allocate ------------------------------------------------------------------
!         Nexmax=maxval( iedim(Nir_0:Nir_1) )
!         allocate( utmp2(Nexmax,Nexmax) )
!         mem=mem+bdmain*size(utmp2) ; memax=max(mem,memax)
!         matmp=(Nexmax*(Nexmax+1))/2
!         allocate( Atmp(matmp,2) )
!         mem=mem+bdreal*size(Atmp) ; memax=max(mem,memax)
!         n=maxval( irdim(Nir_0:Nir_1)*iedim(Nir_0:Nir_1) )
!         allocate( psi_tmp(n1:n2,n) )
!         mem=mem+bdmain*size(psi_tmp) ; memax=max(memax,mem)
!         id=maxval(irdim(Nir_0:Nir_1)
!         allocate( vtmp2(n1:n2,id) ) ; vtmp2=zero
!         allocate( wtmp2(n1:n2,id) ) ; wtmp2=zero
!         mem=mem+bdmain*size(vtmp2)+bdmain*size(wtmp2) ; memax=max(mem,memax)
!-----------------------------------------------------------------------------

         do ir=Nir_0,Nir_1

            N_excited=iedim(ir)
            nd=irdim(ir)
!            fac=dV*nsym/dble(nd)
            fac=dV/dble(nd)
            matmp=(N_excited*(N_excited+1))/2

            allocate( psi_tmp(n1:n2,N_excited*nd) ) ; mem=mem+bdreal*size(psi_tmp)
            allocate( vtmp2(n1:n2,nd) ) ; mem=mem+bdreal*size(vtmp2)
            allocate( wtmp2(n1:n2,nd) ) ; mem=mem+bdreal*size(wtmp2)
            memax=max(mem,memax)

            do ie=1,N_excited
               do id=1,nd
                  n=irlabel2n(ir,id,ie)
                  vtmp2(n1:n2,id)=unk(n1:n2,n,k,s)
               end do
               call hpsi_mol_sym(k,s,vtmp2,wtmp2,n1,n2,nd,ir)
               do id=1,nd
                  n=(ie-1)*nd+id
!                  psi_tmp(n1:n2,n)=wtmp2(n1:n2,id)
!                  psi_tmp(n1:n2,n)=wtmp2(n1:n2,id)*wkgrp(n1:n2)
                  psi_tmp(n1:n2,n)=wtmp2(n1:n2,id)*Nstar(n1:n2)
               end do
            end do

            mem=mem-bdreal*size(wtmp2) ; deallocate( wtmp2 )
            mem=mem-bdreal*size(vtmp2) ; deallocate( vtmp2 )

            allocate( Atmp(matmp,2) ) ; Atmp=0.d0
            mem=mem+bdreal*size(Atmp) ; memax=max(mem,memax)

            i=0
            do je=1,N_excited
            do ie=1,je
               sum0=0.d0
               do id=1,nd
                  n=(je-1)*nd+id
                  m=irlabel2n(ir,id,ie)
!                  sum0=sum0+sum(unk(n1:n2,m,k,s)*psi_tmp(n1:n2,n)*wkgrp(n1:n2))
                  sum0=sum0+sum(unk(n1:n2,m,k,s)*psi_tmp(n1:n2,n))
               end do
               i=i+1
               Atmp(i,1)=sum0*fac
            end do
            end do

            mem=mem-bdreal*size(psi_tmp) ; deallocate( psi_tmp )

            call mpi_allreduce(Atmp(1,1),Atmp(1,2),matmp,mpi_real8,mpi_sum,comm_grid,ierr)

            allocate( utmp2(N_excited,N_excited) ) ; utmp2=0.d0
            mem=mem+bdreal*size(utmp2) ; memax=max(mem,memax)

            i=0
            do je=1,N_excited
            do ie=1,je
               i=i+1
               utmp2(ie,je)=Atmp(i,2)
            end do
            end do

            mem=mem-bdreal*size(Atmp) ; deallocate( Atmp )
!
!- diag -
!
            LWORK =max( LWORK, 1+6*N_excited+2*N_excited**2 )
            LIWORK=max( LIWORK,3+5*N_excited )
            
            allocate( rwork(LWORK),iwork(LIWORK) )
            mem=mem+bdreal*LWORK+bsintg*LIWORK ; memax=max(memax,mem)

            ns=irlabel2n(ir,1,1)
            ne=irlabel2n(ir,1,N_excited)
            nn=ne-ns+1

            call DSYEVD('V','U',N_excited,utmp2,N_excited,esp(ns,k,s),rwork,LWORK,iwork,LIWORK,info)

            LWORK=max( LWORK,nint(rwork(1)) )
            LIWORK=max( LIWORK,iwork(1) )

            mem=mem-bsintg*LIWORK-bdreal*LWORK ; deallocate( iwork,rwork )

            do jd=2,nd
               ms=irlabel2n(ir,jd,1)
               me=irlabel2n(ir,jd,N_excited)
               esp(ms:me,k,s)=esp(ns:ne,k,s)
            end do

            allocate( psi_tmp(n1:n2,N_excited) )
            mem=mem+bdreal*size(psi_tmp) ; memax=max(mem,memax)

            do jd=1,irdim(ir)
               ms=irlabel2n(ir,jd,1)
               me=irlabel2n(ir,jd,N_excited)
               psi_tmp(:,:)=0.d0
               call dgemm(TRANSB,TRANSB,ML0,nn,nn,one,unk(n1,ms,k,s),ML0, &
                    utmp2,N_excited,zero,psi_tmp(n1,1),ML0)
               do n=ms,me
                  unk(n1:n2,n,k,s)=psi_tmp(n1:n2,n-ms+1)
               end do
            end do ! jd

            mem=mem-bdreal*size(psi_tmp) ; deallocate(psi_tmp)
            mem=mem-bdreal*size(utmp2) ; deallocate( utmp2 )

         end do ! ir

      end do ! k
      end do ! s

!      mem=mem-bsintg*np_band*2 ; deallocate( irc,ids )

      call watch(ctime1,etime1)

      s0(29)=ctime1-ctime0
      s0(30)=etime1-etime0
      call mpi_allreduce(s0,s1(1,1),30,mpi_real8,mpi_min,mpi_comm_world,ierr)
      call mpi_allreduce(s0,s1(1,2),30,mpi_real8,mpi_max,mpi_comm_world,ierr)

      if (DISP_SWITCH) then
         write(*,'(1x,"TIME(D_DIAG_MOL_SYM)=",4f10.3)') s1(29,1:2),s1(30,1:2)
!         write(*,'(1x," (MATE)      ",6f10.3)') s1(12,1:2),s1(19,1:2),s1(5,1:2)/1.d6
!         write(*,'(1x,"  mate       ",6f10.3)') s1(8,1:2),s1(15,1:2),s1(1,1:2)/1.d6
!         write(*,'(1x,"  hpsi       ",6f10.3)') s1(9,1:2),s1(16,1:2),s1(2,1:2)/1.d6
!         write(*,'(1x,"  allgatherv ",4f10.3)') s1(10,1:2),s1(17,1:2)
!         write(*,'(1x,"  allreduce  ",4f10.3)') s1(11,1:2),s1(18,1:2)
!         write(*,'(1x," (",a6,")    ",4f10.3)') idiag0,s1(13,1:2),s1(20,1:2)
!         write(*,'(1x," (ROTV)      ",6f10.3)') s1(14,1:2),s1(21,1:2),s1(7,1:2)/1.d6
!         write(*,*) "LWORK,WORK1=",LWORK,WORK1
         write(*,*) "LWORK=",LWORK
         write(*,*) "MEM(MB) =",memax*B2MB,mem
!         write(*,*) "check(nop1-sum(nop(1))=",nop1-s1(22,3)
      end if

!      LWORK=max(LWORK,WORK1)
!      LIWORK=max(LIWORK,WORK2)

      return

      END SUBROUTINE diag_mol_sym
