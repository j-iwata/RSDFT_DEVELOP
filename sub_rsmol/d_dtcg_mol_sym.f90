!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE DTcg_mol_sym(mcg1,mcg2,igs)
      use global_variables
      implicit none
      integer,intent(IN) :: mcg1,mcg2,igs
      integer :: s,k,ns,ne,nn,Mcg,n,m,icg,n1,n2,ML0,Nhpsi,Npc,Ncgtot,ierr
      integer :: mm,icmp,id,ir,nd,ie,ndmax
      real(8),parameter :: ep0=0.d0
      real(8),parameter :: ep1=1.d-15
      real(8) :: rwork(9),W(2),c,d,r,c1,fac
      real(8) :: mem,memax,nop(6),ct(6),et(6),s0(26),s1(26,3)
      real(8) :: ctime0,ctime1,etime0,etime1
      real(8),allocatable :: sb(:),rb(:),E(:),E1(:),gkgk(:),bk(:)
      real(8),allocatable :: esp0(:,:,:),xk(:,:)
      complex(8) :: work(9),zphase,ztmp

      if ( .not.(SYStype/=0 .and. isymmetry==1) ) goto 900

      call watch(ctime0,etime0)

      n1  = idisp(myrank)+1
      n2  = idisp(myrank)+ircnt(myrank)
      ML0 = ircnt(myrank)
      mm  = ML0  ; if (TYPE_MAIN==mpi_complex16) mm=2*ML0
      c1  = 2.d0 ; if (TYPE_MAIN==mpi_complex16) c1=1.d0
      icmp= 1    ; if (TYPE_MAIN==mpi_complex16) icmp=2

      Ncgtot = 0
      Nhpsi  = 0
      Npc    = 0
      mem    = 0.d0
      memax  = 0.d0
      mem_tmp_sub = 0.d0

      nop_pc=0.d0
      nop(:)=0.d0
      ct(:)=0.d0
      et(:)=0.d0

      ctime_hpsi=0.d0
      etime_hpsi=0.d0

      ndmax=maxval(irdim(Nir_0:Nir_1))

      allocate( hxk(n1:n2,ndmax), hpk(n1:n2,ndmax) )
      allocate( gk(n1:n2,ndmax) , Pgk(n1:n2,ndmax) )
      allocate( pk(n1:n2,ndmax) , pko(n1:n2,ndmax) )
      mem=mem+bdmain*ML0*ndmax*6
      allocate( sb(ndmax),rb(ndmax) )
      allocate( E(ndmax),E1(ndmax),gkgk(ndmax),bk(ndmax) )
      mem=mem+bdreal*ndmax*4+bdreal*ndmax*2
      allocate( vtmp2(ndmax,6),wtmp2(ndmax,6) )
      allocate( utmp2(2,2),btmp2(2,2) )
      mem=mem+bdmain*ndmax*6*2+bdmain*4*2
      allocate( xk(n1:n2,ndmax) ) ; xk=0.d0
      mem=mem+bdreal*size(xk)

      memax = max(memax,mem)


      Mcg = mcg1

!
! ---
!
      do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1

         do ir=Nir_0,Nir_1

            nd=irdim(ir)
!            fac=dV*nsym/dble(nd)
            fac=dV/dble(nd)
            N_excited=iedim(ir)

         do ie=1,N_excited

            do id=1,nd
               n=irlabel2n(ir,id,ie)
               xk(n1:n2,id)=unk(n1:n2,n,k,s)
            end do

            E1(:)=1.d10

            call hpsi_mol_sym(k,s,xk,hxk,n1,n2,nd,ir) ; Nhpsi=Nhpsi+1

            call watch1(ct(1),et(1),0)

            call dot_product_mol_sym(xk,hxk,sb,Nstar,fac,nd,mm,nop(1))

            call watch1(ct(1),et(1),1)

            call mpi_allreduce(sb,E,nd,mpi_real8,mpi_sum,comm_grid,ierr)

            call watch1(ct(3),et(3),1)

            do id=1,nd
               gk(n1:n2,id)=-c1*(hxk(n1:n2,id)-E(id)*xk(n1:n2,id))
            end do
            nop(1)=nop(1)+(sop3*ML0+sop3*ML0+sop3*ML0)*nd

            call dot_product_mol_sym(gk,gk,sb,Nstar,fac,nd,mm,nop(1))
!            call dot_product_mol_sym(gk,gk,sb,Nstar,1.d0,nd,mm,nop(1))

            call watch1(ct(1),et(1),1)

            call mpi_allreduce(sb,rb,nd,mpi_real8,mpi_sum,comm_grid,ierr)

            call watch1(ct(3),et(3),1)

            do icg=1,Mcg+1

               Ncgtot=Ncgtot+1

               do id=1,nd
                  Pgk(n1:n2,id)=gk(n1:n2,id)
               end do

               do id=1,nd
                  n=irlabel2n(ir,id,ie)
                  res(n,k,s)=rb(id)/c1**2
!                  res(n,k,s)=rb(id)*fac/c1**2
               end do

!
! --- Convergence check ---
!

!               do id=1,nd
!                  write(*,'(1x,4i5,f20.14,g20.8)') icg,ir,ie,id,E(id),rb(id)
!               end do
!               if(DISP_SWITCH) write(*,*) "ir,ie,icg=",ir,ie,icg
!               call chk_dtcg(s,k,n1,n2,ir,ie,nd,xk)

               if ( all(rb(1:nd)<ep0) ) exit
!               if ( all(rb(1:nd)*fac<ep0) ) exit
               if ( all(abs(E(1:nd)-E1(1:nd))<ep1) ) exit
               if ( icg==Mcg+1 ) exit

!
! --- Preconditioning ---
!

               call watch1(ct(4),et(4),0)

               call precond_cg_mol_sym(E,k,s,nd,ML0,ir,Npc)

               call watch1(ct(4),et(4),1)

!
! --- Orthogonalization ---
!

               select case(igs)
               case(1) ! for all orbitals

               case(2) ! only for unoccupied orbitals

               case(3:)

               case(0)

               end select

               call watch1(ct(5),et(5),1)

!
! ---
!

               call dot_product_mol_sym(Pgk,gk,sb,Nstar,fac,nd,mm,nop(1))
!               call dot_product_mol_sym(Pgk,gk,sb,Nstar,1.d0,nd,mm,nop(1))

               call watch1(ct(1),et(1),1)

               call mpi_allreduce(sb,rb,nd,mpi_real8,mpi_sum,comm_grid,ierr)

               call watch1(ct(3),et(3),1)

               if ( icg==1 ) then
                  pk(n1:n2,1:nd) = Pgk(n1:n2,1:nd)
               else
                  do id=1,nd
                     bk(id)=rb(id)/gkgk(id)
                     pk(n1:n2,id)=Pgk(n1:n2,id)+bk(id)*pk(n1:n2,id)
                  end do
               end if
               gkgk(1:nd)=rb(1:nd)

               call watch1(ct(1),et(1),1)

               call hpsi_mol_sym(k,s,pk,hpk,n1,n2,nd,ir) ; Nhpsi=Nhpsi+1

               call watch1(ct(2),et(2),1)

               vtmp2(:,:)=zero
!               call dot_product_mol_sym(xk,xk ,vtmp2(1,1),Nstar,1.d0,nd,mm,nop(1))
!               call dot_product_mol_sym(pk,xk ,vtmp2(1,2),Nstar,1.d0,nd,mm,nop(1))
!               call dot_product_mol_sym(pk,pk ,vtmp2(1,3),Nstar,1.d0,nd,mm,nop(1))
!               call dot_product_mol_sym(xk,hxk,vtmp2(1,4),Nstar,1.d0,nd,mm,nop(1))
!               call dot_product_mol_sym(pk,hxk,vtmp2(1,5),Nstar,1.d0,nd,mm,nop(1))
!               call dot_product_mol_sym(pk,hpk,vtmp2(1,6),Nstar,1.d0,nd,mm,nop(1))
               call dot_product_mol_sym(xk,xk ,vtmp2(1,1),Nstar,fac,nd,mm,nop(1))
               call dot_product_mol_sym(pk,xk ,vtmp2(1,2),Nstar,fac,nd,mm,nop(1))
               call dot_product_mol_sym(pk,pk ,vtmp2(1,3),Nstar,fac,nd,mm,nop(1))
               call dot_product_mol_sym(xk,hxk,vtmp2(1,4),Nstar,fac,nd,mm,nop(1))
               call dot_product_mol_sym(pk,hxk,vtmp2(1,5),Nstar,fac,nd,mm,nop(1))
               call dot_product_mol_sym(pk,hpk,vtmp2(1,6),Nstar,fac,nd,mm,nop(1))

               call watch1(ct(1),et(1),1)

               call mpi_allreduce(vtmp2,wtmp2,6*ndmax,TYPE_MAIN,mpi_sum,comm_grid,ierr)

               call watch1(ct(3),et(3),1)

               do id=1,nd

                  btmp2(1,1)=wtmp2(id,1)
                  btmp2(2,1)=wtmp2(id,2)
                  btmp2(1,2)=wtmp2(id,2)
                  btmp2(2,2)=wtmp2(id,3)
                  utmp2(1,1)=wtmp2(id,4)
                  utmp2(2,1)=wtmp2(id,5)
                  utmp2(1,2)=wtmp2(id,5)
                  utmp2(2,2)=wtmp2(id,6)
                  call dsygv(1,'V','U',2,utmp2,2,btmp2,2,W,rwork,9,ierr)
                  if ( abs(W(1)-E(id))>1.d-1 .and. abs(W(2)-E(id))<=1.d-1 ) then
                     utmp2(1,1)=utmp2(1,2)
                     utmp2(2,1)=utmp2(2,2)
                     W(1)=W(2)
                  end if
!- Fix the phase -
                  c=utmp2(1,1)
                  if( c<0.d0 ) then
                     utmp2(1,1)=-utmp2(1,1)
                     utmp2(2,1)=-utmp2(2,1)
                  end if

                  E1(id)=E(id)
                  E(id) =W(1)

                  hxk(n1:n2,id)=utmp2(1,1)*hxk(n1:n2,id)+utmp2(2,1)*hpk(n1:n2,id)

                  gk(n1:n2,id)=-c1*( hxk(n1:n2,id)-W(1)*( utmp2(1,1)*xk(n1:n2,id) &
                                                         +utmp2(2,1)*pk(n1:n2,id) ) )

                  nop(1)=nop(1)+sop4*ML0+sop1*ML0
                  nop(1)=nop(1)+2.d0*sop3+sop1*ML0+sop1*ML0 + ML0
                  nop(1)=nop(1)+sop2*ML0+1.d0
               end do

               call dot_product_mol_sym(gk,gk,sb,Nstar,fac,nd,mm,nop(1))
!               call dot_product_mol_sym(gk,gk,sb,Nstar,1.d0,nd,mm,nop(1))

               call watch1(ct(1),et(1),1)

               call mpi_allreduce(sb,rb,nd,mpi_real8,mpi_sum,comm_grid,ierr)

               call watch1(ct(3),et(3),1)

               do id=1,nd
                  m=irlabel2n(ir,id,ie)
                  if ( rb(id)/res(m,k,s)>1.d8 ) then
!                  if ( rb(id)*fac/res(m,k,s)>1.d8 ) then
                     E(id)=E1(id)
                     cycle
                  end if
                  xk(n1:n2,id)=utmp2(1,1)*xk(n1:n2,id)+utmp2(2,1)*pk(n1:n2,id)
                  nop(1)=nop(1)+sop4*ML0+sop1*ML0
               end do

               call watch1(ct(1),et(1),1)

            end do ! icg

            do id=1,nd
               n=irlabel2n(ir,id,ie)
               esp(n,k,s)=E(id)
               unk(n1:n2,n,k,s)=xk(n1:n2,id)
            end do

         end do ! ie
         end do ! ir

      end do ! k-loop
      end do ! s-loop

      mem=mem-bdreal*size(xk) ; deallocate( xk )
      mem=mem-bdmain*size(btmp2)-bdmain*size(utmp2) ; deallocate( btmp2,utmp2  )
      mem=mem-bdmain*size(wtmp2)-bdmain*size(vtmp2) ; deallocate( wtmp2,vtmp2  )
      mem=mem-bdreal*(size(bk)+size(gkgk)+size(E1)+size(E)) ; deallocate( bk,gkgk,E1,E )
      mem=mem-bdreal*(size(rb)+size(sb)) ; deallocate( rb,sb )
      mem=mem-bdmain*size(pko)-bdmain*size(pk)  ; deallocate( pko,pk  )
      mem=mem-bdmain*size(Pgk)-bdmain*size(gk)  ; deallocate( Pgk,gk  )
      mem=mem-bdmain*size(hpk)-bdmain*size(hxk) ; deallocate( hpk,hxk )

      call watch(ctime1,etime1)

      nop(2)=Nhpsi*sum(nop_hpsi(1:7))/MB
      nop(4)=nop_pc
      nop(6)=sum(nop(:))

      ct(2)=ctime_hpsi
      et(2)=etime_hpsi
      ct(6)=ctime1-ctime0
      et(6)=etime1-etime0

      s0(1:6)=ct(1:6)
      s0(7:12)=et(1:6)
      s0(13:18)=nop(1:6)
      where(ct(1:6)>0.d0)
         s0(19:24)=nop(1:6)/ct(1:6)
      end where

      call mpi_allreduce(s0,s1(1,1),24,mpi_real8,mpi_min,mpi_comm_world,ierr)
      call mpi_allreduce(s0,s1(1,2),24,mpi_real8,mpi_max,mpi_comm_world,ierr)
      call mpi_allreduce(s0,s1(1,3),24,mpi_real8,mpi_sum,mpi_comm_world,ierr)

      if (DISP_SWITCH) then
         write(*,'(1x,"TIME(D_DTCG_MOL_SYM)=",6f10.3)') s1(6,1:2),s1(12,1:2),s1(24,1:2)/1.d6
         write(*,'(1x,"  hpsi       ",6f10.3)') s1(2,1:2),s1(8,1:2),s1(20,1:2)/1.d6
         write(*,'(1x,"  pc         ",6f10.3)') s1(4,1:2),s1(10,1:2),s1(22,1:2)/1.d6
         write(*,'(1x,"  gs         ",6f10.3)') s1(5,1:2),s1(11,1:2),s1(23,1:2)/1.d6
         write(*,'(1x,"  op         ",6f10.3)') s1(1,1:2),s1(7,1:2),s1(19,1:2)/1.d6
         write(*,'(1x,"  allreduce  ",4f10.3)') s1(3,1:2),s1(9,1:2)
         write(*,*) " (Ncg,Nhpsi,Npc) ",Ncgtot,Nhpsi,Npc
         write(*,*) " (mem)(MB):",memax*B2MB,mem
      end if

      return

 900  call stop_program

      END SUBROUTINE DTcg_mol_sym


      SUBROUTINE dot_product_mol_sym(a,b,c,w,alpha,n,m,nop)
      implicit none
      integer,intent(IN) :: n,m
      real(8),intent(IN) :: a(m,n),b(m,n),alpha
      integer,intent(IN) :: w(m)
      real(8),intent(INOUT) :: c(n),nop
      real(8) :: sum0
      integer :: i,j

      sum0=0.d0
      do j=1,n
         do i=1,m
            sum0=sum0+a(i,j)*b(i,j)*w(i)
         end do
      end do
      sum0=sum0*alpha
      c(1:n)=sum0

      nop=nop+3.d0*m*n+1.d0

      return
      END SUBROUTINE dot_product_mol_sym



      SUBROUTINE chk_dtcg(s,k,n1,n2,ir0,ie0,nd0,xk)
      use global_variables
      implicit none
      integer,intent(IN) :: s,k,n1,n2,ir0,ie0,nd0
      real(8),intent(IN) :: xk(n1:n2,nd0)
      integer :: n,ir,id,ie,id0,i,isym,lt(3),m,i1,i2,i3,ierr,nd,jr,jd,je
      integer,allocatable :: irc(:),ids(:)
      real(8),allocatable :: w(:,:,:,:),tmp(:,:),tmp0(:,:),sum0(:),sum1(:)
      real(8) :: c,ctime0,ctime1,etime0,etime1

!      call bwatch(ctime0,etime0)

      allocate( w(-ML1:ML1,-ML2:ML2,-ML3:ML3,MB) ) ; w=0.d0

      do n=MB_0,MB_1
         ir=irlabel(1,n)
         id=irlabel(2,n)
         ie=irlabel(3,n)
         nd=irdim(ir)
         do jd=1,nd
            m=irlabel2n(ir,jd,ie)
            do isym=1,nsym
               c=Rir(id,jd,isym,ir)
               do i=n1,n2
                  i1=rga(1,1,isym)*LL(1,i)+rga(1,2,isym)*LL(2,i)+rga(1,3,isym)*LL(3,i)
                  i2=rga(2,1,isym)*LL(1,i)+rga(2,2,isym)*LL(2,i)+rga(2,3,isym)*LL(3,i)
                  i3=rga(3,1,isym)*LL(1,i)+rga(3,2,isym)*LL(2,i)+rga(3,3,isym)*LL(3,i)
                  w(i1,i2,i3,n)=w(i1,i2,i3,n)+c*unk(i,m,k,s)*Nstar(i)
               end do
            end do
         end do
         c=1.d0/nsym
         w(:,:,:,n)=c*w(:,:,:,n)
      end do

      if ( ir0/=0 ) then
      do id0=1,nd0
         n=irlabel2n(ir0,id0,ie0)
         w(:,:,:,n)=zero
         do id=1,nd0
            do isym=1,nsym
               c=Rir(id0,id,isym,ir0)
               do i=n1,n2
                  i1=rga(1,1,isym)*LL(1,i)+rga(1,2,isym)*LL(2,i)+rga(1,3,isym)*LL(3,i)
                  i2=rga(2,1,isym)*LL(1,i)+rga(2,2,isym)*LL(2,i)+rga(2,3,isym)*LL(3,i)
                  i3=rga(3,1,isym)*LL(1,i)+rga(3,2,isym)*LL(2,i)+rga(3,3,isym)*LL(3,i)
                  w(i1,i2,i3,n)=w(i1,i2,i3,n)+c*xk(i,id)*Nstar(i)
               end do
            end do
         end do
         c=1.d0/nsym
         w(:,:,:,n)=c*w(:,:,:,n)
      end do
      end if

      allocate( irc(0:np_band-1),ids(0:np_band-1) )
      irc(0:np_band-1)=size(w(:,:,:,1))*ir_band(0:np_band-1)
      do n=0,np_band-1
         ids(n)=sum(irc(0:n))-irc(n)
      end do
      call mpi_allgatherv(w(:,:,:,MB_0),irc(myrank_b),mpi_real8,w,irc,ids &
           ,mpi_real8,comm_band,ierr)
      deallocate( ids,irc )

      allocate( tmp(MB,MB),tmp0(MB,MB) )
      tmp0=0.d0 ; tmp=0.d0
      do n=1,MB
         do m=1,n
            tmp0(m,n)=sum(w(:,:,:,m)*w(:,:,:,n))*dV
            tmp0(n,m)=tmp0(m,n)
         end do
      end do
      call mpi_allreduce(tmp0,tmp,MB*MB,mpi_real8,mpi_sum,comm_grid,ierr)
      deallocate( tmp0 )

      deallocate( w )

      if ( ir0/=0 ) then
      allocate( sum0(nd0),sum1(nd0) )
      do ir=1,Nir
      do id=1,irdim(ir)
         sum0(:)=0.d0
         sum1(:)=0.d0
         do ie=1,iedim(ir)
            m=irlabel2n(ir,id,ie)
            do id0=1,nd0
               n=irlabel2n(ir0,id0,ie0)
               if ( m==n ) then
                  sum1(id0)=sum1(id0)+abs(tmp(m,n))
               else
                  sum0(id0)=sum0(id0)+abs(tmp(m,n))
               end if
            end do
         end do
         do id0=1,nd0
            if (DISP_SWITCH) write(*,'(1x,4i4,2g24.14)') ir0,id0,ir,id,sum0(id0),sum1(id0)
         end do
      end do
      end do
      deallocate( sum1,sum0 )
      else
         allocate( sum0(4) )
         do jr=1,Nir
         do jd=1,irdim(jr)
         do je=1,iedim(jr)
            n=irlabel2n(jr,jd,je)
            sum0(:)=0.d0
            do ir=1,Nir
            do id=1,irdim(ir)
               do ie=1,iedim(ir)
                  m=irlabel2n(ir,id,ie)
                  if ( m==n ) then
                     sum0(1)=sum0(1)+abs(tmp(m,n))
                  else if ( ir==jr .and. id==jd ) then
                     sum0(2)=sum0(2)+abs(tmp(m,n))
                  else if ( ir==jr ) then
                     sum0(3)=sum0(3)+abs(tmp(m,n))
                  else
                     sum0(4)=sum0(4)+abs(tmp(m,n))
                  end if
               end do
            end do
            end do
            if (DISP_SWITCH) write(*,'(1x,3i4,4g20.10)') jr,jd,je,sum0(1:4)
         end do
         end do
         end do
         deallocate( sum0 )
      end if

      deallocate( tmp )

!      call bwatch(ctime1,etime1)

      return
      END SUBROUTINE chk_dtcg
