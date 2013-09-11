!--------1---------2---------3---------4---------5---------6---------7---------8---------9

      SUBROUTINE gsmol_test
      use global_variables
      implicit none
      integer :: s,k,ir,id,ns,ne,n1,n2,m,n,ierr,nd,m0
      real(8) :: mem,memax,c0,c,fac
      real(8) :: ctime0,ctime1,etime0,etime1
      real(8),allocatable :: tmp1(:),tmp2(:)
      real(8) :: d,sum0
      integer :: lt(3),j1,j2,j3,isym,j,i,i1
      integer :: ie,je,jd

      call bwatch(ctime0,etime0)

      n1=idisp(myrank)+1
      n2=idisp(myrank)+ircnt(myrank)

      mem=0.d0
      memax=0.d0

!      n=maxval( iedim(Nir_0:Nir_1) )
!      allocate( utmp(n) )
!      allocate( vtmp(n) )
!      mem=mem+bdmain*size(utmp)+bdmain*size(vtmp) ; memax=max(mem,memax)
      allocate( tmp1(n1:n2) )
      allocate( tmp2(n1:n2) )
      mem=mem+bdreal*size(tmp1)+bdreal*size(tmp2) ; memax=max(mem,memax)

!      tmp1(n1:n2)=sqrt(wkgrp(n1:n2))
      tmp1(n1:n2)=sqrt(dble(Nstar(n1:n2)))
      tmp2(n1:n2)=1.d0/tmp1(n1:n2)

      do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1

         do ir=Nir_0,Nir_1

            N_excited=iedim(ir)
            nd=irdim(ir)
!            fac=dble(nsym)*dV/dble(nd)
            fac=dV/dble(nd)

            do n=irlabel2n(ir,1,1),irlabel2n(ir,nd,N_excited)
               unk(n1:n2,n,k,s)=unk(n1:n2,n,k,s)*tmp1(n1:n2)
            end do
!
! --- Simple GS ---
!
            goto 1
            do ie=1,N_excited
               if ( ie>1 ) then
                  vtmp(1:ie-1)=zero
                  do jd=1,nd
                     n =irlabel2n(ir,jd,ie)
                     m0=irlabel2n(ir,jd,1)-1
                     do je=1,ie-1
                        vtmp(je)=vtmp(je)+sum(unk(n1:n2,m0+je,k,s)*unk(n1:n2,n,k,s))
                     end do
                  end do
                  vtmp(1:ie-1)=vtmp(1:ie-1)*fac
                  call mpi_allreduce(vtmp,utmp,ie-1,TYPE_MAIN,mpi_sum,comm_grid,ierr)
                  do id=1,nd
                     n =irlabel2n(ir,id,ie)
                     m0=irlabel2n(ir,id,1)-1
                     do je=1,ie-1
                        unk(n1:n2,n,k,s)=unk(n1:n2,n,k,s)-utmp(je)*unk(n1:n2,m0+je,k,s)
                     end do
                  end do
               end if !(ie>1)
               c0=0.d0
               do jd=1,irdim(ir)
                  m=irlabel2n(ir,jd,ie)
                  c0=c0+sum(unk(n1:n2,m,k,s)**2)
               end do
               call mpi_allreduce(c0,c,1,TYPE_MAIN,mpi_sum,comm_grid,ierr)
               c=1.d0/sqrt(c*fac)
               do id=1,irdim(ir)
                  n=irlabel2n(ir,id,ie)
                  unk(n1:n2,n,k,s)=c*unk(n1:n2,n,k,s)
               end do
            end do ! ie
 1          continue

!
!--- Takahashi's GS ---
!

            c0=0.d0
            do id=1,nd
               m=irlabel2n(ir,id,1)
               c0=c0+sum(unk(n1:n2,m,k,s)**2)
            end do
            call mpi_allreduce(c0,c,1,TYPE_MAIN,mpi_sum,comm_grid,ierr)
!            c=1.d0/sqrt(c*fac)
            c=sqrt(dble(nd))/sqrt(c*dV)
            do id=1,nd
               n=irlabel2n(ir,id,1)
               unk(n1:n2,n,k,s)=c*unk(n1:n2,n,k,s)
            end do

            call gsmol_sub(2,N_excited,1,N_excited-1,NBLK,ir,nd,k,s,fac,mem,memax)

            do n=irlabel2n(ir,1,1),irlabel2n(ir,nd,N_excited)
               unk(n1:n2,n,k,s)=unk(n1:n2,n,k,s)*tmp2(n1:n2)
            end do

         end do ! ir

      end do ! k
      end do ! s

      mem=mem-bdreal*size(tmp2) ; deallocate( tmp2 )
      mem=mem-bdreal*size(tmp1) ; deallocate( tmp1 )

!      mem=mem-bdmain*size(vtmp) ; deallocate( vtmp )
!      mem=mem-bdmain*size(utmp) ; deallocate( utmp )

      call bwatch(ctime1,etime1)

      if (DISP_SWITCH) then
         write(*,*) "TIME(D_GSMOL_TEST)",ctime1-ctime0,etime1-etime0
         write(*,*) "MEM(MB)",memax*B2MB,mem
         write(*,*) "MB,NBLK,NBLK1=",MB,NBLK,NBLK1
      end if

      return

 900  call stop_program1("gsmol_test",1)

      END SUBROUTINE gsmol_test

!--------1---------2---------3---------4---------5---------6---------7---------8---------9

      RECURSIVE SUBROUTINE gsmol_sub(b0,b1,a0,a1,MBLK,ir,nd,k,s,fac,mem,memax)
      use global_variables
      implicit none
      integer,intent(IN) :: b0,b1,a0,a1,MBLK,ir,nd,k,s
      real(8),intent(IN) :: fac
      real(8),intent(INOUT) :: mem,memax
      real(8) :: c,c0
      integer :: a,as,ae,b,bs,be,na,nb,n1,n2,ML0,ms,ns,ierr,MBLKH
      integer :: id,m,n,ie,je

      n1 =idisp(myrank)+1
      n2 =idisp(myrank)+ircnt(myrank)
      ML0=n2-n1+1

      do a=a0,a1,MBLK

         as=a
         ae=min(as+MBLK-1,a1)
         na=ae-as+1

         do b=b0,b1,MBLK

            bs=b
            be=min(bs+MBLK-1,b1)
            nb=be-bs+1

            if ( as>=be ) cycle

            if ( ae+1==be ) then

               if ( na<=NBLK1 .and. nb<=NBLK1 ) then

                  allocate( vtmp(as:ae),utmp(as:ae) )
                  mem=mem+bdmain*size(vtmp)*2 ; memax=max(mem,memax)

                  do ie=bs,be

                     vtmp(as:ie-1)=zero
                     do id=1,nd
                        n=irlabel2n(ir,id,ie)
                        do je=as,ie-1
                           m=irlabel2n(ir,id,je)
                           vtmp(je)=vtmp(je)+sum(unk(n1:n2,m,k,s)*unk(n1:n2,n,k,s))
                        end do
                     end do
                     vtmp(as:ie-1)=vtmp(as:ie-1)*fac
                     call mpi_allreduce(vtmp,utmp,ie-as,TYPE_MAIN,mpi_sum,comm_grid,ierr)
                     do id=1,nd
                        n=irlabel2n(ir,id,ie)
                        do je=as,ie-1
                           m=irlabel2n(ir,id,je)
                           unk(n1:n2,n,k,s)=unk(n1:n2,n,k,s)-utmp(je)*unk(n1:n2,m,k,s)
                        end do
                     end do

                     c0=0.d0
                     do id=1,nd
                        m=irlabel2n(ir,id,ie)
                        c0=c0+sum(unk(n1:n2,m,k,s)**2)
                     end do
                     call mpi_allreduce(c0,c,1,TYPE_MAIN,mpi_sum,comm_grid,ierr)
                     c=1.d0/sqrt(c*fac)
                     do id=1,nd
                        n=irlabel2n(ir,id,ie)
                        unk(n1:n2,n,k,s)=c*unk(n1:n2,n,k,s)
!                        write(*,'(1x,8i4)') ir,id,ie,n,bs,be,MBLK
                     end do

                  end do ! ie

                  mem=mem-bdmain*size(utmp)*2 ; deallocate( utmp,vtmp )

               else

                  MBLKH=max( (MBLK+1)/2, NBLK1 )
                  call gsmol_sub(bs,be,as,ae,MBLKH,ir,nd,k,s,fac,mem,memax)
                  
               end if

            else

               allocate( utmp2(na,nb),vtmp2(na,nb) )
               mem=mem+bdmain*size(utmp2)+bdmain*size(vtmp2) ; memax=max(mem,memax)

               utmp2(:,:)=zero
               do id=1,nd
                  ms=irlabel2n(ir,id,as)
                  ns=irlabel2n(ir,id,bs)
                  call dgemm(TRANSA,TRANSB,na,nb,ML0,-fac,unk(n1,ms,k,s),ML0,unk(n1,ns,k,s),ML0,one,utmp2,na)
               end do

               call mpi_allreduce(utmp2,vtmp2,na*nb,TYPE_MAIN,mpi_sum,comm_grid,ierr)

               do id=1,nd
                  ms=irlabel2n(ir,id,as)
                  ns=irlabel2n(ir,id,bs)
                  call dgemm(TRANSB,TRANSB,ML0,nb,na,one,unk(n1,ms,k,s),ML0,vtmp2,na,one,unk(n1,ns,k,s),ML0)
               end do

               mem=mem-bdmain*size(vtmp2)-bdmain*size(utmp2) ; deallocate( vtmp2,utmp2 )

            end if

         end do ! b

      end do ! a

      return
      END SUBROUTINE gsmol_sub
