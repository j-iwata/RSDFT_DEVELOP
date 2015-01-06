!--------1---------2---------3---------4---------5---------6---------7---------8---------9

      SUBROUTINE re_sym
      use global_variables
      implicit none
      integer :: n,k,i,n1,n2,i1,i2,i3,s,m,j,ir,id,ie,jd
      integer :: lt(3),isym,idum,ierr,msz,nd
      real(8) :: c,d,ctime0,ctime1,etime0,etime1,mem,memax
      real(8),allocatable :: w(:,:,:,:),w0(:,:,:),w1(:,:,:)

      call bwatch(ctime0,etime0)
      if (DISP_SWITCH) write(*,'(a60," re_sym")') repeat("-",60)

      n1=idisp(myrank)+1
      n2=idisp(myrank)+ircnt(myrank)
      mem=0.d0
      memax=0.d0

      allocate( w0(-ML1:ML1,-ML2:ML2,-ML3:ML3) ) ; w0=0.d0
      allocate( w1(-ML1:ML1,-ML2:ML2,-ML3:ML3) ) ; w1=0.d0
      msz=size(w0)
      mem=mem+bdreal*msz*2 ; memax=max(mem,memax)

      idum=456214531

      do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1

         do ir=Nir_0,Nir_1

            N_excited=iedim(ir)
            nd=irdim(ir)

            do ie=1,N_excited

               w0(:,:,:)=zero
               do jd=1,nd
                  m=irlabel2n(ir,jd,ie)
                  do isym=1,nsym
                     c=Rir(1,jd,isym,ir)
                     do i=n1,n2
                        i1=rga(1,1,isym)*LL(1,i)+rga(1,2,isym)*LL(2,i)+rga(1,3,isym)*LL(3,i)
                        i2=rga(2,1,isym)*LL(1,i)+rga(2,2,isym)*LL(2,i)+rga(2,3,isym)*LL(3,i)
                        i3=rga(3,1,isym)*LL(1,i)+rga(3,2,isym)*LL(2,i)+rga(3,3,isym)*LL(3,i)
                        w0(i1,i2,i3)=w0(i1,i2,i3)+c*unk(i,m,k,s)*Nstar(i)
                     end do
                  end do
               end do
               c=1.d0/nsym
               w0(:,:,:)=c*w0(:,:,:)
               call mpi_allreduce(w0,w1,msz,mpi_real8,mpi_sum,comm_grid,ierr)

               do id=1,nd
                  n=irlabel2n(ir,id,ie)
                  unk(n1:n2,n,k,s)=zero
               end do
               c=dble(irdim(ir))/dble(nsym)
               do id=1,nd
                  n=irlabel2n(ir,id,ie)
!                  do i=n1,n2
!                     i1=LL(1,i) ; i2=LL(2,i) ; i3=LL(3,i)
!                     do isym=1,nsym
!                        lt(1)=rga(1,1,isym)*i1+rga(2,1,isym)*i2+rga(3,1,isym)*i3
!                        lt(2)=rga(1,2,isym)*i1+rga(2,2,isym)*i2+rga(3,2,isym)*i3
!                        lt(3)=rga(1,3,isym)*i1+rga(2,3,isym)*i2+rga(3,3,isym)*i3
!                        unk(i,n,k,s)=unk(i,n,k,s)+c*Rir(id,1,isym,ir)*w1( lt(1),lt(2),lt(3) )
!                     end do
!                  end do
                  do isym=1,nsym
                     d=c*Rir(id,1,isym,ir)
                     do i=n1,n2
                        i1=LL(1,i) ; i2=LL(2,i) ; i3=LL(3,i)
                        lt(1)=rga(1,1,isym)*i1+rga(2,1,isym)*i2+rga(3,1,isym)*i3
                        lt(2)=rga(1,2,isym)*i1+rga(2,2,isym)*i2+rga(3,2,isym)*i3
                        lt(3)=rga(1,3,isym)*i1+rga(2,3,isym)*i2+rga(3,3,isym)*i3
                        unk(i,n,k,s)=unk(i,n,k,s)+d*w1( lt(1),lt(2),lt(3) )
                     end do
                  end do
               end do

            end do ! ie

         end do ! ir

      end do ! k
      end do ! s

      mem=mem-bdreal*size(w1) ; deallocate( w1 )
      mem=mem-bdreal*size(w0) ; deallocate( w0 )

      call bwatch(ctime1,etime1)
      if ( DISP_SWITCH ) then
         write(*,*) "TIME(RE_SYM)=",ctime1-ctime0,etime1-etime0
         write(*,*) "MEM(MB)",memax*B2MB,mem
      end if

      return

 900  call stop_program

      END SUBROUTINE re_sym
