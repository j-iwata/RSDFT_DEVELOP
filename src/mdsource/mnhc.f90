!-----------------------------------------------------------------------
subroutine mnhc
!   use global_variables
   use cpmd_variables
   implicit none
   real(8) :: sc,ekt,ass
   integer jc,iys,ins,jns,nch1

   call calkin(ekt)
   nch1=nch+1
   sc=1.0d0
   gmns(1)=(ekt*2.0d0-fnkt)/qmns(1)
   do jc=1,ndc
      do iys=1,nys
         vmns(nch)=vmns(nch)+gmns(nch)*wt4(iys)
         do ins=1,nch-1
            ass=exp(-wt8(iys)*vmns(nch1-ins))
            jns=nch-ins
            vmns(jns)=vmns(jns)*ass**2+gmns(jns)*ass*wt4(iys)
         enddo
         ass=exp(-wt2(iys)*vmns(1))
         sc=sc*ass 
         gmns(1)=(sc**2*ekt*2.0d0-fnkt)/qmns(1)
         do ins=1,nch
            xmns(ins)=xmns(ins)+vmns(ins)*wt2(iys) 
         enddo
         do ins=1,nch-1
            ass=exp(-wt8(iys)*vmns(ins+1))
            vmns(ins)  = vmns(ins)*ass**2+gmns(ins)*ass*wt4(iys)
            gmns(ins+1)=(qmns(ins)*vmns(ins)**2-fkt)/qmns(ins+1)
         enddo
         vmns(nch)=vmns(nch)+gmns(nch)*wt4(iys)
      enddo
   enddo
   Velocity=Velocity*sc
   return
end subroutine mnhc

subroutine sicoef
!   use global_variables
   use cpmd_variables
   implicit none
   real(8) :: anc,a2,a4,a8
   real(8) :: w1,w2,w3,w4,w5,c0,c1
   real(8) :: wt(7)
   integer  i

   anc=dfloat(ndc)
   a2 =dt/(anc*2.0d0)
   a4 =dt/(anc*4.0d0)
   a8 =dt/(anc*8.0d0)

   if(nys.eq.1) then
      wt2(1)=a2
      wt4(1)=a4
      wt8(1)=a8
   endif

   if(nys.eq.3) then
      c0    =2.0d0**(1.0d0/3.0d0)
      w1    =1.0d0/(2.0d0-c0)
      w2    =1.0d0-2.0d0*w1
      w3    =w1
      wt2(1)=w1*a2
      wt2(2)=w2*a2
      wt2(3)=w3*a2
      wt4(1)=w1*a4
      wt4(2)=w2*a4
      wt4(3)=w3*a4
      wt8(1)=w1*a8
      wt8(2)=w2*a8
      wt8(3)=w3*a8
   endif

   if(nys.eq.5) then
      c1    =4.0d0**(1.0d0/3.0d0)
      w1    =1.0d0/(4.0d0-c1)
      w2    =w1
      w3    =1.0d0-4.0d0*w1
      w4    =w1
      w5    =w1
      wt2(1)=w1 * a2
      wt2(2)=w2 * a2
      wt2(3)=w3 * a2
      wt2(4)=w4 * a2
      wt2(5)=w5 * a2
      wt4(1)=w1 * a4
      wt4(2)=w2 * a4
      wt4(3)=w3 * a4
      wt4(4)=w4 * a4
      wt4(5)=w5 * a4
      wt8(1)=w1 * a8
      wt8(2)=w2 * a8
      wt8(3)=w3 * a8
      wt8(4)=w4 * a8
      wt8(5)=w5 * a8
   endif

   if(nys.eq.7) then
      wt(1) = 0.784513610477560d0
      wt(2) = 0.235573213359357d0
      wt(3) =-1.17767998417887d0
      wt(4) = 1.0d0 - 2.0d0*(wt(1)+wt(2)+wt(3))
      wt(5) =-1.17767998417887d0
      wt(6) = 0.235573213359357d0
      wt(7) = 0.784513610477560d0
      do i=1,7
         wt2(i)=wt(i)*a2
         wt4(i)=wt(i)*a4
         wt8(i)=wt(i)*a8
      enddo
   endif

   return
end subroutine sicoef


subroutine setmnhc(kine)
!   use global_variables
   use cpmd_variables
   implicit none
   real(8),parameter :: wave=100.0d0       !(1/cm)
   real(8) :: randv(nch*6)
   real(8) :: tunit,kine,bathe
   integer ins,ierr

   if(omegan.eq.0.0d0) omegan=wave

   fkt  =temp
   fnkt =fkt*dble(Mi)*3.0d0
   tunit=tojoul/(clight*planck*omegan)

   allocate(xmns(nch),vmns(nch),gmns(nch),qmns(nch))
!   mem_new = mem_new + bdreal*nch*4
   if(myrank==0) then
      call rnum(randv,nch*6)
      if(inivel) then
         do ins=1,nch
            if(ins.eq.1) then
               qmns(ins)=fnkt*tunit*tunit
               xmns(ins)=0.0d0
               vmns(ins)=sqrt(-2.0d0*fnkt*dlog(randv(2*ins-1))/qmns(ins))*sin(pi2*randv(2*ins))
               vmns(ins)=vmns(ins)*sqrt(fnkt/(qmns(ins)*vmns(ins)**2))
               gmns(ins)=(2.0d0*kine-fnkt)/qmns(ins)
            else
               qmns(ins)=fkt*tunit*tunit
               xmns(ins)=0.0d0
               vmns(ins)=sqrt(-2.0d0*fkt*dlog(randv(2*ins-1))/qmns(ins))*sin(pi2*randv(2*ins))
               vmns(ins)=vmns(ins)*sqrt(fkt/(qmns(ins)*vmns(ins)**2))
               gmns(ins)=(qmns(ins-1)*vmns(ins-1)**2-fkt)/qmns(ins)
            endif
         enddo
      else
         write(*,*) "read bath coordinates from file"
         open(1,file='bathcoordinate.dat0',status='unknown')
         read(1,*) bathe
         read(1,*) (xmns(ins),ins=1,nch)
         read(1,*) (vmns(ins),ins=1,nch)
         read(1,*) (gmns(ins),ins=1,nch)
         read(1,*) (qmns(ins),ins=1,nch)
         close(1)
      endif
      write(*,1000) fkt, fnkt, fkt * tunit * tunit
      write(*,1001) tunit,1.0d0/tunit
   endif
   call mpi_bcast(xmns,nch,mpi_real8,0,mpi_comm_world,ierr)
   call mpi_bcast(vmns,nch,mpi_real8,0,mpi_comm_world,ierr)
   call mpi_bcast(gmns,nch,mpi_real8,0,mpi_comm_world,ierr)
   call mpi_bcast(qmns,nch,mpi_real8,0,mpi_comm_world,ierr)

   return
1000 format("fkt  =",f12.6,2x,"fnkt=",f12.6,2x,"qmns=",f12.6)
1001 format("tunit=",f12.6,2x,"De  =",f12.6)
   contains 
      subroutine rnum(randv,N)
      implicit none
      real(8),intent(inout) :: randv(N)
      integer,intent(in)    :: N
      real(8),parameter :: a=1953125.0d0
      real(8),parameter :: b=33554432.0d0
      real(8),parameter :: r0=177147.0d0
      real(8) r
      integer i
      r=r0
      do i=1,N
        r=mod(a*r,b)
        randv(i)=r/b
      enddo
      return
   end subroutine rnum
end subroutine setmnhc

subroutine unsetmnhc(bathe)
!   use global_variables
   use cpmd_variables
   implicit none
   real(8),intent(in) :: bathe
   integer :: k

   if(myrank==0) then
      write(*,'(2f15.8)') bathe
      write(*,'(8f15.8)') (xmns(k),k=1,nch)
      write(*,'(8f15.8)') (vmns(k),k=1,nch)
      write(*,'(8f15.8)') (gmns(k),k=1,nch)
      write(*,'(8f15.8)') (qmns(k),k=1,nch)
      open(1,file='bathcoordinate.dat1',status='unknown')
      write(1,*) bathe
      write(1,*) (xmns(k),k=1,nch)
      write(1,*) (vmns(k),k=1,nch)
      write(1,*) (gmns(k),k=1,nch)
      write(1,*) (qmns(k),k=1,nch)
      close(1)
   endif
   deallocate(xmns,vmns,gmns,qmns)
!   mem_new = mem_new - bdreal*nch*4

   return
end subroutine unsetmnhc

subroutine emnhc(bathe)
!   use global_variables
   use cpmd_variables
   implicit none
   real(8),intent(inout) :: bathe
   real(8) :: bathk,bathp
   integer ins

   bathk=0.0d0
   do ins=1,nch
      bathk=bathk+vmns(ins)**2*qmns(ins)
   enddo

   bathp=xmns(1)*fnkt
   do ins=2,nch
      bathp=bathp+xmns(ins)*fkt
   enddo

   bathk=bathk*0.5d0
   bathe=bathk+bathp

   return
end subroutine emnhc
