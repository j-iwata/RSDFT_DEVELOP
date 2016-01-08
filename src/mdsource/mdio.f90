!-----------------------------------------------------------------------
!     IO routine for molecular dynamics
!-----------------------------------------------------------------------
subroutine mdio(io,tote0)
   use global_variables
   use cpmd_variables
   implicit none
   integer :: io
   real(8) :: tote0
   integer i,k
   real(8) dum(3)
   character*50 aline
   if(io==0) then
      open(2,file='initial.dat',status='unknown')
      read(2,'(a)') aline
      do i=1,Mi
         read(2,'(3f24.16)') (Rion(k,i),k=1,3)
      enddo
      read(2,'(a)') aline
      do i=1,Mi
         read(2,'(3f24.16)') (Velocity(k,i),k=1,3)
      enddo
      read(2,'(a)') aline
      do i=1,Mi
         read(2,'(3f24.16)') (dum(k),k=1,3)
      enddo
      read(2,*) tote0
      close(2)
   elseif(io==1) then
      open(2,file='final.dat',status='unknown')
      write(2,'(a)') "final coordinate"
      do i=1,Mi
         write(2,'(3g24.16)') (Rion(k,i),k=1,3)
      enddo
      write(2,'(a)') "final velocity"
      do i=1,Mi
         write(2,'(3g24.16)') (Velocity(k,i),k=1,3)
      enddo
      write(2,'(a)') "final force"
      do i=1,Mi
         write(2,'(3g24.16)') (Force(k,i),k=1,3)
      enddo
      write(2,*) tote0
      close(2)
   elseif(io==2) then
      write(*,'(a)') "initial coordinate"
      do i=1,Mi
         write(*,'(3f15.8)') (Rion(k,i),k=1,3)
      enddo
      write(*,'(a)') "initial velocity"
      do i=1,Mi
         write(*,'(3f15.8)') (Velocity(k,i),k=1,3)
      enddo
      write(*,'(a)') "initial force"
      do i=1,Mi
         write(*,'(3f15.8)') (Force(k,i),k=1,3)
      enddo
      write(*,*) tote0
   endif
   return
end subroutine mdio
