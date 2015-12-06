MODULE cpmdio2_module

  implicit none

  PRIVATE
  PUBLIC :: write_data_cpmd_k_para,read_data_cpmd_k_para,gather_wf_data

CONTAINS

!-----------------------------------------------------------------------------------------
      subroutine write_data_cpmd_k_para
      use cpmd_variables
      use array_bound_module
      use wf_module
      implicit none
      integer :: n1,n2,ML0,n,k,i,ispin
      character(len=64) :: filename

      write(filename,"('restart_',i5.5,'.dat')") myrank

      open(1,file=filename,form="unformatted")
      n1  = idisp(myrank)+1
      n2  = idisp(myrank)+ircnt(myrank)
      ML0 = ircnt(myrank)

!      write(*,*) "MSP_0,MSP_1", MSP_0,MSP_1
!      write(*,*) "MBZ_0,MBZ_1",MBZ_0,MBZ_1

      do ispin=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
      do n=MB_0_CPMD,MB_1_CPMD
      do i=n1,n2
        write(1) unk(i,n,k,ispin)
      enddo
      enddo
      enddo
      enddo

      do ispin=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
      do n=MB_0_CPMD,MB_1_CPMD
      do i=n1,n2
        write(1) psi_v(i,n,k,ispin)
      enddo
      enddo
      enddo
      enddo

      close(1)
      return
      end subroutine write_data_cpmd_k_para

!--------------------------------------------------------------------
!--------------------------------------------------------------------
      subroutine read_data_cpmd_k_para
      use cpmd_variables
      use array_bound_module
      use wf_module
      implicit none
      integer :: n1,n2,ML0,n,k,i,ispin
      character(len=64) :: filename

      write(filename,"('restart_',i5.5,'.dat')") myrank

      open(1,file=filename,form="unformatted")
      n1  = idisp(myrank)+1
      n2  = idisp(myrank)+ircnt(myrank)
      ML0 = ircnt(myrank)

      do ispin=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
      do n=MB_0_CPMD,MB_1_CPMD
      do i=n1,n2
        read(1) unk(i,n,k,ispin)
      enddo
      enddo
      enddo
      enddo

      do ispin=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
      do n=MB_0_CPMD,MB_1_CPMD
      do i=n1,n2
        read(1) psi_v(i,n,k,ispin)
      enddo
      enddo
      enddo
      enddo

      close(1)
      return
      end subroutine read_data_cpmd_k_para


subroutine gather_wf_data(iswitch_send)
   use cpmd_variables
   use parallel_module
   use array_bound_module
   use wf_module
   implicit none

   integer,allocatable :: id(:),ir(:)
   integer :: n1,n2,ML0,s,k,mrnk,memax,mem,ierr
   integer :: iswitch_send,MBW_0

   n1    = idisp(myrank)+1
   n2    = idisp(myrank)+ircnt(myrank)
   ML0   = ircnt(myrank)
   memax = 0.d0
   mem   = 0.d0

!- allocate ----------------------------------------------------
   allocate( id(0:np_band-1),ir(0:np_band-1) ) ; id=0 ; ir=0
!---------------------------------------------------------------

   mrnk = id_class(myrank,4)

   if(iswitch_send==0) then

      id(0:np_band-1) = id_band(0:np_band-1)*ML0
      ir(0:np_band-1) = ir_band(0:np_band-1)*ML0
      MBW_0=MB_0
      do s=MSP_0,MSP_1
         do k=MBZ_0,MBZ_1
            call mpi_allgatherv(  unk(n1,MBW_0,k,s),ir(mrnk),TYPE_MAIN,unk(n1,1,k,s),ir,id,TYPE_MAIN,comm_band,ierr)
         enddo
      enddo
   else

      id(0:np_band-1) = id_band_cpmd(0:np_band-1)*ML0
      ir(0:np_band-1) = ir_band_cpmd(0:np_band-1)*ML0
      MBW_0=MB_0_CPMD
      do s=MSP_0,MSP_1
         do k=MBZ_0,MBZ_1
            call mpi_allgatherv(  unk(n1,MBW_0,k,s),ir(mrnk),TYPE_MAIN,unk(n1,1,k,s),ir,id,TYPE_MAIN,comm_band,ierr)
            call mpi_allgatherv(psi_v(n1,MBW_0,k,s),ir(mrnk),TYPE_MAIN,psi_v(n1,1,k,s),ir,id,TYPE_MAIN,comm_band,ierr)
            call mpi_allgatherv(psi_n(n1,MBW_0,k,s),ir(mrnk),TYPE_MAIN,psi_n(n1,1,k,s),ir,id,TYPE_MAIN,comm_band,ierr)
         enddo
      enddo
   endif

!- allocate ----------------------------------------------------
   deallocate( id,ir )
!---------------------------------------------------------------

   return
end subroutine 
end MODULE cpmdio2_module
