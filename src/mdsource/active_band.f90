!-----------------------------------------------------------------------
!     detect active bands
!-----------------------------------------------------------------------
subroutine active_band
  use wf_module, only: occ
  use cpmd_variables,only:mstocck,MBC,MB,MB_0,MB_1,MBZ_0,MBZ_1,MSP_0,MSP_1 &
       ,MB_0_CPMD,MB_1_CPMD,ir_band_cpmd,id_band_cpmd
  use parallel_module, only: myrank,myrank_b,np_band
  implicit none 
  integer :: s,k,n

  allocate( mstocck(MBZ_0:MBZ_1,MSP_0:MSP_1) )
  mstocck=0

  do s=MSP_0,MSP_1
  do k=MBZ_0,MBZ_1
  do n=1,MB
     if ( occ(n,k,s)>1.d-10 ) mstocck(k,s)=mstocck(k,s)+1
  enddo
  enddo
  enddo

  MBC=maxval(mstocck)

  if ( myrank == 0 ) then
     do s=MSP_0,MSP_1
     do k=MBZ_0,MBZ_1
        write(*,*) "occupation number of ",k,"-th BZ=",mstocck(k,s)
     enddo
     enddo
     write(*,*) "maximum occupation number of MB=",mbc
  endif

!--- for band-parallel calculation

  allocate( ir_band_cpmd(0:np_band-1) ) ; ir_band_cpmd= 0
  allocate( id_band_cpmd(0:np_band-1) ) ; id_band_cpmd=-1

  do n=1,MBC
     k=mod(n,np_band)
     ir_band_cpmd(k) = ir_band_cpmd(k) + 1
  end do
  do k=0,np_band-1
     id_band_cpmd(k) = sum( ir_band_cpmd(0:k) ) - ir_band_cpmd(k)
  end do

  MB_0_CPMD = id_band_cpmd(myrank_b) + 1
  MB_1_CPMD = id_band_cpmd(myrank_b) + ir_band_cpmd(myrank_b)

  return
end subroutine active_band
