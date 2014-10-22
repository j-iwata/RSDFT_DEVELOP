MODULE TestModule
  use parallel_module
  implicit none
CONTAINS
  SUBROUTINE export_DensityAndWF
    use wf_module
    use density_module
    use localpot_module
    implicit none
    integer :: s,k,n,i
    write(9000+myrank,*) MS_0_WF,MS_1_WF,MK_0_WF,MK_1_WF,MB_0_WF,MB_1_WF,ML_0_WF,ML_1_WF
    do s=MS_0_WF,MS_1_WF
      do k=MK_0_WF,MK_1_WF
        do n=MB_0_WF,MB_1_WF
          do i=ML_0_WF,ML_1_WF
            write(9000+myrank,*) unk(i,n,k,s)
          enddo
        enddo
      enddo
    enddo
    write(9100+myrank,*) MS_0_RHO,MS_1_RHO,ML_0_RHO,ML_1_RHO
    do s=MS_0_RHO,MS_1_RHO
      do i=ML_0_RHO,ML_1_RHO
        write(9100+myrank,*) rho(i,s)
      enddo
    enddo
  END SUBROUTINE export_DensityAndWF
!----------------------------------------------
  SUBROUTINE import_DensityAndWF
    use wf_module
    use density_module
    use localpot_module
    use array_bound_module
    implicit none
    integer :: s,k,n,i
    logical :: checkInput(1:8)
    read(9000+myrank,*) MS_0_WF,MS_1_WF,MK_0_WF,MK_1_WF,MB_0_WF,MB_1_WF,ML_0_WF,ML_1_WF
    checkInput(1)=MS_0_WF==MSP_0
    checkInput(2)=MS_1_WF==MSP_1
    checkInput(3)=MK_0_WF==MBZ_0
    checkInput(4)=MK_1_WF==MBZ_1
    checkInput(5)=MB_0_WF==MB_0
    checkInput(6)=MB_1_WF==MB_1
    checkInput(7)=ML_0_WF==ML_0
    checkInput(8)=ML_1_WF==ML_1
    if (.not.any(checkInput(1:8))) write(200+myrank,*) 'something wrong with wf input'
    do s=MS_0_WF,MS_1_WF
      do k=MK_0_WF,MK_1_WF
        do n=MB_0_WF,MB_1_WF
          do i=ML_0_WF,ML_1_WF
            read(9000+myrank,*) unk(i,n,k,s)
          enddo
        enddo
      enddo
    enddo
    read(9100+myrank,*) MS_0_RHO,MS_1_RHO,ML_0_RHO,ML_1_RHO
    checkInput(1)=MS_0_RHO==MSP_0
    checkInput(2)=MS_1_RHO==MSP_1
    checkInput(3)=ML_0_RHO==ML_0
    checkInput(4)=ML_1_RHO==ML_1
    if (.not.any(checkInput(1:4))) write(200+myrank,*) 'something wrong with rho input'
    do s=MS_0_RHO,MS_1_RHO
      do i=ML_0_RHO,ML_1_RHO
        read(9100+myrank,*) rho(i,s)
      enddo
    enddo
  END SUBROUTINE import_DensityAndWF
END MODULE TestModule
