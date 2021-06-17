module io_read_wf_general_module

  use parallel_module, only: para, myrank_s,myrank_k,myrank_b,myrank_g
  use watch_module, only: watchb
  use wf_module, only: unk

  implicit none

  private
  public :: read_wf_general

contains

  subroutine read_wf_general( pinfo0,pinfo1,ML1,ML2,Ml3,MB0,MBZ0,MSP0,MB1,MBZ1,MSP1,type_wf )
    implicit none
    type(para),intent(in) :: pinfo0, pinfo1
    integer,intent(in) :: ML1,ML2,ML3,MB0,MBZ0,MSP0,MB1,MBZ1,MSP1
    integer,intent(in) :: type_wf
    real(8) :: ttmp(2),tt(2)
    logical :: disp
    integer :: p1,p2,p3,p4,p5,p6,irank_g,irank_g_old,g0,g1,g2,g3
    integer :: s0,s1,k0,k1,b0,b1,i1_0,i1_1,i2_0,i2_1,i3_0,i3_1
    integer :: i1_2,i1_3,i2_2,i2_3,i3_2,i3_3
    integer :: s,k,b,i1,i2,i3,irank_old,s2,s3,k2,k3,b2,b3,n,i,ng
    integer :: ig, ig_old, n11, n12
    character(13) :: filename
    character(7) :: basename
    character(5) :: crank
    integer,allocatable :: ir_new(:,:), id_new(:,:)
    integer,allocatable :: ir_old(:,:), id_old(:,:)
    real(8),allocatable :: dtmp(:)
    complex(8),allocatable :: utmp(:)
    integer,parameter :: u=3

    call write_border( 0, " read_wf_general(start)" )
    call check_disp_switch( disp, 0 )
    call watchb( ttmp, barrier='on' ); tt=0.0d0

    basename = 'wf.dat1'

    s2 = pinfo1%spin%id(myrank_s)+1
    s3 = s2 + pinfo1%spin%ir(myrank_s)-1
    k2 = pinfo1%bzsm%id(myrank_k)+1
    k3 = k2 + pinfo1%bzsm%ir(myrank_k)-1
    b2 = pinfo1%band%id(myrank_b)+1
    b3 = b2 + pinfo1%band%ir(myrank_b)-1
    g2 = pinfo1%grid%id(myrank_g)+1
    g3 = g2 + pinfo1%grid%ir(myrank_g)-1

    n = maxval( pinfo1%np(1:3) )
    allocate( ir_new(0:n-1,3) ); ir_new=0
    allocate( id_new(0:n-1,3) ); id_new=0
    do i1 = 0, ML1-1
       n = mod( i1, pinfo1%np(1) )
       ir_new(n,1) = ir_new(n,1) + 1
    end do
    do i2 = 0, ML2-1
       n = mod( i2, pinfo1%np(2) )
       ir_new(n,2) = ir_new(n,2) + 1
    end do
    do i3 = 0, ML3-1
       n = mod( i3, pinfo1%np(3) )
       ir_new(n,3) = ir_new(n,3) + 1
    end do
    do i = 1, 3
       do n = 0, pinfo1%np(i)-1
          id_new(n,i) = sum( ir_new(0:n,i) ) - ir_new(n,i)
       end do
    end do

    n = maxval( pinfo0%np(1:3) )
    allocate( ir_old(0:n-1,3) ); ir_old=0
    allocate( id_old(0:n-1,3) ); id_old=0
    do i1 = 0, ML1-1
       n = mod( i1, pinfo0%np(1) )
       ir_old(n,1) = ir_old(n,1) + 1
    end do
    do i2 = 0, ML2-1
       n = mod( i2, pinfo0%np(2) )
       ir_old(n,2) = ir_old(n,2) + 1
    end do
    do i3 = 0, ML3-1
       n = mod( i3, pinfo0%np(3) )
       ir_old(n,3) = ir_old(n,3) + 1
    end do
    do i = 1, 3
       do n = 0, pinfo0%np(i)-1
          id_old(n,i) = sum( ir_old(0:n,i) ) - ir_old(n,i)
       end do
    end do

    loop_i3: do i3 = 0, pinfo1%np(3)-1
             do i2 = 0, pinfo1%np(2)-1
             do i1 = 0, pinfo1%np(1)-1
               irank_g = i1 + i2*pinfo1%np(1) + i3*pinfo1%np(1)*pinfo1%np(2)
               if ( irank_g == myrank_g ) then
                 i1_2 = id_new(i1,1)
                 i1_3 = id_new(i1,1) + ir_new(i1,1) - 1
                 i2_2 = id_new(i2,2)
                 i2_3 = id_new(i2,2) + ir_new(i2,2) - 1
                 i3_2 = id_new(i3,3)
                 i3_3 = id_new(i3,3) + ir_new(i3,3) - 1
                 n11 = i1_3 - i1_2 + 1
                 n12 = n11*(i2_3-i2_2+1)
                 exit loop_i3
               end if
             end do
             end do
             end do loop_i3

    n=maxval(pinfo0%grid%ir)
    if ( type_wf == 1 ) then ! real-wf
       allocate( dtmp(n) ); dtmp=0.0d0
    else ! complex-wf
       allocate( utmp(n) ); utmp=(0.0d0,0.0d0)
    end if

    irank_old = -1

    do p6 = 0, pinfo0%np(6)-1
       s0 = pinfo0%spin%id(p6)+1
       s1 = s0 + pinfo0%spin%ir(p6)-1
    do p5 = 0, pinfo0%np(5)-1
       k0 = pinfo0%bzsm%id(p5)+1
       k1 = k0 + pinfo0%bzsm%ir(p5)-1
    do p4 = 0, pinfo0%np(4)-1
       b0 = pinfo0%band%id(p4)+1
       b1 = b0 + pinfo0%band%ir(p4)-1

       irank_g_old = -1

       do p3 = 0, pinfo0%np(3)-1
          i3_0 = id_old(p3,3)
          i3_1 = i3_0 + ir_old(p3,3) - 1
       do p2 = 0, pinfo0%np(2)-1
          i2_0 = id_old(p2,2)
          i2_1 = i2_0 + ir_old(p2,2) - 1
       do p1 = 0, pinfo0%np(1)-1
          i1_0 = id_old(p1,1)
          i1_1 = i1_0 + ir_old(p1,1) - 1

          irank_g_old = irank_g_old + 1
          g0 = pinfo0%grid%id(irank_g_old)+1
          g1 = g0 + pinfo0%grid%ir(irank_g_old)-1

          irank_old = irank_old + 1

          write(crank,'(i5.5)') irank_old
          filename = basename//'.'//crank
          write(*,*) irank_old, filename

          if ( ( (s0 <= s2 .and. s2 <= s1) .or. (s0 <= s3 .and. s3 <= s1) ) .and. &
               ( (k0 <= k2 .and. k2 <= k1) .or. (k0 <= k3 .and. k3 <= k1) ) .and. &
               ( (b0 <= b2 .and. b2 <= b1) .or. (b0 <= b3 .and. b3 <= b1) ) .and. &
               ( (i1_0 <= i1_2 .and. i1_2 <= i1_1) .or. (i1_0 <= i1_3 .and. i1_3 <= i1_1) ) .and. &
               ( (i2_0 <= i2_2 .and. i2_2 <= i2_1) .or. (i2_0 <= i2_3 .and. i2_3 <= i2_1) ) .and. &
               ( (i3_0 <= i3_2 .and. i3_2 <= i3_1) .or. (i3_0 <= i3_3 .and. i3_3 <= i3_1) ) ) then

             open(u,file=filename,status='old',form='unformatted')

             ng = g1 - g0 + 1

             do s = s0, s1
             do k = k0, k1
             do b = b0, b1

                if ( type_wf == 1 ) then ! real-wf

                   read(u) dtmp(1:ng)

                   if ( s < s2 .or. s3 < s .or. &
                        k < k2 .or. k3 < k .or. &
                        b < b2 .or. b3 < b ) cycle

                   ig_old=0
                   do i3 = i3_0, i3_1
                   do i2 = i2_0, i2_1
                   do i1 = i1_0, i1_1
                      ig_old = ig_old + 1
                      if ( i1 < i1_2 .or. i1_3 < i1 .or. &
                           i2 < i2_2 .or. i2_3 < i2 .or. &
                           i3 < i3_2 .or. i3_3 < i3 ) cycle
                      ig = g2 + i1-i1_2 + n11*(i2-i2_2) + n12*(i3-i3_2)
                      unk(ig,b,k,s)=dtmp(ig_old)
                   end do
                   end do
                   end do

                else ! complex-wf

                   read(u) utmp(1:ng)

                   if ( s < s2 .or. s3 < s .or. &
                        k < k2 .or. k3 < k .or. &
                        b < b2 .or. b3 < b ) cycle

                   ig_old=0
                   do i3 = i3_0, i3_1
                   do i2 = i2_0, i2_1
                   do i1 = i1_0, i1_1
                      ig_old = ig_old + 1
                      if ( i1 < i1_2 .or. i1_3 < i1 .or. &
                           i2 < i2_2 .or. i2_3 < i2 .or. &
                           i3 < i3_2 .or. i3_3 < i3 ) cycle
                      ig = g2 + i1-i1_2 + n11*(i2-i2_2) + n12*(i3-i3_2)
                      unk(ig,b,k,s)=utmp(ig_old)
                   end do
                   end do
                   end do

                end if

             end do !b
             end do !k
             end do !s

             close(u)

          end if

       end do !p1
       end do !p2
       end do !p3

    end do !p4
    end do !p5
    end do !p6

    deallocate( id_old ) 
    deallocate( ir_old ) 
    deallocate( id_new ) 
    deallocate( ir_new ) 
    if ( allocated(dtmp) ) deallocate(dtmp)
    if ( allocated(utmp) ) deallocate(utmp)

    call watchb( ttmp, tt, barrier='on' )
    if ( disp ) then
       write(*,*) "time(read_wf_general)=",tt
    end if
    call write_border( 0, " read_wf_general(end)" )

  end subroutine read_wf_general

end module io_read_wf_general_module
