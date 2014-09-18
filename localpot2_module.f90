MODULE localpot2_module

  use rgrid_module
  use localpot2_variables
  use localpot2_ion_module
  use localpot2_density_module
  use localpot2_vh_module
  use localpot2_xc_module
  use watch_module
  use parallel_module
  use array_bound_module

  implicit none

  PRIVATE
  PUBLIC :: test2_localpot2, Lpot, vloc_nl, MLpot &
           ,read_localpot2, flag_localpot2 &
           ,prep_parallel_localpot2, init_localpot2

  logical :: flag_localpot2 = .false.

  integer :: Nintp_loc

  integer :: MLpot
  integer,allocatable :: Lpot(:,:)
  real(8),allocatable :: vloc_nl(:,:)

CONTAINS


  SUBROUTINE init_localpot2
    implicit none
    integer :: m1,m2,m3,mm1,mm2,mm3

    call read_localpot2(1,myrank)

    if ( flag_localpot2 ) then
       if ( disp_switch_parallel ) then
          write(*,'(a60," Init_LOCALPOT2(START)")') repeat("-",60)
       end if
       call prep_parallel_localpot2
       call initsub_localpot2(Ngrid(1),Ngrid(2),Ngrid(3))
       m1=Ngrid_dense(1)
       m2=Ngrid_dense(2)
       m3=Ngrid_dense(3)
       mm1=Igrid_dense(2,1)-Igrid_dense(1,1)+1
       mm2=Igrid_dense(2,2)-Igrid_dense(1,2)+1
       mm3=Igrid_dense(2,3)-Igrid_dense(1,3)+1
       call test_localpot2
       call test2_localpot2( vion_nl )
    end if
       
  END SUBROUTINE init_localpot2


  SUBROUTINE read_localpot2(unit,rank)
    implicit none
    integer,intent(IN) :: unit,rank
    integer :: i
    character(5) :: cbuf,ckey
    flag_localpot2=.false.
    Ndens_loc=1
    Nintp_loc=1
    fecut_loc=1.0d0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey == "NDLOC" ) then
             backspace(unit)
             read(unit,*) cbuf,Ndens_loc,Nintp_loc,fecut_loc
             exit
          end if
       end do
999    continue
       if ( Ndens_loc > 1 ) flag_localpot2=.true. 
       write(*,*) "Ndens_loc=",Ndens_loc
       write(*,*) "Nintp_loc=",Nintp_loc
       write(*,*) "fecut_loc=",fecut_loc
       write(*,*) "flag_localpot2=",flag_localpot2
    end if
    call send_localpot2
  END SUBROUTINE read_localpot2


  SUBROUTINE send_localpot2
    implicit none
    integer :: ierr
    call mpi_bcast(Ndens_loc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Nintp_loc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(fecut_loc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(flag_localpot2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_localpot2


  SUBROUTINE initsub_localpot2(n1,n2,n3)
    implicit none
    integer,intent(IN) :: n1,n2,n3
    integer :: m1,m2,m3,m1_0,m1_1,m2_0,m2_1,m3_0,m3_1
    Ngrid_dense(1) = n1*Ndens_loc
    Ngrid_dense(2) = n2*Ndens_loc
    Ngrid_dense(3) = n3*Ndens_loc
    Ngrid_dense(0) = Ngrid_dense(1)*Ngrid_dense(2)*Ngrid_dense(3)
    m1_0=Igrid_dense(1,1)
    m1_1=Igrid_dense(2,1)
    m2_0=Igrid_dense(1,2)
    m2_1=Igrid_dense(2,2)
    m3_0=Igrid_dense(1,3)
    m3_1=Igrid_dense(2,3)
    allocate( vion_nl(m1_0:m1_1,m2_0:m2_1,m3_0:m3_1) ) ; vion_nl=0.0d0
    allocate(   vh_nl(m1_0:m1_1,m2_0:m2_1,m3_0:m3_1) ) ;   vh_nl=0.0d0
    allocate(  vxc_nl(m1_0:m1_1,m2_0:m2_1,m3_0:m3_1) ) ;  vxc_nl=0.0d0
    allocate(  rho_nl(m1_0:m1_1,m2_0:m2_1,m3_0:m3_1) ) ;  rho_nl=0.0d0
    allocate( vloc_dense(m1_0:m1_1,m2_0:m2_1,m3_0:m3_1) ) ; vloc_dense=0.0d0
    allocate( vloc_dense_old(m1_0:m1_1,m2_0:m2_1,m3_0:m3_1) )
    vloc_dense_old=0.0d0
  END SUBROUTINE initsub_localpot2


  SUBROUTINE test_localpot2
    implicit none
    integer :: ll_0,ll_1,i1,i2,j1
    real(8) :: x,y

    nitp_0 = min( 0, -Nintp_loc+1 )
    nitp_1 = Nintp_loc

    allocate( Clag1(nitp_0:nitp_1,0:Ndens_loc-1) ) ; Clag1=0.0d0
    allocate( Clag2(nitp_0:nitp_1,0:Ndens_loc-1) ) ; Clag2=0.0d0
    allocate( Clag3(nitp_0:nitp_1,0:Ndens_loc-1) ) ; Clag3=0.0d0

    do j1=0,Ndens_loc-1

       do i1=nitp_0,nitp_1

          Clag1(i1,j1)=1.d0

          do i2=nitp_0,nitp_1

             if ( i2 == i1 ) cycle

             x=dble(j1-i2*Ndens_loc)
             y=dble(i1-i2)*Ndens_loc

             Clag1(i1,j1)=Clag1(i1,j1)*(x/y)

          end do ! i2

       end do ! i1

    end do ! j1

    Clag2(:,:) = Clag1(:,:)
    Clag3(:,:) = Clag1(:,:)

!-

    dV_dense = dV*Ngrid(0)/Ngrid_dense(0)

  END SUBROUTINE test_localpot2


  SUBROUTINE prep_parallel_localpot2
    implicit none
    integer :: i
    do i=1,3
       Igrid_dense(1,i)=Igrid(1,i)*Ndens_loc
       Igrid_dense(2,i)=Igrid(2,i)*Ndens_loc+Ndens_loc-1
    end do
  END SUBROUTINE prep_parallel_localpot2


  SUBROUTINE test2_localpot2( vpot )
    implicit none
    real(8),intent(IN) :: vpot(:,:,:)

    integer :: ic1,ic2,ic3,jd1,jd2,jd3,id1,id2,id3,ML,MK,i,ic,it,jt
    integer :: ML1,ML2,ML3,itp1,itp2,itp3,jtp1,jtp2,jtp3,i1,i2,i3
    integer :: iic1,iic2,iic3,j,jc1,jc2,jc3,k,irank,ierr,ichk0
    real(8) :: ct0,ct1,et0,et1,ct2,ct3,et2,et3
    real(8) :: const,v,vc
    real(8),allocatable :: w(:,:,:),atmp(:)
    integer,allocatable :: LLL(:,:,:),KKK(:,:,:),ichk(:)
    integer,allocatable :: ir(:),id(:)

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    call watch(ct0,et0)

    allocate( LLL(0:ML1-1,0:ML2-1,0:ML3-1) ) ; LLL=0
    irank=-1
    i=0
    do i3=1,node_partition(3)
    do i2=1,node_partition(2)
    do i1=1,node_partition(1)
       irank=irank+1
       do ic3=pinfo_grid(5,irank),pinfo_grid(5,irank)+pinfo_grid(6,irank)-1
       do ic2=pinfo_grid(3,irank),pinfo_grid(3,irank)+pinfo_grid(4,irank)-1
       do ic1=pinfo_grid(1,irank),pinfo_grid(1,irank)+pinfo_grid(2,irank)-1
          i=i+1
          LLL(ic1,ic2,ic3)=i
       end do
       end do
       end do
    end do
    end do
    end do
    ML=i

    allocate( KKK(nitp_0:nitp_1,nitp_0:nitp_1,nitp_0:nitp_1) ) ; KKK=0
    i=0
    do itp3=nitp_0,nitp_1
    do itp2=nitp_0,nitp_1
    do itp1=nitp_0,nitp_1
       i=i+1
       KKK(itp1,itp2,itp3)=i
    end do
    end do
    end do
    MK=i

    if ( disp_switch_parallel ) write(*,*) "ML,MK",ML,MK

    if ( disp_switch_parallel ) then
       call watch(ct1,et1) ; write(*,*) "(1)",ct1-ct0,et1-et0
    end if

    allocate( w(MK,MK,ML) ) ; w=0.0d0

    do ic3=Igrid(1,3),Igrid(2,3)
    do jd3=0,Ndens_loc-1
       id3=ic3*Ndens_loc+jd3-Igrid_dense(1,3)+1

       do ic2=Igrid(1,2),Igrid(2,2)
       do jd2=0,Ndens_loc-1
          id2=ic2*Ndens_loc+jd2-Igrid_dense(1,2)+1

          do ic1=Igrid(1,1),Igrid(2,1)
          do jd1=0,Ndens_loc-1
             id1=ic1*Ndens_loc+jd1-Igrid_dense(1,1)+1

             ic = LLL(ic1,ic2,ic3)
             v = vpot(id1,id2,id3)

             do jtp3=nitp_0,nitp_1
             do jtp2=nitp_0,nitp_1
             do jtp1=nitp_0,nitp_1

                jt = KKK(jtp1,jtp2,jtp3)

                vc = v*Clag1(jtp1,jd1)*Clag2(jtp2,jd2)*Clag3(jtp3,jd3)

             do itp3=nitp_0,nitp_1
             do itp2=nitp_0,nitp_1
             do itp1=nitp_0,nitp_1

                it = KKK(itp1,itp2,itp3)

                w(it,jt,ic) = w(it,jt,ic) &
                     + Clag1(itp1,jd1)*Clag2(itp2,jd2)*Clag3(itp3,jd3) * vc

             end do ! itp1
             end do ! itp2
             end do ! itp3

             end do ! jtp1
             end do ! jtp2
             end do ! jtp3

          end do
          end do

       end do
       end do

    end do
    end do

    if ( disp_switch_parallel ) then
    call watch(ct0,et0) ; write(*,*) "(2)",ct0-ct1,et0-et1
    end if

    allocate( ir(0:np_grid-1) ) ; ir=0
    allocate( id(0:np_grid-1) ) ; id=0

    ir(:) = MK*MK*ir_grid(:)
    id(:) = MK*MK*id_grid(:)

    call mpi_allgatherv( w(1,1,ML_0),ir(myrank_g),mpi_real8 &
            ,w(1,1,1),ir,id,mpi_real8,comm_grid,ierr )

    deallocate( id )
    deallocate( ir )

    if ( disp_switch_parallel ) then
    call watch(ct1,et1) ; write(*,*) "(2-2)",ct1-ct0,et1-et0
    end if

!
! ---
!

    if ( MLpot == 0 ) then

       call watch(ct0,et0)

       allocate( ichk(ML) ) ; ichk=0

       do iic3=0,ML3-1
       do iic2=0,ML2-1
       do iic1=0,ML1-1

          ic = LLL(iic1,iic2,iic3)

          do itp3=nitp_0,nitp_1
          do itp2=nitp_0,nitp_1
          do itp1=nitp_0,nitp_1

             ic1 = mod(iic1+itp1+ML1,ML1)
             ic2 = mod(iic2+itp2+ML2,ML2)
             ic3 = mod(iic3+itp3+ML3,ML3)
             i   = LLL(ic1,ic2,ic3)

             if ( i /= 1 ) cycle

          do jtp3=nitp_0,nitp_1
          do jtp2=nitp_0,nitp_1
          do jtp1=nitp_0,nitp_1

             jc1 = mod(iic1+jtp1+ML1,ML1)
             jc2 = mod(iic2+jtp2+ML2,ML2)
             jc3 = mod(iic3+jtp3+ML3,ML3)
             j   = LLL(jc1,jc2,jc3)

             ichk(j) = ichk(j) + 1

          end do
          end do
          end do

          end do
          end do
          end do

       end do
       end do
       end do

       MLpot = count( ichk /= 0 )

       if ( disp_switch_parallel ) then
       write(*,*) "MLpot=",MLpot
       write(*,*) "sum(ichk)=",sum(ichk)
       end if

       deallocate( ichk )

       allocate( Lpot(MLpot,ML_0:ML_1)    ) ; Lpot=0
       allocate( vloc_nl(MLpot,ML_0:ML_1) ) ; vloc_nl=0.0d0

       if ( disp_switch_parallel ) then
       call watch(ct1,et1) ; write(*,*) "(3)",ct1-ct0,et1-et0
       end if

! ---

       allocate( ichk(ML) ) ; ichk=0

       do iic3=0,ML3-1
       do iic2=0,ML2-1
       do iic1=0,ML1-1

          ic = LLL(iic1,iic2,iic3)

          do itp3=nitp_0,nitp_1
          do itp2=nitp_0,nitp_1
          do itp1=nitp_0,nitp_1

             ic1 = mod(iic1+itp1+ML1,ML1)
             ic2 = mod(iic2+itp2+ML2,ML2)
             ic3 = mod(iic3+itp3+ML3,ML3)
             it  = KKK(itp1,itp2,itp3)
             i   = LLL(ic1,ic2,ic3)

             if ( i < ML_0 .or. ML_1 < i ) cycle

             ichk0 = ichk(i)

          do jtp3=nitp_0,nitp_1
          do jtp2=nitp_0,nitp_1
          do jtp1=nitp_0,nitp_1

             jc1 = mod(iic1+jtp1+ML1,ML1)
             jc2 = mod(iic2+jtp2+ML2,ML2)
             jc3 = mod(iic3+jtp3+ML3,ML3)
             jt  = KKK(jtp1,jtp2,jtp3)
             j   = LLL(jc1,jc2,jc3)

             do k=1,ichk0
                if ( Lpot(k,i) == j ) exit
             end do
             if ( k > ichk0 ) then
                ichk0 = ichk0 + 1
                Lpot(ichk0,i) = j
             end if

          end do ! jtp1
          end do ! jtp2
          end do ! jtp3

          ichk(i) = ichk0

          end do ! itp1
          end do ! itp2
          end do ! itp3

       end do ! ic1
       end do ! ic2
       end do ! ic3

       deallocate( ichk )

       if ( disp_switch_parallel ) then
       call watch(ct0,et0) ; write(*,*) "(4)",ct0-ct1,et0-et1
       end if

! ---

       allocate( ichk(MLpot) ) ; ichk=0
       allocate( atmp(MLpot) ) ; atmp=0.0d0

       do i=ML_0,ML_1
          atmp(:) = Lpot(:,i)
          call indexx( MLpot, atmp, ichk )
          do k=1,MLpot
             Lpot(k,i) = nint( atmp(ichk(k)) )
          end do ! k
       end do ! i

       deallocate( atmp )
       deallocate( ichk )

       if ( disp_switch_parallel ) then
       call watch(ct1,et1) ; write(*,*) "(4-2)",ct1-ct0,et1-et0
       end if

    end if

! ---

    vloc_nl(:,:) = 0.0d0

! ---

    call watch(ct0,et0)

    do iic3=0,ML3-1
    do iic2=0,ML2-1
    do iic1=0,ML1-1

       ic = LLL(iic1,iic2,iic3)

       do itp3=nitp_0,nitp_1
       do itp2=nitp_0,nitp_1
       do itp1=nitp_0,nitp_1

          ic1 = mod(iic1+itp1+ML1,ML1)
          ic2 = mod(iic2+itp2+ML2,ML2)
          ic3 = mod(iic3+itp3+ML3,ML3)
          it  = KKK(itp1,itp2,itp3)
          i   = LLL(ic1,ic2,ic3)

          if ( i < ML_0 .or. ML_1 < i ) cycle

       do jtp3=nitp_0,nitp_1
       do jtp2=nitp_0,nitp_1
       do jtp1=nitp_0,nitp_1

          jc1 = mod(iic1+jtp1+ML1,ML1)
          jc2 = mod(iic2+jtp2+ML2,ML2)
          jc3 = mod(iic3+jtp3+ML3,ML3)
          jt  = KKK(jtp1,jtp2,jtp3)
          j   = LLL(jc1,jc2,jc3)

          call search_index( j, MLpot, Lpot(1,i), k )
          vloc_nl(k,i) = vloc_nl(k,i) + w(it,jt,ic)

       end do ! jtp1
       end do ! jtp2
       end do ! jtp3

       end do ! itp1
       end do ! itp2
       end do ! itp3

    end do ! ic1
    end do ! ic2
    end do ! ic3

    if ( disp_switch_parallel ) then
    call watch(ct1,et1) ; write(*,*) "(5)",ct1-ct0,et1-et0
    end if

! ---

    const=dble(ML)/dble(Ngrid_dense(1)*Ngrid_dense(2)*Ngrid_dense(3))

    vloc_nl(:,:) = vloc_nl(:,:)*const

    deallocate( w )
    deallocate( KKK )
    deallocate( LLL )

    if ( disp_switch_parallel ) then
    call watch(ct0,et0) ; write(*,*) "(6)",ct0-ct1,et0-et1
    end if

  END SUBROUTINE test2_localpot2

  SUBROUTINE search_index( inpt, n, lst, indx )
    implicit none
    integer,intent(IN)  :: inpt, n, lst(n)
    integer,intent(OUT) :: indx
    integer :: i,j,m0,m1,m
    m0=1
    m1=n
    do i=1,n
       if ( m1-m0 == 1 ) then
          if ( lst(m0) == inpt ) then
             indx = m0
             return
          else if ( lst(m1) == inpt ) then
             indx = m1
             return
          end if
          exit
       end if
       m = m0 + (m1-m0)/2
       j = lst(m)
       if ( j == inpt ) then
          indx = m
          return
       else if ( j < inpt ) then
          m0 = m
       else if ( j > inpt ) then
          m1 = m
       end if
    end do ! i
    do i=1,n
       if ( lst(i) == inpt ) then
          write(*,*) i,lst(i)
          exit
       end if
    end do
    write(*,*) "inpt,n=",inpt,n
    write(*,*) "m0,m1,m=",m0,m1,m
    stop "stop@search_index"
  END SUBROUTINE search_index

END MODULE localpot2_module
