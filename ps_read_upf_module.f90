!--------------------------------------------
! Unified Pseudopotential Format (UPF)
! This potential is adopted in QUANTUM ESPRESSO 
! The unit of energy is in Rydberg, and converted to Hartree in this routine
!--------------------------------------------
MODULE ps_read_UPF_module

  implicit none

  PRIVATE
  PUBLIC :: ps_read_UPF, ps_upf

  integer,parameter :: nrmax = 5000
  integer,parameter :: lmax  = 3

  TYPE upf
     integer :: nrr
     integer :: norb
     integer :: NRc(lmax+1)
     integer :: inorm(lmax+1)
     integer :: lo(lmax+1)
     integer :: no(lmax+1)
     real(8) :: znuc
     real(8) :: vps(nrmax,lmax+1)
     real(8) :: ups(nrmax,lmax+1)
     real(8) :: vql(nrmax)
     real(8) :: cdc(nrmax)
     real(8) :: cdd(nrmax)
     real(8) :: rr(nrmax)
     real(8) :: rx(nrmax)
     real(8) :: Rc(lmax+1)
     real(8) :: anorm(lmax+1)
     real(8) :: Dij(lmax+1,lmax+1)
  END TYPE upf

  type(upf) :: ps_upf

CONTAINS

  SUBROUTINE ps_read_UPF(g)
    implicit none
    integer,intent(IN) :: g
    integer,parameter :: max_loop = 100000
    integer :: loop
    character(30) :: cbuf

    ps_upf%nrr     =0
    ps_upf%norb    =0
    ps_upf%znuc    =0.0d0
    ps_upf%vps(:,:)=0.0d0
    ps_upf%ups(:,:)=0.0d0
    ps_upf%vql(:)  =0.0d0
    ps_upf%cdc(:)  =0.0d0
    ps_upf%cdd(:)  =0.0d0
    ps_upf%rr(:)   =0.0d0
    ps_upf%rx(:)   =0.0d0
    ps_upf%Rc(:)   =0.0d0
    ps_upf%NRc(:)  =0
    ps_upf%anorm(:)=0.0d0
    ps_upf%inorm(:)=0
    ps_upf%lo(:)   =0
    ps_upf%no(:)   =0
    ps_upf%Dij(:,:)=0.0d0

    do loop=1,max_loop

       read(g,'(a)',END=10) cbuf
       write(*,*) cbuf

       if ( cbuf(1:21) == '<UPF version="2.0.1">' ) then
          rewind g
          call ps_read_upf_ver201(g)
          return
       else if ( cbuf(1:9) == "<PP_INFO>" ) then
          rewind g
          call ps_read_upf_verorg(g)
          return
       end if

    end do ! loop

10  stop "Format is invalid (stop@ps_read_upf)"

  END SUBROUTINE ps_read_UPF


  SUBROUTINE ps_read_upf_verorg(g)
    implicit none
    integer,intent(IN) :: g
    integer,parameter :: max_loop = 100000
    integer :: i,j,n,i0,i1,loop,nrr,norb,nrc
    real(8) :: tmp
    character(30) :: cbuf

    write(*,'(a40," ps_read_upf_verorg")') repeat("-",40)

! Read

    do loop=1,max_loop

       read(g,*,END=10) cbuf

       if ( cbuf(1:11) == "<PP_HEADER>" ) then

          do i=1,5
             read(g,*)
          end do
          read(g,*) ps_upf%znuc

          do i=1,3
             read(g,*)
          end do
          read(g,*) nrr

       else if ( cbuf(1:11) == "<PP_HEADER " ) then

          do i=1,max_loop
             read(g,*) cbuf
             if ( cbuf(1:9) == "z_valence" ) then
                stop "test"
             end if
          end do

       end if

       if ( cbuf(1:9) == "<PP_MESH>" ) then

          do i=1,max_loop
             read(g,*) cbuf
             if ( cbuf(1:6) == "<PP_R>" ) exit
          end do

          read(g,*) ps_upf%rr(1:nrr)

          do i=1,max_loop
             read(g,*) cbuf
             if ( cbuf(1:8) == "<PP_RAB>" ) exit
          end do ! i

          read(g,*) ps_upf%rx(1:nrr)

       end if ! </PP_MESH>

       if ( cbuf(1:10) == "<PP_LOCAL>" ) then

          read(g,*) ps_upf%vql(1:nrr)

       end if ! </PP_LOCAL>

       if ( cbuf(1:13) == "<PP_NONLOCAL>" ) then

          norb=0
          do i=1,max_loop
             read(g,*) cbuf
             if ( cbuf(1:9) == "<PP_BETA>" ) then
                read(g,*) j, ps_upf%lo(j)
                read(g,*) nrc
                read(g,*) ps_upf%vps(1:nrc,j)
                read(g,*) cbuf
                ps_upf%NRc(j)=nrc
                norb=max( j, norb )
             else if ( cbuf(1:8) == "<PP_DIJ>" ) then
                read(g,*) n
                do j=1,n
                   read(g,*) i0,i1,ps_upf%anorm(j)
                end do
                do j=1,n
                   ps_upf%inorm(j)=nint( sign(1.0d0,ps_upf%anorm(j)) )
                   ps_upf%anorm(j)=abs( ps_upf%anorm(j) )
                end do
                read(g,*) cbuf
             else if ( cbuf(1:8) == "<" ) then
                write(*,*) cbuf,"exit"
                exit
             end if
          end do

       end if

       if ( cbuf(1:10) == "<PP_PSWFC>" ) then
       end if

       if ( cbuf(1:12) == "<PP_RHOATOM>" ) then

          read(g,*) ps_upf%cdd(1:nrr)

       end if

    end do ! loop
10 continue

    do i=nrr+1,2,-1
       ps_upf%rr(i)    = ps_upf%rr(i-1)
       ps_upf%rx(i)    = ps_upf%rx(i-1)
       ps_upf%cdd(i)   = ps_upf%cdd(i-1)
       ps_upf%vql(i)   = ps_upf%vql(i-1)
       ps_upf%vps(i,:) = ps_upf%vps(i-1,:)
       ps_upf%rr(i)    = ps_upf%rr(i-1)
    end do
    ps_upf%rr(1)    = 0.0d0
    ps_upf%rx(1)    = 0.0d0
    ps_upf%cdd(1)   = 0.0d0
    ps_upf%vql(1)   = ps_upf%vql(2)
    ps_upf%vps(1,:) = 0.0d0

    do j=1,norb
       do i=nrr,1,-1
          if ( abs(ps_upf%vps(i,j)) >= 1.d-13 ) then
             ps_upf%Rc(j) = ps_upf%rr(i-1)
             ps_upf%NRc(j)= i-1
             exit
          end if
       end do
    end do

    do j=1,norb
       tmp = sqrt( abs( ps_upf%anorm(j) ) )
       ps_upf%vps(:,j) = sqrt(0.5d0)*ps_upf%vps(:,j)*tmp
    end do

    ps_upf%nrr  = nrr + 1
    ps_upf%norb = norb

    ps_upf%vql(:)   = 0.5d0*ps_upf%vql(:)

    write(*,*) "*** Unified Pseudopotenetial Format (UPF) ***"
    write(*,*) "Znuc=",ps_upf%znuc
    write(*,*) "# of radial mesh points =",nrr
    write(*,*) "# of orbitals =",norb
    write(*,*) "angular momentum =",ps_upf%lo(1:norb)
    write(*,*) "cut off radius =",ps_upf%Rc(1:norb)
    write(*,*) "# of grid points within cut off radius",ps_upf%NRc(1:norb)
    write(*,*) "uVu integral (anorm) ="
    write(*,'(1x,8f10.5)') ( ps_upf%inorm(i)*ps_upf%anorm(i),i=1,norb )
    write(*,*) "sum(rhov)=",sum(ps_upf%cdd*ps_upf%rx)

    write(*,'(a40," ps_read_upf_verorg(end)")') repeat("-",40)

    return
  END SUBROUTINE ps_read_upf_verorg


  SUBROUTINE ps_read_upf_ver201(g)
    implicit none
    integer,intent(IN) :: g
    integer,parameter :: max_loop=1000000
    integer :: loop,i,j,k,l,no(0:5),lo
    character(100) :: cbuf, ckey
    integer :: norb,nrr,nsize
    real(8) :: tmp

    write(*,'(a40," ps_read_upf_ver201")') repeat("-",40)

    no(:) = 0

    do loop=1,max_loop

       read(g,'(a)',END=10) cbuf
       ckey = adjustl( cbuf )

       if ( ckey(1:6) == "</UPF>" ) then
          write(*,*) ckey(1:6)
          exit
       end if

       if ( ckey(1:11) == "<PP_HEADER " ) then

          write(*,*) ckey(1:11)

          backspace(g)

          do i=1,max_loop

             read(g,'(a)') cbuf

             j = index( cbuf, "z_valence=" )
             if ( j /= 0 ) then
                j = j + len_trim( "z_valence=" )
                k = index( cbuf(j+1:), '"' )
                ckey=cbuf(j+1:j+k-1)
                read(ckey,*) ps_upf%znuc
             end if

             j = index( cbuf, "mesh_size=" )
             if ( j /= 0 ) then
                j = j + len_trim( "mesh_size=" )
                k = index( cbuf(j+1:), '"' )
                ckey=cbuf(j+1:j+k-1)
                read(ckey,*) nrr
             end if

             j = index( cbuf, "number_of_proj=" )
             if ( j /= 0 ) then
                j = j + len_trim( "number_of_proj=" )
                k = index( cbuf(j+1:), '"' )
                ckey=cbuf(j+1:j+k-1)
                read(ckey,*) norb
             end if

             j = index( cbuf, "/>" )
             if ( j /= 0 ) exit

          end do ! i

       end if ! </PP_HEADER>

       if ( ckey(1:9) == "<PP_MESH " ) then

          write(*,*) ckey(1:9)

          do i=1,max_loop
             read(g,'(a)') cbuf
             ckey = adjustl( cbuf )
             if ( ckey(1:6) == "<PP_R " ) exit
          end do
          read(g,*) ps_upf%rr(1:nrr)

          do i=1,max_loop
             read(g,'(a)') cbuf
             ckey = adjustl( cbuf )
             if ( ckey(1:8) == "<PP_RAB " ) exit
          end do
          read(g,*) ps_upf%rx(1:nrr)

       end if ! </PP_MESH>

       if ( ckey(1:10) == "<PP_LOCAL " ) then

          write(*,*) ckey(1:10)

          read(g,*) ps_upf%vql(1:nrr)

       end if ! </PP_LOCAL>

       if ( ckey(1:13) == "<PP_NONLOCAL>" ) then

          write(*,*) ckey(1:13)

          j=0
          do i=1,max_loop

             read(g,'(a)') cbuf
             ckey = adjustl( cbuf )

             if ( ckey(1:9) == "<PP_BETA." ) then

                do k=1,max_loop

                   call get_num_from_string( cbuf, "angular_momentum=", lo )

                   l=index( cbuf, ">" )
                   if ( l /= 0 ) then
                      j=j+1
                      read(g,*) ps_upf%vps(1:nrr,j)
                      ps_upf%lo(j) = lo
                      no(lo) = no(lo) + 1
                      ps_upf%no(j) = no(lo)
                      exit
                   end if

                   read(g,'(a)') cbuf

                end do ! k

                if ( j == norb ) exit

             end if

          end do ! i

          do i=1,max_loop

             read(g,'(a)') cbuf
             ckey = adjustl( cbuf )

             if ( ckey(1:8) == "<PP_DIJ " ) then
                call get_num_from_string( cbuf, "size=", nsize )
                write(*,*) "nsize=",nsize
                read(g,*) ps_upf%Dij(1:norb,1:norb)
                j=count( ps_upf%Dij /= 0.0d0 )
                k=0
                do l=1,norb
                   if ( ps_upf%Dij(l,l) /= 0.0d0 ) k=k+1
                end do
                if ( j == k ) then
                   do l=1,norb
                      ps_upf%anorm(l) = ps_upf%Dij(l,l)
                   end do
                   ps_upf%Dij=0.0d0
                end if
                exit
             end if

          end do ! i

       end if ! </PP_NONLOCAL>

       if ( ckey(1:12) == "<PP_RHOATOM " ) then

          write(*,*) ckey(1:12)

          read(g,*) ps_upf%cdd(1:nrr)
          exit

       end if

    end do ! loop
10 continue

    do i=nrr+1,2,-1
       ps_upf%rr(i)    = ps_upf%rr(i-1)
       ps_upf%rx(i)    = ps_upf%rx(i-1)
       ps_upf%cdd(i)   = ps_upf%cdd(i-1)
       ps_upf%vql(i)   = ps_upf%vql(i-1)
       ps_upf%vps(i,:) = ps_upf%vps(i-1,:)
       ps_upf%rr(i)    = ps_upf%rr(i-1)
    end do
    ps_upf%rr(1)    = 0.0d0
    ps_upf%rx(1)    = 0.0d0
    ps_upf%cdd(1)   = 0.0d0
    ps_upf%vql(1)   = ps_upf%vql(2)
    ps_upf%vps(1,:) = 0.0d0

    do j=1,norb
       do i=nrr,1,-1
          if ( abs(ps_upf%vps(i,j)) >=1.d-13 ) then
             ps_upf%Rc(j) = ps_upf%rr(i-1)
             ps_upf%NRc(j)= i-1
             exit
          end if
       end do
    end do

    if ( any( ps_upf%anorm /= 0.0d0 ) ) then
       do j=1,norb
          ps_upf%inorm(j)=1
          if ( ps_upf%anorm(j) < 0.0d0 ) ps_upf%inorm(j)=-1
          tmp = sqrt( abs( ps_upf%anorm(j) ) )
          ps_upf%vps(:,j) = sqrt(0.5d0)*ps_upf%vps(:,j)*tmp
       end do
    else
       do j=1,norb
          ps_upf%vps(:,j) = sqrt(0.5d0)*ps_upf%vps(:,j)
       end do
    end if

    ps_upf%vql(:)   = 0.5d0*ps_upf%vql(:)

    ps_upf%nrr  = nrr + 1
    ps_upf%norb = norb

    write(*,*) "*** Unified Pseudopotenetial Format (UPF) ***"
    write(*,*) "Znuc=",ps_upf%znuc
    write(*,*) "# of radial mesh points =",nrr
    write(*,*) "# of orbitals =",norb
    write(*,*) "angular momentum =",ps_upf%lo(1:norb)
    write(*,*) "cut off radius =",ps_upf%Rc(1:norb)
    write(*,*) "# of grid points within cut off radius",ps_upf%NRc(1:norb)
    write(*,*) "uVu integral (anorm) ="
    write(*,'(1x,8f10.5)') ( ps_upf%inorm(i)*ps_upf%anorm(i),i=1,norb )
    write(*,*) "Dij ="
    write(*,'(1x,9f10.5)') (( ps_upf%Dij(i,j),i=1,norb ),j=1,norb)
    write(*,*) "sum(rhov)=",sum(ps_upf%cdd*ps_upf%rx)

    write(*,'(a40," ps_read_upf_ver201(end)")') repeat("-",40)

    return
  END SUBROUTINE ps_read_upf_ver201


  SUBROUTINE get_num_from_string( buf, str, n )
    implicit none
    character(*),intent(IN) :: buf, str
    integer,intent(OUT) :: n
    integer :: i, l, len_str
    character(8) :: str_out
    l = index( buf, str )
    if ( l /= 0 ) then
       len_str = len_trim( adjustl(str) )
       i=l+len_str+1
       str_out = buf(i:i)
       read(str_out,*) n
    end if
  END SUBROUTINE get_num_from_string


END MODULE ps_read_UPF_module
