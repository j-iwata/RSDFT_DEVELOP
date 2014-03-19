MODULE bc_test_module

  use rgrid_module, only: Ngrid,Igrid
  use parallel_module
  use bc_variables

  implicit none

  PRIVATE
  PUBLIC :: bcset_test1

CONTAINS


  SUBROUTINE bcset_test1(ib1,ib2,ndepth,idir)
    implicit none
    integer,intent(IN) :: ib1,ib2,ndepth,idir
    integer :: a1,a2,a3,b1,b2,b3,nb,ns,ms,mt
    integer :: m,n,ndata,i1,i2,i3,ib,ierr
    integer :: c1,c2,c3,d1,d2,d3,irank,nreq,itags,itagr,ireq(36)
    integer :: i,j,l
    integer :: istatus(mpi_status_size,123)

    a1=Igrid(1,1)
    b1=Igrid(2,1)
    a2=Igrid(1,2)
    b2=Igrid(2,2)
    a3=Igrid(1,3)
    b3=Igrid(2,3)

    nb=ib2-ib1+1

!$OMP master
    nreq=0
!$OMP end master

!(1) [sender][2]-->[1][receiver]

    if ( idir==0 .or. idir==1 .or. idir==2 ) then
       do n=1,n_neighbor(1)
          if ( ndepth==1 ) then
             ndata = fdinfo_recv(10,n,1)*nb
          else
             ndata = fdinfo_recv(9,n,1)*nb
          end if
          if ( ndata<1 ) cycle
          irank = fdinfo_recv(8,n,1)
          itagr = 10
!$OMP barrier
!$OMP master
          nreq  = nreq + 1
          call mpi_irecv(rbuf(1,n,1),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master
!$OMP barrier
       end do
       do n=1,n_neighbor(2)
          if ( ndepth==1 ) then
             c1=b1
             d1=b1
             j=fdinfo_send(10,n,2)*nb
          else
             c1=fdinfo_send(1,n,2)
             d1=fdinfo_send(2,n,2)
             j=fdinfo_send(9,n,2)*nb
          end if
          c2=fdinfo_send(3,n,2) ; d2=fdinfo_send(4,n,2)
          c3=fdinfo_send(5,n,2) ; d3=fdinfo_send(6,n,2)
          if ( j<1 ) cycle
!$OMP do
          do ib=ib1,ib2
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                ndata = 1 + i1-c1 + (i2-c2)*(d1-c1+1) &
                     + (i3-c3)*(d1-c1+1)*(d2-c2+1) &
                     + (ib-ib1)*(d1-c1+1)*(d2-c2+1)*(d3-c3+1) 
                sbuf(ndata,n,2)=www(i1,i2,i3,ib)
             end do
             end do
             end do
          end do
!$OMP end do
          !if ( ndata/=j ) stop
          irank = fdinfo_send(7,n,2)
          itags = 10
!$OMP barrier
!$OMP master
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,2),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master
!$OMP barrier
       end do
    end if

!    call mpi_waitall(nreq,ireq,istatus,ierr) ; nreq=0

!(3) [sender][4]-->[3][receiver]

    if ( idir==0 .or. idir==3 .or. idir==4 ) then
       do n=1,n_neighbor(3)
          if ( ndepth==1 ) then
             ndata = fdinfo_recv(10,n,3)*nb
          else
             ndata = fdinfo_recv(9,n,3)*nb
          end if
          if ( ndata<1 ) cycle
          irank = fdinfo_recv(8,n,3)
          itagr = 30
!$OMP barrier
!$OMP master
          nreq  = nreq + 1
          call mpi_irecv(rbuf(1,n,3),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master
!$OMP barrier
       end do
       do n=1,n_neighbor(4)
          if ( ndepth==1 ) then
             c2=b2
             d2=b2
             j=fdinfo_send(10,n,4)*nb
          else
             c2=fdinfo_send(3,n,4)
             d2=fdinfo_send(4,n,4)
             j=fdinfo_send(9,n,4)*nb
          end if
          c1=fdinfo_send(1,n,4) ; d1=fdinfo_send(2,n,4)
          c3=fdinfo_send(5,n,4) ; d3=fdinfo_send(6,n,4)
          if ( j<1 ) cycle
!$OMP do
          do ib=ib1,ib2
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                ndata = 1 + i1-c1 + (i2-c2)*(d1-c1+1) &
                     + (i3-c3)*(d1-c1+1)*(d2-c2+1) &
                     + (ib-ib1)*(d1-c1+1)*(d2-c2+1)*(d3-c3+1) 
                sbuf(ndata,n,4)=www(i1,i2,i3,ib)
             end do
             end do
             end do
          end do
!$OMP end do
          !if ( ndata/=j ) stop
          irank = fdinfo_send(7,n,4)
          itags = 30
!$OMP barrier
!$OMP master
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,4),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master
!$OMP barrier
       end do
    end if

!    call mpi_waitall(nreq,ireq,istatus,ierr) ; nreq=0

!(5) [sender][6]-->[5][receiver]

    if ( idir==0 .or. idir==5 .or. idir==6 ) then
       do n=1,n_neighbor(5)
          if ( ndepth==1 ) then
             ndata = fdinfo_recv(10,n,5)*nb
          else
             ndata = fdinfo_recv(9,n,5)*nb
          end if
          if ( ndata<1 ) cycle
          irank = fdinfo_recv(8,n,5)
          itagr = 50
!$OMP barrier
!$OMP master
          nreq  = nreq + 1
          call mpi_irecv(rbuf(1,n,5),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master
!$OMP barrier
       end do
       do n=1,n_neighbor(6)
          if ( ndepth==1 ) then
             c3=b3
             d3=b3
             j=fdinfo_send(10,n,6)*nb
          else
             c3=fdinfo_send(5,n,6)
             d3=fdinfo_send(6,n,6)
             j=fdinfo_send(9,n,6)*nb
          end if
          c1=fdinfo_send(1,n,6) ; d1=fdinfo_send(2,n,6)
          c2=fdinfo_send(3,n,6) ; d2=fdinfo_send(4,n,6)
          if ( j<1 ) cycle
!$OMP do
          do ib=ib1,ib2
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                ndata = 1 + i1-c1 + (i2-c2)*(d1-c1+1) &
                     + (i3-c3)*(d1-c1+1)*(d2-c2+1) &
                     + (ib-ib1)*(d1-c1+1)*(d2-c2+1)*(d3-c3+1) 
                sbuf(ndata,n,6)=www(i1,i2,i3,ib)
             end do
             end do
             end do
          end do
!$OMP end do
          !if ( ndata/=j ) stop
          irank = fdinfo_send(7,n,6)
          itags = 50
!$OMP barrier
!$OMP master
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,6),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master
!$OMP barrier
       end do
    end if

!    call mpi_waitall(nreq,ireq,istatus,ierr) ; nreq=0

!(2) [receiver][2]<--[1][sender]

    if ( idir==0 .or. idir==1 .or. idir==2 ) then
       do n=1,n_neighbor(2)
          if ( ndepth==1 ) then
             ndata = fdinfo_recv(10,n,2)*nb
          else
             ndata = fdinfo_recv(9,n,2)*nb
          end if
          if ( ndata<1 ) cycle
          irank = fdinfo_recv(8,n,2)
          itagr = 20
!$OMP barrier
!$OMP master
          nreq  = nreq + 1
          call mpi_irecv(rbuf(1,n,2),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master
!$OMP barrier
       end do
       do n=1,n_neighbor(1)
          if ( ndepth==1 ) then
             c1=a1
             d1=a1
             j=fdinfo_send(10,n,1)*nb
          else
             c1=fdinfo_send(1,n,1)
             d1=fdinfo_send(2,n,1)
             j=fdinfo_send(9,n,1)*nb
          end if
          c2=fdinfo_send(3,n,1) ; d2=fdinfo_send(4,n,1)
          c3=fdinfo_send(5,n,1) ; d3=fdinfo_send(6,n,1)
          if ( j<1 ) cycle
!$OMP do
          do ib=ib1,ib2
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                ndata = 1 + i1-c1 + (i2-c2)*(d1-c1+1) &
                     + (i3-c3)*(d1-c1+1)*(d2-c2+1) &
                     + (ib-ib1)*(d1-c1+1)*(d2-c2+1)*(d3-c3+1) 
                sbuf(ndata,n,1)=www(i1,i2,i3,ib)
             end do
             end do
             end do
          end do
!$OMP end do
          !if ( ndata/=j ) stop
          irank = fdinfo_send(7,n,1)
          itags = 20
!$OMP barrier
!$OMP master
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,1),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master
!$OMP barrier
       end do
    end if

!    call mpi_waitall(nreq,ireq,istatus,ierr) ; nreq=0

!(4) [receiver][4]<--[3][sender]

    if ( idir==0 .or. idir==3 .or. idir==4 ) then
       do n=1,n_neighbor(4)
          if ( ndepth==1 ) then
             ndata = fdinfo_recv(10,n,4)*nb
          else
             ndata = fdinfo_recv(9,n,4)*nb
          end if
          if ( ndata<1 ) cycle
          irank = fdinfo_recv(8,n,4)
          itagr = 40
!$OMP barrier
!$OMP master
          nreq  = nreq + 1
          call mpi_irecv(rbuf(1,n,4),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master
!$OMP barrier
       end do
       do n=1,n_neighbor(3)
          if ( ndepth==1 ) then
             c2=a2
             d2=a2
             j=fdinfo_send(10,n,3)*nb
          else
             c2=fdinfo_send(3,n,3)
             d2=fdinfo_send(4,n,3)
             j=fdinfo_send(9,n,3)*nb
          end if
          c1=fdinfo_send(1,n,3) ; d1=fdinfo_send(2,n,3)
          c3=fdinfo_send(5,n,3) ; d3=fdinfo_send(6,n,3)
          if ( j<1 ) cycle
!$OMP do
          do ib=ib1,ib2
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                ndata = 1 + i1-c1 + (i2-c2)*(d1-c1+1) &
                     + (i3-c3)*(d1-c1+1)*(d2-c2+1) &
                     + (ib-ib1)*(d1-c1+1)*(d2-c2+1)*(d3-c3+1) 
                sbuf(ndata,n,3)=www(i1,i2,i3,ib)
             end do
             end do
             end do
          end do
!$OMP end do
          !if ( ndata/=j ) stop
          irank = fdinfo_send(7,n,3)
          itags = 40
!$OMP barrier
!$OMP master
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,3),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master
!$OMP barrier
       end do
    end if

!    call mpi_waitall(nreq,ireq,istatus,ierr) ; nreq=0

!(6) [receiver][6]<--[5][sender]

    if ( idir==0 .or. idir==5 .or. idir==6 ) then
       do n=1,n_neighbor(6)
          if ( ndepth==1 ) then
             ndata = fdinfo_recv(10,n,6)*nb
          else
             ndata = fdinfo_recv(9,n,6)*nb
          end if
          if ( ndata<1 ) cycle
          irank = fdinfo_recv(8,n,6)
          itagr = 60
!$OMP barrier
!$OMP master
          nreq  = nreq + 1
          call mpi_irecv(rbuf(1,n,6),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master
!$OMP barrier
       end do
       do n=1,n_neighbor(5)
          if ( ndepth==1 ) then
             c3=a3
             d3=a3
             j=fdinfo_send(10,n,5)*nb
          else
             c3=fdinfo_send(5,n,5)
             d3=fdinfo_send(6,n,5)
             j=fdinfo_send(9,n,5)*nb
          end if
          c1=fdinfo_send(1,n,5) ; d1=fdinfo_send(2,n,5)
          c2=fdinfo_send(3,n,5) ; d2=fdinfo_send(4,n,5)
          if ( j<1 ) cycle
!$OMP do
          do ib=ib1,ib2
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                ndata = 1 + i1-c1 + (i2-c2)*(d1-c1+1) &
                     + (i3-c3)*(d1-c1+1)*(d2-c2+1) &
                     + (ib-ib1)*(d1-c1+1)*(d2-c2+1)*(d3-c3+1) 
                sbuf(ndata,n,5)=www(i1,i2,i3,ib)
             end do
             end do
             end do
          end do
!$OMP end do
          !if ( ndata/=j ) stop
          irank = fdinfo_send(7,n,5)
          itags = 60
!$OMP barrier
!$OMP master
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,5),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master
!$OMP barrier
       end do
    end if

!$OMP barrier
!$OMP master
    call mpi_waitall(nreq,ireq,istatus,ierr) ; nreq=0
!$OMP end master
!$OMP barrier

    do m=1,6
       do n=1,n_neighbor(m)

          c1=fdinfo_recv(1,n,m)
          d1=fdinfo_recv(2,n,m)
          c2=fdinfo_recv(3,n,m)
          d2=fdinfo_recv(4,n,m)
          c3=fdinfo_recv(5,n,m)
          d3=fdinfo_recv(6,n,m)

          if ( Md>ndepth .and. ndepth==1 ) then
             if ( fdinfo_recv(10,n,m)<1 ) cycle
             select case(m)
             case(1)
                c1=a1-1
                d1=c1
             case(2)
                c1=b1+1
                d1=c1
             case(3)
                c2=a2-1
                d2=c2
             case(4)
                c2=b2+1
                d2=c2
             case(5)
                c3=a3-1
                d3=c3
             case(6)
                c3=b3+1
                d3=c3
             end select
          else
             if ( fdinfo_recv(9,n,m)<1 ) cycle
          end if

!$OMP do
          do ib=ib1,ib2
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                i = 1 + i1-c1 + (i2-c2)*(d1-c1+1) &
                      + (i3-c3)*(d1-c1+1)*(d2-c2+1) &
                      + (ib-ib1)*(d1-c1+1)*(d2-c2+1)*(d3-c3+1) 
                www(i1,i2,i3,ib)=rbuf(i,n,m)
             end do
             end do
             end do
          end do
!$OMP end do

       end do ! n
    end do ! m

    return

  END SUBROUTINE bcset_test1


END MODULE bc_test_module
