MODULE hsort_module

  implicit none

  PRIVATE
  PUBLIC :: indexx

CONTAINS


  SUBROUTINE indexx( n, arr, indx )
    implicit none
    integer,intent(IN)  :: n
    real(8),intent(IN)  :: arr(n)
    integer,intent(OUT) :: indx(n)
    !call indexx_0( n, arr, indx )
    call indexx_1( n, arr, indx )
  END SUBROUTINE indexx


  SUBROUTINE indexx_0(n,arr,indx)
    implicit none
    integer,intent(IN)  :: n
    real(8),intent(IN)  :: arr(n)
    integer,intent(OUT) :: indx(n)
    integer,parameter :: m=7,nstack=50
    integer :: i,j,k,l,jstack,indxt,ir,itemp,istack(nstack)
    real(8) :: a

    do j=1,n
       indx(j)=j
    end do
    jstack=0
    l=1
    ir=n
1   if ( ir-l < m ) then
       do j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do i=j-1,1,-1
             if ( arr(indx(i))<=a ) goto 2
             indx(i+1)=indx(i)
          end do
          i=0
2         indx(i+1)=indxt
       end do
       if ( jstack==0 ) return
       ir=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
    else
       k=(l+ir)/2
       itemp=indx(k)
       indx(k)=indx(l+1)
       indx(l+1)=itemp
       if ( arr(indx(l+1)) > arr(indx(ir)) ) then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
       end if
       if ( arr(indx(l)) > arr(indx(ir)) ) then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
       end if
       if ( arr(indx(l+1)) > arr(indx(l)) ) then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
       end if
       i=l+1
       j=ir
       indxt=indx(l)
       a=arr(indxt)
3      continue
       i=i+1
       if ( arr(indx(i)) < a ) goto 3
4      continue
       j=j-1
       if ( arr(indx(j)) > a ) goto 4
       if ( j<i ) goto 5
       itemp=indx(i)
       indx(i)=indx(j)
       indx(j)=itemp
       goto 3
5      indx(l)=indx(j)
       indx(j)=indxt
       jstack=jstack+2
       if ( jstack>nstack ) stop "indexx"
       if ( ir-i+1>= j-l ) then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
       else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
       end if
    end if
    goto 1
  END SUBROUTINE indexx_0


  SUBROUTINE indexx_1( ndat, arry, indx )
    implicit none
    integer,intent(IN)  :: ndat
    real(8),intent(IN)  :: arry(ndat)
    integer,intent(OUT) :: indx(ndat)
    integer :: ilayer,iboss,jlayer,jboss,kboss
    integer :: i,i1,i2,max_layer,loop,ndat_max
    integer,allocatable :: irslt(:)
    logical :: flag_boss_change
    real(8) :: arry_max

    arry_max = maxval( arry )

    max_layer = log( dble(ndat+1) )/log( 2.0d0 )

    ndat_max = 2**max_layer - 1
    if ( ndat_max < ndat ) max_layer = max_layer + 1
    ndat_max = 2**max_layer - 1
    if ( ndat_max < ndat ) stop "stop(1)"

    do i=1,ndat
       indx(i)=i
    end do

! ---

    do ilayer=2,max_layer

       loop_boss : do iboss=2**(ilayer-2),2**(ilayer-1)-1

          if ( iboss > ndat ) exit loop_boss

          flag_boss_change = .false.

          do i=2*iboss,2*iboss+1

             if ( i > ndat ) exit loop_boss

             if ( arry(indx(i)) < arry(indx(iboss)) ) then
                call iswap( indx(iboss), indx(i) )
                flag_boss_change = .true.
             end if

             if ( flag_boss_change ) then
                kboss = iboss
                do jlayer=ilayer-1,2,-1
                   do jboss=2**(jlayer-2),2**(jlayer-1)-1
                      if ( arry(indx(kboss)) < arry(indx(jboss)) ) then
                         call iswap( indx(jboss), indx(kboss) )
                         kboss = jboss
                         exit
                      end if
                   end do ! jboss
                end do ! jlayer
                if ( arry(indx(kboss)) < arry(indx(1)) ) then
                   call iswap( indx(1), indx(kboss) )
                end if
             end if

          end do ! i

       end do loop_boss

    end do ! ilayer

! ---

    allocate( irslt(ndat) ) ; irslt=0

    do loop=1,ndat

       irslt(loop)=indx(1)
       indx(1)=-1

       jboss=1
       do ilayer=2,max_layer

          i1=2*jboss
          i2=2*jboss+1

          if ( i1 > ndat ) then
             exit
          else if ( i1 <= ndat .and. i2 > ndat ) then
             call iswap( indx(jboss), indx(i1) )
             exit
          end if

          if ( indx(i1) < 0 .and. indx(i2) < 0 ) then
             exit
          else if ( indx(i1) > 0 .and. indx(i2) < 0 ) then
             call iswap( indx(jboss), indx(i1) )
             jboss=i1
          else if ( indx(i1) < 0 .and. indx(i2) > 0 ) then
             call iswap( indx(jboss), indx(i2) )
             jboss=i2
          else if ( arry(indx(i1)) < arry(indx(i2)) ) then
             call iswap( indx(jboss), indx(i1) )
             jboss=i1
          else
             call iswap( indx(jboss), indx(i2) )
             jboss=i2
          end if

       end do ! ilayer

    end do ! loop

    indx(:)=irslt(:)

    deallocate( irslt )

  END SUBROUTINE indexx_1

  SUBROUTINE iswap( i, j )
    implicit none
    integer,intent(INOUT) :: i,j
    integer :: k
    k=i
    i=j
    j=k
  END SUBROUTINE iswap


END MODULE hsort_module
