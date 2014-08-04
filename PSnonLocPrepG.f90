MODULE PSnonLocPrepG
use parallel_module, only: myrank
  use VarPSMemberG
  use VarPSMember
  use atom_module, only: Natom,ki_atom
  use ps_nloc2_variables
  use electron_module, only: nspin
  ! nzlma

  implicit none
  
  PRIVATE
  PUBLIC :: prepNzqr

CONTAINS
  SUBROUTINE prepNzqr
    implicit none
    integer :: kk1
    integer :: i,j,ik
    integer :: lma1,lma2,a1,a2,l1,l2,m1,m2,i1,i2,a,l,m
    integer,allocatable :: k1a(:)
#ifdef _SHOWALL_INIT_
write(200+myrank,*) ">>>>> inside prepNzqr"
#endif

!----- get N_nzqr -----
    kk1=0
    do lma1=1,nzlma
      a1=amap(lma1)
      l1=lmap(lma1)
      m1=mmap(lma1)
      i1=no(iorbmap(lma1),ki_atom(a1))
      do lma2=1,nzlma
        a2=amap(lma2)
        l2=lmap(lma2)
        m2=mmap(lma2)
        i2=no(iorbmap(lma2),ki_atom(a2))
        if (a2/=a1) cycle
        if (l2>l1) cycle
        if (l1==l2 .and. i2>i1) cycle
        if (l1==l2 .and. i1==i2 .and. m2>m1) cycle
        kk1=kk1+1
      end do
    end do
    N_nzqr=kk1
    if (myrank==0) write(*,*) "N_nzqr= ",N_nzqr
!===== get N_nzqr =====

!----- get N_nlop -----
    kk1=0
    do lma1=1,nzlma
      do lma2=1,nzlma
        if (amap(lma1)/=amap(lma2) .or. lmap(lma1)/=lmap(lma2) .or. mmap(lma1)/=mmap(lma2)) cycle
        kk1=kk1+1
      end do
    end do
    N_nlop=kk1
    if (myrank==0) write(*,*) "N_nlop= ",N_nlop
!===== get N_nlop =====

    call allocateNzqr(nspin,Natom,k1max)

!----- get nzqr_pair, atommap, k1map, kk1map -----
    allocate( k1a(Natom) ) ; k1a(:)=0
    kk1=0
    do lma1=1,nzlma
      a1=amap(lma1)
      l1=lmap(lma1)
      m1=mmap(lma1)
      i1=no(iorbmap(lma1),ki_atom(a1))
      do lma2=1,nzlma
        a2=amap(lma2)
        l2=lmap(lma2)
        m2=mmap(lma2)
        i2=no(iorbmap(lma2),ki_atom(a2))
        if (a2/=a1) cycle
        if (l2>l1) cycle
        if (l1==l2 .and. i2>i1) cycle
        if (l1==l2 .and. i1==i2 .and. m2>m1) cycle
        k1a(a1)=k1a(a1)+1
        kk1=kk1+1
        if (kk1>N_nzqr) then
          write(*,*) 'myrank: kk1>N_nzqr= ',myrank,kk1,N_nzqr
          stop 'inconsistancy in kk1 and N_nzqr'
        end if
        nzqr_pair(kk1,1)  =lma1
        nzqr_pair(kk1,2)  =lma2
        atommap(kk1)      =a1
        k1map(kk1)        =k1a(a1)
        kk1map(k1a(a1),a1)=kk1
      end do
    end do
    deallocate( k1a )
!===== get nzqr_pair, atommap, k1map, kk1map =====

    if (myrank==0) write(*,*) "--- Dij00(1:N_nzqr) ---"
    if (myrank==0) write(*,*) "  [ qij_f(1:N_nzqr) ]  "
    do kk1=1,N_nzqr
      i=nzqr_pair(kk1,1)
      j=nzqr_pair(kk1,2)
      a1=amap(i)
      a2=amap(j)
      l1=lmap(i)
      l2=lmap(j)
      m1=mmap(i)
      m2=mmap(j)
      if (.not.(a1==a2 .and. l1==l2 .and. m1==m2)) cycle
      ik=ki_atom(a1)
      i1=no(iorbmap(i),ik)
      i2=no(iorbmap(j),ik)
      Dij00(kk1)=ddi(i1,i2,l1+1,ik)
      qij_f(kk1)=qqc(i1,i2,l1+1,ik)
    end do
#ifdef _SHOWALL_Q_
    do kk1=1,N_nzqr
      i=nzqr_pair(kk1,1)
      j=nzqr_pair(kk1,2)
      a1=amap(i)
      a2=amap(j)
      l1=lmap(i)
      l2=lmap(j)
      m1=mmap(i)
      m2=mmap(j)
      ik=ki_atom(a1)
      i1=no(iorbmap(i),ik)
      i2=no(iorbmap(j),ik)
      if (myrank==0) then
        if (i==j) then
          write(540,'(1x,i6,2i3,2x,4i4,2x,4i4,f15.10)') kk1,i,j,a1,l1,m1,i1,a2,l2,m2,i2,Dij00(kk1)
        else
          write(540,'(1x,i6,2i3,2x,4i4,2x,4i4,f15.10)') kk1,i,j,a1,l1,m1,i1,a2,l2,m2,i2,Dij00(kk1)
          write(540,'(1x,i6,2i3,2x,4i4,2x,4i4,f15.10,"(oginau)")') kk1,i,j,a1,l1,m1,i1,a2,l2,m2,i2,Dij00(kk1)
        end if
      end if
    end do
#endif

!----- Nlop_type Matrix -----
    if (myrank==0) write(*,*) "--- Dij0 (1:N_nlop) ---"
    if (myrank==0) write(*,*) "   [ qij(1:N_nlop) ]   "
    kk1=0
    do lma1=1,nzlma
      a=amap(lma1)
      l=lmap(lma1)
      m=mmap(lma1)
      ik=ki_atom(a)
      do lma2=1,nzlma
        if (a/=amap(lma2) .or. l/=lmap(lma2) .or. m/=mmap(lma2)) cycle
        kk1=kk1+1
        nlop_pair(1,kk1)=lma1
        nlop_pair(2,kk1)=lma2
        i1=no(iorbmap(lma1),ik)
        i2=no(iorbmap(lma2),ik)
        Dij0(kk1)=ddi(i1,i2,l+1,ik)
#ifdef _SHOWALL_Q_
        if (myrank==0) then
          write(540,*) kk1,lma1,lma2,Dij0(kk1)
        end if
#endif
        qij(kk1)=qqc(i1,i2,l+1,ik)
      end do
    end do
    if (myrank==0) then
      write(*,*) "N_nlop= ",N_nlop,kk1
    end if
!===== Nlop_type Matrix =====
#ifdef _SHOWALL_INIT_
write(200+myrank,*) "<<<<< end of prepNzqr"
#endif

    return
  END SUBROUTINE prepNzqr
  

END MODULE PSnonLocPrepG
