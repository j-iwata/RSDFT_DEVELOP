MODULE ps_nloc1_module

  use pseudopot_module

  implicit none

  PRIVATE
  PUBLIC :: Ndense,Nintp,read_ps_nloc1,read_oldformat_ps_nloc1

  integer :: Ndense, Nintp
  real(8),allocatable :: rad1(:,:)

  INTERFACE
     FUNCTION Ylm(x,y,z,l,m)
       real(8) :: Ylm
       real(8),intent(IN) :: x,y,z
       integer,intent(IN) :: l,m
     END FUNCTION Ylm
  END INTERFACE

CONTAINS


  SUBROUTINE read_ps_nloc1(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(3) :: cbuf,ckey
    Ndense=1
    Nintp =0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:3) == "NL1" ) then
             backspace(unit)
             read(unit,*) cbuf,Ndense,Nintp
          end if
       end do
999    continue
       write(*,*) "Ndense, Nintp =",Ndense,Nintp
    end if
    call send_ps_nloc1(0)
  END SUBROUTINE read_ps_nloc1


  SUBROUTINE read_oldformat_ps_nloc1(rank,unit)
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) Ndense, Nintp
       write(*,*) "Ndense, Nintp =",Ndense,Nintp
    end if
    call send_ps_nloc1(0)
  END SUBROUTINE read_oldformat_ps_nloc1


  SUBROUTINE send_ps_nloc1(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(Ndense ,1,mpi_integer,rank,mpi_comm_world,ierr)
    call mpi_bcast(Nintp  ,1,mpi_integer,rank,mpi_comm_world,ierr)
  END SUBROUTINE send_ps_nloc1


#ifdef TEST
  SUBROUTINE init_pselect_1(MKI)
    integer,intent(IN) :: MKI
    integer :: m
    m=maxval( Mr )
    allocate( rad1(m,MKI) )
    rad1(:,:)=rad(:,:)
  END SUBROUTINE init_pselect_1

  SUBROUTINE prep_pselect_1
    complex(8) :: ztmp0
    integer,allocatable :: icheck_tmp1(:),icheck_tmp2(:),itmp(:,:)
    integer,allocatable :: icheck_tmp3(:,:,:),icheck_tmp4(:,:,:)
    integer,allocatable :: sendmap0(:,:),recvmap0(:,:),ireq(:)
    integer,allocatable :: lma_nsend_tmp(:),maps_tmp(:,:),itmp1(:)
    integer,allocatable :: irad(:,:),nl_rank_map_tmp(:),itmp3(:,:)
    integer,allocatable :: JJ_tmp(:,:),JJ_tmp_0(:,:),itmp2(:)
    integer,allocatable :: jtmp3(:,:,:),mtmp3(:),istatus(:,:)
    integer :: a,i,j,k,L,m,n,mm1,mm2,mm3,m1,m2,ML0,k1,k2,k3
    integer :: i1,i2,i3,j1,j2,j3,ik,ir,m0,iorb,mm,ierr,ir0,irlma
    integer :: ic1,ic2,ic3,id1,id2,id3,ii1,ii2,ii3,iii1,iii2,iii3
    integer :: Nintp_0,nzlma_0,M_irad,NRc,MMJJ_0,lma,lma0,i1_0,i2_0,i3_0
    integer :: nreq,ibuf(3,3),irank
    real(8),parameter :: ep=1.d-8
    real(8) :: x,y,z,r,Rx,Ry,Rz,Rps2,v,v0,d1,d2,d3,r2,kr,pi2
    real(8) :: tmp0,tmp1,tmp2,tmp3,c1,c2,c3,maxerr,err0,err
    real(8) :: ct(10),et(10),ct0,et0,ct1,et1,ct00,et00
    real(8) :: ctime0,ctime1,etime0,etime1,mem,memax
    real(8),allocatable :: uVtmp(:),work(:)
    real(8),allocatable :: rtmp3(:,:),sss(:,:,:,:),rrr(:,:,:,:)


    if ( Ndense <= 0 .or. Nintp <= 0 ) then
       Ndense=1
       Nintp =0
    end if
    if ( Nintp > Md ) then
       Nintp = Md
    end if
    H1d  = H1/Ndense
    H2d  = H2/Ndense
    H3d  = H3/Ndense
    ML1d = ML1*Ndense
    ML2d = ML2*Ndense
    ML3d = ML3*Ndense
    MLd  = ML1d*ML2d*ML3d
    dVd  = (H1d*H2d*H3d)/(H1*H2*H3)
    Nintp_0 = min(0,-Nintp+1)


    allocate( Clag1(Nintp_0:Nintp,0:Ndense-1) ) ; Clag1=0.d0
    allocate( Clag2(Nintp_0:Nintp,0:Ndense-1) ) ; Clag2=0.d0
    allocate( Clag3(Nintp_0:Nintp,0:Ndense-1) ) ; Clag3=0.d0
    do j1=0,Ndense-1
       do i1=Nintp_0,Nintp
          Clag1(i1,j1)=1.d0
          do i2=Nintp_0,Nintp
             if ( i2==i1 ) cycle
             x=dble(j1-i2*Ndense)
             y=dble(i1-i2)*Ndense
             Clag1(i1,j1)=Clag1(i1,j1)*(x/y)
          end do
       end do
    end do
    Clag2(:,:)=Clag1(:,:)
    Clag3(:,:)=Clag1(:,:)


    if ( Mlma < nprocs_g ) then
       nzlma_0 = Mlma
    else
       nzlma_0 = min(Mlma*125/nprocs_g,Mlma)
    end if


    r=maxval(Rps)
    call Make_MinimalBox(r,mm1,mm2,mm3,MMJJ_0)
    mm1_ps = maxval( abs(mcube_grid_ion(:,1)) )
    mm2_ps = maxval( abs(mcube_grid_ion(:,2)) )
    mm3_ps = maxval( abs(mcube_grid_ion(:,3)) )
    mm1 = mm1_ps+Nintp+1
    mm2 = mm2_ps+Nintp+1
    mm3 = mm3_ps+Nintp+1

    allocate( icheck_tmp3(-mm1:mm1,-mm2:mm2,-mm3:mm3) ) ; icheck_tmp3=0
    do i=1,M_grid_ion
       i1=map_grid_ion(1,i)
       i2=map_grid_ion(2,i)
       i3=map_grid_ion(3,i)
       do j3=Nintp_0-1,Nintp+1
       do j2=Nintp_0-1,Nintp+1
       do j1=Nintp_0-1,Nintp+1
          icheck_tmp3(i1+j1,i2+j2,i3+j3)=1
       end do
       end do
       end do
    end do
    MMJJ_0=sum(icheck_tmp3)
    deallocate( icheck_tmp3 )


    L=maxval(lo)
    n=maxval(norb)
    allocate( icheck_tmp3(MI,n,2*L+1)               ) ; icheck_tmp3=0
    allocate( icheck_tmp1(0:nprocs_g-1)             ) ; icheck_tmp1=0
    allocate( JJ_tmp_0(0:9,MMJJ_0)                  ) ; JJ_tmp_0=0
    allocate( rtmp3(MMJJ_0,nzlma_0)                 ) ; rtmp3=0.d0
    allocate( sss(-mm1:mm1,-mm2:mm2,-mm3:mm3,2*L+1) ) ; sss=0.d0
    allocate( rrr(-mm1:mm1,-mm2:mm2,-mm3:mm3,2*L+1) ) ; rrr=0.d0
    allocate( uVtmp(-L:L)                           ) ; uVtmp=0.d0
    allocate( jtmp3(6,MMJJ_0,nzlma_0)               ) ; jtmp3=0
    allocate( mtmp3(nzlma_0)                        ) ; mtmp3=0


    allocate( irad(0:3000,MKI) ) ; irad=0
    M_irad=0
    do ik=1,MKI
       NRc=maxval( NRps(:,ik) )
       NRc=min( 3000, NRc )
       m=0
       irad(0,ik)=1
       do ir=1,NRc
          m=int(100.d0*rad1(ir,ik))+1
          irad( m,ik )=ir
       end do
       ir=irad(0,ik)
       do i=1,m
          if ( irad(i,ik)==0 ) then
             irad(i,ik)=ir
             cycle
          end if
          ir=irad(i,ik)
       end do
       irad(m+1:,ik)=ir
       M_irad=max(M_irad,m)
    end do


    c1                 = 1.d0/ML1
    c2                 = 1.d0/ML2
    c3                 = 1.d0/ML3
    maxerr             = 0
    lma                = 0
    icheck_tmp3(:,:,:) = 0
    MMJJ               = 0
    MMJJ_t             = 0
    nzlma              = 0
    lma0               = 0
    m0                 = size(rrr)

    do a=1,MI

       Rx = aa(1,1)*asi(1,a)+aa(1,2)*asi(2,a)+aa(1,3)*asi(3,a)
       Ry = aa(2,1)*asi(1,a)+aa(2,2)*asi(2,a)+aa(2,3)*asi(3,a)
       Rz = aa(3,1)*asi(1,a)+aa(3,2)*asi(2,a)+aa(3,3)*asi(3,a)

       ic1 = nint( asi(1,a)*ML1 )
       ic2 = nint( asi(2,a)*ML2 )
       ic3 = nint( asi(3,a)*ML3 )

       ik = Kion(a)

       do iorb=1,norb(ik)

          Rps2           = Rps(iorb,ik)**2
          NRc            = NRps(iorb,ik)
          L              = lo(iorb,ik)
          sss(:,:,:,:)   = 0.d0
          rrr(:,:,:,:)   = 0.d0
          icheck_tmp1(:) = 0

          do i=1,M_grid_ion

             i1 = map_grid_ion(1,i)
             i2 = map_grid_ion(2,i)
             i3 = map_grid_ion(3,i)

             id1 = ic1 + i1
             id2 = ic2 + i2
             id3 = ic3 + i3

             do ii3=Nintp_0,Nintp
                iii3=id3+ii3
                k3=iii3/ML3 ; if ( iii3<0 ) k3=(iii3+1)/ML3-1
                i3_0=iii3-k3*ML3
             do ii2=Nintp_0,Nintp
                iii2=id2+ii2
                k2=iii2/ML2 ; if ( iii2<0 ) k2=(iii2+1)/ML2-1
                i2_0=iii2-k2*ML2
             do ii1=Nintp_0,Nintp
                iii1=id1+ii1
                k1=iii1/ML1 ; if ( iii1<0 ) k1=(iii1+1)/ML1-1
                i1_0=iii1-k1*ML1
                irank=PPP2(i1_0,i2_0,i3_0)
                icheck_tmp1(irank)=icheck_tmp1(irank)+1
             end do
             end do
             end do

             k1=id1/ML1 ; if ( id1<0 ) k1=(id1+1)/ML1-1
             k2=id2/ML2 ; if ( id2<0 ) k2=(id2+1)/ML2-1
             k3=id3/ML3 ; if ( id3<0 ) k3=(id3+1)/ML3-1
             i1_0=id1-k1*ML1
             i2_0=id2-k2*ML2
             i3_0=id3-k3*ML3

             if ( PPP2(i1_0,i2_0,i3_0)==myrank_g ) then

                do j3=0,Ndense-1
                   d3=id3*H3+j3*H3d
                do j2=0,Ndense-1
                   d2=id2*H2+j2*H2d
                do j1=0,Ndense-1
                   d1=id1*H1+j1*H1d

                   x  = ee(1,1)*d1+ee(1,2)*d2+ee(1,3)*d3-Rx
                   y  = ee(2,1)*d1+ee(2,2)*d2+ee(2,3)*d3-Ry
                   z  = ee(3,1)*d1+ee(3,2)*d2+ee(3,3)*d3-Rz
                   r2 = x*x+y*y+z*z

                   if ( r2>Rps2+1.d-10 ) cycle

                   r        = sqrt(r2)
                   v0       = 0.d0
                   uVtmp(:) = 0.d0
                   err0     = 0.d0

                   if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 .or. abs(z)>1.d-14 .or. L==0 ) then
                      ir0=irad( int(100.d0*r),ik )
                      do ir=ir0,NRc
                         if ( r<rad1(ir,ik) ) exit
                      end do
                      if ( ir<=2 ) then
                         v0=viod(2,iorb,ik)/rad1(2,ik)
                      else if ( ir<NRc ) then
                         err0=1.d10
                         do mm=1,20
                            m1=max(1,ir-mm)
                            m2=min(ir+mm,NRc)
                            call polint(rad1(m1,ik),viod(m1,iorb,ik),m2-m1+1,r,v,err)
                            if ( abs(err)<err0 ) then
                               v0=v
                               err0=abs(err)
                               if ( err0<ep ) exit
                            end if
                         end do
                         v0=v0/r
                      end if
                      maxerr=max(maxerr,err0)
                      do m=-L,L
                         uVtmp(m)=v0*Ylm(x,y,z,L,m)
                      end do
                   end if

                   if ( v0==0.d0 ) cycle

                   do m=1,2*L+1
                      tmp0=uVtmp(m-1-L)*dVd
                      if ( abs(tmp0)<1.d-10 ) cycle
                      do ii3=Nintp_0,Nintp
                         iii3=i3+ii3
                         tmp3=Clag3(ii3,j3)*tmp0
                      do ii2=Nintp_0,Nintp
                         iii2=i2+ii2
                         tmp2=tmp3*Clag2(ii2,j2)
                      do ii1=Nintp_0,Nintp
                         iii1=i1+ii1
                         sss(iii1,iii2,iii3,m)=sss(iii1,iii2,iii3,m)+tmp2*Clag1(ii1,j1)
                      end do
                      end do
                      end do
                   end do

                end do ! j3
                end do ! j2
                end do ! j1

             end if ! PPP2

          end do ! i ( 1 - M_grid_ion )

          if ( icheck_tmp1(myrank_g)>0 ) then

             call mpi_allreduce_nlpp( sss,rrr,m0,icheck_tmp1,ierr,"R" )

             j1=0
             j2=0
             do i3=-mm3,mm3
                id3=i3+ic3
                k3=id3/ML3 ; if (id3<0) k3=(id3+1)/ML3-1
                ii3=id3-k3*ML3
                if ( ii3<a3b .or. b3b<ii3 ) cycle
             do i2=-mm2,mm2
                id2=i2+ic2
                k2=id2/ML2 ; if (id2<0) k2=(id2+1)/ML2-1
                ii2=id2-k2*ML2
                if ( ii2<a2b .or. b2b<ii2 ) cycle
             do i1=-mm1,mm1
                id1=i1+ic1
                k1=id1/ML1 ; if (id1<0) k1=(id1+1)/ML1-1
                ii1=id1-k1*ML1
                if ( ii1<a1b .or. b1b<ii1 ) cycle
                j1=j1+1
                do m=1,2*L+1
                   if ( abs(rrr(i1,i2,i3,m))<1.d-10 ) cycle
                   j2=j2+1
                   JJ_tmp_0(0,j2)=LLL2(ii1,ii2,ii3)
                   JJ_tmp_0(1,j2)=k1
                   JJ_tmp_0(2,j2)=k2
                   JJ_tmp_0(3,j2)=k3
                   JJ_tmp_0(4,j2)=i1
                   JJ_tmp_0(5,j2)=i2
                   JJ_tmp_0(6,j2)=i3
                   JJ_tmp_0(7,j2)=ii1
                   JJ_tmp_0(8,j2)=ii2
                   JJ_tmp_0(9,j2)=ii3
                   exit
                end do
             end do
             end do
             end do

             MMJJ_t = max(MMJJ_t,j1) 
             MMJJ   = max(MMJJ,j2)
             nzlma  = lma+2*L+1

             do m=1,2*L+1
                lma=lma+1
                icheck_tmp3(a,iorb,m)=lma
             end do

             lma=lma0
             do m=1,2*L+1
                lma=lma+1
                do i=1,j2
                   k1 =JJ_tmp_0(1,i)
                   k2 =JJ_tmp_0(2,i)
                   k3 =JJ_tmp_0(3,i)
                   i1 =JJ_tmp_0(4,i)
                   i2 =JJ_tmp_0(5,i)
                   i3 =JJ_tmp_0(6,i)
                   ii1=JJ_tmp_0(7,i)
                   ii2=JJ_tmp_0(8,i)
                   ii3=JJ_tmp_0(9,i)
                   rtmp3(i,lma)=rrr(i1,i2,i3,m)
                   jtmp3(1,i,lma)=JJ_tmp_0(7,i)
                   jtmp3(2,i,lma)=JJ_tmp_0(8,i)
                   jtmp3(3,i,lma)=JJ_tmp_0(9,i)
                   jtmp3(4,i,lma)=JJ_tmp_0(1,i)
                   jtmp3(5,i,lma)=JJ_tmp_0(2,i)
                   jtmp3(6,i,lma)=JJ_tmp_0(3,i)
                end do
                mtmp3(lma)=j2
             end do

             lma0=lma

          end if ! icheck_tmp1(myrank_g)>0

       end do ! iorb
    end do ! a

    deallocate( JJ_tmp_0 )
    deallocate( irad )
    deallocate( uVtmp )
    deallocate( icheck_tmp1 )
    deallocate( rrr )
    deallocate( sss )

    call mpi_allreduce(nzlma,nzlma_0,1,mpi_integer,mpi_max,comm_grid,ierr)

    nzlma_0 = min(nzlma_0*2,Mlma)
    nrlma   = 0

    allocate( lma_nsend_tmp(0:nprocs_g-1) )
    allocate( icheck_tmp1(0:nprocs_g-1)   )
    allocate( icheck_tmp2(0:nprocs_g-1)   )
    allocate( nl_rank_map_tmp(0:nprocs_g) )
    allocate( maps_tmp(nzlma_0,6) )
    allocate( sendmap0(nzlma_0,0:nprocs_g-1) )
    allocate( recvmap0(nzlma_0,0:nprocs_g-1) )
    n=max(np1,np2,np3)
    allocate( itmp(n,3) ) ; itmp=0

    maps_tmp(:,:)      = 0
    sendmap0(:,:)      = 0
    recvmap0(:,:)      = 0
    icheck_tmp1(:)     = 0
    icheck_tmp2(:)     = 0
    lma_nsend_tmp(:)   = 0
    nl_rank_map_tmp(:) =-1

    do a=1,MI
       ik=Kion(a)
    do iorb=1,norb(ik)
       L=lo(iorb,ik)
    do m=-L,L
       icheck_tmp1(:)=0
       call mpi_allgather(icheck_tmp3(a,iorb,m+L+1),1,mpi_integer &
                         ,icheck_tmp1,1,mpi_integer,comm_grid,ierr)
       itmp(:,:)=0
       do n=0,nprocs_g-1
          if ( icheck_tmp1(n)==0 ) cycle
          itmp(LLp(1,n),1)=LLp(1,n)
          itmp(LLp(2,n),2)=LLp(2,n)
          itmp(LLp(3,n),3)=LLp(3,n)
       end do
       k1=count( itmp(:,1)>0 )
       k2=count( itmp(:,2)>0 )
       k3=count( itmp(:,3)>0 )
       ic1=0
       id1=np1
       do i=1,np1
          if ( ic1==0 .and. itmp(i,1)/=0 ) then
             ic1=i
          else if ( ic1/=0 .and. itmp(i,1)==0 ) then
             id1=i-1
             exit
          end if
       end do
       if ( id1-ic1+1/=k1 ) then
          i1=0
          j1=np1
          do i=id1+1,np1
             if ( i1==0 .and. itmp(i,1)/=0 ) then
                i1=i
             else if ( i1/=0 .and. itmp(i,1)==0 ) then
                j1=i-1
                exit
             end if
          end do
          i1=i1-np1
          j1=j1-np1
          ic1=i1
       end if
       ic2=0
       id2=np2
       do i=1,np2
          if ( ic2==0 .and. itmp(i,2)/=0 ) then
             ic2=i
          else if ( ic2/=0 .and. itmp(i,2)==0 ) then
             id2=i-1
             exit
          end if
       end do
       if ( id2-ic2+1/=k2 ) then
          i2=0
          j2=np2
          do i=id2+1,np2
             if ( i2==0 .and. itmp(i,2)/=0 ) then
                i2=i
             else if ( i2/=0 .and. itmp(i,2)==0 ) then
                j2=i-1
                exit
             end if
          end do
          i2=i2-np2
          j2=j2-np2
          ic2=i2
       end if
       ic3=0
       id3=np3
       do i=1,np3
          if ( ic3==0 .and. itmp(i,3)/=0 ) then
             ic3=i
          else if ( ic3/=0 .and. itmp(i,3)==0 ) then
             id3=i-1
             exit
          end if
       end do
       if ( id3-ic3+1/=k3 ) then
          i3=0
          j3=np3
          do i=id3+1,np3
             if ( i3==0 .and. itmp(i,3)/=0 ) then
                i3=i
             else if ( i3/=0 .and. itmp(i,3)==0 ) then
                j3=i-1
                exit
             end if
          end do
          i3=i3-np3
          j3=j3-np3
          ic3=i3
       end if
       do j3=ic3,id3
       do j2=ic2,id2
       do j1=ic1,id1
          k1=mod(j1+np1-1,np1)+1
          k2=mod(j2+np2-1,np2)+1
          k3=mod(j3+np3-1,np3)+1
          k=LLLp(k1,k2,k3)
          if ( icheck_tmp1(k)==0 ) icheck_tmp1(k)=-1
       end do
       end do
       end do
       do n=0,nprocs_g-1
          if ( icheck_tmp1(n)/=0 ) then
             icheck_tmp2(n)=icheck_tmp2(n)+1
          end if
       end do
       if ( icheck_tmp1(myrank_g)/=0 ) then
          if ( icheck_tmp1(myrank_g)>0 ) then
             maps_tmp(icheck_tmp2(myrank_g),1)=icheck_tmp1(myrank_g)
          end if
          maps_tmp(icheck_tmp2(myrank_g),2)=inorm(iorb,ik)
          maps_tmp(icheck_tmp2(myrank_g),3)=a
          maps_tmp(icheck_tmp2(myrank_g),4)=L
          maps_tmp(icheck_tmp2(myrank_g),5)=m
          maps_tmp(icheck_tmp2(myrank_g),6)=iorb

          do n=0,nprocs_g-1
             if ( n==myrank_g .or. icheck_tmp1(n)==0 ) cycle
             lma_nsend_tmp(n)=lma_nsend_tmp(n)+1  
             sendmap0(lma_nsend_tmp(n),n)=icheck_tmp2(myrank_g)
             recvmap0(lma_nsend_tmp(n),n)=icheck_tmp2(n)
             if ( any(nl_rank_map_tmp(0:nrlma)==n) ) cycle
             nrlma=nrlma+1
             nl_rank_map_tmp(nrlma)=n
          end do
       end if

    end do ! m
    end do ! iorb
    end do ! a


    nzlma = icheck_tmp2(myrank_g)


!    call gv_alloc("uV")
    allocate( uV(MMJJ,nzlma) ) ; uV=0.d0
    allocate( JJ_MAP(6,MMJJ,nzlma) ) ; JJ_MAP=0
    allocate( MJJ_MAP(nzlma) ) ; MJJ_MAP=0

     
    do lma=1,nzlma
       l=maps_tmp(lma,1)
       if (l==0) cycle
       MJJ_MAP(lma)=mtmp3(l)
       do j=1,MJJ_MAP(lma)
          uV(j,lma) = rtmp3(j,l)
          JJ_MAP(1:6,j,lma) = jtmp3(1:6,j,l)
       end do
    end do
     
    deallocate( mtmp3 )
    deallocate( jtmp3 )
    deallocate( rtmp3 )


    allocate( icheck_tmp4(a1b:b1b,a2b:b2b,a3b:b3b) )
    icheck_tmp4=0
    do lma=1,nzlma
       j=0
       icheck_tmp4=0
       do i=1,MJJ_MAP(lma)
          i1=JJ_MAP(1,i,lma)
          i2=JJ_MAP(2,i,lma)
          i3=JJ_MAP(3,i,lma)
          if ( icheck_tmp4(i1,i2,i3)==0 ) then
             j=j+1
             icheck_tmp4(i1,i2,i3)=j
          end if
       end do
       MJJ(lma)=j
    end do
    MAXMJJ = maxval( MJJ(1:nzlma) )
    deallocate( icheck_tmp4 )


    do lma=1,nzlma
       l            = maps_tmp(lma,1)
       if ( l==0 ) cycle
       iuV(lma)     = maps_tmp(lma,2)
       amap(lma)    = maps_tmp(lma,3)
       lmap(lma)    = maps_tmp(lma,4)
       mmap(lma)    = maps_tmp(lma,5)
       iorbmap(lma) = maps_tmp(lma,6)
    end do

    do i=1,nrlma
       nl_rank_map(i)=nl_rank_map_tmp(i)
    end do

    deallocate( itmp )
    deallocate( nl_rank_map_tmp )
    deallocate( icheck_tmp2 )
    deallocate( icheck_tmp1 )
    deallocate( maps_tmp )
    deallocate( icheck_tmp3 )


!    call gv_alloc("lma_nsend")
    allocate( lma_nsend(0:nprocs_g-1) ) ; lma_nsend=0


    do n=0,nprocs_g-1
       sendmap(1:nl_max_send,n) = sendmap0(1:nl_max_send,n)
       lma_nsend(n) = lma_nsend_tmp(n)
    end do


    allocate( ireq(2*nprocs_g) )
    allocate( istatus(MPI_STATUS_SIZE,2*nprocs_g) )
    nreq=0
    do n=0,nprocs_g-1
       if ( lma_nsend(n)<=0 .or. n==myrank_g ) cycle
       nreq=nreq+1
       call mpi_isend(recvmap0(1,n),lma_nsend(n),mpi_integer,n,1,comm_grid,ireq(nreq),ierr)
       nreq=nreq+1
       call mpi_irecv(recvmap(1,n) ,lma_nsend(n),mpi_integer,n,1,comm_grid,ireq(nreq),ierr)
    end do
    call mpi_waitall(nreq,ireq,istatus,ierr)
    deallocate( istatus )
    deallocate( ireq )

    deallocate( recvmap0,sendmap0,lma_nsend_tmp )



    allocate( itmp(3,nrlma) ) ; itmp=0
    allocate( itmp1(nrlma), work(nrlma) )
    allocate( itmp2(nrlma),itmp3(3,nrlma) )

    do irlma=1,nrlma
       n=nl_rank_map(irlma)
       itmp(1,irlma)=LLp(1,n)-LLp(1,myrank_g)
       itmp(2,irlma)=LLp(2,n)-LLp(2,myrank_g)
       itmp(3,irlma)=LLp(3,n)-LLp(3,myrank_g)
    end do

    nrlma_xyz(1:6)=0

    m=0
    n=0
    do i=1,nrlma
       if( itmp(2,i)==0 .and. itmp(3,i)==0 .and. itmp(1,i)>0 )then
          n=n+1
          work(n)=itmp(1,i)
          itmp2(n)=i
       end if
    end do
    if ( n>0 ) then
       call indexx(n,work,itmp1)
       do i=1,n
          j=itmp2( itmp1(i) )
          itmp3(:,m+i)=itmp(:,j)
       end do
    end if
    m=m+n
    nrlma_xyz(1)=nrlma_xyz(1)+n
    n=0
    do i=1,nrlma
       if( itmp(2,i)==0 .and. itmp(3,i)==0 .and. itmp(1,i)<0 )then
          n=n+1
          work(n)=itmp(1,i)
          itmp2(n)=i
       end if
    end do
    if ( n>0 ) then
       call indexx(n,work,itmp1)
       do i=1,n
          j=itmp2(itmp1(i))
          itmp3(:,m+n-i+1)=itmp(:,j)
       end do
    end if
    m=m+n
    nrlma_xyz(2)=nrlma_xyz(2)+n

    n=0
    do i=1,nrlma
       if( itmp(1,i)==0 .and. itmp(3,i)==0 .and. itmp(2,i)>0 )then
          n=n+1
          work(n)=itmp(2,i)
          itmp2(n)=i
       end if
    end do
    if ( n>0 ) then
       call indexx(n,work,itmp1)
       do i=1,n
          j=itmp2( itmp1(i) )
          itmp3(:,m+i)=itmp(:,j)
       end do
    end if
    m=m+n
    nrlma_xyz(3)=nrlma_xyz(3)+n
    n=0
    do i=1,nrlma
       if( itmp(1,i)==0 .and. itmp(3,i)==0 .and. itmp(2,i)<0 )then
          n=n+1
          work(n)=itmp(2,i)
          itmp2(n)=i
       end if
    end do
    if ( n>0 ) then
       call indexx(n,work,itmp1)
       do i=1,n
          j=itmp2(itmp1(i))
          itmp3(:,m+n-i+1)=itmp(:,j)
       end do
    end if
    m=m+n
    nrlma_xyz(4)=nrlma_xyz(4)+n

    n=0
    do i=1,nrlma
       if( itmp(1,i)==0 .and. itmp(2,i)==0 .and. itmp(3,i)>0 )then
          n=n+1
          work(n)=itmp(3,i)
          itmp2(n)=i
       end if
    end do
    if ( n>0 ) then
       call indexx(n,work,itmp1)
       do i=1,n
          j=itmp2( itmp1(i) )
          itmp3(:,m+i)=itmp(:,j)
       end do
    end if
    m=m+n
    nrlma_xyz(5)=nrlma_xyz(5)+n
    n=0
    do i=1,nrlma
       if( itmp(1,i)==0 .and. itmp(2,i)==0 .and. itmp(3,i)<0 )then
          n=n+1
          work(n)=itmp(3,i)
          itmp2(n)=i
       end if
    end do
    if ( n>0 ) then
       call indexx(n,work,itmp1)
       do i=1,n
          j=itmp2(itmp1(i))
          itmp3(:,m+n-i+1)=itmp(:,j)
       end do
    end if
    m=m+n
    nrlma_xyz(6)=nrlma_xyz(6)+n


    n=maxval( nrlma_xyz )
    allocate( num_2_rank(n,6) )
    num_2_rank(:,:)=MPI_PROC_NULL


    m=0
    do i=1,6
       do j=1,nrlma_xyz(i)
          m=m+1
          i1=itmp3(1,m)+LLp(1,myrank_g)
          i2=itmp3(2,m)+LLp(2,myrank_g)
          i3=itmp3(3,m)+LLp(3,myrank_g)
          k=LLLp(i1,i2,i3)
          num_2_rank(j,i)=k
       end do
    end do

    deallocate( itmp,itmp1,itmp2,itmp3,work )

    do i=1,5,2
       n=max( nrlma_xyz(i),nrlma_xyz(i+1) )
       nrlma_xyz(i)=n
       nrlma_xyz(i+1)=n
    end do


!    call gv_alloc("sbufnl")
    n=maxval( lma_nsend )
    allocate( sbufnl(n,0:nprocs_g-1) ) ; sbufnl=zero
    allocate( rbufnl(n,0:nprocs_g-1) ) ; rbufnl=zero


!    call gv_alloc("uVk_p1p2")
    allocate( JJP(MAXMJJ,nzlma) ) ; JJP=0
    allocate( uVk(MAXMJJ,nzlma,MBZ_0:MBZ_1) )


    allocate( icheck_tmp4(a1b:b1b,a2b:b2b,a3b:b3b) )
    icheck_tmp4=0

    c1=1.d0/ML1
    c2=1.d0/ML2
    c3=1.d0/ML3

    do k=MBZ_0,MBZ_1

       d1=pi2*kbb(1,k)
       d2=pi2*kbb(2,k)
       d3=pi2*kbb(3,k)

       do lma=1,nzlma
          j=0
          icheck_tmp4=0
          do i=1,MJJ_MAP(lma)
             i1=JJ_MAP(1,i,lma)
             i2=JJ_MAP(2,i,lma)
             i3=JJ_MAP(3,i,lma)
             k1=JJ_MAP(4,i,lma)
             k2=JJ_MAP(5,i,lma)
             k3=JJ_MAP(6,i,lma)
             j3=icheck_tmp4(i1,i2,i3)
             kr=d1*(c1*i1+k1)+d2*(c2*i2+k2)+d3*(c3*i3+k3)
             ztmp0=dcmplx(cos(kr),-sin(kr))*uV(i,lma)
             if ( j3==0 ) then
                j=j+1
                icheck_tmp4(i1,i2,i3)=j
                uVk(j,lma,k)=ztmp0
                JJP(j,lma)=LLL2(i1,i2,i3)
             else
                uVk(j3,lma,k)=uVk(j3,lma,k)+ztmp0
             end if
          end do
       end do
    end do
    deallocate( icheck_tmp4 )

  END SUBROUTINE prep_pselect_1
#endif

END MODULE ps_nloc1_module
