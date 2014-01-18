MODULE ps_nloc3_module

  use pseudopot_module
  use bb_module
  use aa_module
  use atom_module
  use array_bound_module
  use ggrid_module
  use bz_module
  use rgrid_module
  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: init_ps_nloc3, prep_ps_nloc3, op_ps_nloc3

  integer :: m_mapgk, m_norb
  integer,allocatable :: mapgk(:,:)
  real(8),allocatable :: viodgk2(:,:,:)

  integer,allocatable :: iuV(:),iorbmap(:),amap(:),lmap(:),mmap(:)

  integer :: Mlma
  integer :: nzlma

#ifdef _DRSDFT_
  integer,parameter :: TYPE_MAIN=MPI_REAL8
  real(8),allocatable :: uVk(:,:,:)
  real(8),allocatable :: utmp(:)
  real(8),allocatable :: utmp3(:,:,:)
  real(8),parameter :: zero=0.0d0,one=1.0d0
#else
  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
  complex(8),allocatable :: uVk(:,:,:)
  complex(8),allocatable :: utmp(:)
  complex(8),allocatable :: utmp3(:,:,:)
  complex(8),parameter :: zero=(0.0d0,0.0d0),one=(1.0d0,0.0d0)
#endif

CONTAINS


  SUBROUTINE init_ps_nloc3
    implicit none
    real(8)   ,parameter :: ep=1.d-14, dr=2.d-3
    integer :: ik,ig,i,i0,j,iorb,L,lm,m,iq,k,NRc,MMr,ierr,MKI
    integer :: m0,m1,m2,m3,n1,n2,n3,n,morb,mm,kkk,a,lma,lm0
    integer :: Nqc0,iloc(1),iwork(3),MG
    real(8) :: q,x,qc,Rc,p1,p2,p3,p4,vlong,r,G,G2,sum0,Vcell
    real(8) :: b1,b2,b3,pi4,qc2,const,c1,c2,r1,r2,y,dy,y0,dy0,sb
    real(8) :: c,q1,q2,q3,s,work(3),maxerr,lambda0,eta0
    real(8) :: sb0x,sb0y,sb1x,sb1y,sb2,sb3
    integer,allocatable :: indx(:),itmp(:),itmp2(:)
    real(8),allocatable :: tmp(:),tmp2(:)

    pi4   = 4.0d0*acos(-1.0d0)
    Vcell = abs(va)
    MKI   = Nelement
    MG    = NGgrid(0)

    allocate( mapgk(MG,MBZ_0:MBZ_1) ) ; mapgk=0
    
    call construct_Ggrid(0)

    n=(MBZ_1-MBZ_0+1)*MG
    allocate( tmp(n),indx(n),tmp2(n),itmp(n),itmp2(n) )

    itmp=0
    itmp2=0
    i=0
    do k=MBZ_0,MBZ_1
       do ig=1,MG
          i=i+1
          b1=kbb(1,k)+LLG(1,ig)
          b2=kbb(2,k)+LLG(2,ig)
          b3=kbb(3,k)+LLG(3,ig)
          tmp(i)=( bb(1,1)*b1+bb(1,2)*b2+bb(1,3)*b3 )**2 &
                +( bb(2,1)*b1+bb(2,2)*b2+bb(2,3)*b3 )**2 &
                +( bb(3,1)*b1+bb(3,2)*b2+bb(3,3)*b3 )**2
       end do
    end do
    call indexx(n,tmp,indx)

    m=0
    i=0
    q=1.d50
    do k=MBZ_0,MBZ_1
       do ig=1,MG
          i=i+1
          r=tmp(indx(i))
          if ( abs(q-r)>ep ) then
             m=m+1
             q=r
             tmp2(m)=r
             itmp(m)=1
             itmp2(indx(i))=m
          else
             itmp(m)=itmp(m)+1
             itmp2(indx(i))=m
          end if
       end do
    end do

    i=0
    do k=MBZ_0,MBZ_1
       do ig=1,MG
          i=i+1
          mapgk(ig,k)=itmp2(i)
       end do
    end do

    m_mapgk = m
    m_norb  = maxval(norb)

    deallocate( tmp,itmp,itmp2,indx )

    call destruct_Ggrid

!- allocate --------------------------------
    allocate( viodgk2(m_mapgk,m_norb,MKI) )
    viodgk2=0.0d0
 !------------------------------------------

    NRc=maxval(NRps)
    allocate( tmp(NRc) ) ; tmp=0.d0

    const=pi4/Vcell
    do ik=1,MKI
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          NRc=NRps(iorb,ik)
          select case(L)
          case(0)
             do i=1,m
                q=tmp2(i)
                if ( q==0.d0 ) then
                   do j=1,NRc
                      tmp(j)=rad(j,ik)*viod(j,iorb,ik)*rab(j,ik)
                   end do
                   call simp(tmp(1:NRc),sum0,2)
                else
                   q=sqrt(q)
                   do j=1,NRc
                      r=rad(j,ik)
                      x=q*r
                      if ( x<1.d-1 ) then
                         sb=-(1.d0/39916800.d0*x**10-1.d0/362880.d0*x**8+1.d0/5040.d0*x**6 &
                              -1.d0/120.d0*x**4+1.d0/6.d0*x**2-1.d0)
                      else
                         sb=sin(x)/x
                      end if
                      tmp(j)=r*viod(j,iorb,ik)*sb*rab(j,ik)
                   end do
                   call simp(tmp(1:NRc),sum0,2)
                end if
                viodgk2(i,iorb,ik)=const*sum0
             end do
          case(1)
             do i=1,m
                q=tmp2(i)
                if ( q==0.d0 ) then
                   sum0=0.d0
                else
                   q=sqrt(q)
                   do j=1,NRc
                      x=q*rad(j,ik)
                      if ( x<1.d-1 ) then
                         sb=1.d0/3991680.d0*x**9-1.d0/45360.d0*x**7 &
                              +1.d0/840.d0*x**5-1.d0/30.d0*x**3+1.d0/3.d0*x
                      else
                         sb=sin(x)/x**2-cos(x)/x
                      end if
                      tmp(j)=rad(j,ik)*viod(j,iorb,ik)*sb*rab(j,ik)
                   end do
                   call simp(tmp(1:NRc),sum0,2)
                end if
                viodgk2(i,iorb,ik)=const*sum0
             end do
          case(2)
             do i=1,m
                q=tmp2(i)
                if ( q==0.d0 ) then
                   sum0=0.d0
                else
                   q=sqrt(q)
                   do j=1,NRc
                      x=q*rad(j,ik)
                      if ( x<1.d-1 ) then
                         sb=1.d0/51891840.d0*x**10-1.d0/498960.d0*x**8 &
                              +1.d0/7560.d0*x**6-1.d0/210.d0*x**4+1.d0/15.d0*x**2
                      else
                         sb=(3.d0-1.d0*x**2)/x**3*sin(x)-3.d0/x**2*cos(x)
                      end if
                      tmp(j)=rad(j,ik)*viod(j,iorb,ik)*sb*rab(j,ik)
                   end do
                   call simp(tmp(1:NRc),sum0,2)
                end if
                viodgk2(i,iorb,ik)=const*sum0
             end do
          case default
             write(*,*) "PP for L>2 is not implemented."
             stop
          end select
       end do !iorb
    end do !ik

    deallocate( tmp,tmp2 )

    return
  END SUBROUTINE init_ps_nloc3

  SUBROUTINE simp(f,s,m)
    implicit none
    integer,intent(IN)  :: m
    real(8),intent(IN)  :: f(:)
    real(8),intent(OUT) :: s
    real(8),allocatable :: g(:)
    integer :: i,n,nn,nmax
    n=size(f) ; nmax=int(n/m)*m
    do i=0,m
       nmax=nmax+i ; if ( nmax>=n ) exit
    end do
    allocate( g(nmax) ) ; g(1:n)=f ; if ( nmax>n ) g(n+1:)=0.d0
    select case(m)
    case default
       s = 0.5d0*(f(1)+f(n)) + sum(f(2:n-1))
    case(2)
       s=0.d0
       do i=1,nmax-2,2
          s = s + g(i) + 4.d0*g(i+1) + g(i+2)
       end do
       s=s/3.d0
    case(4)
       s=0.d0
       do i=1,nmax-4,4
          s=s+7*g(i)+32*g(i+1)+12*g(i+2)+32*g(i+3)+7*g(i+4)
       end do
       s=s*2.d0/45.d0
    case(6)
       s=0.d0
       do i=1,nmax-6,6
          s=s+41*g(i)+216*g(i+1)+27*g(i+2)+272*g(i+3)+27*g(i+4)+216*g(i+5)+41*g(i+6)
       end do
       s=s/140.d0
    end select
    deallocate( g )
    return
  END SUBROUTINE simp


  SUBROUTINE prep_ps_nloc3
    implicit none
    complex(8),parameter :: zi=(0.d0,1.d0)
    complex(8),allocatable :: zwork(:,:,:),fftwork(:)
    complex(8) :: zum,ztmp,phase
    integer,allocatable :: icnt(:),idis(:)
    integer :: nn1,nn2,ML0,Mlma_np,l1,l2,l3,m1,m2,m3
    integer :: i1,i2,i3,j1,j2,j3,irank_g,ierr
    integer :: ik,i,j,k,a,L,m,iorb,lma,lma0,lma1
    real(8) :: a1,a2,a3,c1,c2,c3,Gr,x,y,z,pi2
    integer :: k1,kk1,ig,ML1,ML2,ML3,ML,MG
    integer :: ifacx(30),ifacy(30),ifacz(30)
    integer,allocatable :: lx1(:),lx2(:),ly1(:),ly2(:),lz1(:),lz2(:)
    complex(8),allocatable :: wsavex(:),wsavey(:),wsavez(:)

    INTERFACE
       FUNCTION Ylm(x,y,z,l,m)
         real(8) :: Ylm
         real(8),intent(IN) :: x,y,z
         integer,intent(IN) :: l,m
       END FUNCTION Ylm
    END INTERFACE

    nn1     = idisp(myrank)+1
    nn2     = idisp(myrank)+ircnt(myrank)
    ML0     = ircnt(myrank)
    pi2     = 2.d0*acos(-1.0d0)
    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    MG  = NGgrid(0)

    Mlma=0
    do i=1,Natom
       ik=ki_atom(i)
       do iorb=1,norb(ik)
          Mlma=Mlma+2*lo(iorb,ik)+1
       end do
    end do

    if ( Mlma <= 0 ) return

! allocate ---------------------------------------------
    allocate( uVk(nn1:nn2,Mlma,MBZ_0:MBZ_1) ) ; uVk(:,:,:)=zero
    allocate( iuV(Mlma)     ) ; iuV=0
    allocate( amap(Mlma)    ) ; amap=0
    allocate( lmap(Mlma)    ) ; lmap=0
    allocate( mmap(Mlma)    ) ; mmap=0
    allocate( iorbmap(Mlma) ) ; iorbmap=0
!-------------------------------------------------------

!- allocate -----------------------------------------------------
    allocate( zwork(0:ML1-1,0:ML2-1,0:ML3-1),fftwork(ML) )
    allocate( idis(0:nprocs_g-1),icnt(0:nprocs_g-1) )
    allocate( utmp(ML) )
    allocate( utmp3(nn1:nn2,0:nprocs_g-1,0:nprocs_b-1) )
!----------------------------------------------------------------

    Mlma_np = (Mlma+(nprocs_g*nprocs_b)-1)/(nprocs_g*nprocs_b)

    icnt(0:nprocs_g-1)=ML0
    idis(0)=0
    do i=1,nprocs_g-1
       idis(i)=sum( icnt(0:i) )-icnt(i)
    end do

    lma=0
    do a=1,Natom
       ik=ki_atom(a)
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          do m=-L,L
             lma=lma+1
             lmap(lma)=L
             mmap(lma)=m
             amap(lma)=a
             iorbmap(lma)=iorb
             iuV(lma)=inorm(iorb,ik)
          end do
       end do
    end do

    allocate( lx1(ML),lx2(ML),ly1(ML),ly2(ML),lz1(ML),lz2(ML) )
    allocate( wsavex(ML1),wsavey(ML2),wsavez(ML3) )

    call prefft(ML1,ML2,ML3,ML,wsavex,wsavey,wsavez &
         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call construct_Ggrid(0)

    do k=MBZ_0,MBZ_1
       lma=1+myrank_b*nprocs_g+myrank_g-nprocs_g*nprocs_b
       do lma0=1,Mlma_np
          lma=lma+nprocs_g*nprocs_b
          utmp(:)=zero
          if ( lma <= Mlma ) then
             a=amap(lma)
             L=lmap(lma)
             m=mmap(lma)
             iorb=iorbmap(lma)
             a1=aa_atom(1,a)*pi2
             a2=aa_atom(2,a)*pi2
             a3=aa_atom(3,a)*pi2
             ik=ki_atom(a)
             phase=(-zi)**L
             zwork=(0.d0,0.d0)
             do i=1,MG
                i1=mod(ML1+LLG(1,i),ML1)
                i2=mod(ML2+LLG(2,i),ML2)
                i3=mod(ML3+LLG(3,i),ML3)
                Gr=LLG(1,i)*a1+LLG(2,i)*a2+LLG(3,i)*a3
                c1=kbb(1,k)+LLG(1,i)
                c2=kbb(2,k)+LLG(2,i)
                c3=kbb(3,k)+LLG(3,i)
                x=bb(1,1)*c1+bb(1,2)*c2+bb(1,3)*c3
                y=bb(2,1)*c1+bb(2,2)*c2+bb(2,3)*c3
                z=bb(3,1)*c1+bb(3,2)*c2+bb(3,3)*c3
                if ( x==0.d0 .and. y==0.d0 .and. z==0.d0 .and. L/=0 ) cycle
                j=mapgk(i,k)
                zwork(i1,i2,i3)=viodgk2(j,iorb,ik) &
                     *phase*Ylm(x,y,z,L,m)*dcmplx(cos(Gr),-sin(Gr))
             end do
             call fft3bx(ML1,ML2,ML3,ML,zwork,fftwork,wsavex,wsavey,wsavez &
                  ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)
             i=0
             irank_g=-1
             do i3=1,node_partition(3)
             do i2=1,node_partition(2)
             do i1=1,node_partition(1)
                irank_g=irank_g+1
                l1=pinfo_grid(1,irank_g)
                m1=pinfo_grid(2,irank_g)+l1-1
                l2=pinfo_grid(3,irank_g)
                m2=pinfo_grid(4,irank_g)+l2-1
                l3=pinfo_grid(5,irank_g)
                m3=pinfo_grid(6,irank_g)+l3-1
                do j3=l3,m3
                do j2=l2,m2
                do j1=l1,m1
                   i=i+1
                   utmp(i)=zwork(j1,j2,j3)
                end do
                end do
                end do
             end do
             end do
             end do

          end if ! lma<=Mlma

          call mpi_alltoallv(utmp,ir_grid,id_grid,TYPE_MAIN,utmp3(nn1,0,myrank_b) &
               ,icnt,idis,TYPE_MAIN,comm_grid,ierr)

          call mpi_allgather(utmp3(nn1,0,myrank_b),ML0*nprocs_g,TYPE_MAIN &
               ,utmp3(nn1,0,0),ML0*nprocs_g,TYPE_MAIN,comm_band,ierr)

          lma1=lma-(myrank_b*nprocs_g+myrank_g)-1
          do j=0,nprocs_b-1
          do i=0,nprocs_g-1
             lma1=lma1+1
             if ( lma1 <= Mlma ) then
                uVk(nn1:nn2,lma1,k)=utmp3(nn1:nn2,i,j)
             end if
          end do
          end do

       end do ! lma0

    end do ! k

    call destruct_Ggrid

    deallocate( wsavez,wsavey,wsavex )
    deallocate( lz2,lz1,ly2,ly1,lx2,lx1 )

    deallocate( utmp3 )
    deallocate( utmp  )
    deallocate( icnt,idis  )
    deallocate( zwork,fftwork )

    nzlma = Mlma

    return

  END SUBROUTINE prep_ps_nloc3



  SUBROUTINE op_ps_nloc3(k,tpsi,htpsi,n1,n2,ib1,ib2)
    implicit none
    integer,intent(IN) :: k,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    real(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    real(8),allocatable :: uVunk(:,:),uVunk0(:,:)
    character(1),parameter :: TRANSA='T'
    character(1),parameter :: TRANSB='N'
#else
    character(1),parameter :: TRANSA='C'
    character(1),parameter :: TRANSB='N'
    complex(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    complex(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    complex(8),allocatable :: uVunk(:,:),uVunk0(:,:)
#endif
    integer :: ML0,nb,lma,ib,ierr

    if ( Mlma <= 0 ) return

    ML0 = n2-n1+1
    nb  = ib2-ib1+1

    allocate( uVunk(Mlma,ib1:ib2),uVunk0(Mlma,ib1:ib2) )

#ifdef _DRSDFT_
    call dgemm(TRANSA,TRANSB,Mlma,nb,ML0,zdV,uVk(n1,1,k),ML0 &
         ,tpsi(n1,ib1),ML0,zero,uVunk0(1,ib1),Mlma)
#else
    call zgemm(TRANSA,TRANSB,Mlma,nb,ML0,zdV,uVk(n1,1,k),ML0 &
         ,tpsi(n1,ib1),ML0,zero,uVunk0(1,ib1),Mlma)
#endif

    call mpi_allreduce(uVunk0,uVunk,Mlma*nb,TYPE_MAIN,mpi_sum,comm_grid,ierr)

    do ib=ib1,ib2
!$OMP PARALLEL DO
       do lma=1,Mlma
          uVunk(lma,ib)=iuV(lma)*uVunk(lma,ib)
       end do
!$OMP END PARALLEL DO
    end do

#ifdef _DRSDFT_
    call dgemm(TRANSB,TRANSB,ML0,nb,Mlma,one,uVk(n1,1,k),ML0 &
         ,uVunk(1,ib1),Mlma,one,htpsi(n1,ib1),ML0)
#else
    call zgemm(TRANSB,TRANSB,ML0,nb,Mlma,one,uVk(n1,1,k),ML0 &
         ,uVunk(1,ib1),Mlma,one,htpsi(n1,ib1),ML0)
#endif

    deallocate( uVunk0,uVunk )

  END SUBROUTINE op_ps_nloc3


END MODULE ps_nloc3_module
