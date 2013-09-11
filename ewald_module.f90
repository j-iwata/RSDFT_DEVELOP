MODULE ewald_module

  use aa_module
  use bb_module
  use atom_module
  use pseudopot_module, only: Zps
  use parallel_module

  implicit none
  PRIVATE
  PUBLIC :: Eewald,calc_ewald,calc_force_ewald,test_ewald

  real(8) :: Eewald
  real(8) :: eta=0.d0
  real(8) :: qgb=0.d0
  real(8) :: qrb=0.d0

  INTERFACE
     FUNCTION bberf(x)
       real(8) :: bberf,x
     END FUNCTION bberf
  END INTERFACE

CONTAINS

  SUBROUTINE test_ewald(Ewld,disp_switch)
    real(8),intent(IN) :: Ewld
    logical,intent(IN)  :: disp_switch
    if ( disp_switch ) write(*,*) "test_ewald(dummy routine)"
    return
  END SUBROUTINE test_ewald

  SUBROUTINE calc_ewald(Ewld,disp_switch)
    real(8),intent(OUT) :: Ewld
    logical,intent(IN) :: disp_switch
    integer,parameter :: maxloop=20
    integer :: i,j,k,n,i1,i2,loop,ik1,ik2,a,b,ik,MREtmp,MGEtmp
    integer :: MGE,MRE,mmg1,mmg2,mmg3,mmr1,mmr2,mmr3,mmg,mmr
    integer :: ierr,i2_0,i2_1,MI_NP
    integer,allocatable :: ER(:,:),EG(:,:),indx(:),iwork(:,:)
    real(8),parameter :: ep=1.d-15
    real(8) :: Ewldg,Ewldr,Ewldg0,Ewldr0,Ewld0,bbmax,aamax
    real(8) :: x,y,z,ss,qqg,qqa,sum1,t,Qtot,Qtot2,sqeta,qgb0,qrb0
    real(8) :: c1,c2,c3,a1,a2,a3,a1_1,a2_1,a3_1,a1_2,a2_2,a3_2
    real(8) :: const,pi2,pi,Vcell,z1,z2
    real(8),allocatable :: SSR(:),SSG(:),SG2(:),work(:)
    complex(8) :: sg

    pi=acos(-1.d0)
    pi2=2.d0*pi
    Vcell=abs(Va)

    if ( eta == 0.d0 ) eta=pi/Vcell**(2.d0/3.d0)
    if ( qgb == 0.d0 ) qgb=12.d0*sqrt(eta)
    if ( qrb == 0.d0 ) qrb=6.d0/sqrt(eta)

    qgb0=qgb
    qrb0=qrb
    sqeta=sqrt(eta)
    Ewld0 =1.d10
    Ewldg0=1.d10
    Ewldr0=1.d10

    aamax=0.d0
    do i=1,3
       x=sum(aa(:,i)**2)
       if( x > aamax ) aamax=x
    end do
    aamax=sqrt(aamax)+ep
    bbmax=0.d0
    do i=1,3
       x=sum(bb(:,i)**2)
       if( x > bbmax ) bbmax=x
    end do
    bbmax=sqrt(bbmax)+ep

    do loop=1,maxloop

       qqg=qgb*qgb
       qqa=qrb*qrb

       mmg1=int( sqrt(qqg/sum(bb(:,1)**2)) )+1
       mmg2=int( sqrt(qqg/sum(bb(:,2)**2)) )+1
       mmg3=int( sqrt(qqg/sum(bb(:,3)**2)) )+1
       mmg=(2*mmg1+1)*(2*mmg2+1)*(2*mmg3+1)

       mmr1=int( sqrt(qqa/sum(aa(:,1)**2)) )+1
       mmr2=int( sqrt(qqa/sum(aa(:,2)**2)) )+1
       mmr3=int( sqrt(qqa/sum(aa(:,3)**2)) )+1
       mmr=(2*mmr1+1)*(2*mmr2+1)*(2*mmr3+1)

       allocate( EG(3,mmg), SSG(mmg), ER(3,mmr), SSR(mmr) )

       SSG=0.d0
       SSR=0.d0

       n=0
       do k=-mmg3,mmg3
       do j=-mmg2,mmg2
       do i=-mmg1,mmg1
          if( i==0 .and. j==0 .and. k==0 )cycle
          x=i*bb(1,1)+j*bb(1,2)+k*bb(1,3)
          y=i*bb(2,1)+j*bb(2,2)+k*bb(2,3)
          z=i*bb(3,1)+j*bb(3,2)+k*bb(3,3)
          ss=x*x+y*y+z*z
          if ( ss<=qqg-ep ) then
             n=n+1
             EG(1,n)=i
             EG(2,n)=j
             EG(3,n)=k
             SSG(n)=ss
          end if
       end do
       end do
       end do
       MGE=n

       n=0
       do k=-mmr3,mmr3
       do j=-mmr2,mmr2
       do i=-mmr1,mmr1
          x=i*aa(1,1)+j*aa(1,2)+k*aa(1,3)
          y=i*aa(2,1)+j*aa(2,2)+k*aa(2,3)
          z=i*aa(3,1)+j*aa(3,2)+k*aa(3,3)
          ss=x*x+y*y+z*z
          if ( ss<=qqa-ep ) then
             n=n+1
             ER(1,n)=i
             ER(2,n)=j
             ER(3,n)=k
             SSR(n)=ss
          end if
       end do
       end do
       end do
       MRE=n

! Sorting

       n=max(MRE,MGE)
       allocate( indx(n),iwork(3,n),work(n) )

       call indexx(MGE,SSG(1),indx(1))
       iwork(1,1:MGE)=EG(1,1:MGE)
       iwork(2,1:MGE)=EG(2,1:MGE)
       iwork(3,1:MGE)=EG(3,1:MGE)
       work(1:MGE)=SSG(1:MGE)
       do n=1,MGE
          EG(1,n)=iwork(1,indx(n))
          EG(2,n)=iwork(2,indx(n))
          EG(3,n)=iwork(3,indx(n))
          SSG(n)=work(indx(n))
       end do
       call indexx(MRE,SSR(1),indx(1))
       iwork(1,1:MRE)=ER(1,1:MRE)
       iwork(2,1:MRE)=ER(2,1:MRE)
       iwork(3,1:MRE)=ER(3,1:MRE)
       work(1:MRE)=SSR(1:MRE)
       do n=1,MRE
          ER(1,n)=iwork(1,indx(n))
          ER(2,n)=iwork(2,indx(n))
          ER(3,n)=iwork(3,indx(n))
          SSR(n)=work(indx(n))
       end do

       deallocate( indx,iwork,work )

       Qtot=0.d0
       Qtot2=0.d0
       do i=1,Natom
          Qtot =Qtot +Zps(ki_atom(i))
          Qtot2=Qtot2+Zps(ki_atom(i))**2
       end do

! Structure factor

       allocate( SG2(MGE) ) ; SG2=0.d0

       do i=1,MGE
          sg=(0.d0,0.d0)
          x=pi2*EG(1,i)
          y=pi2*EG(2,i)
          z=pi2*EG(3,i)
          do a=1,Natom
             t=x*aa_atom(1,a)+y*aa_atom(2,a)+z*aa_atom(3,a)
             sg=sg+Zps(ki_atom(a))*dcmplx(cos(t),-sin(t))
          end do
          SG2(i)=abs(sg)**2
       end do

! G
       Ewldg=0.d0
       const=1.d0/(4.d0*eta)
       do i=1,MGE
          t=exp(-SSG(i)*const)/SSG(i)
          if ( abs(t)<ep ) then
             MGEtmp=i-1 ; exit
          end if
          Ewldg=Ewldg+SG2(i)*t
          MGEtmp=i
       end do
       Ewldg=Ewldg*4.d0*pi/Vcell
       Ewldg=Ewldg-pi/(eta*Vcell)*Qtot**2

! R
       MI_NP = Natom/nprocs
       if ( MI_NP*nprocs < Natom ) MI_NP=MI_NP+1
       i2_0  = myrank*MI_NP+1
       i2_1  = min( i2_0+MI_NP-1, Natom )
       Ewldr=0.d0
       MREtmp=0
       do i2=i2_0,i2_1
          a1_2=aa_atom(1,i2)
          a2_2=aa_atom(2,i2)
          a3_2=aa_atom(3,i2)
          z2=Zps(ki_atom(i2))
       do i1=i2,Natom
          a1_1=aa_atom(1,i1)
          a2_1=aa_atom(2,i1)
          a3_1=aa_atom(3,i1)
          z1=Zps(ki_atom(i1))
          sum1=0.d0
          do i=1,MRE
             a1=ER(1,i)+a1_1-a1_2
             a2=ER(2,i)+a2_1-a2_2
             a3=ER(3,i)+a3_1-a3_2
             x=aa(1,1)*a1+aa(1,2)*a2+aa(1,3)*a3
             y=aa(2,1)*a1+aa(2,2)*a2+aa(2,3)*a3
             z=aa(3,1)*a1+aa(3,2)*a2+aa(3,3)*a3
             ss=sqrt(x*x+y*y+z*z)
             if ( ss==0.d0 ) cycle
             t=(1.d0-bberf(sqeta*ss))/ss
             if ( abs(t)==0.d0 ) cycle
             sum1=sum1+t
             MREtmp=max(MREtmp,i)
          end do
          if ( i1/=i2 ) sum1=sum1*2.d0
          Ewldr=Ewldr+z1*z2*sum1
       end do
       end do

       sum1=Ewldr
       call mpi_allreduce(sum1,Ewldr,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
       i=MREtmp
       call mpi_allreduce(i,MREtmp,1,mpi_integer,mpi_max,mpi_comm_world,ierr)

       Ewldr=Ewldr-2.d0*sqrt(eta/pi)*Qtot2

       Ewld=0.5d0*(Ewldr+Ewldg)

       qgb=sqrt(SSG(MGEtmp))
       qrb=sqrt(SSR(MREtmp))

       deallocate( SG2,ER,EG,SSG,SSR )

! Check convergence

       if ( disp_switch ) then
          write(*,*) "MGE,MRE",MGEtmp,MGE,MREtmp,MRE
          write(*,*) "qgb,qrb     =",qgb,qrb
          write(*,*) "Ewldg,Ewldr =",Ewldg*0.5d0,Ewldr*0.5d0
          write(*,*) "Ewld        =",Ewld
       end if

       if ( MGEtmp<MGE .and. MREtmp<MRE ) exit
       if ( abs((Ewld -Ewld0)/Ewld)<ep ) then
          qgb=qgb0
          qrb=qrb0
          exit
       end if
       qgb0=qgb
       qrb0=qrb
       if ( abs(Ewldg-Ewldg0)>ep ) qgb=qgb0+bbmax
       if ( abs(Ewldr-Ewldr0)>ep ) qrb=qrb0+aamax
       Ewld0=Ewld
       Ewldg0=Ewldg
       Ewldr0=Ewldr

    end do

    if ( loop>=maxloop ) then
       write(*,*) "Ewald sum is not converged!"
       stop
    end if

    if ( disp_switch ) then
       write(*,*) "Ewald parameters"
       write(*,*) "   eta      =",eta,sqrt(eta)
       write(*,*) "   qgb, qrb =",qgb,qrb
       write(*,*) "   Ewld     =",Ewld
    end if

    return

  END SUBROUTINE calc_ewald


  SUBROUTINE calc_force_ewald(MI,force3)
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force3(3,MI)
    real(8) :: qqg,qqa,sqeta,x,y,z,r,rr,ss,const,fewldg(3),fewldr(3)
    real(8) :: const1,const2,const3
    real(8) :: GR,c,c1,c2,sum0
    real(8) :: Vcell,pi,pi2
    real(8) :: sum1,sum2,sum3
    real(8),allocatable :: SSG(:),SSR(:),dsg2(:,:),work(:,:)
    integer :: mmg1,mmg2,mmg3,mmg,mmr1,mmr2,mmr3,mmr,MGE,MRE
    integer :: i,j,k,n,a,b,ierr
    integer,allocatable :: EG(:,:),ER(:,:)
    integer,allocatable :: id_i(:),ir_i(:)

    force3(:,:) = 0.d0

    pi    = acos(-1.d0)
    pi2   = 2.d0*pi
    Vcell = abs(Va)

    allocate( id_i(0:nprocs-1),ir_i(0:nprocs-1) )
    ir_i(0:nprocs-1)=MI/nprocs
    n=MI-sum(ir_i)
    do i=1,n
       k=mod(i-1,nprocs)
       ir_i(k)=ir_i(k)+1
    end do
    do k=0,nprocs-1
       id_i(k)=sum(ir_i(0:k))-ir_i(k)
    end do

    qqg   = qgb*qgb
    qqa   = qrb*qrb
    sqeta = sqrt(eta)

    mmg1  = int( sqrt(qqg/sum(bb(:,1)**2)) )+1
    mmg2  = int( sqrt(qqg/sum(bb(:,2)**2)) )+1
    mmg3  = int( sqrt(qqg/sum(bb(:,3)**2)) )+1
    mmg   = (2*mmg1+1)*(2*mmg2+1)*(2*mmg3+1)

    mmr1  = int( sqrt(qqa/sum(aa(:,1)**2)) )+1
    mmr2  = int( sqrt(qqa/sum(aa(:,2)**2)) )+1
    mmr3  = int( sqrt(qqa/sum(aa(:,3)**2)) )+1
    mmr   = (2*mmr1+1)*(2*mmr2+1)*(2*mmr3+1)

    allocate( EG(3,mmg), SSG(mmg), ER(3,mmr), SSR(mmr) )

    SSG=0.d0
    SSR=0.d0

    n=0
    do k=-mmg3,mmg3
    do j=-mmg2,mmg2
    do i=-mmg1,mmg1
       if( i==0 .and. j==0 .and. k==0 )cycle
       x=i*bb(1,1)+j*bb(1,2)+k*bb(1,3)
       y=i*bb(2,1)+j*bb(2,2)+k*bb(2,3)
       z=i*bb(3,1)+j*bb(3,2)+k*bb(3,3)
       ss=x*x+y*y+z*z
       if( ss<=qqg )then
          n=n+1
          EG(1,n)=i
          EG(2,n)=j
          EG(3,n)=k
          SSG(n)=ss
       end if
    end do
    end do
    end do
    MGE=n

    n=0
    do k=-mmr3,mmr3
    do j=-mmr2,mmr2
    do i=-mmr1,mmr1
       x=i*aa(1,1)+j*aa(1,2)+k*aa(1,3)
       y=i*aa(2,1)+j*aa(2,2)+k*aa(2,3)
       z=i*aa(3,1)+j*aa(3,2)+k*aa(3,3)
       ss=x*x+y*y+z*z
       if( ss<=qqa )then
          n=n+1
          ER(1,n)=i
          ER(2,n)=j
          ER(3,n)=k
          SSR(n)=ss
       end if
    end do
    end do
    end do
    MRE=n

    allocate( dsg2(3,MGE) ) ; dsg2=0.d0

    const1 = 2.d0*sqeta/sqrt(Pi)
    const2 = 4.d0*pi/Vcell
    const3 = 0.25d0/eta

    do b=id_i(myrank)+1,id_i(myrank)+ir_i(myrank)

       c1=2.d0*Zps(ki_atom(b))

       fewldg(1:3)=0.d0
       do i=1,MGE
          sum0=0.d0
!$OMP parallel do reduction (+:sum0) private(GR)
          do a=1,Natom
             GR=pi2*( EG(1,i)*(aa_atom(1,a)-aa_atom(1,b)) &
                     +EG(2,i)*(aa_atom(2,a)-aa_atom(2,b)) &
                     +EG(3,i)*(aa_atom(3,a)-aa_atom(3,b)) )
             sum0=sum0+Zps(ki_atom(a))*sin(GR)
          end do
!$OMP end parallel do
          dsg2(1,i)=c1*sum(bb(1,1:3)*EG(1:3,i))*sum0
          dsg2(2,i)=c1*sum(bb(2,1:3)*EG(1:3,i))*sum0
          dsg2(3,i)=c1*sum(bb(3,1:3)*EG(1:3,i))*sum0
       end do
       sum1=0.d0
       sum2=0.d0
       sum3=0.d0
!$OMP parallel do reduction(+:sum1,sum2,sum3) private(c)
       do i=1,MGE
          c=exp(-SSG(i)*const3)/SSG(i)
          sum1=sum1+dsg2(1,i)*c
          sum2=sum2+dsg2(2,i)*c
          sum3=sum3+dsg2(3,i)*c
       end do
!$OMP end parallel do
       fewldg(1)=sum1*const2
       fewldg(2)=sum2*const2
       fewldg(3)=sum3*const2

       fewldr(1:3)=0.d0
       do a=1,MI
          if ( a==b ) cycle
          c2=c1*Zps(ki_atom(a))
          sum1=0.d0
          sum2=0.d0
          sum3=0.d0
!$OMP parallel do reduction(+:sum1,sum2,sum3) private(x,y,z,r,rr,c)
          do i=1,MRE
             x=sum( aa(1,1:3)*(aa_atom(1:3,b)-aa_atom(1:3,a)+ER(1:3,i)) )
             y=sum( aa(2,1:3)*(aa_atom(1:3,b)-aa_atom(1:3,a)+ER(1:3,i)) )
             z=sum( aa(3,1:3)*(aa_atom(1:3,b)-aa_atom(1:3,a)+ER(1:3,i)) )
             rr=x*x+y*y+z*z
             r=sqrt(rr)
             c=-c2*( (1.d0-bberf(sqeta*r))+const1*exp(-eta*rr)*r )/(r*rr)
             sum1=sum1+c*x
             sum2=sum2+c*y
             sum3=sum3+c*z
          end do
!$OMP end parallel do
          fewldr(1)=fewldr(1)+sum1
          fewldr(2)=fewldr(2)+sum2
          fewldr(3)=fewldr(3)+sum3
       end do

       force3(1,b)=-0.5d0*( fewldg(1) + fewldr(1) )
       force3(2,b)=-0.5d0*( fewldg(2) + fewldr(2) )
       force3(3,b)=-0.5d0*( fewldg(3) + fewldr(3) )

    end do ! b

    allocate( work(3,MI) )
    work(1:3,1:MI)=force3(1:3,1:MI)
    call mpi_allreduce(work,force3,3*MI,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    deallocate( work )

    deallocate( dsg2 )
    deallocate( SSR, ER, SSG, EG )

    deallocate( ir_i, id_i )

  END SUBROUTINE calc_force_ewald


END MODULE ewald_module
