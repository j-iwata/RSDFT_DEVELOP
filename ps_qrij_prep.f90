MODULE PSQRijPrep

  use parallel_module, only: myrank,myrank_g,nprocs_g,np_grid &
       ,node_partition,COMM_GRID,disp_switch_parallel,MB_d
  use aa_module, only: aa
  use atom_module, only: Natom,ki_atom,aa_atom
  use rgrid_module, only: Igrid,Ngrid,Hgrid
  use VarPSMember
  use VarPSMemberG
  use VarParaPSnonLocG
  use ps_nloc2_module, only: prepMapsTmp
  use ps_nloc2_variables, only: amap,lmap,mmap,iorbmap
  use pseudopot_module, only: pselect
  use minimal_box_module
  use para_rgrid_comm
  use watch_module
  use array_bound_module, only: ML_0
  use polint_module
  use spline_module

  implicit none

  PRIVATE
  PUBLIC :: prepQRijp102

  include 'mpif.h'
  
  complex(8),allocatable :: qaL(:,:) !qaL(k3max,Lrefmax)

  complex(8),parameter :: z0=(0.d0,0.d0),z1=(1.d0,0.d0),zi=(0.d0,1.d0)
  real(8),parameter :: pi4=16.d0*atan(1.d0)
  real(8),allocatable :: y2a(:,:,:,:)

CONTAINS

  SUBROUTINE allocateQaL
    implicit none
    if (allocated(qaL)) deallocate(qaL)
!    allocate( qaL(k2max,max_Lref) ) ; qaL=z0
    allocate( qaL(45,3) ) ; qaL=z0
    return
  END SUBROUTINE allocateQaL


  SUBROUTINE prepQRijp102
    implicit none
    real(8) :: x,y,z
    integer :: l,m

    integer :: ia,ik,ik1,ik2,ik3,il,ir,ir0,ir1
    integer :: i,j,i1,i2,i3,kk1
    integer :: iorb,iorb1,iorb2
    integer :: NRc,n
    real(8) :: Rps2

    integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab2,ab3
    integer :: ML1,ML2,ML3
    integer :: np1,np2,np3

    real(8) :: r
    integer :: mm1,mm2,mm3
    integer :: MMJJ_Q_0,nrqr

    integer :: nl3vmax
    integer,allocatable :: JJ_tmp(:,:,:,:), MJJ_tmp(:,:)

    real(8) :: c1,c2,c3
    real(8) :: maxerr

    real(8) :: Rx,Ry,Rz
    integer :: ic1,ic2,ic3,id1,id2,id3
    integer :: k1,k2,k3,i1_0,i2_0,i3_0
    real(8) :: d1,d2,d3
    real(8) :: r2

    integer :: ll3,mm,m1,m2,iqr,lma1,lma2,a1,a2
    integer :: j3
    real(8) :: v,v0,err,err0,QRtmp

    real(8),allocatable :: QRij_tmp(:,:,:)
    integer,allocatable :: icheck_grid(:,:,:)
    
    real(8) :: ctt(0:9),ett(0:9)
    real(8),parameter :: ep=1.d-8

    call write_border( 40, " prepQRijp102(start)" )

    ctt=0.0d0 ; ett=0.0d0
    call watch(ctt(6),ett(6))

    call allocateQaL

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = Igrid(2,1)-Igrid(1,1)+1
    ab2 = Igrid(2,2)-Igrid(1,2)+1
    ab3 = Igrid(2,3)-Igrid(1,3)+1

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    call watch(ctt(7),ett(7))

    r=maxval(Q_Rps)+maxval(Hgrid(1:3))+1.d-8
    call make_minimal_box(r,mm1,mm2,mm3,MMJJ_Q_0)
    
    call watch(ctt(8),ett(8))

    MMJJ_Q_0 = M_grid_ion

    nl3vmax = maxval( nl3v )

    if ( .not.allocated(JJ_tmp) ) allocate( JJ_tmp(6,MMJJ_Q_0,k1max,Natom) )
    JJ_tmp=0
    if ( .not.allocated(MJJ_tmp) ) allocate( MJJ_tmp(k1max,Natom) )
    MJJ_tmp=0
    if ( .not.allocated(QRij_tmp) ) allocate( QRij_tmp(MMJJ_Q_0,k1max,Natom) )
    QRij_tmp=0.0d0

    call watch(ctt(0),ett(0))

    c1     = 1.d0/Ngrid(1)
    c2     = 1.d0/Ngrid(2)
    c3     = 1.d0/Ngrid(3)
    maxerr = 0.d0

    do ia=1,Natom

       Rx = aa(1,1)*aa_atom(1,ia)+aa(1,2)*aa_atom(2,ia)+aa(1,3)*aa_atom(3,ia)
       Ry = aa(2,1)*aa_atom(1,ia)+aa(2,2)*aa_atom(2,ia)+aa(2,3)*aa_atom(3,ia)
       Rz = aa(3,1)*aa_atom(1,ia)+aa(3,2)*aa_atom(2,ia)+aa(3,3)*aa_atom(3,ia)

       ic1 = nint( aa_atom(1,ia)*Ngrid(1) )
       ic2 = nint( aa_atom(2,ia)*Ngrid(2) )
       ic3 = nint( aa_atom(3,ia)*Ngrid(3) )

       ik = ki_atom(ia)

       do ik1=1,N_k1(ik)

          ik2 = k1_to_k2(ik1,ik)
          ik3 = k1_to_k3(ik1,ik)

          Rps2 = Q_Rps(ik2,ik)*Q_Rps(ik2,ik)
          NRc  = Q_NRps(ik2,ik)

          j=0
          do i=1,M_grid_ion

             i1 = map_grid_ion(1,i)
             i2 = map_grid_ion(2,i)
             i3 = map_grid_ion(3,i)

             id1 = ic1 + i1
             id2 = ic2 + i2
             id3 = ic3 + i3

             k1  = id1/ML1 ; if ( id1<0 ) k1=(id1+1)/ML1-1
             k2  = id2/ML2 ; if ( id2<0 ) k2=(id2+1)/ML2-1
             k3  = id3/ML3 ; if ( id3<0 ) k3=(id3+1)/ML3-1

             i1_0  = id1-k1*ML1
             i2_0  = id2-k2*ML2
             i3_0  = id3-k3*ML3

             if ( Igrid(1,1) <= i1_0 .and. i1_0 <= Igrid(2,1) .and. &
                  Igrid(1,2) <= i2_0 .and. i2_0 <= Igrid(2,2) .and. &
                  Igrid(1,3) <= i3_0 .and. i3_0 <= Igrid(2,3) ) then

                d1 = id1*c1
                d2 = id2*c2
                d3 = id3*c3

                x  = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
                y  = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
                z  = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz

                r2 = x*x+y*y+z*z
                if ( r2 > Rps2+1.d-10 ) cycle

                j=j+1
                r = sqrt(r2)

                call get_qaL( r, x, y, z )
  
                do ll3=1,nl3v(ik2,ik)

                   L     = l3v(ll3,ik2,ik) - 1
                   v0    = 0.d0
                   err0  = 0.d0
                   QRtmp = 0.d0

                   if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 .or. &
                        abs(z)>1.d-14 .or. L==0 ) then

                      ir0 = 1
                      ir1 = NRc
111                   if ( ir1 - ir0 > 1 ) then
                         ir = ( ir0 + ir1 )/2
                         if ( rad1(ir,ik) > r ) then
                            ir1 = ir
                         else
                            ir0 = ir
                         end if
                         goto 111
                      end if
                      ir = ir0

                      if ( ir <= 2 ) then
                         v0  = qrL(2,ll3,ik2,ik)
                         if ( ir < 1 ) stop "prep_QRij_p102(0)"
                      else if ( ir <= NRc ) then
                         err0=1.d10
                         do mm=1,20
                            m1=max(1,ir-mm)
                            m2=min(ir+mm,NRc)
                            call polint(rad1(m1,ik),qrL(m1,ll3,ik2,ik) &
                                 ,m2-m1+1,r,v,err)
                            if ( abs(err) < err0 ) then
                               v0=v
                               err0=abs(err)
                               if ( err0 < ep ) exit
                            end if
                         end do
                      else
                         stop "stop@prep_QRij_p102(1)"
                      end if
                      maxerr=max(maxerr,err0)

                      QRtmp = v0/pi4 * dble( qaL(ik3,ll3)/(-zi)**L )

                   end if ! x,y,z

                   QRij_tmp(j,ik1,ia) = QRij_tmp(j,ik1,ia) + QRtmp

                end do ! ll3

                JJ_tmp(1,j,ik1,ia) = i1_0
                JJ_tmp(2,j,ik1,ia) = i2_0
                JJ_tmp(3,j,ik1,ia) = i3_0
                JJ_tmp(4,j,ik1,ia) = k1
                JJ_tmp(5,j,ik1,ia) = k2
                JJ_tmp(6,j,ik1,ia) = k3

             end if ! Igrid

          end do ! i ( 1 - MMJJ_0 )

          MJJ_tmp(ik1,ia) = j

       end do ! ik1

    end do ! ia

! ---

    allocate( icheck_grid(a1b:b1b,a2b:b2b,a3b:b3b) )
    icheck_grid=0

    MAXMJJ_Q = 0

    do ia=1,Natom
       do ik1=1,N_k1( ki_atom(ia) )

          j=0
          icheck_grid(:,:,:)=0
          do i=1,MJJ_tmp(ik1,ia)
             i1 = JJ_tmp(1,i,ik1,ia)
             i2 = JJ_tmp(2,i,ik1,ia)
             i3 = JJ_tmp(3,i,ik1,ia)
             if ( icheck_grid(i1,i2,i3) == 0 ) then
                j=j+1
                icheck_grid(i1,i2,i3) = j
             end if
          end do

          MAXMJJ_Q = max( MAXMJJ_Q, j )

       end do ! ik1
    end do ! ia

!===================================================================

    if ( allocated(QRij) ) deallocate( QRij )
    allocate( QRij(MAXMJJ_Q,N_nzqr) ) ; QRij=0.0d0

    if ( allocated(JJP_Q) ) deallocate( JJP_Q )
    allocate( JJP_Q(MAXMJJ_Q,N_nzqr) ) ; JJP_Q=0

    if ( allocated(MJJ_Q) ) deallocate( MJJ_Q )
    allocate( MJJ_Q(N_nzqr) ) ; MJJ_Q=0

!===================================================================

    do iqr=1,N_nzqr

       lma1 = nzqr_pair(iqr,1)
       lma2 = nzqr_pair(iqr,2)

       a1 = amap(lma1)
       a2 = amap(lma2)

       iorb1 = iorbmap(lma1)
       iorb2 = iorbmap(lma2)

       m1 = mmap(lma1)
       m2 = mmap(lma2)

       ik = ki_atom(a1)

       do ik1=1,N_k1(ik)
          if ( k1_to_iorb(1,ik1,ik) == iorb1 .and. &
               k1_to_iorb(2,ik1,ik) == iorb2 .and. &
               k1_to_m(1,ik1,ik) == m1 .and. k1_to_m(2,ik1,ik) == m2 ) exit
       end do

       if ( ik1 > N_k1(ik) ) cycle

       j=0
       icheck_grid(:,:,:)=0
       do i=1,MJJ_tmp(ik1,a1)

          i1 = JJ_tmp(1,i,ik1,a1)
          i2 = JJ_tmp(2,i,ik1,a1)
          i3 = JJ_tmp(3,i,ik1,a1)

          if ( icheck_grid(i1,i2,i3) == 0 ) then
             j=j+1
             icheck_grid(i1,i2,i3) = j
             QRij(j,iqr) = QRij_tmp(i,ik1,a1)
             JJP_Q(j,iqr) = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ML_0
          else
             j3 = icheck_grid(i1,i2,i3)
             QRij(j3,iqr) = QRij(j3,iqr) + QRij_tmp(i,ik1,a1)
          end if

       end do ! i

       MJJ_Q(iqr) = j

    end do ! iqr

    deallocate( icheck_grid )
    deallocate( MJJ_tmp )
    deallocate( JJ_tmp )

    call watch(ctt(5),ett(5))

    if ( disp_switch_parallel ) then
       write(*,*) "time(prepQRijp102_1)",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(prepQRijp102_2)",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(prepQRijp102_3)",ctt(3)-ctt(2),ett(3)-ett(2)
       write(*,*) "time(prepQRijp102_4)",ctt(4)-ctt(3),ett(4)-ett(3)
       write(*,*) "time(prepQRijp102_5)",ctt(5)-ctt(4),ett(5)-ett(4)
       write(*,*) "time(prepQRijp102_7)",ctt(7)-ctt(6),ett(7)-ett(6)
       write(*,*) "time(prepQRijp102_8)",ctt(8)-ctt(7),ett(8)-ett(7)
       write(*,*) "time(prepQRijp102_9)",ctt(0)-ctt(8),ett(0)-ett(8)
    end if

    call write_border( 40, " prepQRijp102(en)" )

  END SUBROUTINE prepQRijp102


  SUBROUTINE get_qaL ( r,x,y,z )
    implicit none

    real(8) :: r,x,y,z
    real(8),parameter :: sq3=sqrt(3.d0)
    real(8),parameter :: sq5=sqrt(5.d0)

!----------------------------
    qaL(:,:)   = z0

    qaL(1,1)   = z1
    qaL(3,1)   = z1
    qaL(6,1)   = z1
    qaL(10,1)  = z1
    qaL(15,1)  = z1
    qaL(21,1)  = z1
    qaL(28,1)  = z1
    qaL(36,1)  = z1
    qaL(45,1)  = z1

!--      
    if ( r <= 1.0D-10 ) return
!--      
    x = x/r 
    y = y/r 
    z = z/r 

    qaL(1,1)=1.d0
    qaL(2,1)=sq3*zi*y
    qaL(3,1)=1.d0
    qaL(3,2)=(3.d0*x**2-3.d0*y**2+3.d0*z**2-1.d0)/2.d0
    qaL(4,1)=-sq3*zi*z
    qaL(5,1)=0.d0
    qaL(5,2)=3.d0*y*z
    qaL(6,1)=1.d0
    qaL(6,2)=-(3.d0*z**2-1.d0)
    qaL(7,1)=sq3*zi*x
    qaL(8,1)=0.d0
    qaL(8,2)=-3.d0*x*y
    qaL(9,1)=0.d0
    qaL(9,2)=3.d0*x*z
    qaL(10,1)=1.d0
    qaL(10,2)=-(3.d0*x**2-3.d0*y**2-3.d0*z**2+1.d0)/2.d0
    qaL(11,1)=-sqrt(15.d0)*x*y
    qaL(12,1)=3.d0*sq5*zi*x/5.d0
    qaL(12,2)=(   3.d0*sq5*zi*x * (5.d0*x**2-15.d0*y**2+5.d0*z**2-1.d0) )/20.d0
    qaL(13,1)=0.d0
    qaL(13,2)=3.d0*sq5*zi*x*y*z
    qaL(14,1)=3.d0*sq5*zi*y/5.d0
    qaL(14,2)=-(3.d0*sq5*zi*y*(15.d0*x**2 -5.d0*y**2-5.d0*z**2+1.d0))/(20.d0)
    qaL(15,1)=1.d0
    qaL(15,2)=(5.d0*(3.d0*z**2-1.d0))/7.d0
    qaL(15,3)=-(3.d0*(35.d0*x**4-210.d0*x**2*y**2 +35.d0*y**4-35.d0*z**4+30.d0*z**2-3.d0))/56.d0
    qaL(16,1)=3.d0*sqrt(5.d0/3.0d0)*y*z
    qaL(17,1)=-(3.d0*sq5*zi*z)/5.d0
    qaL(17,2)=-(3.d0*sq5*zi*z *(5.d0*x**2-5.d0*y**2+5.d0*z**2-3.d0))/10.d0
    qaL(18,1)=(3.d0*sq5*zi*y)/5.d0
    qaL(18,2)=-(3.d0*zi*y*(5.d0*z**2-1.d0))/sq5
    qaL(19,1)=0.d0
    qaL(19,2)=3.d0*sq5*zi*x*y*z
    qaL(20,1)=0.d0
    qaL(20,2)=(15.d0*x*z)/7.d0
    qaL(20,3)=(15.d0*x*z *(7.d0*x**2-21.d0*y**2+7.d0*z**2-3.d0))/28.d0
    qaL(21,1)=1.d0
    qaL(21,2)=(5.d0*(3.d0*x**2-3.d0*y**2-3.d0*z**2+1.d0))/14.d0
    qaL(21,3)=-(3.d0*(35.d0*x**2*z**2-5.d0*x**2-35.d0*y**2*z**2 +5.d0*y**2+35.d0*z**4-30.d0*z**2+3.d0))/14.d0
    qaL(22,1)=-(sq5*(3.d0*z**2-1.d0))/2.d0
    qaL(23,1)=-(sqrt(15.d0)*zi*y)/5.d0
    qaL(23,2)=-1.5d0*sqrt(3.d0/5.d0)*zi*y*(5.d0*z**2-1.d0)
    qaL(24,1)=-(2.d0*sqrt(15.d0)*zi*z)/5.d0
    qaL(24,2)=9.d0*zi*z*(5.d0*z**2-3.d0)/(2.d0*sqrt(15.d0))
    qaL(25,1)=-(sqrt(15.d0)*zi*x)/5.d0
    qaL(25,2)=-1.5d0*sqrt(3.d0/5.d0)*zi*x*(5.d0*z**2-1.d0)
    qaL(26,1)=0.d0
    qaL(26,2)=10.d0*sq3*x*y/7.d0
    qaL(26,3)=15.d0*sq3*x*y*(7.d0*z**2-1.d0)/14.d0
    qaL(27,1)=0.d0
    qaL(27,2)=(15.d0*y*z)/(7.d0*sq3)
    qaL(27,3)=-15.d0*sq3*y*z*(7.d0*z**2-3.d0)/14.d0
    qaL(28,1)=1.d0
    qaL(28,2)=-(5.d0*(3.d0*z**2-1.d0))/7.d0
    qaL(28,3)=(9.d0*(35.d0*z**4-30.d0*z**2+3.d0))/28.d0
    qaL(29,1)=(3.d0*sqrt(5.d0/3.d0)*x*z)
    qaL(30,1)=0.d0
    qaL(30,2)=3.d0*sq5*zi*x*y*z
    qaL(31,1)=(3.d0*sq5*zi*x)/5.d0
    qaL(31,2)=-(3.d0*zi*x*(5.d0*z**2-1.d0))/sq5
    qaL(32,1)=-(3.d0*sq5*zi*z)/5.d0
    qaL(32,2)=0.3d0*sq5*zi*z*(5.d0*x**2-5.d0*y**2-5.d0*z**2+3.d0)
    qaL(33,1)=0.d0
    qaL(33,2)=(15.d0*y*z)/7.d0
    qaL(33,3)=-15.d0*y*z *(21.d0*x**2-7.d0*y**2-7.d0*z**2+3.d0)/28.d0
    qaL(34,1)=0.d0
    qaL(34,2)=-(15.d0*x*y)/7.d0
    qaL(34,3)=(15.d0*x*y*(7.d0*z**2-1.d0))/7.d0
    qaL(35,1)=0.d0
    qaL(35,2)=(15.d0*x*z)/(7.d0*sq3)
    qaL(35,3)=-15.d0*sq3*x*z*(7.d0*z**2-3.d0)/14.d0
    qaL(36,1)=1.d0
    qaL(36,2)=-5.d0*(3.d0*x**2-3.d0*y**2+3.d0*z**2-1.d0)/14.d0
    qaL(36,3)=(3.d0*(35.d0*x**2*z**2-5.d0*x**2-35.d0*y**2*z**2 +5.d0*y**2-35.d0*z**4+30.d0*z**2-3.d0))/14.d0
    qaL(37,1)=-(3.d0*sq5*(x**2-y**2))/(2.d0*sq3)
    qaL(38,1)=-(3.d0*sqrt(10.d0)*zi*y)/(5.d0*sqrt(2.d0))
    qaL(38,2)=-0.15d0*sq5*zi*y *(15.d0*x**2-5.d0*y**2+5.d0*z**2-1.d0)
    qaL(39,1)=0.d0
    qaL(39,2)=1.5d0*sq5*zi*z*(x**2-y**2)
    qaL(40,1)=(3.d0*sq5*zi*x)/5.d0
    qaL(40,2)=-0.15d0*sq5*zi*x  *(5.d0*x**2-15.d0*y**2-5.d0*z**2+1.d0)
    qaL(41,1)=0.d0
    qaL(41,2)=0.d0
    qaL(41,3)=7.5d0*x*y*(x**2-y**2)
    qaL(42,1)=0.d0
    qaL(42,2)=-(15.d0*y*z)/7.d0
    qaL(42,3)=-15.d0*y*z *(21.d0*x**2-7.d0*y**2+7.d0*z**2-3.d0)/28.d0
    qaL(43,1)=0.d0
    qaL(43,2)=(15.d0*(x**2-y**2))/(7.d0*sq3)
    qaL(43,3)=15.d0*sq3 *(7.d0*x**2*z**2-x**2-7.d0*y**2*z**2+y**2)/28.d0
    qaL(44,1)=0.d0
    qaL(44,2)=(15.d0*x*z)/7.d0
    qaL(44,3)=-15.d0*x*z *(7.d0*x**2-21.d0*y**2-7.d0*z**2+3.d0)/28.d0
    qaL(45,1)=1.d0
    qaL(45,2)=(5.d0*(3.d0*z**2-1.d0))/7.d0
    qaL(45,3)=3.d0*(35.d0*x**4-210.d0*x**2*y**2 +35.d0*y**4+35.d0*z**4-30.d0*z**2+3.d0)/56.d0

!--    THIS will be done when qaL is used.
!      qaL(:,:)=qaL(:,:)/( (4.d0*Pi)*(-zi)**L )

    return

  END SUBROUTINE get_qaL

END MODULE PSQRijPrep
