MODULE VarPSMemberG
  implicit none

  integer,allocatable :: nl3v(:,:)
  integer,allocatable :: l3v(:,:,:)
  
  real(8),allocatable :: qrL(:,:,:,:)
  real(8),allocatable :: dqrL(:,:,:,:,:)

  complex(8),allocatable :: QG(:,:,:)

  real(8),allocatable :: QRij(:,:)


  real(8),allocatable :: ddi(:,:,:,:)

  real(8),allocatable :: qqr(:,:,:,:)
  real(8),allocatable :: qqc(:,:,:,:)

  integer :: k1max,k2max,k3max
  
  integer,allocatable :: k1_to_k2(:,:)
  integer,allocatable :: k1_to_k3(:,:)
  integer,allocatable :: k1_to_iorb(:,:,:)
  integer,allocatable :: N_k1(:)

  integer,allocatable :: k2_to_iorb(:,:,:)
  integer,allocatable :: N_k2(:)

  integer,allocatable :: icheck_k2(:)

  integer,allocatable :: Q_NRps(:,:)
  real(8),allocatable :: Q_Rps(:,:)

  integer,allocatable :: npq(:)

  integer :: max_Rref=0,max_Lref=0,max_k2=0,max_qgrd=0


CONTAINS

!------------------------------------------
  SUBROUTINE allocateKtoK( k1max,k2max,nki,Rrefmax,Lrefmax )
    implicit none
    integer,intent(IN) :: k1max,k2max,nki,Rrefmax,Lrefmax
    
    if ( allocated(k1_to_k2) ) then
      call deallocateKtoK
    end if

    allocate( k1_to_k2(k1max,nki)     ) ; k1_to_k2(:,:)=0
    allocate( k1_to_k3(k1max,nki)     ) ; k1_to_k3(:,:)=0
    allocate( k1_to_iorb(2,k1max,nki) ) ; k1_to_iorb(:,:,:)=0
    allocate( N_k1(nki)               ) ; N_k1(:)=0

    allocate( k2_to_iorb(2,k2max,nki)   ) ; k2_to_iorb(:,:,:)=0
    allocate( N_k2(nki)               ) ; N_k2(:)=0

    allocate( icheck_k2(nki)          ) ; icheck_k2(:)=0

    allocate( qqc(Rrefmax,Rrefmax,Lrefmax,nki) ) ; qqc(:,:,:,:)=0.d0

    return
  END SUBROUTINE allocateKtoK

!------------------------------------------
  SUBROUTINE deallocateKtoK
    implicit none
    
    deallocate( k1_to_k2 )
    deallocate( k1_to_k3 )
    deallocate( k1_to_iorb )
    deallocate( N_k1 )

    deallocate( k2_to_iorb )
    deallocate( N_k2 )

    deallocate( icheck_k2 )

    deallocate( qqc )

    return
  END SUBROUTINE deallocateKtoK

!------------------------------------------

  SUBROUTINE allocateQRps( k2max,nki )
    implicit none
    integer,intent(IN) :: k2max,nki

    if ( allocated(Q_NRps) ) then
      call deallocateQRps
    end if

    allocate( Q_NRps(k2max,nki) ) ; Q_NRps(:,:)=0
    allocate( Q_Rps(k2max,nki)  ) ; Q_Rps(:,:)=0.d0

    return
  END SUBROUTINE allocateQRps

!------------------------------------------
  SUBROUTINE deallocateQRps
    implicit none

    deallocate( Q_NRps )
    deallocate( Q_Rps  )

    return
  END SUBROUTINE deallocateQRps

!------------------------------------------
  SUBROUTINE allocatePSG( n_l,n_r,n_k,n_g,nki )
    implicit none
    integer,intent(IN) :: n_l,n_r,n_k,n_g,nki
    real(8),allocatable :: ddi_(:,:,:,:),qqr_(:,:,:,:)
    integer,allocatable :: nl3v_(:,:),l3v_(:,:,:)
    real(8),allocatable :: qrL_(:,:,:,:)
    integer :: ml,mr,mk,mg

    if ( max_Rref==0 .or. max_Lref==0 .or. max_k2==0 ) then
        allocate( npq(nki) ) ; npq=0
        allocate( ddi(n_r,n_r,n_l,nki) ) ; ddi=0.d0
        allocate( qqr(n_r,n_r,n_l,nki) ) ; qqr=0.d0
        allocate( qqc(n_r,n_r,n_l,nki) ) ; qqc=0.d0
        allocate( nl3v(n_k,nki) ) ; nl3v=0
        allocate( l3v(n_l,n_k,nki) ) ; l3v=0
        allocate( qrL(max_psgrd,n_l,n_k,nki) ) ; qrL=0.d0
        max_Rref=n_r
        max_Lref=n_l
        max_k2=n_k
        max_qgrd=n_gr
        return
    end if
    ml=max(max_Lref,n_l)
    mr=max(max_Rref,n_r)
    mk=max(max_k2,n_k)
    mg=max(max_qgrd,n_g)
    if ( max_Rref<mr ) then
        allocate( ddi_(mr,mr,ml,nki) ) ; ddi_=0.d0
        allocate( qqr_(mr,mr,ml,nki) ) ; qqr_=0.d0
        ddi_(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)=ddi(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)
        qqr_(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)=qqr(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)
        deallocate( ddi,qqr,qqc )
        allocate( ddi(mr,mr,ml,nki) ) ; ddi=0.d0
        allocate( qqr(mr,mr,ml,nki) ) ; qqr=0.d0
        allocate( qqc(mr,mr,ml,nki) ) ; qqc=0.d0
        ddi(:,:,:,:)=ddi_(:,:,:,:)
        qqr(:,:,:,:)=qqr_(:,:,:,:)
        deallocate( ddi_,qqr_ )
    end if
    if ( max_k2<mk ) then
        allocate( nl3v_(mk,nki) ) ; nl3v_=0
        allocate( l3v_(ml,mk,nki) ) ; l3v_=0
        allocate( qrL_(mg,ml,mk,nki) ) ; qrL_=0.d0
        nl3v_(1:max_k2,1:nki)=nl3v(1:max_k2,1:nki)
        l3v_(1:max_Lref,1:max_k2,1:nki)=l3v(1:max_Lref,1:max_k2,1:nki)
        qrL_(1:max_qgrd,1:max_Lref,max_k2,1:nki)=qrL(1:max_qgrd,1:max_Lref,max_k2,1:nki)
        deallocate( nl3v,l3v,qrL )
        allocate( nl3v(mk,nki) ) ; nl3v=0
        allocate( l3v(ml,mk,nki) ) ; l3v=0
        allocate( qrL(mg,ml,mk,nki) ) ; qrL=0.d0
        nl3v(:,:)=nl3v_(:,:)
        l3v(:,:,:)=l3v_(:,:,:)
        qrL(:,:,:,:)=qrL_(:,:,:,:)
        deallocate( nl3v_,l3v_,qrL_ )
    end if
    if ( max_Lref<ml ) then
        if ( max_k2>=mk ) then
            allocate( l3v_(ml,mk,nki) ) ; l3v_=0
            allocate( qrL_(mg,ml,mk,nki) ) ; qrL_=0.d0
            l3v_(1:max_Lref,1:max_k2,1:nki)=l3v(1:max_Lref,1:max_k2,1:nki)
            qrL_(1:max_qgrd,1:max_Lref,max_k2,1:nki)=qrL(1:max_qgrd,1:max_Lref,max_k2,1:nki)
            deallocate( l3v,qrL )
            allocate( l3v(ml,mk,nki) ) ; l3v=0
            allocate( qrL(mg,ml,mk,nki) ) ; qrL=0.d0
            l3v(:,:,:)=l3v_(:,:,:)
            qrL(:,:,:,:)=qrL_(:,:,:,:)
            deallocate( l3v_,qrL_ )
        end if
        if ( max_Rref>=mr ) then
            allocate( ddi_(mr,mr,ml,nki) ) ; ddi_=0.d0
            allocate( qqr_(mr,mr,ml,nki) ) ; qqr_=0.d0
            ddi_(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)=ddi(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)
            qqr_(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)=qqr(1:max_Rref,1:max_Rref,1:max_Lref,1:nki)
            deallocate( ddi,qqr,qqc )
            allocate( ddi(mr,mr,ml,nki) ) ; ddi=0.d0
            allocate( qqr(mr,mr,ml,nki) ) ; qqr=0.d0
            allocate( qqc(mr,mr,ml,nki) ) ; qqc=0.d0
            ddi(:,:,:,:)=ddi_(:,:,:,:)
            qqr(:,:,:,:)=qqr_(:,:,:,:)
            deallocate( ddi_,qqr_ )
        end if
    end if
    if ( max_qgrd<mg ) then
        if ( max_Lref>=ml .and. max_k2>=mk ) then
            allocate( qrL_(mk,ml,mk,nki) ) ; qrL_=0.d0
            qrL_(1:max_qgrd,1:max_Lref,max_k2,1:nki)=qrL(1:max_qgrd,1:max_Lref,max_k2,1:nki)
            deallocate( qrL )
            allocate( qrL(mg,ml,mk,nki) ) ; qrL=0.d0
            qrL(:,:,:,:)=qrL_(:,:,:,:)
            deallocate( qrL_ )
        end if
    end if
    max_Lref=ml
    max_Rref=mr
    max_k2=mk
    max_qgrd=mg

    return
  END SUBROUTINE allocatePSG

!------------------------------------------
  SUBROUTINE deallocatePSG
    implicit none

    deallocate( nl3v )
    deallocate( l3v  )
    deallocate( qrad )
    deallocate( qrL  )

    return
  END SUBROUTINE deallocatePSG
END MODULE VarPSMemberG
