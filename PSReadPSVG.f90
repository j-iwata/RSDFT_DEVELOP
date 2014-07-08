MODULE PSReadPSVG
    use VarPSMember
    use VarPSMemberG, only: npq,nl3v,l3v,ncf,nrin,rin,coe,qrad,qrL
    implicit none

CONTAINS

    SUBROUTINE readPSVG( unit_ps,ik )
        implicit none
        integer,intent(IN) :: unit_ps,ik
        character(9) :: cbuf9
        integer :: i,j,l1,l2,i1,i2,ic
        integer :: k2,ll1,ll2,ll3,ii1,ii2
        real(8) :: dif,r2
        integer,allocatable :: nl3v_l(:),l3v_l(:,:),ncf(:,:),nrin(:,:)
        real(8),allocatable :: rin(:,:),coe(:,:,:)
        real(8),allocatable :: qrL_l(:,:,:),qrad_l(:,:,:)
        integer :: Rrefmax,Lrefmax,lpsmax,npqmax,ncfmax,nsmpl
        
        Lrefmax=nl(ik)
        Rrefmax=max(nr(:,ik))
        lpsmax=Lrefmax*Rrefmax
        npqmax=lpsmax*(lpsmax+1)/2
        ncfmax=10
        nsmpl=NRps(1,ik)
        
        allocate( nl3v_l(npqmax)       ) ; nl3v_l=0
        allocate( l3v_l(lpsmax,npqmax) ) ; l3v_l =0
        allocate( ncf(lpsmax,npqmax)        ) ; ncf =0
        allocate( nrin(lpsmax,npqmax)       ) ; nrin=0
        allocate( rin(lpsmax,npqmax)        ) ; rin =0.d0
        allocate( coe(ncfmax,lpsmax,npqmax) ) ; coe   =0.d0
        allocate( qrL_l(nsmpl+1,lpsmax,npqmax)  ) ; qrL_l =0.d0
        allocate( qrad_l(nsmpl+1,lpsmax,npqmax) ) ; qrad_l=0.d0
        
        do j=1,10000
            read(unit_ps,'(A)') cbuf9
            if (cbuf9=='#### DATA' ) exit
        end do
        read(unit_ps,*) npq_l
        
        k2=0
        do l1=1,nl(ik)
            do i1=1,nr(l1,ik)
                do l2=1,l1
                    do i2=1,nr(l2,ik)
                        if (l1==l2 .and. i2>i1) cycle
                        k2=k2+1
                        read(unit_ps,*) ll1,ii1,ll2,ii2
                        if ((ll1/=l1) .or. (ll2/=l2) .or. (ii1/=i1) .or. (ii2/=i2)) then
                            stop 'ERROR in pseudization data'
                        end if
                        read(unit_ps,*) nl3v_l(k2)
                        do ll3=1,nl3v_l(k2)
                            read(unit_ps,*) l3v_l(ll3,k2),ncf(ll3,k2),nrin(ll3,k2),rin(ll3,k2)
                            if (ncf(ll3,k2)>0) then
                                read(unit_ps,*) coe(1:ncf(ll3,k2),ll3,k2)
                                dif=abs(rr(nrin(ll3-k2))-rin(ll3,k2))
                                if (dif>1.d-8) then
                                    stop 'ERROR in PSV read space cutoff'
                                end if
                            end if
                            qrL_l(:,ll3,k2)=0.d0
                            do ic=ncf(ll3,k2),1,-1
                                do i=1,nrin(ll3,k2)
                                    r2=rr(i)*rr(i)
                                    qrL_l(i,ll3,k2)=qrL_l(i,ll3,k2)*r2+coe(ic,ll3,k2)
                                end do
                            end do
                            do i=1,nrin(ll3,k2)
                                qrL_l(i,ll3,k2)=qrL_l(i,ll3,k2)*rr(i)**(l3v_l(ll3,k2)+1)
                                qrad_l(i,ll3,k2)=psi(i,i1,l1)*psi(i,i2,l2)-phi(i,i1,l1)*phi(i,i2,l2)
                            end do
                            do i=nrin(ll3,k2)+1,nsmpl
                                qrL_l(i,ll3,k2)=psi(i,i1,l1)*psi(i,i2,l2)-phi(i,i1,l1)*phi(i,i2,l2)
                                qrad_l(i,ll3,k2)=qrL_l(i,ll3,k2)
                            end do
                        end do ! ll3
                    end do ! i2
                end do ! l2
            end do ! i1
        end do ! l1
        
        if (npq_l/=k2) stop 'ERROR npq/=k2'
        npqmax=npq_l

    END SUBROUTINE readPSVG
END MODULE PSReadPSVG
