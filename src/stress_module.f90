MODULE stress_module
  use array_bound_module, only: ML, ML_0, ML_1, MSP, MSP_0, MSP_1, &
       MB, MB_0, MB_1, MBZ, MBZ_0, MBZ_1
  use aa_module, only: aa, Va
  use bb_module, only: bb
  use bz_module, only: kbb
  use rgrid_module, only: Igrid, dV, Ngrid, Hgrid
  use ggrid_module, only: MG_0,MG_1, GG, LLG, MGL, NGgrid, NMGL, &
       construct_Ggrid, destruct_Ggrid, allgatherv_Ggrid
  use atom_module, only: aa_atom, ki_atom, Natom, Nelement
  use density_module, only: rho
  use wf_module, only: unk,occ
  use total_energy_module, only: Ekin, Eion, Enlc, cnst,calc_total_energy
  use hartree_variables, only: E_hartree
  use eion_module, only: Eewald
  use xc_module, only: Exc, Vxc
  use parallel_module, only: comm_grid, MB_d
  use fft_module
  use omp_lib

  implicit none

  PRIVATE
  PUBLIC :: test_stress
  PUBLIC :: calc_stress

  include 'mpif.h'

  real(8), parameter :: delta(3,3) = &
       reshape( (/ &
       1.0d0, 0.0d0, 0.0d0, &
       0.0d0, 1.0d0, 0.0d0, &
       0.0d0, 0.0d0, 1.0d0 /), (/3,3/) )

  real(8), parameter :: M_PI       = 3.14159265358979323846d0
  real(8), parameter :: M_2PI      = M_PI*2.0d0
  real(8), parameter :: M_4PI      = M_PI*4.0d0
  real(8), parameter :: M_SQRTPI   = 1.77245385090551602792981d0
  real(8), parameter :: M_2_SQRTPI = 1.12837916709551257390d0
  complex(8), parameter :: M_CI    = dcmplx(0.0d0,1.0d0)

CONTAINS

  SUBROUTINE test_stress(systype)
    implicit none
    integer,intent(IN) :: systype
    real(8) :: stress(3,3)

    call write_border( 60, " test_stress(start)" )

    select case(systype)
    case default
       call calc_stress(stress)

    case(1)
       ! nothing to do
    end select

    call write_border( 60, " test_stress(end)" )

  END SUBROUTINE test_stress


  SUBROUTINE calc_stress( stress_tot )
    real(8), intent(out) :: stress_tot(3,3)

    real(8) :: stress_kin(3,3)
    real(8) :: stress_ion(3,3)
    real(8) :: stress_har(3,3)
    real(8) :: stress_xc (3,3)
    real(8) :: stress_nlc(3,3)
    real(8) :: stress_ewa(3,3)
    real(8) :: stress_cnt(3,3)
    real(8) :: Etot

    ! calc Ekin, Eion, Enlc, cnst
    call calc_total_energy( .true., Etot )

    ! call each stress
    call calc_stress_kin( stress_kin )
    call calc_stress_ion( stress_ion )
    call calc_stress_har( stress_har )
    call calc_stress_xc ( stress_xc  )
    call calc_stress_nlc( stress_nlc )
    call calc_stress_ewa( stress_ewa )
    call calc_stress_cnt( stress_cnt )

    ! calc total stress
    stress_tot(:,:) = &
         + stress_kin(:,:) &
         + stress_ion(:,:) &
         + stress_har(:,:) &
         + stress_xc(:,:)  &
         + stress_nlc(:,:) &
         + stress_ewa(:,:) &
         + stress_cnt(:,:)

    ! calc total energy corresponding to this stress
    Etot = &
         + Ekin &
         + Eion &
         + E_hartree &
         + Exc &
         + Enlc &
         + Eewald &
         + cnst

    write(*,*) "stress tensor"
    write(*,"(3f12.6)")  real(stress_tot(:,1))
    write(*,"(3f12.6)")  real(stress_tot(:,2))
    write(*,"(3f12.6)")  real(stress_tot(:,3))
    
!!$    write(*,"(a13,f12.6)") "# stress_kin", &
!!$         real(stress_kin(1,1)+stress_kin(2,2)+stress_kin(3,3))
!!$    write(*,"(a13,f12.6)") "# stress_ion", &
!!$         real(stress_ion(1,1)+stress_ion(2,2)+stress_ion(3,3))
!!$    write(*,"(a13,f12.6)") "# stress_har", &
!!$         real(stress_har(1,1)+stress_har(2,2)+stress_har(3,3))
!!$    write(*,"(a13,f12.6)") "# stress_xc ", &
!!$         real(stress_xc (1,1)+stress_xc (2,2)+stress_xc (3,3))
!!$    write(*,"(a13,f12.6)") "# stress_nlc", &
!!$         real(stress_nlc(1,1)+stress_nlc(2,2)+stress_nlc(3,3))
!!$    write(*,"(a13,f12.6)") "# stress_ewa", &
!!$         real(stress_ewa(1,1)+stress_ewa(2,2)+stress_ewa(3,3))
!!$    write(*,"(a13,f12.6)") "# stress_cnt", &
!!$         real(stress_cnt(1,1)+stress_cnt(2,2)+stress_cnt(3,3))
!!$    write(*,"(a13,f12.6)") "# stress_tot", &
!!$         real(stress_tot(1,1)+stress_tot(2,2)+stress_tot(3,3))

!!$    write(*,"(4f12.6,x,a)") real(aa(1,1)), real(Va), real(Etot), &
!!$         real(stress_tot(1,1)+stress_tot(2,2)+stress_tot(3,3)), "# DEBUG STRESS"

  end SUBROUTINE calc_stress


  SUBROUTINE calc_stress_kin( stress_kin )
    real(8), intent(out) :: stress_kin(3,3)

    integer :: i,n,k,s,nb1,nb2,info,la,mu
    complex(8),allocatable :: work(:,:)

    allocate( work(ML_0:ML_1,MB_d) )
    work = 0.0d0

    stress_kin(:,:) = 0.0d0
    do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
          do n=MB_0,MB_1,MB_d
             nb1=n
             nb2=min(nb1+MB_d-1,MB_1)

             do la=1, 3
                do mu=1, 3
                   work = 0.0d0
!$OMP parallel
                   call op_kinetic_sol_eps(k,unk(ML_0,n,k,s),work, &
                        ML_0,ML_1,nb1,nb2,la,mu)
!$OMP end parallel
                   do i=nb1,nb2
                      stress_kin(la,mu) = stress_kin(la,mu) &
                           + occ(i,k,s) &
                           * sum( real(conjg(unk(:,i,k,s))*work(:,i-nb1+1)) )*dV
                   end do
                end do ! mu
             end do ! la
          end do ! n
       end do ! k
    end do ! s

    stress_kin(:,:) = stress_kin(:,:)*(-1.0d0/Va)

    deallocate( work )

    call mpi_allreduce(MPI_IN_PLACE,stress_kin,3*3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,info)

    return
  contains
    SUBROUTINE op_kinetic_sol_eps(k,tpsi,htpsi,n1,n2,ib1,ib2,la,mu)
      use kinetic_variables, only: Md, wk
      use bc_module, only: www, bcset_1, bcset_3
      use fd_module, only: get_coef_lapla_fd, get_coef_nabla_fd
      use omp_variables, only: Igrid_omp

      integer,intent(IN)       :: k,n1,n2,ib1,ib2,la,mu
      complex(8),intent(IN)    ::  tpsi(n1:n2,ib1:ib2)
      complex(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)

      logical :: flag_nab, flag_n12, flag_n23, flag_n31
      real(8) :: coef_lap0, ggg(6), const_k2
      real(8),allocatable :: coef_lap(:,:), coef_nab(:,:)
      complex(8),allocatable :: zcoef_kin(:,:)
      real(8), allocatable :: nab(:), lap(:)

      integer :: i,ib,i1,i2,i3,nb,m,n,j,p,mm
      integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12
      integer :: a1b_omp,b1b_omp,a2b_omp,b2b_omp,a3b_omp,b3b_omp,n1_omp,n2_omp
      integer :: ib1_omp,ib2_omp,nb_omp
      real(8) :: c,d,a(3),kk(3),kG(3)

      a(1) = sqrt(sum(aa(:,1)**2))
      a(2) = sqrt(sum(aa(:,2)**2))
      a(3) = sqrt(sum(aa(:,3)**2))
      ggg(1) = - a(1)*a(1)*(bb(la,1)*bb(mu,1)+bb(mu,1)*bb(la,1))/(4*M_PI**2)
      ggg(2) = - a(2)*a(2)*(bb(la,2)*bb(mu,2)+bb(mu,2)*bb(la,2))/(4*M_PI**2)
      ggg(3) = - a(3)*a(3)*(bb(la,3)*bb(mu,3)+bb(mu,3)*bb(la,3))/(4*M_PI**2)
      ggg(4) = - a(1)*a(2)*(bb(la,1)*bb(mu,2)+bb(mu,1)*bb(la,2))/(4*M_PI**2)
      ggg(5) = - a(2)*a(3)*(bb(la,2)*bb(mu,3)+bb(mu,2)*bb(la,3))/(4*M_PI**2)
      ggg(6) = - a(3)*a(1)*(bb(la,3)*bb(mu,1)+bb(mu,3)*bb(la,1))/(4*M_PI**2)

      allocate( lap(-Md:Md), nab(-Md:Md) ) 
      call get_coef_lapla_fd(Md,lap)
      call get_coef_nabla_fd(Md,nab)

      allocate( coef_lap(3,Md), coef_nab(3,Md) )
      do n=1,Md
         coef_lap(:,n) = -0.5d0*ggg(1:3)*lap(n)/Hgrid(:)**2
         coef_nab(:,n) = nab(n)/Hgrid(:)
      end do
      coef_lap0 = -0.5d0*lap(0)*sum(ggg(1:3)/Hgrid(:)**2)

      deallocate( nab, lap )

      kk(:) = matmul( bb(:,:),kbb(:,k) )
      kG(:) = - kk(la)*bb(mu,:) - kk(mu)*bb(la,:)

      allocate( zcoef_kin(3,-Md:Md) )
      do n=1,Md
         zcoef_kin(:,-n) = coef_lap(:,n) &
              + M_CI*kG(:)/M_2PI*a(:)*coef_nab(:,n)
         zcoef_kin(:, n) = coef_lap(:,n) &
              - M_CI*kG(:)/M_2PI*a(:)*coef_nab(:,n)
      end do
      const_k2 = - kk(la)*kk(mu)

      flag_n12 = ggg(4) /= 0.d0
      flag_n23 = ggg(5) /= 0.d0
      flag_n31 = ggg(6) /= 0.d0
      flag_nab = ANY( kG(:)/=0.d0 )

      a1b = Igrid(1,1)
      b1b = Igrid(2,1)
      a2b = Igrid(1,2)
      b2b = Igrid(2,2)
      a3b = Igrid(1,3)
      b3b = Igrid(2,3)
      ab1 = (b1b-a1b+1)
      ab12= (b1b-a1b+1)*(b2b-a2b+1)

      nb = ib2-ib1+1

      mm=0
!$    mm=omp_get_thread_num()
      n1_omp = Igrid_omp(1,0,mm)
      n2_omp = Igrid_omp(2,0,mm)

      a1b_omp = Igrid_omp(1,1,mm)
      b1b_omp = Igrid_omp(2,1,mm)
      a2b_omp = Igrid_omp(1,2,mm)
      b2b_omp = Igrid_omp(2,2,mm)
      a3b_omp = Igrid_omp(1,3,mm)
      b3b_omp = Igrid_omp(2,3,mm)

      c=coef_lap0+const_k2
      do ib=ib1,ib2
         do i=n1_omp,n2_omp
            htpsi(i,ib) = c*tpsi(i,ib)
         end do
      end do

      do ib=ib1,ib2
         do i3=a3b_omp,b3b_omp
            do i2=a2b_omp,b2b_omp
               do i1=a1b_omp,b1b_omp
                  j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                  www(i1,i2,i3,ib-ib1+1) = tpsi(j,ib)
               end do
            end do
         end do
      end do

!$OMP barrier
      call bcset_3(1,nb,Md,0)
!$OMP barrier

      if ( flag_nab ) then

         do ib=ib1,ib2

            p = ib-ib1+1

            do m=1,Md

               do i3=a3b_omp,b3b_omp
                  do i2=a2b_omp,b2b_omp
                     do i1=a1b_omp,b1b_omp
                        j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                        htpsi(j,ib)=htpsi(j,ib) &
                             +      zcoef_kin(1,m) *www(i1+m,i2,i3,p) &
                             +conjg(zcoef_kin(1,m))*www(i1-m,i2,i3,p) &
                             +      zcoef_kin(2,m) *www(i1,i2+m,i3,p) &
                             +conjg(zcoef_kin(2,m))*www(i1,i2-m,i3,p) &
                             +      zcoef_kin(3,m) *www(i1,i2,i3+m,p) &
                             +conjg(zcoef_kin(3,m))*www(i1,i2,i3-m,p)   
                     end do
                  end do
               end do

            end do ! m

         end do ! ib

      else

         do ib=ib1,ib2

            p=ib-ib1+1

            do m=1,Md

               do i3=a3b_omp,b3b_omp
                  do i2=a2b_omp,b2b_omp
                     do i1=a1b_omp,b1b_omp
                        j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                        htpsi(j,ib)=htpsi(j,ib) &
                             +coef_lap(1,m)*( www(i1-m,i2,i3,p)+www(i1+m,i2,i3,p) ) &
                             +coef_lap(2,m)*( www(i1,i2-m,i3,p)+www(i1,i2+m,i3,p) ) &
                             +coef_lap(3,m)*( www(i1,i2,i3-m,p)+www(i1,i2,i3+m,p) )
                     end do
                  end do
               end do

            end do ! m

         end do ! ib

      end if

!$OMP barrier
      if ( flag_n12 .or. flag_n23 .or. flag_n31 ) then

!$OMP workshare
         wk=www
!$OMP end workshare

         if ( flag_n12 ) then

            do n=1,nb
               do i3=a3b_omp,b3b_omp
                  do i2=a2b_omp,b2b_omp
                     do i1=a1b_omp,b1b_omp
                        www(i1,i2,i3,n) = 0.0d0
                     end do
                  end do
               end do
            end do
            do n=1,nb
               do m=1,Md
                  d=coef_nab(1,m)
                  do i3=a3b_omp,b3b_omp
                     do i2=a2b_omp,b2b_omp
                        do i1=a1b_omp,b1b_omp
                           www(i1,i2,i3,n)=www(i1,i2,i3,n) &
                                -d*( wk(i1-m,i2,i3,n)-wk(i1+m,i2,i3,n) )
                        end do
                     end do
                  end do
               end do
            end do

!$OMP barrier
            call bcset_1(1,nb,Md,3)
!$OMP barrier

            do ib=ib1,ib2
               p=ib-ib1+1
               do m=1,Md
                  d=-ggg(4)*coef_nab(2,m)
                  do i3=a3b_omp,b3b_omp
                     do i2=a2b_omp,b2b_omp
                        do i1=a1b_omp,b1b_omp
                           j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                           htpsi(j,ib)=htpsi(j,ib) &
                                -d*(www(i1,i2-m,i3,p)-www(i1,i2+m,i3,p))
                        end do
                     end do
                  end do
               end do
            end do

!$OMP barrier

         end if

         if ( flag_n23 ) then

            do n=1,nb
               do i3=a3b_omp,b3b_omp
                  do i2=a2b_omp,b2b_omp
                     do i1=a1b_omp,b1b_omp
                        www(i1,i2,i3,n) = 0.0d0
                     end do
                  end do
               end do
            end do

            do n=1,nb
               do m=1,Md
                  d=coef_nab(2,m)
                  do i3=a3b_omp,b3b_omp
                     do i2=a2b_omp,b2b_omp
                        do i1=a1b_omp,b1b_omp
                           www(i1,i2,i3,n)=www(i1,i2,i3,n) &
                                -d*( wk(i1,i2-m,i3,n)-wk(i1,i2+m,i3,n) )
                        end do
                     end do
                  end do
               end do
            end do

!$OMP barrier
            call bcset_1(1,nb,Md,5)
!$OMP barrier

            do ib=ib1,ib2
               p=ib-ib1+1
               do m=1,Md
                  d=-ggg(5)*coef_nab(3,m)
                  do i3=a3b_omp,b3b_omp
                     do i2=a2b_omp,b2b_omp
                        do i1=a1b_omp,b1b_omp
                           j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                           htpsi(j,ib)=htpsi(j,ib) &
                                -d*(www(i1,i2,i3-m,p)-www(i1,i2,i3+m,p))
                        end do
                     end do
                  end do
               end do
            end do

!$OMP barrier

         end if

         if ( flag_n31 ) then

            do n=1,nb
               do i3=a3b_omp,b3b_omp
                  do i2=a2b_omp,b2b_omp
                     do i1=a1b_omp,b1b_omp
                        www(i1,i2,i3,n) = 0.0d0
                     end do
                  end do
               end do
            end do

            do n=1,nb
               do m=1,Md
                  d=coef_nab(3,m)
                  do i3=a3b_omp,b3b_omp
                     do i2=a2b_omp,b2b_omp
                        do i1=a1b_omp,b1b_omp
                           www(i1,i2,i3,n)=www(i1,i2,i3,n) &
                                -d*( wk(i1,i2,i3-m,n)-wk(i1,i2,i3+m,n) )
                        end do
                     end do
                  end do
               end do
            end do

!$OMP barrier
            call bcset_1(1,nb,Md,1)
!$OMP barrier

            do ib=ib1,ib2
               p=ib-ib1+1
               do m=1,Md
                  d=-ggg(6)*coef_nab(1,m)
                  do i3=a3b_omp,b3b_omp
                     do i2=a2b_omp,b2b_omp
                        do i1=a1b_omp,b1b_omp
                           j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                           htpsi(j,ib)=htpsi(j,ib) &
                                -d*(www(i1-m,i2,i3,p)-www(i1+m,i2,i3,p))
                        end do
                     end do
                  end do
               end do
            end do

         end if

      end if

      deallocate( zcoef_kin )
      deallocate( coef_lap, coef_nab )

    END SUBROUTINE op_kinetic_sol_eps

  END SUBROUTINE calc_stress_kin



  SUBROUTINE calc_stress_ion( stress_ion )
    use var_ps_member, only: Mr, parloc, Rps, rad, rab, vql
    use ps_local_module, only: simp
    use pseudopot_module, only: Zps
    use bberf_module, only: bberf
    use simc_module, only: simc

    real(8), intent(out) :: stress_ion(3,3)

    real(8) :: Rc,p1,p2,p3,p4,r,dr,sum0,G(3),G1,G2,Gr
    integer :: s,ir,ia,ie,ig,jg,ig1,ig2,ig3,la,mu,info
    real(8),    allocatable :: vshort(:), tmp(:)
    real(8),    allocatable :: rho1d(:)
    complex(8), allocatable :: rho3d(:,:,:)
    complex(8), allocatable :: dVion1d(:)
    complex(8), allocatable :: work(:,:,:)
    complex(8), allocatable :: SGK(:,:)
    real(8),    allocatable :: dvqlg(:,:)
    real(8) :: const

    allocate( dvqlg(NMGL,Nelement)  )
    dvqlg=0.0d0

    allocate( vshort(maxval(Mr)) )
    allocate( tmp(maxval(Mr)) )
    vshort=0.0d0

    do ie=1,Nelement
       Rc  = maxval( Rps(:,ie) )

       call simc(rad(1,ie),vql(1,ie),Rc,Zps(ie),parloc(1,ie),Mr(ie))
       p1 = parloc(1,ie)
       p2 = sqrt(parloc(2,ie))
       p3 = parloc(3,ie)
       p4 = sqrt(parloc(4,ie))

       do ir=1,Mr(ie)
          r = rad(ir,ie)
          vshort(ir) = vql(ir,ie) &
               - Vlong_r( Zps(ie), p1, p2, p3, p4, r )
       end do

       do ig=1,NMGL
          G1=sqrt(GG(ig))

          do ir=1,Mr(ie)
             r = rad(ir,ie)
             dr = rab(ir,ie)
             tmp(ir) = dr * r*r*r * vshort(ir)*(-1.0d0)*J1(G1*r)
          end do
          call simp(tmp,sum0,2)
          dvqlg(ig,ie) = M_4PI/Va * sum0
       end do

       do ig=1,NMGL
          G2 = GG(ig)
          dvqlg(ig,ie) = dvqlg(ig,ie) &
               + M_4PI/Va * dVlong_g(Zps(ie),p1,p2,p3,p4,G2)
       end do

    end do ! ie

    deallocate( tmp )
    deallocate( vshort )

    call construct_Ggrid(1)

    allocate( SGK(NGgrid(0),Nelement) )
    SGK=(0.d0,0.d0)

    do ia=1,Natom
       ie=ki_atom(ia)
       do ig=MG_0,MG_1
          Gr = M_2PI*sum(LLG(:,ig)*aa_atom(:,ia))
          SGK(ig,ie) = SGK(ig,ie)+dcmplx(cos(Gr),-sin(Gr))
       end do
    end do
    call destruct_Ggrid

    allocate( dVion1d(NGgrid(0)) )
    dVion1d(:) = 0.0d0
    do ie=1,Nelement
       do ig=MG_0,MG_1
          jg=MGL(ig)
          dVion1d(ig) = dVion1d(ig) + dvqlg(jg,ie)*SGK(ig,ie)
       end do
    end do
    deallocate(SGK)
    call allgatherv_Ggrid(dVion1d)

    allocate( rho1d(ML_0:ML_1) )
    rho1d(:) = 0.0d0

    do s=MSP_0,MSP_1
       rho1d(:) = rho1d(:) + rho(:,s)
    end do

    call init_fft
    call d1_to_z3_fft( rho1d, rho3d )
    call forward_fft( rho3d, work )
    if ( allocated(work) ) deallocate( work )
    call finalize_fft
    deallocate(rho1d)

    call construct_Ggrid(2)

    stress_ion(:,:) = 0.0d0
    do ig=1,NGgrid(0)

       ig1=LLG(1,ig)
       ig2=LLG(2,ig)
       ig3=LLG(3,ig)

       G(:) = 0.0d0
       if( ig1>(Ngrid(1)-1)/2 ) then
          G(:) = G(:) + bb(:,1)*(ig1 - Ngrid(1))
       else
          G(:) = G(:) + bb(:,1)*ig1
       end if

       if( ig2>(Ngrid(2)-1)/2 ) then
          G(:) = G(:) + bb(:,2)*(ig2 - Ngrid(2))
       else
          G(:) = G(:) + bb(:,2)*ig2
       end if

       if( ig3>(Ngrid(3)-1)/2 ) then
          G(:) = G(:) + bb(:,3)*(ig3 - Ngrid(3))
       else
          G(:) = G(:) + bb(:,3)*ig3
       end if

       G2=sum(G(:)**2)
       G1=sqrt(G2)
       if( G1 == 0.0d0 ) cycle

       const = real(conjg(rho3d(ig1,ig2,ig3))*dVion1d(ig))/G1
       do la=1, 3
          do mu=1, 3
             stress_ion(la,mu) = stress_ion(la,mu) + G(la)*G(mu) * const
          end do
       end do
    end do

    stress_ion(:,:) = stress_ion(:,:)  + delta(:,:)*Eion/Va

    call destruct_Ggrid

    deallocate( dvqlg )
    deallocate( dVion1d )

    if ( allocated(rho3d) ) deallocate(rho3d)

    return
  contains
    real(8) function Vlong_r( Z, p1, p2, p3, p4, r )
      real(8), intent(in) :: Z, p1, p2, p3, p4
      real(8), intent(in) :: r

      if ( r<1.d-9 ) then
         Vlong_r = -Z*M_2_SQRTPI*(p1*p2+p3*p4)
      else
         Vlong_r = -Z/r*( p1*bberf(p2*r)+p3*bberf(p4*r) )
      end if
    end function Vlong_r

    real(8) function dVlong_g( Z, p1, p2, p3, p4, G2 )
      real(8), intent(in) :: Z, p1, p2, p3, p4
      real(8), intent(in) :: G2
      real(8) :: G

      G = sqrt(G2)
      if ( G2 == 0.d0 ) then
         dVlong_g = 0.0d0
      else
         dVlong_g = 2.0d0*Z/(G*G2)*( p1*exp(-G2/(4*p2**2)) + p3*exp(-G2/(4*p4**2)) ) &
              + Z/G*( p1/(2.0d0*p2**2)*exp(-G2/(4*p2**2)) + p3/(2.0d0*p4**2)*exp(-G2/(4*p4**2)))
      end if
    end function DVlong_g


    real(8) function J1( x )
      real(8), intent(in) :: x

      if ( x<1.d-1 ) then
         J1= 10.d0/39916800.d0*x**9 - 8.d0/362880.d0*x**7 &
              + 6.d0/5040.d0*x**5 - 4.d0/120.d0*x**3 + 2.d0/6.d0*x
      else
         J1=(sin(x)/x-cos(x))/x
      end if
    end function J1

  END SUBROUTINE calc_stress_ion



  SUBROUTINE calc_stress_har( stress_har )
    real(8), intent(out) :: stress_har(3,3)

    integer :: s,ig,ig1,ig2,ig3
    real(8) :: g2,g(3), const
    real(8), allocatable :: rho1d(:)
    complex(8), allocatable :: rho3d(:,:,:)
    complex(8), allocatable :: work(:,:,:)
    integer :: la, mu

    real(8) :: Eh1

    allocate( rho1d(ML_0:ML_1) )
    rho1d(:) = 0.0d0

    do s=MSP_0,MSP_1
       rho1d(:) = rho1d(:) + rho(:,s)
    end do

    call init_fft
    call d1_to_z3_fft( rho1d, rho3d )
    call forward_fft( rho3d, work )
    if ( allocated(work) ) deallocate( work )
    call finalize_fft
    deallocate(rho1d)

    call construct_Ggrid(2)

    stress_har(:,:) = 0.0d0
    do ig=1,NGgrid(0)
       g2=GG(MGL(ig))
       if ( g2 == 0.0d0 ) cycle
       ig1=LLG(1,ig)
       ig2=LLG(2,ig)
       ig3=LLG(3,ig)
       const = abs(rho3d(ig1,ig2,ig3))**2 * M_4PI/g2**2

       G(:) = 0.0d0
       if( ig1>(Ngrid(1)-1)/2 ) then
          G(:) = G(:) + bb(:,1)*(ig1 - Ngrid(1))
       else
          G(:) = G(:) + bb(:,1)*ig1
       end if

       if( ig2>(Ngrid(2)-1)/2 ) then
          G(:) = G(:) + bb(:,2)*(ig2 - Ngrid(2))
       else
          G(:) = G(:) + bb(:,2)*ig2
       end if

       if( ig3>(Ngrid(3)-1)/2 ) then
          G(:) = G(:) + bb(:,3)*(ig3 - Ngrid(3))
       else
          G(:) = G(:) + bb(:,3)*ig3
       end if

       do la=1, 3
          do mu=1, 3
             stress_har(la,mu) = stress_har(la,mu) - G(la)*G(mu)*const
          end do
       end do
    end do
    stress_har(:,:) = stress_har(:,:) + delta(:,:) * (E_hartree/Va)

    call destruct_Ggrid
    if ( allocated(rho3d) ) deallocate(rho3d)

  end SUBROUTINE calc_stress_har


  SUBROUTINE calc_stress_xc( stress_xc )
    real(8), intent(out) :: stress_xc(3,3)

    real(8) :: Evxc
    integer :: info

    Evxc = sum(Vxc(:,:)*rho(:,:))*dV

    call mpi_allreduce(MPI_IN_PLACE,Evxc,1, &
         mpi_real8,mpi_sum,mpi_comm_world,info)

    stress_xc(:,:) = delta(:,:) &
         *(Evxc - Exc)/Va

  end SUBROUTINE calc_stress_xc


  SUBROUTINE calc_stress_nlc( stress_nlc )
    use pseudopot_module, only: NRps,norb,lo
    use ps_nloc2_variables, only: amap,lmap,mmap,iorbmap,nzlma, &
         nrlma_xyz,iuV,lma_nsend,num_2_rank,sendmap,recvmap, &
         MJJ,JJP,uVk,sbufnl,rbufnl,TYPE_MAIN,MMJJ,MJJ_MAP,JJ_MAP
    use ps_nloc2_init_module, only: rad1, dviod
    use force_sub_sub_module, only: gaunt, construct_gaunt_coef_L1
    use ylm_module, only: Ylm
    use polint_module, only: polint

    real(8), intent(out) :: stress_nlc(3,3)

    integer :: i1,i2,i3
    integer :: i,j,k,s,n,ir,iorb,L,L1,L1z,NRc,irank,jrank
    integer :: nreq,max_nreq,ib,ib1,ib2,nnn
    integer :: a,a0,ik,m,lm0,lm1,lma,im,m1,m2
    integer :: ierr,ir0,ir1
    integer,allocatable :: ireq(:),istatus(:,:),ilm1(:,:,:)
    real(8),parameter :: ep=1.d-8
    real(8),save :: Y1(0:3,-3:3,0:4,-4:4)
    real(8),save :: Y2(0:3,-3:3,0:4,-4:4)
    real(8),save :: Y3(0:3,-3:3,0:4,-4:4)
    real(8) :: err,err0
    real(8) :: Ratom(3),d(3)
    real(8) :: R(3),R1,kr,c
    real(8) :: tmp0,tmp1
    real(8) :: yy(3)
    real(8),allocatable :: duVdeps(:,:,:,:) ! (la,mu,grid,lma)
    complex(8) :: ztmp
    complex(8),allocatable :: wtmp5(:,:,:,:,:),vtmp2(:,:,:)
    logical,allocatable :: a_rank(:)
    integer :: k1,k2,k3,a1b,a2b,a3b,ab1,ab2,ab3
    integer :: la, mu
    type(gaunt) :: yyy


    call construct_gaunt_coef_L1( 3, yyy )

    a1b=Igrid(1,1)
    a2b=Igrid(1,2)
    a3b=Igrid(1,3)
    ab1=Igrid(2,1)-Igrid(1,1)+1
    ab2=Igrid(2,2)-Igrid(1,2)+1
    ab3=Igrid(2,3)-Igrid(1,3)+1

    L1=maxval(lo)+1
    n=maxval(norb)
    allocate( ilm1(0:L1,n,Nelement) ) ; ilm1=0
    do ik=1,Nelement
       lm1=0
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          do L1=abs(L-1),L+1
             lm1=lm1+1
             ilm1(L1,iorb,ik)=lm1
          end do
       end do
    end do

    allocate( wtmp5(0:9,nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1) )
    allocate( vtmp2(0:9,nzlma,MB_d) )
    allocate( a_rank(Natom) )
    allocate( duVdeps(3,3,MMJJ,nzlma) )

!$OMP parallel

!$OMP workshare
    wtmp5 = 0.0d0
    a_rank(:)=.false.
!$OMP end workshare

!$OMP do private(i1,i2,i3,k1,k2,k3)
    do a=1,Natom
       i1 = nint( aa_atom(1,a)*Ngrid(1) )
       i2 = nint( aa_atom(2,a)*Ngrid(2) )
       i3 = nint( aa_atom(3,a)*Ngrid(3) )
       k1 = i1/Ngrid(1) ; if ( i1<0 ) k1=(i1+1)/Ngrid(1)-1
       k2 = i2/Ngrid(2) ; if ( i2<0 ) k2=(i2+1)/Ngrid(2)-1
       k3 = i3/Ngrid(3) ; if ( i3<0 ) k3=(i3+1)/Ngrid(3)-1
       i1 = i1 - k1*Ngrid(1)
       i2 = i2 - k2*Ngrid(2)
       i3 = i3 - k3*Ngrid(3)
       if ( Igrid(1,1) <= i1 .and. i1 <= Igrid(2,1) .and. &
            Igrid(1,2) <= i2 .and. i2 <= Igrid(2,2) .and. &
            Igrid(1,3) <= i3 .and. i3 <= Igrid(2,3) ) then
          a_rank(a)=.true.
       end if
    end do
!$OMP end do


!$OMP workshare
    duVdeps = 0.d0
!$OMP end workshare

!$OMP do schedule(dynamic) &
!$OMP    private(a,L,m,iorb,ik,NRc,Ratom,d,R,R1, &
!$OMP            yy,lm1,ir0,ir1,ir,err0,tmp0,m1,m2,tmp1,err)
    do lma=1,nzlma
       a    = amap(lma)
       if ( a <= 0 ) cycle
       L    = lmap(lma)
       m    = mmap(lma)
       iorb = iorbmap(lma)
       ik   = ki_atom(a)
       NRc=NRps(iorb,ik)

       Ratom(:) = matmul( aa(:,:), aa_atom(:,a) )

       do j=1,MJJ_MAP(lma)
          d(1)=(1.0d0/Ngrid(1))*JJ_MAP(1,j,lma)+JJ_MAP(4,j,lma)
          d(2)=(1.0d0/Ngrid(2))*JJ_MAP(2,j,lma)+JJ_MAP(5,j,lma)
          d(3)=(1.0d0/Ngrid(3))*JJ_MAP(3,j,lma)+JJ_MAP(6,j,lma)

          R(:) = matmul( aa(:,:), d(:) ) - Ratom(:)
          R1 = sqrt(sum(R**2))

          yy(:)=0.d0
          do L1=abs(L-1),L+1
             lm1=ilm1(L1,iorb,ik)
             if ( ANY(abs(R(:))>1.d-14) .or. L1==0 ) then
                ir0 = 1
                ir1 = NRc
111             if ( ir1 - ir0 > 1 ) then
                   ir = ( ir0 + ir1 )/2
                   if ( rad1(ir,ik) > R1 ) then
                      ir1 = ir
                   else
                      ir0 = ir
                   end if
                   goto 111
                end if
                ir = ir0
                if ( ir <= 2 ) then
                   err0=0.d0
                   tmp0=dviod(2,lm1,ik)/(rad1(2,ik)**2)
                   if ( ir < 1 ) stop "calc_force_ps_nloc2"
                else if ( ir <= NRc ) then
                   err0=1.d10
                   do im=1,20
                      m1=max(1,ir-im)
                      m2=min(ir+im,NRc)
                      call polint(rad1(m1,ik),dviod(m1,lm1,ik) &
                           ,m2-m1+1,R1,tmp1,err)
                      if ( abs(err)<err0 ) then
                         tmp0=tmp1
                         err0=abs(err)
                         if ( err0<ep ) exit
                      end if
                   end do
                   tmp0=tmp0/(R1*R1)
                else
                   tmp0=0.0d0
                end if

                do L1z=-L1,L1
                   tmp1=tmp0*Ylm(R(1),R(2),R(3),L1,L1z)
                   yy(1)=yy(1)+tmp1*yyy%x(L,m,L1,L1z)
                   yy(2)=yy(2)+tmp1*yyy%y(L,m,L1,L1z)
                   yy(3)=yy(3)-tmp1*yyy%z(L,m,L1,L1z)
                end do
             end if
          end do ! L1

          do la=1, 3
             do mu=1, 3
                duVdeps(la,mu,j,lma)=R(mu)*yy(la)
             end do
          end do

       end do ! j
    end do ! lma
!$OMP end do

!$OMP single
    deallocate( ilm1 )
!$OMP end single

    do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
!$OMP do schedule(dynamic) private(c,i,i1,i2,i3,d,kr,ztmp)
          do n=MB_0,MB_1
             if ( occ(n,k,s) == 0.d0 ) cycle

             c=-2.d0*occ(n,k,s)*dV*dV

             do lma=1,nzlma
                do j=1,MJJ(lma)
                   i=JJP(j,lma)
                   wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s) &
                        +uVk(j,lma,k)*conjg(unk(i,n,k,s))
                end do
                wtmp5(0,lma,n,k,s)=iuV(lma)*c*wtmp5(0,lma,n,k,s)
                do j=1,MJJ_MAP(lma)
                   i1=JJ_MAP(1,j,lma)
                   i2=JJ_MAP(2,j,lma)
                   i3=JJ_MAP(3,j,lma)
                   i = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ML_0
                   d(1)=(1.0d0/Ngrid(1))*i1+JJ_MAP(4,j,lma)
                   d(2)=(1.0d0/Ngrid(2))*i2+JJ_MAP(5,j,lma)
                   d(3)=(1.0d0/Ngrid(3))*i3+JJ_MAP(6,j,lma)

                   kr=M_2PI*sum( kbb(:,k)*d(:) )
                   ztmp=dcmplx(cos(kr),sin(kr))*unk(i,n,k,s)

                   do la=1, 3
                      do mu=1, 3
                         wtmp5(la+(mu-1)*3,lma,n,k,s) &
                              = wtmp5(la+(mu-1)*3,lma,n,k,s) &
                              + duVdeps(la,mu,j,lma)*ztmp
                      end do
                   end do
                end do
             end do ! lma
          end do ! n
!$OMP end do
       end do ! k
    end do ! s
!$OMP end parallel

    deallocate( duVdeps )
    max_nreq=2*maxval( nrlma_xyz )
    allocate( ireq(max_nreq) )
    allocate( istatus(MPI_STATUS_SIZE,max_nreq) )

    stress_nlc(:,:)=0.d0
    do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
          do n=MB_0,MB_1,MB_d

             ib1=n
             ib2=min(ib1+MB_d-1,MB_1)
             nnn=ib2-ib1+1

             if ( occ(n,k,s) == 0.d0 ) cycle

             do i=1,6
                select case(i)
                case(1,3,5)
                   j=i+1
                   vtmp2(:,:,1:nnn)=wtmp5(:,:,ib1:ib2,k,s)
                case(2,4,6)
                   j=i-1
                end select

                do m=1,nrlma_xyz(i)
                   nreq=0
                   irank=num_2_rank(m,i)
                   jrank=num_2_rank(m,j)
                   if( irank >= 0 )then
                      i1=0
                      do ib=1,nnn
                         do i2=1,lma_nsend(irank)
                            do i3=0,9
                               i1=i1+1
                               sbufnl(i1,irank)=vtmp2(i3,sendmap(i2,irank),ib)
                            end do
                         end do
                      end do
                      nreq=nreq+1
                      call mpi_isend(sbufnl(1,irank),10*lma_nsend(irank)*nnn &
                           ,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
                   end if
                   if( jrank >= 0 )then
                      nreq=nreq+1
                      call mpi_irecv(rbufnl(1,jrank),10*lma_nsend(jrank)*nnn &
                           ,TYPE_MAIN,jrank,1,comm_grid,ireq(nreq),ierr)
                   end if
                   call mpi_waitall(nreq,ireq,istatus,ierr)
                   if( jrank >= 0 )then
                      i1=0
                      do ib=ib1,ib2
                         do i2=1,lma_nsend(jrank)
                            do i3=0,9
                               i1=i1+1
                               wtmp5(i3,recvmap(i2,jrank),ib,k,s) &
                                    =wtmp5(i3,recvmap(i2,jrank),ib,k,s)+rbufnl(i1,jrank)
                            end do
                         end do
                      end do
                   end if
                end do ! m
             end do ! i

             do ib=ib1,ib2
                do lma=1,nzlma
                   a=amap(lma)

                   if ( a <= 0 ) cycle
                   if ( a_rank(a) ) then
                      do la=1, 3
                         do mu=1, 3
                            stress_nlc(la,mu) = stress_nlc(la,mu) &
                                 + real(wtmp5(0,lma,ib,k,s)*wtmp5(la+(mu-1)*3,lma,ib,k,s))
                         end do
                      end do
                   end if
                end do
             end do

          end do ! n
       end do ! k
    end do ! s

    deallocate( istatus )
    deallocate( ireq )

    stress_nlc(:,:) = stress_nlc(:,:)*(-1.0d0/Va)

    call mpi_allreduce(MPI_IN_PLACE,stress_nlc,3*3, &
         mpi_real8,mpi_sum,mpi_comm_world,ierr)
    deallocate( a_rank )
    deallocate( vtmp2 )
    deallocate( wtmp5 )

    stress_nlc(:,:) = stress_nlc(:,:) - delta(:,:)*Enlc/Va

  end SUBROUTINE calc_stress_nlc



  SUBROUTINE calc_stress_ewa( stress_ewa )
    use pseudopot_module, only: Zps
    use ewald_variables, only: mg,mr,LG,LR,ipair,mpair
  
    real(8), intent(out) :: stress_ewa(3,3)

    real(8) :: Qtot1, Qtot2
    real(8) :: Ratom(3), t, G(3), GG
    complex(8) :: sg
    integer :: ig, ia, ip, ia1, ia2, il, info
    integer :: la,mu
    real(8) :: a12(3), L(3), L1, sum0(3,3)
    real(8) :: eta
    real(8) :: stress_ewa_g(3,3), stress_ewa_r(3,3)

    ! [note], We use smaller eta then one used for ewald sum of energy
    eta = M_PI/Va**(2.d0/3.d0)/4.0d0

    stress_ewa_g(:,:) = 0.d0
    do ig=1, mg
       if ( all( LG(1:3,ig) == 0 ) ) cycle

       SG=(0.d0,0.d0)
!$OMP parallel do private(t) reduction(+:SG)
       do ia=1, Natom
          t = M_2PI*sum(LG(:,ig)*aa_atom(:,ia))
          SG = SG + Zps(ki_atom(ia))*dcmplx(cos(t),-sin(t))
       end do
!$OMP end parallel do

       G(:) = matmul(bb(:,:),LG(:,ig))
       GG   = sum(G(:)**2)

       stress_ewa_g(:,:) = stress_ewa_g(:,:) &
            + (M_2PI/Va**2 &
            * abs(SG)**2 * exp(-GG/(4.d0*eta))/GG)*delta(:,:)

       do la=1, 3
          do mu=1, 3
             stress_ewa_g(la,mu) = stress_ewa_g(la,mu) &
                  - (M_2PI/Va**2 &
                  * abs(SG)**2 * exp(-GG/(4.d0*eta)) &
                  * 2.0d0*G(la)*G(mu)/GG &
                  * (1.0d0+GG/(4.d0*eta)) )
          end do
       end do
    end do

    call mpi_allreduce(MPI_IN_PLACE,stress_ewa_g,3*3,mpi_real8,mpi_sum, &
         mpi_comm_world,info)

    Qtot1 = sum(Zps(ki_atom(:)))
    stress_ewa_g(:,:) = stress_ewa_g(:,:) &
         - M_PI/(2.0d0*eta*Va**2)*Qtot1**2 *delta(:,:)

    stress_ewa_r(:,:) = 0.d0

    do ip=1, mpair
       ia1 = ipair(1,ip)
       ia2 = ipair(2,ip)

       sum0=0.d0
!$OMP parallel do private(L,L1) reduction(+:sum0)
       do il=1, mr
          L(:) = matmul( aa(:,:), LR(:,il) + aa_atom(:,ia1) - aa_atom(:,ia2) )
          L1 = sqrt(sum(L(:)**2))
          if ( L1 == 0.d0 ) cycle

          do la=1, 3
             do mu=1, 3
                sum0(la,mu) = sum0(la,mu) + L(la)*L(mu)/L1**2*Hprim(sqrt(eta)*L1)
             end do
          end do
       end do
!$OMP end parallel do

       if ( ia1 /= ia2 ) sum0(:,:) = sum0(:,:)*2.d0

       stress_ewa_r(:,:) = stress_ewa_r(:,:) &
            + Zps(ki_atom(ia1))*Zps(ki_atom(ia2))*sum0(:,:)
    end do
    stress_ewa_r(:,:) = stress_ewa_r(:,:)*(-0.5d0*sqrt(eta)/Va)

    call mpi_allreduce(MPI_IN_PLACE,stress_ewa_r,3*3,mpi_real8,mpi_sum, &
         mpi_comm_world,info)

    stress_ewa(:,:) = stress_ewa_g(:,:) + stress_ewa_r(:,:)

    return

  contains
    real(8) function Hprim( x )
      real(8), intent(in) :: x

      Hprim = -M_2_SQRTPI * exp(-x*x) - erfc(x)/x

    end function Hprim

  END SUBROUTINE calc_stress_ewa


  SUBROUTINE calc_stress_cnt( stress_cnt )
    real(8), intent(out) :: stress_cnt(3,3)

    stress_cnt(:,:) = delta(:,:) * (cnst/Va)

  end SUBROUTINE calc_stress_cnt

END MODULE stress_module
