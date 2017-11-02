MODULE atomopt_module

  use parallel_module
  use atom_module
  use total_energy_module
  use bb_module
  use scf_module
  use eion_module, only: calc_eion
  use strfac_module
!  use ps_local_module
!  use ps_pcc_module
!  use pseudopot_module
  use ps_nloc2_module
  use ps_nloc3_module
  use ps_nloc_mr_module
  use force_module
  use kinetic_module, only: SYStype
  use ps_local_mol_module, only: construct_ps_local_mol
  use ps_nloc2_mol_module
  use ps_pcc_mol_module
  use eion_mol_module
  use ps_qrij_prep_module
  use ps_prepNzqr_g_module, only: prepNzqr
  use vdw_grimme_module
  use efield_module
  !--- begin MIZUHO-IR for cellopt
  use aa_module, only: ax, aa, Va
  use stress_module, only: calc_stress
  use rgrid_module, only: Ngrid,Hgrid,Igrid,dV,Init_Rgrid,InitParallel_Rgrid
  use ggrid_module, only: Init_Ggrid,InitParallel_Ggrid,Gcut
  use kinetic_variables, only: Md, ggg
  !use bz_module, only: kbb,Nbzsm,generate_bz
  use bz_module, only: kbb
  use kinetic_module, only: init_kinetic
  use ps_local_module, only: init_ps_local, Vion, construct_ps_local
  use ps_pcc_module, only: init_ps_pcc, construct_ps_pcc
  use ps_nloc_initiate_module, only: ps_nloc_initiate
  use pseudopot_module, only: read_pseudopot, pselect
  use lattice_module, only: lattice, backup_aa_lattice
  use aa_module, only: init_aa
  !--- end MIZUHO-IR for cellopt
  use atomopt_rf_module
  use atomopt_bfgs_module
  use atomopt_diis_module

  implicit none

  PRIVATE
  PUBLIC :: ncycl,most,nrfr,okatom,eeps,feps,decr
  PUBLIC :: read_atomopt
  PUBLIC :: atomopt

  integer :: ncycl,most,nrfr
  real(8) :: okatom,eeps,feps,decr

  logical :: disp_switch_loc
  integer :: diter_opt

  integer :: strlog = 0
  integer,parameter :: unit_atmopt = 87
  integer,parameter :: unit_strlog = 86
  integer,parameter :: unit197 = 197
  integer,parameter :: unit97 = 97
  real(8), parameter :: M_PI       = 3.14159265358979323846d0
  real(8), parameter :: M_2PI      = M_PI*2.0d0

CONTAINS


  SUBROUTINE atomopt( iswitch_opt, iswitch_latopt )
    implicit none
    integer,intent(IN) :: iswitch_opt
    integer,intent(IN) :: iswitch_latopt
    select case( iswitch_opt )
    case( 1, 2 )
       call atomopt_cg( iswitch_opt, iswitch_latopt )
    case( 4 )
       call atomopt_diis( SYStype, feps, diter_opt )
    case( 5 )
       call atomopt_bfgs( SYStype, feps, diter_opt )
    case( 6 )
       call atomopt_rf( SYStype, feps, diter_opt )
    case default
       write(*,*) "Invalid Parameter: iswitch_opt=",iswitch_opt
       call stop_program("atomopt")
    end select
  END SUBROUTINE atomopt


  SUBROUTINE read_atomopt(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(8) :: cbuf,ckey
    ncycl     = 0
    most      = 6
    nrfr      = 5
    diter_opt = 50
    okatom    = 0.5d0
    eeps      = 1.d-10
    feps      = 5.d-4
    decr      = 1.d-1

    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:8) == "ATOMOPT1" ) then
             backspace(unit)
             read(unit,*) cbuf,ncycl,most,nrfr
          else if ( ckey(1:8) == "ATOMOPT2" ) then
             backspace(unit)
             read(unit,*) cbuf,okatom,eeps,feps,decr
          else if ( ckey(1:8) == "ATOMOPT3" ) then
             backspace(unit)
             read(unit,*) cbuf,diter_opt
          else if ( ckey(1:8) == "STRLOG" ) then
             backspace(unit)
             read(unit,*) cbuf,strlog
          end if
       end do
999    continue
       write(*,*) "ncycl, most, nrfr =",ncycl,most,nrfr
       write(*,*) "okatom, eeps      =",okatom,eeps
       write(*,*) "feps, decr        =",feps,decr
       write(*,*) "diter_opt         =",diter_opt
       if ( diter_opt <= 0 ) then
          diter_opt=50
          write(*,*) "diter_opt         =",diter_opt
       end if
       write(*,*) "strlog            =",strlog
    end if
    call send_atomopt
  END SUBROUTINE read_atomopt


  SUBROUTINE send_atomopt
    implicit none
    integer :: ierr
    call mpi_bcast(ncycl,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(most ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nrfr ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(okatom,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(eeps  ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(feps  ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(decr  ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(diter_opt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(strlog,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_atomopt


!-------------------------------------------------- atomopt_cg

  SUBROUTINE atomopt_cg( iswitch_opt, iswitch_latopt )
    implicit none
    integer,intent(IN) :: iswitch_opt
    integer,intent(IN) :: iswitch_latopt ! MIZUHO-IR for cellopt
    integer,parameter :: max_nhist=1000
    integer :: SCF_hist(max_nhist),ICY_hist(max_nhist)
    integer :: LIN_hist(max_nhist)
    integer :: a,nhist,most0,ncycl0,ncycl1,ifar,ierr,icy,nhist0
    integer :: i,itlin,amax,isafe,icflag,iter_final
    real(8) :: Fmax_hist(max_nhist),Etot_hist(max_nhist)
    real(8) :: dmax_hist(max_nhist)
    real(8) :: grad_hist(max_nhist),alpha_hist(max_nhist)
    real(8) :: Fmax,ss,Etot_save,dif,alpha1_0,okstep0
    real(8) :: al(3),all(3),ar(3),arr(3),alp,wd,xmin,emin
    real(8) :: Etot,Etot0,Fmax0_2,alpha0,Fmax0,Fmax00,Etot00,Etsave
    real(8) :: gh,alpha,tmp,ddmax,almax,okstep,alpha0_2,Etot0_2
    real(8) :: alpha2,alpha1,grad0,signh,dif0,gamma,gigi,hh
    real(8) :: alsave,safe,c0,c1,c2,c3,pi,ddmin,grad,safety
    real(8),allocatable :: Force(:,:),aa_atom_0(:,:)
    real(8),allocatable :: gi(:,:),hi(:,:)
    character(22) :: loop_info
    !--- begin MIZUHO-IR for cellopt
    real(8) :: stress(3,3)
    integer :: dim_opt
    type(lattice) :: aa_obj
    !--- end MIZUHO-IR for cellopt

    call write_border( 0, " atomopt_cg(start)" )

    call check_disp_switch( disp_switch_loc, 0 )

    ddmin  = 1.d-8
    safe   = 0.01d0
    safety = 0.01d0
    pi     = acos(-1.d0)

    ! MIZUHO-IR for cellopt
    if( iswitch_latopt >= 1 ) then
       dim_opt = Natom + 3
    else
       dim_opt = Natom
    end if

!- allocate ---------------------------------------
    allocate( Force(3,dim_opt)    ) ; Force=0.d0
    allocate( aa_atom_0(3,dim_opt) ) ; aa_atom_0=0.d0
    allocate( gi(3,dim_opt)       ) ; gi=0.d0
    allocate( hi(3,dim_opt)       ) ; hi=0.d0
!--------------------------------------------------

    if ( iswitch_opt < 2  ) then

       if ( iswitch_opt == 1 ) then
          call calc_total_energy( .false., Etot )
          if ( disp_switch_loc ) write(*,*) "Etot(har)=",Etot
       end if
       if ( disp_switch_loc ) write(*,'(1x,"# Force (total)")')
       call calc_force( Natom, Force )
       ! MIZUHO-IR for cellopt
       if( .not. iswitch_opt >= 1 ) then
          Force(:,1:Natom) = 0.0d0
       end if
       ! MIZUHO-IR for cellopt
       if( iswitch_latopt >= 1 ) then
          call calc_total_energy( .false., Etot )
          call calc_stress( stress )

          Force(:,Natom+1) = Va/M_2PI*matmul( stress(:,:), bb(:,1) )
          Force(:,Natom+2) = Va/M_2PI*matmul( stress(:,:), bb(:,2) )
          Force(:,Natom+3) = Va/M_2PI*matmul( stress(:,:), bb(:,3) )
       end if

       Fmax=0.d0
       do a=1,dim_opt
          ss = sum(Force(:,a)**2)
          Fmax=max(Fmax,ss)
       end do
       Fmax=sqrt(Fmax)

       Etot_save = 0.d0
       dif       = 0.0d0

       nhist                = 1
       Etot_hist(nhist)     = Etot
       alpha_hist(nhist)    = 0.0d0
       Fmax_hist(nhist)     = Fmax
       SCF_hist(nhist)      = 0
       ICY_hist(nhist)      = 0
       LIN_hist(nhist)      = 0
       dmax_hist(nhist)     = 0

       ncycl0   = 1
       most0    = 1

! --- Convergence check 1 ---

       if ( Fmax <= feps ) then
          if ( disp_switch_loc ) then
             write(*,*) 'Fmax,feps =',Fmax,feps
          end if
          goto 999
       end if

! --- Read the previous optimization information ---

    else if ( iswitch_opt >= 2 ) then

       if ( myrank == 0 ) then
          open(1,file="wopt.dat",form='unformatted',status='old')
          read(1) ncycl0,most0,ifar
          read(1) Etot,dif,alpha1_0,okstep0
          read(1) al(1:3),all(1:3),ar(1:3),arr(1:3)
          read(1) Force(:,:)
          read(1) hi(:,:)
          read(1) aa_atom(:,:)
          close(1)
       end if
       call mpi_bcast(ncycl0,1,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(most0,1,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(ifar,1,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(Etot,1,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(dif,1,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(alpha1_0,1,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(okstep0,1,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(al,3,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(all,3,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(ar,3,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(arr,3,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(Force,3*dim_opt,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(hi,3*dim_opt,mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(aa_atom,3*Natom,mpi_real8,0,mpi_comm_world,ierr)

       Fmax=0.d0
       do a=1,dim_opt
          ss = sum(Force(:,a)**2)
          Fmax=max(Fmax,ss)
       end do
       Fmax=sqrt(Fmax)

       nhist = 0

    end if

    call write_atomic_coordinates_log(unit197,0,0,strlog,iswitch_opt)
    if ( strlog == 3 .and. iswitch_opt < 2 ) call write_etot_force_log( 0,0, Etot, Force )

    dif0 = dif

    if ( disp_switch_loc ) then
       write(*,*) "ncycl0,ncycl ",ncycl0,ncycl
       write(*,*) "most  ",most
       write(*,*) "okatom",okatom
    end if

!
! -------------------- CG-loop start --------------------
!

    ncycl1 = ncycl0 + ncycl - 1

    opt_ion : do icy=ncycl0,ncycl1+1

       write(loop_info,'(" ICY    (",i5,")")') icy
       call write_border( 0, loop_info(1:len_trim(loop_info)) )

!
! --- Ion configuration ---
!

       aa_atom_0(1:3,1:Natom) = aa_atom(1:3,1:Natom)
       ! MIZUHO-IR for cellopt
       if( iswitch_latopt >= 1 ) then
          aa_atom_0(1:3,Natom+1:Natom+3) = aa(1:3,1:3)
       end if
!
! gi ---> gradient ( =-grad(Etot)=Force )
! hi ---> search direction
!
       if ( mod(icy-1,nrfr) == 0 ) then

          gamma=0.0d0
          if ( icy > 1 .and. disp_switch_loc ) then
             write(*,*) 'CG-direction is refreshed !!!'
          else
             if ( disp_switch_loc ) write(*,*) 'The first CG step !'
          end if

       else

          gigi=sum(gi(:,:)*gi(:,:))
          if ( gigi > 0.0d0 ) then
             gamma=sum((Force(:,:)-gi(:,:))*Force(:,:))/gigi
          end if

       end if

       gi(:,:)=Force(:,:)

       if ( iswitch_opt >= 2 .and. icy == ncycl0 ) then
       else
          hi(:,:)=gi(:,:)+gamma*hi(:,:)
       end if

       Etot_save = Etot
       hh        = sqrt( sum(hi(:,:)*hi(:,:)) )
       gh        = sum( gi(:,:)*hi(:,:) )
       alpha     = 2.0d0*abs(dif/gh)
       if ( dif == 0.0d0 ) alpha = 0.5d0

!
! --- Check alpha ---
!

       okstep = okatom
!chstep
       tmp    = 0.d0
       amax   = 1
       do a=1,dim_opt
          ss = sum(hi(:,a)*hi(:,a))
          if ( ss > tmp ) then
             tmp =ss
             amax=a
          end if
       end do
       ddmax=sqrt(tmp)*abs(alpha)
       if ( ddmax < 1.d-14 ) then
          if ( disp_switch_loc ) then
             write(*,*) "ddmax is too small : ddmax =",ddmax
             write(*,*) "alpha,amax =",alpha,amax
          end if
          ddmax=1.d-12
       end if
       almax=okstep/ddmax*abs(alpha)
!chstep
       if ( disp_switch_loc ) then
          write(*,'(1x,"Maximum displacement size =",f16.7,3x,"( atom =",i5," )")') ddmax,amax
          write(*,'(1x,"okstep, alpha, almax =",3g16.6)') okstep,alpha,almax
       end if
       if ( ddmax > okstep ) then
          alpha=almax
          if ( disp_switch_loc ) then
             write(*,*) "alpha is too large and replaced by almax."
          end if
       end if

!
! --- hi-direction component of gi, and its sign ---
!

       grad0 = gh
       if ( disp_switch_loc ) write(*,*) "grad0 =",grad0
       if ( grad0 >= 0.d0 ) then
          signh=1.d0
       else
          if ( disp_switch_loc ) write(*,*) 'sign of hi is changed.'
          signh=-1.d0
          hi(:,:)=-hi(:,:)
          grad0=-grad0
       end if

       if ( iswitch_opt>=2 .and. icy==ncycl0 ) then

          alpha1   = alpha1_0
          alpha1_0 = 0.d0
          alpha2   = 0.d0
          okstep   = okstep0

       else

          all(1:3) =-1.d15
          ar(1:3)  = 1.d15
          arr(1:3) = 1.d15

          alpha1 = 0.d0
          alpha2 = 0.d0

          ifar = 0

          al(1) = 0.d0
          al(2) = Etot
          al(3) = grad0

       end if

       Etot00 = Etot
       Fmax00 = Fmax

       if ( disp_switch_loc ) then
          write(*,*) "initial configuration"
          if ( Natom <= 11 ) then
             do a=1,Natom
                write(*,'(1x,i4,3f15.5)') a,aa_atom(:,a)
             end do
          else
             do a=1,min(5,Natom)
                write(*,'(1x,i4,3f15.5)') a,aa_atom(:,a)
             end do
             write(*,'(1x,10x,".")')
             write(*,'(1x,10x,".")')
             write(*,'(1x,10x,".")')
             do a=Natom-5,Natom
                write(*,'(1x,i4,3f15.5)') a,aa_atom(:,a)
             end do
          end if
          ! MIZUHO-IR for cellopt
          if( iswitch_latopt >= 1 ) then
             write(*,'(1x,a4,3f15.5)') "A",aa(:,1)
             write(*,'(1x,a4,3f15.5)') "B",aa(:,2)
             write(*,'(1x,a4,3f15.5)') "C",aa(:,3)
          end if
          write(*,*) "Etot =",Etot
          write(*,*) "grad =",grad0
          write(*,*) "Fmax =",Fmax
       end if

       nhist = nhist + 1

       Etot_hist(nhist)     = Etot
       alpha_hist(nhist)    = 0.d0
       Fmax_hist(nhist)     = Fmax
       grad_hist(nhist)     = grad0
       SCF_hist(nhist)      = 0
       ICY_hist(nhist)      = icy
       LIN_hist(nhist)      = 0
       dmax_hist(nhist)     = sqrt(maxval(hi(1,:)**2+hi(2,:)**2+hi(3,:)**2))

       nhist0 = nhist

       Fmax0  = Fmax
       alpha0 = 0.d0
       Etot0  = Etot

       Fmax0_2  = Fmax
       alpha0_2 = 0.d0
       Etot0_2  = Etot

!
! ---------- Line minimization start (along hi) ----------
!

       linmin : do itlin=most0,most

          write(loop_info,'(" ICY    (",i5,")")') icy
          call write_border( 0, loop_info(1:len_trim(loop_info)) )
          write(loop_info,'(" LINMIN (",i5,")")') itlin
          call write_border( 0, loop_info(1:len_trim(loop_info)) )

          Etsave = Etot
          emin   = 0.d0
          ddmax  = 1000.d0

          if ( disp_switch_loc ) write(*,*) "ifar     =",ifar

          alpha1_0 = alpha1
          okstep0  = okstep

          if ( itlin==1 ) then

             alpha1=alpha

          else

             if ( ifar==0 ) then

                okstep=okstep*2.d0
                alpha2=alpha1
                call get_min_parabola(all,al,xmin,emin)
                alpha1=xmin
                alp=alpha1-alpha2
                if ( disp_switch_loc ) then
                   write(*,*) "alpha1,alpha2=",alpha1,alpha2
                end if
!chstep
                tmp=0.d0
                amax=1
                do a=1,dim_opt
                   ss = sum(hi(:,a)*hi(:,a))
                   if ( ss > tmp ) then
                      tmp=ss
                      amax=a
                   end if
                end do
                ddmax=sqrt(tmp)*abs(alp)
                if ( ddmax < 1.d-14 ) then
                   if ( disp_switch_loc ) then
                      write(*,*) "ddmax is too small : ddmax =",ddmax
                      write(*,*) "alp1,amax =",alp,amax
                   end if
                   ddmax=1.d-12
                end if
                almax=okstep/ddmax*abs(alp)
!chstep
                if ( disp_switch_loc ) then
                   write(*,'(1x,"Maximum displacement =",f16.7,3x,"( atom =",i5," )")') ddmax,amax
                   write(*,'(1x,"okstep, alp, almax =",3g16.6)') &
                        okstep,alp,almax
                end if
                if ( ddmax > okstep ) then
                   if ( disp_switch_loc ) then
                      write(*,*) "alpha1-alpha2 is ","too large and replaced by almax."
                   end if
                   alpha1=alpha2+almax
                   ddmax=okstep
                else if ( alpha1 < alpha2 ) then
                   if ( disp_switch_loc ) then
                      write(*,*) "alpha1 is irregular"
                   end if
                   tmp=abs( al(1)-all(1) )*10.d0
                   alpha1=alpha2+tmp
                   ddmax=tmp/almax*okstep
                end if

             else !( ifar/=0 )

                alpha2=alpha1
                call get_min_parabola(al,ar,xmin,emin)
                alpha1=xmin
                if ( .not.( (al(1)<alpha1).and.(alpha1<ar(1)) ) ) then
                   alpha1=0.5d0*(al(1)+ar(1))
                end if
                alp=alpha1-alpha2
!chstep
                tmp=0.d0
                amax=1
                do a=1,dim_opt
                   ss = sum(hi(:,a)*hi(:,a))
                   if ( ss>tmp ) then
                      tmp=ss
                      amax=a
                   end if
                end do
                ddmax=sqrt(tmp)*abs(alp)
                if ( ddmax < 1.d-14 ) then
                   if ( disp_switch_loc ) then
                      write(*,*) "ddmax is too small : ddmax =",ddmax
                      write(*,*) "alp1,amax =",alp,amax
                   end if
                   ddmax=1.d-12
                end if
                almax=okstep/ddmax*abs(alp)
!chstep
                if ( disp_switch_loc ) then
                   write(*,'(1x,"Maximum displacement =",f16.7,3x,"( atom =",i5," )")') ddmax,amax
                   write(*,'(1x,"okstep, alp, almax =",3g16.6)') &
                        okstep,alp,almax
                end if

!
! --- safety mechanism for stable convergence ---
!

                wd=abs( ar(1)-al(1) )
                isafe=0
                alsave=alpha1
                if ( abs(alpha1-al(1)) < safety*wd ) then
                   isafe=1
                   alpha1=al(1)+safe*(ar(1)-al(1))
                else if ( abs(ar(1)-alpha1) < safety*wd ) then
                   isafe=1
                   alpha1=ar(1)-safe*(ar(1)-al(1))
                end if
                if ( isafe == 1 ) then
                   ddmax=ddmax*abs( (alpha1-alpha2)/alp )
                   if ( disp_switch_loc ) then
                      write(*,*) "Safety mechanism is applied."
                      write(*,*) "safety & safe =",safety,safe
                      write(*,*) "alpha1 is replaced by ",alpha1
                   end if
                end if

             end if

          end if

!
! ---- save the optimization information ---
!

          if ( myrank == 0 ) then
             open(1,file="wopt.dat",form='unformatted')
             write(1) icy,itlin,ifar
             write(1) Etot_save,dif0,alpha1_0,okstep0
             write(1) al,all,ar,arr
             write(1) gi(:,:)
             write(1) hi(:,:)
             write(1) aa_atom_0(:,:)
             close(1)
          end if

          if ( icy == ncycl1 + 1 ) exit opt_ion

!
! --- Trial Configuration ---
!                  
          if ( SYStype == 0 ) then

             c0=alpha1/(2.d0*pi)
             do a=1,Natom
                aa_atom(:,a) = aa_atom_0(:,a) &
                     + alpha1/M_2PI*matmul( hi(:,a),bb(:,:) )
             end do
             ! MIZUHO-IR for cellopt
             if( iswitch_latopt >= 1 ) then
                aa(:,1) = aa_atom_0(:,Natom+1) + alpha1*hi(:,Natom+1)
                aa(:,2) = aa_atom_0(:,Natom+2) + alpha1*hi(:,Natom+2)
                aa(:,3) = aa_atom_0(:,Natom+3) + alpha1*hi(:,Natom+3)
             end if
             
          else if ( SYStype == 1 ) then

             do a=1,Natom
                aa_atom(1,a) = aa_atom_0(1,a) + alpha1*hi(1,a)
                aa_atom(2,a) = aa_atom_0(2,a) + alpha1*hi(2,a)
                aa_atom(3,a) = aa_atom_0(3,a) + alpha1*hi(3,a)
             end do

          end if

          if ( disp_switch_loc ) then
             write(*,*) 'Trial configuration (see fort.97)'
             if ( Natom <= 11 ) then
                do a=1,Natom
                   write(* ,'(1x,i5,3f20.12,i4)') &
                        ki_atom(a),aa_atom(:,a),md_atom(a)
                end do
             else
                do a=1,min(5,Natom)
                   write(* ,'(1x,i5,3f20.12,i4)') &
                        ki_atom(a),aa_atom(:,a),md_atom(a)
                end do
                write(*,'(1x,10x,".")')
                write(*,'(1x,10x,".")')
                write(*,'(1x,10x,".")')
                do a=Natom-5,Natom
                   write(* ,'(1x,i5,3f20.12,i4)') &
                        ki_atom(a),aa_atom(:,a),md_atom(a)
                end do
             end if
             ! MIZUHO-IR for cellopt
             if( iswitch_latopt >= 1 ) then
                write(*,'(1x,a4,3f15.5)') "A",aa(:,1)
                write(*,'(1x,a4,3f15.5)') "B",aa(:,2)
                write(*,'(1x,a4,3f15.5)') "C",aa(:,3)
             end if
          end if

          call write_atomic_coordinates_log(97,icy,itlin,0,iswitch_opt)

!
! --- SCF ---
!

          select case(SYStype)
          case default
             !--- begin MIZUHO-IR for cellopt
             if( iswitch_latopt >= 1 ) then
                ! update Va, volume of the cell
                ax = 1.0d0
                call init_aa( aa_obj )
                ! update aa_backup trough aa_obj
                call backup_aa_lattice( aa_obj )

                ! update Hgrid, grid spacing Hgrid.
                ! update dV, volume of a grid element.
                call Init_Rgrid( SYStype, Md, 2 )

                ! update bb, reciprocal vectors.
                call construct_bb(aa)
                ! update MGL and GG, reciprocal grid.
                call Init_Ggrid( Ngrid, bb, Hgrid, disp_switch_loc )

                ! update kbb, sampling K points
                !call generate_bz

                ! update MG_0 and MG_1, index of G grid.
                call InitParallel_Ggrid( nprocs, myrank )
                ! update ggg, zcoef_kin, etc, FDM coefficients.
                call init_kinetic( aa, bb, size(kbb,2), kbb, Hgrid, Igrid, MB_d, disp_switch_loc )

                ! update vqlg using updated Va.
                call init_ps_local
                ! update cdcg using updated Va.
                call init_ps_pcc
             end if
             !--- end begin MIZUHO-IR for cellopt.

             call calc_eion

             call construct_strfac
             call construct_ps_local
             call construct_ps_pcc
             call destruct_strfac

             !--- begin MIZUHO-IR for cellopt
             if( iswitch_latopt >= 1 ) then
                ! re-read viod and re-construct them using updated Va
                call read_pseudopot( Nelement, myrank )
                call ps_nloc_initiate( Gcut )
             end if
             !--- end begin MIZUHO-IR for cellopt

             select case(pselect)
             case(2)
                call prep_ps_nloc2
             case(3)
                call prep_ps_nloc3
             case(5)
                call prep_ps_nloc_mr
             case(102)
                call prep_ps_nloc2
                call prepNzqr
                call prepQRijp102
             end select

          case(1)

             call calc_eion

             call construct_ps_local_mol
             call construct_ps_pcc_mol
             call prep_ps_nloc2_mol

          end select

          call sawtooth_efield( Vion )

          call calc_E_vdw_grimme( aa_atom )

          write(loop_info,'("( linmin:",i3,", cg:",i3," )")') itlin,icy
          call calc_scf( disp_switch_loc,ierr,diter_opt,feps,loop_info,Etot )

          if ( ierr == -1 .or. ierr == -3 ) then
             exit opt_ion
          else if ( ierr == -2 ) then
             if ( myrank == 0 ) write(*,*) "SCF is not converged"
          end if
          iter_final=ierr

          call calc_force( Natom, Force )
          ! MIZUHO-IR for cellopt
          if( .not. iswitch_opt >= 1 ) then
             Force(:,1:Natom) = 0.0d0
          end if
          ! MIZUHO-IR for cellopt
          if( iswitch_latopt >= 1 ) then
             call calc_stress( stress )
             Force(:,Natom+1) = Va/M_2PI*matmul( stress(:,:), bb(:,1) )
             Force(:,Natom+2) = Va/M_2PI*matmul( stress(:,:), bb(:,2) )
             Force(:,Natom+3) = Va/M_2PI*matmul( stress(:,:), bb(:,3) )
          end if

          if ( disp_switch_loc ) then
             write(*,'(1x,"# Force (total)")')
             if ( Natom <= 11 ) then
                do a=1,Natom
                   write(*,'(1x,i4,i3,3g21.12)') a,ki_atom(a),force(:,a)
                end do
             else
                do a=1,min(5,Natom)
                   write(*,'(1x,i4,i3,3g21.12)') a,ki_atom(a),force(:,a)
                end do
                write(*,'(1x,10x,".")')
                write(*,'(1x,10x,".")')
                write(*,'(1x,10x,".")')
                do a=Natom-5,Natom
                   write(*,'(1x,i4,i3,3g21.12)') a,ki_atom(a),force(:,a)
                end do
             end if
             ! MIZUHO-IR for cellopt
             if( iswitch_latopt >= 1 ) then
                write(*,'(1x,a4,3f15.5)') "A",force(:,Natom+1)
                write(*,'(1x,a4,3f15.5)') "B",force(:,Natom+2)
                write(*,'(1x,a4,3f15.5)') "C",force(:,Natom+3)
             end if
          end if

          grad=sum( Force(:,:)*hi(:,:) )

          Fmax=0.d0
          do a=1,dim_opt
             ss = sum(Force(:,a)**2)
             Fmax=max(Fmax,ss)
          end do
          Fmax=sqrt( Fmax )

          if ( (al(2)<Etot) .or. (grad<=0.d0) ) then
             ifar  = 1
             arr(1:3) = ar(1:3)
             ar(1)    = alpha1
             ar(2)    = Etot
             ar(3)    = grad
          else
             all(1:3) = al(1:3)
             al(1)    = alpha1
             al(2)    = Etot
             al(3)    = grad
          end if

!
! --- history ---
!

          nhist = nhist + 1

          Etot_hist(nhist)     = Etot
          alpha_hist(nhist)    = alpha1
          Fmax_hist(nhist)     = Fmax
          grad_hist(nhist)     = grad
          SCF_hist(nhist)      = iter_final
          ICY_hist(nhist)      = icy
          LIN_hist(nhist)      = itlin
          dmax_hist(nhist)     = &
               abs(alpha1)*sqrt(maxval(hi(1,:)**2+hi(2,:)**2+hi(3,:)**2))

          if ( disp_switch_loc ) then
             write(*,'(1x,3x,1x,a12,1x,a20,1x,3a13)') &
                  'alpha    ','Etot    ','grad    ','Fmax    ','dmax    '
             do i=nhist0,nhist
                write(*,'(1x,i3,1x,g13.6,1x,g20.10,1x,3g13.5,i5)') &
                     i-nhist0,alpha_hist(i),Etot_hist(i) &
                     ,grad_hist(i),Fmax_hist(i),dmax_hist(i),SCF_hist(i)
             end do
          end if

          if ( strlog >= 2 ) then
             call write_atomic_coordinates_log(unit97,icy,itlin,strlog,iswitch_opt)
             if ( strlog == 3 ) call write_etot_force_log(icy,itlin,Etot,Force)
          end if
!
! --- Convergence check 2 ---
!
          dif = Etot - Etsave

          if ( abs(grad) <= abs(grad0*decr) ) then
             icflag=0
             if ( disp_switch_loc ) then
                write(*,*) "--- Convergence achieved in LINMIN ---"
                write(*,*) "decr =",decr
                write(*,*) 'grad, grad0*decr =',grad,grad0*decr
             end if
          else if ( Fmax <= feps ) then
             icflag=0
             if ( disp_switch_loc ) then
                write(*,*) "--- Convergence achieved in LINMIN ---"
                write(*,*) "Fmax =",Fmax
             end if
          else if ( abs(dif) <= eeps ) then
             icflag=0
             if ( disp_switch_loc ) then
                write(*,*) "--- Convergence achieved in LINMIN ---"
                write(*,*) 'dif(Etot) =',Etot
             end if
          else if ( ddmax <= ddmin ) then
             icflag=1
             if ( disp_switch_loc ) then
                write(*,*) "--- ddmax is too small ---"
             end if
          else
             icflag=-1
          end if

          if ( Fmax < Fmax0 ) then
             Fmax0  = Fmax
             alpha0 = alpha1
             Etot0  = Etot
          end if

          if ( Etot < Etot0_2 ) then
             Fmax0_2  = Fmax
             alpha0_2 = alpha1
             Etot0_2  = Etot
          end if

          if ( icflag >= 0 ) then
             if ( Etot > Etot00 ) then
                if ( disp_switch_loc ) then
! why this kind of thing happens?
                   write(*,*) "can not exit linmin !!!"
                   write(*,*) "Etot, Etot00 =",Etot,Etot00
                end if
             else
                exit linmin
             end if
          end if

       end do linmin

       hi(:,:) = signh*hi(:,:)

       most0=1

!  Best structure on a line.

       if ( disp_switch_loc ) then
          write(*,*) 'Best structure on a line (see fort.197)'
       end if  
       call write_atomic_coordinates_log(unit197,icy,itlin,1,iswitch_opt)

!
! --- Convergence check 3 ---
!
       dif  = Etot-Etot_save
       dif0 = dif

       if ( Fmax <= feps ) then

          if ( disp_switch_loc ) then
             write(*,*) "Fmax,feps =",Fmax,feps
          end if
          exit opt_ion

       else if ( abs(dif) < eeps ) then

          if ( disp_switch_loc ) then
             write(*,*) "Etot,Etot_save =",Etot,Etot_save
          end if
          exit opt_ion

       end if

    end do opt_ion

999 continue

    if ( disp_switch_loc ) then
       write(*,*) "histry (Etot,Fmax,SCF)"
       do i=1,nhist
          write(*,'(1x,2i3,f16.8,f14.7,i6,f14.7)') ICY_hist(i) &
               ,LIN_hist(i),Etot_hist(i),Fmax_hist(i),SCF_hist(i),dmax_hist(i)
       end do
    end if


    deallocate( Force )
    deallocate( aa_atom_0, gi, hi )

    call write_border( 1, " atomopt_cg(end)" )

    return

  END SUBROUTINE atomopt_cg


  SUBROUTINE get_min_parabola(g1,g2,xmin,emin)
    implicit none
    real(8),intent(IN)  :: g1(3),g2(3)
    real(8),intent(OUT) :: xmin,emin
    real(8) :: f1(3),f2(3)
    real(8) :: a,b,c,dl,d1,d2,x

    f1(1:3) = (/ g1(1), g1(2), -g1(3) /)
    f2(1:3) = (/ g2(1), g2(2), -g2(3) /)

    a=(f2(3)-f1(3))/(f2(1)-f1(1))
    dl=abs(f2(1)-f1(1))
    if ( abs(a) <= 0.0d0 ) then
       b=(f1(3)+f2(3))*0.5d0
       if ( b >= 0.0d0 ) then
          xmin=f2(1)-1.d6*dl
       else
          xmin=f1(1)+1.d6*dl
       end if
       emin=0.0d0
    else
       b=f1(3)-a*f1(1)
       xmin=-b/a
       d1=abs(xmin-f1(1))
       d2=abs(xmin-f2(1))
       if ( d1 < d2 ) then
          x=f1(1)
          c=f1(2)-0.5d0*a*x*x-b*x
       else
          x=f2(1)
          c=f2(2)-0.5d0*a*x*x-b*x
       end if
       x=xmin
       emin=0.5d0*a*x*x+b*x+c
    end if
    if ( disp_switch_loc ) then
       write(*,*) 'a   =',a
       write(*,*) 'xmin=',xmin
       write(*,*) 'emin=',emin
    end if
    return
  END SUBROUTINE get_min_parabola
  

  SUBROUTINE write_atomic_coordinates_log(unit,icy,itlin,flag,iswitch_opt)
    implicit none
    integer, intent(IN) :: unit,icy,itlin,flag,iswitch_opt

    character(8) :: fname

    if ( myrank == 0 ) then
       write(fname,'(i3)') unit
       fname="fort."//adjustl(fname)
       call write_coordinates_atom( unit_atmopt, 3, fname )
    end if

    if ( strlog /= 0 .and. flag == strlog ) then
       if ( icy == 0 ) then
          if ( myrank == 0 ) open(unit_strlog,file="strlog.dat")
          if ( iswitch_opt >= 2 ) return
          if ( myrank == 0 ) call write_coordinates_atom( unit_strlog, 1 )
       end if
       if ( myrank == 0 ) then
          write(unit_strlog,'("#_STRLOG_",a63," icy, itlin =",2(X,I3))') &
               repeat("-",63),icy,itlin
          call write_coordinates_atom( unit_strlog, 2 )
       end if
    end if

  END SUBROUTINE write_atomic_coordinates_log  


  SUBROUTINE write_etot_force_log( icy, itlin, Etot, Force )
    implicit none
    integer,intent(IN) :: icy, itlin
    real(8),intent(IN) :: etot, force(:,:)
    integer :: i
    if ( myrank == 0 ) then
       write(unit_strlog,'("#_Etot_and_Force_",a63," icy, itlin =",2(X,I3))') &
            repeat("-",63),icy,itlin
       write(unit_strlog,'("Etot(hartree)        :",f24.16)') Etot
       write(unit_strlog,'("Force(Cartesian,hartree/bohr):")')
       do i=1,size(Force,2)
          write(unit_strlog,'(1x,3es22.10)') Force(:,i)
       end do
    end if
  END SUBROUTINE write_etot_force_log

   
END MODULE atomopt_module
