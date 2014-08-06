MODULE force_local_mol_module

  use parallel_module
  use atom_module
  use array_bound_module
  use density_module
  use electron_module, only: Nspin
  use bberf_module
  use pseudopot_module, only: parloc, rad, Zps
  use rgrid_mol_module
  use ps_local_mol_module, only: NRcloc, Rcloc, vqls
  use rgrid_module

  implicit none

  PRIVATE
  PUBLIC :: calc_force_local_mol

CONTAINS


  SUBROUTINE calc_force_local_mol(force)
    implicit none
    real(8),intent(INOUT) :: force(:,:)
    integer :: a,ik,i,i1,i2,i3,ir,mm,m1,m2,ierr
    real(8) :: p1,p2,p3,p4,Rc
    real(8) :: Rx,Ry,Rz,x,y,z,r,vx,vy,vz,v0,v
    real(8) :: tmp0,tmp1,const0,const1,pi
    real(8) :: maxerr,err0,err
    real(8),allocatable :: ftmp(:,:),trho(:)
    integer :: M_irad,NRc,ir0
    integer,allocatable :: irad(:,:)

    pi = acos(-1.0d0)
    const0 = 2.0d0/sqrt(pi)
    const1 = 2.0d0/3.0d0*const0

    allocate( ftmp(3,Natom)   ) ; ftmp=0.0d0
    allocate( trho(ML_0:ML_1) ) ; trho=0.0d0

    do i=1,Nspin
       trho(:) = trho(:) + rho(:,i)
    end do

    allocate( irad(0:3000,Nelement) ) ; irad=0
    M_irad=0
    do ik=1,Nelement
       NRc=min( 3000, NRcloc(ik) )
       mm=0
       irad(0,ik)=1
       do ir=1,NRc
          mm=int(100.d0*rad(ir,ik))+1
          irad( mm,ik )=ir
       end do
       ir=irad(0,ik)
       do i=1,mm
          if ( irad(i,ik)==0 ) then
             irad(i,ik)=ir
             cycle
          end if
          ir=irad(i,ik)
       end do
       irad(mm+1:,ik)=ir
       M_irad=max(M_irad,mm)
    end do

    maxerr = 0.0d0

    do a=1,Natom

       ik = ki_atom(a)

       p1 = -Zps(ik)*parloc(1,ik)
       p2 = sqrt(parloc(2,ik))
       p3 = -Zps(ik)*parloc(3,ik)
       p4 = sqrt(parloc(4,ik))

       Rx = aa_atom(1,a)
       Ry = aa_atom(2,a)
       Rz = aa_atom(3,a)

       Rc = Rcloc(ik)

       do i=ML_0,ML_1

          i1 = LL(1,i)
          i2 = LL(2,i)
          i3 = LL(3,i)

          x = i1*Hsize - Rx
          y = i2*Hsize - Ry
          z = i3*Hsize - Rz

          r = sqrt(x*x+y*y+z*z)

          vx = 0.0d0
          vy = 0.0d0
          vz = 0.0d0

          if ( r <= Rc + 1.d-10 ) then

             ir0=irad( int(100.d0*r),ik )
             do ir=ir0,NRc
                if ( r < rad(ir,ik) ) exit
             end do

             if ( ir <= 2 ) then

                r = rad(2,ik)

                err0 = 1.d10
                v0 = 0.0d0
                do mm=1,20
                   m1=max(1,ir-mm)
                   m2=min(ir+mm,NRc)
                   call dpolint(rad(m1,ik),vqls(m1,ik),m2-m1+1,r,v,err)
                   if ( abs(err) < err0 ) then
                      v0=v
                      err0 = abs(err)
                      if ( err0 < 1.d-8 ) exit
                   end if
                end do

                vx = vx + v0*x/r
                vy = vy + v0*y/r
                vz = vz + v0*z/r

             else

                err0 = 1.d10
                v0 = 0.0d0
                do mm=1,20
                   m1=max(1,ir-mm)
                   m2=min(ir+mm,NRc)
                   call dpolint(rad(m1,ik),vqls(m1,ik),m2-m1+1,r,v,err)
                   if ( abs(err) < err0 ) then
                      v0=v
                      err0 = abs(err)
                      if ( err0 < 1.d-8 ) exit
                   end if
                end do

                vx = vx + v0*x/r
                vy = vy + v0*y/r
                vz = vz + v0*z/r

             end if
             maxerr=max(maxerr,err0)

          end if

          if ( r < 1.d-9 ) then

             tmp1 = -const1*(p1*p2**3+p3*p4**3)

             vx = vx + x*tmp1
             vy = vy + y*tmp1
             vz = vz + z*tmp1

          else

             tmp1 = -p1*bberf(p2*r) - p3*bberf(p4*r) &
                  + r*const0*( p1*p2*exp(-p2*p2*r*r) + p3*p4*exp(-p4*p4*r*r) )

             vx = vx + tmp1*x/r**3
             vy = vy + tmp1*y/r**3
             vz = vz + tmp1*z/r**3

          end if

          ftmp(1,a) = ftmp(1,a) + trho(i)*vx
          ftmp(2,a) = ftmp(2,a) + trho(i)*vy
          ftmp(3,a) = ftmp(3,a) + trho(i)*vz

       end do ! i

    end do ! a

    call mpi_allreduce &
         (MPI_IN_PLACE,ftmp,3*Natom,MPI_REAL8,MPI_SUM,comm_grid,ierr)

    force(:,:) = force(:,:) + ftmp(:,:)*dV

    deallocate( irad )
    deallocate( trho )
    deallocate( ftmp )

  END SUBROUTINE calc_force_local_mol


END MODULE force_local_mol_module
