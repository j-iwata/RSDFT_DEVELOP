!--------1---------2---------3---------4---------5---------6---------7---------8---------9

      SUBROUTINE init_wf_mol_sym
      use global_variables
      implicit none
      integer :: n,k,i,n1,n2,i1,i2,i3,s,m,j,ir,id,ie,jd
      integer :: lt(3),isym,idum,ierr
      real(8) :: c,ctime0,ctime1,etime0,etime1,mem,memax
      real(8),allocatable :: w(:,:,:,:),w0(:,:,:)

      INTERFACE
         FUNCTION ran0(idum)
         real(8) :: ran0
         integer,intent(INOUT) :: idum
         END FUNCTION ran0
      END INTERFACE

      call bwatch(ctime0,etime0)
      if (DISP_SWITCH) write(*,'(a60," init_wf_mol_sym")') repeat("-",60)

      n1=idisp(myrank)+1
      n2=idisp(myrank)+ircnt(myrank)
      mem=0.d0
      memax=0.d0

      unk(:,:,:,:)=zero

      allocate( w0(-ML1:ML1,-ML2:ML2,-ML3:ML3) ) ; w0=0.d0
      mem=mem+bdreal*size(w0) ; memax=max(mem,memax)

      idum=456214531

      do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1

         do ir=1,Nir
            c=dble(irdim(ir))/dble(nsym)
            N_excited=iedim(ir)
            do ie=1,N_excited
!               do i=1,ML_irreducible
!               do isym=1,nsym
!                  lt(1:3)=matmul( rga(1:3,1:3,isym),LL(1:3,i) )
!                  w0(lt(1),lt(2),lt(3))=ran0(idum)-0.5d0
!               end do
!               end do
               do i3=-ML3,ML3
               do i2=-ML2,ML2
               do i1=-ML1,ML1
                  w0(i1,i2,i3)=ran0(idum)-0.5d0
               end do
               end do
               end do
               do id=1,irdim(ir)
                  n=irlabel2n(ir,id,ie) ; if ( n<MB_0 .or. MB_1<n ) cycle
                  do i=n1,n2
                     i1=LL(1,i) ; i2=LL(2,i) ; i3=LL(3,i)
                     do isym=1,nsym
                        lt(1)=rga(1,1,isym)*i1+rga(2,1,isym)*i2+rga(3,1,isym)*i3
                        lt(2)=rga(1,2,isym)*i1+rga(2,2,isym)*i2+rga(3,2,isym)*i3
                        lt(3)=rga(1,3,isym)*i1+rga(2,3,isym)*i2+rga(3,3,isym)*i3
                        unk(i,n,k,s)=unk(i,n,k,s)+c*Rir(id,1,isym,ir)*w0( lt(1),lt(2),lt(3) )
                     end do
                  end do
               end do
            end do ! ie
         end do ! ir

      end do ! k
      end do ! s

      mem=mem-bdreal*size(w0) ; deallocate( w0 )

!check
      goto 100
      allocate( w(-ML1:ML1,-ML2:ML2,-ML3:ML3,MB) ) ; w=0.d0
      mem=mem+bdreal*size(w) ; memax=max(mem,memax)

      do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1

         w(:,:,:,:)=0.d0

         do n=1,MB
            ir=irlabel(1,n)
            id=irlabel(2,n)
            ie=irlabel(3,n)
            do i=n1,n2
               do isym=1,nsym
                  lt(1:3)=matmul( rga(1:3,1:3,isym),LL(1:3,i) )
                  c=0.d0
                  do jd=1,irdim(ir)
                     c=c+Rir(id,jd,isym,ir)*unk(i,irlabel2n(ir,jd,ie),k,s)
                  end do
                  if ( w(lt(1),lt(2),lt(3),n)==0.d0 ) then
                     w(lt(1),lt(2),lt(3),n)=c
                  end if
               end do
            end do
         end do

         if (DISP_SWITCH) then
            write(*,'(1x,2a3,a18,1x,3(a6,1x))') "m","n"," ","irrd","idim","iexc"
            do n=1,MB
               do m=1,n
                  write(*,'(1x,2i3,f18.12,1x,3(2i3,1x),2i8)') &
                       m,n,sum(w(:,:,:,m)*w(:,:,:,n)) &
                       ,irlabel(1,m),irlabel(1,n) &
                       ,irlabel(2,m),irlabel(2,n) &
                       ,irlabel(3,m),irlabel(3,n) &
                       ,count(w(:,:,:,m)/=0.d0),count(w(:,:,:,n)/=0.d0)
               end do
            end do
         end if

      end do ! k
      end do ! s

      mem=mem-bdreal*size(w) ; deallocate( w )
 100  continue

      call bwatch(ctime1,etime1)
      if ( DISP_SWITCH ) then
         write(*,*) "TIME(INIT_WF_MOL_SYM)=",ctime1-ctime0,etime1-etime0
         write(*,*) "MEM(MB)",memax*B2MB,mem
      end if

      return

 900  call stop_program

      END SUBROUTINE init_wf_mol_sym

