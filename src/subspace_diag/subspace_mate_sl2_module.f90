module subspace_mate_sl2_module

  implicit none

  private
  public :: subspace_mate_sl2

  real(8) :: dV=1.0d0
  complex(8) :: zdV=(1.0d0,0.0d0)
  logical :: init_done=.false.

  integer,parameter :: nblk_min = 4
  complex(8),allocatable :: uhu(:,:)

contains

  subroutine init( u )
    use rsdft_allreduce_module, only: rsdft_allreduce
    implicit none
    complex(8),intent(in) :: u(:)
    real(8) :: c
    c = sum( abs(u)**2 )
    call rsdft_allreduce( c, 'grid' )
    dV = 1.0d0/c
    zdV = dV
    init_done = .true.
  end subroutine init

  subroutine subspace_mate_sl2( k, s, u )
    use parallel_module, only: myrank, MB_d, get_range_parallel
    use sl_variables, only: sl0, sl1, Hsub
    use hamiltonian_module, only: hamiltonian
    use rsdft_allreduce_module, only: rsdft_allreduce
    use sl_tools_module, only: slinfo, distribute_matrix
    implicit none
    integer,intent(in) :: k, s
    complex(8),intent(in) :: u(:,:)
    complex(8),allocatable :: hu(:,:)
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    integer :: nb,mg
    integer :: nb_0,nb_1,nblk,ns,ne,nn,ms,me,mm,ib1,ib2,jb1,jb2
    complex(8),allocatable :: echk2(:,:)
    type(slinfo) :: sl

    call write_border( 1, ' z_subspace_mate_sl2(start)' )

    if ( .not.init_done ) call init( u(:,1) )

    mg = size( u, 1 )
    nb = sl0%nband

    call get_range_parallel( nb_0, nb_1, 'band' )

    nblk = sl0%nbsize

    allocate( hu(size(u,1),nblk) ); hu=z0

    allocate( echk2(nb,nb) ); echk2=z0

    do ns = nb_0, nb_1, nblk

      ne = min( ns + nblk - 1, nb_1 )
      nn = ne - ns + 1

      do ib1 = ns, ne, MB_d
        ib2 = min( ib1 + MB_d - 1, ne )
        jb1 = ib1 - ns + 1
        jb2 = ib2 - ns + 1
        call hamiltonian( k,s, u(:,ib1:ib2), hu(:,jb1:jb2)  )
      end do

      ! Diagonal part
      !
      ms = ns
      me = ne
      mm = nn

      allocate( uhu(ms:me,ns:ne) ); uhu=z0
      call mate_sub( u, hu, ms,me ,ns,ne ,nblk )

      echk2(ms:me,ns:ne) = uhu
      deallocate( uhu )

      do ms = ne+1, nb, nblk
        me = min( ms + nblk - 1, nb )
        mm = me - ms + 1
        allocate( uhu(ms:me,ns:ne) ); uhu=z0
        call ZGEMM('C','N',mm,nn,mg,zdV,u(1,ms),mg,hu,mg,z0,uhu,mm)
        echk2(ms:me,ns:ne) = uhu
        deallocate( uhu )
      end do !ms

    end do !ns

    call rsdft_allreduce( echk2, 'g' )
    call rsdft_allreduce( echk2, 'b' )

    if ( allocated(Hsub) ) deallocate(Hsub)

    if ( sl1%idiag /= '' ) then

      sl%nprow  = sl1%nprow
      sl%npcol  = sl1%npcol
      sl%myrow  = sl1%myrow
      sl%mycol  = sl1%mycol
      sl%mbsize = sl1%mbsize
      sl%nbsize = sl1%nbsize

      allocate( Hsub(sl1%ldr,sl1%ldc) ); Hsub=z0
      
      call distribute_matrix( sl, echk2, Hsub )
    
    end if

    call write_border( 1, ' z_subspace_mate_sl2(end)' )

  end subroutine subspace_mate_sl2


  recursive subroutine mate_sub(  u, hu, m1,m2, n1,n2, nblk )
    implicit none
    complex(8),intent(in) :: u(:,:)
    complex(8),intent(in) :: hu(:,:)
    integer,intent(in) :: m1,m2,n1,n2,nblk
    integer :: mg,n,ns,ne,nn,m,ms,me,mm,nblkh
    integer :: ib1,ib2,jb1,jb2,ld
    complex(8),parameter :: z0=(0.0d0,0.0d0)

    mg = size( u, 1 )
    ld = size( uhu, 1 )

    do ns = n1, n2, nblk
      ne = min( ns + nblk - 1, n2 )
      nn = ne - ns + 1
      jb1 = ns - n1 + 1
      jb2 = ne - n1 + 1

      do ms = m1, m2, nblk
        me = min( ms + nblk - 1, m2 )
        mm = me - ms + 1
        ib1 = ms - m1 + 1
        ib2 = me - m1 + 1

        if ( ms == ns ) then

          if ( mm <= nblk_min ) then
            do n = ns, ne
              call ZGEMV('C',mg,ne-n+1,zdV,u(1,n),mg,hu(1,n-n1+1),1,z0,uhu(n,n),1)
            end do
          else
            nblkh = max( nblk/2, nblk_min )
            call mate_sub( u, hu(:,jb1:jb2), ms,me, ns,ne, nblkh )
          end if

        else if ( ms > ne ) then

          call ZGEMM('C','N',mm,nn,mg,zdV,u(1,ms),mg,hu(1,jb1),mg,z0,uhu(ms,ns),ld)

        end if

      end do ! ms
    end do ! ns

    return
  end subroutine mate_sub

end module subspace_mate_sl2_module
