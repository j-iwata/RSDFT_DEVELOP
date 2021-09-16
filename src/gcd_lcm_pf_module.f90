module gcd_lcm_pf_module
  
  implicit none

  private
  public :: gcd
  public :: lcm
  public :: prime_factorization

contains

  integer function gcd( m_in, n_in )
    implicit none
    integer, intent(in) :: m_in, n_in
    integer :: m,n,m_tmp,loop
    if ( m_in >= n_in ) then
      m = m_in
      n = n_in
    else
      m = n_in
      n = m_in
    end if
    do
      if ( n == 0 ) exit
      m_tmp = n
      n = mod( m, n )
      m = m_tmp
    end do
    gcd = m
  end function gcd

  integer function lcm( m_in, n_in )
    implicit none
    integer, intent(in) :: m_in, n_in
    lcm = m_in * n_in / gcd( m_in, n_in )
  end function lcm

  subroutine prime_factorization( n_in, ifact )
    implicit none
    integer,intent(in)  :: n_in
    integer,intent(out) :: ifact(:)
    integer :: n,i,j,m
    ifact(:)=0
    ifact(1)=1
    n = n_in
    m = n
    loop_j : do j = 1, n_in
      loop_i : do i = 2, n
        if ( mod(m,i) == 0 ) then
          m = m/i
          ifact(i) = ifact(i) + 1
          exit loop_i
        else if ( i == n ) then
          exit loop_j
        end if
      end do loop_i
      n = m
    end do loop_j
  end subroutine prime_factorization

end module gcd_lcm_pf_module
