subroutine print_matrix(m, n, a)
  integer, intent(in) :: m, n
  double precision, intent(in) :: a(m, n)
  integer :: i, j
  
  
  do i = 1, m
    do j = 1, n
      write(*, '(1000F14.2)', advance="no")( real(a(i, j)) )
    end do
    print *, ""
  end do
  
  return
end subroutine



program main
  use :: covariance
  
  integer, parameter :: m = 10, n = 3
  double precision :: x(m, n)
  integer :: i, j
  double precision :: ret(n, n)
  
  
  do j = 1, n
    do i = 1, m
      x(i, j) = i + m*(j-1)
    end do
  end do
  
  call cov(m, n, x, ret)
  
  call print_matrix(n, n, ret)
end program
