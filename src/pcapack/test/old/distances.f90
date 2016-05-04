program main
  use :: distances
  use :: float_equals
  integer :: i
  double precision :: x(4) = [1.0d0, 1.2d0, 3.5d0, 2.7d0]
  double precision :: y(4) = [-1.0d0, 2.4d0, 1.5d0, 2.0d0]
  double precision :: truth(6) = [3.1511902513177463d0, 2.0d0, 0.7d0, 5.9d0, 0.0d0, 2.3161589374533507d0]
  character(len=9) :: dist(6) = [character(len=9) :: "euclidean", "supremum", "infimum", "manhattan", "canberra", "minkowski"]
  double precision :: d
  
  
  truth(5) = setinf()
  
  do i = 1, 6
    d = distance(i, 4, x, y, 5)
    print *, dist(i), "=", d, truth(i), fequals(d, truth(i))
  end do
  
  
  contains
  
  function setinf() result(inf)
    double precision :: inf, x
    x = huge(1.0d0)
    inf = x + x
  end function
end program
