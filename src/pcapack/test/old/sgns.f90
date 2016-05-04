program main
  use :: sgns
  
  integer :: x = 5, y = -5
  double precision :: u = 2.345d0, v = -2.345d0
  
  print *, sgn(x) == 1, sgn(y) == -1, sgn(u) == 1, sgn(v) == -1
end program
