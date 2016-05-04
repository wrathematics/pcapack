program main
  use :: means
  
  double precision :: x(5) = (/1.1d0, 2.2d0, 3.3d0, 4.4d0, 5.5d0/)
  
  print *, mean(5, x)
end program
