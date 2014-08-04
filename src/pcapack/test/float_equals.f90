program main
  use :: float_equals
  logical :: ret
  
  
  ret = fequals(2.345d0, 2.34500000000001d0)
  print *, ret
  
  ret = fequals(2.345d0, 2.34500001d0)
  print *, ret
  
  ret = fequals(0.0d0, -0.0d0)
  print *, ret
end program
