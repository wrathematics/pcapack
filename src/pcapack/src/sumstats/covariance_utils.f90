! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module covariance_utils
  implicit none
  
  
  contains
  
  
  ! make symmetric via copying from one triangle to the other.
  ! triang = 'u':  copy from upper; 'l': copy from lower.
  subroutine symmetrize(triang, m, n, x)
    ! in/out
    integer             m, n
    double precision    x(m, n)
    character*1           triang
    ! local
    integer             i, j, k
    
    
    k = min(m, n)
    
    ! copy upper from lower
    if (triang == 'l') then
      do j = 1, n
        do i = j+1, m
          x(j, i) = x(i, j)
        end do
      end do 
    ! copy lower from upper
    else if (triang == 'u') then
      do j = 1, n
        do i = 1, j-1
          x(j, i) = x(i, j)
        end do
      end do 
    else
      stop
    end if
    
    return
  end
  
  
  
  ! x^t * x or x * x^t
  ! trans = 't' :  x^t*x, trans = 'n' : x*x^t
  subroutine crossprod(trans, m, n, alpha, x, c)
    ! in/out
    integer             m, n
    double precision    x( * ), c( * ), alpha
    character*1         trans
    ! local
    integer             ldx, ldc
    character*1         nst
    ! external
    external            dsyrk
    
    
    if (trans == 't') then
      nst = 'n'
      ldx = n
      ldc = m
    else 
      nst = 't'
      ldx = m
      ldc = n
    end if
    
    ! compute upper triangle of x^t*x or x^t*x
    call dsyrk('u', nst, ldc, ldx, alpha, x, m, 0.0d0, c, ldc)
    
    ! fill lower triangle (make symmetric)
    call symmetrize('u', ldc, ldc, c)
    
    return
  end 
  
  
  
  ! matrix inverse
  subroutine inverse(n, x, info)
    ! in/out
    integer             n, info
    double precision    x( * )
    ! local
    integer             lwork, allocerr
    integer, allocatable :: ipiv(:)
    double precision    tmp
    double precision, allocatable :: work(:)
    ! external
    external           dgetrf, dgetri
    
    
    ! factor x=lu
    allocerr = 0
    allocate(ipiv(n), stat=allocerr)
    if (allocerr /= 0) stop "out of memory"
    
    call dgetrf(n, n, x, n, ipiv, info)
    
    if (info /= 0) return
    
    ! invert x
    lwork = -1
    
    call dgetri(n, x, n, ipiv, tmp, lwork, info)
    if (info /= 0) return
    
    lwork = int(tmp)
    allocate(work(lwork), stat=allocerr)
    if (allocerr /= 0) stop "out of memory"
    
    call dgetri(n, x, n, ipiv, work, lwork, info)
    
    1 continue
    deallocate(work)
    deallocate(ipiv)
    
    return
  end
  
end module
