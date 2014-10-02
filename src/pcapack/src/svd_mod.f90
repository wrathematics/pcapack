! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module svd_mod
  use :: lapack, only : dgemm, dscal, dlacpy, dgesdd
  use, intrinsic :: iso_c_binding
  use :: cbool
  implicit none
  
  
  contains
  
  ! m >= n:  u(m, n), vt(n, n)
  ! m < n:  u(m, m), vt(m, n)
  ! i.e., u(m, minmn), vt(minmn, n)
  subroutine svd(nu, nv, m, n, x, s, u, vt, info) &
    bind(C, name='LA_svd_')
    ! in/out
!    logical(kind=c_bool), intent(in) :: retu, retvt
    integer, intent(in) :: nu, nv, m, n
    integer, intent(out) :: info
    double precision, intent(in) :: x(m, n)
    double precision, intent(inout) :: s(*), u(*), vt(*)
    ! local
    character(len=1) :: jobz
    integer :: minmn, lwork
    integer, allocatable :: iwork(:)
    double precision :: tmp(1)
    double precision, allocatable :: work(:)
    double precision, allocatable :: cpx(:,:)
     
    
    minmn = min(m, n)
    
    ! allocations
    if (nu == 0 .and. nv == 0) then
      jobz = 'n'
    else if ((nu == 0 .and. m >= n) .or. (nv == 0 .and. m < n)) then
      jobz = 'o'
    else if (nu <= minmn .and. nv <= minmn) then
      jobz = 's'
    else
      jobz = 'a'
    end if
    
    allocate(cpx(m, n))
    cpx = x
    
    allocate(iwork(8*minmn))
    
    lwork = -1
    call dgesdd(jobz, m, n, cpx, m, s, u, m, vt, minmn, tmp, lwork, iwork, info)
    lwork = int(tmp(1))
    allocate(work(lwork))
    
    call dgesdd(jobz, m, n, cpx, m, s, u, m, vt, minmn, work, lwork, iwork, info)
    
    deallocate(iwork)
    deallocate(work)
    deallocate(cpx)
    
    return
  end subroutine
  
end module
