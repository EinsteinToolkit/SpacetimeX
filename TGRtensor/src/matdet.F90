! $Header$

#include "cctk.h"

module matdet
  use lapack
  implicit none
  private
  
  integer, parameter :: lik = lapack_integer_kind
  
  public calc_symdet3
  public calc_det3
  
  public calc_symdet4
  
contains
  
  subroutine calc_symdet3 (gg, dtg, lerr)
    CCTK_REAL,           intent(in)  :: gg(3,3)
    CCTK_REAL,           intent(out) :: dtg
    logical,   optional, intent(out) :: lerr
    CCTK_REAL    :: tmp(3,3)
    integer(lik) :: ipiv(3)
    integer(lik) :: info
    logical      :: odd
    integer      :: i
    character    :: msg*100
    
    integer(lik), parameter :: lwork = 100
    CCTK_REAL :: work(lwork)
    
    tmp = gg
    call sytrf ('u', 3_lik, tmp, 3_lik, ipiv, work, lwork, info)
    
    if (info < 0) then
       write (msg, '("Error in call to SYTRF, info=",i4)') info
       call CCTK_WARN (1, msg)
    end if
    
    if (present(lerr)) lerr = info /= 0
    
    if (info > 0) then
       dtg = 0
       return
    end if
    
    odd = .false.
    do i=1,3
       if (ipiv(i) /= i) odd = .not. odd
    end do
    
    dtg = 1
    if (odd) dtg = -dtg
    do i=1,3
       dtg = dtg * tmp(i,i)
    end do
  end subroutine calc_symdet3
  
  
  
  subroutine calc_det3 (gg, dtg, lerr)
    CCTK_REAL,           intent(in)  :: gg(3,3)
    CCTK_REAL,           intent(out) :: dtg
    logical,   optional, intent(out) :: lerr
    CCTK_REAL    :: tmp(3,3)
    integer(lik) :: ipiv(3)
    integer(lik) :: info
    logical      :: odd
    integer      :: i
    character    :: msg*100
    
    tmp = gg
    call getrf (3_lik, 3_lik, tmp, 3_lik, ipiv, info)
    
    if (info < 0) then
       write (msg, '("Error in call to GETRF, info=",i4)') info
       call CCTK_WARN (1, msg)
    end if
    
    if (present(lerr)) lerr = info /= 0
    
    if (info > 0) then
       dtg = 0
       return
    end if
    
    odd = .false.
    do i=1,3
       if (ipiv(i) /= i) odd = .not. odd
    end do
    
    dtg = 1
    if (odd) dtg = -dtg
    do i=1,3
       dtg = dtg * tmp(i,i)
    end do
  end subroutine calc_det3
  
  
  
  subroutine calc_symdet4 (g4, dtg4, lerr)
    CCTK_REAL,           intent(in)  :: g4(4,4)
    CCTK_REAL,           intent(out) :: dtg4
    logical,   optional, intent(out) :: lerr
    CCTK_REAL    :: tmp(4,4)
    integer(lik) :: ipiv(4)
    integer(lik) :: info
    logical      :: odd
    integer      :: i
    character    :: msg*100
    
    integer(lik), parameter :: lwork = 100
    CCTK_REAL :: work(lwork)

    tmp = g4
    call sytrf ('u', 4_lik, tmp, 4_lik, ipiv, work, lwork, info)
    
    if (info < 0) then
       write (msg, '("Error in call to SYTRF, info=",i4)') info
       call CCTK_WARN (1, msg)
    end if
    
    if (present(lerr)) lerr = info /= 0
    
    if (info > 0) then
       dtg4 = 0
       return
    end if
    
    odd = .false.
    do i=1,4
       if (ipiv(i) /= i) odd = .not. odd
    end do
    
    dtg4 = 1
    if (odd) dtg4 = -dtg4
    do i=1,4
       dtg4 = dtg4 * tmp(i,i)
    end do
  end subroutine calc_symdet4
  
end module matdet
