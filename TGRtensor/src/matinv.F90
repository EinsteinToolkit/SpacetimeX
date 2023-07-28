! $Header$

#include "cctk.h"

module matinv
  use constants
  use lapack
  implicit none
  private
  
  integer, parameter :: lik = lapack_integer_kind

  public calc_inv2
  
  public calc_posinv3
  public calc_syminv3
  public calc_inv3
  
  public calc_syminv4
  public calc_inv4
  
contains
  
  subroutine calc_inv2 (g2, gu2, lerr)
    CCTK_REAL,           intent(in)  :: g2(2,2)
    CCTK_REAL,           intent(out) :: gu2(2,2)
    logical,   optional, intent(out) :: lerr
    CCTK_REAL    :: tmp(2,2)
    integer(lik) :: ipiv(2)
    integer(lik) :: info
    character    :: msg*100
    
    tmp = g2
    gu2 = delta2
    call gesv (2_lik, 2_lik, tmp, 2_lik, ipiv, gu2, 2_lik, info)
    
    if (.not. present(lerr) .and. info /= 0) then
       write (msg, '("Error in call to GESV, info=",i4)') info
       call CCTK_WARN (1, msg)
    end if
    
    if (present(lerr)) lerr = info /= 0
  end subroutine calc_inv2
  
  
  
  subroutine calc_posinv3 (g3, gu3, lerr)
    CCTK_REAL,           intent(in)  :: g3(3,3)
    CCTK_REAL,           intent(out) :: gu3(3,3)
    logical,   optional, intent(out) :: lerr
    CCTK_REAL    :: tmp(3,3)
    integer(lik) :: info
    character    :: msg*100
    
    tmp = g3
    gu3 = delta3
    call posv ('u', 3_lik, 3_lik, tmp, 3_lik, gu3, 3_lik, info)
    
    if (.not. present(lerr) .and. info /= 0) then
       write (msg, '("Error in call to POSV, info=",i4)') info
       call CCTK_WARN (1, msg)
    end if
    
    if (present(lerr)) lerr = info /= 0
  end subroutine calc_posinv3
  
  subroutine calc_syminv3 (g3, gu3, lerr)
    CCTK_REAL,           intent(in)  :: g3(3,3)
    CCTK_REAL,           intent(out) :: gu3(3,3)
    logical,   optional, intent(out) :: lerr
    CCTK_REAL    :: tmp(3,3)
    integer(lik) :: ipiv(3)
    integer(lik) :: info
    character    :: msg*100
    
    integer(lik), parameter :: lwork = 100
    CCTK_REAL :: work(lwork)
    
    tmp = g3
    gu3 = delta3
    call sysv ('u', 3_lik, 3_lik, tmp, 3_lik, ipiv, gu3, 3_lik, work, lwork, &
         info)
    
    if (.not. present(lerr) .and. info /= 0) then
       write (msg, '("Error in call to SYSV, info=",i4)') info
       call CCTK_WARN (1, msg)
    end if
    
    if (present(lerr)) lerr = info /= 0
  end subroutine calc_syminv3
  
  subroutine calc_inv3 (g3, gu3, lerr)
    CCTK_REAL,           intent(in)  :: g3(3,3)
    CCTK_REAL,           intent(out) :: gu3(3,3)
    logical,   optional, intent(out) :: lerr
    CCTK_REAL    :: tmp(3,3)
    integer(lik) :: ipiv(3)
    integer(lik) :: info
    character    :: msg*100
    
    tmp = g3
    gu3 = delta3
    call gesv (3_lik, 3_lik, tmp, 3_lik, ipiv, gu3, 3_lik, info)
    
    if (.not. present(lerr) .and. info /= 0) then
       write (msg, '("Error in call to GESV, info=",i4)') info
       call CCTK_WARN (1, msg)
    end if
    
    if (present(lerr)) lerr = info /= 0
  end subroutine calc_inv3
  
  
  
  subroutine calc_syminv4 (g4, gu4, lerr)
    CCTK_REAL,           intent(in)  :: g4(0:3,0:3)
    CCTK_REAL,           intent(out) :: gu4(0:3,0:3)
    logical,   optional, intent(out) :: lerr
    CCTK_REAL    :: tmp(0:3,0:3)
    integer(lik) :: ipiv(0:3)
    integer(lik) :: info
    character    :: msg*100
    
    integer(lik), parameter :: lwork = 100
    CCTK_REAL :: work(lwork)
    
    tmp = g4
    gu4 = delta4
    call sysv ('u', 4_lik, 4_lik, tmp, 4_lik, ipiv, gu4, 4_lik, work, lwork, &
         info)
    
    if (.not. present(lerr) .and. info /= 0) then
       write (msg, '("Error in call to SYSV, info=",i4)') info
       call CCTK_WARN (1, msg)
    end if
    
    if (present(lerr)) lerr = info /= 0
  end subroutine calc_syminv4
  
  subroutine calc_inv4 (g4, gu4, lerr)
    CCTK_REAL,           intent(in)  :: g4(0:3,0:3)
    CCTK_REAL,           intent(out) :: gu4(0:3,0:3)
    logical,   optional, intent(out) :: lerr
    CCTK_REAL    :: tmp(0:3,0:3)
    integer(lik) :: ipiv(0:3)
    integer(lik) :: info
    character    :: msg*100
    
    tmp = g4
    gu4 = delta4
    call gesv (4_lik, 4_lik, tmp, 4_lik, ipiv, gu4, 4_lik, info)
    
    if (.not. present(lerr) .and. info /= 0) then
       write (msg, '("Error in call to GESV, info=",i4)') info
       call CCTK_WARN (1, msg)
    end if
    
    if (present(lerr)) lerr = info /= 0
  end subroutine calc_inv4
  
end module matinv
