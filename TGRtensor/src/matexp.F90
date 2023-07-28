! $Header$

#include "cctk.h"

module matexp
  use constants
  implicit none
  private
  
  public calc_exp3
  
contains
  
  subroutine calc_exp3 (h3, g3)
    CCTK_REAL, intent(in)  :: h3(3,3)
    CCTK_REAL, intent(out) :: g3(3,3)
    CCTK_REAL :: nfact
    CCTK_REAL :: tmp(3,3)
    integer   :: n
    
    g3 = delta3
    
    nfact = 1
    tmp = delta3
    do n = 1, 18
       nfact = nfact * n
       tmp = matmul (h3, tmp)
       
       ! exp(x) = sum_n x^n / n!
       g3 = g3 + tmp / nfact
    end do
    
  end subroutine calc_exp3
  
end module matexp
