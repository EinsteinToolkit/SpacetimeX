! $Header$

#include "cctk.h"

module classify
  implicit none
  
  interface
     elemental integer function TAT_isnan (x)
       implicit none
       CCTK_REAL, intent(in) :: x
     end function TAT_isnan
  end interface
  
  interface
     elemental integer function TAT_isinf (x)
       implicit none
       CCTK_REAL, intent(in) :: x
     end function TAT_isinf
  end interface
  
  interface
     elemental integer function TAT_finite (x)
       implicit none
       CCTK_REAL, intent(in) :: x
     end function TAT_finite
  end interface
  
  interface
     pure CCTK_REAL function TAT_nan ()
       implicit none
     end function TAT_nan
  end interface
  
  interface
     pure CCTK_REAL function TAT_inf ()
       implicit none
     end function TAT_inf
  end interface
  
end module classify
