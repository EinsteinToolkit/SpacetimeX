#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_killing_axial (CCTK_ARGUMENTS, hn)
  use cctk
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_REAL :: xi(2), chi
  integer   :: i, j
  
  do j = 1, qlm_nphi(hn)
     do i = 1, qlm_ntheta(hn)
        
        xi(1) = 0
        xi(2) = 1
        
        chi = 0
        
        qlm_xi_t(i,j,hn) = xi(1)
        qlm_xi_p(i,j,hn) = xi(2)
        
        qlm_chi(i,j,hn) = chi
        
     end do
  end do
  
end subroutine qlm_killing_axial
