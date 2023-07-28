#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_init (CCTK_ARGUMENTS)
  use cctk
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  integer :: hn
  
  if (verbose/=0 .or. veryverbose/=0) then
     call CCTK_INFO ("Initialising Quasi-Local Measure calculations")
  end if
  
  do hn = 1, num_surfaces
     
     qlm_calc_error(hn) = 1
     qlm_have_valid_data(hn) = 0
     qlm_have_valid_data_p(hn) = 0
     qlm_have_valid_data_p_p(hn) = 0
     qlm_have_killing_vector(hn) = 0
     qlm_have_killing_vector_p(hn) = 0
     qlm_have_killing_vector_p_p(hn) = 0
     qlm_iteration(hn) = -1
     
  end do
  
end subroutine qlm_init
