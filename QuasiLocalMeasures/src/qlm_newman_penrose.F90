#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_calc_newman_penrose (CCTK_ARGUMENTS, hn)
  use qlm_boundary
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn

  call qlm_calc_newman_penrose1 (CCTK_PASS_FTOF, hn)

  call set_boundary (CCTK_PASS_FTOF, hn, qlm_npkappa  (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_nptau    (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_npsigma  (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_nprho    (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_npepsilon(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_npgamma  (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_npbeta   (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_npalpha  (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_nppi     (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_npnu     (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_npmu     (:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_nplambda (:,:,hn), +1)
  
end subroutine qlm_calc_newman_penrose
