#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_calc_tetrad (CCTK_ARGUMENTS, hn)
  use qlm_boundary
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn

  call qlm_calc_tetrad1 (CCTK_PASS_FTOF, hn)
  
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_l0(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_l1(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_l2(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_l3(:,:,hn), +1)
  
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_n0(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_n1(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_n2(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_n3(:,:,hn), +1)
  
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_m0(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_m1(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_m2(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_m3(:,:,hn), +1)
  
end subroutine qlm_calc_tetrad
