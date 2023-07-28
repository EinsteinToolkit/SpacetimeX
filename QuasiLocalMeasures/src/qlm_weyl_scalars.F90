#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_calc_weyl_scalars (CCTK_ARGUMENTS, hn)
  use qlm_boundary
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn

  call qlm_calc_weyl_scalars1 (CCTK_PASS_FTOF, hn)

  call set_boundary (CCTK_PASS_FTOF, hn, qlm_psi0(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_psi1(:,:,hn), -1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_psi2(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_psi3(:,:,hn), -1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_psi4(:,:,hn), +1)
  
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_i(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_j(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_s(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_sdiff(:,:,hn), +1)
  
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_phi00(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_phi11(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_phi01(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_phi12(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_phi10(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_phi21(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_phi02(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_phi22(:,:,hn), +1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_phi20(:,:,hn), +1)
  
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_lambda(:,:,hn), +1)
  
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_lie_n_theta_l(:,:,hn), +1)
 
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_rsc(:,:,hn), +1)
  
end subroutine qlm_calc_weyl_scalars
