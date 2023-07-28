#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_calc_coordinates (CCTK_ARGUMENTS, hn)
  use cctk
  use constants
  use qlm_boundary
  use tensor2
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_REAL, parameter :: one=1, two=2, half=one/two
  CCTK_REAL :: z0, z0dot, z1, z1dot
  CCTK_REAL :: qq(2,2), dtq
  CCTK_REAL :: integral_z, area, radius
  integer   :: i0,j0
  integer   :: i,j
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Finding invariant coordinates")
  end if
  
  ! latitude of "equator"
  i0 = (qlm_ntheta(hn)+1)/2
  ! longitude of zero meridian
  j0 = 1+qlm_nghostsphi(hn)
  
  
  
  ! calculate area
  area = 0
  do j = 1+qlm_nghostsphi(hn), qlm_nphi(hn)-qlm_nghostsphi(hn)
     do i = 1+qlm_nghoststheta(hn), qlm_ntheta(hn)-qlm_nghoststheta(hn)
        
        ! 2-metric on the horizon
        qq(1,1) = qlm_qtt(i,j,hn)
        qq(1,2) = qlm_qtp(i,j,hn)
        qq(2,2) = qlm_qpp(i,j,hn)
        qq(2,1) = qq(1,2)
        
        call calc_2det (qq, dtq)
        
        area = area + sqrt(dtq) * qlm_delta_theta(hn) * qlm_delta_phi(hn)
        
     end do
  end do
  radius = sqrt(area / (4*pi))
  
  
  
  ! initial value
  qlm_inv_z(i0,j0,hn) = 0
  
  ! d_a z = 1/R^2 xi^b eps_ba
  
  ! transport along equator
  do j = j0+1, qlm_nphi(hn)-qlm_nghostsphi(hn)
     i = i0
     
     z0 = qlm_inv_z(i,j-1,hn)
     z0dot = rhs(i,j-1,(/0,1/))
     z1 = z0 + qlm_delta_phi(hn) * z0dot
     z1dot = rhs(i,j,(/0,1/))
     qlm_inv_z(i,j,hn) = z0 + half * qlm_delta_phi(hn) * (z0dot + z1dot)
  end do
     
  ! transport along meridians
  do j = 1+qlm_nghostsphi(hn), qlm_nphi(hn)-qlm_nghostsphi(hn)
     do i = i0-1, 1+qlm_nghoststheta(hn), -1
        z0 = qlm_inv_z(i+1,j,hn)
        z0dot = rhs(i+1,j,(/-1,0/))
        z1 = z0 + qlm_delta_theta(hn) * z0dot
        z1dot = rhs(i,j,(/-1,0/))
        qlm_inv_z(i,j,hn) = z0 + half * qlm_delta_theta(hn) * (z0dot + z1dot)
     end do
     
     do i = i0+1, qlm_ntheta(hn)-qlm_nghoststheta(hn)
        z0 = qlm_inv_z(i-1,j,hn)
        z0dot = rhs(i-1,j,(/1,0/))
        z1 = z0 + qlm_delta_theta(hn) * z0dot
        z1dot = rhs(i,j,(/1,0/))
        qlm_inv_z(i,j,hn) = z0 + half * qlm_delta_theta(hn) * (z0dot + z1dot)
     end do
  end do
  
  ! normalise
  integral_z = 0
  do j = 1+qlm_nghostsphi(hn), qlm_nphi(hn)-qlm_nghostsphi(hn)
     do i = 1+qlm_nghoststheta(hn), qlm_ntheta(hn)-qlm_nghoststheta(hn)
        
        ! 2-metric on the horizon
        qq(1,1) = qlm_qtt(i,j,hn)
        qq(1,2) = qlm_qtp(i,j,hn)
        qq(2,2) = qlm_qpp(i,j,hn)
        qq(2,1) = qq(1,2)
        
        call calc_2det (qq, dtq)
        
        integral_z = integral_z + qlm_inv_z(i,j,hn) * sqrt(dtq) * qlm_delta_theta(hn) * qlm_delta_phi(hn)
        
     end do
  end do
  qlm_inv_z(:,:,hn) = qlm_inv_z(:,:,hn) - integral_z / area
  
  ! boundary conditions
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_inv_z(:,:,hn), +1)
  
  
  
#if 0
  
  ! initial value
  qlm_inv_phi(i0,j0,hn) = 0
  
  ! xi^a d_a phi = C
  ! z^a d_a phi = 0
  ! z^a = (R^4 / q_bc xi^b xi^c) q^ab d_b z
  ! v^a = A z^a + B xi^a
  ! v^a d_a phi = (A z^a + B xi^a) d_a phi
  !             = B C
  ! (choose C=1, and normalise later)
  
#error "replace z by phi"
  ! transport along equator
  do j = j0+1, qlm_nphi(hn)-qlm_nghostsphi(hn)
     i = i0
     
     z0 = qlm_inv_z(i,j-1,hn)
     z0dot = rhs(i,j-1,(/0,1/))
     z1 = z0 + qlm_delta_phi(hn) * z0dot
     z1dot = rhs(i,j,(/0,1/))
     qlm_inv_z(i,j,hn) = z0 + half * qlm_delta_phi(hn) * (z0dot + z1dot)
  end do
     
  ! transport along meridians
  do j = 1+qlm_nghostsphi(hn), qlm_nphi(hn)-qlm_nghostsphi(hn)
     do i = i0-1, 1+qlm_nghoststheta(hn), -1
        z0 = qlm_inv_z(i+1,j,hn)
        z0dot = rhs(i+1,j,(/-1,0/))
        z1 = z0 + qlm_delta_theta(hn) * z0dot
        z1dot = rhs(i,j,(/-1,0/))
        qlm_inv_z(i,j,hn) = z0 + half * qlm_delta_theta(hn) * (z0dot + z1dot)
     end do
     
     do i = i0+1, qlm_ntheta(hn)-qlm_nghoststheta(hn)
        z0 = qlm_inv_z(i-1,j,hn)
        z0dot = rhs(i-1,j,(/1,0/))
        z1 = z0 + qlm_delta_theta(hn) * z0dot
        z1dot = rhs(i,j,(/1,0/))
        qlm_inv_z(i,j,hn) = z0 + half * qlm_delta_theta(hn) * (z0dot + z1dot)
     end do
  end do
  
#error "normalise"
  
#endif
  
contains
  
  function rhs (i, j, vv) result (zdot)
    integer,   intent(in) :: i, j
    integer,   intent(in) :: vv(2)
    CCTK_REAL             :: zdot
    CCTK_REAL :: qq(2,2), dtq
    
    ! 2-metric on the horizon
    qq(1,1) = qlm_qtt(i,j,hn)
    qq(1,2) = qlm_qtp(i,j,hn)
    qq(2,2) = qlm_qpp(i,j,hn)
    qq(2,1) = qq(1,2)
    
    call calc_2det (qq, dtq)
    
    zdot = vv(1) * (- (1/radius**2) * qlm_xi_p(i,j,hn) * sqrt(dtq)) &
         + vv(2) * (+ (1/radius**2) * qlm_xi_t(i,j,hn) * sqrt(dtq))
  end function rhs
  
end subroutine qlm_calc_coordinates
