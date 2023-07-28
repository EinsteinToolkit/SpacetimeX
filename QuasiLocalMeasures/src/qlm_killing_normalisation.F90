#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



module qlm_killing_normalisation
  use cctk
  use constants
  implicit none
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  private
  public killing_factor
  
contains
  
  subroutine killing_factor (CCTK_ARGUMENTS, hn, theta, factor, nsteps)
    DECLARE_CCTK_ARGUMENTS
    integer,   intent(in)  :: hn
    CCTK_REAL, intent(in)  :: theta
    CCTK_REAL, intent(out) :: factor
    integer,   intent(out) :: nsteps
    
    CCTK_REAL :: lambda0, theta0, phi0
    CCTK_REAL :: lambda1, theta1, phi1
    character :: msg*1000
    
    lambda0 = 0
    theta0  = theta
    phi0    = 0
    phi1    = 2*pi
    call killing_geodesic &
         (CCTK_PASS_FTOF, hn, lambda0, theta0, phi0, phi1, lambda1, theta1, &
         nsteps)
    
    factor = (lambda1 - lambda0) / (2*pi)
    
    if (veryverbose/=0) then
       write (msg, '("Integrated at theta=",g16.6)') theta
       call CCTK_INFO (msg)
       write (msg, '("   Integrated in ",i4," steps")') nsteps
       call CCTK_INFO (msg)
       write (msg, '("   Theta error is ",g16.6)') theta1 - theta0
       call CCTK_INFO (msg)
       write (msg, '("   Normalisation factor is ",g16.6)') factor
       call CCTK_INFO (msg)
    end if
  end subroutine killing_factor
  
  
  
  subroutine killing_geodesic &
       (CCTK_ARGUMENTS, hn, lambda0, theta0, phi0, phi1, lambda1, theta1, nsteps)
    DECLARE_CCTK_ARGUMENTS
    integer,   intent(in)  :: hn
    CCTK_REAL, intent(in)  :: lambda0, theta0, phi0, phi1
    CCTK_REAL, intent(out) :: lambda1, theta1
    integer,   intent(out) :: nsteps
    
    CCTK_REAL :: org_theta, org_phi, del_theta, del_phi
    
    CCTK_REAL :: lambda
    CCTK_REAL :: theta, phi
    CCTK_REAL :: dlambda
    CCTK_REAL :: dtheta, dphi
    CCTK_REAL :: theta2, phi2
    
    integer   :: ierr1, ierr2
    
    org_theta = qlm_origin_theta(hn)
    org_phi   = qlm_origin_phi(hn)
    del_theta = qlm_delta_theta(hn)
    del_phi   = qlm_delta_phi(hn)
    
    nsteps = 0
    lambda = lambda0
    theta  = theta0
    phi    = phi0
    
    dtheta = killing_interp (qlm_xi_t(:,:,hn), &
         org_theta, org_phi, del_theta, del_phi, theta, phi, ierr1)
    dphi   = killing_interp (qlm_xi_p(:,:,hn), &
         org_theta, org_phi, del_theta, del_phi, theta, phi, ierr2)
    
    if (ierr1/=0 .or. ierr2/=0) then
       call CCTK_WARN (2, "Integration path leaves the domain")
       lambda1 = lambda0
       theta1  = theta0
       nsteps  = -1
       return
    end if
    
    if (abs(dphi) < 1.0d-8 .or. abs(dtheta) > abs(dphi)) then
       call CCTK_WARN (2, "Integration path starts out too steep")
       lambda1 = lambda0
       theta1  = theta0
       nsteps  = -1
       return
    end if
    
    dlambda = (qlm_delta_phi(hn) / dphi) / 2
    
    do
       
       dtheta = killing_interp (qlm_xi_t(:,:,hn), &
            org_theta, org_phi, del_theta, del_phi, theta, phi, ierr1)
       dphi   = killing_interp (qlm_xi_p(:,:,hn), &
            org_theta, org_phi, del_theta, del_phi, theta, phi, ierr2)
       
       if (ierr1/=0 .or. ierr2/=0) then
          call CCTK_WARN (2, "Integration path leaves the domain")
          lambda1 = lambda0
          theta1  = theta0
          nsteps  = -1
          return
       end if
       
       theta2 = theta + dlambda * dtheta / 2
       phi2   = phi   + dlambda * dphi   / 2
       
       dtheta = killing_interp (qlm_xi_t(:,:,hn), &
            org_theta, org_phi, del_theta, del_phi, theta2, phi2, ierr1)
       dphi   = killing_interp (qlm_xi_p(:,:,hn), &
            org_theta, org_phi, del_theta, del_phi, theta2, phi2, ierr2)
       
       if (ierr1/=0 .or. ierr2/=0) then
          call CCTK_WARN (2, "Integration path leaves the domain")
          lambda1 = lambda0
          theta1  = theta0
          nsteps  = -1
          return
       end if
       
       if (dphi<=0) then
          call CCTK_WARN (2, "Integration path does not enclose the pole")
          lambda1 = lambda0
          theta1  = theta0
          nsteps  = -1
          return
       end if
       
       theta2 = theta + dlambda * dtheta
       phi2   = phi   + dlambda * dphi
       
       if (phi2 >= phi1) exit
       
       if (nsteps > 100000) then
          call CCTK_WARN (2, "Integration takes too many steps")
          lambda1 = lambda0
          theta1  = theta0
          nsteps  = -1
          return
          exit
       end if
       
       nsteps = nsteps + 1
       lambda = lambda + dlambda
       theta  = theta2
       phi    = phi2
       
    end do
    
    dlambda = (phi1 - phi) / dphi
    
    nsteps = nsteps + 1
    lambda = lambda + dlambda
    theta  = theta  + dlambda * dtheta
    phi    = phi    + dlambda * dphi
    
    lambda1 = lambda
    theta1  = theta
    
  end subroutine killing_geodesic
  
  
  
  function killing_interp &
       (array, origin_x1, origin_x2, delta_x1, delta_x2, x1, x2, ierr)
    CCTK_REAL :: killing_interp
    CCTK_REAL, intent(in)  :: array(:,:)
    CCTK_REAL, intent(in)  :: origin_x1, origin_x2
    CCTK_REAL, intent(in)  :: delta_x1, delta_x2
    CCTK_REAL, intent(in)  :: x1, x2
    integer,   intent(out) :: ierr
    CCTK_REAL, parameter :: eps = 1.0d-10
    CCTK_REAL :: xx1, xx2
    CCTK_REAL :: dx1, dx2
    CCTK_REAL :: f1a, f1b, f1c
    CCTK_REAL :: f2a, f2b, f2c
    CCTK_REAL :: interp2a, interp2b, interp2c
    integer   :: i, j
    xx1 = x1
    xx2 = x2
    i = nint((xx1 - origin_x1) / delta_x1) + 1
    j = nint((xx2 - origin_x2) / delta_x2) + 1
    i = max(2, min(size(array,1)-1, i))
    if (j>size(array,2)-1) then
       ! periodicity in phi-direction
       xx2 = xx2 - (size(array,2)-2) * delta_x2
       j = j - size(array,2)-2
    end if
    j = max(2, min(size(array,2)-1, j))
    dx1 = xx1 - (origin_x1 + (i-1) * delta_x1)
    dx2 = xx2 - (origin_x2 + (j-1) * delta_x2)
    if (abs(dx1)>(1+eps)*delta_x1/2 .or. abs(dx2)>(1+eps)*delta_x2/2) then
       call CCTK_WARN (2, "interpolating out of bounds")
       killing_interp = 0
       ierr = 1
       return
    end if
    f1a =            (-dx1) * (1-dx1) / (  2  *   1 )
    f1b = (-1-dx1) *          (1-dx1) / (  1  * (-1))
    f1c = (-1-dx1) * (-dx1)           / ((-1) * (-2))
    f2a =            (-dx2) * (1-dx2) / (  2  *   1 )
    f2b = (-1-dx2) *          (1-dx2) / (  1  * (-1))
    f2c = (-1-dx2) * (-dx2)           / ((-1) * (-2))
    interp2a = f1a * array(i-1,j-1) + f1b * array(i,j-1) + f1c * array(i+1,j-1)
    interp2b = f1a * array(i-1,j  ) + f1b * array(i,j  ) + f1c * array(i+1,j  )
    interp2c = f1a * array(i-1,j+1) + f1b * array(i,j+1) + f1c * array(i+1,j+1)
    killing_interp = f2a * interp2a + f2b * interp2b + f2c * interp2c
    ierr = 0
  end function killing_interp
  
end module qlm_killing_normalisation
