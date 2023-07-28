! $Header$

#include "cctk.h"

module rotation
  use constants
  implicit none
  private
  
  public make_euler
  
contains
  
  subroutine make_euler1 (phi, tt)
    CCTK_REAL, intent(in)  :: phi
    CCTK_REAL, intent(out) :: tt(0:3,0:3)
    ! y^a = T^a_b x^b
    tt = delta4
    tt(1,1) =  cos(phi)
    tt(1,2) = -sin(phi)
    tt(2,1) =  sin(phi)
    tt(2,2) =  cos(phi)
  end subroutine make_euler1
  
  subroutine make_euler2 (theta, tt)
    CCTK_REAL, intent(in)  :: theta
    CCTK_REAL, intent(out) :: tt(0:3,0:3)
    ! y^a = T^a_b x^b
    tt = delta4
    tt(2,2) =  cos(theta)
    tt(2,3) = -sin(theta)
    tt(3,2) =  sin(theta)
    tt(3,3) =  cos(theta)
  end subroutine make_euler2
  
  subroutine make_euler3 (psi, tt)
    CCTK_REAL, intent(in)  :: psi
    CCTK_REAL, intent(out) :: tt(0:3,0:3)
    ! y^a = T^a_b x^b
    tt = delta4
    tt(1,1) =  cos(psi)
    tt(1,2) = -sin(psi)
    tt(2,1) =  sin(psi)
    tt(2,2) =  cos(psi)
  end subroutine make_euler3
  
  subroutine make_euler (phi, theta, psi, tt)
    CCTK_REAL, intent(in)  :: phi, theta, psi
    CCTK_REAL, intent(out) :: tt(0:3,0:3)
    CCTK_REAL :: tt1(0:3,0:3), tt2(0:3,0:3), tt3(0:3,0:3), tt4(0:3,0:3)
    call make_euler1 (phi, tt1)
    call make_euler2 (theta, tt2)
    call make_euler3 (psi, tt3)
    tt4 = matmul(tt2, tt1)
    tt = matmul(tt3, tt4)
  end subroutine make_euler
  
end module rotation
