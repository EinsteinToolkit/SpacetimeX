! $Header$

#include "cctk.h"

module conversion
  implicit none
  private
  
  public make_cylindrical2cartesian
  public make_spherical2cartesian
  public make_cartesian2spherical
  
contains
  
  subroutine make_cylindrical2cartesian (xx, tt)
    CCTK_REAL, intent(in)  :: xx(3)
    CCTK_REAL, intent(out) :: tt(3,3)
    CCTK_REAL, parameter :: eps = 1.0d-20
    CCTK_REAL :: x,y,z
    CCTK_REAL :: rho
    
    ! cylindrical coordinates: rho phi z
    ! Cartesian   coordinates: x   y   z
    !    rho^2    = x^2 + y^2
    !    tan(phi) = y/x
    
    ! A[cart](p)_i = TT(p)_i^j A[cyl](p)_j
    !    where p = p[cart]
    
    ! Transformation:
    !    x = rho cos(phi)
    !    y = rho sin(phi)
    !    z = z
    !    rho = sqrt(x^2+y^2)
    !    phi = atan(y/x)
    !    z   = z
    
    ! Jacobian:
    !    drho/dx = x/rho
    !    drho/dy = y/rho
    !    drho/dz = 0
    !    dphi/dx = -y/rho^2
    !    dphi/dy = x/rho^2
    !    dphi/dz = 0
    !    dz  /dx = 0
    !    dz  /dy = 0
    !    dz  /dz = 1
    
    x = xx(1)
    y = xx(2)
    z = xx(3)
    
    rho = sqrt(x**2+y**2)
    
    tt(1,1) = x/(rho+eps)
    tt(2,1) = y/(rho+eps)
    tt(3,1) = 0
    tt(1,2) = -y/(rho**2+eps)
    tt(2,2) = x/(rho**2+eps)
    tt(3,2) = 0
    tt(1,3) = 0
    tt(2,3) = 0
    tt(3,3) = 1
    
  end subroutine make_cylindrical2cartesian
  
  subroutine make_spherical2cartesian (xx, tt)
    CCTK_REAL, intent(in)  :: xx(3)
    CCTK_REAL, intent(out) :: tt(3,3)
    CCTK_REAL, parameter :: eps = 1.0d-20
    CCTK_REAL :: x,y,z
    CCTK_REAL :: rho,r
    
    ! Cartesian coordinates: x y z
    ! spherical coordinates: r theta phi
    !    r^2        = x^2 + y^2 + z^2
    !    cos(theta) = z/r
    !    tan(phi)   = y/x
    
    ! A[cart](p)_i = TT(p)_i^j A[spher](p)_j
    !    where p = p[cart]
    
    ! Definition:
    !    rho   = r sin(theta)
    !    z/rho = cos(theta) / sin(theta)
    
    ! Transformation:
    !    x = r sin(theta) cos(phi)
    !    y = r sin(theta) sin(phi)
    !    z = r cos(theta)
    !    r     = sqrt(x^2+y^2+z^2)
    !    theta = acos(z/r)
    !    phi   = atan(y/x)
    
    ! Jacobian:
    !    dr    /dx =  x/r
    !    dr    /dy =  y/r
    !    dr    /dz =  z/r
    !    dtheta/dx =  xz/r^2rho           ! 0
    !    dtheta/dy =  yz/r^2rho           ! 0
    !    dtheta/dz = -rho/r^2             ! 1/rho
    !    dphi  /dx = -y/rho^2
    !    dphi  /dy =  x/rho^2
    !    dphi  /dz =  0
    
    x = xx(1)
    y = xx(2)
    z = xx(3)
    
    rho = sqrt(x**2+y**2)
    r = sqrt(x**2+y**2+z**2)
    
    tt(1,1) =  x/(r+eps)
    tt(2,1) =  y/(r+eps)
    tt(3,1) =  z/(r+eps)
    tt(1,2) =  x*z/(r**2*rho+eps)
    tt(2,2) =  y*z/(r**2*rho+eps)
    tt(3,2) = -rho**2/(r**2*rho+eps)
    tt(1,3) = -y/(rho**2+eps)
    tt(2,3) =  x/(rho**2+eps)
    tt(3,3) =  0
    
  end subroutine make_spherical2cartesian
  
  subroutine make_cartesian2spherical (xx, tt)
    CCTK_REAL, intent(in)  :: xx(3)
    CCTK_REAL, intent(out) :: tt(3,3)
    CCTK_REAL, parameter :: eps = 1.0d-20
    CCTK_REAL :: x,y,z
    CCTK_REAL :: rho,r
    
    ! Cartesian coordinates: x y z
    ! spherical coordinates: r theta phi
    !    r^2        = x^2 + y^2 + z^2
    !    cos(theta) = z/r
    !    tan(phi)   = y/x
    
    ! A[spher](p)_i = TT(p)_i^j A[cart](p)_j
    !    where p = p[cart]
    
    ! Definition:
    !    rho   = r sin(theta)
    !    z/rho = cos(theta) / sin(theta)
    
    ! Transformation:
    !    x = r sin(theta) cos(phi)
    !    y = r sin(theta) sin(phi)
    !    z = r cos(theta)
    !    r     = sqrt(x^2+y^2+z^2)
    !    theta = acos(z/r)
    !    phi   = atan(y/x)
    
    ! Jacobian:
    !    dx/dr     =     sin(theta) cos(phi) =  x/r
    !    dx/dtheta =   r cos(theta) cos(phi) =  x z/rho
    !    dx/dphi   = - r sin(theta) sin(phi) = -y
    !    dy/dr     =     sin(theta) sin(phi) =  y/r
    !    dy/dtheta =   r cos(theta) sin(phi) =  y z/rho
    !    dy/dphi   =   r sin(theta) cos(phi) =  x
    !    dz/dr     =     cos(theta)          =  z/r
    !    dz/dtheta = - r sin(theta)          = -rho
    !    dz/dphi   =   0                     =  0
    
    x = xx(1)
    y = xx(2)
    z = xx(3)
    
    rho = sqrt(x**2+y**2)
    r = sqrt(x**2+y**2+z**2)
    
    tt(1,1) =  x/(r+eps)
    tt(2,1) =  x * z/(rho+eps)
    tt(3,1) = -y
    tt(1,2) =  y/(r+eps)
    tt(2,2) =  y * z/(rho+eps)
    tt(3,2) =  x
    tt(1,3) =  z/(r+eps)
    tt(2,3) = -rho
    tt(3,3) =  0
    
  end subroutine make_cartesian2spherical
  
end module conversion
