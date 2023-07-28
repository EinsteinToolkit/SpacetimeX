#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_multipoles (CCTK_ARGUMENTS, hn)
  use cctk
  use constants
  use qlm_boundary
  use qlm_derivs
  use qlm_variables
  use ricci2
  use tensor2
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_REAL, parameter :: one=1, two=2
  CCTK_REAL, parameter :: half=one/2, fourth=one/4, eighth=one/8
  CCTK_REAL, parameter :: sixteenth=one/16
  CCTK_REAL, parameter :: o128=one/128
  
  CCTK_REAL    :: qq(2,2), dtq, rsc
  CCTK_COMPLEX :: psi2
  CCTK_REAL    :: zz
  CCTK_REAL    :: area, mass, spin
  
  CCTK_REAL    :: delta_space(2)
  
  integer      :: i, j
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Calculating multipole moments")
  end if
  
  delta_space(:) = (/ qlm_delta_theta(hn), qlm_delta_phi(hn) /)
  
  qlm_mp_m0(hn) = 0
  qlm_mp_m1(hn) = 0
  qlm_mp_m2(hn) = 0
  qlm_mp_m3(hn) = 0
  qlm_mp_m4(hn) = 0
  qlm_mp_m5(hn) = 0
  qlm_mp_m6(hn) = 0
  qlm_mp_m7(hn) = 0
  qlm_mp_m8(hn) = 0
  
  qlm_mp_j0(hn) = 0
  qlm_mp_j1(hn) = 0
  qlm_mp_j2(hn) = 0
  qlm_mp_j3(hn) = 0
  qlm_mp_j4(hn) = 0
  qlm_mp_j5(hn) = 0
  qlm_mp_j6(hn) = 0
  qlm_mp_j7(hn) = 0
  qlm_mp_j8(hn) = 0
  
  do j = 1+qlm_nghostsphi(hn), qlm_nphi(hn)-qlm_nghostsphi(hn)
     do i = 1+qlm_nghoststheta(hn), qlm_ntheta(hn)-qlm_nghoststheta(hn)
        
        ! 2-metric on the horizon
        qq(1,1) = qlm_qtt(i,j,hn)
        qq(1,2) = qlm_qtp(i,j,hn)
        qq(2,2) = qlm_qpp(i,j,hn)
        qq(2,1) = qq(1,2)
        
        call calc_2det (qq, dtq)
        
        rsc = qlm_rsc(i,j,hn)
        
        zz = qlm_inv_z(i,j,hn)
        
        area = sqrt(dtq) * qlm_delta_theta(hn) * qlm_delta_phi(hn)
        
        mass = fourth * rsc
        
        qlm_mp_m0(hn) = qlm_mp_m0(hn) + mass * p0(zz) * area
        qlm_mp_m1(hn) = qlm_mp_m1(hn) + mass * p1(zz) * area
        qlm_mp_m2(hn) = qlm_mp_m2(hn) + mass * p2(zz) * area
        qlm_mp_m3(hn) = qlm_mp_m3(hn) + mass * p3(zz) * area
        qlm_mp_m4(hn) = qlm_mp_m4(hn) + mass * p4(zz) * area
        qlm_mp_m5(hn) = qlm_mp_m5(hn) + mass * p5(zz) * area
        qlm_mp_m6(hn) = qlm_mp_m6(hn) + mass * p6(zz) * area
        qlm_mp_m7(hn) = qlm_mp_m7(hn) + mass * p7(zz) * area
        qlm_mp_m8(hn) = qlm_mp_m8(hn) + mass * p8(zz) * area
        
        spin = qlm_spin_density(i,j)
        
        qlm_mp_j0(hn) = qlm_mp_j0(hn) + spin * dp0(zz) * area
        qlm_mp_j1(hn) = qlm_mp_j1(hn) + spin * dp1(zz) * area
        qlm_mp_j2(hn) = qlm_mp_j2(hn) + spin * dp2(zz) * area
        qlm_mp_j3(hn) = qlm_mp_j3(hn) + spin * dp3(zz) * area
        qlm_mp_j4(hn) = qlm_mp_j4(hn) + spin * dp4(zz) * area
        qlm_mp_j5(hn) = qlm_mp_j5(hn) + spin * dp5(zz) * area
        qlm_mp_j6(hn) = qlm_mp_j6(hn) + spin * dp6(zz) * area
        qlm_mp_j7(hn) = qlm_mp_j7(hn) + spin * dp7(zz) * area
        qlm_mp_j8(hn) = qlm_mp_j8(hn) + spin * dp8(zz) * area
        
!!$        spin = 0
!!$        do a=1,2
!!$           do b=1,2
!!$              do c=1,3
!!$                 do d=1,3
!!$                    spin = spin - half * area * epsilon2(a,b) * dzz(b) * ee(d,a) * kk(d,c) * ss(c)
!!$                 end do
!!$              end do
!!$           end do
!!$        end do
!!$        
!!$        qlm_mp_j0(hn) = qlm_mp_j0(hn) + spin * dp0(zz) * area
!!$        qlm_mp_j1(hn) = qlm_mp_j1(hn) + spin * dp1(zz) * area
!!$        qlm_mp_j2(hn) = qlm_mp_j2(hn) + spin * dp2(zz) * area
!!$        qlm_mp_j3(hn) = qlm_mp_j3(hn) + spin * dp3(zz) * area
!!$        qlm_mp_j4(hn) = qlm_mp_j4(hn) + spin * dp4(zz) * area
!!$        qlm_mp_j5(hn) = qlm_mp_j5(hn) + spin * dp5(zz) * area
!!$        qlm_mp_j6(hn) = qlm_mp_j6(hn) + spin * dp6(zz) * area
!!$        qlm_mp_j7(hn) = qlm_mp_j7(hn) + spin * dp7(zz) * area
!!$        qlm_mp_j8(hn) = qlm_mp_j8(hn) + spin * dp8(zz) * area
        
     end do
  end do
  
contains
  
  ! Legendre polynomials
  
  function p0 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: p0
    p0 = 1
  end function p0
  
  function p1 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: p1
    p1 = z
  end function p1
  
  function p2 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: p2
    p2 = 3*half * z**2 - half
  end function p2
  
  function p3 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: p3
    p3 = 5*half * z**3 - 3*half * z
  end function p3
  
  function p4 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: p4
    p4 = 35*eighth * z**4 - 15*fourth * z**2 + 3*eighth
  end function p4
  
  function p5 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: p5
    p5 = 63*eighth * z**5 - 35*fourth * z**3 + 15*eighth * z
  end function p5
  
  function p6 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: p6
    p6 = 231*sixteenth * z**6 - 315*sixteenth * z**4 + 105*sixteenth * z**2 &
         - 5*sixteenth
  end function p6
  
  function p7 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: p7
    p7 = 429*sixteenth * z**7 - 693*sixteenth * z**5 + 315*sixteenth * z**3 &
         - 35*sixteenth * z
  end function p7
  
  function p8 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: p8
    p8 = 6435*o128 * z**8 - 12012*o128 * z**6 + 6930*o128 * z**4 &
         - 1260*o128 * z**2 + 35*o128
  end function p8
  
  ! Derivatives of the Legendre polynomials
  
  function dp0 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: dp0
    dp0 = 0
  end function dp0
  
  function dp1 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: dp1
    dp1 = 1
  end function dp1
  
  function dp2 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: dp2
    dp2 = 3 * z
  end function dp2
  
  function dp3 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: dp3
    dp3 = 15*half * z**2 - 3*half
  end function dp3
  
  function dp4 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: dp4
    dp4 = 35*half * z**3 - 15*half*z
  end function dp4
  
  function dp5 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: dp5
    dp5 = 315*eighth * z**4 - 105*fourth * z**2 + 15*eighth
  end function dp5
  
  function dp6 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: dp6
    dp6 = 693*eighth * z**5 - 315*fourth * z**3 + 105*eighth * z
  end function dp6
  
  function dp7 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: dp7
    dp7 = 3003*sixteenth * z**6 - 3465*sixteenth * z**4 + 945*sixteenth * z**2 &
         - 35*sixteenth
  end function dp7
  
  function dp8 (z)
    CCTK_REAL, intent(in) :: z
    CCTK_REAL             :: dp8
    dp8 = 51480*o128 * z**7 - 72072*o128 * z**5 + 27720*o128 * z**3 &
         - 2520*o128 * z
  end function dp8
  
end subroutine qlm_multipoles



subroutine qlm_multipoles_normalise (CCTK_ARGUMENTS, hn)
  use cctk
  use constants
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Normalising multipole moments")
  end if
  
  ! Normalise
  
!!$  ! This is the normalisation for I_n and L_n
!!$  qlm_mp_m0(hn) = qlm_mp_m0(hn) / sqrt(4*pi/ 1)
!!$  qlm_mp_m1(hn) = qlm_mp_m1(hn) / sqrt(4*pi/ 3)
!!$  qlm_mp_m2(hn) = qlm_mp_m2(hn) / sqrt(4*pi/ 5)
!!$  qlm_mp_m3(hn) = qlm_mp_m3(hn) / sqrt(4*pi/ 7)
!!$  qlm_mp_m4(hn) = qlm_mp_m4(hn) / sqrt(4*pi/ 9)
!!$  qlm_mp_m5(hn) = qlm_mp_m5(hn) / sqrt(4*pi/11)
!!$  qlm_mp_m6(hn) = qlm_mp_m6(hn) / sqrt(4*pi/13)
!!$  qlm_mp_m7(hn) = qlm_mp_m7(hn) / sqrt(4*pi/15)
!!$  qlm_mp_m8(hn) = qlm_mp_m8(hn) / sqrt(4*pi/17)
!!$  
!!$  qlm_mp_j0(hn) = qlm_mp_j0(hn) / sqrt(4*pi/ 1)
!!$  qlm_mp_j1(hn) = qlm_mp_j1(hn) / sqrt(4*pi/ 3)
!!$  qlm_mp_j2(hn) = qlm_mp_j2(hn) / sqrt(4*pi/ 5)
!!$  qlm_mp_j3(hn) = qlm_mp_j3(hn) / sqrt(4*pi/ 7)
!!$  qlm_mp_j4(hn) = qlm_mp_j4(hn) / sqrt(4*pi/ 9)
!!$  qlm_mp_j5(hn) = qlm_mp_j5(hn) / sqrt(4*pi/11)
!!$  qlm_mp_j6(hn) = qlm_mp_j6(hn) / sqrt(4*pi/13)
!!$  qlm_mp_j7(hn) = qlm_mp_j7(hn) / sqrt(4*pi/15)
!!$  qlm_mp_j8(hn) = qlm_mp_j8(hn) / sqrt(4*pi/17)
  
  ! This is the normalisation for M_n and J_n
  qlm_mp_m0(hn) = qlm_mp_m0(hn) * qlm_mass(hn) * qlm_radius(hn)**0 / (2*pi)
  qlm_mp_m1(hn) = qlm_mp_m1(hn) * qlm_mass(hn) * qlm_radius(hn)**1 / (2*pi)
  qlm_mp_m2(hn) = qlm_mp_m2(hn) * qlm_mass(hn) * qlm_radius(hn)**2 / (2*pi)
  qlm_mp_m3(hn) = qlm_mp_m3(hn) * qlm_mass(hn) * qlm_radius(hn)**3 / (2*pi)
  qlm_mp_m4(hn) = qlm_mp_m4(hn) * qlm_mass(hn) * qlm_radius(hn)**4 / (2*pi)
  qlm_mp_m5(hn) = qlm_mp_m5(hn) * qlm_mass(hn) * qlm_radius(hn)**5 / (2*pi)
  qlm_mp_m6(hn) = qlm_mp_m6(hn) * qlm_mass(hn) * qlm_radius(hn)**6 / (2*pi)
  qlm_mp_m7(hn) = qlm_mp_m7(hn) * qlm_mass(hn) * qlm_radius(hn)**7 / (2*pi)
  qlm_mp_m8(hn) = qlm_mp_m8(hn) * qlm_mass(hn) * qlm_radius(hn)**8 / (2*pi)
  
  qlm_mp_j0(hn) = qlm_mp_j0(hn) * qlm_radius(hn)**(-1) / (8*pi)
  qlm_mp_j1(hn) = qlm_mp_j1(hn) * qlm_radius(hn)**0 / (8*pi)
  qlm_mp_j2(hn) = qlm_mp_j2(hn) * qlm_radius(hn)**1 / (8*pi)
  qlm_mp_j3(hn) = qlm_mp_j3(hn) * qlm_radius(hn)**2 / (8*pi)
  qlm_mp_j4(hn) = qlm_mp_j4(hn) * qlm_radius(hn)**3 / (8*pi)
  qlm_mp_j5(hn) = qlm_mp_j5(hn) * qlm_radius(hn)**4 / (8*pi)
  qlm_mp_j6(hn) = qlm_mp_j6(hn) * qlm_radius(hn)**5 / (8*pi)
  qlm_mp_j7(hn) = qlm_mp_j7(hn) * qlm_radius(hn)**6 / (8*pi)
  qlm_mp_j8(hn) = qlm_mp_j8(hn) * qlm_radius(hn)**7 / (8*pi)
  
end subroutine qlm_multipoles_normalise
