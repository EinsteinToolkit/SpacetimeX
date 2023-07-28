#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_analyse (CCTK_ARGUMENTS, hn)
  use adm_metric_simple
  use cctk
  use constants
  use qlm_derivs
  use qlm_variables
  use tensor
  use tensor2
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_REAL, parameter :: zero=0, two=2
  integer,   parameter :: rk = kind(zero)
  
  CCTK_REAL    :: theta, phi
  CCTK_REAL    :: ll(0:3), nn(0:3)
  CCTK_COMPLEX :: mm(0:3)
  CCTK_REAL    :: ss(0:3), tt(0:3)
  CCTK_COMPLEX :: npalpha, npbeta
  CCTK_REAL    :: gg(3,3), dgg(3,3,3), kk(3,3), gg_dot(3,3)
  CCTK_REAL    :: g4(0:3,0:3), dg4(0:3,0:3,0:3)
  CCTK_REAL    :: h4(0:3,0:3), dh4(0:3,0:3,0:3)
  CCTK_REAL    :: dtg, gu(3,3), trk, dgu(3,3,3)
  CCTK_REAL    :: xx(3), ee(3,2)
  CCTK_REAL    :: xi(2), xi1(3)
  CCTK_REAL    :: qq(2,2), dtq
  CCTK_REAL    :: adm_energy, adm_mom(3), adm_amom(3)
  CCTK_REAL    :: w_energy, w_mom(3), w_amom(3,3)
  CCTK_COMPLEX :: ev
  CCTK_REAL    :: spin
  CCTK_REAL    :: npspin
  CCTK_REAL    :: wsspin
  CCTK_REAL    :: tx, ty, tz
  CCTK_REAL    :: xi1_x(3), xi1_y(3), xi1_z(3)
  CCTK_REAL    :: coordspinx, coordspiny, coordspinz
  
  CCTK_REAL    :: delta_space(2)
  
  CCTK_REAL, allocatable  :: weights(:)
  CCTK_REAL               :: sum1

  integer   :: ntheta_inner, nphi_inner

  integer   :: i_eq, j_p0, j_p2
  
  integer   :: i, j, l, m
  integer   :: a, b, c
  
  character :: msg*1000
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Calculating spin")
  end if
  
  delta_space(:) = (/ qlm_delta_theta(hn), qlm_delta_phi(hn) /)
  
  ! Equatorial circumference
  i_eq = (qlm_ntheta(hn) + 1) / 2
  
  ! Polar circumference at phi=0
  j_p0 = 1+qlm_nghostsphi(hn)
  
  ! Polar circumference at phi=pi/2
  j_p2 = 1+qlm_nghostsphi(hn) + (qlm_nphi(hn) - 2*qlm_nghostsphi(hn) - 1) / 4
  
  
  
  ! Initial values
  qlm_equatorial_circumference(hn) = 0
  qlm_polar_circumference_0(hn) = 0
  qlm_polar_circumference_pi_2(hn) = 0
  qlm_area(hn) = 0
  qlm_spin(hn) = 0
  qlm_cvspin(hn) = 0
  qlm_npspin(hn) = 0
  qlm_wsspin(hn) = 0
  qlm_coordspinx(hn) = 0
  qlm_coordspiny(hn) = 0
  qlm_coordspinz(hn) = 0
  
  qlm_adm_energy(hn) = 0
  qlm_adm_momentum_x(hn) = 0
  qlm_adm_momentum_y(hn) = 0
  qlm_adm_momentum_z(hn) = 0
  qlm_adm_angular_momentum_x(hn) = 0
  qlm_adm_angular_momentum_y(hn) = 0
  qlm_adm_angular_momentum_z(hn) = 0

  qlm_w_energy(hn) = 0
  qlm_w_momentum_x(hn) = 0
  qlm_w_momentum_y(hn) = 0
  qlm_w_momentum_z(hn) = 0
  qlm_w_angular_momentum_x(hn) = 0
  qlm_w_angular_momentum_y(hn) = 0
  qlm_w_angular_momentum_z(hn) = 0

  ! Compute weights for spherical integration (see Driscoll and Healy).
  ! These are the correct weights in a Gauss-Legendre-Senc.
  ! Also compare qlm_area with AHFinderDirect's area, they are in much
  ! better agreement than with the previous method using
  ! delta_theta*delta_phi.
  
  ntheta_inner = qlm_ntheta(hn) - 2*qlm_nghoststheta(hn)
  nphi_inner   = qlm_nphi(hn) - 2*qlm_nghostsphi(hn)
  
  allocate (weights(1+qlm_nghoststheta(hn) : qlm_ntheta(hn)-qlm_nghoststheta(hn)))
  
  do i = 1+qlm_nghoststheta(hn), qlm_ntheta(hn)-qlm_nghoststheta(hn)
     theta = qlm_origin_theta(hn) + (i-1)*qlm_delta_theta(hn)
     sum1 = 0
     do l = 0, (ntheta_inner-1)/2
        sum1 = sum1 + sin((2*l+1)*theta)/(2*l+1)
     end do
     weights(i) = 8*pi * sum1 / (nphi_inner * ntheta_inner)
  end do
  
  do j = 1+qlm_nghostsphi(hn), qlm_nphi(hn)-qlm_nghostsphi(hn)
     do i = 1+qlm_nghoststheta(hn), qlm_ntheta(hn)-qlm_nghoststheta(hn)
        
        theta = qlm_origin_theta(hn) + (i-1)*qlm_delta_theta(hn)
        phi   = qlm_origin_phi(hn)   + (j-1)*qlm_delta_phi(hn)
        
        ! 2-metric on the horizon
        qq(1,1) = qlm_qtt(i,j,hn)
        qq(1,2) = qlm_qtp(i,j,hn)
        qq(2,2) = qlm_qpp(i,j,hn)
        qq(2,1) = qq(1,2)
        call calc_2det (qq, dtq)
        
        ll(0) = qlm_l0(i,j,hn)
        ll(1) = qlm_l1(i,j,hn)
        ll(2) = qlm_l2(i,j,hn)
        ll(3) = qlm_l3(i,j,hn)
        
        nn(0) = qlm_n0(i,j,hn)
        nn(1) = qlm_n1(i,j,hn)
        nn(2) = qlm_n2(i,j,hn)
        nn(3) = qlm_n3(i,j,hn)
        
        mm(0) = qlm_m0(i,j,hn)
        mm(1) = qlm_m1(i,j,hn)
        mm(2) = qlm_m2(i,j,hn)
        mm(3) = qlm_m3(i,j,hn)
        
        tt = (ll + nn) / sqrt(two)
        ss = (ll - nn) / sqrt(two)
        
        npalpha = qlm_npalpha(i,j,hn)
        npbeta  = qlm_npbeta(i,j,hn)
        
        gg(1,1) = qlm_gxx(i,j)
        gg(1,2) = qlm_gxy(i,j)
        gg(1,3) = qlm_gxz(i,j)
        gg(2,2) = qlm_gyy(i,j)
        gg(2,3) = qlm_gyz(i,j)
        gg(3,3) = qlm_gzz(i,j)
        gg(2,1) = gg(1,2)
        gg(3,1) = gg(1,3)
        gg(3,2) = gg(2,3)
        
        dgg(1,1,1) = qlm_dgxxx(i,j)
        dgg(1,2,1) = qlm_dgxyx(i,j)
        dgg(1,3,1) = qlm_dgxzx(i,j)
        dgg(2,2,1) = qlm_dgyyx(i,j)
        dgg(2,3,1) = qlm_dgyzx(i,j)
        dgg(3,3,1) = qlm_dgzzx(i,j)
        dgg(2,1,1) = dgg(1,2,1)
        dgg(3,1,1) = dgg(1,3,1)
        dgg(3,2,1) = dgg(2,3,1)
        dgg(1,1,2) = qlm_dgxxy(i,j)
        dgg(1,2,2) = qlm_dgxyy(i,j)
        dgg(1,3,2) = qlm_dgxzy(i,j)
        dgg(2,2,2) = qlm_dgyyy(i,j)
        dgg(2,3,2) = qlm_dgyzy(i,j)
        dgg(3,3,2) = qlm_dgzzy(i,j)
        dgg(2,1,2) = dgg(1,2,2)
        dgg(3,1,2) = dgg(1,3,2)
        dgg(3,2,2) = dgg(2,3,2)
        dgg(1,1,3) = qlm_dgxxz(i,j)
        dgg(1,2,3) = qlm_dgxyz(i,j)
        dgg(1,3,3) = qlm_dgxzz(i,j)
        dgg(2,2,3) = qlm_dgyyz(i,j)
        dgg(2,3,3) = qlm_dgyzz(i,j)
        dgg(3,3,3) = qlm_dgzzz(i,j)
        dgg(2,1,3) = dgg(1,2,3)
        dgg(3,1,3) = dgg(1,3,3)
        dgg(3,2,3) = dgg(2,3,3)
        
        kk(1,1) = qlm_kxx(i,j)
        kk(1,2) = qlm_kxy(i,j)
        kk(1,3) = qlm_kxz(i,j)
        kk(2,2) = qlm_kyy(i,j)
        kk(2,3) = qlm_kyz(i,j)
        kk(3,3) = qlm_kzz(i,j)
        kk(2,1) = kk(1,2)
        kk(3,1) = kk(1,3)
        kk(3,2) = kk(2,3)
        
        
        
        call calc_det (gg, dtg)
        call calc_inv (gg, dtg, gu)
        call calc_trace (gu, kk, trk)
        
        ! Calculate time derivatives
        call calc_3metricdot_simple (kk, gg_dot)
        
        ! Calculate 4-metric
        call calc_4metricderivs_simple (gg, dgg, gg_dot, g4,dg4)

        ! Calculate derivative of 4-metric inverse
        call calc_invderiv (gu, dgg, dgu)
        
        ! "Perturbative" 4-metric
        h4 = g4 - eta4
        dh4 = dg4
        
        xx(1) = qlm_x(i,j,hn)
        xx(2) = qlm_y(i,j,hn)
        xx(3) = qlm_z(i,j,hn)
        
        ee(1,1:2) = deriv (qlm_x(:,:,hn), i, j, delta_space)
        ee(2,1:2) = deriv (qlm_y(:,:,hn), i, j, delta_space)
        ee(3,1:2) = deriv (qlm_z(:,:,hn), i, j, delta_space)
        
        xi(1) = qlm_xi_t(i,j,hn)
        xi(2) = qlm_xi_p(i,j,hn)
        
        if (i == i_eq) then
           qlm_equatorial_circumference(hn) = qlm_equatorial_circumference(hn) &
                + sqrt(qq(2,2)) * qlm_delta_phi(hn)
        end if
        
        if (j == j_p0) then
           qlm_polar_circumference_0(hn) = qlm_polar_circumference_0(hn) &
                + sqrt(qq(1,1)) * qlm_delta_theta(hn)
        end if
        
        if (j == j_p2) then
           qlm_polar_circumference_pi_2(hn) = qlm_polar_circumference_pi_2(hn) &
                + sqrt(qq(1,1)) * qlm_delta_theta(hn)
        end if
        
        qlm_area(hn) = qlm_area(hn) &
             + sqrt(dtq) * weights(i)
        
        ! s^i: outward spacelike normal
        ! K_ij: extrinsic curvature
        ! phi^i: rotational Killing vector
        
        do a=1,3
           xi1(a) = 0
           do b=1,2
              xi1(a) = xi1(a) + ee(a,b) * xi(b)
           end do
        end do
        
        ! phi^i omega_i = - phi^i s^j K_ij
        
        spin = 0
        do a=1,3
           do b=1,3
              spin = spin + xi1(a) * ss(b) * kk(a,b)
           end do
        end do
        qlm_spin_density(i,j) = spin
        
        ! phi^i omega_i = (alpha + ~beta) phi^i m_i + complex conjugate
        npspin = 0
        do a=1,3
           do b=1,3
              npspin = npspin + 2 * real((npalpha + conjg(npbeta)) * xi1(a) * gg(a,b) * mm(b))
           end do
        end do
        
        ! phi^i omega_i = 1/2 f Im Psi_2
        ! (or is it   phi^i omega_i = - f Im Psi_2   ?)
        wsspin = - qlm_inv_z(i,j,hn) * aimag(qlm_psi2(i,j,hn))
        
        ! x = sin theta cos phi
        ! y = sin theta sin phi
        ! z = cos theta
        ! xi_x = (0,-z,y)
        ! xi_y = (z,0,-x)
        ! xi_z = (-y,x,0)
        tx = qlm_x(i,j,hn) - qlm_origin_x(hn)
        ty = qlm_y(i,j,hn) - qlm_origin_y(hn)
        tz = qlm_z(i,j,hn) - qlm_origin_z(hn)
        xi1_x(:) = (/ zero, -tz, ty /)
        xi1_y(:) = (/ tz, zero, -tx /)
        xi1_z(:) = (/ -ty, tx, zero /)
        coordspinx = 0
        coordspiny = 0
        coordspinz = 0
        do a=1,3
           do b=1,3
              coordspinx = coordspinx + xi1_x(a) * ss(b) * kk(a,b)
              coordspiny = coordspiny + xi1_y(a) * ss(b) * kk(a,b)
              coordspinz = coordspinz + xi1_z(a) * ss(b) * kk(a,b)
           end do
        end do
        
        xi1_x(:) = (/ 1, 0, 0 /)
        xi1_y(:) = (/ 0, 1, 0 /)
        xi1_z(:) = (/ 0, 0, 1 /)
        
        qlm_spin(hn) = qlm_spin(hn) &
             + spin * sqrt(dtq) * weights(i)

        qlm_npspin(hn) = qlm_npspin(hn) &
             + npspin * sqrt(dtq) * weights(i)
        
        qlm_wsspin(hn) = qlm_wsspin(hn) &
             + wsspin * sqrt(dtq) * weights(i)
        
        qlm_coordspinx(hn) = qlm_coordspinx(hn) &
             + coordspinx * sqrt(dtq) * weights(i)
        qlm_coordspiny(hn) = qlm_coordspiny(hn) &
             + coordspiny * sqrt(dtq) * weights(i)
        qlm_coordspinz(hn) = qlm_coordspinz(hn) &
             + coordspinz * sqrt(dtq) * weights(i)
        
        ! ADM quantities
        ! Alcubierre, Appendix A, p. 402 ff, eqns. (A.5) - (A.7):
        
        ! ADM energy
        ! E_adm = (1/16 pi) int_S(r) [delta^jk h_ik,j - (trace(h)),i] n^i r^2 dOmega
        adm_energy = 0
        do a=1,3
           do b=1,3
              do c=1,3
                 adm_energy = adm_energy + (delta3(a,b) * dh4(a,c,b) &
             & - gu(a,b) * dh4(a,b,c) - h4(a,b) * dgu(a,b,c) ) * ss(c)
              end do
           end do
        end do
        qlm_adm_energy(hn) = qlm_adm_energy(hn) &
             & + adm_energy / (16*pi) &
             &   * sqrt(dtq) * qlm_delta_theta(hn) * qlm_delta_phi(hn)
        
        ! ADM momentum
        ! P_adm^i = (1/8 pi) int_S(r) [K^i_j - K delta^i_j] n^j r^2 dOmega
        do a=1,3
           adm_mom(a) = 0
           do b=1,3
              do c=1,3
                 adm_mom(a) = adm_mom(a) &
             &        + gu(a,c) * kk(c,b) * ss(b) 
              end do
              adm_mom(a) = adm_mom(a) &
                   + delta3(a,b) * trk * ss(b) 
           end do
        end do
        qlm_adm_momentum_x(hn) = qlm_adm_momentum_x(hn) &
             & + adm_mom(1) / (8*pi) &
             &   * sqrt(dtq) * qlm_delta_theta(hn) * qlm_delta_phi(hn)
        qlm_adm_momentum_y(hn) = qlm_adm_momentum_y(hn) &
             & + adm_mom(2) / (8*pi) &
             &   * sqrt(dtq) * qlm_delta_theta(hn) * qlm_delta_phi(hn)
        qlm_adm_momentum_z(hn) = qlm_adm_momentum_z(hn) &
             & + adm_mom(3) / (8*pi) &
             &   * sqrt(dtq) * qlm_delta_theta(hn) * qlm_delta_phi(hn)
        
        ! ADM angular momentum
        ! J_adm^i = (1/8 pi) int_S(r)
        !   epsilon^ijk x_j (K_kl - delta_kl K) * n^l dS 
        do a=1,3
           adm_amom(a) = 0
           do b=1,3
              do c=1,3
                 do l=1,3
                    do m=1,3
                       adm_amom(a) = adm_amom(a)  &
                &         + epsilon3(a,b,c) * xx(b) * gu(c,m) * kk(m,l) * ss(l)
                    end do
                    adm_amom(a) = adm_amom(a)  &
                &    - epsilon3(a,b,c) * xx(b) *  delta3(c,l) * trk * ss(l)
                 end do
              end do
           end do
        end do
        qlm_adm_angular_momentum_x(hn) = qlm_adm_angular_momentum_x(hn) &
             & + adm_amom(1) / (8*pi) &
             &   * sqrt(dtq) * qlm_delta_theta(hn) * qlm_delta_phi(hn)
        qlm_adm_angular_momentum_y(hn) = qlm_adm_angular_momentum_y(hn) &
             & + adm_amom(2) / (8*pi) &
             &   * sqrt(dtq) * qlm_delta_theta(hn) * qlm_delta_phi(hn)
        qlm_adm_angular_momentum_z(hn) = qlm_adm_angular_momentum_z(hn) &
             & + adm_amom(3) / (8*pi) &
             &   * sqrt(dtq) * qlm_delta_theta(hn) * qlm_delta_phi(hn)

        ! Weinberg pseudotensor quantities
        ! Weinberg, chapter 7.6, pp. 165 ff, eqns. (7.6.22) - (7.6.24):
        
        ! Weinberg energy
        ! E_w = (-1/16 pi) int_S(r) [h_jj,i - hij,j] n_i r^2 dOmega
        w_energy = 0
        do a=1,3
           do b=1,3
              w_energy = w_energy + (dh4(b,b,a) - dh4(a,b,b)) * ss(a)
           end do
        end do
        qlm_w_energy(hn) = qlm_w_energy(hn) &
             & + w_energy / (-16*pi) &
             &   * sqrt(dtq) * weights(i)
        
        ! Weinberg momentum
        ! P_ll^j = (-1/16 pi) int_S(r)
        !    [- h_kk,0 delta_ij + h_k0,k delta_ij
        !     - h_j0,i + h_ij,0                  ] n_i r^2 dOmega
        ! 
        ! P_ll^j = (1/8 pi) int_S(r) [K_ij - K delta_ij] n_i r^2 dOmega
        do a=1,3
           w_mom(a) = 0
           do b=1,3
              w_mom(a) = w_mom(a) &
                   + (- dh4(b,b,0) + dh4(b,0,b)) * ss(a) &
                   + (- dh4(a,0,b) + dh4(b,a,0)) * ss(b)
           end do
        end do
        qlm_w_momentum_x(hn) = qlm_w_momentum_x(hn) &
             & + w_mom(1) / (-16*pi) &
             &   * sqrt(dtq) * weights(i)
        qlm_w_momentum_y(hn) = qlm_w_momentum_y(hn) &
             & + w_mom(2) / (-16*pi) &
             &   * sqrt(dtq) * weights(i)
        qlm_w_momentum_z(hn) = qlm_w_momentum_z(hn) &
             & + w_mom(3) / (-16*pi) &
             &   * sqrt(dtq) * weights(i)
        
        ! Weinberg angular momentum
        ! J_ll^jk = (-1/16 pi) int_S(r)
        !    [- x_j h_0k,i + x_k h_0j,i
        !     + x_j h_ki,0 - x_k h_ji,0
        !     + h_0k delta_ij - h_0j delta_ik] n_i r^2 dOmega
        do a=1,3
           do b=1,3
              w_amom(a,b) = 0
              do c=1,3
                 w_amom(a,b) = w_amom(a,b) + ss(c) * ( &
                      + (- xx(a) * dh4(0,b,c) + xx(b) * dh4(0,a,c)) &
                      + (+ xx(a) * dh4(b,c,0) - xx(b) * dh4(a,c,0)) &
                      + (h4(0,b) * delta4(c,a) - h4(0,a) * delta4(c,b)))
              end do
           end do
        end do
        qlm_w_angular_momentum_x(hn) = qlm_w_angular_momentum_x(hn) &
             & + w_amom(2,3) / (-16*pi) &
             &   * sqrt(dtq) * weights(i)
        qlm_w_angular_momentum_y(hn) = qlm_w_angular_momentum_y(hn) &
             & + w_amom(3,1) / (-16*pi) &
             &   * sqrt(dtq) * weights(i)
        qlm_w_angular_momentum_z(hn) = qlm_w_angular_momentum_z(hn) &
             & + w_amom(1,2) / (-16*pi) &
             &   * sqrt(dtq) * weights(i)
        
     end do
  end do
  
  deallocate(weights)
  
  qlm_polar_circumference_0(hn) = qlm_polar_circumference_0(hn) * 2
  qlm_polar_circumference_pi_2(hn) = qlm_polar_circumference_pi_2(hn) * 2
  
  ! A = 4 pi R^2
  ! R = 2 M
  qlm_radius(hn) = sqrt(qlm_area(hn) / (4*pi))
  qlm_irreducible_mass(hn) = qlm_radius(hn) / 2
  
  qlm_spin(hn) = qlm_spin(hn) / (8*pi)
  qlm_mass(hn) = 1/(2*qlm_radius(hn)) * sqrt(qlm_radius(hn)**4 + 4*qlm_spin(hn)**2)
  qlm_cvspin(hn) = qlm_cvspin(hn) / (8*pi)
  qlm_npspin(hn) = qlm_npspin(hn) / (-8*pi)
  qlm_wsspin(hn) = qlm_wsspin(hn) / (-4*pi)
  qlm_coordspinx(hn) = qlm_coordspinx(hn) / (8*pi)
  qlm_coordspiny(hn) = qlm_coordspiny(hn) / (8*pi)
  qlm_coordspinz(hn) = qlm_coordspinz(hn) / (8*pi)
  
  ! The event horizon is at r = M + sqrt (M^2 - a^2)
  ! with x^2 + y^2 + z^2 = rho^2 = r^2 + a^2 (1 - z^2 / r^2)
  
  call guess_mass_spin &
       (qlm_area(hn), qlm_equatorial_circumference(hn), qlm_mass_guess(hn), qlm_spin_guess(hn))
  
  
  
  if (verbose/=0) then
     
     write (msg, '("Geometric quantities for surface ",i4,":")') hn-1
     call CCTK_INFO (msg)
     write (msg, '("   Area A:                       ",g16.6)') qlm_area(hn)
     call CCTK_INFO (msg)
     write (msg, '("   Irreducible mass M = R/2:     ",g16.6)') qlm_irreducible_mass(hn)
     call CCTK_INFO (msg)
     write (msg, '("   Areal radius R = sqrt(A/4pi): ",g16.6)') qlm_radius(hn)
     call CCTK_INFO (msg)
     
     write (msg, '("Coordinate-dependent quantities for surface ",i4,":")') hn-1
     call CCTK_INFO (msg)
     write (msg, '("   Equatorial circumference:        ",g16.6)') &
          qlm_equatorial_circumference(hn)
     call CCTK_INFO (msg)
     write (msg, '("   Polar circumference at phi=0:    ",g16.6)') &
          qlm_polar_circumference_0(hn)
     call CCTK_INFO (msg)
     write (msg, '("   Polar circumference at phi=pi/2: ",g16.6)') &
          qlm_polar_circumference_pi_2(hn)
     call CCTK_INFO (msg)
     if (qlm_spin_guess(hn) >= 0) then
        write (msg, '("   Spin guess J from distortion:    ",g16.6)') &
             qlm_spin_guess(hn)
        call CCTK_INFO (msg)
     else
        call CCTK_INFO ("   No valid spin guess from distortion.")
        call CCTK_INFO ("   (Spin guess is imaginary.)")
        write (msg, '("   Magnitude of invalid spin guess: ",g16.6)') &
             abs(qlm_spin_guess(hn))
        call CCTK_INFO (msg)
     end if
     write (msg, '("   Mass guess M from distortion:    ",g16.6)') &
          qlm_mass_guess(hn)
     call CCTK_INFO (msg)
     
     write (msg, '("Isolated Horizon quantities for surface ",i4,":")') hn-1
     call CCTK_INFO (msg)
     ev = cmplx(qlm_killing_eigenvalue_re(hn), qlm_killing_eigenvalue_im(hn),rk)
     write (msg, '("   Killing vector field eigenvalue norm:   ",g14.6)') abs(ev)
     call CCTK_INFO (msg)
     write (msg, '("   Spin J:                                 ",g14.6)') qlm_spin(hn)
     call CCTK_INFO (msg)
     write (msg, '("   Kerr spin parameter a = J/M:            ",g14.6)') qlm_spin(hn) / qlm_mass(hn)
     call CCTK_INFO (msg)
     write (msg, '("   Dimensionless spin parameter a = J/M^2: ",g14.6)') qlm_spin(hn) / qlm_mass(hn)**2
     call CCTK_INFO (msg)
     write (msg, '("   Spin J from NP:                         ",g14.6)') qlm_npspin(hn)
     call CCTK_INFO (msg)
     write (msg, '("   Spin J from phi-coordinate-vector:      ",g14.6)') qlm_cvspin(hn)
     call CCTK_INFO (msg)
     write (msg, '("   Mass M:                                 ",g14.6)') qlm_mass(hn)
     call CCTK_INFO (msg)
     
     write (msg, '("Global quantities for surface ",i4,":")') hn-1
     call CCTK_INFO (msg)
     write (msg, '("   ADM energy:             ",g16.6)') qlm_adm_energy(hn)
     call CCTK_INFO (msg)
     write (msg, '("   ADM momentum x:         ",g16.6)') qlm_adm_momentum_x(hn)
     call CCTK_INFO (msg)
     write (msg, '("   ADM momentum y:         ",g16.6)') qlm_adm_momentum_y(hn)
     call CCTK_INFO (msg)
     write (msg, '("   ADM momentum z:         ",g16.6)') qlm_adm_momentum_z(hn)
     call CCTK_INFO (msg)
     write (msg, '("   ADM angular momentum x: ",g16.6)') qlm_adm_angular_momentum_x(hn)
     call CCTK_INFO (msg)
     write (msg, '("   ADM angular momentum y: ",g16.6)') qlm_adm_angular_momentum_y(hn)
     call CCTK_INFO (msg)
     write (msg, '("   ADM angular momentum z: ",g16.6)') qlm_adm_angular_momentum_z(hn)
     call CCTK_INFO (msg)
     write (msg, '("   Weinberg energy:             ",g16.6)') qlm_w_energy(hn)
     call CCTK_INFO (msg)
     write (msg, '("   Weinberg momentum x:         ",g16.6)') qlm_w_momentum_x(hn)
     call CCTK_INFO (msg)
     write (msg, '("   Weinberg momentum y:         ",g16.6)') qlm_w_momentum_y(hn)
     call CCTK_INFO (msg)
     write (msg, '("   Weinberg momentum z:         ",g16.6)') qlm_w_momentum_z(hn)
     call CCTK_INFO (msg)
     write (msg, '("   Weinberg angular momentum x: ",g16.6)') qlm_w_angular_momentum_x(hn)
     call CCTK_INFO (msg)
     write (msg, '("   Weinberg angular momentum y: ",g16.6)') qlm_w_angular_momentum_y(hn)
     call CCTK_INFO (msg)
     write (msg, '("   Weinberg angular momentum z: ",g16.6)') qlm_w_angular_momentum_z(hn)
     call CCTK_INFO (msg)
     
  end if
  
contains
  
  subroutine guess_mass_spin (area, circumference, mass, spin)
    CCTK_REAL, intent(in)  :: area, circumference
    CCTK_REAL, intent(out) :: mass, spin
    CCTK_REAL :: radius, amom2, amom
    
    ! equatorial circumference L, area A
    
    ! L = 2 pi (r^2 + a^2) / r
    ! A = 4 pi (r^2 + a^2)
    ! r = M + sqrt (M^2 - a^2)
    
    ! r = A / (2 L)
    ! a^2 = A / (4 pi) - r^2   ("spin" a = J/M = specific angular momentum)
    ! M = (r^2 + a^2) / (2 r)
    
    ! J = a M   (angular momentum)
    
    radius = area / (2*circumference)
    amom2 = area / (4*pi) - radius**2
    amom = sign(sqrt(abs(amom2)), amom2)
    mass = (radius**2 + amom2) / (2*radius)
    spin = amom * mass
    
    ! equatorial circumference L, area A
    ! mass M, areal radius R, spin a = J/M^2
    
    ! A = 4 pi R^2
    ! L = 4 pi M
    
    ! a^2 = (R/M)^2 - 1/4 (R/M)^4
    ! J = a M^2
    
  end subroutine guess_mass_spin
  
end subroutine qlm_analyse
