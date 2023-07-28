#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_killing_gradient (CCTK_ARGUMENTS, hn)
  use cctk
  use constants
  use qlm_boundary
  use qlm_derivs
  use tensor2
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_REAL, parameter :: two=2, half=1/two
  CCTK_REAL :: qq(2,2), dqq(2,2,2), dtq, qu(2,2), dqu(2,2,2)
  CCTK_REAL :: dpsi2(2), ddpsi2(2,2), ndpsi2, dndpsi2(2)
  CCTK_REAL :: xi(2), dxi(2,2), chi
  integer   :: i, j
  integer   :: a, b
  CCTK_REAL    :: delta_space(2)
  
  delta_space(:) = (/ qlm_delta_theta(hn), qlm_delta_phi(hn) /)
  
  ! Calculate the gradient of a scalar
  do j = 1+qlm_nghostsphi(hn), qlm_nphi(hn)-qlm_nghostsphi(hn)
     do i = 1+qlm_nghoststheta(hn), qlm_ntheta(hn)-qlm_nghoststheta(hn)
        
        ! 2-metric on the horizon
        qq(1,1) = qlm_qtt(i,j,hn)
        qq(1,2) = qlm_qtp(i,j,hn)
        qq(2,2) = qlm_qpp(i,j,hn)
        qq(2,1) = qq(1,2)
        
#if 0
        dqq(1,1,1) = qlm_dqttt(i,j)
        dqq(1,1,2) = qlm_dqttp(i,j)
        dqq(1,2,1) = qlm_dqtpt(i,j)
        dqq(1,2,2) = qlm_dqtpp(i,j)
        dqq(2,2,1) = qlm_dqppt(i,j)
        dqq(2,2,2) = qlm_dqppp(i,j)
        dqq(2,1,:) = dqq(1,2,:)
#endif
        
        call calc_2det (qq, dtq)
#if 0
        call calc_2inv (qq, dtq, qu)
        call calc_2invderiv (qu, dqq, dqu)
#endif
        
        dpsi2(1) = (abs2(qlm_psi2(i+1,j,hn)) - abs2(qlm_psi2(i-1,j,hn))) / (2*qlm_delta_theta(hn))
        dpsi2(2) = (abs2(qlm_psi2(i,j+1,hn)) - abs2(qlm_psi2(i,j-1,hn))) / (2*qlm_delta_phi(hn))
        
#if 0
        ddpsi2(1,1) = (abs2(qlm_psi2(i+1,j,hn)) - 2*abs2(qlm_psi2(i,j,hn)) + abs2(qlm_psi2(i-1,j,hn))) / qlm_delta_theta(hn)**2
        ddpsi2(2,2) = (abs2(qlm_psi2(i,j+1,hn)) - 2*abs2(qlm_psi2(i,j,hn)) + abs2(qlm_psi2(i,j-1,hn))) / qlm_delta_phi(hn)**2
        ddpsi2(1,1) = (abs2(qlm_psi2(i-1,j-1,hn)) - abs2(qlm_psi2(i+1,j-1,hn)) - abs2(qlm_psi2(i-1,j+1,hn)) + abs2(qlm_psi2(i+1,j+1,hn))) / (4*qlm_delta_theta(hn)*qlm_delta_phi(hn))
        ddpsi2(2,1) = ddpsi2(1,2)
        
        ! ndpsi2 = ||grad |Psi_2|^2||
        ndpsi2 = 0
        do a=1,2
           do b=1,2
              ndpsi2 = ndpsi2 + qu(a,b) * dpsi2(a) * dpsi2(b)
           end do
        end do
        ndpsi2 = sqrt(ndpsi2)
        
        ! dndpsi2 = grad ||grad |Psi_2|^2||
        do a=1,2
           dndpsi2(a) = 0
           do b=1,2
              do c=1,2
                 dndpsi2(a) = dndpsi2(a) + 1 / (2*ndpsi2) * (qu(b,c) * ddpsi2(b,a) * dpsi2(c) + qu(b,c) * dpsi2(b) * ddpsi2(c,a) + dqu(b,c,a) * dpsi2(b) * dpsi2(c))
              end do
           end do
        end do
#endif
        
        ! xi^a = eps^ab D_b |Psi_2|^2
        do a=1,2
           xi(a) = 0
           do b=1,2
              xi(a) = xi(a) + sqrt(dtq) * epsilon2(a,b) * dpsi2(b)
           end do
        end do
        
        qlm_xi_t(i,j,hn) = xi(1)
        qlm_xi_p(i,j,hn) = xi(2)
        
#if 0
        ! xi^a = eps^ab D_b ||D_c |Psi_2|^2||
        do a=1,2
           xi(a) = 0
           do b=1,2
              xi(a) = xi(a) + sqrt(dtq) * epsilon2(a,b) * dndpsi2(b)
           end do
        end do
        
        qlm_xi_t(i,j,hn) = xi(1)
        qlm_xi_p(i,j,hn) = xi(2)
#endif
        
     end do
  end do
  
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_xi_t(:,:,hn), -1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_xi_p(:,:,hn), -1)
  
  
  
  ! fix up xi (which must not be zero)
  do j = 1+qlm_nghostsphi(hn), qlm_nphi(hn)-qlm_nghostsphi(hn)
     do i = 1+qlm_nghoststheta(hn), qlm_ntheta(hn)-qlm_nghoststheta(hn)
        
        xi(1) = qlm_xi_t(i,j,hn)
        xi(2) = qlm_xi_p(i,j,hn)
        
        if (sum(xi**2) < 1.0d-4**2) then
           
           qlm_xi_t(i,j,hn) = sum(qlm_xi_t(i:i+1,j:j+1,hn)) / 4
           qlm_xi_p(i,j,hn) = sum(qlm_xi_p(i:i+1,j:j+1,hn)) / 4
           
        end if
        
     end do
  end do
  
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_xi_t(:,:,hn), -1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_xi_p(:,:,hn), -1)
  
  
  
  ! set up the derivative of xi (which is not really needed)
  do j = 1+qlm_nghostsphi(hn), qlm_nphi(hn)-qlm_nghostsphi(hn)
     do i = 1+qlm_nghoststheta(hn), qlm_ntheta(hn)-qlm_nghoststheta(hn)
        
        ! 2-metric on the horizon
        qq(1,1) = qlm_qtt(i,j,hn)
        qq(1,2) = qlm_qtp(i,j,hn)
        qq(2,2) = qlm_qpp(i,j,hn)
        qq(2,1) = qq(1,2)
        call calc_2det (qq, dtq)
        
        dxi(1,1:2) = deriv (qlm_xi_t(:,:,hn), i, j, delta_space)
        dxi(2,1:2) = deriv (qlm_xi_p(:,:,hn), i, j, delta_space)
        
        ! eps_ab sqrt(q) chi = D_b xi_a
        !        sqrt(q) chi = -1/2 eps^ab D_a xi_b
        chi = 0
        do a=1,2
           do b=1,2
              chi = chi - half * sqrt(dtq) * epsilon2(a,b) * dxi(b,a)
           end do
        end do
        
        qlm_chi(i,j,hn) = chi
        
     end do
  end do
  
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_chi (:,:,hn), +1)
  
end subroutine qlm_killing_gradient
