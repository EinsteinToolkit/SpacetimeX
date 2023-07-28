#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



module qlm_killing_transportation
  use cctk
  use constants
  use qlm_variables
  use ricci2
  use tensor2
  implicit none
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  private
  public transport_along_equator
  public transport_along_meridians
  
contains
  
  subroutine transport_along_equator (CCTK_ARGUMENTS, hn, i0, xi, chi)
    DECLARE_CCTK_ARGUMENTS
    integer,   intent(in)    :: hn
    integer,   intent(in)    :: i0
    CCTK_REAL, intent(inout) :: xi(2), chi
    integer :: j0
    integer :: nsteps
    
    j0 = 1+qlm_nghostsphi(hn)
    nsteps = qlm_nphi(hn) - 2*qlm_nghostsphi(hn)
    
    call transport (CCTK_PASS_FTOF, hn, i0, j0, 0, 1, nsteps, xi, chi)
    
  end subroutine transport_along_equator
  
  
  
  subroutine transport_along_meridians (CCTK_ARGUMENTS, hn, i0)
    DECLARE_CCTK_ARGUMENTS
    integer, intent(in) :: hn
    integer, intent(in) :: i0
    CCTK_REAL :: xi(2), chi
    integer   :: j0
    integer   :: dir
    integer   :: nsteps
    
    do j0 = 1+qlm_nghostsphi(hn), qlm_nphi(hn)-qlm_nghostsphi(hn)
       
       do dir=-1,+1,2
          
          if (dir==-1) nsteps = i0 - (1+qlm_nghoststheta(hn))
          if (dir==+1) nsteps = (qlm_ntheta(hn)-qlm_nghoststheta(hn)) - i0
          
          xi(1) = qlm_xi_t(i0,j0,hn)
          xi(2) = qlm_xi_p(i0,j0,hn)
          chi = qlm_chi(i0,j0,hn)
          
          call transport (CCTK_PASS_FTOF, hn, i0, j0, dir, 0, nsteps, xi, chi)
          
       end do
       
    end do
    
  end subroutine transport_along_meridians
  
  
  
  subroutine transport (CCTK_ARGUMENTS, hn, i0, j0, di, dj, nsteps, xi, chi)
    DECLARE_CCTK_ARGUMENTS
    integer,  intent(in)    :: hn
    integer,  intent(in)    :: i0, j0
    integer,  intent(in)    :: di, dj
    integer,  intent(in)    :: nsteps
    CCTK_REAL,intent(inout) :: xi(2), chi
    CCTK_REAL :: vv(2)
    CCTK_REAL :: xi_dot(2), chi_dot
    CCTK_REAL :: xi1(2), chi1
    CCTK_REAL :: xi1_dot(2), chi1_dot
    integer   :: n
    integer   :: i, j
    
    i = i0
    j = j0
    
    vv(1) = di * qlm_delta_theta(hn)
    vv(2) = dj * qlm_delta_phi(hn)
    
    do n=1,nsteps
       
       call transport_rhs (CCTK_PASS_FTOF, hn, i, j, xi, chi, vv, xi_dot, chi_dot)
       xi1 = xi + xi_dot
       chi1 = chi + chi_dot
       i = i + di
       j = j + dj
       if (j <           1+qlm_nghostsphi(hn)) j = j + (qlm_nphi(hn)-2*qlm_nghostsphi(hn))
       if (j > qlm_nphi(hn)-qlm_nghostsphi(hn)) j = j - (qlm_nphi(hn)-2*qlm_nghostsphi(hn))
       
       call transport_rhs &
            (CCTK_PASS_FTOF, hn, i, j, xi1, chi1, vv, xi1_dot, chi1_dot)
       xi = xi + 0.5d0 * (xi_dot + xi1_dot)
       chi = chi + 0.5d0 * (chi_dot + chi1_dot)
       
       qlm_xi_t(i,j,hn) = xi(1)
       qlm_xi_p(i,j,hn) = xi(2)
       qlm_chi(i,j,hn) = chi

    end do
    
  end subroutine transport
  
  
  
  subroutine transport_rhs (CCTK_ARGUMENTS, hn, i, j, xi, chi, vv, xi_dot, chi_dot)
    DECLARE_CCTK_ARGUMENTS
    integer,   intent(in)  :: hn
    integer,   intent(in)  :: i, j
    CCTK_REAL, intent(in)  :: xi(2), chi
    CCTK_REAL, intent(in)  :: vv(2)
    CCTK_REAL, intent(out) :: xi_dot(2), chi_dot
    CCTK_REAL :: qq(2,2), dqq(2,2,2), dtq, qu(2,2), gamma(2,2,2), rsc
    
    if (i<1+qlm_nghoststheta(hn) .or. i>qlm_ntheta(hn)-qlm_nghoststheta(hn) &
         .or. j<1+qlm_nghostsphi(hn) .or. j>qlm_nphi(hn)-qlm_nghostsphi(hn)) then
       call CCTK_WARN (0, "internal error")
    end if
    if (i-1<1 .or. i+1>qlm_ntheta(hn) .or. j-1<1 .or. j+1>qlm_nphi(hn)) then
       call CCTK_WARN (0, "internal error")
    end if
    
    qq(1,1) = qlm_qtt(i,j,hn)
    qq(1,2) = qlm_qtp(i,j,hn)
    qq(2,2) = qlm_qpp(i,j,hn)
    qq(2,1) = qq(1,2)
    
    dqq(1,1,1) = qlm_dqttt(i,j)
    dqq(1,2,1) = qlm_dqtpt(i,j)
    dqq(2,2,1) = qlm_dqppt(i,j)
    dqq(1,1,2) = qlm_dqttp(i,j)
    dqq(1,2,2) = qlm_dqtpp(i,j)
    dqq(2,2,2) = qlm_dqppp(i,j)
    dqq(2,1,:) = dqq(1,2,:)
    
    rsc = qlm_rsc(i,j,hn)
    
    call calc_2det (qq, dtq)
    call calc_2inv (qq, dtq, qu)
    
    call calc_2connections (qu, dqq, gamma)
    
    call killing_transport_rhs &
         (xi, chi, qq, dtq, qu, gamma, rsc, vv, xi_dot, chi_dot)
    
  end subroutine transport_rhs
  
  
  
  subroutine killing_transport_rhs &
       (xi, chi, qq, dtq, qu, gamma2, rsc2, vv, xi_dot, chi_dot)
    CCTK_REAL, intent(in)  :: xi(2), chi
    CCTK_REAL, intent(in)  :: qq(2,2), dtq, qu(2,2), gamma2(2,2,2), rsc2
    CCTK_REAL, intent(in)  :: vv(2)
    CCTK_REAL, intent(out) :: xi_dot(2), chi_dot
    integer :: i, k, l
    
    ! Wald eqn (C.3.6):
    ! D_k D_j xi^i = R^i_jkl xi^l
    
    ! define:
    ! L^i_j = D_j x^i
    
    ! then: see Wald eqns (C.3.7) and (C.3.8):
    ! v^k D_k xi^i = L^i_k v^k
    ! v^k D_k L^i_j = R^i_jkl v^k xi^l
    
    ! in 2D we have:
    ! R_ijkl = 1/2 q R epsilon2_ij epsilon2_kl
    ! L_ij = epsilon2_ij sqrt(q) chi
    
    ! then:
    ! v^k D_k xi^i = epsilon2^i_k sqrt(q) chi v^k
    ! v^k D_k epsilon2^i_j sqrt(q) chi = R^i_jkl v^k xi^l
    ! v^k D_k chi = 1/2 sqrt(q) R epsilon2_kl v^k xi^l
    
    ! define:
    ! X_dot = v^i d/dx^i X   (partial derivatives)
    
    do i=1,2
       xi_dot(i) = 0
       do k=1,2
          do l=1,2
             xi_dot(i) = xi_dot(i) &
                  + qu(i,l) * epsilon2(l,k) * sqrt(dtq) * chi * vv(k) &
                  - vv(k) * gamma2(i,l,k) * xi(l)
          end do
       end do
    end do
    
    chi_dot = 0
    do k=1,2
       do l=1,2
          chi_dot = chi_dot &
               + 0.5d0 * sqrt(dtq) * rsc2 * epsilon2(k,l) * vv(k) * xi(l)
       end do
    end do
    
  end subroutine killing_transport_rhs
  
  
  
#if 0
  subroutine killing_equation (qq, qu, xi, gxi, zeta, trzeta, zetasq)
    CCTK_REAL, intent(in)  :: qq(2,2), qu(2,2)
    CCTK_REAL, intent(in)  :: xi(2), gxi(2,2)
    CCTK_REAL, intent(out) :: zeta(2,2), trzeta, zetasq
    integer :: i, j, k, l
    
    do i=1,2
       do j=1,2
          ! zeta_ij = D_i xi_j + D_j xi_i
          zeta(i,j) = 0
          do k=1,2
             do l=1,2
                zeta(i,j) = zeta(i,j) + qq(j,k) * gxi(k,i) + qq(i,k) * gxi(k,j)
             end do
          end do
       end do
    end do
    
    trzeta = 0
    do i=1,2
       do j=1,2
          trzeta = trzeta + qu(i,j) * zeta(i,j)
       end do
    end do
    
    zetasq = 0
    do i=1,2
       do j=1,2
          do k=1,2
             do l=1,2
                zetasq = zetasq + qu(i,k) * qu(j,l) * zeta(i,j) * zeta(k,l)
             end do
          end do
       end do
    end do
  end subroutine killing_equation
#endif
  
end module qlm_killing_transportation
