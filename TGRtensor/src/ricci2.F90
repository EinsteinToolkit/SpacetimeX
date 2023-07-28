! $Header$

#include "cctk.h"

module ricci2
  implicit none
  private
  public calc_2connections
  public calc_2connectionderivs
  public calc_2ricci
  
contains
  
  subroutine calc_2connections (gu, dgg, gamma)
    CCTK_REAL, intent(in)  :: gu(2,2), dgg(2,2,2)
    CCTK_REAL, intent(out) :: gamma(2,2,2)
    integer :: i,j,k,l
    ! Gamma^i_jk = 1/2 g^il (g_lj,k + g_lk,j - g_jk,l)
    do i=1,2
       do j=1,2
          do k=1,2
             gamma(i,j,k) = 0
             do l=1,2
                gamma(i,j,k) = gamma(i,j,k) + 0.5d0 * gu(i,l) &
                     * (dgg(l,j,k) + dgg(l,k,j) - dgg(j,k,l))
             end do
          end do
       end do
    end do
  end subroutine calc_2connections
  
  subroutine calc_2connectionderivs (gu, dgg, dgu, ddgg, dgamma)
    CCTK_REAL, intent(in)  :: gu(2,2), dgg(2,2,2), dgu(2,2,2), ddgg(2,2,2,2)
    CCTK_REAL, intent(out) :: dgamma(2,2,2,2)
    integer   :: i,j,k,l,m
    ! Gamma^i_jk,l = 1/2 g^im,l (g_mj,k + g_mk,j - g_jk,m)
    !                + 1/2 g^im (g_mj,kl + g_mk,jl - g_jk,ml)
    do i=1,2
       do j=1,2
          do k=1,2
             do l=1,2
                dgamma(i,j,k,l) = 0
                do m=1,2
                   dgamma(i,j,k,l) = dgamma(i,j,k,l) &
                        + 0.5d0 * dgu(i,m,l) * (dgg(m,j,k) + dgg(m,k,j) - dgg(j,k,m)) &
                        + 0.5d0 * gu(i,m) * (ddgg(m,j,k,l) + ddgg(m,k,j,l) - ddgg(j,k,m,l))
                end do
             end do
          end do
       end do
    end do
  end subroutine calc_2connectionderivs
  
  subroutine calc_2ricci (gamma, dgamma, ri)
    CCTK_REAL, intent(in)  :: gamma(2,2,2), dgamma(2,2,2,2)
    CCTK_REAL, intent(out) :: ri(2,2)
    integer :: i,j,k,l
    ! R_ij =   Gamma^k_ij,k - Gamma^k_ik,j
    !        + Gamma^k_lk Gamma^l_ij - Gamma^k_lj Gamma^l_ki
    do i=1,2
       do j=1,2
          ri(i,j) = 0
          do k=1,2
             ri(i,j) = ri(i,j) + dgamma(k,i,j,k) - dgamma(k,i,k,j)
             do l=1,2
                ri(i,j) = ri(i,j) + gamma(k,l,k) * gamma(l,i,j) &
                     &            - gamma(k,l,j) * gamma(l,k,i)
             end do
          end do
       end do
    end do
  end subroutine calc_2ricci
  
end module ricci2
