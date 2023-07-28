! $Header$

#include "cctk.h"

#undef DEBUG

module ricci
  implicit none
  private
  public calc_connections
  public calc_connectionderivs
  public calc_connectionderivs2
  public calc_ricci
  public calc_ricciderivs
  
contains
  
  subroutine calc_connections (gu, dgg, gamma)
    CCTK_REAL, intent(in)  :: gu(3,3), dgg(3,3,3)
    CCTK_REAL, intent(out) :: gamma(3,3,3)
    integer :: i,j,k,l
    ! Gamma^i_jk = 1/2 g^il (g_lj,k + g_lk,j - g_jk,l)
    do i=1,3
       do j=1,3
          do k=1,3
             gamma(i,j,k) = 0
             do l=1,3
                gamma(i,j,k) = gamma(i,j,k) + 0.5d0 * gu(i,l) &
                     * (dgg(l,j,k) + dgg(l,k,j) - dgg(j,k,l))
             end do
          end do
       end do
    end do
  end subroutine calc_connections
  
  subroutine calc_connectionderivs (gu, dgg, dgu, ddgg, dgamma)
    CCTK_REAL, intent(in)  :: gu(3,3), dgg(3,3,3), dgu(3,3,3), ddgg(3,3,3,3)
    CCTK_REAL, intent(out) :: dgamma(3,3,3,3)
    integer   :: i,j,k,l,m
    ! Gamma^i_jk,l = 1/2 g^im,l (g_mj,k + g_mk,j - g_jk,m)
    !                + 1/2 g^im (g_mj,kl + g_mk,jl - g_jk,ml)
    do i=1,3
       do j=1,3
          do k=1,3
             do l=1,3
                dgamma(i,j,k,l) = 0
                do m=1,3
                   dgamma(i,j,k,l) = dgamma(i,j,k,l) &
                        + 0.5d0 * dgu(i,m,l) * (dgg(m,j,k) + dgg(m,k,j) - dgg(j,k,m)) &
                        + 0.5d0 * gu(i,m) * (ddgg(m,j,k,l) + ddgg(m,k,j,l) - ddgg(j,k,m,l))
                end do
             end do
          end do
       end do
    end do
  end subroutine calc_connectionderivs
  
  subroutine calc_connectionderivs2 (gu, dgg, dgu, ddgg, ddgu, dddgg, ddgamma)
    CCTK_REAL, intent(in)  :: gu(3,3), dgg(3,3,3), dgu(3,3,3), ddgg(3,3,3,3), ddgu(3,3,3,3), dddgg(3,3,3,3,3)
    CCTK_REAL, intent(out) :: ddgamma(3,3,3,3,3)
    integer   :: i,j,k,l,m,n
    ! Gamma^i_jk,l = 1/2 g^im,l (g_mj,k + g_mk,j - g_jk,m)
    !                + 1/2 g^im (g_mj,kl + g_mk,jl - g_jk,ml)
    do i=1,3
       do j=1,3
          do k=1,3
             do l=1,3
                do m=1,3
                   ddgamma(i,j,k,l,m) = 0
                   do n=1,3
                      ddgamma(i,j,k,l,m) = ddgamma(i,j,k,l,m) &
                           + 0.5d0 * ddgu(i,n,l,m) * (dgg(n,j,k) + dgg(n,k,j) - dgg(j,k,n)) &
                           + 0.5d0 * dgu(i,n,l) * (ddgg(n,j,k,m) + ddgg(n,k,j,m) - ddgg(j,k,n,m)) &
                           + 0.5d0 * dgu(i,n,m) * (ddgg(n,j,k,l) + ddgg(n,k,j,l) - ddgg(j,k,n,l)) &
                           + 0.5d0 * gu(i,n) * (dddgg(n,j,k,l,m) + dddgg(n,k,j,l,m) - dddgg(j,k,n,l,m))
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine calc_connectionderivs2
  
  subroutine calc_ricci (gamma, dgamma, ri)
    CCTK_REAL, intent(in)  :: gamma(3,3,3), dgamma(3,3,3,3)
    CCTK_REAL, intent(out) :: ri(3,3)
#ifdef DEBUG
    CCTK_REAL :: nrm, cnt
#endif
    integer :: i,j,k,l
    ! R_ij =   Gamma^k_ij,k - Gamma^k_ik,j
    !        + Gamma^k_lk Gamma^l_ij - Gamma^k_lj Gamma^l_ki
    do i=1,3
       do j=1,3
          ri(i,j) = 0
          do k=1,3
             ri(i,j) = ri(i,j) + dgamma(k,i,j,k) - dgamma(k,i,k,j)
             do l=1,3
                ri(i,j) = ri(i,j) + gamma(k,l,k) * gamma(l,i,j) &
                     &            - gamma(k,l,j) * gamma(l,k,i)
             end do
          end do
       end do
    end do
#ifdef DEBUG
    ! check symmetries
    sum = 0
    cnt = 0
    do i=1,3
       do j=1,3
          sum = sum + (ri(i,j) - ri(j,i))**2
          cnt = cnt + ri(i,j)**2
       end do
    end do
    if (sum > 1.0e-12 * cnt) call CCTK_WARN (0, "Ricci tensor is not symmetric")
#endif
  end subroutine calc_ricci
  
  subroutine calc_ricciderivs (gamma, dgamma, ddgamma, dri)
    CCTK_REAL, intent(in)  :: gamma(3,3,3), dgamma(3,3,3,3), ddgamma(3,3,3,3,3)
    CCTK_REAL, intent(out) :: dri(3,3,3)
#if 0
    CCTK_REAL :: nrm, cnt
#endif
    integer :: i,j,k,l,m
    ! R_ij =   Gamma^k_ij,k - Gamma^k_ik,j
    !        + Gamma^k_lk Gamma^l_ij - Gamma^k_lj Gamma^l_ki
    do i=1,3
       do j=1,3
          do k=1,3
             dri(i,j,k) = 0
             do l=1,3
                dri(i,j,k) = dri(i,j,k) + ddgamma(l,i,j,l,k) &
                     &                  - ddgamma(l,i,l,j,k)
                do m=1,3
                   dri(i,j,k) = dri(i,j,k) + dgamma(l,m,l,k) * gamma(m,i,j) &
                        &                  + gamma(l,m,l) * dgamma(m,i,j,k) &
                        &                  - dgamma(l,m,j,k) * gamma(m,l,i) &
                        &                  - gamma(l,m,j) * dgamma(m,l,i,k)
                end do
             end do
          end do
       end do
    end do
#if 0
/* Does not work when dri=0 */
    ! check symmetries
    nrm = 0
    cnt = 0
    do i=1,3
       do j=1,3
          do k=1,3
             nrm = nrm + (dri(i,j,k) - dri(j,i,k))**2
             cnt = cnt + dri(i,j,k)**2
          end do
       end do
    end do
    cnt = cnt + sum(gamma**2)
    if (nrm > 1.0e-12 * cnt) then
       print *,"nrm cnt ",nrm,cnt
       print *,"ddgamma ",ddgamma
       print *,"dri ",dri
    end if
    if (nrm > 1.0e-12 * cnt) call CCTK_WARN (0, "Ricci tensor derivative is not symmetric")
#endif
  end subroutine calc_ricciderivs
  
end module ricci
