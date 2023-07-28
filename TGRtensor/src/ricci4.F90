! $Header$

#include "cctk.h"

#undef DEBUG

module ricci4
  implicit none
  private
  public calc_4connections
  public calc_4connectionderivs
  public calc_4ricci
  public calc_4riemann
  public calc_4weyl
  
contains
  
  subroutine calc_4connections (gu, dg, gamma)
    CCTK_REAL, intent(in)  :: gu(4,4), dg(4,4,4)
    CCTK_REAL, intent(out) :: gamma(4,4,4)
    integer :: i,j,k,l
    ! Gamma^i_jk = 1/2 g^il (g_lj,k + g_lk,j - g_jk,l)
    do i=1,4
       do j=1,4
          do k=1,4
             gamma(i,j,k) = 0
             do l=1,4
                gamma(i,j,k) = gamma(i,j,k) + 0.5d0 * gu(i,l) &
                     * (dg(l,j,k) + dg(l,k,j) - dg(j,k,l))
             end do
          end do
       end do
    end do
  end subroutine calc_4connections
  
  subroutine calc_4connectionderivs (gu, dgg, dgu, ddgg, dgamma)
    CCTK_REAL, intent(in)  :: gu(4,4), dgg(4,4,4), dgu(4,4,4), ddgg(4,4,4,4)
    CCTK_REAL, intent(out) :: dgamma(4,4,4,4)
    integer   :: i,j,k,l,m
    ! Gamma^i_jk,l = 1/2 g^im,l (g_mj,k + g_mk,j - g_jk,m)
    !                + 1/2 g^im (g_mj,kl + g_mk,jl - g_jk,ml)
    do i=1,4
       do j=1,4
          do k=1,4
             do l=1,4
                dgamma(i,j,k,l) = 0
                do m=1,4
                   dgamma(i,j,k,l) = dgamma(i,j,k,l) &
                        + 0.5d0 * dgu(i,m,l) * (dgg(m,j,k) + dgg(m,k,j) - dgg(j,k,m)) &
                        + 0.5d0 * gu(i,m) * (ddgg(m,j,k,l) + ddgg(m,k,j,l) - ddgg(j,k,m,l))
                end do
             end do
          end do
       end do
    end do
  end subroutine calc_4connectionderivs
  
  subroutine calc_4ricci (gamma, dgamma, ri)
    CCTK_REAL, intent(in)  :: gamma(4,4,4), dgamma(4,4,4,4)
    CCTK_REAL, intent(out) :: ri(4,4)
#ifdef DEBUG
    CCTK_REAL :: nrm, cnt
#endif
    integer :: i,j,k,l
#if 1
    ! R_ij =   Gamma^k_ij,k - Gamma^k_ik,j
    !        + Gamma^k_lk Gamma^l_ij - Gamma^k_lj Gamma^l_ki
    do i=1,4
       do j=1,4
          ri(i,j) = 0
          do k=1,4
             ri(i,j) = ri(i,j) + dgamma(k,i,j,k) - dgamma(k,i,k,j)
             do l=1,4
                ri(i,j) = ri(i,j) + gamma(k,l,k) * gamma(l,i,j) &
                     &            - gamma(k,l,j) * gamma(l,k,i)
             end do
          end do
       end do
    end do
#else
    CCTK_REAL :: rmu(4,4,4,4)
    integer   :: m
    do i=1,4
       do j=1,4
          do k=1,4
             do l=1,4
                rmu(i,j,k,l) = dgamma(i,j,l,k) - dgamma(i,j,k,l)
                do m=1,4
                   rmu(i,j,k,l) = rmu(i,j,k,l) &
                        + gamma(m,j,l) * gamma(i,m,k) &
                        - gamma(m,j,k) * gamma(i,m,l)
                end do
             end do
          end do
       end do
    end do
    do i=1,4
       do j=1,4
          ri(i,j) = 0
          do k=1,4
             ri(i,j) = ri(i,j) + rmu(k,i,k,j)
          end do
       end do
    end do
#endif
#ifdef DEBUG
    ! check symmetries
    nrm = 0
    cnt = 0
    do i=1,4
       do j=1,4
          nrm = nrm + (ri(i,j) - ri(j,i))**2
          cnt = cnt + ri(i,j)**2
       end do
    end do
    if (nrm > 1.0e-12 * cnt) call CCTK_WARN (0, "4-Ricci tensor is not symmetric")
#endif
  end subroutine calc_4ricci
  
  subroutine calc_4riemann (gg, gamma, dgamma, rm)
    CCTK_REAL, intent(in)  :: gg(4,4), gamma(4,4,4), dgamma(4,4,4,4)
    CCTK_REAL, intent(out) :: rm(4,4,4,4)
    CCTK_REAL :: rmu(4,4,4,4)
#ifdef DEBUG
    CCTK_REAL :: nrm, cnt
#endif
    integer :: i,j,k,l,m
#if 1
    ! R^i_jkl = + Gamma^i_jl,k - Gamma^i_jk,l
    !           + Gamma^m_jl Gamma^i_mk - Gamma^m_jk Gamma^i_ml
    ! (Wald, 3.4, eqn. (3.4.4), p. 48)
    ! Note that Wald writes R_lkj^i, i.e., we have the order of the
    ! indices exactly reversed.  Due to the Riemann tensor's various
    ! symmetries, this does not change anything.
    do i=1,4
       do j=1,4
          do k=1,4
             do l=1,4
                rmu(i,j,k,l) = dgamma(i,j,l,k) - dgamma(i,j,k,l)
                do m=1,4
                   rmu(i,j,k,l) = rmu(i,j,k,l) &
                        + gamma(m,j,l) * gamma(i,m,k) &
                        - gamma(m,j,k) * gamma(i,m,l)
                end do
             end do
          end do
       end do
    end do
    do i=1,4
       do j=1,4
          do k=1,4
             do l=1,4
                rm(i,j,k,l) = 0
                do m=1,4
                   rm(i,j,k,l) = rm(i,j,k,l) + gg(i,m) * rmu(m,j,k,l)
                end do
             end do
          end do
       end do
    end do
#else
    ! R_ijk^l = + Gamma^l_ik,j - Gamma^l_jk,i
    !           + Gamma^m_ik Gamma^l_mj - Gamma^m_jk Gamma^l_mi
    ! (Wald, 3.4, eqn. (3.4.4), p. 48)
    do i=1,4
       do j=1,4
          do k=1,4
             do l=1,4
                rmu(i,j,k,l) = + dgamma(l,i,k,j) - dgamma(l,j,k,i)
                do m=1,4
                   rmu(i,j,k,l) = rmu(i,j,k,l) &
                        + gamma(m,i,k) * gamma(l,m,j) &
                        - gamma(m,j,k) * gamma(l,m,i)
                end do
             end do
          end do
       end do
    end do
    do i=1,4
       do j=1,4
          do k=1,4
             do l=1,4
                rm(i,j,k,l) = 0
                do m=1,4
                   rm(i,j,k,l) = rm(i,j,k,l) + rmu(i,j,k,m) * gg(m,l)
                end do
             end do
          end do
       end do
    end do
#endif
#ifdef DEBUG
    ! check symmetries
    nrm = 0
    cnt = 0
    do i=1,4
       do j=1,4
          do k=1,4
             do l=1,4
                nrm = nrm + (rm(i,j,k,l) + rm(i,j,l,k))**2
                nrm = nrm + (rm(i,j,k,l) - rm(k,l,i,j))**2
                nrm = nrm + (rm(i,j,k,l) + rm(j,i,k,l))**2
                nrm = nrm + (rm(i,j,k,l) + rm(j,k,i,l) + rm(k,i,j,l))**2
                cnt = cnt + rm(i,j,k,l)**2
             end do
          end do
       end do
    end do
    if (nrm > 1.0e-12 * cnt) call CCTK_WARN (0, "4-Riemann tensor is not symmetric")
#endif
  end subroutine calc_4riemann
  
  subroutine calc_4weyl (gg, rm, ri, rsc, we)
    CCTK_REAL, parameter :: one=1, half=one/2, third=one/3
    CCTK_REAL, intent(in)  :: gg(4,4), rm(4,4,4,4), ri(4,4), rsc
    CCTK_REAL, intent(out) :: we(4,4,4,4)
#ifdef DEBUG
    CCTK_REAL :: nrm, cnt
#endif
    integer :: i,j,k,l
    ! R_ijkl = C_ijkl + 2/(n-2) (g_i[k R_l]j - g_j[k R_l]i)
    !                 - 2/(n-1)(n-2) R g_i[k g_l]j
    ! (Wald, 3.2, eqn. (3.2.28), p. 40)
    do i=1,4
       do j=1,4
          do k=1,4
             do l=1,4
                we(i,j,k,l) = rm(i,j,k,l) &
                     -        (  (gg(i,k) * ri(l,j) - gg(i,l) * ri(k,j)) &
                     &         - (gg(j,k) * ri(l,i) - gg(j,l) * ri(k,i))) &
                     + third * rsc * (gg(i,k) * gg(l,j) - gg(i,l) * gg(k,j))
             end do
          end do
       end do
    end do
#ifdef DEBUG
    ! check symmetries
    nrm = 0
    cnt = 0
    do i=1,4
       do j=1,4
          do k=1,4
             do l=1,4
                nrm = nrm + (we(i,j,k,l) + we(i,j,l,k))**2
                nrm = nrm + (we(i,j,k,l) - we(k,l,i,j))**2
                nrm = nrm + (we(i,j,k,l) + we(j,i,k,l))**2
                nrm = nrm + (we(i,j,k,l) + we(j,k,i,l) + we(k,i,j,l))**2
                cnt = cnt + we(i,j,k,l)**2
             end do
          end do
       end do
    end do
    if (nrm > 1.0e-12 * cnt) call CCTK_WARN (0, "4-Weyl tensor is not symmetric")
#endif
  end subroutine calc_4weyl
  
end module ricci4
