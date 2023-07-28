! $Header$

#include "cctk.h"

module covderivs
  implicit none
  private
  
  public calc_scalargrad
  public calc_covectorgrad
  public calc_vectorgrad
  public calc_tensorgrad
  public calc_scalargradgrad
  public calc_covectorgradgrad
  public calc_longitudinal
  public calc_gradlongitudinal
  
contains
  
  subroutine calc_scalargrad (f, df, gamma, gf)
    CCTK_REAL, intent(in)  :: f
    CCTK_REAL, intent(in)  :: df(3)
    CCTK_REAL, intent(in)  :: gamma(3,3,3)
    CCTK_REAL, intent(out) :: gf(3)
    integer :: i
    ! f;i = f,i
    do i=1,3
       gf(i) = df(i)
    end do
  end subroutine calc_scalargrad
  
  subroutine calc_covectorgrad (f, df, gamma, gf)
    CCTK_REAL, intent(in)  :: f(3)
    CCTK_REAL, intent(in)  :: df(3,3)
    CCTK_REAL, intent(in)  :: gamma(3,3,3)
    CCTK_REAL, intent(out) :: gf(3,3)
    integer :: i,j,k
    ! f_i;j = f_i,j - Gamma^k_ij f_k
    do i=1,3
       do j=1,3
          gf(i,j) = df(i,j)
          do k=1,3
             gf(i,j) = gf(i,j) - gamma(k,i,j) * f(k)
          end do
       end do
    end do
  end subroutine calc_covectorgrad
  
  subroutine calc_vectorgrad (f, df, gamma, gf)
    CCTK_REAL, intent(in)  :: f(3)
    CCTK_REAL, intent(in)  :: df(3,3)
    CCTK_REAL, intent(in)  :: gamma(3,3,3)
    CCTK_REAL, intent(out) :: gf(3,3)
    integer :: i,j,k
    ! f^i;j = f^i,j + Gamma^i_kj f^k
    do i=1,3
       do j=1,3
          gf(i,j) = df(i,j)
          do k=1,3
             gf(i,j) = gf(i,j) + gamma(i,k,j) * f(k)
          end do
       end do
    end do
  end subroutine calc_vectorgrad
  
  subroutine calc_tensorgrad (f, df, gamma, gf)
    CCTK_REAL, intent(in)  :: f(3,3)
    CCTK_REAL, intent(in)  :: df(3,3,3)
    CCTK_REAL, intent(in)  :: gamma(3,3,3)
    CCTK_REAL, intent(out) :: gf(3,3,3)
    integer :: i,j,k,l
    ! f_ij;k = f_ij,k - Gamma^l_ik f_lj - Gamma^l_jk f_il
    do i=1,3
       do j=1,3
          do k=1,3
             gf(i,j,k) = df(i,j,k)
             do l=1,3
                gf(i,j,k) = gf(i,j,k) - gamma(l,i,k) * f(l,j) &
                     &                - gamma(l,j,k) * f(i,l)
             end do
          end do
       end do
    end do
  end subroutine calc_tensorgrad
  
  
  
  subroutine calc_scalargradgrad (f, df, ddf, gamma, ggf)
    CCTK_REAL, intent(in)  :: f
    CCTK_REAL, intent(in)  :: df(3)
    CCTK_REAL, intent(in)  :: ddf(3,3)
    CCTK_REAL, intent(in)  :: gamma(3,3,3)
    CCTK_REAL, intent(out) :: ggf(3,3)
    integer :: i,j,k
    ! f;ij = f,ij - Gamma^k_ij f,k
    do i=1,3
       do j=1,3
          ggf(i,j) = ddf(i,j)
          do k=1,3
             ggf(i,j) = ggf(i,j) - gamma(k,i,j) * df(k)
          end do
       end do
    end do
  end subroutine calc_scalargradgrad
  
  subroutine calc_covectorgradgrad (f, df, ddf, gamma, dgamma, ggf)
    CCTK_REAL, intent(in)  :: f(3)
    CCTK_REAL, intent(in)  :: df(3,3)
    CCTK_REAL, intent(in)  :: ddf(3,3,3)
    CCTK_REAL, intent(in)  :: gamma(3,3,3), dgamma(3,3,3,3)
    CCTK_REAL, intent(out) :: ggf(3,3,3)
    CCTK_REAL :: gf(3,3)
    integer :: i,j,k,l
    ! f_i;j  = f_i,j - Gamma^l_ij f_l
    ! f_i;jk = f_i,jk - Gamma^l_ij,k f_l - Gamma^l_ij f_l,k
    !          - Gamma^l_ik f_l;j - Gamma^l_jk f_i;l
    call calc_covectorgrad (f, df, gamma, gf)
    do i=1,3
       do j=1,3
          do k=1,3
             ggf(i,j,k) = ddf(i,j,k)
             do l=1,3
                ggf(i,j,k) = ggf(i,j,k) &
                     - dgamma(l,i,j,k) * f(l) - gamma(l,i,j) * df(l,k) &
                     - gamma(l,i,k) * gf(l,j) - gamma(l,j,k) * gf(i,l)
             end do
          end do
       end do
    end do
  end subroutine calc_covectorgradgrad
  
  subroutine calc_vectorgradgrad (f, df, ddf, gamma, dgamma, ggf)
    CCTK_REAL, intent(in)  :: f(3)
    CCTK_REAL, intent(in)  :: df(3,3)
    CCTK_REAL, intent(in)  :: ddf(3,3,3)
    CCTK_REAL, intent(in)  :: gamma(3,3,3), dgamma(3,3,3,3)
    CCTK_REAL, intent(out) :: ggf(3,3,3)
    CCTK_REAL :: gf(3,3)
    integer :: i,j,k,l
    ! f^i;j  = f^i,j + Gamma^i_lj f^l
    ! f^i;jk = f^i,jk + Gamma^i_lj,k f^l + Gamma^i_lj f^l,k
    !          - Gamma^i_lk f^l;j - Gamma^l_jk f^i;l
    call calc_vectorgrad (f, df, gamma, gf)
    do i=1,3
       do j=1,3
          do k=1,3
             ggf(i,j,k) = ddf(i,j,k)
             do l=1,3
                ggf(i,j,k) = ggf(i,j,k) &
                     + dgamma(i,l,j,k) * f(l) + gamma(i,l,j) * df(l,k) &
                     - gamma(i,l,k) * gf(l,j) - gamma(l,j,k) * gf(i,l)
             end do
          end do
       end do
    end do
  end subroutine calc_vectorgradgrad
  
  
  
  subroutine calc_longitudinal (gf, gg, gu, lf)
    CCTK_REAL, intent(in)  :: gf(3,3)
    CCTK_REAL, intent(in)  :: gg(3,3), gu(3,3)
    CCTK_REAL, intent(out) :: lf(3,3)
    integer :: i,j,k,l
    ! (Lf)_ij = D_i f_j + D_j f_i - 2/3 g_ij D^l f_l
    do i=1,3
       do j=1,3
          lf(i,j) = gf(i,j) + gf(j,i)
          do k=1,3
             do l=1,3
                lf(i,j) = lf(i,j) - (2.d0/3.d0) * gg(i,j) * gu(k,l) * gf(k,l)
             end do
          end do
       end do
    end do
  end subroutine calc_longitudinal
  
  subroutine calc_gradlongitudinal (ggf, gg, gu, glf)
    CCTK_REAL, intent(in)  :: ggf(3,3,3)
    CCTK_REAL, intent(in)  :: gg(3,3), gu(3,3)
    CCTK_REAL, intent(out) :: glf(3,3,3)
    integer :: i,j,k,l,m
    ! D_k (Lf)_ij = D_k (D_i f_j + D_j f_i - 2/3 g_ij D^l f_l)
    do i=1,3
       do j=1,3
          do k=1,3
             glf(i,j,k) = ggf(j,i,k) + ggf(i,j,k)
             do l=1,3
                do m=1,3
                   glf(i,j,k) = glf(i,j,k) &
                        - (2.d0/3.d0) * gg(i,j) * gu(l,m) * ggf(l,m,k)
                end do
             end do
          end do
       end do
    end do
  end subroutine calc_gradlongitudinal
  
end module covderivs
