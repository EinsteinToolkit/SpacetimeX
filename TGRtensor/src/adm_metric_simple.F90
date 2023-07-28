! $Header$

#include "cctk.h"

module adm_metric_simple
  ! This module assumes that alpha=1 and beta^i=0
  
  use constants
  use tensor
  implicit none
  private
  public calc_4metric_simple
  public calc_4metricderivs_simple
  public calc_4metricderivs2_simple
  
  public calc_3metric_simple
  public calc_3metricderivs_simple
  public calc_extcurv_simple
  public calc_extcurvdot_simple
  
  public calc_3metricdot_simple
  public calc_3metricderivdot_simple
  public calc_3metricdot2_simple
  
contains
  
  subroutine calc_4metric_simple (gg, g4)
    CCTK_REAL, intent(in)  :: gg(3,3)
    CCTK_REAL, intent(out) :: g4(0:3,0:3)
    
    ! ds^2 = - dt^2 + g_ij dx^i dx^j
    g4(0  ,0  ) = -1
    g4(1:3,0  ) = 0
    g4(0  ,1:3) = 0
    g4(1:3,1:3) = gg
  end subroutine calc_4metric_simple
  
  subroutine calc_4metricderivs_simple (gg, dgg, gg_dot, g4,dg4)
    CCTK_REAL, intent(in)  :: gg(3,3)
    CCTK_REAL, intent(in)  :: dgg(3,3,3)
    CCTK_REAL, intent(in)  :: gg_dot(3,3)
    CCTK_REAL, intent(out) :: g4(0:3,0:3),dg4(0:3,0:3,0:3)
    CCTK_REAL :: d4gg(3,3,0:3)
    integer   :: i
    
    ! 4-metric
    g4(0  ,0  ) = -1
    g4(1:3,0  ) = 0
    g4(0  ,1:3) = 0
    g4(1:3,1:3) = gg
    
    ! derivatives
    d4gg  (:,:,0  ) = gg_dot(:,:)
    d4gg  (:,:,1:3) = dgg(:,:,:)
    
    forall (i=0:3)
       dg4(0  ,0  ,i) = 0
       dg4(1:3,0  ,i) = 0
       dg4(0  ,1:3,i) = 0
       dg4(1:3,1:3,i) = d4gg(:,:,i)
    end forall
  end subroutine calc_4metricderivs_simple
  
  subroutine calc_4metricderivs2_simple (gg, dgg, ddgg, gg_dot, gg_dot2, dgg_dot, &
       g4,dg4,ddg4)
    CCTK_REAL, intent(in)  :: gg(3,3)
    CCTK_REAL, intent(in)  :: dgg(3,3,3)
    CCTK_REAL, intent(in)  :: ddgg(3,3,3,3)
    CCTK_REAL, intent(in)  :: gg_dot(3,3)
    CCTK_REAL, intent(in)  :: gg_dot2(3,3)
    CCTK_REAL, intent(in)  :: dgg_dot(3,3,3)
    CCTK_REAL, intent(out) :: g4(0:3,0:3),dg4(0:3,0:3,0:3)
    CCTK_REAL, intent(out) :: ddg4(0:3,0:3,0:3,0:3)
    CCTK_REAL :: d4gg(3,3,0:3)
    CCTK_REAL :: dd4gg(3,3,0:3,0:3)
    integer   :: i,j
    
    ! 4-metric
    g4(0  ,0  ) = -1
    g4(1:3,0  ) = 0
    g4(0  ,1:3) = 0
    g4(1:3,1:3) = gg
    
    ! first derivatives
    d4gg  (:,:,0  ) = gg_dot(:,:)
    d4gg  (:,:,1:3) = dgg(:,:,:)
    
    forall (i=0:3)
       dg4(0  ,0  ,i) = 0
       dg4(1:3,0  ,i) = 0
       dg4(0  ,1:3,i) = 0
       dg4(1:3,1:3,i) = d4gg(:,:,i)
    end forall
    
    ! second derivatives
    dd4gg  (:,:,0  ,0  ) = gg_dot2(:,:)
    dd4gg  (:,:,1:3,0  ) = dgg_dot(:,:,:)
    dd4gg  (:,:,0  ,1:3) = dgg_dot(:,:,:)
    dd4gg  (:,:,1:3,1:3) = ddgg(:,:,:,:)
    
    ! g4(0  ,0  ) = -1
    ! g4(1:3,0  ) = 0
    ! g4(0  ,1:3) = 0
    ! g4(1:3,1:3) = gg
    
    ! dg4(0  ,0  ,i) = 0
    ! dg4(1:3,0  ,i) = 0
    ! dg4(0  ,1:3,i) = 0
    ! dg4(1:3,1:3,i) = d4gg(:,:,i)
    forall (i=0:3, j=0:3)
       ddg4(0  ,0  ,i,j) = 0
       ddg4(1:3,0  ,i,j) = 0
       ddg4(0  ,1:3,i,j) = 0
       ddg4(1:3,1:3,i,j) = dd4gg(:,:,i,j)
    end forall
  end subroutine calc_4metricderivs2_simple
  
  
  
  subroutine calc_3metric_simple (g4, gg)
    CCTK_REAL, intent(in)  :: g4(0:3,0:3)
    CCTK_REAL, intent(out) :: gg(3,3)
    
    ! ds^2 = -alpha^2 dt^2 + g_ij (dx^i + beta^i dt) (dx^j + beta^j dt)
    gg = g4(1:3,1:3)
  end subroutine calc_3metric_simple
  
  subroutine calc_3metricderivs_simple (g4,dg4, gg, dgg, gg_dot)
    CCTK_REAL, intent(in)  :: g4(0:3,0:3),dg4(0:3,0:3,0:3)
    CCTK_REAL, intent(out) :: gg(3,3)
    CCTK_REAL, intent(out) :: dgg(3,3,3)
    CCTK_REAL, intent(out) :: gg_dot(3,3)
    CCTK_REAL :: d4gg(3,3,0:3)
    integer   :: i
    
    ! ds^2 = -alpha^2 dt^2 + g_ij (dx^i + beta^i dt) (dx^j + beta^j dt)
    gg = g4(1:3,1:3)
    
    forall (i=0:3)
       d4gg(:,:,i) = dg4(1:3,1:3,i)
    end forall
    gg_dot = d4gg(:,:,0)
    dgg = d4gg(:,:,1:3)
  end subroutine calc_3metricderivs_simple
  
  subroutine calc_extcurv_simple (gg_dot, kk)
    CCTK_REAL, intent(in)  :: gg_dot(3,3)
    CCTK_REAL, intent(out) :: kk(3,3)
    integer :: i,j
    
    ! d/dt g_ij = -2 K_ij
    
    do i=1,3
       do j=1,3
          kk(i,j) = - gg_dot(i,j)
          kk(i,j) = kk(i,j) / 2
       end do
    end do
    
  end subroutine calc_extcurv_simple
  
  subroutine calc_extcurvdot_simple (gg,gu,ri, kk, tt, kk_dot)
    CCTK_REAL, intent(in)  :: gg(3,3), gu(3,3)
    CCTK_REAL, intent(in)  :: ri(3,3)
    CCTK_REAL, intent(in)  :: kk(3,3)
    CCTK_REAL, intent(in)  :: tt(3,3)
    CCTK_REAL, intent(out) :: kk_dot(3,3)
    integer :: i,j,l,m
    
    ! d/dt K_ij = R_ij - 2 K_il g^lm K_mj + g^lm K_lm K_ij
    !                  - 8 pi (T_ij - 1/2 g_ij T)
    
    do i=1,3
       do j=1,3
          kk_dot(i,j) = ri(i,j) - 8*pi * tt(i,j)
          do l=1,3
             do m=1,3
                kk_dot(i,j) = kk_dot(i,j) &
                     + gu(l,m) * (-2 * kk(i,l) * kk(m,j) + kk(l,m) * kk(i,j)) &
                     + 4*pi * gg(i,j) * gu(l,m) * tt(l,m)
             end do
          end do
       end do
    end do
  end subroutine calc_extcurvdot_simple
  
  
  
  subroutine calc_3metricdot_simple (kk, gg_dot)
    CCTK_REAL, intent(in)  :: kk(3,3)
    CCTK_REAL, intent(out) :: gg_dot(3,3)
    integer :: i,j
    
    ! d/dt g_ij = -2 K_ij
    
    do i=1,3
       do j=1,3
          gg_dot(i,j) = -2 * kk(i,j)
       end do
    end do
  end subroutine calc_3metricdot_simple
  
  subroutine calc_3metricderivdot_simple (dkk, dgg_dot)
    CCTK_REAL, intent(in)  :: dkk(3,3,3)
    CCTK_REAL, intent(out) :: dgg_dot(3,3,3)
    integer :: i,j,k
    
    ! d/dt g_ij,k = -2 K_ij,k
    
    do i=1,3
       do j=1,3
          do k=1,3
             dgg_dot(i,j,k) = -2 * dkk(i,j,k)
          end do
       end do
    end do
  end subroutine calc_3metricderivdot_simple
  
  subroutine calc_3metricdot2_simple (kk_dot, gg_dot2)
    CCTK_REAL, intent(in)  :: kk_dot(3,3)
    CCTK_REAL, intent(out) :: gg_dot2(3,3)
    integer :: i,j
    
    ! d/dt g_ij = -2 K_ij
    ! d^2/dt^2 g_ij = -2 d/dt K_ij
    
    do i=1,3
       do j=1,3
          gg_dot2(i,j) = -2 * kk_dot(i,j)
       end do
    end do
  end subroutine calc_3metricdot2_simple
  
end module adm_metric_simple
