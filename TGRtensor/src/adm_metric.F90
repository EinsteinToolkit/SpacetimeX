! $Header$

#include "cctk.h"

module adm_metric
  use constants
  use tensor
  implicit none
  private
  public calc_4metric
  public calc_4metricderivs
  public calc_4metricderivs2
  public calc_3metric
  public calc_3metricderivs
  
  public calc_extcurv
  public calc_3metricdot
  public calc_3metricderivdot
!!$  public calc_extcurvdot
!!$  public calc_3metricdot2
  
!!$  public simplify_t4
  
contains
  
  subroutine calc_4metric (gg, alfa, beta, g4)
    CCTK_REAL, intent(in)  :: gg(3,3), alfa, beta(3)
    CCTK_REAL, intent(out) :: g4(0:3,0:3)
    CCTK_REAL :: betal(3)
    
    ! ds^2 = -alpha^2 dt^2 + g_ij (dx^i + beta^i dt) (dx^j + beta^j dt)
    betal = matmul(gg, beta)
    g4(0  ,0  ) = -alfa**2 + sum(betal*beta)
    g4(1:3,0  ) = betal
    g4(0  ,1:3) = betal
    g4(1:3,1:3) = gg
  end subroutine calc_4metric
  
  subroutine calc_4metricderivs (gg,alfa,beta, dgg,dalfa,dbeta, &
       gg_dot,alfa_dot,beta_dot, g4,dg4)
    CCTK_REAL, intent(in)  :: gg(3,3),alfa,beta(3)
    CCTK_REAL, intent(in)  :: dgg(3,3,3),dalfa(3),dbeta(3,3)
    CCTK_REAL, intent(in)  :: gg_dot(3,3),alfa_dot,beta_dot(3)
    CCTK_REAL, intent(out) :: g4(0:3,0:3),dg4(0:3,0:3,0:3)
    CCTK_REAL :: d4gg(3,3,0:3),d4alfa(0:3),d4beta(3,0:3)
    CCTK_REAL :: betal(3),d4betal(3,0:3)
    integer   :: i,j
    
    ! 4-metric
    forall (i=1:3)
       betal(i) = sum(gg(i,:) * beta(:))
    end forall
    g4(0  ,0  ) = -alfa**2 + sum(betal(:) * beta(:))
    g4(1:3,0  ) = betal
    g4(0  ,1:3) = betal
    g4(1:3,1:3) = gg
    
    ! derivatives
    d4gg  (:,:,0  ) = gg_dot(:,:)
    d4gg  (:,:,1:3) = dgg(:,:,:)
    d4alfa(    0  ) = alfa_dot
    d4alfa(    1:3) = dalfa(:)
    d4beta(:,  0  ) = beta_dot(:)
    d4beta(:,  1:3) = dbeta(:,:)
    
    forall (i=1:3, j=0:3)
       d4betal(i,j) = sum(d4gg(i,:,j) * beta(:) + gg(i,:) * d4beta(:,j))
    end forall
    
    forall (i=0:3)
       dg4(0  ,0  ,i) = - 2 * alfa * d4alfa(i) &
            &           + sum(d4betal(:,i) * beta(:) + betal(:) * d4beta(:,i))
       dg4(1:3,0  ,i) = d4betal(:,i)
       dg4(0  ,1:3,i) = d4betal(:,i)
       dg4(1:3,1:3,i) = d4gg(:,:,i)
    end forall
  end subroutine calc_4metricderivs
  
  subroutine calc_4metricderivs2 (gg,alfa,beta, dgg,dalfa,dbeta, &
       ddgg,ddalfa,ddbeta, gg_dot,alfa_dot,beta_dot, &
       gg_dot2,alfa_dot2,beta_dot2, dgg_dot,dalfa_dot,dbeta_dot, g4,dg4,ddg4)
    CCTK_REAL, intent(in)  :: gg(3,3),alfa,beta(3)
    CCTK_REAL, intent(in)  :: dgg(3,3,3),dalfa(3),dbeta(3,3)
    CCTK_REAL, intent(in)  :: ddgg(3,3,3,3),ddalfa(3,3),ddbeta(3,3,3)
    CCTK_REAL, intent(in)  :: gg_dot(3,3),alfa_dot,beta_dot(3)
    CCTK_REAL, intent(in)  :: gg_dot2(3,3),alfa_dot2,beta_dot2(3)
    CCTK_REAL, intent(in)  :: dgg_dot(3,3,3),dalfa_dot(3),dbeta_dot(3,3)
    CCTK_REAL, intent(out) :: g4(0:3,0:3),dg4(0:3,0:3,0:3)
    CCTK_REAL, intent(out) :: ddg4(0:3,0:3,0:3,0:3)
    CCTK_REAL :: d4gg(3,3,0:3),d4alfa(0:3),d4beta(3,0:3)
    CCTK_REAL :: dd4gg(3,3,0:3,0:3),dd4alfa(0:3,0:3),dd4beta(3,0:3,0:3)
    CCTK_REAL :: betal(3),d4betal(3,0:3),dd4betal(3,0:3,0:3)
    integer   :: i,j,k
    
    ! 4-metric
    forall (i=1:3)
       betal(i) = sum(gg(i,:) * beta(:))
    end forall
    g4(0  ,0  ) = -alfa**2 + sum(betal(:) * beta(:))
    g4(1:3,0  ) = betal
    g4(0  ,1:3) = betal
    g4(1:3,1:3) = gg
    
    ! first derivatives
    d4gg  (:,:,0  ) = gg_dot(:,:)
    d4gg  (:,:,1:3) = dgg(:,:,:)
    d4alfa(    0  ) = alfa_dot
    d4alfa(    1:3) = dalfa(:)
    d4beta(:,  0  ) = beta_dot(:)
    d4beta(:,  1:3) = dbeta(:,:)
    
    forall (i=1:3, j=0:3)
       d4betal(i,j) = sum(d4gg(i,:,j) * beta(:) + gg(i,:) * d4beta(:,j))
    end forall
    
    forall (i=0:3)
       dg4(0  ,0  ,i) = - 2 * alfa * d4alfa(i) &
            &           + sum(d4betal(:,i) * beta(:) + betal(:) * d4beta(:,i))
       dg4(1:3,0  ,i) = d4betal(:,i)
       dg4(0  ,1:3,i) = d4betal(:,i)
       dg4(1:3,1:3,i) = d4gg(:,:,i)
    end forall
    
    ! second derivatives
    dd4gg  (:,:,0  ,0  ) = gg_dot2(:,:)
    dd4gg  (:,:,1:3,0  ) = dgg_dot(:,:,:)
    dd4gg  (:,:,0  ,1:3) = dgg_dot(:,:,:)
    dd4gg  (:,:,1:3,1:3) = ddgg(:,:,:,:)
    dd4alfa(    0  ,0  ) = alfa_dot2
    dd4alfa(    1:3,0  ) = dalfa_dot(:)
    dd4alfa(    0  ,1:3) = dalfa_dot(:)
    dd4alfa(    1:3,1:3) = ddalfa(:,:)
    dd4beta(:,  0  ,0  ) = beta_dot2(:)
    dd4beta(:,  1:3,0  ) = dbeta_dot(:,:)
    dd4beta(:,  0  ,1:3) = dbeta_dot(:,:)
    dd4beta(:,  1:3,1:3) = ddbeta(:,:,:)
    
    ! betal(i) = gg(i,m) * beta(m)
    ! d4betal(i,j) = d4gg(i,m,j) * beta(m) + gg(i,m) * d4beta(m,j)
    ! dd4betal(i,j,k) =   dd4gg(i,m,j,k) * beta(m) + d4gg(i,m,j) * d4beta(m,k)
    !                   + d4gg(i,m,k) * d4beta(m,j) + gg(i,m) * dd4beta(m,j,k)
    forall (i=1:3, j=0:3, k=0:3)
       dd4betal(i,j,k) = sum(+ dd4gg(i,:,j,k) *    beta(:)      &
            &                +  d4gg(i,:,j)   *  d4beta(:,k)    &
            &                +  d4gg(i,:,k)   *  d4beta(:,j)    &
            &                +    gg(i,:)     * dd4beta(:,j,k))
    end forall
    
    ! g4(0  ,0  ) = -alfa**2 + sum(betal*beta)
    ! g4(1:3,0  ) = betal
    ! g4(0  ,1:3) = betal
    ! g4(1:3,1:3) = gg
    
    ! dg4(0  ,0  ,i) = -2*alfa*d4alfa(i) &
    !      &           + d4betal(m,i)*beta(m) + betal(m)*d4beta(m,i)
    ! dg4(1:3,0  ,i) = d4betal(:,i)
    ! dg4(0  ,1:3,i) = d4betal(:,i)
    ! dg4(1:3,1:3,i) = d4gg(:,:,i)
    forall (i=0:3, j=0:3)
       ddg4(0  ,0  ,i,j) = - 2 * d4alfa(j) * d4alfa(i)        &
            &              - 2 * alfa * dd4alfa(i,j)          &
            &              + sum(+ dd4betal(:,i,j) * beta(:)  &
            &                    + d4betal(:,i) * d4beta(:,j) &
            &                    + d4betal(:,j) * d4beta(:,i) &
            &                    + betal(:) * dd4beta(:,i,j))
       ddg4(1:3,0  ,i,j) = dd4betal(:,i,j)
       ddg4(0  ,1:3,i,j) = dd4betal(:,i,j)
       ddg4(1:3,1:3,i,j) = dd4gg(:,:,i,j)
    end forall
  end subroutine calc_4metricderivs2
  
  
  
  subroutine calc_3metric (g4, gg, alfa, beta)
    CCTK_REAL, intent(in)  :: g4(0:3,0:3)
    CCTK_REAL, intent(out) :: gg(3,3), alfa, beta(3)
    CCTK_REAL :: betal(3)
    CCTK_REAL :: dtg, gu(3,3)
    
    ! ds^2 = -alpha^2 dt^2 + g_ij (dx^i + beta^i dt) (dx^j + beta^j dt)
    betal = g4(1:3,0)
    gg = g4(1:3,1:3)
    call calc_det (gg, dtg)
    call calc_inv (gg, dtg, gu)
    beta = matmul(gu, betal)
    alfa = sqrt(sum(betal*beta) - g4(0,0))
  end subroutine calc_3metric
  
  subroutine calc_3metricderivs (g4,dg4, gg,alfa,beta, dgg,dalfa,dbeta, &
       gg_dot,alfa_dot,beta_dot)
    CCTK_REAL, intent(in)  :: g4(0:3,0:3),dg4(0:3,0:3,0:3)
    CCTK_REAL, intent(out) :: gg(3,3),alfa,beta(3)
    CCTK_REAL, intent(out) :: dgg(3,3,3),dalfa(3),dbeta(3,3)
    CCTK_REAL, intent(out) :: gg_dot(3,3),alfa_dot,beta_dot(3)
    CCTK_REAL :: betal(3),d4betal(3,0:3)
    CCTK_REAL :: dtg,gu(3,3),dgu(3,3,3),gu_dot(3,3)
    CCTK_REAL :: d4gg(3,3,0:3),d4gu(3,3,0:3)
    CCTK_REAL :: d4alfa(0:3),d4beta(3,0:3)
    integer   :: i,j
    
    ! ds^2 = -alpha^2 dt^2 + g_ij (dx^i + beta^i dt) (dx^j + beta^j dt)
    gg = g4(1:3,1:3)
    call calc_det (gg, dtg)
    call calc_inv (gg, dtg, gu)
    betal = g4(1:3,0)
    beta = matmul(gu, betal)
    alfa = sqrt(sum(betal*beta) - g4(0,0))
    
    forall (i=0:3)
       d4gg(:,:,i) = dg4(1:3,1:3,i)
    end forall
    gg_dot = d4gg(:,:,0)
    dgg = d4gg(:,:,1:3)
    
    call calc_invderiv (gu, dgg, dgu)
    call calc_invdot (gu, gg_dot, gu_dot)
    d4gu(:,:,0) = gu_dot
    d4gu(:,:,1:3) = dgu
    
    forall (i=0:3)
       d4betal(:,i) = dg4(1:3,0,i)
    end forall
    forall (i=0:3, j=1:3)
       d4beta(j,i) = sum(d4gu(j,:,i)*betal) + sum(gu(j,:)*d4betal(:,i))
    end forall
    forall (i=0:3)
       d4alfa(i) = 1/(2*alfa) &
            * (sum(d4betal(:,i)*beta) + sum(betal*d4beta(:,i)) - dg4(0,0,i))
    end forall
    
    alfa_dot = d4alfa(0)
    dalfa = d4alfa(1:3)
    beta_dot = d4beta(:,0)
    dbeta = d4beta(:,1:3)
  end subroutine calc_3metricderivs
  
  
  
  subroutine calc_extcurv (gg, dgg, gg_dot, alfa, beta, dbeta, kk)
    CCTK_REAL, intent(in)  :: gg(3,3), dgg(3,3,3), gg_dot(3,3)
    CCTK_REAL, intent(in)  :: alfa
    CCTK_REAL, intent(in)  :: beta(3), dbeta(3,3)
    CCTK_REAL, intent(out) :: kk(3,3)
    integer :: i,j,k
    
    ! d/dt g_ij = -2 alpha K_ij + g_kj beta^k,i + g_ik beta^k,j + beta^k g_ij,k
    
    do i=1,3
       do j=1,3
          kk(i,j) = - gg_dot(i,j)
          do k=1,3
             kk(i,j) = kk(i,j) &
                  + gg(k,j) * dbeta(k,i) + gg(i,k) * dbeta(k,j) &
                  + beta(k) * dgg(i,j,k)
          end do
          kk(i,j) = kk(i,j) / (2*alfa)
       end do
    end do
    
  end subroutine calc_extcurv
  
  subroutine calc_3metricdot (gg, dgg, kk, alfa, beta, dbeta, gg_dot)
    CCTK_REAL, intent(in)  :: gg(3,3), dgg(3,3,3)
    CCTK_REAL, intent(in)  :: kk(3,3)
    CCTK_REAL, intent(in)  :: alfa
    CCTK_REAL, intent(in)  :: beta(3), dbeta(3,3)
    CCTK_REAL, intent(out) :: gg_dot(3,3)
    integer :: i,j,k
    
    ! d/dt g_ij = -2 alpha K_ij + g_kj beta^k,i + g_ik beta^k,j + beta^k g_ij,k
    
    do i=1,3
       do j=1,3
          gg_dot(i,j) = -2*alfa * kk(i,j)
          do k=1,3
             gg_dot(i,j) = gg_dot(i,j) &
                  + gg(k,j) * dbeta(k,i) + gg(i,k) * dbeta(k,j) &
                  + beta(k) * dgg(i,j,k)
          end do
       end do
    end do
  end subroutine calc_3metricdot
  
  subroutine calc_3metricderivdot (gg, dgg, ddgg, kk, dkk, &
       alfa, dalfa, beta, dbeta, ddbeta, dgg_dot)
    CCTK_REAL, intent(in)  :: gg(3,3), dgg(3,3,3), ddgg(3,3,3,3)
    CCTK_REAL, intent(in)  :: kk(3,3), dkk(3,3,3)
    CCTK_REAL, intent(in)  :: alfa, dalfa(3)
    CCTK_REAL, intent(in)  :: beta(3), dbeta(3,3), ddbeta(3,3,3)
    CCTK_REAL, intent(out) :: dgg_dot(3,3,3)
    integer :: i,j,k,l
    
    ! d/dt g_ij = -2 alpha K_ij + g_kj beta^k,i + g_ik beta^k,j + beta^k g_ij,k
    
    ! d/dt g_ij,k = - 2 alpha,k K_ij - 2 alpha K_ij,k
    !               + g_lj,k beta^l,i + g_lj beta^l,ik
    !               + g_il,k beta^l,j + g_il beta^l,jk
    !               + beta^l,k g_ij,l + beta^l g_ij,lk
    
    do i=1,3
       do j=1,3
          do k=1,3
             dgg_dot(i,j,k) = - 2*dalfa(k) * kk(i,j) - 2*alfa * dkk(i,j,k)
             do l=1,3
                dgg_dot(i,j,k) = dgg_dot(i,j,k) &
                     + dgg(l,j,k) * dbeta(l,i) + gg(l,j) * ddbeta(l,i,k) &
                     + dgg(i,l,k) * dbeta(l,j) + gg(i,l) * ddbeta(l,j,k) &
                     + dbeta(l,k) * dgg(i,j,l) + beta(l) * ddgg(i,j,l,k)
             end do
          end do
       end do
    end do
  end subroutine calc_3metricderivdot
  
  
  
  
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
  
  
  
!!$  subroutine simplify_t4 (t4, alfa, beta, t4s)
!!$    CCTK_REAL, intent(in)  :: t4(0:3,0:3)
!!$    CCTK_REAL, intent(in)  :: alfa
!!$    CCTK_REAL, intent(in)  :: beta(3)
!!$    CCTK_REAL, intent(out) :: t4s(0:3,0:3)
!!$    
!!$    ! T_ab   g_ij alfa beta^i
!!$    ! T'_ab   g_ij 1 [0,0,0]
!!$    
!!$    ! x^a
!!$    ! x'^a
!!$    
!!$    ! g'_ab = dx^c/dx'^a dx^d/dx'^b g_cd
!!$    !       = dx^0/dx'^a dx^0/dx'^b g_00
!!$    !         dx^0/dx'^a dx^l/dx'^b g_0l
!!$    !         dx^k/dx'^a dx^0/dx'^b g_k0
!!$    !         dx^k/dx'^a dx^l/dx'^b g_kl
!!$    
!!$    ! g'_ab = J^c_a J^d_b g_cd
!!$    ! T'_ab = J^c_a J^d_b T_cd
!!$    
!!$    
!!$    
!!$    ! delta^i_b = g'^ia g'_ab = J^c_a J^d_b g_cd g'^ia
!!$    
!!$    
!!$    
!!$    ! g' = J^T g J
!!$    ! g' J^(-1) = J^T g
!!$    ! T' = J^T T J
!!$    
!!$    
!!$    
!!$    ! -1 = g'_00 = dx^0/dx'^0 dx^0/dx'^0 g_00
!!$    !              dx^0/dx'^0 dx^l/dx'^0 g_0l
!!$    !              dx^k/dx'^0 dx^0/dx'^0 g_k0
!!$    !              dx^k/dx'^0 dx^l/dx'^0 g_kl
!!$    
!!$    ! 0 = g'_0j = dx^0/dx'^0 dx^0/dx'^j g_00
!!$    !             dx^0/dx'^0 dx^l/dx'^j g_0l   !!! <<< !!!
!!$    !             dx^k/dx'^0 dx^0/dx'^j g_k0
!!$    !             dx^k/dx'^0 dx^l/dx'^j g_kl
!!$    
!!$    ! g_ij = g'_ij = dx^0/dx'^i dx^0/dx'^j g_00
!!$    !                dx^0/dx'^i dx^l/dx'^j g_0l
!!$    !                dx^k/dx'^i dx^0/dx'^j g_k0
!!$    !                dx^k/dx'^i dx^l/dx'^j g_kl
!!$    
!!$    ! dx^0/dx'^0 = 1/sqrt(-g_00)
!!$    ! dx^0/dx'^i = 0
!!$    ! dx^k/dx'^0 = 0
!!$    ! dx^k/dx'^i = delta^k_i
!!$    
!!$    
!!$    
!!$    ! g_00 = -alfa^2 + gama_ij beta^i beta^j
!!$    ! g_0i = gama_ij beta^j
!!$    
!!$    ! alfa^2 = -g_00 + gama_ij beta^i beta^j
!!$    ! beta^i = gama^ij g_0j
!!$    
!!$    ! alfa'^2 = -g'_00 + gama'_ij beta'^i beta'^j
!!$    ! beta'^i = gama'^ij g'_0j
!!$    
!!$    ! -1 = alfa'^2 = -g'_00 + gama'_ij beta'^i beta'^j
!!$    !  0 = beta'^i = gama'^ij g'_0j
!!$    
!!$    ! g'_00 = -1
!!$    ! g'_0i =  0
!!$    
!!$  end subroutine simplify_t4
  
end module adm_metric
