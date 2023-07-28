#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_calc_tetrad1 (CCTK_ARGUMENTS, hn)
  use adm_metric_simple
  use cctk
  use classify
  use matinv
  use pointwise2
  use qlm_derivs
  use qlm_gram_schmidt
  use qlm_variables
  use ricci4
  use tensor
  use tensor4
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_REAL, parameter :: one=1, two=2
  CCTK_REAL, parameter :: gg4(0:3,0:3,0:3) = 0
  CCTK_REAL    :: gg(3,3), dgg(3,3,3), gg_dot(3,3)
  CCTK_REAL    :: kk(3,3)
  CCTK_REAL    :: alfa, beta(3)
  CCTK_REAL    :: g4(0:3,0:3), gu4(0:3,0:3), dg4(0:3,0:3,0:3)
  CCTK_REAL    :: gamma4(0:3,0:3,0:3)
  CCTK_REAL    :: ee(0:3,0:3), ee_p(0:3,1:3), ee_p_p(0:3,1:3)
  CCTK_REAL    :: dee_spher(1:3,1:3,1:3), ee_inv(1:3,1:3)
  CCTK_REAL    :: dee(0:3,0:3,0:3), gee(0:3,0:3,0:3)
  CCTK_REAL    :: m1(0:3), m2(0:3) ! temporary variables to calculate mm
  CCTK_REAL    :: ss(0:3)          ! spacelike outward normal to horizon
  CCTK_REAL    :: tt(0:3)          ! timelike unit normal to hypersurface
  CCTK_REAL    :: ll(0:3)          ! future null vector on the horizon
  CCTK_REAL    :: nn(0:3)          ! future inward null vector
  CCTK_COMPLEX :: mm(0:3)          ! vector on horizon within the hypersurface
  CCTK_REAL    :: gtt(0:3,0:3), gss(0:3,0:3)
  CCTK_REAL    :: gm1(0:3,0:3), gm2(0:3,0:3)
  CCTK_REAL    :: gll(0:3,0:3), gnn(0:3,0:3)
  CCTK_COMPLEX :: gmm(0:3,0:3)
  CCTK_REAL    :: nabla_ll(0:3,0:3), nabla_nn(0:3,0:3)
  CCTK_COMPLEX :: nabla_mm(0:3,0:3)
  
  !CCTK_REAL    :: t0, t1, t2
  !logical      :: ce0, ce1, ce2
  CCTK_REAL    :: delta_space(2)
  
  CCTK_REAL    :: count, accuracy
  
  integer      :: lsh(2)
  integer      :: i, j
  integer      :: a, b, c, d
  CCTK_REAL    :: theta, phi
  
  logical      :: lerr
  
  character*2, parameter :: crlf = achar(13) // achar(10)
  character    :: msg*1000
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Setting tetrad")
  end if
  
  lsh(:) = (/ qlm_ntheta(hn), qlm_nphi(hn) /)
  delta_space(:) = (/ qlm_delta_theta(hn), qlm_delta_phi(hn) /)
  
  count = 0
  accuracy = 0
  
  ! Calculate the coordinates
  do j = 1+qlm_nghostsphi(hn), qlm_nphi(hn)-qlm_nghostsphi(hn)
     do i = 1+qlm_nghoststheta(hn), qlm_ntheta(hn)-qlm_nghoststheta(hn)
        theta = qlm_origin_theta(hn) + (i-1)*qlm_delta_theta(hn)
        phi   = qlm_origin_phi(hn)   + (j-1)*qlm_delta_phi(hn)
        
        ! Get the variables from the arrays
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
        
        alfa = qlm_alpha(i,j)
        
        beta(1) = qlm_betax(i,j)
        beta(2) = qlm_betay(i,j)
        beta(3) = qlm_betaz(i,j)
        
        
        
        ! Calculate 4-metric
        call calc_3metricdot_simple (kk, gg_dot)
        call calc_4metricderivs_simple (gg,dgg,gg_dot, g4,dg4)
        call calc_4inv (g4, gu4)
        call calc_4connections (gu4,dg4, gamma4)
        
        
        
        ee = TAT_nan()
        dee_spher = TAT_nan()
        dee = TAT_nan()
        
        
        
        ! Calculate the future timelike unit normal vector
        ! t^2 = -1
        ee(0,:) = (/ one, -beta /) / alfa
        dee(0,:,:) = 0
        
        
        
        ee(1,0) = 0
        ee(1,1) = qlm_x(i,j,hn) - qlm_origin_x(hn)
        ee(1,2) = qlm_y(i,j,hn) - qlm_origin_y(hn)
        ee(1,3) = qlm_z(i,j,hn) - qlm_origin_z(hn)
        
        !ee_p(1,1) = qlm_x_p(i,j,hn) - qlm_origin_x_p(hn)
        !ee_p(1,2) = qlm_y_p(i,j,hn) - qlm_origin_y_p(hn)
        !ee_p(1,3) = qlm_z_p(i,j,hn) - qlm_origin_z_p(hn)
        
        !ee_p_p(1,1) = qlm_x_p_p(i,j,hn) - qlm_origin_x_p_p(hn)
        !ee_p_p(1,2) = qlm_y_p_p(i,j,hn) - qlm_origin_y_p_p(hn)
        !ee_p_p(1,3) = qlm_z_p_p(i,j,hn) - qlm_origin_z_p_p(hn)
        
        dee(1,0,:) = 0
        !dee(1,1:3,0) = timederiv (ee(1,1:3), ee_p(1,1:3), ee_p_p(1,1:3), t0,t1,t2, ce0,ce1,ce2)
        dee(1,1:3,0) = 0
        dee_spher(1,:,1) = 0    ! this is a choice
        dee_spher(1,1,2:3) = deriv (qlm_x(:,:,hn), i,j, delta_space)
        dee_spher(1,2,2:3) = deriv (qlm_y(:,:,hn), i,j, delta_space)
        dee_spher(1,3,2:3) = deriv (qlm_z(:,:,hn), i,j, delta_space)
        
        
        
        ee(2:3,0) = 0
        ee(2:3,1) = deriv (qlm_x(:,:,hn), i,j, delta_space)
        ee(2:3,2) = deriv (qlm_y(:,:,hn), i,j, delta_space)
        ee(2:3,3) = deriv (qlm_z(:,:,hn), i,j, delta_space)
        
        ee_p(2:3,1) = deriv (qlm_x_p(:,:,hn), i,j, delta_space)
        ee_p(2:3,2) = deriv (qlm_y_p(:,:,hn), i,j, delta_space)
        ee_p(2:3,3) = deriv (qlm_z_p(:,:,hn), i,j, delta_space)
        
        ee_p_p(2:3,1) = deriv (qlm_x_p_p(:,:,hn), i,j, delta_space)
        ee_p_p(2:3,2) = deriv (qlm_y_p_p(:,:,hn), i,j, delta_space)
        ee_p_p(2:3,3) = deriv (qlm_z_p_p(:,:,hn), i,j, delta_space)
        
        dee(2:3,0,:) = 0
        !dee(2:3,1:3,0) = timederiv (ee(2:3,1:3), ee_p(2:3,1:3), ee_p_p(2:3,1:3), t0,t1,t2, ce0,ce1,ce2)
        dee(2:3,1:3,0) = 0
        dee_spher(2:3,:,1) = 0  ! this is a choice
        dee_spher(2:3,1,2:3) = deriv2 (qlm_x(:,:,hn), i,j, delta_space)
        dee_spher(2:3,2,2:3) = deriv2 (qlm_y(:,:,hn), i,j, delta_space)
        dee_spher(2:3,3,2:3) = deriv2 (qlm_z(:,:,hn), i,j, delta_space)
        
        ! ee_a^i
        ! dee_spher_a^i,b
        ! dee_a^i,j = ee_j^b dee_spher_a^i,b
        call calc_inv3 (ee(1:3,1:3), ee_inv, lerr)
        if (lerr) then
           call CCTK_WARN (3, "Could not invert matrix")
        end if
        do a=1,3
           do b=1,3
              do c=1,3
                 dee(a,b,c) = 0
                 do d=1,3
                    dee(a,b,c) = dee(a,b,c) + ee_inv(c,d) * dee_spher(a,b,d)
                 end do
              end do
           end do
        end do
        
        do a=0,3
           do b=0,3
              gee(:,a,b) = dee(:,a,b) 
              do c=0,3
                 gee(:,a,b) = gee(:,a,b) + ee(:,c) * gamma4(a,c,b)
              end do
           end do
        end do
        
        
        
        ! tt
        tt(:) = ee(0,:)
        gtt(:,:) = gee(0,:,:)
        
        ! m1 = ep
        m1(:) = ee(3,:)
        gm1(:,:) = gee(3,:,:)
        call gram_schmidt_normalise (g4,gg4, m1,gm1, one)
        
        ! m2 = et
        m2(:) = ee(2,:)
        gm2(:,:) = gee(2,:,:)
        call gram_schmidt_project (g4,gg4, m1,gm1, one, m2,gm2)
        call gram_schmidt_normalise (g4,gg4, m2,gm2, one)
        
        ! ss = er
        ss(:) = ee(1,:)
        gss(:,:) = gee(1,:,:)
        call gram_schmidt_project (g4,gg4, m1,gm1, one, ss,gss)
        call gram_schmidt_project (g4,gg4, m2,gm2, one, ss,gss)
        call gram_schmidt_normalise (g4,gg4, ss,gss, one)
        
        
        
        ! ll = (tt + ss) / sqrt(two)
        ll = (tt + ss) / sqrt(two)
        gll = (gtt + gss) / sqrt(two)
        
        ! nn = (tt - ss) / sqrt(two)
        nn = (tt - ss) / sqrt(two)
        gnn = (gtt - gss) / sqrt(two)
        
        ! mm = cmplx(m1, m2, kind(mm)) / sqrt(two)
        mm = cmplx(m1, m2, kind(mm)) / sqrt(two)
        gmm = cmplx(gm1, gm2, kind(gmm)) / sqrt(two)
        
        
        
        ! Store the stuff into the arrays
        do a=0,3
           do b=0,3
              nabla_ll(a,b) = 0
              nabla_nn(a,b) = 0
              nabla_mm(a,b) = 0
              do c=0,3
                 nabla_ll(a,b) = nabla_ll(a,b) + g4(a,c) * gll(c,b)
                 nabla_nn(a,b) = nabla_nn(a,b) + g4(a,c) * gnn(c,b)
                 nabla_mm(a,b) = nabla_mm(a,b) + g4(a,c) * gmm(c,b)
              end do
              qlm_tetrad_derivs(i,j)%nabla_ll(a,b) = nabla_ll(a,b)
              qlm_tetrad_derivs(i,j)%nabla_nn(a,b) = nabla_nn(a,b)
              qlm_tetrad_derivs(i,j)%nabla_mm(a,b) = nabla_mm(a,b)
           end do
        end do
        
        qlm_l0(i,j,hn) = ll(0)
        qlm_l1(i,j,hn) = ll(1)
        qlm_l2(i,j,hn) = ll(2)
        qlm_l3(i,j,hn) = ll(3)
        
        qlm_n0(i,j,hn) = nn(0)
        qlm_n1(i,j,hn) = nn(1)
        qlm_n2(i,j,hn) = nn(2)
        qlm_n3(i,j,hn) = nn(3)
        
        qlm_m0(i,j,hn) = mm(0)
        qlm_m1(i,j,hn) = mm(1)
        qlm_m2(i,j,hn) = mm(2)
        qlm_m3(i,j,hn) = mm(3)
        
     end do
  end do
  
  if (count > 0) then
     accuracy = sqrt(accuracy / count)
  end if
  
  if (veryverbose/=0) then
     write (msg, '("Tetrad accuracy L2 norm: ",g12.4)') accuracy
     call CCTK_INFO (msg)
  end if
  
!!$  if (accuracy > 1.0d-12) then
  if (accuracy > 1.0d-8) then
     call CCTK_WARN (1, "Tetrad is not accurate")
  end if
end subroutine qlm_calc_tetrad1
