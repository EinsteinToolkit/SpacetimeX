#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_killing_transport (CCTK_ARGUMENTS, hn)
  use cctk
  use constants
  use lapack
  use qlm_boundary
  use qlm_killing_transportation
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  interface
     integer function TAT_isnan (x)
       implicit none
       CCTK_REAL x
     end function TAT_isnan
  end interface
  
  integer, parameter :: lik = lapack_integer_kind

  CCTK_REAL    :: xi(2,3), chi(3)
  CCTK_REAL    :: vec(3,3)
  CCTK_REAL    :: wr(3), wi(3), vl(3,3), vr(3,3)
  integer      :: i0, j0
  integer      :: n
  integer(lik) :: info
  character    :: msg*1000
  
  integer(lik), parameter :: lwork = 100
  CCTK_REAL :: work(lwork)
  
  
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Transporting Killing vector field")
  end if
  
  
  
  ! latitude of "equator"
  i0 = (qlm_ntheta(hn)+1)/2
  ! longitude of zero meridian
  j0 = 1+qlm_nghostsphi(hn)
  
  vec = delta3
  
  if (veryverbose/=0) then
     write (msg, '("Initial vectors (xi^t xi^p chi):")')
     call CCTK_INFO (msg)
     do n=1,3
        write (msg, '(i2,3g18.8)') n, vec(:,n)
        call CCTK_INFO (msg)
     end do
  end if
  
  xi(1,:) = vec(1,:)
  xi(2,:) = vec(2,:)
  chi(:)  = vec(3,:)
  
  do n=1,3
     call transport_along_equator (CCTK_PASS_FTOF, hn, i0, xi(:,n), chi(n))
  end do
  
  vec(1,:) = xi(1,:)
  vec(2,:) = xi(2,:)
  vec(3,:) = chi(:)
  
  if (veryverbose/=0) then
     write (msg, '("Final vectors (xi^t xi^p chi):")')
     call CCTK_INFO (msg)
     do n=1,3
        write (msg, '(i2,3g18.8)') n, vec(:,n)
        call CCTK_INFO (msg)
     end do
  end if
  
  if (TAT_isnan(sum(vec)) /= 0) then
     ! qlm_calc_error(hn) = 1
     qlm_have_killing_vector(hn) = 0
     call CCTK_WARN (3, "There are nans in the final vectors")
     goto 9999
  end if
  
  call geev ('n', 'v', 3_lik, vec, 3_lik, wr, wi, vl, 3_lik, vr, 3_lik, &
       work, lwork, info)
  if (info/=0) then
     ! qlm_calc_error(hn) = 1
     qlm_have_killing_vector(hn) = 0
     write (msg, '("Error in call to GEEV, info=",i2)') info
     call CCTK_WARN (3, msg)
     goto 9999
  end if
  
  if (veryverbose/=0) then
     write (msg, '("Eigenvectors and eigenvalues:")')
     call CCTK_INFO (msg)
     do n=1,3
        write (msg, '(i2,3g18.8,(" (",g18.8,",",g18.8,")"))') &
             n, vr(:,n), wr(n), wi(n)
        call CCTK_INFO (msg)
     end do
  end if
  
  xi(1,:) = vr(1,:)
  xi(2,:) = vr(2,:)
  chi(:)  = vr(3,:)
  
  ! TODO:
  ! This transport scheme is not ideal.
  ! It leads to large fluctuations in the phi direction.
  n=3
  qlm_killing_eigenvalue_re(hn) = wr(n)
  qlm_killing_eigenvalue_im(hn) = wi(n)
  if (abs(cmplx(wr(n),wi(n),kind(wr)) - (1,0)) > 1.0d-4) then
     call CCTK_WARN (3, "Did not manage to find an eigenvector with the eigenvalue 1")
  end if
  call transport_along_equator (CCTK_PASS_FTOF, hn, i0, xi(:,n), chi(n))
  call transport_along_meridians (CCTK_PASS_FTOF, hn, i0)
  
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_xi_t(:,:,hn), -1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_xi_p(:,:,hn), -1)
  call set_boundary (CCTK_PASS_FTOF, hn, qlm_chi (:,:,hn), +1)
  
9999 continue
end subroutine qlm_killing_transport
