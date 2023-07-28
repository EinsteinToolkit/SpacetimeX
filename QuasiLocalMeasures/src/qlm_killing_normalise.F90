#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_killing_normalise (CCTK_ARGUMENTS, hn)
  use cctk
  use qlm_killing_normalisation
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_REAL, parameter :: one = 1
  
  integer, parameter :: ngeodesics = 7
  
  CCTK_REAL :: theta
  CCTK_REAL :: factor, factors(ngeodesics)
  logical   :: found(ngeodesics)
  
  integer   :: nsteps
  integer   :: ii
  integer   :: i, j
  
  character :: msg*1000
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Normalising Killing vector field")
  end if
  
  
  
  ! Starting meridian
  j = 1+qlm_nghostsphi(hn)
  
  
  
  if (veryverbose/=0) call CCTK_INFO ("Calculating normalisation factor:")
  
  do ii=1,ngeodesics
     i = 1 + qlm_nghoststheta(hn) &
          + ii * (qlm_ntheta(hn) - 2*qlm_nghoststheta(hn) - 1) / (ngeodesics+1)
     if (qlm_xi_p(i,j,hn) < 0) then
        qlm_xi_t(:,:,hn) = -qlm_xi_t(:,:,hn)
        qlm_xi_p(:,:,hn) = -qlm_xi_p(:,:,hn)
        qlm_chi (:,:,hn) = -qlm_chi (:,:,hn)
     end if
     theta = qlm_origin_theta(hn) + (i-1) * qlm_delta_theta(hn)
     call killing_factor (CCTK_PASS_FTOF, hn, theta, factors(ii), nsteps)
     found(ii) = nsteps>0
  end do
  if (any(found)) then
     if (CCTK_EQUALS(killing_vector_normalisation, "average")) then
        if (count(found) >= 3) then
           ! ignore smallest and largest factors
           ii = minloc(factors, 1, found)
           if (ii<=0) call CCTK_WARN (0, "internal error")
           found(ii) = .false.
           ii = maxloc(factors, 1, found)
           if (ii<=0) call CCTK_WARN (0, "internal error")
           found(ii) = .false.
        end if
     else if (CCTK_EQUALS(killing_vector_normalisation, "median")) then
        do while (count(found) >= 3)
           ! ignore smallest and largest factors
           ii = minloc(factors, 1, found)
           if (ii<=0) call CCTK_WARN (0, "internal error")
           found(ii) = .false.
           ii = maxloc(factors, 1, found)
           if (ii<=0) call CCTK_WARN (0, "internal error")
           found(ii) = .false.
        end do
     else
        call CCTK_WARN (0, "internal error")
     end if
     ! average factors
     factor = product(factors, found) ** (one / count(found))
  else
     call CCTK_WARN (1, "Did not manage to integrate along a Killing vector field line loop")
     ! qlm_calc_error(hn) = 1
     qlm_have_killing_vector(hn) = 0
     factor = 1
     goto 9999
  end if
  
  
  
  if (veryverbose/=0) then
     write (msg, '("Normalising xi with the factor ",g16.6)') factor
     call CCTK_INFO (msg)
  end if
  qlm_xi_t(:,:,hn) = qlm_xi_t(:,:,hn) * factor
  qlm_xi_p(:,:,hn) = qlm_xi_p(:,:,hn) * factor
  qlm_chi (:,:,hn) = qlm_chi (:,:,hn) * factor
  
  
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Checking normalisation factors for various Killing vector field line loops:")
     do ii=1,ngeodesics
        i = 1 + qlm_nghoststheta(hn) &
             + ii * (qlm_ntheta(hn) - 2*qlm_nghoststheta(hn) - 1) / (ngeodesics+1)
        theta = qlm_origin_theta(hn) + (i-1) * qlm_delta_theta(hn)
        call killing_factor (CCTK_PASS_FTOF, hn, theta, factor, nsteps)
     end do
  end if
  
  ! qlm_have_killing_vector(hn) = 1
  
9999 continue
  
end subroutine qlm_killing_normalise
