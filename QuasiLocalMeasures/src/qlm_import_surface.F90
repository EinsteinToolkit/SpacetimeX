#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_import_surface (CCTK_ARGUMENTS, hn)
  use cctk
  use qlm_boundary
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  integer   :: i, j, sn, sname_length
  character :: msg*1000, sname*200
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Importing surface shape")
  end if
  
  if ((surface_index(hn) < 0 .or. surface_index(hn) >= nsurfaces) .and. &
      CCTK_EQUALS(surface_name(hn), "")) then
     call CCTK_WARN (0, "Illegal spherical surface index specified")
  end if
  
  sn = sf_IdFromName(surface_index(hn), surface_name(hn)) + 1
  
  if (verbose/=0 .or. veryverbose/=0) then
     call CCTK_FortranString(sname_length, surface_name(hn), sname)
     ! no error checking since FortranString's truncation is sufficient
     write (msg, '("Importing from spherical surface ",i4," ",a)') &
           sn - 1, TRIM(sname)
     call CCTK_INFO (msg)
  end if
  
  
  if (qlm_calc_error(hn) == 0 .and. cctk_iteration > qlm_iteration(hn)) then
     
     ! Cycle time levels
     qlm_have_valid_data_p_p(hn)     = qlm_have_valid_data_p(hn)
     qlm_have_valid_data_p  (hn)     = qlm_have_valid_data  (hn)
     qlm_have_killing_vector_p_p(hn) = qlm_have_killing_vector_p(hn)
     qlm_have_killing_vector_p  (hn) = qlm_have_killing_vector  (hn)
     qlm_iteration(hn)               = cctk_iteration
     
     qlm_time_p_p(hn) = qlm_time_p(hn)
     qlm_time_p  (hn) = qlm_time  (hn)
     qlm_radius_p_p(hn) = qlm_radius_p(hn)
     qlm_radius_p  (hn) = qlm_radius  (hn)
     
     qlm_origin_x_p_p(hn) = qlm_origin_x_p(hn)
     qlm_origin_x_p  (hn) = qlm_origin_x  (hn)
     qlm_origin_y_p_p(hn) = qlm_origin_y_p(hn)
     qlm_origin_y_p  (hn) = qlm_origin_y  (hn)
     qlm_origin_z_p_p(hn) = qlm_origin_z_p(hn)
     qlm_origin_z_p  (hn) = qlm_origin_z  (hn)
     
     qlm_shape_p_p(:,:,hn) = qlm_shape_p(:,:,hn)
     qlm_shape_p  (:,:,hn) = qlm_shape  (:,:,hn)
     
  end if
  
  
  
  ! Check for valid horizon data
  if (sf_valid(sn) <= 0) then
     if (verbose/=0) then
        call CCTK_INFO ("No valid horizon data found")
     end if
     qlm_calc_error(hn) = 1
     qlm_have_valid_data(hn) = 0
     qlm_have_killing_vector(hn) = 0
     return
  end if
  
  
  
  ! Import surface description
  qlm_calc_error(hn) = 0
  qlm_have_valid_data(hn) = 0
  qlm_have_killing_vector(hn) = 1
  
  if (qlm_have_valid_data_p(hn) == 0) then
     qlm_timederiv_order(hn) = 0
  else if (qlm_have_valid_data_p_p(hn) == 0) then
     qlm_timederiv_order(hn) = 1
  else
     qlm_timederiv_order(hn) = 2
  end if
  
  qlm_time(hn) = cctk_time
  
#warning "TODO: Ensure that the surface parameters don't change"
  qlm_nghoststheta(hn) = sf_nghoststheta(sn)
  qlm_nghostsphi  (hn) = sf_nghostsphi(sn)
  qlm_ntheta      (hn) = sf_ntheta(sn)
  qlm_nphi        (hn) = sf_nphi(sn)
  
  qlm_origin_x    (hn) = sf_origin_x(sn)
  qlm_origin_y    (hn) = sf_origin_y(sn)
  qlm_origin_z    (hn) = sf_origin_z(sn)
  qlm_origin_theta(hn) = sf_origin_theta(sn)
  qlm_origin_phi  (hn) = sf_origin_phi(sn)
  qlm_delta_theta (hn) = sf_delta_theta(sn)
  qlm_delta_phi   (hn) = sf_delta_phi(sn)
  
  if (veryverbose /= 0) then
     write (msg, '("calc error      : ",i6)') qlm_calc_error(hn)
     call CCTK_INFO (msg)
     write (msg, '("time deriv order: ",i6)') qlm_timederiv_order(hn)
     call CCTK_INFO (msg)
     write (msg, '("time            : ",g16.6)') qlm_time(hn)
     call CCTK_INFO (msg)
     write (msg, '("nghosts         : ",2i6)') qlm_nghoststheta(hn), qlm_nghostsphi(hn)
     call CCTK_INFO (msg)
     write (msg, '("n               : ",2i6)') qlm_ntheta(hn), qlm_nphi(hn)
     call CCTK_INFO (msg)
     write (msg, '("origin          : ",3g16.6)') &
          qlm_origin_x(hn), qlm_origin_y(hn), qlm_origin_z(hn)
     call CCTK_INFO (msg)
     write (msg, '("origin          : ",2g16.6)') &
          qlm_origin_theta(hn), qlm_origin_phi(hn)
     call CCTK_INFO (msg)
     write (msg, '("delta           : ",2g16.6)') &
          qlm_delta_theta(hn), qlm_delta_phi(hn)
     call CCTK_INFO (msg)
  end if
  
  if (qlm_ntheta(hn) > maxntheta .or. qlm_nphi(hn) > maxnphi) then
     call CCTK_WARN (0, "Surface is too large")
  end if
  
  if (qlm_nghoststheta(hn)<1 .or. qlm_nghostsphi(hn)<1) then
     call CCTK_WARN (0, "Not enough ghost zones for the horizon surface -- need at least 1")
  end if
  
  
  ! Import the surface
  ! Calculate the coordinates
  do j = 1, qlm_nphi(hn)
     do i = 1, qlm_ntheta(hn)
        qlm_shape(i,j,hn) = sf_radius(i,j,sn)
     end do
  end do
  
  
  
  if (mod(int(qlm_ntheta(hn) - 2*qlm_nghoststheta(hn)),2) /= 1) then
     ! We need a grid point on the equator
     call CCTK_WARN (0, "The number of interior grid points in theta direction must be odd")
  end if
  
  if (mod(int(qlm_nphi(hn) - 2*qlm_nghostsphi(hn)),4) /= 0) then
     ! We need grid points on the four major meridians
     call CCTK_WARN (0, "The number of interior grid points in phi direction must be a multiple of four")
  end if
  
end subroutine qlm_import_surface
