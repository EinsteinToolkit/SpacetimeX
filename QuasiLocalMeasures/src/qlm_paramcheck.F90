#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_paramcheck (CCTK_ARGUMENTS)
  use cctk
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer   :: hn
  integer   :: sn
  integer   :: eq_theta, new_ntheta, int_nphi, new_nphi
  character :: msg*1000
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Checking parameters")
  end if
  
  do hn = 1, num_surfaces
     
     if (surface_index(hn) == -1 .and. CCTK_EQUALS(surface_name(hn),"")) then
        ! no surface selected; everything is fine
        goto 9999
     end if
     
     if ((surface_index(hn) < 0 .or. surface_index(hn) >= nsurfaces) .and. &
         CCTK_EQUALS(surface_name(hn), "")) then
        write (msg, '("Illegal surface index specified for surface ",i4," (index is ",i4,", must be less than ",i4,")")') hn-1, surface_index(hn), nsurfaces
        call CCTK_PARAMWARN (msg)
        goto 9999
     end if
     
     sn = sf_IdFromName(surface_index(hn), surface_name(hn)) + 1
     
     ! Import surface description
     qlm_nghoststheta(hn) = nghoststheta(sn)
     qlm_nghostsphi  (hn) = nghostsphi(sn)
     qlm_ntheta      (hn) = ntheta(sn)
     qlm_nphi        (hn) = nphi(sn)
     
     ! Symmetries
     if (symmetric_x(sn) /= 0 .or. &
          symmetric_y(sn) /= 0 .or. &
          symmetric_z(sn) /= 0) then
        call CCTK_WARN (CCTK_WARN_ABORT, "SphericalSurface symmetries are not supported")
     end if
     
     if (auto_res(sn) /= 1 .and. (qlm_ntheta(hn) > maxntheta .or. qlm_nphi(hn) > maxnphi)) then
        write (msg, '("Surface ",i4," is too large: shape is (",2i6,"), maximum is (",2i6,")")') hn-1, qlm_ntheta(hn), qlm_nphi(hn), maxntheta, maxnphi
        call CCTK_PARAMWARN (msg)
     end if
     
     if (qlm_nghoststheta(hn)<1 .or. qlm_nghostsphi(hn)<1) then
        write (msg, '("Not enough ghost zones for surface ",i4,": nghosts=",2i4,", minimum is ",2i4)') hn-1, qlm_nghoststheta(hn), qlm_nghostsphi(hn), 1, 1
        call CCTK_PARAMWARN (msg)
     end if
     
     if (auto_res(sn) /= 1 .and. (mod(int(qlm_ntheta(hn) - 2*qlm_nghoststheta(hn)),2) /= 1)) then
        ! We need a grid point on the equator
        write (msg, '("The number of interior grid points in the theta direction of surface ",i4," must be odd after the symmetries have been removed, but it is ",i6)') hn-1, qlm_ntheta(hn) - 2*qlm_nghoststheta(hn)
        call CCTK_PARAMWARN (msg)
     end if
     
     if (auto_res(sn) /= 1 .and. (mod(int(qlm_nphi(hn) - 2*qlm_nghostsphi(hn)),4) /= 0)) then
        ! We need grid points on the four major meridians
        write (msg, '("The number of interior grid points in the phi direction of surface ",i4," must be a multiple of four after the symmetries have been removed, but it is ",i6)') hn-1, qlm_nphi(hn) - 2*qlm_nghostsphi(hn)
        call CCTK_PARAMWARN (msg)
     end if
     
9999 continue
  end do
  
end subroutine qlm_paramcheck
