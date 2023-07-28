! $Header$

#include "cctk.h"

module tensor2
  implicit none
  private
  
  public calc_2trace
  
  public calc_2det
  public calc_2detderiv
  public calc_2detdot
  
  public calc_2inv
  public calc_2invderiv
  public calc_2invdot
  
contains
  
  subroutine calc_2trace (kk, gu, trk)
    CCTK_REAL, intent(in)  :: kk(2,2)
    CCTK_REAL, intent(in)  :: gu(2,2)
    CCTK_REAL, intent(out) :: trk
    integer :: i,j
    trk = 0
    do i=1,2
       do j=1,2
          trk = trk + gu(i,j) * kk(i,j)
       end do
    end do
  end subroutine calc_2trace
  
  
  
  subroutine calc_2det (gg, dtg)
    CCTK_REAL, intent(in)  :: gg(2,2)
    CCTK_REAL, intent(out) :: dtg
    dtg =  gg(1,1) * gg(2,2) - gg(1,2)**2
  end subroutine calc_2det
  
  subroutine calc_2pdet (gg, pgg, pdtg)
    CCTK_REAL, intent(in)  :: gg(2,2), pgg(2,2)
    CCTK_REAL, intent(out) :: pdtg
    pdtg = pgg(1,1) * gg(2,2) + gg(1,1) * pgg(2,2) - 2 * gg(1,2) * pgg(1,2)
  end subroutine calc_2pdet
  
  subroutine calc_2detderiv (gg, dgg, ddtg)
    CCTK_REAL, intent(in)  :: gg(2,2), dgg(2,2,2)
    CCTK_REAL, intent(out) :: ddtg(2)
    integer :: i
    do i=1,2
       call calc_2pdet (gg, dgg(:,:,i), ddtg(i))
    end do
  end subroutine calc_2detderiv
  
  subroutine calc_2detdot (gg, gg_dot, dtg_dot)
    CCTK_REAL, intent(in)  :: gg(2,2), gg_dot(2,2)
    CCTK_REAL, intent(out) :: dtg_dot
    call calc_2pdet (gg, gg_dot, dtg_dot)
  end subroutine calc_2detdot
  
  
  
  subroutine calc_2inv (gg, dtg, gu)
    CCTK_REAL, intent(in)  :: gg(2,2), dtg
    CCTK_REAL, intent(out) :: gu(2,2)
    gu(1,1) = gg(2,2) / dtg
    gu(2,2) = gg(1,1) / dtg
    gu(1,2) = gg(1,2) / dtg
    gu(2,1) = gu(1,2)
  end subroutine calc_2inv
  
  subroutine calc_2pinv (gu, pgg, pgu)
    CCTK_REAL, intent(in)  :: gu(2,2), pgg(2,2)
    CCTK_REAL, intent(out) :: pgu(2,2)
    integer :: i,j,k,l
    do i=1,2
       do j=1,2
          pgu(i,j) = 0
          do k=1,2
             do l=1,2
                pgu(i,j) = pgu(i,j) - gu(i,k) * gu(j,l) * pgg(k,l)
             end do
          end do
       end do
    end do
  end subroutine calc_2pinv
  
  subroutine calc_2invderiv (gu, dgg, dgu)
    CCTK_REAL, intent(in)  :: gu(2,2), dgg(2,2,2)
    CCTK_REAL, intent(out) :: dgu(2,2,2)
    integer :: i
    do i=1,2
       call calc_2pinv (gu, dgg(:,:,i), dgu(:,:,i))
    end do
  end subroutine calc_2invderiv
  
  subroutine calc_2invdot (gu, gg_dot, gu_dot)
    CCTK_REAL, intent(in)  :: gu(2,2), gg_dot(2,2)
    CCTK_REAL, intent(out) :: gu_dot(2,2)
    call calc_2pinv (gu, gg_dot, gu_dot)
  end subroutine calc_2invdot
  
end module tensor2
