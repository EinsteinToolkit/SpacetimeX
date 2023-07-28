! $Header$

#include "cctk.h"

module tensor
  implicit none
  private
  
  public calc_trace
  
  public calc_det
  public calc_detderiv
  public calc_detdot
  
  public calc_inv
  public calc_invderiv
  public calc_invderiv2
  public calc_invdot
  
contains
  
  subroutine calc_trace (kk, gu, trk)
    CCTK_REAL, intent(in)  :: kk(3,3)
    CCTK_REAL, intent(in)  :: gu(3,3)
    CCTK_REAL, intent(out) :: trk
    integer :: i,j
    trk = 0
    do i=1,3
       do j=1,3
          trk = trk + gu(i,j) * kk(i,j)
       end do
    end do
  end subroutine calc_trace
  
  
  
  subroutine calc_det (gg, dtg)
    CCTK_REAL, intent(in)  :: gg(3,3)
    CCTK_REAL, intent(out) :: dtg
    dtg =      gg(1,1) * gg(2,2) * gg(3,3) &
         & + 2*gg(1,2) * gg(1,3) * gg(2,3) &
         & -   gg(1,1) * gg(2,3)**2        &
         & -   gg(2,2) * gg(1,3)**2        &
         & -   gg(3,3) * gg(1,2)**2
  end subroutine calc_det
  
  subroutine calc_pdet (gg, pgg, pdtg)
    CCTK_REAL, intent(in)  :: gg(3,3), pgg(3,3)
    CCTK_REAL, intent(out) :: pdtg
    pdtg =     pgg(1,1) *   gg(2,2) *  gg(3,3) &
         & +    gg(1,1) *  pgg(2,2) *  gg(3,3) &
         & +    gg(1,1) *   gg(2,2) * pgg(3,3) &
         & + 2*pgg(1,2) *   gg(1,3) *  gg(2,3) &
         & + 2* gg(1,2) *  pgg(1,3) *  gg(2,3) &
         & + 2* gg(1,2) *   gg(1,3) * pgg(2,3) &
         & -   pgg(1,1) *   gg(2,3)**2         &
         & -    gg(1,1) * 2*gg(2,3) * pgg(2,3) &
         & -   pgg(2,2) *   gg(1,3)**2         &
         & -    gg(2,2) * 2*gg(1,3) * pgg(1,3) &
         & -   pgg(3,3) *   gg(1,2)**2         &
         & -    gg(3,3) * 2*gg(1,2) * pgg(1,2)
  end subroutine calc_pdet
  
  subroutine calc_detderiv (gg, dgg, ddtg)
    CCTK_REAL, intent(in)  :: gg(3,3), dgg(3,3,3)
    CCTK_REAL, intent(out) :: ddtg(3)
    integer :: i
    do i=1,3
       call calc_pdet (gg, dgg(:,:,i), ddtg(i))
    end do
  end subroutine calc_detderiv
  
  subroutine calc_detdot (gg, gg_dot, dtg_dot)
    CCTK_REAL, intent(in)  :: gg(3,3), gg_dot(3,3)
    CCTK_REAL, intent(out) :: dtg_dot
    call calc_pdet (gg, gg_dot, dtg_dot)
  end subroutine calc_detdot
  
  
  
  subroutine calc_inv (gg, dtg, gu)
    CCTK_REAL, intent(in)  :: gg(3,3), dtg
    CCTK_REAL, intent(out) :: gu(3,3)
    gu(1,1) = (gg(2,2) * gg(3,3) - gg(2,3) ** 2     ) / dtg
    gu(2,2) = (gg(1,1) * gg(3,3) - gg(1,3) ** 2     ) / dtg
    gu(3,3) = (gg(1,1) * gg(2,2) - gg(1,2) ** 2     ) / dtg
    gu(1,2) = (gg(1,3) * gg(2,3) - gg(1,2) * gg(3,3)) / dtg
    gu(1,3) = (gg(1,2) * gg(2,3) - gg(1,3) * gg(2,2)) / dtg
    gu(2,3) = (gg(1,3) * gg(1,2) - gg(2,3) * gg(1,1)) / dtg
    gu(2,1) = gu(1,2)
    gu(3,1) = gu(1,3)
    gu(3,2) = gu(2,3)
  end subroutine calc_inv
  
  subroutine calc_pinv (gu, pgg, pgu)
    CCTK_REAL, intent(in)  :: gu(3,3), pgg(3,3)
    CCTK_REAL, intent(out) :: pgu(3,3)
    integer :: i,j,k,l
    do i=1,3
       do j=1,3
          pgu(i,j) = 0
          do k=1,3
             do l=1,3
                pgu(i,j) = pgu(i,j) - gu(i,k) * gu(j,l) * pgg(k,l)
             end do
          end do
       end do
    end do
  end subroutine calc_pinv
  
  subroutine calc_invderiv (gu, dgg, dgu)
    CCTK_REAL, intent(in)  :: gu(3,3), dgg(3,3,3)
    CCTK_REAL, intent(out) :: dgu(3,3,3)
    integer :: i
    do i=1,3
       call calc_pinv (gu, dgg(:,:,i), dgu(:,:,i))
    end do
  end subroutine calc_invderiv
  
  subroutine calc_invderiv2 (gu, dgg, dgu, ddgg, ddgu)
    CCTK_REAL, intent(in)  :: gu(3,3), dgg(3,3,3), dgu(3,3,3), ddgg(3,3,3,3)
    CCTK_REAL, intent(out) :: ddgu(3,3,3,3)
    integer :: i,j,k,l,m,n
    do i=1,3
       do j=1,3
          do k=1,3
             do l=1,3
                ddgu(i,j,k,l) = 0
                do m=1,3
                   do n=1,3
                      ddgu(i,j,k,l) = ddgu(i,j,k,l) &
                           - dgu(i,m,l) * gu(j,n) * dgg(m,n,k) &
                           - gu(i,m) * dgu(j,n,l) * dgg(m,n,k) &
                           - gu(i,m) * gu(j,n) * ddgg(m,n,k,l)
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine calc_invderiv2
  
  subroutine calc_invdot (gu, gg_dot, gu_dot)
    CCTK_REAL, intent(in)  :: gu(3,3), gg_dot(3,3)
    CCTK_REAL, intent(out) :: gu_dot(3,3)
    call calc_pinv (gu, gg_dot, gu_dot)
  end subroutine calc_invdot
  
end module tensor
