! $Header$

#include "cctk.h"

module tensor4
  use matinv
  implicit none
  private
  
  public calc_4trace
  
  public calc_4inv
  public calc_4invderiv
  
contains
  
  subroutine calc_4trace (kk, gu, trk)
    CCTK_REAL, intent(in)  :: kk(0:3,0:3)
    CCTK_REAL, intent(in)  :: gu(0:3,0:3)
    CCTK_REAL, intent(out) :: trk
    integer :: i,j
    trk = 0
    do i=0,3
       do j=0,3
          trk = trk + gu(i,j) * kk(i,j)
       end do
    end do
  end subroutine calc_4trace
  
  
  
  subroutine calc_4inv (gg, gu)
    CCTK_REAL, intent(in)  :: gg(0:3,0:3)
    CCTK_REAL, intent(out) :: gu(0:3,0:3)
    call calc_syminv4 (gg, gu)
  end subroutine calc_4inv
  
  subroutine calc_4pinv (gu, pgg, pgu)
    CCTK_REAL, intent(in)  :: gu(0:3,0:3), pgg(0:3,0:3)
    CCTK_REAL, intent(out) :: pgu(0:3,0:3)
    integer :: i,j,k,l
    do i=0,3
       do j=0,3
          pgu(i,j) = 0
          do k=0,3
             do l=0,3
                pgu(i,j) = pgu(i,j) - gu(i,k) * gu(j,l) * pgg(k,l)
             end do
          end do
       end do
    end do
  end subroutine calc_4pinv
  
  subroutine calc_4invderiv (gu, dgg, dgu)
    CCTK_REAL, intent(in)  :: gu(0:3,0:3), dgg(0:3,0:3,0:3)
    CCTK_REAL, intent(out) :: dgu(0:3,0:3,0:3)
    integer :: i
    do i=0,3
       call calc_4pinv (gu, dgg(:,:,i), dgu(:,:,i))
    end do
  end subroutine calc_4invderiv
  
end module tensor4
