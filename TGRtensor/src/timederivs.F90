! $Header$

#include "cctk.h"

module timederivs
  implicit none
  private
  public abs2
  public timederiv
  public timederiv2
  public timederiv_uneven
  public timederiv2_uneven
  public operator(.outer.)
  public operator(.dot.)
  
  interface timederiv
     module procedure rtimederiv
  end interface
  
  interface timederiv2
     module procedure rtimederiv2
  end interface
  
  interface timederiv_uneven
     module procedure rtimederiv_uneven
     module procedure ctimederiv_uneven
  end interface
  
  interface timederiv2_uneven
     module procedure rtimederiv2_uneven
     module procedure ctimederiv2_uneven
  end interface
  
  interface operator(.outer.)
     module procedure router
     module procedure couter
  end interface
  
  interface operator(.dot.)
     module procedure rdot
     module procedure cdot
  end interface
  
contains
  
  ! abs(c)**2 for complex c without a square root
  elemental function abs2 (a)
    CCTK_COMPLEX, intent(in) :: a
    CCTK_REAL :: abs2
    abs2 = a * conjg(a)
  end function abs2
  
  
  
  ! Calculate a time derivate from several time levels
  elemental function rtimederiv (f0, f1, f2, dt) result (fdot)
    CCTK_REAL, intent(in) :: f0, f1, f2
    CCTK_REAL, intent(in) :: dt
    CCTK_REAL             :: fdot
    CCTK_REAL :: dt1, dt2
    CCTK_REAL :: fdot1, fdot2
    
    dt1 = dt
    dt2 = 2*dt
    
    ! f(dt1) = f(0) + dt1 f'(0) + dt1^2/2 f''(0)
    ! f(dt2) = f(0) + dt2 f'(0) + dt2^2/2 f''(0)
    ! f'(0) = [f(dt1) - f(0)] / dt1 - dt1/2 f''(0)
    ! f'(0) = [f(dt2) - f(0)] / dt2 - dt2/2 f''(0)
    fdot1 = (f0 - f1) / dt1
    fdot2 = (f0 - f2) / dt2
    fdot = (fdot1 * dt2 - fdot2 * dt1) / (dt2 - dt1)
  end function rtimederiv
  
  elemental function rtimederiv2 (f0, f1, f2, dt) result (fdotdot)
    CCTK_REAL, intent(in) :: f0, f1, f2
    CCTK_REAL, intent(in) :: dt
    CCTK_REAL             :: fdotdot
    CCTK_REAL :: dt1, dt2
    CCTK_REAL :: fdotdot1, fdotdot2
    
    dt1 = dt
    dt2 = 2*dt
    
    ! f(dt1) = f(0) + dt1 f'(0) + dt1^2/2 f''(0)
    ! f(dt2) = f(0) + dt2 f'(0) + dt2^2/2 f''(0)
    ! f''(0) = [f(dt1) - f(0)] / [dt1^2/2] - f'(0) / [dt1/2]
    ! f''(0) = [f(dt2) - f(0)] / [dt2^2/2] - f'(0) / [dt2/2]
    fdotdot1 = (f1 - f0) / (dt1**2/2)
    fdotdot2 = (f2 - f0) / (dt2**2/2)
    fdotdot = (fdotdot1 * dt1 - fdotdot2 * dt2) / (dt1 - dt2)
  end function rtimederiv2
  
  
  
  ! Calculate a time derivate from several time levels with uneven spacing
  elemental function rtimederiv_uneven (f0, f1, f2, t0, t1, t2, ce0, ce1, ce2) result (fdot)
    CCTK_REAL, intent(in) :: f0, f1, f2
    CCTK_REAL, intent(in) :: t0, t1, t2
    CCTK_INT,  intent(in) :: ce0, ce1, ce2
    CCTK_REAL             :: fdot
    CCTK_REAL :: dt1, dt2
    CCTK_REAL :: fdot1, fdot2
    
!!$    dt1 = ih_time(hn) - ih_time_p(hn)
!!$    dt2 = ih_time(hn) - ih_time_p_p(hn)
    dt1 = t0 - t1
    dt2 = t0 - t2
    
    if (ce1/=0) then
       fdot = 0
    else if (ce2/=0) then
       fdot = (f0 - f1) / dt1
    else
       ! f(dt1) = f(0) + dt1 f'(0) + dt1^2/2 f''(0)
       ! f(dt2) = f(0) + dt2 f'(0) + dt2^2/2 f''(0)
       ! f'(0) = [f(dt1) - f(0)] / dt1 - dt1/2 f''(0)
       ! f'(0) = [f(dt2) - f(0)] / dt2 - dt2/2 f''(0)
       fdot1 = (f0 - f1) / dt1
       fdot2 = (f0 - f2) / dt2
       fdot = (fdot1 * dt2 - fdot2 * dt1) / (dt2 - dt1)
    end if
  end function rtimederiv_uneven
  
  elemental function ctimederiv_uneven (f0, f1, f2, t0, t1, t2, ce0, ce1, ce2) result (fdot)
    CCTK_COMPLEX, intent(in) :: f0, f1, f2
    CCTK_REAL,    intent(in) :: t0, t1, t2
    CCTK_INT,     intent(in) :: ce0, ce1, ce2
    CCTK_COMPLEX             :: fdot
    
    fdot = cmplx(timederiv_uneven(real(f0),real(f1),real(f2), t0,t1,t2, ce0,ce1,ce2), &
         &       timederiv_uneven(aimag(f0),aimag(f1),aimag(f2), t0,t1,t2, ce0,ce1,ce2), &
         &       kind(fdot))
  end function ctimederiv_uneven
  
  elemental function rtimederiv2_uneven (f0, f1, f2, t0, t1, t2, ce0, ce1, ce2) result (fdotdot)
    CCTK_REAL, intent(in) :: f0, f1, f2
    CCTK_REAL, intent(in) :: t0, t1, t2
    CCTK_INT,  intent(in) :: ce0, ce1, ce2
    CCTK_REAL             :: fdotdot
    CCTK_REAL :: dt1, dt2
    CCTK_REAL :: fdotdot1, fdotdot2
    
!!$    dt1 = ih_time(hn) - ih_time_p(hn)
!!$    dt2 = ih_time(hn) - ih_time_p_p(hn)
    dt1 = t0 - t1
    dt2 = t0 - t2
    
    if (ce1/=0) then
       fdotdot = 0
    else if (ce2/=0) then
       fdotdot = 0
    else
       ! f(dt1) = f(0) + dt1 f'(0) + dt1^2/2 f''(0)
       ! f(dt2) = f(0) + dt2 f'(0) + dt2^2/2 f''(0)
       ! f''(0) = [f(dt1) - f(0)] / [dt1^2/2] - f'(0) / [dt1/2]
       ! f''(0) = [f(dt2) - f(0)] / [dt2^2/2] - f'(0) / [dt2/2]
       fdotdot1 = (f1 - f0) / (dt1**2/2)
       fdotdot2 = (f2 - f0) / (dt2**2/2)
       fdotdot = (fdotdot1 * dt1 - fdotdot2 * dt2) / (dt1 - dt2)
    end if
  end function rtimederiv2_uneven
  
  elemental function ctimederiv2_uneven (f0, f1, f2, t0, t1, t2, ce0, ce1, ce2) result (fdotdot)
    CCTK_COMPLEX, intent(in) :: f0, f1, f2
    CCTK_REAL,    intent(in) :: t0, t1, t2
    CCTK_INT,     intent(in) :: ce0, ce1, ce2
    CCTK_COMPLEX             :: fdotdot
    
    fdotdot = cmplx(timederiv2_uneven(real(f0),real(f1),real(f2), t0,t1,t2, ce0,ce1,ce2), &
         &          timederiv2_uneven(aimag(f0),aimag(f1),aimag(f2), t0,t1,t2, ce0,ce1,ce2), &
         &          kind(fdotdot))
  end function ctimederiv2_uneven
  
  
  
  function router (left, right) result (result)
    CCTK_REAL, intent(in) :: left(:), right(:)
    CCTK_REAL             :: result(size(left,1),size(right,1))
    integer :: i, j
    forall (i=1:size(left,1), j=1:size(right,1))
       result(i,j) = left(i) * right(j)
    end forall
  end function router
  
  function couter (left, right) result (result)
    CCTK_COMPLEX, intent(in) :: left(:), right(:)
    CCTK_COMPLEX             :: result(size(left,1),size(right,1))
    integer :: i, j
    forall (i=1:size(left,1), j=1:size(right,1))
       result(i,j) = left(i) * right(j)
    end forall
  end function couter
  
  
  
  function rdot (left, right) result (result)
    CCTK_REAL, intent(in) :: left(:), right(:)
    CCTK_REAL             :: result
    integer :: i
    if (size(left,1) /= size(right,1)) then
       call CCTK_WARN (0, "Array sizes must have the same sizes")
    end if
!!$    result = sum((/( left(i) * right(i), i=1,size(left,1) )/))
    result = 0
    do i=1,size(left,1)
       result = result + left(i) * right(i)
    end do
  end function rdot
  
  function cdot (left, right) result (result)
    CCTK_COMPLEX, intent(in) :: left(:), right(:)
    CCTK_COMPLEX             :: result
    integer :: i
    if (size(left,1) /= size(right,1)) then
       call CCTK_WARN (0, "Array sizes must have the same sizes")
    end if
!!$    result = sum((/( left(i) * right(i), i=1,size(left,1) )/))
    result = 0
    do i=1,size(left,1)
       result = result + left(i) * right(i)
    end do
  end function cdot
  
end module timederivs
