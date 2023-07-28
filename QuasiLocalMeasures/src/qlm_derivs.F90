#include "cctk.h"
#include "cctk_Parameters.h"

module qlm_derivs
  use classify
  implicit none
  private
  public abs2
  public operator(.outer.)
  public operator(.dot.)
  public deriv
  public deriv2
  public timederiv
  
  interface deriv
     module procedure rderiv
  end interface
  
  interface deriv2
     module procedure rderiv2
  end interface
  
  interface timederiv
     module procedure rtimederiv
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
  pure elemental function abs2 (a)
    CCTK_COMPLEX, intent(in) :: a
    CCTK_REAL :: abs2
    abs2 = a * conjg(a)
  end function abs2
  
  
  
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
  
  
  
  ! Calculate spatial derivatives
  
  pure function rderiv (f, i, j, dx) result (df)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(in)           :: f(:,:)
    integer,   intent(in)           :: i, j
    CCTK_REAL, intent(in)           :: dx(2)
    CCTK_REAL                       :: df(2)
    
    select case (spatial_order)
    case (2)
       df(1) = (+ f(i+1,j) - f(i-1,j)) / (2 * dx(1))
       df(2) = (+ f(i,j+1) - f(i,j-1)) / (2 * dx(2))
    case (4)
       df(1) = (- f(i+2,j) + 8*f(i+1,j) - 8*f(i-1,j) + f(i-2,j)) / (12 * dx(1))
       df(2) = (- f(i,j+2) + 8*f(i,j+1) - 8*f(i,j-1) + f(i,j-2)) / (12 * dx(2))
    case default
       ! call CCTK_WARN (0, "internal error")
       df = TAT_nan()
    end select
  end function rderiv
  
  pure function rderiv2 (f, i, j, dx) result (ddf)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(in) :: f(:,:)
    integer,   intent(in) :: i, j
    CCTK_REAL, intent(in) :: dx(2)
    CCTK_REAL             :: ddf(2,2)
    
    select case (spatial_order)
    case (2)
       ddf(1,1) = (+ f(i+1,j) - 2*f(i,j) + f(i-1,j)) / dx(1)**2
       ddf(2,2) = (+ f(i,j+1) - 2*f(i,j) + f(i,j-1)) / dx(2)**2
       ddf(1,2) = (+ f(i-1,j-1) - f(i+1,j-1) - f(i-1,j+1) + f(i+1,j+1)) &
             &    / (4 * dx(1) * dx(2))
       ddf(2,1) = ddf(1,2)
    case (4)
       ddf(1,1) &
            = (- f(i+2,j) + 16*f(i+1,j) - 30*f(i,j) + 16*f(i-1,j) - f(i-2,j)) &
            & / (12 * dx(1)**2)
       ddf(2,2) &
            = (- f(i,j+2) + 16*f(i,j+1) - 30*f(i,j) + 16*f(i,j-1) - f(i,j-2)) &
            & / (12 * dx(2)**2)
       ddf(1,2) &
            = (+   f(i+2,j+2) -  8*f(i+1,j+2) +  8*f(i-1,j+2) -   f(i-2,j+2)  &
            &  - 8*f(i+2,j+1) + 64*f(i+1,j+1) - 64*f(i-1,j+1) + 8*f(i-2,j+1)  &
            &  + 8*f(i+2,j-1) - 64*f(i+1,j-1) + 64*f(i-1,j-1) - 8*f(i-2,j-1)  &
            &  -   f(i+2,j-2) +  8*f(i+1,j-2) -  8*f(i-1,j-2) +   f(i-2,j-2)) &
            & / (144 * dx(1) * dx(2))
       ddf(2,1) = ddf(1,2)
    case default
       ! call CCTK_WARN (0, "internal error")
       ddf = TAT_nan()
    end select
  end function rderiv2
  
  
  
  ! Calculate a time derivate from several time levels
  pure elemental function rtimederiv (f0, f1, f2, t0, t1, t2, ce0, ce1, ce2) result (fdot)
    CCTK_REAL, intent(in) :: f0, f1, f2
    CCTK_REAL, intent(in) :: t0, t1, t2
    logical,   intent(in) :: ce0, ce1, ce2
    CCTK_REAL             :: fdot
    CCTK_REAL :: dt1, dt2
    CCTK_REAL :: fdot1, fdot2
    
!!$    dt1 = qlm_time(hn) - qlm_time_p(hn)
!!$    dt2 = qlm_time(hn) - qlm_time_p_p(hn)
    dt1 = t0 - t1
    dt2 = t0 - t2
    
    if (ce0 .or. ce1) then
       fdot = 0
    else if (ce2) then
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
  end function rtimederiv
  
end module qlm_derivs
