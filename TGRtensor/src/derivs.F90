! $Header$

#ifndef TGR_INCLUDED

#include "cctk.h"
#include "cctk_Parameters.h"

module derivs
  implicit none
  private
  public get_derivs
  public get_derivs2
  public get_derivs3
contains
#endif
  subroutine get_derivs (a, f, pos, off, dx, order)
    CCTK_REAL, intent(in)            :: a(*)
    CCTK_REAL, intent(out)           :: f(3)
    integer,   intent(in)            :: pos, off(3)
    CCTK_REAL, intent(in)            :: dx(3)
    integer,   intent(in),  optional :: order
    integer :: order1
    integer :: i
    order1 = 2
    if (present(order)) order1 = order
    select case (order1)
    case (2)
       do i=1,3
          f(i) = (a(pos+off(i)) - a(pos-off(i))) / (2*dx(i))
       end do
    case (4)
       do i=1,3
          f(i) = (- a(pos+2*off(i)) + 8*a(pos+off(i)) - 8*a(pos-off(i)) + a(pos-2*off(i))) / (12*dx(i))
       end do
    case default
       call CCTK_WARN (0, "Unsupported finite differencing order")
    end select
  end subroutine get_derivs
  subroutine get_derivs2 (a, f, pos, off, dx, order)
    CCTK_REAL, intent(in)            :: a(*)
    CCTK_REAL, intent(out)           :: f(3,3)
    integer,   intent(in)            :: pos, off(3)
    CCTK_REAL, intent(in)            :: dx(3)
    integer,   intent(in),  optional :: order
    integer :: order1
    integer :: i
    order1 = 2
    if (present(order)) order1 = order
    select case (order1)
    case (2)
       do i=1,3
          f(i,i) = (a(pos+off(i)) - 2*a(pos) + a(pos-off(i))) / dx(i)**2
       end do
       f(1,2) = (  a(pos-off(1)-off(2)) - a(pos+off(1)-off(2)) &
            &    - a(pos-off(1)+off(2)) + a(pos+off(1)+off(2))) / (4*dx(1)*dx(2))
       f(2,1) = f(1,2)
       f(1,3) = (  a(pos-off(1)-off(3)) - a(pos+off(1)-off(3)) &
            &    - a(pos-off(1)+off(3)) + a(pos+off(1)+off(3))) / (4*dx(1)*dx(3))
       f(3,1) = f(1,3)
       f(2,3) = (  a(pos-off(2)-off(3)) - a(pos+off(2)-off(3)) &
            &    - a(pos-off(2)+off(3)) + a(pos+off(2)+off(3))) / (4*dx(2)*dx(3))
       f(3,2) = f(2,3)
    case (4)
       do i=1,3
          f(i,i) = (- a(pos-2*off(i)) + 16*a(pos-off(i)) - 30*a(pos) + 16*a(pos+off(i)) - a(pos+2*off(i))) / (12*dx(i)**2)
       end do
       f(1,2) = (    a(pos+2*off(1)+2*off(2)) -  8*a(pos+off(1)+2*off(2)) +  8*a(pos-off(1)+2*off(2)) -   a(pos-2*off(1)+2*off(2)) &
            &    - 8*a(pos+2*off(1)+  off(2)) + 64*a(pos+off(1)+  off(2)) - 64*a(pos-off(1)+  off(2)) + 8*a(pos-2*off(1)+  off(2)) &
            &    + 8*a(pos+2*off(1)-  off(2)) - 64*a(pos+off(1)-  off(2)) + 64*a(pos-off(1)-  off(2)) - 8*a(pos-2*off(1)-  off(2)) &
            &    -   a(pos+2*off(1)-2*off(2)) +  8*a(pos+off(1)-2*off(2)) -  8*a(pos-off(1)-2*off(2)) +   a(pos-2*off(1)-2*off(2))) / (144*dx(1)*dx(2))
       f(2,1) = f(1,2)
       f(1,3) = (    a(pos+2*off(1)+2*off(3)) -  8*a(pos+off(1)+2*off(3)) +  8*a(pos-off(1)+2*off(3)) -   a(pos-2*off(1)+2*off(3)) &
            &    - 8*a(pos+2*off(1)+  off(3)) + 64*a(pos+off(1)+  off(3)) - 64*a(pos-off(1)+  off(3)) + 8*a(pos-2*off(1)+  off(3)) &
            &    + 8*a(pos+2*off(1)-  off(3)) - 64*a(pos+off(1)-  off(3)) + 64*a(pos-off(1)-  off(3)) - 8*a(pos-2*off(1)-  off(3)) &
            &    -   a(pos+2*off(1)-2*off(3)) +  8*a(pos+off(1)-2*off(3)) -  8*a(pos-off(1)-2*off(3)) +   a(pos-2*off(1)-2*off(3))) / (144*dx(1)*dx(3))
       f(3,1) = f(1,3)
       f(2,3) = (    a(pos+2*off(2)+2*off(3)) -  8*a(pos+off(2)+2*off(3)) +  8*a(pos-off(2)+2*off(3)) -   a(pos-2*off(2)+2*off(3)) &
            &    - 8*a(pos+2*off(2)+  off(3)) + 64*a(pos+off(2)+  off(3)) - 64*a(pos-off(2)+  off(3)) + 8*a(pos-2*off(2)+  off(3)) &
            &    + 8*a(pos+2*off(2)-  off(3)) - 64*a(pos+off(2)-  off(3)) + 64*a(pos-off(2)-  off(3)) - 8*a(pos-2*off(2)-  off(3)) &
            &    -   a(pos+2*off(2)-2*off(3)) +  8*a(pos+off(2)-2*off(3)) -  8*a(pos-off(2)-2*off(3)) +   a(pos-2*off(2)-2*off(3))) / (144*dx(2)*dx(3))
       f(3,2) = f(2,3)
    case default
       call CCTK_WARN (0, "Unsupported finite differencing order")
    end select
  end subroutine get_derivs2
  subroutine get_derivs3 (a, f, pos, off, dx, order)
    CCTK_REAL, intent(in)            :: a(*)
    CCTK_REAL, intent(out)           :: f(3,3,3)
    integer,   intent(in)            :: pos, off(3)
    CCTK_REAL, intent(in)            :: dx(3)
    integer,   intent(in),  optional :: order
    integer :: order1
    order1 = 2
    if (present(order)) order1 = order
    select case (order1)
    case (2)
       f(1,1,1) = (a(pos+2*off(1)) - 2*a(pos+off(1)) + 2*a(pos-off(1)) - a(pos-2*off(1))) / (2*dx(1)**3)
       f(2,1,1) = (a(pos+off(1)+off(2)) - 2*a(pos+off(2)) + a(pos-off(1)+off(2)) - a(pos+off(1)-off(2)) + 2*a(pos-off(2)) - a(pos-off(1)-off(2))) / (2*dx(1)*dx(1)*dx(2))
       f(3,1,1) = (a(pos+off(1)+off(3)) - 2*a(pos+off(3)) + a(pos-off(1)+off(3)) - a(pos+off(1)-off(3)) + 2*a(pos-off(3)) - a(pos-off(1)-off(3))) / (2*dx(1)*dx(1)*dx(3))
       f(1,2,1) = f(2,1,1)
       f(2,2,1) = (a(pos+off(1)+off(2)) - 2*a(pos+off(1)) + a(pos+off(1)-off(2)) - a(pos-off(1)+off(2)) + 2*a(pos-off(1)) - a(pos-off(1)-off(2))) / (2*dx(1)*dx(2)*dx(2))
       f(3,2,1) = (a(pos+off(1)+off(2)+off(3)) - a(pos-off(1)+off(2)+off(3)) - a(pos+off(1)-off(2)+off(3)) + a(pos-off(1)-off(2)+off(3)) - a(pos+off(1)+off(2)-off(3)) + a(pos-off(1)+off(2)-off(3)) + a(pos+off(1)-off(2)-off(3)) - a(pos-off(1)-off(2)-off(3))) / (8*dx(1)*dx(2)*dx(3))
       f(1,3,1) = f(3,1,1)
       f(2,3,1) = f(3,2,1)
       f(3,3,1) = (a(pos+off(1)+off(3)) - 2*a(pos+off(1)) + a(pos+off(1)-off(3)) - a(pos-off(1)+off(3)) + 2*a(pos-off(1)) - a(pos-off(1)-off(3))) / (2*dx(1)*dx(3)*dx(3))
       f(:,1,2) = f(:,2,1)
       f(1,2,2) = f(2,2,1)
       f(2,2,2) = (a(pos+2*off(2)) - 2*a(pos+off(2)) + 2*a(pos-off(2)) - a(pos-2*off(2))) / (2*dx(2)**3)
       f(3,2,2) = (a(pos+off(2)+off(3)) - 2*a(pos+off(3)) + a(pos-off(2)+off(3)) - a(pos+off(2)-off(3)) + 2*a(pos-off(3)) - a(pos-off(2)-off(3))) / (2*dx(2)*dx(2)*dx(3))
       f(1,3,2) = f(3,1,2)
       f(2,3,2) = f(3,2,2)
       f(3,3,2) = (a(pos+off(2)+off(3)) - 2*a(pos+off(2)) + a(pos+off(2)-off(3)) - a(pos-off(2)+off(3)) + 2*a(pos-off(2)) - a(pos-off(2)-off(3))) / (2*dx(2)*dx(3)*dx(3))
       f(:,1,3) = f(:,3,1)
       f(:,2,3) = f(:,3,2)
       f(1,3,3) = f(3,1,3)
       f(2,3,3) = f(3,2,3)
       f(3,3,3) = (a(pos+2*off(3)) - 2*a(pos+off(3)) + 2*a(pos-off(3)) - a(pos-2*off(3))) / (2*dx(3)**3)
    case (4)
       ! [+1 -8 +13 0 -13 +8 -1] / 8
       f(1,1,1) = (a(pos-3*off(1)) - 8*a(pos-2*off(1)) + 13*a(pos-off(1)) - 13*a(pos+off(1)) + 8*a(pos+2*off(1)) - a(pos+3*off(1))) / (8*dx(1)**3)
       ! [-1 +8 0 -8 +1] / 12
       ! [-1 +16 -30 +16 -1] / 12
       f(2,1,1) = (+   a(pos-2*off(1)-2*off(2)) -  16*a(pos-off(1)-2*off(2)) +  30*a(pos-2*off(2)) -  16*a(pos+off(1)-2*off(2)) +   a(pos+2*off(1)-2*off(2)) &
            &      - 8*a(pos-2*off(1)-  off(2)) + 128*a(pos-off(1)-  off(2)) - 240*a(pos-  off(2)) + 128*a(pos+off(1)-  off(2)) - 8*a(pos+2*off(1)-  off(2)) &
            &      + 8*a(pos-2*off(1)+  off(2)) - 128*a(pos-off(1)+  off(2)) + 240*a(pos+  off(2)) - 128*a(pos+off(1)+  off(2)) + 8*a(pos+2*off(1)+  off(2)) &
            &      -   a(pos-2*off(1)+2*off(2)) +  16*a(pos-off(1)+2*off(2)) -  30*a(pos+2*off(2)) +  16*a(pos+off(1)+2*off(2)) -   a(pos+2*off(1)+2*off(2))) / (144*dx(1)**2*dx(2))
       f(3,1,1) = (+   a(pos-2*off(1)-2*off(3)) -  16*a(pos-off(1)-2*off(3)) +  30*a(pos-2*off(3)) -  16*a(pos+off(1)-2*off(3)) +   a(pos+2*off(1)-2*off(3)) &
            &      - 8*a(pos-2*off(1)-  off(3)) + 128*a(pos-off(1)-  off(3)) - 240*a(pos-  off(3)) + 128*a(pos+off(1)-  off(3)) - 8*a(pos+2*off(1)-  off(3)) &
            &      + 8*a(pos-2*off(1)+  off(3)) - 128*a(pos-off(1)+  off(3)) + 240*a(pos+  off(3)) - 128*a(pos+off(1)+  off(3)) + 8*a(pos+2*off(1)+  off(3)) &
            &      -   a(pos-2*off(1)+2*off(3)) +  16*a(pos-off(1)+2*off(3)) -  30*a(pos+2*off(3)) +  16*a(pos+off(1)+2*off(3)) -   a(pos+2*off(1)+2*off(3))) / (144*dx(1)**2*dx(3))
       f(1,2,1) = f(2,1,1)
       f(2,2,1) = (+   a(pos-2*off(2)-2*off(1)) -  16*a(pos-off(2)-2*off(1)) +  30*a(pos-2*off(1)) -  16*a(pos+off(2)-2*off(1)) +   a(pos+2*off(2)-2*off(1)) &
            &      - 8*a(pos-2*off(2)-  off(1)) + 128*a(pos-off(2)-  off(1)) - 240*a(pos-  off(1)) + 128*a(pos+off(2)-  off(1)) - 8*a(pos+2*off(2)-  off(1)) &
            &      + 8*a(pos-2*off(2)+  off(1)) - 128*a(pos-off(2)+  off(1)) + 240*a(pos+  off(1)) - 128*a(pos+off(2)+  off(1)) + 8*a(pos+2*off(2)+  off(1)) &
            &      -   a(pos-2*off(2)+2*off(1)) +  16*a(pos-off(2)+2*off(1)) -  30*a(pos+2*off(1)) +  16*a(pos+off(2)+2*off(1)) -   a(pos+2*off(2)+2*off(1))) / (144*dx(2)**2*dx(1))
       f(3,2,1) = (-     a(pos-2*off(1)-2*off(2)-2*off(3)) +   8*a(pos-off(1)-2*off(2)-2*off(3)) -   8*a(pos+off(1)-2*off(2)-2*off(3)) +     a(pos+2*off(1)-2*off(2)-2*off(3)) &
            &      +   8*a(pos-2*off(1)-  off(2)-2*off(3)) -  64*a(pos-off(1)-  off(2)-2*off(3)) +  64*a(pos+off(1)-  off(2)-2*off(3)) -   8*a(pos+2*off(1)-  off(2)-2*off(3)) &
            &      -   8*a(pos-2*off(1)+  off(2)-2*off(3)) +  64*a(pos-off(1)+  off(2)-2*off(3)) -  64*a(pos+off(1)+  off(2)-2*off(3)) +   8*a(pos+2*off(1)+  off(2)-2*off(3)) &
            &      +     a(pos-2*off(1)+2*off(2)-2*off(3)) -   8*a(pos-off(1)+2*off(2)-2*off(3)) +   8*a(pos+off(1)+2*off(2)-2*off(3)) -     a(pos+2*off(1)+2*off(2)-2*off(3)) &
            &      +   8*a(pos-2*off(1)-2*off(2)-  off(3)) -  64*a(pos-off(1)-2*off(2)-  off(3)) +  64*a(pos+off(1)-2*off(2)-  off(3)) -   8*a(pos+2*off(1)-2*off(2)-  off(3)) &
            &      -  64*a(pos-2*off(1)-  off(2)-  off(3)) + 512*a(pos-off(1)-  off(2)-  off(3)) - 512*a(pos+off(1)-  off(2)-  off(3)) +  64*a(pos+2*off(1)-  off(2)-  off(3)) &
            &      +  64*a(pos-2*off(1)+  off(2)-  off(3)) - 512*a(pos-off(1)+  off(2)-  off(3)) + 512*a(pos+off(1)+  off(2)-  off(3)) -  64*a(pos+2*off(1)+  off(2)-  off(3)) &
            &      -   8*a(pos-2*off(1)+2*off(2)-  off(3)) +  64*a(pos-off(1)+2*off(2)-  off(3)) -  64*a(pos+off(1)+2*off(2)-  off(3)) +   8*a(pos+2*off(1)+2*off(2)-  off(3)) &
            &      -   8*a(pos-2*off(1)-2*off(2)+  off(3)) +  64*a(pos-off(1)-2*off(2)+  off(3)) -  64*a(pos+off(1)-2*off(2)+  off(3)) +   8*a(pos+2*off(1)-2*off(2)+  off(3)) &
            &      +  64*a(pos-2*off(1)-  off(2)+  off(3)) - 512*a(pos-off(1)-  off(2)+  off(3)) + 512*a(pos+off(1)-  off(2)+  off(3)) -  64*a(pos+2*off(1)-  off(2)+  off(3)) &
            &      -  64*a(pos-2*off(1)+  off(2)+  off(3)) + 512*a(pos-off(1)+  off(2)+  off(3)) - 512*a(pos+off(1)+  off(2)+  off(3)) +  64*a(pos+2*off(1)+  off(2)+  off(3)) &
            &      +   8*a(pos-2*off(1)+2*off(2)+  off(3)) -  64*a(pos-off(1)+2*off(2)+  off(3)) +  64*a(pos+off(1)+2*off(2)+  off(3)) -   8*a(pos+2*off(1)+2*off(2)+  off(3)) &
            &      +     a(pos-2*off(1)-2*off(2)+2*off(3)) -   8*a(pos-off(1)-2*off(2)+2*off(3)) +   8*a(pos+off(1)-2*off(2)+2*off(3)) -     a(pos+2*off(1)-2*off(2)+2*off(3)) &
            &      -   8*a(pos-2*off(1)-  off(2)+2*off(3)) +  64*a(pos-off(1)-  off(2)+2*off(3)) -  64*a(pos+off(1)-  off(2)+2*off(3)) +   8*a(pos+2*off(1)-  off(2)+2*off(3)) &
            &      +   8*a(pos-2*off(1)+  off(2)+2*off(3)) -  64*a(pos-off(1)+  off(2)+2*off(3)) +  64*a(pos+off(1)+  off(2)+2*off(3)) -   8*a(pos+2*off(1)+  off(2)+2*off(3)) &
            &      -     a(pos-2*off(1)+2*off(2)+2*off(3)) +   8*a(pos-off(1)+2*off(2)+2*off(3)) -   8*a(pos+off(1)+2*off(2)+2*off(3)) +     a(pos+2*off(1)+2*off(2)+2*off(3))) / (512*dx(1)*dx(2)*dx(3))
       f(1,3,1) = f(3,1,1)
       f(2,3,1) = f(3,2,1)
       f(3,3,1) = (+   a(pos-2*off(3)-2*off(1)) -  16*a(pos-off(3)-2*off(1)) +  30*a(pos-2*off(1)) -  16*a(pos+off(3)-2*off(1)) +   a(pos+2*off(3)-2*off(1)) &
            &      - 8*a(pos-2*off(3)-  off(1)) + 128*a(pos-off(3)-  off(1)) - 240*a(pos-  off(1)) + 128*a(pos+off(3)-  off(1)) - 8*a(pos+2*off(3)-  off(1)) &
            &      + 8*a(pos-2*off(3)+  off(1)) - 128*a(pos-off(3)+  off(1)) + 240*a(pos+  off(1)) - 128*a(pos+off(3)+  off(1)) + 8*a(pos+2*off(3)+  off(1)) &
            &      -   a(pos-2*off(3)+2*off(1)) +  16*a(pos-off(3)+2*off(1)) -  30*a(pos+2*off(1)) +  16*a(pos+off(3)+2*off(1)) -   a(pos+2*off(3)+2*off(1))) / (144*dx(3)**2*dx(1))
       f(:,1,2) = f(:,2,1)
       f(1,2,2) = f(2,2,1)
       f(2,2,2) = (a(pos-3*off(2)) - 8*a(pos-2*off(2)) + 13*a(pos-off(2)) - 13*a(pos+off(2)) + 8*a(pos+2*off(2)) - a(pos+3*off(2))) / (8*dx(2)**3)
       f(3,2,2) = (+   a(pos-2*off(2)-2*off(3)) -  16*a(pos-off(2)-2*off(3)) +  30*a(pos-2*off(3)) -  16*a(pos+off(2)-2*off(3)) +   a(pos+2*off(2)-2*off(3)) &
            &      - 8*a(pos-2*off(2)-  off(3)) + 128*a(pos-off(2)-  off(3)) - 240*a(pos-  off(3)) + 128*a(pos+off(2)-  off(3)) - 8*a(pos+2*off(2)-  off(3)) &
            &      + 8*a(pos-2*off(2)+  off(3)) - 128*a(pos-off(2)+  off(3)) + 240*a(pos+  off(3)) - 128*a(pos+off(2)+  off(3)) + 8*a(pos+2*off(2)+  off(3)) &
            &      -   a(pos-2*off(2)+2*off(3)) +  16*a(pos-off(2)+2*off(3)) -  30*a(pos+2*off(3)) +  16*a(pos+off(2)+2*off(3)) -   a(pos+2*off(2)+2*off(3))) / (144*dx(2)**2*dx(3))
       f(1,3,2) = f(3,1,2)
       f(2,3,2) = f(3,2,2)
       f(3,3,2) = (+   a(pos-2*off(3)-2*off(2)) -  16*a(pos-off(3)-2*off(2)) +  30*a(pos-2*off(2)) -  16*a(pos+off(3)-2*off(2)) +   a(pos+2*off(3)-2*off(2)) &
            &      - 8*a(pos-2*off(3)-  off(2)) + 128*a(pos-off(3)-  off(2)) - 240*a(pos-  off(2)) + 128*a(pos+off(3)-  off(2)) - 8*a(pos+2*off(3)-  off(2)) &
            &      + 8*a(pos-2*off(3)+  off(2)) - 128*a(pos-off(3)+  off(2)) + 240*a(pos+  off(2)) - 128*a(pos+off(3)+  off(2)) + 8*a(pos+2*off(3)+  off(2)) &
            &      -   a(pos-2*off(3)+2*off(2)) +  16*a(pos-off(3)+2*off(2)) -  30*a(pos+2*off(2)) +  16*a(pos+off(3)+2*off(2)) -   a(pos+2*off(3)+2*off(2))) / (144*dx(3)**2*dx(2))
       f(:,1,3) = f(:,3,1)
       f(:,2,3) = f(:,3,2)
       f(1,3,3) = f(3,1,3)
       f(2,3,3) = f(3,2,3)
       f(3,3,3) = (a(pos-3*off(3)) - 8*a(pos-2*off(3)) + 13*a(pos-off(3)) - 13*a(pos+off(3)) + 8*a(pos+2*off(3)) - a(pos+3*off(3))) / (8*dx(3)**3)
    case default
       call CCTK_WARN (0, "Unsupported finite differencing order")
    end select
  end subroutine get_derivs3
#ifndef TGR_INCLUDED
end module derivs
#endif
