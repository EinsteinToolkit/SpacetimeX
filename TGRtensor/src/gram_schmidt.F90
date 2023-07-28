! $Header$

#include "cctk.h"

module gram_schmidt
  implicit none
  private
  public gram_schmidt_project4
  public gram_schmidt_normalise4
  
contains
  
  subroutine gram_schmidt_project4 (g, y, dy, ddy, y2, x, dx, ddx)
    CCTK_REAL, intent(in)    :: g(0:3,0:3)
    CCTK_REAL, intent(in)    :: y(0:3), dy(0:3,0:3), ddy(0:3,0:3,0:3), y2
    CCTK_REAL, intent(inout) :: x(0:3), dx(0:3,0:3), ddx(0:3,0:3,0:3)
    CCTK_REAL                :: z(0:3), dz(0:3,0:3), ddz(0:3,0:3,0:3)
    
    integer :: a, b, c, d, e
    
    do a=0,3
       z(a) = x(a)
       do d=0,3
          do e=0,3
             z(a) = z(a) - g(d,e) * x(d) * y(e) * y(a) / y2
          end do
       end do
    end do
    
    do a=0,3
       do b=0,3
          dz(a,b) = dx(a,b)
          do d=0,3
             do e=0,3
                dz(a,b) = dz(a,b) - g(d,e) * ((dx(d,b) * y(e) + x(d) * dy(e,b)) * y(a) + x(d) * y(e) * dy(a,b)) / y2
             end do
          end do
       end do
    end do
    
    do a=0,3
       do b=0,3
          do c=0,3
             ddz(a,b,c) = ddx(a,b,c)
             do d=0,3
                do e=0,3
                   ddz(a,b,c) = ddz(a,b,c) - g(d,e) * (  (ddx(d,b,c) * y(e) + dx(d,b) * dy(e,c) + dx(d,c) * dy(e,b) + x(d) * ddy(e,b,c)) * y(a) + (dx(d,b) * y(e) + x(d) * dy(e,b)) * dy(a,c) &
                        &                              + dx(d,c) * y(e) * dy(a,b) + x(d) * dy(e,c) * dy(a,b) + x(d) * y(e) * ddy(a,b,c)) / y2
                end do
             end do
          end do
       end do
    end do
    
    x = z
    dx = dz
    ddx = ddz
  end subroutine gram_schmidt_project4
  
  subroutine gram_schmidt_normalise4 (g, x, dx, ddx, x2)
    CCTK_REAL, intent(in)    :: g(0:3,0:3)
    CCTK_REAL, intent(inout) :: x(0:3), dx(0:3,0:3), ddx(0:3,0:3,0:3)
    CCTK_REAL, intent(in)    :: x2
    CCTK_REAL, parameter     :: two=2, half=1/two
    CCTK_REAL                :: z2, dz2(0:3), ddz2(0:3,0:3)
    CCTK_REAL                :: z(0:3), dz(0:3,0:3), ddz(0:3,0:3,0:3)
    
    integer :: a, b, c, d, e
    
    z2 = 0
    do d=0,3
       do e=0,3
          z2 = z2 + g(d,e) * x(d) * x(e)
       end do
    end do
    
    do a=0,3
       dz2(a) = 0
       do d=0,3
          do e=0,3
             dz2(a) = dz2(a) + 2 * g(d,e) * dx(d,a) * x(e)
          end do
       end do
    end do
    
    do a=0,3
       do b=0,3
          ddz2(a,b) = 0
          do d=0,3
             do e=0,3
                ddz2(a,b) = ddz2(a,b) + 2 * g(d,e) * ddx(d,a,b) * x(e) + 2 * g(d,e) * dx(d,a) * dx(e,b)
             end do
          end do
       end do
    end do
    
    do a=0,3
       z(a) = x(a) / sqrt(z2 / x2)
    end do
    
    do a=0,3
       do b=0,3
          dz(a,b) = (dx(a,b) * z2 - half * x(a) * dz2(b) / x2) / sqrt(z2 / x2)**3
       end do
    end do
    
    do a=0,3
       do b=0,3
          do c=0,3
             ddz(a,b,c) = ((ddx(a,b,c) * z2 + dx(a,b) * dz2(c) / x2 - half * (dx(a,c) * dz2(b) / x2 + x(a) * ddz2(b,c) / x2**2)) * z2 &
                  &        - 3*half * (dx(a,b) * z2 - half * x(a) * dz2(b) / x2) * dz2(c) / x2) / sqrt(z2 / x2)**5
          end do
       end do
    end do
    
    x = z
    dx = dz
    ddx = ddz
  end subroutine gram_schmidt_normalise4
  
end module gram_schmidt
