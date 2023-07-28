#include "cctk.h"

module qlm_gram_schmidt
  implicit none
  private
  public gram_schmidt_project
  public gram_schmidt_normalise
  
contains
  
  subroutine gram_schmidt_project (g, dg, y, dy, y2, x, dx)
    CCTK_REAL, intent(in)    :: g(0:3,0:3), dg(0:3,0:3,0:3)
    CCTK_REAL, intent(in)    :: y(0:3), dy(0:3,0:3), y2
    CCTK_REAL, intent(inout) :: x(0:3), dx(0:3,0:3)
    CCTK_REAL                :: z(0:3), dz(0:3,0:3)
    
    integer :: a, b, d, e
    
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
                dz(a,b) = dz(a,b) &
                     - dg(d,e,b) * x(d) * y(e) * y(a) / y2 &
                     - g(d,e) * dx(d,b) * y(e) * y(a) / y2 &
                     - g(d,e) * x(d) * dy(e,b) * y(a) / y2 &
                     - g(d,e) * x(d) * y(e) * dy(a,b) / y2
             end do
          end do
       end do
    end do
    
    x = z
    dx = dz
  end subroutine gram_schmidt_project
  
  subroutine gram_schmidt_normalise (g, dg, x, dx, x2)
    CCTK_REAL, intent(in)    :: g(0:3,0:3), dg(0:3,0:3,0:3)
    CCTK_REAL, intent(inout) :: x(0:3), dx(0:3,0:3)
    CCTK_REAL, intent(in)    :: x2
    CCTK_REAL, parameter     :: two=2, half=1/two
    CCTK_REAL                :: z2, dz2(0:3)
    CCTK_REAL                :: z(0:3), dz(0:3,0:3)
    
    integer :: a, b, d, e
    
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
             dz2(a) = dz2(a) &
                  + dg(d,e,a) * x(d) * x(e) &
                  + g(d,e) * dx(d,a) * x(e) &
                  + g(d,e) * x(d) * dx(e,a)
          end do
       end do
    end do
    
    do a=0,3
       z(a) = x(a) / sqrt(z2 / x2)
    end do
    
    do a=0,3
       do b=0,3
          dz(a,b) = (+ dx(a,b) * z2 / x2 &
               &     - half * x(a) * dz2(b) / x2) / sqrt(z2 / x2)**3
       end do
    end do
    
    x = z
    dx = dz
  end subroutine gram_schmidt_normalise
  
end module qlm_gram_schmidt
