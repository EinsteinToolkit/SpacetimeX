! $Header$

#include "cctk.h"

module pointwise2
!!$  use derivs2
  implicit none
  private
  public calc_2position
  public calc_2offsets
  public get_2scalar
  public get_2vector
  public get_2tensor
  public set_2scalar
  public set_2vector
  public set_2tensor
  public get_2scalarderivs
  public get_2vectorderivs
  public get_2tensorderivs
  public get_2scalarderivs2
  public get_2vectorderivs2
  public get_2tensorderivs2
  public get_2scalarderivs3
  public get_2vectorderivs3
  public get_2tensorderivs3
contains
#define TGR_INCLUDED
#include "EinsteinUtils/TGRtensor/src/derivs2.F90"
#undef TGR_INCLUDED
  subroutine calc_2position (shape, i,j, pos)
    integer, intent(in)  :: shape(2)
    integer, intent(in)  :: i,j
    integer, intent(out) :: pos
    pos = 1 + i-1 + shape(1) * (j-1)
  end subroutine calc_2position
  subroutine calc_2offsets (shape, off)
    integer, intent(in)  :: shape(2)
    integer, intent(out) :: off(2)
    off(1) = 1
    off(2) = shape(1)
  end subroutine calc_2offsets
  subroutine get_2scalar (a, f, pos)
    CCTK_REAL, intent(in)  :: a(*)
    CCTK_REAL, intent(out) :: f
    integer,   intent(in)  :: pos
    f = a(pos)
  end subroutine get_2scalar
  subroutine get_2vector (ax,ay, f, pos)
    CCTK_REAL, intent(in)  :: ax(*),ay(*)
    CCTK_REAL, intent(out) :: f(2)
    integer,   intent(in)  :: pos
    f(1) = ax(pos)
    f(2) = ay(pos)
  end subroutine get_2vector
  subroutine get_2tensor (axx,axy,ayy, f, pos)
    CCTK_REAL, intent(in)  :: axx(*),axy(*),ayy(*)
    CCTK_REAL, intent(out) :: f(2,2)
    integer,   intent(in)  :: pos
    f(1,1) = axx(pos)
    f(1,2) = axy(pos)
    f(2,1) = axy(pos)
    f(2,2) = ayy(pos)
  end subroutine get_2tensor
  subroutine set_2scalar (f, a, pos)
    CCTK_REAL, intent(in) :: f
    CCTK_REAL             :: a(*)
    integer,   intent(in) :: pos
    a(pos) = f
  end subroutine set_2scalar
  subroutine set_2vector (f, ax,ay, pos)
    CCTK_REAL, intent(in) :: f(2)
    CCTK_REAL             :: ax(*),ay(*)
    integer,   intent(in) :: pos
    ax(pos) = f(1)
    ay(pos) = f(2)
  end subroutine set_2vector
  subroutine set_2tensor (f, axx,axy,ayy, pos)
    CCTK_REAL, intent(in) :: f(2,2)
    CCTK_REAL             :: axx(*),axy(*),ayy(*)
    integer,   intent(in) :: pos
    axx(pos) = f(1,1)
    axy(pos) = f(1,2)
    ayy(pos) = f(2,2)
  end subroutine set_2tensor
  subroutine get_2scalarderivs (a, f, pos, off, dx, order)
    CCTK_REAL, intent(in)            :: a(*)
    CCTK_REAL, intent(out)           :: f(2)
    integer,   intent(in)            :: pos, off(2)
    CCTK_REAL, intent(in)            :: dx(2)
    integer,   intent(in),  optional :: order
    call get_2derivs (a, f, pos, off, dx, order)
  end subroutine get_2scalarderivs
  subroutine get_2vectorderivs (ax,ay, f, pos, off, dx, order)
    CCTK_REAL, intent(in)            :: ax(*),ay(*)
    CCTK_REAL, intent(out)           :: f(2,2)
    integer,   intent(in)            :: pos, off(2)
    CCTK_REAL, intent(in)            :: dx(2)
    integer,   intent(in),  optional :: order
    CCTK_REAL :: fx(2),fy(2)
    call get_2derivs (ax, fx, pos, off, dx, order)
    call get_2derivs (ay, fy, pos, off, dx, order)
    f(1,:) = fx
    f(2,:) = fy
  end subroutine get_2vectorderivs
  subroutine get_2tensorderivs (axx,axy,ayy, f, pos, off, dx, order)
    CCTK_REAL, intent(in)            :: axx(*),axy(*),ayy(*)
    CCTK_REAL, intent(out)           :: f(2,2,2)
    integer,   intent(in)            :: pos, off(2)
    CCTK_REAL, intent(in)            :: dx(2)
    integer,   intent(in),  optional :: order
    CCTK_REAL :: fxx(2),fxy(2),fyy(2)
    call get_2derivs (axx, fxx, pos, off, dx, order)
    call get_2derivs (axy, fxy, pos, off, dx, order)
    call get_2derivs (ayy, fyy, pos, off, dx, order)
    f(1,1,:) = fxx
    f(1,2,:) = fxy
    f(2,1,:) = fxy
    f(2,2,:) = fyy
  end subroutine get_2tensorderivs
  subroutine get_2scalarderivs2 (a, f, pos, off, dx, order)
    CCTK_REAL, intent(in)            :: a(*)
    CCTK_REAL, intent(out)           :: f(2,2)
    integer,   intent(in)            :: pos, off(2)
    CCTK_REAL, intent(in)            :: dx(2)
    integer,   intent(in),  optional :: order
    call get_2derivs2 (a, f, pos, off, dx, order)
  end subroutine get_2scalarderivs2
  subroutine get_2vectorderivs2 (ax,ay, f, pos, off, dx, order)
    CCTK_REAL, intent(in)            :: ax(*),ay(*)
    CCTK_REAL, intent(out)           :: f(2,2,2)
    integer,   intent(in)            :: pos, off(2)
    CCTK_REAL, intent(in)            :: dx(2)
    integer,   intent(in),  optional :: order
    CCTK_REAL :: fx(2,2),fy(2,2)
    call get_2derivs2 (ax, fx, pos, off, dx, order)
    call get_2derivs2 (ay, fy, pos, off, dx, order)
    f(1,:,:) = fx
    f(2,:,:) = fy
  end subroutine get_2vectorderivs2
  subroutine get_2tensorderivs2 (axx,axy,ayy, f, pos, off, dx, order)
    CCTK_REAL, intent(in)            :: axx(*),axy(*),ayy(*)
    CCTK_REAL, intent(out)           :: f(2,2,2,2)
    integer,   intent(in)            :: pos, off(2)
    CCTK_REAL, intent(in)            :: dx(2)
    integer,   intent(in),  optional :: order
    CCTK_REAL :: fxx(2,2),fxy(2,2),fyy(2,2)
    call get_2derivs2 (axx, fxx, pos, off, dx, order)
    call get_2derivs2 (axy, fxy, pos, off, dx, order)
    call get_2derivs2 (ayy, fyy, pos, off, dx, order)
    f(1,1,:,:) = fxx
    f(1,2,:,:) = fxy
    f(2,1,:,:) = fxy
    f(2,2,:,:) = fyy
  end subroutine get_2tensorderivs2
  subroutine get_2scalarderivs3 (a, f, pos, off, dx, order)
    CCTK_REAL, intent(in)            :: a(*)
    CCTK_REAL, intent(out)           :: f(2,2,2)
    integer,   intent(in)            :: pos, off(2)
    CCTK_REAL, intent(in)            :: dx(2)
    integer,   intent(in),  optional :: order
    call get_2derivs3 (a, f, pos, off, dx, order)
  end subroutine get_2scalarderivs3
  subroutine get_2vectorderivs3 (ax,ay, f, pos, off, dx, order)
    CCTK_REAL, intent(in)            :: ax(*),ay(*)
    CCTK_REAL, intent(out)           :: f(2,2,2,2)
    integer,   intent(in)            :: pos, off(2)
    CCTK_REAL, intent(in)            :: dx(2)
    integer,   intent(in),  optional :: order
    CCTK_REAL :: fx(2,2,2),fy(2,2,2)
    call get_2derivs3 (ax, fx, pos, off, dx, order)
    call get_2derivs3 (ay, fy, pos, off, dx, order)
    f(1,:,:,:) = fx
    f(2,:,:,:) = fy
  end subroutine get_2vectorderivs3
  subroutine get_2tensorderivs3 (axx,axy,ayy, f, pos, off, dx, order)
    CCTK_REAL, intent(in)            :: axx(*),axy(*),ayy(*)
    CCTK_REAL, intent(out)           :: f(2,2,2,2,2)
    integer,   intent(in)            :: pos, off(2)
    CCTK_REAL, intent(in)            :: dx(2)
    integer,   intent(in),  optional :: order
    CCTK_REAL :: fxx(2,2,2),fxy(2,2,2),fyy(2,2,2)
    call get_2derivs3 (axx, fxx, pos, off, dx, order)
    call get_2derivs3 (axy, fxy, pos, off, dx, order)
    call get_2derivs3 (ayy, fyy, pos, off, dx, order)
    f(1,1,:,:,:) = fxx
    f(1,2,:,:,:) = fxy
    f(2,1,:,:,:) = fxy
    f(2,2,:,:,:) = fyy
  end subroutine get_2tensorderivs3
end module pointwise2
