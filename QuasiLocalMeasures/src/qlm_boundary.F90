#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



module qlm_boundary
  use cctk
  implicit none
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  private
  public set_boundary
  
  interface set_boundary
     module procedure set_boundary_real
     module procedure set_boundary_complex
  end interface
  
contains
  
  subroutine set_boundary_real (CCTK_ARGUMENTS, hn, f, parity)
    DECLARE_CCTK_ARGUMENTS
    integer, intent(in) :: hn
    CCTK_REAL           :: f(:,:)
    integer, intent(in) :: parity
    
    CCTK_REAL, allocatable :: ff(:,:)
    
    integer :: ni, nj, gi, gj
    
    integer :: i, j, jj
    
    if (abs(parity) /= 1) then
       call CCTK_WARN (0, "Parity must be +1 or -1")
    end if
    
    
    
    ni = qlm_ntheta(hn)
    nj = qlm_nphi(hn)
    gi = qlm_nghoststheta(hn)
    gj = qlm_nghostsphi(hn)
    
    
    
    f(       :gi, :nj) = -42
    f(ni-gi+1:ni, :nj) = -42
    f(:ni,        :gj) = -42
    f(:ni, nj-gi+1:nj) = -42
    
    allocate (ff(gi, nj))
    
    ! lower theta
    ff(:, :nj) = f(1+gi:2*gi, :nj)
    
    ! polar boundary condition at north pole
    do j=1,nj/2
       do i=1,gi
          jj = j + (nj-2*gj) / 2
          f(i,j) = parity * ff(gi-i+1,jj)
       end do
    end do
    do j=nj/2+1,nj
       do i=1,gi
          jj = j - (nj-2*gj) / 2
          f(i,j) = parity * ff(gi-i+1,jj)
       end do
    end do
    
    ! upper theta
    ff(:, :nj) = f(ni-2*gi+1:ni-gi, :nj)
    
    ! polar boundary condition at south pole
    do j=1,nj/2
       do i=1,gi
          jj = j + (nj-2*gj) / 2
          f(ni-i+1,j) = parity * ff(i,jj)
       end do
    end do
    do j=nj/2+1,nj
       do i=1,gi
          jj = j - (nj-2*gj) / 2
          f(ni-i+1,j) = parity * ff(i,jj)
       end do
    end do
    
    deallocate (ff)
    
    
    
    allocate (ff(ni,gj))
    
    ! lower phi
    ff(:ni, :) = f(:ni, nj-2*gj+1:nj-gj)
    
    ! periodic boundary at null meridian
    do j=1,gj
       do i=1,ni
          f(i,j) = ff(i,j)
       end do
    end do
    
    ! upper phi
    ff(:ni, :) = f(:ni, gj+1:2*gj)
    
    ! periodic boundary at null meridian
    do j=1,gj
       do i=1,ni
          f(i,nj-j+1) = ff(i,gj-j+1)
       end do
    end do
    
    deallocate (ff)
    
    
    
    f(ni+1:, :nj) = 0
    f(:, nj+1:) = 0
    
  end subroutine set_boundary_real
  
  
  
  subroutine set_boundary_complex (CCTK_ARGUMENTS, hn, f, parity)
    DECLARE_CCTK_ARGUMENTS
    integer, intent(in) :: hn
    CCTK_COMPLEX        :: f(:,:)
    integer, intent(in) :: parity
    
    CCTK_REAL :: fre(size(f,1),size(f,2))
    CCTK_REAL :: fim(size(f,1),size(f,2))
    
    integer :: ni, nj, gi, gj
    
    ni = qlm_ntheta(hn)
    nj = qlm_nphi(hn)
    gi = qlm_nghoststheta(hn)
    gj = qlm_nghostsphi(hn)
    
    f(       :gi, :nj) = 0
    f(ni-gi+1:ni, :nj) = 0
    f(:ni,        :gj) = 0
    f(:ni, nj-gj+1:nj) = 0
    
    fre = real(f)
    fim = aimag(f)
    call set_boundary_real (CCTK_PASS_FTOF, hn, fre, parity)
    call set_boundary_real (CCTK_PASS_FTOF, hn, fim, parity)
    f = cmplx(fre, fim, kind(f))
  end subroutine set_boundary_complex
  
end module qlm_boundary
