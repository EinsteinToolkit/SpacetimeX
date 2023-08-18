#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_output_vtk (CCTK_ARGUMENTS, hn, file_name)
  use cctk
  use constants
  use qlm_variables
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  integer,      intent(in) :: hn
  character(*), intent(in) :: file_name

  integer, parameter :: unit = 194

  integer   :: nptheta, npphi
  integer   :: i, j
  CCTK_REAL :: xx, yy, zz

  nptheta = qlm_ntheta(hn) - 2*qlm_nghoststheta(hn) 
  npphi = qlm_nphi(hn) - 2*qlm_nghostsphi(hn)

  open (unit=unit, file=file_name, action='write')

  write (unit, '(A)') '# vtk DataFile Version 2.0' 
  write (unit, '(A)') 'Horizon data' 
  write (unit, '(A)') 'ASCII' 
  write (unit, '(A)') 'DATASET POLYDATA'
  write (unit, '(A,X,I10,X,A)') 'POINTS', nptheta*npphi, 'float'

  do i = 1+qlm_nghoststheta(hn), qlm_ntheta(hn)-qlm_nghoststheta(hn)
     do j = 1+qlm_nghostsphi(hn), qlm_nphi(hn)-qlm_nghostsphi(hn)
        xx = qlm_x(i, j, hn) 
        yy = qlm_y(i, j, hn) 
        zz = qlm_z(i, j, hn) 
        write (unit, *) xx, yy, zz 
     end do
  end do

  write (unit, '()')
  write (unit, '(A,X,I10,X,I10)') &
       'POLYGONS', npphi*(nptheta-1), 5*npphi*(nptheta-1)

  do i = 0, nptheta-2
     do j = 0, npphi-2
        write (unit,'(I1,4(I10))') &
             4, i*npphi+j, (i+1)*npphi+j, (i+1)*npphi+j+1, i*npphi+j+1
     end do
     write (unit,'(I1,4(I10))') &
          4, i*npphi+npphi-1, (i+1)*npphi+npphi-1, (i+1)*npphi, i*npphi
  end do

  write (unit, '()')
  write (unit, '(A,X,I10)') 'POINT_DATA', nptheta*npphi

  call writescalar ('shape', qlm_shape(:,:,hn))
  call writescalar ('l0', qlm_l0(:,:,hn))
  call writescalar ('l1', qlm_l1(:,:,hn))
  call writescalar ('l2', qlm_l2(:,:,hn))
  call writescalar ('l3', qlm_l3(:,:,hn))
  call writescalar ('n0', qlm_n0(:,:,hn))
  call writescalar ('n1', qlm_n1(:,:,hn))
  call writescalar ('n2', qlm_n2(:,:,hn))
  call writescalar ('n3', qlm_n3(:,:,hn))
  call writescalar_complex ('m0', qlm_m0(:,:,hn))
  call writescalar_complex ('m1', qlm_m1(:,:,hn))
  call writescalar_complex ('m2', qlm_m2(:,:,hn))
  call writescalar_complex ('m3', qlm_m3(:,:,hn))
  call writescalar_complex ('npkappa', qlm_npkappa(:,:,hn))
  call writescalar_complex ('nptau', qlm_nptau(:,:,hn))
  call writescalar_complex ('npsigma', qlm_npsigma(:,:,hn))
  call writescalar_complex ('nprho', qlm_nprho(:,:,hn))
  call writescalar_complex ('npepsilon', qlm_npepsilon(:,:,hn))
  call writescalar_complex ('npgamma', qlm_npgamma(:,:,hn))
  call writescalar_complex ('npbeta', qlm_npbeta(:,:,hn))
  call writescalar_complex ('npalpha', qlm_npalpha(:,:,hn))
  call writescalar_complex ('nppi', qlm_nppi(:,:,hn))
  call writescalar_complex ('npnu', qlm_npnu(:,:,hn))
  call writescalar_complex ('npmu', qlm_npmu(:,:,hn))
  call writescalar_complex ('nplambda', qlm_nplambda(:,:,hn))
  call writescalar_complex ('psi0', qlm_psi0(:,:,hn))
  call writescalar_complex ('psi1', qlm_psi1(:,:,hn))
  call writescalar_complex ('psi2', qlm_psi2(:,:,hn))
  call writescalar_complex ('psi3', qlm_psi3(:,:,hn))
  call writescalar_complex ('psi4', qlm_psi4(:,:,hn))
  call writescalar ('xit', qlm_xi_t(:,:,hn))
  call writescalar ('xip', qlm_xi_p(:,:,hn))
  call writescalar ('chi', qlm_chi(:,:,hn))

  close (unit)

contains

  subroutine writescalar (array_name, array)
    character(*), intent(in) :: array_name
    CCTK_REAL,    intent(in) :: array(:, :)

    integer :: i, j

    write (unit, '(/A,X,A,X,A)') 'SCALARS', array_name, 'float 1'
    write (unit, '(A)') 'LOOKUP_TABLE default'
    do i = 1+qlm_nghoststheta(hn), qlm_ntheta(hn)-qlm_nghoststheta(hn)
       do j = 1+qlm_nghostsphi(hn), qlm_nphi(hn)-qlm_nghostsphi(hn)
          write (unit, *) array(i, j)
       end do
    end do
  end subroutine writescalar

  subroutine writescalar_complex (array_name, array)
    character(*), intent(in) :: array_name
    CCTK_COMPLEX, intent(in) :: array(:, :)

    call writescalar ('re' // array_name, real(array))
    call writescalar ('im' // array_name, aimag(array))
  end subroutine writescalar_complex

end subroutine qlm_output_vtk
