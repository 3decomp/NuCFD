! tests/fd-schemes/verify_coeffs/verify_coeff_e.f90
!
!! Part of the fd-schemes test suite.
!
! SPDX-License-Identifier: BSD-3-Clause

program verify_coeff_e
  !! Tests the computation of finite difference coefficient for non-uniform grids acting on f_{i}.

  use nucfd_types
  use nucfd_coeffs

  use nucfd_tests
  
  implicit none
  
  real :: n ! Mesh size
  real :: L ! Domain size
  real :: h ! Average grid spacing.

  type(nucfd_stencil_points) :: stencil ! Stencil of grid points.
  
  real :: e
  real, parameter :: eref = 0.0
  
  call initialise_suite("Verify coefficient E")

  n = 128
  L = 1.0
  h = L / real(n - 1)

  call create_stencil(5, 3, stencil)
  select type(points => stencil%stencil)
  type is(real)
     points(:) = 0.0
     points(-2) = -2.0 * h
     points(-1) = -1.0 * h
     points(0) = 0.0
     points(1) = +1.0 * h
     points(2) = +2.0 * h
  class default
     print *, "Error: Coordinate stencil is misallocated!"
  end select

  e = coeff_e(stencil)
  call test_report("Coefficient E", check_scalar(e, eref / (4.0 * h)))

  call finalise_suite()

end program verify_coeff_e
