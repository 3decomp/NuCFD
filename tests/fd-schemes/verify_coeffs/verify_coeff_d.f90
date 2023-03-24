! tests/fd-schemes/verify_coeffs/verify_coeff_d.f90
!
!! Part of the fd-schemes test suite.
!
! SPDX-License-Identifier: BSD-3-Clause

program verify_coeff_d
  !! Tests the computation of finite difference coefficient for non-uniform grids acting on f_{i-2}.

  use nucfd_types
  use nucfd_coeffs

  use nucfd_tests
  
  implicit none
  
  real :: n ! Mesh size
  real :: L ! Domain size
  real :: h ! Average grid spacing.

  type(nucfd_stencil_points) :: stencil ! Stencil of grid points.
  
  real :: d
  real, parameter :: dref = -1.0 / 9.0

  real :: numerator, denominator, divisor
  real, parameter :: numerator_f2ref =- 2.0 / 3.0
  real, parameter :: denominator_f2ref = 6.0
  real, parameter :: divisor_f2ref = 4.0
  
  call initialise_suite("Verify coefficient D")

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

  call coeff_d_components(points_to_deltas(stencil), numerator, denominator, divisor)
  call test_report("Coefficient D numerator", check_scalar(numerator, numerator_f2ref * (h**3)))
  call test_report("Coefficient D denominator", check_scalar(denominator, denominator_f2ref * (h**3)))
  call test_report("Coefficient D divisor", check_scalar(divisor, divisor_f2ref * h))
  d = coeff_d(stencil)
  call test_report("Coefficient D", check_scalar(d, dref / (4.0 * h)))

  call finalise_suite()

end program verify_coeff_d
