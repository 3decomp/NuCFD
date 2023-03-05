! tests/fd-schemes/verify_coeffs/verify_coeff_b.f90
!
!! Part of the fd-schemes test suite.
!
! SPDX-License-Identifier: BSD-3-Clause

program verify_coeff_b
  !! Tests the computation of finite difference coefficient for non-uniform grids acting on f_{i-1}.

  use nucfd_types
  use nucfd_coeffs

  use nucfd_tests
  
  implicit none
  
  real :: n ! Mesh size
  real :: L ! Domain size
  real :: h ! Average grid spacing.

  type(nucfd_stencil_points) :: stencil ! Stencil of grid points.
  
  real :: b
  real, parameter :: bref = -14.0 / 9.0

  real :: numerator, numerator_corr, denominator, divisor
  real, parameter :: numerator_f1ref = -14.0 / 3.0
  real, parameter :: numerator_corr_f1ref = -0.0
  real, parameter :: denominator_f1ref = 3.0
  real, parameter :: divisor_f1ref = 2.0
  
  call initialise_suite("Verify coefficient B")

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

  call coeff_b_components(points_to_deltas(stencil), numerator, numerator_corr, denominator, divisor)
  call test_report("Coefficient B numerator", check_scalar(numerator, numerator_f1ref * (h**3)))
  call test_report("Coefficient B numerator correction", check_scalar(numerator_corr, numerator_corr_f1ref * (h**3)))
  call test_report("Coefficient B denominator", check_scalar(denominator, denominator_f1ref * (h**3)))
  call test_report("Coefficient B divisor", check_scalar(divisor, divisor_f1ref * h))
  b = coeff_b(stencil)
  call test_report("Coefficient B", check_scalar(b, bref / (2.0 * h)))

  call finalise_suite()

end program verify_coeff_b
