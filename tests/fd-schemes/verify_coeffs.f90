! tests/fd-schemes/verify_coeffs.f90
!
!! Part of the fd-schemes test suite.
!
! SPDX-License-Identifier: BSD-3-Clause

program verify_coeffs
  !! Tests the computation of finite difference coefficients for non-uniform grids.

  use nucfd_types
  use nucfd_coeffs

  use nucfd_tests
  
  implicit none
  
  real :: n ! Mesh size
  real :: L ! Domain size
  real :: h ! Average grid spacing.

  type(nucfd_stencil_points) :: stencil ! Stencil of grid points.
  
  real :: a, b, c, d, e
  real, parameter :: aref = 14.0 / 9.0
  real, parameter :: bref = -aref
  real, parameter :: cref = 1.0 / 9.0
  real, parameter :: dref = -cref
  real, parameter :: eref = 0.0

  real :: numerator, numerator_corr, denominator, divisor
  real, parameter :: numerator_f1ref = 14.0 / 3.0
  real, parameter :: numerator_corr_f1ref = 0.0
  real, parameter :: denominator_f1ref = 3.0
  real, parameter :: divisor_f1ref = 2.0
  real, parameter :: numerator_f2ref = 2.0 / 3.0
  real, parameter :: denominator_f2ref = 6.0
  real, parameter :: divisor_f2ref = 4.0
  
  call initialise_suite("Verify coefficients")

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

  call coeff_a_components(points_to_deltas(stencil), numerator, numerator_corr, denominator, divisor)
  call test_report("Coefficient A numerator", check_scalar(numerator, numerator_f1ref * (h**3)))
  call test_report("Coefficient A numerator correction", check_scalar(numerator_corr, numerator_corr_f1ref * (h**3)))
  call test_report("Coefficient A denominator", check_scalar(denominator, denominator_f1ref * (h**3)))
  call test_report("Coefficient A divisor", check_scalar(divisor, divisor_f1ref * h))
  a = coeff_a(stencil)
  call test_report("Coefficient A", check_scalar(a, aref / (2.0 * h)))

  call coeff_b_components(points_to_deltas(stencil), numerator, numerator_corr, denominator, divisor)
  call test_report("Coefficient B numerator", check_scalar(numerator, -numerator_f1ref * (h**3)))
  call test_report("Coefficient B numerator correction", check_scalar(numerator_corr, -numerator_corr_f1ref * (h**3)))
  call test_report("Coefficient B denominator", check_scalar(denominator, denominator_f1ref * (h**3)))
  call test_report("Coefficient B divisor", check_scalar(divisor, divisor_f1ref * h))
  b = coeff_b(stencil)
  call test_report("Coefficient B", check_scalar(b, bref / (2.0 * h)))

  call coeff_c_components(points_to_deltas(stencil), numerator, denominator, divisor)
  call test_report("Coefficient C numerator", check_scalar(numerator, numerator_f2ref * (h**3)))
  call test_report("Coefficient C denominator", check_scalar(denominator, denominator_f2ref * (h**3)))
  call test_report("Coefficient C divisor", check_scalar(divisor, divisor_f2ref * h))
  c = coeff_c(stencil)
  call test_report("Coefficient C", check_scalar(c, cref / (4.0 * h)))

  call coeff_d_components(points_to_deltas(stencil), numerator, denominator, divisor)
  call test_report("Coefficient C numerator", check_scalar(numerator, -numerator_f2ref * (h**3)))
  call test_report("Coefficient C denominator", check_scalar(denominator, denominator_f2ref * (h**3)))
  call test_report("Coefficient C divisor", check_scalar(divisor, divisor_f2ref * h))
  d = coeff_d(stencil)
  call test_report("Coefficient D", check_scalar(d, dref / (4.0 * h)))

  e = coeff_e(stencil)
  call test_report("Coefficient E", check_scalar(e, eref))

  select type(points => stencil%stencil)
  type is(real)
     points(-2) = -4.0 * h
     points(2)  = +3.0 * h
  class default
     print *, "Error: Coordinate stencil is misallocated!"
  end select
  
  a = coeff_a(stencil)
  b = coeff_b(stencil)
  call test_report("Coefficient A /= B", .not. check_scalar(a, b))
  c = coeff_c(stencil)
  d = coeff_d(stencil)
  call test_report("Coefficient C /= D", .not. check_scalar(c, d))

  call finalise_suite()
  
end program verify_coeffs
