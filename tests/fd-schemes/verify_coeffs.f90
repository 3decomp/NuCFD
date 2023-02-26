program verify_coeffs
  !! Tests the computation of finite difference coefficients for non-uniform grids.
  !!
  !! Part of the fd-schemes test suite.
  !!
  !! SPDX-License-Identifier: BSD-3-Clause

  use nucfd_types
  use nucfd_coeffs

  use nucfd_tests
  
  implicit none

  type(nucfd_stencil_points) :: stencil ! Stencil of grid points.
  
  real :: n ! Mesh size
  real :: L ! Domain size
  real :: h ! Average grid spacing.
  
  real :: a, b, c, d, e
  real, parameter :: aref = 14.0 / 9.0
  real, parameter :: bref = -aref
  real, parameter :: cref = 1.0 / 9.0
  real, parameter :: dref = -cref
  real, parameter :: eref = 0.0
  
  call initialise_suite("Verify coefficients")

  n = 128
  L = 1.0
  h = L / real(n - 1)

  call create_stencil(6, 4, stencil)
  
  stencil%stencil(:) = 0.0
  stencil%stencil(-3) = -3.0 * h
  stencil%stencil(-2) = -2.0 * h
  stencil%stencil(-1) = -1.0 * h
  stencil%stencil(0) = 0.0
  stencil%stencil(1) = +1.0 * h
  stencil%stencil(2) = +2.0 * h
  
  a = coeff_a(stencil)
  call test_report("Coefficient A", check_scalar(a, aref / (2.0 * h)))
  b = coeff_b(stencil)
  call test_report("Coefficient B", check_scalar(b, bref / (2.0 * h)))
  c = coeff_c(stencil)
  call test_report("Coefficient C", check_scalar(c, cref / (4.0 * h)))
  d = coeff_d(stencil)
  call test_report("Coefficient D", check_scalar(d, dref / (4.0 * h)))
  e = coeff_e(stencil)
  call test_report("Coefficient E", check_scalar(e, eref))
  
  stencil%stencil(-3) = -5.0 * h
  stencil%stencil(2)  = +3.0 * h
  a = coeff_a(stencil)
  b = coeff_b(stencil)
  call test_report("Coefficient A /= B", .not. check_scalar(a, b))
  c = coeff_c(stencil)
  d = coeff_d(stencil)
  call test_report("Coefficient C /= D", .not. check_scalar(c, d))

  call finalise_suite()
  
end program verify_coeffs
