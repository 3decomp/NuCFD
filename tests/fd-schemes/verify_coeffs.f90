!!!! tests/fd-schemes/verify_coeffs.f90
!!!
!!!! Description
!!!
!!! Part of the fd-schemes test suite.
!!! Tests the computation of finite difference coefficients for non-uniform grids.
!!!
!!!! LICENSE
!!!
!!! SPDX-License-Identifier: BSD-3-Clause
!!!

program verify_coeffs

  use nucfd_types
  use nucfd_coeffs

  use nucfd_tests
  
  implicit none

  type(nucfd_stencil_points) :: stencil ! Stencil of grid points.
  
  real :: n ! Mesh size
  real :: L ! Domain size
  real :: h ! Average grid spacing.
  
  real :: a, b
  real, parameter :: aref = 14.0 / 9.0
  real, parameter :: bref = -aref
  
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

  stencil%stencil(-3) = -5.0 * h
  stencil%stencil(2)  = +3.0 * h
  a = coeff_a(stencil)
  b = coeff_b(stencil)
  call test_report("Coefficient A /= B", .not. check_scalar(a, b))

  call finalise_suite()
  
end program verify_coeffs
