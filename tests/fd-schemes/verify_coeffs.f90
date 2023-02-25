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

  use nucfd_coeffs

  use nucfd_tests
  
  implicit none

  real :: n
  real :: L
  real :: h
  
  real :: a
  real, parameter :: aref = 14.0 / 9.0
  
  call initialise_suite("Verify coefficients")

  n = 128
  L = 1.0
  h = L / real(n - 1)
  
  a = coeff_a()
  call test_report("Coefficient A", check_scalar(a, aref / (2.0 * h)))
  
  call finalise_suite()
  
end program verify_coeffs
