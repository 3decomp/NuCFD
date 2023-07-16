! tests/fd-schemes/verify_coeffs/verify_coeff_pairsums.f90
!
!! Part of the fd-schemes test suite.
!
! SPDX-License-Identifier: BSD-3-Clause

program verify_coeff_pairsums
  !! Tests the computation of finite difference coefficients for the pairs f_{i+/-1} and f_{i+/-2}.

  use nucfd_types
  use nucfd_coeffs

  use nucfd_tests
  
  implicit none
  
  real :: n ! Mesh size
  real :: L ! Domain size
  real :: h ! Average grid spacing.

  type(nucfd_stencil_points) :: stencil ! Stencil of grid points.
  
  real :: a, b, c, d
  real, parameter :: abref = 0.0, cdref = 0.0
  
  call initialise_suite("Verify coefficient pairsums")

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

  a = coeff_a(stencil)
  b = coeff_b(stencil)
  call test_report("Coefficient A+B", check_scalar(a + b, abref))

  c = coeff_c(stencil)
  d = coeff_d(stencil)
  call test_report("Coefficient A+B", check_scalar(c + d, cdref))

  call finalise_suite()

end program verify_coeff_pairsums
