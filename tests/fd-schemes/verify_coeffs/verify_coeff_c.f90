! tests/fd-schemes/verify_coeffs/verify_coeff_c.f90
!
!! Part of the fd-schemes test suite.
!
! SPDX-License-Identifier: BSD-3-Clause

program verify_coeff_c
  !! Tests the computation of finite difference coefficient for non-uniform grids acting on f_{i+2}.

  use nucfd_types
  use nucfd_coeffs

  use nucfd_tests
  
  implicit none
  
  real :: n ! Mesh size
  real :: L ! Domain size
  real :: h ! Average grid spacing.

  type(nucfd_stencil_points) :: stencil ! Stencil of grid points.
  
  real :: c
  real, parameter :: cref = 1.0 / 9.0

  real :: numerator, denominator, divisor
  real, parameter :: numerator_f2ref = 2.0 / 3.0
  real, parameter :: denominator_f2ref = 6.0
  real, parameter :: divisor_f2ref = 4.0
  
  call initialise_suite("Verify coefficient C")

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

  call coeff_c_components(points_to_deltas(stencil), numerator, denominator, divisor)
  call test_report("Coefficient C numerator", check_scalar(numerator, numerator_f2ref * (h**3)))
  call test_report("Coefficient C denominator", check_scalar(denominator, denominator_f2ref * (h**3)))
  call test_report("Coefficient C divisor", check_scalar(divisor, divisor_f2ref * h))
  c = coeff_c(stencil)
  call test_report("Coefficient C", check_scalar(c, cref / (4.0 * h)))

  call finalise_suite()

end program verify_coeff_c
