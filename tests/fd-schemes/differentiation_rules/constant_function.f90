! tests/fd-schemes/differentiation_rules/constant_function.f90
!
!! Part of the fd-schemes test suite.
!
! SPDX-License-Identifier: BSD-3-Clause

program constant_function
  !! Tests that rules of differentiation are respected for constant functions.

  use nucfd_types
  use nucfd_deriv
  
  use nucfd_tests
  use diff_rules_utils
  
  implicit none

  integer :: n ! Grid size
  real :: L    ! Domain size

  real, dimension(:), allocatable :: x ! Grid
  real, dimension(:), allocatable :: f ! Function
  real :: dfdx ! Derivative
  real :: dgdx ! Derivative

  type(nucfd_index_stencil) :: stencil
  type(test_setup) :: ts
  
  call initialise_suite("Constant function differentiation rules")

  n = 33
  L = 1.0

  ts = initialise_test(n, L)
  f = allocate_test_array(ts)
  x = create_mesh_array(ts)
  stencil = build_stencil(ts)

  f(:) = 1.0
  call deriv_rhs(f, stencil, x, dfdx)
  call deriv_rhs(f + 1.0, stencil, x, dgdx)
  call test_report("Shifted derivative +, const f", check_scalar(dgdx, dfdx))
  call deriv_rhs(f - 1.0, stencil, x, dgdx)
  call test_report("Shifted derivative -, const f", check_scalar(dgdx, dfdx))
  call deriv_rhs(f, stencil, x, dgdx)
  call test_report("Shifted derivative 0, const f", check_scalar(dgdx, dfdx))
  call deriv_rhs(2.0 * f, stencil, x, dgdx)
  call test_report("Scaled derivative 2x, const f", check_scalar(dgdx, 2.0 * dfdx))
  call deriv_rhs(-f, stencil, x, dgdx)
  call test_report("Scaled derivative -1, const f", check_scalar(dgdx, -dfdx))
  
  deallocate(x)
  deallocate(f)
  
  call finalise_suite()

end program constant_function
