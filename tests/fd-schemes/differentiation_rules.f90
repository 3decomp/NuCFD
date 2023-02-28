! tests/fd-schemes/differentiation_rules.f90
!
!! Part of the fd-schemes test suite.
!
! SPDX-License-Identifier: BSD-3-Clause

program differentiation_rules
  !! Tests that rules of differentiation are respected.

  use nucfd_types
  use nucfd_deriv
  
  use nucfd_tests
  
  implicit none

  integer :: n ! Grid size
  real :: L    ! Domain size
  real :: h    ! Grid spacing

  real, dimension(:), allocatable :: x ! Grid
  real, dimension(:), allocatable :: f ! Function
  real :: dfdx ! Derivative
  real :: dgdx ! Derivative

  ! Stencil
  type(nucfd_index_stencil) :: stencil
  integer :: i
  integer, parameter :: width = 5
  integer, parameter :: centre = 3
  
  call initialise_suite("Differentiation rules")

  n = 33
  L = 1.0
  h = L / real(n - 1)

  allocate(x(n))
  allocate(f(n))

  print *, "+++ Initialising stencil +++"
  call create_stencil(width, centre, stencil)
  i = 33 / 2
  select type(indices => stencil%stencil)
  type is(integer)
     indices(-2) = i - 2
     indices(-1) = i - 1
     indices(+0) = i + 0
     indices(+1) = i + 1
     indices(+2) = i + 2
  class default
     print *, "Error: Index stencil is misallocated!"
     error stop
  end select

  print *, "+++ Testing derivative of constant function +++"
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

  print *, "+++ Testing derivative of linearly increasing function +++"
  f(1) = 0.0
  dfdx = 1.0
  do i = 2, n
     f(i) = f(i - 1) + h * dfdx
  end do
  call deriv_rhs(f, stencil, x, dfdx)
  call deriv_rhs(f + 1.0, stencil, x, dgdx)
  call test_report("Shifted derivative +, linear+ f", check_scalar(dgdx, dfdx))
  call deriv_rhs(f - 1.0, stencil, x, dgdx)
  call test_report("Shifted derivative -, linear+ f", check_scalar(dgdx, dfdx))
  call deriv_rhs(f, stencil, x, dgdx)
  call test_report("Shifted derivative 0, linear+ f", check_scalar(dgdx, dfdx))
  call deriv_rhs(2.0 * f, stencil, x, dgdx)
  call test_report("Scaled derivative 2x, linear+ f", check_scalar(dgdx, 2.0 * dfdx))
  call deriv_rhs(-f, stencil, x, dgdx)
  call test_report("Scaled derivative -1, linear+ f", check_scalar(dgdx, -dfdx))

  print *, "+++ Testing derivative of linearly decreasing function +++"
  f(:) = -f(:)
  call deriv_rhs(f, stencil, x, dfdx)
  call deriv_rhs(f + 1.0, stencil, x, dgdx)
  call test_report("Shifted derivative +, linear- f", check_scalar(dgdx, dfdx))
  call deriv_rhs(f - 1.0, stencil, x, dgdx)
  call test_report("Shifted derivative -, linear- f", check_scalar(dgdx, dfdx))
  call deriv_rhs(f, stencil, x, dgdx)
  call test_report("Shifted derivative 0, linear- f", check_scalar(dgdx, dfdx))
  call deriv_rhs(2.0 * f, stencil, x, dgdx)
  call test_report("Scaled derivative 2x, linear- f", check_scalar(dgdx, 2.0 * dfdx))
  call deriv_rhs(-f, stencil, x, dgdx)
  call test_report("Scaled derivative -1, linear- f", check_scalar(dgdx, -dfdx))
  
  deallocate(x)
  deallocate(f)
  
  call finalise_suite()
  
end program differentiation_rules
