! tests/fd-schemes/differentiation_rules/varmesh_constgrad.f90
!
!! Part of the fd-schemes test suite.
!
! SPDX-License-Identifier: BSD-3-Clause

program linear_rising_function
  !! Tests that the derivative of a linear function remains constant on a mesh with local
  !! refinement.

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

  integer :: i

  logical :: passing
  
  call initialise_suite("Linearly increasing function differentiation rules")

  n = 33
  L = 1.0

  ts = initialise_test(n, L)
  f = allocate_test_array(ts)
  x = create_mesh_array(ts)
  stencil = build_stencil(ts)

  f(1) = 0.0
  dfdx = 1.0
  do i = 2, n
     f(i) = f(i - 1) + (x(i) - x(i - 1)) * dfdx
  end do

  i = 3
  call centre_stencil(i, stencil)
  call deriv_rhs(f, stencil, x, dfdx)

  passing = .true.
  do i = 4, n - 2
     print *, i, x(i - 1) - x(i - 2), x(i) - x(i - 1), x(i + 1) - x(i), x(i + 2) - x(i + 1)
     call centre_stencil(i, stencil)
     call deriv_rhs(f, stencil, x, dgdx)
     passing = passing .and. check_scalar(dgdx, dfdx)
  end do
  call test_report("Variable mesh, constant gradient", passing)
  
  deallocate(x)
  deallocate(f)
  
  call finalise_suite()

contains

  subroutine centre_stencil(idx, stencil)
    !! Move a stencil to be centred at specified index.

    integer, intent(in) :: idx                          !! The index
    type(nucfd_index_stencil), intent(inout) :: stencil !! The stencil to be moved

    integer :: shift

    select type(indices => stencil%stencil)
    type is(integer)
       shift = idx - indices(0)
       indices(:) = indices(:) + shift
    end select
    
  end subroutine centre_stencil
 
end program linear_rising_function
