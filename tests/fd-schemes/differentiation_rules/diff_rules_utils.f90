! tests/fd-schemes/differentiation_rules/diff_rules_utils.f90
!
!! Part of the fd-schemes test suite.
!
! SPDX-License-Identifier: BSD-3-Clause

module diff_rules_utils
  !! Utility module for differentiation rules test subsuite.

  use nucfd_types
  
  implicit none

  type test_setup
     integer :: n
     real :: L
     real :: h
  end type test_setup
  
  private
  public :: test_setup
  public :: initialise_test
  public :: allocate_test_array
  public :: create_mesh_array
  public :: build_stencil
  
  integer, parameter :: width = 5
  integer, parameter :: centre = 3
  
contains

  type(test_setup) function initialise_test(n, L)

    integer, intent(in) :: n
    real, intent(in) :: L

    real :: h

    h = L / real(n - 1)
    initialise_test = test_setup(n, L, h)

  end function initialise_test

  function allocate_test_array(ts) result(arr)

    type(test_setup), intent(in) :: ts
    real, dimension(:), allocatable :: arr

    allocate(arr(ts%n))
    
  end function allocate_test_array

  function create_mesh_array(ts) result(x)

    type(test_setup), intent(in) :: ts
    real, dimension(:), allocatable :: x

    integer :: i
    
    x = allocate_test_array(ts)

    associate(n => ts%n, h => ts%h)
      do i = 1, n
         x(i) = real(i - 1) * h
      end do
    end associate

  end function create_mesh_array

  type(nucfd_index_stencil) function build_stencil(ts)

    type(test_setup), intent(in) :: ts

    integer :: i
    
    print *, "+++ Initialising stencil +++"
    call create_stencil(width, centre, build_stencil)
    i = ts%n / 2
    select type(indices => build_stencil%stencil)
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

  end function build_stencil
  
end module diff_rules_utils
