! tests/tridsolver/tridsol_test_utils_mod.f90
!
!! Defines the tridsol_test_utils module.
!!
!! Part of the tridsolver test suite.
!
! SPDX-License-Identifier: BSD-3-Clause

module tridsol_test_utils
  !! Module defining utility functions for the tridiagonal solver test suite.

  implicit none

  private

  public :: compute_rhs
  public :: allocate_system

contains

  subroutine compute_rhs(a, b, c, x, rhs)
    !! Computes the right hand side (i.e. residual) of a tridiagonal system for a given solution.

    real, dimension(:), intent(in) :: a, b, c, x
    real, dimension(:), allocatable, intent(out) :: rhs

    integer :: n
    integer :: i

    n = size(x)

    allocate(rhs(n))

    rhs(1) = b(1) * x(1) + c(1) * x(2)
    do i = 2, n - 1
       rhs(i) = a(i) * x(i - 1) + b(i) * x(i) + c(i) * x(i + 1)
    end do
    rhs(n) = a(1) * x(n - 1) + b(n) * x(n)

  end subroutine compute_rhs

  subroutine allocate_system(n, a, b, c, x, xref)
    !! Allocates a tridiagonal system.

    integer, intent(in) :: n
    real, dimension(:), allocatable, intent(out) :: a, b, c, x, xref

    allocate(a(n), b(n), c(n), x(n), xref(n))

  end subroutine allocate_system

end module tridsol_test_utils
