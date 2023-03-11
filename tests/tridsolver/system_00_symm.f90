program test_system_00_symm
  !! Tests the solution of tridiagonal systems arising from compact finite differences of symmetric
  !! functions on periodic domains.
  !!
  !! Part of the tridsolver test suite.
  !!
  !! SPDX-License-Identifier: BSD-3-Clause

  use nucfd_coeffs
  use nucfd_trid_solver

  use nucfd_tests
  use tridsol_test_utils

  implicit none

  real, parameter :: pi = 4.0 * atan(1.0)

  call initialise_suite("System 00 symmetric")

  call solve_00_system(32)

  call finalise_suite()

contains

  subroutine solve_00_system(n)
    !! Build and solve a periodic system of size n.

    integer, intent(in) :: n
    real, dimension(:), allocatable :: a, b, c, sol, ref, r

    real :: L, x, dx
    integer :: i

    logical :: passing

    passing = .true.

    L = 1.0
    dx = L / real(n)

    call allocate_system(n, a, b, c, sol, ref)

    ! Set bulk coefficients.
    a(:) = alpha
    b(:) = 1.0
    c(:) = alpha

    ! Set reference solution and compute RHS.
    ! For symmetric functions f(-x) = f(x), e.g. cos(x)
    do i = 1, n
       x = real(i - 1) * dx
       ref(i) = cos((x / L) * (2 * pi)) ! Ensure a full period.
    end do
    call compute_rhs(a, b, c, ref, r)
    r(1) = r(1) + a(1) * ref(n)
    r(n) = r(n) + c(n) * ref(1)
    sol(:) = 0.0

    ! Solve the problem and compare solution vs reference.
    call solve_cyclic(a, b, c, r, sol)
    passing = check_rms(sol, ref)

    call test_report("Solve 00 system", passing)

  end subroutine solve_00_system

end program test_system_00_symm
