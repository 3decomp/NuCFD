! tests/tridsolver/system_11_symm.f90
!
!! Defines the test for solving tridiagonal systems of symmetric functions.
!!
!! Part of the tridsolver test suite.
!
! SPDX-License-Identifier: BSD-3-Clause

program test_system_11_symm
  !! Tests the solution of tridiagonal systems arising from compact finite differences of symmetric
  !! functions.

  use nucfd_coeffs
  use nucfd_trid_solver

  use nucfd_tests
  use tridsol_test_utils
  
  implicit none

  real, parameter :: pi = 4.0 * atan(1.0)

  call initialise_suite("System 11 symmetric")

  call solve_11_system(33)
  
  call finalise_suite()

contains

  subroutine solve_11_system(n)
    !! Build and solve a system of size n.

    integer, intent(in) :: n
    real, dimension(:), allocatable :: a, b, c, sol, ref, r

    real :: L, x, dx
    integer :: i

    logical :: passing
    
    passing = .true.

    L = 1.0
    dx = L / real(n - 1)
    
    call allocate_system(n, a, b, c, sol, ref)

    ! Set bulk coefficients.
    a(:) = alpha
    b(:) = 1.0
    c(:) = alpha

    ! Set boundary coefficients.
    ! For symmetric functions f'(-x) = -f'(x).
    a(1) = 0.0; b(1) = 1.0; c(1) = 0.0
    a(n) = 0.0; b(n) = 1.0; c(n) = 0.0

    ! Set reference solution and compute RHS.
    ! For symmetric functions f(-x) = f(x), e.g. cos(x)
    do i = 1, n
       x = real(i - 1) * dx
       ref(i) = cos((x / L) * (2 * pi)) ! Ensure a full period.
    end do
    call compute_rhs(a, b, c, ref, r)
    sol(:) = 0.0

    ! Solve the problem and compare solution vs reference.
    call solve(a, b, c, r, sol)
    passing = check_rms(sol, ref)

    call test_report("Solve 11 system", passing)
    
  end subroutine solve_11_system
  
end program test_system_11_symm
