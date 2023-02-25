!!!! tests/tridsolver.f90
!!!
!!!! Description
!!!
!!! Test suite for the tridiagonal solver.
!!!
!!!! LICENSE
!!!
!!! SPDX-License-Identifier: BSD-3-Clause
!!!

module nucfd_tests

  implicit none

  private
  public :: initialise_suite, finalise_suite
  public :: test_report

  character(len=:), allocatable :: suite_name
  logical, save :: passing
  
contains

  subroutine initialise_suite(test_suite_name)

    character(len=*), intent(in) :: test_suite_name

    suite_name = test_suite_name
    
    print *, "************************************************************************"
    print *, " "//suite_name//" test suite. "
    print *, "------------------------------------------------------------------------"

    passing = .true.
    
  end subroutine initialise_suite

  subroutine finalise_suite()
  
    print *, "------------------------------------------------------------------------"
    call test_report(suite_name)
    print *, "************************************************************************"

    if (.not. passing) then
       error stop 1
    end if
    
  end subroutine finalise_suite
  
  subroutine test_report(test_name, test_status)

    character(len=*), intent(in) :: test_name
    logical, intent(in), optional :: test_status

    logical :: report_status

    if (present(test_status)) then
       report_status = test_status
       passing = test_status .and. passing
    else
       report_status = passing
    end if
    
    if (report_status) then
       print *, " "//test_name//": PASS"
    else
       print *, " "//test_name//": FAIL"
    end if

  end subroutine test_report
  
end module nucfd_tests

program test_tridsolver

  use nucfd_tests
  use nucfd_coeffs
  use nucfd_trid_solver
  
  implicit none

  real, parameter :: pi = 4.0 * atan(1.0)

  call initialise_suite("Tridiagonal solver")

  call solve_11_system(33)
  
  call finalise_suite()

contains

  subroutine solve_11_system(n)

    integer, intent(in) :: n
    real, dimension(:), allocatable :: a, b, c, sol, ref, r

    real :: L, x, dx
    integer :: i
    real :: rms

    logical :: passing
    integer :: fail_ctr
    
    passing = .true.

    L = 1.0
    dx = L / real(n - 1)
    
    call allocate_system(n, a, b, c, sol, ref)
    a(:) = alpha
    b(:) = 1.0
    c(:) = alpha
    do i = 1, n
       x = real(i - 1) * dx
       ref(i) = cos((x / L) * (2 * pi))
    end do
    call compute_rhs(a, b, c, ref, r)
    sol(:) = 0.0

    call solve(a, b, c, r, sol)
    
    rms = sqrt(sum((sol - ref)**2) / real(n))
    if (rms > (2 * epsilon(rms))) then
       passing = .false.

       fail_ctr = 0
       do i = 1, n
          if (fail_ctr >= 10) then
             exit
          end if

          if (abs(sol(i) - ref(i)) > epsilon(rms) * ref(i)) then
             print *, i, ref(i), sol(i)
             fail_ctr = fail_ctr + 1
          end if
       end do
    end if

    call test_report("Solve 11 system", passing)
    
  end subroutine solve_11_system

  subroutine compute_rhs(a, b, c, x, rhs)

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

    integer, intent(in) :: n
    real, dimension(:), allocatable, intent(out) :: a, b, c, x, xref

    allocate(a(n), b(n), c(n), x(n), xref(n))
    
  end subroutine allocate_system
  
end program test_tridsolver
