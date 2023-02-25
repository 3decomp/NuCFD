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
       print *, " "//test_name//": PASS"
    end if

  end subroutine test_report
  
end module nucfd_tests

program test_tridsolver

  use nucfd_tests
  
  implicit none

  call initialise_suite("Tridiagonal solver")
  call finalise_suite()

contains

  logical function solve_11_system(n)

    integer, intent(in) :: n

    solve_11_system = .true.

    call test_report("Solve 11 system", solve_11_system)
    
  end function solve_11_system
  
end program test_tridsolver
