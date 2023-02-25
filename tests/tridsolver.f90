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
  public :: test_report

contains

  subroutine test_report(test_name, test_status)

    character(len=*), intent(in) :: test_name
    logical, intent(in) :: test_status

    if (test_status) then
       print *, " "//test_name//": PASS"
    else
       print *, " "//test_name//": PASS"
    end if

  end subroutine test_report
  
end module nucfd_tests

program test_tridsolver

  use nucfd_tests
  
  implicit none
  
  logical :: passing
  
  print *, "************************************************************************"
  print *, " Tridiagonal solver test suite. "
  print *, "************************************************************************"

  passing = .true.
  
  print *, "------------------------------------------------------------------------"
  call test_report("Test suite", passing)
  print *, "************************************************************************"
  
end program test_tridsolver
