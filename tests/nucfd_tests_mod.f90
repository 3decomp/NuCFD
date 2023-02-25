!!!! tests/tridsolver/system_11_symm.f90
!!!
!!!! Description
!!!
!!! Part of the tridsolver test suite.
!!! Tests the solution of tridiagonal systems arising from compact finite differences of symmetric
!!! functions.
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
  public :: check_rms
  
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

  logical function check_rms(test, ref)
    ! Compute RMS of error and report errors.

    real, dimension(:), intent(in) :: test ! The test data
    real, dimension(:), intent(in) :: ref  ! The reference data

    integer :: n
    integer :: i

    real :: rms
    integer :: fail_ctr
    logical :: test_passing

    n = size(ref)
    
    rms = sqrt(sum((test - ref)**2) / real(n))
    if (rms > (2 * epsilon(rms))) then
       test_passing = .false.
       
       fail_ctr = 0
       do i = 1, n
          if (fail_ctr >= 10) then
             exit
          end if

          if (abs(test(i) - ref(i)) > epsilon(rms) * ref(i)) then
             print *, i, ref(i), test(i)
             fail_ctr = fail_ctr + 1
          end if
       end do
    else
       test_passing = .true.
    end if

    check_rms = test_passing
    
  end function check_rms
  
end module nucfd_tests
