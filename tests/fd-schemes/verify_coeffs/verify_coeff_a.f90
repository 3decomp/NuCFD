! tests/fd-schemes/verify_coeffs/verify_coeff_a.f90
!
!! Part of the fd-schemes test suite.
!
! SPDX-License-Identifier: BSD-3-Clause

program verify_coeff_a
  !! Tests the computation of finite difference coefficient for non-uniform grids acting on f_{i+1}.

  use nucfd_types
  use nucfd_coeffs

  use nucfd_tests
  
  implicit none
  
  real :: n ! Mesh size
  real :: L ! Domain size
  real :: h ! Average grid spacing.

  type(nucfd_stencil_points) :: stencil ! Stencil of grid points.
  
  real :: a

  integer :: scale
  real :: aref
  real :: numerator_f1ref, numerator_corr_f1ref, denominator_f1ref, divisor_f1ref
  
  call initialise_suite("Verify coefficient A")

  n = 128
  L = 1.0
  h = L / real(n - 1)

  call create_stencil(5, 3, stencil)

  scale = 1
  call setup_stencil(scale, scale, scale, scale)
  call compute_uniform_ref_coeffs(scale)
  call test_coefficient("Uniform", stencil, aref, numerator_f1ref, numerator_corr_f1ref, &
       denominator_f1ref, divisor_f1ref)

  !! Check coefficient is independent of scale
  scale = 2
  call setup_stencil(scale, scale, scale, scale)
  call compute_uniform_ref_coeffs(scale)
  call test_coefficient("Uniform (2h)", stencil, aref, numerator_f1ref, numerator_corr_f1ref, &
       denominator_f1ref, divisor_f1ref)
  
  call finalise_suite()

contains

  subroutine setup_stencil(sm2, sm1, sp1, sp2)

    integer, intent(in) :: sm2
    integer, intent(in) :: sm1
    integer, intent(in) :: sp1
    integer, intent(in) :: sp2

    select type(points => stencil%stencil)
    type is(real)
       points(:) = 0.0
       points(-2) = -(sm1 + sm2) * h
       points(-1) = -sm1 * h
       points(0) = 0.0
       points(1) = +sp1 * h
       points(2) = +(sp1 + sp2) * h
    class default
       print *, "Error: Coordinate stencil is misallocated!"
    end select
    
  end subroutine setup_stencil

  subroutine compute_uniform_ref_coeffs(scale)

    integer, intent(in) :: scale
    
    aref = 14.0 / 9.0 / (2.0 * (scale * h))
    numerator_f1ref = 14.0 / 3.0 * ((scale * h)**3)
    numerator_corr_f1ref = 0.0 * ((scale * h)**3)
    denominator_f1ref = 3.0 * ((scale * h)**3)
    divisor_f1ref = 2.0 * (scale * h)
  end subroutine compute_uniform_ref_coeffs
  
  subroutine test_coefficient(grid_desc, stencil, aref, numerator_f1ref, numerator_corr_f1ref, &
       denominator_f1ref, divisor_f1ref)

    character(len=*), intent(in) :: grid_desc
    type(nucfd_stencil_points), intent(in) :: stencil
    real, intent(in) :: aref
    real, intent(in) :: numerator_f1ref
    real, intent(in) :: numerator_corr_f1ref
    real, intent(in) :: denominator_f1ref
    real, intent(in) :: divisor_f1ref

    character(len=:), allocatable :: prefix
    real :: numerator, numerator_corr, denominator, divisor

    prefix = "("//grid_desc//") "
    
    call coeff_a_components(points_to_deltas(stencil), numerator, numerator_corr, &
         denominator, divisor)
    call test_report(prefix//"Coefficient A numerator", &
         check_scalar(numerator, numerator_f1ref))
    call test_report(prefix//"Coefficient A numerator correction", &
         check_scalar(numerator_corr, numerator_corr_f1ref))
    call test_report(prefix//"Coefficient A denominator", &
         check_scalar(denominator, denominator_f1ref))
    call test_report(prefix//"Coefficient A divisor", &
         check_scalar(divisor, divisor_f1ref))
    call test_report(prefix//"Coefficient A den * div", &
         check_scalar(denominator * divisor, &
         (denominator_f1ref * divisor_f1ref)))
    
    a = coeff_a(stencil)
    call test_report("Coefficient A", &
         check_scalar(a, aref))

  end subroutine test_coefficient
  
end program verify_coeff_a
