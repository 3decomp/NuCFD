! src/nucfd_deriv_mod.f90
!
!! Defines the nucfd_deriv module.
!
! SPDX-License-Identifier: BSD-3-Clause

module nucfd_deriv
  !! Module defining the coefficients for compact finite difference schemes on non-uniform grids.

  use nucfd_types
  use nucfd_coeffs
  
  implicit none

  private
  public :: deriv_rhs

contains
  
  subroutine deriv_rhs(f, stencil, x, dfdx)
    !! Compute the RHS of the derivative of a function f based on a stencil of indices.

    real, dimension(:), intent(in) :: f
    type(nucfd_index_stencil), intent(in) :: stencil
    real, dimension(:), intent(in) :: x
    real, intent(out) :: dfdx

    type(nucfd_stencil_points) :: stencil_coordinates

    real :: a, b, c, d, e

    integer :: offset

    offset = lbound(x, 1) - (-1)
    
    call create_stencil(5, 3, stencil_coordinates)

    select type(indices => stencil%stencil)
    type is(integer)
       select type(points => stencil_coordinates%stencil)
       type is(real)
          points(:) = x(indices(:) + offset)
       class default
          print *, "Error: Coordinate stencil is misallocated!"
          error stop
       end select

       a = coeff_a(stencil_coordinates)
       b = coeff_b(stencil_coordinates)
       c = coeff_c(stencil_coordinates)
       d = coeff_d(stencil_coordinates)
       e = coeff_e(stencil_coordinates)

       dfdx = a * (f(indices(+1) + offset) + (b / a) * f(indices(-1) + offset)) &
            + c * (f(indices(+2) + offset) + (d / c) * f(indices(-2) + offset)) &
            + e * f(indices(0) + offset)
    class default
       print *, "Error: Index stencil is misallocated!"
       error stop
    end select
    
  end subroutine deriv_rhs
  
end module nucfd_deriv
