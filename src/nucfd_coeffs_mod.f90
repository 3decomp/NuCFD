!!!! src/nucfd_coeffs_mod.f90
!!!
!!!! Description
!!!
!!! Defines the coefficients for compact finite difference schemes on non-uniform grids.
!!!
!!! Provides the nucfd_coeffs module.
!!!
!!!! LICENSE
!!!
!!! SPDX-License-Identifier: BSD-3-Clause
!!!

module nucfd_coeffs
  !! Module defining the coefficients for compact finite differenceschemes on non-uniform grids.

  use nucfd_types
  
  implicit none

  private
  public :: coeff_a

  interface coeff_a
     module procedure coeff_a_points
     module procedure coeff_a_deltas
  end interface coeff_a
  
  real, parameter, public :: alpha = 1.0 / 3.0 !! Off-diagonal coefficient for first derivative
                                               !! system.
  
contains

  pure real function coeff_a_points(x)
    !! Compute the coefficient acting on f_{i+1} of the finite diference given a stencil of points.

    type(nucfd_stencil_points(6, 4)), intent(in) :: x !! Stencil of points for the finite
                                                      !! difference.

    type(nucfd_stencil_deltas(5, 3)) :: h !! Stencil of grid spacings for the finite difference.

    h = points_to_deltas(x)
    coeff_a_points = coeff_a(h)
    
  end function coeff_a_points

  pure real function coeff_a_deltas(h)
    !! Compute the coefficient acting on f_{i+1} of the finite diference given a stencil of grid
    !! spacings.

    type(nucfd_stencil_deltas(5, 3)), intent(in) :: h !! Stencil of grid spacings for the finite difference.

    coeff_a_deltas = 0.0
    
  end function coeff_a_deltas
  
end module nucfd_coeffs
