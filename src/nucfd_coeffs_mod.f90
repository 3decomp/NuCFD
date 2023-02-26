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
  public :: coeff_b
  public :: coeff_c

  interface coeff_a
     module procedure coeff_a_points
     module procedure coeff_a_deltas
  end interface coeff_a

  interface coeff_b
     module procedure coeff_b_points
     module procedure coeff_b_deltas
  end interface coeff_b

  interface coeff_c
     module procedure coeff_c_points
     module procedure coeff_c_deltas
  end interface coeff_c
  
  real, parameter, public :: alpha = 1.0 / 3.0 !! Off-diagonal coefficient for first derivative
                                               !! system.
  
contains

  real function coeff_a_points(x)
    !! Compute the coefficient acting on f_{i+1} of the finite diference given a stencil of points.

    type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                !! difference.

    type(nucfd_stencil_deltas) :: h  ! Stencil of grid spacings for the finite difference.
    
    h = points_to_deltas(x)
    coeff_a_points = coeff_a(h)
    
  end function coeff_a_points

  pure real function coeff_a_deltas(h)
    !! Compute the coefficient acting on f_{i+1} of the finite diference given a stencil of grid
    !! spacings.

    type(nucfd_stencil_deltas), intent(in) :: h !! Stencil of grid spacings for the finite difference.

    real :: hm2, hm1, h0, hp1, hp2 ! Grid deltas at i -2, -1, 0, +1, +2

    hm2 = h%stencil(-2)
    hm1 = h%stencil(-1)
    h0 = h%stencil(0)
    hp1 = h%stencil(1)
    hp2 = h%stencil(2)

    associate(beta => alpha) ! To match Gamet et al. (1999)
      coeff_a_deltas = hm1 * h0 * hp1 + h0**2 * hp1 + hm1 * h0 * hp2 + h0**2 * hp2 &
           - hm1 * h0**2 * alpha - hm1 * h0 * hp1 * alpha &
           - hm1 * h0 * hp2 * alpha - hm1 * h0 * hp1 * beta - h0**2 * hp1 * beta &
           - hm1 * hp1**2 * beta - 2.0 * h0 * hp1**2 * beta - hp1**3 * beta &
           + hm1 * h0 * hp2 * beta + h0**2 * hp2 * beta + 2.0 * hm1 * hp1 * hp2 * beta &
           + 4.0 * h0 * hp1 * hp2 * beta + 3.0 * hp1**2 * hp2 * beta
      coeff_a_deltas = coeff_a_deltas &
           / (hp1 * (h0 + hp1) * (hm1 + h0 + hp1) * hp2)
    end associate
  end function coeff_a_deltas

  real function coeff_b_points(x)
    !! Compute the coefficient acting on f_{i-1} of the finite diference given a stencil of points.

    type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                !! difference.

    type(nucfd_stencil_deltas) :: h  ! Stencil of grid spacings for the finite difference.
    
    h = points_to_deltas(x)
    coeff_b_points = coeff_b_deltas(h)
    
  end function coeff_b_points

  pure real function coeff_b_deltas(h)
    !! Compute the coefficient acting on f_{i-1} of the finite diference given a stencil of grid
    !! spacings.

    type(nucfd_stencil_deltas), intent(in) :: h !! Stencil of grid spacings for the finite difference.

    real :: hm2, hm1, h0, hp1, hp2 ! Grid deltas at i -2, -1, 0, +1, +2

    hm2 = h%stencil(-2)
    hm1 = h%stencil(-1)
    h0 = h%stencil(0)
    hp1 = h%stencil(1)
    hp2 = h%stencil(2)

    associate(beta => alpha) ! To match Gamet et al. (1999)
      coeff_b_deltas = -hm1 * hp1**2 - h0 * hp1**2 - hm1 * hp1 * hp2 - h0 * hp1 * hp2 &
           - 3.0 * hm1 * h0**2 * alpha - 4.0 * hm1 * h0 * hp1 * alpha & ! End line 1
           + h0**3 * alpha + 2.0 * h0**2 * hp1 * alpha - hm1 * hp1**2 * alpha &
           + h0 * hp1**2 * alpha - 2.0 * hm1 * h0 * hp2 * alpha - hm1 * hp1 * hp2 * alpha & ! End line 2
           + h0**2 * hp2 * alpha + h0 * hp1 * hp2 * alpha + hm1 * hp1 * hp2 * beta &
           + h0 * hp1 * hp2 * beta + hp1**2 * hp2 * beta ! End line 3
      coeff_b_deltas = coeff_b_deltas &
           / (hm1 * h0 * (h0 + hp1) * (h0 + hp1 + hp2))
    end associate
  end function coeff_b_deltas

  real function coeff_c_points(x)
    !! Compute the coefficient acting on f_{i+2} of the finite diference given a stencil of points.

    type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                !! difference.

    type(nucfd_stencil_deltas) :: h  ! Stencil of grid spacings for the finite difference.
    
    h = points_to_deltas(x)
    coeff_c_points = coeff_c_deltas(h)
    
  end function coeff_c_points

  pure real function coeff_c_deltas(h)
    !! Compute the coefficient acting on f_{i+21} of the finite diference given a stencil of grid
    !! spacings.

    type(nucfd_stencil_deltas), intent(in) :: h !! Stencil of grid spacings for the finite difference.

    real :: hm2, hm1, h0, hp1, hp2 ! Grid deltas at i -2, -1, 0, +1, +2

    hm2 = h%stencil(-2)
    hm1 = h%stencil(-1)
    h0 = h%stencil(0)
    hp1 = h%stencil(1)
    hp2 = h%stencil(2)

    associate(beta => alpha) ! To match Gamet et al. (1999)
      coeff_c_deltas = -hm1 * h0 * hp1 - h0**2 * hp1 + hm1 * h0**2 * alpha &
           + hm1 * h0 * hp1 * alpha  + hm1 * h0 * hp1 * beta + h0**2 * hp1 * beta & ! End line 1
           + hm1 * hp1**2 * beta + 2.0 * h0 * hp1**2 * beta + hp1**3 * beta ! End line 2
      coeff_c_deltas = coeff_c_deltas &
           / (hp2 * (hp1 + hp2) * (h0 + hp1 + hp2) * (hm1 + h0 + hp1 + hp2))
    end associate
  end function coeff_c_deltas
  
end module nucfd_coeffs
