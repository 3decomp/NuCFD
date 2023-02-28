module nucfd_coeffs
  !! Module defining the coefficients for compact finite difference schemes on non-uniform grids.
  !!
  !! SPDX-License-Identifier: BSD-3-Clause

  use nucfd_types
  
  implicit none

  private
  public :: coeff_a
  public :: coeff_b
  public :: coeff_c
  public :: coeff_d
  public :: coeff_e

  interface coeff_a
     !! Compute the coefficient acting on f_{i+1} of the finite diference stencil.
     module procedure coeff_a_points
     module procedure coeff_a_deltas
  end interface coeff_a
  
  interface coeff_b
     !! Compute the coefficient acting on f_{i-1} of the finite diference stencil.
     module procedure coeff_b_points
     module procedure coeff_b_deltas
  end interface coeff_b

  interface coeff_c
     !! Compute the coefficient acting on f_{i+2} of the finite diference stencil.
     module procedure coeff_c_points
     module procedure coeff_c_deltas
  end interface coeff_c

  interface coeff_d
     !! Compute the coefficient acting on f_{i-2} of the finite diference stencil.
     module procedure coeff_d_points
     module procedure coeff_d_deltas
  end interface coeff_d

  interface
     module real function coeff_a_points(x)
       type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                   !! difference.
     end function
     module real function coeff_a_deltas(h)
       type(nucfd_stencil_deltas), intent(in) :: h !! Stencil of grid spacings for the finite
                                                   !! difference.
     end function
     module real function coeff_b_points(x)
       type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                   !! difference.
     end function
     module real function coeff_b_deltas(h)
       type(nucfd_stencil_deltas), intent(in) :: h !! Stencil of grid spacings for the finite
                                                   !! difference.
     end function
     module real function coeff_c_points(x)
       type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                   !! difference.
     end function
     module pure real function coeff_c_deltas(h)
       type(nucfd_stencil_deltas), intent(in) :: h !! Stencil of grid spacings for the finite
                                                   !! difference.
     end function
     module real function coeff_d_points(x)
       type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                   !! difference.
     end function
     module pure real function coeff_d_deltas(h)
       type(nucfd_stencil_deltas), intent(in) :: h !! Stencil of grid spacings for the finite
                                                   !! difference.
     end function
     module real function coeff_e(x)
       type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                   !! difference.
     end function
  end interface
  
  real, parameter, public :: alpha = 1.0 / 3.0 !! Off-diagonal coefficient for first derivative
                                               !! system.

end module nucfd_coeffs
