module nucfd_coeffs
  !! Module defining the coefficients for compact finite difference schemes on non-uniform grids.
  !!
  !! SPDX-License-Identifier: BSD-3-Clause

  use nucfd_types
  
  implicit none

  private
  public :: coeff_a
  public :: coeff_b
  public :: coeff_c, coeff_c_components
  public :: coeff_d, coeff_d_components
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
     module real function coeff_c_deltas(h)
       type(nucfd_stencil_deltas), intent(in) :: h !! Stencil of grid spacings for the finite
                                                   !! difference.
     end function
     module subroutine coeff_c_components(h, numerator, denominator, divisor)
       !! Compute the components of the coefficient acting on f_{i+2} of the finite diference given a
       !! stencil of grid spacings.
       !!
       !! For uniform grids the numerator should reduce to (2/3) h^3, the denominator to 6h^3 (for a
       !! coefficient value of (1/9) h^3) and the finite difference divisor to 4h.
       type(nucfd_stencil_deltas), intent(in) :: h !! Stencil of grid spacings for the finite difference.
       real, intent(out) :: numerator   !! The numerator of the coefficient
       real, intent(out) :: denominator !! The denominator of the coefficient
       real, intent(out) :: divisor     !! The finite-difference divisor of the coefficient
     end subroutine coeff_c_components
     
     module real function coeff_d_points(x)
       type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                   !! difference.
     end function
     module real function coeff_d_deltas(h)
       type(nucfd_stencil_deltas), intent(in) :: h !! Stencil of grid spacings for the finite
                                                   !! difference.
     end function
     module subroutine coeff_d_components(h, numerator, denominator, divisor)
       !! Compute the components of the coefficient acting on f_{i-2} of the finite diference given a
       !! stencil of grid spacings.
       !!
       !! For uniform grids the numerator should reduce to -(2/3) h^3, the denominator to 6h^3 (for a
       !! coefficient value of -1/9) and the finite difference divisor to 4h.
       type(nucfd_stencil_deltas), intent(in) :: h !! Stencil of grid spacings for the finite difference.
       real, intent(out) :: numerator   !! The numerator of the coefficient
       real, intent(out) :: denominator !! The denominator of the coefficient
       real, intent(out) :: divisor     !! The finite-difference divisor of the coefficient
     end subroutine coeff_d_components
  
     module real function coeff_e(x)
       type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                   !! difference.
     end function
  end interface
  
  real, parameter, public :: alpha = 1.0 / 3.0 !! Off-diagonal coefficient for first derivative
                                               !! system.

end module nucfd_coeffs
