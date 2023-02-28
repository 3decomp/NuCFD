submodule (nucfd_coeffs) nucfd_coeffs_e

  implicit none

contains

  module real function coeff_e(x)
    !! Compute the coefficient acting on f_{i} of the finite diference given a stencil of points.

    type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                !! difference.

    real :: a, b, c, d ! The neighbouring coefficients

    a = coeff_a(x)
    b = coeff_b(x)
    c = coeff_c(x)
    d = coeff_d(x)

    coeff_e = -(a + b + c + d)
    
  end function coeff_e
  
end submodule nucfd_coeffs_e
