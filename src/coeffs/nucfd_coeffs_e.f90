! src/coeffs/nucfd_coeff_e.f90
!
!! Implements coefficient E of the nucfd_coeffs module.
!
! SPDX-License-Identifier: BSD-3-Clause

submodule (nucfd_coeffs) nucfd_coeffs_e
  !! Submodule defining the coefficient acting on f_i for compact finite difference schemes on
  !! non-uniform grids.

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
