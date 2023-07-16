! src/coeffs/nucfd_coeff_c.f90
!
!! Implements coefficient C of the nucfd_coeffs module.
!
! SPDX-License-Identifier: BSD-3-Clause

submodule (nucfd_coeffs) nucfd_coeffs_c
  !! Submodule defining the coefficient acting on f_{i+2} for compact finite difference schemes on
  !! non-uniform grids.

  implicit none

contains

  module real function coeff_c_points(x)
    !! Compute the coefficient acting on f_{i+2} of the finite diference given a stencil of points.

    type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                !! difference.

    type(nucfd_stencil_deltas) :: h  ! Stencil of grid spacings for the finite difference.

#ifndef NDEBUG
    if (size(x%stencil) /= 5) then
       print *, "Error@coeff_c_points: expecting width 5 stencil, received width = ", &
            size(x%stencil)
       error stop
    end if
    if (1 - lbound(x%stencil, 1) /= 3) then
       print *, "Error@coeff_c_points: expecting centre 3 stencil, received centre = ", &
            1 - lbound(x%stencil, 1)
       print *, size(x%stencil), lbound(x%stencil), ubound(x%stencil)
       error stop
    end if
#endif
    
    h = points_to_deltas(x)
    coeff_c_points = coeff_c_deltas(h)
    
  end function coeff_c_points

  module real function coeff_c_deltas(h)
    !! Compute the coefficient acting on f_{i+2} of the finite diference given a stencil of grid
    !! spacings.

    type(nucfd_stencil_deltas), intent(in) :: h !! Stencil of grid spacings for the finite difference.

    real :: numerator, denominator, divisor

#ifndef NDEBUG
    if (size(h%stencil) /= 4) then
       print *, "Error@coeff_c_deltas: expecting width 4 stencil, received width = ", &
            size(h%stencil)
       error stop
    end if
    if (1 - lbound(h%stencil, 1) /= 2) then
       print *, "Error@coeff_c_deltas: expecting centre 2 stencil, received centre = ", &
            1 - lbound(h%stencil, 1)
       print *, size(h%stencil), lbound(h%stencil), ubound(h%stencil)
       error stop
    end if
#endif
    
    call coeff_c_components(h, numerator, denominator, divisor)
    coeff_c_deltas = (numerator / denominator) / divisor
    
  end function coeff_c_deltas

  module subroutine coeff_c_components(h, numerator, denominator, divisor)
    !! Compute the components of the coefficient acting on f_{i+2} of the finite diference given a
    !! stencil of grid spacings.
    !!
    !! For uniform grids the numerator should reduce to (2/3) h^3, the denominator to 6h^3 (for a
    !! coefficient value of 1/9) and the finite difference divisor to 4h.

    type(nucfd_stencil_deltas), intent(in) :: h !! Stencil of grid spacings for the finite difference.
    real, intent(out) :: numerator   !! The numerator of the coefficient
    real, intent(out) :: denominator !! The denominator of the coefficient
    real, intent(out) :: divisor     !! The finite-difference divisor of the coefficient
    
    real :: hm1, h0, hp1, hp2 ! Grid deltas at i -1, 0, +1, +2

#ifndef NDEBUG
    if (size(h%stencil) /= 4) then
       print *, "Error@coeff_c_components: expecting width 4 stencil, received width = ", &
            size(h%stencil)
       error stop
    end if
    if (1 - lbound(h%stencil, 1) /= 2) then
       print *, "Error@coeff_c_components: expecting centre 2 stencil, received centre = ", &
            1 - lbound(h%stencil, 1)
       print *, size(h%stencil), lbound(h%stencil), ubound(h%stencil)
       error stop
    end if
#endif

    select type(deltas => h%stencil)
    type is(real)
       hm1 = deltas(-1)
       h0 = deltas(0)
       hp1 = deltas(1)
       hp2 = deltas(2)
    class default
       error stop
    end select

    associate(beta => alpha) ! To match Gamet et al. (1999)
      numerator = coeff_numerator(hm1, h0, hp1, alpha, beta)
      denominator = coeff_denominator(h0, hp1, hp2)
      divisor = coeff_divisor(hm1, h0, hp1, hp2)
    end associate
    
  end subroutine coeff_c_components
    
  pure real function coeff_numerator(hm1, h0, hp1, alpha, beta)
    !! Computes the numerator of the coefficient acting on f_{i+2}.
    !!
    !! Reduces to (2/3) h^3 when h=const and alpha=beta=1/3.

    real, intent(in) :: hm1
    real, intent(in) :: h0
    real, intent(in) :: hp1
    real, intent(in) :: alpha
    real, intent(in) :: beta
    
    coeff_numerator = -h0 * hp1 * (hm1 + h0) &
         + hm1 * h0 * (h0 + hp1) * alpha &
         + hp1 * (h0 * (hm1 + h0) + hp1 * (hm1 + 2.0 * h0 + hp1)) * beta
    
  end function coeff_numerator
  
  pure real function coeff_denominator(h0, hp1, hp2)
    !! Computes the denominator of the coefficient acting on f_{i+2}.
    !!
    !! Reduces to 6 h^3 when h=const. Dividing the numerator by this term yields the coefficient
    !! 1/9 when h=const, alpha=beta=1/3.

    real, intent(in) :: h0
    real, intent(in) :: hp1
    real, intent(in) :: hp2

    coeff_denominator = hp2 * (hp1 + hp2) &
         * (h0 + hp1 + hp2)
    
  end function coeff_denominator

  pure real function coeff_divisor(hm1, h0, hp1, hp2)
    !! Computes the non-uniform equivalent to 4h divisor of the coefficient acting on f_{i+2}.
    !!
    !! Dividing the coefficient by this term should reduce to (1/9)/(4h) when h=const,
    !! alpha=beta=1/3.

    real, intent(in) :: hm1
    real, intent(in) :: h0
    real, intent(in) :: hp1
    real, intent(in) :: hp2

    coeff_divisor = (hm1 + h0 + hp1 + hp2) 
  end function coeff_divisor
  
end submodule nucfd_coeffs_c