submodule (nucfd_coeffs) nucfd_coeffs_b
  !! Submodule defining the coefficient acting on f_{i+1} for compact finite difference schemes on
  !! non-uniform grids.
  !!
  !! SPDX-License-Identifier: BSD-3-Clause

  implicit none
  
contains
  
  module real function coeff_b_points(x)
    !! Compute the coefficient acting on f_{i-1} of the finite diference given a stencil of points.

    type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                !! difference.

    type(nucfd_stencil_deltas) :: h  ! Stencil of grid spacings for the finite difference.
    
    h = points_to_deltas(x)
    coeff_b_points = coeff_b_deltas(h)
    
  end function coeff_b_points

  module real function coeff_b_deltas(h)
    !! Compute the coefficient acting on f_{i-1} of the finite diference given a stencil of grid
    !! spacings.

    type(nucfd_stencil_deltas), intent(in) :: h !! Stencil of grid spacings for the finite difference.

    real :: hm2, hm1, h0, hp1, hp2 ! Grid deltas at i -2, -1, 0, +1, +2

    select type(deltas => h%stencil)
    type is(real)
       hm2 = deltas(-2)
       hm1 = deltas(-1)
       h0 = deltas(0)
       hp1 = deltas(1)
       hp2 = deltas(2)
    class default
       error stop
    end select

    associate(beta => alpha) ! To match Gamet et al. (1999)
      coeff_b_deltas = coeff_numerator(hm1, h0, hp1, hp2, alpha)  ! => (14/3) h^3, h=const,
                                                                  ! alpha=1/3
      coeff_b_deltas = coeff_b_deltas &                           ! => Zero correction,
           + coeff_numerator_corr(hm1, h0, hp1, hp2, alpha, beta) ! h=const, alpha=beta
      coeff_b_deltas = coeff_b_deltas &
           /  coeff_denominator(hm1, h0, hp1)                     ! => -14/9, h=const
      coeff_b_deltas = coeff_b_deltas &
           / coeff_denominator_2h(h0, hp1, hp2)                   ! => -(14/9)/(2h), h = const
    end associate
  end function coeff_b_deltas

  pure real function coeff_numerator(hm1, h0, hp1, hp2, alpha)
    !! Computes the numerator of the coefficient acting on f_{i-1}.
    !!
    !! Reduces to -(14/3) h^3 when h=const and alpha=beta.
    
    real, intent(in) :: hm1 
    real, intent(in) :: h0 
    real, intent(in) :: hp1 
    real, intent(in) :: hp2
    real, intent(in) :: alpha
    
    coeff_numerator = -hp1 * ((hm1 + h0) * (hp1 + hp2) + 2.0 * hm1 * h0 * alpha) ! = -(14/3) h^3
  end function coeff_numerator

  pure real function coeff_numerator_corr(hm1, h0, hp1, hp2, alpha, beta)
    !! Computes the non-uniform correction to the numerator acting on f_{i-1}.
    !!
    !! Reduces to zero when h=const and alpha=beta.
    
    real, intent(in) :: hm1 
    real, intent(in) :: h0 
    real, intent(in) :: hp1 
    real, intent(in) :: hp2
    real, intent(in) :: alpha
    real, intent(in) :: beta

    coeff_numerator_corr = -hm1 * hp1 * hp2 * (alpha - beta) & ! Should cancel for case alpha = beta
         - (hp1**2) * (hm1 - h0) * alpha &    ! Should cancel for case h = const
         - h0 * hp2 * (2.0 * hm1 - h0 - hp1) * alpha &      ! Should cancel for case h = const
         - (h0**2) * (3.0 * hm1 - h0 - 2.0 * hp1) * alpha & ! Should cancel for case h = const
         - hp1 * (2.0 * hm1 * h0 * alpha - (h0 + hp1) * hp2 * beta) ! Should cancel for constant h &
                                                                    ! alpha = beta
  end function coeff_numerator_corr

  pure real function coeff_denominator(hm1, h0, hp1)
    !! Computes the denominator of the coefficient acting on f_{i-1}.
    !!
    !! Reduces to 3 h^3 when h=const. Dividing the numerator by this term yields the coefficient
    !! -14/9 when h=const, alpha=beta=1/3.

    real, intent(in) :: hm1
    real, intent(in) :: h0 
    real, intent(in) :: hp1
    
    coeff_denominator = 3.0 * hm1 * h0 * ((h0 + hp1) / 2.0)
  end function coeff_denominator

  pure real function coeff_denominator_2h(h0, hp1, hp2)
    !! Computes the non-uniform equivalent to 2h divisor of the coefficient acting on f_{i-1}.
    !!
    !! Dividing the coefficient by this term should reduce to (14/9)/(2h) when h=const,
    !! alpha=beta=1/3.

    real, intent(in) :: h0
    real, intent(in) :: hp1
    real, intent(in) :: hp2

    coeff_denominator_2h = 2.0 * (h0 + hp1 + hp2) / 3.0
  end function coeff_denominator_2h
  
end submodule nucfd_coeffs_b