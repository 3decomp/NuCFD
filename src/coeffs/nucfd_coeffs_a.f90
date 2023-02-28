submodule (nucfd_coeffs) nucfd_coeffs_a

  implicit none

contains

  module real function coeff_a_points(x)
    !! Compute the coefficient acting on f_{i+1} of the finite diference given a stencil of points.

    type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                !! difference.

    type(nucfd_stencil_deltas) :: h  ! Stencil of grid spacings for the finite difference.
    
    h = points_to_deltas(x)
    coeff_a_points = coeff_a(h)
    
  end function coeff_a_points

  module pure real function coeff_a_deltas(h)
    !! Compute the coefficient acting on f_{i+1} of the finite diference given a stencil of grid
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
      coeff_a_deltas = h0 * ((hm1 + h0) * (hp1 + hp2) + 2.0 * hp1 * hp2 * beta) & ! = (14/3) h^3
           + hm1 * h0 * hp2 * (beta - alpha) & ! Should cancel for case alpha = beta
           + (h0**2) * (hp2 - hp1) * beta    & ! Should cancel for constant h
           + hm1 * hp1 * (2.0 * hp2 - h0 - hp1) * beta &      ! Should cancel for constant h
           + (hp1**2) * (3.0 * hp2 - 2.0 * h0 - hp1) * beta & ! Should cancel for constant h
           + h0 * (2.0 * hp1 * hp2 * beta - hm1 * (h0 + hp1) * alpha) ! Should cancel for constant h
           
      coeff_a_deltas = coeff_a_deltas &
           / (3.0 * hp1 * ((h0 + hp1) / 2.0) * hp2) ! => 14.0 / 9.0 when h = const

      coeff_a_deltas = coeff_a_deltas &
              / (2.0 * ((hm1 + h0 + hp1) / 3.0)) ! => (14.0 / 9.0) / (2h) when h = const
    end associate
  end function coeff_a_deltas
  
end submodule nucfd_coeffs_a
