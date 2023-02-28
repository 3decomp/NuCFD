submodule (nucfd_coeffs) nucfd_coeffs_c

  implicit none

contains

  module real function coeff_c_points(x)
    !! Compute the coefficient acting on f_{i+2} of the finite diference given a stencil of points.

    type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                !! difference.

    type(nucfd_stencil_deltas) :: h  ! Stencil of grid spacings for the finite difference.
    
    h = points_to_deltas(x)
    coeff_c_points = coeff_c_deltas(h)
    
  end function coeff_c_points

  module pure real function coeff_c_deltas(h)
    !! Compute the coefficient acting on f_{i+2} of the finite diference given a stencil of grid
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
      coeff_c_deltas = -hm1 * h0 * hp1 - h0**2 * hp1 + hm1 * h0**2 * alpha &
           + hm1 * h0 * hp1 * alpha  + hm1 * h0 * hp1 * beta + h0**2 * hp1 * beta & ! End line 1
           + hm1 * hp1**2 * beta + 2.0 * h0 * hp1**2 * beta + hp1**3 * beta ! End line 2
      coeff_c_deltas = coeff_c_deltas &
           / (hp2 * (hp1 + hp2) * (h0 + hp1 + hp2) * (hm1 + h0 + hp1 + hp2))
    end associate
  end function coeff_c_deltas
  
end submodule nucfd_coeffs_c
