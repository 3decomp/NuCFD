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

  real function coeff_d_points(x)
    !! Compute the coefficient acting on f_{i-2} of the finite diference given a stencil of points.

    type(nucfd_stencil_points), intent(in) :: x !! Stencil of points for the finite
                                                !! difference.

    type(nucfd_stencil_deltas) :: h  ! Stencil of grid spacings for the finite difference.
    
    h = points_to_deltas(x)
    coeff_d_points = coeff_d_deltas(h)
    
  end function coeff_d_points

  pure real function coeff_d_deltas(h)
    !! Compute the coefficient acting on f_{i-2} of the finite diference given a stencil of grid
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
      coeff_d_deltas = h0 * hp1**2 + h0 * hp1 * hp2 - h0**3 * alpha - 2.0 * h0**2 * hp1 * alpha &
           - h0 * hp1**2 * alpha - h0**2 * hp2 * alpha - h0 * hp1 * hp2 * alpha & ! End line 1
           - h0 * hp1 * hp2 * beta - hp1**2 * hp2 * beta ! End line 2
      coeff_d_deltas = coeff_d_deltas &
           / (hm1 * (hm1 + h0) * (hm1 + h0 + hp1) * (hm1 + h0 + hp1 + hp2))
    end associate
  end function coeff_d_deltas

  real function coeff_e(x)
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
end module nucfd_coeffs
