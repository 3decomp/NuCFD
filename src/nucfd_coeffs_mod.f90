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
  
  implicit none

  private
  public :: coeff_a

  real, parameter, public :: alpha = 1.0 / 3.0 !! Off-diagonal coefficient for first derivative
                                               !! system.
  
contains

  pure real function coeff_a()
    !! Compute the coefficient acting on f_{i+1} of the finite diference.
    
    coeff_a = 0.0
    
  end function coeff_a
  
end module nucfd_coeffs
