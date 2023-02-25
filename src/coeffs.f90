!!!! src/coeffs.f90
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

  implicit none

  private

  real, parameter, public :: alpha = 1.0 / 3.0 ! Off-diagonal coefficient for first-derivative
                                               ! system.
  
end module nucfd_coeffs
