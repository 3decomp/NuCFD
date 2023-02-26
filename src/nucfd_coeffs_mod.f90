module nucfd_coeffs
  !! Module defining the coefficients for compact finite difference schemes on non-uniform grids.
  !!
  !! SPDX-License-Identifier: BSD-3-Clause
  
  implicit none

  private

  real, parameter, public :: alpha = 1.0 / 3.0 !! Off-diagonal coefficient for first derivative
                                               !! system.
  
end module nucfd_coeffs
