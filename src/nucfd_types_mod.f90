!!!! src/nucfd_types_mod.f90
!!!
!!!! Description
!!!
!!! Defines the custom types used by NuCFD.
!!!
!!! Provides the nucfd_types module.
!!!
!!!! LICENSE
!!!
!!! SPDX-License-Identifier: BSD-3-Clause
!!!

module nucfd_types

  implicit none

  type nucfd_stencil(width, centre)
     !! Type representing a stencil.
     integer, kind :: width = 1  !! Stencil width.
     integer, kind :: centre = 1 !! Stencil centre.
     real, dimension(1-centre:(1-centre)+(width-1)) :: stencil !! Stencil data. The lbound, ubound
                                                               !! are set to allow stencil access by
                                                               !! offsets.
  end type nucfd_stencil
  
  type, extends(nucfd_stencil) :: nucfd_stencil_points
     !! Type representing the point coordinates of a stencil.
  end type nucfd_stencil_points

  type, extends(nucfd_stencil) :: nucfd_stencil_deltas
     !! Type representing the grid spacings of a stencil. Following Gamet et al. (1999) these are
     !! defined as h_i = x_i - x_{i-1}.
  end type nucfd_stencil_deltas

contains

  pure function points_to_deltas(x) result(h)
    !! Helper function converting a stencil's point representation to a grid spacing representation.

    type(nucfd_stencil_points(6, 4)), intent(in) :: x !! The stencil's point representation.
    type(nucfd_stencil_deltas(5, 3)) :: h             !! The stencil's grid spacing representation.

   end function points_to_deltas
  
end module nucfd_types
