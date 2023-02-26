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

  type nucfd_stencil
     !! Type representing a stencil.
     real, dimension(:), allocatable :: stencil !! Stencil data
   contains
     final :: free_stencil
  end type nucfd_stencil
  
  type, extends(nucfd_stencil) :: nucfd_stencil_points
     !! Type representing the point coordinates of a stencil.
  end type nucfd_stencil_points

  type, extends(nucfd_stencil) :: nucfd_stencil_deltas
     !! Type representing the grid spacings of a stencil. Following Gamet et al. (1999) these are
     !! defined as h_i = x_i - x_{i-1}.
  end type nucfd_stencil_deltas

contains

  subroutine create_stencil(width, centre, stencil)
    !! Subroutine to create a stencil, specifying the width and centre. The underlying array is
    !! allocated with bounds (1-centre):(1-centre)+(width-1) so that accesses are by offset from the
    !! central position.

    integer, intent(in) :: width                 !! Specifies the width of the stencil
    integer, intent(in) :: centre                !! Specifies the index the stencil should be centred about
    class(nucfd_stencil), intent(out) :: stencil !! The stencil object.

    integer :: lb, ub

    lb = 1 - centre
    ub = (1 - centre) + (width - 1)
    allocate(stencil%stencil(lb:ub))
    
  end subroutine create_stencil
  
  function points_to_deltas(x) result(h)
    !! Helper function converting a stencil's point representation to a grid spacing representation.

    type(nucfd_stencil_points), intent(in) :: x !! The stencil's point representation.
    type(nucfd_stencil_deltas) :: h             !! The stencil's grid spacing representation.

    integer :: i
    real :: xm1, x0 ! Start and end points of grid spacing at i.

    integer :: ub, lb, width, centre ! Upper bound, lower bound, width and centre of grid stencil.

    lb = lbound(x%stencil, 1) + 1
    ub = ubound(x%stencil, 1)
    width = (ub - lb) + 1
    centre = 1 - lb
    call create_stencil(width, centre, h)
    
    do i = lbound(h%stencil, 1) + 1, ubound(h%stencil, 1)
       xm1 = x%stencil(i - 1)
       x0 = x%stencil(i)

       h%stencil(i) = x0 - xm1
    end do
    
   end function points_to_deltas

   subroutine free_stencil(stencil)
     !! Frees a stencil object
     
     type(nucfd_stencil) :: stencil !! The stencil object to be freed.

     if (allocated(stencil%stencil)) then
        deallocate(stencil%stencil)
     end if
     
   end subroutine free_stencil
   
end module nucfd_types
