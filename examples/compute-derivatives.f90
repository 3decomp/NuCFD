! examples/compute-derivatives.f90
!
!! Example program to compute derivatives.
!
! SPDX-License-Identifier: BSD-3-Clause

program compute_derivatives
  !! A program that uses non-uniform compact finite difference schemes to compute derivatives on
  !! (non-)uniform meshes.

  implicit none

  real, dimension(:), allocatable :: x ! 1D "grid"
  
  call build_grid(11, 1.0, (/ 1, 2 /), x)
  
contains

  subroutine build_grid(n, L, spacings, x)
    !! Crude subroutine to build a 1-D grid.
    !! The grid is described by the number of nodes, length and an array of relative grid spacings
    !! for each block (constant within a block), returning a 1-D array of node centres.
    !
    !! @note The number of cells (n - 1) must be divisible by the number of blocks specified by the
    !! spacings array otherwise an error will be raised.
    
    integer, intent(in) :: n                          !! The number of nodes in the mesh
    real, intent(in) :: L                             !! The length of the domain
    integer, dimension(:), intent(in) :: spacings     !! An array of the relative (integer) spacings
                                                      !! for each block
    real, dimension(:), allocatable, intent(out) :: x !! The grid node locations

    integer :: ncells
    integer :: nblocks
    integer :: block_size
    real :: block_length
    real :: h

    integer :: i
    integer :: b
    integer :: idx

    allocate(x(n))

    ncells = n - 1
    nblocks = size(spacings, 1)
    block_size = ncells / nblocks
    if (block_size * nblocks + 1 /= n) then
       print *, "Error: invalid grid specification!"
       print *, " - Grid size n should be specified so that n - 1 is divisible by number of blocks"
       print *, "   e.g. for 2 blocks a grid size of 11 is valid, 10 is not."
       error stop
    end if
    block_length = L / real(sum(spacings))
    h = block_length / real(block_size)

    idx = 1
    x(idx) = 0.0
    idx = 2
    do b = 1, nblocks
       do i = 1, block_size
          x(idx) = x(idx - 1) + real(spacings(b)) * h
          idx = idx + 1
       end do
    end do

    ! Self-check
    if (any(x > L) .or. any(x < 0.0)) then
       print *, "Error: grid exceeds specified length!"
       print *, x
       error stop
    end if
    if (abs(x(1)) > 0.0) then
       print *, "Error: start point is wrong!"
       print *, x(1)
       error stop
    end if
    if (abs(x(n) - L) > (2 * epsilon(L) * L)) then
       print *, "Error: end point is wrong!"
       print *, "- Computed: ", x(n)
       print *, "- Expected: ", L
       print *, "- Relative error: ", abs(x(n) - L) / L
       error stop
    end if
    
  end subroutine build_grid
  
end program compute_derivatives
