! examples/compute-derivatives.f90
!
!! Example program to compute derivatives.
!
! SPDX-License-Identifier: BSD-3-Clause

program compute_derivatives
  !! A program that uses non-uniform compact finite difference schemes to compute derivatives on
  !! (non-)uniform meshes.

  use nucfd_types
  use nucfd_deriv
  use nucfd_trid_solver
  
  implicit none

  real, parameter :: pi = 4.0 * atan(1.0)
  real, parameter :: alpha = 1.0 / 3.0
  
  integer :: n
  real :: L
  integer, dimension(:), allocatable :: block_spacings
  logical :: periodic
  real, dimension(:), allocatable :: x ! 1D "grid"

  type(nucfd_index_stencil) :: stencil
  integer :: width, centre

  integer :: idx, i

  real, dimension(:), allocatable :: f
  real, dimension(:), allocatable :: dfdx

  real, dimension(:), allocatable :: a, b, c, r
  
  n = 66
  L = 1.0
  block_spacings = (/ 1, 2, 1 /)
  periodic = .true.
  call build_grid(n, L, block_spacings, periodic, x)

  print *, "+++ Built grid +++"
  print *, "- lbounds: ", lbound(x)
  print *, "- ubounds: ", ubound(x)
  print *, "- 1st/last node: ", x(1), x(n)
  if (periodic) then
     print *, "- periodic node: ", x(n + 1)
  end if

  allocate(f, mold=x)
  f(:) = sin((x(:) / L) * (2.0 * pi))
  
  allocate(a(n))
  allocate(b(n))
  allocate(c(n))
  allocate(r(n))
  allocate(dfdx(n))
  
  width = 5
  centre = 3
  call create_stencil(width, centre, stencil)

  print *, "+++ Computing RHS +++"
  select type (indices => stencil%stencil)
  type is(integer)
     do idx = 1, n
        ! Move stencil
        do i = lbound(indices, 1), ubound(indices, 1)
           indices(i) = idx + i 
        end do
        
        call deriv_rhs(f, stencil, x, r(idx))
     end do
  end select

  print *, "+++ Solving system +++"
  a(:) = alpha
  b(:) = 1.0
  c(:) = alpha
  call solve_cyclic(a, b, c, r, dfdx)

  open (unit=10, file="dfdx.dat")
  do i = 1, n
     write(10, *) x(i), (2.0 * pi / L) * cos((x(i) / L) * (2.0 * pi)), dfdx(i)
  end do
  close(10)

  deallocate(f, dfdx)
  deallocate(x)

  deallocate(a, b, c, r)
  
contains

  subroutine build_grid(n, L, spacings, periodic, x)
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
    logical, intent(in) :: periodic                   !! Flag indicating periodicity of mesh
    real, dimension(:), allocatable, intent(out) :: x !! The grid node locations

    integer :: ncells
    integer :: nblocks
    integer :: block_size
    real :: block_length
    real :: h

    integer :: i
    integer :: b
    integer :: idx

    integer, parameter :: lghost = 2 ! Ghost node count at left boundary
    integer, parameter :: rghost = 2 ! Ghost node count at right boundary

    allocate(x(1-lghost:n+rghost))

    if (.not. periodic) then
       ncells = n - 1
    else
       ncells = n
    end if
    nblocks = size(spacings, 1)
    if (.not. periodic) then
       block_size = ncells / nblocks
    else
       block_size = ncells / nblocks
    end if
    if (.not. periodic) then
       if (block_size * nblocks + 1 /= n) then
          print *, "Error: invalid grid specification!"
          print *, " - Grid size n should be specified so that n - 1 is divisible by number of blocks"
          print *, "   e.g. for 2 blocks a grid size of 11 is valid, 10 is not."
          error stop
       end if
    else
       if (block_size * nblocks /= n) then
          print *, "Error: invalid grid specification!"
          error stop
       end if
    end if
    block_length = L / real(sum(spacings))
    h = block_length / real(block_size)
    
    idx = 1
    x(idx) = 0.0
    do b = 1, nblocks
       if (periodic .and. (b == nblocks)) then
          block_size = block_size - 1
       end if
       do i = 1, block_size
          idx = idx + 1
          x(idx) = x(idx - 1) + real(spacings(b)) * h
       end do
    end do
    if (idx /= n) then
       print *, "Error: didn't fill all interior nodes!"
       error stop
    end if

    ! Add ghost points
    if (.not. periodic) then
       do i = 1, lghost
          idx = 1 - i
          x(idx) = x(idx + 1) - real(spacings(1)) * h
       end do
       do i = 1, rghost
          idx = n + i
          x(idx) = x(idx - 1) + real(spacings(nblocks)) * h
       end do
    else
       idx = n + 1
       x(idx) = L
       do i = 2, rghost
          idx = n + i
          x(idx) = x(idx - 1) + (x(i) - x(i - 1))
       end do

       do i = 1, lghost
          idx = 1 - i
          x(idx) = x(idx + 1) - (x((n + 1) - (i - 1)) - (x((n + 1) - (i - 1) - 1)))
       end do
    end if

    ! Self-check
    if (any(x(1:n) > L) .or. any(x(1:n) < 0.0)) then
       print *, "Error: grid exceeds specified length!"
       print *, x
       error stop
    end if
    if (abs(x(1)) > 0.0) then
       print *, "Error: start point is wrong!"
       print *, x(1)
       error stop
    end if

    if (.not. periodic) then
       idx = n
    else
       idx = n + 1
    end if
    if (abs(x(idx) - L) > (2 * epsilon(L) * L)) then
       print *, "Error: end point is wrong!"
       print *, "- Computed: ", x(idx)
       print *, "- Expected: ", L
       print *, "- Relative error: ", abs(x(idx) - L) / L
       error stop
    end if
    
  end subroutine build_grid
  
end program compute_derivatives
