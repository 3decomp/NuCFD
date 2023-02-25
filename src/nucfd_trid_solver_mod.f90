!!!! src/nucfd_trid_solver_mod.f90
!!!
!!!! Description
!!!
!!! Implements a simple tridiagonal solver based on the algorithm given at
!!! https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
!!!
!!! Provides the nucfd_trid_solver module.
!!!
!!!! LICENSE
!!!
!!! SPDX-License-Identifier: BSD-3-Clause
!!!

module nucfd_trid_solver

  implicit none

  private

  public :: solve
  public :: solve_cyclic
  
contains
  
  subroutine solve(a, b, c, rhs, x)
    ! Solves a tridiagonal system using the Thomas algorithm.
    
    real, dimension(:), intent(in) :: a   ! The sub-diagonal coefficients vector.
    real, dimension(:), intent(in) :: b   ! The diagonal coefficients vector.
    real, dimension(:), intent(in) :: c   ! The super-diagonal coefficients vector.
    real, dimension(:), intent(in) :: rhs ! The right hand side vector.
    real, dimension(:), intent(out) :: x  ! The solution vector.

    real, dimension(:), allocatable :: bp ! The modified diagonal.

    integer :: n

    n = size(x)
    
    allocate(bp(n))
    call forward_sweep(a, b, c, rhs, bp, x)
    call backward_sweep(bp, c, x)
    deallocate(bp)

  end subroutine solve

  subroutine solve_cyclic(a, b, c, rhs, x)
    ! Solves a cyclic tridiagonal system using the Thomas algorithm.
    
    real, dimension(:), intent(in) :: a   ! The sub-diagonal coefficients vector.
    real, dimension(:), intent(in) :: b   ! The diagonal coefficients vector.
    real, dimension(:), intent(in) :: c   ! The super-diagonal coefficients vector.
    real, dimension(:), intent(in) :: rhs ! The right hand side vector.
    real, dimension(:), intent(out) :: x  ! The solution vector.

    real, dimension(:), allocatable :: bp   ! The modified diagonal
    real, dimension(:), allocatable :: q, u ! Vectors of the augmented system
    real :: v1, vn
    real :: vx, vq
    real :: gamma

    integer :: n
    
    n = size(x)

    allocate(bp(n), q(n), u(n))
    
    ! Create perturbed system
    bp(:) = b(:)
    q(:) = 0.0
    u(:) = 0.0

    gamma = -b(1)

    u(1) = gamma
    u(n) = c(n)

    v1 = 1.0
    vn = a(1) / gamma

    bp(1)= bp(1) - gamma
    bp(n) = bp(n) - a(1) * c(n) / gamma
    
    ! Solve perturbed systems
    call solve(a, bp, c, rhs, x)
    call solve(a, bp, c, u, q)
    
    ! Recontruct the solution
    vx = (v1 * x(1) + vn * x(n))
    vq = (v1 * q(1) + vn * q(n))

    x(:) = x(:) - q(:) * (vx / (1.0 + vq))
    
    deallocate(bp, q, u)
    
  end subroutine solve_cyclic
  
  pure subroutine forward_sweep(a, b, c, rhs, bp, x)
    ! The forward sweep of the Thomas algorithm.
    
    real, dimension(:), intent(in) :: a   ! The sub-diagonal coefficients vector.
    real, dimension(:), intent(in) :: b   ! The diagonal coefficients vector.
    real, dimension(:), intent(in) :: c   ! The super-diagonal coefficients vector.
    real, dimension(:), intent(in) :: rhs ! The right hand side vector.
    real, dimension(:), intent(out) :: bp ! The modified diagonal
    real, dimension(:), intent(out) :: x  ! The solution vector.

    integer :: n
    integer :: i
    real :: w
    
    n = size(x)

    bp(1) = b(1)
    x(1) = rhs(1)
    do i = 2, n
       w = a(i) / bp(i - 1)
       bp(i) = b(i) - w * c(i - 1)
       x(i) = rhs(i) - w * x(i - 1)
    end do
    
  end subroutine forward_sweep

  pure subroutine backward_sweep(bp, c, x)
    ! The backward sweep of the Thomas algorithm.

    real, dimension(:), intent(in) :: bp   ! The modified diagonal
    real, dimension(:), intent(in) :: c    ! The super-diagonal coefficients vector.
    real, dimension(:), intent(inout) :: x ! The solution vector.

    integer :: n
    integer :: i

    n = size(x)

    x(n) = x(n) / bp(n)
    do i = n - 1, 1, -1
       x(i) = (x(i) - c(i) * x(i + 1)) / bp(i)
    end do

  end subroutine backward_sweep
  
end module nucfd_trid_solver
