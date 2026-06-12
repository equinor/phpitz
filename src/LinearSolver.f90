



   !==============================================================================
!module linear_solver
!!==============================================================================
!! Self-contained Gaussian elimination with partial pivoting.
!! Replaces LAPACK dgesv so no external library is needed.
!!==============================================================================
!  use gibbs_data, only: dp
!  implicit none
!
!contains

  subroutine solve_linear(n, A, b, x, info)
   use gibbs_data, only: dp
   implicit none
    !--------------------------------------------------------------------------
    ! Solve A*x = b using Gaussian elimination with partial pivoting.
    ! A is destroyed on output. b is not modified. Result in x.
    ! info = 0 on success, > 0 if matrix is singular.
    !--------------------------------------------------------------------------
    Integer, intent(in)     :: n
    Real(dp), intent(inout) :: A(n, n)
    Real(dp), intent(in)    :: b(n)
    Real(dp), intent(out)   :: x(n)
    Integer, intent(out)    :: info

    Real(dp) :: Ab(n, n+1)   ! augmented matrix
    Real(dp) :: max_val, temp_val, factor, sum_val
    Integer  :: i, j, k, max_row

    info = 0

    ! Build augmented matrix [A|b]
    Ab(1:n, 1:n) = A(1:n, 1:n)
    Ab(1:n, n+1) = b(1:n)

    ! Forward elimination with partial pivoting
    do k = 1, n
      ! Find pivot
      max_val = abs(Ab(k, k))
      max_row = k
      do i = k+1, n
        if (abs(Ab(i, k)) > max_val) then
          max_val = abs(Ab(i, k))
          max_row = i
        end if
      end do

      ! Check for singular matrix
      if (max_val < 1.0e-300_dp) then
        info = k
        x = 0.0_dp
        return
      end if

      ! Swap rows if needed
      if (max_row /= k) then
        do j = k, n+1
          temp_val = Ab(k, j)
          Ab(k, j) = Ab(max_row, j)
          Ab(max_row, j) = temp_val
        end do
      end if

      ! Eliminate below
      do i = k+1, n
        factor = Ab(i, k) / Ab(k, k)
        do j = k+1, n+1
          Ab(i, j) = Ab(i, j) - factor * Ab(k, j)
        end do
        Ab(i, k) = 0.0_dp
      end do
    end do

    ! Back substitution
    do i = n, 1, -1
      sum_val = Ab(i, n+1)
      do j = i+1, n
        sum_val = sum_val - Ab(i, j) * x(j)
      end do
      x(i) = sum_val / Ab(i, i)
    end do

  end subroutine solve_linear

