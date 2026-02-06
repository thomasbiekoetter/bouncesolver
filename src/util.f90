module bouncesolver__util

  use bouncesolver__config, only : wp

  implicit none

  private

  real(wp), public, parameter :: pi = 4.0e0_wp * atan(1.0e0_wp)

  public :: riemann_integrate
  public :: linspace
  public :: is_equal

contains

  function riemann_integrate(x, f) result(I)

    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: f(:)
    real(wp) :: I

    integer :: n
    real(wp) :: dx
    real(wp) :: df

    I = 0.0e0_wp
    do n = 1, size(x) - 1
      dx = x(n + 1) - x(n)
      df = 0.5e0_wp * (f(n + 1) + f(n))
      I = I + dx * df
    end do

  end function

  pure function linspace(a, b, num) result(y)

    real(wp), intent(in) :: a
    real(wp), intent(in) :: b
    integer, intent(in) :: num
    real(wp), dimension(num) :: y

    integer :: i
    real(wp) :: step

    step = (b - a) / real(num - 1, wp)
    y(1) = a
    do i=2,num
      y(i) = a + (i - 1) * step
    end do

  end function linspace

  function is_equal(a, b, eps) result(y)

    real(wp), intent(in) :: a
    real(wp), intent(in) :: b
    real(wp), intent(in), optional :: eps
    logical :: y

    real(wp) :: tol

    tol = 1.0e-4_wp
    if (present(eps)) tol = eps

    y = (abs(a - b) <= tol)

  end function is_equal

end module bouncesolver__util
