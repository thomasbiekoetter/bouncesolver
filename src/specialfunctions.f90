module bouncesolver__specialfunction

  use bouncesolver__config, only : wp
  use bouncesolver__util, only : pi

  implicit none

  private

  public :: bessel_mod_1
  public :: bessel_mod_onehalf

  real(wp), parameter :: p1 = 0.5e0_wp
  real(wp), parameter :: p2 = 0.87890594e0_wp
  real(wp), parameter :: p3 = 0.51498869e0_wp
  real(wp), parameter :: p4 = 0.15084934e0_wp
  real(wp), parameter :: p5 = 0.2658733e-1_wp
  real(wp), parameter :: p6 = 0.301532e-2_wp
  real(wp), parameter :: p7 = 0.32411e-3_wp

  real(wp), parameter :: q1 = 0.39894228e0_wp
  real(wp), parameter :: q2 = -0.3988024e-1_wp
  real(wp), parameter :: q3 = -0.362018e-2_wp
  real(wp), parameter :: q4 = 0.163801e-2_wp
  real(wp), parameter :: q5 = -0.1031555e-1_wp
  real(wp), parameter :: q6 = 0.2282967e-1_wp
  real(wp), parameter :: q7 = -0.2895312e-1_wp
  real(wp), parameter :: q8 = 0.1787654e-1_wp
  real(wp), parameter :: q9 = -0.420059e-2_wp

contains

  ! Source: NUMERICAL RECIPES IN FORTRAN 77
  pure function bessel_mod_1(x) result(res)

    real(wp), intent(in) :: x
    real(wp) :: res

    real(wp) :: ax
    real(wp) :: y

    if (abs(x) < 3.75e0_wp) then
      y = (x / 3.75e0_wp) ** 2
      res = x * (p1 + y * (p2 + y * (p3 +  &
        y * (p4 + y * (p5 + y * (p6 + y * p7))))))
    else
      ax = abs(x)
      y = 3.75e0_wp / ax
      res = (exp(ax) / sqrt(ax)) * (q1 +  &
        y * (q2 + y * (q3 + y * (q4 +  &
        y * (q5 + y * (q6 + y * (q7 + &
        y * (q8 + y * q9))))))))
      if (x < 0.0e0_wp) res = -res
    endif

  end function bessel_mod_1

  pure function bessel_mod_onehalf(x) result(res)

    real(wp), intent(in) :: x
    real(wp) :: res

    if (abs(x) > 1.0e-20_wp) then
      res = sqrt(2.0e0_wp / (pi * x)) * sinh(x)
    else
      res = sqrt(2.0e0_wp * x / pi)
    end if

  end function bessel_mod_onehalf

end module bouncesolver__specialfunction
