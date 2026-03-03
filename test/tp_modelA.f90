program bouncesolver__test_tp_modelA

  use bouncesolver__config, only : wp
  use bouncesolver__config, only : fmt2
  use bouncesolver__tp, only : solver
  use gradmin__descent, only : minimize

  implicit none

  real(wp), parameter :: vev = 246.0e0_wp
  real(wp), parameter :: musq = 0.05e0_wp * vev ** 2
  real(wp), parameter :: lam = 0.1e0_wp
  real(wp), parameter :: lamb = 0.1e0_wp

  type(solver) :: tp
  integer, parameter :: d = 2
  real(wp) :: phi_F(d) = [0.0e0_wp, vev + 1.0e0_wp]
  real(wp) :: phi_T(d) = [vev + 1.0e0_wp, 0.0e0_wp]

  call get_minima(phi_F, phi_T)

  tp = solver(  &
    d, V, phi_F, phi_T,  &
    alpha=2,  &
    nx=400)

contains

  function V(phi) result(y)

    real(wp), intent(in) :: phi(:)
    real(wp) :: y

    real(wp) :: p
    real(wp) :: s

    p = phi(1)
    s = phi(2)

    y = lam * (p ** 2 + s ** 2 - vev ** 2) ** 2 +  &
      lamb * p ** 2 * s ** 2 - musq * p ** 2

  end function V

  subroutine get_minima(phi_F, phi_T)

    real(wp), intent(inout) :: phi_F(d)
    real(wp), intent(inout) :: phi_T(d)

    real(wp) :: phi_min(d)
    real(wp) :: V_min

    call minimize(  &
      V, phi_T, phi_min, V_min,  &
      maxiter=10000, mode=1)

    if (phi_min(1) > 0.0e0_wp) then
      phi_T = phi_min
    else
      phi_T = -phi_min
    end if

    write(*,"(a)", advance="no")  &
      " True minimum at phi / v ="
    write(*,fmt2) phi_T / vev
    write(*,*) "                V / v^4 =",  &
      V_min / vev ** 4

    call minimize(  &
      V, phi_F, phi_min, V_min,  &
      maxiter=10000, mode=1)

    phi_F = phi_min

    write(*,"(a)", advance="no")  &
      "False minimum at phi / v ="
    write(*,fmt2) phi_F
    write(*,*) "                V / v^4 =",  &
      V_min / vev ** 4

  end subroutine get_minima

end program bouncesolver__test_tp_modelA
