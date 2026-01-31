program bouncesolver__app_pdm_thickwalled

  use gradmin__config, only : wp
  use gradmin__descent, only : minimize
  use bouncesolver__pdm, only : solver
  use bouncesolver__util, only : linspace
  use csv_module

  implicit none

  real(wp), parameter :: a = 1.8e0_wp
  real(wp), parameter :: b = 0.2e0_wp
  real(wp), parameter :: cmin = 0.1e0_wp
  real(wp), parameter :: cmax = 0.6e0_wp
  integer, parameter :: num = 101
  real(wp) :: clist(num)
  real(wp) :: c

  type(solver) :: pd1
  type(solver) :: pd2
  integer, parameter :: d = 2
  real(wp) :: phi_false(d) = [0.0e0_wp, 0.0e0_wp]
  real(wp) :: phi_true(d) = [1.0e0_wp, 1.0e0_wp]
  real(wp) :: phi_min(d)
  real(wp) :: V_min

  real(wp) :: S
  integer :: i

  type(csv_file) :: f
  logical :: status_ok

  call create_datafile()

  clist = linspace(cmin, cmax, num)

  do i = 1, num

    c = clist(i)

    write(*,*)
    write(*,*)
    write(*,*) "RUN: ", i
    write(*,*) "  c =", c
    write(*,*)

    call minimize(  &
      V, phi_true, phi_min, V_min, maxiter=10000, mode=1)
    phi_true = phi_min

    call minimize(  &
      V, phi_false, phi_min, V_min, maxiter=10000, mode=1)
    phi_false = phi_min

    pd1 = solver(  &
      d, V,  &
      phi_false, phi_true,  &
      deform_eps=2.0e-2_wp,  &
      num_spline_knots=50,  &
      max_iter=100,  &
      verbose_level=1)

    call write_line(  &
      1, c, pd1%S, pd1%Sp, pd1%Sk, pd1%S0)

    pd2 = solver(  &
      d, V,  &
      phi_false, phi_true,  &
      deform_eps=5.0e-2_wp,  &
      num_spline_knots=20,  &
      max_iter=100,  &
      verbose_level=1)

    call write_line(  &
      2, c, pd2%S, pd2%Sp, pd2%Sk, pd2%S0)

  end do

  call f%close(status_ok)

contains

  function V(x) result(y)

    real(wp), intent(in) :: x(:)
    real(wp) :: y

    y = (x(1) ** 2 + x(2) ** 2) * (  &
      a * (x(1) - 1.0e0_wp) ** 2 +  &
      b * (x(2) - 1.0e0_wp) ** 2 - c)

  end function V

  subroutine create_datafile()

    call f%initialize(verbose=.true.)
    call f%open(  &
      "plots/pdm/thick_walled/data.csv",  &
      n_cols=6,  &
      status_ok=status_ok)
    call f%add(["id", "c ", "S ", "Sp", "Sk", "S0"])
    call f%next_row()

  end subroutine create_datafile

  subroutine write_line(id, c, S, Sp, Sk, S0)

    integer, intent(in) :: id
    real(wp), intent(in) :: c
    real(wp), intent(in) :: S
    real(wp), intent(in) :: Sp
    real(wp), intent(in) :: Sk
    real(wp), intent(in) :: S0

    call f%add(  &
      [real(id, wp), c, S, Sp, Sk, S0],  &
      real_fmt="(es15.5)")
    call f%next_row()

    write(*,*)

  end subroutine write_line

end program bouncesolver__app_pdm_thickwalled
