program bouncesolver__test_pdm_thickwalledstraight

  use bouncesolver__config, only : wp
  use bouncesolver__config, only : fmt
  use bouncesolver__pdm, only : solver
  use bouncesolver__util, only : linspace
  use bouncesolver__export, only : csv_x_of_rho
  use bouncesolver__export, only : csv_pot_of_rho
  use bouncesolver__export, only : csv_pot_of_x
  use bouncesolver__export, only : csv_phi_of_rho
  use bouncesolver__export, only : csv_phi_of_x
  use bouncesolver__export, only : csv_forces_of_x
  use gradmin__descent, only : minimize
  use csv_module

  implicit none

  real(wp), parameter :: a = 1.8e0_wp
  real(wp), parameter :: b = 0.2e0_wp
  real(wp), parameter :: c = 0.06e0_wp

  type(solver) :: pd
  integer, parameter :: d = 2
  real(wp) :: phi_false(d) = [0.0e0_wp, 0.0e0_wp]
  real(wp) :: phi_true(d) = [1.0e0_wp, 1.0e0_wp]

  type(csv_file) :: f
  logical :: status_ok

  call get_minima(phi_false, phi_true)

  call export_potential(  &
    [-0.5e0_wp, 1.5e0_wp],  &
    [-0.5e0_wp, 2.0e0_wp],  &
    200)

  pd = solver(  &
    d, V,  &
    phi_false, phi_true,  &
    alpha=2,  &
    n_odeint=2000,  &
    verbose_level=1,  &
!   deform_eps=1.0e-1_wp,  &
    smoothing=.true.,  &
    rho_max_fac=40.0e0_wp,  &
    init_path_mode=1,  &
    maxiter=40)

  call pd%print_exit_status()

  call csv_x_of_rho(pd, "plots/tests/thick_walled/straight/x_of_rho.csv")

  call csv_pot_of_rho(pd, "plots/tests/thick_walled/straight/pot_of_rho.csv")

  call csv_pot_of_x(pd, "plots/tests/thick_walled/straight/pot_of_x.csv")

  call csv_phi_of_rho(pd, "plots/tests/thick_walled/straight/phi_of_rho.csv")

  call csv_phi_of_x(pd,  "plots/tests/thick_walled/straight/phi_of_x.csv")

  call csv_forces_of_x(pd,  "plots/tests/thick_walled/straight/forces_of_x.csv")

contains

  function V(phi) result(y)

    real(wp), intent(in) :: phi(:)
    real(wp) :: y

    y = (phi(1) ** 2 + phi(2) ** 2) * (  &
      a * (phi(1) - 1.0e0_wp) ** 2 +  &
      b * (phi(2) - 1.0e0_wp) ** 2 - c)

  end function V

  subroutine get_minima(phi_false, phi_true)

    real(wp), intent(inout) :: phi_false(d)
    real(wp), intent(inout) :: phi_true(d)

    real(wp) :: phi_min(d)
    real(wp) :: V_min

    call minimize(V, phi_true, phi_min, V_min, maxiter=10000, mode=1)
    phi_true = phi_min

    call minimize(V, phi_false, phi_min, V_min, maxiter=10000, mode=1)
    phi_false = phi_min

!   write(*,*) phi_true
!   write(*,*) phi_false
!   call exit

  end subroutine get_minima

  subroutine export_potential(  &
    phi1_lims, phi2_lims, number_points)

    real(wp), intent(in) :: phi1_lims(2)
    real(wp), intent(in) :: phi2_lims(2)
    integer, intent(in) :: number_points

    real(wp), allocatable :: phi1(:)
    real(wp), allocatable :: phi2(:)

    integer :: n
    integer :: i
    integer :: j

    n = number_points

    phi1 = linspace(phi1_lims(1), phi1_lims(2), n)
    phi2 = linspace(phi2_lims(1), phi2_lims(2), n)

    call f%initialize(verbose=.true.)

    call f%open(  &
      "plots/tests/thick_walled/straight/pot_of_phi1_phi2.csv",  &
      n_cols=3,  &
      status_ok=status_ok)

    call f%add(["phi1", "phi2", "pot "])
    call f%next_row()

    do i = 1, n
      do j = 1, n
        call f%add(  &
          [phi1(i), phi2(j), V([phi1(i), phi2(j)])],  &
          real_fmt=fmt)
        call f%next_row()
      end do
    end do

    call f%close(status_ok)

  end subroutine export_potential

end program bouncesolver__test_pdm_thickwalledstraight
