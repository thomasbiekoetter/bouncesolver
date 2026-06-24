program bouncesolver__test_pdm_thickwalled_initpaths
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
  real(wp), parameter :: c = 0.1e0_wp

  type(solver) :: pd
  integer, parameter :: d = 2
  real(wp) :: phi_false(d) = [0.0e0_wp, 0.0e0_wp]
  real(wp) :: phi_true(d) = [1.0e0_wp, 1.0e0_wp]

  call get_minima(phi_false, phi_true)

  call export_potential(  &
    [-0.5e0_wp, 1.5e0_wp],  &
    [-0.5e0_wp, 2.0e0_wp],  &
    200)

  call calc_bounce(maxiter=1)

  call calc_bounce(maxiter=40)

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

      call minimize(  &
        V, phi_true, phi_min, V_min,  &
        maxiter=10000, mode=1)

      phi_true = phi_min

      call minimize(  &
        V, phi_false, phi_min, V_min,  &
        maxiter=10000, mode=1)

      phi_false = phi_min

    end subroutine get_minima

    subroutine export_potential(  &
      phi1_lims, phi2_lims, number_points)

      real(wp), intent(in) :: phi1_lims(2)
      real(wp), intent(in) :: phi2_lims(2)
      integer, intent(in) :: number_points

      type(csv_file) :: f
      logical :: status_ok
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
        "plots/tests/thick_walled/initpaths/" //  &
          "pot_of_phi1_phi2.csv",  &
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

    subroutine calc_bounce(maxiter)

      integer, intent(in) :: maxiter

      character(len=2) :: sm
      character(len=100) :: prefix

      if (maxiter < 10) then
        write(sm, '(I1)') maxiter
      else
        write(sm, '(I2)') maxiter
      end if

      pd = solver(  &
        d, V,  &
        phi_false, phi_true,  &
        alpha=2,  &
        n_odeint=2000,  &
        verbose_level=1,  &
        smoothing=.true.,  &
        rho_max_fac=30.0e0_wp,  &
        init_path_mode=3,  &
        maxiter=maxiter)

      if (pd%exit_status >= 0) then
        prefix = "plots/tests/thick_walled/initpaths/" // sm
        call csv_x_of_rho(pd, trim(prefix) // "_x_of_rho.csv")
        call csv_pot_of_rho(pd, trim(prefix) // "_pot_of_rho.csv")
        call csv_pot_of_x(pd, trim(prefix) // "_pot_of_x.csv")
        call csv_phi_of_rho(pd, trim(prefix) // "_phi_of_rho.csv")
        call csv_phi_of_x(pd, trim(prefix) // "_phi_of_x.csv")
        call csv_forces_of_x(pd, trim(prefix) // "_forces_of_x.csv")
      end if

    end subroutine calc_bounce

end program bouncesolver__test_pdm_thickwalled_initpaths
