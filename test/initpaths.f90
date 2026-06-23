program bouncesolver__test_initpaths
  use gradmin__descent, only : minimize
  use bouncesolver__config, only : wp
  use bouncesolver__config, only : fmt
  use bouncesolver__pdm, only : solver
  use bouncesolver__export, only : csv_phi_of_rho
  use bouncesolver__util, only : linspace
  use cxsm__potential_threedim, only : cxsm3D
  use cxsm__vacstruc_phases, only: phase
  use cxsm__vacstruc_phases, only: get_phases
  use csv_module, only : csv_file
  implicit none

  real(wp), parameter :: ms = 55.0e0_wp
  real(wp), parameter :: lams = 1.0e0_wp
  real(wp), parameter :: lamhs = 0.8e0_wp
  real(wp), parameter :: T = 60.0e0_wp

  type(solver) :: pd
  type(cxsm3D) :: pot
  real(wp) :: xF(2)
  real(wp) :: xT(2)

  pot = cxsm3D(ms, lams, lamhs)

  call export_potential(  &
    [-1.0e0_wp, 40.0e0_wp],  &
    [-1.0e0_wp, 40.0e0_wp],  &
    200)

  call get_minima(xF, xT)
  call calc_bounce_apprx(1)
  call calc_bounce_apprx(2)
  call calc_bounce_apprx(3)
  call calc_bounce(1)
  call calc_bounce(3)

  contains

    function V(phi) result(y)

      real(wp), intent(in) :: phi(:)
      real(wp) :: y

      y = pot%V(phi, T)

    end function V

    subroutine export_potential(  &
      phi1_lims, phi2_lims, number_points)

      real(wp), intent(in) :: phi1_lims(2)
      real(wp), intent(in) :: phi2_lims(2)
      integer, intent(in) :: number_points

      real(wp), allocatable :: phi1(:)
      real(wp), allocatable :: phi2(:)
      type(csv_file) :: f
      logical :: status_ok
      integer :: n
      integer :: i
      integer :: j

      n = number_points

      phi1 = linspace(phi1_lims(1), phi1_lims(2), n)
      phi2 = linspace(phi2_lims(1), phi2_lims(2), n)

      call f%initialize(verbose=.true.)

      call f%open(  &
        "plots/tests/initpaths/pot_of_phi1_phi2.csv",  &
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

    subroutine get_minima(xF, xT)

      real(wp), intent(out) :: xF(2)
      real(wp), intent(out) :: xT(2)

      type(phase) :: ewsb_phase
      type(phase) :: singlet_phase
      real(wp) :: xmin(2)
      real(wp) :: Vmin

      call get_phases(pot, ewsb_phase, singlet_phase)
      call ewsb_phase%get_x_at_T(T, xT)
      call minimize(V, xT, xmin, Vmin, maxiter=10000, mode=1)
      xT = xmin
      call singlet_phase%get_x_at_T(T, xF)
      call minimize(V, xF, xmin, Vmin, maxiter=10000, mode=1)
      xF = xmin

      write(*, *) "  True minimum: T =", T
      write(*, *) "       xT =", xT
      write(*, *) "    V(xT) =", V(xT)
      write(*, *)

      write(*, *) "  False minimum: T =", T
      write(*, *) "          xF =", xF
      write(*, *) "    V(xF, T) =", V(xF)
      write(*, *)

    end subroutine get_minima

    subroutine calc_bounce_apprx(init_path_mode)

      integer, intent(in) :: init_path_mode

      character(len=1) :: sf
      character(len=2) :: si

      write(*, *)
      write(*, *) "Compute bounce without path deformation:"

      write(sf, '(I1)') init_path_mode

      pd = solver(  &
        2, V, xF, xT,  &
        rho_max_fac=40.0e0_wp,  &
        n_odeint=2000,  &
        deform_eps=2.0e-2_wp,  &
        Rforce_threshold=5.0e-2_wp,  &
        maxiter=1,  &
        init_path_mode=init_path_mode,  &
        barrier_buffer=1.0e-5_wp,  &
        verbose_level=1)

      if (pd%exit_status >= 0)  &
        call csv_phi_of_rho(pd,  &
          "plots/tests/initpaths/phi_of_rho_" // sf // "_apprx.csv")

      write(*, *)

    end subroutine calc_bounce_apprx

    subroutine calc_bounce(init_path_mode)

      integer, intent(in) :: init_path_mode

      character(len=1) :: sf

      write(*, *)
      write(*, *) "Compute bounce with deformation:"

      write(sf, '(I1)') init_path_mode

      pd = solver(  &
        2, V, xF, xT,  &
        rho_max_fac=40.0e0_wp,  &
        n_odeint=2000,  &
        deform_eps=2.0e-2_wp,  &
        Rforce_threshold=5.0e-2_wp,  &
        maxiter=40,  &
        init_path_mode=init_path_mode,  &
        barrier_buffer=1.0e-5_wp,  &
        verbose_level=1)

      if (pd%exit_status >= 0)  &
        call csv_phi_of_rho(pd,  &
          "plots/tests/initpaths/phi_of_rho_" // sf // "_deform.csv")

      write(*, *)

    end subroutine calc_bounce

end program bouncesolver__test_initpaths
