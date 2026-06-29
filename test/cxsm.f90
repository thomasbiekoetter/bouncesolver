program bouncesolver__test_cxsm

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
  use cxsm__potential_threedim, only : cxsm3D
  use cxsm__vacstruc_temperatures, only : get_Tcrit
  use cxsm__vacstruc_phases, only : phase
  use cxsm__vacstruc_temperatures, only : get_Tnuc
  use csv_module

  implicit none

  real(wp), parameter :: ms = 60.0e0_wp
  real(wp), parameter :: lams = 1.0e0_wp
  real(wp), parameter :: lamhs = 0.8e0_wp

  type(solver) :: pd
  type(cxsm3D) :: pot
  type(phase) :: ewsb_phase
  type(phase) :: singlet_phase
  real(wp) :: Tcrit
  real(wp) :: T
  real(wp) :: xT(2)
  real(wp) :: xF(2)
  real(wp) :: xmin(2)
  real(wp) :: Vmin
  real(wp) :: S
  real(wp) :: Tnuc

  pot = cxsm3D(ms, lams, lamhs)

  call get_Tcrit(  &
    pot, Tcrit,  &
    ewsb_phase=ewsb_phase,  &
    singlet_phase=singlet_phase)

  write(*,*)
  write(*,*) "Tcrit =", Tcrit
  write(*,*)

  T = 63.0e0_wp

  call export_potential(  &
    [-1.0e0_wp, 40.0e0_wp],  &
    [-1.0e0_wp, 40.0e0_wp],  &
    200)

  write(*,*) "Computing bounce at T =", T
  write(*,*)

  call ewsb_phase%get_x_at_T(T, xT)
  call minimize(V, xT, xmin, Vmin, maxiter=10000, mode=1)
  xT = xmin
  write(*,*) "  True minimum:"
  write(*,*) "       xT =", xT
  write(*,*) "    V(xT) =", V(xT)
  write(*,*)

  call singlet_phase%get_x_at_T(T, xF)
  call minimize(V, xF, xmin, Vmin, maxiter=10000, mode=1)
  xF = xmin
  write(*,*) "  False minimum:"
  write(*,*) "          xF =", xF
  write(*,*) "    V(xF, T) =", V(xF)
  write(*,*)

  write(*,*) "Computing bounce with straight initial path."
  pd = solver(  &
    2, V, xF, xT,  &
    rho_max_fac=80.0e0_wp,  &
    n_odeint=2000,  &
    deform_eps=2.0e-2_wp,  &
    Rforce_threshold=5.0e-2_wp,  &
    maxiter=40,  &
    init_path_mode=1,  &
    barrier_buffer=1.0e-5_wp,  &
    verbose_level=0)
  if (pd%exit_status >= 0) then
    call csv_x_of_rho(pd, "plots/tests/cxsm/straight/x_of_rho.csv")
    call csv_pot_of_rho(pd, "plots/tests/cxsm/straight/pot_of_rho.csv")
    call csv_pot_of_x(pd, "plots/tests/cxsm/straight/pot_of_x.csv")
    call csv_phi_of_rho(pd, "plots/tests/cxsm/straight/phi_of_rho.csv")
    call csv_phi_of_x(pd,  "plots/tests/cxsm/straight/phi_of_x.csv")
    call csv_forces_of_x(pd,  "plots/tests/cxsm/straight/forces_of_x.csv")
    write(*,*) "  -> S =", pd%S
  else
    write(*,*) "  FAILED"
  end if
  write(*,*)


  write(*,*) "Computing bounce with shallowest initial path."
  pd = solver(  &
    2, V, xF, xT,  &
    rho_max_fac=40.0e0_wp,  &
    n_odeint=2000,  &
    deform_eps=2.0e-3_wp,  &
    Rforce_threshold=5.0e-2_wp,  &
    maxiter=40,  &
    init_path_mode=3,  &
    barrier_buffer=1.0e-5_wp,  &
    verbose_level=0)
  if (pd%exit_status >= 0) then
    call csv_x_of_rho(pd, "plots/tests/cxsm/shallow/x_of_rho.csv")
    call csv_pot_of_rho(pd, "plots/tests/cxsm/shallow/pot_of_rho.csv")
    call csv_pot_of_x(pd, "plots/tests/cxsm/shallow/pot_of_x.csv")
    call csv_phi_of_rho(pd, "plots/tests/cxsm/shallow/phi_of_rho.csv")
    call csv_phi_of_x(pd,  "plots/tests/cxsm/shallow/phi_of_x.csv")
    call csv_forces_of_x(pd,  "plots/tests/cxsm/shallow/forces_of_x.csv")
    write(*,*) "  -> S =", pd%S
  else
    write(*,*) "  FAILED"
  end if
  write(*,*)

  write(*,*) "Computing bounce with straight initial path and no path deformation."
  pd = solver(  &
    2, V, xF, xT,  &
    rho_max_fac=40.0e0_wp,  &
    n_odeint=2000,  &
    deform_eps=2.0e-3_wp,  &
    Rforce_threshold=5.0e-2_wp,  &
    maxiter=1,  &
    init_path_mode=1,  &
    barrier_buffer=1.0e-5_wp,  &
    verbose_level=0)
  if (pd%exit_status >= 0) then
    write(*,*) "  -> S =", pd%S
  else
    write(*,*) "  FAILED"
  end if
  write(*,*)

  write(*,*) "Computing bounce with shallowest initial path and no path deformation."
  pd = solver(  &
    2, V, xF, xT,  &
    rho_max_fac=40.0e0_wp,  &
    n_odeint=2000,  &
    deform_eps=2.0e-3_wp,  &
    Rforce_threshold=5.0e-2_wp,  &
    maxiter=1,  &
    init_path_mode=3,  &
    barrier_buffer=1.0e-5_wp,  &
    verbose_level=0)
  if (pd%exit_status >= 0) then
    write(*,*) "  -> S =", pd%S
  else
    write(*,*) "  FAILED"
  end if
  write(*,*)

  write(*,*) "Finding nucleation tempreature:"
  call get_Tnuc(  &
    pot, Tnuc,  &
    S=S,  &
    Tcrit=Tcrit,  &
    ewsb_phase=ewsb_phase,  &
    singlet_phase=singlet_phase)
  write(*,*) "  -> Tnuc =", Tnuc

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
      "plots/tests/cxsm/pot_of_phi1_phi2.csv",  &
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

end program bouncesolver__test_cxsm
