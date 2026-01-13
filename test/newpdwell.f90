program bouncesolver__test_newpdwell

  use bouncesolver__config, only : wp
  use bouncesolver__pathdeformation, only : solver
  use gradmin__descent, only : minimize
  use bouncesolver__util, only : linspace
  use csv_module

  implicit none

  real(wp), parameter :: lam = 1.0e0_wp
  real(wp), parameter :: vev = 1.0e0_wp
  real(wp), parameter :: eps = 0.03e0_wp

  type(solver) :: pd
  integer, parameter :: num_fields = 1
  real(wp) :: phi_false(num_fields) = [1.0e0_wp]
  real(wp) :: phi_true(num_fields) = [-1.0e0_wp]

  integer :: i
  real(wp), allocatable :: x(:)
  real(wp) :: phimin(num_fields)
  real(wp) :: Vmin

  type(csv_file) :: f
  logical :: status_ok

  real(wp), allocatable :: rho(:)
  real(wp), allocatable :: pot(:)

  ! Get true minimum
  call minimize(V, [phi_true], phimin, Vmin, mode=1, maxiter=10000)
  write(*,*) "phi(true) =", phimin
  write(*,*) "V(true) =", Vmin
  write(*,*)
  phi_true = phimin

  ! Get false minimum
  call minimize(V, [phi_false], phimin, Vmin, mode=1, maxiter=10000)
  write(*,*) "x(false) =", phimin
  write(*,*) "V(false) =", Vmin
  write(*,*)
  phi_false = phimin

  pd = solver(  &
    num_fields,  &
    V,  &
    phi_false,  &
    phi_true,  &
    alpha=2,  &
    rho_max_fac=35.0e0_wp,  &
    xmin=1.0e-20_wp)

! x = linspace(pd%x_min, pd%x_max, 40)

! x(1) = pd%dV_dx(1.0e-8_wp)
! write(*,*) x(1)
! call exit

! write(*,*) "x        phi        V        dV        d2V"
! do i = 1, 40
!   write(*,*) x(i), pd%phi_of_x(x(i)), pd%V_of_x(x(i)), pd%dV_dx(x(i)), pd%d2V_dx2(x(i))
! end do

! write(*,*) pd%x_barrier

! write(*,*) pd%x_min, pd%x_max, maxval(pd%x_bounce(:))
  call f%initialize(verbose=.true.)
  call f%open(  &
    "x_of_rho.csv",  &
    n_cols=4,  &
    status_ok=status_ok)
  call f%add(["rho", "x1 ", "x2 ", "x1i"])
  call f%next_row()
  do i = 1, pd%nsteps_odeint
    call f%add([  &
      pd%rho_bounce(i),  &
      pd%x_bounce(i),  &
      pd%xdot_bounce(i),  &
      pd%x_of_rho(pd%rho_bounce(i))],  &
      real_fmt="(F15.10)")
    call f%next_row()
  end do
  call f%close(status_ok)

  rho = linspace(pd%rho_min, pd%rho_max, 100)
  allocate(pot(size(rho)))
  do i = 1, size(rho)
    pot(i) = pd%V_of_rho(rho(i))
  end do
  call f%initialize(verbose=.true.)
  call f%open(  &
    "pot_of_rho.csv",  &
    n_cols=2,  &
    status_ok=status_ok)
  call f%add(["rho", "pot"])
  call f%next_row()
  do i = 1, size(rho)
    call f%add([rho(i), pot(i)], real_fmt="(F15.10)")
    call f%next_row()
  end do
  call f%close(status_ok)

  write(*,*) "Action = ", pd%S

contains

  function V(x) result(y)

    real(wp), intent(in) :: x(:)
    real(wp) :: y

    y = lam * (x(1) ** 2 - vev ** 2) ** 2 / 4.0e0_wp + eps * x(1)

  end function V

end program bouncesolver__test_newpdwell
