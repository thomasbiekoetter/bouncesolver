module bouncesolver__pdm

  use bouncesolver__config, only : wp
  use bouncesolver__util, only : linspace
  use bouncesolver__util, only : pi
  use bouncesolver__util, only : riemann_integrate
  use bouncesolver__util, only : is_equal
  use odeint__rk4, only : integrate
  use gradmin__derivatives, only : gradient
  use gradmin__descent, only : minimize
  use bspline_module, only : bspline_1d
  use bspline_module, only : get_status_message
  use cubicsplines__smoothing, only : add_viscosity

  implicit none

  private

  real(wp), parameter :: eps_gradient = 1.0e-6_wp
  real(wp), parameter :: eps_hessian = 1.0e-4_wp

  integer, parameter :: alpha_default = 2
  real(wp), parameter :: x_min_start = 1.0e-18_wp
  real(wp), parameter :: rho_max_fac_default = 10.0e0_wp

  abstract interface
    function V_abstract(x) result(y)
      import wp
      real(wp), intent(in) :: x(:)
      real(wp) :: y
    end function V_abstract
  end interface

  type, public :: solver
    procedure(V_abstract), pointer, nopass :: V => null()
    integer :: d
    integer :: alpha
    integer :: n = 1000
    integer :: n_spls = 10
    integer :: maximum_iterations = 100
    integer :: verbose_level = 1
    real(wp) :: Rforce_threshold = 0.1e0_wp
    real(wp) :: deform_eps = 2.0e-2_wp
    real(wp) :: Rforce
    real(wp) :: Rforce_previous = 1.0e10_wp
    real(wp) :: Rforce_initial = 1.0e10_wp
    real(wp) :: Rforce_best = 1.0e10_wp
    real(wp), allocatable :: phi_false(:)
    real(wp), allocatable :: phi_true(:)
    real(wp), allocatable :: x(:)
    real(wp), allocatable :: xdot(:)
    real(wp), allocatable :: phi(:, :)
    real(wp), allocatable :: rho(:)
    real(wp), allocatable :: pot(:)
    real(wp), allocatable :: xb(:)
    real(wp), allocatable :: xbdot(:)
    real(wp) :: x_min
    real(wp) :: rho_min
    real(wp) :: rho_max
    real(wp) :: rho_max_fac
    integer :: ib_max
    real(wp) :: S
    real(wp) :: Sp
    real(wp) :: Sk
    real(wp) :: S0
    real(wp) :: S_best = 1.0e10_wp
    real(wp) :: Sp_best
    real(wp) :: Sk_best
    real(wp), allocatable :: Nforce(:, :)
    real(wp), allocatable :: Pforce(:, :)
    real(wp), allocatable :: xforce(:)
    real(wp) :: grad_V_max_on_path
    integer :: current_iteration
    integer, allocatable :: ixs(:)
    real(wp) :: gradV_max
    real(wp), allocatable :: phi_deformed(:, :)
    type(bspline_1d), allocatable :: bspls(:)
    integer, allocatable :: iflag_bspls(:)
    integer :: kx_bspls
    logical :: flag_too_thin = .false.
    logical :: flag_trivial_solution = .false.
  contains
    procedure, private :: allocate_arrays
    procedure, private :: construct_starting_path
    procedure, private :: check_minima
    procedure, private :: phi_of_x
    procedure, public :: dphi_dx
    procedure, public :: d2phi_dx2
    procedure, private :: V_of_x
    procedure, private :: find_x_barrier
    procedure, private :: d2V_dx2
    procedure, private :: estimate_rho_max
    procedure, private :: bounce_on_path
    procedure, private :: dV_dx
    procedure, private :: clip_bounce
    procedure, private :: calc_action
    procedure, private :: calc_forces
    procedure, private :: gradV
    procedure, private :: deform_path
    procedure, public :: check_convergence
  end type solver

  interface solver
    module procedure create_solver
  end interface solver

contains

  function create_solver(  &
    d, V, phi_false, phi_true,  &
    alpha, rho_max_fac, x_min,  &
    deform_eps, max_iter,  &
    num_odeint,  &
    num_spline_knots,  &
    verbose_level) result(this)

    integer, intent(in) :: d
    procedure(V_abstract) :: V
    real(wp), intent(in) :: phi_false(d)
    real(wp), intent(in) :: phi_true(d)
    integer, intent(in), optional :: alpha
    real(wp), intent(in), optional :: rho_max_fac
    real(wp), intent(in), optional :: x_min
    real(wp), intent(in), optional :: deform_eps
    integer, intent(in), optional :: max_iter
    integer, intent(in), optional :: num_odeint
    integer, intent(in), optional :: num_spline_knots
    integer, intent(in), optional :: verbose_level
    type(solver) :: this

    logical :: converged
    integer :: iteration

    this%d = d
    this%V => V
    this%phi_false = phi_false
    this%phi_true = phi_true
    if (present(alpha)) then
      this%alpha = alpha
    else
      this%alpha = alpha_default
    end if
    if (present(rho_max_fac)) then
      this%rho_max_fac = rho_max_fac
    else
      this%rho_max_fac = rho_max_fac_default
    end if
    if (present(x_min)) then
      this%x_min = x_min
    else
      this%x_min = x_min_start
    end if

    if (present(verbose_level)) then
      this%verbose_level = verbose_level
    end if

    if (present(deform_eps)) then
      this%deform_eps = deform_eps
    end if

    if (present(max_iter)) then
      this%maximum_iterations = max_iter
    end if

    if (present(num_odeint)) then
      this%n = num_odeint
    end if

    if (present(num_spline_knots)) then
      this%n_spls = num_spline_knots
    end if

    this%kx_bspls = 5

    call this%check_minima()
    call this%allocate_arrays()
    call this%construct_starting_path()
    call this%estimate_rho_max()
    converged = .false.
    iteration = 0
    do while (.not. converged)
      iteration = iteration + 1
      call this%bounce_on_path()
      call this%calc_action()
      if (this%flag_too_thin) exit
      if (this%flag_trivial_solution) exit
      if (this%verbose_level >= 1) then
        write(*,*) "Action = ", this%S, this%Sp, this%Sk
      end if
      if (iteration == 1) this%S0 = this%S
      if (abs(this%S / this%Sp) > 1.0e2_wp) then
        write(*,*) "Bounce wrong. Exiting."
        exit
      end if
      call this%calc_forces(iteration)
      call this%deform_path()
      call this%check_convergence(iteration, converged)
    end do

    this%S = this%S_best
    this%Sp = this%Sp_best
    this%Sk = this%Sk_best

    if (this%verbose_level >= 1) then
      write(*,*)
      write(*,*) "Result:"
      write(*,*) "Action = ", this%S, this%Sp, this%Sk
      write(*,*) "Force ratio =", this%Rforce_best
    end if

  end function create_solver

  subroutine allocate_arrays(this)

    class(solver), intent(inout) :: this

    integer :: n
    integer :: d

    n = this%n
    d = this%d

    allocate(this%phi(n, d))
    allocate(this%x(n))
    allocate(this%xdot(n))
    allocate(this%rho(n))
    allocate(this%pot(n))
    allocate(this%xb(n))
    allocate(this%xbdot(n))
    allocate(this%Nforce(n, d))
    allocate(this%Pforce(n, d))
    allocate(this%xforce(n))
    allocate(this%ixs(n))
    allocate(this%phi_deformed(n, d))
    allocate(this%bspls(d))
    allocate(this%iflag_bspls(d))

  end subroutine allocate_arrays

  subroutine check_minima(this)

    class(solver), intent(inout) :: this

    real(wp) :: phi_min(this%d)
    real(wp) :: V_min

    call minimize(  &
      this%V, this%phi_true, phi_min, V_min,  &
      maxiter=10000, mode=1)

    if (norm2(phi_min - this%phi_true) > 10 * eps_gradient) then
      write(*,*) "Warning: check true minimum."
      write(*,*) "Using refined phi_true =", phi_min
      write(*,*) "Given phi_true =", this%phi_true
      this%phi_true = phi_min
    end if

    call minimize(  &
      this%V, this%phi_false, phi_min, V_min,  &
      maxiter=10000, mode=1)

    if (norm2(phi_min - this%phi_false) > 10 * eps_gradient) then
      write(*,*) "Warning: check false minimum."
      write(*,*) "Using refined phi_false =", phi_min
      write(*,*) "Given phi_false =", this%phi_false
      this%phi_false = phi_min
    end if

  end subroutine check_minima

  subroutine construct_starting_path(this)

    class(solver), intent(inout) :: this

    real(wp) :: x(this%n)
    real(wp) :: phi(this%n, this%d)
    real(wp) :: dx(this%n - 1)
    real(wp) :: pot(this%n)
    integer :: i
    integer :: j
    integer :: n
    integer :: d
    real(wp) :: dphidx
    real(wp), allocatable :: x_spline(:)
    real(wp), allocatable :: phi_spline(:, :)
    integer :: idx_bspls

    n = this%n
    d = this%d

    x = linspace(0.0e0_wp, 1.0e0_wp, n)

    do i = 1, n
      phi(i, :) = this%phi_true +  &
        x(i) * (this%phi_false - this%phi_true)
    end do

    do i = 1, n - 1
      dx(i) = norm2(phi(i + 1, :) - phi(i, :))
    end do
    do i = 2, n
      x(i) = x(i - 1) + dx(i - 1)
    end do

    do i = 1, n - 1
      dphidx = norm2(phi(i + 1, :) - phi(i, :)) /  &
        (x(i + 1) - x(i))
      if (.not. is_equal(abs(dphidx), 1.0e0_wp)) then
        write(*,*) "init path went wrong"
        call exit
      end if
    end do

    allocate(x_spline(this%n_spls))
    allocate(phi_spline(this%n_spls, d))

    do i = 1, this%n_spls
      if (i == 1) then
        x_spline(i) = x(i)
        phi_spline(i, :) = phi(i, :)
      else if (i == this%n_spls) then
        x_spline(i) = x(this%n)
        phi_spline(i, :) = this%phi_false
      else
        j = 1 + (i - 1) * n / this%n_spls
        x_spline(i) = x(j)
        phi_spline(i, :) = phi(j, :)
      end if
    end do

    do j = 1, d
      call this%bspls(j)%initialize(  &
        x_spline,  &
        phi_spline(:, j), this%kx_bspls, this%iflag_bspls(j))
      if (this%iflag_bspls(j) /= 0) then
        write(*,*) 'Error initializing ', j, 'D spline: ' //  &
          get_status_message(this%iflag_bspls(j))
        call exit
      end if
    end do

    idx_bspls = 0
    do j = 1, d
      do i = 1, n
        call this%bspls(j)%evaluate(  &
          x(i), idx_bspls, phi(i, j), this%iflag_bspls(j))
      end do
    end do

    do i = 1, n
      pot(i) = this%V(phi(i, :))
    end do

    this%x = x
    this%phi = phi
    this%pot = pot

  end subroutine construct_starting_path

  function phi_of_x(this, x) result(phi)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: x
    real(wp) :: phi(this%d)

    integer :: j
    integer :: idx_bspls

    idx_bspls = 0

    do j = 1, this%d
      call this%bspls(j)%evaluate(  &
        x, idx_bspls, phi(j), this%iflag_bspls(j))
    end do

  end function phi_of_x

  function dphi_dx(this, x) result(dphi)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: x
    real(wp) :: dphi(this%d)

    integer :: j
    integer :: idx_bspls

    idx_bspls = 1

    do j = 1, this%d
      call this%bspls(j)%evaluate(  &
        x, idx_bspls, dphi(j), this%iflag_bspls(j))
    end do

  end function dphi_dx

  function d2phi_dx2(this, x) result(d2phi)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: x
    real(wp) :: d2phi(this%d)

    integer :: j
    integer :: idx_bspls

    idx_bspls = 2

    do j = 1, this%d
      call this%bspls(j)%evaluate(  &
        x, idx_bspls, d2phi(j), this%iflag_bspls(j))
    end do

  end function d2phi_dx2

  function V_of_x(this, x) result(V)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: x
    real(wp) :: V

    real(wp) :: phi(this%d)

    phi = this%phi_of_x(x)
    V = this%V(phi)

  end function V_of_x

  function d2V_dx2(this, x) result(d2V)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: x
    real(wp) :: d2V

    real(wp) :: V_m2
    real(wp) :: V_m1
    real(wp) :: V
    real(wp) :: V_p1
    real(wp) :: V_p2
    real(wp) :: eps
    logical :: lower_ok
    logical :: upper_ok
    integer :: n
    integer :: i

    n = this%n

    V = abs(this%V_of_x(x))
    eps = maxval([eps_hessian * V, eps_hessian])

    lower_ok = (x - 2.0e0_wp * eps) > this%x(1)
    upper_ok = (x + 2.0e0_wp * eps) < this%x(n)

    if (lower_ok .and. upper_ok) then
      ! central
      V_m2 = this%V_of_x(x - 2.0e0_wp * eps)
      V_m1 = this%V_of_x(x - 1.0e0_wp * eps)
      V = this%V_of_x(x)
      V_p1 = this%V_of_x(x + 1.0e0_wp * eps)
      V_p2 = this%V_of_x(x + 2.0e0_wp * eps)
      d2V = (-V_m2 + 16.0e0_wp * V_m1 -  &
        30.0e0_wp * V + 16.0e0_wp * V_p1 -  &
        V_p2) / (12.e0_wp * eps ** 2)
    else if (upper_ok) then
      ! forward
      V_m2 = this%V_of_x(x)
      V_m1 = this%V_of_x(x + 1.0e0_wp * eps)
      V = this%V_of_x(x + 2.0e0_wp * eps)
      V_p1 = this%V_of_x(x + 3.0e0_wp * eps)
      V_p2 = this%V_of_x(x + 4.0e0_wp * eps)
      d2V = (35.0e0_wp * V_m2 - 104.0e0_wp * V_m1 +  &
        114.0e0_wp * V - 56.0e0_wp * V_p1 +  &
        11.0e0_wp * V_p2) / (12.0e0_wp * eps ** 2)
    else if (lower_ok) then
      V_m2 = this%V_of_x(x - 4.0e0_wp * eps)
      V_m1 = this%V_of_x(x - 3.0e0_wp * eps)
      V = this%V_of_x(x - 2.0e0_wp * eps)
      V_p1 = this%V_of_x(x - 1.0e0_wp * eps)
      V_p2 = this%V_of_x(x)
      d2V = (35.0e0_wp * V_p2 - 104.0e0_wp * V_p1 +  &
        114.0e0_wp * V - 56.0e0_wp * V_m1 +  &
        11.0e0_wp * V_m2) / (12.0e0_wp * eps ** 2)
    else
      write(*,*) "Problem in d2V_dx2"
      call exit
    end if

  end function d2V_dx2

  function dV_dx(this, x) result(dV)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: x
    real(wp) :: dV

    real(wp) :: V_m2
    real(wp) :: V_m1
    real(wp) :: V
    real(wp) :: V_p1
    real(wp) :: V_p2
    real(wp) :: eps
    logical :: lower_ok
    logical :: upper_ok
    integer :: n
    integer :: i
!   real(wp) :: blend
!   real(wp) :: dV_blend
!   real(wp) :: delta_phi_true(this%d)

    n = this%n

    V = abs(this%V_of_x(x))
    eps = maxval([eps_gradient * V, eps_gradient])

    lower_ok = (x - 2.0e0_wp * eps) > this%x(1)
    upper_ok = (x + 2.0e0_wp * eps) < this%x(n)

    if (lower_ok .and. upper_ok) then
      ! central
      V_m2 = this%V_of_x(x - 2.0e0_wp * eps)
      V_m1 = this%V_of_x(x - 1.0e0_wp * eps)
      V_p1 = this%V_of_x(x + 1.0e0_wp * eps)
      V_p2 = this%V_of_x(x + 2.0e0_wp * eps)
      dV = (-V_p2 + 8.0e0_wp * V_p1 -  &
        8.0e0_wp * V_m1 + V_m2) / (12.0e0_wp * eps)
    else if (upper_ok) then
      ! forward
      V_m2 = this%V_of_x(x)
      V_m1 = this%V_of_x(x + 1.0e0_wp * eps)
      V_p1 = this%V_of_x(x + 2.0e0_wp * eps)
      V_p2 = this%V_of_x(x + 3.0e0_wp * eps)
      dV = (-11.0e0_wp * V_m2 + 18.0e0_wp * V_m1 -  &
        9.0e0_wp * V_p1 + 2.0e0_wp * V_p2) /  &
        (6.0e0_wp * eps)
    else if (lower_ok) then
      V_m2 = this%V_of_x(x - 4.0e0_wp * eps)
      V_m1 = this%V_of_x(x - 3.0e0_wp * eps)
      V = this%V_of_x(x - 2.0e0_wp * eps)
      V_p1 = this%V_of_x(x - 1.0e0_wp * eps)
      V_p2 = this%V_of_x(x)
      dV = (25.0e0_wp * V_p2 - 48.0e0_wp * V_p1 +  &
        36.0e0_wp * V - 16.0e0_wp * V_m1 +  &
        3.0e0_wp * V_m2) / (12.0e0_wp * eps)
    else
      dV = 0.0e0_wp
      this%flag_trivial_solution = .true.
    end if

!   delta_phi_true = this%phi_of_x(x) - this%phi_true
!   if (norm2(delta_phi_true) < eps) then
!     dV_blend = this%d2V_dx2(x) * (x - this%x(1))
!     write(*,*) dV, dV_blend
!     blend = exp(-(norm2(delta_phi_true) / eps) ** 2)
!     dV = dV * (1.0e0_wp - blend) + dV_blend * blend
!   end if

  end function dV_dx

  subroutine estimate_rho_max(this)

    class(solver), intent(inout) :: this
    real(wp) :: msq_false

    msq_false = this%d2V_dx2(this%x(this%n))

    this%rho_max = this%rho_max_fac * sqrt(msq_false)

  end subroutine estimate_rho_max

  subroutine bounce_on_path(this)

    class(solver), intent(inout) :: this

    logical :: shooting_converged
    real(wp) :: x_min
    real(wp) :: x_max
    real(wp) :: rho_min
    real(wp) :: rho_max
    real(wp) :: x1
    real(wp) :: x0(2)
    integer :: over_under_flag
    logical :: overshot
    integer :: n
    real(wp), allocatable :: rho(:)
    real(wp), allocatable :: x(:, :)
    integer :: i

    real(wp) :: temp
    real(wp) :: temp2(this%d)

    n = this%n
    x_min = this%x_min
    x_max = this%find_x_barrier()

    rho_min = 1.0e-4_wp
    rho_max = this%rho_max

    over_under_flag = 0
    shooting_converged = .false.

    do while (.not. shooting_converged)

      select case (over_under_flag)
        case (0) ! First iteration
          x1 = (x_max + x_min) / 100.0e0_wp
        case (1) ! Overshoot
          x_min = x1
          x1 = (x_max + x_min) / 2.0e0_wp
        case (-1) ! Undershoot
          x_max = x1
          x1 = (x_max + x_min) / 2.0e0_wp
      end select

      x0(1) = x1
      x0(2) = 0.0e0_wp

      call integrate(  &
        dxdrho,  &
        x0,  &
        rho_min,  &
        rho_max,  &
        n,  &
        rho, x)

      if (overshot) then
        over_under_flag = 1
      else
        over_under_flag = -1
      end if

      if (is_equal(  &
            x_min / x_max,  &
            1.0e0_wp,  &
            eps=1.0e-10_wp)) then
        shooting_converged = .true.
      else
        shooting_converged = .false.
      end if

      if (this%verbose_level >= 2) then
        write(*,*) x_min, x_max, over_under_flag
      end if

    end do

    if (is_equal(x_max / this%x_min, 1.0e0_wp, eps=1.0e-10_wp)) then
      this%flag_too_thin = .true.
    end if

    this%rho = rho
    this%xb = x(:, 1)
    this%xbdot = x(:, 2)
    this%x = this%xb
    do i = 1, n
      this%phi(i, :) = this%phi_of_x(this%x(i))
    end do
    do i = 1, n
      this%pot(i) = this%V_of_x(this%x(i))
    end do

    contains

      function dxdrho(x, rho) result(xdot)

        real(wp), intent(in) :: x(:)
        real(wp), intent(in) :: rho
        real(wp), allocatable :: xdot(:)

        integer :: ix

        allocate(xdot(2))

        if (x(1) > this%x(this%n)) then
          overshot = .true.
        else
          overshot = .false.
        end if

        if (.not. overshot) then
          xdot(1) = x(2)
          xdot(2) = this%dV_dx(x(1)) -  &
            this%alpha * x(2) / rho
        else
          xdot(1) = 0.0e0_wp
          xdot(2) = 0.0e0_wp
        end if

      end function dxdrho

  end subroutine bounce_on_path

  function find_x_barrier(this) result(x_barrier)

    class(solver), intent(inout) :: this
    real(wp) :: x_barrier

    real(wp) :: V1
    real(wp) :: V2
    integer :: i

    V1 = this%V(this%phi(1, :))
    ! Not start close to phi_true here
    do i = this%n / 10, this%n
      V2 = this%V(this%phi(i, :))
      if (V2 > V1) then
        V1 = V2
      else
        x_barrier = this%x(i - 1)
        exit
      end if
    end do

    if (i == this%n + 1) then
      write(*,*) "No potential barrier on path."
      call exit
    end if

  end function find_x_barrier

  subroutine clip_bounce(this, iclip)

    class(solver), intent(inout) :: this
    integer, intent(inout) :: iclip

    integer :: i
    integer :: n
    integer :: clip_range
    real(wp) :: xb(this%n)
    real(wp) :: xbdot(this%n)
    real(wp) :: x_max

    n = this%n
    xb = this%x
    xbdot = this%xbdot
    x_max = maxval(this%x)
    clip_range = 0

    do i = 1, n
      if (is_equal(xb(i) / x_max, 1.0e0_wp, eps=1.0e-3_wp)) then
        clip_range = clip_range + 1
      else
        clip_range = 0
      end if
      if (clip_range == 10) then
        if (this%verbose_level >= 2) then
          write(*,"(a,I5,a,ES12.3)")  &
            "Clipped solution" // &
            "at rho(",  i - 1, ") = ",  &
            this%rho(i - 1)
        end if
        exit
      end if
    end do

    iclip = minval([i - 1, n])
    this%ib_max = iclip

  end subroutine clip_bounce

  subroutine calc_action(this)

    class(solver), intent(inout) :: this

    real(wp) :: dimen
    integer :: n
    integer :: i
    real(wp), allocatable :: fpre(:)
    real(wp), allocatable :: f(:)
    real(wp), allocatable :: rho(:)
    real(wp), allocatable :: pot_term(:)
    real(wp), allocatable :: kin_term(:)
    real(wp) :: dphi
    real(wp) :: drho

    dimen = this%alpha + 1.0e0_wp

    call this%clip_bounce(n)
    allocate(rho(n))
    allocate(fpre(n))
    allocate(f(n))
    allocate(pot_term(n))
    allocate(kin_term(n))

    rho = this%rho(1:n)

    if (this%alpha == 2) then
      fpre = 4.0e0_wp * pi * rho ** 2
    else if (this%alpha == 3) then
      fpre = 2.0e0_wp * pi ** 2 * rho ** 3
    else
      write(*,*) "Wrong alpha value."
      call exit
    end if

    pot_term = this%pot(1:n) - this%pot(n)

    ! d phi / d rho = dphi / dx * dx / drho = dx / drho
    do i = 1, n - 1
      dphi = this%xb(i + 1) - this%xb(i)
      drho = rho(i + 1) - rho(i)
      kin_term(i) = (dphi / drho) ** 2 / 2.0e0_wp
    end do
    dphi = this%xb(n) - this%xb(n - 1)
    drho = rho(n) - rho(n - 1)
    kin_term(n) = (dphi / drho) ** 2 / 2.0e0_wp

    ! From potential term only
    f = fpre * 2.0e0_wp * pot_term / (2.0e0_wp - dimen)
    this%Sp = riemann_integrate(rho, f)

    ! From kinetic term only
    f = fpre * 2.0e0_wp * kin_term / dimen
    this%Sk = riemann_integrate(rho, f)

    ! All terms
    f = fpre * (kin_term + pot_term)
    this%S = riemann_integrate(rho, f)

  end subroutine calc_action

  function gradV(this, phi) result(dV)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: phi(this%d)
    real(wp) :: dV(this%d)

    real(wp) :: V
    real(wp) :: h

    V = abs(this%V(phi))
    h = maxval([V * eps_gradient, eps_gradient])

    dV = gradient(this%V, phi, eps=h)

  end function gradV

  subroutine calc_forces(this, iteration)

    class(solver), intent(inout) :: this
    integer, intent(in) :: iteration

    integer :: i
    integer :: j
    integer :: n
    real(wp) :: rho(this%ib_max)
    real(wp) :: dx_drho(this%ib_max)
    real(wp) :: phi(this%ib_max, this%d)
    real(wp) :: x(this%ib_max)
    real(wp) :: d2phi_dx2(this%ib_max, this%d)
    real(wp) :: gradV(this%ib_max, this%d)
    real(wp) :: pathdir(this%ib_max, this%d)
    real(wp) :: gradV_perp(this%ib_max, this%d)
    real(wp) :: gradV_para(this%ib_max, this%d)
    real(wp) :: Nforce(this%ib_max, this%d)
    real(wp) :: Pforce(this%ib_max, this%d)
    real(wp) :: gradV_len(this%ib_max)

    n = this%ib_max
    rho = this%rho(1:n)
    dx_drho = this%xbdot(1:n)

    x = this%x(1:n)
    phi = this%phi(1:n, :)

    if (iteration == 0) then
      d2phi_dx2 = 0.0e0_wp
      dx_drho = 0.0e0_wp
    else
      do i = 1, n
        d2phi_dx2(i, :) = this%d2phi_dx2(x(i))
      end do
    end if

    do i = 1, n
      gradV(i, :) = this%gradV(phi(i, :))
    end do

    do i = 1, n
      gradV_len(i) = norm2(gradV(i, :))
    end do
    this%gradV_max = maxval(gradV_len)

    do i = 1, n
      pathdir(i, :) = this%dphi_dx(x(i))
      if (this%verbose_level >= 3) then
        if (.not. is_equal(norm2(pathdir(i, :)), 1.0e0_wp, eps=1.0e-2_wp)) then
          write(*,*) "Problem with pathdir in calc_forces."
          write(*,*) i, pathdir(i, :), norm2(pathdir(i, :))
        end if
      end if
    end do

    do i = 1, n
      gradV_perp(i, :) = gradV(i, :) -  &
        dot_product(gradV(i, :), pathdir(i, :)) *  &
        pathdir(i, :)
    end do

    do i = 1, n
      gradV_para(i, :) =  &
        dot_product(gradV(i, :), pathdir(i, :)) *  &
        pathdir(i, :)
    end do

    do i = 1, n
      Nforce(i, :) = d2phi_dx2(i, :) * dx_drho(i) ** 2 -  &
        gradV_perp(i, :)
    end do

    do i = 1, n
      Pforce(i, :) = gradV_para(i, :)
    end do

    this%xforce = 0.0e0_wp
    this%xforce(1:n) = x

    this%Nforce = 0.0e0_wp
    this%Nforce(1:n, :) = Nforce

    this%Pforce = 0.0e0_wp
    this%Pforce(1:n, :) = Pforce

    contains

      function dT_step(Tl, T0, Tu, xl, x0, xu) result(y)

        real(wp), intent(in) :: Tl
        real(wp), intent(in) :: T0
        real(wp), intent(in) :: Tu
        real(wp), intent(in) :: xl
        real(wp), intent(in) :: x0
        real(wp), intent(in) :: xu
        real(wp) :: y

        y = 2.0e0_wp * (  &
          (Tu - T0) / (xu - x0) -  &
          (T0 - Tl) / (x0 - xl)) / (xu - xl)

      end function dT_step

  end subroutine calc_forces

  subroutine deform_path(this)

    class(solver), intent(inout) :: this

    integer :: i
    integer :: j
    integer :: n
    integer :: d
    real(wp) :: phi_old(this%d)
    real(wp) :: phi_new(this%d)
    real(wp) :: phi_deformed(this%ib_max, this%d)
    real(wp) :: rescale
    real(wp) :: eps
    real(wp), allocatable :: x_spline(:)
    real(wp), allocatable :: phi_spline(:, :)
    integer :: idx_bspls
    real(wp) :: phi(this%n, this%d)
    real(wp) :: x(this%n)
    real(wp) :: dx(this%n)
    real(wp) :: dphidx

    d = this%d

    n = this%ib_max

    rescale = norm2(this%phi_false - this%phi_true)
    rescale = rescale / this%gradV_max

    eps = this%deform_eps * rescale

    do i = 1, n
      phi_old = this%phi(i, :)
      phi_new = phi_old + eps * this%Nforce(i, :)
      phi_deformed(i, :) = phi_new
    end do

    allocate(x_spline(this%n_spls))
    allocate(phi_spline(this%n_spls, d))

    x_spline = linspace(0.0e0_wp, 1.0e0_wp, this%n_spls)
    do i = 1, this%n_spls
      if (i == 1) then
        phi_spline(i, :) = phi_deformed(1, :)
      else if (i == this%n_spls) then
        phi_spline(i, :) = this%phi_false
      else
        j = 1 + (i - 1) * n / this%n_spls
        phi_spline(i, :) = phi_deformed(j, :)
      end if
    end do

    do j = 1, d
      phi_spline(:, j) = add_viscosity(  &
        x_spline,  &
        phi_spline(:, j),  &
        this%n_spls,  &
        eps=2.0e0_wp,  &
        iterations=100)
    end do

    ! Construct B-spline approximation to get path with
    ! this%n points and for more precise calculation
    ! of d2phi_dx2
    do j = 1, d
      call this%bspls(j)%initialize(  &
        x_spline,  &
        phi_spline(:, j), this%kx_bspls, this%iflag_bspls(j))
      if (this%iflag_bspls(j) /= 0) then
        write(*,*) 'Error initializing ', j, 'D spline: ' //  &
          get_status_message(this%iflag_bspls(j))
        call exit
      end if
    end do

    n = this%n
    x = linspace(0.0e0_wp, 1.0e0_wp, n)

    idx_bspls = 0
    do j = 1, d
      do i = 1, n
        call this%bspls(j)%evaluate(  &
          x(i), idx_bspls, phi(i, j), this%iflag_bspls(j))
      end do
    end do

    do i = 1, n - 1
      dx(i) = norm2(phi(i + 1, :) - phi(i, :))
    end do
    do i = 2, n
      x(i) = x(i - 1) + dx(i - 1)
    end do

    do i = 1, n - 1
      dphidx = norm2(phi(i + 1, :) - phi(i, :)) /  &
        (x(i + 1) - x(i))
      if (.not. is_equal(abs(dphidx), 1.0e0_wp)) then
        write(*,*) "spline path went wrong"
        write(*,*) i, dphidx
!       call exit
      end if
    end do

    this%x = x
    this%phi = phi
    do i = 1, n
      this%pot(i) = this%V(this%phi(i, :))
    end do

    ! Above I changed x such that |d phi / d x| = 1
    !   -> Construct spline again with new x
    do i = 1, this%n_spls
      if (i == 1) then
        x_spline(i) = x(i)
        phi_spline(i, :) = phi(i, :)
      else if (i == this%n_spls) then
        x_spline(i) = x(this%n)
        phi_spline(i, :) = this%phi_false
      else
        j = 1 + (i - 1) * n / this%n_spls
        x_spline(i) = x(j)
        phi_spline(i, :) = phi(j, :)
      end if
    end do
    do j = 1, d
      call this%bspls(j)%initialize(  &
        x_spline,  &
        phi_spline(:, j), this%kx_bspls, this%iflag_bspls(j))
      if (this%iflag_bspls(j) /= 0) then
        write(*,*) 'Error initializing ', j, 'D spline: ' //  &
          get_status_message(this%iflag_bspls(j))
        call exit
      end if
    end do

    ! Reduce rho_max if initial guess
    ! overestimates the wall width
    this%rho_max = minval([  &
      this%rho_max,  &
      2.0e0_wp * this%rho(this%ib_max)])

  end subroutine deform_path

  subroutine check_convergence(this, iteration, conv)

    class(solver), intent(inout) :: this
    integer, intent(in) :: iteration
    logical, intent(inout) :: conv

    real(wp) :: N_max
    real(wp) :: P_max

    N_max = maxval(abs(this%Nforce))
    P_max = maxval(abs(this%Pforce))

    this%Rforce = N_max / P_max
    if (this%verbose_level >= 1) then
      write(*,*) "Force ratio =", this%Rforce
    end if

    if (iteration == 1) then
      this%Rforce_initial = this%Rforce
    end if

    if (this%Rforce < this%Rforce_best) then
      this%Rforce_best = this%Rforce
      this%S_best = this%S
      this%Sp_best = this%Sp
      this%Sk_best = this%Sk
    end if

    if (this%Rforce < this%Rforce_threshold) then
      conv = .true.
      if (this%verbose_level >= 1) then
        write(*,*) "Converged succesfully. Force ratio:", this%Rforce
      end if
      return
    else
      conv = .false.
    end if

    if (this%Rforce > 2 * this%Rforce_previous) then
      conv = .true.
      if (this%verbose_level >= 1) then
        write(*,*) "Stopped before reaching Rforce_threshold"
        write(*,*) "because Rforce started to increase."
      end if
    end if

    if ((this%Rforce > 1.5e0_wp * this%Rforce_initial) .and.  &
      (iteration > 10)) then
      conv = .true.
      if (this%verbose_level >= 1) then
        write(*,*) "Stopped before reaching Rforce_threshold"
        write(*,*) "because Rforce increase above the one of"
        write(*,*) "straight-path approximation."
      end if
    end if

    if (iteration == this%maximum_iterations) then
      conv = .true.
      if (this%verbose_level >= 1) then
        write(*,*) "Maximum iterations reached:", this%maximum_iterations
      end if
    end if

    if (is_equal(this%Rforce / this%Rforce_previous, 1.0e0_wp, eps=1.0e-2_wp)) then
      this%deform_eps = 2.0e0_wp * this%deform_eps
      if (this%verbose_level >= 1) then
        write(*,*) "Slow convergence: Increased deform_eps to", this%deform_eps
      end if
    end if

    if (this%Rforce / this%Rforce_previous > 1.0e0_wp) then
      this%deform_eps = 0.5e0_wp * this%deform_eps
      if (this%verbose_level >= 1) then
        write(*,*) "Bad deformation: Decrease deform_eps to", this%deform_eps
      end if
    end if

    this%Rforce_previous = this%Rforce

  end subroutine check_convergence

end module bouncesolver__pdm
