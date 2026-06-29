module bouncesolver__pdm

  use bouncesolver__config, only : wp
  use bouncesolver__config, only : fmt
  use bouncesolver__config, only : fmt2
  use bouncesolver__config, only : fmt3
  use bouncesolver__util, only : linspace
  use bouncesolver__util, only : pi
  use bouncesolver__util, only : riemann_integrate
  use bouncesolver__util, only : is_equal
  use bouncesolver__specialfunction, only : bessel_mod_1
  use bouncesolver__specialfunction, only : bessel_mod_onehalf
  use odeint__integrate, only : integrate_ode => integrate
  use gradmin__descent, only : fmin_descent => minimize
  use gradmin__derivatives, only : gradient
  use bspline_module, only : bspline_1d
  use bspline_module, only : get_status_message
  use cubicsplines__smoothing, only : add_viscosity
  use linalg_matrices_diagonalize, only : svd_square
  use ieee_arithmetic, only : ieee_is_nan

  implicit none

  private

  real(wp), parameter :: eps_gradient = 1.0e-8_wp

  integer, parameter :: alpha_default = 2
  real(wp), parameter :: x_min_start = 1.0e-8_wp
  real(wp), parameter :: rho_max_fac_default = 20.0e0_wp

  integer, parameter :: solver_status_ok = 0
  integer, parameter :: fail_construct_starting_path = -11
  integer, parameter :: fail_estimate_rho_max = -12
  integer, parameter :: fail_calc_msq_true = -13
  integer, parameter :: fail_bounce_on_path_thin_init = -14
  integer, parameter :: fail_bounce_on_path_thin_stitch = -15
  integer, parameter :: fail_solver_alpha = -16
  integer, parameter :: fail_deform_path_initsplines = -17
  integer, parameter :: fail_deform_path_evalsplines = -18
  integer, parameter :: fail_deform_path_evalsplines_nans = -19
  integer, parameter :: fail_deform_path_evalsplines_unitlen = -20
  integer, parameter :: fail_d2V_dx2 = -21
  integer, parameter :: fail_dV_dx = -22
  integer, parameter :: fail_find_x_barrier = -23
  integer, parameter :: fail_apply_blending = -24
  integer, parameter :: fail_xdot_wiggles = -25
! integer, parameter :: fail_bounce_on_path_thin_notfit = -26
  integer, parameter :: fail_bounce_on_path_thin_shooting = -27
  integer, parameter :: fail_check_minima = -28
  integer, parameter :: status_maxiter = 1

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
    integer :: n = 2000
    integer :: n_spls = 20
    integer :: maximum_iterations = 100
    integer :: verbose_level = 1
    real(wp) :: Rforce_threshold = 0.1e0_wp
    real(wp) :: Rforce
    real(wp) :: Rforce_previous = 1.0e10_wp
    real(wp) :: Rforce_best = 1.0e10_wp
    real(wp) :: deform_eps = 0.02e0_wp
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
    integer, allocatable :: ixs(:)
    real(wp) :: gradV_max
    real(wp), allocatable :: phi_deformed(:, :)
    type(bspline_1d), allocatable :: bspls(:)
    integer, allocatable :: iflag_bspls(:)
    integer :: kx_bspls
    real(wp), allocatable :: x_best(:)
    real(wp), allocatable :: xb_best(:)
    real(wp), allocatable :: xbdot_best(:)
    real(wp), allocatable :: rho_best(:)
    real(wp), allocatable :: phi_best(:, :)
    real(wp), allocatable :: pot_best(:)
    real(wp), allocatable :: Nforce_best(:, :)
    real(wp), allocatable :: Pforce_best(:, :)
    real(wp), allocatable :: xforce_best(:)
    integer :: iter_best
    integer :: ib_max_best
    logical :: smoothing = .false.
    real(wp) :: msq_false
    real(wp) :: msq_true
    integer :: exit_status
    integer :: init_path_mode = 1
    real(wp) :: barrier_buffer = 1.0e-3_wp
  contains
    procedure, private :: allocate_arrays
    procedure, private :: construct_starting_path
    procedure, private :: find_x_barrier
    procedure, private :: d2V_dx2
    procedure, private :: estimate_rho_max
    procedure, private :: get_index_from_x
    procedure, private :: bounce_on_path
    procedure, private :: bounce_on_path_thin
    procedure, private :: dV_dx
    procedure, private :: clip_bounce
    procedure, private :: calc_action
    procedure, private :: calc_forces
    procedure, private :: gradV
    procedure, private :: calc_phi_on_bounce
    procedure, private :: deform_path
    procedure, public :: check_convergence
    procedure, private :: calc_msq_true
    procedure, private :: reparam_path
    procedure, private :: apply_blending
    procedure, public :: print_exit_status
    procedure, private :: check_minima
  end type solver

  interface solver
    module procedure create_solver
  end interface solver

contains

  function create_solver(  &
    d, V, phi_false, phi_true,  &
    alpha, rho_max_fac, x_min,  &
    n_odeint,  &
    maxiter, deform_eps,  &
    Rforce_threshold,  &
    smoothing,  &
    init_path_mode,  &
    barrier_buffer,  &
    verbose_level) result(this)

    integer, intent(in) :: d
    procedure(V_abstract) :: V
    real(wp), intent(in) :: phi_false(d)
    real(wp), intent(in) :: phi_true(d)
    integer, intent(in), optional :: alpha
    real(wp), intent(in), optional :: rho_max_fac
    real(wp), intent(in), optional :: x_min
    integer, intent(in), optional :: n_odeint
    integer, intent(in), optional :: maxiter
    real(wp), intent(in), optional :: deform_eps
    real(wp), intent(in), optional :: Rforce_threshold
    logical, intent(in), optional :: smoothing
    integer, intent(in), optional :: init_path_mode
    integer, intent(in), optional :: verbose_level
    real(wp), intent(in), optional :: barrier_buffer
    type(solver) :: this

    logical :: converged
    integer :: iteration
    logical :: too_thin

    this%exit_status = solver_status_ok

    this%d = d
    this%V => V
    this%phi_false = phi_false
    this%phi_true = phi_true

    call this%check_minima()
    if (this%exit_status /= solver_status_ok) then
      call this%print_exit_status()
      return
    end if

    if (present(alpha)) then
      if ((alpha < 2) .or. (alpha > 3)) then
        this%exit_status = fail_solver_alpha
      end if
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

    if (present(n_odeint)) this%n = n_odeint

    if (present(maxiter)) this%maximum_iterations = maxiter

    if (present(deform_eps)) this%deform_eps = deform_eps

    if (present(Rforce_threshold)) this%Rforce_threshold = Rforce_threshold

    if (present(smoothing)) this%smoothing = smoothing

    if (present(init_path_mode)) this%init_path_mode = init_path_mode

    if (present(barrier_buffer)) this%barrier_buffer = barrier_buffer

    if (present(verbose_level)) this%verbose_level = verbose_level

    this%kx_bspls = 5

    call this%allocate_arrays()
    if (this%exit_status /= solver_status_ok) then
      call this%print_exit_status()
      return
    end if

    call this%construct_starting_path(this%init_path_mode)
    if (this%exit_status /= solver_status_ok) then
      call this%print_exit_status()
      return
    end if

    call this%estimate_rho_max()
    if (this%exit_status /= solver_status_ok) then
      call this%print_exit_status()
      return
    end if

    call this%calc_msq_true()
    if (this%exit_status /= solver_status_ok) then
      call this%print_exit_status()
      return
    end if

    converged = .false.
    iteration = 0

    do while (.not. converged)

      iteration = iteration + 1

      if (this%verbose_level > 0) then
        write(*,*) "========================================"
        write(*,*) "Iteration = ", iteration
      end if

      call this%bounce_on_path(too_thin)
      if (this%exit_status /= solver_status_ok) then
        call this%print_exit_status()
        return
      end if

      if (too_thin) then

        if (this%verbose_level > 0) then
          write(*,*) "  Wall is thin: approximate solution."
        end if

        call this%bounce_on_path_thin()
        if (this%exit_status /= solver_status_ok) then
          call this%print_exit_status()
          return
        end if

      end if

      call this%reparam_path()

      call this%calc_action()
      if (this%exit_status /= solver_status_ok) then
        call this%print_exit_status()
        return
      end if

      if (this%verbose_level > 0) then
        write(*,'(a)', advance='no') "   Action = "
        write(*,fmt3) this%S, this%Sp, this%Sk
      end if

      if (iteration == 1) this%S0 = this%S

      call this%calc_forces(iteration)

      call this%check_convergence(iteration, too_thin, converged)

      if (.not. converged) call this%deform_path()
      if (this%exit_status < solver_status_ok) then
        call this%print_exit_status()
        return
      end if

      if (this%verbose_level > 0) then
        write(*,*)
      end if

    end do

    this%S = this%S_best
    this%Sp = this%Sp_best
    this%Sk = this%Sk_best
    this%x = this%x_best
    this%xb = this%xb_best
    this%xbdot = this%xbdot_best
    this%rho = this%rho_best
    this%phi = this%phi_best
    this%pot = this%pot_best
    this%Nforce = this%Nforce_best
    this%Pforce = this%Pforce_best
    this%xforce = this%xforce_best
    this%ib_max = this%ib_max_best

    if (this%verbose_level >= 1) then
      write(*,*)
      write(*,*) "========================================"
      write(*,*) "Result:"
      write(*,'(a)', advance='no') "       Action = "
      write(*,fmt3) this%S, this%Sp, this%Sk
      write(*,'(a)', advance='no') "  Force ratio ="
      write(*,fmt) this%Rforce_best
      write(*,*) "   Iteration =", this%iter_best, "/", iteration
      write(*,*) "========================================"
      write(*,*)
    end if

  end function create_solver

  subroutine check_minima(this)

    class(solver), intent(inout) :: this

    integer :: i
    logical :: are_same

    ! Put here tests for minima

    ! Are they the same?
    are_same = .true.
    do i = 1, this%d
      if (.not. is_equal(this%phi_true(i), this%phi_false(i), eps=1.0e-8_wp)) then
        are_same = .false.
        exit
      end if
    end do
    if (are_same) this%exit_status = fail_check_minima

  end subroutine check_minima

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
    allocate(this%x_best(n))
    allocate(this%xb_best(n))
    allocate(this%xbdot_best(n))
    allocate(this%rho_best(n))
    allocate(this%phi_best(n, d))
    allocate(this%pot_best(n))
    allocate(this%Nforce_best(n, d))
    allocate(this%Pforce_best(n, d))

  end subroutine allocate_arrays

  subroutine construct_starting_path(this, mode)

    class(solver), intent(inout) :: this
    integer, intent(in) :: mode

    real(wp) :: x(this%n)
    real(wp) :: phi(this%n, this%d)
    real(wp) :: dx(this%n - 1)
    real(wp) :: pot(this%n)
    integer :: i
    integer :: n
    real(wp) :: dphidx
    real(wp) :: phi1_pivot
    real(wp) :: vmin
    real(wp) :: pmin(this%d - 1)
    real(wp) :: phi_pivot(this%d)
    real(wp) :: t1(this%d)
    real(wp) :: t2(this%d)
    real(wp) :: t(this%d)
    real(wp) :: tanvec(this%d)
    real(wp) :: proj(this%d, this%d)
    real(wp) :: evproj(this%d)
    real(wp) :: leftproj(this%d, this%d)
    real(wp) :: rightproj(this%d, this%d)
    integer :: j
    integer :: k
    integer :: svd_kont
    integer :: svd_order
    real(wp) :: orthobasis(this%d - 1, this%d)
    real(wp) :: restrictedx0(this%d - 1)
    integer :: minimize_status

    n = this%n

    x = linspace(0.0e0_wp, 1.0e0_wp, n)

    svd_order = 1
    svd_kont = 1
    restrictedx0 = 0.0e0_wp

    if (mode == 1) then ! Straight path
      do i = 1, n
        phi(i, :) = this%phi_true +  &
          x(i) * (this%phi_false - this%phi_true)
      end do
    else if (mode == 2) then ! Path along well
      do i = 1, n
        phi(i, :) = this%phi_true +  &
          x(i) * (this%phi_false - this%phi_true)
        phi1_pivot = phi(i, 1)
        call fmin_descent(  &
          V, phi(i, 2 : this%d), pmin, vmin,  &
          status=minimize_status,  &
          maxiter=10000, mode=1)
        phi(i, 2 : this%d) = pmin
      end do
    else if (mode == 3) then ! Path along orthogonal well
      do i = 1, n
        phi_pivot = this%phi_true +  &
          x(i) * (this%phi_false - this%phi_true)
        ! Compute tangent unit vector
        if (i < n) then
          t1 = this%phi_true +  &
            x(i) * (this%phi_false - this%phi_true)
          t2 = this%phi_true +  &
            x(i + 1) * (this%phi_false - this%phi_true)
        else
          t1 = this%phi_true +  &
            x(i - 1) * (this%phi_false - this%phi_true)
          t2 = this%phi_true +  &
            x(i) * (this%phi_false - this%phi_true)
        end if
        t = t2 - t1
        t = t / norm2(t)
        ! Compute projector
        proj = 0.0e0_wp
        do j = 1, this%d
          do k = 1, this%d
            if (j == k) proj(j, k) = 1.0e0_wp
            proj(j, k) = proj(j, k) - t(j) * t(k)
          end do
        end do
        call svd_square(  &
          this%d,  &
          proj,  &
          evproj,  &
          leftproj,  &
          rightproj,  &
          svd_kont,  &
          svd_order)
        ! Check eigenvalues
        do j = 1, this%d - 1
          if (.not. is_equal(evproj(j), 1.0e0_wp)) then
            write(*, *) "cadsscsd"
            call exit
          end if
        end do
        if (.not. is_equal(evproj(this%d), 0.0e0_wp)) then
          write(*, *) "adadxsda"
          call exit
        end if
        ! First d - 1 columns of leftproj (the columns belonging
        ! to the eigenvalues equal to 1) are the eigenvectors
        ! of the orthogonal basis
        do j = 1, this%d - 1
          orthobasis(j, :) = leftproj(j, :)
        end do
        ! Check that eigenvectors are orthogonal to t
        do j = 1, this%d - 1
          if (.not. is_equal(dot_product(orthobasis(j, :), t), 0.0e0_wp)) then
            write(*,*) "crtythrt"
            call exit
          end if
        end do
        call fmin_descent(  &
          restrictedV, restrictedx0,  &
          pmin, vmin,  &
          status=minimize_status,  &
          maxiter=10000, mode=1)
        phi(i, :) = phi_pivot
        do j = 1, this%d - 1
          phi(i, :) = phi(i, :) + pmin(j) * orthobasis(j, :)
        end do
      end do
    end if

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
        this%exit_status = fail_construct_starting_path
        return
      end if
    end do

    do i = 1, n
      pot(i) = this%V(phi(i, :))
    end do

    this%x = x
    this%phi = phi
    this%pot = pot

    contains

      function restrictedV(x) result(y)

        real(wp), intent(in) :: x(:)
        real(wp) :: y

        integer :: i
        real(wp) :: p(this%d)

        p = phi_pivot
        do i = 1, this%d - 1
          p = p + x(i) * orthobasis(i, :)
        end do
        y = this%V(p)

      end function restrictedV

      function V(x) result(y)

        real(wp), intent(in) :: x(:)
!       real(wp), intent(in) :: x(this%d - 1)
        real(wp) :: y

        real(wp) :: x_pivot(this%d)

        if (this%d > 2) then
          x_pivot(2 : this%d) = x
          x_pivot(1) = phi1_pivot
        else
          x_pivot(2) = x(1)
          x_pivot(1) = phi1_pivot
        end if

        y = this%V(x_pivot)

      end function V

  end subroutine construct_starting_path

  function d2V_dx2(this, i) result(d2V)

    class(solver), intent(inout) :: this
    integer, intent(in) :: i
    real(wp) :: d2V

    real(wp) :: Vm
    real(wp) :: V0
    real(wp) :: Vp
    real(wp) :: hm
    real(wp) :: hp
    real(wp) :: pot(this%n)
    real(wp) :: x(this%n)
    integer :: n
    integer :: j
    integer :: dj

    n = this%n
    dj = 1

    x = this%x
    pot = this%pot

    ! For now always use central differences.
    ! Should be fine if n large because
    ! d2V(i + 1) ~ d2V(i) ~ d2V(i - 1)

    if ((i > dj) .and. (i < n - dj)) then
      j = i
    else if (i < n - 2 * dj) then
      j = i + dj
    else if (i > 2 * dj) then
      j = i - dj
    else
      this%exit_status = fail_d2V_dx2
      return
    end if

    Vm = pot(j - dj)
    V0 = pot(j)
    Vp = pot(j + dj)
    hm = x(j) - x(j - dj)
    hp = x(j + dj) - x(j)
    d2V = 2.0e0_wp * (  &
      (Vp - V0) / hp -  &
      (V0 - Vm) / hm) / (hm + hp)

  end function d2V_dx2

  function dV_dx(this, i) result(dV)

    class(solver), intent(inout) :: this
    integer, intent(in) :: i
    real(wp) :: dV

    real(wp) :: Vm
    real(wp) :: Vp
    real(wp) :: h
    real(wp) :: pot(this%n)
    real(wp) :: x(this%n)
    integer :: n
    integer :: di

    n = this%n
    di = 1

    x = this%x
    pot = this%pot

    if ((i > di) .and. (i < n - di)) then
      Vm = pot(i - di)
      Vp = pot(i + di)
      h = x(i + di) - x(i - di)
    else if (i < n - 2 * di) then
      Vm = pot(i)
      Vp = pot(i + di)
      h = x(i + di) - x(i)
    else if (i > 2 * di) then
      Vm = pot(i - di)
      Vp = pot(i)
      h = x(i) - x(i - di)
    else
      this%exit_status = fail_dV_dx
      return
    end if

    dV = (Vp - Vm) / h

  end function dV_dx

  subroutine estimate_rho_max(this)

    class(solver), intent(inout) :: this
    real(wp) :: msq_false

    msq_false = this%d2V_dx2(this%n)

    if (msq_false <= 0.0e0_wp) then
      this%exit_status = fail_estimate_rho_max
      return
    end if

    this%msq_false = msq_false

    this%rho_max = this%rho_max_fac / sqrt(msq_false)

  end subroutine estimate_rho_max

  subroutine get_index_from_x(this, a, ia, overshot)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: a
    integer, intent(out) :: ia
    logical, optional :: overshot

    real(wp) :: x(this%n)
    integer :: n
    integer :: low
    integer :: high
    integer :: mid
    logical :: os

    n = this%n
    x = this%x

    os = .false.

    if (a <= x(1)) then

      ia = 1

    else if (a >= x(n)) then

      ia = n
      os = .true.

    else

      low = 1
      high = n

      do while (high - low > 1)
        mid = (low + high) / 2
        if (x(mid) == a) then
          ia = mid
          return
        else if (x(mid) < a) then
          low = mid
        else
          high = mid
        end if
      end do

      if (abs(x(low) - a) <= abs(x(high) -a)) then
        ia = low
      else
        ia = high
      end if
    end if

    if (present(overshot)) then
      overshot = os
    end if

  end subroutine get_index_from_x

  subroutine bounce_on_path(this, too_thin)

    class(solver), intent(inout) :: this
    logical, intent(out) :: too_thin

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

    n = this%n
    x_min = this%x_min
    x_max = this%find_x_barrier()
    if (x_max < 0.0e0_wp) then
      this%exit_status = fail_find_x_barrier
      return
    end if

    rho_min = 1.0e-4_wp
    rho_max = this%rho_max

    over_under_flag = 0
    shooting_converged = .false.

    too_thin = .false.

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

      call integrate_ode(  &
        dxdrho,  &
        x0,  &
        rho_min,  &
        rho_max,  &
        n,  &
        rho, x,  &
        method="rk4")

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
      too_thin = .true.
      return
    end if

    this%rho = rho
    this%xb = x(:, 1)
    this%xbdot = x(:, 2)

    contains

      function dxdrho(x, rho) result(xdot)

        real(wp), intent(in) :: x(:)
        real(wp), intent(in) :: rho
        real(wp), allocatable :: xdot(:)

        integer :: ix

        allocate(xdot(2))

        call this%get_index_from_x(x(1), ix, overshot)

        if (.not. overshot) then
          xdot(1) = x(2)
          xdot(2) = this%dV_dx(ix) -  &
            this%alpha * x(2) / rho
        else
          xdot(1) = 0.0e0_wp
          xdot(2) = 0.0e0_wp
        end if

      end function dxdrho

  end subroutine bounce_on_path

  subroutine reparam_path(this)

    class(solver), intent(inout) :: this

    integer :: i
    real(wp) :: phi(this%d)

    do i = 1, this%n
      call this%get_index_from_x(  &
        this%xb(i),  &
        this%ixs(i))
    end do

    call this%calc_phi_on_bounce()

    do i = 1, this%n
      phi = this%phi(i, :)
      this%pot(i) = this%V(phi)
    end do

  end subroutine reparam_path

  subroutine bounce_on_path_thin(this)

    class(solver), intent(inout) :: this

    integer :: n
    integer :: n_apprx
    integer :: n_num
    real(wp) :: x0(2)
    real(wp) :: xB
    real(wp) :: xT
    real(wp) :: del
    real(wp) :: del_min
    real(wp) :: del_max
    real(wp) :: b
    real(wp) :: nu
    integer ::  D
    real(wp) :: x_eps
    real(wp) :: rho_eps
    integer :: over_under_flag
    logical :: shooting_converged
    logical :: overshot
    logical :: bubble_not_fitting
    real(wp) :: r
    real(wp) :: r_min
    real(wp) :: r_max
    real(wp) :: rho_apprx(this%n)
    real(wp), allocatable :: rho_mask(:)
    real(wp) :: drho
    real(wp) :: xb_apprx(this%n)
    real(wp) :: xbdot_apprx(this%n)
    real(wp) :: rho_num_min
    real(wp) :: rho_num_max
    real(wp), allocatable :: rho_num(:)
    real(wp), allocatable :: x_num(:, :)
    integer :: nuinv

    integer :: i
    integer :: k
    real(wp) :: del_max_init

    n = this%n
    k = 40 * int(n / 2000)

    b = sqrt(this%msq_true)

    nu = (this%alpha - 1.0e0_wp) / 2.0e0_wp
    D = this%alpha + 1

    nuinv = 2 / (this%alpha - 1)

    del_max_init = 50.0e0_wp
    del_max = del_max_init ! -> x(0) = xT - 0 (xB - xT) ~ xT
    del_min = 3.0e0_wp   ! -> x(0) = xT - 1 (xB - xT) ~ xB

    xB = this%find_x_barrier()
    if (xB < 0.0e0_wp) then
      this%exit_status = fail_find_x_barrier
      return
    end if

    xT = this%x(1)

    r_min = 1.0e-2_wp * real(this%rho_max / n, wp)
    r_max = this%rho_max
    x_eps = xT + (xB - xT) / 50.0e0_wp ! To check
    rho_mask = linspace(r_min, r_max, n)
    drho = rho_mask(2) - rho_mask(1)
    rho_num_max = this%rho_max

    shooting_converged = .false.
    over_under_flag = 0

    ! Just to avoid that elements beyond n_apprx
    ! are not initialized
    xb_apprx = 0.0e0_wp
    xbdot_apprx = 0.0e0_wp

    do while (.not. shooting_converged)

      select case (over_under_flag)
        case (0) ! First iteration
          del = 10.0e0_wp ! 14.0e0_wp ! 12.94e0_wp
        case (1) ! Overshoot -> del smaller -> closer to xB
          del_max = del
          del = (del_max + del_min) / 2.0e0_wp
        case (-1) ! Undershoot -> del larger -> closer to xT
          del_min = del
          del = (del_max + del_min) / 2.0e0_wp
      end select

      ! Add check if del = del_max
      if (del >= del_max_init - 1.0e-4_wp) then
        if (this%verbose_level >= 1) write(*,*) del, del_max
        this%exit_status = fail_bounce_on_path_thin_shooting
        return
      end if

      ! Determine rho-value where we switch to numerical
      r = this%rho_max / 1.0e8_wp
      r_min = 1.0e-2_wp * real(this%rho_max / n, wp)
      r_max = this%rho_max
      do while (abs(r_max - r_min) > 1.0e-10_wp)
        if (x_approx(r) < x_eps) then
          r_min = r
        else
          r_max = r
        end if
        r = (r_max + r_min) / 2.0e0_wp
      end do

!     write(*,*) r, this%rho_max, del, x_eps, x_approx(r)
!     if (r >= this%rho_max * 0.8e0_wp) then
!       if (this%verbose_level >= 1) write(*,*) r, this%rho_max
!       this%exit_status = fail_bounce_on_path_thin_notfit
!       return
!     end if

      ! Check if x_eps = x_approx(r), otherwise
      ! already x(r_min) > x_eps
      if (.not. is_equal(x_eps, x_approx(r), eps=1.0e-4_wp)) then
        if (this%verbose_level >= 1) write(*,*) x_eps, x_approx(r)
        this%exit_status = fail_bounce_on_path_thin_init
        return
      end if


      ! Avoid integrating until rho_max analytically
      if (r > 0.9e0_wp * this%rho_max) then
        r = 0.6 * this%rho_max
        bubble_not_fitting = .true.
        if (this%verbose_level >= 2) write(*,*) "Probably rho_max too small."
      else
        bubble_not_fitting = .false.
      end if

      ! Take fraction r / rho_max of total points for analytic part
      n_apprx = floor((r / this%rho_max) * n)
      do i = 1, n_apprx
        rho_apprx(i) = rho_mask(i)
        xb_apprx(i) = x_approx(rho_apprx(i))
      end do

      xbdot_apprx(1) = (xb_apprx(2) - xb_apprx(1)) / drho
      do i = 2, n_apprx - 1
        xbdot_apprx(i) = (xb_apprx(i + 1) - xb_apprx(i - 1)) /  &
          (2.0e0_wp * drho)
      end do
      xbdot_apprx(n_apprx) = (xb_apprx(n_apprx) -  &
        xb_apprx(n_apprx - 1)) / drho

      ! now integrate the rest with odeint
      ! overlap region over k indices (check that n_aprrx > k)
      if (n_apprx <= k) n_apprx = k + 10
      x0(1) = xb_apprx(n_apprx - k)
      x0(2) = xbdot_apprx(n_apprx - k)
      rho_num_min = rho_mask(n_apprx  - k + 1)
      n_num = n - (n_apprx - k)
      if (n_apprx + n_num /= n + k) then
        if (this%verbose_level >= 2) write(*,*) "n_apprx + n_num /= n + k"
        this%exit_status = fail_bounce_on_path_thin_stitch
        return
      end if

      call integrate_ode(  &
        dxdrho,  &
        x0,  &
        rho_num_min,  &
        rho_num_max,  &
        n_num,  &
        rho_num, x_num,  &
        method="rk4")

      if (overshot) then
        over_under_flag = 1
      else
        over_under_flag = -1
      end if
      if (bubble_not_fitting) then
        over_under_flag = 1
      end if

      if (is_equal(  &
            del_min / del_max,  &
            1.0e0_wp,  &
            eps=1.0e-8_wp)) then
        shooting_converged = .true.
      else
        shooting_converged = .false.
      end if

      if (this%verbose_level >= 2) then
        write(*,*) del_min, del, del_max, over_under_flag
      end if

    end do

    call this%apply_blending(  &
      n_apprx, k,  &
      rho_apprx(1:n_apprx), rho_num,  &
      xb_apprx(1:n_apprx), x_num(:, 1),  &
      x_num(:, 2))

  contains

    function x_approx(rho) result(y)

      real(wp), intent(in) :: rho
      real(wp) :: y

      real(wp) :: x0
      real(wp) :: br
      real(wp) :: c

      x0 = xT + exp(-del) * (xB - xT)
      br = rho * b

      if (D == 3) then
        c = sqrt(pi / (2.0e0_wp * br))
      else if (D == 4) then
        c = 2.0e0_wp / br
      end if

      y = xT - (xT - x0) * c * Inu(br)

    end function x_approx

    function Inu(br) result(y)

      real(wp), intent(in) :: br
      real(wp) :: y

      if (nuinv == 2) then
        y = bessel_mod_onehalf(br) ! = sqrt(2.0e0_wp / (pi * br)) * sinh(br)
      else if (nuinv == 1) then
        y = bessel_mod_1(br)
      end if

    end function Inu

      function dxdrho(x, rho) result(xdot)

        real(wp), intent(in) :: x(:)
        real(wp), intent(in) :: rho
        real(wp), allocatable :: xdot(:)

        integer :: ix

        allocate(xdot(2))

        call this%get_index_from_x(x(1), ix, overshot)

        if (.not. overshot) then
          xdot(1) = x(2)
          xdot(2) = this%dV_dx(ix) -  &
            this%alpha * x(2) / rho
        else
          xdot(1) = 0.0e0_wp
          xdot(2) = 0.0e0_wp
        end if

      end function dxdrho

  end subroutine bounce_on_path_thin

  subroutine apply_blending(  &
    this, nstich, k,  &
    r1, r2,  &
    x1, x2,  &
    x2dot)

    class(solver), intent(inout) :: this
    integer, intent(in) :: nstich
    integer, intent(in) :: k
    real(wp), intent(in) :: r1(:)
    real(wp), intent(in) :: r2(:)
    real(wp), intent(in) :: x1(:)
    real(wp), intent(in) :: x2(:)
    real(wp), intent(in) :: x2dot(:)

    integer :: n
    integer :: n1
    integer :: n2
    integer :: i
    integer :: i1
    integer :: i2
    real(wp) :: kr
    real(wp) :: t
    real(wp) :: w
    real(wp) :: x(this%n)
    real(wp) :: xdot(this%n)
    real(wp) :: r(this%n)
    integer :: j1
    integer :: j2
    real(wp) :: dr

    n = this%n
    n1 = size(r1)
    n2 = size(r2)

    this%rho = 0
    this%xb = 0
    this%xbdot = 0

    if (n1 + n2 - k /= n) then
      this%exit_status = fail_apply_blending
      return
    end if

    if ((size(x1) /= n1) .or. (size(x2) /= n2)) then
      this%exit_status = fail_apply_blending
      return
    end if

    ! First part: analytic solution
    j1 = 1
    j2 = nstich - k
    x(j1:j2) = x1(j1:j2)
    r(j1:j2) = r1(j1:j2)

    ! Intermediate part: blend solution
    i1 = nstich - k
    i2 = nstich - 1
    kr = real(k, wp)
    do i = i1, i2
      t = real(i - i1, wp) / kr
      w = 1.0e0_wp - (6.0e0_wp * t ** 5 - 15.0e0_wp * t ** 4 + 10.0e0_wp * t ** 3)
      j1 = i + 1
      j2 = i - i1 + 1
      x(j1) = w * x1(j1) + (1.0e0_wp - w) * x2(j2)
      r(j1) = r1(j1)
    end do

    ! Last part: numeric solution
    j1 = k + 1
    j2 = nstich + 1
    x(j2:n) = x2(j1:n2)
    r(j2:n) = r2(j1:n2)

    ! xdot: 1 -- nstich: compute here
    ! 1 -- nstich: compute here
    dr = r(2) - r(1)
    do i = 1, nstich
      xdot(i) = (x(i + 1) - x(i)) / dr
    end do
    ! nstich + 1 -- n: take from odeint
    j1 = nstich + 1
    j2 = k + 1
    xdot(j1:n) = x2dot(j2:n2)

    this%rho = r
    this%xb = x
    this%xbdot = xdot

  end subroutine apply_blending

  subroutine calc_msq_true(this)

    class(solver), intent(inout) :: this

    this%msq_true = this%d2V_dx2(1)

    if (this%msq_true <= 0.0e0_wp) then
      this%exit_status = fail_calc_msq_true
    end if

  end subroutine calc_msq_true

  subroutine calc_phi_on_bounce(this)

    class(solver), intent(inout) :: this

    real(wp) :: x(this%n)
    real(wp) :: xb(this%n)
    real(wp) :: phi(this%n, this%d)
    real(wp) :: phib(this%n, this%d)
    integer :: ixs(this%n)
    integer :: n
    integer :: i
    integer :: j
    real(wp) :: x1
    real(wp) :: x2
    real(wp) :: dx
    real(wp) :: phi1(this%d)
    real(wp) :: phi2(this%d)
    real(wp) :: dphi(this%d)

    n = this%n
    x = this%x
    xb = this%xb
    phi = this%phi
    ixs = this%ixs

    do i = 1, n
      j = ixs(i)
      if (j == n) then
        phib(i, :) = phi(j, :)
      else if (j == 1) then
        if (xb(i) > x(j)) then
          x1 = x(j)
          x2 = x(j + 1)
          dx = x2 - x1
          dphi = phi(j + 1, :) - phi(j, :)
          phib(i, :) = phi(j, :) +  &
            ((xb(i) - x1) / dx) * dphi
        else
          phib(i, :) = phi(j, :)
        end if
      else
        if (xb(i) > x(j)) then
          x1 = x(j)
          x2 = x(j + 1)
          dx = x2 - x1
          dphi = phi(j + 1, :) - phi(j, :)
          phib(i, :) = phi(j, :) +  &
            ((xb(i) - x1) / dx) * dphi
        else if (xb(i) < x(j)) then
          x1 = x(j - 1)
          x2 = x(j)
          dx = x2 - x1
          dphi = phi(j, :) - phi(j - 1, :)
          phib(i, :) = phi(j - 1, :) +  &
            ((xb(i) - x1) / dx) * dphi
        else
          phib(i, :) = phi(i, :)
        end if
      end if
    end do

    this%x = xb
    this%phi = phib

  end subroutine calc_phi_on_bounce

  function find_x_barrier(this) result(x_barrier)

    class(solver), intent(inout) :: this
    real(wp) :: x_barrier

    real(wp) :: V1
    real(wp) :: V2
    real(wp) :: V3
    real(wp) :: V4
    integer :: i
    real(wp) :: buffer
    real(wp) :: phi(this%d)
    real(wp) :: phi_true(this%d)
    real(wp) :: phi_false(this%d)

    phi_false = this%phi(this%n, :)
    phi_true = this%phi(1, :)
    buffer = this%barrier_buffer * (this%V(phi_false) - this%V(phi_true))

    V1 = this%V(phi_true)
    ! Not start close to phi_true here
    do i = this%n / 10, this%n
      phi = this%phi(i, :)
      V2 = this%V(phi)
      if (V2 > V1) then
        V1 = V2
      else
        if (V2 < V1 - buffer) then
          x_barrier = this%x(i - 1)
          exit
        end if
      end if
    end do

    if (i == this%n + 1) then
      ! If no barrier found than retourn negative x_barrier
      ! Catch this error in bounce_on_path().
      x_barrier = -1.0e0_wp
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
    integer :: buffer
    logical :: was_clipped
    integer :: xbdot_zeros
    logical :: was_positive

    n = this%n
    xb = this%x
    xbdot = this%xbdot
    x_max = maxval(this%x)
    clip_range = 0
    buffer = 4
    was_clipped = .false.

    do i = n / 10, n
      if (is_equal(xb(i) / x_max, 1.0e0_wp, eps=1.0e-3_wp)) then
        clip_range = clip_range + 1
      else
        clip_range = 0
      end if
      if (clip_range == buffer) then
        was_clipped = .true.
        if (this%verbose_level > 0) then
          write(*,"(a,I5,a)", advance="no")  &
            "   Clipped solution" // &
            "at rho(",  i - 1, ") = "
          write(*,fmt) this%rho(i - 1)
        end if
        exit
      end if
    end do

    iclip = minval([i - 2 * buffer, n])
    this%ib_max = iclip

    if ((iclip / real(n, wp)) < 0.4e0_wp) then
      if (this%verbose_level > 0) then
        write(*,*) "Clipped early. Reduce rho_max by a factor of two."
      end if
      this%rho_max = this%rho_max / 2
    end if

    ! Check bounce after clipping
    ! by counting zeros of xbdot
    !   -> should be zero
    was_positive = .true.
    xbdot_zeros = 0
    do i = n / 10, this%ib_max
      if ((xbdot(i) < 0.0e0_wp) .and. (was_positive)) then
        xbdot_zeros = xbdot_zeros + 1
        was_positive = .false.
      end if
      if ((xbdot(i) > 0.0e0_wp) .and. (.not. was_positive)) then
        xbdot_zeros = xbdot_zeros + 1
        was_positive = .true.
      end if
    end do
    if (xbdot_zeros > 0) then
      if (this%verbose_level >= 1) write(*,*) xbdot_zeros
      this%exit_status = fail_xdot_wiggles
      return
    end if

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
    h = maxval([V * eps_gradient, 1.0e-12_wp])

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
    real(wp) :: dx
    real(wp) :: p(this%d)

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
        do j = 1, this%d
          if (i == 1) then
            d2phi_dx2(i, j) = dT_step(  &
              phi(1, j), phi(2, j), phi(3, j),  &
              x(1), x(2), x(3))
          else if (i == n) then
            d2phi_dx2(i, j) = dT_step(  &
              phi(n - 2, j), phi(n - 1, j), phi(n, j),  &
              x(n - 2), x(n - 1), x(n))
          else
            d2phi_dx2(i, j) = dT_step(  &
              phi(i - 1, j), phi(i, j), phi(i + 1, j),  &
              x(i - 1), x(i), x(i + 1))
          end if
        end do
      end do
    end if

    do i = 1, n
      p = phi(i, :) ! Temporary copy to avoid copy runtime warning
      gradV(i, :) = this%gradV(p)
    end do

    do i = 1, n
      gradV_len(i) = norm2(gradV(i, :))
    end do
    this%gradV_max = maxval(gradV_len)

    do i = 1, n
      if (i == 1) then
        dx = x(2) - x(1)
        if (dx == 0.0e0_wp) dx = dx + 1.0e-30_wp
        pathdir(i, :) = (phi(2, :) - phi(1, :)) / dx
      else if (i == n) then
        dx = x(n) - x(n - 1)
        if (dx == 0.0e0_wp) dx = dx + 1.0e-30_wp
        pathdir(i, :) = (phi(n, :) - phi(n - 1, :)) / dx
      else
        dx = x(i + 1) - x(i - 1)
        if (dx == 0.0e0_wp) dx = dx + 1.0e-30_wp
        pathdir(i, :) = (phi(i + 1, :) - phi(i - 1, :)) / dx
      end if
      if ((.not. is_equal(norm2(pathdir(i, :)), 1.0e0_wp, eps=1.0e-3_wp)) .and.  &
          this%verbose_level >= 3) then
        write(*,*) "Problem with pathdir in calc_forces."
        write(*,*) i, pathdir(i, :), norm2(pathdir(i, :))
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

        real(wp) :: dxl
        real(wp) :: dxu
        real(wp) :: dxlu

        dxu = xu - x0
        dxl = x0 - xl

        ! Avoid 1 / 0 exceptions
        if (dxu == 0.0e0_wp) dxu = dxu + 1.0e-10_wp
        if (dxl == 0.0e0_wp) dxl = dxl + 1.0e-10_wp

        dxlu = dxu + dxl

        y = 2.0e0_wp * (  &
          (Tu - T0) / dxu -  &
          (T0 - Tl) / dxl) / dxlu

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
    real(wp) :: p(this%d)

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
        ! Extend path slightly beyond release point towards
        ! true minimum to avoid boundary artifacts in the
        ! spline interpolation
        phi_spline(i, :) = phi_deformed(1, :) -  &
          1.0e2_wp * (phi_deformed(2, :) - phi_deformed(1, :))
      else if (i == this%n_spls) then
        phi_spline(i, :) = this%phi_false
      else
        j = 1 + (i - 1) * n / this%n_spls
        phi_spline(i, :) = phi_deformed(j, :)
      end if
    end do

!   do j = 1, d
!     phi_spline(:, j) = add_viscosity(  &
!       x_spline,  &
!       phi_spline(:, j),  &
!       this%n_spls,  &
!       eps=1.0e0_wp,  &
!       iterations=10)
!   end do

    ! Construct B-spline approximation to get path with
    ! this%n points and for more precise calculation
    ! of d2phi_dx2
    do j = 1, d
      call this%bspls(j)%initialize(  &
        x_spline,  &
        phi_spline(:, j), this%kx_bspls, this%iflag_bspls(j))
      if (this%iflag_bspls(j) /= 0) then
        if (this%verbose_level >= 2) then
          write(*,*) 'Error initializing ', j, 'D spline: ' //  &
            get_status_message(this%iflag_bspls(j))
        end if
        this%exit_status = fail_deform_path_initsplines
        return
      end if
    end do

    n = this%n
    x = linspace(0.0e0_wp, 1.0e0_wp, n)

    idx_bspls = 0
    do j = 1, d
      do i = 1, n
        call this%bspls(j)%evaluate(  &
          x(i), idx_bspls, phi(i, j), this%iflag_bspls(j))
        if (ieee_is_nan(phi(i, j))) then
          if (this%verbose_level >= 2) then
            write(*,*) i, j, phi(i, j), this%iflag_bspls(j)
          end if
          this%exit_status = fail_deform_path_evalsplines_nans
          return
        end if
      end do
      if (this%iflag_bspls(j) /= 0) then
        if (this%verbose_level >= 2) then
          write(*,*) 'Error evaluating ', j, 'D spline: ' //  &
            get_status_message(this%iflag_bspls(j))
        end if
        this%exit_status = fail_deform_path_evalsplines
        return
      end if
    end do

    if (this%smoothing) then
      do j = 1, d
        phi(:, j) = add_viscosity(  &
          x,  &
          phi(:, j),  &
          n,  &
          eps=2.0e0_wp,  &
          iterations=200)
      end do
    end if

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
!       if (i == 1) then
!         write(*,*) "spline path went wrong"
!         write(*,*) i, dphidx
!         write(*,*) x(i + 1), x(i)
!         write(*,*) phi(i + 1, :)
!         write(*,*) phi(i, :)
!         write(*,*) x(i + 1) - x(i)
!         write(*,*) norm2(phi(i + 1, :) - phi(i, :))
!       end if
        this%exit_status = fail_deform_path_evalsplines_unitlen
        return
      end if
    end do

    this%x = x
    this%phi = phi
    do i = 1, n
      p = this%phi(i, :)
      this%pot(i) = this%V(p)
    end do

    ! Reduce rho_max if initial guess
    ! overestimates the wall width
    this%rho_max = minval([  &
      this%rho_max,  &
      2.0e0_wp * this%rho(this%ib_max)])

  end subroutine deform_path

  subroutine check_convergence(this, iteration, too_thin, conv)

    class(solver), intent(inout) :: this
    integer, intent(in) :: iteration
    logical, intent(in) :: too_thin
    logical, intent(inout) :: conv

    real(wp) :: N_max
    real(wp) :: P_max

    N_max = maxval(abs(this%Nforce))
    P_max = maxval(abs(this%Pforce))

    this%Rforce = N_max / P_max
    if (this%verbose_level > 0) then
      write(*,'(a)', advance='no') "   Force ratio ="
      write(*,fmt) this%Rforce
    end if

    if ((this%Rforce < this%Rforce_best) .or. (iteration == 1)) then
      this%Rforce_best = this%Rforce
      this%S_best = this%S
      this%Sp_best = this%Sp
      this%Sk_best = this%Sk
      this%x_best = this%x
      this%xb_best = this%xb
      this%xbdot_best = this%xbdot
      this%rho_best = this%rho
      this%phi_best = this%phi
      this%pot_best = this%pot
      this%Nforce_best = this%Nforce
      this%Pforce_best = this%Pforce
      this%xforce_best = this%xforce
      this%ib_max_best = this%ib_max
      this%iter_best = iteration
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

!   if (this%Rforce > 2 * this%Rforce_previous) then
!     conv = .true.
!     write(*,*) "Stopped before reaching Rforce_threshold"
!     write(*,*) "because Rforce started to increase."
!   end if

    if (iteration == this%maximum_iterations) then
      conv = .true.
      this%exit_status = status_maxiter
      if (this%verbose_level >= 1) then
        write(*,*) "Maximum iterations reached:", this%maximum_iterations
      end if
    end if

!   if (is_equal(this%Rforce / this%Rforce_previous, 1.0e0_wp, eps=1.0e-2_wp)) then
!     this%deform_eps = 1.2e0_wp * this%deform_eps
!     if (this%verbose_level >= 1) then
!       write(*,*) "Slow convergence: Increased deform_eps to", this%deform_eps
!     end if
!   end if

    if (.not. too_thin) then
      if (this%Rforce / this%Rforce_previous > 1.2e0_wp) then
        this%deform_eps = 0.8e0_wp * this%deform_eps
        if (this%verbose_level >= 1) then
          write(*,*) "Bad deformation: Decrease deform_eps to", this%deform_eps
        end if
      end if
    end if

    this%Rforce_previous = this%Rforce

  end subroutine check_convergence

  subroutine print_exit_status(this)

    class(solver), intent(inout) :: this

    select case (this%exit_status)
      case (solver_status_ok)
        write(*,*) "No errors. Converged corectly."
      case (fail_construct_starting_path)
        write(*,*) "Constructing initial straight path"
        write(*,*) "guess went wrong."
      case (fail_estimate_rho_max)
        write(*,*) "Could not estimate max rho value"
        write(*,*) "because curvature in false minimum"
        write(*,*) "was computed to be negative."
        write(*,*) "Check argument phi_false."
      case (fail_calc_msq_true)
        write(*,*) "Could not calculate msq_true"
        write(*,*) "because curvature in true minimum"
        write(*,*) "was computed to be negative."
        write(*,*) "Check argument phi_true."
      case (fail_bounce_on_path_thin_init)
        write(*,*) "Initializing analytic approximation for"
        write(*,*) "small rho in bounce_on_path_thin went wrong."
      case (fail_bounce_on_path_thin_stitch)
        write(*,*) "A problem appeared in determining the"
        write(*,*) "lengths of the analytic and numeric parts"
        write(*,*) "of the solution in bounce_on_path_thin went wrong."
      case (fail_solver_alpha)
        write(*,*) "Wrong alpha value. Possible values are 2 or 3."
      case (fail_deform_path_initsplines)
        write(*,*) "An error occured when constructing the B-spline"
        write(*,*) "interpolations of the path before deformation."
      case (fail_deform_path_evalsplines_nans)
        write(*,*) "Some of the B-spline interpolations return nans."
      case (fail_deform_path_evalsplines)
        write(*,*) "During the evaluation of the B-spline interpolations"
        write(*,*) "an error occured."
      case (fail_deform_path_evalsplines_unitlen)
        write(*,*) "The B-spline approximation of the deformed path"
        write(*,*) "does not satisfy |d phi / dx| = 1."
      case (fail_d2V_dx2)
        write(*,*) "d2V_dx2 could not be evaluated for some index"
        write(*,*) "along the path."
      case (fail_dV_dx)
        write(*,*) "dV_dx could not be evaluated for some index"
        write(*,*) "along the path."
      case (fail_find_x_barrier)
        write(*,*) "Could not determine peak of potential barrier"
        write(*,*) "along the path. This could mean that either"
        write(*,*) "phi_true or phi_min are not the exact locations"
        write(*,*) "of the minima in the potential, or that there"
        write(*,*) "either of those is a saddle point."
        write(*,*) "If you are sure that there should be a barrier,"
        write(*,*) "you can try to reduce barrier_buffer (default is 1e-3)."
      case (fail_apply_blending)
        write(*,*) "Problem with sizes of arrays in apply_blending."
      case (fail_xdot_wiggles)
        write(*,*) "xdot is not monothonically increasing, but wiggles"
        write(*,*) "around. Usually this happens if the clipping of the"
        write(*,*) "bounce was unsuccessful. Try to decrease rho_max_fac."
!     case (fail_bounce_on_path_thin_notfit)
!       write(*,*) "The rho-value of the stich point would be close to or"
!       write(*,*) "larger than rho_max. Cannot construct the analytic"
!       write(*,*) "part of bounce because does not fit in rho range."
!       write(*,*) "You can try to increase rho_max_fac."
      case (fail_bounce_on_path_thin_shooting)
        write(*,*) "Always undershot in _thin. This shouldn't happen."
      case (fail_check_minima)
        write(*,*) "The true and false minima given to the solver are"
        write(*,*) "identical. Cannot find bounce."
      case (status_maxiter)
        write(*,*) "Maximum number of iterations have been performed"
        write(*,*) "before the convergence condition was satisfied."
        write(*,*) "The path might not yet be optimal."
      case default
        write(*,*) "An unknown error occured."
    end select

  end subroutine print_exit_status

end module bouncesolver__pdm
