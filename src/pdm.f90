module bouncesolver__pdm

  use bouncesolver__config, only : wp
  use bouncesolver__util, only : linspace
  use bouncesolver__util, only : pi
  use bouncesolver__util, only : riemann_integrate
  use bouncesolver__util, only : is_equal
  use bouncesolver__specialfunction, only : bessel_mod_1
  use odeint__rk4, only : integrate
  use gradmin__derivatives, only : gradient
  use bspline_module, only : bspline_1d
  use bspline_module, only : get_status_message
  use cubicsplines__smoothing, only : add_viscosity
  use ieee_arithmetic, only : ieee_is_nan

  implicit none

  private

  real(wp), parameter :: eps_gradient = 1.0e-8_wp

  integer, parameter :: alpha_default = 2
  real(wp), parameter :: x_min_start = 1.0e-18_wp
  real(wp), parameter :: rho_max_fac_default = 20.0e0_wp

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
    logical :: spline_fail = .false.
    real(wp), allocatable :: x_best(:)
    real(wp), allocatable :: xb_best(:)
    real(wp), allocatable :: xbdot_best(:)
    real(wp), allocatable :: rho_best(:)
    real(wp), allocatable :: phi_best(:, :)
    real(wp), allocatable :: pot_best(:)
    real(wp), allocatable :: Nforce_best(:, :)
    real(wp), allocatable :: Pforce_best(:, :)
    integer :: ib_max_best
    logical :: smoothing = .false.
    real(wp) :: msq_false
    real(wp) :: msq_true
  contains
    procedure, private :: allocate_arrays
    procedure, private :: construct_starting_path
    procedure, private :: find_x_barrier
    procedure, private :: d2V_dx2
    procedure, private :: estimate_rho_max
    procedure, private :: get_index_from_x
    procedure, private :: bounce_on_path
    procedure, private :: bounce_on_path_thin
    procedure, private :: bounce_on_path_thin2
    procedure, private :: dV_dx
    procedure, private :: clip_bounce
    procedure, private :: calc_action
    procedure, private :: calc_forces
    procedure, private :: gradV
    procedure, private :: calc_phi_on_bounce
    procedure, private :: deform_path
    procedure, public :: check_convergence
    procedure, private :: calc_msq_true
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
    integer, intent(in), optional :: verbose_level
    type(solver) :: this

    logical :: converged
    integer :: iteration
    logical :: too_thin

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

    if (present(n_odeint)) this%n = n_odeint

    if (present(maxiter)) this%maximum_iterations = maxiter

    if (present(deform_eps)) this%deform_eps = deform_eps

    if (present(Rforce_threshold)) this%Rforce_threshold = Rforce_threshold

    if (present(smoothing)) this%smoothing = smoothing

    if (present(verbose_level)) this%verbose_level = verbose_level

    this%kx_bspls = 5

    call this%allocate_arrays()
    call this%construct_starting_path()
    call this%estimate_rho_max()
    call this%calc_msq_true
    converged = .false.
    iteration = 0
    do while (.not. converged)
      iteration = iteration + 1
      call this%bounce_on_path(too_thin)
!     if (too_thin) then
!       ! Only for iteration > 1 because for first
!       ! iteration not dramatic because for the
!       ! straight initial path the deformation only
!       ! depends on the gradient of the potential
!       ! NO: deform_path sets rho_max to zero then...
!       call this%bounce_on_path_thin()
!     end if
!     if (too_thin) then
!       write(*,*) "Wall too thin: approximate try."
!       call this%bounce_on_path_thin()
!     end if
      if (too_thin) then
        write(*,*) "Wall too thin: approximate try."
        call this%bounce_on_path_thin2()
      end if
!     call this%bounce_on_path_thin()
      call this%calc_action()
      write(*,*) "Action = ", this%S, this%Sp, this%Sk
      if (iteration == 1) this%S0 = this%S
      call this%calc_forces(iteration)
      call this%check_convergence(iteration, converged)
      if (.not. converged) call this%deform_path()
      if (this%spline_fail) exit
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
    this%ib_max = this%ib_max_best

    if (this%verbose_level >= 1) then
      write(*,*)
      write(*,*) "Result:"
      write(*,*) "  Action = ", this%S, this%Sp, this%Sk
      write(*,*) "  Force ratio =", this%Rforce_best
      write(*,*) "  Iterations =", iteration
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
    allocate(this%x_best(n))
    allocate(this%xb_best(n))
    allocate(this%xbdot_best(n))
    allocate(this%rho_best(n))
    allocate(this%phi_best(n, d))
    allocate(this%pot_best(n))
    allocate(this%Nforce_best(n, d))
    allocate(this%Pforce_best(n, d))

  end subroutine allocate_arrays

  subroutine construct_starting_path(this)

    class(solver), intent(inout) :: this

    real(wp) :: x(this%n)
    real(wp) :: phi(this%n, this%d)
    real(wp) :: dx(this%n - 1)
    real(wp) :: pot(this%n)
    integer :: i
    integer :: n
    real(wp) :: dphidx

    n = this%n

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

    do i = 1, n
      pot(i) = this%V(phi(i, :))
    end do

    this%x = x
    this%phi = phi
    this%pot = pot

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
      write(*,*) "Problem in d2V_dx2"
      call exit
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
      write(*,*) "Problem in dV_dx"
      call exit
    end if

    dV = (Vp - Vm) / h

  end function dV_dx

  subroutine estimate_rho_max(this)

    class(solver), intent(inout) :: this
    real(wp) :: msq_false

    msq_false = this%d2V_dx2(this%n)
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
      too_thin = .true.
      return
    end if

    this%rho = rho
    this%xb = x(:, 1)
    this%xbdot = x(:, 2)
    do i = 1, n
      call this%get_index_from_x(  &
        this%xb(i),  &
        this%ixs(i))
    end do

    call this%calc_phi_on_bounce()

    do i = 1, n
      this%pot(i) = this%V(this%phi(i, :))
    end do

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

! subroutine bounce_on_path_thin(this)

!   class(solver), intent(inout) :: this

!   real(wp) :: x_match
!   real(wp) :: msq_true
!   real(wp) :: d2V
!   logical :: shooting_converged
!   real(wp) :: x_min
!   real(wp) :: x_max
!   real(wp) :: rho_min
!   real(wp) :: rho_max
!   real(wp) :: x1
!   real(wp) :: x0(2)
!   integer :: over_under_flag
!   logical :: overshot
!   integer :: n
!   real(wp), allocatable :: rho(:)
!   real(wp), allocatable :: x(:, :)
!   integer :: i

!   n = this%n

!   msq_true = this%d2V_dx2(1)
!   do i = 1, n
!     d2V = this%d2V_dx2(i)
!     ! At the bubble wall d2V = 0, abs(d2V - msq_true) ~ msq_true
!     ! I want to use thin-walled approximation not until bubble wall,
!     ! but only a subrange below bubble wall, so added factor 4.
!     if (abs(d2V - msq_true) > msq_true) then
!       x_match = this%x(i)
!       exit
!     end if
!   end do
!   if ((i - 1) == n) then
!     write(*,*) "x_match problem"
!     call exit
!   end if

!   x_min = this%x_min
!   x_max = this%find_x_barrier()

!   rho_min = 1.0e-4_wp
!   rho_max = this%rho_max

!   over_under_flag = 0
!   shooting_converged = .false.

!   do while (.not. shooting_converged)

!     select case (over_under_flag)
!       case (0) ! First iteration
!         x1 = (x_max + x_min) / 100.0e0_wp
!       case (1) ! Overshoot
!         x_min = x1
!         x1 = (x_max + x_min) / 2.0e0_wp
!       case (-1) ! Undershoot
!         x_max = x1
!         x1 = (x_max + x_min) / 2.0e0_wp
!     end select

!     x0(1) = x1
!     x0(2) = 0.0e0_wp

!     call integrate(  &
!       dxdrho,  &
!       x0,  &
!       rho_min,  &
!       rho_max,  &
!       n,  &
!       rho, x)

!     if (overshot) then
!       over_under_flag = 1
!     else
!       over_under_flag = -1
!     end if

!     if (is_equal(  &
!           x_min / x_max,  &
!           1.0e0_wp,  &
!           eps=1.0e-10_wp)) then
!       shooting_converged = .true.
!     else
!       shooting_converged = .false.
!     end if

!     if (this%verbose_level >= 2) then
!       write(*,*) 'xxx', x_min, x_max, over_under_flag
!     end if

!   end do

!   this%rho = rho
!   this%xb = x(:, 1)
!   this%xbdot = x(:, 2)
!   do i = 1, n
!     call this%get_index_from_x(  &
!       this%xb(i),  &
!       this%ixs(i))
!   end do

!   call this%calc_phi_on_bounce()

!   do i = 1, n
!     this%pot(i) = this%V(this%phi(i, :))
!   end do

!   contains

!     function dxdrho(x, rho) result(xdot)

!       real(wp), intent(in) :: x(:)
!       real(wp), intent(in) :: rho
!       real(wp), allocatable :: xdot(:)

!       integer :: ix

!       allocate(xdot(2))

!       call this%get_index_from_x(x(1), ix, overshot)

!       if (.not. overshot) then
!         if (x(1) > x_match) then
!           xdot(1) = x(2)
!           xdot(2) = this%dV_dx(ix) -  &
!             this%alpha * x(2) / rho
!         else
!           xdot(1) = x(2)
!           xdot(2) = msq_true * x(1) -  &
!             this%alpha * x(2) / rho
!         end if
!       else
!         xdot(1) = 0.0e0_wp
!         xdot(2) = 0.0e0_wp
!       end if

!     end function dxdrho

! end subroutine bounce_on_path_thin

  subroutine bounce_on_path_thin2(this)

    class(solver), intent(inout) :: this

    integer :: n
!   real(wp) :: msq_true
    real(wp) :: x0
    real(wp) :: V0
    real(wp) :: dV
    real(wp) :: d2V
    real(wp) :: nu
    real(wp) :: gam
!   real(wp) :: V_true
    integer :: i
    integer :: j
    integer :: i_switch
    real(wp) :: x_switch
    real(wp) :: rho_switch
    real(wp) :: beta
    real(wp) :: br
    real(wp) :: rho(this%n)
    real(wp) :: drho
    real(wp) :: x1
    real(wp) :: x_start(2)
    real(wp) :: x_min
    real(wp) :: x_max
    real(wp) :: rho_min
    real(wp) :: rho_max
    integer :: over_under_flag
    logical :: shooting_converged
    logical :: too_thin
    real(wp) :: x(this%n, 2)
    real(wp), allocatable :: rho_ode(:)
    real(wp), allocatable :: x_ode(:, :)
    logical :: overshot

    n = this%n
    nu = 0.5e0_wp * (this%alpha - 1.0e0_wp)
    gam = gamma(nu + 1.0e0_wp)

    x_min = this%x_min
    x_max = this%find_x_barrier()

    rho_min = 1.0e-4_wp
    rho_max = this%rho_max
    rho = linspace(1.0e-4_wp, rho_max, n)
    drho = rho(2) - rho(1)

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

!     x_start(1) = x1
!     x_start(2) = 0.0e0_wp

      ! Construct harmonic approximation of V
!     msq_true = this%msq_true
!     V_true = this%V(this%phi_true)
      call this%get_index_from_x(x1, j)
      x0 = this%x(j)
      V0 = this%pot(j)
      dV = this%dV_dx(j)
      d2V = this%d2V_dx2(j)

      ! Check until which x-value is good approximation
      do i = j, n
        if (this%d2V_dx2(i) < 0.0e0_wp) then
!       if (.not. is_equal(Vapprx(this%x(i)) / this%pot(i), 1.0e0_wp, eps=2.0e2_wp)) then
          x_switch = this%x(i)
          exit
        end if
      end do
!     write(*,*) dV, d2V, x_switch

      ! Use analytic solution in valid x-range
!     write(*,*) 'd2V', d2V, x0
      beta = sqrt(d2V)
      do i = 1, n
        br = beta * rho(i)
        x(i, 1) = (dV / d2V) * (gam * (br / 2.0e0_wp) ** (-nu) * Inu(br) - 1.0e0_wp)
        if (x(i, 1) > x_switch) then
          i_switch = i - 1
          exit
        end if
      end do
!     write(*,*) i_switch, rho(i_switch), x(i_switch, 1)

      ! Get xdot for analytic part
      x(1, 2) = (x(2, 1) - x(1, 1)) / drho
      do i = 2, i_switch - 1
        x(i, 2) = (x(i + 1, 1) - x(i - 1, 1)) / (2.0e0_wp * drho)
      end do
      x(i_switch, 2) = (x(i_switch, 1) - x(i_switch - 1, 1)) / drho

      ! Integrate the remaining rho region with odeint
      x_start = x(i_switch, :)
      rho_switch = rho(i_switch)
      call integrate(  &
        dxdrho,  &
        x_start,  &
        rho_switch,  &
        rho_max,  &
        n - i_switch + 1,  &
        rho_ode, x_ode)

      ! Stitching together thin + odeint parts
      rho(i_switch:n) = rho_ode
      x(i_switch:n, :) = x_ode

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

    this%rho = rho
    this%xb = x(:, 1)
    this%xbdot = x(:, 2)
    do i = 1, n
      call this%get_index_from_x(  &
        this%xb(i),  &
        this%ixs(i))
    end do

    call this%calc_phi_on_bounce()

    do i = 1, n
      this%pot(i) = this%V(this%phi(i, :))
    end do

  contains

    function Vapprx(x) result(y)

      real(wp), intent(in) :: x
      real(wp) :: y

!     y = msq_true * x ** 2 + V_true
      y = V0 + dV * (x - x0) + 0.5e0_wp * d2V * (x - x0) ** 2

    end function

    function Inu(br) result(y)

      real(wp), intent(in) :: br
      real(wp) :: y

      if (this%alpha == 2) then
        ! Source: 1. https://ameli.github.io/special_functions/api/besseli.html
        !         2. https://archive.lib.msu.edu/crcmath/math/math/m/m312.htm
        y = sqrt(2.0e0_wp / (pi * br)) * sinh(br)
      else if (this%alpha == 3) then
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

  end subroutine bounce_on_path_thin2

  subroutine bounce_on_path_thin(this)

    class(solver), intent(inout) :: this

    real(wp) :: x_match
    real(wp) :: msq_true
    real(wp) :: d2V
    logical :: shooting_converged
    real(wp) :: nu
!   real(wp) :: x_min
!   real(wp) :: x_max
    real(wp) :: rho_min
    real(wp) :: rho_max
!   real(wp) :: x1
    real(wp) :: x0(2)
    integer :: over_under_flag
    logical :: overshot
    integer :: n
    real(wp), allocatable :: rho(:)
    real(wp), allocatable :: x(:, :)
    integer :: i
    real(wp) :: r
    real(wp) :: rho_match
    integer :: i_match
    integer :: i_match_min
    integer :: i_match_max
    real(wp) :: C_min
    real(wp) :: C_max
    real(wp) :: C
    real(wp) :: b
    real(wp) :: drho
    real(wp), allocatable :: rho_ode(:)
    real(wp), allocatable :: x_ode(:, :)

    n = this%n

    nu = 0.5e0_wp * (this%alpha - 1.0e0_wp)

    msq_true = this%msq_true ! this%d2V_dx2(1)
!   do i = 1, n
!     d2V = this%d2V_dx2(i)
!     ! At the bubble wall d2V = 0, abs(d2V - msq_true) ~ msq_true
!     ! I want to use thin-walled approximation not until bubble wall,
!     ! but only a subrange below bubble wall, so added factor 2.
!     if (4.0e0_wp * abs(d2V - msq_true) > msq_true) then
!       x_match = this%x(i)
!       exit
!     end if
!   end do
!   if ((i - 1) == n) then
!     write(*,*) "x_match problem"
!     call exit
!   end if
!   write(*,*) "asdfasd", i

!   do i = 1, n
!     if (norm2(this%phi(i, :) - this%phi_true) >  &
!         1.0e-2_wp * norm2(this%phi_false - this%phi_true)) then
!       x_match = this%x(i)
!       i_match = i
!       exit
!     end if
!   end do
!   write(*,*) "asdfasd", i, x_match

!   if (i_match == 0) then
!     write(*,*) "Couldn't determine i_match"
!     write(*,*) i_match, x_match
!     call exit
!   end if

!     write(*,*) i_match, rho_match, x_match, C

    shooting_converged = .false.
    over_under_flag = 0
    C_min = 1.0e-30_wp
    C_max = 1.0e-1_wp ! Determine from position of barrier
    b = sqrt(msq_true)

    rho_min = 1.0e-4_wp
    rho_max = 1.0e0_wp * this%rho_max

    rho = linspace(1.0e-4_wp, rho_max, n)
!   rho_match = rho(i_match)
    drho = rho(2) - rho(1)
    allocate(x(n, 2))

    C = 1.0e-5_wp
    i_match_min = 2
    i_match_max = n / 2

    do while (.not. shooting_converged)

!     select case (over_under_flag)
!       case (0) ! First iteration
!         C = (C_max + C_min) / 100.0e0_wp
!       case (1) ! Overshoot
!         C_max = C
!         C = (C_max + C_min) / 2.0e0_wp
!       case (-1) ! Undershoot
!         C_min = C
!         C = (C_max + C_min) / 2.0e0_wp
!     end select
!     C = 0.00000000001e0_wp

      select case (over_under_flag)
        case (0) ! First iteration
          i_match = (i_match_max + i_match_min) / 2.0e0_wp
        case (1) ! Overshoot
          i_match_max = i_match
          i_match = (i_match_max + i_match_min) / 2.0e0_wp
        case (-1) ! Undershoot
          i_match_min = i_match
          i_match = (i_match_max + i_match_min) / 2.0e0_wp
      end select

      ! Analytic solution until x = x_match
!     do i = 1, i_match
!       x(i, 1) = C * Inu(rho(i))
!       if (x(i, 1) > x_match) then
!         rho_match = rho(i)
!         exit
!       end if
!     end do
      do i = 1, i_match
        x(i, 1) = C * Inu(rho(i))
      end do

      ! Get xdot for analytic part
      x(1, 2) = (x(2, 1) - x(1, 1)) / drho
      do i = 2, i_match - 1
        x(i, 2) = (x(i + 1, 1) - x(i - 1, 1)) / (2.0e0_wp * drho)
      end do
      x(i_match, 2) = (x(i_match, 1) - x(i_match - 1, 1)) / drho
!     write(*,*) x(1:i_match, 2)

      ! TODO: make tighter cut on rho if needed

      ! Integrate the remaining rho region with odeint
      x0 = x(i_match, :)
      rho_match = rho(i_match)
      call integrate(  &
        dxdrho,  &
        x0,  &
        rho_match + drho,  &
        rho_max,  &
        n - i_match,  &
        rho_ode, x_ode)

      ! Stitching together thin + odeint parts
      rho(i_match+1:n) = rho_ode
      x(i_match+1:n, :) = x_ode

!     do i = 570, 600
!       write(*,*) i, rho(i), x(i, 2)
!     end do
!     call exit

      if (overshot) then
        over_under_flag = 1
      else
        over_under_flag = -1
      end if

!     if (is_equal(  &
!           C_min / C_max,  &
!           1.0e0_wp,  &
!           eps=1.0e-10_wp)) then
!       shooting_converged = .true.
!     else
!       shooting_converged = .false.
!     end if

      if (abs(i_match_min -i_match_max) <= 1) then
        shooting_converged = .true.
      else
        shooting_converged = .false.
      end if

!     if (this%verbose_level >= 2) then
!       write(*,*) C_min, C_max, over_under_flag
!     end if

      if (this%verbose_level >= 2) then
        write(*,*) i_match_min, i_match_max, over_under_flag
      end if

!     shooting_converged = .true.

    end do

    this%rho = rho
    this%xb = x(:, 1)
    this%xbdot = x(:, 2)
    do i = 1, n
      call this%get_index_from_x(  &
        this%xb(i),  &
        this%ixs(i))
    end do

    call this%calc_phi_on_bounce()

    do i = 1, n
      this%pot(i) = this%V(this%phi(i, :))
    end do

  contains

    function Inu(x) result(y)

      real(wp), intent(in) :: x
      real(wp) :: y

      if (this%alpha == 2) then
        y = sinh(x * b) / x
      else if (this%alpha == 3) then
        y = bessel_mod_1(x * b) / x
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

  subroutine calc_msq_true(this)

    class(solver), intent(inout) :: this

    this%msq_true = this%d2V_dx2(1)

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
        write(*,"(a,I5,a,ES12.3)")  &
          "Clipped solution" // &
          "at rho(",  i - 1, ") = ",  &
          this%rho(i - 1)
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
      gradV(i, :) = this%gradV(phi(i, :))
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
!       call exit
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
        if (ieee_is_nan(phi(i, j))) then
          write(*,*) i, j, phi(i, j), this%iflag_bspls(j)
          call exit
        end if
      end do
      if (this%iflag_bspls(j) /= 0) then
        write(*,*) 'Error evaluating ', j, 'D spline: ' //  &
          get_status_message(this%iflag_bspls(j))
        call exit
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
        if (i == 1) then
          write(*,*) "spline path went wrong"
          write(*,*) i, dphidx
          write(*,*) x(i + 1), x(i)
          write(*,*) phi(i + 1, :)
          write(*,*) phi(i, :)
          write(*,*) x(i + 1) - x(i)
          write(*,*) norm2(phi(i + 1, :) - phi(i, :))
        end if
        this%spline_fail = .true.
!       call exit
      end if
    end do

    this%x = x
    this%phi = phi
    do i = 1, n
      this%pot(i) = this%V(this%phi(i, :))
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
    write(*,*) "Force ratio =", this%Rforce

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
      this%ib_max_best = this%ib_max
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
      write(*,*) "Maximum iterations reached:", this%maximum_iterations
    end if

!   if (is_equal(this%Rforce / this%Rforce_previous, 1.0e0_wp, eps=1.0e-2_wp)) then
!     this%deform_eps = 1.2e0_wp * this%deform_eps
!     if (this%verbose_level >= 1) then
!       write(*,*) "Slow convergence: Increased deform_eps to", this%deform_eps
!     end if
!   end if

    if (this%Rforce / this%Rforce_previous > 1.2e0_wp) then
      this%deform_eps = 0.5e0_wp * this%deform_eps
      if (this%verbose_level >= 1) then
        write(*,*) "Bad deformation: Decrease deform_eps to", this%deform_eps
      end if
    end if

    this%Rforce_previous = this%Rforce

  end subroutine check_convergence

end module bouncesolver__pdm
