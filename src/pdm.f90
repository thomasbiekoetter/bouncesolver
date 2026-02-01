module bouncesolver__pdm

  use bouncesolver__config, only : wp
  use bouncesolver__util, only : linspace
  use bouncesolver__util, only : pi
  use bouncesolver__util, only : riemann_integrate
  use bouncesolver__util, only : is_equal
  use odeint__rk4, only : integrate
  use gradmin__derivatives, only : gradient
  use bspline_module, only : bspline_1d

  implicit none

  private

  real(wp), parameter :: eps_gradient = 1.0e-8_wp

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
    integer :: current_iteration
    integer, allocatable :: ixs(:)
    real(wp) :: gradV_max
    real(wp), allocatable :: phi_deformed(:, :)
    type(bspline_1d), allocatable :: bspls(:)
    integer, allocatable :: iflag_bspls(:)
    integer :: kx_bspls
    logical :: spline_fail = .false.
  contains
    procedure, private :: allocate_arrays
    procedure, private :: construct_starting_path
    procedure, private :: find_x_barrier
    procedure, private :: d2V_dx2
    procedure, private :: estimate_rho_max
    procedure, private :: get_index_from_x
    procedure, private :: bounce_on_path
    procedure, private :: dV_dx
    procedure, private :: clip_bounce
    procedure, private :: calc_action
    procedure, private :: calc_forces
    procedure, private :: gradV
    procedure, private :: calc_phi_on_bounce
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
    maxiter, deform_eps,  &
    verbose_level) result(this)

    integer, intent(in) :: d
    procedure(V_abstract) :: V
    real(wp), intent(in) :: phi_false(d)
    real(wp), intent(in) :: phi_true(d)
    integer, intent(in), optional :: alpha
    real(wp), intent(in), optional :: rho_max_fac
    real(wp), intent(in), optional :: x_min
    integer, intent(in), optional :: maxiter
    real(wp), intent(in), optional :: deform_eps
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

    if (present(maxiter)) this%maximum_iterations = maxiter

    if (present(deform_eps)) this%deform_eps = deform_eps

    if (present(verbose_level)) this%verbose_level = verbose_level

    this%kx_bspls = 5

    call this%allocate_arrays()
    call this%construct_starting_path()
    call this%estimate_rho_max()
    converged = .false.
    iteration = 0
    do while (.not. converged)
      iteration = iteration + 1
      call this%bounce_on_path()
      call this%calc_action()
      write(*,*) "Action = ", this%S, this%Sp, this%Sk
      if (iteration == 1) this%S0 = this%S
      call this%calc_forces(iteration)
      call this%deform_path()
      if (this%spline_fail) exit
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

    this%rho_max = this%rho_max_fac * sqrt(msq_false)

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

!     write(*,*) x_min, x_max, over_under_flag

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
        phib(i, :) = phi(i, :)
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
        pathdir(i, :) = (phi(2, :) - phi(1, :)) / (x(2) - x(1))
      else if (i == n) then
        pathdir(i, :) = (phi(n, :) - phi(n - 1, :)) / (x(n) - x(n - 1))
      else
        pathdir(i, :) = (phi(i + 1, :) - phi(i - 1, :)) / (x(i + 1) - x(i - 1))
      end if
      if (.not. is_equal(norm2(pathdir(i, :)), 1.0e0_wp, eps=1.0e-3_wp)) then
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

    ! Construct B-spline approximation to get path with
    ! this%n points and for more precise calculation
    ! of d2phi_dx2
    do j = 1, d
      call this%bspls(j)%initialize(  &
        x_spline,  &
        phi_spline(:, j), this%kx_bspls, this%iflag_bspls(j))
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
      write(*,*) "Stopped before reaching Rforce_threshold"
      write(*,*) "because Rforce started to increase."
    end if

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
