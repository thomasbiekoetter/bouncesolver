module bouncesolver__pathdeformation

  use bouncesolver__config, only : wp
  use bouncesolver__util, only : linspace
  use bouncesolver__util, only : pi
  use bouncesolver__util, only : riemann_integrate
  use bouncesolver__util, only : is_equal
  use odeint__rk4, only : integrate
  use evortran__util_interp_spline, only :  &
    spline_construct
  use evortran__util_interp_spline, only :  &
    spline_getval
  use gradmin__derivatives, only : gradient

  implicit none

  private

  real(wp), parameter :: eps_gradient = 1.0e-6_wp
  real(wp), parameter :: eps_hessian = 1.0e-4_wp

  real(wp), parameter :: xmin_start = 1.0e-18_wp

  abstract interface
    function V_abstract(x) result(y)
      import wp
      real(wp), intent(in) :: x(:)
      real(wp) :: y
    end function V_abstract
  end interface

  type, public :: solver
    procedure(V_abstract), pointer, nopass :: V => null()
    integer :: num_fields
    integer :: alpha
    integer :: nsteps_odeint = 10000
    integer :: nnodes_phi_of_x_init = 10000
    integer :: nnodes_phi_of_x
    real(wp), allocatable :: phi_false(:)
    real(wp), allocatable :: phi_true(:)
    real(wp), allocatable :: x_init(:)
    real(wp), allocatable :: phi_init(:, :)
    real(wp), allocatable :: phi_of_x_init_b(:, :)
    real(wp), allocatable :: phi_of_x_init_c(:, :)
    real(wp), allocatable :: phi_of_x_init_d(:, :)
    real(wp) :: x_init_min
    real(wp) :: x_init_max
    real(wp), allocatable :: phi(:, :)
    real(wp), allocatable :: x(:)
    real(wp), allocatable :: phi_of_x_b(:, :)
    real(wp), allocatable :: phi_of_x_c(:, :)
    real(wp), allocatable :: phi_of_x_d(:, :)
    real(wp) :: x_min
    real(wp) :: x_max
    real(wp) :: rho_min
    real(wp) :: rho_max
    real(wp) :: rho_max_fac
    real(wp) :: x_barrier
    real(wp) :: B_approx
    real(wp) :: xmin_bounce
    real(wp), allocatable :: rho_bounce(:)
    real(wp), allocatable :: x_bounce(:)
    real(wp), allocatable :: x_of_rho_bounce_b(:)
    real(wp), allocatable :: x_of_rho_bounce_c(:)
    real(wp), allocatable :: x_of_rho_bounce_d(:)
    real(wp), allocatable :: xdot_bounce(:)
    real(wp), allocatable :: xdot_of_rho_bounce_b(:)
    real(wp), allocatable :: xdot_of_rho_bounce_c(:)
    real(wp), allocatable :: xdot_of_rho_bounce_d(:)
    integer :: i_bounce_clip    ! Action is integrated up to here
    real(wp) :: rho_bounce_clip
    real(wp) :: S
    real(wp) :: Sp
    real(wp) :: Sk
    real(wp), allocatable :: normal_force(:, :)
  contains
    procedure, private :: allocate_arrays
    procedure, private :: construct_starting_path
    procedure, public :: phi_of_x_init
    procedure, public :: phi_of_x
    procedure, public :: V_of_x
    procedure, public :: dV_dx
    procedure, public :: d2V_dx2
    procedure, private :: estimate_rho_max
    procedure, private :: set_rho_min
    procedure, public :: bounce_on_path
    procedure, private :: find_x_barrier
    procedure, private :: calc_B_approx
    procedure, private :: check_over_under
    procedure, public :: x_of_rho
    procedure, public :: xdot_of_rho
    procedure, public :: d2x_of_drho2
    procedure, public :: V_of_rho
    procedure, public :: phi_of_rho
    procedure, public :: dphi_drho
    procedure, private :: clip_solution
    procedure, private :: calc_action
    procedure, private :: calc_normal_forces
    procedure, public :: grad_V
    procedure, private :: dphi_dx
  end type solver

  interface solver
    module procedure create_solver
  end interface solver

contains

  function create_solver(  &
    num_fields, V, phi_false, phi_true,  &
    alpha, rho_max_fac, xmin) result(this)

    integer, intent(in) :: num_fields
    procedure(V_abstract) :: V
    real(wp), intent(in) :: phi_false(num_fields)
    real(wp), intent(in) :: phi_true(num_fields)
    integer, intent(in), optional :: alpha
    real(wp), intent(in), optional :: rho_max_fac
    real(wp), intent(in), optional :: xmin
    type(solver) :: this

    this%num_fields = num_fields
    this%V => V
    this%phi_false = phi_false
    this%phi_true = phi_true
    if (present(alpha)) then
      this%alpha = alpha
    else
      this%alpha = 2
    end if
    if (present(rho_max_fac)) then
      this%rho_max_fac = rho_max_fac
    else
      this%rho_max_fac = 20.0e0_wp
    end if
    if (present(xmin)) then
      this%xmin_bounce = xmin
    else
      this%xmin_bounce = xmin_start
    end if

    call this%allocate_arrays()

    call this%construct_starting_path()

    call this%estimate_rho_max()

    call this%set_rho_min()

    call this%find_x_barrier()

    call this%calc_B_approx()

    ! This will have to be a loop
    ! do while (normal forces not zero)
        call this%bounce_on_path()
        call this%clip_solution()
        call this%calc_normal_forces()
    ! call this%deform_path
    ! end do

    call this%calc_action()

  end function create_solver

  subroutine allocate_arrays(this)

    class(solver), intent(inout) :: this

    ! Allocate data arrays for cubic
    ! interpolation of initial path
    allocate(this%x_init(  &
      this%nnodes_phi_of_x_init))
    allocate(this%phi_init(  &
      this%nnodes_phi_of_x_init,  &
      this%num_fields))
    allocate(this%phi_of_x_init_b(  &
      this%nnodes_phi_of_x_init,  &
      this%num_fields))
    allocate(this%phi_of_x_init_c(  &
      this%nnodes_phi_of_x_init,  &
      this%num_fields))
    allocate(this%phi_of_x_init_d(  &
      this%nnodes_phi_of_x_init,  &
      this%num_fields))

    ! Allocate data array for x of rho
    allocate(this%rho_bounce(this%nsteps_odeint))
    allocate(this%x_bounce(this%nsteps_odeint))
    allocate(this%x_of_rho_bounce_b(this%nsteps_odeint))
    allocate(this%x_of_rho_bounce_c(this%nsteps_odeint))
    allocate(this%x_of_rho_bounce_d(this%nsteps_odeint))

    ! Allocate data array for xdot of rho
    allocate(this%xdot_bounce(this%nsteps_odeint))
    allocate(this%xdot_of_rho_bounce_b(this%nsteps_odeint))
    allocate(this%xdot_of_rho_bounce_c(this%nsteps_odeint))
    allocate(this%xdot_of_rho_bounce_d(this%nsteps_odeint))

  end subroutine allocate_arrays

  ! Constructs an initial path as a straight
  ! line phi_i(x) that connects the true
  ! to the false minimum under the condition
  ! that x is the length of the tunneling path
  ! in field space: |dphi_i / dx| = 1
  subroutine construct_starting_path(this)

    class(solver), intent(inout) :: this

    integer :: num_nodes
    real(wp) :: x_min
    real(wp) :: x_max
    real(wp), allocatable :: x(:)
    real(wp), allocatable :: dx(:)
    real(wp), allocatable :: t(:)
    real(wp), allocatable :: phi(:, :)
    integer :: i
    real(wp) :: dphidx

    ! Construct first straight path
    ! without requiring |dphi / dx| = 1.
    ! s will later be x after inteprolation
    num_nodes = this%nnodes_phi_of_x_init
    x_min = 0.0e0_wp
    x_max = 1.0e0_wp
    allocate(phi(num_nodes, this%num_fields))

    x = linspace(x_min, x_max, num_nodes)

    t = (x - x_min) / (x_max - x_min)
    do i = 1, num_nodes
      phi(i, :) = this%phi_true +  &
        t(i) * (this%phi_false - this%phi_true)
    end do

    ! Rescale to have |dphi / dx | = 1.
    allocate(dx(num_nodes - 1))
    do i = 1, num_nodes - 1
      dx(i) = sqrt(sum((phi(i + 1, :) - phi(i, :)) ** 2))
    end do
    do i = 2, num_nodes
      x(i) = x(i - 1) + dx(i - 1)
    end do
    ! Check
    do i = 1, num_nodes - 1
      dphidx = norm2(phi(i + 1, :) - phi(i, :)) /  &
        (x(i + 1) - x(i))
      if (.not. is_equal(abs(dphidx), 1.0e0_wp)) then
        write(*,*) "init path went wrong"
        call exit
      end if
    end do

    ! Make cubic spline interpolation
    this%x_init = x
    this%x_init_min = x(1)
    this%x_init_max = x(num_nodes)
    this%phi_init = phi
    do i = 1, this%num_fields
      call spline_construct(  &
        this%x_init, this%phi_init(:, i),  &
        this%phi_of_x_init_b(:, i),  &
        this%phi_of_x_init_c(:, i),  &
        this%phi_of_x_init_d(:, i),  &
        num_nodes)
    end do

    ! Set current path to init path
    this%x = this%x_init
    this%x_min = this%x_init_min
    this%x_max = this%x_init_max
    this%phi = this%phi_init
    this%phi_of_x_b = this%phi_of_x_init_b
    this%phi_of_x_c = this%phi_of_x_init_c
    this%phi_of_x_d = this%phi_of_x_init_d
    this%nnodes_phi_of_x = this%nnodes_phi_of_x_init

  end subroutine construct_starting_path

  function phi_of_x_init(this, x) result(phi)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: x
    real(wp) :: phi(this%num_fields)

    integer :: i

    do i = 1, this%num_fields
      phi(i) = spline_getval(  &
        x,  &
        this%x_init,  &
        this%phi_init(:, i),  &
        this%phi_of_x_init_b(:, i),  &
        this%phi_of_x_init_c(:, i),  &
        this%phi_of_x_init_d(:, i),  &
        this%nnodes_phi_of_x_init)
    end do

  end function phi_of_x_init

  function phi_of_x(this, x) result(phi)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: x
    real(wp) :: phi(this%num_fields)

    integer :: i

    do i = 1, this%num_fields
      phi(i) = spline_getval(  &
        x,  &
        this%x,  &
        this%phi(:, i),  &
        this%phi_of_x_b(:, i),  &
        this%phi_of_x_c(:, i),  &
        this%phi_of_x_d(:, i),  &
        this%nnodes_phi_of_x)
    end do

  end function phi_of_x

  function V_of_x(this, x) result(V)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: x
    real(wp) :: V

    real(wp) :: phi(this%num_fields)

    phi = this%phi_of_x(x)
    V = this%V(phi)

  end function V_of_x

  function dV_dx(this, x, h) result(dV)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: x
    real(wp), optional :: h
    real(wp) :: dV

    real(wp) :: phi(this%num_fields)
    real(wp) :: phi_m2(this%num_fields)
    real(wp) :: phi_m1(this%num_fields)
    real(wp) :: phi_p1(this%num_fields)
    real(wp) :: phi_p2(this%num_fields)
    real(wp) :: V
    real(wp) :: V_m2
    real(wp) :: V_m1
    real(wp) :: V_p1
    real(wp) :: V_p2
    real(wp) :: eps
    logical :: lower_ok
    logical :: upper_ok
    real(wp) :: blend
    real(wp) :: dV_blend
    real(wp) :: delta_phi_true(this%num_fields)

    phi = this%phi_of_x(x)
    V = maxval([abs(this%V(phi)), 1.0e0_wp])

    if (present(h)) then
      eps = h * V
    else
      eps = eps_gradient * V
    end if

    lower_ok = (x - 2.0e0_wp * eps) > this%x_min
    upper_ok = (x + 2.0e0_wp * eps) < this%x_max

    if (lower_ok .and. upper_ok) then
      ! central
      phi_m2 = this%phi_of_x(x - 2.0e0_wp * eps)
      phi_m1 = this%phi_of_x(x - 1.0e0_wp * eps)
      phi_p1 = this%phi_of_x(x + 1.0e0_wp * eps)
      phi_p2 = this%phi_of_x(x + 2.0e0_wp * eps)
      V_m2 = this%V(phi_m2)
      V_m1 = this%V(phi_m1)
      V_p1 = this%V(phi_p1)
      V_p2 = this%V(phi_p2)
      dV = (-V_p2 + 8.0e0_wp * V_p1 -  &
        8.0e0_wp * V_m1 + V_m2) / (12.0e0_wp * eps)
    else if (upper_ok) then
      ! forward
      phi_m2 = this%phi_of_x(x)
      phi_m1 = this%phi_of_x(x + 1.0e0_wp * eps)
      phi_p1 = this%phi_of_x(x + 2.0e0_wp * eps)
      phi_p2 = this%phi_of_x(x + 3.0e0_wp * eps)
      V_m2 = this%V(phi_m2)
      V_m1 = this%V(phi_m1)
      V_p1 = this%V(phi_p1)
      V_p2 = this%V(phi_p2)
      dV = (-11.0e0_wp * V_m2 + 18.0e0_wp * V_m1 -  &
        9.0e0_wp * V_p1 + 2.0e0_wp * V_p2) /  &
        (6.0e0_wp * eps)
    else if (lower_ok) then
      phi_m2 = this%phi_of_x(x - 4.0e0_wp * eps)
      phi_m1 = this%phi_of_x(x - 3.0e0_wp * eps)
      phi = this%phi_of_x(x - 2.0e0_wp * eps)
      phi_p1 = this%phi_of_x(x - 1.0e0_wp * eps)
      phi_p2 = this%phi_of_x(x)
      V_m2 = this%V(phi_m2)
      V_m1 = this%V(phi_m1)
      V = this%V(phi)
      V_p1 = this%V(phi_p1)
      V_p2 = this%V(phi_p2)
      dV = (25.0e0_wp * V_p2 - 48.0e0_wp * V_p1 +  &
        36.0e0_wp * V - 16.0e0_wp * V_m1 +  &
        3.0e0_wp * V_m2) / (12.0e0_wp * eps)
    else
      write(*,*) "Problem in dV_dx"
      call exit
    end if

    delta_phi_true = this%phi_of_x(x) - this%phi_true
    if (norm2(delta_phi_true) < eps * 100.0e0_wp) then
      dV_blend = this%d2V_dx2(x) * (x - this%x(1))
      blend = exp(-(norm2(delta_phi_true) / eps) ** 2)
      dV = dV * (1.0e0_wp - blend) + dV_blend * blend
    end if

  end function dV_dx

  function d2V_dx2(this, x, h) result(d2V)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: x
    real(wp), optional :: h
    real(wp) :: d2V

    real(wp) :: phi_m2(this%num_fields)
    real(wp) :: phi_m1(this%num_fields)
    real(wp) :: phi(this%num_fields)
    real(wp) :: phi_p1(this%num_fields)
    real(wp) :: phi_p2(this%num_fields)
    real(wp) :: V_m2
    real(wp) :: V_m1
    real(wp) :: V
    real(wp) :: V_p1
    real(wp) :: V_p2
    real(wp) :: eps
    logical :: lower_ok
    logical :: upper_ok

    phi = this%phi_of_x(x)
    V = maxval([abs(this%V(phi)), 1.0e0_wp])

    if (present(h)) then
      eps = h * V
    else
      eps = eps_hessian * V
    end if

    lower_ok = (x - 2.0e0_wp * eps) > this%x_min
    upper_ok = (x + 2.0e0_wp * eps) < this%x_max

    if (lower_ok .and. upper_ok) then
      ! central
      phi_m2 = this%phi_of_x(x - 2.0e0_wp * eps)
      phi_m1 = this%phi_of_x(x - 1.0e0_wp * eps)
      phi = this%phi_of_x(x)
      phi_p1 = this%phi_of_x(x + 1.0e0_wp * eps)
      phi_p2 = this%phi_of_x(x + 2.0e0_wp * eps)
      V_m2 = this%V(phi_m2)
      V_m1 = this%V(phi_m1)
      V = this%V(phi)
      V_p1 = this%V(phi_p1)
      V_p2 = this%V(phi_p2)
      d2V = (-V_m2 + 16.0e0_wp * V_m1 -  &
        30.0e0_wp * V + 16.0e0_wp * V_p1 -  &
        V_p2) / (12.e0_wp * eps ** 2)
    else if (upper_ok) then
      ! forward
      phi_m2 = this%phi_of_x(x)
      phi_m1 = this%phi_of_x(x + 1.0e0_wp * eps)
      phi = this%phi_of_x(x + 2.0e0_wp * eps)
      phi_p1 = this%phi_of_x(x + 3.0e0_wp * eps)
      phi_p2 = this%phi_of_x(x + 4.0e0_wp * eps)
      V_m2 = this%V(phi_m2)
      V_m1 = this%V(phi_m1)
      V = this%V(phi)
      V_p1 = this%V(phi_p1)
      V_p2 = this%V(phi_p2)
      d2V = (35.0e0_wp * V_m2 - 104.0e0_wp * V_m1 +  &
        114.0e0_wp * V - 56.0e0_wp * V_p1 +  &
        11.0e0_wp * V_p2) / (12.0e0_wp * eps ** 2)
    else if (lower_ok) then
      phi_m2 = this%phi_of_x(x - 4.0e0_wp * eps)
      phi_m1 = this%phi_of_x(x - 3.0e0_wp * eps)
      phi = this%phi_of_x(x - 2.0e0_wp * eps)
      phi_p1 = this%phi_of_x(x - 1.0e0_wp * eps)
      phi_p2 = this%phi_of_x(x)
      V_m2 = this%V(phi_m2)
      V_m1 = this%V(phi_m1)
      V = this%V(phi)
      V_p1 = this%V(phi_p1)
      V_p2 = this%V(phi_p2)
      d2V = (35.0e0_wp * V_p2 - 104.0e0_wp * V_p1 +  &
        114.0e0_wp * V - 56.0e0_wp * V_m1 +  &
        11.0e0_wp * V_m2) / (12.0e0_wp * eps ** 2)
    else
      write(*,*) "Problem in dV_dx"
      call exit
    end if

  end function d2V_dx2

  ! Estimate rho_max as 15 * m_false with
  ! m_false ** 2 = d2V_true (curvature =
  ! mass of particle in false minimum)
  subroutine estimate_rho_max(this)

    class(solver), intent(inout) :: this

    real(wp) :: msq_false

    msq_false = this%d2V_dx2(this%x_max)
    this%rho_max = this%rho_max_fac * sqrt(msq_false)

  end subroutine estimate_rho_max

  ! For now just make small
  subroutine set_rho_min(this)

    class(solver), intent(inout) :: this

    this%rho_min = 1.0e-4_wp

  end subroutine set_rho_min

  subroutine bounce_on_path(this)

    class(solver), intent(inout) :: this

    logical :: shooting_converged

    real(wp), allocatable :: x(:, :)
    real(wp), allocatable :: rho(:)
    real(wp) :: xmin
    real(wp) :: xmax
    real(wp) :: x1
    real(wp) :: x0(2)
    integer :: num
    integer :: over_under_flag

    xmin = this%xmin_bounce
    xmax = this%x_barrier
    num = this%nsteps_odeint
    shooting_converged = .false.
    over_under_flag = 0

    do while (.not. shooting_converged)

      select case (over_under_flag)
        case (0) ! First iteration
          x1 = (xmax + xmin) / 100.0e0_wp
        case (1) ! Overshoot
          xmin = x1
          x1 = (xmax + xmin) / 2.0e0_wp
        case (-1) ! Undershoot
          xmax = x1
          x1 = (xmax + xmin) / 2.0e0_wp
      end select

      x0(1) = x1
      x0(2) = 0.0e0_wp

      call integrate(  &
        dxdrho,  &
        x0,  &
        this%rho_min,  &
        this%rho_max,  &
        num,  &
        rho, x)

      call this%check_over_under(  &
        rho, x,  &
        over_under_flag)

      if (is_equal(xmin / xmax, 1.0e0_wp, eps=1.0e-10_wp)) then
        shooting_converged = .true.
      else
        shooting_converged = .false.
      end if

!     write(*,*) xmin, xmax

    end do

    this%rho_bounce = rho
    this%x_bounce = x(:, 1)
    call spline_construct(  &
      this%rho_bounce, this%x_bounce,  &
      this%x_of_rho_bounce_b,  &
      this%x_of_rho_bounce_c,  &
      this%x_of_rho_bounce_d,  &
      num)

    this%xdot_bounce = x(:, 2)
    call spline_construct(  &
      this%rho_bounce, this%xdot_bounce,  &
      this%xdot_of_rho_bounce_b,  &
      this%xdot_of_rho_bounce_c,  &
      this%xdot_of_rho_bounce_d,  &
      num)

    if (is_equal(xmax / xmin_start, 1.0e0_wp, eps=1.0e-10_wp)) then
      write(*,*) "Warning: Not converged correctly. Wall too thin."
      write(*,*) "You can try again with larger rho_max_fac (default = 20)."
    end if

  contains

    function dxdrho(xi, rho) result(xdot)

      real(wp), intent(in) :: xi(:)
      real(wp), intent(in) :: rho
      real(wp), allocatable :: xdot(:)

      allocate(xdot(size(xi)))

      xdot(1) = xi(2)

      if ((rho < 1.0e-1_wp) .and. (xi(1) < 1.0e-1_wp)) then
        ! Not sure if this helps, maybe remove
        xdot(2) = this%B_approx * xi(1) -  &
          this%alpha * xi(2) / rho
      else
        xdot(2) = this%dV_dx(xi(1)) -  &
          this%alpha * xi(2) / rho
      end if

    end function dxdrho

  end subroutine bounce_on_path

  subroutine find_x_barrier(this)

    class(solver), intent(inout) :: this

    real(wp), allocatable :: x(:)
    real(wp) :: V1
    real(wp) :: V2
    integer :: i

    x = linspace(this%x_min, this%x_max, this%nsteps_odeint)

    V1 = this%V_of_x(x(1))
    do i = 2, this%nsteps_odeint
      V2 = this%V_of_x(x(i))
      if (V2 > V1) then
        V1 = V2
      else
        this%x_barrier = x(i - 1)
        exit
      end if
    end do

    if (i == this%nsteps_odeint) then
      write(*,*) "No potential barrier on path."
      call exit
    end if

  end subroutine find_x_barrier

  subroutine calc_B_approx(this)

    class(solver), intent(inout) :: this

    real(wp) :: dx
    real(wp) :: x0
    real(wp) :: x1
    real(wp) :: x2
    real(wp) :: x3
    real(wp) :: x4
    real(wp) :: phi0(this%num_fields)
    real(wp) :: phi1(this%num_fields)
    real(wp) :: phi2(this%num_fields)
    real(wp) :: phi3(this%num_fields)
    real(wp) :: phi4(this%num_fields)
    real(wp) :: V0
    real(wp) :: V1
    real(wp) :: V2
    real(wp) :: V3
    real(wp) :: V4
    real(wp) :: dphi
    integer :: i
    real(wp) :: d2V

    dx = this%x(2) - this%x(1)
    x0 = 0.0e0_wp
    x1 = x0 + dx
    x2 = x1 + dx
    x3 = x2 + dx
    x4 = x3 + dx

    phi0 = this%phi_of_x(x0)
    phi1 = this%phi_of_x(x1)
    phi2 = this%phi_of_x(x2)
    phi3 = this%phi_of_x(x3)
    phi4 = this%phi_of_x(x4)

    dphi = norm2(phi1 - phi0) ! sign doesn't matter because
                              ! h ** 2 in denominator of 2nd derivative

    ! Since |dphi / dx| = 1 the dphis should be
    ! equidistant if the dxs are. Check:
    if (.not. is_equal(norm2(phi2 - phi1), dphi)) then
      write(*,*) "Problem with dphi in calc_B_approx."
      call exit
    end if
    if (.not. is_equal(norm2(phi3 - phi2), dphi)) then
      write(*,*) "Problem with dphi in calc_B_approx."
      call exit
    end if
    if (.not. is_equal(norm2(phi4 - phi3), dphi)) then
      write(*,*) "Problem with dphi in calc_B_approx."
      call exit
    end if

    V0 = this%V(phi0)
    V1 = this%V(phi1)
    V2 = this%V(phi2)
    V3 = this%V(phi3)
    V4 = this%V(phi4)

    ! 2nd derivative of V wrt. phi along path
    d2V = (35.0e0_wp * V0 - 104.0e0_wp * V1 +  &
      114e0_wp * V2 - 56.0e0_wp * V3 +  &
      11.0e0_wp * V4) / (12.0e0_wp * dphi ** 2)

    this%B_approx = d2V

  end subroutine calc_B_approx

  subroutine check_over_under(this, rho, x, flag)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: rho(this%nsteps_odeint)
    real(wp), intent(in) :: x(this%nsteps_odeint, 2)
    integer, intent(out) :: flag

    ! Overshoot if x goes beyond the false minimum,
    ! undershoot otherwise
    if (maxval(x(:, 1)) > this%x_max + 1.0e-4) then
      flag = 1
    else
      flag = -1
    end if

  end subroutine check_over_under

  function x_of_rho(this, rho) result(x)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: rho
    real(wp) :: x

    x = spline_getval(  &
      rho,  &
      this%rho_bounce,  &
      this%x_bounce,  &
      this%x_of_rho_bounce_b,  &
      this%x_of_rho_bounce_c,  &
      this%x_of_rho_bounce_d,  &
      this%nsteps_odeint)

  end function x_of_rho

  function xdot_of_rho(this, rho) result(xdot)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: rho
    real(wp) :: xdot

    xdot = spline_getval(  &
      rho,  &
      this%rho_bounce,  &
      this%xdot_bounce,  &
      this%xdot_of_rho_bounce_b,  &
      this%xdot_of_rho_bounce_c,  &
      this%xdot_of_rho_bounce_d,  &
      this%nsteps_odeint)

  end function xdot_of_rho

  function d2x_of_drho2(this, rho) result(dxdot)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: rho
    real(wp) :: dxdot

    real(wp) :: h
    real(wp) :: xdot_m2
    real(wp) :: xdot_m1
    real(wp) :: xdot_p1
    real(wp) :: xdot_p2

    h = maxval(this%rho_bounce) * 1.0e-6_wp

    xdot_m2 = this%xdot_of_rho(rho)
    xdot_m1 = this%xdot_of_rho(rho + 1.0e0_wp * h)
    xdot_p1 = this%xdot_of_rho(rho + 2.0e0_wp * h)
    xdot_p2 = this%xdot_of_rho(rho + 3.0e0_wp * h)

    dxdot = (-11.0e0_wp * xdot_m2 + 18.0e0_wp * xdot_m1 -  &
      9.0e0_wp * xdot_p1 + 2.0e0_wp * xdot_p2) /  &
      (6.0e0_wp * h)

  end function d2x_of_drho2

  function V_of_rho(this, rho) result(V)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: rho
    real(wp) :: V

    real(wp) :: x

    x = this%x_of_rho(rho)
    V = this%V_of_x(x)

  end function V_of_rho

  function phi_of_rho(this, rho) result(phi)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: rho
    real(wp) :: phi(this%num_fields)

    real(wp) :: x

    x = this%x_of_rho(rho)
    phi = this%phi_of_x(x)

  end function phi_of_rho

  function dphi_dx(this, x) result(y)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: x
    real(wp) :: y(this%num_fields)

    real(wp) :: phi1(this%num_fields)
    real(wp) :: phi2(this%num_fields)
    real(wp) :: h(this%num_fields)
    integer :: i

    phi1 = this%phi_of_x(x)

    do i = 1, this%num_fields
      h(i) = maxval([  &
        abs(phi1(i)) * eps_gradient,  &
        1.0e-6_wp])
      phi2 = this%phi_of_x(x + h(i))
      y(i) = (phi2(i) - phi1(i)) / h(i)
    end do

  end function dphi_dx

  function dphi_drho(this, rho) result(y)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: rho
    real(wp) :: y(this%num_fields)

    real(wp) :: phi1(this%num_fields)
    real(wp) :: phi2(this%num_fields)
    real(wp) :: h

    h = this%rho_max / 1.0e5_wp

    phi1 = this%phi_of_rho(rho)
    phi2 = this%phi_of_rho(rho + h)

    y = (phi2 - phi1) / h

  end function dphi_drho

  subroutine calc_action(this)

    class(solver), intent(inout) :: this

    real(wp) :: D
    integer :: n
    integer :: i
    real(wp), allocatable :: rho(:)
    real(wp), allocatable :: pot_term(:)
    real(wp), allocatable :: kin_term(:)
    real(wp), allocatable :: fpre(:)
    real(wp), allocatable :: f(:)

    D = this%alpha + 1

    rho = this%rho_bounce(1:this%i_bounce_clip)
    n = size(rho)

    if (this%alpha == 2) then
      fpre = 4.0e0_wp * pi * rho ** 2
    else if (this%alpha == 3) then
      fpre = 2.0e0_wp * pi ** 2 * rho ** 3
    else
      write(*,*) "Wrong alpha value."
      call exit
    end if

    allocate(pot_term(n))
    do i = 1, n
      pot_term(i) = this%V_of_rho(rho(i))
    end do
    pot_term = pot_term - this%V_of_x(this%x_max)

    allocate(kin_term(n))
    do i = 1, n
      kin_term(i) = sum(this%dphi_drho(rho(i)) ** 2)
    end do
    kin_term = kin_term * 0.5e0_wp

    ! From potential term only
    f = fpre * 2.0e0_wp * pot_term / (2.0e0_wp - D)
    this%Sp = riemann_integrate(rho, f)

    ! From kinetic term only
    f = fpre * 2.0e0_wp * kin_term / D
    this%Sk = riemann_integrate(rho, f)

    ! All terms
    f = fpre * (kin_term + pot_term)
    this%S = riemann_integrate(rho, f)

    if ((.not. is_equal(this%Sp / this%S, 1.0e0_wp, eps=1.0e-2_wp)) .or.  &
        (.not. is_equal(this%Sk / this%S, 1.0e0_wp, eps=1.0e-2_wp))) then
      write(*,"(a,ES12.3,a,ES12.3,a,ES12.3,a)")  &
        "Problem with action: (S, Sp, Sk) = (",  &
        this%S, ", ", this%Sp, ", ", this%Sk, ")"
    end if

  end subroutine calc_action

  subroutine clip_solution(this)

    class(solver), intent(inout) :: this

    integer :: i
    integer :: i_clip
    real(wp) :: rho_clip
    integer :: clip_range
    integer :: n
    real(wp) :: xmax

    n = size(this%rho_bounce)
    clip_range = 0

    xmax = maxval(this%x_bounce)

    do i = 1, n
      if (is_equal(this%x_bounce(i), xmax, eps=1.0e-2_wp) .and.  &
          is_equal(this%xdot_bounce(i), 0.0e0_wp, eps=1.0e-3_wp)) then
        clip_range = clip_range + 1
      else
        clip_range = 0
      end if
      ! If for ten steps in a row the boundary condition
      ! at rho_end is satisfied, clip the solution at this
      ! index to avoid undershoot from good solution
      if (clip_range == 10)  then
        exit
      end if
    end do
    i_clip = minval([i, n])
    rho_clip =this%rho_bounce(i_clip)
    this%i_bounce_clip = i_clip
    this%rho_bounce_clip = rho_clip
    if (i_clip < n) then
      write(*,"(a,I5,a,ES12.3)") "Warning: clipped solution at rho(",  &
        i_clip, ") = ", rho_clip
    end if

  end subroutine clip_solution

  subroutine calc_normal_forces(this)

    class(solver), intent(inout) :: this

    integer :: i
    integer :: n
    real(wp), allocatable :: rho(:)
    real(wp), allocatable :: dx_drho(:)
    real(wp), allocatable :: d2phi_dx2(:, :)
    real(wp), allocatable :: grad_perp_V(:, :)
    real(wp), allocatable :: x(:)
    real(wp) :: phi(this%num_fields)
    real(wp) :: grad_V(this%num_fields)
    real(wp) :: pathdir(this%num_fields)

    n = this%i_bounce_clip
    rho = this%rho_bounce(1:n)

    ! Compute dx / drho along path
    allocate(dx_drho(n))
    do i = 1, n
      dx_drho(i) = this%xdot_of_rho(rho(i))
    end do

    ! Compute d2 phi_i / dx dx along path
    allocate(d2phi_dx2(n, this%num_fields))
    allocate(x(n))
    do i = 1, n
      x(i) = this%x_of_rho(rho(i))
      d2phi_dx2(i, :) = d2p_dx2(x(i))
    end do

    ! Compute grad_T V(phi_i) along path
    allocate(grad_perp_V(n, this%num_fields))
    do i = 1, n
      phi = this%phi_of_rho(rho(i))
      grad_V = this%grad_V(phi)
      pathdir = this%dphi_dx(x(i))
      grad_perp_V(i, :) = grad_V - dot_product(grad_V, pathdir) * pathdir
    end do

    ! Compute normal force along path
    allocate(this%normal_force(n, this%num_fields))
    do i = 1, n
      this%normal_force(i, :) = d2phi_dx2(i, :) * dx_drho(i) ** 2 -  &
        grad_perp_V(i, :)
    end do

  contains

    function d2p_dx2(x) result(d2p)

      real(wp), intent(in) :: x
      real(wp) :: d2p(this%num_fields)

      real(wp) :: h
      real(wp) :: p_m2(this%num_fields)
      real(wp) :: p_m1(this%num_fields)
      real(wp) :: p(this%num_fields)
      real(wp) :: p_p1(this%num_fields)
      real(wp) :: p_p2(this%num_fields)

      h = maxval(this%x_bounce) * 1.0e-4_wp

      p_m2 = this%phi_of_x(x)
      p_m1 = this%phi_of_x(x + 1.0e0_wp * h)
      p = this%phi_of_x(x + 2.0e0_wp * h)
      p_p1 = this%phi_of_x(x + 3.0e0_wp * h)
      p_p2 = this%phi_of_x(x + 4.0e0_wp * h)
      d2p = (35.0e0_wp * p_m2 - 104.0e0_wp * p_m1 +  &
        114.0e0_wp * p - 56.0e0_wp * p_p1 +  &
        11.0e0_wp * p_p2) / (12.0e0_wp * h ** 2)

    end function d2p_dx2

  end subroutine calc_normal_forces

  function grad_V(this, phi) result(dV)

    class(solver), intent(inout) :: this
    real(wp), intent(in) :: phi(this%num_fields)
    real(wp) :: dV(this%num_fields)

    real(wp) :: V
    real(wp) :: eps

    V = this%V(phi)
    eps = V * eps_gradient

    dV = gradient(this%V, phi, eps=eps)

  end function grad_V

end module bouncesolver__pathdeformation
