module bouncesolver__tp

  use bouncesolver__config, only : wp
  use bouncesolver__util, only : linspace
  use bouncesolver__util, only : pi
  use bouncesolver__util, only : is_equal
  use gradmin__descent, only : minimize
! use gradmin__newton, only : newton => minimize
! use gradmin__newton, only : set_eps_newton => set_eps
  use evortran__individuals_float, only : individual
  use evortran__evolutions_float, only : evolve_population
  use evortran__prng_rand, only : initialize_rands

  implicit none

  private

  real(wp) :: vev = 246.0e0_wp

  abstract interface
    function V_abstract(x) result(y)
      import wp
      real(wp), intent(in) :: x(:)
      real(wp) :: y
    end function V_abstract
  end interface

  type, public :: solver
    procedure(V_abstract), pointer,  &
      nopass :: V => null()
    integer :: d
    real(wp), allocatable :: phi_false(:)
    real(wp), allocatable :: phi_true(:)
    integer :: alpha = 3
    integer :: nx = 1000
    real(wp), allocatable :: x(:)
    real(wp), allocatable :: Vt(:)
    real(wp) :: S_min
    real(wp), allocatable :: x_min(:)
    real(wp), allocatable :: phi_min(:, :)
    real(wp), allocatable :: Vt_min(:)
    real(wp), allocatable :: V_min(:)
    real(wp), allocatable :: varphi_min(:)
  contains
    procedure, private :: allocate_arrays
    procedure, private :: calc_S
    procedure, private :: minimize_action
    procedure, private :: save_solution
  end type

  interface solver
    module procedure create_solver
  end interface solver

contains

  function create_solver(  &
    d, V, phi_false, phi_true,  &
    alpha, nx) result(this)

    integer, intent(in) :: d
    procedure(V_abstract) :: V
    real(wp), intent(in) :: phi_false(d)
    real(wp), intent(in) :: phi_true(d)
    integer, intent(in), optional :: alpha
    integer, intent(in), optional :: nx
    type(solver) :: this

    this%d = d
    this%V => V
    this%phi_false = phi_false
    this%phi_true = phi_true
    if (present(alpha)) this%alpha = alpha
    if (present(nx)) this%nx = nx

    call this%allocate_arrays()
    call this%minimize_action()
    call this%save_solution()

  end function create_solver

  subroutine allocate_arrays(this)

    class(solver), intent(inout) :: this

    integer :: nx
    integer :: d

    nx = this%nx
    d = this%d

    this%x = linspace(0.0e0_wp, 1.0e0_wp, nx)
    allocate(this%Vt(nx))
    allocate(this%x_min(nx * d))
    allocate(this%phi_min(nx, d))
    allocate(this%Vt_min(nx))
    allocate(this%V_min(nx))
    allocate(this%varphi_min(nx))

  end subroutine allocate_arrays

  function calc_S(  &
      this, p,  &
      Vt_out, V_out, varphi_out) result(y)

    class(solver), intent(in) :: this
    real(wp), intent(in) :: p(this%nx, this%d)
    real(wp), intent(out), optional :: Vt_out(this%nx)
    real(wp), intent(out), optional :: V_out(this%nx)
    real(wp), intent(out), optional :: varphi_out(this%nx)
    real(wp) :: y

    real(wp) :: S(this%nx - 1)
    real(wp) :: Vt(this%nx)
    real(wp) :: V(this%nx)
    real(wp) :: varphi(this%nx)
    real(wp) :: pc(this%d)
    real(wp) :: p0(this%d)
    real(wp) :: dp(this%nx - 1, this%d)
    integer :: i
    real(wp), parameter :: prefac =  &
      27.0e0_wp * pi ** 2 / 2.0e0_wp

    pc = this%phi_false
    p0 = p(this%nx, :)

    do i = 1, this%nx
      Vt(i) = this%V(pc) +  &
        this%x(i) ** 3 *  &
        (3.0e0_wp - 2.0e0_wp * this%x(i)) *  &
        (this%V(p0) - this%V(pc))
      if (i > 1) then
        if (Vt(i) > Vt(i - 1)) then
!         write(*,*) "Vt not monothonically decreasing."
          y = 1.0e10_wp
          return
        end if
      end if
    end do

    varphi(1) = 0.0e0_wp
    do i = 2, this%nx
      dp(i - 1, :) = p(i, :) - p(i - 1, :)
      varphi(i) = varphi(i - 1) +  &
        sqrt(dot_product(dp(i - 1, :), dp(i - 1, :)))
    end do

    do i = 1, this%nx
      V(i) = this%V(p(i, :))
    end do

    do i = 1, this%nx - 1
      ! Extra minus sign compared to Eq. (31)
      ! because otherwise S would be negative
      ! because Vt monothonically decreasing
      ! and negative
      S(i) = prefac *  &
        (V(i) + V(i + 1) - Vt(i) - Vt(i + 1)) ** 2 *  &
        dot_product(dp(i, :), dp(i, :)) ** 2 /  &
        (-(Vt(i + 1) - Vt(i))) ** 3
    end do

    y = sum(S)

    if (present(Vt_out)) Vt_out = Vt
    if (present(V_out)) V_out = V
    if (present(varphi_out)) varphi_out = varphi

  end function calc_S

  subroutine minimize_action(this)

    class(solver), intent(inout) :: this

    real(wp) :: phi_init(this%nx, this%d)
    real(wp) :: x_init(this%nx * this%d)
    real(wp) :: x_min(this%nx * this%d)
    real(wp) :: S_min
    real(wp) :: phi_min(this%nx, this%d)
    integer :: i
    integer :: j
!   integer :: newton_status
    type(individual) :: best_ind
    type(individual) :: init_ind

    do i = 1, this%d
      phi_init(:, i) = linspace(  &
        this%phi_false(i),  &
        this%phi_true(i),  &
        this%nx)
    end do

    do i = 1, this%nx
      do j = 0, this%d - 1
        x_init(i + j * this%nx) = phi_init(i, j + 1)
      end do
    end do

    call minimize(  &
      f, x_init, x_min, S_min,  &
      maxiter=100000, mode=1)

!   newton_status = 0
!   call set_eps_newton(1.0e-3_wp)
!   call newton(  &
!     f, x_init, x_min, S_min, maxiter=3, status=newton_status)
!   write(*,*) S_min, newton_status

!   call initialize_rands(mode="twister")
!   init_ind = individual(this%nx * this%d, fitness)
!   init_ind%genes = x_init / (1.2e0_wp * vev)
!   best_ind = evolve_population(  &
!     100, this%nx * this%d, fitness,  &
!     max_generations=100000,  &
!     selection="rank",  &
!     selection_size=50,  &
!     mating="two-point",  &
!     sbx_eta_c=3.0e0_wp,  &
!     elite_size=10,  &
!     verbose=.true.,  &
!     add_ind=init_ind,  &
!     mutate_prob=0.6e0_wp,  &
!     mutate_gene_prob=0.1e0_wp,  &
!     mutate_gaussian_sigma=0.1e0_wp,  &
!     mutate="gaussian")
!   x_min = 1.2e0_wp * vev * best_ind%genes

    this%x_min = x_min
    phi_min = reshape(x_min, [this%nx, this%d])
    this%phi_min = phi_min
    this%S_min = this%calc_S(phi_min)

  contains

    subroutine fitness(ind, fit)

      class(individual), intent(in) :: ind
      real(wp), intent(out) :: fit

      real(wp) :: x(size(ind%genes))
      real(wp) :: regu
      integer :: i
      integer :: j
      real(wp) :: phi(this%nx, this%d)

      x = 1.2e0_wp * vev * ind%genes

      phi = reshape(x, [this%nx, this%d])
      regu = 0.0e0_wp
      do i = 1, this%nx - 1
        do j  = 1, this%d
          regu = regu + (phi(i, j) - phi(i + 1, j)) ** 2
        end do
      end do

      fit = f(x) + regu

    end subroutine fitness

    function f(x) result(y)

      real(wp), intent(in) :: x(:)
      real(wp) :: y

      real(wp) :: phi(this%nx, this%d)
      integer :: i

      phi = reshape(x, [this%nx, this%d])

      y = this%calc_S(phi)

      ! boundary condition
       y = y + 1.0e2_wp * sum(phi(1, :) - this%phi_false) ** 2

      ! avoid trivial solution
      ! y = y + 1.0e0_wp / (sum(abs(x)) ** 2 + 1.0e-10_wp) + 1.0e0_wp


    end function f

  end subroutine minimize_action

  subroutine save_solution(this)

    class(solver), intent(inout) :: this

    real(wp) :: S_min
    real(wp) :: Vt_out(this%nx)
    real(wp) :: V_out(this%nx)
    real(wp) :: varphi_out(this%nx)

    S_min = this%calc_S(this%phi_min ,  &
      Vt_out, V_out, varphi_out)
    write(*,*) S_min

    if (.not. is_equal(S_min, this%S_min)) then
      write(*,*) "Kaputt.", S_min, this%S_min
      call exit
    end if

!   this%Vt_min = Vt_out
!   this%V_min = V_out
!   this%varphi_min = varphi_out

  end subroutine save_solution

end module bouncesolver__tp
