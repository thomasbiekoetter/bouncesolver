module bouncesolver__tp

  use bouncesolver__config, only : wp
  use bouncesolver__util, only : linspace

  implicit none

  private

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
  contains
    procedure, private :: allocate_arrays
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

  end function create_solver

  subroutine allocate_arrays(this)

    class(solver), intent(inout) :: this

    integer :: nx

    nx = this%nx

    this%x = linspace(0.0e0_wp, 1.0e0_wp, nx)
    allocate(this%Vt(nx))

  end subroutine allocate_arrays

end module bouncesolver__tp
