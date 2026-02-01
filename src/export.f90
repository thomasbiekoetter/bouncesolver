module bouncesolver__export

  use bouncesolver__config, only : wp
  use bouncesolver__config, only : fmt
  use bouncesolver__pdm, only : solver
  use csv_module

  implicit none

  private

  type(csv_file) :: f
  logical :: status_ok

  public :: csv_x_of_rho
  public :: csv_pot_of_rho
  public :: csv_pot_of_x
  public :: csv_phi_of_rho
  public :: csv_phi_of_x
  public :: csv_forces_of_x

contains

  subroutine csv_x_of_rho(pd, filename)

    type(solver), intent(in) :: pd
    character(len=*), intent(in) :: filename

    integer :: i

    call f%initialize(verbose=.true.)

    call f%open(  &
      trim(filename),  &
      n_cols=3,  &
      status_ok=status_ok)

    call f%add(["rho ", "x   ", "xdot"])
    call f%next_row()

    do i = 1, pd%n
      call f%add(  &
        [pd%rho(i), pd%xb(i), pd%xbdot(i)],  &
        real_fmt=fmt)
      call f%next_row()
    end do

    call f%close(status_ok)

  end subroutine csv_x_of_rho

  subroutine csv_pot_of_rho(pd, filename)

    type(solver), intent(in) :: pd
    character(len=*), intent(in) :: filename

    integer :: i

    call f%initialize(verbose=.true.)

    call f%open(  &
      trim(filename),  &
      n_cols=2,  &
      status_ok=status_ok)

    call f%add(["rho", "pot"])
    call f%next_row()

    do i = 1, pd%ib_max
      call f%add(  &
        [pd%rho(i), pd%pot(i)],  &
        real_fmt=fmt)
      call f%next_row()
    end do

    call f%close(status_ok)

  end subroutine csv_pot_of_rho

  subroutine csv_pot_of_x(pd, filename)

    type(solver), intent(in) :: pd
    character(len=*), intent(in) :: filename

    integer :: i

    call f%initialize(verbose=.true.)

    call f%open(  &
      trim(filename),  &
      n_cols=2,  &
      status_ok=status_ok)

    call f%add(["x", "p"])
    call f%next_row()

    do i = 1, pd%ib_max
      call f%add( &
        [pd%x(i), pd%pot(i)],  &
        real_fmt=fmt)
      call f%next_row()
    end do

    call f%close(status_ok)

  end subroutine csv_pot_of_x

  subroutine csv_phi_of_rho(pd, filename)

    type(solver), intent(in) :: pd
    character(len=*), intent(in) :: filename

    integer :: i

    call f%initialize(verbose=.true.)

    call f%open(  &
      trim(filename),  &
      n_cols=3,  &
      status_ok=status_ok)

    call f%add(["rho ", "phi1", "phi2"])
    call f%next_row()

    do i = 1, pd%n
      call f%add(  &
        [pd%rho(i), pd%phi(i, 1), pd%phi(i, 2)],  &
        real_fmt=fmt)
      call f%next_row()
    end do

    call f%close(status_ok)

  end subroutine csv_phi_of_rho

  subroutine csv_phi_of_x(pd, filename)

    type(solver), intent(in) :: pd
    character(len=*), intent(in) :: filename

    integer :: i

    call f%initialize(verbose=.true.)

    call f%open(  &
      trim(filename),  &
      n_cols=3,  &
      status_ok=status_ok)

    call f%add(["x   ", "phi1", "phi2"])
    call f%next_row()

    do i = 1, pd%ib_max
      call f%add(  &
        [pd%x(i), pd%phi(i, 1), pd%phi(i, 2)],  &
        real_fmt=fmt)
      call f%next_row()
    end do

    call f%close(status_ok)

  end subroutine csv_phi_of_x

  subroutine csv_forces_of_x(pd, filename)

    type(solver), intent(in) :: pd
    character(len=*), intent(in) :: filename

    integer :: i

    call f%initialize(verbose=.true.)

    call f%open(  &
      trim(filename),  &
      n_cols=5,  &
      status_ok=status_ok)

    call f%add(["x ", "N1", "N2", "P1", "P2"])
    call f%next_row()

    do i = 1, pd%ib_max
      call f%add([  &
        pd%xforce(i),  &
        pd%Nforce(i, 1), pd%Nforce(i, 2),  &
        pd%Pforce(i, 1), pd%Pforce(i, 2)], &
        real_fmt=fmt)
      call f%next_row()
    end do

    call f%close(status_ok)

  end subroutine csv_forces_of_x

end module bouncesolver__export
