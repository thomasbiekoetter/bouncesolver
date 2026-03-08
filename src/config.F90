module bouncesolver__config

  implicit none

  private

  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: qp = selected_real_kind(30,4931)

#ifdef QUAD
  integer, parameter, public :: wp = qp
#else
  integer, parameter, public :: wp = dp
#endif

  character(len=*), parameter, public :: fmt = "(G20.10)"
  character(len=*), parameter, public :: fmt2 = "(G20.10,G20.10)"
  character(len=*), parameter, public :: fmt3 = "(G20.10,G20.10,G20.10)"

end module bouncesolver__config
