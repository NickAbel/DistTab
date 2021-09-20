module communicator
  use :: mpi
  use :: kind_params

  implicit none

  private
  public :: comm

  type :: comm

    integer(i4) :: ierror, rank, nprocs, window
    integer(i4) :: int_size_mpi, dbl_size_mpi
    integer(kind=mpi_address_kind) :: target_displacement

  contains

    !procedure, public, pass(this) :: read_in
    final :: comm_destructor

  end type comm

  interface comm
    module procedure :: comm_constructor
  end interface comm

contains

!> Constructor for the MPI communicator object.
!!
!! @result this the created comm object
  type(comm) function comm_constructor() result(this)

    call mpi_init(this % ierror)
    call mpi_type_size(mpi_integer, this % int_size_mpi, this % ierror)
    call mpi_type_size(mpi_double, this % dbl_size_mpi, this % ierror)
    call mpi_comm_rank(mpi_comm_world, this % rank, this % ierror)
    call mpi_comm_size(mpi_comm_world, this % nprocs, this % ierror)

  end function comm_constructor

!> destructor for the comm type
!!
!! @param this the comm object to destruct
  subroutine comm_destructor(this)
    type(comm) :: this

    call mpi_finalize(this % ierror)

  end subroutine comm_destructor

end module communicator
