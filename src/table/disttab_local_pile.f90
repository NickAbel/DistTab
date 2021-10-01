!> This module contains the local pile object, for storing recently
!! obtained blocks from other subtables in a FIFO queue.
module disttab_local_pile
  use :: mpi
  use :: kind_params

  implicit none

  private
  public :: local_pile

  type :: local_pile

    ! The size of a block
    integer(i4), allocatable, dimension(:) :: block_dims
    ! An array (size equal to all blocks) of integers
    ! that is <= -1 if the block corresponding to the entry
    ! is not local, else gives the location of the block in
    ! the local pile FIFO queue
    integer(i4), allocatable, dimension(:) :: block_locator
    ! Number of state variables
    integer(i4) :: nvar
    ! Total number of blocks on the entire table
    integer(i4) :: total_blocks

  contains
    procedure, public, pass(this) :: push

    final :: local_pile_destructor

  end type local_pile

  interface local_pile
    module procedure :: local_pile_constructor
  end interface local_pile

contains

!> Constructor for the local pile object.
  type(local_pile) function local_pile_constructor() result(this)

  end function local_pile_constructor

!> destructor for the local pile type
  subroutine local_pile_destructor(this)
    type(local_pile) :: this

  end subroutine local_pile_destructor

!> Push a block onto the local pile
  subroutine push(this)
    class(local_pile), intent(in) :: this

  end subroutine push

end module disttab_local_pile
