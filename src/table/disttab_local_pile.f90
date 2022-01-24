!> This module contains the local pile object, for storing recently
!! obtained blocks from other subtables in a FIFO queue.
module disttab_local_pile
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
    ! Size of the queue
    integer(i4) :: queue_size

    real(sp), allocatable, dimension(:, :) :: pile

  contains
    procedure, public, pass(this) :: push
    procedure, public, pass(this) :: resize

    final :: local_pile_destructor

  end type local_pile

  interface local_pile
    module procedure :: local_pile_constructor
  end interface local_pile

contains

!> Constructor for the local pile object.
  type(local_pile) function local_pile_constructor(queue_size, total_blocks, block_dims, nvar) result(this)
    integer(i4), dimension(:), intent(in) :: block_dims
    integer(i4), intent(in) :: total_blocks, nvar, queue_size

    !write (*, '(A,I0)') "Block size on table: ", total_blocks
    !write (*, '(A,I0)') "# of state vars: ", nvar
    !write (*, *) "Intra-block dimensions: ", block_dims

    allocate (this % block_dims(size(block_dims)))

    ! Starting with 1-entry blocks, todo
    this % block_dims = block_dims
    this % nvar = nvar
    this % queue_size = queue_size
    this % total_blocks = total_blocks

    allocate (this % pile(nvar, product(this % block_dims) * this % queue_size))

    allocate (this % block_locator(this % total_blocks))
    this % block_locator = -1

  end function local_pile_constructor

!> destructor for the local pile type
  subroutine local_pile_destructor(this)
    type(local_pile) :: this

  end subroutine local_pile_destructor

!> Push a block onto the local pile
  subroutine push(this, new_block, new_block_loc)
    class(local_pile) :: this
    real(sp), dimension(this % nvar, product(this % block_dims)), intent(in) :: new_block
    integer(i4) :: new_block_loc

    print *, "NVAR = ", this % nvar, " PILE SIZE = ", size(this % pile), " PRODUCT BLOCK DIMS = ", &
      & product(this % block_dims), " PILE UBOUND2 = ", ubound(this % pile, dim=2)

    ! Move the old blocks back
    this % pile(:, (product(this % block_dims) + 1):ubound(this % pile, dim=2)) = &
      & this % pile(:, 1:(ubound(this % pile, dim=2) - product(this % block_dims)))

    ! Place the new block at the front
    this % pile(:, 1:product(this % block_dims)) = new_block

    ! Indicate the removed block is gone, in the block locator
    if (maxval(this % block_locator) .ge. this % queue_size) then
      this % block_locator(maxloc(this % block_locator)) = -1
    end if

    ! Indicate the locations of the old blocks in the block locator
    this % block_locator = this % block_locator + sign(1, this % block_locator)

    ! Indicate the new block is at the front of the pile
    this % block_locator(new_block_loc) = 1

  end subroutine push

!> Resize the local pile
!! TODO an MPI fence is necessary when using this to prevent accessing deallocated memory, etc.
  subroutine resize(this, queue_size, total_blocks, block_dims, nvar)
    class(local_pile) :: this
    integer(i4), dimension(:), intent(in) :: block_dims
    integer(i4), intent(in) :: total_blocks, nvar, queue_size

    !write (*, '(A,I0)') "Block size on table: ", total_blocks
    !write (*, '(A,I0)') "# of state vars: ", nvar
    !write (*, *) "Intra-block dimensions: ", block_dims

    deallocate (this % block_dims)
    allocate (this % block_dims(size(block_dims)))

    ! Starting with 1-entry blocks, todo
    this % block_dims = block_dims
    this % nvar = nvar
    this % queue_size = queue_size
    this % total_blocks = total_blocks

    deallocate (this % pile)
    allocate (this % pile(nvar, product(this % block_dims) * this % queue_size))

    deallocate (this % block_locator)
    allocate (this % block_locator(this % total_blocks))
    this % block_locator = -1

  end subroutine resize


end module disttab_local_pile
