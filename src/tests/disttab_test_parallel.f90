module disttab_test_parallel
  use :: disttab_table
  use :: kind_params
  use :: mpi

  implicit none
  private
  public :: parallel_test

  type :: parallel_test

    private

    type(table) :: lookup

    integer(i4), allocatable, dimension(:) :: tile_dims
    integer(i4), allocatable, dimension(:) :: subtable_dims
    integer(i4), allocatable, dimension(:) :: table_dims

  contains

    procedure, pass(this) :: run_parallel_get_test
    procedure, pass(this) :: run_local_pile_test

    procedure, pass(this), private :: parallel_get_test
    procedure, pass(this), private :: local_pile_test

    final :: parallel_test_destructor

  end type parallel_test

  interface parallel_test
    module procedure :: parallel_test_constructor
  end interface parallel_test

contains

!> The constructor for the parallel_test type.
!! Initializes the table and partition dimensions and creates a table
!! object lookup for testing.
!!
!! @param table_dimensions specifies the size of the table in the two control variable and one state variable dimension
!! @param tile_dimensions specifies the size of the tiles in each direction
!! @return this the parallel_test object which parallel_test_constructor instantiates
  type(parallel_test) function parallel_test_constructor(table_dimensions, subtable_dimensions, &
     & tile_dimensions) result(this)
    integer(i4), dimension(:), intent(in) :: table_dimensions
    integer(i4), dimension(:), intent(in) :: subtable_dimensions
    integer(i4), dimension(:), intent(in) :: tile_dimensions

    allocate (this % table_dims(size(table_dimensions)))
    allocate (this % subtable_dims(size(table_dimensions) - 1))
    allocate (this % tile_dims(size(table_dimensions) - 1))

    this % table_dims = table_dimensions
    this % subtable_dims = subtable_dimensions
    this % tile_dims = tile_dimensions
    this % lookup = table(this % table_dims, this % subtable_dims)

    ! Impose tile dimensions
    this % lookup % part_dims = this % tile_dims

  end function parallel_test_constructor

!> destructor for the parallel test type
!!
!! @param this the table object to destruct
!! @todo what is this exactly doing? is deallocate_table call necessary?
  subroutine parallel_test_destructor(this)
    type(parallel_test) :: this

  end subroutine parallel_test_destructor

  subroutine parallel_get_test(this)
    class(parallel_test), intent(inout) :: this
    integer(i4) :: rank, nprocs, real_size, i, j, k, ierror, ind, tlb_index, tub_index
    integer(i4), dimension(size(this % tile_dims)) :: global_coords
    integer(i4), dimension(size(this % tile_dims)) :: tile_lower_bound, tile_upper_bound
    real(sp) :: r
    real(sp), dimension(this % lookup % nvar, product(this % tile_dims)) :: tile_buffer

    ! MPI variables we'll need
    call mpi_comm_size(mpi_comm_world, nprocs, ierror)
    call mpi_type_size(mpi_real, real_size, ierror)
    call mpi_comm_rank(mpi_comm_world, rank, ierror)

    do i = lbound(this % lookup % elems, dim=2), ubound(this % lookup % elems, dim=2)
      this % lookup % elems(:, i) = i
    end do

    call random_number(r)
    r = r * product(this % table_dims(1:size(this % tile_dims)))

    ind = ceiling(r)
    print *, "rank ", rank, " unreorganized table parallel_get_test requesting index ", ind

    if (any(abs(ind - this % lookup % index_to_value(ind)) .ne. 0)) then
      print *, "!!DIFF in parallel_get_test getting index ", ind, ", value ", this % lookup % index_to_value(ind), &
        & " should equal index"
    end if

  end subroutine parallel_get_test

  subroutine local_pile_test(this)
    class(parallel_test), intent(inout) :: this
    integer(i4) :: rank, nprocs, real_size, ierror, ind, i, j, k, blk
    real(sp), dimension(:,:), allocatable :: gold_tile, push_tile
    real(sp) :: r

    ! MPI variables we'll need
    call mpi_comm_size(mpi_comm_world, nprocs, ierror)
    call mpi_type_size(mpi_real, real_size, ierror)
    call mpi_comm_rank(mpi_comm_world, rank, ierror)

    print *, "local_pile_test, rank ", rank

    ! print subtable bounds on each rank
    print *, "rank ", rank, " subtable bounds: ", lbound(this % lookup % elems, dim=2), &
      & ubound(this % lookup % elems, dim=2)

    do i = lbound(this % lookup % elems, dim=2), ubound(this % lookup % elems, dim=2)
      this % lookup % elems(:, i) = i
    end do

    print *, "part dims: ", this % lookup % part_dims

    do i = 1, 100

      ! Check for duplicate entries in the local piles
      do j = 1, size(this % lookup % pile % pile) 
        do k = 1, size(this % lookup % pile % pile) 
          if (j .ne. k .and. all(this % lookup % pile % pile(:, j) .eq. this % lookup % pile % pile(:, k)) &
              .and. any(this % lookup % pile % pile(:,j) .ne. 0)) then 
            print *, "ERROR in local_pile_test: duplicated entries found in the pile on rank ", rank, &
                   " at entries ", j, " and ", k, this % lookup % pile % pile(:,j) 
          end if
        end do
      end do

      call random_number(r)
      r = r * product(this % table_dims)
      ind = ceiling(r)

      !write (*, '(A, I0, A, I0)') "Rank ", rank, ": Request index ", ind

      if (ind .ge. lbound(this % lookup % elems, dim=2) .and. ind .le. &
        & ubound(this % lookup % elems, dim=2)) then
      !  write (*, *) "Obtaining index ", ind, " which is local on rank ", rank, &
      !    & " --> ", this % lookup % index_to_value(ind)
      else if (ind .lt. 1 .or. ind .gt. product(this % table_dims)) then
        write (*, '(A, I0, A)') "ERROR in index_to_value call: Requested index (", ind, &
          & ") is not within lookup table bounds."
      else
      !  write (*, '(A, I0, A, I0, A, F8.2)') "Rank ", rank, ": Requested index (", ind, &
      !    & ") is on another sub-table.", this % lookup % index_to_value(ind)
      blk = floor(1.0 * (ind - 1) / product(this % lookup % pile % block_dims)) + 1
        if (this % lookup % pile % block_locator(blk) .gt. 0) then
          write (*, *) "Rank ", rank, ": Requested index (", ind, &
          & " is on another sub-table ", this % lookup % index_to_value(ind), &
          & " but is found on the local pile slot ", this % lookup % pile % block_locator(blk), &
          & this % lookup % pile % pile(:, &
          & product(this % lookup % pile % block_dims)*(this % lookup % pile % block_locator(blk) - 1) + 1 : &
          & product(this % lookup % pile % block_dims)*(this % lookup % pile % block_locator(blk) - 1) + &
          & product(this % lookup % pile % block_dims))
        else
          allocate (push_tile(this % lookup % pile % nvar, product(this % lookup % pile % block_dims)))
          do j = 1, product(this % lookup % pile % block_dims)
            k = floor(1.0 * (ind - 1) / product(this % lookup % pile % block_dims)) * &
            & product(this % lookup % pile % block_dims) + j
            push_tile(:,j) = this % lookup % index_to_value(k)
          end do
          call this % lookup % pile % push(push_tile, blk)
          deallocate (push_tile)
        end if
      end if
    end do

    if (rank .eq. 0) print *, this % lookup % pile % block_locator
    if (rank .eq. 0) print *, this % lookup % pile % pile

    allocate (gold_tile(this % lookup % pile % nvar, product(this % lookup % pile % block_dims)))

!    do i = 1, this % lookup % pile % total_blocks
!      gold_tile = (i-1)*product(this % lookup % pile % block_dims) + 1
!      if (this % lookup % pile % block_locator(i) .gt. 0) then
!        print *, gold_tile
!        if (any(gold_tile .ne. this % lookup % pile % pile(:, & 
!        & this % lookup % pile % block_locator(i)*product(this % lookup % pile % block_dims) : &
!        & (this % lookup % pile % block_locator(i))*product(this % lookup % pile % block_dims) + &
!        & product(this % lookup % pile % block_dims) - 1))) then
!          print *, "!!!!!!FAIL in local pile test", &
!        & gold_tile, this % lookup % pile % pile(:, & 
!        & this % lookup % pile % block_locator(i)*product(this % lookup % pile % block_dims): &
!        & (this % lookup % pile % block_locator(i))*product(this % lookup % pile % block_dims) + &
!        & product(this % lookup % pile % block_dims) - 1)
!        end if
!      end if
!    end do

    deallocate (gold_tile)

  end subroutine local_pile_test

!> Runs the partition mapping test.
!! @param this the parallel_test object
  subroutine run_parallel_get_test(this)
    class(parallel_test), intent(inout) :: this

    call this % parallel_get_test()

  end subroutine run_parallel_get_test

!> Runs the local pile test.
!! @param this the parallel_test object
  subroutine run_local_pile_test(this)
    class(parallel_test), intent(inout) :: this

    call this % local_pile_test()

  end subroutine run_local_pile_test

end module disttab_test_parallel
