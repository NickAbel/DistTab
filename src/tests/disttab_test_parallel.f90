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

    integer(i4), allocatable, dimension(:) :: part_dims
    integer(i4), allocatable, dimension(:) :: subtable_dims
    integer(i4), allocatable, dimension(:) :: table_dims

  contains

    procedure, pass(this) :: run_parallel_get_test

    procedure, pass(this), private :: parallel_get_test

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
!! @param partition_dimensions specifies the size of the partition blocks in each direction
!! @return this the parallel_test object which parallel_test_constructor instantiates
  type(parallel_test) function parallel_test_constructor(table_dimensions, subtable_dimensions, &
     & partition_dimensions) result(this)
    integer(i4), dimension(:), intent(in) :: table_dimensions
    integer(i4), dimension(:), intent(in) :: subtable_dimensions
    integer(i4), dimension(:), intent(in) :: partition_dimensions

    allocate (this % table_dims(size(table_dimensions)))
    allocate (this % subtable_dims(size(table_dimensions) - 1))
    allocate (this % part_dims(size(table_dimensions) - 1))

    this % table_dims = table_dimensions
    this % subtable_dims = subtable_dimensions
    this % part_dims = partition_dimensions
    this % lookup = table(this % table_dims, this % subtable_dims)

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
    integer(i4) :: rank, nprocs, real_size, i, ierror, ind
    real(sp) :: r

    call mpi_comm_size(mpi_comm_world, nprocs, ierror)
    call mpi_type_size(mpi_real, real_size, ierror)
    call mpi_comm_rank(mpi_comm_world, rank, ierror)

    if (rank .eq. 0) write (*, *) "table size ", this % table_dims, "with ", nprocs, " ranks"

    write (*, *) "on rank ", rank, " we have a table of dimensions ", &
      & this % table_dims(1:size(this % table_dims) - 1)

    if (mod(product(this % table_dims(1:size(this % table_dims) - 1)), nprocs) .ne. 0) then
      write (*, '(A,I0,A,I0,A)') "however the table size ", &
      & product(this % table_dims(1:size(this % table_dims) - 1)), " won't divide evenly over ", nprocs, " ranks"
    end if

    if (rank .eq. 0) print *, "sub-table dimensions supplied (", this % subtable_dims, &
      & ") give total ", product(ceiling(1.0 * this % table_dims(1:size(this % table_dims) - 1) &
      & / this % subtable_dims)), " sub-tables on table size ", &
      & this % table_dims(1:size(this % table_dims) - 1)

    print *, "rank ", rank, " subtable bounds: ", lbound(this % lookup % elems, dim=2), ubound(this % lookup % elems, dim=2)

    do i = lbound(this % lookup % elems, dim=2), ubound(this % lookup % elems, dim=2)
      this % lookup % elems(:, i) = i * 1.0
    end do

    call random_number(r)
    r = r * product(this % table_dims)
    ind = ceiling(r)
    write (*, '(A, I0, A, I0)') "Rank ", rank, ": Request index ", ind
    if (ind .ge. lbound(this % lookup % elems, dim=2) .and. ind .le. ubound(this % lookup % elems, dim=2)) then
      write (*, *) "Obtaining index ", ind, " which is local on rank ", rank, " --> ", this % lookup % index_to_value(ind)
    else if (ind .lt. 1 .or. ind .gt. product(this % table_dims)) then
      write (*, '(A, I0, A)') "ERROR in index_to_value call: Requested index (", ind, ") is not within lookup table bounds."
    else
      write (*, '(A, I0, A, I0, A, F8.2)') "Rank ", rank, ": Requested index (", ind, &
        & ") is on another sub-table. --> ", this % lookup % index_to_value(ind)
    end if

  end subroutine parallel_get_test

!> Runs the partition mapping test.
!! @param this the parallel_test object to which run_test belongs
  subroutine run_parallel_get_test(this)
    class(parallel_test), intent(inout) :: this

    call this % parallel_get_test()

  end subroutine run_parallel_get_test

end module disttab_test_parallel
