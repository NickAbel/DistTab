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
    procedure, pass(this) :: run_parallel_partition_map_test
    procedure, pass(this) :: run_parallel_partition_map_unmap_test

    procedure, pass(this), private :: create_test_tables_padded_parallel
    procedure, pass(this), private :: create_test_tables_unpadded_parallel

    procedure, pass(this), private :: parallel_get_test
    procedure, pass(this), private :: local_pile_test
    procedure, pass(this), private :: parallel_partition_map_test
    procedure, pass(this), private :: parallel_partition_map_unmap_test

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
    integer(i4), dimension(:), allocatable :: ones

    allocate (this % table_dims(size(table_dimensions)))
    allocate (this % subtable_dims(size(table_dimensions) - 1))
    allocate (this % tile_dims(size(table_dimensions) - 1))
    allocate (ones(size(table_dimensions) - 1))

    ones = 1

    this % table_dims = table_dimensions
    this % subtable_dims = subtable_dimensions
    this % tile_dims = tile_dimensions
    this % lookup = table(this % table_dims, this % subtable_dims, mpi_comm_world)

    !call this % lookup % partition_remap_subtable(tile_dimensions, ones)

    ! Impose tile dimensions
    this % lookup % part_dims = this % tile_dims
    !print *, "total blocks: ", product(this % table_dims(1:size(this % table_dims) - 1)) / product(this % tile_dims)
    call this % lookup % pile % resize(10, product(this % table_dims(1:size(this % table_dims) - 1)) / &
      & product(this % tile_dims), this % tile_dims, this % lookup % nvar)

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
    call mpi_comm_size(this % lookup % communicator, nprocs, ierror)
    call mpi_type_size(mpi_real, real_size, ierror)
    call mpi_comm_rank(this % lookup % communicator, rank, ierror)

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
    integer(i4) :: rank, nprocs, real_size, ierror, ind, i, j, k, blk, displacement_cv
    real(sp), dimension(:, :), allocatable :: gold_tile, push_tile
    real(sp), dimension(this % lookup % nvar) :: table_value, pile_value
    real(sp) :: r

    ! MPI variables we'll need
    call mpi_comm_size(this % lookup % communicator, nprocs, ierror)
    call mpi_type_size(mpi_real, real_size, ierror)
    call mpi_comm_rank(this % lookup % communicator, rank, ierror)

    ! print subtable bounds on each rank
    print *, "rank ", rank, " subtable bounds: ", lbound(this % lookup % elems, dim=2), &
      & ubound(this % lookup % elems, dim=2), " nvar = ", this % lookup % nvar

    ! Impose tile dimensions
    this % lookup % pile % block_dims = this % tile_dims

    ! Impose table entries
    do i = lbound(this % lookup % elems, dim=2), ubound(this % lookup % elems, dim=2)
      do j = lbound(this % lookup % elems, dim=1), ubound(this % lookup % elems, dim=1)
        this % lookup % elems(j, i) = i + 0.01 * j
      end do
    end do

    ! Zero out the pile (for checking duplicate entries)
    this % lookup % pile % pile = 0.0

    do i = 1, 1000

      call mpi_win_fence(0, this % lookup % window, ierror)

      ! Check for duplicate entries in the local piles
      !do j = lbound(this % lookup % pile % pile, dim=2), ubound(this % lookup % pile % pile, dim=2)
      !  do k = lbound(this % lookup % pile % pile, dim=2), ubound(this % lookup % pile % pile, dim=2)
      !    if (j .ne. k .and. all(this % lookup % pile % pile(:, j) .eq. this % lookup % pile % pile(:, k)) &
      !        .and. any(this % lookup % pile % pile(:, j) .ne. 0)) then
      !      print *, "ERROR in local_pile_test: duplicated entries found in the pile on rank ", rank, &
      !        " at entries ", j, " and ", k, this % lookup % pile % pile(:, j), this % lookup % pile % pile(:, k)
      !    end if
      !  end do
      !end do

      ! Random target index to retrieve
      call random_number(r)
      r = r * product(this % table_dims(1:size(this % table_dims) - 1))
      ind = ceiling(r)

      !write (*, *) "Rank ", rank, ": Request index ", ind, lbound(this % lookup % elems, dim=2), &
      !  &  ubound(this % lookup % elems, dim=2)

      ! Check if index is on the table
      if (ind .ge. lbound(this % lookup % elems, dim=2) .and. ind .le. &
        & ubound(this % lookup % elems, dim=2)) then

        !write (*, *) "Obtaining index ", ind, " which is local on rank ", rank, &
        !  & " --> ", this % lookup % index_to_value(ind)

      else if (ind .lt. 1 .or. ind .gt. product(this % table_dims)) then

        !write (*, '(A, I0, A)') "ERROR in index_to_value call: Requested index (", ind, &
        !  & ") is not within lookup table bounds."

      else

        !write (*, '(A, I0, A, I0, A, F8.2)') "Rank ", rank, ": Requested index (", ind, &
        !  & ") is on another sub-table.", this % lookup % index_to_value(ind)

        blk = floor(1.0 * (ind - 1) / product(this % lookup % pile % block_dims)) + 1

        if (this % lookup % pile % block_locator(blk) .gt. 0) then

          displacement_cv = merge(mod(ind, product(this % lookup % part_dims)), &
          & product(this % lookup % part_dims), &
          & mod(ind, product(this % lookup % part_dims)) .ne. 0)

          table_value = this % lookup % index_to_value(ind)
          pile_value = this % lookup % pile % pile(:, product(this % lookup % pile % block_dims) * &
          & (this % lookup % pile % block_locator(blk) - 1) + displacement_cv)

          if (any(table_value .ne. pile_value)) print *, "FAIL in local_pile_test: table value and pile value not equal"

          !write (*, *) "Rank ", rank, ": Requested index ", ind, &
          !& " is on another sub-table ", this % lookup % index_to_value(ind), &
          !& " but is found on the local pile slot ", this % lookup % pile % block_locator(blk) &
          !&, ":", displacement_cv, this % lookup % pile % pile(:, &
          !& product(this % lookup % pile % block_dims) * (this % lookup % pile % block_locator(blk) - 1) + 1: &
          !& product(this % lookup % pile % block_dims) * (this % lookup % pile % block_locator(blk) - 1) + &
          !& product(this % lookup % pile % block_dims))

        else

          allocate (push_tile(this % lookup % pile % nvar, product(this % lookup % pile % block_dims)))

          do j = 1, product(this % lookup % pile % block_dims)
            k = floor(1.0 * (ind - 1) / product(this % lookup % pile % block_dims)) * &
            & product(this % lookup % pile % block_dims) + j
            push_tile(:, j) = this % lookup % index_to_value(k)
          end do

          call this % lookup % pile % push(push_tile, blk)

          deallocate (push_tile)

        end if
      end if
    end do

    print *, this % lookup % pile % block_locator
    print *, this % lookup % pile % pile

    !allocate (gold_tile(this % lookup % pile % nvar, product(this % lookup % pile % block_dims)))

    !do i = 1, this % lookup % pile % total_blocks
    !  do j = 1, product(this % lookup % pile % block_dims)
    !    do k = 1, this % lookup % pile % nvar
    !      gold_tile(k, j) = (i - 1) * product(this % lookup % pile % block_dims) + j + (k - 1) / 10.0
    !    end do
    !  end do
    !  if (this % lookup % pile % block_locator(i) .gt. 0) then
    !    print *, gold_tile
    !    !if (any(gold_tile .ne. this % lookup % pile % pile(:, &
    !    !& this % lookup % pile % block_locator(i)*product(this % lookup % pile % block_dims) : &
    !    !& (this % lookup % pile % block_locator(i))*product(this % lookup % pile % block_dims) + &
    !    !& product(this % lookup % pile % block_dims) - 1))) then
    !    !  print *, "!!!!!!FAIL in local pile test", &
    !    !& "GOLD: [", gold_tile, "] ", "PILE SLOT: [", this % lookup % pile % pile(:, &
    !    !& this % lookup % pile % block_locator(i-1)*product(this % lookup % pile % block_dims): &
    !    !& (this % lookup % pile % block_locator(i))*product(this % lookup % pile % block_dims) - 1), "] "
    !    !end if
    !  end if
    !end do

    !deallocate (gold_tile)

  end subroutine local_pile_test

  subroutine create_test_tables_padded_parallel(this)
    class(parallel_test), intent(inout) :: this

  end subroutine create_test_tables_padded_parallel

  subroutine create_test_tables_unpadded_parallel(this)
    class(parallel_test), intent(inout) :: this
    integer(i4) :: rank, nprocs, i, j, k, ierror, N, ind_b
    character :: a
    character(len=:), allocatable :: str, str_sorted
    integer(i4), dimension(size(this % lookup % part_dims)) :: coord, coord_reord, tile_dims, coord_p, coord_b, &
    & t_b, s_b, p_b, t_1d, t_flat, coord_gs, coord_re

    N = size(this % lookup % part_dims)

    ! MPI variables we'll need
    call mpi_comm_size(this % lookup % communicator, nprocs, ierror)
    call mpi_comm_rank(this % lookup % communicator, rank, ierror)

    ! Find total partitions in each dimension
    tile_dims = this % lookup % subtable_dims(1:N) / this % lookup % part_dims

    ! Edit string to reflect rank number
    str = 'partition_test_table.nopad.DistTab.tmp.dat.rank0'
    a = str(len(str):len(str))
    str(len(str):len(str)) = char(ichar(a) + rank)

    open (52, file=str, action='write')

    t_b = this % lookup % table_dims(1:N) / this % lookup % part_dims
    s_b = this % lookup % table_dims(1:N) / this % lookup % subtable_dims(1:N)
    p_b = this % lookup % subtable_dims(1:N) / this % lookup % part_dims
    t_1d = 1
    t_1d(1) = product((this % lookup % table_dims(1:N) / this % lookup % part_dims) / &
            & (this % lookup % subtable_dims(1:N) / this % lookup % part_dims))
    t_flat = p_b * t_1d

    do i = product(this % lookup % subtable_dims) * rank + 1, (rank + 1) * product(this % lookup % subtable_dims)
      call this % lookup % index_to_local_coord(i, this % lookup % part_dims, tile_dims, coord_p, coord_b)
      coord = this % lookup % local_coord_to_global_coord(coord_p, coord_b, tile_dims)

      ind_b = this % lookup % global_coord_to_index(coord_p, p_b, t_flat)
      coord_gs = this % lookup % index_to_global_coord(ind_b, s_b, p_b)
      coord_re = this % lookup % local_coord_to_global_coord(coord_gs, coord_b, tile_dims)

      write (52, *) (coord_re(j), j=1, N), (i, k=1, this % lookup % nvar)
    end do

    close (52)

    ! Edit string to reflect rank number
    str_sorted = 'partition_test_table_sorted.nopad.DistTab.tmp.dat.rank0'
    a = str_sorted(len(str_sorted):len(str_sorted))
    str_sorted(len(str_sorted):len(str_sorted)) = char(ichar(a) + rank)

    ! This gawk loop will sort the coordinate columns 1 to N in order, which
    ! emulates the storage format of Alya tables.
    call execute_command_line("gawk -F ',' '                                &
                              &   {                                         &
                              &      for(i=1;i<=NF;i++){sorter[i][NR]=$i}   &
                              &   }                                         &
                              &   END{                                      &
                              &      for(i=1;i<=NF;i++){asort(sorter[i])}   &
                              &      for(j=1;j<=NR;j++){                    &
                              &          for(i=1;i<NF;i++){                 &
                              &              printf ""%s,"",sorter[i][j]    &
                              &          }                                  &
                              &          print sorter[i][j]                 &
                              &      }                                      &
                              &  }' "//str//" > "//str_sorted)

    open (53, file=str_sorted, action='read')

    ! After sorting to the "Alya format," load the table back into the elements array.
    do i = product(this % lookup % subtable_dims) * rank + 1, (rank + 1) * product(this % lookup % subtable_dims)
      read (53, *) (coord(j), j=1, N), this % lookup % elems(:, i)
    end do

    close (53)

  end subroutine create_test_tables_unpadded_parallel

  subroutine parallel_partition_map_test(this)
    class(parallel_test), intent(inout) :: this
    integer(i4) :: rank, nprocs, real_size, i_old, i, j, k, ndim, n_subtable, ierror, head, tail, temp, ind_r
    integer(i4), dimension(size(this % tile_dims)) :: coord, coord_p, coord_b, coord_r, subtable_blks, ones, &
    & part_blks, table_dims_map
    double precision, allocatable, dimension(:, :) :: elems_gold_std
    character :: a
    character(len=:), allocatable :: str

    ndim = size(this % tile_dims)
    ones = 1

    ! MPI variables we'll need
    call mpi_comm_size(this % lookup % communicator, nprocs, ierror)
    call mpi_type_size(mpi_real, real_size, ierror)
    call mpi_comm_rank(this % lookup % communicator, rank, ierror)

    do i = product(this % subtable_dims) * rank + 1, (rank + 1) * product(this % subtable_dims)
      this % lookup % elems(:, i) = i
    end do

    call this % lookup % partition_remap_subtable(this % tile_dims, this % lookup % subtable_dims)

    ! Below to the end of the subroutine is the real parallel mapping test
    ! Create gold standard and fill table with Alya-format sorted version of table
    !call this % create_test_tables_unpadded_parallel()

    !call this % lookup % partition_remap_subtable(this % tile_dims, this % lookup % subtable_dims)

    ! ! Edit string to reflect rank number
    ! str = 'partition_test_table.nopad.DistTab.tmp.dat.rank0'
    ! a = str(len(str):len(str))
    ! str(len(str):len(str)) = char(ichar(a) + rank)

    ! ! Open, read, and close the file to the gold standard
    ! allocate (elems_gold_std(this % lookup % nvar, &
    ! & product(this % lookup % subtable_dims) * rank + 1:(rank + 1) * product(this % lookup % subtable_dims)))
    ! open (54, file=str, action='read')
    ! do i = product(this % lookup % subtable_dims) * rank + 1, (rank + 1) * product(this % lookup % subtable_dims)
    !   read (54, *) (coord(j), j=1, N), elems_gold_std(:, i)
    ! end do
    ! close (54)

    ! ! Check if the partition mapping of the elements is equal to the gold standard generated in create_test_tables
    ! if (all(this % lookup % elems .eq. elems_gold_std)) then
    !   write (*, *) "n = [", this % lookup % table_dims, "], q = [", &
    !     this % lookup % part_dims, "]: ", "parallel_partition_map_test passed! Tables will be deleted. Rank ", rank, "."
    !   if (rank .eq. 0) call execute_command_line('rm *.DistTab.tmp.dat.rank*')
    ! else if (all(abs(this % lookup % elems - elems_gold_std) .lt. 0.0005)) then
    !   write (*, *) "n = [", this % lookup % table_dims, "], q = [", &
    !     this % lookup % part_dims, "]: ", &
    !     & "trivially small diff in parallel_partition_map_test. most likely a pass. Tables not deleted. Rank ", rank, "."
    ! else
    !   write (*, *) "n = [", this % lookup % table_dims, "], q = [", &
    !     this % lookup % part_dims, "]: ", &
    !     & "parallel_partition_map_test conditional not working right. Tables not deleted. Rank ", rank, "."
    ! end if

    !deallocate (elems_gold_std)

  end subroutine parallel_partition_map_test

  subroutine parallel_partition_map_unmap_test(this)
    class(parallel_test), intent(inout) :: this
    integer(i4) :: rank, nprocs, real_size, i, j, k, ierror, ind, tlb_index, tub_index
    integer(i4), dimension(size(this % tile_dims)) :: global_coords
    integer(i4), dimension(size(this % tile_dims)) :: tile_lower_bound, tile_upper_bound
    real(sp) :: r
    real(sp), dimension(this % lookup % nvar, product(this % tile_dims)) :: tile_buffer

    ! MPI variables we'll need
    call mpi_comm_size(this % lookup % communicator, nprocs, ierror)
    call mpi_type_size(mpi_real, real_size, ierror)
    call mpi_comm_rank(this % lookup % communicator, rank, ierror)

    print *, "parallel partition mapping-unmapping test"

  end subroutine parallel_partition_map_unmap_test

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

!> Runs the parallel partition mapping test.
!! @param this the parallel_test object
  subroutine run_parallel_partition_map_test(this)
    class(parallel_test), intent(inout) :: this

    call this % parallel_partition_map_test()

  end subroutine run_parallel_partition_map_test

!> Runs the parallel partition mapping-unmapping test.
!! @param this the parallel_test object
  subroutine run_parallel_partition_map_unmap_test(this)
    class(parallel_test), intent(inout) :: this

    call this % parallel_partition_map_unmap_test()

  end subroutine run_parallel_partition_map_unmap_test

end module disttab_test_parallel
