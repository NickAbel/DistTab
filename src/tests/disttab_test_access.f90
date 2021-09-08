module disttab_test_access
  use :: disttab_table
  use :: iso_fortran_env

  implicit none
  private
  public :: access_test

  type :: access_test
    private
    type(table) :: lookup

  contains

    procedure, pass(this) :: run_value_test
    procedure, pass(this) :: run_value_cloud_test
    procedure, pass(this) :: run_get_map_get_test
    procedure, pass(this) :: run_get_perf_test
    procedure, pass(this) :: run_locality_test

    procedure, pass(this), private :: get_value_test
    procedure, pass(this), private :: get_value_cloud_test
    procedure, pass(this), private :: get_map_get_test
    procedure, pass(this), private :: get_perf_test
    procedure, pass(this), private :: locality_test

    procedure, pass(this), private :: fill_table_ascending_integers
    procedure, pass(this), private :: fill_cvars_linspace

  end type access_test

  interface access_test
    module procedure :: access_test_constructor
  end interface access_test

contains

  !> The constructor for the access_test type.
  !! Initializes the table given table dimensions.
  !!
  !! @param table_dimensions specifies the size of the table
  !! @return this the access_test object created
  type(access_test) function access_test_constructor(table_dimensions) result(this)
    integer (kind = int64), dimension(:), intent(in) :: table_dimensions

    this%lookup = table(table_dimensions)
    call this%fill_table_ascending_integers()
    call this%fill_cvars_linspace()

  end function access_test_constructor

  !> Test the ability to get correct values from the table from index,
  !! local coordinate, and global coordinates.
  !! @param this access_test object
  subroutine get_value_test(this)
    class(access_test), intent(inout) :: this
    real (kind = real64) :: r
    integer (kind = int64) :: ind, N
    integer (kind = int64), dimension(size(this%lookup%part_dims)) :: coord, coord_p, coord_b, box_dims
    real (kind = real64), dimension(this%lookup%table_dim_svar) :: val_ind, val_local_coord, val_global_coord

    N = size(this%lookup%part_dims)

    call random_number(r)
    r = r * this%lookup%table_dims_flat
    ind = ceiling(r)

    print *, "get_value_test begin, ind = ", ind

    ! Via index
    val_ind = this%lookup%index_to_value(ind)

    ! Via local coord
    box_dims = this%lookup%table_dims_padded(1:N) / this%lookup%part_dims
    call this%lookup%index_to_local_coord(ind, this%lookup%part_dims, box_dims, coord_p, coord_b)
    val_local_coord = this%lookup%local_coord_to_value(coord_p, coord_b)

    ! Via global coord
    coord = this%lookup%local_coord_to_global_coord(coord_p, coord_b, box_dims)
    val_global_coord = this%lookup%global_coord_to_value(coord)

    if (any(val_ind .ne. val_local_coord) .or. any(val_local_coord .ne. val_global_coord)) then
      print *, "<FAIL> get_value_test not working correctly..."
    else
      print *, "get_value_test passed!"
    end if

  end subroutine get_value_test

  !> Test the ability to get correct value clouds from the table from index,
  !! local coordinate, and global coordinates.
  !! Example in 3D case:
  !! In general, when the value at coordinate (i,j,k) is desired, return
  !! {i,i+1}x{j,j+1}x{k,k+1} (returning err if i+1,j+1,k+1, etc > Z).
  !! @param this access_test object
  subroutine get_value_cloud_test(this)
    class(access_test), intent(inout) :: this
    real (kind = real64) :: real_val(size(this%lookup%part_dims)), r
    integer (kind = int64) :: ind, N
    integer (kind = int64), dimension(size(this%lookup%part_dims)) :: coord, coord_p, coord_b, box_dims
    real (kind = real64), dimension(this%lookup%table_dim_svar, 2**size(this%lookup%part_dims)) :: val_cloud_ind, &
                                                     & val_cloud_local_coord, val_cloud_global_coord, &
                                                     & val_cloud_real

    N = size(this%lookup%part_dims)

    call random_number(r)
    r = r * this%lookup%table_dims_flat
    ind = ceiling(r)

    box_dims = this%lookup%table_dims_padded(1:N) / this%lookup%part_dims

    ! Re-roll if the value cloud would spill off the table's end
    do while (any(this%lookup%index_to_global_coord(ind, this%lookup%part_dims, box_dims) &
      & .ge. this%lookup%table_dims(1:N)))
      print *, "Random cloud will fall off table! Re-rolling..."
      call random_number(r)
      r = r * this%lookup%table_dims_flat
      ind = ceiling(r)
    end do

    print *, "get_value_cloud_test begin, ind = ", ind

    ! Get value cloud from integerized index
    val_cloud_ind = this%lookup%index_to_value_cloud(ind)

    ! Get value cloud from local coordinate decomposition
    call this%lookup%index_to_local_coord(ind, this%lookup%part_dims, box_dims, coord_p, coord_b)
    val_cloud_local_coord = this%lookup%local_coord_to_value_cloud(coord_p, coord_b)

    ! Get value cloud from global coordinate index
    coord = this%lookup%local_coord_to_global_coord(coord_p, coord_b, box_dims)
    val_cloud_global_coord = this%lookup%global_coord_to_value_cloud(coord)

    ! Get value cloud from normalized control variable location
    real_val = 0.0
    val_cloud_real = this%lookup%real_to_value_cloud(real_val)

    ! Compare value clouds (except the real value which is not related)
    if (any(val_cloud_ind .ne. val_cloud_local_coord)   &
      & .or. any(val_cloud_local_coord .ne. val_cloud_global_coord)) then
      print *, "<FAIL> get_value_cloud_test not working correctly..."
      print *, "Index --> Val Cloud: "
      print *, val_cloud_ind
      print *, "Local Coord --> Val Cloud: "
      print *, val_cloud_local_coord
      print *, "Global Coord --> Val Cloud: "
      print *, val_cloud_global_coord
      print *, "Real --> Val Cloud (not supposed to match): "
      print *, val_cloud_real
    else
      print *, "get_value_cloud_test passed!"
    end if

  end subroutine get_value_cloud_test

  !> Get a table value from random global coordinates, stash the result,
  !! change the partitioning scheme, get a table value from the same
  !! global coordinates, and check that the two values remain the same.
  !! @param this access_test object
  !! @param partition_dims partition dimensions to remap to before re-obtaining values
  subroutine get_map_get_test(this, partition_dims)
    class(access_test), intent(inout) :: this
    real (kind = real64) :: real_coords_rand(size(this%lookup%part_dims)), r
    integer (kind = int64) :: ind, N
    integer (kind = int64), dimension(size(this%lookup%part_dims)) :: coord, box_dims
    integer (kind = int64), dimension(size(this%lookup%part_dims)), intent(in) :: partition_dims
    real (kind = real64), dimension(this%lookup%table_dim_svar) :: val_real, val_real_map
    real (kind = real64), dimension(this%lookup%table_dim_svar, 2**size(this%lookup%part_dims)) :: val_cloud_global_coord, &
                                                                                  & val_cloud_global_coord_map

    N = size(this%lookup%part_dims)

    call random_number(r)
    r = r * this%lookup%table_dims_flat
    ind = ceiling(r)

    box_dims = this%lookup%table_dims_padded(1:N) / this%lookup%part_dims

    ! Re-roll if the value cloud would spill off the table's end
    do while (any(this%lookup%index_to_global_coord(ind, this%lookup%part_dims, box_dims) &
      & .ge. this%lookup%table_dims(1:N)))
      print *, "Random cloud will fall off table! Re-rolling..."
      call random_number(r)
      r = r * this%lookup%table_dims_flat
      ind = ceiling(r)
    end do

    print *, "get_map_get_test start", this%lookup%part_dims

    ! Get value cloud from global coordinate index
    coord = this%lookup%index_to_global_coord(ind, this%lookup%part_dims, box_dims)
    val_cloud_global_coord = this%lookup%global_coord_to_value_cloud(coord)

    ! Get value cloud from normalized control variable location
    call random_number(real_coords_rand)
    val_real = this%lookup%real_to_value(real_coords_rand)

    ! Remap
    call this%lookup%partition_remap(partition_dims, this%lookup%table_dims)

    ! Get value cloud from global coordinate index
    val_cloud_global_coord_map = this%lookup%global_coord_to_value_cloud(coord)
    val_real_map = this%lookup%real_to_value(real_coords_rand)
    if (any(val_cloud_global_coord .ne. val_cloud_global_coord_map)) then 
      print *, "<FAIL> get_map_get_test not working correctly...", this%lookup%part_dims
    else
      print *, "get_map_get_test passed!", this%lookup%part_dims
    end if

  end subroutine get_map_get_test

  !> Get a table value from random global coordinates, stash the result,
  !! change the partitioning scheme, get a table value from the same
  !! global coordinates, and check that the two values remain the same.
  !! @param this access_test object
  !! @param runs the number of times to fetch a random real coordinate
  subroutine get_perf_test(this, runs, M)
    class(access_test), intent(inout) :: this
    real (kind = real64) :: real_coords_rand(size(this%lookup%part_dims))
    real (kind = real64) :: t1, t2, time_total, time_total_opt
    integer (kind = int64) :: i, runs
    integer (kind = int64), dimension(size(this%lookup%part_dims)) :: coord, coord_opt
    integer (kind = int64), dimension(size(this%lookup%part_dims)), intent(in) :: M
    integer (kind = int64), dimension(:), allocatable :: buckets

    !print *, "starting get_perf_test"

    time_total = 0.0
    time_total_opt = 0.0

    allocate (buckets(sum(M)))
    buckets = this%lookup%real_to_global_coord_opt_preprocessor(M)

    do i = 1, runs
      call random_number(real_coords_rand)

      call cpu_time(t1)
      coord = this%lookup%real_to_global_coord(real_coords_rand)
      call cpu_time(t2)

      time_total = time_total + (t2 - t1)

      call cpu_time(t1)
      coord_opt = this%lookup%real_to_global_coord_opt(real_coords_rand, M, buckets)
      call cpu_time(t2)

      time_total_opt = time_total_opt + (t2 - t1)

      if (any(coord .ne. coord_opt)) then
        print *, "<FAIL> optimized real_to_global_coord_opt"
        print *, "requested value: ", real_coords_rand
        print *, "coord: ", coord
        print *, "coord_opt: ", coord_opt
        print *, "coord vars", this%lookup%ctrl_vars(coord(1)), this%lookup%ctrl_vars(coord(2) +&
          &this%lookup%table_dims(1)), "next: ",&
          &this%lookup%ctrl_vars(coord(1) + 1), &
          this%lookup%ctrl_vars(coord(2) + this%lookup%table_dims(1) + 1)
        print *, this%lookup%ctrl_vars(coord_opt(1)), this%lookup%ctrl_vars(coord_opt(2) + this%lookup%table_dims(1)), "next: ", &
          this%lookup%ctrl_vars(coord_opt(1) + 1), &
          this%lookup%ctrl_vars(coord_opt(2) + this%lookup%table_dims(1) + 1)
        print *, "table dims: ", this%lookup%table_dims
      end if

    end do

    deallocate (buckets)

    !write (*, *) "------------ TOTALS -----------"
    write (*, *) "time_total = ", time_total
    write (*, *) "time_total_opt = ", time_total_opt
  end subroutine get_perf_test

  subroutine locality_test(this, runs)
    class(access_test), intent(inout) :: this
    integer (kind = int64), intent(in) :: runs
    integer (kind = int64), dimension(size(this%lookup%table_dims)) :: coord
    integer (kind = int64), dimension(size(this%lookup%part_dims)) :: box_dims, coord_p, coord_b, part_dims
    real (kind = real64), dimension(size(this%lookup%table_dims)) :: coord_r
    real (kind = real64), dimension(this%lookup%table_dim_svar) :: s
    real (kind = real64) :: t1, t2, time_total
    integer (kind = int64) :: i, table_size, ind, N, k, ind_p, ind_b, div_p, div_b, box_size
    integer (kind = int64), dimension(:), allocatable :: accesses
    N = size(this%lookup%part_dims)

    allocate (accesses(runs))
    table_size = product(this%lookup%table_dims(1:N))
    print *, "locality_test"
    time_total = 0.0

    this%lookup%part_dims = 10
    box_dims = this%lookup%table_dims(1:N) / this%lookup%part_dims
    box_size = product(box_dims)

    print *, "partitioning: ", this%lookup%part_dims
    print *, "table size in memory: ", sizeof(this%lookup%elems) / 1000, "kilobytes"
    print *, "elem size in memory: ", sizeof(this%lookup%elems(:, 1))
    print *, "tile size in memory: ", 4.d0 * box_size / 1000.d0, "kilobytes"

    ind = 1
    print *, "Total accesses ", runs
    print *, "Table dims ", this%lookup%table_dims
    part_dims = this%lookup%part_dims

    do i = 1, runs

      ! Generate random global coordinate to grab
      ! Formula: coord_r(i) * range + offset
      ! for example, coord_r(i) * 10.0 + 20 gives a range of 10 and
      ! offset of 20, so coord_r(i) will be in [20, 30]
      call random_number(coord_r)
      coord_r(1) = coord_r(1) * 99.0 + 1
      coord_r(2) = coord_r(2) * 99.0 + 1
      coord_r(3) = coord_r(3) * 99.0 + 1
      coord = ceiling(coord_r)
      !coord(2) = mod(i * 100, this%lookup%table_dims(2)) + 1
      coord(4) = 1

      ! For test: Generate all the coordinates outside this timed loop
      ! store in search-coord array, loop through that array see how that goes
      do k = 1, N
        coord_p(k) = ceiling(dfloat(coord(k)) / box_dims(k))
        coord_b(k) = mod(coord(k), box_dims(k))
        if (coord_b(k) .eq. 0) coord_b(k) = box_dims(k)
      end do

      div_p = box_size
      ind_p = coord_p(N)
      do k = 1, N - 1
        div_p = div_p / part_dims(k)
        ind_p = ind_p + (coord_p(k) - 1) * div_p
      end do

      div_b = box_size
      ind_b = coord_b(N)
      do k = 1, N - 1
        div_b = div_b / box_dims(k)
        ind_b = ind_b + (coord_b(k) - 1) * div_b
      end do

      accesses(i) = (ind_p - 1) * box_size + ind_b

    end do

    call cpu_time(t1)
    do i = 1, runs

      ! Sanity check
      !if (any(coord .gt. this%lookup%table_dims) .or. any(coord .lt. 1)) then
      !  print *, "ERROR in locality_test: Generated coords off table!"
      !  print *, coord_r, coord
      !end if

      ! Time value retrieval
      !val = this%lookup%global_coord_to_value(coord)

      !ind = mod(ind + stride, table_size) + 1

      s = s + this%lookup%elems(1, accesses(i))

    end do
    call cpu_time(t2)
    time_total = time_total + t2 - t1

    print *, "Total time = ", time_total

    deallocate (accesses)

  end subroutine locality_test

  !> Fill lookup's table with ascending integers.
  !! @param this access_test object
  subroutine fill_table_ascending_integers(this)
    class(access_test), intent(inout) :: this
    integer (kind = int64) :: i

    do i = 1, this%lookup%table_dims_flat
      this%lookup%elems(:, i) = i
    end do

  end subroutine fill_table_ascending_integers

  !> Fill lookup's control variables with linspace values in [0, 1].
  !! @param this access_test object
  subroutine fill_cvars_linspace(this)
    class(access_test), intent(inout) :: this
    integer (kind = int64) :: i, N, j

    N = size(this%lookup%part_dims)

    do i = 1, N
      do j = 1, this%lookup%table_dims(i)
        this%lookup%ctrl_vars(j + sum(this%lookup%table_dims(:i - 1))) = &
                        & 1.d0 * (j - 1) / (this%lookup%table_dims(i) - 1)
      end do
    end do

  end subroutine fill_cvars_linspace

  !> Runs the get value test.
  !! @param this access_test object
  subroutine run_value_test(this)
    class(access_test), intent(inout) :: this

    call this%get_value_test()

  end subroutine run_value_test

  !> Runs the get value cloud test.
  !! @param this access_test object
  subroutine run_value_cloud_test(this)
    class(access_test), intent(inout) :: this

    call this%get_value_cloud_test()

  end subroutine run_value_cloud_test

  !> Runs the get value-remap-get value test.
  !! @param this access_test object
  !! @param partition_dims partition dimensions to remap to before re-obtaining values
  subroutine run_get_map_get_test(this, partition_dims)
    class(access_test), intent(inout) :: this
    integer (kind = int64), dimension(size(this%lookup%part_dims)) :: partition_dims

    call this%get_map_get_test(partition_dims)

  end subroutine run_get_map_get_test

  !> Runs the performance test for get functionality.
  !! @param this access_test object
  !! @param runs number of times to access a random real-valued coordinate
  subroutine run_get_perf_test(this, runs, M)
    class(access_test), intent(inout) :: this
    integer (kind = int64), intent(in) :: runs
    integer (kind = int64), intent(in), dimension(size(this%lookup%part_dims)) :: M

    call this%get_perf_test(runs, M)

  end subroutine run_get_perf_test

  !> Runs the performance test for get functionality.
  !! @param this access_test object
  !! @param runs number of times to access a random real-valued coordinate
  subroutine run_locality_test(this, runs)
    class(access_test), intent(inout) :: this
    integer (kind = int64), intent(in) :: runs

    call this%locality_test(runs)

  end subroutine run_locality_test

end module disttab_test_access
