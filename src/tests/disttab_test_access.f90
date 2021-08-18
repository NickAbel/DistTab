module disttab_test_access
  use :: disttab_table

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

    procedure, pass(this), private :: get_value_test
    procedure, pass(this), private :: get_value_cloud_test
    procedure, pass(this), private :: get_map_get_test

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
    integer, dimension(:), intent(in) :: table_dimensions

    this%lookup = table(table_dimensions)
    call this%fill_table_ascending_integers()

  end function access_test_constructor

  !> Test the ability to get correct values from the table from index,
  !! local coordinate, and global coordinates.
  !! @param this access_test object
  subroutine get_value_test(this)
    class(access_test), intent(inout) :: this
    real :: r
    integer :: ind, N
    integer, dimension(size(this%lookup%part_dims)) :: coord, coord_p, coord_b, box_dims
    real, dimension(this%lookup%table_dim_svar) :: val_ind, val_local_coord, val_global_coord

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
      print *, "get_value_test not working correctly..."
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
    real :: real_val(size(this%lookup%part_dims)), r
    integer :: ind, N
    integer, dimension(size(this%lookup%part_dims)) :: coord, coord_p, coord_b, box_dims
    real, dimension(this%lookup%table_dim_svar, 2**size(this%lookup%part_dims)) :: val_cloud_ind, &
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
    call this%fill_cvars_linspace()
    real_val = 0.0
    val_cloud_real = this%lookup%real_to_value_cloud(real_val)

    ! Compare value clouds (except the real value which is not related)
    if (any(val_cloud_ind .ne. val_cloud_local_coord)   &
      & .or. any(val_cloud_local_coord .ne. val_cloud_global_coord)) then
      print *, "get_value_cloud_test not working correctly..."
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
    real :: real_coords_rand(size(this%lookup%part_dims)), r
    integer :: ind, N
    integer, dimension(size(this%lookup%part_dims)) :: coord, coord_p, coord_b, box_dims
    integer, dimension(size(this%lookup%part_dims)), intent(in) :: partition_dims
    real, dimension(this%lookup%table_dim_svar) :: val_real, val_real_map
    real, dimension(this%lookup%table_dim_svar, 2**size(this%lookup%part_dims)) :: val_cloud_global_coord, &
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
    call this%fill_cvars_linspace()
    call random_number(real_coords_rand)
    val_real = this%lookup%real_to_value(real_coords_rand)

    ! Remap
    call this%lookup%partition_remap(partition_dims, this%lookup%table_dims)

    ! Get value cloud from global coordinate index
    val_cloud_global_coord_map = this%lookup%global_coord_to_value_cloud(coord)
    val_real_map = this%lookup%real_to_value(real_coords_rand)
    if (any(val_cloud_global_coord .ne. val_cloud_global_coord_map) .or. &
        & any(val_real .ne. val_real_map)) then
      print *, "<FAIL> get_map_get_test", this%lookup%part_dims
    else
      print *, "get_map_get_test passed!", this%lookup%part_dims
    end if

  end subroutine get_map_get_test

  !> Fill lookup's table with ascending integers.
  !! @param this access_test object
  subroutine fill_table_ascending_integers(this)
    class(access_test), intent(inout) :: this
    integer :: i

    do i = 1, this%lookup%table_dims_flat
      this%lookup%elems(:, i) = i
    end do

  end subroutine fill_table_ascending_integers

  !> Fill lookup's control variables with linspace values in [0, 1].
  !! @param this access_test object
  subroutine fill_cvars_linspace(this)
    class(access_test), intent(inout) :: this
    integer :: i, N, j

    N = size(this%lookup%part_dims)

    do i = 1, N
      do j = 1, this%lookup%table_dims(i)
        this%lookup%ctrl_vars(j + sum(this%lookup%table_dims(:i - 1))) = 1.0 * (j - 1) / (this%lookup%table_dims(i) - 1)
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
    integer, dimension(size(this%lookup%part_dims)) :: partition_dims

    call this%get_map_get_test(partition_dims)

  end subroutine run_get_map_get_test

end module disttab_test_access
