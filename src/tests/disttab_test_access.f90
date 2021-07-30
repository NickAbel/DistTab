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

    procedure, pass(this), private :: get_value_test
    procedure, pass(this), private :: get_value_cloud_test
    procedure, pass(this), private :: fill_table_ascending_integers

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
    val_ind = this%lookup%index_to_value(ind)
    box_dims = this%lookup%table_dims_padded(1:N) / this%lookup%part_dims
    call this%lookup%index_to_local_coord(ind, this%lookup%part_dims, box_dims, coord_p, coord_b)
    val_local_coord = this%lookup%local_coord_to_value(coord_p, coord_b)
    coord = this%lookup%local_coord_to_global_coord(coord_p, coord_b, this%lookup%part_dims, box_dims)
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
  subroutine get_value_cloud_test(this)
    class(access_test), intent(inout) :: this
    real :: r
    integer :: ind, N
    integer, dimension(size(this%lookup%part_dims)) :: coord, coord_p, coord_b, box_dims
    real, dimension(this%lookup%table_dim_svar, 2**size(this%lookup%part_dims)) :: val_cloud_ind, &
                                                     & val_cloud_local_coord, val_cloud_global_coord

    N = size(this%lookup%part_dims)
    call random_number(r)
    !r = r*this%lookup%table_dims_flat
    !ind = ceiling(r)
    ind = 1
    print *, "get_value_cloud_test begin, ind = ", ind
    val_cloud_ind = this%lookup%index_to_value_cloud(ind)
    box_dims = this%lookup%table_dims_padded(1:N) / this%lookup%part_dims
    call this%lookup%index_to_local_coord(ind, this%lookup%part_dims, box_dims, coord_p, coord_b)
    val_cloud_local_coord = this%lookup%local_coord_to_value_cloud(coord_p, coord_b)
    coord = this%lookup%local_coord_to_global_coord(coord_p, coord_b, this%lookup%part_dims, box_dims)
    val_cloud_global_coord = this%lookup%global_coord_to_value_cloud(coord)
    if (any(val_cloud_ind .ne. val_cloud_local_coord) .or. any(val_cloud_local_coord .ne. val_cloud_global_coord)) then
      print *, "get_value_cloud_test not working correctly..."
    else
      print *, "get_value_cloud_test passed!"
    end if
    print *, val_cloud_ind
    print *, val_cloud_local_coord
    print *, val_cloud_global_coord

  end subroutine get_value_cloud_test

  !> Fill lookup's table with ascending integers.
  subroutine fill_table_ascending_integers(this)
    class(access_test), intent(inout) :: this
    integer :: i

    do i = 1, this%lookup%table_dims_flat
      this%lookup%elems(:, i) = i
    end do

  end subroutine fill_table_ascending_integers

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

end module disttab_test_access
