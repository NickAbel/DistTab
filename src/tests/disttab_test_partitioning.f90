module disttab_test_partitioning
  use :: disttab_table

  implicit none
  private
  public :: partitioning_test

  type :: partitioning_test
    private
    type(table) :: lookup
    integer, allocatable, dimension(:) :: table_dims
    integer, allocatable, dimension(:) :: part_dims
  contains
    procedure, pass(this) :: run_map_test
    procedure, pass(this) :: run_map_unmap_test

    procedure, pass(this), private :: partition_map_test
    procedure, pass(this), private :: partition_map_unmap_test
    procedure, pass(this), private :: create_test_tables_padded
    procedure, pass(this), private :: create_test_tables_unpadded

  end type partitioning_test

  interface partitioning_test
    module procedure :: partitioning_test_constructor
  end interface partitioning_test

contains

  !> The constructor for the partitioning_test type.
  !! Initializes the table and partition dimensions and creates a table
  !! object lookup for testing.
  !!
  !! @param table_dimensions specifies the size of the table in the two control variable and one state variable dimension
  !! @param partition_dimensions specifies the size of the partition blocks in each direction
  !! @return this the partitioning_test object which partitioning_test_constructor instantiates
  type(partitioning_test) function partitioning_test_constructor(table_dimensions, partition_dimensions) result(this)
    integer, dimension(:), intent(in) :: table_dimensions
    integer, dimension(:), intent(in) :: partition_dimensions

    allocate (this%table_dims(size(table_dimensions)))
    allocate (this%part_dims(size(table_dimensions) - 1))

    this%table_dims = table_dimensions
    this%part_dims = partition_dimensions
    this%lookup = table(this%table_dims)

  end function partitioning_test_constructor

  !> A verification test for the table object's partition mapping algorithm.
  !! First, calls create_test_tables to fill the elements of this%lookup
  !! with pre-partitioned entries.
  !! Then, calls this%lookup%partition_remap to get the result of the
  !! partition mapping algorithm.
  !! Next creates the array elems_gold_std, which is a reference
  !! solution to the partition mapping problem with table size of
  !! this%lookup and the partition dimensions of this%part_dims.
  !!
  !! The entries of elems_gold_std are ascending integers.
  !!
  !! elems_gold_std is then compared to this%lookup%elems. If all
  !! entries are equal, a pass is reported. If there is no difference in
  !! absolute value greater than 0.0005 between all elements of the two tables,
  !! a different message is given (as this could be due to precision issues
  !! that I am not sure don't exist.) Else, the test is failed.
  !!
  !! @param this the partition_test object to which partition_map_test belongs
  !! @todo be sure there is no precision issue
  subroutine partition_map_test(this)
    class(partitioning_test), intent(inout) :: this
    double precision, allocatable, dimension(:, :) :: elems_gold_std
    integer :: i, j, N
    integer :: coord(size(this%part_dims))

    print *, "partition_map_test begin"

      !! Compute padded dimensions for test tables
    this%lookup%table_dims_padded = this%table_dims
    do i = lbound(this%lookup%table_dims_padded, dim=1), ubound(this%lookup%table_dims_padded, dim=1) - 1
      do while (mod(this%lookup%table_dims_padded(i), this%part_dims(i)) .ne. 0)
        this%lookup%table_dims_padded(i) = this%lookup%table_dims_padded(i) + 1
      end do
    end do
    this%lookup%table_dims_padded_flat = &
      product(this%lookup%table_dims_padded(1:ubound(this%lookup%table_dims_padded, dim=1) - 1))
    this%lookup%table_dims_padded_flat = product(this%lookup%table_dims_padded)
    ! Create a new padded table
    deallocate (this%lookup%elems)
    allocate (this%lookup%elems(this%lookup%table_dim_svar, this%lookup%table_dims_padded_flat))
    this%lookup%elems = 0.d0

      !! Create gold standard and fill table with Alya-format sorted version of table
    call this%create_test_tables_padded()
      !! Run partition mapping algorithm
    call this%lookup%partition_remap(this%part_dims, this%lookup%table_dims_padded)

    ! Open, read, and close the file to the gold standard
    allocate (elems_gold_std(this%lookup%table_dim_svar, this%lookup%table_dims_padded_flat))
    open (36, file='partition_test_table.pad.tmp.dat', action='read')
    N = size(this%part_dims)
    do i = 1, this%lookup%table_dims_padded_flat
      read (36, *) (coord(j), j=1, N), elems_gold_std(:, i)
    end do
    close (36)

      !! Check if the partition mapping of the elements is equal to the gold standard generated in create_test_tables
    if (all(this%lookup%elems .eq. elems_gold_std)) then
      write (*, *) "n = [", this%lookup%table_dims, "], q = [", &
        this%part_dims, "]: ", "partition_map_test passed! Tables will be deleted."
      !call execute_command_line('rm partition_test_table.tmp.dat partition_test_table_sorted.tmp.dat')
    else if (all(abs(this%lookup%elems - elems_gold_std) .lt. 0.0005)) then
      write (*, *) "n = [", this%lookup%table_dims, "], q = [", &
        this%part_dims, "]: ", "trivially small diff in partition_map_test. most likely a pass. Tables not deleted."
    else
      write (*, *) "n = [", this%lookup%table_dims, "], q = [", &
        this%part_dims, "]: ", "partition_map_test conditional not working right. Tables not deleted."
    end if

    deallocate (elems_gold_std)

  end subroutine partition_map_test

  !> Test to map and unmap a flamelet table to a partitioning scheme.
  !! @param this the partition_test object
  subroutine partition_map_unmap_test(this)
    class(partitioning_test), intent(inout) :: this
    double precision, allocatable, dimension(:, :) :: elems_gold_std
    integer :: i, j, N
    integer :: coord(size(this%part_dims))

    print *, "partition_map_unmap_test begin"

    N = size(this%part_dims)
      !! Create gold standard and fill table with Alya-format sorted version of table
    call this%create_test_tables_unpadded()
      !! Run partition mapping algorithm
    call this%lookup%partition_remap(this%part_dims, this%lookup%table_dims)
    call this%lookup%partition_remap(this%lookup%table_dims, this%part_dims)

      !! Open, read, and close the file to the gold standard
    allocate (elems_gold_std(this%lookup%table_dim_svar, this%lookup%table_dims_flat))
    open (36, file='partition_test_table_sorted.nopad.tmp.dat', action='read')
    do i = 1, this%lookup%table_dims_flat
      read (36, *) (coord(j), j=1, N), elems_gold_std(:, i)
    end do
    close (36)

      !! Check if the partition mapping of the elements is equal to the gold standard generated in create_test_tables
    if (all(this%lookup%elems .eq. elems_gold_std)) then
      write (*, *) "n = [", this%lookup%table_dims, "], q = [", &
        this%part_dims, "]: ", "partition_map_unmap_test passed!"
    else if (all(abs(this%lookup%elems - elems_gold_std) .lt. 0.0005)) then
      write (*, *) "n = [", this%lookup%table_dims, "], q = [", &
        this%part_dims, "]: ", "trivially small diff in partition_map_unmap_test. most likely a pass."
    else
      write (*, *) "n = [", this%lookup%table_dims, "], q = [", &
        this%part_dims, "]: ", "partition_map_unmap_test conditional not working right."
    end if

    deallocate (elems_gold_std)

  end subroutine partition_map_unmap_test

  !> Writes to file the gold-standard table for the partitioning test,
  !! which is in the correct partition-major ordering, with easy-to-
  !! understand values for the control and state variables:
  !!
  !! phi_i1 phi_i2 ... phi_iN Eta_j
  !!
  !! where:
  !!
  !! phi_i1 is the coordinate in direction 1
  !! phi_i2 is the coordinate in direction 2
  !! ... ... ...
  !! phi_iN is the coordinate in direction N
  !! Eta_j is an ascending integer.
  !!
  !! The file that has been written is sorted with a GNU awk command,
  !! which yields the order in which Alya-formatted flamelet tables are stored.
  !! @param this the partition_test object to which create_test_tables belongs
  subroutine create_test_tables_padded(this)
    class(partitioning_test), intent(inout) :: this
    integer, allocatable, dimension(:) :: box_dims
    integer :: svar, N, i, j, k
    integer, dimension(size(this%part_dims)) :: coord

    N = size(this%part_dims)
      !! Find total partitions in each dimension
    allocate (box_dims(N))
    box_dims = this%lookup%table_dims_padded(1:N) / this%part_dims

    open (34, file='partition_test_table.pad.tmp.dat', action='write')

    do i = 1, this%lookup%table_dims_padded_flat
      coord = this%lookup%index_to_global_coord(i, this%part_dims, box_dims)
      svar = i
      if (any(coord .gt. this%lookup%table_dims(1:N))) svar = 0
      write (34, *) (coord(j), j=1, N), (svar, k=1, this%lookup%table_dim_svar)
    end do

    close (34)
    deallocate (box_dims)

      !! This gawk loop will sort the coordinate columns 1 to N in order, which
      !! emulates the storage format of Alya tables.
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
                              &  }' partition_test_table.pad.tmp.dat >      &
                              &  partition_test_table_sorted.pad.tmp.dat")

    open (35, file='partition_test_table_sorted.pad.tmp.dat', action='read')

      !! After sorting to the "Alya format," load the table back into the elements array.
    do i = 1, this%lookup%table_dims_padded_flat
      read (35, *) (coord(j), j=1, N), this%lookup%elems(:, i)
    end do

    close (35)

  end subroutine create_test_tables_padded

  !> Writes to file the gold-standard table for the partitioning test,
  !! which is in the correct partition-major ordering, with easy-to-
  !! understand values for the control and state variables:
  !!
  !! phi_i1 phi_i2 ... phi_iN Eta_j
  !!
  !! where:
  !!
  !! phi_i1 is the coordinate in direction 1
  !! phi_i2 is the coordinate in direction 2
  !! ... ... ...
  !! phi_iN is the coordinate in direction N
  !! Eta_j is an ascending integer.
  !!
  !! The file that has been written is sorted with a GNU awk command,
  !! which yields the order in which Alya-formatted flamelet tables are stored.
  !! @param this the partition_test object to which create_test_tables belongs
  subroutine create_test_tables_unpadded(this)
    class(partitioning_test), intent(inout) :: this
    integer, allocatable, dimension(:) :: box_dims
    integer :: svar, N, i, j, k
    integer, dimension(size(this%part_dims)) :: coord

    N = size(this%part_dims)
      !! Find total partitions in each dimension
    allocate (box_dims(N))
    box_dims = this%lookup%table_dims_padded(1:N) / this%part_dims

    open (34, file='partition_test_table.nopad.tmp.dat', action='write')

    do i = 1, this%lookup%table_dims_flat
      coord = this%lookup%index_to_global_coord(i, this%part_dims, box_dims)
      svar = i
      write (34, *) (coord(j), j=1, N), (svar, k=1, this%lookup%table_dim_svar)
    end do

    close (34)
    deallocate (box_dims)

      !! This gawk loop will sort the coordinate columns 1 to N in order, which
      !! emulates the storage format of Alya tables.
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
                              &  }' partition_test_table.nopad.tmp.dat >    &
                              &  partition_test_table_sorted.nopad.tmp.dat")

    open (35, file='partition_test_table_sorted.nopad.tmp.dat', action='read')

      !! After sorting to the "Alya format," load the table back into the elements array.
    do i = 1, this%lookup%table_dims_flat
      read (35, *) (coord(j), j=1, N), this%lookup%elems(:, i)
    end do

    close (35)

  end subroutine create_test_tables_unpadded

  !> Runs the partition mapping test.
  !! @param this the partitioning_test object to which run_test belongs
  subroutine run_map_test(this)
    class(partitioning_test), intent(inout) :: this

    call this%partition_map_test()

  end subroutine run_map_test

  !> Runs the partition mapping-unmapping test.
  !! @param this the partitioning_test object to which run_test belongs
  subroutine run_map_unmap_test(this)
    class(partitioning_test), intent(inout) :: this

    call this%partition_map_unmap_test()

  end subroutine run_map_unmap_test

end module disttab_test_partitioning
