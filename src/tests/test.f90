!> A driver program for the DistTab partitioning test.
program test
  use :: disttab_test_access
  use :: disttab_test_partitioning
  use :: disttab_table
  use :: iso_fortran_env

  implicit none
  real (kind = real64), allocatable, dimension(:) :: table_dims_real, part_dims_real
  integer (kind = int64), allocatable, dimension(:) :: table_dims, part_dims, segments
  integer (kind = int64) :: dims, i, j
  integer (kind = int64) :: total_runs
  type(partitioning_test) :: test_partition
  type(access_test) :: test_access
  type(table) :: lookup

  ! Built-in tests that verify mapping and padding
  call square_test()
  call rand_test_full()

  ! Other tests
  !call read_test()
  call locality_test()
  !call access_perf_test()

contains

  !> Simple test of 4x4 table with 2x2 partitions.
  !! Good for quick verification of functionality.
  subroutine square_test()
    allocate (table_dims(3))
    allocate (part_dims(2))

    table_dims = (/4, 4, 1/)
    part_dims = (/2, 2/)

    test_partition = partitioning_test(table_dims, part_dims)
    call test_partition%run_map_test()

    test_partition = partitioning_test(table_dims, part_dims)
    call test_partition%run_map_unmap_test()

    test_access = access_test(table_dims)
    call test_access%run_value_test()
    call test_access%run_value_cloud_test()
    call test_access%run_get_map_get_test(part_dims)

    deallocate (part_dims)
    deallocate (table_dims)
  end subroutine square_test

  !> Test of 2-dim to 5-dim partitioning with randomly generated
  !! table and partition sizes on each dimension.
  subroutine rand_test_full()
    do dims = 2, 5
        allocate (table_dims_real(dims + 1))
        allocate (part_dims_real(dims))
        allocate (table_dims(dims + 1))
        allocate (part_dims(dims))

        call random_number(part_dims_real)
        part_dims_real = part_dims_real * 12.0
        part_dims = ceiling(part_dims_real)

        call random_number(table_dims_real)
        table_dims_real = table_dims_real * 20.0
        table_dims = ceiling(table_dims_real)
        table_dims(1:dims) = part_dims + table_dims(1:dims)
        table_dims(dims + 1) = 1

        test_partition = partitioning_test(table_dims, part_dims)
        call test_partition%run_map_test()

        test_partition = partitioning_test(table_dims, part_dims)
        call test_partition%run_map_unmap_test()

        test_access = access_test(table_dims)
        call test_access%run_value_test()
        call test_access%run_value_cloud_test()
        call test_access%run_get_map_get_test(part_dims)

        deallocate (table_dims_real)
        deallocate (part_dims_real)
        deallocate (table_dims)
        deallocate (part_dims)
    end do
  end subroutine rand_test_full

  !> Read in a table of state variables (more than 1 per line) into table object
  !! todo: Read table in Alya format properly that stores table size, control variables, etc.
  !! todo: Line 495 error when reading in
  subroutine read_test()
    character(len=120) :: file_id
    integer (kind = int64) :: part_dims_prev(3)

    allocate (table_dims_real(3))
    allocate (part_dims_real(2))
    allocate (table_dims(3))
    allocate (part_dims(2))

    file_id = "../results/cdf-mesh-0-1.000000.dat"

    print *, "the read_test"

    table_dims = (/600, 600, 1/)
    part_dims_prev = 1
    part_dims = (/60, 60, 1/)
    lookup = table(table_dims)
    call lookup%read_in(file_id)
    call lookup%partition_remap(part_dims, part_dims_prev)
    !call lookup%partition_remap(part_dims_prev, part_dims)

    deallocate (table_dims_real)
    deallocate (part_dims_real)
    deallocate (table_dims)
    deallocate (part_dims)
  end subroutine read_test

  subroutine locality_test()
    allocate (table_dims(4))

    table_dims = (/1000, 1000, 1000, 1/)
    test_access = access_test(table_dims)
    total_runs = 100000000
    call test_access%run_locality_test(total_runs)

    deallocate (table_dims)
  end subroutine locality_test

  !> Run bucketed real -> coordinate search performance tests.
  subroutine access_perf_test()
    do dims = 2, 2
      do j = 1, 2
        allocate (table_dims_real(dims + 1))
        allocate (part_dims_real(dims))
        allocate (table_dims(dims + 1))
        allocate (part_dims(dims))
        allocate (segments(dims))

        i = 10
        table_dims(1:dims) = 10 * (2**i)
        table_dims(dims + 1) = 1
        segments = 2**j

        if (segments(1) .gt. table_dims(1)) exit

        print *, "table dims = ", table_dims(1)
        print *, "number of segments = ", segments(1)

        test_access = access_test(table_dims)
        call test_access%run_get_perf_test(total_runs, segments)

        deallocate (segments)
        deallocate (table_dims_real)
        deallocate (part_dims_real)
        deallocate (table_dims)
        deallocate (part_dims)
      end do
    end do
  end subroutine access_perf_test


end program test
