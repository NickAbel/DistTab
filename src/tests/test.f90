!> A driver program for DistTab testing.
program test
  use :: disttab_table
  use :: disttab_test_access
  use :: disttab_test_parallel
  use :: disttab_test_partitioning
  use :: kind_params
  use :: mpi

  implicit none
  real(sp), allocatable, dimension(:) :: table_dims_real, part_dims_real
  integer(i4), allocatable, dimension(:) :: table_dims, part_dims, subtable_dims, segments
  integer(i4) :: dims, i, j, ierror
  integer(i4) :: total_runs
  type(table) :: lookup
  type(access_test) :: test_access
  type(parallel_test) :: test_parallel
  type(partitioning_test) :: test_partition

  call mpi_init(ierror)

  ! Built-in tests that verify mapping and padding
  call square_test()
  call rand_test_fast()
  !call rand_test_full()
  !call test_mpi()
  !call reshape_test()

  ! Other tests
  !call read_test() ! This won't work for the time being as MPI is introduced in table % read_in()
  !call locality_test()
  !call access_perf_test()

contains

!> Simple test of 4x4 table with 2x2 partitions.
!! Good for quick verification of functionality.
  subroutine square_test()
    block
      type(partitioning_test) :: test_part
      type(access_test) :: test_acc

      allocate (table_dims(3))
      allocate (part_dims(2))

      table_dims = (/4, 4, 1/)
      part_dims = (/2, 2/)

      test_part = partitioning_test(table_dims, part_dims)
      call test_part % run_map_test()

      test_part = partitioning_test(table_dims, part_dims)
      call test_part % run_map_unmap_test()

      test_acc = access_test(table_dims)
      call test_acc % run_value_test()
      call test_acc % run_value_cloud_test()
      call test_acc % run_get_map_get_test(part_dims)

      deallocate (part_dims)
      deallocate (table_dims)
    end block
  end subroutine square_test

!> Test of 2-dim to 3-dim partitioning with randomly generated
!! table and partition sizes on each dimension.
  subroutine rand_test_fast()
    do dims = 2, 3
      allocate (table_dims_real(dims + 1))
      allocate (part_dims_real(dims))
      allocate (table_dims(dims + 1))
      allocate (part_dims(dims))

      call random_number(part_dims_real)
      part_dims_real = part_dims_real * 15.0
      part_dims = ceiling(part_dims_real)

      call random_number(table_dims_real)
      table_dims_real = table_dims_real * 50.0
      table_dims = ceiling(table_dims_real)
      table_dims(1:dims) = part_dims + table_dims(1:dims)
      table_dims(dims + 1) = 1

      test_partition = partitioning_test(table_dims, part_dims)
      call test_partition % run_map_test()

      test_partition = partitioning_test(table_dims, part_dims)
      call test_partition % run_map_unmap_test()

      test_access = access_test(table_dims)
      call test_access % run_value_test()
      call test_access % run_value_cloud_test()
      call test_access % run_get_map_get_test(part_dims)

      deallocate (table_dims_real)
      deallocate (part_dims_real)
      deallocate (table_dims)
      deallocate (part_dims)
    end do
  end subroutine rand_test_fast

!> Test of 2-dim to 6-dim partitioning with randomly generated
!! table and partition sizes on each dimension.
  subroutine rand_test_full()
    do dims = 2, 4
      allocate (table_dims_real(dims + 1))
      allocate (part_dims_real(dims))
      allocate (table_dims(dims + 1))
      allocate (part_dims(dims))

      call random_number(part_dims_real)
      part_dims_real = part_dims_real * 12.0
      part_dims = ceiling(part_dims_real)

      call random_number(table_dims_real)
      table_dims_real = table_dims_real * 100.0
      table_dims = ceiling(table_dims_real)
      table_dims(1:dims) = part_dims + table_dims(1:dims)
      table_dims(dims + 1) = 1

      test_partition = partitioning_test(table_dims, part_dims)
      call test_partition % run_map_test()

      test_partition = partitioning_test(table_dims, part_dims)
      call test_partition % run_map_unmap_test()

      test_access = access_test(table_dims)
      call test_access % run_value_test()
      call test_access % run_value_cloud_test()
      call test_access % run_get_map_get_test(part_dims)

      deallocate (table_dims_real)
      deallocate (part_dims_real)
      deallocate (table_dims)
      deallocate (part_dims)
    end do
  end subroutine rand_test_full

!> Test MPI capabilities somehow.
  subroutine test_mpi()
    character(len=120) :: file_id

    allocate (table_dims(3))
    allocate (subtable_dims(2))
    allocate (part_dims(2))

    table_dims = (/4, 16, 1/)
    subtable_dims = (/4, 4/)
    part_dims = (/2, 2/)

    test_parallel = parallel_test(table_dims, subtable_dims, part_dims)
    !call test_parallel % run_parallel_get_test()
    call test_parallel % run_parallel_partition_map_test()
    !call test_parallel % run_local_pile_test()
    !call test_parallel % run_parallel_partition_map_unmap_test()

    ! WIP todo parallel file I/O with MPI in lookup % read_in()
    !file_id = "../tables/4x4_1var.dat"
    !call lookup % read_in(file_id)

    deallocate (table_dims)
    deallocate (subtable_dims)
    deallocate (part_dims)

    call mpi_finalize(ierror)

  end subroutine test_mpi

!> Test of 2-dim to 3-dim partitioning with randomly generated
!! table and partition sizes on each dimension.
  subroutine reshape_test()
    type(table) :: lookup_shaped
    real(sp) :: nvar_real
    integer(i4) :: nvar
    print *, "reshape_test: "
    do dims = 2, 3
      lookup_shaped = table((/1/))
      print *, "generated dims: ", dims, " unshaped initialized table dimensions (should be 1): ", lookup_shaped % table_dims
      allocate (table_dims_real(dims + 1))
      allocate (table_dims(dims + 1))

      call random_number(nvar_real)
      nvar_real = nvar_real * 18.0
      nvar = ceiling(nvar_real)

      call random_number(table_dims_real)
      table_dims_real = table_dims_real * 20.0
      table_dims = ceiling(table_dims_real)
      table_dims(1:dims) = 1 + table_dims(1:dims)
      table_dims(dims + 1) = nvar

      print *, "reshaping to (last dim is nvar): ", table_dims

      call lookup_shaped % reshape_table(dims, nvar, table_dims)

      print *, size(lookup_shaped % elems)

      deallocate (table_dims_real)
      deallocate (table_dims)
    end do
  end subroutine reshape_test

!> Read in a table of state variables (more than 1 per line) into table object
!! todo: Read table in Alya format properly that stores table size, control variables, etc.
!! todo: Line 495 error when reading in
  subroutine read_test()
    character(len=120) :: file_id
    integer(i4) :: part_dims_prev(3)

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
    call lookup % read_in(file_id)
    call lookup % partition_remap(part_dims, part_dims_prev)
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
    call test_access % run_locality_test(total_runs)

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
        call test_access % run_get_perf_test(total_runs, segments)

        deallocate (segments)
        deallocate (table_dims_real)
        deallocate (part_dims_real)
        deallocate (table_dims)
        deallocate (part_dims)
      end do
    end do
  end subroutine access_perf_test

end program test
