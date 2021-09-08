!> A driver program for the DistTab partitioning test.
program test
  use :: disttab_test_access
  use :: disttab_test_partitioning
  use :: disttab_table
  use :: iso_fortran_env

  implicit none
  real (kind = real64), allocatable, dimension(:) :: nr, qr
  integer (kind = int64), allocatable, dimension(:) :: n, q, M
  integer (kind = int64) :: dims, i, j
  integer (kind = int64) :: runs
  type(partitioning_test) :: test_partition
  type(access_test) :: test_access
  type(table) :: lookup

  ! Built-in tests that verify mapping and padding
  call square_test()
  call rand_test_full()
  ! Other tests
  !call read_test()
  !call locality_test()

contains

  !> Simple test of 4x4 table with 2x2 partitions.
  !! Good for quick verification of functionality.
  subroutine square_test()
    allocate (n(3))
    allocate (q(2))
    n = (/4, 4, 1/)
    q = (/2, 2/)
    test_partition = partitioning_test(n, q)
    call test_partition%run_map_test()
    test_partition = partitioning_test(n, q)
    call test_partition%run_map_unmap_test()
    test_access = access_test(n)
    call test_access%run_value_test()
    call test_access%run_value_cloud_test()
    call test_access%run_get_map_get_test(q)
    !call test_access%run_get_perf_test(runs)
    deallocate (q)
    deallocate (n)
  end subroutine square_test

  !> Test of 2-dim to 5-dim partitioning with randomly generated
  !! table and partition sizes on each dimension.
  subroutine rand_test_full()
    do dims = 2, 5
      do j = 1, 2
        allocate (nr(dims + 1))
        allocate (qr(dims))
        allocate (n(dims + 1))
        allocate (q(dims))
        call random_number(qr)
        qr = qr * 12.0
        q = ceiling(qr)
        call random_number(nr)
        nr = nr * 20.0
        n = ceiling(nr)
        n(1:dims) = q + n(1:dims)
        n(dims + 1) = 1
        test_partition = partitioning_test(n, q)
        call test_partition%run_map_test()
        test_partition = partitioning_test(n, q)
        call test_partition%run_map_unmap_test()
        test_access = access_test(n)
        call test_access%run_value_test()
        call test_access%run_value_cloud_test()
        call test_access%run_get_map_get_test(q)
        !allocate (M(dims))
        !i = 10
        !n(1:dims) = 10 * (2**i)
        !n(dims + 1) = 1
        !M = 2**j
        !if (M(1) .gt. n(1)) exit
        !print *, "table dims = ", n(1)
        !print *, "number of segments = ", M(1)
        !call test_access%run_get_perf_test(runs, M)
        !deallocate (M)
        deallocate (nr)
        deallocate (qr)
        deallocate (n)
        deallocate (q)
      end do
    end do
  end subroutine rand_test_full

  !> Read in a table of state variables (more than 1 per line) into table object
  !! todo: Read table in Alya format properly that stores table size, control variables, etc.
  subroutine read_test()
    character(len=50) :: file_id
    integer (kind = int64) :: q_prev(3)
    real (kind = real64) :: interp_random_cvars(2)
    allocate (nr(3))
    allocate (qr(2))
    allocate (n(3))
    allocate (q(2))
    file_id = "../results/cdf-mesh-0-1.000000.dat"
    print *, "the read_test"
    n = (/600, 600, 1/)
    q_prev = 1
    !call random_number(qr)
    !call random_number(interp_random_cvars)
    !qr = qr * 9.0
    !q = ceiling(qr)
    q = (/60, 60, 1/)
    lookup = table(n)
    call lookup%read_in(file_id)
    call lookup%partition_remap(q, q_prev)
    !call lookup%partition_remap(q_prev, q)
    deallocate (nr)
    deallocate (qr)
    deallocate (n)
    deallocate (q)
  end subroutine read_test

  subroutine locality_test()
    allocate (n(4))
    n = (/1000, 1000, 1000, 1/)
    test_access = access_test(n)
    runs = 100000000
    call test_access%run_locality_test(runs)
    deallocate (n)
  end subroutine locality_test

end program test
