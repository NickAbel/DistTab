!> A driver program for the DistTab partitioning test.
program test
   use disttab_test_access
   use disttab_test_partitioning
   use disttab_table
   implicit none
   real, allocatable, dimension(:) :: nr, qr
   integer, allocatable, dimension(:) :: n, q
   integer :: dims, i
   type(partitioning_test) :: test_partition
   type(table) :: lookup

   !! Built-in tests that verify mapping and padding
   call square_test()
   call rand_test_full()
   !! Other tests
   !call read_test()

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
      deallocate (q)
      deallocate (n)
   end subroutine square_test

   !> Test of 2-dim to 5-dim partitioning with randomly generated
  !! table and partition sizes on each dimension.
   subroutine rand_test_full()
      do dims = 2, 3
         do i = 1, 5
            allocate (nr(dims + 1))
            allocate (qr(dims))
            allocate (n(dims + 1))
            allocate (q(dims))
            call random_number(qr)
            qr = qr*10.0
            q = ceiling(qr)
            call random_number(nr)
            nr = nr*20.0
            n = ceiling(nr)
            n(1:dims) = q + n(1:dims)
            n(dims + 1) = 1
            test_partition = partitioning_test(n, q)
            call test_partition%run_map_test()
            test_partition = partitioning_test(n, q)
            call test_partition%run_map_unmap_test()
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
      integer :: q_prev(3)
      allocate (nr(3))
      allocate (qr(2))
      allocate (n(3))
      allocate (q(2))
      file_id = "2d_sv.dat"
      print *, "the read_test"
      n = (/121, 11, 15/)
      q_prev = 1
      call random_number(qr)
      qr = qr*9.0
      q = ceiling(qr)
      lookup = table(n)
      call lookup%read_in(file_id)
      call lookup%partition_remap(q, q_prev)
      call lookup%partition_remap(q_prev, q)
      deallocate (nr)
      deallocate (qr)
      deallocate (n)
      deallocate (q)
   end subroutine read_test

end program test
