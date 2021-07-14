!> A driver program for the DistTab partitioning test.
program test
   use disttab_test_partitioning
   implicit none
   real, allocatable, dimension(:) :: nr, qr
   integer, allocatable, dimension(:) :: n, q
   integer :: dims
   type(partitioning_test) :: test_partition

   call square_test()
   call rand_test_full()

contains

   !> Simple test of 4x4 table with 2x2 partitions.
  !! Good for quick verification of functionality.
   subroutine square_test()
      allocate (n(3))
      allocate (q(2))
      n = (/4, 4, 1/)
      q = (/2, 2/)
      test_partition = partitioning_test(n, q)
      call test_partition%run_test()
      deallocate (n)
      deallocate (q)
   end subroutine square_test

   !> Test of 2-dim to 5-dim partitioning with randomly generated
  !! table and partition sizes on each dimension.
   subroutine rand_test_full()
      do dims = 2, 5
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
         call test_partition%run_test()
         deallocate (nr)
         deallocate (qr)
         deallocate (n)
         deallocate (q)
      end do
   end subroutine rand_test_full

end program test
