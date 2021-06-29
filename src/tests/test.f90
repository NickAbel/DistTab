!> A driver program for the DistTab partitioning test.
program test
  use disttab_test_partitioning
  implicit none
  integer :: n(3), q(2)
  type(partitioning_test) :: test_partition_2d
  type(partitioning_test) :: test_partition_2d_0pad
  type(partitioning_test) :: test_partition_nd
  type(partitioning_test) :: test_partition_nd_0pad

  n = (/25,10,1/)
  q = (/5,5/)

  test_partition_2d = partitioning_test(n, q)
  call test_partition_2d%run_test_2d()

  n = (/3,3,1/)
  q = (/2,2/)

  test_partition_2d_0pad = partitioning_test(n, q)
  call test_partition_2d_0pad%run_test_2d_0pad()

  !n = (/4,4,4,1/)
  !q = (/2,2,2/)

  !test_partition_nd = partitioning_test(n, q)
  !call test_partition_nd%run_test()

  !n = (/3,3,3,1/)
  !q = (/2,2,2/)

  !test_partition_nd_0pad = partitioning_test(n, q)
  !call test_partition_nd_0pad%run_test()

end program test
