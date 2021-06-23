!> A driver program for the DistTab partitioning test.
program test
  use disttab_test_partitioning
  implicit none
  integer :: n(3), q(2)
  type(partitioning_test) :: test_partition

  n = (/4,4,1/)
  q = (/2,2/)

  test_partition = partitioning_test(n, q)
end program test
