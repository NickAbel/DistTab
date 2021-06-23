module disttab_test_partitioning
   use :: disttab_table

   implicit none
   private
   public :: partitioning_test

   type :: partitioning_test
      private
      type(table) :: lookup
      integer :: table_dims(3), partition_dims(2)
   contains
      procedure, pass(this) :: run_test
      procedure, pass(this), private :: partition_mapping_test
      procedure, pass(this), private :: part_test_fill_table

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
      integer, intent(in) :: table_dimensions(3), partition_dimensions(2)
      this%table_dims = table_dimensions
      this%partition_dims = partition_dimensions
      this%lookup = table(this%table_dims)

   end function partitioning_test_constructor

   !> A verification test for the table object's partition mapping algorithm.
   !! First, calls part_test_fill_table to fill the elements of this%lookup
   !! with pre-partitioned entries.
   !! Then, calls this%lookup%partition_mapping to get the result of the
   !! partition mapping algorithm.
   !! Next creates the array elements_gold_std, which is a reference
   !! solution to the partition mapping problem with table size of 
   !! this%lookup and the partition dimensions of this%partition_dims.
   !! 
   !! The entries of elements_gold_std are selected to be easily understandable
   !! values Phi_k. In particular,
   !! Phi_k is the sum of the number of the partition which the entry is located,
   !! and the local number of the entry in the partition divided by 1000.
   !! For example, for block number 2 and entry number 3, Phi_k = 2.003.
   !! 
   !! elements_gold_std is then compared to this%lookup%elements. If all
   !! entries are equal, a pass is reported. If there is no difference in
   !! absolute value greater than 0.0005 between all elements of the two tables,
   !! a different message is given (as this could be due to precision issues
   !! that I am not sure don't exist.) Else, the test is failed.
   !!
   !! @param this the partition_test object to which partition_mapping_test belongs
   !! @todo be sure there is no precision issue
   subroutine partition_mapping_test(this)
      class(partitioning_test), intent(inout) :: this
      integer, dimension(2) :: total_partitions
      double precision, allocatable, dimension(:, :) :: elements_gold_std
      integer :: i, j, k, l, m, n, part_ctr

      call part_test_fill_table(this)
      call this%lookup%partition_mapping(this%partition_dims)

      allocate (elements_gold_std(this%lookup%table_dims_cvar_flat, this%lookup%table_dim_svar))

      total_partitions = this%lookup%table_dims(1:2)/this%partition_dims(1:2)

      m = 0
      n = 1
      part_ctr = 1
      do i = 1, total_partitions(2)
         do j = 1, total_partitions(1)
            do k = 1, this%partition_dims(2)
               do l = 1, this%partition_dims(1)
                  m = m + 1
                  elements_gold_std(n, 1) = DBLE(part_ctr) + (DBLE(m)/DBLE(1000))
                  n = n + 1
               end do
            end do
            part_ctr = part_ctr + 1
            m = 0
         end do
      end do

      if (all(this%lookup%elements .eq. elements_gold_std)) then
         write (*, '(a5, 3i3, a8, 2i3, a3, a30) ') "n = [", this%lookup%table_dims, "], q = [", &
            this%partition_dims, "]: ", "partition_mapping_test passed!"
      else if (all(abs(this%lookup%elements - elements_gold_std) .lt. 0.0005)) then
         write (*, '(a5, 3i3, a8, 2i3, a3, a70) ') "n = [", this%lookup%table_dims, "], q = [", &
            this%partition_dims, "]: ", "trivially small diff in partition_mapping_test. most likely a pass."
      else
         write (*, '(a5, 3i3, a8, 2i3, a3, a55) ') "n = [", this%lookup%table_dims, "], q = [", &
            this%partition_dims, "]: ", "partition_mapping_test conditional not working right..."
         !write (*, fmt='(f9.3)') abs(elements_gold_std - this%lookup%elements)
      end if

      deallocate (elements_gold_std)

   end subroutine partition_mapping_test

   !> Writes to file the gold-standard table for the partitioning test,
   !! which is in the correct partition-major ordering, with easy-to-
   !! understand values for the control and state variables:
   !! 
   !! phi_i phi_j Phi_k
   !! 
   !! where:
   !! 
   !! phi_i is the coordinate in direction 1
   !! phi_j is the coordinate in direction 2
   !! Phi_k is the sum of the number of the partition which the entry is located,
   !! and the local number of the entry in the partition divided by 1000.
   !! For example, for block number 2 and entry number 3, Phi_k = 2.003.
   !! 
   !! The file that has been written is sorted according to the phi_j column first,
   !! and then by the phi_i column second, which yields the order in which Alya-
   !! formatted flamelet tables are stored.
   !! @param this the partition_test object to which part_test_fill_table belongs
   subroutine part_test_fill_table(this)
      class(partitioning_test), intent(inout) :: this
      integer, dimension(2) :: total_partitions
      integer :: i, j, k, l, m, part_ctr, c1, c2

      total_partitions = this%lookup%table_dims(1:2)/this%partition_dims(1:2)

      m = 0
      part_ctr = 1
      open (34, file='partition_test_table.tmp.dat', action='write')

      do i = 1, total_partitions(2)
         do j = 1, total_partitions(1)
            do k = 1, this%partition_dims(2)
               do l = 1, this%partition_dims(1)
                  m = m + 1
                  write (unit=34, fmt='(2(i4.3), f9.3)') l + (j - 1)*this%partition_dims(1), &
                     k + (i - 1)*this%partition_dims(2), &
                     part_ctr + (DBLE(m)/DBLE(1000))
               end do
            end do
            part_ctr = part_ctr + 1
            m = 0
         end do
      end do

      close (34)

      call execute_command_line('sort -k 2,2 partition_test_table.tmp.dat >&
  &       partition_test_table_sorted.tmp.dat')

      open (35, file='partition_test_table_sorted.tmp.dat', action='read')

      do i = 1, this%lookup%table_dims_cvar_flat
         read (unit=35, fmt='(2(i4.3), f9.3)') c1, c2, this%lookup%elements(i, 1)
      end do

      close (35)
      call execute_command_line('rm partition_test_table.tmp.dat partition_test_table_sorted.tmp.dat')

   end subroutine part_test_fill_table

   !> Runs the partitioning test.
   !! @param this the partitioning_test object to which run_test belongs
   subroutine run_test(this)
      class(partitioning_test), intent(inout) :: this

      call partition_mapping_test(this)

   end subroutine run_test

end module disttab_test_partitioning
