module disttab_test_partitioning
   use :: disttab_table

   implicit none
   private
   public :: partitioning_test

   type :: partitioning_test
      private
      type(table) :: lookup
      integer, allocatable, dimension(:) :: table_dims
      integer, allocatable, dimension(:) :: partition_dims
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
      integer, dimension(:), intent(in) :: table_dimensions
      integer, dimension(:), intent(in) :: partition_dimensions

      allocate (this%table_dims(size(table_dimensions)))
      allocate (this%partition_dims(size(table_dimensions)))

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
      integer, allocatable, dimension(:) :: counters, uppers
      double precision, allocatable, dimension(:, :) :: elements_gold_std
      integer :: i, j, k, l, idx, offset, m, n, part_ctr, h, o, coords(size(this%partition_dims))

      print *, "partition_mapping_test begin"

      ! Pad out table to maintain shape
      this%lookup%table_dims_padded = this%table_dims
      do i = lbound(this%lookup%table_dims_padded, dim=1), ubound(this%lookup%table_dims_padded, dim=1) - 1
         do while (mod(this%lookup%table_dims_padded(i), this%partition_dims(i)) .ne. 0)
            this%lookup%table_dims_padded(i) = this%lookup%table_dims_padded(i) + 1
         end do
      end do
      this%lookup%table_dims_padded_cvar_flat = &
         product(this%lookup%table_dims_padded(1:ubound(this%lookup%table_dims_padded, dim=1) - 1))
      this%lookup%table_dims_padded_flat = product(this%lookup%table_dims_padded)
      deallocate (this%lookup%elements)
      allocate (this%lookup%elements(this%lookup%table_dims_padded_cvar_flat, this%lookup%table_dim_svar))

      N = size(this%partition_dims)
      call part_test_fill_table(this)
      call this%lookup%partition_mapping(this%partition_dims)

      allocate (elements_gold_std(this%lookup%table_dims_padded_cvar_flat, this%lookup%table_dim_svar))

      open (36, file='partition_test_table.tmp.dat', action='read')

      do i = 1, this%lookup%table_dims_padded_cvar_flat
         read (36, *) (coords(j), j=1, N), elements_gold_std(i, 1)
      end do

      close (36)

      !call execute_command_line('rm partition_test_table.tmp.dat partition_test_table_sorted.tmp.dat')

      if (all(this%lookup%elements .eq. elements_gold_std)) then
         write (*, *) "n = [", this%lookup%table_dims, "], q = [", &
            this%partition_dims, "]: ", "partition_mapping_test passed!"
      else if (all(abs(this%lookup%elements - elements_gold_std) .lt. 0.0005)) then
         write (*, *) "n = [", this%lookup%table_dims, "], q = [", &
            this%partition_dims, "]: ", "trivially small diff in partition_mapping_test. most likely a pass."
      else
         write (*, *) "n = [", this%lookup%table_dims, "], q = [", &
            this%partition_dims, "]: ", "partition_mapping_test conditional not working right..."
!         write (*, *) "Sorted table    Gold std"
!         write (*, fmt='(2(f9.3))') abs(this%lookup%elements), abs(elements_gold_std)
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
      integer, allocatable, dimension(:) :: total_partitions
      integer :: h, N, i, ind_local, j, k, l, m, part_ctr
      integer, dimension(size(this%partition_dims)) :: coords, coords_p, coords_b
      character(len=400) :: sort_cmd

      N = size(this%partition_dims)
      allocate (total_partitions(N))
      total_partitions = this%lookup%table_dims_padded(1:N)/this%partition_dims

      open (34, file='partition_test_table.tmp.dat', action='write')

      do i = 1, this%lookup%table_dims_padded_cvar_flat
         part_ctr = ceiling(real(i)/product(this%partition_dims))

         coords_b = this%lookup%flat2coords(part_ctr, total_partitions)

         !! Localized intra-partition index
         if (mod(i, product(this%partition_dims)) .ne. 0) then
            ind_local = mod(i, product(this%partition_dims))
         else
            ind_local = product(this%partition_dims)
         end if

         coords_p = this%lookup%flat2coords(ind_local, this%partition_dims)

         coords = coords_p + this%partition_dims*(coords_b - 1)
         h = i
         if (any(coords .gt. this%table_dims)) h = 0
         write (34, *) (coords(j), j=1, N), h
      end do

      deallocate (total_partitions)

      close (34)

      call execute_command_line("gawk -F ',' '                                &
                                    {                                         &
                                       for(i=1;i<=NF;i++){sorter[i][NR]=$i}   &
                                    }                                         &
                                    END{                                      &
                                       for(i=1;i<=NF;i++){asort(sorter[i])}   &
                                       for(j=1;j<=NR;j++){                    &
                                           for(i=1;i<NF;i++){                 &
                                               printf ""%s,"",sorter[i][j]    &
                                           }                                  &
                                           print sorter[i][j]                 &
                                       }                                      &
                                   }' partition_test_table.tmp.dat > partition_test_table_sorted.tmp.dat")

      open (35, file='partition_test_table_sorted.tmp.dat', action='read')

      do i = 1, this%lookup%table_dims_padded_cvar_flat
         read (35, *) (coords(j), j=1, N), this%lookup%elements(i, 1)
      end do

      close (35)

   end subroutine part_test_fill_table

   !> Runs the partitioning test.
  !! @param this the partitioning_test object to which run_test belongs
   subroutine run_test(this)
      class(partitioning_test), intent(inout) :: this

      print *, "calling partition_mapping_test"
      call partition_mapping_test(this)

   end subroutine run_test

end module disttab_test_partitioning
