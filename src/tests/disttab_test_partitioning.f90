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

   type(partitioning_test) function partitioning_test_constructor(n, q) result(this)
      integer, intent(in) :: n(3), q(2)
      this%table_dims = n
      this%partition_dims = q
      this%lookup = table(n)

   end function partitioning_test_constructor

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

   subroutine run_test(this)
      class(partitioning_test), intent(inout) :: this

      call partition_mapping_test(this)

   end subroutine run_test

end module disttab_test_partitioning
