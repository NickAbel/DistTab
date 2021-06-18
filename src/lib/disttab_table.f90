module disttab_table

   use :: mpi

   implicit none

   private
   public :: table

   type :: table

      integer, allocatable, dimension(:) :: table_dims
      integer, allocatable, dimension(:) :: partition_dims
      integer, allocatable, dimension(:) :: total_partitions
      double precision, allocatable, dimension(:, :) :: elements
      integer :: table_dims_cvar_flat
      integer :: table_dims_flat
      integer :: partition_dims_flat
      integer :: table_dim_svar

   contains

      procedure, public, pass(this)  :: read_in                ! Read in a lookup table from file
      procedure, public, pass(this)  :: fill_example           ! Fill table with example entries
      procedure, public, pass(this)  :: partition_mapping      ! Map table to partition-major order
      procedure, public, pass(this)  :: partition_mapping_test ! Run test to verify block reorganization functionality
      procedure, public, pass(this)  :: print_elements         ! Print the table elements
      procedure, public, pass(this)  :: table_deallocate       ! Deallocate table member variables (this is not the dtor!)

      procedure, private, pass(this) :: part_test_fill_table   ! Fill table for the partition mapping test (uses file I/O)
      procedure, private, pass(this) :: coords2flat            ! Convert table coordinates to flat-major index
      procedure, private, pass(this) :: get_partition_bounds   ! Get the bounds of the partition containing the argument
      final :: table_destructor

   end type table

   interface table
      module procedure :: table_constructor
   end interface table

contains

   type(table) function table_constructor(table_dimensions, partition_dimensions) result(this)
      integer, intent(in) :: table_dimensions(3), partition_dimensions(2)

      allocate (this%table_dims(size(table_dimensions)))
      allocate (this%partition_dims(size(partition_dimensions)))
      allocate (this%total_partitions(size(partition_dimensions)))

      this%table_dims = table_dimensions
      this%partition_dims = partition_dimensions

      this%total_partitions = this%table_dims(1:2)/this%partition_dims(1:2)

      this%table_dims_cvar_flat = product(this%table_dims(1:ubound(this%table_dims, dim=1) - 1))
      this%table_dims_flat = product(this%table_dims)
      this%partition_dims_flat = product(this%partition_dims)
      this%table_dim_svar = this%table_dims(ubound(this%table_dims, dim=1))

      allocate (this%elements(this%table_dims_cvar_flat, this%table_dim_svar))

   end function table_constructor

   subroutine table_destructor(this)
      type(table) :: this

      call table_deallocate(this)
      print *, "If you are reading this, table_destructor is running automatically."

   end subroutine table_destructor

   subroutine read_in(this, file_id)
      class(table), intent(inout) :: this
      character(len=*), intent(in) :: file_id
      integer :: dims(2), i
      double precision :: c1, c2

      open (1, file=file_id, action='read')
      read (unit=1, fmt=*) dims
      read (unit=1, fmt=*) dims
      do i = 1, this%table_dims_cvar_flat
         read (unit=1, fmt=*) c1, c2, this%elements(i, 1)
      end do
      close (1)

   end subroutine read_in

   !recursive subroutine fill_example(this,loop_counters,loop_uppers,idx)
   subroutine fill_example(this)
      class(table), intent(inout) :: this
      integer :: i, j, k
      integer :: rank, ierror

      call mpi_comm_rank(mpi_comm_world, rank, ierror)

      !integer, dimension(size(this%table_dims)-1) :: loop_counters, loop_uppers
      !integer, dimension(size(this%table_dims)-1) :: loop_counters_copy

      !if (idx .eq. 1) then
      !  this%elements(coords2flat(this,(/loop_counters(1),loop_counters(2)/)) + this%table_dims_cvar_flat*rank,1) = &
      !    coords2flat(this,(/loop_counters(1),loop_counters(2)/)) + 1.d0*this%table_dims_cvar_flat*rank
      !  print *, "?", coords2flat(this,(/loop_counters(1),loop_counters(2)/)), loop_counters
      !  do while (loop_counters(idx) .lt. loop_uppers(idx))
      !    loop_counters(idx) = loop_counters(idx) + 1
      !    this%elements(coords2flat(this,(/loop_counters(1),loop_counters(2)/)) + this%table_dims_cvar_flat*rank,1) = &
      !      coords2flat(this,(/loop_counters(1),loop_counters(2)/)) + 1.d0*this%table_dims_cvar_flat*rank
      !    print *, "?", coords2flat(this,(/loop_counters(1),loop_counters(2)/)), loop_counters
      !  enddo
      !else if (idx .gt. 1) then
      !  do while (loop_counters(idx) .le. loop_uppers(idx))
      !    loop_counters_copy = loop_counters
      !    call fill_example(this,loop_counters_copy,loop_uppers,idx-1,rank)
      !    loop_counters(idx) = loop_counters(idx) + 1
      !  enddo
      !endif

      k = 1
      do i = 1, this%table_dims(1)
         do j = 1, this%table_dims(2)
            this%elements(k, 1) = k + &
                                  this%table_dims_cvar_flat*rank*1.d0
            k = k + 1
         end do
      end do

   end subroutine fill_example

   subroutine partition_mapping(this)
      class(table), intent(inout) :: this
      integer :: offset, i, j, k, l
      integer, dimension(4) :: partition_bounds
      integer :: rank, ierror
      double precision, allocatable, dimension(:, :) :: elements_old

      call mpi_comm_rank(mpi_comm_world, rank, ierror)

      allocate (elements_old(this%table_dims_cvar_flat, this%table_dim_svar))

      elements_old = this%elements

      offset = 1

      do i = 1, this%table_dims(2), this%partition_dims(2)
         do j = 1, this%table_dims(1), this%partition_dims(1)
            partition_bounds = get_partition_bounds(this, (/j, i/))
            do k = partition_bounds(2), partition_bounds(4)
               do l = partition_bounds(1), partition_bounds(3)
                  this%elements(offset, 1) = elements_old(coords2flat(this, (/l, k/)), 1)
                  offset = offset + 1
               end do
            end do
         end do
      end do

      deallocate (elements_old)

   end subroutine partition_mapping

   subroutine partition_mapping_test(this)
      class(table), intent(inout) :: this
      double precision, allocatable, dimension(:, :) :: elements_gold_std
      integer :: i, j, k, l, m, n, part_ctr

      call part_test_fill_table(this)
      call partition_mapping(this)

      allocate (elements_gold_std(this%table_dims_cvar_flat, this%table_dim_svar))

      m = 0
      n = 1
      part_ctr = 1
      do i = 1, this%total_partitions(2)
         do j = 1, this%total_partitions(1)
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

      if (all(this%elements .eq. elements_gold_std)) then
         write (*, '(a5, 3i3, a8, 2i3, a3, a30) ') "n = [", this%table_dims, "], q = [", &
            this%partition_dims, "]: ", "partition_mapping_test passed!"
      else if (all(abs(this%elements - elements_gold_std) .lt. 0.0005)) then
         write (*, '(a5, 3i3, a8, 2i3, a3, a70) ') "n = [", this%table_dims, "], q = [", &
            this%partition_dims, "]: ", "trivially small diff in partition_mapping_test. most likely a pass."
      else
         write (*, '(a5, 3i3, a8, 2i3, a3, a55) ') "n = [", this%table_dims, "], q = [", &
            this%partition_dims, "]: ", "partition_mapping_test conditional not working right..."
         !write (*, fmt='(f9.3)') abs(elements_gold_std - this%elements)
      end if

      deallocate (elements_gold_std)

   end subroutine partition_mapping_test

   subroutine print_elements(this)
      class(table), intent(inout) :: this
      integer :: rank, ierror

      call mpi_comm_rank(mpi_comm_world, rank, ierror)

      write (*, fmt='(f9.3)') this%elements

   end subroutine print_elements

   subroutine table_deallocate(this)
      class(table), intent(inout) :: this

      deallocate (this%table_dims)
      deallocate (this%partition_dims)
      deallocate (this%total_partitions)
      deallocate (this%elements)

   end subroutine table_deallocate

   subroutine part_test_fill_table(this)
      class(table), intent(inout) :: this
      integer :: i, j, k, l, m, part_ctr, c1, c2

      m = 0
      part_ctr = 1
      open (34, file='partition_test_table.tmp.dat', action='write')

      do i = 1, this%total_partitions(2)
         do j = 1, this%total_partitions(1)
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

      do i = 1, this%table_dims_cvar_flat
         read (unit=35, fmt='(2(i4.3), f9.3)') c1, c2, this%elements(i, 1)
      end do

      close (35)
      call execute_command_line('rm partition_test_table.tmp.dat partition_test_table_sorted.tmp.dat')

   end subroutine part_test_fill_table

   function coords2flat(this, coords) result(flat)
      class(table), intent(inout) :: this
      integer, dimension(2), intent(in)  :: coords
      integer :: flat

      flat = (coords(2) - 1)*this%table_dims(1) + coords(1)

   end function

   function get_partition_bounds(this, coords) result(partition_bounds)
      class(table), intent(inout) :: this
      integer, dimension(2), intent(in)  :: coords
      integer, dimension(4) :: partition_bounds

      partition_bounds(1) = coords(1)
      partition_bounds(2) = coords(2)
      partition_bounds(3) = coords(1) + this%partition_dims(1) - 1
      partition_bounds(4) = coords(2) + this%partition_dims(2) - 1

    !! This, while more general, is broken here, so I use the above now
    !! which is more simple but only works on partition bases.
    !! Lower bound direction 1
      !partition_bounds(1) = coords(1) - &
      !  merge(mod(coords(1),this%partition_dims(1))-1,this%partition_dims(1)-1, &
      !  mod(coords(1),this%partition_dims(1)) .ne. 0)

    !! Lower bound direction 2
      !partition_bounds(2) = coords(2) - &
      !  merge(mod(coords(2),this%partition_dims(2))-1,this%partition_dims(2)-1, &
      !  mod(coords(2),this%partition_dims(2)) .ne. 0)

    !! Upper bound direction 1
      !partition_bounds(3) = coords(1) - &
      !  merge(this%partition_dims(1)-mod(coords(1),this%partition_dims(1)),0, &
      !  mod(coords(1),this%partition_dims(1)) .ne. 0)

    !! Upper bound direction 2
      !partition_bounds(4) = coords(2) - &
      !  merge(this%partition_dims(2)-mod(coords(2),this%partition_dims(2)),0, &
      !  mod(coords(2),this%partition_dims(2)) .ne. 0)

   end function

end module disttab_table
