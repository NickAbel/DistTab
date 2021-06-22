module disttab_table

   !use :: mpi

   implicit none

   private
   public :: table

   type :: table

      double precision, allocatable, dimension(:, :) :: elements
      integer, allocatable, dimension(:) :: table_dims
      integer, allocatable, dimension(:) :: table_dims_padded

      integer :: table_dims_cvar_flat
      integer :: table_dims_flat      
      integer :: table_dims_padded_cvar_flat
      integer :: table_dims_padded_flat
      integer :: table_dim_svar

   contains

      procedure, public, pass(this)  :: read_in                ! Read in a lookup table from file
      procedure, public, pass(this)  :: partition_mapping      ! Map table to partition-major order
      procedure, public, pass(this)  :: print_elements         ! Print the table elements
      procedure, public, pass(this)  :: table_deallocate       ! Deallocate table member variables (this is not the dtor!)

      procedure, private, pass(this) :: coords2flat            ! Convert table coordinates to flat-major index
      procedure, private, pass(this) :: get_partition_bounds   ! Get the bounds of the partition containing the argument
      final :: table_destructor

   end type table

   interface table
      module procedure :: table_constructor
   end interface table

contains

   type(table) function table_constructor(table_dimensions) result(this)
      integer, intent(in) :: table_dimensions(3)
      integer :: i

      allocate (this%table_dims(size(table_dimensions)))
      allocate (this%table_dims_padded(size(table_dimensions)))

      this%table_dims = table_dimensions
      this%table_dims_cvar_flat = product(this%table_dims(1:ubound(this%table_dims, dim=1) - 1))
      this%table_dims_flat = product(this%table_dims)
      this%table_dims_padded_cvar_flat = product(this%table_dims(1:ubound(this%table_dims, dim=1) - 1))
      this%table_dims_padded_flat = product(this%table_dims)
      this%table_dim_svar = this%table_dims(ubound(this%table_dims, dim=1))

      allocate (this%elements(this%table_dims_padded_cvar_flat, this%table_dim_svar))

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

      !call mpi_comm_rank(mpi_comm_world, rank, ierror)

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

   subroutine partition_mapping(this, partition_dims)
      class(table), intent(inout) :: this
      integer, dimension(2), intent(in) :: partition_dims
      integer, dimension(4) :: partition_bounds
      integer :: offset, i, j, k, l
      integer :: rank, ierror
      double precision, allocatable, dimension(:, :) :: elements_old

      !call mpi_comm_rank(mpi_comm_world, rank, ierror)

      allocate (elements_old(this%table_dims_cvar_flat, this%table_dim_svar))

      elements_old = this%elements

      ! Pad out table to maintain shape
      !this%table_dims_padded = this%table_dims
      !do i = lbound(this%table_dims_padded, dim=1), ubound(this%table_dims_padded, dim=1) - 1
      !   do while (mod(this%table_dims_padded(i), partition_dims(i)) .ne. 0)
      !      this%table_dims_padded(i) = this%table_dims_padded(i) + 1
      !   end do
      !end do

      offset = 1

      do i = 1, this%table_dims(2), partition_dims(2)
         do j = 1, this%table_dims(1), partition_dims(1)
            partition_bounds = get_partition_bounds(this, (/j, i/), partition_dims)
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

   subroutine print_elements(this)
      class(table), intent(inout) :: this
      integer :: rank, ierror

      !call mpi_comm_rank(mpi_comm_world, rank, ierror)

      write (*, fmt='(f9.3)') this%elements

   end subroutine print_elements

   subroutine table_deallocate(this)
      class(table), intent(inout) :: this

      if (allocated(this%table_dims))        deallocate (this%table_dims)
      if (allocated(this%table_dims_padded)) deallocate (this%table_dims_padded)
      if (allocated(this%elements))          deallocate (this%elements)

   end subroutine table_deallocate

   function coords2flat(this, coords) result(flat)
      class(table), intent(inout) :: this
      integer, dimension(2), intent(in)  :: coords
      integer :: flat

      flat = (coords(2) - 1)*this%table_dims(1) + coords(1)

   end function

   function get_partition_bounds(this, coords, partition_dims) result(partition_bounds)
      class(table), intent(inout) :: this
      integer, dimension(2), intent(in)  :: coords
      integer, dimension(2) :: partition_dims
      integer, dimension(4) :: partition_bounds

      partition_bounds(1) = coords(1)
      partition_bounds(2) = coords(2)
      partition_bounds(3) = coords(1) + partition_dims(1) - 1
      partition_bounds(4) = coords(2) + partition_dims(2) - 1

    !! This, while more general, is broken here, so I use the above now
    !! which is more simple but only works on partition bases.
    !! Lower bound direction 1
      !partition_bounds(1) = coords(1) - &
      !  merge(mod(coords(1),partition_dims(1))-1,partition_dims(1)-1, &
      !  mod(coords(1),partition_dims(1)) .ne. 0)

    !! Lower bound direction 2
      !partition_bounds(2) = coords(2) - &
      !  merge(mod(coords(2),partition_dims(2))-1,partition_dims(2)-1, &
      !  mod(coords(2),partition_dims(2)) .ne. 0)

    !! Upper bound direction 1
      !partition_bounds(3) = coords(1) - &
      !  merge(partition_dims(1)-mod(coords(1),partition_dims(1)),0, &
      !  mod(coords(1),partition_dims(1)) .ne. 0)

    !! Upper bound direction 2
      !partition_bounds(4) = coords(2) - &
      !  merge(partition_dims(2)-mod(coords(2),partition_dims(2)),0, &
      !  mod(coords(2),partition_dims(2)) .ne. 0)

   end function

end module disttab_table
