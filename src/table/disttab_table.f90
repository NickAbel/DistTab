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
      !procedure, public, pass(this)  :: get_svar_by_coords    ! Given integer coordinates in CV space, return corresponding SV's

      procedure, private, pass(this) :: coords2flat            ! Convert table coordinates to flat-major index
      procedure, private, pass(this) :: get_partition_bounds   ! Get the bounds of the partition containing the argument
      final :: table_destructor

   end type table

   interface table
      module procedure :: table_constructor
   end interface table

contains

   !> Constructor for the table object.  
   !! Allocates the elements array according to the table dimensions.
   !! Initializes the member variables table_dims, table_dims_padded.
   !! Computes the parameters table_dims_cvar_flat, table_dims_flat,
   !! table_dims_padded_cvar_flat, table_dims_padded_flat, table_dim_svar,
   !! which are all a property of the table_dimensions argument, necessary
   !! for partition mapping functionality.
   !! 
   !! @param table_dimensions the length of the control variable space in each dimension, and the number of state variables
   !! @result this table object
   !! @todo if partition dimensions are a property of a lookup table which
   !! may change, the partition_dimensions can be stored as a member variable
   !! but perhaps not initialized here. Then, the same is true of the "padded"
   !! variables which generally will change as the partition_dimensions change.
   !! Organize the code as such. Perhaps use an optional argument for partition_dimensions
   !! here instead.
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

   !> destructor for the table type
   !! 
   !! @param this the table object to destruct
   !! @todo what is this exactly doing? is table_deallocate call necessary?
   subroutine table_destructor(this)
      type(table) :: this

      call table_deallocate(this)
      print *, "If you are reading this, table_destructor is running automatically."

   end subroutine table_destructor

   !> Reads in a lookup table to the elements array of the table object.
   !! Currently reads in tables with the following format:
   !!
   !! phi_1,1 phi_2,1 Phi_1
   !! phi_1,2 phi_2,2 Phi_2
   !! ...     ...     ...
   !!
   !! Where phi_i,j are control variables, and Phi_k are state variables.
   !! @param this table object to which read_in is a member.
   !! @param file_id the lookup table filename to read in.
   !! @todo Should read in a table in the "Alya format" (see tables/2d.dat)
   subroutine read_in(this, file_id)
      class(table), intent(inout) :: this
      character(len=*), intent(in) :: file_id
      integer :: dims(2), i
      double precision :: c1, c2

      open (1, file=file_id, action='read')
      do i = 1, this%table_dims_cvar_flat
         read (unit=1, fmt=*) c1, c2, this%elements(i, 1)
      end do
      close (1)

   end subroutine read_in


   !> A quick and dirty table filler that populates the table with integers 1, 2, ...
   !! 
   !! @param this the table object to fill
   !! @todo for N-dimensional lookup tables this is the first functionality
   !! to develop. Note the commented function prototype below.
   !recursive subroutine fill_example(this,loop_counters,loop_uppers,idx)
   subroutine fill_example(this)
      class(table), intent(inout) :: this
      integer :: i, j, k

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

   !> Remaps the partition from assumed Alya-format ordering to given partition
   !! ordering.
   !! 
   !! @param this table object to perform partition mapping
   !! @param partition_dims partition dimensions to use
   !! @todo zero padding
   !! @todo n-dimensionalization
   !! @todo assumed currently that table is in Alya-format. Add functionality
   !! to remap from one partition size to another. This can be achieved by
   !! writing and sorting the table as in test_partitioning%part_test_fill_table
   !! to recover the Alya format and then reading that sorted file back in.
   !! @todo partition_dims should also be stored as a member variable, initialized
   !! as [1,1,...,1] when in Alya format.
   subroutine partition_mapping(this, partition_dims)
      class(table), intent(inout) :: this
      integer, dimension(2), intent(in) :: partition_dims
      integer, dimension(4) :: partition_bounds
      integer :: offset, i, j, k, l
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

   !> Write elements to stdout. 
   !! 
   !! @param this table object whose elements are to be written
   subroutine print_elements(this)
      class(table), intent(inout) :: this

      write (*, fmt='(f9.3)') this%elements

   end subroutine print_elements

   !> Deallocates the allocated (allocatable) member variables of the table object
   !! 
   !! @param this the table object whose allocatable member variables are to be deallocated
   !! @todo This is called by the destructor, but is it necessary?
   subroutine table_deallocate(this)
      class(table), intent(inout) :: this

      if (allocated(this%table_dims))        deallocate (this%table_dims)
      if (allocated(this%table_dims_padded)) deallocate (this%table_dims_padded)
      if (allocated(this%elements))          deallocate (this%elements)

   end subroutine table_deallocate

   !> Converts entry (i, j) in coordinate indexing to flat index n
   !! 
   !! @param this table object to which coords2flat belongs
   !! @param coords the coordinates of the entry to be returned in flat index
   !! @result flat the flat index of entry located at coordinates coords
   function coords2flat(this, coords) result(flat)
      class(table), intent(inout) :: this
      integer, dimension(2), intent(in)  :: coords
      integer :: flat

      flat = (coords(2) - 1)*this%table_dims(1) + coords(1)

   end function

   !> Return the dimensional coordinate bounds of the block beginning at
   !! given coordinates coords.
   !! Currently only works on partition bases. The commented code is broken.
   !!
   !! @param this the table object to which get_partition_bounds belongs
   !! @param coords the entry (partition base only) whose four bounds is to be returned
   !! @param partition_dims dimensions of the partition in each direction
   !! @result partition_bounds coordinates of the entries which bound the partition containing coords
   !! @todo fix the broken commented code which works on entries which are not partition bases
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
