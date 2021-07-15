module disttab_table

   !use :: mpi

   implicit none

   private
   public :: table

   type :: table

      double precision, allocatable, dimension(:, :) :: elements
      integer, allocatable, dimension(:) :: table_dims
      integer, allocatable, dimension(:) :: table_dims_padded
      integer, allocatable, dimension(:) :: partition_dims

      integer :: table_dims_cvar_flat
      integer :: table_dims_flat
      integer :: table_dims_padded_cvar_flat
      integer :: table_dims_padded_flat
      integer :: table_dim_svar
      integer :: partition_dims_flat

   contains

      procedure, public, pass(this)  :: read_in                ! Read in a lookup table from file
      procedure, public, pass(this)  :: partition_mapping      ! Map table to partition-major order
      procedure, public, pass(this)  :: print_elements         ! Print the table elements
      procedure, public, pass(this)  :: table_deallocate       ! Deallocate table member variables (this is not the dtor!)
      !procedure, public, pass(this)  :: get_svar_by_coords    ! Given integer coordinates in CV space, return corresponding SV's

      procedure, public, pass(this) :: coords2flat             ! Convert table coordinates to flat-major index
      procedure, public, pass(this) :: flat2coords             ! Convert flat-major index to table coordinates
      procedure, public, pass(this) :: get_global_coords       ! ravel a linear index to partitioned index
      procedure, public, pass(this) :: get_local_coords        ! Unravel a partitioned index to linear index
      procedure, public, pass(this) :: get_index_global_coords ! ravel a linear index to partitioned index
      procedure, public, pass(this) :: get_index_local_coords  ! Unravel a partitioned index to linear index
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
      integer, dimension(:), intent(in) :: table_dimensions

      allocate (this%table_dims(size(table_dimensions)))
      allocate (this%table_dims_padded(size(table_dimensions)))
      allocate (this%partition_dims(size(table_dimensions) - 1))

      this%table_dims = table_dimensions
      this%table_dims_cvar_flat = product(this%table_dims(1:ubound(this%table_dims, dim=1) - 1))
      this%table_dims_flat = product(this%table_dims)
      this%table_dims_padded_cvar_flat = product(this%table_dims(1:ubound(this%table_dims, dim=1) - 1))
      this%table_dims_padded_flat = product(this%table_dims)
      this%table_dim_svar = this%table_dims(ubound(this%table_dims, dim=1))

      this%partition_dims = 1 ! Table initially considered to have 1x1 partitions, i.e. 'unpartitioned'
      this%partition_dims_flat = product(this%partition_dims)

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
      integer :: dims(size(this%table_dims) - 1), i
      double precision :: c1, c2

      open (1, file=file_id, action='read')
      do i = 1, this%table_dims_cvar_flat
         read (unit=1, fmt=*) c1, c2, this%elements(i, 1)
      end do
      close (1)

   end subroutine read_in

   !> Remaps the partition from assumed Alya-format ordering to given partition
   !! ordering.
   !!
   !! @param this table object to perform partition mapping
   !! @param partition_dimensions partition dimensions to use
   !! @todo zero padding
   !! @todo n-dimensionalization
   !! @todo assumed currently that table is in Alya-format. Add functionality
   !! to remap from one partition size to another. This can be achieved by
   !! writing and sorting the table as in test_partitioning%part_test_fill_table
   !! to recover the Alya format and then reading that sorted file back in.
   subroutine partition_mapping(this, partition_dimensions)
      class(table), intent(inout) :: this
      integer, dimension(size(this%table_dims) - 1), intent(in) :: partition_dimensions
      integer, dimension(size(this%table_dims) - 1) :: coords, coords_b, coords_p, partition_dimensions_prev
      integer, dimension(size(this%table_dims) - 1) :: total_partitions
      integer :: i, i_new, N
      double precision, allocatable, dimension(:, :) :: elements_old

      !call mpi_comm_rank(mpi_comm_world, rank, ierror)
      N = size(this%table_dims) - 1
      allocate (elements_old(this%table_dims_cvar_flat, this%table_dim_svar))

      elements_old = this%elements
      partition_dimensions_prev = this%partition_dims
      this%partition_dims = partition_dimensions
      this%partition_dims_flat = product(this%partition_dims)

      ! Pad out table to maintain shape
      this%table_dims_padded = this%table_dims
      do i = lbound(this%table_dims_padded, dim=1), ubound(this%table_dims_padded, dim=1) - 1
         do while (mod(this%table_dims_padded(i), partition_dimensions(i)) .ne. 0)
            this%table_dims_padded(i) = this%table_dims_padded(i) + 1
         end do
      end do

      this%table_dims_padded_cvar_flat = &
         product(this%table_dims_padded(1:ubound(this%table_dims_padded, dim=1) - 1))
      this%table_dims_padded_flat = product(this%table_dims_padded)

      deallocate (this%elements)
      allocate (this%elements(this%table_dims_padded_cvar_flat, this%table_dim_svar))
      this%elements = 0.d0

      do i = 1, N
         total_partitions(i) = this%table_dims_padded(i)/this%partition_dims(i)
      end do

      do i = 1, this%table_dims_padded_cvar_flat
         call this%get_local_coords(i, this%partition_dims, total_partitions, coords_p, coords_b)
         coords = coords_p + (coords_b - 1)*this%partition_dims
         !print *, "get_index_global_coords at ", coords, " = ", this%get_index_global_coords(coords)
         !print *, "get_index_local_coords at ", coords_p, " + ", coords_b, " = ", &
         !  this%get_index_local_coords(coords_p, coords_b)
         i_new = this%get_index_global_coords(coords, this%table_dims_padded(1:N))
         this%elements(i, 1) = elements_old(i_new, 1)
         if (any(coords .gt. this%table_dims)) then
            this%elements(i, 1) = 0
         end if
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

      if (allocated(this%table_dims)) deallocate (this%table_dims)
      if (allocated(this%table_dims_padded)) deallocate (this%table_dims_padded)
      if (allocated(this%partition_dims)) deallocate (this%partition_dims)
      if (allocated(this%elements)) deallocate (this%elements)

   end subroutine table_deallocate

   !> Converts flat index n entry (i_1, i_2, ..., i_N) in coordinate indexing
   !! using the dimensions given in part_dims
   !!
   !! @param this table object to which flat2coords belongs
   !! @param flat the flat index to be returned in coordinate index
   !! @param part_dims table dimensions for coordinate indexing
   !! @result coords the coordinates of the entry at flat index flat
   function flat2coords(this, flat, part_dims) result(coords)
      class(table), intent(inout) :: this
      integer :: flat, flat_cpy, part
      integer, dimension(size(this%partition_dims)) :: coords, part_dims
      integer :: k, div, N

      N = size(part_dims)
      div = product(part_dims(2:N))
      flat_cpy = flat

      do k = 1, N
         coords(k) = ceiling(real(flat_cpy)/real(div))
         if (mod(flat_cpy, div) .ne. 0) then
            flat_cpy = mod(flat_cpy, div)
         else
            flat_cpy = div
         end if
         if (k .lt. N) div = div/part_dims(k + 1)
      end do

   end function flat2coords

   !> Return coordinates from a global index on the object padded table dimensions.
   !!
   !! @param this table object
   !! @param ind the index whose coordinates are to be found
   !! @result coords global coordinates on the object padded table dimensions
   function get_global_coords(this, ind) result(coords)
      class(table), intent(inout) :: this
      integer :: ind
      integer, dimension(size(this%partition_dims)) :: coords

      coords = this%flat2coords(ind, this%table_dims_padded)

   end function get_global_coords

   !> Given an index, intra-partition dimensions part_dims, and inter-partition dimensions box_dims,
   !! return the coordinates of index under the partitioning scheme defined by
   !! part_dims and box_dims.
   !!
   !! @param this table object
   !! @param ind the partitioned index
   !! @param part_dims the intra-partition dimensions of the partitioning scheme
   !! @param box_dims the inter-partition dimensions of the partitioning scheme
   !! @param coords_part the array to return the intra-partition dimensions in
   !! @param coords_box the array to return the inter-partition dimensions in
   subroutine get_local_coords(this, ind, part_dims, box_dims, coords_part, coords_box)
      class(table), intent(inout) :: this
      integer :: ind, ind_local, N, k, box
      integer, dimension(size(this%partition_dims)) :: coords_part, coords_box
      integer, dimension(size(this%partition_dims)) :: part_dims, box_dims

      N = size(this%partition_dims)

      box = ceiling(real(ind)/product(part_dims))

      !! Compute inter-partition contribution to index
      coords_box = this%flat2coords(box, box_dims)

      !! Localized intra-partition index
      if (mod(ind, product(part_dims)) .ne. 0) then
         ind_local = mod(ind, product(part_dims))
      else
         ind_local = product(part_dims)
      end if

      !! Compute intra-partition contribution to index
      coords_part = this%flat2coords(ind_local, part_dims)

   end subroutine get_local_coords

   !> Converts entry (i_1, i_2, ..., i_N) in coordinate indexing to flat index n
   !! in global order, according to the given partition size part_dims.
   !!
   !! @param this table object to which coords2flat belongs
   !! @param coords the coordinates of the entry to be returned in flat index
   !! @param part_dims the partition dimensions on which to find the flat index.
   !! @result flat the flat index of entry located at coordinates coords
   function coords2flat(this, coords, part_dims) result(flat)
      class(table), intent(inout) :: this
      integer, dimension(size(this%partition_dims)) :: coords, part_dims
      integer :: flat, k, N, div

      N = size(this%partition_dims)
      div = product(part_dims)
      flat = coords(N)

      do k = 1, N - 1
         div = div/part_dims(k)
         flat = flat + (coords(k) - 1)*div
      end do

   end function coords2flat

   !> Return the global flat index, given coordinates and control variable space dimensions.
   !!
   !! @param this table object
   !! @param coords the coordinates
   !! @param table_dims_cvar the dimensions of the control variable space
   !! @result ind global flat index
   function get_index_global_coords(this, coords, table_dims_cvar) result(ind)
      class(table), intent(inout) :: this
      integer :: ind
      integer, dimension(size(this%partition_dims)) :: coords, table_dims_cvar

      ind = this%coords2flat(coords, table_dims_cvar)

   end function get_index_global_coords

   !> Return partitioned flat index, given coordinates and an associated two-level partitioning scheme
   !! described by the intra-partition dimensions part_dims and inter-partition dimensions box_dims.
   !!
   !! @param this table object
   !! @param coords_part the intra-partition term of the coordinates
   !! @param coords_box the inter-partition term of the coordinates
   !! @param part_dims the intra-partition dimensions
   !! @param box_dims the inter-partition box dimensions
   !! @result ind the partitioned index
   function get_index_local_coords(this, coords_part, coords_box, part_dims, box_dims) result(ind)
      class(table), intent(inout) :: this
      integer :: ind, ind_local, N, k, part
      integer, dimension(size(this%partition_dims)) :: coords_part, coords_box, part_dims, box_dims

      N = size(part_dims)

      do k = 1, N
         if (coords_part(k) .gt. part_dims(k)) then
            print '("get_index_local_coords ERROR: coords_part larger than part_dims at dimension ", I0)', k
         end if
         if (coords_box(k) .gt. box_dims(k)) then
            print '("get_index_local_coords ERROR: coords_box larger than box_dims at dimension ", I0)', k
         end if
         if (mod(part_dims(k), box_dims(k)) .ne. 0) then
            print '("get_index_local_coords ERROR: box_dims not divisible by part_dims at dimension ", I0)', k
         end if
      end do

      ind = this%coords2flat(coords_part, part_dims) + &
            (this%coords2flat(coords_box, box_dims) - 1)*(product(part_dims))

   end function get_index_local_coords

end module disttab_table
