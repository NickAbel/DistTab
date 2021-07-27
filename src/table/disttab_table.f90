module disttab_table

   implicit none

   private
   public :: table

   type :: table

      real, allocatable, dimension(:, :) :: elems
      integer, allocatable, dimension(:) :: table_dims
      integer, allocatable, dimension(:) :: table_dims_padded
      integer, allocatable, dimension(:) :: part_dims

      integer :: table_dims_flat
      integer :: table_dims_padded_flat
      integer :: part_dims_flat
      integer :: table_dim_svar

   contains
      ! Read in a lookup table from file
      procedure, public, pass(this)  :: read_in
      !! Pad out table for new partitioning scheme
      !procedure, public, pass(this)  :: reshape_table
      ! Map table to partition-major order
      procedure, public, pass(this)  :: partition_remap
      ! Print the table elements
      procedure, public, pass(this)  :: print_elems
      ! Deallocate table member variables (this is not the dtor!)
      procedure, public, pass(this)  :: deallocate_table

      !! Given integer coordinates in CV space, return corresponding SV value
      procedure, public, pass(this)  :: index_to_value
      procedure, public, pass(this)  :: local_coord_to_value
      procedure, public, pass(this)  :: global_coord_to_value

      !! Given integer coordinates in CV space, return corresponding SV value cloud
      procedure, public, pass(this)  :: index_to_value_cloud
      procedure, public, pass(this)  :: local_coord_to_value_cloud
      procedure, public, pass(this)  :: global_coord_to_value_cloud

      ! Convert index to table coordinates
      procedure, public, pass(this) :: index_to_coord
      ! Find global coordinate of index given partitioning scheme
      procedure, public, pass(this) :: index_to_global_coord
      ! Find two-level decomposition of coordinates of index given partitioning scheme
      procedure, public, pass(this) :: index_to_local_coord

      ! Convert coordinates to index
      procedure, public, pass(this) :: coord_to_index
      ! Get index from global coordinate and partitioning scheme
      procedure, public, pass(this) :: global_coord_to_index
      ! Get index from two-level decomposition of coordinates and partitioning scheme
      procedure, public, pass(this) :: local_coord_to_index

      ! Convert global coordinate two two-level decomposition of coordinates in partitioning scheme
      procedure, public, pass(this) :: global_coord_to_local_coord
      ! Convert two-level decomposition of coordinates and partitioning scheme to global coordinate
      procedure, public, pass(this) :: local_coord_to_global_coord

      final :: table_destructor

   end type table

   interface table
      module procedure :: table_constructor
   end interface table

!   interface get_value
!     module procedure :: index_to_value, local_coord_to_value, global_coord_to_value
!   end interface get_value

contains

   function index_to_value(this, ind) result(val)
      class(table), intent(inout) :: this
      integer, intent(in) :: ind
      real, dimension(this%table_dim_svar) :: val
      print *, "index_to_value"
      val = this%elems(:, ind)
   end function index_to_value

   function local_coord_to_value(this, coord_p, coord_b) result(val)
      class(table), intent(inout) :: this
      integer, dimension(:), intent(in)  :: coord_p, coord_b
      integer :: ind, N
      integer, dimension(size(this%part_dims)) :: box_dims
      real, dimension(this%table_dim_svar) :: val
      N = size(this%part_dims)
      print *, "local_coord_to_value"
      box_dims = this%table_dims_padded(1:N)/this%part_dims
      ind = this%local_coord_to_index(coord_p, coord_b, this%part_dims, box_dims)
      val = this%elems(:, ind)
   end function local_coord_to_value

   function global_coord_to_value(this, coord) result(val)
      class(table), intent(inout) :: this
      integer, dimension(:), intent(in) :: coord
      integer :: ind, N
      integer, dimension(size(this%part_dims)) :: box_dims
      real, dimension(this%table_dim_svar) :: val
      N = size(this%part_dims)
      print *, "global_coord_to_value"
      box_dims = this%table_dims_padded(1:N)/this%part_dims
      ind = this%global_coord_to_index(coord, this%part_dims, box_dims)
      val = this%elems(:, ind)
   end function global_coord_to_value

   function index_to_value_cloud(this, ind) result(val_cloud)
      class(table), intent(inout) :: this
      integer, intent(in) :: ind
      integer :: i, N, ind_cpy
      integer, dimension(size(this%part_dims)) :: box_dims, coord
      real, dimension(this%table_dim_svar, size(this%part_dims) + 1) :: val_cloud
      N = size(this%part_dims)
      print *, "index_to_value_cloud"
      box_dims = this%table_dims_padded(1:N)/this%part_dims
      coord = this%index_to_global_coord(ind, this%part_dims, box_dims)
      val_cloud(:, 1) = this%elems(:, ind)
      do i = 2, size(this%part_dims) + 1
         coord(i-1) = coord(i-1) + 1
         ind_cpy = this%global_coord_to_index(coord, this%part_dims, box_dims)
         val_cloud(:, i) = this%elems(:, ind_cpy)
         coord(i-1) = coord(i-1) - 1
      end do
   end function index_to_value_cloud

   function local_coord_to_value_cloud(this, coord_p, coord_b) result(val_cloud)
      class(table), intent(inout) :: this
      integer, dimension(:), intent(in)  :: coord_p, coord_b
      integer :: ind, N, i
      integer, dimension(size(this%part_dims)) :: box_dims, coord
      real, dimension(this%table_dim_svar, size(this%part_dims) + 1) :: val_cloud
      N = size(this%part_dims)
      print *, "local_coord_to_value_cloud"
      box_dims = this%table_dims_padded(1:N)/this%part_dims
      ind = this%local_coord_to_index(coord_p, coord_b, this%part_dims, box_dims)
      val_cloud(:, 1) = this%elems(:, ind)
      coord = this%local_coord_to_global_coord(coord_p, coord_b, this%part_dims, box_dims)
      do i = 2, size(this%part_dims) + 1
         coord(i-1) = coord(i-1) + 1
         ind = this%global_coord_to_index(coord, this%part_dims, box_dims)
         val_cloud(:, i) = this%elems(:, ind)
         coord(i-1) = coord(i-1) - 1
      end do
      val_cloud(:, i) = this%elems(:, ind)
   end function local_coord_to_value_cloud

   function global_coord_to_value_cloud(this, coord) result(val_cloud)
      class(table), intent(inout) :: this
      integer, dimension(:), intent(in) :: coord
      integer :: ind, N, i
      integer, dimension(size(this%part_dims)) :: box_dims, coord_cpy
      real, dimension(this%table_dim_svar, size(this%part_dims) + 1) :: val_cloud
      N = size(this%part_dims)
      print *, "global_coord_to_value_cloud"
      box_dims = this%table_dims_padded(1:N)/this%part_dims
      ind = this%global_coord_to_index(coord, this%part_dims, box_dims)
      val_cloud(:, 1) = this%elems(:, ind)
      coord_cpy = coord
      do i = 2, size(this%part_dims) + 1
         coord_cpy(i-1) = coord_cpy(i-1) + 1
         ind = this%global_coord_to_index(coord_cpy, this%part_dims, box_dims)
         val_cloud(:, i) = this%elems(:, ind)
         coord_cpy(i-1) = coord_cpy(i-1) - 1
      end do
      val_cloud(:, i) = this%elems(:, ind)
   end function global_coord_to_value_cloud

   !> Constructor for the table object.
   !! Allocates the elements array according to the table dimensions.
   !! Initializes the member variables table_dims, table_dims_padded.
   !! Computes the parameters table_dims_flat,
   !! table_dims_padded_flat, table_dim_svar,
   !! which are all a property of the table_dims argument, necessary
   !! for partition mapping functionality.
   !!
   !! @param table_dims the length of the control variable space
   !! in each dimension, and the number of state variables
   !! @result this table object
   !! @todo if partition dimensions are a property of a lookup table which
   !! may change, the part_dims can be stored as a member variable
   !! but perhaps not initialized here. Then, the same is true of the "padded"
   !! variables which generally will change as the part_dims change.
   !! Organize the code as such. Perhaps use an optional argument for part_dims
   !! here instead.
   type(table) function table_constructor(table_dims) result(this)
      integer, dimension(:), intent(in) :: table_dims

      allocate (this%table_dims(size(table_dims)))
      allocate (this%table_dims_padded(size(table_dims)))
      allocate (this%part_dims(size(table_dims) - 1))

      this%table_dims = table_dims
      this%table_dims_padded = table_dims
      this%table_dims_flat = product(this%table_dims(1:ubound(this%table_dims, dim=1) - 1))
      this%table_dims_padded_flat = product(this%table_dims(1:ubound(this%table_dims, dim=1) - 1))
      this%table_dim_svar = this%table_dims(ubound(this%table_dims, dim=1))

      this%part_dims = this%table_dims_padded(1:ubound(this%table_dims, dim=1) - 1) ! Table initially considered to have one table-sized partition, i.e. 'unpartitioned'
      this%part_dims_flat = product(this%part_dims)

      allocate (this%elems(this%table_dim_svar, this%table_dims_padded_flat))

   end function table_constructor

   !> destructor for the table type
   !!
   !! @param this the table object to destruct
   !! @todo what is this exactly doing? is deallocate_table call necessary?
   subroutine table_destructor(this)
      type(table) :: this

      call deallocate_table(this)

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
      integer :: i

      open (1, file=file_id, action='read')
      do i = 1, this%table_dims_flat
         read (unit=1, fmt=*) this%elems(:, i)
      end do
      close (1)

   end subroutine read_in

   !> Remaps the partition from assumed Alya-format ordering to given partition
   !! ordering.
   !!
   !! @param this table object to perform partition mapping
   !! @param part_dims partition dims to use
   !! @todo assumed currently that table is in Alya-format. Add functionality
   !! to remap from one partition size to another. This can be achieved by
   !! writing and sorting the table as in test_partitioning%part_test_fill_table
   !! to recover the Alya format and then reading that sorted file back in.
   !! @todo move nasty reshaping code to another function
   subroutine partition_remap(this, part_dims, part_dims_prev)
      class(table), intent(inout) :: this
      integer, dimension(size(this%table_dims) - 1), intent(in) :: part_dims, part_dims_prev

      integer, dimension(size(this%table_dims) - 1) :: coord, coord_b, coord_p
      integer, dimension(size(this%table_dims) - 1) :: box_dims, box_dims_prev
      integer :: i, i_old, N
      real, allocatable, dimension(:, :) :: elems_old

      N = size(this%table_dims) - 1
      allocate (elems_old(this%table_dim_svar, this%table_dims_padded_flat))
      elems_old = this%elems

      box_dims_prev = this%table_dims_padded(1:N)/part_dims_prev
      !! Pad out table to maintain shape
      ! Find padded table dims
      this%table_dims_padded = this%table_dims
      do i = lbound(this%table_dims_padded, dim=1), ubound(this%table_dims_padded, dim=1) - 1
         do while (mod(this%table_dims_padded(i), part_dims(i)) .ne. 0)
            this%table_dims_padded(i) = this%table_dims_padded(i) + 1
         end do
      end do
      this%table_dims_padded_flat = &
         product(this%table_dims_padded(1:ubound(this%table_dims_padded, dim=1) - 1))

      ! Create a new padded table
      deallocate (this%elems)
      allocate (this%elems(this%table_dim_svar, this%table_dims_padded_flat))
      this%elems = 0.d0

      this%part_dims = part_dims
      box_dims = this%table_dims_padded(1:N)/this%part_dims

      do i = 1, this%table_dims_padded_flat
         call this%index_to_local_coord(i, this%part_dims, box_dims, coord_p, coord_b)
         coord = this%local_coord_to_global_coord(coord_p, coord_b, this%part_dims, box_dims)
         if (any(coord .gt. this%table_dims(1:N))) then
            this%elems(:, i) = 0
         else
            i_old = this%global_coord_to_index(coord, part_dims_prev, box_dims_prev)
            this%elems(:, i) = elems_old(:, i_old)
         end if
      end do

      deallocate (elems_old)

   end subroutine partition_remap

   !> Write elements to stdout.
   !!
   !! @param this table object whose elements are to be written
   subroutine print_elems(this)
      class(table), intent(inout) :: this
      integer :: i
      !write (*, fmt='(f9.3)') this%elems
      do i = 1, this%table_dims_flat
         write (*, fmt='(*(e16.8))') this%elems(:, i)
      end do

   end subroutine print_elems

   !> Deallocates the allocated (allocatable) member variables of the table object
   !!
   !! @param this the table object whose allocatable member variables are to be deallocated
   !! @todo This is called by the destructor, but is it necessary?
   subroutine deallocate_table(this)
      class(table), intent(inout) :: this

      if (allocated(this%table_dims)) deallocate (this%table_dims)
      if (allocated(this%table_dims_padded)) deallocate (this%table_dims_padded)
      if (allocated(this%part_dims)) deallocate (this%part_dims)
      if (allocated(this%elems)) deallocate (this%elems)

   end subroutine deallocate_table

   !> Converts flat index n entry (i_1, i_2, ..., i_N) in coordinate indexing
   !! using the dimensions given in part_dims
   !!
   !! @param this table object to which flat2coord belongs
   !! @param ind the flat index to be returned in coordinate index
   !! @param dims dimensions for coordinate indexing
   !! @result coord the coordinates of the entry at flat index flat
   function index_to_coord(this, ind, dims) result(coord)
      class(table), intent(inout) :: this
      integer :: ind, ind_cpy
      integer, dimension(size(this%part_dims)) :: coord
      integer, dimension(size(this%part_dims)), intent(in) :: dims
      integer :: k, div, N

      N = size(dims)
      div = product(dims(2:N))
      ind_cpy = ind

      outer: do k = 1, N
         coord(k) = ceiling(real(ind_cpy)/real(div))
         if (mod(ind_cpy, div) .ne. 0) then
            ind_cpy = mod(ind_cpy, div)
         else
            ind_cpy = div
            coord(k + 1:N) = dims(k + 1:N)
            exit outer
         end if
         if (k .lt. N) div = div/dims(k + 1)
      end do outer

   end function index_to_coord

   !> Return coordinates from a global index on the object padded table dimensions.
   !!
   !! @param this table object
   !! @param ind the index whose coordinates are to be found
   !! @result coord global coordinates on the object padded table dimensions
   function index_to_global_coord(this, ind, part_dims, box_dims) result(coord)
      class(table), intent(inout) :: this
      integer :: ind
      integer, dimension(size(this%part_dims)) :: coord, coord_p, coord_b
      integer, dimension(size(this%part_dims)) :: part_dims, box_dims

      call this%index_to_local_coord(ind, part_dims, box_dims, coord_p, coord_b)
      coord = this%local_coord_to_global_coord(coord_p, coord_b, part_dims, box_dims)

   end function index_to_global_coord

   !> Given an index, intra-partition dimensions part_dims, and inter-partition dimensions box_dims,
   !! return the coordinates of index under the partitioning scheme defined by
   !! part_dims and box_dims.
   !!
   !! @param this table object
   !! @param ind the partitioned index
   !! @param part_dims the inter-partition dimensions of the partitioning scheme
   !! @param box_dims the intra-partition dimensions of the partitioning scheme
   !! @param coord_p the array to return the inter-partition dimensions in
   !! @param coord_b the array to return the intra-partition dimensions in
   subroutine index_to_local_coord(this, ind, part_dims, box_dims, coord_p, coord_b)
      class(table), intent(inout) :: this
      integer :: ind, ind_p, ind_b, N
      integer, dimension(size(this%part_dims)) :: coord_p, coord_b
      integer, dimension(size(this%part_dims)) :: part_dims, box_dims

      N = size(part_dims)

      ind_p = ceiling(real(ind)/product(box_dims))

      !! Compute inter-partition contribution to index
      coord_p = this%index_to_coord(ind_p, part_dims)

      !! Localized intra-partition index
      ind_b = mod(ind, product(box_dims))
      if (ind_b .ne. 0) then
         coord_b = this%index_to_coord(ind_b, box_dims)
      else
         coord_b = box_dims
      end if

   end subroutine index_to_local_coord

   !> Converts entry (i_1, i_2, ..., i_N) in coordinate indexing to flat index n
   !! in global order, according to the given partition size part_dims.
   !!
   !! @param this table object
   !! @param coord the coordinates of the entry to be returned in flat index
   !! @param dims the dimensions on which to find the flat index.
   !! @result ind the flat index of entry located at coordinates coord
   function coord_to_index(this, coord, dims) result(ind)
      class(table), intent(inout) :: this
      integer, dimension(size(this%part_dims)) :: coord, dims
      integer :: ind, k, N, div

      N = size(dims)
      div = product(dims)
      ind = coord(N)

      do k = 1, N - 1
         div = div/dims(k)
         ind = ind + (coord(k) - 1)*div
      end do

   end function coord_to_index

   !> Return the global flat index, given coordinates and control variable space dimensions.
   !!
   !! @param this table object
   !! @param coord the coordinates
   !! @param table_dims_cvar the dimensions of the control variable space
   !! @result ind global flat index
   function global_coord_to_index(this, coord, part_dims, box_dims) result(ind)
      class(table), intent(inout) :: this
      integer :: ind
      integer, dimension(size(this%part_dims)) :: coord, coord_p, coord_b
      integer, dimension(size(this%part_dims)) :: part_dims, box_dims

      call this%global_coord_to_local_coord(coord, part_dims, box_dims, coord_p, coord_b)
      ind = this%local_coord_to_index(coord_p, coord_b, part_dims, box_dims)

   end function global_coord_to_index

   !> Return partitioned flat index, given coordinates and an associated two-level partitioning scheme
   !! described by the intra-partition dimensions part_dims and inter-partition dimensions box_dims.
   !!
   !! @param this table object
   !! @param coord_part the intra-partition term of the coordinates
   !! @param coord_box the inter-partition term of the coordinates
   !! @param part_dims the intra-partition dimensions
   !! @param box_dims the inter-partition box dimensions
   !! @result ind the partitioned index
   function local_coord_to_index(this, coord_p, coord_b, part_dims, box_dims) result(ind)
      class(table), intent(inout) :: this
      integer :: ind, N
      integer, dimension(size(this%part_dims)), intent(in) :: coord_p, coord_b, part_dims, box_dims
      N = size(part_dims)

      ind = (this%coord_to_index(coord_p, part_dims) - 1)*product(box_dims) + &
            this%coord_to_index(coord_b, box_dims)

   end function local_coord_to_index

   subroutine global_coord_to_local_coord(this, coord, part_dims, box_dims, coord_p, coord_b)
      class(table), intent(inout) :: this
      integer :: N, k
      integer, dimension(size(this%part_dims)) :: part_dims, box_dims
      integer, dimension(size(this%part_dims)) :: coord_p, coord_b, coord
      N = size(part_dims)

      do k = 1, N
         coord_p(k) = ceiling(real(coord(k))/box_dims(k))
         coord_b(k) = mod(coord(k), box_dims(k))
         if (coord_b(k) .eq. 0) coord_b(k) = box_dims(k)
      end do

   end subroutine global_coord_to_local_coord

   function local_coord_to_global_coord(this, coord_p, coord_b, part_dims, box_dims) result(coord)
      class(table), intent(inout) :: this
      integer :: N
      integer, dimension(size(this%part_dims)) :: part_dims, box_dims
      integer, dimension(size(this%part_dims)) :: coord_p, coord_b, coord
      N = size(part_dims)

      coord = (coord_p - 1)*box_dims + coord_b

   end function local_coord_to_global_coord

end module disttab_table
