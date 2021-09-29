!> This module contains the table object, which contains member functions
!! and variables pertaining to the creation, storing, spatial tiling, and
!! access to entries of the lookup table itself.
module disttab_table
  use :: mpi
  use :: kind_params

  implicit none

  private
  public :: table

  type :: table

    real(sp), allocatable, dimension(:) :: ctrl_vars
    real(sp), allocatable, dimension(:, :) :: elems
    integer(i4), allocatable, dimension(:) :: part_dims
    integer(i4), allocatable, dimension(:) :: table_dims
    integer(i4), allocatable, dimension(:) :: table_dims_padded
    integer(i4), allocatable, dimension(:) :: subtable_dims
    integer(i4), allocatable, dimension(:) :: subtable_dims_padded

    integer(i4) :: nvar
    integer(i4) :: table_dims_flat
    integer(i4) :: table_dims_padded_flat

    ! MPI RMA memory view window
    integer(i4) :: window

  contains

! Read in a lookup table from file
    procedure, public, pass(this) :: read_in

! Map table to partition-major order
    procedure, public, pass(this) :: partition_remap

! Print the table elements
    procedure, public, pass(this) :: print_elems

! Deallocate table member variables (this is not the dtor!)
    procedure, public, pass(this) :: deallocate_table

! Given integer coordinates in CV space, return corresponding SV value
    procedure, public, pass(this) :: index_to_value
    procedure, public, pass(this) :: local_coord_to_value
    procedure, public, pass(this) :: global_coord_to_value

! Given integer coordinates in CV space, return corresponding SV value cloud
    procedure, public, pass(this) :: index_to_value_cloud
    procedure, public, pass(this) :: local_coord_to_value_cloud
    procedure, public, pass(this) :: global_coord_to_value_cloud

! Given real coordinates in CV space, return SV value cloud for interpolation
    procedure, public, pass(this) :: real_to_value_cloud

! Find integerized coordinates or indices given real coordinates in the CV space
    procedure, public, pass(this) :: real_to_global_coord
    procedure, public, pass(this) :: real_to_global_coord_opt
    procedure, public, pass(this) :: real_to_global_coord_opt_preprocessor
    procedure, public, pass(this) :: real_to_index
    procedure, public, pass(this) :: real_to_value

    procedure, public, pass(this) :: gather_value_cloud

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

!> Given index, return corresponding value.
!!
!! @param this the table object
!! @param ind linear index in lookup table
!! @result val state variable values located at ind
  function index_to_value(this, ind) result(val)
    class(table), intent(in) :: this
    integer(i4), intent(in) :: ind
    integer(i4) :: target_rank, ierror
    integer(kind=mpi_address_kind) :: target_displacement
    real(sp), dimension(this % nvar) :: val

    if (ind .ge. lbound(this % elems, dim=2) .and. ind .le. ubound(this % elems, dim=2)) then
      val = this % elems(1, ind)
    else if (ind .lt. 1 .or. ind .gt. product(this % table_dims)) then
      write (*, '(A, I0, A)') "ERROR in index_to_value call: Requested index (", ind, ") is not within lookup table bounds."
    else
      target_rank = (ind - 1) / product(this % subtable_dims)
      target_displacement = merge(mod(ind, product(this % subtable_dims)), &
                              & product(this % subtable_dims), &
                              & mod(ind, product(this % subtable_dims)) .ne. 0) - 1
      call mpi_get(val, this % nvar, mpi_real, target_rank, target_displacement, this % nvar, &
                 & mpi_real, this % window, ierror)
    end if

  end function index_to_value

!> Given local coordinate decomposition, return corresponding value.
!!
!! @param this the table object
!! @param coord_p the intra-partition term of the coordinates
!! @param coord_b the inter-partition term of the coordinates
!! @result val state variable values located at coordinate given in (coord_p, coord_b)
  pure function local_coord_to_value(this, coord_p, coord_b) result(val)
    class(table), intent(in) :: this
    integer(i4), dimension(:), intent(in) :: coord_p, coord_b
    integer(i4) :: ind, N
    integer(i4), dimension(size(this % part_dims)) :: tile_dims
    real(sp), dimension(this % nvar) :: val

    N = size(this % part_dims)

    tile_dims = this % table_dims_padded(1:N) / this % part_dims
    ind = this % local_coord_to_index(coord_p, coord_b, this % part_dims, tile_dims)
    val = this % elems(:, ind)

  end function local_coord_to_value

!> Given global coordinate, return corresponding value.
!!
!! @param this the table object
!! @param coord global coordinate
!! @result val state variable values located at coordinate coord
  pure function global_coord_to_value(this, coord) result(val)
    class(table), intent(in) :: this
    integer(i4), dimension(:), intent(in) :: coord
    integer(i4) :: ind, N
    integer(i4), dimension(size(this % part_dims)) :: tile_dims
    real(sp), dimension(this % nvar) :: val

    N = size(this % part_dims)

    tile_dims = this % table_dims_padded(1:N) / this % part_dims
    ind = this % global_coord_to_index(coord, this % part_dims, tile_dims)
    val = this % elems(:, ind)

  end function global_coord_to_value

!> Given coordinates [0, 1]x[0, 1]x...x[0, 1] in the normalized state space, return
!! corresponding global coordinates.
!!
!! @param this the table object
!! @param real_val Coordinates of length N [0, 1]x...x[0, 1] in normalized state space
!! @result global_coord resultant global coordinates
  pure function real_to_global_coord(this, real_val) result(global_coord)
    class(table), intent(in) :: this
    real(sp), dimension(size(this % part_dims)), intent(in) :: real_val
    integer(i4) :: N, i, j, offset_cv
    integer(i4), dimension(size(this % part_dims)) :: global_coord

    N = size(this % part_dims)

    do i = 1, N
      offset_cv = sum(this % table_dims(:i - 1))
      j = 1
      do while ((real_val(i) .lt. this % ctrl_vars(j + offset_cv)) &
                & .or. (real_val(i) .ge. this % ctrl_vars(j + 1 + offset_cv)))
        j = j + 1
      end do
      global_coord(i) = j
    end do

  end function real_to_global_coord

!> Preprocessor for bucketed real to global coord function.
!! Populates the bucket array used.
!! The argument M must be equal to the argument M passed to the real_to_global_coord_opt
!! function.
!!
!! @param this table object
!! @param segments desired number of segments per dimension in the bucket array
!! @result buckets the bucket array
  pure function real_to_global_coord_opt_preprocessor(this, segments) result(buckets)
    class(table), intent(in) :: this
    integer(i4), dimension(size(this % part_dims)), intent(in) :: segments
    integer(i4), dimension(sum(segments)) :: buckets
    integer(i4) :: i, j, k, l, N
    real(sp) :: offset, delta

    N = size(this % part_dims)

    l = 1
    do i = 1, N
      delta = 1.0 / (segments(i))
      do j = 0, segments(i) - 1
        offset = delta * j
        k = sum(this % table_dims(:i - 1)) + 1
        do while (offset .gt. this % ctrl_vars(k))
          k = k + 1
        end do
        buckets(l) = max(sum(this % table_dims(:i - 1)) + 1, k - 1)
        l = l + 1
      end do
    end do

  end function real_to_global_coord_opt_preprocessor

!> Given coordinates [0, 1]x[0, 1]x...x[0, 1] in the normalized state space, return
!! corresponding global coordinates.
!! Uses optimized bucketed preprocessed technique.
!!
!! @param this the table object
!! @param real_val Coordinates of length N [0, 1]x...x[0, 1] in normalized state space
!! @param segments number of segments in the pre-processing array
!! @param buckets bucket array, corresponds to coarse linear table discretization
!! @result global_coord resultant global coordinates
  pure function real_to_global_coord_opt(this, real_val, segments, buckets) result(global_coord)
    class(table), intent(in) :: this
    integer(i4) :: N, i, j, offset_cv
    real(sp) :: delta
    real(sp), dimension(size(this % part_dims)), intent(in) :: real_val
    integer(i4), dimension(size(this % part_dims)), intent(in) :: segments
    integer(i4), dimension(:), intent(in) :: buckets
    integer(i4), dimension(size(this % part_dims)) :: global_coord, coord_start_indices

    N = size(this % part_dims)

    do i = 1, N
      delta = 1.0 / (segments(i))
      offset_cv = sum(this % table_dims(:i - 1))
      coord_start_indices(i) = buckets(ceiling(real_val(i) / delta) + sum(segments(:i - 1)))
      j = coord_start_indices(i)
      do while ((real_val(i) .lt. this % ctrl_vars(j)) &
        & .or. (real_val(i) .ge. this % ctrl_vars(j + 1)))
        j = j + 1
      end do
      global_coord(i) = j - offset_cv
    end do

  end function real_to_global_coord_opt

!> Given coordinates [0, 1]x[0, 1]x...x[0, 1] in the normalized state space, return
!! corresponding linear index
!!
!! @param this the table object
!! @param real_val Coordinates of length N [0, 1]x...x[0, 1] in normalized state space
!! @result resultant index
  pure function real_to_index(this, real_val) result(ind)
    class(table), intent(in) :: this
    real(sp), dimension(size(this % part_dims)), intent(in) :: real_val
    integer(i4), dimension(size(this % part_dims)) :: tile_dims, global_coord
    integer(i4) :: ind, N

    N = size(this % part_dims)

    global_coord = this % real_to_global_coord(real_val)
    tile_dims = this % table_dims_padded(1:N) / this % part_dims
    ind = this % global_coord_to_index(global_coord, this % part_dims, tile_dims)

  end function real_to_index

!> Given coordinates [0, 1]x[0, 1]x...x[0, 1] in the normalized state space, return
!! corresponding state variable values.
!!
!! @param this the table object
!! @param real_val Coordinates of length N [0, 1]x...x[0, 1] in normalized state space
!! @result val resultant SV values
  pure function real_to_value(this, real_val) result(val)
    class(table), intent(in) :: this
    real(sp), dimension(size(this % part_dims)), intent(in) :: real_val
    real(sp), dimension(this % nvar) :: val
    integer(i4) :: ind

    ind = this % real_to_index(real_val)
    val = this % elems(:, ind)

  end function real_to_value

!> Return corresponding 2**N cloud of neighbors in the state space
!! beginning at the given index.
!!
!! @param this the table object
!! @param ind linearized index of point
!! @result val_cloud cloud of values based at coordinate corresponding to
!! linearized index ind.
  function index_to_value_cloud(this, ind) result(val_cloud)
    class(table), intent(inout) :: this
    integer(i4), intent(in) :: ind
    integer(i4) :: N, j
    integer(i4), dimension(size(this % part_dims)) :: tile_dims, coord
    real(sp), dimension(this % nvar, 2**size(this % part_dims)) :: val_cloud

    N = size(this % part_dims)

    tile_dims = this % table_dims_padded(1:N) / this % part_dims
    coord = this % index_to_global_coord(ind, this % part_dims, tile_dims)

    if (any(coord + 1 .gt. this % table_dims_padded(1:N))) then
      print *, "ERROR: value cloud off table"
      print *, "coord = ", coord, "table_dims_padded = ", this % table_dims_padded(1:N)
    end if

    j = 1
    call this % gather_value_cloud(N, coord, coord + 1, val_cloud, j, tile_dims)

  end function index_to_value_cloud

!> Return corresponding 2**N cloud of neighbors in the state space
!! beginning at the point indicated by the local coordinate decomposition
!! given.
!!
!! @param this the table object
!! @param coord_p the intra-partition term of the coordinates
!! @param coord_b the inter-partition term of the coordinates
!! @result val_cloud resultant value cloud of 2**N [0, 1]x...x[0, 1]
  function local_coord_to_value_cloud(this, coord_p, coord_b) result(val_cloud)
    class(table), intent(inout) :: this
    integer(i4), dimension(:), intent(in) :: coord_p, coord_b
    integer(i4) :: ind, N, j
    integer(i4), dimension(size(this % part_dims)) :: tile_dims, coord
    real(sp), dimension(this % nvar, 2**size(this % part_dims)) :: val_cloud

    N = size(this % part_dims)

    tile_dims = this % table_dims_padded(1:N) / this % part_dims
    ind = this % local_coord_to_index(coord_p, coord_b, this % part_dims, tile_dims)
    coord = this % local_coord_to_global_coord(coord_p, coord_b, tile_dims)

    if (any(coord + 1 .gt. this % table_dims_padded(1:N))) then
      print *, "ERROR: value cloud off table"
      print *, "coord = ", coord, "table_dims_padded = ", this % table_dims_padded(1:N)
    end if

    j = 1
    call this % gather_value_cloud(N, coord, coord + 1, val_cloud, j, tile_dims)

  end function local_coord_to_value_cloud

!> Return corresponding 2**N cloud of neighbors in the state space
!! beginning at the point indicated by the global coordinate given.
!!
!! @param this the table object
!! @param coord global coordinate
!! @result val_cloud resultant value cloud of 2**N [0, 1]x...x[0, 1]
  function global_coord_to_value_cloud(this, coord) result(val_cloud)
    class(table), intent(inout) :: this
    integer(i4), dimension(:), intent(in) :: coord
    integer(i4) :: ind, N, j
    integer(i4), dimension(size(this % part_dims)) :: tile_dims, coord_cpy
    real(sp), dimension(this % nvar, 2**size(this % part_dims)) :: val_cloud

    N = size(this % part_dims)

    if (any(coord + 1 .gt. this % table_dims_padded(1:N))) then
      print *, "ERROR: value cloud off table"
      print *, "coord = ", coord, "table_dims_padded = ", this % table_dims_padded(1:N)
    end if

    tile_dims = this % table_dims_padded(1:N) / this % part_dims
    ind = this % global_coord_to_index(coord, this % part_dims, tile_dims)
    coord_cpy = coord
    j = 1
    call this % gather_value_cloud(N, coord_cpy, coord_cpy + 1, val_cloud, j, tile_dims)

  end function global_coord_to_value_cloud

!> Given coordinates [0, 1]x[0, 1]x...x[0, 1] in the normalized state space, return
!! a value cloud.
!!
!! @param this the table object
!! @param real_val Coordinates of length N [0, 1]x...x[0, 1] in normalized state space
!! @result val_cloud resultant value cloud of 2**N [0, 1]x...x[0, 1]
  function real_to_value_cloud(this, real_val) result(val_cloud)
    class(table), intent(inout) :: this
    integer(i4) :: N, j
    real(sp), dimension(size(this % part_dims)), intent(in) :: real_val
    integer(i4), dimension(size(this % part_dims)) :: coord_base, tile_dims
    real(sp), dimension(this % nvar, 2**size(this % part_dims)) :: val_cloud

    N = size(this % part_dims)

    coord_base = this % real_to_global_coord(real_val)

    tile_dims = this % table_dims_padded(1:N) / this % part_dims

    j = 1
    call this % gather_value_cloud(N, coord_base, coord_base + 1, val_cloud, j, tile_dims)

  end function real_to_value_cloud

!> Fills a 2^N-sized array with a value cloud of neighbors in state space.
!!
!!
!! @param this the table object
!! @param idx index of the array to increment
!! @param ctrs counter array
!! @param uppers when the counter array reaches this value, stop incrementing
!! @param val_cloud output 2**N value cloud
!! @param j incrementer for value cloud array
!! @param tile_dims intra-partition box dimension size in partitioning scheme
  recursive subroutine gather_value_cloud(this, idx, ctrs, uppers, val_cloud, j, tile_dims)
    class(table), intent(inout) :: this
    integer(i4), intent(in) :: idx
    integer(i4) :: j, ind, N
    integer(i4), dimension(size(this % part_dims)) :: ctrs, uppers, tile_dims
    integer(i4), dimension(size(this % part_dims)) :: ctrs_copy
    real(sp), dimension(this % nvar, 2**size(this % part_dims)) :: val_cloud

    N = size(this % part_dims)

    if (idx .eq. 1) then
      ind = this % global_coord_to_index(ctrs, this % part_dims, tile_dims)
      val_cloud(:, j) = this % elems(:, ind)
      j = j + 1
      do while (ctrs(N - idx + 1) .lt. uppers(N - idx + 1))
        ctrs(N - idx + 1) = ctrs(N - idx + 1) + 1
        ind = this % global_coord_to_index(ctrs, this % part_dims, tile_dims)
        val_cloud(:, j) = this % elems(:, ind)
        j = j + 1
      end do
    else if (idx .gt. 1) then
      do while (ctrs(N - idx + 1) .le. uppers(N - idx + 1))
        ctrs_copy = ctrs
        call this % gather_value_cloud(idx - 1, ctrs_copy, uppers, val_cloud, j, tile_dims)
        ctrs(N - idx + 1) = ctrs(N - idx + 1) + 1
      end do
    end if

  end subroutine gather_value_cloud

!> Constructor for the table object.
!! Allocates the elements array according to the table dimensions.
!! Initializes the member variables table_dims, table_dims_padded.
!! Computes the parameters table_dims_flat,
!! table_dims_padded_flat, nvar,
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
  type(table) function table_constructor(table_dims, subtable_dims) result(this)
    integer(i4), dimension(:), intent(in) :: table_dims
    integer(i4), dimension(:), intent(in), optional :: subtable_dims
    integer(i4) :: nprocs, ierror, rank, real_size
    integer(kind=mpi_address_kind) :: subtable_size

    ! Parallel
    if (present(subtable_dims)) then
      call mpi_comm_size(mpi_comm_world, nprocs, ierror)
      call mpi_type_size(mpi_real, real_size, ierror)
      call mpi_comm_rank(mpi_comm_world, rank, ierror)

      allocate (this % table_dims(size(table_dims)))
      allocate (this % table_dims_padded(size(table_dims)))
      allocate (this % part_dims(size(table_dims) - 1))
      allocate (this % subtable_dims(size(table_dims) - 1))

      this % table_dims = table_dims
      this % table_dims_padded = table_dims
      this % table_dims_flat = product(this % table_dims(1:ubound(this % table_dims, dim=1) - 1))
      this % table_dims_padded_flat = product(this % table_dims(1:ubound(this % table_dims, dim=1) - 1))
      this % nvar = this % table_dims(ubound(this % table_dims, dim=1))
      this % subtable_dims = subtable_dims
      this % subtable_dims_padded = subtable_dims
      this % part_dims = this % table_dims_padded(1:ubound(this % table_dims, dim=1) - 1)

      ! Create subtable elems with universal integer bounds
      allocate (this % elems(this % nvar, &
        & rank * product(this % subtable_dims) + 1:(rank + 1) * product(this % subtable_dims)))

      ! Create the MPI window
      subtable_size = product(this % subtable_dims) * real_size
      call mpi_win_create(this % elems, subtable_size, &
        & real_size, mpi_info_null, mpi_comm_world, this % window, ierror)
      call mpi_win_fence(0, this % window, ierror)

      ! TODO control vars in parallel
      !allocate (this % ctrl_vars(sum(this % table_dims(1:ubound(this % table_dims, dim=1) - 1))))

      ! Serial
    else

      allocate (this % table_dims(size(table_dims)))
      allocate (this % table_dims_padded(size(table_dims)))
      allocate (this % part_dims(size(table_dims) - 1))

      this % table_dims = table_dims
      this % table_dims_padded = table_dims
      this % table_dims_flat = product(this % table_dims(1:ubound(this % table_dims, dim=1) - 1))
      this % table_dims_padded_flat = product(this % table_dims(1:ubound(this % table_dims, dim=1) - 1))
      this % subtable_dims = table_dims
      this % subtable_dims_padded = table_dims
      this % nvar = this % table_dims(ubound(this % table_dims, dim=1))

! Table initially considered to have one table-sized partition, i.e. 'unpartitioned'
! Note there are other partitioning schemes which are identical to this scheme.
      this % part_dims = this % table_dims_padded(1:ubound(this % table_dims, dim=1) - 1)

      allocate (this % elems(this % nvar, this % table_dims_padded_flat))
      allocate (this % ctrl_vars(sum(this % table_dims(1:ubound(this % table_dims, dim=1) - 1))))
    end if

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
!! Assumes that the first size(this%part_dims) lines are normalized
!! control variable discretizations, and the following lines are
!! state variable values.
!! @param this table object to which read_in is a member.
!! @param file_id the lookup table filename to read in.
  subroutine read_in(this, file_id)
    class(table), intent(inout) :: this
    character(len=*), intent(in) :: file_id
    integer(i4) :: access_mode, rank, ierror, handle, N, &
    & cnt
    integer :: statu(mpi_status_size)
    integer(kind=mpi_offset_kind) :: offset

    N = size(this % part_dims)

    offset = 0
    cnt = product(this % table_dims(1:N))

    ! Sequential
    !integer(i4) :: i
    !open (1, file=file_id, action='read')
    !read (unit=1, fmt=*) this % ctrl_vars
    !print *, size(this % ctrl_vars)
    !print *, this % table_dims_flat
    !do i = 1, this % table_dims_flat
    !  read (unit=1, fmt=*) this % elems(:, i)
    !end do
    !close (1)

    call mpi_comm_rank(mpi_comm_world, rank, ierror)
    write (*, '(A,I0,A)') '[MPI process ', rank, '] read_in'

    access_mode = mpi_mode_rdonly
    call mpi_file_open(mpi_comm_world, file_id, access_mode, mpi_info_null, handle, ierror)
    if (ierror .ne. mpi_success) then
      write (*, '(A,I0,A)') '[MPI process ', rank, '] Failure in opening the file.'
      call mpi_abort(mpi_comm_world, -1, ierror)
    end if

    call mpi_file_sync(handle, ierror)
    call mpi_barrier(mpi_comm_world, ierror)

    call mpi_file_read(handle, this % elems, cnt, &
       & mpi_float, statu, ierror)

    print *, this % elems

    call mpi_file_close(handle, ierror)
    if (ierror .ne. mpi_success) then
      write (*, '(A,I0,A)') '[MPI process ', rank, '] Failure in closing the file.'
      call mpi_abort(mpi_comm_world, -1, ierror)
    end if
  end subroutine read_in

!> Remaps the partition from a given previous partition ordering to given new partition
!! ordering.
!!
!! @param this table object to perform partition mapping
!! @param part_dims partition dims to use
!! @param part_dims_prev partition dims in previous partition scheme
!! @todo move nasty reshaping code to another function
  subroutine partition_remap(this, part_dims, part_dims_prev)
    class(table), intent(inout) :: this
    integer(i4), dimension(size(this % table_dims) - 1), intent(in) :: part_dims, part_dims_prev
    integer(i4), dimension(size(this % table_dims) - 1) :: coord, coord_b, coord_p
    integer(i4), dimension(size(this % table_dims) - 1) :: tile_dims, tile_dims_prev
    integer(i4) :: i, i_old, N
    real(sp), allocatable, dimension(:, :) :: elems_old

    N = size(this % table_dims) - 1

    allocate (elems_old(this % nvar, this % table_dims_padded_flat))
    elems_old = this % elems

    tile_dims_prev = this % table_dims_padded(1:N) / part_dims_prev
! Pad out table to maintain shape
! Find padded table dims
    this % table_dims_padded = this % table_dims
    do i = lbound(this % table_dims_padded, dim=1), ubound(this % table_dims_padded, dim=1) - 1
      do while (mod(this % table_dims_padded(i), part_dims(i)) .ne. 0)
        this % table_dims_padded(i) = this % table_dims_padded(i) + 1
      end do
    end do
    this % table_dims_padded_flat = &
      product(this % table_dims_padded(1:ubound(this % table_dims_padded, dim=1) - 1))

! Create a new padded table
    deallocate (this % elems)
    allocate (this % elems(this % nvar, this % table_dims_padded_flat))
    this % elems = 0.d0

    this % part_dims = part_dims
    tile_dims = this % table_dims_padded(1:N) / this % part_dims

    do i = 1, this % table_dims_padded_flat
      call this % index_to_local_coord(i, this % part_dims, tile_dims, coord_p, coord_b)
      coord = this % local_coord_to_global_coord(coord_p, coord_b, tile_dims)
      if (any(coord .gt. this % table_dims(1:N))) then
        this % elems(:, i) = 0
      else
        i_old = this % global_coord_to_index(coord, part_dims_prev, tile_dims_prev)
        this % elems(:, i) = elems_old(:, i_old)
      end if
    end do

    deallocate (elems_old)

  end subroutine partition_remap

!> Write elements to stdout.
!!
!! @param this table object whose elements are to be written
  subroutine print_elems(this)
    class(table), intent(inout) :: this
    integer(i4) :: i
!write (*, fmt='(f9.3)') this%elems
    do i = 1, this % table_dims_flat
      write (*, fmt='(*(e16.8))') this % elems(:, i)
    end do

  end subroutine print_elems

!> Deallocates the allocated (allocatable) member variables of the table object
!!
!! @param this the table object whose allocatable member variables are to be deallocated
!! @todo This is called by the destructor, but is it necessary?
  subroutine deallocate_table(this)
    class(table), intent(inout) :: this

    if (allocated(this % table_dims)) deallocate (this % table_dims)
    if (allocated(this % table_dims_padded)) deallocate (this % table_dims_padded)
    if (allocated(this % part_dims)) deallocate (this % part_dims)
    if (allocated(this % elems)) deallocate (this % elems)
    if (allocated(this % ctrl_vars)) deallocate (this % ctrl_vars)

  end subroutine deallocate_table

!> Converts flat index n entry (i_1, i_2, ..., i_N) in coordinate indexing
!! using the dimensions given in part_dims
!!
!! @param this table object to which flat2coord belongs
!! @param ind the flat index to be returned in coordinate index
!! @param dims dimensions for coordinate indexing
!! @result coord the coordinates of the entry at flat index flat
  pure function index_to_coord(this, ind, dims) result(coord)
    class(table), intent(in) :: this
    integer(i4), intent(in) :: ind
    integer(i4), dimension(size(this % part_dims)), intent(in) :: dims
    integer(i4), dimension(size(this % part_dims)) :: coord
    integer(i4) :: k, div, N, ind_cpy

    N = size(this % part_dims)

    div = product(dims(2:N))

    ind_cpy = ind

    outer: do k = 1, N
      coord(k) = ceiling(real(ind_cpy, sp) / real(div, sp))
      if (mod(ind_cpy, div) .ne. 0) then
        ind_cpy = mod(ind_cpy, div)
      else
        ind_cpy = div
        coord(k + 1:N) = dims(k + 1:N)
        exit outer
      end if
      if (k .lt. N) div = div / dims(k + 1)
    end do outer

  end function index_to_coord

!> Return coordinates from a global index on the object padded table dimensions.
!!
!! @param this table object
!! @param ind the index whose coordinates are to be found
!! @param part_dims the inter-partition dimensions of the partitioning scheme
!! @param tile_dims the intra-partition dimensions of the partitioning scheme
!! @result coord global coordinates on the object padded table dimensions
  pure function index_to_global_coord(this, ind, part_dims, tile_dims) result(coord)
    class(table), intent(in) :: this
    integer(i4), intent(in) :: ind
    integer(i4), dimension(size(this % part_dims)), intent(in) :: part_dims, tile_dims
    integer(i4) :: ind_p, ind_b
    integer(i4), dimension(size(this % part_dims)) :: coord, coord_p, coord_b

    ind_p = ceiling(real(ind, sp) / product(tile_dims))

! Compute inter-partition contribution to index
    coord_p = this % index_to_coord(ind_p, part_dims)

! Localized intra-partition index
    ind_b = mod(ind, product(tile_dims))
    if (ind_b .ne. 0) then
      coord_b = this % index_to_coord(ind_b, tile_dims)
    else
      coord_b = tile_dims
    end if

    coord = this % local_coord_to_global_coord(coord_p, coord_b, tile_dims)

  end function index_to_global_coord

!> Given an index, intra-partition dimensions part_dims, and inter-partition dimensions tile_dims,
!! return the coordinates of index under the partitioning scheme defined by
!! part_dims and tile_dims.
!!
!! @param this table object
!! @param ind the partitioned index
!! @param part_dims the inter-partition dimensions of the partitioning scheme
!! @param tile_dims the intra-partition dimensions of the partitioning scheme
!! @param coord_p the array to return the inter-partition dimensions in
!! @param coord_b the array to return the intra-partition dimensions in
  subroutine index_to_local_coord(this, ind, part_dims, tile_dims, coord_p, coord_b)
    class(table), intent(inout) :: this
    integer(i4) :: ind, ind_p, ind_b
    integer(i4), dimension(size(this % part_dims)) :: coord_p, coord_b
    integer(i4), dimension(size(this % part_dims)) :: part_dims, tile_dims

    ind_p = ceiling(real(ind, sp) / product(tile_dims))

! Compute inter-partition contribution to index
    coord_p = this % index_to_coord(ind_p, part_dims)

! Localized intra-partition index
    ind_b = mod(ind, product(tile_dims))
    if (ind_b .ne. 0) then
      coord_b = this % index_to_coord(ind_b, tile_dims)
    else
      coord_b = tile_dims
    end if

  end subroutine index_to_local_coord

!> Converts entry (i_1, i_2, ..., i_N) in coordinate indexing to flat index n
!! in global order, according to the given partition size part_dims.
!!
!! @param this table object
!! @param coord the coordinates of the entry to be returned in flat index
!! @param dims the dimensions on which to find the flat index.
!! @result ind the flat index of entry located at coordinates coord
  pure function coord_to_index(this, coord, dims) result(ind)
    class(table), intent(in) :: this
    integer(i4), dimension(size(this % part_dims)), intent(in) :: coord, dims
    integer(i4) :: ind, k, N, div

    N = size(dims)

    div = product(dims)
    ind = coord(N)

    do k = 1, N - 1
      div = div / dims(k)
      ind = ind + (coord(k) - 1) * div
    end do

  end function coord_to_index

!> Return the global flat index, given coordinates and control variable space dimensions.
!!
!! @param this table object
!! @param coord the coordinates
!! @param part_dims the intra-partition dimensions
!! @param tile_dims the inter-partition box dimensions
!! @result ind global flat index
  pure function global_coord_to_index(this, coord, part_dims, tile_dims) result(ind)
    class(table), intent(in) :: this
    integer(i4) :: ind, N, k
    integer(i4), dimension(size(this % part_dims)) :: coord_p, coord_b
    integer(i4), dimension(size(this % part_dims)), intent(in) :: coord, part_dims, tile_dims

    N = size(this % part_dims)

    do k = 1, N
      coord_p(k) = ceiling(real(coord(k), sp) / tile_dims(k))
      coord_b(k) = mod(coord(k), tile_dims(k))
      if (coord_b(k) .eq. 0) coord_b(k) = tile_dims(k)
    end do
    ind = this % local_coord_to_index(coord_p, coord_b, part_dims, tile_dims)

  end function global_coord_to_index

!> Return partitioned flat index, given coordinates and an associated two-level partitioning scheme
!! described by the intra-partition dimensions part_dims and inter-partition dimensions tile_dims.
!!
!! @param this table object
!! @param coord_p the intra-partition term of the coordinates
!! @param coord_b the inter-partition term of the coordinates
!! @param part_dims the intra-partition dimensions
!! @param tile_dims the inter-partition box dimensions
!! @result ind the partitioned index
  pure function local_coord_to_index(this, coord_p, coord_b, part_dims, tile_dims) result(ind)
    class(table), intent(in) :: this
    integer(i4), dimension(size(this % part_dims)), intent(in) :: coord_p, coord_b, part_dims, tile_dims
    integer(i4) :: ind

    ind = (this % coord_to_index(coord_p, part_dims) - 1) * product(tile_dims) + &
          this % coord_to_index(coord_b, tile_dims)

  end function local_coord_to_index

!> Decompose global coordinate to local coordinate pair in given partition and box dimension
!! scheme.
!!
!! @param this table object
!! @param coord input global coordinate
!! @param part_dims intra-partition dimension size
!! @param tile_dims inter-partition box dimension size
!! @param coord_p output intra-partition coordinate term
!! @param coord_b output inter-partition box coordinate term
  subroutine global_coord_to_local_coord(this, coord, part_dims, tile_dims, coord_p, coord_b)
    class(table), intent(in) :: this
    integer(i4) :: N, k
    integer(i4), dimension(size(this % part_dims)) :: part_dims, tile_dims
    integer(i4), dimension(size(this % part_dims)) :: coord_p, coord_b, coord

    N = size(part_dims)

    do k = 1, N
      coord_p(k) = ceiling(real(coord(k), sp) / tile_dims(k))
      coord_b(k) = mod(coord(k), tile_dims(k))
      if (coord_b(k) .eq. 0) coord_b(k) = tile_dims(k)
    end do

  end subroutine global_coord_to_local_coord

!> Convert a local coordinates decomposition pair to global coordinates.
!!
!! @param this table object
!! @param coord_p intra-partition coordinate
!! @param coord_b inter-partition box coordinate
!! @param tile_dims size of inter-partition box
!! @result coord global coordinate
  pure function local_coord_to_global_coord(this, coord_p, coord_b, tile_dims) result(coord)
    class(table), intent(in) :: this
    integer(i4), dimension(size(this % part_dims)), intent(in) :: tile_dims, coord_p, coord_b
    integer(i4), dimension(size(this % part_dims)) :: coord

    coord = (coord_p - 1) * tile_dims + coord_b

  end function local_coord_to_global_coord

end module disttab_table
