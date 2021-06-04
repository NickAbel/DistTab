program get_ex
  use mpi
  implicit none
  integer :: ierror, rank, nprocs
  integer, parameter, dimension(3) :: n = (/ 4,4,1 /)
  integer, dimension(3) :: npad
  integer :: ncv_flat, n_flat, ncv_padded, n_padded
  integer, parameter, dimension(2) :: q = (/ 2,2 /)
  integer, dimension(2) :: coords
  integer :: q_flat
  double precision, dimension(:), allocatable :: partition
  integer, dimension(:), allocatable :: nloop_counters, nloop_uppers
  integer :: x(2), y ! Coordinate index (x(1),x(2)), flat index y
  integer :: a, b, c, d, i, j, k, l, m, bmi
  double precision :: z
  double precision, dimension(:,:), allocatable :: lookup_table, lookup_table_flat
  integer :: window
  integer :: integer_size, dbl_size
  integer(kind=mpi_address_kind) :: lookup_table_size
  integer(kind=mpi_address_kind) :: target_displacement

  ! Setup MPI communicator, retrieve MPI type sizes
  call mpi_init(ierror)
  call mpi_type_size(mpi_integer, integer_size, ierror)
  call mpi_type_size(mpi_double, dbl_size, ierror)
  call mpi_comm_rank(mpi_comm_world, rank, ierror)
  call mpi_comm_size(mpi_comm_world, nprocs, ierror)

  ncv_flat = product(n(1:ubound(n,dim=1)-1))
  n_flat = product(n)
  q_flat = product(q)

  ! Pad out table to maintain shape
  npad = n
  do i = lbound(npad, dim=1), ubound(npad, dim=1) - 1
    do while (mod(npad(i),q(i)) .ne. 0)
      npad(i) = npad(i) + 1
    enddo
  enddo

  ncv_padded = product(npad(1:ubound(n,dim=1)-1))
  n_padded = product(npad)

  ! Allocate lookup table
  allocate(lookup_table(((ncv_padded*rank) + 1):(ncv_padded*(rank+1)), n(ubound(n,dim=1))))
  allocate(lookup_table_flat(((ncv_padded*rank) + 1):(ncv_padded*(rank+1)), n(ubound(n,dim=1))))

  lookup_table = 0.d0
  lookup_table_flat = 0.d0

  ! Open lookup table file
  !lookup_table_filename = 'tables/2d.dat'
  !open(unit = 46, file = lookup_table_filename, status='old', action='read')

  ! Compute size of lookup table in bytes for window creation
  lookup_table_size = n_padded*dbl_size

  ! Populate lookup table with self-explanatory entries
  i=1+ncv_flat*rank
  do j = 1,n(3)
    do k = 1, n(2)
      do l = 1, n(1)
        lookup_table_flat(coords2flat((/l,k/)) + ncv_padded*rank,j) = i
        i=i+1
      enddo
    enddo
  enddo

  ! Create memory window covering whole sub-table on each rank
  call mpi_win_create(lookup_table, lookup_table_size, dbl_size, mpi_info_null, mpi_comm_world, window, ierror)
  call mpi_win_fence(0, window, ierror)

  call block_major_order()
  deallocate(lookup_table_flat)

  !Pick a random number in the global lookup table range
  call random_number(z)
  y = ceiling(z*ncv_flat)
  x = flat2coords(y)
  coords = find_local_partition(x)
  ! Print information about the point requested x
  print *, 'Flat entry y = ', y, ' has coordinates: ', x, ' at ', coords

  y = coords2flat(x)

  print *, 'Base of partition ', coords, ' is ', bm_partition_base(find_local_partition(x))

  ! Dynamically allocate one partition
  allocate(partition(1:q_flat))

  target_displacement = npad(2)/q(2)*(coords(1)-1)*q_flat + (coords(2)-1)*q_flat + 1

  print *, "RANK = ", rank, " DISPLACEMENT = ", target_displacement

  if (find_block_rank_flat(y) .ne. rank) then
  !  print *, 'Fetch block of ', q_flat, ' from rank ', find_block_rank_flat(y), ' to rank ', rank
  !  print *, 'The block containing ', y, ' begins with entry ', bm_partition_base(find_local_partition(x))
  !  print *, 'Which is located at index ', target_displacement
  !  call mpi_get(partition, q_flat, mpi_double, find_block_rank_flat(y), &
  !    target_displacement, q_flat, mpi_double, window, ierror)
  !  print *, partition
  endif

  ! Clean up, exit
  !close(46)
  deallocate(partition)
  deallocate(lookup_table)
  call mpi_win_free(window, ierror)
  call mpi_finalize(ierror)

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! A bunch of helper functions containing ugly modular arithmetic. !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !recursive function nloop_increment(idx,ctrs,uppers) result(ctrs_incr)
  !  integer, intent(in) :: idx
  !  integer, dimension(ubound(nloop_counters)), intent(in) :: ctrs, uppers
  !  integer, dimension(ubound(nloop_counters)) :: ctrs_incr
  !  print *, ctrs
  !  if (idx .eq. 1) then
  !    do while (ctrs(idx) .lt. uppers(idx))
  !      ctrs(idx) = ctrs(idx) + 1
  !    enddo
  !  else if (idx .gt. 1) then
  !    do while (ctrs(idx) .le. uppers(idx)) 
  !      call nloop_increment(idx - 1)
  !      ctrs(idx) = ctrs(idx) + 1
  !    enddo
  !  endif
  !  ctrs_incr = ctrs
  !end function


  ! Reorder the chemistry table into block-major organization
  ! so that each block is contiguous in memory for MPI GET
  subroutine block_major_order()
    integer :: offset, i, j, k, l
    integer, dimension(4) :: blk_bounds, base
    offset = 1 + ncv_padded*rank
    do i = 1 + npad(2)*rank, npad(2)*(rank + 1), q(2)
      do j = 1, npad(1), q(1) 
        blk_bounds = find_partition_bounds_coord((/j,i/))
        do k = blk_bounds(2), blk_bounds(4)
          do l = blk_bounds(1), blk_bounds(3)
            if (rank .eq. 0) then
            endif
            lookup_table(offset,1) = lookup_table_flat(coords2flat((/l,k/)),1)
            offset = offset + 1
          enddo
        enddo
      enddo
    enddo
  end subroutine

  ! Find the base of the partition in block-major organization
  function bm_partition_base(partition) result(base)
    integer, dimension(2), intent(in) :: partition
    integer :: base
    base = npad(2)*q(1)*(partition(1)-1) + q(2)*(partition(2)-1) + 1
  end function

  ! Return upper and lower bounds of rank in global flat indices
  ! Assuming the tables on each rank are NOT indexed starting at 1
  !function find_rank_bounds(bloc) result(bounds)
  !  integer, intent(in)  :: bloc
  !  integer, dimension(2) :: bounds
  !  bounds(1) = lbound(lookup_table,dim=1)!+ncv_flat*rank
  !  bounds(2) = ubound(lookup_table,dim=1)!+ncv_flat*rank
  !end function

  ! Return rank where a flat index f is located
  ! TODO Using the max() function is a hack-y fix for the problem
  ! ie when n1*n2 = 250, 500/250 = 2 but entry 500 is on table 1.
  ! This might not work in general cases
  function find_block_rank_flat(f) result(bloc)
    integer, intent(in)  :: f
    integer :: bloc
    bloc = max(0,(f-1)/(ncv_flat))
  end function

  ! Return rank where a coordinate index (c(1),c(2)) is located
  !function find_block_rank_coord(c) result(bloc)
  !  integer, dimension(2), intent(in)  :: c
  !  integer :: bloc
  !  bloc = coords2flat(c)/ncv_flat
  !end function

  ! Return local partition containing coordinate index (c(1),c(2))
  function find_local_partition(c) result(partition_index)
    integer, dimension(2) :: c
    integer, dimension(2) :: partition_index
    partition_index(1) = ceiling(real(c(1))/real(q(1)))
    partition_index(2) = ceiling(real(c(2))/real(q(2)))
  end function

  ! Find a global coordinate (c(1),c(2))'s local (on partition) coordinates
  !function find_local_partition_entry(c) result(local_entry)
  !  integer, dimension(2) :: c
  !  integer, dimension(2) :: local_entry
  !  local_entry(1) = merge(mod(c(1),q(1)),q(1),mod(c(1),q(1)) .ne. 0)
  !  local_entry(2) = merge(mod(c(2),q(2)),q(2),mod(c(2),q(2)) .ne. 0)
  !end function

  ! Return coordinate bounds of flat index f's partition
  !function find_partition_bounds_flat(f) result(bounds)
  !  integer, intent(in)  :: f
  !  integer, dimension(2) :: c
  !  integer, dimension(4) :: bounds
  !  c = flat2coords(f)
  !  bounds = find_partition_bounds_coord(c)
  !end function

  !! Return coordinate bounds of coordinate index (c(1),c(2))'s partition
  function find_partition_bounds_coord(c) result(bounds)
    integer, dimension(2), intent(in)  :: c
    integer, dimension(4) :: bounds
    ! Lower bound direction 1
    bounds(1) = c(1)-merge(mod(c(1),q(1))-1,q(1)-1,mod(c(1),q(1)) .ne. 0)
    ! Lower bound direction 2
    bounds(2) = c(2)-merge(mod(c(2),q(2))-1,q(2)-1,mod(c(2),q(2)) .ne. 0)
    ! Upper bound direction 1
    bounds(3) = c(1)+merge(q(1)-mod(c(1),q(1)),0,mod(c(1),q(1)) .ne. 0)
    ! Upper bound direction 2
    bounds(4) = c(2)+merge(q(2)-mod(c(2),q(2)),0,mod(c(2),q(2)) .ne. 0)
  end function

  ! Given coordinate indices, return flat index
  function coords2flat(c) result(f)
    integer, dimension(2), intent(in)  :: c
    integer :: f
    f = (c(2)-1)*npad(1) + c(1)
  end function

  ! Given flat index, return coordinate indices
  function flat2coords(f) result(c)
    integer, intent(in)   :: f
    integer, dimension(2) :: c
    c(1) = ceiling(real(f)/real(n(2)))
    c(2) = merge(npad(2), mod(f,npad(2)), mod(f,npad(2)) .eq. 0)
  end function

end program get_ex

