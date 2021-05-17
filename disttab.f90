program get_ex
  use mpi
  implicit none
  character (len=255) :: lookup_table_filename
  integer :: ierror, rank, nprocs
  integer, parameter, dimension(3) :: n = (/ 25, 10, 1 /)
  integer, parameter, dimension(2) :: q = (/ 5, 2 /), p = (/ 3, 3 /)
  integer, dimension(3) :: bd
  integer, dimension(:,:), allocatable :: partition
  integer :: x(2), y ! Coordinate index (x(1),x(2)), flat index y
  integer :: block_bounds(4)
  integer :: i, j, k
  double precision :: value_fetched, z
  double precision, dimension(:,:), allocatable :: lookup_table
  integer :: window
  integer :: integer_size, dbl_size
  integer(kind=mpi_address_kind) :: win_size, lookup_table_size
  integer(kind=mpi_address_kind) :: target_displacement

  ! Setup MPI communicator, retrieve MPI type sizes
  call mpi_init(ierror)
  call mpi_type_size(mpi_integer, integer_size, ierror)
  call mpi_type_size(mpi_double, dbl_size, ierror)
  call mpi_comm_rank(mpi_comm_world, rank, ierror)
  call mpi_comm_size(mpi_comm_world, nprocs, ierror)

  ! Allocate lookup table
  allocate(lookup_table(((n(1)*n(2)*rank) + 1):(n(1)*n(2)*(rank+1)), n(3)))

  ! Open lookup table file
  lookup_table_filename = 'tables/2d.dat'
  open(unit = 46, file = lookup_table_filename, status='old', action='read')

  ! Compute size of lookup table in bytes for window creation
  lookup_table_size = n(1)*n(2)*n(3)*dbl_size

  ! Insist on some comm size
  if (nprocs .ne. 2) then
    print *, "use 2 mpi processes for now, we detect ", nprocs
    call mpi_abort(mpi_comm_world, -1, ierror)
  endif

  ! Block division per-rank
  do i = 1, 3
    bd(i) = n(i)/nprocs
  enddo

  ! Populate lookup table with self-explanatory entries
  do j = 1,n(3)
    do k = lbound(lookup_table,dim=1), ubound(lookup_table,dim=1)
      !y = coords2flat((/ k,j /))
      !x = flat2coords(y)
      !!! Debug output for indexing conversion
      !if (rank .eq. 1) then
      !  print *, k-x(1), j-x(2), n(1)*n(2)*rank+y
      !endif 
      lookup_table(k,j) = k + 0.d0
    enddo
  enddo
  !print *, lookup_table

  ! Create memory window covering whole sub-table on each rank
  call mpi_win_create(lookup_table, lookup_table_size, dbl_size, mpi_info_null, mpi_comm_world, window, ierror)
  call mpi_win_fence(0, window, ierror)

  ! Pick a random number in the global lookup table range
  call random_number(z)
  y = ceiling(z*n(1)*n(2)*nprocs)
  x = flat2coords(y)

  ! Print information about the point requested x
  print *, 'Flat entry y = ', y, ' has coordinates: ', x

  y = coords2flat(x)

  print *,  'Verify: (', lookup_table(y,1), ')', &
    ' in partition block: ', find_local_partition(x), &
    'with local partition entry ', find_local_partition_entry(x)

  ! Print information about x's partition block
  block_bounds = find_partition_bounds_coord(x)

  print *, "lookup table on rank ", rank, " has range ", find_rank_bounds(rank)

  ! Dynamically allocate one partition
  allocate(partition(block_bounds(1):block_bounds(2), block_bounds(3):block_bounds(4)))

  ! Determine if the partition entries are on-table. MPI_GET any entries that are not.
  do i = block_bounds(1), block_bounds(2)
    do j = block_bounds(3), block_bounds(4)
      y = coords2flat((/ i,j /))
      partition(i,j) = merge(lookup_table(y,1), -1.0d0, find_block_rank_flat(y) .eq. rank)
      print *, i,j, " corresponds to x = ", y, merge("     (On table at rank ", " (OUT OF TABLE ON RANK ", &
        find_block_rank_flat(y) .eq. rank), rank, ") "
      if (partition(i,j) .eq. -1.0) then
        print *, "MPI GETting entry ", y, " from table ", find_block_rank_flat(y)
        target_displacement = mod(y-1,250)
        call mpi_get(value_fetched, 1, mpi_double, find_block_rank_flat(y), &
          target_displacement, 1, mpi_double, window, ierror)
        partition(i,j) = value_fetched
      endif
    enddo
  enddo
  print *, partition

  ! Clean up, exit
  close(46)
  deallocate(partition)
  call mpi_win_free(window, ierror)
  call mpi_finalize(ierror)

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! A bunch of helper functions containing ugly modular arithmetic. !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!
  ! Return upper and lower bounds of rank in global flat indices
  ! Assuming the tables on each rank are NOT indexed starting at 1
  function find_rank_bounds(bloc) result(bounds)
    integer, intent(in)  :: bloc
    integer, dimension(2) :: bounds
    bounds(1) = lbound(lookup_table,dim=1)!+n(1)*n(2)*rank
    bounds(2) = ubound(lookup_table,dim=1)!+n(1)*n(2)*rank
  end function

  !!!!
  ! Return rank where a flat index f is located
  ! TODO Using the max() function is a hack-y fix for the problem
  ! ie when n1*n2 = 250, 500/250 = 2 but entry 500 is on table 1.
  ! This might not work in general cases
  function find_block_rank_flat(f) result(bloc)
    integer, intent(in)  :: f
    integer :: bloc
    bloc = max(0,(f-1)/(n(1)*n(2)))
  end function

  ! Return rank where a coordinate index (c(1),c(2)) is located
  function find_block_rank_coord(c) result(bloc)
    integer, dimension(2), intent(in)  :: c
    integer :: bloc
    bloc = coords2flat(c)/(n(1)*n(2))
  end function

  ! Return local partition containing coordinate index (c(1),c(2))
  function find_local_partition(c) result(partition_index)
    integer, dimension(2) :: c
    integer, dimension(2) :: partition_index
    partition_index(1) = ceiling(real(c(1))/real(q(1)))
    partition_index(2) = ceiling(real(c(2))/real(q(2)))
  end function

  ! Find a global coordinate (c(1),c(2))'s local (on partition) coordinates
  function find_local_partition_entry(c) result(local_entry)
    integer, dimension(2) :: c
    integer, dimension(2) :: local_entry
    local_entry(1) = merge(mod(c(1),q(1)),q(1),mod(c(1),q(1)) .ne. 0)
    local_entry(2) = merge(mod(c(2),q(2)),q(2),mod(c(2),q(2)) .ne. 0)
  end function


  ! Return coordinate bounds of flat index f's partition
  function find_partition_bounds_flat(f) result(bounds)
    integer, intent(in)  :: f
    integer, dimension(2) :: c
    integer, dimension(4) :: bounds
    c = flat2coords(f)
    bounds = find_partition_bounds_coord(c)
  end function

  ! Return coordinate bounds of coordinate index (c(1),c(2))'s partition
  function find_partition_bounds_coord(c) result(bounds)
    integer, dimension(2), intent(in)  :: c
    integer, dimension(4) :: bounds
    print*, 'So, the bounds of the partition block containing x are:'
    ! Lower bound direction 1
    bounds(1) = c(1)-merge(mod(c(1),q(1))-1,q(1)-1,mod(c(1),q(1)) .ne. 0)
    ! Upper bound direction 1
    bounds(2) = c(1)+merge(q(1)-mod(c(1),q(1)),0,mod(c(1),q(1)) .ne. 0)
    ! Lower bound direction 2
    bounds(3) = c(2)-merge(mod(c(2),q(2))-1,q(2)-1,mod(c(2),q(2)) .ne. 0)
    ! Upper bound direction 2
    bounds(4) = c(2)+merge(q(2)-mod(c(2),q(2)),0,mod(c(2),q(2)) .ne. 0)
    print*, 'LLB = (',bounds(1),', ',bounds(3),')'
    print*, 'ULB = (',bounds(1),', ',bounds(4),')'
    print*, 'LRB = (',bounds(2),', ',bounds(3),')'
    print*, 'URB = (',bounds(2),', ',bounds(4),')'
    print*, 'Per-rank block division is: ', bd(1), ' * ', bd(2)
  end function

  !!!!
  ! Given coordinate indices, return flat index
  function coords2flat(c) result(f)
    integer, dimension(2), intent(in)  :: c
    integer :: f
    f = (c(1)-1)*n(2) + c(2)
  end function

  ! Given flat index, return coordinate indices
  function flat2coords(f) result(c)
    integer, intent(in)   :: f
    integer, dimension(2) :: c
    c(1) = ceiling(real(f)/real(n(2)))
    c(2) = merge(n(2), mod(f,n(2)), mod(f,n(2)) .eq. 0)
  end function

end program get_ex

