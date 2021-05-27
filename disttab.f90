program get_ex
  use mpi
  implicit none
  character (len=255) :: lookup_table_filename
  integer :: ierror, rank, nprocs
  integer, parameter, dimension(3) :: n = (/ 25, 10, 1 /)
  integer :: ncv_flat, n_flat
  integer, parameter, dimension(2) :: q = (/ 5, 2 /), p = (/ 3, 3 /)
  integer :: q_flat
  integer, dimension(3) :: bd
  double precision, dimension(:), allocatable :: partition
  integer :: x(2), y ! Coordinate index (x(1),x(2)), flat index y
  integer :: block_bounds(4)
  integer :: a, b, c, i, j, k, l, m, ro, co, bmi
  double precision :: value_fetched, z, d
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

  ncv_flat = product(n(1:ubound(n,dim=1)-1))
  n_flat = product(n)
  q_flat = product(q)

  ! Allocate lookup table
  allocate(lookup_table(((ncv_flat*rank) + 1):(ncv_flat*(rank+1)), n(3)))

  ! Open lookup table file
  !lookup_table_filename = 'tables/2d.dat'
  !open(unit = 46, file = lookup_table_filename, status='old', action='read')

  ! Compute size of lookup table in bytes for window creation
  lookup_table_size = n_flat*dbl_size

  ! Block division per-rank
  do i = 1, 3
    bd(i) = n(i)/nprocs
  enddo

  ! Populate lookup table with self-explanatory entries
  do j = 1,n(3)
    do k = lbound(lookup_table,dim=1), ubound(lookup_table,dim=1)
      lookup_table(k,j) = k + 0.d0
    enddo
  enddo

  ! Rearrange lookup table into block-major ordering
  call block_major_order(lookup_table)

  !print *, lookup_table

  ! Create memory window covering whole sub-table on each rank
  call mpi_win_create(lookup_table, lookup_table_size, dbl_size, mpi_info_null, mpi_comm_world, window, ierror)
  call mpi_win_fence(0, window, ierror)

  ! Pick a random number in the global lookup table range
  call random_number(z)
  y = ceiling(z*ncv_flat*nprocs)
  x = flat2coords(y)

  ! Print information about the point requested x
  print *, 'Flat entry y = ', y, ' has coordinates: ', x

  y = coords2flat(x)

  !print *,  'Verify: (', lookup_table(y,1), ')', &
  !  ' in partition block: ', find_local_partition(x), &
  !  'with local partition entry ', find_local_partition_entry(x)

  print *, 'Base of partition is ', bm_partition_base(find_local_partition(x))

  ! Print information about x's partition block
  block_bounds = find_partition_bounds_coord(x)

  !print *, "lookup table on rank ", rank, " has range ", find_rank_bounds(rank)

  ! Dynamically allocate one partition
  allocate(partition(1:q_flat))

  target_displacement = fm2bm(bm_partition_base(find_local_partition(x))) - ncv_flat*rank - 1

  print *, "RANK = ", rank, " Y = ", y, " DISPLACEMENT = ", target_displacement

  if (find_block_rank_flat(y) .ne. rank) then
    print *, 'Fetch block of ', q_flat, ' from rank ', find_block_rank_flat(y), ' to rank ', rank
    print *, 'The block containing ', y, ' begins with entry ', bm_partition_base(find_local_partition(x))
    print *, 'Which is located at index ', target_displacement
    call mpi_get(partition, q_flat, mpi_double, find_block_rank_flat(y), &
      target_displacement, q_flat, mpi_double, window, ierror)
    print *, partition
  endif


  ! Clean up, exit
  !close(46)
  deallocate(partition)
  call mpi_win_free(window, ierror)
  call mpi_finalize(ierror)

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! A bunch of helper functions containing ugly modular arithmetic. !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Reorder the chemistry table into block-major organization
  ! so that each block is contiguous in memory for MPI GET
  subroutine block_major_order(table)
    double precision, dimension(:,:) :: table
    m = ncv_flat*rank + 1
    do i = 1, n(1)/q(1)
      do j = 1, n(2)/q(2)
        x = flat2coords(bm_partition_base((/ i, j /)))
        do k = 0, q(1)-1
          do l = 0, q(2)-1
            bmi = coords2flat(x + (/ k,l /)) + ncv_flat*rank
            !print *, m, bmi
            if (m .gt. bmi) then
              z = lookup_table(m,1)
              lookup_table(m,1) = lookup_table(bmi,1)
              lookup_table(bmi,1) = z
            endif
            m = m + 1
          enddo
        enddo
      enddo
    enddo
  end subroutine

  ! Slowest way to convert flat end-major index to block-major index
  function fm2bm(fm) result(bm)
    integer, intent(in) :: fm
    integer :: bm
    m = 1
    do a = 0, nprocs-1
      do i = 1, n(1)/q(1)
        do j = 1, n(2)/q(2)
          x = flat2coords(bm_partition_base((/ i, j /)))
          do k = 0, q(1)-1
            do l = 0, q(2)-1
              bm = coords2flat(x + (/ k,l /)) + ncv_flat*rank
              if (m .eq. fm) return
              m = m + 1
            enddo
          enddo
        enddo
      enddo
    enddo
  end function

  ! Find the base of the partition in block-major organization
  function bm_partition_base(partition) result(base)
    integer, dimension(2), intent(in) :: partition
    integer :: base
    base = n(2)*q(1)*(partition(1)-1) + q(2)*(partition(2)-1) + 1
  end function

  ! Return upper and lower bounds of rank in global flat indices
  ! Assuming the tables on each rank are NOT indexed starting at 1
  function find_rank_bounds(bloc) result(bounds)
    integer, intent(in)  :: bloc
    integer, dimension(2) :: bounds
    bounds(1) = lbound(lookup_table,dim=1)!+ncv_flat*rank
    bounds(2) = ubound(lookup_table,dim=1)!+ncv_flat*rank
  end function

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
    ! Lower bound direction 1
    bounds(1) = c(1)-merge(mod(c(1),q(1))-1,q(1)-1,mod(c(1),q(1)) .ne. 0)
    ! Upper bound direction 1
    bounds(2) = c(1)+merge(q(1)-mod(c(1),q(1)),0,mod(c(1),q(1)) .ne. 0)
    ! Lower bound direction 2
    bounds(3) = c(2)-merge(mod(c(2),q(2))-1,q(2)-1,mod(c(2),q(2)) .ne. 0)
    ! Upper bound direction 2
    bounds(4) = c(2)+merge(q(2)-mod(c(2),q(2)),0,mod(c(2),q(2)) .ne. 0)
    !print*, 'LLB = (',bounds(1),', ',bounds(3),')'
    !print*, 'ULB = (',bounds(1),', ',bounds(4),')'
    !print*, 'LRB = (',bounds(2),', ',bounds(3),')'
    !print*, 'URB = (',bounds(2),', ',bounds(4),')'
  end function

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

