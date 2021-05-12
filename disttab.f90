program get_ex
  use mpi
  implicit none
  character (len=255) :: lookup_table_filename
  integer :: ierror, rank, nprocs
  integer, parameter :: n1 = 25, n2 = 10, n3 = 1, p1 = 3, p2 = 3
  integer, parameter :: q1 = 5, q2 = 2
  integer :: x1, x2, x, y ! Coordinate index (x1,x2), flat index x
  integer :: block_bounds(4)
  integer :: l, m, n, bd1, bd2, flat_index, coord_index(2)
  double precision :: lookup_table(n1*n2,n3), values_fetched(q1), k
  integer :: window
  integer :: integer_size, dbl_size
  integer(kind=mpi_address_kind) :: win_size, lookup_table_size
  integer(kind=mpi_address_kind) :: target_displacement
  interface find_block_rank

    ! Return rank where a flat index x is located
    subroutine find_block_rank_flat(x)
      integer, intent(in) :: x
    end subroutine

    ! Return rank where a coordinate index (x1,x2) is located
    subroutine find_block_rank_coord(x1,x2)
      integer, intent(in) :: x1,x2
    end subroutine

  end interface find_block_rank

  interface find_partition_flat

    ! Return all flat indices of flat index x's partition
    subroutine find_partition_flat2flat(x)
      integer, intent(in) :: x
    end subroutine

    ! Return all flat indices of coordinate index (x1,x2)'s partition
    subroutine find_partition_coord2flat(x1,x2)
      integer, intent(in) :: x1,x2
    end subroutine

  end interface find_partition_flat

  interface find_partition_coord

    ! Return all coordinate indices of flat index x's partition
    subroutine find_partition_flat2coord(x)
      integer, intent(in) :: x
    end subroutine

    ! Return all coordinate indices of coordinate index (x1,x2)'s partition
    subroutine find_partition_coord2coord(x1,x2)
      integer, intent(in) :: x1,x2
    end subroutine

  end interface find_partition_coord

  ! Setup MPI communicator, retrieve MPI type sizes
  call mpi_init(ierror)
  call mpi_type_size(mpi_integer, integer_size, ierror)
  call mpi_type_size(mpi_double, dbl_size, ierror)
  call mpi_comm_rank(mpi_comm_world, rank, ierror)
  call mpi_comm_size(mpi_comm_world, nprocs, ierror)

  ! Open lookup table file
  lookup_table_filename = 'tables/2d.dat'
  open(unit = 46, file = lookup_table_filename, status='old', action='read')

  ! Compute size of lookup table and window block size
  lookup_table_size = n1*n2*n3*dbl_size
  win_size = q1*dbl_size

  ! Insist on some comm size
  if (nprocs .ne. 2) then
    print *, "use 2 mpi processes for now, we detect ", nprocs
    call mpi_abort(mpi_comm_world, -1, ierror)
  endif

  ! Block division per-rank
  bd1 = n1/nprocs
  bd2 = n2/nprocs

  ! Populate lookup table with self-explanatory entries
  do l = 1,n3
    do m = 1,n2
      do n = 1,n1
        call coords2flat(n,m,n2,flat_index)
        coord_index(1) = ceiling(real(flat_index)/real(n2))
        coord_index(2) = merge(m, mod(flat_index,n2), m .eq. n2)
        !!! Debug output for indexing conversion
        !if (rank .eq. 0) then
        !  print *, n-coord_index(1), m-coord_index(2), flat_index,n2
        !endif 
        lookup_table(flat_index,l) = n1*n2*rank+flat_index !1000*coord_index(1) + coord_index(2)
      enddo
    enddo
  enddo

  ! Create memory window covering whole sub-table on each rank
  call mpi_win_create(lookup_table, q1*win_size, dbl_size, mpi_info_null, mpi_comm_world, window, ierror)
  call mpi_win_fence(0, window, ierror)

  ! If block is available on own rank... print it, for now
  call random_number(k)
  x = ceiling(k*500)
  call flat2coords(x,n2,x1,x2)
  call coords2flat(x1,x2,n2,y)
  ! Print information about the point requested x
  print*, 'Flat entry x = ', x, ' has coordinates: ', x1, x2, 'Verify: (', lookup_table(y,1), ')', &
    ' in partition block: ', ceiling(real(x1)/real(q1)), ceiling(real(x2)/real(q2)), &
    'with local partition entry ', merge(mod(x1,q1),q1,mod(x1,q1) .ne. 0), merge(mod(x2,q2),q2,mod(x2,q2) .ne. 0)
  ! Print information about x's partition block
  print*, 'So, the bounds of the partition block containing x are:'
  ! Lower bound direction 1
  block_bounds(1) = x1-merge(mod(x1,q1)-1,q1-1,mod(x1,q1) .ne. 0)
  ! Upper bound direction 1
  block_bounds(2) = x1+merge(q1-mod(x1,q1),0,mod(x1,q1) .ne. 0)
  ! Lower bound direction 2
  block_bounds(3) = x2-merge(mod(x2,q2)-1,q2-1,mod(x2,q2) .ne. 0)
  ! Upper bound direction 2
  block_bounds(4) = x2+merge(q2-mod(x2,q2),0,mod(x2,q2) .ne. 0)
  print*, 'LLB = (',block_bounds(1),', ',block_bounds(3),')'
  print*, 'ULB = (',block_bounds(1),', ',block_bounds(4),')'
  print*, 'LRB = (',block_bounds(2),', ',block_bounds(3),')'
  print*, 'URB = (',block_bounds(2),', ',block_bounds(4),')'
  print*, 'Per-rank block division is: ', bd1, ' x ', bd2

  print *, "lookup table on rank ", rank, " starts at ", lbound(lookup_table,dim=1)+n1*n2*rank, &
    ", ends at ", ubound(lookup_table,dim=1)+n1*n2*rank

  do m = block_bounds(1), block_bounds(2)
    do n = block_bounds(3), block_bounds(4)
      call coords2flat(m,n,n2,y)
      print *, m,n, " corresponds to x = ", y, merge("     (On table at rank ", " (OUT OF TABLE ON RANK ", y >= &
        lbound(lookup_table,dim=1)+n1*n2*rank  .and. y <= ubound(lookup_table,dim=1)+n1*n2*rank), rank, ") "
    enddo
  enddo

  ! If block is available on other rank, get it and... print it, for now
  if (rank .eq. 0) then
    target_displacement = 0
    call mpi_get(values_fetched, q1, mpi_double, 1, target_displacement, q1, mpi_double, window, ierror)
  else if (rank .eq. 1) then
    target_displacement = 0
    call mpi_get(values_fetched, q1, mpi_double, 0, target_displacement, q1, mpi_double, window, ierror)
  endif
  call mpi_win_fence(0, window, ierror)

  ! Print, clean up, exit
  !print *, "Fetched value = ", values_fetched, rank
  call mpi_win_free(window, ierror)
  close(46)
  call mpi_finalize(ierror)
end program get_ex

subroutine coords2flat(y1,y2,n2,y)
  integer, intent(in) :: y1,y2,n2
  integer, intent(out) :: y
  y = (y1-1)*n2 + y2
end subroutine

subroutine flat2coords(y,n2,y1,y2)
  integer, intent(in) :: y,n2
  integer, intent(out) :: y1,y2
  y1 = ceiling(real(y)/real(n2))
  y2 = merge(n2, mod(y,n2), mod(y,n2) .eq. 0)
end subroutine

