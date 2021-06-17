program disttab
  use :: mpi
  use :: disttab_table

  implicit none

  integer :: ierror, rank, nprocs
  integer :: ncv_flat, n_flat, ncv_padded, n_padded
  integer :: q_flat
  double precision, dimension(:), allocatable :: partition
  integer, dimension(:), allocatable :: n, q, nloop_counters, nloop_uppers
  integer, dimension(:), allocatable :: npad, coords, x
  character(len=32), dimension(:), allocatable :: nstr, qstr
  integer :: y ! Coordinate index (x(1),x(2)), flat index y
  integer :: a, b, c, d, i, j, k, l, m, bmi, n_len, q_len
  ! Note that option 2 is used for the reordering unit test
  double precision :: z
  double precision, dimension(:,:), allocatable :: lookup_table, lookup_table_flat
  integer :: window
  integer :: integer_size, dbl_size
  integer(kind=mpi_address_kind) :: lookup_table_size
  integer(kind=mpi_address_kind) :: target_displacement
  character(len=32) :: arg
  character(len=120) :: file_id

  type(table) :: lookup

  ! Read in the dimensions for table and block sizes, allocate relevant arrays
  ! according to the input sizes
  call getarg(1, arg)
  read (arg, *) n_len

  allocate(nstr(n_len))
  allocate(npad(n_len))
  allocate(n(n_len))
  allocate(coords(n_len-1))
  allocate(x(n_len-1))

  do i = 2, n_len + 1
    call getarg(i, arg)
    nstr(i-1) = arg
  enddo

  call getarg(n_len + 2, arg)
  read (arg, *) q_len

  allocate(qstr(q_len))
  allocate(q(q_len))

  do i = n_len + 3, n_len + q_len + 2
    call getarg(i, arg)
    qstr(i - (n_len + 2)) = arg
  enddo

  call getarg(i, file_id)

  print *, i, file_id
  read (nstr, *) n
  read (qstr, *) q

  call mpi_init(ierror)

  lookup = table(n, q)

  !call lookup % fill_example()
  call lookup % read_in(file_id)
  call lookup % partition_mapping_test()
  write (*,*) " "
  call lookup % partition_mapping()
  call lookup % table_deallocate()

  call mpi_finalize(ierror)

end program disttab 
