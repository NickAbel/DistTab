program disttab
   use :: mpi
   use :: disttab_table

   implicit none

   integer :: ierror
   integer, dimension(:), allocatable :: n, q

   type(table) :: lookup

   call read_command_args()

   call mpi_init(ierror)

   lookup = table(n, q)

   !call lookup % fill_example()
   !call lookup % read_in(file_id)
   call lookup%partition_mapping_test()
   !call lookup % partition_mapping()
   call lookup%table_deallocate()

   call mpi_finalize(ierror)

contains

   ! Read in the dimensions for table and block sizes, allocate relevant arrays
   ! according to the input sizes
   subroutine read_command_args()
      character(len=32) :: arg
      !character(len=120) :: file_id
      character(len=32), dimension(:), allocatable :: nstr, qstr, npad, coords, x
      integer :: n_len, q_len, i

      call getarg(1, arg)
      read (arg, *) n_len

      allocate (nstr(n_len))
      allocate (npad(n_len))
      allocate (n(n_len))
      allocate (coords(n_len - 1))
      allocate (x(n_len - 1))

      do i = 2, n_len + 1
         call getarg(i, arg)
         nstr(i - 1) = arg
      end do

      call getarg(n_len + 2, arg)
      read (arg, *) q_len

      allocate (qstr(q_len))
      allocate (q(q_len))

      do i = n_len + 3, n_len + q_len + 2
         call getarg(i, arg)
         qstr(i - (n_len + 2)) = arg
      end do

      !  call getarg(i, file_id)
      !
      !  print *, i, file_id
      read (nstr, *) n
      read (qstr, *) q

   end subroutine read_command_args

end program disttab
