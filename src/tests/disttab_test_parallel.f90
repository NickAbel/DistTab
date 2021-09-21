!> This module contains functions pertaining to testing the parallel
!! capabilities of the DistTab library.
module disttab_test_parallel
  use :: communicator
  use :: disttab_table
  use :: kind_params

  implicit none
  private
  public :: parallel_test

  type :: parallel_test
    private
    type(table) :: lookup
    type(comm) :: disttab_comm

  contains

    procedure, pass(this) :: run_comm_test

    procedure, pass(this), private :: comm_test
    procedure, pass(this), private :: fill_table_ascending_integers
    procedure, pass(this), private :: fill_cvars_linspace

    final :: parallel_test_destructor

  end type parallel_test

  interface parallel_test
    module procedure :: parallel_test_constructor
  end interface parallel_test

contains

!> The constructor for the parallel_test type.
!!
!! @param table_dimensions specifies the size of the table
!! @return this the parallel_test object created
  type(parallel_test) function parallel_test_constructor(table_dimensions) result(this)
    integer(i4), dimension(:), intent(in) :: table_dimensions

    this % lookup = table(table_dimensions)
    call this % fill_table_ascending_integers()

  end function parallel_test_constructor

!> destructor for the parallel test type
!!
!! @param this the parallel test object to destruct
  subroutine parallel_test_destructor(this)
    type(parallel_test) :: this

    call this % lookup % deallocate_table()

  end subroutine parallel_test_destructor

!> Runs the get value test.
!! @param this parallel_test object
  subroutine comm_test(this)
    class(parallel_test), intent(inout) :: this

    write (*, *) "comm_test"

  end subroutine comm_test

!> Fill lookup's table with ascending integers.
!! @param this parallel_test object
  subroutine fill_table_ascending_integers(this)
    class(parallel_test), intent(inout) :: this
    integer(i4) :: i

    do i = 1, this % lookup % table_dims_flat
      this % lookup % elems(:, i) = i
    end do

  end subroutine fill_table_ascending_integers

!> Fill lookup's control variables with linspace values in [0, 1].
!! @param this parallel_test object
  subroutine fill_cvars_linspace(this)
    class(parallel_test), intent(inout) :: this
    integer(i4) :: i, N, j

    N = size(this % lookup % part_dims)

    do i = 1, N
      do j = 1, this % lookup % table_dims(i)
        this % lookup % ctrl_vars(j + sum(this % lookup % table_dims(:i - 1))) = &
                        & 1.0 * (j - 1) / (this % lookup % table_dims(i) - 1)
      end do
    end do

  end subroutine fill_cvars_linspace

!> Runs the get value test.
!! @param this parallel_test object
  subroutine run_comm_test(this)
    class(parallel_test), intent(inout) :: this

    call this % comm_test()

  end subroutine run_comm_test

end module disttab_test_parallel
