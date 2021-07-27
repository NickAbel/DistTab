module disttab_test_access
   use :: disttab_table

   implicit none
   private
   public :: access_test

   type :: access_test
      private
      type(table) :: lookup
   contains
      procedure, pass(this) :: run_value_test
      procedure, pass(this) :: run_value_cloud_test

      procedure, pass(this), private :: get_value_test
      procedure, pass(this), private :: get_value_cloud_test

   end type access_test

   interface access_test
      module procedure :: access_test_constructor
   end interface access_test

contains

  !> The constructor for the access_test type.
  !! Initializes the table given table dimensions.
  !!
  !! @param table_dimensions specifies the size of the table
  !! @return this the access_test object created
   type(access_test) function access_test_constructor(table_dimensions) result(this)
      integer, dimension(:), intent(in) :: table_dimensions

      this%lookup = table(table_dimensions)

   end function access_test_constructor

  !> Test the ability to get correct values from the table from index,
  !! local coordinate, and global coordinates.
   subroutine get_value_test(this)
      class(access_test), intent(inout) :: this

      print *, "get_value_test begin"

   end subroutine get_value_test 

  !> Test the ability to get correct value clouds from the table from index,
  !! local coordinate, and global coordinates.
  !! Example in 3D case:
  !! In general, when the value at coordinate (i,j,k) is desired, return
  !! {i,i+1}x{j,j+1}x{k,k+1} (returning err if i+1,j+1,k+1, etc > Z).
   subroutine get_value_cloud_test(this)
      class(access_test), intent(inout) :: this

      print *, "get_value_cloud_test begin"

   end subroutine get_value_cloud_test 


  !> Runs the get value test.
  !! @param this access_test object
   subroutine run_value_test(this)
      class(access_test), intent(inout) :: this

      call this%get_value_test()

   end subroutine run_value_test

  !> Runs the get value cloud test.
  !! @param this access_test object
   subroutine run_value_cloud_test(this)
      class(access_test), intent(inout) :: this

      call this%get_value_cloud_test()

   end subroutine run_value_cloud_test

end module disttab_test_access
