!> \brief Support routines for dl_mg_defco module.
!!
!> These routines replace routines formerly provided by ONETEP and are mostly
!> based on the corresponding ONETEP routines, but make use of internal
!> DL_MG data and routines.
!>
!> Written by J. C. Womack, 2016
module dl_mg_defco_utils
  use dl_mg_params, only: wp
#ifdef MPI
  use dl_mg_mpi_header, only: MPI_ADDRESS_KIND
#endif
  implicit none

  private

  public :: dl_mg_defco_utils_assert
  public :: dl_mg_defco_utils_abort
  public :: dl_mg_defco_utils_alloc_check
  public :: dl_mg_defco_utils_dealloc_check
  public :: dl_mg_defco_utils_get_unit
  public :: dl_mg_defco_utils_integrate_product_on_grid

  ! TODO Remove need to additional MPI communicator and used DL_MG-specific
  !      communicator.
  ! MPI communicator to use for parallel operations in DL_MG defect correction
  ! wrapper
  !integer, public :: dl_mg_defco_comm
  ! TODO Is it necessary to have a upper bound for tags? This replicates
  !      functionality in ONETEP's comms module.
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND), public :: dl_mg_defco_tag_ub
#else
  ! Don't use MPI_ADDRESS_KIND if we do not have the mpi module
  integer, public :: dl_mg_defco_tag_ub
#endif

  contains

    !> Replaces ONETEP routine utils_assert with a routine that uses dl_mg
    !> error-handling routines. Written by JCW, June 2016.
    subroutine dl_mg_defco_utils_assert(condition,msg)
      use dl_mg_errors, only: handle_error, DL_MG_ERR_ASSERTION
      implicit none
      ! Arguments
      logical, intent(in) :: condition
      character(len=*), intent(in) :: msg

      ! If condition satisfied, return to calling unit:
      if(condition) return

      ! Else, abort with error:
      call handle_error(DL_MG_ERR_ASSERTION,msg=msg)
      ! NB. Do not pass ierror argument, so handle_error aborts regardless of
      !     value of abort_on_failure
    end subroutine dl_mg_defco_utils_assert

    !> Replaces ONETEP routine utils_abort with a routine that uses dl_mg
    !! error-handling routines. This routine will always raise an
    !! "unspecified" error, and the error should be described in the "msg"
    !! string. Written by JCW, June 2016.
    subroutine dl_mg_defco_utils_abort(msg)
      use dl_mg_errors, only: handle_error, DL_MG_ERR_UNSPECIFIED
      implicit none
      character(len=*), intent(in) :: msg
      ! Abort with error:
      call handle_error(DL_MG_ERR_UNSPECIFIED,msg=msg)
      ! NB. Do not pass ierror argument, so handle_error aborts regardless of
      !     value of abort_on_failure
    end subroutine dl_mg_defco_utils_abort

    !> Replaces ONETEP routine utils_alloc_check with a routine that uses dl_mg
    !> error-handling routines. Written by JCW, June 2016.
    subroutine dl_mg_defco_utils_alloc_check(status,routine,array)
      use dl_mg_errors, only: handle_error, DL_MG_ERR_ALLOCATION
      implicit none
      integer,intent(in) :: status !< status integer returned by allocate()
      character(len=*), intent(in) :: routine !< name of calling routine
      character(len=*), intent(in) :: array   !< name of array

      ! If ierr==0, return to calling unit:
      if (status==0) return
      ! Else, abort with error
      call handle_error(DL_MG_ERR_ALLOCATION,msg=routine//" "//array)
      ! NB. Do not pass ierror argument, so handle_error aborts regardless of
      !     value of abort_on_failure
    end subroutine dl_mg_defco_utils_alloc_check

    !> Replaces ONETEP routine utils_dealloc_check with a routine that uses dl_mg
    !> error-handling routines. Written by JCW, June 2016.
    subroutine dl_mg_defco_utils_dealloc_check(status,routine,array)
      use dl_mg_errors, only: handle_error, DL_MG_ERR_DEALLOCATION
      implicit none
      integer,intent(in) :: status !< status integer returned by deallocate()
      character(len=*), intent(in) :: routine !< name of calling routine
      character(len=*), intent(in) :: array   !< name of array

      ! If ierr==0, return to calling unit:
      if (status==0) return
      ! Else, abort with error
      call handle_error(DL_MG_ERR_DEALLOCATION,msg=routine//" "//array)
      ! NB. Do not pass ierror argument, so handle_error aborts regardless of
      !     value of abort_on_failure
    end subroutine dl_mg_defco_utils_dealloc_check

    !> Finds a free unit for I/O
    !! Based on utils_unit from ONETEP utils module, which was an amalgamation of
    !! restart_find_unit by Chris-Kriton Skylaris and esdf_unit by Chris Pickard.
    !! Modified by James Womack to make use of internal DL_MG objects.
    integer function dl_mg_defco_utils_get_unit()
      use dl_mg_errors, only: handle_error, DL_MG_ERR_IO

      ! Local variables
      integer :: trial_unit
      integer :: ierr
      logical :: ex
      logical :: op
      character(len=2) :: unit_str

      dl_mg_defco_utils_get_unit = -1
      do trial_unit=10,99
         inquire(unit=trial_unit,exist=ex,opened=op,iostat=ierr)
         if (ierr /= 0) then
            write(unit_str,'(i2)') trial_unit
            call handle_error(DL_MG_ERR_IO, &
                 msg='Error in dl_mg_defco_utils_get_unit: inquiring about unit '//&
                 trim(unit_str)//'failed.')
         end if
         if (ex .and. (.not. op)) then
            dl_mg_defco_utils_get_unit = trial_unit
            exit
         end if
      end do

      if (dl_mg_defco_utils_get_unit == -1) then
            call handle_error(DL_MG_ERR_IO, &
                 msg='Error in dl_mg_defco_utils_get_unit: no I/O units available')
      end if
    end function dl_mg_defco_utils_get_unit

    !> Integrate the functions represented on a 3D grid. Integration occurs over
    !! slices array slices (1:m1, 1:m2, 1:m3).
    !!
    !! If x2 is present, integration is over the product of x1 and x2.
    !! If x2 is not present, integration is over x1.
    !!
    !! Adapted from integrals_product_on_grid from ONETEP's integrals module,
    !! for inclusion in DL_MG. See original documentation in source.
    !! J. C. Womack, 2016
    real(kind=wp) function dl_mg_defco_utils_integrate_product_on_grid(weight, x1, x2,&
         m1, m2, m3, comm)
      !======================================================================!
      ! Returns the integral of a product of two functions                   !
      ! If x2 is omitted, integrates x1 only.                                !
      !----------------------------------------------------------------------!
      ! Arguments:                                                           !
      !   grid (input): The grid on which the integration takes place.       !
      !   x1, x2 (input): The quantities to integrate. x2 is optional.       !
      !   m1, m2, m3 (input): Integration region --- integration occurs over !
      !                 array slices 1:m1,1:m2,1:m3                          !
      ! Return value: Integral of x1(r)*x2(r) on the grid.                   !
      !----------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in June 2010.                              !
      !======================================================================!

      ! JCW: "use [...], only" syntax is problematic with mpi headers, so
      ! JCW: just import entire dl_mg_mpi_header module
      use dl_mg_mpi_header
      !use dl_mg_mpi_header, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_SUM, &
      !     MPI_DOUBLE_PRECISION, MPI_SUCCESS

      ! jd: Arguments
      integer, intent(in)                 :: m1, m2, m3
      real(kind=wp), intent(in)           :: x1(:, :, :)
      real(kind=wp), intent(in), optional :: x2(:, :, :)
      real(kind=wp), intent(in)     :: weight
      integer, intent(in), optional       :: comm !< MPI communicator for ALLREDUCE

      ! jd: Internal variables
      integer :: i1, i2, i3
      real(kind=wp) :: integral,aux
      integer :: ierr
      character(len=*), parameter :: myself = "integrate_product_on_grid"

      !------------------------------------------------------------------------

      ! Sanity checks
      call dl_mg_defco_utils_assert(m1>0 .and. m2>0 .and. m3>=0,&
           'Error in '//myself//': Bad integration bounds.')
      call dl_mg_defco_utils_assert(size(x1,1)>=m1.and.size(x1,2)>=m2.and.size(x1,3)>=m3,&
           'Error in '//myself//': Integration bounds greater than extents &
           &of x1 array')
      if (present(x2)) then
         call dl_mg_defco_utils_assert(size(x2,1)>=m1.and.size(x2,2)>=m2.and.size(x2,3)>=m3,&
              'Error in '//myself//': Integration bounds greater than extents &
              &of x2 array')
      endif

      integral = 0.0_wp
      ! jd: The following duplicates code in order not to check the
      !     if(present) condition repeatedly in a tight loop
      if(present(x2)) then
         do i3=1, m3
            do i2=1, m2
               do i1=1, m1
                  integral = integral + x1(i1,i2,i3) * x2(i1,i2,i3)
               end do
            end do
         end do
      else
         do i3=1, m3
            do i2=1, m2
               do i1=1, m1
                  integral = integral + x1(i1,i2,i3)
               end do
            end do
         end do
      end if

#ifdef MPI
      call dl_mg_defco_utils_assert(present(comm),"Error in "//myself//":&
           &No MPI communicator specified")
      aux = integral
      call MPI_ALLREDUCE(aux,integral,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
              comm,ierr)
      if (ierr /= MPI_SUCCESS) call dl_mg_defco_utils_abort("Error in "//myself//&
           "MPI_ALLREDUCE failed")
#endif
      ! jd: NB, the weight remains unchanged even if integration bounds overridden
      dl_mg_defco_utils_integrate_product_on_grid = integral * weight

    end function dl_mg_defco_utils_integrate_product_on_grid


end module dl_mg_defco_utils
