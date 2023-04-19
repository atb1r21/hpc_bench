!==============================================================================!
!                                                                              !
!                                 <Module name>                                !
!                                                                              !
!                       <Author name(s)>, <Date of creation>                   !
!                                                                              !
!------------------------------------------------------------------------------!
! # MODULE DESCRIPTION                                                         !
! Write a general description of the functionality offered by the module here. !
!------------------------------------------------------------------------------!
! # BACKGROUND                                                                 !
! Optionally provide a more detailed description of any theoretical or         !
! implementation details which may help developers working with the module.    !
! This is useful, for example, if the module implements complicated or subtle  !
! theory, the implementation is complex or difficult to follow, or there are   !
! important data structures which require further explanation.                 !
!------------------------------------------------------------------------------!
! # NOTES                                                                      !
! Write notes about the module which may be relevant to other developers, e.g. !
! known issues or constraints.                                                 !
! ## TODO LIST                                                                 !
! * List of known tasks that remain to be completed                            !
!==============================================================================!

!------------------ DELETE THIS WHEN CREATING A NEW MODULE --------------------!
! # USING THIS TEMPLATE                                                        !
! ## General guidelines                                                        !
! (These are adapted from `CONTRIBUTING.markdown`. See this document for more  !
! guidance on ONETEP coding conventions and guidelines.)                       !
!                                                                              !
! Before creating a new module, you should consider carefully whether your new !
! functionality fits within the framework of an existing module, or is generic !
! enough to be part of a multi-purpose module, such as `utils` or `services`.  !
! If a new module is needed to encapsulate some new functionality, then you    !
! should follow the following guidelines:                                      !
!                                                                              !
! * Give your module a name which indicates the functionality it contains.     !
! * The filename for the module should have the form `<module_name>_mod.F90`   !
!   where `<module_name>` is the name you have given the module.               !
! * By default, variables and procedures in your module should be private      !
!   (i.e. they should have the `private` attribute).                           !
! * Global variables (private or public variables declared at the level of     !
!   the module itself rather than within its routines) constitute "hidden      !
!   state". Think carefully before declaring any module-level global variables !
!   and avoid doing so where it is reasonably possible.                        !
! * Variables and procedures which do have to be public (accessible outside    !
!   the module) should be explicitly specified (i.e. they should have the      !
!   `public` attribute).                                                       !
! * In general, public variable and procedure names should be prepended by a   !
!   standard prefix (typically the module name, or a shortened version of the  !
!   name).                                                                     !
!                                                                              !
! ## Optional elements                                                         !
! This template contains some elements which may not be relevant to all        !
! modules. For example, the functionality in the module may not lend itself to !
! unit testing or the use of initialization or de-initialization routines.     !
! Please remove elements that are not suited for your module.                  !
!                                                                              !
! ## Markdown                                                                  !
! It is recommended that long-form comment text (e.g. in the module or         !
! procedure headers is written with Markdown formatting. This enables text to  !
! be clearly formatted into sections and bullet-point lists, but should also   !
! facilitate automated extraction and formatting of documentation in the       !
! future.                                                                      !
!                                                                              !
! Markdown is simple to learn and easy to apply, see e.g.                      !
! [CommonMark help page](https://commonmark.org/help/)                         !
!                                                                              !
! Note that, for clarity, Markdown conventions are not fully followed in the   !
! template. For example in line with usual man-page syntax, we use angle       !
! brackets to denote text which should be replaced in the actual module (e.g.  !
! <Module name>). This conflicts with the interpretation of text between angle !
! brackets as HTML tags in Markdown. Since these are to be replaced when the   !
! module is written, this should not be an issue.                              !
!                                                                              !
! In general, it is recommended to use Markdown syntax on a "best effort"      !
! basis, to provide documentation with a readable and clear structure. Where   !
! using Markdown conflicts with clarity (as with the case of text within angle !
! brackets), priority should be placed on clarity.                             !
!                                                                              !
! ## AUTHOR(S) & CHANGELOG                                                     !
! Author(s):        James C. Womack                                            !
! Date of creation: 09/2018                                                    !
! List of major changes:                                                       !
! 11/09/18, Initial creation of template module, J. C. Womack                  !
!------------------------------------------------------------------------------!

module template

  ! # MODULE-WIDE INCLUDES
  use comms, only : pub_on_root
  use constants, only: DP, stdout, stderr

  implicit none
  private

  ! # PUBLIC PROCEDURES
  ! Initialize module
  public :: template_initialise
  ! De-initialize module (reset module state to default)
  public :: template_exit
  ! Example public procedure
  public :: template_pub_proc
  ! Unit tests for module procedures
  public :: template_run_unit_tests

  ! # PUBLIC MODULE VARIABLES AND TYPES
  ! Example public variable
  logical, public :: template_example_pub_var
  !   ^-- Module global variables should be avoided where possible to avoid
  !       introducing 'hidden state'.

  ! Example public type
  public :: EXAMPLE_TYPE
  type EXAMPLE_TYPE
    integer :: a
    integer :: b
    integer :: c
  end type EXAMPLE_TYPE

  ! # PRIVATE MODULE VARIABLES AND TYPES
  ! Example private variable
  logical, save :: example_priv_var
  !   ^-- Module global variables should be avoided where possible to avoid
  !       introducing 'hidden state'.

  ! Short code for reporting (useful for grepping output)
  ! This is also the block name for devel codes
  character(len=*), parameter :: module_short_code = "TMPLT"

  ! # INTERFACE BLOCKS
  ! Interface blocks may be included here, if needed (e.g. for overloading
  ! procedure interfaces)

contains
  ! ## PUBLIC PROCEDURES

  subroutine template_initialise()
    !==========================================================================!
    ! # DESCRIPTION                                                            !
    ! Initialize the <template> module:                                        !
    ! * List of steps involved in module initialization                        !
    ! * ...                                                                    !
    !--------------------------------------------------------------------------!
    ! # ARGUMENTS                                                              !
    ! <name> <in/out/inout> <arg descrption>                                   !
    ! None                                                                     !
    !--------------------------------------------------------------------------!
    ! # AUTHORS & CHANGELOG                                                    !
    ! Author(s):        <Author name(s)>                                       !
    ! Date of creation: <Date of creation>                                     !
    ! List of major changes:                                                   !
    ! <date> <change description> <author>                                     !
    !==========================================================================!
    use constants, only: NORMAL
    use rundat, only: pub_devel_code, pub_output_detail
    use utils, only: utils_trace_in, utils_trace_out, &
         utils_abort, utils_devel_code

    implicit none

    ! Parameters
    character(len=*), parameter :: myself = "template_initialise"

    ! Local variables
    logical :: run_unit_tests
    integer :: ierror

    call utils_trace_in(myself)

    ! Insert initialization code here. If initialization is not required, this
    ! routine can be deleted.

    !------- Example code for calling unit tests (delete if not needed) -------!
    ! This is a good place to run unit tests, if this module has any.
    ! It is recommended to activate these with a devel code value. In this
    ! case, the devel_code `TMPLT:UNIT_TEST=[T|F]:TMPLT` is checked, where
    ! TMPLT is the module parameter `module_short_code`.
    run_unit_tests = &
         utils_devel_code(.false.,module_short_code,"UNIT_TEST",pub_devel_code)
    if (run_unit_tests) then
       if (pub_on_root) then
          write(stdout,"(a)") module_short_code//": Running unit tests."
       end if
       call template_run_unit_tests(ierror)
       if (ierror /= 0) then
          call utils_abort("Error in "//myself//": unit tests failed.")
       else
          write(stdout,"(a)") module_short_code//": Unit tests passed."

       end if
    end if
    !--------------------------------------------------------------------------!

    call utils_trace_out(myself)

  end subroutine template_initialise

  subroutine template_exit()
    !==========================================================================!
    ! # DESCRIPTION                                                            !
    ! De-initialize / clean up the <template> module:                          !
    ! * List of steps involved in module initialization                        !
    ! * ...                                                                    !
    !--------------------------------------------------------------------------!
    ! # ARGUMENTS                                                              !
    ! <name> <in/out/inout> <arg descrption>                                   !
    ! None                                                                     !
    !--------------------------------------------------------------------------!
    ! # AUTHORS & CHANGELOG                                                    !
    ! Author(s):        <Author name(s)>                                       !
    ! Date of creation: <Date of creation>                                     !
    ! List of major changes:                                                   !
    ! <date> <change description> <author>                                     !
    !==========================================================================!
    use utils, only: utils_trace_in, utils_trace_out

    implicit none

    ! Parameters
    character(len=*), parameter :: myself = "template_exit"

    ! Local variables

    call utils_trace_in(myself)

    ! Insert de-initialization code here. If initialization is not required,
    ! this routine can be deleted.

    call utils_trace_out(myself)

  end subroutine template_exit

  subroutine template_pub_proc(example_arg1,example_arg2,example_arg3)
    !==========================================================================!
    ! # DESCRIPTION                                                            !
    ! Write a description of the purpose and functionality of this procedure   !
    ! here.                                                                    !
    !--------------------------------------------------------------------------!
    ! # ARGUMENTS                                                              !
    ! <name> <in/out/inout> <arg descrption>                                   !
    ! example_arg1   out    An example argument with intent out.               !
    ! example_arg1   inout  An example argument with intent in.                !
    ! example_arg2   in     An example argument with intent in (optional)      !
    !--------------------------------------------------------------------------!
    ! # AUTHORS & CHANGELOG                                                    !
    ! Author(s):        <Author name(s)>                                       !
    ! Date of creation: <Date of creation>                                     !
    ! List of major changes:                                                   !
    ! <date> <change description> <author>                                     !
    !==========================================================================!
    use utils, only: utils_trace_in, utils_trace_out

    implicit none

    ! Parameters
    character(len=*), parameter :: myself = "template_pub_proc"

    ! Arguments
    real(kind=DP), intent(out)    :: example_arg1(:)
    integer, intent(inout)        :: example_arg2
    logical, intent(in), optional :: example_arg3

    ! Local variables
    logical :: example_arg3_loc

    call utils_trace_in(myself)

    ! Insert procedure code here

    !------ Example code for copying optional arg (delete if not needed) ------!
    if (present(example_arg3)) then
       example_arg3_loc = example_arg3
    else
       example_arg3_loc = .false.
    end if
    !--------------------------------------------------------------------------!

    call utils_trace_out(myself)

  end subroutine template_pub_proc

  subroutine template_run_unit_tests(ierror)
    !==========================================================================!
    ! # DESCRIPTION                                                            !
    ! Run unit tests for <template> module and output results to stdout.       !
    !                                                                          !
    ! An error code is returned. If this is zero, all tests passed and no      !
    ! discrepancies in outcomes were detected between MPI ranks. A non-zero    !
    ! error code indicates test failure and/or discrepancies between ranks.    !
    !--------------------------------------------------------------------------!
    ! # ARGUMENTS                                                              !
    ! <name>    <in/out/inout> <arg descrption>                                !
    ! ierror    out            error code (== 0 success, != 0 failure)         !
    !--------------------------------------------------------------------------!
    ! # AUTHORS & CHANGELOG                                                    !
    ! Author(s):        <Author name(s)>                                       !
    ! Date of creation: <Date of creation>                                     !
    ! List of major changes:                                                   !
    ! <date> <change description> <author>                                     !
    !==========================================================================!
    use comms, only: pub_total_num_nodes, comms_allgather
    use utils, only: utils_trace_in, utils_trace_out, &
         utils_alloc_check, utils_dealloc_check

    implicit none

    ! Parameters
    character(len=*), parameter :: myself = "template_run_unit_tests"
    integer, parameter :: IFAIL = 1
    integer, parameter :: IPASS = 2
    integer, parameter :: ISKIP = 3
    integer, parameter :: TEST_UNDEFINED = 1
    integer, parameter :: TEST_PASSED = 0
    integer, parameter :: TEST_FAILED = -1
    integer, parameter :: TEST_SKIPPED = -2
    integer, parameter :: ERROR_FAILED_TEST  = 1000

    ! Arguments
    integer, intent(out) :: ierror ! 0 for success

    ! Local variables
    ! - General
    integer :: itest   ! count total number of tests
    integer :: outcome ! integer representation of test outcome

    integer :: my_results(3) ! test results
      ! results(IFAIL) counts fails
      ! results(IPASS) counts passes
      ! results(ISKIP) counts skips

    call utils_trace_in(myself)

    !-- Example implementation of unit test framework (delete if not needed) --!
    ! This is a template which provides a suggested structure for implementing
    ! unit tests for a module.
    !
    ! In this example, the individual unit tests are implemented as internal
    ! procedures within the `template_run_unit_tests` subroutine. Each test
    ! returns the integer `outcome` (which has the value of one of the
    ! parameters TEST_UNDEFINED, TEST_PASSED, TEST_FAILED, TEST_SKIPPED), and
    ! overall results (fails, passes, skips) are counted in the array
    ! `my_results`.
    !
    ! Once all the tests are completed, the results are reported to `stdout`
    ! and the `template_run_unit_tests` routine returns a single integer error
    ! code indicating which indicates the outcome of the overall tests.
    !
    ! This is a serial unit testing framework. If you need tests to make use
    ! of MPI or OpenMP parallelism, you will need to take care to ensure the
    ! results are correctly communicated between multiple ranks.

    ! Assume the best, unless otherwise informed
    ierror = 0

    ! Initialize counters
    itest = 0
    my_results(1:3) = 0

    ! START OF UNIT TESTS
    call internal_test_example_priv_proc(outcome)
    call internal_update_results(itest,my_results,outcome)
    ! END OF UNIT TESTS

    ! REPORT RESULTS
    if (pub_on_root) then
       write(stdout,"(a)")      module_short_code//": [Unit test report]"
       write(stdout,"(a,i0,a)") module_short_code//": ", itest, " tests run"
       write(stdout,"(a,i0,a)") module_short_code//": * ", &
            my_results(IPASS), " passed"
       write(stdout,"(a,i0,a)") module_short_code//": * ", &
            my_results(IFAIL), " failed"
       write(stdout,"(a,i0,a)") module_short_code//": * ", &
            my_results(ISKIP), " skipped"
    end if

    ! Non-zero error code in event of failed tests
    if (my_results(IFAIL) > 0) ierror = ierror + ERROR_FAILED_TEST
    !--------------------------------------------------------------------------!

    call utils_trace_out(myself)

    contains

      subroutine internal_update_results(itest,my_results,outcome)
        !======================================================================!
        ! # INTERNAL PROCEDURE DESCRIPTION                                     !
        ! Host:      template_run_unit_tests                                  !
        ! Function:  Update test and results counters based on outcome of      !
        !            a unit test                                               !
        !----------------------------------------------------------------------!
        ! # ARGUMENTS                                                          !
        ! <name>      <in/out/inout> <arg descrption>                          !
        ! itest       inout          test counter (increment by 1)             !
        ! my_results  inout          array counting of passes, fails and skips !
        ! outcome     inout          integer outcome code of preceding test    !
        !                            set to TEST_UNDEFINED on exit             !
        !----------------------------------------------------------------------!
        ! # AUTHORS & CHANGELOG                                                !
        ! Author(s):        <Author name(s)>                                   !
        ! Date of creation: <Date of creation>                                 !
        ! List of major changes:                                               !
        ! <date> <change description> <author>                                 !
        !======================================================================!
        use utils, only: utils_abort

        implicit none

        ! Parameters
        character(len=*), parameter :: myself = &
             "template_run_unit_tests::internal_update_results"

        ! Arguments
        integer,intent(inout) :: itest
        integer,intent(inout) :: my_results(3)
        integer,intent(inout) :: outcome

        ! Increment test counter
        itest = itest + 1

        select case (outcome)

        case (TEST_PASSED)
           my_results(IPASS) = my_results(IPASS) + 1

        case (TEST_FAILED)
           my_results(IFAIL) = my_results(IFAIL) + 1

        case (TEST_SKIPPED)
           my_results(ISKIP) = my_results(ISKIP) + 1

        case default
           call utils_abort("Error in "//myself//": &
                &outcome code not recognized.")
        end select

        ! Reset outcome to undefined state
        outcome = TEST_UNDEFINED

      end subroutine internal_update_results

      subroutine internal_test_example_priv_proc(outcome)
        !======================================================================!
        ! # INTERNAL PROCEDURE DESCRIPTION                                     !
        ! Host:      template_run_unit_tests                                   !
        ! Function:  Unit test for <example_priv_proc> routine                 !
        !----------------------------------------------------------------------!
        ! ## TESTS                                                             !
        ! * List of tests being performed in this procedure                    !
        ! * ...                                                                !
        !                                                                      !
        !----------------------------------------------------------------------!
        ! # ARGUMENTS                                                          !
        ! <name>      <in/out/inout> <arg descrption>                          !
        ! outcome     out            integer outcome code for this unit test   !
        !----------------------------------------------------------------------!
        ! # AUTHORS & CHANGELOG                                                !
        ! Author(s):        <Author name(s)>                                   !
        ! Date of creation: <Date of creation>                                 !
        ! List of major changes:                                               !
        ! <date> <change description> <author>                                 !
        !======================================================================!
        implicit none

        ! Parameters
        character(len=*), parameter :: unit_test_routine = &
             "template_run_unit_tests::internal_test_example_priv_proc"

        ! Arguments
        integer,intent(out) :: outcome

        ! Local variables
        ! ...

        ! Assume everything is fine, unless it's not
        outcome = TEST_PASSED

        ! Insert unit test code here for the routine being tested. The value of
        ! `outcome` should be set to TEST_PASSED, TEST_FAILED, or TEST_SKIPPED
        ! based on the outcome of your tests

      end subroutine internal_test_example_priv_proc

      subroutine internal_test_template_pub_proc(outcome)
        !======================================================================!
        ! # INTERNAL PROCEDURE DESCRIPTION                                     !
        ! Host:      template_run_unit_tests                                   !
        ! Function:  Unit test for <template_pub_proc> routine                 !
        !----------------------------------------------------------------------!
        ! ## TESTS                                                             !
        ! * List of tests being performed in this procedure                    !
        ! * ...                                                                !
        !                                                                      !
        !----------------------------------------------------------------------!
        ! # ARGUMENTS                                                          !
        ! <name>      <in/out/inout> <arg descrption>                          !
        ! outcome     out            integer outcome code for this unit test   !
        !----------------------------------------------------------------------!
        ! # AUTHORS & CHANGELOG                                                !
        ! Author(s):        <Author name(s)>                                   !
        ! Date of creation: <Date of creation>                                 !
        ! List of major changes:                                               !
        ! <date> <change description> <author>                                 !
        !======================================================================!
        implicit none

        ! Parameters
        character(len=*), parameter :: unit_test_routine = &
             "template_run_unit_tests::internal_test_template_pub_proc"

        ! Arguments
        integer,intent(out) :: outcome

        ! Local variables
        ! ...

        ! Assume everything is fine, unless it's not
        outcome = TEST_PASSED

        ! Insert unit test code here for the routine being tested. The value of
        ! `outcome` should be set to TEST_PASSED, TEST_FAILED, or TEST_SKIPPED
        ! based on the outcome of your tests

      end subroutine internal_test_template_pub_proc

  end subroutine template_run_unit_tests

  ! ## PRIVATE PROCEDURES

  subroutine example_priv_proc(example_arg1,example_arg2,example_arg3)
    !==========================================================================!
    ! # DESCRIPTION                                                            !
    ! Write a description of the purpose and functionality of this procedure   !
    ! here.                                                                    !
    !--------------------------------------------------------------------------!
    ! # ARGUMENTS                                                              !
    ! <name> <in/out/inout> <arg descrption>                                   !
    ! example_arg1   out    An example argument with intent out.               !
    ! example_arg1   inout  An example argument with intent in.                !
    ! example_arg2   in     An example argument with intent in (optional)      !
    !--------------------------------------------------------------------------!
    ! # AUTHORS & CHANGELOG                                                    !
    ! Author(s):        <Author name(s)>                                       !
    ! Date of creation: <Date of creation>                                     !
    ! List of major changes:                                                   !
    ! <date> <change description> <author>                                     !
    !==========================================================================!
    use utils, only: utils_trace_in, utils_trace_out

    implicit none

    ! Parameters
    character(len=*), parameter :: myself = "example_priv_proc"

    ! Arguments
    real(kind=DP), intent(out)    :: example_arg1(:)
    integer, intent(inout)        :: example_arg2
    logical, intent(in), optional :: example_arg3

    ! Local variables
    logical :: example_arg3_loc

    call utils_trace_in(myself)

    ! Insert procedure code here

    !------ Example code for copying optional arg (delete if not needed) ------!
    if (present(example_arg3)) then
       example_arg3_loc = example_arg3
    else
       example_arg3_loc = .false.
    end if
    !--------------------------------------------------------------------------!

    call utils_trace_out(myself)

  end subroutine example_priv_proc

end module template
