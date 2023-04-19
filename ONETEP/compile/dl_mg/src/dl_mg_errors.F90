!> data types and subroutines for error handling
module dl_mg_errors
  use dl_mg_params
  implicit none
  public

  !> Warning : If a new error code is added nerror must be increased by one and
  !! a corressponding error message must be added in error_list array.
  !! This is awkward but I cannot think of a better solution.
  !! error codes are defined in dl_mg_common_params.inc
  
  integer, parameter :: nerror = 34

  type errm
     character(len=DL_MG_MAX_ERROR_STRING) m
  end type errm

  type(errm), parameter :: error_list(nerror) = (/ &
       !! [error-codes]
       ! DL_MG_ERR_UNSPECIFIED = 1
       errm("Unspecified error:"), &
       ! Dl_MG_ERR_INIT = 2
       errm("dl_mg_init called before dl_mg_free"), &
       ! DL_MG_ERR_NEUMANM_BC = 3
       errm("Neumann boundary condition is not implemented"), &
       !DL_MG_ERR_UNKNOWN_BC = 4
       errm("dl_mg_init: wrong boundary condition specification &
       &in input. It must be one of: DL_MG_BC_PERIODIC, DL_MG_BC_DIRICHLET"),&
       !DL_MG_ERR_MPI_TOPO = 5
       errm("input communicator does not have a cartesian topology"), &
       !DL_MG_ERR_MPI_COMM_DUP = 6
       errm("duplicate communicator does not have a &
       &cartesian topology"), &
       !DL_MG_ERR_MPI_TOPO_BC = 7
       errm("boundary conditions inconsistent with MPI topology"),&
       ! DL_MG_ERR_GRID_EMPTY = 8
       errm("empty inner grid in set_mg_grids"),&
       ! DL_MG_ERR_PROLONG_MAP = 9
       errm("build_prolong_map: inactive site left behind"), &
       ! DL_MG_ERR_DEBYE = 10
       errm("dl_mg_init_nonlin: Debye length squared is negative"), &
       !DL_MG_ERR_NOINIT_POISSON = 11
       errm("dl_mg_solver_poisson is called before dl_mg_init"), &
       ! DL_MG_ERR_NOINIT_PBE = 12
       errm("dl_mg_solver_pbe called before initialisation subroutines"), &
       ! DL_MG_ERR_DERPOT = 13
       errm("dl_mg_solver_pbe: not an error, functionality added"), &
       ! DL_MG_ERR_MOD_BC = 14
       errm("Dirichlet boundary values were modified by the solver"), &
       ! DL_MG_ERR_RES_RATIO = 15
       errm("Consecutive residuals ratio is too large.&
       & It is possible that the round-off limit was reached.&
       & Check the residual ratio distribution histogram (enabled with DL_MG_RES_HIST)."), &
       ! DL_MG_ERR_NITER = 16
       errm("Computation failed to converge after the prescribed number of iterations"), &
       ! DL_MG_ERR_NEWTON_TYPE = 17
       errm("wrong eq_type in Newton's MG"), &
       ! DL_MG_ERR_NEWTON_DAMP = 18
       errm("compute_damping_parameter: failed to find a damping parameter for&
       & Newton method"), &
       ! DL_MG_ERR_NEWTON_DAMP_DERFUN = 19
       errm("compute_damping_parameter: functional derivative is non-negative"), &
       ! DL_MG_ERR_ASSERTION          = 20    ! JCW
       errm("assertion error"), &
       ! DL_MG_ERR_ALLOCATION         = 21    ! JCW
       errm("error in memory allocation"), &
       ! DL_MG_ERR_DEALLOCATION       = 22    ! JCW
       errm("error in memory deallocation"), &
       ! DL_MG_ERR_IO                 = 23    ! JCW
       errm("error in file I/O"),&
       ! DL_MG_ERR_GRID_GEOM          = 24    ! JCW
       errm("error in setting up grid geometry"), &
       ! DL_MG_ERR_MPI_FAIL           = 25    ! JCW
       errm("MPI routine returned a non-zero error code"), &
       ! DL_MG_ERR_EQTYPE = 26
       errm("Wrong equation type"), &
       ! DL_MG_ERR_FD_ORDER           = 27    ! JCW
       errm("requested FD order not supported"), &
       ! DL_MG_ERR_DEFCO_ITER         = 28    ! JCW
       errm("convergence was not achieved within the maximum number of &
       &defect correction iterations"), &
       ! DL_MG_ERR_DEFCO_UNPHYSICAL   = 29    ! JCW
       errm("unphysical solution obtained in defect correction"), &
       ! DL_MG_ERR_DEFCO_DAMPING      = 30    ! JCW
       errm("error damping for defect correction failed"), &
       ! DL_MG_ERR_DEFCO_DERPOT       = 31    ! JCW
       errm("der_pot should not be passed to the solver when the defect &
       &correction is used, since the potential at which the derivative is &
       &evaluated is generated internally&
       & (within the defect correction loop)"), &
       ! DL_MG_ERR_NOMPI              = 32
       errm("Using serial build code in MPI run?"), &
       ! DL_MG_ERR_CION_ZERO          = 33
       errm("Ion concentration is zero in PBE solver. &
       &Use Poisson or linearised PBE"), &
       ! DL_MG_ERR_CION_NEG           = 34
       errm("Shifted ion concentration is negative.") &
       !! [error-codes]
       /)

  logical, save :: abort_on_failure = .true. !< control wether
  !! to abort or return control to the calling up when error conditions are detected.
  !! Control is return to the calling routine only
  !! if handle_error is called with ierror argument
  integer, save :: first_error = -1 !< keeps the first error code reported by handle_error()
  !! this useful because the ierror code could be changed as dlmg is quitting and the control
  !! is returned to the higher level subroutines
  !! NB the error code are positive integers

  private error_list, nerror, error_abort, errm

contains

  subroutine handle_error(err_code, ierror, msg)
    implicit none
    integer, intent(in) :: err_code !< error code
    integer, optional, intent(out) :: ierror !< value to be return to the caller
    character(len=*), optional, intent(in) :: msg !< supplementary message
                                                  !! used only if the execution aborts
                                                  !! on error

    if (abort_on_failure) then
       if (err_code /= DL_MG_SUCCESS) then
          call error_abort(err_code, msg)
       endif
    else
       if (.not. present(ierror)) then
          write(0,*) "!!! Error_returns is true but this comes from a call&
               & without ierror argument."
          call error_abort(err_code, msg)
       else
          ierror = err_code
          if (first_error == -1) then
             first_error = ierror
          end if
       end if
    end if

  end subroutine handle_error

  
  subroutine error_abort(error_code, msg)
    use dl_mg_mpi_header
    use dl_mg_common_data, only : mg_comm, report_unit
    implicit none
    integer, intent(in) :: error_code
    character(len=*), optional, intent(in) :: msg

    integer ierr
    character(len=DL_MG_MAX_ERROR_STRING) s

    call get_errmsg(error_code, s)
    if (present(msg)) then
       if (error_code == DL_MG_ERR_UNSPECIFIED .and. &
            len_trim(msg) > 0) then
          ! no need to print "unspecified error" if a message is provided,
          ! it could be connfusing
          s=msg
       else
          s = trim(s)//" "//msg
       end if
    endif

    write(report_unit,*) "DL_MG fatal error: "//s
    write(0,*) "DL_MG fatal error: "//s
#ifdef USE_BACKTRACE
    call backtrace
#endif

#ifdef MPI
    call MPI_abort(mg_comm, error_code,ierr)
#else
    stop
#endif

   end subroutine error_abort


   subroutine get_errmsg(err_code, msg)
     implicit none
     integer, intent(in) :: err_code
     character(len=DL_MG_MAX_ERROR_STRING), intent(out) :: msg

     integer i
     character(len=15) c

     if ( err_code < 0 .or. err_code > nerror ) then
        write(c,*) err_code
        msg=trim(adjustl(c))//": wrong error code?"
     else if (err_code == DL_MG_SUCCESS) then
        msg = "successful execution"
     else
        msg = error_list(err_code)%m
     end if
   end subroutine get_errmsg

   !> can be used to trigger assertion failure.
   !! useful for testing error handling
   function check_assertion(assertion, err_code)
     implicit none

     logical, intent(in) :: assertion !< logical expression whose value is
                                      !! is returned
     integer, optional, intent(in) :: err_code !< used for testing the error reporting;
                                               !! if missing nothing happens

     logical check_assertion

#ifdef DL_MG_TEST_ERR_CODES
     integer  len, tvals(8)
     character(len=63) val
     integer, save :: err_code_to_test=0
     logical, save :: firstcall=.true.
#endif

     check_assertion = assertion

#ifdef DL_MG_TEST_ERR_CODES
     ! this is used to simulate an error at a selected point in
     ! the code for error testing
     if (firstcall) then
        firstcall=.false.
        call get_environment_variable("DL_MG_TEST_ERR_CODE", value=val, &
             length = len)
        if (len > 0) then
           read(val, *) err_code_to_test
        else
           write(0,*) "WARNING: DL_MG_TEST_ERR_CODE length is 0. !!!"
        endif
     end if
        ! if err_code_to_test < 0 pick at random a failure point
        ! in this case -err_code_to_test is the percentual
        ! probability to fail
     if (err_code_to_test < 0) then
        call date_and_time(values=tvals)
        if (mod(tvals(8),100) < -err_code_to_test ) check_assertion = .false.
     else if ( present(err_code)) then
        if (err_code == err_code_to_test)  check_assertion = .false.
     endif
#endif
   end function check_assertion


end module dl_mg_errors
