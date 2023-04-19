! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                                 O N E T E P                                 !
!=============================================================================!
!                                                                             !
! The main ONETEP program.                                                    !
!                                                                             !
!-----------------------------------------------------------------------------!
!                                                                             !
!                 The ONETEP code is written and maintained by                !
!         Chris-Kriton Skylaris, Arash A. Mostofi, Nicholas D.M. Hine,        !
!                     Jacek Dziedzic and Peter D. Haynes                      !
!                                                                             !
!-----------------------------------------------------------------------------!
! Arash A. Mostofi, Version 0.01, 01/11/2004                                  !
! Originally written by Chris-Kriton Skylaris in 2000.                        !
! Improved and parallelised by Chris-Kriton Skylaris in November 2003.        !
! Improvements by Peter D. Haynes in 2004.                                    !
!-----------------------------------------------------------------------------!

program onetep
  use alloc_checker, only: alloc_init,alloc_report
  use anharmonic,   only: anharmonic_main
  use bibliography, only: bibliography_cite, bibliography_output
  use comms, only: comms_barrier, comms_bcast, comms_exit, &
       comms_init, comms_groups_init, pub_on_root, pub_root_proc_id
  use image_comms, only: image_comms_init, image_comms_finalize, &
       pub_on_imroots_root
  use constants, only: stdout, DP, file_maxsize, file_maxsize_rundat
  use couplings, only: couplings_exit
  use energy_and_force, only: energy_and_force_calculate, &
       energy_and_force_init_cell, energy_and_force_exit_cell
  use unit_test, only: unit_test_exec
  use eda_main, only: eda_task
  use esdf, only: esdf_init, esdf_close, esdf_errorout, &
       esdf_dump_input_file, esdf_stdout_dump, esdf_warnout
  use forcetest, only: forcetest_test_forces
  use geometry_optimiser, only: geometry_optimise
  use hubbard_init, only: hubbard_init_species_exit
  ! agrecokpt: include k_points_mod as well
  use k_points, only: kpoints_init
#ifdef ACCELRYS
  use license
#endif
  use md, only: md_main
  use model_type, only: MODEL
  use parallel_strategy, only: pub_par
  use phonon, only: phonon_main
  use rundat, only: get_rundat, rundat_exit, rundat_check_inputs, &
       rundat_threads_init, pub_all_tasks, pub_do_properties, &
       pub_rootname, pub_write_params, pub_write_xyz, &
       pub_task, pub_geom_continuation, pub_hubbard, &
       pub_anharmonic_calculate, pub_check_stack_size, pub_eda, &
       pub_cmplx_ngwfs, pub_num_images, pub_emft, &
       mix_dkn_type ! jd: Do not remove, used in #ifdef ACCELRYS
  use rundat_blocks, only: rundat_blocks_exec, rundat_blocks_exit
  use services, only: services_flush, services_write_xyz
  use spherical_wave, only: sw_exit
  use timer, only: timer_clock, timer_pin_subroutines_to_top
  use tssearch, only: tssearch_run
  use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
       utils_open_unit_check, utils_close_unit_check, utils_abort, &
       utils_check_stack_size, utils_init_array_checker, &
       utils_array_checker_report, utils_banner, &
       utils_filename_length_check, utils_update_rootname, &
       utils_feature_not_supported
  use vdwcorrection, only: vdwcorrection_override_dealloc
#ifdef ITC_TRACE
  use vt
#endif

  implicit none

  ! Local Variables
  type(MODEL), target                        :: mdl
  real(kind=DP), dimension(:,:), allocatable :: total_forces
  real(kind=DP)                              :: total_energy
  character(len=file_maxsize)                :: input_file,output_file
  integer                                    :: ierr,iunit
  integer                                    :: iidx
#ifdef ACCELRYS
  integer                                    :: num_licenses
  character (len=16)                         :: msversion="2016.0"
  character(len=20)                          :: prev_mix_dkn_type
#endif
  character (len=4)                          :: copyrightyear="2017"
  character (len=12)                         :: compilation_date
  logical                                    :: final_properties
  character(len=file_maxsize), dimension(:), allocatable, save :: block_args
  ! agrecokpt: single k-point in cartesian coordinates for testing
  real(kind=DP) :: kpt_cart(3)

  ! ----------------------------------------------------------------------------

  ! jd: Set up alloc checker
  call alloc_init

  ! vv: Call help feature
  call internal_help

  ! Initialise communications module
  call comms_init

  ! pdh: extract input file name from command line
  call internal_arguments

  ! jd: Let utils_mod know pub_rootname, so that it doesn't depend on rundat.
  ! This also deletes any potential .error_message files.
  call utils_update_rootname(pub_rootname, delete_err_file = .true.)

#ifdef __DATE__
    compilation_date=__DATE__
#ifdef _CRAYFTN
    ! jme: In CRAY Fortran, the date format is mm/dd/yy
    copyrightyear='20'//compilation_date(7:8)
#else
    copyrightyear=compilation_date(8:11)
#endif
#endif

#ifdef ACCELRYS
#ifdef MS_VERSION_STRING
    msversion=MS_VERSION_STRING
#endif
#endif

  ! pdh: print out something friendly so you know it's working...
  if (pub_on_root) call internal_welcome_banner


#ifdef ACCELRYS
  ! checkout licenses
  ! from 4.2 release there is no increase in number of checked out
  ! licenses with the number of processes; one license per job
  num_licenses = 1
  call license_checkout(num_licenses,ierr)
  if (ierr /= LIC_SUCCESS) call utils_abort('Error in ONETEP: &
          &failed to check out ONETEP licences')

  ! Set traps for signals
  call license_setup_traps()

#endif
#ifdef ITC_TRACE
  call vt_init_symbols()
#endif
  call utils_init_array_checker

  ! initialisation of the timer subroutine
  call timer_clock('total_time',0)

  ! cite standard ONETEP papers
  call bibliography_cite('JCP_special_2020')
  call bibliography_cite('ONETEP')

  ! pdh: read input file
  if (pub_on_root) call esdf_init(input_file,ierr)
  call comms_bcast(pub_root_proc_id,ierr)
  if (ierr /= 0) call utils_abort(&
       'Reading file "'//trim(input_file)//'" failed')

  ! cks: read run-time parameters from file
  call comms_barrier
  if (pub_on_root) write(stdout,'(3a)',advance='no') &
       'Reading parameters from file "',trim(input_file),'" ...'
  call get_rundat
  if (pub_on_root) write(stdout,'(a/)') '... done'

  ! jcap: if EMFT is defined alongside ScaLAPACK, abort gracefully
#ifdef SCALAPACK
  if (pub_emft) then
     call utils_feature_not_supported('EMFT with ScaLAPACK')
  end if
#endif

  ! kkbd: Setup ONETEP images if we're using them.
  call image_comms_init

  ! kkbd: If using images, get the others up to speed
  if(pub_num_images > 1)then
    if(pub_on_root .and. (.not.pub_on_imroots_root))then
      call esdf_init(input_file,ierr)
      if(ierr .ne. 0) &
        & call utils_abort("Error initializing esdf stuff on new image roots")
    end if
  end if

  ! jd: Ensure we have reasonable stack space.
  if(pub_check_stack_size) then
     call utils_check_stack_size
  end if

  ! jme: initialise OpenMP threading
  if (pub_on_root) write(stdout,'(a)') 'Checking processes and threads...'
  call rundat_threads_init
  call comms_barrier
  if (pub_on_root) write(stdout,'(a/)') '... done'

  ! jd: Reduce the timer overhead of commonly used, very short subroutines
  !     by pinning their names at the top of the timer name list.
  call timer_pin_subroutines_to_top

  ! ndmh: find first task out of pub_all_tasks list
  if (index(pub_all_tasks,' ')==0) pub_all_tasks = trim(pub_all_tasks)//' '
  iidx = index(pub_all_tasks,' ')
  pub_task = trim(pub_all_tasks(1:iidx))

  ! ndmh: perform input variable checks
  if (pub_on_root) write(stdout,'(a)',advance='no') 'Basic input checks...'
  call rundat_check_inputs
  call comms_barrier
  if (pub_on_root) write(stdout,'(a/)') '... done'
  call comms_barrier
  if (pub_on_root) write(stdout,'(3a)',advance='no') &
       'Reading geometry and species blocks from file "',trim(input_file), &
       '" ...'
  ! aam: Initialise elements and most of cell
  call rundat_blocks_exec(mdl)
  ! rc2013: maintain for now to avoid allocation issues in hubbard etc.
  pub_par=>mdl%regions(1)%par
  mdl%par=>mdl%regions(1)%par
  call comms_barrier
  if (pub_on_root) write(stdout,'(a)') '... done'

  ! cks: write coordinates in xyz file
  ! qoh: but not if we are doing a geometry optimistation continuation
  if (pub_write_xyz .and. .not. pub_geom_continuation) call services_write_xyz( &
       mdl%elements, pub_rootname, "Initial atomic coordinates")

  if (pub_on_root) then

     ! dump the input file to stdout
     call esdf_dump_input_file

     if (pub_write_params) then
        write(stdout,'(a)') REPEAT('-',80)
        write(stdout,'(a)') utils_banner('-', 'RUN-TIME PARAMETERS')
        write(stdout,'(a)') REPEAT('-',80)

        ! ndmh: write parameters to stdout
        call esdf_stdout_dump

        write(stdout,'(a)') REPEAT('-',80)
        write(stdout,'(a/)') REPEAT('-',80)
     else
        write(stdout,*)
     end if

     ! cks: Close ESDF subroutines
     call esdf_warnout
     call esdf_errorout
     call esdf_close

  end if

  call services_flush

  ! Initialise groups of cores in communications module
  call comms_groups_init()

  ! Allocate forces
  allocate(total_forces(1:3,1:mdl%nat),stat=ierr)
  call utils_alloc_check('ONETEP (main program)','total_forces', ierr)

  call services_flush

  call comms_barrier

  call energy_and_force_init_cell(mdl)

  call kpoints_init(kpt_cart, mdl%cell)

  do

     ! ndmh: find next task out of pub_all_tasks list
     if (index(pub_all_tasks,' ')==0) pub_all_tasks = trim(pub_all_tasks)//' '
     iidx = index(pub_all_tasks,' ')
     pub_task = trim(adjustl(pub_all_tasks(1:iidx)))

     ! Re-check inputs based on task
     call rundat_check_inputs

     ! aam: Determine what to do next based on value of "task"
     select case (pub_task)

     case ('SINGLEPOINT')

        ! agreco: added optional argument to specify if using complex NGWFs
        ! and single k-point
        call energy_and_force_calculate(total_energy,total_forces,mdl, &
            is_cmplx=pub_cmplx_ngwfs, kpt=kpt_cart)

     case ('PROPERTIES')

        if(pub_anharmonic_calculate) then
           call anharmonic_main(mdl%elements,mdl%cell)
        else

        call energy_and_force_calculate(total_energy,total_forces,mdl)

        end if

     case ('PHONON')

        call phonon_main(total_energy,total_forces,mdl)

     case ('GEOMETRYOPTIMIZATION')

        write(output_file,'(a,a)') trim(pub_rootname),'.geom'
        if (pub_do_properties) then
           final_properties = .true.
           pub_do_properties = .false.
        else
           final_properties = .false.
        end if
        ! agrecocmplx: allow for complex NGWFs to be used
        call geometry_optimise(total_energy,total_forces,mdl,output_file, &
             is_cmplx=pub_cmplx_ngwfs)
        if (final_properties) then
           pub_do_properties = .true.
           call energy_and_force_calculate(total_energy,total_forces,mdl, &
                is_cmplx=pub_cmplx_ngwfs)
        end if

     case ('MOLECULARDYNAMICS')

        write(output_file,'(a,a)') trim(pub_rootname),'.md'
        call md_main(total_energy,total_forces,mdl)
#ifdef ACCELRYS
        pub_task = 'SINGLEPOINT'         ! jd: Disable electronic_history
        prev_mix_dkn_type = mix_dkn_type !     for Biovia's e&f calculation
        mix_dkn_type = 'NONE'            !     below.
        call energy_and_force_calculate(total_energy,total_forces,mdl, &
             properties_only=.true.)
        pub_task = 'MOLECULARDYNAMICS'   ! jd: Restore state in case further
        mix_dkn_type = prev_mix_dkn_type !     tasks need it
#endif

     case ('TRANSITIONSTATESEARCH')

        write(output_file,'(a,a)') trim(pub_rootname),'.ts'
        call tssearch_run(mdl,output_file)
#ifdef ACCELRYS
        call energy_and_force_calculate(total_energy,total_forces,mdl, &
             properties_only=.true.)
#endif

     case ('FORCETEST')

        write(output_file,'(a,a)') trim(pub_rootname),'.forcetest'
        call forcetest_test_forces(total_energy,total_forces, &
             mdl,output_file)

     case ('HUBBARDSCF')

        write(output_file,'(a,a)') trim(pub_rootname),'.hubbardscf'
        call energy_and_force_calculate(total_energy,total_forces,mdl, &
             is_cmplx=pub_cmplx_ngwfs)

     case ('TDDFT')

        write(output_file,'(a,a)') trim(pub_rootname),'.tddft'
        call energy_and_force_calculate(total_energy,total_forces,mdl)

     case ('LR_TDDFT')

        write(output_file,'(a,a)') trim(pub_rootname), '.lr_tddft'
        call energy_and_force_calculate(total_energy,total_forces,mdl)

     case ('LR_PHONONS') ! gcc32
! gcc32 TEST
!       call utils_abort('LR_PHONONS IS NOT YET SUITABLE FOR RUNNING')
        write(output_file,'(a,a)') trim(pub_rootname), '.lr_phonons'
        call energy_and_force_calculate(total_energy,total_forces,mdl, &
             is_cmplx=pub_cmplx_ngwfs)

     case ('COND','PROPERTIES_COND')

        write(output_file,'(a,a)') trim(pub_rootname),'.cond'
        call energy_and_force_calculate(total_energy,total_forces,mdl, &
             is_cmplx=pub_cmplx_ngwfs)

     case ('COUPLINGS')

        call energy_and_force_calculate(total_energy,total_forces,mdl)

     case ('LUMOSEARCH')

        call energy_and_force_calculate(total_energy,total_forces,mdl)

     case ('UNIT_TEST')

        call unit_test_exec(mdl)

     case ('EDA', 'EDA_PREP')

        write(output_file,'(a,a)') trim(pub_rootname),'.eda'
        call eda_task(mdl)

     case default

        call utils_abort('Error in ONETEP: illegal value for task: '//trim(pub_task))

     end select

     call timer_clock('onetep.F90 main loop imbalance',1)
     call comms_barrier
     call timer_clock('onetep.F90 main loop imbalance',2)

     ! trim down pub_all_tasks list and exit if no more tasks are left
     pub_all_tasks = pub_all_tasks(iidx+1:)
     if (len(trim(pub_all_tasks))==0) exit

  end do

  if (.not. pub_eda) call energy_and_force_exit_cell(mdl)

  call rundat_exit
  call couplings_exit
  call vdwcorrection_override_dealloc
  call sw_exit

  ! Deallocate forces and elements
  deallocate(total_forces,stat=ierr)
  call utils_dealloc_check('ONETEP (main program)','total_forces', ierr)

  if (pub_hubbard) then
     ! ddor: deallocate hubbard_species type array in DFT+U
     call hubbard_init_species_exit
  end if

  ! deallocate arrays read in by rundat_blocks
  call rundat_blocks_exit(mdl)

  ! write the bibliography
  write(output_file,'(a80)') trim(pub_rootname)//'.bib'
  output_file = adjustl(output_file)
  iunit = utils_unit()
  if (pub_on_root)then
     open(unit=iunit,file=output_file,iostat=ierr)
     call utils_open_unit_check('ONETEP (main program)',output_file,ierr)
  end if
  call bibliography_output(iunit)
  if (pub_on_root) then
     close(iunit,iostat=ierr)
     call utils_close_unit_check('ONETEP (main program)',output_file,ierr)
  end if

  ! jd: Measure how long a pair of timer_clock calls takes, to get an idea of
  !     how big the overhead is.
  call timer_clock('#timer_clock overhead for 1000 calls',1)
  do iidx = 1, 1000
     call timer_clock('#ignore_me',1)
     call timer_clock('#ignore_me',2)
  end do
  call timer_clock('#timer_clock overhead for 1000 calls',2)

  call timer_clock('onetep.F90 final imbalance',1)
  call comms_barrier
  call timer_clock('onetep.F90 final imbalance',2)

  ! shutdown of the timer subroutine
  call timer_clock('total_time',3)

  call utils_array_checker_report

  call alloc_report

  ! say goodbye
  if (pub_on_root) call internal_farewell

#ifdef ACCELRYS
  ! check the license back in
  call license_checkin(ierr)
#endif

  ! close down images if we're using them
  call image_comms_finalize

  ! clean up comms module
  call comms_exit

contains

  subroutine internal_help

    use esdf, only: esdf_help
    use utils, only: utils_assert
    implicit none

    ! Local variables
    integer :: num_args
    logical :: halt
    integer :: iargs,ierr
    character(file_maxsize):: arg(1:3)
    character(file_maxsize) :: helpflags(1:3)
    character(file_maxsize) :: searchflags(1:2)

    helpflags   =["-h    ","--help","-help "]
    searchflags =["-s      ","--search"]

    num_args = command_argument_count()

    allocate(block_args(0:num_args),stat=ierr)
    call utils_alloc_check('ONETEP (main program)','block_args',ierr)

    do iargs = 0, num_args
       call get_command_argument(iargs, value=block_args(iargs), status=ierr)
       call utils_assert(ierr == 0, 'ERROR in onetep (internal_help): &
            &get_command_argument returned an error status for one &
            &of the arguments. (argument, status) =', iargs, ierr)
    end do

    halt = .false.

    arg(1:3)="-"

    ! Copy to our local cmd_args variable
    if(num_args==0) then

    elseif(num_args==1) then
       arg(1)=block_args(1)
    elseif(num_args==2) then
       arg(1:2)=block_args(1:2)
    elseif(num_args>=3) then
       arg(1:3)=block_args(1:3)
    endif

    ! --help search <searchword>
    if(any(helpflags==arg(1)).and. (trim(adjustl(arg(2)))=="SEARCH" .or. &
         trim(adjustl(arg(2)))=='search') .and.arg(3)(1:1)/="-") then
       call esdf_help('search', arg(3))
       halt=.true.
       ! --help <keyowrd>
    elseif(any(helpflags==arg(1)).and. arg(2)(1:1)/="-" .and. &
         arg(2)/="SEARCH" .and. arg(2)/="search") then
       call esdf_help(arg(2), arg(2))
       halt=.true.
       ! --search <searchword>
    elseif(any(searchflags==arg(1)).and. arg(2)(1:1)/="-" .and. arg(3)=="-")&
         then
       call esdf_help('search', arg(2))
       halt=.true.
       ! <inputname>
    elseif(arg(1)(1:1)/="-")then ! It must be an inputname
       halt=.false.
    elseif(arg(1)=="-" .and. arg(2)=="-" .and. arg(3)=="-") then
       halt=.false.
       ! -help
    elseif(any(helpflags==arg(1)) .and. arg(2)(1:1)=="-" .and. &
         arg(3)(1:1)=="-")then
       halt=.true.
       call esdf_help('help','')
       ! <none of the above!>
    elseif(arg(1)=="-search") then
       write(6,'(a,a,a)') 'Did you mean --search ', trim(adjustl(arg(2))),' ?'
       halt=.true.
    elseif(arg(1)=="--h") then
       write(6,'(a,a,1x,a,a)') 'Did you mean [-help|--help|-h] ', &
            trim(adjustl(arg(2))), trim(adjustl(arg(3))),' ?'
       halt=.true.
    else ! Couldn't make any sense of it. Call help.
       call esdf_help('help',arg(2))
       halt=.true.
    endif

    if (halt) then
       stop!!cannot use utils_abort() here, too early in the initialisation.
    end if

  end subroutine internal_help

  subroutine internal_arguments

    !===================================================!
    ! This subroutine deals with command line arguments !
    !===================================================!

    implicit none

    ! Local variables
    integer :: num_args, len_root

    if (pub_on_root) then
       num_args = command_argument_count()
    end if
    call comms_bcast(pub_root_proc_id,num_args)

    ! JCW: Each block_args element is a character array of length file_maxsize
    ! JCW: input_file is a character array of length file_maxsize
    if (num_args == 0) then
       input_file = 'onetep.dat'
    else
       input_file = block_args(1)
       call comms_bcast(pub_root_proc_id,input_file)
    end if

    ! JCW: Check if length of input_file exceeds string length which causes
    ! JCW: failure in get_rundat. If so, abort with informative error message.
    ! JCW: Set file_output=.false., since we have not yet established whether
    ! JCW: the filename length is acceptable and have not initialized utils
    ! JCW: module's version of pub_rootname.
    call utils_filename_length_check(trim(input_file),&
         "onetep, internal_arguments",file_maxsize_rundat,file_output=.false.)

    len_root = index(input_file,'.dat') - 1
    ! JCW: pub_rootname is a character array of length 80 (!)
    ! JCW: input_file could therefore be 80+4 in length, but in many places
    ! JCW: filenames of length 80 are created as pub_rootname+extension.
    ! JCW: It is (slightly) safer to make 80 the limit for input_file.
    ! JCW: In future, we should standardize allowed filename lengths across
    ! JCW: the code to avoid danger where files approach 80 characters in length.
    if (len_root > 0) then
       pub_rootname = input_file(1:len_root)
    else
       pub_rootname = input_file
    end if
    write(input_file,'(a80)') trim(pub_rootname)//'.dat'
    input_file = adjustl(input_file)
#ifdef ACCELRYS
    write(output_file,'(a80)') trim(pub_rootname)//'.onetep'
    ! JCW: output_file is a character array of length file_maxsize
    output_file = adjustl(output_file)
    if (pub_on_root)then
       open (unit=stdout,file=output_file)
    end if
#endif

    deallocate(block_args,stat=ierr)
    call utils_dealloc_check('ONETEP (main program)','block_args',ierr)

  end subroutine internal_arguments

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine internal_welcome_banner

    !========================================================!
    ! Banner printed at the beginning of ONETEP calculation. !
    !--------------------------------------------------------!
    ! Updated by Chris-Kriton Skylaris on 14/10/2004.        !
    !========================================================!

    implicit none

    character(len=8)  :: date
    character(len=10) :: time
    character(len=5)  :: zone

#ifdef ACCELRYS
    character(len=128) :: mathlib_version_string
    integer :: mkl_version_end
#endif

    call date_and_time(date,time,zone)

    write(stdout,*) '+---------------------------------------------------------------+'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|        ####### #     # ####### ####### ####### ######         |'
    write(stdout,*) '|        #     # ##    # #          #    #       #     #        |'
    write(stdout,*) '|        #     # # #   # #          #    #       #     #        |'
    write(stdout,*) '|        #     # #  #  # #####      #    #####   ######         |'
    write(stdout,*) '|        #     # #   # # #          #    #       #              |'
    write(stdout,*) '|        #     # #    ## #          #    #       #              |'
    write(stdout,*) '|        ####### #     # #######    #    ####### #              |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|        Linear-Scaling Ab Initio Total Energy Program          |'
    write(stdout,*) '|                                                               |'
#ifdef ACCELRYS
    write(stdout,*) '|                Materials Studio version '//msversion//'      |'
#else
    write(stdout,*) '|          Release for academic collaborators of ODG            |'
    write(stdout,*) '|                                             Version 6.1.12.0  |'
#endif
    write(stdout,*) '|                                                               |'
    write(stdout,*) '+---------------------------------------------------------------+'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  Authors:                                                     |'
    write(stdout,*) '|  Jacek Dziedzic, Peter D. Haynes, Nicholas D. M. Hine,        |'
    write(stdout,*) '|  Arash A. Mostofi, Mike C. Payne and Chris-Kriton Skylaris    |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  Contributors:                                                |'
    write(stdout,*) '|  J. Aarons, L. Andrinopoulos, P. W. Avraam, R. A. Bell,       |'
    write(stdout,*) '|  A. Bhandari, G. A. Bramley, R. J. Charlton, S. J. Clark,     |'
    write(stdout,*) '|  R. J. Clements, G. C. Constantinescu, F. Corsetti,           |'
    write(stdout,*) '|  N. Corsini, O. Dieguez, S. M. M. Dubois, K. K. B. Duff,      |'
    write(stdout,*) '|  J. M. Escartin, A. Greco, H. Helal, Q. O. Hill, L. P. Lee,   |'
    write(stdout,*) '|  J.-H. Li, E. B. Linscott, G. Moynihan, D. D. O`Regan,        |'
    write(stdout,*) '|  O. K. Okan, M. J. S. Phipps, C. J. Pickard,                  |'
    write(stdout,*) '|  J. C. A. Prentice, M. I. J. Probert, L. E. Ratcliff,         |'
    write(stdout,*) '|  M. Robinson, A. Ruiz Serrano, J. S. Spencer, E. W. Tait,     |'
    write(stdout,*) '|  G. Teobaldi, D. Turban, V. Vitale, K. A. Wilkinson,          |'
    write(stdout,*) '|  C. Weber, J. C. Womack, N. Yeung and T. J. Zuehlsdorff       |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|                                   Copyright (c) 2004-'//&
         & copyrightyear//'     |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  Please cite:                                                 |'
    write(stdout,*) '|  "The ONETEP linear-scaling density functional theory program"|'
    write(stdout,*) '|   J. Chem. Phys. 152(17) 174111 (2020)                        |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|          in all publications arising from your use of ONETEP. |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '+---------------------------------------------------------------+'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|   ONETEP is based on developments described in the following  |'
    write(stdout,*) '|   publications:                                               |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  "Nonorthogonal generalized Wannier function pseudopotential  |'
    write(stdout,*) '|   plane-wave method".                                         |'
    write(stdout,*) '|   C.-K. Skylaris, A. A. Mostofi, P. D. Haynes, O. Dieguez,    |'
    write(stdout,*) '|   M. C. Payne.                                                |'
    write(stdout,*) '|   Phys. Rev. B 66 035119 (2002).                              |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  "Preconditioned iterative minimization for linear-scaling    |'
    write(stdout,*) '|   electronic structure calculations".                         |'
    write(stdout,*) '|   A. A. Mostofi, P. D. Haynes, C.-K. Skylaris, M. C. Payne.   |'
    write(stdout,*) '|   J. Chem. Phys. 119(17), pp.8842-8848 (2003).                |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  "Total-energy calculations on a real space grid with         |'
    write(stdout,*) '|   localized functions and a plane-wave basis".                |'
    write(stdout,*) '|   A. A. Mostofi, C.-K. Skylaris, P. D. Haynes, M. C. Payne.   |'
    write(stdout,*) '|   Comput. Phys. Commun. 147, pp.788-802 (2002).               |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  "Accurate kinetic energy evaluation in electronic structure  |'
    write(stdout,*) '|   calculations with localized functions on real space grids"  |'
    write(stdout,*) '|   C.-K. Skylaris, A. A. Mostofi, P. D. Haynes, C. J. Pickard, |'
    write(stdout,*) '|   M. C. Payne.                                                |'
    write(stdout,*) '|   Comput. Phys. Commun. 140, pp.315-322 (2001).               |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  "Accurate ionic forces and geometry optimization in linear-  |'
    write(stdout,*) '|   scaling density-functional theory with local orbitals"      |'
    write(stdout,*) '|   N. D. M. Hine, M. Robinson, P. D. Haynes, C.-K. Skylaris,   |'
    write(stdout,*) '|   M. C. Payne, and A. A. Mostofi.                             |'
    write(stdout,*) '|   Phys. Rev. B 83 195102 (2011).                              |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '+---------------------------------------------------------------+'

#ifdef ACCELRYS
#ifdef __DATE__
#ifndef PLATFORM
#define PLATFORM "MS Windows"
#endif
#ifdef debug
#define DEBUG " DEBUG "
#else
#define DEBUG " "
#endif
    write(stdout,*)
    write(stdout,*) "This",DEBUG, "version was compiled for ",PLATFORM, &
             & " on ", __DATE__
    write(stdout,*)
    call mkl_get_version_string(mathlib_version_string)
    mkl_version_end = index(mathlib_version_string,'Product ') - 1
    write(stdout,*)  mathlib_version_string(1:mkl_version_end)

    write(stdout,*)
#endif
#endif

    write(stdout,'(/a,2(a2,a1),a4,1x,a2,a1,a2,a2,a5,a1)') 'Job started: ', &
         date(7:8),'-',date(5:6),'-',date(1:4),time(1:2),':',time(3:4),&
         ' (',zone,')'

    write(stdout,*)

    call services_flush

  end subroutine internal_welcome_banner

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine internal_farewell

    implicit none

    character(len=8)  :: date
    character(len=10) :: time
    character(len=5)  :: zone

    call date_and_time(date,time,zone)

    write(stdout,'(/a,2(a2,a1),a4,1x,a2,a1,a2,a2,a5,a1/)') 'Job completed: ', &
         date(7:8),'-',date(5:6),'-',date(1:4),time(1:2),':',time(3:4),&
         ' (',zone,')'

    call services_flush

  end subroutine internal_farewell

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end program onetep
