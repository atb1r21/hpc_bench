! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                                Utils module                                 !
!=============================================================================!
! The module for miscellaneous utilities in ONETEP. Routines in this module   !
! can only depend on the constants, comms and io modules.                     !
!-----------------------------------------------------------------------------!
!                                                                             !
!                 The ONETEP code is written and maintained by                !
!          Chris-Kriton Skylaris, Arash A. Mostofi and Peter D. Haynes        !
!                                                                             !
!-----------------------------------------------------------------------------!
! This file written by Peter Haynes, 22 March 2005.                           !
! Subsequent modifications by Peter Haynes, Jacek Dziedzic,                   !
! Alvaro Ruiz Serrano, and Jose M Escartin.                                   !
!=============================================================================!

module utils

  use constants, only: DP, LONG, MAX_WORD_LENGTH

  implicit none

  private

! jme: if we know that the compiler is buggy, abort compilation here
#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER == 1800)
#if defined(__INTEL_COMPILER_UPDATE) && (__INTEL_COMPILER_UPDATE < 5)
! Bug (Intel issue #03143573): intent(out) pointers mess up with the target
! of the pointer at entry.  Appeared in 18.0.0, was fixed in 18.0.2.
! However, other suspected bugs (segfaults in for_dealloc_all_nocheck)
! seem to persist in 18.0.2 and 18.0.3, while 18.0.5 does not appear to
! have this issue. We thus exclude 18.0.[1234].
#error "ONETEP compilation ERROR: Intel Compiler versions 18.0.0-18.0.4 are not supported due to known bugs."
#endif
#endif

  ! jd: Cannot use 'pub_debug' here, utils cannot depend on rundat
#ifdef DEBUG
  logical, parameter :: loc_debug = .True.
#else
  logical, parameter :: loc_debug = .False.
#endif

#ifdef DEBUG_ARRAYS
  logical, parameter :: loc_debug_arrays = .True.
#else
  logical, parameter :: loc_debug_arrays = .False.
#endif
  logical :: debug_arrays_alloc = .False.

  type array_info
     character(len=36) :: name
     character(len=72) :: routine
     integer :: allocated
     integer :: max_allocated
  end type array_info

  ! rc2013: need a lot more arrays for embedding
  integer, parameter :: max_arrays = 20000
  integer :: num_arrays = 0
  type(array_info), allocatable :: arrays(:)

  public :: utils_init_array_checker
  public :: utils_array_checker_report

  interface utils_heapsort
     module procedure utils_heapsort_integer
     module procedure utils_heapsort_real
  end interface

  interface utils_assert
     module procedure utils_assert_with_string
     module procedure utils_assert_with_integers
     module procedure utils_assert_with_reals
     module procedure utils_assert_with_logicals
  end interface

  interface utils_use_var
     module procedure utils_use_var_character
     module procedure utils_use_var_integer
     module procedure utils_use_var_logical
     module procedure utils_use_var_real
  end interface

  interface utils_qc_print
     module procedure utils_qc_print_integer
     module procedure utils_qc_print_real
     module procedure utils_qc_print_string
  end interface

  interface utils_devel_code
     module procedure utils_devel_code_logical
     module procedure utils_devel_code_integer
     module procedure utils_devel_code_real
  end interface

  interface utils_compare_vectors
     module procedure utils_compare_vectors_integer
     module procedure utils_compare_vectors_real
     module procedure utils_compare_vectors_complex
  end interface utils_compare_vectors

  interface utils_sanity_check
     module procedure utils_sanity_check_real_0
     module procedure utils_sanity_check_real_1
     module procedure utils_sanity_check_real_2
     module procedure utils_sanity_check_real_3
  end interface

  public :: utils_alloc_check
  public :: utils_dealloc_check
  public :: utils_flush
  public :: utils_unit
  public :: utils_binary_copy
  public :: utils_heapsort
  public :: utils_open_unit_check
  public :: utils_close_unit_check
  public :: utils_read_check
  public :: utils_write_check
  public :: utils_assert
  public :: utils_abort
  public :: utils_trace_in
  public :: utils_trace_out
  public :: utils_sanity_check
  public :: utils_erf
  public :: utils_erfc
  public :: utils_custom_erfc
  public :: utils_dump_array1D_to_file
  public :: utils_dump_array2D_to_file
  public :: utils_dump_array3D_to_file
  public :: utils_is_finite
  public :: utils_isnan
  public :: utils_use_var
  public :: utils_qc_print
  public :: utils_write_csv
  public :: utils_read_csv
  public :: utils_devel_code
  public :: utils_busy_wait
  public :: utils_check_stack_size
  public :: utils_report_memory_estimate
  public :: utils_fnv1_hash
  public :: utils_fnv1_hash_file
  public :: utils_rand_gen
  public :: utils_feature_not_supported
  public :: utils_compare_vectors
  public :: utils_flushed_string_output
  public :: utils_logical_to_real
  public :: utils_character_to_real
  public :: utils_real_to_logical
  public :: utils_real_to_character
  public :: utils_long_int_to_str
  public :: utils_int_to_str
  public :: utils_real_to_str
  public :: utils_point_to_str
  public :: utils_str_to_int
  public :: utils_strip_char_from_string
  public :: utils_lock_ping
  public :: utils_lock_block
  public :: utils_nth_word_of
  public :: utils_banner
  public :: utils_erase_file
  public :: utils_filename_length_check
!  public :: utils_stream_copy
  public :: utils_safe_int
  public :: utils_safe_nint
  public :: utils_parse_bc
  public :: utils_postfix_to_ngwf_set_name
  public :: utils_update_rootname

#ifdef TRACE
  integer, private, save :: nindent = 0
#endif

  ! jd: Do not change this string
  character(len=20), public, save :: utils_random_seed = '_%++%_7cd70ff7_%++%_'

  ! jd: A copy of rundat::pub_rootname, so that we don't depend on rundat
  ! ebl: to modify this copy outside of utils_mod, use utils_update_rootname
  character(len=80), private, save :: pub_rootname_local

contains

  subroutine utils_init_array_checker()

    implicit none
    integer :: ierr

    if (loc_debug_arrays) then
       allocate(arrays(max_arrays), stat=ierr)
       if (ierr /=0) then
          call utils_abort('Error in utils_init_array_checker: allocating &
               &arrays failed with code ', ierr)
       end if
       debug_arrays_alloc = .True.
       num_arrays = 0
       arrays(:)%name = ''
       arrays(:)%routine = ''
       arrays(:)%allocated = 0
       arrays(:)%max_allocated = 0
    end if

  end subroutine utils_init_array_checker

  subroutine utils_array_checker_report()
  !----------------------------------------------------------------------------!
  ! @docme                                                                     !
  !----------------------------------------------------------------------------!
  ! Caveats:                                                                   !
  !   jd: The array checker is limited to thread 0 (cf. utils_alloc_check(),   !
  !       utils_dealloc_check(). Thus, allocations that happen inside OMP      !
  !       regions on shared data structures protected with a critical will be  !
  !       missed by the array checker on all threads except 0. This will lead  !
  !       to the array checker complaining if these datastructures are         !
  !       subsequently deallocated on thread 0 -- the checker will think there !
  !       are more deallocs than allocs. This happens for hash tables in HFx.  !
  !----------------------------------------------------------------------------!

    use comms, only: pub_on_root
    use constants, only: stdout
    implicit none

    integer :: iarray, ierr

    if ( loc_debug_arrays .and. pub_on_root) then

       write(stdout,'(a)') utils_banner('=', 'ARRAY (DE)ALLOCATION CHECKS - REPORT')
       if ( ANY(arrays(1:num_arrays)%allocated /= 0 ) ) then
          write(stdout,'(a)')  '       Allocated'
          write(stdout,'(3a)') '   Id Final   Max Array tag', &
               repeat(' ',27), 'Subprogram'
          write(stdout,'(a)') repeat('-',80)
          do iarray=1,num_arrays
             if (arrays(iarray)%allocated /= 0) &
                  write(stdout,'(i5,i6,i6,1x,a36,a64)') iarray, &
                       arrays(iarray)%allocated, arrays(iarray)%max_allocated, &
                       arrays(iarray)%name, arrays(iarray)%routine
          end do
       else
          write(stdout,'(a)') 'All allocation and deallocation tag counts match.'
       end if
       write(stdout,'(a)') repeat('=',80)

       deallocate(arrays, stat=ierr)
       if (ierr /=0) then
          call utils_abort('Error in utils_array_checker_report: deallocating &
               &arrays failed with code ', ierr)
       end if
       debug_arrays_alloc = .False.

    end if

  end subroutine utils_array_checker_report

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !============================================================================!
  ! This subroutine checks whether an array has been allocated successfully.   !
  ! On failure it writes an error message and aborts.                          !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   routine (in) : Name of the subroutine which called this routine          !
  !   array   (in) : Name of the array allocated                               !
  !   status  (in) : Status flag returned by allocate                          !
  !----------------------------------------------------------------------------!
  ! Caveats:                                                                   !
  !   jd: The array checker part is limited to thread 0. Thus, allocations that!
  !       happen inside OMP regions on shared data structures protected with a !
  !       critical will be missed by the array checker on all threads except 0.!
  !       This will lead to the array checker complaining if these datastruc-  !
  !       tures are subsequently deallocated on thread 0 -- the checker will   !
  !       think there are more deallocs than allocs. This happens for hash     !
  !       tables in HFx.                                                       !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2005.                                       !
  !============================================================================!

  subroutine utils_alloc_check(routine,array,status)

    implicit none

    ! Arguments
    character(len=*), intent(in) :: routine    ! Routine name
    character(len=*), intent(in) :: array      ! Array name
    integer, intent(in)          :: status     ! Status flag

    ! Locals
    integer :: iarray
!$  integer, external :: omp_get_thread_num

    if (status /= 0) then
       call utils_abort('Error in '//trim(routine)//': allocating '//&
            trim(array)//' failed with code ',status)
    end if

    if (.not. (loc_debug_arrays .and. debug_arrays_alloc) ) return
!$  if (omp_get_thread_num()==0) then
    do iarray=1,num_arrays
       ! ndmh: if current array name and parent has been used before
       if ((array == arrays(iarray)%name) .and. &
            (routine == arrays(iarray)%routine)) then
          ! ndmh: increment allocated count for this array
          arrays(iarray)%allocated = arrays(iarray)%allocated + 1
          arrays(iarray)%max_allocated = max(arrays(iarray)%allocated, &
               arrays(iarray)%max_allocated)
          return
       end if
    end do

    ! ndmh: make sure you don't exceed allocated memory
    call utils_assert(num_arrays /= max_arrays, 'Error in utils_alloc_check(): &
         &maximum number of arrays exceeded: ', max_arrays)

    num_arrays = num_arrays + 1
    arrays(num_arrays)%name = array
    arrays(num_arrays)%routine = routine
    arrays(num_arrays)%allocated = 1
    arrays(iarray)%max_allocated = max(arrays(iarray)%allocated, &
         arrays(iarray)%max_allocated)

!$  end if

  end subroutine utils_alloc_check

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !============================================================================!
  ! This subroutine checks whether an array has been deallocated successfully. !
  ! On failure it writes an error message and aborts.                          !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   routine (in) : Name of the subroutine which called this routine          !
  !   array   (in) : Name of the array deallocated                             !
  !   status  (in) : Status flag returned by deallocate                        !
  !   allocated_in (in, opt): Name of the subroutine where the array was       !
  !                           allocated (helps the array checker find a match).!
  !----------------------------------------------------------------------------!
  ! Caveats:                                                                   !
  !   jd: The array checker part is limited to thread 0. Thus, allocations that!
  !       happen inside OMP regions on shared data structures protected with a !
  !       critical will be missed by the array checker on all threads except 0.!
  !       This will lead to the array checker complaining if these datastruc-  !
  !       tures are subsequently deallocated on thread 0 -- the checker will   !
  !       think there are more deallocs than allocs. This happens for hash     !
  !       tables in HFx.                                                       !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2005.                                       !
  ! Extended with 'allocated_in' by Jacek Dziedzic, November 2019.             !
  !============================================================================!

  subroutine utils_dealloc_check(routine,array,status,allocated_in)

    use comms, only: pub_on_root
    use constants, only: stdout
    implicit none

    ! Arguments
    character(len=*), intent(in) :: routine    ! Routine name
    character(len=*), intent(in) :: array      ! Array name
    integer, intent(in)          :: status     ! Status flag
    character(len=*), intent(in), optional :: allocated_in ! Routine name

    ! Locals
    integer :: iarray
    logical :: any_routine
    character(72) :: routine_no_de
    character(72) :: loc_allocated_in

!$  integer, external :: omp_get_thread_num

    if (status /= 0) then
       call utils_abort('Error in '//trim(routine)//': deallocating '//&
            trim(array)//' failed with code ',status)
    end if

    if (.not. (loc_debug_arrays .and. debug_arrays_alloc) ) return

    if(present(allocated_in)) then
       loc_allocated_in = allocated_in
    else
       loc_allocated_in = ''
    end if

!$  if (omp_get_thread_num()==0) then
    routine_no_de = ''
    iarray = index(routine,'dealloc')
    if (iarray>0) then
       routine_no_de(1:iarray) = routine(1:iarray)
       routine_no_de(iarray:(len_trim(routine)-2)) = &
            routine((iarray+2):len_trim(routine))
    else
       routine_no_de = routine
    end if
    any_routine = .false.
    iarray = index(routine,'_exit')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'_free')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'_destroy')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'_remove')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'_close')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'_[close]')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'_purge')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'classical_pot_dealloc')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'function_basis_deallocate')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'basis_sphere_deallocate')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'sparse_segments_dealloc')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'kernel_workspace_dealloc')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'hash_table_add_data') ! jd: as it potentially reallocs
    if (iarray>0) any_routine = .true.

    do iarray=1,num_arrays
       ! ndmh: find array in list of allocated arrays
       if ((array == arrays(iarray)%name) .and. &
            ((routine_no_de == arrays(iarray)%routine).or. &
            any_routine .or. loc_allocated_in == arrays(iarray)%routine)) then
          ! ndmh: decrement number of allocated versions
          if (arrays(iarray)%allocated > 0) then
             arrays(iarray)%allocated = arrays(iarray)%allocated - 1
             return
          end if
       end if
    end do

    if (pub_on_root) write(stdout,'(3a/3a/l1)') 'WARNING in &
         &utils_dealloc_check: matching tag "',trim(array),'" not found', &
         'when deallocating in routine "',trim(routine),'"', any_routine

!$  end if

  end subroutine utils_dealloc_check

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !============================================================================!
  ! Flushes output (stdout by default)                                         !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes 12 November 2004                                   !
  ! Modified by Jacek Dziedzic on 08/06/2010 to accept an optional parameter,  !
  ! which, if .true., abandons the MPI reduction at the end. This allows one to!
  ! call this routine on certain procs only, from places where services_flush  !
  ! is not accessible due to circular dependencies.                            !
  ! Default behaviour remains unchanged.                                       !
  ! Simplified to use Fortran2003 intrinsic by JM Escartin 5 June 2015         !
  !============================================================================!

  subroutine utils_flush(unit, dont_reduce_error_flag)

    use comms, only : pub_on_root, comms_reduce
    use constants, only: stdout

    implicit none

    ! Arguments
    integer, optional, intent(in) :: unit
    logical, optional, intent(in) :: dont_reduce_error_flag

    ! Local variables
    integer :: unit_local     ! Local copy of argument
    integer :: ierr           ! Error flag

    if (present(unit)) then
       unit_local = unit
    else
       unit_local = stdout
    end if

    flush(unit_local, iostat=ierr)

    ! jd: Prevent reduction, if asked to, so that flushing on selected
    !     procs becomes possible
    if (present(dont_reduce_error_flag)) then
       if (dont_reduce_error_flag) return
    end if

    ! jme: Make sure that iostat values don't cancel (we only want to know
    ! if any of them was nonzero)
    ierr = abs(ierr)

    ! Compare error flags
    call comms_reduce('SUM', ierr)
    if (ierr /= 0 .and. pub_on_root) write(stdout,'(a)') &
         'WARNING in utils_flush: flush failed on at least one processor'

  end subroutine utils_flush

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !============================================================================!
  ! Finds a free unit for I/O                                                  !
  !----------------------------------------------------------------------------!
  ! Amalgamation of restart_find_unit by Chris-Kriton Skylaris and esdf_unit   !
  ! by Chris Pickard.                                                          !
  !============================================================================!

  integer function utils_unit()

    ! Local variables
    integer :: trial_unit
    integer :: ierr
    logical :: ex
    logical :: op

    utils_unit = -1
    do trial_unit=10,99
       inquire(unit=trial_unit,exist=ex,opened=op,iostat=ierr)
       call utils_assert(ierr == 0, 'Error in utils_unit(): inquiring about uni&
            &t failed. Unit number and error were as follows: ',trial_unit,ierr)
       if (ex .and. (.not. op)) then
          utils_unit = trial_unit
          exit
       end if
    end do

    call utils_assert(utils_unit /= -1, 'Error in utils_unit(): &
         &no I/O units available.')

  end function utils_unit

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !============================================================================!
  ! Copies data from/to binary files in the extensible format used for the     !
  ! density kernel.                                                            !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, February 2007.                                    !
  !============================================================================!

  subroutine utils_binary_copy(iunit,ounit)

    use constants, only: DP

    implicit none

    ! Arguments
    integer, intent(in) :: iunit              ! Input I/O unit
    integer, intent(in) :: ounit              ! Output I/O unit

    ! Local variables
    integer :: ierr                           ! Error flag
    integer :: ilen                           ! Record length
    character(len=12) :: tag                  ! Record tag
    integer, allocatable :: ibuf(:)           ! I/O buffer for integer data
    real(kind=DP), allocatable :: dbuf(:)     ! I/O buffer for real data
    complex(kind=DP), allocatable :: zbuf(:)  ! I/O buffer for complex data


    ! Indefinite loop until endfile record found
    do
       ! Read tag for record and check for endfile
       read(iunit) tag
       if (tag(1:7) == 'ENDFILE') then
          write(ounit) tag
          write(ounit) 1
          write(ounit) 0
          exit
       end if


       ! Check record type is valid
       if (tag(12:12) /= 'I' .and. tag(12:12) /= 'D' .and. &
            tag(12:12) /= 'Z') then
          call utils_abort('Error in utils_binary_copy(): malformed record tag.')
       end if
       write(ounit) tag


       ! Read record length
       read(iunit) ilen
       if (ilen <= 0) then
          call utils_abort('Error in utils_binary_copy(): malformed record length.')
       end if
       write(ounit) ilen


       ! Allocate workspace for this type as required
       if (tag(12:12) == 'I') then
          if (allocated(ibuf)) then
             if (size(ibuf) < ilen) then
                deallocate(ibuf,stat=ierr)
                call utils_dealloc_check('utils_binary_copy','ibuf',ierr)
             end if
          end if
          if (.not. allocated(ibuf)) then
             allocate(ibuf(ilen),stat=ierr)
             call utils_alloc_check('utils_binary_copy','ibuf',ierr)
          end if
          read(iunit) ibuf(1:ilen)
          write(ounit) ibuf(1:ilen)
          deallocate(ibuf,stat=ierr)
          call utils_dealloc_check('utils_binary_copy','ibuf',ierr)
       else if (tag(12:12) == 'D') then
          if (allocated(dbuf)) then
             if (size(dbuf) < ilen) then
                deallocate(dbuf,stat=ierr)
                call utils_dealloc_check('utils_binary_copy','dbuf',ierr)
             end if
          end if
          if (.not. allocated(dbuf)) then
             allocate(dbuf(ilen),stat=ierr)
             call utils_alloc_check('utils_binary_copy','dbuf',ierr)
          end if
          read(iunit) dbuf(1:ilen)
          write(ounit) dbuf(1:ilen)
          deallocate(dbuf,stat=ierr)
          call utils_dealloc_check('utils_binary_copy','dbuf',ierr)
       else if (tag(12:12) == 'Z') then
          if (allocated(zbuf)) then
             if (size(zbuf) < ilen) then
                deallocate(zbuf,stat=ierr)
                call utils_dealloc_check('utils_binary_copy','zbuf',ierr)
             end if
          end if
          if (.not. allocated(zbuf)) then
             allocate(zbuf(ilen),stat=ierr)
             call utils_alloc_check('utils_binary_copy','zbuf',ierr)
          end if
          read(iunit) zbuf(1:ilen)
          write(ounit) zbuf(1:ilen)
          deallocate(zbuf,stat=ierr)
          call utils_dealloc_check('utils_binary_copy','zbuf',ierr)

       end if

       ! End loop over file
    end do

    return

  end subroutine utils_binary_copy

  !--------------------------------------------------------------------------

  subroutine utils_heapsort_integer(n,key,idx)

    !=======================================================================!
    ! This subroutine creates an index based on a list of keys, which are   !
    ! sorted by the heap sort algorithm, which is guaranteed to complete    !
    ! after O(n log n) operations.                                          !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! n (input)          : Number of items in list to be indexed            !
    ! key (input/output) : Key according to which the items are sorted      !
    ! idx (output)       : The resulting index                              !
    !-----------------------------------------------------------------------!
    ! Written by Peter Haynes, 10/7/03                                      !
    !=======================================================================!

    use constants, only: LONG

    implicit none

    ! Arguments
    integer, intent(in) :: n
    integer(kind=LONG), intent(inout) :: key(n)
    integer, intent(out) :: idx(n)

    ! Local variables
    integer :: i,j                 ! Loop variables over list
    integer :: ic,is               ! Counters for heap creation/selection
    integer :: temp_idx            ! Temporary copy of index
    integer(kind=LONG) :: temp_key ! Temporary copy of key

    ! If the list is empty, there's nothing to do!
    if (n == 0) return

    ! If there is only one item in the list, there's no sorting to be done!
    if (n == 1) then
       idx(1) = 1
       return
    end if

    ! Initialise the sorting index
    do i=1,n
       idx(i) = i
    end do

    ! Do the heap sort
    ic = n/2+1
    is = n
    do
       if (ic > 1) then ! in heap creation phase
          ic = ic - 1
          temp_key = key(ic) ; temp_idx = idx(ic)
       else ! in heap selection phase
          temp_key = key(is) ; temp_idx = idx(is)
          key(is) = key(1) ; idx(is) = idx(1)
          is = is - 1
          if (is == 1) then
             key(1) = temp_key ; idx(1) = temp_idx
             exit
          end if
       end if
       ! Sift down temporary copy to correct level in heap
       i = ic
       j = 2 * ic
       do while (j <= is)
          if (j < is) then
             if (key(j) < key(j+1)) j = j + 1
          end if
          if (temp_key < key(j)) then ! demote temporary copy
             key(i) = key(j) ; idx(i) = idx(j)
             i = j
             j = 2 * j
          else ! found level for temporary copy
             j = is + 1
          end if
       end do
       key(i) = temp_key ; idx(i) = temp_idx
    end do

  end subroutine utils_heapsort_integer

  !--------------------------------------------------------------------------

  subroutine utils_heapsort_real(n,key,idx)

    !=======================================================================!
    ! This subroutine creates an index based on a list of keys, which are   !
    ! sorted by the heap sort algorithm, which is guaranteed to complete    !
    ! after O(n log n) operations.                                          !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! n (input)          : Number of items in list to be indexed            !
    ! key (input/output) : Key according to which the items are sorted      !
    ! idx (output)       : The resulting index                              !
    !-----------------------------------------------------------------------!
    ! Written by Peter Haynes, 10/7/03                                      !
    !=======================================================================!

    use constants, only: DP

    implicit none

    ! Arguments
    integer, intent(in) :: n
    real(kind=DP), intent(inout) :: key(n)
    integer, intent(out) :: idx(n)

    ! Local variables
    integer :: i,j                 ! Loop variables over list
    integer :: ic,is               ! Counters for heap creation/selection
    integer :: temp_idx            ! Temporary copy of index
    real(kind=DP) :: temp_key      ! Temporary copy of key

    ! If the list is empty, there's nothing to do!
    if (n == 0) return

    ! If there is only one item in the list, there's no sorting to be done!
    if (n == 1) then
       idx(1) = 1
       return
    end if

    ! Initialise the sorting index
    do i=1,n
       idx(i) = i
    end do

    ! Do the heap sort
    ic = n/2+1
    is = n
    do
       if (ic > 1) then ! in heap creation phase
          ic = ic - 1
          temp_key = key(ic) ; temp_idx = idx(ic)
       else ! in heap selection phase
          temp_key = key(is) ; temp_idx = idx(is)
          key(is) = key(1) ; idx(is) = idx(1)
          is = is - 1
          if (is == 1) then
             key(1) = temp_key ; idx(1) = temp_idx
             exit
          end if
       end if
       ! Sift down temporary copy to correct level in heap
       i = ic
       j = 2 * ic
       do while (j <= is)
          if (j < is) then
             if (key(j) < key(j+1)) j = j + 1
          end if
          if (temp_key < key(j)) then ! demote temporary copy
             key(i) = key(j) ; idx(i) = idx(j)
             i = j
             j = 2 * j
          else ! found level for temporary copy
             j = is + 1
          end if
       end do
       key(i) = temp_key ; idx(i) = temp_idx
    end do

  end subroutine utils_heapsort_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_open_unit_check(routine,file_name,ios_status)

    !==========================================================================!
    ! This subroutine checks whether an unit has been opened successfully.     !
    ! On failure it writes an error message and aborts.                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   routine   (in) : Name of the subroutine which called this routine      !
    !   file_name (in) : Name of the file to be opened                         !
    !   ios_status(in) : Status flag returned by open                          !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano, July 2009.                               !
    !==========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: routine    ! Routine name
    character(len=*), intent(in) :: file_name  ! File name
    integer, intent(in)          :: ios_status ! Status flag

    call utils_assert(ios_status == 0, 'Error in '//trim(routine)//': failed &
         &to open file "'//trim(file_name)//'" with code ',ios_status)

  end subroutine utils_open_unit_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_close_unit_check(routine,unit_name,ios_status)

    !==========================================================================!
    ! This subroutine checks whether an unit has been closed successfully.     !
    ! On failure it writes an error message and aborts.                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   routine   (in) : Name of the subroutine which called this routine      !
    !   unit_name (in) : Name of the unit to be closed                         !
    !   ios_status(in) : Status flag returned by open                          !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano, July 2009.                               !
    !==========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: routine    ! Routine name
    character(len=*), intent(in) :: unit_name  ! File name
    integer, intent(in)          :: ios_status ! Status flag

    call utils_assert(ios_status == 0, 'Error in '//trim(routine)//': failed &
         &to close unit "'//trim(unit_name)//'" with code ',ios_status)

  end subroutine utils_close_unit_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_read_check(routine,var_name,ios_status)

    !==========================================================================!
    ! This subroutine checks whether a value has been readed successfully.     !
    ! On failure it writes an error message and aborts.                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   routine   (in) : Name of the subroutine which called this routine      !
    !   var_name  (in) : Name of the variable to be read                       !
    !   ios_status(in) : Status flag returned by open                          !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano, July 2009.                               !
    !==========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: routine    ! Routine name
    character(len=*), intent(in) :: var_name   ! Variable name
    integer, intent(in)          :: ios_status ! Status flag

    call utils_assert(ios_status == 0, 'Error in '//trim(routine)//': failed &
         &to read variable "'//trim(var_name)//'" with code ',ios_status)

  end subroutine utils_read_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_write_check(routine,var_name,ios_status)

    !==========================================================================!
    ! This subroutine checks whether a value has been written successfully.    !
    ! On failure it writes an error message and aborts.                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   routine   (in) : Name of the subroutine which called this routine      !
    !   var_name  (in) : Name of the variable to be read                       !
    !   ios_status(in) : Status flag returned by open                          !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano, July 2009.                               !
    !==========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: routine    ! Routine name
    character(len=*), intent(in) :: var_name   ! Variable name
    integer, intent(in)          :: ios_status ! Status flag

    call utils_assert(ios_status == 0, 'Error in '//trim(routine)//': failed &
         &to write variable "'//trim(var_name)//'" with code ',ios_status)

  end subroutine utils_write_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function utils_custom_erfc(x)
    !=========================================================================!
    ! Calculate an accurage approximation to the complemenary error function, !
    !    erfc(x)=(2/sqrt(pi))*integral(x->infinity)[exp(-t^2)]dt.             !
    !    Based upon parameterization given in NSWC Mathematics Library.       !
    !    NB This is machine portable and not dependent on a system function.  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   x, intent=in, argument of erfc to evaluate.                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v0.1, 25/09/2001                               !
    !=========================================================================!

    use constants, only: DP

    implicit none

    real(kind=dp), intent(in) :: x
    real(kind=dp)             :: utils_custom_erfc

    !expansion parameters
    real(kind=dp), dimension(1:5) :: a=(/  0.771058495001320E-04_dp, -0.133733772997339E-02_dp,  &
         &  0.323076579225834E-01_dp,  0.479137145607681E-01_dp,  0.128379167095513E+00_dp /)
    real(kind=dp), dimension(1:3) :: b=(/  0.301048631703895E-02_dp,  0.538971687740286E-01_dp,  &
         &  0.375795757275549E+00_dp /)
    real(kind=dp), dimension(1:8) :: p=(/ -1.36864857382717E-07_dp,  5.64195517478974E-01_dp, &
         &  7.21175825088309E+00_dp,  4.31622272220567E+01_dp,  1.52989285046940E+02_dp,  &
         &  3.39320816734344E+02_dp,  4.51918953711873E+02_dp,  3.00459261020162E+02_dp /)
    real(kind=dp), dimension(1:8) :: q=(/  1.00000000000000E+00_dp,  1.27827273196294E+01_dp, &
         &  7.70001529352295E+01_dp, 2.77585444743988E+02_dp,  6.38980264465631E+02_dp, &
         &  9.31354094850610E+02_dp, 7.90950925327898E+02_dp,  3.00459260956983E+02_dp /)
    real(kind=dp), dimension(1:5) :: r=(/  2.10144126479064E+00_dp,  2.62370141675169E+01_dp, &
         &  2.13688200555087E+01_dp, 4.65807828718470E+00_dp,  2.82094791773523E-01_dp /)
    real(kind=dp), dimension(1:4) :: s=(/ 9.41537750555460E+01_dp, 1.87114811799590E+02_dp, &
         &  9.90191814623914E+01_dp, 1.80124575948747E+01_dp /)
    real(kind=dp)                 :: c = .564189583547756_dp

    ! Local variables
    real(kind=dp) :: ax, x2, t, top, bot

    ax=abs(x)
    x2=x*x
    if (ax<0.5_dp) then
       t=x2
       top=((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0_dp
       bot=((  b(1)*t + b(2))*t + b(3)) * t + 1.0_dp
       utils_custom_erfc=0.5_dp + (0.5_dp-x*(top/bot))

    else if (ax<4.0_dp) then
       top=((((((p(1) *ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax + p(6))*ax  &
            &  + p(7))*ax + p(8)
       bot=((((((q(1) *ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax + q(6))*ax  &
            &  + q(7))*ax + q(8)
       utils_custom_erfc=exp(-x2)*top/bot
       if (x<0.0_dp) utils_custom_erfc=2.0_dp - utils_custom_erfc

    else
       if (x<=-5.6_dp) then
          utils_custom_erfc=2.0_dp    !large negative x limit

       else if (x>100.0_dp) then
          utils_custom_erfc=0.0_dp    !large positive x limit

       else
          t=1.0_dp/x2
          top=(((r(1)*t + r(2))*t + r(3))*t + r(4)) * t + r(5)
          bot=(((s(1)*t + s(2))*t + s(3))*t + s(4)) * t + 1.0_dp
          utils_custom_erfc=exp(-x2)*(c-t*top/bot)/ax
          if (x<0.0_dp) utils_custom_erfc=2.0_dp - utils_custom_erfc

       end if
    end if

    return
  end function utils_custom_erfc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function utils_erf(x)
    !==========================================================================!
    ! Returns erf(x), calculated from derf(x) which is either vendor-provided, !
    ! or linked from libc. Alternatively, iff CUSTOM_ERF is defined, it uses   !
    ! utils_custom_erfc, defined above, to compute erf(x).                     !
    !--------------------------------------------------------------------------!
    ! Both versions are accurate to at least 3D-15, utils_custom_erfc() is     !
    ! several times slower.                                                    !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, November 2010.                                !
    !==========================================================================!

    use constants, only: DP

    implicit none

    real(kind=DP) :: utils_erf

    ! jd: PGI and FLANG require this to be declared
#if defined(__PGI) || defined __FLANG
    real(kind=DP) :: derf
#endif

    ! Arguments
    real(kind=DP), intent(in) :: x

    ! -----------------------------------------------------------------------

#ifndef CUSTOM_ERF
#if defined(__PGI) || defined __FLANG
    utils_erf = derf(x)
#else
    utils_erf = erf(x)
#endif
#else
    utils_erf = 1.0_DP - utils_custom_erfc(x)
#endif

  end function utils_erf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function utils_erfc(x)
    !==========================================================================!
    ! Returns erfc(x), calculated from derfc(x), which is either               !
    ! vendor-provided, or linked from libc. Alternatively, iff CUSTOM_ERF is   !
    ! defined, it uses utils_custom_erfc, defined above.                       !
    !--------------------------------------------------------------------------!
    ! Both versions are accurate to at least 3D-15, utils_custom_erfc() is     !
    ! several times slower.                                                    !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, November 2010.                                !
    !==========================================================================!

    use constants, only: DP

    implicit none

    real(kind=DP) :: utils_erfc

    ! jd: PGI and FLANG require this to be declared
#if defined(__PGI) || defined __FLANG
    real(kind=DP) :: derfc
#endif

    ! Arguments
    real(kind=DP), intent(in) :: x

    ! -----------------------------------------------------------------------

#ifndef CUSTOM_ERF
#if defined(__PGI) || defined __FLANG
    utils_erfc = derfc(x)
#else
    utils_erfc = erfc(x)
#endif
#else
    utils_erfc = utils_custom_erfc(x)
#endif

  end function utils_erfc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine utils_abort(error_message, opt_int_to_print1, opt_int_to_print2, &
       opt_int_to_print3, opt_int_to_print4, opt_int_to_print5, &
       opt_real_to_print1, opt_real_to_print2, &
       opt_logical_to_print1, opt_logical_to_print2, opt_logical_to_print3, &
       file_output, avoid_mpi_calls)
    !=========================================================================!
    ! Aborts ONETEP, making a best effort to produce a single copy of the     !
    ! error message on standard output and in a .error_message file.          !
    !
    ! Can be called on any combination of procs, including, or not, the root  !
    ! proc. If after about 10 seconds the procs cannot agree on who should be !
    ! the spokesman that should produce the error message (e.g. due to heavy  !
    ! comms delays or huge numbers of procs), everyone has a go at producing  !
    ! the error message ('graceless abort').                                  !
    !                                                                         !
    ! Aborting (as in calling comms_abort eventually) is guaranteed in all    !
    ! cases.                                                                  !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   error_message (in): The error message to be displayed before ONETEP   !
    !                       is terminated.                                    !
    !   opt_int_to_print1..5 (in, optional): Optional int values to be        !
    !                                        printed after the error message. !
    !   opt_real_to_print1..2 (in, optional): Optional real value to be       !
    !                                         printed after the error message.!
    !   file_output (in, optional): if present, and .false., do not output    !
    !                               to file (necessary if root filename is    !
    !                               too long or filesystem is known to be     !
    !                               inaccessible, default is .true.)          !
    !   avoid_mpi_calls (in, optional): will not elect a spokesman to avoid   !
    !                                   MPI calls. This is useful when called !
    !                                   from OMP regions. The abort will be   !
    !                                   messy (multiple outputs), but al least!
    !                                   it won't crash by trying to do MPI    !
    !                                   from multiple OMP contexts.           !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 07/2012                                      !
    ! Updated by Edward Linscott 05/2019, moving the code for updating        !
    ! pub_rootname_local to a separate routine (utils_update_rootname)        !
    !=========================================================================!

    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_send, &
         comms_test, comms_irecv, comms_abort_for_use_in_utils_abort

    implicit none

    ! Arguments
    character(len=*), intent(in)            :: error_message
    integer, intent(in), optional           :: opt_int_to_print1
    integer, intent(in), optional           :: opt_int_to_print2
    integer, intent(in), optional           :: opt_int_to_print3
    integer, intent(in), optional           :: opt_int_to_print4
    integer, intent(in), optional           :: opt_int_to_print5
    real(kind=DP), intent(in), optional     :: opt_real_to_print1
    real(kind=DP), intent(in), optional     :: opt_real_to_print2
    logical, intent(in), optional           :: opt_logical_to_print1
    logical, intent(in), optional           :: opt_logical_to_print2
    logical, intent(in), optional           :: opt_logical_to_print3
    logical, intent(in), optional           :: file_output
    logical, intent(in), optional           :: avoid_mpi_calls

    ! jd: Internal variables
    integer :: val_to_send
    integer :: send_handles(0:pub_total_num_procs-1)
    integer :: recv_handles(0:pub_total_num_procs-1)
    integer :: src, dest
    integer :: recv(0:pub_total_num_procs-1)
    logical :: test(0:pub_total_num_procs-1)
    integer :: sleep_counter
    integer :: max_sleep_milliseconds
    character(len=1024) :: error_message_header
    character(len=1024) :: error_message_file
    integer :: error_unit
    integer :: spokesman
!$  integer, external :: omp_get_thread_num
    integer           :: tid
    character(len=256) :: int_1_as_string, int_2_as_string, int_3_as_string, &
         int_4_as_string, int_5_as_string, real_1_as_string, real_2_as_string, &
         logical_1_as_string, logical_2_as_string, logical_3_as_string
    logical :: loc_file_output
    logical :: loc_avoid_mpi_calls

    integer, parameter :: sleep_tick = 50

    !------------------------------------------------------------------------

    ! JCW: loc_file_output controls whether file operations are done.
    ! JCW: By default, they are, but this may be disabled by passing
    ! JCW: file_output = .false. as an optional argument
    loc_file_output = .true.
    if (present(file_output)) loc_file_output = file_output

    loc_avoid_mpi_calls = .false.
    if (present(avoid_mpi_calls)) loc_avoid_mpi_calls = avoid_mpi_calls

    tid = 0
!$  tid = omp_get_thread_num()

    error_message_header=''

    ! Convert optional values to be printed to strings
    int_1_as_string=''
    int_2_as_string=''
    int_3_as_string=''
    int_4_as_string=''
    int_5_as_string=''
    real_1_as_string=''
    real_2_as_string=''
    logical_1_as_string=''
    logical_2_as_string=''
    logical_3_as_string=''
    if(present(opt_int_to_print1)) then
       write(int_1_as_string,'(i32)') opt_int_to_print1
       int_1_as_string = trim(adjustl(int_1_as_string))
    end if
    if(present(opt_int_to_print2)) then
       write(int_2_as_string,'(i32)') opt_int_to_print2
       int_2_as_string = trim(adjustl(int_2_as_string))
    end if
    if(present(opt_int_to_print3)) then
       write(int_3_as_string,'(i32)') opt_int_to_print3
       int_3_as_string = trim(adjustl(int_3_as_string))
    end if
    if(present(opt_int_to_print4)) then
       write(int_4_as_string,'(i32)') opt_int_to_print4
       int_4_as_string = trim(adjustl(int_4_as_string))
    end if
    if(present(opt_int_to_print5)) then
       write(int_5_as_string,'(i32)') opt_int_to_print5
       int_5_as_string = trim(adjustl(int_5_as_string))
    end if
    if(present(opt_real_to_print1)) &
         write(real_1_as_string,'(f16.7)') opt_real_to_print1
    if(present(opt_real_to_print2)) &
         write(real_2_as_string,'(f16.7)') opt_real_to_print2
    if(present(opt_logical_to_print1)) &
         write(logical_1_as_string,'(l1)') opt_logical_to_print1
    if(present(opt_logical_to_print2)) &
         write(logical_2_as_string,'(l1)') opt_logical_to_print2
    if(present(opt_logical_to_print3)) &
         write(logical_3_as_string,'(l1)') opt_logical_to_print3

    ! Try to elect a spokesman, unless there is only one proc (which makes the
    ! selection trivial) and unless we're rank 0 (which automatically makes us
    ! the spokesman)
    if(pub_total_num_procs > 1 .and. pub_my_proc_id > 0 .and. &
         .not. loc_avoid_mpi_calls .and. tid == 0) then
       val_to_send = pub_my_proc_id
       test(:) = .false.

       ! Attempt to send a trivial message to everyone, non-blocking
       do dest = 0, pub_total_num_procs-1
          call comms_send(dest, pub_my_proc_id, 1, tag = pub_my_proc_id, &
               return_handle = send_handles(dest), add_to_stack = .false.)
       end do

       ! Post the receives for same
       do src = 0, pub_total_num_procs-1
          call comms_irecv(src, recv(src), 1, tag = src, &
               handle = recv_handles(src))
       end do

       ! Wait until everyone is pinged successfully, or timeout after at most
       ! a few seconds
       max_sleep_milliseconds = int(1000.0_DP * &
            log10(real(pub_total_num_procs,kind=DP)))
       do sleep_counter = 1, max_sleep_milliseconds/sleep_tick
          do src = 0, pub_total_num_procs-1
             call comms_test(test(src),recv_handles(src))
          end do
          if(all(test)) exit
          call utils_busy_wait(sleep_tick)
       end do

       ! Find the proc with the smallest index that was succesfully pinged,
       ! this will be the spokesman that outputs the error message
       do spokesman = 0, pub_total_num_procs-1
          if(test(spokesman)) exit
       end do
    else
       spokesman = 0
    end if

    ! Have the spokesman output the error message
    if(pub_my_proc_id == spokesman .or. loc_avoid_mpi_calls) then
       call internal_output_error_message

       ! Wait for a short while in the hope that the error message will get to
       ! the screen in its entirety before the abort message (interactive jobs)
       call utils_busy_wait(200)

       ! Go down in flames
       call comms_abort_for_use_in_utils_abort

       ! In the unlikely case MPI_Abort returns
       stop!!cannot use utils_abort() here, obviously

    end if

    ! This is what non-spokesmen do -- patiently wait to be aborted from the
    ! outside, through the comms_abort() executed by the spokesman.

    ! This also handles the corner case where a spokesman could not be elected,
    ! e.g. thousands of procs called utils_abort and could not communicate
    ! within the timeout window. Each of them might have a different idea of who
    ! the spokesman is, and, if we're unlucky, everyone will think it's someone
    ! else. If so, abort gracelessly (everyone has a go at outputting the error
    ! message and aborting).

    ! Wait for 10 seconds
    do sleep_counter = 1, 10
       call utils_busy_wait(1000)

       ! Dummy MPI call, might help MPI realize an abort is in progress
       if(tid /= 0) call comms_test(test(0),recv_handles(0))
    end do

    ! Try to output the message from any ranks that got here
    error_message_header = &
         'NB: There were problems performing a graceful abort. &
         &ONETEP will now try to abort with no regard for grace.'

    call internal_output_error_message

    ! And go down in flames
    call comms_abort_for_use_in_utils_abort

    ! In the unlikely case MPI_Abort returns
    stop!!cannot use utils_abort() here, obviously

contains

    subroutine internal_output_error_message

       use constants, only: stdout

       implicit none

       ! ----------------------------------------------------------------------

       if (loc_file_output) then
          ! JCW: Only give error_message_file a name if file output is
          ! JCW: allowed
          write(error_message_file,'(a,a)') trim(pub_rootname_local), &
               '.error_message'
       end if

       write(stdout,'(a)') ''
       write(stdout,'(a)') repeat('+',80)
       write(stdout,'(a)') trim(error_message_header)
       write(stdout,'(a/)') &
            'ONETEP terminated abnormally due to the following error:'
       write(stdout,'(a)') trim(error_message)//' '//&
            trim(int_1_as_string)//' '//&
            trim(int_2_as_string)//' '//&
            trim(int_3_as_string)//' '//&
            trim(int_4_as_string)//' '//&
            trim(int_5_as_string)//' '//&
            trim(real_1_as_string)//' '//&
            trim(real_2_as_string)//' '//&
            trim(logical_1_as_string)//' '//&
            trim(logical_2_as_string)//' '//&
            trim(logical_3_as_string)
       write(stdout,'(a,i0,a)') '... detected on MPI rank #', pub_my_proc_id,'.'
       write(stdout,'(a)') repeat('+',80)
       write(stdout,'(a)') ''
       call utils_flush(stdout,.true.)

       ! JCW: Only output to .error_message file if allowed by loc_file_output
       if (loc_file_output) then
          error_unit = utils_unit()
          open(error_unit, file=error_message_file, status="new", err=1000) ! (*)
          write(error_unit,'(a)') repeat('+',80)
          write(error_unit,'(a)') trim(error_message_header)
          write(error_unit,'(a/)') &
               'ONETEP terminated abnormally due to the following error:'
          write(error_unit,'(a)') trim(error_message)//' '//&
               trim(int_1_as_string)//' '//&
               trim(int_2_as_string)//' '//&
               trim(int_3_as_string)//' '//&
               trim(int_4_as_string)//' '//&
               trim(int_5_as_string)//' '//&
               trim(real_1_as_string)//' '//&
               trim(real_2_as_string)//' '//&
               trim(logical_1_as_string)//' '//&
               trim(logical_2_as_string)//' '//&
               trim(logical_3_as_string)
          write(error_unit,'(a,i0,a)') '... detected on MPI rank #', pub_my_proc_id,'.'
          write(error_unit,'(a)') repeat('+',80)

          call utils_flush(error_unit,.true.)
          close(error_unit)
       end if

       return

       ! (*) When aborting gracelessly, with many procs racing to output the
       !     error message, only allow one of them to open the file to have
       !     a semblance of order. The late procs will remain silent.
1000   continue

    end subroutine internal_output_error_message

  end subroutine utils_abort

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_assert_with_integers(condition,error_message, &
       value_1_to_print,value_2_to_print,value_3_to_print,value_4_to_print)

    !==========================================================================!
    ! Asserts that 'condition' is .true., if not -- forwards 'error_message'   !
    ! followed a set of 'values_to_print' to utils_assert_with_string.         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   condition (input): The condition to check.                             !
    !   error_message (input): The error message to print if assertion fails.  !
    !   value_[1..4]_to_print (input): The values to print if assertion fails. !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Avoid using conditions with side effects (such as function calls).     !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2010.                                   !
    !==========================================================================!

    ! Arguments
    logical, intent(in) :: condition
    character(len=*), intent(in) :: error_message
    integer, intent(in) :: value_1_to_print
    integer, intent(in), optional :: value_2_to_print
    integer, intent(in), optional :: value_3_to_print
    integer, intent(in), optional :: value_4_to_print

    ! Local variables
    character(len=256) :: value_1_as_string, value_2_as_string, &
         value_3_as_string, value_4_as_string

    ! -----------------------------------------------------------------------

    if(condition) then
       return ! jd: Incur no extra overhead
    else
       write(value_1_as_string,'(i32)') value_1_to_print
       value_1_as_string = trim(adjustl(value_1_as_string))
       value_2_as_string=''
       value_3_as_string=''
       value_4_as_string=''
       if(present(value_2_to_print)) then
          write(value_2_as_string,'(i32)') value_2_to_print
          value_2_as_string = trim(adjustl(value_2_as_string))
       end if
       if(present(value_3_to_print)) then
          write(value_3_as_string,'(i32)') value_3_to_print
          value_3_as_string = trim(adjustl(value_3_as_string))
       end if
       if(present(value_4_to_print)) then
          write(value_4_as_string,'(i32)') value_4_to_print
          value_4_as_string = trim(adjustl(value_4_as_string))
       end if
       call utils_assert_with_string(condition,error_message, &
            trim(value_1_as_string)//' '//&
            trim(value_2_as_string)//' '//&
            trim(value_3_as_string)//' '//&
            trim(value_4_as_string))
    end if

  end subroutine utils_assert_with_integers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_assert_with_reals(condition, error_message, &
       value_1_to_print, value_2_to_print, value_3_to_print, value_4_to_print)

    !==========================================================================!
    ! Asserts that 'condition' is .true., if not -- forwards 'error_message'   !
    ! followed a set of 'values_to_print' to utils_assert_with_string.         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   condition (input): The condition to check.                             !
    !   error_message (input): The error message to print if assertion fails.  !
    !   value_[1..4]_to_print (input): The values to print if assertion fails. !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Avoid using conditions with side effects (such as function calls).     !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2010.                                   !
    !==========================================================================!

    use constants, only: DP

    ! Arguments
    logical, intent(in)                 :: condition
    character(len=*), intent(in)        :: error_message
    real(kind=DP), intent(in)           :: value_1_to_print
    real(kind=DP), intent(in), optional :: value_2_to_print
    real(kind=DP), intent(in), optional :: value_3_to_print
    real(kind=DP), intent(in), optional :: value_4_to_print

    ! Local variables
    character(len=256) :: value_1_as_string, value_2_as_string, &
         value_3_as_string, value_4_as_string

    ! -----------------------------------------------------------------------

    if(condition) then
       return ! jd: Incur no extra overhead
    else
       write(value_1_as_string,'(e20.12)') value_1_to_print
       value_2_as_string=''
       value_3_as_string=''
       value_4_as_string=''
       if(present(value_2_to_print)) &
            write(value_2_as_string,'(e20.12)') value_2_to_print
       if(present(value_3_to_print)) &
            write(value_3_as_string,'(e20.12)') value_3_to_print
       if(present(value_4_to_print)) &
            write(value_4_as_string,'(e20.12)') value_4_to_print
       call utils_assert_with_string(condition,error_message, &
            trim(value_1_as_string)//' '//&
            trim(value_2_as_string)//' '//&
            trim(value_3_as_string)//' '//&
            trim(value_4_as_string))
    end if

  end subroutine utils_assert_with_reals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_assert_with_logicals(condition, error_message, &
       value_1_to_print, value_2_to_print, value_3_to_print, value_4_to_print)

    !==========================================================================!
    ! Asserts that 'condition' is .true., if not -- forwards 'error_message'   !
    ! followed a set of 'values_to_print' to utils_assert_with_string.         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   condition (input): The condition to check.                             !
    !   error_message (input): The error message to print if assertion fails.  !
    !   value_[1..4]_to_print (input): The values to print if assertion fails. !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Avoid using conditions with side effects (such as function calls).     !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2010.                                   !
    !==========================================================================!

    ! Arguments
    logical, intent(in)                 :: condition
    character(len=*), intent(in)        :: error_message
    logical, intent(in)                 :: value_1_to_print
    logical, intent(in), optional       :: value_2_to_print
    logical, intent(in), optional       :: value_3_to_print
    logical, intent(in), optional       :: value_4_to_print

    ! Local variables
    character(len=256) :: value_1_as_string, value_2_as_string, &
         value_3_as_string, value_4_as_string

    ! -----------------------------------------------------------------------

    if(condition) then
       return ! jd: Incur no extra overhead
    else
       write(value_1_as_string,'(l1)') value_1_to_print
       value_2_as_string=''
       value_3_as_string=''
       value_4_as_string=''
       if(present(value_2_to_print)) &
            write(value_2_as_string,'(l1)') value_2_to_print
       if(present(value_3_to_print)) &
            write(value_3_as_string,'(l1)') value_3_to_print
       if(present(value_4_to_print)) &
            write(value_4_as_string,'(l1)') value_4_to_print
       call utils_assert_with_string(condition,error_message, &
            trim(value_1_as_string)//' '//&
            trim(value_2_as_string)//' '//&
            trim(value_3_as_string)//' '//&
            trim(value_4_as_string))
    end if

  end subroutine utils_assert_with_logicals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_assert_with_string(condition,error_message,optional_string)

    !==========================================================================!
    ! Asserts that 'condition' is .true., if not -- forwards 'error_message'   !
    ! optionally followed by 'optional_string' to utils_abort.                 !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   condition (input): The condition to check.                             !
    !   error_message (input): The error message to print if assertion fails.  !
    !   optional_string (input): The optional string to print if assertion     !
    !                            fails.                                        !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Avoid using conditions with side effects (such as function calls).     !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2010.                                   !
    !==========================================================================!

    implicit none

    ! Arguments
    logical, intent(in) :: condition
    character(len=*), intent(in) :: error_message
    character(len=*), intent(in), optional :: optional_string

    ! -----------------------------------------------------------------------

    if(condition) return

    if(present(optional_string)) then
       call utils_abort(error_message//' '//optional_string)
    else
       call utils_abort(error_message)
    end if

  end subroutine utils_assert_with_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_trace_in(msg)
    !=======================================================================!
    ! Enters a proc in the call graph, outputs the caller name to stdout.   !
    !-----------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 06/2010.                                 !
    !=======================================================================!

#ifdef TRACE
#ifdef MPI
    use comms, only: pub_my_proc_id
#endif
#endif
    use constants, only: stdout

    implicit none

    character(len=*), intent(in) :: msg

    ! -----------------------------------------------------------------------

    call utils_use_var(msg)

#ifdef TRACE
    write(stdout,*)
    write(stdout,'(a)',advance='no') repeat(' ', nindent)

#ifdef MPI
    write(stdout,'(3a,i7,a)',advance='no') '-->> ', msg, ' (proc #', &
         pub_my_proc_id,')'
#else
    write(stdout,*) '-->> ', msg
#endif

    nindent = nindent+4
    call utils_flush(stdout,.true.)

#endif

  end subroutine utils_trace_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_trace_out(msg)
    !=======================================================================!
    ! Exits a proc in the call graph, outputs the caller name to stdout.    !
    !-----------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 06/2010.                                 !
    !=======================================================================!

#ifdef TRACE
    use constants, only: stdout
    use comms, only: pub_my_proc_id
#endif

    implicit none

    character(len=*), intent(in) :: msg

    ! -----------------------------------------------------------------------

#ifdef TRACE
    nindent = nindent-4
    call utils_assert(nindent >= 0, 'utils_trace_out: negative indentation. &
         &This means you either have mispaired calls to utils_trace_in() and &
         &utils_trace_out() (e.g. early function return skipping over utils_&
         &trace_out()) or you attempted to use utils_trace_{in,out}() from &
         &an OMP region. This is not allowed because nindent is a module-wide &
         &variable vulnerable to race conditions.', nindent)
    write(stdout,*)
    write(stdout,'(a)',advance='no') repeat(' ', nindent)

#ifdef MPI
    write(stdout,'(3a,i7,a)',advance='no') '<<-- ', msg, ' (proc #', &
         pub_my_proc_id,')'
#else
    write(stdout,*) '<<-- ', msg
#endif

    call utils_flush(stdout,.true.)

#else
    call utils_use_var(msg)
#endif

  end subroutine utils_trace_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function utils_is_finite(x)
    !=========================================================================!
    ! Returns true if the double precision argument is an IEEE finite number. !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   x (input) : the double precision value to check                       !
    !-------------------------------------------------------------------------!
    ! Written by JM Escartin in September 2016.                               !
    !=========================================================================!

    ! If you find this to be the case, add "-DNO_IS_FINITE" to your compile options.
    ! This will cause is_infty to always return true.

    use constants, only: DP

#if defined(NO_IS_FINITE) || defined(__GFORTRAN__) && __GNUC__ < 5 || defined(__PGI) && __PGIC__ < 12
    implicit none
    real(kind=DP), intent(in) :: x

    utils_is_finite = .true.
    return
#else
    use, intrinsic :: IEEE_arithmetic, only: IEEE_is_finite !! External dependency
    ! gfortran note: For full compliance with the Fortran standards, code using
    ! the IEEE_EXCEPTIONS or IEEE_ARITHMETIC modules should be compiled with
    ! the following options: -fno-unsafe-math-optimizations -frounding-math
    ! -fsignaling-nans.

    implicit none

    real(kind=DP), intent(in) :: x

    utils_is_finite = IEEE_is_finite(x)
#endif

  end function utils_is_finite

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function utils_isnan(x)
    !=========================================================================!
    ! Returns .true. iff the double precision argument x is NaN.              !
    ! Allows abstracting away from the platform-specific mess.                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   x (input) : the double precision value to check against NaN.          !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 01/2011.                                   !
    !=========================================================================!

    ! jd: This implementation was found to work with the following compilers:
    !     - Intel's ifort 10.1, 11.1
    !     - GNU fortran 4.3, 4.4
    !     - PGI fortran 7.1, 10.9, 11.4
    !     - Sun Pro fortran 95 8.5
    !
    ! jd: It does not work (fails to compile) with
    !     - GNU fortran 4.1.2
    !
    ! If you find this to be the case, add "-DNO_ISNAN" to your compile options.
    ! This will cause isnan to always return false.

    use constants, only: DP, stderr

#if defined(NO_ISNAN)
    implicit none
    real(kind=DP), intent(in) :: x

    utils_isnan = .false.
    return
#else
#if defined(__GFORTRAN__) && __GNUC__ < 5 || defined(__PGI) && __PGIC__ < 12
#else
    use, intrinsic :: IEEE_arithmetic, only: IEEE_support_nan, IEEE_is_nan
    ! gfortran note: For full compliance with the Fortran standards, code using
    ! the IEEE_EXCEPTIONS or IEEE_ARITHMETIC modules should be compiled with
    ! the following options: -fno-unsafe-math-optimizations -frounding-math
    ! -fsignaling-nans.
#endif

    implicit none
    real(kind=DP), intent(in) :: x

#if defined(__PGI) && __PGIC__ < 12
    ! jd: PGI uses an externally linked isnand function
    ! jme: at least versions 12 and above have the new intrinsic function
    logical :: isnand
    utils_isnan = isnand(x)
#elif defined(__GFORTRAN__) && __GNUC__ < 5
    ! use C's isnan
    utils_isnan = isnan(x)
#else
    ! jme: new default
    if (IEEE_support_nan(x)) then
       utils_isnan = ieee_is_nan(x)
    else
       ! jme: this processor does not support IEEE NaNs
       write(stderr,'(a)') 'Warning: this processor does not support IEEE NaNs &
       &for real variables of the provided kind.'
       utils_isnan = .false.
    end if
#endif
#endif

  end function utils_isnan

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_sanity_check_real_3(arr, arrname, l1, l2, l3, excessive)
    !=========================================================================!
    ! Scans a 3D array of double precision numbers for NaNs, Inf, -Inf and    !
    ! excessively large values. If any are found, the location and value of   !
    ! the first one is printed and utils_abort is called.                     !
    !                                                                         !
    ! Exits immediately with no action if pub_debug is off.                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   arr (input)     : The array to check                                  !
    !   arrname (input) : The name of the array, written out if check fails.  !
    !   l1,l2,l3 (input, optional): Lower bounds of arr, only needed if the   !
    !                               array is not dimensioned from 1.          !
    !   excessive (input, opt): A threshold for deeming a value excessive.    !
    !                           If omitted, only NaNs and +/- Infs will be    !
    !                           detected.                                     !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 06/2010.                                   !
    !=========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: stderr

    implicit none

    ! jd: Input arguments
    real(kind=DP), intent(in)     :: arr(:,:,:)
    character(len=*), intent(in)  :: arrname
    integer, intent(in), optional :: l1, l2, l3
    real(kind=DP), intent(in), optional :: excessive

    integer :: i,j,k       ! jd: Indices
    integer :: ii,jj,kk    ! jd: Indices
    real(kind=DP) :: value ! jd: arr(i,j,k)
    logical :: failed      ! jd: .true. iff arr(i,j,k) failed sanity check
    integer :: myrank      ! jd: Alias for pub_my_proc_id

    !------------------------------------------------------------------------

    if (.not. loc_debug) return

    myrank = pub_my_proc_id ! jd: Local copy is accessible from debugger
    do k=lbound(arr,3), ubound(arr,3)
       do j=lbound(arr,2), ubound(arr,2)
          do i=lbound(arr,1), ubound(arr,1)

             value = arr(i,j,k)
             failed = .false.

             if(.not. utils_is_finite(value)) then
                write(stderr,*) 'Non-finite number detected in ', arrname
                failed = .true.
             end if
             if(utils_isnan(value)) then
                write(stderr,*) 'NaN detected in ', arrname
                failed = .true.
             end if
             if(present(excessive)) then
                if(abs(value) > excessive) then
                   write(stderr,*) 'An excessively large value detected in ', &
                        arrname
                   failed = .true.
                end if
             end if

             if(failed) then
                write(stderr,*) 'Element at ', i,' ',j,' ',k,' is :',arr(i,j,k)
                if(present(l1) .or. present(l2) .or. present(l3)) then
                   ii=i
                   jj=j
                   kk=k
                   if(present(l1)) ii=ii+l1-1
                   if(present(l2)) jj=jj+l2-1
                   if(present(l3)) kk=kk+l3-1
                   write(stderr,*) 'NB: The array might not be dimensioned from&
                        & 1 -- actually the element is at ',ii,' ',jj,' ',kk
                end if
                write(stderr,*) 'proc: ',myrank
                call utils_flush(stderr,.true.)
                call utils_abort('Sanity check failed for '//arrname)
             end if

          end do
       end do
    end do

  end subroutine utils_sanity_check_real_3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_sanity_check_real_2(arr, arrname, l1, l2, excessive)
    !=========================================================================!
    ! Scans a 2D array of double precision numbers for NaNs, Inf, -Inf and    !
    ! excessively large values. If any are found, the location and value of   !
    ! the first one is printed and utils_abort is called.                     !
    !                                                                         !
    ! Exits immediately with no action if pub_debug is off.                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   arr (input)     : The array to check                                  !
    !   arrname (input) : The name of the array, written out if check fails.  !
    !   l1,l2 (input, optional): Lower bounds of arr, only needed if the      !
    !                            array is not dimensioned from 1.             !
    !   excessive (input, opt): A threshold for deeming a value excessive.    !
    !                           If omitted, only NaNs and +/- Infs will be    !
    !                           detected.                                     !
    !-------------------------------------------------------------------------!
    ! Cloned from 3D version by Jacek Dziedzic in December 2015.              !
    !=========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: stderr

    implicit none

    ! jd: Input arguments
    real(kind=DP), intent(in)     :: arr(:,:)
    character(len=*), intent(in)  :: arrname
    integer, intent(in), optional :: l1, l2
    real(kind=DP), intent(in), optional :: excessive

    integer :: i,j         ! jd: Indices
    integer :: ii,jj       ! jd: Indices
    real(kind=DP) :: value ! jd: arr(i,j)
    logical :: failed      ! jd: .true. iff arr(i,j) failed sanity check
    integer :: myrank      ! jd: Alias for pub_my_proc_id

    !------------------------------------------------------------------------

    if (.not. loc_debug) return

    myrank = pub_my_proc_id ! jd: Local copy is accessible from debugger
    do j=lbound(arr,2), ubound(arr,2)
       do i=lbound(arr,1), ubound(arr,1)

          value = arr(i,j)
          failed = .false.

          if(.not. utils_is_finite(value)) then
             write(stderr,*) 'Non-finite number detected in ', arrname
             failed = .true.
          end if
          if(utils_isnan(value)) then
             write(stderr,*) 'NaN detected in ', arrname
             failed = .true.
          end if
          if(present(excessive)) then
             if(abs(value) > excessive) then
                write(stderr,*) 'An excessively large value detected in ', &
                     arrname
                failed = .true.
             end if
          end if

          if(failed) then
             write(stderr,*) 'Element at ', i,' ',j,' is :',arr(i,j)
             if(present(l1) .or. present(l2)) then
                ii=i
                jj=j
                if(present(l1)) ii=ii+l1-1
                if(present(l2)) jj=jj+l2-1
                write(stderr,*) 'NB: The array might not be dimensioned from&
                     & 1 -- actually the element is at ',ii,' ',jj
             end if
             write(stderr,*) 'proc: ',myrank
             call utils_flush(stderr,.true.)
             call utils_abort('Sanity check failed for '//arrname)
          end if

       end do
    end do

  end subroutine utils_sanity_check_real_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_sanity_check_real_1(arr, arrname, l1, excessive)
    !=========================================================================!
    ! Scans a 1D array of double precision numbers for NaNs, Inf, -Inf and    !
    ! excessively large values. If any are found, the location and value of   !
    ! the first one is printed and utils_abort is called.                     !
    !                                                                         !
    ! Exits immediately with no action if pub_debug is off.                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   arr (input)     : The array to check                                  !
    !   arrname (input) : The name of the array, written out if check fails.  !
    !   l1 (input, optional): Lower bound of arr, only needed if the          !
    !                               array is not dimensioned from 1.          !
    !   excessive (input, opt): A threshold for deeming a value excessive.    !
    !                           If omitted, only NaNs and +/- Infs will be    !
    !                           detected.                                     !
    !-------------------------------------------------------------------------!
    ! Cloned from 3D version by Jacek Dziedzic in December 2015.              !
    !=========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: stderr

    implicit none

    ! jd: Input arguments
    real(kind=DP), intent(in)     :: arr(:)
    character(len=*), intent(in)  :: arrname
    integer, intent(in), optional :: l1
    real(kind=DP), intent(in), optional :: excessive

    integer :: i           ! jd: Indices
    integer :: ii          ! jd: Indices
    real(kind=DP) :: value ! jd: arr(i)
    logical :: failed      ! jd: .true. iff arr(i) failed sanity check
    integer :: myrank      ! jd: Alias for pub_my_proc_id

    !------------------------------------------------------------------------

    if (.not. loc_debug) return

    myrank = pub_my_proc_id ! jd: Local copy is accessible from debugger
    do i=lbound(arr,1), ubound(arr,1)

       value = arr(i)
       failed = .false.

       if(.not. utils_is_finite(value)) then
          write(stderr,*) 'Non-finite number detected in ', arrname
          failed = .true.
       end if
       if(utils_isnan(value)) then
          write(stderr,*) 'NaN detected in ', arrname
          failed = .true.
       end if
       if(present(excessive)) then
          if(abs(value) > excessive) then
             write(stderr,*) 'An excessively large value detected in ', arrname
             failed = .true.
          end if
       end if

       if(failed) then
          write(stderr,*) 'Element at ', i,' is :',arr(i)
          if(present(l1)) then
             ii=i
             if(present(l1)) ii=ii+l1-1
             write(stderr,*) 'NB: The array might not be dimensioned from&
                  & 1 -- actually the element is at ',ii
          end if
          write(stderr,*) 'proc: ',myrank
          call utils_flush(stderr,.true.)
          call utils_abort('Sanity check failed for '//arrname)
       end if

    end do

  end subroutine utils_sanity_check_real_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_sanity_check_real_0(r, rname, excessive)
    !=========================================================================!
    ! Checks if a double precision number is a  NaNs, Inf, -Inf or a value    !
    ! with an excessive magnitude. If so, the value is printed and utils_abort!
    ! is called.                                                              !
    !                                                                         !
    ! Exits immediately with no action if pub_debug is off.                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   r (input): The value to check.                                        !
    !   rname (input): The name of the value, written out if the check fails. !
    !   excessive (input, opt): A threshold for deeming a value excessive.    !
    !                           If omitted, only NaNs and +/- Infs will be    !
    !                           detected.                                     !
    !-------------------------------------------------------------------------!
    ! Cloned from 1D version by Jacek Dziedzic in Octoberr 2016.              !
    !=========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: stderr

    implicit none

    ! jd: Input arguments
    real(kind=DP), intent(in)     :: r
    character(len=*), intent(in)  :: rname
    real(kind=DP), intent(in), optional :: excessive

    logical :: failed      ! jd: .true. iff arr(i) failed sanity check
    integer :: myrank      ! jd: Alias for pub_my_proc_id

    !------------------------------------------------------------------------

    if (.not. loc_debug) return

    myrank = pub_my_proc_id ! jd: Local copy is accessible from debugger

    failed = .false.

    if(.not. utils_is_finite(r)) then
       write(stderr,*) 'The value of ', rname,' is infinite: ', r
       failed = .true.
    end if
    if(utils_isnan(r)) then
       write(stderr,*) 'The value of ', rname, ' is NaN: ', r
       failed = .true.
    end if
    if(present(excessive)) then
       if(abs(r) > excessive) then
          write(stderr,*) 'The magnitude of ',rname, ' is excessively large: ',r
          failed = .true.
       end if
    end if

    if(failed) then
       write(stderr,*) 'proc: ',myrank
       call utils_flush(stderr,.true.)
       call utils_abort('Sanity check failed for '//rname)
    end if

  end subroutine utils_sanity_check_real_0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_dump_array3D_to_file(arr, filename, output_format)
    !==========================================================================!
    ! Writes the contents of a 3D local (non-distributed) array to a file      !
    ! for debugging purposes. The output filename is constructed by suffixing  !
    ! the passed filename with the proc number and a global running counter.   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   arr (input): the 3D array to write out.                                !
    !   filename (input): The basename of the file where arr will be written.  !
    !   output_format (input, optional): The format to use for output.         !
    !                                    If omitted, it defaults to (e20.12).  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, November 2010.                                !
    !==========================================================================!

    use constants, only: DP
    use comms, only: pub_my_proc_id

    implicit none

    ! jd: Input arguments
    real(kind=DP), intent(in)              :: arr(:,:,:)
    character(len=*), intent(in)           :: filename
    character(len=*), intent(in), optional :: output_format

    !------------------------------------------------------------------------

    ! jd: Local variables
    integer :: outunit
    integer, save :: file_counter = 0
    character(len=200) :: outfilename
    character(len=200) :: suffix1, suffix2
    integer :: i,j,k

    ! jd: Construct the full filename
    write(suffix1,'(i10)') pub_my_proc_id
    suffix1 = adjustl(suffix1)
    write(suffix2,'(i10)') file_counter
    suffix2 = adjustl(suffix2)
    outfilename=trim(trim(trim(filename)//'_proc'//suffix1)//'_'//suffix2)

    ! jd: Open the file
    outunit = utils_unit()
    open(outunit, file=outfilename, err=10)

    ! jd: Output the contents of the array
    do k=lbound(arr,3), ubound(arr,3)
       do j=lbound(arr,2), ubound(arr,2)
          do i=lbound(arr,1), ubound(arr,1)
             if(present(output_format)) then
                write(outunit,output_format) arr(i,j,k)
             else
                write(outunit,'(e20.12)') arr(i,j,k)
             end if
          end do
       end do
    end do

    ! jd: Close the file, increment running counter
    close(outunit, err=20)
    file_counter = file_counter + 1

    return

    ! jd: I/O error handling
10  call utils_abort('Error during creation of file: '//outfilename)
20  call utils_abort('Error during closing of file: '//outfilename)

  end subroutine utils_dump_array3D_to_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_dump_array2D_to_file(arr, filename, output_format)
    !==========================================================================!
    ! Writes the contents of a 2D local (non-distributed) array to a file      !
    ! for debugging purposes. The output filename is constructed by suffixing  !
    ! the passed filename with the proc number and a global running counter.   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   arr (input): the 2D array to write out.                                !
    !   filename (input): The basename of the file where arr will be written.  !
    !   output_format (input, optional): The format to use for output.         !
    !                                    If omitted, it defaults to (e20.12).  !
    !--------------------------------------------------------------------------!
    ! Cloned from utils_dump_array3D_to_file by Jacek Dziedzic, April 2012.    !
    !==========================================================================!

    use constants, only: DP
    use comms, only: pub_my_proc_id

    implicit none

    ! jd: Input arguments
    real(kind=DP), intent(in)              :: arr(:,:)
    character(len=*), intent(in)           :: filename
    character(len=*), intent(in), optional :: output_format

    !------------------------------------------------------------------------

    ! jd: Local variables
    integer :: outunit
    integer, save :: file_counter = 0
    character(len=200) :: outfilename
    character(len=200) :: suffix1, suffix2
    integer :: i,j

    ! jd: Construct the full filename
    write(suffix1,'(i10)') pub_my_proc_id
    suffix1 = adjustl(suffix1)
    write(suffix2,'(i10)') file_counter
    suffix2 = adjustl(suffix2)
    outfilename=trim(trim(trim(filename)//'_proc'//suffix1)//'_'//suffix2)

    ! jd: Open the file
    outunit = utils_unit()
    open(outunit, file=outfilename, err=10)

    ! jd: Output the contents of the array
    do j=lbound(arr,2), ubound(arr,2)
       do i=lbound(arr,1), ubound(arr,1)
          if(present(output_format)) then
             write(outunit,output_format) arr(i,j)
          else
             write(outunit,'(e20.12)') arr(i,j)
          end if
       end do
    end do

    ! jd: Close the file, increment running counter
    close(outunit, err=20)
    file_counter = file_counter + 1

    return

    ! jd: I/O error handling
10  call utils_abort('Error during creation of file: '//outfilename)
20  call utils_abort('Error during closing of file: '//outfilename)

  end subroutine utils_dump_array2D_to_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_dump_array1D_to_file(arr, filename, output_format, no_suffix)
    !==========================================================================!
    ! Writes the contents of a 1D local (non-distributed) array to a file      !
    ! for debugging purposes. The output filename is constructed by suffixing  !
    ! the passed filename with the proc number and a global running counter.   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   arr (input): the 2D array to write out.                                !
    !   filename (input): The basename of the file where arr will be written.  !
    !   output_format (input, optional): The format to use for output.         !
    !                                    If omitted, it defaults to (e20.12).  !
    !   no_suffix (input, optional): If present and true, the proc number and  !
    !                                global running counter are not suffixed.  !
    !--------------------------------------------------------------------------!
    ! Cloned from utils_dump_array2D_to_file by Jacek Dziedzic, August 2012.   !
    !==========================================================================!

    use constants, only: DP
    use comms, only: pub_my_proc_id

    implicit none

    ! jd: Input arguments
    real(kind=DP), intent(in)              :: arr(:)
    character(len=*), intent(in)           :: filename
    character(len=*), intent(in), optional :: output_format
    logical, intent(in), optional          :: no_suffix

    ! jd: Local variables
    integer            :: outunit
    integer, save      :: file_counter = 0
    logical            :: want_suffix
    character(len=200) :: outfilename
    character(len=200) :: suffix1, suffix2
    integer            :: i

    !------------------------------------------------------------------------

    ! jd: Construct the full filename
    want_suffix = .true.
    if(present(no_suffix)) then
       if(no_suffix) then
          want_suffix = .false.
       end if
    end if
    if(want_suffix) then
       write(suffix1,'(i10)') pub_my_proc_id
       suffix1 = adjustl(suffix1)
       write(suffix2,'(i10)') file_counter
       suffix2 = adjustl(suffix2)
       outfilename=trim(trim(trim(filename)//'_proc'//suffix1)//'_'//suffix2)
    else
       outfilename=trim(filename)
    end if

    ! jd: Open the file
    outunit = utils_unit()
    open(outunit, file=outfilename, err=10)

    ! jd: Output the contents of the array
    do i=lbound(arr,1), ubound(arr,1)
       if(present(output_format)) then
          write(outunit,output_format) arr(i)
       else
          write(outunit,'(e20.12)') arr(i)
       end if
    end do

    ! jd: Close the file, increment running counter
    close(outunit, err=20)
    file_counter = file_counter + 1

    return

    ! jd: I/O error handling
10  call utils_abort('Error during creation of file: '//outfilename)
20  call utils_abort('Error during closing of file: '//outfilename)

  end subroutine utils_dump_array1D_to_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_use_var_integer(x)
    !==========================================================================!
    ! Pretends to use its argument, so that compiler warnings about unused     !
    ! variables can be suppressed.                                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   x (input): An integer value to be 'used'.                              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2011.                                   !
    !==========================================================================!

    implicit none

    integer, intent(in) :: x
    integer :: dummy
    dummy = x

  end subroutine utils_use_var_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_use_var_logical(x)
    !==========================================================================!
    ! Pretends to use its argument, so that compiler warnings about unused     !
    ! variables can be suppressed.                                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   x (input): A logical value to be 'used'.                               !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2011.                                   !
    !==========================================================================!

    implicit none

    logical, intent(in) :: x
    logical :: dummy
    dummy = x

  end subroutine utils_use_var_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_use_var_real(x)
    !==========================================================================!
    ! Pretends to use its argument, so that compiler warnings about unused     !
    ! variables can be suppressed.                                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   x (input): An integer value to be 'used'.                              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2011.                                   !
    !==========================================================================!

    use constants, only: DP

    implicit none

    real(kind=DP), intent(in) :: x
    real(kind=DP) :: dummy
    dummy = x

  end subroutine utils_use_var_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_use_var_character(x)
    !==========================================================================!
    ! Pretends to use its argument, so that compiler warnings about unused     !
    ! variables can be suppressed.                                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   x (input): A character value to be 'used'.                             !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2011.                                   !
    !==========================================================================!

    implicit none

    character(len=1), intent(in) :: x
    character(len=1) :: dummy
    dummy = x

  end subroutine utils_use_var_character

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_qc_print_integer(tag, value)
    !==========================================================================!
    ! Prints a line of QC output to stdout.                                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   tag (input): A tag (string) describing the QC output.                  !
    !   value (input): The value corresponding to the QC output.               !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, November 2012.                                !
    !==========================================================================!

    use constants, only: stdout

    implicit none

    character(len=*), intent(in) :: tag
    integer, intent(in)    :: value

    character(len=24) :: outtag

    ! --------------------------------------------------------------------------

    outtag = '['//trim(tag)//']'

    write(stdout,'(a30,i23)') '<QC> '//adjustr(outtag)//':', value

  end subroutine utils_qc_print_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_qc_print_real(tag, value, minimum_magnitude)
    !==========================================================================!
    ! Prints a line of QC output to stdout.                                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   tag (input): A tag (string) describing the QC output.                  !
    !   value (input): The value corresponding to the QC output.               !
    !   minimum_magnitude (input, opt): If provided, sets a threshold for the  !
    !                                   absolute value of quantities to be     !
    !                                   printed.                               !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, November 2012.                                !
    ! Extended with minimum_magnitude by Jacek Dziedzic, October 2019.         !
    !==========================================================================!

    use constants, only: stdout, DP

    implicit none

    character(len=*), intent(in)        :: tag
    real(kind=DP), intent(in)           :: value
    real(kind=DP), intent(in), optional :: minimum_magnitude

    character(len=24) :: outtag
    real(kind=DP) :: loc_minimum_magnitude

    ! --------------------------------------------------------------------------

    if(present(minimum_magnitude)) then
       loc_minimum_magnitude = minimum_magnitude
    else
       loc_minimum_magnitude = 0.0_DP
    end if

    if(abs(value) >= loc_minimum_magnitude) then
       outtag = '['//trim(tag)//']'

       write(stdout,'(a30,f23.12)') '<QC> '//adjustr(outtag)//':', value
    end if

  end subroutine utils_qc_print_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_qc_print_string(tag, value)
    !==========================================================================!
    ! Prints a line of QC output to stdout.                                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   tag (input): A tag (string) describing the QC output.                  !
    !   value (input): The value corresponding to the QC output.               !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, May 2014.                                     !
    !==========================================================================!

    use constants, only: stdout

    implicit none

    character(len=*), intent(in) :: tag
    character(len=*), intent(in) :: value

    character(len=24) :: outtag

    ! --------------------------------------------------------------------------

    outtag = '['//trim(tag)//']'

    write(stdout,'(a30,a50)') '<QC> '//adjustr(outtag)//':', trim(value)

  end subroutine utils_qc_print_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_write_csv(mat, filename, nelems, format, proc_id)

    !========================================================================!
    ! ars: This subroutine takes a 2D array and writes a csv spreadsheet.    !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                           !
    !========================================================================!

    use constants, only: DP, stdout

    implicit none

    ! ars: Arguments
    integer, intent(in) :: nelems
    real(kind=DP) :: mat(nelems,nelems)
    character(len=*), intent(in) :: filename
    character(len=*) :: format
    integer, optional, intent(in) :: proc_id

    ! ars: local variables
    integer :: myunit, jjj
    character(len=80) :: fmt, pname, loc_filename


    if (present(proc_id)) then
       write(pname, '(i10)') proc_id
       loc_filename=trim(adjustl(filename))//"_"//trim(adjustl(pname))//".csv"
    else
       loc_filename=trim(adjustl(filename))//".csv"
    end if

    write(stdout,*) "Writing file ", loc_filename

    ! ars: get format
    write(fmt,'(a1)') '('
    write(fmt,'(a1,i5)') trim(adjustl(fmt)), nelems
    write(fmt,'(a,a1, a, a6)') trim(adjustl(fmt)), '(', trim(adjustl(format)), '";") )'

    myunit = utils_unit()
    open(unit=myunit, file=trim(adjustl(loc_filename)), action="write")

    do jjj = 1, nelems
       write(myunit, fmt) mat(jjj,:)
    end do

    close(myunit)


  end subroutine utils_write_csv


  !.............................................................................


  subroutine utils_read_csv(mat, filename, nelems, format)

    !========================================================================!
    ! ars: This subroutine reads a 2D csv spreadsheet and places the data in !
    !      a Fortran 2D array.                                               !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                           !
    !========================================================================!

    use constants, only: DP, stdout

    implicit none

    ! ars: Arguments
    integer, intent(in) :: nelems
    real(kind=DP), intent(inout) :: mat(:,:)
    character(len=*), intent(in) :: filename
    character(len=*) :: format

    ! ars: local variables
    integer :: myunit, jjj
    character(len=80) :: fmt

    write(stdout,*) "Reading file ", filename

    ! ars: get format
    write(fmt,'(a1)') '('
    write(fmt,'(a1,i5)') trim(adjustl(fmt)), nelems
    write(fmt,'(a,a1, a, a6)') trim(adjustl(fmt)), '(', trim(adjustl(format)), '";") )'

    myunit = utils_unit()
    open(unit=myunit, file=trim(adjustl(filename)), action="read")

    do jjj = 1, nelems
       read(myunit, fmt) mat(jjj,:)
    end do

    close(myunit)

  end subroutine utils_read_csv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function utils_devel_code_logical(default, blockname, keyword, &
       devel_code, no_bcast, no_warn)

    !========================================================================!
    ! ars: This function reads a logical variable included in the devel_code !
    !      input character array.                                            !
    ! ars: Copied from similar code by Simon Dubois.                         !
    !------------------------------------------------------------------------!
    ! ars: syntax: devel_code BLOCKNAME:KEYWORD=T/F:BLOCKNAME                !
    !------------------------------------------------------------------------!
    ! ars: usage notes:                                                      !
    !      default = default output value of the function if the keyword is  !
    !                not specified in devel_code.                            !
    !      blockname must be passed into this function with CAPITAL LETTERS  !
    !      keyword must be passed into this function with CAPITAL LETTERS    !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                           !
    !========================================================================!

    use comms, only: pub_on_root, pub_root_proc_id, comms_bcast
    use constants, only: stdout

    implicit none

    logical, intent(in) :: default
    character(len=*), intent(in) :: blockname, keyword, devel_code
    logical, optional, intent(in) :: no_bcast
    logical, optional, intent(in) :: no_warn

    integer :: start_pos, stop_pos, test_pos
    logical :: loc_no_bcast, loc_no_warn

    if (present(no_bcast)) then
       loc_no_bcast = no_bcast
    else
       loc_no_bcast = .false.
    end if

    if (present(no_warn)) then
       loc_no_warn = no_warn
    else
       loc_no_warn = .false.
    end if

    if (pub_on_root .or. loc_no_bcast) then

       utils_devel_code_logical = default

       if (len_trim(devel_code)>0) then
          start_pos=index(devel_code,trim(adjustl(blockname))//':')
          stop_pos=index(devel_code,':'//trim(adjustl(blockname)))
          if (stop_pos<=0) stop_pos=len_trim(devel_code) !missing end so go to end of string

          if (start_pos>0) then

             ! ars: set finite differences scheme
             test_pos=index(devel_code,trim(adjustl(keyword))//'=')

             if (test_pos>start_pos.and.test_pos<stop_pos) then
                test_pos=test_pos+len(trim(adjustl(keyword))//'=')

                read(devel_code(test_pos:test_pos+ &
                     & index(devel_code(test_pos:stop_pos),':')-2),*) &
                     utils_devel_code_logical

                ! ars: write warning on screen
                if(pub_on_root .and. .not. loc_no_warn) then
                   write(stdout,'(a,a,l1)') "Warning: devel_code: ", &
                        trim(adjustl(blockname)//':')//trim(adjustl(keyword))&
                        //' = ', utils_devel_code_logical
                end if

             end if
          end if
       end if
    end if

    if (.not.loc_no_bcast) then
       call comms_bcast(pub_root_proc_id, utils_devel_code_logical)
    end if

  end function utils_devel_code_logical

  !.............................................................................

  integer function utils_devel_code_integer(idefault, blockname, keyword, &
       devel_code, no_bcast, no_warn)

    !=========================================================================!
    ! ars: This function reads an integer variable included in the devel_code !
    !      input character array.                                             !
    ! ars: Copied from similar code by Simon Dubois.                          !
    !-------------------------------------------------------------------------!
    ! ars: syntax: devel_code BLOCKNAME:KEYWORD=I:BLOCKNAME                   !
    !-------------------------------------------------------------------------!
    ! ars: usage notes:                                                       !
    !      default = default output value of the function if the keyword is   !
    !                not specified in devel_code.                             !
    !      blockname must be passed into this function with CAPITAL LETTERS   !
    !      keyword must be passed into this function with CAPITAL LETTERS     !
    !-------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                            !
    !=========================================================================!

    use comms, only: pub_on_root, pub_root_proc_id, comms_bcast
    use constants, only: stdout

    implicit none

    integer, intent(in) :: idefault
    character(len=*), intent(in) :: blockname, keyword, devel_code
    logical, optional, intent(in) :: no_bcast
    logical, optional, intent(in) :: no_warn

    integer :: start_pos, stop_pos, test_pos
    logical :: loc_no_bcast, loc_no_warn

    if (present(no_bcast)) then
       loc_no_bcast = no_bcast
    else
       loc_no_bcast = .false.
    end if

    if (present(no_warn)) then
       loc_no_warn = no_warn
    else
       loc_no_warn = .false.
    end if

    if (pub_on_root .or. loc_no_bcast) then

       utils_devel_code_integer = idefault

       if (len_trim(devel_code)>0) then
          start_pos=index(devel_code,trim(adjustl(blockname))//':')
          stop_pos=index(devel_code,':'//trim(adjustl(blockname)))
          if (stop_pos<=0) stop_pos=len_trim(devel_code) !missing end so go to end of string
          if (start_pos>0) then

             ! ars: set finite differences scheme
             test_pos=index(devel_code,trim(adjustl(keyword))//'=')
             if (test_pos>start_pos.and.test_pos<stop_pos) then
                test_pos=test_pos+len(trim(adjustl(keyword))//'=')
                read(devel_code(test_pos:test_pos+ &
                     & index(devel_code(test_pos:stop_pos),':')-2),*) &
                     utils_devel_code_integer

                ! ars: write warning on screen
                if(pub_on_root .and. .not. loc_no_warn) then
                   write(stdout,'(a,a,i10)') "Warning: devel_code: ", &
                        trim(adjustl(blockname)//':')//trim(adjustl(keyword))&
                        //' = ', utils_devel_code_integer
                end if

             end if
          end if
       end if
    end if

    if (.not.loc_no_bcast) then
       call comms_bcast(pub_root_proc_id,utils_devel_code_integer)
    end if

  end function utils_devel_code_integer

  !.............................................................................

  real(kind=DP) function utils_devel_code_real(default, blockname, keyword, &
       devel_code, no_bcast, no_warn)

    !=========================================================================!
    ! ars: This function reads an real variable included in the devel_code    !
    !      input character array.                                             !
    ! ars: Copied from similar code by Simon Dubois.                          !
    !-------------------------------------------------------------------------!
    ! ars: syntax: devel_code BLOCKNAME:KEYWORD=R:BLOCKNAME                   !
    !-------------------------------------------------------------------------!
    ! ars: usage notes:                                                       !
    !      default = default output value of the function if the keyword is   !
    !                not specified in devel_code.                             !
    !      blockname must be passed into this function with CAPITAL LETTERS   !
    !      keyword must be passed into this function with CAPITAL LETTERS     !
    !-------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                            !
    !=========================================================================!

    use comms, only: pub_on_root, pub_root_proc_id, comms_bcast
    use constants, only: stdout, DP

    implicit none

    real(kind=DP), intent(in) :: default
    character(len=*), intent(in) :: blockname, keyword, devel_code
    logical, optional, intent(in) :: no_bcast
    logical, optional, intent(in) :: no_warn

    integer :: start_pos, stop_pos, test_pos
    logical :: loc_no_bcast, loc_no_warn

    if (present(no_bcast)) then
       loc_no_bcast = no_bcast
    else
       loc_no_bcast = .false.
    end if

    if (present(no_warn)) then
       loc_no_warn = no_warn
    else
       loc_no_warn = .false.
    end if

    if (pub_on_root .or. loc_no_bcast) then

       utils_devel_code_real = default

       if (len_trim(devel_code)>0) then
          start_pos=index(devel_code,trim(adjustl(blockname))//':')
          stop_pos=index(devel_code,':'//trim(adjustl(blockname)))
          if (stop_pos<=0) stop_pos=len_trim(devel_code) !missing end so go to end of string
          if (start_pos>0) then

             ! ars: set finite differences scheme
             test_pos=index(devel_code,trim(adjustl(keyword))//'=')
             if (test_pos>start_pos.and.test_pos<stop_pos) then
                test_pos=test_pos+len(trim(adjustl(keyword))//'=')
                read(devel_code(test_pos:test_pos+ &
                     & index(devel_code(test_pos:stop_pos),':')-2),*) &
                     utils_devel_code_real

                ! ars: write warning on screen
                if(pub_on_root .and. .not. loc_no_warn) then
                   write(stdout,'(a,a,f22.12)') "Warning: devel_code: ", &
                        trim(adjustl(blockname)//':')//trim(adjustl(keyword))&
                        //' = ', utils_devel_code_real
                end if

             end if
          end if
       end if
    end if

    if (.not.loc_no_bcast) then
       call comms_bcast(pub_root_proc_id,utils_devel_code_real)
    end if

  end function utils_devel_code_real

!.............................................................................

  subroutine utils_busy_wait(delta_t)
    !==========================================================================!
    ! Busy-waits for a number of milliseconds.                                 !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   delta_t (in): The number of milliseconds to busy wait.                 !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Likely misbehaves if a leap second happens while waiting.              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, July 2012.                                    !
    ! Largely rewritten by Jose M Escartin, 11 February 2017.                  !
    !==========================================================================!

    use constants, only: LONG
    implicit none

    ! jd: Arguments
    integer, intent(in) :: delta_t

    ! jd: Internal variables
    integer, dimension(8) :: vals1, vals2
    integer(kind=LONG)    :: time1, time2

    !------------------------------------------------------------------------

    call date_and_time(values=vals1)
    time1 = internal_vals_to_int(vals1)
    time2 = time1

    do while(time2-time1 <= delta_t)
       call date_and_time(values=vals2)
       time2 = internal_vals_to_int(vals2)
       if(time2 < time1) exit ! jd: overflow, better end the wait than hang
    end do

    contains

       function internal_vals_to_int(vals) result(time)
          !============================================!
          ! Returns number of milliseconds since EPOCH !
          ! (1st January 2017 00:00:00 +00:00)         !
          ! excluding leap seconds.                    !
          !============================================!
          implicit none

          ! Argument
          integer, dimension(8), intent(in) :: vals
          ! Result
          integer(kind=LONG)                :: time

          ! Local variables
          integer, parameter :: EPOCH = 2017
          integer, parameter, dimension(12) :: days_in_month = &
                                           [31,28,31,30,31,30,31,31,30,31,30,31]
          integer :: y, day_of_year, days_since_EPOCH

          ! Day of the year (New Year is day zero).
          day_of_year = sum(days_in_month(1:(vals(2)-1))) + (vals(3)-1)
          if (internal_year_is_leap(vals(1)) .and. (vals(2)>2)) then
             day_of_year = day_of_year + 1
          end if

          ! Count complete days since beginning of EPOCH.
          days_since_EPOCH = day_of_year
          do y = EPOCH, vals(1)-1
             if (internal_year_is_leap(y)) then
                days_since_EPOCH = days_since_EPOCH + 366
             else
                days_since_EPOCH = days_since_EPOCH + 365
             end if
          end do

          ! Convert to milliseconds and add milliseconds within day.
          ! (takes into account timezone).
          time =  vals(8) + &
               1000 * (vals(7) + 60 * (vals(6)-vals(4)) + 3600 * vals(5)) + &
               int(days_since_EPOCH,kind=LONG) * (24 * 3600 * 1000)

       end function internal_vals_to_int

       function internal_year_is_leap(year) result(is_leap)
          !===========================================!
          ! Returns whether year is leap year or not. !
          !===========================================!
          implicit none
          integer, intent(in) :: year     ! Argument
          logical             :: is_leap  ! Result
          is_leap = (mod(year,4)==0) .and. ( .not. ( &
               mod(year,100)==0 .and. (.not. (mod(year,400)==0 ) ) ) )
       end function

  end subroutine utils_busy_wait

!.............................................................................

  recursive subroutine utils_check_stack_size(dummy)
    !==========================================================================!
    ! Attempts to create a stack-based temporary to force a crash early in the !
    ! run if the stack is too small. Currently assumes 100MB per MPI rank is   !
    ! sufficient.                                                              !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   Do not use the dummy argument.                                         !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, October 2013.                                 !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout

    implicit none

    real(kind=DP), intent(in), optional :: dummy(:)

    real(kind=DP), allocatable :: heap_dummy(:)
    integer :: ierr
    integer :: error_unit
    character(len=256) :: error_message_file
    integer, parameter :: reasonable_stack_size = 100000000 ! 100MB

    !------------------------------------------------------------------------
#ifndef WIN32
    if(.not. present(dummy)) then ! call from onetep.F90
       allocate(heap_dummy(reasonable_stack_size/8),stat=ierr)
       call utils_alloc_check('utils_check_stack_size','heap_dummy',ierr)
       heap_dummy(:) = 0.0_DP ! even if not referenced, an intent(in) arg needs init
       if(pub_on_root) then
          write(stdout,'(a)') 'If your calculation crashes here, before &
               &"Checking processes and threads..."'
          write(stdout,'(a)') 'is displayed, then your stack size is insufficient.'
          write(stdout,'(a)') 'If so, use "ulimit -s unlimited" at runtime.'
          write(stdout,'(a/)') 'You can disable this check with "check_stack_size F".'
       end if
       call utils_flush

       ! Create an .error_message file, it will be immediately deleted
       ! if the stack is of sufficient size.
       if(pub_on_root) then
          write(error_message_file,'(a,a)') trim(pub_rootname_local), &
               '.error_message'
          error_unit = utils_unit()
          open(error_unit, file=error_message_file, status="new", err=1000)
          write(error_unit,'(a)') 'Your stack size is insufficient. &
               &Use "ulimit -s unlimited".'
          close(error_unit)
1000   continue
       end if

       ! This creates the stack-based temporary and will crash if the
       ! stack size is insufficient
       call utils_check_stack_size(heap_dummy + heap_dummy)

       ! Delete the .error_message file
       if(pub_on_root) then
          error_unit = utils_unit()
          open(unit=error_unit, file=error_message_file, status="old", err=500)
          close(unit=error_unit, status="delete", err=500)
500    continue ! no-op if file was not there to begin with
       end if

       deallocate(heap_dummy,stat=ierr)
       call utils_dealloc_check('utils_check_stack_size','heap_dummy',ierr)

    else ! pretend to use the data so that the call is not optimised away
       call utils_use_var(dummy(1))
    end if
#endif
  end subroutine utils_check_stack_size

!.............................................................................

  subroutine utils_report_memory_estimate(functionality_name, &
       names, memory_requirements, opt_computed_total)
    !==========================================================================!
    ! Create a consistent pretty report of the estimated memory useage of a    !
    ! module. Memory useage must be inputted as number of bytes. This is then  !
    ! automatically formatted into kB/MB/GB etc.                               !
    !                                                                          !
    ! Arguments:                                                               !
    !    functionality_name: the header of the table to be printed.            !
    !    names, memory_requirements: entries in the table.                     !
    !    opt_computed_total (opt, out): If passed, will be set to the sum of   !
    !                                   memory requirements. Note that this is !
    !                                   only computed and returned on root.    !
    !                                   Other procs return -1.                 !
    !--------------------------------------------------------------------------!
    ! Caveat:                                                                  !
    !    If you use an array constructor (the syntax (/.../)) for passing the  !
    !    names, be aware that all entries *must* have the same length. It's up !
    !    to you to choose the length, but be sure to pad all the names that    !
    !    are shorter. The intel compiler will let you get away with this, but  !
    !    gfortran and NAG will (rightfully) complain.                          !
    !--------------------------------------------------------------------------!
    ! Written by Robert Bell, November 2013 using code by Louis Lee            !
    ! Generalized to any number of entries by Jacek Dziedzic in April 2014.    !
    ! Extended with 'opt_computed_total' by Jacek Dziedzic in May 2020.        !
    !==========================================================================!
     use comms, only: pub_on_root
     use constants, only: stdout, LONG

     implicit none

     ! arguments
     character(len=*), intent(in)              :: functionality_name
     character(len=*), intent(in)              :: names(:)
     integer(kind=LONG), intent(in)            :: memory_requirements(:)
     integer(kind=LONG), intent(out), optional :: opt_computed_total

     ! internal
     integer(kind=LONG) :: total
     integer            :: n_entries, i
     integer, parameter :: fmt_str_len = 35
     ! formatted names
     character(len=49) :: fmt_func_name
     character(len=fmt_str_len) :: fmt_str

     write(fmt_func_name,*) functionality_name

     total = -1_LONG
     if (pub_on_root) then
        write(stdout,'(1x,a)') &
             & '+----------------------------------------------------+'
        write(stdout,'(1x,a,a49,a)') &
             & '|  ', adjustl(fmt_func_name),                     ' |'
        write(stdout,'(1x,a)')  &
             & '|  Estimated memory requirement per MPI rank         |'
        write(stdout,'(1x,a)') &
             & '+----------------------------------------------------+'

        call utils_assert(.true., &
             'utils_report_memory_requirements: Inconsistent arguments passed')

        n_entries = ubound(names,1)
        total = 0_LONG

        do i = 1, n_entries
           call utils_assert(memory_requirements(i) >= 0,&
                'utils_report_memory_estimate: Supplied memory estimate for "'&
                //trim(names(i))//'" is negative. It''s likely you encountered &
                &integer wraparound *before* calling utils_report_memory_&
                &estimate (which itself uses LONGs).')

           if(memory_requirements(i) > 0) then
              call utils_assert(len(names(i)) < fmt_str_len, &
                   'utils_report_memory_requirements: Supplied name too long.',&
                   len(names(i)),fmt_str_len)
              write(fmt_str,*) names(i)
              write(stdout,'(1x,a,a35,3a)') &
                   '|  ',adjustl(fmt_str), ' : ', &
                   format_estimate(memory_requirements(i)),'  |'
           end if

           total = total + memory_requirements(i)
        end do

        write(stdout,'(1x,a)') &
             & '+----------------------------------------------------+'
        ! overall summary
        if(n_entries > 1) then
          fmt_str = 'Estimated peak total per MPI rank  '
          write(stdout,'(1x,a,a35,3a)') &
               & '|  ', adjustl(fmt_str), ' : ', format_estimate(total),'  |'
          write(stdout,'(1x,a)') &
               & '+----------------------------------------------------+'
        endif

     endif

     if(present(opt_computed_total)) then
        opt_computed_total = total
     end if

     ! ------------------------------------------------------------------------

     contains
        character(len=10) function format_estimate(memory)

           use constants, only: DP, LONG

           implicit none

           ! arguments
           integer(kind=LONG), intent(in) :: memory

           ! internal
           character(len=2) :: byte_mult(0:6)
           data byte_mult /' B','kB','MB','GB','TB','PB','EB'/
           integer, parameter :: two_pow_ten = 1024
           integer :: temp
           real(kind=DP) :: loc_memory

           loc_memory = memory

           temp = 0
           do while (loc_memory >= two_pow_ten .and. temp <= 6)
             loc_memory = loc_memory/real(two_pow_ten,kind=DP)
             temp = temp + 1
           end do
           write(format_estimate,'(f7.2,1x,2a)') loc_memory, byte_mult(temp)
        end function format_estimate

  end subroutine utils_report_memory_estimate

!.............................................................................

   integer(kind=LONG) function utils_fnv1_hash(string,previous_hash)
    !==========================================================================!
    ! Computes the Fowler-Noll-Vo (FNV-1) hash of the string. Multiple strings !
    ! can be hashed together by providing the result of the previous hash.     !
    ! An initial hash can be generated by calling it on the empty string.      !
    !--------------------------------------------------------------------------!
    ! Written by Robert Bell, 23/07/2014.                                      !
    !==========================================================================!

      use constants, only: LONG

      implicit none

      ! arguments
      character(len=*), intent(in) :: string
      integer(kind=LONG), intent(in), optional :: previous_hash

      ! internal
      integer(kind=LONG), parameter :: FNV_OFFSET_BASIS = 2166136261_LONG
      integer(kind=LONG), parameter :: FNV_PRIME = 16777619_LONG
      integer :: ii
      character(len=1) :: c, prev_c


      if (present(previous_hash)) then
         utils_fnv1_hash = previous_hash
      else
         utils_fnv1_hash = FNV_OFFSET_BASIS
      endif

      prev_c = ' ' ! init to remove any whitespace at start of string
      do ii = 1, len(string)
         c = string(ii:ii)

         ! skip multiple whitespace
         if (prev_c == ' ' .and. c == ' ') cycle
         prev_c = c

         utils_fnv1_hash = ieor(utils_fnv1_hash, ichar(c,kind=LONG))
         utils_fnv1_hash = utils_fnv1_hash * FNV_PRIME
      enddo


   end function utils_fnv1_hash

!.............................................................................

   integer(kind=LONG) function utils_fnv1_hash_file(filename)
    !==========================================================================!
    ! Computes the Fowler-Noll-Vo (FNV-1) hash of the contents of a file.      !
    !--------------------------------------------------------------------------!
    ! Written by Robert Bell, 23/07/2014.                                      !
    !==========================================================================!

      use constants, only: LONG

      implicit none

      ! arguments
      character(len=*), intent(in) :: filename

      ! internal
      integer :: iunit, ierr
      character(len=128) :: line

      ! initialise hash
      utils_fnv1_hash_file = utils_fnv1_hash('')

      iunit = utils_unit()

      open(unit=iunit,file=trim(filename),action='READ',iostat=ierr)
      call utils_open_unit_check('utils_fnv1_hash_file',trim(filename),ierr)

      ! calculate the hash for each line of the file
      do
         read(iunit,'(a)',end=101,err=103,iostat=ierr) line
         utils_fnv1_hash_file = &
              utils_fnv1_hash(line,utils_fnv1_hash_file)
      enddo

101   close(iunit,iostat=ierr)
      call utils_close_unit_check('utils_fnv1_hash_file',trim(filename),ierr)

      return

103   call utils_abort('Error in utils_fnv1_hash: error reading from file '//&
                       &trim(filename)//'; returned code ', ierr)

   end function utils_fnv1_hash_file

!.............................................................................

  subroutine utils_rand_gen(rand_array,num_rand, seed)
    !============================================================================!
    ! Very simple Linear congruental random number generator that will generate  !
    ! num_rand random numbers and store them in rand_array. The string of random !
    ! numbers is initalised with seed "seed", the multplyer used is 214013, the  !
    ! modulus is 2^31 and the increment is 2531011                               !
    ! The random numbers are normalised to lie between 0 and 1. Note that this   !
    ! is a very simple random number generator with a periodicity of around      !
    ! 31000, so it should not be used anywhere where the quality of the random   !
    ! numbers is very important.                                                 !
    !----------------------------------------------------------------------------!
    ! Written by Tim Zuehlsdorff in July 2014                                    !
    !============================================================================!

    use constants, only: LONG

    implicit none

    ! Arguments
    integer, intent(in) :: num_rand
    real(kind=DP), intent(inout) :: rand_array(num_rand)
    integer, intent(in) :: seed

    ! local variables:
    integer(kind=LONG) :: temp_val
    integer(kind=LONG) :: temp_val_prev
    integer(kind=LONG), parameter :: a=214013_LONG
    integer(kind=LONG), parameter :: c=2531011_LONG
    integer(kind=LONG), parameter :: modulus=2_LONG**31_LONG
    integer(kind=LONG), parameter :: factor=2_LONG**16_LONG
    integer(kind=LONG), parameter :: normalisation=32767_LONG
    integer :: icount

    temp_val=a*seed+c
    temp_val_prev=mod(temp_val,modulus)
    rand_array(1)=1.0_DP*(temp_val_prev/factor)/(1.0_DP*normalisation)

    do icount=2, num_rand
       temp_val=a*temp_val_prev+c
       temp_val_prev=mod(temp_val,modulus)
       rand_array(icount)=1.0_DP*(temp_val_prev/factor)/(1.0_DP*normalisation)
    enddo

  end subroutine utils_rand_gen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !==========================================================================!
  ! Aborts with a 'feature not supported' error message.                     !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   f (input): The name of the feature that is not there.                  !
  !--------------------------------------------------------------------------!
  ! Moved from multigrid_methods and generalized a bit                       !
  ! by Jacek Dziedzic, 2011/11/03.                                           !
  !--------------------------------------------------------------------------!
  subroutine utils_feature_not_supported(f)

    implicit none

    character(len=*), intent(in) :: f

    call utils_abort(trim(f)//' not present in this version.')

  end subroutine utils_feature_not_supported

!.............................................................................
!.............................................................................

    !============================================================================!
    ! The following three functions compare two input vectors, vecA and vecB,    !
    ! which must be of the same type (integer, real, or complex, respectively).  !
    ! The possible results are:                                                  !
    !    0   both vectors are equal                                              !
    !    1   vectors have different dimensions                                   !
    !   -1   vectors have the same dimensions, but different values              !
    ! The tolerance of the comparison can be set with the optional argument eps, !
    ! except for integer values, where the hardcoded tolerance is 0.             !
    !----------------------------------------------------------------------------!
    ! A possible use of the integer function is to check whether the shapes of   !
    ! 2 arrays conform: shapes of array1 and array2 conform if and only if       !
    !    utils_compare_vectors(shape(array1), shape(array2)) == 0                !
    !----------------------------------------------------------------------------!
    ! Written by Jose M Escartin in November 2014.                               !
    !============================================================================!

!.............................................................................

  integer function utils_compare_vectors_integer(vecA, vecB) result(return_code)

  implicit none

  ! Arguments
  integer, intent(in) :: vecA(:), vecB(:)

  if ( size(vecA) /= size(vecB) ) then
     ! vectors have different dimensions
     return_code = 1
  else
     if ( all (vecA == vecB) ) then
        return_code = 0
     else
        return_code = -1
     end if
  end if

  end function utils_compare_vectors_integer

!.............................................................................

  integer function utils_compare_vectors_real(vecA, vecB, eps) result(return_code)

  implicit none

  ! Arguments
  real(DP), intent(in) :: vecA(:), vecB(:)
  real(DP), intent(in), optional :: eps

  ! Local variables
  real(DP), parameter :: zero = 0.0_DP
  real(DP) :: loc_eps
  loc_eps = zero

  if ( present(eps) ) then
     if (eps < zero) then
        call utils_abort ('Error in utils_compare_vectors_real: negative eps.')
     else
        loc_eps = eps
     end if
  end if

  if ( size(vecA) /= size(vecB) ) then
     ! vectors have different dimensions
     return_code = 1
  else
     if ( loc_eps == zero ) then
        if ( all( vecA == vecB ) ) then
           return_code = 0
        else
           return_code = -1
        end if
     else
        if ( all( abs(vecA-vecB) < eps ) ) then
           return_code = 0
        else
           return_code = -1
        end if
     end if
  end if

  end function utils_compare_vectors_real

!.............................................................................

  integer function utils_compare_vectors_complex(vecA, vecB, eps) result(return_code)

  implicit none

  ! Arguments
  complex(DP), intent(in) :: vecA(:), vecB(:)
  real(DP), intent(in), optional :: eps

  ! Local variables
  real(DP), parameter :: zero = 0.0_DP
  real(DP) :: loc_eps
  loc_eps = zero

  if ( present(eps) ) then
     if (eps < zero) then
        call utils_abort ('Error in utils_compare_vectors_complex: negative eps.')
     else
        loc_eps = eps
     end if
  end if

  if ( size(vecA) /= size(vecB) ) then
     ! vectors have different dimensions
     return_code = 1
  else
     if ( loc_eps == zero ) then
        if ( all( vecA == vecB ) ) then
           return_code = 0
        else
           return_code = -1
        end if
     else
        if ( all( abs(vecA-vecB) < eps ) ) then
           return_code = 0
        else
           return_code = -1
        end if
     end if
  end if

  end function utils_compare_vectors_complex

!.............................................................................

  subroutine utils_flushed_string_output(string_to_output, &
       dont_reduce_error_flag)
    !==========================================================================!
    ! Outputs the string on root, then calls utils_flush.
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   string_to_output (input): The string that is to be output.             !
    !   dont_reduce_error_flag (input, optional): If present, it is passed to  !
    !                                             utils_flush.                 !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   - Carriage is not advanced. Attach CRLF to your string for newline.    !
    !   - When calling on a subset of ranks, pass .true. for 2nd argument, or  !
    !     else utils_flush is going to deadlock.                               !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, February 2015.                                !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout

    implicit none

    ! jd: Arguments
    character(len=*), intent(in)  :: string_to_output
    logical, optional, intent(in) :: dont_reduce_error_flag

    !------------------------------------------------------------------------

    if(pub_on_root) then
       write(stdout,'(a)',advance='no') string_to_output
    end if
    call utils_flush(stdout,dont_reduce_error_flag)

  end subroutine utils_flushed_string_output

!.............................................................................

  function utils_int_to_str(i, format_string) result(result)
    !==========================================================================!
    ! Returns a string corresponding to an int.                                !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   i (in): The int to convert.                                            !
    !   format_string (in, opt): Optional formatting string. Defaults to 'i0'. !
    ! Return value:                                                            !
    !   The string.                                                            !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic somewhere in 2015.                             !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    integer, intent(in) :: i
    character(len=*), intent(in), optional :: format_string

    ! jd: Local variables
    character(len=32) :: result

    !------------------------------------------------------------------------

    if(present(format_string)) then
       write(result, '('//format_string//')') i
    else
       write(result, '(i0)') i
    end if

  end function utils_int_to_str

!.............................................................................

  function utils_long_int_to_str(i, format_string) result(result)
    !==========================================================================!
    ! Returns a string corresponding to a long int.                            !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   i (in): The long int to convert.                                       !
    !   format_string (in, opt): Optional formatting string. Defaults to 'i0'. !
    ! Return value:                                                            !
    !   The string.                                                            !
    !--------------------------------------------------------------------------!
    ! Cloned by Jacek Dziedzic from utils_int_to_str() in May 2020.            !
    !==========================================================================!

    use constants, only: LONG

    implicit none

    ! jd: Arguments
    integer(kind=LONG), intent(in) :: i
    character(len=*), intent(in), optional :: format_string

    ! jd: Local variables
    character(len=32) :: result

    !------------------------------------------------------------------------

    if(present(format_string)) then
       write(result, '('//format_string//')') i
    else
       write(result, '(i0)') i
    end if

  end function utils_long_int_to_str

!.............................................................................

  character(len=32) function utils_real_to_str(r, format_string)
    !==========================================================================!
    ! Returns a string corresponding to a real. Five digits after the decimal  !
    ! point are used by default. The real is formatted as 'G' and is           !
    ! left-justified, which, when suitably trimmed afterwards, automatically   !
    ! takes care of the field width. Supply 'format_string' to override.       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r (in): The real to convert.                                           !
    !   format_string (in, opt): Opt. formatting string. Defaults to 'g32.5'.  !
    ! Return value:                                                            !
    !   The string.                                                            !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2015.                             !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in) :: r
    character(len=*), intent(in), optional :: format_string

    ! jd: Local variables
    character(len=32) :: buffer

    !------------------------------------------------------------------------

    if(present(format_string)) then
       write(buffer, '('//format_string//')') r
    else
       write(buffer, '(g32.5)') r
    end if
    utils_real_to_str = adjustl(buffer)

  end function utils_real_to_str

!.............................................................................

  character(len=100) function utils_point_to_str(r, format_string)
    !==========================================================================!
    ! Returns a string corresponding to a POINT.                               !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r (in): The POINT to convert.                                          !
    !   format_string (in): Format string (e.g. 'f12.6') for each component.   !
    ! Return value:                                                            !
    !   The string.                                                            !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2015.                             !
    !==========================================================================!

    use geometry, only: POINT

    implicit none

    ! jd: Arguments
    type(POINT), intent(in) :: r
    character(len=*), intent(in) :: format_string

    ! jd: Local variables
    character(len=100) :: buffer

    !------------------------------------------------------------------------

    write(buffer, '('//format_string//',a,'//format_string//',a,'//&
         format_string//')') r%x,',',r%y,',',r%z
!@    utils_point_to_str = adjustl(buffer)
    utils_point_to_str = buffer

  end function utils_point_to_str

!.............................................................................

  integer function utils_str_to_int(s)
    !==========================================================================!
    ! Returns an corresponding to a string.
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   s (in): The string to convert.                                         !
    ! Return value:                                                            !
    !   The int.                                                               !
    ! If the conversion fails, utils_assert is called.                         !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2015.                                  !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    character(len=*), intent(in) :: s

    ! jd: Locals
    integer :: ierr

    !------------------------------------------------------------------------

    read(s,*,iostat=ierr) utils_str_to_int
    call utils_assert(ierr == 0, &
         'utils_str_to_int: Expected an integer. Got the following: "'//&
         trim(s)//'"')

  end function utils_str_to_int

!.............................................................................

  subroutine utils_strip_char_from_string(string,ch)
    !==========================================================================!
    ! Strips instances of character ch from string                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   string (inout): The string that is to be stripped.                     !
    !   ch (input): The character to strip.                                    !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in March 2015.                                  !
    !==========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(inout) :: string
    character, intent(in) :: ch

    ! Local Variables
    integer :: ic
    logical :: keep_going
    character(len=3000) :: ctemp

    !------------------------------------------------------------------------
    keep_going = .true.
    do while(keep_going)
       ic = index(string,ch)
       if (ic>0) then
          ctemp = string
          string(ic:)=ctemp(ic+1:)//' '
       else
          keep_going=.false.
       end if
    end do


  end subroutine utils_strip_char_from_string

!.............................................................................

  subroutine utils_lock_ping(lock_name, ping_value, called_from_utils_abort)
    !==========================================================================!
    ! Pings a lock, by writing an integer to a named POSIX pipe.               !
    ! The pipe is opened beforehand and closed afterwards. The idea here is    !
    ! that whoever (another process on the same node) is doing a blocking read !
    ! on the pipe will become unblocked.                                       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   lock_name (in): A string with the filename of the lock (named pipe).   !
    !   ping_value (in, opt): An optional override of the value, the default   !
    !                         being 1.                                         !
    !   called_from_utils_abort (in, opt): Please don't use this. It will      !
    !                                      become handy when we make an attempt!
    !                                      to ping all active locks on abort.  !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   - Only run this, and utils_lock_block, on root. The assumption is that !
    !     the 'other side' is a single process and we don't want race condi-   !
    !     tions. An assert checks against calling this subroutine on non-root  !
    !     ranks.                                                               !
    !   - If the pipe cannot be opened, written to, or closed abort is called. !
    !   - No one ensures that the filename actually refers to a named pipe.    !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2015.                                   !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: CRLF

    implicit none

    ! jd: Arguments
    character(len=*), intent(in)  :: lock_name
    integer, intent(in), optional :: ping_value
    logical, intent(in), optional :: called_from_utils_abort

    ! jd: Local variables
    integer :: fileunit
    integer :: ierr

    ! -------------------------------------------------------------------------

    if(.not. present(called_from_utils_abort)) then
       call utils_flushed_string_output('ONETEP is waiting for an external &
            &consumer ('//trim(lock_name)//').'//CRLF,.true.)
    end if

    if(pub_on_root) then
       fileunit = utils_unit()
       open(unit=fileunit, file=trim(lock_name), action="write", status="old", &
            err=100, iostat=ierr)
       if(present(ping_value)) then
          write(fileunit,'(i0)',err=200, iostat=ierr) ping_value
       else
          write(fileunit,'(i0)',err=200, iostat=ierr) 1
       end if
       close(fileunit, err=300, iostat=ierr)
    else
       call utils_assert(present(called_from_utils_abort), &
            'Race condition in utils_lock_ping().')
    end if

    return

100 if(.not. present(called_from_utils_abort)) then
       call utils_abort('utils_lock_ping: Could not open lock file: '//&
            trim(lock_name), ierr)
    end if
200 if(.not. present(called_from_utils_abort)) then
       call utils_abort('utils_lock_ping: Could not write to lock file: '//&
            trim(lock_name), ierr)
    end if
300 if(.not. present(called_from_utils_abort)) then
       call utils_abort('utils_lock_ping: Could not close lock file: '//&
            trim(lock_name), ierr)
    end if

  end subroutine utils_lock_ping

!.............................................................................

  subroutine utils_lock_block(lock_name)
    !==========================================================================!
    ! Blocks on a lock, by reading an integer from a named POSIX pipe.         !
    ! The pipe is opened beforehand and closed afterwards. The idea here is    !
    ! that by doing a blocking read on a pipe we can wait until someone        !
    ! (another process on the same node) signals us to continue by writing to  !
    ! the pipe, thus unblocking us.                                            !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   lock_name (in): A string with the filename of the lock (named pipe).   !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   - Only run this, and utils_lock_block, on root. The assumption is that !
    !     the 'other side' is a single process and we don't want race condi-   !
    !     tions. An assert checks against calling this subroutine on non-root  !
    !     ranks.                                                               !
    !   - If the pipe cannot be opened, read from, or closed abort is called.  !
    !   - No one ensures that the filename actually refers to a named pipe.    !
    !   - The value that is read is not returned, as I repeated tests did not  !
    !     uncover the right settings of Fortran flags and formats that would   !
    !     return the actual value written to the pipe in the other process.    !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2015.                                   !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: CRLF

    implicit none

    ! jd: Arguments
    character(len=*), intent(in) :: lock_name

    ! jd: Local variables
    integer :: fileunit
    integer :: dummy, ierr

    ! -------------------------------------------------------------------------

    call utils_flushed_string_output('ONETEP is waiting for an &
         &external producer ('//trim(lock_name)//').'//CRLF,.true.)

    if(pub_on_root) then
       fileunit = utils_unit()
       open(unit=fileunit, file=trim(lock_name), action="read", &
            access="stream", form="formatted", iostat=ierr, &
            err=100)
       read(fileunit,'(i1)', iostat=ierr) dummy
       close(fileunit, err=300, iostat=ierr)
    else
       call utils_abort('Race condition in utils_lock_block().')
    end if

    return

100  call utils_abort('utils_lock_block: Could not open lock file: '//&
          trim(lock_name),ierr)
200  call utils_abort('utils_lock_block: Could not read from lock file: '//&
          trim(lock_name),ierr)
300  call utils_abort('utils_lock_block: Could not close lock file: '//&
          trim(lock_name),ierr)

  end subroutine utils_lock_block

!.............................................................................

  character(len=MAX_WORD_LENGTH) function utils_nth_word_of(string, n)
    !==========================================================================!
    ! Returns n-th word of a string.                                           !
    ! Words are deemed to be delimited by whitespace, currently defined as     !
    ! spaces and tabs. Multiple whitespace characters are coalesced, leading   !
    ! and trailing whitespace is ignored.                                      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   string (in): The string to extract the word from.                      !
    !   n (in): The number of the word wanted (1-based).                       !
    ! Return value:                                                            !
    !   n-th word of the string. If the string has fewer words, '' is returned.!
    ! Limitations:                                                             !
    !   There is no restriction on the length of the string.                   !
    !   The returned words are strings MAX_WORD_LENGTH-characters long.        !
    !   Any word longer than that will be truncated when returned.             !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2015.                                  !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    character(len=*), intent(in) :: string
    integer, intent(in)          :: n

    ! jd: Locals
    character(len=2), parameter :: whitespace = char(9)//' '
    integer :: word, curpos, blankpos, leadblank, j
    integer :: maxlen

    ! -------------------------------------------------------------------------

    maxlen = len(string)

    ! Skip leading whitespace
    do leadblank = 1, maxlen
      if(scan(string(leadblank:leadblank),whitespace) == 0) exit
    end do
    curpos = leadblank

    ! Go over tokens (words)
    do word = 1, n
       ! Look for the next whitespace
       blankpos=scan(string(curpos:),whitespace)

       if(blankpos /= 0) then
          ! Found more whitespace
          if(word == n) then
             ! If we counted to n by now, return everything up to this
             ! whitespace, exclusive.
             utils_nth_word_of = string(curpos:curpos+blankpos-2)
             return
          end if
          ! Still looking for more words. Skip to the next non-whitespace,
          ! in case there are multiple consecutive whitespaces.
          do j = blankpos, maxlen
             if(curpos+j > maxlen) then
                ! Oops, ran out of string, not enough words
                utils_nth_word_of = ''
                return
             end if
             if(scan(string(curpos+j:curpos+j),whitespace) == 0) exit
          end do
          ! Position at the next non-whitespace for next iteration
          curpos = curpos + j
       else
          ! No more whitespace.
          if(word == n) then
             ! Return what's left over, if we counted to n.
             utils_nth_word_of = string(curpos:)
          else
             ! Or an empty string, if we didn't
             utils_nth_word_of = ''
          end if
          return
       end if
    end do

    ! We should never get here.
    call utils_abort('Internal error in utils_nth_word_of')
    utils_nth_word_of = 'cromulent' ! silence 'no return value' warning

  end function utils_nth_word_of

!.............................................................................

  function utils_banner(adornment, title, delimiter, opt_terminal_width) result(string)
     !=======================================================================!
     ! Prints a centred title with optional adornments and delimiters        !
     ! The structure of the banner is:                                       !
     ! delimiter (repeated adornment) title (repeated adornment) delimiter   !
     ! Spaces may be added around the title to improve readability.          !
     !-----------------------------------------------------------------------!
     ! Written by J.M. Escartin in April 2015.                               !
     ! Trivial extension to shorter banners (via opt_terminal_width)         !
     ! by Jacek Dziedzic, February 2019.                                     !
     !=======================================================================!

!!!  Temporary fixes due to PGI bug.
     implicit none

     ! Output
     character(len=80) :: string

     ! Arguments
!!!     character(len=*), intent(in), optional :: adornment, title, delimiter
     character(len=*), intent(in), optional :: title, delimiter  !!! temporary fix
     character(len=*), intent(in) :: adornment                   !!! temporary fix
     integer, intent(in), optional :: opt_terminal_width ! jd: for shorter banners

     ! Local Variables
     integer, parameter :: terminal_width = 80
     integer :: loc_terminal_width
!!!     character, parameter :: default_adornment = ' '
     character(len=:), allocatable :: loc_adornment
!!!     character(len=:), allocatable :: reversed_delimiter
     character(len=len(adornment)) :: reversed_adornment        !!! temporary fix
     character(len=:), allocatable :: reversed_delimiter
     logical :: title_spaces
     integer :: len_title, len_adornment, len_delimiter
     integer :: len_both_adornments, len_left_adornments, len_right_adornments
     integer :: num_left_adornments, num_right_adornments, frac_adornment
     integer :: min_string_len
     integer :: i0, i1

     ! Check present variables and get their lengths.
     if (present(title)) then
        len_title = len(title)
     else
        len_title = 0
     end if
     if (present(delimiter)) then
        len_delimiter = len(delimiter)
        ! jd: FLANG v12.0.0 ICE's on this line.
        !     Workaround: delimiters will not be reversed with FLANG.
#ifndef __FLANG
        reversed_delimiter = internal_reverse_string(delimiter)
#else
        reversed_delimiter = delimiter
#endif
     else
        len_delimiter = 0
     end if
!!!     if (present(adornment)) then
        loc_adornment = adornment
!!!     else
!!!        loc_adornment = default_adornment
!!!     end if
     len_adornment = len(loc_adornment)
        ! jd: FLANG v12.0.0 ICE's on this line.
        !     Workaround: adornments will not be reversed with FLANG.
#ifndef __FLANG
     reversed_adornment = internal_reverse_string(loc_adornment)
#else
     reversed_adornment = loc_adornment
#endif

     if(present(opt_terminal_width)) then
        loc_terminal_width = opt_terminal_width
     else
        loc_terminal_width = terminal_width
     end if

     ! Check if banner fits in terminal.
     min_string_len = len_title + 2 * ( len_delimiter + len_adornment )
!!!     if ( present(title) .and. (present(adornment).or.present(delimiter)) ) then
     if ( present(title) ) then                              !!! temporary fix
        ! Spaces will be needed.
        title_spaces = .True.
        min_string_len = min_string_len + 2
     else
        title_spaces = .False.
     end if
     if (min_string_len > loc_terminal_width) then
        call utils_abort ('Error in utils_banner: &
             &banner length exceeds default terminal width.')
     end if

     ! Set the number of characters to be filled with adornments.
     len_both_adornments = loc_terminal_width - min_string_len + 2 * len_adornment
     ! Split them into left and right.
     num_left_adornments = len_both_adornments / (2 * len_adornment)
     len_left_adornments = len_adornment * num_left_adornments
     len_right_adornments = len_both_adornments - len_left_adornments
     num_right_adornments = len_right_adornments / len_adornment
     ! Characters left.
     frac_adornment = len_both_adornments - &
          (len_right_adornments + len_left_adornments)

     ! Build the string
     i0 = 1
     i1 = 0
     if (present(delimiter)) then
        i1 = len_delimiter
        string(i0:i1) = delimiter
        i0 = i1 + 1
     end if
     i1 = i1 + num_left_adornments * len_adornment
!!!     string(i0:i1) = repeat(loc_adornment, num_left_adornments)
     string(i0:i1) = repeat(adornment, num_left_adornments)  !!! temporary fix
     i0 = i1 + 1
     if (present(title)) then
        if (title_spaces) then
           i1 = i1 + 1
           string(i0:i1) = ' '
           i0 = i1 + 1
        end if
        i1 = i1 + len_title
        string(i0:i1) = title
        i0 = i1 + 1
        if (title_spaces) then
           i1 = i1 + 1
           string(i0:i1) = ' '
           i0 = i1 + 1
        end if
     end if
     i1 = i1 + num_right_adornments * len_adornment
     string(i0:i1) = repeat(reversed_adornment, num_right_adornments)
     i0 = i1 + 1
     if (frac_adornment > 0) then
        i1 = i1 + frac_adornment
        string(i0:i1) = reversed_adornment(1:frac_adornment)
        i0 = i1 + 1
     end if
     if (present(delimiter)) then
        i1 = i1 + len_delimiter
        string(i0:i1) = reversed_delimiter
     end if

     if (i1 /= loc_terminal_width) then
        call utils_abort ("Error in utils_banner: &
             &final banner length doesn't match default terminal width." // string, i1)
     end if

     ! jd: Ensure we pad with spaces to 80 chars
     if(present(opt_terminal_width)) then
        string(i1+1:) = ' '
     end if

  contains

     function internal_reverse_string(my_string) result(reverse_string)
        implicit none

        character(len=*), intent(in) :: my_string
        character(len=len(my_string)) :: reverse_string

        ! Local variables
        integer, parameter :: nchar = 4
        logical :: mask(len(my_string),nchar)
        character(len=1), dimension(nchar) :: ch, rch
        integer :: length
        integer :: i, j

        ! Reverse string keeping characters.
        length = len(my_string)
        do j = 1, length
           i = length - j + 1
           reverse_string(j:j) = my_string(i:i)
        end do

        ! Exchange '<' and '>', '/' and '\', etc.
        ! jd: Different compilers have different ideas about
        !     escaping '\'.
#ifdef __FLANG
        ch  = [ '<' , '>' , '/', '\\' ]
        rch = [ '>' , '<' , '\\', '/' ]
#else
        ch  = [ '<' , '>' , '/', '\' ]
        rch = [ '>' , '<' , '\', '/' ]
#endif
        do i = 1, nchar
           mask(:,i) = ( reverse_string == ch(i) )
        end do
        do i = 1, nchar
            do j = 1, length
               if ( mask(j,i) ) reverse_string(j:j) = rch(i)
            end do
        end do

     end function internal_reverse_string

  end function utils_banner

!.............................................................................

  subroutine utils_erase_file(filename)
    !=======================================================================!
    ! Erases a file with a given filename.                                  !
    ! If the file did not exist, or could not be opened for any other       !
    ! reason -- nothing happens.                                            !
    ! If the file could not be deleted -- nothing happens.                  !
    !-----------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2016.                            !
    !=======================================================================!

    implicit none

    ! jd: Arguments
    character(len=*), intent(in) :: filename

    ! jd: Local variables
    integer :: erase_unit

    ! -------------------------------------------------------------------------

    erase_unit = utils_unit()
    open(unit=erase_unit, file=filename, status="old", err=500)
    close(unit=erase_unit, status="delete", err=500)
500 continue ! jd: Ignore errors

  end subroutine utils_erase_file

!.............................................................................

  real(kind=DP) elemental function utils_logical_to_real(l)
    !==========================================================================!
    ! Returns 0D0 if the argument is .true., and something else if .false..    !
    ! Trivial, but it makes packing arrays of logicals in arrays of reals      !
    ! bulletproof.                                                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   l (in): The logical to convert.                                        !
    ! Return value:                                                            !
    !   The real.                                                              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in August 2016.                                !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    logical, intent(in) :: l

    !------------------------------------------------------------------------

    utils_logical_to_real = merge(1D0,0D0,l)

  end function utils_logical_to_real

!.............................................................................

  logical elemental function utils_real_to_logical(r)
    !==========================================================================!
    ! Returns .true. if r is anything but exactly 0D0.                         !
    ! Trivial, but it makes unpacking arrays of logicals from arrays of reals  !
    ! bulletproof.                                                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r (in): The real to convert.                                           !
    ! Return value:                                                            !
    !   The logical.                                                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in August 2016.                                !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in) :: r

    !------------------------------------------------------------------------

    utils_real_to_logical = (r /= 0D0)

  end function utils_real_to_logical

!.............................................................................

  real(kind=DP) elemental function utils_character_to_real(c)
    !==========================================================================!
    ! Returns some kind of opaque representation of a character as a real.     !
    ! Trivial, but it makes packing arrays of characters in arrays of reals    !
    ! bulletproof.                                                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   c (in): The character to convert.                                      !
    ! Return value:                                                            !
    !   The real.                                                              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2019.                                 !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    character, intent(in) :: c

    !------------------------------------------------------------------------

    utils_character_to_real = real(ichar(c),kind=DP)

  end function utils_character_to_real

!.............................................................................

  character elemental function utils_real_to_character(r)
    !==========================================================================!
    ! Returns a character from a real representation previously obtained by    !
    ! utils_character_to_real(). No guarantees are made on the internals of    !
    ! the representation (opaque).                                             !
    ! Trivial, but it makes unpacking arrays of logicals from arrays of reals  !
    ! bulletproof.                                                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r (in): The real to convert.                                           !
    ! Return value:                                                            !
    !   The character.                                                         !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2019.                                 !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in) :: r

    !------------------------------------------------------------------------

    utils_real_to_character = char(int(r))

  end function utils_real_to_character

!.............................................................................

  subroutine utils_filename_length_check(filename,routine,max_length,file_output)
    !==========================================================================!
    ! Check the length of a filename and return an error if it is too long     !
    ! If no maximum length is specified, the routine will check against the    !
    ! maximum length specified in the constants module                         !
    !==========================================================================!
    ! filename (in)   : filename, ideally trimmed to remove trailing blanks    !
    ! routine (in)    : name of calling routine, for error reporting           !
    ! max_length (in) : optional maximum length, if global value from          !
    !                   constants module is not used                           !
    ! file_output (in, optional) : if present, and .false., do not output      !
    !                              error message to file (passed to            !
    !                              utils_abort, necessary if utils module's    !
    !                              local pub_rootname has not been             !
    !                              initialized, default is .true.)             !
    !==========================================================================!
    ! James C. Womack, 09/2016
    !==========================================================================!
    use constants, only: file_maxsize
    implicit none

    ! Arguments
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: routine
    integer, optional, intent(in) :: max_length
    logical, optional, intent(in) :: file_output

    ! Parameters
    character(len=*), parameter :: myself = "utils_filename_length_check"

    ! Local variables
    integer :: loc_max_length
    character(len=6) :: int_str ! <-- should be enough to hold any integer
                                ! character string length anyone could want!
    logical :: loc_file_output

    if ( present(max_length) ) then
       loc_max_length = max_length
    else
       loc_max_length = file_maxsize ! <-- from constants module
    end if

    loc_file_output = .true.
    if ( present(file_output) ) loc_file_output = file_output

    if (len(trim(filename)) <= loc_max_length) then
      ! JCW: No problem
      return
    else
      ! JCW: Overlong filename, abort
      write(int_str,'(i6)') loc_max_length
      call utils_abort("Error in "//myself//": &
           &Filename "//trim(filename)//" used in routine "//trim(routine)//" &
           &is longer than "//trim(adjustl(int_str))//" characters. Please &
           &shorten the filename (including any directory path included) and &
           &run ONETEP again.",file_output=loc_file_output)
    end if

  end subroutine utils_filename_length_check

!.............................................................................

! jme: This subroutine was written in order to copy arbitrary files, but its use
!      introduced a performance regression (see CCPForge bug #1702), and since
!      Nov 2017, when the fix was committed (v4.5.14.10, r1645), it is no longer
!      called from anywhere in the code.
!
! subroutine utils_stream_copy(original_file,target_file)
!   !==========================================================================!
!   ! Makes an exact copy of a file, using unformatted stream I/O              !
!   ! (Fortran 2003). Works on formatted or unformatted files. Original and    !
!   ! target files must not already be open. If the target file already exists !
!   ! it is replaced.                                                          !
!   !--------------------------------------------------------------------------!
!   ! Arguments:                                                               !
!   !   original_file (in): Filename of the file to copy.                      !
!   !   target_file (in): Filename of the target file.                         !
!   !--------------------------------------------------------------------------!
!   ! Written by Fabiano Corsetti in September 2016.                           !
!   !==========================================================================!

!   use constants, only: LONG
!   implicit none

!   character(*), intent(in) :: original_file
!   character(*), intent(in) :: target_file

!   integer, parameter :: one_byte_int_kind = selected_int_kind(2)

!   integer :: i, ierr, ipos, extra, original_unit, target_unit
!   integer(kind=LONG) :: i8
!   integer(kind=one_byte_int_kind) :: i1

!   ! jme: check that i1 has 8 bits
!   call utils_assert(storage_size(i1)==8, 'Error in utils_stream_copy: &
!        &the processor didn''t provide a 1 Byte integer.')

!   original_unit=utils_unit()
!   open(original_unit,file=trim(original_file), &
!        form='unformatted',access='stream',status='old',iostat=ierr)
!   call utils_open_unit_check('stream_copy','original_file',ierr)

!   target_unit=utils_unit()
!   open(target_unit,file=trim(target_file), &
!        form='unformatted',access='stream',status='replace',iostat=ierr)
!   call utils_open_unit_check('stream_copy','target_file',ierr)

!   do
!      read(original_unit,iostat=ierr) i8
!      if (ierr/=0) then
!         inquire(original_unit,pos=ipos)
!         exit
!      end if
!      write(target_unit) i8
!   end do
!   extra=mod(ipos-1,8)
!   ipos=ipos-extra

!   do i=0,extra-1
!   jme: this prevents CRAY bug #848598 (fixed in CCE 8.5.8).
!#ifdef _CRAYFTN
!!dir$ suppress i
!#endif
!      read(original_unit,pos=ipos+i) i1
!      write(target_unit) i1
!   end do

!   close(target_unit,iostat=ierr)
!   call utils_close_unit_check('stream_copy','target_file',ierr)

!   close(original_unit,iostat=ierr)
!   call utils_close_unit_check('stream_copy','original_file',ierr)

! end subroutine utils_stream_copy

!.............................................................................

  function utils_safe_int(x) result(n)
     !=============================================!
     ! Take integer part, with safety.             !
     !---------------------------------------------!
     ! Written by JM Escartin on 4th October 2016. !
     !=============================================!
     implicit none

     ! Arguments and result
     real(kind=DP), intent(in) :: x
     integer :: n

     ! Local variables
     character(len=25) :: str

     n = int(x) ! Always defined

     ! If x was too big, abort.
     if (abs(x) > real(huge(n), kind=DP)) then
        write(str, '(1P,G25.15)') x
        if (x > 0.0_DP) then
           call utils_abort('Error in utils_safe_int: input real value (' // &
                trim(adjustl(str)) // ') is bigger than the largest positive &
                &integer of the current default kind', huge(n))
        else ! x < 0
           call utils_abort('Error in utils_safe_int: input real value (' // &
                trim(adjustl(str)) // ') is bigger in magnitude than the &
                &largest positive integer of the current default kind', huge(n))
        end if
     end if

  end function utils_safe_int

!-------------------------------------------------------------------------------

  function utils_safe_nint(x) result(n)
     !=============================================!
     ! Compute nearest integer, with safety.       !
     !---------------------------------------------!
     ! Written by JM Escartin on 4th October 2016. !
     !=============================================!
     implicit none

     ! Arguments and result
     real(kind=DP), intent(in) :: x
     integer :: n

     ! Local variables
     character(len=25) :: str

     n = nint(x) ! Always defined

     ! If x was too big, abort.
     if (abs(x) > real(huge(n), kind=DP)) then
        write(str, '(1P,G25.15)') x
        if (x > 0.0_DP) then
           call utils_abort('Error in utils_safe_nint: input real value (' // &
                trim(adjustl(str)) // ') is bigger than the largest positive &
                &integer of the current default kind', huge(n))
        else ! x < 0
           call utils_abort('Error in utils_safe_nint: input real value (' // &
                trim(adjustl(str)) // ') is bigger in magnitude than the &
                &largest positive integer of the current default kind', huge(n))
        end if
     end if

  end function utils_safe_nint

  subroutine utils_parse_bc(bc_string,bc_is_zero,bc_is_periodic)
    !=========================================================================!
    ! Parses a character string containing BC designators, as formatted       !
    ! in the pub_multigrid_bc parameter. This is a string which contains 3    !
    ! characters indicating the BC along each lattice vector direction        !
    !     "O", "o" : open BC (open with Dirichlet BCs computed at cell faces) !
    !     "P", "p" : periodic BC                                              !
    !     "Z", "z" : zero BC (open but potential is zero at cell faces)       !
    ! We assume that bc_string has already passed through some checks in      !
    ! rundat, so that it only contains the allowed characters (" OoPpZz")     !
    ! and has exactly three of the non-whitespace characters.                 !
    !                                                                         !
    ! The result of the parsing can be returned in bc_is_zero, which is       !
    ! .true. for zero BCs and .false. otherwise for each direction.           !
    !                                                                         !
    ! A logical array, bc_is_periodic, can also be output, indicating whether !
    ! the grid is periodic (T) or non-periodic (F) along each lattice vector. !
    !                                                                         !
    ! A separate bc_is_zero array is useful since DL_MG does not              !
    ! distinguish between computed Dirichlet BCs and zero BCs, but ONETEP     !
    ! needs to use this to determine whether to compute the BCs.              !
    !                                                                         !
    ! The separate bc_is_periodic array can be directly passed to             !
    ! comms_cart_create, and used to create the MPI communicator with         !
    ! Cartesian topology required by DL_MG.                                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  bc_string           intent=in,  the string of BC designators           !
    !  bc_is_zero     intent=out, opt, array of logicals indicating whether   !
    !                                  zero BCs are used in each direction    !
    !  bc_is_periodic intent=out, opt, array of logicals indicating whether   !
    !                                  BCs are periodic in each direction     !
    !-------------------------------------------------------------------------!
    ! Written by James C. Womack, 05/2017                                     !
    !=========================================================================!

    implicit none

    ! Arguments
    character(len=*),      intent(in)  :: bc_string
    logical, dimension(3), optional, intent(out) :: bc_is_zero
    logical, dimension(3), optional, intent(out) :: bc_is_periodic

    ! Parameters
    character(len=*), parameter :: myself = "utils_parse_bc"

    ! Local variables
    logical :: bc_is_zero_loc(3)
    logical :: bc_is_periodic_loc(3)
    integer :: ibc
    integer :: ichar

    ibc = 1

    do ichar = 1, len(bc_string)

       select case ( bc_string(ichar:ichar) )

       case("O","o")
         bc_is_zero_loc(ibc)     = .false.
         bc_is_periodic_loc(ibc) = .false.
         ibc = ibc + 1
       case("P","p")
         bc_is_zero_loc(ibc)     = .false.
         bc_is_periodic_loc(ibc) = .true.
         ibc = ibc + 1
       case("Z","z")
         bc_is_zero_loc(ibc)     = .true.
         bc_is_periodic_loc(ibc) = .false.
         ibc = ibc + 1
       case(" ")
         ! Do nothing
       case default
         call utils_abort("Error in "//myself//": Unrecognized BC designator '"//&
              bc_string(ichar:ichar)//"'.")
       end select

       ! Only 3 BC designators allowed
       if (ibc > 3) exit

    end do

    if (ibc-1 < 3) then
         call utils_abort("Error in "//myself//": Too few BC designators in &
              &bc_string = '"//trim(bc_string)//"'.")
    end if

    ! Copy local variables to optional arguments, if present
    if (present(bc_is_zero)) bc_is_zero(1:3) = bc_is_zero_loc(1:3)
    if (present(bc_is_periodic)) bc_is_periodic(1:3) = bc_is_periodic_loc(1:3)

  end subroutine utils_parse_bc

!.............................................................................

  character(len=MAX_WORD_LENGTH) function utils_postfix_to_ngwf_set_name(rep_postfix)
    !==========================================================================!
    ! Parses rep_postfix, returns a corresponding human-readable NGWF set name.!
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in June 2018.                                  !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    character(len=*), intent(in) :: rep_postfix

    ! jd: Local variables
    character(len=*), parameter :: myself = 'utils_postfix_to_ngwf_set_name'

    ! -----------------------------------------------------------------------

    if (trim(rep_postfix)=='c') then
       utils_postfix_to_ngwf_set_name = 'cond'
    else if (trim(rep_postfix)=='j') then
       utils_postfix_to_ngwf_set_name = 'joint'
    else if (trim(rep_postfix)=='a') then
       utils_postfix_to_ngwf_set_name = 'aux'
    else if (trim(rep_postfix)=='z') then
       utils_postfix_to_ngwf_set_name = 'jointPH'
    else if (trim(rep_postfix)=='#') then ! swex not attached to any NGWF_REP
       utils_postfix_to_ngwf_set_name = 'free-floating'
    else if (trim(rep_postfix)=='') then
       utils_postfix_to_ngwf_set_name = 'val'
    else
       call utils_abort(trim(myself)//': Unrecognized rep_postfix: "'//&
            trim(rep_postfix)//'".')
    end if

  end function utils_postfix_to_ngwf_set_name

!.............................................................................

  subroutine utils_update_rootname(rootfilename, delete_err_file, file_output)
     !=========================================================================!
     ! # DESCRIPTION                                                           !
     ! As utils cannot be dependent on rundat.mod, we must have a local        !
     ! version of pub_rootname ("pub_rootame_local"). This routine updates     !
     ! this local version. There are very few instances where this variable    !
     ! should be modified, so this is done via a routine (rather than simply   !
     ! making pub_rootname_local public)                                       !
     !-------------------------------------------------------------------------!
     ! # ARGUMENTS                                                             !
     ! rootfilename       in    the new value for pub_rootname_local           !
     ! delete_err_file    in    at the same time, delete any pre-existing      !
     !                          .error_message file of that name               !
     ! file_output        in    allow output to file                           !
     !-------------------------------------------------------------------------!
     ! # AUTHORS & CHANGELOG                                                   !
     ! Author(s):        Edward Linscott                                       !
     ! Date of creation: May 2019                                              !
     ! Split from code previously in utils_abort                               !
     !=========================================================================!

     use comms, only: pub_on_root

     implicit none

     character(len=80), intent(in) :: rootfilename
     logical, optional, intent(in) :: delete_err_file
     logical, optional, intent(in) :: file_output

     logical             :: loc_delete_err_file
     logical             :: loc_file_output
     character(len=1024) :: error_message_file

     ! By default, delete pre-existing .error_message file
     loc_delete_err_file = .true.
     if (present(delete_err_file)) loc_delete_err_file = delete_err_file

     ! By default, permit output to file
     loc_file_output = .true.
     if (present(file_output)) loc_file_output = file_output

     ! JCW: Only store pub_rootname_local if file output is allowed
     if (loc_file_output) then
        pub_rootname_local = rootfilename
        ! Delete pub_rootname_local.error_message
        if (pub_on_root .and. delete_err_file) then
           write(error_message_file,'(a,a)') trim(pub_rootname_local), &
                 '.error_message'
           call utils_erase_file(error_message_file)
        end if
     end if

  end subroutine utils_update_rootname
!-------------------------------------------------------------------------------

end module utils
