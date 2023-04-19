! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                               alloc_check module                            !
!=============================================================================!
! This module is far from being completed yet. Treat it as a stub for now.    !
!@docme                                                                       !
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
!=============================================================================!

module alloc_checker

  use constants, only: CRLF, DP, stdout

  implicit none

  private

  public :: alloc_init
  public :: alloc_fail
  public :: alloc_in
  public :: alloc_out
  public :: alloc_report

  public :: ah_, as_, ai_, ax_, ay_

  integer, save :: ah_, as_, ai_

  integer, save :: aa_ = 0 ! n alloc calls
  integer, save :: ad_ = 0 ! n dealloc calls
  integer, save :: ap_ = 0 ! alloc/dealloc pool (bytes)
  integer, save :: aw_ = 0 ! total alloc watermark (max ap_ ever)

  integer, parameter :: max_id_len = 64
  integer, parameter :: max_n_ids = 256
  integer, parameter :: max_n_threads = 128
  integer :: mem_n_ids(max_n_threads) = 0
  character(len=max_id_len) :: mem_routine_names(max_n_ids, max_n_threads)
  character(len=max_id_len) :: mem_var_names(max_n_ids, max_n_threads)
  integer                   :: mem_bytes(max_n_ids, max_n_threads)
  integer                   :: mem_call_balance(max_n_ids, max_n_threads)

  character(len=max_id_len), save :: ax_, ay_

  real(kind=DP), parameter :: mibibyte = 1048576.0_DP

!$OMP THREADPRIVATE(aa_, ad_, ah_, ai_, as_, ap_, aw_, ax_, ay_)

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_init
    !==========================================================================!
    !--------------------------------------------------------------------------!
    !--------------------------------------------------------------------------!
    !==========================================================================!

     implicit none

     mem_n_ids(:) = 0
     mem_routine_names(:,:) = ''
     mem_var_names(:,:) = ''
     mem_bytes(:,:) = 0
     mem_call_balance(:,:) = 0

  end subroutine alloc_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_fail
    !==========================================================================!
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   routine (in) : Name of the subroutine which called this routine        !
    !   array   (in) : Name of the array allocated                             !
    !   ierr    (in) : Flag returned by allocate                               !
    !--------------------------------------------------------------------------!
    !==========================================================================!

    use utils, only: utils_abort, utils_int_to_str

    implicit none

    ! jd: Local variables

    ! -------------------------------------------------------------------------

    call utils_abort('Error '//trim(utils_int_to_str(ai_))//&
         ' allocating array "'//trim(ay_)//"'"//CRLF//&
         ' in subroutine or function "'//trim(ax_)//'".')

  end subroutine alloc_fail

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_in(footprint)
    !==========================================================================!
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   routine (in) : Name of the subroutine which called this routine        !
    !   array   (in) : Name of the array deallocated                           !
    ! @
    !--------------------------------------------------------------------------!
    !==========================================================================!

    use utils, only: utils_abort, utils_alloc_check

    implicit none

    ! jd: Arguments
    integer, intent(in)          :: footprint

    ! jd: Local variables

    ! -------------------------------------------------------------------------

!$OMP CRITICAL(SEC_ALLOC_CHECKER_GLOBAL_STATE)
!    write(stdout,*) trim(ax_)//'@:@'//trim(ay_)//': +', footprint

    ! Check the allocation succeeded
    call utils_alloc_check(trim(ax_), trim(ay_), 0)

    ! Remember that alloc size has not been stored (yet)
    ah_ = 0
    as_ = -1

    ! Add alloc to call balance
    aa_ = aa_ + 1

    ! Add alloc footprint to current balance
    ap_ = ap_ + footprint

    ! Update high-water mark, if exceeded
    aw_ = max(aw_,ap_)

    ! Check for leaks
#ifdef ALLOC_LEAK_CHECK
    call alloc_leak_check_in(ax_,ay_,footprint)
#endif

!$OMP END CRITICAL(SEC_ALLOC_CHECKER_GLOBAL_STATE)

  end subroutine alloc_in

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_out
    !==========================================================================!
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   routine (in) : Name of the subroutine which called this routine        !
    !   array   (in) : Name of the array deallocated                           !
    ! @
    !--------------------------------------------------------------------------!
    !==========================================================================!

    use utils, only: utils_abort, utils_dealloc_check

    implicit none

    ! jd: Local variables

    ! -------------------------------------------------------------------------

!$OMP CRITICAL(SEC_ALLOC_CHECKER_GLOBAL_STATE)
!    write(stdout,*) trim(ax_)//'@:@'//trim(ay_)//': -', as_

    ! If allocation status has been verified, judge if it worked
    if(ah_ == 1) then
       if(as_ == -1) then
          call utils_abort('Attempt to deallocate an array "'//trim(ay_)//&
               '" whose allocation status is already'//CRLF//'.false. in &
               &subroutine or function "'//trim(ax_)//'".')
       end if
    end if

    if(ai_ /= 0) then
       call utils_abort('Error in '//trim(ax_)//': deallocating '//&
            trim(ay_)//' failed with code ',ai_)
    end if

    call utils_dealloc_check(trim(ax_), trim(ay_), 0)

    ! Subtract alloc footprint from current balance
    ap_ = ap_ - as_

    ! Remember that alloc size has not been stored (yet)
    as_ = -1
    ah_ = 0

    ! Add dealloc to call balance
    ad_ = ad_ + 1

    ! Check for leaks
#ifdef ALLOC_LEAK_CHECK
    call alloc_leak_check_out(ax_,ay_,as_)
#endif

!$OMP END CRITICAL(SEC_ALLOC_CHECKER_GLOBAL_STATE)

  end subroutine alloc_out

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef ALLOC_COUNT
  subroutine alloc_report
    !==========================================================================!
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   routine (in) : Name of the subroutine which called this routine        !
    !   array   (in) : Name of the array deallocated                           !
    ! @
    !--------------------------------------------------------------------------!
    !==========================================================================!

    use comms, only: pub_on_root, pub_my_proc_id

    implicit none

    ! jd: Arguments

    ! jd: Local variables

    ! -------------------------------------------------------------------------

    if(pub_on_root) then
       write(stdout,'(a,i6,a)') CRLF//&
            '+------------------ MEMORY USE REPORT (proc ',pub_my_proc_id, &
            ') -------------------+'
       write(stdout,'(a,i10,a)') '| Total number of allocations:           ', &
            aa_, '                    |'
       write(stdout,'(a,i10,a)') '| Total number of deallocations:         ', &
            ad_, '                    |'
       write(stdout,'(a,i10,a)') '| Balance:                               ', &
            aa_-ad_, '                    |'
       write(stdout,'(a,f13.3,a)') '| Total leaked memory (MiB):          ', &
            ap_/mibibyte, '                    |'
       write(stdout,'(a,f13.3,a)') '| High-water mark (MiB):              ', &
            aw_/mibibyte, '                    |'
       write(stdout,'(a)') &
            '+-------------------------------------------------------------&
            &---------+'

    end if

  end subroutine alloc_report
#else

  subroutine alloc_report

    use comms, only: pub_on_root

    implicit none

    ! -------------------------------------------------------------------------

    if(pub_on_root) then
#if 0
       write(stdout,'(a)') CRLF//&
            'Compile ONETEP with -DALLOC_COUNT to see a memory use report.'
#endif
    end if

  end subroutine alloc_report

#endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_leak_check_in(routine_name, var_name, footprint)
    !==========================================================================!
    !--------------------------------------------------------------------------!
    !--------------------------------------------------------------------------!
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use utils, only: utils_int_to_str, utils_assert

    implicit none

    ! jd: Arguments
    character(len=max_id_len), intent(in) :: routine_name
    character(len=max_id_len), intent(in) :: var_name
    integer, intent(in) :: footprint

    ! jd: Local variables
    character(len=max_id_len) :: stored_var_name, stored_routine_name
    integer :: thread_id
    integer :: i_name
    logical :: found_it
!$  integer, external :: omp_get_thread_num

    ! -------------------------------------------------------------------------

    thread_id = 1
!$  thread_id = omp_get_thread_num() + 1

    ! See if we have this name in the list
    found_it = .false.
    do i_name = 1, mem_n_ids(thread_id)
       stored_routine_name = mem_routine_names(i_name, thread_id)
       stored_var_name = mem_var_names(i_name, thread_id)

       if(fuzzy_string_match(stored_var_name, var_name) .and. &
          fuzzy_string_match(stored_routine_name, routine_name)) then

          ! Already have it, increment counters
          mem_bytes(i_name, thread_id) = &
               mem_bytes(i_name, thread_id) + footprint
          mem_call_balance(i_name, thread_id) = &
               mem_call_balance(i_name, thread_id) + 1
          found_it = .true.
          exit
       end if

    end do

    if(found_it) return

    ! Don't have this name in the list, add it
    mem_n_ids(thread_id) = mem_n_ids(thread_id) + 1
    call utils_assert(mem_n_ids(thread_id) <= max_n_ids, &
         'alloc_leak_check_in: Maximum number of tracked allocatable &
         &entities ('//trim(utils_int_to_str(max_n_ids))//&
         ') exceeded -- thread '//trim(utils_int_to_str(thread_id))//&
         ', rank '//trim(utils_int_to_str(pub_my_proc_id)))

    mem_routine_names(mem_n_ids(thread_id), thread_id) = routine_name
    mem_var_names(mem_n_ids(thread_id), thread_id) = var_name
    mem_bytes(mem_n_ids(thread_id), thread_id) = footprint
    mem_call_balance(mem_n_ids(thread_id), thread_id) = 1

  end subroutine alloc_leak_check_in

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_leak_check_out(routine_name, var_name, footprint)
    !==========================================================================!
    !--------------------------------------------------------------------------!
    !--------------------------------------------------------------------------!
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use utils, only: utils_assert, utils_int_to_str

    implicit none

    ! jd: Arguments
    character(len=max_id_len), intent(in) :: routine_name
    character(len=max_id_len), intent(in) :: var_name
    integer, intent(in) :: footprint

    ! jd: Local variables
    character(len=max_id_len) :: stored_var_name, stored_routine_name
    integer :: thread_id
    integer :: i_name
    integer :: j_name
    logical :: found_it
!$  integer, external :: omp_get_thread_num

    ! -------------------------------------------------------------------------

    thread_id = 1
!$  thread_id = omp_get_thread_num() + 1

    ! See if we have this name in the list
    found_it = .false.
    do i_name = 1, mem_n_ids(thread_id)
       stored_routine_name = mem_routine_names(i_name, thread_id)
       stored_var_name = mem_var_names(i_name, thread_id)

       if(fuzzy_string_match(stored_var_name, var_name) .and. &
          fuzzy_string_match(stored_routine_name, routine_name)) then

          ! Already have it, decrement counters
          mem_bytes(i_name, thread_id) = &
               mem_bytes(i_name, thread_id) - footprint
          mem_call_balance(i_name, thread_id) = &
               mem_call_balance(i_name, thread_id) - 1

          ! If the balance is now at zero, remove the entry
          if(mem_call_balance(i_name, thread_id) == 0) then
             ! Shift each subsequent routine up by one
             do j_name = i_name + 1, mem_n_ids(thread_id)
                mem_routine_names(j_name-1, thread_id) = &
                     mem_routine_names(j_name, thread_id)
                mem_var_names(j_name-1, thread_id) = &
                     mem_var_names(j_name, thread_id)
                mem_bytes(j_name-1, thread_id) = &
                     mem_bytes(j_name, thread_id)
                mem_call_balance(j_name-1, thread_id) = &
                     mem_call_balance(j_name, thread_id)
             end do
             ! Overwrite the last entry (no longer used)
             mem_routine_names(mem_n_ids(thread_id),thread_id) = '***DELETED***'
             mem_var_names(mem_n_ids(thread_id),thread_id) = '***DELETED***'
             mem_bytes(mem_n_ids(thread_id),thread_id) = 0
             mem_call_balance(mem_n_ids(thread_id),thread_id) = 0
             ! Decrement the number of entries
             mem_n_ids(thread_id) = mem_n_ids(thread_id) - 1
             call utils_assert(mem_n_ids(thread_id) >= 0, &
                  'alloc_leak_check_out: Logic error')
          end if

          found_it = .true.
          exit
       end if

    end do

    ! If not found, look for a "partial match".
    ! e.g. esdf_close::llist will match esdf_init::llist.

    call utils_assert(found_it, &
         'alloc_leak_check_out: You tried to deallocate "'//&
         trim(routine_name)//'::'//trim(var_name)//&
         '", which has either never been allocated or has been deallocated &
         &already -- thread '//trim(utils_int_to_str(thread_id))//&
         ', rank '//trim(utils_int_to_str(pub_my_proc_id)))

  end subroutine alloc_leak_check_out

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function fuzzy_string_match(s1, s2)

    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    character(len=max_id_len), intent(in) :: s1
    character(len=max_id_len), intent(in) :: s2

    ! jd: Local variables
    integer :: i1, i2
    character :: c1, c2
    character(len=*), parameter :: myself = 'fuzzy_string_match'

    ! ------------------------------------------------------------------------

    fuzzy_string_match = .true.
    if(trim(s1) == trim(s2)) return

    i1 = 1
    i2 = 1
    do while(i1 <= max_id_len .or. i2 <= max_id_len)

       if(i1 <= max_id_len) then
          c1 = s1(i1:i1)
       else
          c1 = ' '
       end if

       if(i2 <= max_id_len) then
          c2 = s2(i2:i2)
       else
          c2 = ' '
       end if

       do while(c1 == '[')
          do while (c1 /= ']')
             i1 = i1 + 1
             call utils_assert(i1 <= max_id_len, &
                  myself//': Unterminated [ in s1="'//trim(s1)//'"')
             c1 = s1(i1:i1)
          end do
          i1 = i1 + 1
          if(i1 <= max_id_len) then
             c1 = s1(i1:i1)
          else
             c1 = ' '
          end if
       end do

       do while(c2 == '[')
          do while (c2 /= ']')
             i2 = i2 + 1
             call utils_assert(i2 <= max_id_len, &
                  myself//': Unterminated [ in s2="'//trim(s2)//'"')
             c2 = s2(i2:i2)
          end do
          i2 = i2 + 1
          if(i2 <= max_id_len) then
             c2 = s2(i2:i2)
          else
             c2 = ' '
          end if
       end do

       if(c1 /= c2) then
          fuzzy_string_match = .false.
          return
       end if

       i1 = i1 + 1
       i2 = i2 + 1

    end do

  end function fuzzy_string_match

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module alloc_checker
