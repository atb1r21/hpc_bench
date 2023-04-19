! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
#ifdef ACCELRYS
!=============================================================================!
! This module contains a dummy interface to a timer routine which could be    !
! used in future to provide a breakdown of the time spent in various routines !
! in the code.                                                                !
!-----------------------------------------------------------------------------!
! This version by Peter Haynes, 2 March 2005                                  !
!=============================================================================!

module timer

  use constants
  use rundat, only : pub_timings_level

  implicit none

  private

  public :: timer_clock
  public :: timer_check_iteration_time
  public :: wrappers_etime
  public :: timer_pin_subroutines_to_top

  real(kind=DP) :: start_time

contains

  subroutine timer_clock(name,option,work,multiplier)

    use comms
    use utils, only: utils_abort

    implicit none

    ! Arguments
    character(len=*), intent(in) :: name
    integer, intent(in) :: option
    real(kind=DP), optional, intent(in) :: work
    integer, intent(in), optional :: multiplier

    ! Local variables
    real(kind=DP) :: tot_time
    character(len=80) :: tot_str
!$  integer, external :: omp_get_max_threads

    if (pub_timings_level == 0 .and. (option == 1 .or. option == 2)) return

    select case (option)

    case (0)

       start_time = wrappers_etime()

    case (1:2)

       ! Do nothing

    case (3)

       tot_time = wrappers_etime() - start_time
       call comms_reduce('MAX',tot_time)

       if (pub_on_root) then
          write(tot_str,'(a,f12.3,a,i6,a)') 'TOTAL TIME:', tot_time, &
            's on ', pub_total_num_procs, ' proc(s)'
!$        write(tot_str,'(a,a,i4,a)') trim(adjustl(tot_str)), ' each with ', &
!$             omp_get_max_threads(),' thread(s)'
          write(stdout,'(a)')
          write(stdout,'(a)') tot_str
       end if

    case default

       call utils_abort('Error in timer_clock(): unknown option ',option)

    end select

  end subroutine timer_clock

  subroutine timer_check_iteration_time(loop_name,option,break_loop,global_mpi_loop)
    implicit none

    ! arguments
    character(len=*), intent(in) :: loop_name
    character(len=*), intent(in) :: option
    logical, optional, intent(out) :: break_loop
    logical, optional, intent(in) :: global_mpi_loop ! this is only optional at start of loop

    if (present(break_loop)) break_loop = .false.

  end subroutine timer_check_iteration_time

  real(kind=DP) function wrappers_etime()

#ifdef MPI
    use comms, only: MPI_WTIME
#endif
    implicit none

#ifdef MPI
    wrappers_etime = MPI_WTIME()
#else
    ! Local variable
    real(kind=kind(0.0)) :: tmp_cpu_time   ! Intrinsic call is single-precision

    call cpu_time(tmp_cpu_time)
    wrappers_etime = real(tmp_cpu_time,DP)
#endif

  end function wrappers_etime


  subroutine timer_pin_subroutines_to_top
    !==========================================================================!
    ! No-op in ACCELRYS build.                                                 !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2019.                              !
    !==========================================================================!

    implicit none

    ! -------------------------------------------------------------------------

  end subroutine timer_pin_subroutines_to_top

end module timer

#else
!######################################################################!
!#                                                                    #!
!#                                                                    #!
!#              This module belongs to the package "ONES"             #!
!#                                                                    #!
!#                                                                    #!
!#                          Oswaldo Dieguez                           #!
!#                 The Theory of Condensed Matter Group               #!
!#            Cavendish Laboratory (University of Cambridge)          #!
!#                           Madingley Road                           #!
!#                         Cambridge  CB3 0HE                         #!
!#                           United Kingdom                           #!
!#                                                                    #!
!#                        odl22@phy.cam.ac.uk                         #!
!#                                                                    #!
!#                                                                    #!
!#                          27 November 2001                          #!
!#                                                                    #!
!#                                                                    #!
!######################################################################!

! USAGE:
! This module can be used to measure the computation time for different
! parts of a program.
! Follow these steps:
! 1) Include the line "USE timer" in the subroutines where timings are
!    to be performed.
! 2) Include the line "CALL clock('TOTAL TIME:',0)" before the first
!    executable order in the main program. Include the line
!    "CALL clock('TOTAL TIME:',3)" after the last executable order in
!    the main program.
! 3) Include the lines "CALL clock(string,1)" and "CALL clock(string,2)"
!    delimiting the beginning and end of a part of the program to be
!    timed. "string" is a string of characters chosen by the user that
!    will identify on output the part of the program timed.

!=====================================================================!
! Adapted for the parallel version of ONETEP by Chris-Kriton Skylaris !
! in November 2003.                                                   !
! Revised for MPI by Peter Haynes, July 2004                          !
! Now also counts amount of work (floating points operations or bytes !
! transferred).                                                       !
! Revised for OpenMP by Karl Wilkinson, July 2013                     !
!=====================================================================!

module timer

  use constants, only: DP
  use rundat, only: pub_timings_level

  implicit none

  private

  integer, parameter :: name_len = 50

  type timing_info
     real(kind=DP) :: num_calls     ! no. of calls
     real(kind=DP) :: cum_work      ! cumulative amount of work
     real(kind=DP) :: cum_time      ! cumulative time
     real(kind=DP) :: self_time     ! cum. time - cum. time of timed children
     real(kind=DP) :: last_time     ! last time
     integer       :: num_procs     ! number of procs that executed this subr.
     integer       :: thread_count
     integer       :: print_index
     integer       :: multiplier
     character(len=name_len) :: name      ! name label
  end type timing_info

  integer, parameter :: max_threads = 128     ! no. of possible threads used
  integer, parameter :: max_names = 200  ! no. of different timings possible
  integer, parameter :: max_call_stack_depth = 1000 ! max. depth of call stack
  integer :: num_names(max_threads) = 0
  ! jd: max_threads+1 used by kaw to reduce over threads
  ! jd: max_threads+2 used by jd for storing the final, ordered version
  type(timing_info) :: routines(max_threads+2,max_names)
  integer :: call_stack_depth(max_threads) = 0
  integer :: call_stack(max_threads,max_call_stack_depth)
  real(kind=DP) :: loop_runtime
  logical :: loop_runtime_in_use = .false.
  character(len=80) :: loop_runtime_current_name

  public :: timer_clock
  public :: timer_check_iteration_time
  public :: wrappers_etime
  public :: timer_pin_subroutines_to_top

  ! jd: Names of subroutines that take a very short time to execute and are
  !     executed very many times. By pinning their names to the top of the
  !     timer name list we reduce timing overheads.
  character(len=*), parameter :: pinned_subroutines(11) = (/&
       'swri_swop_calc_all_in_ppd_fast  ', &
       'remote_unpack_ngwf              ', &
       'remote_obtain_dkn_block         ', &
       'swri_metric_mat_cheb_int_product', &
       'swri_metric_mat_cheb_expand     ', &
       'swri_metric_mat_eval_swpot      ', &
       'basis_add_function_to_box       ', &
       'basis_dot_function_with_box     ', &
       'remote_unpack_ngwf_no_sphere    ', &
       'swri_swop_calc_all_in_ppd       ', &
       'swri_metric_mat_dispatch        '/)

contains

  subroutine timer_clock(name,option,work,multiplier)

    use comms, only: comms_barrier, pub_on_root, comms_reduce, &
         pub_my_proc_id, pub_total_num_procs, pub_root_proc_id, comms_send, &
         comms_irecv, comms_wait
    use constants, only: DP, stdout, CRLF
    use rundat, only: pub_timings_level, pub_timings_order
    use utils, only: utils_abort, utils_assert, utils_alloc_check, &
         utils_dealloc_check, utils_banner

    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in) :: option
    real(kind=DP), optional, intent(in) :: work
    integer, intent(in), optional :: multiplier

    ! cks: internal declarations
    integer :: iname
    real(kind=DP) :: tot_time, av_tot_time, max_tot_time
    real(kind=DP) :: percent
    real(kind=DP) :: sum_time, sum_self
    real(kind=DP) :: time_passed
    real(kind=DP) :: sum_work
    integer :: proc_count
    real(kind=DP) :: sum_calls
    real(kind=DP) :: calls_per_proc      ! no. of calls per proc
    integer :: parent
    character(len=6) :: calls_string
    character(len=8) :: rate_string
    character(len=name_len) :: local_name
    character(len=name_len) :: remote_name
    character(len=80) :: tot_str
    integer :: thread_id, num_threads, out, num_out, thread, jname
    real(kind=DP) :: inv_num_threads
    logical :: match
    logical :: already_have_it
    integer :: in
    integer :: local_i, remote_i
    real(kind=DP) :: remote_num_calls, remote_cum_time, remote_cum_work, &
         remote_self_time
    integer :: remote_multiplier
    integer :: sender
    integer :: sum_procs
    integer :: rh, all_rh
    integer :: recv_handles(pub_total_num_procs * 6), send_handles(7)
    integer :: print_idx, best_idx
    real(kind=DP), allocatable :: remote_num_calls_all(:)
    real(kind=DP), allocatable :: remote_cum_time_all(:)
    real(kind=DP), allocatable :: remote_cum_work_all(:)
    real(kind=DP), allocatable :: remote_self_time_all(:)
    character, allocatable :: remote_name_all(:)
    integer, allocatable :: remote_multiplier_all(:)
    character :: all_my_names(max_names*name_len)
    real(kind=DP) :: all_my_num_calls(max_names)
    real(kind=DP) :: all_my_cum_time(max_names)
    real(kind=DP) :: all_my_cum_work(max_names)
    real(kind=DP) :: all_my_self_time(max_names)
    integer :: all_my_multiplier(max_names)
    integer :: nums_out_on_remote(pub_total_num_procs-1)
    integer :: max_nums_out_on_remote, total_nums_out_on_remote
    integer :: send_handle_idx
    real(kind=DP) :: best_val, cur_val
    logical :: timer_overhead
    real(kind=DP) :: timer_call_overhead
    integer :: ierr
    integer :: s
    integer :: loc_multiplier
    integer :: final_multiplier
    real(kind=DP) :: t1,t2,t3,t4,t5
!$  integer, external :: omp_get_thread_num, omp_get_num_threads
!$  logical, external :: omp_in_parallel

    if (pub_timings_level == 0 .and. (option == 1 .or. option == 2)) return

#ifdef BARRIER_IN_TIMER
    call comms_barrier()
#endif

    thread_id = 1
    num_threads = 1
!$  thread_id = omp_get_thread_num() + 1
!$  num_threads = omp_get_num_threads()
    inv_num_threads = 1.0_DP/real(num_threads,kind=DP)

    if(present(multiplier)) then
       loc_multiplier = multiplier
    else
       loc_multiplier = 1
    end if

    select case (option)

    case (0)
       ! cks: initialise variables corresponding to the start of time

       num_names(:) = 1
       routines(:,:)%cum_work = 0.0_DP
       routines(:,:)%cum_time = 0.0_DP
       routines(:,:)%self_time = 0.0_DP
       routines(:,1)%num_calls = 1.0_DP
       routines(:,1)%name = name
       routines(:,1)%last_time = wrappers_etime()
       routines(:,:)%thread_count = 1
       routines(:,:)%multiplier = 1

       ! jd: Manually put 'total_time' as the first entry in the call stack
       call_stack(:,1) = 1
       call_stack_depth(:) = 1

       out = max_threads + 1

       return

    case (1)

          ! cks: Loop over all names used so far, looking for this one
          do iname=1,num_names(thread_id)
             ! cks: break out if current time "tag" has been used before
             if (name == routines(thread_id,iname)%name) exit
          end do

          ! jd: If the loop ran to completion, iname = num_names + 1 and we
          !     don't have this "tag" yet. In both cases, update this tag,
          !     but first make sure we don't go past the end of the array

          ! cks: make sure you don't exceed allocated memory
          if (iname > max_names) then
          call utils_abort('Error in timer_clock: &
               &maximum number of timing tags exeeded: ', max_names)
          end if

          routines(thread_id,iname)%name = name
          routines(thread_id,iname)%last_time = wrappers_etime()
          routines(thread_id,iname)%num_calls = &
             routines(thread_id,iname)%num_calls + 1 * loc_multiplier
          routines(thread_id,iname)%thread_count = num_threads
          routines(thread_id,iname)%multiplier = loc_multiplier

          ! jd: Put the index of this "tag" on top of the call stack
          call_stack_depth(thread_id) = call_stack_depth(thread_id) + 1
          if (call_stack_depth(thread_id) > max_call_stack_depth) then
          call utils_abort('Error in timer_clock: maximum &
               &call stack depth exceeded. Check for unbounded recursion and &
               &erroneous uses of "1" instead of "2" in subroutine exits')
          end if
          call_stack(thread_id,call_stack_depth(thread_id)) = iname

          ! jd: If this was a new tag, increase count
          if (iname == num_names(thread_id) + 1) num_names(thread_id) = &
               num_names(thread_id) + 1

    case (2)

           if (routines(thread_id,call_stack(&
                thread_id,call_stack_depth(thread_id)))%name /= name) then
              call utils_abort('Error in timer_clock: &
                  &call stack inconsistent. Check for mismatched calls to &
                  &timer_clock.'//CRLF//'Timer indicates exit from '//&
                  trim(name)//' while the most recently entered routine was '//&
                  trim(routines(thread_id,call_stack(&
                  thread_id,call_stack_depth(thread_id)))%name))
           end if

           ! cks: Loop over all names used so far, looking for this one
           do iname=1,num_names(thread_id)
              if (name == routines(thread_id,iname)%name) then

                 ! cks: add time to the time already stored for the current "tag"
                 time_passed = wrappers_etime() - routines(thread_id,iname)%last_time
                 time_passed = time_passed * inv_num_threads
                 time_passed = time_passed * loc_multiplier
                 routines(thread_id,iname)%cum_time = &
                    routines(thread_id,iname)%cum_time + time_passed
                 routines(thread_id,iname)%self_time = &
                    routines(thread_id,iname)%self_time + time_passed

                 ! jd: Add the work for this "tag"
                 if (present(work)) routines(thread_id,iname)%cum_work = &
                      routines(thread_id,iname)%cum_work + work

                 ! jd: Pop this finished "tag" from the top of the call stack
                 call_stack(thread_id,call_stack_depth(thread_id)) = -1
                 call_stack_depth(thread_id) = call_stack_depth(thread_id) - 1
                 if (call_stack_depth(thread_id) < 0) then
                    call utils_abort('Error in timer_clock: &
                         &call stack inconsistent. Check for mismatched calls to &
                         &timer_clock.')
                 end if

                 ! jd: Subtract the timing of this routine from the self time
                 !     of the routine which is now the top of the call stack
                 !     (and which is the parent)
                 parent = call_stack(thread_id,call_stack_depth(thread_id))
                 routines(thread_id,parent)%self_time = &
                     routines(thread_id,parent)%self_time - time_passed / inv_num_threads &
                     / real(routines(thread_id,parent)%thread_count,kind=DP)

                 return

              end if
           end do

           if (pub_on_root) write(stdout,'(a)') 'WARNING in timer_clock: &
                &matching tag "',trim(name),'" not found for option 2'

    case (3)
       ! cks: total time the program was running
       tot_time = wrappers_etime() - routines(1,1)%last_time
       av_tot_time = tot_time
       max_tot_time = tot_time
       call comms_reduce('SUM',av_tot_time)
       call comms_reduce('MAX',max_tot_time)
       av_tot_time = av_tot_time / pub_total_num_procs
       if (pub_on_root) then
          write(stdout,'(/a)') utils_banner('-', 'TIMING INFORMATION')
          write(tot_str,'(a,f12.3,a,i6,a)') 'AVERAGE TIME:', av_tot_time, &
            's on ', pub_total_num_procs, ' proc(s)'
          write(stdout,'(a)') tot_str
          write(tot_str,'(a,f12.3,a,i6,a)') 'TOTAL TIME:  ', max_tot_time, &
            's on ', pub_total_num_procs, ' proc(s)'
          write(stdout,'(a)') tot_str
       end if
       routines(1,1)%cum_time = tot_time
       routines(1,1)%self_time = routines(1,1)%self_time + tot_time

       if (pub_timings_level == 0) return

       call comms_barrier
       t1 = wrappers_etime()
       ! kaw: Reduction of thread timings for a proc. This maintains the
       !      original data as we will probably want it later.

       ! Copy data from thread 1.
       out = max_threads + 1
       routines(out,:) = routines(1,:)
       num_out = num_names(1)

       ! Accumulate data from all threads into seperate dimension of array
       if (max_threads.gt.1) then ! jd: @Should this not be 'pub_threads_max'?
          do thread = 2, max_threads ! loop over threads
             ! loop over names of routines for this thread
             do iname = 2, num_names(thread)
                local_name = routines(thread,iname)%name
                match = .false.
                ! loop over names of routines currently stored in out
                do jname = 1, num_out
                   if (local_name == routines(out,jname)%name) then
                      ! kaw: add data to the that stored for the current "tag"
                      routines(out,jname)%num_calls = &
                               routines(out,jname)%num_calls+ &
                               routines(thread,iname)%num_calls
                      routines(out,jname)%cum_time  = &
                               routines(out,jname)%cum_time  + &
                               routines(thread,iname)%cum_time
                      routines(out,jname)%self_time = &
                               routines(out,jname)%self_time + &
                               routines(thread,iname)%self_time
                      routines(out,jname)%cum_work  = &
                               routines(out,jname)%cum_work  + &
                               routines(thread,iname)%cum_work
                      routines(out,jname)%thread_count  = &
                               routines(out,jname)%thread_count  + 1
                      ! jd: %multiplier is not summed over threads
                      match = .true.
                   end if
                end do
                if (.not.match) then
                   ! kaw: add data for the current "tag"
                   num_out = num_out + 1
                   routines(out,num_out)%name = routines(thread,iname)%name
                   routines(out,num_out)%num_calls = &
                      routines(thread,iname)%num_calls
                   routines(out,num_out)%cum_time      = &
                      routines(thread,iname)%cum_time
                   routines(out,num_out)%self_time     = &
                      routines(thread,iname)%self_time
                   routines(out,num_out)%cum_work      = &
                      routines(thread,iname)%cum_work
                   routines(out,num_out)%thread_count  = &
                      routines(thread,iname)%thread_count
                   routines(out,num_out)%multiplier  = &
                      routines(thread,iname)%multiplier
                end if
             end do
          end do

          ! kaw: Average the timings for the threaded routines
          ! loop over names of routines currently stored in out
          do iname = 2, num_out
             !if (routines(out,iname)%thread_count.gt.1) then
                !routines(out,iname)%cum_time  = &
                !   routines(out,iname)%cum_time/routines(out,iname)%thread_count
                !routines(out,iname)%self_time = &
                !   routines(out,iname)%self_time/routines(out,iname)%thread_count
                !routines(out,iname)%cum_work  = &
                !   routines(out,iname)%cum_work/routines(out,iname)%thread_count
             !end if
          end do
       end if

       call comms_barrier
       t2 = wrappers_etime()

       ! ********************************
       ! *** Timing breakdown by proc ***
       ! ********************************
       if (iand(pub_timings_level,4) /= 0) then

          ! cks: print name of each "tag", number of calls, times and
          !      percent of total time
          proc_loop: do proc_count=0,pub_total_num_procs-1

             call comms_barrier

             if (pub_my_proc_id == proc_count) then

                ! cks: print timings regarding pub_my_proc_id
                write(stdout,'(a)') '-----------------------------------------&
                     &---------------------------------------'
                write(stdout,'(a)') '| TAG                                     &
                     &#calls    cpu time   %total     proc |'
                do iname=1,num_names(1)
                   ! jd: Do not show names with '#'
                   if(routines(out,iname)%name(1:1)=='#') cycle
                   ! jd: Do not show names with zero calls (pinned subroutines
                   !     that were never actually called)
                   if(routines(out,iname)%num_calls == 0) cycle

                   percent = 100.0_DP * routines(out,iname)%cum_time / tot_time
                   call internal_write_int(int(routines(out,iname)%num_calls), &
                        calls_string)
                   write(stdout,'(a2,a38,a2,a6,1x,f10.3,a2,f7.3,a2,i8,a2)') &
                        '| ',routines(out,iname)%name,': ',calls_string, &
                        routines(out,iname)%cum_time,'s ',percent,'% ', &
                        pub_my_proc_id,' |'

                end do
                write(stdout,'(a)') '-----------------------------------------&
                     &---------------------------------------'

             end if

          end do proc_loop

       end if

       ! jd: Now ensure all MPI tasks have a consistent view of the routines
       !     called and that they are ordered in the same fashion. This takes
       !     care of scenarios where certain routines are not called on
       !     certain procs or are called in a different order.
       in = max_threads+1
       out = max_threads+2
       if(.not. pub_on_root) then
          ! jd: Tell root how many names we have on this proc
          call comms_send(pub_root_proc_id,num_out, &
               return_handle = send_handles(1), add_to_stack = .false.)
          ! jd: 'Pack' all relevant data into arrays that can be sent in one go
          do iname = 1, num_out
             do s = 1, name_len
                all_my_names((iname-1)*name_len+s) = &
                     routines(in,iname)%name(s:s)
             end do
             all_my_num_calls(iname) = routines(in,iname)%num_calls
             all_my_cum_work(iname) = routines(in,iname)%cum_work
             all_my_cum_time(iname) = routines(in,iname)%cum_time
             all_my_self_time(iname) = routines(in,iname)%self_time
             all_my_multiplier(iname) = routines(in,iname)%multiplier
          end do

          ! jd: Send every local name followed by num_calls, cum_time, etc.
          !     to root.
          call comms_send(pub_root_proc_id, all_my_names(1), &
               length = num_out * name_len, return_handle = send_handles(2), &
               add_to_stack = .false.)
          call comms_send(pub_root_proc_id, all_my_num_calls, length = num_out, &
               return_handle = send_handles(3), add_to_stack = .false.)
          call comms_send(pub_root_proc_id, all_my_cum_work, length = num_out, &
               return_handle = send_handles(4), add_to_stack = .false.)
          call comms_send(pub_root_proc_id, all_my_cum_time, length = num_out, &
               return_handle = send_handles(5), add_to_stack = .false.)
          call comms_send(pub_root_proc_id, all_my_self_time, length = num_out, &
               return_handle = send_handles(6), add_to_stack = .false.)
          call comms_send(pub_root_proc_id, all_my_multiplier, length = num_out, &
               return_handle = send_handles(7), add_to_stack = .false.)

          do send_handle_idx = 1, 6
             call comms_wait(send_handles(send_handle_idx))
          end do

       else
          ! jd: On root receive the names and corresponding data from all procs
          !     and put them in routines(out).
          ! jd: But first copy root's routines to out.
          routines(out,:) = routines(in,:)
          routines(out,:)%num_procs = 1

          ! jd: Receive the numbers of entries from all procs first
          do sender = 1, pub_total_num_procs-1
             call comms_irecv(sender, nums_out_on_remote(sender), &
                  handle = recv_handles(sender))
          end do
          do sender = 1, pub_total_num_procs-1
             call comms_wait(recv_handles(sender))
          end do

          ! jd: Longest list of timings across all procs
          max_nums_out_on_remote = &
               maxval(nums_out_on_remote(1:pub_total_num_procs-1))

          ! jd: Upper bound on all lists collated
          total_nums_out_on_remote = &
               max_nums_out_on_remote * pub_total_num_procs

          allocate(remote_name_all(total_nums_out_on_remote*name_len),stat=ierr)
          call utils_alloc_check('timer_clock','remote_name_all',ierr)
          allocate(remote_num_calls_all(total_nums_out_on_remote),stat=ierr)
          call utils_alloc_check('timer_clock','remote_num_calls_all',ierr)
          allocate(remote_cum_work_all(total_nums_out_on_remote),stat=ierr)
          call utils_alloc_check('timer_clock','remote_cum_work_all',ierr)
          allocate(remote_cum_time_all(total_nums_out_on_remote),stat=ierr)
          call utils_alloc_check('timer_clock','remote_cum_time_all',ierr)
          allocate(remote_self_time_all(total_nums_out_on_remote),stat=ierr)
          call utils_alloc_check('timer_clock','remote_self_time_all',ierr)
          allocate(remote_multiplier_all(total_nums_out_on_remote),stat=ierr)
          call utils_alloc_check('timer_clock','remote_multiplier_all',ierr)

          ! jd: Post receives for data from all remote procs
          rh = 1
          do sender = 1, pub_total_num_procs-1
             call comms_irecv(sender, &
                  remote_name_all(max_nums_out_on_remote*(sender-1)*name_len+1),&
                  length = nums_out_on_remote(sender)*name_len, &
                  handle=recv_handles(rh))
             rh = rh + 1
             call comms_irecv(sender, &
                  remote_num_calls_all(max_nums_out_on_remote * (sender-1)+1), &
                  length=nums_out_on_remote(sender), handle=recv_handles(rh))
             rh = rh + 1
             call comms_irecv(sender, &
                  remote_cum_work_all(max_nums_out_on_remote * (sender-1)+1), &
                  length=nums_out_on_remote(sender), handle=recv_handles(rh))
             rh = rh + 1
             call comms_irecv(sender, &
                  remote_cum_time_all(max_nums_out_on_remote * (sender-1)+1), &
                  length=nums_out_on_remote(sender), handle=recv_handles(rh))
             rh = rh + 1
             call comms_irecv(sender, &
                  remote_self_time_all(max_nums_out_on_remote * (sender-1)+1), &
                  length=nums_out_on_remote(sender), handle=recv_handles(rh))
             rh = rh + 1
             call comms_irecv(sender, &
                  remote_multiplier_all(max_nums_out_on_remote * (sender-1)+1), &
                  length=nums_out_on_remote(sender), handle=recv_handles(rh))
             rh = rh + 1
          end do
          all_rh = rh - 1

          ! jd: Wait for completion
          do rh = 1, all_rh
             call comms_wait(recv_handles(rh))
          end do

          ! jd: Process
          do sender = 1, pub_total_num_procs-1
             do remote_i = 1, nums_out_on_remote(sender)

                ! jd: Extract remote data
                do s = 1, name_len
                   remote_name(s:s) = &
                        remote_name_all((max_nums_out_on_remote*(sender-1) + &
                        remote_i-1) * name_len + s)
                end do
                remote_num_calls = &
                     remote_num_calls_all(max_nums_out_on_remote*(sender-1) + &
                     remote_i)
                remote_cum_work = &
                     remote_cum_work_all(max_nums_out_on_remote*(sender-1) + &
                     remote_i)
                remote_cum_time = &
                     remote_cum_time_all(max_nums_out_on_remote*(sender-1) + &
                     remote_i)
                remote_self_time = &
                     remote_self_time_all(max_nums_out_on_remote*(sender-1) + &
                     remote_i)
                remote_multiplier = &
                     remote_multiplier_all(max_nums_out_on_remote*(sender-1) + &
                     remote_i)

                ! jd: See if we have this routine on root yet
                already_have_it = .false.
                do local_i = 1, num_out

                   if(trim(routines(out,local_i)%name) == trim(remote_name)) then
                      ! jd: Yes. Add remote's data.
                      routines(out,local_i)%num_calls = &
                           routines(out,local_i)%num_calls + remote_num_calls
                      routines(out,local_i)%cum_work = &
                           routines(out,local_i)%cum_work + remote_cum_work
                      routines(out,local_i)%cum_time = &
                           routines(out,local_i)%cum_time + remote_cum_time
                      routines(out,local_i)%self_time = &
                           routines(out,local_i)%self_time + remote_self_time
                      routines(out,local_i)%num_procs = &
                           routines(out,local_i)%num_procs + 1
                      routines(out,local_i)%multiplier = &
                           routines(out,local_i)%multiplier + remote_multiplier
                      already_have_it = .true.
                      exit
                   end if
                end do
                if(.not. already_have_it) then
                   num_out = num_out + 1
                   routines(out,num_out)%name = remote_name
                   routines(out,num_out)%num_calls = remote_num_calls
                   routines(out,num_out)%cum_work = remote_cum_work
                   routines(out,num_out)%cum_time = remote_cum_time
                   routines(out,num_out)%self_time = remote_self_time
                   routines(out,num_out)%num_procs = 1
                   routines(out,num_out)%multiplier = remote_multiplier
                end if
             end do
          end do

          deallocate(remote_multiplier_all,stat=ierr)
          call utils_dealloc_check('timer_clock','remote_multiplier_all',ierr)
          deallocate(remote_self_time_all,stat=ierr)
          call utils_dealloc_check('timer_clock','remote_self_time_all',ierr)
          deallocate(remote_cum_time_all,stat=ierr)
          call utils_dealloc_check('timer_clock','remote_cum_time_all',ierr)
          deallocate(remote_cum_work_all,stat=ierr)
          call utils_dealloc_check('timer_clock','remote_cum_work_all',ierr)
          deallocate(remote_num_calls_all,stat=ierr)
          call utils_dealloc_check('timer_clock','remote_num_calls_all',ierr)
          deallocate(remote_name_all,stat=ierr)
          call utils_dealloc_check('timer_clock','remote_name_all',ierr)

       end if

       ! cks: print maximum timings across all procs

       call comms_barrier
       t3 = wrappers_etime()

       ! jd: The below is no longer necessary.
       ! cks: see if every proc has timings for the same tasks
!       min_global_num_names = num_names(1)
!       call comms_reduce('MIN', min_global_num_names)
!       max_global_num_names = num_names(1)
!       call comms_reduce('MAX', max_global_num_names)
!       if (min_global_num_names /= max_global_num_names .and. pub_on_root) &
!            write(stdout,'(a)') 'WARNING in timer_clock: inconsistent timing &
!                 &data'

       ! **********************************
       ! *** Usual (cumulative) timings ***
       ! **********************************
       if (iand(pub_timings_level,1) /= 0) then

          if (pub_on_root) then

             write(stdout,'(a)') utils_banner('=', &
                  'AVERAGE TIMINGS FROM ALL PROCS (CUMULATIVE)')
             write(stdout,'(a)') '|| TAG                                    &
                  &#calls    cpu time   %total  Gflops ||'

             ! jd: Determine how to order the values
             routines(out,:)%print_index = 99999
             do print_idx=1, num_out
                best_idx = -1
                best_val = huge(1.0_DP)
                do iname=1, num_out
                   if(routines(out,iname)%print_index /= 99999) cycle

                   sum_calls = routines(out,iname)%num_calls
                   sum_time = routines(out,iname)%cum_time

                   if(trim(pub_timings_order) == trim('NCALLS')) then
                      cur_val = sum_calls
                   else if(trim(pub_timings_order) == trim('TIME')) then
                      cur_val = sum_time
                   else if(trim(pub_timings_order) == trim('NONE')) then
                      cur_val = iname
                   else if(trim(pub_timings_order) == trim('NAME')) then
                      cur_val = real(ichar(routines(out,iname)%name(1:1)), &
                           kind=DP)
                   else
                      call utils_abort('Unrecognized timings_order.')
                   end if

                   if(cur_val < best_val) then
                      best_val = cur_val
                      best_idx = iname
                   end if
                end do
                call utils_assert(best_idx > 0, &
                     'Internal error [1] in timer_mod.')
                routines(out,best_idx)%print_index = print_idx
             end do

             ! jd: Print in order
             do print_idx=1, num_out
                iname = minloc(routines(out,:)%print_index,1)
                routines(out,iname)%print_index = 99999

                sum_calls = routines(out,iname)%num_calls
                sum_time = routines(out,iname)%cum_time
                sum_work = routines(out,iname)%cum_work
                sum_self = routines(out,iname)%self_time
                sum_procs = routines(out,iname)%num_procs
                local_name = routines(out,iname)%name
                if(local_name(1:1)=='#') cycle ! Do not show names with '#'
                ! jd: Do not show names with zero calls (pinned subroutines
                !     that were never actually called)
                if(routines(out,iname)%num_calls == 0) cycle

                percent = 100.0_DP * sum_time / (av_tot_time * sum_procs)

                if (sum_time > epsilon(1.0_DP) .and. sum_work > 0.0_DP) then
                   !write(rate_string,'(f8.3)') 1.0e-9 * sum_work * &
                   !     pub_total_num_procs / sum_time
                   call internal_write_real(1.0e-9 * sum_work * &
                        sum_procs / sum_time,rate_string)
                else
                   rate_string = '  ------'
                end if

                calls_per_proc = sum_calls/sum_procs
                call internal_write_real_as_int(calls_per_proc,calls_string)
                write(stdout,'(a3,a37,a2,a6,1x,f10.2,a2,f7.3,a1,a8,a3)', &
                     advance='no') &
                     '|| ',local_name,': ',calls_string,&
                        sum_time/sum_procs,'s ',percent,'%',&
                        rate_string, ' ||'
                if(routines(out,iname)%num_procs /= pub_total_num_procs) then
                   write(stdout,'(a,i0,a)') ' (', routines(out,iname)%num_procs,&
                        ' procs only)'
                else
                   write(stdout,'(a)') ''
                end if
             end do
          end if

          if (pub_on_root) write(stdout,'(a)') repeat('=',80)

       end if

       call comms_barrier
       t4 = wrappers_etime()

       ! ******************
       ! *** Self-times ***
       ! ******************
       if (iand(pub_timings_level,2) /= 0) then

          call comms_barrier

          if (pub_on_root) then

             write(stdout,'(a)') utils_banner('=', &
                  'AVERAGE TIMINGS FROM ALL PROCS (SELF)')
             write(stdout,'(a)') '++ TAG                                    &
                  &#calls    cpu time   %total  Gflops ++'

             routines(out,1)%name = 'main_program_(onetep.F90)'

             ! jd: Determine how to order the values
             routines(out,:)%print_index = 99999
             do print_idx=1, num_out
                best_idx = -1
                best_val = huge(1.0_DP)
                if(trim(routines(out,print_idx)%name) == &
                     trim('#timer_clock overhead for 1000 calls')) then
                   timer_call_overhead = &
                        routines(out,print_idx)%cum_time / 1000.0_DP / &
                        routines(out,print_idx)%num_procs
                end if

                do iname=1, num_out
                   if(routines(out,iname)%print_index /= 99999) cycle

                   sum_calls = routines(out,iname)%num_calls
                   sum_self = routines(out,iname)%self_time

                   if(trim(pub_timings_order) == trim('NCALLS')) then
                      cur_val = sum_calls
                   else if(trim(pub_timings_order) == trim('TIME')) then
                      cur_val = sum_self
                   else if(trim(pub_timings_order) == trim('NONE')) then
                      cur_val = iname
                   else if(trim(pub_timings_order) == trim('NAME')) then
                      cur_val = real(ichar(routines(out,iname)%name(1:1)), &
                           kind=DP)
                   else
                      call utils_abort('Unrecognized timings_order.')
                   end if

                   if(cur_val < best_val) then
                      best_val = cur_val
                      best_idx = iname
                   end if
                end do
                call utils_assert(best_idx > 0, &
                     'Internal error [2] in timer_mod.')
                routines(out,best_idx)%print_index = print_idx
             end do

             timer_overhead = .false.

             ! jd: Print in order
             do print_idx=1, num_out
                iname = minloc(routines(out,:)%print_index,1)
                routines(out,iname)%print_index = 99999

                sum_calls = routines(out,iname)%num_calls
                sum_time = routines(out,iname)%cum_time
                sum_work = routines(out,iname)%cum_work
                sum_self = routines(out,iname)%self_time
                sum_procs = routines(out,iname)%num_procs
                local_name = routines(out,iname)%name
                if(local_name(1:1)=='#') cycle ! Do not show names with '#'
                ! jd: Do not show names with zero calls (pinned subroutines
                !     that were never actually called)
                if(routines(out,iname)%num_calls == 0) cycle

                percent = 100.0_DP * sum_self / (av_tot_time * sum_procs)

                if (sum_time > epsilon(1.0_DP) .and. sum_work > 0.0_DP) then
                   call internal_write_real(1.0e-9 * sum_work * &
                        sum_procs / sum_time,rate_string)
                else
                   rate_string = '  ------'
                end if

                calls_per_proc = sum_calls/sum_procs
                call internal_write_real_as_int(calls_per_proc,calls_string)
                write(stdout,'(a3,a37,a2,a6,1x,f10.2,a1,f8.3,a1,a8,a3)', &
                     advance='no') &
                     '++ ',local_name,': ',calls_string,&
                     sum_self/sum_procs,'s',percent,'%',&
                     rate_string, ' ++'
                final_multiplier = routines(out,iname)%multiplier/sum_procs
                ! jd: Each timer_clock call pair costs about 1-5 us. If the
                !     overhead of measuring short subroutines is >2% and this is
                !     for a non-trivial routine (measured selftime per proc>1s),
                !     warn the developer they're wasting too many cycles on
                !     measurements.
                if(timer_call_overhead / final_multiplier > &
                     0.02_DP * sum_self/sum_procs/calls_per_proc .and. &
                     sum_self/sum_procs > 1.0_DP) then
                   write(stdout,'(a,f6.1,a)',advance='no') ' (!) overhead: ~', &
                        100.0_DP * timer_call_overhead / final_multiplier / &
                        (sum_self/sum_procs/calls_per_proc),'%'
                   timer_overhead = .true.
                end if
                if(routines(out,iname)%num_procs /= pub_total_num_procs) then
                   write(stdout,'(a,i0,a)') ' (',routines(out,iname)%num_procs,&
                        ' procs only)'
                else
                   write(stdout,'(a)') ''
                end if
             end do

             write(stdout,'(a)') repeat('=',80)

             if(timer_overhead) then
                write(stdout,'(a,f6.2,a)') '(!) WARNING: You''re wasting&
                     & precious cycles calling timer. '//CRLF//'            &
                     & Each pair of timer_clock calls itself cost ', &
                     timer_call_overhead*1E6_DP,' microseconds!'
             end if

          end if


       end if

       call comms_barrier
       t5 = wrappers_etime()

       if(pub_on_root) then
          write(stdout,'(a)') 'Overhead of generating timer report: '
          write(stdout,'(a,i0,a)') '- ', int((t2-t1)*1000.0_DP), &
               ' ms averaging over threads,'
          write(stdout,'(a,i0,a)') '- ', int((t3-t2)*1000.0_DP), &
               ' ms ensuring consistency across MPI ranks,'
          write(stdout,'(a,i0,a)') '- ', int((t5-t3)*1000.0_DP), &
               ' ms printing in order.'
       end if

    case default

       call utils_abort('Error in timer_clock(): unknown option ',option)

    end select

  contains

    ! pdh: compressed format for large integers
    subroutine internal_write_int(n,text)

      implicit none

      integer, intent(in) :: n
      character(len=6), intent(out) :: text

      if (n < 0) then  ! integer overflow has occurred
         write(text,'(a5,1x)') 'OVFLW'
      else if (n < 100000) then
         write(text,'(i5,1x)') n
      else if (n < 999950) then
         write(text,'(f5.1,a1)') n*1.0e-3_DP,'k'
      else if (n < 9999500) then
         write(text,'(f5.3,a1)') n*1.0e-6_DP,'M'
      else if (n < 99995000) then
         write(text,'(f5.2,a1)') n*1.0e-6_DP,'M'
      else if (n < 999950000) then
         write(text,'(f5.1,a1)') n*1.0e-6_DP,'M'
      else
         write(text,'(i5,a1)') nint(n*1.0e-6_DP),'M'
      end if

    end subroutine internal_write_int

    ! ndmh: compressed format for large reals
    subroutine internal_write_real(rn,text)

      implicit none

      real(kind=DP), intent(in) :: rn
      character(len=8), intent(out) :: text

      if (rn < 0.0_DP) then  ! integer overflow has occurred
         write(text,'(a5,3x)') 'OVFLW'
      else if (rn < 999.9995_DP) then
         write(text,'(1x,f7.3)') rn
      else if (rn < 99999.5_DP) then
         write(text,'(f7.3,a1)') rn*1.0e-3_DP,'k'
      else if (rn < 999999.5_DP) then
         write(text,'(f7.2,a1)') rn*1.0e-3_DP,'k'
      else if (rn < 999995.0_DP) then
         write(text,'(f7.3,a1)') rn*1.0e-6_DP,'M'
      else if (rn < 9999500.0_DP) then
         write(text,'(f7.2,a1)') rn*1.0e-6_DP,'M'
      else
         write(text,'(f7.1,a1)') rn*1.0e-6_DP,'M'
      end if

    end subroutine internal_write_real

    ! jd: compressed format for large integers represented as reals
    subroutine internal_write_real_as_int(n,text)

      implicit none

      real(kind=DP), intent(in) :: n
      character(len=6), intent(out) :: text

      if (n < 100000.0_DP) then
         write(text,'(i5,1x)') int(n)
      else if (n < 999950.0_DP) then
         write(text,'(f5.1,a1)') n*1.0e-3_DP,'k'
      else if (n < 9999500.0_DP) then
         write(text,'(f5.3,a1)') n*1.0e-6_DP,'M'
      else if (n < 99995000.0_DP) then
         write(text,'(f5.2,a1)') n*1.0e-6_DP,'M'
      else if (n < 999950000.0_DP) then
         write(text,'(f5.1,a1)') n*1.0e-6_DP,'M'
      else if (n < 9999500000.0_DP) then
         write(text,'(f5.3,a1)') n*1.0e-9_DP,'G'
      else if (n < 99995000000.0_DP) then
         write(text,'(f5.2,a1)') n*1.0e-9_DP,'G'
      else if (n < 999950000000.0_DP) then
         write(text,'(f5.1,a1)') n*1.0e-9_DP,'G'
      else
         write(text,'(a6)') 'OVRFLW'
      end if

    end subroutine internal_write_real_as_int


  end subroutine timer_clock


  subroutine timer_check_iteration_time(loop_name,option,break_loop,global_mpi_loop)
  !=====================================================================!
  ! Routine for checking if the next iteration of a computationally     !
  ! intensive loop will have chance to finish. Works by comparing the   !
  ! timing of the previous iteration to the remaining time.             !
  !                                                                     !
  ! Be careful to ensure that a check is made that this is not the last !
  ! iteration making a warning message unnecessary.                     !
  !                                                                     !
  ! The routine does not stop the calculation, it instead returns a     !
  ! flag so that necessary data can be saved before runtime expires.    !
  !                                                                     !
  ! Flag global_mpi_loop indicates that loop is across all procs and    !
  ! if one proc detects that execution should stop, then all procs must !
  ! stop. A comms_reduce is used to compare break_loop on all procs.    !
  ! This flag must *NOT* be set if not all the procs will complete this !
  ! loop iteration else comms_reduce will hang.                         !
  !                                                                     !
  ! Only one loop may use the timer at a time. It should therefore be   !
  ! put around only the most computationally intensive i.e. outer loop. !
  !                                                                     !
  !=====================================================================!
  !                                                                     !
  ! To use, insert:                                                     !
  !   >>  call timer_check_iteration_time('loop_name', 'start')         !
  ! inside and at the start of the loop, then use:                      !
  !   >>  call timer_check_iteration_time('loop_name', 'stop', &        !
  !                                  & global_mpi_loop, break_loop)     !
  ! at the end of the loop.                                             !
  ! Outside (and after) the loop, use:                                  !
  !   >>  call timer_check_iteration_time('loop_name', 'cancel')        !
  ! to reset the loop. This is a clean way to stop the timer if         !
  ! convergence is reached.                                             !
  !                                                                     !
  ! The flag break_loop is set true if insufficient runtime remains.    !
  !                                                                     !
  ! Written by Robert Bell Aug 2013                                     !
  !=====================================================================!


    use comms, only: comms_reduce, pub_on_root
    use constants, only: DP, stdout
    use rundat, only: pub_max_runtime
    use utils, only: utils_abort

    implicit none

    ! arguments
    character(len=*), intent(in) :: loop_name
    character(len=*), intent(in) :: option
    logical, optional, intent(out) :: break_loop
    logical, optional, intent(in) :: global_mpi_loop ! this is only optional at start of loop

    ! internal
    real(kind=DP) :: current_runtime

    select case (option)

      case ('start')
        if (loop_runtime_in_use) then
          ! can only use timer over one loop at once
          write(stdout,'(a,/4x,a,/a)') 'WARNING in timer_check_iteration_time():', &
               & 'timer already in use for loop '//trim(loop_name)//' but set also for '//&
               & trim(loop_runtime_current_name),'. Cannot use in multiple nested loops.'
        endif

        loop_runtime_current_name = trim(loop_name)

        ! initialise the timer with current time
        loop_runtime = wrappers_etime()
        loop_runtime_in_use = .true.

      case ('cancel')
        loop_runtime_in_use = .false.

      case ('stop')
        if (.not. loop_runtime_in_use) then
          call utils_abort('Error in timer_check_iteration_time(): timer not started. &
              & Timer must be initialised at the start of the loop.')
        endif

        ! get loop iteration runtime
        loop_runtime = wrappers_etime() - loop_runtime
        loop_runtime_in_use = .false.

        current_runtime = wrappers_etime() - routines(1,1)%last_time

        ! compare to iteration time to runtime left
        ! (ignore if pub_max_runtime less than zero)
        if (pub_max_runtime > 0.0_DP .and. (pub_max_runtime - current_runtime < loop_runtime) ) then
          break_loop = .true.
        else
          break_loop = .false.
        endif

        ! only break the loop if all procs enter the same loop
        ! (if there is any race condition between procs, must be close to max time:
        ! break loop anyway)
        if (present(global_mpi_loop)) then
          if (global_mpi_loop) then
            call comms_reduce('OR',break_loop)

            ! only write out message if all procs are involved
            if (break_loop .and. pub_on_root) then
              write(stdout,'(//a)') repeat('x',80)
              write(stdout,'(a/)') repeat('x',80)
              write(stdout,'(a)') "WARNING:"
              write(stdout,'(2x,a,f10.0,a)') "Last "//loop_name//" iteration took ", &
                & loop_runtime, " seconds;"
              write(stdout,'(2x,a,f10.0,a)') "Remaining runtime is ", &
                & pub_max_runtime-current_runtime, " seconds."

              if (pub_max_runtime-current_runtime < 0.0_DP) then
                write(stdout,'(2x,a)') "User allocated runtime expired."
              else
                write(stdout,'(2x,a)') "User allocated runtime likely to expire before next &
                  &iteration is completed."
              endif
              write(stdout,'(/a)') repeat('x',80)
              write(stdout,'(a//)') repeat('x',80)
            endif
          endif

        ! global_mpi_loop must be defined! It's not optional!
        else
          call utils_abort('Error in timer_check_iteration_time: you must indicate &
            &whether the loop is MPI safe!')
        endif

      case default

       call utils_abort('Error in timer_check_iteration_time(): unknown option '//option)

    end select


  end subroutine timer_check_iteration_time

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine timer_pin_subroutines_to_top
    !==========================================================================!
    ! Ensures subroutines that are known to take little time and are called    !
    ! very many times wind up at the top of timer's list of subroutines. This  !
    ! decreases timer overhead, as they are then found faster in timer_clock() !
    ! calls.                                                                   !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   We do this on all threads, so this must be called *after*              !
    !   rundat_threads_init(). Also we must be called after get_rundat() for   !
    !   pub_timings_level to be set (needed in timer_clock()).                 !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2019.                              !
    !==========================================================================!

    use rundat, only: pub_threads_max

    implicit none

    integer :: i
    integer :: thread_id
!$  integer, external :: omp_get_thread_num

    ! -------------------------------------------------------------------------

    ! jd: This is not a parallelised loop. All iterations are meant to be
    !     executed on all OMP threads.
!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(NONE)                      &
!$OMP PRIVATE(i, thread_id)                                                    &
!$OMP SHARED(routines, num_names)
    do i = 1, size(pinned_subroutines)
       call timer_clock(trim(pinned_subroutines(i)),1)
       call timer_clock(trim(pinned_subroutines(i)),2)

       thread_id = 1
!$     thread_id = omp_get_thread_num() + 1

       ! jd: Reset num calls for the pinned subroutines to 0 so the ones that
       !     do not actually get called are not shown in timer output.
       routines(thread_id,num_names(thread_id))%num_calls = 0

    end do
!$OMP END PARALLEL

  end subroutine timer_pin_subroutines_to_top

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function wrappers_etime()

#ifdef MPI
    use comms, only: MPI_WTIME
#endif

    implicit none

#ifdef MPI
    wrappers_etime = MPI_WTIME()
#else
    ! internal variables
    integer, parameter :: SP  = KIND(1.0)
    real(kind=SP) :: delta,tarray(2),etime

    delta = etime(tarray)
    wrappers_etime = real(delta,DP)
#endif

  end function wrappers_etime

end module timer
#endif




