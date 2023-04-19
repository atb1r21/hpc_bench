! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!-----------------------------------------------------------------------!
!                                                                       !
! Neighbour list module                                                 !
!                                                                       !
! Written by Jacek Dziedzic in November 2013.                           !
!                                                                       !
!-----------------------------------------------------------------------!
!                                                                       !
! This module implements a Verlet neighbour list.                       !
!                                                                       !
!-----------------------------------------------------------------------!

module neighbour_list

  use constants, only: DP

  implicit none
  private

  ! Public types
  public :: NL_NEIGHBOUR_LIST

  ! Public subroutines
  public :: neighbour_list_init
  public :: neighbour_list_free
  public :: neighbour_list_populate
  public :: neighbour_list_init_from_sparse
  public :: neighbour_list_for_atom_list
  public :: neighbour_list_for_local_atoms
  public :: neighbour_list_procs_of_neighbours
  public :: neighbour_list_are_neighbours

  ! Example:
  !   Neighbour list for 1000 atoms, each can have up to 200 neighbours.
  !   n_elements: 1000, max_total_neighbours = 200000.
  type NL_NEIGHBOUR_LIST
     logical :: populated
     integer :: n_elements
     integer :: max_total_neighbours
     integer :: listunit = 0
     character(len=16) :: list_name
     integer, allocatable :: first_idx(:)
     integer, allocatable :: last_idx(:)
     integer, allocatable :: n_neighbours(:)
     integer, allocatable :: neighbours(:)
  end type NL_NEIGHBOUR_LIST

  ! ---------------------------------------------------------------------------

contains

  ! ---------------------------------------------------------------------------
  ! ----- P U B L I C   S U B R O U T I N E S   A N D   F U N C T I O N S -----
  ! ---------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine neighbour_list_init(list, list_name, n_elements, &
       max_neighbours_per_elem, listunit)
    !==========================================================================!
    ! Initialises a fresh neighbour list object.                               !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   list  (out):     The neighbour list object that is initialised on      !
    !                    return.                                               !
    !   list_name (in):  User-defined name associated with the list.           !
    !   n_elements (in): Number of entries (e.g. atoms) in the index.          !
    !   max_neighbours_per_elem (in): Maximum average number of neighbours per !
    !                                 element. Used for allocation.            !
    !   listunit (in, opt): An open unit to where the output of                !
    !                       neighbour_list_list is to be directed.             !
    !                       Defaults to stdout if omitted.                     !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, November 2013                                 !
    !==========================================================================!

    use constants, only: stdout
    use utils, only: utils_alloc_check, utils_assert

    implicit none

    ! Arguments
    type(NL_NEIGHBOUR_LIST), intent(out) :: list
    character(len=*)                 :: list_name
    integer, intent(in)              :: n_elements
    integer, intent(in)              :: max_neighbours_per_elem
    integer, intent(in), optional    :: listunit

    ! Local variables
    integer :: ierr

    ! ------------------------------------------------------------------------

    call utils_assert(n_elements>=1, &
         'Cannot have neighbour_lists with fewer than 1 element.')
    call utils_assert(max_neighbours_per_elem > 0, &
         'Cannot have neighbour lists with max_neighbours_per_elem < 1')

    list%populated = .false.
    list%n_elements = n_elements
    list%max_total_neighbours = n_elements * max_neighbours_per_elem
    list%list_name = list_name
    if(present(listunit)) then
       list%listunit = listunit
    else
       list%listunit = stdout
    end if

    ! Allocate index
    allocate(list%first_idx(n_elements),stat=ierr)
    call utils_alloc_check('neighbour_list_init','list%first_idx',ierr)
    allocate(list%last_idx(n_elements),stat=ierr)
    call utils_alloc_check('neighbour_list_init','list%last_idx',ierr)
    allocate(list%n_neighbours(n_elements),stat=ierr)
    call utils_alloc_check('neighbour_list_init','list%n_neighbours',ierr)

    ! Allocate data
    allocate(list%neighbours(list%max_total_neighbours),stat=ierr)
    call utils_alloc_check('neighbour_list_init','list%neighbours',ierr)

  end subroutine neighbour_list_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine neighbour_list_free(list)
    !==========================================================================!
    ! Frees storage associated with a neighbour list object.                   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   list  (in/out):     The neighbour list object that is freed on return. !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, November 2013                                 !
    !==========================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(NL_NEIGHBOUR_LIST), intent(inout) :: list

    ! Local variables
    integer :: ierr

    ! ------------------------------------------------------------------------

    ! Deallocate data
    deallocate(list%neighbours,stat=ierr)
    call utils_dealloc_check('neighbour_list_init','list%neighbours',ierr)

    ! Deallocate index
    deallocate(list%n_neighbours,stat=ierr)
    call utils_dealloc_check('neighbour_list_init','list%n_neighbours',ierr)
    deallocate(list%last_idx,stat=ierr)
    call utils_dealloc_check('neighbour_list_init','list%last_idx',ierr)
    deallocate(list%first_idx,stat=ierr)
    call utils_dealloc_check('neighbour_list_init','list%first_idx',ierr)

    list%n_elements = -1
    list%max_total_neighbours = -1
    list%populated = .false.
    list%list_name = "! STALE !"

  end subroutine neighbour_list_free

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine neighbour_list_populate(list, n_elements, elements, &
       neighbour_counts, neighbours)
    !==========================================================================!
    ! Populates a neighbour list with data passed in the input parameters.     !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   list (inout): The NL to populate, assumed to have been initialised.    !
    !   n_elements (in): The number of elements on this MPI rank               !
    !                    (so, typically, n_my_atoms).                          !
    !   elements (in): Array with the indices of elements on this MPI rank     !
    !                  so, typically, my_atoms)                                !
    !   neighbour_counts (in): Numbers of neighbours on this MPI rank.         !
    !   neighbours (in): Neighbours themselves. 1st index: neighbour,          !
    !                    2nd index: element (atom).                            !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, November 2013                                 !
    ! Bugfixed by Jacek Dziedzic, November 2021 -- could deadlock, see [*].    !
    !==========================================================================!

    use comms, only: pub_root_proc_id, pub_on_root, pub_any_source, &
         comms_send, comms_irecv, comms_bcast, comms_barrier, comms_wait, &
         pub_total_num_procs, pub_my_proc_id
    use constants, only: CRLF, NORMAL
    use rundat, only: pub_output_detail
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_flushed_string_output

    implicit none

    ! Arguments
    type(NL_NEIGHBOUR_LIST), intent(inout) :: list
    integer, intent(in)                    :: n_elements
    integer, intent(in)                    :: elements(:)
    integer, intent(in)                    :: neighbour_counts(:)
    integer, intent(in)                    :: neighbours(:,:)

    ! Local variables
    logical :: comms_needed
    integer :: elem_i, elem_dest_i
    integer :: elem
    integer :: ierr
    integer, allocatable                   :: send_handles(:)
    integer, allocatable                   :: recv_handles(:)
    character(len=*), parameter            :: myself = 'neighbour_list_populate'

    ! ------------------------------------------------------------------------

    call comms_barrier

    call timer_clock(myself,1)

    if(pub_output_detail >= NORMAL) then
       call utils_flushed_string_output(&
            'NL: Populating neighbour list: '//trim(list%list_name)//'... ')
    end if

    ! jd: Take care of the tails that will not be filled
    list%n_neighbours(:) = -1
    list%first_idx(:) = -1
    list%last_idx(:) = -1
    list%neighbours(:) = -1

    ! jd: Elide the comms stage if MPI absent or MPI present, but only 1 rank
    comms_needed = .false.
#ifdef MPI
    if(pub_total_num_procs > 1) then
      comms_needed = .true.
    end if
#endif

    ! jd: [*] Note on comms in code below.
    !     We need to be careful to avoid deadlocks from the textbook scenario,
    !     where a lot of irecvs are posted, and one of them finally blocks
    !     because MPI buffers are exhausted. This would happen here on root,
    !     who posts as many irecvs as there are atoms, and some of them are for
    !     sends that would later come from itself. If any irecv blocks (and on
    !     my MPI implementation this happens after 1024 irecvs), this would
    !     deadlock. To avoid this, we jump through a few hoops to ensure
    !     information transfer between root and root happens not via MPI, but
    !     via local copy mechanisms -- similarly to what we do when comms_needed
    !     is .false., just for individual elements.
    if(comms_needed) then
       allocate(send_handles(n_elements),stat=ierr)
       call utils_alloc_check('neighbour_list_populate','send_handles',ierr)

       ! --- 1. transmit the neighbour counts ---

       ! jd: On root, post receives for neighbour counts of each element
       if(pub_on_root) then
          allocate(recv_handles(list%n_elements),stat=ierr)
          call utils_alloc_check('neighbour_list_populate','recv_handles',ierr)

          do elem_dest_i = 1, list%n_elements
             if(any(elem_dest_i == elements(:))) then
                ! [*] On root just do a local copy. Here, conveniently
                !     elem_i == elem_dest_i, so RHS does not need to find the
                !     corresponding elem_i in elements(:).
                list%n_neighbours(elem_dest_i) = neighbour_counts(elem_dest_i)
             else
                ! [*] Elsewhere just post the receives.
                call comms_irecv(pub_any_source, list%n_neighbours(elem_dest_i),&
                     1, tag = elem_dest_i, handle = recv_handles(elem_dest_i))
             end if
          end do
       end if

       ! jd: Blocking-send the neighbour counts
       do elem_i = 1, n_elements
          elem = elements(elem_i)
          if(.not. pub_on_root) then ! [*]
             call comms_send(pub_root_proc_id, neighbour_counts(elem_i), 1, &
                  tag = elem, return_handle = send_handles(elem_i), &
                  add_to_stack = .false.)
             call comms_wait(send_handles(elem_i))
          end if
       end do
    else
       ! jd: In serial mode comms_irecv and comms_send are no-ops. Transfer the
       !     data manually.
       list%n_neighbours(1:list%n_elements) = neighbour_counts(1:list%n_elements)
    end if

    ! jd: Wait for the receives to complete
    if(pub_on_root) then
       if(comms_needed) then
          do elem_dest_i = 1, list%n_elements
             ! [*] Do not post waits for irecvs that were elided earlier.
             if(.not. any(elem_dest_i == elements(:))) then
                call comms_wait(recv_handles(elem_dest_i))
             end if
          end do
       end if

       ! jd: Fill first_idx and last_idx based on n_neighbours
       do elem_dest_i = 1, list%n_elements
          if(elem_dest_i == 1) then
             list%first_idx(elem_dest_i) = 1
          else
             list%first_idx(elem_dest_i) = list%first_idx(elem_dest_i-1) + &
                  list%n_neighbours(elem_dest_i-1)
          end if
          list%last_idx(elem_dest_i) = &
                list%first_idx(elem_dest_i) + list%n_neighbours(elem_dest_i)-1
       end do
    end if

    ! --- 2. transmit the neighbour lists themselves ---

    if(comms_needed) then
       ! jd: On root, post receives for neighbour lists of each element
       if(pub_on_root) then
          do elem_dest_i = 1, list%n_elements
             if(any(elem_dest_i == elements(:))) then
                ! [*] On root just do a local copy.
                list%neighbours(list%first_idx(elem_dest_i):&
                     list%first_idx(elem_dest_i)+list%n_neighbours(elem_dest_i)-1) = &
                     neighbours(1:neighbour_counts(elem_dest_i),elem_dest_i)
             else
                ! [*] Elsewhere just post the receives.
                call comms_irecv(pub_any_source, &
                     list%neighbours(list%first_idx(elem_dest_i)), &
                     list%n_neighbours(elem_dest_i), &
                     tag = elem_dest_i, handle = recv_handles(elem_dest_i))
             end if
          end do
       end if

       ! jd: Blocking-send the neighbour lists
       do elem_i = 1, n_elements
          elem = elements(elem_i)
          if(.not. pub_on_root) then
             call comms_send(pub_root_proc_id, neighbours(:,elem_i), &
                  neighbour_counts(elem_i), &
                  tag = elem, return_handle = send_handles(elem_i), &
                  add_to_stack = .false.)
             call comms_wait(send_handles(elem_i))
          end if
       end do
    else
       ! jd: In serial mode comms_irecv and comms_send are no-ops. Transfer the
       !     data manually.
       do elem_dest_i = 1, list%n_elements
          list%neighbours(list%first_idx(elem_dest_i):&
               list%first_idx(elem_dest_i)+list%n_neighbours(elem_dest_i)-1) = &
               neighbours(1:neighbour_counts(elem_dest_i),elem_dest_i)
       end do
    end if

    if(comms_needed) then
       ! jd: Wait for the receives to complete
       if(pub_on_root) then
          do elem_dest_i = 1, list%n_elements
             ! [*] Do not post waits for irecvs that were elided earlier.
             if(.not. any(elem_dest_i == elements(:))) then
                call comms_wait(recv_handles(elem_dest_i))
             end if
          end do
       end if

       ! --- 3. broadcast results from root ---
       call comms_bcast(pub_root_proc_id,list%n_neighbours(:),list%n_elements)
       call comms_bcast(pub_root_proc_id,list%first_idx(:),list%n_elements)
       call comms_bcast(pub_root_proc_id,list%last_idx(:),list%n_elements)
       call comms_bcast(pub_root_proc_id,list%neighbours(:), &
            list%max_total_neighbours)

       if(pub_on_root) then
          deallocate(recv_handles,stat=ierr)
          call utils_dealloc_check('neighbour_list_populate','recv_handles',ierr)
       end if

       deallocate(send_handles,stat=ierr)
       call utils_dealloc_check('neighbour_list_populate','send_handles',ierr)

       call comms_barrier
    end if

    list%populated = .true.

    if(pub_output_detail >= NORMAL) then
       call utils_flushed_string_output('done.'//CRLF)
    end if

    call timer_clock(myself,2)

  end subroutine neighbour_list_populate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine neighbour_list_init_from_sparse(list, &
       list_name, sparse_mat, listunit)
    !==========================================================================!
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2013.                              !
    ! Modified to remove pub_par by Joseph Prentice, June 2018                 !
    !==========================================================================!

    use comms, only: pub_my_proc_id, comms_reduce
    use parallel_strategy, only: PARAL_INFO
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index, &
         sparse_get_par
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    type(NL_NEIGHBOUR_LIST), intent(out) :: list
    character(len=*), intent(in)         :: list_name
    type(SPAM3), intent(in)              :: sparse_mat
    integer, optional                    :: listunit

    integer :: n_my_atoms
    integer :: idxlen_a
    integer :: fstnzbridx_cola, lstnzbridx_cola
    integer :: local_a
    integer :: global_a
    integer :: highest_neighbour_count
    integer :: ierr

    integer, allocatable :: idx_a(:)
    integer, allocatable :: my_atoms(:)
    integer, allocatable :: my_neighbour_counts(:)
    integer, allocatable :: my_neighbours(:,:)

    character(len=*), parameter :: myself = 'neighbour_list_init_from_sparse'
    type(PARAL_INFO), pointer :: par

    ! -------------------------------------------------------------------------

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, sparse_mat)

    n_my_atoms = par%num_atoms_on_proc(pub_my_proc_id)

    allocate(my_atoms(n_my_atoms),stat=ierr)
    call utils_alloc_check(myself,'my_atoms',ierr)
    allocate(my_neighbour_counts(n_my_atoms),stat=ierr)
    call utils_alloc_check(myself,'my_neighbour_counts',ierr)

    ! jd: Find global indices of my atoms A, put them in an array
    do local_a = 1, n_my_atoms
       global_a = par%first_atom_on_proc(pub_my_proc_id) + local_a-1
       my_atoms(local_a) = global_a
    end do

    ! jd: Find neighbours (in the sense of sparsity pattern of 'sparse_mat')
    !     of all my atoms A
    idxlen_a = sparse_index_length(sparse_mat)
    allocate(idx_a(idxlen_a), stat=ierr)
    call utils_alloc_check(myself, 'idx_a', ierr)
    call sparse_generate_index(idx_a, sparse_mat)

    ! jd: Determine atoms B that are s-neighbours with each A, count them
    highest_neighbour_count = -1
    do local_a = 1, n_my_atoms
       fstnzbridx_cola = idx_a(local_a)
       lstnzbridx_cola = idx_a(local_a+1)-1
       my_neighbour_counts(local_a) = idx_a(local_a+1) - idx_a(local_a)
       highest_neighbour_count = &
            max(highest_neighbour_count, my_neighbour_counts(local_a))
    end do

    allocate(my_neighbours(highest_neighbour_count,n_my_atoms),stat=ierr)
    call utils_alloc_check(myself,'my_neighbours',ierr)

    do local_a = 1, n_my_atoms
       fstnzbridx_cola = idx_a(local_a)
       lstnzbridx_cola = idx_a(local_a+1)-1
       my_neighbours(1:my_neighbour_counts(local_a),local_a) = &
            idx_a(fstnzbridx_cola:lstnzbridx_cola)
    end do

    ! jd: Each proc can have a different value for highest_neighbour_count,
    !     use the global maximum for allocation.
    call comms_reduce('MAX', highest_neighbour_count)

    call neighbour_list_init(list, list_name, par%nat, &
         highest_neighbour_count, listunit)

    call neighbour_list_populate(list, n_my_atoms, my_atoms, &
         my_neighbour_counts, my_neighbours)

    deallocate(my_neighbours,stat=ierr)
    call utils_dealloc_check(myself,'my_neighbours',ierr)
    deallocate(idx_a, stat=ierr)
    call utils_dealloc_check(myself, 'idx_a', ierr)
    deallocate(my_neighbour_counts,stat=ierr)
    call utils_dealloc_check(myself,'my_neighbour_counts',ierr)
    deallocate(my_atoms,stat=ierr)
    call utils_dealloc_check(myself,'my_atoms',ierr)

  end subroutine neighbour_list_init_from_sparse

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine neighbour_list_for_atom_list(out_atom_list, n_out_atom_list, &
       in_atom_list, n_in_atom_list, nl, global_nat)
    !==========================================================================!
    ! Determines a list of atoms (array of global atom indices) that are       !
    ! neighbours of any atom in a given list of atoms (array of global atom    !
    ! indices) according to a given neighbour list.                            !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   out_atom_list (out): The list of atoms that are neighbours is returned !
    !                        here. This is allocated here, caller is responsi- !
    !                        ble for freeing it later.                         !
    !   n_out_atom_list (out): On exit it will contain the number of elements  !
    !                          in out_atom_list.                               !
    !   in_atom_list (in): The list of atoms whose neighbours we want to find. !
    !   n_in_atom_list (in): The number of elements in in_atom_list.           !
    !   nl (in): The neighbour list which decides who is a neighbour of who.   !
    !   global_nat (in): Global number of atoms. Use size(elements) or par%nat.!
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Storage for out_atom_list is allocated here. Caller is responsible for !
    !   freeing it later.                                                      !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2019.                                   !
    !==========================================================================!

    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    ! jd: Arguments
    integer, intent(out), allocatable   :: out_atom_list(:)
    integer, intent(out)                :: n_out_atom_list
    integer, intent(in)                 :: in_atom_list(:)
    integer, intent(in)                 :: n_in_atom_list
    type(NL_NEIGHBOUR_LIST), intent(in) :: nl
    integer, intent(in)                 :: global_nat

    ! jd: Local variables
    integer :: global_a
    integer :: global_b
    integer :: a_idx
    integer :: b_idx
    integer :: my_a_atom_idx
    logical, allocatable :: a_atoms(:)
    integer :: ierr
    character(len=*), parameter :: myself = 'neighbour_list_for_atom_list'

    ! -----------------------------------------------------------------------

    ! jd: Sanity check on neighbour list
    call utils_assert(nl%populated, myself//': Neighbour list '//&
         trim(nl%list_name)//' not populated.')

    ! --------------------------------------------------------------------------
    ! jd: [1] Have a bitmap of all possible A atoms.
    ! --------------------------------------------------------------------------
    allocate(a_atoms(global_nat),stat=ierr)
    call utils_alloc_check(myself,'a_atoms',ierr)
    a_atoms(:) = .false.

    ! --------------------------------------------------------------------------
    ! jd: [2] Go over all A's and tag them in the bitmap.
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    ! jd: Loop over B's in in_atom_list.
    ! --------------------------------------------------------------------------
    loop_B:                                                                    &
    do b_idx = 1, n_in_atom_list
       global_b = in_atom_list(b_idx)

       ! -----------------------------------------------------------------------
       ! jd: Loop over A's that are neighbours of B
       ! -----------------------------------------------------------------------
       loop_A:                                                                 &
       do a_idx = nl%first_idx(global_b), nl%last_idx(global_b)
          global_a = nl%neighbours(a_idx)
          a_atoms(global_a) = .true.
       end do loop_A
    end do loop_B

    ! --------------------------------------------------------------------------
    ! jd: [3] Now that we know which A atoms are in the bitmap, and how many
    !         there are, we can construct a list.
    ! --------------------------------------------------------------------------
    n_out_atom_list = count(a_atoms)
    allocate(out_atom_list(n_out_atom_list),stat=ierr)
    call utils_alloc_check(myself,'out_atom_list',ierr)

    ! --------------------------------------------------------------------------
    ! jd: Loop over *all* global atoms A
    ! --------------------------------------------------------------------------
    my_a_atom_idx = 1
    loop_A2:                                                                   &
    do global_a = 1, global_nat
       if(a_atoms(global_a)) then
          out_atom_list(my_a_atom_idx) = global_a
          my_a_atom_idx = my_a_atom_idx + 1
       end if
    end do loop_A2

    ! jd: Check we've indeed stored all A atoms from the bitmap
    call utils_assert(my_a_atom_idx-1 == n_out_atom_list, myself//&
         ': Logic error: ', my_a_atom_idx-1, n_out_atom_list, global_nat)

    ! --------------------------------------------------------------------------
    ! jd: [4] Bitmap is no longer needed.
    ! --------------------------------------------------------------------------
    deallocate(a_atoms,stat=ierr)
    call utils_dealloc_check(myself,'a_atoms',ierr)

  end subroutine neighbour_list_for_atom_list

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine neighbour_list_for_local_atoms(out_atom_list, n_out_atom_list, &
       nl, par)
    !==========================================================================!
    ! Returns a list of atoms (array of global atom indices) that are          !
    ! neighbours of any sparse-local atom on this proc, according to a given   !
    ! neighbour list.                                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   out_atom_list (out): The list of atoms that are neighbours is returned !
    !                        here. This is allocated here, caller is responsi- !
    !                        ble for freeing it later.                         !
    !   n_out_atom_list (out): On exit it will contain the number of elements  !
    !                          in out_atom_list.                               !
    !   nl (in): The neighbour list which decides who is a neighbour of who.   !
    !   par (in): PARAL_INFO needed to establish local atoms and par%nat.      !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Storage for out_atom_list is allocated here. Caller is responsible for !
    !   freeing it later.                                                      !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2019.                                   !
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! jd: Arguments
    integer, intent(out), allocatable   :: out_atom_list(:)
    integer, intent(out)                :: n_out_atom_list
    type(NL_NEIGHBOUR_LIST), intent(in) :: nl
    type(PARAL_INFO), intent(in)        :: par

    ! jd: Local variables
    integer, allocatable :: in_atom_list(:)
    integer :: n_in_atom_list
    integer :: local_a, global_a
    integer :: ierr
    character(len=*), parameter :: myself = 'neighbour_list_for_local_atoms'

    ! -------------------------------------------------------------------------

    n_in_atom_list = par%num_atoms_on_proc(pub_my_proc_id)
    allocate(in_atom_list(n_in_atom_list),stat=ierr)
    call utils_alloc_check(myself,'in_atom_list',ierr)

    do local_a = 1, n_in_atom_list
       global_a = par%first_atom_on_proc(pub_my_proc_id) + local_a - 1
       in_atom_list(local_a) = global_a
    end do

    call neighbour_list_for_atom_list(out_atom_list, n_out_atom_list, &
         in_atom_list, n_in_atom_list, nl, par%nat)

    deallocate(in_atom_list,stat=ierr)
    call utils_dealloc_check(myself,'in_atom_list',ierr)

  end subroutine neighbour_list_for_local_atoms

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine neighbour_list_procs_of_neighbours(proc_list, n_proc_list, &
       global_at, nl, ngwf_basis)
    !==========================================================================!
    ! Establishes a list of MPI procs that sparse-own neighbours of a given    !
    ! atom, according to a given neighbour list.                               !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   proc_list (out): List of MPI ranks that sparse-own at least one        !
    !                    neighbour of global_at. The first n_proc_list         !
    !                    elements are set, the rest is -1.                     !
    !   n_proc_list (out): On exit, stores the number of elements in proc_list.!
    !   global_at (in): The global atom index whose neighbours' procs we seek. !
    !   nl (in): The neighbour list that determines who is a neighbour of who. !
    !   ngwf_basis (in): NGWF basis in which the *neighbours* live (sic!).     !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   ngwf_basis refers to the neighbours, *not* the atom global_at.         !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in June 2019.                                  !
    !==========================================================================!

    use comms, only: pub_total_num_procs
    use function_basis, only: FUNC_BASIS

    implicit none

    ! jd: Arguments
    integer, intent(out)                :: proc_list(1:pub_total_num_procs)
    integer, intent(out)                :: n_proc_list
    integer, intent(in)                 :: global_at
    type(NL_NEIGHBOUR_LIST), intent(in) :: nl
    type(FUNC_BASIS), intent(in)        :: ngwf_basis

    ! jd Local variables
    integer :: neigh_idx
    integer :: global_neigh
    integer :: first_ngwf_idx_of_neigh
    integer :: proc_of_neigh

    ! -----------------------------------------------------------------------

    n_proc_list = 0

    proc_list(:) = -1

    do neigh_idx = nl%first_idx(global_at), nl%last_idx(global_at)
       global_neigh = nl%neighbours(neigh_idx)
       first_ngwf_idx_of_neigh = ngwf_basis%first_on_atom(global_neigh)
       proc_of_neigh = ngwf_basis%proc_of_func(first_ngwf_idx_of_neigh)

       if(.not. any(proc_list(1:n_proc_list) == proc_of_neigh)) then
          n_proc_list = n_proc_list + 1
          proc_list(n_proc_list) = proc_of_neigh
       end if

    end do

  end subroutine neighbour_list_procs_of_neighbours

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  elemental logical function neighbour_list_are_neighbours(global_at1, &
       global_at2, nl)
    !==========================================================================!
    ! Ascertains if two atoms are neighbours according to a given neighbour    !
    ! list. Or, strictly speaking, if global_at2 is a neighbour of global_at1. !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   global_at1 (in): Global atom index of the first atom.                  !
    !   global_at2 (in): Global atom index of the second atom.                 !
    !   nl (in): The neighbour list that determines who is a neighbour of who. !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   See the comment [#] at the top of hf_exchange_mod to understand that   !
    !   some neighbouring relations are not commutative.                       !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2019.                                  !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    integer, intent(in)                 :: global_at1
    integer, intent(in)                 :: global_at2
    type(NL_NEIGHBOUR_LIST), intent(in) :: nl

    ! -------------------------------------------------------------------------

    neighbour_list_are_neighbours = any(nl%neighbours(&
         nl%first_idx(global_at1):nl%last_idx(global_at1)) == global_at2)

  end function neighbour_list_are_neighbours

end module neighbour_list
