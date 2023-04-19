! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!-----------------------------------------------------------------------!
!                                                                       !
! Hash table module v1.5.                                               !
!                                                                       !
! Written by Jacek Dziedzic in May 2012                                 !
!                                                                       !
! v1.0 (May 2012): Original version.                                    !
! v1.1 (April 2015): Added cleanup prevention.                          !
! v1.2 (January 2017): Added stashes, dont_add_if_full attribute.       !
!                      5 keys.                                          !
! v1.3 (May 2017): Added iterate.                                       !
! v1.4 (March 2019): Generalised cleanup prevention to overfill_strategy!
!                    Fixed integer wraparound on counters.              !
! v1.5 (March 2019): Axpy, iterate by pointer.                          !
!                                                                       !
!-----------------------------------------------------------------------!
!                                                                       !
! This module implements a hash table for storing batches of            !
! real(kind=DP) values that can be retrieved in near-constant time.     !
! Between 1 and 5 keys (integers) can be used to index the hash table.  !
! The batch sizes can differ across entries in the table.               !
!                                                                       !
!-----------------------------------------------------------------------!
! Glossary:                                                             !
!   key: an integer value.                                              !
!                                                                       !
!   keyset: an ordered set of 1-5 keys, used to index the hash table.   !
!                                                                       !
!   batch AKA element: a 1D array of real(kind=DP) values, a unit of    !
!                      storage in the hash table.                       !
!                                                                       !
!   node AKA record: A stored batch, with associated metadata, such as  !
!                    batch size, corresponding keyset, or flag          !
!                    indicating if in use.                              !
!                    Corresponding type: HT_STORAGE_NODE.               !
!                                                                       !
!   primary storage: An array of nodes, part of the hash table.         !
!                    Corresponding field: slots(:).                     !
!                                                                       !
!   slot: A node that is an entry in primary storage.                   !
!                                                                       !
!   chain link: A node that is not in primary storage, but rather in    !
!               a, presumably short, linked list attached to a slot.    !
!                                                                       !
!   chain: A linked list of chain links.                                !
!                                                                       !
!-----------------------------------------------------------------------!
! Currently supported operations include:                               !
! - init: for initialising the hash table.                              !
! - free: for freeing all storage associated with the hash table.       !
! - purge: removes all elements, but keeps the hash table.              !
! - add: for adding a batch of real(kind=DP) values to the hash table.  !
! - axpy: for incrementing a batch of real(kind=DP) values already      !
!         present in the hash table, or adding them if they are absent. !
! - lookup: for retrieving a batch matching a keyset from the table.    !
!           (performs a deep-copy into a supplied buffer).              !
! - lookup_nocount: like lookup, but does not count hits and misses and !
!                   hence allows an intent(in) table.                   !
! - lookup_ptr: for retrieving a batch matching a keyset from the table.!
!               (returns a pointer, does not deep-copy).                !
! - lookup_ptr_nocount.                                                 !
! - probe: for checking if the hash table contains a batch              !
!          corresponding to a given keyset.                             !
! - probe_nocount.                                                      !
! - remove: for removing a batch matching a keyset from the hash table. !
! - list: for showing the statistics about the hash table (memory load, !
!         collision statistics, hit ratio, etc) and, optionally, the    !
!         headers of the stored records.                                !
! - cleanup_at_will: for performing a deferred trimming of the hash     !
!                    table back to maximum size. Deferring cleanup is   !
!                    useful in OMP environments.                        !
! - iterate: for iterating over the hash table and examining the keys   !
!            and, optionally, data.                                     !
! - stash_prepare: initialises a stash associated with a hash table.    !
! - stash: adds data to a stash.                                        !
! - stash_commit: commits data from a stash to associated hash table.   !
! - stash_free: frees a stash.                                          !
!                                                                       !
! Elements can be overwritten (by adding a batch with the same keyset   !
! and setting overwrite = .true.).                                      !
!                                                                       !
! The maximum number of slots is supplied during the initialisation of  !
! the hash table. Approximately 120 bytes are needed to store the       !
! metadata for every slot. Elements added to the hash table fill the    !
! slots, with storage for the actual data allocated when the elements   !
! are added (as the size of the storage depends on the size of the      !
! added elements (batches)). Occasionally, when two keysets yield the   !
! same hash and a collision occurs, the newly added element cannot      !
! be added to the slot table ('primary storage') and will be stored in  !
! chains (short linked lists) attached to the original slot. In the     !
! absence of collisions the addition and lookup times are constant.     !
! Unless the size of the primary storage is ridiculously small compared !
! to the number of actually stored elements, the chance of collisions   !
! is very small and the vast majority of elements are stored in slots,  !
! rather than in chains, leading to a near-constant addition and lookup !
! times. As the cost of empty slots in primary storage is very small,   !
! large hash tables can be constructed cheaply, leading to extremely    !
! good performance.                                                     !
!                                                                       !
! The maximum number of elements the hash table can ever hold is also   !
! supplied during initialisation. This refers to the combined number of !
! elements in primary storage *and* in chains. The hash table is        !
! guaranteed *never* [1] to hold more than this number of elements, how-!
! ever there is no storage cost to any elements before they are added,  !
! and the small node-overhead cost is only paid for elements in primary !
! storage. What happens once the hash table becomes full (the total     !
! number of elements reaches the maximum number of elements), is        !
! determined by 'overfill_strategy'. There are several modes:           !
!                                                                       !
! 'F' - allow overfilling (essentially ignoring max_elements).          !
! 'N' - do not add the element (ignore the add request).                !
! 'B' - break (abort program).                                          !
! 'C' - cleanup (some elements will be removed before the current one   !
!       is added). This means that a certain fraction of elements will  !
! be automatically discarded during the next add operation.             !
! This is known as cleanup. By default 10% of the contents of the       !
! hash table will be discarded, this fraction can be changed by setting !
! the cleanup_fraction field. Cleanup can be ensured never to happen by !
! always keeping the number of elements below the maximum. During       !
! cleanup, the elements to be discarded first are the ones occupying the!
! slots, and the first chain links are moved into the slots. This is an !
! approximation to the "Least-recently added" strategy. Cleanup can be  !
! temporarily prevented by setting prevent_cleanup in hash_table_add(). !
! This is useful in OMP contexts, as cleanups invalidate pointers to    !
! hash table data. [1] In this scenario the hash table *can* take up    !
! more memory that is allowed to -- until the next call to              !
! hash_table_cleanup_at_will().                                         !
!                                                                       !
! Stashes allow temporarily adding new records to a separate object     !
! (stash) that is merged with the hash table later. This is useful in   !
! OMP environments, because it lifts the requirement of using a critical!
! section for read accesses (lookups) to the hash table while simulta-  !
! neously storing new elements. This is best understood by considering  !
! two scenarios:                                                        !
!                                                                       !
! 1) Not using stashes. A single hash table is accessed from multiple   !
!    threads, using lookup and add operations. Add operations, must of  !
!    course be protected with a critical section to avoid race condi-   !
!    tions where two threads update the same data structure simulta-    !
!    neously. Lookups have to be protected with a critical section too, !
!    because if thread B updates the HT during thread A's lookup,       !
!    thread A will see an inconsistent data structure. This is subopti- !
!    mal, because it serialises read accesses. Finally, even with criti-!
!    cal sections around lookups and adds care must be taken when adds  !
!    can potentially lead to cleanups. Following a cleanup all pointers !
!    to data in a HT are invalidated. Thus any data obtained with       !
!    lookup_by_ptr is invalid. This can be worked around using          !
!    prevent_cleanup or dont_add_if_full.                               !
!                                                                       !
! 2) Using stashes. All add operations are replaced with stash opera-   !
!    tions. A stash is thread-local storage. Thus any new elements wind !
!    up in the stash, and no critical sections are needed in the main   !
!    loop. The critical around lookups is not needed, because the HT is !
!    not modified in the loop. The critical around stashes is not needed!
!    since each thread has a separate stash. Pointers are not invalida- !
!    ted, since there are no cleanups. After the OMP loop is complete,  !
!    each thread can merge its stash into the table (from within a      !
!    critical section) -- this is called a commit, and then free the    !
!    stash. There are two drawbacks. First, the data added to the stash !
!    is not accessible until the stash is committed, so will normally   !
!    yield hits only the next time around. This is typically not a prob-!
!    lem -- e.g. in a loop over PPDs we won't be looking at the same PPD!
!    twice in one loop execution, only the next execution will benefit  !
!    from the cached data. Second, the stashes have to be dimensioned   !
!    suitably, and the cost is paid per thread. Stash sizes are not     !
!    reported when listing hash tables.                                 !
!-----------------------------------------------------------------------!
! Planned additions include:                                            !
! - elision of storage for batches that contain only zeroes.            !
! - using a separate datastructure for metadata (linked lists of chain  !
!   links) and dropping pointers, replacing them with integer indices.  !
!   Using a pool instead of allocating every node each time. This will  !
!   improve 'seek' times when chains are long and avoid the use of      !
!   Fortran pointers, allowing dropping 'target' in most instances.     !
!   Difficulties include pool fragmentation and dealing with it.        !
! - strategies other than FIFO for cleaning up the table. Counting      !
!   accesses per node should not cost much, and never-hit nodes could   !
!   be eliminated first.                                                !
! - node pinning.                                                       !
! - serialisation of stashes to files, deserialisation.                 !
! - conversion of a HT to a stash, MPI ops for stashes, akin to packed  !
!   NGWFs.                                                              !
! - separation of counters from main HT objects and storing them in     !
!   a module-wide library. This will enable getting rid of _nocount()   !
!   versions of operations, and all trivial state would be decoupled    !
!   from objects, which could be passed as intent(in) then.             !
! - better hashing functions.                                           !
!-----------------------------------------------------------------------!

module hash_table

  use constants, only: DP, LONG

  implicit none
  private

  ! Public types
  public :: HT_HASH_TABLE
  public :: HT_STORAGE_NODE

  ! Public subroutines
  public :: hash_table_init
  public :: hash_table_free
  public :: hash_table_purge
  public :: hash_table_add
  public :: hash_table_axpy
  public :: hash_table_lookup
  public :: hash_table_lookup_nocount
  public :: hash_table_lookup_ptr
  public :: hash_table_lookup_ptr_nocount
  public :: hash_table_probe
  public :: hash_table_probe_nocount
  public :: hash_table_remove
  public :: hash_table_cleanup_at_will
  public :: hash_table_list
  public :: hash_table_iterate
  public :: hash_table_stash_prepare
  public :: hash_table_stash
  public :: hash_table_stash_commit
  public :: hash_table_stash_free
  public :: hash_table_self_test
  public :: hash_table_calc_cache_size
  public :: hash_table_mem_usage

  ! Public type definitions
  type HT_STORAGE_NODE
     integer(kind=LONG)             :: n_hits
     logical                        :: in_use
     integer                        :: key1, key2, key3, key4, key5
     integer                        :: batch_size
     real(kind=DP), allocatable     :: node_data(:)
     logical                        :: contains_only_zeroes
     type(HT_STORAGE_NODE), pointer :: next
  end type HT_STORAGE_NODE

  type HT_HASH_TABLE
     integer :: n_keys = -1       ! number of keys
     integer :: n_slots = 0       ! number of occupied slots in slots(:)
     integer :: max_slots = 0     ! size of the slots(:) array
     integer :: max_elements = 0  ! maximum number of nodes in total
     integer :: listunit = 0      ! unit for output from hash_table_list
     integer :: n_chained = 0     ! number of nodes stored in chains
     integer :: longest_chain = 0 ! length of the longest chain
     integer(kind=LONG) :: n_hits = 0           ! total number of hits so far
     integer(kind=LONG) :: n_misses = 0         ! total number of misses so far
     integer(kind=LONG) :: n_adds_ok = 0        ! total number of records added
     integer(kind=LONG) :: n_adds_when_full = 0 ! total n. recs added, overfilling
     integer(kind=LONG) :: n_adds_rejected = 0  ! total n. recs not added b/c table full
     integer(kind=LONG) :: n_removed_recs = 0   ! total number of removed records
     integer(kind=LONG) :: total_slot_data_size = 0  ! sum of all batch_size in slots
     integer(kind=LONG) :: total_chain_data_size = 0 ! sum of all batch_size in chains
     integer :: n_cleanups = 0          ! total number of cleanups so far
     integer :: n_purges = 0            ! total number of purges so far
     integer :: first_slot_for_cleanup = 1 ! cleanup starts from this slot
     character(len=16) :: table_name       ! user-defined name
     real(kind=DP) :: cleanup_fraction = 0.1_DP ! fraction of max_elements to
                                                ! remove when table full.
     type(HT_STORAGE_NODE), allocatable :: slots(:) ! primary storage
  end type HT_HASH_TABLE

  integer, parameter, public :: node_size = 128 ! Educated guesstimate from measurement

  integer, parameter :: C_STASH_HEADER_SIZE = 4
  integer, parameter :: C_STASH_OVERHEAD = 6 ! five keys + ndata

  ! ---------------------------------------------------------------------------

contains

  ! ---------------------------------------------------------------------------
  ! ----- P U B L I C   S U B R O U T I N E S   A N D   F U N C T I O N S -----
  ! ---------------------------------------------------------------------------

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_init(table, table_name, n_keys, max_slots, &
       max_elements, listunit)
    !==========================================================================!
    ! Initialises a fresh hash table object.                                   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table (out):    The hash table object that is initialised on return.   !
    !   table_name (in): User-defined name associated with the table.
    !   n_keys (in):    The number of keys that'll be used to index the table. !
    !                   Currently between 1 and 4 keys are supported.          !
    !   max_slots (in): The maximum number of entries in primary storage, or   !
    !                   the total number of hashes that can be stored. Note    !
    !                   that the hash table CAN store more elements than there !
    !                   are slots, but performance will slowly degrade as more !
    !                   elements are stored in chains.                         !
    !   max_elements (in): Maximum number of entries in all storage. The table !
    !                      will undergo cleanup if an attempt to add more      !
    !                      elements will be made.                              !
    !   listunit (in, opt): An open unit to where the output of hash_table_list!
    !                       is to be directed. Defaults to stdout if omitted.  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, May 2012                                      !
    !==========================================================================!

    use constants, only: stdout
    use utils, only: utils_alloc_check, utils_assert

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(out) :: table
    character(len=*)                 :: table_name
    integer, intent(in)              :: n_keys
    integer, intent(in)              :: max_slots
    integer, intent(in)              :: max_elements
    integer, intent(in), optional    :: listunit

    ! Local variables
    integer :: ierr
    integer :: slot
    character(len=*), parameter :: myself = 'hash_table_init'

    ! ------------------------------------------------------------------------

    call utils_assert(n_keys>=1 .and. n_keys<=5, myself//': '//&
         trim(table_name)//': Only hash tables with 1 to 5 keys are supported.')
    call utils_assert(max_slots>=0, myself//': '//trim(table_name)//&
         ': Cannot have hash tables with fewer than 0 slots.')
    call utils_assert(max_elements >= max_slots, myself//': '//&
         trim(table_name)//': Cannot have hash tables with max_elements &
         &smaller than max_slots.', max_elements, max_slots)

    table%max_slots = max_slots
    table%max_elements = max_elements
    table%n_keys = n_keys
    table%table_name = table_name
    if(present(listunit)) then
       table%listunit = listunit
    else
       table%listunit = stdout
    end if

    ! Corner case of dummy, empty hash_table
    if(max_slots == 0) return

    ! Allocate primary storage (metadata only)
    allocate(table%slots(max_slots),stat=ierr)
    call utils_alloc_check('hash_table_init','table%slots',ierr)

    ! Hash table is initially empty -- no slots or chains are in use
    do slot = 1, max_slots
       table%slots(slot)%n_hits = 0
       table%slots(slot)%in_use = .false.
       nullify(table%slots(slot)%next)
    end do

    table%total_slot_data_size = 0
    table%total_chain_data_size = 0

  end subroutine hash_table_init
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_free(table)
    !==========================================================================!
    ! Frees all storage associated with a hash table. The data contained in    !
    ! the table is freed, so is metadata describing primary storage. Chains    !
    ! are truncated to zero. The resulting hash table remains valid, but is    !
    ! tagged by setting n_keys to -1, so that it cannot be meaningfully        !
    ! accessed any more, unless it is initialised anew.                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table (inout): The hash table object that is to be freed.              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, May 2012                                      !
    !==========================================================================!

    use timer, only: timer_clock
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(inout) :: table

    ! Local variables
    character(len=*), parameter    :: myself = 'hash_table_free'
    integer :: ierr

    ! ------------------------------------------------------------------------

    call timer_clock(myself,1)

    ! Delete slots and chains
    call hash_table_purge(table)

    ! Deallocate the slot table itself
    if(table%max_slots /= 0) then
       deallocate(table%slots,stat=ierr)
       call utils_dealloc_check('hash_table_free','table%slots',ierr)
    end if

    ! Tag the hash table so that it's not accidentally accessed
    table%n_keys = -1
    table%max_slots = 0

    call timer_clock(myself,2)

  end subroutine hash_table_free
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_purge(table)
    !==========================================================================!
    ! Removes all elements of a hash_table, but retains the table itself.      !
    ! Hit/miss counters are preserved.                                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table (inout): The hash table object that is to be purged.             !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, October 2016.                                 !
    !==========================================================================!

    use constants, only: garbage_real, CRLF
    use rundat, only: pub_debug_on_root, pub_debug
    use timer, only: timer_clock
    use utils, only: utils_dealloc_check, utils_assert, utils_abort, &
         utils_flushed_string_output

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(inout) :: table

    ! Local variables
    character(len=*), parameter    :: myself = 'hash_table_purge'
    type(HT_STORAGE_NODE), pointer :: ptr, next_ptr
    integer :: slot
    integer :: ierr

    ! ------------------------------------------------------------------------

    table%n_purges = table%n_purges + 1

    ! Corner case of dummy, empty hash_table
    if(table%max_slots == 0) return

    call timer_clock(myself//'_['//trim(table%table_name)//']',1)

    if(pub_debug_on_root) then
       call utils_flushed_string_output('HT: Purging '//&
            trim(table%table_name)//'.'//CRLF,.true.)
    end if

    call utils_assert(table%n_keys /= -1, &
         myself//': Attempt to purge a hash table that has been freed already.')

    ! Clear all used slots, deallocate their chains. Do not deallocate slots.
    do slot = 1, table%max_slots
       if(table%slots(slot)%in_use) then

          table%slots(slot)%n_hits = 0

          if(pub_debug) then
             table%slots(slot)%node_data(:) = garbage_real
             table%slots(slot)%key1 = -1
             table%slots(slot)%key2 = -1
             table%slots(slot)%key3 = -1
             table%slots(slot)%key4 = -1
             table%slots(slot)%key5 = -1
          end if

          ! Remove storage associated with slot's data
          if(.not. table%slots(slot)%contains_only_zeroes) then
             deallocate(table%slots(slot)%node_data,stat=ierr)
             call utils_dealloc_check(myself,'node_data', ierr)
          else
             call utils_abort(myself//': contains_only_zeroes not yet supported.')
          end if

          ! Remove the chain
          ptr => table%slots(slot)%next
          do while (associated(ptr))
             next_ptr => ptr%next
             ! Deallocate data storage of this node
             deallocate(ptr%node_data,stat=ierr)
             call utils_dealloc_check(myself,'node_data', ierr)
             ! Deallocate the node itself
             deallocate(ptr,stat=ierr)
             call utils_dealloc_check(myself,'node', ierr)
             ! Update stats
             table%n_chained = table%n_chained - 1
             ! Continue to the next chain link
             ptr => next_ptr
          end do

          ! Mark as not in use
          table%slots(slot)%in_use = .false.

          ! Update stats
          table%n_slots = table%n_slots - 1

       end if
    end do

    table%first_slot_for_cleanup = 1
    table%total_slot_data_size = 0
    table%total_chain_data_size = 0

    call timer_clock(myself//'_['//trim(table%table_name)//']',2)

  end subroutine hash_table_purge
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_add(table, data_to_add, ndata, &
       key1, key2, key3, key4, key5, &
       overwrite, overfill_strategy)
    !==========================================================================!
    ! Adds data (a batch of real(kind=DP) values) to the hash table.           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table (inout): The hash table object where the data is to be added.    !
    !   data_to_add (in): The data to be added.                                !
    !   ndata (in): The number of values in the batch to be added.             !
    !   key1..key5 (in): The keys (constituting a 'keyset') which will be used !
    !                    to retrieve the data. Key1 is mandatory, further keys !
    !                    are optional. The number of supplied keys must match  !
    !                    the number of keys declared when initialising the     !
    !                    hash table.                                           !
    !   overwrite (inout, optional): If omitted or .false., an attempt to add  !
    !                             data to a keyset already present in the      !
    !                             table will cause an abort. If .true., the    !
    !                             new data will replace the old data. If the   !
    !                             new batch size is different, storage will be !
    !                             automatically re-allocated to accommodate    !
    !                             the new data. In case of an overwrite, it    !
    !                             will be .true. on exit. If there was no      !
    !                             overwrite, it will be .false. on exit.       !
    !   overfill_strategy (in, opt): What to do if the HT is full already:     !
    !      'F' - allow overfilling (ignore max_elements).                      !
    !      'C' - cleanup (remove some elements before adding this one).        !
    !      'N' - do not add the element (ignore request).                      !
    !      'B' - break (abort).                                                !
    !      Default (if omitted) is 'C'.                                        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, May 2012                                      !
    ! Extended with prevent_cleanup by Jacek Dziedzic, April 2015.             !
    ! Extended with dont_add_if_full by Jacek Dziedzic, Jan 2017.              !
    ! Generalised to overfill_strategy by Jacek Dziedzic, March 2019.          !
    !==========================================================================!

    use rundat, only: pub_debug
    use utils, only: utils_isnan, utils_abort

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(inout) :: table
    real(kind=DP), intent(in)          :: data_to_add(:)
    integer, intent(in)                :: ndata
    integer, intent(in)                :: key1
    integer, intent(in), optional      :: key2
    integer, intent(in), optional      :: key3
    integer, intent(in), optional      :: key4
    integer, intent(in), optional      :: key5
    logical, intent(inout), optional   :: overwrite
    character, intent(in), optional    :: overfill_strategy

    ! Local variables
    integer   :: local_key1, local_key2, local_key3, local_key4, local_key5
    integer   :: hash
    integer   :: i
    character :: loc_overfill_strategy

    ! ------------------------------------------------------------------------

    ! NB. This subroutine is called from OMP contexts. Do not utils_trace() it.

    ! Corner case of dummy, empty hash_table
    if(table%max_slots == 0) return

    ! Sanity check keys, put in zeroes where optional args are not supplied
    call hash_table_process_keys(table, &
         local_key1, local_key2, local_key3, local_key4, local_key5, & ! out
         key1, key2, key3, key4, key5)                                 ! in

    ! Generate hash
    hash = hash_table_generate_hash(&
         local_key1, local_key2, local_key3, local_key4, local_key5, &
         table%max_slots)

    if(present(overfill_strategy)) then
       loc_overfill_strategy = overfill_strategy
    else
       loc_overfill_strategy = 'C'
    end if

    ! Add to table
    call hash_table_add_data(table, data_to_add, ndata, hash, &
         local_key1, local_key2, local_key3, local_key4, local_key5, &
         overwrite, loc_overfill_strategy)

    if(pub_debug) then
       do i = 1, ndata
          if(utils_isnan(data_to_add(i))) then
             call utils_abort('Sanity check failed. NaN detected when adding &
             &data to '//trim(table%table_name),i,ndata);
          end if
       end do
    end if

  end subroutine hash_table_add
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_axpy(table, x, ndata, key1, key2, key3, key4, key5, a, &
       overfill_strategy, dont_axpy_first_nelems)
    !==========================================================================!
    ! Performs y := y + a*x on a record in a hash table.                       !
    ! In the event the element is not found in the table, it is added afresh.  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table (inout): The hash table object where the data is to be axpied in.!
    !   x (in): The data to be added.                                          !
    !   ndata (in): The number of values in the batch to be added.             !
    !   key1..key5 (in): The keys (constituting a 'keyset') which will be used !
    !                    to retrieve the data. Key1 is mandatory, further keys !
    !                    are optional. The number of supplied keys must match  !
    !                    the number of keys declared when initialising the     !
    !                    hash table.                                           !
    !   a (in, opt): The scalar for optionally multiplying data in x.          !
    !   overfill_strategy (in): See documentation of hash_table_add().         !
    !   dont_axpy_first_nelems (in, opt): If given, this many elements at the  !
    !                                     beginning of 'x' will be stored, not !
    !                                     axpied. Useful for metadata.         !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, March 2019.                                   !
    !==========================================================================!

    use rundat, only: pub_debug
    use utils, only: utils_isnan, utils_alloc_check, utils_dealloc_check, &
         utils_abort, utils_assert

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(inout), target :: table
    real(kind=DP), intent(in)                  :: x(:)
    integer, intent(in)                        :: ndata
    integer, intent(in)                        :: key1
    integer, intent(in), optional              :: key2
    integer, intent(in), optional              :: key3
    integer, intent(in), optional              :: key4
    integer, intent(in), optional              :: key5
    real(kind=DP), intent(in), optional        :: a
    character, intent(in)                      :: overfill_strategy
    integer, intent(in), optional              :: dont_axpy_first_nelems

    ! Local variables
    real(kind=DP), pointer :: y(:)
    integer       :: local_key1, local_key2, local_key3, local_key4, local_key5
    integer       :: hash
    integer       :: i
    integer       :: ysize
    integer       :: loc_dont_axpy_first_nelems
    real(kind=DP) :: loc_a
    integer       :: ierr
    logical       :: overwrite
    character(len=*), parameter :: myself = 'hash_table_axpy'

    ! ------------------------------------------------------------------------

    ! Corner case of dummy, empty hash_table
    if(table%max_slots == 0) return

    ! Sanity check keys, put in zeroes where optional args are not supplied
    call hash_table_process_keys(table, &
         local_key1, local_key2, local_key3, local_key4, local_key5, & ! out
         key1, key2, key3, key4, key5)                                 ! in

    ! Generate hash
    hash = hash_table_generate_hash(&
         local_key1, local_key2, local_key3, local_key4, local_key5, &
         table%max_slots)

    if(present(a)) then
       loc_a = a
    else
       loc_a = 1.0_DP
    end if

    if(present(dont_axpy_first_nelems)) then
       loc_dont_axpy_first_nelems = dont_axpy_first_nelems
    else
       loc_dont_axpy_first_nelems = 0
    end if

    ! Retrieve record from table
    call hash_table_lookup_ptr_nocount(y,ysize,table,key1,key2,key3,key4,key5)
    if(ysize == -1) then
       ! jd: If not present, imitate a record of zeroes
       allocate(y(ndata),stat=ierr)
       call utils_alloc_check(myself,'y',ierr)
       y(:) = 0.0_DP
    else
       if(ndata /= ysize) then
          call utils_assert(.false., myself//': Data sizes incompatible', &
               ndata, ysize)
       end if
    end if

    ! Elements to simply store, if any
    y(1:loc_dont_axpy_first_nelems) = x(1:loc_dont_axpy_first_nelems)

    ! Elements to axpy (usually all)
    y(1+loc_dont_axpy_first_nelems:ndata) = &
         y(1+loc_dont_axpy_first_nelems:ndata) + &
         loc_a * x(1+loc_dont_axpy_first_nelems:ndata)

    ! Add to table
    overwrite = .true.
    call hash_table_add_data(table, y, ndata, hash, &
         local_key1, local_key2, local_key3, local_key4, local_key5, &
         overwrite, overfill_strategy)

    if(pub_debug) then
       do i = 1, ndata
          if(utils_isnan(y(i))) then
             call utils_abort(myself//': Sanity check failed. NaN detected &
                  &axpying a record in '//trim(table%table_name),i,ndata)
          end if
       end do
    end if

    if(ysize == -1) then
       ! jd: Clean up our potential imitation of a record
       deallocate(y,stat=ierr)
       call utils_dealloc_check(myself,'y',ierr)
    end if

  end subroutine hash_table_axpy
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  integer function hash_table_probe(table, key1, key2, key3, key4, key5)
    !==========================================================================!
    ! Returns the size of the batch stored for the given keyset, or -1 if the  !
    ! hash table does not contain an element corresponding to the keyset.      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table (inout): The hash table object. Declared inout, so that hits and !
    !                  misses can be counted and stored inside the object.     !
    !   key1..key5 (in): The keys (constituting a 'keyset') which correspond   !
    !                    to the sought data. Key1 is mandatory, further keys   !
    !                    are optional. The number of supplied keys must match  !
    !                    the number of keys declared when initialising the     !
    !                    hash table.                                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, May 2012                                      !
    !==========================================================================!

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(inout), target :: table
    integer, intent(in)                        :: key1
    integer, intent(in), optional              :: key2
    integer, intent(in), optional              :: key3
    integer, intent(in), optional              :: key4
    integer, intent(in), optional              :: key5

    ! Local variables
    real(kind=DP) :: dummy_result_never_referenced(1)
    integer       :: retval

    ! ------------------------------------------------------------------------

    call hash_table_lookup(dummy_result_never_referenced, &
         retval, table, key1, key2, key3, key4, key5, want_result = .false.)

    hash_table_probe = retval

  end function hash_table_probe
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  integer function hash_table_probe_nocount(table, key1, key2, key3, key4, key5)
    !==========================================================================!
    ! Returns the size of the batch stored for the given keyset, or -1 if the  !
    ! hash table does not contain an element corresponding to the keyset.      !
    ! Does not count hits and misses, allowing for the table to be intent(in)  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table (in): The hash table object.                                     !
    !   key1..key5 (in): The keys (constituting a 'keyset') which correspond   !
    !                    to the sought data. Key1 is mandatory, further keys   !
    !                    are optional. The number of supplied keys must match  !
    !                    the number of keys declared when initialising the     !
    !                    hash table.                                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, October 2016                                  !
    !==========================================================================!

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(in), target :: table
    integer, intent(in)                     :: key1
    integer, intent(in), optional           :: key2
    integer, intent(in), optional           :: key3
    integer, intent(in), optional           :: key4
    integer, intent(in), optional           :: key5

    ! Local variables
    real(kind=DP) :: dummy_result_never_referenced(1)
    integer       :: retval

    ! ------------------------------------------------------------------------

    call hash_table_lookup_nocount(dummy_result_never_referenced, &
         retval, table, key1, key2, key3, key4, key5, want_result = .false.)

    hash_table_probe_nocount = retval

  end function hash_table_probe_nocount
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_lookup(lookup_result, batch_size, table, &
       key1, key2, key3, key4, key5, want_result)
    !==========================================================================!
    ! Retrieves the data (batch of real(kind=DP)'s) corresponding to the given !
    ! keyset. A deep copy from the hash table to lookup_result is performed.   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   lookup_result (out): The buffer where the data will be retrieved to.   !
    !                        The caller is responsible for supplying a buffer  !
    !                        that is big enough to hold the data.              !
    !                        hash_table_probe() can be used to query for the   !
    !                        size of the stored data. If want_result is .false.!
    !                        this buffer will NOT BE REFERENCED, thus giving   !
    !                        another way of querying for the data length or    !
    !                        just for the presence of the data in the table.   !
    !   batch_size (out): Will contain, on return, the number of values in the !
    !                     returned batch. The argument 'lookup_result' must be !
    !                     at least that big, unless want_result is .false.     !
    !   table (inout): The hash table object. Declared inout, so that hits and !
    !                  misses can be counted and stored inside the object.     !
    !   key1..key5 (in): The keys (constituting a 'keyset') which correspond   !
    !                    to the sought data. Key1 is mandatory, further keys   !
    !                    are optional. The number of supplied keys must match  !
    !                    the number of keys declared when initialising the     !
    !                    hash table.                                           !
    !   want_result (in, optional): If .false. is passed, lookup_result will   !
    !                               not be written to or even referenced and   !
    !                               only batch_size will be set on return.     !
    !                               Defaults to .true. when omitted.           !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Only the first batch_size elements in lookup_result will be written to.!
    !   This is done for the sake of performance. As lookup_result is an       !
    !   intent(out) parameter, this may trigger debug-mode checks with         !
    !   overzealous compilers -- technically it's our responsibility to write  !
    !   to the entire array.                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, May 2012                                      !
    !==========================================================================!

    implicit none

    ! Arguments
    real(kind=DP), intent(out)                 :: lookup_result(:)
    integer, intent(out)                       :: batch_size
    type(HT_HASH_TABLE), intent(inout), target :: table
    integer, intent(in)                        :: key1
    integer, intent(in), optional              :: key2
    integer, intent(in), optional              :: key3
    integer, intent(in), optional              :: key4
    integer, intent(in), optional              :: key5
    logical, intent(in), optional              :: want_result

    ! Local variables
    logical                         :: local_want_result
    type(HT_STORAGE_NODE), pointer  :: ptr, prev_ptr

    ! ------------------------------------------------------------------------
    ! jd: Do not time this subroutine, you will get overheads.

    ! Corner case of dummy, empty hash_table
    if(table%max_slots == 0) then
       batch_size = -1
       return
    end if

    if(present(want_result)) then
       local_want_result = want_result
    else
       local_want_result = .true.
    end if

    ! Locate the record
    call hash_table_locate(ptr, prev_ptr, table, key1, key2, key3, key4, key5)

    if(associated(ptr)) then
       ! Record found, return data
       batch_size = ptr%batch_size
       if(local_want_result) then
          lookup_result(1:batch_size) = ptr%node_data(1:batch_size)
       end if
    else
       ! No such record, notify caller
       batch_size = -1
    end if

  end subroutine hash_table_lookup
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_lookup_nocount(lookup_result, batch_size, table, &
       key1, key2, key3, key4, key5, want_result)
    !==========================================================================!
    ! Retrieves the data (batch of real(kind=DP)'s) corresponding to the given !
    ! keyset. A deep copy from the hash table to lookup_result is performed.   !
    ! Hit/miss counters in the hash table are not updated, which allows it to  !
    ! be intent(in).                                                           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   lookup_result (out): The buffer where the data will be retrieved to.   !
    !                        The caller is responsible for supplying a buffer  !
    !                        that is big enough to hold the data.              !
    !                        hash_table_probe() can be used to query for the   !
    !                        size of the stored data. If want_result is .false.!
    !                        this buffer will NOT BE REFERENCED, thus giving   !
    !                        another way of querying for the data length or    !
    !                        just for the presence of the data in the table.   !
    !   batch_size (out): Will contain, on return, the number of values in the !
    !                     returned batch. The argument 'lookup_result' must be !
    !                     at least that big, unless want_result is .false.     !
    !   table (in): The hash table object.                                     !
    !   key1..key5 (in): The keys (constituting a 'keyset') which correspond   !
    !                    to the sought data. Key1 is mandatory, further keys   !
    !                    are optional. The number of supplied keys must match  !
    !                    the number of keys declared when initialising the     !
    !                    hash table.                                           !
    !   want_result (in, optional): If .false. is passed, lookup_result will   !
    !                               not be written to or even referenced and   !
    !                               only batch_size will be set on return.     !
    !                               Defaults to .true. when omitted.           !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Only the first batch_size elements in lookup_result will be written to.!
    !   This is done for the sake of performance. As lookup_result is an       !
    !   intent(out) parameter, this may trigger debug-mode checks with         !
    !   overzealous compilers -- technically it's our responsibility to write  !
    !   to the entire array.                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, October 2016                                  !
    !==========================================================================!

    implicit none

    ! Arguments
    real(kind=DP), intent(out)                 :: lookup_result(:)
    integer, intent(out)                       :: batch_size
    type(HT_HASH_TABLE), intent(in), target    :: table
    integer, intent(in)                        :: key1
    integer, intent(in), optional              :: key2
    integer, intent(in), optional              :: key3
    integer, intent(in), optional              :: key4
    integer, intent(in), optional              :: key5
    logical, intent(in), optional              :: want_result

    ! Local variables
    logical                         :: local_want_result
    type(HT_STORAGE_NODE), pointer  :: ptr, prev_ptr

    ! ------------------------------------------------------------------------
    ! jd: Do not time this subroutine, you will get overheads.

    ! Corner case of dummy, empty hash_table
    if(table%max_slots == 0) then
       batch_size = -1
       return
    end if

    if(present(want_result)) then
       local_want_result = want_result
    else
       local_want_result = .true.
    end if

    ! Locate the record
    call hash_table_locate_nocount(ptr, prev_ptr, table, &
         key1, key2, key3, key4, key5)

    if(associated(ptr)) then
       ! Record found, return data
       batch_size = ptr%batch_size
       if(local_want_result) then
          lookup_result(1:batch_size) = ptr%node_data(1:batch_size)
       end if
    else
       ! No such record, notify caller
       batch_size = -1
    end if

  end subroutine hash_table_lookup_nocount
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_lookup_ptr(ret_ptr, batch_size, table, &
       key1, key2, key3, key4, key5)
    !==========================================================================!
    ! Retrieves the data (batch of real(kind=DP)'s) corresponding to the given !
    ! keyset. The data is returned by pointer (ret_ptr), to avoid deep-copying !
    ! the contents. The size of the batch is returned in batch_size. If the    !
    ! element is absent in the hash table, a batch_size of -1 is returned and  !
    ! ret_ptr remains untouched.                                               !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   ret_ptr (inout): (Pointer to) the buffer where the data resides.       !
    !                    Intent(inout) is used, because this is not touched    !
    !                    when the element is not present in the hash table.    !
    !                    The pointer remains valid ONLY until the next call to !
    !                    hash_table_remove(), hash_table_cleanup(),            !
    !                    hash_table_free() or hash_table_add(overwrite=.true.).!
    !   batch_size (out): Will contain, on return, the number of values in the !
    !                     returned batch. The argument 'lookup_result' must be !
    !                     at least that big, unless want_result is .false.     !
    !   table (inout): The hash table object. Declared inout, so that hits and !
    !                  misses can be counted and stored inside the object.     !
    !   key1..key5 (in): The keys (constituting a 'keyset') which correspond   !
    !                    to the sought data. Key1 is mandatory, further keys   !
    !                    are optional. The number of supplied keys must match  !
    !                    the number of keys declared when initialising the     !
    !                    hash table.                                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, June 2012                                     !
    !==========================================================================!

    implicit none

    ! Arguments
    real(kind=DP), intent(inout), pointer      :: ret_ptr(:)
    integer, intent(out)                       :: batch_size
    type(HT_HASH_TABLE), intent(inout), target :: table
    integer, intent(in)                        :: key1
    integer, intent(in), optional              :: key2
    integer, intent(in), optional              :: key3
    integer, intent(in), optional              :: key4
    integer, intent(in), optional              :: key5

    ! Local variables
!   character(len=*), parameter     :: myself = 'hash_table_lookup_ptr'
    type(HT_STORAGE_NODE), pointer  :: ptr, prev_ptr

    ! ------------------------------------------------------------------------

    ! Corner case of dummy, empty hash_table
    if(table%max_slots == 0) then
       batch_size = -1
       return
    end if

    ! Locate the record
    call hash_table_locate(ptr, prev_ptr, table, key1, key2, key3, key4, key5)

    if(associated(ptr)) then
       ! Record found, return data
       batch_size = ptr%batch_size
       ret_ptr => ptr%node_data(:)
    else
       ! No such record, notify caller
       batch_size = -1
    end if

  end subroutine hash_table_lookup_ptr
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_lookup_ptr_nocount(ret_ptr, batch_size, table, &
       key1, key2, key3, key4, key5)
    !==========================================================================!
    ! Retrieves the data (batch of real(kind=DP)'s) corresponding to the given !
    ! keyset. The data is returned by pointer (ret_ptr), to avoid deep-copying !
    ! the contents. The size of the batch is returned in batch_size. If the    !
    ! element is absent in the hash table, a batch_size of -1 is returned and  !
    ! ret_ptr remains untouched. No counting of hits and misses is performed,  !
    ! allowing the table to be intent(in).                                     !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   ret_ptr (inout): (Pointer to) the buffer where the data resides.       !
    !                    Intent(inout) is used, because this is not touched    !
    !                    when the element is not present in the hash table.    !
    !                    The pointer remains valid ONLY until the next call to !
    !                    hash_table_remove(), hash_table_cleanup(),            !
    !                    hash_table_free() or hash_table_add(overwrite=.true.).!
    !   batch_size (out): Will contain, on return, the number of values in the !
    !                     returned batch. The argument 'lookup_result' must be !
    !                     at least that big, unless want_result is .false.     !
    !   table (in): The hash table object.                                     !
    !   key1..key5 (in): The keys (constituting a 'keyset') which correspond   !
    !                    to the sought data. Key1 is mandatory, further keys   !
    !                    are optional. The number of supplied keys must match  !
    !                    the number of keys declared when initialising the     !
    !                    hash table.                                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, June 2012                                     !
    !==========================================================================!

    implicit none

    ! Arguments
    real(kind=DP), intent(inout), pointer      :: ret_ptr(:)
    integer, intent(out)                       :: batch_size
    type(HT_HASH_TABLE), intent(in), target    :: table
    integer, intent(in)                        :: key1
    integer, intent(in), optional              :: key2
    integer, intent(in), optional              :: key3
    integer, intent(in), optional              :: key4
    integer, intent(in), optional              :: key5

    ! Local variables
    type(HT_STORAGE_NODE), pointer  :: ptr, prev_ptr

    ! ------------------------------------------------------------------------

    ! NB. This subroutine is called from OMP contexts. Do not utils_trace() it.

    ! Corner case of dummy, empty hash_table
    if(table%max_slots == 0) then
       batch_size = -1
       return
    end if

    ! Locate the record
    call hash_table_locate_nocount(ptr, prev_ptr, table, &
         key1, key2, key3, key4, key5)

    if(associated(ptr)) then
       ! Record found, return data
       batch_size = ptr%batch_size
       ret_ptr => ptr%node_data(:)
    else
       ! No such record, notify caller
       batch_size = -1
    end if

  end subroutine hash_table_lookup_ptr_nocount
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_remove(table, key1, key2, key3, key4, key5)
    !==========================================================================!
    ! Removes an entry from the hash table.                                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table (inout): The table on which to operate.                          !
    !   key1..key5:    The keys defining the entry to be removed.              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, June 2012                                     !
    !==========================================================================!

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(inout), target :: table
    integer, intent(in)                        :: key1
    integer, intent(in), optional              :: key2
    integer, intent(in), optional              :: key3
    integer, intent(in), optional              :: key4
    integer, intent(in), optional              :: key5

    ! Local variables
    type(HT_STORAGE_NODE), pointer  :: ptr, prev_ptr

    ! ------------------------------------------------------------------------

    ! Corner case of dummy, empty hash_table
    if(table%max_slots == 0) return

    ! Locate the record
    call hash_table_locate(ptr, prev_ptr, table, key1, key2, key3, key4, key5)

    ! Remove it
    call hash_table_remove_record(ptr, prev_ptr, table)

  end subroutine hash_table_remove
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_list(table, verbosity_level)
    !==========================================================================!
    ! Gives some statistics on the contents of the hash table.                 !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table (in): The hash table object.                                     !
    !   verbosity_level (in): 0 - basic stats only,                            !
    !                         1 - also show metadata (keys, hash, batch_size,  !
    !                             pointers, etc.)                              !
    !                         2 - also show data itself (bulky!).              !
    !                         The default is 0.                                !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, May 2012                                      !
    ! Reduced complexity from O(N) to O(1) (for verbosity_level 0)             !
    ! by Jacek Dziedzic, April 2019.                                           !
    !==========================================================================!

    use constants, only: SIZEOF_DOUBLE
    use timer, only: wrappers_etime, timer_clock

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(in), target :: table
    integer, intent(in), optional           :: verbosity_level

    ! Local variables
    integer                                 :: slot
    integer                                 :: i_element
    integer                                 :: i_slot
    integer                                 :: i_chained
    type(HT_STORAGE_NODE), pointer          :: ptr
    real(kind=DP)                           :: n_total
    real(kind=DP)                           :: mem_size_slots_mib
    real(kind=DP)                           :: mem_size_chains_mib
    real(kind=DP)                           :: mem_overhead_mib
    integer                                 :: local_verbosity_level
    character(len=21)                       :: s_tag
    character(len=*), parameter             :: myself = 'hash_table_list'

    ! ------------------------------------------------------------------------

    call timer_clock(myself,1)

    if(present(verbosity_level)) then
       local_verbosity_level = verbosity_level
    else
       local_verbosity_level = 0
    end if

    if(table%max_slots == 0) then
       s_tag = ' #T: DUMMY'
    else
       s_tag = ' #T: '//table%table_name
    end if

    write(table%listunit,'(a,a)') &
         '+---------------------------------------------------+', s_tag

    if(table%max_slots == 0) then
       write(table%listunit,'(a,a)') &
            '|         A dummy hash table (holds 0 slots).       |',s_tag
       write(table%listunit,'(a,a)') &
            '+---------------------------------------------------+',s_tag
       call timer_clock(myself,2)
       return
    end if

    if(table%n_keys == -1) then
       write(table%listunit,'(a,a)') &
            '|         A hash table that has been freed.         |',s_tag
       write(table%listunit,'(a,a)') &
            '+---------------------------------------------------+',s_tag
       call timer_clock(myself,2)
       return
    end if

    write(table%listunit,'(a,i1,a,a)') '|        A hash table indexed by ', &
         table%n_keys,' keys.            |', s_tag
    write(table%listunit,'(a,a)') &
         '+---------------------------------------------------+', s_tag
    write(table%listunit,'(a,f24.6,a,a)') '| Timestamp: ', wrappers_etime(), &
         '               |', s_tag
    write(table%listunit,'(a,i11,a,i11,a,f5.1,a,a)') &
         '| Occupied slots: ',table%n_slots, '/', table%max_slots,' (', &
         real(table%n_slots,kind=DP)/real(table%max_slots,kind=DP)*100.0_DP, &
         '%). |', s_tag
    write(table%listunit,'(a,a)') &
         '| Number of stored elements:                        |', s_tag
    n_total = (real(table%n_slots,kind=DP) + real(table%n_chained,kind=DP))
    if(n_total == 0.0_DP) n_total = 1.0_DP ! corner case with no elements
    write(table%listunit,'(a,i11,a,f6.2,a,a)') '| - in slots:  ', &
         table%n_slots, ' (', &
         real(table%n_slots,kind=DP) / n_total * 100.0_DP, &
         '%).               |', s_tag
    write(table%listunit,'(a,i11,a,f6.2,a,a)') '| - in chains: ', &
         table%n_chained, ' (', &
         real(table%n_chained,kind=DP) / n_total * 100.0_DP,&
         '%).               |',s_tag
    write(table%listunit,'(a,i11,a,i11,a,a)') '| = in toto:   ', &
         table%n_slots+table%n_chained,' (max: ',table%max_elements, &
         ').      |',s_tag
    write(table%listunit,'(a,a)') &
         '| Number of accesses:                               |', s_tag
    n_total = table%n_hits + table%n_misses
    if(n_total == 0.0_DP) n_total = 1.0_DP ! corner case with no accesses
    write(table%listunit,'(a,i13,a,f6.2,a,a)') '| - hits:     ', &
         table%n_hits, &
         ' (', &
         real(table%n_hits,kind=DP) / n_total * 100.0_DP, '%).              |', &
         s_tag
    write(table%listunit,'(a,i13,a,f6.2,a,a)') '| - misses:   ', &
         table%n_misses, ' (', &
         real(table%n_misses,kind=DP) / n_total * 100.0_DP, '%).              |',&
         s_tag
    write(table%listunit,'(a,i13,a,a)') '| = in toto:  ', &
         table%n_hits + table%n_misses, '.                        |', s_tag
    write(table%listunit,'(a,a)') &
         '| Number of record adds:                            |', s_tag
    n_total = real(table%n_adds_ok,kind=DP) + &
         real(table%n_adds_when_full,kind=DP) + &
         real(table%n_adds_rejected,kind=DP)
    if(n_total == 0.0_DP) n_total = 1.0_DP ! corner case with no adds
    write(table%listunit,'(a,i15,a,f6.2,a,a)') '| - successful:          ', &
         table%n_adds_ok, ' (', &
         real(table%n_adds_ok,kind=DP) / n_total * 100.0_DP, '%). |', &
         s_tag
    write(table%listunit,'(a,i15,a,f6.2,a,a)') '| - succesf. w/overflow: ', &
         table%n_adds_when_full, ' (', &
         real(table%n_adds_when_full,kind=DP) / n_total * 100.0_DP, '%). |',&
         s_tag
    write(table%listunit,'(a,i15,a,f6.2,a,a)') '| - rejected (ht full):  ', &
         table%n_adds_rejected, ' (', &
         real(table%n_adds_rejected,kind=DP) / n_total * 100.0_DP, '%). |',&
         s_tag
    write(table%listunit,'(a,i15,a,a)') '| = in toto:             ', &
         table%n_adds_ok + table%n_adds_when_full + table%n_adds_rejected,&
         '.           |', s_tag
    write(table%listunit,'(a,i5,a,a)') '| Longest chain (ever):  ', &
         table%longest_chain, '.                     |', s_tag
    write(table%listunit,'(a,i15,a,a)') '| Removed records:       ', &
         table%n_removed_recs, '.           |', s_tag
    write(table%listunit,'(a,i15,a,a)') '| Cleanups:              ', &
         table%n_cleanups, '.           |', s_tag
    write(table%listunit,'(a,i15,a,a)') '| Purges:                ', &
         table%n_purges, '.           |', s_tag

    if(local_verbosity_level > 0) then
       write(table%listunit,'(a,a)') &
       '+----------- Contents of the hash table: -----------+',s_tag
       i_element = 1
       i_slot = 1
       do slot = 1, table%max_slots
          if(table%slots(slot)%in_use) then
             ptr => table%slots(slot)
             i_chained = 1
             if(local_verbosity_level > 0) then
                write(table%listunit,'(a,a)') &
                     '|===================================================|', s_tag
             end if
             do while(associated(ptr))
                if(ptr%in_use) then
                   if(local_verbosity_level > 0) then
                      if(i_chained > 1) then
                         write(table%listunit,'(a,a)') &
                         '| ------------------------------------------------- |', &
                         s_tag
                      end if
                      write(table%listunit,'(a,i13,a,i13,a,a)') '| Used slot: ', &
                           i_slot, ', element: ', i_element,'  |', s_tag
                      write(table%listunit,'(a,i13,a,a)') '| Hash: ', slot, &
                           '                               |', s_tag
#ifndef DONT_HAVE_LOC
                      write(table%listunit,'(a,i24,a,a)') '| Node @    ', &
                           loc(ptr), '                |', s_tag
                      write(table%listunit,'(a,i24,a,a)') '| Storage @ ', &
                           loc(ptr%node_data), '                |', s_tag
#endif
                      write(table%listunit,'(a,i10,i10,i10,i10,i10,a,a)') &
                           '| - Keys: ', &
                           ptr%key1, ptr%key2, ptr%key3, ptr%key4, ptr%key5, &
                           '  |', s_tag
                      write(table%listunit,'(a,i11,a,a)') '| - Batch size: ', &
                           ptr%batch_size, '                         |', s_tag
                      write(table%listunit,'(a,i11,a,a)') '| - # hits:     ', &
                           ptr%n_hits, '                         |', s_tag
                      write(table%listunit,'(a,l1,a,a)') &
                           '| - contains_only_zeroes: ', ptr%contains_only_zeroes,&
                           '                         |', s_tag
                      if(local_verbosity_level > 1 .and. &
                           .not. ptr%contains_only_zeroes) then
                         write(table%listunit,*) ptr%node_data(1:ptr%batch_size)
                      end if
                   end if

                   i_slot = i_slot + 1
                   i_element = i_element + 1
                   ptr => ptr%next
                   if(associated(ptr)) then
                      if(local_verbosity_level > 0) then
                         write(table%listunit,'(a,a)') &
                              '| ^                               &
                              &                  |', s_tag
                         write(table%listunit,'(a,i8,a,a)') '| | chained entry', &
                              i_chained, '                           |', s_tag
                         write(table%listunit,'(a,a)') &
                              '| |                               &
                              &                  |', s_tag
                      end if
                      i_chained = i_chained + 1
                      i_slot = i_slot - 1 ! inside a chain, does not count as slots
                   end if
                 end if ! ptr%in_use
             end do ! while associated(ptr)
          end if ! slot in use
       end do ! over slots
    end if

    mem_size_slots_mib = table%total_slot_data_size * &
         SIZEOF_DOUBLE/1024.0_DP/1024.0_DP
    mem_size_chains_mib = table%total_chain_data_size * &
         SIZEOF_DOUBLE/1024.0_DP/1024.0_DP

    write(table%listunit,'(a,a)') &
         '+---------------------------------------------------+', s_tag
    write(table%listunit,'(a,a)') &
         '| Estimated memory usage:                           |', s_tag
    mem_overhead_mib = real(node_size,kind=DP)*real(table%max_slots,kind=DP) / &
         1024.0_DP/1024.0_DP
    write(table%listunit,'(a,f11.3,a,a)') '| - node overhead: ', &
         mem_overhead_mib, ' MiB.                 |', s_tag
    write(table%listunit,'(a,f11.3,a,a)') '| - data (slots):  ', &
         mem_size_slots_mib, ' MiB.                 |', s_tag
    write(table%listunit,'(a,f11.3,a,a)') '| - data (chains): ', &
         mem_size_chains_mib, ' MiB.                 |', s_tag
    write(table%listunit,'(a,f11.3,a,a)') '| = in toto:       ', &
         mem_overhead_mib + mem_size_slots_mib + mem_size_chains_mib, &
         ' MiB.                 |', s_tag
    write(table%listunit,'(a,a)') &
         '+---------------------------------------------------+', s_tag

    call timer_clock(myself,2)

  end subroutine hash_table_list
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_mem_usage(table, mem_total, mem_slots, mem_chains, &
       mem_overhead)
    !==========================================================================!
    ! Returns the memory usage of a given hash table. All values are in MiB    !
    ! and correspond to memory used on current MPI rank.                       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table (in): The hash table object.                                     !
    !   mem_total (out): Total memory used by the HT is returned here.         !
    !   mem_slots (out, opt): Memory used by data in slots is returned here.   !
    !   mem_chains (out, opt): Memory used by data in chains is returned here. !
    !   mem_overhead (out, opt): Memory lost to HT overheads is returned here. !
    !                                                                          !
    !   mem_total = mem_slots + mem_chains + mem_overhead.                     !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2020.                              !
    !==========================================================================!

    use constants, only: SIZEOF_DOUBLE

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(in), target :: table
    real(kind=DP), intent(out)              :: mem_total
    real(kind=DP), intent(out), optional    :: mem_slots
    real(kind=DP), intent(out), optional    :: mem_chains
    real(kind=DP), intent(out), optional    :: mem_overhead

    ! Local variables
    real(kind=DP) :: mem_size_slots_mib
    real(kind=DP) :: mem_size_chains_mib
    real(kind=DP) :: mem_overhead_mib

    ! -------------------------------------------------------------------------

    mem_size_slots_mib = table%total_slot_data_size * &
         SIZEOF_DOUBLE/1024.0_DP/1024.0_DP
    mem_size_chains_mib = table%total_chain_data_size * &
         SIZEOF_DOUBLE/1024.0_DP/1024.0_DP
    mem_overhead_mib = real(node_size,kind=DP)*real(table%max_slots,kind=DP) / &
         1024.0_DP/1024.0_DP

    mem_total = mem_size_slots_mib + mem_size_chains_mib + mem_overhead_mib
    if(present(mem_slots)) then
       mem_slots = mem_size_slots_mib
    end if
    if(present(mem_chains)) then
       mem_chains = mem_size_chains_mib
    end if
    if(present(mem_overhead)) then
       mem_overhead = mem_overhead_mib
    end if

  end subroutine hash_table_mem_usage
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_cleanup_at_will(table, fraction_to_remove)
    !==========================================================================!
    ! If the table is overfull, calls cleanup.                                 !
    ! This is useful if cleanup has been prevented earlier in hash_table_add   !
    ! in OMP contexts. A cleanup invalidates pointers to hash table entries,   !
    ! so being able to temporarily prevent cleanups is beneficial.             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table (inout): The HT to clean up.                                     !
    !   fraction_to_remove (in, opt): If provided, gets passed to cleanup().   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2015                                    !
    !==========================================================================!

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(inout), target  :: table
    real(kind=DP), intent(in), optional         :: fraction_to_remove

    ! ------------------------------------------------------------------------

    if(table%n_slots + table%n_chained >= table%max_elements) then
       call hash_table_cleanup(table, fraction_to_remove)
    end if

  end subroutine hash_table_cleanup_at_will
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_stash_prepare(stash, max_stash_size, &
      overwrite, overfill_strategy)
    !==========================================================================!
    ! Sets up a stash of a specified size, with specified values for the       !
    ! 'overwrite' and 'overfill_strategy' attributes.                          !
    !                                                                          !
    ! Stashes can be used as temporary thread-local storage for hash tables,   !
    ! to avoid having to lock a table during read accesses -- as all updates   !
    ! can be arranged to only happen to the stash and not the table.           !
    ! Each thread's stash can then be merged with the hash table ('committed') !
    ! within a critical section.                                               !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   stash (inout): An allocatable array that will be allocated by this     !
    !                  routine. This will constitute the stash.                !
    !   max_stash_size (in): The size of the stash (in MiB).                   !
    !   overwrite (in):        } The attributes will be stored in the stash and!
    !   overfill_strategy (in): } will be enforced when the stash is committed.!
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in Jan 2017.                                   !
    ! Generalised to overfill_strategy by Jacek Dziedzic, March 2019.          !
    !==========================================================================!

    use constants, only: SIZEOF_DOUBLE
    use utils, only: utils_alloc_check, utils_character_to_real, &
         utils_logical_to_real

    implicit none

    ! Arguments
    real(kind=DP), intent(inout), allocatable :: stash(:)
    integer, intent(in)   :: max_stash_size
    logical, intent(in)   :: overwrite
    character, intent(in) :: overfill_strategy

    ! Local variables
    integer :: max_stash_size_elems
    integer :: ierr
    integer, parameter :: MiB = 1048576

    ! ------------------------------------------------------------------------

    max_stash_size_elems = max_stash_size * MiB / SIZEOF_DOUBLE

    allocate(stash(max_stash_size_elems), stat=ierr)
    call utils_alloc_check('hash_table_stash_prepare','stash', ierr)

    stash(1) = real(max_stash_size_elems, kind=DP)
    stash(2) = C_STASH_HEADER_SIZE+1 ! points to 5, first available record
    stash(3) = utils_logical_to_real(overwrite)
    stash(4) = utils_character_to_real(overfill_strategy)
  ! stash(5) ... data goes here

  end subroutine hash_table_stash_prepare
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_stash(stash, table, data_to_add, ndata, &
       key1, key2, key3, key4, key5)
    !==========================================================================!
    ! Adds data (a batch of real(kind=DP) values) to a stash.                  !
    ! See comment in hash_table_stash_prepare() for general info on stashes.   !
    ! If the stash is full, it will ignore adds that would overflow it.        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   stash (inout): The stash to which data is to be added.                 !
    !   table (in): The hash table to which the stash will be later committed. !
    !               This is only used for debug checks (table_name and so on). !
    !   data_to_add (in): The data to be added to the stash.                   !
    !   ndata (in): The number of values in the batch to be added.             !
    !   key1..key5 (in): The keys (constituting a 'keyset') which will be used !
    !                    to retrieve the data. Key1 is mandatory, further keys !
    !                    are optional. The number of supplied keys must match  !
    !                    the number of keys declared when initialising the     !
    !                    hash table.                                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in Jan 2017, using hash_table_add as template. !
    !==========================================================================!

    use rundat, only: pub_debug
    use utils, only: utils_isnan, utils_abort

    implicit none

    ! Arguments
    real(kind=DP), intent(inout)       :: stash(:)
    type(HT_HASH_TABLE), intent(in), target :: table
    real(kind=DP), intent(in)          :: data_to_add(:)
    integer, intent(in)                :: ndata
    integer, intent(in)                :: key1
    integer, intent(in), optional      :: key2
    integer, intent(in), optional      :: key3
    integer, intent(in), optional      :: key4
    integer, intent(in), optional      :: key5

    ! Local variables
    integer :: local_key1, local_key2, local_key3, local_key4, local_key5
    integer :: i

    ! ------------------------------------------------------------------------

    ! Corner case of dummy, empty hash_table
    if(table%max_slots == 0) return

    ! Sanity check keys, put in zeroes where optional args are not supplied
    call hash_table_process_keys(table, &
         local_key1, local_key2, local_key3, local_key4, local_key5, & ! out
         key1, key2, key3, key4, key5)                                 ! in

    ! Add to stash
    call hash_table_stash_data(stash, table, data_to_add, ndata, &
         local_key1, local_key2, local_key3, local_key4, local_key5)

    if(pub_debug) then
       do i = 1, ndata
          if(utils_isnan(data_to_add(i))) then
             call utils_abort('Sanity check failed. NaN detected when stashing &
             &data associated with '//trim(table%table_name),i,ndata);
          end if
       end do
    end if

  end subroutine hash_table_stash
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_stash_data(stash, table, data_to_add, ndata, &
       key1, key2, key3, key4, key5)
    !==========================================================================!
    ! Worker routine for stashing a single record.                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   stash (inout): The stash to which data is to be added.                 !
    !   table (in): The hash table to which the stash will be later committed. !
    !               This is only used for debug checks (table_name and so on). !
    !   data_to_add(in): The data added.                                       !
    !   ndata(in):       The number of elements in data_to_add.                !
    !   key1..key5:      The keys.                                             !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in Jan 2017, using hash_table_add_data as      !
    ! template.                                                                !
    !==========================================================================!

    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(inout)       :: stash(:)
    type(HT_HASH_TABLE), intent(in), target :: table
    real(kind=DP), intent(in)          :: data_to_add(:)
    integer, intent(in)                :: ndata
    integer, intent(in)                :: key1, key2, key3, key4, key5

    ! Local variables
    integer :: max_stash_size_elems
    integer :: first_free
    character(len=*), parameter    :: myself = 'hash_table_stash_data'

    ! ------------------------------------------------------------------------

    max_stash_size_elems = stash(1)
    first_free = stash(2)

    ! jd: Exit if we don't have enough space for stashing ndata elements
    if(max_stash_size_elems - first_free - C_STASH_OVERHEAD <= ndata) return

    stash(first_free) = key1
    first_free = first_free + 1
    stash(first_free) = key2
    first_free = first_free + 1
    stash(first_free) = key3
    first_free = first_free + 1
    stash(first_free) = key4
    first_free = first_free + 1
    stash(first_free) = key5
    first_free = first_free + 1
    stash(first_free) = ndata
    first_free = first_free + 1
    stash(first_free:first_free + ndata-1) = data_to_add(1:ndata)
    first_free = first_free + ndata
    call utils_assert(first_free <= max_stash_size_elems, &
         myself//': stash overflow')
    stash(2) = first_free

  end subroutine hash_table_stash_data
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_stash_commit(stash, table)
    !==========================================================================!
    ! Commits a stash to a hash table.                                         !
    ! See comment in hash_table_stash_prepare() for general info on stashes.   !
    ! The behaviour in case of overwrite and in case of a full table depends on!
    ! the values passed in 'overwrite' and 'dont_add_if_full' in the call to   !
    ! hash_table_stash_prepare().                                              !
    !                                                                          !
    ! This subroutine is not protected from race conditions, so it's up to the !
    ! caller to ensure the hash table is only updated from one thread at a     !
    ! time.                                                                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   stash (in): The stash to be committed.                                 !
    !   table (in/out): The table to which the stash is to be committed.       !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2017.                               !
    !==========================================================================!

    use utils, only: utils_abort, utils_assert, utils_real_to_character, &
         utils_real_to_logical

    implicit none

    ! Arguments
    real(kind=DP), intent(in), allocatable :: stash(:)
    type(HT_HASH_TABLE), intent(inout), target :: table

    ! Local variables
    integer   :: max_stash_size_elems
    integer   :: first_free
    logical   :: overwrite
    character :: overfill_strategy
    integer   :: cur_pos
    integer   :: ndata
    integer   :: key1, key2, key3, key4, key5
    character(len=*), parameter :: myself = 'hash_table_stash_commit'

    ! ------------------------------------------------------------------------

    max_stash_size_elems = stash(1)
    first_free = stash(2)
    overwrite = utils_real_to_logical(stash(3))
    overfill_strategy = utils_real_to_character(stash(4))
    cur_pos = 5

    if(overfill_strategy /= 'F' .and. overfill_strategy /= 'C' .and. &
         overfill_strategy /= 'N' .and. overfill_strategy /= 'B') then
       call utils_abort(myself//': Unrecognised overfill_strategy: '//&
            overfill_strategy)
    end if

    ! Add all stashed records in sequence
    do while(cur_pos < first_free)
       call utils_assert(cur_pos + C_STASH_OVERHEAD-1 <= max_stash_size_elems, &
            myself//': stash bounds overrun (overhead)', cur_pos, &
            max_stash_size_elems)
       key1 = stash(cur_pos)
       cur_pos = cur_pos + 1
       key2 = stash(cur_pos)
       cur_pos = cur_pos + 1
       key3 = stash(cur_pos)
       cur_pos = cur_pos + 1
       key4 = stash(cur_pos)
       cur_pos = cur_pos + 1
       key5 = stash(cur_pos)
       cur_pos = cur_pos + 1
       ndata = stash(cur_pos)
       cur_pos = cur_pos + 1
       call utils_assert(cur_pos + ndata <= max_stash_size_elems, &
            myself//': stash bounds overrun (data)', cur_pos, ndata, &
            max_stash_size_elems)
       ! jd: Add to table. The ugliness below is a consequence of the fact that
       !     the stash does not hold info on how many keys there were, it
       !     already has zeroes for missing keys, which are indistinguishable
       !     from valid keys equal to zero. process_keys() could figure that
       !     out, but that would add overheads in a common routine, I prefer
       !     the slight ugliness here.
       select case(table%n_keys)
       case(5)
          call hash_table_add(table, stash(cur_pos:cur_pos+ndata-1), ndata, &
               key1, key2, key3, key4, key5, overwrite = overwrite, &
               overfill_strategy = overfill_strategy)
       case(4)
          call hash_table_add(table, stash(cur_pos:cur_pos+ndata-1), ndata, &
               key1, key2, key3, key4, overwrite = overwrite, &
               overfill_strategy = overfill_strategy)
       case(3)
          call hash_table_add(table, stash(cur_pos:cur_pos+ndata-1), ndata, &
               key1, key2, key3, overwrite = overwrite, &
               overfill_strategy = overfill_strategy)
       case(2)
          call hash_table_add(table, stash(cur_pos:cur_pos+ndata-1), ndata, &
               key1, key2, overwrite = overwrite, &
               overfill_strategy = overfill_strategy)
       case(1)
          call hash_table_add(table, stash(cur_pos:cur_pos+ndata-1), ndata, &
               key1, overwrite = overwrite, &
               overfill_strategy = overfill_strategy)
       case default
          call utils_abort(myself//': Unsupported number of keys in stash', &
               table%n_keys)
       end select
       cur_pos = cur_pos + ndata
    end do

  end subroutine hash_table_stash_commit
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_stash_free(stash)
    !==========================================================================!
    ! Deallocates memory associated with a stash.                              !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !--------------------------------------------------------------------------!
    !  stash (inout): The stash whose memory is to be deallocated.             !
    !==========================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(inout), allocatable :: stash(:)

    ! Local variables
    integer :: ierr
    character(len=*), parameter :: myself = 'hash_table_stash_free'

    ! ------------------------------------------------------------------------

    deallocate(stash, stat=ierr)
    call utils_dealloc_check(myself,'stash',ierr)

  end subroutine hash_table_stash_free
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  ! ---------------------------------------------------------------------------
  ! ---- P R I V A T E   S U B R O U T I N E S   A N D   F U N C T I O N S ----
  ! ---------------------------------------------------------------------------

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  integer function hash_table_generate_hash(key1, key2, key3, key4, key5, &
       max_slots)
    !==========================================================================!
    ! The hashing function itself. The returned hash is always between 1 and   !
    ! max_slots.                                                               !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   key1..key5 (in): The keyset from which the hash will be calculated.    !
    !                    All five keys are mandatory (use zeroes for keys      !
    !                    which are not used.                                   !
    !   max_slots (in):  The same value that was used to initialise the table. !
    !--------------------------------------------------------------------------!
    ! Return value:                                                            !
    !   The calculated hash.                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, May 2012                                      !
    !==========================================================================!

    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: key1, key2, key3, key4, key5
    integer, intent(in) :: max_slots

    ! Local variables
    integer, parameter :: HT_CONST_0 = 17
    integer, parameter :: HT_PRIME_1 = 1399
    integer, parameter :: HT_PRIME_2 = 193741
    integer, parameter :: HT_PRIME_3 = 21928757
    integer, parameter :: HT_PRIME_4 = 1378477297
    integer, parameter :: HT_PRIME_5 = 31

    ! ------------------------------------------------------------------------

    if(max_slots <= 0) then
       call utils_abort('hash_table_generate_hash: Bad max_slots.')
    end if

    hash_table_generate_hash = mod(iabs((((((HT_CONST_0 + key1) * HT_PRIME_1 + &
         HT_CONST_0 + key2) * HT_PRIME_2 + &
         HT_CONST_0 + key3) * HT_PRIME_3 + &
         HT_CONST_0 + key4) * HT_PRIME_4 + &
         HT_CONST_0 + key5) * HT_PRIME_5 ),max_slots) + 1

    ! jd: Trap extreme corner case where the value under iabs is -2147483648.
    !     The absolute value of this has no representation in the integer type
    !     (as maxint is 2147483647), and iabs returns a negative value (!), its
    !     own argument. This then leads to a negative hash.
    if(hash_table_generate_hash < 0) then
       hash_table_generate_hash = 1
    end if

  end function hash_table_generate_hash
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_locate(found_ptr, prev_ptr, table, &
       key1, key2, key3, key4, key5)
    !==========================================================================!
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table (inout): The hash table object. Declared inout, so that hits and !
    !                  misses can be counted and stored inside the object.     !
    !   key1..key5 (in): The keys (constituting a 'keyset') which correspond   !
    !                    to the sought data. Key1 is mandatory, further keys   !
    !                    are optional. The number of supplied keys must match  !
    !                    the number of keys declared when initialising the     !
    !                    hash table.                                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, May 2012                                      !
    !==========================================================================!

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(inout), target  :: table
    type(HT_STORAGE_NODE), intent(out), pointer :: found_ptr
    type(HT_STORAGE_NODE), intent(out), pointer :: prev_ptr
    ! ^note: ifort 18.0.[01], due to issue #03143573, required a
    ! @workaround: declaration of intent(out) pointers as intent(inout)
    integer, intent(in)                         :: key1
    integer, intent(in), optional               :: key2
    integer, intent(in), optional               :: key3
    integer, intent(in), optional               :: key4
    integer, intent(in), optional               :: key5

    ! Local variables
    integer                         :: hash
    integer                         :: local_key1, local_key2
    integer                         :: local_key3, local_key4, local_key5
    type(HT_STORAGE_NODE), pointer  :: ptr, prev
    logical                         :: found

    ! ------------------------------------------------------------------------

    ! Sanity check keys, put in zeroes where optional args are not supplied
    call hash_table_process_keys(table, &
         local_key1, local_key2, local_key3, local_key4, local_key5, & ! out
         key1, key2, key3, key4, key5)                                 ! in

    ! Generate hash
    hash = hash_table_generate_hash(local_key1, local_key2, local_key3, &
         local_key4, local_key5, table%max_slots)

    if(.not. table%slots(hash)%in_use) then
       ! No entry for this hash, thus element absent
       table%n_misses = table%n_misses + 1
       nullify(found_ptr)
       nullify(prev_ptr)
    else
       ! There is entry for this hash, check if the keyset matches
       table%n_hits = table%n_hits + 1
       if(local_key1 == table%slots(hash)%key1 .and. &
            local_key2 == table%slots(hash)%key2 .and. &
            local_key3 == table%slots(hash)%key3 .and. &
            local_key4 == table%slots(hash)%key4 .and. &
            local_key5 == table%slots(hash)%key5) then
          ! Keys in the slot match, return what we found
          table%slots(hash)%n_hits = table%slots(hash)%n_hits + 1
          found_ptr => table%slots(hash)
          nullify(prev_ptr)
       else
          ! Keys in the slot do not match, go along the chain
          found = .false.
          prev => table%slots(hash)
          ptr => table%slots(hash)%next
          do while (associated(ptr))
             if(local_key1 == ptr%key1 .and. &
                  local_key2 == ptr%key2 .and. &
                  local_key3 == ptr%key3 .and. &
                  local_key4 == ptr%key4 .and. &
                  local_key5 == ptr%key5) then

                ! Keys in the chain match, return what we found
                ptr%n_hits = ptr%n_hits + 1
                found_ptr => ptr
                prev_ptr => prev
                found = .true.
                exit
             end if
             prev => ptr
             ptr => ptr%next
          end do
          ! If we get here and haven't found the record, this means that the
          ! keyset was not found in the chain (or there was no chain at all),
          ! which means it's a miss after all -- we only got something with the
          ! same hash.
          if(.not. found) then
             table%n_hits = table%n_hits - 1             ! Miss, after all
             table%n_misses = table%n_misses + 1         ! Miss, after all
             nullify(found_ptr)
             nullify(prev_ptr)
          end if
       end if
    end if

  end subroutine hash_table_locate
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_locate_nocount(found_ptr, prev_ptr, table, &
       key1, key2, key3, key4, key5)
    !==========================================================================!
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table (in): The hash table object.                                     !
    !   key1..key5 (in): The keys (constituting a 'keyset') which correspond   !
    !                    to the sought data. Key1 is mandatory, further keys   !
    !                    are optional. The number of supplied keys must match  !
    !                    the number of keys declared when initialising the     !
    !                    hash table.                                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, May 2012                                      !
    !==========================================================================!

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(in), target     :: table
    type(HT_STORAGE_NODE), intent(out), pointer :: found_ptr
    type(HT_STORAGE_NODE), intent(out), pointer :: prev_ptr
    ! ^note: ifort 18.0.[01], due to issue #03143573, required a
    ! @workaround: declaration of intent(out) pointers as intent(inout)
    integer, intent(in)                         :: key1
    integer, intent(in), optional               :: key2
    integer, intent(in), optional               :: key3
    integer, intent(in), optional               :: key4
    integer, intent(in), optional               :: key5

    ! Local variables
    integer                         :: hash
    integer                         :: local_key1, local_key2
    integer                         :: local_key3, local_key4, local_key5
    type(HT_STORAGE_NODE), pointer  :: ptr, prev
    logical                         :: found

    ! ------------------------------------------------------------------------

    ! Sanity check keys, put in zeroes where optional args are not supplied
    call hash_table_process_keys(table, &
         local_key1, local_key2, local_key3, local_key4, local_key5, & ! out
         key1, key2, key3, key4, key5)                                 ! in

    ! Generate hash
    hash = hash_table_generate_hash(local_key1, local_key2, local_key3, &
         local_key4, local_key5, table%max_slots)

    if(.not. table%slots(hash)%in_use) then
       ! No entry for this hash, thus element absent
       nullify(found_ptr)
       nullify(prev_ptr)
    else
       ! There is entry for this hash, check if the keyset matches
       if(local_key1 == table%slots(hash)%key1 .and. &
            local_key2 == table%slots(hash)%key2 .and. &
            local_key3 == table%slots(hash)%key3 .and. &
            local_key4 == table%slots(hash)%key4 .and. &
            local_key5 == table%slots(hash)%key5) then
          ! Keys in the slot match, return what we found
          found_ptr => table%slots(hash)
          nullify(prev_ptr)
       else
          ! Keys in the slot do not match, go along the chain
          found = .false.
          prev => table%slots(hash)
          ptr => table%slots(hash)%next
          do while (associated(ptr))
             if(local_key1 == ptr%key1 .and. &
                  local_key2 == ptr%key2 .and. &
                  local_key3 == ptr%key3 .and. &
                  local_key4 == ptr%key4 .and. &
                  local_key5 == ptr%key5) then

                ! Keys in the chain match, return what we found
                found_ptr => ptr
                prev_ptr => prev
                found = .true.
                exit
             end if
             prev => ptr
             ptr => ptr%next
          end do
          ! If we get here and haven't found the record, this means that the
          ! keyset was not found in the chain (or there was no chain at all),
          ! which means it's a miss after all -- we only got something with the
          ! same hash.
          if(.not. found) then
             nullify(found_ptr)
             nullify(prev_ptr)
          end if
       end if
    end if

  end subroutine hash_table_locate_nocount
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_add_data(table, data_to_add, ndata, hash, &
       key1, key2, key3, key4, key5, &
       overwrite, overfill_strategy)
    !==========================================================================!
    ! Worker routine for adding an entry to the table.                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table(inout):    The table operated on.                                !
    !   data_to_add(in): The data added.                                       !
    !   ndata(in):       The number of elements in data_to_add.                !
    !   hash(in):        The hash corresponding to key1..key5.                 !
    !   key1..key5:      The keys.                                             !
    !   overwrite(inout), optional. If not specified or specified as .false.,  !
    !                               the presence of an entry matching the      !
    !                               keyset will trigger an abort. If specified !
    !                               and .true., a potentially duplicate entry  !
    !                               will be overwritten. If this happens, on   !
    !                               exit 'overwrite' will remain .true..       !
    !   overfill_strategy (in): What to do if the HT is full already:          !
    !      'F' - allow overfilling (ignore max_elements).                      !
    !      'C' - cleanup (remove some elements before adding this one).        !
    !      'N' - do not add the element (ignore request).                      !
    !      'B' - break (abort).                                                !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, May 2012                                      !
    ! Extended with prevent_cleanup by Jacek Dziedzic, April 2015.             !
    ! Extended with dont_add_if_full by Jacek Dziedzic, January 2017.          !
    ! Generalised with overfill_strategy by Jacek Dziedzic, March 2019.        !
    !==========================================================================!

    use utils, only: utils_abort, utils_assert, utils_alloc_check, &
         utils_dealloc_check

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(inout), target :: table
    real(kind=DP), intent(in)          :: data_to_add(:)
    integer, intent(in)                :: ndata
    integer, intent(in)                :: hash
    integer, intent(in)                :: key1, key2, key3, key4, key5
    logical, intent(inout), optional   :: overwrite
    character, intent(in)              :: overfill_strategy

    ! Local variables
    character(len=*), parameter    :: myself = 'hash_table_add_data'
    integer                        :: ierr
    integer                        :: cur_chain_length
    type(HT_STORAGE_NODE), pointer :: ptr, next_ptr
    logical                        :: local_overwrite
    logical                        :: already_overwritten
    logical                        :: was_full

    ! ------------------------------------------------------------------------

    call utils_assert(hash >= 1 .and. hash <= table%max_slots, &
         myself//': Invalid hash.',hash)

    if(present(overwrite)) then
       local_overwrite = overwrite
       overwrite = .false. ! From now on acts as a returned value
    else
       local_overwrite = .false.
    end if

    if(overfill_strategy /= 'F' .and. overfill_strategy /= 'C' .and. &
         overfill_strategy /= 'N' .and. overfill_strategy /= 'B') then
       call utils_abort(myself//': Unrecognised overfill_strategy: '//&
            overfill_strategy)
    end if

    was_full = .false.

    ! Check if hash table full. If so, depending on flags, either
    ! remove 10% or so records, do nothing and continue adding record if
    ! cleanup expressly forbidden (useful in OMP contexts), do nothing and
    ! don't add the record (useful if FIFO eviction is undesired) or abort
    if(table%n_slots + table%n_chained >= table%max_elements) then
       ! jd: 'B': Break with an error
       if(overfill_strategy == 'B') then
          call utils_abort(myself//': Hash table '//trim(table%table_name)//&
               ' full and overfill_strategy is "BREAK ON FULL".')
       end if
       ! jd: 'N': Ignore request
       if(overfill_strategy == 'N') then
          table%n_adds_rejected = table%n_adds_rejected + 1
          return
       end if
       ! jd: 'C': Cleanup and continue adding
       if(overfill_strategy == 'C') then
          call hash_table_cleanup(table)
       end if
       ! jd: 'F': Allow overfilling (no-op, just add record)
       table%n_adds_when_full = table%n_adds_when_full + 1
       table%n_adds_ok = table%n_adds_ok - 1 ! counter the subsequent increment
       was_full = .true.
    end if

    table%n_adds_ok = table%n_adds_ok + 1

    if(.not. table%slots(hash)%in_use) then
       ! Simple case of no collision and no overwrite
       ! - allocate storage
       allocate(table%slots(hash)%node_data(ndata),stat=ierr)
       call utils_alloc_check(myself,'node_data',ierr)
       ! - copy data
       table%slots(hash)%node_data(1:ndata) = data_to_add(1:ndata)
       ! - set metadata
       table%slots(hash)%n_hits = 0
       table%slots(hash)%key1 = key1
       table%slots(hash)%key2 = key2
       table%slots(hash)%key3 = key3
       table%slots(hash)%key4 = key4
       table%slots(hash)%key5 = key5
       table%slots(hash)%batch_size = ndata
       table%slots(hash)%contains_only_zeroes = .false. ! @fixthis
       nullify(table%slots(hash)%next)
       table%n_slots = table%n_slots + 1
       table%slots(hash)%in_use = .true.
       table%total_slot_data_size = table%total_slot_data_size + ndata

    else
       ! Hash collision
       ! - Go through the chain, inspecting keys in it in case this is an
       !   overwrite. Start from the slot itself.
       ptr => table%slots(hash)
       next_ptr => table%slots(hash)%next
       cur_chain_length = 1
       already_overwritten = .false.
       do while (associated(ptr))

          if(ptr%key1 == key1 .and. &
               ptr%key2 == key2 .and. &
               ptr%key3 == key3 .and. &
               ptr%key4 == key4 .and. &
               ptr%key5 == key5) then
             ! Keyset already present: overwrite or abort
             if(local_overwrite) then
                ! De-allocate and reallocate if batch size is different
                if(ptr%batch_size /= ndata) then
                   deallocate(ptr%node_data,stat=ierr)
                   call utils_dealloc_check(myself, &
                        'node_data', ierr)
                   allocate(ptr%node_data(ndata),stat=ierr)
                   call utils_alloc_check(myself,'node_data',&
                        ierr)
                end if
                ! Overwrite with new data
                ptr%node_data(1:ndata) = data_to_add(1:ndata)
                ptr%batch_size = ndata
                ptr%contains_only_zeroes = .false. ! @fixthis
                already_overwritten = .true.
                overwrite = .true.
                ! Don't count overwrites as adds
                if(.not. was_full) then
                   table%n_adds_ok = table%n_adds_ok - 1
                else
                   table%n_adds_when_full = table%n_adds_when_full - 1
                end if
                exit
             else
                ! No permission to overwrite. Abort.
                call utils_abort(myself//': '//trim(table%table_name)//&
                     ': Keyset already present and overwrite is .false.',&
                     key1,key2,key3,key4,key5)
             end if
          end if

          ! Move along in the chain
          cur_chain_length = cur_chain_length + 1
          if(.not. associated(next_ptr)) then
             exit
          end if
          ptr => next_ptr
          next_ptr => ptr%next
       end do
       if (.not. already_overwritten) then
          ! Reached the end of the chain, thus it's a hash collision with no
          ! overwrite. Add the node at the end of the chain.
          ! - allocate storage at the end of the chain
          allocate(next_ptr,stat=ierr)
          call utils_alloc_check(myself,'node',ierr)
          allocate(next_ptr%node_data(ndata),stat=ierr)
          call utils_alloc_check(myself,'node_data',ierr)
          ! - copy data
          next_ptr%node_data(1:ndata) = data_to_add(1:ndata)
          ! - set metadata
          next_ptr%n_hits = 0
          next_ptr%key1 = key1
          next_ptr%key2 = key2
          next_ptr%key3 = key3
          next_ptr%key4 = key4
          next_ptr%key5 = key5
          next_ptr%batch_size = ndata
          next_ptr%contains_only_zeroes = .false. ! @fixthis
          next_ptr%in_use = .true.
          table%n_chained = table%n_chained + 1
          table%total_chain_data_size = table%total_chain_data_size + ndata
          ptr%next => next_ptr
          nullify(next_ptr%next)
          table%longest_chain = max(table%longest_chain,cur_chain_length)
       end if
    end if

  end subroutine hash_table_add_data
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_process_keys(table, &              ! in
       out_key1, out_key2, out_key3, out_key4, out_key5, & ! out
       in_key1, in_key2, in_key3, in_key4, in_key5)        ! in
    !==========================================================================!
    ! Internal routine for pre-processing keys. Keys from 'in' are copied to   !
    ! keys from 'out', absent 'in' keys are replaced with 0 in 'out'.          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table(in):               The table in question.                        !
    !   in_key1..in_key5(in):    Input keys, last four are optional.           !
    !   out_key1..out_key5(out): Output keys.                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, May 2012                                      !
    !==========================================================================!

    use rundat, only: pub_debug
    use utils, only: utils_assert, utils_abort, utils_int_to_str

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(in) :: table
    integer, intent(out)            :: out_key1
    integer, intent(out)            :: out_key2
    integer, intent(out)            :: out_key3
    integer, intent(out)            :: out_key4
    integer, intent(out)            :: out_key5
    integer, intent(in)             :: in_key1
    integer, intent(in), optional   :: in_key2
    integer, intent(in), optional   :: in_key3
    integer, intent(in), optional   :: in_key4
    integer, intent(in), optional   :: in_key5

    ! Local variables
    character(len=*), parameter     :: myself = 'hash_table_process_keys'
    integer                         :: n_keys_passed

    ! ------------------------------------------------------------------------

    if(pub_debug) then
       ! Sanity check optional args
       if(present(in_key5)) then
          call utils_assert(present(in_key4) .and. present(in_key3) .and. &
               present(in_key2), &
               myself//': Inconsistency in the number of passed keys for '&
               //trim(table%table_name))
       end if
       if(present(in_key4)) then
          call utils_assert(present(in_key3) .and. present(in_key2), &
               myself//': Inconsistency in the number of passed keys for '&
               //trim(table%table_name))
       end if
       if(present(in_key3)) then
          call utils_assert(present(in_key2), &
               myself//': Inconsistency in the number of passed keys for ' &
               //trim(table%table_name))
       end if
    end if

    ! Count keys, use zeroes where optional args were not supplied
    n_keys_passed = 1
    out_key1 = in_key1
    if(present(in_key2)) then
       n_keys_passed = n_keys_passed + 1
       out_key2 = in_key2
    else
       out_key2 = 0
    end if
    if(present(in_key3)) then
       n_keys_passed = n_keys_passed + 1
       out_key3 = in_key3
    else
       out_key3 = 0
    end if
    if(present(in_key4)) then
       n_keys_passed = n_keys_passed + 1
       out_key4 = in_key4
    else
       out_key4 = 0
    end if
    if(present(in_key5)) then
       n_keys_passed = n_keys_passed + 1
       out_key5 = in_key5
    else
       out_key5 = 0
    end if

    if(pub_debug .and. table%n_keys /= n_keys_passed) then
       call utils_abort(myself//': Number of passed keys ('//&
            trim(utils_int_to_str(n_keys_passed))//&
            ') is inconsistent with the number of keys &
            &the hash table was initialised with ('//&
            trim(utils_int_to_str(table%n_keys))//&
            ') (or, likely, you missed initialisation entirely) for '//&
            trim(table%table_name))
    end if

  end subroutine hash_table_process_keys
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_remove_record(ptr, prev_ptr, table)
    !==========================================================================!
    ! Removes a record pointed by 'ptr' from 'table'.                          !
    ! The need to pass 'prev_ptr' stems from the fact that the remainder of the!
    ! chain may need to be-reattached to the previous node.                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table(inout): The table operated on.                                   !
    !   ptr(inout):   Pointer to the record that is to be removed.             !
    !   prev_ptr(inout): Pointer to a record whose 'next' is ptr.              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, June 2012                                     !
    !==========================================================================!

    use rundat, only: pub_debug
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(HT_STORAGE_NODE), intent(inout), pointer :: ptr
    type(HT_STORAGE_NODE), intent(inout), pointer :: prev_ptr
    type(HT_HASH_TABLE), intent(inout), target    :: table

    ! Local variables
    integer :: ndata
    character(len=*), parameter    :: myself = 'hash_table_remove_record'
    integer                        :: ierr
    type(HT_STORAGE_NODE), pointer :: ptr_next_next
    real(kind=DP), parameter :: garbage = 9D99

    ! ------------------------------------------------------------------------

    if(.not. associated(ptr)) then
       call utils_abort(myself//': Attempt to remove a non-existent record.')
    end if

    if(.not. ptr%in_use) then
       call utils_abort(myself//&
            ': Attempt to remove a record that is not in use.')
    end if

    table%n_removed_recs = table%n_removed_recs + 1

    ndata = ptr%batch_size

    ! OK, the record is definitely present.

    ! Deallocate storage, if allocated
    if(.not. ptr%contains_only_zeroes) then
       !write(*,*) '@DE data (remove record) ', loc(ptr%node_data)
       if(pub_debug) then
          ! Overwrite deallocated data with garbage for easier race detection
          ptr%node_data(:) = garbage
       endif
       deallocate(ptr%node_data,stat=ierr)
       call utils_dealloc_check(myself,'node_data', ierr)
    else
       call utils_abort('contains_only_zeroes not implemented')
    end if

    ! Are we removing a record in a slot?
    if(.not. associated(prev_ptr)) then
       ! YES
       ! - Check if there was a chain attached to the slot
       if(associated(ptr%next)) then
          ! There was a chain, move the first chain link into the slot by
          ! deep-copying
          ptr%in_use = .true.
          ptr%n_hits = ptr%next%n_hits
          ptr%key1 = ptr%next%key1
          ptr%key2 = ptr%next%key2
          ptr%key3 = ptr%next%key3
          ptr%key4 = ptr%next%key4
          ptr%key5 = ptr%next%key5
          ptr%batch_size = ptr%next%batch_size
          ! Data has just been deallocated, must reallocate (batch size can
          ! differ).
          allocate(ptr%node_data(ptr%batch_size),stat=ierr)
          call utils_alloc_check(myself,'node_data', ierr)
          ptr%node_data = ptr%next%node_data
          if(pub_debug) then
             ! Overwrite deallocated data with garbage for easier race detection
             ptr%next%node_data(:) = garbage
          end if
          deallocate(ptr%next%node_data,stat=ierr)
          call utils_dealloc_check(myself,'node_data',ierr)
          ptr%contains_only_zeroes = ptr%next%contains_only_zeroes
          ptr_next_next => ptr%next%next
          deallocate(ptr%next, stat=ierr)
          call utils_dealloc_check(myself,'node',ierr)
          ptr%next => ptr_next_next
          ! We now have one fewer chained record
          table%n_chained = table%n_chained - 1
          table%total_chain_data_size = table%total_chain_data_size - ndata
       else
          ! There was no chain, mark the slot as unused
          ptr%in_use = .false.
          ptr%n_hits = 0
          ! Also, we now have one fewer occupied slot
          table%n_slots = table%n_slots - 1
          table%total_slot_data_size = table%total_slot_data_size - ndata
       end if
    else
       ! NO, we are removing a record in a chain
       ! - Tie the previous link to the next one
       prev_ptr%next => ptr%next
       ! - Get rid of this link. Storage has been deallocated already.
       deallocate(ptr,stat=ierr)
       call utils_dealloc_check(myself,'node', ierr)
       ! - We now have one fewer chained record
       table%n_chained = table%n_chained - 1
       table%total_chain_data_size = table%total_chain_data_size - ndata
    end if

  end subroutine hash_table_remove_record
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_cleanup(table, fraction_to_remove)
    !==========================================================================!
    ! Removes certain percentage of records ('fraction_to_remove') from table. !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   table (inout): The HT to clean up.                                     !
    !   fraction_to_remove (in, opt): If provided, specifies the fraction of   !
    !                                 elements to remove (e.g. 0.05).          !
    !                                 If omitted, it uses the cleanup_fraction !
    !                                 stored in the HT.                        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, June 2012                                     !
    !==========================================================================!

    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(inout), target  :: table
    real(kind=DP), intent(in), optional         :: fraction_to_remove

    ! Local variables
    character(len=*), parameter    :: myself = 'hash_table_cleanup'
    type(HT_STORAGE_NODE), pointer :: ptr
    type(HT_STORAGE_NODE), pointer :: nullptr
    real(kind=DP)                  :: local_fraction_to_remove
    integer                        :: n_elements_to_remove
    integer                        :: slot, first_slot

    ! ------------------------------------------------------------------------

    call timer_clock(myself//'_['//trim(table%table_name)//']',1)

    table%n_cleanups = table%n_cleanups + 1

    ! Corner case of dummy, empty hash_table
    if(table%max_slots == 0) then
       call timer_clock(myself//'_['//trim(table%table_name)//']',2)
       return
    end if

    ! Use provided fraction to remove, or default to the one stored in the table
    if(present(fraction_to_remove)) then
       local_fraction_to_remove = fraction_to_remove
    else
       local_fraction_to_remove = table%cleanup_fraction
    end if

    call utils_assert(local_fraction_to_remove >= 0.0_DP .and. &
         local_fraction_to_remove <= 1.0_DP, myself// &
         ': Fraction of elements to remove must be in [0, 1].')

    n_elements_to_remove = &
         ceiling(real(table%max_elements,kind=DP) * local_fraction_to_remove)

    ! Never try to remove more elements than there are actually there
    n_elements_to_remove = min(n_elements_to_remove, &
         table%n_slots + table%n_chained)

    call utils_assert(n_elements_to_remove > 0, myself// &
         ': n_elements_to_remove must be positive.')

    nullify(nullptr)

    do while (n_elements_to_remove > 0)
       first_slot = table%first_slot_for_cleanup

       do slot = first_slot, table%max_slots
          ptr => table%slots(slot)
          if(ptr%in_use) then
             call hash_table_remove_record(ptr, nullptr, table)
             n_elements_to_remove = n_elements_to_remove - 1
             if(n_elements_to_remove == 0) exit
          end if
       end do

       table%first_slot_for_cleanup = mod(slot-1,table%max_slots) + 1
    end do

    call timer_clock(myself//'_['//trim(table%table_name)//']',2)

  end subroutine hash_table_cleanup
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_iterate(batch_size, ptr, slot, &
       table, key1, key2, key3, key4, key5, eot, data_ptr, want_result, &
       lookup_result)
    !--------------------------------------------------------------------------!
    ! Iterates over the hash table without the need to provide keys.           !
    ! Each call returns the contents of a single node in the hash table, and   !
    ! adjusts the arguments 'ptr' and 'slot' so that the next call will return !
    ! the subsequent element.                                                  !
    !                                                                          !
    ! To begin iterating, pass NULL for 'ptr', and something that is not -1    !
    ! for 'slot'. For subsequent iterations leave 'ptr' and 'slot' be.         !
    !                                                                          !
    ! If a node is succesfully found, the stored data is returned in           !
    ! 'lookup_result' (unless .false. was passed in 'want_result'),            !
    ! 'batch_size' is set to the size of the stored data, 'key1..5' are set to !
    ! the values of the keys indexing the node, and .false. is returned for    !
    ! 'eot'. 'ptr' and 'slot' are adjusted in such a way that the next call    !
    ! will start at the next node that is in use.                              !
    !                                                                          !
    ! If a node is not found (table is empty, or all nodes have been iterated  !
    ! over), 'eot' is returned as .true, 'batch_size' is 0, 'key1..5' are -1,  !
    ! and 'lookup_result' remains untouched.                                   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   batch_size (out): Will contain, on return, the number of values in the !
    !                     returned data. The argument 'lookup_result' must be  !
    !                     at least that big, unless 'want_result' is .false.   !
    !   ptr (inout): Pass NULL to start iterating. Pass a previously returned  !
    !                value to continue iterating.                              !
    !   slot (inout): Pass anything but -1 to start iterating. Pass a          !
    !                 previously returned value to continue iterating.         !
    !   table (in): The hash table object.                                     !
    !   key1..key5 (out): The keys (constituting a 'keyset') which correspond  !
    !                     to the returned data. All keys are mandatory. Zeros  !
    !                     will be returned for parts of keyset that are absent.!
    !                     -1 will be returned for all keys if there are no more!
    !                     nodes.                                               !
    !   eot (out): .false. will be returned if a node was found. .true. will   !
    !              be returned if there are no more nodes.                     !
    !   data_ptr (out, opt): If passed, on exit it will point to the data.     !
    !                        This can be useful if you do not want to allocate !
    !                        storage for it. Be mindful that the pointer       !
    !                        becomes invalid after any cleanup, add, etc.      !
    !                        If eot is .true. on exit, data_ptr will be NULL.  !
    !   want_result (in, optional): If .false. is passed, lookup_result will   !
    !                               not be written to or even referenced.      !
    !                               Defaults to .true. when omitted.           !
    !   lookup_result (out, opt): The buffer where the data will be retrieved  !
    !                        to. The caller is responsible for supplying a     !
    !                        buffer that is big enough to hold the data.       !
    !                        If want_result is .false. this buffer will        !
    !                        NOT BE REFERENCED.                                !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   - Only the first 'batch_size' elements in lookup_result will be        !
    !     written to. This is done for the sake of performance.                !
    !     As 'lookup_result' is an intent(out) parameter, this may trigger     !
    !     debug-mode checks with overzealous compilers -- technically it's our !
    !     responsibility to write to the entire array.                         !
    !   - Note the semantics of 'eot' which acts like C++'s eof flag. It is    !
    !     only returned as .true. if a node was sought, but not found. When    !
    !     the last node is successfully returned, 'eot' is NOT .true..         !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, May 2017                                      !
    !==========================================================================!

    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    integer, intent(out)                          :: batch_size
    type(HT_STORAGE_NODE), intent(inout), pointer :: ptr
    integer, intent(inout)                        :: slot
    type(HT_HASH_TABLE), intent(in), target       :: table
    integer, intent(out)                          :: key1
    integer, intent(out)                          :: key2
    integer, intent(out)                          :: key3
    integer, intent(out)                          :: key4
    integer, intent(out)                          :: key5
    logical, intent(out)                          :: eot
    real(kind=DP), intent(out), pointer, optional :: data_ptr(:)
    logical, intent(in), optional                 :: want_result
    real(kind=DP), intent(out), optional          :: lookup_result(:)

    ! Local variables
    logical :: loc_want_result
    integer :: cur_slot
    character(len=*), parameter :: myself = 'hash_table_iterate'

    ! -------------------------------------------------------------------------

    call utils_assert(table%max_slots /=0, &
         myself//': Hash table not yet initialised or previously freed.')

    if(slot == -1) goto 100

    if(present(want_result)) then
       loc_want_result = want_result
    else
       loc_want_result = .true.
    end if

    if(loc_want_result) then
       call utils_assert(present(lookup_result), myself//': Inconsistent args')
    end if

    ! jd: If starting from the beginning, find first slot in use
    if(.not. associated(ptr)) then
       do cur_slot = 1, table%max_slots
          if(table%slots(cur_slot)%in_use) then
             ptr => table%slots(cur_slot)
             slot = cur_slot
             exit
          end if
       end do
    end if

    ! jd: If all slots are unused, table is empty and we're done
    if(.not. associated(ptr)) goto 100

    ! jd: OK, we located a node that is in use
    key1 = ptr%key1
    key2 = ptr%key2
    key3 = ptr%key3
    key4 = ptr%key4
    key5 = ptr%key5
    batch_size = ptr%batch_size
    eot = .false.
    if(present(data_ptr)) data_ptr => ptr%node_data
    if(loc_want_result) then
       lookup_result(1:batch_size) = ptr%node_data(1:batch_size)
    end if

    ! jd: Aim for the next element in chain
    ptr => ptr%next

    ! jd: If there is none, scan for the next slot in use
    if(.not. associated(ptr)) then
       do cur_slot = slot+1, table%max_slots
          if(table%slots(cur_slot)%in_use) then
             ptr => table%slots(cur_slot)
             slot = cur_slot
             return
          end if
       end do
       ! jd: If we're here, we just found the last node. Return it, but make
       !     sure we return eot on next call
       slot = -1
       return
    else
       return
    end if

100 continue ! jd: We've reached the end
    key1 = -1
    key2 = -1
    key3 = -1
    key4 = -1
    key5 = -1
    batch_size = 0
    slot = -1
    ptr => NULL()
    if(present(data_ptr)) then
       data_ptr => NULL()
    end if
    eot = .true.
    ! jd: lookup_result remains untouched
    return

  end subroutine hash_table_iterate
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine hash_table_self_test
    !==========================================================================!
    ! Internal use only.                                                       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   NB: Doesn't test cleanup, because it would be hard to figure out       !
    !       records removed in cleanup.
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, June 2012                                     !
    !==========================================================================!

    use utils, only: utils_assert

    implicit none

    ! Local variables
    integer, parameter          :: MAX_KEY = 5000
    integer, parameter          :: NDATA = 50000
    integer, parameter          :: NRUNS = 4000
    type(HT_HASH_TABLE)         :: table
    integer                     :: run, prev_run, key1, key2, batch_size
    real(kind=DP)               :: r
    real(kind=DP)               :: mydata(1:NDATA)
    integer                     :: n_keys_added
    integer                     :: keys(2,NRUNS)
    logical                     :: overwrite
    integer                     :: key_pair_loc

    ! ------------------------------------------------------------------------

    write(*,*) 'Intializing a new hash table...'
    ! NB: The test only works fine when there is no cleanup involved, thus
    !     the maximum number of elements is set to NRUNS.
    call hash_table_init(table, 'TEST', 2, 2000, NRUNS)

    write(*,*) 'Dump of the empty hash table...'
    call hash_table_list(table, 0)

    write(*,*) 'Adding records...'
    ! Add some records
    do run = 1, NRUNS
       call random_number(r)
       key1 = int(r * real(MAX_KEY,kind=DP))
       call random_number(r)
       key2 = int(r * real(MAX_KEY,kind=DP))
       keys(1,run) = key1
       keys(2,run) = key2

       mydata(:) = real(run*run,kind=DP)

       write(*,*) 'Run ',run,': Adding ', key1, ' ', key2
       overwrite = .true.
       call hash_table_add(table, mydata, NDATA, key1, key2, &
            overwrite = overwrite)
       if(overwrite) then
          write(*,*) '!!! There was an overwrite.'
          do prev_run = 1, run-1
             if(keys(1,prev_run) == key1 .and. keys(2,prev_run) == key2) then
                write(*,*) 'Previous keys removed from run ', prev_run
                keys(1,prev_run) = -1
                keys(2,prev_run) = -1
             end if
          end do
       end if
    end do

    call hash_table_list(table, 0)

    write(*,*) 'Looking up records...'
    ! Look up all the records
    do run = 1, NRUNS
       key1 = keys(1,run)
       key2 = keys(2,run)
       if(key1 /= -1 .and. key2 /= -1) then
          write(*,*) 'Run ',run,': Looking up ', key1, ' ', key2
          call hash_table_lookup(mydata, batch_size, table, key1, key2)
          call utils_assert(batch_size == NDATA, &
               'Hash table self-test failed [1].', batch_size,NDATA)
          call utils_assert(abs(mydata(10) - real(run*run,kind=DP)) < 1D-5, &
               'Hash table self-test failed [2].', &
               mydata(10),real(run*run,kind=DP))
       end if
    enddo

    ! Not necessarily == NRUNS, because duplicates might have been randomized
    n_keys_added = table%n_slots + table%n_chained
    write(*,*) n_keys_added,' keys added after ',NRUNS,' hash_table_add''s.'

    write(*,*) 'Deleting records...'
    do run = 1, n_keys_added
       ! Pick a key to delete
       do while(.true.)
          call random_number(r)
          key_pair_loc = int(r * real(NRUNS,kind=DP)) + 1
          ! write(*,*) 'Perhaps deleting key_pair ', key_pair_loc
          key1 = keys(1,key_pair_loc)
          key2 = keys(2,key_pair_loc)
          if(key1 /= -1 .and. key2 /= -1) exit
       end do
       write(*,*) 'Keys to be deleted: ', key1, ' ', key2
       call hash_table_remove(table, key1, key2)
       ! Mark as removed
       keys(1,key_pair_loc) = -1
       keys(2,key_pair_loc) = -1
    end do

    write(*,*) 'All records removed.'

    call hash_table_list(table, 1)

    write(*,*) 'Freeing hash table...'
    call hash_table_free(table)

    write(*,*) 'Hash table after freeing...'
    call hash_table_list(table, 1)

  end subroutine hash_table_self_test
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  integer function hash_table_calc_cache_size(limit_in_mib, element_size)
    !==========================================================================!
    ! Returns the maximum number of elements in a cache, given the size of an  !
    ! element and the maximum size of the cache in MiB.                        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in June 2012.                                  !
    ! Improved by Jacek Dziedzic in February 2019 to take node overhead into   !
    ! account. This required moving it from utils, as the overhead constant    !
    ! is a part of this module.                                                !
    !==========================================================================!

    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    integer, intent(in) :: limit_in_mib
    integer, intent(in) :: element_size

    ! jd: Local variables
    real(kind=DP) :: element_size_in_mib
    real(kind=DP) :: node_overhead_in_mib

    ! ------------------------------------------------------------------------

    call utils_assert(element_size > 0, &
         'hash_table_calc_cache_size: Element size must be positive.')

    element_size_in_mib = real(element_size,kind=DP) / 1024.0_DP / 1024.0_DP
    node_overhead_in_mib = real(node_size,kind=DP) / 1024.0_DP / 1024.0_DP

    hash_table_calc_cache_size = floor(limit_in_mib / &
         (element_size_in_mib + node_overhead_in_mib))

  end function hash_table_calc_cache_size
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module hash_table
