!================================================================!
!                                                                !
!               Spherical wave expansion container               !
!                                                                !
!----------------------------------------------------------------!
! Defines a container for a spherical-wave expansion of NGWF     !
! pairs (SW_EX), and basic operations on it.                     !
!                                                                !
! Higher-level operations are defined in sw_expansion_mod.       !
!----------------------------------------------------------------!
! Split off from sw_expansion by Jacek Dziedzic in October 2016. !
!================================================================!

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

module sw_expansion_type

  use constants, only: DP, CRLF, stdout
  use hash_table, only: HT_HASH_TABLE
  use sw_resolution_of_identity, only: SW_QUALITY
  use utils, only: utils_trace_in, utils_trace_out

  implicit none

  private

  ! --------------------------------
  ! --- Major public subroutines ---
  ! --------------------------------

  public :: swx_init_container
  public :: swx_destroy_container
  public :: swx_invalidate_expansion
  public :: swx_merge_swex_pair

  ! --------------------------------
  ! --- Minor public subroutines ---
  ! --------------------------------
  public :: swx_coeffs_n_slots
  public :: swx_merge_coeffs_pair
  public :: swx_merge_dcoeffs_pair

  ! ----------------------------------
  ! --- Public types and variables ---
  ! ----------------------------------

  public :: SW_EX

  type SW_EX

     logical           :: initialised = .false.
     logical           :: uses_full_swri = .true.
     character(len=32) :: swex_name = '(*I have not been set up yet!*) '
     type(SW_QUALITY)  :: quality ! size of the aux. basis set (min_l, max_l, max_q)

     ! jd: xlat table from swex SW index to swri SW index
     integer, allocatable :: swex_sw_to_swri_sw(:) ! max_sws_per_centre elements

     ! jd: Stores expansion coefficients
     type(HT_HASH_TABLE) :: coeffs_ht

     ! jd: Stores expansion coefficients q-contracted with J_lq
     !     Only present with polarisable embedding
     type(HT_HASH_TABLE) :: dcoeffs_ht

     ! jd: Hash table file output control
     integer             :: hash_table_info_unit = -1
     character(len=256)  :: hash_table_info_filename = '_not_set_up_'

  end type SW_EX

  ! jd: Number of *slots* for coeffs_ht per atom. This will be scaled by the
  !     number of MPI ranks. Value comes from a fit to large calculations.
  !     Overheads are small -- 500000 slots corresponds to 60MB. The value of
  !     300 leads to typical numbers of slots in the order of 100000.
  !     This is *not* the total maximum number of coeffs, which is unbounded
  !     (HT is allowed to overfill).
  integer, parameter, public :: swx_slots_coeffs_per_atom = 300

  character(len=32), parameter :: hash_table_info_rootname = 'hash_table_swx'

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !---------------------------------------------------------------------------!
  ! ****               P U B L I C          R O U T I N E S              **** !
  !---------------------------------------------------------------------------!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_init_container(swex, &
       swri, swex_name, min_l, max_l, max_q, metric, rep_postfix)
    !==========================================================================!
    ! This subroutine initialises a SW_EX containter.                          !
    !  - initialises fields based on arguments.                                !
    !  - pre-calculates swex_sw_to_swri_sw xlat table.                         !
    !  - does sanity checks on the state of the SW_EX.                         !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 17/02/2015.                                 !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: SW_V, SW_O, METRIC_ELECTROSTATIC, METRIC_OVERLAP, &
         CRLF, NORMAL
    use rundat, only: pub_swx_output_detail
    use sw_resolution_of_identity, only: SW_RI
    use utils, only: utils_flushed_string_output, utils_alloc_check, &
         utils_assert, utils_int_to_str, utils_postfix_to_ngwf_set_name

    implicit none

    ! jd: Arguments
    type(SW_EX), intent(out)      :: swex
    type(SW_RI), intent(in)       :: swri
    character(len=*), intent(in)  :: swex_name
    integer, intent(in)           :: min_l
    integer, intent(in)           :: max_l
    integer, intent(in)           :: max_q
    integer, intent(in)           :: metric
    character(len=*), intent(in)  :: rep_postfix

    ! jd: Local variables
    integer :: ierr
    integer :: swri_sw, swex_sw
    integer :: swri_l, swri_m, swri_b, swri_q
    character(len=32) :: ngwf_set_name
    character(len=*), parameter :: myself = 'swx_init_container'

    ! -------------------------------------------------------------------------

    ngwf_set_name = utils_postfix_to_ngwf_set_name(rep_postfix)

    if(pub_swx_output_detail >= NORMAL) then
       call utils_flushed_string_output(&
            'SWX: Initialising '//trim(ngwf_set_name)//&
            '%['//trim(swex_name)//']->')
    end if

    swex%quality%min_l = min_l
    swex%quality%max_l = max_l
    swex%quality%max_q = max_q
    swex%quality%max_sws_per_centre = max_q * (1+max_l)**2
    swex%swex_name = swex_name

    ! jd: Check if the attached SWRI is compatible
    if(metric == METRIC_ELECTROSTATIC) then
       call utils_assert(swri%has_metric(SW_V), &
            'The SW resolution of identity ['//trim(swri%swri_name)//&
            '] attached to the SW expansion ['//trim(swex%swex_name)//&
            '] does not define a V metric matrix, which the expansion needs.')
    end if
    if(metric == METRIC_OVERLAP) then
       call utils_assert(swri%has_metric(SW_O), &
            'The SW resolution of identity ['//trim(swri%swri_name)//&
            '] attached to the SW expansion ['//trim(swex%swex_name)//&
            '] does not define an O metric matrix, which the expansion needs.')
    end if
    call utils_assert(swri%quality%min_l <= swex%quality%min_l, &
         'The SW resolution of identity ['//trim(swri%swri_name)//']''s value &
         &of min_l ('//trim(utils_int_to_str(swri%quality%min_l))//&
         ') is larger than min_l specified in the SW expansion ['//&
         trim(swex%swex_name)//'] ('//&
         trim(utils_int_to_str(swex%quality%min_l))//&
         '). The resolution of identity must encompass &
         &all values of l that will be used in the expansion.')
    call utils_assert(swri%quality%max_l >= swex%quality%max_l, &
         'The SW resolution of identity ['//trim(swri%swri_name)//']''s value &
         &of max_l ('//trim(utils_int_to_str(swri%quality%max_l))//&
         ') is smaller than max_l specified in the SW expansion ['//&
         trim(swex%swex_name)//'] ('//&
         trim(utils_int_to_str(swex%quality%max_l))//&
         '). The resolution of identity must encompass &
         &all values of l that will be used in the expansion.')
    call utils_assert(swri%quality%max_q >= swex%quality%max_q, &
         'The SW resolution of identity ['//trim(swri%swri_name)//']''s value &
         &of max_q ('//trim(utils_int_to_str(swri%quality%max_q))//&
         ') is smaller than max_q specified in the SW expansion ['//&
         trim(swex%swex_name)//'] ('//&
         trim(utils_int_to_str(swex%quality%max_q))//&
         '). The resolution of identity must encompass &
         &all values of q that will be used in the expansion.')

    allocate(swex%swex_sw_to_swri_sw(swex%quality%max_sws_per_centre),stat=ierr)
    call utils_alloc_check(myself,'swex%swex_sw_to_swri_sw',ierr)

    ! jd: Fill in swex_sw_to_swri_sw, count SWs in swex
    swex%swex_sw_to_swri_sw(:) = -1
    swex_sw = 1
    do swri_sw = 1, swri%quality%num_sws_per_centre
       swri_b = swri%sw_idx_to_bessel_idx(swri_sw)
       swri_l = swri%sw_idx_to_l(swri_sw)
       swri_m = swri%sw_idx_to_m(swri_sw)
       swri_q = swri%sw_idx_to_q_idx(swri_sw)
       if(swri_q > swex%quality%max_q) cycle
       if(swri_l < swex%quality%min_l) cycle
       if(swri_l > swex%quality%max_l) cycle

       swex%swex_sw_to_swri_sw(swex_sw) = swri_sw
       swex_sw = swex_sw + 1
    end do

    swex%quality%num_sws_per_centre = swex_sw-1

    ! jd: We use a full swri if and only if two conditions are satisfied
    !     - max_l = max_l
    !     - num_sws_per_centre = num_sws_per_centre.
    !     In a scenario where there is q truncation, we might end up with a
    !     full-RI even if different max_q's were specified, ie. in the end
    !     we have the same set in swri and swex. But just checking if the
    !     number of SWs is the same is insufficient -- say max_l = 1, max_q = 2
    !     vs max_l = 0, max_q = 8 -- both yield 8 SWs.
    swex%uses_full_swri = (swri%quality%max_l == max_l .and. &
         swex%quality%num_sws_per_centre == swri%quality%num_sws_per_centre)

    if(pub_on_root .and. pub_swx_output_detail >= NORMAL) then
       write(stdout,'(a,a,a,a)') &
            merge('full-RI','part-RI',swex%uses_full_swri),'->[', &
            trim(swri%swri_name), '].'
       write(stdout,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') &
            'SWX: - min_l: ', swex%quality%min_l, &
            ', max_l: ', swex%quality%max_l, &
            ', max_q: ', swex%quality%max_q, &
            ', num_sws: ', swex%quality%num_sws_per_centre, &
            ', max_sws: ', swex%quality%max_sws_per_centre,'.'
    end if

    call utils_assert(swex%quality%max_sws_per_centre /= 0, myself//&
         ': Cannot have zero spherical waves in an expansion.')
    call utils_assert(swex%quality%num_sws_per_centre /= 0, myself//&
         ': Cannot effectively have zero spherical waves in an expansion.')

    ! jd: Initialise persistent caches
    call swx_init_persistent_caches(swex, ngwf_set_name, swri%natoms)

    if(pub_swx_output_detail >= NORMAL) then
       call utils_flushed_string_output('SWX: Done.'//CRLF)
    end if

    swex%initialised = .true.

  end subroutine swx_init_container

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_destroy_container(swex, rep_postfix)
    !==========================================================================!
    ! Destroys a SW_EX container.                                              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 18/06/2012.                                 !
    ! Modified by Jacek Dziedzic on 08/04/2014.                                !
    ! Generalised for SW_EX by Jacek Dziedzic on 17/02/2015.                   !
    !==========================================================================!

    use constants, only: NORMAL
    use rundat, only: pub_swx_output_detail
    use utils, only: utils_flushed_string_output, utils_dealloc_check, &
         utils_postfix_to_ngwf_set_name

    implicit none

    ! jd: Arguments
    type(SW_EX), intent(inout)             :: swex
    character(len=*), intent(in), optional :: rep_postfix

    ! jd: Local variables
    integer :: ierr
    character(len=*), parameter :: myself = 'swx_destroy_container'

    ! -------------------------------------------------------------------------

    if(pub_swx_output_detail >= NORMAL) then
       if(present(rep_postfix)) then
          call utils_flushed_string_output(&
               'SWX: Destroying ['//&
               trim(utils_postfix_to_ngwf_set_name(rep_postfix))//'%'//&
               trim(swex%swex_name)//']... ')
       else
          call utils_flushed_string_output(&
               'SWX: Destroying ['//trim(swex%swex_name)//']... ')
       end if
    end if

    deallocate(swex%swex_sw_to_swri_sw,stat=ierr)
    call utils_dealloc_check(myself,'swex%swex_sw_to_swri_sw',ierr)

    ! jd: Free persistent caches
    call swx_free_persistent_caches(swex)

    swex%initialised = .false.

    if(pub_swx_output_detail >= NORMAL) then
       call utils_flushed_string_output('done'//CRLF)
    end if

  end subroutine swx_destroy_container

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_invalidate_expansion(swex, postfix)
    !==========================================================================!
    ! Invalidates the expansion coefficients (coeffs, dcoeffs) of an SW_EX     !
    ! container. Hash tables are purged.                                       !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2016.10.14.                                 !
    !==========================================================================!

    use constants, only: NORMAL
    use hash_table, only: hash_table_purge
    use rundat, only: pub_pol_emb_qmstar, pub_swx_output_detail
    use utils, only: utils_flushed_string_output, utils_postfix_to_ngwf_set_name

    implicit none

    ! jd: Arguments
    type(SW_EX), intent(inout)   :: swex
    character(len=*), intent(in) :: postfix

    ! jd: Local variables
    character(len=*), parameter :: myself = 'swx_invalidate_expansion'

    ! -------------------------------------------------------------------------

    if(pub_swx_output_detail >= NORMAL) then
       call utils_flushed_string_output(&
            'SWX: Invalidating '//trim(utils_postfix_to_ngwf_set_name(postfix))&
            //'%['//trim(swex%swex_name)//'].'//CRLF)
    end if

    if(pub_pol_emb_qmstar) call hash_table_purge(swex%dcoeffs_ht)
    call hash_table_purge(swex%coeffs_ht)

  end subroutine swx_invalidate_expansion

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !---------------------------------------------------------------------------!
  ! ****               P R I V A T E        R O U T I N E S              **** !
  !---------------------------------------------------------------------------!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_init_persistent_caches(swex, ngwf_set_name, natoms)
    !==================================================================!
    ! Initialises the cache structures that need to persist across     !
    ! calls. Currently these are COEFFS and DCOEFFS.                   !
    !==================================================================!

    use comms, only: pub_my_proc_id
    use hash_table, only: hash_table_init
    use rundat, only: pub_pol_emb_qmstar
    use utils, only: utils_unit, utils_abort

    implicit none

    ! Arguments
    type(SW_EX), intent(inout)   :: swex
    character(len=9), intent(in) :: ngwf_set_name
    integer, intent(in)          :: natoms

    ! Local variables
    integer :: n_slots

    ! -------------------------------------------------------------------------

    ! jd: Open a file for debugging hash-tables, if desired
    swex%hash_table_info_unit = utils_unit()
    write(swex%hash_table_info_filename,'(a,a,a,a,a,a,i0,a)') &
         trim(hash_table_info_rootname), '_', trim(ngwf_set_name), &
         '_', trim(swex%swex_name), &
         '_proc_', pub_my_proc_id, '.log'
    open(swex%hash_table_info_unit, file=swex%hash_table_info_filename, err=10)

    n_slots = swx_coeffs_n_slots(natoms)
    call hash_table_init(swex%coeffs_ht, 'COEFFS', 5, &
         n_slots, n_slots, swex%hash_table_info_unit)

    if(pub_pol_emb_qmstar) then
       call hash_table_init(swex%dcoeffs_ht, 'DCOEFFS', 4, &
            n_slots, n_slots, swex%hash_table_info_unit)
    end if

    return

10  call utils_abort('Could not create file: '//&
         trim(swex%hash_table_info_filename))

  end subroutine swx_init_persistent_caches

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_free_persistent_caches(swex)
    !==================================================================!
    ! Frees the persistent caches.                                     !
    !==================================================================!

    use hash_table, only: hash_table_free
    use rundat, only: pub_pol_emb_qmstar
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    type(SW_EX), intent(inout) :: swex

    ! -------------------------------------------------------------------------

    if(pub_pol_emb_qmstar) call hash_table_free(swex%dcoeffs_ht)
    call hash_table_free(swex%coeffs_ht)

    close(swex%hash_table_info_unit, err=10)

    return

10  call utils_abort('Could not close file: '//swex%hash_table_info_filename)

  end subroutine swx_free_persistent_caches

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_merge_swex_pair(swexes3, &              ! out
       swexes1, swexes2, swri, metric, swex3_name, &     ! in
       swex1_factor, swex2_factor)                       ! in, opt
    !==========================================================================!
    ! Merges two swex pairs (if Bessel averaging in effect) or two swexes      !
    ! (in the absence of Bessel averaging). The merging results in a new swex  !
    ! (or pair) whose min_l and max_l encompass the two input swexes (or pairs)!
    ! and whose coeffs and dcoeffs are summed. Optional factors can be provided!
    ! to achieve e.g. swex coeff subtraction. A metric for the resultant       !
    ! swex (pair) needs to be provided.                                        !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Caller is responsible for freeing the resultant swexes3.               !
    !   Both swexes (or pairs) must be tied to the same swri.                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2017.                                   !
    !==========================================================================!

    use hash_table, only: HT_STORAGE_NODE, hash_table_iterate, &
         hash_table_lookup_nocount, hash_table_add
    use rundat, only: pub_dma_bessel_averaging, pub_pol_emb_qmstar
    use sw_resolution_of_identity, only: SW_RI
    use utils, only: utils_abort, utils_assert, utils_int_to_str, &
         utils_alloc_check, utils_dealloc_check

    implicit none

    ! jd: Arguments
    type(SW_EX), intent(out)     :: swexes3(2)
    type(SW_EX), intent(in)      :: swexes1(2)
    type(SW_EX), intent(in)      :: swexes2(2)
    type(SW_RI), intent(in)      :: swri
    integer, intent(in)          :: metric
    character(len=*), intent(in) :: swex3_name
    real(kind=DP), intent(in), optional :: swex1_factor
    real(kind=DP), intent(in), optional :: swex2_factor

    ! jd: Local variables
    integer :: n_swexes
    integer :: i_swex
    integer :: min_l, max_l, min_l1, max_l1, min_l2, max_l2, max_q
    integer :: i_centre, n_centres
    integer :: num_coeffs1, num_coeffs2
    integer :: ht1_slot
    integer :: centre_offset1, centre_offset2, centre_offset3
    integer :: num_sws1, num_sws2, num_sws3
    integer :: num_lms1, num_lms2, num_lms3
    integer :: coeffs_kind, dcoeffs_kind
    integer :: ngwf1, ngwf2, dummy_key4, dummy_key5
    type(HT_STORAGE_NODE), pointer :: ht1_ptr
    real(kind=DP), pointer :: dummy_data_ptr(:)
    logical :: eot
    integer :: ierr
    real(kind=DP), allocatable :: coeffs1(:), coeffs2(:), coeffs3(:)
    real(kind=DP), allocatable :: dcoeffs1(:), dcoeffs2(:), dcoeffs3(:)
    character(len=*), parameter :: myself = 'swx_merge_swex_pair'

    ! -------------------------------------------------------------------------

    call utils_abort(myself//': Four-index coeffs not supported')

    n_swexes = merge(2,1,pub_dma_bessel_averaging)

    call utils_assert(swexes1(1)%quality%max_q == swexes2(1)%quality%max_q, &
         myself//': max_q mismatch (1)', &
         swexes1(1)%quality%max_q, swexes2(1)%quality%max_q)
    if(n_swexes == 2) then
       call utils_assert(swexes1(1)%quality%min_l == swexes1(2)%quality%min_l, &
            myself//': min_l mismatch', &
            swexes1(1)%quality%min_l, swexes1(2)%quality%min_l)
       call utils_assert(swexes1(1)%quality%max_l == swexes1(2)%quality%max_l, &
            myself//': max_l mismatch', &
            swexes1(1)%quality%max_l, swexes1(2)%quality%max_l)
       call utils_assert(swexes1(2)%quality%max_q == swexes2(2)%quality%max_q, &
            myself//': max_q mismatch (2)', &
            swexes1(2)%quality%max_q, swexes2(2)%quality%max_q)
    end if

    min_l1 = swexes1(1)%quality%min_l
    max_l1 = swexes1(1)%quality%max_l
    min_l2 = swexes2(1)%quality%min_l
    max_l2 = swexes2(1)%quality%max_l
    min_l = min(min_l1,min_l2)
    max_l = max(max_l1,max_l2)

    ! jd: Go over both swexes if Bessel averaging
    do i_swex = 1, n_swexes
       max_q = swexes1(i_swex)%quality%max_q
       call swx_init_container(swexes3(i_swex), swri, &
            trim(swex3_name)//'-'//trim(utils_int_to_str(i_swex)), &
            min_l, max_l, max_q, metric, '#')

       num_sws1 = swexes1(i_swex)%quality%num_sws_per_centre
       num_sws2 = swexes2(i_swex)%quality%num_sws_per_centre
       num_sws3 = swexes3(i_swex)%quality%num_sws_per_centre
       num_lms1 = (1+max_l1-min_l1)*(1+max_l1+min_l1)
       num_lms2 = (1+max_l2-min_l2)*(1+max_l2+min_l2)
       num_lms3 = (1+max_l-min_l)*(1+max_l+min_l)

       allocate(coeffs1(2*num_sws1), stat=ierr)
       call utils_alloc_check(myself,'coeffs1',ierr)
       allocate(coeffs2(2*num_sws2), stat=ierr)
       call utils_alloc_check(myself,'coeffs2',ierr)
       allocate(coeffs3(2*num_sws3), stat=ierr)
       call utils_alloc_check(myself,'coeffs3',ierr)

       ! -----------------------------------------------------------------------
       ! jd: Iterate through the COEFF hash tables of swex-1.
       !     Find corresponding entries in swex-2. Merge them.
       ! -----------------------------------------------------------------------
       ht1_ptr => NULL()
       ht1_slot = 0
       do
          ! jd: Acquire coeffs from both swexes
          call hash_table_iterate(num_coeffs1, ht1_ptr, ht1_slot, &
               swexes1(i_swex)%coeffs_ht, ngwf1, ngwf2, coeffs_kind, &
               dummy_key4, dummy_key5, eot, dummy_data_ptr, .true., coeffs1)
          if(eot) exit
          call hash_table_lookup_nocount(coeffs2, num_coeffs2, &
               swexes2(i_swex)%coeffs_ht, ngwf1, ngwf2, coeffs_kind)
          call utils_assert(num_coeffs2 /= -1, myself//&
               ': Internal error locating coeffs', ngwf1, ngwf2, coeffs_kind)

          ! jd: Figure out if this is a 1-centre or 2-centre expansion
          if(num_coeffs1 == num_sws1) then
             n_centres = 1
             call utils_assert(num_coeffs2 == num_sws2, &
                  myself//': n_centres mismatch between two swexes (C1)')
          elseif(num_coeffs1 == 2*num_sws1) then
             n_centres = 2
             call utils_assert(num_coeffs2 == 2*num_sws2, &
                  myself//': n_centres mismatch between two swexes (C2)')
          else
             call utils_abort(myself//': Expansion is neither 1-centre nor &
                  &2-centre (C)', num_coeffs1, num_sws1, num_coeffs2, num_sws2)
          end if

          ! jd: Merge coeffs, translating SW indices
          coeffs3(:) = 0.0_DP
          do i_centre = 1, n_centres
             centre_offset1 = (i_centre - 1) * num_sws1 + 1
             centre_offset2 = (i_centre - 1) * num_sws2 + 1
             centre_offset3 = (i_centre - 1) * num_sws1 + 1
             call swx_merge_coeffs_pair(coeffs3(centre_offset3:), &
                  coeffs1(centre_offset1:), coeffs2(centre_offset2:), &
                  swri, swexes1(i_swex), swexes2(i_swex), &
                  swex1_factor, swex2_factor)
          end do

          call hash_table_add(swexes3(i_swex)%coeffs_ht, coeffs3, &
               max(num_coeffs1, num_coeffs2), ngwf1, ngwf2, coeffs_kind)

       end do ! HT elements (ngwf pairs)

       deallocate(coeffs1, stat=ierr)
       call utils_dealloc_check(myself,'coeffs1',ierr)
       deallocate(coeffs2, stat=ierr)
       call utils_dealloc_check(myself,'coeffs2',ierr)
       deallocate(coeffs3, stat=ierr)
       call utils_dealloc_check(myself,'coeffs3',ierr)

       ! -----------------------------------------------------------------------
       ! jd: Iterate through the DCOEFF hash tables of swex-1.
       !     Find corresponding entries in swex-2. Merge them.
       !     This needs to be in a separate loop from the one above, as the two
       !     HTs (coeffs_ht vs dcoeffs_ht) may iterate differently.
       ! -----------------------------------------------------------------------
       if(pub_pol_emb_qmstar) then
          ht1_ptr => NULL()
          ht1_slot = 0

          allocate(dcoeffs1(2*num_lms1), stat=ierr)
          call utils_alloc_check(myself,'dcoeffs1',ierr)
          allocate(dcoeffs2(2*num_lms2), stat=ierr)
          call utils_alloc_check(myself,'dcoeffs2',ierr)
          allocate(dcoeffs3(2*num_lms3), stat=ierr)
          call utils_alloc_check(myself,'dcoeffs3',ierr)

          do
             ! jd: Acquire dcoeffs from both swexes
             call hash_table_iterate(num_coeffs1, ht1_ptr, ht1_slot, &
                  swexes1(i_swex)%dcoeffs_ht, ngwf1, ngwf2, dcoeffs_kind, &
                  dummy_key4, dummy_key5, eot, dummy_data_ptr, .true., dcoeffs1)
             if(eot) exit
             call hash_table_lookup_nocount(dcoeffs2, num_coeffs2, &
                  swexes2(i_swex)%dcoeffs_ht, ngwf1, ngwf2, dcoeffs_kind)
             call utils_assert(num_coeffs2 /= -1, myself//&
                  ': Internal error locating dcoeffs', ngwf1, ngwf2, dcoeffs_kind)

             ! jd: Figure out if this is a 1-centre or 2-centre expansion
             if(num_coeffs1 * max_q == num_sws1) then
                n_centres = 1
                call utils_assert(num_coeffs2 * max_q == num_sws2, &
                     myself//': n_centres mismatch between two swexes (D1)')
             elseif(num_coeffs1 * max_q == 2*num_sws1) then
                n_centres = 2
                call utils_assert(num_coeffs2 * max_q == 2 * num_sws2, &
                     myself//': n_centres mismatch between two swexes (D2)')
             else
                call utils_abort(myself//': Expansion is neither 1-centre nor &
                     &2-centre (D)', num_coeffs1 * max_q, num_sws1, &
                     num_coeffs2 * max_q, num_sws2)
             end if

             ! jd: Merge dcoeffs, translating SW indices
             dcoeffs3(:) = 0.0_DP
             do i_centre = 1, n_centres
                centre_offset1 = (i_centre - 1) * num_lms1 + 1
                centre_offset2 = (i_centre - 1) * num_lms2 + 1
                centre_offset3 = (i_centre - 1) * num_lms3 + 1

                call swx_merge_dcoeffs_pair(dcoeffs3(centre_offset3:), &
                     dcoeffs1(centre_offset1:), dcoeffs2(centre_offset2:), &
                     swri, swexes1(i_swex), swexes2(i_swex), &
                     swex1_factor, swex2_factor)
             end do

             call hash_table_add(swexes3(i_swex)%dcoeffs_ht, dcoeffs3, &
                  max(num_coeffs1, num_coeffs2), ngwf1, ngwf2, dcoeffs_kind)

          end do ! HT elements (ngwf pairs)

          deallocate(dcoeffs1, stat=ierr)
          call utils_dealloc_check(myself,'dcoeffs1',ierr)
          deallocate(dcoeffs2, stat=ierr)
          call utils_dealloc_check(myself,'dcoeffs2',ierr)
          deallocate(dcoeffs3, stat=ierr)
          call utils_dealloc_check(myself,'dcoeffs3',ierr)

      end if ! pol_emb_qmstar?
    end do ! swexes (Bessel averaging)

  end subroutine swx_merge_swex_pair

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function swx_xlat_sw_idx(src_sw_idx, swri, src_swex, dest_swex_quality)
    !==========================================================================!
    ! Translates a SW index from src_swex to dest_swex. Both swexes must share !
    ! the same swri.                                                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2017.                                   !
    !==========================================================================!

    use sw_resolution_of_identity, only: SW_RI, SW_QUALITY, swri_sw_to_ablmq
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    integer, intent(in)     :: src_sw_idx
    type(SW_RI), intent(in) :: swri
    type(SW_EX), intent(in) :: src_swex
    type(SW_QUALITY), intent(in) :: dest_swex_quality

    ! jd: Local variables
    integer :: src_swri_sw_idx
    integer :: src_l, src_m, src_q
    integer :: dest_min_l, dest_max_l, dest_max_q
    integer :: src_min_l, src_max_q
    integer :: minloffset
    integer :: dummy_b
    real(kind=DP) :: dummy_a, dummy_q
    character(len=*), parameter :: myself = 'swx_xlat_sw_idx'

    ! -------------------------------------------------------------------------

    src_swri_sw_idx = src_swex%swex_sw_to_swri_sw(src_sw_idx)
    call swri_sw_to_ablmq(swri, src_swri_sw_idx, dummy_a, dummy_b, &
         src_l, src_m, dummy_q)
    src_q = swri%sw_idx_to_q_idx(src_swri_sw_idx)
    src_max_q = src_swex%quality%max_q
    dest_max_q = dest_swex_quality%max_q
    dest_min_l = dest_swex_quality%min_l
    dest_max_l = dest_swex_quality%max_l
    src_min_l = src_swex%quality%min_l

    if(src_max_q /= dest_max_q) then
       call utils_abort(myself//': max_q not equal', src_max_q, dest_max_q)
    end if
    if(src_q > dest_max_q) then
       call utils_abort(myself//': q out of range', src_q, dest_max_q)
    end if
    if(src_l > dest_max_l .or. src_l < dest_min_l) then
       call utils_abort(myself//': l out range', src_l, dest_min_l, dest_max_l)
    end if
    minloffset = (dest_min_l - src_min_l)*(dest_min_l + src_min_l)*dest_max_q
    swx_xlat_sw_idx = src_m + (src_q-1) * (2*src_l+1) + (src_l-1) * (src_l+1) *&
         dest_max_q + src_l + dest_max_q + 1 !- minloffset

  end function swx_xlat_sw_idx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function swx_xlat_lm_idx(src_lm_idx, swri, src_swex, dest_swex_quality)
    !==========================================================================!
    ! Translates an LM index from src_swex to dest_swex. Both swexes must share!
    ! the same swri. An LM index indexes pairs (l,m), ie. it runs like a SW    !
    ! index, except it is contracted over q's.                                 !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2017.                                   !
    !==========================================================================!

    use sw_resolution_of_identity, only: SW_RI, SW_QUALITY
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    integer, intent(in)     :: src_lm_idx
    type(SW_RI), intent(in) :: swri
    type(SW_EX), intent(in) :: src_swex
    type(SW_QUALITY), intent(in) :: dest_swex_quality

    ! jd: Local variables
    integer :: src_max_q
    character(len=*), parameter :: myself = 'swx_xlat_lm_idx'

    ! -------------------------------------------------------------------------

    src_max_q = src_swex%quality%max_q

    if(src_max_q /= dest_swex_quality%max_q) then
       call utils_abort(myself//': max_q differs between src and dest swexes', &
            src_max_q, dest_swex_quality%max_q)
    end if

    swx_xlat_lm_idx = swx_xlat_sw_idx(src_lm_idx*src_max_q, swri, src_swex, &
         dest_swex_quality) / src_max_q

  end function swx_xlat_lm_idx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_merge_coeffs_pair(dest_coeffs, &
       src_coeffs1, src_coeffs2, swri, src_swex1, src_swex2,&
       factor1, factor2)
    !==========================================================================!
    ! Merges two sets of coeffs (src_coeffs1 from src_swex1, and src_coeffs2   !
    ! from swrc_swex2) into a set of coeffs 'dest_coeffs', by summing the      !
    ! coefficients. Optional weights for the sum can be specified via factor1  !
    ! and factor2.                                                             !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2017.                                   !
    !==========================================================================!

    use sw_resolution_of_identity, only: SW_RI, SW_QUALITY

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)          :: dest_coeffs(:)
    real(kind=DP), intent(in)           :: src_coeffs1(:)
    real(kind=DP), intent(in)           :: src_coeffs2(:)
    type(SW_RI), intent(in)             :: swri
    type(SW_EX), intent(in)             :: src_swex1
    type(SW_EX), intent(in)             :: src_swex2
    real(kind=DP), intent(in), optional :: factor1
    real(kind=DP), intent(in), optional :: factor2

    ! jd: Local variables
    integer          :: num_sws1, num_sws2
    integer          :: swex1_sw, swex2_sw, dest_swex_sw
    type(SW_QUALITY) :: dest_swex_quality
    real(kind=DP)    :: loc_factor1, loc_factor2

    ! -------------------------------------------------------------------------

    if(present(factor1)) then
       loc_factor1 = factor1
    else
       loc_factor1 = 1.0_DP
    end if
    if(present(factor2)) then
       loc_factor2 = factor2
    else
       loc_factor2 = 1.0_DP
    end if

    dest_swex_quality%min_l = min(src_swex1%quality%min_l, src_swex2%quality%min_l)
    dest_swex_quality%max_l = max(src_swex1%quality%max_l, src_swex2%quality%max_l)
    dest_swex_quality%max_q = max(src_swex1%quality%max_q, src_swex2%quality%max_q)

    num_sws1 = src_swex1%quality%num_sws_per_centre
    num_sws2 = src_swex2%quality%num_sws_per_centre

    do swex1_sw = 1, num_sws1
       dest_swex_sw = swx_xlat_sw_idx(swex1_sw, swri, src_swex1, dest_swex_quality)
       dest_coeffs(dest_swex_sw) = loc_factor1 * src_coeffs1(swex1_sw)
    end do
    do swex2_sw = 1, num_sws2
       dest_swex_sw = swx_xlat_sw_idx(swex2_sw, swri, src_swex2, dest_swex_quality)
       dest_coeffs(dest_swex_sw) = dest_coeffs(dest_swex_sw) + &
            loc_factor2 * src_coeffs2(swex2_sw)
    end do

  end subroutine swx_merge_coeffs_pair

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_merge_dcoeffs_pair(dest_dcoeffs, &
       src_dcoeffs1, src_dcoeffs2, swri, src_swex1, src_swex2, &
       factor1, factor2)
    !==========================================================================!
    ! Merges two sets of dcoeffs (src_dcoeffs1 from src_swex1, and src_coeffs2 !
    ! from swrc_swex2) into a set of dcoeffs 'dest_dcoeffs', by summing the    !
    ! d-coefficients. Optional weights for the sum can be specified via factor1!
    ! and factor2.                                                             !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2017.                                   !
    !==========================================================================!

    use sw_resolution_of_identity, only: SW_RI, SW_QUALITY

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)          :: dest_dcoeffs(:)
    real(kind=DP), intent(in)           :: src_dcoeffs1(:)
    real(kind=DP), intent(in)           :: src_dcoeffs2(:)
    type(SW_RI), intent(in)             :: swri
    type(SW_EX), intent(in)             :: src_swex1
    type(SW_EX), intent(in)             :: src_swex2
    real(kind=DP), intent(in), optional :: factor1
    real(kind=DP), intent(in), optional :: factor2

    ! jd: Local variables
    integer          :: num_lms1, num_lms2
    integer          :: swex1_lm, swex2_lm, dest_swex_lm
    type(SW_QUALITY) :: dest_swex_quality
    real(kind=DP)    :: loc_factor1, loc_factor2

    ! -------------------------------------------------------------------------

    if(present(factor1)) then
       loc_factor1 = factor1
    else
       loc_factor1 = 1.0_DP
    end if
    if(present(factor2)) then
       loc_factor2 = factor2
    else
       loc_factor2 = 1.0_DP
    end if

    dest_swex_quality%min_l = min(src_swex1%quality%min_l, src_swex2%quality%min_l)
    dest_swex_quality%max_l = max(src_swex1%quality%max_l, src_swex2%quality%max_l)
    dest_swex_quality%max_q = max(src_swex1%quality%max_q, src_swex2%quality%max_q)

    num_lms1 = (1+src_swex1%quality%max_l-src_swex1%quality%min_l) * &
         (1+src_swex1%quality%max_l+src_swex1%quality%min_l)
    num_lms2 = (1+src_swex2%quality%max_l-src_swex2%quality%min_l) * &
         (1+src_swex2%quality%max_l+src_swex2%quality%min_l)

    do swex1_lm = 1, num_lms1
       dest_swex_lm = swx_xlat_lm_idx(swex1_lm, swri, src_swex1, dest_swex_quality)
       dest_dcoeffs(dest_swex_lm) = loc_factor1 * src_dcoeffs1(swex1_lm)
    end do
    do swex2_lm = 1, num_lms2
       dest_swex_lm = swx_xlat_lm_idx(swex2_lm, swri, src_swex2, dest_swex_quality)
       dest_dcoeffs(dest_swex_lm) = dest_dcoeffs(dest_swex_lm) + &
            loc_factor2 * src_dcoeffs2(swex2_lm)
    end do

  end subroutine swx_merge_dcoeffs_pair

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function swx_coeffs_n_slots(natoms)
    !==========================================================================!
    ! Returns a reasonable value for the number slots in the coeffs HT.        !
    ! Calculations on large systems were used to find a good fit to how this   !
    ! changes with the number of atoms and MPI ranks.                          !
    ! Details do not matter much -- this HT is allowed to overfill, we just    !
    ! want to balance overheads with chain length.                             !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in April 2020.                                 !
    !==========================================================================!

    use comms, only: pub_total_num_procs

    implicit none

    ! jd: Arguments
    integer, intent(in) :: natoms

    ! -------------------------------------------------------------------------

    swx_coeffs_n_slots = &
         max(int(real(swx_slots_coeffs_per_atom * natoms,kind=DP) / &
         sqrt(real(pub_total_num_procs,kind=DP))),1)

  end function swx_coeffs_n_slots

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module sw_expansion_type
