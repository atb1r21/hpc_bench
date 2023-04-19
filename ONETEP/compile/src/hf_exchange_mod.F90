!================================================================!
!                                                                !
!                   Hartree-Fock exchange module                 !
!                                                                !
! This module contains subroutines for calculating a Hartree-Fock!
! exchange contribution to the energy using a resolution of      !
! identity approach.                                             !
!----------------------------------------------------------------!
! Groundwork laid by Quintin Hill in 2008/9 with supervision by  !
! Chris-Kriton Skylaris. Taken over by Jacek Dziedzic from 2011. !
! Several massive generalizations and expansions since.          !
!================================================================!
! References:                                                    !
! [1] Linear-scaling calculation of Hartree-Fock exchange energy !
!     with non-orthogonal generalised Wannier functions,         !
!     J. Dziedzic, Q. Hill, and C.-K. Skylaris,                  !
!     J. Chem. Phys. 139, 214103 (2013); doi: 10.1063/1.4832338. !
! [2] Massively parallel linear-scaling Hartree-Fock exchange    !
!     and hybrid exchange-correlation functionals with plane     !
!     wave basis set accuracy,                                   !
!     J. Dziedzic, J. C. Womack, R. Ali, and C.-K. Skylaris,     !
!     J. Chem. Phys. 155 224106 (2021).                          !
!================================================================!

!@TODO:
!
!Correct and clean code:
!- Clean up 4PI factor -- there is unnecessary * and / in places.
!- Get rid of old convention where NGWF radii could be different.
!  Look for @IDENTICAL_RADII

! ------------------------------------------------------------------------------
! [#] *** A very important note on neighbourhoods ***
! We call a pair of atoms 'S-neighbours' when their NGWFs overlap. This means
! the pair is within the sparsity pattern of the overlap matrix. This means
! that the spheres of radius r_NGWF centred on the two atoms intersect.
! Note that this does not necessarily mean that the PPD grids of the two atoms
! share a point or a PPD. The intersection may be so small that there are zero
! grid points within it.
!
! In HFx D atoms are S-neighbours of A atoms (and vice versa), and
! C atoms are S-neighbours of B atoms (and vice versa).
!
! The relation of being S-neighbours is commutative.
!
! Similarly, we call a pair of atoms 'X-neighbours' when they are not further
! away from one another than the HFx (X) cutoff. In HFx B atoms are X-neighbours
! of A atoms (and vice versa).
!
! The relation of being X-neighbours is commutative.
!
! Now, for the important point. Atom I is an XS-neighbour of atom J when there
! exists an atom K such that J and K are X-neighbours, and K and I are
! S-neighbours. In the context of HFx this always boils down to finding D atoms
! that are XS-neighbours of some B atom. To check for this, either examine the
! XS neighbour list of B and look for D's index, or:
!
!       if(neighbour_list_are_neighbours(global_b, global_d, xs_atoms_nl)) then
!          D is an XS-neighbour of B
!       end if
!
! It is imperative to realise that the XS relation is *NOT* commutative (!).
! This is easy to understand on a 1D example. Imagine r_NGWF = 7, r_X = 20.
! x_D = 0, x_A = 6, x_B = 20.
! D is an XS-neighbour of B, because D is within r_NGWF of A and A is within
! r_X of B. But B is *not* an XS-neighbour of D, because x_A is not within
! r_ngwf of B.
!
! Finally, SXS-neighbours are the atom pairs D-C. This one is commutative.
! ------------------------------------------------------------------------------

! ==========================
! HFx hash tables explained:
! ==========================
!
!-------------------------------------------------------------------------------
! Now residing in other modules:
!-------------------------------------------------------------------------------
!                  | # of | used in    | used in    | used in | independent of |
! Name             | keys | expansion? | calculate? | gradnt? |  DKN |  NGWFs  |
!-------------------------------------------------------------------------------
! coeffs           |   5  |    yes     |    yes     |    no   |  yes |   no    | Lives in SW_EX now, part of NGWF_REP.
! ngwfs_hts        |   1  |    yes     |    yes     |   yes   |  yes |   no    | Lives in remote now, as remote_ngwf_cache_hts.
! swops_in_ppds    |   3  |    yes     |    yes     |   yes   |  yes |  yes    | Lives in SW_RI now.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Residing in HFX_STATE:
!-------------------------------------------------------------------------------
!                  | # of | used in    | used in    | used in | independent of |
! Name             | keys | expansion? | calculate? | gradnt? |  DKN |  NGWFs  |
!-------------------------------------------------------------------------------
! expansions       |   3  |     no     |    yes     |    no   |  yes |   no    |
! ppd_products     |   3  |     no     |    yes     |    no   |  yes |   no    |
! ppd_ngwfs_dd     |   2  |     no     |     no     |   yes   |  yes |   no    |
! dlists           |   2  |     no     |    yes     |    no   |  yes |  yes    |
! alists           |   1  |     no     |    yes     |   yes   |  yes |  yes    |
! packed_metric... |   2  |    yes     |     no     |    no   |  yes |  yes    |
!-------------------------------------------------------------------------------
! NB: The expansions HT cannot be easily moved to NGWF_REP, where it would
!     perhaps belong, because it is modified in hfx_main(), where the rep is
!     intent(in). It cannot be pre-populated in hfx_expand_ngwf_pairs(), because
!     we might not be able to store all of it, and elements might have to be
!     recycled (expansions recomputed on the fly).
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Local:
!-----------------------------------------------------------------------------------------------------------
!                      | # of | used in | used in | used in | independent of |                             |
! Name                 | keys | expan?  | calc.?  | gradnt? |  DKN |  NGWFs  | Lives in                    |
!----------------------------------------------------------------------------------------------------------|
! packed_dkn_cd_blocks |   2  |   no    |   yes   |   yes   |   no |  yes    | hf_exchange_calculate()     |
! packed_dkn_ab_blocks |   2  |   no    |    no   |   yes   |   no |  yes    | hf_exchange_calculate()     |
! x_ab_contributions   |   2  |   no    |   yes   |    no   |   no |   no    | hfx_main()                  |
! grad_akets_dkn_loc   |   4  |   no    |    no   |   yes   |   no |   no    | hfx_main()                  |
! grad_akets_dkn_my    |   4  |   no    |    no   |   yes   |   no |   no    | hfx_main()                  |
! swop_hits            |   2  |   no    |    no   |    no   |  yes |  yes    | hfx_dkn_indep_preparation() |
!-----------------------------------------------------------------------------------------------------------
!
! *** coeffs_ht: ***
! Indexed by (Bb, Cc, coeffs_kind, swex_h, 1). The last key is always 1 and
! reflects the fact that the SW coefficients are DKN-agnostic, but formally
! DKN-dependent. 'swex_h' will always be REP_SWEX_HFX in HFx (in single NGWF basis
! mode) or REP_SWEX_HFX_OTHER (when mixing NGWF bases).
! The third key is either SW_O or SW_V, depending on metric.
! Stores coefficients of expansion of \phi_Bb * 'phi_Cc in terms of SWOPs.
! SWOPs are centred on B and C. This hash table is not allowed to evict elements.
! Size of 1 entry: num_sws_in_expansion. Typically: 500 doubles.
!
! *** ngwfs_hts *** (aka remote_ngwf_cache_hts(:)):
! Indexed by global NGWF index. Stores NGWFs in packed form. Used for caching
! remote NGWFs. This duplicates as a convenient device for getting the PPD
! indices (as opposed to actual ngwfs_on_grid) of remote NGWFs.
! Size of 1 entry: max_n_ppds_sphere * (cell%n_pts + 2) + 6.
! Different NGWF sets live in different indices of this HT-array.
! Typically: 20k doubles.

! *** swops_in_ppds_ht: ***
! Indexed by (src_atom, ppd, sw_or_swpot).
! Stores all SWOPs generated by src_atom, in a ppd.
! Size of 1 entry: num_sws * cell%n_pts. Typically: 60k doubles.
! SWOPs are stored with swriside quality. If swexside quality is different,
! filtering is done upon lookup.

! *** expansions_ht: ***
! Indexed by (ppd, Bb, Cc). Stores SWOP-expanded version of potential due to
! an NGWF pair (Bb, Cc), in a ppd. Bb and Cc are ordered (Bb<=Cc) (iff there
! is Bb-Cc symmetry, as happens when bases are not mixed).
! SWOPs are centred on B and C.
! Size of 1 entry: cell%n_pts. Typically: 125 doubles.
! NB: This cannot be moved to NGWF_REP -- we don't calculate all expansions in
!     advance, as there is typically not enough RAM to keep them all. They are
!     cached and purged on the fly in places where rep is not necessarily
!     in/out. Purging cannot happen in ngwf_representation, as this would lead
!     to a circular dependency between the two modules.
! NB: It might seem tempting to get rid of 'ppd' as the key and store the
!     expansion for all PPDs that are relevant to BC. This will not work because
!     of X truncation, at least in the case of valence calculations, where
!     symmetry dictates we only keep one expansion for BbCc and CcBb. Note that
!     the set of relevant PPDs differs between these two. The 1st set will
!     include PPDs of atoms that are X-neighbours with B, the 2nd set will
!     include PPDs of atoms that are X-neighbours with C.
! *** IMPORTANT ***
!     Note that we can get away with having only *one* hash table despite
!     multiple NGWF sets, because we *never* make subsequent calls to
!     hf_exchange_calculate(), with different NGWF-set-pairs, without an
!     intervening purge of expansions_ht, which happens in expand_ngwf_pairs().
!     That is, we always expand first, then make a call (in conduction) or
!     a number of calls (for LNV iterations, in valence), without the NGWFs
!     changing and *with the same NGWF-set-pairs* in the call. When we
!     switch to a new NGWF-set-pair, we always re-expand first.

! *** ppd_products_ht: ***
! Indexed by (A, D, ppd). Stores \phi_Aa * \phi_Dd for a given PPD.
! Uses a packed format. Includes all combinations of {a,d}.
! Size of 1 entry: num_ngwfs_a * num_ngwfs_d * cell%n_pts.

! *** ppd_ngwfs_dd_ht: ***
! Indexed by (D, ppd). Stores \phi_Dd for a given PPD.
! Uses a packed format. Includes all d's.
! Size of 1 entry: num_ngwfs_d * cell%n_pts.

! *** dlists_ht: ***
! Indexed by (B, ppd). Stores indices of D atoms that share this PPD *and*
! are XS-neighbours of this B. This makes it easier to iterate over relevant D.

! *** alists_ht: ***
! Indexed by (ppd). Stores indices of A atoms that share this PPD.

! *** packed_metric_matrix_blocks_ht: ***
! Indexed by (B,C). Stores metric matrix atomblocks that are my-local.

! *** packed_dkn_cd_blocks_ht: ***
! Indexed by (Cc, Dd). Stores DKN blocks. Used for caching remote DKN elements.
! Size of 1 entry: max_ngwfs_on_atom * max_ngwfs_on_atom * num_spins.
! Typically: 32 doubles.

! *** packed_dkn_ab_blocks_hts(3): ***
! Indexed by (Aa, Bb). Stores DKN blocks. Used for caching remote DKN elements
! for the NGWF gradient calculation. (1) stores K^AB, (2) stores tcK^A_B,
! (3) stores tcK^B_A -- tcK is not symmetric.
! Size of 1 entry: max_ngwfs_on_atom * max_ngwfs_on_atom * num_spins.
! Typically: 32 doubles.

! *** x_ab_contributions_ht: ***
! Indexed by (A,B). Stores contributions (blocks) to the X matrix from different
! MPI ranks.
! Size of 1 entry: num_ngwfs_on_atom * num_ngwfs_on_atom * pub_num_spins

! *** grad_akets_dkn_loc_ht: ***
! Indexed by (Aa, ppd, is, which_kernel). Like grad_akets_dkn_my_ht,
! but redistributed to sparse-local representation and summed over all B.
! Size of 1 entry: cell%n_pts.

! *** grad_akets_dkn_my_ht: ***
! Indexed by (Aa, ppd, is, which_kernel). NGWF gradient term 1 kets multiplied
! by K or tcK elements and summed over Bb. B only goes over my B atoms.
! Also includes a constant multiplicative factor. Term 2 kets are also
! converted to this form.
! Size of 1 entry: cell%n_pts.

! *** swop_hits_ht: ***
! Indexed by (PPD, src_atom). Trivial HT storing single elements -- number of
! times the (PPD, src_atom) pair is going to generate SWOPs in one
! NGWF iteration. This is only used for estimating which SWOPs are going to be
! used more heavily, and is short-lived.
! Size of 1 entry: 1 real.
!-------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

module hf_exchange

  use constants, only: DP, garbage_int
  use function_basis, only: FUNC_BASIS
  use hash_table, only: HT_HASH_TABLE
  use neighbour_list, only: NL_NEIGHBOUR_LIST
  use ngwf_representation, only: NGWF_REP, NGWF_HAM

  implicit none

  private

  public :: hf_exchange_init
  public :: hf_exchange_dkn_indep_stage
  public :: hf_exchange_calculate
  public :: hf_exchange_free

  public :: HFX_STATE

  type HFX_STATE

     ! ========================================================================
     ! State that is independent of everything (electrons and ions)
     ! ========================================================================

     ! jd: Handle to HFx's SW_RI
     integer            :: hfx_swri_h = -1
     logical            :: swri_is_shared_with_dma = .false.

     ! jd: Counters for coeff statistics
     real(kind=DP)      :: n_coeffs_processed = 0.0_DP ! jd: Currently
     real(kind=DP)      :: n_coeffs_elided = 0.0_DP    !     unused

     ! jd: ... and expansion HT and AD products HT statistics
     logical            :: expansions_feedback_printed = .false.
     logical            :: prods_feedback_printed = .false.

     ! jd: Control of debug files
     integer            :: hash_table_info_unit
     character(len=256) :: hash_table_info_filename


     ! ========================================================================
     ! State that is independent of NGWFs and DKN, only changing when ions move
     ! ========================================================================

     logical                     :: swop_cache_populated
     logical                     :: all_swops_cached
     logical                     :: all_expansions_cached

     ! jd: Neighbour lists
     type(NL_NEIGHBOUR_LIST)     :: s_atoms_nl   ! S-neighbours (A<->D, B<->C)
     type(NL_NEIGHBOUR_LIST)     :: x_atoms_nl   ! X-neighbours (A<->B)
     type(NL_NEIGHBOUR_LIST)     :: xs_atoms_nl  ! x-S-neighbours (B->D, A->C)
     type(NL_NEIGHBOUR_LIST)     :: sxs_atoms_nl ! s-x-S-neighbours (C<->D)

     ! jd: Indices of atoms and PPDs involved in exchange on this rank
     integer, allocatable        :: my_bc_pairs(:,:) ! BC pairs assigned to
     integer                     :: n_my_bc_pairs    ! this rank
     integer, allocatable        :: my_a_atoms(:)    ! X-neighbours of my B
     integer                     :: n_my_a_atoms
     integer, allocatable        :: my_b_atoms(:)    ! Unique B's in my BC pairs
     integer                     :: n_my_b_atoms
     integer, allocatable        :: my_c_atoms(:)    ! Unique C's in my BC pairs
     integer                     :: n_my_c_atoms
     integer, allocatable        :: my_d_atoms(:)    ! s-X-neighbours of my B
     integer                     :: n_my_d_atoms

     integer, allocatable        :: my_a_ppd_list(:) ! PPD indices of my A
     integer                     :: n_my_a_ppd_list

     integer, allocatable        :: my_bb_ngwf_offsets(:) ! idx by my_B
     integer                     :: n_my_bb_ngwfs    ! # of NGWFs on my B atoms

     ! jd: Lists of D's that share PPDs with A's of each B we own
     type(HT_HASH_TABLE)         :: dlists_ht

     ! jd: Lists of A's that share PPDs
     type(HT_HASH_TABLE)         :: alists_ht

     ! jd: Metric matrices distributed over my BC rather than local B
     type(HT_HASH_TABLE)         :: packed_metric_matrix_blocks_hts(2) ! V, O

     ! ========================================================================
     ! State that is independent of DKN, but not NGWFs.
     ! We want that to persist for the entire LNV loop.
     ! ========================================================================

     type(HT_HASH_TABLE)         :: expansions_ht
     type(HT_HASH_TABLE)         :: ppd_products_ht
     type(HT_HASH_TABLE)         :: ppd_ngwfs_dd_ht
     logical                     :: expansions_ht_in_use
     logical                     :: ppd_products_ht_in_use

     ! ========================================================================
     ! State that changes within a single call to hfx_calculate().
     ! ========================================================================

     ! jd: Nothing here

  end type HFX_STATE

  type HFX_NGWF_INDEX
     ! JCW: Derived type encapsulating pointers to FUNC_BASIS and NGWF_REP
     !      instances associated with an NGWF index in HFX evaluation. Using
     !      the type as a container for pointers allows array-of-pointers-like
     !      functionality in Fortran.
     type(FUNC_BASIS), pointer :: basis => NULL()
     type(NGWF_REP), pointer   :: rep => NULL()
  end type HFX_NGWF_INDEX

  ! jd: Module-wide parameters

  ! jd: Maximum number of atoms that can share a PPD. Used to dimension index
  !     arrays. Protected by an assert against overrun.
  integer, parameter :: max_num_at_per_ppd = 200

  ! jd: The below max* values do not really matter much, these HTs are
  !     permitted to overfill.
  integer, parameter :: maxslots_ppd_products_per_atom = 500
  ! ^ We need roughly 5000 entries per atom, at least for sytems in the 200-3000
  !   atom range (and it seems to be a constant). An 'entry' here is a PPD with
  !   ad products for a given AD pair (as opposed to the user feedback that
  !   counts AD pairs and so is about N_PPD_per_at ~ 50 times lower).
  !   Our strategy is to allocate 10x fewer slots to get a good balance between
  !   memory overheads and chain length. Remember that in this case we actually
  !   stop populating the HT once we hit the user-defined limit.
  integer, parameter :: maxslots_ppd_ngwfs_dd_per_atom = 200
  ! ^ This HT does not occupy a lot of space -- about 1.5 GB for largest doable
  !   systems (3000 atoms). Overheads are tiny. The above aims at providing
  !   a reasonably large number of slots so we don't waste time going along
  !   chains. This HT overfills if necessary.
  integer, parameter :: maxslots_x_ab_contributions = 10000
  integer, parameter :: max_x_ab_contributions = 20000
  ! ^ Entire HT is negligible
  integer, parameter :: maxslots_dkn_blocks_per_atom = 400
  integer, parameter :: max_cached_dkn_blocks_default = 800000
  ! ^ Typical requirements are about 550/atom on 32 MPI ranks, practically
  !   independent of system size (checked up to 3000 atoms). With fewer MPI
  !   ranks this increases sublinearly (so e.g. on 8 MPI ranks its 950/atom).
  !   However, overheads are significant here -- 400/atom leads to 100 MiB of
  !   overhead at 2000 atoms. Since this HT is permitted to overfill, we settle
  !   on 400/atom (unless user overrides with cache_limit_for_dknblks). But we
  !   never permit exceeding 800000, so the default then costs 100 MiB of
  !   overhead *max* and starts to chain a little if very few ranks are used.
  integer, parameter :: maxslots_dlists = 20000
  integer, parameter :: max_dlists = 20000
  ! ^ Negligible overheads
  integer, parameter :: maxslots_alists = 2000
  integer, parameter :: max_alists = 2000
  ! ^ Negligible overheads
  integer, parameter :: maxslots_swop_hits = 1000000
  integer, parameter :: max_swop_hits = 1000000
  ! ^ Roughly 130 MB, and it only lives briefly during initialisation.
  !   For 1132 atoms we need about 5M elements, so this strikes a good balance
  !   between overheads and chain length.
  integer, parameter :: maxslots_metric = 10000
  ! ^ Roughly 1MB of overhead.
  integer, parameter :: maxslots_grad_kets_dkn_loc_per_atom = 500
  ! ^ Roughly 0.06MB/atom of overhead
  integer, parameter :: maxslots_grad_kets_dkn_my_per_atom = 200
  ! ^ Roughly 0.025MB/atom of overhead

  character(len=32), parameter :: hash_table_info_rootname = 'hash_table_hfx'

  ! JCW: Integers for NGWF integral indices in arrays of NGWF_REP and
  !      FUNC_BASIS instances
  integer, parameter :: A_NGWF = 1 ! alpha
  integer, parameter :: B_NGWF = 2 ! beta
  integer, parameter :: C_NGWF = 3 ! gamma
  integer, parameter :: D_NGWF = 4 ! delta

contains

  ! ---------------------------------------------------------------------------
  ! ----- P U B L I C   S U B R O U T I N E S   A N D   F U N C T I O N S -----
  ! ---------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine hf_exchange_init(hfxstate, rep, mat, cell, natoms, init_rep_only)
    !==========================================================================!
    ! This subroutine initialises Hartree-Fock exchange.                       !
    ! - prints out banners,                                                    !
    ! - outputs .bib file,                                                     !
    ! - does sanity checks on set-up,                                          !
    ! - opens hash table log files,                                            !
    ! - initialises S, X, XS, and S-X-S neighbour lists,                       !
    ! - initialises hash tables:                                               !
    !   - EXPANSIONS,                                                          !
    !   - PPD_PRODUCTS,                                                        !
    !   - PPD_NGWFS_DD,                                                        !
    !   - DLISTS,                                                              !
    !   - ALISTS.                                                              !
    ! - initialises minor hfxstate details,                                    !
    ! - initialises SWEX containers for HFx.                                   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (inout): The HFX_STATE container to initialise.               !
    !   rep (inout): The NGWF_REP associated with HFx. Its overlap is used to  !
    !                generate the S-neighbour list. Its SW_EX structure will   !
    !                be initialised.                                           !
    !   mat (in): X matrix, the sparsity of which is used to generate the      !
    !             X-neighbour list. For HFx evaluation, this should be the     !
    !             exchange matrix of an NGWF_HAM.                              !
    !   cell (in): Needed to dimension caches operating on PPDs.               !
    !   natoms (in): Number of atoms in the HFx system (size of elements(:)).  !
    !   init_rep_only (in, opt): If present and .true., only the rep will be   !
    !                            initialised.                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 17/04/2012.                                 !
    ! Modified by Jacek Dziedzic in November 2013, April 2014.                 !
    ! Modified for SWRI/SWEX by Jacek Dziedzic in February 2015.               !
    ! Cleaned up by Jacek Dziedzic in October 2016.                            !
    ! Modified for hfxstate by Jacek Dziedzic in February 2019.                !
    ! Modified to handle subsystem calculations by Robert Charlton, 13/09/2018.!
    !==========================================================================!

    use bibliography, only: bibliography_cite
    use comms, only: pub_my_proc_id
    use constants, only: SIZEOF_DOUBLE, REP_SWEX_HFX, REP_SWEX_HFX_OTHER, &
         SW_V, SW_O
    use hash_table, only: hash_table_init, hash_table_calc_cache_size
    use neighbour_list, only: neighbour_list_init_from_sparse
    use rundat, only: pub_use_swx, pub_hfx_use_ri, pub_hfx_max_l, &
         pub_hfx_max_q, pub_hfx_metric, pub_cache_limit_for_expansions, &
         pub_cond_calculate, pub_use_activeswx, pub_use_activehfx, &
         pub_active_region, pub_dma_use_ri, pub_cache_limit_for_prods
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_create, sparse_destroy
    use sw_resolution_of_identity, only: swri_get_handle_to, swri_library
    use sw_expansion_type, only: swx_init_container
    use utils, only: utils_abort, utils_assert, utils_unit

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target :: hfxstate
    type(NGWF_REP), intent(inout)          :: rep
    type(SPAM3), intent(in)                :: mat
    type(CELL_INFO), intent(in)            :: cell
    integer, intent(in)                    :: natoms
    logical, intent(in), optional          :: init_rep_only

    ! jd: Local variables
    type(SPAM3) :: xs, sxs ! only for structures to get XS and SXS NL's
    integer     :: max_cached_expansions
    logical     :: filtering_needed
    logical     :: local_init_rep_only
    character(len=*), parameter :: myself = 'hf_exchange_init'

    ! -------------------------------------------------------------------------

    if(present(init_rep_only)) then
       local_init_rep_only = init_rep_only
    else
       local_init_rep_only = .false.
    end if

    hfxstate%swop_cache_populated = .false.
    hfxstate%all_swops_cached = .false.
    hfxstate%all_expansions_cached = .false.

    ! jd: Prepare SW_RI
    hfxstate%hfx_swri_h = swri_get_handle_to(pub_hfx_use_ri)

    ! jd: Check if SW_RI is shared with DMA. If so, we won't be destroying the
    !     full metric matrices later, once they are my-distributed, because DMA
    !     uses the old sparse-local distribution.
    if(trim(pub_dma_use_ri) /= '<unset>') then
       hfxstate%swri_is_shared_with_dma = &
            (swri_get_handle_to(pub_dma_use_ri) == hfxstate%hfx_swri_h)
    end if

    filtering_needed = .not. (&
         swri_library(hfxstate%hfx_swri_h)%quality%max_l == pub_hfx_max_l .and.&
         swri_library(hfxstate%hfx_swri_h)%quality%max_q == pub_hfx_max_q)

    if(filtering_needed) then
       call utils_abort('This version of Hartree-Fock exchange does not support&
            & swri-to-swex filtering. That is, the quality of the SW resolution&
            & of identity (lmax and qmax specified in the block selected using&
            & hfx_use_ri) must match the actual quality of the SW expansion&
            & (hfx_max_l and hfx_max_q). This is done for the sake of&
            & efficiency.')
    end if

    if(.not. local_init_rep_only) then

       call bibliography_cite('HFX')

       ! rc2013: check if this is an active subsystem calculation
       if(pub_use_activehfx) then
          call utils_assert(pub_use_activeswx, myself//': Spherical wave&
               & expansion (pub_use_activeswx) must be turned on for &
               & subsystem Hartree-Fock exchange.')
       else
          call utils_assert(pub_use_swx, myself//': Spherical wave expansion&
               & (pub_use_swx) must be turned on for Hartree-Fock exchange.')
       end if

       if(trim(pub_hfx_use_ri) == trim('<unset>')) then
          call utils_abort(myself//&
               ': ''hfx_use_ri'' is mandatory when using Hartree-Fock exchange.')
       end if

       if(pub_hfx_max_l == -1) then
          call utils_abort(myself//&
               ': ''hfx_max_l'' is mandatory when using Hartree-Fock exchange.')
       end if

       if(pub_hfx_max_q == -1) then
          call utils_abort(myself//&
               ': ''hfx_max_q'' is mandatory when using Hartree-Fock exchange.')
       end if

       ! jd: Open a file for debugging hash-tables
       hfxstate%hash_table_info_unit = utils_unit()
       write(hfxstate%hash_table_info_filename,'(a,a,a,a,i0,a)') &
            trim(hash_table_info_rootname), '_', '', &
            'proc_', pub_my_proc_id, '.log'
       open(hfxstate%hash_table_info_unit, &
            file=hfxstate%hash_table_info_filename, err=10)

       ! jd: Initialise neighbour lists: S, X
       ! rc2013: using the active subsystem overlap matrix
       ! -- by default pub_active_region = 1
       call neighbour_list_init_from_sparse(hfxstate%s_atoms_nl, &
            'HFx_S_ATOMS', rep%overlap%m(pub_active_region,pub_active_region))
       call neighbour_list_init_from_sparse(hfxstate%x_atoms_nl, &
            'HFx_X_ATOMS', mat)
       ! jd: Initialise neighbour lists:
       !     - XS (between B and D, B -> D),
       !     - SXS (between C and D).
       !     - SX (between D and B, D -> B),
       call sparse_create(xs, &
            rep%overlap%m(pub_active_region,pub_active_region), &
            mat, rep%overlap%m(pub_active_region,pub_active_region)%iscmplx)
       call neighbour_list_init_from_sparse(hfxstate%xs_atoms_nl, &
            'HFx_XS_ATOMS', xs)
       call sparse_create(sxs, xs, &
            rep%overlap%m(pub_active_region,pub_active_region), &
            rep%overlap%m(pub_active_region,pub_active_region)%iscmplx)
       call neighbour_list_init_from_sparse(hfxstate%sxs_atoms_nl, &
            'HFx_SXS_ATOMS', sxs)
       call sparse_destroy(sxs)
       call sparse_destroy(xs)
       ! jd: @IDENTICAL_RADII between val and cond implicitly assumed
       !     in that only the val call populates the NLs.

       ! jd: Prepare expansions hash table
       max_cached_expansions = hash_table_calc_cache_size(&
            pub_cache_limit_for_expansions, cell%n_pts * SIZEOF_DOUBLE)
       hfxstate%expansions_ht_in_use = .not. (max_cached_expansions == 0)
       if(hfxstate%expansions_ht_in_use) then
          call hash_table_init(hfxstate%expansions_ht, 'EXPANSIONS', 3, &
               max_cached_expansions/5, max_cached_expansions, &
               hfxstate%hash_table_info_unit) ! /5 to reduce overhd of procs
       end if

       ! jd: Prepare PPD products hash table
       hfxstate%ppd_products_ht_in_use = &
            .not. (pub_cache_limit_for_prods == 0)
       if(hfxstate%ppd_products_ht_in_use) then
          call hash_table_init(hfxstate%ppd_products_ht, 'PPD_PRODUCTS', 3, &
               maxslots_ppd_products_per_atom * natoms, &
               maxslots_ppd_products_per_atom * natoms, &
               hfxstate%hash_table_info_unit)
       end if

       ! jd: Prepare PPD NGWFs Dd hash table
       call hash_table_init(hfxstate%ppd_ngwfs_dd_ht, 'PPD_NGWFS_DD', 2, &
            maxslots_ppd_ngwfs_dd_per_atom * natoms, &
            maxslots_ppd_ngwfs_dd_per_atom * natoms, &
            hfxstate%hash_table_info_unit)

       ! jd: Prepare DLISTS hash table
       call hash_table_init(hfxstate%dlists_ht, 'DLISTS', 2, &
            maxslots_dlists, max_dlists, hfxstate%hash_table_info_unit)

       ! jd: Prepare ALISTS hash table
       call hash_table_init(hfxstate%alists_ht, 'ALISTS', 1, &
            maxslots_alists, max_alists, hfxstate%hash_table_info_unit)

       ! jd: Prepare METRIC MATRIX hash table
       if(swri_library(hfxstate%hfx_swri_h)%has_metric(SW_V)) then
          call hash_table_init(hfxstate%packed_metric_matrix_blocks_hts(1), &
               'VMETRIC', 2, maxslots_metric, maxslots_metric, &
               hfxstate%hash_table_info_unit)
       end if
       if(swri_library(hfxstate%hfx_swri_h)%has_metric(SW_O)) then
          call hash_table_init(hfxstate%packed_metric_matrix_blocks_hts(2), &
               'OMETRIC', 2, maxslots_metric, maxslots_metric, &
               hfxstate%hash_table_info_unit)
       end if

    end if

    ! jd: Prepare SW_EX
    call swx_init_container(rep%swexes(REP_SWEX_HFX), &
         swri_library(hfxstate%hfx_swri_h), 'hfx_swex', &
         0, pub_hfx_max_l, pub_hfx_max_q, pub_hfx_metric, rep%postfix)

    ! jd: Prepare SW_EX for conduction, if needed
    if(pub_cond_calculate) then
       call swx_init_container(rep%swexes(REP_SWEX_HFX_OTHER), &
            swri_library(hfxstate%hfx_swri_h), 'hfx_swex_other', &
            0, pub_hfx_max_l, pub_hfx_max_q, pub_hfx_metric, rep%postfix)
    end if

    return

10  call utils_abort(myself//': Could not create file: '//&
         trim(hfxstate%hash_table_info_filename))

  end subroutine hf_exchange_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_free(hfxstate)
    !==========================================================================!
    ! Cleans up after Hartree-Fock exchange.                                   !
    ! Currently this amounts to freeing the neighbour lists, the hash tables   !
    ! (EXPANSIONS and PPD_PRODUCTS), allocatables in hfxstate and closing the  !
    ! hash table log file.                                                     !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 18/06/2012.                                 !
    ! Modified for hfxstate by Jacek Dziedzic in February 2019.                !
    !==========================================================================!

    use constants, only: SW_V, SW_O
    use hash_table, only: hash_table_free
    use neighbour_list, only: neighbour_list_free
    use sw_resolution_of_identity, only: swri_library
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target :: hfxstate

    ! jd: Local variables
    character(len=*), parameter :: myself = 'hf_exchange_free'

    ! -------------------------------------------------------------------------

    ! jd: Free remaining HFX_STATE caches
    call hash_table_free(hfxstate%alists_ht)
    call hash_table_free(hfxstate%dlists_ht)
    call hash_table_free(hfxstate%ppd_ngwfs_dd_ht)
    if(hfxstate%ppd_products_ht_in_use) then
       call hash_table_free(hfxstate%ppd_products_ht)
    end if
    if(hfxstate%expansions_ht_in_use) then
       call hash_table_free(hfxstate%expansions_ht)
    end if
    if(swri_library(hfxstate%hfx_swri_h)%has_metric(SW_V)) then
       call hash_table_free(hfxstate%packed_metric_matrix_blocks_hts(1))
    end if
    if(swri_library(hfxstate%hfx_swri_h)%has_metric(SW_O)) then
       call hash_table_free(hfxstate%packed_metric_matrix_blocks_hts(2))
    end if
    close(hfxstate%hash_table_info_unit, err=10)

    ! jd: Free persistent neighbour lists
    call neighbour_list_free(hfxstate%sxs_atoms_nl)
    call neighbour_list_free(hfxstate%xs_atoms_nl)
    call neighbour_list_free(hfxstate%x_atoms_nl)
    call neighbour_list_free(hfxstate%s_atoms_nl)

    call hfx_free_arrays(hfxstate)

    return

10  call utils_abort(myself//': Could not close file: '//&
         hfxstate%hash_table_info_filename)

  end subroutine hf_exchange_free

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_dkn_indep_stage(hfxstate, mdl, ireg, par, &
       rep, ngwf_basis, rep2, ngwf_basis2, basis_selector)
    !==========================================================================!
    ! Performs the DKN-independent stage of HFx:                               !
    ! (1) Expands NGWF pairs, storing the coefficients in rep (or in rep2, if  !
    !     basis_selector is present and says so). Since this means NGWFs       !
    !     changed, this also invalidates all NGWF-dependent state.             !
    ! (2) Does preparations for the actual HFx calculation that can be done    !
    !     once per every NGWF change (before the LNV loop). Some of these      !
    !     do not even depend on NGWFs, only ionic positions, but these are     !
    !     cheap and we can re-do them in the NGWF SCF loop. For details see    !
    !     hfx_dkn_indep_preparation().                                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (inout): The results of the preparations are stored here.     !
    !   mdl (in): The model. Needed for ionic positions, cell, regions.        !
    !   ireg (in): Region ID for embedding.                                    !
    !   par (in): Needed to figure out who rank-owns which atoms.              !
    !   rep (inout): The coefficients of NGWF pair expansions go here unless   !
    !                basis_selector(5) says otherwise.                         !
    !   ngwf_basis (in): The first NGWF basis.                                 !
    !                                                                          !
    ! Optional arguments for mixed-basis calculations:                         !
    !   rep2, ngwf_basis2, basis_selector.                                     !
    !   - See banner for hf_exchange_calculate().                              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2019.                              !
    !==========================================================================!

    use model_type, only: MODEL
    use parallel_strategy, only: PARAL_INFO
    use sw_resolution_of_identity, only: swri_library

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target         :: hfxstate
    type(NGWF_REP), intent(inout), target          :: rep
    type(FUNC_BASIS), intent(in), target           :: ngwf_basis
    type(MODEL), intent(in)                        :: mdl
    integer, intent(in)                            :: ireg
    type(PARAL_INFO), intent(in)                   :: par
    ! jd: Arguments for mixed NGWF bases
    type(NGWF_REP), optional, intent(in), target   :: rep2
    type(FUNC_BASIS), optional, intent(in), target :: ngwf_basis2
    integer, optional, intent(in)                  :: basis_selector(5)

    ! -------------------------------------------------------------------------

    ! jd: Invalidate all state that depends on NGWFs (expansions, NGWF products,
    !     Dd NGWFs).
    call hfx_invalidate_ngwf_dependent_state(hfxstate)

    ! jd: Do all the book-keeping and comms that are DKN-independent.
    !     Populate SWOP caches.
    call hfx_dkn_indep_preparation(hfxstate, &
         swri_library(hfxstate%hfx_swri_h), ireg, par, mdl, &
         rep, ngwf_basis, rep2, ngwf_basis2, basis_selector)

    ! jd: Expand NGWF pairs
    call hfx_expand_ngwf_pairs(hfxstate, rep, ngwf_basis, mdl, ireg, &
         rep2, ngwf_basis2, basis_selector)

  end subroutine hf_exchange_dkn_indep_stage

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_calculate(hfxstate, mat, &  ! in/out
       rep, ireg, denskern_ab, denskern_cd, ngwf_basis, fftbox, cell, elements, & ! in
       calcgradient, &                               ! input
       tc_denskern_ab, cov_grad, &                   ! opt args for gradient
       contra_grad, precond_func_recip, &            ! opt args for gradient
       rep2, ngwf_basis2, basis_selector, &          ! opt args for mixed NGWFs
       energy_prefactor, grad_prefactor, &           ! opt args for mixed NGWFs
       symmetrise_xmatrix)                           ! opt arg
    !==========================================================================!
    ! Calculates the Hartree-Fock exchange matrix and optionally the exchange  !
    ! contribution to the NGWF gradient. It also calculates the HFx energy,    !
    ! but this is only used for diagnostics -- it's printed out, but never     !
    ! referenced (instead it comes form the X matrix contribution to H).       !
    !--------------------------------------------------------------------------!
    ! Requirements:                                                            !
    ! - hf_exchange_dkn_indep_stage() must have been called first.             !
    !--------------------------------------------------------------------------!
    ! For notation (meaning of abcd indices) see [1].                          !
    ! hfxstate (inout): The HFX_STATE container that stores all HFx data that  !
    !                   need to persist across calls.                          !
    ! mat (inout): The SPAM3_EMBED matrix into which the resulting matrix (X)  !
    !              will be placed. For evaluating HFX, this should be the      !
    !              hfexchange component of an NGWF_HAM instance.               !
    ! rep (in):        } Describe the NGWF basis.                              !
    ! ngwf_basis (in): }                                                       !
    ! ireg (in): Region identifier for embedding.                              !
    ! denskern_ab (in): The density kernel, or in general, some contravariant  !
    !                   tensor with which the X matrix is contracted to obtain !
    !                   the HFx energy. In conduction calcs: the P matrix.     !
    !                   In valence calcs: the usual DKN.                       !
    ! denskern_cd (in): The density kernel, or in general, some contravariant  !
    !                   tensor that is part of the exchange operator. In both  !
    !                   the valence and conduction calcs: the usual DKN.       !
    ! fftbox (in): The usual. Needed for shaving NGWFs, grid weight and for    !
    !              KE preconditioning.                                         !
    ! cell (in): The usual. Needed for n_pts and weight.                       !
    ! elements (in): The usual.                                                !
    ! calcgradient (in): Pass .true. is HFx NGWF gradient is to be calculated. !
    !                                                                          !
    ! *** Optional arguments for NGWF gradient calculation, needed only when   !
    ! calcgradient is .true.:                                                  !
    !                                                                          !
    ! tc_denskern_ab (in): Tensorially-correct version of denskern_ab.         !
    ! cov_grad (in/out): Covariant NGWF gradient, HFx contribution will        !
    !                    be *added* here.                                      !
    ! contra_grad (in/out): Contrvariant NGWF gradient, HFx contribution will  !
    !                    be *added* here.                                      !
    ! precond_func_recip (in): Needed for gradient preconditioning.            !
    !                                                                          !
    ! *** Optional arguments for mixed NGWF bases. Either pass none, or all:   !
    !                                                                          !
    ! rep2 (in):        } Describe the 2nd NGWF basis.                         !
    ! ngwf_basis2 (in): }                                                      !
    ! basis_selector (in): Chooses the NGWF bases for the calculation of       !
    !                      TEFCIs. Quintin's alpha-beta-gamma-delta convention !
    !                      convention is used (cf. [1], eq. 5). Alpha is in the!
    !                      bra, beta is in the ket. These two are the indices  !
    !                      of the exchange matrix. Delta is in the bra,        !
    !                      gamma  is in the ket. K^{\delta,\gamma}, or whatever!
    !                      is passed for denskern_cd is contracted over when   !
    !                      building the X matrix. This is different from Tim's !
    !                      convention in JCP 2013, and different from the CCP9 !
    !                      Flagship meeting convention.                        !
    !                      basis_selector(1) chooses basis for alpha (1 or 2). !
    !                      basis_selector(2) chooses basis for beta (1 or 2).  !
    !                      basis_selector(3) chooses basis for gamma (1 or 2). !
    !                      basis_selector(4) chooses basis for delta (1 or 2). !
    !                      basis_selector(5) chooses the *source* swex in rep, !
    !                                        from which the expansion coeffs   !
    !                                        will be retrieved.                !
    !                      (/1,1,1,1,REP_SWEX_HFX/) is for valence.            !
    !                      (/1,1,2,2,REP_SWEX_HFX_OTHER/) is for conduction,   !
    !                      where rep = cond_rep, rep2 = val_rep, so cond is 1, !
    !                      val is 2.                                           !
    !                                                                          !
    ! *** Optional argument for pre-multiplying the exchange matrix:           !
    ! energy_prefactor (in): Defaults to 1.0. In conduction a different value  !
    !                        is passed to account for the different way that   !
    !                        conduction handles spin degeneracy.               !
    ! *** Optional argument for pre-multiplying the NGWF gradient:             !
    ! grad_prefactor (in): Defaults to 1.0. In conduction 0.5 can be passed to !
    !                      take into account the fact that only one set of     !
    !                      NGWFs is optimised, and the other one is fixed      !
    !                      (I think).                                          !
    !                      This also takes care of spin degeneracy conventions !
    !                      in conduction (see the call to hf_exchange_calculate!
    !                      from conduction).                                   !
    ! *** Optional argument for controlling symmetrisation                     !
    ! symmetrise_xmatrix (in) : If present and false, then disable             !
    !                           symmetrising the X-matrix by averaging with    !
    !                           its transpose. In this case, the X matrix      !
    !                           returned will be approximated in an unbalanced !
    !                           way, with only the ket (beta-gamma) expanded   !
    !                           in SWs.                                        !
    ! The X matrix is symmetrised by default by taking the average of the      !
    ! matrix with its transpose. This is equivalent to averaging over matrix   !
    ! elements where either the bra or ket NGWF products are fitted using the  !
    ! local 2-centre scheme described in II.D.3 of [1]. This relies on the     !
    ! symmetry of the 4-index electron repulsion integrals when the alpha-beta !
    ! and gamma-delta NGWF pairs share the same NGWF basis.                    !
    !                                                                          !
    ! When the alpha and beta NGWF bases differ (as occurs when evaluating     !
    ! the HF-contribution to V_SCF in LR-TDDFT), symmetrisation by averaging   !
    ! using the transpose of X is no longer valid and symmetrise_xmatrix       !
    ! should be .false.                                                        !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in 2008 and January 2009.                        !
    ! Modified by Quintin Hill on 12/02/2009 to make hfexchange symmetric.     !
    ! Rewritten by Jacek Dziedzic in July 2012 to use 2-centre 'BC' expansion. !
    ! Split by Jacek Dziedzic in October 2016 to separate NGWF expansion from  !
    ! the calculation of exchange.                                             !
    ! Extended by Jacek Dziedzic in June 2018 to support mixed NGWF bases.     !
    ! Modified by Jacek Dziedzic in Feb-Mar 2019 for new parallel scheme.      !
    ! Generalized by James C. Womack in Jul 2019 to allow evaluation of        !
    ! HF-contribution of SCF response potential matrix.                        !
    !==========================================================================!

    use constants, only: stdout, REP_SWEX_HFX
    use datatypes, only: FUNCTIONS
    use comms, only: pub_on_root, comms_barrier
    use constants, only: VERBOSE, PROLIX, CRLF
    use fft_box, only: FFTBOX_INFO
    use hash_table, only: HT_HASH_TABLE, hash_table_init, hash_table_free, &
         hash_table_list
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use remote, only: remote_ngwf_cache_hts
    use rundat, only: pub_num_spins, pub_hhf_factor, pub_hhf_nstates, &
         pub_hfx_output_detail, pub_precond_recip
    use sparse, only: SPAM3, sparse_scale, sparse_create, sparse_copy, &
         sparse_destroy
    use sparse_embed, only: SPAM3_EMBED
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use sparse, only: SPAM3, sparse_show_matrix, sparse_trace, sparse_get_par
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check, &
         utils_postfix_to_ngwf_set_name, utils_trace_in, utils_trace_out, &
         utils_abort

    implicit none

    ! jd: Arguments
    ! hfxstate cannot be intent(in), hfx_main() modifies expansions
    type(HFX_STATE), intent(inout), target :: hfxstate
    type(SPAM3_EMBED), intent(inout)       :: mat(:)
    type(NGWF_REP), intent(in), target     :: rep
    type(SPAM3), intent(in)                :: denskern_ab(pub_num_spins) ! left K
    type(SPAM3), intent(in)                :: denskern_cd(pub_num_spins) ! right K
    type(FUNC_BASIS), intent(in), target   :: ngwf_basis
    type(FFTBOX_INFO), intent(in)          :: fftbox
    type(CELL_INFO), intent(in)            :: cell
    type(ELEMENT), intent(in)              :: elements(:)
    logical, intent(in)                    :: calcgradient ! Calculate NGWF grad?
    integer, intent(in)                    :: ireg  ! region
    ! qoh: Arguments for calculating NGWF gradient
    ! qoh: (required if calcgradient is true)
    type(SPAM3), optional, intent(in)        :: tc_denskern_ab(pub_num_spins)
    type(FUNCTIONS), optional, intent(inout) :: cov_grad
    type(FUNCTIONS), optional, intent(inout) :: contra_grad
    real(kind=DP), optional, intent(in)      :: precond_func_recip(:,:,:)
    ! jd: Arguments for mixed NGWF bases
    type(NGWF_REP), optional, intent(in), target   :: rep2
    type(FUNC_BASIS), optional, intent(in), target :: ngwf_basis2
    integer, optional, intent(in)                  :: basis_selector(5)
    real(kind=DP), optional, intent(in)            :: energy_prefactor
    real(kind=DP), optional, intent(in)            :: grad_prefactor
    ! JCW: Control symmetrisation
    logical, optional, intent(in) :: symmetrise_xmatrix

    ! jd: Local variables
    integer                     :: src_swex_h
    type(HT_HASH_TABLE), target :: packed_dkn_cd_blocks_ht
    type(HT_HASH_TABLE), target :: packed_dkn_ab_blocks_hts(3) ! K^AB, tcK^A_B, tcK^B_A
    integer                     :: max_cached_dkn_blocks_cd
    integer                     :: max_cached_dkn_blocks_ab
    integer                     :: is
    integer                     :: natoms
    integer                     :: i_ht
    real(kind=DP)               :: hfx_energy_for_spin
    real(kind=DP)               :: hfx_energy
    real(kind=DP)               :: trKS_ab, trKS_cd

    ! jd: Local variables for mixed NGWF basis sets
    type(HFX_NGWF_INDEX) :: ngwf_indices(A_NGWF:D_NGWF)
    ! agrecocmplx
    logical :: loc_cmplx
    logical :: loc_symmetrise_xmatrix

    ! jcap: temporary arrays for hf_exchange
    type(SPAM3), allocatable :: xmatrix(:)
    type(PARAL_INFO), pointer :: par ! rc2013: parallel strategy
    integer :: ierr

    ! JCW: Pointers to NGWF overlap matrices (within rep_ab, rep_cd) that
    !      correspond to the density kernels denskern_ab, denskern_cd
    type(SPAM3), pointer :: overlap_ab
    type(SPAM3), pointer :: overlap_cd
    logical :: x_overlap_rep(A_NGWF:D_NGWF) ! used to identify whether overlap
                                            ! is available

    character(len=*), parameter :: myself = 'hf_exchange_calculate'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    ! JCW: Determine whether to symmetrise X matrix by averaging with its
    !      transpose or not
    if (present(symmetrise_xmatrix)) then
       ! Symmetrise X by averaging with transpose if optional argument is true
       loc_symmetrise_xmatrix = symmetrise_xmatrix
    else
       ! Default: symmetrise X by averaging with transpose
       loc_symmetrise_xmatrix = .true.
    end if

    ! rc2013: get the parallel strategy from the diagonal matrix
    call sparse_get_par(par, denskern_ab(1))
    call utils_assert(par%nat == size(elements), 'Error in '//myself//': &
         &allocated parallel strategy is incompatible with elements.')

    natoms = size(elements)

    allocate(xmatrix(pub_num_spins),stat=ierr)
    call utils_alloc_check(myself,'xmatrix',ierr)

    do is=1,pub_num_spins
       call sparse_create(xmatrix(is),mat(is)%m(ireg,ireg))
       call sparse_copy(xmatrix(is),mat(is)%m(ireg,ireg))
    end do


    if(pub_on_root .and. pub_hfx_output_detail >= VERBOSE) then
       if(.not. present(rep2)) then
          write(stdout,'(a,a,a)') 'HFx: Calculation for ', &
               trim(utils_postfix_to_ngwf_set_name(rep%postfix)), &
               ' (single NGWF basis).'
       else
          write(stdout,'(a,a,a,a,a,a,a,i0,i0,i0,i0,i0,a)') 'HFx: Calculation for ', &
               trim(utils_postfix_to_ngwf_set_name(rep%postfix)), '+', &
               trim(utils_postfix_to_ngwf_set_name(rep2%postfix)), &
               ', structure: ', trim(xmatrix(1)%structure), ', selector: ', &
               basis_selector, ' (mixed NGWF basis).'
       end if
    end if

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(ireg)%iscmplx

    ! agrecocmplx
    call utils_assert(.not. loc_cmplx, 'Error in&
         & hf_exchange_calculate: not ready yet for complex NGWFs.')

    ! jd: Ensure consistency between optional parameters for mixed NGWF bases
    call utils_assert((present(rep2) .eqv. present(ngwf_basis2)) .and. &
         (present(rep2) .eqv. present(basis_selector)), &
         myself//': Optional arguments for HFx with mixed NGWF bases must &
         &either all be present or all be absent.')

    ! jd: Ensure consistency between optional parameters for gradient
    if(calcgradient) then
       call utils_assert(present(tc_denskern_ab) .and. present(cov_grad) .and. &
            present(contra_grad) .and. &
            (present(precond_func_recip) .or. .not. pub_precond_recip), &
            myself//': When calcgradient is .true., certain optional arguments&
            & must be passed.')
    end if
    call utils_assert((present(rep2) .eqv. present(ngwf_basis2)) .and. &
         (present(rep2) .eqv. present(basis_selector)), &
         myself//': Optional arguments for HFx with mixed NGWF bases must &
         &either all be present or all be absent.')

    ! jd: Select the right source NGWF_REP (whose SW_EX holds the coeffs),
    !     and right NGWF bases.
    if(present(basis_selector)) then
       !jmecmplx
       call utils_assert(.not.rep2%ngwfs_on_grid(ireg)%iscmplx, myself//&
            ': not ready for complex NGWFs yet (rep2).')
       call utils_assert(rep%ngwfs_on_grid(ireg)%iscmplx .eqv. &
            rep2%ngwfs_on_grid(ireg)%iscmplx, myself//&
            ': Cannot mix real and complex NGWFs between REPs.')
       !call utils_assert(basis_selector(1) == basis_selector(2), myself//&
       !     ': Mixing bases *within* a density kernel (AB) not supported yet.')
       if (.not.basis_selector(1) == basis_selector(2).and.pub_on_root) then
          write(stdout,'(a)') "HFx WARNING: mixed basis density kernel (AB) --&
               & experimental functionality"
       end if
       !call utils_assert(basis_selector(3) == basis_selector(4), myself//&
       !     ': Mixing bases *within* a density kernel (CD) not supported yet.')
       if (.not.basis_selector(3) == basis_selector(4).and.pub_on_root) then
          write(stdout,'(a)') "HFx WARNING: mixed basis density kernel (CD) --&
               & experimental functionality"
       end if

       src_swex_h = basis_selector(5)

       if(basis_selector(1) == 1) then
          ngwf_indices(A_NGWF)%basis => ngwf_basis
          ngwf_indices(A_NGWF)%rep => rep
       else if(basis_selector(1) == 2) then
          ngwf_indices(A_NGWF)%basis => ngwf_basis2
          ngwf_indices(A_NGWF)%rep => rep2
       else
          call utils_abort(myself//': Illegal basis_selector(1)')
       end if
       if(basis_selector(2) == 1) then
          ngwf_indices(B_NGWF)%basis => ngwf_basis
          ngwf_indices(B_NGWF)%rep => rep
       else if(basis_selector(2) == 2) then
          ngwf_indices(B_NGWF)%basis => ngwf_basis2
          ngwf_indices(B_NGWF)%rep => rep2
       else
          call utils_abort(myself//': Illegal basis_selector(2)')
       end if
       if(basis_selector(3) == 1) then
          ngwf_indices(C_NGWF)%basis => ngwf_basis
          ngwf_indices(C_NGWF)%rep => rep
       else if(basis_selector(3) == 2) then
          ngwf_indices(C_NGWF)%basis => ngwf_basis2
          ngwf_indices(C_NGWF)%rep => rep2
       else
          call utils_abort(myself//': Illegal basis_selector(3)')
       end if
       if(basis_selector(4) == 1) then
          ngwf_indices(D_NGWF)%basis => ngwf_basis
          ngwf_indices(D_NGWF)%rep => rep
       else if(basis_selector(4) == 2) then
          ngwf_indices(D_NGWF)%basis => ngwf_basis2
          ngwf_indices(D_NGWF)%rep => rep2
       else
          call utils_abort(myself//': Illegal basis_selector(4)')
       end if

       ! JCW: Select the overlap matrices associated with the NGWF sets
       !      for the pairs of indices (alpha-beta, gamma-delta)
       ! JCW: rep%cross_overlap is only sparse_embed_create'd by
       !      ngwf_rep_create when rep%postfix is 'c', 'j', 'a' or 'z'
       x_overlap_rep(A_NGWF) = &
            any(ngwf_indices(A_NGWF)%rep%postfix == [ 'c', 'j', 'a', 'z'])
       x_overlap_rep(B_NGWF) = &
            any(ngwf_indices(B_NGWF)%rep%postfix == [ 'c', 'j', 'a', 'z'])
       x_overlap_rep(C_NGWF) = &
            any(ngwf_indices(C_NGWF)%rep%postfix == [ 'c', 'j', 'a', 'z'])
       x_overlap_rep(D_NGWF) = &
            any(ngwf_indices(D_NGWF)%rep%postfix == [ 'c', 'j', 'a', 'z'])
       if (.not.basis_selector(1) == basis_selector(2)) then
          ! JCW: Expect NGWF overlap matrix to be cross_overlap when
          !      NGWF bases for alpha and beta differ and one of the bases
          !      is valence

          ! JCW: cross_overlap should only be available for cond or joint
          !      NGWF representations, so find which of the reps for alpha
          !      and beta NGWF indices are cond or joint
          ! JCW: cross_overlap evaluated in hamiltonian_dens_indep_matrices is
          !      <val|cond> and denskern_ab is expected to be P^{cond,val}
          if (all(x_overlap_rep(A_NGWF:B_NGWF).eqv.[.true.,.true.])) then
             ! Both are non-val, overlap (probably) not available
             overlap_ab => NULL()
          else if (all(x_overlap_rep(A_NGWF:B_NGWF).eqv.[.true.,.false.])) then
             ! A is non-val, expect <val|cond> cross overlap and P^{cond,val}
             overlap_ab => ngwf_indices(A_NGWF)%rep%cross_overlap%p
          else if (all(x_overlap_rep(A_NGWF:B_NGWF).eqv.[.false.,.true.])) then
             ! B is non-val, expect <cond|val> cross overlap and P^{val,cond}
             ! but cross_overlap_tr does not seem to be evaluated by default
             !overlap_ab => ngwf_indices(B_NGWF)%rep%cross_overlap_tr%p
             overlap_ab => NULL()
          else
             ! reps differ and are both not in the set of reps for which
             ! cross_overlap is evaluated
             overlap_ab => NULL()
          end if
       else
          ! JCW: Expect NGWF overlap matrix to be standard overlap when
          !      NGWF bases for alpha and beta are identical
          overlap_ab => ngwf_indices(A_NGWF)%rep%overlap%p
       end if

       if (.not.basis_selector(3) == basis_selector(4)) then
          ! JCW: Expect NGWF overlap matrix to be cross_overlap when
          !      NGWF bases for gamma and delta differ and one of the bases
          !      is valence

          ! JCW: cross_overlap should only be available for cond or joint
          !      NGWF representations, so find which of the reps for gamma
          !      and delta NGWF indices are cond or joint
          ! JCW: cross_overlap evaluated in hamiltonian_dens_indep_matrices is
          !      <val|cond> and denskern_ab is expected to be P^{cond,val}
          if (all(x_overlap_rep(C_NGWF:D_NGWF).eqv.[.true.,.true.])) then
             ! Both are non-val, overlap (probably) not available
             overlap_cd => NULL()
          else if (all(x_overlap_rep(C_NGWF:D_NGWF).eqv.[.true.,.false.])) then
             ! C is non-val, expect <val|cond> cross overlap and P^{cond,val}
             overlap_cd => ngwf_indices(C_NGWF)%rep%cross_overlap%p
          else if (all(x_overlap_rep(C_NGWF:D_NGWF).eqv.[.false.,.true.])) then
             ! D is non-val, expect <cond|val> cross overlap and P^{val,cond}
             ! but cross_overlap_tr does not seem to be evaluated by default
             !overlap_cd => ngwf_indices(D_NGWF)%rep%cross_overlap_tr%p
             overlap_cd => NULL()
          else
             ! reps differ and are both not in the set of reps for which
             ! cross_overlap is evaluated
             overlap_cd => NULL()
          end if
       else
          ! JCW: Expect NGWF overlap matrix to be standard overlap when
          !      NGWF bases for gamma and delta are identical
          overlap_cd => ngwf_indices(C_NGWF)%rep%overlap%p
       end if
    else
       src_swex_h = REP_SWEX_HFX
       ngwf_indices(A_NGWF)%basis => ngwf_basis
       ngwf_indices(B_NGWF)%basis => ngwf_basis
       ngwf_indices(C_NGWF)%basis => ngwf_basis
       ngwf_indices(D_NGWF)%basis => ngwf_basis
       ngwf_indices(A_NGWF)%rep => rep
       ngwf_indices(B_NGWF)%rep => rep
       ngwf_indices(C_NGWF)%rep => rep
       ngwf_indices(D_NGWF)%rep => rep
       overlap_ab => ngwf_indices(A_NGWF)%rep%overlap%p
       overlap_cd => ngwf_indices(C_NGWF)%rep%overlap%p
    end if

    if(pub_hfx_output_detail >= VERBOSE) then
       do is = 1, pub_num_spins
          if (associated(overlap_ab)) then
             trKS_ab = sparse_trace(denskern_ab(is), overlap_ab)
             if(pub_on_root) then
                write(stdout,'(a,f12.3,a,i1,a)') 'HFx: Density kernel AB &
                     &normalised to ', trKS_ab, ' electrons (spin ',is,').'
             end if
          else
             if(pub_on_root) then
                write(stdout,'(a)') 'HFx: Density kernel AB &
                     &normalisation not reported as overlap AB not available.'
             end if
          end if
       end do
       do is = 1, pub_num_spins
          if (associated(overlap_cd)) then
             trKS_cd = sparse_trace(denskern_cd(is), overlap_cd)
             if(pub_on_root) then
                write(stdout,'(a,f12.3,a,i1,a)') 'HFx: Density kernel CD &
                     &normalised to ', trKS_cd, ' electrons (spin ',is,').'
             end if
          else
             if(pub_on_root) then
                write(stdout,'(a)') 'HFx: Density kernel CD &
                     &normalisation not reported as overlap CD not available.'
             end if
          end if
       end do
       if(pub_on_root) then
          if(present(energy_prefactor)) then
             write(stdout,'(a,f8.4)') 'HFx: Energy prefactor: ', energy_prefactor
          end if
          if(present(grad_prefactor)) then
             write(stdout,'(a,f8.4)') 'HFx: Grad prefactor: ', grad_prefactor
          end if
       end if
    end if

    ! jd: Initialise hash tables that depend on DKN. Assume DKN always changes
    !     across calls to hf_exchange_calculate().
    call hfx_dkn_blocks_ht_size(max_cached_dkn_blocks_cd, &
         max_cached_dkn_blocks_ab, ngwf_indices, natoms)

    call hash_table_init(packed_dkn_cd_blocks_ht, 'DKNBLKS_CD', 2, &
         max_cached_dkn_blocks_cd, max_cached_dkn_blocks_cd, &
         hfxstate%hash_table_info_unit)
    if(calcgradient) then
       call hash_table_init(packed_dkn_ab_blocks_hts(1), 'DKNBLKS_AB', 2, &
            max_cached_dkn_blocks_ab, max_cached_dkn_blocks_ab, &
            hfxstate%hash_table_info_unit)
       call hash_table_init(packed_dkn_ab_blocks_hts(2), 'TCDKNBLKS_AB', 2, &
            max_cached_dkn_blocks_ab, max_cached_dkn_blocks_ab, &
            hfxstate%hash_table_info_unit)
       call hash_table_init(packed_dkn_ab_blocks_hts(3), 'TCDKNBLKS_BA', 2, &
            max_cached_dkn_blocks_ab, max_cached_dkn_blocks_ab, &
            hfxstate%hash_table_info_unit)
    end if

    ! ========================================================================

    ! jd: Communicate DKN elements K^{CcDd}. Store them in HTs.
    call hfx_comms_dkn(packed_dkn_cd_blocks_ht, hfxstate%my_c_atoms, &
         hfxstate%n_my_c_atoms, hfxstate%sxs_atoms_nl, denskern_cd, &
         ngwf_indices(C_NGWF)%basis%max_on_atom, &
         ngwf_indices(D_NGWF)%basis%max_on_atom)

    call hash_table_list(packed_dkn_cd_blocks_ht,0)

    ! jd: Communicate DKN elements K^{AaBb} and tcK^{Aa}_{Bb}. Store them in HTs.
    !     This is only needed in the gradient calc.
    if(calcgradient) then
       call hfx_comms_dkn(packed_dkn_ab_blocks_hts(1), &
            hfxstate%my_b_atoms, hfxstate%n_my_b_atoms, hfxstate%x_atoms_nl,&
            denskern_ab, &
            ngwf_indices(B_NGWF)%basis%max_on_atom, &
            ngwf_indices(A_NGWF)%basis%max_on_atom, opt_transpose = .true.)
       call hash_table_list(packed_dkn_ab_blocks_hts(1),0)
       call hfx_comms_dkn(packed_dkn_ab_blocks_hts(2), &
            hfxstate%my_b_atoms, hfxstate%n_my_b_atoms, hfxstate%x_atoms_nl,&
            tc_denskern_ab, &
            ngwf_indices(B_NGWF)%basis%max_on_atom, &
            ngwf_indices(A_NGWF)%basis%max_on_atom, opt_transpose = .true.)
       call hash_table_list(packed_dkn_ab_blocks_hts(2),0)
       call hfx_comms_dkn(packed_dkn_ab_blocks_hts(3), &
            hfxstate%my_b_atoms, hfxstate%n_my_b_atoms, hfxstate%x_atoms_nl,&
            tc_denskern_ab, &
            ngwf_indices(B_NGWF)%basis%max_on_atom, &
            ngwf_indices(A_NGWF)%basis%max_on_atom)
       call hash_table_list(packed_dkn_ab_blocks_hts(3),0)
    end if

    ! jd: Main part of X matrix calculation (actual accumulation).
    call hfx_main(xmatrix, &
         hfxstate, remote_ngwf_cache_hts, rep, par, elements, cell, fftbox, &
         packed_dkn_cd_blocks_ht, packed_dkn_ab_blocks_hts, &
         calcgradient, src_swex_h, ngwf_indices, cov_grad, contra_grad, &
         precond_func_recip, grad_prefactor)

    ! ========================================================================

    ! jd: Clean up DKN hash tables
    call hash_table_free(packed_dkn_cd_blocks_ht)
    if(calcgradient) then
       do i_ht = 1, 3
          call hash_table_free(packed_dkn_ab_blocks_hts(i_ht))
       end do
    end if

    ! jd: Print relevant matrices if verbosity set to maximum
    if(pub_hfx_output_detail > PROLIX) then
       do is = 1, pub_num_spins
          if(pub_on_root) then
             write(stdout,'(a,i0,a)') CRLF//"-XMAT:--(spin ", is, ")-----------"
          end if
          call sparse_show_matrix(xmatrix(is))
          if(pub_on_root) then
             write(stdout,'(a,i0,a)') CRLF//"-DKN_AB:--(spin ", is, ")---------"
          end if
          call sparse_show_matrix(denskern_ab(is))
          if(pub_on_root) then
             write(stdout,'(a,i0,a)') CRLF//"-DKN_CD:--(spin ", is, ")---------"
          end if
          call sparse_show_matrix(denskern_cd(is))
          if(calcgradient) then
             if(pub_on_root) then
                write(stdout,'(a,i0,a)') CRLF//"-TC_DKN_AB:--(spin ", is, ")------"
             end if
             call sparse_show_matrix(tc_denskern_ab(is))
          end if
       end do
    end if

    ! jd:  Symmetrise X
    ! JCW: ...by averaging with its transpose, unless disabled via optional
    !      argument
    if (loc_symmetrise_xmatrix) call hfx_symmetrise_xmatrix(xmatrix)

    ! jd: Apply energy_prefactor
    if(present(energy_prefactor)) then
       if(energy_prefactor /= 1.0_DP) then
         do is = 1, pub_num_spins
            call sparse_scale(xmatrix(is), energy_prefactor)
         enddo
       end if
    endif

    ! cks: apply hyper Hartree Fock factor
    if (pub_hhf_nstates > 0) then
       do is = 1, pub_num_spins
          call sparse_scale(xmatrix(is), pub_hhf_factor)
       enddo
    endif

    ! jd: Print symmetrised (and possibly hHF-scaled) X
    if(pub_hfx_output_detail > PROLIX) then
       do is = 1, pub_num_spins
          if(pub_on_root) then
             write(stdout,'(a,i0,a)') CRLF//"-XMAT_SYMM:--(spin ", is, ")-------"
          end if
          call sparse_show_matrix(xmatrix(is))
          if(pub_on_root) then
             write(stdout,'(a)') CRLF
          end if
       end do
    end if

    ! JCW: Expect (based on call to sparse_put_block in hfx_comms_xmatrix) that
    !      xmatrix has NGWF index beta (B) as row index and NGWF index alpha
    !      (A) as column index. This should therefore be contracted with a
    !      denskern_ab that has NGWF index alpha (A) as row index and NGWF
    !      index (B) as column index.

    if(pub_hfx_output_detail >= VERBOSE .and. rep%postfix /= 'j') then
       hfx_energy = 0.0_DP              ! jd: ^ For joint calc dkn is zero
       do is = 1, pub_num_spins         !       so no point in calculating energy
          hfx_energy_for_spin = - 0.5_DP * sparse_trace(denskern_ab(is), &
               xmatrix(is))
          ! jd: Assumes that when 'rep2' is passed, we're in cond+val mode
          !     Verify if this is correct for joint+val and, in the future,
          !     for hybrid TDDFT.
          if(present(rep2)) then
             ! jd: Fixes reported energy for hybrid conduction to be consistent
             !     with the actual contribution from X.
             hfx_energy_for_spin = &
                  hfx_energy_for_spin * 4.0_DP / real(pub_num_spins, kind=DP)
          end if
          hfx_energy = hfx_energy + hfx_energy_for_spin
          if(pub_on_root) then
             ! jd: spin-degenerate case uses different DKN normalisations
             !     between val and cond.
             if(calcgradient .and. pub_num_spins == 1) then
                ! jd: Assumes that when 'rep2' is passed, we're in cond+val mode
                !     Verify if this is correct for joint+val and, in the future,
                !     for hybrid TDDFT.
                if(.not. present(rep2)) then
                   write(stdout,'(a,i0,a,f19.12,a)') &
                        'HFx: Hartree-Fock exchange energy (spin ',is,'): ', &
                        hfx_energy_for_spin*2.0_DP, ' Ha.'
                   write(stdout,'(a)') 'HFx: ^ (accounting for spin degeneracy &
                        &in density kernel in val gradient calc)'
                else
                   write(stdout,'(a,i0,a,f19.12,a)') &
                        'HFx: Hartree-Fock exchange energy (spin ',is,'): ', &
                        hfx_energy_for_spin*0.5_DP, ' Ha.'
                   write(stdout,'(a)') 'HFx: ^ (accounting for spin degeneracy &
                        &in density kernel in non-val gradient calc)'
                end if
             else ! non-gradient case, no issues
                write(stdout,'(a,i0,a,f19.12,a)') &
                     'HFx: Hartree-Fock exchange energy (spin ',is,'): ', &
                     hfx_energy_for_spin, ' Ha.'
             end if
          end if
       end do
       if(pub_on_root .and. pub_num_spins==2) then
         write(stdout,'(a,f19.12,a)') &
              'HFx: Hartree-Fock exchange energy (total):  ', &
              hfx_energy, ' Ha.'
       end if

    end if

    do is=1,pub_num_spins
       call sparse_copy(mat(is)%m(ireg,ireg),xmatrix(is))
       call sparse_destroy(xmatrix(is))
    end do
    deallocate(xmatrix,stat=ierr)
    call utils_dealloc_check(myself,'xmatrix',ierr)

    call timer_clock(myself//'_imbalance',1)
    call comms_barrier
    call timer_clock(myself//'_imbalance',2)
    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hf_exchange_calculate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! ---------------------------------------------------------------------------
  ! ---- P R I V A T E   S U B R O U T I N E S   A N D   F U N C T I O N S ----
  ! ---------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_main(xmatrix, &
       hfxstate, ngwfs_hts, rep, par, elements, cell, fftbox, &
       packed_dkn_cd_blocks_ht, packed_dkn_ab_blocks_hts, calcgradient, &
       src_swex_h, ngwf_indices, cov_grad, contra_grad, precond_func_recip, &
       grad_prefactor)
    !==========================================================================!
    ! Main stage of DKN-dependent part of HFx calculation.                     !
    ! - Expansions from all my B-C pairs are generated for all relevant PPDs.  !
    ! - Expansions are multiplied by K^CcDd, acted on \phi_Aa \phi_Dd and      !
    !   summed over PPDs.                                                      !
    ! - Contributions to X matrix are reduced over MPI ranks and the X matrix  !
    !   is assembled.                                                          !
    ! - NGWF gradient term 1 is calculated iff calcgradient == .true..         !
    ! - NGWF gradient term 2 is calculated iff calcgradient == .true..         !
    ! Intermediate results (contributions to X_AB from different procs) are    !
    ! stored in x_ab_contributions_ht. This is freed once X is assembled.      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   xmatrix (inout): Results are stored here.                              !
    !   hfxstate (inout): Contains most book-keeping.                          !
    !   ngwfs_hts (in): NGWF caches, only needed for extracting A PPD lists.   !
    !   rep (in): Expansion coefficients are extracted from here. Coefficients !
    !             always live in 'rep', regardless of whether or not we are    !
    !             doing mixed bases. Also used for swex quality.               !
    !   par (in): Usual parallel info. Needed for determining atom owners.     !
    !   elements (in): Used for centres.                                       !
    !   cell (in): The usual.                                                  !
    !   fftbox (in): Needed for weight and for NGWF gradient preconditioning.  !
    !   packed_dkn_cd_blocks_ht (in): Pre-communicated DKN blocks.             !
    !   packed_dkn_ab_blocks_hts (in): Pre-communicated DKN and tcK blocks.    !
    !   calcgradient (in): If .true., NGWF gradient will be calculated.        !
    !   src_swex_h (in): Handle to SW_EX associated with HFx.                  !
    !   ngwf_indices (in): Array of HFX_NGWF_INDEX instances containing        !
    !                      pointers to the NGWF basis and the NGWF_REP used    !
    !                      by atoms A, B, C and D.                             !
    !                                                                          !
    ! Optional arguments only needed if calcgradient == .true.:                !
    !   cov_grad (in/out): Covariant NGWF gradient, HFx contribution will      !
    !                      be *added* here.                                    !
    !   contra_grad (in/out): Contravariant NGWF gradient, HFx contribution    !
    !                      will be *added* here.                               !
    !   precond_func_recip (in): Needed for gradient preconditioning.          !
    !   grad_prefactor (in):      Defaults to 1.0. In conduction 0.5 can be    !
    !                             passed to take into account the fact that    !
    !                             only one set NGWFs is optimised, and the     !
    !                             other one is fixed (I think).                !
    !                             This also takes care of spin degeneracy      !
    !                             conventions in conduction (see the call to   !
    !                             hf_exchange_calculate() from conduction).    !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2019.                              !
    !==========================================================================!

    use comms, only: pub_total_num_procs, pub_on_root, pub_my_proc_id, &
         comms_barrier, comms_allgather
    use constants, only: SWEX_BB_CC_SYMMETRIC, SW_V, SW_O, METRIC_OVERLAP, &
         METRIC_ELECTROSTATIC, VERBOSE, PROLIX, stdout
    use datatypes, only: FUNCTIONS
    use fft_box, only: FFTBOX_INFO
    use hash_table, only: HT_HASH_TABLE, hash_table_lookup_nocount, &
         hash_table_list, hash_table_init, hash_table_free, hash_table_purge
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_hfx_metric, pub_num_spins, pub_hfx_output_detail, &
         pub_precond_recip
    use simulation_cell, only: CELL_INFO, minimum_image_number_to_displacement
    use sparse, only: SPAM3, sparse_scale
    use sw_resolution_of_identity, only: ATOM_CENTRE, swri_expansion_centres, &
         swri_init_centre, swri_library
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_assert, &
         utils_dealloc_check, utils_postfix_to_ngwf_set_name, utils_abort
    use visual, only: visual_ngwfs

    implicit none

    ! jd: Arguments
    type(SPAM3), intent(inout)               :: xmatrix(:) ! :spins
    type(HFX_STATE), intent(inout), target   :: hfxstate
    type(HT_HASH_TABLE), intent(in), target  :: ngwfs_hts(:) ! :handles
    type(NGWF_REP), intent(in)               :: rep
    type(PARAL_INFO), intent(in), pointer    :: par
    type(ELEMENT), intent(in)                :: elements(:)
    type(CELL_INFO), intent(in)              :: cell
    type(FFTBOX_INFO), intent(in)            :: fftbox
    type(HT_HASH_TABLE), intent(in), target  :: packed_dkn_cd_blocks_ht
    type(HT_HASH_TABLE), intent(in), target  :: packed_dkn_ab_blocks_hts(3)
    logical, intent(in)                      :: calcgradient
    integer, intent(in)                      :: src_swex_h
    type(HFX_NGWF_INDEX), intent(in)         :: ngwf_indices(A_NGWF:D_NGWF)
    type(FUNCTIONS), intent(inout), optional :: cov_grad
    type(FUNCTIONS), intent(inout), optional :: contra_grad
    real(kind=DP), intent(in), optional      :: precond_func_recip(:,:,:)
    real(kind=DP), intent(in), optional      :: grad_prefactor

    ! jd: Local variables
    type(HT_HASH_TABLE), target :: x_ab_contributions_ht
    type(HT_HASH_TABLE), target :: grad_akets_dkn_my_ht
    type(HT_HASH_TABLE), target :: grad_akets_dkn_loc_ht
    real(kind=DP), allocatable  :: f(:,:,:) ! fused PPD-Dd-pt, is, ngwf_b
    real(kind=DP), allocatable  :: grad_ket_in_ppds(:,:,:,:)
    integer, allocatable        :: a_ppd_list(:)
    integer :: maxslots_grad_kets_dkn
    integer :: n_a_ppd_list
    integer :: pair_idx
    integer :: b_idx
    integer :: prev_global_b
    integer :: orig_global_b, orig_global_c
    integer :: global_b, global_c
    integer :: ngwf_b, ngwf_c
    integer :: global_bb_ngwf_idx, global_cc_ngwf_idx
    integer :: first_ngwf_idx_of_b, first_ngwf_idx_of_c
    integer :: num_sws_in_expansion
    integer :: atoms(2), expansion_atoms(2)
    type(ATOM_CENTRE) :: centres(2), expansion_centres(2)
    real(kind=DP) :: coeffs(2*rep%swexes(src_swex_h)%quality%max_sws_per_centre)
    integer :: num_of_cached_coeffs
    integer :: coeffs_kind
    integer :: is
    logical :: table_full
    real(kind=DP) :: n_expansions_hits, n_expansions_misses
    real(kind=DP) :: n_prods_hits, n_prods_misses
    real(kind=DP) :: n_hits_all(0:pub_total_num_procs-1)
    real(kind=DP) :: n_misses_all(0:pub_total_num_procs-1)
    real(kind=DP) :: n_hits, n_misses, n_tot
    real(kind=DP) :: frac_hit
    integer :: proc
    integer :: mode
    integer :: img_at, img_ppd, img_pt
    integer :: ierr
    character(len=16) :: string1, string2, string3
    character(len=*), parameter :: myself = 'hfx_main'
    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    ! jd: Ensure consistency between optional parameters for gradient
    if(calcgradient) then
       call utils_assert(present(contra_grad) .and. &
            (present(precond_func_recip) .or. .not. pub_precond_recip), &
            myself//': When calcgradient is .true., certain optional arguments&
            & must be passed.')
    end if

    if(pub_hfx_metric == METRIC_ELECTROSTATIC) coeffs_kind = SW_V
    if(pub_hfx_metric == METRIC_OVERLAP) coeffs_kind = SW_O

    n_expansions_hits = 0.0_DP
    n_expansions_misses = 0.0_DP
    n_prods_hits = 0.0_DP
    n_prods_misses = 0.0_DP

    call hash_table_init(x_ab_contributions_ht, 'XAB', 2, &
         maxslots_x_ab_contributions, max_x_ab_contributions, &
         hfxstate%hash_table_info_unit)

    if(calcgradient) then
       maxslots_grad_kets_dkn = maxslots_grad_kets_dkn_loc_per_atom * par%nat
       call hash_table_init(grad_akets_dkn_loc_ht, 'AKETS_DKN_LOC', 4, &
            maxslots_grad_kets_dkn, maxslots_grad_kets_dkn, &
            hfxstate%hash_table_info_unit)
       maxslots_grad_kets_dkn = maxslots_grad_kets_dkn_my_per_atom * par%nat
       call hash_table_init(grad_akets_dkn_my_ht, 'AKETS_DKN_MY', 4, &
            maxslots_grad_kets_dkn, maxslots_grad_kets_dkn, &
            hfxstate%hash_table_info_unit)
    end if

    do is = 1, pub_num_spins
       call sparse_scale(xmatrix(is),0.0_DP)
    end do

    ! *************************************************************************
    ! *************************************************************************
    ! *** EXPANSIONS -> HT
    ! *************************************************************************
    ! *************************************************************************

    ! jd: Only expand if we don't have expansions from the previous LNV
    !     iteration. They may or may not be all of them (cachelimit constraint).
    !     Also, do not bother if the cachelimit is zero (then all are calculated
    !     on the fly later).
    if(hfxstate%expansions_ht%n_slots + hfxstate%expansions_ht%n_chained == 0 &
         .and. hfxstate%expansions_ht_in_use) then

       table_full = .false.

       ! -----------------------------------------------------------------------
       ! jd: Two outer 'mode' iterations.                               MODE 1,2
       ! -----------------------------------------------------------------------
       !     In mode == 1 we only do expansions that will prove the most costly
       !     to recalculate, ie. where both SWOPs are not in the cache. IOW,
       !     we skip the ones where either or both SWOPs are in the cache.
       !     In mode == 2 we do all remaining expansions, until the cache is
       !     full.
       !     In this way if expansion memory is limited (which is the case in
       !     practice, we focus on the expansions that will be the most costly
       !     to recalculate).
       loop_mode:                                                              &
       do mode = 1, 2

          ! jd: Skip mode 2 if table already full after mode 1.
          if(table_full) cycle

          ! --------------------------------------------------------------------
          ! jd: Loop over B-C's given to this rank                        MY B-C
          ! --------------------------------------------------------------------
          prev_global_b = garbage_int
          loop_BC:                                                             &
          do pair_idx = 1, hfxstate%n_my_bc_pairs
             global_b = hfxstate%my_bc_pairs(1,pair_idx)
             global_c = hfxstate%my_bc_pairs(2,pair_idx)
             orig_global_b = par%orig_atom(global_b)
             orig_global_c = par%orig_atom(global_c)
             first_ngwf_idx_of_b = ngwf_indices(B_NGWF)%basis%first_on_atom(global_b)
             first_ngwf_idx_of_c = ngwf_indices(C_NGWF)%basis%first_on_atom(global_c)

             ! jd: Ignore atoms not represented in this SWRI
             if(.not. elements(orig_global_b)%in_swri(hfxstate%hfx_swri_h)) cycle
             if(.not. elements(orig_global_c)%in_swri(hfxstate%hfx_swri_h)) cycle

             if(prev_global_b /= global_b .and. pub_hfx_output_detail >= VERBOSE) then
                write(stdout,'(a,i7,a,i0,a,i0,a)') '  - B: ', global_b, &
                     ' [', pub_my_proc_id, '] (', mode, ')'
             end if

             ! *****************************************************************
             ! jd: Find the PPDs of all A's that are X-neighbours with B. These
             !     are the PPDs where we are going to need expansions from BbCc.
             !     If this B is the same as in the last pair, no need to update.
             ! *****************************************************************
             if(global_b /= prev_global_b) then
                ! jd: Deallocate previous PPD list, if any
                if(allocated(a_ppd_list)) then
                   deallocate(a_ppd_list,stat=ierr)
                   call utils_dealloc_check(myself,'a_ppd_list',ierr, &
                        allocated_in='hfx_determine_a_ppd_list_for_b')
                end if
                call hfx_determine_a_ppd_list_for_b(a_ppd_list, n_a_ppd_list, &
                     ngwfs_hts, hfxstate, ngwf_indices(A_NGWF), global_b)
                if(hfxstate%expansions_ht_in_use) then
                   call hash_table_list(hfxstate%expansions_ht,0)
                end if
             end if

             ! jd: Figure out the expansion centres, sort them,
             !     pad with -1 if fewer than 2 are needed
             call swri_init_centre(centres(1), par, elements, global_b)
             call swri_init_centre(centres(2), par, elements, global_c)
             atoms = (/global_b, global_c/)
             call swri_expansion_centres(swri_library(hfxstate%hfx_swri_h), &
                  expansion_centres, expansion_atoms, num_sws_in_expansion, & ! out
                  centres, atoms, rep%swexes(src_swex_h)%quality, &           ! in
                  SWEX_BB_CC_SYMMETRIC(src_swex_h))                           ! in

             ! -----------------------------------------------------------------
             ! jd: All NGWFs b of B                                          bbb
             ! -----------------------------------------------------------------
             loop_ngwf_b:                                                      &
             do ngwf_b = 1, ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)
                global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

                ! --------------------------------------------------------------
                ! jd: All NGWFs c of C                                       ccc
                ! --------------------------------------------------------------
                loop_ngwf_c:                                                   &
                do ngwf_c = 1, ngwf_indices(C_NGWF)%basis%num_on_atom(global_c)
                   global_cc_ngwf_idx = first_ngwf_idx_of_c + ngwf_c - 1

                   ! ***********************************************************
                   ! jd: Look up expansion coefficients for this Bb Cc
                   ! ***********************************************************
                   if(.not. SWEX_BB_CC_SYMMETRIC(src_swex_h)) then
                      ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                      call hash_table_lookup_nocount(coeffs, num_of_cached_coeffs, &
                           rep%swexes(src_swex_h)%coeffs_ht, &
                           global_bb_ngwf_idx, global_cc_ngwf_idx, coeffs_kind, &
                           src_swex_h, 1) ! 1 is for formal spin-dependence
                   else
                      ! jd: Here there's Bb-Cc symmetry and we only store one triangle
                      call hash_table_lookup_nocount(coeffs, num_of_cached_coeffs, &
                           rep%swexes(src_swex_h)%coeffs_ht, &
                           min(global_bb_ngwf_idx,global_cc_ngwf_idx), &
                           max(global_bb_ngwf_idx,global_cc_ngwf_idx), coeffs_kind, &
                           src_swex_h, 1) ! 1 is for formal spin-dependence
                   end if

                   if(num_of_cached_coeffs == -1) then
                      call utils_abort(myself//': Internal error: coeffs hash &
                           &table too small or you forgot to call hf_exchange_&
                           &dkn_indep_stage() on '//&
                           trim(utils_postfix_to_ngwf_set_name(rep%postfix))//&
                           '%['//trim(rep%swexes(src_swex_h)%swex_name)//']: &
                           &Bb, Cc, coeffs_kind follow: ', &
                           global_bb_ngwf_idx, global_cc_ngwf_idx, coeffs_kind)
                   end if

                   ! ***********************************************************
                   ! jd: Calculate expanded potential from BbCc in all relevant
                   !     PPDs. Store results in expansion cache.
                   !     Uses OMP internally.
                   ! ***********************************************************
                   if(.not. table_full) then
                      call hfx_calc_expansions_bb_cc(hfxstate%expansions_ht, & ! in/out
                           table_full, &                                       ! out
                           SWEX_BB_CC_SYMMETRIC(src_swex_h), &                 ! in
                           swri_library(hfxstate%hfx_swri_h), &                ! in
                           rep%swexes(src_swex_h)%quality, &                   ! in
                           cell, a_ppd_list, n_a_ppd_list, &                   ! in
                           expansion_atoms, expansion_centres, &               ! in
                           coeffs, global_bb_ngwf_idx, global_cc_ngwf_idx, &   ! in
                           mode)                                               ! in
                   end if

                end do loop_ngwf_c

             end do loop_ngwf_b

             prev_global_b = global_b

          end do loop_BC

       end do loop_mode

       hfxstate%all_expansions_cached = .not. table_full

    end if

    if(.not. hfxstate%expansions_ht_in_use) then
       hfxstate%all_expansions_cached = .false.
    end if

    ! *************************************************************************
    ! *** END EXPANSIONS -> HT
    ! *************************************************************************

    ! jd: Deallocate final PPD list, if any (corner case of BC list empty)
    if(allocated(a_ppd_list)) then
       deallocate(a_ppd_list,stat=ierr)
       call utils_dealloc_check(myself,'a_ppd_list',ierr, &
            allocated_in='hfx_determine_a_ppd_list_for_b')
    end if

    ! --------------------------------------------------------------------------
    ! jd: Loop over B's given to this rank                                MY BBB
    ! --------------------------------------------------------------------------
    loop_B:                                                                    &
    do b_idx = 1, hfxstate%n_my_b_atoms
       global_b = hfxstate%my_b_atoms(b_idx)

       if(pub_hfx_output_detail >= VERBOSE) then
          write(stdout,'(a,i7,a,i0,a)') '  + B: ', global_b, ' [', &
               pub_my_proc_id, ']'
       end if

       ! jd: Deallocate previous PPD list, if any
       if(allocated(a_ppd_list)) then
          deallocate(a_ppd_list,stat=ierr)
          call utils_dealloc_check(myself,'a_ppd_list',ierr, &
               allocated_in='hfx_determine_a_ppd_list_for_b')
       end if

       ! jd: Figure out all PPDs on all A's relevant to current B
       call hfx_determine_a_ppd_list_for_b(a_ppd_list, n_a_ppd_list, &
            ngwfs_hts, hfxstate, ngwf_indices(A_NGWF), global_b)

       ! jd: f_bDdsi;B = \sum_Cc ( e_bCci;B * K^CcDds )
       !     This eliminates index Cc.
       call hfx_contract_over_dkn_cc(hfxstate, f, a_ppd_list, n_a_ppd_list, &
            packed_dkn_cd_blocks_ht, SWEX_BB_CC_SYMMETRIC(src_swex_h), &
            global_b, ngwf_indices, cell, rep, par, src_swex_h, elements, &
            n_expansions_hits, n_expansions_misses)

       ! jd: Calculate actual contributions to X from current B.
       !     This eliminates index Dd.
       call hfx_all_ad_for_this_b(x_ab_contributions_ht, &
            hfxstate, f, a_ppd_list, n_a_ppd_list, global_b, cell, &
            size(elements), fftbox%weight, ngwf_indices, n_prods_hits, &
            n_prods_misses, calcgradient, grad_ket_in_ppds)

       if(calcgradient) then
          call hfx_grad_ket_dkn_baccum(grad_akets_dkn_my_ht, &
               grad_ket_in_ppds, packed_dkn_ab_blocks_hts, &
               a_ppd_list, n_a_ppd_list, global_b, cell%n_pts, ngwf_indices, &
               hfxstate, grad_prefactor)
       end if

       deallocate(f, stat=ierr)
       call utils_dealloc_check(myself,'f',ierr, &
            allocated_in='hfx_contract_over_dkn_cc')

    end do loop_B

    ! jd: Deallocate final PPD list, if any (corner case of B list empty)
    if(allocated(a_ppd_list)) then
       deallocate(a_ppd_list,stat=ierr)
       call utils_dealloc_check(myself,'a_ppd_list',ierr, &
            allocated_in='hfx_determine_a_ppd_list_for_b')
    end if

    call timer_clock(myself//'_work_imbalance',1)
    call comms_barrier
    call timer_clock(myself//'_work_imbalance',2)

    ! jd: Communicate the contributions to X matrix to where they sparse-belong
    call hfx_comms_xmatrix(xmatrix, x_ab_contributions_ht, ngwf_indices)

    call timer_clock(myself//'_comms_imbalance',1)
    call comms_barrier
    call timer_clock(myself//'_comms_imbalance',2)

    if(calcgradient) then

       call timer_clock(myself//'_ngwf_grad1',1)

       ! -----------------------------------------------------------------------
       ! jd: NGWF gradient term 1
       ! -----------------------------------------------------------------------
       ! jd: It has been calculated already -- first stage took place in
       !     hfx_calc_all_ad_for_this_b(), then it got summed over my B's in
       !     hfx_grad_ket_dkn_baccum(). What is left is to communicate the
       !     contributions to where they sparse-belong (my -> loc). This also
       !     collates my-B sums from different procs into one sum over B.
       call hash_table_list(grad_akets_dkn_my_ht,0)
       call hfx_comms_gradient(grad_akets_dkn_loc_ht, &
            grad_akets_dkn_my_ht, cell%n_pts, ngwf_indices(A_NGWF))
       call hash_table_list(grad_akets_dkn_loc_ht,0)
       call hash_table_purge(grad_akets_dkn_my_ht)
       ! jd: Accumulate gradient PPDs proper for sparse-local A's.
       call hfx_gradient(cov_grad, contra_grad, grad_akets_dkn_loc_ht, &
            precond_func_recip, cell, fftbox, par, ngwf_indices(A_NGWF))
       call hash_table_purge(grad_akets_dkn_loc_ht)
       call timer_clock(myself//'_ngwf_grad1',2)

       ! -----------------------------------------------------------------------
       ! jd: NGWF gradient term 2
       ! -----------------------------------------------------------------------
       call timer_clock(myself//'_ngwf_grad2',1)
       call hfx_gradient_term2(hfxstate, grad_akets_dkn_my_ht, &
            packed_dkn_cd_blocks_ht, packed_dkn_ab_blocks_hts, ngwfs_hts, &
            rep%swexes(src_swex_h)%quality, src_swex_h, &
            swri_library(hfxstate%hfx_swri_h), cell, par, elements, &
            ngwf_indices, grad_prefactor)
       call hash_table_list(grad_akets_dkn_my_ht,0)
       ! jd: Communicate the contributions to gradient term 2 to where they
       !     sparse-belong (my -> loc), collate my-B sums from different procs
       !     into one B sum.
       call hfx_comms_gradient(grad_akets_dkn_loc_ht, &
            grad_akets_dkn_my_ht, cell%n_pts, ngwf_indices(A_NGWF))
       call hash_table_list(grad_akets_dkn_loc_ht,0)
       call hash_table_free(grad_akets_dkn_my_ht)

       if(.false.) then
          ! jd: Dump of NGWF gradient terms on a grid. Requires changing the
          !     'if' above and adding a %block species_ngwf_plot to the input
          !     file, plus dx_format T.
          call visual_ngwfs(cov_grad, ngwf_indices(A_NGWF)%basis, &
               'ngwf_cov_grad_pre_term2', elements, cell, fftbox, par)
          call visual_ngwfs(cov_grad, ngwf_indices(A_NGWF)%basis, &
               'ngwf_contra_grad_pre_term2', elements, cell, fftbox, par)
       end if

       ! jd: Accumulate gradient PPDs proper for sparse-local A's.
       call hfx_gradient(cov_grad, contra_grad, grad_akets_dkn_loc_ht, &
            precond_func_recip, cell, fftbox, par, ngwf_indices(A_NGWF))

       if(.false.) then
          ! jd: Dump of NGWF gradient terms on a grid. Requires changing the
          !     'if' above and adding a %block species_ngwf_plot to the input
          !     file, plus dx_format T.
          call visual_ngwfs(cov_grad, ngwf_indices(A_NGWF)%basis, &
               'ngwf_cov_grad_post_term2', elements, cell, fftbox, par)
          call visual_ngwfs(cov_grad, ngwf_indices(A_NGWF)%basis, &
               'ngwf_contra_grad_post_term2', elements, cell, fftbox, par)
       end if

       call hash_table_free(grad_akets_dkn_loc_ht)
       call timer_clock(myself//'_ngwf_grad2',2)

    end if

    call hash_table_list(x_ab_contributions_ht,0)
    call hash_table_free(x_ab_contributions_ht)

    ! -------------------------------------------------------------------------
    ! jd: User feedback regarding caching of expansions and products
    ! -------------------------------------------------------------------------
    if(pub_hfx_output_detail >= PROLIX .or. (pub_hfx_output_detail >= VERBOSE &
         .and. .not. hfxstate%expansions_feedback_printed)) then

       if(hfxstate%expansions_ht_in_use) then
          call comms_allgather(n_hits_all, n_expansions_hits, 1, &
               gather_not_allgather = .true.)
          call comms_allgather(n_misses_all, n_expansions_misses, 1, &
               gather_not_allgather = .true.)

          if(pub_on_root) then
             write(stdout,'(a)') 'HFx: +'//repeat('-',71)//'+'
             write(stdout,'(a)') 'HFx: |    MPI |                          Expansi&
                  &on cache                     |'
             write(stdout,'(a)') 'HFx: |   rank |           hits |         misses &
                  &|          total | hit ratio |'
             write(stdout,'(a)') 'HFx: +'//repeat('-',71)//'+'
             do proc = 0, pub_total_num_procs-1
                n_hits = n_hits_all(proc)
                n_misses = n_misses_all(proc)
                n_tot = n_hits + n_misses
                if(n_tot > 0.0_DP) then
                   frac_hit = n_hits / n_tot * 100.0_DP
                else
                   frac_hit = 1.0_DP
                end if
                write(string1,'(f16.0)') n_hits
                write(string2,'(f16.0)') n_misses
                write(string3,'(f16.0)') n_tot
                write(stdout,'(a,i7,a,a15,a,a15,a,a15,a,f8.2,a)') 'HFx: |', proc, &
                     ' |', string1(1:15), ' |', string2(1:15), &
                     ' |', string3(1:15), ' |', frac_hit ,' % |'
             end do
             write(stdout,'(a)') 'HFx: +'//repeat('-',71)//'+'
          end if
       else
          if(pub_on_root) then
             write(stdout,'(a)') 'HFx: +'//repeat('-',71)//'+'
             write(stdout,'(a)') 'HFx: | Expansion cache is not used &
                  &(cache_limit_for_expansions == 0).        |'
             write(stdout,'(a)') 'HFx: +'//repeat('-',71)//'+'
          end if
       end if

       hfxstate%expansions_feedback_printed = .true.

    end if

    if(pub_hfx_output_detail >= PROLIX .or. (pub_hfx_output_detail >= VERBOSE &
         .and. .not. hfxstate%prods_feedback_printed)) then

       if(hfxstate%expansions_ht_in_use) then
          call comms_allgather(n_hits_all, n_prods_hits, 1, &
               gather_not_allgather = .true.)
          call comms_allgather(n_misses_all, n_prods_misses, 1, &
               gather_not_allgather = .true.)

          if(pub_on_root) then
             write(stdout,'(a)') 'HFx: +'//repeat('-',71)//'+'
             write(stdout,'(a)') 'HFx: |    MPI |                    AD NGWF produ&
                  &ct cache                     |'
             write(stdout,'(a)') 'HFx: |   rank |           hits |         misses &
                  &|          total | hit ratio |'
             write(stdout,'(a)') 'HFx: +'//repeat('-',71)//'+'
             do proc = 0, pub_total_num_procs-1
                n_hits = n_hits_all(proc)
                n_misses = n_misses_all(proc)
                n_tot = n_hits + n_misses
                if(n_tot > 0.0_DP) then
                   frac_hit = n_hits / n_tot * 100.0_DP
                else
                   frac_hit = 1.0_DP
                end if
                write(string1,'(f16.0)') n_hits
                write(string2,'(f16.0)') n_misses
                write(string3,'(f16.0)') n_tot
                write(stdout,'(a,i7,a,a15,a,a15,a,a15,a,f8.2,a)') 'HFx: |', proc, &
                     ' |', string1(1:15), ' |', string2(1:15), &
                     ' |', string3(1:15), ' |', frac_hit ,' % |'
             end do
             write(stdout,'(a)') 'HFx: +'//repeat('-',71)//'+'
          else
             if(pub_on_root) then
                write(stdout,'(a)') 'HFx: +'//repeat('-',71)//'+'
                write(stdout,'(a)') 'HFx: | AD product cache is not used.      &
                     &                                 |'
                write(stdout,'(a)') 'HFx: +'//repeat('-',71)//'+'
             end if
          end if
       end if

       hfxstate%prods_feedback_printed = .true.

    end if

#if 0
    write(*,*) '@IMAGE MAP'
    do img_at = 1, par%nat
       do img_ppd = 1, cell%n_ppds
          do img_pt = 1, cell%n_pts
             if(swri_library(hfxstate%hfx_swri_h)%&
                  image_map(img_pt,img_ppd,img_at) /= garbage_int) then
                write(*,'(i5,i14,f14.8,f14.8,f14.8)') img_at, &
                     swri_library(hfxstate%hfx_swri_h)%&
                     image_map(img_pt,img_ppd,img_at), &
                     minimum_image_number_to_displacement(&
                     swri_library(hfxstate%hfx_swri_h)%&
                     image_map(img_pt,img_ppd,img_at),cell)
             end if

          end do
       end do
    end do

#endif

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hfx_main

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_determine_bc_pairs(hfxstate, & ! inout
       par, ngwf_indices)                       ! in
    !==========================================================================!
    ! Determines all (system-wide) B-C atom pairs (C's are X-neighbours of B). !
    ! Divides the list of pairs across MPI ranks. Results get stored in        !
    ! hfxstate's my_bc_pairs.                                                  !
    !                                                                          !
    ! Attempts to balance work by trying to assign similar load, measured by   !
    ! n_ngwfs_b * n_ngwfs_c per each pair. In any case, the number of pairs is !
    ! bound to be much larger than the number of MPI ranks, so load should be  !
    ! fairly well balanced.                                                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! hfxstate (inout): The HFX_STATE container that stores all HFx data that  !
    !                   need to persist across calls.                          !
    ! par (in): Needed to figure out who rank-owns which atoms.                !
    ! ngwf_indices (in): Array of HFX_NGWF_INDEX instances containing          !
    !                    pointers to the NGWF basis and the NGWF_REP used      !
    !                    by atoms A, B, C and D (only FUNC_BASIS instances     !
    !                    for atoms B and C used).                              !
    !--------------------------------------------------------------------------!
    ! Notes and caveats:                                                       !
    !  - Allocates my_bc_pairs to desired size on each proc. Caller is         !
    !    responsible for freeing my_bc_pairs later (happens in hfx_free_arrays)!
    !  - The global list of BC pairs is allocated on and communicated to all   !
    !    ranks (allgather). This is currently not needed, but might prove      !
    !    useful if one day we'd like to know the my-owners of sparse-local     !
    !    atoms. The 'owner' field is filled, but not used ATM.                 !
    !  @unsure about ignoring pairs outside of SWRI, possibly needed.          !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2019.                              !
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_root_proc_id, pub_on_root, &
         pub_total_num_procs, comms_reduce, comms_allgather, comms_irecv, &
         comms_free, comms_wait, comms_barrier
    use constants, only: VERBOSE, PROLIX
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_hfx_output_detail
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_assert, &
         utils_abort, utils_alloc_check, utils_dealloc_check, utils_int_to_str

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target :: hfxstate
    type(PARAL_INFO), intent(in)           :: par
    type(HFX_NGWF_INDEX), intent(in)       :: ngwf_indices(A_NGWF:D_NGWF)

    ! jd: Local variables
    ! idx1: 1:B, 2:C, 3:n_ngwfs_B*n_ngwf_C (~effort), 4: unused
    !       5: owner (rank) of this pair, so always pub_my_proc_id
    ! stores data only for sparse-local B's
    integer, allocatable              :: local_bc_pairs(:,:)
    ! idx1: 1:B, 2:C, 3:n_ngwfs_B*n_ngwf_C (~effort), 4: accum of #3
    !       5: owner (rank) of this pair
    ! stores data for all B's, but only on root
    integer, allocatable              :: global_bc_pairs(:,:)
    integer :: global_b, local_b
    integer :: global_c, c_idx, n_atoms_c
    integer :: n_pairs_local, n_pairs_global, pair_idx
    integer :: n_ngwfs_on_b, n_ngwfs_on_c
    integer :: accum_n_ngwf_prods
    integer :: n_ngwf_prods_per_rank
    integer :: cur_rank
    integer :: start_pair
    integer :: ierr
    integer :: list_length
    integer :: recv_handle
    character(len=*), parameter :: myself = 'hfx_determine_bc_pairs'

    ! -----------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    ! jd: Sanity check on neighbour list
    call utils_assert(hfxstate%s_atoms_nl%populated, &
         myself//': s_atoms neighbour list not populated.')

    ! --------------------------------------------------------------------------
    ! [1] Count the number of B-C pairs, where B's are sparse-local and C's are
    !     B's S-neighbours. This will be used to dimension the array of pairs.
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    ! jd: Loop over B's -- all atoms on this rank                            BBB
    ! --------------------------------------------------------------------------
    n_pairs_local = 0
    loop_B1:                                                                   &
    do local_b=1, par%num_atoms_on_proc(pub_my_proc_id)
       global_b = par%first_atom_on_proc(pub_my_proc_id) + local_b-1
       n_atoms_c = hfxstate%s_atoms_nl%n_neighbours(global_b)
       n_pairs_local = n_pairs_local + n_atoms_c
    end do loop_B1

    allocate(local_bc_pairs(5,n_pairs_local),stat=ierr)
    call utils_alloc_check(myself,'local_bc_pairs',ierr)

    ! --------------------------------------------------------------------------
    ! [2] Find and store the B-C pairs, where B's are sparse-local and C's are
    !     B's S-neighbours.
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    ! jd: Loop over B's -- all atoms on this rank                            BBB
    ! --------------------------------------------------------------------------
    pair_idx = 1
    loop_B2:                                                                   &
    do local_b=1, par%num_atoms_on_proc(pub_my_proc_id)
       global_b = par%first_atom_on_proc(pub_my_proc_id) + local_b-1
       n_ngwfs_on_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)

       ! -----------------------------------------------------------------------
       ! jd: Loop over C's that are S-neighbours of B                        CCC
       ! -----------------------------------------------------------------------
       loop_C2:                                                                &
       do c_idx = hfxstate%s_atoms_nl%first_idx(global_b), &
            hfxstate%s_atoms_nl%last_idx(global_b)
          global_c = hfxstate%s_atoms_nl%neighbours(c_idx)
          n_ngwfs_on_c = ngwf_indices(C_NGWF)%basis%num_on_atom(global_c)

          local_bc_pairs(1,pair_idx) = global_b
          local_bc_pairs(2,pair_idx) = global_c
          local_bc_pairs(3,pair_idx) = n_ngwfs_on_b * n_ngwfs_on_c
          local_bc_pairs(4,pair_idx) = -1
          local_bc_pairs(5,pair_idx) = -1

          pair_idx = pair_idx + 1

       end do loop_C2
    end do loop_B2

    call utils_assert(pair_idx-1 == n_pairs_local, myself//': Logic error [1]',&
         pair_idx-1, n_pairs_local)

    ! --------------------------------------------------------------------------
    ! [3] Combine results from all procs on root
    ! --------------------------------------------------------------------------
    n_pairs_global = n_pairs_local
    call comms_reduce('SUM', n_pairs_global)
    allocate(global_bc_pairs(5,n_pairs_global),stat=ierr)
    call utils_alloc_check(myself,'global_bc_pairs',ierr)

    call comms_allgather(global_bc_pairs, local_bc_pairs, 5*n_pairs_local, &
         (/-1/), (/-1/))

    deallocate(local_bc_pairs,stat=ierr)
    call utils_dealloc_check(myself,'local_bc_pairs',ierr)

    ! --------------------------------------------------------------------------
    ! [4] Divide across ranks, trying to give every rank similar effort
    !     (measured by total n_NGWF_b*n_ngwf_c).
    ! --------------------------------------------------------------------------
#ifdef MPI
    ! jd: Have all procs posts receives for the sizes of the actual data
    call comms_irecv(pub_root_proc_id, list_length, handle = recv_handle)

    ! jd: On root divvy up the list and initiate sends
    if(pub_on_root) then
       ! jd: Calculate cumulative n_NGWF_b*n_NGWF_c (index 4)
       accum_n_ngwf_prods = 0
       do pair_idx = 1, n_pairs_global
          accum_n_ngwf_prods = accum_n_ngwf_prods + global_bc_pairs(3,pair_idx)
          global_bc_pairs(4,pair_idx) = accum_n_ngwf_prods
       end do

       if(pub_hfx_output_detail > PROLIX) then
          call internal_show_pairs(global_bc_pairs, n_pairs_global, &
               'Global B-C pairs (before distribution)')
       end if

       ! jd: Divide the global number of NGWF products in B-C's over MPI ranks,
       !     rounding up
       n_ngwf_prods_per_rank = accum_n_ngwf_prods / pub_total_num_procs
       if(mod(accum_n_ngwf_prods,pub_total_num_procs) /= 0) then
          n_ngwf_prods_per_rank = n_ngwf_prods_per_rank + 1
       end if

       ! jd: Find out each rank's portion
       cur_rank = 0
       start_pair = 1
       do pair_idx = 1, n_pairs_global
          accum_n_ngwf_prods = global_bc_pairs(4,pair_idx)
          if(accum_n_ngwf_prods >= (cur_rank+1) * n_ngwf_prods_per_rank) then
             ! jd: This rank's portion is [start_pair,pair_idx]
             global_bc_pairs(5,start_pair:pair_idx) = cur_rank
             call internal_send_pairs(global_bc_pairs, start_pair, pair_idx, &
                  cur_rank)
             ! jd: Next rank's portion starts after the end of this one's
             start_pair = pair_idx + 1
             cur_rank = cur_rank + 1
          end if
       end do
       ! jd: If n_pairs_global does not divide evenly between ranks, then
       !     last rank gets all that's left
       if(cur_rank < pub_total_num_procs) then
          pair_idx = n_pairs_global
          ! jd: This rank's portion is [start_pair,pair_idx]
          global_bc_pairs(5,start_pair:pair_idx) = cur_rank
          call internal_send_pairs(global_bc_pairs, start_pair, pair_idx, &
               cur_rank)
          cur_rank = cur_rank + 1
       end if
       call utils_assert(cur_rank == pub_total_num_procs, &
            myself//': Something went wrong dividing B-C pairs across ranks', &
            cur_rank, pub_total_num_procs, n_pairs_global)
    end if

    ! jd: Wait for the recv on the size of the actual data
    call comms_wait(recv_handle)

    ! jd: Now list_length is known and we're ready to alloc and receive
    hfxstate%n_my_bc_pairs = list_length / 5
#else
    ! jd: Non-MPI (commless) version
    hfxstate%n_my_bc_pairs = n_pairs_global
#endif
    allocate(hfxstate%my_bc_pairs(5,hfxstate%n_my_bc_pairs),stat=ierr)
    call utils_alloc_check(myself,'my_bc_pairs',ierr)

#ifdef MPI
    call comms_irecv(pub_root_proc_id, hfxstate%my_bc_pairs, handle = recv_handle)
    ! jd: Wait for the recv on the actual data
    call comms_wait(recv_handle)

    ! jd: Fix the 'accum' column to count in a proc-local fashion
    hfxstate%my_bc_pairs(4,1) = hfxstate%my_bc_pairs(3,1)
    do pair_idx = 2, hfxstate%n_my_bc_pairs
       hfxstate%my_bc_pairs(4,pair_idx) = hfxstate%my_bc_pairs(3,pair_idx) + &
            hfxstate%my_bc_pairs(4,pair_idx-1)
    end do

    ! jd: Wait for the sends on root to complete before the send buffer is nuked
    if(pub_on_root) then
       call comms_free
    end if
#else
    ! jd: Non-MPI (commless) version
    hfxstate%my_bc_pairs(:,:) = global_bc_pairs(:,:)
    hfxstate%my_bc_pairs(5,:) = 0 ! jd: Root owns every pair
#endif

    ! jd: User feedback -- detailed
    if(pub_hfx_output_detail >= PROLIX) then
       if(pub_on_root) then
          call internal_show_pairs(global_bc_pairs, n_pairs_global, &
               'Global B-C pairs (after distribution)')
       end if
       if(pub_hfx_output_detail > PROLIX) then
          call internal_show_pairs(hfxstate%my_bc_pairs, hfxstate%n_my_bc_pairs, &
               'B-C pairs of rank #'//trim(utils_int_to_str(pub_my_proc_id)))
       end if
    end if

    ! jd: User feedback -- shortened
    if(pub_hfx_output_detail >= VERBOSE .and. pub_on_root) then
       call internal_show_stats(global_bc_pairs, n_pairs_global, par%nat)
    end if

    deallocate(global_bc_pairs,stat=ierr)
    call utils_dealloc_check(myself,'global_bc_pairs',ierr)

    ! jd: Sanity check. I should be the owner of all my BC pairs
    do pair_idx = 1, hfxstate%n_my_bc_pairs
       if(hfxstate%my_bc_pairs(5,pair_idx) /= pub_my_proc_id) then
          call utils_abort(myself//': Unexpected owner of BC pair.', &
               hfxstate%my_bc_pairs(5,pair_idx), pub_my_proc_id, pair_idx)
       end if
    end do

    call comms_barrier
    call timer_clock(myself,2)
    call utils_trace_out(myself)

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_show_stats(pair_list, n_pairs, global_nat)
      !========================================================================!
      ! Outputs a digest of BC distribution to stdout in human-readable format.!
      !------------------------------------------------------------------------!
      ! Caveats:                                                               !
      !  - pair_list is global (over all procs).                               !
      !------------------------------------------------------------------------!
      ! written by Jacek Dziedzic in November 2019.                            !
      !========================================================================!

      use comms, only: pub_total_num_procs
      use constants, only: stdout, CRLF
      use utils, only: utils_banner, utils_alloc_check, utils_dealloc_check

      implicit none

      ! jd: Arguments
      integer, intent(in)          :: pair_list(:,:)
      integer, intent(in)          :: n_pairs
      integer, intent(in)          :: global_nat

      ! jd: Local variables
      integer :: i_pair
      integer :: i_proc
      integer :: n_b_atoms, n_c_atoms, n_b_c_atoms
      integer :: sum_bc, sum_ngwf_prod, sum_load
      integer :: min_n_b_atoms, min_n_c_atoms, min_n_b_c_atoms
      integer :: min_sum_bc, min_sum_ngwf_prod, min_sum_load
      integer :: max_n_b_atoms, max_n_c_atoms, max_n_b_c_atoms
      integer :: max_sum_bc, max_sum_ngwf_prod, max_sum_load
      real(kind=DP) :: avg_n_b_atoms, avg_n_c_atoms, avg_n_b_c_atoms
      real(kind=DP) :: avg_sum_bc, avg_sum_ngwf_prod, avg_sum_load
      logical, allocatable :: b_atoms(:)
      logical, allocatable :: c_atoms(:)
      character(len=*), parameter :: myself = &
           'hfx_determine_bc_pairs::internal_show_stats'

      ! -----------------------------------------------------------------------

      write(stdout,'(a)') utils_banner('=', 'Load distribution', &
           opt_terminal_width = 77)
      write(stdout,'(a)') '|     MPI |      # BC |     # B |     # C |    #B+C &
           &|    # NGWF | estimated |'//CRLF//'|    rank |     pairs |   atoms &
           &|   atoms |   atoms |  products |      load |'//CRLF//&
           '+-----------------------------------------------------------&
           &----------------+'

      ! jd: Have two bitmaps: all possible B atoms, C atoms.
      allocate(b_atoms(global_nat),stat=ierr)
      call utils_alloc_check(myself,'b_atoms',ierr)
      allocate(c_atoms(global_nat),stat=ierr)
      call utils_alloc_check(myself,'c_atoms',ierr)

      min_sum_bc = huge(1)
      min_sum_ngwf_prod = huge(1)
      min_sum_load = huge(1)
      min_n_b_atoms = huge(1)
      min_n_c_atoms = huge(1)
      min_n_b_c_atoms = huge(1)

      max_sum_bc = -1
      max_sum_ngwf_prod = -1
      max_sum_load = -1
      max_n_b_atoms = -1
      max_n_c_atoms = -1
      max_n_b_c_atoms = -1

      avg_sum_bc = 0.0_DP
      avg_sum_ngwf_prod = 0.0_DP
      avg_sum_load = 0.0_DP
      avg_n_b_atoms = 0.0_DP
      avg_n_c_atoms = 0.0_DP
      avg_n_b_c_atoms = 0.0_DP

      do i_proc = 0, pub_total_num_procs-1

         sum_bc = 0
         sum_ngwf_prod = 0
         sum_load = 0

         b_atoms(:) = .false.
         c_atoms(:) = .false.

         do i_pair = 1, n_pairs

            if(pair_list(5,i_pair) /= i_proc) cycle

            global_b = pair_list(1,i_pair)
            global_c = pair_list(2,i_pair)
            b_atoms(global_b) = .true.
            c_atoms(global_c) = .true.

            sum_bc = sum_bc + 1
            sum_ngwf_prod = sum_ngwf_prod + pair_list(3,i_pair)
            sum_load = sum_ngwf_prod

         end do

         n_b_atoms = count(b_atoms)
         n_c_atoms = count(c_atoms)
         n_b_c_atoms = n_b_atoms + n_c_atoms

         min_sum_bc = min(sum_bc,min_sum_bc)
         min_sum_ngwf_prod = min(sum_ngwf_prod,min_sum_ngwf_prod)
         min_sum_load = min(sum_load,min_sum_load)
         min_n_b_atoms = min(n_b_atoms,min_n_b_atoms)
         min_n_c_atoms = min(n_c_atoms,min_n_c_atoms)
         min_n_b_c_atoms = min(n_b_c_atoms,min_n_b_c_atoms)

         max_sum_bc = max(sum_bc,max_sum_bc)
         max_sum_ngwf_prod = max(sum_ngwf_prod,max_sum_ngwf_prod)
         max_sum_load = max(sum_load,max_sum_load)
         max_n_b_atoms = max(n_b_atoms,max_n_b_atoms)
         max_n_c_atoms = max(n_c_atoms,max_n_c_atoms)
         max_n_b_c_atoms = max(n_b_c_atoms,max_n_b_c_atoms)

         avg_sum_bc = avg_sum_bc + real(sum_bc,kind=DP) / real(pub_total_num_procs,kind=DP)
         avg_sum_ngwf_prod = avg_sum_ngwf_prod + real(sum_ngwf_prod,kind=DP) / real(pub_total_num_procs,kind=DP)
         avg_sum_load = avg_sum_load + real(sum_load,kind=DP) / real(pub_total_num_procs,kind=DP)
         avg_n_b_atoms = avg_n_b_atoms + real(n_b_atoms,kind=DP) / real(pub_total_num_procs,kind=DP)
         avg_n_c_atoms = avg_n_c_atoms + real(n_c_atoms,kind=DP) / real(pub_total_num_procs,kind=DP)
         avg_n_b_c_atoms = avg_n_b_c_atoms + real(n_b_c_atoms,kind=DP) / real(pub_total_num_procs,kind=DP)

         write(stdout,'(a,i8,a,i10,a,i8,a,i8,a,i8,a,i10,a,i10,a)') &
              '|', i_proc, ' |', sum_bc, ' |', n_b_atoms, ' |', n_c_atoms, &
              ' |', n_b_c_atoms, ' |', sum_ngwf_prod, &
              ' |', sum_load, ' |'

      end do
      write(stdout,'(a)') '+---------------------------------------------------&
           &------------------------+'

      write(stdout,'(a,i10,a,i8,a,i8,a,i8,a,i10,a,i10,a)') &
           '|     min |', min_sum_bc, ' |', min_n_b_atoms, ' |', min_n_c_atoms,&
           ' |', min_n_b_c_atoms, ' |', min_sum_ngwf_prod, &
           ' |', min_sum_load, ' |'
      write(stdout,'(a,i10,a,i8,a,i8,a,i8,a,i10,a,i10,a)') &
           '|     max |', max_sum_bc, ' |', max_n_b_atoms, ' |', max_n_c_atoms,&
           ' |', max_n_b_c_atoms, ' |', max_sum_ngwf_prod, &
           ' |', max_sum_load, ' |'

      write(stdout,'(a,f10.1,a,f8.1,a,f8.1,a,f8.1,a,f10.1,a,f10.1,a)') &
           '| average |', avg_sum_bc, ' |', avg_n_b_atoms, ' |', avg_n_c_atoms,&
           ' |', avg_n_b_c_atoms, ' |', avg_sum_ngwf_prod, &
           ' |', avg_sum_load, ' |'

      deallocate(b_atoms,stat=ierr)
      call utils_dealloc_check(myself,'b_atoms',ierr)
      deallocate(c_atoms,stat=ierr)
      call utils_dealloc_check(myself,'c_atoms',ierr)

      write(stdout,'(a)') repeat('=',77)

    end subroutine internal_show_stats

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_show_pairs(pair_list, n_pairs, list_name)
      !========================================================================!
      ! Outputs a list of pairs to stdout in human-readable format.            !
      ! When called in parallel for many MPI ranks simultaneously, the ranks   !
      ! will race, so this is best used in serial or for debugging only.       !
      !========================================================================!

      use comms, only: pub_my_proc_id
      use constants, only: stdout
      use utils, only: utils_banner

      implicit none

      ! jd: Arguments
      integer, intent(in)          :: pair_list(:,:)
      integer, intent(in)          :: n_pairs
      character(len=*), intent(in) :: list_name

      ! jd: Local variables
      integer :: i_pair

      ! -----------------------------------------------------------------------

      write(stdout,'(a)') utils_banner('=', trim(list_name), &
           opt_terminal_width = 66)
      write(stdout,'(a)') &
           '|  Pair  | B atom | C atom | NGWF cost |  accum  | owner | rank |'
      do i_pair = 1, n_pairs
         write(stdout,'(a,i8,a,i8,a,i8,a,i11,a,i9,a,i7,a,i6,a)') &
              '|', i_pair, '|', &
              pair_list(1,i_pair), '|', pair_list(2,i_pair), '|', &
              pair_list(3,i_pair), '|', pair_list(4,i_pair), '|', &
              pair_list(5,i_pair), '|', pub_my_proc_id, '|'
      end do
      write(stdout,'(a)') repeat('=',66)

    end subroutine internal_show_pairs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_send_pairs(pair_list, start_pair, end_pair, dest_rank)
      !========================================================================!
      ! Initiates a send of a range from a list of pairs to a remote MPI rank. !
      ! The caller is expected to wait for completion using comms_free().      !
      !========================================================================!

      use comms, only: comms_send

      implicit none

      ! jd: Arguments
      integer, intent(in) :: pair_list(:,:)
      integer, intent(in) :: start_pair, end_pair
      integer, intent(in) :: dest_rank

      ! jd: Local variables
      integer :: length

      ! -----------------------------------------------------------------------

      length = (end_pair - start_pair + 1) * 5
      call comms_send(dest_rank, length)
      call comms_send(dest_rank, pair_list(1,start_pair), length)

    end subroutine internal_send_pairs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine hfx_determine_bc_pairs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_split_pair_list(hfxstate, global_nat, ngwf_indices)
    !==========================================================================!
    ! Splits a list of BC atom pairs that are mine into lists of B's and C's.  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (in): List of BC pairs is obtained from here, results go here.!
    !   global_nat (in): size(elements).                                       !
    !   ngwf_indices (in): Array of HFX_NGWF_INDEX instances containing        !
    !                      pointers to the NGWF basis and the NGWF_REP used    !
    !                      by atoms A, B, C and D (only FUNC_BASIS instances   !
    !                      for atoms B are used -- to count Bb NGWFs).         !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2019.                              !
    !==========================================================================!

    use comms, only: comms_barrier
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target :: hfxstate
    integer                                :: global_nat
    type(HFX_NGWF_INDEX), intent(in)       :: ngwf_indices(A_NGWF:D_NGWF)

    ! jd: Local variables
    logical, allocatable :: b_atoms(:), c_atoms(:)
    integer :: pair_idx
    integer :: global_b, global_c
    integer :: my_b_atom_idx, my_c_atom_idx
    integer :: n_ngwfs_b
    integer :: ierr
    character(len=*), parameter :: myself = 'hfx_split_pair_list'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    ! --------------------------------------------------------------------------
    ! jd: [1] Have two bitmaps: all possible B atoms, C atoms.
    ! --------------------------------------------------------------------------
    allocate(b_atoms(global_nat),stat=ierr)
    call utils_alloc_check(myself,'b_atoms',ierr)
    b_atoms(:) = .false.
    allocate(c_atoms(global_nat),stat=ierr)
    call utils_alloc_check(myself,'c_atoms',ierr)
    c_atoms(:) = .false.

    ! --------------------------------------------------------------------------
    ! jd: [2] Go over all B's and C's and tag them in the bitmaps.
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    ! jd: Loop over B-C's given to this rank                              MY B-C
    ! --------------------------------------------------------------------------
    loop_BC:                                                                   &
    do pair_idx = 1, hfxstate%n_my_bc_pairs
       global_b = hfxstate%my_bc_pairs(1,pair_idx)
       global_c = hfxstate%my_bc_pairs(2,pair_idx)
       b_atoms(global_b) = .true.
       c_atoms(global_c) = .true.
    end do loop_BC

    ! --------------------------------------------------------------------------
    ! jd: [3] Now that we know which B, C atoms are in the bitmaps, and how many
    !         there are, we can construct the lists.
    ! --------------------------------------------------------------------------
    hfxstate%n_my_b_atoms = count(b_atoms)
    allocate(hfxstate%my_b_atoms(hfxstate%n_my_b_atoms),stat=ierr)
    call utils_alloc_check(myself,'my_b_atoms',ierr)
    allocate(hfxstate%my_bb_ngwf_offsets(hfxstate%n_my_b_atoms),stat=ierr)
    call utils_alloc_check(myself,'my_bb_ngwf_offsets',ierr)
    hfxstate%n_my_c_atoms = count(c_atoms)
    allocate(hfxstate%my_c_atoms(hfxstate%n_my_c_atoms),stat=ierr)
    call utils_alloc_check(myself,'my_c_atoms',ierr)

    ! --------------------------------------------------------------------------
    ! jd: Loop over *all* global atoms B                                     BBB
    ! --------------------------------------------------------------------------
    my_b_atom_idx = 1
    loop_B:                                                                    &
    do global_b = 1, global_nat
       if(b_atoms(global_b)) then
          hfxstate%my_b_atoms(my_b_atom_idx) = global_b
          my_b_atom_idx = my_b_atom_idx + 1
       end if
    end do loop_B

    ! jd: Check we've indeed stored all B atoms from the bitmap
    call utils_assert(my_b_atom_idx-1 == hfxstate%n_my_b_atoms, &
         myself//': Logic error [B]', my_b_atom_idx-1, hfxstate%n_my_b_atoms)

    ! jd: Establish n_my_bb_ngwfs
    hfxstate%n_my_bb_ngwfs = 0
    do my_b_atom_idx = 1, hfxstate%n_my_b_atoms
       n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(&
            hfxstate%my_b_atoms(my_b_atom_idx))
       hfxstate%my_bb_ngwf_offsets(my_b_atom_idx) = hfxstate%n_my_bb_ngwfs + 1
       hfxstate%n_my_bb_ngwfs = hfxstate%n_my_bb_ngwfs + n_ngwfs_b
    end do

    ! --------------------------------------------------------------------------
    ! jd: Loop over *all* global atoms C                                     CCC
    ! --------------------------------------------------------------------------
    my_c_atom_idx = 1
    loop_C:                                                                    &
    do global_c = 1, global_nat
       if(c_atoms(global_c)) then
          hfxstate%my_c_atoms(my_c_atom_idx) = global_c
          my_c_atom_idx = my_c_atom_idx + 1
       end if
    end do loop_C

    ! jd: Check we've indeed stored all C atoms from the bitmap
    call utils_assert(my_c_atom_idx-1 == hfxstate%n_my_c_atoms, &
         myself//': Logic error [C]', my_c_atom_idx-1, hfxstate%n_my_c_atoms)

    ! --------------------------------------------------------------------------
    ! jd: [4] Bitmaps are no longer needed.
    ! --------------------------------------------------------------------------
    deallocate(c_atoms,stat=ierr)
    call utils_dealloc_check(myself,'c_atoms',ierr)
    deallocate(b_atoms,stat=ierr)
    call utils_dealloc_check(myself,'b_atoms',ierr)

    call comms_barrier
    call utils_trace_out(myself)
    call timer_clock(myself,2)

  end subroutine hfx_split_pair_list

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_determine_my_a_atoms(hfxstate, global_nat)
    !==========================================================================!
    ! Determines the list of A atoms for a proc, that is, a list of atoms A    !
    ! who are X-neighbours with my atoms B. Atoms B for the current proc are   !
    ! inferred from my_b_atoms. Results (global indices of A atoms) are        !
    ! stored in my_a_atoms.                                                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (inout): All book-keeping happens here.                       !
    !   global_nat (in): size(elements).                                       !
    !--------------------------------------------------------------------------!
    ! Caller is responsible for freeing my_a_atoms later (this happens in      !
    ! hfx_free_arrays).                                                        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2019.                              !
    !==========================================================================!

    use comms, only: comms_barrier
    use neighbour_list, only: neighbour_list_for_atom_list
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target :: hfxstate
    integer, intent(in)                    :: global_nat

    ! jd: Local variables
    character(len=*), parameter :: myself = 'hfx_determine_my_a_atoms'

    ! -----------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    call neighbour_list_for_atom_list(hfxstate%my_a_atoms, &
         hfxstate%n_my_a_atoms, hfxstate%my_b_atoms, &
         hfxstate%n_my_b_atoms, hfxstate%x_atoms_nl, global_nat)

    call comms_barrier
    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hfx_determine_my_a_atoms

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_determine_my_d_atoms(hfxstate, global_nat)
    !==========================================================================!
    ! Determines the list of D atoms for a proc, that is, a list of atoms D who!
    ! are S-neighbours with my atoms A. My atoms A are inferred from           !
    ! my_a_atoms. Results (global indices of D atoms) are stored in my_d_atoms.!
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (inout): All book-keeping happens here.                       !
    !   global_nat (in): size(elements).                                       !
    !--------------------------------------------------------------------------!
    ! Caller is responsible for freeing my_d_atoms later (this happens in      !
    ! hfx_free_arrays).                                                        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2019.                              !
    !==========================================================================!

    use comms, only: comms_barrier
    use neighbour_list, only: neighbour_list_for_atom_list
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target :: hfxstate
    integer, intent(in)                    :: global_nat

    ! jd: Local variables
    character(len=*), parameter :: myself = 'hfx_determine_my_d_atoms'

    ! -----------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    call neighbour_list_for_atom_list(hfxstate%my_d_atoms, &
         hfxstate%n_my_d_atoms, hfxstate%my_a_atoms, &
         hfxstate%n_my_a_atoms, hfxstate%s_atoms_nl, global_nat)

    call comms_barrier
    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hfx_determine_my_d_atoms

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_determine_dlists(hfxstate, ngwfs_hts, ngwf_indices )
    !==========================================================================!
    ! Determines which D atoms share which PPDs.                               !
    ! For all my atoms B:                                                      !
    !   - determine PPDs in the union of all corresponding A,                  !
    !   - determine which D's share each PPD.                                  !
    ! B atoms are taken from my_b_atoms. D atoms are taken from my_d_atoms.    !
    !                                                                          !
    ! As a result, for every occuring pair of (B, ppd#) a list of all D's that !
    ! share it is stored.                                                      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (inout): All book-keeping happens here.                       !
    !   ngwfs_hts (in): Needed to determine PPD lists of A's for current B.    !
    !                   This is because atoms A are, in gen., not sparse local.!
    !   ngwf_indices (in): Array of HFX_NGWF_INDEX instances containing        !
    !                      pointers to the NGWF basis and the NGWF_REP used    !
    !                      by atoms A, B, C and D (only A and D indices used). !
    !                      for atoms A and B used).                            !
    !--------------------------------------------------------------------------!
    ! Caveat/note:                                                             !
    ! It may seem strange that dlists are indexed with B in addition to a      !
    ! PPD index. This is needed because some A's will be out of X cutoff from  !
    ! B and will therefore not contribute their D's.                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2019.                                 !
    !==========================================================================!

    use comms, only: comms_barrier
    use hash_table, only: hash_table_add, hash_table_purge, hash_table_list
    use remote, only: remote_ppd_list_of_atom
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_abort, &
         utils_alloc_check, utils_dealloc_check

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target  :: hfxstate
    type(HT_HASH_TABLE), intent(in), target :: ngwfs_hts(:)
    type(HFX_NGWF_INDEX), intent(in)        :: ngwf_indices(A_NGWF:D_NGWF)

    ! jd: Local variables
    integer, allocatable :: a_ppd_list(:) ! PPDs on all A's *relevant* to
    integer :: n_a_ppd_list               ! current B. Not "all my A PPDs".
    integer :: atom_d_ppd_list(ngwf_indices(D_NGWF)%basis%max_n_ppds_sphere)
    integer, allocatable :: delems(:)   ! number of D's in corresponding ielems
    integer, allocatable :: dlist(:,:)  ! PPD indices of D's sharing PPDs with given B
    integer :: n_atom_d_ppds
    integer :: b_idx
    integer :: global_b
    integer :: d_idx
    integer :: global_d
    integer :: i_ppd
    integer :: cur_ppd
    integer :: i_d
    integer :: ierr
    character(len=*), parameter :: myself = 'hfx_determine_dlists'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    ! jd: Purge any prior DLISTS
    call hash_table_list(hfxstate%dlists_ht,0)
    call hash_table_purge(hfxstate%dlists_ht)

    ! --------------------------------------------------------------------------
    ! jd: Loop over B's given to this rank                                MY BBB
    ! --------------------------------------------------------------------------
    loop_B:                                                                   &
    do b_idx = 1, hfxstate%n_my_b_atoms
       global_b = hfxstate%my_b_atoms(b_idx)

       call hfx_determine_a_ppd_list_for_b(a_ppd_list, n_a_ppd_list, &
            ngwfs_hts, hfxstate, ngwf_indices(A_NGWF), global_b)

       allocate(delems(n_a_ppd_list),stat=ierr)
       call utils_alloc_check(myself,'delems',ierr)
       allocate(dlist(n_a_ppd_list,max_num_at_per_ppd),stat=ierr)
       call utils_alloc_check(myself,'dlist',ierr)

       delems(:) = 0
       dlist(:,:) = garbage_int

       ! --------------------------------------------------------------------
       ! jd: For all my D's                                               DDD
       ! --------------------------------------------------------------------
       loop_D:                                                              &
       do d_idx = 1, hfxstate%n_my_d_atoms
          global_d = hfxstate%my_d_atoms(d_idx)

          ! jd: Obtain PPD indices on D (not necessarily local)
          call remote_ppd_list_of_atom(atom_d_ppd_list, n_atom_d_ppds, &
               global_d, ngwf_indices(D_NGWF)%basis, &
               ngwf_indices(D_NGWF)%rep%ngwf_cache_handle)

          ! ********************************************************************
          ! * For each of A PPDs determine in which D it's going to appear
          ! ********************************************************************
          ! --------------------------------------------------------------------
          ! jd: For all PPDs in union of A's relevant to current B           PPD
          ! --------------------------------------------------------------------
          loop_PPD:                                                            &
          do i_ppd = 1, n_a_ppd_list
             cur_ppd = a_ppd_list(i_ppd)

             if(any(atom_d_ppd_list(1:n_atom_d_ppds) == cur_ppd)) then
                ! jd: This D shares cur_ppd with some A that is X-neighbour of this B
                delems(i_ppd) = delems(i_ppd) + 1
                i_d = delems(i_ppd)
                if(i_d > max_num_at_per_ppd) then
                   call utils_abort(&
                        myself//': Insufficient max_num_at_per_ppd', i_d)
                end if
                dlist(i_ppd,i_d) = global_d
             end if

          end do loop_PPD

       end do loop_D

       ! -----------------------------------------------------------------------
       ! jd: For all PPDs in union of A's relevant to current B              PPD
       ! -----------------------------------------------------------------------
       loop_PPD2:                                                              &
       do i_ppd = 1, n_a_ppd_list
          cur_ppd = a_ppd_list(i_ppd)
          call hash_table_add(hfxstate%dlists_ht, &
               real(dlist(i_ppd,1:delems(i_ppd)),kind=DP), delems(i_ppd),&
               global_b, cur_ppd, overfill_strategy = 'F')
       end do loop_PPD2

       deallocate(a_ppd_list,stat=ierr)
       call utils_dealloc_check(myself,'a_ppd_list',ierr, &
            allocated_in='hfx_determine_a_ppd_list_for_b')
       deallocate(dlist,stat=ierr)
       call utils_dealloc_check(myself,'dlist',ierr)
       deallocate(delems,stat=ierr)
       call utils_dealloc_check(myself,'delems',ierr)

    end do loop_B

    call comms_barrier
    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hfx_determine_dlists

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_determine_alists(hfxstate, ngwfs_hts, alpha_ngwf_index)
    !==========================================================================!
    ! Determines which A atoms share which PPDs.                               !
    ! For all my atoms A:                                                      !
    !   - determine PPDs in the union of all corresponding A,                  !
    !   - determine which A's share each PPD.                                  !
    ! A atoms are taken from my_a_atoms.                                       !
    !                                                                          !
    ! As a result, for every ppd# that occurs in my A, a list of all A's that  !
    ! share it is stored.                                                      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (inout): All book-keeping happens here.                       !
    !   ngwfs_hts (in): Needed to determine PPD lists of A's.                  !
    !                   This is because atoms A are, in gen., not sparse local.!
    !   alpha_ngwf_index (in): The HFX_NGWF_INDEX instance containing the      !
    !                          basis and the NGWF_REP used by A atoms.         !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2019.                                 !
    !==========================================================================!

    use comms, only: comms_barrier
    use hash_table, only: hash_table_add, hash_table_purge, hash_table_list
    use remote, only: remote_ppd_list_of_atom
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_abort, &
         utils_alloc_check, utils_dealloc_check

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target  :: hfxstate
    type(HT_HASH_TABLE), intent(in), target :: ngwfs_hts(:)
    type(HFX_NGWF_INDEX), intent(in)        :: alpha_ngwf_index

    ! jd: Local variables
    integer :: atom_a_ppd_list(alpha_ngwf_index%basis%max_n_ppds_sphere)
    integer, allocatable :: aelems(:)   ! number of A's in corresponding alist
    integer, allocatable :: alist(:,:)  ! PPD indices of A's sharing PPDs
    integer :: n_atom_a_ppds
    integer :: a_idx
    integer :: global_a
    integer :: i_ppd
    integer :: cur_ppd
    integer :: i_a
    integer :: ierr
    character(len=*), parameter :: myself = 'hfx_determine_alists'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    ! jd: Purge any prior ALISTS
    call hash_table_list(hfxstate%alists_ht,0)
    call hash_table_purge(hfxstate%alists_ht)

    allocate(aelems(hfxstate%n_my_a_ppd_list),stat=ierr)
    call utils_alloc_check(myself,'aelems',ierr)
    allocate(alist(hfxstate%n_my_a_ppd_list,max_num_at_per_ppd),stat=ierr)
    call utils_alloc_check(myself,'alist',ierr)

    aelems(:) = 0
    alist(:,:) = garbage_int

    ! --------------------------------------------------------------------------
    ! jd: Loop over A's given to this rank                                MY AAA
    ! --------------------------------------------------------------------------
    loop_A:                                                                    &
    do a_idx = 1, hfxstate%n_my_a_atoms
       global_a = hfxstate%my_a_atoms(a_idx)

       ! jd: Obtain PPD indices on A (not necessarily local)
       call remote_ppd_list_of_atom(atom_a_ppd_list, n_atom_a_ppds, &
            global_a, alpha_ngwf_index%basis, &
            alpha_ngwf_index%rep%ngwf_cache_handle)

       ! -----------------------------------------------------------------------
       ! jd: For all PPDs in union of my A's                                 PPD
       ! -----------------------------------------------------------------------
       loop_PPD:                                                               &
       do i_ppd = 1, hfxstate%n_my_a_ppd_list
          cur_ppd = hfxstate%my_a_ppd_list(i_ppd)

          if(any(atom_a_ppd_list(1:n_atom_a_ppds) == cur_ppd)) then
             ! jd: This A has cur_ppd
             aelems(i_ppd) = aelems(i_ppd) + 1
             i_a = aelems(i_ppd)
             if(i_a > max_num_at_per_ppd) then
                call utils_abort(&
                     myself//': Insufficient max_num_at_per_ppd', i_a)
             end if
             alist(i_ppd,i_a) = global_a
          end if

       end do loop_PPD

    end do loop_A

    ! --------------------------------------------------------------------------
    ! jd: For all PPDs in union of my A's                                    PPD
    ! --------------------------------------------------------------------------
    loop_PPD2:                                                                 &
    do i_ppd = 1, hfxstate%n_my_a_ppd_list
       cur_ppd = hfxstate%my_a_ppd_list(i_ppd)
       call hash_table_add(hfxstate%alists_ht, &
            real(alist(i_ppd,1:aelems(i_ppd)),kind=DP), aelems(i_ppd),&
            cur_ppd, overfill_strategy = 'F')
    end do loop_PPD2

    deallocate(alist,stat=ierr)
    call utils_dealloc_check(myself,'alist',ierr)
    deallocate(aelems,stat=ierr)
    call utils_dealloc_check(myself,'aelems',ierr)

    call comms_barrier
    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hfx_determine_alists

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_determine_ppd_list(hfxstate, ngwfs_hts, &
       alpha_ngwf_index)
    !==========================================================================!
    ! Determines the list of PPD inidices for all my A atoms.                  !
    ! 'My' A atoms are the atoms that are X-neighbours to my B atoms. They are !
    ! not, in general, sparse-local to us.                                     !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (inout): List of A atoms is read from here, my_a_ppd_list(:)  !
    !                     member is allocated and populated.                   !
    !   ngwfs_hts (in): NGWF caches. PPD lists of atoms that are in general    !
    !                   remote to this proc are obtained from here.            !
    !   alpha_ngwf_index (in): The HFX_NGWF_INDEX instance containing the      !
    !                          basis and the NGWF_REP used by A atoms.         !
    !                          Needed only to get max_n_ppds_sphere and to     !
    !                          find out who sparse owns what and get the       !
    !                          NGWF set cache handle.                          !
    !--------------------------------------------------------------------------!
    ! Caller is responsible for freeing my_a_ppd_list later (happens in        !
    ! hfx_free_arrays).                                                        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2019.                              !
    !==========================================================================!

    use comms, only: comms_barrier
    use hash_table, only: HT_HASH_TABLE
    use ppd_ops, only: MAX_PPDS_OF_INTEREST, ppd_set_union
    use remote, only: remote_ppd_list_of_atom
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_alloc_check

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target  :: hfxstate
    type(HT_HASH_TABLE), intent(in), target :: ngwfs_hts(:)
    type(HFX_NGWF_INDEX), intent(in)        :: alpha_ngwf_index

    ! jd: Local variables
    integer :: my_a_atom_idx
    integer :: global_a
    integer :: atom_a_ppd_list(alpha_ngwf_index%basis%max_n_ppds_sphere)
    integer :: n_atom_a_ppds
    integer :: local_ppd_list(1:MAX_PPDS_OF_INTEREST)
    integer :: n_local_ppds
    integer :: ierr
    character(len=*), parameter :: myself = 'hfx_determine_ppd_list'

    ! -----------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    local_ppd_list(:) = garbage_int
    n_local_ppds = 0

    ! --------------------------------------------------------------------------
    ! jd: Loop over A's                                                   MY AAA
    ! --------------------------------------------------------------------------
    loop_A:                                                                    &
    do my_a_atom_idx = 1, hfxstate%n_my_a_atoms
       global_a = hfxstate%my_a_atoms(my_a_atom_idx)

       ! jd: Obtain PPD indices on A (not necessarily local)
       call remote_ppd_list_of_atom(atom_a_ppd_list, n_atom_a_ppds, &
            global_a, alpha_ngwf_index%basis, &
            alpha_ngwf_index%rep%ngwf_cache_handle)

       ! jd: Accumulate the PPD indices in a local workspace with fixed size.
       call ppd_set_union(local_ppd_list, n_local_ppds, &
            atom_a_ppd_list, n_atom_a_ppds)

    end do loop_A

    ! jd: Copy over the filled PPD indices to result
    hfxstate%n_my_a_ppd_list = n_local_ppds
    allocate(hfxstate%my_a_ppd_list(hfxstate%n_my_a_ppd_list),stat=ierr)
    call utils_alloc_check(myself,'my_a_ppd_list',ierr)
    hfxstate%my_a_ppd_list(:) = local_ppd_list(1:hfxstate%n_my_a_ppd_list)

    call comms_barrier
    call utils_trace_out(myself)
    call timer_clock(myself,2)

  end subroutine hfx_determine_ppd_list

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_determine_a_ppd_list_for_b(a_ppd_list, n_a_ppd_list, &
       ngwfs_hts, hfxstate, alpha_ngwf_index, global_b)
    !==========================================================================!
    ! Determines the list of PPDs that is a union over all A's that are        !
    ! X-neighbours with a given B.                                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   a_ppd_list (out): The list of PPD indices is returned here. This array !
    !                     is allocated here.                                   !
    !   n_a_ppd_list (out): The number of elements in the returned a_ppd_list. !
    !   ngwfs_hts (in): NGWF caches. Needed, because atom A is, in general,    !
    !                   not sparse-local.                                      !
    !   hfxstate (in): Needed for its X neighbour list.                        !
    !   alpha_ngwf_index (in): The HFX_NGWF_INDEX instance containing the      !
    !                          basis and the NGWF_REP used by A atoms.         !
    !   global_b (in): Global atom index for B whose A's' PPDs we want.        !
    !--------------------------------------------------------------------------!
    ! Caller is responsible for freeing a_ppd_list later.                      !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2019.                                 !
    !==========================================================================!

    use hash_table, only: HT_HASH_TABLE
    use ppd_ops, only: MAX_PPDS_OF_INTEREST, ppd_set_union
    use remote, only: remote_ppd_list_of_atom
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_assert, &
         utils_alloc_check

    implicit none

    ! jd: Arguments
    integer, allocatable, intent(out)       :: a_ppd_list(:)
    integer, intent(out)                    :: n_a_ppd_list
    type(HT_HASH_TABLE), intent(in), target :: ngwfs_hts(:)
    type(HFX_STATE), intent(in), target     :: hfxstate
    type(HFX_NGWF_INDEX), intent(in)        :: alpha_ngwf_index
    integer, intent(in)                     :: global_b

    ! jd: Local variables
    integer :: a_idx
    integer :: global_a
    integer :: atom_a_ppd_list(alpha_ngwf_index%basis%max_n_ppds_sphere)
    integer :: n_atom_a_ppds
    integer :: local_ppd_list(1:MAX_PPDS_OF_INTEREST)
    integer :: n_local_ppds
    integer :: ierr
    character(len=*), parameter :: myself = 'hfx_determine_a_ppd_list_for_b'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)
    ! jd: No barriers here, each rank calls this a different number of times
    call timer_clock(myself,1)

    ! jd: Sanity check on neighbour list
    call utils_assert(hfxstate%x_atoms_nl%populated, &
         myself//': x_atoms neighbour list not populated.')

    local_ppd_list(:) = garbage_int
    n_local_ppds = 0

    ! --------------------------------------------------------------------------
    ! jd: Loop over A's that are X-neighbours of B                           AAA
    ! --------------------------------------------------------------------------
    loop_A:                                                                    &
    do a_idx = hfxstate%x_atoms_nl%first_idx(global_b), &
         hfxstate%x_atoms_nl%last_idx(global_b)
       global_a = hfxstate%x_atoms_nl%neighbours(a_idx)

       ! jd: Obtain PPD indices on A (not necessarily local)
       call remote_ppd_list_of_atom(atom_a_ppd_list, n_atom_a_ppds, &
            global_a, alpha_ngwf_index%basis, &
            alpha_ngwf_index%rep%ngwf_cache_handle)

       ! jd: Accumulate the PPD indices in a local workspace with fixed size.
       call ppd_set_union(local_ppd_list, n_local_ppds, &
            atom_a_ppd_list, n_atom_a_ppds)

    end do loop_A

    ! jd: Copy over the filled PPD indices to result
    n_a_ppd_list = n_local_ppds
    allocate(a_ppd_list(n_a_ppd_list),stat=ierr)
    call utils_alloc_check(myself,'a_ppd_list',ierr)
    a_ppd_list(:) = local_ppd_list(1:n_a_ppd_list)

    ! jd: No barriers here, each rank calls this a different number of times
    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hfx_determine_a_ppd_list_for_b

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_comms_ngwfs(ngwfs_hts, wanted_atoms, n_wanted_atoms, &
       wanted_ngwf_index, ireg)
    !==========================================================================!
    ! Communicates NGWFs across MPI ranks. Every MPI rank gets NGWFs for atoms !
    ! whose (global) indices are provided in 'wanted_atoms'. Every MPI rank    !
    ! serves NGWFs of atoms it owns locally.                                   !
    ! Results (packed NGWFs) are stored in 'ngwfs_hts', at the index corres-   !
    ! ponding to the NGWF set (cache_handle) being exchanged, determined by    !
    ! rep's ngwf_cache_handle. Everything happens in embedding region ireg.    !
    ! 'ngwfs_hts' would normally be remote's NGWF caches.                      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   ngws_hts (inout): Caches for NGWFs. Results stored here.               !
    !   wanted_atoms (in): Array of global atom indices whose NGWF we'd like   !
    !                      to obtain from whoever has them.                    !
    !   n_wanted_atoms (in): Number of elements in the above.                  !
    !   wanted_ngwf_index (in): The HFX_NGWF_INDEX instance containing the     !
    !                           basis from which we want the NGWFs and the     !
    !                           NGWF_REP in which the NGWFs live.              !
    !   ireg (in): Embedding region identifier.                                !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2019.                              !
    !==========================================================================!

    use comms, only: pub_total_num_procs, pub_null_handle, comms_barrier, &
         comms_send, comms_wait
    use hash_table, only: HT_HASH_TABLE
    use remote, only: NGWF_REQUEST_TAG, remote_obtain_ngwf, remote_serve_ngwfs
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    type(HT_HASH_TABLE), target, intent(inout) :: ngwfs_hts(:)
    integer, intent(in)                        :: wanted_atoms(:)
    integer, intent(in)                        :: n_wanted_atoms
    type(HFX_NGWF_INDEX), intent(in)           :: wanted_ngwf_index
    integer, intent(in)                        :: ireg

    ! jd: Local variables
    integer :: wanted_idx
    integer :: wanted_atom_e
    integer :: ngwf_e
    integer :: first_ngwf_idx_of_e
    integer :: global_ee_ngwf_idx
    integer :: proc
    logical :: done
    real(kind=DP), pointer :: dummy_ngwf_ee(:) ! remote
    integer :: send_handles(0:pub_total_num_procs-1)
    logical :: who_is_done(0:pub_total_num_procs-1)
    character(len=*), parameter :: myself = 'hfx_comms_ngwfs'

    ! -----------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    send_handles(:) = pub_null_handle
    who_is_done(:) = .false.

    ! --------------------------------------------------------------------------
    ! jd: All wanted atoms E                                                 EEE
    ! --------------------------------------------------------------------------
    loop_E:                                                                    &
    do wanted_idx = 1, n_wanted_atoms
       wanted_atom_e = wanted_atoms(wanted_idx)
       first_ngwf_idx_of_e = &
           wanted_ngwf_index%basis%first_on_atom(wanted_atom_e)

       ! -----------------------------------------------------------------------
       ! jd: All NGWFs e of the wanted atom E                                eee
       ! -----------------------------------------------------------------------
       loop_ngwf_e:                                                            &
       do ngwf_e = 1, wanted_ngwf_index%basis%num_on_atom(wanted_atom_e)
          global_ee_ngwf_idx = first_ngwf_idx_of_e + ngwf_e - 1
          ! jd: Ignore the actual NGWF, just make sure it becomes cached
          call remote_obtain_ngwf(dummy_ngwf_ee, ngwfs_hts, who_is_done, &
               global_ee_ngwf_idx, &
               wanted_ngwf_index%rep%ngwfs_on_grid(ireg), &
               wanted_ngwf_index%basis, &
               wanted_ngwf_index%rep%ngwf_cache_handle, &
               user_tag = wanted_ngwf_index%rep%ngwf_cache_handle, &
               overfill_strategy = 'F')
       end do loop_ngwf_e

    end do loop_E

    ! jd: Finish NGWF comms
    if(pub_total_num_procs > 1) then
       ! jd: After we're done, first notify everyone, then keep serving
       !     other procs with NGWFs, until everyone is done
       do proc = 0, pub_total_num_procs-1
          call comms_send(proc, -1, 1, tag = NGWF_REQUEST_TAG + 1, & ! [x]
               return_handle = send_handles(proc), add_to_stack = .false.)
               ! [x]: For a 'done' request it's OK to only send it for
               !      NGWF set 1 -- it's understood that *all* sets finished
               !      participating.
       end do

       done = .false.
       do while(.not. done)
          call remote_serve_ngwfs(done, who_is_done, &
               wanted_ngwf_index%rep%ngwfs_on_grid(ireg), &
               wanted_ngwf_index%basis, &
               wanted_ngwf_index%rep%ngwf_cache_handle)
       end do

       do proc = 0, pub_total_num_procs-1
          call comms_wait(send_handles(proc))
       end do

       call comms_barrier
    end if

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hfx_comms_ngwfs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_comms_dkn(packed_dkn_blocks_ht, &
       atom_list, n_atoms, nl, denskern, row_max_on_atom, col_max_on_atom, &
       opt_transpose)
    !==========================================================================!
    ! Communicates necessary DKN blocks between MPI ranks.                     !
    ! 'atom_list' is an array 'n_atoms' long, which specifies rows that are    !
    ! going to be needed on the calling MPI proc. Needed columns are establi-  !
    ! shed from the neighbour list 'nl'. Source data is taken from 'denskern', !
    ! results wind up in packed_dkn_blocks_ht.                                 !
    !                                                                          !
    ! This routine is used to ensure all K^CcDd where C are my atoms and D are !
    ! their s-x-S-neighbours are available locally.                            !
    ! Similarly, we use it for K^BbAa where B are my atoms, and A are their    !
    ! X-neighbours.                                                            !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   packed_dkn_cd_blocks_ht (inout): Results go here.                      !
    !   atom_list (in): Global indices of atoms whose DKN rows we need.        !
    !   n_atoms (in): Number of elements in the above.                         !
    !   nl (in): Neighbour list specifying which columns are going to be       !
    !            needed for each row.                                          !
    !   denskern (in) : Density kernel matrix to be communicated.              !
    !   row_max_on_atom (in): Maximum number of NGWFs in a row.                !
    !   col_max_on_atom (in): Maximum number of NGWFs in a column.             !
    !   opt_transpose (in, opt): If present and .true., dimensions and indices !
    !                            for rows and columns will be transposed.      !
    !                            Mind that it will still be rows that are de-  !
    !                            fined by 'atom_list' and it will still be     !
    !                            columns that are defined thy 'nl'.            !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2019.                              !
    ! Generalised by Jacek Dziedzic in April 2020.                             !
    !==========================================================================!

    use comms, only: pub_total_num_procs, pub_null_handle, comms_send, &
         comms_wait, comms_barrier
    use hash_table, only: HT_HASH_TABLE
    use neighbour_list, only: NL_NEIGHBOUR_LIST
    use remote, only: DKN_BLOCK_REQUEST_TAG, remote_obtain_dkn_block, &
         remote_serve_dkn_blocks
    use rundat, only: pub_num_spins
    use sparse, only: SPAM3
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_assert

    implicit none

    ! jd: Arguments
    type(HT_HASH_TABLE), intent(inout), target :: packed_dkn_blocks_ht
    integer, intent(in)                        :: atom_list(:)
    integer, intent(in)                        :: n_atoms
    type(NL_NEIGHBOUR_LIST), intent(in)        :: nl
    type(SPAM3), intent(in)                    :: denskern(pub_num_spins)
    integer, intent(in)                        :: row_max_on_atom
    integer, intent(in)                        :: col_max_on_atom
    logical, intent(in), optional              :: opt_transpose

    ! jd: Local variables
    real(kind=DP), pointer   :: dummy_packed_dkn_block(:)
    logical                  :: who_is_done(0:pub_total_num_procs-1)
    integer                  :: send_handles(0:pub_total_num_procs-1,2)
    integer :: my_row_atom_idx
    integer :: global_row
    integer :: global_col
    integer :: col_idx
    logical :: done
    logical :: loc_transpose
    integer :: proc
    character(len=*), parameter :: myself = 'hfx_comms_dkn'

    ! -----------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    ! jd: Sanity check on neighbour list
    call utils_assert(nl%populated, myself//': neighbour list not populated.')

    if(present(opt_transpose)) then
       loc_transpose = opt_transpose
    else
       loc_transpose = .false.
    end if

    who_is_done(:) = .false.
    send_handles(:,:) = pub_null_handle

    ! --------------------------------------------------------------------------
    ! jd: Loop over row atoms                                                ROW
    ! --------------------------------------------------------------------------
    loop_row:                                                                    &
    do my_row_atom_idx = 1, n_atoms
       global_row = atom_list(my_row_atom_idx)

       ! -----------------------------------------------------------------------
       ! jd: Loop over col atoms that are neighbours of row                  COL
       ! -----------------------------------------------------------------------
       loop_col:                                                               &
       do col_idx = nl%first_idx(global_row), nl%last_idx(global_row)
          global_col = nl%neighbours(col_idx)

          ! jd: We don't really care about the data passed by pointer, only
          !     about the DKN block getting to the cache
          call remote_obtain_dkn_block(dummy_packed_dkn_block, & ! out
               packed_dkn_blocks_ht, who_is_done, &
               merge(global_col, global_row, loc_transpose), &
               merge(global_row, global_col, loc_transpose), &
               denskern, &
               merge(col_max_on_atom, row_max_on_atom, loc_transpose), &
               merge(row_max_on_atom, col_max_on_atom, loc_transpose), &
               overfill_strategy = 'F')

       end do loop_col
    end do loop_row

    ! jd: After we're done, first notify everyone, then keep serving
    !     other procs with DKN blocks, until everyone is done
    if(pub_total_num_procs > 1) then
       do proc = 0, pub_total_num_procs-1
          ! jd: Two notifications are needed for DKN blocks, they expect 2 sends
          call comms_send(proc, -1, 1, tag = DKN_BLOCK_REQUEST_TAG + 1, &
               return_handle = send_handles(proc,1), add_to_stack = .false.)
          call comms_send(proc, -1, 1, tag = DKN_BLOCK_REQUEST_TAG + 1, &
               return_handle = send_handles(proc,2), add_to_stack = .false.)
       end do

       done = .false.
       do while(.not. done)
          call remote_serve_dkn_blocks(done, who_is_done, denskern, &
               merge(col_max_on_atom, row_max_on_atom, loc_transpose), &
               merge(row_max_on_atom, col_max_on_atom, loc_transpose))
       end do

       do proc = 0, pub_total_num_procs-1
          call comms_wait(send_handles(proc,1))
          call comms_wait(send_handles(proc,2))
       end do

    end if

    call comms_barrier
    call utils_trace_out(myself)
    call timer_clock(myself,2)

  end subroutine hfx_comms_dkn

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_comms_metric(hfxstate, full_metric_matrix, which_metric)
    !==========================================================================!
    ! Communicates necessary metric matrix blocks between MPI ranks.           !
    ! Ensures all {V,O}_BbCc where B are my atoms and C are their S-neighbours !
    ! wind up in packed_metric_matrix_blocks_ht. In general neither of these   !
    ! are sparse local, so we request them from remote procs.                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (inout): List of my BC pairs is extracted from here.          !
    !                     Results go into packed_metric_matrix_blocks_ht.      !
    !   full_metric_matrices (in) : Metric matrices V and O in SPAM3 rep.      !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in April 2019 using hfx_comms_dkn() as         !
    ! template.                                                                !
    !==========================================================================!

    use comms, only: pub_total_num_procs, pub_null_handle, comms_send, &
         comms_wait, comms_barrier
    use hash_table, only: hash_table_list
    use remote, only: MATRIX_BLOCK_REQUEST_TAG, remote_obtain_matrix_block, &
         remote_serve_matrix_blocks
    use sparse, only: SPAM3
    use sw_resolution_of_identity, only: swri_library
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target     :: hfxstate
    type(SPAM3), intent(in)                    :: full_metric_matrix
    integer, intent(in)                        :: which_metric

    ! jd: Local variables
    real(kind=DP), pointer :: dummy_packed_matrix_block(:)
    logical                :: who_is_done(0:pub_total_num_procs-1)
    integer                :: send_handles(0:pub_total_num_procs-1,2)
    integer :: global_b
    integer :: global_c
    integer :: pair_idx
    logical :: done
    integer :: proc
    character(len=*), parameter :: myself = 'hfx_comms_metric'

    ! -----------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    who_is_done(:) = .false.
    send_handles(:,:) = pub_null_handle

    ! --------------------------------------------------------------------------
    ! jd: Loop over my BC pairs                                               BC
    ! --------------------------------------------------------------------------
    loop_BC:                                                                   &
    do pair_idx = 1, hfxstate%n_my_bc_pairs
       global_b = hfxstate%my_bc_pairs(1,pair_idx)
       global_c = hfxstate%my_bc_pairs(2,pair_idx)

       ! jd: We don't really care about the data passed by pointer, only
       !     about the matrix block getting to the cache.
       call remote_obtain_matrix_block(dummy_packed_matrix_block, & ! out
            hfxstate%packed_metric_matrix_blocks_hts(which_metric), &
            who_is_done, global_b, global_c, full_metric_matrix, &
            swri_library(hfxstate%hfx_swri_h)%quality%num_sws_per_centre, &
            swri_library(hfxstate%hfx_swri_h)%quality%num_sws_per_centre, &
            overfill_strategy = 'F')

    end do loop_BC

    ! jd: After we're done, first notify everyone, then keep serving
    !     other procs with metric matrix blocks, until everyone is done
    if(pub_total_num_procs > 1) then
       do proc = 0, pub_total_num_procs-1
          ! jd: Two notifications are needed for matrix blocks, they expect 2 sends
          call comms_send(proc, -1, 1, tag = MATRIX_BLOCK_REQUEST_TAG, &
               return_handle = send_handles(proc,1), add_to_stack = .false.)
          call comms_send(proc, -1, 1, tag = MATRIX_BLOCK_REQUEST_TAG, &
               return_handle = send_handles(proc,2), add_to_stack = .false.)
       end do

       done = .false.
       do while(.not. done)
          call remote_serve_matrix_blocks(done, who_is_done, full_metric_matrix, &
               swri_library(hfxstate%hfx_swri_h)%quality%num_sws_per_centre, &
               swri_library(hfxstate%hfx_swri_h)%quality%num_sws_per_centre)
       end do

       do proc = 0, pub_total_num_procs-1
          call comms_wait(send_handles(proc,1))
          call comms_wait(send_handles(proc,2))
       end do

    end if

    call hash_table_list(hfxstate%packed_metric_matrix_blocks_hts(which_metric),0)

    call comms_barrier
    call utils_trace_out(myself)
    call timer_clock(myself,2)

  end subroutine hfx_comms_metric

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_comms_xmatrix(xmatrix, &                   ! inout, adds to
       x_ab_contributions_ht, ngwf_indices)                 ! in
    !==========================================================================!
    ! Assembles the X matrix from contributions on different procs.            !
    ! Contributions are taken from x_ab_contributions_ht, which contains       !
    ! atomblocks indexed by global atom indices A and B. Multiple MPI ranks    !
    ! can own contributions with the same A and B (they have to be summed).    !
    ! The idea is to get each contribution to the proc that sparse-owns A and  !
    ! to sum over potentially multiple contributions to the same atomblock.    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   xmatrix (inout): Results are *added* to here.                          !
    !   x_ab_contributions_ht (in): My contributions to X are taken from here. !
    !   ngwf_indices (in): Array of HFX_NGWF_INDEX instances containing        !
    !                      pointers to the NGWF basis and the NGWF_REP used    !
    !                      by atoms A, B, C and D (only FUNC_BASIS instances   !
    !                      for atoms A and B used).                            !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   This subroutine *adds* to xmatrix, so caller is responsible for        !
    !   zeroing xmatrix first.                                                 !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2019.                                 !
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_any_source, pub_total_num_procs, &
         pub_null_handle, comms_irecv, comms_send, comms_barrier, comms_wait, &
         comms_test, comms_cancel
    use hash_table, only: HT_HASH_TABLE, HT_STORAGE_NODE, hash_table_iterate
    use rundat, only: pub_num_spins
    use sparse, only: SPAM3, sparse_put_block
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_alloc_check, &
         utils_dealloc_check, utils_assert

    implicit none

    ! jd: Arguments
    type(SPAM3), intent(inout)              :: xmatrix(:) ! :spins
    type(HT_HASH_TABLE), intent(in), target :: x_ab_contributions_ht
    type(HFX_NGWF_INDEX), intent(in)        :: ngwf_indices(A_NGWF:D_NGWF)

    ! jd: Local variables
    real(kind=DP), allocatable :: recved_xatomblocks(:,:)
    integer :: first_ngwf_idx_of_a
    integer :: global_a
    integer :: n_ngwfs_a
    integer :: n_ngwfs_b
    integer :: global_b
    integer :: recved_global_a
    integer :: recved_global_b
    integer :: i_batch
    integer :: ndata
    type(HT_STORAGE_NODE), pointer :: ptr
    real(kind=DP), pointer :: data_ptr(:)
    integer :: slot
    integer :: key3, key4, key5
    logical :: eot
    logical :: recv_finished
    logical :: send_finished
    logical :: i_am_done_sending
    integer :: offset
    integer :: recved_nelements
    integer :: source
    integer :: max_msg_size
    integer :: sparse_owner
    integer :: i_rank
    integer :: is
    integer :: ierr
    real(kind=DP)      :: done(3,0:pub_total_num_procs-1)
    integer, parameter :: recv_batch_size = 1000
    integer, parameter :: send_batch_size = 1000
    integer            :: recv_handles(1:recv_batch_size)
    integer            :: send_handles(1:send_batch_size)
    integer            :: send_done_handles(0:pub_total_num_procs-1)
    integer            :: nsent(0:pub_total_num_procs-1)
    integer            :: nrecv(0:pub_total_num_procs-1)
    integer            :: nrecv_expected(0:pub_total_num_procs-1)
    logical            :: recvs_in_progress(1:recv_batch_size)
    logical            :: sends_in_progress(1:send_batch_size)
    logical            :: done_procs(0:pub_total_num_procs-1)
    character(len=*), parameter :: myself = 'hfx_comms_xmatrix'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    ! jd: Starts iteration from the 1st record of the HT
    ptr => NULL()
    slot = 0

#ifdef MPI
    max_msg_size = 3 + &
         ngwf_indices(A_NGWF)%basis%max_on_atom*&
         ngwf_indices(B_NGWF)%basis%max_on_atom*pub_num_spins
    allocate(recved_xatomblocks(max_msg_size, recv_batch_size), stat=ierr)
    call utils_alloc_check(myself,'recved_xatomblocks',ierr)

    recv_handles(:) = pub_null_handle
    send_handles(:) = pub_null_handle
    send_done_handles(:) = pub_null_handle
    i_am_done_sending = .false.
    done_procs(:) = .false.
    nsent(:) = 0
    nrecv(:) = 0
    nrecv_expected(:) = -1
    done(:,:) = -1
    recvs_in_progress(:) = .false.
    sends_in_progress(:) = .false.

    ! jd: Loop until everyone is done, ie. we get a 'done sending' notification
    !     from every sender and we received as many messages from every sender
    !     as they claim they sent.
    do while(.not. (all(done_procs(:)) .and. all(nrecv(:) == nrecv_expected(:))))

       ! jd: Post a batch of receives
       ! --- recv --- recv --- recv --- recv --- recv --- recv --- recv ---
       do i_batch = 1, recv_batch_size
          if(.not. recvs_in_progress(i_batch)) then
             call comms_irecv(pub_any_source, recved_xatomblocks(:,i_batch), &
                  handle=recv_handles(i_batch))
             recvs_in_progress(i_batch) = .true.
          end if
       end do

       ! --- send --- send --- send --- send --- send --- send --- send ---
       if(.not. i_am_done_sending) then
          ! jd: Send a batch of records from our HT
          do i_batch = 1, send_batch_size
             if(sends_in_progress(i_batch)) cycle ! this send slot is being used
             ! jd: Look for next record to send
             call hash_table_iterate(ndata, ptr, slot, x_ab_contributions_ht, &
                  global_a, global_b, key3, key4, key5, eot, data_ptr, .false.)
             if(.not. eot) then
                ! jd: Have a record to send
                first_ngwf_idx_of_a = &
                     ngwf_indices(A_NGWF)%basis%first_on_atom(global_a)
                sparse_owner = &
                     ngwf_indices(A_NGWF)%basis%proc_of_func(first_ngwf_idx_of_a)
                call comms_send(sparse_owner, data_ptr, ndata, &
                     return_handle = send_handles(i_batch), add_to_stack = .false.)
                sends_in_progress(i_batch) = .true.
                nsent(sparse_owner) = nsent(sparse_owner) + 1
             else
                ! jd: I am done sending
                call utils_assert(sum(nsent(:)) == x_ab_contributions_ht%n_slots + &
                     x_ab_contributions_ht%n_chained, myself//&
                     ': Mismatch between number of stripes to be sent and &
                     & number of stripes actually sent.', sum(nsent(:)), &
                     x_ab_contributions_ht%n_slots, &
                     x_ab_contributions_ht%n_chained, pub_my_proc_id)
                ! jd: No more records, we iterated through entire HT.
                !     Send a 'done' message to everyone, tell them how many
                !     messages were sent.
                i_am_done_sending = .true.
                do i_rank = 0, pub_total_num_procs-1
                   done(3,i_rank) = nsent(i_rank)
                   call comms_send(i_rank, done(:,i_rank), 3, return_handle = &
                        send_done_handles(i_rank), add_to_stack = .false.)
                end do
                exit
             end if
          end do
       end if

       ! jd: Check the receives for completion
       ! --- recv --- recv --- recv --- recv --- recv --- recv --- recv ---
       do i_batch = 1, recv_batch_size
          if(.not. recvs_in_progress(i_batch)) cycle
          call comms_test(recv_finished, recv_handles(i_batch), source)
          if(recv_finished) then
             ! jd: Receive completed
             recved_global_a = nint(recved_xatomblocks(1,i_batch))
             recved_global_b = nint(recved_xatomblocks(2,i_batch))
             recvs_in_progress(i_batch) = .false.
             if(recved_global_a == -1 .and. recved_global_b == -1) then
                ! jd: Received a done notification. Remember this proc is done
                !     sending, remember how many messages they say they sent
                done_procs(source) = .true.
                nrecv_expected(source) = int(recved_xatomblocks(3,i_batch))
             else
                ! jd: Received an X atomblock, insert into X
                recved_nelements = nint(recved_xatomblocks(3,i_batch))
                first_ngwf_idx_of_a = &
                     ngwf_indices(A_NGWF)%basis%first_on_atom(recved_global_a)
                sparse_owner = &
                     ngwf_indices(A_NGWF)%basis%proc_of_func(first_ngwf_idx_of_a)
                call utils_assert(pub_my_proc_id == sparse_owner, myself//&
                     ': Received an atomblock that is not sparse-local to me',&
                     sparse_owner, pub_my_proc_id, first_ngwf_idx_of_a, &
                     recved_global_a)
                do is = 1, pub_num_spins
                   n_ngwfs_a = &
                        ngwf_indices(A_NGWF)%basis%num_on_atom(recved_global_a)
                   n_ngwfs_b = &
                        ngwf_indices(B_NGWF)%basis%num_on_atom(recved_global_b)
                   offset = 4 + (is-1) * n_ngwfs_a * n_ngwfs_b
                   call sparse_put_block(&
                        reshape(recved_xatomblocks(offset:,i_batch), &
                        (/n_ngwfs_b,n_ngwfs_a/)), xmatrix(is), &
                        recved_global_b, recved_global_a, accum = .true.)
                end do
                nrecv(source) = nrecv(source) + 1
             end if
          end if
       end do

       ! jd: Check the sends for completion
       ! --- send --- send --- send --- send --- send --- send --- send ---
       do i_batch = 1, send_batch_size
          if(sends_in_progress(i_batch)) then
             call comms_test(send_finished, send_handles(i_batch))
             if(send_finished) then
                sends_in_progress(i_batch) = .false.
             end if
          end if
       end do

    end do ! not done with the exchange yet

    ! jd: Everyone is done. Wait on remaining sends. Cancel unneeded receives.
    !     Wait until cancellations take effect.

    ! jd: Wait on sends
    do i_batch = 1, send_batch_size
       if(sends_in_progress(i_batch)) then
          call comms_wait(send_handles(i_batch))
       end if
    end do

    ! jd: Since we exited the loop, we must have already received all records.
    !     We only need to cancel any posted irecvs and wait for the cancels to
    !     finish before we can safely get rid of buffers.
    do i_batch = 1, recv_batch_size
       if(recvs_in_progress(i_batch)) then
          call comms_cancel(recv_handles(i_batch))
       end if
    end do

    ! jd: Wait for the cancellations to take effect
    do i_batch = 1, recv_batch_size
       if(recvs_in_progress(i_batch)) then
          call comms_wait(recv_handles(i_batch))
       end if
    end do

    ! jd: Wait on the sends for the done notifications
    do i_batch = 0, pub_total_num_procs-1
       call comms_wait(send_done_handles(i_batch))
    end do

    deallocate(recved_xatomblocks, stat=ierr)
    call utils_dealloc_check(myself,'recved_xatomblocks',ierr)

    call timer_clock(myself//'_imbalance',1)
    call comms_barrier
    call timer_clock(myself//'_imbalance',2)
#else
    ! jd: Non-MPI (commless) version
    eot = .false.
    do while(.not. eot)
       call hash_table_iterate(ndata, ptr, slot, x_ab_contributions_ht, &
            global_a, global_b, key3, key4, key5, eot, data_ptr, .false.)
       if(.not. eot) then
          ! jd: Have a record to store
          do is = 1, pub_num_spins
             n_ngwfs_a = ngwf_indices(A_NGWF)%basis%num_on_atom(global_a)
             n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)
             offset = 4 + (is-1) * n_ngwfs_a * n_ngwfs_b
             call sparse_put_block(reshape(data_ptr(offset:), &
                  (/n_ngwfs_b,n_ngwfs_a/)), xmatrix(is), &
                  global_b, global_a, accum = .true.)
          end do
       end if
    end do
#endif

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hfx_comms_xmatrix

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_comms_gradient(grad_akets_dkn_out_ht, &
       grad_akets_dkn_in_ht, cell_npts, alpha_ngwf_index)
    !==========================================================================!
    ! Redistributes NGWF gradient term kets from 'my B-C' distribution to the  !
    ! usual 'sparse-local A' distribution. On exit all ranks have the terms    !
    ! corresponding to their sparse-local atoms A.                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   grad_akets_dkn_out_ht (inout): Results are stored here. On entry it is !
    !                                  expected to have been initialised and   !
    !                                  empty.                                  !
    !   grad_akets_dkn_in_ht (in): Gradient dkn-contracted kets for my Bb are  !
    !                              taken from here.                            !
    !   cell_npts (in): cell%n_pts, needed to dimension the PPDs.              !
    !   alpha_ngwf_index (in): The HFX_NGWF_INDEX instance containing the      !
    !                          basis and the NGWF_REP used by A atoms.         !
    !                          Needed to know who sparse-owns which A's.       !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in April 2020.                                 !
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_any_source, pub_total_num_procs, &
         pub_null_handle, pub_on_root, comms_irecv, comms_send, comms_barrier, &
         comms_wait, comms_test, comms_cancel
    use constants, only: stdout, VERBOSE
    use hash_table, only: HT_STORAGE_NODE, hash_table_iterate, hash_table_axpy
    use rundat, only: pub_hfx_output_detail
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_assert, &
         utils_alloc_check, utils_dealloc_check, utils_flush

    implicit none

    ! jd: Arguments
    type(HT_HASH_TABLE), intent(inout), target :: grad_akets_dkn_out_ht
    type(HT_HASH_TABLE), intent(in), target    :: grad_akets_dkn_in_ht
    integer, intent(in)                        :: cell_npts
    type(HFX_NGWF_INDEX), intent(in)           :: alpha_ngwf_index

    ! jd: Local variables
    real(kind=DP), allocatable :: recved_records(:,:)
    integer :: global_aa_ngwf_idx
    integer :: recved_global_aa
    integer :: recved_cur_ppd
    integer :: recved_is
    integer :: recved_wk
    integer :: wk
    integer :: i_batch
    integer :: ndata
    type(HT_STORAGE_NODE), pointer :: ptr
    real(kind=DP), pointer :: data_ptr(:)
    integer :: slot
    integer :: key5
    logical :: eot
    logical :: recv_finished
    logical :: send_finished
    logical :: i_am_done_sending
    integer :: source
    integer :: msg_size
    integer :: sparse_owner
    integer :: i_rank
    integer :: cur_ppd
    integer :: is
    integer :: ierr
    real(kind=DP)      :: done(4,0:pub_total_num_procs-1)
    integer, parameter :: recv_batch_size = 1000
    integer, parameter :: send_batch_size = 1000
    integer            :: recv_handles(1:recv_batch_size)
    integer            :: send_handles(1:send_batch_size)
    integer            :: send_done_handles(0:pub_total_num_procs-1)
    integer            :: nsent(0:pub_total_num_procs-1)
    integer            :: nrecv(0:pub_total_num_procs-1)
    integer            :: nrecv_expected(0:pub_total_num_procs-1)
    logical            :: recvs_in_progress(1:recv_batch_size)
    logical            :: sends_in_progress(1:send_batch_size)
    logical            :: done_procs(0:pub_total_num_procs-1)
    character(len=*), parameter :: myself = 'hfx_comms_gradient'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    if(pub_hfx_output_detail >= VERBOSE .and. pub_on_root) then
       write(stdout,'(a)') 'HFx: Communicating NGWF gradient kets.'
       call utils_flush(stdout,.true.)
    end if

    call utils_assert(grad_akets_dkn_out_ht%n_slots + &
         grad_akets_dkn_out_ht%n_chained == 0, &
         myself//': Destination HT is meant to be empty on entry.', &
         grad_akets_dkn_out_ht%n_slots, &
         grad_akets_dkn_out_ht%n_chained)

    ! jd: Starts iteration from the 1st record of the HT
    ptr => NULL()
    slot = 0

#ifdef MPI
    msg_size = 4 + cell_npts
    allocate(recved_records(msg_size, recv_batch_size), stat=ierr)
    call utils_alloc_check(myself,'recved_records',ierr)

    recv_handles(:) = pub_null_handle
    send_handles(:) = pub_null_handle
    send_done_handles(:) = pub_null_handle
    i_am_done_sending = .false.
    done_procs(:) = .false.
    nsent(:) = 0
    nrecv(:) = 0
    nrecv_expected(:) = -1
    done(:,:) = -1
    recvs_in_progress(:) = .false.
    sends_in_progress(:) = .false.

    ! jd: Loop until everyone is done, ie. we get a 'done sending' notification
    !     from every sender and we received as many messages from every sender
    !     as they claim they sent.
    do while(.not. (all(done_procs(:)) .and. all(nrecv(:) == nrecv_expected(:))))

       ! jd: Post a batch of receives
       ! --- recv --- recv --- recv --- recv --- recv --- recv --- recv ---
       do i_batch = 1, recv_batch_size
          if(.not. recvs_in_progress(i_batch)) then
             call comms_irecv(pub_any_source, recved_records(:,i_batch), &
                  handle=recv_handles(i_batch))
             recvs_in_progress(i_batch) = .true.
          end if
       end do

       ! --- send --- send --- send --- send --- send --- send --- send ---
       if(.not. i_am_done_sending) then
          ! jd: Send a batch of records from our HT
          do i_batch = 1, send_batch_size
             if(sends_in_progress(i_batch)) cycle ! this send slot is being used
             ! jd: Look for next record to send
             call hash_table_iterate(ndata, ptr, slot, &
                  grad_akets_dkn_in_ht, global_aa_ngwf_idx, cur_ppd, &
                  is, wk, key5, eot, data_ptr, .false.)

             if(.not. eot) then
                ! jd: Have a record to send
                sparse_owner = &
                     alpha_ngwf_index%basis%proc_of_func(global_aa_ngwf_idx)
                call comms_send(sparse_owner, data_ptr, ndata, &
                     return_handle = send_handles(i_batch), &
                     add_to_stack = .false.)
                sends_in_progress(i_batch) = .true.
                nsent(sparse_owner) = nsent(sparse_owner) + 1
             else
                ! jd: I am done sending
                call utils_assert(sum(nsent(:)) == &
                     grad_akets_dkn_in_ht%n_slots + &
                     grad_akets_dkn_in_ht%n_chained, myself//&
                     ': Mismatch between number of grad_akets_dkn_in to be sent&
                     & and number of grad_akets_dkn_in_ppds_in actually sent.',&
                     sum(nsent(:)), grad_akets_dkn_in_ht%n_slots, &
                     grad_akets_dkn_in_ht%n_chained, pub_my_proc_id)
                ! jd: No more records, we iterated through entire HT.
                !     Send a 'done' message to everyone, tell them how many
                !     messages were sent.
                i_am_done_sending = .true.
                do i_rank = 0, pub_total_num_procs-1
                   done(3,i_rank) = nsent(i_rank)
                   call comms_send(i_rank, done(:,i_rank), 4, return_handle = &
                        send_done_handles(i_rank), add_to_stack = .false.)
                end do
                exit
             end if
          end do
       end if

       ! jd: Check the receives for completion
       ! --- recv --- recv --- recv --- recv --- recv --- recv --- recv ---
       do i_batch = 1, recv_batch_size
          if(.not. recvs_in_progress(i_batch)) cycle
          call comms_test(recv_finished, recv_handles(i_batch), source)
          if(recv_finished) then
             ! jd: Receive completed
             recved_global_aa = nint(recved_records(1,i_batch))
             recvs_in_progress(i_batch) = .false.
             if(recved_global_aa == -1) then
                ! jd: Received a done notification
                done_procs(source) = .true.
                nrecv_expected(source) = nint(recved_records(3,i_batch))
             else
                ! jd: Received a record, insert into dest HT
                recved_cur_ppd = nint(recved_records(2,i_batch))
                recved_is = nint(recved_records(3,i_batch))
                recved_wk = nint(recved_records(4,i_batch))
                sparse_owner = &
                     alpha_ngwf_index%basis%proc_of_func(recved_global_aa)
                call utils_assert(pub_my_proc_id == sparse_owner, myself//&
                     ': Received a gradient record that is not sparse-local &
                     &to me', sparse_owner, pub_my_proc_id, recved_global_aa)
                ! jd: Need to axpy, as multiple procs might have sent data from
                !     the same Aa (as [in index rename land] B-C pairs, not B's,
                !     are distributed).
                call hash_table_axpy(grad_akets_dkn_out_ht, &
                     recved_records(5:,i_batch), cell_npts, recved_global_aa, &
                     recved_cur_ppd, recved_is, recved_wk, &
                     overfill_strategy = 'F')
                nrecv(source) = nrecv(source) + 1
             end if
          end if
       end do

       ! jd: Check the sends for completion
       ! --- send --- send --- send --- send --- send --- send --- send ---
       do i_batch = 1, send_batch_size
          if(sends_in_progress(i_batch)) then
             call comms_test(send_finished, send_handles(i_batch))
             if(send_finished) then
                sends_in_progress(i_batch) = .false.
             end if
          end if
       end do

    end do ! not done with the exchange yet

    ! jd: Everyone is done. Wait on remaining sends. Cancel unneeded receives.
    !     Wait until cancellations take effect.

    ! jd: Wait on sends
    do i_batch = 1, send_batch_size
       if(sends_in_progress(i_batch)) then
          call comms_wait(send_handles(i_batch))
       end if
    end do

    ! jd: Since we exited the loop, we must have already received all records.
    !     We only need to cancel any posted irecvs and wait for the cancels to
    !     finish before we can safely get rid of buffers.
    do i_batch = 1, recv_batch_size
       if(recvs_in_progress(i_batch)) then
          call comms_cancel(recv_handles(i_batch))
       end if
    end do

    ! jd: Wait for the cancellations to take effect
    do i_batch = 1, recv_batch_size
       if(recvs_in_progress(i_batch)) then
          call comms_wait(recv_handles(i_batch))
       end if
    end do

    ! jd: Wait on the sends for the done notifications
    do i_batch = 0, pub_total_num_procs-1
       call comms_wait(send_done_handles(i_batch))
    end do

    deallocate(recved_records, stat=ierr)
    call utils_dealloc_check(myself,'recved_records',ierr)

    call timer_clock(myself//'_imbalance',1)
    call comms_barrier
    call timer_clock(myself//'_imbalance',2)
#else
    ! jd: Non-MPI (commless) version
    eot = .false.
    do while(.not. eot)
       call hash_table_iterate(ndata, ptr, slot, grad_akets_dkn_in_ht, &
            global_aa_ngwf_idx, cur_ppd, is, wk, key5, eot, data_ptr, .false.)
       if(.not. eot) then
          ! jd: Have a record to store
          call hash_table_axpy(grad_akets_dkn_out_ht, data_ptr(5:), &
               cell_npts, global_aa_ngwf_idx, cur_ppd, is, wk, &
               overfill_strategy = 'F')
       end if
    end do
#endif
    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hfx_comms_gradient

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_calc_aa_dd_ngwf_prods_in_ppds(hfxstate, ngwfs_hts, cell, &
       ngwf_indices, max_n_ppds_sphere)
    !==========================================================================!
    ! Calculates and stores in a HT all \phi_Aa \phi_Dd necessary on this proc !
    ! (using 'my' scheme, not 'sparse-local'). IOW, for all my A's it          !
    ! calculates the products of \phi_Aa with all \phi_Dd that s-overlap it.   !
    !                                                                          !
    ! However, given that these products can occupy a lot of RAM, particularly !
    ! in conduction calculations, we limit our appetite to what the user per-  !
    ! mitted via cache_limit_for_prods. In the end then, not *all* products are!
    ! stored, only the ones that fit. No attempt is made (yet) to judge if some!
    ! products are more valuable (more reusable). This is because the actual   !
    ! cost of re-doing these products over and over again is, currently, far   !
    ! from being a bottleneck.                                                 !
    !                                                                          !
    ! Any products that did not fit will be subsequently found absent in the   !
    ! HT and will be recalculated by hfx_calc_aa_dd_ngwf_prod_in_ppd().        !
    !                                                                          !
    ! This subroutine also stores \phi_Dd that will be needed for all my A's.  !
    ! These do not cost much and are all stored, there is no capping mechanism.!
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (inout): Results are stored in ppd_products_ht and in         !
    !                     ppd_ngwfs_dd_ht.                                     !
    !   ngwfs_hts (in): NGWF caches. NGWF data is taken from here.             !
    !   cell (in): Needed for n_pts.                                           !
    !   ngwf_indices (in): Array of HFX_NGWF_INDEX instances containing        !
    !                      pointers to the NGWF basis and the NGWF_REP used    !
    !                      by atoms A, B, C and D (only used for A and D).     !
    !   max_n_ppds_sphere (in): Bigger of the two. Used to dimension arrays.   !
    !--------------------------------------------------------------------------!
    ! Caveat:                                                                  !
    !   If A and D do not sparse-overlap, no entry will be added for this pair.!
    !   Be careful later not to look for A-D pairs that share a PPD, but       !
    !   do not actually overlap.                                               !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2019.                                 !
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_total_num_procs, pub_on_root, &
         comms_allgather
    use constants, only: garbage_real, VERBOSE, stdout
    use hash_table, only: hash_table_lookup_ptr_nocount, hash_table_add, &
         hash_table_list, hash_table_mem_usage
    use ppd_ops, only: ppd_product
    use rundat, only: pub_threads_max, pub_cache_limit_for_prods, &
         pub_hfx_output_detail
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_assert, &
         utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target  :: hfxstate
    type(HT_HASH_TABLE), intent(in), target :: ngwfs_hts(:)
    type(CELL_INFO), intent(in)             :: cell
    type(HFX_NGWF_INDEX), intent(in)        :: ngwf_indices(A_NGWF:D_NGWF)
    integer, intent(in)                     :: max_n_ppds_sphere

    ! jd: Local variables
    integer :: ppd_indices_in_product_of_aa_dd(1:max_n_ppds_sphere)
    integer :: ppd_indices_in_ngwf_dd(1:max_n_ppds_sphere)
    real(kind=DP) :: ppd_product_of_aa_dd(max_n_ppds_sphere * cell%n_pts)
    real(kind=DP) :: ppd_ngwf_dd(max_n_ppds_sphere * cell%n_pts)
    real(kind=DP), allocatable :: ppd_product_ad(:,:,:)
    real(kind=DP), allocatable :: ppd_ngwfs_d(:,:)
    integer :: n_ppds_in_product_of_aa_dd
    integer :: n_ppds_in_ngwf_dd
    integer :: my_a_atom_idx
    integer :: my_d_atom_idx
    integer :: d_idx
    integer :: global_a
    integer :: global_d
    real(kind=DP), pointer :: packed_ngwf_aa(:) ! ptr to ht's internals
    real(kind=DP), pointer :: packed_ngwf_dd(:) ! ptr to ht's internals
    integer :: packed_ngwf_aa_size, packed_ngwf_dd_size
    integer :: first_ngwf_idx_of_a, first_ngwf_idx_of_d
    integer, allocatable :: ad_pairs(:,:)
    integer :: ngwf_a, ngwf_d
    integer :: n_ngwfs_a, n_ngwfs_d
    integer :: global_aa_ngwf_idx
    integer :: global_dd_ngwf_idx
    integer :: ndata_product
    integer :: ndata_ngwf_dd
    integer :: max_n_ppds_points
    integer :: npoints
    integer :: i_ppd
    integer :: cur_ppd
    integer :: ppd_start
    integer :: ppd_end
    integer :: i_pair
    integer :: n_pairs
    real(kind=DP) :: mem_used
    logical :: exit_now
    integer :: n_ad_pairs_added
    integer :: n_pairs_added_all(0:pub_total_num_procs-1)
    integer :: n_pairs_all(0:pub_total_num_procs-1)
    integer :: proc
    real(kind=DP) :: frac_hit
    integer :: ierr
    character(len=*), parameter :: myself = 'hfx_calc_aa_dd_ngwf_prods_in_ppds'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)
    call timer_clock(myself,1)
    ! jd: Note, can't have barriers in here now that this is called a different
    !     number of times from each proc.

    if(.not. hfxstate%ppd_products_ht_in_use) goto 100

    ! *************************************************************************
    ! jd: 1. Generate a list of A-D pairs -- a fused loop will OMP-parallelise
    !        better
    ! *************************************************************************

    n_pairs = 0
    ! --------------------------------------------------------------------------
    ! jd: Loop over A's whose S-neighbours D we're looking for            MY AAA
    ! --------------------------------------------------------------------------
    loop_A0:                                                                   &
    do my_a_atom_idx = 1, hfxstate%n_my_a_atoms
       global_a = hfxstate%my_a_atoms(my_a_atom_idx)

       ! -----------------------------------------------------------------------
       ! jd: Loop over D's that are S-neighbours of A                     MY DDD
       ! -----------------------------------------------------------------------
       loop_D0:                                                                &
       do d_idx = hfxstate%s_atoms_nl%first_idx(global_a), &
            hfxstate%s_atoms_nl%last_idx(global_a)
          global_d = hfxstate%s_atoms_nl%neighbours(d_idx)

          n_pairs = n_pairs + 1

       end do loop_D0

    end do loop_A0

    allocate(ad_pairs(2,n_pairs),stat=ierr)
    call utils_alloc_check(myself,'ad_pairs',ierr)

    i_pair = 1
    ! --------------------------------------------------------------------------
    ! jd: Loop over A's whose S-neighbours D we're looking for            MY AAA
    ! --------------------------------------------------------------------------
    loop_A1:                                                                   &
    do my_a_atom_idx = 1, hfxstate%n_my_a_atoms
       global_a = hfxstate%my_a_atoms(my_a_atom_idx)

       ! -----------------------------------------------------------------------
       ! jd: Loop over D's that are S-neighbours of A                     MY DDD
       ! -----------------------------------------------------------------------
       loop_D1:                                                                &
       do d_idx = hfxstate%s_atoms_nl%first_idx(global_a), &
            hfxstate%s_atoms_nl%last_idx(global_a)
          global_d = hfxstate%s_atoms_nl%neighbours(d_idx)

          ad_pairs(1,i_pair) = global_a
          ad_pairs(2,i_pair) = global_d
          i_pair = i_pair + 1

       end do loop_D1

    end do loop_A1

    call utils_assert(i_pair-1 == n_pairs,myself//': Internal error.')

    ! --------------------------------------------------------------------------
    ! jd: Fused loop over A-D pairs.                                     AAA-DDD
    ! --------------------------------------------------------------------------
    exit_now = .false.
    n_ad_pairs_added = 0
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(DYNAMIC) &
!$OMP PRIVATE(i_pair, first_ngwf_idx_of_a, first_ngwf_idx_of_d, global_a,      &
!$OMP      global_d, n_ngwfs_a, n_ngwfs_d, ppd_product_ad, ierr, ndata_product,&
!$OMP      ngwf_d, global_dd_ngwf_idx, packed_ngwf_dd, packed_ngwf_dd_size,    &
!$OMP      ngwf_a, global_aa_ngwf_idx, packed_ngwf_aa, packed_ngwf_aa_size,    &
!$OMP      ppd_indices_in_product_of_aa_dd, ppd_product_of_aa_dd, npoints,     &
!$OMP      i_ppd, n_ppds_in_product_of_aa_dd, cur_ppd, ppd_start, ppd_end,     &
!$OMP      mem_used)                                                           &
!$OMP SHARED(n_pairs, max_n_ppds_points, cell, ngwfs_hts, hfxstate, ad_pairs,  &
!$OMP      ngwf_indices, max_n_ppds_sphere, pub_cache_limit_for_prods,         &
!$OMP      exit_now, pub_my_proc_id)                                           &
!$OMP REDUCTION(+:n_ad_pairs_added)
    loop_AD:                                                                   &
    do i_pair = 1, n_pairs

       if(exit_now) cycle ! jd: Cannot 'exit' from a parallel DO

       global_a = ad_pairs(1,i_pair)
       global_d = ad_pairs(2,i_pair)

       first_ngwf_idx_of_a = ngwf_indices(A_NGWF)%basis%first_on_atom(global_a)
       n_ngwfs_a = ngwf_indices(A_NGWF)%basis%num_on_atom(global_a)

       first_ngwf_idx_of_d = ngwf_indices(D_NGWF)%basis%first_on_atom(global_d)
       n_ngwfs_d = ngwf_indices(D_NGWF)%basis%num_on_atom(global_d)

       max_n_ppds_points = max_n_ppds_sphere * cell%n_pts
       allocate(ppd_product_ad(max_n_ppds_points,n_ngwfs_a,n_ngwfs_d),stat=ierr)
       call utils_alloc_check(myself,'ppd_product_ad',ierr)

       ! -----------------------------------------------------------------------
       ! jd: for all d on D                                                  ddd
       ! -----------------------------------------------------------------------
       ndata_product = 0
       loop_ngwf_d:                                                            &
       do ngwf_d = 1, n_ngwfs_d
          global_dd_ngwf_idx = first_ngwf_idx_of_d + ngwf_d - 1

          ! jd: Get NGWF Dd
          call hash_table_lookup_ptr_nocount(packed_ngwf_dd, &
               packed_ngwf_dd_size, &
               ngwfs_hts(ngwf_indices(D_NGWF)%rep%ngwf_cache_handle), &
               global_dd_ngwf_idx)
          if(packed_ngwf_dd_size == -1) then
             call utils_abort(myself//': NGWF Dd absent in cache [1].', &
                  global_dd_ngwf_idx, global_d, ngwf_d, global_a, &
                  avoid_mpi_calls = .true.)
          end if

          ! --------------------------------------------------------------------
          ! jd: for all a on A                                               aaa
          ! --------------------------------------------------------------------
          loop_ngwf_a:                                                         &
          do ngwf_a = 1, n_ngwfs_a
             global_aa_ngwf_idx = first_ngwf_idx_of_a + ngwf_a - 1

             ! jd: Get NGWF Aa
             call hash_table_lookup_ptr_nocount(packed_ngwf_aa, &
                  packed_ngwf_aa_size, &
                  ngwfs_hts(ngwf_indices(A_NGWF)%rep%ngwf_cache_handle), &
                  global_aa_ngwf_idx)
             if(packed_ngwf_aa_size == -1) then
                call utils_abort(myself//': NGWF Aa absent in cache [1].', &
                     global_aa_ngwf_idx, global_a, ngwf_a, &
                     avoid_mpi_calls = .true.)
             end if

             ! jd: Calculate PPD product
             call ppd_product(packed_ngwf_aa, packed_ngwf_dd, cell%n_pts, &
                  ppd_product_of_aa_dd, ppd_indices_in_product_of_aa_dd, &
                  n_ppds_in_product_of_aa_dd)
             npoints = n_ppds_in_product_of_aa_dd * cell%n_pts

             ppd_product_ad(1:npoints,ngwf_a,ngwf_d) = &
                  ppd_product_of_aa_dd(1:npoints)
             ppd_product_ad(npoints+1:,ngwf_a,ngwf_d) = garbage_real
             ! jd: Safety net, might want to remove     ^

             ! jd: # of reals in {all a, all d, one PPD}.
             ndata_product = ndata_product + cell%n_pts

          end do loop_ngwf_a

       end do loop_ngwf_d

       ! jd: Split this into PPDs and add each one to the HT
       ! jd: Note-1: Moving the critical into the i_ppd loop is slower.
       !     Note-2: Using stashes instead of the critical is slower too.
       !             Stashes seem to help only when there are lookups concurrent
       !             to stores -- because they permit eliding the critical in
       !             the lookup.
       ! @possible_optimisation:
       !              Do not use a HT for storing PPD products. Use a giant
       !              array and manually do the book-keeping for indices.
       !              Then we will be able to populate this without a critical
       !              here. Also lookups and especially purges will be faster.
       !              However, since this subroutine has negligible cost, it's
       !              probably not worth the effort. Maybe for conduction?
!$OMP CRITICAL(SEC_PPD_PRODUCTS_HT_WRITE)
       do i_ppd = 1, n_ppds_in_product_of_aa_dd

          if(exit_now) exit

          cur_ppd = ppd_indices_in_product_of_aa_dd(i_ppd)
          ppd_start = 1+(i_ppd-1) * cell%n_pts
          ppd_end = ppd_start + cell%n_pts - 1
          call hash_table_add(hfxstate%ppd_products_ht, &
               reshape(ppd_product_ad(ppd_start:ppd_end,:,:),&
               (/ndata_product/)), ndata_product, global_a, global_d, &
               cur_ppd, overfill_strategy = 'F')

          call hash_table_mem_usage(hfxstate%ppd_products_ht, mem_used)

          if(mem_used >= pub_cache_limit_for_prods) then
             exit_now = .true.
             exit
          end if

       end do
!$OMP END CRITICAL(SEC_PPD_PRODUCTS_HT_WRITE)

       deallocate(ppd_product_ad,stat=ierr)
       call utils_dealloc_check(myself,'ppd_product_ad',ierr)

       n_ad_pairs_added = n_ad_pairs_added + 1

    end do loop_AD
!$OMP END PARALLEL DO

    if(pub_hfx_output_detail >= VERBOSE) then

       call comms_allgather(n_pairs_added_all, n_ad_pairs_added, 1, &
            gather_not_allgather = .true.)
       call comms_allgather(n_pairs_all, n_pairs, 1, &
            gather_not_allgather = .true.)

       if(pub_on_root) then
          write(stdout,'(a)') '+'//repeat('-',78)//'+'
          write(stdout,'(a)') '|    MPI |                          AD NGWF prod&
               &uct cache                      |'
          write(stdout,'(a)') '|   rank | products stored | products left &
               &out |           total | store ratio |'
          write(stdout,'(a)') '+'//repeat('-',78)//'+'
          do proc = 0, pub_total_num_procs-1
             n_ad_pairs_added = n_pairs_added_all(proc)
             n_pairs = n_pairs_all(proc)
             if(n_pairs > 0) then
                frac_hit = real(n_ad_pairs_added,kind=DP) / &
                     real(n_pairs,kind=DP) * 100.0_DP
             else
                frac_hit = 1.0_DP
             end if
             write(stdout,'(a,i7,a,i16,a,i18,a,i16,a,f11.2,a)') '|', proc, &
                  ' |', n_ad_pairs_added, ' |', n_pairs - n_ad_pairs_added, &
                  ' |', n_pairs, ' |', frac_hit ,'% |'
          end do
          write(stdout,'(a)') '+'//repeat('-',78)//'+'
       end if

    end if

    deallocate(ad_pairs,stat=ierr)
    call utils_dealloc_check(myself,'ad_pairs',ierr)

    call hash_table_list(hfxstate%ppd_products_ht,0)

100 continue

    if(pub_hfx_output_detail >= VERBOSE .and. pub_on_root .and. &
         .not. hfxstate%ppd_products_ht_in_use) then
       write(stdout,'(a)') '+'//repeat('-',78)//'+'
       write(stdout,'(a)') '| AD NGWF products are not cached (cache_limit_for_&
            &prods == 0).                |'
       write(stdout,'(a)') '+'//repeat('-',78)//'+'
    end if

    ! *************************************************************************
    ! --- D NGWFs ---
    ! *************************************************************************

    ! --------------------------------------------------------------------------
    ! jd: Loop over my D's                                                MY DDD
    ! --------------------------------------------------------------------------
    loop_D2:                                                                   &
    do my_d_atom_idx = 1, hfxstate%n_my_d_atoms

       global_d = hfxstate%my_d_atoms(my_d_atom_idx)

       first_ngwf_idx_of_d = ngwf_indices(D_NGWF)%basis%first_on_atom(global_d)
       n_ngwfs_d = ngwf_indices(D_NGWF)%basis%num_on_atom(global_d)

       max_n_ppds_points = max_n_ppds_sphere * cell%n_pts
       allocate(ppd_ngwfs_d(max_n_ppds_points,n_ngwfs_d),stat=ierr)
       call utils_alloc_check(myself,'ppd_ngwfs_d',ierr)

       ! -----------------------------------------------------------------------
       ! jd: for all d on D                                                  ddd
       ! -----------------------------------------------------------------------
       ndata_ngwf_dd = 0
       loop_ngwf_d2:                                                           &
       do ngwf_d = 1, n_ngwfs_d
          global_dd_ngwf_idx = first_ngwf_idx_of_d + ngwf_d - 1

          ! jd: Get NGWF Dd
          call hash_table_lookup_ptr_nocount(packed_ngwf_dd, &
               packed_ngwf_dd_size, &
               ngwfs_hts(ngwf_indices(D_NGWF)%rep%ngwf_cache_handle), &
               global_dd_ngwf_idx)
          call utils_assert(packed_ngwf_dd_size /= -1, &
               myself//': NGWF Dd absent in cache [2].', &
               global_dd_ngwf_idx, global_d, ngwf_d)

          ! jd: Calculate NGWF Dd, used later in gradient
          call ppd_product(packed_ngwf_dd, packed_ngwf_dd, cell%n_pts, &
               ppd_ngwf_dd, ppd_indices_in_ngwf_dd, &
               n_ppds_in_ngwf_dd, use_1_for_bra = .true.)
          npoints = n_ppds_in_ngwf_dd * cell%n_pts

          ppd_ngwfs_d(1:npoints,ngwf_d) = ppd_ngwf_dd(1:npoints)
          ppd_ngwfs_d(npoints+1:,ngwf_d) = garbage_real
          ! jd: Safety net, might want to remove ^

          ! jd: # of reals in {all d, one PPD}.
          ndata_ngwf_dd = ndata_ngwf_dd + cell%n_pts

       end do loop_ngwf_d2

       ! jd: Split this into PPDs and add each one to the HT
       do i_ppd = 1, n_ppds_in_ngwf_dd
          cur_ppd = ppd_indices_in_ngwf_dd(i_ppd)
          ppd_start = 1+(i_ppd-1) * cell%n_pts
          ppd_end = ppd_start + cell%n_pts - 1
          call hash_table_add(hfxstate%ppd_ngwfs_dd_ht, &
               reshape(ppd_ngwfs_d(ppd_start:ppd_end,:),(/ndata_ngwf_dd/)), &
               ndata_ngwf_dd, global_d, cur_ppd, overfill_strategy = 'F')
       end do

       deallocate(ppd_ngwfs_d,stat=ierr)
       call utils_dealloc_check(myself,'ppd_ngwfs_d',ierr)

    end do loop_D2

    call hash_table_list(hfxstate%ppd_ngwfs_dd_ht,0)

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hfx_calc_aa_dd_ngwf_prods_in_ppds

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_calc_aa_dd_ngwf_prod_in_ppd(ppd_product_ad, &
       ngwfs_hts, cell, ngwf_indices, global_a, global_d, cur_ppd, &
       sample, sample_period)
    !==========================================================================!
    ! Calculates a product of all pairs of NGWFs centred on two atoms in a PPD.!
    ! The NGWFs are taken from NGWF HTs. Results are stored in ppd_product_ad. !
    !                                                                          !
    ! This subroutine is called billions of times for larger systems and as    !
    ! such needs to be reasonably efficient. Hence we only call timer every    !
    ! once in a while to reduce the overhead.                                  !
    !                                                                          !
    ! There is still plenty of room for optimisation -- for instance we look up!
    ! both NGWFs in the HT every time we are called (so, once for every PPD    !
    ! they share). We could do better by fetching entire NGWFs once, but this  !
    ! is made complicated, as elsewhere (in the caller) it is important to     !
    ! iterate *by PPDs* and not by NGWFs. So, optimisation is possible, but has!
    ! not been attempted yet.                                                  !
    !                                                                          !
    ! This subroutine needs to be (and is) threadsafe. It is called from an    !
    ! OMP loop over PPDs.                                                      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   ppd_product_ad: Results go here. Caller is responsible for passing an  !
    !                   array that is max_on_atom * max_on_atom * n_pts big.   !
    !                   This is sufficiently small to fit on the (OMP) stack.  !
    !   ngwfs_hts (in): NGWF caches. NGWF data is taken from here.             !
    !   cell (in): Needed for n_pts.                                           !
    !   ngwf_indices (in): Array of HFX_NGWF_INDEX instances containing        !
    !                      pointers to the NGWF basis and the NGWF_REP used    !
    !                      by atoms A, B, C and D (only used for A and D).     !
    !   global_a (in): The global index of the 1st NGWF in the product.        !
    !   global_d (in): The global index of the 2nd NGWF in the product.        !
    !   cur_ppd (in): The index of the PPD.                                    !
    !   sample, sample_period (in): Timer control. See caller.                 !
    !--------------------------------------------------------------------------!
    ! Caveat:                                                                  !
    !   Behaviour when passing a PPD that does not belong to the intersection  !
    !   of A and D (or when A and D do not overlap) is undefined.              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2020.                              !
    !==========================================================================!

    use constants, only: garbage_real
    use hash_table, only: hash_table_lookup_ptr_nocount
    use ppd_ops, only: ppd_product_one_ppd
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    type(HFX_NGWF_INDEX), intent(in)        :: ngwf_indices(A_NGWF:D_NGWF)
    type(CELL_INFO), intent(in)             :: cell
    real(kind=DP), intent(out)              :: ppd_product_ad(cell%n_pts* &
         ngwf_indices(A_NGWF)%basis%max_on_atom * &
         ngwf_indices(D_NGWF)%basis%max_on_atom)
    integer, intent(in)                     :: global_a
    integer, intent(in)                     :: global_d
    type(HT_HASH_TABLE), intent(in), target :: ngwfs_hts(:)
    integer, intent(in)                     :: cur_ppd
    integer, intent(in)                     :: sample
    integer, intent(in)                     :: sample_period

    ! jd: Local variables
    real(kind=DP), pointer :: packed_ngwf_aa(:) ! ptr to ht's internals
    real(kind=DP), pointer :: packed_ngwf_dd(:) ! ptr to ht's internals
    integer :: packed_ngwf_aa_size, packed_ngwf_dd_size
    integer :: first_ngwf_idx_of_a, first_ngwf_idx_of_d
    integer :: ngwf_a, ngwf_d
    integer :: n_ngwfs_a, n_ngwfs_d
    integer :: global_aa_ngwf_idx
    integer :: global_dd_ngwf_idx
    integer :: cache_offset_a, cache_offset_d
    integer :: dest
    character(len=*), parameter :: myself = 'hfx_calc_aa_dd_ngwf_prod_in_ppd'

    ! -------------------------------------------------------------------------

    ! jd: Do not use utils_trace_{in,out} -- we are called from OMP regions.
    ! jd: Do not use MPI barriers here -- we are called a different number of
    !     times on different procs.

    if(mod(sample,sample_period) == 0) then
       call timer_clock(myself,1,multiplier = sample_period)
    end if

    first_ngwf_idx_of_a = ngwf_indices(A_NGWF)%basis%first_on_atom(global_a)
    n_ngwfs_a = ngwf_indices(A_NGWF)%basis%num_on_atom(global_a)

    first_ngwf_idx_of_d = ngwf_indices(D_NGWF)%basis%first_on_atom(global_d)
    n_ngwfs_d = ngwf_indices(D_NGWF)%basis%num_on_atom(global_d)

    dest = 1

    ! -------------------------------------------------------------------------
    ! jd: for all d on D                                                    ddd
    ! -------------------------------------------------------------------------
    cache_offset_d = -1
    loop_ngwf_d:                                                               &
    do ngwf_d = 1, n_ngwfs_d
       global_dd_ngwf_idx = first_ngwf_idx_of_d + ngwf_d - 1

       ! jd: Get NGWF Dd
       call hash_table_lookup_ptr_nocount(packed_ngwf_dd, &
            packed_ngwf_dd_size, &
            ngwfs_hts(ngwf_indices(D_NGWF)%rep%ngwf_cache_handle), &
            global_dd_ngwf_idx)
       if(packed_ngwf_dd_size == -1) then
          call utils_abort(myself//': NGWF Dd absent in cache [1].', &
               global_dd_ngwf_idx, global_d, ngwf_d, global_a, &
               avoid_mpi_calls = .true.)
       end if

       ! ----------------------------------------------------------------------
       ! jd: for all a on A                                                 aaa
       ! ----------------------------------------------------------------------
       cache_offset_a = -1
       loop_ngwf_a:                                                            &
       do ngwf_a = 1, n_ngwfs_a
          global_aa_ngwf_idx = first_ngwf_idx_of_a + ngwf_a - 1

          ! jd: Get NGWF Aa
          call hash_table_lookup_ptr_nocount(packed_ngwf_aa, &
               packed_ngwf_aa_size, &
               ngwfs_hts(ngwf_indices(A_NGWF)%rep%ngwf_cache_handle), &
               global_aa_ngwf_idx)
          if(packed_ngwf_aa_size == -1) then
             call utils_abort(myself//': NGWF Aa absent in cache [1].', &
                  global_aa_ngwf_idx, global_a, ngwf_a, &
                  avoid_mpi_calls = .true.)
          end if

          ! jd: Calculate PPD product.
          call ppd_product_one_ppd(packed_ngwf_aa, packed_ngwf_dd, &
               cell%n_pts, cur_ppd, ppd_product_ad(dest:dest+cell%n_pts-1), &
               cache_offset_a, cache_offset_d)
          !    ^ cache_offsets allow finding the offset to a and d in the
          !      packed NGWFs in O(1) time. Otherwise we'd have to scan the
          !      packed structure for this PPD multiple times.

          dest = dest + cell%n_pts

       end do loop_ngwf_a

    end do loop_ngwf_d

    ppd_product_ad(dest:) = garbage_real

    if(mod(sample,sample_period) == 0) then
       call timer_clock(myself,2, multiplier = sample_period)
    end if

  end subroutine hfx_calc_aa_dd_ngwf_prod_in_ppd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_dkn_indep_preparation(hfxstate, &
       swri, ireg, par, mdl, &
       rep, ngwf_basis, rep2, ngwf_basis2, basis_selector)
    !==========================================================================!
    ! Does all the necessary preparation of HFx datastructures that are        !
    ! density kernel independent.                                              !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (inout): Container for the HFx datastructures.                !
    !   swri (in/out): SW_RI container, holds metric matrices. In/out because  !
    !                  its SWOP cache is populated here.                       !
    !   ireg (in): Embedding region identifier.                                !
    !   par (in): Needed to come up with parallel distribution of B-C pairs.   !
    !   mdl (in): Needed for cell, elements.                                   !
    !   rep (in):        } Describe the NGWF basis.                            !
    !   ngwf_basis (in): }                                                     !
    !                                                                          !
    ! Optional arguments for mixed NGWF bases. Either pass none, or all.       !
    !   See banner for hf_exchange_calculate() for meaning and convention.     !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2019, cannibalising                !
    ! hf_exchange_calculate() by Quintin Hill (2008-2009), rewritten by JD to  !
    ! use 2-centre 'BC' expansion (2012), split by JD to separate NGWF exp'n   !
    ! from calculation of exchange (2016), extended by JD (2018) to support    !
    ! mixed NGWF bases.                                                        !
    !==========================================================================!

    use constants, only: stdout, SW_O, SW_V, REP_SWEX_HFX
    use comms, only: pub_on_root
    use constants, only: VERBOSE
    use hash_table, only: hash_table_purge
    use model_type, only: MODEL
    use remote, only: remote_ngwf_cache_hts
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_hfx_output_detail, pub_cond_calculate
    use sparse, only: sparse_destroy
    use sw_resolution_of_identity, only: SW_RI
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_abort, utils_postfix_to_ngwf_set_name,&
         utils_trace_in, utils_trace_out, utils_flush

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target         :: hfxstate
    type(SW_RI), intent(inout)                     :: swri
    integer, intent(in)                            :: ireg  ! embed. region
    type(PARAL_INFO), intent(in)                   :: par
    type(MODEL), intent(in)                        :: mdl
    type(NGWF_REP), intent(in), target             :: rep
    type(FUNC_BASIS), intent(in), target           :: ngwf_basis
    ! jd: Arguments for mixed NGWF bases
    type(NGWF_REP), optional, intent(in), target   :: rep2
    type(FUNC_BASIS), optional, intent(in), target :: ngwf_basis2
    integer, optional, intent(in)                  :: basis_selector(5)

    ! jd: Local variables
    type(HFX_NGWF_INDEX) :: ngwf_indices(A_NGWF:D_NGWF)
    logical :: loc_cmplx ! agrecocmplx
    logical :: destroying_full_metric
    integer :: region_nat
    integer :: swex_handle
    character(len=*), parameter :: myself = 'hfx_dkn_indep_preparation'

    ! -----------------------------------------------------------------------

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    if(pub_on_root .and. pub_hfx_output_detail >= VERBOSE) then
       if(.not. present(rep2)) then
          write(stdout,'(a,a,a)') 'HFx: Stage-1 for ', &
               trim(utils_postfix_to_ngwf_set_name(rep%postfix)), &
               ' (single NGWF basis).'
       else
          write(stdout,'(a,a,a,a,a,i0,i0,i0,i0,i0,a)') 'HFx: Stage-1 for ', &
               trim(utils_postfix_to_ngwf_set_name(rep%postfix)), '+', &
               trim(utils_postfix_to_ngwf_set_name(rep2%postfix)), &
               ', selector: ', &
               basis_selector, ' (mixed NGWF basis).'
       end if
    end if

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(ireg)%iscmplx

    ! agrecocmplx
    call utils_assert(.not. loc_cmplx, &
         myself//': not ready yet for complex NGWFs.')

    call utils_assert(par%nat == size(mdl%regions(ireg)%elements), myself//&
         ': Allocated parallel strategy is incompatible with elements.', &
         par%nat, size(mdl%regions(ireg)%elements))

    region_nat = par%nat

    ! jd: Ensure consistency between optional parameters for mixed NGWF bases
    call utils_assert((present(rep2) .eqv. present(ngwf_basis2)) .and. &
         (present(rep2) .eqv. present(basis_selector)), &
         myself//': Optional arguments for HFx with mixed NGWF bases must &
         &either all be present or all be absent.')

    call utils_assert((present(rep2) .eqv. present(ngwf_basis2)) .and. &
         (present(rep2) .eqv. present(basis_selector)), &
         myself//': Optional arguments for HFx with mixed NGWF bases must &
         &either all be present or all be absent.')

    ! jd: Select the right source NGWF_REP (whose SW_EX holds the coeffs),
    !     and right NGWF bases.
    if(present(basis_selector)) then
       !jmecmplx
       call utils_assert(.not.rep2%ngwfs_on_grid(ireg)%iscmplx, myself//&
            ': not ready for complex NGWFs yet (rep2).')
       call utils_assert(rep%ngwfs_on_grid(ireg)%iscmplx .eqv. &
            rep2%ngwfs_on_grid(ireg)%iscmplx, myself//&
            ': Cannot mix real and complex NGWFs between REPs.')
       !call utils_assert(basis_selector(1) == basis_selector(2), myself//&
       !     ': Mixing bases *within* a density kernel (AB) not supported yet.')
       if (.not.basis_selector(1) == basis_selector(2).and.pub_on_root) then
          write(stdout,'(a)') "HFx WARNING: mixed basis density kernel (AB) --&
               & experimental functionality"
       end if
       !call utils_assert(basis_selector(3) == basis_selector(4), myself//&
       !     ': Mixing bases *within* a density kernel (CD) not supported yet.')
       if (.not.basis_selector(3) == basis_selector(4).and.pub_on_root) then
          write(stdout,'(a)') "HFx WARNING: mixed basis density kernel (CD) --&
               & experimental functionality"
       end if
       if(basis_selector(1) == 1) then
          ngwf_indices(A_NGWF)%basis => ngwf_basis
          ngwf_indices(A_NGWF)%rep => rep
       else if(basis_selector(1) == 2) then
          ngwf_indices(A_NGWF)%basis => ngwf_basis2
          ngwf_indices(A_NGWF)%rep => rep2
       else
          call utils_abort(myself//': Illegal basis_selector(1)')
       end if
       if(basis_selector(2) == 1) then
          ngwf_indices(B_NGWF)%basis => ngwf_basis
          ngwf_indices(B_NGWF)%rep => rep
       else if(basis_selector(2) == 2) then
          ngwf_indices(B_NGWF)%basis => ngwf_basis2
          ngwf_indices(B_NGWF)%rep => rep2
       else
          call utils_abort(myself//': Illegal basis_selector(2)')
       end if
       if(basis_selector(3) == 1) then
          ngwf_indices(C_NGWF)%basis => ngwf_basis
          ngwf_indices(C_NGWF)%rep => rep
       else if(basis_selector(3) == 2) then
          ngwf_indices(C_NGWF)%basis => ngwf_basis2
          ngwf_indices(C_NGWF)%rep => rep2
       else
          call utils_abort(myself//': Illegal basis_selector(3)')
       end if
       if(basis_selector(4) == 1) then
          ngwf_indices(D_NGWF)%basis => ngwf_basis
          ngwf_indices(D_NGWF)%rep => rep
       else if(basis_selector(4) == 2) then
          ngwf_indices(D_NGWF)%basis => ngwf_basis2
          ngwf_indices(D_NGWF)%rep => rep2
       else
          call utils_abort(myself//': Illegal basis_selector(4)')
       end if
    else
       ngwf_indices(A_NGWF)%basis => ngwf_basis
       ngwf_indices(B_NGWF)%basis => ngwf_basis
       ngwf_indices(C_NGWF)%basis => ngwf_basis
       ngwf_indices(D_NGWF)%basis => ngwf_basis
       ngwf_indices(A_NGWF)%rep => rep
       ngwf_indices(B_NGWF)%rep => rep
       ngwf_indices(C_NGWF)%rep => rep
       ngwf_indices(D_NGWF)%rep => rep
    end if

    ! jd: If arrays in hfxstate are allocated (just check one), we've been
    !     called already (in previous NGWF steps). Need to free previous arrays
    !     before we start allocating new ones.
    if(allocated(hfxstate%my_bc_pairs)) then
       call hfx_free_arrays(hfxstate)
    end if

    ! ========================================================================
    ! This part is independent of NGWFs. Could be optimised to only be done
    ! when atoms move.
    ! ========================================================================
    if(pub_on_root .and. pub_hfx_output_detail >= VERBOSE) then
       write(stdout,'(a)') 'HFx: - Parallel distribution... '
       call utils_flush(stdout, .true.)
    end if

    ! jd: Determine B-C pairs to process. The array gets allocated here.
    call hfx_determine_bc_pairs(hfxstate, par, ngwf_indices)

    ! jd: Determine B and C atoms to process. The arrays get allocated here.
    call hfx_split_pair_list(hfxstate, region_nat, ngwf_indices)

    ! jd: Determine A atoms to process. The array gets allocated here.
    !     A's are X-neighbours with my B.
    call hfx_determine_my_a_atoms(hfxstate, region_nat)

    ! jd: Determine D atoms to process. The array gets allocated here.
    !     D's are s-X-neighbours with my B.
    call hfx_determine_my_d_atoms(hfxstate, region_nat)

    if(pub_on_root .and. pub_hfx_output_detail >= VERBOSE) then
       write(stdout,'(a)') 'HFx: - Communicating metric matrix (ket: BC)...'
       call utils_flush(stdout, .true.)
    end if

    ! jd: Communicate metric matrices (turn sparse-local order into my-order).
    !     This is for kets(BbCc). If possible, destroy the original SPAM3
    !     structures for full metric matrices -- they are no longer needed
    !     now that we switched to my-BC distribution. However, there are
    !     exceptions:
    !     - if the swri is shared with DMA. Since DMA uses the old sparse-local
    !       distribution, if we share the swri, we don't destroy the matrices
    !       so that swx_expand_local_ngwf_pairs() keeps working,
    !     - if doing conduction. This is because B-C pairs will be rebalanced
    !       across procs (because # of NGWFs changes when moving from val to
    !       cond), and my-owners of matrix blocks will change. This will require
    !       re-communicating the metric. We'd either have to recalculate the V
    !       matrix or to keep a copy. We choose the latter.
    destroying_full_metric = &
         .not. hfxstate%swri_is_shared_with_dma .and. .not. pub_cond_calculate

    if(swri%has_metric(SW_V)) then
       call hfx_comms_metric(hfxstate, swri%full_metric_matrices(SW_V), SW_V)
       if(destroying_full_metric) then
          call sparse_destroy(swri%full_metric_matrices(SW_V))
       end if
    end if
    if(swri%has_metric(SW_O)) then
       call hfx_comms_metric(hfxstate, swri%full_metric_matrices(SW_O), SW_O)
       if(destroying_full_metric) then
          call sparse_destroy(swri%full_metric_matrices(SW_O))
       end if
    end if

    ! jd: Make sure SW_RI's destructor doesn't try to destroy again if we
    !     did destroy it here.
    swri%full_metric_matrices_destroyed = .not. destroying_full_metric

    ! ========================================================================

    if(pub_on_root .and. pub_hfx_output_detail >= VERBOSE) then
       write(stdout,'(a)') 'HFx: - Communicating NGWFs...'
       call utils_flush(stdout, .true.)
    end if

    ! jd: *** Communicate necessary NGWFs. ***
    !     - A NGWFs are needed for \phi_Aa \phi_Dd products.
    !     - B NGWFs are needed for NGWF expansion (of \phi_Bb \phi_Cc).
    !     - C NGWFs are needed for NGWF expansion (of \phi_Bb \phi_Cc).
    !     - D NGWFs are needed for \phi_Aa \phi_Dd products.
    ! jd: Communicate NGWFs of my A atoms. Store them in HTs.
    call hfx_comms_ngwfs(remote_ngwf_cache_hts, &      ! A's are not necessarily
         hfxstate%my_a_atoms, hfxstate%n_my_a_atoms, & ! sparse-local, hence [*]
         ngwf_indices(A_NGWF), ireg)                   ! involvement of remote

    ! jd: Communicate NGWFs of my B atoms. Store them in HTs.
    call hfx_comms_ngwfs(remote_ngwf_cache_hts, &      ! B's are not necessarily
         hfxstate%my_b_atoms, hfxstate%n_my_b_atoms, & ! sparse-local, hence [*]
         ngwf_indices(B_NGWF), ireg)                   ! involvement of remote

    ! jd: Communicate NGWFs of my C atoms. Store them in HTs.
    call hfx_comms_ngwfs(remote_ngwf_cache_hts, &      ! C's are not necessarily
         hfxstate%my_c_atoms, hfxstate%n_my_c_atoms, & ! sparse-local, hence [*]
         ngwf_indices(C_NGWF), ireg)                   ! involvement of remote

    ! jd: Communicate NGWFs of my D atoms. Store them in HTs.
    call hfx_comms_ngwfs(remote_ngwf_cache_hts, &      ! D's are not necessarily
         hfxstate%my_d_atoms, hfxstate%n_my_d_atoms, & ! sparse-local, hence [*]
         ngwf_indices(D_NGWF), ireg)                   ! involvement of remote

    if(pub_on_root .and. pub_hfx_output_detail >= VERBOSE) then
       write(stdout,'(a)') 'HFx: - Index book-keeping...'
       call utils_flush(stdout, .true.)
    end if

    ! jd: For all my B determine lists of D's for every A PPD,
    !     where A is x-neigbhour of B
    call hfx_determine_dlists(hfxstate, remote_ngwf_cache_hts, &
         ngwf_indices)

    ! jd: Determine all PPDs where we need the expansions (PPDs of A's).
    !     The array gets allocated here.
    call hfx_determine_ppd_list(hfxstate, remote_ngwf_cache_hts, & ! same [*]
         ngwf_indices(A_NGWF))                                     ! here

    ! jd: Determine lists of A's for every my A PPD,
    call hfx_determine_alists(hfxstate, remote_ngwf_cache_hts, &
         ngwf_indices(A_NGWF))

    ! jd: Populate SWOP caches with the SWOPs deemed most intensely used.
    if(.not. hfxstate%swop_cache_populated) then

       ! jd: In conduction things get a little tricky:
       !     - init() gets called for val,
       !     - we get here for val calc and populate SWOP caches,
       !     - init() with init_rep_only gets called again from conduction,
       !       setting swop_cache_populated to .false., which makes sense as
       !       they will have to be recalculated, since the change in n_ngwfs
       !       affects HFx distribution,
       !     - we get here again for cond calc. We should populate the caches,
       !       but they contain the old SWOPs. Purge them.
       !     In valence the HT is empty, so it won't matter.
       call hash_table_purge(swri%swops_in_ppds_ht)

       ! jd: How SWOP hits are counted depends on Bb-Cc symmetry, so establish
       !     that first. The handle is also needed for establishing quality.
       if(present(basis_selector)) then
          swex_handle = basis_selector(5)
       else
          swex_handle = REP_SWEX_HFX
       end if

       if(pub_on_root .and. pub_hfx_output_detail >= VERBOSE) then
          write(stdout,'(a)') 'HFx: - Estimating memory use...'
          call utils_flush(stdout, .true.)
       end if

       call hfx_populate_swop_cache(hfxstate, swri, par, mdl, ireg, &
            ngwf_indices, swex_handle)

    end if

    ! jd: Calculate \phi_Aa * \phi_Dd in PPDs. Due to memory constraints, it's
    !     likely not all will fit
    if(pub_on_root .and. pub_hfx_output_detail >= VERBOSE) then
       write(stdout,'(a)') 'HFx: - Calculating NGWF products...'
       call utils_flush(stdout, .true.)
    end if
    call hfx_calc_aa_dd_ngwf_prods_in_ppds(hfxstate, remote_ngwf_cache_hts, &
         mdl%cell, ngwf_indices, &
         max(ngwf_indices(A_NGWF)%basis%max_n_ppds_sphere, &
         ngwf_indices(D_NGWF)%basis%max_n_ppds_sphere))
    ! =========================================================================

    if(pub_on_root .and. pub_hfx_output_detail >= VERBOSE) then
       write(stdout,'(a)') 'HFx: ...done.'
       call utils_flush(stdout, .true.)
    end if

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hfx_dkn_indep_preparation

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_expand_ngwf_pairs(hfxstate, &
       rep, ngwf_basis, mdl, ireg, &
       rep2, ngwf_basis2, basis_selector) ! opt args for mixed NGWF bases
    !==========================================================================!
    ! Expands products of NGWFs in spherical waves. This is a wrapper for      !
    ! swx_expand_ngwf_pairs(). Suitable initialisation is ensured first.       !
    ! Expansions are invalidated (as NGWFs must have changed).                 !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (inout): Container for HFx state. Needed for SWRI handle and  !
    !                     X neighbour list.                                    !
    !   rep (in/out): The NGWF rep whose NGWF pairs will be expanded. Also,    !
    !                 the results of the expansion are stored here (in the     !
    !                 'swexes' member). In a single NGWF basis scenario, the   !
    !                 coeffs are stored at REP_SWEX_HFX. With mixed NGWF bases,!
    !                 basis_selector(5) decides.                               !
    !   ngwf_basis (in): The NGWF basis corresponding to rep.                  !
    !   mdl (in): Needed for cell, elements, regions.                          !
    !   ireg (in): Embedding region ID.                                        !
    !                                                                          !
    ! *** Optional arguments for mixed NGWF bases only.                        !
    !   rep2 (in): Second rep to use when mixing NGWF bases. Can be intent-in, !
    !              because results (coeffs) are always stored in the first rep.!
    !   ngwf_basis2 (in): Second NGWF basis.                                   !
    !   basis_selector (in): Chooses the NGWF bases for the calculation of     !
    !                      TEFCIs. Quintin's alpha-beta-gamma-delta convention !
    !                      convention is used (cf. [1], eq. 5). Alpha is in the!
    !                      bra, beta is in the ket. These two are the indices  !
    !                      of the exchange matrix. Delta is in the bra,        !
    !                      gamma  is in the ket. K^{\delta,\gamma}, or whatever!
    !                      is passed for denskern_cd is contracted over when   !
    !                      building the X matrix. This is different from Tim's !
    !                      convention in JCP 2013, and different from the CCP9 !
    !                      Flagship meeting convention.                        !
    !                      basis_selector(1) chooses basis for alpha (1 or 2). !
    !                      basis_selector(2) chooses basis for beta (1 or 2).  !
    !                      basis_selector(3) chooses basis for gamma (1 or 2). !
    !                      basis_selector(4) chooses basis for delta (1 or 2). !
    !                      basis_selector(5) chooses the *source* swex in rep, !
    !                                        from which the expansion coeffs   !
    !                                        will be retrieved.                !
    !                      Here only the ket part is referenced.               !
    !                      (/1,1,1,1,REP_SWEX_HFX/) is for valence.            !
    !                      (/1,1,2,2,REP_SWEX_HFX_OTHER/) is for conduction,   !
    !                      where rep = cond_rep, rep2 = val_rep, so cond is 1, !
    !                      val is 2.                                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October 2016.                               !
    ! Extended by Jacek Dziedzic in June-August 2018 for mixed NGWF basis sets.!
    ! Modified by Robert Charlton to allow embedding in a subsystem, 13/09/18. !
    !==========================================================================!

    use comms, only: comms_barrier
    use constants, only: METRIC_ELECTROSTATIC, METRIC_OVERLAP, &
         SW_V, SW_O, REP_SWEX_HFX
    use model_type, only: MODEL
    use rundat, only: pub_hfx_metric, pub_active_region, pub_emft, &
         pub_use_activehfx, pub_emft_follow
    use sw_expansion, only: swx_expand_my_ngwf_pairs
    use sw_resolution_of_identity, only: swri_library
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_abort, utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target         :: hfxstate
    type(NGWF_REP), intent(inout), target          :: rep
    type(FUNC_BASIS), intent(in), target           :: ngwf_basis
    type(MODEL), intent(in)                        :: mdl
    integer, intent(in)                            :: ireg
    ! jd: Arguments for mixed NGWF bases
    type(NGWF_REP), optional, intent(in), target   :: rep2
    type(FUNC_BASIS), optional, intent(in), target :: ngwf_basis2
    integer, optional, intent(in)                  :: basis_selector(5)

    ! jd: --- Local variables ---
    integer :: coeffs_kind
    integer :: ket_basis_selector(3)
    ! agrecocmplx
    logical :: loc_cmplx
    character(len=*), parameter :: myself = 'hfx_expand_ngwf_pairs'

    ! -----------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    if(present(basis_selector)) then
       ket_basis_selector(1) = basis_selector(2)
       ket_basis_selector(2) = basis_selector(3)
       ket_basis_selector(3) = basis_selector(5)
    else
       ket_basis_selector(:) = 1
    end if

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(ireg)%iscmplx

    ! agrecocmplx
    call utils_assert(.not. loc_cmplx, &
         myself//': not ready yet for complex NGWFs.')

    ! jd: Ensure we're initialised
    call utils_assert(hfxstate%hfx_swri_h /= -1, &
         myself//': Need to call hf_exchange_init() first')

    ! rc2013: If we're doing embedding check that the SWRI covers all
    ! atoms in the active region
    if((pub_emft .or. pub_emft_follow) .and. pub_use_activehfx) then
       if(count(mdl%regions(pub_active_region)%elements(:)%in_swri(hfxstate%hfx_swri_h))&
            /= mdl%regions(pub_active_region)%par%nat) then
          call utils_abort('Hartree-Fock exchange cannot be used on just a &
               &subsystem of the active region. Ensure that your SWRI ['//&
               trim(swri_library(hfxstate%hfx_swri_h)%swri_name)//'] includes all &
               &species listed in the first line of %BLOCK species_ngwf_regions.')
       end if
    else
       ! jd: Ensure the HFx SWRI covers the entire system
       if(count(mdl%elements(:)%in_swri(hfxstate%hfx_swri_h)) /= mdl%nat) then
           call utils_abort('Your SWRI does not include all atoms in the input &
               &file. If you want to do a subsystem Hartree-Fock exchange &
               &calculation include %BLOCK species_ngwf_regions. &
               &Otherwise ensure that your SWRI ['//&
               trim(swri_library(hfxstate%hfx_swri_h)%swri_name)//&
               '] includes all atoms.')
       end if
    end if

    ! jd: Sanity check on X-neighbour list
    call utils_assert(hfxstate%x_atoms_nl%populated, &
         myself//': x_atoms neighbour list not populated.')

    ! jd: Ensure consistency between optional parameters for mixed NGWF bases
    call utils_assert((present(rep2) .eqv. present(ngwf_basis2)) .and. &
         (present(rep2) .eqv. present(basis_selector)), &
         myself//': Optional arguments for HFx with mixed NGWF bases must &
         &either all be present or all be absent.')

    ! =========================================================================

    if(pub_hfx_metric == METRIC_ELECTROSTATIC) coeffs_kind = SW_V
    if(pub_hfx_metric == METRIC_OVERLAP) coeffs_kind = SW_O

    ! jd: Perform SW expansion for all the products of NGWFs, where B is mine
    !     (in my_b_atoms). This is distinct from rank-local B's.
    !     B's S-neighbours C are also not in general rank-local. No comms are
    !     needed. The results are stored in hfx_swex%coeffs_ht.
    if(.not. present(rep2)) then
       ! jd: --- DEFAULT CASE, NO MIXED BASIS ---
       call swx_expand_my_ngwf_pairs(rep, swri_library(hfxstate%hfx_swri_h), & ! in/out
            hfxstate%packed_metric_matrix_blocks_hts, &                 ! in
            hfxstate%my_bc_pairs, hfxstate%n_my_bc_pairs, &             ! in
            REP_SWEX_HFX, ngwf_basis, mdl%cell, mdl%regions(ireg)%par, &! input
            mdl%regions(ireg)%elements, coeffs_kind)                    ! input
    else
       ! jd: --- MIXED BASIS CASE ---
       call swx_expand_my_ngwf_pairs(rep, swri_library(hfxstate%hfx_swri_h), & ! in/out
            hfxstate%packed_metric_matrix_blocks_hts, &                ! in
            hfxstate%my_bc_pairs, hfxstate%n_my_bc_pairs, &            ! in
            ket_basis_selector(3), &                                   ! in
            ngwf_basis, mdl%cell, mdl%regions(ireg)%par, &             ! input
            mdl%regions(ireg)%elements, coeffs_kind, &                 ! in
            rep2=rep2, ngwf_basis2=ngwf_basis2, &                      ! in opt
            ket_basis_selector=ket_basis_selector)                     ! in opt
    end if

    call comms_barrier
    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hfx_expand_ngwf_pairs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_calc_expansions_bb_cc(expansions_ht, & ! in/out (result)
       table_full, &                                                 ! out
       bb_cc_symmetry, swri, swex_quality, cell, ppd_list, n_ppds, & ! in
       expansion_atoms, expansion_centres, coeffs, &                 ! in
       global_bb_ngwf_idx, global_cc_ngwf_idx, mode)                 ! in
    !==========================================================================!
    ! Generates a potential coming from an SW expansion of an NGWF product.    !
    ! The potential is generated in PPDs and stored in expansions_ht.          !
    ! This operation is costly, and the aim is to store all needed expansions  !
    ! in the HT and to reuse them later.                                       !
    !                                                                          !
    ! SWOPs in PPDs are cached behind the scenes. NGWF pairs must have been    !
    ! expanded in advance (coeffs must be available).                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   expansions_ht (inout): Results go here.                                !
    !   table_full (out): Will be set to .T. if expansions HT became full.     !
    !   bb_cc_symmetry (in): If .true. Bb-Cc symmetry will be exploited.       !
    !                        When mixed NGWF bases are used, this symmetry is  !
    !                        no longer present.                                !
    !   swri (inout): SW_RI object in which we do the resolution of identiy.   !
    !                 Needed for the SWOP cache.                               !
    !   swex_quality (in): Determines the quality of the expansion.            !
    !   cell (in): Needed for generating SWOPs.                                !
    !   ppd_list (in): The list of PPDs for which we want the expansions.      !
    !   n_ppds (in): Number of elements in ppd_list.                           !
    !   expansion_atoms (in): Atoms in the expansion.                          !
    !   expansion_centres (in): Centres in the expansion.                      !
    !   coeffs (in): Expansion coefficients.                                   !
    !   global_bb_ngwf_idx (in): } Global NGWF indices for the source pair of  !
    !   global_cc_ngwf_idx (in): } the expansion.                              !
    !   mode (in): If mode == 1, expansions where SWOPs are cached for either  !
    !              centre (or the centre, if there is one) are not calculated. !
    !              This helps populate the HT with expansions that will be more!
    !              costly to obtain later (by calling this routine with mode 1 !
    !              first, and then mode 2).                                    !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Note that this subroutine assumes all PPDs belong to the same sphere   !
    !   (cf. array dimensions).                                                !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in August 2012, as hfx_calculate_expanded_pot. !
    ! Extended by Jacek Dizedzic in June 2018 to handle mixed NGWF basis sets. !
    ! Adapted by Jacek Dziedzic in March 2019 to the new parallel scheme.      !
    !==========================================================================!

    use basis, only: basis_ppd_location
    use constants, only: CRLF
    use geometry, only: POINT
    use hash_table, only: HT_HASH_TABLE, hash_table_probe, hash_table_add, &
         hash_table_probe_nocount
    use rundat, only: pub_threads_max, pub_hfx_debug
    use simulation_cell, only: CELL_INFO, minimum_image_number_to_displacement,&
         minimum_image_number_to_indices
    use sw_resolution_of_identity, only: SW_RI, ATOM_CENTRE, SW_QUALITY, &
         PBC_IMAGE_INFO, swri_obtain_swops_in_ppd_ptr_fast
    use utils, only: utils_int_to_str, utils_point_to_str, utils_real_to_str, &
         utils_abort
    use timer, only: timer_clock

    implicit none

    ! jd: Arguments
    type(HT_HASH_TABLE), intent(inout), target :: expansions_ht
    logical, intent(out)          :: table_full
    logical, intent(in)           :: bb_cc_symmetry
    type(SW_RI), intent(in)       :: swri
    type(SW_QUALITY), intent(in)  :: swex_quality
    type(CELL_INFO), intent(in)   :: cell
    integer, intent(in)           :: expansion_atoms(2)
    type(ATOM_CENTRE), intent(in) :: expansion_centres(2)
    real(kind=DP), intent(in)     :: coeffs(2*swex_quality%max_sws_per_centre)
    integer, intent(in)           :: n_ppds
    integer, intent(in)           :: ppd_list(:) ! 1:n_ppds
    integer, intent(in)           :: global_bb_ngwf_idx
    integer, intent(in)           :: global_cc_ngwf_idx
    integer, intent(in)           :: mode

    ! jd: Local variables
    integer           :: cur_ppd, cur_ppd_n
    integer           :: result
    logical           :: expansion_is_cached
    integer           :: size_of_cached_data
    integer           :: sw1_idx
    integer           :: e
    integer           :: global_ngwf_idx_1, global_ngwf_idx_2
    integer           :: centre_offset
    integer           :: src_atom
    integer           :: num_sws_on_src, num_sws
!$  integer, external :: omp_get_thread_num
    integer           :: tid
    integer           :: ndata1, ndata2
    integer           :: ipt
    real(kind=DP)     :: coeff
    type(ATOM_CENTRE) :: src_centre
    type(POINT)       :: ppd_corner
    type(POINT)       :: im1_disp, im2_disp
    integer           :: im1_a1_neighbour, im2_a1_neighbour
    integer           :: im1_a2_neighbour, im2_a2_neighbour
    integer           :: im1_a3_neighbour, im2_a3_neighbour
    real(kind=DP), pointer :: swops_in_ppd_ptr(:)
    real(kind=DP), target  :: swops_in_ppd_workspace(&
         cell%n_pts * swri%quality%num_sws_per_centre, 0:pub_threads_max-1)
    real(kind=DP)     :: expansion_in_ppd(cell%n_pts)
    type(PBC_IMAGE_INFO) :: which_images(cell%n_pts,2)
    character(len=*), parameter :: myself = 'hfx_calc_expansions_bb_cc'

    ! --------------------------------------------------------------------------

    call timer_clock(myself,1)

    if(.not. bb_cc_symmetry) then
       global_ngwf_idx_1 = global_bb_ngwf_idx
       global_ngwf_idx_2 = global_cc_ngwf_idx
    else
       global_ngwf_idx_1 = min(global_bb_ngwf_idx,global_cc_ngwf_idx)
       global_ngwf_idx_2 = max(global_bb_ngwf_idx,global_cc_ngwf_idx)
    end if

    ! --------------------------------------------------------------------------
    ! Go over relevant PPDs                                                  PPD
    ! --------------------------------------------------------------------------
    ! Optimisation note. I initially thought that it would be better to use
    ! swri_obtain_swops_in_ppd_set_omp() to reduce contention on HT accesses
    ! here, and to store expansions for entire PPD sets in one go (dropping the
    ! PPD as key). The latter won't work because of Bb/Cc symmetry. The set of
    ! PPDs is different for Bb-Cc than for Cc-Bb, so this breaks down as soon
    ! as X is truncated. Also, it wasn't actually faster. Finally, it seems that
    ! the CRITICALs here do not cost much, less than 1% of the execution time.
    ! Perhaps using a stash would help a little.
    !
    ! The entire loop takes:
    !  - ~10-40ms (8-20 threads) when calculating the expansions.
    !  - 0.4-2.0ms (8-20 threads) when they are retrieved from cache.
    !                             ^ this we don't do anymore, we just don't
    !                               call this subroutine when we have expansions
    ! Type of SCHEDULE does not make practically any difference.
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(DYNAMIC) &
!$OMP PRIVATE(cur_ppd_n, cur_ppd, expansion_is_cached, result, tid,            &
!$OMP      size_of_cached_data, e, src_atom, src_centre, num_sws, ndata1,      &
!$OMP      num_sws_on_src, centre_offset, coeff, expansion_in_ppd, ndata2,     &
!$OMP      swops_in_ppd_ptr, which_images, ppd_corner, im1_disp, im2_disp,     &
!$OMP      im1_a1_neighbour, im1_a2_neighbour, im1_a3_neighbour,               &
!$OMP      im2_a1_neighbour, im2_a2_neighbour, im2_a3_neighbour)               &
!$OMP SHARED(n_ppds, ppd_list, cell, expansions_ht, global_ngwf_idx_1, swri,   &
!$OMP      global_ngwf_idx_2, expansion_atoms, expansion_centres, swex_quality,&
!$OMP      coeffs, swops_in_ppd_workspace, mode, pub_hfx_debug)
    all_ppds:                                                                  &
    do cur_ppd_n = 1, n_ppds

       cur_ppd = ppd_list(cur_ppd_n)
       tid = 0
!$     tid = omp_get_thread_num()

       ! Check if the entire expansion of BbCc into SWs is cached for this PPD.
       ! If so, nothing to do.
!$OMP CRITICAL(SEC_EXPANSIONS_HT_ACCESS)
       result = hash_table_probe(expansions_ht, &
            cur_ppd, global_ngwf_idx_1, global_ngwf_idx_2)
!$OMP END CRITICAL(SEC_EXPANSIONS_HT_ACCESS)
       if(result /= -1) cycle

       ! ************************************************
       ! jd: Expansion not cached, need to calculate it.
       ! ************************************************

       ! jd: But first check if we should skip this expansion in this turn.
       !     If mode == 1 we only do the ones that will prove the most costly
       !     to recalculate, ie. where both SWOPs are not in the cache. IOW,
       !     we skip the ones where either or both SWOPs are in the cache.
       if(mode == 1) then
          ndata1 = hash_table_probe_nocount(swri%swops_in_ppds_ht, &
               key1 = expansion_atoms(1), key2 = cur_ppd, key3 = ichar('P'))

          if(ndata1 > 0) cycle ! swop 1 in cache, cycle

          if(expansion_atoms(2) /= -1) then
             ndata2 = hash_table_probe_nocount(swri%swops_in_ppds_ht, &
                  key1 = expansion_atoms(2), key2 = cur_ppd, key3 = ichar('P'))
          else
             ndata2 = ndata1
          end if

          if(ndata2 > 0) cycle ! swop 2 in cache, cycle

       end if

       !     If mode == 2 we do every expansion left

       expansion_in_ppd(:) = 0.0_DP

       ! jd: Loop over unique centres in the expansion
       two_centres:                                                         &
       do e = 1, 2
          src_atom = expansion_atoms(e)
          src_centre = expansion_centres(e)
          if(src_atom == -1) exit

          num_sws_on_src = swex_quality%num_sws_per_centre

          ! jd: Get these SWOPs from hash table or calculate, if absent.
          !     This assumes the HT has been populated (although possibly
          !     not exhaustively). It doesn't update the HT. No locks or
          !     criticals are required.
          if(.not. pub_hfx_debug) then
             call swri_obtain_swops_in_ppd_ptr_fast(swops_in_ppd_ptr, num_sws, & ! out
                  src_centre, src_atom, cur_ppd, 'P', cell, swri, &              ! in
                  swops_in_ppd_workspace(:,tid)) !^ we always use SWpots once OVO
          else                                   !  metric got dropped
             call swri_obtain_swops_in_ppd_ptr_fast(swops_in_ppd_ptr, num_sws, & ! out
                  src_centre, src_atom, cur_ppd, 'P', cell, swri, &              ! in
                  swops_in_ppd_workspace(:,tid), which_images(:,e))
          end if

          if(e == 1) centre_offset = 0
          if(e == 2) centre_offset = num_sws_on_src
          ! --------------------------------------------------------------------
          ! jd: for all SWs on src atom                                       SW
          ! --------------------------------------------------------------------
          sws_on_one:                                                          &
          do sw1_idx = 1, num_sws_on_src

             coeff = coeffs(sw1_idx + centre_offset)
!            call utils_sanity_check((/coeff/),'single_coeff',excessive=1D12)

             expansion_in_ppd(:) = expansion_in_ppd(:) + &
                  coeff * swops_in_ppd_ptr(&
                  (sw1_idx-1)*cell%n_pts+1:sw1_idx*cell%n_pts)

!            jd: ^ Equivalent BLAS call, practically identical performance
!            call linalg_daxpy(cell%n_pts, coeff, swops_in_ppd_ptr(&
!                 (sw1_idx-1)*cell%n_pts+1:sw1_idx*cell%n_pts), 1, expansion_in_ppd, 1)

          end do sws_on_one

       end do two_centres

       ! jd: Debug two-centre expansions
       if(pub_hfx_debug .and. expansion_atoms(2) /= -1) then
          do ipt = 1, cell%n_pts
             if(which_images(ipt,1)%image /= which_images(ipt,2)%image) then

                ppd_corner = basis_ppd_location(cur_ppd,cell)

                im1_disp = minimum_image_number_to_displacement(&
                     which_images(ipt,1)%image, cell)
                im2_disp = minimum_image_number_to_displacement(&
                     which_images(ipt,2)%image, cell)

                call minimum_image_number_to_indices(which_images(ipt,1)%image, &
                     im1_a1_neighbour, im1_a2_neighbour, im1_a3_neighbour)
                call minimum_image_number_to_indices(which_images(ipt,2)%image, &
                     im2_a1_neighbour, im2_a2_neighbour, im2_a3_neighbour)

#if 1
             write(*,'(a,i3,a,i2,i2,i2,a,i3,a,i2,i2,i2,a,f14.6,f14.6,f14.6,a,&
                  &f14.6,f14.6,f14.6,a)') '[2] Image 1: '//&
                  trim(utils_int_to_str(which_images(ipt,1)%image))//&
                  '['//trim(utils_int_to_str(im1_a1_neighbour))//&
                  ' '//trim(utils_int_to_str(im1_a2_neighbour))//&
                  ' '//trim(utils_int_to_str(im1_a3_neighbour))//&
                  '], Image 2: '//&
                  trim(utils_int_to_str(which_images(ipt,2)%image))//&
                  '['//trim(utils_int_to_str(im2_a1_neighbour))//&
                  ' '//trim(utils_int_to_str(im2_a2_neighbour))//&
                  ' '//trim(utils_int_to_str(im2_a3_neighbour))//'].'//&
                  ' Centre 1 at '//trim(utils_point_to_str(&
                  expansion_centres(1)%incell,'f14.6'))//'. '//&
                  'Centre 2 at '//trim(utils_point_to_str(&
                  expansion_centres(2)%incell,'f14.6'))//'.'
#endif

#if 1
                call utils_abort(myself//&
                     ': [2] Inconsistent images for point '//&
                     trim(utils_int_to_str(ipt))//' in PPD '//&
                     trim(utils_int_to_str(cur_ppd))//'.'//&
                     CRLF//'Image 1: '//&
                     trim(utils_int_to_str(which_images(ipt,1)%image))//&
                     '['//trim(utils_int_to_str(im1_a1_neighbour))//&
                     ' '//trim(utils_int_to_str(im1_a2_neighbour))//&
                     ' '//trim(utils_int_to_str(im1_a3_neighbour))//&
                     '], Image 2: '//&
                     trim(utils_int_to_str(which_images(ipt,2)%image))//&
                     '['//trim(utils_int_to_str(im2_a1_neighbour))//&
                     ' '//trim(utils_int_to_str(im2_a2_neighbour))//&
                     ' '//trim(utils_int_to_str(im2_a3_neighbour))//'].'//&
                     CRLF//'Centre 1 at '//trim(utils_point_to_str(&
                     expansion_centres(1)%incell,'f14.6'))//'.'//CRLF//&
                     'Centre 2 at '//trim(utils_point_to_str(&
                     expansion_centres(2)%incell,'f14.6'))//'.'//CRLF//&
                     'PPD starts at '//trim(utils_point_to_str(&
                     ppd_corner,'f14.6'))//'.'//CRLF//&
                     'Point at '//trim(utils_point_to_str(&
                     which_images(ipt,1)%curpoint,'f14.6'))//'.'//CRLF//&
                     'Direct displacement to centre 1: '//trim(utils_point_to_str(&
                     which_images(ipt,1)%disp0,'f14.6'))//'.'//CRLF//&
                     'Direct displacement to centre 2: '//trim(utils_point_to_str(&
                     which_images(ipt,2)%disp0,'f14.6'))//'.'//CRLF//&
                     'Displacement to image 1: '//trim(utils_point_to_str(&
                     which_images(ipt,1)%disp,'f14.6'))//'.'//CRLF//&
                     'Displacement to image 2: '//trim(utils_point_to_str(&
                     which_images(ipt,2)%disp,'f14.6'))//'.'//CRLF//&
                     'Direct distance to centre 1: '//trim(utils_real_to_str(&
                     which_images(ipt,1)%rad0))//'.'//CRLF//&
                     'Direct distance to centre 2: '//trim(utils_real_to_str(&
                     which_images(ipt,2)%rad0))//'.'//CRLF//&
                     'Distance to image 1: '//trim(utils_real_to_str(&
                     which_images(ipt,1)%rad))//'.'//CRLF//&
                     'Distance to image 2: '//trim(utils_real_to_str(&
                     which_images(ipt,2)%rad))//'.', avoid_mpi_calls = .true.)
#endif

             end if
          end do
       end if

       ! jd: Add the expansion for this PPD to cache.
!$OMP CRITICAL(SEC_EXPANSIONS_HT_ACCESS)
       call hash_table_add(expansions_ht, expansion_in_ppd, cell%n_pts, &
            cur_ppd, global_ngwf_idx_1, global_ngwf_idx_2, &
            overfill_strategy = 'N')
!$OMP END CRITICAL(SEC_EXPANSIONS_HT_ACCESS)

    end do all_ppds
!$OMP END PARALLEL DO

    ! jd: If table full, remember it
    table_full = (expansions_ht%n_slots + expansions_ht%n_chained >= &
         expansions_ht%max_elements)

    call timer_clock(myself,2)

  end subroutine hfx_calc_expansions_bb_cc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_calculate_single_expansion(expansion_in_ppd, &
       expansion_atoms, expansion_centres, cell, swri, coeffs, cur_ppd)
    !==========================================================================!
    ! Generates a potential coming from an SW expansion of an NGWF product.    !
    ! The potential is generated in a PPD and stored in expansion_in_ppd.      !
    !                                                                          !
    ! This subroutine is only used to calculate expansions on the fly, when    !
    ! some did not fit in the expansion cache (expansions_ht). It is then      !
    ! called from an OMP region.                                               !
    !                                                                          !
    ! This almost duplicates some of the code in hfx_calc_expansions_bb_cc(),  !
    ! but this is for reasons of efficiency.                                   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   expansion_in_ppd (out): The expansion is returned here. The buffer     !
    !                           must be exactly cell%n_pts elements long.      !
    !   expansion_atoms (in): Atoms in the expansion.                          !
    !   expansion_centres (in): Centres in the expansion.                      !
    !   cell (in): The simulation cell.                                        !
    !   swri (in): Needed for SWOPs, quality.                                  !
    !   coeffs (in): Expansion coefficients.                                   !
    !   cur_ppd (in): The index of the PPD in which to expand.                 !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Keep in mind that this subroutine is called from an OMP region.        !
    !   Do not do any comms here (no barriers!), do not output (utils_trace).  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October 2019.                               !
    !==========================================================================!

    use basis, only: basis_ppd_location
    use constants, only: CRLF
    use geometry, only: POINT
    use rundat, only: pub_hfx_debug
    use simulation_cell, only: CELL_INFO, minimum_image_number_to_displacement,&
         minimum_image_number_to_indices
    use sw_resolution_of_identity, only: SW_RI, ATOM_CENTRE, PBC_IMAGE_INFO, &
         swri_obtain_swops_in_ppd_ptr_fast
    use utils, only: utils_int_to_str, utils_point_to_str, utils_real_to_str, &
         utils_abort

    implicit none

    ! jd: Arguments
    type(CELL_INFO), intent(in)   :: cell
    real(kind=DP), intent(out)    :: expansion_in_ppd(cell%n_pts)
    integer, intent(in)           :: expansion_atoms(2)
    type(ATOM_CENTRE), intent(in) :: expansion_centres(2)
    type(SW_RI), intent(in)       :: swri
    real(kind=DP), intent(in)     :: coeffs(:)
    integer, intent(in)           :: cur_ppd

    ! jd: Local variables
    type(ATOM_CENTRE)      :: src_centre
    real(kind=DP), pointer :: swops_in_ppd_ptr(:)
    real(kind=DP), target  :: swops_in_ppd_workspace(&
         cell%n_pts * swri%quality%num_sws_per_centre)
    real(kind=DP)          :: coeff
    integer                :: sw1_idx
    integer                :: e
    integer                :: src_atom
    integer                :: centre_offset
    integer                :: num_sws
    type(POINT)            :: ppd_corner
    type(POINT)            :: im1_disp, im2_disp
    integer                :: ipt
    integer                :: im1_a1_neighbour, im2_a1_neighbour
    integer                :: im1_a2_neighbour, im2_a2_neighbour
    integer                :: im1_a3_neighbour, im2_a3_neighbour
    type(PBC_IMAGE_INFO)   :: which_images(cell%n_pts,2)
    character(len=*), parameter :: myself = 'hfx_calculate_single_expansion'

    ! -------------------------------------------------------------------------

    expansion_in_ppd(:) = 0.0_DP

    ! jd: Loop over unique centres in the expansion
    two_centres:                                                               &
    do e = 1, 2
       src_atom = expansion_atoms(e)
       src_centre = expansion_centres(e)
       if(src_atom == -1) exit

       ! jd: Get these SWOPs from hash table or calculate, if absent.
       !     This assumes the HT has been populated (although possibly
       !     not exhaustively). It doesn't update the HT. No locks or
       !     criticals are required.
       if(.not. pub_hfx_debug) then
          call swri_obtain_swops_in_ppd_ptr_fast(swops_in_ppd_ptr, num_sws, & ! out
               src_centre, src_atom, cur_ppd, 'P', cell, swri, &              ! in
               swops_in_ppd_workspace)        !^ we always use SWpots once OVO
       else                                   !  metric got dropped
          call swri_obtain_swops_in_ppd_ptr_fast(swops_in_ppd_ptr, num_sws, & ! out
               src_centre, src_atom, cur_ppd, 'P', cell, swri, &              ! in
               swops_in_ppd_workspace, which_images(:,e))
       end if

       if(e == 1) centre_offset = 0
       if(e == 2) centre_offset = num_sws
       ! --------------------------------------------------------------------
       ! jd: for all SWs on src atom                                       SW
       ! --------------------------------------------------------------------
       sws_on_one:                                                          &
       do sw1_idx = 1, num_sws

          coeff = coeffs(sw1_idx + centre_offset)

          expansion_in_ppd(:) = expansion_in_ppd(:) + &
               coeff * swops_in_ppd_ptr(&
               (sw1_idx-1)*cell%n_pts+1:sw1_idx*cell%n_pts)
          ! ^ NB: we rely on the fact that expansion_in_ppd(:) is correctly
          !       dimensioned. Otherwise there'd be a bounds mismatch.

       end do sws_on_one

    end do two_centres

    ! jd: Debug two-centre expansions
    if(pub_hfx_debug .and. expansion_atoms(2) /= -1) then
       do ipt = 1, cell%n_pts
          if(which_images(ipt,1)%image /= which_images(ipt,2)%image) then

             ppd_corner = basis_ppd_location(cur_ppd,cell)

             im1_disp = minimum_image_number_to_displacement(&
                  which_images(ipt,1)%image, cell)
             im2_disp = minimum_image_number_to_displacement(&
                  which_images(ipt,2)%image, cell)

             call minimum_image_number_to_indices(which_images(ipt,1)%image, &
                  im1_a1_neighbour, im1_a2_neighbour, im1_a3_neighbour)
             call minimum_image_number_to_indices(which_images(ipt,2)%image, &
                  im2_a1_neighbour, im2_a2_neighbour, im2_a3_neighbour)

#if 1
             write(*,'(a,i3,a,i2,i2,i2,a,i3,a,i2,i2,i2,a,f14.6,f14.6,f14.6,a,&
                  &f14.6,f14.6,f14.6,a)') '[2E] Image 1: '//&
                  trim(utils_int_to_str(which_images(ipt,1)%image))//&
                  '['//trim(utils_int_to_str(im1_a1_neighbour))//&
                  ' '//trim(utils_int_to_str(im1_a2_neighbour))//&
                  ' '//trim(utils_int_to_str(im1_a3_neighbour))//&
                  '], Image 2: '//&
                  trim(utils_int_to_str(which_images(ipt,2)%image))//&
                  '['//trim(utils_int_to_str(im2_a1_neighbour))//&
                  ' '//trim(utils_int_to_str(im2_a2_neighbour))//&
                  ' '//trim(utils_int_to_str(im2_a3_neighbour))//'].'//&
                  ' Centre 1 at '//trim(utils_point_to_str(&
                  expansion_centres(1)%incell,'f14.6'))//'. '//&
                  'Centre 2 at '//trim(utils_point_to_str(&
                  expansion_centres(2)%incell,'f14.6'))//'.'
#endif

#if 1
             call utils_abort(myself//': [2E] Inconsistent images for point '//&
                  trim(utils_int_to_str(ipt))//' in PPD '//&
                  trim(utils_int_to_str(cur_ppd))//'.'//&
                  CRLF//'Image 1: '//&
                  trim(utils_int_to_str(which_images(ipt,1)%image))//&
                  '['//trim(utils_int_to_str(im1_a1_neighbour))//&
                  ' '//trim(utils_int_to_str(im1_a2_neighbour))//&
                  ' '//trim(utils_int_to_str(im1_a3_neighbour))//&
                  '], Image 2: '//&
                  trim(utils_int_to_str(which_images(ipt,2)%image))//&
                  '['//trim(utils_int_to_str(im2_a1_neighbour))//&
                  ' '//trim(utils_int_to_str(im2_a2_neighbour))//&
                  ' '//trim(utils_int_to_str(im2_a3_neighbour))//'].'//&
                  CRLF//'Centre 1 at '//trim(utils_point_to_str(&
                  expansion_centres(1)%incell,'f14.6'))//'.'//CRLF//&
                  'Centre 2 at '//trim(utils_point_to_str(&
                  expansion_centres(2)%incell,'f14.6'))//'.'//CRLF//&
                  'PPD starts at '//trim(utils_point_to_str(&
                  ppd_corner,'f14.6'))//'.'//CRLF//&
                  'Point at '//trim(utils_point_to_str(&
                  which_images(ipt,1)%curpoint,'f14.6'))//'.'//CRLF//&
                  'Direct displacement to centre 1: '//trim(utils_point_to_str(&
                  which_images(ipt,1)%disp0,'f14.6'))//'.'//CRLF//&
                  'Direct displacement to centre 2: '//trim(utils_point_to_str(&
                  which_images(ipt,2)%disp0,'f14.6'))//'.'//CRLF//&
                  'Displacement to image 1: '//trim(utils_point_to_str(&
                  which_images(ipt,1)%disp,'f14.6'))//'.'//CRLF//&
                  'Displacement to image 2: '//trim(utils_point_to_str(&
                  which_images(ipt,2)%disp,'f14.6'))//'.'//CRLF//&
                  'Direct distance to centre 1: '//trim(utils_real_to_str(&
                  which_images(ipt,1)%rad0))//'.'//CRLF//&
                  'Direct distance to centre 2: '//trim(utils_real_to_str(&
                  which_images(ipt,2)%rad0))//'.'//CRLF//&
                  'Distance to image 1: '//trim(utils_real_to_str(&
                  which_images(ipt,1)%rad))//'.'//CRLF//&
                  'Distance to image 2: '//trim(utils_real_to_str(&
                  which_images(ipt,2)%rad))//'.', avoid_mpi_calls = .true.)
#endif

          end if
       end do
    end if

  end subroutine hfx_calculate_single_expansion

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_calc_expansion_for_q(expansion_in_ppds, &                 ! out
       ppd_list, n_ppds, swri, cell, expansion_atoms, expansion_centres, & ! in
       q_coeffs)                                                           ! in
    !==========================================================================!
    ! Generates a potential coming from an SW expansion of Q term.             !
    ! The potential is generated in PPDs.                                      !
    !                                                                          !
    ! SWOPs in PPDs are cached behind the scenes. Q must have been expanded    !
    ! in advance (q_coeffs must be available).                                 !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   expansion_in_ppds (out): Result goes here.                             !
    !   ppd_list (in): The list of PPDs for which we want the expansion.       !
    !   n_ppds (in): Number of elements in ppd_list.                           !
    !   swri (inout): SW_RI object in which we do the resolution of identiy.   !
    !                 Needed for the SWOP cache.                               !
    !   swex_quality (in): Determines the quality of the expansion.            !
    !   cell (in): Needed for generating SWOPs.                                !
    !   expansion_atoms (in): Atoms in the expansion.                          !
    !   expansion_centres (in): Centres in the expansion.                      !
    !   q_coeffs (in): Expansion coefficients.                                 !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2019 using hfx_calc_expansions_bb_cc   !
    ! as a template.
    !==========================================================================!

    use basis, only: basis_ppd_location
    use constants, only: CRLF
    use geometry, only: POINT
    use rundat, only: pub_threads_max, pub_hfx_debug
    use simulation_cell, only: CELL_INFO, minimum_image_number_to_displacement,&
         minimum_image_number_to_indices
    use sw_resolution_of_identity, only: SW_RI, ATOM_CENTRE, PBC_IMAGE_INFO, &
         swri_obtain_swops_in_ppd_ptr_fast
    use timer, only: timer_clock
    use utils, only: utils_int_to_str, utils_point_to_str, utils_real_to_str, &
         utils_abort, utils_sanity_check

    implicit none

    ! jd: Arguments
    type(CELL_INFO), intent(in)   :: cell
    integer, intent(in)           :: n_ppds
    real(kind=DP), intent(out)    :: expansion_in_ppds(n_ppds*cell%n_pts)
    type(SW_RI), intent(in)       :: swri
    integer, intent(in)           :: expansion_atoms(2)
    type(ATOM_CENTRE), intent(in) :: expansion_centres(2)
    real(kind=DP), intent(in)     :: q_coeffs(2*swri%quality%max_sws_per_centre)
    integer, intent(in)           :: ppd_list(1:n_ppds)

    ! jd: Local variables
    real(kind=DP), pointer :: swops_in_ppd_ptr(:)
    real(kind=DP), target  :: swops_in_ppd_workspace(&
         cell%n_pts * swri%quality%num_sws_per_centre,0:pub_threads_max-1)
    integer           :: cur_ppd, cur_ppd_n
    integer           :: offs_to_cur_ppd
    integer           :: result
    logical           :: expansion_is_cached
    integer           :: size_of_cached_data
    integer           :: sw1_idx
    integer           :: e
    integer           :: centre_offset
    integer           :: src_atom
    integer           :: num_sws_on_src, num_sws
!$  integer, external :: omp_get_thread_num
    integer           :: tid
    real(kind=DP)     :: coeff
    type(ATOM_CENTRE) :: src_centre
    type(POINT)          :: ppd_corner
    type(POINT)          :: im1_disp, im2_disp
    integer              :: ipt
    integer              :: im1_a1_neighbour, im2_a1_neighbour
    integer              :: im1_a2_neighbour, im2_a2_neighbour
    integer              :: im1_a3_neighbour, im2_a3_neighbour
    type(PBC_IMAGE_INFO) :: which_images(cell%n_pts,2)
    character(len=*), parameter :: myself = 'hfx_calc_expansion_for_q'

    ! --------------------------------------------------------------------------

    call timer_clock(myself,1)

    expansion_in_ppds(:) = 0.0_DP

    ! --------------------------------------------------------------------------
    ! Go over relevant PPDs                                                  PPD
    ! --------------------------------------------------------------------------

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(DYNAMIC) &
!$OMP PRIVATE(cur_ppd_n, cur_ppd, offs_to_cur_ppd, expansion_is_cached, result,&
!$OMP      size_of_cached_data, e, src_atom, src_centre, tid,                  &
!$OMP      num_sws_on_src, num_sws, centre_offset, coeff, swops_in_ppd_ptr,    &
!$OMP      which_images, ppd_corner, im1_disp, im2_disp,                       &
!$OMP      im1_a1_neighbour, im1_a2_neighbour, im1_a3_neighbour,               &
!$OMP      im2_a1_neighbour, im2_a2_neighbour, im2_a3_neighbour)               &
!$OMP SHARED(n_ppds, ppd_list, cell, swri, expansion_atoms, expansion_centres, &
!$OMP      q_coeffs, expansion_in_ppds, swops_in_ppd_workspace, pub_hfx_debug)
    all_ppds:                                                                  &
    do cur_ppd_n = 1, n_ppds

       cur_ppd = ppd_list(cur_ppd_n)
       offs_to_cur_ppd = (cur_ppd_n-1)*cell%n_pts+1

       tid = 0
!$     tid = omp_get_thread_num()

       ! jd: Loop over unique centres in the expansion
       two_centres:                                                            &
       do e = 1, 2
          src_atom = expansion_atoms(e)
          src_centre = expansion_centres(e)
          if(src_atom == -1) exit

          num_sws_on_src = swri%quality%num_sws_per_centre

          ! jd: Get these SWOPs from hash table or calculate, if absent.
          !     This populates the hash table (from a critical section) and
          !     temporarily locks it against evictions.
          if(.not. pub_hfx_debug) then
             call swri_obtain_swops_in_ppd_ptr_fast(swops_in_ppd_ptr, num_sws, & ! out
                  src_centre, src_atom, cur_ppd, 'P', cell, swri, &              ! in
                  swops_in_ppd_workspace(:,tid)) !^ we always use SWpots once OVO
          else                                   !  metric got dropped
             call swri_obtain_swops_in_ppd_ptr_fast(swops_in_ppd_ptr, num_sws, & ! out
                  src_centre, src_atom, cur_ppd, 'P', cell, swri, &              ! in
                  swops_in_ppd_workspace(:,tid), which_images(:,e))
          end if

          if(e == 1) centre_offset = 0
          if(e == 2) centre_offset = num_sws_on_src
          ! --------------------------------------------------------------------
          ! jd: for all SWs on src atom                                       SW
          ! --------------------------------------------------------------------
          sws_on_one:                                                          &
          do sw1_idx = 1, num_sws_on_src

             coeff = q_coeffs(sw1_idx + centre_offset)

#if 0
             call utils_sanity_check((/coeff/),'single_q_coeff',excessive=1D12)
             call utils_sanity_check(swops_in_ppd_ptr(&
                  (sw1_idx-1)*cell%n_pts+1:sw1_idx*cell%n_pts),&
                  'swops_in_ppd_ptr',excessive=1D12)
#endif

             expansion_in_ppds(offs_to_cur_ppd:offs_to_cur_ppd+cell%n_pts-1) = &
                  expansion_in_ppds(offs_to_cur_ppd:offs_to_cur_ppd+cell%n_pts-1) + &
                  coeff * swops_in_ppd_ptr(&
                  (sw1_idx-1)*cell%n_pts+1:sw1_idx*cell%n_pts)

!            jd: ^ Equivalent BLAS call, practically identical performance
!            call linalg_daxpy(cell%n_pts, coeff, swops_in_a_ppd_ptr(&
!                 (sw1_idx-1)*cell%n_pts+1:sw1_idx*cell%n_pts), 1, expansion_in_ppd, 1)

          end do sws_on_one

       end do two_centres

       ! jd: Debug two-centre expansions
       if(pub_hfx_debug .and. expansion_atoms(2) /= -1) then
          do ipt = 1, cell%n_pts
             if(which_images(ipt,1)%image /= which_images(ipt,2)%image) then

                ppd_corner = basis_ppd_location(cur_ppd,cell)

                im1_disp = minimum_image_number_to_displacement(&
                     which_images(ipt,1)%image, cell)
                im2_disp = minimum_image_number_to_displacement(&
                     which_images(ipt,2)%image, cell)

                call minimum_image_number_to_indices(which_images(ipt,1)%image, &
                     im1_a1_neighbour, im1_a2_neighbour, im1_a3_neighbour)
                call minimum_image_number_to_indices(which_images(ipt,2)%image, &
                     im2_a1_neighbour, im2_a2_neighbour, im2_a3_neighbour)

#if 1
             write(*,'(a,i3,a,i2,i2,i2,a,i3,a,i2,i2,i2,a,f14.6,f14.6,f14.6,a,&
                  f14.6,f14.6,f14.6,a)') '[3B] Image 1: '//&
                  trim(utils_int_to_str(which_images(ipt,1)%image))//&
                  '['//trim(utils_int_to_str(im1_a1_neighbour))//&
                  ' '//trim(utils_int_to_str(im1_a2_neighbour))//&
                  ' '//trim(utils_int_to_str(im1_a3_neighbour))//&
                  '], Image 2: '//&
                  trim(utils_int_to_str(which_images(ipt,2)%image))//&
                  '['//trim(utils_int_to_str(im2_a1_neighbour))//&
                  ' '//trim(utils_int_to_str(im2_a2_neighbour))//&
                  ' '//trim(utils_int_to_str(im2_a3_neighbour))//'].'//&
                  ' Centre 1 at '//trim(utils_point_to_str(&
                  expansion_centres(1)%incell,'f14.6'))//'. '//&
                  'Centre 2 at '//trim(utils_point_to_str(&
                  expansion_centres(2)%incell,'f14.6'))//'.'
#endif

#if 1
                call utils_abort(myself//&
                     ': [2] Inconsistent images for point '//&
                     trim(utils_int_to_str(ipt))//' in PPD '//&
                     trim(utils_int_to_str(cur_ppd))//'.'//&
                     CRLF//'Image 1: '//&
                     trim(utils_int_to_str(which_images(ipt,1)%image))//&
                     '['//trim(utils_int_to_str(im1_a1_neighbour))//&
                     ' '//trim(utils_int_to_str(im1_a2_neighbour))//&
                     ' '//trim(utils_int_to_str(im1_a3_neighbour))//&
                     '], Image 2: '//&
                     trim(utils_int_to_str(which_images(ipt,2)%image))//&
                     '['//trim(utils_int_to_str(im2_a1_neighbour))//&
                     ' '//trim(utils_int_to_str(im2_a2_neighbour))//&
                     ' '//trim(utils_int_to_str(im2_a3_neighbour))//'].'//&
                     CRLF//'Centre 1 at '//trim(utils_point_to_str(&
                     expansion_centres(1)%incell,'f14.6'))//'.'//CRLF//&
                     'Centre 2 at '//trim(utils_point_to_str(&
                     expansion_centres(2)%incell,'f14.6'))//'.'//CRLF//&
                     'PPD starts at '//trim(utils_point_to_str(&
                     ppd_corner,'f14.6'))//'.'//CRLF//&
                     'Point at '//trim(utils_point_to_str(&
                     which_images(ipt,1)%curpoint,'f14.6'))//'.'//CRLF//&
                     'Direct displacement to centre 1: '//trim(utils_point_to_str(&
                     which_images(ipt,1)%disp0,'f14.6'))//'.'//CRLF//&
                     'Direct displacement to centre 2: '//trim(utils_point_to_str(&
                     which_images(ipt,2)%disp0,'f14.6'))//'.'//CRLF//&
                     'Displacement to image 1: '//trim(utils_point_to_str(&
                     which_images(ipt,1)%disp,'f14.6'))//'.'//CRLF//&
                     'Displacement to image 2: '//trim(utils_point_to_str(&
                     which_images(ipt,2)%disp,'f14.6'))//'.'//CRLF//&
                     'Direct distance to centre 1: '//trim(utils_real_to_str(&
                     which_images(ipt,1)%rad0))//'.'//CRLF//&
                     'Direct distance to centre 2: '//trim(utils_real_to_str(&
                     which_images(ipt,2)%rad0))//'.'//CRLF//&
                     'Distance to image 1: '//trim(utils_real_to_str(&
                     which_images(ipt,1)%rad))//'.'//CRLF//&
                     'Distance to image 2: '//trim(utils_real_to_str(&
                     which_images(ipt,2)%rad))//'.', avoid_mpi_calls = .true.)
#endif

             end if
          end do
       end if

    end do all_ppds
!$OMP END PARALLEL DO

    call timer_clock(myself,2)

  end subroutine hfx_calc_expansion_for_q

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_contract_over_dkn_cc(hfxstate, f, a_ppd_list, n_a_ppd_list, &
       packed_dkn_cd_blocks_ht, bb_cc_symmetry, global_b, ngwf_indices, &
       cell, rep, par, src_swex_h, elements, &
       n_expansions_hits, n_expansions_misses) ! adds to
    !==========================================================================!
    ! Contracts the SWOP expansions with K^CcDd, removing the Cc index from    !
    ! further consideration.                                                   !
    !                                                                          !
    ! For all Cc that are S-neighbours with current B and are in my BC pairs,  !
    ! for all PPDs that are of interest (that is, belonging to A's that are    !
    ! X-neighbours with current B), calculate                                  !
    !                                                                          !
    ! f_bDdsi;B = \sum_Cc ( e_bCci;B * K^CcDds ).                              !
    !                                                                          !
    ! e are SWOP expansions, f is the result, i is the PPD index, s is spin.   !
    ! ";B" indicates parametric-like dependence on B -- i.e. we only perform   !
    ! the operation for a "current B" provided as an argument.                 !
    !                                                                          !
    ! In a perfect world with a lot of RAM this subroutine extracts all e from !
    ! expansions_ht as they have been calculated in advance and cached. This   !
    ! enables reusing e over inner loop iterations (they are DKN-independent). !
    ! In the sad reality of 2019 we find that for systems with hundreds of     !
    ! atoms we cannot cache all the expansions. All the missing expansions are !
    ! calculated here on the fly.
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (in): Needed for neighbour list, atom lists.                  !
    !   f (out): Result goes here. It gets allocated here.                     !
    !            1st index - fused PPD-D-d-pt_in_ppd.                          !
    !            2nd index - spin.                                             !
    !            3rd index - b.                                                !
    !   a_ppd_list (in): Indices of PPDs that are of interest (see above).     !
    !   n_a_ppd_list (in): Number of elements in the above.                    !
    !   packed_dkn_cd_blocks_ht (in): HT from which DKN^CcDds are read.        !
    !   bb_cc_symmetry (in): If .true., it will be assumed e is invariant wrt  !
    !                        swapping Bb and Cc (like in valence calcs).       !
    !   global_b (in): Global index of B atom for which this is calculated.    !
    !   ngwf_indices (in): Array of HFX_NGWF_INDEX instances containing        !
    !                      pointers to the NGWF basis and the NGWF_REP used    !
    !                      by atoms A, B, C and D (only FUNC_BASIS instances   !
    !                      for atoms B, C and D used).                         !
    !   cell (in): Needed mostly for n_pts, but also in call to                !
    !              hfx_calculate_single_expansion() that becomes necessary if  !
    !              we have to expansions on the fly.                           !
    !   rep (in): Needed for swexes, where coeffs_ht lives, needed for         !
    !             expansions on the fly.                                       !
    !   par (in): Needed for swri_init_centre() for expansions on the fly.     !
    !   src_swex_h (in): Needed for indexing swexes when doing expansions on   !
    !                    the fly.                                              !
    !   elements (in): Needed for swri_init_centre() for expansions on the fly.!
    !   n_expansions_hits (inout):   } performance counters for expansions HT, !
    !   n_expansions_misses (inout): } this subroutine adds to them.           !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   It is the responsibility of the caller to deallocate f.                !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2019.                                 !
    !==========================================================================!

    use constants, only: garbage_real, METRIC_ELECTROSTATIC, METRIC_OVERLAP, &
         SW_V, SW_O, CRLF
    use hash_table, only: HT_HASH_TABLE, hash_table_lookup_ptr_nocount, &
         hash_table_lookup_nocount
    use ion, only: ELEMENT
    use ngwf_representation, only: NGWF_REP
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_num_spins, pub_threads_max, pub_hfx_metric, &
         pub_cache_limit_for_expansions
    use simulation_cell, only: CELL_INFO
    use sw_resolution_of_identity, only: ATOM_CENTRE, swri_expansion_centres, &
         swri_init_centre, swri_library
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_abort, utils_int_to_str, &
         utils_alloc_check, utils_dealloc_check, utils_postfix_to_ngwf_set_name

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(in), target     :: hfxstate
    real(kind=DP), intent(out), allocatable :: f(:,:,:)
    integer, intent(in)                     :: a_ppd_list(:)
    integer, intent(in)                     :: n_a_ppd_list
    type(HT_HASH_TABLE), intent(in), target :: packed_dkn_cd_blocks_ht
    logical, intent(in)                     :: bb_cc_symmetry
    integer, intent(in)                     :: global_b
    type(CELL_INFO), intent(in)             :: cell
    type(NGWF_REP), intent(in)              :: rep
    type(PARAL_INFO), intent(in), pointer   :: par
    integer, intent(in)                     :: src_swex_h
    type(ELEMENT), intent(in)               :: elements(:)
    type(HFX_NGWF_INDEX), intent(in)        :: ngwf_indices(A_NGWF:D_NGWF)
    real(kind=DP), intent(inout)            :: n_expansions_hits
    real(kind=DP), intent(inout)            :: n_expansions_misses

    ! jd: Local variables
    real(kind=DP), allocatable :: ppd_a_swpot_sum(:,:,:)
    real(kind=DP), allocatable :: dkn_blocks_cd(:,:,:,:) ! c, d, is, d_idx
    real(kind=DP), allocatable :: dkn_blocks_cd_reordered(:,:,:) ! is, global_Dd, c
    integer                    :: dlist(max_num_at_per_ppd)
    integer                    :: dactual(n_a_ppd_list, max_num_at_per_ppd)
    integer                    :: ndactual(n_a_ppd_list)
    integer                    :: offsets_to_ppds(n_a_ppd_list)
    real(kind=DP), pointer     :: dlist_real_ptr(:)
    real(kind=DP)              :: coeffs(&
         2*swri_library(hfxstate%hfx_swri_h)%quality%max_sws_per_centre)
    real(kind=DP) :: dkn_cc_dd_is
    integer :: offset, offset_save_b, offset_save_c
    integer :: packed_dkn_block_size
    integer :: global_c, global_d
    integer :: d_idx
    integer :: d_ord
    integer :: num_d
    integer :: ngwf_b, ngwf_c, ngwf_d
    integer :: global_bb_ngwf_idx, global_cc_ngwf_idx, global_dd_ngwf_idx
    integer :: first_ngwf_idx_of_b, first_ngwf_idx_of_c, first_ngwf_idx_of_d
    integer :: n_ngwfs_b, n_ngwfs_c, n_ngwfs_d
    integer :: global_ngwf_idx_1, global_ngwf_idx_2
    integer :: i_ppd
    integer :: cur_ppd
    integer :: size_of_cached_data
    integer :: num_of_cached_coeffs
    integer :: fused_index_range
    integer :: npts
    integer :: is
    integer :: pair_idx
    integer :: i_d_actual
    integer :: num_sws_in_expansion
    integer :: atoms(2), expansion_atoms(2)
    type(ATOM_CENTRE) :: centres(2), expansion_centres(2)
    integer :: coeffs_kind
    integer :: ierr
    integer :: tid
!$  integer, external :: omp_get_thread_num
    character(len=*), parameter :: myself = 'hfx_contract_over_dkn_cc'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    ! jd: Sanity check on neighbour lists
    call utils_assert(hfxstate%xs_atoms_nl%populated, &
         myself//': xs_atoms neighbour list not populated.')

    packed_dkn_block_size = pub_num_spins * &
         ngwf_indices(C_NGWF)%basis%max_on_atom * &
         ngwf_indices(D_NGWF)%basis%max_on_atom

    first_ngwf_idx_of_b = ngwf_indices(B_NGWF)%basis%first_on_atom(global_b)
    n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)

    npts = cell%n_pts

    ! **************************************************************************
    ! * Establish range of fused index PPD-D-d-pt.
    ! * Relevant D's change depending on PPD, so we don't know in advance.
    ! * Also, figure out which D's share this PPD, but do not actually intersect
    ! * with any A that is X-neighbour of this B (D's that are not x-S-neighbours
    ! * of B). These have to be cycled over. We store 'good' D's in dactual.
    ! **************************************************************************

    fused_index_range = 0
    offsets_to_ppds(:) = garbage_int
    ! --------------------------------------------------------------------------
    ! jd: All my A PPD *relevant to current B*                             r PPD
    ! --------------------------------------------------------------------------
    loop_relev_PPD:                                                            &
    do i_ppd = 1, n_a_ppd_list
       cur_ppd = a_ppd_list(i_ppd)

       ! jd: Retrieve list of all D's sharing this PPD
       call hash_table_lookup_ptr_nocount(dlist_real_ptr, num_d, &
            hfxstate%dlists_ht, global_b, cur_ppd)
       call utils_assert(num_d <= max_num_at_per_ppd, myself//&
            ': Insufficient max_num_at_per_ppd [1]', &
            max_num_at_per_ppd, num_d)
       call utils_assert(num_d /= -1, myself//': dlist not found [1]', &
            global_b, cur_ppd)
       call utils_assert(num_d /= 0, myself//': dlist empty [1]', &
            global_b, cur_ppd)

       ! jd: Convert to ints, pad
       dlist(1:num_d) = nint(dlist_real_ptr(1:num_d))
       dlist(num_d+1:) = garbage_int

       ! -----------------------------------------------------------------------
       ! jd: For all D's that are relevant to this PPD                     r DDD
       ! -----------------------------------------------------------------------
       i_d_actual = 1
       dactual(i_ppd,:) = -1
       loop_relev_D:                                                           &
       do d_idx = 1, num_d
          global_d = dlist(d_idx)

          ! jd: Careful. There are corner cases where D shares a PPD with
          !     an A that is an X-neighbour of B, but A and D *do not*
          !     intersect (from the sparsity point of view). For every
          !     candidate D we must double-check if it's an x-S-neighbour of B.
          if(.not. any(hfxstate%xs_atoms_nl%neighbours(&
               hfxstate%xs_atoms_nl%first_idx(global_b):&
               hfxstate%xs_atoms_nl%last_idx(global_b)) == global_d)) then
             cycle
          end if

          n_ngwfs_d = ngwf_indices(D_NGWF)%basis%num_on_atom(global_d)
          ! jd: Store offset to PPD when we hit the 1st non-cycled-over D
          !     Can't do this for d_idx == 1, because this D can get cycled over
          if(offsets_to_ppds(i_ppd) == garbage_int) then
             offsets_to_ppds(i_ppd) = fused_index_range + 1
          end if
          fused_index_range = fused_index_range + n_ngwfs_d * npts

          if(fused_index_range < 0) then
             call utils_abort(myself//': Overflow in auxiliary quantity F.'//&
                  CRLF//'Increase the number of MPI ranks (processors) in &
                  &your calculation (preferably), or reduce your PPD size. &
                  &As a last resort, reduce the exchange cutoff (hfx_cutoff).')
          end if

          dactual(i_ppd,i_d_actual) = global_d
          i_d_actual = i_d_actual + 1
       end do loop_relev_D
       ndactual(i_ppd) = i_d_actual - 1

    end do loop_relev_PPD

    ! jd: ~1GB for 443 atoms on 8 ranks with 4 NGWFs max.
    !     If this turns out too big (conduction?), consider moving n_ngwfs_b
    !     out of the subroutine.
    allocate(f(fused_index_range, pub_num_spins, n_ngwfs_b),stat=ierr)
    call utils_alloc_check(myself,'f',ierr)
    f(:,:,:) = 0.0_DP

    ! --------------------------------------------------------------------------
    ! jd: Loop over C's that are in my B-C pair list for this B              CCC
    ! --------------------------------------------------------------------------
    loop_C:                                                                    &
    do pair_idx = 1, hfxstate%n_my_bc_pairs
       if(hfxstate%my_bc_pairs(1,pair_idx) /= global_b) cycle ! not this B
       global_c = hfxstate%my_bc_pairs(2,pair_idx)

       n_ngwfs_c = ngwf_indices(C_NGWF)%basis%num_on_atom(global_c)
       first_ngwf_idx_of_c = ngwf_indices(C_NGWF)%basis%first_on_atom(global_c)

       allocate(ppd_a_swpot_sum(n_a_ppd_list * npts, n_ngwfs_c, n_ngwfs_b), &
            stat=ierr)
       call utils_alloc_check(myself,'ppd_a_swpot_sum',ierr)

       ! **********************************************************************
       ! *** EXPANSIONS ON THE FLY
       ! **********************************************************************

       ! jd: If we couldn't fit all the expansions and will have to recalculate
       !     some, some preparation is in order
       if(.not. hfxstate%all_expansions_cached) then
          ! jd: Figure out the expansion centres, sort them,
          !     pad with -1 if fewer than 2 are needed
          call swri_init_centre(centres(1), par, elements, global_b)
          call swri_init_centre(centres(2), par, elements, global_c)
          atoms = (/global_b, global_c/)
          call swri_expansion_centres(swri_library(hfxstate%hfx_swri_h), &
               expansion_centres, expansion_atoms, num_sws_in_expansion, & ! out
               centres, atoms, swri_library(hfxstate%hfx_swri_h)%quality, &! in
               bb_cc_symmetry)                                             ! in
       end if

       ! **********************************************************************
       ! *** END EXPANSIONS ON THE FLY
       ! **********************************************************************

       ! ***********************************************************************
       ! * Retrieve or calculate on the fly all necessary expansions for
       ! * all my A PPDs.
       ! ***********************************************************************
       call timer_clock(myself//'_expa_get',1)

       ! -----------------------------------------------------------------------
       ! jd: for all b on B                                                  bbb
       ! -----------------------------------------------------------------------
       loop_ngwf_b:                                                            &
       do ngwf_b = 1, n_ngwfs_b
          global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

          ! --------------------------------------------------------------------
          ! jd: for all c on C                                               ccc
          ! --------------------------------------------------------------------
          loop_ngwf_c:                                                         &
          do ngwf_c = 1, n_ngwfs_c
             global_cc_ngwf_idx = first_ngwf_idx_of_c + ngwf_c - 1

             if(.not. bb_cc_symmetry) then
                global_ngwf_idx_1 = global_bb_ngwf_idx
                global_ngwf_idx_2 = global_cc_ngwf_idx
             else
                global_ngwf_idx_1 = min(global_bb_ngwf_idx,global_cc_ngwf_idx)
                global_ngwf_idx_2 = max(global_bb_ngwf_idx,global_cc_ngwf_idx)
             end if

             ! ****************************************************************
             ! *** EXPANSIONS ON THE FLY
             ! ****************************************************************

             ! jd: If we couldn't fit all the expansions and will have to
             !     recalculate some, some preparation is in order
             if(.not. hfxstate%all_expansions_cached) then

                if(pub_hfx_metric == METRIC_ELECTROSTATIC) coeffs_kind = SW_V
                if(pub_hfx_metric == METRIC_OVERLAP) coeffs_kind = SW_O

                ! **************************************************************
                ! jd: Look up expansion coefficients for this Bb Cc
                ! **************************************************************
                if(.not. bb_cc_symmetry) then
                   ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                   call hash_table_lookup_nocount(coeffs, num_of_cached_coeffs,&
                        rep%swexes(src_swex_h)%coeffs_ht, &
                        global_bb_ngwf_idx, global_cc_ngwf_idx, coeffs_kind, &
                        src_swex_h, 1) ! 1 is for formal spin-dependence
                else
                   ! jd: Here there's Bb-Cc symmetry and we only store one triangle
                   call hash_table_lookup_nocount(coeffs, num_of_cached_coeffs,&
                        rep%swexes(src_swex_h)%coeffs_ht, &
                        min(global_bb_ngwf_idx,global_cc_ngwf_idx), &
                        max(global_bb_ngwf_idx,global_cc_ngwf_idx), coeffs_kind, &
                        src_swex_h, 1) ! 1 is for formal spin-dependence
                end if

                call utils_assert(num_of_cached_coeffs /= -1, myself//&
                     ': Internal error: coeffs hash table too small or you &
                     &forgot to call hf_exchange_dkn_indep_stage() on '//&
                     trim(utils_postfix_to_ngwf_set_name(rep%postfix))//'%['//&
                     trim(rep%swexes(src_swex_h)%swex_name)//']: &
                     &Bb, Cc, coeffs_kind follow: ', &
                     global_bb_ngwf_idx, global_cc_ngwf_idx, coeffs_kind)
             end if

             ! ****************************************************************
             ! *** END EXPANSIONS ON THE FLY
             ! ****************************************************************

             ! -----------------------------------------------------------------
             ! jd: Get expansions for all my A PPDs for current Bb, Cc       PPD
             ! -----------------------------------------------------------------
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(DYNAMIC) &
!$OMP PRIVATE(i_ppd, cur_ppd, size_of_cached_data) &
!$OMP SHARED(hfxstate, ppd_a_swpot_sum, global_ngwf_idx_1, global_ngwf_idx_2,  &
!$OMP      npts, ngwf_b, ngwf_c, a_ppd_list, n_a_ppd_list, expansion_atoms,    &
!$OMP      expansion_centres, cell, pub_cache_limit_for_expansions, coeffs,    &
!$OMP      swri_library)                                                       &
!$OMP REDUCTION(+:n_expansions_hits,n_expansions_misses)
             loop_PPD:                                                         &
             do i_ppd = 1, n_a_ppd_list
                cur_ppd = a_ppd_list(i_ppd)

                ! jd: Retrieve from HT to buffer. Threads access different
                !     slices of 1st idx of ppd_a_swpot_sum, so no race.
                !     If we did not cache *any* expansions, don't bother looking.
                if(hfxstate%expansions_ht_in_use) then
                   call hash_table_lookup_nocount(&
                        ppd_a_swpot_sum((i_ppd-1)*npts+1:,ngwf_c,ngwf_b), &
                        size_of_cached_data, hfxstate%expansions_ht, cur_ppd, &
                        global_ngwf_idx_1, global_ngwf_idx_2)
                else
                   size_of_cached_data = -1
                end if

                ! jd: If we didn't have enough RAM for all the expansions,
                !     we will need to recalculate some.
                if(size_of_cached_data == -1) then
                   if(.not. hfxstate%all_expansions_cached) then
                      call hfx_calculate_single_expansion(ppd_a_swpot_sum(&
                           (i_ppd-1)*npts+1:i_ppd*npts,ngwf_c,ngwf_b), &
                           expansion_atoms, expansion_centres, cell, &
                           swri_library(hfxstate%hfx_swri_h), coeffs, cur_ppd)
                      n_expansions_misses = n_expansions_misses + 1.0_DP
                   else
                      call utils_abort(myself//': Expansion not found in HT.', &
                           cur_ppd, global_ngwf_idx_1, global_ngwf_idx_2)
                   end if
                else if(size_of_cached_data /= npts) then
                   call utils_abort(myself//': Mismatch between the size of &
                        &cached expansion ('//trim(utils_int_to_str(&
                        size_of_cached_data))//') and PPD size ('//&
                        trim(utils_int_to_str(npts))//').', &
                        avoid_mpi_calls = .true.)
                else
                   n_expansions_hits = n_expansions_hits + 1.0_DP
                end if

             end do loop_PPD
!$OMP END PARALLEL DO

          end do loop_ngwf_c

       end do loop_ngwf_b

       call timer_clock(myself//'_expa_get',2)

       call timer_clock(myself//'_dkn_get',1)

       ! ***********************************************************************
       ! * Get all necessary DKN elements (all c, all s, all my Dd) from HTs
       ! ***********************************************************************
       ! @room for opt: this is not directly dependent on B (but indirectly,
       ! since C=f(B)).

       allocate(dkn_blocks_cd(ngwf_indices(C_NGWF)%basis%max_on_atom, &
            ngwf_indices(D_NGWF)%basis%max_on_atom, pub_num_spins, &
            hfxstate%n_my_d_atoms),stat=ierr)
       call utils_alloc_check(myself,'dkn_block_cd',ierr)
       dkn_blocks_cd(:,:,:,:) = garbage_real ! helps uncover padding bugs
       allocate(dkn_blocks_cd_reordered(pub_num_spins,&
            ngwf_indices(D_NGWF)%basis%num,&
            ngwf_indices(C_NGWF)%basis%max_on_atom),stat=ierr)
       call utils_alloc_check(myself,'dkn_block_cd_reordered',ierr)
       dkn_blocks_cd_reordered(:,:,:) = garbage_real ! helps uncover padding bugs

       ! -----------------------------------------------------------------------
       ! jd: Loop over D's that are x-S-neighbours of B                      DDD
       ! -----------------------------------------------------------------------
       loop_D:                                                                 &
       do d_idx = hfxstate%xs_atoms_nl%first_idx(global_b), &
            hfxstate%xs_atoms_nl%last_idx(global_b)
          global_d = hfxstate%xs_atoms_nl%neighbours(d_idx)
          d_ord = d_idx - hfxstate%xs_atoms_nl%first_idx(global_b) + 1
          first_ngwf_idx_of_d = ngwf_indices(D_NGWF)%basis%first_on_atom(global_d)
          n_ngwfs_d = ngwf_indices(D_NGWF)%basis%num_on_atom(global_d)

          ! jd: Get DKN^CD (all c's, all d's, spins)
          call hfx_get_dkn_atomblock(dkn_blocks_cd(:,:,:,d_ord), & ! out
               packed_dkn_cd_blocks_ht, global_c, global_d, &
               ngwf_indices(C_NGWF)%basis%max_on_atom, &
               ngwf_indices(D_NGWF)%basis%max_on_atom)

          ! jd: Reorder DKN elements for cache efficiency.
          do ngwf_c = 1, n_ngwfs_c
             do ngwf_d = 1, n_ngwfs_d
                global_dd_ngwf_idx = first_ngwf_idx_of_d + ngwf_d - 1
                do is = 1, pub_num_spins
                   dkn_blocks_cd_reordered(is,global_dd_ngwf_idx,ngwf_c) = &
                        dkn_blocks_cd(ngwf_c,ngwf_d,is,d_ord)
                end do
             end do
          end do

       end do loop_D

       deallocate(dkn_blocks_cd,stat=ierr)
       call utils_dealloc_check(myself,'dkn_block_cd',ierr)

       call timer_clock(myself//'_dkn_get',2)
       call timer_clock(myself//'_main',1)

       ! -----------------------------------------------------------------------
       ! jd: All my A PPD *relevant to current B*                          r PPD
       ! -----------------------------------------------------------------------
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(DYNAMIC) &
!$OMP PRIVATE(i_ppd, cur_ppd, first_ngwf_idx_of_d, global_d, n_ngwfs_d, &
!$OMP      dkn_cc_dd_is, global_dd_ngwf_idx, tid, offset_save_b, &
!$OMP      offset_save_c) &
!$OMP FIRSTPRIVATE(npts, n_ngwfs_b, n_ngwfs_c) &
!$OMP LASTPRIVATE(offset) &
!$OMP SHARED(hfxstate, a_ppd_list, n_a_ppd_list, ngwf_indices, f, &
!$OMP      ndactual, dactual, dkn_blocks_cd_reordered, pub_num_spins, &
!$OMP      ppd_a_swpot_sum, offsets_to_ppds)
       loop_relev_PPD2:                                                        &
       do i_ppd = 1, n_a_ppd_list
          cur_ppd = a_ppd_list(i_ppd)

          tid = 0
!$        tid = omp_get_thread_num()

          offset = offsets_to_ppds(i_ppd)

          ! --------------------------------------------------------------------
          ! jd: For all D's that are relevant to this PPD                  r DDD
          !     ^ with D's that share a PPD with A, but are not S-X-neighbours
          !       of current B eliminated.
          ! --------------------------------------------------------------------
          loop_actual_D2:                                                      &
          do d_idx = 1, ndactual(i_ppd)
             global_d = dactual(i_ppd,d_idx)
             first_ngwf_idx_of_d = &
                  ngwf_indices(D_NGWF)%basis%first_on_atom(global_d)
             n_ngwfs_d = ngwf_indices(D_NGWF)%basis%num_on_atom(global_d)

             ! -----------------------------------------------------------------
             ! jd: For all d                                                 ddd
             ! -----------------------------------------------------------------
             loop_ngwf_d:                                                      &
             do ngwf_d = 1, n_ngwfs_d
                global_dd_ngwf_idx = first_ngwf_idx_of_d + ngwf_d - 1

                ! --------------------------------------------------------------
                ! jd: For all b                                              bbb
                ! --------------------------------------------------------------
                loop_ngwf_b2:                                                  &
                do ngwf_b = 1, n_ngwfs_b

                   if(ngwf_b == 1) then
                      offset_save_b = offset
                   else
                      offset = offset_save_b
                   end if

                   ! -----------------------------------------------------------
                   ! jd: For all c                                           ccc
                   ! -----------------------------------------------------------
                   loop_ngwf_c2:                                               &
                   do ngwf_c = 1, n_ngwfs_c

                      if(ngwf_c == 1) then
                         offset_save_c = offset
                      else
                         offset = offset_save_c
                      end if

                      ! --------------------------------------------------------
                      ! jd: For all spins                                    sss
                      ! --------------------------------------------------------
                      loop_spins:                                              &
                      do is = 1, pub_num_spins
                         dkn_cc_dd_is = dkn_blocks_cd_reordered(&
                              is, global_dd_ngwf_idx, ngwf_c)

                         f(offset:offset+npts-1,is,ngwf_b) = &
                              f(offset:offset+npts-1,is,ngwf_b) + &
                              ppd_a_swpot_sum((i_ppd-1)*npts+1:i_ppd*npts,&
                              ngwf_c,ngwf_b) * dkn_cc_dd_is

                      end do loop_spins

                      offset = offset + npts

                   end do loop_ngwf_c2

                end do loop_ngwf_b2

             end do loop_ngwf_d

             ! *************************************************************

          end do loop_actual_D2

          if(i_ppd < n_a_ppd_list) then
             if(offset /= offsets_to_ppds(i_ppd+1)) then
                call utils_abort(myself//&
                     ': Book-keeping error writing f, thread', tid, &
                     offset, offsets_to_ppds(i_ppd+1), avoid_mpi_calls = .true.)
             end if
          end if

       end do loop_relev_PPD2
!$OMP END PARALLEL DO

       call timer_clock(myself//'_main',2)

       call utils_assert(offset-1 == ubound(f,1), &
            myself//': Book-keeping error writing f.', offset-1, ubound(f,1))

       deallocate(dkn_blocks_cd_reordered,stat=ierr)
       call utils_dealloc_check(myself,'dkn_block_cd_reordered',ierr)

       deallocate(ppd_a_swpot_sum,stat=ierr)
       call utils_dealloc_check(myself,'ppd_a_swpot_sum',ierr)

    end do loop_C

    call timer_clock(myself,2)

  end subroutine hfx_contract_over_dkn_cc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_all_ad_for_this_b(x_ab_contributions_ht, &
       hfxstate, f, a_ppd_list, n_a_ppd_list, global_b, cell, global_nat, &
       fftbox_weight, ngwf_indices, n_prods_hits, n_prods_misses, calcgradient,&
       grad_ket_in_ppds)
    !==========================================================================!
    ! Calculates an X matrix stripe for a given B. The columns in the stripe   !
    ! count all NGWFs b on current B. Rows are all rows (Aa), but values in    !
    ! the stripe will be zero where A is not an X-neighbour of B.              !
    !                                                                          !
    ! This stage of the calculation eliminates the final indices, Dd and i_PPD.!
    !                                                                          !
    ! What we do is:                                                           !
    ! for all PPDs relevant to current B {                                     !
    !   for all D's that share this PPD {                                      !
    !     for all A's that share this PPD {                                    !
    !       for all d, b, a {                                                  !
    !         for all spins {                                                  !
    !            Xstripe(spin,Aa,b) += sum_pts (f_bDdsi;B * ppd_product(Aa,Dd))!
    !         }                                             ^ for current PPD  !
    !       }                                                                  !
    !     }                                                                    !
    !   }                                                                      !
    ! }                                                                        !
    !                                                                          !
    ! By "relevant PPDs" we mean PPDs in the union of all A's that are         !
    ! X-neighbours with current B.                                             !
    !                                                                          !
    ! Results are stored in x_ab_contributions_ht, they are only merged into   !
    ! X in hfx_comms_matrix().                                                 !
    !                                                                          !
    ! If calcgradient is .true., we also calculate NGWF gradient term 1 kets:  !
    ! k_{s,b,i_PPD;B} = \sum_{Dd} \phi_{Dd,i_PPD} * f_{bDdsi_PPD;B}.           !
    ! These are calculated for the set of PPDs shared by all relevant A's.     !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   x_ab_contributions_ht (inout): Results go here. They are reordered     !
    !                                  slightly and packed for later conve-    !
    !                                  nience.                                 !
    !   hfxstate (in): Book-keeping of neighbour lists, alists, dlists and     !
    !                  PPD products HT.                                        !
    !   f (in): Expansions contracted with DKN, calculated in                  !
    !           hfx_contract_over_dkn_cc().                                    !
    !   a_ppd_list (in): PPD indices of relevant PPDs.                         !
    !   n_a_ppd_list (in): Number of elements in the above.                    !
    !   global_b (in): Global index of current B atom.                         !
    !   cell (in): Needed only for n_pts and to give more helpful error msgs.  !
    !   global_nat (in): size(elements). Needed to dimension array.            !
    !   fftbox_weight (in): Needed in prefactor.                               !
    !   ngwf_indices (in): Array of HFX_NGWF_INDEX instances containing        !
    !                      pointers to the NGWF basis and the NGWF_REP used    !
    !                      by atoms A, B, C and D (only FUNC_BASIS instances   !
    !                      for atoms A, B and D used).                         !
    !   n_prods_hits (inout): Hit counter for the Aa-Dd NGWF product cache.    !
    !   n_prods_misses (inout): Miss counter for the Aa-Dd NGWF product cache. !
    !   calcgradient (in): Pass .true. if this is an NGWF gradient calculation.!
    !                      In the energy part, this affects the spin degeneracy!
    !                      factor. Also, it instructs this subroutine to       !
    !                      calculate grad_ket_in_ppds (see below).             !
    !   grad_ket_in_ppds (out): If calcgradient is .true., we allocate and     !
    !                           return the gradient ket k_{s,b,i_PPD;B} here.  !
    !--------------------------------------------------------------------------!
    ! Notes:                                                                   !
    !   The main part is OMP-parallelised over PPDs. It seems OMP reduction of !
    !   arrays is flaky (I got occasional zeroes in the result despite avoiding!
    !   the reduction-races-against-initialisation trap). Instead I'm reducing !
    !   manually -- each OMP thread has its own stripe, and I sum over threads !
    !   in the end.                                                            !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   If calcgradient is .true., we allocate grad_ket_in_ppds here. Caller   !
    !   is responsible for freeing it later.                                   !
    !                                                                          !
    !   Other MPI ranks may own X matrix stripes that overlap with the stripe  !
    !   on current proc -- this happens if two or more procs get B-C pairs that!
    !   share a B. Care must be taken to sum over MPI ranks. This happens in   !
    !   hfx_comms_matrix().                                                    !
    !                                                                          !
    !   Extra care must be taken to avoid:                                     !
    !     - D's that share a PPD with an A, but A,D do not sparse-intersect    !
    !       ("D cycling"),                                                     !
    !     - A's that share this PPD, but A,D do not sparse-intesect            !
    !       ("A cycling").                                                     !
    !     - A's that share this PPD, do sparse-intersect with some D, but are  !
    !       not X-neighbours with current B. Another way to put it:            !
    !       We go over PPDs that are relevant to current B, meaning there is   !
    !       some A that is an X-neighbour of B that has this PPD. Then we go   !
    !       over a list of A's sharing this PPD, but not all of these A's must !
    !       be X-neighbours of B. Some of them will be outside the X cutoff.   !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2019.                                 !
    ! Revised by Jacek Dziedzic in April 2020 to move further gradient         !
    ! operations elsewhere.                                                    !
    !==========================================================================!

    use constants, only: PI
    use hash_table, only: hash_table_lookup_ptr_nocount, hash_table_axpy
    use remote, only: remote_ngwf_cache_hts
    use rundat, only: pub_num_spins, pub_threads_max, pub_debug
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_abort, utils_alloc_check, &
         utils_dealloc_check, utils_sanity_check
    use xc, only: pub_hfxfraction

    implicit none

    ! jd: Arguments
    type(HT_HASH_TABLE), intent(inout), target :: x_ab_contributions_ht
    type(HFX_STATE), intent(in), target        :: hfxstate
    real(kind=DP), intent(in)                  :: f(:,:,:)
    integer, intent(in)                        :: a_ppd_list(:)
    integer, intent(in)                        :: n_a_ppd_list
    integer, intent(in)                        :: global_b
    type(CELL_INFO), intent(in)                :: cell
    integer, intent(in)                        :: global_nat
    real(kind=DP), intent(in)                  :: fftbox_weight
    type(HFX_NGWF_INDEX), intent(in)           :: ngwf_indices(A_NGWF:D_NGWF)
    real(kind=DP), intent(inout)               :: n_prods_hits
    real(kind=DP), intent(inout)               :: n_prods_misses
    logical, intent(in)                        :: calcgradient
    real(kind=DP), intent(out), optional, allocatable :: grad_ket_in_ppds(:,:,:,:)
    ! Indices are: pts in PPD, spins, b, PPD. PPDs go over *all* A's relevant
    !              to current B.

    ! jd: Local variables
    integer                    :: alist(max_num_at_per_ppd)
    integer                    :: dlist(max_num_at_per_ppd)
    real(kind=DP), pointer     :: alist_real_ptr(:)
    real(kind=DP), pointer     :: dlist_real_ptr(:)
    real(kind=DP), pointer     :: ppd_ad_product(:)
    integer                    :: dactual(n_a_ppd_list, max_num_at_per_ppd)
    integer                    :: ndactual(n_a_ppd_list)
    integer                    :: offsets_to_ppds(n_a_ppd_list)
    logical                    :: touched_a_atoms(global_nat,0:pub_threads_max-1)
    real(kind=DP), allocatable :: xstripe(:,:,:,:)  ! is, global_Aa, ngwf_b, tid
    real(kind=DP)              :: data_to_add(&
         ngwf_indices(A_NGWF)%basis%max_on_atom*&
         ngwf_indices(B_NGWF)%basis%max_on_atom*pub_num_spins+3)
    real(kind=DP), target      :: ppd_ad_product_workspace(cell%n_pts * &
         ngwf_indices(A_NGWF)%basis%max_on_atom * &
         ngwf_indices(D_NGWF)%basis%max_on_atom)
    real(kind=DP), pointer     :: ppd_ngwfs_dd(:) ! fused: d,pt
    real(kind=DP) :: common_factor
    real(kind=DP) :: spin_grad_fac
    integer :: fused_index_range
    integer :: num_a, num_d
    integer :: i_d_actual
    integer :: global_a, global_d
    integer :: a_idx, d_idx
    integer :: is
    integer :: n_ngwfs_a, n_ngwfs_b, n_ngwfs_d
    integer :: first_ngwf_idx_of_a, first_ngwf_idx_of_b, first_ngwf_idx_of_d
    integer :: global_aa_ngwf_idx
    integer :: ngwf_a, ngwf_b, ngwf_d
    integer :: n_cycled_a
    integer :: i_ppd
    integer :: cur_ppd
    integer :: npts
    integer :: offset, offset_save_atom_a, offset_save_a, offset_save_b
    integer :: prod_offset, prod_offset_save_b
    integer :: ngwf_dd_offset
    integer :: size_of_cached_data
    integer :: nelements ! = n_ngwfs_a * n_ngwfs_b * spins
    integer :: ierr
    integer :: sample
    integer :: tid
!$  integer, external :: omp_get_thread_num
    integer, parameter :: sample_period = 100
    character(len=*), parameter :: myself = 'hfx_all_ad_for_this_b'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    call utils_assert(hfxstate%xs_atoms_nl%populated, &
         myself//': xs_atoms neighbour list not populated.')

    first_ngwf_idx_of_b = ngwf_indices(B_NGWF)%basis%first_on_atom(global_b)
    n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)

    npts = cell%n_pts
    common_factor = 1.0_DP / (4.0_DP*PI) * fftbox_weight
    touched_a_atoms(:,:) = .false.

    call timer_clock(myself//'_main',1)

    ! jd: Allocate and zero a stripe that will be axpied into X
    allocate(xstripe(pub_num_spins,ngwf_indices(A_NGWF)%basis%num, &
         ngwf_indices(B_NGWF)%basis%max_on_atom,0:pub_threads_max-1),stat=ierr)
    call utils_alloc_check(myself,'xstripe',ierr)
    xstripe(:,:,:,:) = 0.0_DP

    fused_index_range = 0
    offsets_to_ppds(:) = garbage_int
    ! --------------------------------------------------------------------------
    ! jd: All my A PPD *relevant to current B*                             r PPD
    ! --------------------------------------------------------------------------
    loop_relev_PPD:                                                            &
    do i_ppd = 1, n_a_ppd_list
       cur_ppd = a_ppd_list(i_ppd)

       ! jd: Retrieve list of all D's sharing this PPD
       call hash_table_lookup_ptr_nocount(dlist_real_ptr, num_d, &
            hfxstate%dlists_ht, global_b, cur_ppd)
       call utils_assert(num_d <= max_num_at_per_ppd, myself//&
            ': Insufficient max_num_at_per_ppd [1]', &
            max_num_at_per_ppd, num_d)
       call utils_assert(num_d /= -1, myself//': dlist not found [1]', &
            global_b, cur_ppd)
       call utils_assert(num_d /= 0, myself//': dlist empty [1]', &
            global_b, cur_ppd)

       ! jd: Convert to ints, pad
       dlist(1:num_d) = nint(dlist_real_ptr(1:num_d))
       dlist(num_d+1:) = garbage_int

       ! -----------------------------------------------------------------------
       ! jd: For all D's that are relevant to this PPD                     r DDD
       ! -----------------------------------------------------------------------
       i_d_actual = 1
       dactual(i_ppd,:) = -1
       loop_relev_D:                                                           &
       do d_idx = 1, num_d

          global_d = dlist(d_idx)

          ! jd: Careful. There are corner cases where D shares a PPD with
          !     an A that is an X-neighbour of B, but A and D *do not*
          !     intersect (from the sparsity point of view). For every
          !     candidate D we must double-check if it's an x-S-neighbour of B.
          if(.not. any(hfxstate%xs_atoms_nl%neighbours(&
               hfxstate%xs_atoms_nl%first_idx(global_b):&
               hfxstate%xs_atoms_nl%last_idx(global_b)) == global_d)) then
             cycle
          end if

          n_ngwfs_d = ngwf_indices(D_NGWF)%basis%num_on_atom(global_d)
          ! jd: Store offset to PPD when we hit the 1st non-cycled-over D
          !     Can't do this for d_idx == 1, because this D can get cycled over
          if(offsets_to_ppds(i_ppd) == garbage_int) then
             offsets_to_ppds(i_ppd) = fused_index_range + 1
          end if
          fused_index_range = fused_index_range + n_ngwfs_d * npts
          dactual(i_ppd,i_d_actual) = global_d
          i_d_actual = i_d_actual + 1
       end do loop_relev_D
       ndactual(i_ppd) = i_d_actual - 1

    end do loop_relev_PPD

    if(calcgradient) then
       ! jd: ~20MB (443 atoms, 8 ranks)
       allocate(grad_ket_in_ppds(npts,pub_num_spins,n_ngwfs_b,n_a_ppd_list), &
            stat=ierr)
       call utils_alloc_check(myself,'grad_ket_in_ppds',ierr)
       grad_ket_in_ppds(:,:,:,:) = 0.0_DP
    end if

    ! --------------------------------------------------------------------------
    ! jd: All my A PPD *relevant to current B*                             r PPD
    ! --------------------------------------------------------------------------
    sample=0
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(DYNAMIC) &
!$OMP PRIVATE(i_ppd, cur_ppd, first_ngwf_idx_of_d, global_d, n_ngwfs_d,        &
!$OMP      tid, offset_save_b, offset_save_a, offset_save_atom_a, num_a,       &
!$OMP      alist_real_ptr, global_aa_ngwf_idx, size_of_cached_data, alist,     &
!$OMP      prod_offset, prod_offset_save_b, first_ngwf_idx_of_a, global_a,     &
!$OMP      n_ngwfs_a, ppd_ad_product, n_cycled_a, ngwf_dd_offset, ppd_ngwfs_dd,&
!$OMP      ppd_ad_product_workspace)                                           &
!$OMP FIRSTPRIVATE(npts, n_ngwfs_b, sample)                                    &
!$OMP LASTPRIVATE(offset)                                                      &
!$OMP SHARED(hfxstate, a_ppd_list, n_a_ppd_list, ngwf_indices, cell, xstripe,  &
!$OMP      ndactual, dactual, pub_num_spins, offsets_to_ppds, global_b, f,     &
!$OMP      calcgradient, grad_ket_in_ppds, touched_a_atoms,                    &
!$OMP      remote_ngwf_cache_hts)                                              &
!$OMP REDUCTION(+:n_prods_hits,n_prods_misses)

    loop_relev_PPD2:                                                           &
    do i_ppd = 1, n_a_ppd_list
       cur_ppd = a_ppd_list(i_ppd)

       offset = offsets_to_ppds(i_ppd)

       tid = 0
!$     tid = omp_get_thread_num()

       ! -----------------------------------------------------------------------
       ! jd: For all D's that are relevant to this PPD                     r DDD
       !     ^ with D's that share a PPD with A, but are not s-X-neighbours
       !       of current B eliminated.
       ! -----------------------------------------------------------------------
       loop_actual_D2:                                                         &
       do d_idx = 1, ndactual(i_ppd)

          global_d = dactual(i_ppd,d_idx)

          first_ngwf_idx_of_d = &
              ngwf_indices(D_NGWF)%basis%first_on_atom(global_d)
          n_ngwfs_d = ngwf_indices(D_NGWF)%basis%num_on_atom(global_d)

          ! jd: Get a list of all my A's that share this PPD
          call hash_table_lookup_ptr_nocount(alist_real_ptr, num_a, &
               hfxstate%alists_ht, cur_ppd)
          if(num_a > max_num_at_per_ppd) then
             call utils_abort(myself//': Insufficient max_num_a_per_ppd',&
                  max_num_at_per_ppd, num_a)
          end if
          if(num_a == -1) then
             call utils_abort(myself//': alist not found', cur_ppd)
          end if
          if(num_a == 0) then
             call utils_abort(myself//': alist empty', cur_ppd)
          end if
          ! jd: Convert to ints, pad
          alist(1:num_a) = nint(alist_real_ptr(1:num_a))
          alist(num_a+1:) = garbage_int

          ! --------------------------------------------------------------------
          ! jd: For all A's that are relevant to this PPD                    AAA
          ! --------------------------------------------------------------------
          n_cycled_a = 0
          loop_relev_A:                                                        &
          do a_idx = 1, num_a
             global_a = alist(a_idx)
             n_ngwfs_a = ngwf_indices(A_NGWF)%basis%num_on_atom(global_a)
             first_ngwf_idx_of_a = &
                 ngwf_indices(A_NGWF)%basis%first_on_atom(global_a)

             if(a_idx == 1) then
                offset_save_atom_a = offset
             end if

             if(hfxstate%ppd_products_ht_in_use) then
                ! jd: Retrieve PPD products for all {a,d} for this PPD
                call hash_table_lookup_ptr_nocount(ppd_ad_product, &
                     size_of_cached_data, hfxstate%ppd_products_ht, &
                     global_a, global_d, cur_ppd)
             else
                size_of_cached_data = -1
             end if

             ! jd: If a PPD product is absent, this is either:
             !     - a cache miss (it did not fit in the HT),
             !     - a corner case, where A and D do share a PPD, but
             !       do not s-overlap,
             !     - a bug.
             if(size_of_cached_data == -1) then
                if(any(hfxstate%s_atoms_nl%neighbours(&
                     hfxstate%s_atoms_nl%first_idx(global_a):&
                     hfxstate%s_atoms_nl%last_idx(global_a)) == global_d)) then
                   ! jd: They do s-overlap, but the product is missing.
                   !     We assume this means a cache miss. Bugs are unlikely,
                   !     since this has been tested extensively in a previous
                   !     scheme where *all* products were cached, which exposed
                   !     bugs easily.
                   ! jd: Calculate the product for this PPD on the fly.
                   call hfx_calc_aa_dd_ngwf_prod_in_ppd(&
                        ppd_ad_product_workspace, remote_ngwf_cache_hts, cell, &
                        ngwf_indices, global_a, global_d, cur_ppd, sample, &
                        sample_period)
                   ppd_ad_product => ppd_ad_product_workspace
                   size_of_cached_data = npts * n_ngwfs_a * n_ngwfs_d
                   sample = sample + 1
                   n_prods_misses = n_prods_misses + 1.0_DP
                else
                   ! jd: They don't s-overlap, meaning corner case. Skip this A.
                   n_cycled_a = n_cycled_a + 1
                   cycle
                end if
             else
                n_prods_hits = n_prods_hits + 1.0_DP
             end if
             if(size_of_cached_data /= npts * n_ngwfs_a * n_ngwfs_d) then
                call utils_abort(myself//': Unexpected ndata for cached &
                     &PPD product.', size_of_cached_data, npts, &
                     n_ngwfs_a, n_ngwfs_d)
             end if

             ! jd: Another corner case where we're in a PPD that is relevant
             !     to current B (so *some* A that is X-neighbours with B has
             !     it), and there's another A (us) that shares this PPD, but
             !     even though the two A's share a D (so we did not get cycled
             !     above), we are not an X-neighbour of B.
             if (.not. any(hfxstate%x_atoms_nl%neighbours(&
                  hfxstate%x_atoms_nl%first_idx(global_b):&
                  hfxstate%x_atoms_nl%last_idx(global_b)) == global_a)) then
                n_cycled_a = n_cycled_a + 1
                cycle
             end if

             touched_a_atoms(global_a,tid) = .true.

             ! jd: Crucial that this comes *after* the above potential cycles
             if (a_idx > 1) then
                offset = offset_save_atom_a
             end if

             prod_offset = 1

             ! -----------------------------------------------------------------
             ! jd: For all d                                                 ddd
             ! -----------------------------------------------------------------
             ngwf_dd_offset = 1
             loop_ngwf_d:                                                      &
             do ngwf_d = 1, n_ngwfs_d

                ! jd: Retrieve NGWFs Dd for all d for this PPD if doing gradient
                !     Do this only for the 1st A to save effort.
                if(calcgradient .and. (a_idx-n_cycled_a == 1)) then
                   call hash_table_lookup_ptr_nocount(ppd_ngwfs_dd, &
                        size_of_cached_data, hfxstate%ppd_ngwfs_dd_ht, &
                        global_d, cur_ppd)
                   if(size_of_cached_data /= npts * n_ngwfs_d) then
                      call utils_abort(myself//': Unexpected ndata for &
                           &cached NGWF Dd.', size_of_cached_data, npts, &
                           n_ngwfs_d)
                   end if
                end if

                ! --------------------------------------------------------------
                ! jd: For all b                                              bbb
                ! --------------------------------------------------------------
                loop_ngwf_b:                                                   &
                do ngwf_b = 1, n_ngwfs_b

                   if(ngwf_b == 1) then
                      offset_save_b = offset
                      prod_offset_save_b = prod_offset
                   else
                      offset = offset_save_b
                      prod_offset = prod_offset_save_b
                   end if

                   ! -----------------------------------------------------------
                   ! jd: For all a                                           aaa
                   ! -----------------------------------------------------------
                   loop_ngwf_a:                                                &
                   do ngwf_a = 1, n_ngwfs_a
                      global_aa_ngwf_idx = first_ngwf_idx_of_a + ngwf_a - 1

                      if(ngwf_a == 1) then
                         offset_save_a = offset
                      else
                         offset = offset_save_a
                      end if

                      ! --------------------------------------------------------
                      ! jd: For all spins                                    sss
                      ! --------------------------------------------------------
                      loop_spins:                                              &
                      do is = 1, pub_num_spins

                         if(pub_debug) then
                            call utils_sanity_check(xstripe(is, &
                                 global_aa_ngwf_idx,ngwf_b,tid), 'stripe', &
                                 excessive=1D10)
                            call utils_sanity_check(f(offset:offset+npts-1,&
                                 is,ngwf_b), 'f', excessive=1D10)
                            call utils_sanity_check(ppd_ad_product(&
                                 prod_offset:prod_offset+npts-1), 'product', &
                                 excessive=1D10)
                         end if

                         ! -----------------------------------------------------
                         xstripe(is,global_aa_ngwf_idx,ngwf_b,tid) = &
                              xstripe(is,global_aa_ngwf_idx,ngwf_b,tid) + &
                              sum(f(offset:offset+npts-1,is,ngwf_b) * &
                              ppd_ad_product(prod_offset:prod_offset+npts-1))
                         ! -----------------------------------------------------

                         ! jd: Calculate gradient term 1 kets when necessary.
                         !     Do not recalculate it when only A or a change.
                         if(calcgradient .and. (a_idx-n_cycled_a == 1) .and. &
                              ngwf_a == 1) then
                            grad_ket_in_ppds(1:npts,is,ngwf_b,i_ppd) = &
                                 grad_ket_in_ppds(1:npts,is,ngwf_b,i_ppd) + &
                                 f(offset:offset+npts-1,is,ngwf_b) * &
                                 ppd_ngwfs_dd(ngwf_dd_offset:&
                                 ngwf_dd_offset+npts-1)

                            if(pub_debug) then
                               call utils_sanity_check(f(&
                                    offset:offset+npts-1,is,ngwf_b),'f')
                               call utils_sanity_check(ppd_ngwfs_dd(&
                                    ngwf_dd_offset:ngwf_dd_offset+npts-1),&
                                    'ppd_ngwfs_dd')
                               call utils_sanity_check(grad_ket_in_ppds(&
                                    1:npts,is,ngwf_b,i_ppd),'grad_ket_in_ppds')
                            end if

                         end if

                      end do loop_spins

                      offset = offset + npts
                      prod_offset = prod_offset + npts

                   end do loop_ngwf_a

                end do loop_ngwf_b

                if(calcgradient) ngwf_dd_offset = ngwf_dd_offset + npts

             end do loop_ngwf_d

          end do loop_relev_A

          ! jd: If all A's were skipped for this D, need to fast-forward offset
          if(n_cycled_a == num_a) then
             offset = offset + npts * n_ngwfs_d
          end if

       end do loop_actual_D2

       if(i_ppd < n_a_ppd_list) then
          if(offset /= offsets_to_ppds(i_ppd+1)) then
             call utils_abort(myself//': Book-keeping error reading f, thread',&
                  tid, offset, offsets_to_ppds(i_ppd+1), avoid_mpi_calls = .true.)
          end if
       end if

    end do loop_relev_PPD2
!$OMP END PARALLEL DO

    call utils_assert(offset-1 == ubound(f,1), &
         myself//': Book-keeping error reading f.', offset-1, ubound(f,1))

    call timer_clock(myself//'_main',2)
    call timer_clock(myself//'_store',1)

    ! **************************************************************************
    ! * 5. Scale and store results
    ! **************************************************************************

    ! qoh: Scale by spin factor
    if (calcgradient) then
       spin_grad_fac = 2.0_DP
    else
       spin_grad_fac = real(pub_num_spins,kind=DP)
    end if

    ! qoh: Scale by 0.5 (for Hartree-2-exchange ratio) and pub_hfxfraction
    xstripe = xstripe * 0.5_DP * pub_hfxfraction * spin_grad_fac * common_factor

    ! --------------------------------------------------------------------------
    ! jd: For all A that were of interest above                              AAA
    ! --------------------------------------------------------------------------
    loop_A2:                                                                   &
    do global_a = 1, global_nat

       if(.not. any(touched_a_atoms(global_a,:))) cycle

       n_ngwfs_a = ngwf_indices(A_NGWF)%basis%num_on_atom(global_a)
       first_ngwf_idx_of_a = ngwf_indices(A_NGWF)%basis%first_on_atom(global_a)

       ! jd: xstripe is dimensioned for max_on_atom, so must extract
       !     the relevant part. Also we want to store global_[ab] there.
       !     We also need to reduce across threads (last index).
       nelements = n_ngwfs_a * n_ngwfs_b * pub_num_spins
       offset = 4
       do is = 1, pub_num_spins
          do ngwf_a = 1, n_ngwfs_a
             global_aa_ngwf_idx = first_ngwf_idx_of_a + ngwf_a - 1
             do ngwf_b = 1, n_ngwfs_b
                data_to_add(offset) = sum(xstripe(is,global_aa_ngwf_idx,ngwf_b,:))
                offset = offset + 1
             end do
          end do
       end do

       ! jd: Add metadata: global_a, global_b, nelements. This makes our job
       !     easier during comms later.
       data_to_add(1) = real(global_a,kind=DP)
       data_to_add(2) = real(global_b,kind=DP)
       data_to_add(3) = real(nelements,kind=DP)

       call hash_table_axpy(x_ab_contributions_ht, data_to_add, &
            nelements+3, global_a, global_b, overfill_strategy = 'F', &
            dont_axpy_first_nelems = 3)

    end do loop_A2

    deallocate(xstripe,stat=ierr)
    call utils_dealloc_check(myself,'xstripe',ierr)

    call timer_clock(myself//'_store',2)

    call timer_clock(myself,2)

  end subroutine hfx_all_ad_for_this_b

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_free_arrays(hfxstate)
    !==========================================================================!
    ! Frees the allocatables in hfxstate.                                      !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2019.                              !
    !==========================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target :: hfxstate

    ! jd: Local variables
    integer :: ierr
    character(len=*), parameter :: myself = 'hfx_free_arrays'

    ! -------------------------------------------------------------------------

    ! jd: NB: B and C atom lists are allocated via allocate in this module.
    !         A and D atom lists are returned by neighbour_list_for_atom_list().
    !         local B and local D atom lists are returned (indirectly) by
    !         neighbour_list_for_atom_list(). Hence the different treatment below.

    deallocate(hfxstate%my_bc_pairs,stat=ierr)
    call utils_dealloc_check(myself,'my_bc_pairs',ierr)
    deallocate(hfxstate%my_a_atoms,stat=ierr)
    call utils_dealloc_check(myself,'out_atom_list',ierr, &
         allocated_in='neighbour_list_for_atom_list')
    deallocate(hfxstate%my_bb_ngwf_offsets,stat=ierr)
    call utils_dealloc_check(myself,'my_bb_ngwf_offsets',ierr)
    deallocate(hfxstate%my_b_atoms,stat=ierr)
    call utils_dealloc_check(myself,'my_b_atoms',ierr)
    deallocate(hfxstate%my_c_atoms,stat=ierr)
    call utils_dealloc_check(myself,'my_c_atoms',ierr)
    deallocate(hfxstate%my_d_atoms,stat=ierr)
    call utils_dealloc_check(myself,'out_atom_list',ierr, &
         allocated_in='neighbour_list_for_atom_list')
    deallocate(hfxstate%my_a_ppd_list,stat=ierr)
    call utils_dealloc_check(myself,'my_a_ppd_list',ierr)

  end subroutine hfx_free_arrays

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_symmetrise_xmatrix(hfexchange)
    !==========================================================================!
    ! Symmetrises the exchange matrix.                                         !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2012 basing on Quintin Hill's code.    !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: VERBOSE, stdout
    use rundat, only: pub_num_spins, pub_spin_fac, pub_hfx_output_detail
    use sparse, only: sparse_axpy, sparse_create, sparse_copy, sparse_destroy, &
         sparse_rms_element, sparse_scale, sparse_transpose, SPAM3
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(SPAM3), intent(inout) :: hfexchange(pub_num_spins)

    ! jd: Local variables
    type(SPAM3)   :: transhfx  ! Transposed HF exchange matrix
    type(SPAM3)   :: hfx_temp  ! Temporary matrix to calculate asymmetry
    real(kind=DP) :: rms_diff
    integer       :: is

    ! ------------------------------------------------------------------------

    ! qoh: Make hfexchange matrix symmetric
    ! jd:  ... and calculate the rms asymmetry
    call sparse_create(transhfx,hfexchange(1))
    call sparse_create(hfx_temp,hfexchange(1))
    rms_diff = 0.0_DP
    do is = 1, pub_num_spins

       call sparse_copy(hfx_temp,hfexchange(is))
       call sparse_scale(hfexchange(is),0.5_DP)
       call sparse_transpose(transhfx,hfexchange(is))
       call sparse_axpy(hfexchange(is),transhfx,1.0_DP)
       call sparse_axpy(hfx_temp,transhfx,-2.0_DP) ! (already scaled by 0.5)
       rms_diff = rms_diff + sparse_rms_element(hfx_temp) * 0.5_DP * pub_spin_fac
    end do
    call sparse_destroy(transhfx)
    call sparse_destroy(hfx_temp)

    if(pub_on_root) then
       if(pub_hfx_output_detail >= VERBOSE) then
          if(rms_diff < 1D-16) then
             write(stdout,'(a)') &
                  'HFx: Exchange matrix symmetric to machine precision.'
          else
             write(stdout,'(a,f5.2,a)') 'HFx: Exchange matrix symmetric to ', &
                  -log10(rms_diff),' digits (rms).'
          end if
       end if
       ! jd: Don't let calculations continue if exchange matrix appears broken
       !     The threshold corresponds to 3.4 correct decimals.
       call utils_assert(rms_diff < 3.98107D-4, 'Exchange matrix not deemed &
            &accurate enough for a stable calculation',rms_diff, &
            -log10(rms_diff))
    end if

  end subroutine hfx_symmetrise_xmatrix

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_grad_ket_dkn_baccum(grad_akets_dkn_my_ht, &
       grad_ket_in_ppds, packed_dkn_ab_blocks_hts, a_ppd_list, n_a_ppd_list, &
       global_b, cell_npts, ngwf_indices, hfxstate, grad_prefactor)
    !==========================================================================!
    ! Multiplies NGWF gradient term 1 ket in PPDs by K^{Aa,Bb} or tcK^{Aa}_{Bb}!
    ! and sums over b. The ket (grad_ket_in_PPDs) is for a given B:            !
    ! k_{s,b,i_PPD;B} = \sum_{Dd} \phi_{Dd,i_PPD} * f_{bDdsi_PPD;B}, and has   !
    ! been calculated in hfx_all_ad_for_this_b().                              !
    ! B is a my-atom, nothing is sparse-local here. The PPDs are the set of    !
    ! all PPDs belonging to A atoms that are relevant to this B.               !
    !                                                                          !
    ! The result is accumulated in 'grad_aket_dkn' ('z'):                      !
    !                                                                          !
    ! z_{Aa,ppd,s,wk} = factor * \sum_B \sum_b K^{Aa,Bb,s,wk} * k_{s,b,i_ppd;B}!
    !                                                                          !
    ! This subroutine only does the sum over b. Multiple calls are issued, one !
    ! for each B, to obtain the total sum. The grad_aket_dkn's are axpied into !
    ! the HT grad_akets_dkn_my_ht. The HT only stores partial sums over my B   !
    ! These are later combined, and redistributed to the sparse-local repre-   !
    ! sentation in hfx_comms_gradient().                                       !
    ! The keys in the HT are (Aa, ppd, is).                                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   grad_akets_dkn_my_ht (inout): Results are accumulated here.            !
    !   grad_ket_in_ppds (inout): 'k' is read from here. We deallocate it here.!
    !   packed_dkn_ab_blocks_hts (in): K and tcK elements are taken from here. !
    !   a_ppd_list (in): List of PPDs in the union of all relevant A's.        !
    !   n_a_ppd_list (in): Number of elements in the above.                    !
    !   global_b (in): Global index of atom B over whose b we sum here.        !
    !   cell_npts (in): Pass cell%n_pts.                                       !
    !   ngwf_indices (in): Array of HFX_NGWF_INDEX instances containing        !
    !                      pointers to the NGWF basis and the NGWF_REP used    !
    !                      by atoms A, B, C and D (we only use A and B here).  !
    !   hfxstate (in): Needed for alists and X-neighbour list.                 !
    !   grad_prefactor (in, opt): Defaults to 1.0. In conduction 0.5 can be    !
    !                             passed to take into account the fact that    !
    !                             only one set NGWFs is optimised, and the     !
    !                             other one is fixed (I think).                !
    !                             This also takes care of spin degeneracy      !
    !                             conventions in conduction (see the call to   !
    !                             hf_exchange_calculate() from conduction).    !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in April 2020.                                 !
    !==========================================================================!

    use constants, only: garbage_int, PI
    use hash_table, only: HT_HASH_TABLE, hash_table_lookup_ptr_nocount, &
         hash_table_axpy
    use rundat, only: pub_num_spins
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_dealloc_check, utils_unit

    implicit none

    ! jd: Arguments
    type(HT_HASH_TABLE), intent(inout), target :: grad_akets_dkn_my_ht
    real(kind=DP), intent(inout), allocatable  :: grad_ket_in_ppds(:,:,:,:)
    type(HT_HASH_TABLE), intent(in), target    :: packed_dkn_ab_blocks_hts(3)
    integer, intent(in)                        :: a_ppd_list(:) ! union over relev A
    integer, intent(in)                        :: n_a_ppd_list
    integer                                    :: global_b
    integer                                    :: cell_npts
    type(HFX_NGWF_INDEX), intent(in)           :: ngwf_indices(A_NGWF:D_NGWF)
    type(HFX_STATE), intent(in)                :: hfxstate
    real(kind=DP), intent(in), optional        :: grad_prefactor

    ! jd: Local variables
    integer       :: alist(max_num_at_per_ppd)
    real(kind=DP) :: dkn_ba(&
         ngwf_indices(B_NGWF)%basis%max_on_atom, &
         ngwf_indices(A_NGWF)%basis%max_on_atom, pub_num_spins, 2)
    real(kind=DP) :: ket(cell_npts+4)
    real(kind=DP), pointer :: alist_real_ptr(:)
    real(kind=DP) :: loc_grad_prefactor
    real(kind=DP) :: factor
    integer :: i_ppd
    integer :: cur_ppd
    integer :: num_a
    integer :: a_idx
    integer :: global_a
    integer :: n_ngwfs_a, n_ngwfs_b
    integer :: ngwf_a, ngwf_b
    integer :: first_ngwf_idx_of_a, first_ngwf_idx_of_b
    integer :: global_aa, global_bb
    integer :: is
    integer :: ierr
    logical :: first_time
    character(len=*), parameter :: myself = 'hfx_grad_ket_dkn_baccum'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    ! jd: Say hi, but only for the first B on root, which we recognize by the
    !     HT being empty.
    first_time = (grad_akets_dkn_my_ht%n_slots + &
         grad_akets_dkn_my_ht%n_chained == 0)

    n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)
    first_ngwf_idx_of_b = ngwf_indices(B_NGWF)%basis%first_on_atom(global_b)

    ! jd: Honestly have no idea of the reason for the 4 PI here.
    !     It came about when PPDs were introduced instead of TBs.
    !     The 1/2 is due to the fact that the bra-ket symmetry is broken
    !     (kets are expanded, bras are not) and the two terms in the
    !     NGWF gradient are calculated separately.
    if(present(grad_prefactor)) then
       loc_grad_prefactor = grad_prefactor
    else
       loc_grad_prefactor = 1.0_DP
    end if
    ! jd: loc_grad_prefactor allows us to scale the NGWF gradient by a factor
    !     of 0.5 in conduction to account for the fact that half of the NGWFs
    !     (the valence ones) are fixed.

    factor = 0.5_DP / (4.0_DP * PI) * loc_grad_prefactor

    ! --------------------------------------------------------------------------
    ! jd: All my A PPD *relevant to current B*                             r PPD
    ! --------------------------------------------------------------------------
    loop_relev_PPD:                                                            &
    do i_ppd = 1, n_a_ppd_list
       cur_ppd = a_ppd_list(i_ppd)

       ket(2) = real(cur_ppd,kind=DP)

       ! jd: Get a list of all my A's that share this PPD
       call hash_table_lookup_ptr_nocount(alist_real_ptr, num_a, &
            hfxstate%alists_ht, cur_ppd)
       if(num_a > max_num_at_per_ppd) then
          call utils_abort(myself//': Insufficient max_num_a_per_ppd [2]',&
               max_num_at_per_ppd, num_a)
       end if
       if(num_a == -1) then
          call utils_abort(myself//': alist not found [2]', cur_ppd)
       end if
       if(num_a == 0) then
          call utils_abort(myself//': alist empty [2]', cur_ppd)
       end if
       ! jd: Convert to ints, pad
       alist(1:num_a) = nint(alist_real_ptr(1:num_a))
       alist(num_a+1:) = garbage_int

       ! -----------------------------------------------------------------------
       ! jd: For all A's that are relevant to this PPD                       AAA
       ! -----------------------------------------------------------------------
       loop_relev_A:                                                           &
       do a_idx = 1, num_a
          global_a = alist(a_idx)
          n_ngwfs_a = ngwf_indices(A_NGWF)%basis%num_on_atom(global_a)
          first_ngwf_idx_of_a = &
               ngwf_indices(A_NGWF)%basis%first_on_atom(global_a)

          ! jd: Another corner case where we're in a PPD that is relevant
          !     to current B (so *some* A that is X-neighbours with B has
          !     it), and there's another A (us) that shares this PPD, but
          !     even though the two A's share a D (so we did not get cycled
          !     above), we are not an X-neighbour of B.
          if (.not. any(hfxstate%x_atoms_nl%neighbours(&
               hfxstate%x_atoms_nl%first_idx(global_b):&
               hfxstate%x_atoms_nl%last_idx(global_b)) == global_a)) then
             cycle
          end if

          ! jd: Get K^{BA} and tcK^B_A (sic!) (all b's, all a's, spins).
          !     This not very efficient, since we're inside the PPD loop.
          !     For K we rely on symmetry to re-use the dkn_ab HT.
          !     For tcK, since it's not symmetric, we have a separate HT.
          call hfx_get_dkn_atomblock(dkn_ba(:,:,:,1), & ! out
               packed_dkn_ab_blocks_hts(1), &
               global_a, global_b, &
               ngwf_indices(A_NGWF)%basis%max_on_atom, &
               ngwf_indices(B_NGWF)%basis%max_on_atom)

          call hfx_get_dkn_atomblock(dkn_ba(:,:,:,2), & ! out
               packed_dkn_ab_blocks_hts(3), &
               global_b, global_a, &
               ngwf_indices(B_NGWF)%basis%max_on_atom, &
               ngwf_indices(A_NGWF)%basis%max_on_atom)

          dkn_ba(:,:,:,2) = -dkn_ba(:,:,:,2)
                            ! ^ equivalent to using -tcK in old code

          ! --------------------------------------------------------------------
          ! jd: For all b                                                    bbb
          ! --------------------------------------------------------------------
          loop_ngwf_b:                                                         &
          do ngwf_b = 1, n_ngwfs_b
             global_bb = first_ngwf_idx_of_b + ngwf_b - 1

             ! -----------------------------------------------------------------
             ! jd: For all spins                                             sss
             ! -----------------------------------------------------------------
             loop_spins:                                                       &
             do is = 1, pub_num_spins

                ket(3) = real(is,kind=DP)

                ! --------------------------------------------------------------
                ! jd: For all a                                              aaa
                ! --------------------------------------------------------------
                loop_ngwf_a:                                                   &
                do ngwf_a = 1, n_ngwfs_a
                   global_aa = first_ngwf_idx_of_a + ngwf_a - 1

                   ket(1) = real(global_aa,kind=DP)

#if 0
    ! jd: Useful when debugging
    do which_kernel = 1, 2
       write(filename,'(a,i0,a,i0,a,i0,a,i0)') &
            'k_', global_aa, '_', global_bb, '_', cur_ppd, '_', which_kernel
       u = utils_unit()
       open(u, file=filename)
       write(u,*) dkn_ba(ngwf_b,ngwf_a,1,which_kernel)
       close(u)
    end do
#endif

                   ket(4) = real(1.0_DP,kind=DP)
                   ket(5:5+cell_npts-1) = &
                        grad_ket_in_ppds(:,is,ngwf_b,i_ppd) * &
                        dkn_ba(ngwf_a,ngwf_b,is,1) * factor
                        ! ^ uses symmetry

                   ! jd: Axpy result into output HT
                   call hash_table_axpy(grad_akets_dkn_my_ht, &
                        ket(:), cell_npts+4, &
                        global_aa, cur_ppd, is, 1, &
                        overfill_strategy = 'F', dont_axpy_first_nelems = 4)

                   ket(4) = real(2.0_DP,kind=DP)
                   ket(5:5+cell_npts-1) = &
                        grad_ket_in_ppds(:,is,ngwf_b,i_ppd) * &
                        dkn_ba(ngwf_b,ngwf_a,is,2) * factor

                   ! jd: Axpy result into output HT
                   call hash_table_axpy(grad_akets_dkn_my_ht, &
                        ket(:), cell_npts+4, &
                        global_aa, cur_ppd, is, 2, &
                        overfill_strategy = 'F', dont_axpy_first_nelems = 4)


                end do loop_ngwf_a

             end do loop_spins

          end do loop_ngwf_b

       end do loop_relev_A

    end do loop_relev_PPD

    deallocate(grad_ket_in_ppds, stat=ierr)
    call utils_dealloc_check(myself, 'grad_ket_in_ppds', ierr, &
         allocated_in = 'hfx_all_ad_for_this_b')

    call timer_clock(myself,2)

  end subroutine hfx_grad_ket_dkn_baccum

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_gradient(cov_grad, contra_grad, grad_akets_dkn_loc_ht,&
       precond_func_recip, cell, fftbox, par, alpha_ngwf_index)
    !==========================================================================!
    ! Calculates the HFx contribution to the NGWF gradient for all NGWFs a on  !
    ! all sparse-local atoms A (so we use the usual ONETEP distribution) from  !
    ! term 1. Reciprocal space preconditioning is applied if necessary.        !
    ! HFx contributions are *added* to cov_grad and contra_grad.               !
    !--------------------------------------------------------------------------!
    ! Contains vestiges of the subroutine:                                     !
    !   - originally written by Quintin Hill in late 2008 and January 2009,    !
    !   - adapted by Jacek Dziedzic to new O(N) HF exchange in July 2012,      !
    !   - reworked by Jacek Dziedzic for PPDs in April 2015,                   !
    !   - extended by Jacek Dziedzic to handle mixed NGWF bases in June 2018.  !
    ! Rehauled by Jacek Dziedzic for new paralellisation scheme in May 2019.   !
    !==========================================================================!

    use basis, only: basis_clean_function, basis_copy_function_to_box, &
         basis_extract_function_from_box
    use comms, only: pub_my_proc_id, comms_barrier
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero, &
         data_functions_alloc, data_functions_dealloc, data_fftbox_alloc, &
         data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box_pair
    use function_ops, only: function_ops_batch_col_start
    use hash_table, only: HT_HASH_TABLE, hash_table_lookup_ptr_nocount
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_spin_fac, pub_num_spins, pub_precond_recip, &
         pub_hhf_nstates, pub_hhf_factor, pub_debug
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_abort, &
         utils_alloc_check, utils_dealloc_check, utils_sanity_check, &
         utils_int_to_str, utils_assert
    use xc, only: pub_hfxfraction

    implicit none

    ! jd: Arguments
    type(FUNCTIONS), intent(inout)  :: cov_grad
    type(FUNCTIONS), intent(inout)  :: contra_grad
    type(HT_HASH_TABLE), intent(in) :: grad_akets_dkn_loc_ht
    real(kind=DP), intent(in), optional :: precond_func_recip(:,:,:)
    type(CELL_INFO), intent(in)     :: cell
    type(FFTBOX_INFO), intent(in)   :: fftbox
    type(PARAL_INFO), intent(in)    :: par
    type(HFX_NGWF_INDEX), intent(in) :: alpha_ngwf_index

    ! jd: Local variables
    real(kind=DP),    allocatable  :: work_fftbox(:,:,:)
    complex(kind=DP), allocatable  :: zwork_fftbox(:,:,:)
    type(FFTBOX_DATA), allocatable :: kets_in_fftbox(:,:)
    type(FUNCTIONS), allocatable   :: kets_in_ppds_on_a(:,:,:)
    type(FUNCTIONS)                :: aket_dkn_in_ppds(2)
    real(kind=DP), pointer         :: aket_dkn_in_ppd_ptr(:)
    integer :: fa_box_start(3,1)
    integer :: fa_start_in_box(3,1)
    integer :: global_a
    integer :: local_a
    integer :: n_ngwfs_a
    integer :: ngwf_a
    integer :: local_ngwf_aa
    integer :: global_ngwf_aa
    integer :: i_ppd
    integer :: cur_ppd
    integer :: n_ppds
    integer :: ndata
    integer :: ipt
    integer :: dot_npts
    integer :: bra_start
    integer :: is
    integer :: which_kernel
    integer :: offset
    integer :: ierr
    real(kind=DP) :: factor
    character(len=*), parameter :: myself = 'hfx_gradient'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

#if 0
    print *, "-NGWF_CONTRA_TOT_PRE:-============================"
    print *, contra_grad%d

    print *, "-NGWF_COV_TOT_PRE:-============================"
    print *, cov_grad%d
#endif

    ! qoh: Apply scaling for exchange to all density kernel elements
    ! qoh: Scaling is 4 (for the 4 NGWFs) * 0.5 (Hartree:Exchange ratio)
    ! qoh: * 0.5 (double counting of electron pairs) * (2/num_spins)^2
    ! qoh: (scale both density kernels by spin factor) * hfxfraction *
    ! qoh: grid point weight * num_spins
    factor = pub_hfxfraction * cell%weight * 2.0_DP * pub_spin_fac

    ! cks: apply hyper Hartree-Fock factor
    if (pub_hhf_nstates > 0) then
       factor = factor * pub_hhf_factor
    endif

    ! jme: allocate buffers
    do which_kernel = 1,2
       call data_functions_alloc(aket_dkn_in_ppds(which_kernel), &
            alpha_ngwf_index%basis%max_n_ppds_sphere * cell%n_pts, &
            iscmplx = contra_grad%iscmplx) !jmecmplx
    end do

    ! --------------------------------------------------------------------------
    ! jd: Loop over local A's                                                AAA
    ! --------------------------------------------------------------------------
    loop_A:                                                                    &
    do local_a = 1, par%num_atoms_on_proc(pub_my_proc_id)
       global_a = par%first_atom_on_proc(pub_my_proc_id) + local_a - 1
       n_ngwfs_a = alpha_ngwf_index%basis%num_on_atom(global_a)

       if(pub_precond_recip) then
          allocate(kets_in_ppds_on_a(2,n_ngwfs_a,pub_num_spins), stat=ierr)
          call utils_alloc_check(myself,'kets_in_ppds_on_a',ierr)
          !jmecmplx
          do is = 1, pub_num_spins
             do ngwf_a = 1, n_ngwfs_a
                do which_kernel = 1, 2
                   call data_functions_alloc(kets_in_ppds_on_a(&
                        which_kernel,ngwf_a,is), &
                        alpha_ngwf_index%basis%max_n_ppds_sphere * cell%n_pts, &
                        iscmplx = contra_grad%iscmplx)
                end do
             end do
          end do
       end if

       ! -----------------------------------------------------------------------
       ! jd: for all spins                                                   sss
       ! -----------------------------------------------------------------------
       loop_spins:                                                             &
       do is = 1, pub_num_spins

          ! --------------------------------------------------------------------
          ! jd: for all a on A                                               aaa
          ! --------------------------------------------------------------------
          loop_ngwf_a:                                                         &
          do ngwf_a = 1, n_ngwfs_a

             global_ngwf_aa = &
                  alpha_ngwf_index%basis%first_on_atom(global_a)+ngwf_a-1
             local_ngwf_aa = global_ngwf_aa - &
                  alpha_ngwf_index%basis%first_on_proc(pub_my_proc_id) + 1
             n_ppds = alpha_ngwf_index%basis%&
                  spheres(local_ngwf_aa)%n_ppds_sphere

             ! -----------------------------------------------------------------
             ! jd: for all PPDs in Aa                                        PPD
             ! -----------------------------------------------------------------
             offset = 1
             loop_PPD:                                                         &
             do i_ppd = 1, n_ppds
                cur_ppd = alpha_ngwf_index%basis%&
                     spheres(local_ngwf_aa)%ppd_list(1,i_ppd)

                ! --------------------------------------------------------------
                ! jd: for K and tcK                                        K/tcK
                ! --------------------------------------------------------------
                loop_which_kernel:                                             &
                do which_kernel = 1, 2
                   call hash_table_lookup_ptr_nocount(aket_dkn_in_ppd_ptr, &
                        ndata, grad_akets_dkn_loc_ht, &
                        global_ngwf_aa, cur_ppd, is, which_kernel)
                   if(ndata /= cell%n_pts) then
                      call utils_abort(myself//': grad_aket_dkn not found or &
                           &size mismatch: Aa='//&
                           trim(utils_int_to_str(global_ngwf_aa))//&
                           ', PPD='//trim(utils_int_to_str(cur_ppd))//&
                           ', is='//trim(utils_int_to_str(is))//&
                           ', wk='//trim(utils_int_to_str(which_kernel)),&
                           ndata, cell%n_pts)
                   end if
                   if(pub_debug) then
                      call utils_sanity_check(aket_dkn_in_ppd_ptr(1:cell%n_pts), &
                           'aket_dkn_in_ppd_ptr')
                   end if

                   ! jd: Store in FUNCTIONS structure
                   aket_dkn_in_ppds(which_kernel)%d(offset:offset+cell%n_pts-1) = &
                        aket_dkn_in_ppd_ptr(1:cell%n_pts)

                   if(pub_precond_recip) then
                      kets_in_ppds_on_a(which_kernel,ngwf_a,is)%d(&
                           offset:offset+cell%n_pts-1) = &
                           aket_dkn_in_ppd_ptr(1:cell%n_pts)
                   end if

                end do loop_which_kernel

                offset = offset + cell%n_pts

             end do loop_PPD

             ! cks: shaving: zero points outside NGWF sphere in PPD rep.
             if(cell%n_pts > 1) then
                call basis_clean_function(aket_dkn_in_ppds(1),&
                     alpha_ngwf_index%basis%spheres(local_ngwf_aa), &
                     cell, fftbox, 1)
                call basis_clean_function(aket_dkn_in_ppds(2),&
                     alpha_ngwf_index%basis%spheres(local_ngwf_aa), &
                     cell, fftbox, 1)
             end if

             dot_npts = alpha_ngwf_index%basis%spheres(local_ngwf_aa)%&
                  n_ppds_sphere * cell%n_pts
             bra_start = alpha_ngwf_index%basis%spheres(local_ngwf_aa)%offset

             ! jd: Add contribution to {contra,cov}_grad.
             do ipt = 0, dot_npts - 1
                contra_grad%d(bra_start+ipt) = contra_grad%d(bra_start+ipt)- & !jmecmplx
                     factor * aket_dkn_in_ppds(1)%d(ipt+1)                     !jmecmplx
                if(.not. pub_precond_recip) then
                   cov_grad%d(bra_start+ipt) = cov_grad%d(bra_start+ipt)- &    !jmecmplx
                        factor * aket_dkn_in_ppds(2)%d(ipt+1)                  !jmecmplx
                end if
             end do

          end do loop_ngwf_a

       end do loop_spins

       ! -------------------------------------------
       ! jd: Apply reciprocal-space preconditioning
       ! -------------------------------------------

       if(pub_precond_recip) then

          call utils_assert(present(precond_func_recip), myself//': precond_&
               &func_recip must be passed when pub_precond_recip is .true.')
          ! jd: Allocate workspace for ket in FFTbox around Aa
          n_ngwfs_a = alpha_ngwf_index%basis%num_on_atom(global_a)
          allocate(kets_in_fftbox(n_ngwfs_a, pub_num_spins), stat=ierr)
          call utils_alloc_check(myself,'kets_in_fftbox',ierr)
          do is = 1, pub_num_spins
             do ngwf_a = 1, alpha_ngwf_index%basis%num_on_atom(global_a)
                call data_fftbox_alloc(kets_in_fftbox(ngwf_a, is), &
                     fftbox%total_ld1, fftbox%total_ld2, fftbox%total_pt3, &
                     iscmplx=contra_grad%iscmplx)      !jmecmplx
             end do
          end do

          ! jme: allocate work arrays
          allocate(work_fftbox(fftbox%total_ld1, fftbox%total_ld2, &
               fftbox%total_pt3), stat=ierr)
          call utils_alloc_check(myself, 'work_fftbox', ierr)
          allocate(zwork_fftbox(fftbox%total_ld1, fftbox%total_ld2, &
               fftbox%total_pt3), stat=ierr)
          call utils_alloc_check(myself, 'zwork_fftbox', ierr)

          ! --------------------------------------------------------------------
          ! jd: for all a on A                                               aaa
          ! --------------------------------------------------------------------
          loop_ngwf_a_2:                                                          &
          do ngwf_a = 1, alpha_ngwf_index%basis%num_on_atom(global_a)

             global_ngwf_aa = &
                  alpha_ngwf_index%basis%first_on_atom(global_a) + ngwf_a - 1
             local_ngwf_aa = global_ngwf_aa - &
                  alpha_ngwf_index%basis%first_on_proc(pub_my_proc_id) + 1

             ! jd: Determine the position of the FFTbox around Aa and the position
             !     of Aa in the FFTbox
             call function_ops_batch_col_start(fa_box_start, fa_start_in_box, &
                  1, local_ngwf_aa, local_ngwf_aa, fftbox, cell, &
                  alpha_ngwf_index%basis)

             ! -----------------------------------------------------------------
             ! jd: for all spins                                             sss
             ! -----------------------------------------------------------------
             loop_spins_2:                                                     &
             do is = 1, pub_num_spins

#ifdef HFX_EXTRA_SHAVE
  ! jd: We could be additionally shaving at this point, before we
  !     FFT to recip space. This helps remove the tiny dependence
  !     of energy on position wrt PPD (so, there is a separate
  !     small eggbox effect even when we move the system by an
  !     integer number of psincs), but otherwise yields worse NGWF
  !     gradients and worse final ground state energies.
  !     Rough idea about the magnitude of the eggbox effect (for acetonitrile):
  ! PBE:
  ! eggbox amplitude: 0.1 kcal/mol
  ! eggbox amplitude when moving by an integer multiple of psincs: 1E-8 kcal/mol
  ! HF:
  ! eggbox amplitude: 0.4 kcal/mol
  ! eggbox amplitude when moving by an integer multiple of psincs: 5E-4 kcal/mol
  ! HF with this additional shaving:
  ! eggbox amplitude: 0.4 kcal/mol
  ! eggbox amplitude when moving by an integer multiple of psincs: 2E-6 kcal/mol
  !
  !  Total energy ~= -15000 kcal/mol

                if (cell%n_pts > 1) then
                   call basis_clean_function(kets_in_ppds_on_a(2,ngwf_a,is), &
                        alpha_ngwf_index%basis%spheres(local_ngwf_aa), &
                        cell, fftbox, 1)
                end if
#endif

                ! jd: Transfer ket in PPDs on A to an FFTbox around Aa. Only the
                !     tcK component is needed (which_kernel==2).
                call basis_copy_function_to_box(kets_in_fftbox(ngwf_a,is), & ! target
                     fa_start_in_box(1,1), &
                     fa_start_in_box(2,1), &
                     fa_start_in_box(3,1), &
                     alpha_ngwf_index%basis%tight_boxes(local_ngwf_aa), &
                     kets_in_ppds_on_a(2,ngwf_a,is), &
                     alpha_ngwf_index%basis%spheres(local_ngwf_aa), &
                     cell, fftbox, 1)
             end do loop_spins_2

             work_fftbox = 0.0_DP
             zwork_fftbox = (0.0_DP,0.0_DP)

             if (pub_num_spins == 2) then
                ! qoh: Forward FFT the covariant gradient to reciprocal space
                call fourier_apply_box_pair('C','Forward', &
                     kets_in_fftbox(ngwf_a,1)%d, &               !jmecmplx
                     kets_in_fftbox(ngwf_a,2)%d, zwork_fftbox)   !jmecmplx
                ! qoh: Apply kinetic energy preconditioning to covariant gradient
                zwork_fftbox = zwork_fftbox * precond_func_recip
                ! qoh: Backward FFT the covariant gradient to real space
                call fourier_apply_box_pair('C','Backward', &
                     kets_in_fftbox(ngwf_a,1)%d,  &              !jmecmplx
                     kets_in_fftbox(ngwf_a,2)%d, zwork_fftbox)   !jmecmplx
             else
                ! qoh: Forward FFT the covariant gradient to reciprocal space
                call fourier_apply_box_pair('C','Forward', &
                     kets_in_fftbox(ngwf_a,1)%d, &               !jmecmplx
                     work_fftbox, zwork_fftbox)
                ! qoh: Apply kinetic energy preconditioning to covariant gradient
                zwork_fftbox = zwork_fftbox * precond_func_recip
                ! qoh: Backward FFT the covariant gradient to real space
                call fourier_apply_box_pair('C','Backward', &
                     kets_in_fftbox(ngwf_a,1)%d,  &              !jmecmplx
                     work_fftbox, zwork_fftbox)
             end if

          end do loop_ngwf_a_2

          ! --------------------------------------------------------------------
          ! jd: for all spins                                                sss
          ! --------------------------------------------------------------------
          loop_spins_3:                                                        &
          do is = 1, pub_num_spins

             ! -----------------------------------------------------------------
             ! jd: for all a on A                                            aaa
             ! -----------------------------------------------------------------
             loop_ngwf_a_3:                                                    &
             do ngwf_a = 1, alpha_ngwf_index%basis%num_on_atom(global_a)

                global_ngwf_aa = &
                     alpha_ngwf_index%basis%first_on_atom(global_a) + ngwf_a-1
                local_ngwf_aa = global_ngwf_aa - &
                     alpha_ngwf_index%basis%first_on_proc(pub_my_proc_id)+1

                call data_set_to_zero(aket_dkn_in_ppds(2))
                ! cks: extract ppds belonging to bra function from ket fftbox
                call basis_extract_function_from_box(aket_dkn_in_ppds(2), &
                     kets_in_fftbox(ngwf_a,is), &                 !jmecmplx
                     alpha_ngwf_index%basis%spheres(local_ngwf_aa), &
                     alpha_ngwf_index%basis%tight_boxes(local_ngwf_aa), &
                     fa_start_in_box(1,1), &
                     fa_start_in_box(2,1), &
                     fa_start_in_box(3,1), 1, cell, fftbox)

                ! cks: shaving: zero points outside NGWF sphere in PPD rep.
                if (cell%n_pts > 1) then
                   call basis_clean_function(aket_dkn_in_ppds(2), &
                        alpha_ngwf_index%basis%spheres(local_ngwf_aa), &
                        cell, fftbox, 1)
                end if

                dot_npts = alpha_ngwf_index%basis%spheres(local_ngwf_aa)%&
                     n_ppds_sphere * cell%n_pts

                bra_start = alpha_ngwf_index%basis%spheres(local_ngwf_aa)%offset
                do ipt = 0, dot_npts - 1
                   cov_grad%d(bra_start+ipt) = cov_grad%d(bra_start+ipt)+ & ! sic, plus  !jmecmplx
                        factor * aket_dkn_in_ppds(2)%d(ipt+1)                            !jmecmplx
                end do

             end do loop_ngwf_a_3

          end do loop_spins_3

          ! jme: deallocate work arrays
          deallocate(zwork_fftbox, stat=ierr)
          call utils_dealloc_check(myself, 'zwork_fftbox', ierr)
          deallocate(work_fftbox, stat=ierr)
          call utils_dealloc_check(myself, 'work_fftbox', ierr)

          do is = pub_num_spins, 1, -1
             do ngwf_a = alpha_ngwf_index%basis%num_on_atom(global_a), 1, -1
                call data_fftbox_dealloc(kets_in_fftbox(ngwf_a, is))
             end do
          end do

          deallocate(kets_in_fftbox, stat=ierr)
          call utils_dealloc_check(myself,'kets_in_fftbox',ierr)

          !jmecmplx
          do is = 1, pub_num_spins
             do ngwf_a = 1, n_ngwfs_a
                do which_kernel = 1, 2
                   call data_functions_dealloc(kets_in_ppds_on_a(&
                        which_kernel,ngwf_a,is))
                end do
             end do
          end do

          deallocate(kets_in_ppds_on_a, stat=ierr)
          call utils_dealloc_check(myself,'kets_in_ppds_on_a',ierr)

       end if ! pub_precond_recip

    end do loop_A

    ! jme: deallocate buffers
    do which_kernel = 2, 1, -1
       call data_functions_dealloc(aket_dkn_in_ppds(which_kernel))
    end do

#if 0
    print *, "-NGWF_CONTRA_TOT_POST:-============================"
    print *, contra_grad%d

    print *, "-NGWF_COV_TOT_POST:-============================"
    print *, cov_grad%d
#endif

    call comms_barrier
    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hfx_gradient

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_gradient_term2(hfxstate, grad_akets_dkn_my_ht, &
       packed_dkn_cd_blocks_ht, packed_dkn_ab_blocks_hts, ngwfs_hts, &
       swex_quality, src_swex_h, swri, cell, par, elements, &
       ngwf_indices, grad_prefactor)
    !==========================================================================!
    ! Calculates term 2 (second line of (34) in [1]) of the NGWF gradient.     !
    ! This now follows the new parallelisation scheme, so we are only dealing  !
    ! with our B-C pairs. To get better alignment with this distribution, and  !
    ! to re-use B and C SWOPs that we already have, we rename indices in (34)  !
    ! so that we don't work with A and D, but with B and C. Indices are        !
    ! renamed as follows:                                                      !
    !     A <-> B                                                              !
    !     C <-> D                                                              !
    ! In the end we obtain DKN-reduced gradient kets for all my B atoms. These !
    ! are then redistributed by hfx_comms_gradient_term2() to the sparse-local !
    ! scheme, and that's where we back-rename the last remaining index (B->A). !
    !                                                                          !
    ! If any loop ordering appears unnatural or inconvenient, this is because  !
    ! we go to great lengths to minimise re-calculation of quantities without  !
    ! introducing a significant memory load. For instance, the outermost loop  !
    ! goes over an inconvenient index C because we do not want to recalculate  !
    ! P_AC needlessly. If the loop goes over C, we only store P_AC for current !
    ! C, not for all C's. The loop cannot go over A because then we'd be       !
    ! recalculating and re-expanding Q_BC. We use similar reasoning for all    !
    ! indices, gradually eliminating them by doing summations.                 !
    !                                                                          !
    ! The general idea is this:                                                !
    ! for all C {                                                              !
    !   Calculate P_{Aac;C} for all relevant A. This eliminates index Dd.      !
    !   for all c {                                                            !
    !     for relevant B {                                                     !
    !       for relevant A {                                                   !
    !         accumulate in Q_{Bb,wk,s;Cc} the sum_a K^BbAa P_{Aac;C}.         !
    !       }                                                                  !
    !     }                                                                    !
    !     Index Aa is eliminated.                                              !
    !     for relevant B {                                                     !
    !       for all spins {                                                    !
    !         for which_kernel (K/tcK) {                                       !
    !           for b {
    !             Expand Q_{Bb,wk,s;Cc}.
    !             Calculate potential of expansion in PPDs common to B and C.
    !             Accumulate \phi_Cc * pot. in PPDs on B (kets_in_ppds_on_bs)
    !           } b
    !         } wk
    !       } s
    !     } B
    !   } c
    ! } C
    ! Index Cc is eliminated.
    ! Convert from array to HT, rename to akets_dkn_my.                        !
    !                                                                          !
    ! Results are stored in grad_akets_dkn_my_ht. They are for all my B atoms  !
    ! (that will become local A atoms after comms_gradient(), hence the name   !
    ! 'akets').                                                                !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (in): HFx state. Needed for neighbour lists and my atom lists.!
    !   grad_akets_dkn_my_ht (inout): Results go here.                         !
    !   packed_dkn_cd_blocks_ht (in): HT that stores required K^CcDd blocks.   !
    !                                 Needed for calculating P terms.          !
    !   packed_dkn_ab_blocks_hts (in): HT that stores required K^AaBB and tcK  !
    !                                  blocks. Needed for calculating Q terms. !
    !   ngwfs_hts (in): NGWF caches. Only needed for B's PPD indices in the    !
    !                   calc of Q.                                             !
    !   swex_quality (in): Quality of the expansion.                           !
    !   src_swex_h (in): Needed in the expansion.                              !
    !   swri (inout): Needed in the expansion. Caches SWOPs behind the scenes. !
    !   cell (in): Needed in the expansion and for n_pts.                      !
    !   par (in): Needed to establish who sparse-owns which atoms.             !
    !   elements (in): Needed for expansion centres.                           !
    !   ngwf_indices (in): Array of HFX_NGWF_INDEX instances containing        !
    !                      pointers to the NGWF basis and the NGWF_REP used    !
    !                      by atoms A, B, C and D.                             !
    !   grad_prefactor (in): Defaults to 1.0. In conduction 0.5 can be passed  !
    !   [optional]           to take into account the fact that only one set of!
    !                        NGWFs is optimised, and the other one is fixed    !
    !                        (I think).                                        !
    !                        This also takes care of spin degeneracy conven-   !
    !                        tions in conduction (see the call to              !
    !                        hf_exchange_calculate from conduction).           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May-July 2019.                              !
    !==========================================================================!

    use comms, only: comms_barrier, pub_my_proc_id
    use constants, only: PI, SW_O, SW_V, METRIC_ELECTROSTATIC, METRIC_OVERLAP, &
         SWEX_BB_CC_SYMMETRIC, VERBOSE, CRLF, stdout, LONG
    use hash_table, only: HT_HASH_TABLE, hash_table_add
    use ion, only: ELEMENT
    use neighbour_list, only: neighbour_list_are_neighbours, &
         neighbour_list_for_atom_list
    use parallel_strategy, only: PARAL_INFO
    use ppd_ops, only: ppd_set_intersection_unpacked, ppd_accumulate, &
         ppd_extract_contents_of_subset
    use remote, only: remote_ppd_list_of_atom, remote_n_ppds_of_atom, &
         remote_obtain_unpacked_ngwfs
    use rundat, only: pub_num_spins, pub_hfx_metric, pub_swx_dbl_grid, &
         pub_threads_max, pub_hfx_output_detail, pub_debug
    use simulation_cell, only: CELL_INFO
    use sw_expansion, only: swx_expansion_2_fast, swx_rms_fitting_error, &
         swx_calculate_rms_fitting_error
    use sw_resolution_of_identity, only: ATOM_CENTRE, SW_RI, SW_QUALITY, &
         swri_init_centre, swri_expansion_centres, swri_build_matrix_2
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_alloc_check, &
         utils_dealloc_check, utils_assert, utils_sanity_check, utils_abort, &
         utils_unit

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(in), target        :: hfxstate
    type(HT_HASH_TABLE), intent(inout), target :: grad_akets_dkn_my_ht
    type(HT_HASH_TABLE), intent(in), target    :: packed_dkn_cd_blocks_ht
    type(HT_HASH_TABLE), intent(in), target    :: packed_dkn_ab_blocks_hts(3)
    type(HT_HASH_TABLE), intent(in), target    :: ngwfs_hts(:)
    type(SW_QUALITY), intent(in)               :: swex_quality
    integer, intent(in)                        :: src_swex_h
    type(SW_RI), intent(in)                    :: swri
    type(CELL_INFO), intent(in)                :: cell
    type(PARAL_INFO), intent(in)               :: par
    type(ELEMENT), intent(in)                  :: elements(:)
    type(HFX_NGWF_INDEX), intent(in)           :: ngwf_indices(A_NGWF:D_NGWF)
    real(kind=DP), optional, intent(in)        :: grad_prefactor

    ! jd: Local allocatables
    real(kind=DP), allocatable :: aa_ngwfs_in_ppds(:) ! idx: PPD-pt-relev_Aa
    real(kind=DP), allocatable :: dd_ngwfs_in_ppds(:) ! idx: PPD-pt-relev_Dd
    real(kind=DP), allocatable :: c_ngwfs_in_ppds(:)  ! idx: PPD-pt-ngwf_c
    real(kind=DP), allocatable :: qs_in_ppds(:) ! fused index: relevB-ngwf_b-wk-is-PPD-pt
    real(kind=DP), allocatable :: p_ac(:) ! fused index: PPDs_on_A-a-c-is
    real(kind=DP), allocatable :: kets_in_ppds_on_bs(:,:,:,:) ! pt-PPD, wk, is, my_ngwf_b
    integer, allocatable       :: offsets_to_aa_ngwfs(:) ! idx: relev_A
    integer, allocatable       :: offsets_to_dd_ngwfs(:) ! idx: relev_D
    integer, allocatable       :: dummy(:)
    integer, allocatable       :: ppd_indices_in_qs(:,:)
    integer, allocatable       :: n_ppds_in_qs(:)
    integer, allocatable       :: relev_a_atoms(:) ! relev_a -> global_a
    integer, allocatable       :: relev_d_atoms(:) ! relev_d -> global_d

    ! jd: Local stack-based arrays
    real(kind=DP) :: q_expansion_in_bc_ppds(min(&
         ngwf_indices(B_NGWF)%basis%max_n_ppds_sphere, &
         ngwf_indices(C_NGWF)%basis%max_n_ppds_sphere) * cell%n_pts)
    real(kind=DP) :: ngwf_cc_in_bc_ppds(min(&
         ngwf_indices(B_NGWF)%basis%max_n_ppds_sphere, &
         ngwf_indices(C_NGWF)%basis%max_n_ppds_sphere) * cell%n_pts)
    real(kind=DP) :: data_to_add(cell%n_pts+4)
    integer :: relev_b_atoms(hfxstate%n_my_b_atoms) ! relev_b -> global_b
    integer :: relev_b_atoms_my(hfxstate%n_my_b_atoms) ! relev_b -> my_b
    integer(kind=LONG) :: qs_offsets(hfxstate%n_my_b_atoms) ! indexed by relev_b
    integer(kind=LONG) :: p_offsets(hfxstate%n_my_a_atoms) ! indexed by relev_a
    integer :: ppd_indices_on_as(ngwf_indices(A_NGWF)%basis%max_n_ppds_sphere, &
         hfxstate%n_my_a_atoms) ! indexed by relev_a
    integer :: ppd_indices_bb_cc_intersection(min(&
         ngwf_indices(B_NGWF)%basis%max_n_ppds_sphere, &
         ngwf_indices(C_NGWF)%basis%max_n_ppds_sphere))
    integer :: ppd_indices_on_b(ngwf_indices(B_NGWF)%basis%max_n_ppds_sphere)
    integer :: ppd_indices_on_c(ngwf_indices(C_NGWF)%basis%max_n_ppds_sphere)

    integer :: n_ppds_bb_cc_intersection
    integer :: is
    integer :: which_kernel
    integer :: ndata_in_q_one_ngwf

    ! jd: Atom and NGWF indices
    integer :: c_idx
    integer :: my_b
    integer :: relev_a, relev_b, relev_d
    integer :: global_a, global_b, global_c, global_d
    integer :: first_ngwf_idx_of_a, first_ngwf_idx_of_b, first_ngwf_idx_of_c
    integer :: global_bb_ngwf_idx, global_cc_ngwf_idx
    integer :: n_ngwfs_a, n_ngwfs_b, n_ngwfs_c
    integer :: ngwf_a, ngwf_b, ngwf_c ! for current atom
    integer :: my_ngwf_b              ! for all my B atoms
    integer :: n_relev_a_atoms, n_relev_b_atoms, n_relev_d_atoms
    integer :: n_relev_aa_ngwfs

    ! jd: PPD counting
    integer(kind=LONG) :: n_points_in_p
    integer :: n_points_in_bb_cc_intersection
    integer :: n_ppds_on_as(hfxstate%n_my_a_atoms) ! indexed by relev_a
    integer :: n_ppds_on_ds(hfxstate%n_my_d_atoms) ! indexed by relev_d
    integer :: n_ppds_on_b
    integer :: n_ppds_on_c
    integer :: n_ppds_in_p_from_this_a
    integer :: n_ppds_in_p
    integer :: n_ppds_in_qset
    integer :: n_ppds_in_q_from_this_b
    integer(kind=LONG) :: n_points_in_qset
    integer :: spin_wk_ngwf_b_offset
    integer :: dump_offset
    integer :: offset_to_ngwf_c
    integer :: offset_to_ppd
    integer :: i_ppd
    integer :: cur_ppd

    ! jd: Expansion stuff
    type(ATOM_CENTRE) :: expansion_centres(2), centres(2)
    integer :: coeffs_kind
    integer :: expansion_atoms(2), atoms(2)
    integer :: num_sws_in_expansion
    real(kind=DP), allocatable :: vmatrix(:,:)
    real(kind=DP), allocatable :: omatrix(:,:)
    real(kind=DP), allocatable :: wmatrix(:,:)
    real(kind=DP)              :: q_coeffs(2*swri%quality%max_sws_per_centre)

    ! jd: Varia
    real(kind=DP) :: loc_grad_prefactor
    real(kind=DP) :: factor
    real(kind=DP) :: rms
    integer :: n_bc_pairs_processed
    integer :: ierr
    character(len=*), parameter :: myself = 'hfx_gradient_term2'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    if(pub_hfx_metric == METRIC_OVERLAP) coeffs_kind = SW_O
    if(pub_hfx_metric == METRIC_ELECTROSTATIC) coeffs_kind = SW_V

    ! *************************
    ! jd: Index renaming zone:
    !     A <-> B
    !     C <-> D
    ! *************************

    allocate(kets_in_ppds_on_bs(ngwf_indices(B_NGWF)%basis%max_n_ppds_sphere*&
         cell%n_pts, 2, pub_num_spins, hfxstate%n_my_bb_ngwfs), stat=ierr)
    call utils_alloc_check(myself, 'kets_in_ppds_on_bs', ierr)
    kets_in_ppds_on_bs(:,:,:,:) = 0.0_DP

    n_bc_pairs_processed = 0

    ! --------------------------------------------------------------------------
    ! jd: Loop over my C atoms                                               CCC
    ! --------------------------------------------------------------------------
    loop_C:                                                                    &
    do c_idx = 1, hfxstate%n_my_c_atoms
       global_c = hfxstate%my_c_atoms(c_idx)
       first_ngwf_idx_of_c = ngwf_indices(C_NGWF)%basis%first_on_atom(global_c)
       n_ngwfs_c = ngwf_indices(C_NGWF)%basis%num_on_atom(global_c)

       if(pub_hfx_output_detail >= VERBOSE) then
          write(stdout,'(a,i7,a,i0,a)') '  * C: ', global_c, ' [', &
               pub_my_proc_id, ']'
       end if

       ! *******************************************************************
       ! jd: Start a new set of Q (current C, all my B).
       ! *******************************************************************

       call timer_clock(myself//'_part0',1)

       ! *******************************************************************
       ! jd: Just book-keeping
       ! *******************************************************************

       ! jd: Obtain PPD indices on C (not necessarily local)
       call remote_ppd_list_of_atom(ppd_indices_on_c, n_ppds_on_c, &
            global_c, ngwf_indices(C_NGWF)%basis, &
            ngwf_indices(C_NGWF)%rep%ngwf_cache_handle)

       ! jd: Figure out PPDs that go into this set of Q
       call hfx_determine_ppd_indices_in_qs(ppd_indices_in_qs, n_ppds_in_qs,&
            global_c, hfxstate, ngwf_indices, cell%n_ppds)

       ! jd: Figure out relevant B atoms -- those, who have >0 PPDs and feature
       !     in our BC pair list, so B atoms that are present in this Q set.
       n_relev_b_atoms = 0
       relev_b_atoms(:) = garbage_int
       relev_b_atoms_my(:) = garbage_int
       do my_b = 1, hfxstate%n_my_b_atoms
          if(n_ppds_in_qs(my_b) > 0) then
             n_relev_b_atoms = n_relev_b_atoms + 1
             global_b = hfxstate%my_b_atoms(my_b)
             relev_b_atoms(n_relev_b_atoms) = global_b
             relev_b_atoms_my(n_relev_b_atoms) = my_b
          end if
       end do

       ! jd: Figure out the dimensions of the Q set.
       n_ppds_in_qset = 0
       n_points_in_qset = 0_LONG
       qs_offsets(:) = garbage_int
       do relev_b = 1, n_relev_b_atoms
          global_b = relev_b_atoms(relev_b)
          my_b = relev_b_atoms_my(relev_b)
          n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)
          n_ppds_in_q_from_this_b = &
               n_ppds_in_qs(my_b) * n_ngwfs_b * 2 * pub_num_spins
          n_ppds_in_qset = n_ppds_in_qset + n_ppds_in_q_from_this_b
          qs_offsets(relev_b) = n_points_in_qset + 1
          n_points_in_qset = n_points_in_qset + &
               int(n_ppds_in_q_from_this_b * cell%n_pts, kind=LONG)
       end do

       ! jd: Figure out relevant A's -- s-X-neighbours of C that are also
       !     my-A atoms. The 2nd condition ([*]) is needed to exclude A's that
       !     are C's s-X-neighbours but only via a B that is not my B. These are
       !     not our responsibility.
       !     Determine the dimension in P, offsets to different A's.
       n_relev_aa_ngwfs = 0
       n_ppds_in_p = 0
       n_points_in_p = 0_LONG
       p_offsets(:) = garbage_int
       n_ppds_on_as(:) = garbage_int
       ppd_indices_on_as(:,:) = garbage_int

       call neighbour_list_for_atom_list(relev_a_atoms, n_relev_a_atoms, &
            relev_b_atoms, n_relev_b_atoms, hfxstate%x_atoms_nl, par%nat)

       do relev_a = 1, n_relev_a_atoms
          global_a = relev_a_atoms(relev_a)
          n_ngwfs_a = ngwf_indices(A_NGWF)%basis%num_on_atom(global_a)
          n_relev_aa_ngwfs = n_relev_aa_ngwfs + n_ngwfs_a

          ! jd: Obtain PPD indices on A (not necessarily local)
          call remote_ppd_list_of_atom(ppd_indices_on_as(:,relev_a), &
               n_ppds_on_as(relev_a), global_a, ngwf_indices(A_NGWF)%basis, &
               ngwf_indices(A_NGWF)%rep%ngwf_cache_handle)

          n_ppds_in_p_from_this_a = n_ppds_on_as(relev_a) * &
               n_ngwfs_a * n_ngwfs_c * pub_num_spins
          n_ppds_in_p = n_ppds_in_p + n_ppds_in_p_from_this_a
          p_offsets(relev_a) = n_points_in_p + 1
          n_points_in_p = n_points_in_p + &
               int(n_ppds_in_p_from_this_a * cell%n_pts, kind=LONG)
          if(n_points_in_p < 0) then
             call utils_abort(myself//': Overflow in auxiliary quantity P.'//&
                  CRLF//'Increase the number of MPI ranks (processors) in &
                  &your calculation (preferably), or reduce your PPD size. &
                  &As a last resort, reduce the exchange cutoff (hfx_cutoff).')
          end if
       end do

       ! jd: Unpack all relevant Aa NGWFs
       call remote_obtain_unpacked_ngwfs(aa_ngwfs_in_ppds, offsets_to_aa_ngwfs,&
            n_relev_a_atoms, relev_a_atoms, n_ppds_on_as, &
            ngwf_indices(A_NGWF)%basis, &
            ngwf_indices(A_NGWF)%rep%ngwf_cache_handle, cell%n_pts)

       ! jd: Figure out relevant D atoms. These are S-neighbours of all relevant
       !     A atoms.
       call neighbour_list_for_atom_list(relev_d_atoms, n_relev_d_atoms, &
            relev_a_atoms, n_relev_a_atoms, hfxstate%s_atoms_nl, par%nat)

       ! jd: Obtain number of PPDs on relevant D's (not necessarily local)
       do relev_d = 1, n_relev_d_atoms
          global_d = relev_d_atoms(relev_d)
          call remote_n_ppds_of_atom(n_ppds_on_ds(relev_d), global_d, &
               ngwf_indices(D_NGWF)%basis, &
               ngwf_indices(D_NGWF)%rep%ngwf_cache_handle)
       end do

       ! jd: Unpack all relevant Dd NGWFs
       call remote_obtain_unpacked_ngwfs(dd_ngwfs_in_ppds, offsets_to_dd_ngwfs,&
            n_relev_d_atoms, relev_d_atoms, n_ppds_on_ds, &
            ngwf_indices(D_NGWF)%basis, &
            ngwf_indices(D_NGWF)%rep%ngwf_cache_handle, cell%n_pts)

       ! jd: Unpack all c NGWFs
       call remote_obtain_unpacked_ngwfs(c_ngwfs_in_ppds, dummy, &
            1, [global_c], [n_ppds_on_c], ngwf_indices(C_NGWF)%basis, &
            ngwf_indices(C_NGWF)%rep%ngwf_cache_handle, cell%n_pts)

       call timer_clock(myself//'_part0',2)

       ! *******************************************************************
       ! *** PART 1: Calculating P                                       ***
       ! *******************************************************************

       call timer_clock(myself//'_part1',1)

       ! jd: Dimension of P: Relevant A atoms * their ngwfs * n_ngwfs_c * spins
       allocate(p_ac(max(n_points_in_p,1_LONG)),stat=ierr)
       call utils_alloc_check(myself,'p_ac',ierr)
       p_ac(:) = 0.0_DP

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(DYNAMIC) &
!$OMP PRIVATE(relev_a, global_a, dump_offset, ngwf_a, ngwf_c, n_ngwfs_a)       &
!$OMP SHARED(hfxstate, ngwf_indices, global_c, cell, packed_dkn_cd_blocks_ht,  &
!$OMP      relev_a_atoms, n_relev_a_atoms, p_ac, p_offsets, n_ppds_on_as,      &
!$OMP      ppd_indices_on_as, n_ngwfs_c, pub_num_spins, aa_ngwfs_in_ppds,      &
!$OMP      dd_ngwfs_in_ppds, offsets_to_aa_ngwfs, offsets_to_dd_ngwfs,         &
!$OMP      relev_d_atoms, n_relev_d_atoms)
       ! -----------------------------------------------------------------------
       ! jd: Loop over A's that are relevant                                 AAA
       ! -----------------------------------------------------------------------
       loop_A1:                                                                &
       do relev_a = 1, n_relev_a_atoms
          global_a = relev_a_atoms(relev_a)

          ! *************************************************************
          ! jd: Calc P_Aac;C for all {a,c}. This eliminates index Dd.
          !     Each thread writes to a different slice (different A).
          ! *************************************************************
          call hfx_calc_p_ac(p_ac, hfxstate%s_atoms_nl, &
               packed_dkn_cd_blocks_ht, cell%n_pts, global_a, global_c, &
               n_ppds_on_as(relev_a), ppd_indices_on_as(:,relev_a), &
               p_offsets(relev_a), aa_ngwfs_in_ppds, &
               offsets_to_aa_ngwfs(relev_a), dd_ngwfs_in_ppds, &
               offsets_to_dd_ngwfs, relev_d_atoms, n_relev_d_atoms, &
               ngwf_indices)

#if 0
          ! jd: Debug dump of P's for comparison against old implementation
          n_ngwfs_a = ngwf_indices(A_NGWF)%basis%num_on_atom(global_a)
          do ngwf_a = 1, n_ngwfs_a
             do ngwf_c = 1, n_ngwfs_c
                dump_offset = p_offsets(relev_a) + &
                     ((ngwf_a-1) * n_ngwfs_c * pub_num_spins + &
                     (ngwf_c-1) * pub_num_spins + &
                     (1-1)) * n_ppds_on_as(relev_a) * cell%n_pts

                write(filename,'(a,i0,a,i0,a,i0,a,i0)') &
                     'p_', global_a, '_', ngwf_a, '_', global_c, '_', ngwf_c
                u = utils_unit()
                open(u, file=filename)
                write(u,*) p_ac(dump_offset:dump_offset+&
                     n_ppds_on_as(relev_a)*cell%n_pts-1)
                close(u)
             end do
          end do
#endif

       end do loop_A1
!$OMP END PARALLEL DO

       deallocate(offsets_to_dd_ngwfs, stat=ierr)
       call utils_dealloc_check(myself, 'offsets_to_ngwfs', ierr, &
            allocated_in = 'remote_obtain_unpacked_ngwfs')
       deallocate(dd_ngwfs_in_ppds, stat=ierr)
       call utils_dealloc_check(myself, 'ngwfs_in_ppds', ierr, &
            allocated_in = 'remote_obtain_unpacked_ngwfs')
       deallocate(offsets_to_aa_ngwfs, stat=ierr)
       call utils_dealloc_check(myself, 'offsets_to_ngwfs', ierr, &
            allocated_in = 'remote_obtain_unpacked_ngwfs')
       deallocate(aa_ngwfs_in_ppds, stat=ierr)
       call utils_dealloc_check(myself, 'ngwfs_in_ppds', ierr, &
            allocated_in = 'remote_obtain_unpacked_ngwfs')
       deallocate(relev_d_atoms, stat=ierr)
       call utils_dealloc_check(myself, 'relev_d_atoms', ierr, &
            allocated_in = 'neighbour_list_for_atom_list')

       ! *******************************************************************
       ! *** END OF PART 1                                               ***
       ! *******************************************************************

       ! ####################################
       ! ### Index Dd has been eliminated ###
       ! ####################################

       call timer_clock(myself//'_part1',2)

       ! -----------------------------------------------------------------------
       ! jd: All NGWFs c of C                                                ccc
       ! -----------------------------------------------------------------------
       offset_to_ngwf_c = 1
       loop_ngwf_c:                                                            &
       do ngwf_c = 1, n_ngwfs_c
          global_cc_ngwf_idx = first_ngwf_idx_of_c + ngwf_c - 1

          ! *******************************************************************
          ! *** PART 2: Calculating Q                                       ***
          ! *******************************************************************

          call timer_clock(myself//'_part2',1)

          ! jd: This reaches ~300MB in a 171-atom system. This is with 4 NGWFs
          !     and 33 relevant B's. This will not become much bigger for larger
          !     systems, as ~30 B's on a rank is already plenty. For a 701-atom
          !     system on 16 MPI ranks I got 800MB, but systems that big should
          !     be using more ranks already. It will, however become bigger in:
          !     - spin-polarised systems (by a factor of 2),
          !     - conduction calculations (by a factor of maybe 15/4 from NGWFs
          !       and maybe (12/8)^3 from NGWF radii). This yields ~3800MB,
          !       which is largish, but bearable.
          !     However, it means there is no hope of
          !     - spawning a copy of Q per OMP thread,
          !     - promoting 'c' from 'current' index to storing Q for all 'c'.
          !     If the memory load from Q ever becomes problematic, workable
          !     solutions include:
          !     - demoting 'b' from storing all 'b' to 'current' index.
          !     - ditto for spin.
          !     - batching over B (permitting any desired saving factor in mem
          !       load for Q at the price of repeating the calculation of P_AC
          !       that many times).
          call utils_assert(n_points_in_qset >= 0, &
               'Overflow in HFx auxiliary quantity Q')
          allocate(qs_in_ppds(max(n_points_in_qset,1_LONG)), stat=ierr)
          call utils_alloc_check(myself,'qs_in_ppds',ierr)
          qs_in_ppds(:) = 0.0_DP

          ! --------------------------------------------------------------------
          ! jd: Loop over relevant B's                                       BBB
          ! --------------------------------------------------------------------
          ! jd: This loop OMP-parallelises very poorly. This is because there
          !     are insufficient B's -- often fewer than the number of OMP
          !     threads. For instance in a 171-atom system on 16 MPI ranks,
          !     there were 1-31 B's, with an average of 11.5. Moreover, the cost
          !     for each B varies by quite a lot depending on the # of A's and
          !     # of NGWFs, so balancing is poor. Fortunately, this stage is
          !     comparably short, part 3 is much more intensive.
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(DYNAMIC) &
!$OMP PRIVATE(relev_b, global_b, first_ngwf_idx_of_b, n_ngwfs_b, relev_a,      &
!$OMP      global_a, first_ngwf_idx_of_a, n_ngwfs_a)                           &
!$OMP SHARED(n_relev_a_atoms, n_relev_b_atoms, ngwf_indices, hfxstate, cell,   &
!$OMP      qs_in_ppds, qs_offsets, p_offsets, p_ac, global_c, ngwf_c,          &
!$OMP      ppd_indices_in_qs, n_ppds_in_qs, ppd_indices_on_as, n_ppds_on_as,   &
!$OMP      packed_dkn_ab_blocks_hts, relev_a_atoms, relev_b_atoms)
          loop_B:                                                              &
          do relev_b = 1, n_relev_b_atoms
             global_b = relev_b_atoms(relev_b)
             first_ngwf_idx_of_b = &
                  ngwf_indices(B_NGWF)%basis%first_on_atom(global_b)
             n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)

             ! -----------------------------------------------------------------
             ! jd: Loop over A's that are relevant                           AAA
             ! -----------------------------------------------------------------
             loop_A2:                                                          &
             do relev_a = 1, n_relev_a_atoms
                global_a = relev_a_atoms(relev_a)
                first_ngwf_idx_of_a = ngwf_indices(A_NGWF)%basis%first_on_atom(global_a)
                n_ngwfs_a = ngwf_indices(A_NGWF)%basis%num_on_atom(global_a)

                ! jd: Exclude A's that are relevant, but not through this
                !     particular B. This ensures X sparsity is maintained,
                !     that is when calculating Q_BbCc = \sum_A ..., which then
                !     determines the NGWF gradient on Bb, we only go over A's
                !     that actually interact (are not outside of X sparsity).
                !     An alternative would be to go over X-neighbours of B in
                !     the A2 loop above, but then we'd have to eliminate A's
                !     that are not S-neighbours of current C.
                if(.not. neighbour_list_are_neighbours(global_b, global_a, &
                     hfxstate%x_atoms_nl)) cycle

                ! **********************************************************
                ! jd: Accumulate into Q_{Bb,wk,is};Cc
                ! **********************************************************
                ! Each thread writes to a different B-slice.
                call hfx_accum_q_beta(qs_in_ppds(:), &
                     qs_offsets(relev_b), p_offsets(relev_a), p_ac, &
                     global_b, global_a, global_c, ngwf_c, ppd_indices_in_qs, &
                     n_ppds_in_qs, ppd_indices_on_as(:,relev_a), &
                     n_ppds_on_as(relev_a), cell, ngwf_indices, hfxstate, &
                     packed_dkn_ab_blocks_hts)

             end do loop_A2

          end do loop_B
!$OMP END PARALLEL DO

          call timer_clock(myself//'_part2',2)

          ! ####################################
          ! ### Index Aa has been eliminated ###
          ! ####################################

          ! *******************************************************************
          ! *** PART 3: Expanding Q_{Bb,wk,is};Cc                           ***
          ! *******************************************************************

          call timer_clock(myself//'_part3',1)

          ! --------------------------------------------------------------------
          ! jd: Loop over relevant B's                                       BBB
          ! --------------------------------------------------------------------
          ! jd: This could be easily made OMP-parallel, but it doesn't help,
          !     (in fact it makes things worse), because the costly part
          !     (expand) is already internally OMP-parallelised. Similarly,
          !     parallelising here, and turning off OMP there is also counter-
          !     productive. See also comments before the loop in part 2 for why
          !     it would parallelise poorly.
          loop_B2:                                                             &
          do relev_b = 1, n_relev_b_atoms
             global_b = relev_b_atoms(relev_b)
             my_b = relev_b_atoms_my(relev_b)
             first_ngwf_idx_of_b = &
                  ngwf_indices(B_NGWF)%basis%first_on_atom(global_b)
             n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)

             if(ngwf_c == 1) n_bc_pairs_processed = n_bc_pairs_processed + 1

             ! jd: Obtain PPD indices on B (not necessarily local)
             call remote_ppd_list_of_atom(ppd_indices_on_b, &
                  n_ppds_on_b, global_b, ngwf_indices(B_NGWF)%basis, &
                  ngwf_indices(B_NGWF)%rep%ngwf_cache_handle)

             ! jd: Find PPDs common to B and C
             call ppd_set_intersection_unpacked(&
                  ppd_indices_bb_cc_intersection, n_ppds_bb_cc_intersection, & ! out
                  ppd_indices_on_b, n_ppds_on_b, ppd_indices_on_c, n_ppds_on_c)! in

             ! jd: Init expansion centres for B and C
             call swri_init_centre(centres(1), par, elements, global_b)
             call swri_init_centre(centres(2), par, elements, global_c)

             ! jd: Figure out the expansion centres, sort them, pad with -1
             !     if fewer than 2 are needed
             atoms = (/global_b, global_c/) ! jd: Avoid 'temporary created' warning
             call swri_expansion_centres(swri, &
                  expansion_centres, expansion_atoms, num_sws_in_expansion, &     ! out
                  centres, atoms, swex_quality, SWEX_BB_CC_SYMMETRIC(src_swex_h)) ! in

             ! jd: Allocate the 2-centre matrices (V, O).
             if(pub_hfx_metric == METRIC_ELECTROSTATIC) then
                allocate(vmatrix(num_sws_in_expansion,num_sws_in_expansion), &
                     stat=ierr)
                call utils_alloc_check(myself,'vmatrix',ierr)

                ! jd: Prepare the 2-centre V matrix
                call swri_build_matrix_2(vmatrix, &                              ! out
                     hfxstate%packed_metric_matrix_blocks_hts(SW_V), swri, &     ! in
                     num_sws_in_expansion, expansion_centres, expansion_atoms, & ! in
                     SW_V)
             end if
             if(pub_hfx_metric == METRIC_OVERLAP) then
                allocate(omatrix(num_sws_in_expansion,num_sws_in_expansion), &
                     stat=ierr)
                call utils_alloc_check(myself,'omatrix',ierr)

                ! jd: Prepare the 2-centre O matrix
                call swri_build_matrix_2(omatrix, &                              ! out
                     hfxstate%packed_metric_matrix_blocks_hts(SW_O), swri, &     ! in
                     num_sws_in_expansion, expansion_centres, expansion_atoms, & ! in
                     SW_O)
             end if

             ndata_in_q_one_ngwf = n_ppds_in_qs(my_b) * cell%n_pts

             ! -----------------------------------------------------------------
             ! jd: for all spins                                             sss
             ! -----------------------------------------------------------------
             spin_wk_ngwf_b_offset = 0
             loop_spins2:                                                      &
             do is = 1, pub_num_spins

                ! --------------------------------------------------------------
                ! jd: K and tcK                                               wk
                ! --------------------------------------------------------------
                loop_which_kernel2:                                            &
                do which_kernel = 1, 2

                   my_ngwf_b = hfxstate%my_bb_ngwf_offsets(my_b)

                   ! -----------------------------------------------------------
                   ! jd: for all b on B                                      bbb
                   ! -----------------------------------------------------------
                   loop_ngwf_b:                                                &
                   do ngwf_b = 1, n_ngwfs_b
                      global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

                      ! jd: Extract relevant fragment (PPDs common to B and C) of NGWF Cc
                      call ppd_extract_contents_of_subset(&
                           ngwf_cc_in_bc_ppds, &                 ! out
                           ppd_indices_bb_cc_intersection, &     ! } wanted
                           n_ppds_bb_cc_intersection, &          ! } subset
                           c_ngwfs_in_ppds(offset_to_ngwf_c:), & ! } mother
                           ppd_indices_on_c, n_ppds_on_c, &      ! } set
                           cell%n_pts)
                      if(pub_debug) then
                         call utils_sanity_check(ngwf_cc_in_bc_ppds, &
                              'ngwf_cc_in_ad_ppds',excessive=1D12)
                      end if

                      call utils_assert(.not. pub_swx_dbl_grid, myself//&
                              ': Q term calc only supported on coarse grid &
                              &(please disable swx_dbl_grid).')

                      ! *******************************************************
                      ! jd: Expand Q_Bb,wk,is;Cc in terms of SWs.
                      ! *******************************************************
                      ! swx_expansion_2_fast() calls
                      ! swx_swop_quantity_overlap_ppd_fast(), which uses OMP
                      call timer_clock(myself//'_qexpand',1)
                      if(pub_debug) then
                         call utils_sanity_check(qs_in_ppds(qs_offsets(relev_b):&
                              qs_offsets(relev_b) + int(ndata_in_q_one_ngwf-1,kind=LONG)), &
                              'qs_in_ppds', excessive=1D12)
                      end if
                      call swx_expansion_2_fast(q_coeffs, &            ! out
                           swri, expansion_centres, expansion_atoms, & ! in
                           num_sws_in_expansion, coeffs_kind, cell, &  ! in
                           vmatrix, wmatrix, omatrix, &                ! in
                           qs_in_ppds(&                                ! in
                           qs_offsets(relev_b)+int(spin_wk_ngwf_b_offset,kind=LONG):& ! in
                           qs_offsets(relev_b)+int(spin_wk_ngwf_b_offset+& ! in
                           ndata_in_q_one_ngwf-1,kind=LONG)), &            ! in
                           n_ppds_in_qs(my_b), ppd_indices_in_qs(:,my_b)) ! in
                      call timer_clock(myself//'_qexpand',2)

                      if(swx_calculate_rms_fitting_error) then
                         rms = swx_rms_fitting_error(&
                              ngwf_indices(A_NGWF)%rep%swexes(src_swex_h), &
                              ! ^ Doesn't matter which NGWF index is used, the
                              !   swex will only be used for '%quality', and
                              !   they are all identical
                              swri, n_ppds_in_qs(my_b), &
                              ppd_indices_in_qs(:,my_b), &
                              qs_in_ppds(qs_offsets(relev_b) + &
                              int(spin_wk_ngwf_b_offset,kind=LONG): &
                              qs_offsets(relev_b)+int(spin_wk_ngwf_b_offset + &
                              ndata_in_q_one_ngwf-1,kind=LONG)), &
                              expansion_atoms, expansion_centres, &
                              -1, -1, cell, coeffs_kind, src_swex_h, is, &
                              q_coeffs)

                         write(stdout,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,e10.4)') &
                              'fit_rms_error [Q]: ',&
                              global_b,':',ngwf_b,',',global_c,':',ngwf_c, &
                              ', spin ', is, ', wk ',which_kernel, ': ', rms
                      end if

#if 0
                      ! jd: Debug dump of Q's for comparison against old
                      !     implementation

                      write(filename,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                           'qcoeffs_', global_b, '_', ngwf_b, '_', global_c, '_', &
                           ngwf_c, '_', which_kernel, '_', is
                      u = utils_unit()
                      open(u, file=filename, status='new')
                      write(u,*) q_coeffs(1:num_sws_in_expansion)
                      close(u)
#endif
#if 0
                      write(filename,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                           'vdata_', global_b, '_', ngwf_b, '_', global_c, '_', &
                           ngwf_c, '_', which_kernel, '_', is
                      u = utils_unit()
                      open(u, file=filename, status='new')
                      write(u,*) expansion_atoms(:), num_sws_in_expansion, coeffs_kind
                      write(u,*) vmatrix
                      close(u)
#endif
#if 0
                      write(filename,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                           'qdata_', global_b, '_', ngwf_b, '_', global_c, '_', &
                           ngwf_c, '_', which_kernel, '_', is
                      u = utils_unit()
                      open(u, file=filename, status='new')
                       write(u,*) qs_in_ppds(&
                           qs_offsets(relev_b)+int(spin_wk_ngwf_b_offset,kind=LONG):&
                           qs_offsets(relev_b)+int(spin_wk_ngwf_b_offset+&
                           ndata_in_q_one_ngwf-1,kind=LONG))
                      close(u)
#endif

                      ! jd: Calculate potential due to the expansion of Q in
                      !     PPDs common to B and C. Uses OMP internally.
                      call timer_clock(myself//'_qexpand2',1)
                      call hfx_calc_expansion_for_q(q_expansion_in_bc_ppds, & ! out
                           ppd_indices_bb_cc_intersection, &                  ! in
                           n_ppds_bb_cc_intersection, swri, cell, &           ! in
                           expansion_atoms, expansion_centres, q_coeffs)      ! in
                      n_points_in_bb_cc_intersection = &
                           n_ppds_bb_cc_intersection * cell%n_pts
                      if(pub_debug) then
                         call utils_sanity_check(q_expansion_in_bc_ppds(&
                              1:n_points_in_bb_cc_intersection), &
                              'q_expansion_in_bc_ppds',excessive=1D12)
                      end if
                      call timer_clock(myself//'_qexpand2',2)

                      ! jd: Accumulate ngwf_cc_in_bc_ppds * q_expansion_in_bc_ppds,
                      !     transfering result to PPDs on B
                      call timer_clock(myself//'_qexpand3',1)
                      if(pub_debug) then
                         call utils_sanity_check(ngwf_cc_in_bc_ppds(&
                              1:n_points_in_bb_cc_intersection), &
                              'ngwf_cc_in_bc_ppds',excessive=1D12)
                         call utils_sanity_check(q_expansion_in_bc_ppds(&
                              1:n_points_in_bb_cc_intersection), &
                              'q_expansion_in_bc_ppds',excessive=1D12)
                      end if
                      call ppd_accumulate(&
                           kets_in_ppds_on_bs(:,which_kernel,is,my_ngwf_b), & ! adds to
                           ppd_indices_on_b, n_ppds_on_b, &      ! target PPDs
                           ngwf_cc_in_bc_ppds(&                  ! }
                           1:n_points_in_bb_cc_intersection) * & ! } data to add
                           q_expansion_in_bc_ppds(&              ! }
                           1:n_points_in_bb_cc_intersection), &  ! }
                           ppd_indices_bb_cc_intersection, &     ! src PPDs
                           n_ppds_bb_cc_intersection, &          ! src PPDs
                           .true., cell%n_pts)
                      if(pub_debug) then
                         call utils_sanity_check(&
                              kets_in_ppds_on_bs(:,which_kernel,is,my_ngwf_b),&
                              'kets_in_ppds_on_bs',excessive=1D12)
                      end if
                      call timer_clock(myself//'_qexpand3',2)

                      my_ngwf_b = my_ngwf_b + 1
                      spin_wk_ngwf_b_offset = spin_wk_ngwf_b_offset + ndata_in_q_one_ngwf

                   end do loop_ngwf_b

                end do loop_which_kernel2

             end do loop_spins2

             if(pub_hfx_metric == METRIC_OVERLAP) then
                deallocate(omatrix,stat=ierr)
                call utils_dealloc_check(myself,'omatrix',ierr)
             end if
             if(pub_hfx_metric == METRIC_ELECTROSTATIC) then
                deallocate(vmatrix,stat=ierr)
                call utils_dealloc_check(myself,'vmatrix',ierr)
             end if

          end do loop_B2

          call timer_clock(myself//'_part3',2)

          ! *******************************************************************
          ! jd: Finish with this Q set
          ! *******************************************************************
          deallocate(qs_in_ppds,stat=ierr)
          call utils_dealloc_check(myself,'qs_in_ppds',ierr)

          offset_to_ngwf_c = offset_to_ngwf_c + n_ppds_on_c * cell%n_pts

       end do loop_ngwf_c

       deallocate(relev_a_atoms,stat=ierr)
       call utils_dealloc_check(myself, 'relev_a_atoms', ierr, &
            allocated_in = 'neighbour_list_for_atom_list')

       deallocate(p_ac,stat=ierr)
       call utils_dealloc_check(myself,'p_ac',ierr)
       deallocate(ppd_indices_in_qs, stat=ierr)
       call utils_dealloc_check(myself,'ppd_indices_in_qs', ierr, &
            allocated_in = 'hfx_determine_ppd_indices_in_qs')
       deallocate(n_ppds_in_qs, stat=ierr)
       call utils_dealloc_check(myself,'n_ppds_in_qs',ierr, &
            allocated_in = 'hfx_determine_ppd_indices_in_qs')

       deallocate(dummy, stat=ierr)
       call utils_dealloc_check(myself, 'offsets_to_ngwfs', ierr, &
            allocated_in = 'remote_obtain_unpacked_ngwfs')
       deallocate(c_ngwfs_in_ppds, stat=ierr)
       call utils_dealloc_check(myself, 'ngwfs_in_ppds', ierr, &
            allocated_in = 'remote_obtain_unpacked_ngwfs')

    end do loop_C

    call utils_assert(n_bc_pairs_processed == hfxstate%n_my_bc_pairs, myself//&
         ': Bookkeeping error for BC pairs.', n_bc_pairs_processed, &
         hfxstate%n_my_bc_pairs)

    ! ####################################
    ! ### Index Cc has been eliminated ###
    ! ####################################

    ! jd: loc_grad_prefactor allows us to scale the NGWF gradient by a factor
    !     of 0.5 in conduction to account for the fact that half of the NGWFs
    !     (the valence ones) are fixed.
    if(present(grad_prefactor)) then
       loc_grad_prefactor = grad_prefactor
    else
       loc_grad_prefactor = 1.0_DP
    end if

    ! jd: Honestly have no idea of the reason for the 4 PI here.
    !     It came about when PPDs were introduced instead of TBs.
    !     The 1/2 is due to the fact that the bra-ket symmetry is broken
    !     (kets are expanded, bras are not) and the two terms in the
    !     NGWF gradient are calculated separately.
    factor = 0.5_DP / (4.0_DP * PI) * loc_grad_prefactor

    kets_in_ppds_on_bs = kets_in_ppds_on_bs * factor

    ! jd: Convert kets_in_ppds_on_bs to a HT understood by hfx_gradient().

    ! --------------------------------------------------------------------------
    ! jd: Loop over my B's                                                   BBB
    ! --------------------------------------------------------------------------
    my_ngwf_b = 1
    loop_B3:                                                                   &
    do my_b = 1, hfxstate%n_my_b_atoms
       global_b = hfxstate%my_b_atoms(my_b)
       n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)
       first_ngwf_idx_of_b = ngwf_indices(B_NGWF)%basis%first_on_atom(global_b)

       ! jd: Obtain PPD indices on B (not necessarily local)
       call remote_ppd_list_of_atom(ppd_indices_on_b, &
            n_ppds_on_b, global_b, ngwf_indices(B_NGWF)%basis, &
            ngwf_indices(B_NGWF)%rep%ngwf_cache_handle)

       ! -----------------------------------------------------------------------
       ! jd: for all b on B                                                  bbb
       ! -----------------------------------------------------------------------
       loop_ngwf_b3:                                                           &
       do ngwf_b = 1, n_ngwfs_b

          global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

          ! --------------------------------------------------------------------
          ! jd: for all spins                                                sss
          ! --------------------------------------------------------------------
          loop_spins3:                                                         &
          do is = 1, pub_num_spins

             ! -----------------------------------------------------------------
             ! jd: K and tcK
             ! -----------------------------------------------------------------
             loop_which_kernel3:                                               &
             do which_kernel = 1, 2

                ! --------------------------------------------------------------
                !                                                            PPD
                ! --------------------------------------------------------------
                offset_to_ppd = 1
                loop_ppds_on_b:                                                &
                do i_ppd = 1, n_ppds_on_b

                   cur_ppd = ppd_indices_on_b(i_ppd)

                   data_to_add(1) = real(global_bb_ngwf_idx,kind=DP)
                   data_to_add(2) = real(cur_ppd,kind=DP)
                   data_to_add(3) = real(is,kind=DP)
                   data_to_add(4) = real(which_kernel,kind=DP)
                   data_to_add(5:5+cell%n_pts-1) = kets_in_ppds_on_bs(&
                        offset_to_ppd:offset_to_ppd+cell%n_pts-1,&
                        which_kernel,is,my_ngwf_b)

                   call hash_table_add(grad_akets_dkn_my_ht, &
                        data_to_add, cell%n_pts + 4, global_bb_ngwf_idx, &
                        cur_ppd, is, which_kernel, overfill_strategy = 'F')

                   offset_to_ppd = offset_to_ppd + cell%n_pts

                end do loop_ppds_on_b

             end do loop_which_kernel3

          end do loop_spins3

          my_ngwf_b = my_ngwf_b + 1

       end do loop_ngwf_b3

    end do loop_B3

    deallocate(kets_in_ppds_on_bs, stat=ierr)
    call utils_dealloc_check(myself, 'kets_in_ppds_on_bs', ierr)

    call timer_clock(myself//'_imbal',1) ! ~15% effort
    call comms_barrier
    call timer_clock(myself//'_imbal',2)
    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hfx_gradient_term2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_calc_p_ac(p_ac, &
       s_atoms_nl, packed_dkn_cd_blocks_ht, cell_npts, global_a, global_c, &
       n_ppds_on_a, ppd_indices_on_a, p_offset, aa_ngwfs_in_ppds, &
       offset_to_aa_ngwf, dd_ngwfs_in_ppds, offsets_to_dd_ngwfs, &
       relev_d_atoms, n_relev_d_atoms, ngwf_indices, opt_effort_spent)
    !==========================================================================!
    ! Calculates the auxiliary term P_Aacs;C = \sum_Dd \phi_Aa \phi_Dd K^CcDds !
    ! for a given A and a given C. The sum over D runs over relevant D's that  !
    ! are S-neighbours with A.                                                 !
    !                                                                          !
    ! It is possible for more than one MPI rank to calculate the same P_AC.    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   p_ac (in/out): Results will be accumulated here. Caller is responsible !
    !                  for passing a suitably large buffer and zeroing it.     !
    !   s_atoms_nl (in): S-neighbour list.                                     !
    !   packed_dkn_cd_blocks_ht (in): Relevant DKN atomblocks are grabbed from !
    !                                 here.                                    !
    !   cell_npts (in): Needed for dimensioning.                               !
    !   global_a (in): The global atom index of A.                             !
    !   global_c (in): The global atom index of C.                             !
    !   n_ppds_on_a (in): Number of PPDs on atom A.                            !
    !   ppd_indices_on_a (in): Indices of PPDs on A.                           !
    !   p_offset (in): Keeps track of which part of p_ac(:) we're writing to.  !
    !   aa_ngwfs_in_ppds (in): Unpacked Aa NGWFs in PPDs prepared in advance   !
    !                          for efficiency.                                 !
    !   offset_to_aa_ngwf (in): Offset to 1st NGWF of A in the above.          !
    !   dd_ngwfs_in_ppds (in): Unpacked Dd NGWFs in PPDs prepared in advance   !
    !                          for efficiency.                                 !
    !   offsets_to_dd_ngwf (in): Offsets to 1st NGWF of D's in the above.      !
    !   relev_d_atoms (in): Global indices of relevant D atoms.                !
    !   n_relev_d_atoms (in): Number of elements in the above.                 !
    !   ngwf_indices (in): Array of HFX_NGWF_INDEX instances containing        !
    !                      pointers to the NGWF basis and the NGWF_REP used    !
    !                      by atoms A, B, C and D (only A, C and D are used).  !
    !   opt_effort_spent (out, opt): If provided, it will be set to the number !
    !                                of calls to ppd_accumulate() issued.      !
    !                                This is the number of iterations in       !
    !                                D-d-a-c-s.                                !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Keep in mind that this subroutine is called from an OMP region.        !
    !   Do not do any comms here (no barriers!), do not output (utils_trace).  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May-June 2019.                              !
    ! Adapted for OMP thread safety by Jacek Dziedzic in October 2019.         !
    ! Rehauled by Jacek Dziedzic in March-April 2020.                          !
    !==========================================================================!

    use constants, only: LONG
    use hash_table, only: HT_HASH_TABLE
    use neighbour_list, only: NL_NEIGHBOUR_LIST
    use ppd_ops, only: ppd_accumulate, ppd_product_unpacked
    use remote, only: remote_ppd_list_of_atom
    use rundat, only: pub_num_spins, pub_debug
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_sanity_check

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(inout), allocatable  :: p_ac(:)
    type(NL_NEIGHBOUR_LIST), intent(in)        :: s_atoms_nl
    type(HT_HASH_TABLE), intent(in), target    :: packed_dkn_cd_blocks_ht
    integer, intent(in)                        :: cell_npts
    integer, intent(in)                        :: global_a
    integer, intent(in)                        :: global_c
    integer, intent(in)                        :: n_ppds_on_a
    integer, intent(in)                        :: ppd_indices_on_a(:)
    integer(kind=LONG), intent(in)             :: p_offset
    real(kind=DP), intent(in)                  :: aa_ngwfs_in_ppds(:)
    integer, intent(in)                        :: offset_to_aa_ngwf
    real(kind=DP), intent(in)                  :: dd_ngwfs_in_ppds(:)
    integer, intent(in)                        :: offsets_to_dd_ngwfs(:)
    integer, intent(in)                        :: relev_d_atoms(:)
    integer, intent(in)                        :: n_relev_d_atoms
    type(HFX_NGWF_INDEX), intent(in)           :: ngwf_indices(A_NGWF:D_NGWF)
    integer, intent(out), optional             :: opt_effort_spent

    ! jd: Local variables
    integer :: ppd_indices_in_product_of_aa_dd(&
         ngwf_indices(A_NGWF)%basis%max_n_ppds_sphere) ! @IDENTICAL_RADII
    real(kind=DP) :: ppd_product_of_aa_dd(& ! @IDENTICAL_RADII
         ngwf_indices(A_NGWF)%basis%max_n_ppds_sphere * cell_npts)
    real(kind=DP) :: dkn_blocks_cd(&
         ngwf_indices(C_NGWF)%basis%max_on_atom, &
         ngwf_indices(D_NGWF)%basis%max_on_atom, pub_num_spins)
    integer :: ppd_indices_on_d(ngwf_indices(D_NGWF)%basis%max_n_ppds_sphere)
    real(kind=DP) :: dkn_cc_dd_is
    integer :: n_ppds_on_d
    integer :: n_ppds_in_product_of_aa_dd
    integer :: packed_dkn_block_size
    integer :: global_dd_ngwf_idx
    integer :: global_d
    integer :: d_idx
    integer :: relev_d
    integer :: ngwf_a, ngwf_c, ngwf_d
    integer :: n_ngwfs_a, n_ngwfs_c, n_ngwfs_d
    integer :: is
    integer(kind=LONG) :: dest_offset
    integer :: src_aa_offset, src_dd_offset
    logical :: found
    character(len=*), parameter :: myself = 'hfx_calc_p_ac'

    ! -------------------------------------------------------------------------

    ! jd: Do not trace or use barriers -- this is now called from OMP regions.

    call timer_clock(myself,1)

    if(present(opt_effort_spent)) opt_effort_spent = 0

    packed_dkn_block_size = pub_num_spins * &
         ngwf_indices(C_NGWF)%basis%max_on_atom * &
         ngwf_indices(D_NGWF)%basis%max_on_atom

    n_ngwfs_a = ngwf_indices(A_NGWF)%basis%num_on_atom(global_a)
    n_ngwfs_c = ngwf_indices(C_NGWF)%basis%num_on_atom(global_c)

    ! --------------------------------------------------------------------------
    ! jd: Loop over D's that are S-neighbours of A                           DDD
    !     This set of D's is a subset of relevant D's.
    ! --------------------------------------------------------------------------
    loop_D:                                                                    &
    do d_idx = s_atoms_nl%first_idx(global_a), s_atoms_nl%last_idx(global_a)
       global_d = s_atoms_nl%neighbours(d_idx)
       n_ngwfs_d = ngwf_indices(D_NGWF)%basis%num_on_atom(global_d)

       ! jd: Obtain PPD indices on D (not necessarily local)
       call remote_ppd_list_of_atom(ppd_indices_on_d, n_ppds_on_d, &
            global_d, ngwf_indices(D_NGWF)%basis, &
            ngwf_indices(D_NGWF)%rep%ngwf_cache_handle)

       ! jd: Xlat global_d to relev_d
       found = .false.
       do relev_d = 1, n_relev_d_atoms
          if(relev_d_atoms(relev_d) == global_d) then
             found = .true.
             exit
          end if
       end do
       call utils_assert(found, myself//': D not relevant.', global_d, d_idx)

       ! ***********************************************************************
       ! * Get all necessary DKN elements (all d, all s, all c) from HTs
       ! ***********************************************************************
       ! @room for opt: this is not directly dependent on A (but indirectly,
       ! since D=f(A)).

       ! Mind symmetry trickiness C<->D. We renamed indices, so we actually need
       ! K^{DdCc}, but we still call the array 'cd', and that's what it actually
       ! contains. So we implicitly assume K^{CD} is symmetric. So far this is
       ! fine, because this is the case both in valence and conduction
       ! calculations, and in hybrid TDDFT, where this is no longer the case,
       ! we do not calculate the NGWF gradient and never get here.

       ! jd: Get DKN^CD (all c's, all d's, spins)
       call hfx_get_dkn_atomblock(dkn_blocks_cd, & ! out
            packed_dkn_cd_blocks_ht, global_c, global_d, &
            ngwf_indices(C_NGWF)%basis%max_on_atom, &
            ngwf_indices(D_NGWF)%basis%max_on_atom)

       src_dd_offset = offsets_to_dd_ngwfs(relev_d)

       ! --------------------------------------------------------------------
       ! jd: for all d on D                                               ddd
       ! --------------------------------------------------------------------
       loop_ngwf_d:                                                            &
       do ngwf_d = 1, n_ngwfs_d
          global_dd_ngwf_idx = &
               ngwf_indices(D_NGWF)%basis%first_on_atom(global_d) + ngwf_d - 1

          dest_offset = p_offset ! jd: Points to 1st PPD on this A
          src_aa_offset = offset_to_aa_ngwf

          ! --------------------------------------------------------------------
          ! jd: for all a on A                                               aaa
          ! --------------------------------------------------------------------
          loop_ngwf_a2:                                                        &
          do ngwf_a = 1, n_ngwfs_a

             ! jd: Compute product of NGWFs Aa and Dd
             call ppd_product_unpacked(ppd_product_of_aa_dd, &
                  ppd_indices_in_product_of_aa_dd, n_ppds_in_product_of_aa_dd, &
                  aa_ngwfs_in_ppds(src_aa_offset:), ppd_indices_on_a, n_ppds_on_a, &
                  dd_ngwfs_in_ppds(src_dd_offset:), ppd_indices_on_d, n_ppds_on_d, &
                  cell_npts)
             if(pub_debug) then
                call utils_sanity_check(ppd_product_of_aa_dd, &
                     'ppd_product_of_aa_dd')
             end if

             ! -----------------------------------------------------------------
             ! jd: for all c on C                                            ccc
             ! -----------------------------------------------------------------
             loop_ngwf_c:                                                      &
             do ngwf_c = 1, n_ngwfs_c

                ! --------------------------------------------------------------
                ! jd: for all spins                                          sss
                ! --------------------------------------------------------------
                loop_spins:                                                    &
                do is = 1, pub_num_spins

                   ! jd: Retrieve K^{Cc,Dd,is} from K^{CD} obtained earlier
                   dkn_cc_dd_is = dkn_blocks_cd(ngwf_c,ngwf_d,is)

                   call ppd_accumulate(p_ac, & ! out
                        ppd_indices_on_a, n_ppds_on_a, ppd_product_of_aa_dd,&
                        ppd_indices_in_product_of_aa_dd, n_ppds_in_product_of_aa_dd, &
                        .true., cell_npts, dkn_cc_dd_is, dest_shift = dest_offset)

                   if(present(opt_effort_spent)) then
                      opt_effort_spent = opt_effort_spent + 1
                   end if

                   dest_offset = dest_offset + int(n_ppds_on_a * cell_npts,kind=LONG)

                end do loop_spins

             end do loop_ngwf_c

             src_aa_offset = src_aa_offset + n_ppds_on_a * cell_npts

          end do loop_ngwf_a2

          src_dd_offset = src_dd_offset + n_ppds_on_d * cell_npts

      end do loop_ngwf_d

    end do loop_D

    call timer_clock(myself,2)

  end subroutine hfx_calc_p_ac

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_determine_ppd_indices_in_qs(ppd_indices_in_qs, n_ppds_in_qs, &
       global_c, hfxstate, ngwf_indices, cell_nppds)
    !==========================================================================!
    ! Determines the PPD indices in a set of Q's. A set of Q's spans the       !
    ! indices Bb,wk,s for a current Cc. Atoms B are 'relevant atoms to current !
    ! C' -- my B atoms that share a PPD with C, where the BC pair is also mine.!
    ! This is explained in more detail below.                                  !
    !                                                                          !
    ! This subroutine serves as a definition of relevancy (we later examine    !
    ! the returned array to determine relevant B atoms).                       !
    !                                                                          !
    ! Returns PPD indices in the set of Qs for every B and the numbers of PPDs !
    ! (number of meaningful elements in the previous array) for every B.       !
    ! Unlike in the previous implementation, the list of PPDs is *sorted*.     !
    ! It spans more than one sphere.                                           !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   ppd_indices_in_qs (out): PPD indices in the set of Q's are returned    !
    !                            here. 1st index: ordinal of PPD, 2nd index:   !
    !                            ord_b (see 'Caveats'). The leading dimension  !
    !                            of 1st index is MAX_PPDS_OF_INTEREST.         !
    !   n_ppds_in_qs (out): Number of actual meaningful elements in each row   !
    !                       of ppd_indices_in_qs. Also indexed by ord_b.       !
    !   global_c (in): Global C atom whose S-neighbours B we are considering.  !
    !   hfxstate (in): Needed for neighbour lists and my atom lists.           !
    !   ngwf_indices (in): Only index A needed.                                !
    !   cell_nppds (in): Number of PPDs in the cell.                           !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    ! - In the returned arrays B is indexed with 'ord_b', that is, a value     !
    !   between 1 and n_my_b_atoms, that counts my B atoms, and not, as one    !
    !   might think, relev_b.                                                  !
    ! - Results are allocated here, it is someone else's responsibility to     !
    !   deallocated them when no longer needed.                                !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2020.                                 !
    !==========================================================================!

    use constants, only: garbage_int, LONG
    use neighbour_list, only: neighbour_list_are_neighbours
    use ppd_ops, only: MAX_PPDS_OF_INTEREST, ppd_set_union_unordered_fast
    use remote, only: remote_ppd_list_of_atom
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_heapsort

    implicit none

    ! jd: Arguments
    integer, intent(out), allocatable   :: ppd_indices_in_qs(:,:) ! PPD_idx, ord-my-B
    integer, intent(out), allocatable   :: n_ppds_in_qs(:) ! ord-my-B
    integer, intent(in)                 :: global_c
    type(HFX_STATE), intent(in), target :: hfxstate
    type(HFX_NGWF_INDEX), intent(in)    :: ngwf_indices(A_NGWF:D_NGWF)
    integer, intent(in)                 :: cell_nppds

    ! jd: Local variables
    integer :: ppd_indices_in_q_unsorted(MAX_PPDS_OF_INTEREST)
    integer :: indices_xlat(MAX_PPDS_OF_INTEREST)
    integer(kind=LONG) :: ppd_indices_in_q_unsorted_long(MAX_PPDS_OF_INTEREST)
    integer, allocatable :: ppd_indices_on_as(:,:)
    integer, allocatable :: n_ppds_on_as(:)
    logical :: ppd_index_lookup(cell_nppds)
    integer :: n_a_atoms
    integer :: a_idx
    integer :: global_a, global_b
    integer :: a_ord, my_b
    integer :: i_pair
    logical :: found
    integer :: ierr
    character(len=*), parameter :: myself = 'hfx_determine_ppd_indices_in_qs'

    ! -------------------------------------------------------------------------

    allocate(ppd_indices_in_qs(MAX_PPDS_OF_INTEREST,hfxstate%n_my_b_atoms), &
         stat=ierr)
    call utils_alloc_check(myself,'ppd_indices_in_qs',ierr)
    allocate(n_ppds_in_qs(hfxstate%n_my_a_atoms), stat=ierr)
    call utils_alloc_check(myself,'n_ppds_in_qs',ierr)

    n_ppds_in_qs(:) = 0
    ppd_indices_in_qs(:,:) = garbage_int
    ppd_indices_in_q_unsorted(:) = garbage_int

    ! --------------------------------------------------------------------------
    ! jd: Loop over B's given to this rank                                MY BBB
    ! --------------------------------------------------------------------------
    loop_B:                                                                    &
    do my_b = 1, hfxstate%n_my_b_atoms
       global_b = hfxstate%my_b_atoms(my_b)

       ! jd: Exclude my B's that are not S-neighbours of the current C.
       !     IOW, there's no such Q_BC. Not all my B's are relevant for this C.
       if(.not. neighbour_list_are_neighbours(&
            global_c, global_b, hfxstate%s_atoms_nl)) then
          cycle
       end if

       ! jd: Exclude my B's that, although thy are S-neighbours with C, are
       !     not in my BC-pairs. IOW, we have this B, we have this C, they
       !     are neighbours, but this B-C pair is someone else's responsibility.
       found = .false.
       do i_pair = 1, hfxstate%n_my_bc_pairs
          if(hfxstate%my_bc_pairs(1,i_pair) == global_b .and. &
               hfxstate%my_bc_pairs(2,i_pair) == global_c) then
             found = .true.
             exit
          end if
       end do
       ! *******************
       if(.not. found) cycle
       ! *******************

       n_a_atoms = hfxstate%x_atoms_nl%last_idx(global_b) - &
            hfxstate%x_atoms_nl%first_idx(global_b) + 1

       allocate(ppd_indices_on_as(&
            ngwf_indices(A_NGWF)%basis%max_n_ppds_sphere, n_a_atoms), stat=ierr)
       call utils_alloc_check(myself,'ppd_indices_on_as',ierr)
       allocate(n_ppds_on_as(n_a_atoms), stat=ierr)
       call utils_alloc_check(myself,'n_ppds_on_as',ierr)

       ppd_indices_on_as(:,:) = garbage_int
       n_ppds_on_as(:) = garbage_int
       ppd_index_lookup(:) = .false.

       ! -----------------------------------------------------------------------
       ! jd: Loop over A's that are X-neighbours of B                        AAA
       ! -----------------------------------------------------------------------
       loop_A:                                                                 &
       do a_idx = hfxstate%x_atoms_nl%first_idx(global_b), &
            hfxstate%x_atoms_nl%last_idx(global_b)
          global_a = hfxstate%x_atoms_nl%neighbours(a_idx)
          a_ord = a_idx - hfxstate%x_atoms_nl%first_idx(global_b) + 1

          ! *****************************************************
          if(.not. any(hfxstate%my_a_atoms(:) == global_a)) cycle
          ! *****************************************************

          ! jd: Obtain PPD indices on A (not necessarily local)
          call remote_ppd_list_of_atom(ppd_indices_on_as(:,a_ord), &
               n_ppds_on_as(a_ord), global_a, ngwf_indices(A_NGWF)%basis, &
               ngwf_indices(A_NGWF)%rep%ngwf_cache_handle)

          ! jd: Add these PPDs to those of Q_{Bb,wk,s;Cc}. Indices only depend
          !     {B;C}. *CAVEAT*: The indices in q are *NOT* sorted (yet) -- they
          !     are only in ascending order for each sphere, but there are
          !     multiple spheres in Q. We will sort them in a momend. A PPD
          !     index lookup is used to make this operation faster.
          call ppd_set_union_unordered_fast(&
               ppd_indices_in_q_unsorted(:), n_ppds_in_qs(my_b), &
               ppd_indices_on_as(:,a_ord), n_ppds_on_as(a_ord), &
               ppd_index_lookup)

       end do loop_A

       deallocate(n_ppds_on_as, stat=ierr)
       call utils_dealloc_check(myself,'n_ppds_on_as',ierr)
       deallocate(ppd_indices_on_as, stat=ierr)
       call utils_dealloc_check(myself,'ppd_indices_on_as',ierr)

       ppd_indices_in_q_unsorted_long(1:n_ppds_in_qs(my_b)) = &
            int(ppd_indices_in_q_unsorted(1:n_ppds_in_qs(my_b)),kind=LONG)

       ! jd: Sort the indices.
       call utils_heapsort(n_ppds_in_qs(my_b), &
            ppd_indices_in_q_unsorted_long(:), indices_xlat)
       ppd_indices_in_qs(1:n_ppds_in_qs(my_b),my_b) = &
            int(ppd_indices_in_q_unsorted_long(1:n_ppds_in_qs(my_b)))

    end do loop_B

  end subroutine hfx_determine_ppd_indices_in_qs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_accum_q_beta(qs_in_ppds, qs_offset, p_offset, &
       p_ac, global_b, global_a, global_c, ngwf_c, ppd_indices_in_qs, &
       n_ppds_in_qs, ppd_indices_on_a, n_ppds_on_a, cell, ngwf_indices, &
       hfxstate, packed_dkn_ab_blocks_hts)
    !==========================================================================!
    ! Accumulates a contribution from all ngwfs a of a single atom A to        !
    ! Q_{Bb,wk,s;Cc}, for a given atom B and a given Cc, for all ngwfs b,      !
    ! for all spins and both K and tcK.                                        !
    !                                                                          !
    ! The auxiliary term Q_{Bb,wk,s;Cc} = \sum_Aa K^BbAas P_{Aacs;C}.          !
    !                                                                          !
    ! This subroutine is threadsafe. It is called from an OMP region. Each     !
    ! thread writes to a different B, so no synchronisation primitives are     !
    ! needed.                                                                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   qs_in_ppds (inout): Results are accumulated here.                      !
    !   qs_offset (in): Keeps track of which part of the array we write to.    !
    !   p_offset (in): Keeps track of which part of P we read from.            !
    !   p_ac (in): Contains P_{Aacs;C}.                                        !
    !   global_b (in): Global index of atom B.                                 !
    !   global_a (in): Global index of atom A.                                 !
    !   global_c (in): Global index of atom C.                                 !
    !   ngwf_c (in): NGWF c on atom C.                                         !
    !   ppd_indices_in_qs (in): Indices of PPDs in a set of Q's for current C, !
    !                           for different my-B.                            !
    !   n_ppds_in_qs (in): Number of meaningful elements in the above.         !
    !   ppd_indices_on_a (in): Indices of PPD on atom A.                       !
    !   n_ppds_on_a (in): Number of elements in the above.                     !
    !   cell (in): Needed for n_pts, n_ppds.                                   !
    !   ngwf_indices (in): Array of HFX_NGWF_INDEX instances containing        !
    !                      pointers to the NGWF basis and the NGWF_REP used    !
    !                      by atoms A, B, C and D (only A, B and C used).      !
    !   hfxstate (in): Needed for the list of my B atoms.                      !
    !   packed_dkn_ab_blocks_hts (in): Distributed K^{AaBb}, tcK^{Aa}_{Bb} and !
    !                                  tcK^{Bb}_{Aa} blocks.                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in June 2019.                                  !
    !==========================================================================!

    use constants, only: LONG
    use hash_table, only: HT_HASH_TABLE
    use ppd_ops, only: ppd_accumulate
    use rundat, only: pub_num_spins, pub_debug
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_assert, &
         utils_sanity_check, utils_unit

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(inout)            :: qs_in_ppds(:) ! fused index
    integer(kind=LONG), intent(in)          :: qs_offset
    integer(kind=LONG), intent(in)          :: p_offset
    real(kind=DP), intent(in)               :: p_ac(:)
    integer, intent(in)                     :: global_b
    integer, intent(in)                     :: global_a
    integer, intent(in)                     :: global_c
    integer, intent(in)                     :: ngwf_c
    integer, intent(in), allocatable        :: ppd_indices_in_qs(:,:)
    integer, intent(in), allocatable        :: n_ppds_in_qs(:)
    integer, intent(in)                     :: ppd_indices_on_a(:)
    integer, intent(in)                     :: n_ppds_on_a
    type(CELL_INFO), intent(in)             :: cell
    type(HFX_NGWF_INDEX), intent(in)        :: ngwf_indices(A_NGWF:D_NGWF)
    type(HFX_STATE), intent(in), target     :: hfxstate
    type(HT_HASH_TABLE), intent(in), target :: packed_dkn_ab_blocks_hts(3)

    ! jd: Local variables
    real(kind=DP) :: contrib_from_this_a(n_ppds_on_a * cell%n_pts)
    real(kind=DP) :: dkn_ab(ngwf_indices(A_NGWF)%basis%max_on_atom, &
         ngwf_indices(B_NGWF)%basis%max_on_atom, pub_num_spins, 2)
    integer :: n_ppds_in_q
    integer :: my_b
    integer :: ngwf_a, ngwf_b
    integer :: n_ngwfs_a, n_ngwfs_b, n_ngwfs_c
    integer :: is
    integer :: which_kernel
    integer :: ndata_per_ngwf
    integer :: ndata_in_ppd_q
    integer :: saved_start_ppd
    integer(kind=LONG) :: q_offset
    integer(kind=LONG) :: src_offset
    logical :: found
    character(len=*), parameter :: myself = 'hfx_accum_q_beta'

    ! -------------------------------------------------------------------------

    ! jd: Do not add a barrier here, this routine is called a different number
    !     of times from different ranks, possibly even zero times.
    call utils_trace_in(myself)
    call timer_clock(myself,1)

    n_ngwfs_a = ngwf_indices(A_NGWF)%basis%num_on_atom(global_a)
    n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)
    n_ngwfs_c = ngwf_indices(C_NGWF)%basis%num_on_atom(global_c)

    ! jd: Xlat global_b to my_b (and so, an index in n_ppds_in_qs).
    found = .false.
    do my_b = 1, hfxstate%n_my_b_atoms
       if(hfxstate%my_b_atoms(my_b) == global_b) then
          found = .true.
          exit
       end if
    end do
    call utils_assert(found,myself//': Not my B.', global_b)

    n_ppds_in_q = n_ppds_in_qs(my_b)
    ndata_per_ngwf = n_ppds_in_q*cell%n_pts
    ndata_in_ppd_q = ndata_per_ngwf * n_ngwfs_b * 2 * pub_num_spins

    ! jd: Get DKN^AB (all a's, all b's, spins)
    do which_kernel = 1, 2
       call hfx_get_dkn_atomblock(dkn_ab(:,:,:,which_kernel), & ! out
            packed_dkn_ab_blocks_hts(which_kernel), global_a, global_b, &
            ngwf_indices(A_NGWF)%basis%max_on_atom, &
            ngwf_indices(B_NGWF)%basis%max_on_atom)

    end do
    dkn_ab(:,:,:,2) = -dkn_ab(:,:,:,2)
                    ! ^ equivalent to using -tcK in old code

    q_offset = qs_offset
    saved_start_ppd = 1

    ! -----------------------------------------------------------------------
    ! jd: For all spins                                                   sss
    ! -----------------------------------------------------------------------
    loop_spins1:                                                            &
    do is = 1, pub_num_spins

       ! -----------------------------------------------------------------
       ! jd: for both kernels                                        K/tcK
       ! ------------------------------------------------------------------
       loop_which_kernel:                                                &
       do which_kernel = 1, 2

          ! --------------------------------------------------------------------
          ! jd: for all b on B                                               bbb
          ! --------------------------------------------------------------------
          loop_ngwf_b:                                                         &
          do ngwf_b = 1, n_ngwfs_b

             contrib_from_this_a(:) = 0.0_DP
             do ngwf_a = 1, n_ngwfs_a

                ! jd: p_offset points to first PPD of this A. Offset to current
                !     a-c-is.
                src_offset = p_offset + int(&
                     ((ngwf_a-1) * n_ngwfs_c * pub_num_spins + &
                     (ngwf_c-1) * pub_num_spins + &
                     (is-1)) * n_ppds_on_a * cell%n_pts,kind=LONG)

                contrib_from_this_a(:) = contrib_from_this_a(:) + &
                     dkn_ab(ngwf_a,ngwf_b,is,which_kernel) * &
                     p_ac(src_offset:src_offset+int(n_ppds_on_a*cell%n_pts-1,kind=LONG))
             end do

             if(pub_debug) then
                call utils_sanity_check(contrib_from_this_a(&
                     1:n_ppds_on_a*cell%n_pts), 'contrib_from_this_a', &
                     excessive=1D12)
             end if

             call ppd_accumulate(qs_in_ppds, ppd_indices_in_qs(:,my_b), &
                  n_ppds_in_qs(my_b), contrib_from_this_a(:), &
                  ppd_indices_on_a, n_ppds_on_a, .true., cell%n_pts, &
                  sum_start_ppd = saved_start_ppd, dest_shift = q_offset)

             q_offset = q_offset + int(ndata_per_ngwf,kind=LONG)

          end do loop_ngwf_b

       end do loop_which_kernel

    end do loop_spins1

#if 0
    ! jd: Useful when debugging
    write(filename,'(a,i0,a,i0,a,i0,a,i0)') &
         'k_', global_a, '_', global_b, '_', global_c, '_', ngwf_c
    u = utils_unit()
    open(u, file=filename)
    write(u,*) dkn_ab(:,:,:,:)
    close(u)
#endif

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine hfx_accum_q_beta

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_populate_swop_cache(hfxstate, swri, par, mdl, ireg, &
       ngwf_indices, swex_h)
    !==========================================================================!
    ! On each proc, predicts which SWOPs will be needed in the calculation and !
    ! how many times. Ranks the SWOPs in the order of usefulness, from those   !
    ! that will be used most to those that will be used least. Precalculates   !
    ! all the SWOPs that fit in the cache, honouring pub_cache_limit_for_swops.!
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (in/out): Needed for B-C pairs, HT units, among other things. !
    !   swri (in/out): Contains the SWOP cache.                                !
    !   par (in): Needed for orig_atom, proc_of_atom.                          !
    !   mdl (in): Needed for elements, cell.                                   !
    !   ireg (in): Embedding region ID.                                        !
    !   ngwf_indices (in): Array of HFX_NGWF_INDEX instances containing        !
    !                      pointers to the NGWF basis and the NGWF_REP used    !
    !                      by atoms A, B, C and D.                             !
    !   swex_h (in): Handle to the SW_EX in use. Needed e.g. for SW quality.   !
    !--------------------------------------------------------------------------!
    ! More info/caveats:                                                       !
    !   This subroutine goes through a dry-run of an energy calculation, and,  !
    !   if it is deemed that NGWF gradient will be calculated, the gradient    !
    !   calculation. Determining if we're going to do the NGWF gradient is     !
    !   tricky, we do this in rundat (pub_ngwf_gradient_needed). There are     !
    !   three stages for the dry-run:                                          !
    !   1) Expanding my B-C pairs, this happens once per NGWF iteration.       !
    !   2) Generating potentials from B and C on all A. This, in principle,    !
    !      happens in LNV iterations, but the assumption is that because       !
    !      expansions are cached, this will only happen once per NGWF iter.    !
    !      This assumption will no longer be true if we permit running in      !
    !      RAM-starved scenarios where expansions are re-calculated, and it    !
    !      might shift the balance (SWOPs needed for this stage will become    !
    !      more important).                                                    !
    !   3) (gradient only) 3A: Expanding Q,                                    !
    !                      3B: Generating the potentials from Q's expansion.   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in August 2019.                                !
    ! Revised by Jacek Dziedzic for rehauled gradient calc in April 2020.      !
    !==========================================================================!

    use comms, only: pub_total_num_procs, pub_on_root, comms_allgather, &
         comms_barrier, comms_reduce
    use constants, only: garbage_int, stdout, NORMAL, SWEX_BB_CC_SYMMETRIC, &
         METRIC_ELECTROSTATIC, METRIC_OVERLAP, VERBOSE
    use hash_table, only: hash_table_init, hash_table_free, hash_table_list, &
         hash_table_lookup_ptr_nocount
    use model_type, only: MODEL
    use neighbour_list, only: neighbour_list_for_atom_list
    use parallel_strategy, only: PARAL_INFO
    use ppd_ops, only: ppd_set_intersection_unpacked
    use remote, only: remote_ppd_list_of_atom, remote_n_ppds_of_atom, &
         remote_ngwf_cache_hts
    use rundat, only: pub_ngwf_gradient_needed, pub_hfx_output_detail, &
         pub_hfx_metric, pub_num_spins, pub_cache_limit_for_swops, &
         pub_threads_max
    use sw_resolution_of_identity, only: SW_RI, ATOM_CENTRE, SW_QUALITY, &
         swri_init_centre, swri_qualities_equal, swri_obtain_swops_in_ppd_ptr
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_dealloc_check, utils_trace_in, &
         utils_trace_out, utils_flush

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout)   :: hfxstate
    type(SW_RI), intent(inout)       :: swri
    type(PARAL_INFO), intent(in)     :: par
    type(MODEL), intent(in)          :: mdl
    integer, intent(in)              :: ireg
    type(HFX_NGWF_INDEX), intent(in) :: ngwf_indices(A_NGWF:D_NGWF)
    integer, intent(in)              :: swex_h

    ! jd: Local variables
    type(HT_HASH_TABLE), target :: swop_hits_ht
    integer, allocatable :: a_ppd_list(:)
    integer, allocatable :: swop_hits(:,:)
    integer, allocatable :: ppd_indices_in_qs(:,:)
    integer, allocatable :: n_ppds_in_qs(:)
    integer, allocatable :: relev_a_atoms(:) ! relev_a -> global_a
    real(kind=DP), pointer :: swops_in_ppd_ptr(:)
    real(kind=DP), pointer :: dlist_real_ptr(:)
    integer :: n_swops_possible_on_procs(0:pub_total_num_procs-1)
    integer :: n_swops_needed_on_procs(0:pub_total_num_procs-1)
    integer :: pub_cache_limit_for_swops_on_procs(0:pub_total_num_procs-1)
    real(kind=DP) :: n_hits_stored_on_procs(0:pub_total_num_procs-1)
    real(kind=DP) :: n_hits_total_on_procs(0:pub_total_num_procs-1)
    integer :: ppd_indices_on_a(ngwf_indices(A_NGWF)%basis%max_n_ppds_sphere)
    integer :: ppd_indices_on_b(ngwf_indices(B_NGWF)%basis%max_n_ppds_sphere)
    integer :: ppd_indices_on_c(ngwf_indices(C_NGWF)%basis%max_n_ppds_sphere)
    integer :: ppd_indices(max(&
         ngwf_indices(B_NGWF)%basis%max_n_ppds_sphere, &
         ngwf_indices(C_NGWF)%basis%max_n_ppds_sphere))
    integer :: ppd_indices_bb_cc_intersection(min(&
         ngwf_indices(B_NGWF)%basis%max_n_ppds_sphere, &
         ngwf_indices(C_NGWF)%basis%max_n_ppds_sphere))
    integer :: relev_b_atoms(hfxstate%n_my_b_atoms) ! relev_b -> global_b
    integer :: relev_b_atoms_my(hfxstate%n_my_b_atoms) ! relev_b -> my_b
    integer :: n_ppds_on_as(hfxstate%n_my_a_atoms) ! indexed by relev_a
    integer :: dlist(max_num_at_per_ppd)

    type(ATOM_CENTRE) :: src_centre
    type(SW_QUALITY) :: swex_quality
    integer :: n_ppds_on_a, n_ppds_on_b, n_ppds_on_c
    integer :: n_a_ppd_list
    integer :: n_common_ppds
    integer :: n_hits
    integer :: pair_idx
    integer :: global_a, global_b, global_c, global_d
    integer :: prev_global_b
    integer :: orig_global_b, orig_global_c
    integer :: n_ngwfs_a, n_ngwfs_b, n_ngwfs_c, n_ngwfs_d
    integer :: b_idx, c_idx
    integer :: relev_a, relev_b
    integer :: my_a, my_b
    integer :: n_relev_a_atoms, n_relev_b_atoms
    integer :: n_ppds_in_qset
    integer :: n_ppds_in_q_from_this_b
    integer :: n_points_in_qset
    integer :: n_ppds_in_p
    integer :: n_ppds_in_p_from_this_a
    integer :: n_points_in_p
    integer :: n_ppds_bb_cc_intersection
    integer :: n_swops
    integer :: proc
    integer :: elem
    integer :: src_atom
    integer :: ppd_idx
    integer :: n_sws
    integer :: n_points_in_f_max, n_points_in_p_max, n_points_in_qset_max
    integer :: i_ppd
    integer :: d_idx
    integer :: num_d
    integer :: cur_ppd
    integer :: fused_index_range
    integer :: n_records_in_grad_dkn_my_term1, n_records_in_grad_dkn_my_term2
    integer :: n_records_in_grad_dkn_my
    integer :: n_hits_here
    real(kind=DP) :: n_hits_total
    real(kind=DP) :: n_hits_stored
    integer :: ierr
!$  integer, external :: omp_get_thread_num
    integer :: tid
    logical :: bb_cc_symmetry
    real(kind=DP) :: fraction_cacheable
    real(kind=DP) :: effort_saved
    character :: sw_or_swpot
    character(len=26) :: buffer
    character(len=*), parameter :: myself = 'hfx_populate_swop_cache'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    call timer_clock(myself//'_bookkeeping',1)

    call utils_assert(.not. hfxstate%swop_cache_populated, &
         myself//': SWOP cache already populated.')

    call hash_table_init(swop_hits_ht, 'SWOP_HITS', 2, maxslots_swop_hits, &
         max_swop_hits, hfxstate%hash_table_info_unit)

    bb_cc_symmetry = SWEX_BB_CC_SYMMETRIC(swex_h)

    call utils_assert(swri_qualities_equal(&
         ngwf_indices(A_NGWF)%rep%swexes(swex_h)%quality, &
         ngwf_indices(B_NGWF)%rep%swexes(swex_h)%quality), myself//&
         ': SW_EX qualities do not match between two NGWF bases (A and B).')
    call utils_assert(swri_qualities_equal(&
         ngwf_indices(A_NGWF)%rep%swexes(swex_h)%quality, &
         ngwf_indices(C_NGWF)%rep%swexes(swex_h)%quality), myself//&
         ': SW_EX qualities do not match between two NGWF bases (A and C).')
    call utils_assert(swri_qualities_equal(&
         ngwf_indices(A_NGWF)%rep%swexes(swex_h)%quality, &
         ngwf_indices(D_NGWF)%rep%swexes(swex_h)%quality), myself//&
         ': SW_EX qualities do not match between two NGWF bases (A and D).')

    swex_quality = ngwf_indices(A_NGWF)%rep%swexes(swex_h)%quality
    ! --------------------------------------------------------------------------
    ! This is a dry-run version of stage 1 (expansion). We emulate this chain: !
    ! - hf_exchange_dkn_indep_stage()                                          !
    !   - hfx_expand_ngwf_pairs()                                              !
    !     - swx_expand_my_ngwf_pairs()                 S T A G E  1            !
    !       - swx_my_ngwf_pair_fitting()                                       !
    !         - swx_expansion_2_fast()                                         !
    !           - swx_swop_quantity_overlap_ppd_fast()                         !
    !             - swri_obtain_swops_in_ppd_ptr_fast()                        !
    ! --------------------------------------------------------------------------
    ! jd: Loop over all my BC pairs                                           BC
    ! --------------------------------------------------------------------------
    prev_global_b = garbage_int
    loop_BC1:                                                                  &
    do pair_idx = 1, hfxstate%n_my_bc_pairs
       global_b = hfxstate%my_bc_pairs(1,pair_idx)
       global_c = hfxstate%my_bc_pairs(2,pair_idx)
       orig_global_b = par%orig_atom(global_b)
       orig_global_c = par%orig_atom(global_c)
       n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)
       n_ngwfs_c = ngwf_indices(C_NGWF)%basis%num_on_atom(global_c)

       ! jd: Ignore atoms not represented in this SWRI
       if(.not. mdl%regions(ireg)%elements(orig_global_b)%&
            in_swri(hfxstate%hfx_swri_h)) cycle
       if(.not. mdl%regions(ireg)%elements(orig_global_c)%&
            in_swri(hfxstate%hfx_swri_h)) cycle

       ! jd: Get the PPD list for atom B (in general, not sparse-local)
       !     Do not repeat unnecessarily if B did not change.
       if(global_b /= prev_global_b) then
          call remote_ppd_list_of_atom(ppd_indices_on_b, n_ppds_on_b, global_b,&
               ngwf_indices(B_NGWF)%basis, &
               ngwf_indices(B_NGWF)%rep%ngwf_cache_handle)
       end if

       ! jd: If there is Bb-Cc symmetry, then it will be exploited later to
       !     elide re-doing expansions where C>B, *but only* if we also own
       !     the pair (C,B). If the pair (C,B) is not in the list of my pairs,
       !     someone else gets the mirror expansion and it is effectively
       !     redone. This is the reason for the more complicated check here.
       if(bb_cc_symmetry .and. .not. global_b <= global_c) then
          if(internal_is_pair_mine(global_c, global_b, &
               hfxstate%my_bc_pairs, hfxstate%n_my_bc_pairs)) cycle
       end if

       ! jd: Get the PPD list for atom C (in general, not sparse-local)
       call remote_ppd_list_of_atom(ppd_indices_on_c, n_ppds_on_c, global_c,&
            ngwf_indices(C_NGWF)%basis, &
            ngwf_indices(C_NGWF)%rep%ngwf_cache_handle)

       ! jd: Determine the PPD indices in the intersection
       call ppd_set_intersection_unpacked(ppd_indices, n_common_ppds, &   ! out
            ppd_indices_on_b , n_ppds_on_b, ppd_indices_on_c, n_ppds_on_c)

       n_hits = n_ngwfs_b * n_ngwfs_c

       ! jd: If symmetry in effect, for same-centre only include the lower
       !     triangle and diagonal.
       if(bb_cc_symmetry .and. global_b == global_c) then
          n_hits = (n_ngwfs_b + 1) * n_ngwfs_b / 2
       end if

       call internal_hit_swop(ppd_indices, n_common_ppds, n_hits, global_b)
       if(global_c /= global_b) then
          call internal_hit_swop(ppd_indices, n_common_ppds, n_hits, global_c)
       end if

       prev_global_b = global_b

    end do loop_BC1

    ! --------------------------------------------------------------------------
    ! This is a dry-run version of stage 2 (generation of potentials from the  !
    ! expansion. We emulate this chain:                                        !
    ! - hf_exchange_calculate()                                                !
    !   - hfx_main()                                                           !
    !     - hfx_calc_expansions_bb_cc()                                        !
    !       - swri_obtain_swops_in_ppd_ptr_ptr_fast()  S T A G E  2            !
    ! --------------------------------------------------------------------------
    ! jd: Loop over all my BC pairs                                           BC
    ! --------------------------------------------------------------------------
    prev_global_b = garbage_int
    loop_BC2:                                                                  &
    do pair_idx = 1, hfxstate%n_my_bc_pairs
       global_b = hfxstate%my_bc_pairs(1,pair_idx)
       global_c = hfxstate%my_bc_pairs(2,pair_idx)
       orig_global_b = par%orig_atom(global_b)
       orig_global_c = par%orig_atom(global_c)
       n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)
       n_ngwfs_c = ngwf_indices(C_NGWF)%basis%num_on_atom(global_c)

       ! jd: Ignore atoms not represented in this SWRI
       if(.not. mdl%regions(ireg)%elements(orig_global_b)%&
            in_swri(hfxstate%hfx_swri_h)) cycle
       if(.not. mdl%regions(ireg)%elements(orig_global_c)%&
            in_swri(hfxstate%hfx_swri_h)) cycle

       if(global_b /= prev_global_b) then
          ! jd: Deallocate previous PPD list, if any
          if(allocated(a_ppd_list)) then
             deallocate(a_ppd_list,stat=ierr)
             call utils_dealloc_check(myself,'a_ppd_list',ierr, &
                  allocated_in='hfx_determine_a_ppd_list_for_b')
          end if
          call hfx_determine_a_ppd_list_for_b(a_ppd_list, n_a_ppd_list, &
               remote_ngwf_cache_hts, hfxstate, ngwf_indices(A_NGWF), global_b)
       end if

       ! jd: If there is Bb-Cc symmetry, then it will be exploited later to
       !     elide re-doing expansions where C>B, *but only* if we also own
       !     the pair (C,B). If the pair (C,B) is not in the list of my pairs,
       !     someone else gets the mirror expansion and it is effectively
       !     redone. This is the reason for the more complicated check here.
       if(bb_cc_symmetry .and. .not. global_b <= global_c) then
          if(internal_is_pair_mine(global_c, global_b, &
               hfxstate%my_bc_pairs, hfxstate%n_my_bc_pairs)) cycle
       end if

       n_hits = n_ngwfs_b * n_ngwfs_c

       ! jd: If symmetry in effect, for same-centre only include the lower
       !     triangle and diagonal.
       if(bb_cc_symmetry .and. global_b == global_c) then
          n_hits = (n_ngwfs_b + 1) * n_ngwfs_b / 2
       end if

       call internal_hit_swop(a_ppd_list, n_a_ppd_list, n_hits, global_b)
       if(global_b /= global_c) then
          call internal_hit_swop(a_ppd_list, n_a_ppd_list, n_hits, global_c)
       end if

       prev_global_b = global_b

    end do loop_BC2

    ! jd: Deallocate final PPD list, if any (corner case of BC list empty)
    if(allocated(a_ppd_list)) then
       deallocate(a_ppd_list,stat=ierr)
       call utils_dealloc_check(myself,'a_ppd_list',ierr, &
            allocated_in='hfx_determine_a_ppd_list_for_b')
    end if

    if(pub_ngwf_gradient_needed) then
       ! -----------------------------------------------------------------------
       ! This is a dry-run version of stage 3. We do                           !
       ! - 3A (gradient term 2, expanding Q) and                               !
       ! - 3B (gradient term 2, potential due to expanded Q) simultaneously    !
       ! We emulate this chain:                                                !
       ! - hf_exchange_calculate()                                             !
       !   - hfx_main()                                                        !
       !     - hfx_gradient_term2()                                            !
       !       - swx_expansion_2_fast()                                        !
       !           - swx_swop_quantity_overlap_ppd_fast()                      !
       !             - swri_obtain_swops_in_ppd_ptr_fast()   S T A G E  3 A    !
       !       - hfx_calc_expansion_for_q()                                    !
       !         - swri_obtain_swops_in_ppd_ptr_fast()       S T A G E  3 B    !
       ! -----------------------------------------------------------------------

       n_points_in_p_max = 0
       n_points_in_qset_max = 0

       ! -----------------------------------------------------------------------
       ! jd: Loop over my C atoms                                            CCC
       ! -----------------------------------------------------------------------
       loop_C:                                                                 &
       do c_idx = 1, hfxstate%n_my_c_atoms
          global_c = hfxstate%my_c_atoms(c_idx)
          n_ngwfs_c = ngwf_indices(C_NGWF)%basis%num_on_atom(global_c)

          ! jd: This section duplicates some 60 lines of hfx_gradient_term2(),
          !     but putting this in a 19-param subroutine would be even worse.

          ! jd: Obtain PPD indices on C (not necessarily local)
          call remote_ppd_list_of_atom(ppd_indices_on_c, n_ppds_on_c, &
               global_c, ngwf_indices(C_NGWF)%basis, &
               ngwf_indices(C_NGWF)%rep%ngwf_cache_handle)

          ! jd: Figure out PPDs that go into this set of Q
          call hfx_determine_ppd_indices_in_qs(ppd_indices_in_qs, n_ppds_in_qs,&
               global_c, hfxstate, ngwf_indices, mdl%cell%n_ppds)

          ! jd: Figure out relevant B atoms -- those, who have >0 PPDs and feature
          !     in our BC pair list, so B atoms that are present in this Q set.
          n_relev_b_atoms = 0
          relev_b_atoms(:) = garbage_int
          relev_b_atoms_my(:) = garbage_int
          do my_b = 1, hfxstate%n_my_b_atoms
             if(n_ppds_in_qs(my_b) > 0) then
                n_relev_b_atoms = n_relev_b_atoms + 1
                global_b = hfxstate%my_b_atoms(my_b)
                relev_b_atoms(n_relev_b_atoms) = global_b
                relev_b_atoms_my(n_relev_b_atoms) = my_b
             end if
          end do

          ! jd: Figure out the dimensions of the Q set -- to estimate memory use
          n_ppds_in_qset = 0
          n_points_in_qset = 0
          do relev_b = 1, n_relev_b_atoms
             global_b = relev_b_atoms(relev_b)
             my_b = relev_b_atoms_my(relev_b)
             n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)
             n_ppds_in_q_from_this_b = &
                  n_ppds_in_qs(my_b) * n_ngwfs_b * 2 * pub_num_spins
             n_ppds_in_qset = n_ppds_in_qset + n_ppds_in_q_from_this_b
             n_points_in_qset = n_points_in_qset + &
                  n_ppds_in_q_from_this_b * mdl%cell%n_pts
          end do

          n_points_in_qset_max =  max(n_points_in_qset_max, n_points_in_qset)

          ! jd: Figure out relevant A's -- s-X-neighbours of C that are also
          !     my-A atoms. The 2nd condition ([*]) is needed to exclude A's that
          !     are C's s-X-neighbours but only via a B that is not my B. These are
          !     not our responsibility.
          !     Determine the dimension in P, offsets to different A's.
          n_ppds_in_p = 0
          n_points_in_p = 0
          n_ppds_on_as(:) = garbage_int
          call neighbour_list_for_atom_list(relev_a_atoms, n_relev_a_atoms, &
               relev_b_atoms, n_relev_b_atoms, hfxstate%x_atoms_nl, par%nat)

          do relev_a = 1, n_relev_a_atoms
             global_a = relev_a_atoms(relev_a)
             n_ngwfs_a = ngwf_indices(A_NGWF)%basis%num_on_atom(global_a)

             ! jd: Obtain # of PPDs on A (not necessarily local)
             call remote_n_ppds_of_atom(n_ppds_on_as(relev_a), global_a, &
                  ngwf_indices(A_NGWF)%basis, &
                  ngwf_indices(A_NGWF)%rep%ngwf_cache_handle)

             n_ppds_in_p_from_this_a = n_ppds_on_as(relev_a) * &
                  n_ngwfs_a * n_ngwfs_c * pub_num_spins
             n_ppds_in_p = n_ppds_in_p + n_ppds_in_p_from_this_a
             n_points_in_p = &
                  n_points_in_p + n_ppds_in_p_from_this_a * mdl%cell%n_pts
          end do

          n_points_in_p_max =  max(n_points_in_p_max, n_points_in_p)

          deallocate(relev_a_atoms,stat=ierr)
          call utils_dealloc_check(myself, 'relev_a_atoms', ierr, &
               allocated_in = 'neighbour_list_for_atom_list')

          ! --------------------------------------------------------------------
          ! jd: Loop over relevant B's                                       BBB
          ! --------------------------------------------------------------------
          loop_B:                                                              &
          do relev_b = 1, n_relev_b_atoms
             global_b = relev_b_atoms(relev_b)
             my_b = relev_b_atoms_my(relev_b)
             n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)

             ! jd: Obtain PPD indices on B (not necessarily local)
             call remote_ppd_list_of_atom(ppd_indices_on_b, &
                  n_ppds_on_b, global_b, ngwf_indices(B_NGWF)%basis, &
                  ngwf_indices(B_NGWF)%rep%ngwf_cache_handle)

             ! jd: Find PPDs common to B and C
             call ppd_set_intersection_unpacked(ppd_indices_bb_cc_intersection,&
                  n_ppds_bb_cc_intersection, ppd_indices_on_b, n_ppds_on_b, &
                  ppd_indices_on_c, n_ppds_on_c)

             n_hits = n_ngwfs_b * n_ngwfs_c * pub_num_spins * 2 ! 2 is for wk

             ! --- 3A ---
             call internal_hit_swop(ppd_indices_in_qs(:,my_b), &
                  n_ppds_in_qs(my_b), n_hits, global_b)
             if(global_b /= global_c) then
             call internal_hit_swop(ppd_indices_in_qs(:,my_b), &
                  n_ppds_in_qs(my_b), n_hits, global_c)
             end if

             ! --- 3B ---
             call internal_hit_swop(ppd_indices_bb_cc_intersection, &
                  n_ppds_bb_cc_intersection, n_hits, global_b)
             if(global_b /= global_c) then
                call internal_hit_swop(ppd_indices_bb_cc_intersection, &
                     n_ppds_bb_cc_intersection, n_hits, global_c)
             end if

          end do loop_B

          deallocate(ppd_indices_in_qs, stat=ierr)
          call utils_dealloc_check(myself,'ppd_indices_in_qs', ierr, &
               allocated_in = 'hfx_determine_ppd_indices_in_qs')
          deallocate(n_ppds_in_qs, stat=ierr)
          call utils_dealloc_check(myself,'n_ppds_in_qs',ierr, &
               allocated_in = 'hfx_determine_ppd_indices_in_qs')

       end do loop_C

       call comms_reduce('MAX',n_points_in_p_max)
       call comms_reduce('MAX',n_points_in_qset_max)

    else
       n_points_in_p_max = 0
       n_points_in_qset_max = 0
    end if ! NGWF gradient needed?

    ! -----------------------------------------------------------------------
    ! Now we'd like to find out how big the most important temporaries are  !
    ! going to be. These do not involve SWOPs. We do a dry-run version of   !
    ! - hf_exchange_calculate()                                             !
    !   - hfx_main()                                                        !
    !     - hfx_contract_over_dkn_cc()                                      !
    ! -----------------------------------------------------------------------

    n_points_in_f_max = 0
    ! --------------------------------------------------------------------------
    ! jd: Loop over B's given to this rank                                MY BBB
    ! --------------------------------------------------------------------------
    loop_B2:                                                                   &
    do b_idx = 1, hfxstate%n_my_b_atoms
       global_b = hfxstate%my_b_atoms(b_idx)
       n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)

       ! jd: Deallocate previous PPD list, if any
       if(allocated(a_ppd_list)) then
          deallocate(a_ppd_list,stat=ierr)
          call utils_dealloc_check(myself,'a_ppd_list',ierr, &
               allocated_in='hfx_determine_a_ppd_list_for_b')
       end if
       call hfx_determine_a_ppd_list_for_b(a_ppd_list, n_a_ppd_list, &
            remote_ngwf_cache_hts, hfxstate, ngwf_indices(A_NGWF), global_b)

       ! jd: Mimic relevant part of hfx_contract_over_dkn_cc()
       fused_index_range = 0
       ! -----------------------------------------------------------------------
       ! jd: All my A PPD *relevant to current B*                          r PPD
       ! -----------------------------------------------------------------------
       loop_relev_PPD:                                                          &
       do i_ppd = 1, n_a_ppd_list
          cur_ppd = a_ppd_list(i_ppd)

          ! jd: Retrieve list of all D's sharing this PPD
          call hash_table_lookup_ptr_nocount(dlist_real_ptr, num_d, &
               hfxstate%dlists_ht, global_b, cur_ppd)
          call utils_assert(num_d <= max_num_at_per_ppd, myself//&
               ': Insufficient max_num_at_per_ppd [1]', &
               max_num_at_per_ppd, num_d)
          call utils_assert(num_d /= -1, myself//': dlist not found [1]', &
               global_b, cur_ppd)
          call utils_assert(num_d /= 0, myself//': dlist empty [1]', &
               global_b, cur_ppd)

          ! jd: Convert to ints, pad
          dlist(1:num_d) = nint(dlist_real_ptr(1:num_d))
          dlist(num_d+1:) = garbage_int

          ! --------------------------------------------------------------------
          ! jd: For all D's that are relevant to this PPD                  r DDD
          ! --------------------------------------------------------------------
          loop_relev_D:                                                        &
          do d_idx = 1, num_d
             global_d = dlist(d_idx)

             ! jd: Careful. There are corner cases where D shares a PPD with
             !     an A that is an X-neighbour of B, but A and D *do not*
             !     intersect (from the sparsity point of view). For every
             !     candidate D we must double-check if it's an x-S-neighbour of B.
             if(.not. any(hfxstate%xs_atoms_nl%neighbours(&
                  hfxstate%xs_atoms_nl%first_idx(global_b):&
                  hfxstate%xs_atoms_nl%last_idx(global_b)) == global_d)) then
                cycle
             end if

             n_ngwfs_d = ngwf_indices(D_NGWF)%basis%num_on_atom(global_d)
             fused_index_range = fused_index_range + n_ngwfs_d * mdl%cell%n_pts

          end do loop_relev_D

       end do loop_relev_PPD
       n_points_in_f_max = &
            max(n_points_in_f_max,fused_index_range * n_ngwfs_b * pub_num_spins)

    end do loop_B2

    ! jd: Deallocate final PPD list, if any (corner case of B list empty)
    if(allocated(a_ppd_list)) then
       deallocate(a_ppd_list,stat=ierr)
       call utils_dealloc_check(myself,'a_ppd_list',ierr, &
            allocated_in='hfx_determine_a_ppd_list_for_b')
    end if

    call comms_reduce('MAX',n_points_in_f_max)

    ! -----------------------------------------------------------------------
    ! Now we'd like to find out how big the most important temporaries are  !
    ! going to be for the NGWF gradient part, if any. These do not involve  !
    ! SWOPs. We do a dry-run version of                                     !
    ! For term 1:                                                           !
    ! - hf_exchange_calculate()                                             !
    !   - hfx_main()                                                        !
    !     - hfx_grad_ket_dkn_baccum()                                       !
    ! For term 2:                                                           !
    ! - hf_exchange_calculate()                                             !
    !   - hfx_main()                                                        !
    !     - hfx_gradient_term2()                                            !
    !                                                                       !
    ! The logic is practically identical, except term 2 goes over all my B  !
    ! atoms, while term 1 goes over all my A atoms. We need to pick the     !
    ! larger of the two values, since the HT gets purged between term 1 and !
    ! term 2, so they are never used simultaneously. Of course there is     !
    ! always going to be more my-A atoms than my-B atoms, but in the        !
    ! interest of future-proofing the code we still calculate both, in case !
    ! NGWF radii or bases ever differ between A and B.                      !
    ! -----------------------------------------------------------------------

    if(pub_ngwf_gradient_needed) then

       n_records_in_grad_dkn_my_term1 = 0
       n_records_in_grad_dkn_my_term2 = 0

       ! ********* TERM1 *********

       ! -----------------------------------------------------------------------
       ! jd: Loop over my A's                                                AAA
       ! -----------------------------------------------------------------------
       loop_A3:                                                                &
       do my_a = 1, hfxstate%n_my_a_atoms
          global_a = hfxstate%my_a_atoms(my_a)
          n_ngwfs_a = ngwf_indices(A_NGWF)%basis%num_on_atom(global_a)

          ! jd: Obtain PPD indices on A (not necessarily local)
          call remote_ppd_list_of_atom(ppd_indices_on_a, &
               n_ppds_on_a, global_a, ngwf_indices(A_NGWF)%basis, &
               ngwf_indices(A_NGWF)%rep%ngwf_cache_handle)

          ! jd: Emulate loop over NGWFs a, spins, which_kernel, PPDs on A.
          !     Each innermost iteration adds a record of npts+4 doubles.
          !     These are summed over my A, but unreduced across MPI ranks.
          !     We are only counting records, not points, to avoid LONG ints
          !     (we need comms_reduce() later, it doesn't work for LONGs).
          n_records_in_grad_dkn_my_term1 = n_records_in_grad_dkn_my_term1 + &
               n_ngwfs_a * n_ppds_on_a * pub_num_spins * 2

       end do loop_A3

       ! ********* TERM2 *********

       ! -----------------------------------------------------------------------
       ! jd: Loop over my B's                                                BBB
       ! -----------------------------------------------------------------------
       loop_B3:                                                                 &
       do my_b = 1, hfxstate%n_my_b_atoms
          global_b = hfxstate%my_b_atoms(my_b)
          n_ngwfs_b = ngwf_indices(B_NGWF)%basis%num_on_atom(global_b)

          ! jd: Obtain PPD indices on B (not necessarily local)
          call remote_ppd_list_of_atom(ppd_indices_on_b, &
               n_ppds_on_b, global_b, ngwf_indices(B_NGWF)%basis, &
               ngwf_indices(B_NGWF)%rep%ngwf_cache_handle)

          ! jd: Emulate loop over NGWFs b, spins, which_kernel, PPDs on B.
          !     Each innermost iteration adds a record of npts+4 doubles.
          !     These are summed over my B, but unreduced across MPI ranks.
          !     We are only counting records, not points, to avoid LONG ints
          !     (we need comms_reduce() later, it doesn't work for LONGs).
          n_records_in_grad_dkn_my_term2 = n_records_in_grad_dkn_my_term2 + &
               n_ngwfs_b * n_ppds_on_b * pub_num_spins * 2

       end do loop_B3

       n_records_in_grad_dkn_my = &
            max(n_records_in_grad_dkn_my_term1, n_records_in_grad_dkn_my_term2)

    end if

    call timer_clock(myself//'_bookkeeping',2)
    call timer_clock(myself//'_sort',1)

    call internal_crude_sort(swop_hits_ht, swop_hits, n_swops, n_hits_total)
    call hash_table_list(swop_hits_ht,0)
    call hash_table_free(swop_hits_ht)

    call timer_clock(myself//'_sort',2)
    call timer_clock(myself//'_calc',1)

    ! jd: Barrier, or else imbalances show up in 'calc' because of
    !     the comms_allgather() below.
    call timer_clock(myself//'_imbalance',1)
    call comms_barrier
    call timer_clock(myself//'_imbalance',2)

    call internal_memory_use_estimate_and_adjust

    ! --------------------------------------------------------------------------
    ! jd: Dry-run for the actual generation of cached SWOPS so that we can
    !     show the user feedback first.
    ! --------------------------------------------------------------------------
    n_hits_stored = 0.0_DP
    do elem = 1, min(swri%swops_in_ppds_ht%max_elements,n_swops)
       n_hits_here = swop_hits(1,elem)
       n_hits_stored = n_hits_stored + n_hits_here
    end do

    ! jd: User feedback on how we're doing with the caching
    if(pub_hfx_output_detail >= NORMAL) then
       call comms_allgather(n_swops_possible_on_procs, &
            swri%swops_in_ppds_ht%max_elements, 1, gather_not_allgather = .true.)
       call comms_allgather(n_swops_needed_on_procs, &
            n_swops, 1, gather_not_allgather = .true.)
       call comms_allgather(n_hits_stored_on_procs, &
            n_hits_stored, 1, gather_not_allgather = .true.)
       call comms_allgather(n_hits_total_on_procs, &
            n_hits_total, 1, gather_not_allgather = .true.)
       call comms_allgather(pub_cache_limit_for_swops_on_procs, &
            pub_cache_limit_for_swops, 1, gather_not_allgather = .true.)
       if(pub_on_root) then
          write(stdout,'(a)') '+'//repeat('-',78)//'+'
          write(stdout,'(a)') '|    MPI |                           |          &
               &    |           |              |'
          write(stdout,'(a)') '|   rank |       SWOP cache capacity | SWOPs nee&
               &ded | Cacheable |    Hit ratio |'
          write(stdout,'(a)') '+'//repeat('-',78)//'+'
          do proc = 0, pub_total_num_procs-1
             if(n_swops_needed_on_procs(proc) /= 0) then
                fraction_cacheable = min(&
                     real(n_swops_possible_on_procs(proc),kind=DP) / &
                     real(n_swops_needed_on_procs(proc),kind=DP),1.0_DP) * &
                     100.0_DP
             else
                ! jd: Corner case of zero swops needed, avoid division by zero
                fraction_cacheable = 100.0_DP
             end if
             if(n_hits_total_on_procs(proc) /= 0) then
                effort_saved = n_hits_stored_on_procs(proc) / &
                     n_hits_total_on_procs(proc) * 100.0_DP
             else
                ! jd: Corner case of zero swops needed, avoid division by zero
                effort_saved = 100.0_DP
             end if
             write(buffer,'(i0,a,i0,a)') &
                  pub_cache_limit_for_swops_on_procs(proc), ' MiB (', &
                  n_swops_possible_on_procs(proc), ' el)'
             write(stdout,'(a,i7,a,a26,a,i10,a,f9.2,a,f12.2,a)') '|', proc,' |', &
                  adjustr(trim(buffer)), ' |', &
                  n_swops_needed_on_procs(proc), ' el |', fraction_cacheable, &
                  '% |', effort_saved, '% |'
          end do
          write(stdout,'(a)') '+'//repeat('-',78)//'+'
       end if
    end if

    if(pub_on_root .and. pub_hfx_output_detail >= VERBOSE) then
       write(stdout,'(a)') 'HFx: - Populating SWOP cache...'
       call utils_flush(stdout, .true.)
    end if

    if(pub_hfx_metric == METRIC_ELECTROSTATIC) sw_or_swpot = 'P'
    if(pub_hfx_metric == METRIC_OVERLAP) sw_or_swpot = 'S'

    ! jd: Starting from the most heavily used SWOPs, get them into cache.
    !     Stop when cache limit is hit, or when we add all we needed.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    n_hits_stored = 0
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(STATIC)  &
!$OMP PRIVATE(elem, n_hits_here, src_atom, ppd_idx, swops_in_ppd_ptr, n_sws,   &
!$OMP      src_centre, tid)                                                    &
!$OMP SHARED(swri, n_swops, swop_hits, par, mdl, sw_or_swpot, swex_quality,    &
!$OMP      ireg)                                                               &
!$OMP REDUCTION(+:n_hits_stored)
     do elem = 1, min(swri%swops_in_ppds_ht%max_elements,n_swops)
       tid = 0
!$     tid = omp_get_thread_num()

       n_hits_here = swop_hits(1,elem)
       src_atom = swop_hits(2,elem)
       ppd_idx = swop_hits(3,elem)
       call swri_init_centre(src_centre, par, mdl%regions(ireg)%elements, src_atom)

       call swri_obtain_swops_in_ppd_ptr(swri, swops_in_ppd_ptr, n_sws, &
            src_centre, src_atom, ppd_idx, sw_or_swpot, mdl%cell, swex_quality,&
            dont_lookup = .true., use_fast = .true.)

       if(tid == 0 .and. mod(elem,10) == 0) then
          call hash_table_list(swri%swops_in_ppds_ht,0)
       end if

       n_hits_stored = n_hits_stored + n_hits_here
    end do
!$OMP END PARALLEL DO
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    deallocate(swop_hits,stat=ierr)
    call utils_dealloc_check(myself,'swop_hits',ierr, &
         allocated_in='hfx_populate_swop_cache::internal_crude_sort')

    hfxstate%swop_cache_populated = .true.
    hfxstate%all_swops_cached = (n_hits_stored == n_hits_total)

    call timer_clock(myself//'_calc',2)

    call timer_clock(myself//'_imbalance',1)
    call comms_barrier
    call timer_clock(myself//'_imbalance',2)
    call timer_clock(myself,2)
    call utils_trace_out(myself)

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    logical function internal_is_pair_mine(global_b, global_c, my_bc_pairs, &
         n_my_bc_pairs)
      !========================================================================!
      ! Returns .true. if pair B-C is in the list of my BC pairs.              !
      ! Uses a linear scan, but this is not expected to be problematic.        !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   global_b (in): Global index of B atom.                               !
      !   global_c (in): Global index of C atom.                               !
      !   my_bc_pairs (in): Array of my B-C pairs.                             !
      !   n_my_bc_pairs (in): Number of entries in the above.                  !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in August 2019.                              !
      !========================================================================!

      implicit none

      ! jd: Arguments
      integer, intent(in) :: global_b
      integer, intent(in) :: global_c
      integer, intent(in) :: my_bc_pairs(:,:)
      integer, intent(in) :: n_my_bc_pairs

      ! jd: Local variables
      integer :: i_pair

      ! -----------------------------------------------------------------------

      do i_pair = 1, n_my_bc_pairs
         if(my_bc_pairs(1,i_pair) == global_b .and. &
              my_bc_pairs(2,i_pair) == global_c) then
            internal_is_pair_mine = .true.
            return
         end if
      end do

      internal_is_pair_mine = .false.
      return

    end function internal_is_pair_mine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_hit_swop(ppd_indices, n_ppds, n_hits, src_atom)
      !========================================================================!
      ! Increases the number of hits for SWOPs originating on src_atom needed  !
      ! in PPDs from ppd_indices. That is, we use swops_hits_ht to simply store!
      ! usage counts (hits) for particular SWOP cache entries (src_atom-PPD    !
      ! pairs).                                                                !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   ppd_indices (in): Array of PPD indices where the SWOPs are needed.   !
      !   n_ppds (in): The number of elements in the above.                    !
      !   n_hits (in): Number of hits to increase the count by.                !
      !   src_atom (in): Global index of the atom on which the SWOP is centred.!
      !                                                                        !
      ! Side effects:                                                          !
      !   Modifies swop_hits_ht in parent.                                     !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in August 2019.                              !
      !========================================================================!

      use hash_table, only: hash_table_axpy

      implicit none

      ! jd: Arguments
      integer, intent(in) :: ppd_indices(:)
      integer, intent(in) :: n_ppds
      integer, intent(in) :: n_hits
      integer, intent(in) :: src_atom

      ! jd: Local variables
      integer :: i_ppd
      integer :: cur_ppd

      ! -------------------------------------------------------------------------

      do i_ppd = 1, n_ppds
         cur_ppd = ppd_indices(i_ppd)
         call hash_table_axpy(swop_hits_ht, (/real(n_hits,kind=DP)/), 1, &
              cur_ppd, src_atom, overfill_strategy = 'F')
      end do

    end subroutine internal_hit_swop

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_crude_sort(swop_hits_ht, swop_hits, n_swops, n_hits_total)
      !========================================================================!
      ! Converts an unsorted hash table of SWOP cache entry hits, so a map     !
      ! (src_atom, PPD) -> n_hits, into a "mostly sorted" simple array of      !
      ! triples: (n_hits, src_atom, PPD). "Mostly sorted" means that n_hits    !
      ! will, in general, be decreasing for subsequent elements, but           !
      ! not strictly. Instead there will be groups of elements with similar    !
      ! magnitude of n_hits, followed by groups of elements with a smaller     !
      ! magnitude, in decrements of 0.75 (factor). This is done for the sake of!
      ! efficiency -- an exact sort would require O(N^2) or fancy coding.      !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   swop_hits_ht (in): The input HT containing number of hits (counts)   !
      !                      for SWOP cache entries.                           !
      !   swop_hits (out): Array of triples (n_hits, src_atom, PPD), "mostly   !
      !                    sorted". Allocated here.                            !
      !   n_swops (out): Number of elements in the above.                      !
      !   n_hits_total (out): Sum of the 1st index in the above.               !
      !------------------------------------------------------------------------!
      ! Caveats:                                                               !
      !   swop_hits is allocated here, caller is responsible for deallocation. !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in August 2019.                              !
      !========================================================================!

      use hash_table, only: HT_HASH_TABLE, HT_STORAGE_NODE, hash_table_iterate
      use utils, only: utils_alloc_check, utils_assert

      implicit none

      ! jd: Arguments
      type(HT_HASH_TABLE), intent(in), target :: swop_hits_ht
      integer, intent(out), allocatable       :: swop_hits(:,:)
      integer, intent(out)                    :: n_swops
      real(kind=DP), intent(out)              :: n_hits_total

      ! jd: Local variables
      type(HT_STORAGE_NODE), pointer :: ptr
      integer :: slot
      integer :: cur_ppd, src_atom
      integer :: ndata
      integer :: key3, key4, key5
      logical :: eot
      integer :: n_elements, n_elements_done
      integer :: n_hits, n_hits_max
      integer :: n_hits_up, n_hits_dn
      real(kind=DP) :: n_hits_real(1)
      integer :: ierr
      real(kind=DP), pointer :: dummy_data_ptr(:)
      real(kind=DP), parameter :: factor = 0.75_DP               ! [*]
      character(len=*), parameter :: myself = &
           'hfx_populate_swop_cache::internal_crude_sort'
      ! [*] 'factor' determines the relative magnitude differences between
      !     subsequent groups of elements in the output. For the default value
      !     of 0.75 this means that the first group of elements will contain
      !     the largest n_hits and all n_hits not smaller than 0.75 of the
      !     largest n_hits. The second group will contain elements smaller than
      !     0.75 of largest n_hits, but not smaller than 0.75**2 of largest
      !     n_hits, and so on. Increasing 'factor' towards 1.0 will tend towards
      !     a perfect sort at the cost of increased number of iterate scans of
      !     the HT. Decreasing 'factor' towards 0.0 will reduce the quality of
      !     the sort, but decrease the time spent sorting.

      ! -----------------------------------------------------------------------

      ! jd: Allocate swop_hits to the requisite size
      n_elements = swop_hits_ht%n_slots + swop_hits_ht%n_chained
      allocate(swop_hits(3,n_elements),stat=ierr)
      call utils_alloc_check(myself,'swop_hits',ierr)

      ! jd: Find maximum value in the HT
      n_hits_max = -1
      ptr => NULL()
      slot = 0
      eot = .false.
      do while (.not. eot)
         call hash_table_iterate(ndata, ptr, slot, swop_hits_ht, &
              cur_ppd, src_atom, key3, key4, key5, eot, dummy_data_ptr, .true.,&
              n_hits_real)
         if(.not. eot) then
            n_hits = nint(n_hits_real(1))
            n_hits_max = max(n_hits_max,n_hits)
         end if
      end do

!      call utils_assert(n_hits_max /= -1, myself//': SWOP hits table empty.')
!      ^ SWOP hits table can legitimately be empty in embedding when some ranks
!        only get pairs from outside the SWRI, so don't uncomment that.

      ! jd: Iterate many times through the HT, picking elements from the largest
      !     to the smallest in ranges where the lower bound is smaller than the
      !     upper bound by 'factor'.
      n_elements_done = 0
      n_hits_up = n_hits_max
      do while (n_elements_done < n_elements)
         n_hits_dn = int(real(n_hits_up,kind=DP) * factor)
         ptr => NULL()
         slot = 0
         eot = .false.
         do while (.not. eot)
            call hash_table_iterate(ndata, ptr, slot, swop_hits_ht, &
                 cur_ppd, src_atom, key3, key4, key5, eot, dummy_data_ptr, &
                 .true., n_hits_real)
            if(.not. eot) then
               n_hits = nint(n_hits_real(1))
               if(n_hits > n_hits_dn .and. n_hits <= n_hits_up) then
                  n_elements_done = n_elements_done + 1
                  call utils_assert(n_elements_done <= n_elements, myself//&
                       ': Range check error.', n_elements_done, n_elements)
                  swop_hits(1,n_elements_done) = n_hits
                  swop_hits(2,n_elements_done) = src_atom
                  swop_hits(3,n_elements_done) = cur_ppd
               end if
            end if
         end do

         n_hits_up = n_hits_dn

      end do

      n_swops = n_elements
      n_hits_total = sum(real(swop_hits(1,:),kind=DP))

    end subroutine internal_crude_sort

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_memory_use_estimate_and_adjust
      !========================================================================!
      ! Calculates and prints out a memory use estimate for HFx's TEFCI engine.!
      ! If 'hfx_memory_limit' is specified, adjustments to cache limits are    !
      ! made.                                                                  !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   None.                                                                !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in April 2020.                               !
      !========================================================================!

      use comms, only: pub_my_proc_id, pub_on_root, pub_root_proc_id, &
           comms_reduce, comms_bcast
      use constants, only: SW_V, SW_O, SIZEOF_DOUBLE, BRIEF, SAFE_DIV_EPS, &
           LONG, garbage_real
      use hash_table, only: node_size, hash_table_mem_usage, hash_table_init, &
           hash_table_calc_cache_size
      use remote, only: remote_ngwf_cache_hts
      use rundat, only: pub_cache_limit_for_expansions, pub_hfx_output_detail, &
           pub_cache_limit_for_prods, pub_cache_limit_for_swops, pub_task, &
           pub_hfx_memory_limit, pub_hfx_memory_weights, &
           pub_hfx_bessel_rad_nptsx
      use utils, only: utils_report_memory_estimate, utils_int_to_str, &
           utils_real_to_str, utils_abort

      implicit none

      ! jd: Local variables
      real(kind=DP) :: mem_swops_mib
      real(kind=DP) :: mem_prods_mib
      real(kind=DP) :: mem_expansions_mib
      real(kind=DP) :: mem_ngwfs_dd_mib
      real(kind=DP) :: mem_ngwfs_mib, mem_ngwfs_tot_mib
      real(kind=DP) :: mem_dlists_mib
      real(kind=DP) :: mem_vmetric_mib
      real(kind=DP) :: mem_ometric_mib
      real(kind=DP) :: minimum_required_memory_mib
      real(kind=DP) :: mem_to_spare_mib
      real(kind=DP) :: dummy
      integer(kind=LONG) :: mem_radial_bessel_long
      integer(kind=LONG) :: mem_swops_long
      integer(kind=LONG) :: mem_prods_long
      integer(kind=LONG) :: mem_expansions_long
      integer(kind=LONG) :: mem_ngwfs_dd_long
      integer(kind=LONG) :: mem_dlists_long
      integer(kind=LONG) :: mem_ngwfs_long
      integer(kind=LONG) :: mem_vmetric_long
      integer(kind=LONG) :: mem_ometric_long
      integer(kind=LONG) :: overhead_from_slots, overhead_from_records
      integer(kind=LONG) :: n_points_in_grad_dkn_loc
      integer(kind=LONG) :: n_bytes_in_grad_kets_dkn_my
      integer(kind=LONG) :: minimum_required_memory
      integer(kind=LONG) :: user_adjustable_memory
      real(kind=DP)      :: n_bytes_in_grad_dkn_loc_real
      real(kind=DP)      :: n_bytes_in_dkn_cd
      real(kind=DP)      :: n_bytes_in_dkn_ab
      real(kind=DP)      :: n_bytes_in_coeffs
      real(kind=DP)      :: coeffs_overhead_mib
      real(kind=DP)      :: scaling_frac
      integer :: max_cached_expansions
      integer :: max_cached_swops_in_ppds
      integer :: n_max_records_in_grad_dkn_my
      integer :: ngwf_cache_handles(4)
      integer :: unique(4)
      integer :: n_unique
      integer :: i_unique
      integer :: ngwf_index
      integer :: c_idx
      integer :: global_c
      integer :: n_d_atoms, n_a_atoms
      integer :: n_cd_blocks, n_ab_blocks
      integer :: max_cached_dkn_blocks_cd, max_cached_dkn_blocks_ab
      integer :: dkn_blocks_cd_overhead, dkn_blocks_ab_overhead
      integer :: n_my_bb_cc_ngwf_pairs
      integer :: n_swops_max
      integer :: num_sws
      integer :: actual_n_swops, actual_n_swops_max
      real(kind=DP) :: adjusted_cache_weights(3)
      real(kind=DP) :: cache_weights(3)
      real(kind=DP) :: weight_sum
      integer :: auto_mem_limit_for_swops
      integer :: auto_mem_limit_for_expansions
      integer :: auto_mem_limit_for_prods

      ! Used to force promotion of ints in memory estimate i/o/2 avoid wraparound
      integer(kind=LONG), parameter :: L = 1_LONG

      ! jd: Default memory weights for valence calculations:
      !     - mostly SWOPs and expansions.
      real(kind=DP), parameter :: default_hfx_memory_weights_val(3) = &
           [10.0,5.0,1.0]
      ! jd: Default memory weights for conduction calculations:
      !     - Only SWOPs. Expansions are not re-used (LNV iterations do not
      !       invoke HFx, so NGWFs change between every call to HFx). There are
      !       so many products that we wind up caching a tiny minority and
      !       lookup costs exceed the gains. Better not to cache at all.
      real(kind=DP), parameter :: default_hfx_memory_weights_cond(3) = &
           [1.0,0.0,0.0]

      character(len=*), parameter :: myself = &
           'hfx_populate_swop_cache::internal_memory_use_estimate_and_adjust'

      ! -----------------------------------------------------------------------

      ! jd: Radial Bessel lookup in SWRI
      mem_radial_bessel_long = pub_hfx_bessel_rad_nptsx * swri%max_num_bessels * &
           SIZEOF_DOUBLE

      ! jd: Dd NGWFs HT
      call hash_table_mem_usage(hfxstate%ppd_ngwfs_dd_ht, mem_ngwfs_dd_mib)
      call comms_reduce('MAX', mem_ngwfs_dd_mib)
      mem_ngwfs_dd_long = mem_ngwfs_dd_mib * 1024.0_DP * 1024.0_DP

      ! jd: Remote NGWFs HT
      ! jd: Find unique cache handles so that we don't count the same basis more
      !     than once, but work correctly in the presence of multiple bases.
      do ngwf_index = A_NGWF, D_NGWF
         ngwf_cache_handles(ngwf_index) = &
              ngwf_indices(ngwf_index)%rep%ngwf_cache_handle
      end do

      unique(:) = -1
      unique(1) = ngwf_cache_handles(A_NGWF)
      n_unique = 1
      do ngwf_index = B_NGWF, D_NGWF
         if(all(ngwf_cache_handles(ngwf_index) /= unique(1:n_unique))) then
            n_unique = n_unique + 1
            unique(n_unique) = ngwf_cache_handles(ngwf_index)
         end if
      end do

      mem_ngwfs_tot_mib = 0
      do i_unique = 1, n_unique
         call hash_table_mem_usage(remote_ngwf_cache_hts(unique(i_unique)), &
              mem_ngwfs_mib)
         mem_ngwfs_tot_mib = mem_ngwfs_tot_mib + mem_ngwfs_mib
      end do
      call comms_reduce('MAX', mem_ngwfs_tot_mib)
      mem_ngwfs_long = mem_ngwfs_tot_mib * 1024.0_DP * 1024.0_DP

      ! jd: DLISTS HT
      call hash_table_mem_usage(hfxstate%dlists_ht, mem_dlists_mib)
      call comms_reduce('MAX', mem_dlists_mib)
      mem_dlists_long = mem_dlists_mib * 1024.0_DP * 1024.0_DP

      ! jd: [VO]METRIC HT
      if(swri%has_metric(SW_V)) then
         call hash_table_mem_usage(hfxstate%packed_metric_matrix_blocks_hts(1),&
              mem_vmetric_mib)
      else
         mem_vmetric_mib = 0.0_DP
      end if
      if(swri%has_metric(SW_O)) then
         call hash_table_mem_usage(hfxstate%packed_metric_matrix_blocks_hts(2),&
              mem_ometric_mib)
      else
         mem_ometric_mib = 0.0_DP
      end if
      call comms_reduce('MAX', mem_vmetric_mib)
      call comms_reduce('MAX', mem_ometric_mib)
      mem_vmetric_long = mem_vmetric_mib * 1024.0_DP * 1024.0_DP
      mem_ometric_long = mem_ometric_mib * 1024.0_DP * 1024.0_DP

      ! jd: KETS_DKN_MY HT
      ! jd: First, account for the HT overhead. The HT itself is unavailable
      !     (it will only be constructed in hfx_main()), so we do this manually.
      !     We slightly overestimate the cost, by assuming there is overhead
      !     from all slots (just overhead, not data in slots), but all the
      !     records are stored in chains. It's hard to be more accurate at this
      !     stage.
      overhead_from_slots = &
           L*maxslots_grad_kets_dkn_my_per_atom * par%nat * node_size
      overhead_from_records = L*n_records_in_grad_dkn_my * node_size

      ! jd: Find a max over all procs, but keep the individual values
      n_max_records_in_grad_dkn_my = n_records_in_grad_dkn_my
      call comms_reduce('MAX',n_max_records_in_grad_dkn_my)

      ! jd: Get the total -- for the max value
      n_bytes_in_grad_kets_dkn_my = L*n_max_records_in_grad_dkn_my * &
           (mdl%cell%n_pts+4) * SIZEOF_DOUBLE + &
           overhead_from_slots + overhead_from_records

      ! jd: KETS_DKN_LOC HT
      ! jd: Estimate by multiplying the value for kets_dkn_my-term 2 by the
      !     factor n_local_a_atoms/n_my_b_atoms. This is done on proc-unreduced
      !     values, since the factor varies across procs. Then we reduce.
      scaling_frac = real(par%num_atoms_on_proc(pub_my_proc_id),kind=DP) / &
           real(max(1,hfxstate%n_my_b_atoms),kind=DP)
      overhead_from_slots = &
           L*maxslots_grad_kets_dkn_loc_per_atom * par%nat * node_size
      overhead_from_records = int(real(n_records_in_grad_dkn_my_term2,kind=DP)*&
           scaling_frac,kind=LONG) * node_size

      n_bytes_in_grad_dkn_loc_real = &
           real(n_records_in_grad_dkn_my_term2,kind=DP) * &
           (mdl%cell%n_pts+4) * SIZEOF_DOUBLE * scaling_frac + &
           overhead_from_slots + overhead_from_records
      call comms_reduce('MAX',n_bytes_in_grad_dkn_loc_real)

      ! jd: DKN blocks CD and AB.
      !     First determine HT overheads.
      call hfx_dkn_blocks_ht_size(max_cached_dkn_blocks_cd, &
           max_cached_dkn_blocks_ab, ngwf_indices, par%nat)
      dkn_blocks_cd_overhead = real(max_cached_dkn_blocks_cd * node_size,kind=DP)
      dkn_blocks_ab_overhead = real(max_cached_dkn_blocks_ab * node_size,kind=DP)

      ! jd: Count CD blocks
      n_cd_blocks = 0
      do c_idx = 1, hfxstate%n_my_c_atoms
         global_c = hfxstate%my_c_atoms(c_idx)
         n_d_atoms = hfxstate%sxs_atoms_nl%n_neighbours(global_c)
         n_cd_blocks = n_cd_blocks + n_d_atoms
      end do

      n_bytes_in_dkn_cd = real(n_cd_blocks,kind=DP) * pub_num_spins * &
           ngwf_indices(C_NGWF)%basis%max_on_atom * &
           ngwf_indices(D_NGWF)%basis%max_on_atom * SIZEOF_DOUBLE + &
           dkn_blocks_cd_overhead
      call comms_reduce('MAX',n_bytes_in_dkn_cd)

      ! jd: Count AB blocks
      n_ab_blocks = 0
      do b_idx = 1, hfxstate%n_my_b_atoms
         global_b = hfxstate%my_b_atoms(b_idx)
         n_a_atoms = hfxstate%x_atoms_nl%n_neighbours(global_b)
         n_ab_blocks = n_ab_blocks + n_a_atoms
      end do

      n_bytes_in_dkn_ab = real(n_ab_blocks,kind=DP) * pub_num_spins * &
           ngwf_indices(A_NGWF)%basis%max_on_atom * &
           ngwf_indices(B_NGWF)%basis%max_on_atom * SIZEOF_DOUBLE + &
           dkn_blocks_ab_overhead

      call comms_reduce('MAX',n_bytes_in_dkn_ab)

      ! jd: COEFFS HT. This is an estimate.
      !     We count my Bb-Cc pairs, multiply by number of coeffs,
      !     and take into account node_size bytes of overhead.
      !     We add overhead from the HT procs.
      !     Simplifications:
      !     - ignores the fact that some pairs are on-site (smaller matrix),
      !       which leads to a slight overestimate.
      !     - ignores Bb-Cc symmetry. Taking this into account is tricker than
      !       just dividing by 2 -- we only have distributed fragments of the
      !       lower triangle. IOW, the total number of stored coeffs is well-
      !       approximated by half the number of total coeffs (modulo the dia-
      !       gonal terms), but we only get sub-rectangles of this. This leads
      !       to a non-negligible overestimate, but so be it.
      !     - ignores the fact that many coeffs will wind up in procs and thus
      !       double-counts their overhead (slight overestimate).
      n_my_bb_cc_ngwf_pairs = sum(hfxstate%my_bc_pairs(3,:))
      call hash_table_mem_usage(ngwf_indices(B_NGWF)%rep%swexes(swex_h)%&
           coeffs_ht, dummy, mem_overhead = coeffs_overhead_mib)
      n_bytes_in_coeffs = real(n_my_bb_cc_ngwf_pairs,kind=DP) * &
           (2.0_DP * swri%quality%num_sws_per_centre * SIZEOF_DOUBLE + &
           node_size)
      n_bytes_in_coeffs = n_bytes_in_coeffs + &
           coeffs_overhead_mib * 1024.0_DP * 1024.0_DP

      call comms_reduce('MAX',n_bytes_in_coeffs)

      ! jd: Do not restrict this output to higher output details.
      !     It doesn't only print, it also calculates the required sum.
      if(pub_hfx_output_detail >= BRIEF) then
         call utils_report_memory_estimate(&
              'HFx TEFCI engine: minimum requirements',&
              [&
              'Radial Bessel lookup              ',&
              'Dd NGWFs hash table               ',&
              'All remote NGWFs hash table       ',&
              'dlists hash table                 ',&
              'coeffs hash table (estimate)      ',&
              'V metric matrix hash table        ',&
              'O metric matrix hash table        ',&
              'f auxiliary term                  ',&
              'P term in NGWF gradient           ',&
              'Q term in NGWF gradient           ',&
              'My kets in NGWF grad. (estimate)  ',&
              'Local kets in NGWF grad. (estim.) ',&
              'K^{CD} hash table                 ', &
              'K^{AB} hash table                 ', &
              'tcK^A_B hash table                ', &
              'tcK^B_A hash table                '&
              ], &
              [&
              ! Radial Bessel lookup in SWRI
              mem_radial_bessel_long, &
              ! PPD NGWFs Dd hash table
              mem_ngwfs_dd_long, &
              ! Remote NGWFs hash table
              mem_ngwfs_long, &
              ! dlists hash table
              mem_dlists_long, &
              ! coeffs hash table
              int(n_bytes_in_coeffs,kind=LONG), &
              ! V metric hash table
              merge(mem_vmetric_long, 0_LONG,swri%has_metric(SW_V)),&
              ! O metric hash table
              merge(mem_ometric_long, 0_LONG,swri%has_metric(SW_O)),&
              ! f auxiliary term in hfx_contract_over_dkn_cc()
              L*n_points_in_f_max*SIZEOF_DOUBLE, &
              ! P_AC term in gradient term2
              merge(L*n_points_in_p_max*SIZEOF_DOUBLE,0_LONG,&
              pub_ngwf_gradient_needed), &
              ! Q auxiliary term in gradient term2
              merge(L*n_points_in_qset_max*SIZEOF_DOUBLE,0_LONG,&
              pub_ngwf_gradient_needed), &
              ! grad_kets_dkn_my HT
              merge(n_bytes_in_grad_kets_dkn_my,0_LONG,&
              pub_ngwf_gradient_needed), &
              ! grad_kets_dkn_loc HT
              merge(int(n_bytes_in_grad_dkn_loc_real,kind=LONG),0_LONG,&
              pub_ngwf_gradient_needed), &
              ! DKN^CD
              int(n_bytes_in_dkn_cd,kind=LONG), &
              ! DKN^AB
              merge(int(n_bytes_in_dkn_ab,kind=LONG), 0_LONG, &
              pub_ngwf_gradient_needed), &
              ! tcDKN^AB
              merge(int(n_bytes_in_dkn_ab,kind=LONG),0_LONG, &
              pub_ngwf_gradient_needed), &
              ! tcDKN^BA
              merge(int(n_bytes_in_dkn_ab,kind=LONG),0_LONG, &
              pub_ngwf_gradient_needed) &
              ], minimum_required_memory)
      end if

      ! ***********************************************
      ! *** Adjust user-controlled part if asked to ***
      ! ***********************************************

      if(pub_hfx_memory_limit /= -1) then

         if(pub_on_root) then
            ! jd: Check to see if this is even doable, i.e. if we haven't
            !     exceeded the maximum with the parts we have no control of
            minimum_required_memory_mib = &
                 real(minimum_required_memory,kind=DP) / 1024.0_DP / 1024.0_DP
            if(pub_hfx_memory_limit < minimum_required_memory_mib) then
               call utils_abort(myself//&
                    ': Specified maximum memory limit for HFx ('//&
                    trim(utils_int_to_str(pub_hfx_memory_limit))//&
                    ' MB per MPI rank) cannot be satisfied. The minimum memory &
                    &requirement of this system is '//trim(utils_real_to_str(&
                    minimum_required_memory_mib,'f16.1'))//' MB per MPI rank.')
            end if

            ! jd: Choose suitable default for hfx_memory_weights,
            !     if not specified
            if(all(pub_hfx_memory_weights == -1.0_DP)) then
               if(pub_task == 'COND' .or. pub_task == 'PROPERTIES_COND') then
                  pub_hfx_memory_weights = default_hfx_memory_weights_cond
               else
                  pub_hfx_memory_weights = default_hfx_memory_weights_val
               end if
            end if

            adjusted_cache_weights(:) = pub_hfx_memory_weights(:)

            ! jd: If expansions and/or products not in use (cache limit set to
            !     zero explicitly), do not take them into account.
            if(.not. hfxstate%expansions_ht_in_use) then
               adjusted_cache_weights(2) = 0.0_DP
            end if
            if(.not. hfxstate%ppd_products_ht_in_use) then
               adjusted_cache_weights(3) = 0.0_DP
            end if

            ! jd: Normalise weights so that they add up to 1.
            weight_sum = sum(adjusted_cache_weights)
            call utils_assert(all(adjusted_cache_weights >= 0.0_DP), &
                 myself//': cache weights must all be non-negative.')
            call utils_assert(weight_sum > SAFE_DIV_EPS, &
                 myself//': cache weights cannot all be zero.')
            cache_weights(:) = adjusted_cache_weights(:) / weight_sum

            if(pub_hfx_output_detail >= VERBOSE) then
               write(stdout,'(a)') 'HFx: - Adjusting cache sizes according to &
                    &weights: '//&
                    trim(utils_real_to_str(cache_weights(1),'f12.4'))//', '//&
                    trim(utils_real_to_str(cache_weights(2),'f12.4'))//', '//&
                    trim(utils_real_to_str(cache_weights(3),'f12.4'))//'.'
            end if

            ! jd: Divide the memory to spare accordingly. First take care of
            !     SWOPs.
            mem_to_spare_mib = pub_hfx_memory_limit-minimum_required_memory_mib

            auto_mem_limit_for_swops = real(mem_to_spare_mib,kind=DP) * &
                 cache_weights(1)
         end if

         call comms_bcast(pub_root_proc_id,auto_mem_limit_for_swops)

         ! jd: Adjust SWOPs
         !     ... this is what we plan according to the memory limit.
         pub_cache_limit_for_swops = auto_mem_limit_for_swops

      end if

      ! jd: SWOP HT is adjusted down in case we can cache 100% of SWOPs,
      !     regardless of whether we have hfx_memory_limit or not.
      max_cached_swops_in_ppds = hash_table_calc_cache_size( &
           pub_cache_limit_for_swops, &
           swri%quality%max_sws_per_centre * mdl%cell%n_pts * SIZEOF_DOUBLE)

      !     ... if this is more than what is needed to hold everything,
      !     scale down. NB: n_swops can be, and usually is, different across
      !     procs. So the scalings might be different too, and, consequently,
      !     the memory pools assigned to expansions and products too. The
      !     banner shows what happens on root proc.
      if(max_cached_swops_in_ppds > n_swops) then
         max_cached_swops_in_ppds = n_swops
      end if
      ! jd: Reinitialise HT to the new size.
      call hash_table_free(swri%swops_in_ppds_ht)
      call hash_table_init(swri%swops_in_ppds_ht, 'SWOPS_IN_PPDS', 3, &
           max_cached_swops_in_ppds/4, max_cached_swops_in_ppds, &
           swri%hash_table_info_unit)
      ! jd: The table is not yet populated, so we calculate the size manually
      !     -- we know how many elements there will be.
      ! jd: Actual data
      mem_swops_mib = real(max_cached_swops_in_ppds,kind=DP) * &
           real(swri%quality%max_sws_per_centre * mdl%cell%n_pts, kind=DP) * &
           SIZEOF_DOUBLE/1024.0_DP/1024.0_DP
      ! jd: plus HT overheads
      mem_swops_mib = mem_swops_mib + real(node_size,kind=DP) * &
           real(swri%swops_in_ppds_ht%max_slots,kind=DP)/1024.0_DP/1024.0_DP
      mem_swops_long = mem_swops_mib * 1024.0_DP * 1024.0_DP
      pub_cache_limit_for_swops = mem_swops_mib

      ! jd: Expansions and products are adjusted if hfx_memory_limit set.
      if(pub_hfx_memory_limit /= -1) then

         if(pub_on_root) then ! mem_to_spare_mib is only initialised on root
            ! jd: Divide the rest across expansions and prods
            mem_to_spare_mib = mem_to_spare_mib - mem_swops_mib

            ! jd: Normalise weights again so that they add up to 1.
            cache_weights(1) = garbage_real
            weight_sum = sum(adjusted_cache_weights(2:3))
            if(weight_sum /= 0.0_DP) then
               call utils_assert(all(adjusted_cache_weights(2:3) >= 0.0_DP), &
                    myself//': cache weights(2:3) must all be non-negative.')
               call utils_assert(weight_sum > SAFE_DIV_EPS, &
                    myself//': cache weights(2:3) cannot all be zero.')
               cache_weights(2:3) = adjusted_cache_weights(2:3) / weight_sum

               auto_mem_limit_for_expansions = real(mem_to_spare_mib,kind=DP) * &
                    cache_weights(2)
               auto_mem_limit_for_prods = real(mem_to_spare_mib,kind=DP) * &
                    cache_weights(3)
            else
               ! jd: Corner case where both cache_limit_for_expansions and
               !     cache_limit_for_prods are 0. This means user doesn't want
               !     them.
               auto_mem_limit_for_expansions = 0
               auto_mem_limit_for_prods = 0
            end if
         end if

         call comms_bcast(pub_root_proc_id,auto_mem_limit_for_expansions)
         call comms_bcast(pub_root_proc_id,auto_mem_limit_for_prods)

         ! jd: Adjust expansions and products.
         pub_cache_limit_for_expansions = auto_mem_limit_for_expansions
         pub_cache_limit_for_prods = auto_mem_limit_for_prods

         ! jd: Reinitialise expansions HT to the new size.
         if(hfxstate%expansions_ht_in_use) then
            call hash_table_free(hfxstate%expansions_ht)
            max_cached_expansions = hash_table_calc_cache_size(&
                 pub_cache_limit_for_expansions, mdl%cell%n_pts * SIZEOF_DOUBLE)
            if(max_cached_expansions > 0) then
               call hash_table_init(hfxstate%expansions_ht, 'EXPANSIONS', 3, &
                    max_cached_expansions/5, max_cached_expansions, &
                    hfxstate%hash_table_info_unit) ! /5 to reduce overhd of procs
            end if
            hfxstate%expansions_ht_in_use = .not. (max_cached_expansions == 0)
         end if

         ! jd: Products HT does not need to be reinitialised (it has been
         !     reasonably auto-sized in init() and pub_cache_limit_for_prods
         !     will honoured when populating). We only need to make sure
         !     ppd_products_ht_in_use is consistent. If it was F, it stayed F
         !     and nothing to do. If it was T, but we asked for a zero weight
         !     for products, it now must become F so that lookups are not
         !     needlesly performed, and we need to free the HT.
         if(hfxstate%ppd_products_ht_in_use) then
            hfxstate%ppd_products_ht_in_use = &
                 .not. (pub_cache_limit_for_prods == 0)
            if(.not. hfxstate%ppd_products_ht_in_use) then
               call hash_table_free(hfxstate%ppd_products_ht)
            end if
         end if

      end if

      ! jd: Estimate memory use for expansions HT. The table is not yet
      !     populated. We assume it will be filled entirely.
      if(hfxstate%expansions_ht_in_use) then
         mem_expansions_mib = pub_cache_limit_for_expansions
      else
         mem_expansions_mib = 0.0_DP
      end if
      mem_expansions_long = mem_expansions_mib * 1024.0_DP * 1024.0_DP

      ! jd: Estimate memory use for PPD products HT. The table is not yet
      !     populated. We assume it will be filled entirely.
      if(hfxstate%ppd_products_ht_in_use) then
         mem_prods_mib = pub_cache_limit_for_prods
      else
         mem_prods_mib = 0.0_DP
      end if
      mem_prods_long = mem_prods_mib * 1024.0_DP * 1024.0_DP

      if(pub_hfx_output_detail >= BRIEF) then
         call utils_report_memory_estimate(&
              'HFx TEFCI engine: user-adjustable requirements',&
              [&
              'SWOP hash table                   ',&
              'Expansions hash table             ',&
              'AD NGWF products hash table       ' &
              ], &
              [&
              ! SWOP hash table
              mem_swops_long, &
              ! Expansions hash table
              mem_expansions_long, &
              ! PPD products hash table
              mem_prods_long &
              ], user_adjustable_memory)
         if(pub_on_root) then
            if(pub_hfx_memory_limit /= -1) then
               write(stdout,'(a,a,a)') 'HFx: - Peak memory use capped at ', &
                    trim(utils_real_to_str(real(minimum_required_memory + &
                    user_adjustable_memory,kind=DP)/1024.0_DP/1024.0_DP, &
                    'f16.1')), ' MB per MPI rank.'
            else
               write(stdout,'(a,a,a)') 'HFx: - Peak memory use estimated at ', &
                    trim(utils_real_to_str(real(minimum_required_memory + &
                    user_adjustable_memory,kind=DP)/1024.0_DP/1024.0_DP, &
                    'f16.1')), ' MB per MPI rank.'
            end if
         end if
      end if

    end subroutine internal_memory_use_estimate_and_adjust

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine hfx_populate_swop_cache

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_invalidate_ngwf_dependent_state(hfxstate)
    !==========================================================================!
    ! Invalidates all NGWF-dependent state in HFx.                             !
    ! This should be called every time NGWFs change -- in practice it's called !
    ! from hf_exchange_dkn_indep_stage().                                      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   hfxstate (inout): HFx state object. It's NGWF-dependent HTs will be    !
    !                     purged.                                              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2020.                                 !
    !==========================================================================!

    use hash_table, only: hash_table_purge

    implicit none

    ! jd: Arguments
    type(HFX_STATE), intent(inout), target :: hfxstate

    ! -------------------------------------------------------------------------

    ! jd: Invalidate old NGWF expansions, since NGWFs changed. This can't be
    !     done from ngwf_rep_mod, as it would lead to a circular dependency
    !     between hf_exchange_mod and ngwf_rep_mod. [footnote 1].
    ! jd: Note that we can get away with having only *one* hash table despite
    !     multiple NGWF sets, because we *never* make subsequent calls to
    !     hf_exchange_calculate(), with different NGWF-set-pairs, without an
    !     intervening purge of expansions_ht, which happens in expand_ngwf_pairs().
    !     That is, we always expand first, then make a call (in conduction) or
    !     a number of calls (for LNV iterations, in valence), without the NGWFs
    !     changing and *with the same NGWF-set-pairs* in the call. When we
    !     switch to a new NGWF-set-pair, we always re-expand first.
    if(hfxstate%expansions_ht_in_use) then
       call hash_table_purge(hfxstate%expansions_ht)
    end if

    ! jd: Aa-Dd NGWF products in PPDs and Dd NGWFs in PPDs must also be purged
    !     since NGWFs changed.
    if(hfxstate%ppd_products_ht_in_use) then
       call hash_table_purge(hfxstate%ppd_products_ht)
    end if
    call hash_table_purge(hfxstate%ppd_ngwfs_dd_ht)

  end subroutine hfx_invalidate_ngwf_dependent_state

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_get_dkn_atomblock(dkn_atomblock, packed_dkn_blocks_ht, &
       global_row, global_col, row_max_on_atom, col_max_on_atom)
    !==========================================================================!
    ! Retrieves a DKN atomblock from the DKN hash table.                       !
    ! It is assumed the atomblock is present, we abort if it is not found.     !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   dkn_atomblock (out): Result is returned here. The dimensions are       !
    !                        (row_max_on_atom, col_max_on_atom, pub_num_spins).!
    !   packed_dkn_blocks_ht (in): Hash table with DKN blocks from which we    !
    !                              extract the atomblock.                      !
    !   global_row (in): Global atom index for the row (typically A or C).     !
    !   global_col (in): Global atom index for the column (typically B or D).  !
    !   row_max_on_atom (in): Maximum number of NGWFs per atom in row basis.   !
    !   col_max_on_atom (in): Maximum number of NGWFs per atom in col basis.   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in April 2020.                                 !
    !==========================================================================!

    use constants, only: garbage_real
    use hash_table, only: HT_HASH_TABLE, hash_table_lookup_ptr_nocount
    use rundat, only: pub_num_spins, pub_debug
    use utils, only: utils_int_to_str, utils_abort, utils_sanity_check

    implicit none

    ! jd: Arguments
    integer, intent(in) :: row_max_on_atom ! } jd: Defines the size of the
    integer, intent(in) :: col_max_on_atom ! }     atomblock
    real(kind=DP), intent(out) :: dkn_atomblock(row_max_on_atom, &
         col_max_on_atom, pub_num_spins)
    type(HT_HASH_TABLE), intent(in) :: packed_dkn_blocks_ht
    integer, intent(in) :: global_row ! jd: Global atom index
    integer, intent(in) :: global_col ! jd: Global atom index

    ! jd: Local variables
    real(kind=DP), pointer :: packed_dkn_block_ptr(:) ! points to HT's internals
    integer :: packed_dkn_block_size_in_ht
    integer :: packed_dkn_block_size
    character(len=*), parameter :: myself = 'hfx_get_dkn_atomblock'

    ! -------------------------------------------------------------------------

    if(pub_debug) then
       dkn_atomblock(:,:,:) = garbage_real ! helps uncover padding bugs
    end if

    packed_dkn_block_size = pub_num_spins * row_max_on_atom * col_max_on_atom

    call hash_table_lookup_ptr_nocount(packed_dkn_block_ptr, &
         packed_dkn_block_size_in_ht, packed_dkn_blocks_ht, &
         global_row, global_col)
    if(packed_dkn_block_size_in_ht /= packed_dkn_block_size) then
       call utils_abort(myself//': DKN packed block size mismatch. Got: '//&
            trim(utils_int_to_str(packed_dkn_block_size_in_ht))//&
            ', expected: '//trim(utils_int_to_str(packed_dkn_block_size))//'.',&
            global_row, global_col)
    end if
    dkn_atomblock(:,:,:) = reshape(packed_dkn_block_ptr, &
         [row_max_on_atom, col_max_on_atom, pub_num_spins])

    if(pub_debug) then
       call utils_sanity_check(dkn_atomblock, 'dkn_atomblock', excessive=1D10)
    end if

  end subroutine hfx_get_dkn_atomblock

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hfx_dkn_blocks_ht_size(max_cached_dkn_blocks_cd, &
       max_cached_dkn_blocks_ab, ngwf_indices, natoms)
    !==========================================================================!
    ! Calculates reasonable values for the sizes of DKN hash tables.           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   max_cached_dkn_blocks_cd (out): Calculated size for K^CD blocks is     !
    !                                   returned here.                         !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in April 2020.                                 !
    !==========================================================================!

    use constants, only: SIZEOF_DOUBLE
    use hash_table, only: hash_table_calc_cache_size
    use rundat, only: pub_cache_limit_for_dknblks, pub_num_spins

    implicit none

    ! jd: Arguments
    integer, intent(out)             :: max_cached_dkn_blocks_cd
    integer, intent(out)             :: max_cached_dkn_blocks_ab
    type(HFX_NGWF_INDEX), intent(in) :: ngwf_indices(A_NGWF:D_NGWF)
    integer, intent(in)              :: natoms

    ! -------------------------------------------------------------------------

    if(pub_cache_limit_for_dknblks == -1) then
       ! jd: No user-specified value, use default. The default is '400/atom, but
       !     no more than 1M slots'.
       max_cached_dkn_blocks_cd = min(maxslots_dkn_blocks_per_atom * natoms, &
            max_cached_dkn_blocks_default)
       max_cached_dkn_blocks_ab = max_cached_dkn_blocks_cd
    else
       ! jd: User overrode default, use their value
       max_cached_dkn_blocks_cd = hash_table_calc_cache_size( &
            pub_cache_limit_for_dknblks, &
            ngwf_indices(C_NGWF)%basis%max_on_atom * SIZEOF_DOUBLE * &
            ngwf_indices(D_NGWF)%basis%max_on_atom * pub_num_spins)
       max_cached_dkn_blocks_ab = hash_table_calc_cache_size( &
            pub_cache_limit_for_dknblks, &
            ngwf_indices(A_NGWF)%basis%max_on_atom * SIZEOF_DOUBLE * &
            ngwf_indices(B_NGWF)%basis%max_on_atom * pub_num_spins * 3)
       !                                      (K^AB, tcK^AB, tcK^BA) ^
    end if

  end subroutine hfx_dkn_blocks_ht_size

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module hf_exchange
