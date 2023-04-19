! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!-----------------------------------------------------------------------!
!               P O L A R I S A B L E   E M B E D D I N G               !
!-----------------------------------------------------------------------!
!                                                                       !
! This module implements embedding classical charges, dipoles and       !
! quadrupoles (termed 'embedding multipoles' from now on) that can      !
! change over the course of an SCF cycle.                               !
!                                                                       !
! The values of the embedding multipole components would normally vary  !
! in response to changes in the electronic density during the SCF cycle,!
! while their positions would be fixed during the SCF cycle but would   !
! change as we move QM atoms (in geometry optimisation, MD, etc).       !
!                                                                       !
! This functionality has been described in                              !
! [1] J. Dziedzic, Y. Mao, Y. Shao, J. Ponder, T. Head-Gordon,          !
!     M. Head-Gordon, C.-K. Skylaris, J. Chem. Phys. 145, 124106 (2016).!
!                                                                       !
! Currently supported are:                                              !
! - Permanent classical monopoles, dipoles and quadrupoles, acting on   !
!   the full ONETEP density (and cores) through Coulombic tensors with  !
!   crude singularity masking or (default) Thole-like smearing,         !
! - Inducible (polarisable) dipoles, acting on the ONETEP DMA-decomposed!
!   density through Coulombic or Thole-damped tensors. The approach     !
!   where the energy expression does not involve integration of the QM  !
!   density with an external potential, but rather involves QM multi-   !
!   poles obtained from DMA is termed the QM* (qmstar) representation.  !
!   This representation, although vastly more complex, allows using     !
!   Thole damping, which operates on pairs of atoms and cannot be used  !
!   for the entire density.                                             !
! - Since Feb 2017 both LNV and NGWF gradients are supported for the    !
!   QM* representation. Gradients for the full density representation   !
!   are trivial (only involve an addition of an external potential to   !
!   lhxc). In the QM* representation, the LNV gradient is calculated in !
!   the matrix representation (in the NGWF basis) -- a 'P' or 'polemb'  !
!   matrix is added to the Hamiltonian. The NGWF gradient is composed of!
!   two terms. The first term is handled automatically by ONETEP once   !
!   the P matrix is included in the Hamiltonian. The second term is more!
!   tricky and is evaluated in polarisable_embedding_ngwf_gradient().   !
! - Currently the NGWF gradient does not support DMA charge scaling,    !
!   an assert checks against trying to do that.                         !
! - Permanent classical monopoles, dipoles and quadrupoles, acting on   !
!   the ONETEP DMA-decomposed density through Coulombic or Thole-damped !
!   tensors. Same notes on gradients as for induced dipoles. This is    !
!   coded up correctly, but will likely lead to excessive charge sucking!
!   by the MM atoms, reflecting a catastrophic breakdown of the multipo-!
!   le expansion within the NGWF sphere.                                !
! - Both mutual and direct induction work. With mutual induction the    !
!   LNV gradients depend on tight convergence of SCF dipoles. With      !
!   direct induction they are exact.                                    !
! - Since February 2017 the NGWF gradient can be computed on the double !
!   grid (pol_emb_dbl_grid T). By default it is computed on the coarse  !
!   grid, like in HFx. Technically, to prevent aliasing, the double grid!
!   version should be used, but it is unbearably inefficient and has not!
!   been optimised much, as I plan this to remain an experimental faci- !
!   lity at most. To ensure consistency, pol_emb_dbl_grid T should only !
!   be used together with swx_dbl_grid T, which ensures NGWF-swop       !
!   overlaps are calculated on the double grid too. I strongly recommend!
!   using the defaults (pol_emb_dbl_grid F, swx_dbl_grid F).            !
!                                                                       !
! Written by Jacek Dziedzic in 2015.05-2017.06.                         !
!                                                                       !
!-----------------------------------------------------------------------!
#define FFTBOX_DIMS mdl%fftbox%total_ld1,mdl%fftbox%total_ld2,mdl%fftbox%total_pt3
#define FFTBOX_DBL_DIMS mdl%fftbox%total_ld1_dbl,mdl%fftbox%total_ld2_dbl,mdl%fftbox%total_pt3_dbl
#define PPDS_DIMS ngwf_basis%max_n_ppds_sphere*cell%n_pts
#define MDL_PPDS_DIMS ngwf_basis%max_n_ppds_sphere*mdl%cell%n_pts

module polarisable_embedding

  use constants, only: CRLF, DP, VERBOSE, stdout
  use kernel, only: DKERN
  use multipole_ops, only: SPHERICAL_MULTIPOLE_SET
  use ngwf_representation, only: NGWF_REP

  implicit none

  public :: polarisable_embedding_init
  public :: polarisable_embedding_init_mm_rep_params
  public :: polarisable_embedding_expand_ngwf_pairs
  public :: polarisable_embedding_calculate
  public :: polarisable_embedding_ngwf_gradient
  public :: polarisable_embedding_exit

  private

  ! jd: The maximum number of multipole sets that can be used at any time.
  !     This is used to dimension arrays. Don't be afraid to increase this
  !     if need be.
  integer, parameter :: MAX_N_MULTIPOLE_SETS = 2

  ! jd: An array for caching potential of sets of fixed multipoles
  real(kind=DP), allocatable, save :: fixed_sets_polemb_locpot(:,:,:,:)
  real(kind=DP), allocatable, save :: stored_local_rep_mm_pol_emb_pot(:,:,:)
  ! jd: ... and their core interaction energies.
  real(kind=DP), save :: fixed_core_energy_terms(MAX_N_MULTIPOLE_SETS)
  ! jd: Tracks if the above have been filled
  logical, save :: fixed_sets_calculated = .false.
  logical, save :: repulsive_mm_pot_calculated = .false.

  ! jd: MM repulsive potential params for all MM species.
  !     1st: param# (1 or 2), 2nd: MM species index.
  real(kind=DP), allocatable, save :: mm_rep_params(:,:)   ! } only allocated if
  ! jd: MM species names for all MM species.               ! } mm_rep_params
  character(len=4), allocatable, save :: mm_all_species(:) ! } block present

  type VACUUM_STATE
     type(NGWF_REP)                :: rep
     type(DKERN)                   :: kdenskern
     type(SPHERICAL_MULTIPOLE_SET) :: vacuum_multipoles(3)
     type(SPHERICAL_MULTIPOLE_SET) :: vacuum_aux_multipoles(3)
     logical                       :: initialised = .false.
  end type VACUUM_STATE

  type(VACUUM_STATE), save, public :: polarisable_embedding_vacuum_state

  ! jd: These reflect the fact that we're working with up-to-quadrupolar
  !     MM force fields.
  integer, parameter :: n_qterm_sph_components = 9    ! 1(l=0) + 3(l=1) + 5(l=2)
  integer, parameter :: n_qterm_cart_components = 13  ! 1(l=0) + 3(l=1) + 9(l=2)

  integer, save                :: pol_emb_scf_iter = 0

  integer, save                :: hash_table_info_unit
  character(len=256), save     :: hash_table_info_filename
  character(len=32), parameter :: hash_table_info_rootname = 'hash_table_polemb'

contains

  ! ---------------------------------------------------------------------------
  ! ----- P U B L I C   S U B R O U T I N E S   A N D   F U N C T I O N S -----
  ! ---------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine polarisable_embedding_init(ngwf_basis, mdl, n_occ, overlap)!, denskern)
    !==========================================================================!
    ! Sets up the module, currently this amounts to:                           !
    ! - setting fixed_sets_calculated to .false.,                              !
    ! - preparing the Gamma lookup,                                            !
    ! - allocating caches for locpots due to fixed sets.                       !
    ! - allocating caches for MM reppots due to fixed sets.                    !
    ! - issuing a relevant citation to the .bib file                           !
    ! - initialising the vacuum state, iff pol_emb_vacuum_qmstar T.            !
    !--------------------------------------------------------------------------!
    !  Arguments:                                                              !
    !    overlap and denskern are only used for their sparsity pattern.        !
    !    n_occ only needed to set up *vacuum* state's n_occ.                   !
    ! JCAP: Is denskern actually used? It doesn't seem to be.                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2015.07.                                    !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,         !
    ! September 2018                                                           !
    !==========================================================================!

    use bibliography, only: bibliography_cite
    use comms, only: pub_my_proc_id
    use constants, only: REP_N_SWEXES
    use function_basis, only: FUNC_BASIS
    !use kernel, only: DKERN
    use model_type, only: MODEL
    use multipole_ops, only: multipole_prepare_gamma_lookup
    use rundat, only: pub_aug, pub_dma_calculate, pub_pol_emb_vacuum_qmstar, &
         pub_use_swx
    use sparse, only: SPAM3
    use sparse_embed, only: SPAM3_EMBED
    use sw_expansion_type, only: SW_EX
    use utils, only: utils_assert, utils_alloc_check, utils_abort, utils_unit

    implicit none

    ! jd: Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    type(MODEL), intent(in)      :: mdl
    integer, intent(in)          :: n_occ(:,:)
    type(SPAM3_EMBED), intent(in)      :: overlap
!    type(DKERN), intent(in)      :: denskern

    ! jd: Local variables
    integer :: n_qterm_components
    integer :: ierr
    character(len=*), parameter :: myself = 'polarisable_embedding_init'

    ! -------------------------------------------------------------------------

    call utils_assert(.not. pub_aug, &
         'Polarisable embedding cannot be used with charge augmentation')
    call utils_assert(pub_dma_calculate, 'Polarisable embedding requires DMA to&
         & be enabled (dma_calculate T + further set-up)')
    call utils_assert(pub_use_swx, 'Polarisable embedding requires SWx.&
         & This indicates an internal error in ONETEP''s logic')

    call utils_assert(.not. allocated(fixed_sets_polemb_locpot), &
         myself//': fixed_sets_polemb_locpot already allocated.')

    call bibliography_cite('POLEMB')

    ! jd: Allocate locpots for fixed multipole sets
    allocate(fixed_sets_polemb_locpot(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, MAX_N_MULTIPOLE_SETS), stat=ierr)
    call utils_alloc_check(myself,'fixed_sets_polemb_locpot',ierr)
    allocate(stored_local_rep_mm_pol_emb_pot(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12), stat=ierr)
    call utils_alloc_check(myself,'stored_local_rep_mm_pol_emb_pot',ierr)
    fixed_sets_calculated = .false.

    call multipole_prepare_gamma_lookup()

    ! jd: Ensure vacuum state is not initialised to boot
    polarisable_embedding_vacuum_state%initialised = .false.

    ! jd: ... and initialise it, if needed
    if(pub_pol_emb_vacuum_qmstar) then
       call polarisable_embedding_init_vacuum_state(&
            polarisable_embedding_vacuum_state, ngwf_basis, mdl, n_occ, overlap)
!            denskern)
    end if

    ! jd: Open a file for debugging hash-tables
    hash_table_info_unit = utils_unit()
    write(hash_table_info_filename,'(a,a,a,a,i0,a)') &
         trim(hash_table_info_rootname), '_', '', &
         '_proc_', pub_my_proc_id, '.log'
    open(hash_table_info_unit, file=hash_table_info_filename, err=10)

    return

10  call utils_abort('Could not create file: '//trim(hash_table_info_filename))

  end subroutine polarisable_embedding_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine polarisable_embedding_init_mm_rep_params(n_mm_species, cbuf, dbuf)
    !==========================================================================!
    ! @docme
    !--------------------------------------------------------------------------!
    !  Arguments:                                                              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2018.03.                                    !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout, CRLF
    use utils, only: utils_assert, utils_alloc_check

    implicit none

    ! jd: Arguments
    integer, intent(in)           :: n_mm_species
    character(len=10), intent(in) :: cbuf(n_mm_species)
    real(kind=DP), intent(in)     :: dbuf(2,n_mm_species)

    ! jd: Local variables
    integer :: ispecies
    integer :: ierr
    character(len=*), parameter :: myself = 'polarisable_embedding_init'

    ! -------------------------------------------------------------------------

    call utils_assert(n_mm_species > 0, myself//': Invalid n_mm_species')

    allocate(mm_all_species(n_mm_species), stat=ierr)
    call utils_alloc_check(myself,'mm_all_species',ierr)
    allocate(mm_rep_params(2,n_mm_species), stat=ierr)
    call utils_alloc_check(myself,'mm_rep_params',ierr)

    do ispecies = 1, n_mm_species
       mm_all_species(ispecies) = trim(cbuf(ispecies)(1:4))
       mm_rep_params(:,ispecies) = dbuf(:,ispecies)
    end do

    if(pub_on_root) then
       write(stdout,'(/a)') &
            '----------------------------------------------------------'&
            //CRLF//&
            '|           MM repulsive potential parameters:           |'&
            //CRLF//&
            '+--------------------------------------------------------+'&
            //CRLF//&
            '|    Species | Amplitude (Ha/e) |           Zeta (a0^-1) |'
       do ispecies = 1, n_mm_species
          write(stdout,'(a2,a10,a2,f17.3,a6,f19.6,a2)') '| ', &
               adjustr(mm_all_species(ispecies)), ' |', &
               mm_rep_params(1, ispecies), ' |    ',&
               mm_rep_params(2, ispecies),' |'
       end do
       write(stdout,'(a)') '---------------------------------------------------&
            &-------'
    end if


  end subroutine polarisable_embedding_init_mm_rep_params


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine polarisable_embedding_exit
    !==========================================================================!
    ! Cleans up after the module, deallocating the array holding the cached    !
    ! local potential of fixed multipole sets, and destroying the Gamma lookup.!
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2015.07.                                    !
    !==========================================================================!

    use datatypes, only: data_functions_dealloc
    use dma, only: dma_free_multipoles
    use kernel, only: kernel_destroy
    use multipole_ops, only: multipole_destroy_gamma_lookup
    use ngwf_representation, only: ngwf_rep_destroy
    use utils, only: utils_assert, utils_dealloc_check, utils_abort

    implicit none

    ! jd: Local variables
    integer :: ierr
    character(len=*), parameter :: myself = 'polarisable_embedding_exit'

    ! -------------------------------------------------------------------------

    call utils_assert(allocated(fixed_sets_polemb_locpot), &
         myself//': fixed_sets_polemb_locpot already allocated.')

    deallocate(stored_local_rep_mm_pol_emb_pot, stat=ierr)
    call utils_dealloc_check(myself,'stored_local_rep_mm_pol_emb_pot',ierr)
    deallocate(fixed_sets_polemb_locpot, stat=ierr)
    call utils_dealloc_check(myself,'fixed_sets_polemb_locpot',ierr)

    fixed_sets_calculated = .false.
    call multipole_destroy_gamma_lookup

    if(polarisable_embedding_vacuum_state%initialised) then
       call data_functions_dealloc(polarisable_embedding_vacuum_state%rep%ngwfs_on_grid(1))
       call ngwf_rep_destroy(polarisable_embedding_vacuum_state%rep)
       call kernel_destroy(polarisable_embedding_vacuum_state%kdenskern)
       call dma_free_multipoles(polarisable_embedding_vacuum_state%vacuum_multipoles)
       call dma_free_multipoles(polarisable_embedding_vacuum_state%vacuum_aux_multipoles)
       polarisable_embedding_vacuum_state%initialised = .false.
    end if

    ! jd: MM rep params and MM species were only allocated if an mm_rep_params
    !     block was present.
    if(allocated(mm_rep_params)) then
       deallocate(mm_rep_params, stat=ierr)
       call utils_dealloc_check(myself,'mm_rep_params',ierr)
    end if
    if(allocated(mm_all_species)) then
       deallocate(mm_all_species, stat=ierr)
       call utils_dealloc_check(myself,'mm_all_species',ierr)
    end if

    close(hash_table_info_unit, err=10)

    return

10  call utils_abort('Could not close file: '//hash_table_info_filename)

  end subroutine polarisable_embedding_exit

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine polarisable_embedding_calculate(polemb, vacuum_polemb, & ! in/out
       vacuum_aux_polemb, &                                    ! in/out
       polemb_energy_mm_valence, &                             ! out
       polemb_energy_mm_vdw, &                                 ! out
       polemb_energy_mm_es, &                                  ! out
       polemb_energy_qmstar_elec_mm_perm, &                    ! out
       polemb_energy_qmstar_elec_mm_ind, &                     ! out
       polemb_energy_qmfull_elec_mm_perm, &                    ! out
       polemb_energy_qmfull_elec_mm_ind, &                     ! out
       polemb_energy_qm_core_mm_perm, &                        ! out
       polemb_energy_qm_core_mm_ind, &                         ! out
       polemb_energy_qmmm_vdw, &                               ! out
       polemb_energy_other, &                                  ! out
       lhxc_fine, lhxc_energy, &                               ! inout
       density_on_fine, mdl, rep, ngwf_basis, denskern, &      ! in
       energy_term_names, energy_term_values, n_energy_terms)  ! out, opt

    !==========================================================================!
    ! Driver routine for polarisable embedding.                                !
    ! 1) Performs DMA analysis of current electronic density, outputting the   !
    !    multipoles to a GDMA-like file.                                       !
    ! 2) Pings a FIFO lock to let the MM side know the multipoles are ready.   !
    ! 3) Waits on a FIFO lock until the MM side acks their readiness.          !
    ! 4) Parses a file generated by the MM side, containing positions and      !
    !    values of permanent and SCF-induced multipoles, and also any external !
    !    energy terms that are to be included.                                 !
    ! 5) Unless told otherwise, *adds* the smeared-Coulombic potential         !
    !    of MM multipoles in non-QM* sets to lhxc_fine, and the energy of their!
    !    interaction with electrons to lhxc_energy. The energy of their inter- !
    !    action with cores is returned in qm_core_mm_{perm,ind}.               !
    ! 6) Constructs the polarisation ('P') matrices for MM multipoles in QM*   !
    !    sets and stores them in polemb. Index 1 stores contributions from     !
    !    MM induced dipoles, index 2 stores contributions from MM permanent    !
    !    multipoles. If either set is non-QM*, the corresponding matrix is 0.  !
    ! 7) Returns a number of energy terms associated with the embedding.       !
    ! 8) Returns MM valence terms (bond, angle, etc) in mm_valence.            !
    ! 9) Returns MM vdW energies and QM/MM vdW energies.                       !
    ! A) Returns any other MM-side terms in energy_other.                      !
    ! B) If the last three optional arguments are specified, returns a compre- !
    !    hensive breakdown of energy terms in two arrays (names, values), each !
    !    n_energy_terms-long. This is useful in hamiltonian_energy_components. !
    !    It is then, and only then, the responsibility of the caller to deallo-!
    !    cate these two arrays. Energy terms with '#' in the name should be    !
    !    ignored.                                                              !
    !--------------------------------------------------------------------------!
    ! Note:                                                                    !
    ! - The matrix-pairs polemb(:), vacuum_polemb(:) and vacuum_aux_polemb(:)  !
    !   are only passed as in/out for their sparsity pattern. They will be     !
    !   zeroed first, they are *not* added to.                                 !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    ! - Returns MM *and* QM/MM vdW terms in mm_vdw.                            !
    ! - Careful with double-counting: polemb_energy_qmstar_elec_mm_ind and     !
    !   polemb_energy_qmstar_elec_mm_perm will be *added* to lhxc_fine.        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2015.11-2015.12.                            !
    ! Reworked by Jacek Dziedzic in 2016.09 for possibility of QM* perm mpoles.!
    !==========================================================================!

    use comms, only: comms_barrier, pub_on_root
    use constants, only: REP_SWEX_POL_EMB_DMA_1, REP_SWEX_POL_EMB_DMA_2, &
         REP_SWEX_POL_EMB_AUX_DMA_1, REP_SWEX_POL_EMB_AUX_DMA_2, SIZEOF_DOUBLE
    use dma, only: dma_calculate, dma_free_multipoles, dma_write_multipoles, &
         dma_expand_pair_densities
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE, hash_table_init, hash_table_free, &
         hash_table_calc_cache_size
    use integrals, only: integrals_locpot
    use model_type, only: MODEL
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET, multipole_set_merge, &
         multipole_set_info, multipole_set_scale_monopoles, &
         multipole_set_scale_quadrupoles
    use rundat, only: PUB_1K, pub_pol_emb_pot, pub_num_spins, pub_dma_use_ri, &
         pub_devel_code, pub_inner_loop_iteration, pub_pol_emb_vacuum_qmstar, &
         pub_dma_bessel_averaging
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_scale
    use sparse_array, only: SPAM3_ARRAY
    use sparse_embed, only: sparse_embed_extract_from_array, &
         sparse_embed_destroy_extracted_array
    use sw_expansion_type, only: swx_coeffs_n_slots
    use sw_resolution_of_identity, only: swri_library, swri_get_handle_to
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_lock_ping, utils_lock_block, &
         utils_alloc_check, utils_dealloc_check, utils_int_to_str, &
         utils_devel_code
    use visual, only: visual_scalarfield, visual_scalarfield_ngwf_basis

    implicit none

    ! jd: Arguments
    type(SPAM3), intent(inout)       :: polemb(:) ! 1: induced, 2: permanent
    type(SPAM3), intent(inout)       :: vacuum_polemb(:) ! 1: induced, 2: permanent
    type(SPAM3), intent(inout)       :: vacuum_aux_polemb(:) ! 1: induced, 2: permanent
    real(kind=DP), intent(out)       :: polemb_energy_mm_valence
    real(kind=DP), intent(out)       :: polemb_energy_mm_vdw
    real(kind=DP), intent(out)       :: polemb_energy_mm_es
    real(kind=DP), intent(out)       :: polemb_energy_qmstar_elec_mm_perm
    real(kind=DP), intent(out)       :: polemb_energy_qmstar_elec_mm_ind
    real(kind=DP), intent(out)       :: polemb_energy_qmfull_elec_mm_perm
    real(kind=DP), intent(out)       :: polemb_energy_qmfull_elec_mm_ind
    real(kind=DP), intent(out)       :: polemb_energy_qm_core_mm_perm
    real(kind=DP), intent(out)       :: polemb_energy_qm_core_mm_ind
    real(kind=DP), intent(out)       :: polemb_energy_qmmm_vdw
    real(kind=DP), intent(out)       :: polemb_energy_other
    real(kind=DP), intent(inout)     :: lhxc_fine(:,:,:,:)
    real(kind=DP), intent(inout)     :: lhxc_energy
    real(kind=DP), intent(in)        :: density_on_fine(:,:,:,:)
    type(MODEL), intent(in)          :: mdl
    type(NGWF_REP), intent(in)       :: rep
    type(FUNC_BASIS), intent(in)     :: ngwf_basis(1)
    type(SPAM3_ARRAY), intent(in)    :: denskern
    character(len=27), intent(out), optional, allocatable :: energy_term_names(:)
    real(kind=DP), intent(out), optional, allocatable   :: energy_term_values(:)
    integer, intent(out), optional                      :: n_energy_terms

    ! jd: Local variables
    character(len=27), allocatable   :: loc_energy_term_names(:)
    real(kind=DP), allocatable       :: loc_energy_term_values(:)
    integer                          :: loc_n_energy_terms
    real(kind=DP), allocatable       :: loc_pol_emb_pot_fine(:,:,:)
    real(kind=DP), allocatable       :: loc_rep_mm_pol_emb_pot_fine(:,:,:)
    integer                          :: i_term
    ! jd: DMA multipoles in real spherical representation for the total density
    !     at this SCF step. There are three sets -- one for max_q, one for
    !     max_q-1 (Bessel averaging), the last one for the averaging
    type(SPHERICAL_MULTIPOLE_SET)    :: scf_dma_polarisation_multipoles(3)
    type(SPHERICAL_MULTIPOLE_SET)    :: scf_dma_total_multipoles(3)
    integer                          :: is
    integer                          :: ierr
    integer                          :: dma_swri_h
    type(SPAM3)                      :: tmp1, tmp2, tmp3, tmp_spam3
    type(HT_HASH_TABLE), target      :: coeffs_ht
    integer                          :: i_set
    integer                          :: natoms
    integer, save                    :: pot_output_iter = 1
    logical, save                    :: pot_output_done = .false.
    character(len=*), parameter :: myself = 'polarisable_embedding_calculate'

    ! jcap: temporary arrays
    type(SPAM3), allocatable  :: kern_array(:)

    ! -------------------------------------------------------------------------

    call utils_assert(pub_pol_emb_pot, myself//&
         ': Cannot call this subroutine without "pol_emb_pot_filename"')
    call utils_assert(mdl%nat_classical == 0, &
         'Polarisable embedding cannot be used simultaneously with &
         &classical atoms.')

    dma_swri_h = swri_get_handle_to(pub_dma_use_ri)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! ---------------
    ! --- SCF DMA ---
    ! ---------------
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    natoms = size(mdl%elements)

    call hash_table_init(coeffs_ht, 'POLEMB_COEFFS', 5, &
         swx_coeffs_n_slots(natoms), swx_coeffs_n_slots(natoms), &
         rep%swexes(REP_SWEX_POL_EMB_DMA_1)%hash_table_info_unit)

    if(pub_pol_emb_vacuum_qmstar) then
       allocate(kern_array(pub_num_spins),stat=ierr)
       call utils_alloc_check(myself,'kern_array',ierr)

       ! jd: Expand the polarisation density
       ! jcap: the 1 here corresponds to a region - only set up for 1
       ! region at the moment

       call sparse_embed_extract_from_array(kern_array,&
            polarisable_embedding_vacuum_state%kdenskern%kern%m(:,PUB_1K))
       call dma_expand_pair_densities(coeffs_ht, &
            rep, 1, denskern%m(:,PUB_1K), ngwf_basis(1), mdl, REP_SWEX_POL_EMB_DMA_1,&
            vacuum_rep = polarisable_embedding_vacuum_state%rep, &
            vacuum_denskern = kern_array)
       call sparse_embed_destroy_extracted_array(kern_array)
       call dma_calculate(scf_dma_polarisation_multipoles, 'scf_dma_polarisation', &
            rep%swexes(REP_SWEX_POL_EMB_DMA_1:REP_SWEX_POL_EMB_DMA_2), &
            rep%n_occ, rep%overlap, ngwf_basis, mdl, REP_SWEX_POL_EMB_DMA_1, &
            'D', coeffs_ht = coeffs_ht)

       deallocate(kern_array,stat=ierr)
       call utils_dealloc_check(myself,'kern_array',ierr)

!       write(*,*) '@@@ HACK FOR MASKING CHARGES AND QPOLES in polarisation multipoles @@@'
!       do i_set = 1, merge(3,1,pub_dma_bessel_averaging)
!          call multipole_set_scale_monopoles(scf_dma_polarisation_multipoles(i_set), 0D0)
!          call multipole_set_scale_quadrupoles(scf_dma_polarisation_multipoles(i_set), 0D0)
!       end do

!       call multipole_set_info(scf_dma_polarisation_multipoles(1), stdout, .false.)
       ! jd: Combine multipoles from polarisation density with vacuum multipoles
       !     obtaining total multipoles that drive MM
       do i_set = 1, merge(3,1,pub_dma_bessel_averaging)
          ! jd: Merge delta with 0
          call multipole_set_merge(scf_dma_total_multipoles(i_set), &
               'tmp', scf_dma_polarisation_multipoles(i_set), &
               polarisable_embedding_vacuum_state%vacuum_multipoles(i_set), &
               1D0, 1D0)
       end do
!       call multipole_set_info(polarisable_embedding_vacuum_state%vacuum_multipoles(1), stdout, .false.)
!       call multipole_set_info(scf_dma_total_multipoles(1), stdout, .false.)
       call dma_write_multipoles('scf_dma_total', scf_dma_total_multipoles, &
            mdl%elements, mdl%par, dma_swri_h, gather_first = .true.)
    else
       ! jd: Calculate and output total multipoles [these drive the MM subsystem]
       ! jcap: the 1 here corresponds to a region - only set up for 1
       ! region at the moment
       call dma_expand_pair_densities(coeffs_ht, &
            rep, 1, denskern%m(:,PUB_1K), ngwf_basis(1), mdl, REP_SWEX_POL_EMB_DMA_1)
       call dma_calculate(scf_dma_total_multipoles, 'scf_dma_total', &
            rep%swexes(REP_SWEX_POL_EMB_DMA_1:REP_SWEX_POL_EMB_DMA_2), &
            rep%n_occ, rep%overlap, ngwf_basis, mdl, REP_SWEX_POL_EMB_DMA_1, &
            'D', coeffs_ht = coeffs_ht)
!       call multipole_set_info(scf_dma_total_multipoles(1), stdout, .false.)
    end if

    call hash_table_free(coeffs_ht)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    call timer_clock('work imbalance (pol emb)',1)
    call comms_barrier
    call timer_clock('work imbalance (pol emb)',2)
    call timer_clock('ONETEP_waiting_for_MM',1)
    ! jd: Let the outside world know the multipoles are ready
    if(pub_on_root) call utils_lock_ping('$QM2MM.lock')
    ! jd: Wait until the outside world acknowledged
    if(pub_on_root) call utils_lock_block('$MM2QM.lock')
    call comms_barrier
    call timer_clock('ONETEP_waiting_for_MM',2)
    pol_emb_scf_iter = pol_emb_scf_iter + 1

    ! -------------------------------------------------------
    ! --- POL EMB: Calculation of matrix and energy terms ---
    ! -------------------------------------------------------

    polemb_energy_mm_valence = 0D0
    polemb_energy_mm_vdw = 0D0
    polemb_energy_mm_es = 0D0
    polemb_energy_qmstar_elec_mm_perm = 0D0
    polemb_energy_qmstar_elec_mm_ind = 0D0
    polemb_energy_qmfull_elec_mm_perm = 0D0
    polemb_energy_qmfull_elec_mm_ind = 0D0
    polemb_energy_qm_core_mm_perm = 0D0
    polemb_energy_qm_core_mm_ind = 0D0
    polemb_energy_qmmm_vdw = 0D0
    polemb_energy_other = 0D0

    allocate(loc_pol_emb_pot_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12), stat=ierr)
    call utils_alloc_check(myself,'loc_pol_emb_pot_fine', ierr)
    allocate(loc_rep_mm_pol_emb_pot_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12), stat=ierr)
    call utils_alloc_check(myself,'loc_rep_mm_pol_emb_pot_fine', ierr)

    ! jd: Calculate embedding potential and energy terms.
    !     The energy terms are:
    !     - interaction of external potential with the electronic density,
    !     - interaction of external potential with the cores,
    !     - energy terms calculated outside of ONETEP, if any (these are
    !       communicated through a file).
    !     The above terms are stored in pol_emb_energy_term_values(:) with
    !     a corresponding array for names. Whether a particular term is
    !     included, depends on the keywords. The variable pol_emb_n_energy
    !     returns the number of these energy terms, which can be zero.
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------
    if(pub_pol_emb_vacuum_qmstar) then
       call pol_emb_potential_and_energy(polemb, vacuum_polemb, &      ! in/out
            vacuum_aux_polemb, &                                       ! in/out
            loc_energy_term_names, loc_energy_term_values, &           ! out
            loc_n_energy_terms, &                                      ! out
            loc_pol_emb_pot_fine, loc_rep_mm_pol_emb_pot_fine, &       ! out
            density_on_fine, mdl%fine_grid, ngwf_basis(1), &           ! in
            swri_library(dma_swri_h)%s_atoms_nl, &                     ! in
            scf_dma_total_multipoles, &                                ! in
            rep%swexes(REP_SWEX_POL_EMB_DMA_1:&                        ! in
            REP_SWEX_POL_EMB_DMA_2)%dcoeffs_ht, &                      ! in
            rep%overlap%p, denskern, mdl, &                            ! in
            scf_dma_polarisation_multipoles, &                         ! OPT
            polarisable_embedding_vacuum_state%vacuum_multipoles, &    ! OPT
            polarisable_embedding_vacuum_state%rep%&
            swexes(REP_SWEX_POL_EMB_DMA_1:REP_SWEX_POL_EMB_DMA_2)%&
            dcoeffs_ht, &                                                ! OPT
            polarisable_embedding_vacuum_state%rep%&
            swexes(REP_SWEX_POL_EMB_AUX_DMA_1:REP_SWEX_POL_EMB_AUX_DMA_2)%&
            dcoeffs_ht)                                                ! OPT

    else
       call pol_emb_potential_and_energy(polemb, vacuum_polemb, &      ! in/out
            vacuum_aux_polemb, &                                       ! in/out
            loc_energy_term_names, loc_energy_term_values, &           ! out
            loc_n_energy_terms, &                                      ! out
            loc_pol_emb_pot_fine, loc_rep_mm_pol_emb_pot_fine, &       ! out
            density_on_fine, mdl%fine_grid, ngwf_basis(1), &              ! in
            swri_library(dma_swri_h)%s_atoms_nl, &                     ! in
            scf_dma_total_multipoles, &                                ! in
            rep%swexes(REP_SWEX_POL_EMB_DMA_1:&                        ! in
            REP_SWEX_POL_EMB_DMA_2)%dcoeffs_ht, &                      ! in
            rep%overlap%p, denskern, mdl)                              ! in
    end if
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------

    ! jd: Devel code hook to see what polemb pot looks like when generated
    !     in the NGWF basis. This can reveal tricky issues with representing
    !     external potentials that are difficult to capture in the NGWF basis
    if(utils_devel_code(.false.,'DEBUG','DUMP_POLEMB', pub_devel_code, &
         no_warn = .true.)) then
       if(pub_inner_loop_iteration == 1 .and. .not. pot_output_done) then
          ! jcap: the 1s here represents the region - only set up for
          ! one region at the moment
          call visual_scalarfield_ngwf_basis(polemb(1), rep, mdl, 1, &
               ngwf_basis(1), 'polemb pot (QM* induced) when represented in NGWF &
               &basis (internal units)', '_polemb_qmstar_pot_ind_'//&
               trim(utils_int_to_str(pot_output_iter)))
          call visual_scalarfield_ngwf_basis(polemb(2), rep, mdl, 1, &
               ngwf_basis(1), 'polemb pot (QM* perm) when represented in NGWF &
               &basis (internal units)', '_polemb_qmstar_pot_perm_'//&
               trim(utils_int_to_str(pot_output_iter)))
          call visual_scalarfield(loc_pol_emb_pot_fine, mdl%fine_grid, mdl%cell, &
               'total polemb pot contrib to LHXC (internal units)', &
               '_polemb_local_pot_'//trim(utils_int_to_str(pot_output_iter)))
          call visual_scalarfield(loc_rep_mm_pol_emb_pot_fine, mdl%fine_grid, &
               mdl%cell, 'total rep MM pot contrib to LHXC (internal units)', &
               '_polemb_rep_MM_pot_'//trim(utils_int_to_str(pot_output_iter)))

          call sparse_create(tmp_spam3, polemb(1))
          call sparse_scale(tmp_spam3, 0D0)
          call integrals_locpot(tmp_spam3, &
               rep%ngwfs_on_grid(1),ngwf_basis(1),rep%ngwfs_on_grid(1),ngwf_basis(1), &
               mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, &
               loc_pol_emb_pot_fine)
          call visual_scalarfield_ngwf_basis(tmp_spam3, rep, mdl, 1, &
               ngwf_basis(1), 'polemb local pot when represented in NGWF &
               &basis (internal units)', '_polemb_reconstructed_local_pot_'//&
               trim(utils_int_to_str(pot_output_iter)))
          call sparse_scale(tmp_spam3, 0D0)
          call integrals_locpot(tmp_spam3, &
               rep%ngwfs_on_grid(1),ngwf_basis(1),rep%ngwfs_on_grid(1),ngwf_basis(1), &
               mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, &
               loc_rep_mm_pol_emb_pot_fine)
          call visual_scalarfield_ngwf_basis(tmp_spam3, rep, mdl, 1, &
               ngwf_basis(1), 'polemb rep MM pot when represented in NGWF &
               &basis (internal units)', '_polemb_reconstructed_rep_MM_pot_'//&
               trim(utils_int_to_str(pot_output_iter)))
          call sparse_destroy(tmp_spam3)

          pot_output_iter = pot_output_iter + 1
          pot_output_done = .true.
       else
          pot_output_done = .false.
       end if
    end if

    call utils_assert(&
         (present(energy_term_names) .eqv. present(energy_term_values)) .and. &
         (present(energy_term_names) .eqv. present(n_energy_terms)), &
         myself//': Either all optional arguments must be provided or none')

    ! jd: If called from hamiltonian_energy_components(), return a deep copy of
    !     all components before we zero some of them
    if(present(energy_term_names)) then

       n_energy_terms = loc_n_energy_terms

       allocate(energy_term_names(n_energy_terms),stat=ierr)
       call utils_alloc_check(myself,'energy_term_names',ierr)
       allocate(energy_term_values(n_energy_terms),stat=ierr)
       call utils_alloc_check(myself,'energy_term_values',ierr)

       energy_term_names = loc_energy_term_names
       energy_term_values = loc_energy_term_values
    end if

    ! jd: Go through the polarisable embedding energy terms and partition them.
    !     Leave the original list unaltered.
    do i_term = 1, loc_n_energy_terms
       ! jd: Ensure hashed terms have been zeroed beforehand
       if(index(loc_energy_term_names(i_term),'#') /= 0) then
          call utils_assert(loc_energy_term_values(i_term) == 0D0, &
               myself//': Non-zero energy term that should be ignored: '//&
               trim(loc_energy_term_names(i_term)))
       end if

       ! jd: lhxc_energy needs to include terms that have representation in the
       !     lhxc potential, that is everything with QM elec in it. Technically
       !     it should also include MM+ pol iff polarisation is done via
       !     locpot. We keep this term separately.
       if(index(loc_energy_term_names(i_term),'QM elec <-> MM') == 1) then
          lhxc_energy = lhxc_energy + loc_energy_term_values(i_term)
          ! Note we don't zero loc_energy_term_values here, since these
          ! contributions are returned in separate energy terms
       end if

       ! jd: Include the energy due to the repulsive MM potential in lhxc
       if(index(loc_energy_term_names(i_term),'QM elec <-> rep MM') == 1) then
!         lhxc_energy = lhxc_energy + loc_energy_term_values(i_term)
!         loc_energy_term_values(i_term) = 0D0 ! sic
!         ^ Leave both commented out!
       end if

       ! jd: QM cores <-> MM perm (regardless of scheme)
       if(index(loc_energy_term_names(i_term),'QM core <-> MM perm') == 1) then
          polemb_energy_qm_core_mm_perm = &
               polemb_energy_qm_core_mm_perm + loc_energy_term_values(i_term)
          loc_energy_term_values(i_term) = 0D0
       end if
       ! jd: QM cores <-> MM ind (regardless of scheme)
       if(index(loc_energy_term_names(i_term),'QM core <-> MM ind') == 1) then
          polemb_energy_qm_core_mm_ind = &
               polemb_energy_qm_core_mm_ind + loc_energy_term_values(i_term)
          loc_energy_term_values(i_term) = 0D0
       end if

       ! --- Non-QM* (full) QM/MM electrostatics ---
       ! jd: QM elec <-> MM perm (via locpot)
       if(index(loc_energy_term_names(i_term),'QM elec <-> MM perm') == 1) then
          polemb_energy_qmfull_elec_mm_perm = &
               polemb_energy_qmfull_elec_mm_perm + loc_energy_term_values(i_term)
          loc_energy_term_values(i_term) = 0D0
       end if
       ! jd: QM elec <-> MM ind (via locpot)
       if(index(loc_energy_term_names(i_term),'QM elec <-> MM ind') == 1) then
          polemb_energy_qmfull_elec_mm_ind = &
               polemb_energy_qmfull_elec_mm_ind + loc_energy_term_values(i_term)
          loc_energy_term_values(i_term) = 0D0
       end if

       ! --- QM*/MM electrostatics ---
       ! jd: QM* elec <-> MM perm
       if(index(loc_energy_term_names(i_term),'QM* elec <-> MM perm') == 1) then
          polemb_energy_qmstar_elec_mm_perm = &
               polemb_energy_qmstar_elec_mm_perm + loc_energy_term_values(i_term)
          loc_energy_term_values(i_term) = 0D0
       end if
       ! jd: QM* elec <-> MM ind
       if(index(loc_energy_term_names(i_term),'QM* elec <-> MM ind') == 1) then
          polemb_energy_qmstar_elec_mm_ind = &
               polemb_energy_qmstar_elec_mm_ind + loc_energy_term_values(i_term)
          loc_energy_term_values(i_term) = 0D0
       end if

       ! jd: MM permanent electrostatics
       !                            don't remove this space V
       if(index(loc_energy_term_names(i_term),'MM perm mpole ') == 1) then
          polemb_energy_mm_es = &
               polemb_energy_mm_es + loc_energy_term_values(i_term)
          loc_energy_term_values(i_term) = 0D0
       end if
       ! jd: MM+ polarisation calculated externally. We lump it with mm_es, but
       !     note its presence in DFT gradients (polemb or locpot).
       if(index(loc_energy_term_names(i_term),'MM+ pol') == 1) then
          polemb_energy_mm_es = &
               polemb_energy_mm_es + loc_energy_term_values(i_term)
          loc_energy_term_values(i_term) = 0D0
       end if

       ! jd: MM bond, angle and other valence terms
       if(index(loc_energy_term_names(i_term),'MMv') == 1) then
          polemb_energy_mm_valence = &
               polemb_energy_mm_valence + loc_energy_term_values(i_term)
          loc_energy_term_values(i_term) = 0D0
       end if
       ! jd: MM vdW (calculated tinker-side)
       if(index(loc_energy_term_names(i_term),'MM vdW') == 1) then
          polemb_energy_mm_vdw = &
               polemb_energy_mm_vdw + loc_energy_term_values(i_term)
          loc_energy_term_values(i_term) = 0D0
       end if
       ! jd: QM/MM vdW (calculated tinker-side, for now)
       if(index(loc_energy_term_names(i_term),'QM/MM vdW') == 1) then
          polemb_energy_qmmm_vdw = &
               polemb_energy_qmmm_vdw + loc_energy_term_values(i_term)
          loc_energy_term_values(i_term) = 0D0
       end if

    end do

    ! jd: All remaining energy terms are lumped into 'other'.
    !     E.g. externally calculated electrostatics, tinker-side solvation, etc.
    !     MM repulsive potential too.
    polemb_energy_other = sum(loc_energy_term_values(1:loc_n_energy_terms))

    ! jd: Add the permanent embedding potential to lhxc
    ! jd: Add the repulsive MM potential to lhxc
    do is = 1, pub_num_spins
       lhxc_fine(:,:,:,is) = lhxc_fine(:,:,:,is) + &
            loc_pol_emb_pot_fine(:,:,:) + loc_rep_mm_pol_emb_pot_fine(:,:,:)
    end do

    ! jd: Deallocate the energy terms allocated in
    !     pol_emb_potential_and_energy() or calls therein
    deallocate(loc_energy_term_values, stat=ierr)
    call utils_dealloc_check(myself,'energy_term_values', ierr)
    deallocate(loc_energy_term_names, stat=ierr)
    call utils_dealloc_check(myself,'energy_term_names', ierr)
    ! jd: ... and the local potential allocated here
    deallocate(loc_rep_mm_pol_emb_pot_fine, stat=ierr)
    call utils_dealloc_check(myself,'loc_rep_mm_pol_emb_pot_fine', ierr)
    deallocate(loc_pol_emb_pot_fine, stat=ierr)
    call utils_dealloc_check(myself,'loc_pol_emb_pot_fine', ierr)

    ! jd: Destroy multipole sets, temporary swexes
    call dma_free_multipoles(scf_dma_total_multipoles)
    if(pub_pol_emb_vacuum_qmstar) then
       call dma_free_multipoles(scf_dma_polarisation_multipoles)
    end if

  end subroutine polarisable_embedding_calculate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine polarisable_embedding_expand_ngwf_pairs(rep, ngwf_basis, mdl, &
       aux_swex_h)
    !==========================================================================!
    ! Expands the NGWF pairs in 'rep' in terms of spherical waves. Expansions  !
    ! are stored in rep%swexes(:) -- in REP_SWEX_POL_EMB_DMA_1 (and, if Bessel !
    ! averaging is in effect, in REP_SWEX_POL_EMB_DMA_2) by default. This can  !
    ! be overridden via aux_swex_h.
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   rep (in/out): The NGWF_REP from which the NGWFs to expand are read, and!
    !                 where the expansion coefficients are stored.             !
    !   ngwf_basis, mdl (in): The usual.                                       !
    !   aux_swex_h (in, opt): Allows overriding the destination swex(es) in    !
    !                         rep to something else.                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2017.                                   !
    ! Modified for embedding by Joseph Prentice, September 2018                !
    !==========================================================================!

    use constants, only: REP_SWEX_POL_EMB_DMA_1
    use dma, only: dma_expand_ngwf_pairs
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_pol_emb_qmstar
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(NGWF_REP), intent(inout) :: rep
    type(FUNC_BASIS), intent(in)  :: ngwf_basis
    type(MODEL), intent(in)       :: mdl
    integer, intent(in), optional :: aux_swex_h

    ! jd: Local variables
    integer :: swex_h
    character(len=*), parameter :: myself = &
         'polarisable_embedding_expand_ngwf_pairs'

    ! -------------------------------------------------------------------------

    call utils_assert(pub_pol_emb_qmstar, myself//': Spurious call')

    if(present(aux_swex_h)) then
       swex_h = aux_swex_h
    else
       swex_h = REP_SWEX_POL_EMB_DMA_1
    end if

    call dma_expand_ngwf_pairs(rep, ngwf_basis, mdl, 1, swex_h)

!    if(present(aux_swex_h)) then
!       write(*,*) '@@@ HACK FOR MASKING CHARGES AND QPOLES in coeffs @@@'
!    do i_set = 1, merge(2,1,pub_dma_bessel_averaging)
!       call multipole_set_scale_monopoles(scf_dma_polarisation_multipoles(i_set), 0D0)
!       call multipole_set_scale_quadrupoles(scf_dma_polarisation_multipoles(i_set), 0D0)
!    end do


  end subroutine polarisable_embedding_expand_ngwf_pairs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine polarisable_embedding_ngwf_gradient(rep, &                ! in/out
       denskern, ngwf_basis, mdl, tc_denskern, &                       ! input
       cov_grad, contra_grad, &                               ! in/out (adds to)
       precond_func_recip)                                             ! in
    !==========================================================================!
    ! Driver routine for calculating the NGWF gradient contributions from      !
    ! polarisable embedding.                                                   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   rep (in/out): The NGWF_REP that contributes to the gradient. In/out    !
    !                 because NGWF caching happens behind the scenes.          !
    !   denskern, ngwf_basis, mdl (in): The usual.                             !
    !   tc_denskern (in): Tensorially-correct denskern.                        !
    !   cov_grad, cov_grad (in/out): Results are accumulated here.             !
    !   precond_func_recip (in): Needed for recip-space precionditioning.      !
    !--------------------------------------------------------------------------!
    ! Notes:                                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2017.                                   !
    ! Modified for embedding by Joseph Prentice, September 2018                !
    !==========================================================================!

    use constants, only: REP_SWEX_POL_EMB_DMA_1, REP_SWEX_POL_EMB_DMA_2
    use datatypes, only: FUNCTIONS
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use remote, only: remote_ngwf_cache_hts
    use rundat, only: pub_dma_use_ri, pub_pol_emb_vacuum_qmstar
    use sparse_array, only: SPAM3_ARRAY
    use sw_expansion_type, only: SW_EX
    use sw_resolution_of_identity, only: swri_get_handle_to
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(NGWF_REP), intent(in)             :: rep
    type(SPAM3_ARRAY), intent(in)          :: denskern
    type(FUNC_BASIS), intent(in)           :: ngwf_basis
    type(MODEL), intent(in)                :: mdl
    type(SPAM3_ARRAY), intent(in)          :: tc_denskern
    type(FUNCTIONS), intent(inout)         :: cov_grad
    type(FUNCTIONS), intent(inout)         :: contra_grad
    real(kind=DP), intent(in)              :: precond_func_recip(:,:,:)

    ! jd: Local variables
    type(SPAM3_ARRAY) :: vacuum_kdenskern_in_new_ngwfs
    type(SPAM3_ARRAY) :: tc_vacuum_kdenskern_in_new_ngwfs
    integer           :: dma_swri_h
    type(SW_EX)       :: merged_swexes(2)
    character(len=*), parameter :: myself = 'polarisable_embedding_ngwf_gradient'

    ! ------------------------------------------------------------------------

    dma_swri_h = swri_get_handle_to(pub_dma_use_ri)

    ! jd: Term (1) in NGWF gradient, due to P matrix, current K, and
    !     polarisation coefficients (if vacuum QM*) or total coefficients (otherwise)
    ! jcap: The 1 in rep%ngwfs_on_grid is a region number - only set up for one region
    call pol_emb_ngwf_grad_contribution(&
         rep%swexes(REP_SWEX_POL_EMB_DMA_1:REP_SWEX_POL_EMB_DMA_2), &
         rep%ngwfs_on_grid(1), remote_ngwf_cache_hts, &
         denskern, ngwf_basis, mdl, &
         tc_denskern, cov_grad, contra_grad, precond_func_recip, &
         rep%ngwf_cache_handle)

    ! jd: There are *no* additional terms in vacuum QM* mode once we'd abandoned
    !     the kernel rotation scheme.

  end subroutine polarisable_embedding_ngwf_gradient

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine polarisable_embedding_init_vacuum_state(vs, &
       ngwf_basis, mdl, n_occ, overlap)!, denskern)
    !==========================================================================!
    ! Initialises the vacuum state, used only if pol_emb_vacuum_qmstar T.      !
    !--------------------------------------------------------------------------!
    ! Note:                                                                    !
    !   overlap and denskern are only needed for their sparsity patterns.      !
    !   n_occ, in practice, is copied into the vacuum state from current state,!
    !   reflecting the fact that both states have the same total charge.       !
    ! JCAP: Is denskern actually used? It doesn't seem to be...                !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2017.06.                                    !
    ! Modified for embedding by Joseph Prentice, September 2018                !
    !==========================================================================!

    use constants, only: REP_N_SWEXES, SIZEOF_DOUBLE, &
         REP_SWEX_POL_EMB_DMA_1, REP_SWEX_POL_EMB_DMA_2, &
         REP_SWEX_POL_EMB_AUX_DMA_1, REP_SWEX_POL_EMB_AUX_DMA_2
    use datatypes, only: data_functions_alloc
    use dma, only: dma_expand_pair_densities, dma_calculate
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE, hash_table_init, hash_table_free, &
         hash_table_calc_cache_size
    use kernel, only: kernel_create
    use model_type, only: MODEL
    use ngwf_representation, only: ngwf_rep_create
    use restart, only: restart_ngwfs_tightbox_input
    use rundat, only: pub_rootname, pub_debug_on_root, pub_spin_fac, PUB_1K, &
         pub_dma_max_q, pub_dma_bessel_averaging, &
         pub_pol_emb_vacuum_dma_min_l, pub_pol_emb_vacuum_dma_max_l, &
         pub_pol_emb_dma_min_l, pub_pol_emb_dma_max_l, pub_dma_metric, &
         pub_dma_use_ri, pub_num_spins
    use sparse, only: SPAM3, sparse_read
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_array_read, &
         sparse_embed_array_scale, sparse_embed_extract_from_array, &
         sparse_embed_destroy_extracted_array
    use sw_expansion_type, only: SW_EX, swx_init_container, swx_coeffs_n_slots
    use sw_resolution_of_identity, only: swri_get_handle_to, swri_library
    use utils, only: utils_assert, utils_int_to_str, utils_alloc_check, &
         utils_dealloc_check

    implicit none

    ! jd: Arguments
    type(VACUUM_STATE), intent(inout) :: vs
    type(FUNC_BASIS), intent(in)      :: ngwf_basis(1)
    type(MODEL), intent(in)           :: mdl
    integer, intent(in)               :: n_occ(:,:)
    type(SPAM3_EMBED), intent(in)     :: overlap

    ! jd: Local variables
    type(HT_HASH_TABLE)         :: coeffs_ht
    integer                     :: dma_swri_h
    integer                     :: swex_h
    integer                     :: n_dma_runs, dma_run
    integer                     :: max_q
    integer                     :: natoms
    character(len=*), parameter :: myself = &
         'polarisable_embedding_init_vacuum_state'

    ! jcap: temporary arrays
    type(SPAM3), allocatable  :: kern_array(:)
    integer :: ierr

    ! -------------------------------------------------------------------------

    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering '//myself

    call utils_assert(.not. vs%initialised, &
         myself//': Vacuum state already initialised')

    ! jd ------------------- restore vacuum DKN and NGWFs -------------------

    call kernel_create(vs%kdenskern, 'K', is_cmplx=.false.)
    ! rc2013: EMBED_FIX!
    call sparse_embed_array_read(vs%kdenskern%kern, &
         trim(pub_rootname)//'.pol_emb_kdenskern')
    call ngwf_rep_create(vs%rep,'', mdl, is_cmplx = .false.)
    vs%rep%ngwf_cache_handle = 6
    call data_functions_alloc(vs%rep%ngwfs_on_grid(1), ngwf_basis(1)%size_on_grid, &
         iscmplx=.false.) ! jd: since ngwf_rep_create does not alloc ngwfs_on_grid
    call restart_ngwfs_tightbox_input(vs%rep%ngwfs_on_grid(1), &
         ngwf_basis(1), mdl%cell, mdl%fftbox, mdl%elements, &
         'pol_emb_vac_tightbox_'//ngwf_basis(1)%name, mdl%regions(1))

    call sparse_embed_array_scale(vs%kdenskern%kern, pub_spin_fac)

    ! jd: Copy n_occ from current state to vacuum state
    vs%rep%n_occ(:,:) = n_occ(:,:)

    ! rc2013: allocate temporary kernel array
    allocate(kern_array(pub_num_spins),stat=ierr)
    call utils_alloc_check(myself,'kern_array',ierr)

    ! jd ------------------- initialise vs%rep%swexes -------------------
    dma_swri_h = swri_get_handle_to(pub_dma_use_ri)
    n_dma_runs = merge(2,1,pub_dma_bessel_averaging)
    swex_h = REP_SWEX_POL_EMB_DMA_1
    max_q = pub_dma_max_q

    do dma_run = 1, n_dma_runs
       call utils_assert(.not. vs%rep%swexes(swex_h)%initialised, &
            myself//': Attempt to initialise a SW_EX that has already &
            &been initialised (VACUUM POLEMB DMA)')
       call swx_init_container(vs%rep%swexes(swex_h), &
            swri_library(dma_swri_h), &
            'vacuum_polemb_dma_swex-'//trim(utils_int_to_str(dma_run)), &
            pub_pol_emb_vacuum_dma_min_l, pub_pol_emb_vacuum_dma_max_l, &
            max_q, pub_dma_metric, vs%rep%postfix)
       max_q = max_q - 1
       swex_h = swex_h + 1
    end do

    ! jd ---------------- initialise vs%rep%swexes -- aux_swex ----------------
    swex_h = REP_SWEX_POL_EMB_AUX_DMA_1
    max_q = pub_dma_max_q

    do dma_run = 1, n_dma_runs
       call utils_assert(.not. vs%rep%swexes(swex_h)%initialised, &
            myself//': Attempt to initialise a SW_EX that has already &
            &been initialised (VACUUM POLEMB DMA -- AUX)')
       call swx_init_container(vs%rep%swexes(swex_h), &
            swri_library(dma_swri_h), &
            'vacuum_polemb_aux_dma_swex-'//trim(utils_int_to_str(dma_run)), &
            pub_pol_emb_dma_min_l, pub_pol_emb_dma_max_l, max_q, &
            pub_dma_metric, vs%rep%postfix)
       max_q = max_q - 1
       swex_h = swex_h + 1
    end do

    ! jd ------------------- expand vacuum state -------------------
    natoms = size(mdl%elements)
    call hash_table_init(coeffs_ht, 'VACUUM_COEFFS', 5, &
         swx_coeffs_n_slots(natoms), swx_coeffs_n_slots(natoms), &
         vs%rep%swexes(REP_SWEX_POL_EMB_DMA_1)%hash_table_info_unit)

    ! jcap: the 1 here corresponds to a region - only set up for 1
    ! region at the moment
    ! rc2013: EMBED_FIX!
    call sparse_embed_extract_from_array(kern_array,vs%kdenskern%kern%m(:,PUB_1K))
    call dma_expand_pair_densities(coeffs_ht, &
         vs%rep, 1, kern_array, ngwf_basis(1), mdl, &
         REP_SWEX_POL_EMB_DMA_1)
    call sparse_embed_destroy_extracted_array(kern_array)
    call dma_calculate(vs%vacuum_multipoles, 'vacuum_dma', &
         vs%rep%swexes(REP_SWEX_POL_EMB_DMA_1:REP_SWEX_POL_EMB_DMA_2), &
         vs%rep%n_occ, overlap, ngwf_basis(1), mdl, REP_SWEX_POL_EMB_DMA_1, &
         'D', coeffs_ht = coeffs_ht) ! @overlap fishy

    call hash_table_free(coeffs_ht)

    ! jd: Vacuum state dcoeffs -- full expansion
    call polarisable_embedding_expand_ngwf_pairs(vs%rep, ngwf_basis(1), mdl)

    ! jd ------------------- expand vacuum state -- aux expansion -------------------

    call hash_table_init(coeffs_ht, 'VACAUX_COEFFS', 5, &
         swx_coeffs_n_slots(natoms), swx_coeffs_n_slots(natoms), &
         vs%rep%swexes(REP_SWEX_POL_EMB_AUX_DMA_1)%hash_table_info_unit)

    ! jcap: the 1 here corresponds to a region - only set up for 1
    ! region at the moment
    call sparse_embed_extract_from_array(kern_array,vs%kdenskern%kern%m(:,PUB_1K))
    call dma_expand_pair_densities(coeffs_ht, &
         vs%rep, 1, kern_array, ngwf_basis(1), mdl, &
         REP_SWEX_POL_EMB_AUX_DMA_1)
    call sparse_embed_destroy_extracted_array(kern_array)
    call dma_calculate(vs%vacuum_aux_multipoles, 'vacuum_aux_dma', &
         vs%rep%swexes(REP_SWEX_POL_EMB_AUX_DMA_1:REP_SWEX_POL_EMB_AUX_DMA_2), &
         vs%rep%n_occ, overlap, ngwf_basis, mdl, REP_SWEX_POL_EMB_AUX_DMA_1, &
         'D', coeffs_ht = coeffs_ht) ! @overlap fishy

    call hash_table_free(coeffs_ht)

    ! jd: Vacuum state dcoeffs -- auxiliary expansion
    call polarisable_embedding_expand_ngwf_pairs(vs%rep, ngwf_basis(1), mdl, &
         REP_SWEX_POL_EMB_AUX_DMA_1)

    write(*,*) '@@@VACUUM STATE EXPANDED-2@@@'

    vs%initialised = .true.

    ! rc2013: deallocate temporary array
    deallocate(kern_array,stat=ierr)
    call utils_dealloc_check(myself,'kern_array',ierr)

    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving '//myself

  end subroutine polarisable_embedding_init_vacuum_state

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! ---------------------------------------------------------------------------
  ! ---- P R I V A T E   S U B R O U T I N E S   A N D   F U N C T I O N S ----
  ! ---------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pol_emb_potential_and_energy(polemb, vacuum_polemb, &     ! in/out
       vacuum_aux_polemb, &                                            ! in/out
       energy_term_names, energy_term_values, n_energy_terms, &        ! out
       local_pol_emb_pot, &                                            ! out
       local_rep_mm_pol_emb_pot, &                                     ! out
       density_on_grid, grid, ngwf_basis, s_atoms_nl, &                ! in
       scf_dma_driving_multipoles, driving_dcoeffs_ht, &               ! in
       overlap, denskern, mdl, &                                       ! in
       scf_dma_polarisation_multipoles, vacuum_multipoles, &           ! OPT
       vacuum_dcoeffs_ht, vacuum_aux_dcoeffs_ht)                       ! OPT
    !==========================================================================!
    ! Go away and use the driver routine, polarisable_embedding_calculate().   !
    !--------------------------------------------------------------------------!
    ! Calculates the electrostatic potential and electrostatic energy resulting!
    ! from the interaction of ONETEP density (elec+core) with an external      !
    ! embedding potential, represented as sets of multipoles. The multipoles   !
    ! are assumed to be changing values (technically changing positions are    !
    ! supported too) between invocations of this subroutine, allowing for the  !
    ! embedding to respond to changes in ONETEP density.                       !
    ! The positions and values of the multipoles are read from a file. So are  !
    ! extra ("external") energy terms that can later be added to the ONETEP    !
    ! Hamiltonian. These can include, e.g. classical bonded terms, vdW, or     !
    ! externally computed approximations to the interaction of the embedding   !
    ! with ONETEP density.                                                     !
    !                                                                          !
    ! The potential of all permanent multipole sets is returned in             !
    ! local_pol_emb_pot.                                                       !
    ! The energy terms read from the file are returned in energy_term_values,  !
    ! and their names -- in energy_term_names. Both arrays are n_energy_terms  !
    ! long. Two extra energy terms are *always* placed in these arrays for     !
    ! every multipole set -- these are:                                        !
    ! 1) the interaction energy between the multipole set and ONETEP electrons,!
    ! 2) the interaction energy between the multipole set and ONETEP cores.    !
    !                                                                          !
    ! The extra (3*nsets) energy terms are always the first entries in the     !
    ! arrays and are included in the n_energy_terms counter. However, the user !
    ! can wish these not to be included in the Hamiltonian, by specifying the  !
    ! damping strategy of some, or all, multipole sets as DAMPING_ZERO. In this!
    ! case these terms are still calculated, and the entries are present, but  !
    ! are set to zero.                                                         !
    !                                                                          !
    ! All external energy terms whose names contain '#' are *NOT* added to the !
    ! Hamiltonian.                                                             !
    !                                                                          !
    ! The arrays energy_term_names and energy_term_values are allocated here,  !
    ! it is the responsibility of the caller to later destroy them.            !
    ! Local_emb_pot, in contrast, is assumed to have been allocated beforehand.!
    !                                                                          !
    ! The argument 'polemb' is expected to contain a pair of S-sparsity pattern!
    ! SPAM3's on entry, and on exit it will contain the polarisability ('P')   !
    ! matrix needed in gradients. Element 1 is for induced interactions.       !
    ! Element 2 is for permanent interactions if using QM* representation for  !
    ! these. If not, then it is zero.                                          !
    !                                                                          !
    ! The role of arguments density_on_grid, grid, ngwf_basis, overlap,        !
    ! denskern and elements is obvious, they are all intent(in).               !
    !                                                                          !
    ! s_atoms_nl (in) is used to determine S-neighbours.                       !
    !                                                                          !
    ! scf_dma_multipoles (in), as prepared by dma_evaluate_multipoles(), is    !
    ! needed to calculate the interaction energy of polarisable dipoles with   !
    ! DMA multipoles. Also, in the calculation of the P matrix the DMA monopo- !
    ! les are needed and their scaling factor, if any. It is stored along with !
    ! all multipoles in scf_dma_multipoles.                                    !
    !                                                                          !
    ! pol_dcoeffs_ht is a hash table containing the d coefficients needed for the  !
    ! calculation of P. This hash table is prepared in DMA.                    !
    ! @vastly updatedoc (also inv_overlap) vacuum_polemb
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2015.05.                                    !
    ! Gradient correction term added by Jacek Dziedzic in 2015.09.             !
    ! Reworked to get rid of the correction term, Jacek Dziedzic 2015.11.      !
    ! Modified to remove pub_par by Joseph Prentice, September 2018            !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root, comms_reduce
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE
    use integrals, only: integrals_product_on_grid
    use model_type, only: MODEL
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET, &
         multipole_init_spherical_set, multipole_free_spherical_set, &
         DAMPING_COULOMBIC_MASKED, DAMPING_COULOMBIC_SMEARED, DAMPING_THOLE, &
         DAMPING_ZERO, DAMPING_FROM_POTENTIAL
    use neighbour_list, only: NL_NEIGHBOUR_LIST
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_pol_emb_pot_filename, pub_pol_emb_pot, &
         pub_num_spins, pub_output_detail, pub_pol_emb_perm_scaling, &
         pub_dma_multipole_scaling, pub_pol_emb_fixed_charge, &
         pub_pol_emb_vacuum_qmstar, pub_pol_emb_repulsive_mm_pot_alpha, &
         pub_pol_emb_repulsive_mm_pot_verbose
    use sparse, only: SPAM3, sparse_create, sparse_axpy, sparse_destroy, &
         sparse_scale, sparse_get_par
    use sparse_array, only: SPAM3_ARRAY
    use sw_expansion_type, only: SW_EX
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_abort, utils_dealloc_check, &
         utils_alloc_check, utils_flushed_string_output, utils_banner, &
         utils_real_to_str

    implicit none

    ! jd: Arguments
    type(SPAM3), intent(inout)                  :: polemb(:)
    type(SPAM3), intent(inout)                  :: vacuum_polemb(:)
    type(SPAM3), intent(inout)                  :: vacuum_aux_polemb(:)
    character(len=27), intent(out), allocatable :: energy_term_names(:)
    real(kind=DP), intent(out), allocatable     :: energy_term_values(:)
    integer, intent(out)                        :: n_energy_terms
    type(GRID_INFO), intent(in)                 :: grid
    real(kind=DP), intent(out)                  :: local_pol_emb_pot(&
         grid%ld1, grid%ld2, grid%max_slabs12)
    real(kind=DP), intent(out)                  :: local_rep_mm_pol_emb_pot(&
         grid%ld1, grid%ld2, grid%max_slabs12)
    real(kind=DP), intent(in)                   :: density_on_grid(:,:,:,:)
    type(FUNC_BASIS), intent(in)                :: ngwf_basis
    type(NL_NEIGHBOUR_LIST), intent(in)         :: s_atoms_nl
    type(SPHERICAL_MULTIPOLE_SET), intent(in)   :: scf_dma_driving_multipoles(:)
    type(HT_HASH_TABLE), intent(in), target     :: driving_dcoeffs_ht(:)
    type(SPAM3), intent(in)                     :: overlap
    type(SPAM3_ARRAY), intent(in)               :: denskern
    type(MODEL), intent(in)                     :: mdl
    type(SPHERICAL_MULTIPOLE_SET), intent(in), optional :: scf_dma_polarisation_multipoles(:)
    type(SPHERICAL_MULTIPOLE_SET), intent(in), optional :: vacuum_multipoles(:)
    type(HT_HASH_TABLE), intent(in), target, optional   :: vacuum_dcoeffs_ht(:)
    type(HT_HASH_TABLE), intent(in), target, optional   :: vacuum_aux_dcoeffs_ht(:)

    ! jd: Local variables
    type(SPHERICAL_MULTIPOLE_SET), allocatable :: sph_mm_multipoles(:)
    type(SPAM3)       :: polemb_cur_set
    type(SPAM3)       :: vacuum_polemb_cur_set
    type(SPAM3)       :: vacuum_aux_polemb_cur_set
    integer           :: n_sets
    integer           :: mpole_set
    integer           :: ierr
    integer           :: is
    integer           :: last_repulsive_set
    real(kind=DP)     :: pot_max, pot_min, pot_norm
    real(kind=DP)     :: external_permanent, external_induced
    real(kind=DP)     :: internal_permanent, internal_induced
    real(kind=DP)     :: difference_permanent, difference_induced
    real(kind=DP), allocatable :: local_pol_emb_pot_cur_set(:,:,:)
    real(kind=DP), allocatable :: temp_rep_mm_pol_emb_pot(:,:,:)
    real(kind=DP)     :: scal_factor
    real(kind=DP)     :: ind_factor
    logical           :: set_is_fixed
    logical           :: set_is_permanent
    logical           :: set_is_induced
    logical           :: set_is_qmstar
    logical           :: set_is_repulsive
    logical           :: induction_absent, permanent_absent
    logical           :: induction_damped, permanent_damped
    logical           :: induction_from_qmstar, permanent_from_qmstar
    integer           :: i_term
    logical           :: dont_compare_perm_with_tinker
    logical           :: dont_compare_ind_with_tinker
    integer           :: mm_atom
    real(kind=DP)     :: mm_rep_energy
    type(PARAL_INFO), pointer :: par

    ! Pretty printing
    character(len=6)  :: use_string
    character(len=6)  :: detail_string
    character(len=80) :: outstr
    real(kind=DP), parameter :: difference_induced_threshold = 1.5936D-5 ! 0.01 kcal/mol
    real(kind=DP), parameter :: difference_permanent_threshold = 1.5936D-5 ! 0.01 kcal/mol
    character(len=*), parameter :: myself = 'pol_emb_potential_and_energy'

    ! ---------------------------------------------------------------------------

    call timer_clock(myself,1)

    call utils_assert(pub_pol_emb_pot, myself//&
         ': Cannot call this subroutine without "pol_emb_pot_filename"')

    call utils_assert(pub_pol_emb_vacuum_qmstar .eqv. &
         present(scf_dma_polarisation_multipoles), &
         myself//': optional arg scf_dma_polarisation_multipoles &
         &must be passed iff pub_pol_emb_vacuum_qmstar')
    call utils_assert(pub_pol_emb_vacuum_qmstar .eqv. &
         present(vacuum_multipoles), &
         myself//': optional arg vacuum_multipoles &
         &must be passed iff pub_pol_emb_vacuum_qmstar')
    call utils_assert(pub_pol_emb_vacuum_qmstar .eqv. &
         present(vacuum_dcoeffs_ht), &
         myself//': optional arg vacuum_dcoeffs_ht &
         &must be passed iff pub_pol_emb_vacuum_qmstar')
    call utils_assert(pub_pol_emb_vacuum_qmstar .eqv. &
         present(vacuum_aux_dcoeffs_ht), &
         myself//': optional arg vacuum_aux_dcoeffs_ht &
         &must be passed iff pub_pol_emb_vacuum_qmstar')

    call sparse_get_par(par, overlap)

    ! jd: Read Cartesian multipoles from file, convert to sphericals
    call pol_emb_read_from_file(energy_term_names, energy_term_values, &
         n_energy_terms, sph_mm_multipoles, n_sets, pub_pol_emb_pot_filename)

    ! jd: Must have at least the 3 energy terms per set that are calculated
    !     internally in ONETEP
    call utils_assert(n_energy_terms >= 3*n_sets, myself//': Internal error')

    allocate(local_pol_emb_pot_cur_set(grid%ld1, grid%ld2, grid%max_slabs12), &
         stat=ierr)
    call utils_alloc_check(myself,'local_pol_emb_pot_cur_set',ierr)

    if(pub_pol_emb_repulsive_mm_pot_verbose) then
       allocate(temp_rep_mm_pol_emb_pot(grid%ld1, grid%ld2, grid%max_slabs12), &
            stat=ierr)
       call utils_alloc_check(myself,'temp_rep_mm_pol_emb_pot',ierr)
    end if

    ! jd: Local potentials and polemb matrices start with zero
    local_pol_emb_pot(:,:,:) = 0D0
    local_rep_mm_pol_emb_pot(:,:,:) = 0D0

    call sparse_scale(polemb(1),0D0)
    call sparse_scale(polemb(2),0D0)
    call sparse_create(polemb_cur_set,polemb(1))

    if(pub_pol_emb_vacuum_qmstar) then
       call sparse_scale(vacuum_polemb(1),0D0)
       call sparse_scale(vacuum_polemb(2),0D0)
       call sparse_scale(vacuum_aux_polemb(1),0D0)
       call sparse_scale(vacuum_aux_polemb(2),0D0)
       call sparse_create(vacuum_polemb_cur_set,vacuum_polemb(1))
       call sparse_create(vacuum_aux_polemb_cur_set,vacuum_polemb(1))
    end if

    induction_absent = .true.
    permanent_absent = .true.
    induction_from_qmstar = .false.
    permanent_from_qmstar = .false.
    induction_damped = .false.
    permanent_damped = .false.

    ! ==========================================================================
    ! jd: Find out the index of the last repulsive set
    ! ==========================================================================
    last_repulsive_set = -1
    do mpole_set = 1, n_sets
       set_is_repulsive = &
            (index(sph_mm_multipoles(mpole_set)%set_name,'rep')/=0)
       if(set_is_repulsive) then
          last_repulsive_set = mpole_set
       end if
    end do

    ! ==========================================================================
    ! jd: Go over all multipole sets and calculate energies and potentials
    ! ==========================================================================
    do mpole_set = 1, n_sets

       ! jd: Each set starts with zero local potential and zero P matrix
       local_pol_emb_pot_cur_set = 0D0
       call sparse_scale(polemb_cur_set,0D0)
       if(pub_pol_emb_vacuum_qmstar) then
          call sparse_scale(vacuum_polemb_cur_set,0D0)
          call sparse_scale(vacuum_aux_polemb_cur_set,0D0)
       end if

       ! jd: Figure out the attributes of current set
       set_is_permanent = &
            (index(sph_mm_multipoles(mpole_set)%set_name,'perm')/=0)
       set_is_induced = &
            (index(sph_mm_multipoles(mpole_set)%set_name,'ind')/=0)
       set_is_fixed = &
            (index(sph_mm_multipoles(mpole_set)%set_name,'fix')/=0)
       set_is_qmstar = &
            (index(sph_mm_multipoles(mpole_set)%set_name,'qmstar')/=0)
       set_is_repulsive = &
            (index(sph_mm_multipoles(mpole_set)%set_name,'rep')/=0)

       ! jd: This will break down if there is ever more than one permanent
       !     or more than one induced set, hence the assert
       if(set_is_permanent .and. set_is_qmstar) permanent_from_qmstar = .true.
       if(set_is_induced .and. set_is_qmstar) induction_from_qmstar = .true.
       call utils_assert(MAX_N_MULTIPOLE_SETS <= 2, &
            myself//': Logic not entirely ready for more than 2 multipole sets')
       call utils_assert(set_is_induced .neqv. set_is_permanent, &
            myself//': Each multipole set must be either permanent or induced')

       if(set_is_induced) then
          call utils_assert(sph_mm_multipoles(mpole_set)%potential_damping /= &
               DAMPING_COULOMBIC_MASKED .and. &
               sph_mm_multipoles(mpole_set)%potential_damping /= &
               DAMPING_COULOMBIC_SMEARED, &
               myself//': Cannot have Coulombic induced dipoles, &
               &that would lead to polarisation catastrophe. Also TINKER''s &
               &MM+ term always responds to *Thole-damped* QM multipoles, so &
               &that would be inconsistent too.')
       end if

       ! jd: Induced dipoles have a 1/2 factor from the work cost of inducing
       !     the dipoles (cf. Simmonett paper or TINKTEP paper).
       if(set_is_induced) then
          ind_factor = 0.5_DP
       else
          ind_factor = 1.0_DP
       end if

       ! jd: Optional additional scaling factors.
       scal_factor = pub_dma_multipole_scaling
       if(set_is_permanent) scal_factor = scal_factor * pub_pol_emb_perm_scaling

       if(set_is_qmstar) then
          ! -----------------------------------------------------------------
          ! jd: QM* sets need to calculate P matrix every time, even if set
          !     is fixed.
          ! -----------------------------------------------------------------
          if(pub_pol_emb_vacuum_qmstar) then
             call pol_emb_qmstar_energy_and_lnv_gradient(&
                  energy_term_values(3*mpole_set-1), polemb_cur_set,  &
                  driving_dcoeffs_ht, sph_mm_multipoles(mpole_set), &
                  scf_dma_driving_multipoles, &
                  ngwf_basis, mdl, denskern, overlap, s_atoms_nl, &
                  scal_factor, ind_factor, &
                  scf_dma_polarisation_multipoles, vacuum_multipoles, &
                  vacuum_polemb_cur_set, vacuum_aux_polemb_cur_set, &
                  vacuum_dcoeffs_ht, vacuum_aux_dcoeffs_ht)
             ! jd: Accumulate total vacuum_polemb
             call sparse_axpy(vacuum_polemb(merge(1,2,set_is_induced)), &
                  vacuum_polemb_cur_set,1D0)
             call sparse_axpy(vacuum_aux_polemb(merge(1,2,set_is_induced)), &
                  vacuum_aux_polemb_cur_set,1D0)
          else
             call pol_emb_qmstar_energy_and_lnv_gradient(&
                  energy_term_values(3*mpole_set-1), polemb_cur_set, &
                  driving_dcoeffs_ht, sph_mm_multipoles(mpole_set), &
                  scf_dma_driving_multipoles, &
                  ngwf_basis, mdl, denskern, overlap, s_atoms_nl, &
                  scal_factor, ind_factor)
          end if
          ! jd: Accumulate total polemb
          call sparse_axpy(polemb(merge(1,2,set_is_induced)),polemb_cur_set,1D0)
       else
          ! -----------------------------------------------------------------
          ! jd: Non-QM* sets calculate local potential, or look it up if set
          !     is fixed.
          ! -----------------------------------------------------------------
          if(set_is_fixed .and. fixed_sets_calculated) then
             ! jd: If this is a fixed set, and it has been dealt with already,
             !     just look up the local potential.
             local_pol_emb_pot_cur_set(:,:,:) = &
                  fixed_sets_polemb_locpot(:,:,:,mpole_set)
          else
             ! jd: Else calculate it.
             call pol_emb_potential_rspc(local_pol_emb_pot_cur_set, &
                  sph_mm_multipoles(mpole_set), grid)
          end if
          ! jd: And store it, if fixed.
          if(set_is_fixed) then
             fixed_sets_polemb_locpot(:,:,:,mpole_set) = &
                  local_pol_emb_pot_cur_set(:,:,:)
          end if

          ! jd: apply factors to local potential.
          local_pol_emb_pot_cur_set = &
               local_pol_emb_pot_cur_set * scal_factor * ind_factor

          ! jd: add local potential from this set to total
          local_pol_emb_pot(:,:,:) = local_pol_emb_pot(:,:,:) + &
               local_pol_emb_pot_cur_set(:,:,:)

          ! -----------------------------------------------------------------
          ! jd: Calculate the energy of the multipole set interacting with
          !     the electrons, from integrating the local potential with
          !     density. This needs to be recalculated even for fixed sets,
          !     as obviously electronic dofs (density) change.
          ! -----------------------------------------------------------------
          energy_term_values(3*mpole_set-1) = 0D0
          do is = 1, pub_num_spins
             energy_term_values(3*mpole_set-1) = &
                  energy_term_values(3*mpole_set-1) + &
                  integrals_product_on_grid(grid, &
                  local_pol_emb_pot_cur_set(:,:,:), density_on_grid(:,:,:,is))
          end do

       end if ! QM* or non-QM*

       ! ---------------------------------------------------------------------
       ! jd: If a set uses the MM repulsive pot, calculate it or look it up
       ! ---------------------------------------------------------------------
       energy_term_values(3*mpole_set-2) = 0D0
       if(set_is_repulsive) then
          ! jd: @Assumes MM ions do not move@
          if(.not. repulsive_mm_pot_calculated) then
             if(mpole_set == last_repulsive_set) then
                ! jd: Last set. Calculate and store repulsive potential
                if(pub_pol_emb_repulsive_mm_pot_alpha > 0D0) then
                   call pol_emb_repulsive_mm_pot(local_rep_mm_pol_emb_pot, &
                        sph_mm_multipoles(mpole_set), grid, mdl%cell)
                else
                   local_rep_mm_pol_emb_pot(:,:,:) = 0D0
                end if
                stored_local_rep_mm_pol_emb_pot(:,:,:) = local_rep_mm_pol_emb_pot(:,:,:)
                repulsive_mm_pot_calculated = .true.
             end if
          else
             ! jd: Use the stored repulsive pot, attach to last repulsive set
             if(mpole_set == last_repulsive_set) then
                local_rep_mm_pol_emb_pot(:,:,:) = stored_local_rep_mm_pol_emb_pot(:,:,:)
             end if
          end if ! calculate or restore repulsive MM pot

          if(mpole_set == last_repulsive_set) then
             ! -----------------------------------------------------------------
             ! jd: Calculate the energy of the repulsive MM pot interacting with
             !     the electrons, from integrating the local potential with
             !     density, regardless if the potential has just been calculated
             !     or reused. Electronic DOFs changed in any case.
             ! -----------------------------------------------------------------
             do is = 1, pub_num_spins
                energy_term_values(3*mpole_set-2) = &
                     energy_term_values(3*mpole_set-2) + &
                     integrals_product_on_grid(grid, &
                     local_rep_mm_pol_emb_pot(:,:,:), &
                     density_on_grid(:,:,:,is))
             end do

             ! jd: If verbose MM potentials used, calculate and output contribu-
             !     tions from each MM atom separately
             if(pub_pol_emb_repulsive_mm_pot_verbose) then
                do mm_atom = 1, sph_mm_multipoles(mpole_set)%n_sites
                   ! jd: Skip dummy QM multipoles masquerading as MM multipoles
                   if(sph_mm_multipoles(mpole_set)%max_l_mask(mm_atom) == -1) cycle
                   ! jd: Calculate MM pot from *just this MM atom*
                   call pol_emb_repulsive_mm_pot(temp_rep_mm_pol_emb_pot, &
                        sph_mm_multipoles(mpole_set), grid, mdl%cell, mm_atom)
                   ! jd: Integrate with density, to get energy
                   mm_rep_energy = 0D0
                   do is = 1, pub_num_spins
                      mm_rep_energy = &
                           mm_rep_energy + &
                           integrals_product_on_grid(grid, &
                           temp_rep_mm_pol_emb_pot(:,:,:), &
                           density_on_grid(:,:,:,is))
                   end do
                   if(pub_on_root) then
                      write(*,'(a,i7,a5,f10.3,f10.3,f10.3,f16.9,a)') &
                           '| => MM rep contrib ', mm_atom, &
                           sph_mm_multipoles(mpole_set)%species(mm_atom), &
                           sph_mm_multipoles(mpole_set)%centres(mm_atom)%x, &
                           sph_mm_multipoles(mpole_set)%centres(mm_atom)%y, &
                           sph_mm_multipoles(mpole_set)%centres(mm_atom)%z, &
                           mm_rep_energy, ' |'
                   end if
                end do
             end if
          end if

       end if ! repulsive MM set

       ! --------------------------------------------------------------------
       ! jd: Calculate the energy of the multipole set interacting with cores
       !     regardless of whether set is QM* or not.
       ! --------------------------------------------------------------------
       if(set_is_fixed .and. fixed_sets_calculated) then
          ! jd: If this is a fixed set, and we have the energy already,
          !     just look it up. This is done not only for efficiency, but
          !     is needed for correct operation if the multipoles in the file
          !     actually do change, and we want them fixed.
          energy_term_values(3*mpole_set) = &
               fixed_core_energy_terms(mpole_set)
       else
          ! jd: Else calculate it
          energy_term_values(3*mpole_set) = &
               pol_emb_core_energy_rspc(sph_mm_multipoles(mpole_set), mdl%elements, par) * &
               scal_factor * ind_factor

          ! jd: and store, if set is fixed
          if(set_is_fixed) then
             fixed_core_energy_terms(mpole_set) = &
                  energy_term_values(3*mpole_set)
          end if
       end if

       ! --------------------------------------------------------------------
       ! jd: User feedback
       ! --------------------------------------------------------------------

       ! jd: Information on the values of the internally calculated terms
       if(pub_output_detail >= VERBOSE .and. pub_on_root) then
          select case(sph_mm_multipoles(mpole_set)%energy_damping)
          case (DAMPING_FROM_POTENTIAL)
             select case(sph_mm_multipoles(mpole_set)%potential_damping)
             case (DAMPING_COULOMBIC_MASKED)
                use_string = 'COUL-M'
             case (DAMPING_COULOMBIC_SMEARED)
                use_string = 'COUL-S'
             case (DAMPING_ZERO)
                use_string = '   NO!'
             case (DAMPING_THOLE)
                use_string = ' THOLE'
             case default
                call utils_abort(myself//': Unrecognized potential damping')
             end select
          case (DAMPING_ZERO)
             use_string = '    NO'
          case default
             call utils_abort(myself//': Unrecognized energy damping')
          end select

          ! jd: Take note if some interactions are absent
          if(.not. (sph_mm_multipoles(mpole_set)%energy_damping == DAMPING_ZERO&
               .or. (sph_mm_multipoles(mpole_set)%energy_damping == &
               DAMPING_FROM_POTENTIAL .and. &
               sph_mm_multipoles(mpole_set)%potential_damping == DAMPING_ZERO &
               ))) then
             if(set_is_permanent) permanent_absent = .false.
             if(set_is_induced) induction_absent = .false.
          end if

          ! jd: Take note if some interactions are damped
          if(sph_mm_multipoles(mpole_set)%energy_damping == DAMPING_THOLE&
               .or. (sph_mm_multipoles(mpole_set)%energy_damping == &
               DAMPING_FROM_POTENTIAL .and. &
               sph_mm_multipoles(mpole_set)%potential_damping == DAMPING_THOLE &
               )) then
             if(set_is_permanent) permanent_damped = .true.
             if(set_is_induced) induction_damped = .true.
          end if

          detail_string = &
               merge('P','I',set_is_permanent)//&
               merge('*',' ',set_is_qmstar)//&
               merge('F',' ',set_is_fixed)//&
               merge(merge('r','1',fixed_sets_calculated),' ',set_is_fixed &
                    .and. .not. set_is_qmstar)//&
               merge('S',' ',pub_dma_multipole_scaling /= 1D0)//&
               merge('s',' ',set_is_permanent .and. &
                    pub_pol_emb_perm_scaling /= 1D0)

          ! jd: elec <-> rep MM
          if(set_is_repulsive) then
             write(outstr,'(a,f18.9,a)') &
                  '| - '//energy_term_names(3*mpole_set-2)//' ', &
                  energy_term_values(3*mpole_set-2), ' Ha  ONETEP   '//&
                  'REPULS'//'  '//detail_string//' |'
             write(stdout,'(a80)') outstr
          else
             call utils_assert(energy_term_values(3*mpole_set-2) == 0D0, &
                  myself//': Unexpected MM rep term for a set that is not &
                  &repulsive')
          end if

          ! jd: elec <-> MM
          write(outstr,'(a,f18.9,a)') &
               '| - '//energy_term_names(3*mpole_set-1)//' ', &
               energy_term_values(3*mpole_set-1), ' Ha  ONETEP   '//&
               use_string//'  '//detail_string//' |'
          write(stdout,'(a80)') outstr

          ! jd: Cores interacting with fixed sets always reuse energy,
          !     not just potential.
          if(set_is_fixed) then
             if(fixed_sets_calculated) then
                detail_string(4:4) = 'R'
             else
                detail_string(4:4) = '1'
             end if
          end if

          ! jd: core <-> MM
          write(outstr,'(a,f18.9,a)') &
               '| - '//energy_term_names(3*mpole_set)//' ', &
               energy_term_values(3*mpole_set), ' Ha  ONETEP   '//&
               use_string//'  '//detail_string//' |'
          write(stdout,'(a80)') outstr

       end if ! verbose output on root

    end do ! mpole sets

    fixed_sets_calculated = .true.

    ! jd: Provide comparison of the QM/MM energy terms between the full
    !     density and the QMstar model
    call utils_flushed_string_output(utils_banner('-', delimiter='|') // CRLF)

    if(pub_on_root) then
       external_permanent = internal_eterm_by_name('QM/MM perm mpole')
       internal_permanent = & ! 'QM' dropped here, as this must also match 'QM*'
            internal_eterm_by_name('elec <-> MM perm') + &
            internal_eterm_by_name('core <-> MM perm')
       external_induced = internal_eterm_by_name('QM/MM polarisation')
       if(.not. pub_pol_emb_fixed_charge) then
          internal_induced = & ! 'QM' dropped here, as this must also match 'QM*'
               internal_eterm_by_name('elec <-> MM ind') + &
               internal_eterm_by_name('core <-> MM ind')
       else
          internal_induced = 0D0
       end if

       if(pub_output_detail >= VERBOSE) then
          write(stdout,'(a80)') '|                      External |   &
               &        Internal |         Difference      |'
          write(outstr, '(a13,f18.9,a3,f18.9,a3,f18.12,a7)') &
               '| Permanent: ', external_permanent, ' | ', internal_permanent, &
               ' | ', external_permanent-internal_permanent,' (Ha) |'
          write(stdout,'(a80)') outstr
          if(pub_dma_multipole_scaling /= 1D0) then
             write(outstr, '(a13,f18.9,a3,f18.9,a3,f18.12,a7)') &
                  '| ^with DMS: ', &
                  external_permanent/pub_dma_multipole_scaling, ' | ', &
                  internal_permanent, ' | ', &
                  external_permanent/pub_dma_multipole_scaling-&
                  internal_permanent,' (Ha) |'
             write(stdout,'(a80)') outstr
          end if
          if(pub_pol_emb_perm_scaling /= 1D0) then
             write(outstr, '(a13,f18.9,a3,f18.9,a3,f18.12,a7)') &
                  '| ^with PEPS:', &
                  external_permanent*pub_pol_emb_perm_scaling, ' | ', &
                  internal_permanent, ' | ', &
                  external_permanent*pub_pol_emb_perm_scaling-&
                  internal_permanent,' (Ha) |'
             write(stdout,'(a80)') outstr
          end if
          if(pub_dma_multipole_scaling /= 1D0 .and. &
               pub_pol_emb_perm_scaling /= 1D0) then
             write(outstr, '(a13,f18.9,a3,f18.9,a3,f18.12,a7)') &
                  '| ^with both:', &
                  external_permanent*pub_pol_emb_perm_scaling/&
                  pub_dma_multipole_scaling, ' | ', &
                  internal_permanent, ' | ', &
                  external_permanent*pub_pol_emb_perm_scaling/&
                  pub_dma_multipole_scaling-&
                  internal_permanent,' (Ha) |'
             write(stdout,'(a80)') outstr
          end if
          difference_permanent = &
               external_permanent*pub_pol_emb_perm_scaling-internal_permanent
          difference_induced = external_induced-internal_induced
          write(outstr, '(a13,f18.9,a3,f18.9,a3,f18.12,a7)') &
               '| Induced:   ', external_induced, ' | ', internal_induced, &
               ' | ', difference_induced,' (Ha) |'
          write(stdout,'(a80)') outstr
       end if

       if(permanent_absent) then
          write(stdout,'(a80)') '| (perm term absent: energy_zero or &
               &potential_zero+energy_from_potential used) |'
       else if(permanent_damped) then
          write(stdout,'(a80)') '| (perm term damped: will not exactl&
               &y match TINKER''s idea)                     |'
       else if(permanent_from_qmstar) then
          if(abs(difference_permanent) > difference_permanent_threshold) then
             call utils_abort(myself//': Serious discrepancy detected betwee&
                  &n TINKER''s calculation of QM/MM permanent interacti&
                  &on energy and ONETEP''s idea. The discrepancy is '//&
                  trim(utils_real_to_str(difference_permanent))//' Ha')
          end if
       end if

       ! jd: If induction turned off, notify user. If turned on, check if there
       !     is no discrepancy between us and TINKER.
       if(induction_absent) then
          if(pub_pol_emb_fixed_charge) then
             write(stdout,'(a80)') '| (induction absent: pol_emb_fixed_charge &
                  &used)                                      |'
          else
             write(stdout,'(a80)') '| (induction absent: energy_zero or &
                  &potential_zero+energy_from_potential used) |'
          end if
          if(internal_eterm_by_name('MM+ pol',.true.) /= 0D0 .and. &
               internal_eterm_by_name('#MM+ pol',.true.) == 0D0) then
             call utils_abort('MM+ polarisation NOT disabled (tinker_mm_pol_&
                  &energy 1), but induction ignored in ONETEP. This is not &
                  &physically correct. The MM+ pol term must be accounted for &
                  &in the gradient in ONETEP.')
          end if
       else
          if(induction_from_qmstar) then
             if(abs(difference_induced) > difference_induced_threshold) then
                call utils_abort(myself//': Serious discrepancy detected betwee&
                     &n TINKER''s calculation of QM/MM induced dipole interacti&
                     &on energy and ONETEP''s idea. The discrepancy is '//&
                     trim(utils_real_to_str(difference_induced))//' Ha')
             end if
          end if
       end if

    end if ! pub_on_root

    ! jd: Clean up local buffers used for verbose calcs
    if(pub_pol_emb_repulsive_mm_pot_verbose) then
       deallocate(temp_rep_mm_pol_emb_pot, stat=ierr)
       call utils_dealloc_check(myself,'temp_rep_mm_pol_emb_pot',ierr)
    end if

    ! jd: Clean up locpot array and temp P matrix
    deallocate(local_pol_emb_pot_cur_set, stat=ierr)
    call utils_dealloc_check(myself,'local_pol_emb_pot_cur_set',ierr)
    call sparse_destroy(polemb_cur_set)
    if(pub_pol_emb_vacuum_qmstar) then
       call sparse_destroy(vacuum_polemb_cur_set)
       call sparse_destroy(vacuum_aux_polemb_cur_set)
    end if

    ! -----------------------------------------------------------------------
    ! jd: Calculate some metrics of the total potential
    ! -----------------------------------------------------------------------
    if(pub_output_detail >= VERBOSE) then
       pot_max = maxval(local_pol_emb_pot(:,:,:))
       pot_min = minval(local_pol_emb_pot(:,:,:))
       call comms_reduce('MAX',pot_max)
       call comms_reduce('MIN',pot_min)
       pot_norm = sqrt(integrals_product_on_grid(grid, &
            local_pol_emb_pot,local_pol_emb_pot) / grid%weight / &
            grid%n1 / grid%n2 / grid%n3)
       if(pub_on_root) then
          write(stdout,'(a)') '|-----------------------------&
               &-------------------------------------------------|'
          write(stdout,'(a,e11.4,a,e11.4,a,e11.4,a)') &
               '| Perm embed. potential  min: ', pot_min,'  max: ', pot_max, &
               '  norm: ', pot_norm,' |'
       end if
    end if

    ! jd: Zero terms that are not meant to be included. This is needed for the
    !     case where we generate the potential (to drive the electronic dof),
    !     but throw away the energies calculated, using the external values
    !     instead. This happens when energy_damping == DAMPING_ZERO, and
    !     potential_damping is /= DAMPING_ZERO.
    do mpole_set = 1, n_sets
       if(sph_mm_multipoles(mpole_set)%energy_damping == DAMPING_ZERO) then
          ! jd: Note that the repulsive term (3*mpole_set-2) is NOT zeroed.
          energy_term_values(3*mpole_set-1) = 0D0
          energy_term_values(3*mpole_set) = 0D0
       end if
    end do

    ! jd: Zero terms that are not meant to be included and are denoted with #
    do i_term = 1, n_energy_terms
       if(index(energy_term_names(i_term),'#') /= 0) then
          energy_term_values(i_term) = 0D0
       end if
    end do

    ! jd: Destroy spherical representation of each set
    do mpole_set = 1, n_sets
       call multipole_free_spherical_set(sph_mm_multipoles(mpole_set))
    end do

    ! jd: Destroy the array of sets itself
    deallocate(sph_mm_multipoles,stat=ierr)
    call utils_dealloc_check(myself,'sph_multipoles',ierr)

    call utils_flushed_string_output(utils_banner('~', delimiter='\\') // CRLF)

    call timer_clock(myself,2)

    return

contains

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  real(kind=DP) function internal_eterm_by_name(name, allowed_to_be_absent)
    !==========================================================================!
    ! Returns a value from energy_term_values corresponding to a name in       !
    ! energy_term_names. If the energy term is absent, aborts, unless          !
    ! 'allowed_to_be_absent' is present and .true., then 0D0 is returned.      !
    ! Note that the energy term name only has to *cotain* name for a match.    !
    ! This is needed to match 'hashed' energy terms.                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2015.06.                                    !
    !==========================================================================!

    implicit none

    character(len=*), intent(in)  :: name
    logical, intent(in), optional :: allowed_to_be_absent

    integer :: found, term, found_at
    logical :: allowed_to_be_absent_local

    ! -------------------------------------------------------------------------
    if(present(allowed_to_be_absent)) then
       allowed_to_be_absent_local = allowed_to_be_absent
    else
       allowed_to_be_absent_local = .false.
    end if

    found = 0
    do term = 1, n_energy_terms
       if(index(trim(energy_term_names(term)), trim(name)) /= 0) then
          found = found + 1
          found_at = term
       end if
    end do
    call utils_assert(found > 0 .or. allowed_to_be_absent_local, &
         'internal_eterm_by_name: Expected energy term "'//&
         trim(name)//'" not found')
    call utils_assert(found <= 1, &
         'internal_eterm_by_name: Energy term "'//&
         trim(name)//'" occurs more than once')

    if(found == 0) then
       internal_eterm_by_name = 0D0
    else
       internal_eterm_by_name = energy_term_values(found_at)
    end if

  end function internal_eterm_by_name

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  end subroutine pol_emb_potential_and_energy

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pol_emb_potential_rspc(local_pol_emb_pot, &
       sph_multipoles, grid)
    !==========================================================================!
    ! Calculates, in real space, the potential coming from a set of spherical  !
    ! multipoles. The potential is returned in local_pol_emb_pot. Damping is   !
    ! applied, depending on the settings in sph_multipoles.                    !
    ! Usual slab-based distribution is used, with OMP used for points in the   !
    ! slab. The multipoles are replicated on all ranks.                        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2015.05.                                    !
    !==========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_real_pt
    use geometry, only: POINT
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET, DAMPING_ZERO, &
         multipole_pot_of_spherical_set
    use rundat, only: pub_threads_max
    use timer, only: timer_clock

    implicit none

    ! jd: Arguments
    type(SPHERICAL_MULTIPOLE_SET), intent(in) :: sph_multipoles
    type(GRID_INFO), intent(in)               :: grid
    real(kind=DP), intent(out)                :: local_pol_emb_pot(&
         grid%ld1, grid%ld2, grid%max_slabs12)

    ! jd: Local variables
    real(kind=DP) :: r_point(3)
    integer       :: i1, i2, islab12, ipt
    character(len=*), parameter :: myself='pol_emb_potential_rspc'

    ! ---------------------------------------------------------------------------

    call timer_clock(myself,1)

    local_pol_emb_pot(:,:,:) = 0D0 ! jd: Take care of pt-ld padding, and zero

    if(sph_multipoles%potential_damping == DAMPING_ZERO) goto 999

    ! jd: Generate the multipole potential in real space, for my slabs
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(none) SCHEDULE(dynamic) &
!$OMP PRIVATE(ipt, i1, i2, islab12, r_point) &
!$OMP SHARED(grid, pub_threads_max, local_pol_emb_pot, sph_multipoles)
    do ipt=1, grid%num_my_slabs12 * grid%n2 * grid%n1
       i1 = modulo(ipt-1,grid%n1) + 1
       i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
       islab12 = (ipt-(i2-1)*grid%n1-i1) / &
            (grid%n1*grid%n2) + 1

       call cell_grid_real_pt(r_point,i1,i2,islab12,grid)

       local_pol_emb_pot(i1,i2,islab12) = &
           multipole_pot_of_spherical_set(sph_multipoles, &
           POINT(r_point(1),r_point(2),r_point(3)))
    end do
!$OMP END PARALLEL DO

999 call timer_clock(myself,2)

  end subroutine pol_emb_potential_rspc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function pol_emb_core_energy_rspc(sph_multipoles, elements, par)
    !==========================================================================!
    ! Calculates, in real space, the energy of interaction of a set of sph.    !
    ! multipoles with ONETEP ionic cores. Classical atoms do not count as cores!
    ! here. Damping is applied, depending on the settings in sph_multipoles.   !
    ! Each MPI rank calculates the interaction of the spherical multipoles with!
    ! its cores (as defined by 'par'), the results are then reduced.           !
    ! OMP is not used. The multipoles are replicated on all ranks.             !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2015.05.                                    !
    !==========================================================================!

    use comms, only: pub_my_proc_id, comms_reduce
    use geometry, only: POINT
    use ion, only: ELEMENT
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET, &
         DAMPING_THOLE, DAMPING_ZERO, DAMPING_FROM_POTENTIAL, &
         DAMPING_COULOMBIC_MASKED, DAMPING_COULOMBIC_SMEARED, &
         multipole_pot_of_spherical_set
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_pol_emb_vacuum_qmstar, pub_pol_emb_dma_min_l, &
         pub_pol_emb_vacuum_dma_min_l
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    type(PARAL_INFO), intent(in)              :: par
    type(SPHERICAL_MULTIPOLE_SET), intent(in) :: sph_multipoles
    type(ELEMENT), intent(in)                 :: elements(par%nat)

    ! jd: Local variables
    integer       :: i_atom, global_i_atom
    integer       :: n_my_atoms
    real(kind=DP) :: energy
    real(kind=DP) :: total_mpole_pot
    real(kind=DP) :: alpha_core
    type(POINT)   :: r_core
    integer       :: min_l
    character(len=*), parameter :: myself = 'pol_emb_core_energy_rspc'

    ! ---------------------------------------------------------------------------

    call timer_clock(myself,1)

    energy = 0D0

    ! jd: Calculate the energy even when DAMPING_ZERO is set, for informational
    !     purposes. The energy will be zeroed elsewhere, after the report
    !     comparing it with the external energy has been shown. Only if
    !     POTENTIAL_DAMPING_ZERO is set, zero it.
    select case(sph_multipoles%energy_damping)
    case(DAMPING_ZERO)
    case(DAMPING_FROM_POTENTIAL)
       select case(sph_multipoles%potential_damping)
       case(DAMPING_ZERO)
          goto 999 ! return 0D0
       case(DAMPING_COULOMBIC_MASKED)
       case(DAMPING_COULOMBIC_SMEARED)
       case(DAMPING_THOLE)
          continue ! proceed as usual
       case default
          call utils_abort(myself//': Unrecognized potential damping setting')
       end select
    case default
       call utils_abort(myself//': Unrecognized energy damping setting')
    end select

    ! jd: Handle corner case where no DMA set includes charges (min_l > 0).
    !     This is very exotic (expanding the molecule in terms of dipoles only),
    !     but technically possible. In that case monopoles are not written to
    !     the GDMA file and MM does not see the cores at all.
    if(pub_pol_emb_vacuum_qmstar) then
       min_l = min(pub_pol_emb_dma_min_l, pub_pol_emb_vacuum_dma_min_l)
    else
       min_l = pub_pol_emb_dma_min_l
    end if
    if(min_l > 0) goto 999 ! return 0D0

    n_my_atoms = par%num_atoms_on_proc(pub_my_proc_id)

    ! jd: Multipoles are replicated, parallelisation is done over QM ion cores
    do i_atom = 1, n_my_atoms
       global_i_atom = par%first_atom_on_proc(pub_my_proc_id) + i_atom -1
       r_core = elements(par%orig_atom(global_i_atom))%centre
       alpha_core = elements(par%orig_atom(global_i_atom))%thole_polarisability
       total_mpole_pot = &
            - multipole_pot_of_spherical_set(sph_multipoles, r_core, alpha_core)
       energy = energy + &
            elements(par%orig_atom(global_i_atom))%ion_charge * total_mpole_pot
    end do

    call comms_reduce('SUM', energy)

999 pol_emb_core_energy_rspc = energy

    call timer_clock(myself,2)

  end function pol_emb_core_energy_rspc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pol_emb_read_from_file(energy_term_names, &    ! out
       energy_term_values, n_energy_terms, sph_multipoles, n_sets, &      ! out
       filename, silent)                                                  ! in
    !==========================================================================!
    ! Reads a definition of a polarisable embedding (a set of sets of mdqpoles,!
    ! each mdqpole located at a site) from a file. Each set of mdqpoles is     !
    ! converted a spherical multipole set and returned in this format, the     !
    ! result being an allocatable array of type(SPHERICAL_MULTIPOLE_SET), with !
    ! the array length of n_sets. It is the responsibility of the caller to    !
    ! deallocate this structure, calling multipole_free_spherical_set() on each!
    ! of the array elements first.
    !                                                                          !
    ! The data is replicated -- everyone gets a copy of the entire set of mul- !
    ! tipoles. The assumption is that slab-based parallelism will be used and  !
    ! the embedding will not need to be distributed.                           !
    !                                                                          !
    ! Also, a number of energy terms are potentially read. These will be termed!
    ! 'external' energy terms. These should go in the first lines of the file  !
    ! and are identified by the presence of a ':', like that:                  !
    ! valence: 0.00034                                                         !
    ! permanent multipole: -0.0074                                             !
    ! The names and values (in au) will be stored in energy_term_names(:) and  !
    ! energy_term_values(:), respectively. Additionally, two energy terms,     !
    ! calculated internally (and hence termed 'internal') will be  added for   !
    ! every multipole set, as described in the documentation to                !
    ! pol_emb_potential_and_energy(). The internal energy terms                !
    ! are the first (3*n_sets) ones in the returned arrays. The arrays are     !
    ! allocated to the dimension of n_energy_terms here, where n_energy_terms  !
    ! includes the (3*n_sets) extra energy terms, and is thus at least         !
    ! (3*n_sets). The caller will be responsible for deallocating these arrays.!
    !                                                                          !
    ! The energy term values, and names are replicated, this becomes handy in  !
    ! hamiltonian_mod.                                                         !
    !                                                                          !
    ! We follow a convention where external energy terms containing '#' in the !
    ! name are later *NOT* included in the Hamiltonian. Whether the internal   !
    ! described in pol_emb_potential_and_energy() are added or not is          !
    ! controlled by the header of each multipole set. The headers are          !
    ! constructed thus:                                                        !
    ! multipole <n> <name> <directive_1> <directive_2> <directive_3> ...       !
    ! where <n> is the number of multipoles in the set, <name> is a name by    !
    ! which the multipole set will be identified and <directive_n> are optional!
    ! directives for controlling the behaviour of the multipole set. The       !
    ! following directives are understood so far:                              !
    ! - potential_coulombic: Generate a Coulombic potential for the mpole set  !
    !                        and include it in the Hamiltonian.                !
    ! - potential_thole: Generate a Thole-damped potential for the mpole set   !
    !                    and include it in the Hamiltonian.                    !
    !                    (not implemented at this stage, asserts).             !
    ! - potential_zero: Generate no potential for the mpole set, essentially   !
    !                   ignoring it.                                           !
    ! - energy_from_potential: Include the internal energy contribution from   !
    !                          the multipole set interacting with ONETEP cores !
    !                          and electrons, using the same damping strategy  !
    !                          that is used for the potential, thus ensuring   !
    !                          consistency.                                    !
    ! - energy_zero: Do not include the internal energy contribution, with the !
    !                assumption that one of the external energy terms does.    !
    ! Unrecognised directives cause asserts.                                   !
    ! An example multipole header could look like:                             !
    ! multipoles 5 permanent potential_coulombic energy_from_potential         !
    ! The 5 multipoles themselves would then follow the format:                !
    ! <species> <x> <y> <z> <q> <dx> <dy> <dz> <Qxx> <Qxy> <Qxz> <Qyx> \       !
    ! <Qyy> <Qyz> <Qzx> <Qzy> <Qzz>                                            !
    ! where <species> is ignored; <x>, <y>, <z> define the site's position,    !
    ! <q> is the charge (electrons negative conventin), <dx>, <dy> <dz> are    !
    ! the Cartesian components of the dipole, and <Qxx>..<Qzz> are the         !
    ! Cartesian primitive components of the quadrupole. Be careful with the    !
    ! conventions here.                                                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   energy_term_names (out): Energy term names will be returned here.      !
    !   energy_term_values (out): Energy term values will be returned here.    !
    !   sph_multipoles (out): External multipoles, converted to spherical      !
    !                         representation will be returned here.            !
    !   n_sets (out): The number of sets of multipoles (currently 1 or 2) will !
    !                 be returned here.                                        !
    !   filename (in): Input filename from which the multipoles and energy     !
    !                  terms are to be read.                                   !
    !   silent (in, opt): If provided and .true., no user feedback will be     !
    !                     produced. This is useful when this subroutine is     !
    !                     called from the NGWF gradient.                       !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   All the multipole sets and energy terms are replicated.                !
    !   The caller is responsible for deallocating energy_term_names later on. !
    !   The caller is responsible for deallocating energy_term_values later on.!
    !   The caller is responsible for deinitialising, and then deallocating    !
    !   sph_multipoles when no longer needed.                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2015.05-07.                                 !
    ! Updated by Jacek Dziedzic for NGWF gradient, 2016.11.                    !
    !==========================================================================!

    use comms, only: pub_on_root, pub_root_proc_id, comms_bcast
    use constants, only: garbage_real
    use multipole_ops, only: CART_PRIM_MDQPOLE, CART_PRIM_MDQPOLE_NDATA, &
         SPHERICAL_MULTIPOLE_SET, multipole_mdq_cart_to_spherical_set
    use rundat, only: pub_output_detail
    use timer, only: timer_clock
    use utils, only: utils_flushed_string_output, utils_unit, utils_abort, &
         utils_int_to_str, utils_alloc_check, utils_dealloc_check, &
         utils_assert, utils_int_to_str, utils_str_to_int, utils_nth_word_of, &
         utils_banner, utils_int_to_str, utils_point_to_str

    implicit none

    ! jd: Arguments
    character(len=27), intent(out), allocatable       :: energy_term_names(:)
    real(kind=DP), intent(out), allocatable           :: energy_term_values(:)
    type(SPHERICAL_MULTIPOLE_SET), intent(out), allocatable :: sph_multipoles(:)
    integer, intent(out)                              :: n_energy_terms
    integer, intent(out)                              :: n_sets
    character(len=*), intent(in)                      :: filename
    logical, intent(in), optional                     :: silent

    ! jd: Local variables
    integer            :: n_sites(MAX_N_MULTIPOLE_SETS)
    character(len=80)  :: set_names(MAX_N_MULTIPOLE_SETS)
    character(len=39)  :: filename_truncated
    integer            :: mpole_file_unit
    character(len=512) :: dummy_string
    character(len=78)  :: outstr
    character(len=80)  :: ignorestr
    logical            :: local_silent
    integer            :: end_of_name
    integer            :: ierr
    integer            :: i
    integer            :: offs
    integer            :: n_energy_terms_to_ignore
    integer            :: n_sites_total
    integer            :: mpole_set
    logical            :: set_is_qmstar, set_is_repulsive
    integer            :: potential_damping, energy_damping
    real(kind=DP)      :: Qxx, Qxy, Qxz, Qyx, Qyy, Qyz, Qzx, Qzy, Qzz
    real(kind=DP)      :: alpha ! polarisability
    real(kind=DP)      :: vdw_rmin
    real(kind=DP)      :: vdw_eps
    real(kind=DP), allocatable :: bcast_buf(:)
    ! jd: Cartesian {mono-di-quadru-}poles read from file for current set.
    !     Indexed by site.
    type(CART_PRIM_MDQPOLE), allocatable :: mdqpoles_cart(:)

    logical, save      :: mpole_file_continue = .false.
    character(len=*), parameter :: myself = 'pol_emb_read_from_file'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    local_silent = .false.
    if(present(silent)) then
       local_silent = silent
    end if

    if(pub_output_detail >= VERBOSE .and. .not. local_silent) then
       filename_truncated = filename
       call utils_flushed_string_output(CRLF // &
            utils_banner('~', delimiter='/') // CRLF // &
            '| Polarisable embedding potential from '//filename_truncated//&
            ' |'//CRLF)
    end if

    if(pub_on_root) then
       ! jd: Open the multipole file
       mpole_file_unit = utils_unit()
       open(unit=mpole_file_unit, file=filename, status="old", err=100)
       call internal_seek_to_record()

       ! jd: Count energy terms and, multipole sets and sites in each set
       !     so that we know how big the arrays need to be
       n_energy_terms = 0
       n_energy_terms_to_ignore = 0
       n_sets = 0
       n_sites(:) = 0
       set_names(:) = "????????"
       n_sites_total = 0
       do while(.true.)
          read(mpole_file_unit,'(a)',iostat=ierr) dummy_string
          if(ierr > 0) goto 200 ! read error
          if(ierr < 0) exit     ! EOF
          if(index(dummy_string,':') /= 0) then
             ! jd: Energy terms contain ':' and they come first
             n_energy_terms = n_energy_terms + 1
             call utils_assert(n_sites_total == 0, &
                  myself//': Syntax error in "'//trim(filename)//&
                  '". All energy terms must precede multipole sites. Also, if &
                  &using blank lines to delineate records, make sure none are &
                  &missing')
             if(dummy_string(1:1) == '#') then
                n_energy_terms_to_ignore = n_energy_terms_to_ignore + 1
             end if
          else if(index(dummy_string,'multipoles') /= 0) then
             ! jd: Multipole headers contain 'multipoles'
             n_sets = n_sets + 1
             n_sites(n_sets) = &
                  utils_str_to_int(utils_nth_word_of(dummy_string,2))
             set_names(n_sets) = trim(utils_nth_word_of(dummy_string,3))

             if(pub_output_detail >= VERBOSE .and. .not. local_silent) then
                outstr = &
                     adjustl('| Multipole set "'//&
                     trim(set_names(n_sets))//&
                     '": '//trim(utils_int_to_str(n_sites(n_sets)))//&
                     ' sites.')
                write(stdout,'(a80)') outstr//' |'
             end if
          else if(len(trim(dummy_string)) == 0) then
             ! jd: Blank line indicates end of record
             mpole_file_continue = .true.
             exit
          else
             ! jd: All other lines are multipole definitions
             n_sites_total = n_sites_total + 1
          end if
       end do ! lines in mpole input file
       call utils_assert(n_sites_total == sum(n_sites(:)), &
            'Mismatch in the number of multipole sites in "'//trim(filename)//&
            '": Expecting '//trim(utils_int_to_str(sum(n_sites)))//&
            ' sites, got '//trim(utils_int_to_str(n_sites_total))//&
            ' sites instead.')
    end if ! on root
    if(pub_output_detail >= VERBOSE .and. pub_on_root .and. &
         .not. local_silent) then
       if(n_energy_terms_to_ignore == 0) then
          ignorestr = ' (all included).'
       else
          ignorestr = ' (out of which '//&
               trim(utils_int_to_str(n_energy_terms_to_ignore))//' #ignored).'
       end if
       outstr = adjustl('| All in all '//&
            trim(utils_int_to_str(n_sites_total))//' sites, and '//&
            trim(utils_int_to_str(n_energy_terms))//' external energy terms'//&
            ignorestr)
       write(stdout,'(a80)') outstr//' |'
       write(stdout,'(a80)') '| Energy term                           &
            &   Energy      Source   Included?      |'
    end if

    ! jd: Add three extra energy terms per set that are calculated internally
    !     in ONETEP (interaction of the set with cores and with electrons,
    !     elecs <-> MM repulsive potential)
    if(pub_on_root) n_energy_terms = n_energy_terms + 3*n_sets

    call comms_bcast(pub_root_proc_id, n_sets)
    call comms_bcast(pub_root_proc_id, n_sites)
    call comms_bcast(pub_root_proc_id, n_sites_total)
    call comms_bcast(pub_root_proc_id, n_energy_terms)

    ! jd: Allocate spherical multipole sets
    allocate(sph_multipoles(n_sets),stat=ierr)
    call utils_alloc_check(myself,'sph_multipoles',ierr)

    ! jd: Allocate datastructure
    allocate(energy_term_names(n_energy_terms),stat=ierr)
    call utils_alloc_check(myself,'energy_term_names',ierr)
    allocate(energy_term_values(n_energy_terms),stat=ierr)
    call utils_alloc_check(myself,'energy_term_values',ierr)

    ! jd: Set dummy values for energy terms on all procs, they will be
    !     later overwritten.
    energy_term_names(:) = ""
    energy_term_values(:) = 0D0

    if(pub_on_root) then

       ! jd: Rewind and read the records
       rewind(mpole_file_unit,err=300)
       call internal_seek_to_record()

       ! ---------------------------------------
       ! --- Energy terms external to ONETEP ---
       ! ---------------------------------------
       do i = 3*n_sets+1, n_energy_terms
          read(mpole_file_unit,'(a)',iostat=ierr) dummy_string
          end_of_name = index(dummy_string,':')
          read(dummy_string(end_of_name+1:),*) energy_term_values(i)
          energy_term_names(i) = dummy_string(1:end_of_name-1)
          if(pub_output_detail >= VERBOSE .and. .not. local_silent) then
             write(outstr,'(a,a,a,f18.9,a)') '| - ', &
                  energy_term_names(i), ' ', energy_term_values(i),' Ha  '//&
                  'TINKER   '//merge('    NO','   YES',dummy_string(1:1)=='#')
             write(stdout,'(a80)') outstr//' |'
          end if
       end do

       ! -----------------------------------------
       ! --- Energy terms calculated in ONETEP ---
       ! -----------------------------------------
       do mpole_set = 1, n_sets
          set_is_qmstar = (index(set_names(mpole_set),'qmstar')/=0)
          if(set_is_qmstar) then
             energy_term_names(3*mpole_set-1) = &
                  'QM* elec <-> MM '//trim(set_names(mpole_set))
          else
             energy_term_names(3*mpole_set-1) = &
                  'QM elec <-> MM '//trim(set_names(mpole_set))
          end if
          energy_term_names(3*mpole_set-2) = &
               'QM elec <-> rep MM '//trim(set_names(mpole_set))
          energy_term_names(3*mpole_set) = &
               'QM core <-> MM '//trim(set_names(mpole_set))
          ! jd: Initialise values to garbage, so that comms_bcast doesn't touch
          !     uninitialised data. These terms will be set in
          !     pol_emb_potential_and_energy().
          energy_term_values(3*mpole_set-2) = garbage_real
          energy_term_values(3*mpole_set-1) = garbage_real
          energy_term_values(3*mpole_set) = garbage_real
       end do
    end if

    ! jd: Broadcast energy term names (needed, as in hamiltonian we're going
    !     to need to exclude energy terms with hashes in the name).
    !     The values are also needed everywhere. Set names too, for keeping
    !     track of what's permanent.
    call comms_bcast(pub_root_proc_id, energy_term_names)
    call comms_bcast(pub_root_proc_id, energy_term_values)
    call comms_bcast(pub_root_proc_id, set_names)

    ! --------------------------------------------------
    ! --- Loop over multipole sets and read each one ---
    ! --------------------------------------------------
    do mpole_set = 1, n_sets

       allocate(mdqpoles_cart(n_sites(mpole_set)),stat=ierr)
       call utils_alloc_check(myself,'mdqpoles_cart',ierr)

       if(pub_on_root) then

          ! First the header for a set
          read(mpole_file_unit,'(a)',iostat=ierr) dummy_string

          ! jd: Parse 'dummy_string', set potential_damping, energy_damping
          call internal_parse_set_properties
       end if

       call comms_bcast(pub_root_proc_id, potential_damping)
       call comms_bcast(pub_root_proc_id, energy_damping)

       if(pub_on_root) then

          set_is_repulsive = (index(set_names(mpole_set),'rep')/=0)

          ! --- Multipole sets ---
          do i = 1, n_sites(mpole_set)

             read(mpole_file_unit,*,iostat=ierr) &
                  mdqpoles_cart(i)%species, &
                  mdqpoles_cart(i)%centre%x, &
                  mdqpoles_cart(i)%centre%y, &
                  mdqpoles_cart(i)%centre%z, &
                  mdqpoles_cart(i)%charge, &
                  mdqpoles_cart(i)%dipole(1), &
                  mdqpoles_cart(i)%dipole(2), &
                  mdqpoles_cart(i)%dipole(3), &
                  Qxx, Qxy, Qxz, Qyx, Qyy, Qyz, Qzx, Qzy, Qzz, alpha, &
                  vdw_rmin, vdw_eps
             call utils_assert(ierr == 0,myself//': Error reading '//&
                  trim(filename)//': the file is likely misformatted.',ierr)
             call utils_assert(Qxy == Qyx, myself//': Qxy != Qyx for site '//&
                  trim(utils_int_to_str(i))//' at '//&
                  utils_point_to_str(mdqpoles_cart(i)%centre,'f12.6'))
             call utils_assert(Qxz == Qzx, myself//': Qxz != Qzx for site '//&
                  trim(utils_int_to_str(i))//' at '//&
                  utils_point_to_str(mdqpoles_cart(i)%centre,'f12.6'))
             call utils_assert(Qyz == Qzy, myself//': Qyz != Qzy for site '//&
                  trim(utils_int_to_str(i))//' at '//&
                  utils_point_to_str(mdqpoles_cart(i)%centre,'f12.6'))
             call utils_assert(alpha >= 0D0, &
                  myself//': Negative polarisability for site '//&
                  trim(utils_int_to_str(i))//' at '//&
                  utils_point_to_str(mdqpoles_cart(i)%centre,'f12.6'))

             ! @ The input monopoles are negated due to the charge sign
             !   convention being different in ONETEP. Dipoles are not, as
             !   we have (-1)^l in the potential for now. Quadrupoles are
             !   negated for the same reason.
             mdqpoles_cart(i)%charge = -mdqpoles_cart(i)%charge

             mdqpoles_cart(i)%quadrupole(1) = -Qxx
             mdqpoles_cart(i)%quadrupole(2) = -Qxy
             mdqpoles_cart(i)%quadrupole(3) = -Qxz
             mdqpoles_cart(i)%quadrupole(4) = -Qyy
             mdqpoles_cart(i)%quadrupole(5) = -Qyz
             mdqpoles_cart(i)%quadrupole(6) = -Qzz
             mdqpoles_cart(i)%polarisability = alpha
             mdqpoles_cart(i)%masked = (index(mdqpoles_cart(i)%species,'#') /= 0)
             mdqpoles_cart(i)%userdata(:) = garbage_real
             if(set_is_repulsive .and. .not. mdqpoles_cart(i)%masked) then
                mdqpoles_cart(i)%userdata(1) = vdw_rmin ! jd: for old scheme
                mdqpoles_cart(i)%userdata(2) = vdw_eps  ! jd: for old scheme
             end if

             ! jd: New MM rep scheme. Store MM rep params in the multipole struct.
             if(set_is_repulsive .and. .not. mdqpoles_cart(i)%masked) then
                mdqpoles_cart(i)%userdata(1) = &
                     pol_emb_mm_species_to_rep_param(1,mdqpoles_cart(i)%species)
                mdqpoles_cart(i)%userdata(2) = &
                     pol_emb_mm_species_to_rep_param(2,mdqpoles_cart(i)%species)
             else
                mdqpoles_cart(i)%userdata(1) = garbage_real
                mdqpoles_cart(i)%userdata(2) = garbage_real
             endif

          end do

       end if ! on root

       ! jd: Broadcast the site data to everyone
       if(n_sites(mpole_set) > 0) then
          allocate(bcast_buf(n_sites(mpole_set) * CART_PRIM_MDQPOLE_NDATA), &
               stat=ierr)
          call utils_alloc_check(myself,'bcast_buf',ierr)
          if(pub_on_root) then
             ! jd: On root pack data and send through a bcast
             do i = 1, n_sites(mpole_set)
                offs = (i-1)*CART_PRIM_MDQPOLE_NDATA
                bcast_buf(offs+1) = mdqpoles_cart(i)%centre%x
                bcast_buf(offs+2) = mdqpoles_cart(i)%centre%y
                bcast_buf(offs+3) = mdqpoles_cart(i)%centre%z
                bcast_buf(offs+4) = mdqpoles_cart(i)%charge
                bcast_buf(offs+5) = mdqpoles_cart(i)%dipole(1)
                bcast_buf(offs+6) = mdqpoles_cart(i)%dipole(2)
                bcast_buf(offs+7) = mdqpoles_cart(i)%dipole(3)
                bcast_buf(offs+8) = mdqpoles_cart(i)%quadrupole(1)
                bcast_buf(offs+9) = mdqpoles_cart(i)%quadrupole(2)
                bcast_buf(offs+10) = mdqpoles_cart(i)%quadrupole(3)
                bcast_buf(offs+11) = mdqpoles_cart(i)%quadrupole(4)
                bcast_buf(offs+12) = mdqpoles_cart(i)%quadrupole(5)
                bcast_buf(offs+13) = mdqpoles_cart(i)%quadrupole(6)
                bcast_buf(offs+14) = mdqpoles_cart(i)%polarisability
                bcast_buf(offs+15:offs+17) = mdqpoles_cart(i)%userdata(1:3)
                if(mdqpoles_cart(i)%masked) then
                   bcast_buf(offs+18) = 1.0_DP
                else
                   bcast_buf(offs+18) = 0.0_DP
                end if
             end do
             call comms_bcast(pub_root_proc_id, bcast_buf)
          else
             ! jd: On slaves recv through a bcast and unpack data
             call comms_bcast(pub_root_proc_id, bcast_buf)
             do i = 1, n_sites(mpole_set)
                offs = (i-1)*CART_PRIM_MDQPOLE_NDATA
                mdqpoles_cart(i)%centre%x = bcast_buf(offs+1)
                mdqpoles_cart(i)%centre%y = bcast_buf(offs+2)
                mdqpoles_cart(i)%centre%z = bcast_buf(offs+3)
                mdqpoles_cart(i)%charge = bcast_buf(offs+4)
                mdqpoles_cart(i)%dipole(1) = bcast_buf(offs+5)
                mdqpoles_cart(i)%dipole(2) = bcast_buf(offs+6)
                mdqpoles_cart(i)%dipole(3) = bcast_buf(offs+7)
                mdqpoles_cart(i)%quadrupole(1) = bcast_buf(offs+8)
                mdqpoles_cart(i)%quadrupole(2) = bcast_buf(offs+9)
                mdqpoles_cart(i)%quadrupole(3) = bcast_buf(offs+10)
                mdqpoles_cart(i)%quadrupole(4) = bcast_buf(offs+11)
                mdqpoles_cart(i)%quadrupole(5) = bcast_buf(offs+12)
                mdqpoles_cart(i)%quadrupole(6) = bcast_buf(offs+13)
                mdqpoles_cart(i)%polarisability = bcast_buf(offs+14)
                mdqpoles_cart(i)%userdata(1:3) = bcast_buf(offs+15:offs+17)
                if(bcast_buf(offs+18) == 1.0_DP) then
                   mdqpoles_cart(i)%masked = .true.
                else
                   mdqpoles_cart(i)%masked = .false.
                end if

             end do
          end if
          deallocate(bcast_buf, stat=ierr)
          call utils_dealloc_check(myself,'bcast_buf',ierr)
       end if

       ! jd: Convert Cartesian multipoles to sphericals
       call multipole_mdq_cart_to_spherical_set(sph_multipoles(mpole_set), &
            mdqpoles_cart, n_sites(mpole_set), set_names(mpole_set), &
            potential_damping, energy_damping)

       ! jd: Cartesian multipoles no longer needed, destroy them
       deallocate(mdqpoles_cart,stat=ierr)
       call utils_dealloc_check(myself,'mdqpoles_cart',ierr)

    end do ! over mpole sets

    if(pub_on_root) then
       close(unit=mpole_file_unit, err=400)
    end if

    call timer_clock(myself,2)

    return

100 call utils_abort('Error opening file: '//trim(filename))
200 call utils_abort('Error reading file: '//trim(filename))
300 call utils_abort('Error rewinding file: '//trim(filename))
400 call utils_abort('Error closing file: '//trim(filename))

contains

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine internal_seek_to_record

    use utils, only: utils_abort

    implicit none

    integer :: records_skipped

    ! -------------------------------------------------------------------------

    if(.not. mpole_file_continue) return
    if(pol_emb_scf_iter == 1) return

    records_skipped = 0
    do while(.true.)
       read(mpole_file_unit,'(a)',iostat=ierr) dummy_string
       if(ierr > 0) goto 200 ! read error
       if(ierr < 0) goto 300 ! EOF
       if(len(trim(dummy_string)) == 0) then
          records_skipped = records_skipped + 1
          if(records_skipped == pol_emb_scf_iter-1) then
             exit
          end if
       end if
    end do

    return

200 call utils_abort('Error reading file: "'//trim(filename)//&
         '" when seeking to record '//trim(utils_int_to_str(pol_emb_scf_iter)))
300 call utils_abort('Unexpected EOF in file: "'//trim(filename)//&
         '" when seeking to record '//trim(utils_int_to_str(pol_emb_scf_iter))&
         //'. Managed to fast-forward over '//&
         trim(utils_int_to_str(records_skipped))//' records.')

  end subroutine internal_seek_to_record

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine internal_parse_set_properties

    use constants, only: MAX_WORD_LENGTH
    use multipole_ops, only: DAMPING_THOLE, DAMPING_FROM_POTENTIAL, &
         DAMPING_COULOMBIC_MASKED, DAMPING_COULOMBIC_SMEARED, DAMPING_ZERO

    implicit none

    integer :: n
    character(len=MAX_WORD_LENGTH) :: word

    ! -------------------------------------------------------------------------
    do n = 4, 999
       word = utils_nth_word_of(dummy_string,n)

       select case(trim(word))
          case ('potential_coulombic_masked')
             potential_damping = DAMPING_COULOMBIC_MASKED
          case ('potential_coulombic_smeared')
             potential_damping = DAMPING_COULOMBIC_SMEARED
          case ('potential_zero')
             potential_damping = DAMPING_ZERO
          case ('potential_thole_damped')
             potential_damping = DAMPING_THOLE
          case ('energy_zero')
             energy_damping = DAMPING_ZERO
          case ('energy_from_potential')
             energy_damping = DAMPING_FROM_POTENTIAL
          case ('') ! ran out of keywords
             exit
          case default
             call utils_abort('Unrecognised keyword "'//trim(word)//'" in "'//&
                  trim(filename)//'"')
       end select

    end do

  end subroutine internal_parse_set_properties

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  end subroutine pol_emb_read_from_file

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pol_emb_qmstar_energy_and_lnv_gradient(elec_energy, polemb, & ! out
       driving_dcoeffs_ht, mm_multipoles, qm_driving_multipoles, ngwf_basis, & ! in
       mdl, denskern, overlap, s_atoms_nl, scal_factor, ind_factor, & ! in
       qm_scf_dma_polarisation_multipoles, qm_vacuum_multipoles, &    ! OPT
       vacuum_polemb, vacuum_aux_polemb, vacuum_dcoeffs_ht, vacuum_aux_dcoeffs_ht) ! OPT
    !==========================================================================!
    ! Calculates the P matrix and associated energy term Tr[KP]. The P matrix  !
    ! contains contributions to the Hamiltonian from the polarisable embedding !
    ! terms that use the QM* approach.                                         !
    !                                                                          !
    ! If charge scaling is in effect, there is an extra term due to the        !
    ! dependence of total number of electrons on the density kernel. This is   !
    ! slightly tricky and cannot be evaluated in the same loop -- the correc-  !
    ! tion needs a term that needs to be summed across all contributions first.!
    ! Only then a second loop weights the correction with S matrix elements and!
    ! includes it. Note that the first loop deals with all combinations of     !
    ! QM {charges, dipoles, quadrupoles} <-> MM {charges, dipoles, quadupoles},!
    ! while the charge correction only involves QM charges (as only they are   !
    ! scaled), interacting with MM {charges, dipoles, quadupoles}.             !
    ! This has been tested with kernel FDs to ensure gradients are consistent. !
    ! Any potential breakage in charge scaling correction is sneaky and is     !
    ! nigh impossible to detect in standard LNV convergence loops -- i.e. the  !
    ! commutator does go down nicely even if the charge correction is omitted. !
    ! Only an LNV FD test uncovers problems there. Now there are none.         !
    !--------------------------------------------------------------------------!
    ! elec_energy (out): The calculated energy term is returned here.          !
    ! polemb (inout): A SPAM3 where the P matrix will be returned.             !
    !                 Must be sparse_created on entry and should be zeroed by  !
    !                 the caller, since we add to it.                          !
    ! dcoeffs_ht (in): The D coefficients from which the P matrix will be      !
    !                     calculated.                                          !
    ! mm_multipoles (in): Spherical multipoles describing the MM side.         !
    ! qm_multipoles (in): Spherical multipoles (two sets: bessel runs)         !
    !                     describing the QM side.                              !
    ! ngwf_basis, elements, denskern, overlap (in): The usual.                 !
    ! s_atoms_nl (in): Neighbour list for S sparsity.                          !
    ! scal_factor (in): Factor with which the elements of the P matrix are     !
    !                   pre-multiplied -- depends on pub_dma_multipole_scaling !
    !                   and pub_pol_emb_perm_scaling.                          !
    ! ind_factor (in): Factor accounting for the different constant in energy  !
    !                  expressions for permanent (1) interactions and induced  !
    !                  interactions (1/2).                                     !
    ! @updatedoc: args
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2016.09, based on code extracted from       !
    ! pol_emb_potential_and_energy().                                          !
    ! Modified to remove pub_par by Joseph Prentice, September 2018            !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: REP_SWEX_POL_EMB_DMA_1, REP_SWEX_POL_EMB_AUX_DMA_1
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE
    use model_type, only: MODEL
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET, DAMPING_ZERO
    use neighbour_list, only: NL_NEIGHBOUR_LIST
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_dma_bessel_averaging, pub_num_spins, PUB_1K, &
         pub_pol_emb_dma_min_l, pub_pol_emb_dma_max_l, pub_debug, &
         pub_pol_emb_vacuum_dma_min_l, pub_pol_emb_vacuum_dma_max_l, &
         pub_pol_emb_vacuum_qmstar
    use sparse, only: SPAM3, sparse_trace, sparse_create, sparse_axpy, &
         sparse_destroy, sparse_show_matrix, sparse_get_par
    use sparse_array, only: SPAM3_ARRAY
    use sw_expansion_type, only: SW_EX
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)                :: elec_energy ! jd: Tr[KP] @notreally
    type(SPAM3), intent(inout)                :: polemb
    type(HT_HASH_TABLE), intent(in), target   :: driving_dcoeffs_ht(:)
    type(SPHERICAL_MULTIPOLE_SET), intent(in) :: mm_multipoles
    type(SPHERICAL_MULTIPOLE_SET), intent(in) :: qm_driving_multipoles(:)
    type(FUNC_BASIS), intent(in)              :: ngwf_basis
    type(MODEL), intent(in)                   :: mdl
    type(SPAM3_ARRAY), intent(in)             :: denskern
    type(SPAM3), intent(in)                   :: overlap
    type(NL_NEIGHBOUR_LIST), intent(in)       :: s_atoms_nl
    real(kind=DP), intent(in)                 :: scal_factor
    real(kind=DP), intent(in)                 :: ind_factor
    type(SPHERICAL_MULTIPOLE_SET), intent(in), optional :: qm_scf_dma_polarisation_multipoles(:)
    type(SPHERICAL_MULTIPOLE_SET), intent(in), optional :: qm_vacuum_multipoles(:)
    type(SPAM3), intent(inout), optional                :: vacuum_polemb
    type(SPAM3), intent(inout), optional                :: vacuum_aux_polemb
    type(HT_HASH_TABLE), intent(in), target, optional   :: vacuum_dcoeffs_ht(:)
    type(HT_HASH_TABLE), intent(in), target, optional   :: vacuum_aux_dcoeffs_ht(:)

    ! jd: Local variables
    type(SPAM3)                 :: cur_polemb
    type(SPAM3)                 :: cur_polemb_vac
    type(SPAM3)                 :: cur_polemb_vac_aux
    real(kind=DP)               :: polemb_energy_pol
    real(kind=DP)               :: polemb_energy_vac
    real(kind=DP)               :: polemb_energy_vac_aux
    integer                     :: bessel_run, n_bessel_runs
    integer                     :: is
    integer                     :: one_atom_only
    character(len=*), parameter :: myself = &
         'pol_emb_qmstar_energy_and_lnv_gradient'
    type(PARAL_INFO), pointer   :: par

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    if(pub_dma_bessel_averaging) then
       n_bessel_runs = 2
    else
       n_bessel_runs = 1
    end if

    call utils_assert(pub_pol_emb_vacuum_qmstar .eqv. &
         present(qm_scf_dma_polarisation_multipoles), &
         myself//': optional arg qm_scf_dma_polarisation_multipoles &
         &must be passed iff pub_pol_emb_vacuum_qmstar')
    call utils_assert(pub_pol_emb_vacuum_qmstar .eqv. &
         present(qm_vacuum_multipoles), &
         myself//': optional arg qm_vacuum_multipoles &
         &must be passed iff pub_pol_emb_vacuum_qmstar')
    call utils_assert(pub_pol_emb_vacuum_qmstar .eqv. &
         present(vacuum_polemb), &
         myself//': optional arg vacuum_polemb &
         &must be passed iff pub_pol_emb_vacuum_qmstar')
    call utils_assert(pub_pol_emb_vacuum_qmstar .eqv. &
         present(vacuum_aux_polemb), &
         myself//': optional arg vacuum_aux_polemb &
         &must be passed iff pub_pol_emb_vacuum_qmstar')
    call utils_assert(pub_pol_emb_vacuum_qmstar .eqv. &
         present(vacuum_dcoeffs_ht), &
         myself//': optional arg vacuum_dcoeffs_ht &
         &must be passed iff pub_pol_emb_vacuum_qmstar')
    call utils_assert(pub_pol_emb_vacuum_qmstar .eqv. &
         present(vacuum_aux_dcoeffs_ht), &
         myself//': optional arg vacuum_aux_dcoeffs_ht &
         &must be passed iff pub_pol_emb_vacuum_qmstar')

    call sparse_get_par(par, overlap)

    ! jd: If potential not needed, leave P untouched,
    !     just use zero runs to elide the AVG loop entirely.
    if(mm_multipoles%potential_damping == DAMPING_ZERO) n_bessel_runs = 0

    ! jd: Similarly of MM set is effectively empty (zero non-ignorable
    !     multipoles or all multipoles are exactly zero).
    if(maxval(mm_multipoles%max_l_mask(:)) == -1) n_bessel_runs = 0

    ! *************************************************************************
    !                ---***--- Calculate the P matrix ---***---
    ! *************************************************************************

    ! --------------------------------------------------------------------------
    ! jd: Loop over Bessel averaging runs                                   AVG
    ! --------------------------------------------------------------------------
    do bessel_run = 1, n_bessel_runs

       call sparse_create(cur_polemb, polemb)

       if(pub_pol_emb_vacuum_qmstar) then

!          write(*,*) '@ -------------- THIS IS PMAT-POL ------------------ @'
          ! jd: P matrix contribution from current state
          call pol_emb_calc_pmatrix(cur_polemb, overlap, &
               qm_scf_dma_polarisation_multipoles(bessel_run), mm_multipoles, &
               driving_dcoeffs_ht(bessel_run), ngwf_basis, mdl%elements, &
               s_atoms_nl, scal_factor, &
               pub_pol_emb_dma_min_l, pub_pol_emb_dma_max_l, &
               REP_SWEX_POL_EMB_DMA_1 + bessel_run - 1)

!          write(*,*) '@ -------------- THIS IS PMAT-VAC ------------------ @'
          ! jd: P matrix contribution from vacuum state
          ! @@@uses wrong overlap, charge scaling will be borked
          call sparse_create(cur_polemb_vac, polemb)
          call pol_emb_calc_pmatrix(cur_polemb_vac, &
               overlap, qm_vacuum_multipoles(bessel_run), mm_multipoles, &
               vacuum_dcoeffs_ht(bessel_run), ngwf_basis, mdl%elements, &
               s_atoms_nl, scal_factor, &
               pub_pol_emb_vacuum_dma_min_l, pub_pol_emb_vacuum_dma_max_l, &
               REP_SWEX_POL_EMB_DMA_1 + bessel_run - 1)
          call sparse_axpy(vacuum_polemb, cur_polemb_vac, &
               1.0_DP/real(n_bessel_runs,kind=DP))
          call sparse_destroy(cur_polemb_vac)

!          write(*,*) '@ -------------- THIS IS PMAT-AUX ------------------ @'
          ! jd: P matrix contribution from vacuum state -- auxiliary expansion
          ! @@@uses wrong overlap, charge scaling will be borked
          call sparse_create(cur_polemb_vac_aux, polemb)
          call pol_emb_calc_pmatrix(cur_polemb_vac_aux, &
               overlap, qm_scf_dma_polarisation_multipoles(bessel_run), & ! <-- pol multipoles only for quality setting
               mm_multipoles, &
               vacuum_aux_dcoeffs_ht(bessel_run), ngwf_basis, mdl%elements, &
               s_atoms_nl, scal_factor, &
               pub_pol_emb_dma_min_l, pub_pol_emb_dma_max_l, &
               REP_SWEX_POL_EMB_AUX_DMA_1 + bessel_run - 1)
          call sparse_axpy(vacuum_aux_polemb, cur_polemb_vac_aux, &
               1.0_DP/real(n_bessel_runs,kind=DP))
          call sparse_destroy(cur_polemb_vac_aux)
       else
          ! @@@ HACK -- split P over atoms -- HACK @@@
!          do one_atom_only = 1, par%nat
!             ! jd: P matrix contribution from current (total or polarisation) state
!             call pol_emb_calc_pmatrix(cur_polemb, overlap, &
!                  qm_driving_multipoles(bessel_run), mm_multipoles, &
!                  driving_dcoeffs_ht(bessel_run), ngwf_basis, mdl%elements, s_atoms_nl, &
!                  scal_factor, pub_pol_emb_dma_min_l, pub_pol_emb_dma_max_l, &
!                  REP_SWEX_POL_EMB_DMA_1 + bessel_run - 1, one_atom_only = one_atom_only)
!             do is = 1, pub_num_spins
!                write(*,'(a,i0,a,i0,a,i0,a,f20.10)') '@CONTR atom ', one_atom_only, ' sp ', is, ' bessr ', bessel_run, ' is ', &
!                   ind_factor * sparse_trace(denskern%m(is, PUB_1K), cur_polemb)
!             end do
!
!          end do

          ! jd: P matrix contribution from current (total or polarisation) state
          call pol_emb_calc_pmatrix(cur_polemb, overlap, &
               qm_driving_multipoles(bessel_run), mm_multipoles, &
               driving_dcoeffs_ht(bessel_run), ngwf_basis, mdl%elements, s_atoms_nl, &
               scal_factor, pub_pol_emb_dma_min_l, pub_pol_emb_dma_max_l, &
               REP_SWEX_POL_EMB_DMA_1 + bessel_run - 1)
       end if

       ! jd: Accumulate average over Bessel runs in polemb, vacuum_polemb
       call sparse_axpy(polemb, cur_polemb, 1.0_DP/real(n_bessel_runs,kind=DP))
       call sparse_destroy(cur_polemb)


    end do ! 1 or 2 Bessel runs
    ! -------------------------------------------------------------------------

!@@
#if 0
    if(pub_debug) then
       if(pub_on_root) then
          print *, ""
          print *, "-DKN:------------------------------"
       end if
       call sparse_show_matrix(denskern%m(1,PUB_1K))
       if(pub_pol_emb_vacuum_qmstar) then
          if(pub_on_root) then
             print *, ""
             print *, "-VACDKN:------------------------------"
          end if
          call sparse_show_matrix(vacuum_kdenskern_in_new_ngwfs%m(1,PUB_1K))
       end if
       if(pub_on_root) then
          print *, ""
          print *, "-OVERLAP:------------------------------"
       end if
       call sparse_show_matrix(overlap)
    end if
#endif

    ! --------------------------------------------------------------------
    ! jd: Calculate the energy of the MM multipole set interacting
    !     with electrons, from tr[KP]. Induced dipoles get a 1/2.
    ! --------------------------------------------------------------------
    elec_energy = 0D0
    polemb_energy_pol = 0D0
    polemb_energy_vac = 0D0
    polemb_energy_vac_aux = 0D0
    do is = 1, pub_num_spins
       polemb_energy_pol = polemb_energy_pol + &
            ind_factor * sparse_trace(denskern%m(is, PUB_1K), polemb)
       if(pub_pol_emb_vacuum_qmstar) then
          polemb_energy_vac = polemb_energy_vac + &
               ind_factor * sparse_trace(&
               polarisable_embedding_vacuum_state%kdenskern%kern%m(is, PUB_1K)%p, &
               vacuum_polemb)
          polemb_energy_vac_aux = polemb_energy_vac_aux - &
               ind_factor * sparse_trace(&
               polarisable_embedding_vacuum_state%kdenskern%kern%m(is, PUB_1K)%p, &
               vacuum_aux_polemb)
       end if
    end do

    elec_energy = elec_energy + polemb_energy_pol + polemb_energy_vac + polemb_energy_vac_aux
!    write(*,*) '@POL ', polemb_energy_pol
!    write(*,*) '@VAC ', polemb_energy_vac
!    write(*,*) '@AUX ', polemb_energy_vac_aux
!    write(*,*) '@TOT ', elec_energy

    call timer_clock(myself,2)

  end subroutine pol_emb_qmstar_energy_and_lnv_gradient

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pol_emb_qterm(mm_multipoles, elements, par, global_b, & ! in
       qm_min_l, qm_max_l, &                                         ! in
       M00_B, &                                                      ! in, opt
       qterm_sph, sum_corr_terms)                                    ! out, opt
    !==========================================================================!
    ! Calculates the Q term for a single QM site B in polarisable embedding.   !
    ! Q_B is an auxiliary quantity in the calculation of the P matrix.         !
    ! Also calculates S_corr, the charge scaling correction term auxiliary     !
    ! quantity. Both of Q term and sph_corr_terms are optional -- they will    !
    ! not be provided if not passed.                                           !
    !                                                                          !
    ! The elements in the qterm array are arranged according to their l and m  !
    ! values: qterm_sph(1) is the QM charge term, qterm_sph(2:4) are the QM    !
    ! dipole terms, qterm_sph (5:9) are the QM quadrupole terms. Spherical     !
    ! multipoles are used, quadrupoles are traceless and use the Gray-Gubbins  !
    ! convention.                                                              !
    ! If dma_max_l is lower than 2, the qterm_sph array has a correspondingly  !
    ! narrower range.                                                          !
    !                                                                          !
    ! The units of Q are potential (for charges), potential/length   (for      !
    ! dipoles), potential/length^2 (for quadrupoles). Q can be identified with !
    ! the scalar electric potential, electric field vector and electric field  !
    ! derivative tensor at r_B. Multiplying Q by the corresponding DMA d       !
    ! coefficients (unitless, charge * m, charge * m^2, respectively) yields   !
    ! P matrix elements, which being energy derivatives wrt DKN (charge), have !
    ! the units of potential.                                                  !
    !                                                                          !
    ! The type of underlying electrostatics (Coulombic, Thole,...) is defined  !
    ! by the MM multipole set.                                                 !
    !                                                                          !
    ! If charge-scaling is in effect, at the same time a second auxiliary      !
    ! quantity is calculated and returned in sum_corr_terms. This term is      !
    ! similar in nature -- it is the sum of potentials, fields and field deri- !
    ! vatives due to unscaled charges of all QM atoms, summed over all MM      !
    ! sites. For the calculation of this term the unscaled monopole of B must  !
    ! be passed as the last argument.                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   mm_multipoles (in): Describes the MM multipole set.                    !
    !   elements (in): Needed for QM site positions and polarisabilities.      !
    !   global_b (in): Global index B of the QM site for which Q_B is calc'ed. !
    !                  B does *not* need to be an atom local to this rank,     !
    !                  owing to the replication of elements(:).                !
    !   qterm_sph (out, opt): Result for Q term is returned here, if provided. !
    !   M00_B (in, opt): The unscaled DMA monopole on atom B. Ignored if charge!
    !                    scaling is not in effect. Only needed if              !
    !                    sum_corr_terms is also passed.                        !
    !   sum_corr_terms (in/out, opt): The auxiliary quantity is *added* to this!
    !                                 argument. This will need to be suitably  !
    !                                 comms_reduced later. Ignored if omitted. !
    !                                 Filled with garbage if charge scaling    !
    !                                 is not in effect.                        !
    !                                                                          !
    ! @doc qm_min_l, qm_max_l -- needed because in split-dma qterms with different
    !                            masks are needed
    ! @doc that cartesians are set to (remain at) zero if an l is absent in QM,
    !      but conversion to sph takes care to align it
    !--------------------------------------------------------------------------!
    ! Notes:
    ! - for the origin of ONE_THIRD, cf. eq. (2.63) in CG Gray and KE Gubbins  !
    !   "Theory of Molecular Fluids, Volume1: Fundamentals".                   !
    ! - minuses for charge and quadrupole terms are a consequence of a combina-!
    !   tion of ONETEP using a reverse charge convention, and having (-)^l in  !
    !   the definition of the potential, and not the multipole. Consequently,  !
    !   we need '-' for every odd multipole.                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October 2016.                               !
    ! Re-arranged to use spherical multipoles by Jacek Dziedzic in Nov 2016.   !
    ! Re-arranged again by Jacek Dziedzic in Feb 2017.                         !
    ! Modified to remove pub_par by Joseph Prentice, September 2018            !
    !==========================================================================!

    use constants, only: garbage_real, ONE_THIRD
    use geometry, only: POINT
    use ion, only: ELEMENT
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET, &
         multipole_0_sph_real_to_cart, multipole_1_sph_real_to_cart, &
         multipole_2_sph_real_to_cart_traceless, multipole_T0, multipole_T1, &
         multipole_T2, multipole_T3, multipole_T4
    use parallel_strategy, only: PARAL_INFO
    use multipole_ops, only: multipole_1_cart_to_sph_real, &
         multipole_2_cart_traceless_to_sph_real, &
         multipole_2_cart_primitive_to_traceless
    use rundat, only: pub_pol_emb_mpole_exclusion_radius, pub_debug, &
         pub_pol_emb_thole_a, pub_dma_scale_charge
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(PARAL_INFO), intent(in)              :: par
    type(SPHERICAL_MULTIPOLE_SET), intent(in) :: mm_multipoles
    type(ELEMENT), intent(in)                 :: elements(:)
    integer, intent(in)                       :: global_b
    integer, intent(in)                       :: qm_min_l
    integer, intent(in)                       :: qm_max_l
    real(kind=DP), intent(out), optional      :: qterm_sph(:)
    real(kind=DP), intent(in), optional       :: M00_B
    real(kind=DP), intent(inout), optional    :: sum_corr_terms

    ! jd: Local variables
    real(kind=DP) :: qterm_cart(n_qterm_cart_components)
    real(kind=DP) :: qterm_cart_traceless(6)
    integer       :: orig_global_b
    integer       :: site_l
    type(POINT)   :: vec_b, vec_l
    real(kind=DP) :: alpha_b, alpha_l
    integer       :: mm_min_l, mm_max_l
    integer       :: min_ll, max_ll, ll ! mmm
    integer       :: mpole_offs ! offset in multipole array
    integer       :: sph_offset ! offset in spherical representation of result
    real(kind=DP) :: T0_BL, T1_BL(3), T2_BL(3,3), T3_BL(3,3,3), T4_BL(3,3,3,3)
    real(kind=DP) :: q_l       ! MM charge at site L
    real(kind=DP) :: mu_l(3)   ! MM dipole at site L
    real(kind=DP) :: QQ_l(3,3) ! MM quadrupole at site L
    real(kind=DP) :: qpole_fac
    real(kind=DP), parameter :: epstol = 1D-9
    character(len=*), parameter :: myself = 'pol_emb_qterm'

    ! -------------------------------------------------------------------------

    call utils_assert(present(M00_B) .eqv. present(sum_corr_terms), myself//&
         'M00_B needs to be present iff sum_corr_terms is present')

    ! jd: Thole polarisability and position of atom B
    orig_global_b = par%orig_atom(global_b)
    alpha_b = elements(orig_global_b)%thole_polarisability
    vec_b = elements(orig_global_b)%centre

    ! jd: Range of angular momenta and their pairs we'll be interested in
    mm_min_l = minval(mm_multipoles%min_l_mask(:))
    mm_max_l = maxval(mm_multipoles%max_l_mask(:))
    min_ll = mm_min_l + qm_min_l
    max_ll = mm_max_l + qm_max_l

    qterm_cart(:) = 0D0
    qterm_sph(:) = 0D0 ! jd: Crucial, some components may remain untouched after that

    ! --------------------------------------------------------------------------
    ! jd: Loop over MM multipole sites L                                    LLL
    ! --------------------------------------------------------------------------
    loop_L:                                                                   &
    do site_l = 1, mm_multipoles%n_sites
       ! jd: Only take into account sites with unmasked dipoles.
       !     The ones that are masked are actually QM sites, with
       !     zero polarisabilities, an artifact of Jay's approach.
       !     They should be dropped here.
       if(mm_multipoles%max_l_mask(site_l) /= -1) then
          vec_l = mm_multipoles%centres(site_l)
          alpha_l = mm_multipoles%polarisabilities(site_l)

          if(pub_debug) then
             T0_BL = garbage_real
             T1_BL = garbage_real
             T2_BL = garbage_real
             T3_BL = garbage_real
             T4_BL = garbage_real
             q_l = garbage_real
             mu_l = garbage_real
             QQ_l = garbage_real
          end if

          ! --------------------------------------------------------------------
          ! jd: Loop over angular momentum pairs ll                          ll
          ! --------------------------------------------------------------------
          loop_ll:                                                            &
          do ll = min_ll, max_ll
             if(ll == 0) then
                T0_BL = multipole_T0(vec_b, vec_l, &
                    mm_multipoles%potential_damping, &
                    pub_pol_emb_mpole_exclusion_radius, &
                    alpha_b, alpha_l, pub_pol_emb_thole_a)
             end if
             if(ll == 1) then
                T1_BL(:) = multipole_T1(vec_b, vec_l, &
                    mm_multipoles%potential_damping, &
                    pub_pol_emb_mpole_exclusion_radius, &
                    alpha_b, alpha_l, pub_pol_emb_thole_a)
             end if
             if(ll == 2) then
                T2_BL(:,:) = multipole_T2(vec_b, vec_l, &
                    mm_multipoles%potential_damping, &
                    pub_pol_emb_mpole_exclusion_radius, &
                    alpha_b, alpha_l, pub_pol_emb_thole_a)
             end if
             if(ll == 3) then
                T3_BL(:,:,:) = multipole_T3(vec_b, vec_l, &
                    mm_multipoles%potential_damping, &
                    pub_pol_emb_mpole_exclusion_radius, &
                    alpha_b, alpha_l, pub_pol_emb_thole_a)
             end if
             if(ll == 4) then
                T4_BL(:,:,:,:) = multipole_T4(vec_b, vec_l, &
                    mm_multipoles%potential_damping, &
                    pub_pol_emb_mpole_exclusion_radius, &
                    alpha_b, alpha_l, pub_pol_emb_thole_a)
             end if
          end do loop_ll
       end if ! MM multipole not entirely masked

       ! **********************
       ! jd: --- MM charges ---
       ! **********************
       if(mm_multipoles%max_l_mask(site_l) >= 0 .and. &
            mm_multipoles%min_l_mask(site_l) <= 0) then
          ! jd: Get \q_L
          mpole_offs = (site_l-1) * mm_multipoles%n_multipoles_per_site + 0
          call multipole_0_sph_real_to_cart(&
               q_l, mm_multipoles%mpoles(mpole_offs+1), .true.)
          ! jd: QM charge <-> MM charge
          if(qm_min_l <= 0 .and. qm_max_l >= 0) then
             qterm_cart(1) = qterm_cart(1) + (-q_l * T0_BL)
          end if
          ! jd: QM dipole <-> MM charge
          if(qm_min_l <= 1 .and. qm_max_l >= 1) then
             qterm_cart(2:4) = qterm_cart(2:4) + (-q_l * T1_BL(:))
          end if
          ! jd: QM quadrupole <-> MM charge
          if(qm_min_l <=1 .and. qm_max_l >= 2) then
             qterm_cart(5:13) = qterm_cart(5:13) + &
                  ONE_THIRD * reshape((-q_l * T2_BL(:,:)),(/9/))
          end if
       end if

       ! **********************
       ! jd: --- MM dipoles ---
       ! **********************
       if(mm_multipoles%max_l_mask(site_l) >= 1 .and. &
            mm_multipoles%min_l_mask(site_l) <= 1) then
          ! jd: Get \mu_L
          mpole_offs = (site_l-1) * &
               mm_multipoles%n_multipoles_per_site + 1
          call multipole_1_sph_real_to_cart(&
               mu_l(1), mu_l(2), mu_l(3), &
               mm_multipoles%mpoles(mpole_offs+1), &
               mm_multipoles%mpoles(mpole_offs+2), &
               mm_multipoles%mpoles(mpole_offs+3), .true.)
          ! jd: QM charge <-> MM dipole
          if(qm_min_l <= 0 .and. qm_max_l >= 0) then
             qterm_cart(1) = qterm_cart(1) + sum(mu_l(:) * T1_BL(:))
          end if
          ! jd: QM dipole <-> MM dipole
          if(qm_min_l <= 1 .and. qm_max_l >= 1) then
             qterm_cart(2:4) = qterm_cart(2:4) + (matmul(mu_l(:),T2_BL(:,:)))
          end if
          ! jd: QM quadrupole <-> MM dipole
          if(qm_min_l <= 2 .and. qm_max_l >= 2) then
             qterm_cart(5:7)   = qterm_cart(5:7) + &
                  ONE_THIRD * matmul(mu_l(:),T3_BL(1,:,:))
             qterm_cart(8:10)  = qterm_cart(8:10) + &
                  ONE_THIRD * matmul(mu_l(:),T3_BL(2,:,:))
             qterm_cart(11:13) = qterm_cart(11:13) + &
                  ONE_THIRD * matmul(mu_l(:),T3_BL(3,:,:))
          end if
       end if

       ! **************************
       ! jd: --- MM quadrupoles ---
       ! **************************
       if(mm_multipoles%max_l_mask(site_l) >= 2 .and. &
            mm_multipoles%min_l_mask(site_l) <= 2) then
          ! jd: Get \Q_L
          mpole_offs = (site_l-1) * &
               mm_multipoles%n_multipoles_per_site + 4
          call multipole_2_sph_real_to_cart_traceless(&
               QQ_l(1,1), QQ_l(1,2), QQ_l(1,3), QQ_l(2,2), QQ_l(2,3), QQ_l(3,3),&
               mm_multipoles%mpoles(mpole_offs+1), &
               mm_multipoles%mpoles(mpole_offs+2), &
               mm_multipoles%mpoles(mpole_offs+3), &
               mm_multipoles%mpoles(mpole_offs+4), &
               mm_multipoles%mpoles(mpole_offs+5), .true.)
          QQ_l(2,1) = QQ_l(1,2)
          QQ_l(3,1) = QQ_l(1,3)
          QQ_l(3,2) = QQ_l(2,3)
          ! jd: QM charge <-> MM quadrupole
          if(qm_min_l <= 0 .and. qm_max_l >= 0) then
             qterm_cart(1) = &
                  qterm_cart(1) + (-ONE_THIRD) * sum(QQ_l(:,:) * T2_BL(:,:))
          end if
          ! jd: QM dipole <-> MM quadrupole
          if(qm_min_l <= 1 .and. qm_max_l >= 1) then
             qterm_cart(2) = &
                  qterm_cart(2) + (-ONE_THIRD) * sum(QQ_l(:,:) * T3_BL(1,:,:))
             qterm_cart(3) = &
                  qterm_cart(3) + (-ONE_THIRD) * sum(QQ_l(:,:) * T3_BL(2,:,:))
             qterm_cart(4) = &
                  qterm_cart(4) + (-ONE_THIRD) * sum(QQ_l(:,:) * T3_BL(3,:,:))
          end if
          ! jd: QM quadrupole <-> MM quadrupole
          if(qm_min_l <= 2 .and. qm_max_l >= 2) then
             qpole_fac = (-ONE_THIRD) * ONE_THIRD
             qterm_cart(5) = &
                  qterm_cart(5)   + qpole_fac * sum(QQ_l(:,:) * T4_BL(1,1,:,:))
             qterm_cart(6) = &
                  qterm_cart(6)   + qpole_fac * sum(QQ_l(:,:) * T4_BL(2,1,:,:))
             qterm_cart(7) = &
                  qterm_cart(7)   + qpole_fac * sum(QQ_l(:,:) * T4_BL(3,1,:,:))
             qterm_cart(8) = &
                  qterm_cart(8)   + qpole_fac * sum(QQ_l(:,:) * T4_BL(1,2,:,:))
             qterm_cart(9) = &
                  qterm_cart(9)   + qpole_fac * sum(QQ_l(:,:) * T4_BL(2,2,:,:))
             qterm_cart(10) = &
                  qterm_cart(10) + qpole_fac * sum(QQ_l(:,:) * T4_BL(3,2,:,:))
             qterm_cart(11) = &
                  qterm_cart(11) + qpole_fac * sum(QQ_l(:,:) * T4_BL(1,3,:,:))
             qterm_cart(12) = &
                  qterm_cart(12) + qpole_fac * sum(QQ_l(:,:) * T4_BL(2,3,:,:))
             qterm_cart(13) = &
                  qterm_cart(13) + qpole_fac * sum(QQ_l(:,:) * T4_BL(3,3,:,:))
          end if
       end if

       ! jd: Sanity check: Cartesian quadrupole components should have the
       !     following symmetries: xy<->yx, xz<->zx, yz<->zy.
       if(abs(qterm_cart(6)-qterm_cart(8)) > epstol) then
          call utils_assert(.false., 'pol_emb_qterm: Spontaneous XY symmetry &
               &breaking in Cartesian quadrupole.',qterm_cart(6)-qterm_cart(8))
       end if
       if(abs(qterm_cart(7)-qterm_cart(11)) > epstol) then
          call utils_assert(.false., 'pol_emb_qterm: Spontaneous XZ symmetry &
               &breaking in Cartesian quadrupole.',qterm_cart(11)-qterm_cart(7))
       end if
       if(abs(qterm_cart(10)-qterm_cart(12)) > epstol) then
          call utils_assert(.false., 'pol_emb_qterm: Spontaneous XY symmetry &
               &breaking in Cartesian quadrupole.',qterm_cart(10)-qterm_cart(12))
       end if

       ! jd: Convert primitive       5   8   11      1  2  3
       !     Cartesian quadrupole    6   9   12  ->  2  4  5
       !     to traceless            7  10   13      3  5  6
       call multipole_2_cart_primitive_to_traceless(&
            qterm_cart_traceless(1), qterm_cart_traceless(2), &
            qterm_cart_traceless(3), qterm_cart_traceless(4), &
            qterm_cart_traceless(5), qterm_cart_traceless(6), &
            qterm_cart(5), qterm_cart(6), qterm_cart(7), qterm_cart(9), &
            qterm_cart(10), qterm_cart(13), convention='GrayGubbins')
       if(abs(qterm_cart_traceless(1) + qterm_cart_traceless(4) + &
            qterm_cart_traceless(6)) > epstol) then
          call utils_assert(.false., 'pol_emb_qterm: Cartesian quadrupole not &
               &traceless.', qterm_cart_traceless(1)+qterm_cart_traceless(4) + &
               qterm_cart_traceless(6))
       end if

       ! jd: Convert qterm from Cartesian to spherical
       if(present(qterm_sph)) then
          sph_offset = 1
          if(qm_min_l <= 0 .and. qm_max_l >=0) then
             qterm_sph(sph_offset) = qterm_cart(sph_offset)
             sph_offset = sph_offset + 1
          end if
          if(qm_min_l <= 1 .and. qm_max_l >=1) then
             call multipole_1_cart_to_sph_real(qterm_sph(sph_offset), &
                  qterm_sph(sph_offset+1), qterm_sph(sph_offset+2), &
                  qterm_cart(2), qterm_cart(3), qterm_cart(4), .true.)
             sph_offset = sph_offset + 3
          end if
          if(qm_min_l <= 2 .and. qm_max_l >=2) then
             call multipole_2_cart_traceless_to_sph_real(qterm_sph(sph_offset),&
                  qterm_sph(sph_offset+1), qterm_sph(sph_offset+2), &
                  qterm_sph(sph_offset+3), qterm_sph(sph_offset+4), &
                  qterm_cart_traceless(1), qterm_cart_traceless(2), &
                  qterm_cart_traceless(3), qterm_cart_traceless(4), &
                  qterm_cart_traceless(5), qterm_cart_traceless(6), .true.)
          end if
       end if

       if(present(sum_corr_terms)) then
          if(pub_dma_scale_charge) then
             call utils_assert(qm_min_l <= 0, &
                 myself//': Cannot charge-scale a set with no charges')
             ! jd: Only take into account sites with unmasked mpoles.
             !     Partial masking is used to avoid calculating zeroes.
             !     Full masking excludes mpoles that are actually QM sites,
             !     with zero polarisabilities, an artifact of Jay's approach.
             !     They are dropped here.
             if(mm_multipoles%max_l_mask(site_l) >= 0 .and. &
                  mm_multipoles%min_l_mask(site_l) <= 0) then
                sum_corr_terms = sum_corr_terms + M00_B * T0_BL * (-q_l)
             end if
             if(mm_multipoles%max_l_mask(site_l) >= 1 .and. &
                  mm_multipoles%min_l_mask(site_l) <= 1) then
                sum_corr_terms = sum_corr_terms + M00_B * sum(T1_BL(:) * mu_l(:))
             end if
             if(mm_multipoles%max_l_mask(site_l) >= 2 .and. &
                  mm_multipoles%min_l_mask(site_l) <= 2) then
                sum_corr_terms = sum_corr_terms + &
                     M00_B * (-ONE_THIRD * sum(T2_BL(:,:) * QQ_l(:,:)))
             end if
          else ! jd: Must write something there: intent(out)
             sum_corr_terms = garbage_real
          end if
       end if
    end do loop_L

  end subroutine pol_emb_qterm

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pol_emb_jqterm(jqterm_sph, qterm_sph, swri, dma_swex)
    !==========================================================================!
    ! Calculates the "JQ term", that is Q_A^lm * J_lq, where the first term    !
    ! is the "Q term" on an atom A, and the second term is the Bessel radial   !
    ! integral. This is evaluated for all SWs originating on A, that is for    !
    ! s_A = {l,m,q}.                                                           !
    !                                                                          !
    ! If Bessel averaging is in effect, this subroutine needs to be called     !
    ! separately for each of the two swexes.                                   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   jqterm_sph (out): The result is returned here.                         !
    !   qterm_sph (in): Values of Q term calculated in advance, used as input. !
    !   swri (in): SWRI in which we operate. Needed for index bookkeeping.     !
    !   dma_swex (in): The SWEX describing the expansion.                      !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2016.                              !
    !==========================================================================!

    use constants, only: SQRT_PI, PI
    use rundat, only: pub_pol_emb_dma_max_l, pub_pol_emb_vacuum_dma_max_l, &
         pub_pol_emb_dma_min_l, pub_pol_emb_vacuum_dma_min_l
    use sw_resolution_of_identity, only: SW_RI, swri_sph_bess_solidharm_int, &
         swri_sw_to_ablmq
    use sw_expansion_type, only: SW_EX
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: jqterm_sph(:)
    real(kind=DP), intent(in)  :: qterm_sph(:)
    type(SW_RI), intent(in)    :: swri
    type(SW_EX), intent(in)    :: dma_swex

    ! jd: Local variables
    real(kind=DP) :: qterm_cart(n_qterm_cart_components)
    integer       :: swex_sw, swri_sw
    integer       :: num_sws_per_centre
    integer       :: i_b, i_l, i_m, qterm_idx
    integer       :: min_l
    real(kind=DP) :: a, q
    real(kind=DP) :: J_lq
    real(kind=DP) :: factor
    character(len=*), parameter :: myself = 'pol_emb_jqterm'

    ! -------------------------------------------------------------------------

    num_sws_per_centre = dma_swex%quality%num_sws_per_centre
    min_l = dma_swex%quality%min_l
    jqterm_sph(:) = 0D0 ! jd: As not necessarily all elements are filled
                        !     (when Bessel averaging is used)
    do swex_sw = 1, num_sws_per_centre
       swri_sw = dma_swex%swex_sw_to_swri_sw(swex_sw)
       call swri_sw_to_ablmq(swri, swri_sw, a, i_b ,i_l ,i_m, q)
       qterm_idx = i_l*i_l - min_l*min_l + i_l + 1 + i_m ! jd: Indexes lm's
!       write(*,'(a,i5,i5,i5,i5,i5)') '@@@XLAT: ', swex_sw, swri_sw, i_l, i_m,  qterm_idx

       call utils_assert(i_l <= max(pub_pol_emb_dma_max_l, pub_pol_emb_vacuum_dma_max_l), &
            myself//': Internal error (l too big)',i_l,max(pub_pol_emb_dma_max_l, pub_pol_emb_vacuum_dma_max_l))
       call utils_assert(i_l >= min(pub_pol_emb_dma_min_l, pub_pol_emb_vacuum_dma_min_l), &
            myself//': Internal error (l too small)',i_l,min(pub_pol_emb_dma_min_l, pub_pol_emb_vacuum_dma_min_l))
!       call utils_assert(i_l <= pub_pol_emb_dma_max_l, &
!            myself//': Internal error (l)',i_l,pub_pol_emb_dma_max_l)
       call utils_assert(i_m >= -i_l .and. i_m <= i_l, &
            myself//': Internal error (m)',i_m)

       ! jd: Room for @optimisation: J_lq does not need to be
       !     recalculated if only m changes.
       ! jd: Always keep sure that this is consistent with
       !     its counterpart in compute_multipole_elms()
       J_lq = swri_sph_bess_solidharm_int(i_l,q,a)
       ! jd: Employ Racah normalisation
       factor=2D0*SQRT_PI/sqrt(2D0*i_l+1D0)
       ! jd: Required fudge factor, origin unknown
       factor=factor/(2.0_DP*PI)

       jqterm_sph(swex_sw) = J_lq * factor * qterm_sph(qterm_idx)
    end do

  end subroutine pol_emb_jqterm

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pol_emb_uabterm(swri,                                     & ! inout
       ngwf_basis, cell, elements, par, local_a, global_b, jqterms_ht, & ! in
       swex, bessel_run, u00_only, &                                     ! in
       uabterm_in_ppds_on_a,               & ! adds to } PPD version
       ppd_indices_on_a, n_ppds_on_a,      & ! out     } PPD version
       uabterm_in_fftbox_dbl_on_a, &         ! adds to } FFTbox version
       fftbox, hack_mask)                    ! in      } FFTbox version
    !==========================================================================!
    ! Calculates the "U_AB" term:                                              !
    ! U_AB(\vec{r}) =                                                          !
    ! \sum_{s_A} JQ_A^{s_A} \sum_{t_AB} v_{t_AB} (\vec(r)) V^{t_AB,s_A} +      !
    ! \sum_{s_B} JQ_B^{s_B} \sum_{t_AB} v_{t_AB} (\vec(r)) V^{t_AB,s_B},       !
    ! where s_A, s_B are spherical waves on atoms A, B. JQ_A and JQ_B are the  !
    ! "JQ terms" on atoms A, B. t_AB are the spherical waves on atoms A, B,    !
    ! v are spherical wave potentials, V is the inverse V matrix block between !
    ! atoms A, B.                                                              !
    !                                                                          !
    ! If u00_only is .true., limits the sums over s_A and s_B (but not t_AB)   !
    ! to charge terms. This is useful when calculating the charge-scaling      !
    ! correction to the NGWF gradient.                                         !
    !                                                                          !
    ! Atom A must be owned by the rank on which this is executed.              !
    !                                                                          !
    ! In a typical scenario the PPD arguments are passed, and U_AB is calcula- !
    ! ted in PPDs on atom A, on the coarse grid.                               !
    !                                                                          !
    ! If pol_emb_dbl_grid T instead, U_AB is calculated in a double-grid FFT   !
    ! box around atom A, with uabterm_in_ppds_on_a unreferenced, *but* the PPDs!
    ! are still determined, so ppd_indices_on_a and n_ppds_on_a are always     !
    ! filled.                                                                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   uabterm_in_ppds_on_a (inout): Output is accumulated here, unless       !
    !                                 pol_emb_fine_grid T (then this remains   !
    !                                 unreferenced).                           !
    !   ppd_indices_on_a (out): Will be filled with the indices of PPDs where  !
    !                           the output is returned.                        !
    !   n_ppds_on_a (out): Will be set to the number of PPDs in output.        !
    !   swri (inout): SWRI in which we operate. In/out, because SWOP caching   !
    !                 takes place behind the scenes.                           !
    !   ngwf_basis, cell (in): Needed to dimension PPD arrays.                 !
    !   elements (in): Needed for ionic positions.                             !
    !   local_a (in): Local index of atom A.                                   !
    !   global_b (in): Global index of atom B.                                 !
    !   jqterms_ht (in): Stores values of JQ_A and JQ_B calculated in advance. !
    !   swex (in): SWEX in which we operate.                                   !
    !   bessel_run (in): Needed to pick the right jqvterm from the HT.         !
    !   u00_only (in): If .true., returns U^{00}_AB. If .false., returns U_AB. !
    !   uabterm_in_fftbox_dbl_on_a (inout): Output is accumulated here iff     !
    !                                       pol_emb_dbl_grid T.                !
    !   fftbox (in): Needed iff pol_emb_dbl_grid T.                            !
    !--------------------------------------------------------------------------!
    ! Notes:                                                                   !
    ! Values are *accumulated* in uabterm_in_ppds_on_a, so it's the responsi-  !
    ! bility of the caller to zero this in advance.                            !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2016.                              !
    ! Extended by Jacek Dziedzic in January 2017 with the double grid version. !
    ! Modified to remove pub_par by Joseph Prentice, September 2018            !
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: SW_V, SW_O, METRIC_ELECTROSTATIC, VERBOSE
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_batch_col_start
    use geometry, only: POINT, operator(*), operator(-)
! jme: workaround for cray compiler bug 843178
#if _CRAYFTN
    use geometry, only: add_points
#else
    use geometry, only: operator(+)
#endif
    use hash_table, only: hash_table_stash_prepare, hash_table_stash_commit, &
         hash_table_stash_free, hash_table_lookup_nocount, HT_HASH_TABLE
    use ion, only: ELEMENT
    use linalg, only: linalg_invert_sym_matrix
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_dma_metric, pub_output_detail, pub_pol_emb_dbl_grid,&
         pub_threads_max, pub_dbl_grid_scale, pub_dma_use_ri, pub_devel_code
    use simulation_cell, only: CELL_INFO
    use sw_resolution_of_identity, only: SW_RI, ATOM_CENTRE, swri_init_centre, &
         swri_expansion_centres, swri_extract_matrix_2, swri_sw_to_ablmq, &
         swri_obtain_swops_in_ppd, swri_obtain_swops_at_point, &
         swri_get_handle_to
    use sw_expansion_type, only: SW_EX
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check, &
         utils_abort, utils_devel_code

    implicit none

    ! jd: Arguments
    type(PARAL_INFO), intent(in)            :: par
    type(FUNC_BASIS), intent(in)            :: ngwf_basis
    type(CELL_INFO), intent(in)             :: cell
    type(ELEMENT), intent(in)               :: elements(:)
    integer, intent(in)                     :: local_a
    integer, intent(in)                     :: global_b
    type(HT_HASH_TABLE), intent(in), target :: jqterms_ht
    type(SW_RI), intent(inout)              :: swri
    type(SW_EX), intent(in)                 :: swex
    integer, intent(in)                     :: bessel_run
    logical, intent(in)                     :: u00_only
    logical, intent(in), optional           :: hack_mask

    real(kind=DP), intent(inout), optional :: uabterm_in_ppds_on_a(PPDS_DIMS)
    integer, intent(out) :: ppd_indices_on_a(ngwf_basis%max_n_ppds_sphere)
    integer, intent(out) :: n_ppds_on_a

    type(FFTBOX_INFO), intent(in), optional :: fftbox
    real(kind=DP), intent(inout), optional  :: uabterm_in_fftbox_dbl_on_a(:,:,:)

    ! jd: Local variables
    integer :: orig_global_b
    integer :: first_ngwf_idx_of_a, local_idx_of_1st_ngwf_on_a
    integer :: global_a, orig_global_a
    integer :: num_sws_in_expansion
    integer :: expansion_atoms_out(2)
    integer :: expansion_atoms_in(2)
    type(ATOM_CENTRE) :: centres(2)
    type(ATOM_CENTRE) :: expansion_centres(2)
    type(ATOM_CENTRE) :: src_centre
    real(kind=DP), allocatable  :: metricmatrix(:,:)
    real(kind=DP) :: cond
    integer :: row, col
    integer :: swex_sw
    integer :: ierr
    integer :: swex_sw_t_ab
    integer :: ec
    integer :: src_atom
    integer :: n_sws, num_sws
    integer :: ppd_on_a
    integer :: ppd_idx
    integer :: src_offs, dest_offs
    integer :: swri_sw_s_a, swex_sw_s_a, swex_sw_s_a_min, swex_sw_s_a_max
    integer :: i_l, i_q, i_m, i_b
    real(kind=DP) :: q, a
    character :: sw_or_swpot
    integer :: didx, d1idx, d2idx, d3idx
    integer :: npts_dbl
    integer :: dma_swri_h
    real(kind=DP) :: jqterm_a(swex%quality%num_sws_per_centre)
    real(kind=DP) :: jqterm_b(swex%quality%num_sws_per_centre)
    integer :: n_elems
    integer :: fa_box_start(3), fa_start_in_box(3), fftbox_start_dbl(3)
    type(POINT) :: fftbox_dbl_start_in_cell_vec
    type(POINT)   :: a1_dbl, a2_dbl, a3_dbl
    real(kind=DP) :: common_factor
    real(kind=DP) :: jqvterms(2*swex%quality%max_sws_per_centre)
    real(kind=DP) :: swops_in_ppd(cell%n_pts *2*swex%quality%max_sws_per_centre)
    integer       :: points_used
    type(POINT)   :: curpoint
    integer       :: centre_sw
    real(kind=DP) :: swops_at_point(swex%quality%num_sws_per_centre)
    real(kind=DP), allocatable :: swops_stash(:)
    logical :: debug_show_skipped
    integer, parameter :: max_stash_size = 1024 ! == 1GiB, generalizeme, perthread
    real(kind=DP), parameter :: epstol = 1D-10
    character(len=*), parameter :: myself = 'pol_emb_uabterm'

    ! -------------------------------------------------------------------------

    debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
         pub_devel_code, no_warn=.true., no_bcast=.true.)
    dma_swri_h = swri_get_handle_to(pub_dma_use_ri)

    if(pub_dma_metric == METRIC_ELECTROSTATIC) then
       sw_or_swpot = 'P'
    else
       sw_or_swpot = 'S'
    end if

    global_a = par%first_atom_on_proc(pub_my_proc_id) + local_a-1
    first_ngwf_idx_of_a = ngwf_basis%first_on_atom(global_a)
    local_idx_of_1st_ngwf_on_a = first_ngwf_idx_of_a + 1 - &
         ngwf_basis%first_on_proc(pub_my_proc_id)
    orig_global_a = par%orig_atom(global_a)

    if(pub_pol_emb_dbl_grid) then

       ! jd: Determine the position of the FFTbox around Aa and the position
       !     of Aa in the FFTbox
       call function_ops_batch_col_start(fa_box_start, fa_start_in_box, &
            1, local_idx_of_1st_ngwf_on_a, local_idx_of_1st_ngwf_on_a, &
            fftbox, cell, ngwf_basis)

       ! jd: Determine where the double grid FFT box starts
       if (abs(pub_dbl_grid_scale-2.0_DP) < epstol) then
          fftbox_start_dbl(:) = 2*fa_box_start(:) - 1
       else
          call utils_abort(myself//': only double grid scale of 2 is supported')
       end if

! jme: workaround for cray compiler bug 843178
#if _CRAYFTN
       fftbox_dbl_start_in_cell_vec = 1D0 / pub_dbl_grid_scale * ( &
            add_points(add_points( &
            ((fftbox_start_dbl(1)-1) * fftbox%d1) * fftbox%a1_unit, &
            ((fftbox_start_dbl(2)-1) * fftbox%d2) * fftbox%a2_unit), &
            ((fftbox_start_dbl(3)-1) * fftbox%d3) * fftbox%a3_unit) )
#else
       fftbox_dbl_start_in_cell_vec = 1D0 / pub_dbl_grid_scale * (&
            real(fftbox_start_dbl(1)-1,kind=DP) * fftbox%d1 * fftbox%a1_unit+&
            real(fftbox_start_dbl(2)-1,kind=DP) * fftbox%d2 * fftbox%a2_unit+&
            real(fftbox_start_dbl(3)-1,kind=DP) * fftbox%d3 * fftbox%a3_unit)
#endif
    end if

    ! jd: Determine the PPDs that belong to A
    n_ppds_on_a = ngwf_basis%spheres(local_idx_of_1st_ngwf_on_a)%n_ppds_sphere
    ppd_indices_on_a(1:n_ppds_on_a) = ngwf_basis%spheres(&
         local_idx_of_1st_ngwf_on_a)%ppd_list(1,1:n_ppds_on_a)
    ppd_indices_on_a(n_ppds_on_a+1:) = -1

    ! jd: Extract V matrix atomblock
    expansion_atoms_in = (/global_a, global_b/)
    call swri_init_centre(centres(1), par, elements, global_a)
    call swri_init_centre(centres(2), par, elements, global_b)
    call swri_expansion_centres(swri, &
         expansion_centres, expansion_atoms_out, &  ! out
         num_sws_in_expansion, &                    ! out
         centres, expansion_atoms_in, swex%quality, & ! in
         .true.) ! @ <-This assumes the SW_EX to be Bb-Cc symmetric
                 ! @   Not so for mixed NGWF bases, but this doesn't happen in
                 ! @   polemb yet.
    allocate(metricmatrix(num_sws_in_expansion,num_sws_in_expansion), stat=ierr)
    call utils_alloc_check(myself,'metricmatrix',ierr)
    call swri_extract_matrix_2(metricmatrix, swri, num_sws_in_expansion, &
         expansion_centres, expansion_atoms_out, merge(SW_V, SW_O, &
         pub_dma_metric == METRIC_ELECTROSTATIC), swex%quality)

    ! jd: Invert V atomblock
    call linalg_invert_sym_matrix(metricmatrix, num_sws_in_expansion, cond)
    if(pub_output_detail > VERBOSE) then
       write(stdout,'(a,e8.2)') '_Condition number for matrix inversion: ',cond
    end if

    num_sws = swex%quality%num_sws_per_centre

    ! jd: Go over all s_A by default, unless asked to calculate U^{00} only.
    !     Then have s_A only span {l=0, m=0, q=any}.
    if(u00_only) then
       ! jd: NB that Using 'max_q' is not viable, because we don't know if last
       !     Bessels were dropped. And so we count.
       do swex_sw_s_a = 1, num_sws
          swri_sw_s_a = swex%swex_sw_to_swri_sw(swex_sw_s_a)
          call swri_sw_to_ablmq(swri, swri_sw_s_a, a, i_b ,i_l ,i_m, q)
          if(i_l /= 0) exit
       end do
       swex_sw_s_a_max = swex_sw_s_a - 1
    else
       swex_sw_s_a_max = num_sws
    end if

    ! jd: Go over all s_A by default, unless asked to calculate U^{00} only.
    !     Then have s_A only span {l=0, m=0, q=any}.
    if(present(hack_mask)) then
       if(hack_mask) then
          call utils_assert(.not. u00_only, myself//': u00_only not supported with hack_mask')
          ! jd: NB that Using 'max_q' is not viable, because we don't know if last
          !     Bessels were dropped. And so we count.
          ! jd: Find the first l==1
          do swex_sw_s_a = 1, num_sws
             swri_sw_s_a = swex%swex_sw_to_swri_sw(swex_sw_s_a)
             call swri_sw_to_ablmq(swri, swri_sw_s_a, a, i_b ,i_l ,i_m, q)
             if(i_l == 1) exit
          end do
          swex_sw_s_a_min = swex_sw_s_a ! sic, no -1
          ! jd: Find the first l==2, subtract one, get last l==1
          do swex_sw_s_a = 1, num_sws
             swri_sw_s_a = swex%swex_sw_to_swri_sw(swex_sw_s_a)
             call swri_sw_to_ablmq(swri, swri_sw_s_a, a, i_b ,i_l ,i_m, q)
             if(i_l == 2) exit
          end do
          swex_sw_s_a_max = swex_sw_s_a - 1 ! sic,-1
       else
          swex_sw_s_a_min = 1
          swex_sw_s_a_max = num_sws
       end if
    else
       swex_sw_s_a_min = 1
       swex_sw_s_a_max = num_sws
    end if

!    if(present(hack_mask)) then
!       write(*,*) '@HACK-UAB: ', present(hack_mask), swex_sw_s_a_min, swex_sw_s_a_max
!    end if

    ! jd: Ignore atoms not represented in this SWRI
    if(.not. elements(orig_global_a)%in_swri(dma_swri_h)) then
       if(debug_show_skipped) then
          write(stdout,'(a,i0,a,i0,a,a,a)') &
          '"Skipping" A ('//trim(myself)//' orig/SFC: ', orig_global_a, '/', &
          global_a, ' (', trim(elements(orig_global_a)%species_id),')'
       end if
       jqterm_a(:) = 0D0
    else
       call hash_table_lookup_nocount(jqterm_a, n_elems, jqterms_ht, global_a, &
            bessel_run)
       call utils_assert(n_elems == swex%quality%num_sws_per_centre, &
            myself//': Inconsistent number of elements in jqterm_ht (A)', &
            global_a, n_elems, swex%quality%num_sws_per_centre, bessel_run)
    end if

    orig_global_b = par%orig_atom(global_b)
    if(global_a /= global_b) then
       ! jd: Ignore atoms not represented in this SWRI
       if(.not. elements(orig_global_b)%in_swri(dma_swri_h)) then
          if(debug_show_skipped) then
             write(stdout,'(a,i0,a,i0,a,a,a)') &
             '"Skipping" B ('//trim(myself)//' orig/SFC: ', orig_global_b, '/',&
             global_b, ' (', trim(elements(orig_global_b)%species_id),')'
          end if
          jqterm_b(:) = 0D0
       else
          call hash_table_lookup_nocount(jqterm_b, n_elems, jqterms_ht, &
               global_b, bessel_run)
          call utils_assert(n_elems == swex%quality%num_sws_per_centre, &
               myself//': Inconsistent number of elements in jqterm_ht (B)', &
               global_b, n_elems, swex%quality%num_sws_per_centre, bessel_run)
       end if
    else
       ! jd: if global_a == global_b, jqterm_b(:) is not referenced
    end if

    ! jd: Calculate J_lq * Q_AB * V^-1
    if(global_a == global_b) then
       ! jd: Single-centre expansion
       do swex_sw_t_ab = 1, num_sws_in_expansion
          jqvterms(swex_sw_t_ab) = &
               sum(metricmatrix(swex_sw_s_a_min:swex_sw_s_a_max,swex_sw_t_ab) * &
               jqterm_a(swex_sw_s_a_min:swex_sw_s_a_max))
       end do
    elseif(global_b > global_a) then
       ! jd: Two-centre expansion with Bb>Cc
       do swex_sw_t_ab = 1, num_sws_in_expansion
          jqvterms(swex_sw_t_ab) = &
               sum(metricmatrix(swex_sw_s_a_min:swex_sw_s_a_max,swex_sw_t_ab) * &
               jqterm_a(swex_sw_s_a_min:swex_sw_s_a_max)) + &
               sum(metricmatrix(num_sws+swex_sw_s_a_min:&
               num_sws+swex_sw_s_a_max,swex_sw_t_ab) * &
               jqterm_b(swex_sw_s_a_min:swex_sw_s_a_max))
       end do
    else
       ! jd: Two-centre expansion with Bb<Cc
       do swex_sw_t_ab = 1, num_sws_in_expansion
          jqvterms(swex_sw_t_ab) = &
              sum(metricmatrix(swex_sw_s_a_min:swex_sw_s_a_max,swex_sw_t_ab) * &
               jqterm_b(swex_sw_s_a_min:swex_sw_s_a_max)) + &
               sum(metricmatrix(num_sws+swex_sw_s_a_min:&
               num_sws+swex_sw_s_a_max,swex_sw_t_ab) * &
               jqterm_a(swex_sw_s_a_min:swex_sw_s_a_max))
       end do
    end if

    deallocate(metricmatrix, stat=ierr)
    call utils_dealloc_check(myself,'metricmatrix',ierr)

    if(.not. pub_pol_emb_dbl_grid) then
       ! - COARSE - COARSE - COARSE - COARSE - COARSE - COARSE - COARSE - COARSE
       ! -----------------------------------------------------------------------
       ! jd: Loop over PPDs on A                                            PPD
       ! -----------------------------------------------------------------------
       ppd_loop:                                                              &
       do ppd_on_a = 1, n_ppds_on_a
          ppd_idx = ppd_indices_on_a(ppd_on_a)
          dest_offs = (ppd_on_a-1)*cell%n_pts+1
          ! --------------------------------------------------------------------
          ! jd: Loop over expansion centres                                 EEE
          ! --------------------------------------------------------------------
          src_offs = 1
          two_centres_ppd_branch:                                                &
          do ec = 1, 2
             src_atom = expansion_atoms_out(ec)
             src_centre = expansion_centres(ec)
             if(src_atom == -1) exit

             call swri_obtain_swops_in_ppd(swri, &                     ! inout
                  swops_in_ppd(src_offs:), n_sws, &                    ! out
                  src_centre, src_atom, ppd_idx, sw_or_swpot, cell, &  ! in
                  swex%quality)                                        ! in

             if(n_sws /= num_sws) call utils_assert(.false., &
                  myself//': Logic error', n_sws, num_sws)

             ! jd: Skip to second centre, if any
             src_offs = src_offs + n_sws * cell%n_pts

          end do two_centres_ppd_branch
          ! --------------------------------------------------------------------

          src_offs = 1
          do swex_sw_t_ab = 1, num_sws_in_expansion
             uabterm_in_ppds_on_a(dest_offs:dest_offs+cell%n_pts-1) = &
                  uabterm_in_ppds_on_a(dest_offs:dest_offs+cell%n_pts-1) + &
                  swops_in_ppd(src_offs:src_offs+cell%n_pts-1) * &
                  jqvterms(swex_sw_t_ab)
             src_offs = src_offs + cell%n_pts
          end do

       end do ppd_loop
       ! -----------------------------------------------------------------------
    else
       ! - FINE - FINE - FINE - FINE - FINE - FINE - FINE - FINE - FINE - FINE
       common_factor = 1.0_DP/8.0_DP ! accounts for double grid
       npts_dbl = &
            fftbox%total_pt1_dbl * fftbox%total_pt2_dbl * fftbox%total_pt3_dbl

       a1_dbl = 0.5_DP * fftbox%d1 * fftbox%a1_unit
       a2_dbl = 0.5_DP * fftbox%d2 * fftbox%a2_unit
       a3_dbl = 0.5_DP * fftbox%d3 * fftbox%a3_unit

       points_used = 0
       src_offs = 1
       ! -----------------------------------------------------------------------
       ! jd: Loop over expansion centres                                    EEE
       ! -----------------------------------------------------------------------
       two_centres_fftbox_branch:                                             &
       do ec = 1, 2
          src_atom = expansion_atoms_out(ec)
          src_centre = expansion_centres(ec)
          if(src_atom == -1) exit

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP FIRSTPRIVATE(swops_stash) &
!$OMP PRIVATE(didx, curpoint, swops_at_point, d1idx, d2idx, d3idx, centre_sw) &
!$OMP SHARED(fftbox_dbl_start_in_cell_vec, npts_dbl, a1_dbl, a2_dbl, a3_dbl, &
!$OMP      swex, fftbox, src_centre, src_atom, src_offs, jqvterms, &
!$OMP      sw_or_swpot, swri, cell, uabterm_in_fftbox_dbl_on_a) &
!$OMP REDUCTION(+:points_used)

          ! jd: Have each thread prepare its own stash
          call hash_table_stash_prepare(swops_stash, max_stash_size, &
               overwrite = .false., overfill_strategy = 'N')

          call timer_clock('pol_emb_uabterm_main_loop',1)
!$OMP DO SCHEDULE(DYNAMIC)
          do didx = 1, npts_dbl

             d1idx = modulo(didx-1,fftbox%total_pt1_dbl) + 1
             d2idx = modulo((didx-d1idx)/fftbox%total_pt1_dbl,fftbox%total_pt2_dbl) + 1
             d3idx = (didx-(d2idx-1)*fftbox%total_pt1_dbl-d1idx) / &
                  (fftbox%total_pt1_dbl*fftbox%total_pt2_dbl) + 1

! jme: workaround for cray compiler bug 843178
#if _CRAYFTN
             curpoint = add_points( &
                  add_points(fftbox_dbl_start_in_cell_vec, (d1idx-1) * a1_dbl), &
                  add_points((d2idx-1) * a2_dbl, (d3idx-1) * a3_dbl) )
#else
             curpoint = fftbox_dbl_start_in_cell_vec + &
                  (d1idx-1) * a1_dbl + (d2idx-1) * a2_dbl + (d3idx-1) * a3_dbl
#endif

             points_used = points_used + 1

             ! jd: This performs read accesses to swops_at_points_ht and
             !     writes to the stash.
             call swri_obtain_swops_at_point(swri, cell, swops_at_point, &
                  curpoint, src_centre, src_atom, sw_or_swpot, swex%quality, &
                  swops_stash)

             ! jd: Accumulate
             centre_sw = 1
             do swex_sw_t_ab = src_offs, src_offs + swex%quality%num_sws_per_centre - 1
                uabterm_in_fftbox_dbl_on_a(d1idx,d2idx,d3idx) = &
                     uabterm_in_fftbox_dbl_on_a(d1idx,d2idx,d3idx) + &
                     swops_at_point(centre_sw) * jqvterms(swex_sw_t_ab)
                centre_sw = centre_sw + 1
             end do

          end do
!$OMP END DO
          call timer_clock('pol_emb_uabterm_main_loop',2)

          ! jd: Merge each thread's stash into the shared hash table now that
          !     it's not being accessed anymore. Destroy stashes.
          call timer_clock('pol_emb_uabterm_stash_commit',1)
!$OMP CRITICAL(SEC_SWOPS_AT_POINTS_HT_ACCESS)
          call hash_table_stash_commit(swops_stash, swri%swops_at_points_ht)
!$OMP END CRITICAL(SEC_SWOPS_AT_POINTS_HT_ACCESS)
          call timer_clock('pol_emb_uabterm_stash_commit',2)

          call hash_table_stash_free(swops_stash)

!$OMP END PARALLEL

          ! jd: Skip to second centre, if any
          src_offs = src_offs + swex%quality%num_sws_per_centre

       end do two_centres_fftbox_branch
       ! ----------------------------------------------------------------------

    end if

  end subroutine pol_emb_uabterm

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pol_emb_tau0term(tau0term_in_ppd, &
       local_a, s_atoms_nl, ngwf_basis, rep, denskern, cell, elements)
    !==========================================================================!
    ! Calculates sum_Bb K^{Aa,Bb} \phi_{Bb} in PPDs on A. The sum goes over all!
    ! overlapping NGWFs, including those on A. The result is calculated for    !
    ! all NGWFs a on A, for all spins.                                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   tau0term_in_ppd (out): The result is returned here. The array must be  !
    !                          allocated by caller in advance.                 !
    !   local_a (in): The local index of atom A.                               !
    !   s_atoms_nl (in): Neighbour list.                                       !
    !   ngwf_basis (in): The NGWF basis in which we operate.                   !
    !   rep (inout): The NGWF_REP container of the NGWF basis. Note that it is !
    !                the responsibility of the caller for all NGWFs Bb to be   !
    !                cached there, so that we avoid comms here.                !
    !   denskern (in): The density kernel.                                     !
    !   cell (in): The simulation cell.                                        !
    !   elements (in).                                                         !
    !--------------------------------------------------------------------------!
    ! Notes:                                                                   !
    !   Mind the requirements listed in the Arguments section.                 !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2017.                              !
    ! Modified to remove pub_par by Joseph Prentice, September 2018            !
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use datatypes, only: FUNCTIONS, data_functions_dealloc
    use function_basis, only: FUNC_BASIS
    use hash_table, only: hash_table_lookup
    use ion, only: ELEMENT
    use neighbour_list, only: NL_NEIGHBOUR_LIST
    use parallel_strategy, only: PARAL_INFO
    use ppd_ops, only: ppd_extract_contents_of_subset, ppd_accumulate, &
         ppd_set_intersection_unpacked
    use remote, only: remote_unpack_ngwf_no_sphere, packed_ngwf_size, &
         remote_ngwf_cache_hts
    use rundat, only: pub_num_spins, pub_dma_use_ri, pub_devel_code
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_get_element, sparse_get_par
    use sw_resolution_of_identity, only: swri_get_handle_to
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_devel_code

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)          :: tau0term_in_ppd(:,:,:) ! ppdpt, is, ngwf_a
    integer, intent(in)                 :: local_a
    type(NL_NEIGHBOUR_LIST), intent(in) :: s_atoms_nl
    type(FUNC_BASIS), intent(in)        :: ngwf_basis
    type(NGWF_REP), intent(in)          :: rep
    type(SPAM3), intent(in)             :: denskern(pub_num_spins)
    type(CELL_INFO), intent(in)         :: cell
    type(ELEMENT), intent(in)           :: elements(:)

    ! jd: Local variables
    integer :: global_a
    integer :: n_ngwfs_on_a, ngwf_a, global_aa_ngwf_idx, local_aa_ngwf_idx
    integer :: local_first_ngwf_idx_of_a, first_ngwf_idx_of_a
    integer :: b_idx, global_b, orig_global_b
    integer :: ngwf_b, first_ngwf_idx_of_b, global_bb_ngwf_idx
    integer :: is
    real(kind=DP) :: dkn_aa_bb_is
    integer :: ppd_indices_on_a(ngwf_basis%max_n_ppds_sphere)
    integer :: ppd_indices_on_b(ngwf_basis%max_n_ppds_sphere)
    integer :: ppd_indices_aa_bb_intersection(ngwf_basis%max_n_ppds_sphere)
    integer :: n_ppds_on_a, n_ppds_on_b
    integer :: n_ppds_aa_bb_intersecion
    type(FUNCTIONS) :: ngwf_bb_in_ppds
    real(kind=DP) :: ngwf_bb_in_ppds_on_a(PPDS_DIMS) ! Phi_Bb in PPDs common to A, B
    real(kind=DP) :: packed_ngwf_bb(packed_ngwf_size)
    integer :: ndata
    integer :: swri_h
    logical :: debug_show_skipped
    character(len=*), parameter :: myself = 'pol_emb_tau0term'

    type(PARAL_INFO), pointer :: par

    ! ------------------------------------------------------------------------

    call timer_clock(myself,1)

    call sparse_get_par(par, denskern(1))

    global_a = par%first_atom_on_proc(pub_my_proc_id) + local_a-1
    n_ngwfs_on_a = ngwf_basis%num_on_atom(global_a)
    first_ngwf_idx_of_a = ngwf_basis%first_on_atom(global_a)
    local_first_ngwf_idx_of_a = ngwf_basis%first_on_atom(global_a) - &
         ngwf_basis%first_on_proc(pub_my_proc_id)+1

    swri_h = swri_get_handle_to(pub_dma_use_ri)
    debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
         pub_devel_code, no_bcast = .true., no_warn=.true.)

    tau0term_in_ppd(:,:,:) = 0D0

    ! jd: Determine the PPDs that belong to A ('s first NGWF)
    n_ppds_on_a = ngwf_basis%spheres(local_first_ngwf_idx_of_a)%n_ppds_sphere
    ppd_indices_on_a(1:n_ppds_on_a) = &
         ngwf_basis%spheres(local_first_ngwf_idx_of_a)%ppd_list(1,1:n_ppds_on_a)
    ppd_indices_on_a(n_ppds_on_a+1:) = -1

    ! -----------------------------------------------------------------------
    ! jd: Loop over B's that are s-neighbours with A                     BBB
    ! -----------------------------------------------------------------------
    loop_B:                                                                &
    do b_idx = s_atoms_nl%first_idx(global_a), s_atoms_nl%last_idx(global_a)
       global_b = s_atoms_nl%neighbours(b_idx)
       orig_global_b = par%orig_atom(global_b)

       ! jd: Ignore atoms not represented in this SWRI
       if(.not. elements(orig_global_b)%in_swri(swri_h)) then
          if(debug_show_skipped) then
             write(stdout,'(a,i0,a,i0,a,a,a)') &
             'Skipping ('//trim(myself)//' orig/SFC: ', orig_global_b, '/', &
             global_b, ' (', trim(elements(orig_global_b)%species_id),')'
          end if
          cycle
       end if

       ! --------------------------------------------------------------------
       ! jd: for all NGWFs b on B                                        bbb
       ! --------------------------------------------------------------------
       loop_ngwf_b:                                                        &
       do ngwf_b = 1, ngwf_basis%num_on_atom(global_b)
          first_ngwf_idx_of_b = ngwf_basis%first_on_atom(global_b)
          global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

          ! jd: Obtain NGWF_Bb and its PPD indices
          call hash_table_lookup(packed_ngwf_bb, ndata, &
               remote_ngwf_cache_hts(rep%ngwf_cache_handle), global_bb_ngwf_idx)
          call utils_assert(ndata /= -1, myself//': NGWF not found in NGWF &
               &cache. This indicates a bug in the logic.')
          call remote_unpack_ngwf_no_sphere(packed_ngwf_bb, &
               ppd_indices_on_b, n_ppds_on_b, ngwf_bb_in_ppds)

          ! --------------------------------------------------------------------
          ! jd: for all NGWFs a on A                                        aaa
          ! --------------------------------------------------------------------
          loop_ngwf_a:                                                        &
          do ngwf_a = 1, ngwf_basis%num_on_atom(global_a)
             global_aa_ngwf_idx = first_ngwf_idx_of_a + ngwf_a - 1

             ! jd: Figure out PPDs common between A and B, once per AB pair.
             if(ngwf_a == 1 .and. ngwf_b == 1) then
                call ppd_set_intersection_unpacked(&
                     ppd_indices_aa_bb_intersection, & ! out
                     n_ppds_aa_bb_intersecion, &       ! out
                     ppd_indices_on_a, n_ppds_on_a, &  ! in
                     ppd_indices_on_b, n_ppds_on_b)    ! in
             end if

             ! -----------------------------------------------------------------
             ! jd: for all spins                                            sss
             ! -----------------------------------------------------------------
             loop_spins:                                                      &
             do is = 1, pub_num_spins

                ! jd: Get K^{Aa,Bb}
                call sparse_get_element(dkn_aa_bb_is, denskern(is), &
                     global_bb_ngwf_idx, global_aa_ngwf_idx)

                ! jd: Generate (relevant fragment of) NGWF Bb in PPDs
                !     common to A and B by extracting
                call ppd_extract_contents_of_subset(ngwf_bb_in_ppds_on_a, & !out
                     ppd_indices_aa_bb_intersection, &     ! wanted subset
                     n_ppds_aa_bb_intersecion, &           ! wanted subset
                     packed_ngwf_bb, ppd_indices_on_b, &   ! mother set
                     n_ppds_on_b, &                        ! mother set
                     cell%n_pts)

                ! jd: Accumulate \tau = K^{Aa,Bb} \phi_{Bb} in PPDs on A.
                call ppd_accumulate(tau0term_in_ppd(:,is,ngwf_a), &
                     ppd_indices_on_a, n_ppds_on_a, &
                     ngwf_bb_in_ppds%d, ppd_indices_aa_bb_intersection, &
                     n_ppds_aa_bb_intersecion, .true., cell%n_pts, dkn_aa_bb_is)

             end do loop_spins

          end do loop_ngwf_a

          call data_functions_dealloc(ngwf_bb_in_ppds)

       end do loop_ngwf_b

    end do loop_B

    call timer_clock(myself,2)

  end subroutine pol_emb_tau0term

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pol_emb_tau1term(tau1term_in_ppd, &
       local_a, s_atoms_nl, ngwf_basis, rep, denskern, cell, &
       jqterms_ht, swri, elements)
    !==========================================================================!
    ! Calculates sum_B U^00_AB sum_b K^{Aa,Bb} \phi_{Bb} in PPDs on A. The sum !
    ! goes over all overlapping NGWFs, including those on A.                   !
    !--------------------------------------------------------------------------!
    ! @docme
    ! Arguments:                                                               !
    !   tau1term_in_ppd: Must be allocated by caller.                          !
    !--------------------------------------------------------------------------!
    ! Notes:                                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2017.                              !
    ! Modified to remove pub_par by Joseph Prentice, September 2018            !
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: REP_SWEX_POL_EMB_DMA_1, REP_SWEX_POL_EMB_DMA_2
    use datatypes, only: FUNCTIONS, data_functions_dealloc
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE, hash_table_lookup_nocount, &
         hash_table_lookup
    use ion, only: ELEMENT
    use neighbour_list, only: NL_NEIGHBOUR_LIST
    use parallel_strategy, only: PARAL_INFO
    use ppd_ops, only: ppd_extract_contents_of_subset, ppd_accumulate, &
         ppd_set_intersection_unpacked, ppd_product
    use remote, only: remote_unpack_ngwf_no_sphere, remote_pack_data_real, &
         packed_ngwf_size, remote_ngwf_cache_hts
    use rundat, only: pub_num_spins, pub_dma_bessel_averaging, pub_dma_use_ri, &
         pub_devel_code
    use simulation_cell, only: CELL_INFO
    use sw_expansion_type, only: SW_EX
    use sw_resolution_of_identity, only: SW_RI, swri_get_handle_to
    use sparse, only: SPAM3, sparse_get_element, sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_sanity_check, utils_assert, utils_alloc_check, &
         utils_dealloc_check, utils_devel_code

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)              :: tau1term_in_ppd(:,:,:) ! ppdpt, is, ngwf_a
    integer, intent(in)                     :: local_a
    type(NL_NEIGHBOUR_LIST), intent(in)     :: s_atoms_nl
    type(FUNC_BASIS), intent(in)            :: ngwf_basis
    type(NGWF_REP), intent(in), target      :: rep
    type(SPAM3), intent(in)                 :: denskern(pub_num_spins)
    type(CELL_INFO), intent(in)             :: cell
    type(HT_HASH_TABLE), intent(in), target :: jqterms_ht
    type(SW_RI), intent(inout)              :: swri ! (caching in uabterm)
    type(ELEMENT), intent(in)               :: elements(:)

    ! jd: Local variables
    integer :: global_a
    integer :: n_ngwfs_on_a, ngwf_a, global_aa_ngwf_idx, local_aa_ngwf_idx
    integer :: local_first_ngwf_idx_of_a, first_ngwf_idx_of_a
    integer :: b_idx, global_b, orig_global_b
    integer :: ngwf_b, first_ngwf_idx_of_b, global_bb_ngwf_idx
    integer :: is
    real(kind=DP) :: dkn_aa_bb_is
    integer :: ppd_indices_on_a(ngwf_basis%max_n_ppds_sphere)
    integer :: ppd_indices_on_b(ngwf_basis%max_n_ppds_sphere)
    integer :: ppd_indices_aa_bb_intersection(ngwf_basis%max_n_ppds_sphere)
    integer :: n_ppds_on_a, n_ppds_on_b
    integer :: n_ppds_aa_bb_intersecion
    type(FUNCTIONS) :: ngwf_bb_in_ppds
    real(kind=DP) :: ngwf_bb_in_ppds_on_a(PPDS_DIMS) ! Phi_Bb in PPDs common to A, B
    real(kind=DP) :: packed_ngwf_bb(packed_ngwf_size)
    integer :: n_elems
    real(kind=DP), allocatable :: jqterm_a(:)
    real(kind=DP), allocatable :: jqterm_b(:)
    real(kind=DP) :: uab00term_in_ppds_on_a(PPDS_DIMS)
    real(kind=DP) :: uab00_ngwf_bb_prod_in_ppds_on_a(PPDS_DIMS)
    integer :: ppd_indices_in_prod(ngwf_basis%max_n_ppds_sphere)
    integer :: n_ppds_in_prod
    real(kind=DP) :: packed_uab00_in_ppds_on_a(packed_ngwf_size)
    type(SW_EX), pointer :: swex
    integer :: bessel_run, n_bessel_runs
    integer :: ierr
    integer :: ndata
    integer :: swri_h
    logical :: debug_show_skipped

    ! jd: Varia
    character(len=*), parameter :: myself = 'pol_emb_tau1term'

    type(PARAL_INFO), pointer :: par

    ! ------------------------------------------------------------------------

    call timer_clock(myself,1)

    call sparse_get_par(par, denskern(1))

    global_a = par%first_atom_on_proc(pub_my_proc_id) + local_a-1
    n_ngwfs_on_a = ngwf_basis%num_on_atom(global_a)
    first_ngwf_idx_of_a = ngwf_basis%first_on_atom(global_a)
    local_first_ngwf_idx_of_a = ngwf_basis%first_on_atom(global_a) - &
         ngwf_basis%first_on_proc(pub_my_proc_id)+1

    swri_h = swri_get_handle_to(pub_dma_use_ri)
    debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
         pub_devel_code, no_bcast = .true., no_warn=.true.)

    tau1term_in_ppd(:,:,:) = 0D0

    ! jd: Determine the PPDs that belong to A ('s first NGWF)
    n_ppds_on_a = ngwf_basis%spheres(local_first_ngwf_idx_of_a)%n_ppds_sphere
    ppd_indices_on_a(1:n_ppds_on_a) = &
         ngwf_basis%spheres(local_first_ngwf_idx_of_a)%ppd_list(1,1:n_ppds_on_a)
    ppd_indices_on_a(n_ppds_on_a+1:) = -1

    if(pub_dma_bessel_averaging) then
       n_bessel_runs = 2
    else
       n_bessel_runs = 1
    end if

    ! -----------------------------------------------------------------
    ! jd: Loop over Bessel averaging runs                          AVG
    ! -----------------------------------------------------------------
    loop_bessel_avg_run:                                           &
    do bessel_run = 1, n_bessel_runs

       if(bessel_run == 1) swex => rep%swexes(REP_SWEX_POL_EMB_DMA_1)
       if(bessel_run == 2) swex => rep%swexes(REP_SWEX_POL_EMB_DMA_2)

       allocate(jqterm_a(swex%quality%num_sws_per_centre),stat=ierr)
       call utils_alloc_check(myself,'jqterm_a',ierr)
       allocate(jqterm_b(swex%quality%num_sws_per_centre),stat=ierr)
       call utils_alloc_check(myself,'jqterm_b',ierr)

       call hash_table_lookup_nocount(jqterm_a, n_elems, jqterms_ht, global_a, &
            bessel_run)
       call utils_assert(n_elems == swex%quality%num_sws_per_centre, &
            myself//': Inconsistent number of elements in jqterm_ht (A)', &
            global_a, n_elems, swex%quality%num_sws_per_centre)

       ! -----------------------------------------------------------------------
       ! jd: Loop over B's that are s-neighbours with A                     BBB
       ! -----------------------------------------------------------------------
       loop_B:                                                                &
       do b_idx = s_atoms_nl%first_idx(global_a), s_atoms_nl%last_idx(global_a)
          global_b = s_atoms_nl%neighbours(b_idx)
          orig_global_b = par%orig_atom(global_b)

          ! jd: Ignore atoms not represented in this SWRI
          if(.not. elements(orig_global_b)%in_swri(swri_h)) then
             if(debug_show_skipped) then
                write(stdout,'(a,i0,a,i0,a,a,a)') &
                'Skipping ('//trim(myself)//' orig/SFC: ', orig_global_b, '/', &
                global_b, ' (', trim(elements(orig_global_b)%species_id),')'
             end if
             cycle
          end if

          ! jd: Calculate uab00, re-determine PPDs on A (even though unneeded)
          uab00term_in_ppds_on_a(:) = 0D0 ! uabterm adds to
          call hash_table_lookup_nocount(jqterm_b, n_elems, jqterms_ht, &
               global_b, bessel_run)
          call utils_assert(n_elems == swex%quality%num_sws_per_centre, &
               myself//': Inconsistent number of elements in jqterm_ht (B)', &
               global_b, n_elems, swex%quality%num_sws_per_centre)
          call pol_emb_uabterm(swri, ngwf_basis, cell, elements, par, &
               local_a, global_b, jqterms_ht, swex, bessel_run, .true., &
               uabterm_in_ppds_on_a = uab00term_in_ppds_on_a, &
               ppd_indices_on_a = ppd_indices_on_a, &
               n_ppds_on_a = n_ppds_on_a)
          call utils_sanity_check(uab00term_in_ppds_on_a, &
               'U_AB00 (ppd) term in '//myself, excessive=1D20)
          call remote_pack_data_real(packed_uab00_in_ppds_on_a, &
               uab00term_in_ppds_on_a, ngwf_basis, local_first_ngwf_idx_of_a)

          ! --------------------------------------------------------------------
          ! jd: for all NGWFs b on B                                        bbb
          ! --------------------------------------------------------------------
          loop_ngwf_b:                                                        &
          do ngwf_b = 1, ngwf_basis%num_on_atom(global_b)
             first_ngwf_idx_of_b = ngwf_basis%first_on_atom(global_b)
             global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

             ! jd: Retrieve packed NGWF Bb from cache
             call hash_table_lookup(packed_ngwf_bb, ndata, &
                  remote_ngwf_cache_hts(rep%ngwf_cache_handle), global_bb_ngwf_idx)
             call utils_assert(ndata /= -1, myself//': NGWF not found in NGWF &
                  &cache. This indicates a bug in the logic.')
             call remote_unpack_ngwf_no_sphere(packed_ngwf_bb, &
                  ppd_indices_on_b, n_ppds_on_b, ngwf_bb_in_ppds)

             if(n_bessel_runs == 2) then
                ngwf_bb_in_ppds%d(:) = ngwf_bb_in_ppds%d * 0.5D0
             end if

             ! --------------------------------------------------------------------
             ! jd: for all NGWFs a on A                                        aaa
             ! --------------------------------------------------------------------
             loop_ngwf_a:                                                        &
             do ngwf_a = 1, ngwf_basis%num_on_atom(global_a)
                global_aa_ngwf_idx = first_ngwf_idx_of_a + ngwf_a - 1

                ! jd: Figure out PPDs common between A and B, once per AB pair.
                if(ngwf_a == 1 .and. ngwf_b == 1) then
                   call ppd_set_intersection_unpacked(&
                        ppd_indices_aa_bb_intersection, & ! out
                        n_ppds_aa_bb_intersecion, &       ! out
                        ppd_indices_on_a, n_ppds_on_a, &  ! in
                        ppd_indices_on_b, n_ppds_on_b)    ! in
                end if

                ! -----------------------------------------------------------------
                ! jd: for all spins                                            sss
                ! -----------------------------------------------------------------
                loop_spins:                                                      &
                do is = 1, pub_num_spins

                   ! jd: Get K^{Aa,Bb}
                   call sparse_get_element(dkn_aa_bb_is, denskern(is), &
                        global_bb_ngwf_idx, global_aa_ngwf_idx)

                   ! jd: Generate (relevant fragment of) NGWF Bb in PPDs
                   !     common to A and B by extracting
                   call ppd_extract_contents_of_subset(ngwf_bb_in_ppds_on_a, & !out
                        ppd_indices_aa_bb_intersection, &     ! wanted subset
                        n_ppds_aa_bb_intersecion, &           ! wanted subset
                        packed_ngwf_bb, ppd_indices_on_b, &   ! mother set
                        n_ppds_on_b, &                        ! mother set
                        cell%n_pts)

                   call ppd_product(packed_uab00_in_ppds_on_a, packed_ngwf_bb, &
                        cell%n_pts, uab00_ngwf_bb_prod_in_ppds_on_a)
                   call utils_sanity_check(uab00_ngwf_bb_prod_in_ppds_on_a, &
                        'U_AB00 * \Phi_Bb term in '//myself,excessive=1D20)

                   ! jd: Accumulate \tau = @@@ K^{Aa,Bb} \phi_{Bb} in PPDs on A.
                   call ppd_accumulate(tau1term_in_ppd(:,is,ngwf_a), &
                        ppd_indices_on_a, n_ppds_on_a, &
                        uab00_ngwf_bb_prod_in_ppds_on_a, &
                        ppd_indices_aa_bb_intersection, &
                        n_ppds_aa_bb_intersecion, .true., cell%n_pts, dkn_aa_bb_is)

                end do loop_spins

             end do loop_ngwf_a

             call data_functions_dealloc(ngwf_bb_in_ppds)

          end do loop_ngwf_b

       end do loop_B

       deallocate(jqterm_b,stat=ierr)
       call utils_dealloc_check(myself,'jqterm_b',ierr)
       deallocate(jqterm_a,stat=ierr)
       call utils_dealloc_check(myself,'jqterm_a',ierr)

    end do loop_bessel_avg_run

    call timer_clock(myself,2)

  end subroutine pol_emb_tau1term

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pol_emb_qsc1term(ngwf_basis, cell, local_a, global_b, &   ! in
       qterms_ht, bessel_run, overlap, tau0term_in_ppds_on_a, &        ! in
       tau1term_in_ppds_on_a, n_sets, sum_corr_terms, lambda, q_DMA, & ! in
       qsc1term_in_ppds_on_a)               ! out/accum depending on bessel_run
    !==========================================================================!
    ! @docme
    ! Work in progress.
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !--------------------------------------------------------------------------!
    ! Notes:                                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2017.                              !
    ! Modified to remove pub_par by Joseph Prentice, September 2018            !
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: SW_V, SW_O, METRIC_OVERLAP
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE, hash_table_lookup_nocount
    use parallel_strategy, only: PARAL_INFO
    use ppd_ops, only: ppd_product
    use remote, only: packed_ngwf_size, remote_pack_data_real
    use rundat, only: pub_num_spins, pub_dma_scale_charge, pub_pol_emb_dbl_grid
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_get_element, sparse_get_par
    use rundat, only: pub_dma_scale_charge, pub_pol_emb_dbl_grid
    use utils, only: utils_assert, utils_abort

    implicit none

    ! jd: Arguments
    type(FUNC_BASIS), intent(in)    :: ngwf_basis
    type(CELL_INFO), intent(in)     :: cell
    integer, intent(in)             :: local_a
    integer, intent(in)             :: global_b
    type(HT_HASH_TABLE), intent(in), target :: qterms_ht
    real(kind=DP), intent(in)       :: tau0term_in_ppds_on_a(PPDS_DIMS)
    real(kind=DP), intent(in)       :: tau1term_in_ppds_on_a(PPDS_DIMS)
    integer, intent(in)             :: bessel_run
    type(SPAM3), intent(in)         :: overlap
    integer, intent(in)             :: n_sets
    real(kind=DP)                   :: q_DMA(2) ! idx: Bessel runs
    real(kind=DP)                   :: sum_corr_terms(2,MAX_N_MULTIPOLE_SETS)
    real(kind=DP)                   :: lambda(2) ! idx: Bessel runs
    real(kind=DP), intent(inout)    :: qsc1term_in_ppds_on_a(PPDS_DIMS, &
         ngwf_basis%max_on_atom, ngwf_basis%max_on_atom, pub_num_spins)

    ! jd: Local variables
    integer       :: first_ngwf_idx_of_a, global_aa_ngwf_idx, ngwf_a, global_a
    integer       :: first_ngwf_idx_of_b, global_bb_ngwf_idx, ngwf_b
    integer       :: n_ngwfs_on_a, n_ngwfs_on_b
    real(kind=DP) :: qterm_a(n_qterm_sph_components) ! Q_A^{lm}
    real(kind=DP) :: QA00
    integer       :: ndata
    integer       :: is
    real(kind=DP) :: prefactor
    real(kind=DP) :: qsc1term_in_ppds_on_a_cur_bessel(PPDS_DIMS, &
         ngwf_basis%max_on_atom, ngwf_basis%max_on_atom, pub_num_spins)
    real(kind=DP) :: S_aa_bb

    character(len=*), parameter :: myself = 'pol_emb_qsc1term'

    type(PARAL_INFO), pointer :: par

    ! -------------------------------------------------------------------------

    call utils_assert(pub_dma_scale_charge, &
         myself//': This subroutine can only be used if dma_scale_charge is T')

    call utils_assert(.not. pub_pol_emb_dbl_grid, &
         myself//': pol_emb_dbl_grid T not supported with DMA charge scaling')

    call sparse_get_par(par, overlap)

    global_a = par%first_atom_on_proc(pub_my_proc_id) + local_a-1
    first_ngwf_idx_of_a = ngwf_basis%first_on_atom(global_a)
    n_ngwfs_on_a = ngwf_basis%num_on_atom(global_a)
    first_ngwf_idx_of_b = ngwf_basis%first_on_atom(global_b)
    n_ngwfs_on_b = ngwf_basis%num_on_atom(global_b)

    call hash_table_lookup_nocount(qterm_a, ndata, qterms_ht, global_a)
    call utils_assert(ndata /= -1, myself//': qterm_a not there', global_a)
    QA00 = qterm_a(1)

    prefactor = QA00 - sum(sum_corr_terms(bessel_run,1:n_sets)) / q_DMA(bessel_run)

    ! --------------------------------------------------------------------------
    ! jd: for all spins                                                     sss
    ! --------------------------------------------------------------------------
    loop_spins:                                                               &
    do is = 1, pub_num_spins

       ! -----------------------------------------------------------------------
       ! jd: for all NGWFs a on A                                           aaa
       ! -----------------------------------------------------------------------
       loop_ngwf_a:                                                           &
       do ngwf_a = 1, n_ngwfs_on_a
          global_aa_ngwf_idx = first_ngwf_idx_of_a + ngwf_a - 1

          ! -----------------------------------------------------------------
          ! jd: for all NGWFs b on B                                     bbb
          ! -----------------------------------------------------------------
          loop_ngwf_b:                                                     &
          do ngwf_b = 1, n_ngwfs_on_b
             global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

             call sparse_get_element(S_aa_bb, overlap, &
                  global_bb_ngwf_idx, global_aa_ngwf_idx)

             ! jd: Look up d coefficients for this Aa Bb
!             call hash_table_lookup_nocount(d_alpha_beta, &          ! out
!                  num_coeffs, &                                      ! out
!                  dcoeffs_ht(bessel_run), &                          ! in
!                  global_aa_ngwf_idx, global_bb_ngwf_idx, &          ! in
!                  merge(SW_O,SW_V,pub_dma_metric == METRIC_OVERLAP)) ! in

!             call utils_assert(num_coeffs == &
!                  (pub_dma_max_l+1)**2, myself//& !@@min_l
!                  ': Internal problem with dcoeff lookup.', &
!                  global_bb_ngwf_idx, global_cc_ngwf_idx, coeffs_kind, &
!                  num_coeffs)

!             qsc1term_in_ppds_on_a_cur_bessel(:,ngwf_b,ngwf_a,is) =
!                  prefactor * &
!                  (d_alpha_beta * tau0term_in_ppds_on_a(:) + &
!                  2D0 * (S_aa_bb - lambda(bessel_run) * d_alpha_beta) * &
!                  tau1term_in_ppds_on_a(:))

          end do loop_ngwf_b

       end do loop_ngwf_a

    end do loop_spins

    if(bessel_run == 1) then
!       qsc1term_in_ppds_on_a(:,:,:,:) = &
!            qsc1term_in_ppds_on_a_cur_bessel(:,:,:,:)
    else
!       qsc1term_in_ppds_on_a(:) = &
!            qsc1term_in_ppds_on_a(:,:,:,:) + &
!            qsc1term_in_ppds_on_a_cur_bessel(:,:,:,:)
    end if

  end subroutine pol_emb_qsc1term

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pol_emb_qsc2term(&
       ngwf_basis, cell, par, local_a, global_b, qterms_ht, ngwfs_ht, & ! in
       uab00term_in_ppds_on_a, bessel_run, n_sets, sum_corr_terms, &    ! in
       lambda, q_DMA, &                                                 ! in
       qsc2term_in_ppds_on_a)               ! out/accum depending on bessel_run
    !==========================================================================!
    ! @docme
    ! Work in progress.
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !--------------------------------------------------------------------------!
    ! Notes:                                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2017.                              !
    ! Modified to remove pub_par by Joseph Prentice, September 2018            !
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE, hash_table_lookup_nocount
    use parallel_strategy, only: PARAL_INFO
    use ppd_ops, only: ppd_product
    use remote, only: packed_ngwf_size, remote_pack_data_real
    use simulation_cell, only: CELL_INFO
    use rundat, only: pub_dma_scale_charge, pub_pol_emb_dbl_grid
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(PARAL_INFO), intent(in)    :: par
    type(FUNC_BASIS), intent(in)    :: ngwf_basis
    type(CELL_INFO), intent(in)     :: cell
    integer, intent(in)             :: local_a
    integer, intent(in)             :: global_b
    type(HT_HASH_TABLE), intent(in), target :: qterms_ht
    type(HT_HASH_TABLE), intent(in), target :: ngwfs_ht
    real(kind=DP), intent(in)       :: uab00term_in_ppds_on_a(PPDS_DIMS)
    integer, intent(in)             :: bessel_run
    integer, intent(in)             :: n_sets
    real(kind=DP)                   :: q_DMA(2) ! idx: Bessel runs
    real(kind=DP)                   :: sum_corr_terms(2,MAX_N_MULTIPOLE_SETS)
    real(kind=DP)                   :: lambda(2) ! idx: Bessel runs
    real(kind=DP), intent(inout)    :: qsc2term_in_ppds_on_a(PPDS_DIMS)

    ! jd: Local variables
    integer       :: global_a
    real(kind=DP) :: qterm_a(n_qterm_sph_components) ! Q_A^{lm}
    real(kind=DP) :: QA00
    integer       :: ndata
    real(kind=DP) :: bracket_term(PPDS_DIMS)
    real(kind=DP) :: packed_ngwf_bb(packed_ngwf_size)
    real(kind=DP) :: packed_bracket_term(packed_ngwf_size)
    real(kind=DP) :: qsc2term_in_ppds_on_a_cur_bessel(PPDS_DIMS)

    character(len=*), parameter :: myself = 'pol_emb_qsc2term'

    ! -------------------------------------------------------------------------

    call utils_assert(pub_dma_scale_charge, &
         myself//': This subroutine can only be used if dma_scale_charge is T')

    call utils_assert(.not. pub_pol_emb_dbl_grid, &
         myself//': pol_emb_dbl_grid T not supported with DMA charge scaling')

    global_a = par%first_atom_on_proc(pub_my_proc_id) + local_a-1

    call hash_table_lookup_nocount(qterm_a, ndata, qterms_ht, global_a)
    call utils_assert(ndata /= -1, myself//': qterm_a not there', global_a)
    QA00 = qterm_a(1)

    bracket_term(:) = uab00term_in_ppds_on_a(:) * &
          ((lambda(bessel_run) - 1D0) * QA00 - &
          lambda(bessel_run) * &
          sum(sum_corr_terms(bessel_run,1:n_sets))/q_DMA(bessel_run)) + &
          sum(sum_corr_terms(bessel_run,1:n_sets))/(2D0*q_DMA(bessel_run))
    call remote_pack_data_real(packed_bracket_term, bracket_term, ngwf_basis, 1)

    call hash_table_lookup_nocount(packed_ngwf_bb, ndata, ngwfs_ht, global_b)
    call utils_assert(ndata == packed_ngwf_size, myself//': NGWF_Bb not there &
         &or length mismatch', global_b, ndata, packed_ngwf_size)

    call ppd_product(packed_bracket_term, packed_ngwf_bb, cell%n_pts, &
         qsc2term_in_ppds_on_a_cur_bessel)

    if(bessel_run == 1) then
       qsc2term_in_ppds_on_a(:) = qsc2term_in_ppds_on_a_cur_bessel(:)
    else
       qsc2term_in_ppds_on_a(:) = &
            qsc2term_in_ppds_on_a(:) + qsc2term_in_ppds_on_a_cur_bessel(:)
    end if

  end subroutine pol_emb_qsc2term

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pol_emb_qterms_and_sum_corr_terms(&
       ngwf_basis, qm_multipoles, mm_multipoles, elements, par, & ! in
       sum_corr_terms, lambda, q_DMA, qterms)
    !==========================================================================!
    ! Calculates Q terms (see pol_emb_qterm()) for all atoms local to this     !
    ! MPI rank. Calculates the charge-scaling correction auxiliary term        !
    ! s_corr (reduced over all atoms in the system). Calculates the charge-    !
    ! scaling factor \lambda and the total DMA charge q_DMA.                   !
    !                                                                          !
    ! This is calculated for *one* MM multipole set (mm_multipoles), and one   !
    ! Bessel run (corresponding to the passed qm_multipoles). Suitable summa-  !
    ! tion or averaging must be done externally.                               !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   ngwf_basis (in): The NGWF basis in which we operate.                   !
    !   qm_multipoles (in): Multipoles obtained from DMA for which results are !
    !                       calculated. Only needed for their %{min,max}_l and !
    !                       to get unscaled monopole for charge scaling.       !
    !   mm_multipoles (in): Multipoles describing the MM system.               !
    !   elements (in): Needed for QM site positions and polarisabilities.      !
    !   sum_corr_terms (out): s_corr_term, reduced over all atoms in the system!
    !                         is returned here.                                !
    !   lambda (out): Charge-scaling factor is returned here, if charge-scaling!
    !                 is in effect. Otherwise returns 1.0.                     !
    !   q_DMA (out): Total DMA charge is returned here.                        !
    !   qterms (out, opt): If provided, Q terms for atoms local to this MPI    !
    !                      rank are returned here.                             !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    ! - Be careful to call it from all MPI ranks, as there's a reduce in here. !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2017.                              !
    ! Modified to remove pub_par by Joseph Prentice, September 2018            !
    !==========================================================================!

    use comms, only: pub_my_proc_id, comms_reduce
    use constants, only: SAFE_DIV_EPS
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_dma_use_ri, pub_dma_scale_charge, pub_devel_code
    use utils, only: utils_assert, utils_devel_code
    use sw_resolution_of_identity, only: swri_get_handle_to

    implicit none

    ! jd: Arguments
    type(PARAL_INFO), intent(in)              :: par
    type(FUNC_BASIS), intent(in)              :: ngwf_basis
    type(SPHERICAL_MULTIPOLE_SET), intent(in) :: qm_multipoles
    type(SPHERICAL_MULTIPOLE_SET), intent(in) :: mm_multipoles
    type(ELEMENT), intent(in)                 :: elements(:)
    real(kind=DP), intent(out)                :: sum_corr_terms
    real(kind=DP), intent(out)                :: lambda
    real(kind=DP), intent(out)                :: q_DMA
    real(kind=DP), intent(out), optional      :: qterms(:,:) ! 2nd: loc atom idx

    ! jd: Local variables
    integer :: local_b, global_b, orig_global_b
    integer :: mp_begin_loc
    real(kind=DP) :: M00_B ! unscaled monopole on atom B
    integer :: swri_h
    logical :: debug_show_skipped
    character(len=*), parameter :: myself = 'pol_emb_qterms_and_sum_corr_terms'

    ! -------------------------------------------------------------------------

    swri_h = swri_get_handle_to(pub_dma_use_ri)

    debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
         pub_devel_code, no_bcast = .true., no_warn = .true.)

    q_DMA = 0D0
    sum_corr_terms = 0D0

    if(present(qterms)) then
       qterms(:,:) = 0D0 ! jd: since some B's may be cycled over
    end if

    ! -----------------------------------------------------------------------
    ! jd: Loop over B's -- all atoms on this core                        BBB
    ! -----------------------------------------------------------------------
    loop_B0:                                                               &
    do local_b = 1, par%num_atoms_on_proc(pub_my_proc_id)
       global_b = par%first_atom_on_proc(pub_my_proc_id) + local_b-1
       orig_global_b = par%orig_atom(global_b)

       ! jd: Ignore atoms not represented in this SWRI
       if(.not. elements(orig_global_b)%in_swri(swri_h)) then
          if(debug_show_skipped) then
             write(stdout,'(a,i0,a,i0,a,a,a)') &
             'Skipping ('//trim(myself)//' orig/SFC: ', orig_global_b, '/', &
             global_b, ' (', trim(elements(orig_global_b)%species_id),')'
          end if
          cycle
       end if

       mp_begin_loc = (local_b-1) * qm_multipoles%n_multipoles_per_site+1
       M00_B = qm_multipoles%mpoles(mp_begin_loc)
       ! jd: Recover unscaled monopole, accumulate in q_DMA
       if(pub_dma_scale_charge) then
          lambda = qm_multipoles%dma_monopole_scaling_factor
          call utils_assert(lambda > SAFE_DIV_EPS, &
               myself//': Bad dma_monopole_scaling_factor', lambda)
       else
          lambda = 1D0
       end if
       M00_B = M00_B / lambda
       q_DMA = q_DMA + M00_B

       if(present(qterms)) then
          call pol_emb_qterm(mm_multipoles, elements, par, global_b, &
               qm_multipoles%min_l, qm_multipoles%max_l, &
               M00_B, qterms(:,local_b),  sum_corr_terms)
       else
          call pol_emb_qterm(mm_multipoles, elements, par, global_b, &
               qm_multipoles%min_l, qm_multipoles%max_l, &
               M00_B = M00_B, sum_corr_terms = sum_corr_terms)
       end if

    end do loop_B0

    ! jd: Reduce previously accumulated values (QM sites are distributed).
    call comms_reduce('SUM',sum_corr_terms)
    call comms_reduce('SUM',q_DMA)
!@@@ not valid now that we can have qterms for polarisation multipoles
!    call utils_assert(q_DMA > SAFE_DIV_EPS, &
!         myself//': Total DMA charge unreasonable', q_DMA)

  end subroutine pol_emb_qterms_and_sum_corr_terms

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pol_emb_repulsive_mm_pot(repulsive_mm_pot, &
       mm_multipoles, grid, cell, only_one_mm_atom)
    !=========================================================================!
    ! @docme
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   only_one_mm_atom (in,opt) -- if specified, the MM rep pot originating !
    !                                from just one (this) MM atom is calc'ed, !
    !                                not from all mm_multipoles.              !
    ! ----------------------------------------------------------------------- !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic in March 2017.                                !
    ! Rewritten by Jacek Dziedzic in September 2017.                          !
    ! 'only_one_mm_atom' added by Jacek Dziedzic in September 2018.           !
    !=========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_real_pt
    use comms, only: pub_on_root, pub_my_proc_id
    use constants, only: stdout, SAFE_DIV_EPS, PI
    use geometry, only: POINT, magnitude, OPERATOR(-), OPERATOR(+), OPERATOR(*), operator(.dot.)
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET, multipole_1_sph_real_to_cart
    use rundat, only: pub_pol_emb_repulsive_mm_pot_a, &
         pub_pol_emb_repulsive_mm_pot_b, pub_pol_emb_repulsive_mm_pot_c, &
         pub_pol_emb_repulsive_mm_pot_alpha, pub_pol_emb_repulsive_mm_pot_beta,&
         pub_pol_emb_repulsive_mm_pot_r0, pub_pol_emb_repulsive_mm_pot_cutoff, &
         pub_pol_emb_repulsive_mm_pot_write, pub_pol_emb_repulsive_mm_pot_verbose
!$  use rundat, only: pub_threads_max
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_erf, utils_real_to_str, utils_int_to_str, &
         utils_abort, utils_assert
    use visual, only: visual_scalarfield

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)    :: repulsive_mm_pot(:,:,:)
    type(SPHERICAL_MULTIPOLE_SET) :: mm_multipoles
    type(GRID_INFO), intent(in)   :: grid
    type(CELL_INFO), intent(in)   :: cell
    integer, intent(in), optional :: only_one_mm_atom

    ! jd: Local variables
    integer       :: islab12, i1, i2, i3, ipt, L
    integer       :: n_mm_at
    integer       :: mm_atom_first, mm_atom_last
    type(POINT)   :: r, R_L, dvec
    real(kind=DP) :: r_pt(3)
    real(kind=DP) :: dvec_magn, accum
    real(kind=DP) :: eta
    real(kind=DP) :: zeta, magn
    character(len=256) :: output_filename_suffix
    integer, save :: pot_output_iter = 1
    character(len=*), parameter :: myself = 'pol_emb_repulsive_mm_pot'

    !------------------------------------------------------------------------

    call timer_clock(myself,1)

    ! jd: Take care of padding between pt and ld
    repulsive_mm_pot = 0.0_DP

    n_mm_at = mm_multipoles%n_sites

    if(present(only_one_mm_atom)) then
       mm_atom_first = only_one_mm_atom
       mm_atom_last = only_one_mm_atom
    else
       mm_atom_first = 1
       mm_atom_last = n_mm_at
    end if

    call utils_assert(mm_atom_first >= 1, &
         myself//': Illegal mm_atom_first: ', mm_atom_first)
    call utils_assert(mm_atom_last <= mm_multipoles%n_sites, &
         myself//': Illegal mm_atom_last: ', mm_atom_last)

    ! jd: Generate MM repulsive potential on every point of fine grid
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt, i1, i2, islab12, i3, r, r_pt, accum, L, R_L, dvec_magn, &
!$OMP      dvec, eta, zeta, magn) &
!$OMP SHARED(pub_pol_emb_repulsive_mm_pot_cutoff, n_mm_at, grid, mm_multipoles, &
!$OMP      repulsive_mm_pot, pub_my_proc_id, pub_threads_max, mm_atom_first, &
!$OMP      mm_atom_last)
    do ipt = 1, grid%num_my_slabs12 * grid%n1 * grid%n2
       i1 = modulo(ipt-1,grid%n1) + 1
       i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
       islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1
       i3 = grid%first_slab12(pub_my_proc_id) + islab12 - 1

       call cell_grid_real_pt(r_pt,i1,i2,islab12,grid)
       r%x=r_pt(1); r%y=r_pt(2); r%z=r_pt(3)

       accum = 0.0_DP
       do L = mm_atom_first, mm_atom_last

          ! jd: Skip dummy QM multipoles masquerading as MM multipoles
          if(mm_multipoles%max_l_mask(L) == -1) cycle

          R_L = mm_multipoles%centres(L)
          dvec = r-R_L ! @MIC
          dvec_magn = magnitude(dvec)

#if 0
          ! jd: This is for the old scheme that tried using MM-side classical
          !     vdW params to determine suitable zeta and magnitude. Leaving
          !     snippets of the code in case this needs to be resurrected someday.
          vdw_rmin = mm_multipoles%userdata(1,L)
          vdw_eps = mm_multipoles%userdata(2,L)
          if(vdw_rmin > SAFE_DIV_EPS) then
             denom = a * (beta/vdw_rmin)*(beta/vdw_rmin) + b*(beta/vdw_rmin) + c
          else
             call utils_abort(myself//': vdw_rmin negative or extremely close &
                  &to zero, atom #'//trim(utils_int_to_str(L)))
          end if
          if(abs(denom) < SAFE_DIV_EPS) then
             call utils_abort(myself//': denom is zero or extremely close &
                  &to zero')
          end if

          ! jd: If R0 is -1, override eta to 1
          if(r0 == -1D0) then
             eta = 1D0
          else
             eta = num / denom
          end if
#endif

          ! --------------------------------------------------------------
          ! --- New scheme for MM rep params: user-provided via %block ---
          ! --------------------------------------------------------------
          ! jd: Look up magn and zeta based on MM index. This needs to be
          !     efficient, as we're inside a N_MMat x N_QMfinegridpts loop.
          magn = mm_multipoles%userdata(1,L)
          zeta = mm_multipoles%userdata(2,L)
          eta = 1D0 ! unused in the new scheme

          ! jd: Only calculate steric potential within a cutoff around each
          !     atom, to improve efficiency
          if(dvec_magn <= pub_pol_emb_repulsive_mm_pot_cutoff) then
             accum = accum + magn / PI * eta * &
                  exp(-2D0 * zeta * dvec_magn)
          end if

       end do ! over MM atoms L

       repulsive_mm_pot(i1,i2,islab12) = accum

    end do ! ipt
!$OMP END PARALLEL DO

    if(pub_pol_emb_repulsive_mm_pot_write) then
       if(pub_on_root) then
          write(stdout,'(a)') ''
       end if
       ! jd: Output filename includes extra details if verbose MM pots used
       if(pub_pol_emb_repulsive_mm_pot_verbose) then
          write(output_filename_suffix,'(a,i0,a,i0,a,i0)') &
               '_repulsive_MM_pot_'//trim(mm_multipoles%set_name)//'_atoms_', &
               mm_atom_first, '_to_', mm_atom_last, '_', pot_output_iter
       else
          write(output_filename_suffix,'(a,i0)') &
               '_repulsive_MM_pot_'//trim(mm_multipoles%set_name)//'_', &
               pot_output_iter
       end if
       call visual_scalarfield(repulsive_mm_pot, grid, cell, &
            'repulsive MM potential, atomic units', &
            output_filename_suffix)
       pot_output_iter = pot_output_iter + 1
    end if

    call timer_clock(myself,2)

  end subroutine pol_emb_repulsive_mm_pot

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pol_emb_calc_pmatrix(cur_polemb, &                         ! in/out
       overlap, qm_multipoles, mm_multipoles, dcoeffs_ht, ngwf_basis, & ! in
       elements, s_atoms_nl, scal_factor, min_l, max_l, swex_h, hack_mask, one_atom_only)         ! in
    !==============================================================!
    !--------------------------------------------------------------!
    ! @docme
    ! Caller is responsible for providing a created cur_polemb on
    ! entry. It is zeroed here.                                    !
    ! Only one bessel run here.
    ! @updatedoc: swex_h -- only for HT key
    ! @qm_multipoles only needed for Qterms, but only for their
    ! {min,max}_l and to get unscaled monopole for charge scaling.
    ! @overlap only needed for charge scaling
    !--------------------------------------------------------------!
    ! Extracted from pol_emb_energy_and_lnv_gradient by Jacek      !
    ! Dziedzic on 2017/05/22.                                      !
    ! Modified to remove pub_par by Joseph Prentice, September 2018!
    !==============================================================!

    use comms, only: pub_my_proc_id, pub_on_root
    use constants, only: SW_V, SW_O, SW_D, METRIC_OVERLAP
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(-)
    use hash_table, only: HT_HASH_TABLE, hash_table_lookup_nocount
    use ion, only: ELEMENT
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET
    use neighbour_list, only: NL_NEIGHBOUR_LIST
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_dma_scale_charge, pub_dma_metric, pub_devel_code, &
         pub_dma_use_ri
    use sparse, only: SPAM3, sparse_put_element, &
         sparse_get_element, sparse_scale, sparse_create, sparse_transpose, &
         sparse_axpy, sparse_destroy, sparse_show_matrix, sparse_get_par
    use sw_resolution_of_identity, only: swri_get_handle_to
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_devel_code, &
         utils_assert, utils_int_to_str

    implicit none

    ! jd: Arguments
    type(SPAM3), intent(inout)                :: cur_polemb
    type(SPAM3), intent(in)                   :: overlap
    type(SPHERICAL_MULTIPOLE_SET), intent(in) :: qm_multipoles
    type(SPHERICAL_MULTIPOLE_SET), intent(in) :: mm_multipoles
    type(HT_HASH_TABLE), intent(in), target   :: dcoeffs_ht
    type(FUNC_BASIS), intent(in)              :: ngwf_basis
    type(ELEMENT), intent(in)                 :: elements(:)
    type(NL_NEIGHBOUR_LIST), intent(in)       :: s_atoms_nl
    real(kind=DP), intent(in)                 :: scal_factor
    integer, intent(in)                       :: min_l
    integer, intent(in)                       :: max_l
    integer, intent(in)                       :: swex_h
    logical, intent(in), optional             :: hack_mask
    integer, intent(in), optional             :: one_atom_only

    ! jd: Local variables
    !     'L' are polarisable sites, 'B' and 'C' are QM atoms.
    integer           :: site_l     ! site index
    type(POINT)       :: vec_l      ! position vector
    real(kind=DP)     :: alpha_l    ! polarisability

    ! jd: DMA coefficients for current MM atom, for all l
    real(kind=DP)     :: d_beta_gamma(n_qterm_sph_components)

    ! P matrix auxiliaries
    type(SPAM3)       :: cur_tpolemb  ! Transpose of current P, for symmetrising
    real(kind=DP)     :: P_beta_gamma      ! P matrix element
    real(kind=DP)     :: P_beta_gamma_corr ! P matrix element correction term

    ! Varia
    integer           :: swri_h
    real(kind=DP)     :: S_gamma_beta ! Overlap matrix element
    real(kind=DP)     :: centre_factor! scaling factor
    integer           :: num_coeffs
    integer           :: coeffs_kind
    real(kind=DP)     :: lambda ! charge scaling factor
    integer           :: ierr
    logical           :: debug_show_skipped

    ! jd: Qterm and charge scaling
    real(kind=DP), allocatable :: qterms(:,:)
    real(kind=DP)     :: sum_corr_terms       ! Total charge-scaling correction
    real(kind=DP)     :: M00_B                ! Unscaled monopole on site B
    real(kind=DP)     :: q_DMA                ! Total unscaled DMA monopole

    ! Atom and NGWF counters
    integer           :: local_b, global_b, first_ngwf_idx_of_b, orig_global_b
    integer           :: c_idx, global_c, first_ngwf_idx_of_c, orig_global_c
    integer           :: ngwf_b, ngwf_c
    integer           :: global_bb_ngwf_idx, global_cc_ngwf_idx

    character(len=*), parameter :: myself = 'pol_emb_calc_pmatrix'

    type(PARAL_INFO), pointer :: par

    ! -------------------------------------------------------------------------

!    write(*,*) '@PMATRIX ', min_l, max_l

    swri_h = swri_get_handle_to(pub_dma_use_ri)

    coeffs_kind = merge(SW_O,SW_V,pub_dma_metric == METRIC_OVERLAP)

    debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
         pub_devel_code,no_warn=.true.)

    call sparse_get_par(par, overlap)

    allocate(qterms(n_qterm_sph_components, &
         par%num_atoms_on_proc(pub_my_proc_id)), stat=ierr)
    call utils_alloc_check(myself,'qterms',ierr)

    call sparse_scale(cur_polemb, 0D0)

    d_beta_gamma = 0D0 ! Later we only fill entries up to DMA's lmax

    ! jd: Fill qterms for all local atoms B, calculate sum_corr_term,
    !     charge-scaling factor, total DMA monopole (scaled) -- all for the
    !     current Bessel run, current MM set.
!@    write(*,*) '@@@QM_MULTIPOLES: '
!@    call multipole_set_info(qm_multipoles, stdout)
!@    write(*,*) '@@@MM_MULTIPOLES: '
!@    call multipole_set_info(mm_multipoles, stdout)
!@@    mmm_multipoles = mm_multipoles
!    unit = utils_unit()
!    open(unit, file='mpoles')
!    read(unit,*) mmm_multipoles%mpoles
!    write(*,*) '@@@MM_MULTIPOLES_OVERWRITTEN: '
!    call multipole_set_info(mmm_multipoles, stdout)
!    close(unit)

    call pol_emb_qterms_and_sum_corr_terms(ngwf_basis, qm_multipoles, &
         mm_multipoles, elements, par, sum_corr_terms, lambda, q_DMA, qterms)

    ! -------------------------------------------------------------------------
    ! jd: Loop over B's -- all atoms on this core                        BBB
    ! -------------------------------------------------------------------------
    loop_B:                                                                   &
    do local_b = 1, par%num_atoms_on_proc(pub_my_proc_id)
       global_b = par%first_atom_on_proc(pub_my_proc_id) + local_b-1
       orig_global_b = par%orig_atom(global_b)
       first_ngwf_idx_of_b = ngwf_basis%first_on_atom(global_b)

       if(present(one_atom_only)) then
          if(one_atom_only /= orig_global_b) then
             write(*,*) '@@@ IGNORING ATOM ', orig_global_b, ', ONLY CALCULATING P for ', one_atom_only
             cycle
          end if
       end if

       ! jd: Ignore atoms not represented in this SWRI
       if(.not. elements(orig_global_b)%in_swri(swri_h)) then
          if(debug_show_skipped) then
             write(stdout,'(a,i0,a,i0,a,a,a)') &
             'Skipping B ('//trim(myself)//' orig/SFC: ', orig_global_b, '/', &
             global_b, ' (', trim(elements(orig_global_b)%species_id),')'
          end if
          cycle
       end if

       ! ----------------------------------------------------------------------
       ! jd: Loop over C's that are s-neighbours with B                  CCC
       ! ----------------------------------------------------------------------
       loop_C:                                                                &
       do c_idx = s_atoms_nl%first_idx(global_b), s_atoms_nl%last_idx(global_b)
          global_c = s_atoms_nl%neighbours(c_idx)
          orig_global_c = par%orig_atom(global_c)
          first_ngwf_idx_of_c = ngwf_basis%first_on_atom(global_c)

          ! jd: No skipping of C atoms that are not in SWRI.
          !     C's outside of SWRI do contribute to non-diagonal elements
          !     of P. Only C-C (diagonal) elements are zero.

          ! jd: 1/2 from the fact that we later add the transpose of P
          !     rather than averaging with it.
          !     Another 1/2 from ???
          if(global_b == global_c) then
             centre_factor = 0.25_DP
          else
             centre_factor = 0.25_DP
          end if

          ! -------------------------------------------------------------------
          ! jd: for all b on B                                           bbb
          ! -------------------------------------------------------------------
          loop_ngwf_b:                                                        &
          do ngwf_b = 1, ngwf_basis%num_on_atom(global_b)
             global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

             ! ----------------------------------------------------------------
             ! jd: for all c on C                                        ccc
             ! ----------------------------------------------------------------
             loop_ngwf_c:                                                     &
             do ngwf_c = 1, ngwf_basis%num_on_atom(global_c)
                global_cc_ngwf_idx = first_ngwf_idx_of_c + ngwf_c - 1

                ! jd: Look up d coefficients for this Bb Cc
                call hash_table_lookup_nocount(d_beta_gamma, num_coeffs, & ! out
                     dcoeffs_ht, global_bb_ngwf_idx, global_cc_ngwf_idx, & ! in
                     coeffs_kind + SW_D, swex_h)                           ! in

                call utils_assert(num_coeffs == &
                     (1+max_l-min_l) * (1+max_l+min_l), &
                     myself//': Internal problem with dcoeff lookup for '//&
                     trim(utils_int_to_str(global_bb_ngwf_idx))//'-'//&
                     trim(utils_int_to_str(global_cc_ngwf_idx)), &
                     coeffs_kind + SW_D, swex_h, num_coeffs, &
                     (1+max_l-min_l) * (1+max_l+min_l))

                if(present(hack_mask)) then
                   if(hack_mask) then
                      write(*,*) '@@@HACK-MASK in DCOEFFS @@@'
                      d_beta_gamma(1) = 0D0
                      d_beta_gamma(5:9) = 0D0
                   end if
                end if

                if(pub_dma_scale_charge) then
                   d_beta_gamma(1) = d_beta_gamma(1) * lambda
                end if

                ! *** Matrix element of the polemb (P) matrix ***
                P_beta_gamma = - scal_factor * centre_factor * &
                     sum(d_beta_gamma(:) * qterms(:,local_b))

!                write(*,*) '@D: ', d_beta_gamma
!                write(*,*) '@Q: ', qterms(:,local_b)

                ! ***********************************************************
                ! -*- Calculate charge scaling correction terms.          -*-
                ! -*- These come from DMA's charges only, which interact  -*-
                ! -*- with MM charges, dipoles and quadrupoles.           -*-
                ! ***********************************************************
                if(pub_dma_scale_charge) then
                   ! jd: Retrieve S_gamma_beta. Assume overlap symmetry,
                   !     as S_beta_gamma is easier to get, being local.
                   call sparse_get_element(S_gamma_beta, overlap, &
                        global_cc_ngwf_idx, global_bb_ngwf_idx)

                   P_beta_gamma_corr = &
                        - scal_factor * centre_factor * sum_corr_terms * &
                        (2.0_DP*S_gamma_beta - d_beta_gamma(1)) / q_DMA

                   P_beta_gamma = P_beta_gamma + P_beta_gamma_corr
                end if

                ! jd: Store calculated final P matrix element
                call sparse_put_element(P_beta_gamma, cur_polemb, &
                     global_cc_ngwf_idx, global_bb_ngwf_idx)

             end do loop_ngwf_c
          end do loop_ngwf_b
       end do loop_C
    end do loop_B

    deallocate(qterms,stat=ierr)
    call utils_dealloc_check(myself,'qterms',ierr)

    ! @@@comment is untrue
    ! jd: Take care of the lower triangle, fill by symmetry. Elements
    !     have been scaled in advance, so just sum with the tranpose.
    call sparse_create(cur_tpolemb, cur_polemb)
    call sparse_transpose(cur_tpolemb, cur_polemb)
    call sparse_axpy(cur_polemb, cur_tpolemb, 1.0_DP)
    call sparse_destroy(cur_tpolemb)

!    if(pub_on_root) then
!       print *, ""
!       print *, "-PMAT:-------------------------------"
!    end if
!    call sparse_show_matrix(cur_polemb)

  end subroutine pol_emb_calc_pmatrix

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pol_emb_ngwf_grad_contribution(swexes, &
       ngwfs_on_grid, packed_ngwf_cache_hts, &
       denskern, ngwf_basis, mdl, tc_denskern, &              ! input
       cov_grad, contra_grad, &                               ! in/out (adds to)
       precond_func_recip, ngwf_cache_handle, prefactor)      ! in
    !==========================================================================!
    ! This subroutine calculates the additional contribution to the NGWF gra-  !
    ! dient from polarisable embedding treated in the QM* representation.      !
    ! Embeddings handled through the 'full QM' scheme lack this contribution,  !
    ! and are handled through a usual contribution to the LHXC potential.      !
    !                                                                          !
    ! Contributions can arise from either or both sets of MM multipoles (perma-!
    ! nent or/and induced).                                                    !
    !                                                                          !
    ! Reciprocal-space preconditioning is applied, if in effect.               !
    !                                                                          !
    ! Depending on swx_dbl_grid, the gradient is calculated on the coarse grid !
    ! (default), or, painfully, on the double grid. In the latter case, NGWFs  !
    ! are interpolated to the double grid, U_AB is calculated in an entire     !
    ! FFTbox around atom A, multiplication happens on the double grid, and the !
    ! result is coarsened back to the coarse grid. This is experimental        !
    ! functionality and far from being clean. The default coarse grid version  !
    ! is clean.                                                                !
    !                                                                          !
    ! DMA charge scaling leads to an additional term in the gradient and is    !
    ! not yet implemented.                                                     !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! @updatedoc
    !   rep (in/out): The NGWF_REP that contributes to the gradient. In/out    !
    !                 because NGWF caching happens behind the scenes.          !
    !   denskern, ngwf_basis, mdl (in): The usual.                             !
    !   tc_denskern (in): Tensorially-correct denskern.                        !
    !   cov_grad, cov_grad (in/out): Results are accumulated here.             !
    !   precond_func_recip (in): Needed for recip-space precionditioning.      !
    !   @factor
    !--------------------------------------------------------------------------!
    ! Notes:                                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October-November 2016.                      !
    ! Modified to remove pub_par by Joseph Prentice, September 2018            !
    !==========================================================================!

    use basis, only: basis_clean_function, basis_copy_function_to_box, &
         basis_extract_function_from_box
    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_send, &
         comms_barrier, pub_on_root, comms_wait
    use constants, only: PI, stdout, VERBOSE
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_functions_alloc, &
         data_functions_copy, data_functions_dealloc, data_fftbox_alloc, &
         data_fftbox_dealloc, data_set_to_zero
    use dma, only: dma_calculate, dma_free_multipoles
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_batch_col_start
    use fourier, only: fourier_apply_box_pair, fourier_interpolate, &
         fourier_filter
    use geometry, only: POINT, operator(-)
    use hash_table, only: HT_HASH_TABLE, hash_table_init, hash_table_free, &
         hash_table_add, hash_table_probe, hash_table_lookup_nocount
    use model_type, only: MODEL
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET, &
         multipole_free_spherical_set
    use parallel_strategy, only: PARAL_INFO
    use ppd_ops, only: ppd_product, ppd_accumulate
    use remote, only: remote_obtain_ngwf_by_copy, packed_ngwf_size, &
         remote_serve_ngwfs, remote_pack_data_real, NGWF_REQUEST_TAG
    use rundat, only: pub_dma_use_ri, pub_pol_emb_pot_filename, pub_num_spins, &
         pub_dma_bessel_averaging, pub_dma_scale_charge, pub_spin_fac, &
         pub_precond_recip, pub_output_detail, pub_pol_emb_dbl_grid, PUB_1K, &
         pub_devel_code, pub_pol_emb_vacuum_qmstar
    use sparse, only: sparse_get_element
    use sparse_array, only: SPAM3_ARRAY
    use sw_expansion_type, only: SW_EX
    use sw_resolution_of_identity, only: SW_RI, swri_library, swri_get_handle_to
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check, &
         utils_sanity_check, utils_int_to_str, utils_devel_code, utils_abort

    implicit none

    ! jd: Arguments
    type(SW_EX), intent(in), target        :: swexes(2)
    type(FUNCTIONS), intent(in)            :: ngwfs_on_grid
    type(HT_HASH_TABLE), intent(inout)     :: packed_ngwf_cache_hts(:)
    type(SPAM3_ARRAY), intent(in)          :: denskern
    type(FUNC_BASIS), intent(in)           :: ngwf_basis
    type(MODEL), intent(in), target        :: mdl
    type(SPAM3_ARRAY), intent(in)          :: tc_denskern
    type(FUNCTIONS), intent(inout)         :: cov_grad
    type(FUNCTIONS), intent(inout)         :: contra_grad
    real(kind=DP), intent(in)              :: precond_func_recip(:,:,:)
    integer, intent(in)                    :: ngwf_cache_handle
    real(kind=DP), intent(in), optional    :: prefactor

    ! jd: --- Local variables ---

    type(SW_EX), pointer :: swex_2
    type(SW_RI), pointer :: swri

    ! jd: Atom bookkeeping
    integer :: local_a, global_a, orig_global_a
    integer :: local_b, global_b, b_idx, orig_global_b

    ! jd: NGWF bookkeeping
    integer :: first_ngwf_idx_of_a, ngwf_a, global_aa_ngwf_idx
    integer :: first_ngwf_idx_of_b, ngwf_b, global_bb_ngwf_idx
    integer :: local_idx_of_1st_ngwf_on_a, local_idx_of_1st_ngwf_on_b
    integer :: local_aa_ngwf_idx
    integer :: n_ngwfs_on_a
    real(kind=DP) :: packed_ngwf_bb(packed_ngwf_size)
    real(kind=DP) :: packed_initial_ngwf_bb(packed_ngwf_size)
    real(kind=DP) :: packed_ngwf_bb_delta(packed_ngwf_size)
    integer :: ndata

    ! jd: DKN elements
    real(kind=DP) :: dkn_aa_bb_is, tc_dkn_aa_bb_is

    ! jd: PPD bookkeeping
    integer :: n_ppds_on_a, n_ppds_on_b, n_ppds_in_prod
    integer :: ppd_indices_on_a(ngwf_basis%max_n_ppds_sphere)
    integer :: ppd_indices_on_b(ngwf_basis%max_n_ppds_sphere)
    integer :: ppd_indices_in_prod(ngwf_basis%max_n_ppds_sphere)
    integer :: cur_ppd_n, cur_ppd_a
    integer :: offs_to_cur_ppd_a
    integer :: dot_npts
    integer :: bra_start
    integer :: ipt
    type(FUNCTIONS)  :: ket_on_bragrid_buffer(2) ! 1st: K, 2nd: tcK

    ! jd: Bessel averaging bookkeeping
    integer :: bessel_run, n_bessel_runs

    ! jd: Comms
    integer :: proc
    logical :: who_is_done_ngwf_xchg(0:pub_total_num_procs-1)
    integer :: send_handles(0:pub_total_num_procs-1)
    logical :: done_ngwf

    ! jd: Auxiliary terms
    type(HT_HASH_TABLE), target :: qterms_ht
    type(HT_HASH_TABLE), target :: jqterms_ht

    real(kind=DP) :: qterm(n_qterm_sph_components) ! Q_B^{lm}
    real(kind=DP) :: qterm_cur_set(n_qterm_sph_components) ! Q_B^{lm}
    real(kind=DP), allocatable :: jqterm(:) ! JQ_A^{s_A}.
    real(kind=DP) :: uabterm_in_ppds_on_a(MDL_PPDS_DIMS)
    real(kind=DP) :: uab00term_in_ppds_on_a(MDL_PPDS_DIMS)
    real(kind=DP) :: uab_ngwf_bb_prod_in_ppds_on_a(MDL_PPDS_DIMS)
    real(kind=DP) :: packed_uab_in_ppds_on_a(packed_ngwf_size)
    real(kind=DP) :: qsc1term_in_ppds_on_a(MDL_PPDS_DIMS)
    real(kind=DP) :: qsc2term_in_ppds_on_a(MDL_PPDS_DIMS)
    real(kind=DP) :: qsc2_ngwf_bb_prod_in_ppds_on_a(MDL_PPDS_DIMS)
    real(kind=DP) :: packed_qsc2_in_ppds_on_a(packed_ngwf_size)
    real(kind=DP), allocatable :: tau0term_in_ppds_on_a(:,:,:) ! 1st: ppdpt, is, ngwf
    real(kind=DP), allocatable :: tau1term_in_ppds_on_a(:,:,:) ! 1st: ppdpt, is, ngwf
    real(kind=DP) :: q_DMA(2) ! idx: Bessel runs
    real(kind=DP) :: sum_corr_terms(2,MAX_N_MULTIPOLE_SETS) ! idx1: Bessel runs, MM sets
    real(kind=DP) :: lambda(2) ! idx: Bessel runs

    ! jd: Gradient kets in PPDs
    type(FUNCTIONS), allocatable :: kets_in_ppds_on_a(:,:,:)

    ! jd: FFTbox-based double grid version
    real(kind=DP) :: uabterm_in_fftbox_dbl_on_a(FFTBOX_DBL_DIMS)
    ! jd: sum_B K^AaBb \Phi Bb in FFTbox around Aa: idx1: a idx2: spin idx3: k/tck
    type(FFTBOX_DATA), allocatable :: dkn_phi_bb_in_fftbox_on_a(:,:,:)
    type(FFTBOX_DATA), allocatable :: uab_dkn_phi_bb_in_fftbox_on_a(:,:,:)
    type(FFTBOX_DATA), allocatable :: sum_b_uab_dkn_phi_bb_in_fftbox_on_a(:,:,:)

    ! jd: Work arrays for the double grid version
    real(kind=DP)    :: dwork_box_dbl(FFTBOX_DBL_DIMS,2)
    complex(kind=DP) :: dbl_work(FFTBOX_DBL_DIMS)
    complex(kind=DP) :: coarse_work(FFTBOX_DIMS)

    ! jd: MM side
    integer                                    :: mpole_set
    integer                                    :: n_sets
    type(SPHERICAL_MULTIPOLE_SET), allocatable :: mm_multipoles(:)
    integer                                    :: n_energy_terms
    character(len=27), allocatable             :: dummy_energy_term_names(:)
    real(kind=DP), allocatable                 :: dummy_energy_term_values(:)

    ! jd: Reciprocal-space preconditioning
    integer          :: fa_box_start(3,1)
    integer          :: fa_start_in_box(3,1)
    type(FFTBOX_DATA), allocatable :: kets_in_fftbox(:,:)
    real(kind=DP), allocatable :: work_fftbox(:,:,:)
    complex(kind=DP), allocatable :: zwork_fftbox(:,:,:)

    ! jd: Varia
    integer       :: is
    integer       :: k_tck
    real(kind=DP) :: alpha_a
    type(POINT)   :: vec_a
    real(kind=DP) :: factor
    integer       :: dma_swri_h
    integer       :: ierr
    logical       :: debug_show_skipped
    type(SPHERICAL_MULTIPOLE_SET) :: scf_dma_multipoles(3)
    character(len=*), parameter :: myself = 'pol_emb_ngwf_grad_contribution'
    type(PARAL_INFO), pointer :: par

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    call utils_sanity_check(contra_grad%d, &
         'contra_grad on entry to '//myself, excessive=1D20)
    call utils_sanity_check(cov_grad%d,&
         'contra_grad on entry to '//myself, excessive=1D20)

#if 0
       print *, "-NGWF_CONTRA_PRE:-============================@"
       print *, contra_grad%d
       print *, "-NGWF_COV_PRE:-============================@"
       print *, cov_grad%d
#endif

    dma_swri_h = swri_get_handle_to(pub_dma_use_ri)

    debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
         pub_devel_code, no_warn=.true.)

    par=>mdl%regions(1)%par

    ! jd: Read Cartesian multipoles from file, convert to sphericals
    call pol_emb_read_from_file(&
         dummy_energy_term_names, dummy_energy_term_values, n_energy_terms, &
         mm_multipoles, n_sets, pub_pol_emb_pot_filename, silent = .true.)

    ! jd: If no sets with QM* gradient, bail out. Clean up and stop timer first.
    if(all(index(mm_multipoles(:)%set_name,'qmstar')==0)) then
       goto 999
    end if

    call utils_assert(.not. (pub_dma_scale_charge .and. pub_pol_emb_dbl_grid), &
         myself//': Charge scaling in pol_emb NGWF gradient is only supported &
         &on the coarse grid (do not use pol_emb_dbl_grid T)')

    call utils_assert(.not. ngwfs_on_grid%iscmplx, &
         myself//' not ready for complex NGWFs')

    who_is_done_ngwf_xchg(:) = .false.

    swri => swri_library(dma_swri_h)
    if(pub_dma_bessel_averaging) then
       n_bessel_runs = 2
    else
       n_bessel_runs = 1
    end if

    ! jd: Do a DMA calculation, needed only for charge-scaling correction,
    !     which needs the charge scaling factor and the monopoles. Alternative
    !     would be to store those in SW_EX rather than in scf_dma_multipoles.
    if(pub_dma_scale_charge) then
       call utils_abort(myself//': Charge-scaling not implemented in NGWF grad')
#if 0
       call dma_calculate(scf_dma_multipoles, &
            merged_swexes, &
            denskern, rep%n_occ, rep%overlap, ngwf_basis, mdl, &
            dkn_scale_factor = pub_spin_fac*0.5_DP)
#endif
    end if
    ! jd: Create hash tables for Q and JQ. Q is indexed with the atom number,
    !     so there will never be more than nat entries. JQ is indexed with
    !     atom number and Bessel averaging run, so there will never be more than
    !     2*nat entries. We use these upper bounds to ensure the HT does not
    !     overfill, but only pay for 1000 entries of primary storage, which is
    !     next to nothing.
    call hash_table_init(qterms_ht, 'QTERMS', 1, par%nat, 1000, stdout)
    call hash_table_init(jqterms_ht, 'JQTERMS', 2, n_bessel_runs * par%nat, &
         1000, stdout)

    do k_tck = 1, 2
       call data_functions_alloc(ket_on_bragrid_buffer(k_tck), &
            MDL_PPDS_DIMS, iscmplx=.false.)
    end do

    allocate(jqterm(swexes(1)%quality%num_sws_per_centre), stat=ierr)
    call utils_alloc_check(myself,'jqterm',ierr)

    ! **************************************************************************
    ! **************************************************************************
    ! Preparation-1: set up Q terms, JQ terms and NGWFs BB for all atoms of
    !                revelance -- all A's and their S-neighbours B.
    ! **************************************************************************
    ! **************************************************************************

    ! --------------------------------------------------------------------------
    ! jd: Loop over A's -- all atoms on this rank                           AAA
    ! --------------------------------------------------------------------------
    loop_A1:                                                                  &
    do local_a = 1, par%num_atoms_on_proc(pub_my_proc_id)
       global_a = par%first_atom_on_proc(pub_my_proc_id) + local_a-1
       orig_global_a = par%orig_atom(global_a)

       ! jd: Atoms outside this SWRI should *not* be skipped, or else the NGWF
       !     gradient will be wrong for subsystem DMA. This is a result of
       !     elements of the P matrix being generated through symmetries.
       !     cf. page (4B) in gradient notes.

       ! -----------------------------------------------------------------------
       ! jd: Loop over B's that are s-neighbours with A                     BBB
       ! -----------------------------------------------------------------------
       loop_B1:                                                               &
       do b_idx = swri%s_atoms_nl%first_idx(global_a), &
            swri%s_atoms_nl%last_idx(global_a)
          global_b = swri%s_atoms_nl%neighbours(b_idx)
          orig_global_b = par%orig_atom(global_b)
          first_ngwf_idx_of_b = ngwf_basis%first_on_atom(global_b)

          ! jd: Atoms outside this SWRI should *not* be skipped, or else the NGWF
          !     gradient will be wrong for subsystem DMA. This is a result of
          !     elements of the P matrix being generated through symmetries.
          !     cf. page (4B) in gradient notes.

          ! jd: If quantities of interest have not been calculated for this atom
          !     yet, calculate them
          if(-1 == hash_table_probe(qterms_ht, global_b)) then
             ! -----------------------------------------------------------------
             ! jd: Loop over multipole sets                                  MS
             ! -----------------------------------------------------------------
             qterm(:) = 0D0
             loop_set_1:                                                      &
             do mpole_set = 1, n_sets
                if(index(mm_multipoles(mpole_set)%set_name,'qmstar')/=0) then
                   ! jd: Calculate auxiliary term Q_B
                   call pol_emb_qterm(mm_multipoles(mpole_set), mdl%elements, &
                        par, global_b, swexes(1)%quality%min_l, &
                        swexes(1)%quality%max_l, qterm_sph = qterm_cur_set)
                   qterm = qterm + qterm_cur_set

                   if(pub_output_detail >= VERBOSE .and. pub_on_root .and. &
                        local_a == 1 .and. &
                        b_idx == swri%s_atoms_nl%first_idx(global_a)) then
                      if(mpole_set == 1) write(stdout,'(/a)')
                      write(stdout,'(a,i0,a,i0)') &
                           'Polarisable embedding QM* NGWF gradient for set '//&
                           trim(mm_multipoles(mpole_set)%set_name)//'. l: ', &
                           swexes(1)%quality%min_l, '-', swexes(1)%quality%max_l
                   end if
                end if
             end do loop_set_1
             ! -----------------------------------------------------------------

             call utils_sanity_check(qterm,'qterm in '//myself, excessive=1D20)
 !            write(*,*) '@@@STORE: ', qterm
             call hash_table_add(qterms_ht, qterm, n_qterm_sph_components, &
                  global_b)

             ! -----------------------------------------------------------------
             ! jd: for all NGWFs b on B                                     bbb
             ! -----------------------------------------------------------------
             loop_ngwf_b_1:                                                   &
             do ngwf_b = 1, ngwf_basis%num_on_atom(global_b)
                global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

                ! jd: Obtain all NGWFs Bb that will be required. These are
                !     stored in the HT.
                call remote_obtain_ngwf_by_copy(&
                     packed_ngwf_bb, packed_ngwf_cache_hts, &
                     who_is_done_ngwf_xchg, global_bb_ngwf_idx, &
                     ngwfs_on_grid, ngwf_basis, ngwf_cache_handle, &
                     overfill_strategy='F')

             end do loop_ngwf_b_1
             ! -----------------------------------------------------------------

             ! -----------------------------------------------------------------
             ! jd: Loop over Bessel averaging runs                          AVG
             ! -----------------------------------------------------------------
             loop_bessel_avg_run_1:                                           &
             do bessel_run = 1, n_bessel_runs
                swex_2 => swexes(bessel_run)
                ! jd: Calculate auxiliary term JQ
                call pol_emb_jqterm(jqterm, qterm, swri, swex_2)
                call utils_sanity_check(jqterm, 'jqterm in '//myself, &
                     excessive=1D20)
                call hash_table_add(jqterms_ht, jqterm, &
                     swex_2%quality%num_sws_per_centre, global_b, bessel_run)

             end do loop_bessel_avg_run_1
             ! -----------------------------------------------------------------

          end if ! New atom B

       end do loop_B1
       ! -----------------------------------------------------------------------

    end do loop_A1
    ! --------------------------------------------------------------------------

    ! jd: Ensure all NGWF comms have finished
    if(pub_total_num_procs > 1) then
       do proc = 0, pub_total_num_procs-1
          call comms_send(proc, -1, 1, tag = NGWF_REQUEST_TAG + 1, & ! [x]
               return_handle = send_handles(proc), add_to_stack = .false.)
               ! [x]: For a 'done' request it's OK to only send it for
               !      NGWF set 1 -- it's understood that *all* sets finished
               !      participating.
       end do

       done_ngwf = .false.
       if(pub_total_num_procs > 1) then
          do while(.not. done_ngwf)
             call remote_serve_ngwfs(done_ngwf, who_is_done_ngwf_xchg, &
                  ngwfs_on_grid, ngwf_basis, ngwf_cache_handle)
          end do
       end if

       do proc = 0, pub_total_num_procs-1
          call comms_wait(send_handles(proc))
       end do
    end if

    ! **************************************************************************
    ! **************************************************************************
    ! Preparation-2: calculate sum_corr_terms, scaling lambda, q_DMA
    ! **************************************************************************
    ! **************************************************************************
    if(pub_dma_scale_charge) then
       ! -----------------------------------------------------------------------
       ! jd: Loop over multipole sets                                        MS
       ! -----------------------------------------------------------------------
       loop_set_2:                                                            &
       do mpole_set = 1, n_sets
          if(index(mm_multipoles(mpole_set)%set_name,'qmstar')/=0) then

             ! -----------------------------------------------------------------
             ! jd: Loop over Bessel averaging runs                          AVG
             ! -----------------------------------------------------------------
             loop_bessel_avg_run_2:                                           &
             do bessel_run = 1, n_bessel_runs
                call pol_emb_qterms_and_sum_corr_terms(ngwf_basis, &
                     scf_dma_multipoles(bessel_run), mm_multipoles(mpole_set), &
                     mdl%elements, par, sum_corr_terms(bessel_run, mpole_set), &
                     lambda(bessel_run), q_DMA(bessel_run))
             end do loop_bessel_avg_run_2
          end if

       end do loop_set_2
       ! -----------------------------------------------------------------------
    end if

    ! **************************************************************************
    ! **************************************************************************
    ! 4: Main calculation
    ! **************************************************************************
    ! **************************************************************************

    ! --------------------------------------------------------------------------
    ! jd: Loop over A's -- all atoms on this rank                           AAA
    ! --------------------------------------------------------------------------
    loop_A4:                                                                  &
    do local_a = 1, par%num_atoms_on_proc(pub_my_proc_id)
       global_a = par%first_atom_on_proc(pub_my_proc_id) + local_a-1
       first_ngwf_idx_of_a = ngwf_basis%first_on_atom(global_a)
       local_idx_of_1st_ngwf_on_a = first_ngwf_idx_of_a + 1 - &
            ngwf_basis%first_on_proc(pub_my_proc_id)
       orig_global_a = par%orig_atom(global_a)
       n_ngwfs_on_a = ngwf_basis%num_on_atom(global_a)

       ! jd: Atoms outside this SWRI should *not* be skipped, or else the NGWF
       !     gradient will be wrong for subsystem DMA. This is a result of
       !     elements of the P matrix being generated through symmetries.
       !     cf. page (4B) in gradient notes.

       ! jd: Thole polarisability and position of atom A
       alpha_a = mdl%elements(orig_global_a)%thole_polarisability
       vec_a = mdl%elements(orig_global_a)%centre

       ! jd: Calculate auxiliary quantities tau0 and tau1 iff charge-scaling
       !     is in effect
       if(pub_dma_scale_charge) then
#if 0
          allocate(tau0term_in_ppds_on_a(MDL_PPDS_DIMS, &
               pub_num_spins, n_ngwfs_on_a), stat=ierr)
          call utils_alloc_check(myself, 'tau0term_in_ppds_on_a', ierr)
          allocate(tau1term_in_ppds_on_a(MDL_PPDS_DIMS, &
               pub_num_spins, n_ngwfs_on_a), stat=ierr)
          call utils_alloc_check(myself, 'tau1term_in_ppds_on_a', ierr)
          call pol_emb_tau0term(tau0term_in_ppds_on_a, local_a, swri%s_atoms_nl, &
               ngwf_basis, rep, denskern%m(:,PUB_1K), mdl%cell, mdl%elements)
          call pol_emb_tau1term(tau1term_in_ppds_on_a, local_a, swri%s_atoms_nl, &
               ngwf_basis, rep, denskern%m(:,PUB_1K), mdl%cell, &
               jqterms_ht, swri, mdl%elements)
#endif
       end if

       allocate(kets_in_ppds_on_a(n_ngwfs_on_a,pub_num_spins,2), stat=ierr)
       call utils_alloc_check(myself,'kets_in_ppds_on_a',ierr)

       ! jd: Allocations for dbl grid branch
       if(pub_pol_emb_dbl_grid) then
          allocate(dkn_phi_bb_in_fftbox_on_a(n_ngwfs_on_a,pub_num_spins,2), &
               stat=ierr)
          call utils_alloc_check(myself,'dkn_phi_bb_in_fftbox_on_a',ierr)
          allocate(uab_dkn_phi_bb_in_fftbox_on_a(n_ngwfs_on_a,pub_num_spins,2),&
               stat=ierr)
          call utils_alloc_check(myself,'uab_dkn_phi_bb_in_fftbox_on_a',ierr)
          allocate(sum_b_uab_dkn_phi_bb_in_fftbox_on_a(&
               n_ngwfs_on_a,pub_num_spins, 2), stat=ierr)
          call utils_alloc_check(myself,'sum_b_uab_dkn_phi_bb_in_fftbox_on_a', &
               ierr)
       end if

       ! jd: Allocate work datastructures
       do k_tck = 1, 2
          do is = 1, pub_num_spins
             do ngwf_a = 1, n_ngwfs_on_a
                call data_functions_alloc(kets_in_ppds_on_a(ngwf_a, is, k_tck),&
                     MDL_PPDS_DIMS, iscmplx=ngwfs_on_grid%iscmplx)
                if(pub_pol_emb_dbl_grid) then
                   call data_fftbox_alloc(dkn_phi_bb_in_fftbox_on_a(&
                        ngwf_a,is,k_tck), FFTBOX_DIMS, &
                        iscmplx=contra_grad%iscmplx)      !jmecmplx
                   call data_fftbox_alloc(uab_dkn_phi_bb_in_fftbox_on_a(&
                        ngwf_a,is,k_tck), FFTBOX_DIMS, &
                        iscmplx=contra_grad%iscmplx)      !jmecmplx
                   call data_fftbox_alloc(sum_b_uab_dkn_phi_bb_in_fftbox_on_a(&
                        ngwf_a,is,k_tck), FFTBOX_DIMS, &
                        iscmplx=contra_grad%iscmplx)      !jmecmplx
                end if
             end do
          end do
       end do

       ! jd: Start new sum over B
       call data_set_to_zero(kets_in_ppds_on_a)
       if(pub_pol_emb_dbl_grid) then
          call data_set_to_zero(sum_b_uab_dkn_phi_bb_in_fftbox_on_a(:,:,:))
       end if

       ! -----------------------------------------------------------------------
       ! jd: Loop over B's that are s-neighbours with A                     BBB
       ! -----------------------------------------------------------------------
       loop_B4:                                                               &
       do b_idx = swri%s_atoms_nl%first_idx(global_a), &
            swri%s_atoms_nl%last_idx(global_a)
          global_b = swri%s_atoms_nl%neighbours(b_idx)
          orig_global_b = par%orig_atom(global_b)
          first_ngwf_idx_of_b = ngwf_basis%first_on_atom(global_b)
          local_idx_of_1st_ngwf_on_b = first_ngwf_idx_of_b + 1 - &
               ngwf_basis%first_on_proc(pub_my_proc_id)

          ! jd: No skipping of atoms outside SWRI here

          if(pub_pol_emb_dbl_grid) then ! jd: Since pol_emb_uabterm adds to
             uabterm_in_fftbox_dbl_on_a(:,:,:) = 0D0
          else
             uabterm_in_ppds_on_a(:) = 0D0
          end if

          ! --------------------------------------------------------------------
          ! jd: Loop over Bessel averaging runs                             AVG
          ! --------------------------------------------------------------------
          loop_bessel_avg_run_4:                                              &
          do bessel_run = 1, n_bessel_runs
             swex_2 => swexes(bessel_run)

             ! DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL
             if(pub_pol_emb_dbl_grid) then
                ! jd: Calculate uabterm in a double FFTbox around A
                !     Also determine the PPDs on A
                call pol_emb_uabterm(swri, ngwf_basis, mdl%cell, mdl%elements, par, &
                     local_a, global_b, jqterms_ht, swex_2, bessel_run, .false., &
                     uabterm_in_fftbox_dbl_on_a = uabterm_in_fftbox_dbl_on_a, &
                     fftbox = mdl%fftbox, &
                     ppd_indices_on_a = ppd_indices_on_a, &
                     n_ppds_on_a = n_ppds_on_a)
!@@@                     hack_mask = pub_pol_emb_vacuum_qmstar)
                call utils_sanity_check(uabterm_in_fftbox_dbl_on_a, &
                     'U_AB (fftbox) term in '//myself, excessive=1D20)
             ! COARSE - COARSE - COARSE - COARSE - COARSE - COARSE - COARSE
             else
                ! jd: Calculate uabterm, determine PPDs on A
                call pol_emb_uabterm(swri, ngwf_basis, mdl%cell, mdl%elements, par, &
                     local_a, global_b, jqterms_ht, swex_2, bessel_run, .false., &
                     uabterm_in_ppds_on_a = uabterm_in_ppds_on_a, &
                     ppd_indices_on_a = ppd_indices_on_a, &
                     n_ppds_on_a = n_ppds_on_a)
!@@@                     hack_mask = pub_pol_emb_vacuum_qmstar)
                call utils_sanity_check(uabterm_in_ppds_on_a, &
                     'U_AB (ppd) term in '//myself, excessive=1D20)

#if 0
                if(pub_dma_scale_charge) then
                   ! jd: Calculate uab00term, unnecessarily redetermine PPDs on A
                   call pol_emb_uabterm(swri, ngwf_basis, mdl%cell, mdl%elements, par, &
                        local_a, global_b, jqterms_ht, swex_2, bessel_run, .true., &
                        uabterm_in_ppds_on_a = uab00term_in_ppds_on_a, &
                        ppd_indices_on_a = ppd_indices_on_a, &
                        n_ppds_on_a = n_ppds_on_a)
                        hack_mask ??
                   call utils_sanity_check(uabterm_in_ppds_on_a, &
                        'U_AB (ppd) term in '//myself, excessive=1D20)

                   ! jd: Calculate qscterm1 (averages over bessel runs internally)
                   call pol_emb_qsc1term(ngwf_basis, mdl%cell, local_a, global_b, &
                        qterms_ht, bessel_run, rep%overlap, tau0term_in_ppds_on_a, &
                        tau1term_in_ppds_on_a, n_sets, sum_corr_terms, &
                        lambda, q_DMA, qsc1term_in_ppds_on_a)
                   call utils_sanity_check(qsc1term_in_ppds_on_a, &
                        'QSC1 (ppd) term in '//myself, excessive=1D20)

                   ! jd: Calculate qscterm2 (averages over bessel runs internally)
                   call pol_emb_qsc2term(ngwf_basis, mdl%cell, par, local_a, global_b, &
                        qterms_ht, rep%remote_ngwf_cache_ht, &
                        uab00term_in_ppds_on_a, bessel_run, &
                        n_sets, sum_corr_terms, lambda, q_DMA, qsc2term_in_ppds_on_a)
                   call utils_sanity_check(qsc2term_in_ppds_on_a, &
                        'QSC2 (ppd) term in '//myself, excessive=1D20)
                end if ! charge scaling?
#endif
             end if ! dbl/coarse

          end do loop_bessel_avg_run_4
          ! --------------------------------------------------------------------

          ! DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL
          if(pub_pol_emb_dbl_grid) then
             ! jd: Scale sums by 1/2 if Bessel averaging.
             if(pub_dma_bessel_averaging) then
                uabterm_in_fftbox_dbl_on_a = uabterm_in_fftbox_dbl_on_a * 0.5_DP
             end if ! [1]
             ! jd: Slightly unsure why we scale twice ([1] and [2]), but
             !     FD tests seem to indicate this is the right thing to do.
             ! jd: @Likely because we use use_1_for_bra later, the scaling in [2]
             !     does not really matter.
             ! jd: Use 1.0 as U_AB when using double grid (or 0.5 when averaging)
             call remote_pack_data_real(packed_uab_in_ppds_on_a, &
                  uabterm_in_ppds_on_a, ngwf_basis, local_idx_of_1st_ngwf_on_a, &
                  override_value = merge(0.5D0,1D0,pub_dma_bessel_averaging)) ! [2]
          ! COARSE - COARSE - COARSE - COARSE - COARSE - COARSE - COARSE
          else
             ! jd: Scale sums by 1/2 if Bessel averaging.
             if(pub_dma_bessel_averaging) then
                uabterm_in_ppds_on_a = uabterm_in_ppds_on_a * 0.5_DP
! @scaled internally?
!                if(pub_dma_scale_charge) then
!                   qsc1term_in_ppds_on_a = qsc1term_in_ppds_on_a * 0.5_DP
!                   qsc2term_in_ppds_on_a = qsc2term_in_ppds_on_a * 0.5_DP
!                end if
             end if

             ! jd: Pack U_AB for ppd_product
             call remote_pack_data_real(packed_uab_in_ppds_on_a, &
                  uabterm_in_ppds_on_a, ngwf_basis, local_idx_of_1st_ngwf_on_a)

             ! jd: Pack QSCterm2 for ppd_product
             if(pub_dma_scale_charge) then
                call remote_pack_data_real(packed_qsc2_in_ppds_on_a, &
                     qsc2term_in_ppds_on_a, ngwf_basis, &
                     local_idx_of_1st_ngwf_on_a)
             end if
          end if

          ! jd: The meaning of kets_in_ppds_on_a changes for the double grid
          !     case, and it must be re-initialised for every atom B.
          if(pub_pol_emb_dbl_grid) then
             call data_set_to_zero(kets_in_ppds_on_a)
          end if

          ! --------------------------------------------------------------------
          ! jd: for all NGWFs b on B                                        bbb
          ! --------------------------------------------------------------------
          loop_ngwf_b4:                                                       &
          do ngwf_b = 1, ngwf_basis%num_on_atom(global_b)
             global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

             ! jd: Retrieve packed NGWF Bb from cache
             call hash_table_lookup_nocount(packed_ngwf_bb, ndata, &
                  packed_ngwf_cache_hts(ngwf_cache_handle), global_bb_ngwf_idx)
             call utils_assert(ndata /= -1, myself//': NGWF not found in NGWF &
                  &cache. This indicates a bug in the logic.')

             ! jd: Calculate U_AB * \Phi_Bb in PPDs on A (if coarse grid)
             !     or               \Phi_Bb in PPDs on A (if double grid)
             call ppd_product(packed_uab_in_ppds_on_a, packed_ngwf_bb, &
                  mdl%cell%n_pts, uab_ngwf_bb_prod_in_ppds_on_a, &
                  ppd_indices_in_prod, n_ppds_in_prod, &
                  use_1_for_bra = pub_pol_emb_dbl_grid)
             call utils_sanity_check(uab_ngwf_bb_prod_in_ppds_on_a, &
                  'U_AB * \Phi_Bb term in '//myself,excessive=1D20)

             ! jd: Calculate charge-scaling correction term (notes, p. 14)
             if(pub_dma_scale_charge) then
                call ppd_product(packed_qsc2_in_ppds_on_a, packed_ngwf_bb, &
                     mdl%cell%n_pts, qsc2_ngwf_bb_prod_in_ppds_on_a)
                call utils_sanity_check(uab_ngwf_bb_prod_in_ppds_on_a, &
                     'charge scaling term2 term in '//myself,excessive=1D20)
             end if

             ! -----------------------------------------------------------------
             ! jd: for all spins                                            sss
             ! -----------------------------------------------------------------
             loop_spins4:                                                     &
             do is = 1, pub_num_spins

                ! -------------------------------------------------------------
                ! jd: for all NGWFs a on A                                  aaa
                ! -------------------------------------------------------------
                loop_ngwf_a4:                                                 &
                do ngwf_a = 1, n_ngwfs_on_a
                   global_aa_ngwf_idx = first_ngwf_idx_of_a + ngwf_a - 1

                   call sparse_get_element(dkn_aa_bb_is, denskern%m(is,PUB_1K),&
                        global_bb_ngwf_idx, global_aa_ngwf_idx)
                   call sparse_get_element(tc_dkn_aa_bb_is, tc_denskern%m(is,PUB_1K),&
                        global_bb_ngwf_idx, global_aa_ngwf_idx)

                   ! jd: Accumulate U_AB * K^AaBb * \Phi_Bb
                   !            and U_AB * tcK^AaBb * \Phi_Bb
                   !            (directly on coarse grid)
                   !     or
                   ! jd: Accumulate        K^AaBb * \Phi_Bb
                   !                       tcK^AaBb * \Phi_Bb
                   !   (to be multiplied by U_AB later on double grid)
                   call ppd_accumulate(kets_in_ppds_on_a(ngwf_a,is,1)%d, &
                        ppd_indices_on_a, n_ppds_on_a, &
                        uab_ngwf_bb_prod_in_ppds_on_a(:), &
                        ppd_indices_in_prod, n_ppds_in_prod, .true., &
                        mdl%cell%n_pts, -dkn_aa_bb_is) ! sic, -

                   call ppd_accumulate(kets_in_ppds_on_a(ngwf_a,is,2)%d, &
                        ppd_indices_on_a, n_ppds_on_a, &
                        uab_ngwf_bb_prod_in_ppds_on_a(:), &
                        ppd_indices_in_prod, n_ppds_in_prod, .true., &
                        mdl%cell%n_pts, -tc_dkn_aa_bb_is) ! sic, -

                end do loop_ngwf_a4
                ! --------------------------------------------------------------
             end do loop_spins4
             ! -----------------------------------------------------------------
          end do loop_ngwf_b4
          ! --------------------------------------------------------------------

          ! DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL
          if(pub_pol_emb_dbl_grid) then

             ! jd: Determine the position of the FFTbox around A's 1st NGWF
             !     and its position in the FFTbox
             call function_ops_batch_col_start(fa_box_start, fa_start_in_box, &
                  1, local_idx_of_1st_ngwf_on_a, local_idx_of_1st_ngwf_on_a, &
                  mdl%fftbox, mdl%cell, ngwf_basis)

             ! jd: Transfer sum_b K^AaBb   * \Phi_Bb and
             !              sum_b tcK^AaBb * \Phi_Bb from PPDs to A's FFTbox
             do is = 1, pub_num_spins
                do ngwf_a = 1, ngwf_basis%num_on_atom(global_a)

                   global_aa_ngwf_idx = &
                        ngwf_basis%first_on_atom(global_a) + ngwf_a - 1
                   local_aa_ngwf_idx = global_aa_ngwf_idx - &
                        ngwf_basis%first_on_proc(pub_my_proc_id)+1

                   do k_tck = 1, 2
                      ! jd: Transfer data from PPDs to an FFTbox around Aa.
                      call basis_copy_function_to_box(&
                           dkn_phi_bb_in_fftbox_on_a(ngwf_a,is,k_tck), & ! target
                           fa_start_in_box(1,1), &
                           fa_start_in_box(2,1), &
                           fa_start_in_box(3,1), &
                           ngwf_basis%tight_boxes(local_aa_ngwf_idx), &
                           kets_in_ppds_on_a(ngwf_a,is,k_tck), &
                           ngwf_basis%spheres(local_aa_ngwf_idx), &
                           mdl%cell, mdl%fftbox, 1)
                   end do

                   ! jd: Multiply sum_b (tc)K^AaBb * \Phi_Bb with U_AB potential
                   !     on double grid
                   call fourier_interpolate(coarse_work, dbl_work, &
                        dkn_phi_bb_in_fftbox_on_a(ngwf_a,is,1)%d, &
                        dkn_phi_bb_in_fftbox_on_a(ngwf_a,is,2)%d, &
                        dwork_box_dbl(:,:,:,1),dwork_box_dbl(:,:,:,2))

                   dwork_box_dbl(:,:,:,1) = &
                        dwork_box_dbl(:,:,:,1) * uabterm_in_fftbox_dbl_on_a
                   dwork_box_dbl(:,:,:,2) = &
                        dwork_box_dbl(:,:,:,2) * uabterm_in_fftbox_dbl_on_a

                   call fourier_filter(coarse_work, dbl_work, &
                        dwork_box_dbl(:,:,:,1), dwork_box_dbl(:,:,:,2), &
                        uab_dkn_phi_bb_in_fftbox_on_a(ngwf_a,is,1)%d, &
                        uab_dkn_phi_bb_in_fftbox_on_a(ngwf_a,is,2)%d)

                end do ! ngwf_a
             end do ! is

             ! jd: Accumulate sum over atoms B
             do k_tck = 1, 2
                do is = 1, pub_num_spins
                   do ngwf_a = 1, ngwf_basis%num_on_atom(global_a)
                      sum_b_uab_dkn_phi_bb_in_fftbox_on_a(ngwf_a,is,k_tck)%d(:,:,:) = &
                           sum_b_uab_dkn_phi_bb_in_fftbox_on_a(ngwf_a,is,k_tck)%d + &
                           uab_dkn_phi_bb_in_fftbox_on_a(ngwf_a,is,k_tck)%d
                   end do
                end do
             end do

          end if ! double grid
          ! DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL

       end do loop_B4
       ! -----------------------------------------------------------------------

       ! *******************************************
       ! jd: Calculate the gradient proper
       ! *******************************************

       ! -----------------------------------------------------------------------
       ! jd: for all spins                                                  sss
       ! -----------------------------------------------------------------------
       loop_spins_4b:                                                         &
       do is = 1, pub_num_spins

          ! --------------------------------------------------------------------
          ! jd: for all NGWFs a on A                                        aaa
          ! --------------------------------------------------------------------
          loop_ngwf_a4b:                                                       &
          do ngwf_a = 1, ngwf_basis%num_on_atom(global_a)

             global_aa_ngwf_idx = ngwf_basis%first_on_atom(global_a) + ngwf_a - 1
             local_aa_ngwf_idx = global_aa_ngwf_idx - &
                  ngwf_basis%first_on_proc(pub_my_proc_id)+1

             ! DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL - DBL
             if(pub_pol_emb_dbl_grid) then

                ! jd: Extract from FFTbox to PPDs
                do k_tck = 1, 2
                   call basis_extract_function_from_box(&
                        ket_on_bragrid_buffer(k_tck), &
                        sum_b_uab_dkn_phi_bb_in_fftbox_on_a(ngwf_a,is,k_tck), &  !jmecmplx
                        ngwf_basis%spheres(local_aa_ngwf_idx), &
                        ngwf_basis%tight_boxes(local_aa_ngwf_idx), &
                        fa_start_in_box(1,1), fa_start_in_box(2,1), &
                        fa_start_in_box(3,1), 1, mdl%cell, mdl%fftbox)
                end do

                ! jd: Keep a copy in kets_in_ppds_on_a, as it is needed in
                !     pub_precond_recip below. Only the tck component (2) is needed.
                if(pub_precond_recip) then
                   call data_functions_copy(kets_in_ppds_on_a(ngwf_a,is,2), &
                        ket_on_bragrid_buffer(2))
                end if
             ! COARSE - COARSE - COARSE - COARSE - COARSE - COARSE - COARSE
             else
                ! jd: Copy ket in PPDs on A to a temp buffer for shaving
                do k_tck = 1, 2
                   call data_functions_copy(ket_on_bragrid_buffer(k_tck), &
                        kets_in_ppds_on_a(ngwf_a,is,k_tck))
                end do
             end if

             ! cks: shaving: zero points outside NGWF sphere in PPD rep.
             if (mdl%cell%n_pts > 1) then
                call basis_clean_function(ket_on_bragrid_buffer(1), &
                     ngwf_basis%spheres(local_aa_ngwf_idx), mdl%cell, mdl%fftbox, 1)
                call basis_clean_function(ket_on_bragrid_buffer(2), &
                     ngwf_basis%spheres(local_aa_ngwf_idx), mdl%cell, mdl%fftbox, 1)
             end if

             dot_npts = ngwf_basis%spheres(local_aa_ngwf_idx)%n_ppds_sphere * &
                  mdl%cell%n_pts

             factor = -mdl%cell%weight * pub_spin_fac * 4.0_DP * PI
             if(present(prefactor)) factor = factor * prefactor

             bra_start = ngwf_basis%spheres(local_aa_ngwf_idx)%offset
             do ipt = 0, dot_npts-1
                contra_grad%d(bra_start+ipt) = contra_grad%d(bra_start+ipt)- &   !jmecmplx
                     factor * ket_on_bragrid_buffer(1)%d(ipt+1)                  !jmecmplx
                if(.not. pub_precond_recip) then
                   cov_grad%d(bra_start+ipt) = cov_grad%d(bra_start+ipt)- &      !jmecmplx
                        factor * ket_on_bragrid_buffer(2)%d(ipt+1)               !jmecmplx
                end if
             end do

          end do loop_ngwf_a4b

       end do loop_spins_4b

#if 0
       print *,global_a
       print *, "-KIP_1:-============================@"
       print *, kets_in_ppds_on_a(1,1,1)%d
       print *, "-NGWF_CONTRA_1:-============================@"
       print *, contra_grad%d
       print *, "-NGWF_COV_1:-============================@"
       print *, cov_grad%d
#endif

       ! -------------------------------------------
       ! jd: Apply reciprocal-space preconditioning
       ! -------------------------------------------

       if(pub_precond_recip) then

          ! jd: Allocate workspace for ket in FFTbox around Aa
          allocate(kets_in_fftbox(ngwf_basis%num_on_atom(global_a), &
               pub_num_spins), stat=ierr)
          call utils_alloc_check(myself,'kets_in_fftbox',ierr)
          do is = 1, pub_num_spins
             do ngwf_a = 1, ngwf_basis%num_on_atom(global_a)
                call data_fftbox_alloc(kets_in_fftbox(ngwf_a, is), FFTBOX_DIMS,&
                     iscmplx=contra_grad%iscmplx)      !jmecmplx
             end do
          end do

          allocate(work_fftbox(FFTBOX_DIMS),stat=ierr)
          call utils_alloc_check(myself,'work_fftbox',ierr)
          allocate(zwork_fftbox(FFTBOX_DIMS),stat=ierr)
          call utils_alloc_check(myself,'zwork_fftbox',ierr)

          ! --------------------------------------------------------------------
          ! jd: for all a on A                                              aaa
          ! --------------------------------------------------------------------
          loop_ngwf_a4c:                                                      &
          do ngwf_a = 1, ngwf_basis%num_on_atom(global_a)

             global_aa_ngwf_idx = ngwf_basis%first_on_atom(global_a) + ngwf_a - 1
             local_aa_ngwf_idx = global_aa_ngwf_idx - &
                  ngwf_basis%first_on_proc(pub_my_proc_id)+1

             ! jd: Determine the position of the FFTbox around Aa and the
             !     position of Aa in the FFTbox
             call function_ops_batch_col_start(fa_box_start, fa_start_in_box, &
                  1, local_aa_ngwf_idx, local_aa_ngwf_idx, mdl%fftbox, &
                  mdl%cell, ngwf_basis)

             ! -----------------------------------------------------------------
             ! jd: for all spins                                            sss
             ! -----------------------------------------------------------------
             loop_spins_4c:                                                   &
             do is = 1, pub_num_spins
                ! jd: Transfer ket in PPDs on A to an FFTbox around Aa. Only the
                !     tcK component is needed (which_kernel==2).
                call basis_copy_function_to_box(&
                     kets_in_fftbox(ngwf_a,is), & ! target
                     fa_start_in_box(1,1), &
                     fa_start_in_box(2,1), &
                     fa_start_in_box(3,1), &
                     ngwf_basis%tight_boxes(local_aa_ngwf_idx), &
                     kets_in_ppds_on_a(ngwf_a,is,2), &
                     ngwf_basis%spheres(local_aa_ngwf_idx), mdl%cell, mdl%fftbox, 1)

#if 0
               call visual_output_dx(kets_in_fftbox(ngwf_a,is)%d, '', '', &
                    'grad_'//trim(utils_int_to_str(global_a))//'_'//&
                    trim(utils_int_to_str(ngwf_a))//'_'//&
                    trim(utils_int_to_str(is))//'.dx', &
                    0D0, 0D0, 0D0, fftbox%a1, fftbox%a2, fftbox%a3, &
                    fftbox%total_ld1, fftbox%total_ld2, fftbox%total_pt3, &
                    fftbox%total_ld1, .true., .true.)
#endif

             end do loop_spins_4c


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

          end do loop_ngwf_a4c

          deallocate(zwork_fftbox,stat=ierr)
          call utils_dealloc_check(myself,'zwork_fftbox',ierr)
          deallocate(work_fftbox,stat=ierr)
          call utils_dealloc_check(myself,'work_fftbox',ierr)

          ! --------------------------------------------------------------------
          ! jd: for all spins                                               sss
          ! --------------------------------------------------------------------
          loop_spins_4d:                                                      &
          do is = 1, pub_num_spins

             ! -----------------------------------------------------------------
             ! jd: for all a on A                                           aaa
             ! -----------------------------------------------------------------
             loop_ngwf_a4d:                                                   &
             do ngwf_a = 1, ngwf_basis%num_on_atom(global_a)

                global_aa_ngwf_idx = ngwf_basis%first_on_atom(global_a) + ngwf_a-1
                local_aa_ngwf_idx = global_aa_ngwf_idx - &
                     ngwf_basis%first_on_proc(pub_my_proc_id)+1

                call data_set_to_zero(ket_on_bragrid_buffer(2))
                ! cks: extract ppds belonging to bra function from ket fftbox
                call basis_extract_function_from_box(ket_on_bragrid_buffer(2), &
                     kets_in_fftbox(ngwf_a,is), &                 !jmecmplx
                     ngwf_basis%spheres(local_aa_ngwf_idx), &
                     ngwf_basis%tight_boxes(local_aa_ngwf_idx), &
                     fa_start_in_box(1,1), &
                     fa_start_in_box(2,1), &
                     fa_start_in_box(3,1), &
                     1, mdl%cell, mdl%fftbox)

                ! cks: shaving: zero points outside NGWF sphere in PPD rep.
                if (mdl%cell%n_pts > 1) then
                   call basis_clean_function(ket_on_bragrid_buffer(2), &
                        ngwf_basis%spheres(local_aa_ngwf_idx), mdl%cell, &
                        mdl%fftbox, 1)
                end if

                dot_npts = ngwf_basis%spheres(local_aa_ngwf_idx)%n_ppds_sphere*&
                     mdl%cell%n_pts

                bra_start = ngwf_basis%spheres(local_aa_ngwf_idx)%offset
                do ipt = 0, dot_npts-1
                   cov_grad%d(bra_start+ipt) = cov_grad%d(bra_start+ipt)- & ! sic, minus !jmecmplx
                        factor * ket_on_bragrid_buffer(2)%d(ipt+1)       !jmecmplx
                end do

             end do loop_ngwf_a4d

          end do loop_spins_4d

          do is = pub_num_spins, 1, -1
             do ngwf_a = ngwf_basis%num_on_atom(global_a), 1, -1
                call data_fftbox_dealloc(kets_in_fftbox(ngwf_a, is))
             end do
          end do
          deallocate(kets_in_fftbox, stat=ierr)
          call utils_dealloc_check(myself,'kets_in_fftbox',ierr)

       end if ! pub_precond_recip

#if 0
       print *,global_a
       print *, "-KIP_2:-============================@"
       print *, kets_in_ppds_on_a(1,1,1)%d
       print *, "-NGWF_CONTRA_2:-============================@"
       print *, contra_grad%d

       print *, "-NGWF_COV_2:-============================@"
       print *, cov_grad%d
#endif

       if(pub_pol_emb_dbl_grid) then
          deallocate(sum_b_uab_dkn_phi_bb_in_fftbox_on_a, stat=ierr)
          call utils_dealloc_check(myself,'sum_b_uab_dkn_phi_bb_in_fftbox_on_a', &
               ierr)
          deallocate(uab_dkn_phi_bb_in_fftbox_on_a, stat=ierr)
          call utils_dealloc_check(myself,'uab_dkn_phi_bb_in_fftbox_on_a',ierr)
          deallocate(dkn_phi_bb_in_fftbox_on_a, stat=ierr)
          call utils_dealloc_check(myself,'dkn_phi_bb_in_fftbox_on_a',ierr)
       end if

       deallocate(kets_in_ppds_on_a,stat=ierr)
       call utils_dealloc_check(myself,'kets_in_ppds_on_a',ierr)

       if(pub_dma_scale_charge) then
          deallocate(tau1term_in_ppds_on_a, stat=ierr)
          call utils_dealloc_check(myself, 'tau1term_in_ppds_on_a', ierr)
          deallocate(tau0term_in_ppds_on_a, stat=ierr)
          call utils_dealloc_check(myself, 'tau0term_in_ppds_on_a', ierr)
       end if

    end do loop_A4
    ! --------------------------------------------------------------------------

    call hash_table_free(qterms_ht)
    call hash_table_free(jqterms_ht)

    ! jd: Free multipole sets allocated in pol_emb_read_from_file()
    do mpole_set = 1, n_sets
       call multipole_free_spherical_set(mm_multipoles(mpole_set))
    end do

    deallocate(jqterm, stat=ierr)
    call utils_dealloc_check(myself,'jqterm',ierr)

    do k_tck = 1, 2
       call data_functions_dealloc(ket_on_bragrid_buffer(k_tck))
    end do

    call timer_clock(myself//'_imbalance',1)

    if(pub_dma_scale_charge) then
       call utils_abort(myself//': Charge-scaling not implemented in NGWF grad')
#if 0
       call dma_free_multipoles(scf_dma_multipoles)
#endif
    end if

    call comms_barrier

    call timer_clock(myself//'_imbalance',2)

999 continue

    ! jd: Clean up what pol_emb_read_from_file() set up
    deallocate(dummy_energy_term_values, stat=ierr)
    call utils_dealloc_check(myself, 'dummy_energy_term_values', ierr)
    deallocate(dummy_energy_term_names, stat=ierr)
    call utils_dealloc_check(myself, 'dummy_energy_term_names', ierr)
    deallocate(mm_multipoles,stat=ierr)
    call utils_dealloc_check(myself,'mm_multipoles',ierr)

    call utils_sanity_check(contra_grad%d,'contra_grad on exit from '//myself, &
         excessive=1D20)
    call utils_sanity_check(cov_grad%d,'contra_grad on exit from '//myself, &
         excessive=1D20)

    call timer_clock(myself,2)

  end subroutine pol_emb_ngwf_grad_contribution

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function pol_emb_mm_species_to_rep_param(iparam, mm_species)
    !==========================================================================!
    ! Returns the value of a given MM repulsive potential parameter for a      !
    ! MM species. The value is found in a table prepared in advance by         !
    ! rundat_blocks.                                                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2018.                                 !
    !==========================================================================!

    use constants, only: DP
    use utils, only: utils_int_to_str, utils_assert, utils_abort

    implicit none

    ! jd: Arguments
    integer, intent(in)          :: iparam
    character(len=4), intent(in) :: mm_species

    ! jd: Local variables
    integer :: ispecies
    character(len=*), parameter :: myself = 'pol_emb_mm_species_to_rep_param'

    ! -------------------------------------------------------------------------

    call utils_assert(allocated(mm_all_species) .and. allocated(mm_rep_params),&
         myself//': At least one of your MM multipole sets is tagged "rep", but&
         & %block mm_rep_params was absent. You need to specify MM repulsive&
         & potential parameters for every MM species in this block.')

    call utils_assert(iparam > 0 .and. &
         iparam <= ubound(mm_rep_params,1), &
         myself//': iparam out of range: '//trim(utils_int_to_str(iparam)))

    do ispecies = 1, ubound(mm_all_species,1)
       if(trim(mm_all_species(ispecies)) == trim(mm_species)) then
          pol_emb_mm_species_to_rep_param = mm_rep_params(iparam,ispecies)
          return
       end if
    end do

    ! jd: Not found in list
    call utils_abort(myself//': species "'//trim(mm_species)//'" not found in &
         &mm_all_species. You need to specify MM repulsive potential parameters&
         & for every MM species in %block mm_rep_params.')

  end function pol_emb_mm_species_to_rep_param

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module polarisable_embedding
