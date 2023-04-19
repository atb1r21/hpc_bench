! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!-----------------------------------------------------------------------!
!                              D   M   A                                !
!-----------------------------------------------------------------------!
!                                                                       !
! Distributed multipole analysis module.                                !
!                                                                       !
! The current version of this module (using SW_EX and SW_RI) was written!
! by Jacek Dziedzic in October-November 2016.                           !
! Generalised to expand either via NGWF pairs or via atom-pair densities!
! by Jacek Dziedzic in June 2017.                                       !
!-----------------------------------------------------------------------!

!@TODO: - Periodicity does not honour global settings, OBCs are assumed
!         everywhere. Look for @MIC.
!       - Identical radii of all NGWFs are assumed, or things get complicated.
!         Look for @IDENTICAL_RADII

module dma

  use constants, only: DP, CRLF

  implicit none

  private

  public :: dma_expand_ngwf_pairs
  public :: dma_expand_pair_densities
  public :: dma_calculate
  public :: dma_output_potential
  public :: dma_write_multipoles
  public :: dma_free
  public :: dma_free_multipoles
  public :: dma_free_swexes

  integer, parameter, public :: N_DMA_FILES = 8

contains

  ! ---------------------------------------------------------------------------
  ! ----- P U B L I C   S U B R O U T I N E S   A N D   F U N C T I O N S -----
  ! ---------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dma_expand_ngwf_pairs(rep, ngwf_basis, mdl, ireg, main_swex_h)
    !==========================================================================!
    ! Performs an expansion of NGWF pairs into spherical waves. The resultant  !
    ! expansion coefficients are stored in rep's dma_swexes(:). If polarisable !
    ! embedding with QM* representation is in use, a second set of coefficients!
    ! needed for LNV gradients ('dcoeffs') is also computed and stored.        !
    !                                                                          !
    ! Both sets of coefficients are communicated across procs so that all      !
    ! S-neighbours of each pair have access to that pair's coefficients. This  !
    ! is needed in the DMA calculation proper later on. S-neighbourhood is     !
    ! determined from the s_atoms_nl neighbour list of the parent SWRI.        !
    !                                                                          !
    ! This routine is a wrapper around swx_expand_ngwf_pairs(), taking care of !
    ! initialisation of containers, sanity checks, Bessel averaging loop,      !
    ! and the call to a slave routine for calculating dcoeffs.                 !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   rep (inout): an NGWF_REP that is expanded.                             !
    !   ngwf_basis (in): NGWF basis corresponding to rep.                      !
    !   mdl (in): Needed for the ionic positions and cell.                     !
    !   main_swex_h (in): Index to rep%swexes(:) where the expansion will be   !
    !                     stored. Note that if Bessel averaging is used, also  !
    !                     the next index will be written to.                   !
    !   ireg (in):        Index for the region of this NGWF basis.             !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Even if charge-scaling is in effect, neither set of coefficients is    !
    !   scaled! This is because the scaling factor is unknown at the time --   !
    !   -- after all it depends on DKN. During the actual calculation of the   !
    !   multipoles in dma_calculate() the scaling factors are determined, and  !
    !   the monopoles *are* scaled. Still, the coefficients in rep%swexes(:)   !
    !   remain unscaled.                                                       !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October 2016.                               !
    !==========================================================================!

    use bibliography, only: bibliography_cite
    use comms, only: pub_on_root
    use constants, only: SW_O, SW_V, stdout, VERBOSE, METRIC_OVERLAP, &
         REP_SWEX_POL_EMB_DMA_1, REP_SWEX_POL_EMB_AUX_DMA_1
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_dma_metric, pub_dma_use_ri, pub_dma_bessel_averaging,&
         pub_devel_code, pub_output_detail, pub_dma_max_q, &
         pub_num_kpoints, PUB_1K
    use sw_expansion, only: swx_expand_local_ngwf_pairs
    use sw_resolution_of_identity, only: swri_get_handle_to, swri_library
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_flush, utils_devel_code

    implicit none

    ! jd: Arguments
    type(NGWF_REP), intent(inout)              :: rep
    type(FUNC_BASIS), intent(in)               :: ngwf_basis
    type(MODEL), intent(in)                    :: mdl
    integer, intent(in)                        :: ireg
    integer, intent(in)                        :: main_swex_h

    ! jd: --- Local variables ---

    ! jd: Expansion business
    integer     :: dma_swri_h  ! jd: Handle to container for SW res-of-id
    integer     :: swex_h
    integer     :: dma_run, n_dma_runs

    ! jd: Atom index book-keeping
    integer :: global_b, orig_global_b
    integer :: n_atoms_in_dma

    ! jd: SW index book-keeping
    integer :: max_q
    integer :: min_l, max_l

    ! jd: Varia
    logical       :: debug_show_skipped
    character(len=*), parameter :: myself = 'dma_expand_ngwf_pairs'
    ! agrecocmplx
    logical :: loc_cmplx

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)
    call bibliography_cite('ANHARM_AND_DMA')

    debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
         pub_devel_code, no_bcast = .true., no_warn=.true.)

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(ireg)%iscmplx

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    !   ireg (in):        Index for the region of this NGWF basis.             !
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine '//myself//' not ready yet for more than one k-point.')

    ! agrecokpt: TO-DO: check/implement k-points here as well
    call utils_assert(loc_cmplx.eqv..false., &
         'Subroutine '//myself//' not ready yet for complex NGWFs.')

    ! jd: Locate SW_RI in the library
    dma_swri_h = swri_get_handle_to(pub_dma_use_ri)

    n_atoms_in_dma = count(mdl%elements(:)%in_swri(dma_swri_h))
    call utils_assert(n_atoms_in_dma > 0, &
         myself//': Cannot do DMA with zero atoms.')

    ! jd: Verbose output for subsystem DMA, if devel code
    if(debug_show_skipped .and. pub_on_root) then
       write(stdout,'(a)') 'Atoms considered for DMA:'
       write(stdout,'(a)') ' considered? |  atom #  | orig atom # | species'
       ! rc2013: this forces 1 region
       do global_b = 1, mdl%regions(ireg)%par%nat
          orig_global_b = mdl%regions(ireg)%par%orig_atom(global_b)
          write(stdout,'(a,i14,i15,a,a)') &
               merge('     yes','      no', &
               mdl%regions(ireg)%elements(orig_global_b)%in_swri(dma_swri_h)), &
               global_b, orig_global_b, '   ', &
               trim(mdl%regions(ireg)%elements(orig_global_b)%species_id)
       end do
    end if

    if(pub_dma_bessel_averaging) then
       n_dma_runs = 2
    else
       n_dma_runs = 1
    end if

    max_q = pub_dma_max_q
    ! -------------------------------------------------------------------------
    ! Loop over DMA runs (for Bessel averaging)
    ! -------------------------------------------------------------------------
    do dma_run = 1, n_dma_runs

       ! jd: Current SW_EX
       swex_h = main_swex_h + dma_run - 1
       min_l = rep%swexes(swex_h)%quality%min_l
       max_l = rep%swexes(swex_h)%quality%max_l

       ! { --------------------------------------------------------------------
       ! { Actual DMA expansion -- runs 1 and 2.
       ! { --------------------------------------------------------------------

       if (pub_on_root .and. pub_output_detail >= VERBOSE) then
          write(stdout,'(a,i0,a,i0,a,i0,a,i0,a)') &
               'DMA: Performing SW expansion of NGWF pairs (multipole orders ',&
               min_l, '-', max_l, ', ',max_q, &
               ' Bessels) for ', n_atoms_in_dma, ' atoms... '
       end if

       ! jd: Perform SW expansion of NGWF products
       call swx_expand_local_ngwf_pairs(rep, ireg, swri_library(dma_swri_h), &
            swex_h, ngwf_basis, mdl%cell, mdl%fftbox, mdl%regions(ireg)%par, &
            mdl%regions(ireg)%elements, &
            merge(SW_O,SW_V,pub_dma_metric == METRIC_OVERLAP), &
            swri_library(dma_swri_h)%s_atoms_nl)

       ! jd: Ensure all output is flushed and does not interfere with QC output
       !     @Remove this once debug statements from non-root ranks are removed
       call utils_flush(stdout)

       ! jd: *************************************************************
       !     *** Calculate the D coefficients if polarisable embedding ***
       !     *************************************************************
       if(main_swex_h == REP_SWEX_POL_EMB_DMA_1 .or. &
            main_swex_h == REP_SWEX_POL_EMB_AUX_DMA_1) then
          call dma_calculate_dcoeffs_ngwf_pairs(rep%swexes(swex_h), &
               ngwf_basis, mdl, dma_swri_h, &
               merge(SW_O,SW_V,pub_dma_metric == METRIC_OVERLAP), swex_h)
       end if

       if (pub_on_root .and. pub_output_detail == VERBOSE) then
          write(stdout,'(a)') '... done.'
       end if

       max_q = max_q - 1

    end do ! dma runs

    call timer_clock(myself,2)

  end subroutine dma_expand_ngwf_pairs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dma_expand_pair_densities(coeffs_ht, rep, ireg, denskern, &
       ngwf_basis, mdl, main_swex_h, vacuum_rep, vacuum_denskern)
    !==========================================================================!
    ! Performs an expansion of atom-pairs densities into SWs. In contrast to   !
    ! dma_expand_ngwf_pairs(), the resultant expansion coefficients are stored !
    ! in coeffs_ht (*not* in rep!), and D coefficients are *not* calculated.   !
    !                                                                          !
    ! The coefficients are communicated across procs so that all S-neighbours  !
    ! of each pair have access to that pair's coefficients. This is needed in  !
    ! the DMA calculation proper later on. S-neighbourhood is determined from  !
    ! the s_atoms_nl neighbour list of the parent SWRI.                        !
    !                                                                          !
    ! This routine is a wrapper around swx_expand_local_pair_densities(),      !
    ! taking care of initialisation of containers, sanity checks and the       !
    ! Bessel averaging loop.                                                   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   rep (inout): an NGWF_REP that is expanded.                             !
    !   ngwf_basis (in): NGWF basis corresponding to rep.                      !
    !   mdl (in): Needed for the ionic positions and cell.                     !
    !   main_swex_h (in): Index to rep%swexes(:) for getting the quality of the!
    !                     expansion. Note that if Bessel averaging is used,    !
    !                     also the next index in rep%swexes will be accessed.  !
    !   ireg (in):        Index for the region of this NGWF basis.             !
    !--------------------------------------------------------------------------!
    ! Expert mode:                                                             !
    !   When arguments 'vacuum_rep' and 'vacuum denskern' are specified, the   !
    !   expanded quantity will be the density difference, ie.                  !
    !   \phi_Aa K^Aa,Bb \phi_Bb - \phi0_Aa K0^Aa,Bb \phi_Bb, where \phi0 refer !
    !   to vacuum_rep, and K0 refers to vacuum_denskern.                       !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Even if charge-scaling is in effect, the coefficients are *not* scaled.!
    !   This is a vestige of the NGWF-product fitting, where the scaling factor!
    !   would not be known at this stage (as it depends on DKN).               !
    !   During the actual calculation of the multipoles in dma_calculate() the !
    !   scaling factors are determined, and the monopoles *are* scaled. Still, !
    !   the coefficients in coeffs_ht remain unscaled.                         !
    !   Also, charge-scaling has not been properly tested in density fitting.  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in June 2017.                                  !
    ! Modified for embedding by Joseph Prentice, September 2018                !
    !==========================================================================!

    use bibliography, only: bibliography_cite
    use comms, only: pub_on_root
    use constants, only: SW_O, SW_V, stdout, VERBOSE, METRIC_OVERLAP, &
         REP_SWEX_POL_EMB_DMA_1
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_dma_metric, pub_dma_use_ri, pub_dma_bessel_averaging,&
         pub_devel_code, pub_output_detail, pub_dma_max_q, pub_num_spins, &
         pub_num_kpoints, PUB_1K
    use sparse, only: SPAM3, sparse_get_par
    use sw_expansion, only: swx_expand_local_pair_densities
    use sw_resolution_of_identity, only: swri_get_handle_to, swri_library
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_flush, utils_devel_code

    implicit none

    ! jd: Arguments
    type(HT_HASH_TABLE), intent(inout), target :: coeffs_ht
    type(NGWF_REP), intent(in)                 :: rep
    type(SPAM3), intent(in)                    :: denskern(pub_num_spins)
    type(FUNC_BASIS), intent(in)               :: ngwf_basis
    type(MODEL), intent(in)                    :: mdl
    integer, intent(in)                        :: main_swex_h
    integer, intent(in)                        :: ireg
    type(NGWF_REP), intent(inout), optional    :: vacuum_rep
    type(SPAM3), intent(in), optional          :: vacuum_denskern(pub_num_spins)

    ! jd: --- Local variables ---

    ! jd: Expansion business
    integer     :: dma_swri_h  ! jd: Handle to container for SW res-of-id
    integer     :: swex_h
    integer     :: dma_run, n_dma_runs

    ! jd: Atom index book-keeping
    integer :: global_b, orig_global_b
    integer :: n_atoms_in_dma

    ! jd: SW index book-keeping
    integer :: max_q
    integer :: min_l, max_l

    ! jd: Varia
    logical :: debug_show_skipped
    character(len=*), parameter :: myself = 'dma_expand_densities'
    ! agrecocmplx
    logical :: loc_cmplx
    type(PARAL_INFO), pointer :: par

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)
    call bibliography_cite('ANHARM_AND_DMA')

    debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
         pub_devel_code, no_bcast = .true., no_warn=.true.)

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(ireg)%iscmplx

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine '//myself//' not ready yet for more than one k-point.')

    ! agrecokpt: TO-DO: check/implement k-points here as well
    call utils_assert(loc_cmplx.eqv..false., &
         'Subroutine '//myself//' not ready yet for complex NGWFs.')

    call utils_assert(present(vacuum_rep) .eqv. present(vacuum_denskern), &
         myself//'Optional arguments must all be specified or all be absent')

    ! rc2013: get the parallel strategy from the sparse matrix
    call sparse_get_par(par, denskern(1))

    ! jd: Locate SW_RI in the library
    dma_swri_h = swri_get_handle_to(pub_dma_use_ri)

    n_atoms_in_dma = count(mdl%regions(ireg)%elements(:)%in_swri(dma_swri_h))
    call utils_assert(n_atoms_in_dma > 0, &
         myself//': Cannot do DMA with zero atoms.')

    ! jd: Verbose output for subsystem DMA, if devel code
    if(debug_show_skipped .and. pub_on_root) then
       write(stdout,'(a)') 'Atoms considered for DMA:'
       write(stdout,'(a)') ' considered? |  atom #  | orig atom # | species'
       do global_b = 1, par%nat
          orig_global_b = par%orig_atom(global_b)
          write(stdout,'(a,i14,i15,a,a)') &
               merge('     yes','      no', &
               mdl%regions(ireg)%elements(orig_global_b)%in_swri(dma_swri_h)), &
               global_b, orig_global_b, '   ', &
               trim(mdl%regions(ireg)%elements(orig_global_b)%species_id)
       end do
    end if

    if(pub_dma_bessel_averaging) then
       n_dma_runs = 2
    else
       n_dma_runs = 1
    end if

    max_q = pub_dma_max_q
    ! -------------------------------------------------------------------------
    ! Loop over DMA runs (for Bessel averaging)
    ! -------------------------------------------------------------------------
    do dma_run = 1, n_dma_runs

       ! jd: Current SW_EX
       swex_h = main_swex_h + dma_run - 1
       min_l = rep%swexes(swex_h)%quality%min_l
       max_l = rep%swexes(swex_h)%quality%max_l

       ! { --------------------------------------------------------------------
       ! { Actual DMA expansion -- runs 1 and 2.
       ! { --------------------------------------------------------------------

       if (pub_on_root .and. pub_output_detail >= VERBOSE) then
          write(stdout,'(a,i0,a,i0,a,i0,a,i0,a)') &
               'DMA: Performing SW expansion of pair densities (multipole orders ',&
               min_l, '-', max_l, ', ',max_q, &
               ' Bessels) for ', n_atoms_in_dma, ' atoms... '
       end if

       call swx_expand_local_pair_densities(coeffs_ht, rep, ireg, &
            swri_library(dma_swri_h), swex_h, denskern, ngwf_basis, mdl%cell, &
            mdl%fftbox, mdl%regions(ireg)%elements, &
            merge(SW_O,SW_V,pub_dma_metric == METRIC_OVERLAP), &
            swri_library(dma_swri_h)%s_atoms_nl, vacuum_rep, vacuum_denskern)

       ! jd: Ensure all output is flushed and does not interfere with QC output
       !     @Remove this once debug statements from non-root ranks are removed
       call utils_flush(stdout)

       ! *******************************************************************
       ! jd: NB: D coefficients are not calculated here, as in D mode they
       !         require a DKN-independent, ngwf-pair fit.
       ! *******************************************************************

       if (pub_on_root .and. pub_output_detail == VERBOSE) then
          write(stdout,'(a)') '... done.'
       end if

       max_q = max_q - 1

    end do ! dma runs

    call timer_clock(myself,2)

  end subroutine dma_expand_pair_densities

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dma_calculate(my_multipoles, &                            ! out
       output_name, dma_swexes, &                                      ! in
       n_occ, overlap, ngwf_basis, mdl, swex_h, &                      ! in
       mode, denskern, coeffs_ht, dkn_scale_factor)                    ! in, opt
    !==========================================================================!
    ! Computes a Distributed Multipole Analysis of the density.                !
    ! Multipoles for each centre are determined either from pre-expanded NGWF  !
    ! pairs obtained earlier by dma_expand_ngwf_pairs() and the current kernel,!
    ! known as mode 'P', or from pre-expanded atom-pair densities obtained     !
    ! earlier by dma_expand_pair_densities(), this is mode 'D'. For mode 'D'   !
    ! do not pass a density kernel.                                            !
    !                                                                          !
    ! Multipoles are output to files determined by 'output_name'.              !
    !                                                                          !
    ! This subroutine also works with 'subsystem DMA', ie. in the case where   !
    ! only some atoms belong to the DMA SWRI. This becomes slightly tricky when!
    ! charge scaling is in effect. Either a desired number of electrons in the !
    ! subsystem must be specified, or Mulliken analysis is performed on the fly!
    ! and the desired charge (a 'partial trace' of KS) is obtained from it.    !
    !--------------------------------------------------------------------------!
    ! my_multipoles (out): Contains the calculated multipoles as spherical     !
    !                     multipole sets. If there is no Bessel averaging, only!
    !                     element 1 is meaningful. If there is Bessel avera-   !
    !                     ging, the first element contains the multipoles for  !
    !                     the nominal max_q, the second -- for (max_q-1), and  !
    !                     element 3 -- the averages. Averaging is performed    !
    !                     after possible scaling of monopoles.                 !
    !                     Only stores multipoles for atoms local to this proc. !
    ! output_name (in): Added as an infix to obtain the output filename.       !
    ! dma_swexes (in): SW expansions of NGWF pairs.                            !
    ! n_occ, ngwf_basis, mdl and elements (in) -- obvious.                     !
    ! overlap (in): Only needed for Mulliken analysis (when doing subsystem DMA!
    !               with auto-determination of target charge for subsystem).   !
    ! swex_h (in): User-defined identifier used to distinguish coeffs from     !
    !              different expansions (when Bessel averaging).               !
    ! mode (in): 'P' for calculating multipoles from NGWF-pair expansions      !
    !            precomputed in dma_swexes(:);                                 !
    !            'D' for calculating multipoles from pair-densities precomputed!
    !            in coeffs_ht.                                                 !
    ! denskern (in, opt): Density kernel. Needed only in mode 'P'.             !
    ! coeffs_ht (in, opt): Precomputed expansion coefficients. Needed only in  !
    !                      mode 'P'.                                           !
    ! dkn_scale_factor(input, opt): The kernel is scaled differently between   !
    !                               properties and hamiltonian. This scalar    !
    !                               factor can be used to make DMA work in both!
    !                               contexts.                                  !
    !--------------------------------------------------------------------------!
    ! Current version of this subroutine is due to Jacek Dziedzic.             !
    ! Cleaned up by Jacek Dziedzic in May 2014.                                !
    ! Generalised for max_q averaging by Jacek Dziedzic in February 2015.      !
    ! Extended to calculate derivative terms ('dcoeffs') by Jacek Dziedzic in  !
    ! December 2015.                                                           !
    ! Carved out the NGWF expansion business to a separate routine, cleaned up !
    ! -- Jacek Dziedzic, October 2016.                                         !
    ! Extended with Mulliken analysis for subsystem DMA when used in gradients,!
    ! -- Jacek Dziedzic, May 2017.                                             !
    ! Extended to support two modes (P, D) by Jacek Dziedzic in June 2017.     !
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_on_root, pub_root_proc_id, comms_send,&
         comms_recv, comms_wait, comms_reduce, pub_any_source, comms_bcast
    use constants, only: stdout, VERBOSE, CRLF, SW_O, SW_V, METRIC_OVERLAP, &
         SAFE_DIV_EPS, SWEX_BB_CC_SYMMETRIC, garbage_real
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE, hash_table_lookup_nocount
    use model_type, only: MODEL
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET, multipole_set_gather, &
         multipole_set_scale_monopoles
    use parallel_strategy, only: PARAL_INFO
    use population, only: population_analysis_mulliken
    use rundat, only: pub_dma_max_q, &
         pub_dma_scale_charge, pub_dma_use_ri, pub_dma_bessel_averaging, &
         pub_num_spins, pub_devel_code, pub_output_detail, pub_spin_fac, &
         pub_dma_target_num_val_elec, pub_dma_output_potential, &
         pub_dma_output_potential_reference, pub_dma_multipole_scaling, &
         pub_dma_metric, pub_pol_emb_qmstar, PUB_1K
    use sparse, only: sparse_get_element, SPAM3, sparse_get_par
    use sparse_array, only: SPAM3_ARRAY
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY
    use sw_expansion_type, only: SW_EX
    use sw_resolution_of_identity, only: swri_get_handle_to, swri_library
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_abort, utils_alloc_check, utils_unit, &
         utils_flushed_string_output, utils_int_to_str, utils_devel_code, &
         utils_dealloc_check, utils_real_to_str

    implicit none

    ! jd: Arguments
    !     my_multipoles and all_multipoles are indexed by the number of DMA run
    !     (for Bessel averaging): 1 or 2. Entry 3 holds the average of the two
    !     runs. My_multipoles contains multipoles for MPI-rank-local atoms,
    !     all_multipoles contains all multipoles (on root) is dummy (elsewhere).
    type(SPHERICAL_MULTIPOLE_SET), intent(out)        :: my_multipoles(3)
    character(len=*)                                  :: output_name
    type(SW_EX), intent(in)                           :: dma_swexes(:)
    integer, intent(in)                               :: n_occ(pub_num_spins)
    type(SPAM3_EMBED), intent(in)                     :: overlap
    type(FUNC_BASIS), intent(in)                      :: ngwf_basis(1)
    type(MODEL), intent(in)                           :: mdl
    integer, intent(in)                               :: swex_h
    character, intent(in)                             :: mode
    type(SPAM3_EMBED_ARRAY), intent(in), optional     :: denskern
    type(HT_HASH_TABLE), intent(in), optional, target :: coeffs_ht
    real(kind=DP), intent(in), optional               :: dkn_scale_factor

    ! jd: --- Local variables ---
    type(SPHERICAL_MULTIPOLE_SET)                     :: all_multipoles(3)

    ! jd: Expansion business
    integer     :: dma_swri_h  ! jd: Handle to container for SW res-of-id
    integer     :: num_sws_per_centre
    integer     :: max_q
    integer     :: min_l, max_l
    integer     :: dma_run, n_dma_runs

    ! jd: Atom number book-keeping
    integer :: local_b
    integer :: global_b, global_c
    integer :: orig_global_b
    integer :: c_idx

    ! jd: NGWF number book-keeping
    integer :: global_bb_ngwf_idx, global_cc_ngwf_idx
    integer :: first_ngwf_idx_of_b, first_ngwf_idx_of_c
    integer :: ngwf_b, ngwf_c

    ! jd: Subsequent levels of coefficients
    real(kind=DP), allocatable  :: sphw_coeffs(:,:) ! 1st idx: spin
    real(kind=DP), allocatable  :: c_coeffs(:,:)    ! 1st idx: spin
    real(kind=DP), allocatable  :: g_coeffs(:)
    integer                     :: coeffs_kind

    ! jd: All multipoles (on root only) for printing to file
    real(kind=DP)               :: total_monopole_charge ! sum of l=0 components

    ! jd: Comms-related
    integer :: n_multipoles_per_atom

    ! jd: Kernel element scale factor, to tackle inconsistency of DKN scaling
    !     between properties and calls to hamiltonian_lhxc_calculate.
    real(kind=DP) :: local_dkn_scale_factor

    ! jd: Varia
    integer       :: file_units(N_DMA_FILES)
    integer       :: is
    integer       :: coeffs_is
    integer       :: loc_num_spins
    integer       :: cur_swex_h
    real(kind=DP), allocatable :: q_atom(:,:) ! Mulliken pop analysis
    real(kind=DP) :: total_core_charge, total_wanted_charge, total_dma_charge
    real(kind=DP) :: dkn_el
    integer       :: num_of_sphw_coeffs
    integer       :: n_atoms_in_dma
    real(kind=DP) :: scaling_factor
    ! jd: Don't change these to integers or truncation issues will be upon you
    real(kind=DP) :: nelecs, nelecs_inactive, nelecs_total, nelecs_expected
    real(kind=DP), parameter :: nelecs_fraction_threshold = 0.04_DP
    integer       :: ierr
    character(len=32)           :: swex_name
    character(len=128)          :: output_type
    logical                     :: debug_show_skipped
    character(len=*), parameter :: myself = 'dma_calculate'
    type(PARAL_INFO), pointer   :: par

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    call utils_assert(mode == 'P' .or. mode == 'D', myself//': Invalid mode', &
         mode)

    call utils_assert((mode == 'P' .and. present(denskern) .and. &
         .not. present(coeffs_ht)) .or. &
         (mode == 'D' .and. .not. present(denskern) .and. &
         present(coeffs_ht)), myself//&
              ": Mode 'P' requires passing the denskern argument and &
              &not passing the coeffs_ht argument. Mode 'D' requires *not* &
              &passing denskern, and passing coeffs_ht.", mode)

    debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
         pub_devel_code, no_warn=.true.)

    ! rc2013: get the parallel strategy from the sparse matrix
    call sparse_get_par(par, overlap%p)

    if(present(dkn_scale_factor)) then
       local_dkn_scale_factor = dkn_scale_factor
    else
       local_dkn_scale_factor = 0.5_DP
    end if

    ! jd: Obtain swex name without trailing '-1' or '-2'
    swex_name = dma_swexes(1)%swex_name(1:index(dma_swexes(1)%swex_name,'-')-1)

    ! jd: Create and open output files for printing out the multipoles
    call dma_write_multipoles_open_outputs(output_name, file_units)

    ! jd: Locate SW_RI in the library
    dma_swri_h = swri_get_handle_to(pub_dma_use_ri)

    n_atoms_in_dma = count(mdl%elements(:)%in_swri(dma_swri_h))
    call utils_assert(n_atoms_in_dma > 0, &
         myself//': Cannot do DMA with zero atoms.')

    ! jd: Determine the total ionic charge of the subsystem
    total_core_charge = sum(mdl%elements(:)%ion_charge, &
         mask = mdl%elements(:)%in_swri(dma_swri_h))

    ! jd: Determine the total electronic charge of the subsystem
    if(n_atoms_in_dma == par%nat) then
       ! jd: DMA for entire system
       nelecs = pub_spin_fac * real(sum(n_occ),kind=DP) * &
            pub_dma_multipole_scaling
    else
       ! jd: DMA for a subsystem
       if(pub_dma_target_num_val_elec /= -999999.0_DP) then
          ! jd: User-specified number of valence electrons in effect
          nelecs = pub_dma_target_num_val_elec * pub_dma_multipole_scaling
          call utils_assert(.not. pub_pol_emb_qmstar, &
               "Cannot combine (1) polarisable embedding using the QM* represen&
               &tation with (2) subsystem DMA that declares 'dma_target_num_val&
               &_elec'. (1) implies DMA analysis will be used in LNV and NGWF&
               &gradients and energy expressions. (2) implies approximate treat&
               &ment, since Tr[KS] (used in charge-scaling expressions) over &
               &the subsystem will not be constant and so will not be strictly &
               &equal to 'dma_target_num_val_elec'. Either do DMA on the entire&
               &system, or drop 'dma_target_num_val_elec' (Mulliken analysis&
               &will then be used to determine the target number of electrons&
               &in the subsystem, or do not do (1).")
       else
          ! jd: Number of valence electrons from Mulliken analysis for the
          !     active subset of atoms
          nelecs = 0D0
          nelecs_inactive = 0D0
          call utils_assert(mode == 'P', myself//": mode 'D' incompatible with &
               &Mulliken analysis for active subset -- kernel and overlap needed")
          call population_analysis_mulliken(overlap, denskern%m(:,PUB_1K), &
               mdl, ngwf_basis, -1, q_atom)

          if(pub_on_root) then
             do global_b = 1, par%nat
                orig_global_b = par%orig_atom(global_b)
                if(mdl%elements(orig_global_b)%in_swri(dma_swri_h)) then
                   nelecs = nelecs + sum(q_atom(orig_global_b,1:pub_num_spins))
                else
                   nelecs_inactive = nelecs_inactive + &
                        sum(q_atom(orig_global_b,1:pub_num_spins))
                end if
             end do
             nelecs_total = nelecs + nelecs_inactive
             nelecs_expected = pub_spin_fac * real(sum(n_occ),kind=DP)
             ! jd: Work around the fact that ngwf_gradient_lnv() does not
             !     account for spin degeneracy correctly.
             if(abs(nelecs_total-nelecs_expected)/nelecs_expected > &
                  nelecs_fraction_threshold) then
                ! jd: Number of electrons far from expected total
                if(abs(2D0*nelecs_total-nelecs_expected)/nelecs_expected > &
                     nelecs_fraction_threshold) then
                   ! jd: ... and far from half of it, something's off
                   call utils_abort(myself//": The total # of electrons ("//&
                        trim(utils_real_to_str(nelecs_total))//&
                        ") is suspicious. It's neither close enough to the &
                        &expected number ("//&
                        trim(utils_real_to_str(nelecs_expected))//&
                        "), nor to half of it (which happens when called from &
                        &ngwf_gradient_lnv() due to differing spin normalisatio&
                        &n conventions). There were "//&
                        trim(utils_real_to_str(nelecs))//&
                        " electrons in the DMA region, and "//&
                        trim(utils_real_to_str(nelecs_inactive))//&
                        " electrons outside it.")
                else
                   ! jd: ... but close to half of it, we must have been called
                   !     from ngwf_gradient_lnv().
                   nelecs = 2D0*nelecs
                   nelecs_inactive = 2D0*nelecs_inactive
                   nelecs_expected = 2D0*nelecs_expected
                   nelecs_total = 2D0*nelecs_total
                   write(stdout,'(a)') CRLF//'DMA: Looks like a gradient step, &
                        &fixing spin degeneracy normalisation.'
                end if
             else
                ! jd: No-op, number of electrons close to expected.
             end if

             nelecs = nelecs * pub_dma_multipole_scaling
             nelecs_inactive = nelecs_inactive * pub_dma_multipole_scaling
             write(stdout,'(a,f20.8)') &
                  'DMA: Target number of electrons in DMA region:      ', nelecs
             write(stdout,'(a,f20.8)') &
                  'DMA: Target number of electrons outside DMA region: ', &
                  nelecs_inactive
             write(stdout,'(a,f20.8)') &
                  'DMA: Total number of electrons accounted for:       ', &
                  nelecs_total
          end if

          call comms_bcast(pub_root_proc_id, nelecs)
          call comms_bcast(pub_root_proc_id, nelecs_inactive)
          call comms_bcast(pub_root_proc_id, nelecs_total)

          deallocate(q_atom, stat=ierr)
          call utils_dealloc_check(myself,'q_atom',ierr)

       end if
    end if
    if(pub_dma_scale_charge) then
       call utils_assert(nelecs/=0D0, myself//&
            ': Cannot have zero electrons as target for charge scaling')
    end if

    total_wanted_charge = -total_core_charge + nelecs

    ! jd: Warn against outputting DMA potential and reference potential if
    !     DMA potential only includes a subset of atoms
    if(n_atoms_in_dma /= par%nat .and. pub_dma_output_potential .and. &
         pub_dma_output_potential_reference) then
       call utils_flushed_string_output('WARNING: DMA in ['//&
            trim(swri_library(dma_swri_h)%swri_name)//&
            '] is performed for a *subsystem*. The multipole'//CRLF//&
            '         potential will &
            &be due to this subsystem, while the reference potential'//CRLF//&
            '         will be &
            &due to the entire system. Do not expect these to match!'//CRLF)
    end if

    if(pub_dma_bessel_averaging) then
       n_dma_runs = 3
    else
       n_dma_runs = 1
    end if

    max_q = pub_dma_max_q
    min_l = dma_swexes(1)%quality%min_l
    max_l = dma_swexes(1)%quality%max_l

    ! jd: Initialise multipole datastructure
    call init_multipoles(my_multipoles, mdl%elements, min_l, max_l, par)

    ! -------------------------------------------------------------------------
    ! Loop over DMA runs (for Bessel averaging)
    ! -------------------------------------------------------------------------
    do dma_run = 1, n_dma_runs

       if(dma_run < 3) then
          output_type = 'Values for max_q = '//trim(utils_int_to_str(max_q))//'.'
       else
          output_type = 'Values from averages over max_q in {'//&
               trim(utils_int_to_str(max_q+2))//','//&
               trim(utils_int_to_str(max_q+1))//'}.'
       end if

       if(dma_run == 3) goto 1000 ! only do output of averages
       cur_swex_h = swex_h + dma_run - 1

       ! { --------------------------------------------------------------------
       ! { Actual DMA calculation -- runs 1 and 2.
       ! { --------------------------------------------------------------------

       if (pub_on_root .and. pub_output_detail == VERBOSE) then
          write(stdout,'(a,i0,a,i0,a,i0,a,i0,a)') &
               CRLF//'DMA: Using SWRI ['//&
               trim(swri_library(dma_swri_h)%swri_name)//'] and SWEX ['//&
               trim(dma_swexes(dma_run)%swex_name)//'] to generate multipoles'//&
               CRLF//'DMA: (mode '//mode//', multipole orders ', min_l, '-', &
               max_l, ', ',max_q, ' Bessels) for ', n_atoms_in_dma, &
               ' atoms: '//trim(output_name)
       end if

       ! jd: ******************************************************************
       !     *** Calculate the G coefficients first (CKS's notes, eq. (21)) ***
       !     ******************************************************************

       num_sws_per_centre = dma_swexes(dma_run)%quality%num_sws_per_centre

       allocate(c_coeffs(pub_num_spins,num_sws_per_centre),stat=ierr)
       call utils_alloc_check(myself,'c_coeffs',ierr)
       allocate(sphw_coeffs(pub_num_spins,2*num_sws_per_centre),stat=ierr)
       call utils_alloc_check(myself,'sphw_coeffs',ierr)
       allocate(g_coeffs(num_sws_per_centre),stat=ierr)
       call utils_alloc_check(myself,'g_coeffs',ierr)

       total_monopole_charge = 0D0
       ! -----------------------------------------------------------------------
       ! jd: Loop over B's -- all atoms on this core                         BBB
       ! -----------------------------------------------------------------------
       loop_B:                                                                 &
       do local_b=1, par%num_atoms_on_proc(pub_my_proc_id)
          global_b = par%first_atom_on_proc(pub_my_proc_id) + local_b-1
          first_ngwf_idx_of_b = ngwf_basis(1)%first_on_atom(global_b)
          orig_global_b = par%orig_atom(global_b)

          ! jd: Ignore atoms not represented in this SWRI
          if(.not. mdl%elements(orig_global_b)%in_swri(dma_swri_h)) then
             if(debug_show_skipped) then
                write(stdout,'(a,i0,a,i0,a,a,a)') &
                     'Skipping (dma_calculate) B:orig/SFC: ',&
                     orig_global_b, '/', global_b, ' (', &
                     trim(mdl%elements(orig_global_b)%species_id),')'
             end if
             cycle
          end if

          ! jd: Start a fresh centre
          g_coeffs = 0D0

          ! --------------------------------------------------------------------
          ! jd: Loop over C's that are s-neighbours with B                   CCC
          ! --------------------------------------------------------------------
          loop_C:                                                              &
          do c_idx = swri_library(dma_swri_h)%s_atoms_nl%first_idx(global_b), &
               swri_library(dma_swri_h)%s_atoms_nl%last_idx(global_b)
             global_c = swri_library(dma_swri_h)%s_atoms_nl%neighbours(c_idx)
             first_ngwf_idx_of_c = ngwf_basis(1)%first_on_atom(global_c)

             ! -----------------------------------------------------------------
             ! jd: for all b on B                                            bbb
             ! -----------------------------------------------------------------
             loop_ngwf_b:                                                      &
             do ngwf_b = 1, ngwf_basis(1)%num_on_atom(global_b)
                global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

                ! --------------------------------------------------------------
                ! jd: for all c on C                                         ccc
                ! --------------------------------------------------------------
                loop_ngwf_c:                                                   &
                do ngwf_c = 1, ngwf_basis(1)%num_on_atom(global_c)
                   global_cc_ngwf_idx = first_ngwf_idx_of_c + ngwf_c - 1

                   coeffs_kind = merge(SW_O,SW_V,pub_dma_metric==METRIC_OVERLAP)

                   ! jd: In mode P coeffs are DKN-agnostic, in mode D not so
                   if(mode == 'D') then
                      loc_num_spins = pub_num_spins
                   else
                      loc_num_spins = 1
                   end if

                   ! -----------------------------------------------------------
                   ! jd: for all spins (if in mode D)                        sss
                   ! -----------------------------------------------------------
                   loop_spins_1:                                               &
                   do is = 1, loc_num_spins
                      ! jd: Look up expansion coefficients for this Bb Cc
                      if(present(coeffs_ht)) then
                         ! jd: Mode D -- spin-dependent
                         if(.not. SWEX_BB_CC_SYMMETRIC(cur_swex_h)) then
                            ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                            call hash_table_lookup_nocount(&
                                 sphw_coeffs(is,:), num_of_sphw_coeffs, &  ! out
                                 coeffs_ht, &                              ! in
                                 global_bb_ngwf_idx, global_cc_ngwf_idx, & ! in
                                 coeffs_kind, cur_swex_h, is)              ! in
                         else
                            ! jd: Here there's Bb-Cc symmetry and we only store one triangle
                            call hash_table_lookup_nocount(&
                                 sphw_coeffs(is,:), num_of_sphw_coeffs, &      ! out
                                 coeffs_ht, &                                  ! in
                                 min(global_bb_ngwf_idx,global_cc_ngwf_idx), & ! in
                                 max(global_bb_ngwf_idx,global_cc_ngwf_idx), & ! in
                                 coeffs_kind, cur_swex_h, is)                  ! in
                         end if
                      else
                         ! jd: Mode P -- spin dependence is only formal (is=1)
                         if(.not. SWEX_BB_CC_SYMMETRIC(cur_swex_h)) then
                            ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                            call hash_table_lookup_nocount(&
                                 sphw_coeffs(is,:), num_of_sphw_coeffs, &  ! out
                                 dma_swexes(dma_run)%coeffs_ht, &          ! in
                                 global_bb_ngwf_idx, global_cc_ngwf_idx, & ! in
                                 coeffs_kind, cur_swex_h, is)              ! in
                         else
                            ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                            call hash_table_lookup_nocount(&
                                 sphw_coeffs(is,:), num_of_sphw_coeffs, &      ! out
                                 dma_swexes(dma_run)%coeffs_ht, &              ! in
                                 min(global_bb_ngwf_idx,global_cc_ngwf_idx), & ! in
                                 max(global_bb_ngwf_idx,global_cc_ngwf_idx), & ! in
                                 coeffs_kind, cur_swex_h, is)                  ! in
                         end if
                      end if

                      if(num_of_sphw_coeffs == -1) then
                         call utils_abort(myself//': ['//&
                              trim(dma_swexes(dma_run)%swex_name)//&
                              '] Internal problem with coeff lookup [1]. Possibly &
                              &you forgot to call dma_expand_ngwf_pairs() or &
                              &dma_expand_expand_pair_densities() first. Indices: ', &
                              global_bb_ngwf_idx, global_cc_ngwf_idx, &
                              coeffs_kind, cur_swex_h)
                      end if

                      ! jd: Normalise sphw_coeffs to get c_coeffs, deal with
                      !     b-c ordering.
                      if (num_of_sphw_coeffs == num_sws_per_centre) then
                         ! jd: One-centre expansion
                         c_coeffs(is,1:num_sws_per_centre) = &
                              sphw_coeffs(is,1:num_sws_per_centre)
                      else ! jd: Two-centre, B first
                         if (global_b <= global_c) then
                            c_coeffs(is,1:num_sws_per_centre) = &
                                 2.0_DP*sphw_coeffs(is,1:num_sws_per_centre)

                         else ! jd: Two-centre, C first
                            c_coeffs(is,1:num_sws_per_centre) = &
                                 2.0_DP*sphw_coeffs(is,num_sws_per_centre+1:&
                                 num_of_sphw_coeffs)
                         end if
                      endif

                   end do loop_spins_1

                   ! -----------------------------------------------------------
                   ! jd: for all spins                                       sss
                   ! -----------------------------------------------------------
                   loop_spins_2:                                               &
                   do is = 1, pub_num_spins ! even in mode P -- we need the right dkn_el

                      if(mode == 'P') then
                         call sparse_get_element(dkn_el, denskern%m(is,PUB_1K)%p, &
                              global_cc_ngwf_idx, global_bb_ngwf_idx)
                         coeffs_is = 1
                      else
                         dkn_el = 1D0 ! jd: In mode P DKN is included in coeffs
                         coeffs_is = is
                      end if

                      ! jd: This takes care of the discrepancy between calls
                      !     to hamiltonian_lhxc_calculate, where each call is
                      !     preceded by kernel scaling; and properties, where
                      !     this is absent. Note that the default scale factor
                      !     is 0.5, not 1.0, which *may* mean that sw_expansion
                      !     does return values that are twice as large.
                      dkn_el = dkn_el * local_dkn_scale_factor

                      g_coeffs(:) = g_coeffs(:) + c_coeffs(coeffs_is,:) * dkn_el

                   end do loop_spins_2
                end do loop_ngwf_c
             end do loop_ngwf_b
          end do loop_C

          ! jd: Compute the multipoles themselves for atom B
          call compute_multipole_elms(my_multipoles(dma_run), &  ! in/out
               total_monopole_charge, &                          ! out
               dma_swexes(dma_run), swri_library(dma_swri_h), &  ! input
               g_coeffs, local_b, global_b, mdl%elements, par)   ! input

       end do loop_B
       call comms_reduce('SUM',total_monopole_charge)
       total_dma_charge = -total_core_charge + total_monopole_charge

       ! jd: Deallocate working arrays
       deallocate(sphw_coeffs,stat=ierr)
       call utils_dealloc_check(myself, 'sphw_coeffs', ierr)
       deallocate(c_coeffs,stat=ierr)
       call utils_dealloc_check(myself, 'c_coeffs', ierr)
       deallocate(g_coeffs,stat=ierr)
       call utils_dealloc_check(myself, 'g_coeffs', ierr)

       ! } --------------------------------------------------------------------
       ! } End of actual DMA calculation -- runs 1 and 2.
       ! } --------------------------------------------------------------------

       max_q = max_q - 1

1000   n_multipoles_per_atom = (1+max_l-min_l) * (1+max_l+min_l) ! == sum_min_l^max_l (2l+1)

       ! ----------------------------------------------------------------------
       ! Charge scaling (in runs 1 and 2)
       ! ----------------------------------------------------------------------

       if(dma_run < 3) then ! No scaling in dummy run
          ! jd: If charge scaling in effect, find the scaling factor
          if(pub_dma_scale_charge) then
             call utils_assert(abs(total_monopole_charge)>SAFE_DIV_EPS, &
                  trim(myself)//': Total monopole charge too close to zero: ', &
                  total_monopole_charge)
             scaling_factor = real(nelecs,kind=DP)/total_monopole_charge
          else
             scaling_factor = 1D0
          end if
       end if

       ! jd: User feedback on charge scaling
       if(pub_on_root) then
          if(dma_run < 3) then
             if(pub_dma_multipole_scaling /= 1D0) then
                write(stdout,'(a, f12.6,a)') &
                     'DMA: total charge contribution from the monopoles: ', &
                     total_monopole_charge, ' (following dma_multipole_scaling)'
             else
                write(stdout,'(a, f12.6)') &
                     'DMA: total charge contribution from the monopoles: ', &
                     total_monopole_charge
             end if

             if(pub_dma_scale_charge) then
                write(stdout,'(a, f12.6)') &
                  'DMA: expected charge:                              ', &
                  real(nelecs,kind=DP)
                write(stdout,'(a, f12.6)') &
                     'DMA: charge scaling factor:                        ', &
                     scaling_factor
             else
                write(stdout,'(a)') 'DMA: no charge scaling applied.'
             end if
          end if
       end if

       ! jd: Effect actual scaling. Do not do this in dma_run 3 or when there
       !     are no monopoles.
       if(pub_dma_scale_charge .and. dma_run < 3 .and. min_l <= 0) then
          call multipole_set_scale_monopoles(my_multipoles(dma_run), &
               scaling_factor)
          my_multipoles(dma_run)%dma_nelecs = nelecs
       end if

       ! jd: Gather multipoles to root for output to file
       call multipole_set_gather(all_multipoles(dma_run), my_multipoles(dma_run))

       ! jd: Output to file
       call dma_write_multipoles_main(all_multipoles(dma_run), mdl%elements, &
            par, dma_swri_h, file_units, output_type, &
            want_gdma_output = (dma_run == n_dma_runs))

       ! ----------------------------------------------------------------------
       ! Bessel averaging of multipoles (happens in run 2).
       ! ----------------------------------------------------------------------

       ! jd: Average the atomic multipoles into 3rd component after 2nd run
       if(dma_run == 2) then
          my_multipoles(3)%mpoles = &
               0.5_DP * (my_multipoles(1)%mpoles + my_multipoles(2)%mpoles)
          my_multipoles(3)%dma_monopole_scaling_factor = garbage_real
          ! jd: NB. When averaging scaled monopoles, there is no unique
          !         scalar scaling factor for the averaged ones, hence
          !         we tag it as garbage.
          my_multipoles(3)%dma_nelecs = nelecs
       end if

       ! jd: The third loop is dummy and only gathers at root and outputs

    end do ! over Bessel averaging runs

    ! jd: Close output files
    call dma_write_multipoles_close_outputs(file_units)

    call timer_clock(myself,2)

  end subroutine dma_calculate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dma_free(my_multipoles, dma_swexes) ! inout
    !==========================================================================!
    ! Deallocates DMA multipoles. Frees the the spherical sets and destroys    !
    ! the SWEX containers.                                                     !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2015.                                   !
    !==========================================================================!

    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET
    use sw_expansion_type, only: SW_EX

    implicit none

    ! jd: Arguments
    type(SPHERICAL_MULTIPOLE_SET), intent(inout) :: my_multipoles(3)
    type(SW_EX), intent(inout)                   :: dma_swexes(:)

    ! -------------------------------------------------------------------------

    call dma_free_multipoles(my_multipoles)
    call dma_free_swexes(dma_swexes)

  end subroutine dma_free

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dma_free_swexes(dma_swexes) ! inout
    !==========================================================================!
    ! Destroys DMA's SWEX containers.                                          !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2016.                                   !
    !==========================================================================!

    use sw_expansion_type, only: SW_EX, swx_destroy_container

    implicit none

    ! jd: Arguments
    type(SW_EX), intent(inout)                   :: dma_swexes(:)

    ! -------------------------------------------------------------------------

    ! jd: Destroy SW_EXes
    if(dma_swexes(1)%initialised) call swx_destroy_container(dma_swexes(1))
    if(dma_swexes(2)%initialised) call swx_destroy_container(dma_swexes(2))

  end subroutine dma_free_swexes

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dma_free_multipoles(my_multipoles) ! inout
    !==========================================================================!
    ! Deallocates DMA multipoles.                                              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2016.                                   !
    !==========================================================================!

    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET, &
         multipole_free_spherical_set
    use rundat, only: pub_dma_bessel_averaging

    implicit none

    ! jd: Arguments
    type(SPHERICAL_MULTIPOLE_SET), intent(inout) :: my_multipoles(3)

    ! jd: Local variables
    integer :: dma_run

    ! -------------------------------------------------------------------------

    do dma_run = 1, merge(3,1,pub_dma_bessel_averaging)
       call multipole_free_spherical_set(my_multipoles(dma_run))
    end do

  end subroutine dma_free_multipoles

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dma_output_potential(my_multipoles, & ! in
       mdl, reference_too, ngwf_basis, &           ! in
       denskern, overlap, ngwfs_on_grid)           ! only needed if reference_too
    !==========================================================================!
    ! Computes, if 'reference_too' is true, the potential of the electronic    !
    ! density through the pointwise approach and, in all cases, the potential  !
    ! from the multipole expansion on the faces of the cell (currently only    !
    ! Z=0 and Z=max) and outputs them to files.                                !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   my_multipoles (in): The array holding the multipole expansion.         !
    !   reference_too (in): If true, the reference potential is also calculated!
    !                       and output (this is costly).                       !
    !   ngwf_basis (in)                                                        !
    !   denskern, overlap, ngwfs_on_grid (in, opt):                            !
    !   These only need to be supplied if reference_too is true.               !
    !--------------------------------------------------------------------------!
    ! Rewritten from scratch by Jacek Dziedzic in May 2014.                    !
    ! Made reference potential output optional by Jacek Dziedzic in Feb 2015.  !
    ! Modified to use SPHERICAL_MULTIPOLE_SET by Jacek Dziedzic in May 2015.   !
    ! Generalised to only need most params for the reference calculation, by   !
    ! Jacek Dziedzic in July 2015.                                             !
    !==========================================================================!

    use datatypes, only: FUNCTIONS
    use cell_grid, only: GRID_INFO, cell_grid_real_pt
    use comms, only: comms_reduce, pub_my_proc_id
    use constants, only: DP
    use density, only: density_on_grid
    use finite_differences, only: finite_difference_set_geometry
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, OPERATOR(+), OPERATOR(*), operator(-), magnitude
    use model_type, only: MODEL
    use multigrid_methods, only: multigrid_prepare_bound_cond, mg
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET, &
         multipole_pot_of_spherical_set
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_num_spins, pub_spin_fac, pub_dma_bessel_averaging, &
         pub_dma_max_q, pub_num_kpoints, PUB_1K
    use sparse, only: SPAM3, sparse_get_par
    use sparse_array, only: SPAM3_ARRAY
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_abort, utils_assert

    implicit none

    ! jd: arguments
    type(SPHERICAL_MULTIPOLE_SET), intent(in)  :: my_multipoles(3)
    type(MODEL), intent(in)                    :: mdl
    logical, intent(in)                        :: reference_too
    type(FUNC_BASIS), intent(in)               :: ngwf_basis
    type(SPAM3_ARRAY), intent(in), optional    :: denskern
    type(SPAM3), intent(in), optional          :: overlap
    type(FUNCTIONS), intent(in), optional      :: ngwfs_on_grid

    ! jd: Local variables
    integer :: n_dma_runs, dma_run
    integer :: ierr
    real(kind=DP), allocatable :: density_fine(:,:,:,:)
    ! jd: potential due to electronic density calculated through
    !     multigrid_prepare_bound_cond(). Zeroes inside, values on the faces.
    real(kind=DP), allocatable :: ref_potential_on_grid(:,:,:)
    ! jd: potential due to electronic density from the multipoles, on the
    !     top face of the simulation cell only.
    real(kind=DP), allocatable :: mp_potential_on_face_top(:,:)
    ! jd: potential due to electronic density from the multipoles, on the
    !     bottom face of the simulation cell only.
    real(kind=DP), allocatable :: mp_potential_on_face_bot(:,:)
    integer         :: file_unit
    integer         :: i1, i2, islab12 ! jd: Gridpoint indices
    type(POINT)     :: r_point, R_I
    real(kind=DP)   :: r_pt(3)
    integer         :: hairy_shaved
    integer         :: is
    integer         :: ipt
    integer         :: I
    real(kind=DP)   :: d
    logical         :: shave_me
    type(GRID_INFO) :: grid
    integer     :: owner_of_first_slab, owner_of_last_slab ! jd: Proc numbers
    character(len=*), parameter :: ref_pot_filename_template(2) = &
         (/'ref_elec_potential_for_z=       ', &
         'ref_elec_shaved_potential_for_z='/)
    character(len=*), parameter :: mp_pot_filename_template = &
         'mp_elec_potential_for_z='
    character(len=256)          :: pot_filename
    character(len=80)           :: str_coord
    character(len=80)           :: str_coord2
    character(len=80)           :: str_max_q
    character(len=*), parameter :: myself = 'dma_output_potential'
    type(PARAL_INFO), pointer   :: par

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine dma_output_potential not ready yet for more&
         & than one k-point.')

    if(reference_too) then
       call utils_assert(present(denskern) .and. present(overlap) .and. &
            present(ngwfs_on_grid),myself//&
            ': All optional arguments must be specified when "reference_too"&
            & is .true.)')
    end if

    ! rc2013: get the parallel strategy from the sparse matrix
    call sparse_get_par(par, overlap)

    grid = mdl%fine_grid

    owner_of_first_slab = grid%proc_slab12(1)
    owner_of_last_slab = grid%proc_slab12(grid%n3)

    ! jd: Ensure we have a cuboid cell
    call utils_assert(&
         grid%da1%y == 0D0 .and. grid%da1%z == 0D0 .and. &
         grid%da2%z == 0D0 .and. grid%da2%x == 0D0 .and. &
         grid%da3%x == 0D0 .and. grid%da3%y == 0D0, &
         trim(myself)//' only supports cuboid cells, sorry.')

    ! jd: -----------------------------------------------------------------
    ! jd: --- Obtain the reference potential directly from the density  ---
    ! jd: -----------------------------------------------------------------

    if(reference_too) then

       ! jd: Prepare electronic density on fine grid
       allocate(density_fine(grid%ld1, grid%ld2, &
            grid%max_slabs12, pub_num_spins),stat=ierr)
       call utils_alloc_check(myself,'density_fine',ierr)
       call density_on_grid(density_fine, grid, mdl%dbl_grid, &
            mdl%cell, mdl%fftbox, denskern%m(:,PUB_1K), overlap, &
            ngwfs_on_grid, ngwf_basis, ngwfs_on_grid, ngwf_basis)

       ! jd: Account for spin degeneracy
       density_fine = density_fine * pub_spin_fac

       ! jd: Allocate array for potential, set up MG
       allocate(ref_potential_on_grid(grid%ld1, grid%ld2, &
            grid%max_slabs12),stat=ierr)
       call utils_alloc_check(myself,'ref_potential_on_grid',ierr)
       mg = finite_difference_set_geometry(grid,1)

       ! jd: Calculate the potential due to this density: hairy and shaved.
       do hairy_shaved = 1, 2

          ! jd: Calculate potential
          call multigrid_prepare_bound_cond(ref_potential_on_grid, &
               density_fine, grid, mdl%cell, &
               bc_coarseness=1, bc_surface_coarseness=1)

          ! jd: Output the Z=0 face of the reference potential on the grid
          if(pub_my_proc_id == owner_of_first_slab) then
             file_unit = utils_unit()
             write(str_coord,'(f9.4)') 0D0
             pot_filename = trim(ref_pot_filename_template(hairy_shaved))//&
                  trim(adjustl(str_coord))//trim('.txt')
             open(unit=file_unit, file=pot_filename, status="unknown", err=100)
             do i1= 1, grid%n1
                do i2 = 1, grid%n2
                   r_point = &
                        real((i1-1),kind=DP) * grid%da1 + &
                        real((i2-1),kind=DP) * grid%da2
                   write(file_unit,'(f12.4,f12.4,e18.9)') &
                        r_point%X, r_point%Y, ref_potential_on_grid(i1,i2,1)
                end do
             end do
             close(file_unit,err=200)
          end if

          ! jd: Output the Z=max face of the reference potential on the grid
          if(pub_my_proc_id == owner_of_last_slab) then
             file_unit = utils_unit()
             write(str_coord2,'(f9.4)') &
                  real((grid%n3-1),kind=DP) * grid%da3%z
             pot_filename = trim(ref_pot_filename_template(hairy_shaved))//&
                  trim(adjustl(str_coord2))//trim('.txt')
             open(unit=file_unit,file=pot_filename, status="unknown", err=100)
             do i1= 1, grid%n1
                do i2 = 1, grid%n2
                   r_point = &
                        real((i1-1),kind=DP) * grid%da1 + &
                        real((i2-1),kind=DP) * grid%da2
                   write(file_unit,'(f12.4,f12.4,e18.9)') r_point%X, r_point%Y, &
                        ref_potential_on_grid(i1,i2,grid%num_my_slabs12)
                end do
             end do
             close(file_unit,err=200)
          end if

          ! jd: Shave density between iterations 1 and 2
          if(hairy_shaved == 1) then
             do is = 1, pub_num_spins
                do ipt=1, grid%num_my_slabs12 * grid%n2 * &
                     grid%n1
                   i1 = modulo(ipt-1,grid%n1) + 1
                   i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
                   islab12 = (ipt-(i2-1)*grid%n1-i1) / &
                        (grid%n1*grid%n2) + 1

                   call cell_grid_real_pt(r_pt,i1,i2,islab12,grid)

                   shave_me = .true.
                   do I = 1, par%nat
                      R_I = mdl%elements(I)%centre

                      d = magnitude(POINT(r_pt(1),r_pt(2),r_pt(3))-R_I) ! @MIC@

                      if(d <= mdl%elements(I)%radius) then
                         shave_me = .false.
                         exit
                      end if
                   end do

                   if(shave_me) density_fine(i1,i2,islab12,is) = 0D0

                end do ! points to shave
             end do ! spins
          end if ! if shaving

       end do ! hairy_shaved

       ! jd: Density no longer needed
       deallocate(density_fine,stat=ierr)
       call utils_dealloc_check(myself,'density_fine',ierr)

       ! jd: Reference potential no longer needed
       deallocate(ref_potential_on_grid,stat=ierr)
       call utils_dealloc_check(myself,'ref_potential_on_grid',ierr)

    end if ! reference_too

    ! jd: ---------------------------------------------------------------------
    ! jd: --- Obtain potential from the multipoles on top and bottom faces  ---
    ! jd: ---------------------------------------------------------------------

    allocate(mp_potential_on_face_top(grid%n1, grid%n2), &
         stat=ierr)
    call utils_alloc_check(myself,'mp_potential_on_face_top',ierr)
    allocate(mp_potential_on_face_bot(grid%n1, grid%n2), &
         stat=ierr)
    call utils_alloc_check(myself,'mp_potential_on_face_bot',ierr)

    if(pub_dma_bessel_averaging) then
       n_dma_runs = 3
    else
       n_dma_runs = 1
    end if

    do dma_run = 1, n_dma_runs

       ! jd: This calculates the potential from the multipoles of atoms local to
       !     this proc only, need to reduce later.
       do i2 = 1, grid%n2
          do i1 = 1, grid%n1
             ! jd: Points on Z=0 face
             r_point = &
                  real((i1-1),kind=DP) * grid%da1 + &
                  real((i2-1),kind=DP) * grid%da2
                  mp_potential_on_face_top(i1,i2) = &
                       multipole_pot_of_spherical_set(my_multipoles(dma_run), &
                       r_point)
             ! jd: Points on Z=max face
             r_point = &
                  real((i1-1),kind=DP) * grid%da1 + &
                  real((i2-1),kind=DP) * grid%da2 + &
                  real((grid%n3-1),kind=DP) * grid%da3

                  mp_potential_on_face_bot(i1,i2) = &
                       multipole_pot_of_spherical_set(my_multipoles(dma_run), &
                       r_point)
          end do
       end do

       ! jd: Reduce over multipoles from all procs
       call comms_reduce('SUM',mp_potential_on_face_top)
       call comms_reduce('SUM',mp_potential_on_face_bot)

       if(dma_run == 1) write(str_max_q,'(i0)') pub_dma_max_q
       if(dma_run == 2) write(str_max_q,'(i0)') pub_dma_max_q-1
       if(dma_run == 3) write(str_max_q,'(a3)') 'avg'

       ! jd: Output the Z=0 face of the multipole potential on the grid
       if(pub_my_proc_id == owner_of_first_slab) then
          file_unit = utils_unit()
          write(str_coord,'(f9.4)') 0D0
          pot_filename = trim(mp_pot_filename_template)//&
               trim(adjustl(str_coord))//'_'//trim(adjustl(str_max_q))//'.txt'
          open(unit=file_unit,file=pot_filename, status="unknown", err=100)
          do i1= 1, grid%n1
             do i2 = 1, grid%n2
                r_point = &
                     real((i1-1),kind=DP) * grid%da1 + &
                     real((i2-1),kind=DP) * grid%da2
                write(file_unit,'(f12.4,f12.4,e18.9)') &
                     r_point%X, r_point%Y, mp_potential_on_face_top(i1,i2)
             end do
          end do
          close(file_unit,err=200)
       end if

       ! jd: Output the Z=max face of the multipole potential on the grid
       if(pub_my_proc_id == owner_of_last_slab) then
          file_unit = utils_unit()
          write(str_coord2,'(f9.4)') &
               real((grid%n3-1),kind=DP) * grid%da3%z
          pot_filename = trim(mp_pot_filename_template)//&
               trim(adjustl(str_coord2))//'_'//trim(adjustl(str_max_q))//'.txt'
          open(unit=file_unit,file=pot_filename, status="unknown", err=100)
          do i1= 1, grid%n1
             do i2 = 1, grid%n2
                r_point = &
                     real((i1-1),kind=DP) * grid%da1 + &
                     real((i2-1),kind=DP) * grid%da2
                write(file_unit,'(f12.4,f12.4,e18.9)') &
                     r_point%X, r_point%Y, mp_potential_on_face_bot(i1,i2)
             end do
          end do
          close(file_unit,err=200)
       end if
    end do ! over dma_runs

    deallocate(mp_potential_on_face_bot,stat=ierr)
    call utils_dealloc_check(myself,'mp_potential_on_face_bot',ierr)
    deallocate(mp_potential_on_face_top,stat=ierr)
    call utils_dealloc_check(myself,'mp_potential_on_face_top',ierr)

    call timer_clock(myself,2)

    return

    ! jd: I/O error handling
100 call utils_abort('Error during creation of file: '//pot_filename)
200 call utils_abort('Error during closing of file: '//pot_filename)

  end subroutine dma_output_potential

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dma_write_multipoles(output_name, multipoles, elements, par, &
       dma_swri_h, gather_first)
    !==========================================================================!
    ! Writes DMA multipoles to user-readable text files and, optionally,       !
    ! a GDMA-like file.                                                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   output_name (in): Infix to be used in generating output file names.    !
    !   multipoles (in): DMA-generated multipole sets.                         !
    !   elements (in): The usual. Needed to potentially skip atoms outside of  !
    !                  DMA swri.                                               !
    !   dma_swri_h (in): Handle to DMA swri, needed to potentially skip atoms  !
    !                    outside of DMA swri.                                  !
    !   gather_first (in): If .true., the multipole set will be internally     !
    !                      gathered across procs to a temporary. If .false.,   !
    !                      caller is responsible for gathering the multipoles  !
    !                      to root first.                                      !
    !   par (in):  Parallel strategy for region associated with elements.      !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2017.                                  !
    !==========================================================================!

    use ion, only: ELEMENT
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET, multipole_set_gather, &
         multipole_free_spherical_set
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_dma_bessel_averaging, pub_dma_max_q
    use utils, only: utils_int_to_str

    implicit none

    ! jd: Arguments
    character(len=*), intent(in)                      :: output_name
    type(SPHERICAL_MULTIPOLE_SET), intent(in), target :: multipoles(3)
    type(PARAL_INFO), intent(in)                      :: par
    type(ELEMENT), intent(in)                         :: elements(:)
    integer, intent(in)                               :: dma_swri_h
    logical, intent(in)                               :: gather_first

    ! jd: Local variables
    type(SPHERICAL_MULTIPOLE_SET), target  :: all_multipoles(3)
    type(SPHERICAL_MULTIPOLE_SET), pointer :: cur_multipoles
    integer                                :: file_units(N_DMA_FILES)
    integer                                :: n_dma_runs
    integer                                :: dma_run
    integer                                :: max_q
    character(len=128)                     :: output_type

    ! -------------------------------------------------------------------------

    if(pub_dma_bessel_averaging) then
       n_dma_runs = 3
    else
       n_dma_runs = 1
    end if

    max_q = pub_dma_max_q

    call dma_write_multipoles_open_outputs(output_name, file_units)

    do dma_run = 1, n_dma_runs

       if(gather_first) then
          call multipole_set_gather(all_multipoles(dma_run), multipoles(dma_run))
          cur_multipoles => all_multipoles(dma_run)
       else
          cur_multipoles => multipoles(dma_run)
       end if

       if(dma_run < 3) then
          output_type = 'Values for max_q = '//trim(utils_int_to_str(max_q))//'.'
       else
          output_type = 'Values from averages over max_q in {'//&
               trim(utils_int_to_str(max_q+2))//','//&
               trim(utils_int_to_str(max_q+1))//'}.'
       end if

       call dma_write_multipoles_main(cur_multipoles, &
            elements, par, dma_swri_h, file_units, output_type, &
            want_gdma_output = (dma_run == n_dma_runs))

       max_q = max_q - 1

       if(gather_first) then
          call multipole_free_spherical_set(all_multipoles(dma_run))
       end if

    end do

    call dma_write_multipoles_close_outputs(file_units)

  end subroutine dma_write_multipoles

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! ---------------------------------------------------------------------------
  ! ---- P R I V A T E   S U B R O U T I N E S   A N D   F U N C T I O N S ----
  ! ---------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dma_write_multipoles_open_outputs(output_name, file_units)
    !==========================================================================!
    ! Opens DMA output files in preparation for writing out the multipoles.    !
    !--------------------------------------------------------------------------!
    ! Extracted by Jacek Dziedzic in July 2017 from dma_calculate().           !
    !==========================================================================!

    use comms, only: pub_on_root
    use rundat, only: pub_rootname
    use utils, only: utils_unit, utils_abort

    implicit none

    ! jd: Arguments
    character(len=*), intent(in) :: output_name
    integer, intent(out)         :: file_units(N_DMA_FILES)

    ! jd: Local variables
    integer :: ifile
    character(len=*), parameter :: file_name_suffixes(N_DMA_FILES) = (/ &
         'dma_multipoles_atomcentred_elec_spherical_real.txt     ', &
         'dma_multipoles_atomcentred_elec_cartesian_traceless.txt', &
         'dma_multipoles_atomcentred_elec_cartesian_primitive.txt', &
         'dma_multipoles_refpt_elec_cartesian_primitive.txt      ', &
         'dma_multipoles_refpt_total_cartesian_primitive.txt     ', &
         'dma_multipoles_refpt_elec_cartesian_traceless.txt      ', &
         'dma_multipoles_refpt_total_cartesian_traceless.txt     ', &
         'dma_multipoles_gdma_like.txt                           ' /)
    character(len=512)          :: file_names(N_DMA_FILES)

    ! -------------------------------------------------------------------------

    if(pub_on_root) then
       do ifile = 1, N_DMA_FILES
          write(file_names(ifile),'(a)') trim(pub_rootname)//'_'//&
               trim(output_name)//'_'//trim(file_name_suffixes(ifile))
          file_units(ifile) = utils_unit()
          open(unit=file_units(ifile), file=trim(file_names(ifile)), &
               status="unknown", err=1000)
       end do
    end if

    return

    ! jd: I/O error handling
1000 call utils_abort('Error during creation of file: '//file_names(ifile))

  end subroutine dma_write_multipoles_open_outputs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dma_write_multipoles_close_outputs(file_units)
    !==========================================================================!
    ! Closes DMA output files.                                                 !
    !--------------------------------------------------------------------------!
    ! Extracted by Jacek Dziedzic in July 2017 from dma_calculate().           !
    !==========================================================================!

    use comms, only: pub_on_root
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    integer, intent(in) :: file_units(N_DMA_FILES)

    ! jd: Local variables
    integer :: ifile

    ! -------------------------------------------------------------------------

    if(pub_on_root) then
       do ifile = 1, N_DMA_FILES
          close(file_units(ifile), err=1000)
       end do
    end if

    return

    ! jd: I/O error handling
1000 call utils_abort('Error closing of DMA output file #',ifile)

  end subroutine dma_write_multipoles_close_outputs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dma_write_multipoles_main(all_multipoles, elements, par, dma_swri_h, &
       file_units, output_type, want_gdma_output)
    !==========================================================================!
    ! Outputs DMA multipoles to files.                                         !
    ! Caller is responsible for calling dma_write_multipoles_open_outputs()    !
    ! first, and dma_write_multipoles_close_outputs() subsequently.            !
    !--------------------------------------------------------------------------!
    ! Caveat:                                                                  !
    !   'all_multipoles' is expected to have been gathered to root in advance. !
    !--------------------------------------------------------------------------!
    ! Extracted by Jacek Dziedzic in July 2017 from dma_calculate().           !
    !==========================================================================!

    use comms, only: pub_on_root
    use ion, only: ELEMENT
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_devel_code
    use utils, only: utils_devel_code, utils_assert

    implicit none

    ! jd: Arguments
    type(SPHERICAL_MULTIPOLE_SET), intent(in) :: all_multipoles
    type(ELEMENT), intent(in)                 :: elements(:)
    type(PARAL_INFO), intent(in)              :: par
    integer, intent(in)                       :: dma_swri_h
    integer, intent(in)                       :: file_units(N_DMA_FILES)
    character(len=*), intent(in)              :: output_type
    logical, intent(in)                       :: want_gdma_output

    ! jd: Local variables
    logical :: debug_show_skipped
    integer :: global_b_orig_idx, global_b
    integer :: mp_begin, mp_end
    character(len=*), parameter :: myself = 'dma_write_multipoles_main'

    ! -------------------------------------------------------------------------

    debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
         pub_devel_code, no_bcast = .true., no_warn=.true.)

    if(pub_on_root) then
       ! jd: Output in input-ordering, but all_multipoles(:) uses indices
       !     after SF-reordering, so take care to undo that.
       do global_b_orig_idx = 1, par%nat
          global_b = par%distr_atom(global_b_orig_idx)
          mp_begin = (global_b-1) * all_multipoles%n_multipoles_per_site+1
          mp_end = global_b * all_multipoles%n_multipoles_per_site
          call utils_assert(mp_begin >= lbound(all_multipoles%mpoles,1) .and. &
               mp_begin <= ubound(all_multipoles%mpoles,1) .and. &
               mp_end >= lbound(all_multipoles%mpoles,1) .and. &
               mp_end <= ubound(all_multipoles%mpoles,1), myself//': Multipole &
               &index out of range. Forgot to gather multipole set?', &
               mp_begin, mp_end, lbound(all_multipoles%mpoles,1), &
               ubound(all_multipoles%mpoles,1))
          call print_mp_elms_for_atom(all_multipoles%mpoles(mp_begin:mp_end), &
               all_multipoles%min_l, all_multipoles%max_l, &
               dma_swri_h, debug_show_skipped, global_b_orig_idx, &
               elements, file_units, output_type, want_gdma_output)
       end do
    end if

  end subroutine dma_write_multipoles_main

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dma_calculate_dcoeffs_ngwf_pairs(dma_swex, ngwf_basis, mdl, &
       dma_swri_h, coeffs_kind, swex_h)
    !==========================================================================!
    ! Calculates the 'D-coefficients' needed for LNV gradient in the presence  !
    ! of polarisable embedding with the QM* representation, when fitting       !
    ! NGWF pairs. See                                                          !
    ! [1] J. Dziedzic, Y. Mao, Y. Shao, J. Ponder, T. Head-Gordon,             !
    !     M. Head-Gordon, C.-K. Skylaris, J. Chem. Phys. 145, 124106 (2016).   !
    !                                                                          !
    ! The D-coefficients are communicated across procs, destinations are       !
    ! determined from SWRI's S-atoms neighbour list. The D-coefficients are    !
    ! stored in dma_swex%coeffs_ht.                                            !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   dma_swex (inout): The SW_EX container for which we calculate D-coeffs. !
    !                     dma_expand_ngwf_pairs() or dma_expand_pair_densities !
    !                     must have been called first on it.                   !
    !   ngwf_basis (in): The NGWF basis on which we operate, used for atom     !
    !                    counting.                                             !
    !   mdl (in): The model, needed for elements.                              !
    !   dma_swri_h (in): Handle to the DMA SWRI container.                     !
    !   coeffs_kind (in): The type of coefficients (SW_O, SW_V, ...).          !
    !   swex_h (in): Used as a key for locating coefficients.                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October 2016.                               !
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: SQRT_PI, PI, STDOUT, SW_D, SWEX_BB_CC_SYMMETRIC
    use function_basis, only: FUNC_BASIS
    use hash_table, only: hash_table_lookup_nocount, hash_table_add
    use model_type, only: MODEL
    use rundat, only: pub_devel_code
    use sw_expansion_type, only: SW_EX
    use sw_expansion, only: swx_communicate_coeffs
    use sw_resolution_of_identity, only: swri_sph_bess_solidharm_int, &
         swri_sw_to_ablmq, swri_library
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check, &
         utils_devel_code, utils_abort

    implicit none

    ! jd: --- Arguments ---
    type(SW_EX), intent(inout)                 :: dma_swex
    type(FUNC_BASIS), intent(in)               :: ngwf_basis
    type(MODEL), intent(in)                    :: mdl
    integer, intent(in)                        :: dma_swri_h
    integer, intent(in)                        :: coeffs_kind
    integer, intent(in)                        :: swex_h

    ! jd: --- Local variables ---
    ! jd: Counters for atoms and NGWFs
    integer :: local_b, global_b, orig_global_b
    integer :: c_idx, global_c, orig_global_c
    integer :: first_ngwf_idx_of_b, first_ngwf_idx_of_c
    integer :: global_bb_ngwf_idx, global_cc_ngwf_idx
    integer :: ngwf_b, ngwf_c

    ! jd: SW-counting and indices
    integer :: num_sws_per_centre
    integer :: swex_sw, swri_sw
    integer :: i_l, i_m, i_b, species
    integer :: offset_of_l, offset_of_lm
    integer :: ndata
    real(kind=DP) :: a, q

    ! jd: Bessel integral J_lq and normalisation factors
    real(kind=DP) :: J_lq
    real(kind=DP) :: factor

    ! jd: Subsequent levels of coefficients
    real(kind=DP), allocatable  :: sphw_coeffs(:)
    real(kind=DP), allocatable  :: c_coeffs(:)
    real(kind=DP), allocatable  :: d_beta_gamma(:)

    ! jd: Varia
    logical       :: debug_show_skipped
    logical       :: overwrite
    integer       :: num_of_sphw_coeffs
    integer       :: ierr
    character(len=*), parameter :: myself = 'dma_calculate_dcoeffs_ngwf_pairs'

    ! -------------------------------------------------------------------------

    num_sws_per_centre = dma_swex%quality%num_sws_per_centre

    debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
         pub_devel_code, no_bcast = .true., no_warn=.true.)

    allocate(c_coeffs(num_sws_per_centre),stat=ierr)
    call utils_alloc_check(myself,'c_coeffs',ierr)

    allocate(sphw_coeffs(2*num_sws_per_centre),stat=ierr)
    call utils_alloc_check(myself,'sphw_coeffs',ierr)

    allocate(d_beta_gamma((1+dma_swex%quality%max_l-dma_swex%quality%min_l)* &
         (1+dma_swex%quality%max_l+dma_swex%quality%min_l)),stat=ierr)
    call utils_alloc_check(myself,'d_beta_gamma',ierr)

    ! --------------------------------------------------------------------------
    ! jd: Loop over B's -- all atoms on this core                            BBB
    ! --------------------------------------------------------------------------
    ! rc2013: PAR WARNING!
    loop_B:                                                                    &
    do local_b=1, mdl%par%num_atoms_on_proc(pub_my_proc_id)
       global_b = mdl%par%first_atom_on_proc(pub_my_proc_id) + local_b-1
       first_ngwf_idx_of_b = ngwf_basis%first_on_atom(global_b)
       orig_global_b = mdl%par%orig_atom(global_b)

       ! jd: No ignoring of B atoms that are outside SWRI happens here, because
       !     we need D-coeffs for all pairs of atoms, except outsider-outsider
       !     pairs.

       ! -----------------------------------------------------------------------
       ! jd: Loop over C's that are s-neighbours with B                      CCC
       ! -----------------------------------------------------------------------
       loop_C:                                                                 &
       do c_idx = swri_library(dma_swri_h)%s_atoms_nl%first_idx(global_b), &
            swri_library(dma_swri_h)%s_atoms_nl%last_idx(global_b)
          global_c = swri_library(dma_swri_h)%s_atoms_nl%neighbours(c_idx)
          orig_global_c = mdl%par%orig_atom(global_c)
          first_ngwf_idx_of_c = ngwf_basis%first_on_atom(global_c)

          ! jd: Only ignore atom pairs where both atoms are outside this SWRI
          if(.not. mdl%elements(orig_global_b)%in_swri(dma_swri_h) .and. &
               .not. mdl%elements(orig_global_c)%in_swri(dma_swri_h)) then
             if(debug_show_skipped) then
                write(stdout,'(a,i0,a,i0,a,a,a,i0,a,i0,a,i0,a)') &
                     'Skipping pair ('//trim(myself)//') B:orig/SFC: ',&
                     orig_global_b, '/', global_b, ' (', &
                     trim(mdl%elements(orig_global_b)%species_id),&
                     ') C:orig/SFC: ', orig_global_c, '/', global_c, ' (', &
                     trim(mdl%elements(orig_global_c)%species_id),')'
             end if
             cycle
          end if

          ! --------------------------------------------------------------------
          ! jd: for all b on B                                               bbb
          ! --------------------------------------------------------------------
          loop_ngwf_b:                                                         &
          do ngwf_b = 1, ngwf_basis%num_on_atom(global_b)
             global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

             ! -----------------------------------------------------------------
             ! jd: for all c on C                                            ccc
             ! -----------------------------------------------------------------
             loop_ngwf_c:                                                      &
             do ngwf_c = 1, ngwf_basis%num_on_atom(global_c)
                global_cc_ngwf_idx = first_ngwf_idx_of_c + ngwf_c - 1

                ! jd: Look up expansion coefficients for this Bb Cc
                if(.not. SWEX_BB_CC_SYMMETRIC(swex_h)) then
                   ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                   call hash_table_lookup_nocount(&
                        sphw_coeffs, num_of_sphw_coeffs, &                 ! out
                        dma_swex%coeffs_ht, &                              ! in
                        global_bb_ngwf_idx, global_cc_ngwf_idx, &          ! in
                        coeffs_kind, swex_h, 1)                            ! in
                else
                   ! jd: Here there's Bb-Cc symmetry and we only store one triangle
                   call hash_table_lookup_nocount(&
                        sphw_coeffs, num_of_sphw_coeffs, &                 ! out
                        dma_swex%coeffs_ht, &                              ! in
                        min(global_bb_ngwf_idx,global_cc_ngwf_idx), &      ! in
                        max(global_bb_ngwf_idx,global_cc_ngwf_idx), &      ! in
                        coeffs_kind, swex_h, 1)                            ! in
                end if

                if(num_of_sphw_coeffs == -1) then
                   call utils_abort(myself//': ['//&
                        trim(dma_swex%swex_name)//&
                     '] Internal problem with coeff lookup [2]. &
                     &Possibly forgot to call dma_expand_ngwf_pairs() or &
                     &dma_expand_pair_densities() first.', &
                     global_bb_ngwf_idx,global_cc_ngwf_idx, num_of_sphw_coeffs)
                end if

                ! jd: Normalise sphw_coeffs to get c_coeffs, deal with
                !     b-c ordering.
                if (num_of_sphw_coeffs == num_sws_per_centre) then
                   ! jd: One-centre expansion
                   c_coeffs(1:num_sws_per_centre) = &
                        sphw_coeffs(1:num_sws_per_centre)
                else ! jd: Two-centre, B first
                   if (global_b <= global_c) then
                      c_coeffs(1:num_sws_per_centre) = &
                           2.0_DP*sphw_coeffs(1:num_sws_per_centre)

                   else ! jd: Two-centre, C first
                      c_coeffs(1:num_sws_per_centre) = &
                           2.0_DP*sphw_coeffs(num_sws_per_centre+1:&
                           num_of_sphw_coeffs)
                   end if
                endif

                d_beta_gamma(:) = 0D0
                do swex_sw = 1, num_sws_per_centre
                   species = 1 ! @IDENTICAL RADII
                   swri_sw = dma_swex%swex_sw_to_swri_sw(swex_sw)
                   call swri_sw_to_ablmq(swri_library(dma_swri_h), &
                        swri_sw, a, i_b ,i_l ,i_m, q)
                   offset_of_l = (i_l+1)**2-i_l-dma_swex%quality%min_l**2
                   offset_of_lm = offset_of_l + i_m

                   call utils_assert(i_l <= dma_swex%quality%max_l, &
                        myself//': Internal error (l > max_l)',i_l, &
                        dma_swex%quality%max_l)
                   call utils_assert(i_l >= dma_swex%quality%min_l, &
                        myself//': Internal error (l < min_l)',i_l, &
                        dma_swex%quality%min_l)
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

                   d_beta_gamma(offset_of_lm) = &
                        d_beta_gamma(offset_of_lm) + &
                        c_coeffs(swex_sw) * factor * J_lq
                end do ! over SWs

                ndata = (1 + dma_swex%quality%max_l - dma_swex%quality%min_l) *&
                     (1 + dma_swex%quality%max_l + dma_swex%quality%min_l)
                ! jd: Store d_beta_gamma in hash table
                overwrite = .true. ! (potentially overwrite old data)

                call hash_table_add(dma_swex%dcoeffs_ht, d_beta_gamma, &
                     ndata, global_bb_ngwf_idx, global_cc_ngwf_idx, &
                     coeffs_kind + SW_D, swex_h, &
                     overwrite = overwrite)

             end do loop_ngwf_c
          end do loop_ngwf_b
       end do loop_C
    end do loop_B

    ! jd: Deallocate working arrays
    deallocate(d_beta_gamma,stat=ierr)
    call utils_dealloc_check(myself, 'd_beta_gamma', ierr)
    deallocate(sphw_coeffs,stat=ierr)
    call utils_dealloc_check(myself, 'sphw_coeffs', ierr)
    deallocate(c_coeffs,stat=ierr)
    call utils_dealloc_check(myself, 'c_coeffs', ierr)

    ! jd: Communicate dcoeffs across procs
    ! rc2013: EMBED_FIX!
    call swx_communicate_coeffs(dma_swex%dcoeffs_ht, &
         ngwf_basis, swri_library(dma_swri_h)%s_atoms_nl, &
         swri_library(dma_swri_h)%s_atoms_nl, swri_library(dma_swri_h), &
         mdl%par, mdl%elements, coeffs_kind + SW_D, swex_h, 'P')
         ! NB: Aliasing of arguments 3-5 is OK as all are intent in.

  end subroutine dma_calculate_dcoeffs_ngwf_pairs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_multipoles(my_multipoles, & ! out
       elements, min_l, max_l, par)           ! in
    !==========================================================================!
    ! Sets up a datastructure for keeping DMA multipoles for atoms local to    !
    ! this proc. Spherical real representation is used.                        !
    ! The internal array storing multipole values is allocated and the centres !
    ! are set up (values are copied from elements). Values of the multipoles   !
    ! are set to zero.                                                         !
    !--------------------------------------------------------------------------!
    ! Current version of this subroutine is due to Jacek Dziedzic.             !
    ! Cleaned up by Jacek Dziedzic in May 2014.                                !
    ! Cleaned up by Jacek Dziedzic in February 2015.                           !
    ! Largely rewritten by Jacek Dziedzic in May 2015 to use multipole_ops.    !
    ! Modified to remove pub_par by Joseph Prentice, September 2018            !
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use ion, only: ELEMENT
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET, &
         multipole_init_spherical_set
    use parallel_strategy, only: PARAL_INFO

    implicit none

    ! jd: Arguments
    type(SPHERICAL_MULTIPOLE_SET), intent(out) :: my_multipoles(3)
    type(PARAL_INFO), intent(in)               :: par
    type(ELEMENT), intent(in)                  :: elements(par%nat)
    integer, intent(in)                        :: min_l
    integer, intent(in)                        :: max_l

    ! jd: Local variables
    integer :: n_atoms
    integer :: i_atom
    integer :: global_i_atom
    integer :: dma_run

    ! -------------------------------------------------------------------------

    do dma_run = 1, 3
       n_atoms = par%num_atoms_on_proc(pub_my_proc_id)
       call multipole_init_spherical_set(my_multipoles(dma_run), min_l, max_l, &
            n_atoms, 'DMA multipoles')

       my_multipoles(dma_run)%mpoles(:) = 0D0

       do i_atom = 1, n_atoms
          global_i_atom = par%first_atom_on_proc(pub_my_proc_id) + i_atom -1
          my_multipoles(dma_run)%centres(i_atom) = &
               elements(par%orig_atom(global_i_atom))%centre
          my_multipoles(dma_run)%safe_radii(i_atom) = &
               elements(par%orig_atom(global_i_atom))%radius
       end do
    end do

  end subroutine init_multipoles

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_multipole_elms(my_multipoles, &            ! in/out
       total_monopole, &                                        ! out, added to
       swex, swri, g_coeffs, local_atom_idx, global_atom_idx, & ! in
       elements, par)                                           ! in
    !==========================================================================!
    ! Computes elements of a multipole vector.                                 !
    ! The multipoles are atom-centered, electronic-only, and in spherical rep. !
    ! total_monopole is rank-local and is added to.                            !
    ! my_multipoles stores only multipoles local to this rank.                 !
    !--------------------------------------------------------------------------!
    ! Current version of this subroutine is due to Jacek Dziedzic.             !
    ! Cleaned up by Jacek Dziedzic in May 2014.                                !
    ! Largely rewritten by Jacek Dziedzic in February 2015 for generalising    !
    ! towards SW_EX/SW_RI.                                                     !
    ! Modified to remove pub_par by Joseph Prentice, September 2018            !
    !==========================================================================!

    use constants, only: PI, SQRT_PI
    use ion, only: ELEMENT
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_dma_multipole_scaling, &
         pub_dma_dipole_scaling, pub_dma_quadrupole_scaling
    use sparse, only: sparse_get_par
    use sw_expansion_type, only: SW_EX
    use sw_resolution_of_identity, only: SW_RI, swri_sph_bess_solidharm_int
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(SPHERICAL_MULTIPOLE_SET), intent(inout) :: my_multipoles
    real(kind=DP), intent(inout)                 :: total_monopole ! rank-local!
    type(SW_EX), intent(in)                      :: swex
    type(SW_RI), intent(in)                      :: swri
    real(kind=DP), intent(in)                    :: g_coeffs(:)
    integer, intent(in)                          :: local_atom_idx
    integer, intent(in)                          :: global_atom_idx
    type(PARAL_INFO), intent(in)                 :: par
    type(ELEMENT), intent(in)                    :: elements(par%nat)

    ! jd: internal variables:
    real(kind=DP) :: J_lq, q, a
    real(kind=DP) :: factor
    integer       :: i_m, i_l, i_b, offset_of_l
    integer       :: swex_sw, swri_sw
    integer       :: species
    integer       :: offs_to_atom
    integer       :: n_mpoles_per_atom
    real(kind=DP) :: charge_on_atom
    character(len=*), parameter :: myself = 'compute_multipole_elms'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    species = elements(par%orig_atom(global_atom_idx))%species_number

    n_mpoles_per_atom = my_multipoles%n_multipoles_per_site
    offs_to_atom = (local_atom_idx-1)*n_mpoles_per_atom+1

    my_multipoles%mpoles(offs_to_atom:offs_to_atom+n_mpoles_per_atom-1) = 0D0

    do swex_sw = 1, swex%quality%num_sws_per_centre
       swri_sw = swex%swex_sw_to_swri_sw(swex_sw)
       i_b = swri%sw_idx_to_bessel_idx(swri_sw)
       i_l = swri%sphbessels(species,i_b)%lval
       q = swri%sphbessels(species,i_b)%qval
       a = swri%sphbessels(species,i_b)%aval
       i_m = swri%sw_idx_to_m(swri_sw)
       offset_of_l = (i_l+1)**2-i_l-swex%quality%min_l**2

       call utils_assert(i_l <= swex%quality%max_l, &
            myself//': Internal error (l > max_l)',i_l,swex%quality%max_l)
       call utils_assert(i_l >= swex%quality%min_l, &
            myself//': Internal error (l < min_l)',i_l,swex%quality%min_l)
       call utils_assert(i_m >= -i_l .and. i_m <= i_l, &
            myself//': Internal error (m)',i_m)

       ! jd: Room for optimisation: J_lq does not need to be
       !     recalculated if only m changes.

       ! jd: Always keep sure that this is consistent with
       !     its counterpart in dma_evaluate_multipoles().
       J_lq = swri_sph_bess_solidharm_int(i_l,q,a)
       ! jd: Employ Racah normalisation
       factor=2D0*SQRT_PI/sqrt(2D0*i_l+1D0)
       ! jd: Required fudge factor, origin unknown
       factor=factor/(2.0_DP*PI)

       ! jd: Some corrections are applied a priori
       factor = factor * pub_dma_multipole_scaling ! (default is 1.0)
       if(i_l == 1) then
          factor = factor * pub_dma_dipole_scaling
       end if
       if(i_l == 2) then
          factor = factor * pub_dma_quadrupole_scaling
       end if

       ! jd: Main formula
       my_multipoles%mpoles(offs_to_atom+offset_of_l+i_m-1) = &
            my_multipoles%mpoles(offs_to_atom+offset_of_l+i_m-1) + &
            g_coeffs(swex_sw) * factor * J_lq

    end do

    ! jd: Accum the total monopole, if there are any monopoles in the expansion
    if(my_multipoles%min_l <= 0 .and. my_multipoles%max_l >= 0) then
       charge_on_atom = my_multipoles%mpoles(offs_to_atom)
       total_monopole = total_monopole + charge_on_atom
    end if

    call timer_clock(myself,2)

  end subroutine compute_multipole_elms

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_mp_elms_for_atom(a_multipole, min_l, max_l, dma_swri_h, &
       debug_show_skipped, orig_idx, elements, &
       file_units, output_type, want_gdma_like)
    !==========================================================================!
    ! Prints all elements of a multipole vector for a single atom.             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   a_multipole (in): Multipole data for this atom.                        !
    !   min_l (in): Minimum angular momentum for the output multipoles.        !
    !   max_l (in): Maximum angular momentum for the output multipoles.        !
    !   dma_swri_h (in): Needed to skip multipoles that are outside the SWRI.  !
    !   debug_show_skipped (in): T if atoms outside of SWRI are to be printed  !
    !                            (to stdout, not to multipole files).          !
    !   orig_idx (in): Global atom index, with SF curve reordering undone.     !
    !   elements (in): The usual.                                              !
    !   file_unit_* (in): An open unit to which the data is to be written.     !
    !   output_type (in): A string to be written as header.                    !
    !   want_gdma_like (in): .true. is GDMA-like output is desired.            !
    !--------------------------------------------------------------------------!
    ! Notes (jd):                                                              !
    !   Changed from the type(MULTIPOLE) representation to a straightforward   !
    !   array (with the same indexing pattern), to make it easier for comms.   !
    !--------------------------------------------------------------------------!
    ! Current version of this subroutine is due to Jacek Dziedzic.             !
    ! Cleaned up by Jacek Dziedzic in May 2014.                                !
    ! Rewritten by Jacek Dziedzic in September 2014.                           !
    ! GDMA-like output by Jacek Dziedzic, added in March 2015.                 !
    ! Modified to remove pub_par by Joseph Prentice, September 2018            !
    !==========================================================================!

    use constants, only: ANGSTROM, stdout
    use geometry, only: POINT
    use ion, only: ELEMENT
    use multipole_ops, only: multipole_2_sph_real_to_cart_traceless, &
         multipole_2_cart_traceless_to_primitive, &
         multipole_2_cart_primitive_to_traceless, multipole_1_sph_real_to_cart,&
         multipole_1_translate_cartesian, &
         multipole_2_translate_cartesian_primitive, &
         multipole_1_cart_to_sph_real, multipole_2_cart_traceless_to_sph_real,&
         multipole_1_magnitude_cartesian, multipole_2_magnitude_cartesian, &
         multipole_1_magnitude_sph_real, multipole_2_magnitude_sph_real
    use rundat, only: pub_print_qc, pub_polarisation_simcell_refpt, &
         pub_rootname, pub_dma_precise_gdma_output, pub_dma_multipole_scaling
    use utils, only: utils_qc_print, utils_abort, utils_int_to_str

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in)      :: a_multipole(:)
    integer, intent(in)            :: min_l
    integer, intent(in)            :: max_l
    integer, intent(in)            :: dma_swri_h
    logical, intent(in)            :: debug_show_skipped
    integer, intent(in)            :: orig_idx
    type(ELEMENT), intent(in)      :: elements(:)
    integer, intent(in)            :: file_units(N_DMA_FILES)
    character(len=128)             :: output_type
    logical, intent(in)            :: want_gdma_like

    ! jd: Internal variables:
    integer                        :: i_m, i_l, offset_of_l
    real(kind=DP)                  :: refpt(3) ! =pub_polarisation_simcell_refpt
    real(kind=DP)                  :: centre_pt(3) ! origin of current multipole
    type(POINT)                    :: centre       ! same, in POINT represn'n

    ! ******* ELECTRONIC ********
    ! jd: Atom-centred spherical multipoles of electronic density
    real(kind=DP) :: Q00
    real(kind=DP) :: Q1m1, Q10, Q11
    real(kind=DP) :: Q2m2, Q2m1, Q20, Q21, Q22
    ! jd: Atom-centred Cartesian multipoles of electronic density
    real(kind=DP) :: elec_q
    real(kind=DP) :: elec_mu_X, elec_mu_Y, elec_mu_Z
    real(kind=DP) :: elec_Qxx, elec_Qxy, elec_Qxz, elec_Qyy, elec_Qyz, elec_Qzz
    real(kind=DP) :: elec_Txx, elec_Txy, elec_Txz, elec_Tyy, elec_Tyz, elec_Tzz
    ! jd: Refpt-relative Cartesian multipoles of electronic density
    real(kind=DP) :: elec_mu_X_refpt, elec_mu_Y_refpt, elec_mu_Z_refpt
    real(kind=DP) :: elec_Qxx_refpt, elec_Qxy_refpt, elec_Qxz_refpt, &
         elec_Qyy_refpt, elec_Qyz_refpt, elec_Qzz_refpt
    real(kind=DP) :: elec_Txx_refpt, elec_Txy_refpt, elec_Txz_refpt, &
         elec_Tyy_refpt, elec_Tyz_refpt, elec_Tzz_refpt
    ! jd: Refpt-relative Cartesian multipoles of electronic density, accumulated
    real(kind=DP), save :: elec_q_total, elec_mu_X_refpt_total, &
         elec_mu_Y_refpt_total, elec_mu_Z_refpt_total, elec_Qxx_refpt_total, &
         elec_Qxy_refpt_total, elec_Qxz_refpt_total, elec_Qyy_refpt_total, &
         elec_Qyz_refpt_total, elec_Qzz_refpt_total, elec_Txx_refpt_total, &
         elec_Txy_refpt_total, elec_Txz_refpt_total, elec_Tyy_refpt_total, &
         elec_Tyz_refpt_total, elec_Tzz_refpt_total

    ! ******* IONIC (CORE) ********
    real(kind=DP)       :: core_q
    real(kind=DP)       :: core_mu_X_refpt, core_mu_Y_refpt, core_mu_Z_refpt
    real(kind=DP)       :: core_Qxx_refpt, core_Qxy_refpt, core_Qxz_refpt, &
         core_Qyy_refpt, core_Qyz_refpt, core_Qzz_refpt
    real(kind=DP)       :: core_Txx_refpt, core_Txy_refpt, core_Txz_refpt, &
         core_Tyy_refpt, core_Tyz_refpt, core_Tzz_refpt

    ! ******* (ELEC+CORE) ********
    real(kind=DP)       :: elec_and_core_q
    real(kind=DP)       :: elec_and_core_mu_X_refpt, elec_and_core_mu_Y_refpt, &
         elec_and_core_mu_Z_refpt
    real(kind=DP)       :: elec_and_core_Qxx_refpt, elec_and_core_Qxy_refpt, &
         elec_and_core_Qxz_refpt, elec_and_core_Qyy_refpt, &
         elec_and_core_Qyz_refpt, elec_and_core_Qzz_refpt
    real(kind=DP)       :: elec_and_core_Txx_refpt, elec_and_core_Txy_refpt, &
         elec_and_core_Txz_refpt, elec_and_core_Tyy_refpt, &
         elec_and_core_Tyz_refpt, elec_and_core_Tzz_refpt
    real(kind=DP)       :: magnitude
    ! jd: Refpt-relative Cartesian multipoles of elec+core density, accumulated
    real(kind=DP), save :: elec_and_core_q_total, &
         elec_and_core_mu_X_refpt_total, elec_and_core_mu_Y_refpt_total, &
         elec_and_core_mu_Z_refpt_total, elec_and_core_Qxx_refpt_total, &
         elec_and_core_Qxy_refpt_total, elec_and_core_Qxz_refpt_total, &
         elec_and_core_Qyy_refpt_total, elec_and_core_Qyz_refpt_total, &
         elec_and_core_Qzz_refpt_total, elec_and_core_Txx_refpt_total, &
         elec_and_core_Txy_refpt_total, elec_and_core_Txz_refpt_total, &
         elec_and_core_Tyy_refpt_total, elec_and_core_Tyz_refpt_total, &
         elec_and_core_Tzz_refpt_total

    ! jd: Refpt-relative spherical-real multipoles of elec+core density,
    !     calculated *from* Cartesian values
    real(kind=DP) :: elec_and_core_Q10_refpt_total, &
         elec_and_core_Q11_refpt_total, elec_and_core_Q1m1_refpt_total, &
         elec_and_core_Q20_refpt_total, elec_and_core_Q21_refpt_total, &
         elec_and_core_Q2m1_refpt_total, elec_and_core_Q22_refpt_total, &
         elec_and_core_Q2m2_refpt_total

    ! jd: Varia
    integer                        :: ifile
    logical                        :: skip_this_atom
    character(len=32)              :: qc_tag
    character(len=10)              :: padded_symbol
    character(len=*), parameter    :: myself = 'print_mp_elms_for_atom'
    ! rc2013: number of atoms in system
    integer                        :: nat

    ! -------------------------------------------------------------------------

    nat=size(elements)

    refpt = pub_polarisation_simcell_refpt

    centre = elements(orig_idx)%centre
    core_q = elements(orig_idx)%ion_charge * pub_dma_multipole_scaling
    ! ^ jd: Multipole scaling is 1.0 by default. If in use, elec mpoles have
    !       already been scaled when computed (along with dcoeffs). Now cores
    !       must be scaled in the same fashion.

    ! jd: Initialize totals before first atom is output
    if(orig_idx == 1) then
       elec_q_total = 0D0
       elec_mu_X_refpt_total = 0D0
       elec_mu_Y_refpt_total = 0D0
       elec_mu_Z_refpt_total = 0D0
       elec_Qxx_refpt_total = 0D0
       elec_Qxy_refpt_total = 0D0
       elec_Qxz_refpt_total = 0D0
       elec_Qyy_refpt_total = 0D0
       elec_Qyz_refpt_total = 0D0
       elec_Qzz_refpt_total = 0D0
       elec_Txx_refpt_total = 0D0
       elec_Txy_refpt_total = 0D0
       elec_Txz_refpt_total = 0D0
       elec_Tyy_refpt_total = 0D0
       elec_Tyz_refpt_total = 0D0
       elec_Tzz_refpt_total = 0D0
       elec_and_core_q_total = 0D0
       elec_and_core_mu_X_refpt_total = 0D0
       elec_and_core_mu_Y_refpt_total = 0D0
       elec_and_core_mu_Z_refpt_total = 0D0
       elec_and_core_Qxx_refpt_total = 0D0
       elec_and_core_Qxy_refpt_total = 0D0
       elec_and_core_Qxz_refpt_total = 0D0
       elec_and_core_Qyy_refpt_total = 0D0
       elec_and_core_Qyz_refpt_total = 0D0
       elec_and_core_Qzz_refpt_total = 0D0
       elec_and_core_Txx_refpt_total = 0D0
       elec_and_core_Txy_refpt_total = 0D0
       elec_and_core_Txz_refpt_total = 0D0
       elec_and_core_Tyy_refpt_total = 0D0
       elec_and_core_Tyz_refpt_total = 0D0
       elec_and_core_Tzz_refpt_total = 0D0
    end if

    ! jd: Write output headers before first atom is output
    if(orig_idx == 1) then
       do ifile = 1, N_DMA_FILES-1 ! jd: -1: absent for GDMA file
          write(file_units(ifile),'(a)')  CRLF//trim(output_type)
       end do

       write(file_units(1),'(a)') '# Real-valued, atom-centred spherical&
            & multipoles of electronic density only (atomic units, &
            &charge convention: electrons positive).'
       write(file_units(2),'(a)') '# Atom-centred, Cartesian traceless &
            &(Gray-Gubbins convention) multipoles of electronic density only &
            &(atomic units, charge convention: electrons negative).'
       write(file_units(3),'(a)') '# Atom-centred, Cartesian primitive multipoles&
            & of electronic density only (atomic units, charge convention: &
            &electrons negative).'
       write(file_units(3),'(a)') '# NB: The resulting potential is invariant to &
            &addition of \lambda * I to the primitive quadrupole tensor, so Qxx&
            &, Qyy, Qzz may not be directly comparable with other calculations.'
       write(file_units(4),'(a,f16.8,f16.8,f16.8)') '# Cartesian primitive multi&
            &poles of electronic density only (atomic units, charge convention:&
            & electrons negative), calculated at &
            &refpt = ', refpt(1), refpt(2), refpt(3)
       write(file_units(4),'(a)') '# NB: The resulting potential is invariant to &
            &addition of \lambda * I to the primitive quadrupole tensor, so Qxx&
            &, Qyy, Qzz may not be directly comparable with other calculations.'
       write(file_units(5),'(a,f16.8,f16.8,f16.8)') '# Cartesian primitive multi&
            &poles of total (elec+core) density (atomic units, charge &
            &convention: electrons negative), calculated&
            & at refpt = ', refpt(1), refpt(2), refpt(3)
       write(file_units(5),'(a)')'# NB: The resulting potential is invariant to &
            &addition of \lambda * I to the primitive quadrupole tensor, so Qxx&
            &, Qyy, Qzz may not be directly comparable with other calculations.'
       write(file_units(6),'(a,f16.8,f16.8,f16.8)') '# Cartesian traceless &
            &(Gray-Gubbins convention, atomic units, charge convention: &
            &electrons negative) multipoles of electronic &
            &density only, calculated at refpt = ', refpt(1), refpt(2), refpt(3)
       write(file_units(7),'(a,f16.8,f16.8,f16.8)') '# Cartesian traceless &
            &(Gray-Gubbins convention, atomic units, charge convention: &
            &electrons negative) multipoles of total (elec+core) density, &
            &calculated at refpt = ', refpt(1), refpt(2), refpt(3)
       if(want_gdma_like) then
          write(file_units(8),'(a)') &
               '                                     G D M A - like output'//&
               CRLF//''//CRLF//&
               '                                by ONETEP''s implementation of'&
               //CRLF//&
               ''//CRLF//&
               '                            the Distributed Multipole Analysis &
               &method'//CRLF//&
               ''//CRLF//&
               ''//CRLF//&
               ''//CRLF//&
               ''//CRLF//&
               ''//CRLF//&
               ''//CRLF//&
               trim(pub_rootname)//&
               ''//CRLF//&
               ''//CRLF//&
               ''//CRLF//&
               ''//CRLF//&
               ''//CRLF//&
               ''//CRLF//&
               ''//CRLF//&
               '                         Distributed Multipole Analysis'//CRLF&
               //''//CRLF//&
               ''//CRLF//&
               ''//CRLF//&
               ''//CRLF//&
               'Radii have no meaning in ONETEP and will all be given as&
               & 0.500 A'//CRLF//&
               ''//CRLF//&
               'Positions and radii in angstrom'//CRLF//&
               'Multipole moments in atomic units, ea_0^k for rank k'//CRLF
       end if
    end if

    ! jd: Ignore atoms not represented in this SWRI
    if(.not. elements(orig_idx)%in_swri(dma_swri_h)) then
       if(debug_show_skipped) then
          write(stdout,'(a,i0,a,a,a)') &
               'Skipping (print_mp_elms_for_atom) orig: ', orig_idx, ' (', &
               trim(elements(orig_idx)%species_id),') -- outside of SWRI'
       end if
       skip_this_atom = .true.
    else
       skip_this_atom = .false.
    end if

    ! jd: Ignore dummy atoms (dma_dummy_atom_symbol as the 1st character of
    !     species). This allows hiding these atoms from TINKER.
#if 0
    ! --- This is now disabled, as we are not hiding dummy atoms from TINKER ---
    ! --- to preserve the correctness of the gradients.                      ---
!    if(elements(orig_idx)%symbol(1:1) == dma_dummy_atom_symbol) then
!       if(debug_show_skipped) then
!          write(stdout,'(a,i0,a,a,a)') &
!               'Skipping (print_mp_elms_for_atom) orig: ', orig_idx, ' (', &
!               trim(elements(orig_idx)%species_id),') -- a dummy atom'
!       end if
!       skip_this_atom = .true.
!    end if
#endif
    ! jd: Skip next bit of output, but ensure totals are written out
    if(skip_this_atom) goto 100

    ! --------------------------------------------------------------------------
    ! jd: Output atom-centred spherical multipoles (raw DMA values)
    ! --------------------------------------------------------------------------
    loop_l:                                                                    &
    do i_l = min_l, max_l
       offset_of_l = (i_l+1)**2-i_l-min_l**2

       loop_coeffs_m:                                                          &
       do i_m = -i_l, i_l
          write(file_units(1),'(i7, a3, f16.8, f16.8, f16.8, i2, i3, f16.8)') &
               orig_idx, elements(orig_idx)%symbol, &
               centre%X, centre%Y, centre%Z, i_l, i_m, &
               a_multipole(offset_of_l+i_m)
          if(pub_print_qc) then
             if(i_m >= 0) then
                write(qc_tag,'(i0,a,i0,a,i0)') orig_idx, '_', i_l, '_', i_m
             else
                ! jd: Work around buildbot's difficulties with tags containing '-'.
                write(qc_tag,'(i0,a,i0,a,i0)') orig_idx, '_', i_l, '_m', -i_m
             end if
             call utils_qc_print('DMA_mp_'//trim(qc_tag), &
                  a_multipole(offset_of_l+i_m),minimum_magnitude = 0.001_DP)
          end if
       end do loop_coeffs_m
    end do loop_l

    ! --------------------------------------------------------------------------
    ! jd: Output atom-centred spherical multipoles (GDMA-like output)
    ! --------------------------------------------------------------------------
    if(want_gdma_like) then
       padded_symbol = elements(orig_idx)%symbol
       write(file_units(8),'(a10,a4,f10.6,a5,f10.6,a5,f10.6,a9)') &
            padded_symbol, ' x =', &
            centre%X/ANGSTROM, '  y =', centre%Y/ANGSTROM, '  z =', &
            centre%Z/ANGSTROM, ' angstrom'
       write(file_units(8),'(a27,i1,a27)') '           Maximum rank =  ', &
            max_l,'  Radius =  0.500 angstrom'

       ! jd: Warning against monopoles corresponding to negative # of electrons
       !     The check is here, because
       !     a) it's on root, and
       !     b) it's only executed in the last dma_run, when want_gdma_like is T
       if(min_l <=0 .and. max_l >= 0) then
          offset_of_l = (0+1)**2-0-min_l**2
          if(a_multipole(offset_of_l) < 0D0) then
             write(stdout,*) &
               'SERIOUS WARNING: DMA charge on atom #'//&
               trim(utils_int_to_str(orig_idx))//&
               ' (original numbering), species '''//&
               trim(elements(orig_idx)%species_id)//&
               ''' is negative. Since this is an *electronic* rather than total &
               &charge, something must have gone very wrong, unless your system co&
               &ntains positrons. If so, ONETEP does not support antimatter calcul&
               &ations yet. Seriously, though, consider using dma_bessel_averaging&
               & T (if not using it already) and/or increasing dma_max_q (to ~15-2&
               &0).'//CRLF
          end if

          ! jd: Take core charges into account. Core dipoles, quadrupoles etc.
          !     vanish, since we're in the atom-centred frame of reference.
          Q00 = a_multipole(offset_of_l) - core_q
       else
          Q00 = -core_q
       end if
       if(abs(Q00) >= 1000.0_DP) then
          call utils_abort(myself//': Atom-centred monopole exceeded +/-1000. &
               &This would produce a non-conforming DMA-like output file and &
               &is likely highly unphysical. This might indicate serious &
               &problems with your metric matrix of the SWRI procedure in &
               &general. The value of the monopole follows',&
               opt_real_to_print1=Q00)
       end if
       if(min_l <= 0 .and. max_l >= 0) then
          if(pub_dma_precise_gdma_output) then
             write(file_units(8),'(19x, a, f20.15)') 'Q00  =', -Q00 ! ONETEP uses
          else                                                   ! reverse charge
             write(file_units(8),'(19x, a, f11.6)') 'Q00  =', -Q00 ! convention
          end if
       end if
       if(min_l <= 1 .and. max_l >= 1) then
          offset_of_l = (1+1)**2-1-min_l**2
          Q1m1 = a_multipole(offset_of_l-1)
          Q10 =  a_multipole(offset_of_l  )
          Q11 =  a_multipole(offset_of_l+1)
          magnitude = multipole_1_magnitude_sph_real(Q10,Q1m1,Q11)
          if(pub_dma_precise_gdma_output) then
             if(magnitude >= 1D3) then
                write(file_units(8),'(a6,1p,e11.3,a8,e20.12,a8,e20.12,a8,e20.12)') &
                '|Q1| =', magnitude, &
                '  Q10  =', Q10, '  Q11c =', Q11, '  Q11s =', Q1m1
             else
                 write(file_units(8),'(a6,f11.6,a8,f20.15,a8,f20.15,a8,f20.15)') &
                      '|Q1| =', magnitude, &
                      '  Q10  =', Q10, '  Q11c =', Q11, '  Q11s =', Q1m1
             end if
          else
             if(magnitude >= 1D3) then
                write(file_units(8),'(a6,1p,e11.3,a8,e11.3,a8,e11.3,a8,e11.3)') &
                '|Q1| =', magnitude, &
                '  Q10  =', Q10, '  Q11c =', Q11, '  Q11s =', Q1m1
             else
                write(file_units(8),'(a6,f11.6,a8,f11.6,a8,f11.6,a8,f11.6)') &
                '|Q1| =', magnitude, &
                '  Q10  =', Q10, '  Q11c =', Q11, '  Q11s =', Q1m1
             end if
          end if
       end if
       if(min_l <= 2 .and. max_l >= 2) then
          offset_of_l = (2+1)**2-2-min_l**2

          Q2m2 = -a_multipole(offset_of_l-2) ! @ In the future these minuses will disappear
          Q2m1 = -a_multipole(offset_of_l-1) !   as (-)^l is moved from the potential
          Q20  = -a_multipole(offset_of_l  ) !   to the multipole definition. For now
          Q21  = -a_multipole(offset_of_l+1) !   the minuses remain here to satisfy
          Q22  = -a_multipole(offset_of_l+2) !   TINKER
          magnitude = multipole_2_magnitude_sph_real(Q20,Q2m1,Q21,Q2m2,Q22)

          if(pub_dma_precise_gdma_output) then
             if(magnitude >= 1D3) then
                write(file_units(8),'(a6,1p,e11.3,a8,e20.12,a8,e20.12,a8,e20.12)') &
                '|Q2| =', magnitude, &
                '  Q20  =', Q20, '  Q21c =', Q21, '  Q21s =', Q2m1
                write(file_units(8),'(17x,a8,1p,e20.12,a8,e20.12)') &
                '  Q22c =', Q22, '  Q22s =', Q2m2
             else
                write(file_units(8),'(a6,f11.6,a8,f20.15,a8,f20.15,a8,f20.15)') &
                '|Q2| =', magnitude, &
                '  Q20  =', Q20, '  Q21c =', Q21, '  Q21s =', Q2m1
                write(file_units(8),'(17x,a8,f20.15,a8,f20.15)') &
                '  Q22c =', Q22, '  Q22s =', Q2m2
             end if
          else
             if(magnitude >= 1D3) then
                write(file_units(8),'(a6,1p,e11.3,a8,e11.3,a8,e11.3,a8,e11.3)') &
                '|Q2| =', magnitude, &
                '  Q20  =', Q20, '  Q21c =', Q21, '  Q21s =', Q2m1
                write(file_units(8),'(17x,a8,1p,e11.3,a8,e11.3)') &
                '  Q22c =', Q22, '  Q22s =', Q2m2
             else
                write(file_units(8),'(a6,f11.6,a8,f11.6,a8,f11.6,a8,f11.6)') &
                '|Q2| =', magnitude, &
                '  Q20  =', Q20, '  Q21c =', Q21, '  Q21s =', Q2m1
                write(file_units(8),'(17x,a8,f11.6,a8,f11.6)') &
                '  Q22c =', Q22, '  Q22s =', Q2m2
             end if
          end if
       end if
       write(file_units(8),'(a)') ''
    end if

    ! jd: Calculate dipole moment from core, at refpt
    if(min_l <= 1 .and. max_l >= 1) then
       core_mu_X_refpt = core_q * (centre%X-refpt(1))
       core_mu_Y_refpt = core_q * (centre%Y-refpt(2))
       core_mu_Z_refpt = core_q * (centre%Z-refpt(3))
    else
       core_mu_X_refpt = 0D0
       core_mu_Y_refpt = 0D0
       core_mu_Z_refpt = 0D0
    end if

    ! jd: Calculate quadrupole moment (Cartesian primitive) from core, at refpt
    if(min_l <= 2 .and. max_l >= 2) then
       core_Qxx_refpt = core_q * (centre%X-refpt(1)) * (centre%X-refpt(1))
       core_Qxy_refpt = core_q * (centre%X-refpt(1)) * (centre%Y-refpt(2))
       core_Qxz_refpt = core_q * (centre%X-refpt(1)) * (centre%Z-refpt(3))
       core_Qyy_refpt = core_q * (centre%Y-refpt(2)) * (centre%Y-refpt(2))
       core_Qyz_refpt = core_q * (centre%Y-refpt(2)) * (centre%Z-refpt(3))
       core_Qzz_refpt = core_q * (centre%Z-refpt(3)) * (centre%Z-refpt(3))
    else
       core_Qxx_refpt = 0D0
       core_Qxy_refpt = 0D0
       core_Qxz_refpt = 0D0
       core_Qyy_refpt = 0D0
       core_Qyz_refpt = 0D0
       core_Qzz_refpt = 0D0
    end if

    ! jd: Calculate Racah-normalized real-spherical, atom-centred multipoles
    !     of electronic density. Set higher multipoles to zero, if DMA order
    !     was too low.
    if(min_l <= 0 .and. max_l >= 0) then
       offset_of_l = (0+1)**2-0-min_l**2
       Q00 = a_multipole(offset_of_l)
    else
       Q00 = 0D0
    end if

    if(min_l <= 1 .and. max_l >= 1) then
       offset_of_l = (1+1)**2-1-min_l**2
       Q1m1 = a_multipole(offset_of_l-1)
       Q10  = a_multipole(offset_of_l  )
       Q11  = a_multipole(offset_of_l+1)
    else
       Q1m1 = 0D0
       Q10 = 0D0
       Q11 = 0D0
    end if

    if(min_l <= 2 .and. max_l >= 2) then
       offset_of_l = (2+1)**2-2-min_l**2
       Q2m2 = a_multipole(offset_of_l-2)
       Q2m1 = a_multipole(offset_of_l-1)
       Q20  = a_multipole(offset_of_l  )
       Q21  = a_multipole(offset_of_l+1)
       Q22  = a_multipole(offset_of_l+2)
    else
       Q2m2 = 0D0
       Q2m1 = 0D0
       Q20 = 0D0
       Q21 = 0D0
       Q22 = 0D0
    end if

    ! jd: Internally we keep elec_q in the electrons-positive convention,
    !     and only flip for output.
    elec_q = Q00
    elec_and_core_q = elec_q - core_q

    ! jd: Convert electronic multipoles to traceless Cartesian,
    !     still atom-centred
    call multipole_1_sph_real_to_cart(elec_mu_X, elec_mu_Y, elec_mu_Z, &
         Q1m1, Q10, Q11, .true.)
    call multipole_2_sph_real_to_cart_traceless(&
         elec_Txx, elec_Txy, elec_Txz, elec_Tyy, elec_Tyz, elec_Tzz, &
         Q2m2, Q2m1, Q20, Q21, Q22, .true.)

    ! --------------------------------------------------------------------------
    ! jd: Output atom-centred traceless Cartesian multipoles
    ! --------------------------------------------------------------------------
    if(min_l <= 0 .and. max_l >= 0) then
       write(file_units(2),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' q =   ', -elec_q
    end if
    if(min_l <= 1 .and. max_l >= 1) then
       write(file_units(2),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' dx =  ', elec_mu_X
       write(file_units(2),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' dy =  ', elec_mu_Y
       write(file_units(2),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' dz =  ', elec_mu_Z
    end if
    if(min_l <= 2 .and. max_l >= 2) then
       write(file_units(2),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qxx = ', elec_Txx
       write(file_units(2),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qxy = ', elec_Txy
       write(file_units(2),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qxz = ', elec_Txz
       write(file_units(2),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qyy = ', elec_Tyy
       write(file_units(2),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qyz = ', elec_Tyz
       write(file_units(2),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qzz = ', elec_Tzz
    end if

    ! Convert to Cartesian primitive, still atom-centred
    call multipole_2_cart_traceless_to_primitive(&
         elec_Qxx, elec_Qxy, elec_Qxz, elec_Qyy, elec_Qyz, elec_Qzz, &
         elec_Txx, elec_Txy, elec_Txz, elec_Tyy, elec_Tyz, elec_Tzz, 0.0_DP)

    ! --------------------------------------------------------------------------
    ! jd: Output atom-centred primitive Cartesian multipoles
    ! --------------------------------------------------------------------------
    if(min_l <= 0 .and. max_l >= 0) then
       write(file_units(3),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' q =   ', -elec_q
    end if
    if(min_l <= 1 .and. max_l >= 1) then
       write(file_units(3),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' dx =  ', elec_mu_X
       write(file_units(3),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' dy =  ', elec_mu_Y
       write(file_units(3),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' dz =  ', elec_mu_Z
    end if
    if(min_l <= 2 .and. max_l >= 2) then
       write(file_units(3),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qxx = ', elec_Qxx
       write(file_units(3),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qxy = ', elec_Qxy
       write(file_units(3),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qxz = ', elec_Qxz
       write(file_units(3),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qyy = ', elec_Qyy
       write(file_units(3),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qyz = ', elec_Qyz
       write(file_units(3),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qzz = ', elec_Qzz
    end if

    ! jd: Translate Cartesian primitive to refpt
    centre_pt(1) = centre%X
    centre_pt(2) = centre%Y
    centre_pt(3) = centre%Z

    ! - dipoles
    if(min_l <= 1 .and. max_l >= 1) then
       call multipole_1_translate_cartesian(elec_mu_X_refpt, &
            elec_mu_X,elec_q,centre_pt,refpt,1)
       call multipole_1_translate_cartesian(elec_mu_Y_refpt, &
            elec_mu_Y,elec_q,centre_pt,refpt,2)
       call multipole_1_translate_cartesian(elec_mu_Z_refpt, &
            elec_mu_Z,elec_q,centre_pt,refpt,3)
    else
       elec_mu_X_refpt = 0D0
       elec_mu_Y_refpt = 0D0
       elec_mu_Z_refpt = 0D0
    end if

    ! - quadrupoles
    if(max_l <= 2 .and. max_l >= 2) then
       call multipole_2_translate_cartesian_primitive(elec_Qxx_refpt, &
            elec_Qxx, elec_mu_X_refpt, elec_mu_X_refpt, elec_q, &
            centre_pt, refpt, 1, 1)
       call multipole_2_translate_cartesian_primitive(elec_Qxy_refpt, &
            elec_Qxy, elec_mu_X_refpt, elec_mu_Y_refpt, elec_q, &
            centre_pt, refpt, 1, 2)
       call multipole_2_translate_cartesian_primitive(elec_Qxz_refpt, &
            elec_Qxz, elec_mu_X_refpt, elec_mu_Z_refpt, elec_q, &
            centre_pt, refpt, 1, 3)
       call multipole_2_translate_cartesian_primitive(elec_Qyy_refpt, &
            elec_Qyy, elec_mu_Y_refpt, elec_mu_Y_refpt, elec_q, &
            centre_pt, refpt, 2, 2)
       call multipole_2_translate_cartesian_primitive(elec_Qyz_refpt, &
            elec_Qyz, elec_mu_Y_refpt, elec_mu_Z_refpt, elec_q, &
            centre_pt, refpt, 2, 3)
       call multipole_2_translate_cartesian_primitive(elec_Qzz_refpt, &
            elec_Qzz, elec_mu_Z_refpt, elec_mu_Z_refpt, elec_q, &
            centre_pt, refpt, 3, 3)
    else
       elec_Qxx_refpt = 0D0
       elec_Qxy_refpt = 0D0
       elec_Qxz_refpt = 0D0
       elec_Qyy_refpt = 0D0
       elec_Qyz_refpt = 0D0
       elec_Qzz_refpt = 0D0
    end if

    ! jd: Convert electronic quadrupoles to traceless (Gray-Gubbins convention)
    call multipole_2_cart_primitive_to_traceless(&
         elec_Txx_refpt, elec_Txy_refpt, elec_Txz_refpt, &
         elec_Tyy_refpt, elec_Tyz_refpt, elec_Tzz_refpt, &
         elec_Qxx_refpt, elec_Qxy_refpt, elec_Qxz_refpt, &
         elec_Qyy_refpt, elec_Qyz_refpt, elec_Qzz_refpt, &
         'GrayGubbins')

    ! jd: Calculate traceless Cartesian quadrupoles of the ions
    call multipole_2_cart_primitive_to_traceless(&
         core_Txx_refpt, core_Txy_refpt, core_Txz_refpt, &  ! output
         core_Tyy_refpt, core_Tyz_refpt, core_Tzz_refpt, &  ! output
         core_Qxx_refpt, core_Qxy_refpt, core_Qxz_refpt, &  ! input
         core_Qyy_refpt, core_Qyz_refpt, core_Qzz_refpt, &  ! input
         & 'GrayGubbins')

    ! jd: Calculate elec+core multipoles
    elec_and_core_mu_X_refpt = elec_mu_X_refpt + core_mu_X_refpt
    elec_and_core_mu_Y_refpt = elec_mu_Y_refpt + core_mu_Y_refpt
    elec_and_core_mu_Z_refpt = elec_mu_Z_refpt + core_mu_Z_refpt
    elec_and_core_Qxx_refpt = elec_Qxx_refpt + core_Qxx_refpt
    elec_and_core_Qxy_refpt = elec_Qxy_refpt + core_Qxy_refpt
    elec_and_core_Qxz_refpt = elec_Qxz_refpt + core_Qxz_refpt
    elec_and_core_Qyy_refpt = elec_Qyy_refpt + core_Qyy_refpt
    elec_and_core_Qyz_refpt = elec_Qyz_refpt + core_Qyz_refpt
    elec_and_core_Qzz_refpt = elec_Qzz_refpt + core_Qzz_refpt
    elec_and_core_Txx_refpt = elec_Txx_refpt + core_Txx_refpt
    elec_and_core_Txy_refpt = elec_Txy_refpt + core_Txy_refpt
    elec_and_core_Txz_refpt = elec_Txz_refpt + core_Txz_refpt
    elec_and_core_Tyy_refpt = elec_Tyy_refpt + core_Tyy_refpt
    elec_and_core_Tyz_refpt = elec_Tyz_refpt + core_Tyz_refpt
    elec_and_core_Tzz_refpt = elec_Tzz_refpt + core_Tzz_refpt

    ! --------------------------------------------------------------------------
    ! jd: Output refpt-relative primitive Cartesian multipoles
    ! --------------------------------------------------------------------------
    if(min_l <= 0 .and. max_l >= 0) then
       write(file_units(4),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' q =   ', -elec_q
    end if
    if(min_l <= 1 .and. max_l >= 1) then
       write(file_units(4),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' dx =  ', elec_mu_X_refpt
       write(file_units(4),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' dy =  ', elec_mu_Y_refpt
       write(file_units(4),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' dz =  ', elec_mu_Z_refpt
    end if
    if(min_l <= 2 .and. max_l >= 2) then
       write(file_units(4),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qxx = ', elec_Qxx_refpt
       write(file_units(4),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qxy = ', elec_Qxy_refpt
       write(file_units(4),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qxz = ', elec_Qxz_refpt
       write(file_units(4),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qyy = ', elec_Qyy_refpt
       write(file_units(4),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qyz = ', elec_Qyz_refpt
       write(file_units(4),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qzz = ', elec_Qzz_refpt
    end if

    if(min_l <= 0 .and. max_l >= 0) then
       write(file_units(5),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' q =   ', elec_and_core_q
    end if
    if(min_l <= 1 .and. max_l >= 1) then
       write(file_units(5),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' dx =  ', elec_and_core_mu_X_refpt
       write(file_units(5),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' dy =  ', elec_and_core_mu_Y_refpt
       write(file_units(5),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' dz =  ', elec_and_core_mu_Z_refpt
    end if
    if(min_l <= 2 .and. max_l >= 2) then
       write(file_units(5),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qxx = ', elec_and_core_Qxx_refpt
       write(file_units(5),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qxy = ', elec_and_core_Qxy_refpt
       write(file_units(5),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qxz = ', elec_and_core_Qxz_refpt
       write(file_units(5),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qyy = ', elec_and_core_Qyy_refpt
       write(file_units(5),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qyz = ', elec_and_core_Qyz_refpt
       write(file_units(5),'(i7, a3, f16.8, f16.8, f16.8, a7, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, ' Qzz = ', elec_and_core_Qzz_refpt
    end if

    ! jd: Accumulate total electronic monopole, dipole, primitive quadrupole
    elec_q_total = elec_q_total + elec_q
    elec_mu_X_refpt_total = elec_mu_X_refpt_total + elec_mu_X_refpt
    elec_mu_Y_refpt_total = elec_mu_Y_refpt_total + elec_mu_Y_refpt
    elec_mu_Z_refpt_total = elec_mu_Z_refpt_total + elec_mu_Z_refpt
    elec_Qxx_refpt_total = elec_Qxx_refpt_total + elec_Qxx_refpt
    elec_Qxy_refpt_total = elec_Qxy_refpt_total + elec_Qxy_refpt
    elec_Qxz_refpt_total = elec_Qxz_refpt_total + elec_Qxz_refpt
    elec_Qyy_refpt_total = elec_Qyy_refpt_total + elec_Qyy_refpt
    elec_Qyz_refpt_total = elec_Qyz_refpt_total + elec_Qyz_refpt
    elec_Qzz_refpt_total = elec_Qzz_refpt_total + elec_Qzz_refpt

    ! jd: Accumulate total elec+core monopole, dipole, primitive quadrupole
    elec_and_core_q_total = elec_and_core_q_total + elec_and_core_q
    elec_and_core_mu_X_refpt_total = &
         elec_and_core_mu_X_refpt_total + elec_and_core_mu_X_refpt
    elec_and_core_mu_Y_refpt_total = &
         elec_and_core_mu_Y_refpt_total + elec_and_core_mu_Y_refpt
    elec_and_core_mu_Z_refpt_total = &
         elec_and_core_mu_Z_refpt_total + elec_and_core_mu_Z_refpt
    elec_and_core_Qxx_refpt_total = &
         elec_and_core_Qxx_refpt_total + elec_and_core_Qxx_refpt
    elec_and_core_Qxy_refpt_total = &
         elec_and_core_Qxy_refpt_total + elec_and_core_Qxy_refpt
    elec_and_core_Qxz_refpt_total = &
         elec_and_core_Qxz_refpt_total + elec_and_core_Qxz_refpt
    elec_and_core_Qyy_refpt_total = &
         elec_and_core_Qyy_refpt_total + elec_and_core_Qyy_refpt
    elec_and_core_Qyz_refpt_total = &
         elec_and_core_Qyz_refpt_total + elec_and_core_Qyz_refpt
    elec_and_core_Qzz_refpt_total = &
         elec_and_core_Qzz_refpt_total + elec_and_core_Qzz_refpt

100 continue

    ! jd: Output totals for refpt-relative primitive Cartesian multipoles
    if(orig_idx == nat) then
       if(min_l <= 0 .and. max_l >= 0) then
          write(file_units(4),'(a5, a7, f16.8)') &
               'TOTAL', ' q =   ', -elec_q_total
       end if
       if(min_l <= 1 .and. max_l >= 1) then
          write(file_units(4),'(a5, a7, f16.8)') &
               'TOTAL', ' dx =  ', elec_mu_X_refpt_total
          write(file_units(4),'(a5, a7, f16.8)') &
               'TOTAL', ' dy =  ', elec_mu_Y_refpt_total
          write(file_units(4),'(a5, a7, f16.8)') &
               'TOTAL', ' dz =  ', elec_mu_Z_refpt_total
          write(file_units(4),'(a5, a7, f16.8)') &
               'TOTAL', ' |d| = ', multipole_1_magnitude_cartesian(&
               elec_mu_X_refpt_total, elec_mu_Y_refpt_total, elec_mu_Z_refpt_total)
       end if
       if(min_l <= 2 .and. max_l >= 2) then
          write(file_units(4),'(a5, a7, f16.8)') &
               'TOTAL', ' Qxx = ', elec_Qxx_refpt_total
          write(file_units(4),'(a5, a7, f16.8)') &
               'TOTAL', ' Qxy = ', elec_Qxy_refpt_total
          write(file_units(4),'(a5, a7, f16.8)') &
               'TOTAL', ' Qxz = ', elec_Qxz_refpt_total
          write(file_units(4),'(a5, a7, f16.8)') &
               'TOTAL', ' Qyy = ', elec_Qyy_refpt_total
          write(file_units(4),'(a5, a7, f16.8)') &
               'TOTAL', ' Qyz = ', elec_Qyz_refpt_total
          write(file_units(4),'(a5, a7, f16.8)') &
               'TOTAL', ' Qzz = ', elec_Qzz_refpt_total
          write(file_units(4),'(a5, a7, f16.8)') &
               'TOTAL', ' |Q| = ', multipole_2_magnitude_cartesian(&
               elec_Qxx_refpt_total,elec_Qxy_refpt_total,elec_Qxz_refpt_total,&
               elec_Qyy_refpt_total,elec_Qyz_refpt_total,elec_Qzz_refpt_total)
       end if

       if(min_l <= 0 .and. max_l >= 0) then
          write(file_units(5),'(a5, a7, f16.8)') &
               'TOTAL', ' q =   ', -elec_and_core_q_total
       end if
       if(min_l <= 1 .and. max_l >= 1) then
          write(file_units(5),'(a5, a7, f16.8)') &
               'TOTAL', ' dx =  ', elec_and_core_mu_X_refpt_total
          write(file_units(5),'(a5, a7, f16.8)') &
               'TOTAL', ' dy =  ', elec_and_core_mu_Y_refpt_total
          write(file_units(5),'(a5, a7, f16.8)') &
               'TOTAL', ' dz =  ', elec_and_core_mu_Z_refpt_total
          write(file_units(5),'(a5, a7, f16.8)') &
               'TOTAL', ' |d| = ', multipole_1_magnitude_cartesian(&
               elec_and_core_mu_X_refpt_total, elec_and_core_mu_Y_refpt_total, &
               elec_and_core_mu_Z_refpt_total)
       end if
       if(min_l <= 2 .and. max_l >= 2) then
          write(file_units(5),'(a5, a7, f16.8)') &
               'TOTAL', ' Qxx = ', elec_and_core_Qxx_refpt_total
          write(file_units(5),'(a5, a7, f16.8)') &
               'TOTAL', ' Qxy = ', elec_and_core_Qxy_refpt_total
          write(file_units(5),'(a5, a7, f16.8)') &
               'TOTAL', ' Qxz = ', elec_and_core_Qxz_refpt_total
          write(file_units(5),'(a5, a7, f16.8)') &
               'TOTAL', ' Qyy = ', elec_and_core_Qyy_refpt_total
          write(file_units(5),'(a5, a7, f16.8)') &
               'TOTAL', ' Qyz = ', elec_and_core_Qyz_refpt_total
          write(file_units(5),'(a5, a7, f16.8)') &
               'TOTAL', ' Qzz = ', elec_and_core_Qzz_refpt_total
          write(file_units(5),'(a5, a7, f16.8)') &
               'TOTAL', ' |Q| = ', multipole_2_magnitude_cartesian(&
               elec_and_core_Qxx_refpt_total,elec_and_core_Qxy_refpt_total,&
               elec_and_core_Qxz_refpt_total,elec_and_core_Qyy_refpt_total,&
               elec_and_core_Qyz_refpt_total,elec_and_core_Qzz_refpt_total)
       end if
    end if ! last atom

    ! jd: Skip next bit of output, but ensure totals are written out
    if(skip_this_atom) goto 200

    ! --------------------------------------------------------------------------
    ! jd: Output refpt-relative traceless Cartesian multipoles
    ! --------------------------------------------------------------------------
    if(min_l <= 0 .and. max_l >= 0) then
       write(file_units(6),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  q =   ', -elec_q
    end if
    if(min_l <= 1 .and. max_l >= 1) then
       write(file_units(6),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  dx =  ', elec_mu_X_refpt
       write(file_units(6),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  dy =  ', elec_mu_Y_refpt
       write(file_units(6),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  dz =  ', elec_mu_Z_refpt
    end if
    if(min_l <= 2 .and. max_l >= 2) then
       write(file_units(6),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  Qxx = ', elec_Txx_refpt
       write(file_units(6),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  Qxy = ', elec_Txy_refpt
       write(file_units(6),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  Qxz = ', elec_Txz_refpt
       write(file_units(6),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  Qyy = ', elec_Tyy_refpt
       write(file_units(6),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  Qyz = ', elec_Tyz_refpt
       write(file_units(6),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  Qzz = ', elec_Tzz_refpt
    end if

    if(min_l <= 0 .and. max_l >= 0) then
       write(file_units(7),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  q =   ', -elec_and_core_q
    end if
    if(min_l <= 1 .and. max_l >= 1) then
       write(file_units(7),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  dx =  ', elec_and_core_mu_X_refpt
       write(file_units(7),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  dy =  ', elec_and_core_mu_Y_refpt
       write(file_units(7),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  dz =  ', elec_and_core_mu_Z_refpt
    end if
    if(min_l <= 2 .and. max_l >= 2) then
       write(file_units(7),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  Qxx = ', elec_and_core_Txx_refpt
       write(file_units(7),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  Qxy = ', elec_and_core_Txy_refpt
       write(file_units(7),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  Qxz = ', elec_and_core_Txz_refpt
       write(file_units(7),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  Qyy = ', elec_and_core_Tyy_refpt
       write(file_units(7),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  Qyz = ', elec_and_core_Tyz_refpt
       write(file_units(7),'(i7, a3, f16.8, f16.8, f16.8, a8, f16.8)') &
            orig_idx, elements(orig_idx)%symbol, &
            centre%X, centre%Y, centre%Z, '  Qzz = ', elec_and_core_Tzz_refpt
    end if

    ! jd: Accumulate total electronic traceless quadrupole
    elec_Txx_refpt_total = elec_Txx_refpt_total + elec_Txx_refpt
    elec_Txy_refpt_total = elec_Txy_refpt_total + elec_Txy_refpt
    elec_Txz_refpt_total = elec_Txz_refpt_total + elec_Txz_refpt
    elec_Tyy_refpt_total = elec_Tyy_refpt_total + elec_Tyy_refpt
    elec_Tyz_refpt_total = elec_Tyz_refpt_total + elec_Tyz_refpt
    elec_Tzz_refpt_total = elec_Tzz_refpt_total + elec_Tzz_refpt

    ! jd: Accumulate total elec+core traceless quadrupole
    elec_and_core_Txx_refpt_total = &
         elec_and_core_Txx_refpt_total + elec_and_core_Txx_refpt
    elec_and_core_Txy_refpt_total = &
         elec_and_core_Txy_refpt_total + elec_and_core_Txy_refpt
    elec_and_core_Txz_refpt_total = &
         elec_and_core_Txz_refpt_total + elec_and_core_Txz_refpt
    elec_and_core_Tyy_refpt_total = &
         elec_and_core_Tyy_refpt_total + elec_and_core_Tyy_refpt
    elec_and_core_Tyz_refpt_total = &
         elec_and_core_Tyz_refpt_total + elec_and_core_Tyz_refpt
    elec_and_core_Tzz_refpt_total = &
         elec_and_core_Tzz_refpt_total + elec_and_core_Tzz_refpt

200 continue

    ! jd: Output totals for refpt-relative traceless Cartesian multipoles
    if(orig_idx == nat) then
       if(min_l <= 0 .and. max_l >= 0) then
          write(file_units(6),'(a5, a7, f16.8)') &
               'TOTAL', ' q =   ', -elec_q_total
       end if
       if(min_l <= 1 .and. max_l >= 1) then
          write(file_units(6),'(a5, a7, f16.8)') &
               'TOTAL', ' dx =  ', elec_mu_X_refpt_total
          write(file_units(6),'(a5, a7, f16.8)') &
               'TOTAL', ' dy =  ', elec_mu_Y_refpt_total
          write(file_units(6),'(a5, a7, f16.8)') &
               'TOTAL', ' dz =  ', elec_mu_Z_refpt_total
          write(file_units(6),'(a5, a7, f16.8)') &
               'TOTAL', ' |d| = ', multipole_1_magnitude_cartesian(&
               elec_mu_X_refpt_total,elec_mu_Y_refpt_total,elec_mu_Z_refpt_total)
       end if
       if(min_l <= 2 .and. max_l >= 2) then
          write(file_units(6),'(a5, a7, f16.8)') &
               'TOTAL', ' Qxx = ', elec_Txx_refpt_total
          write(file_units(6),'(a5, a7, f16.8)') &
               'TOTAL', ' Qxy = ', elec_Txy_refpt_total
          write(file_units(6),'(a5, a7, f16.8)') &
               'TOTAL', ' Qxz = ', elec_Txz_refpt_total
          write(file_units(6),'(a5, a7, f16.8)') &
               'TOTAL', ' Qyy = ', elec_Tyy_refpt_total
          write(file_units(6),'(a5, a7, f16.8)') &
               'TOTAL', ' Qyz = ', elec_Tyz_refpt_total
          write(file_units(6),'(a5, a7, f16.8)') &
               'TOTAL', ' Qzz = ', elec_Tzz_refpt_total
          write(file_units(6),'(a5, a7, f16.8)') &
               'TOTAL', ' |Q| = ', multipole_2_magnitude_cartesian(&
               elec_Txx_refpt_total,elec_Txy_refpt_total,elec_Txz_refpt_total,&
               elec_Tyy_refpt_total,elec_Tyz_refpt_total,elec_Tzz_refpt_total)
       end if

       if(min_l <= 0 .and. max_l >= 0) then
          write(file_units(7),'(a5, a7, f16.8)') &
               'TOTAL', ' q =   ', -elec_and_core_q_total
       end if
       if(min_l <= 1 .and. max_l >= 1) then
          write(file_units(7),'(a5, a7, f16.8)') &
               'TOTAL', ' dx =  ', elec_and_core_mu_X_refpt_total
          write(file_units(7),'(a5, a7, f16.8)') &
               'TOTAL', ' dy =  ', elec_and_core_mu_Y_refpt_total
          write(file_units(7),'(a5, a7, f16.8)') &
               'TOTAL', ' dz =  ', elec_and_core_mu_Z_refpt_total
          write(file_units(7),'(a5, a7, f16.8)') &
               'TOTAL', ' |d| = ', multipole_1_magnitude_cartesian(&
               elec_and_core_mu_X_refpt_total,elec_and_core_mu_Y_refpt_total, &
               elec_and_core_mu_Z_refpt_total)
       end if
       if(min_l <= 2 .and. max_l >= 2) then
          write(file_units(7),'(a5, a7, f16.8)') &
               'TOTAL', ' Qxx = ', elec_and_core_Txx_refpt_total
          write(file_units(7),'(a5, a7, f16.8)') &
               'TOTAL', ' Qxy = ', elec_and_core_Txy_refpt_total
          write(file_units(7),'(a5, a7, f16.8)') &
               'TOTAL', ' Qxz = ', elec_and_core_Txz_refpt_total
          write(file_units(7),'(a5, a7, f16.8)') &
               'TOTAL', ' Qyy = ', elec_and_core_Tyy_refpt_total
          write(file_units(7),'(a5, a7, f16.8)') &
               'TOTAL', ' Qyz = ', elec_and_core_Tyz_refpt_total
          write(file_units(7),'(a5, a7, f16.8)') &
               'TOTAL', ' Qzz = ', elec_and_core_Tzz_refpt_total
          write(file_units(7),'(a5, a7, f16.8)') &
               'TOTAL', ' |Q| = ', multipole_2_magnitude_cartesian(&
               elec_and_core_Txx_refpt_total,elec_and_core_Txy_refpt_total, &
               elec_and_core_Txz_refpt_total,elec_and_core_Tyy_refpt_total, &
               elec_and_core_Tyz_refpt_total,elec_and_core_Tzz_refpt_total)
       end if

       call multipole_1_cart_to_sph_real(&
            elec_and_core_Q1m1_refpt_total, &
            elec_and_core_Q10_refpt_total, &
            elec_and_core_Q11_refpt_total, &
            elec_and_core_mu_X_refpt_total, &
            elec_and_core_mu_Y_refpt_total, &
            elec_and_core_mu_Z_refpt_total, .true.)

       if(min_l <= 2 .and. max_l >= 1) then
          write(file_units(7),'(a)') &
               CRLF//'The above in spherical-real (GDMA-like) representation:'
       end if

       if(min_l <= 1 .and. max_l >= 1) then
          write(file_units(7),'(a5, a9, f16.8)') &
               'TOTAL', ' Q10  =  ', elec_and_core_Q10_refpt_total
          write(file_units(7),'(a5, a9, f16.8)') &
               'TOTAL', ' Q11c  =  ', elec_and_core_Q11_refpt_total
          write(file_units(7),'(a5, a9, f16.8)') &
               'TOTAL', ' Q11s  =  ', elec_and_core_Q1m1_refpt_total
          write(file_units(7),'(a5, a9, f16.8)') &
               'TOTAL', ' |Q1|  =  ', multipole_1_magnitude_sph_real(&
               elec_and_core_Q10_refpt_total, &
               elec_and_core_Q11_refpt_total, &
               elec_and_core_Q1m1_refpt_total)
       end if

       call multipole_2_cart_traceless_to_sph_real(&
            elec_and_core_Q2m2_refpt_total, &
            elec_and_core_Q2m1_refpt_total, &
            elec_and_core_Q20_refpt_total, &
            elec_and_core_Q21_refpt_total, &
            elec_and_core_Q22_refpt_total, &
            elec_and_core_Txx_refpt_total, &
            elec_and_core_Txy_refpt_total, &
            elec_and_core_Txz_refpt_total, &
            elec_and_core_Tyy_refpt_total, &
            elec_and_core_Tyz_refpt_total, &
            elec_and_core_Tzz_refpt_total, .true.)

       if(min_l <= 2 .and. max_l >= 2) then
          write(file_units(7),'(a5, a9, f16.8)') &
               'TOTAL', ' Q20  =  ', elec_and_core_Q20_refpt_total
          write(file_units(7),'(a5, a9, f16.8)') &
               'TOTAL', ' Q21c  =  ', elec_and_core_Q21_refpt_total
          write(file_units(7),'(a5, a9, f16.8)') &
               'TOTAL', ' Q21s  =  ', elec_and_core_Q2m1_refpt_total
          write(file_units(7),'(a5, a9, f16.8)') &
               'TOTAL', ' Q22c  =  ', elec_and_core_Q22_refpt_total
          write(file_units(7),'(a5, a9, f16.8)') &
               'TOTAL', ' Q22s  =  ', elec_and_core_Q2m2_refpt_total
          write(file_units(7),'(a5, a9, f16.8)') &
               'TOTAL', ' |Q2|  =  ', multipole_2_magnitude_sph_real(&
               elec_and_core_Q20_refpt_total, &
               elec_and_core_Q21_refpt_total, &
               elec_and_core_Q2m1_refpt_total, &
               elec_and_core_Q22_refpt_total, &
               elec_and_core_Q2m2_refpt_total)
       end if

       if(want_gdma_like) then
          write(file_units(8),'(a38)') 'Total multipoles referred to origin at'
          write(file_units(8),'(10x,3(a,f11.6),a9)') " x =", refpt(1)/ANGSTROM, &
               ',  y = ', refpt(2)/ANGSTROM, ',  z = ', refpt(3)/ANGSTROM, &
               ' angstrom'

          if(min_l <= 0 .and. max_l >= 0) then
             write(file_units(8),'(19x, a, f11.6)') 'Q00  =', -elec_and_core_q_total
          end if
          if(min_l <=1 .and. max_l >= 1) then
             Q1m1 = elec_and_core_Q1m1_refpt_total
             Q10 = elec_and_core_Q10_refpt_total
             Q11 = elec_and_core_Q11_refpt_total
             magnitude = multipole_1_magnitude_sph_real(Q10,Q1m1,Q11)

             if(magnitude >=1D3) then
                write(file_units(8),'(a6,1p,e11.3,a8,e11.3,a8,e11.3,a8,e11.3)') &
                '|Q1| =', magnitude, &
                '  Q10  =', Q10, '  Q11c =', Q11, '  Q11s =', Q1m1
             else
                write(file_units(8),'(a6,f11.6,a8,f11.6,a8,f11.6,a8,f11.6)') &
                '|Q1| =', magnitude, &
                '  Q10  =', Q10, '  Q11c =', Q11, '  Q11s =', Q1m1
             end if
          end if
          if(min_l <= 2 .and. max_l >= 2) then
             Q2m2 = elec_and_core_Q2m2_refpt_total
             Q2m1 = elec_and_core_Q2m1_refpt_total
             Q20 = elec_and_core_Q20_refpt_total
             Q21 = elec_and_core_Q21_refpt_total
             Q22 = elec_and_core_Q22_refpt_total
             magnitude = multipole_2_magnitude_sph_real(Q20,Q2m1,Q21,Q2m2,Q22)

             if(magnitude >=1D3) then
                write(file_units(8),'(a6,1p,e11.3,a8,e11.3,a8,e11.3,a8,e11.3)') &
                '|Q2| =', magnitude, &
                '  Q20  =', Q20, '  Q21c =', Q21, '  Q21s =', Q2m1
                write(file_units(8),'(17x,a8,1p,e11.3,a8,e11.3)') &
                '  Q22c =', Q22, '  Q22s =', Q2m2
             else
                write(file_units(8),'(a6,f11.6,a8,f11.6,a8,f11.6,a8,f11.6)') &
                '|Q2| =', multipole_2_magnitude_sph_real(Q20,Q2m1,Q21,Q2m2,Q22), &
                '  Q20  =', Q20, '  Q21c =', Q21, '  Q21s =', Q2m1
                write(file_units(8),'(17x,a8,f11.6,a8,f11.6)') &
                '  Q22c =', Q22, '  Q22s =', Q2m2
             end if
          end if

       end if

    end if ! last atom

    ! jd: Skip next bit of output, but ensure cleanup happens
    if(skip_this_atom) goto 300

    if(pub_print_qc) then
       write(qc_tag,'(i0)') orig_idx
       if(min_l <= 0 .and. max_l >= 0) then
          call utils_qc_print('DMA_mp_'//trim(qc_tag)//'_0', -elec_q, &
               minimum_magnitude = 0.001_DP)
       end if
       if(min_l <= 1 .and. max_l >= 1) then
          call utils_qc_print('DMA_mp_'//trim(qc_tag)//'_1_x', elec_mu_X_refpt,&
               minimum_magnitude = 0.001_DP)
          call utils_qc_print('DMA_mp_'//trim(qc_tag)//'_1_y', elec_mu_Y_refpt,&
               minimum_magnitude = 0.001_DP)
          call utils_qc_print('DMA_mp_'//trim(qc_tag)//'_1_z', elec_mu_Z_refpt,&
               minimum_magnitude = 0.001_DP)
       end if
       if(min_l <= 2 .and. max_l >= 2) then
          call utils_qc_print('DMA_mp_'//trim(qc_tag)//'_2_xx', elec_Txx_refpt,&
               minimum_magnitude = 0.001_DP)
          call utils_qc_print('DMA_mp_'//trim(qc_tag)//'_2_xy', elec_Txy_refpt,&
               minimum_magnitude = 0.001_DP)
          call utils_qc_print('DMA_mp_'//trim(qc_tag)//'_2_xz', elec_Txz_refpt,&
               minimum_magnitude = 0.001_DP)
          call utils_qc_print('DMA_mp_'//trim(qc_tag)//'_2_yy', elec_Tyy_refpt,&
               minimum_magnitude = 0.001_DP)
          call utils_qc_print('DMA_mp_'//trim(qc_tag)//'_2_yz', elec_Tyz_refpt,&
               minimum_magnitude = 0.001_DP)
          call utils_qc_print('DMA_mp_'//trim(qc_tag)//'_2_zz', elec_Tzz_refpt,&
               minimum_magnitude = 0.001_DP)
       end if
    end if

300 continue

  end subroutine print_mp_elms_for_atom

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module dma
