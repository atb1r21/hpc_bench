!================================================================!
!                                                                !
!                 Spherical wave expansion module                !
!                                                                !
! This module contains subroutines for expanding products of     !
! NGWFs in an auxiliary basis set of truncated spherical waves.  !
!                                                                !
! Only things that relate to the actual expansion of NGWFs live  !
! here.                                                          !
!                                                                !
! Everything that pertains to the metric matrices, generation of !
! SWs and SWpots, and associated caches should go into           !
! sw_resolution_of_identity_mod.                                 !
!                                                                !
! Since June 2017 also pair-densities (NGWF pairs with a DKN     !
! element sandwiched in) can be fitted. The resultant expansion  !
! coefficient are then spin-dependent. Density fitting (as       !
! opposed to NGWF-pair fitting) is currently only used in an     !
! exotic variant of polarisable embedding (where vacuum density  !
! and polarisation density are split).                           !
!                                                                !
!----------------------------------------------------------------!
! - The original module (as hf_exchange_mod) was written by      !
!   Quintin Hill in 2008/9 with supervision by Chris-Kriton      !
!   Skylaris, then extensively modified by Jacek Dziedzic in     !
!   2012 and 2013.                                               !
! - sw_expansion_mod was then extracted in April 2014 from       !
!   hf_exchange_mod by Jacek Dziedzic, in preparation for DMA.   !
! - This module was then pruned by moving some of the things to  !
!   sw_resolution_of_identity_mod in February 2015 by            !
!   Jacek Dziedzic.                                              !
! - Extended to density fitting (as opposed to ngwf_pair fitting)!
!   in June 2017 by Jacek Dziedzic.                              !
! - Extended to handle mixed NGWF sets (for hybrid conduction    !
!   and in preparation for hybrid TDDFT) in July/August 2018 by  !
!   Jacek Dziedzic.                                              !
!================================================================!

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

#define PPDS_DIMS ngwf_basis%max_n_ppds_sphere*cell%n_pts

module sw_expansion

  use constants, only: DP, CRLF, stdout
  use utils, only: utils_trace_in, utils_trace_out

  implicit none

  private

  ! --------------------------------
  ! --- Major public subroutines ---
  ! --------------------------------

  public :: swx_init_module
  public :: swx_cleanup_module
  public :: swx_expand_local_ngwf_pairs
  public :: swx_expand_my_ngwf_pairs
  public :: swx_expand_local_pair_densities
  public :: swx_communicate_coeffs

  ! --------------------------------
  ! --- Minor public subroutines ---
  ! --------------------------------

  public :: swx_expansion_2
  public :: swx_expansion_2_fast
  public :: swx_rms_fitting_error
  public :: swx_swop_quantity_overlap_ppd
  public :: swx_swop_quantity_overlap_fftbox_dbl

  ! ----------------------------------
  ! --- Public types and variables ---
  ! ----------------------------------

  ! jd: --- Public state ---
  public :: swx_calculate_rms_fitting_error

  ! ----------------------------------------
  ! --- Module-wide state and parameters ---
  ! ----------------------------------------

  logical, save :: swx_calculate_rms_fitting_error ! Should we calc fitting error?
  integer, parameter :: VACUUM_TAG = 1

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !---------------------------------------------------------------------------!
  ! ****               P U B L I C          R O U T I N E S              **** !
  !---------------------------------------------------------------------------!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_init_module(rep)
    !==========================================================================!
    ! This subroutine initialises the swx module.                              !
    !  - DMA swexes are initialised for properties and for pol_emb.            !
    !  - Some trivial module-wide constants are initialised.                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2015.                              !
    ! Extended by Jacek Dziedzic in May 2017.                                  !
    !==========================================================================!

    use constants, only: REP_SWEX_PROPERTIES_DMA_1, REP_SWEX_POL_EMB_DMA_1
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_devel_code, pub_dma_max_q, pub_dma_metric, &
         pub_dma_bessel_averaging, pub_dma_use_ri, pub_dma_max_l, &
         pub_pol_emb_qmstar, pub_pol_emb_dma_min_l, pub_pol_emb_dma_max_l
    use sw_expansion_type, only: swx_init_container
    use sw_resolution_of_identity, only: swri_get_handle_to, swri_library
    use utils, only: utils_devel_code, utils_flushed_string_output, &
         utils_int_to_str, utils_assert

    implicit none

    ! jd: Arguments
    type(NGWF_REP), intent(inout) :: rep

    ! jd: Local variables
    integer :: max_q
    integer :: dma_run, n_dma_runs
    integer :: swex_h
    integer :: dma_swri_h
    character(len=*), parameter :: myself = 'swx_init_module'

    ! -------------------------------------------------------------------------

    call utils_flushed_string_output(&
         'SWX: Initialising module... '//CRLF)

    ! jd: Initialise module-wide constants
    swx_calculate_rms_fitting_error = utils_devel_code(.false., 'SWX', &
         'CALCULATE_RMS_FITTING_ERROR', pub_devel_code)

    ! jd: Initialise DMA SW_EX or SW_EXes if needed
    if(pub_dma_max_l > -1) then
       dma_swri_h = swri_get_handle_to(pub_dma_use_ri)
       n_dma_runs = merge(2,1,pub_dma_bessel_averaging)
       swex_h = REP_SWEX_PROPERTIES_DMA_1
       max_q = pub_dma_max_q

       do dma_run = 1, n_dma_runs
          call utils_assert(.not. rep%swexes(swex_h)%initialised, &
               myself//': Attempt to initialise a SW_EX that has already &
               &been initialised (properties DMA)')
          call swx_init_container(rep%swexes(swex_h), &
               swri_library(dma_swri_h), &
               'properties_dma_swex-'//trim(utils_int_to_str(dma_run)), &
               0, pub_dma_max_l, max_q, pub_dma_metric, rep%postfix)
          max_q = max_q - 1
          swex_h = swex_h + 1
       end do

    end if

    ! jd: Initialise POLEMB SW_EX or SW_EXes if needed
    if(pub_pol_emb_qmstar) then
       dma_swri_h = swri_get_handle_to(pub_dma_use_ri)
       n_dma_runs = merge(2,1,pub_dma_bessel_averaging)
       swex_h = REP_SWEX_POL_EMB_DMA_1
       max_q = pub_dma_max_q

       do dma_run = 1, n_dma_runs
          call utils_assert(.not. rep%swexes(swex_h)%initialised, &
               myself//': Attempt to initialise a SW_EX that has already &
               &been initialised (POLEMB DMA)')
          call swx_init_container(rep%swexes(swex_h), &
               swri_library(dma_swri_h), &
               'polemb_dma_swex-'//trim(utils_int_to_str(dma_run)), &
               pub_pol_emb_dma_min_l, pub_pol_emb_dma_max_l, &
               max_q, pub_dma_metric, rep%postfix)
          max_q = max_q - 1
          swex_h = swex_h + 1
       end do

    end if

    call utils_flushed_string_output('done.'//CRLF)

  end subroutine swx_init_module

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_cleanup_module
    !==========================================================================!
    ! This subroutine cleans up after the swx module.                          !
    ! Currently nothing happens here, as SW_EX containers are freed via        !
    ! ngwf_rep_destroy().                                                      !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2015.                              !
    !==========================================================================!

    use utils, only: utils_flushed_string_output

    implicit none

    ! -------------------------------------------------------------------------

    call utils_flushed_string_output('SWX: Clean-up... ')

    call utils_flushed_string_output('done.'//CRLF)

  end subroutine swx_cleanup_module

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_expand_local_ngwf_pairs(rep, ireg, swri, &           ! in/out
      swex_h, ngwf_basis, cell, fftbox, par, elements, coeffs_kind, & ! input
      nl, &                                                           ! in-opt
      rep2, ngwf_basis2, ket_basis_selector)  ! opts for mixed NGWF bases
    !==========================================================================!
    ! This subroutine performs a spherical wave expansion of all overlapping   !
    ! NGWF pairs, storing the resulting expansion coefficients in coeffs_ht of !
    ! a SW_EX container specified via swex_h, of the NGWF_REP rep.             !
    !                                                                          !
    ! The pairs to expand are rank-local-B's with their s-neighbours C's.      !
    ! C's are not necessarily rank-local.                                      !
    !                                                                          !
    ! Typically, the coefficients are then communicated between procs according!
    ! to the neighbour list nl (which would be the X-neighbour list in HFx or  !
    ! the S-neighbour list in DMA). If 'nl' is omitted, the coefficients are   !
    ! *not* communicated (each proc then owns the coefficients for NGWF        !
    ! products of its atoms and their S-neighbours, modulo some complications  !
    ! when the SWRI does not span all atoms).                                  !
    !                                                                          !
    ! Two levels of caching are in effect. The coefficients themselves are     !
    ! cached, and actucally returned, in coeffs_ht. These are automatically    !
    ! invalidated when necessary. NGWFs are cached in via the remote module,   !
    ! and rep%ngwf_cache_handle selects the appropriate cache (so that we      !
    ! do not confuse different NGWF bases). These caches are also invalidated  !
    ! when necessary (=> ngwf_rep_register_change()). Also, ngwf_pair_fitting()!
    ! that is called from here caches SWOPs behind the scenes, when it calls   !
    ! expansion_2(), which calls swop_quantity_overlap_ppd(), which calls      !
    ! swri_obain_swops_in_ppd(), where the actual caching happens (in SWRI).   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   rep (in/out):  The NGWF_REP to expand. Results of the expansion also   !
    !                  go here.                                                !
    !   swri (in/out): The SWRI in which we expand. In/out because it caches   !
    !                  SWOPs behind the scenes.                                !
    !   swex_h (in):   Handle to the SW_EX container in rep. The SW_EX itself  !
    !                  is not passed, because it would alias with rep.         !
    !                  Overridden by ket_basis_selector(3), if it is passed.   !
    !   ngwf_basis, cell, fftbox, elements: The usual.                         !
    !   coeffs_kind (in): Specifies the expansion kind from which the coeffs   !
    !                     originate (SW_V, SW_O, ...).                         !
    !   nl (in, opt): Neighbour list for communicating coeffs, if they are     !
    !                 to be communicated.                                      !
    !                                                                          !
    ! *Optional* arguments for mixed NGWF bases only.                          !
    !   rep2 (in): Second rep to use when mixing NGWF bases. Can be intent-in, !
    !              because results (coeffs) are always stored in the first rep.!
    !   ngwf_basis2 (in): Second NGWF basis.                                   !
    !   ket_basis_selector (in): Selects the NGWF basis (1 or 2) in which      !
    !                            - Bb is indexed (ket_basis_selector(1)),      !
    !                            - Cc is indexed (ket_basis_selector(2)).      !
    !                            ket_basis_selector(3) supplies an override    !
    !                            for swex_h (allows choosing the destination   !
    !                            SW_EX in a rep for the coeffs).               !
    !--------------------------------------------------------------------------!
    ! Note:                                                                    !
    !   - Regardless of what is passed for 'nl', if anything, the routine      !
    !     swx_density_fitting called from here uses swri%s_atoms_nl            !
    !     to find out S-neighbours of our atoms. This is intended.             !
    !--------------------------------------------------------------------------!
    ! Cornerstones written by Quintin Hill in 2008 and January 2009.           !
    ! Modified by Quintin Hill on 12/02/2009 to make hfexchange symmetric.     !
    ! Rewritten by Jacek Dziedzic in July 2012 to use 2-centre 'BC' expansion. !
    ! Heavily modified by Jacek Dziedzic in April 2014 to separate from HFx.   !
    ! Generalised for SW_RI/SW_EX in February 2015 by Jacek Dziedzic.          !
    ! Extended for mixed NGWF basis sets in June-August 2018 by Jacek Dziedzic.!
    ! Modified for embedding by Joseph Prentice, September 2018                !
    !==========================================================================!

    use comms, only: pub_total_num_procs
    use constants, only: SW_LETTERS, NORMAL
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use neighbour_list, only: NL_NEIGHBOUR_LIST
    use ngwf_representation, only: NGWF_REP
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_swx_output_detail
    use timer, only: timer_clock
    use simulation_cell, only: CELL_INFO
    use sw_resolution_of_identity, only: SW_RI
    use utils, only: utils_assert, utils_abort, utils_flushed_string_output, &
         utils_postfix_to_ngwf_set_name

    implicit none

    ! jd: Arguments
    type(NGWF_REP), intent(inout)                   :: rep
    type(SW_RI), intent(inout)                      :: swri
    integer, intent(in)                             :: swex_h
    type(FUNC_BASIS), intent(in), target            :: ngwf_basis
    type(FFTBOX_INFO), intent(in)                   :: fftbox
    type(CELL_INFO), intent(in)                     :: cell
    type(PARAL_INFO), intent(in)                    :: par
    type(ELEMENT), intent(in)                       :: elements(par%nat)
    integer, intent(in)                             :: coeffs_kind
    integer, intent(in)                             :: ireg
    type(NL_NEIGHBOUR_LIST), intent(in), optional   :: nl
    ! jd: Optional arguments for mixed NGWF bases
    type(NGWF_REP), optional, intent(in)            :: rep2
    type(FUNC_BASIS), optional, intent(in), target  :: ngwf_basis2
    integer, optional, intent(in)                   :: ket_basis_selector(3)

    ! jd: Local variables
    integer                     :: dest_swex_h
    type(FUNC_BASIS), pointer   :: beta_ngwf_basis
    type(FUNC_BASIS), pointer   :: gamma_ngwf_basis
    character(len=*), parameter :: myself = 'swx_expand_local_ngwf_pairs'

    ! -----------------------------------------------------------------------

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    !jmecmplx
    call utils_assert(.not. rep%ngwfs_on_grid(ireg)%iscmplx, 'Error in '//myself//&
         ': not ready for complex NGWFs yet.')

    ! jd: Ensure consistency between optional parameters for mixed NGWF bases
    call utils_assert((present(rep2) .eqv. present(ngwf_basis2)) .and. &
         (present(rep2) .eqv. present(ket_basis_selector)), &
         myself//': Optional arguments for SWx with mixed NGWF bases must &
         &either all be present or all be absent.')

    ! jd: Select the right NGWF basis and SW_EX
    if(present(ket_basis_selector)) then
       dest_swex_h = ket_basis_selector(3)
       if(ket_basis_selector(1) == 1) then
          beta_ngwf_basis => ngwf_basis
       else if(ket_basis_selector(1) == 2) then
          beta_ngwf_basis => ngwf_basis2
       else
          call utils_abort(myself//': Illegal ket_basis_selector(1)')
       end if
       if(ket_basis_selector(2) == 1) then
          gamma_ngwf_basis => ngwf_basis
       else if(ket_basis_selector(2) == 2) then
          gamma_ngwf_basis => ngwf_basis2
       else
          call utils_abort(myself//': Illegal ket_basis_selector(2)')
       end if
    else
       dest_swex_h = swex_h
       beta_ngwf_basis => ngwf_basis
       gamma_ngwf_basis => ngwf_basis
    end if

    ! jd: The majority of times that we get here is because NGWFs changed and
    !     need to be re-expanded. Other scenarios, where NGWFs have not changed
    !     since last call are: a) Second run of Bessel averaging, where we need
    !     a second expansion for the same NGWFs, b) Mixing HFx and DMA with
    !     different coeff_kind's.
    if(pub_swx_output_detail >= NORMAL) then
       call utils_flushed_string_output('SWX: '//&
            trim(utils_postfix_to_ngwf_set_name(rep%postfix))//'%['//&
            trim(rep%swexes(dest_swex_h)%swex_name)//']: Fitting local NGWF pairs &
            &('//SW_LETTERS(coeffs_kind)//').'//CRLF)
    end if

    call utils_assert(rep%swexes(dest_swex_h)%initialised, &
         myself//': SW_EX container not initialised')

    ! jd: Perform pair fitting on all products of NGWFs local to this core
    !     and their neighbours (not necessarily local).
    call swx_local_ngwf_pair_fitting(rep, &                             ! in/out
         ireg, swri, swex_h, ngwf_basis, cell, fftbox, par, elements, &
         coeffs_kind, rep2, ngwf_basis2, ket_basis_selector)

    ! jd: Engage in the exchange of the obtained expansion coefficients, they
    !     will be needed on procs that are, in general, different from the
    !     ones where the expansion took place.
    if(pub_total_num_procs > 1 .and. present(nl)) then
       call swx_communicate_coeffs(rep%swexes(dest_swex_h)%coeffs_ht, &
            ngwf_basis, nl, swri%s_atoms_nl, swri, par, elements, coeffs_kind, &
            dest_swex_h, 'P', ngwf_basis2, ket_basis_selector)
    end if

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine swx_expand_local_ngwf_pairs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_expand_my_ngwf_pairs(rep, swri, &                    ! in/out
      packed_metric_matrix_blocks_hts, my_bc_pairs, n_my_bc_pairs, &  ! in
      swex_h, ngwf_basis, cell, par, elements, coeffs_kind, &         ! input
      rep2, ngwf_basis2, ket_basis_selector)  ! opts for mixed NGWF bases
    !==========================================================================!
    ! This subroutine performs a spherical wave expansion of all overlapping   !
    ! NGWF pairs, storing the resulting expansion coefficients in coeffs_ht of !
    ! a SW_EX container specified via swex_h, of the NGWF_REP rep.             !
    !                                                                          !
    ! The pairs to expand are specified via arguments. Neither B atoms nor C's !
    ! are in general rank-local.                                               !
    !                                                                          !
    ! This subroutine is only called from HFx. DMA and polemb rely on the      !
    ! usual ONETEP parallelisation scheme and use swx_expand_local_ngwf_pairs()!
    ! instead.                                                                 !
    !                                                                          !
    ! Unlike in swx_expand_local_ngwf_pairs(), there are NO comms of           !
    ! coefficients taking place.                                               !
    !                                                                          !
    ! Two levels of caching are in effect. The coefficients themselves are     !
    ! cached, and actually returned, in rep%coeffs_ht. These are automatically !
    ! invalidated when necessary. NGWFs are cached via the remote module,      !
    ! and rep%ngwf_cache_handle selects the appropriate cache (so that we      !
    ! do not confuse different NGWF bases). These caches are also invalidated  !
    ! when necessary (=> ngwf_rep_register_change()).                          !
    !                                                                          !
    ! Finallly, swx_my_ngwf_pair_fitting() that is called from here relies on  !
    ! SWOP caches having been populated (in hfx_populate_swop_cache()) when it !
    ! calls swx_expansion_2_fast(), which                                      !
    ! calls swx_swop_quantity_overlap_ppd_fast(), which                        !
    ! calls swri_obain_swops_in_ppd_ptr_fast(), which actually references the  !
    ! SWOP cache. Elements absent in the cache are recalculated on the fly, but!
    ! not added (since the cache is full already).                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   rep (in/out):  The NGWF_REP to expand. Results of the expansion also   !
    !                  go here.                                                !
    !   swri (in/out): The SWRI in which we expand. In/out because it caches   !
    !                  SWOPs behind the scenes.                                !
    !   packed_metric_matrix_blocks_hts (in): Provides necessary metric matrix !
    !                                         blocks (since B-C pairs are not  !
    !                                         rank-local).                     !
    !   my_bc_pairs (in): List of BC pairs for which expansions will be calced.!
    !   n_my_bc_pairs (in): Number of elements in the above.                   !
    !   swex_h (in):   Handle to the SW_EX container in rep. The SW_EX itself  !
    !                  is not passed, because it would alias with rep.         !
    !                  Overridden by ket_basis_selector(3), if it is passed.   !
    !   ngwf_basis, cell, par, elements: The usual.                            !
    !   coeffs_kind (in): Specifies the expansion kind from which the coeffs   !
    !                     originate (SW_V, SW_O, ...).                         !
    !                                                                          !
    ! *Optional* arguments for mixed NGWF bases only.                          !
    !   rep2 (in): Second rep to use when mixing NGWF bases. Can be intent-in, !
    !              because results (coeffs) are always stored in the first rep.!
    !   ngwf_basis2 (in): Second NGWF basis.                                   !
    !   ket_basis_selector (in): Selects the NGWF basis (1 or 2) in which      !
    !                            - Bb is indexed (ket_basis_selector(1)),      !
    !                            - Cc is indexed (ket_basis_selector(2)).      !
    !                            ket_basis_selector(3) supplies an override    !
    !                            for swex_h (allows choosing the destination   !
    !                            SW_EX in a rep for the coeffs).               !
    !--------------------------------------------------------------------------!
    ! Cannibalised from swx_expand_local_ngwf_pairs() by Jacek Dziedzic in     !
    ! April 2019.                                                              !
    !==========================================================================!

    use constants, only: SW_LETTERS, NORMAL
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE
    use ion, only: ELEMENT
    use ngwf_representation, only: NGWF_REP
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_swx_output_detail
    use timer, only: timer_clock
    use simulation_cell, only: CELL_INFO
    use sw_resolution_of_identity, only: SW_RI
    use utils, only: utils_assert, utils_abort, utils_flushed_string_output, &
         utils_postfix_to_ngwf_set_name

    implicit none

    ! jd: Arguments
    type(NGWF_REP), intent(inout), target   :: rep
    type(HT_HASH_TABLE), intent(in), target :: packed_metric_matrix_blocks_hts(:)
    integer, intent(in)                     :: my_bc_pairs(:,:)
    integer, intent(in)                     :: n_my_bc_pairs
    type(SW_RI), intent(inout)              :: swri
    integer, intent(in)                     :: swex_h
    type(FUNC_BASIS), intent(in), target    :: ngwf_basis
    type(CELL_INFO), intent(in)             :: cell
    type(PARAL_INFO), intent(in)            :: par
    type(ELEMENT), intent(in)               :: elements(par%nat)
    integer, intent(in)                     :: coeffs_kind
    ! jd: Optional arguments for mixed NGWF bases
    type(NGWF_REP), optional, intent(in), target    :: rep2
    type(FUNC_BASIS), optional, intent(in), target  :: ngwf_basis2
    integer, optional, intent(in)                   :: ket_basis_selector(3)

    ! jd: Local variables
    integer                     :: dest_swex_h
    type(FUNC_BASIS), pointer   :: beta_ngwf_basis
    type(FUNC_BASIS), pointer   :: gamma_ngwf_basis
    type(NGWF_REP), pointer     :: beta_rep
    type(NGWF_REP), pointer     :: gamma_rep
    character(len=*), parameter :: myself = 'swx_expand_my_ngwf_pairs'

    ! -----------------------------------------------------------------------

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    ! jd: Ensure consistency between optional parameters for mixed NGWF bases
    call utils_assert((present(rep2) .eqv. present(ngwf_basis2)) .and. &
         (present(rep2) .eqv. present(ket_basis_selector)), &
         myself//': Optional arguments for SWx with mixed NGWF bases must &
         &either all be present or all be absent.')

    ! jd: Select the right NGWF basis and SW_EX
    if(present(ket_basis_selector)) then
       dest_swex_h = ket_basis_selector(3)
       if(ket_basis_selector(1) == 1) then
          beta_rep => rep
          beta_ngwf_basis => ngwf_basis
       else if(ket_basis_selector(1) == 2) then
          beta_rep => rep2
          beta_ngwf_basis => ngwf_basis2
       else
          call utils_abort(myself//': Illegal ket_basis_selector(1)')
       end if
       if(ket_basis_selector(2) == 1) then
          gamma_rep => rep
          gamma_ngwf_basis => ngwf_basis
       else if(ket_basis_selector(2) == 2) then
          gamma_rep => rep2
          gamma_ngwf_basis => ngwf_basis2
       else
          call utils_abort(myself//': Illegal ket_basis_selector(2)')
       end if
    else
       dest_swex_h = swex_h
       beta_rep => rep
       gamma_rep => rep
       beta_ngwf_basis => ngwf_basis
       gamma_ngwf_basis => ngwf_basis
    end if

    ! jd: The majority of times that we get here is because NGWFs changed and
    !     need to be re-expanded. Other scenarios, where NGWFs have not changed
    !     since last call are: a) Second run of Bessel averaging, where we need
    !     a second expansion for the same NGWFs, b) Mixing HFx and DMA with
    !     different coeff_kind's.
    if(pub_swx_output_detail >= NORMAL) then
       call utils_flushed_string_output('SWX: '//&
            trim(utils_postfix_to_ngwf_set_name(rep%postfix))//'%['//&
            trim(rep%swexes(dest_swex_h)%swex_name)//']: Fitting my NGWF pairs &
            &('//SW_LETTERS(coeffs_kind)//').'//CRLF)
    end if

    call utils_assert(rep%swexes(dest_swex_h)%initialised, &
         myself//': SW_EX container not initialised')

    ! jd: Perform density fitting on all products of NGWFs specified in the
    !     arguments. Neither B's nor C's are in general rank-local.
    call swx_my_ngwf_pair_fitting(beta_rep, gamma_rep, &           ! in/out, in
         my_bc_pairs, n_my_bc_pairs, packed_metric_matrix_blocks_hts, & ! in
         swri, swex_h, cell, par, elements, coeffs_kind, &              ! in
         beta_ngwf_basis, gamma_ngwf_basis)                             ! in

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine swx_expand_my_ngwf_pairs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_expand_local_pair_densities(coeffs_ht, rep, ireg, swri, &
      swex_h, denskern, ngwf_basis, cell, fftbox, elements, coeffs_kind, &
      nl, vacuum_rep, vacuum_denskern)
    !==========================================================================!
    ! This subroutine performs a spherical wave expansion of all atom-pair     !
    ! densities, i.e. of all \phi_Aa K^Aa,Bb \phi_Bb where A and B have overla-!
    ! pping NGWFs. This differs from swx_expand_nwgf_pairs() in the following: !
    ! - the coefficients are DKN-dependent;                                    !
    ! - the coefficients are stored in a separate coeffs_ht (first argument),  !
    !   and not in rep%swexes(:)%coeffs_ht, which means rep can be intent(in). !
    !                                                                          !
    ! The pairs to expand are rank-local-B's with their s-neighbours C's.      !
    ! C's are not necessarily rank-local.                                      !
    !                                                                          !
    ! Typically, the coefficients are then communicated between procs according!
    ! to the neighbour list nl (which would be the X-neighbour list in HFx or  !
    ! the S-neighbour list in DMA). If 'nl' is omitted, the coefficients are   !
    ! *not* communicated (each proc then owns the coefficients for NGWF        !
    ! products of its atoms and their S-neighbours, modulo some complications  !
    ! when the SWRI does not span all atoms).                                  !
    !                                                                          !
    ! Two levels of caching are in effect. The coefficients themselves are     !
    ! cached in coeffs_ht. Since these are disassociated from any rep, it is   !
    ! up to the caller to track invalidations. The second level of caching is  !
    ! via remote_ngwf_cache_hts in remote_mod. These  are invalidated when     !
    ! necessary, via ngwf_rep_register_change(). Finally, swx_density_fitting()!
    ! that is called from here caches SWOPs behind the scenes, when it calls   !
    ! expansion_2(), which calls swop_quantity_overlap_ppd(), which calls      !
    ! swri_obain_swops_in_ppd(), where the actual caching happens (in SWRI).   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   coeffs_ht (in/out): The expansion coefficients are stored here. It is  !
    !                       the responsibility of the caller to pass a suitably!
    !                       initialised hash table.                            !
    !   rep (in):      The NGWF_REP whose NGWFs to expand.                     !
    !   swri (in/out): The SWRI in which we expand. In/out because it caches   !
    !                  SWOPs behind the scenes.                                !
    !   swex_h (in):   Handle to the SW_EX container in rep. The SW_EX itself  !
    !                  is not passed, because it would alias with rep.         !
    !   denskern (in): The density kernel.                                     !
    !   ngwf_basis, cell, fftbox, elements: The usual.                         !
    !   coeffs_kind (in): Specifies the expansion kind from which the coeffs   !
    !                     originate.                                           !
    !   nl (in, opt): Neighbour list for communicating coeffs.                 !
    !--------------------------------------------------------------------------!
    ! Expert mode:                                                             !
    !   When arguments 'vacuum_rep' and 'vacuum denskern' are specified, the   !
    !   expanded quantity will be the density difference, ie.                  !
    !   \phi_Aa K^Aa,Bb \phi_Bb - \phi0_Aa K0^Aa,Bb \phi_Bb, where \phi0 refer !
    !   to vacuum_rep, and K0 refers to vacuum_denskern.                       !
    !--------------------------------------------------------------------------!
    ! Note:                                                                    !
    !   - Regardless of what is passed for 'nl', if anything, the routine      !
    !     swx_density_fitting called from here uses swri%s_atoms_nl            !
    !     to find out S-neighbours of our atoms. This is intended.             !
    !--------------------------------------------------------------------------!
    ! Caveat:                                                                  !
    !   This subroutine is *not* ready for mixed NGWF bases. I don't envision  !
    !   hybrid polarisable embedding conduction with vacuum reference states   !
    !   to ever be needed. If it turns out I'm wrong, modify this subroutine   !
    !   in the same fashion swx_expand_ngwf_pairs() was modified in June-August!
    !   2018.                                                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in June 2017 using swx_expand_ngwfs() as a     !
    ! template.                                                                !
    ! Modified for embedding by Joseph Prentice, September 2018                !
    !==========================================================================!

    use comms, only: pub_total_num_procs
    use constants, only: SW_LETTERS, NORMAL
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE
    use ion, only: ELEMENT
    use neighbour_list, only: NL_NEIGHBOUR_LIST
    use ngwf_representation, only: NGWF_REP
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_num_spins, pub_swx_output_detail
    use sparse, only: SPAM3
    use timer, only: timer_clock
    use simulation_cell, only: CELL_INFO
    use sparse, only: sparse_get_par
    use sw_resolution_of_identity, only: SW_RI
    use utils, only: utils_assert, utils_flushed_string_output

    implicit none

    ! jd: Arguments
    type(HT_HASH_TABLE), intent(inout), target    :: coeffs_ht
    type(NGWF_REP), intent(in)                    :: rep
    type(SW_RI), intent(inout)                    :: swri
    integer, intent(in)                           :: swex_h
    type(SPAM3), intent(in)                       :: denskern(pub_num_spins)
    type(FUNC_BASIS), intent(in)                  :: ngwf_basis
    type(FFTBOX_INFO), intent(in)                 :: fftbox
    type(CELL_INFO), intent(in)                   :: cell
    type(ELEMENT), intent(in)                     :: elements(:)
    integer, intent(in)                           :: coeffs_kind
    integer, intent(in)                           :: ireg
    type(NL_NEIGHBOUR_LIST), intent(in), optional :: nl
    type(NGWF_REP), intent(inout), optional       :: vacuum_rep
    type(SPAM3), intent(in), optional             :: vacuum_denskern(pub_num_spins)

    ! jd: Local variables
    character(len=*), parameter :: myself = 'swx_expand_local_pair_densities'
    type(PARAL_INFO), pointer   :: par

    ! -----------------------------------------------------------------------

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    !jmecmplx
    call utils_assert(.not. rep%ngwfs_on_grid(ireg)%iscmplx, 'Error in '//myself//&
         ': not ready for complex NGWFs yet.')

    call utils_assert(present(vacuum_rep) .eqv. present(vacuum_denskern), &
         myself//'Optional arguments must all be specified or all be absent')

    ! rc2013: get the parallel strategy from the diagonal matrix
    call sparse_get_par(par, denskern(1))
    call utils_assert(par%nat == size(elements), 'Error in '//myself//': &
         &allocated parallel strategy is incompatible with elements.')

    ! jd: The majority of times that we get here is because NGWFs changed and
    !     need to be re-expanded. Other scenarios, where NGWFs have not changed
    !     since last call are: a) Second run of Bessel averaging, where we need
    !     a second expansion for the same NGWFs, b) Mixing HFx and DMA with
    !     different coeff_kind's.
    if(pub_swx_output_detail >= NORMAL) then
       call utils_flushed_string_output('SWX: ['//&
            trim(rep%swexes(swex_h)%swex_name)//']: Performing density fitting &
            &('//SW_LETTERS(coeffs_kind)//').'//CRLF)
    end if

    ! jd: Perform density fitting on all products of DKN*NGWFs local to this
    !     proc, and their neighbours (not necessarily local).
    call swx_local_density_fitting(coeffs_ht, rep, ireg, swri, swex_h, &
         denskern, ngwf_basis, cell, fftbox, elements, coeffs_kind, &
         vacuum_rep, vacuum_denskern)

    ! jd: Engage in the exchange of the obtained expansion coefficients, they
    !     will be needed on procs that are, in general, different from the
    !     ones where the expansion took place.
    if(pub_total_num_procs > 1 .and. present(nl)) then
       call swx_communicate_coeffs(coeffs_ht, ngwf_basis, &
            nl, swri%s_atoms_nl, swri, par, elements, coeffs_kind, swex_h, 'D')
    end if

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine swx_expand_local_pair_densities

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_expansion_2(swri, &
       swex_quality, swex_h, is, &
       expansion_centres, expansion_atoms, &
       num_sws_in_expansion, coeffs_kind, cell, vmatrix, wmatrix, omatrix, &
       global_bb_ngwf_idx, global_cc_ngwf_idx, &
       coeffs_ht, q_coeffs, &
       ppd_quantity, n_ppds_in_quantity, ppd_indices_in_quantity, &
       fftbox_quantity, fftbox, fftbox_start_in_cell_vec)
    !==========================================================================!
    ! This subroutine expands a quantity in a PPD representation in terms of   !
    ! SWs/SWpots originating on up to 2 centres.                               !
    ! The 'quantity' is typically either a) an NGWF product (in DMA, mode 'P'),!
    ! b) an atom-pair density \phi_Aa K^Aa,Bb \phi_Bb) (in DMA, mode 'D'), or  !
    ! c) something more involved (Q_AaDd in HFx gradients). As of 2020.04      !
    ! (c) is not used anymore -- HFx gradients use swx_expansion_2_fast().     !
    !                                                                          !
    ! The results (expansion coeffs) are returned in coeffs_ht in cases a) and !
    ! b), or in q_coeffs (not a hash table!) in case c).                       !
    !                                                                          !
    ! We do not store the coefficients in a SWEX, since swx_expansion_2() has  !
    ! to be callable from hfx_gradient_term_atomblock(), where 'rep', and      !
    ! hence its swexes, are intent(in). It is still OK to pass                 !
    ! rep%swexes(idx)%coeffs_ht as coeffs_ht if desired.                       !
    !                                                                          !
    ! In cases a) and b) global_bb_ngwf_idx and global_cc_ngwf_idx need to be  !
    ! provided, because there is opportunity to cache SWOPs in PPDs.           !
    !                                                                          !
    ! In case c) ppd_quantity is a sum of weighted products originating on     !
    ! multiple spheres (Q_AaDd). Since there are no meaningful NGWF indices    !
    ! under which this could be stored in coeffs_ht, the coeffs are returned   !
    ! in q_coeffs.                                                             !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   swri (inout): The SW_RI underlying the expansion.                      !
    !                 Inout because the SWOP cache is updated.                 !
    !   swex_quality (in): Describes the quality of the expansion.             !
    !   swex_h (in): A user-defined index to be stored in the fourth key of    !
    !                the hash table. Typically a swex handle like REP_SWEX_HFX !
    !                would be passed here. This is used to distinguish coeffs  !
    !                from different SWEXes, e.g. in Bessel averaging.          !
    !   is (in): Spin index which will be used as the fifth key of the hash    !
    !            table.                                                        !
    !   expansion_centres (in): } Describe the SW centres where the SWs        !
    !   expansion_atoms (in):   } originate.                                   !
    !   num_sws_in_expansion (in): Number of SWs in the expansion, yes.        !
    !   coeffs_kind (in): See below.                                           !
    !   cell (in): The usual.                                                  !
    !   vmatrix, wmatrix, omatrix (in): Relevant metric matrix blocks.         !
    !   global_bb_ngwf_idx, global_cc_ngwf_idx (in): Identify the stored coeffs!
    !   coeffs_ht (inout, opt): The results go here in cases a) and b).        !
    !   q_coeffs (out, opt): The results go here in case c).                   !
    !                                                                          !
    !   For calculations on the coarse grid (swx_dbl_grid F) the quantity to   !
    !   be expanded is in PPDs, and the following arguments are relevant:      !
    !                                                                          !
    !   ppd_quantity (in):            } Describes the quantity                 !
    !   n_ppds_in_quantity (in):      } to be expanded, in terms               !
    !   ppd_indices_in_quantity (in): } of a set of PPDs.                      !
    !                                                                          !
    !   Expert mode:                                                           !
    !   ------------                                                           !
    !                                                                          !
    !   For calculations on the double grid (swx_dbl_grid T) the quantity is   !
    !   in a double-grid FFT box, and the following arguments are relevant:    !
    !                                                                          !
    !   fftbox_quantity (in):          } Describes the quantity                !
    !   fftbox (in):                   } to be expanded, in terms              !
    !   fftbox_start_in_cell_vec (in): } of a double-grid FFT box.             !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! coeffs_kind specifies what kind of expansion coefficients are desired:   !
    ! - SW_V - expansion of quantity's potential in SWpots, via                !
    !   electrostatic metric,                                                  !
    ! - SW_W - expansion of quantity's potential in SWs, via                   !
    !   overlap metric ('OVO'),                                                !
    ! - SW_O - expansion of quantity in SWs (overlap metric).                  !
    !--------------------------------------------------------------------------!
    ! - Written by Jacek Dziedzic in February 2012, basing on a similar        !
    !   subroutine by Quintin Hill.                                            !
    ! - Extended in August 2012 by Jacek Dziedzic to handle the second case.   !
    ! - Simplified in February 2014 by Jacek Dziedzic to elide caching of      !
    !   swop-ngwf-product overlaps (SNPOs), as these are now never re-used.    !
    ! - Generalised to SW_RI/SW_EX by Jacek Dziedzic in February 2015.         !
    ! - Reworked for SW_RI/SW_EX responsibility reshuffle by JD in Oct 2016.   !
    ! - Generalised to support double grid by Jacek Dziedzic in January 2017.  !
    ! - Generalised to support density expansions (case b) by Jacek Dziedzic   !
    !   in June 2017.                                                          !
    !==========================================================================!

    use constants, only: SW_O, SW_V, SW_W, SWEX_BB_CC_SYMMETRIC
    use datatypes, only: FFTBOX_DATA
    use fft_box, only: FFTBOX_INFO
    use geometry, only: POINT
    use hash_table, only: HT_HASH_TABLE, hash_table_add, hash_table_list
    use rundat, only: pub_swri_overlap_indirect, pub_use_swx, pub_swx_dbl_grid, &
         pub_pol_emb_dbl_grid, pub_use_activeswx
    use linalg, only: linalg_dgesv
    use simulation_cell, only: CELL_INFO
    use sw_resolution_of_identity, only: ATOM_CENTRE, SW_RI, SW_QUALITY
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(inout)                   :: swri
    type(SW_QUALITY), intent(in)                 :: swex_quality
    integer, intent(in)                          :: swex_h
    integer, intent(in)                          :: is
    integer, intent(in)                          :: num_sws_in_expansion
    type(ATOM_CENTRE), intent(in)                :: expansion_centres(2)
    integer, intent(in)                          :: expansion_atoms(2)
    integer, intent(in)                          :: coeffs_kind
    type(CELL_INFO), intent(in)                  :: cell
    real(kind=DP), intent(in), allocatable       :: vmatrix(:,:)
    real(kind=DP), intent(in), allocatable       :: wmatrix(:,:)
    real(kind=DP), intent(in), allocatable       :: omatrix(:,:)
    integer, intent(in), optional                :: global_bb_ngwf_idx
    integer, intent(in), optional                :: global_cc_ngwf_idx
    type(HT_HASH_TABLE), intent(inout), optional, target :: coeffs_ht
    real(kind=DP), intent(out), optional         :: q_coeffs(:)
    real(kind=DP), intent(in), optional          :: ppd_quantity(:)
    integer, intent(in), optional                :: n_ppds_in_quantity
    integer, intent(in), optional                :: ppd_indices_in_quantity(:)

    ! jd: Arguments for expert mode
    type(FFTBOX_DATA), intent(in), optional      :: fftbox_quantity
    type(FFTBOX_INFO), intent(in), optional      :: fftbox
    type(POINT), intent(in), optional            :: fftbox_start_in_cell_vec

    ! Local variables
    real(kind=DP)              :: swcoeff(num_sws_in_expansion)
    integer                    :: ierr ! error flag
    real(kind=DP), allocatable :: swop_quantity_overlap(:)
    real(kind=DP)              :: accum
    integer                    :: row, k
    integer                    :: swex_max_sws_per_centre
    character                  :: sw_or_swpot
    character(len=*), parameter :: myself = 'swx_expansion_2'

    ! -----------------------------------------------------------------------

    call utils_trace_in(myself)

    call utils_assert(pub_use_swx .or. pub_use_activeswx, &
         trim(myself)//' should only be used with pub_use_swx or pub_use_activeswx.')

    if(coeffs_kind == SW_V) then
       call utils_assert(allocated(vmatrix) .and. .not. allocated(wmatrix) &
            .and. .not. allocated(omatrix), trim(myself)//&
            ': SW_V requires V matrix and no W matrix or O matrix.', &
            allocated(vmatrix),allocated(wmatrix),allocated(omatrix))
    end if

    if(coeffs_kind == SW_W) then
       call utils_assert(allocated(vmatrix) .and. allocated(wmatrix) &
            .and. allocated(omatrix), trim(myself)//&
            ': SW_W requires V matrix and W matrix and O matrix.', &
            allocated(vmatrix),allocated(wmatrix),allocated(omatrix))
    end if

    if(coeffs_kind == SW_O) then
       call utils_assert(.not. allocated(vmatrix) .and. &
            .not. allocated(wmatrix) &
            .and. allocated(omatrix), trim(myself)//&
            ': SW_O requires no V matrix or W matrix but O matrix.', &
            allocated(vmatrix),allocated(wmatrix),allocated(omatrix))
    end if

    ! jd: All three of Bb/Cc/coeffs_ht optional parameters should be supplied
    !     or omitted.
    call utils_assert(&
         (present(global_bb_ngwf_idx) .and. present(global_cc_ngwf_idx) .and. &
         present(coeffs_ht)) .eqv. &
         (present(global_bb_ngwf_idx) .or. present(global_cc_ngwf_idx) .or. &
         present(coeffs_ht)), &
         trim(myself)//': Inconsistency in optional Bb/Cc/coeffs_ht parameters')

    ! jd: Consistency between Bb/Cc/coeffs_ht and q_coeffs optional parameters
    call utils_assert(&
         present(global_bb_ngwf_idx) .neqv. present(q_coeffs), trim(myself)//&
         ': If Bb/Cc opt parameters are specified, q_coeffs cannot be, and vice&
         & versa')

    if(present(global_bb_ngwf_idx)) then
       if(SWEX_BB_CC_SYMMETRIC(swex_h)) then
          call utils_assert(global_bb_ngwf_idx <= global_cc_ngwf_idx, &
               trim(myself)//': For Bb-Cc symmetric SW_EX''s, NGWF indices must&
               & be sorted on entry.')
       end if
    end if

    ! jd: All three of PPD-related optional parameters should be supplied
    !     or omitted.
    call utils_assert((present(ppd_quantity) .and. &
         present(n_ppds_in_quantity) .and. &
         present(ppd_indices_in_quantity)) .eqv. &
         present(ppd_quantity) .or. present(n_ppds_in_quantity) .or. &
         present(ppd_indices_in_quantity), &
         myself//': Inconsistency in PPD-related parameters')

    ! jd: All three of FFT box-related optional parameters should be supplied
    !     or omitted.
    call utils_assert((present(fftbox_quantity) .and. present(fftbox) .and. &
         present(fftbox_start_in_cell_vec)) .eqv. (present(fftbox_quantity) .or. &
         present(fftbox) .or. present(fftbox_start_in_cell_vec)), &
         myself//': Inconsistency in FFT box-related parameters')

    call utils_assert(present(ppd_quantity) .neqv. present(fftbox_quantity), &
         myself//': Quantity to be expanded must be provided *either* in PPDs &
         &or in an FFT box')

    call utils_assert(present(fftbox_quantity) .eqv. pub_swx_dbl_grid, &
         myself//': FFT box quantity must only be provided iff swx_dbl_grid T')

    swex_max_sws_per_centre = swex_quality%max_q * &
         (1+swex_quality%max_l-swex_quality%min_l) * &
         (1+swex_quality%max_l+swex_quality%min_l)

    allocate(swop_quantity_overlap(num_sws_in_expansion),stat=ierr)
    call utils_alloc_check(myself,'swop_quantity_overlap',ierr)

    swcoeff = 0.0_DP

    if(coeffs_kind == SW_O .or. coeffs_kind == SW_W) then
       sw_or_swpot = 'S'
    else
       sw_or_swpot = 'P'
    end if

    ! jd: Calculate or retrieve from cache the overlaps of SWOPs with
    !     the quantity in PPDs.
    if(pub_swx_dbl_grid) then
       call swx_swop_quantity_overlap_fftbox_dbl(swri, &
            swex_quality, &
            fftbox_start_in_cell_vec, fftbox_quantity, fftbox, &
            swop_quantity_overlap, expansion_centres, expansion_atoms, &
            num_sws_in_expansion, sw_or_swpot, cell)
    else
       call swx_swop_quantity_overlap_ppd(swri, &
            swex_quality, &
            ppd_quantity, n_ppds_in_quantity, ppd_indices_in_quantity, &
            swop_quantity_overlap, expansion_centres, expansion_atoms, &
            num_sws_in_expansion, sw_or_swpot, cell)
    end if

    if(pub_swx_dbl_grid .or. pub_pol_emb_dbl_grid) then
       call hash_table_list(swri%swops_at_points_ht,0)
    end if

    if(coeffs_kind == SW_V) then
       ! jd: If the electrostatic metric is used we solve Vc = b for c,
       !     where V = vmatrix, c = swcoeff, b = swop_quantity_overlap
       call linalg_dgesv(&
            swcoeff(1:num_sws_in_expansion), &
            vmatrix, swop_quantity_overlap, num_sws_in_expansion, 1)
    end if

    if(coeffs_kind == SW_O) then
       ! jd: If the overlap coefficients are needed we solve Oc = b for c,
       !     where O = omatrix, c = swcoeff, b = swop_quantity_overlap
       call linalg_dgesv(&
            swcoeff(1:num_sws_in_expansion), &
            omatrix, swop_quantity_overlap, num_sws_in_expansion, 1)
    end if

    if(coeffs_kind == SW_W) then
       ! jd: If the overlap metric is used
       !     - for overlap_indirect F we calculate Wb = c
       !       where W = wmatrix (O^-1.V.O^-1), c = swcoeff,
       !       b = swop_quantity_overlap
       !     - for overlap_indirect T we solve Wc = b for c, just like
       !       with the electrostatic metric, except W = O.V^-1.O
       if(pub_swri_overlap_indirect) then
          call linalg_dgesv(&
               swcoeff(1:num_sws_in_expansion), &
               wmatrix, swop_quantity_overlap, num_sws_in_expansion, 1)
       else
          do row = 1, num_sws_in_expansion
             accum = 0.0_DP
             do k = 1, num_sws_in_expansion
                accum = accum + wmatrix(k,row) * swop_quantity_overlap(k)
                !@optimize me: mat*vec mul
             end do
             swcoeff(row) = accum
          end do
       end if
    end if

    if(present(global_bb_ngwf_idx)) then
       ! jd: Cases a) and b)
       call hash_table_add(coeffs_ht, swcoeff, num_sws_in_expansion, &
            global_bb_ngwf_idx, global_cc_ngwf_idx, coeffs_kind, swex_h, is, &
            overfill_strategy = 'F')
    else
       ! jd: Case c)
       q_coeffs(:) = 0.0_DP ! jd: With 1-centre expansion q_coeffs has a different
                            !     dimension than swcoeff (it's bigger).
       q_coeffs(1:num_sws_in_expansion) = swcoeff(1:num_sws_in_expansion)
    end if

    deallocate(swop_quantity_overlap,stat=ierr)
    call utils_dealloc_check(myself,'swop_quantity_overlap',ierr)

    call utils_trace_out(myself)

  end subroutine swx_expansion_2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_expansion_2_fast(coeffs, swri, &
       expansion_centres, expansion_atoms, &
       num_sws_in_expansion, coeffs_kind, cell, vmatrix, wmatrix, omatrix, &
       ppd_quantity, n_ppds_in_quantity, ppd_indices_in_quantity)
    !==========================================================================!
    ! This subroutine expands a quantity in a PPD representation in terms of   !
    ! SWs/SWpots originating on up to 2 centres. It is currently only used in  !
    ! HFx, but not DMA or polemb.                                              !
    !                                                                          !
    ! The 'quantity' is either an NGWF product or Q_AaDd in HFx NGWF gradient, !
    ! but this subroutine does not care. It just expands 'things in PPDs'.     !
    !                                                                          !
    ! The results (expansion coeffs) are returned in coeffs.                   !
    !                                                                          !
    ! In contrast to the non-fast version, it does not modify any HTs, so it's !
    ! safe to use in OMP with no CRITICALs needed. It relies on the caching of !
    ! SWOPs, but it does not add to the cache. Another difference is that here !
    ! we do not support swri-to-swex filtering (for efficiency). Finally, this !
    ! routine works on the coarse grid only.                                   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   coeffs(out): The first num_sws_in_expansion elements are populated with!
    !                the expansion coefficients, rest is set to garbage_real.  !
    !   swri (in): The SW_RI underlying the expansion.                         !
    !   expansion_centres (in): } Describe the SW centres where the SWs        !
    !   expansion_atoms (in):   } originate.                                   !
    !   num_sws_in_expansion (in): Number of SWs in the expansion, yes.        !
    !   coeffs_kind (in): See below.                                           !
    !   cell (in): The usual.                                                  !
    !   vmatrix, wmatrix, omatrix (in): Relevant metric matrix blocks.         !
    !   ppd_quantity (in):            } Describes the quantity                 !
    !   n_ppds_in_quantity (in):      } to be expanded, in terms               !
    !   ppd_indices_in_quantity (in): } of a set of PPDs.                      !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! coeffs_kind specifies what kind of expansion coefficients are desired:   !
    ! - SW_V - expansion of quantity's potential in SWpots, via                !
    !   electrostatic metric,                                                  !
    ! - SW_W - expansion of quantity's potential in SWs, via                   !
    !   overlap metric ('OVO'),                                                !
    ! - SW_O - expansion of quantity in SWs (overlap metric).                  !
    !--------------------------------------------------------------------------!
    ! - Written by Jacek Dziedzic in February 2012, basing on a similar        !
    !   subroutine by Quintin Hill.                                            !
    ! - Extended in August 2012 by Jacek Dziedzic to handle the second case.   !
    ! - Simplified in February 2014 by Jacek Dziedzic to elide caching of      !
    !   swop-ngwf-product overlaps (SNPOs), as these are now never re-used.    !
    ! - Generalised to SW_RI/SW_EX by Jacek Dziedzic in February 2015.         !
    ! - Reworked for SW_RI/SW_EX responsibility reshuffle by JD in Oct 2016.   !
    ! - Generalised to support double grid by Jacek Dziedzic in January 2017.  !
    ! - Generalised to support density expansions (case b) by Jacek Dziedzic   !
    !   in June 2017.                                                          !
    ! - Converted to 'fast' version by Jacek Dziedzic in August 2019.          !
    !==========================================================================!

    use constants, only: SW_O, SW_V, SW_W, garbage_real
    use rundat, only: pub_swri_overlap_indirect, pub_use_swx, pub_use_activeswx
    use linalg, only: linalg_dgesv
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use sw_resolution_of_identity, only: ATOM_CENTRE, SW_RI
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)                   :: coeffs(:)
    type(SW_RI), intent(in)                      :: swri
    integer, intent(in)                          :: num_sws_in_expansion
    type(ATOM_CENTRE), intent(in)                :: expansion_centres(2)
    integer, intent(in)                          :: expansion_atoms(2)
    integer, intent(in)                          :: coeffs_kind
    type(CELL_INFO), intent(in)                  :: cell
    real(kind=DP), intent(in), allocatable       :: vmatrix(:,:)
    real(kind=DP), intent(in), allocatable       :: wmatrix(:,:)
    real(kind=DP), intent(in), allocatable       :: omatrix(:,:)
    real(kind=DP), intent(in)                    :: ppd_quantity(:)
    integer, intent(in)                          :: n_ppds_in_quantity
    integer, intent(in)                          :: ppd_indices_in_quantity(:)

    ! Local variables
    real(kind=DP)              :: swcoeff(num_sws_in_expansion)
    integer                    :: ierr ! error flag
    real(kind=DP), allocatable :: swop_quantity_overlap(:)
    real(kind=DP)              :: accum
    integer                    :: row, k
    character                  :: sw_or_swpot
    character(len=*), parameter :: myself = 'swx_expansion_2_fast'

    ! -----------------------------------------------------------------------

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    call utils_assert(pub_use_swx .or. pub_use_activeswx, trim(myself)//&
         ' should only be used with pub_use_swx or pub_use_activeswx.')

    call utils_assert(coeffs_kind == SW_V .or. coeffs_kind == SW_O .or. &
         coeffs_kind == SW_W, myself//': Illegal coeffs_kind')

    if(coeffs_kind == SW_V) then
       call utils_assert(allocated(vmatrix) .and. .not. allocated(wmatrix) &
            .and. .not. allocated(omatrix), trim(myself)//&
            ': SW_V requires V matrix and no W matrix or O matrix.', &
            allocated(vmatrix),allocated(wmatrix),allocated(omatrix))
    end if

    if(coeffs_kind == SW_W) then
       call utils_assert(allocated(vmatrix) .and. allocated(wmatrix) &
            .and. allocated(omatrix), trim(myself)//&
            ': SW_W requires V matrix and W matrix and O matrix.', &
            allocated(vmatrix),allocated(wmatrix),allocated(omatrix))
    end if

    if(coeffs_kind == SW_O) then
       call utils_assert(.not. allocated(vmatrix) .and. &
            .not. allocated(wmatrix) &
            .and. allocated(omatrix), trim(myself)//&
            ': SW_O requires no V matrix or W matrix but O matrix.', &
            allocated(vmatrix),allocated(wmatrix),allocated(omatrix))
    end if

    allocate(swop_quantity_overlap(num_sws_in_expansion),stat=ierr)
    call utils_alloc_check(myself,'swop_quantity_overlap',ierr)

    swcoeff = 0D0

    if(coeffs_kind == SW_O .or. coeffs_kind == SW_W) then
       sw_or_swpot = 'S'
    else
       sw_or_swpot = 'P'
    end if

    ! jd: Calculate the overlaps of SWOPs with the quantity in PPDs.
    !     Uses OMP internally.
    call swx_swop_quantity_overlap_ppd_fast(swri, &
         ppd_quantity, n_ppds_in_quantity, ppd_indices_in_quantity, &
         swop_quantity_overlap, expansion_centres, expansion_atoms, &
         num_sws_in_expansion, sw_or_swpot, cell)

    if(coeffs_kind == SW_V) then
       ! jd: If the electrostatic metric is used we solve Vc = b for c,
       !     where V = vmatrix, c = swcoeff, b = swop_quantity_overlap
       call linalg_dgesv(swcoeff(1:num_sws_in_expansion), &
            vmatrix, swop_quantity_overlap, num_sws_in_expansion, 1)
    end if

    if(coeffs_kind == SW_O) then
       ! jd: If the overlap coefficients are needed we solve Oc = b for c,
       !     where O = omatrix, c = swcoeff, b = swop_quantity_overlap
       call linalg_dgesv(swcoeff(1:num_sws_in_expansion), &
            omatrix, swop_quantity_overlap, num_sws_in_expansion, 1)
    end if

    if(coeffs_kind == SW_W) then
       ! --- obsoleted 'OVO' way to do overlap ---
       ! jd: If the overlap metric is used
       !     - for overlap_indirect F we calculate Wb = c
       !       where W = wmatrix (O^-1.V.O^-1), c = swcoeff,
       !       b = swop_quantity_overlap
       !     - for overlap_indirect T we solve Wc = b for c, just like
       !       with the electrostatic metric, except W = O.V^-1.O
       if(pub_swri_overlap_indirect) then
          call linalg_dgesv(swcoeff(1:num_sws_in_expansion), &
               wmatrix, swop_quantity_overlap, num_sws_in_expansion, 1)
       else
          do row = 1, num_sws_in_expansion
             accum = 0.0_DP
             do k = 1, num_sws_in_expansion
                accum = accum + wmatrix(k,row) * swop_quantity_overlap(k)
                !@optimize me: mat*vec mul
             end do
             swcoeff(row) = accum
          end do
       end if
    end if

    coeffs(1:num_sws_in_expansion) = swcoeff(1:num_sws_in_expansion)
    coeffs(num_sws_in_expansion+1:) = garbage_real

    deallocate(swop_quantity_overlap,stat=ierr)
    call utils_dealloc_check(myself,'swop_quantity_overlap',ierr)

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine swx_expansion_2_fast

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function swx_rms_fitting_error(swex, &                ! in
       swri, &                                                        ! in/out
       n_ppds_in_quantity, ppd_indices_in_quantity, ppd_quantity, &   ! in
       expansion_atoms, expansion_centres, global_bb_ngwf_idx, &      ! in
       global_cc_ngwf_idx, cell, coeffs_kind, swex_h, is, &           ! in
       q_coeffs)                                                      ! in, opt
    !==========================================================================!
    ! Calculates the RMS error in a SW fit of a quantity in PPDs.              !
    !                                                                          !
    ! There are three cases of interest -- a) the quantity is a simple product !
    ! of two NGWFs (Bb, Cc), b) the quantity is an atom-pair density, and      !
    ! c) a sum of weighted products originating on multiple spheres (Q_AaDd).  !
    ! The values of global_bb_ngwf_idx and global_cc_ngwf_idx need to be -1 in !
    ! this last case.                                                          !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in August 2012.                                !
    ! Generalised for SW_EX/SW_RI in February 2015.                            !
    ! Generalised for case b) by Jacek Dziedzic in June 2017.                  !
    ! Fixed case c), by Jacek Dziedzic in August 2018.                         !
    !==========================================================================!

    use constants, only: PI, SWEX_BB_CC_SYMMETRIC
    use hash_table, only: hash_table_lookup_nocount
    use rundat, only: pub_swx_c_threshold, pub_debug
    use simulation_cell, only: CELL_INFO
    use sw_expansion_type, only: SW_EX
    use sw_resolution_of_identity, only: SW_RI, ATOM_CENTRE, &
         swri_swop_calc_all_in_ppd
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(SW_EX), intent(in)              :: swex
    type(SW_RI), intent(in)              :: swri
    integer, intent(in)                  :: n_ppds_in_quantity
    integer, intent(in)                  :: ppd_indices_in_quantity(:)
    real(kind=DP), intent(in)            :: ppd_quantity(:)
    integer, intent(in)                  :: expansion_atoms(2)
    type(ATOM_CENTRE), intent(in)        :: expansion_centres(2)
    integer, intent(in)                  :: global_bb_ngwf_idx
    integer, intent(in)                  :: global_cc_ngwf_idx
    type(CELL_INFO), intent(in)          :: cell
    integer, intent(in)                  :: coeffs_kind
    integer, intent(in)                  :: swex_h
    integer, intent(in)                  :: is
    real(kind=DP), intent(in), optional  :: q_coeffs(2*swex%quality%max_sws_per_centre)

    ! jd: Local variables
    real(kind=DP) :: rms
    integer :: sw1_idx, sw2_idx
    integer :: src_atom
    type(ATOM_CENTRE):: src_centre
    integer :: num_of_cached_coeffs
    integer :: num_sws
    integer :: ec
    integer :: pt
    integer :: cur_ppd_n, cur_ppd, offs_to_cur_ppd
    real(kind=DP) :: coeff
    real(kind=DP) :: coeffs(2*swex%quality%max_sws_per_centre)
    real(kind=DP) :: ppd_sw_sum(n_ppds_in_quantity * cell%n_pts)
    real(kind=DP) :: swops_in_a_ppd(cell%n_pts * swex%quality%max_sws_per_centre)
    character(len=*), parameter :: myself = 'swx_rms_fitting_error'

    ! -------------------------------------------------------------------------

    if(present(q_coeffs)) then
       coeffs(:) = q_coeffs(:)
    else
       ! --- Default case ---
       ! jd: Look up freshly calculated expansion coefficients for this Bb Cc
       if(.not. SWEX_BB_CC_SYMMETRIC(swex_h)) then
          ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
          call hash_table_lookup_nocount(coeffs, num_of_cached_coeffs, &
               swex%coeffs_ht, &
               global_bb_ngwf_idx, global_cc_ngwf_idx, coeffs_kind, swex_h, is)
       else
          ! jd: Here there's Bb-Cc symmetry and we only store one triangle
          call hash_table_lookup_nocount(coeffs, num_of_cached_coeffs, &
               swex%coeffs_ht, &
               min(global_bb_ngwf_idx,global_cc_ngwf_idx), &
               max(global_bb_ngwf_idx,global_cc_ngwf_idx), coeffs_kind, swex_h, is)
       end if

       call utils_assert(num_of_cached_coeffs /= -1, &
            myself//': coeffs hash table too small or coeffs missing', &
            global_bb_ngwf_idx, global_cc_ngwf_idx, coeffs_kind, swex_h)
    end if

    ppd_sw_sum(:) = 0.0_DP

    ! --------------------------------------------------------------
    ! Go over relevant PPDs (those on in BbCc)                   PPD
    ! --------------------------------------------------------------
    all_ppds:                                                      &
    do cur_ppd_n = 1, n_ppds_in_quantity

       cur_ppd = ppd_indices_in_quantity(cur_ppd_n)
       offs_to_cur_ppd = (cur_ppd_n-1)*cell%n_pts+1

       ! jd: Loop over unique centres in the expansion
       sw2_idx = 1
       two_centres:                                          &
       do ec = 1, 2
          src_atom = expansion_atoms(ec)
          src_centre = expansion_centres(ec)
          if(src_atom == -1) exit

          ! - obtain the SWOPs for this PPD
          call swri_swop_calc_all_in_ppd(swops_in_a_ppd, num_sws, &   ! out
               src_centre, cur_ppd, 'S', cell, swri, swex%quality)    ! in

          if(pub_debug) call utils_assert(swex%quality%num_sws_per_centre == &
               num_sws, myself//': Inconsistent number of SWOPs')

          ! --------------------------------------------------------------------
          ! jd: for all SWs on src atom                                       SW
          ! --------------------------------------------------------------------
          sws_on_one:                                                          &
          do sw1_idx = 1, swex%quality%num_sws_per_centre

             ! elide when expansion coefficient is zero or very small
             ! zero coefficients happen e.g. when Bb and Cc spheres barely
             ! overlap and there are no points inside. Then |BbCc) is zero
             ! and all coeffs are zero.

             coeff = coeffs(sw2_idx)
             if(abs(coeff) > pub_swx_c_threshold) then
                ppd_sw_sum(&
                     offs_to_cur_ppd:offs_to_cur_ppd+cell%n_pts-1) = &
                     ppd_sw_sum(&
                     offs_to_cur_ppd:offs_to_cur_ppd+cell%n_pts-1) + &
                     coeff * swops_in_a_ppd(&
                     (sw1_idx-1)*cell%n_pts+1:sw1_idx*cell%n_pts)
             end if

             sw2_idx = sw2_idx + 1

          end do sws_on_one

       end do two_centres

    end do all_ppds

    ppd_sw_sum = ppd_sw_sum / (4.0_DP * PI)

    rms = 0.0_DP
    do pt = 1, n_ppds_in_quantity * cell%n_pts
       rms = rms + (ppd_quantity(pt)-ppd_sw_sum(pt))**2
    end do
    swx_rms_fitting_error = sqrt(rms / (n_ppds_in_quantity * cell%n_pts))

  end function swx_rms_fitting_error

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

  subroutine swx_local_ngwf_pair_fitting(rep, ireg, swri, &             ! in/out
       swex_h, ngwf_basis, cell, fftbox, par, elements, coeffs_kind, &  ! in
       rep2, ngwf_basis2, ket_basis_selector)
    !==========================================================================!
    ! Performs ngwf-product-fitting, that is, for all pairs of S-overlapping   !
    ! atoms all NGWF products are expanded in terms of SWs originating on      !
    ! those two atoms. The resulting expansion coefficients are stored in      !
    ! coeffs_ht of a SW_EX container in rep. The container is selected with a  !
    ! handle swex_h.                                                           !
    ! Caching of SWOPs in PPDs is performed behind the scenes. The metric used !
    ! for expansion is selected with coeffs_kind.                              !
    !                                                                          !
    ! The pairs to expand are rank-local-B's with their s-neighbours C's.      !
    ! C's are not necessarily rank-local.                                      !
    !                                                                          !
    ! If the underlying SW_RI does not encompass the entire system, only the   !
    ! NGWF products between in-region and in-region and in-region and their    !
    ! overlapping atoms are expanded (ie. all products except those between    !
    ! an outsider and an outsider).                                            !
    !                                                                          !
    ! The calculated coefficients are *not* communicated across procs, this is !
    ! done in swx_expand_local_ngwf_pairs().                                   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   rep (inout): The NGWF_REP for which to perform density fitting.        !
    !                Results are also stored here.                             !
    !   swri (inout): The SW_RI underlying the expansion. SWOPs are cached     !
    !                 behind the scenes.                                       !
    !   swex_h (in): Handle for choosing which SW_EX container in rep will     !
    !                be used to store results. The container itself is not     !
    !                passed to avoid aliasing with rep.                        !
    !   ngwf_basis, cell, elements (in): The usual.                            !
    !   fftbox (in): Only used when swx_double_grid T.                         !
    !   coeffs_kind (in): See swx_expansion_2().                               !
    !                                                                          !
    ! *Optional* arguments for mixed NGWF bases only.                          !
    !   rep2 (in): Second rep to use when mixing NGWF bases. Can be intent-in, !
    !              because results (coeffs) are always stored in the first rep.!
    !   ngwf_basis2 (in): Second NGWF basis.                                   !
    !   ket_basis_selector (in): Selects the NGWF basis (1 or 2) in which      !
    !                            - Bb is indexed (ket_basis_selector(1)),      !
    !                            - Cc is indexed (ket_basis_selector(2)).      !
    !                            ket_basis_selector(3) supplies an override    !
    !                            for swex_h (allows choosing the destination   !
    !                            SW_EX in a rep for the coeffs).               !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2012.                                  !
    ! Generalised for SW_EX/SW_RI in February 2015.                            !
    ! Cleaned up by Jacek Dziedzic in June 2017.                               !
    ! Extended for mixed NGWF basis sets in June-August 2018 by Jacek Dziedzic.!
    ! Modified for embedding by Joseph Prentice, September 2018                !
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_send, &
         comms_barrier, comms_wait
    use constants, only: SW_V, SW_W, SW_O, VERBOSE, SWEX_BB_CC_SYMMETRIC
    use datatypes, only: FFTBOX_DATA, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT
    use hash_table, only: hash_table_list, hash_table_probe
    use ion, only: ELEMENT
    use ngwf_representation, only: NGWF_REP
    use parallel_strategy, only: PARAL_INFO
    use ppd_ops, only: ppd_product
    use remote, only: packed_ngwf_size, NGWF_REQUEST_TAG, remote_obtain_ngwf, &
         remote_serve_ngwfs, remote_pack_ngwf, remote_unpack_ngwf_to_fftbox, &
         remote_ngwf_cache_hts
    use rundat, only: pub_devel_code, pub_swx_output_detail, pub_swx_dbl_grid, &
         pub_pol_emb_dbl_grid
    use simulation_cell, only: CELL_INFO
    use sw_resolution_of_identity, only: SW_RI, ATOM_CENTRE, swri_init_centre, &
         swri_extract_matrix_2, swri_expansion_centres, swri_make_wmatrix_2, &
         swri_get_handle_to
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_alloc_check, &
         utils_dealloc_check, utils_assert, utils_devel_code, utils_abort

    implicit none

    ! jd: Arguments
    type(NGWF_REP), intent(inout), target :: rep
    type(SW_RI), intent(inout)            :: swri
    integer, intent(in)                   :: swex_h
    type(FUNC_BASIS), intent(in), target  :: ngwf_basis
    type(CELL_INFO), intent(in)           :: cell
    type(FFTBOX_INFO), intent(in)         :: fftbox
    type(PARAL_INFO), intent(in)          :: par
    type(ELEMENT), intent(in)             :: elements(par%nat)
    integer, intent(in)                   :: coeffs_kind
    integer, intent(in)                   :: ireg
    ! jd: Arguments for mixed NGWF bases
    type(NGWF_REP), optional, intent(in), target    :: rep2
    type(FUNC_BASIS), optional, intent(in), target  :: ngwf_basis2
    integer, optional, intent(in)                   :: ket_basis_selector(3)

    ! jd: --- Local variables ---
    ! jd: Atom-counting and indices
    integer :: local_b
    integer :: global_b, global_c, c_idx, orig_global_b
    type(ATOM_CENTRE) :: centres(2)

    ! jd: NGWF-counting and indices
    integer :: local_bb_ngwf_idx
    integer :: global_bb_ngwf_idx, global_cc_ngwf_idx
    integer :: first_ngwf_idx_of_b, first_ngwf_idx_of_c
    integer :: local_idx_of_1st_ngwf_on_b, local_idx_of_1st_ngwf_on_c
    integer :: ngwf_b, ngwf_c

    ! jd: Centres and atoms in expansion
    integer :: num_sws_in_expansion
    integer :: expansion_atoms_out(2)
    integer :: expansion_atoms_in(2)
    type(ATOM_CENTRE) :: expansion_centres(2)

    ! jd: Matrices in expansion
    real(kind=DP), allocatable :: vmatrix(:,:)
    real(kind=DP), allocatable :: omatrix(:,:)
    real(kind=DP), allocatable :: wmatrix(:,:)

    ! jd: PPDs and PPD indices
    integer :: n_ppds_in_product_of_bb_cc
    integer :: ppd_indices_in_product_of_bb_cc(ngwf_basis%max_n_ppds_sphere)
    real(kind=DP) :: ppd_product_of_bb_cc(PPDS_DIMS)

    ! jd: NGWFs in FFT boxes
    type(FFTBOX_DATA) :: fftbox_ngwf_bb, fftbox_ngwf_cc, &
         product_of_bb_cc_in_b_fftbox
    type(POINT) :: fftbox_bb_start_in_cell_vec
    integer :: delta_mid_complement(3)

    ! jd: Packed NGWFs and their comms
    real(kind=DP)          :: packed_ngwf_bb(packed_ngwf_size) ! local
    real(kind=DP), pointer :: packed_ngwf_cc(:) ! remote
    logical :: done
    integer :: proc
    integer :: send_handles(0:pub_total_num_procs-1)
    logical :: who_is_done(0:pub_total_num_procs-1)

    ! jd: Varia
    integer :: probe_result
    integer :: ierr
    real(kind=DP) :: rms
    integer :: swri_handle
    logical :: debug_show_skipped

    ! jd: For mixed NGWF basis sets
    integer                   :: dest_swex_h
    type(NGWF_REP), pointer   :: beta_rep
    type(NGWF_REP), pointer   :: gamma_rep
    type(FUNC_BASIS), pointer :: beta_ngwf_basis
    type(FUNC_BASIS), pointer :: gamma_ngwf_basis

    character(len=*), parameter :: myself = 'swx_local_ngwf_pair_fitting'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)

    call comms_barrier
    call timer_clock(myself,1)

    ! jd: Ensure consistency between optional parameters for mixed NGWF bases
    call utils_assert((present(rep2) .eqv. present(ngwf_basis2)) .and. &
         (present(rep2) .eqv. present(ket_basis_selector)), &
         myself//': Optional arguments for SWx with mixed NGWF bases must &
         &either all be present or all be absent.')

    ! jd: Select the right NGWF basis and SW_EX
    if(present(ket_basis_selector)) then
       dest_swex_h = ket_basis_selector(3)
       if(ket_basis_selector(1) == 1) then
          beta_ngwf_basis => ngwf_basis
          beta_rep => rep
       else if(ket_basis_selector(1) == 2) then
          beta_ngwf_basis => ngwf_basis2
          beta_rep => rep2
       else
          call utils_abort(myself//': Illegal ket_basis_selector(1)')
       end if
       if(ket_basis_selector(2) == 1) then
          gamma_ngwf_basis => ngwf_basis
          gamma_rep => rep
       else if(ket_basis_selector(2) == 2) then
          gamma_ngwf_basis => ngwf_basis2
          gamma_rep => rep2
       else
          call utils_abort(myself//': Illegal ket_basis_selector(2)')
       end if
    else
       dest_swex_h = swex_h
       beta_ngwf_basis => ngwf_basis
       gamma_ngwf_basis => ngwf_basis
       beta_rep => rep
       gamma_rep => rep
    end if

    swri_handle = swri_get_handle_to(swri%swri_name)

    debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
         pub_devel_code, no_bcast=.true., no_warn=.true.)

    ! jd: Sanity check on neighbour lists
    call utils_assert(swri%s_atoms_nl%populated, &
         'Internal error in '//trim(myself)&
         //': s_atoms neighbour list not populated for '//trim(swri%swri_name))

    who_is_done(:) = .false.

    ! --------------------------------------------------------------------------
    ! jd: Loop over B's -- all atoms on this core                            BBB
    ! --------------------------------------------------------------------------
    loop_B:                                                                    &
    do local_b=1, par%num_atoms_on_proc(pub_my_proc_id)
       global_b = par%first_atom_on_proc(pub_my_proc_id) + local_b-1
       orig_global_b = par%orig_atom(global_b)

       ! jd: Ignore atoms not represented in this SWRI
       if(.not. elements(orig_global_b)%in_swri(swri_handle)) then
          if(debug_show_skipped) then
             write(stdout,'(a,i0,a,i0,a,a,a)') &
             'Skipping ('//myself//') orig/SFC: ', orig_global_b, '/', &
             global_b, ' (', trim(elements(orig_global_b)%species_id),')'
          end if
          cycle ! [*]
       end if

       if(pub_swx_output_detail >= VERBOSE) write(stdout,*) ' - B: ', global_b

       call swri_init_centre(centres(1), par, elements, global_b)

       ! jd: List hash table usage for SWOPS, coeffs, NGWFs.
       call hash_table_list(remote_ngwf_cache_hts(rep%ngwf_cache_handle),0)
       if(present(rep2)) then
          if(rep2%ngwf_cache_handle /= rep%ngwf_cache_handle) then
             call hash_table_list(remote_ngwf_cache_hts(rep2%ngwf_cache_handle),0)
          end if
       end if
       call hash_table_list(rep%swexes(dest_swex_h)%coeffs_ht,0)
       call hash_table_list(swri%swops_in_ppds_ht,0)
       if(pub_swx_dbl_grid .or. pub_pol_emb_dbl_grid) then
          call hash_table_list(swri%swops_at_points_ht,0)
       end if

       first_ngwf_idx_of_b = beta_ngwf_basis%first_on_atom(global_b)
       local_idx_of_1st_ngwf_on_b = first_ngwf_idx_of_b + 1 - &
           beta_ngwf_basis%first_on_proc(pub_my_proc_id)

       ! -----------------------------------------------------------------------
       ! jd: Loop over C's that are s-neighbours with B                      CCC
       ! -----------------------------------------------------------------------
       loop_C:                                                                 &
       do c_idx = swri%s_atoms_nl%first_idx(global_b), &
            swri%s_atoms_nl%last_idx(global_b)

          global_c = swri%s_atoms_nl%neighbours(c_idx)
          ! jd: No elimination of atoms outside SWRI here, we do want the
          !     input from BC if B is in SWRI and C is outside

          call swri_init_centre(centres(2), par, elements, global_c)

          first_ngwf_idx_of_c = gamma_ngwf_basis%first_on_atom(global_c)
          local_idx_of_1st_ngwf_on_c = first_ngwf_idx_of_c + 1 - &
              gamma_ngwf_basis%first_on_proc(pub_my_proc_id)

          ! jd: Work around 'temporary was created' debug mode warning
          expansion_atoms_in = (/global_b, global_c/)

          ! jd: Figure out the expansion centres, sort them, pad with -1
          !     if fewer than 2 are needed
          call swri_expansion_centres(swri, &
               expansion_centres, expansion_atoms_out, &  ! out
               num_sws_in_expansion, &                    ! out
               centres, expansion_atoms_in, &             ! in
               rep%swexes(dest_swex_h)%quality, &         ! in
               SWEX_BB_CC_SYMMETRIC(dest_swex_h))         ! in

          ! jd: Allocate the required 2-centre matrices (O, V, W).
          if(coeffs_kind == SW_V .or. coeffs_kind == SW_W) then
             allocate(vmatrix(num_sws_in_expansion,num_sws_in_expansion), &
                  stat=ierr)
             call utils_alloc_check(myself,'vmatrix',ierr)
          end if
          if(coeffs_kind == SW_W) then
             allocate(wmatrix(num_sws_in_expansion,num_sws_in_expansion), &
                  stat=ierr)
             call utils_alloc_check(myself,'wmatrix',ierr)
          end if
          if(coeffs_kind == SW_O .or. coeffs_kind == SW_W) then
             allocate(omatrix(num_sws_in_expansion,num_sws_in_expansion), &
                  stat=ierr)
             call utils_alloc_check(myself,'omatrix',ierr)
          end if

          ! jd: Prepare the 2-centre V matrix
          if(coeffs_kind == SW_V .or. coeffs_kind == SW_W) then
             call swri_extract_matrix_2(vmatrix, &                ! out
                  swri, num_sws_in_expansion, &                   ! in
                  expansion_centres, expansion_atoms_out, SW_V, & ! in
                  rep%swexes(dest_swex_h)%quality)                ! in
          end if

          ! jd: Prepare the 2-centre O matrix
          if(coeffs_kind == SW_O .or. coeffs_kind == SW_W) then
             call swri_extract_matrix_2(omatrix, &                ! out
                  swri, num_sws_in_expansion, &                   ! in
                  expansion_centres, expansion_atoms_out, SW_O, & ! in
                  rep%swexes(dest_swex_h)%quality)                ! in
          end if

          ! jd: Prepare the 2-centre W matrix from O and V
          if(coeffs_kind == SW_W) then
             call swri_make_wmatrix_2(wmatrix, vmatrix, omatrix, &
                  num_sws_in_expansion)
          end if

          ! --------------------------------------------------------------------
          ! jd: for all c on C                                               ccc
          ! --------------------------------------------------------------------
          loop_ngwf_c:                                                         &
          do ngwf_c = 1, gamma_ngwf_basis%num_on_atom(global_c)
             global_cc_ngwf_idx = first_ngwf_idx_of_c + ngwf_c - 1

             ! @optimizeme: check first if we perhaps have all coefficients
             ! for this c, ie for all b's. If so, no need to obtain and expand Cc
             ! since no-one does comms during swx_expansion_2, there are likely
             ! convoy effects here now.
             call remote_obtain_ngwf(packed_ngwf_cc, &
                  remote_ngwf_cache_hts, who_is_done, &
                  global_cc_ngwf_idx, gamma_rep%ngwfs_on_grid(ireg), &
                  gamma_ngwf_basis, gamma_rep%ngwf_cache_handle, &
                  overfill_strategy='F')
             if(pub_swx_dbl_grid) then
                ! jd: Unpack Cc NGWF to an FFT box midway between Aa and Bb
                call remote_unpack_ngwf_to_fftbox(fftbox_ngwf_cc, &
                     packed_ngwf_cc, fftbox, cell, &
                     centre_for_target_fftbox = &
                     beta_ngwf_basis%spheres(local_idx_of_1st_ngwf_on_b)%centre, &
                     delta_mid_complement = delta_mid_complement)
             end if

             ! -----------------------------------------------------------------
             ! jd: for all b on B                                            bbb
             ! -----------------------------------------------------------------
             loop_ngwf_b:                                                         &
             do ngwf_b = 1, beta_ngwf_basis%num_on_atom(global_b)
                global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1
                local_bb_ngwf_idx = global_bb_ngwf_idx - &
                      beta_ngwf_basis%first_on_proc(pub_my_proc_id)+1

                ! jd: In a single NGWF basis scenario, expansions for (Bb,Cc)
                !     and (Cc,Bb) are identical (symmetry under label exchange).
                !     But we cannot just expand for Bb <= Cc, when some Bs can
                !     be skipped due to not being in the SWRI region. We need to
                !     handle overlaps of in-region Bb's with out-of-region Cc's
                !     regardless of ordering. We thus explicitly check if the
                !     coeff is there, and expand if not.

                ! jd: Check if we have coefficients for this Bb Cc

                if(.not. SWEX_BB_CC_SYMMETRIC(dest_swex_h)) then
                   ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                   probe_result = &
                        hash_table_probe(rep%swexes(dest_swex_h)%coeffs_ht, & ! in
                        global_bb_ngwf_idx, global_cc_ngwf_idx, &             ! in
                        coeffs_kind, dest_swex_h, 1)
                else
                   ! jd: Here there's Bb-Cc symmetry and we only store Bb<=Cc triangle
                   probe_result = &
                        hash_table_probe(rep%swexes(dest_swex_h)%coeffs_ht, & ! in
                        min(global_bb_ngwf_idx,global_cc_ngwf_idx), &         ! in
                        max(global_bb_ngwf_idx,global_cc_ngwf_idx), &         ! in
                        coeffs_kind, dest_swex_h, 1)
                end if
                if(probe_result == -1) then

                   ! *** If absent, expand ***

                   ! jd: Pack NGWF Bb for ppd_product
                   call remote_pack_ngwf(packed_ngwf_bb, &
                        beta_rep%ngwfs_on_grid(ireg), beta_ngwf_basis, &
                        local_bb_ngwf_idx)

                   ! ==========================================
                   ! ============ EXPANSION PROPER ============
                   ! ==========================================
                   if(pub_swx_dbl_grid) then
                      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      ! ~~~ double grid version ~~~
                      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      ! jd: Compute product of NGWFs Bb and Cc -- double grid
                      call remote_unpack_ngwf_to_fftbox(fftbox_ngwf_bb, &
                           packed_ngwf_bb, fftbox, cell, &
                           fftbox_start_in_cell_vec = fftbox_bb_start_in_cell_vec, &
                           gridpoint_offset = delta_mid_complement)
                      call swx_ngwf_product(product_of_bb_cc_in_b_fftbox, & !out
                           fftbox_ngwf_bb, fftbox_ngwf_cc, fftbox)
                      call data_fftbox_dealloc(fftbox_ngwf_bb)
                      ! jd: Expand product in SWOPs -- double grid

                      if(.not. SWEX_BB_CC_SYMMETRIC(dest_swex_h)) then
                         ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                         call swx_expansion_2(swri, &                       ! in/out
                              rep%swexes(dest_swex_h)%quality, dest_swex_h, 1, & ! in
                              expansion_centres, expansion_atoms_out, & ! in
                              num_sws_in_expansion, coeffs_kind, &      ! in
                              cell, vmatrix, wmatrix, omatrix, &        ! in
                              global_bb_ngwf_idx, global_cc_ngwf_idx, & !} for caching
                              rep%swexes(dest_swex_h)%coeffs_ht, &      ! in/out
                              fftbox_quantity = product_of_bb_cc_in_b_fftbox, &
                              fftbox = fftbox, &
                              fftbox_start_in_cell_vec = fftbox_bb_start_in_cell_vec)
                      else
                         ! jd: Here there's Bb-Cc symmetry and we only store one triangle
                         call swx_expansion_2(swri, &                       ! in/out
                              rep%swexes(dest_swex_h)%quality, dest_swex_h, 1, &      ! in
                              expansion_centres, expansion_atoms_out, &     ! in
                              num_sws_in_expansion, coeffs_kind, &          ! in
                              cell, vmatrix, wmatrix, omatrix, &            ! in
                              min(global_bb_ngwf_idx,global_cc_ngwf_idx), & !} for
                              max(global_bb_ngwf_idx,global_cc_ngwf_idx), & !} caching
                              rep%swexes(dest_swex_h)%coeffs_ht, &     ! in/out
                              fftbox_quantity = product_of_bb_cc_in_b_fftbox, &
                              fftbox = fftbox, &
                              fftbox_start_in_cell_vec = fftbox_bb_start_in_cell_vec)
                      end if
                      call data_fftbox_dealloc(product_of_bb_cc_in_b_fftbox)
                   else
                      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      ! ~~~ coarse grid (usual) version ~~~
                      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      ! jd: Compute product of NGWFs Bb and Cc -- coarse grid
                      call ppd_product(packed_ngwf_bb, packed_ngwf_cc, cell%n_pts, & ! in
                           ppd_product_of_bb_cc, &                         ! out
                           ppd_indices_in_product_of_bb_cc, &              ! out
                           n_ppds_in_product_of_bb_cc)                     ! out
                      ! jd: Expand product in SWOPs -- coarse grid
                      if(.not. SWEX_BB_CC_SYMMETRIC(dest_swex_h)) then
                         ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                         call swx_expansion_2(swri, &                       ! in/out
                              rep%swexes(dest_swex_h)%quality, dest_swex_h, 1, &      ! in
                              expansion_centres, expansion_atoms_out, & ! in
                              num_sws_in_expansion, coeffs_kind, &      ! in
                              cell, vmatrix, wmatrix, omatrix, &        ! in
                              global_bb_ngwf_idx, global_cc_ngwf_idx, & !} for caching
                              rep%swexes(dest_swex_h)%coeffs_ht, &      ! in/out
                              ppd_quantity = ppd_product_of_bb_cc, &    ! in
                              n_ppds_in_quantity = n_ppds_in_product_of_bb_cc, & ! in
                              ppd_indices_in_quantity = ppd_indices_in_product_of_bb_cc) ! in
                      else
                         ! jd: Here there's Bb-Cc symmetry and we only store one triangle
                         call swx_expansion_2(swri, &                       ! in/out
                              rep%swexes(dest_swex_h)%quality, dest_swex_h, 1, &      ! in
                              expansion_centres, expansion_atoms_out, &     ! in
                              num_sws_in_expansion, coeffs_kind, &          ! in
                              cell, vmatrix, wmatrix, omatrix, &            ! in
                              min(global_bb_ngwf_idx,global_cc_ngwf_idx), & !} for
                              max(global_bb_ngwf_idx,global_cc_ngwf_idx), & !} caching
                              rep%swexes(dest_swex_h)%coeffs_ht, &     ! in/out
                              ppd_quantity = ppd_product_of_bb_cc, &        ! in
                              n_ppds_in_quantity = n_ppds_in_product_of_bb_cc, & ! in
                              ppd_indices_in_quantity = ppd_indices_in_product_of_bb_cc) ! in
                      end if

                      if(swx_calculate_rms_fitting_error) then
                         rms = swx_rms_fitting_error(rep%swexes(dest_swex_h), swri, &
                              n_ppds_in_product_of_bb_cc, &
                              ppd_indices_in_product_of_bb_cc, &
                              ppd_product_of_bb_cc, &
                              expansion_atoms_out, expansion_centres, &
                              global_bb_ngwf_idx, global_cc_ngwf_idx, cell, &
                              coeffs_kind, dest_swex_h, 1)
                         write(stdout,'(a,i0,a,i0,a,i0,a,i0,a,e11.4)') &
                              'SWX: fit_rms_error: ', global_b,':',ngwf_b,',', &
                              global_c,':',ngwf_c,' ', rms
                      end if

                   end if ! fine/coarse grid product
                   ! ==========================================
                   ! ==========================================

                else
                   ! If present, then already expanded, so no-op
                end if

             end do loop_ngwf_b

             if(pub_swx_dbl_grid) call data_fftbox_dealloc(fftbox_ngwf_cc)

          end do loop_ngwf_c

          if(coeffs_kind == SW_W) then
             deallocate(wmatrix,stat=ierr)
             call utils_dealloc_check(myself,'wmatrix',ierr)
          end if
          if(coeffs_kind == SW_V .or. coeffs_kind == SW_W) then
             deallocate(vmatrix,stat=ierr)
             call utils_dealloc_check(myself,'vmatrix',ierr)
          end if
          if(coeffs_kind == SW_O .or. coeffs_kind == SW_W) then
             deallocate(omatrix,stat=ierr)
             call utils_dealloc_check(myself,'omatrix',ierr)
          end if

       end do loop_C

       ! cycle [*] gets here

    end do loop_B

    ! NB: Only gamma NGWFs are communicated
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
       if(pub_total_num_procs > 1) then
          do while(.not. done)
             if(.not. present(rep2)) then
                call remote_serve_ngwfs(done, who_is_done, &
                     rep%ngwfs_on_grid(ireg), ngwf_basis, rep%ngwf_cache_handle)
             else
                call remote_serve_ngwfs(done, who_is_done, &
                     rep%ngwfs_on_grid(ireg), ngwf_basis, rep%ngwf_cache_handle, &
                     ngwfs_on_grid2 = rep2%ngwfs_on_grid(ireg), &
                     ngwf_basis2 = ngwf_basis2, &
                     cache_handle2 = rep2%ngwf_cache_handle)
             end if
          end do
       end if

       do proc = 0, pub_total_num_procs-1
          call comms_wait(send_handles(proc))
       end do

       call comms_barrier
    end if

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine swx_local_ngwf_pair_fitting

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_my_ngwf_pair_fitting(beta_rep, gamma_rep, &        ! in/out, in
       my_bc_pairs, n_my_bc_pairs, packed_metric_matrix_blocks_hts, & ! in
       swri, dest_swex_h, cell, par, elements, coeffs_kind, &         ! in
       beta_ngwf_basis, gamma_ngwf_basis)                             ! in
    !==========================================================================!
    ! Performs ngwf-product-fitting, that is, for all pairs of S-overlapping   !
    ! atoms all NGWF products are expanded in terms of SWs originating on      !
    ! those two atoms. The resulting expansion coefficients are stored in      !
    ! coeffs_ht of a SW_EX container in beta_rep.                              !                             !
    !                                                                          !
    ! Caching of SWOPs in PPDs is performed behind the scenes (the cache is    !
    ! only read, not written to or cleaned up). The metric used for expansion  !
    ! is selected with coeffs_kind.                                            !
    !                                                                          !
    ! The pairs to expand are specified in the arguments. Neither B's nor C's  !
    ! are, in general, rank local. They are expected to be in the NGWF cache.  !
    !                                                                          !
    ! The double-grid version is not supported here, neither is swri-to-swex   !
    ! filtering.                                                               !
    !                                                                          !
    ! If the underlying SW_RI does not encompass the entire system, only the   !
    ! NGWF products between in-region and in-region and in-region and their    !
    ! overlapping atoms are expanded (ie. all products except those between    !
    ! an outsider and an outsider).                                            !
    !                                                                          !
    ! The calculated coefficients are *not* communicated across procs.         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   beta_rep (inout): The NGWF_REP containing B atoms. Results (coeffs)    !
    !                     are also stored here.                                !
    !   gamma_rep (in): The NGWF_REP containing C atoms. Not modified.         !
    !   my_bc_pairs (in): List of B-C pairs given to this proc in HFx paralle- !
    !                     lisation scheme.                                     !
    !   n_my_bc_pairs (in): The number of elements in the above.               !
    !   packed_metric_matrix_blocks_hts (in): Provides necessary metric matrix !
    !                                         blocks (since B-C pairs are not  !
    !                                         rank-local).                     !
    !   swri (in): The SW_RI underlying the expansion. SWOPs are read-cached   !
    !              behind the scenes.                                          !
    !   dest_swex_h (in): Handle for choosing which SW_EX container in beta_rep!
    !                will be used to store results. The container itself is not!
    !                passed to avoid aliasing with beta_rep.                   !
    !   cell, par, elements (in): The usual.                                   !
    !   coeffs_kind (in): See swx_expansion_2_fast().                          !
    !   beta_ngwf_basis (in): NGWF basis for B atoms.                          !
    !   gamma_ngwf_basis (in): NGWF basis for C atoms.                         !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in April 2019 using swx_local_pair_fitting()   !
    ! as a template.
    !==========================================================================!

    use constants, only: SW_V, SW_W, SW_O, SWEX_BB_CC_SYMMETRIC
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE, hash_table_list, hash_table_probe, &
         hash_table_lookup_ptr_nocount, hash_table_add
    use ion, only: ELEMENT
    use ngwf_representation, only: NGWF_REP
    use parallel_strategy, only: PARAL_INFO
    use ppd_ops, only: ppd_product
    use remote, only: remote_ngwf_cache_hts
    use rundat, only: pub_devel_code
    use simulation_cell, only: CELL_INFO
    use sw_resolution_of_identity, only: SW_RI, ATOM_CENTRE, swri_init_centre, &
         swri_build_matrix_2, swri_expansion_centres, swri_make_wmatrix_2, &
         swri_get_handle_to
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_alloc_check, &
         utils_dealloc_check, utils_assert, utils_devel_code

    implicit none

    ! jd: Arguments
    type(NGWF_REP), intent(inout), target   :: beta_rep
    type(NGWF_REP), intent(in), target      :: gamma_rep
    type(HT_HASH_TABLE), intent(in), target :: packed_metric_matrix_blocks_hts(:)
    integer, intent(in)                     :: my_bc_pairs(:,:)
    integer, intent(in)                     :: n_my_bc_pairs
    type(SW_RI), intent(inout)              :: swri
    integer, intent(in)                     :: dest_swex_h
    type(FUNC_BASIS), intent(in), target    :: beta_ngwf_basis
    type(FUNC_BASIS), intent(in), target    :: gamma_ngwf_basis
    type(CELL_INFO), intent(in)             :: cell
    type(PARAL_INFO), intent(in)            :: par
    type(ELEMENT), intent(in)               :: elements(par%nat)
    integer, intent(in)                     :: coeffs_kind

    ! jd: Local variables
    ! jd: Atoms and pairs
    integer :: global_b, global_c
    integer :: orig_global_b
    integer :: prev_global_b
    integer :: pair_idx

    ! jd: Centres and atoms in expansion
    integer :: num_sws_in_expansion
    integer :: expansion_atoms_out(2)
    integer :: expansion_atoms_in(2)
    type(ATOM_CENTRE) :: expansion_centres(2)
    type(ATOM_CENTRE) :: centres(2)

    ! jd: NGWFs
    integer :: n_ngwfs_b, n_ngwfs_c
    integer :: ngwf_b, ngwf_c
    integer :: first_ngwf_idx_of_b, first_ngwf_idx_of_c
    integer :: global_bb_ngwf_idx, global_cc_ngwf_idx
    real(kind=DP), pointer :: packed_ngwf_bb(:) ! ptr to ht's internals
    real(kind=DP), pointer :: packed_ngwf_cc(:) ! ptr to ht's internals
    integer :: packed_ngwf_bb_size, packed_ngwf_cc_size

    ! jd: Matrices in expansion
    real(kind=DP), allocatable :: vmatrix(:,:)
    real(kind=DP), allocatable :: omatrix(:,:)
    real(kind=DP), allocatable :: wmatrix(:,:)

    ! jd: BbCc product
    integer :: n_ppds_in_product_of_bb_cc
    integer :: ppd_indices_in_product_of_bb_cc(beta_ngwf_basis%max_n_ppds_sphere)
    real(kind=DP) :: ppd_product_of_bb_cc(beta_ngwf_basis%max_n_ppds_sphere*cell%n_pts)

    ! jd: Other
    real(kind=DP) :: coeffs(2*swri%quality%max_sws_per_centre)
    integer :: probe_result
    integer :: swri_handle
    logical :: debug_show_skipped
    real(kind=DP) :: rms
    integer :: ierr
    character(len=*), parameter :: myself = 'swx_my_ngwf_pair_fitting'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
         pub_devel_code, no_bcast=.true., no_warn=.true.)

    swri_handle = swri_get_handle_to(swri%swri_name)

    ! --------------------------------------------------------------------------
    ! jd: Loop over all my BC pairs                                           BC
    ! --------------------------------------------------------------------------
    prev_global_b = -1
    loop_BC:                                                                   &
    do pair_idx = 1, n_my_bc_pairs
       global_b = my_bc_pairs(1,pair_idx)
       global_c = my_bc_pairs(2,pair_idx)
       orig_global_b = par%orig_atom(global_b)
       first_ngwf_idx_of_b = beta_ngwf_basis%first_on_atom(global_b)
       first_ngwf_idx_of_c = gamma_ngwf_basis%first_on_atom(global_c)
       n_ngwfs_b = beta_ngwf_basis%num_on_atom(global_b)
       n_ngwfs_c = gamma_ngwf_basis%num_on_atom(global_c)

       ! jd: Ignore atoms not represented in this SWRI
       if(.not. elements(orig_global_b)%in_swri(swri_handle)) then
          if(debug_show_skipped) then
             write(stdout,'(a,i0,a,i0,a,a,a)') &
             'Skipping ('//myself//') orig/SFC: ', orig_global_b, '/', &
             global_b, ' (', trim(elements(orig_global_b)%species_id),')'
          end if
          cycle ! [*]
       end if

       call swri_init_centre(centres(1), par, elements, global_b)
       call swri_init_centre(centres(2), par, elements, global_c)

       ! jd: List hash table usage for SWOPS, coeffs, NGWFs.
       if(prev_global_b /= global_b) then
          call hash_table_list(remote_ngwf_cache_hts(beta_rep%ngwf_cache_handle),0)
          call hash_table_list(beta_rep%swexes(dest_swex_h)%coeffs_ht,0)
          call hash_table_list(swri%swops_in_ppds_ht,0)
          if(gamma_rep%ngwf_cache_handle /= beta_rep%ngwf_cache_handle) then
             call hash_table_list(remote_ngwf_cache_hts(gamma_rep%ngwf_cache_handle),0)
          end if
       end if

       ! jd: Work around 'temporary was created' debug mode warning
       expansion_atoms_in = (/global_b, global_c/)

       ! jd: Figure out the expansion centres, sort them, pad with -1
       !     if fewer than 2 are needed
       call swri_expansion_centres(swri, &
            expansion_centres, expansion_atoms_out, &  ! out
            num_sws_in_expansion, &                    ! out
            centres, expansion_atoms_in, &             ! in
            beta_rep%swexes(dest_swex_h)%quality, &    ! in
            SWEX_BB_CC_SYMMETRIC(dest_swex_h))         ! in

       ! jd: Allocate the required 2-centre matrices (O, V, W).
       if(coeffs_kind == SW_V .or. coeffs_kind == SW_W) then
          allocate(vmatrix(num_sws_in_expansion,num_sws_in_expansion),stat=ierr)
          call utils_alloc_check(myself,'vmatrix',ierr)
       end if
       if(coeffs_kind == SW_W) then
          allocate(wmatrix(num_sws_in_expansion,num_sws_in_expansion),stat=ierr)
          call utils_alloc_check(myself,'wmatrix',ierr)
       end if
       if(coeffs_kind == SW_O .or. coeffs_kind == SW_W) then
          allocate(omatrix(num_sws_in_expansion,num_sws_in_expansion),stat=ierr)
          call utils_alloc_check(myself,'omatrix',ierr)
       end if

       ! jd: Prepare the 2-centre V matrix
       if(coeffs_kind == SW_V .or. coeffs_kind == SW_W) then
          call swri_build_matrix_2(vmatrix, &                                  ! out
               packed_metric_matrix_blocks_hts(SW_V), swri, &                  ! in
               num_sws_in_expansion, expansion_centres, expansion_atoms_out, & ! in
               SW_V)
       end if

       ! jd: Prepare the 2-centre O matrix
       if(coeffs_kind == SW_O .or. coeffs_kind == SW_W) then
          call swri_build_matrix_2(omatrix, &                                 ! out
               packed_metric_matrix_blocks_hts(SW_O), swri, &                 ! in
               num_sws_in_expansion, expansion_centres, expansion_atoms_out,& !in
               SW_O)
       end if

       ! jd: Prepare the 2-centre W matrix from O and V
       if(coeffs_kind == SW_W) then
          call swri_make_wmatrix_2(wmatrix, vmatrix, omatrix, &
               num_sws_in_expansion)
       end if

       ! -----------------------------------------------------------------------
       ! jd: for all b on B                                                  bbb
       ! -----------------------------------------------------------------------
       loop_ngwf_b:                                                            &
       do ngwf_b = 1, n_ngwfs_b
          global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

          ! jd: Get NGWF Bb
          call hash_table_lookup_ptr_nocount(packed_ngwf_bb, &
               packed_ngwf_bb_size, &
               remote_ngwf_cache_hts(beta_rep%ngwf_cache_handle), &
               global_bb_ngwf_idx)
          call utils_assert(packed_ngwf_bb_size /= -1, &
               myself//': NGWF Bb absent in cache.', &
               global_bb_ngwf_idx, global_b, ngwf_b)

          ! --------------------------------------------------------------------
          ! jd: for all c on C                                               ccc
          ! --------------------------------------------------------------------
          loop_ngwf_c:                                                         &
          do ngwf_c = 1, n_ngwfs_c
             global_cc_ngwf_idx = first_ngwf_idx_of_c + ngwf_c - 1

             ! jd: Get NGWF Cc
             call hash_table_lookup_ptr_nocount(packed_ngwf_cc, &
                  packed_ngwf_cc_size, &
                  remote_ngwf_cache_hts(gamma_rep%ngwf_cache_handle), &
                  global_cc_ngwf_idx)
             call utils_assert(packed_ngwf_cc_size /= -1, &
                  myself//': NGWF Cc absent in cache.', &
                  global_cc_ngwf_idx, global_c, ngwf_c)

             ! jd: In a single NGWF basis scenario, expansions for (Bb,Cc)
             !     and (Cc,Bb) are identical (symmetry under label exchange).
             !     But we cannot just expand for Bb <= Cc, when some Bs can
             !     be skipped due to not being in the SWRI region. We need to
             !     handle overlaps of in-region Bb's with out-of-region Cc's
             !     regardless of ordering. We thus explicitly check if the
             !     coeff is there, and expand if not.

             ! jd: Check if we have coefficients for this Bb Cc
             if(.not. SWEX_BB_CC_SYMMETRIC(dest_swex_h)) then
                ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                probe_result = hash_table_probe(&
                     beta_rep%swexes(dest_swex_h)%coeffs_ht, & ! in
                     global_bb_ngwf_idx, global_cc_ngwf_idx, & ! in
                     coeffs_kind, dest_swex_h, 1)
             else
                ! jd: Here there's Bb-Cc symmetry and we only store Bb<=Cc triangle
                probe_result = hash_table_probe(&
                     beta_rep%swexes(dest_swex_h)%coeffs_ht, &       ! in
                     min(global_bb_ngwf_idx,global_cc_ngwf_idx), &   ! in
                     max(global_bb_ngwf_idx,global_cc_ngwf_idx), &   ! in
                     coeffs_kind, dest_swex_h, 1)
             end if

             if(probe_result == -1) then

                ! *** If absent, expand ***
                ! ==========================================
                ! ============ EXPANSION PROPER ============
                ! ==========================================
                ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ! ~~~ coarse grid (usual) version ~~~
                ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ! jd: Compute product of NGWFs Bb and Cc -- coarse grid
                call ppd_product(packed_ngwf_bb, packed_ngwf_cc, cell%n_pts, & ! in
                     ppd_product_of_bb_cc, &                         ! out
                     ppd_indices_in_product_of_bb_cc, &              ! out
                     n_ppds_in_product_of_bb_cc)                     ! out

                ! jd: Expand product in SWOPs -- coarse grid
                call swx_expansion_2_fast(coeffs, &            ! out
                     swri, expansion_centres, expansion_atoms_out, & ! in
                     num_sws_in_expansion, coeffs_kind, &      ! in
                     cell, vmatrix, wmatrix, omatrix, &        ! in
                     ppd_quantity = ppd_product_of_bb_cc, &    ! in
                     n_ppds_in_quantity = n_ppds_in_product_of_bb_cc, & ! in
                     ppd_indices_in_quantity = ppd_indices_in_product_of_bb_cc) ! in

                if(.not. SWEX_BB_CC_SYMMETRIC(dest_swex_h)) then
                   ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                   call hash_table_add(beta_rep%swexes(dest_swex_h)%coeffs_ht, &
                        coeffs, num_sws_in_expansion, &
                        global_bb_ngwf_idx, global_cc_ngwf_idx, &
                        coeffs_kind, dest_swex_h, 1, overfill_strategy = 'F')
                else
                   call hash_table_add(beta_rep%swexes(dest_swex_h)%coeffs_ht, &
                        coeffs, num_sws_in_expansion, &
                        min(global_bb_ngwf_idx,global_cc_ngwf_idx), &
                        max(global_bb_ngwf_idx,global_cc_ngwf_idx), &
                        coeffs_kind, dest_swex_h, 1, overfill_strategy = 'F')
                end if

                if(swx_calculate_rms_fitting_error) then
                   rms = swx_rms_fitting_error(beta_rep%swexes(dest_swex_h), &
                        swri, n_ppds_in_product_of_bb_cc, &
                        ppd_indices_in_product_of_bb_cc, ppd_product_of_bb_cc, &
                        expansion_atoms_out, expansion_centres, &
                        global_bb_ngwf_idx, global_cc_ngwf_idx, cell, &
                        coeffs_kind, dest_swex_h, 1)
                   write(stdout,'(a,i0,a,i0,a,i0,a,i0,a,e11.4)') &
                        'SWX: fit_rms_error: ', global_b,':',ngwf_b,',', &
                        global_c,':',ngwf_c,' ', rms
                end if
                ! ==========================================
                ! ==========================================

             else
                ! If present, then already expanded, so no-op
             end if

          end do loop_ngwf_c

       end do loop_ngwf_b

       if(coeffs_kind == SW_W) then
          deallocate(wmatrix,stat=ierr)
          call utils_dealloc_check(myself,'wmatrix',ierr)
       end if
       if(coeffs_kind == SW_V .or. coeffs_kind == SW_W) then
          deallocate(vmatrix,stat=ierr)
          call utils_dealloc_check(myself,'vmatrix',ierr)
       end if
       if(coeffs_kind == SW_O .or. coeffs_kind == SW_W) then
          deallocate(omatrix,stat=ierr)
          call utils_dealloc_check(myself,'omatrix',ierr)
       end if

       prev_global_b = global_b

    end do loop_BC

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine swx_my_ngwf_pair_fitting

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_local_density_fitting(coeffs_ht, rep, ireg, swri, swex_h, &
       denskern, ngwf_basis, cell, fftbox, elements, coeffs_kind, &
       vacuum_rep, vacuum_denskern)
    !==========================================================================!
    ! Performs density-fitting, that is, for all pairs of S-overlapping atoms  !
    ! all NGWF products multiplied by a kernel element are expanded in terms   !
    ! of SWs originating on those two atoms. The resulting expansion           !
    ! coefficients are stored the coeffs_ht argument (*not* in rep's SW_EX!).  !
    ! The handle swex_h is used to distinguish multiple sets (e.g. from Bessel !
    ! averaging). The coefficients are spin-dependent, and spin will be the    !
    ! fifth index in the returned HT.                                          !
    !                                                                          !
    ! The pairs to expand are rank-local-B's with their s-neighbours C's.      !
    ! C's are not necessarily rank-local.                                      !
    !                                                                          !
    ! Caching of SWOPs in PPDs is performed behind the scenes. The metric used !
    ! for expansion is selected with coeffs_kind.                              !
    !                                                                          !
    ! If the underlying SW_RI does not encompass the entire system, only the   !
    ! NGWF products between in-region and in-region and in-region and their    !
    ! overlapping atoms are expanded (ie. all products except those between    !
    ! an outsider and an outsider).                                            !
    !                                                                          !
    ! The calculated coefficients are not communicated across procs, this is   !
    ! done in swx_expand_pair_densities().                                     !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   coeffs_ht (in/out): The expansion coefficients are stored here. It is  !
    !                       the responsibility of the caller to pass a suitably!
    !                       initialised hash table.                            !
    !   rep (in): The NGWF_REP for which to perform density fitting.           !
    !   swri (inout): The SW_RI underlying the expansion. SWOPs are cached     !
    !                 behind the scenes.                                       !
    !   swex_h (in): This will be stored in the fourth key of the hash table.  !
    !                In contrast to swx_ngwf_pair_fitting(), this will also be !
    !                used to choose a SW_EX from rep, so only pass appropriate !
    !                SW_EX handles like REP_SWEX_POL_EMB_DMA_1. The quality of !
    !                the expansion will be determined by                       !
    !                rep%swexes(swex_h)%quality.                               !
    !   ngwf_basis, cell, elements (in): The usual.                            !
    !   fftbox (in): Only used when swx_double_grid T.                         !
    !   coeffs_kind (in): See swx_expansion_2() or swx_expand_pair_densities().!
    !--------------------------------------------------------------------------!
    ! Expert mode:                                                             !
    !   When arguments 'vacuum_rep' and 'vacuum_denskern' are specified, the   !
    !   expanded quantity will be the density difference, ie.                  !
    !   \phi_Aa K^Aa,Bb \phi_Bb - \phi0_Aa K0^Aa,Bb \phi_Bb, where \phi0 refer !
    !   to vacuum_rep, and K0 refers to vacuum_denskern.                       !
    !   CAVEAT: The quality of the expansion will be determined from           !
    !           rep%swexes(swex_h)%quality, vacuum_rep's quality will be       !
    !           ignored.
    !--------------------------------------------------------------------------!
    ! Caveat:                                                                  !
    !   This subroutine is *not* ready for mixed NGWF bases. I don't envision  !
    !   hybrid polarisable embedding conduction with vacuum reference states   !
    !   to ever be needed. If it turns out I'm wrong, modify this subroutine   !
    !   in the same fashion swx_ngwf_pair_fitting() was modified in June-August!
    !   2018.                                                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in June 2017 using swx_ngwf_pair_fitting() as  !
    ! a template.                                                              !
    ! Modified for embedding by Joseph Prentice, September 2018                !
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_send, &
         comms_barrier, comms_wait
    use constants, only: SW_V, SW_W, SW_O, VERBOSE, SWEX_BB_CC_SYMMETRIC
    use datatypes, only: FFTBOX_DATA, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT
    use hash_table, only: HT_HASH_TABLE, hash_table_list, hash_table_probe
    use ion, only: ELEMENT
    use ngwf_representation, only: NGWF_REP
    use parallel_strategy, only: PARAL_INFO
    use ppd_ops, only: ppd_product
    use remote, only: packed_ngwf_size, NGWF_REQUEST_TAG, remote_obtain_ngwf, &
         remote_serve_ngwfs, remote_pack_ngwf, remote_unpack_ngwf_to_fftbox, &
         remote_ngwf_cache_hts
    use rundat, only: pub_devel_code, pub_swx_output_detail, pub_swx_dbl_grid, &
         pub_pol_emb_dbl_grid, pub_num_spins
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_get_element, sparse_get_par
    use sw_resolution_of_identity, only: SW_RI, ATOM_CENTRE, swri_init_centre, &
         swri_extract_matrix_2, swri_expansion_centres, swri_make_wmatrix_2, &
         swri_get_handle_to
    use timer, only: timer_clock
    use utils, only: utils_trace_in, utils_trace_out, utils_alloc_check, &
         utils_dealloc_check, utils_assert, utils_devel_code, utils_abort

    implicit none

    ! jd: Arguments
    type(HT_HASH_TABLE), intent(inout), target :: coeffs_ht
    type(NGWF_REP), intent(in)                 :: rep
    type(SW_RI), intent(inout)                 :: swri
    integer, intent(in)                        :: swex_h
    type(SPAM3), intent(in)                    :: denskern(pub_num_spins)
    type(FUNC_BASIS), intent(in)               :: ngwf_basis
    type(CELL_INFO), intent(in)                :: cell
    type(FFTBOX_INFO), intent(in)              :: fftbox
    type(ELEMENT), intent(in)                  :: elements(:)
    integer, intent(in)                        :: coeffs_kind
    integer, intent(in)                        :: ireg
    type(NGWF_REP), intent(inout), optional    :: vacuum_rep
    type(SPAM3), intent(in), optional          :: vacuum_denskern(pub_num_spins)

    ! jd: --- Local variables ---

    ! jd: Atom-counting and indices
    integer :: local_b
    integer :: global_b, global_c, c_idx, orig_global_b
    type(ATOM_CENTRE) :: centres(2)

    ! jd: NGWF-counting and indices
    integer :: local_bb_ngwf_idx
    integer :: global_bb_ngwf_idx, global_cc_ngwf_idx
    integer :: first_ngwf_idx_of_b, first_ngwf_idx_of_c
    integer :: local_idx_of_1st_ngwf_on_b, local_idx_of_1st_ngwf_on_c
    integer :: ngwf_b, ngwf_c

    ! jd: Centres and atoms in expansion
    integer :: num_sws_in_expansion
    integer :: expansion_atoms_out(2)
    integer :: expansion_atoms_in(2)
    type(ATOM_CENTRE) :: expansion_centres(2)

    ! jd: Matrices in expansion
    real(kind=DP), allocatable :: vmatrix(:,:)
    real(kind=DP), allocatable :: omatrix(:,:)
    real(kind=DP), allocatable :: wmatrix(:,:)

    ! jd: PPDs and PPD indices
    integer :: n_ppds_in_product_of_bb_cc
    integer :: ppd_indices_in_product_of_bb_cc(ngwf_basis%max_n_ppds_sphere)
    real(kind=DP) :: ppd_product_of_bb_cc(PPDS_DIMS)
    real(kind=DP) :: ppd_product_of_bb_cc_dkn(PPDS_DIMS)
    integer :: vacuum_n_ppds_in_product_of_bb_cc
    integer :: vacuum_ppd_indices_in_product_of_bb_cc(ngwf_basis%max_n_ppds_sphere)
    real(kind=DP) :: vacuum_ppd_product_of_bb_cc(PPDS_DIMS)
    real(kind=DP) :: vacuum_ppd_product_of_bb_cc_dkn(PPDS_DIMS)
    real(kind=DP) :: ppd_density(PPDS_DIMS)

    ! jd: NGWFs in FFT boxes
    type(FFTBOX_DATA) :: fftbox_ngwf_bb, fftbox_ngwf_cc, &
         product_of_bb_cc_in_b_fftbox
    type(FFTBOX_DATA) :: vacuum_fftbox_ngwf_bb, vacuum_fftbox_ngwf_cc, &
         vacuum_product_of_bb_cc_in_b_fftbox
    type(POINT) :: fftbox_bb_start_in_cell_vec
    integer :: delta_mid_complement(3)

    ! jd: DKN stuff
    real(kind=DP) :: dkn_el, vacuum_dkn_el
    integer :: is

    ! jd: Packed NGWFs and their comms
    real(kind=DP)          :: packed_ngwf_bb(packed_ngwf_size) ! local
    real(kind=DP)          :: packed_vacuum_ngwf_bb(packed_ngwf_size) ! local
    real(kind=DP), pointer :: packed_ngwf_cc(:) ! remote
    real(kind=DP), pointer :: packed_vacuum_ngwf_cc(:) ! remote
    logical :: done
    integer :: proc
    integer :: send_handles(0:pub_total_num_procs-1)
    integer :: vacuum_send_handles(0:pub_total_num_procs-1)
    logical :: who_is_done(0:pub_total_num_procs-1)
    logical :: vacuum_who_is_done(0:pub_total_num_procs-1)

    ! jd: Varia
    integer :: probe_result
    integer :: ierr
    real(kind=DP) :: rms
    integer :: swri_handle
    logical :: debug_show_skipped
    character(len=*), parameter :: myself = 'swx_local_density_fitting'
    type(PARAL_INFO), pointer   :: par

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    call utils_assert(present(vacuum_rep) .eqv. present(vacuum_denskern), &
         myself//'Optional arguments must all be specified or all be absent')

    ! rc2013: get the parallel strategy from the sparse matrix
    call sparse_get_par(par, denskern(1))
    call utils_assert(par%nat == size(elements), 'Error in '//myself//': &
         &allocated parallel strategy is incompatible with elements.')

    swri_handle = swri_get_handle_to(swri%swri_name)

    debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
         pub_devel_code, no_bcast=.true., no_warn=.true.)

    ! jd: Sanity check on neighbour lists
    call utils_assert(swri%s_atoms_nl%populated, &
         'Internal error in '//trim(myself)&
         //': s_atoms neighbour list not populated for '//trim(swri%swri_name))

    vacuum_who_is_done(:) = .false.
    who_is_done(:) = .false.

    ! --------------------------------------------------------------------------
    ! jd: Loop over B's -- all atoms on this core                            BBB
    ! --------------------------------------------------------------------------
    loop_B:                                                                    &
    do local_b=1, par%num_atoms_on_proc(pub_my_proc_id)
       global_b = par%first_atom_on_proc(pub_my_proc_id) + local_b-1
       orig_global_b = par%orig_atom(global_b)

       ! jd: Ignore atoms not represented in this SWRI
       if(.not. elements(orig_global_b)%in_swri(swri_handle)) then
          if(debug_show_skipped) then
             write(stdout,'(a,i0,a,i0,a,a,a)') &
             'Skipping ('//myself//') orig/SFC: ', orig_global_b, '/', &
             global_b, ' (', trim(elements(orig_global_b)%species_id),')'
          end if
          cycle ! [*]
       end if

       if(pub_swx_output_detail >= VERBOSE) write(stdout,*) ' - B: ', global_b

       call swri_init_centre(centres(1), par, elements, global_b)

       call hash_table_list(remote_ngwf_cache_hts(rep%ngwf_cache_handle),0)
       call hash_table_list(coeffs_ht,0)
       if(present(vacuum_rep)) then
          call hash_table_list(remote_ngwf_cache_hts(vacuum_rep%ngwf_cache_handle),0)
       end if
       call hash_table_list(swri%swops_in_ppds_ht,0)
       if(pub_swx_dbl_grid .or. pub_pol_emb_dbl_grid) then
          call hash_table_list(swri%swops_at_points_ht,0)
       end if

       first_ngwf_idx_of_b = ngwf_basis%first_on_atom(global_b)
       local_idx_of_1st_ngwf_on_b = first_ngwf_idx_of_b + 1 - &
           ngwf_basis%first_on_proc(pub_my_proc_id)

       ! -----------------------------------------------------------------------
       ! jd: Loop over C's that are s-neighbours with B                      CCC
       ! -----------------------------------------------------------------------
       loop_C:                                                                 &
       do c_idx = swri%s_atoms_nl%first_idx(global_b), &
            swri%s_atoms_nl%last_idx(global_b)

          global_c = swri%s_atoms_nl%neighbours(c_idx)
          ! jd: No elimination of atoms outside SWRI here, we do want the
          !     input from BC if B is in SWRI and C is outside

          call swri_init_centre(centres(2), par, elements, global_c)

          first_ngwf_idx_of_c = ngwf_basis%first_on_atom(global_c)
          local_idx_of_1st_ngwf_on_c = first_ngwf_idx_of_c + 1 - &
              ngwf_basis%first_on_proc(pub_my_proc_id)

          ! jd: Work around 'temporary was created' debug mode warning
          expansion_atoms_in = (/global_b, global_c/)

          ! jd: Figure out the expansion centres, sort them, pad with -1
          !     if fewer than 2 are needed
          call swri_expansion_centres(swri, &
               expansion_centres, expansion_atoms_out, &  ! out
               num_sws_in_expansion, &                    ! out
               centres, expansion_atoms_in, rep%swexes(swex_h)%quality, & ! in
               SWEX_BB_CC_SYMMETRIC(swex_h))              ! in

          ! jd: Allocate the required 2-centre matrices (O, V, W).
          if(coeffs_kind == SW_V .or. coeffs_kind == SW_W) then
             allocate(vmatrix(num_sws_in_expansion,num_sws_in_expansion), &
                  stat=ierr)
             call utils_alloc_check(myself,'vmatrix',ierr)
          end if
          if(coeffs_kind == SW_W) then
             allocate(wmatrix(num_sws_in_expansion,num_sws_in_expansion), &
                  stat=ierr)
             call utils_alloc_check(myself,'wmatrix',ierr)
          end if
          if(coeffs_kind == SW_O .or. coeffs_kind == SW_W) then
             allocate(omatrix(num_sws_in_expansion,num_sws_in_expansion), &
                  stat=ierr)
             call utils_alloc_check(myself,'omatrix',ierr)
          end if

          ! jd: Prepare the 2-centre V matrix
          if(coeffs_kind == SW_V .or. coeffs_kind == SW_W) then
             call swri_extract_matrix_2(vmatrix, &                ! out
                  swri, num_sws_in_expansion, &                   ! in
                  expansion_centres, expansion_atoms_out, SW_V, & ! in
                  rep%swexes(swex_h)%quality)                     ! in
          end if

          ! jd: Prepare the 2-centre O matrix
          if(coeffs_kind == SW_O .or. coeffs_kind == SW_W) then
             call swri_extract_matrix_2(omatrix, &                ! out
                  swri, num_sws_in_expansion, &                   ! in
                  expansion_centres, expansion_atoms_out, SW_O, & ! in
                  rep%swexes(swex_h)%quality)                     ! in
          end if

          ! jd: Prepare the 2-centre W matrix from O and V
          if(coeffs_kind == SW_W) then
             call swri_make_wmatrix_2(wmatrix, vmatrix, omatrix, &
                  num_sws_in_expansion)
          end if

          ! --------------------------------------------------------------------
          ! jd: for all c on C                                               ccc
          ! --------------------------------------------------------------------
          loop_ngwf_c:                                                         &
          do ngwf_c = 1, ngwf_basis%num_on_atom(global_c)
             global_cc_ngwf_idx = first_ngwf_idx_of_c + ngwf_c - 1

             ! @optimizeme: check first if we perhaps have all coefficients
             ! for this c, ie for all b's. If so, no need to obtain and expand Cc
             ! since no-one does comms during swx_expansion_2, there are likely
             ! convoy effects here now.
             call remote_obtain_ngwf(packed_ngwf_cc, &
                  remote_ngwf_cache_hts, &
                  who_is_done, global_cc_ngwf_idx, rep%ngwfs_on_grid(ireg), &
                  ngwf_basis, rep%ngwf_cache_handle, overfill_strategy='F')
             if(present(vacuum_rep)) then
                 call utils_abort(myself//': Vacuum rep functionality &
                      &temporarily disabled due to changes in remote')
!                call remote_obtain_ngwf(packed_vacuum_ngwf_cc, &
!                     remote_ngwf_cache_hts, &
!                     vacuum_who_is_done, global_cc_ngwf_idx, &
!                     vacuum_rep%ngwfs_on_grid, ngwf_basis, &
!                     vacuum_rep%ngwf_cache_handle, )
             end if
             if(pub_swx_dbl_grid) then
                ! jd: Unpack Cc NGWF to an FFT box midway between Aa and Bb
                call remote_unpack_ngwf_to_fftbox(fftbox_ngwf_cc, &
                     packed_ngwf_cc, fftbox, cell, &
                     centre_for_target_fftbox = &
                     ngwf_basis%spheres(local_idx_of_1st_ngwf_on_b)%centre, &
                     delta_mid_complement = delta_mid_complement)
                if(present(vacuum_rep)) then
                   call remote_unpack_ngwf_to_fftbox(vacuum_fftbox_ngwf_cc, &
                        packed_vacuum_ngwf_cc, fftbox, cell, &
                        centre_for_target_fftbox = &
                        ngwf_basis%spheres(local_idx_of_1st_ngwf_on_b)%centre, &
                        delta_mid_complement = delta_mid_complement)
                end if
             end if

             ! -----------------------------------------------------------------
             ! jd: for all b on B                                            bbb
             ! -----------------------------------------------------------------
             loop_ngwf_b:                                                         &
             do ngwf_b = 1, ngwf_basis%num_on_atom(global_b)
                global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1
                local_bb_ngwf_idx = global_bb_ngwf_idx - &
                      ngwf_basis%first_on_proc(pub_my_proc_id)+1

                ! jd: Expansions for (Bb,Cc) and (Cc,Bb) are identical. But we
                !     cannot just expand for Bb <= Cc, when some Bs can be
                !     skipped due to not being in the SWRI region. We need to
                !     handle overlaps of in-region Bb's with out-of-region Cc's
                !     regardless of ordering. We thus explicitly check if the
                !     coeff is there, and expand if not.
                ! jd: Also, when using mixed NGWF bases there's no Bb-Cc symmetry.

                ! jd: Check if we have coefficients for this Bb Cc. It is
                !     sufficient to check spin 1.
                if(.not. SWEX_BB_CC_SYMMETRIC(swex_h)) then
                   ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                   probe_result = hash_table_probe(coeffs_ht, &     ! in
                        global_bb_ngwf_idx, global_cc_ngwf_idx, &   ! in
                        coeffs_kind, swex_h, 1)
                else
                   ! jd: Here there's Bb-Cc symmetry and we only store one triangle
                   probe_result = hash_table_probe(coeffs_ht, &  ! in
                        min(global_bb_ngwf_idx,global_cc_ngwf_idx), &   ! in
                        max(global_bb_ngwf_idx,global_cc_ngwf_idx), &   ! in
                        coeffs_kind, swex_h, 1)
                end if
                if(probe_result == -1) then

                   ! *** If absent, expand ***

                   ! jd: Pack NGWF Bb for ppd_product
                   call remote_pack_ngwf(packed_ngwf_bb, &
                        rep%ngwfs_on_grid(ireg), ngwf_basis, local_bb_ngwf_idx)

                   if(present(vacuum_rep)) then
                      call remote_pack_ngwf(packed_vacuum_ngwf_bb, &
                           vacuum_rep%ngwfs_on_grid(ireg), ngwf_basis, &
                           local_bb_ngwf_idx)
                   end if

                   ! ==========================================
                   ! ============ EXPANSION PROPER ============
                   ! ==========================================
                   if(pub_swx_dbl_grid) then
                      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      ! ~~~ double grid version ~~~
                      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      call utils_abort(myself//': density fitting not &
                           &implemented with swx_dbl_grid yet')
#if 0
                      ! jd: @must loop over spins here like in the coarse grid
                      !     case
                      !     @must take into account vacuum rep like in the
                      !     coarse grid case
                      ! jd: Compute product of NGWFs Bb and Cc -- double grid
                      call remote_unpack_ngwf_to_fftbox(fftbox_ngwf_bb, &
                           packed_ngwf_bb, fftbox, cell, &
                           fftbox_start_in_cell_vec = fftbox_bb_start_in_cell_vec, &
                           gridpoint_offset = delta_mid_complement)
                      call swx_ngwf_product(product_of_bb_cc_in_b_fftbox, & !out
                           fftbox_ngwf_bb, fftbox_ngwf_cc, fftbox)
                      call data_fftbox_dealloc(fftbox_ngwf_bb)

                      ! jd: Expand product in SWOPs -- double grid
                      if(.not. SWEX_BB_CC_SYMMETRIC(swex_h)) then
                         ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                         call swx_expansion_2(swri, &                   ! in/out
                              rep%swexes(swex_h)%quality, swex_h, is, & ! in
                              expansion_centres, expansion_atoms_out, & ! in
                              num_sws_in_expansion, coeffs_kind, &      ! in
                              cell, vmatrix, wmatrix, omatrix, &        ! in
                              global_bb_ngwf_idx, global_cc_ngwf_idx, & !} for caching
                              coeffs_ht, &                              ! in/out
                              fftbox_density = product_of_bb_cc_in_b_fftbox, &
                              fftbox = fftbox, &
                              fftbox_start_in_cell_vec = fftbox_bb_start_in_cell_vec)
                       else
                         ! jd: Here there's Bb-Cc symmetry and we only store one triangle
                         call swx_expansion_2(swri, &                       ! in/out
                              rep%swexes(swex_h)%quality, swex_h, is, &     ! in
                              expansion_centres, expansion_atoms_out, &     ! in
                              num_sws_in_expansion, coeffs_kind, &          ! in
                              cell, vmatrix, wmatrix, omatrix, &            ! in
                              min(global_bb_ngwf_idx,global_cc_ngwf_idx), & !} for
                              max(global_bb_ngwf_idx,global_cc_ngwf_idx), & !} caching
                              coeffs_ht, &                                  ! in/out
                              fftbox_density = product_of_bb_cc_in_b_fftbox, &
                              fftbox = fftbox, &
                              fftbox_start_in_cell_vec = fftbox_bb_start_in_cell_vec)
                       end if
                       call data_fftbox_dealloc(product_of_bb_cc_in_b_fftbox)
#endif
                   else
                      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      ! ~~~ coarse grid (usual) version ~~~
                      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      ! jd: Compute product of NGWFs Bb and Cc -- coarse grid
                      call ppd_product(packed_ngwf_bb, packed_ngwf_cc, &
                           cell%n_pts, &                       ! in
                           ppd_product_of_bb_cc, &             ! out
                           ppd_indices_in_product_of_bb_cc, &              ! out
                           n_ppds_in_product_of_bb_cc)                     ! out
                      ! jd: Compute product of vacuum NGWFs Bb and Cc -- coarse grid
                      if(present(vacuum_rep)) then
                         call ppd_product(packed_vacuum_ngwf_bb, &
                              packed_vacuum_ngwf_cc, cell%n_pts, &   ! in
                              vacuum_ppd_product_of_bb_cc, &         ! out
                              vacuum_ppd_indices_in_product_of_bb_cc, &       ! out
                              vacuum_n_ppds_in_product_of_bb_cc)              ! out
                         call utils_assert(vacuum_n_ppds_in_product_of_bb_cc == &
                              n_ppds_in_product_of_bb_cc, myself//': n_ppds_in_&
                              &product_of_bb_cc inconsistent between vacuum and&
                              & current state ',vacuum_n_ppds_in_product_of_bb_cc, &
                              n_ppds_in_product_of_bb_cc)
                         call utils_assert(all(&
                              vacuum_ppd_indices_in_product_of_bb_cc(1:&
                              vacuum_n_ppds_in_product_of_bb_cc) == &
                              ppd_indices_in_product_of_bb_cc(1:&
                              n_ppds_in_product_of_bb_cc)), myself//&
                              ': ppd_indices_in_product_of_bb_cc inconsistent &
                              &between vacuum and current state')
                      end if

                      ! --------------------------------------------------------
                      ! jd: for all spins                                    sss
                      ! --------------------------------------------------------
                      loop_spin:                                               &
                      do is = 1, pub_num_spins
                         ! jd: Obtain DKN element
                         call sparse_get_element(dkn_el, denskern(is), &
                              global_cc_ngwf_idx, global_bb_ngwf_idx)
                         if(present(vacuum_rep)) then
                            call sparse_get_element(vacuum_dkn_el, &
                                 vacuum_denskern(is), &
                                 global_cc_ngwf_idx, global_bb_ngwf_idx)
                         end if

                         ! jd: Obtain atom-pair density
                         ppd_product_of_bb_cc_dkn = ppd_product_of_bb_cc * dkn_el
                         if(present(vacuum_rep)) then
                            vacuum_ppd_product_of_bb_cc_dkn = &
                                 vacuum_ppd_product_of_bb_cc * vacuum_dkn_el
                            ppd_density = ppd_product_of_bb_cc_dkn - &
                                 vacuum_ppd_product_of_bb_cc_dkn
                         else
                            ppd_density = ppd_product_of_bb_cc_dkn
                         end if

                         ! jd: Expand density difference in SWOPs -- coarse grid
                         if(.not. SWEX_BB_CC_SYMMETRIC(swex_h)) then
                            ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                            call swx_expansion_2(swri, &                    ! in/out
                                 rep%swexes(swex_h)%quality, swex_h, is, &  ! in
                                 expansion_centres, expansion_atoms_out, &  ! in
                                 num_sws_in_expansion, coeffs_kind, &       ! in
                                 cell, vmatrix, wmatrix, omatrix, &         ! in
                                 global_bb_ngwf_idx, global_cc_ngwf_idx, &  !} for caching
                                 coeffs_ht, &                               ! in/out
                                 ppd_quantity = ppd_density, &
                                 n_ppds_in_quantity = n_ppds_in_product_of_bb_cc, & ! in
                                 ppd_indices_in_quantity = ppd_indices_in_product_of_bb_cc) ! in
                         else
                            ! jd: Here there's Bb-Cc symmetry and we only store one triangle
                            call swx_expansion_2(swri, &                       ! in/out
                                 rep%swexes(swex_h)%quality, swex_h, is, &     ! in
                                 expansion_centres, expansion_atoms_out, &     ! in
                                 num_sws_in_expansion, coeffs_kind, &          ! in
                                 cell, vmatrix, wmatrix, omatrix, &            ! in
                                 min(global_bb_ngwf_idx,global_cc_ngwf_idx), & !} for
                                 max(global_bb_ngwf_idx,global_cc_ngwf_idx), & !} caching
                                 coeffs_ht, &                                  ! in/out
                                 ppd_quantity = ppd_density, &
                                 n_ppds_in_quantity = n_ppds_in_product_of_bb_cc, & ! in
                                 ppd_indices_in_quantity = ppd_indices_in_product_of_bb_cc) ! in
                         end if

                         if(swx_calculate_rms_fitting_error) then
                            rms = swx_rms_fitting_error(rep%swexes(swex_h), swri, &
                                 n_ppds_in_product_of_bb_cc, &
                                 ppd_indices_in_product_of_bb_cc, &
                                 ppd_product_of_bb_cc, &
                                 expansion_atoms_out, expansion_centres, &
                                 global_bb_ngwf_idx, global_cc_ngwf_idx, cell, &
                                 coeffs_kind, swex_h, is)
                            write(stdout,'(a,i0,a,i0,a,i0,a,i0,a,e11.4)') &
                                 'SWX: fit_rms_error: ', global_b,':',ngwf_b,',', &
                                 global_c,':',ngwf_c,' ', rms
                         end if ! calculate_rms_fitting_error
                     end do loop_spin ! spins

                   end if ! fine/coarse grid product
                   ! ==========================================
                   ! ==========================================

                else
                   ! If present, then already expanded, so no-op
                end if

             end do loop_ngwf_b

             if(pub_swx_dbl_grid) call data_fftbox_dealloc(fftbox_ngwf_cc)

          end do loop_ngwf_c

          if(coeffs_kind == SW_W) then
             deallocate(wmatrix,stat=ierr)
             call utils_dealloc_check(myself,'wmatrix',ierr)
          end if
          if(coeffs_kind == SW_V .or. coeffs_kind == SW_W) then
             deallocate(vmatrix,stat=ierr)
             call utils_dealloc_check(myself,'vmatrix',ierr)
          end if
          if(coeffs_kind == SW_O .or. coeffs_kind == SW_W) then
             deallocate(omatrix,stat=ierr)
             call utils_dealloc_check(myself,'omatrix',ierr)
          end if

       end do loop_C

       ! cycle [*] gets here

    end do loop_B

    if(pub_total_num_procs > 1) then
       ! jd: After we're done, first notify everyone, then keep serving
       !     other procs with NGWFs, until everyone is done
       do proc = 0, pub_total_num_procs-1
          call comms_send(proc, -1, 1, tag = NGWF_REQUEST_TAG + 1, &
               return_handle = send_handles(proc), add_to_stack = .false.)
               ! [x]: For a 'done' request it's OK to only send it for
               !      NGWF set 1 -- it's understood that *all* sets finished
               !      participating.

          if(present(vacuum_rep)) then
              call utils_abort(myself//': Vacuum rep functionality &
                   &temporarily disabled due to changes in remote')
!             call comms_send(proc, -1, 1, tag = NGWF_REQUEST_TAG + VACUUM_TAG, &
!                  return_handle = vacuum_send_handles(proc), &
!                  add_to_stack = .false.)
          end if
       end do

       done = .false.
       if(pub_total_num_procs > 1) then
          do while(.not. done)
             call remote_serve_ngwfs(done, who_is_done, &
                  rep%ngwfs_on_grid(ireg), ngwf_basis, rep%ngwf_cache_handle)
             if(present(vacuum_rep)) then
                 call utils_abort(myself//': Vacuum rep functionality &
                      &temporarily disabled due to changes in remote')
!                call remote_serve_ngwfs_and_dkn_blocks(done, &
!                     vacuum_who_is_done, vacuum_rep%ngwfs_on_grid, ngwf_basis, &
!                     user_tag = VACUUM_TAG)
             end if
          end do
       end if

       do proc = 0, pub_total_num_procs-1
          call comms_wait(send_handles(proc))
          if(present(vacuum_rep)) then
             call comms_wait(vacuum_send_handles(proc))
          end if
       end do

       call comms_barrier
    end if

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine swx_local_density_fitting

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_swop_quantity_overlap_ppd(swri, &                   ! in/out
       swex_quality, &                                               ! in
       ppd_quantity, n_ppds_in_quantity, ppd_indices_in_quantity, &  ! in
       swop_quantity_overlap, &                                      ! out
       expansion_centres, expansion_atoms, &                         ! in
       num_sws_in_expansion, sw_or_swpot, cell)                      ! in
    !==========================================================================!
    ! This subroutine evaluates the electrostatic integrals that result from   !
    ! a quantity in PPDs interacting with a set of spherical waves or          !
    ! potentials thereof originating on up to 2 atoms. The quantity is usually !
    ! either a product of NGWFs, an atom-pair density, or a sum of weighted    !
    ! products originating on multiple spheres (Q_AaDd).                       !
    !                                                                          !
    ! SWOPs in PPDs are cached in each case.                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2012 basing on a similar subroutine!
    ! by Quintin Hill. Extended in August 2012 by Jacek Dziedzic to work with  !
    ! the second case. Simplified in February 2014 by Jacek Dziedzic to elide  !
    ! the snpos cache, no longer needed.                                       !
    ! Generalised for SW_RI/SW_EX by Jacek Dziedzic in February 2015.          !
    ! Cleaned up by Jacek Dziedzic in June 2017.                               !
    !==========================================================================!

    use constants, only: PI
    use geometry, only: POINT, operator(+)
    use hash_table, only: hash_table_cleanup_at_will
    use rundat, only: pub_hfx_debug
!$  use rundat, only: pub_threads_max
    use simulation_cell, only: CELL_INFO, minimum_image_number_to_displacement,&
         minimum_image_number_to_indices
    use sw_resolution_of_identity, only: ATOM_CENTRE, SW_RI, SW_QUALITY, &
         PBC_IMAGE_INFO, swri_obtain_swops_in_ppd
    use timer, only: timer_clock
    use utils, only: utils_real_to_str, utils_point_to_str, utils_int_to_str, &
         utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(inout)   :: swri
    type(SW_QUALITY), intent(in) :: swex_quality
    integer, intent(in)          :: num_sws_in_expansion
    real(kind=DP), intent(out)   :: swop_quantity_overlap(num_sws_in_expansion)
    real(kind=DP), intent(in)    :: ppd_quantity(:)
    integer, intent(in)          :: n_ppds_in_quantity
    integer, intent(in)          :: ppd_indices_in_quantity(:)
    type(ATOM_CENTRE), intent(in):: expansion_centres(2)
    integer, intent(in)          :: expansion_atoms(2)
    character, intent(in)        :: sw_or_swpot
    type(CELL_INFO), intent(in)  :: cell

    ! jd: Local variables
    real(kind=DP)          :: common_factor
    type(ATOM_CENTRE)      :: src_centre
    real(kind=DP)          :: accum(swex_quality%max_sws_per_centre) ! for cur. centre
    real(kind=DP)          :: swops_in_ppd_workspace(&
         cell%n_pts * swex_quality%max_sws_per_centre)
    type(PBC_IMAGE_INFO), allocatable :: which_images(:,:)
    type(POINT)            :: im1_disp, im2_disp
    integer                :: im1_a1_neighbour, im2_a1_neighbour
    integer                :: im1_a2_neighbour, im2_a2_neighbour
    integer                :: im1_a3_neighbour, im2_a3_neighbour
    integer                :: cur_centre_n
    integer                :: offset_to_src
    integer                :: src_atom
    integer                :: sw_idx
    integer                :: size_of_cached_data
    integer                :: npts
    integer                :: ppd_offset
    integer                :: cur_ppd_n
    integer                :: cur_ppd
    integer                :: num_sws
    integer                :: ipt
    integer                :: ierr

    character(len=*), parameter :: myself = 'swx_swop_quantity_overlap_ppd'

    ! ------------------------------------------------------------------------

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    common_factor = 4.0_DP * PI * cell%weight
    npts = cell%n_pts

    swop_quantity_overlap = 0.0_DP
    offset_to_src = 1

    if(pub_hfx_debug) then
       allocate(which_images(n_ppds_in_quantity*cell%n_pts,2),stat=ierr)
       call utils_alloc_check(myself,'which_images',ierr)
    end if

    all_centres:                                                               &
    do cur_centre_n = 1, 2

       src_atom = expansion_atoms(cur_centre_n)
       src_centre = expansion_centres(cur_centre_n)

       if(src_atom == -1) exit

       ! jd: Recalculate the overlap using (hopefully) cached swops in PPDs.
       accum(:) = 0.0_DP
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(STATIC) &
!$OMP PRIVATE(cur_ppd_n, cur_ppd, size_of_cached_data, num_sws, &
!$OMP      ppd_offset, sw_idx, cur_centre_n) &
!$OMP FIRSTPRIVATE(swops_in_ppd_workspace) &
!$OMP SHARED(n_ppds_in_quantity, swri, src_atom, sw_or_swpot, cell, &
!$OMP      src_centre, npts, ppd_indices_in_quantity, ppd_quantity, &
!$OMP      swex_quality, pub_threads_max, pub_hfx_debug, which_images) &
!$OMP REDUCTION(+:accum)
       ! -----------------------------------------------------------------------
       ! Go over relevant PPDs                                               PPD
       ! -----------------------------------------------------------------------
       all_ppds:                                                               &
       do cur_ppd_n = 1, n_ppds_in_quantity

          cur_ppd = ppd_indices_in_quantity(cur_ppd_n)

          ppd_offset = (cur_ppd_n-1)*npts

          ! jd: Get these SWOPs from hash table or calculate, if absent.
          !     This populates the hash table (from a critical section) and
          !     temporarily locks it against evictions.
          if(.not. pub_hfx_debug) then
             call swri_obtain_swops_in_ppd(swri, &    ! in/out (SWOP caching)
                  swops_in_ppd_workspace, num_sws, &  ! out
                  src_centre, src_atom, cur_ppd, &    ! in
                  sw_or_swpot, cell, swex_quality)    ! in
          else
             call swri_obtain_swops_in_ppd(swri, &    ! in/out (SWOP caching)
                  swops_in_ppd_workspace, num_sws, &  ! out
                  src_centre, src_atom, cur_ppd, &    ! in
                  sw_or_swpot, cell, swex_quality, &  ! in
                  which_images(ppd_offset+1,cur_centre_n)) ! optional out (debug)
          end if

          ! We now have all SWOPs from the current centre for the current PPD
          ! in swops_in_ppd. Accumulate products.
          do sw_idx = 1, swex_quality%num_sws_per_centre
             accum(sw_idx) = accum(sw_idx) + &
                  sum(swops_in_ppd_workspace((sw_idx-1)*npts+1 : &
                  (sw_idx-1)*npts+npts) * &
                  ppd_quantity(ppd_offset+1 : ppd_offset+npts))
          end do

       end do all_ppds
!$OMP END PARALLEL DO

       swop_quantity_overlap(&
            offset_to_src:offset_to_src+swex_quality%num_sws_per_centre-1) = &
            common_factor * accum(1:swex_quality%num_sws_per_centre)

       offset_to_src = offset_to_src + swex_quality%num_sws_per_centre

    end do all_centres

    ! jd: Debug two-centre expansions
    ! NB: This code path is not taken in HFx (which uses swx_expansion_2_fast(),
    !     and thus swx_swop_quantity_overlap_ppd_fast()), it's DMA only.
    if(pub_hfx_debug .and. expansion_atoms(2) /= -1) then

       do ipt = 1, n_ppds_in_quantity * cell%n_pts
          if(which_images(ipt,1)%image /= which_images(ipt,2)%image) then

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
                  f14.6,f14.6,f14.6,a)') '[D] Image 1: '//&
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
             call utils_abort(myself//': [D] Inconsistent images for point '//&
                  trim(utils_int_to_str(ipt))//' in a multiple-PPD quantity.'//&
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
                  'Image 1 at  '//trim(utils_point_to_str(&
                  expansion_centres(1)%incell + im1_disp,'f14.6'))//'.'//CRLF//&
                  'Image 2 at  '//trim(utils_point_to_str(&
                  expansion_centres(2)%incell + im2_disp,'f14.6'))//'.'//CRLF//&
                  'Point at    '//trim(utils_point_to_str(&
                  which_images(ipt,1)%curpoint,'f14.6'))//'.'//CRLF//&
                  'Direct displacement to centre 1: '//trim(utils_point_to_str(&
                  which_images(ipt,1)%disp0,'f14.6'))//'.'//CRLF//&
                  'Direct displacement to centre 2: '//trim(utils_point_to_str(&
                  which_images(ipt,2)%disp0,'f14.6'))//'.'//CRLF//&
                  'Displacement to image 1:         '//trim(utils_point_to_str(&
                  which_images(ipt,1)%disp,'f14.6'))//'.'//CRLF//&
                  'Displacement to image 2:         '//trim(utils_point_to_str(&
                  which_images(ipt,2)%disp,'f14.6'))//'.'//CRLF//&
                  'Direct distance to centre 1: '//trim(utils_real_to_str(&
                  which_images(ipt,1)%rad0))//'.'//CRLF//&
                  'Direct distance to centre 2: '//trim(utils_real_to_str(&
                  which_images(ipt,2)%rad0))//'.'//CRLF//&
                  'Distance to image 1:         '//trim(utils_real_to_str(&
                  which_images(ipt,1)%rad))//'.'//CRLF//&
                  'Distance to image 2:         '//trim(utils_real_to_str(&
                  which_images(ipt,2)%rad))//'.')
#endif
          end if
       end do

    end if

    if(pub_hfx_debug) then
       deallocate(which_images,stat=ierr)
       call utils_dealloc_check(myself,'which_images',ierr)
    end if

    ! Perform deferred cleanups if necessary
    call hash_table_cleanup_at_will(swri%swops_in_ppds_ht)

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine swx_swop_quantity_overlap_ppd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_swop_quantity_overlap_ppd_fast(swri, &              ! in
       ppd_quantity, n_ppds_in_quantity, ppd_indices_in_quantity, &  ! in
       swop_quantity_overlap, &                                      ! out
       expansion_centres, expansion_atoms, &                         ! in
       num_sws_in_expansion, sw_or_swpot, cell)                      ! in
    !==========================================================================!
    ! This subroutine evaluates the electrostatic integrals that result from   !
    ! a quantity in PPDs interacting with a set of spherical waves or          !
    ! potentials thereof originating on up to 2 atoms. The quantity is usually !
    ! either a product of NGWFs, an atom-pair density, or a sum of weighted    !
    ! products originating on multiple spheres (Q_AaDd).                       !
    !                                                                          !
    ! This version is only used in HFx. DMA and polemb use the 'non-fast' ver. !
    !                                                                          !
    ! The differences between the 'non-fast' version and this one are:         !
    ! - We do not modify any HTs here, particularly we do not add to the SWOP  !
    !   cache, only reference it, and we do not do any HT cleanups. This allows!
    !   us to not have any OMP CRICITAL regions and be thread-safe.            !
    ! - swri-to-swex filtering is not supported for efficiency.                !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2012 basing on a similar subroutine!
    ! by Quintin Hill. Extended in August 2012 by Jacek Dziedzic to work with  !
    ! the second case. Simplified in February 2014 by Jacek Dziedzic to elide  !
    ! the snpos cache, no longer needed.                                       !
    ! Generalised for SW_RI/SW_EX by Jacek Dziedzic in February 2015.          !
    ! Cleaned up by Jacek Dziedzic in June 2017. Converted to 'fast' version   !
    ! by Jacek Dziedzic in August 2019.                                        !
    !==========================================================================!

    use constants, only: PI
    use geometry, only: POINT, operator(+)
    use rundat, only: pub_threads_max, pub_hfx_debug
    use simulation_cell, only: CELL_INFO, minimum_image_number_to_displacement,&
         minimum_image_number_to_indices
    use sw_resolution_of_identity, only: ATOM_CENTRE, SW_RI, PBC_IMAGE_INFO, &
         swri_obtain_swops_in_ppd_ptr_fast
    use timer, only: timer_clock
    use utils, only: utils_real_to_str, utils_point_to_str, utils_int_to_str, &
         utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(in)      :: swri
    integer, intent(in)          :: num_sws_in_expansion
    real(kind=DP), intent(out)   :: swop_quantity_overlap(num_sws_in_expansion)
    real(kind=DP), intent(in)    :: ppd_quantity(:)
    integer, intent(in)          :: n_ppds_in_quantity
    integer, intent(in)          :: ppd_indices_in_quantity(:)
    type(ATOM_CENTRE), intent(in):: expansion_centres(2)
    integer, intent(in)          :: expansion_atoms(2)
    character, intent(in)        :: sw_or_swpot
    type(CELL_INFO), intent(in)  :: cell

    ! jd: Local variables
    real(kind=DP)          :: common_factor
    type(ATOM_CENTRE)      :: src_centre
    real(kind=DP)          :: accum(swri%quality%max_sws_per_centre) ! for cur. centre
    real(kind=DP), pointer :: swops_in_ppd_ptr(:)
    real(kind=DP), target  :: swops_in_ppd_workspace(&
         cell%n_pts * swri%quality%num_sws_per_centre, 0:pub_threads_max-1)
    type(PBC_IMAGE_INFO), allocatable :: which_images(:,:)
    type(POINT)            :: im1_disp, im2_disp
    integer                :: im1_a1_neighbour, im2_a1_neighbour
    integer                :: im1_a2_neighbour, im2_a2_neighbour
    integer                :: im1_a3_neighbour, im2_a3_neighbour
    integer                :: cur_centre_n
    integer                :: offset_to_src
    integer                :: src_atom
    integer                :: sw_idx
    integer                :: size_of_cached_data
    integer                :: npts
    integer                :: ppd_offset
    integer                :: cur_ppd_n
    integer                :: cur_ppd
    integer                :: num_sws
!$  integer, external      :: omp_get_thread_num
    integer                :: tid
    integer                :: ipt
    integer                :: ierr

    character(len=*), parameter :: myself = 'swx_swop_quantity_overlap_ppd_fast'

    ! ------------------------------------------------------------------------

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    common_factor = 4.0_DP * PI * cell%weight
    npts = cell%n_pts

    swop_quantity_overlap = 0.0_DP
    offset_to_src = 1

    if(pub_hfx_debug) then
       allocate(which_images(n_ppds_in_quantity*cell%n_pts,2),stat=ierr)
       call utils_alloc_check(myself,'which_images',ierr)
    end if

    all_centres:                                                               &
    do cur_centre_n = 1, 2

       src_atom = expansion_atoms(cur_centre_n)
       src_centre = expansion_centres(cur_centre_n)

       if(src_atom == -1) exit

       ! jd: Recalculate the overlap using (hopefully) cached swops in PPDs.
       accum(:) = 0.0_DP
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(DYNAMIC) &
!$OMP PRIVATE(cur_ppd_n, cur_ppd, size_of_cached_data, num_sws, & ! do NOT change
!$OMP      tid, ppd_offset, sw_idx, swops_in_ppd_ptr)           & ! schedule from
!$OMP SHARED(n_ppds_in_quantity, swri, src_atom, sw_or_swpot,   & ! DYNAMIC
!$OMP      cell, src_centre, npts, ppd_indices_in_quantity, ppd_quantity,      &
!$OMP      pub_threads_max, pub_hfx_debug, swops_in_ppd_workspace, &
!$OMP      which_images, cur_centre_n) &
!$OMP REDUCTION(+:accum)
       ! -----------------------------------------------------------------------
       ! Go over relevant PPDs                                               PPD
       ! -----------------------------------------------------------------------
       all_ppds:                                                               &
       do cur_ppd_n = 1, n_ppds_in_quantity

          cur_ppd = ppd_indices_in_quantity(cur_ppd_n)
          tid = 0
!$        tid = omp_get_thread_num()

          ppd_offset = (cur_ppd_n-1)*npts

          ! jd: Get these SWOPs from hash table or calculate, if absent.
          !     This populates the hash table (from a critical section) and
          !     temporarily locks it against evictions.
          if(.not. pub_hfx_debug) then
             call swri_obtain_swops_in_ppd_ptr_fast(swops_in_ppd_ptr, num_sws, & ! out
                  src_centre, src_atom, cur_ppd, sw_or_swpot, cell, swri, &      ! in
                  swops_in_ppd_workspace(:,tid))
          else
             call swri_obtain_swops_in_ppd_ptr_fast(swops_in_ppd_ptr, num_sws, & ! out
                  src_centre, src_atom, cur_ppd, sw_or_swpot, cell, swri, &      ! in
                  swops_in_ppd_workspace(:,tid), &
                  which_images(ppd_offset+1,cur_centre_n))
          end if

          ! We now have all SWOPs from the current centre for the current PPD
          ! in swops_in_ppd. Accumulate products.
          do sw_idx = 1, swri%quality%num_sws_per_centre
             accum(sw_idx) = accum(sw_idx) + &
                  sum(swops_in_ppd_ptr((sw_idx-1)*npts+1 : &
                  (sw_idx-1)*npts+npts) * &
                  ppd_quantity(ppd_offset+1 : ppd_offset+npts))
          end do

       end do all_ppds
!$OMP END PARALLEL DO

       swop_quantity_overlap(&
            offset_to_src:offset_to_src+swri%quality%num_sws_per_centre-1) = &
            common_factor * accum(1:swri%quality%num_sws_per_centre)

       offset_to_src = offset_to_src + swri%quality%num_sws_per_centre

    end do all_centres

    ! jd: Debug two-centre expansions
    if(pub_hfx_debug .and. expansion_atoms(2) /= -1) then

       do ipt = 1, n_ppds_in_quantity * cell%n_pts
          if(which_images(ipt,1)%image /= which_images(ipt,2)%image) then

             im1_disp = minimum_image_number_to_displacement(&
                  which_images(ipt,1)%image, cell)
             im2_disp = minimum_image_number_to_displacement(&
                  which_images(ipt,2)%image, cell)

             call minimum_image_number_to_indices(which_images(ipt,1)%image, &
                  im1_a1_neighbour, im1_a2_neighbour, im1_a3_neighbour)
             call minimum_image_number_to_indices(which_images(ipt,2)%image, &
                  im2_a1_neighbour, im2_a2_neighbour, im2_a3_neighbour)

#if 1
             write(*,'(a,i3,a,i2,i2,i2,a,i3,a,i2,i2,i2,a,&
                  f14.6,f14.6,f14.6,a,f14.6,f14.6,f14.6,a)') '[1] Image 1: '//&
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
             call utils_abort(myself//': [1] Inconsistent images for point '//&
                  trim(utils_int_to_str(ipt))//' in a multiple-PPD quantity.'//&
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
                  'Image 1 at  '//trim(utils_point_to_str(&
                  expansion_centres(1)%incell + im1_disp,'f14.6'))//'.'//CRLF//&
                  'Image 2 at  '//trim(utils_point_to_str(&
                  expansion_centres(2)%incell + im2_disp,'f14.6'))//'.'//CRLF//&
                  'Point at    '//trim(utils_point_to_str(&
                  which_images(ipt,1)%curpoint,'f14.6'))//'.'//CRLF//&
                  'Direct displacement to centre 1: '//trim(utils_point_to_str(&
                  which_images(ipt,1)%disp0,'f14.6'))//'.'//CRLF//&
                  'Direct displacement to centre 2: '//trim(utils_point_to_str(&
                  which_images(ipt,2)%disp0,'f14.6'))//'.'//CRLF//&
                  'Displacement to image 1:         '//trim(utils_point_to_str(&
                  which_images(ipt,1)%disp,'f14.6'))//'.'//CRLF//&
                  'Displacement to image 2:         '//trim(utils_point_to_str(&
                  which_images(ipt,2)%disp,'f14.6'))//'.'//CRLF//&
                  'Direct distance to centre 1: '//trim(utils_real_to_str(&
                  which_images(ipt,1)%rad0))//'.'//CRLF//&
                  'Direct distance to centre 2: '//trim(utils_real_to_str(&
                  which_images(ipt,2)%rad0))//'.'//CRLF//&
                  'Distance to image 1:         '//trim(utils_real_to_str(&
                  which_images(ipt,1)%rad))//'.'//CRLF//&
                  'Distance to image 2:         '//trim(utils_real_to_str(&
                  which_images(ipt,2)%rad))//'.')

#endif

          end if
       end do

    end if

    if(pub_hfx_debug) then
       deallocate(which_images,stat=ierr)
       call utils_dealloc_check(myself,'which_images',ierr)
    end if

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine swx_swop_quantity_overlap_ppd_fast

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_communicate_coeffs(coeffs_ht, ngwf_basis, nl, s_atoms_nl, &
       swri, par, elements, coeffs_kind, swex_h, mode, ngwf_basis2, ket_basis_selector)
    !==========================================================================!
    ! Communicates expansion coefficients between procs.                       !
    !                                                                          !
    ! During NGWF-pair fitting or atom-pair density fitting every proc expands !
    ! products between NGWFs on local atoms (Bb) and their remote S-neighbours !
    ! (Cc) (multiplied by a kernel element when doing density fitting).        !
    ! During the acting of the potentials in HFx, local NGWFs are Aa, and both !
    ! Bb and Cc are then, in general, remote. B's are then X-neighbours of A,  !
    ! while C are S-neighbours of B. Similarly in DMA coefficients are needed  !
    ! on remote procs, but in an S sparsity pattern, instead of X.             !
    !                                                                          !
    ! This subroutine copies the expansion coefficients to where they will be  !
    ! needed. A general 'nl' neighbour list is passed as an argument to make   !
    ! this more general, e.g. for non-HFx applications. For HFx, simply pass   !
    ! x_atoms_nl for nl. For DMA, simply pass s_atoms_nl for nl.               !
    !                                                                          !
    ! *Regardless* of what is passed for nl, S-neighbours must be figured out  !
    ! between B and C, both remote. Thus the need for s_atoms_nl.              !
    !                                                                          !
    ! The resultant (and input) expansion coefficients are stored in coeffs_ht.!
    ! The coefficients that we started with ('owned' during the density        !
    ! fitting) on a proc and that may no longer be needed during the acting of !
    ! the potentials are *not* evicted from coeffs_ht for simplicity, but also !
    ! for efficiency. If memory consumption becomes an issue, this might need  !
    ! to be rethought.                                                         !
    !                                                                          !
    ! 'coeffs_ht' can contain either of:                                       !
    ! a) spin-agnostic, but formally spin-dependent SW coefficients (in HFx    !
    !    in DMA in mode 'P');                                                  !
    ! b) truly spin-dependent SW coefficients (in DMA in mode 'D');            !
    ! c) truly spin-agnostic SW coefficients (D coefficients in HFx and in DMA !
    !    regardless of mode).                                                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   coeffs_ht (in/out): The coefficients that we own (on input).           !
    !                       The coefficients that we want or own (on output).  !
    !   ngwf_basis (in): Needed to figure out how many NGWFs on each atom.     !
    !   nl (in): Describes who is going to need which coefficients.            !
    !   s_atoms_nl (in): Needed to figure out S-neighbours. Can be identical   !
    !                    to nl.                                                !
    !   swri (in): Needed to skip atoms that are not in the expansion.         !
    !   elements (in): Needed for the same reason.                             !
    !   coeffs_kind (in): Selects the type of coefficients (SW_O, SW_V,        !
    !                     SW_O+SW_D, SW_V+SW_D).                               !
    !   swex_h (in): Used to choose the right set of coefficients.             !
    !   mode (in): 'P' if coeffs resulted from an expansion of NGWF pairs, or  !
    !              'D' if coeffs resulted from an expansion of pair densities. !
    !                                                                          !
    ! *Optional* arguments for mixed NGWF bases only.                          !
    !   rep2 (in): Second rep to use when mixing NGWF bases. Can be intent-in, !
    !              because results (coeffs) are always stored in the first rep.!
    !   ngwf_basis2 (in): Second NGWF basis.                                   !
    !   ket_basis_selector (in): Selects the NGWF basis (1 or 2) in which      !
    !                            - Bb is indexed (ket_basis_selector(1)),      !
    !                            - Cc is indexed (ket_basis_selector(2)).      !
    !                            ket_basis_selector(3) is ignored here, since  !
    !                            we're working with a free coeffs_ht, and not  !
    !                            an NGWF_REP.                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2014.                               !
    ! Generalised for SW_RI/SW_EX by Jacek Dziedzic in February 2015.          !
    ! Generalised to support spin-dependent coefficients by Jacek Dziedzic in  !
    ! June 2017.                                                               !
    ! Extended for mixed NGWF basis sets in June-August 2018 by Jacek Dziedzic.!
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_send, &
         comms_barrier
    use constants, only: VERBOSE, SW_D, SWEX_BB_CC_SYMMETRIC
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE
    use ion, only: ELEMENT
    use neighbour_list, only: NL_NEIGHBOUR_LIST
    use parallel_strategy, only: PARAL_INFO
    use remote, only: COEFFS_REQUEST_TAG, remote_obtain_coeffs, &
         remote_serve_coeffs
    use rundat, only: pub_devel_code, pub_swx_output_detail, pub_num_spins
    use sw_resolution_of_identity, only: SW_RI, swri_get_handle_to
    use timer, only: timer_clock
    use utils, only: utils_devel_code, utils_assert

    implicit none

    ! Arguments
    type(HT_HASH_TABLE), intent(inout), target     :: coeffs_ht
    type(FUNC_BASIS), intent(in)                   :: ngwf_basis
    type(NL_NEIGHBOUR_LIST), intent(in)            :: nl
    type(NL_NEIGHBOUR_LIST), intent(in)            :: s_atoms_nl
    type(SW_RI), intent(in)                        :: swri
    type(PARAL_INFO), intent(in)                   :: par
    type(ELEMENT), intent(in)                      :: elements(par%nat)
    integer, intent(in)                            :: coeffs_kind
    integer, intent(in)                            :: swex_h
    character, intent(in)                          :: mode
    type(FUNC_BASIS), intent(in), optional         :: ngwf_basis2
    integer, intent(in), optional                  :: ket_basis_selector(3)

    ! Local variables
    integer :: local_a
    integer :: global_a
    integer :: b_idx
    integer :: global_b, orig_global_b
    integer :: c_idx
    integer :: global_c, orig_global_c
    logical :: done
    integer :: proc
    logical :: who_is_done(0:pub_total_num_procs-1)
    integer :: send_handle
    integer :: swri_handle
    integer :: loc_num_spins
    integer :: is
    logical :: debug_show_skipped
    logical :: spin_dependent
    character(len=*), parameter :: myself = 'swx_communicate_coeffs'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier
    call timer_clock(myself,1)

    call utils_assert(mode == 'P' .or. mode == 'D', myself//': Invalid mode', &
         mode)

    ! jd: Ensure consistency between optional parameters for mixed NGWF bases
    call utils_assert(present(ngwf_basis2) .eqv. present(ket_basis_selector), &
         myself//': Optional arguments for SWx with mixed NGWF bases must &
         &either all be present or all be absent.')

    swri_handle = swri_get_handle_to(swri%swri_name)

    ! jd: beta_ngwf_basis and gamma_ngwf_basis are not needed here, we simply
    !     pass ngwf_basis2 and the selector down the call chain

    debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
         pub_devel_code, no_bcast=.true., no_warn=.true.)

    if(coeffs_kind >= SW_D) then
       ! jd: D coeffs are spin-agnostic
       spin_dependent = .false.
       loc_num_spins = 1
    else
       spin_dependent = .true.
       ! jd: C coeffs are truly spin-dependent in mode D and
       !     formally spin-dependent in mode P.
       if(mode == 'D') then
          loc_num_spins = pub_num_spins
       else
          loc_num_spins = 1
       end if
    end if

    ! --------------------------------------------------------------------------
    ! jd: Loop over spins (potentially)                                      sss
    ! --------------------------------------------------------------------------
    loop_spins:                                                                &
    do is = 1, loc_num_spins

       who_is_done(:) = .false.

       ! -----------------------------------------------------------------------
       ! jd: Loop over A's -- all atoms on this core                         AAA
       ! -----------------------------------------------------------------------
       loop_A:                                                                 &
       do local_a=1, par%num_atoms_on_proc(pub_my_proc_id)
          global_a = par%first_atom_on_proc(pub_my_proc_id) + local_a-1

          if(pub_swx_output_detail >= VERBOSE) then
             write(stdout,*) ' / A: ', global_a
          end if

          ! --------------------------------------------------------------------
          ! jd: Loop over B's that are nl-neighbours of A                    BBB
          ! --------------------------------------------------------------------
          loop_B:                                                              &
          do b_idx = nl%first_idx(global_a), nl%last_idx(global_a)
             global_b = nl%neighbours(b_idx)
             orig_global_b = par%orig_atom(global_b)

             ! jd: No ignoring of atoms outside SWRI happens here, we only want
             !     to ignore atom pairs where both atoms are outsiders (below).

             ! -----------------------------------------------------------------
             ! jd: Loop over C's that are s-neighbours with B                CCC
             ! -----------------------------------------------------------------
             loop_C:                                                           &
             do c_idx = &
                  s_atoms_nl%first_idx(global_b), s_atoms_nl%last_idx(global_b)
                global_c = s_atoms_nl%neighbours(c_idx)
                orig_global_c = par%orig_atom(global_c)

                ! jd: Ignore atom pairs where neither atom is in this SWRI
                if(.not. elements(orig_global_b)%in_swri(swri_handle) .and. &
                     .not. elements(orig_global_c)%in_swri(swri_handle)) then
                   if(debug_show_skipped) then
                      write(stdout,'(a,i0,a,i0,a,a,a,i0,a,i0,a,a,a)') &
                           'Skipping (swx_communicate_coeffs) B:orig/SFC: ',&
                           orig_global_b, '/', global_b, ' (', &
                           trim(elements(orig_global_b)%species_id),&
                           '): C:orig/SFC', orig_global_c, '/', global_c, ' (',&
                           trim(elements(orig_global_c)%species_id),')'
                   end if
                   cycle
                end if

                ! jd: Obtain the expansion coefficients of the product
                !     of all NGWFs Bb and Cc on atom pair B-C from remote.
                !     Store into coeffs_ht.
                if(.not. spin_dependent) then
                   ! jd: D coeffs are spin-agnostic
                   if(.not. SWEX_BB_CC_SYMMETRIC(swex_h)) then
                      ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                      call remote_obtain_coeffs(coeffs_ht, who_is_done, &
                           global_b, global_c, &
                           coeffs_kind, swex_h, ngwf_basis, elements, swri_handle, par, &
                           ngwf_basis2 = ngwf_basis2, &
                           ket_basis_selector = ket_basis_selector, &
                           overfill_strategy='F')
                   else
                      ! jd: Here there's Bb-Cc symmetry and we only store one triangle
                      call remote_obtain_coeffs(coeffs_ht, who_is_done, &
                           min(global_b,global_c), max(global_b,global_c), &
                           coeffs_kind, swex_h, ngwf_basis, elements, swri_handle, par, &
                           ngwf_basis2 = ngwf_basis2, &
                           ket_basis_selector = ket_basis_selector, &
                           overfill_strategy='F')
                   end if
                else
                   if(.not. SWEX_BB_CC_SYMMETRIC(swex_h)) then
                      ! jd: Here there's no Bb-Cc symmetry and we store all coeffs
                      call remote_obtain_coeffs(coeffs_ht, who_is_done, &
                           global_b, global_c, &
                           coeffs_kind, swex_h, ngwf_basis, elements, swri_handle, par, is, &
                           ngwf_basis2 = ngwf_basis2, &
                           ket_basis_selector = ket_basis_selector, &
                           overfill_strategy='F')
                   else
                      ! jd: Here there's Bb-Cc symmetry and we only store one triangle
                      call remote_obtain_coeffs(coeffs_ht, who_is_done, &
                           min(global_b,global_c), max(global_b,global_c), &
                           coeffs_kind, swex_h, ngwf_basis, elements, swri_handle, par, is, &
                           ngwf_basis2 = ngwf_basis2, &
                           ket_basis_selector = ket_basis_selector, &
                           overfill_strategy='F')
                   end if
                end if

             end do loop_C

          end do loop_B

       end do loop_A

       ! jd: After we're done, first notify everyone, then keep serving
       !     other procs with coeffs, until everyone is done
       do proc = 0, pub_total_num_procs-1
         ! jd: NB, two notifications are needed for coeffs
          call comms_send(proc, -1, 1, tag = COEFFS_REQUEST_TAG, &
               return_handle = send_handle, add_to_stack = .false.)
          call comms_send(proc, -1, 1, tag = COEFFS_REQUEST_TAG, &
               return_handle = send_handle, add_to_stack = .false.)
       end do

       done = .false.
       if(pub_total_num_procs > 1) then
          do while(.not. done)
             if(spin_dependent) then
                call remote_serve_coeffs(done, who_is_done, ngwf_basis, &
                     coeffs_ht, coeffs_kind, swex_h, is, &
                     ngwf_basis2 = ngwf_basis2, &
                     ket_basis_selector = ket_basis_selector)
             else
                call remote_serve_coeffs(done, who_is_done, ngwf_basis, &
                     coeffs_ht, coeffs_kind, swex_h, &
                     ngwf_basis2 = ngwf_basis2, &
                     ket_basis_selector = ket_basis_selector)
             end if
          end do
       end if

    end do loop_spins

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine swx_communicate_coeffs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_ngwf_product(prod_in_fftbox_dbl_on_bb, &
       fftbox_ngwf_bb, fftbox_ngwf_cc, fftbox)
    !==========================================================================!
    ! Calculates a product of two NGWFs supplied in FFTboxes, interpolating    !
    ! first to the double grid. The result is returned in a double-grid FFtbox !
    ! on the first of the NGWFs.                                               !
    !                                                                          !
    ! The result is allocated here, it's up to the caller to free it later.    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   prod_in_fftbox_dbl_on_bb (out): The result is returned here.           !
    !   fftbox_ngwf_bb (in): The first of the NGWFs to multiply.               !
    !   fftbox_ngwf_cc (in): The second of the NGWFs to multiply.              !
    !   fftbox (in): The definition of an FFTbox.                              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2017                                !
    !==========================================================================!

    use datatypes, only: FFTBOX_DATA, data_fftbox_alloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_interpolate_product
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FFTBOX_DATA), intent(out) :: prod_in_fftbox_dbl_on_bb
    type(FFTBOX_DATA), intent(in)  :: fftbox_ngwf_bb
    type(FFTBOX_DATA), intent(in)  :: fftbox_ngwf_cc
    type(FFTBOX_INFO), intent(in)  :: fftbox

    ! Local Variables
    complex(kind=DP), allocatable, dimension(:,:,:) :: coarse_work, fine_work
    real(kind=DP), allocatable, dimension(:,:,:)    :: prod_in_fftbox_on_bb_real
    integer :: ierr
    character(len=*), parameter :: myself = 'swx_ngwf_product'

    ! -------------------------------------------------------------------------

    ! jd: Find interpolated product of Bb and Cc in CC's FFT box
    allocate(coarse_work(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3),stat=ierr)
    call utils_alloc_check(myself,'coarse_work',ierr)
    allocate(fine_work(fftbox%total_ld1*2, &
         fftbox%total_ld2*2,fftbox%total_pt3*2),stat=ierr)
    call utils_alloc_check(myself,'fine_work',ierr)
    allocate(prod_in_fftbox_on_bb_real(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl),stat=ierr)
    call utils_alloc_check(myself,'prod_in_fftbox_on_bb_real',ierr)

    call fourier_interpolate_product(coarse_work, fine_work, &
         fftbox_ngwf_bb%d, fftbox_ngwf_cc%d, prod_in_fftbox_on_bb_real)

    call data_fftbox_alloc(prod_in_fftbox_dbl_on_bb, fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl, fftbox%total_pt3_dbl, .false.)
    prod_in_fftbox_dbl_on_bb%d(:,:,:) = prod_in_fftbox_on_bb_real

    deallocate(prod_in_fftbox_on_bb_real, stat=ierr)
    call utils_dealloc_check(myself,'prod_in_fftbox_on_bb_real',ierr)
    deallocate(fine_work,stat=ierr)
    call utils_dealloc_check(myself,'fine_work',ierr)
    deallocate(coarse_work, stat=ierr)
    call utils_dealloc_check(myself,'coarse_work',ierr)

  end subroutine swx_ngwf_product

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swx_swop_quantity_overlap_fftbox_dbl(swri, &               ! in/out
       swex_quality, &                                                  ! in
       fftbox_start_in_cell_vec, fftbox_dbl_quantity, fftbox, &         ! in
       swop_quantity_overlap, &                                         ! out
       expansion_centres, expansion_atoms, &                            ! in
       num_sws_in_expansion, sw_or_swpot, cell)                         ! in
    !==========================================================================!
    ! Calculates the overlap of SWs or SWpots originating on up to two atoms   !
    ! with a quantity specified in a double-grid FFTbox.                       !
    ! This is a first, ugly attempt at an alternative to                       !
    ! swx_swop_quantity_overlap_ppd, needed for double-grid SW expansions.     !
    !                                                                          !
    ! SWOPs at points are cached (inefficiently). Stashes are used to improve  !
    ! the performance of caching.                                              !
    !--------------------------------------------------------------------------!
    ! This routine will likely be retired, as double-grid SW expansions do not !
    ! seem to offer much benefit, only headaches.                              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2017.                               !
    !==========================================================================!

    use constants, only: PI
    use datatypes, only: FFTBOX_DATA
    use fft_box, only: FFTBOX_INFO
    use geometry, only: POINT, operator(*), operator(-)
! jme: workaround for cray compiler bug 843178
#if _CRAYFTN
    use geometry, only: add_points
#else
    use geometry, only: operator(+)
#endif
    use hash_table, only: hash_table_stash_prepare, hash_table_stash_commit, &
         hash_table_stash_free
    use rundat, only: pub_threads_max
    use simulation_cell, only: CELL_INFO
    use sw_resolution_of_identity, only: ATOM_CENTRE, SW_RI, SW_QUALITY, &
         swri_obtain_swops_at_point
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(inout)    :: swri
    type(SW_QUALITY), intent(in)  :: swex_quality
    integer, intent(in)           :: num_sws_in_expansion
    real(kind=DP), intent(out)    :: swop_quantity_overlap(num_sws_in_expansion)
    type(POINT), intent(in)       :: fftbox_start_in_cell_vec
    type(FFTBOX_DATA), intent(in) :: fftbox_dbl_quantity
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(ATOM_CENTRE), intent(in) :: expansion_centres(2)
    integer, intent(in)           :: expansion_atoms(2)
    character, intent(in)         :: sw_or_swpot
    type(CELL_INFO), intent(in)   :: cell

    ! jd: Local variables
    real(kind=DP)     :: common_factor
    type(ATOM_CENTRE) :: src_centre
    type(ATOM_CENTRE) :: other_centre
    real(kind=DP)     :: accum(swex_quality%num_sws_per_centre) ! for cur. centre
    real(kind=DP)     :: buffer(swex_quality%num_sws_per_centre)
    integer           :: cur_centre_n
    integer           :: offset_to_src
    integer           :: src_atom
    integer           :: npts_dbl
    type(POINT)       :: a1_dbl, a2_dbl, a3_dbl
    integer           :: d1idx, d2idx, d3idx, didx
    integer           :: points_used
    type(POINT)       :: curpoint
    real(kind=DP)     :: rad, rad_other
    type(POINT)       :: disp, disp_other
    real(kind=DP), allocatable :: swops_stash(:)
    integer, parameter :: max_stash_size = 1024 ! == 1GiB, generalizeme, perthread
    character(len=*), parameter :: myself = 'swx_swop_quantity_overlap_fftbox_dbl'

    ! ------------------------------------------------------------------------

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    call utils_assert(.not. fftbox_dbl_quantity%iscmplx, 'Error in '//myself//&
         ': not ready for complex NGWFs yet.')

    common_factor = 4.0_DP * PI * fftbox%weight * 1.0_DP/8.0_DP ! (dbl grid)
    npts_dbl = fftbox%total_pt1_dbl * fftbox%total_pt2_dbl * fftbox%total_pt3_dbl

    a1_dbl = 0.5_DP * fftbox%d1 * fftbox%a1_unit
    a2_dbl = 0.5_DP * fftbox%d2 * fftbox%a2_unit
    a3_dbl = 0.5_DP * fftbox%d3 * fftbox%a3_unit

    swop_quantity_overlap = 0.0_DP
    offset_to_src = 1
    points_used = 0

    ! jd: Loop over the two SW centres
    all_centres: do cur_centre_n = 1, 2

       src_atom = expansion_atoms(cur_centre_n)
       src_centre = expansion_centres(cur_centre_n)

       if(src_atom == -1) exit

       other_centre = expansion_centres(3-cur_centre_n)

       accum(:) = 0.0_DP

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP FIRSTPRIVATE(swops_stash) &
!$OMP PRIVATE(didx, curpoint, buffer, d1idx, d2idx, d3idx, rad, disp, &
!$OMP      rad_other, disp_other) &
!$OMP SHARED(fftbox_start_in_cell_vec, npts_dbl, a1_dbl, a2_dbl, a3_dbl, &
!$OMP      swex_quality, fftbox, fftbox_dbl_quantity, src_centre, src_atom, &
!$OMP      sw_or_swpot, swri, cell, other_centre) &
!$OMP REDUCTION(+:accum, points_used)

       ! jd: Have each thread prepare its own stash
       call hash_table_stash_prepare(swops_stash, max_stash_size, &
            overwrite = .false., overfill_strategy = 'N')

       call timer_clock('swx_swop_quantity...main_loop',1)
!$OMP DO SCHEDULE(DYNAMIC)
       do didx = 1, npts_dbl

          d1idx = modulo(didx-1,fftbox%total_pt1_dbl) + 1
          d2idx = modulo((didx-d1idx)/fftbox%total_pt1_dbl,fftbox%total_pt2_dbl) + 1
          d3idx = (didx-(d2idx-1)*fftbox%total_pt1_dbl-d1idx) / &
               (fftbox%total_pt1_dbl*fftbox%total_pt2_dbl) + 1

          ! jd: === Truncation based on value ===
          ! jd: Simple shaving to zero outside the NGWF radius is not stable.
          if(abs(fftbox_dbl_quantity%d(d1idx,d2idx,d3idx)) < 1D-25) cycle

! jme: workaround for cray compiler bug 843178
#if _CRAYFTN
          curpoint = add_points( &
               add_points(fftbox_start_in_cell_vec, (d1idx-1) * a1_dbl), &
               add_points((d2idx-1) * a2_dbl, (d3idx-1) * a3_dbl) )
#else
          curpoint = fftbox_start_in_cell_vec + &
               (d1idx-1) * a1_dbl + (d2idx-1) * a2_dbl + (d3idx-1) * a3_dbl
#endif

          points_used = points_used + 1

          ! jd: This performs read accesses to swops_at_points_ht and
          !     writes to the stash.
          call swri_obtain_swops_at_point(swri, cell, buffer, curpoint, &
               src_centre, src_atom, sw_or_swpot, swex_quality, swops_stash)

          accum(:) = accum(:) + &
               buffer(:) * fftbox_dbl_quantity%d(d1idx,d2idx,d3idx)

       end do
!$OMP END DO
       call timer_clock('swx_swop_quantity...main_loop',2)

       ! jd: Merge each thread's stash into the shared hash table now that
       !     it's not being accessed anymore. Destroy stashes.
       call timer_clock('swx_swop_quantity...stash_commit',1)
!$OMP CRITICAL(SEC_SWOPS_AT_POINTS_HT_ACCESS)
       call hash_table_stash_commit(swops_stash, swri%swops_at_points_ht)
!$OMP END CRITICAL(SEC_SWOPS_AT_POINTS_HT_ACCESS)
       call timer_clock('swx_swop_quantity...stash_commit',2)

       call hash_table_stash_free(swops_stash)

!$OMP END PARALLEL

       swop_quantity_overlap(&
            offset_to_src:offset_to_src+swex_quality%num_sws_per_centre-1) = &
            common_factor * accum(1:swex_quality%num_sws_per_centre)

       offset_to_src = offset_to_src + swex_quality%num_sws_per_centre

    end do all_centres

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine swx_swop_quantity_overlap_fftbox_dbl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module sw_expansion
