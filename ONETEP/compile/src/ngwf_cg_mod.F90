! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!            NGWF Conjugate Gradients optimisation module        !
!                                                                !
! This module performs the optimisation of the electronic energy !
! with respect to the NGWFs by Conjugate Gradients Minimisation. !
!----------------------------------------------------------------!
! This module was created by Nicholas Hine on 28/04/2010 out of  !
! part of the previous version of electronic_mod.                !
! Originally written by Chris-Kriton Skylaris in 2000.           !
! Subsequent modifications by Chris-Kriton Skylaris,             !
! Arash A. Mostofi, Peter D. Haynes and Nicholas D.M. Hine.      !
! Modified by Nicholas Hine to use the NGWF_REP container for    !
! NGWF-related data and matrices, October 2010.                  !
! Further modifications by Laura Ratcliff, David O'Regan and     !
! Alvaro Ruiz Serrano.                                           !
!================================================================!


module ngwf_cg

  implicit none

  private

  public :: ngwf_cg_optimise
  public :: ngwf_cg_emft_reset

  ! jcap: module-wide variables required to ensure correct resetting
  ! of NGWF optimisation after an EMFT follow calculation
  logical, save :: orig_pub_emft
  logical, save :: orig_pub_build_bo
  integer, save :: orig_pub_maxit_ngwf_cg

contains

  subroutine ngwf_cg_optimise(total_energy, converged, out_of_runtime, & ! out
       ham, denskern, edft, rep, ngwf_nonsc_forces, lhxc_fine, &         ! inout
       ngwf_basis, proj_basis, nl_projectors, hub_proj_basis, hub, mdl,& ! in
       hfxstate, &                                                       ! in
       val_rep, val_ngwf_basis, val_dkn, val_ham, kpt, & ! lr408: optional cond args
       dfdtau_fine) ! JCW: optional tau-dependent XC functional args

    !==========================================================================!
    ! This subroutine minimises the total energy with respect to the NGWF      !
    ! expansion coefficients and the density kernel elements subject to the    !
    ! constraint that the density matrix is idempotent and integrates to       !
    ! the correct number of electrons.                                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! total_energy    (output): the total energy                               !
    ! converged       (output): whether the total energy was converged         !
    ! out_of_runtime  (output): whether the calculation aborted as out of time !
    ! rep              (inout): NGWF Representation (functions and matrices)   !
    ! denskern         (inout): density kernel                                 !
    ! ngwf_basis       (input): Function basis type describing the NGWFs       !
    ! proj_basis       (input): Function basis type describing the nonlocal    !
    !                           pseudopotential projectors                     !
    ! hub_proj_basis   (input): Function basis type describing the Hubbard     !
    !                           projectors                                     !
    ! lhxc_fine        (inout): Local-Hartree-Exhange-Correlation potential    !
    ! ngwf_nonsc_forces(inout): Outer loop non self-consistent force correction!
    ! For Conduction NGWF optimisation only:                                   !
    ! val_rep          (input): Valence NGWF representation (optional)         !
    ! val_ngwf_basis   (input): Valence NGWF basis (optional)                  !
    ! val_dkn          (input): Valence density kernel (optional)              !
    ! dfdtau_fine     (output): Gradient of XC energy per unit volume wrt      !
    !                        KE density energy (optional)                      !
    !--------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in June 2001.                           !
    ! Modified by Arash Mostofi.                                               !
    ! Modified by Chris-Kriton Skylaris in June 2003.                          !
    ! Modified by Chris-Kriton Skylaris October-December 2003 so that          !
    ! it runs in parallel.                                                     !
    ! Modified by Chris-Kriton Skylaris on 22/02/2004 to reduce                !
    ! memory requirements.                                                     !
    ! Modified by Peter Haynes on 1/07/2004 to use fourier parallelisation     !
    ! Fix of bug in calculation summary by Chris-Kriton Skylaris on 16/07/2004.!
    ! Addition of qc printout and improvements in tests for energy convergence !
    ! by Chris-Kriton Skylaris on 12/03/2005.                                  !
    ! Modification to increase number of lnv iterations when NGWF convergence  !
    ! stagnates by Chris-Kriton Skylaris on 22/03/2005.                        !
    ! Modified for spin polarisation by Peter Haynes, July 2006                !
    ! Modified for Nonlinear Core Corrections by Nicholas Hine, January 2009   !
    ! Modified by David O'Regan for DFT+U, April 2009                          !
    ! Adapted for SPAM3 and function basis by Nicholas Hine, May-July 2009.    !
    ! Renamed from electronic_energy_minimise_tc to ngwf_cg_optimise and moved !
    ! to new module by Nicholas Hine on 28/04/2010.                            !
    ! Modifications for NGWF_REP and NGWF_HAM types by Nicholas Hine in        !
    ! October 2010.                                                            !
    ! Modified by Laura Ratcliff for conduction calculations, Oct 2010.        !
    ! Modified for Christoffel Kernel updates by David O'Regan and Nicholas    !
    ! Hine in November 2010.                                                   !
    ! New convergence criteria by Nicholas Hine, 11/11/2010.                   !
    ! Non self-consistent forces correction by Alvaro Ruiz Serrano, 19/11/2010.!
    ! Return flag indicating convergence added by Nicholas Hine 20/04/2011.    !
    ! Return flag if insufficient runtime left for NGWF CGs. Rob Bell, Aug 2013!
    ! Modified by Chris-Kriton Skylaris to make sure that NGWFs are always     !
    ! consistent with the kernel 22/7/2014.                                    !
    ! Modified for embedding by Robert Charlton, August 2017.                  !
    !==========================================================================!

    use cdft_intermediate_cg, only: cdft_intermediate_u_cg_optimise
    use comms, only: pub_on_root, pub_root_proc_id, comms_bcast, comms_reduce
    use constants, only: DP, verbose, max_spins, normal, stdout, brief, LONG, &
         EDA_ISOLATED, EDA_POLSIMUL, EDA_POLSIMUL_OFF, &
         EDA_POLFRAGLOC_DEVEL, EDA_POLFRAGLOC_DEVEL_OFF, &
         EDA_CTFRAGLOC_DEVEL, EDA_CTFRAGLOC_DEVEL_OFF, k_B
    use datatypes, only: COEF, FUNCTIONS,  &
         data_functions_alloc, data_functions_dealloc, data_functions_dot, &
         data_functions_copy, data_functions_axpy, data_set_to_zero, &
         data_functions_scale, data_size
    use dense, only: DEM, dense_create, dense_destroy, dense_copy
    use electronic, only: electronic_energy, electronic_lagrangian
    use electronic_history, only: elec_history_store_ngwfs, &
         elec_history_store_dkn, elec_history_backup_history, &
         elec_history_update_dkn, elec_history_compute_temp
    use ensemble_dft_type, only: EDFT_MODEL
    use fragment_data, only: fragment_data_add_fragment, eda_dfdtau_fine
    use fragment_scfmi, only: denskern_R
    use function_basis, only: FUNC_BASIS, function_basis_est_num_psincs
    ! agrecokpt: needed for kpt argument in projectors_func_ovlp_box
    use geometry, only: POINT
    use hamiltonian, only: hamiltonian_build_matrix, &
         hamiltonian_dens_indep_matrices, hamiltonian_energy_components, &
         hamiltonian_dens_dep_nonsc
    use hf_exchange, only: HFX_STATE, hf_exchange_dkn_indep_stage
    use hubbard_build, only: HUBBARD_MODEL, hubbard_energy_info, &
         cdft_energy_info, dft_nu_energy_info
    use kernel, only: DKERN, kernel_create, kernel_destroy, kernel_validate_ks,&
         kernel_workspace_create, kernel_workspace_destroy
    use mermin, only: mermin_switch_mu
    use model_type, only: MODEL
    use ngwf_gradient, only: ngwf_gradient_lnv, &
         ngwf_gradient_paw_precond_init, ngwf_gradient_paw_precond_exit
    use ngwf_representation, only: NGWF_REP, NGWF_HAM, &
         ngwf_rep_register_change, ngwf_ham_register_change
    use nonsc_forces, only: nonsc_forces_ngwfs_calc
    use parallel_strategy, only: PARAL_INFO
    use polarisable_embedding, only: polarisable_embedding_expand_ngwf_pairs
    use projectors, only: PROJECTOR_SET
    use restart, only: restart_kernel_write, restart_ngwfs_tightbox_output, &
         restart_hamiltonian_write, restart_overlap_write, &
         restart_sph_waves_output
    use rundat, only: pub_kernel_update, pub_write_tightbox_ngwfs, pub_ngwf_cg_type, &
         pub_output_detail, pub_nnho, pub_precond_real, pub_precond_recip, &
         pub_maxit_ngwf_cg, minit_lnv, ngwf_threshold_orig, pub_elec_cg_max, &
         pub_write_ngwf_plot, pub_cube_format, pub_dx_format, pub_dftb, &
         pub_grd_format, pub_hubbard, pub_hub_on_the_fly, &
         pub_print_qc, pub_k_zero, lnv_threshold_orig, &
         pub_write_sw_ngwfs, cond_firstit_lnv, pub_firstit_lnv, &
         pub_task, pub_ngwf_cg_max_step, pub_write_denskern, pub_write_hamiltonian, &
         pub_cdft, pub_maxit_cdft_u_cg,pub_perturbative_soc, &
         pub_cdft_tight, pub_cdft_cg_threshold, pub_cdft_max_grad, pub_edft, &
         pub_write_overlap, pub_write_converged_dk_ngwfs, &
         pub_aug, pub_cond_calculate, pub_kernel_christoffel_update, &
         pub_elec_force_tol, pub_nonsc_forces, &
         pub_devel_code, pub_write_ngwf_radial, &
         pub_ngwf_cg_rotate, cond_minit_lnv, pub_forces_needed,&
         pub_debug_on_root, pub_num_spins, pub_dft_nu, pub_foe, &
         pub_confined_ngwfs, pub_maxit_ngwf_cg_confined, &
         pub_confined_ngwfs_barrier, md_write_history, md_iter_global, &
         mix_dkn_type, pub_num_kpoints, PUB_1K, &
         pub_kpoint_method, pub_pol_emb_qmstar, &
         pub_eda, pub_eda_mode, pub_eda_scfmi, pub_frag_counter, &
         pub_dmft_spoil_kernel, pub_dmft_fully_sc, pub_dmft_fully_sc_h, &
         pub_xc_ke_density_required, mix_dkn_init_num, md_aux_rep, pub_use_hfx,&
         pub_energy_components_interval, pub_pol_emb_vacuum_qmstar, &
         pub_freeze_switch_steps, pub_quit_region, pub_do_fandt, pub_emft, &
         pub_emft_follow, pub_emft_lnv_only, pub_freeze_envir_ngwfs, &
         pub_active_region, pub_emft_lnv_steps, pub_build_bo, &
         pub_mermin, pub_mermin_smearing_width, pub_spin_fac, pub_mermin_cheb, &
         pub_mermin_temp, pub_firstit_mermin, maxit_mermin, pub_mermin_mu_sq, &
         mermin_threshold_orig
    use services, only: services_flush, services_cubic_fit_minimum, &
        services_line_search_parabola
    use smearing_operator, only: smearing_matrix, lowdin_transformation, &
        calculate_fermi_entropy_mermin
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_copy, &
         sparse_product
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_copy, sparse_embed_product, sparse_embed_axpy, &
         sparse_embed_array_create, sparse_embed_array_destroy, &
         sparse_embed_array_copy, SPAM3_EMBED_ARRAY, sparse_embed_trace, &
         sparse_embed_scale
    use timer, only: timer_clock, timer_check_iteration_time
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_devel_code, &
         utils_assert, utils_abort, utils_banner, utils_flushed_string_output
    use visual, only: visual_ngwfs, visual_ngwfs_radial

    implicit none

    ! Arguments
    type(NGWF_REP), intent(inout) :: rep
    type(FUNC_BASIS), intent(in) :: proj_basis(rep%nsub)
    type(FUNC_BASIS), intent(inout) :: hub_proj_basis(rep%nsub)
    type(MODEL), intent(inout), target :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    type(FUNC_BASIS), intent(inout) :: ngwf_basis(rep%nsub)
    type(DKERN), intent(inout)   :: denskern
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(rep%nsub)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(EDFT_MODEL), intent(inout) :: edft
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, rep%nsub)
    real(kind=DP), intent(out)   :: total_energy
    type(NGWF_HAM), intent(inout) :: ham
    ! ars: non self-consistent forces correction due to NGWF loop
    real(kind=DP), intent(inout) :: ngwf_nonsc_forces(:,:)
    logical, intent(out) :: converged  ! pdh: convergence flag
    logical, intent(out) :: out_of_runtime  ! rab207
    ! JCW: Optional gradient of XC energy per unit volume wrt kinetic energy
    ! JCW: density.
    real(kind=DP), optional, intent(out) :: dfdtau_fine(:,:,:,:)
    ! JCW: In general, we have
    ! JCW:    dfdtau_fine(mdl%fine_grid%ld1,&
    ! JCW:    mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    ! JCW: or
    ! JCW:    dfdtau_fine(1,1,1,1)
    ! JCW: when it is unneeded.

    ! lr408: Optional arguments needed for conduction calculation
    type(NGWF_REP), optional, intent(in) :: val_rep
    type(NGWF_HAM), optional, intent(in) :: val_ham
    type(FUNC_BASIS), optional, intent(in) :: val_ngwf_basis(:) ! size(val_rep%nsub)
    type(SPAM3_EMBED_ARRAY), optional, intent(in)      :: val_dkn

    ! agrecokpt: optional argument to specify non-Gamma k-point
    real(kind=DP), optional, intent(in) :: kpt(3)

    ! Local Variables
    type(DKERN) :: pur_denskern
    type(SPAM3_EMBED_ARRAY) :: start_denskern
    type(SPAM3_EMBED) :: start_sp_overlap
    type(FUNCTIONS) :: cov_grad_on_grid(mdl%nsub)
    type(FUNCTIONS) :: prev_contra_grad_on_grid(mdl%nsub)
    type(FUNCTIONS) :: contra_grad_on_grid(mdl%nsub)
    type(FUNCTIONS) :: direction_on_grid(mdl%nsub)
    type(FUNCTIONS) :: start_ngwfs_on_grid(mdl%nsub)
    type(FUNCTIONS) :: prev_direction_on_grid(mdl%nsub)
    type(smearing_matrix) :: tmp_overlap, lwd_overlap
    real(kind=DP), allocatable :: total_forces(:,:,:)
    real(kind=DP) :: last_n_energies(3) ! space to store the 3 most recent energies
    real(kind=DP) :: mu(max_spins)
    real(kind=DP) :: rms_gradient
    real(kind=DP) :: max_gradient
    real(kind=DP) :: previous_rms_gradient ! RMS NGWF grad of previous iteration
    real(kind=DP) :: line_search_coeff
    real(kind=DP) :: lnv_threshold
    real(kind=DP) :: mermin_threshold
    real(kind=DP) :: ngwf_threshold
    real(kind=DP) :: F0,F1,F2,trial_length
    type(COEF) :: G_init
    type(COEF) :: cg_coeff   !jmecmplx
    real(kind=DP) :: quadratic_coeff,cubic_coeff,rejected_quadratic_coeff
    real(kind=DP) :: predicted_functional
    type(COEF) :: previous_g_dot_g
    type(COEF) :: current_g_dot_g
    type(COEF) :: previous_dir_dot_g
    character(len=88), allocatable, dimension(:) :: summary_lines
    integer :: cvstat            ! =0 if lnv_threshold, /=0 otherwise
    integer :: iteration         ! current iteration
    integer :: cg_count          ! current number of steps since CG reset
    integer(kind=LONG) :: est_num_psincs ! estimate number of psincs in all NGWF spheres
    integer :: current_maxit_lnv ! number of lnv iterations to do for current NGWF step
    integer :: current_maxit_mermin ! number of mermin iterations to do for current NGWF step
    integer :: is         ! pdh: spin loop counter
    integer :: ierr       ! error flag
    integer :: minit      ! qoh: Minimum number of CG iterations
    logical :: trial2     ! ndmh: flag to perform second trial step
    logical :: retrial1   ! ndmh: flag to perform repeat first trial step
    logical :: reversing  ! ndmh: line search is going uphill
    logical :: line_search_success ! ndmh: line search fit success flag
    logical :: check_conjugacy     ! ndmh: flag for doing cg conjugacy check
    logical :: updated_shift ! lr408: Flag needed for conduction calculations
    logical :: cdft_converged, cdft_coarse_converged ! gibo: cDFT-convergence flags
    type(NGWF_REP) :: vacuum_rep ! jd: For polarisable embedding
    ! agrecocmplx
    logical :: loc_cmplx
    real(kind=DP) :: G_init_real
    ! agrecokpt
    type(POINT) :: loc_kpt
    real(kind=DP) :: loc_kpt_cart(3)
    ! rc2013: F+T convergence checks
    logical :: ft_all_converged
    logical, allocatable :: ft_converged(:)

    logical :: nonscenergy ! JCW: Logical value associated with devel_code
                           ! NGWFCG:NONSCENERGY=T/F:NGWFCG

    ! ndmh: FD variables
    integer :: ifd
    ! jd: More closely-spaced points where it matters.
    integer, parameter :: nfd=31
    real(kind=DP), parameter :: logd=1.58489_DP
    real(kind=DP), parameter :: start = 1D-6
    real(kind=DP) :: fd_trial_step(nfd)=(/&
         start*logd**0, start*logd**1, start*logd**2, start*logd**3, &
         start*logd**4, start*logd**5, start*logd**6, start*logd**7, &
         start*logd**8, start*logd**9, start*logd**10, start*logd**11, &
         start*logd**12, start*logd**13, start*logd**14, start*logd**15, &
         start*logd**16, start*logd**17, start*logd**18, start*logd**19, &
         start*logd**20, start*logd**21, start*logd**22, start*logd**23, &
         start*logd**24, start*logd**25, start*logd**26, start*logd**27, &
         start*logd**28, start*logd**29, start*logd**30/)
    real(kind=DP) :: FFD(nfd)

    !ep: FAST DEBUG MERMIN
    !real(kind=DP) :: fd_trial_step(nfd)=(/0.000001_DP, 0.00001_DP,0.0001_DP/)
    !ep: FAST DEBUG MERMIN

    !ep: LONG STEP DEBUG FOR LOW MERMIN CHEB EXPANSION
    !real(kind=DP) :: fd_trial_step(nfd)=(/0.000001_DP, 0.00001_DP,0.0001_DP, &
    !     0.001_DP,0.01_DP,0.1_DP,0.2_DP,0.3_DP,0.4_DP,0.5_DP, 0.6_DP,  &
    !     0.7_DP,0.8_DP, 0.9_DP, 1.0_DP, 1.5_DP, 2.0_DP, 2.5_DP, 3.0_DP, 3.5_DP,
    !     4.0_DP, &
    !     4.5_DP, 5.0_DP, 6.0_DP, 7.0_DP, 8.0_DP, 9.0_DP, 10.0_DP, 20.0_DP,
    !     30.0_DP, &
    !     40.0_DP, 50.0_DP /)
    !ep: LONG STEP DEBUG FOR LOW MERMIN CHEB EXPANSION


    ! cks: FOR PLOTTING FUNCTIONS
    character(len=64) :: fun_name, sub_fun_name, iter_number

    ! ars: force the actual evaluation of the functional after polynomial fit
    logical :: evaluate_functional, force_cubic_step
    real(kind=DP) :: Ffinal

    ! ars: buffer to store edft MOs
    type(DEM), allocatable :: start_mo(:)
    ! ars: buffer to store the initial inv_overlap
    type(SPAM3_EMBED) :: start_inv_overlap
    ! ars: inner product of the gradient
    type(COEF) :: inner_grad
    ! ars: flags for devel_code
    logical :: rotate_prev_dir, calculate_inner_grad, orthogonalise_dir, old_cg_logic

    ! ndmh: flag for whether current gradient is consistent with kernel
    logical :: gradient_valid

    ! ja531 --> When not doing FOE, enable this to still compute the FOE rotation, for testing and debug only.
    logical :: foe_test_rot=.false.

    ! vv: flag for whether the current kernel workspace needs to be allocated
    logical :: deallocate_workspace
    type(SPAM3_EMBED) :: ks_half, overlap_half
    type(SPAM3_EMBED_ARRAY) :: Qmat
    real(kind=DP) :: normS
    ! rc2013: parallel strategy
    type(PARAL_INFO), pointer :: par
    integer :: mrows, ncols
    ! rc2013: temporary structure to hold DKERN%kern%m(:,PUB_1K) info
    integer        :: ireg, ift, isub
    ! jcap: temporary variables for calculation of nonsc forces
    integer :: iat
    real(kind=DP),allocatable :: sub_nonsc_forces(:,:)
    character(len=3) :: suffix
    ! rc2013: subsystem counter
    character(len=1) :: isub_str
    character(len=30) :: file_suffix
    ! rc2013: have we switched to the active functional yet?
    logical :: fnl_switched
    !ep: mermin variables
    real(kind=DP) :: smermin(pub_num_spins)
    integer :: ik, jsub
    real(kind=DP) :: Ftmp
    real(kind=DP) :: muext(max_spins,pub_num_kpoints)
    type(COEF) :: pre_dir_dot_g
    type(COEF) :: now_dir_dot_g
    !ep: mermin variables


    if (pub_debug_on_root) write(stdout,'(a)')'DEBUG: Entering ngwf_cg_optimise'

    ! Flush output
    call services_flush

    ! Start timer
    call timer_clock('ngwf_cg_optimise',1)

    if (.not.pub_xc_ke_density_required) then
       ! JCW: Non-tau-dependent XC functional
       ! JCW: dfdtau_fine is either not present, or has a length of 1
       if ( present(dfdtau_fine ) ) then
          call utils_assert( &
               size(dfdtau_fine) == 1,&
               "Error in ngwf_cg_optimise: &
               &pub_xc_ke_density_required is false, but dfdtau_fine is present &
               &and has a size > 1.")
       end if
    else
       ! JCW: tau-dependent XC functional
       ! JCW: Check dfdtau_fine is present
       call utils_assert( present(dfdtau_fine),&
            "Error in ngwf_cg_optimise: &
            &pub_xc_ke_density_required is true, but dfdtau_fine is not &
            &present.")
    end if

    mrows=rep%overlap%mrows
    ncols=rep%overlap%ncols

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(1)%iscmplx

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine ngwf_cg_optimise not ready yet for more&
         & than one k-point.')

    ! agrecocmplx: should be ok now, need to check everything is
    ! implemented correctly
    !call utils_assert(.not. loc_cmplx, 'Error in&
    !     & ngwf_cg_optimise: not ready yet for complex NGWFs.')  !jmecmplx

    ! agrecokpt: currently only KP method is being implemented
    call utils_assert(pub_kpoint_method == 'KP', &
         'Subroutine ngwf_cg_optimise currently supports&
         & only KP method for BZ sampling')

    ! agrecokpt: default is Gamma point
    if (present(kpt)) then
       loc_kpt%x = kpt(1)
       loc_kpt%y = kpt(2)
       loc_kpt%z = kpt(3)
       loc_kpt_cart(:) = kpt(:)
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
       loc_kpt_cart(:) = 0.0_DP
    end if

    ! ars: get evaluate_functional
    evaluate_functional = &
         utils_devel_code(.false., 'NGWFCG', 'EVALUATEFUNCTIONAL', pub_devel_code)

    ! ars: get evaluate_functional
    force_cubic_step = &
         utils_devel_code(.false., 'NGWFCG', 'FORCECUBICSTEP', pub_devel_code)

    ! ars: define flag for check_conjugacy
    check_conjugacy = &
         utils_devel_code(.false., 'NGWFCG', 'CHECKCONJUGACY', pub_devel_code)

    ! ars: define flag for rotating the previous search direction
    rotate_prev_dir = &
         utils_devel_code(.false., 'NGWFCG', 'ROTATEPREVDIR', pub_devel_code)

    ! ars: orthogonalise search direction
    orthogonalise_dir = &
         utils_devel_code(.false., 'NGWFCG', 'ORTHOGONALISE', pub_devel_code)

    ! ars: define flag for calculating the inner product of the gradient
    calculate_inner_grad = &
         utils_devel_code(.false., 'NGWFCG', 'INNERGRAD', pub_devel_code)

    ! ndmh: define flag for using the old logic for NGWF CG (ie do not evaluate
    ! ndmh: the energy again after a step is taken at the end - useful for
    ! ndmh: comparison to old timings)
    old_cg_logic = &
         utils_devel_code(.false., 'NGWFCG', 'OLD_CG_LOGIC', pub_devel_code)

    ! cks: write initial NGWFs in plotting formats
    if (pub_write_ngwf_plot .and. &
         (pub_cube_format .or. pub_grd_format .or. pub_dx_format) .and. &
         (.not.(pub_eda))) then
       do isub=1,mdl%nsub
          if (mdl%nsub.gt.1) then
             write(isub_str,'(i1)') isub
             call visual_ngwfs(rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                  'initial_'//isub_str, mdl%regions(isub)%elements, &
                  mdl%cell, mdl%fftbox, mdl%regions(isub)%par)
          else
             call visual_ngwfs(rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                  'initial', mdl%regions(isub)%elements, &
                  mdl%cell, mdl%fftbox, mdl%regions(isub)%par)
          end if
       end do
    endif

    ! smmd: write initial NGWFs radial distribution
    if (pub_write_ngwf_radial.ge.1) then
       do isub=1,mdl%nsub
          if (mdl%nsub.gt.1) then
             write(isub_str,'(i1)') isub
             call visual_ngwfs_radial(rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                  'initial_'//isub_str,'ngwf',mdl%fftbox,mdl%cell, &
                  mdl%regions(isub)%par)
          else
             call visual_ngwfs_radial(rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                  'initial','ngwf',mdl%fftbox,mdl%cell, &
                  mdl%regions(isub)%par)
          end if
       end do
    endif

    ! pdh: allocate local sparse matrices
    ! agrecocmplx
    call kernel_create(pur_denskern, denskern%core_struc(1), is_cmplx=loc_cmplx)

    ! rc2013: allocate structures
    call internal_allocate

    if (pub_kernel_christoffel_update) then
       call sparse_embed_array_create(start_denskern, denskern%kern)
       if (pub_aug) call sparse_embed_create(start_sp_overlap,rep%sp_overlap)
    end if
    if (pub_ngwf_cg_rotate) then
       if (.not.pub_kernel_christoffel_update) then
          call sparse_embed_array_create(start_denskern, denskern%kern)
          if (pub_aug) call sparse_embed_create(start_sp_overlap,rep%sp_overlap)
       end if
       if (pub_edft.and..not.pub_foe) then
          allocate(start_mo(1:pub_num_spins), stat=ierr)
          call utils_alloc_check('ngwf_cg_optimise', 'start_mo', ierr)
          do is = 1, pub_num_spins
             ! agrecocmplx: does this need to be complex? I guess so...
             call dense_create(start_mo(is), edft%num, edft%num, &
                  iscmplx=loc_cmplx)
          end do
       end if
       if (((pub_elec_cg_max.gt.0).and.(rotate_prev_dir)).or.pub_foe.or.foe_test_rot) &
            call sparse_embed_create(start_inv_overlap, rep%inv_overlap)
    end if
    allocate(summary_lines(0:max(pub_maxit_ngwf_cg+1,2)),stat=ierr)
    call utils_alloc_check('ngwf_cg_optimise','summary_lines',ierr)

    ! cks: <<< parameter initialisations >>>
    ngwf_threshold     = ngwf_threshold_orig
    lnv_threshold      = lnv_threshold_orig
    mermin_threshold      = mermin_threshold_orig
    est_num_psincs     = function_basis_est_num_psincs(ngwf_basis,mdl%cell)

    ! cks: <<< variable intialisations >>>
    previous_rms_gradient = huge(1.0_DP)
    max_gradient          = huge(1.0_DP)
    line_search_coeff  = 0.15_DP
    !ep : adaptable step with mermin method due to high chebyshev expansions
    if (pub_mermin) then
       trial_length       = 0.1_DP/(pub_mermin_cheb*0.5_DP)
    else
       trial_length       = 0.1_DP
    end if
    F0=0.0_DP ; F1=0.0_DP ; F2=0.0_DP
    rms_gradient=1.0_DP ; mu=0.0_DP ; total_energy=0.0_DP ; cg_count=0
    predicted_functional =0.0_DP
    line_search_success = .true.
    trial2 = .false.
    retrial1 = .false.
    reversing = .false.
    last_n_energies(:) = huge(1.0_DP)
    quadratic_coeff = 0.0_DP
    rejected_quadratic_coeff = 0.0_DP
    cubic_coeff = 0.0_DP
    converged = .false.
    cvstat = 0 ! smmd
    cdft_coarse_converged = .false. ! gibo: cDFT initialistion
    out_of_runtime = .false.
    gradient_valid = .false.
    ! rc2013: allocate and initialise F+T convergence checks
    ft_all_converged = .false.
    allocate(ft_converged(mdl%nsub), stat=ierr)
    call utils_alloc_check('ngwf_cg_optimise','ft_converged',ierr)
    ! rc2013: ngwf_gradient loop
    ift = 1
    ireg = pub_active_region
    fnl_switched = .false.
    !ep mermin outer chemical potential var
    muext=0.0_DP

    ! lr408: Default is false
    updated_shift = .false.

    ! ars: nonsc_forces zero by default
    ngwf_nonsc_forces(:,:) = 0.0_DP

    ! cks: prepare to output summary of NGWF optimisation
    ! rc2013: say which region we're optimising if doing freeze-and-thaw
    if(pub_do_fandt .and. pub_on_root) then
       write(summary_lines(0),'(a88)') '|ITER| REGION |    RMS GRADIENT   &
            &|     TOTAL ENERGY    |   step   |     Epredicted  '
    else
       if (pub_on_root) write(summary_lines(0),'(a80)') '|ITER|    RMS GRADIENT   &
            &|     TOTAL ENERGY    |   step   |     Epredicted  '
    end if
    if (pub_on_root) summary_lines(0) = adjustl(summary_lines(0))

    ! cks: First message for brief output level
    if (pub_on_root .and. pub_output_detail == BRIEF) then
       write(stdout,'(/a)') repeat('#',80)
       write(stdout,'(a)') utils_banner('#','NGWF self-consistent optimisation')
       write(stdout,'(a)') repeat('#',80)
       write(stdout,'(a80)') summary_lines(0)
    end if

    ! qoh: Allow blank calculations for timings purposes if pub_maxit_ngwf_cg < 0
    minit = 1
    if (pub_maxit_ngwf_cg < 0) minit = 0

    ! cks: NGWF iterations loop
    iteration = 0 ! ndmh: prevent unitialised variable when pub_maxit_ngwf_cg = 0

    ! gcc32: initialize preconditioning for PAW
    do isub=1,rep%nsub
       call ngwf_gradient_paw_precond_init(nl_projectors(isub), &
            mdl%regions(isub)%paw_sp,mdl%cell,mdl%fftbox)
    end do

    ! set number of iterations
    if (pub_cond_calculate) then
       current_maxit_lnv     = cond_firstit_lnv
    else if (pub_mermin) then
       current_maxit_mermin  = pub_firstit_mermin
    else
       current_maxit_lnv     = pub_firstit_lnv
    end if

    ! lr408: Update conduction Hamiltonian
    if (pub_cond_calculate) call hamiltonian_dens_dep_nonsc(ham,rep, &
         ngwf_basis,mdl, hfxstate, lhxc_fine,hub,val_rep,val_ham, &
         val_dkn%m(:,PUB_1K), updated_shift, &
         val_ngwf_basis = val_ngwf_basis, &    ! jd: Needed for HFx
         cond_dkn = denskern%kern%m(:,PUB_1K)) ! jd: Needed for HFx

    if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then
       write(stdout,'(/a)') '>>> Optimising kernel for current NGWFs:'
    endif

    ! cks: +++++++++++++++++++++++ OPTIMISE DENSITY KERNEL ++++++++++++++++++
    !gom
    if ( (((.not.pub_cdft).AND.(.not.pub_dft_nu)) .OR. (pub_maxit_cdft_u_cg==0)).OR. &
         ((pub_cdft.OR.pub_dft_nu).AND.(((pub_task=='HUBBARDSCF').AND. & !1stHUBBARDSCF-cDFT iter.
         (hub%consistency_iteration==1)) .OR. pub_hub_on_the_fly)) ) then

       ! agrecokpt: k-point dependence here needed for TB method
       total_energy = electronic_energy(denskern, pur_denskern, ham, &
            lhxc_fine, mu, edft, rep, ngwf_basis, hub_proj_basis, hub, &
            mdl, hfxstate, lnv_threshold, current_maxit_lnv, &
            kernel_update=.true., conv_status=cvstat, &
            dfdtau_fine = dfdtau_fine, kpt=loc_kpt_cart, &
            mermin_threshold=mermin_threshold, &
            current_maxit_mermin = current_maxit_mermin)
    else
       ! JCW: Abort if tau-dependent XC functional requested, as combination of
       ! JCW: this with cDFT not implemented/tested
       call utils_assert(.not.pub_xc_ke_density_required,"Error in &
            &ngwf_cg_optimise: combination of cDFT and tau-dependent XC &
            &functionals has not yet been implemented/tested.")
       ! rc2013: incompatible with embedding
       call utils_assert(mdl%nsub .eq. 1,"Error in &
            &ngwf_cg_optimise: combination of cDFT and embedding &
            &has not yet been implemented/tested.")
       ! gibo: modified for cDFT simultanous NGWFs and Ucdft potential optimisation
       ! agrecokpt: k-point dependence needed because we call electronic_energy?
       call cdft_intermediate_u_cg_optimise(total_energy, cdft_converged, &
            denskern, pur_denskern, ham, &
            lhxc_fine, mu, edft, rep, ngwf_basis, hub_proj_basis, hub, &
            mdl, hfxstate, lnv_threshold, current_maxit_lnv, conv_stat=cvstat, &
            dfdtau_fine = dfdtau_fine, kpt=loc_kpt_cart)
    endif
    ! cks: ++++++++++++++++++++ END OPTIMISE DENSITY KERNEL +++++++++++++++++

   ! set number of iterations
    if (pub_cond_calculate) then
       current_maxit_lnv     = cond_minit_lnv
    else
       if (pub_mermin) then
          current_maxit_mermin  = maxit_mermin
       else
          current_maxit_lnv     = minit_lnv
       end if
    end if

    do iteration=1,max(pub_maxit_ngwf_cg,minit)
       ! rc2013: assign parallel strategy
       par=>mdl%regions(ireg)%par
       ! rc2013: only need a subset of psincs for subsystem optimisations
       if(pub_do_fandt) est_num_psincs = &
            function_basis_est_num_psincs(ngwf_basis(ireg),mdl%cell)

       ! gcc32: set confined ngwfs barrier, if needed
       if(pub_confined_ngwfs) then
         if(iteration<=pub_maxit_ngwf_cg_confined) then
            pub_confined_ngwfs_barrier = pub_confined_ngwfs_barrier* &
                (1.0_DP-real(iteration-1,kind=DP) / &
                real(pub_maxit_ngwf_cg_confined,kind=DP))
         else
            pub_confined_ngwfs = .false.
            pub_confined_ngwfs_barrier = 0.0_DP
            pub_maxit_ngwf_cg_confined = -1
         end if
       end if

       ! rc2013: set all functions to zero for F+T switch
       ! jcap: this also sets the functions to zero at the beginning
       ! of a normal calculation. The link to F+T is that in a F+T
       ! calculation, ift will get reset to 1, so this will be called
       ! again, whilst in other cases it will not
       if(ift == 1) then
          if(pub_debug_on_root) write(stdout,'(a, i4)') &
               'DEBUG: resetting functions for iteration ', iteration
          call data_set_to_zero(G_init)
          call data_set_to_zero(cg_coeff)
          call data_set_to_zero(current_g_dot_g)
          call data_set_to_zero(previous_g_dot_g)
          call data_set_to_zero(prev_direction_on_grid(:))
          call data_set_to_zero(inner_grad)
          !ep : mermin var
          call data_set_to_zero(pre_dir_dot_g)
          call data_set_to_zero(now_dir_dot_g)
          !ep
       endif

       ! rab207: start loop timer
       call timer_check_iteration_time('ngwf_cg','start')

       ! cks: ++++++++++ ITERATION HEADER +++++++++++++++++++++++++++++
       if (pub_on_root .and. pub_output_detail >= NORMAL .and. &
            pub_maxit_ngwf_cg > 0) then
          write(stdout,'(///a)') repeat('#',80)
          if (pub_cond_calculate) then
             write(stdout,'(a,i4.3,a)')'######################### COND NGWF &
                  &CG iteration ',iteration,' ##########################'
          else
             write(stdout,'(a,i4.3,a)')'########################### NGWF &
                  &CG iteration ',iteration,' #############################'
          end if
          write(stdout,'(a)') repeat('#',80)
          if (pub_output_detail >= NORMAL .and. pub_do_fandt) &
               write(stdout,'(a,i4,a)') &
               '*************** Optimising NGWFs in Region', ireg, &
               ' ***************'

          if (pub_output_detail >= VERBOSE) then
             if (pub_precond_recip) then
                write(stdout,'(a,f7.4,a)') &
                     '****** Reciprocal space K.E. preconditioning with &
                     &k_zero = ', pub_k_zero, ' a0^-1 *******'
             end if
             if (pub_precond_real) then
                write(stdout,'(a,f7.4,a)') &
                     '************ Real space K.E. preconditioning with &
                     &k_zero = ', pub_k_zero, ' a0^-1 *******'
             end if
          end if
       end if
       call services_flush
       ! cks: ++++++ END ITERATION HEADER +++++++++++++++++++++++++++++



       ! cks: =========== HYBRID NGWFs==========================================
       if (pub_nnho .and. mdl%nsub==1) then
          call internal_hybridize_ngwfs ! takes care to update rep's matrices
          ! jd: *** NGWFS changed ***
          call ngwf_rep_register_change(rep, &
               'ngwf_cg: internal_hybridize_ngwfs', updated_by_nnho=.true.)
          call ngwf_ham_register_change(ham,rep,updated_by_nnho=.true.)
          ! jd: Re-expand new NGWFs in SWs, if needed
          if(pub_pol_emb_qmstar) then
             ! jcap: the 1 here corresponds to a region - only set up
             ! for one region
             call polarisable_embedding_expand_ngwf_pairs(rep, ngwf_basis(1), mdl)
          end if
          if(pub_use_hfx) then
             if(pub_cond_calculate) then
                call utils_abort('The combination of NNHO, conduction and HFx &
                     &is not supported yet')
             end if
             call hf_exchange_dkn_indep_stage(hfxstate, mdl, 1, &
                  mdl%regions(1)%par, rep, ngwf_basis(1))
          end if
       endif
       ! cks: ======= END HYBRID NGWFs==========================================

       ! lr408: Update Hamiltonian (in conduction case this is unnecessary)
       ! ebl: updating Hamiltonian is also unnecessary for some DMFT calculations
       if (.not.(pub_cond_calculate .or. pub_dmft_spoil_kernel .or. pub_dftb)) then
          if(.not.pub_dmft_fully_sc.or.(pub_dmft_fully_sc.and.pub_dmft_fully_sc_h)) then

             ! mjsp: SCF-MI NGWF gradient calculated without Hamiltonian projection
             if (pub_eda) then
                if (pub_eda_mode == EDA_POLSIMUL) &
                     pub_eda_mode = EDA_POLSIMUL_OFF
                if (pub_eda_mode == EDA_POLFRAGLOC_DEVEL) &
                     pub_eda_mode = EDA_POLFRAGLOC_DEVEL_OFF
                if (pub_eda_mode == EDA_CTFRAGLOC_DEVEL) &
                     pub_eda_mode = EDA_CTFRAGLOC_DEVEL_OFF
             end if

             ! agrecokpt: each component of the hamiltonian has
             ! the proper k-point dependence included at this stage
             call hamiltonian_build_matrix(ham, rep)

             if (pub_eda) then
                if (pub_eda_mode == EDA_POLSIMUL_OFF) &
                     pub_eda_mode = EDA_POLSIMUL
                if (pub_eda_mode == EDA_POLFRAGLOC_DEVEL_OFF) &
                     pub_eda_mode = EDA_POLFRAGLOC_DEVEL
                if (pub_eda_mode == EDA_CTFRAGLOC_DEVEL_OFF) &
                     pub_eda_mode = EDA_CTFRAGLOC_DEVEL
             end if

          end if
       end if

       ! ars: calculate gradient if nonsc forces required
       if ((pub_maxit_ngwf_cg == 0).and..not.pub_nonsc_forces) then
          converged = .false.
          exit
       end if

       if (pub_on_root.and.pub_output_detail>=NORMAL) then
          write(stdout,'(/a)',advance='NO') '>>> Checking for convergence of &
               &NGWFs: '
          call services_flush()
       endif

       ! cks: ~~~~~~~~~~~~~~~~~~~~~~~ NGWF GRADIENT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       if (.not.pub_cond_calculate) then

          ! mjsp: if EDA polarisation calculation, pass the density kernel
          ! represented in the full overlap
          if (pub_eda .and. ((pub_eda_mode == EDA_POLSIMUL).or.&
               (pub_eda_mode == EDA_POLFRAGLOC_DEVEL).or.&
               (pub_eda_mode == EDA_CTFRAGLOC_DEVEL))) then

             ! rc2013: or if doing embedding
             call utils_assert(mdl%nsub .eq. 1,"Error in &
                  &ngwf_cg_optimise: combination of EDA and embedding &
                  &has not yet been implemented/tested.")

             ! agrecokpt: call with kpt argument
             ! jd: ~~~~~~~~~~~~~~~ EDA branch ~~~~~~~~~~~~~~~
             ! rc2013: structures not suitable for embedding
             call ngwf_gradient_lnv(contra_grad_on_grid, cov_grad_on_grid, & ! in/out
                  denskern_R, rep, ngwf_basis, proj_basis, hub_proj_basis, & ! in
                  nl_projectors, lhxc_fine, ham, hub, mu, mdl, hfxstate, &   ! in
                  ireg, kpt=loc_kpt, &                                       ! in
                  dfdtau_fine = eda_dfdtau_fine)
             !     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else
             if(pub_pol_emb_vacuum_qmstar) then
                ! jd: ~~~~~~~~~~~~ polarisable embedding branch ~~~~~~~~~~~~
                call ngwf_gradient_lnv(contra_grad_on_grid, cov_grad_on_grid, & ! in/out
                     denskern, rep, ngwf_basis, proj_basis, hub_proj_basis, &   ! in
                     nl_projectors, lhxc_fine, ham, hub, mu, mdl, hfxstate, &   ! in
                     ireg, kpt=loc_kpt, &                                       ! in
                     dfdtau_fine = dfdtau_fine)                                 ! in
                !     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             else if (pub_mermin) then
                if (pub_mermin_mu_sq) then
                  ! poli: ~~~~~~~~~~~ mermin branch mu_square~~~~~~~~~
                       call ngwf_gradient_lnv(contra_grad_on_grid,  &
                            cov_grad_on_grid, denskern, rep, ngwf_basis, &
                            proj_basis, hub_proj_basis, nl_projectors, &
                            lhxc_fine, ham, hub, mu, mdl, hfxstate, ireg, &
                            muext, kpt=loc_kpt, step=trial_length, &
                            dfdtau_fine = dfdtau_fine)
                  !     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                else
                  ! poli: ~~~~~~~~~~~~~~~ mermin branch ~~~~~~~~~~~~~~~
                       call ngwf_gradient_lnv(contra_grad_on_grid, &
                            cov_grad_on_grid, denskern, rep, ngwf_basis, &
                            proj_basis, hub_proj_basis, nl_projectors, lhxc_fine, &
                            ham, hub, mu, mdl, hfxstate, ireg, muext, &
                            kpt=loc_kpt, dfdtau_fine = dfdtau_fine)
                   !     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                end if
             else
                ! jd: ~~~~~~~~~~~~~~~ default branch ~~~~~~~~~~~~~~~
                call ngwf_gradient_lnv(contra_grad_on_grid, cov_grad_on_grid, & ! in/out
                     denskern, rep, ngwf_basis, proj_basis, hub_proj_basis, &   ! in
                     nl_projectors, lhxc_fine, ham, hub, mu, mdl, hfxstate, &   ! in
                     ireg, kpt=loc_kpt, &                                       ! in
                     dfdtau_fine = dfdtau_fine)                                 ! in
                !     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             end if

          end if

       else
          ! JCW: Abort if tau-dependent XC functional requested
          call utils_assert(.not.pub_xc_ke_density_required,"Error in &
               &ngwf_cg_optimise: combination of conduction NGWFs and &
               &tau-dependent XC functionals has not yet been &
               &implemented/tested.")
          ! lr408: Need extra arguments for conduction calculation
          ! agrecokpt: call with kpt argument
          ! jd: ~~~~~~~~~~~~~~~ conduction branch ~~~~~~~~~~~~~~~
          call ngwf_gradient_lnv(contra_grad_on_grid, cov_grad_on_grid, & ! in/out
               denskern, rep, ngwf_basis, proj_basis, hub_proj_basis, &   ! in
               nl_projectors, lhxc_fine, ham, hub, mu, mdl, hfxstate, &   ! in
               ireg, val_dkn=val_dkn, val_rep=val_rep,                &   ! in
               val_ham=val_ham%ham, val_ngwf_basis=val_ngwf_basis, &      ! in
               cond_shift=ham%cond_shift, kpt=loc_kpt)                   ! in
          !     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       end if
       gradient_valid = .true.

       ! cks: ~~~~~~~~~~~~~~~~~~~~ END NGWF GRADIENT ~~~~~~~~~~~~~~~~~~~~~~~~~~~

       ! ars: orthogonalise SD search direction
       do isub=1,rep%nsub
          ! rc2013: EMBED_FIX!
          if (orthogonalise_dir) &
               call internal_orthogonalise_direction(cov_grad_on_grid(isub))
       end do

       ! ndmh: check conjugacy condition on search direction
       if (check_conjugacy) then
          ! agrecocmplx: I think it is correct to consider the complex dot
          ! product here, hence taking the conjugate of the first argument.
          ! This is done automatically in data_functions_dot
          previous_dir_dot_g = data_functions_dot(contra_grad_on_grid, &  !jmecmplx
               prev_direction_on_grid)                                     !jmecmplx
          if (previous_dir_dot_g%iscmplx) then
             call comms_reduce('SUM', previous_dir_dot_g%z)
          else
             call comms_reduce('SUM', previous_dir_dot_g%d)
          end if
       end if

       ! cks: ###################### TEST CONVERGENCE ##########################
       converged = internal_test_convergence()

       ! rc2013: only exit if convergence satisfied in all regions
       ft_all_converged = internal_ft_convergence()
       call comms_bcast(pub_root_proc_id,ft_all_converged)
       if (ft_all_converged) then

          exit

       else if(pub_quit_region .and. converged) then
          ! rc2013: we're not really converged...
          converged = .false.
          ! rc2013: but we do need the summary info!
          if (pub_on_root) then
             if(abs(total_energy)<100000.0_DP) then
                if(pub_do_fandt) then
                   write(summary_lines(iteration),'(i4,i8,f21.14,f22.14,f11.6,f22.14)')&
                        iteration-1,ireg,rms_gradient,total_energy,line_search_coeff, &
                        predicted_functional
                else
                   write(summary_lines(iteration),'(i4,f21.14,f22.14,f11.6,f22.14)')&
                        iteration-1,rms_gradient,total_energy,line_search_coeff, &
                        predicted_functional
                end if
             else
                if(pub_do_fandt) then
                   write(summary_lines(iteration),'(i4,i8,f21.14,f22.12,f11.6,f22.12)')&
                        iteration-1,ireg,rms_gradient,total_energy,line_search_coeff, &
                        predicted_functional
                else
                   write(summary_lines(iteration),'(i4,f21.14,f22.12,f11.6,f22.12)')&
                        iteration-1,rms_gradient,total_energy,line_search_coeff, &
                        predicted_functional
                end if
             end if
             if (pub_output_detail == BRIEF) then
                write(stdout,'(a80)') summary_lines(iteration)
             end if
          end if
          ! rc2013: move into the next region and skip the rest of this iteration
          ! rc2013: unless we're already doing EMFT
          if(pub_do_fandt .or. (pub_emft_follow .and. .not. fnl_switched)) &
               call internal_region_reset

          ! rab207: check timings -- if likely to run out of time, quit now
          call timer_check_iteration_time('ngwf_cg','stop',out_of_runtime, &
               global_mpi_loop=.true.)

          cycle
       end if
       ! cks: #################### END TEST CONVERGENCE ########################

       ! pdh: exit if no NGWF optimisation required
       if (pub_maxit_ngwf_cg == 0) then
          converged = .false.
          exit
       end if

       ! ndmh: kernel will now be adjusted, so gradient will no longer be
       ! ndmh: consistent with kernel
       gradient_valid = .false.

       !ep : Mermin adjustment to change between mu-internal
       !ep : Mermin adjustment to change between mu-external (real chem. pot.)
       if (pub_mermin) then
            call mermin_switch_mu(total_energy, rep%overlap, &
                denskern%kern%m(:,PUB_1K), mu(:), rep%n_occ)
       end if

       ! cks: **************** FUNCTIONAL AT INITIAL POINT *********************
       F0 = electronic_lagrangian(total_energy, rep%overlap, &
            denskern%kern%m(:,PUB_1K), &
            ham%ham, mu, rep%n_occ, edft, muext=muext)
       ! cks: ************** END FUNCTIONAL AT INITIAL POINT *******************


       ! cks: ************************ LINE SEARCH *****************************
       ! the density kernel is not updated during the line search

       call internal_find_direction

       ! rc2013: EMBED_FIX!
       do isub=1,mdl%nsub
          ! cks: store direction if doing conjugate gradients
          if (pub_elec_cg_max > 0) then
             ! ars: orthogonalise CG direction
             if (orthogonalise_dir) &
                  call internal_orthogonalise_direction(direction_on_grid(isub))
             call data_functions_copy(prev_direction_on_grid(isub),direction_on_grid(isub))
          endif

          if (pub_ngwf_cg_type == 'NGWF_POLAK') &
               call data_functions_copy(prev_contra_grad_on_grid(isub),contra_grad_on_grid(isub))

          previous_g_dot_g = current_g_dot_g

          ! ndmh_pointerfun
          call data_functions_copy(start_ngwfs_on_grid(isub), &
               rep%ngwfs_on_grid(isub))
       end do

       ! ep: for mermin search (to be checked)
       if (pub_mermin) then
          do isub=1,mdl%nsub
             if (iteration .gt. 1) then
                call data_functions_copy(prev_contra_grad_on_grid(isub),contra_grad_on_grid(isub))
             end if
          end do
       end if


       ! ndmh: if doing Christoffel updates to the kernel, save the denskern at
       ! the initial point now
       ! ars: or if rotating to new NGWF representation
       if (pub_kernel_christoffel_update.or.pub_ngwf_cg_rotate) then
          call sparse_embed_array_copy(start_denskern, denskern%kern)
          do is=1,pub_num_spins
             if (pub_aug) call sparse_embed_copy(start_sp_overlap,rep%sp_overlap)
             if (pub_edft.and..not.pub_foe) call dense_copy(start_mo(is), edft%mo(is))
             if (((pub_elec_cg_max.gt.0).and.(rotate_prev_dir)).or.pub_foe.or.foe_test_rot) &
                  call sparse_embed_copy(start_inv_overlap, rep%inv_overlap)
          end do
       end if


       if (index(pub_devel_code,'NGWF_FD')>0) then
          do ifd=1,nfd
             do isub=1,mdl%nsub
                call data_functions_copy(rep%ngwfs_on_grid(isub), &
                     start_ngwfs_on_grid(isub))
                call data_functions_axpy(rep%ngwfs_on_grid(isub), &
                     direction_on_grid(isub),fd_trial_step(ifd))
             end do
             ! jd: *** NGWFS changed ***
             call ngwf_rep_register_change(rep,'ngwf_cg: NGWF_FD')
             call ngwf_ham_register_change(ham,rep)

             ! ndmh: %%%%%%%%%%%%% FUNCTIONAL AT FD STEP (DEBUG) %%%%%%%%%%%%%%%
             FFD(ifd) = internal_step_energy()
             if ((pub_output_detail>=VERBOSE).and.(.not.pub_cond_calculate)) then
                if (.not.(pub_eda_scfmi) .and. .not. (pub_mermin)) then
                   call hamiltonian_energy_components(  &
                        pur_denskern%kern%m(:,PUB_1K), rep, mdl, &
                        ngwf_basis, hub_proj_basis, hub, &
                        ham%hfexchange, edft=edft)
                else if (pub_mermin) then
                   !ep: call Hamiltonian component with Mermin params
                   call hamiltonian_energy_components(  &
                        denskern%kern%m(:,PUB_1K), rep, mdl, &
                        ngwf_basis, hub_proj_basis, hub, &
                        ham%hfexchange,muext(:,PUB_1K))
                else
                   ! mjsp: if SCF-MI then use the supermolecule density kernel
                   ! mjsp: representation
                   call hamiltonian_energy_components(  &
                        denskern_R%kern%m(:,PUB_1K), rep, mdl, &
                        ngwf_basis, hub_proj_basis, hub, &
                        ham%hfexchange)
                end if
             end if
             ! ndmh: %%%%%%%%%%%%% FUNCTIONAL AT FD STEP (DEBUG) %%%%%%%%%%%%%%%
          end do
       end if

       ! Update functions by trial step * search direction
       do isub=1,mdl%nsub
          call data_functions_copy(rep%ngwfs_on_grid(isub), &
               start_ngwfs_on_grid(isub))
          call data_functions_axpy(rep%ngwfs_on_grid(isub), &
               direction_on_grid(isub),trial_length)
       end do
       ! jd: *** NGWFS changed ***
       call ngwf_rep_register_change(rep,'ngwf_cg: moved by trial length 1')
       call ngwf_ham_register_change(ham,rep)

       ! cks: %%%%%%%%%%%%%%%%%% FUNCTIONAL AT TRIAL LENGTH 1 %%%%%%%%%%%%%%%%%%
       F1 = internal_step_energy()
       ! cks: %%%%%%%%%%%%%%%% END FUNCTIONAL AT TRIAL LENGTH 1 %%%%%%%%%%%%%%%%

       ! cks: ===================== SEARCH BY FITTING PARABOLA =================
       ! agrecocmplx: G_init is the slope of energy in search direction, consider
       ! only the real part in complex case? add a check on imaginary part?
       if (loc_cmplx) then
          G_init_real = real(G_init%z,kind=DP)
       else
          G_init_real = G_init%d
       end if

       call services_line_search_parabola(&
            quadratic_coeff, predicted_functional,line_search_success, & !output
            G_init_real, F0, F1, trial_length, pub_ngwf_cg_max_step)         !input  !jmecmplx

       line_search_coeff = quadratic_coeff

       ! cks: ================= END SEARCH BY FITTING PARABOLA =================


       ! cks: CUBIC CUBIC CUBIC ----- SEARCH BY FITTING CUBIC ---- CUBIC CUBIC
       if (force_cubic_step.or.(quadratic_coeff * G_init_real > 0.0_DP)) then !jmecmplx

          trial2 = .true.

          ! Update functions by 2 * trial step * search direction
          do isub=1,mdl%nsub
             call data_functions_copy(rep%ngwfs_on_grid(isub), &
                  start_ngwfs_on_grid(isub))
             call data_functions_axpy(rep%ngwfs_on_grid(isub), &
                  direction_on_grid(isub),2.0_DP*trial_length)
          end do
          ! jd: *** NGWFS changed ***
          call ngwf_rep_register_change(rep, 'ngwf_cg: moved by trial length 2')
          call ngwf_ham_register_change(ham,rep)

          ! cks: %%%%%%%%%%%%%%%% FUNCTIONAL AT TRIAL LENGTH 2 %%%%%%%%%%%%%%%%%
          F2 = internal_step_energy()
          ! cks: %%%%%%%%%%%%% END FUNCTIONAL AT TRIAL LENGTH 2 %%%%%%%%%%%%%%%%

          ! Cubic Fit
          call services_cubic_fit_minimum( &
               cubic_coeff, predicted_functional, line_search_success, & !output
               F0, F1, F2, G_init_real, trial_length, 2.0_DP*trial_length, &  !input !jmecmplx
               pub_ngwf_cg_max_step)                                         !input

          line_search_coeff = cubic_coeff

       end if
       ! cks: CUBIC CUBIC --------- END SEARCH BY FITTING CUBIC -------- CUBIC

       ! ndmh: SECOND TRIAL STEP ----- SECOND TRIAL STEP ----- SECOND TRIAL STEP
       ! ndmh: protection against bad line search results: redo trial step if
       ! ndmh: line search result is too much bigger or smaller than trial step
       ! ndmh: or if we are searching 'uphill' and went past trial step position
       if (.not.force_cubic_step.and.( ((line_search_coeff < 0.05_DP*trial_length) .or. &
            (line_search_coeff > 10.0_DP*trial_length) .or. &
            (.not. line_search_success) .or. &
            (reversing .and. (line_search_coeff > trial_length))) &
            .and. .not. trial2 ) ) then

          retrial1 = .true.
          rejected_quadratic_coeff = quadratic_coeff
          if ((rejected_quadratic_coeff == 0.15_DP).and. &
               (.not.line_search_success)) then
             ! Try fit again, relaxing maximum CG step, to find sensible trial step
             call services_line_search_parabola(&
                  quadratic_coeff, predicted_functional,line_search_success, & !output
                  G_init_real, F0, F1, trial_length, 10.0_DP*pub_ngwf_cg_max_step) !input
          end if

          do isub=1,mdl%nsub
             call data_functions_copy(rep%ngwfs_on_grid(isub), &
                  start_ngwfs_on_grid(isub))
             call data_functions_axpy(rep%ngwfs_on_grid(isub), &
                  direction_on_grid(isub),quadratic_coeff)
          end do
          ! jd: *** NGWFS changed ***
          call ngwf_rep_register_change(rep,'ngwf_cg: moved by retrial step')
          call ngwf_ham_register_change(ham,rep)

          ! cks: %%%%%%%%%%%%%%%% FUNCTIONAL AT TRIAL LENGTH 2 %%%%%%%%%%%%%%%%%
          F2 = internal_step_energy()
          ! cks: %%%%%%%%%%%%% END FUNCTIONAL AT TRIAL LENGTH 2 %%%%%%%%%%%%%%%%

          ! ndmh: quadratic fit at new trial length
          call services_line_search_parabola(&
               quadratic_coeff, predicted_functional, &               !output
               line_search_success, G_init_real, F0, F2, &            !input !jmecmplx
               line_search_coeff, pub_ngwf_cg_max_step)                   !input

          line_search_coeff = quadratic_coeff

       end if
       ! ndmh: SECOND TRIAL STEP ----- SECOND TRIAL STEP ----- SECOND TRIAL STEP

       if (.not.(trial2 .or. retrial1)) F2 = 0.0_DP

       ! cks: set the new values of the NGWFs
       do isub=1,mdl%nsub
          call data_functions_copy(rep%ngwfs_on_grid(isub), &
               start_ngwfs_on_grid(isub))
          call data_functions_axpy(rep%ngwfs_on_grid(isub), &
               direction_on_grid(isub), line_search_coeff)
       end do
       ! jd: *** NGWFS changed ***
       call ngwf_rep_register_change(rep, &
            'ngwf_cg: moved to final NGWFs following search')
       call ngwf_ham_register_change(ham,rep)

       Ffinal = internal_step_energy(evaluate_functional)

       if (pub_on_root .and. pub_output_detail >= NORMAL) &
            call internal_print_search_info

       ! cks: write new NGWFs to file if required
       if (pub_write_tightbox_ngwfs.and.(.not.pub_write_converged_dk_ngwfs)) then
          ! cks: in universal tightbox representation
          do isub=1,mdl%nsub
             write(file_suffix,'(a30)') ngwf_basis(isub)%name
             if(mdl%nsub .gt. 1) then
                write(isub_str,'(i1)') isub
                write(file_suffix,'(a30)') trim(file_suffix)//"_"//isub_str
             endif
             file_suffix = adjustl(trim(file_suffix))
             call restart_ngwfs_tightbox_output(rep%ngwfs_on_grid(isub),&
                  ngwf_basis(isub),mdl%cell,mdl%fftbox,&
                  mdl%regions(isub)%elements,'tightbox_'//file_suffix,&
                  mdl%regions(isub))
          end do
       endif
       ! cks: ************************* END LINE SEARCH ***********************

       ! ndmh: reset flags
       trial2 = .false.
       retrial1 = .false.
       reversing = .false.
       F2 = 0.0_DP

       ! cks: write line with summary info for curent iteration
       if (pub_on_root) then
          if(abs(total_energy)<100000.0_DP) then
             if(pub_do_fandt) then
                write(summary_lines(iteration),'(i4,i8,f21.14,f22.14,f11.6,f22.14)')&
                     iteration-1,ireg,rms_gradient,total_energy,line_search_coeff, &
                     predicted_functional
             else
                write(summary_lines(iteration),'(i4,f21.14,f22.14,f11.6,f22.14)')&
                     iteration-1,rms_gradient,total_energy,line_search_coeff, &
                     predicted_functional
             end if
          else
             if(pub_do_fandt) then
                write(summary_lines(iteration),'(i4,i8,f21.14,f22.12,f11.6,f22.12)')&
                     iteration-1,ireg,rms_gradient,total_energy,line_search_coeff, &
                     predicted_functional
             else
                write(summary_lines(iteration),'(i4,f21.14,f22.12,f11.6,f22.12)')&
                     iteration-1,rms_gradient,total_energy,line_search_coeff, &
                     predicted_functional
             end if
          end if
          if (pub_output_detail == BRIEF) then
             ! rc2013: make sure we print all F+T info!
             write(stdout,'(a88)') summary_lines(iteration)
          end if
       end if

       if (old_cg_logic.and.(iteration==max(pub_maxit_ngwf_cg,minit))) cycle

       if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then
          write(stdout,'(/a)') '>>> Optimising kernel for current NGWFs:'
       endif

       ! cks: +++++++++++++++++++++++ OPTIMISE DENSITY KERNEL ++++++++++++++++++
       !gom
       if ( (((.not.pub_cdft).AND.(.not.pub_dft_nu)) .OR. (pub_maxit_cdft_u_cg==0)).OR. &
            ((pub_cdft.OR.pub_dft_nu).AND.(((pub_task=='HUBBARDSCF').AND. & !1stHUBBARDSCF-cDFT iter.
            (hub%consistency_iteration==1)) .OR. pub_hub_on_the_fly)) ) then

          total_energy = electronic_energy(denskern, pur_denskern, ham, &
               lhxc_fine, mu, edft, rep, ngwf_basis, hub_proj_basis, hub, &
               mdl, hfxstate, lnv_threshold, current_maxit_lnv, &
               kernel_update=.true., conv_status=cvstat, &
               dfdtau_fine = dfdtau_fine, &
               mermin_threshold=mermin_threshold, &
               current_maxit_mermin=current_maxit_mermin)
       else
          ! JCW: Abort if tau-dependent XC functional requested, as combination of
          ! JCW: this with cDFT not implemented/tested
          call utils_assert(.not.pub_xc_ke_density_required,"Error in &
               &ngwf_cg_optimise: cDFT and tau-dependent XC functionals have not &
               &yet been implemented/tested.")
          ! rc2013: incompatible with embedding
          call utils_assert(mdl%nsub .eq. 1,"Error in &
               &ngwf_cg_optimise: combination of cDFT and embedding &
               &has not yet been implemented/tested.")
          !gibo: modified for cDFT simultanous NGWFs and Ucdft potential optimisation
          call cdft_intermediate_u_cg_optimise(total_energy, cdft_converged, &
               denskern, pur_denskern, ham, &
               lhxc_fine, mu, edft, rep, ngwf_basis, hub_proj_basis, hub, &
               mdl, hfxstate, lnv_threshold, current_maxit_lnv, &
               conv_stat=cvstat, dfdtau_fine = dfdtau_fine)
       endif
       ! cks: ++++++++++++++++++++ END OPTIMISE DENSITY KERNEL +++++++++++++++++


       !ep : Mermin adjustment to change between mu
       if (pub_mermin) then
            call mermin_switch_mu(total_energy, rep%overlap, &
                denskern%kern%m(:,PUB_1K), mu(:), rep%n_occ)
            Ftmp = electronic_lagrangian(total_energy, rep%overlap, &
                 denskern%kern%m(:,PUB_1K), &
                 ham%ham, mu, rep%n_occ, edft, muext=muext)
            total_energy=Ftmp
       endif

       ! cks: ===== write out different energy components every so often =======
       ! lr408: Don't print out components if this is a conduction calculation
       ! rc2013: But do print if we're swapping regions
       if (((mod(iteration, 5) == 0 .and. pub_output_detail >= VERBOSE) .or. &
            (iteration == pub_maxit_ngwf_cg).or.(ift==pub_freeze_switch_steps))&
            .and. (.not.pub_cond_calculate)) then

          if (.not. (pub_eda_scfmi) .and. .not.(pub_mermin)) then
             call hamiltonian_energy_components(pur_denskern%kern%m(:,PUB_1K), &
                  rep, mdl, ngwf_basis, hub_proj_basis, hub, ham%hfexchange, &
                  edft=edft)

          else if (pub_mermin) then
             call hamiltonian_energy_components(denskern%kern%m(:,PUB_1K), &
                  rep, mdl, ngwf_basis, hub_proj_basis, hub, ham%hfexchange, &
                  muext(:,PUB_1K))
          else
             ! mjsp: if SCF-MI then use the supermolecule density kernel
             ! mjsp: representation
             call hamiltonian_energy_components(denskern_R%kern%m(:,PUB_1K), &
                  rep, mdl, ngwf_basis, hub_proj_basis, hub, ham%hfexchange)
          end if

       end if

       ! cks: ==================================================================

       ! rab207: check timings -- if likely to run out of time, quit now
       call timer_check_iteration_time('ngwf_cg','stop',out_of_runtime,global_mpi_loop=.true.)
       if (out_of_runtime) exit

       ! rc2013: move onto the next region if we've done enough steps here
       ! Or the gradient in this region is below the threshold
       ! Or if we're at the end of the NGWF optimisation and doing EMFT
       if( ift == pub_freeze_switch_steps .or. converged .or. &
            iteration==max(pub_maxit_ngwf_cg,minit)) then
          call internal_region_reset
       else
          ift = ift + 1
       endif
    end do ! rc2013: end of CG iteration

    ! rab207: cancel timer
    call timer_check_iteration_time('ngwf_cg','cancel')

    ! write NGWFs if required and minimisation ends in first iteration
    ! ndmh: or if writing final converged results only
    do isub=1,mdl%nsub
       ! rc2013: attach regional string to file name if needed
       write(file_suffix,'(a30)') ngwf_basis(isub)%name
       if(mdl%nsub .gt. 1) then
          write(isub_str,'(i1)') isub
          write(file_suffix,'(a30)') trim(file_suffix)//"_"//isub_str
       endif
       file_suffix = adjustl(trim(file_suffix))
       if (pub_write_tightbox_ngwfs .and. ((iteration.eq.1).or. &
            pub_write_converged_dk_ngwfs)) then
          call restart_ngwfs_tightbox_output(rep%ngwfs_on_grid(isub),&
               ngwf_basis(isub),mdl%cell,mdl%fftbox,&
               mdl%regions(isub)%elements, 'tightbox_'//trim(file_suffix),&
               mdl%regions(isub))
       endif
       ! ars: write in spherical waves representation
       if (pub_write_sw_ngwfs) then! .and. iteration.eq.1) then
          ! rc2013: use basis name to write these NGWFs instead
          call restart_sph_waves_output(rep%ngwfs_on_grid(isub),&
               ngwf_basis(isub),mdl%fftbox,mdl%uni_tightbox,mdl%cell,&
               mdl%regions(isub)%elements,'sw_'//trim(file_suffix), &
               mdl%regions(isub))
       endif
    end do

    ! ndmh: write denskern here if we have not been writing it as we go along
    if (pub_write_denskern .and. pub_write_converged_dk_ngwfs) &
         call restart_kernel_write(denskern%kern, write_cond=pub_cond_calculate)

    ! ndmh: write hamiltonian in sparse format if requested
    if (pub_write_hamiltonian) then
       if(.not.pub_edft) then ! JA: The correct Hamiltonian for EDFT was written in EDFT.
          call restart_hamiltonian_write(ham%ham, write_cond=pub_cond_calculate)
       end if
       if (pub_perturbative_soc) then
           call internal_write_soc_ham()
       end if
    end if
    ! ndmh: write overlap in sparse format if requested
    if (pub_write_overlap) then
       call restart_overlap_write(rep%overlap, write_cond=pub_cond_calculate)
    end if


    ! ars: calculate ngwf_nonsc_forces if required
    if (pub_nonsc_forces .and. pub_forces_needed &
         .and.(.not.(pub_cond_calculate))) then

       if (.not.gradient_valid) then

          ! recalculate gradient
          call ngwf_gradient_lnv(contra_grad_on_grid, cov_grad_on_grid, & ! in/out
               denskern, rep, ngwf_basis, proj_basis, hub_proj_basis, &   ! in
               nl_projectors, lhxc_fine, ham, hub, mu, mdl, hfxstate, &   ! in
               ireg, kpt=loc_kpt, &                                       ! in
               dfdtau_fine = dfdtau_fine)                                 ! in

       end if

       ! jcap: calculate nonsc forces subregion by subregion, as
       ! we only need the diagonal terms anyway
       do isub=1,mdl%nsub
          ! jcap: allocate temporary array
          allocate(sub_nonsc_forces(1:3,mdl%regions(isub)%par%nat),stat=ierr)
          call utils_alloc_check('energy_and_force_calculate', &
               'sub_nonsc_forces', ierr)
          sub_nonsc_forces=0.0_DP

          if((pub_output_detail >= VERBOSE).and.(mdl%nsub.gt.1)) then
             write(suffix,'(i0)')isub
             call utils_flushed_string_output('Region '//trim(suffix)//':')
          end if
          ! jcap: subregion suffix
          ! rc2013: only if we actually need it
          suffix = trim(rep%postfix)
          if(mdl%nsub .gt. 1) &
               write(suffix,'(a,i0,i0)') trim(rep%postfix),isub,isub
          call nonsc_forces_ngwfs_calc(sub_nonsc_forces,&
               rep%ngwfs_on_grid(isub), contra_grad_on_grid(isub), &
               &ngwf_basis(isub), mdl%cell, mdl%fftbox, &
               mdl%regions(isub)%par, suffix=trim(suffix))

          ! jcap: copy data into full array
          do iat=1,mdl%regions(isub)%par%nat
             ngwf_nonsc_forces(:,mdl%regions(isub)%elements(iat)%global_atom_number)&
                  = sub_nonsc_forces(:,iat)
          end do
          deallocate(sub_nonsc_forces,stat=ierr)
          call utils_dealloc_check('energy_and_force_calculate', &
              'sub_nonsc_forces',  ierr)
       end do

    end if

    ! smmd: store NGWFs for extrapolation at subsequent MD/geom steps
    do isub=1,mdl%nsub
       par=>mdl%regions(isub)%par
       call elec_history_store_ngwfs(rep%ngwfs_on_grid(isub),ngwf_basis(isub),&
            mdl%cell,mdl%fftbox,mdl%regions(isub)%elements,isub,'scf', &
            par = par)
    end do

    deallocate_workspace = .false.
    if ( mix_dkn_type == 'XLI' .or. mix_dkn_type == 'XLD' .or. &
         mix_dkn_type == 'NAIVE' .or. mix_dkn_type == 'XLIS') then

       if ( md_aux_rep == 'ORTHO' .or. md_aux_rep == 'ortho') then
          !vv: Initialise workspace
          lwd_overlap%matrix_type = 4
          tmp_overlap%matrix_type = 4
          lwd_overlap%dataSPAM3_EMBED%structure = 'K'
          tmp_overlap%dataSPAM3_EMBED%structure = 'K'
          ks_half%structure = 'K'
          call sparse_embed_create(ks_half)
          ! agrecocmplx: do these need to be complex? if so, we need to set
          ! standard_is_cmplx to true, and call routines with iscmplx=.true.
          call sparse_embed_create(lwd_overlap%dataSPAM3_EMBED)
          call sparse_embed_create(tmp_overlap%dataSPAM3_EMBED)
          call sparse_embed_array_create(Qmat, n_kpoints=PUB_1K, n_spins=pub_num_spins, &
               structure='K')

          ! agrecocmplx: danger of unsafe complex to real conversion
          ! if tmp_overlap is kept real
          call sparse_embed_copy(tmp_overlap%dataSPAM3_EMBED,rep%overlap)

          ! vv: Use Lowdin transformation to compute S^(1/2)
          call lowdin_transformation(tmp_overlap,sqrt_S_out=lwd_overlap)
          do is = 1,pub_num_spins
             call sparse_embed_product(ks_half,denskern%kern%m(is,PUB_1K), &
                  lwd_overlap%dataSPAM3_EMBED)
             call sparse_embed_product(Qmat%m(is,PUB_1K),lwd_overlap%dataSPAM3_EMBED,ks_half)
          end do

          call elec_history_store_dkn(Qmat,'scf')

          call sparse_embed_destroy(lwd_overlap%dataSPAM3_EMBED)
          call sparse_embed_destroy(tmp_overlap%dataSPAM3_EMBED)
          call sparse_embed_destroy(ks_half)
          call sparse_embed_array_destroy(Qmat)

       else if (md_aux_rep == 'ASYM' .or. md_aux_rep == 'asym') then
          if(.not. denskern%workspace_created) then
             call kernel_workspace_create(denskern,rep%overlap)
             call kernel_validate_ks(denskern,rep%overlap)
             deallocate_workspace = .true.
          end if
          call elec_history_store_dkn(denskern%ks,'scf')
          if (deallocate_workspace) call kernel_workspace_destroy(denskern)
       end if

    else
       ! smmd: store denskern for extrapolation at subsequent MD/geom steps
       call elec_history_store_dkn(denskern%kern,'scf')
    end if

    ! vv : If Berendsen thermostat is applied to the auxiliary dof
    if ( pub_task == 'MOLECULARDYNAMICS' .and. (mix_dkn_type == 'XLI' .or. &
         mix_dkn_type == 'NAIVE' .or. mix_dkn_type == 'XLIS')) then
       call elec_history_update_dkn()
    end if

    ! write temperature of KS if needed
    if ( pub_task == 'MOLECULARDYNAMICS' .and. &
         md_iter_global >= mix_dkn_init_num + 1 .AND. &
        (mix_dkn_type == 'XLD' .or.     &
         mix_dkn_type == 'XLI'  .or. &
         mix_dkn_type == 'XLIS' )) then
       call elec_history_compute_temp(mix_dkn_type,'real', &
            pub_num_kpoints, pub_num_spins)
       !call elec_history_compute_temp(mix_dkn_type,'aux', &
       !     pub_num_kpoints, pub_num_spins)
    end if

    !   vv: backup electronic history
    ! agrecocmplx
    if (pub_task == 'MOLECULARDYNAMICS' .and.  (md_iter_global >= md_write_history) .and. &
         md_write_history > 0 .and. mod(md_iter_global,md_write_history) == 0) then
       call elec_history_backup_history(mdl%nat, mrows, ncols, loc_cmplx, &
            mdl%cell, mdl%fftbox)
    end if

    ! gibo: Successful NGWFs-cDFT-optimisation!
    !  print out cDFT-projector (for cDFT-CI runs)
    ! [for HUBBARDSCF+cDFT, let sbrtne hubbard_projector_consistency do it]
    ! ddor: write new Hubbard NGWF projectors to file in tightbox rep.
    !       write out all NGWFs, so we won't need projector tightboxes,
    !       or write out hub%consistency_overlap.
    ! rab207: TODO does this need storing if run out of runtime?
    ! gom
    if (converged .AND. (pub_task/='HUBBARDSCF') .AND. (pub_cdft.OR.pub_dft_nu)) then
        call restart_ngwfs_tightbox_output(&
             hub%consistency_ngwfs,ngwf_basis(1),mdl%cell,mdl%fftbox, &
             mdl%regions(1)%elements, 'tightbox_hub_projs',mdl%regions(1))
    end if

    ! gibo: prevent progressive convergence tigthening at each NGWFs-opt call
    if (pub_cdft.and.pub_cdft_tight)then
       if (cdft_coarse_converged) then
          pub_cdft_cg_threshold = 10._DP*pub_cdft_cg_threshold
          pub_cdft_max_grad = 10._DP*pub_cdft_max_grad
          !=== (possibly redundant) make all procs aware of the changes
          call comms_bcast(pub_root_proc_id, pub_cdft_max_grad)
          call comms_bcast(pub_root_proc_id, pub_cdft_cg_threshold)
          !=== (possibly redundant) make all procs aware of the changes
       endif
    endif


    if (.not. converged .and. .not. out_of_runtime) then
       ! cks: reset iteration number just for storing final line of calculation
       !      summary
       iteration = iteration - 1

       ! cks: print warning that calculation failed to converge
       if (pub_on_root) write(stdout,'(a,i4,a)') &
            'WARNING: maximum number of NGWF CG iterations (',pub_maxit_ngwf_cg, &
            ') exceeded!'
    end if

    ! JCW: If non-self-consistent energy calculation (zero iterations), simply
    ! JCW: calculate the energy components and print
    nonscenergy = utils_devel_code(.false.,'NGWFCG','NONSCENERGY',pub_devel_code)
    if (pub_debug_on_root) then
       write(stdout,'(a,l1)') "DEBUG: NONSCENERGY = ", nonscenergy
    end if
    if (nonscenergy) then
       call hamiltonian_energy_components(  &
            pur_denskern%kern%m(:,PUB_1K), rep, mdl, ngwf_basis, &
            hub_proj_basis, hub, ham%hfexchange, edft=edft)
    end if

    if (pub_mermin) then
         call mermin_switch_mu(total_energy, rep%overlap, &
             denskern%kern%m(:,PUB_1K), mu(:), rep%n_occ)
         Ftmp = electronic_lagrangian(total_energy, rep%overlap, &
              denskern%kern%m(:,PUB_1K), &
              ham%ham, mu, rep%n_occ, edft, muext=muext)
    endif

    if (pub_eda .and. (pub_eda_mode == EDA_ISOLATED)) &
       ! mjsp: Store fragment data internally for
       ! mjsp: later supermolecule driver stage
       ! mjsp: (Done in NGWF CG as pur_denskern is required)
       call fragment_data_add_fragment(mdl, pur_denskern%kern, rep, &
            ngwf_basis, pub_frag_counter)

    ! cks: print calculation summary
    call internal_calculation_summary

    ! cks: print quality control information
    if (pub_print_qc) call internal_qc_output

    ! ddor: Write out DFT+U occupancies if it hasn't already been done
    if (pub_hubbard.and.(pub_output_detail .ne. VERBOSE)) then
       !gibo: for cDFT use the cDFT-labelling
       if (pub_cdft) then
          call cdft_energy_info(hub,hub_proj_basis(1))
       else if (pub_dft_nu) then
          call dft_nu_energy_info(hub,hub_proj_basis(1))
       else if (pub_hubbard) then
          call hubbard_energy_info(hub,hub_proj_basis(1))
       endif
    endif

    ! gcc32: exit preconditioning for PAW
    do isub=1,rep%nsub
       call ngwf_gradient_paw_precond_exit(nl_projectors(isub))
    end do

    ! ndmh: for force convergence testing
    if (pub_elec_force_tol > 0.0_DP) then
       deallocate(total_forces,stat=ierr)
       call utils_dealloc_check('ngwf_cg_optimise','total_forces',ierr)
    end if
    deallocate(summary_lines,stat=ierr)
    call utils_dealloc_check('ngwf_cg_optimise','summary_lines',ierr)
    if (pub_ngwf_cg_rotate) then
       if (.not.pub_kernel_christoffel_update) then
          if (pub_aug) call sparse_embed_destroy(start_sp_overlap)
          call sparse_embed_array_destroy(start_denskern)
       end if
       if (pub_edft.and..not.pub_foe) then
          do is = 1, pub_num_spins
             call dense_destroy(start_mo(is))
          end do
          deallocate(start_mo, stat=ierr)
          call utils_dealloc_check('ngwf_cg_optimise', 'start_mo', ierr)
       end if
       if (((pub_elec_cg_max.gt.0).and.(rotate_prev_dir)).or.pub_foe.or.foe_test_rot) &
            call sparse_embed_destroy(start_inv_overlap)
    end if
    if (pub_kernel_christoffel_update) then
       if (pub_aug) call sparse_embed_destroy(start_sp_overlap)
       call sparse_embed_array_destroy(start_denskern)
    end if

    do isub=1,mdl%nsub
       if (pub_ngwf_cg_type == 'NGWF_POLAK') then
          call data_functions_dealloc(prev_contra_grad_on_grid(isub))
       end if
       call data_functions_dealloc(cov_grad_on_grid(isub))
       call data_functions_dealloc(prev_direction_on_grid(isub))
       call data_functions_dealloc(start_ngwfs_on_grid(isub))
       call data_functions_dealloc(direction_on_grid(isub))
       call data_functions_dealloc(contra_grad_on_grid(isub))
    end do
    deallocate(ft_converged, stat=ierr)
    call utils_dealloc_check('ngwf_cg_optimise','ft_converged',ierr)

    call kernel_destroy(pur_denskern)

    call timer_clock('ngwf_cg_optimise',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving ngwf_cg_optimise'

    ! Flush output
    call services_flush

    return

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    real(kind=DP) function internal_step_energy(evaluate_functional)

      !==============================================================!
      ! Update all relevant matrices
      !--------------------------------------------------------------!
      ! Written by Nicholas Hine on 22/08/17                         !
      ! Generalised by Jacek Dziedzic on 10/07/2018 for conduction   !
      ! with HFx.                                                    !
      ! Changes made to make compatible with embedding by Joseph     !
      ! Prentice, summer 2018                                        !
      !==============================================================!

      use constants, only: REP_SWEX_HFX_OTHER
      use rundat, only: pub_active_region, pub_use_activehfx

      logical, intent(in), optional :: evaluate_functional

      logical :: loc_evaluate_functional
      real(kind=DP) :: F

      loc_evaluate_functional = .true.
      if (present(evaluate_functional)) loc_evaluate_functional = evaluate_functional

      ! jcap: loop over regions and check if we need to expand the NGWFs in each region
      do isub=1,mdl%nsub
         ! jd: Re-expand new NGWFs in SWs, if needed
         if(pub_pol_emb_qmstar) then
            call polarisable_embedding_expand_ngwf_pairs(rep, ngwf_basis(isub), &
                 mdl)
         end if
         if(pub_use_hfx.or.(pub_use_activehfx.and.(isub==pub_active_region))) then
            if(.not. pub_cond_calculate) then
               ! jd: Usually we expand valence NGWFs
               call hf_exchange_dkn_indep_stage(hfxstate, mdl, isub, &
                    mdl%regions(isub)%par, &
                    rep, ngwf_basis(isub))
            else
               ! jd: In conduction, we expand mixed NGWF products.
               !     Then rep, ngwf_basis refer to COND.
               call hf_exchange_dkn_indep_stage(hfxstate, mdl, isub, &
                    mdl%regions(isub)%par, &
                    rep, ngwf_basis(isub), &
                    rep2 = val_rep, ngwf_basis2 = val_ngwf_basis(isub), &
                    basis_selector = (/1,1,2,2,REP_SWEX_HFX_OTHER/))
               !                      ^ qoh's alpha-beta-gamma-delta conv'n
               ! basis for beta: 1 (cond),
               ! basis for gamma: 2 (val).
               ! swex selector: HFX_OTHER (of cond_rep)
            end if
         end if
      end do

      ! cks: calculate matrix elements with trial NGWFs
      ! agrecokpt: at specified k-point kpt
      call hamiltonian_dens_indep_matrices(rep, ngwf_basis, proj_basis, &
           nl_projectors, hub_proj_basis, hub, mdl, val_rep, val_ngwf_basis, &
           kpt=loc_kpt_cart)

      ! ars: rotate kernel to the new NGWF representation
      if (pub_ngwf_cg_rotate) then
         if (.not.pub_foe.and..not.foe_test_rot) then
            call internal_transform_kernel(rep%ngwfs_on_grid,&
                 start_ngwfs_on_grid, start_sp_overlap, start_denskern, &
                 old_mo=start_mo)
         else if (pub_foe) then
            call internal_transform_kernel(rep%ngwfs_on_grid,&
                 start_ngwfs_on_grid, start_sp_overlap, start_denskern, &
                 old_inv_overlap=start_inv_overlap)
         else if (foe_test_rot) then
            call internal_transform_kernel(rep%ngwfs_on_grid,&
                 start_ngwfs_on_grid, start_sp_overlap, start_denskern, &
                 old_mo=start_mo, old_inv_overlap=start_inv_overlap)
         end if
      end if

      if (pub_kernel_christoffel_update) call internal_kernel_christoffel( &
           rep%ngwfs_on_grid,start_ngwfs_on_grid,denskern%kern,start_denskern, &
           start_sp_overlap)

      ! lr408: Update conduction Hamiltonian
      if (pub_cond_calculate) call hamiltonian_dens_dep_nonsc(ham,rep, &
           ngwf_basis,mdl, hfxstate, lhxc_fine,hub,val_rep,val_ham, &
           val_dkn%m(:,PUB_1K),updated_shift, &
           val_ngwf_basis = val_ngwf_basis, &      ! jd: Needed for HFx
           cond_dkn = denskern%kern%m(:,PUB_1K))   ! jd: Needed for HFx

      if (pub_on_root.and.pub_kernel_update.and.(pub_output_detail>=NORMAL) &
           .and.loc_evaluate_functional) then
         write(stdout,'(a)') 'Updating kernel for this line step:'
      endif

      ! Now evaluate the energy functional, unless asked not to
      if (.not.loc_evaluate_functional) then
         internal_step_energy = 0.0_DP
      else
         F = electronic_energy(denskern, pur_denskern, ham, &
              lhxc_fine, mu, edft, rep, ngwf_basis, hub_proj_basis, hub, &
              mdl, hfxstate, lnv_threshold, current_maxit_lnv, &
              pub_kernel_update, dfdtau_fine = dfdtau_fine, &
              mermin_threshold=mermin_threshold , &
              current_maxit_mermin=current_maxit_mermin )

         ! ep
         if (pub_mermin .and. pub_kernel_update) then
            call mermin_switch_mu(total_energy, rep%overlap, &
                denskern%kern%m(:,PUB_1K), mu(:), rep%n_occ)
         end if

         internal_step_energy = electronic_lagrangian(F, rep%overlap, &
              denskern%kern%m(:,PUB_1K),ham%ham,mu,rep%n_occ,edft, muext=muext)
      end if

    end function internal_step_energy


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_transform_kernel(new_ngwfs_on_grid,&
         old_ngwfs_on_grid, old_sp_overlap, old_denskern, old_mo, old_inv_overlap)

      !==============================================================!
      ! Transforms the density kernel to its representation in terms !
      ! of new NGWFs                                                 !
      !--------------------------------------------------------------!
      ! Written by Nicholas Hine on 22/10/09                         !
      ! PAW added by David O'Regan 29/3/11                           !
      ! EDFT added by Alvaro Ruiz Serrano on 27/12/2012.             !
      ! Fixed to rotate old_denskern (and not another                !
      ! previously-rotated denskern), Alvaro Ruiz Serrano 27/12/2012 !
      ! Modified for embedding by Joseph Prentice, August 2018       !
      !==============================================================!

      use augmentation, only: augmentation_overlap
      use datatypes, only: FUNCTIONS
      use dense, only: DEM, dense_create, dense_destroy, &
           dense_convert, dense_product
      use function_ops, only: function_ops_brappd_ketppd
      use rundat, only: pub_aug, pub_edft, pub_debug_on_root
      use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
           sparse_embed_create, sparse_embed_transpose, sparse_embed_product, &
           sparse_embed_destroy, sparse_embed_copy, sparse_embed_array_product

      implicit none

      ! Arguments
      type(FUNCTIONS), intent(in) :: new_ngwfs_on_grid(:)
      type(FUNCTIONS), intent(in) :: old_ngwfs_on_grid(:)
      type(SPAM3_EMBED), intent(in) :: old_sp_overlap
      type(SPAM3_EMBED_ARRAY), intent(in) :: old_denskern
      type(DEM), optional, intent(in) :: old_mo(:)
      type(SPAM3_EMBED),optional, intent(in) :: old_inv_overlap

      ! Local Variables
      type(SPAM3_EMBED) :: overlap_old_new
      type(SPAM3_EMBED) :: sk, ks
      type(SPAM3_EMBED_ARRAY) :: ksk
      type(DEM) :: ks_dense
      type(SPAM3_EMBED) :: ham_rot, loc_old_inv_overlap
      integer :: isub,jsub

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering internal_transform_kernel'

      ! Create temporary matrix structures
      call sparse_embed_create(overlap_old_new, rep%overlap)
      call sparse_embed_create(sk, rep%overlap, denskern%kern%m(1,1))
      call sparse_embed_create(ks, denskern%kern%m(1,1), rep%overlap)
      call sparse_embed_array_create(ksk, ks, denskern%kern)

      ! Calculate overlap of old NGWFs with new NGWFs
      ! jcap: loop over regions
      do isub=1,mdl%nsub
         do jsub=1,mdl%nsub
            call function_ops_brappd_ketppd(overlap_old_new%m(isub,jsub), &
                 old_ngwfs_on_grid(isub),ngwf_basis(isub),&
                 new_ngwfs_on_grid(jsub),ngwf_basis(jsub),mdl%cell)

            if (pub_aug) then
               ! ddor: Calculate the augmentation of the overlap matrix
               call augmentation_overlap(overlap_old_new%m(isub,jsub),&
                    mdl%regions(isub)%pseudo_sp,mdl%regions(jsub)%paw_sp, &
                    old_sp_overlap%m(isub,jsub),rep%sp_overlap%m(isub,jsub))
            end if
         end do
      end do

      ! Calculate sk_a^b = <f_a|f'_c>.(S^-1)^cb
      call sparse_embed_product(sk,overlap_old_new,rep%inv_overlap)
      ! Transpose to find ks^a_b = (S^-1)^ac.<f'_c|f_b>
      call sparse_embed_transpose(ks, sk)

      ! Calculate density kernel transformed to new NGWF basis
      ! K'^ah = (S^-1)^ab.<f'_b|f_g> K^ge <f_e|f'_z>.(S^-1)^zh
      call sparse_embed_array_product(ksk, ks, old_denskern)
      call sparse_embed_array_product(denskern%kern, ksk, sk)

      ! ars: rotate edft%mo
      if (pub_edft) then
         ! agrecocmplx
         call dense_create(ks_dense, edft%num, edft%num, iscmplx=loc_cmplx)
         call dense_convert(ks_dense, ks)
         if(present(old_mo)) then
            do is = 1, pub_num_spins
               call dense_product(edft%mo(is), ks_dense, old_mo(is))
            end do
         end if

         if(pub_foe.or.foe_test_rot) then
            if(present(old_inv_overlap)) then
               ! ja531 -> Set up the Hamiltonian rotation matrix.
               call sparse_embed_create(ham_rot,old_inv_overlap,overlap_old_new)
               ! Calculate ham_rot^a_b = (S^-1)^ac.<f_c|f'_b> (left)
               call sparse_embed_product(ham_rot,old_inv_overlap,overlap_old_new)
               ! Transpose to find ham_rot_a^b = <f'_a|f_c>.(S^-1)^cb (right)
            else
               call sparse_embed_create(loc_old_inv_overlap,rep%inv_overlap)

               ! Calculate overlap of old NGWFs with new NGWFs
               ! jcap: loop over regions
               do isub=1,mdl%nsub
                  do jsub=1,mdl%nsub
                     call function_ops_brappd_ketppd(loc_old_inv_overlap%m(isub,jsub), &
                          old_ngwfs_on_grid(isub),ngwf_basis(isub),&
                          old_ngwfs_on_grid(jsub),ngwf_basis(jsub),mdl%cell)

                     if (pub_aug) then
                        call augmentation_overlap(loc_old_inv_overlap%m(isub,jsub), &
                             mdl%regions(isub)%pseudo_sp,mdl%regions(jsub)%paw_sp,&
                             old_sp_overlap%m(isub,jsub),old_sp_overlap%m(isub,jsub))
                     end if
                  end do
               end do

               ! ja531 -> Set up the Hamiltonian rotation matrix.
               call sparse_embed_create(ham_rot,loc_old_inv_overlap,overlap_old_new)
               ! Calculate ham_rot^a_b = (S^-1)^ac.<f_c|f'_b>
               call sparse_embed_product(ham_rot,loc_old_inv_overlap,overlap_old_new)
               call sparse_embed_destroy(loc_old_inv_overlap)
            end if
            ! Transpose to find ham_rot_a^b = <f'_a|f_c>.(S^-1)^cb
            call sparse_embed_copy(edft%rot,ham_rot)
            call sparse_embed_destroy(ham_rot)
         end if

         call dense_destroy(ks_dense)
      end if

      ! Destroy temporary matrix structures
      call sparse_embed_array_destroy(ksk)
      call sparse_embed_destroy(ks)
      call sparse_embed_destroy(sk)
      call sparse_embed_destroy(overlap_old_new)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving internal_transform_kernel'

    end subroutine internal_transform_kernel

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_kernel_christoffel(new_ngwfs_on_grid, &
         old_ngwfs_on_grid,new_denskern,old_denskern,old_sp_overlap)

      !==============================================================!
      ! Updates the density kernel with terms from the Christoffel   !
      ! symbols which preserve the completeness of the basis.        !
      !--------------------------------------------------------------!
      ! Written by David O'Regan in June 2010                        !
      ! Modifications by Nicholas Hine, November 2010.               !
      ! PAW added by David O'Regan 29/3/11                           !
      ! Modified for embedding by Joseph Prentice, August 2018       !
      !==============================================================!

      use augmentation, only: augmentation_overlap
      use datatypes, only: FUNCTIONS, data_functions_alloc, data_functions_dealloc
      use bibliography, only: bibliography_cite
      use function_ops, only: function_ops_brappd_ketppd
      use rundat, only: pub_aug
      use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
           sparse_embed_create, sparse_embed_product, sparse_embed_destroy, &
           sparse_embed_axpy, sparse_embed_array_product, &
           sparse_embed_array_scale, sparse_embed_array_axpy, &
           sparse_embed_array_transpose

      implicit none

      ! Arguments
      type(FUNCTIONS), intent(in) :: new_ngwfs_on_grid(:)
      type(FUNCTIONS), intent(in) :: old_ngwfs_on_grid(:)
      type(SPAM3_EMBED_ARRAY), intent(inout) :: new_denskern
      type(SPAM3_EMBED_ARRAY), intent(in) :: old_denskern
      type(SPAM3_EMBED), intent(in) :: old_sp_overlap

      ! Local Variables
      type(SPAM3_EMBED) :: step, tc_step, step_sp_overlap
      type(SPAM3_EMBED_ARRAY) :: ktmp
      type(FUNCTIONS) :: step_on_grid(mdl%nsub)
      ! agrecocmplx
      logical :: loc_cmplx

      integer :: isub,jsub

      call bibliography_cite('CHRISTOFFEL')

      ! agrecocmplx
      loc_cmplx = old_ngwfs_on_grid(1)%iscmplx

      ! agrecocmplx
      ! jcap: loop over regions
      do isub=1,mdl%nsub
         call data_functions_alloc(step_on_grid(isub), &
              ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
      end do

      ! Create temporary matrix structures
      call sparse_embed_create(step,rep%overlap)
      call sparse_embed_create(tc_step,rep%overlap,rep%inv_overlap)
      call sparse_embed_array_create(ktmp, old_denskern)

      ! Calculate change |dphi> in NGWFs
      ! jcap: loop over regions
      do isub=1,mdl%nsub
         call data_functions_copy(step_on_grid(isub),new_ngwfs_on_grid(isub))
         call data_functions_axpy(step_on_grid(isub),old_ngwfs_on_grid(isub),&
              -1.0_DP)
         !step_on_grid = new_ngwfs_on_grid - old_ngwfs_on_grid
      end do

      ! jcap: create the augmentated overlap operator
      if (pub_aug) then

         ! ndmh: <proj_i|dphi_a> = <proj_i|new_phi_a> - <proj_i|old_phi_a>
         call sparse_embed_create(step_sp_overlap,old_sp_overlap)
         call sparse_embed_copy(step_sp_overlap,rep%sp_overlap)
         call sparse_embed_axpy(step_sp_overlap,old_sp_overlap,-1.0_DP)

      end if

      ! Calculate overlap of the change with the old NGWFs
      ! jcap: loop over regions
      do isub=1,mdl%nsub
         do jsub=1,mdl%nsub
            call function_ops_brappd_ketppd(step%m(isub,jsub), &
                 step_on_grid(isub),ngwf_basis(isub),&
                 old_ngwfs_on_grid(jsub),ngwf_basis(jsub),mdl%cell)

            ! ndmh: Apply the augmentated overlap operator to this matrix
            if (pub_aug) then
               ! ndmh: Calculate the augmentation of the "step" matrix
               call augmentation_overlap(step%m(isub,jsub),&
                    mdl%regions(isub)%pseudo_sp,mdl%regions(jsub)%paw_sp, &
                    old_sp_overlap%m(isub,jsub),step_sp_overlap%m(isub,jsub))
            end if
         end do
      end do

      if (pub_aug) then
         ! ndmh: Clean up <proj_i|dphi_a> (no longer needed)
         call sparse_embed_destroy(step_sp_overlap)
      end if

      ! Raise an index in the step-NGWF overlap
      call sparse_embed_product(tc_step,step,rep%inv_overlap)

      ! Compute -K <g|phi> S^^
      call sparse_embed_array_product(new_denskern, old_denskern, tc_step)
      call sparse_embed_array_scale(new_denskern, -1.0_DP)
      ! Transpose to get -S^^ <phi|g> K
      call sparse_embed_array_transpose(ktmp, new_denskern)
      ! Add this to what we got before
      call sparse_embed_array_axpy(new_denskern, ktmp, 1.0_DP)
      ! Add the previous kernel to the calculated change to get new kernel
      call sparse_embed_array_axpy(new_denskern, old_denskern, 1.0_DP)

      ! Destroy temporary matrix structures
      call sparse_embed_array_destroy(ktmp)
      call sparse_embed_destroy(tc_step)
      call sparse_embed_destroy(step)

      ! jcap: loop over regions
      do isub=1,mdl%nsub
         call data_functions_dealloc(step_on_grid(isub))
      end do

    end subroutine internal_kernel_christoffel

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_hybridize_ngwfs

      !================================================!
      ! Construct atomic hybrids from current NGWFs.   !
      !------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 17/02/2006 !
      !================================================!

      use nho, only: nho_generate
      use projectors, only: projectors_func_ovlp_box
      use rundat, only: pub_any_nl_proj, pub_paw, pub_num_kpoints, &
           PUB_1K, pub_pol_emb_qmstar, pub_use_hfx, pub_eda_scfmi_any, &
           pub_kpoint_method
      use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create, sparse_embed_scale, &
           sparse_embed_axpy, sparse_embed_transpose, sparse_embed_product, &
           sparse_embed_destroy, sparse_embed_array_create, &
           sparse_embed_array_product, sparse_embed_array_destroy
      use utils, only: utils_assert, utils_abort

      implicit none

      type(SPAM3_EMBED) :: ao_to_hy
      type(SPAM3_EMBED) :: trans_ao_to_hy
      type(SPAM3_EMBED) :: hy_to_ao
      type(SPAM3_EMBED) :: trans_hy_to_ao
      type(SPAM3_EMBED_ARRAY) :: k_hao
      type(SPAM3_EMBED) :: h_aoh
      type(SPAM3_EMBED) :: kernel
      integer :: jsub

      ! jme: KPOINTS_DANGER
      ! Some parts of this subroutine enforce a single k-point.
      call utils_assert(pub_num_kpoints == PUB_1K, &
           'Subroutine internal_hybridize_ngwfs (ngwf_cg) not ready yet for more&
           & than one k-point.')

      ao_to_hy%structure ='D'
      ! agrecocmplx
      call sparse_embed_create(ao_to_hy,iscmplx=loc_cmplx)
      trans_ao_to_hy%structure ='D'
      ! agrecocmplx
      call sparse_embed_create(trans_ao_to_hy,iscmplx=loc_cmplx)

      hy_to_ao%structure ='D'
      ! agrecocmplx
      call sparse_embed_create(hy_to_ao,iscmplx=loc_cmplx)
      trans_hy_to_ao%structure ='D'
      ! agrecocmplx
      call sparse_embed_create(trans_hy_to_ao,iscmplx=loc_cmplx)


      ! KPOINTS_DANGER: subroutine only considers 1 k-point
      call sparse_embed_array_create(k_hao, denskern%kern)
      call sparse_embed_create(h_aoh, ham%ham(1))

      ! pdh: take average of up and down kernels for spin polarisation
      call sparse_embed_create(kernel, denskern%kern%m(1,PUB_1K))
      do is=1,pub_num_spins
         call sparse_embed_axpy(kernel, denskern%kern%m(is,PUB_1K), 1.0_DP)
      end do
      call sparse_embed_scale(kernel, 1.0_DP/pub_num_spins)

      ! cks: convert the NGWFs to NHOs
      do isub=1,mrows
         call nho_generate(ao_to_hy%m(isub,isub), hy_to_ao%m(isub,isub), & ! output
              rep%ngwfs_on_grid(isub), prev_direction_on_grid(isub), &     ! in-out
              kernel%m(isub,isub), rep%overlap%m(isub,isub), &
              mdl%regions(isub)%elements, mdl%cell, ngwf_basis(isub), isub)! input
      end do

      call sparse_embed_transpose(trans_ao_to_hy, ao_to_hy)
      call sparse_embed_transpose(trans_hy_to_ao, hy_to_ao)

      ! cks: transform denskern
      call sparse_embed_array_product(k_hao, denskern%kern, trans_hy_to_ao)
      call sparse_embed_array_product(denskern%kern, hy_to_ao, k_hao)
      ! ndmh: mark sk and ks workspaces in kernel mod as invalid
      !call kernel_workspace_invalidate()

      ! cks: transform pur_denskern
      call sparse_embed_array_product(k_hao, pur_denskern%kern, trans_hy_to_ao)
      call sparse_embed_array_product(pur_denskern%kern, hy_to_ao, k_hao)

      ! cks: transform inverse overlap
      call sparse_embed_product(k_hao%m(1,PUB_1K), rep%inv_overlap, trans_hy_to_ao)
      call sparse_embed_product(rep%inv_overlap, hy_to_ao, k_hao%m(1,PUB_1K))

      ! cks: transform overlap
      call sparse_embed_product(h_aoh, rep%overlap, ao_to_hy)
      call sparse_embed_product(rep%overlap, trans_ao_to_hy, h_aoh)

      ! cks: transform kinetic
      call sparse_embed_product(h_aoh, rep%kinet, ao_to_hy)
      call sparse_embed_product(rep%kinet, trans_ao_to_hy, h_aoh)

      ! cks: transform nonlocpot
      if (pub_any_nl_proj) then
         call sparse_embed_product(h_aoh, rep%nonlocpot, ao_to_hy)
         call sparse_embed_product(rep%nonlocpot, trans_ao_to_hy, h_aoh)
      end if
      if (pub_aug) then
         do is=1,pub_num_spins
            call sparse_embed_product(h_aoh, ham%nonlocpot(is), ao_to_hy)
            call sparse_embed_product(ham%nonlocpot(is), trans_ao_to_hy, h_aoh)
         end do
      end if

      ! cks: transform lhxc
      do is=1,pub_num_spins
         call sparse_embed_product(h_aoh, ham%lhxc(is), ao_to_hy)
         call sparse_embed_product(ham%lhxc(is), trans_ao_to_hy, h_aoh)
      end do

      do jsub=1,ncols
         do isub=1,mrows
            ! cks: temporary measure - recompute sp_overlap matrix
            ! agrecokpt: take into account k-point dependence here as well
            ! only needed in KP method
            if (pub_any_nl_proj.or.pub_paw) then
               select case (pub_kpoint_method)
                  ! KP method
                  case('KP')
                     call projectors_func_ovlp_box(rep%sp_overlap%m(isub,jsub), &
                          rep%ngwfs_on_grid(isub),ngwf_basis(isub),proj_basis(jsub),nl_projectors(jsub), &
                          mdl%fftbox,mdl%cell, kshift=loc_kpt)
                  ! TB method: no need to include k-point dependence in projectors
                  case('TB')
                     call projectors_func_ovlp_box(rep%sp_overlap%m(isub,jsub), &
                          rep%ngwfs_on_grid(isub),ngwf_basis(isub),proj_basis(jsub),nl_projectors(jsub), &
                          mdl%fftbox,mdl%cell)
                  case default
                     call utils_abort('Illegal k-point method specified')
               end select
            end if

            if (pub_hubbard) then
               ! ddor: recompute the DFT+U on-site ngwf-projector overlap matrix
               !       since the NGWFs have been modified.
               ! agrecokpt: take into account k-point dependence here as well
               ! need to check it is correct for KP method to have k-point dependent
               ! Hubbard projectors as well?
               select case(pub_kpoint_method)
                  case('KP')
                     ! rc2013: EMBED_FIX!
                     call projectors_func_ovlp_box(rep%hub_overlap%p, & ! out
                          rep%ngwfs_on_grid(1),ngwf_basis(1), &
                          hub_proj_basis(1),hub%projectors, &
                          mdl%fftbox,mdl%cell, kshift=loc_kpt)
                  case('TB')
                     call projectors_func_ovlp_box(rep%hub_overlap%p, & ! out
                          rep%ngwfs_on_grid(1),ngwf_basis(1), &
                          hub_proj_basis(1),hub%projectors, &
                          mdl%fftbox,mdl%cell)
                  case default
                     call utils_abort('Illegal k-point method specified')
               end select
            endif
         end do
      end do

      call utils_assert(.not. pub_use_hfx, 'HFx is not compatible with nnho -- &
           &that would require transforming the X matrix in internal_hybridize_&
           &ngwfs() and some invalidations')

      call utils_assert(.not. pub_pol_emb_qmstar, 'Polarisable embedding, &
           &at least using the QM* representation, is not compatible with nnho &
           &-- that would require transforming the X matrix in internal_hybridi&
           &ze_ngwfs() and some invalidations')

      call utils_assert(.not. pub_eda_scfmi_any, 'EDA is incompatible with &
           &nnho -- that would require transforming EDA overlap matrices in &
           &internal_hybridize_ngwfs()')


      call sparse_embed_destroy(kernel)
      call sparse_embed_destroy(h_aoh)
      call sparse_embed_destroy(ao_to_hy)
      call sparse_embed_destroy(trans_ao_to_hy)
      call sparse_embed_destroy(hy_to_ao)
      call sparse_embed_destroy(trans_hy_to_ao)
      call sparse_embed_array_destroy(k_hao)

    end subroutine internal_hybridize_ngwfs

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_qc_output

      !====================================================!
      ! This subroutine prints out quality control info in !
      ! a form that can be easily accessed and compared.   !
      !----------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 12/03/2005.    !
      !====================================================!

      use comms, only: comms_barrier
      use constants, only: stdout
      use utils, only: utils_qc_print

      if (pub_on_root) then
         write(stdout,'(a1)')' '
         call utils_qc_print('NGWF iterations',iteration)
         call utils_qc_print('total_energy',total_energy)
         call utils_qc_print('rms_gradient',rms_gradient)
      endif

      call comms_barrier


    end subroutine internal_qc_output


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    logical function internal_test_convergence()

      use datatypes, only: data_coef_abs
      use datatypes, only: COEF, data_coef_abs
      use comms, only: comms_reduce
      use forces, only: forces_calculate
      use hubbard_build, only: hubbard_spin_splitting_zero
      use rundat, only: pub_delta_e_conv, pub_hub_ngwf_spin_thr, &
           pub_write_ngwf_grad_plot, pub_elec_energy_tol, &
           pub_ngwf_max_grad, pub_write_forces, maxit_lnv, &
           pub_write_ngwf_grad_radial, pub_kernel_force_conv, cond_maxit_lnv, &
           pub_ngwf_plot_every_it, pub_num_kpoints, PUB_1K, pub_show_overlap, &
           pub_rootname, pub_edft, pub_edft_spin_fix, maxit_mermin
      use sparse, only: sparse_show_matrix
      use sparse_embed, only: sparse_embed_array_is_dense
      use utils, only: utils_unit, utils_sanity_check

      implicit none

      logical :: energy_increasing
      logical :: energy_change_converged
      logical :: ngwf_grad_converged
      logical :: ngwf_max_grad_converged
      logical :: forces_converged
      logical :: dkn_converged !smmd
      logical :: all_converged
      logical :: orig_write_forces
      logical :: test_e_conv, test_f_conv
      logical :: test_dkn_conv !smmd
      logical :: test_rms_grad, test_max_grad
      logical :: maxit_lnv_reached
      logical :: maxit_mermin_reached !ep
      logical :: edft_spin_fix_done
      real(kind=DP) :: max_force_diff
      real(kind=DP) :: max_energy_diff
      real(kind=DP) :: diff12, diff23
      integer :: iat
      integer :: global_iat
      ! agreco: added to write overlap matrix to file
      integer :: fileunit
      integer :: isub, jsub
      real(kind=DP) :: max_gradient_arr(mdl%nsub), tmp
      type(COEF) :: tmp_coef
      integer :: ii

      ! ddor: Used in the case of DFT+U with self-consistent projectors
      character(len=64) :: our_name, hub_proj_iter_number

      ! jme: KPOINTS_DANGER
      ! Some parts of this subroutine enforce a single k-point.
      call utils_assert(pub_num_kpoints == PUB_1K, &
           'Function internal_test_convergence (ngwf_cg) not ready yet for more&
           & than one k-point.')

      ! ndmh: If required, recalculate forces to see if they are converged
      forces_converged = .false.
      test_f_conv = (pub_elec_force_tol > 0.0_DP).and.(iteration>1).and. &
           (.not.pub_cond_calculate).and..not.pub_dmft_spoil_kernel
      if ((pub_elec_force_tol > 0.0_DP).and.(.not.pub_cond_calculate).and..not.pub_dmft_spoil_kernel) then

         ! ndmh: Update storage of most recent three forces
         total_forces(:,:,1) = total_forces(:,:,2)
         total_forces(:,:,2) = total_forces(:,:,3)

         ! ndmh: temporarily override write_forces (we do not want to display
         ! ndmh: the forces here every iteration)
         orig_write_forces = pub_write_forces
         pub_write_forces = .false.

         ! ars: calculate nonsc_forces if required
         if (pub_nonsc_forces) then
            ! jcap: calculate nonsc forces subregion by subregion, as
            ! we only need the diagonal terms anyway
            allocate(sub_nonsc_forces(1:3,mdl%regions(ireg)%par%nat),stat=ierr)
            call utils_alloc_check('energy_and_force_calculate','sub_nonsc_forces', &
                 ierr)
            sub_nonsc_forces=0.d0
            ! jcap: subregion suffix
            suffix = trim(rep%postfix)
            if(mdl%nsub .gt. 1) &
                 write(suffix,'(a,i0,i0)') trim(rep%postfix),ireg,ireg
            call nonsc_forces_ngwfs_calc(sub_nonsc_forces,&
                 rep%ngwfs_on_grid(ireg), contra_grad_on_grid(ireg), &
                 &ngwf_basis(ireg), mdl%cell, mdl%fftbox, &
                 mdl%regions(ireg)%par, suffix=trim(suffix))
            ! jcap: copy data into full array
            do iat=1,mdl%regions(ireg)%par%nat
               ngwf_nonsc_forces(:,mdl%regions(ireg)%elements(iat)%global_atom_number)&
                    = sub_nonsc_forces(:,iat)
            end do
            deallocate(sub_nonsc_forces,stat=ierr)
            call utils_dealloc_check('energy_and_force_calculate','sub_nonsc_forces', &
                 ierr)
         end if

         ! agrecokpt: at specified k-point
         call forces_calculate(total_forces(:,:,3), &
              denskern,ham,lhxc_fine, &
              rep,ngwf_basis,proj_basis,nl_projectors,hub_proj_basis,hub, &
              mdl, hfxstate, ngwf_nonsc_forces,kpt=loc_kpt, &
              dfdtau_fine=dfdtau_fine)
         pub_write_forces = orig_write_forces

         diff12 = 0.0_DP
         diff23 = 0.0_DP
         ! Loop over atoms to find greatest change in total_forces
         do iat=1,mdl%regions(ireg)%par%nat
            global_iat=mdl%regions(ireg)%elements(iat)%global_atom_number
            diff12 = max(diff12, sqrt(sum((total_forces(:,global_iat,1) - &
                 total_forces(:,global_iat,2))**2)))
            diff23 = max(diff23, sqrt(sum((total_forces(:,global_iat,2) - &
                 total_forces(:,global_iat,3))**2)))
         end do
         if (iteration<2) diff12 = 0.0_DP
         max_force_diff = max(diff12,diff23)

         if ((diff23<pub_elec_force_tol).and.(diff12<pub_elec_force_tol)) then
            forces_converged = .true.
         end if

         ! tjz21: bugfix provided by Edward Linscott.
         ! forces_converged needs to be set consistently across all procs. If
         ! one proc reports forces_converged=.false., set variable to .false.
         ! accross all procs.
         call comms_reduce("AND",forces_converged)

      end if

      ! cks: Update storage of most recent three energies
      last_n_energies(1) = last_n_energies(2)
      last_n_energies(2) = last_n_energies(3)
      last_n_energies(3) = total_energy

      ! ndmh: Check energy change per atom vs tolerance
      energy_change_converged = .false.
      test_e_conv = (pub_elec_energy_tol > 0.0_DP).and.(iteration>1)
      if (test_e_conv) then
         ! mjsp: Modified to avoid bug in which "Maximum change in energy" is
         ! 1E+306 on the first iteration.
         diff12 = abs(last_n_energies(1)-last_n_energies(2)) &
              / real(mdl%nat,kind=DP)
         diff23 = abs(last_n_energies(2)-last_n_energies(3)) &
              / real(mdl%nat,kind=DP)
         if (iteration==2) diff12 = 0.0_DP
         max_energy_diff = max(diff12,diff23)
         if (max_energy_diff<pub_elec_energy_tol) then
            energy_change_converged = .true.
         end if
      end if

      ! smmd: Check convergence of density kernel
      test_dkn_conv = pub_kernel_force_conv
      if (test_dkn_conv) then
         if (pub_cond_calculate) then
            maxit_lnv_reached = (current_maxit_lnv==cond_maxit_lnv)
         else
           if (pub_mermin) then
               maxit_mermin_reached = (current_maxit_mermin==maxit_mermin)
            else
               maxit_lnv_reached = (current_maxit_lnv==maxit_lnv)
            end if
         end if
         if ((cvstat.gt.0 .and. (.not.maxit_lnv_reached)) .or. &
            (cvstat.gt.0 .and. (.not.maxit_mermin_reached))) then
            dkn_converged = .false.
         else
            dkn_converged = .true.
         endif
      endif

      ! calculate rms gradient
      ! cks: parallel calculation of NGWF rms_gradient
      call data_set_to_zero(current_g_dot_g)
      current_g_dot_g = &
           data_functions_dot(contra_grad_on_grid, cov_grad_on_grid)
      if (current_g_dot_g%iscmplx) then
         call comms_reduce('SUM', current_g_dot_g%z)
      else
         call comms_reduce('SUM', current_g_dot_g%d)
      end if

      ! cks: store for Fletcher-Reeves calculation of CG coeff
      rms_gradient = sqrt(data_coef_abs(current_g_dot_g) / real(est_num_psincs, kind=DP))
      ! rc2013: print some debugging info -- can probably remove
      if(pub_debug_on_root) write(stdout,'(a,f15.9)') 'DEBUG: g_dot_g = ', &
           data_coef_abs(current_g_dot_g)
      if(pub_debug_on_root) write(stdout,'(a,f20.9)') 'DEBUG: est_num_psincs = ', &
           real(est_num_psincs, kind=DP)
      if(pub_debug_on_root) write(stdout,'(a,f15.9)') 'DEBUG: rms_gradient = ', &
           rms_gradient

      ! ars: calculate inner product of the gradient
      if (calculate_inner_grad) then
         ! rc2013: EMBED_FIX
         inner_grad = internal_inner_grad(contra_grad_on_grid(ireg), &
              rep%overlap%m(ireg,mdl%regions(ireg)%reg_count), ngwf_basis(ireg))
      end if

      ! ndmh: Calculate max value of <g^a|g_a> over all NGWFs on all procs
      max_gradient = 0.0_DP
      test_max_grad = (pub_ngwf_max_grad > 0.0_DP)
      if (test_max_grad) then
         do isub=1,mdl%nsub
            ! ndmh_pointerfun
            ! agrecocmplx: take conjugate of left hand side in complex case
            if (loc_cmplx) then
                max_gradient_arr = maxval(abs( &
                    conjg(contra_grad_on_grid(isub)%z) * cov_grad_on_grid(isub)%z))
            else
                max_gradient_arr = maxval(abs( &
                    contra_grad_on_grid(isub)%d * cov_grad_on_grid(isub)%d))
            end if
         end do
         max_gradient = maxval(max_gradient_arr)
         call comms_reduce('MAX', max_gradient)
         max_gradient = sqrt(max_gradient)
      end if

      ! cks: Test for if allowed to use energy gain as convergence criterion
      energy_increasing = .false.
      energy_increasing = pub_delta_e_conv .and. &
           ( last_n_energies(3) > last_n_energies(2)) .and. &
           ( last_n_energies(2) > last_n_energies(1))

      ! gibo: for cDFT-optimisation, do not care if energy is going up
      !       for single-point cDFT [Fixed cDFT-potential], do care...
      !gom
      if ((pub_cdft.OR.pub_dft_nu).AND.(pub_maxit_cdft_u_cg>0)) energy_increasing=.false.

      !ep: TO LEAVE UNTIL PROBLEM WITH LOW ORDER CHEBYSHEV SOLVED
      !      IN PRINCIPLE SHOUD NOT BE HERE
      if ((pub_mermin)) energy_increasing=.false.
      !ep

      ! cks: NGWF RMS gradient convergence criterion
      test_rms_grad = (ngwf_threshold > 0.0_DP)
      if (test_rms_grad) then
         ngwf_grad_converged = (rms_gradient < ngwf_threshold)
      else
         ngwf_grad_converged = .false.
      end if

      ! ndmh: NGWF MAX gradient convergence criterion
      ngwf_max_grad_converged = .false.
      if (test_max_grad) then
         ngwf_max_grad_converged = (max_gradient < pub_ngwf_max_grad)
      end if

      ! kkbd: EDFT spin fix criterion: ban convergence if we're holding the spin
      ! fixed for the first few iterations of a free-spin run.
      edft_spin_fix_done = .true.
      if (pub_edft .and. (pub_edft_spin_fix > 0)) then
         edft_spin_fix_done = .false.
      end if


      ! ndmh: Check all relevant criteria for convergence
      all_converged = .true.
      if ((pub_elec_force_tol > 0.0_DP).and.(iteration<2)) &
           all_converged = .false.
      if ((pub_elec_energy_tol > 0.0_DP).and.(iteration<2)) &
           all_converged = .false.
      if (test_f_conv) &
           all_converged = all_converged.and.forces_converged
      if (test_e_conv) &
           all_converged = all_converged.and.energy_change_converged
      if (test_rms_grad) &
           all_converged = all_converged.and.ngwf_grad_converged
      if (test_max_grad) &
           all_converged = all_converged.and.ngwf_max_grad_converged
      if (pub_delta_e_conv) &
           all_converged = all_converged.or.energy_increasing
      if (test_dkn_conv) then
           if (all_converged .and. .not. dkn_converged) then
              if (pub_cond_calculate) then
                 current_maxit_lnv = cond_maxit_lnv
              else
                 if (pub_mermin) then
                    current_maxit_mermin = maxit_mermin
                 else
                    current_maxit_lnv = maxit_lnv
                 end if
              end if
           end if
           all_converged = all_converged.and.dkn_converged
      endif
      if (pub_edft) &
         all_converged = all_converged.and.edft_spin_fix_done

      ! gibo: cDFT-optimisation additions
      !      [cDFT-opt. is different from single-point, fixed potential cDFT]
      !gom
      if ((pub_cdft .or. pub_dft_nu).AND. (pub_maxit_cdft_u_cg>0)) &
           all_converged = all_converged.and.cdft_converged
      ! gibo: if cdft_tight=T, perform extra NGWFs steps with
      !       10-times tigther cDFT-optimisation criteria
      if (pub_cdft.and.pub_cdft_tight.and.all_converged.and.&
          (.not.cdft_coarse_converged)) then
              cdft_coarse_converged = cdft_converged
              cdft_converged = .false.
              all_converged = .false.
              pub_cdft_cg_threshold = 0.1_DP*pub_cdft_cg_threshold
              pub_cdft_max_grad = 0.1_DP*pub_cdft_max_grad
              !=== (possibly redundant) make all procs aware of the changes
              call comms_bcast(pub_root_proc_id, pub_cdft_max_grad)
              call comms_bcast(pub_root_proc_id, pub_cdft_cg_threshold)
              !=== (possibly redundant) make all procs aware of the changes
      endif

      if (all_converged) then

         ! cks: print details only when output is not set to brief
         if (pub_output_detail >= NORMAL) then

            if (pub_on_root) then
               write(stdout,'(/a)')'           ............................&
                    &............................'
               write(stdout,'(a)')'           |         *** NGWF optimisation &
                    &converged ***          |'
               if (test_rms_grad) write(stdout,'(a,f19.14,a)') &
                    '           | RMS NGWF gradient = ',rms_gradient,'              |'
               if (test_max_grad) write(stdout,'(a,f19.14,a)') &
                    '           | MAX NGWF gradient = ',max_gradient,'              |'
               if (test_f_conv) write(stdout,'(a,f17.14,a)') &
                    '           | Maximum change in forces = ',max_force_diff,' au      |'
               if (test_e_conv) write(stdout,'(a,f17.14,a)') &
                    '           | Maximum change in energy = ',max_energy_diff,' / atom  |'
               write(stdout,'(a)')'           | Criteria satisfied: &
                    &                                 |'
               if (test_rms_grad.and.ngwf_grad_converged) write(stdout,'(a)') &
                    '           | -> RMS NGWF gradient lower than set threshold.       |'
               if (test_max_grad.and.ngwf_max_grad_converged) write(stdout,'(a)') &
                    '           | -> MAX NGWF gradient lower than set threshold.       |'
               if (test_e_conv.and.energy_change_converged) write(stdout,'(a)') &
                    '           | -> Energy change per atom lower than set threshold.  |'
               if (test_f_conv.and.forces_converged) write(stdout,'(a)') &
                    '           | -> Maximum force change lower than set threshold.    |'
               if (test_dkn_conv.and.cvstat.le.0) write(stdout,'(a)') &
                    '           | -> Density kernel threshold reached                  |'
               if (test_dkn_conv.and.cvstat.gt.0 .and. .not. pub_mermin) then
                  write(stdout,'(a)') &
                    '           | -> Density kernel unconverged but                    |'
                  write(stdout,'(a)') &
                    '           |             # of lnv iterations = maxit_lnv          |'
               endif
              if (test_dkn_conv.and.cvstat.gt.0 .and. pub_mermin) then
                  write(stdout,'(a)') &
                    '           | -> Density kernel unconverged but                    |'
                  write(stdout,'(a)') &
                    '           |             # of mermin iterations = maxit_mermin    |'
               endif
               if (energy_increasing) then
                  write(stdout,'(a)') &
                       '           | -> Maximum degree of convergence for applied level   |'
                  write(stdout,'(a)') &
                       '           |    of density kernel truncation has been reached.    | '
               endif
               ! rc2013: if this is a freeze-and-thaw calc then we're not done yet...
               if(pub_do_fandt) then
                  write(stdout,'(a)') &
                       '           |    Proceeding with freeze-and-thaw calculation:      | '
                  write(stdout,'(a)') &
                       '           |        swapping regions.                             | '
               endif
               ! rc2013: nor are we done if it's time to swap XC functionals!
               if(pub_emft_follow .and. .not. fnl_switched) then
                  write(stdout,'(a)') &
                       '           |    Proceeding with embedding calculation:            | '
                  write(stdout,'(a)') &
                       '           |        swapping exchange-correlation functionals.    | '
               endif
               write(stdout,'(a)') '           ===========================&
                    &============================='
            end if
         end if

         ! ndmh: if there is no kernel truncation but energy gain convergence
         ! ndmh: has been set due to rising energy, something must be wrong
         ! ndmh: with the kernel (unless unreasonable demands on NGWF
         ! ndmh: convergence have been requested), so warn user
         if (pub_on_root.and.energy_increasing.and. &
              (sparse_embed_array_is_dense(denskern%kern))) then
            write(stdout,'(a,f8.4,a)')  'WARNING: No kernel truncation, &
                 &yet energy has risen on last 2 iterations.'
            write(stdout,'(a,f8.4,a)')  'WARNING: This most likely indicates &
                 &a problem with the density kernel.'
            write(stdout,'(a)') 'WARNING: Either kernel occupation numbers &
                 &may be unreliable, or'
            write(stdout,'(a)') 'WARNING: requested NGWF gradient tolerance &
                 &may be unachievable.'
         end if

         if (.not. (pub_cond_calculate .or. pub_do_fandt)) then

            if (.not.(pub_eda_scfmi) .and. .not.(pub_mermin)) then
               call hamiltonian_energy_components(  &
                    pur_denskern%kern%m(:,PUB_1K), rep, mdl, &
                    ngwf_basis, hub_proj_basis, hub, ham%hfexchange, edft=edft)
            !ep: Print energy components with mermin args
            else if (pub_mermin) then
               call hamiltonian_energy_components(  &
                    denskern%kern%m(:,PUB_1K), rep, mdl, &
                    ngwf_basis, hub_proj_basis, hub, ham%hfexchange, &
                    muext(:,PUB_1K))
            else
               ! mjsp: if SCF-MI then use the supermolecule density kernel
               ! mjsp: representation
               call hamiltonian_energy_components(  &
                    denskern_R%kern%m(:,PUB_1K), rep, mdl, &
                    ngwf_basis, hub_proj_basis, hub, ham%hfexchange)
            end if

         end if

         ! cks: write final NGWFs in plotting formats
         if (pub_write_ngwf_plot .and. &
              ( pub_cube_format .or. pub_grd_format .or. pub_dx_format) .and. &
              (.not.(pub_eda))) then

            ! rc2013: write out all regions
            do isub=1,mdl%nsub
               write(isub_str,'(i1)') isub
               ! ddor: If carrying out a DFT+U calculation with self-consistent
               ! ddor: projectors, then write out the final NGWFs on each projector
               ! ddor: optimisation iteration.
               if ( pub_task == 'HUBBARDSCF' ) then
                  write(hub_proj_iter_number,*) hub%consistency_iteration
                  if (mdl%nsub.gt.1) then
                     write(our_name,'(a64)') 'hub_consistency_iteration'// &
                          trim(adjustl(hub_proj_iter_number))//'_'//isub_str
                  else
                     write(our_name,'(a64)') 'hub_consistency_iteration'// &
                          trim(adjustl(hub_proj_iter_number))
                  end if
                  call visual_ngwfs(rep%ngwfs_on_grid(isub), &
                       ngwf_basis(isub),our_name, &
                       mdl%regions(isub)%elements,mdl%cell,mdl%fftbox, &
                       mdl%regions(isub)%par)
               else
                  if (mdl%nsub.gt.1) then
                     call visual_ngwfs(rep%ngwfs_on_grid(isub), &
                          ngwf_basis(isub),'final_'//isub_str, &
                          mdl%regions(isub)%elements,mdl%cell,mdl%fftbox, &
                          mdl%regions(isub)%par)
                  else
                     call visual_ngwfs(rep%ngwfs_on_grid(isub), &
                          ngwf_basis(isub),'final', &
                          mdl%regions(isub)%elements,mdl%cell,mdl%fftbox, &
                          mdl%regions(isub)%par)
                  end if
               endif
            end do

         end if

         ! smmd: Write NGWFs radial distribution for final iteration
         if (pub_write_ngwf_radial.ge.1) then
            ! rc2013: write out all regions to separate files
            do isub=1,mdl%nsub
               if (mdl%nsub.gt.1) then
                  write(isub_str,'(i1)') isub
                  call visual_ngwfs_radial(rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                       'final_'//isub_str, 'ngwf',mdl%fftbox,mdl%cell, &
                       mdl%regions(isub)%par)
               else
                  call visual_ngwfs_radial(rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                       'final', 'ngwf',mdl%fftbox,mdl%cell, &
                       mdl%regions(isub)%par)
               end if
            end do
         end if

         ! smmd: Write covariant gradient radial distribution for final iteration
         if (pub_write_ngwf_grad_radial.ge.1) then
            ! rc2013: write out all regions to separate files
            do isub=1,mdl%nsub
               if (mdl%nsub.gt.1) then
                  write(isub_str,'(i1)') isub
                  call visual_ngwfs_radial(cov_grad_on_grid(isub), ngwf_basis(isub), &
                       'final_'//isub_str, 'cov_grad',mdl%fftbox,mdl%cell, &
                       mdl%regions(isub)%par)
                  call visual_ngwfs_radial(contra_grad_on_grid(isub), ngwf_basis(isub), &
                       'final_'//isub_str, 'contra_grad',mdl%fftbox,mdl%cell, & ! ars
                       mdl%regions(isub)%par)
               else
                  call visual_ngwfs_radial(cov_grad_on_grid(isub), ngwf_basis(isub), &
                       'final', 'cov_grad',mdl%fftbox,mdl%cell, &
                       mdl%regions(isub)%par)
                  call visual_ngwfs_radial(contra_grad_on_grid(isub), ngwf_basis(isub), &
                       'final', 'contra_grad',mdl%fftbox,mdl%cell, & ! ars
                       mdl%regions(isub)%par)
               end if
            end do
         end if

         ! agreco: Write overlap matrix for final iteration
         if (pub_show_overlap) then
            if(pub_on_root) then
               fileunit = utils_unit()
               open(unit=fileunit, file=trim(pub_rootname)//'.sparse_overlap')
            end if

            do jsub=1,ncols
               do isub=1,mrows
                  call sparse_show_matrix(rep%overlap%m(isub,jsub),fileunit, &
                       matlab_format=.true.)
               enddo
            enddo

            if(pub_on_root) then
               close(fileunit)
            end if

         end if

         internal_test_convergence = .true.
      else

         ! print un-converged values on NORMAL or above
         if (pub_output_detail >= NORMAL .and. pub_on_root) then
            write(stdout,'(2x,a)') 'NOT CONVERGED'

            write(stdout,'(2x,a)') repeat('=',72)

            if (test_rms_grad) call internal_print_param_status( &
                 'NGWF RMS gradient', rms_gradient, ngwf_threshold_orig)
            if (test_max_grad) call internal_print_param_status( &
                 'MAX NGWF gradient', max_gradient, pub_ngwf_max_grad)
            if (test_e_conv) call internal_print_param_status( &
                 'Maximum change in energy', max_energy_diff, pub_elec_energy_tol)
            if (test_f_conv) call internal_print_param_status( &
                 'Maximum change in forces', max_force_diff,pub_elec_force_tol)
            if (test_dkn_conv.and. (cvstat.gt.0)) &
               write(stdout,'(2x,a)') 'Density kernel threshold unsatisfied'
            if (.not.edft_spin_fix_done) then
               write(stdout,*) "Net spin fixed for",&
                        pub_edft_spin_fix,"more iterations."
            end if
            if (energy_increasing) then
               write(stdout,*)
               write(stdout,'(a)') 'Energy is increasing: Maximum degree of &
                    &convergence may have been reached'
               write(stdout,'(a)') 'for applied level of kernel truncation.'
            end if
            write(stdout,'(2x,a/)') repeat('=',72)
         endif

         internal_test_convergence = .false.
      end if

      ! cks: set number of lnv iterations for next NGWF step
      if (rms_gradient > 0.9_DP*previous_rms_gradient) then
         if (pub_mermin) then
             current_maxit_mermin = current_maxit_mermin+1
         else
             current_maxit_lnv = current_maxit_lnv +1
         end if
      end if

      if (pub_cond_calculate) then
         if (current_maxit_lnv > cond_maxit_lnv) &
              current_maxit_lnv = cond_maxit_lnv
      else
        if (pub_mermin) then
            if (current_maxit_mermin > maxit_mermin) &
                 current_maxit_mermin = maxit_mermin
         else
            if (current_maxit_lnv > maxit_lnv) &
                 current_maxit_lnv = maxit_lnv
         end if
      end if
      ! cks: store current rms grad for next NGWF step
      previous_rms_gradient = rms_gradient

      if (.not.all_converged) then
         !ddor: Return DFT+U spin-splitting to zero if we have
         !      proceeded far enough in the calculation.
         if ((pub_hubbard) .and. (rms_gradient .lt. pub_hub_ngwf_spin_thr)) then
            call hubbard_spin_splitting_zero(hub,hub_proj_basis(1))
         endif
      end if

      ! cks: Write NGWFs of current iteration in plotting formats
      ! agreco: write at every iteration only if requested
      if ((pub_write_ngwf_plot .and. pub_ngwf_plot_every_it) .and. &
           (pub_cube_format .or. pub_grd_format .or. pub_dx_format) .and. &
           (.not.(pub_eda))) then
         ! rc2013: write all regions to separate files
         do isub=1,mdl%nsub
            write(isub_str,'(i1)') isub
            write(iter_number,*) iteration
            write(fun_name,'(a64)') 'iteration'// &
                 trim(adjustl(iter_number))
            if (mdl%nsub.gt.1) then
               write(isub_str,'(i1)') isub
               write(sub_fun_name,'(a64)') trim(fun_name)//'_'//isub_str
            else
               sub_fun_name = fun_name
            end if
            call visual_ngwfs(rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                 sub_fun_name, mdl%regions(isub)%elements, &
                 mdl%cell, mdl%fftbox, mdl%regions(isub)%par)
         end do
      end if

      ! ndmh: write gradient on grid to file for visualisation
      if (pub_write_ngwf_grad_plot .and. &
           (pub_cube_format .or. pub_grd_format .or. pub_dx_format)) then
         ! rc2013: write all regions to separate files
         do isub=1,mdl%nsub
            if (mdl%nsub.gt.1) then
               write(isub_str,'(i1)') isub
               call visual_ngwfs(cov_grad_on_grid(isub), ngwf_basis(isub), &
                    'ngwf_cov_grad_'//isub_str, mdl%regions(isub)%elements, &
                    mdl%cell, mdl%fftbox, mdl%regions(isub)%par)
               call visual_ngwfs(contra_grad_on_grid(isub), ngwf_basis(isub), &
                    'ngwf_contra_grad_'//isub_str, mdl%regions(isub)%elements, &
                    mdl%cell, mdl%fftbox, mdl%regions(isub)%par)
            else
               call visual_ngwfs(cov_grad_on_grid(isub), ngwf_basis(isub), &
                    'ngwf_cov_grad', mdl%regions(isub)%elements, &
                    mdl%cell, mdl%fftbox, mdl%regions(isub)%par)
               call visual_ngwfs(contra_grad_on_grid(isub), ngwf_basis(isub), &
                    'ngwf_contra_grad', mdl%regions(isub)%elements, &
                    mdl%cell, mdl%fftbox, mdl%regions(isub)%par)
            end if
         end do
      end if

      ! smmd: Write NGWFs radial distribution for current iteration
      if (pub_write_ngwf_radial.ge.2) then
         write(iter_number,*) iteration
         write(fun_name,'(a64)') 'iter'//trim(adjustl(iter_number))
         ! rc2013: write all regions to separate files
         do isub=1,mdl%nsub
            if (mdl%nsub.gt.1) then
               write(isub_str,'(i1)') isub
               write(sub_fun_name,'(a64)') trim(fun_name)//'_'//isub_str
            else
               sub_fun_name = fun_name
            end if
            call visual_ngwfs_radial(rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                 sub_fun_name, 'ngwf',mdl%fftbox,mdl%cell, &
                 mdl%regions(isub)%par)
         end do
      end if

      ! smmd: Write covariant gradient radial distribution for current iteration
      if (pub_write_ngwf_grad_radial.ge.2) then
         write(iter_number,*) iteration
         write(fun_name,'(a64)') 'iter'//trim(adjustl(iter_number))
         ! rc2013: write all regions to separate files
         do isub=1,mdl%nsub
            if (mdl%nsub.gt.1) then
               write(isub_str,'(i1)') isub
               write(sub_fun_name,'(a64)') trim(fun_name)//'_'//isub_str
            else
               sub_fun_name = fun_name
            end if
            call visual_ngwfs_radial(cov_grad_on_grid(isub), ngwf_basis(isub), &
                 sub_fun_name, 'cov_grad',mdl%fftbox,mdl%cell, &
                 mdl%regions(isub)%par)
         end do
      end if

      call utils_sanity_check(rms_gradient,'NGWF RMS gradient',excessive=1D5)

    end function internal_test_convergence

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_print_param_status(name,value,tol)

       implicit none

       ! arguments
       character(len=*), intent(in) :: name
       real(kind=DP), intent(in)    :: value
       real(kind=DP), intent(in)    :: tol

       ! internal
       logical :: param_converged

       param_converged = (value < tol)

       if (param_converged) then
          write(stdout,'(a26,a,es11.4,a,es11.4,a)') &
               adjustl(name),' = ', value, ' < ', tol, '  | CONVERGED'
       else
          write(stdout,'(a26,a,es11.4,a,es11.4,a)') &
               adjustl(name),' = ', value, ' > ', tol, '  | above tolerance'
       endif

    end subroutine internal_print_param_status

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! cks: print a concise summary of the calculation
    subroutine internal_calculation_summary

      use utils, only: utils_banner
      implicit none

      integer :: sumrow
      integer :: lastrow

      lastrow = max(iteration,1)

      if (pub_on_root) then
         ! ndmh: adapt number of digits in total energy to suit its magnitude
         if(abs(total_energy)<100000.0_DP) then
            if(pub_do_fandt) then
               write(summary_lines(lastrow),'(i4, a8, f21.14, f22.14, a)') &
                    iteration-1, '      ', rms_gradient, total_energy, '  <-- CG'
            else
               write(summary_lines(lastrow),'(i4, f21.14, f22.14, a)') &
                    iteration-1, rms_gradient, total_energy, '  <-- CG'
            end if
         else
            if(pub_do_fandt) then
               write(summary_lines(lastrow),'(i4, a8, f21.14, f22.12, a)') &
                    iteration-1, '      ', rms_gradient, total_energy, '  <-- CG'
            else
               write(summary_lines(lastrow),'(i4, f21.14, f22.12, a)') &
                    iteration-1, rms_gradient, total_energy, '  <-- CG'
            end if
         end if

         if (pub_output_detail >= NORMAL) then
            write(stdout,'(/a)') utils_banner(' ', &
                 '<<<<< CALCULATION SUMMARY >>>>>')
            do sumrow=0,lastrow
               write(stdout,'(a88)') summary_lines(sumrow)
            end do
         else
            write(stdout,'(a88)') summary_lines(iteration)
         endif

      end if

    end subroutine internal_calculation_summary

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_find_direction

      use constants, only: cmplx_0
      use datatypes, only: data_coef_abs, data_functions_axpy, &
           data_set_to_zero, data_functions_scale, data_coef_scale
      use services, only: services_polak_cg_coeff

      implicit none

      ! agrecocmplx
      logical :: slope_is_positive
      ! rc2013: placeholder for slope
      type(COEF) :: tmp_G

      tmp_G%iscmplx = loc_cmplx
      call data_set_to_zero(tmp_G)

      if (.not.line_search_success) then
         trial_length = trial_length * 0.5_DP
      else if (line_search_coeff > 0.0_DP) then
         trial_length = &
              max(sqrt(trial_length * line_search_coeff),epsilon(1.0_DP))
      end if
      trial_length = max(trial_length,0.0001_DP)

      if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then
         write(stdout,'(a)') '>>> Improving NGWFs using line search:'
      endif


      ! calculate CG coefficient
      if ( iteration > 1 ) then
         if ((cg_count >= pub_elec_cg_max) .or. (.not.line_search_success) &
              .or. (updated_shift)) then
            ! cks: reset CG after "cg_max" steps
            ! ndmh: or after a fitting failure
            ! lr408: or if conduction shift has just been updated
            ! rc2013: or if we've just switched regions
            call data_set_to_zero(cg_coeff)
            cg_count = 0
            if (pub_on_root  .and. (.not.line_search_success) .and. &
                 (pub_output_detail >= NORMAL) ) write(stdout,'(a)') &
                 'Resetting NGWF CG'
         else
            ! cks <<< FLETCHER >>>
            if (pub_ngwf_cg_type == 'NGWF_FLETCHER') then

               ! cks: original Fletcher-Reeves formula (cheaper in memory)
               if (data_coef_abs(previous_g_dot_g) > epsilon(1.0_DP)) then
                  ! agrecocmplx: in complex case, take ratio of abs values
                  if (loc_cmplx) then
                     cg_coeff%z = cmplx(abs(current_g_dot_g%z)/abs(previous_g_dot_g%z),&
                             0.0_DP, kind=DP)
                  else
                     cg_coeff%d = current_g_dot_g%d / previous_g_dot_g%d
                  end if
                  ! cks: protection from crazy coefficients
                  if (data_coef_abs(cg_coeff) > 2.0_DP ) then                 !jmecmplx
                     if (pub_on_root  .and.(pub_output_detail >= VERBOSE) ) then
                        ! agrecocmplx
                        if (loc_cmplx) then
                           write(stdout,'(a,f8.4,1x,f8.4,a)') 'WARNING: NGWF&
                                & Fletcher-Reeves CG coeff too large (',  &
                                cg_coeff%z, ') - setting to zero'
                        else
                           write(stdout,'(a,f8.4,a)') 'WARNING: NGWF&
                                & Fletcher-Reeves CG coeff too large (',  &
                                cg_coeff%d, ') - setting to zero'
                        end if
                     end if
                     call data_set_to_zero(cg_coeff)                          !jmecmplx
                     cg_count = 0
                  end if
               else
                  ! agrecocmplx
                  if (pub_on_root .and.(pub_output_detail >= NORMAL) ) then
                       if (loc_cmplx) then
                          write(stdout,*) ' previous_g_dot_g=', &
                          previous_g_dot_g%z, &
                          'CG coeffient set to zero in ngwf_cg_optimise'

                       else
                          write(stdout,*) ' previous_g_dot_g=', &
                          previous_g_dot_g%d, &
                          'CG coeffient set to zero in ngwf_cg_optimise'
                       end if
                  end if
                  call data_set_to_zero(cg_coeff)
                  cg_count = 0
               end if

               ! cks: <<< POLAK >>>
            else if (pub_ngwf_cg_type == 'NGWF_POLAK') then

               ! POLAK FORMULA
               do isub=1,mdl%nsub
                  cg_coeff = services_polak_cg_coeff(prev_direction_on_grid(isub), &
                       cov_grad_on_grid(isub),contra_grad_on_grid(isub), &
                       prev_contra_grad_on_grid(isub), &
                       ngwf_basis(isub)%size_on_grid)
               end do

               ! FLETCHER-REEVES FORMULA
               !        cg_coeff=services_fr2_cg_coeff(prev_direction_on_grid, &
               !             cov_grad_on_grid,contra_grad_on_grid, &
               !             prev_contra_grad_on_grid,n_ngwf_ppds*cell%n_pts)

            end if


            ! cks: re-initialise the periodic reset process if cg_coeff was zero
            ! cks: otherwise increase cg_count
            if(loc_cmplx) then
               if(cg_coeff%z == cmplx_0) then
                  cg_count = 0
               else
                  cg_count = cg_count +1
               end if
            else
               if(cg_coeff%d == 0.0_DP) then
                  cg_count = 0
               else
                  cg_count = cg_count + 1
               end if
            end if

         end if
      else
         call data_set_to_zero(cg_coeff)
      endif

      ! cks: HACK
      !      cg_coeff = 0.0_DP      ! STEEPEST DESCENTS ONLY
      !      write(stdout,*)'Steepest descents enforced!'

      do isub=1,mdl%nsub
        ! Find search direction
        if (pub_elec_cg_max > 0) then
            ! ndmh_pointerfun
            call data_functions_copy(direction_on_grid(isub),cov_grad_on_grid(isub))
            call data_functions_scale(direction_on_grid(isub),-1.0_DP)
            call data_functions_axpy(direction_on_grid(isub),prev_direction_on_grid(isub),cg_coeff)
            !direction_on_grid = cg_coeff * prev_direction_on_grid - cov_grad_on_grid
        else
            ! cks:steepest descents
            call data_set_to_zero(direction_on_grid(isub))
            call data_functions_axpy(direction_on_grid(isub), cov_grad_on_grid(isub), -1.0_DP)
        endif


      end do
      ! Slope of energy in search direction
      G_init = data_functions_dot(contra_grad_on_grid, direction_on_grid)

      ! cks: collect the work of each proc
      if (loc_cmplx) then
         call comms_reduce('SUM', G_init%z)
      else
         call comms_reduce('SUM', G_init%d)
      end if

      ! take action in case of positive slope along search direction
      slope_is_positive = .false.
      ! agrecocmplx: consider sign of real part in complex case?
      if (loc_cmplx) then
         G_init_real = real(G_init%z,kind=DP)
      else
         G_init_real = G_init%d
      end if

      if (G_init_real > 0.0_DP) &
            slope_is_positive = .true.

      ! agrecocmplx
      if (slope_is_positive) then                  !jmecmplx
         if (pub_on_root .and. (pub_output_detail >= VERBOSE) ) then
            if (loc_cmplx) then
               write(stdout,'(a,2(e16.6))') &                                      !jmecmplx
                    &'WARNING: slope along search direction is positive:', G_init%z !jmecmplx
            else
               write(stdout,'(a,e16.6)') &                                         !jmecmplx
                    &'WARNING: slope along search direction is positive:', G_init%d !jmecmplx
            end if
            write(stdout,'(a)') '         Resetting conjugate gradients!'
         end if
         do isub=1,mdl%nsub
            call data_set_to_zero(direction_on_grid(isub))
            call data_functions_axpy(direction_on_grid(isub), &
                 cov_grad_on_grid(isub), -1.0_DP)
         end do
         G_init = data_functions_dot(contra_grad_on_grid, direction_on_grid)
         cg_count = 0

         ! cks: collect the work of each proc
         ! agrecocmplx
         if (loc_cmplx) then
            call comms_reduce('SUM', G_init%z)
         else
            call comms_reduce('SUM', G_init%d)
         end if

         ! take action in case of positive slope along search direction
         slope_is_positive = .false.
         ! agrecocmplx: consider sign of real part in complex case?
         if (loc_cmplx) then
            G_init_real = real(G_init%z,kind=DP)
         else
            G_init_real = G_init%d
         end if

         if (G_init_real > 0.0_DP) &
               slope_is_positive = .true.

         ! agrecocmplx
         if (slope_is_positive) then                  !jmecmplx
            if (pub_on_root .and. (pub_output_detail >= VERBOSE) ) then
               write(stdout,'(a)') 'WARNING: slope along search direction is still &
                    &positive.'
               write(stdout,'(a)') '         Reversing search direction!!'
            end if
            do isub=1,mdl%nsub
                call data_functions_scale(direction_on_grid(isub), -1.0_DP)
            end do
            call data_coef_scale(G_init, -1.0_DP)
            ! ndmh: if searching 'uphill', always re-check final step
            ! ndmh: before accepting, to avoid accepting very bad steps.
            reversing = .true.
         end if

      end if

      call services_flush

    end subroutine internal_find_direction


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    subroutine internal_print_search_info

      use rundat, only: pub_ngwf_max_grad

      implicit none

      if (pub_output_detail == NORMAL) then
         write(stdout,'(2x,a)')  repeat('=',72)
         write(stdout,'(3x,a,es26.14,a)') 'Predicted total energy    = ', &
                 predicted_functional, ' Eh'
         write(stdout,'(3x,a,es26.14,a)') 'Predicted gain in energy  = ', &
                 predicted_functional - F0, ' Eh'
         write(stdout,'(2x,a/)') repeat('=',72)
      else if (pub_output_detail >= VERBOSE) then
         write(stdout,'(a)') &
              '------------------------------- NGWF line search &
              &-------------------------------'
         write(stdout,'(a,f22.14)')   'RMS gradient                = ',rms_gradient
         if (pub_ngwf_max_grad>0.0_DP) write(stdout,'(a,f22.14)') &
              'MAX gradient                = ',max_gradient
         if (index(pub_devel_code,'NGWF_FD')>0) then
            do ifd=1,nfd
               write(stdout,'(a,f22.6)')    'FD trial step length        = ', &
                    fd_trial_step(ifd)
            end do
         end if
         write(stdout,'(a,f22.6)')    'Trial step length           = ',trial_length

         ! agrecocmplx
         if (loc_cmplx) then
            write(stdout,'(a,2(f22.11))')   'Gradient along search dir.  = ',G_init%z    !jmecmplx
         else
            write(stdout,'(a,f22.11)')   'Gradient along search dir.  = ',G_init%d
         end if
         if (index(pub_devel_code,'NGWF_FD')>0) then
            do ifd=1,nfd
               write(stdout,'(a,f22.11)')   'Gradient by FD              = ',&
                    (FFD(ifd)-F0)/fd_trial_step(ifd)
            end do
         end if
         if (check_conjugacy) then
            if (previous_dir_dot_g%iscmplx)  then
               write(stdout,'(a,2(f22.11))')   'Conjugacy Test              = ', &
                    previous_dir_dot_g%z
            else
               write(stdout,'(a,f22.11)')   'Conjugacy Test              = ', &
                    previous_dir_dot_g%d
            end if
         end if
         if (calculate_inner_grad) then
            if (inner_grad%iscmplx) then
               write(stdout,'(a,f22.11,1x,f22.11)')   'Inner gradient              = ', &
                    inner_grad%z
            else
               write(stdout,'(a,f22.11)')   'Inner gradient              = ', &
                    inner_grad%d
            end if
         end if
         write(stdout,'(a,f22.14)')   'Functional at step 0        = ',F0
         if (index(pub_devel_code,'NGWF_FD')>0) then
            do ifd=1,nfd
               write(stdout,'(a,f22.14)')   'Functional at FD step       = ',FFD(ifd)
            end do
         end if
         write(stdout,'(a,f22.14)')   'Functional at step 1        = ',F1
         if (trial2) then
            write(stdout,'(a,f22.14)')'Functional at step 2        = ',F2
            write(stdout,'(a,f22.14)')'Functional predicted        = ', &
                 predicted_functional
            write(stdout,'(a,f22.6)') 'Rejected quadratic step     = ', &
                 quadratic_coeff
            write(stdout,'(a,f22.6)') 'Selected cubic step         = ',cubic_coeff
         else if (retrial1) then
            write(stdout,'(a,f22.6)') 'Rejected quadratic step     = ', &
                 rejected_quadratic_coeff
            write(stdout,'(a,f22.14)')'Functional at new step 1    = ',F2
            write(stdout,'(a,f22.14)')'Functional predicted        = ', &
                 predicted_functional
            write(stdout,'(a,f22.6)') 'Selected quadratic step     = ', &
                 quadratic_coeff
         else
            write(stdout,'(a,f22.14)')'Functional predicted        = ', &
                 predicted_functional
            write(stdout,'(a,f22.6)') 'Selected quadratic step     = ', &
                 quadratic_coeff
         endif
         if (Ffinal/=0.0_DP) then
            write(stdout,'(a,f22.14)')'Final Functional            = ', Ffinal
         end if
         if ( pub_elec_cg_max > 0 ) then
            ! agrecocmplx
            if (loc_cmplx) then
               write(stdout,'(a,2(f22.6))')    'Conjugate gradients coeff.  = ',cg_coeff%z
            else
               write(stdout,'(a,f22.6)')    'Conjugate gradients coeff.  = ',cg_coeff%d  !jmecmplx
            end if
         endif
         write(stdout,'(a)')'--------------------------- NGWF line search &
              &finished --------------------------'
         write(stdout,'(a)')'                                               '
      endif

    end subroutine internal_print_search_info

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    type(COEF) function internal_inner_grad(contra_on_grid, overlap, fun_basis)

      ! Calculates the inner product of the gradient with itself, using the
      ! overlap as the metric.
      ! Written by Alvaro Ruiz Serrano in January 2013.

      use datatypes, only: COEF, FUNCTIONS, data_set_to_zero, &
           data_functions_alloc, data_functions_dealloc
      use function_ops, only: function_ops_sum_ppd_funcs

      implicit none

      ! ars: arguments
      type(FUNC_BASIS), intent(in) :: fun_basis
      type(FUNCTIONS), intent(in) :: contra_on_grid
      type(SPAM3), intent(in) :: overlap

      ! ars: local variables
      type(FUNCTIONS) :: cov_on_grid(1)
      type(SPAM3) :: overlap_matrix(1)

      ! ars: allocate
      call data_functions_alloc(cov_on_grid(1), fun_basis%size_on_grid, &
           iscmplx=contra_on_grid%iscmplx)
      call sparse_create(overlap_matrix(1), overlap)

      ! ars: init
      call data_set_to_zero(cov_on_grid(1))
      call sparse_copy(overlap_matrix(1), overlap)

      ! ars: calculate cov_on_grid = contra_on_grid * S_\alpha\beta
      call function_ops_sum_ppd_funcs(cov_on_grid(1:1), &        ! inout
           fun_basis, overlap_matrix(1:1), 1, 1, overlap, &           ! input
           contra_on_grid, fun_basis)

      ! ars: calculate g_dot_g
      internal_inner_grad = data_functions_dot(contra_on_grid, cov_on_grid(1))

      ! ars: deallocate
      call sparse_destroy(overlap_matrix(1))
      call data_functions_dealloc(cov_on_grid(1))

    end function internal_inner_grad


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_rotate_prev_search_dir(old_dir_on_grid, fun_basis, &
         old_fun_on_grid, new_fun_on_grid, old_inv_overlap)

      ! Rotates the previous NGWF search direction to the new NGWF representation
      ! Written by Alvaro Ruiz Serrano in January 2013.

      use basis, only: basis_clean_function, basis_extract_function_from_box
      use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_functions_copy, &
           data_functions_alloc, data_functions_dealloc, &
           data_fftbox_alloc, data_fftbox_dealloc, data_set_to_zero
      use function_ops, only: function_ops_sum_fftbox_batch, &
           function_ops_batch_col_start, function_ops_brappd_ketppd
      use rundat, only: pub_fftbox_batch_size
      use sparse, only: sparse_product, sparse_index_length, sparse_generate_index

      implicit none

      ! ars: arguments
      type(FUNC_BASIS), intent(in) :: fun_basis
      type(FUNCTIONS), intent(inout) :: old_dir_on_grid
      type(FUNCTIONS), intent(in) :: old_fun_on_grid
      type(FUNCTIONS), intent(in) :: new_fun_on_grid
      type(SPAM3), intent(in) :: old_inv_overlap

      ! ars: local variables
      type(FUNCTIONS) :: buffer_on_grid
      type(SPAM3) :: overlap_new_old, os(1)
      type(FFTBOX_DATA), allocatable :: fftbox_batch(:,:,:)
      integer :: local_fa, local_start, local_end, local_len
      integer :: batch_size, max_current_size
      integer :: idx_len, batch_count, n_batches, batch_index
      integer :: is, ib
      integer, allocatable, dimension(:) :: overlap_idx
      integer, allocatable :: fa_box_start(:,:)
      integer, allocatable :: fa_start_in_box(:,:)
      ! agrecocmplx
      logical :: loc_cmplx


      if (pub_on_root) write(stdout,*) "Entering internal_rotate_prev_search_dir"

      ! agrecocmplx
      loc_cmplx = old_fun_on_grid%iscmplx

      ! ars: set batch_sizes and n_batches
      batch_size = pub_fftbox_batch_size
      n_batches = ngwf_basis(ireg)%max_on_proc / batch_size
      if (mod(ngwf_basis(ireg)%max_on_proc,batch_size) > 0) n_batches = n_batches + 1

      ! ars: allocate buffer
      call data_functions_alloc(buffer_on_grid, &
           ngwf_basis(ireg)%size_on_grid, iscmplx=loc_cmplx)

      ! ars: create sparse structures
      overlap_new_old%structure = 'S'
      ! agrecocmplx
      call sparse_create(overlap_new_old, iscmplx=loc_cmplx)
      call sparse_create(os(1), overlap_new_old, old_inv_overlap)

      ! ars: calculate new-old overlap
      call function_ops_brappd_ketppd(overlap_new_old, &
           new_fun_on_grid,fun_basis,old_fun_on_grid,fun_basis,mdl%cell)

      ! ars: calculate OS product
      call sparse_product(os(1), overlap_new_old, old_inv_overlap)

      ! pdh: overlap matrix index
      idx_len = sparse_index_length(rep%ngwf_overlap%m(ireg,ireg))
      allocate(overlap_idx(idx_len),stat=ierr)
      call utils_alloc_check('internal_rotate_prev_search_dir','overlap_idx',ierr)

      ! ars: calculate \sum_\beta \phi_\beta (r) * <d|phi>S^
      call sparse_generate_index(overlap_idx,rep%ngwf_overlap%m(ireg,ireg))

      ! ndmh: fb start positions in box and start positions of FFT boxes
      allocate(fa_box_start(3, batch_size), stat=ierr)
      call utils_alloc_check('internal_rotate_prev_search_dir','fa_box_start',ierr)
      allocate(fa_start_in_box(3, batch_size), stat=ierr)
      call utils_alloc_check('internal_rotate_prev_search_dir','fa_start_in_box',ierr)

      ! ars: allocate storage for fftboxes for this batch
      allocate(fftbox_batch(pub_num_spins,1,batch_size), stat=ierr)
      call utils_alloc_check('internal_rotate_prev_search_dir','fftbox_batch',ierr)
      do batch_count = 1, batch_size
         do is = 1, pub_num_spins
            ! agrecocmplx
            call data_fftbox_alloc(fftbox_batch(is,1,batch_count), &
                 mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, &
                 mdl%fftbox%total_pt3, iscmplx=loc_cmplx)
         end do
      end do

      ! ars: loop over batches of NGWFs
      local_start = 1
      do batch_count=1,n_batches

         local_end = min(local_start+batch_size-1,ngwf_basis(ireg)%proc_num)
         local_len = local_end - local_start + 1

         ! cks: maximum size of current batch over all procs
         max_current_size = local_len
         call comms_reduce('MAX', max_current_size)

         call function_ops_batch_col_start(fa_box_start,fa_start_in_box, &
              batch_size,local_start,local_end,mdl%fftbox,mdl%cell, &
              ngwf_basis(ireg))

         ! cks: zero before accumulation
         do is = 1, pub_num_spins
            do ib = 1, local_len
               call data_set_to_zero(fftbox_batch(is,1,ib))
            end do
         end do

         call function_ops_sum_fftbox_batch(fftbox_batch, mdl%fftbox, mdl%cell, &
              old_dir_on_grid, ngwf_basis(ireg), fa_box_start, &
              batch_size, local_start, local_end, overlap_idx, idx_len, os(1), &
              1, 1, 1.0_DP)

         ! ars: shave
         batch_index = 0
         do local_fa = local_start, local_end

            batch_index = batch_index + 1

            if (pub_on_root) write(stdout,*) "Point O"

            ! ars: FFTbox to PPD
            call basis_extract_function_from_box(buffer_on_grid, &
                 fftbox_batch(1,1,batch_index), &
                 ngwf_basis(ireg)%spheres(local_fa),&
                 ngwf_basis(ireg)%tight_boxes(local_fa), &
                 fa_start_in_box(1,batch_index), &
                 fa_start_in_box(2,batch_index), &
                 fa_start_in_box(3,batch_index), &
                 ngwf_basis(ireg)%spheres(local_fa)%offset,mdl%cell,mdl%fftbox)

            ! ars: PPD to sphere
            call basis_clean_function(buffer_on_grid, &
                 ngwf_basis(ireg)%spheres(local_fa), mdl%cell,mdl%fftbox)

         end do

         local_start = local_start + batch_size

      end do


      do batch_count = batch_size, 1, -1
         do is = pub_num_spins, 1, -1
            call data_fftbox_dealloc(fftbox_batch(is,1,batch_count))
         end do
      end do
      deallocate(fftbox_batch, stat=ierr)
      call utils_dealloc_check('internal_rotate_prev_search_dir','fftbox_batch',ierr)
      deallocate(fa_start_in_box,stat=ierr)
      call utils_dealloc_check('internal_rotate_prev_search_dir','fa_start_in_box',ierr)
      deallocate(fa_box_start,stat=ierr)
      call utils_dealloc_check('internal_rotate_prev_search_dir','fa_box_start',ierr)
      deallocate(overlap_idx,stat=ierr)
      call utils_dealloc_check('internal_rotate_prev_search_dir','overlap_idx',ierr)

      ! ars: destroy sparse_structures
      call sparse_destroy(os(1))
      call sparse_destroy(overlap_new_old)

      ! ars: copy to old_dir_on_grid
      call data_functions_copy(old_dir_on_grid, buffer_on_grid)

      ! ars: deallocate buffer
      call data_functions_dealloc(buffer_on_grid)

      if (pub_on_root) write(stdout,*) "internal_rotate_prev_search_dir"


    end subroutine internal_rotate_prev_search_dir

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_orthogonalise_direction(array_on_grid)

      use basis, only: basis_clean_function, basis_extract_function_from_box
      use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero, &
           data_functions_alloc, data_functions_dealloc, &
           data_fftbox_alloc, data_fftbox_dealloc, &
           data_functions_axpy
      use function_ops, only: function_ops_sum_fftbox_batch, &
           function_ops_batch_col_start, function_ops_brappd_ketppd
      use rundat, only: pub_fftbox_batch_size
      use sparse, only: sparse_product, sparse_index_length, sparse_generate_index

      implicit none

      ! ars: arguments
      type(FUNCTIONS), intent(inout) :: array_on_grid

      ! ars: local variables
      type(FUNCTIONS) :: projection_on_grid
      type(FFTBOX_DATA), allocatable :: fftbox_batch(:,:,:)
      type(SPAM3) :: overlap_dir_ngwfs, cc(1)
      integer :: local_fa, local_start, local_end, local_len
      integer :: batch_size, max_current_size
      integer :: is, ib
      integer :: idx_len, batch_count, n_batches, batch_index
      integer, allocatable, dimension(:) :: overlap_idx
      integer, allocatable :: fa_box_start(:,:)
      integer, allocatable :: fa_start_in_box(:,:)
      ! agrecocmplx
      logical :: loc_cmplx

      if (pub_on_root) write(stdout,*) "Entering internal_orthogonalise_direction"

      ! agrecocmplx
      loc_cmplx = rep%ngwfs_on_grid(ireg)%iscmplx

      ! ars: set batch_sizes and n_batches
      batch_size = pub_fftbox_batch_size
      n_batches = ngwf_basis(ireg)%max_on_proc / batch_size
      if (mod(ngwf_basis(ireg)%max_on_proc,batch_size) > 0) n_batches = n_batches + 1

      ! ars: allocate buffer
      call data_functions_alloc(projection_on_grid, ngwf_basis(ireg)%size_on_grid, &
           iscmplx=loc_cmplx)

      ! ars: create sparse structures
      call sparse_create(overlap_dir_ngwfs, rep%overlap%m(ireg,ireg))
      call sparse_create(cc(1), overlap_dir_ngwfs, rep%inv_overlap%m(ireg,ireg))

      ! ars: calculate overlap between array_on_grid and the NGWFs
      call function_ops_brappd_ketppd(overlap_dir_ngwfs, &
           array_on_grid, ngwf_basis(ireg), &
           rep%ngwfs_on_grid(mdl%regions(ireg)%reg_count), &
           ngwf_basis(mdl%regions(ireg)%reg_count), mdl%cell)

      ! ars: calculate product C=<d|phi>S^
      call sparse_product(cc(1), overlap_dir_ngwfs, rep%inv_overlap%m(ireg,ireg))

      ! pdh: overlap matrix index
      idx_len = sparse_index_length(rep%ngwf_overlap%m(ireg,ireg))
      allocate(overlap_idx(idx_len),stat=ierr)
      call utils_alloc_check('internal_orthogonalise_direction','overlap_idx',ierr)

      ! ars: calculate \sum_\beta \phi_\beta (r) * <d|phi>S^
      call sparse_generate_index(overlap_idx,rep%ngwf_overlap%m(ireg,ireg))

      ! ndmh: fb start positions in box and start positions of FFT boxes
      allocate(fa_box_start(3, batch_size), stat=ierr)
      call utils_alloc_check('internal_orthogonalise_direction','fa_box_start',ierr)
      allocate(fa_start_in_box(3, batch_size), stat=ierr)
      call utils_alloc_check('internal_orthogonalise_direction','fa_start_in_box',ierr)

      ! ars: allocate storage for fftboxes for this batch
      allocate(fftbox_batch(pub_num_spins,1,batch_size), stat=ierr)
      call utils_alloc_check('internal_orthogonalise_direction','fftbox_batch',ierr)
      do batch_count = 1, batch_size
         do is = 1, pub_num_spins
            call data_fftbox_alloc(fftbox_batch(is,1,batch_count), &
                 mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, &
                 mdl%fftbox%total_pt3, iscmplx=loc_cmplx)
         end do
      end do

      ! ars: loop over batches of NGWFs
      local_start = 1
      do batch_count=1,n_batches

         local_end = min(local_start+batch_size-1,ngwf_basis(ireg)%proc_num)
         local_len = local_end - local_start + 1

         ! cks: maximum size of current batch over all procs
         max_current_size = local_len
         call comms_reduce('MAX', max_current_size)

         call function_ops_batch_col_start(fa_box_start,fa_start_in_box, &
              batch_size,local_start,local_end,mdl%fftbox,mdl%cell, &
              ngwf_basis(ireg))

         ! cks: zero before accumulation
         do is = 1, pub_num_spins
            do ib = 1, local_len
               call data_set_to_zero(fftbox_batch(is,1,ib))
            end do
         end do

         call function_ops_sum_fftbox_batch(fftbox_batch, mdl%fftbox, mdl%cell, &
              rep%ngwfs_on_grid(ireg), ngwf_basis(ireg), fa_box_start, &
              batch_size, local_start, local_end, overlap_idx, idx_len, cc(1), &
              1, 1, 1.0_DP)

         ! ars: shave
         batch_index = 0
         do local_fa = local_start, local_end

            batch_index = batch_index + 1

            ! ars: FFTbox to PPD
            call basis_extract_function_from_box(projection_on_grid, &
                 fftbox_batch(1,1,batch_index), &
                 ngwf_basis(ireg)%spheres(local_fa),&
                 ngwf_basis(ireg)%tight_boxes(local_fa), &
                 fa_start_in_box(1,batch_index), &
                 fa_start_in_box(2,batch_index), &
                 fa_start_in_box(3,batch_index), &
                 ngwf_basis(ireg)%spheres(local_fa)%offset,mdl%cell,mdl%fftbox)

            ! ars: PPD to sphere
            call basis_clean_function(array_on_grid, &
                 ngwf_basis(ireg)%spheres(local_fa), mdl%cell,mdl%fftbox)

         end do

         local_start = local_start + batch_size
      end do

      do batch_count = batch_size, 1, -1
         do is = pub_num_spins, 1, -1
            call data_fftbox_dealloc(fftbox_batch(is,1,batch_count))
         end do
      end do
      deallocate(fftbox_batch, stat=ierr)
      call utils_dealloc_check('internal_orthogonalise_direction','fftbox_batch',ierr)
      deallocate(fa_start_in_box,stat=ierr)
      call utils_dealloc_check('internal_orthogonalise_direction','fa_start_in_box',ierr)
      deallocate(fa_box_start,stat=ierr)
      call utils_dealloc_check('internal_orthogonalise_direction','fa_box_start',ierr)
      deallocate(overlap_idx,stat=ierr)
      call utils_dealloc_check('internal_orthogonalise_direction','overlap_idx',ierr)

      ! ars: destroy sparse_structures
      call sparse_destroy(cc(1))
      call sparse_destroy(overlap_dir_ngwfs)

      ! ars: orthogonalise direction
      call data_functions_axpy(array_on_grid, projection_on_grid, -1.0_DP)

      ! ars: deallocate buffer
      call data_functions_dealloc(projection_on_grid)

      if (pub_on_root) write(stdout,*) "Leaving internal_orthogonalise_direction"

    end subroutine internal_orthogonalise_direction

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_write_soc_ham()

      use augmentation, only: aug_projector_denskern
      use hamiltonian, only: hamiltonian_soc_matrices
      use rundat, only: pub_num_spins, pub_cond_calculate, pub_paw
      use sparse_embed, only: sparse_embed_extract_from_array, &
           sparse_embed_destroy_extracted_array
      use utils, only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Arguments

      ! Local Variables
      integer :: is
      type(SPAM3_EMBED), allocatable :: rhoij(:), hamso(:)
      type(SPAM3), allocatable :: aug_array(:), kern_array(:)

      ! ndmh: return if not PAW
      if (.not.pub_paw) then
         write(stdout,*) 'WARNING: Perturbative SOC hamiltonian &
              &is only available in PAW.'
         return
      end if

      ! Set up matrices
      allocate(hamso(4),stat=ierr)
      call utils_alloc_check('internal_write_soc_ham','hamso',ierr)
      allocate(rhoij(pub_num_spins),stat=ierr)
      call utils_alloc_check('internal_write_soc_ham','rhoij',ierr)
      allocate(aug_array(pub_num_spins),stat=ierr)
      call utils_alloc_check('internal_write_soc_ham','aug_array',ierr)
      allocate(kern_array(pub_num_spins),stat=ierr)
      call utils_alloc_check('internal_write_soc_ham','kern_array',ierr)

      ! Build rhoij
      do is=1,pub_num_spins
         rhoij(is)%structure = 'E'
         call sparse_embed_create(rhoij(is))
      end do
      call sparse_embed_extract_from_array(aug_array,rhoij)
      if (pub_cond_calculate) then
         call sparse_embed_extract_from_array(kern_array,val_dkn%m(:,PUB_1K))
         call aug_projector_denskern(aug_array,kern_array,val_rep%sp_overlap%p)
      else
         call sparse_embed_extract_from_array(kern_array,denskern%kern%m(:,PUB_1K))
         call aug_projector_denskern(aug_array,kern_array,rep%sp_overlap%p)
      end if
      call sparse_embed_destroy_extracted_array(aug_array,rhoij,.true.)
      call sparse_embed_destroy_extracted_array(kern_array)

      ! Allocate 4 components of Hamiltonian
      do is=1,4
         call sparse_embed_create(hamso(is),ham%ham(1),iscmplx=.true.)
      end do

      ! Calculate S.O. terms
      call hamiltonian_soc_matrices(hamso(:),ham,rhoij,rep,mdl)
      ! Write to binary file
      call restart_hamiltonian_write(hamso,write_soc=.true.)
      ! Clean up matrices
      do is=4,1,-1
         call sparse_embed_destroy(hamso(is))
      end do
      do is=pub_num_spins,1,-1
         call sparse_embed_destroy(rhoij(is))
      end do
      deallocate(rhoij,stat=ierr)
      call utils_dealloc_check('internal_write_soc_ham','rhoij',ierr)
      deallocate(aug_array,stat=ierr)
      call utils_dealloc_check('internal_write_soc_ham','aug_array',ierr)
      deallocate(kern_array,stat=ierr)
      call utils_dealloc_check('internal_write_soc_ham','kern_array',ierr)
      deallocate(hamso,stat=ierr)
      call utils_dealloc_check('internal_write_soc_ham','hamso',ierr)

    end subroutine internal_write_soc_ham

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! rc2013: do various checks for convergence with embedding
    logical function internal_ft_convergence()

      implicit none

      if(pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: call internal_ft_convergence'
      if(converged .and. pub_do_fandt) then
         ft_converged(ireg) = .true.
         ift = 1
         ! rc2013: this region is converged, check if the others are too
         conv_loop: do isub=1,mdl%nsub
            if(ft_converged(isub)) then
               internal_ft_convergence = .true.
            else
               internal_ft_convergence = .false.
               exit conv_loop
            endif
         enddo conv_loop
         if(pub_freeze_envir_ngwfs) internal_ft_convergence = .true.
      else if(pub_do_fandt) then
         ! rc2013: only declare convergence if all regions converge in succession
         ft_converged = .false.
         internal_ft_convergence = .false.
      else if(pub_emft_follow .and. .not. fnl_switched) then
         ! rc2013: if we're doing an EMFT follow-up and haven't swapped
         ! functionals yet then we're not really converged
         internal_ft_convergence = .false.
      else
         ! rc2013: default case
         internal_ft_convergence = converged
      endif
      if(pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: leaving internal_ft_convergence'

    end function internal_ft_convergence

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_region_reset()

      use rundat, only: pub_build_bo, pub_block_orthogonalise

      implicit none

      if(pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering internal_region_reset'

      ! rc2013: reset the F+T counter
      ift = 1
      ! rc2013: if we've gone through all regions, start over
      if(pub_do_fandt .and. .not. pub_freeze_envir_ngwfs) then
         if(ireg == mdl%nsub) then
            ireg = 1
         else
            ireg = ireg + 1
         endif
      end if

      ! rc2013: reset search parameters to default values (new SD search)
      line_search_coeff   = 0.15_DP
      trial_length        = 0.1_DP
      line_search_success = .true.
      cg_count            = 0

      ! rc2013: if we're doing EMFT after the NGWF optimisation, set it up now
      if(pub_emft_follow) then
         orig_pub_emft = pub_emft
         pub_emft = .true.
         if (pub_block_orthogonalise) then
             ! rc2013: recalculate the overlap matrices etc in BO form
             orig_pub_build_bo = pub_build_bo
             pub_build_bo = .true.
             call hamiltonian_dens_indep_matrices(rep, ngwf_basis, proj_basis, &
                  nl_projectors, hub_proj_basis, hub, mdl, val_rep, val_ngwf_basis, &
                  kpt=loc_kpt_cart)
         end if
         current_maxit_lnv = pub_emft_lnv_steps ! run extra-long DK optimisation
         ! rc2013: we must reset the kernel and the Hamiltonian before the next iteration
         total_energy = electronic_energy(denskern, pur_denskern, ham, &
              lhxc_fine, mu, edft, rep, ngwf_basis, hub_proj_basis, hub, &
              mdl, hfxstate, lnv_threshold, current_maxit_lnv, &
              kernel_update=.true., conv_status=cvstat, &
              dfdtau_fine = dfdtau_fine )

         ! rc2013: if it's an LNV-only run then reset maxit_ngwf_cg
         if (pub_emft_lnv_only) then
            orig_pub_maxit_ngwf_cg = pub_maxit_ngwf_cg
            pub_maxit_ngwf_cg = 0
         end if
         fnl_switched = .true.
      end if

      if(pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving internal_region_reset'

    end subroutine internal_region_reset

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_allocate

        implicit none

        ! Allocate workspace
        cg_coeff%iscmplx = loc_cmplx   !jmecmplx
        G_init%iscmplx = loc_cmplx
        inner_grad%iscmplx = loc_cmplx
        previous_g_dot_g%iscmplx = loc_cmplx
        current_g_dot_g%iscmplx = loc_cmplx
        !ep: mermin var
        pre_dir_dot_g%iscmplx = loc_cmplx
        now_dir_dot_g%iscmplx = loc_cmplx


        ! rc2013: create structures for all subsystems
        do isub=1,mdl%nsub
           call data_functions_alloc(contra_grad_on_grid(isub), &
                ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
           call data_functions_alloc(direction_on_grid(isub), &
                ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
           call data_functions_alloc(start_ngwfs_on_grid(isub), &
                ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
           !ep: alloc mermin vars
           if (pub_elec_cg_max > 0) then
              ! cks: allocate properly only when not doing steepest descents
              call data_functions_alloc(prev_direction_on_grid(isub), &
                   ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
           else
              ! cks: otherwise allocate only token memory to keep compiler happy
              call data_functions_alloc(prev_direction_on_grid(isub), 1, &
                   iscmplx=loc_cmplx)
           end if

           call data_functions_alloc(cov_grad_on_grid(isub), &
                ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)

           if (pub_ngwf_cg_type == 'NGWF_POLAK') then
              call data_functions_alloc(prev_contra_grad_on_grid(isub), &
                   ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
           end if

           !ep: alloc mermin vars
           if (pub_mermin) then
              call data_functions_alloc(prev_contra_grad_on_grid(isub), &
                   ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
           end if


           ! ndmh: for force convergence testing
           if (pub_elec_force_tol > 0.0_DP) then
              allocate(total_forces(3,mdl%nat,3),stat=ierr)
              call utils_alloc_check('ngwf_cg_optimise','total_forces',ierr)
           end if
        end do

    end subroutine internal_allocate

  end subroutine ngwf_cg_optimise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ngwf_cg_emft_reset

    !======================================================================!
    ! This subroutine resets NGWF optimisation global variables to         !
    ! their original values if we have been doing an EMFT follow           !
    ! calculation, in case another set of NGWF optimisations is            !
    ! required (e.g. conduction)                                           !
    !======================================================================!
    ! Written by Joseph Prentice, May 2019                                 !
    !======================================================================!

    use rundat, only: pub_emft, pub_build_bo, pub_maxit_ngwf_cg, &
         pub_block_orthogonalise, pub_emft_lnv_only

    implicit none

    pub_emft = orig_pub_emft
    if (pub_block_orthogonalise) pub_build_bo = orig_pub_build_bo
    if (pub_emft_lnv_only) pub_maxit_ngwf_cg = orig_pub_maxit_ngwf_cg

  end subroutine ngwf_cg_emft_reset


end module ngwf_cg
