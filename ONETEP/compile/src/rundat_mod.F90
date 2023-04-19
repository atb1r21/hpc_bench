! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!======================================================================!
!                                                                      !
!            The ONETEP code is written and maintained by              !
!     Chris-Kriton Skylaris, Arash A. Mostofi, Nicholas D.M. Hine,     !
!                 Jacek Dziedzic and Peter D. Haynes.                  !
!                                                                      !
!======================================================================!
module rundat

  use constants, only: DP, PI, HARTREE_IN_EVS, METRIC_ELECTROSTATIC, &
       METRIC_OVERLAP, SAFE_DIV_EPS, mu_ref_SHE

  implicit none

  ! ndmh: General properties
  character(len=80)            :: pub_rootname  ! set in onetep.F90, modified in e&f
  character(len=80)            :: pub_task          ! set in onetep.F90
  character(len=80)            :: pub_all_tasks ! modified in onetep.F90
  real(kind=DP), protected     :: pub_cutoff_energy
  real(kind=DP), protected     :: pub_kernel_cutoff
  character(len=80)            :: pub_xc_functional
  integer, protected           :: pub_libxc_x_func_id
  integer, protected           :: pub_libxc_c_func_id
  real(kind=DP)                :: pub_charge  ! modified in eda_main_mod
  integer                      :: pub_spin   ! modified in energy_and_force_mod
  real(kind=DP)                :: pub_real_spin ! modified in ensemble_dft_mod and energy_and_force_mod
  logical                      :: pub_spin_polarised ! modified in energy_and_force_mod
  logical                      :: pub_confined_ngwfs ! gcc32
  integer                      :: pub_maxit_ngwf_cg_confined ! gcc32
  real(kind=DP)                :: pub_confined_ngwfs_barrier
  integer                      :: pub_num_spins   ! expected: set in energy_and_force_mod
  real(kind=DP)                :: pub_spin_fac   ! expected: set in energy_and_force_mod
  integer, parameter           :: pub_num_kpoints=1
  integer, parameter           :: PUB_1K=1
  real(kind=DP), protected     :: pub_constant_efield(3)
  integer, protected           :: pub_fftbox_pref(3)
  integer, protected           :: pub_augbox_pref(3)
  integer, protected           :: pub_ppd_npoints(3)
  real(kind=DP)                :: pub_psinc_spacing(3) ! modified in rundat_blocks_mod
  character(len=80), protected :: pub_devel_code
  logical, protected           :: pub_nnho
  real(kind=DP), protected     :: pub_ngwf_halo
  logical, protected           :: pub_nonsc_forces
  real(kind=DP), protected     :: pub_external_pressure
  real(kind=DP), protected     :: pub_smoothing_factor
  real(kind=DP), protected     :: pub_isosurface_cutoff
  character(len=80), protected :: pub_input_xyz_file
  character(len=80), protected :: pub_input_xyz_file_intermediate
  character(len=80), protected :: pub_input_xyz_file_product
  ! agreco: directions along which NGWFs are extended
  logical, protected           :: pub_extend_ngwf(3)
  ! agreco: initialise NGWFs to completely random values
  logical, protected           :: pub_full_rand_ngwf
  ! agreco: user specified seed for pseudo-random generator
  integer, protected           :: pub_rand_seed_ngwf
  ! agreco: sigma of normal distribution of random numbers
  real(kind=DP), protected     :: pub_rand_sigma
  ! agreco: convolute fully random NGWFs to smooth oscillations
  ! at the edges of localisation region
  logical, protected           :: pub_rand_conv
  ! agreco: convolution function
  character(len=80)            :: pub_conv_func
  ! agreco: write NGWFs for plotting at every step
  logical, protected           :: pub_ngwf_plot_every_it
  ! agreco: use complex valued NGWFs
  logical                      :: pub_cmplx_ngwfs
  ! agreco: read real NGWFs from file into complex NGWFs
  logical                      :: pub_read_real_ngwfs
  ! agreco: origin of the sawtooth potential associated to e-field
  real(kind=DP), protected     :: pub_efield_origin(3)
  ! agreco: store overlap matrix to file at the end of calculation
  logical, protected           :: pub_show_overlap
  ! agreco: specify width of region to which apply convolution
  real(kind=DP), protected     :: pub_conv_region_width
  character(len=40)            :: conv_region_width_unit
  ! agreco: specify single k-point for testing
  ! k-point is specified in fractional coordinates
  real(kind=DP), protected     :: pub_single_kpt(3)
  ! agreco: method to be used for k-point sampling of BZ
  character(len=40), protected :: pub_kpoint_method
  ! agrecocmplx: threshold to check imaginary part of supposedly
  ! real quantities
  real(kind=DP), protected     :: pub_imag_thr
  ! agrecocmplx: phase to add to complex NGWFs (testing)
  real(kind=DP), protected     :: pub_ngwfs_phase
  ! agrecocmplx: if applying a random phase to complex NGWFs (testing)
  logical, protected           :: pub_ngwfs_rand_phase
  ! agrecocmplx: if checking hermitian character of complex relevant
  ! matrices S/H/K and so on
  logical, protected           :: pub_check_hermitian
  ! agrecocmplx: check density is real when using complex NGWFs
  logical, protected           :: pub_check_density
  ! agrecokpt: size of k-point grid
  integer, protected           :: pub_kp_grid_size(3)
  ! agrecokpt: shif of k-point grid
  integer, protected           :: pub_kp_grid_shift(3)

  ! ja531-> contraco Ham cutoff multiplier
  real(kind=dp), protected     :: pub_contracoham_radmult

  ! jd: Debug?
#ifdef  DEBUG
  !jd:^^ Please keep the double space in the line above
  logical, parameter  :: pub_debug = .true.
#else
  logical, parameter  :: pub_debug = .false.
#endif
  logical             :: pub_debug_on_root
  logical, protected  :: pub_check_stack_size

  ! ndmh: Kernel Optimisation parameters
  integer, protected           :: pub_maxit_palser_mano
  integer, protected           :: pub_maxit_kernel_occ_check
  integer                      :: pub_maxit_pen ! modified in energy_and_force_mod
  real(kind=DP), protected     :: pub_pen_param
  integer                      :: maxit_lnv ! modified in md_mod
  integer                      :: minit_lnv ! modified in md_mod
  integer                      :: pub_firstit_lnv
  real(kind=DP)                :: lnv_threshold_orig ! modified in md_mod
  character(len=80), protected :: pub_lnv_cg_type
  real(kind=DP), protected     :: pub_lnv_cg_max_step
  logical                      :: pub_exact_lnv ! modified in eda_driver_supermol_mod
  logical, protected           :: pub_old_lnv
  logical, protected           :: pub_lnv_check_trial_steps
  integer                      :: pub_kerfix ! modified in eda_driver_supermol_mod
  integer, protected           :: pub_maxit_kernel_fix
  logical, protected           :: pub_initial_dens_realspace
  logical, protected           :: pub_kernel_track_mid_occ
  logical, protected           :: pub_kernel_check_all

  logical, protected           :: pub_kernel_diis          !use kernel DIIS (not a keyword)
  integer                      :: pub_kernel_diis_size      !max # of kernels to store in memory   ! modified in kernel_diis_mod
  integer, protected           :: pub_kernel_diis_maxit    !max # of DIIS iterations
  real(kind=DP)                :: pub_kernel_diis_threshold    !DIIS convergence threshold             ! modified in kernel_diis_mod
  integer                      :: pub_kernel_diis_linear_iter    !#linear DIIS iter before Pulay DIIS  ! modified in kernel_diis_mod
  real(kind=dp)                :: pub_kernel_diis_coeff    !coefficients for linear DIIS               ! modified in kernel_diis_mod
  character(len=4), protected  :: pub_kernel_diis_conv_criteria !convergence criteria
  real(kind=DP), protected     :: pub_kernel_diis_lshift   ! cks: value of the level shifter parameter
  real(kind=DP), protected     :: pub_kernel_diis_ls_iter  !number of level shifting iteration
  character(len=80), protected :: pub_kernel_diis_scheme   ! ars: scheme of dkn optimisation in inner loop

  logical, protected           :: pub_foe
  logical, protected           :: pub_dense_foe
  logical, protected           :: pub_H2denskern_sparsity
  real(kind=dp), protected     :: pub_foe_mu_tol
  logical, protected           :: pub_foe_test_sparsity
  real(kind=dp), protected     :: pub_foe_cheby_thres
  logical, protected           :: pub_foe_avoid_inversions
  logical, protected           :: pub_foe_check_entropy

  logical, protected           :: pub_edft
  integer                      :: pub_edft_maxit ! modified in ensemble_dft_mod
  real(kind=DP)                :: pub_edft_max_step ! modified in ensemble_dft_mod
  real(kind=DP), protected     :: pub_edft_smearing_width
  real(kind=DP), protected     :: pub_edft_free_energy_thres
  real(kind=DP), protected     :: pub_edft_energy_thres
  real(kind=DP), protected     :: pub_edft_entropy_thres
  real(kind=DP), protected     :: pub_edft_rms_gradient_thres
  real(kind=DP), protected     :: pub_edft_commutator_thres
  real(kind=DP), protected     :: pub_edft_fermi_thres
  real(kind=DP), protected     :: pub_edft_trial_step
  real(kind=DP), protected     :: pub_edft_nelec_thres
  logical, protected           :: pub_edft_write_occ
  integer                      :: pub_edft_extra_bands ! modified in ensemble_dft_mod
  integer                      :: pub_edft_spin_fix ! kkbd - Number of NGWF CG iterations to hold the spin fixed.
                                                    !        If negative, hold forever. (Default: -1)
  integer                      :: pub_edft_spin_fix_orig
  character(len=80), protected :: pub_edft_update_scheme   ! ars: scheme of dkn optimisation in inner loop
  integer, protected           :: pub_edft_ham_diis_size        ! max # of ham to store in memory
  real(kind=DP), protected     :: pub_eigensolver_orfac
  real(kind=DP), protected     :: pub_eigensolver_abstol
  integer, protected           :: pub_edft_round_evals
  integer                      :: pub_edft_init_maxit
  logical, protected           :: pub_edft_grand_canonical   ! ab: Grand Canonical ensemble DFT
  real(kind=DP), protected     :: pub_edft_reference_potential
  real(kind=DP), protected     :: pub_edft_electrode_potential
  logical                      :: pub_pbc_smeared_ion_rep ! ja531-> whether we are doing PBC smeared ions in
                                                          ! properties / electrostatic potential output
  logical                      :: pub_chemical_softness

  integer                      :: pub_inner_loop_iteration = 0
  !                   modified in lnv, kernel_diis, ensemble_dft, ngwf_cg
  ! jd:                           ^ Globally accessible copy of inner loop
  !                                 iteration number, only used for controlling
  !                                 devel_code dumps of grid quantities from
  !                                 various modules. Works with LNV, KDIIS, EDFT.
  !                                 Iterations are numbered from 1.
  !                                 0 indicates we're in properties
  !                                -1 indicates we're in the NGWF gradient

  ! ndmh: NGWF optimisation parameters
  logical                      :: pub_ngwf_gradient_needed = .false. ! modified in e&f
  integer                      :: pub_maxit_ngwf_cg ! modified in hubbard_build_mod
  character(len=80), protected :: pub_ngwf_cg_type
  real(kind=DP), protected     :: pub_ngwf_cg_max_step
  integer, protected           :: pub_elec_cg_max
  character(len=80), protected :: pub_precond_scheme
  real(kind=DP), protected     :: pub_k_zero
  real(kind=DP), protected     :: pub_r_precond
  logical, protected           :: pub_precond_recip
  logical, protected           :: pub_precond_real
  character(len=80), protected :: pub_smooth_scheme
  real(kind=DP), protected     :: pub_r_smooth
  real(kind=DP), protected     :: pub_k_smooth
  real(kind=DP), protected     :: pub_occ_mix
  logical, protected           :: pub_kernel_update
  logical, protected           :: pub_kernel_christoffel_update
  integer, protected           :: pub_maxit_hotelling
  real(kind=DP), protected     :: pub_max_resid_hotelling
  logical                      :: pub_use_aux_ngwfs ! modified in rundat_blocks_mod
  integer                      :: pub_energy_components_interval

  ! ndmh: Convergence criteria
  real(kind=DP)                :: ngwf_threshold_orig  ! modified in md_mod
  real(kind=DP)                :: pub_elec_energy_tol   ! Max energy/atom change per iteration ! modified in md_mod
  real(kind=DP)                :: pub_elec_force_tol    ! Max force change per iteration      ! modified in md_mod
  real(kind=DP)                :: pub_ngwf_max_grad     ! Maximum permissable NGWF gradient   ! modified in md_mod
  logical, protected           :: pub_kernel_force_conv  ! Force denskern convergence on last NGWF iteration
  logical, protected           :: pub_delta_e_conv

  ! rc2013: embedding parameters -- CHECK!
  integer                      :: pub_freeze_switch_steps
  integer                      :: pub_ngwf_regions_ngroups ! modified in rundat_blocks_mod
  integer,allocatable          :: pub_ngwf_regions_group_nsp(:) ! modified in rundat_blocks_mod
  character(len=4), allocatable:: pub_ngwf_regions_groups(:,:) ! modified in rundat_blocks_mod
  ! rc2013: print out matrix information (for debugging)
  logical                      :: pub_embed_debug
  logical                      :: pub_quit_region = .false. ! modified in energy_and_force
  ! rc2013: control for projector embedding method
  logical                      :: pub_project_embed
  logical                      :: pub_do_fandt ! rc2013: do embedding calc
  logical                      :: pub_freeze_envir_ngwfs ! rc2013: do embedding calc


  ! ndmh: Calculation parameters not directly related to physics
  integer, protected       :: pub_fftbox_batch_size
  real(kind=DP), protected :: pub_dense_threshold
  logical, protected       :: pub_ovlp_for_nonlocal
  logical, protected       :: pub_use_space_filling_curve
  logical, protected       :: pub_coreham_denskern_guess
  logical, protected       :: pub_check_atoms
  character(len=80), protected :: pub_locpot_scheme
  real(kind=DP), protected :: pub_smooth_projectors
  real(kind=DP), protected :: pub_smooth_loc_pspot
  logical, protected       :: pub_odd_psinc_grid
  logical, protected       :: pub_even_psinc_grid
  logical, protected       :: pub_realspace_projectors
  logical, protected       :: pub_projectors_precalculate
  character(len=80), protected :: pub_parallel_scheme

  ! ndmh: Output/Visualisation parameters
  integer, protected           :: pub_output_detail
  integer, protected           :: pub_paw_output_detail
  integer, protected           :: pub_md_output_detail
  integer, protected           :: pub_timings_level
  character(len=80), protected :: pub_timings_order
  real(kind=DP), protected     :: pub_max_runtime
  logical, protected           :: pub_write_params
  logical, protected           :: pub_esdf_dump
  logical                      :: pub_write_forces  ! modified in ngwf_cg_mod
  logical, protected           :: pub_write_positions
  logical, protected           :: pub_write_velocities
  logical, protected           :: pub_write_xyz
  logical, protected           :: pub_print_qc
  logical, protected           :: pub_cube_format
  logical, protected           :: pub_dx_format
  logical, protected           :: pub_dx_coarse
  integer, protected           :: pub_dx_sig_digits
  logical, protected           :: pub_grd_format
  logical, protected           :: pub_write_density_plot
  logical, protected           :: pub_write_polarisation_plot
  logical, protected           :: pub_write_ngwf_plot
  logical, protected           :: pub_write_ngwf_grad_plot
  integer, protected           :: pub_write_ngwf_radial
  integer, protected           :: pub_write_ngwf_grad_radial
  real(kind=DP), protected     :: pub_write_radial_step
  real(kind=DP), protected     :: pub_write_radial_smear
  logical, protected           :: pub_write_initial_radial_ngwfs

  ! NGWF/kernel I/O parameters
  logical                      :: pub_read_denskern ! modified in energy_and_force
  logical                      :: pub_write_denskern ! modified in energy_and_force
  logical                      :: pub_read_sub_denskern ! modified in energy_and_force
  logical                      :: pub_write_tightbox_ngwfs ! modified in energy_and_force
  logical                      :: pub_read_tightbox_ngwfs ! modified in energy_and_force
  logical, protected           :: pub_write_converged_dk_ngwfs
  logical, protected           :: pub_write_sw_ngwfs
  logical                      :: pub_read_sw_ngwfs ! modified in md_mod
  integer, protected           :: pub_write_max_l
  integer                      :: pub_read_max_l ! modified in restart_mod
  integer, protected           :: pub_extra_n_sw
  logical                      :: pub_write_hamiltonian ! modified in forcetest_mod
  logical                      :: pub_read_hamiltonian ! modified in forcetest_mod
  logical, protected           :: pub_write_overlap
  logical, protected           :: pub_ngwf_cg_rotate

  ! ndmh: System parameters relating to grids, augmentation and pseudopotentials
  logical, protected           :: pub_paw
  real(kind=DP), protected     :: pub_fine_grid_scale
  real(kind=DP), protected     :: pub_dbl_grid_scale
  logical, protected           :: pub_aug_funcs_recip
  logical                      :: pub_usp         ! not input  ! modified in pseudopotentials_mod and paw_mod
  logical                      :: pub_aug         ! not input  ! modified in pseudopotentials_mod
  logical                      :: pub_any_nl_proj ! not input  ! modified in pseudopotentials_mod and energy_and_force
  logical                      :: pub_nlcc        ! not input  ! modified in pseudopotentials_mod and paw_mod
  logical                      :: pub_fine_is_dbl ! not input  ! modified in energy_and_force
  logical                      :: pub_dbl_is_std  ! not input  ! modified in energy_and_force
  integer                      :: pub_aug_den_dim ! not input  ! modified in energy_and_force
  logical                      :: pub_nhat_in_xc  ! not input  ! modified in pseudopotentials_mod and paw_mod
  character(len=500), protected :: pub_pseudo_path

  ! qoh: for dispersion correction
  character(len=80), protected :: pub_dispersion  ! damping function to use
  real(kind=DP), protected     :: pub_vdw_dcoeff  ! override of damping coefficient
  real(kind=DP), protected     :: pub_vdw_radial_cutoff ! radial cutoff for van der Waals interactions

  ! ddor: parameters for DFT+U with self-consistent projectors
  logical                  :: pub_hubbard ! modified in rundat_blocks_mod
  integer, protected       :: pub_hub_max_iter      ! Maximum number of DFT+U projector optimisation steps
  real(kind=DP), protected :: pub_hub_energy_tol    ! Energy tolerance when using DFT+U projector optimisation
  integer, protected       :: pub_hub_conv_win      ! Energy convergence window when using DFT+U projector optimisation
  real(kind=DP)            :: pub_hub_proj_mixing   ! Proportion of old Hubbard proj to mix with new ! modified in rundat_blocks_mod
  integer, protected       :: pub_hub_functional    ! DFT+U correction functional to use
  integer, protected       :: pub_hub_tensor_corr   ! DFT+U tensorial correction to use
  real(kind=DP), protected :: pub_hub_ngwf_spin_thr ! NGWF RMS gradient at which to switch off DFT+U spin-splitting
  logical, protected       :: pub_hub_on_the_fly    ! Carry out on-the-fly HUBBARDSCF with new projectors for new NGWFs
  logical                  :: pub_hubbard_restart   ! If restarting a HUBBARDSCF calc (not input) ! modified in rundat_blocks_mod
  logical                  :: pub_hubbard_atomsolve ! If using atomic guesses as DFT+U projs (not input) ! modif. in rundat_blocks
  logical, protected       :: pub_hub_calculating_u ! If calculating a bare static response (fixed Hamiltonian) (not input)
  logical, protected       :: pub_hub_read_projectors   ! Read hubbard-projectors from file (if .true.)
  logical, protected       :: pub_hubbard_compute_u_or_j  ! Compute U or J hubbard corrections
  logical, protected       :: pub_hubbard_j_minority_term  ! Include minority-only energy term in DFT+U+J
  logical, protected       :: pub_hub_tensor_forces    ! Calculate forces due to Hubbard metric tensor changing
  ! ebl: hubbard parameter needed for DMFT calculations
  logical                  :: pub_hub_proj_read_only

  ! ebl: parameters for DMFT
  logical           :: pub_dmft_fully_sc
  logical           :: pub_dmft_fully_sc_h
  integer           :: pub_dmft_kernel
  logical           :: pub_dmft_nbo
  integer           :: pub_dmft_nkpoints
  logical           :: pub_dmft_optics
  real(kind=DP)     :: pub_dmft_order_proj
  logical           :: pub_dmft_plot_real_space
  integer           :: pub_dmft_points
  logical           :: pub_dmft_purify_sc
  logical           :: pub_dmft_sc
  logical           :: pub_dmft_spoil_kernel
  logical           :: pub_dmft_switch_off_proj_order

  logical           :: pub_dmft_complex_freq
  real(kind=DP)     :: pub_dmft_cutoff_small
  real(kind=DP)     :: pub_dmft_dos_max
  real(kind=DP)     :: pub_dmft_dos_min
  real(kind=DP)     :: pub_dmft_emax
  real(kind=DP)     :: pub_dmft_emin
  real(kind=DP)     :: pub_dmft_kernel_mix
  logical           :: pub_dmft_kpoints_sym
  logical           :: pub_dmft_ks_shift
  real(kind=DP)     :: pub_dmft_mu_diff_max
  integer           :: pub_dmft_mu_order
  integer           :: pub_dmft_nmu_loop
  integer           :: pub_dmft_nval
  integer           :: pub_dmft_optics_i1
  integer           :: pub_dmft_optics_i2
  real(kind=DP)     :: pub_dmft_optics_window
  logical           :: pub_dmft_paramagnetic
  logical           :: pub_dmft_rotate_green
  real(kind=DP)     :: pub_dmft_scaling_cutoff
  integer           :: pub_dmft_scaling_meth
  integer           :: pub_dmft_scaling_nmpi
  real(kind=DP)     :: pub_dmft_scaling_tail
  logical           :: pub_dmft_skip_energy
  real(kind=DP)     :: pub_dmft_smear
  real(kind=DP)     :: pub_dmft_smear_eta
  real(kind=DP)     :: pub_dmft_smear_shift
  real(kind=DP)     :: pub_dmft_smear_T
  real(kind=DP)     :: pub_dmft_smear_w
  real(kind=DP)     :: pub_dmft_temp
  real(kind=DP)     :: pub_dmft_win
  logical           :: pub_dmft_write
  logical           :: pub_dmft_read

  ! ebl: parameters for DMFT that are currently devel_code
  ! real(kind=DP)     :: pub_dmft_chem_shift
  ! real(kind=DP)     :: pub_dmft_cutoff_tail
  ! real(kind=DP)     :: pub_dmft_doping
  ! integer           :: pub_dmft_embed_iter
  ! real(kind=DP)     :: pub_dmft_embed_mix
  ! real(kind=DP)     :: pub_dmft_free_green_frequ
  ! integer           :: pub_dmft_gpu_num
  ! logical           :: pub_dmft_impose_same_coeffs
  ! logical           :: pub_dmft_impose_chem_spin
  ! logical           :: pub_dmft_in_bohr
  ! logical           :: pub_dmft_integrate_green
  ! logical           :: pub_dmft_invert_overlap
  ! logical           :: pub_dmft_kpoints_kernel_gamma
  ! logical           :: pub_dmft_lin_scaling
  ! logical           :: pub_dmft_local_scratch
  ! integer           :: pub_dmft_norm_proj
  ! real(kind=DP)     :: pub_dmft_optics_x1
  ! real(kind=DP)     :: pub_dmft_optics_y1
  ! real(kind=DP)     :: pub_dmft_optics_z1
  ! logical           :: pub_dmft_plot_all_proj
  ! logical           :: pub_dmft_plot_real_space_sigma
  ! real(kind=DP)     :: pub_dmft_scaling_cutoff_h
  ! integer           :: pub_dmft_scaling_iter
  ! integer           :: pub_dmft_scaling_maxspace
  ! real(kind=DP)     :: pub_dmft_scaling_tol
  ! logical           :: pub_dmft_split
  ! logical           :: pub_dmft_splitk

  !gom: parameters for DFT+nu CG
  logical                  :: pub_turn_off_hartree     ! Omit Hartree terms in energy and potential
  logical                  :: pub_hubbard_unify_sites  ! Combine all projectors into one Hubbard site
  logical                  :: pub_dft_nu               ! Activate DFT+nu mode
  logical, protected       :: pub_dft_nu_opt_u1_only   ! Optimise onlyU1 potential
  logical, protected       :: pub_dft_nu_opt_u2_only   ! Optimise onlyU2 potential
  logical, protected       :: pub_dft_nu_continuation  ! RestartDFT_nu optimisation

  logical, protected       :: pub_turn_off_ewald ! jd: Omit Ewald energy & force
                                                 !     Saves time in pub_task properties in huge systems

  ! mjsp: parameters for EDA
  logical, protected   :: pub_eda              ! if this is an EDA type calculation
  logical, protected   :: pub_eda_continuation ! if this is an EDA calculation continuation
  logical              :: pub_eda_continuation_loaded_dkn   ! if EDA continuation data has been loaded
  logical              :: pub_eda_continuation_loaded_ngwfs ! if EDA continuation data has been loaded
  logical, protected   :: pub_eda_write        ! if EDA calculation continuation data to be written
  logical, protected   :: pub_eda_preptool     ! if this is an EDA supermolecule preparation task
  integer :: num_supermolecules
  logical, protected   :: pub_eda_read_frags   ! whether to load in fragment      dkn and NGWFs from file
  logical, protected   :: pub_eda_read_super   ! whether to load in supermolecule dkn and NGWFs from file
  character(len=80), allocatable    :: pub_frag_file_prefix(:) ! fragment filename prefixes
                                                         ! e.g. 'frag1' for frag1.dkn and frag1.tightbox_ngwfs
  character(len=80)                 :: pub_super_file_prefix ! supermolecule filename prefix
  integer              :: pub_eda_mode
  integer              :: pub_frag_counter
  integer              :: pub_frag_counter2
  logical              :: pub_eda_scfmi     ! if this is an SCF MI calculation (as explicit task or as part of EDA calculation)
  logical              :: pub_eda_scfmi_any ! if this calculation involves (at any point) an SCF MI calculation
  logical              :: pub_eda_have_sp_frags   ! if we have any spin polarised fragments (updated after fragment calculations)
  logical, protected   :: pub_eda_deltadens       ! Perform delta density calculations on the EDA energy components
  logical, protected   :: pub_eda_frag_isol_pol   ! Perform isolated fragment polarisation analysis as part of EDA
  logical, protected   :: pub_eda_frag_isol_ct    ! Perform fragment-pair delocalisation analysis as part of EDA
  integer, allocatable :: pub_frag_deloc(:,:)     ! Fragments to calculate delocalisations for
  logical, protected   :: pub_eda_reset_ngwfs_pol ! Resets the NGWFs after frozen stage of the EDA calculation
  logical, protected   :: pub_eda_reset_ngwfs_ct  ! Resets the NGWFs after polarisation stage of the EDA calculation
!  logical              :: pub_eda_scfmi_restrict_ngwfs ! Restrict NGWFs during polarisation calculation (in development)
  integer, allocatable       :: pub_frag_iatm(:)  ! Number of atoms in each fragment
  real(kind=DP), allocatable :: pub_frag_charge(:)  ! Charge on this fragment
  logical                    :: pub_eda_split_atoms ! Whether to partition fragments
  integer, allocatable       :: pub_frag_atoms(:,:) ! Atoms used to partition fragments
  real(kind=DP), allocatable :: pub_eda_isol(:)     ! Isolated fragment energies (idx=0 is total)
  real(kind=DP), allocatable :: pub_eda_isol_x(:)   ! Intrafragmental exchange
  real(kind=DP), allocatable :: pub_eda_isol_c(:)   ! Intrafragmental correlation
  real(kind=DP)   :: pub_eda_frz_x      ! Intra+interfragmental exchange (FRZ)
  real(kind=DP)   :: pub_eda_frz_c      ! Intra+interfragmental correlation (FRZ)
  real(kind=DP)   :: pub_eda_rep_x      ! Intra+interfragmental exchange (REP)
  real(kind=DP)   :: pub_eda_rep_c      ! Intra+interfragmental correlation (REP)
  real(kind=DP)   :: pub_eda_intrafrag_x_frz      ! Intrafragmental dE exchange contrib (isol-frz)
  real(kind=DP)   :: pub_eda_intrafrag_c_frz      ! Intrafragmental dE correlation contrib (isol-frz)
  real(kind=DP)   :: pub_eda_dE_x_rep      ! Intra+interfragmental dE exchange contrib (rep-frz)
  real(kind=DP)   :: pub_eda_dE_c_rep      ! Intra+interfragmental dE correlation contrib (rep-frz)
  integer         :: pub_eda_nodiag     ! Avoid diagonalisation (1=Frozen, 2=Polarisation, 3=Frozen and polarisation)

  ! ep: parameters for MERMIN MODE
  logical, protected           :: pub_mermin
  logical, protected           :: pub_check_mermin
  logical, protected           :: pub_mermin_temp
  real(kind=DP), protected     :: pub_mermin_smearing_width
  integer      , protected     :: pub_mermin_cheb
  logical, protected           :: pub_mermin_mu_sq
  integer                      :: maxit_mermin
  integer                      :: pub_firstit_mermin
  real(kind=DP)                :: mermin_threshold_orig
  character(len=80), protected :: pub_mermin_cg_type
  real(kind=DP), protected     :: pub_mermin_cg_max_step
  logical, protected           :: pub_mermin_check_trial_steps
  integer, protected           :: pub_mermin_round_evals
  real(kind=DP), protected     :: pub_mermin_free_energy_thres
  ! ep: parameters for MERMIN MODE

  ! gibo: parameters for Constrained_DFT (cDFT) ======= START
  logical            :: pub_cdft                          ! Constrained DFT run (if .true.) ! modified in rundat_blocks_mod
  logical, protected :: pub_ci_cdft                       ! Configuration-Interaction cDFT run (if .true.)
  logical, protected :: pub_cdft_hubbard                  ! Constrained DFT+U run (if .true.)
  logical, protected :: pub_cdft_atom_charge              ! ATOM-CHARGE-constrained CDFT simulation (if .true.)
  logical, protected :: pub_cdft_atom_spin                ! ATOM-SPIN-constrained CDFT simulation (if .true.)
  logical, protected :: pub_cdft_group_charge_diff        ! GROUP-CHARGE-DIFFERENCE-constrained CDFT (if .true.)
  logical, protected :: pub_cdft_group_spin_diff          ! GROUP-SPIN-DIFFERENCE-constrained CDFT (if .true.)
  logical, protected :: pub_cdft_group_charge_acceptor    ! GROUP-CHARGE-ACCEPTOR constrained CDFT simulation (if .true.)
  logical, protected :: pub_cdft_group_charge_donor       ! GROUP-CHARGE-DONOR constrained CDFT simulation (if .true.)
  logical, protected :: pub_cdft_group_spin_acceptor      ! GROUP-SPIN-ACCEPTOR constrained CDFT simulation (if .true.)
  logical, protected :: pub_cdft_group_spin_donor         ! GROUP-SPIN-DONOR constrained CDFT simulation (if .true.)
  logical, protected :: pub_cdft_print_all_occ            ! Print occupancies for all cDFT-atoms (if .true.)
  logical, protected :: pub_cdft_read_projectors          ! Read cDFT-projectors from file (if .true.)
  logical, protected :: pub_cdft_group_charge_up_only     ! logical to constrain only UP-electrons
  logical, protected :: pub_cdft_group_charge_down_only   ! logical to constrain only DOWN-electrons
  logical, protected :: pub_cdft_multi_proj               ! logical to use multiple-ang.mom projector on one cDFT-site
  logical, protected :: pub_cdft_write_potentials         ! Write cDFT-potentials into (.cdft) file (if .true.)

  integer, protected :: pub_ci_cdft_num_conf              ! Number of cDFT configurations for CI_CDFT run
  ! number of active cdft-modes [1: Charge OR Spin only, 2: Charge AND Spin]
  integer, protected :: pub_cdft_modes
  integer, protected :: pub_maxit_cdft_u_cg                   ! maximum number of cDFT-U CG iterations

  real(kind=DP)                :: pub_cdft_max_grad                 ! Maximum cDFT U-gradient component    ! modified in ngwf_cg_mod
  ! Targeted (acceptor-donor) CHARGE difference for  GROUP-CHARGE-DIFFERENCE-constrained cDFT
  real(kind=DP), protected     :: pub_cdft_group_charge_diff_target
  ! Targeted (acceptor-donor) SPIN difference for  GROUP-SPIN-DIFFERENCE-constrained cDFT
  real(kind=DP), protected     :: pub_cdft_group_spin_diff_target
  ! Constraining potential (eV) for GROUP-CHARGE-DIFFERENCE-constrained CDFT simulation
  real(kind=DP)                :: pub_cdft_group_charge_diff_u  ! modified in cdft_intermediate_cg_mod
  ! Constraining potential (eV) for GROUP-SPIN-DIFFERENCE-constrained CDFT simulation
  real(kind=DP)                :: pub_cdft_group_spin_diff_u ! modified in cdft_intermediate_cg_mod
  ! Targeted group-CHARGE for GROUP-CHARGE-ACCEPTOR-constrained cDFT
  real(kind=DP), protected     :: pub_cdft_group_charge_acceptor_target
  ! Targeted group-CHARGE for GROUP-CHARGE-DONOR-constrained cDFT
  real(kind=DP), protected     :: pub_cdft_group_charge_donor_target
  ! Targeted group-MAGMOM for GROUP-SPIN-ACCEPTOR-constrained cDFT
  real(kind=DP), protected     :: pub_cdft_group_spin_acceptor_target
  ! Targeted group-MAGMOM for GROUP-SPIN-DONOR-constrained cDFT
  real(kind=DP), protected     :: pub_cdft_group_spin_donor_target
  ! Constraining potential (eV) for GROUP-CHARGE-ACCEPTOR-constrained CDFT simulation
  real(kind=DP)                :: pub_cdft_group_charge_acceptor_u ! modified in cdft_intermediate_cg_mod
  ! Constraining potential (eV) for GROUP-CHARGE-DONOR-constrained CDFT simulation
  real(kind=DP)                :: pub_cdft_group_charge_donor_u ! modified in cdft_intermediate_cg_mod
  ! Constraining potential (eV) for GROUP-SPIN-ACCEPTOR-constrained CDFT simulation
  real(kind=DP)                :: pub_cdft_group_spin_acceptor_u ! modified in cdft_intermediate_cg_mod
  ! Constraining potential (eV) for GROUP-SPIN-DONOR-constrained CDFT simulation
  real(kind=DP)                :: pub_cdft_group_spin_donor_u ! modified in cdft_intermediate_cg_mod

  real(kind=DP), protected     :: pub_cdft_trial_length              ! initial trial_length for cDFT-line search
  real(kind=DP)                :: pub_cdft_cg_threshold              ! U-opt RMS gradient convg threshold ! modified in ngwf_cg_mod
  character(len=80), protected :: pub_cdft_cg_type                   ! type for CG CDFT U-optimisation
  integer, protected           :: pub_cdft_cg_max                ! max number iterations before U-opt CG reset
  real(kind=DP), protected     :: pub_cdft_elec_energy_tol       ! Max energy/atom change per cDFT iteration
  real(kind=DP), protected     :: pub_cdft_cg_max_step               ! Max trial-step in CG U-opt

  logical, protected           :: pub_cdft_guru                  ! The user is a cDFT_guru, let her/him free to mess around...
  logical, protected           :: pub_cdft_continuation          ! Continuate a previous cDFT-optimisation
  logical, protected           :: pub_cdft_tight                 ! Activate extra NGWFs-cDFT tight optimisation
  ! gibo: parameters for Constrained_DFT (cDFT) ======= END

  ! ndmh: parameters for Correction of PBCs
  character(len=80), protected :: pub_coulomb_cutoff_type
  real(kind=DP), protected     :: pub_coulomb_radius
  real(kind=DP), protected     :: pub_coulomb_length
  logical, protected           :: pub_coulomb_cutoff
  logical, protected           :: pub_coulomb_cutoff_write_int
  real(kind=DP), protected     :: pub_mt_cutoff

  ! ndmh: Properties calculation parameters
  logical                      :: pub_do_properties ! modified in onetep.F90
  integer, protected           :: pub_num_eigenvalues
  integer                      :: pub_homo_dens_plot ! modified in conduction_properties_mod
  integer                      :: pub_lumo_dens_plot ! modified in conduction_properties_mod
  integer                      :: pub_homo_plot ! modified in conduction_properties_mod
  integer                      :: pub_lumo_plot ! modified in conduction_properties_mod
  real(kind=DP), protected     :: pub_dos_smear
  integer                      :: pub_ldos_ngroups ! modified in rundat_blocks_mod
  logical, protected           :: pub_popn_calculate
  logical, protected           :: pub_lowdin_popn_calculate
  real(kind=DP), protected     :: pub_popn_bond_cutoff
  logical, protected           :: pub_popn_mulliken_partial
  logical, protected           :: pub_ngwf_analysis
  logical, protected           :: pub_polarisation_calculate
  logical, protected           :: pub_polarisation_local
  logical, protected           :: pub_polarisation_berry
  logical, protected           :: pub_polarisation_simcell_calculate
  real(kind=DP), protected     :: pub_polarisation_simcell_refpt(3)
  logical, protected           :: pub_efield_calculate
  logical, protected           :: pub_spread_calculate
  integer, allocatable         :: pub_ldos_group_nsp(:) ! not input ! modified in rundat_blocks_mod
  character(len=4), allocatable:: pub_ldos_groups(:,:) ! modified in rundat_blocks_mod
  logical, protected           :: pub_anharmonic_calculate
  logical, protected           :: pub_print_potential_noxc ! print Ion+Hartree (NO XC) local potential

  logical, protected           :: pub_pdos_optados_output
  integer, protected           :: pub_pdos_max_l
  integer                      :: pub_pdos_ngroups ! modified in rundat_blocks_mod
  integer, allocatable         :: pub_pdos_group_nsp(:) ! not input ! modified in rundat_blocks_mod
  character(len=4), allocatable:: pub_pdos_groups(:,:) ! modified in rundat_blocks_mod
  logical, protected           :: pub_pdos_reduce_sws
  integer, protected           :: pub_pdos_max_n
  integer, protected           :: pub_pdos_stride_n
  logical, protected           :: pub_pdos_construct_basis
  logical, protected           :: pub_pdos_lowdin
  logical, protected           :: pub_pdos_lcao_optimize
  logical, protected           :: pub_pdos_orth_atom_blocks
  logical, protected           :: pub_pdos_output_basis
  logical, protected           :: pub_pdos_output_swopt_kernham
  logical, protected           :: pub_pdos_pseudoatomic
  logical, protected           :: pub_pdos_sum_mag


  ! aam: geometry optimiser parameters
  real(kind=DP)                :: pub_geom_modulus_est ! modified in geometry_optimiser_mod
  real(kind=DP)                :: pub_geom_frequency_est ! modified in geometry_optimiser_mod
  real(kind=DP), protected     :: pub_geom_energy_tol
  real(kind=DP), protected     :: pub_geom_force_tol
  real(kind=DP), protected     :: pub_geom_disp_tol
  integer, protected           :: pub_geom_max_iter
  integer                      :: pub_geom_convergence_win ! modified in geometry_optimiser_mod
  integer, protected           :: pub_geom_backup_iter
  integer, protected           :: pub_geom_reset_dk_ngwfs_iter
  integer, protected           :: pub_geom_lbfgs_max_updates
  integer, protected           :: pub_geom_lbfgs_block_length
  logical, protected           :: pub_geom_continuation
  logical, protected           :: pub_geom_print_inv_hessian
  logical, protected           :: pub_geom_reuse_dk_ngwfs
  character(len=80)            :: pub_geom_method ! modified in geometry_optimiser_mod
  logical, protected           :: pub_geom_lbfgs
  ! lk: added pub_geom_output_detail variable
  !     to control the level of ouput detail
  !     of geometry optimization.
  integer, protected           :: pub_geom_output_detail
  ! lam81: preconditioners
  character(len=80)            :: pub_geom_precond_type           ! type of preconditioner
  logical, protected           :: pub_geom_precond_scale_cell     ! scaling cell
  real(kind=DP)                :: pub_geom_precond_exp_c_stab     ! stabilisation constant for Exp preconditioner
  real(kind=DP)                :: pub_geom_precond_exp_A          ! A value for Exp preconditioner
  real(kind=DP)                :: pub_geom_precond_exp_r_NN       ! nearest neighbor distance for Exp precond
  real(kind=DP)                :: pub_geom_precond_exp_r_cut      ! cutoff distance for Exp precond
  real(kind=DP)                :: pub_geom_precond_exp_mu         ! mu value for Exp precond
  real(kind=DP)                :: pub_geom_precond_ff_c_stab      ! stabilisation constant for FF preconditioner
  real(kind=DP)                :: pub_geom_precond_ff_r_cut       ! cutoff distance for FF precond

  ! ny: Species dependent scissor parameters
  integer                      :: pub_scissor_ngroups ! modified in rundat_blocks_mod
  integer, allocatable         :: pub_scissor_group_nsp(:) ! not input ! modified in rundat_blocks_mod
  character(len=4), allocatable:: pub_scissor_groups(:,:) ! modified in rundat_blocks_mod
  real(kind=DP), allocatable   :: pub_scissor_group_shifts(:,:) ! modified in rundat_blocks_mod

  ! vm: transition state search parameters
  character(len=8), protected  :: tssearch_method          ! transition state search method (LSTQST for now)
  character(len=20), protected :: tssearch_lstqst_protocol ! which protocol to use
  integer, protected           :: tssearch_qst_max_iter    !
  integer, protected           :: tssearch_cg_max_iter     !
  real (kind=dp), protected    :: tssearch_force_tol       ! force tolerance
  real (kind=dp), protected    :: tssearch_disp_tol        ! displacement tolerance
  real (kind=dp), protected    :: tssearch_energy_tol      ! energy tolerance

  ! kkbd: NEB-specific transition state search parameters
  integer :: pub_neb_ci_delay ! Number of NEB chain LBFGS iterations before we enable a climbing image.
                              ! Negative number to disable climbing image
  logical, protected :: pub_neb_print_summary ! if .true., print out convergence info and reaction pathway info to stdout
  real(kind=dp), protected :: pub_neb_spring_constant ! Spring constant for the NEB chain
  logical, protected :: pub_neb_read_xyz
  logical, protected :: pub_neb_converge_all
  character(len=80), protected :: pub_neb_update_method
  integer, protected :: pub_neb_max_iter
  integer, protected :: pub_neb_glbfgs_history_size
  logical, protected :: pub_neb_continuation

  real(kind=DP), protected :: pub_reactant_energy
  real(kind=DP), protected :: pub_product_energy
  character(len=80), protected :: pub_reactant_rootname
  character(len=80), protected :: pub_product_rootname

  ! aam: molecular dymanics
  real(kind=dp), protected     :: md_delta_t
  integer                      :: md_num_iter ! It is overwritten in md_mod
  integer                      :: md_iter_global ! Not a keyword
  logical                      :: md_restart ! Modified if md_global_restart=T
  logical                      :: md_restart_thermo ! Modified if md_global_restart=T
  logical                      :: md_global_restart ! Modified in
                                                    ! ngwfs_initialise
  logical                      :: md_write_out
  logical                      :: md_init_velocities ! modified in rundat_blocks_mod
  integer, protected           :: md_reset_history
  logical, protected           :: pub_md_properties
  integer                      :: md_write_history
  character(len=5), protected  :: md_aux_rep
  real(kind=dp), protected     :: md_aux_dkn_t
  real(kind=dp), protected     :: md_aux_beren_tc
  real(kind=dp), protected     :: md_lnv_threshold
  real(kind=dp)                :: md_ngwf_threshold
  ! vv
  logical, protected           :: md_autocorr

  logical, protected           :: mts_xi
  integer                      :: mts_nstep ! modified in md_mod
  real(kind=dp), protected     :: mts_delta_t
  real(kind=dp)                :: mts_ngwf_threshold ! modified in md_mod
  real(kind=dp)                :: mts_lnv_threshold ! modified in md_mod
  real(kind=DP)                :: mts_elec_energy_tol ! modified in md_mod
  real(kind=DP)                :: mts_elec_force_tol ! modified in md_mod
  real(kind=DP)                :: mts_ngwf_max_grad ! modified in md_mod
  integer                      :: mts_maxit_ngwf_cg ! modified in md_mod
  integer                      :: mts_minit_lnv ! modified in md_mod
  integer                      :: mts_maxit_lnv ! modified in md_mod
  integer                      :: mts_maxit_pen ! modified in md_mod

  logical, protected           :: mts_mix_inc

  character(len=20)            :: mix_dkn_type ! modified in onetep.F90
  integer, protected           :: mix_dkn_num
  integer, protected           :: mix_dkn_reset
  character(len=20), protected :: mix_dkn_init_type
  integer, protected           :: mix_dkn_init_num

  character(len=20), protected :: mix_ngwfs_type
  integer, protected           :: mix_ngwfs_num
  integer, protected           :: mix_ngwfs_reset
  character(len=20), protected :: mix_ngwfs_init_type
  integer, protected           :: mix_ngwfs_init_num

  real(kind=dp), protected     :: mix_ngwfs_coeff
  real(kind=dp), protected     :: mix_local_length
  real(kind=dp), protected     :: mix_local_smear

  ! lr408: Conduction parameters
  logical                      :: pub_cond_calculate ! modified in rundat_blocks
  logical                      :: pub_cond_calculate_any_task ! modified in rundat_blocks
  logical                      :: cond_read_denskern ! modified in conduction_mod
  logical, protected           :: cond_read_tightbox_ngwfs
  logical, protected           :: cond_fixed_shift
  logical, protected           :: cond_calc_max_eigen
  integer, protected           :: cond_num_states
  integer, protected           :: cond_maxit_lnv
  integer, protected           :: cond_minit_lnv
  integer, protected           :: cond_firstit_lnv
  real(kind=DP), protected     :: cond_kernel_cutoff
  real(kind=DP), protected     :: cond_init_shift
  real(kind=DP), protected     :: cond_shift_buffer
  integer                      :: cond_num_extra_states ! modified in conduction_mod
  integer, protected           :: cond_num_extra_its
  logical, protected           :: cond_plot_joint_orbitals
  logical, protected           :: cond_plot_vc_orbitals
  real(kind=DP), protected     :: cond_energy_range
  real(kind=DP), protected     :: cond_energy_gap

  ! gcc32: LR_PHONONS PARAMETERS
  logical                      :: pub_lr_phonons_calculate
  logical                      :: pub_lr_phonons_zero_dim
  logical                      :: pub_lr_phonons_restart
  real(kind=DP), protected     :: pub_lr_phonons_kernel_cutoff

  ! tjz07: LR_TDDFT parameters
  logical                      :: pub_lr_tddft_calculate ! modified in onetep.F90
  integer, protected           :: pub_lr_tddft_num_states
  real(kind=DP), protected     :: pub_lr_tddft_cg_threshold
  integer, protected           :: pub_lr_tddft_maxit_cg
  integer, protected           :: pub_lr_tddft_maxit_pen
  integer, protected           :: pub_lr_tddft_reset_cg
  real(kind=DP), protected     :: pub_lr_tddft_penalty_tol
  real(kind=DP), protected     :: pub_lr_tddft_kernel_cutoff
  logical, protected           :: pub_lr_tddft_write_densities
  logical, protected           :: pub_lr_tddft_write_kernels
  logical, protected           :: pub_lr_tddft_restart
  logical, protected           :: pub_lr_tddft_triplet
  logical, protected           :: pub_lr_tddft_projector
  logical, protected           :: pub_lr_tddft_joint_set
  logical                      :: pub_lr_tddft_preopt ! modified in lr_tddft_mod
  integer                      :: pub_lr_tddft_preopt_iter ! modified in lr_tddft_mod
  logical, protected           :: pub_lr_tddft_analysis
  integer                      :: pub_lr_tddft_num_conv_states ! modified in lr_tddft_mod
  integer, protected           :: pub_lr_tddft_check_conv_iter
  logical, protected           :: pub_lr_tddft_precond
  real(kind=DP), protected     :: pub_lr_tddft_precond_tol
  integer, protected           :: pub_lr_tddft_precond_iter
  logical, protected           :: pub_lr_tddft_RPA
  logical, protected           :: pub_lr_tddft_init_random
  logical, protected           :: pub_lr_tddft_init_max_overlap
  logical, protected           :: pub_lr_tddft_mlwf_analysis
  logical, protected           :: pub_lr_tddft_mom_mat_els
  logical, protected           :: pub_lr_tddft_penalty_func
  logical                      :: pub_lr_tddft_sparse_region ! defines   ! modified in rundat_blocks_mod
                       ! whether sparsity region for kernel is used
  logical, protected           :: pub_lr_tddft_properties_only
  real(kind=DP), protected     :: pub_lr_tddft_ct_length
  integer                      :: pub_tddft_kernel_ngroups ! modified in rundat_blocks_mod
  integer,allocatable          :: pub_tddft_kernel_group_nsp(:) ! modified in rundat_blocks_mod
  character(len=4), allocatable:: pub_tddft_kernel_groups(:,:) ! modified in rundat_blocks_mod
  integer                      :: pub_tddft_ct_ngroups ! modified in rundat_blocks_mod
  integer,allocatable          :: pub_tddft_ct_group_nsp(:) ! modified in rundat_blocks_mod
  character(len=4), allocatable:: pub_tddft_ct_groups(:,:) ! modified in rundat_blocks_mod
  integer, protected           :: pub_lr_tddft_HOMO_num
  integer, protected           :: pub_lr_tddft_LUMO_num
  logical, protected           :: pub_lr_tddft_xc_finite_diff
  logical, protected           :: pub_lr_tddft_subsystem_coupling
  logical, protected           :: pub_lr_tddft_restart_from_TDA
  real(kind=DP)                :: pub_lr_optical_permittivity ! modified in lr_tddft_mod
  real(kind=DP)                :: pub_lr_tddft_spectrum_smear
  real(kind=DP), private       :: lr_tddft_default_permittivity ! only defined in here

  ! jhl52: QNTO parameters
  logical, protected           :: pub_qnto_analysis
  integer, protected           :: pub_qnto_num_transition
  logical, protected           :: pub_qnto_write_orbitals
  integer, protected           :: pub_qnto_svd_method
  logical, protected           :: pub_qnto_nbo_proj
  character(len=80), protected :: pub_qnto_ref_dir
  integer, protected           :: pub_qnto_num_ref_states
  integer, protected           :: pub_qnto_num_core_atoms

  ! smmd: Transmission coefficients
  logical, protected           :: pub_etrans_bulk
  logical, protected           :: pub_etrans_lcr
  logical, protected           :: pub_etrans_same_leads
  logical, protected           :: pub_etrans_write_xyz
  logical, protected           :: pub_etrans_write_hs
  real(kind=dp), protected     :: pub_etrans_ecmplx

  real(kind=dp), protected     :: pub_etrans_emax
  real(kind=dp), protected     :: pub_etrans_emin
  integer, protected           :: pub_etrans_enum
  real(kind=dp), protected     :: pub_etrans_eref
  character(len=80), protected :: pub_etrans_eref_method
  logical, protected           :: pub_etrans_calc_lead_pot
  integer, protected           :: pub_etrans_lead_nkpoints
  real(kind=dp), protected     :: pub_etrans_lead_disp_tol
  integer, protected           :: pub_etrans_num_eigchan
  logical, protected           :: pub_etrans_plot_eigchan
  logical, protected           :: pub_etrans_lead_lcheck
  integer, protected           :: pub_etrans_seed_lead

  ! lr408: Conduction / optical spectra parameters
  logical, protected :: pub_spectra_calculate
  ! lr408: if true calculate momentum matrix elements for spectra,
  ! lr408: otherwise calculate position matrix elements  (for molecules)
  logical, protected :: pub_calc_mom_mat_els
  logical, protected :: pub_calc_nonloc_comm
  logical, protected :: pub_spec_cont_deriv
  real(kind=dp), protected :: pub_spec_nonloc_fin_diff
  logical, protected :: pub_spectra_print_mat_els
  real(kind=dp), protected :: pub_spec_scissor_op
  real(kind=dp), protected :: pub_opt_smear

  logical, protected :: pub_eels_calculate
  logical, protected :: pub_eels_fine_projectors
  logical, protected :: pub_eels_realspace
  ! fc: phonon parameters
  logical, protected :: pub_have_phonon_disp_list, pub_have_phonon_except_list
  logical, protected :: pub_phonon_energy_check, pub_phonon_write_eigenvecs
  logical, protected :: pub_phonon_DOS, pub_phonon_SK
  logical, protected :: pub_have_phonon_animate_list, pub_have_supercell
  logical, protected :: pub_have_phonon_grid, pub_have_phonon_qpoints
  integer, protected :: pub_num_disp, pub_num_except, pub_num_vib, pub_num_ani
  integer, protected :: pub_phonon_farming_task, pub_phonon_sampling
  integer, protected :: pub_phonon_vib_free, pub_supercell(1:3), pub_phonon_grid(1:3)
  integer, protected :: pub_nat_unit, pub_num_grid, pub_num_qpoints
  integer, allocatable :: pub_phonon_disp_list(:)         ! Removed protection
  integer, allocatable :: pub_phonon_animate_list(:)      !     as these are
  integer, allocatable :: pub_phonon_iexcept_list(:,:)    !     deallocated in
  integer, allocatable :: pub_supercell_unit_list(:)      !     phonon_main
  real(kind=DP), allocatable :: pub_phonon_qpoints(:,:)   ! [jd, retreat2015]
  real(kind=DP), allocatable :: pub_phonon_dexcept_list(:)!
  real(kind=DP), protected :: pub_phonon_disp, pub_phonon_fmax
  real(kind=DP), protected :: pub_phonon_tmin, pub_phonon_tmax
  real(kind=DP), protected :: pub_phonon_deltat, pub_phonon_min_freq
  real(kind=DP), protected :: pub_phonon_animate_scale
  real(kind=DP), protected :: pub_phonon_DOS_min, pub_phonon_DOS_max
  real(kind=DP), protected :: pub_phonon_DOS_delta

  ! vv: anharmonic parameters
  character(len=80), protected :: pub_anh_qc_factor
  character(len=80), protected :: pub_anh_acf_factor
  integer, protected           :: pub_anh_first_iter
  integer, protected           :: pub_anh_last_iter
  real(kind=DP), protected     :: pub_anh_plot_firstfreq
  real(kind=DP), protected     :: pub_anh_plot_lastfreq
  real(kind=DP), protected     :: pub_anh_md_temp
  logical, protected           :: pub_anh_plot_all
  logical, protected           :: pub_anh_apply_filter
  character(len=80), protected :: pub_anh_type

  ! pdh: bandstructure parameters
  integer, protected           :: pub_bs_kpoint_path_length
  integer                      :: pub_bs_num_eigenvalues ! modified in bandstructure_mod
  integer                      :: pub_bs_unfold(3) ! modified in bandstructure_mod
  logical, protected           :: pub_do_bandstructure
  character(len=80), protected :: pub_bs_method
  real(kind=DP), protected     :: pub_bs_kpoint_path_spacing
  real(kind=DP), allocatable, protected :: pub_bs_kpoint_path_start(:,:)
  real(kind=DP), allocatable, protected :: pub_bs_kpoint_path_end(:,:)
  logical, protected           :: pub_perturbative_soc

  ! gcc32: bandstructure unfolding parameters
  integer, protected :: pub_bsunfld_kpoint_path_length
  integer, protected :: pub_bsunfld_num_kpts_path
  integer, protected :: pub_bsunfld_num_atoms_prim
  logical, protected :: pub_bsunfld_restart
  logical, protected :: pub_bsunfld_calculate
  integer            :: pub_bsunfld_num_eigenvalues !modif in bandstructure_mod
  real(kind=DP), protected :: pub_bsunfld_transfo(9)
  real(kind=DP), allocatable, protected :: pub_bsunfld_kpoint_path_start(:,:)
  real(kind=DP), allocatable, protected :: pub_bsunfld_kpoint_path_end(:,:)
  integer :: pub_bsunfld_ngroups, pub_bsunfld_nprojatoms ! modified in
                                                         ! rundat_blocks_mod
  integer, allocatable :: pub_bsunfld_group_nsp(:),pub_bsunfld_projatoms_nsp(:)
                                                   ! not input, modified in
                                                   ! rundat_blocks_mod
  character(len=4), allocatable:: pub_bsunfld_groups(:,:), &
      pub_bsunfld_projatoms(:,:) ! modified in rundat_blocks_mod


  ! agreco: k-points parameters for BZ sampling in scf calculation
  ! agreco: copy of pub_num_kpoints, to be merged when
  ! everything is working properly for multiple k-points
  integer :: pub_num_kpoints_temp
  ! agreco: whether the k-point list for BZ sampling is specified
  ! in the input file
  logical, protected :: pub_kpoints_specified_list
  ! agreco: array to store k-points with specified weights
  real(kind=DP), allocatable, protected :: pub_kpoint_list(:,:)
  ! agreco: real variable to store k-point weight as in input file
  real(kind=DP) :: kpoint_weight
  ! agreco: global flag whether we are doing Gamma point only or not
  logical :: pub_gamma_point_only

  ! jd: --- SWRI & SWX ----
  ! jd: Batch sizes for products of Chebyshev expansions
  integer, protected :: pub_swri_cheb_batchsize
  ! jd: Apply SW/SWpot stencil smoothing for calculating ovlp on coarse grids?
  logical, protected :: pub_swri_swop_smoothing
  ! jd: Specifies how matrix inversions are to be done in SWX w/ overlap metric
  logical, protected :: pub_swri_overlap_indirect
  character(len=1024), protected :: pub_swri_assembly_prefix
  ! jd: Should SWRI output be extra verbose?
  logical, protected :: pub_swri_verbose
  real(kind=DP), protected     :: pub_swri_proximity_sort_point(3)
  ! jd: Should matrix inversions in SWX be improved upon with Hotelling?
  logical, protected :: pub_swri_improve_inverse
  ! jd: Print metric matrix eigenvalues to assess linear dependence of SWs?
  logical, protected :: pub_swri_print_eigenvalues
  real(kind=DP), protected :: pub_swx_c_threshold
  integer, protected       :: pub_swx_output_detail
  logical, protected       :: pub_swx_dbl_grid
  ! -----------------------

  ! JCW: --- SPH HARM ROTATION ----
  logical, protected       :: pub_use_sph_harm_rot
  ! -------------------------------

  ! qoh: Hartree-Fock exchange parameters:
  ! qoh: Whether Hartree-Fock exchange is being used. This is not an input
  ! jd:  parameter, and is set in xc_hfxinit, from energy_and_force_init_cell.
  ! jd:  It cannot be set here, as rundat would have to depend on xc.
  logical, protected :: pub_use_hfx ! set in xc_hfxinit, through rundat_set_pub_use_hfx
  logical, protected :: pub_hfx_nlpp_for_exchange
  ! jd: X matrix read/write for restarts
  logical, protected :: pub_hfx_read_xmatrix, pub_hfx_write_xmatrix
  ! jd: Electrostatic or overlap metric
  character(len=80), private   :: hfx_metric_string
  character(len=32), protected :: pub_hfx_use_ri
  real(kind=DP), protected     :: pub_hfx_cutoff
  integer, protected           :: pub_hfx_metric
  integer, protected           :: pub_hfx_max_l
  integer, protected           :: pub_hfx_max_q
  integer, protected           :: pub_hfx_output_detail
  integer, protected           :: pub_hfx_bessel_rad_nptsx
  integer, protected           :: pub_hfx_memory_limit
  real(kind=DP)                :: pub_hfx_memory_weights(3) ! adjusted in hfx
  logical, protected           :: pub_hfx_debug

  integer            :: pub_cache_limit_for_swops      ! these are
  integer            :: pub_cache_limit_for_expansions ! adjusted in
  integer            :: pub_cache_limit_for_prods      ! hf_exchange_mod
  integer, protected :: pub_cache_limit_for_dknblks
  integer, protected :: pub_cache_limit_for_swops2

  integer, protected :: pub_ht_stash_size

  ! jd: If true (default), total ionic force is zeroed by subtracting the
  !     average force from the force experienced by (non-classical) atoms
  logical            :: pub_zero_total_force
  ! vv: If true, ionic forces are corrected with a mass-weighted factor from
  !     quantum atoms only
  logical            :: pub_mw_total_force

  ! qoh: Whether we are using tightboxes to do FFTs
  logical, protected :: pub_tightbox_fft_coarse ! (not input)

  ! ddor: TDDFT flags and parameters
  logical, protected           :: pub_do_tddft
  real(kind=DP), protected     :: pub_tddft_maximum_energy
  real(kind=DP), protected     :: pub_tddft_resolution
  character(len=20), protected :: pub_tddft_propagation_method
  integer, protected           :: pub_tddft_sparsity_level
  logical, protected           :: pub_tddft_tammdancoff
  real(kind=DP), protected     :: pub_tddft_dipole_kick_strength(3)
  character(len=80)            :: pub_tddft_xc_functional ! modified in xc_mod
  integer, protected           :: pub_tddft_hamiltonian_mixing
  real(kind=DP), protected     :: pub_tddft_damping
  logical, protected           :: pub_tddft_enforced_idempotency
  integer, protected           :: pub_tddft_maxit_hotelling
  real(kind=DP), protected     :: pub_tddft_max_resid_hotelling
  logical, protected           :: pub_tddft_inv_overlap_exact

  ! jd: For QM/MM
  logical, protected           :: pub_pol_emb_pot ! not a keyword
  logical, protected           :: pub_pol_emb_qmstar
  logical, protected           :: pub_pol_emb_vacuum_qmstar
  logical, protected           :: pub_pol_emb_write_vacuum_restart
  character(len=80), protected :: pub_pol_emb_pot_filename
  real(kind=DP), protected     :: pub_pol_emb_mpole_exclusion_radius
  real(kind=DP), protected     :: pub_pol_emb_thole_a
  real(kind=DP), protected     :: pub_pol_emb_pairwise_polarisability
  real(kind=DP), protected     :: pub_pol_emb_smearing_a
  integer, protected           :: pub_pol_emb_dma_min_l
  integer, protected           :: pub_pol_emb_dma_max_l
  integer, protected           :: pub_pol_emb_vacuum_dma_min_l
  integer, protected           :: pub_pol_emb_vacuum_dma_max_l
  real(kind=DP), protected     :: pub_pol_emb_perm_scaling
  logical, protected           :: pub_pol_emb_fixed_charge
  logical, protected           :: pub_pol_emb_dbl_grid
  real(kind=DP), protected     :: pub_pol_emb_repulsive_mm_pot_a
  real(kind=DP), protected     :: pub_pol_emb_repulsive_mm_pot_b
  real(kind=DP), protected     :: pub_pol_emb_repulsive_mm_pot_c
  real(kind=DP), protected     :: pub_pol_emb_repulsive_mm_pot_alpha
  real(kind=DP), protected     :: pub_pol_emb_repulsive_mm_pot_beta
  real(kind=DP), protected     :: pub_pol_emb_repulsive_mm_pot_r0
  real(kind=DP), protected     :: pub_pol_emb_repulsive_mm_pot_cutoff
  logical, protected           :: pub_pol_emb_repulsive_mm_pot_write
  logical, protected           :: pub_pol_emb_repulsive_mm_pot_verbose
  real(kind=DP), protected     :: pub_pol_emb_polscal

  ! jd: for implicit solvent
  ! - Fattebert-Gygi-Scherlis model and MPSM
  real(kind=DP), protected     :: pub_is_density_threshold
  real(kind=DP), protected     :: pub_is_solvation_beta
  ! - Andreussi (SCCS) model
  real(kind=DP), protected     :: pub_is_density_min_threshold
  real(kind=DP), protected     :: pub_is_density_max_threshold
  ! - Soft-sphere model
  real(kind=DP), protected     :: pub_is_soft_sphere_delta
  real(kind=DP), protected     :: pub_is_soft_sphere_scale
  real(kind=DP)                :: pub_is_bulk_permittivity ! modified in lr_tddft_rpa_mod
  logical, protected           :: pub_is_multigrid_verbose
  real(kind=DP), protected     :: pub_is_multigrid_verbose_y
  real(kind=DP), protected     :: pub_is_multigrid_verbose_z
  integer, protected           :: pub_is_bc_coarseness
  integer, protected           :: pub_is_bc_surface_coarseness
  real(kind=DP), protected     :: pub_is_bc_threshold
  logical                      :: pub_is_bc_allow_frac_charge
  real(kind=DP), protected     :: pub_is_smeared_ion_width
  real(kind=DP), protected     :: pub_is_core_width
  real(kind=DP), protected     :: pub_is_surface_thickness
  real(kind=DP), protected     :: pub_is_solvent_surf_tension
  real(kind=DP), protected     :: pub_is_solvent_pressure
  real(kind=DP), protected     :: pub_is_apolar_scaling_factor
  logical                      :: pub_is_implicit_solvent ! modified in energy_and_force
  logical, protected           :: pub_is_check_solv_energy_grad
  logical, protected           :: pub_is_smeared_ion_rep
  logical                      :: pub_is_include_apolar ! modified in energy_and_force
  character(len=80), protected :: pub_is_apolar_sasa_definition
  character(len=80), protected :: pub_is_apolar_method
  logical, protected           :: pub_is_separate_restart_files
  logical                      :: pub_is_auto_solvation ! modified in energy_and_force
  real(kind=DP), protected     :: pub_is_sc_steric_magnitude
  real(kind=DP), protected     :: pub_is_sc_steric_cutoff
  real(kind=DP), protected     :: pub_is_sc_steric_smoothing_alpha
  character(len=80)            :: pub_is_steric_pot_type
  real(kind=DP), protected     :: pub_is_hc_steric_smearing
  real(kind=DP), protected     :: pub_is_hc_steric_dens_isovalue
  logical, protected           :: pub_is_steric_write
  character(len=80), protected :: pub_is_solvation_method
  character(len=80), protected :: pub_is_dielectric_model
  character(len=80), protected :: pub_is_dielectric_function
  character(len=80), protected :: pub_is_solvation_output_detail
  real(kind=DP), protected     :: pub_is_dielectric_exclusions_smearing
  logical, protected           :: pub_is_solvation_properties
  character(len=80)            :: pub_is_pbe ! modified in energy_and_force
  logical                      :: pub_is_pbe_bc_debye_screening ! modified in e&f
  real(kind=DP), protected     :: pub_is_pbe_energy_tolerance
  real(kind=DP), protected     :: pub_is_pbe_exp_cap
  real(kind=DP), protected     :: pub_is_pbe_temperature
  character(len=80)            :: pub_is_pbe_neutralisation_scheme ! modif. in check_inputs
  logical, protected           :: pub_is_restart_vac_from_vac
  logical, protected           :: pub_is_emft_cavity

  real(kind=DP), private :: is_default_bulk_permittivity ! jd: not used outside here
  real(kind=DP), private :: obsolete_solvent_surface_tension


  ! jd: for open BCs in local pseudo and ion-ion energies
  !   - is ion-ion energy calculated with open BC (by direct summation)
  logical, protected :: pub_ii_energy_direct = .false.
  !   - is hartree calculated with MG? (currently this implies open BC)
  logical, protected :: pub_multigrid_hartree = .false.
  !   - is multgrid used? (affects cell_grid_distribute)
  logical, protected :: pub_multigrid_in_use = .false.
  !   - is local pseudo calculated with open BC?
  logical, protected :: pub_open_localpseudo = .false.
  !   - did the user force open BCs in ion-ion energy?
  logical, protected :: pub_openbc_ion_ion
  !   - did the user force open BCs in hartree calculation?
  logical, protected :: pub_openbc_hartree
  !   - did the user force open BCs in local pseudo calculation?
  logical, protected :: pub_openbc_pspot
  !   - parameter 'npts_x' in the open BC local pseudo calculation
  integer, protected :: pub_openbc_pspot_finetune_nptsx
  !   - parameter 'alpha' in the open BC local pseudo calculation
  real(kind=DP), protected :: pub_openbc_pspot_finetune_alpha
  !   - parameter 'fineness' in the open BC local pseudo calculation
  integer, protected :: pub_openbc_pspot_finetune_f


  ! JCW: Updated parameters for controlling the multigrid solver
  ! JCW:
  ! JCW: (These may be a good candidate to consolidate into a block specifically
  ! JCW: for controlling the solver)
  ! JCW:
  ! JCW: These keywords control DL_MG using the new v2.0 API. Some of the old
  ! JCW: keywords which are MG-solver specific have been replaced with updated
  ! JCW: versions and new keywords have been added as new parameters are available
  ! JCW: for external control.
  integer, protected       :: pub_mg_granularity_power    ! replaces pub_is_multigrid_nlevels
  integer, protected       :: pub_mg_defco_fd_order ! replaces pub_is_discretization_order
                                                    ! for DL_MG's high-order FD
  real(kind=DP), protected :: pub_mg_tol_res_rel    ! relative tolerance in residual
  real(kind=DP), protected :: pub_mg_tol_res_abs    ! absolute tolerance in residual
    ! pub_mg_tol_res_{rel,abs} apply to the "top" iteration, which is defect correction if
    ! fd_order > 2, otherwise Newton or V-cycle. However, since DL_MG is always called from ONETEP
    ! with the Newton and V-cycle tolerances explicitly set, when fd_order == 2 the tolerances for
    ! the Newton and V-cycle iterations are determined by their specific keyword value.
  real(kind=DP), protected :: pub_mg_tol_pot_rel    ! relative tolerance in potential (defect corr)
  real(kind=DP), protected :: pub_mg_tol_pot_abs    ! absolute tolerance in potential (defect corr)
  real(kind=DP), protected :: pub_mg_tol_vcyc_rel   ! relative tolerance in V-cycle iterations
  real(kind=DP), protected :: pub_mg_tol_vcyc_abs   ! absolute tolerance in V-cycle iterations
  real(kind=DP), protected :: pub_mg_tol_newton_rel   ! relative tolerance in Newton iterations
  real(kind=DP), protected :: pub_mg_tol_newton_abs   ! absolute tolerance in Newton iterations
  real(kind=DP), protected :: pub_mg_max_res_ratio    ! threshold for giving up on convergence
  integer, protected       :: pub_mg_vcyc_smoother_iter_pre  ! V cycle smoother iterations pre-smoothing
  integer, protected       :: pub_mg_vcyc_smoother_iter_post ! V cycle smoother iterations post-smoothing
  integer, protected       :: pub_mg_max_iters_defco   ! maximum iterations for defect correction
  integer, protected       :: pub_mg_max_iters_vcycle  ! maximum iteration for V-cycle
  integer, protected       :: pub_mg_max_iters_newton ! maximum iterations for Newton method
  logical, protected       :: pub_mg_use_error_damping ! use error damping in defect corr loop
  integer, protected       :: pub_mg_max_iters_cg      ! maximum iterations for conjugate gradients
  logical, protected       :: pub_mg_pbe_use_fas       ! use Full Approximation Scheme for non-linear
                                                       ! PBE (default is Newton method)
  real(kind=DP), protected :: pub_mg_tol_mu_rel   ! relative tolerance for chemical potential
  real(kind=DP), protected :: pub_mg_tol_mu_abs   ! absolute tolerance for chemical potential

  logical, protected       :: pub_mg_continue_on_error ! do not abort on errors

  real(kind=DP), protected :: pub_mg_tol_cg_res_rel    ! relative tolerance in residual for conjugate gradients
  real(kind=DP), protected :: pub_mg_tol_cg_res_abs    ! absolute tolerance in residual for conjugate gradients

  logical, protected       :: pub_mg_use_cg            ! use conjugate gradients

  ! JCW: Has a keyword that was made obsolete in DL_MG v2.0 update been found?
  logical, private         :: obsolete_parameter_found_DLMGv2

  ! JCW: Updated parameters for controlling finite differences
  ! JCW: These only affect ONETEP's internal finite_differences module. DL_MG's
  ! JCW: finite differences should be controlled with DL_MG-specific keywords
  ! JCW: (see below).
  integer, protected       :: pub_finite_difference_order ! replaces pub_is_discretization order
                                                          ! for non-DL_MG usage of high-order FD
  ! JCW: Additional parameters for control of BCs
  ! JCW:  - allow the user to select BCs for multigrid, ion-ion and local pseudo
  ! JCW:    along each Cartesian direction using a short string. This string has
  ! JCW:    three characters representing the BC along X, Y and Z directions
  ! JCW:    (in this order) and the characters may be P (periodic), O (open) or Z
  ! JCW:    (zero, only allowed for multigrid). White space is ignored, e.g.
  ! JCW:    "P P P" (all periodic)
  ! JCW:    "P P O" (periodic along X and Y, i.e. a surface)
  ! JCW:    "O O P" (periodic along Z direction)
  ! JCW: (this is overidden by pub_openbc_hartree, which is equivalent to
  ! JCW: pub_multigrid_bc = "O O O")
  character(len=80), protected :: pub_multigrid_bc
  logical, protected           :: pub_multigrid_bc_is_periodic(3) ! <-- set in rundat_check_inputs
  logical, protected           :: pub_multigrid_bc_is_zero(3)     ! <-- set in rundat_check_inputs
  character(len=80), protected :: pub_ion_ion_bc
  logical, protected           :: pub_ion_ion_bc_is_periodic(3)   ! <-- set in rundat_check_inputs
  character(len=80), protected :: pub_pspot_bc
  logical, protected           :: pub_pspot_bc_is_periodic(3)     ! <-- set in rundat_check_inputs
  character(len=80), protected :: pub_smeared_ion_bc
  logical, protected           :: pub_smeared_ion_bc_is_periodic(3) ! <-- set in rundat_check_inputs
  ! JCW: NOTE: To avoid confusion over BCs, the pub_*_is_periodic and pub_*_is_zero
  ! JCW: logical arrays should be used to determine the current BCs in other
  ! JCW: modules. References to pub_open_localpseudo, pub_ii_energy_direct
  ! JCW: and other keywords should be replaced with checks based on these.
  ! JCW: DO NOT USE pub_*_bc (CHARACTER VARIABLES) TO CHECK FOR BCS OUTSIDE
  ! JCW: OF THE RUNDAT MODULE.
  ! ab: van der Waals BCs
  character(len=80), protected :: pub_vdw_bc
  logical, protected           :: pub_vdw_bc_is_periodic(3) ! <-- set in rundat_check_inputs
  ! ab: DFTB BCs
  character(len=80), protected :: pub_dftb_bc
  logical, protected           :: pub_dftb_bc_is_periodic(3) ! <-- set in rundat_check_inputs
  ! jd: HFx BCs
  character(len=80), protected :: pub_hfx_bc
  logical, protected           :: pub_hfx_bc_is_periodic(3) ! <-- set in rundat_check_inputs

  logical, protected           :: pub_external_bc_from_cube ! Read in boundary conditions from cube file

  ! vv: Not a keyword. Check whether mixed bcs are used for aborting MD
  ! calculations
  logical            :: pub_mix_bcs = .false.  ! vv: modified in rundat
  ! vv: Not a keyword
  logical            :: pub_isthermo = .false.      ! modified in rundat_blocks_mod

  ! jd: Not a keyword. Set automatically depending on tasks and pub_write_forces
  logical, protected :: pub_ions_will_move = .false.

  ! jd: Not a keyword. Set automatically depending on tasks and pub_write_forces
  logical, protected :: pub_forces_needed = .false.

  logical, protected :: pub_permit_unusual_ngwf_count

  ! lpl: DDEC parameters
  logical, protected :: pub_ddec
  logical, protected :: pub_ddec_write
  logical, protected :: pub_ddec_hirshfeld

  integer, protected :: pub_ddec_rad_npts
  real(kind=DP), protected :: pub_ddec_rad_rcut

  integer, protected :: pub_ddec_maxit
  integer, protected :: pub_ddec_core_maxit
  real(kind=DP), protected :: pub_ddec_conv_thresh
  real(kind=DP), protected :: pub_ddec_IH_frac
  real(kind=DP), protected :: pub_ddec_zero_thresh
  logical, protected :: pub_ddec_use_coredens

  logical, protected :: pub_ddec_core_correction
  logical, protected :: pub_ddec_reshape_dens
  integer, protected :: pub_ddec_core_corr_maxit

  integer, protected :: pub_ddec_ionic_range
  !integer, protected :: pub_ddec_atomsolve_maxit
  !integer, protected :: pub_ddec_rcomp_maxit
  !real(kind=DP), protected :: pub_ddec_target_radius
  !real(kind=DP), protected :: pub_ddec_rcomp_econv
  !logical, protected :: pub_ddec_write_rcomp

  character(len=80), protected :: pub_ddec_refdens_path
  logical, protected :: pub_ddec_refdens_init
  logical, protected :: pub_ddec_renormalize_refdens

  logical :: pub_ddec_c3_refdens
  character(len=80), protected :: pub_ddec_rad_shell_mode

  logical, protected :: pub_ddec_multipole
  integer, protected :: pub_ddec_moment
  integer, protected :: pub_ddec_moment_order

  logical, protected :: pub_ddec_interp_rad_dens
  logical, protected :: pub_ddec_avg_rad
  real(kind=DP), protected :: pub_ddec_min_shell_dens
  integer, protected :: pub_ddec_ref_shell_mode

  logical, protected :: pub_ddec_format_dens

  logical, protected :: pub_ddec_eff_decay_exp
  real(kind=DP), protected :: pub_ddec_eff_decay_rmin
  ! lpl: This is related to the still-dodgy subroutine 'ddec_rmse'
  !logical, protected :: pub_ddec_rmse

  ! lpl: Not inputs
  character(len=512), allocatable :: pub_ddec_rcomp(:) ! modified in ddec_mod
  ! lpl: This is related ro the still-dodgy subroutine 'ddec_rmse'
  !character(len=256), allocatable :: pub_ddec_rmse_vdW(:) ! modified in rundat_blocks_mod
  ! lpl: END DDEC parameters

  ! aeaa: DDEC anisotropy parameters
  logical, protected :: pub_ddec_aniso
  real(kind=DP), protected :: pub_ddec_aniso_max_dis
  real(kind=DP), protected :: pub_ddec_aniso_max_dis_halogen
  real(kind=DP), protected :: pub_ddec_aniso_error_thres
  real(kind=DP), protected :: pub_ddec_aniso_error_reduce

  ! lpl: NBO stuff (part of properties_calculate)
  logical, protected :: pub_write_nbo
  logical, protected :: pub_nbo_init_lclowdin
  logical, protected :: pub_nbo_write_lclowdin
  logical, protected :: pub_nbo_write_npacomp
  logical, protected :: pub_nbo_scale_dm
  logical, protected :: pub_nbo_write_dipole
  logical, protected :: pub_nbo_scale_spin

  logical, allocatable            :: pub_nbo_write_species(:)      ! not input  ! modified in rundat_blocks_mod
  character(len=256), allocatable            :: pub_nbo_ngwf_label(:) ! not input ! modified in rundat_blocks_mod

  logical, protected :: pub_nbo_pnao_analysis
  character(len=80), protected :: pub_nbo_aopnao_scheme

  logical, protected :: pub_plot_nbo
  character(len=8)            :: pub_nbo_plotorbtype ! modified in properties_mod
  integer, allocatable            :: pub_nbo_list_plotnbo(:) ! not input ! modified in rundat_blocks_mod
  ! lpl: NBO stuff (part of properties_calculate)

  ! kkbd: Image-comms related variables
  integer, protected :: pub_num_images
  character(len=80), protected :: pub_image_sizes

  ! kaw: OpenMP related variables
  integer, protected :: pub_threads_max
  integer, protected :: pub_threads_max_possible

  ! Number of "Loop parallel" threads
  integer, protected :: pub_threads_num_fftboxes

  ! Number of "Data parallel" threads
  integer, protected :: pub_threads_per_fftbox
  integer, protected :: pub_threads_per_cellfft

  ! "Data parallel" switches, activated if more than one thread used.
  logical, protected :: pub_threads_fftbox
  logical, protected :: pub_threads_cellfft

  ! smmd: Number of threads available for MKL
  integer, protected :: pub_threads_num_mkl

  ! ddor: Use root mean squared measure of kernel-Hamiltonian commutator and delta K
  logical, protected :: pub_rms_kernel_measure

  ! az, jd: DMA-related parameters
  character(len=80), private   :: dma_metric_string
  character(len=32), protected :: pub_dma_use_ri
  integer, protected           :: pub_dma_metric
  logical, protected           :: pub_dma_calculate
  logical, protected           :: pub_dma_output_potential
  logical, protected           :: pub_dma_output_potential_reference
  logical, protected           :: pub_dma_scale_charge
  logical, protected           :: pub_dma_bessel_averaging
  integer, protected           :: pub_dma_max_l
  integer, protected           :: pub_dma_max_q
  logical, protected           :: pub_dma_precise_gdma_output
  real(kind=DP), protected     :: pub_dma_target_num_val_elec
  real(kind=DP), protected     :: pub_dma_multipole_scaling
  real(kind=DP), protected     :: pub_dma_dipole_scaling
  real(kind=DP), protected     :: pub_dma_quadrupole_scaling

  ! jd: Whether the spherical wave expansion module is used
  logical, protected :: pub_use_swx

  ! jd: Whether the remote NGWF facility is used
  logical :: pub_use_remote_ngwfs ! not protected, set from e&f

  ! JCW: Additional exchange-correlation functional options
  real(kind=DP), protected :: pub_xc_min_tau ! minimum tau threshold (tautol)
  character(len=80), protected :: pub_xc_initial_functional ! alternative functional used
                             ! to do initial guess -- only affects meta-GGAs
  ! JCW: Whether kinetic energy density is required
  logical, save, public :: pub_xc_ke_density_required   ! <-- set in xc_init
  ! JCW: Whether Laplacian of density is required
  logical, save, public :: pub_xc_lapl_density_required ! <-- set in xc_init

  ! vv: Seed for generating velocities in md
  integer, save :: pub_rand_seed

  ! cks: HHF parameters
  integer, protected       :: pub_hhf_nstates
  real(kind=DP)            :: pub_hhf_factor  ! jd: set in e&f

  ! jcap: Parameters for embedding mean field theory
  character(len=80)        :: pub_active_xc_functional
  logical                  :: pub_emft
  integer                  :: pub_active_region
  logical, protected :: pub_use_activehfx ! set in xc_activehfxinit, through rundat_set_pub_use_activehfx
  logical, protected :: pub_use_activeswx
  logical                  :: pub_emft_follow ! rc2013: do EMFT after a regular optimisation
  logical                  :: pub_emft_lnv_only ! rc2013: do LNV only
  logical                  :: pub_block_orthogonalise
  logical                  :: pub_singlet_triplet_split ! rc2013: calculate singlet-triplet splitting
  logical                  :: pub_build_bo = .false. ! modified in E+F, ngwf_cg
  integer                  :: pub_emft_lnv_steps ! for EMFT LNV-only runs

  ! rjc: Electron localisation descriptors and Kinetic Energy Density output
  logical, protected, public :: pub_eld_calculate ! enable/disable eld
  character(len=80), protected, public :: pub_eld_function ! choose which eld function to use
  logical, protected, public :: pub_ke_density_calculate ! enable/disable ke density

  ! jd: Parameters for DFTB
  logical, protected            :: pub_dftb
  integer, protected            :: pub_dftb_method
  character(len=80), private    :: dftb_method_string
  character(len=512), protected :: pub_dftb_method_param_file
  character(len=512), protected :: pub_dftb_common_param_file
  real(kind=DP), protected      :: pub_dftb_coord_cutoff
  real(kind=DP), protected      :: pub_dftb_rep_cutoff
  real(kind=DP), protected      :: pub_dftb_srb_cutoff
  real(kind=DP), protected      :: pub_dftb_overlap_cutoff
  logical, protected            :: pub_dftb_overlap_analytical
  logical, protected            :: pub_dftb_cartesian_ngwfs
  logical, protected            :: pub_dftb_ewald_replicate_xtb
  real(kind=DP), protected      :: pub_dftb_ewald_parameter
  logical, protected            :: pub_dftb_calculate_force

contains

  subroutine get_rundat

    use comms, only: comms_bcast, pub_on_root, pub_root_proc_id, &
         pub_comms_group_size
    use constants, only: DP, stdout, ANGSTROM, NORMAL, BRIEF, VERBOSE, PROLIX, &
         MAXIMUM, CRLF, garbage_real, garbage_int, DFTB_GFN0, DFTB_GFN1, DFTB_GFN2, &
         DFTB_DEFAULT_COMMON_PARAM_FILE, DFTB_DEFAULT_METHOD_PARAM_FILES
    use esdf, only: esdf_double, esdf_string, esdf_integer, esdf_boolean,&
         esdf_block, block_data, esdf_physical, esdf_convfac, esdf_defined
    use utils, only: utils_alloc_check, utils_assert, utils_abort, &
                     utils_assert, utils_strip_char_from_string

    character(len=80) :: buf, tdbuf ! ddor: tdbuf for TDDFT
    character(len=80) :: output_detail
    character(len=80) :: paw_output_detail
    character(len=80) :: md_output_detail
    character(len=80) :: swx_output_detail
    character(len=80) :: hfx_output_detail
    character(len=80) :: geom_output_detail
    character(len=40) :: psinc_spacing_unit
    character(len=40) :: efield_origin_unit
    character(len=40) :: efield_unit
    character(len=40) :: tmp_unit
    character(len=2000) :: errmsg
    integer :: ierr
    integer :: nlines,iline,ipath
    character(len=500) :: default_pseudo_path

    logical :: bs_start, bsunfld_start
    real(kind=DP) :: dummy_real(3)

    character(len=*), parameter :: myself = 'get_rundat'

    ! -------------------------------------------------------------------------

    pub_debug_on_root = pub_debug .and. pub_on_root

    ! Root proc reads from input file...

    if (pub_on_root) then
       pub_check_stack_size    = esdf_boolean('check_stack_size',.true.)
       pub_all_tasks           = esdf_string('task', 'SINGLEPOINT')
       if (index(pub_all_tasks, 'UNIT_TEST') == 0) &
          pub_cutoff_energy           = esdf_physical('cutoff_energy', -20.0_DP, &
              'hartree') ! jd: The minus will be flipped if psinc_spacing is
                         !     not provided.

       pub_kernel_cutoff       = esdf_physical('kernel_cutoff', 1000.0_DP, &
            'bohr')
       pub_xc_functional       = esdf_string('xc_functional','LDA')
       pub_libxc_x_func_id     = esdf_integer('libxc_x_func_id', 0)
       pub_libxc_c_func_id     = esdf_integer('libxc_c_func_id', 0)
       pub_charge              = esdf_double('charge', 0.0_DP)
       ! kkbd: see internal_setup_spin in energy_and_force_mod for why we must initialize the spin like this.
       pub_real_spin           = esdf_double('spin', real(garbage_int, kind=dp))
       pub_spin                = nint(pub_real_spin)
       pub_spin_polarised      = esdf_boolean('spin_polarized',.false.)
       pub_confined_ngwfs          = esdf_boolean('confined_ngwfs',.false.)
       pub_maxit_ngwf_cg_confined  = esdf_integer('maxit_ngwf_cg_confined',5)
       pub_confined_ngwfs_barrier  = esdf_double('confined_ngwfs_barrier',2000.0_DP)
       pub_spin_polarised      = esdf_boolean('spin_polarised', &
            pub_spin_polarised)
       ! rab207: append "Ha/bohr/e" to constant_efield as default units
       buf                     = esdf_string('constant_efield','0.0 0.0 0.0')
       buf = trim(adjustl(buf))//' "Ha/bohr/e"'
       read(buf,*) pub_constant_efield(:), efield_unit
       pub_constant_efield(:) = pub_constant_efield(:) * &
            esdf_convfac(efield_unit,'Ha/bohr/e')
       ! agreco: by default set the origin of the sawtooth potential at (0,0,0)
       buf                     = esdf_string('efield_origin','0.0 0.0 0.0')
       ! append "bohr" to efield_origin as default units
       buf = trim(adjustl(buf))//' bohr'
       read(buf,*) pub_efield_origin(:),efield_origin_unit
       pub_efield_origin(:) = pub_efield_origin(:) * &
           esdf_convfac(efield_origin_unit,'bohr')
       buf                     = esdf_string('fftbox_pref','0 0 0')
       read(buf,*) pub_fftbox_pref(:)
       buf                     = esdf_string('augbox_pref','0 0 0')
       read(buf,*) pub_augbox_pref(:)
       buf                     = esdf_string('ppd_npoints', '0 0 1')
       read(buf,*) pub_ppd_npoints(:)
       ! agreco: by default do not use extended NGWFs
       buf                     = esdf_string('extend_ngwf', 'F F F')
       read(buf,*) pub_extend_ngwf(:)
       ! agreco: by default NGWFs are not initialised to random values
       pub_full_rand_ngwf          = esdf_boolean('full_rand_ngwf',.false.)
       ! agreco: by default set seed to negative number
       pub_rand_seed_ngwf      = esdf_integer('rand_seed_ngwf',-1)
       ! agreco: by default not using normal distribution
       pub_rand_sigma          = esdf_double('rand_normal_sigma',0.0_DP)
       ! agreco: by default, convolution is turned off
       pub_rand_conv           = esdf_boolean('convolute_rand',.false.)
       ! agreco: use erfc as default convolution function
       pub_conv_func           = esdf_string('convolute_func','erfc')
       ! agreco: specify width of convolution region
       buf                     = esdf_string('conv_region_width','-1.0')
       ! append "bohr" to efield_origin as default units
       buf = trim(adjustl(buf))//' bohr'
       read(buf,*) pub_conv_region_width, conv_region_width_unit
       pub_conv_region_width = pub_conv_region_width * &
                esdf_convfac(conv_region_width_unit,'bohr')
       ! agreco: if true, write NGWFs for plotting at every step
       pub_ngwf_plot_every_it  = esdf_boolean('write_ngwf_plot_every_it',.false.)
       ! agreco: if true, use complex valued NGWFs
       pub_cmplx_ngwfs         = esdf_boolean('use_cmplx_ngwfs',.false.)
       ! agreco: if true, read real NGWFs from file into complex NGWFs
       pub_read_real_ngwfs     = esdf_boolean('read_real_ngwfs',.false.)
       ! agreco: by default do not store overlap matrix in a file
       pub_show_overlap        = esdf_boolean('show_overlap',.false.)
       ! agreco: single k-point for testing
       buf                     = esdf_string('single_kpoint','0.0 0.0 0.0')
       read(buf,*) pub_single_kpt(:)
       ! agreco: k-point method for BZ sampling (use KP by default)
       pub_kpoint_method       = esdf_string('kpoint_method','KP')
       ! agreco: threshold to check imaginary part of density with complex NGWFs
       pub_imag_thr            = esdf_double('imag_thr', 1.0e-6_dp)
       ! agreco: multiply all initial complex NGWFs by exp(i * phase)
       pub_ngwfs_phase         = esdf_double('mult_ngwf_by_phase', 0.0_DP)
       ! agreco: multiply all initial complex NGWFs by a random phase in [0,2PI]
       pub_ngwfs_rand_phase    = esdf_boolean('mult_ngwf_by_random_phase',.false.)
       ! agreco: check hermitian character of relevant complex matrices
       pub_check_hermitian     = esdf_boolean('check_hermitian_mats',.false.)
       ! agreco: check density is real when using complex NGWFs
       pub_check_density       = esdf_boolean('check_density',.false.)
       ! agreco: by default use 1x1x1 k-point grid
       buf                     = esdf_string('kpoint_grid_size','1 1 1')
       read(buf,*) pub_kp_grid_size(:)
       ! agreco: by default do not shift the grid
       buf                     = esdf_string('kpoint_grid_shift','0 0 0')
       read(buf,*) pub_kp_grid_shift(:)


       pub_contracoham_radmult = esdf_double('contracoham_radmult', 1.0_dp)

       ! mjsp: Force cutoff=523.78eV for unit testing
       if (index(pub_all_tasks, 'UNIT_TEST') > 0) then
          if (pub_on_root) write(stdout,'(//a)') 'WARNING: Cutoff energy overridden to &
             &523.78 eV for unit testing.'
          pub_psinc_spacing =  [0.628205128205_dp,0.628205128205_dp,0.628205128205_dp]
          psinc_spacing_unit = 'bohr'
          pub_cutoff_energy = -(((6.0_DP*(PI**2))**(2.0_DP/3.0_DP)) / 2.0_DP) / &
               (maxval(pub_psinc_spacing)**2)
       else
          buf                     = esdf_string('psinc_spacing','0.0 0.0 0.0')
          ! ndmh: append "bohr" to pub_psinc_spacing as default units
          buf = trim(adjustl(buf))//' bohr'
          read(buf,*) pub_psinc_spacing(:),psinc_spacing_unit
          pub_psinc_spacing(:) = pub_psinc_spacing(:) * &
               esdf_convfac(psinc_spacing_unit,'bohr')
       end if

       pub_devel_code          = esdf_string('devel_code','')
       ! rab207: third term in esdf_string prevents case conversion
       pub_input_xyz_file      = esdf_string('positions_xyz_file','',.true.)
       pub_input_xyz_file_intermediate = esdf_string('positions_intermediate_xyz_file','',.true.)
       pub_input_xyz_file_product      = esdf_string('positions_product_xyz_file','',.true.)
       pub_nnho                = esdf_boolean('nnho',.false.)
       pub_ngwf_halo           = esdf_physical('ngwf_halo', -1.0_DP, 'bohr')
#ifdef OLD_DEFAULTS
       pub_nonsc_forces        = esdf_boolean('nonsc_forces',.false.)
#else
       pub_nonsc_forces        = esdf_boolean('nonsc_forces',.true.)
#endif
       pub_external_pressure   = esdf_physical('external_pressure',0.0_DP, &
            'ha/bohr**3')
       pub_smoothing_factor    = esdf_double('smoothing_factor',5.0_DP)
       pub_isosurface_cutoff   = esdf_double('isosurface_cutoff',0.0005_DP)
       pub_maxit_palser_mano   = esdf_integer('maxit_palser_mano', 200)
       pub_maxit_kernel_occ_check = esdf_integer('maxit_kernel_occ_check',0)
       pub_maxit_pen           = esdf_integer('maxit_pen', 0)
       pub_pen_param               = esdf_double('pen_param', 4.0_DP)
#ifdef OLD_DEFAULTS
       minit_lnv               = esdf_integer('minit_lnv',3)
       maxit_lnv               = esdf_integer('maxit_lnv',8)
#else
       minit_lnv               = esdf_integer('minit_lnv',5)
       maxit_lnv               = esdf_integer('maxit_lnv',5)
#endif
       lnv_threshold_orig      = esdf_double('lnv_threshold_orig',1.0e-9_dp)
       pub_lnv_cg_type         = esdf_string('lnv_cg_type','LNV_FLETCHER')
       pub_lnv_cg_max_step     = esdf_double('lnv_cg_max_step', 3.00_DP)
       pub_exact_lnv           = esdf_boolean('exact_lnv',.true.)
       pub_old_lnv             = esdf_boolean('old_lnv',.false.)
#ifdef OLD_DEFAULTS
       pub_lnv_check_trial_steps = esdf_boolean('lnv_check_trial_steps', &
            .false.)
       pub_kerfix              = esdf_integer('kerfix', 1)
       pub_maxit_kernel_fix    = esdf_integer('maxit_kernel_fix', 3)
       pub_initial_dens_realspace = esdf_boolean('initial_dens_realspace', &
            .false.)
#else
       pub_lnv_check_trial_steps = esdf_boolean('lnv_check_trial_steps', &
            .true.)
       pub_kerfix              = esdf_integer('kerfix', 1)
       pub_maxit_kernel_fix    = esdf_integer('maxit_kernel_fix', 3)
       pub_initial_dens_realspace = esdf_boolean('initial_dens_realspace', &
            .true.)
#endif
       pub_kernel_track_mid_occ= esdf_boolean('kernel_track_mid_occ',.false.)
       pub_kernel_check_all    = esdf_boolean('kernel_check_all',.false.)
       pub_kernel_diis_size    = esdf_integer('kernel_diis_size', 10)
       pub_kernel_diis_maxit   = esdf_integer('kernel_diis_maxit', 25)
       pub_kernel_diis_threshold   = esdf_double('kernel_diis_threshold', 1.0e-9_DP)
       pub_kernel_diis_linear_iter   = esdf_integer('kernel_diis_linear_iter', 5)
       pub_kernel_diis_coeff   = esdf_double('kernel_diis_coeff', 0.1_DP)
       buf = esdf_string('kernel_diis_conv_criteria', '1000')
       pub_kernel_diis_conv_criteria = buf(1:4)
       pub_kernel_diis_lshift  = esdf_physical('kernel_diis_lshift',1.0_DP,'hartree')
       pub_kernel_diis_ls_iter = esdf_integer('kernel_diis_ls_iter', 0)
       pub_kernel_diis_scheme  = esdf_string('kernel_diis_scheme','NONE')

       pub_foe                     = esdf_boolean('foe', .false.)
       pub_dense_foe               = esdf_boolean('dense_foe', .false.)
       pub_H2denskern_sparsity     = esdf_boolean('H2denskern_sparsity', .false.)
       pub_foe_mu_tol              = esdf_physical('foe_mu_tol',1.0e-7_dp,'hartree')
       pub_foe_test_sparsity       = esdf_boolean('foe_test_sparsity', .false.)
       pub_foe_cheby_thres         = esdf_double('foe_cheby_thres',1.0e-9_dp)
       pub_foe_avoid_inversions    = esdf_boolean('foe_avoid_inversions', .false.)
       pub_foe_check_entropy       = esdf_boolean('foe_check_entropy', .true.)


       ! jd: DFTB
       pub_dftb  = esdf_boolean('dftb',.false.)
       dftb_method_string = esdf_string('dftb_method','GFN0')

       select case(trim(dftb_method_string))
       case ('GFN0')
          pub_dftb_method = DFTB_GFN0
       case ('GFN1')
          pub_dftb_method = DFTB_GFN1
       case ('GFN2')
          pub_dftb_method = DFTB_GFN2
       case default
          call utils_abort("Only GFN0, GFN1 and GFN2 are supported dftb_method.")
       end select

       ! ab: use EDFT procedures in DFTB
       if (pub_dftb) then
          pub_edft                 = .true.
       else
          pub_edft                 = esdf_boolean('edft', .false.)
       end if
       pub_edft_maxit              = esdf_integer('edft_maxit', 10)
       pub_edft_max_step           = esdf_double('edft_max_step', 1.0_DP)
       pub_edft_smearing_width     = esdf_physical('edft_smearing_width', &
           3.166811429e-3_DP, 'hartree')
       pub_edft_free_energy_thres  = esdf_physical('edft_free_energy_thres', 1.0e-6_DP,'hartree')
       pub_edft_energy_thres       = esdf_physical('edft_energy_thres', 1.0e-4_DP,'hartree')
       pub_edft_entropy_thres      = esdf_physical('edft_entropy_thres', 1.0e-4_DP,'hartree')
       pub_edft_rms_gradient_thres = esdf_double('edft_rms_gradient_thres', 0.0_DP)
       pub_edft_commutator_thres   = esdf_physical('edft_commutator_thres', 1.0e-5_DP,'hartree')
       pub_edft_fermi_thres        = esdf_physical('edft_fermi_thres', 1.0e-3_DP,'hartree')
       pub_edft_trial_step         = esdf_double('edft_trial_step', -1.0_DP)
       pub_edft_nelec_thres        = esdf_double('edft_nelec_thres',1.0e-6_DP)
       pub_edft_write_occ          = esdf_boolean('edft_write_occ', .false.)
       pub_edft_extra_bands        = esdf_integer('edft_extra_bands', -1)
       pub_edft_init_maxit         = esdf_integer('edft_init_maxit', -1)
       pub_edft_spin_fix           = esdf_integer('edft_spin_fix', -1)
       pub_edft_spin_fix_orig      = pub_edft_spin_fix
       pub_edft_update_scheme      = esdf_string('edft_update_scheme','DAMP_FIXPOINT')
       pub_edft_ham_diis_size      = esdf_integer('edft_ham_diis_size', 10)
       pub_eigensolver_orfac       = esdf_double('eigensolver_orfac', 1.0e-4_DP)
       pub_eigensolver_abstol      = esdf_double('eigensolver_abstol', 1.0e-9_DP)
       pub_write_hamiltonian       = esdf_boolean('write_hamiltonian', .false.)
       pub_read_hamiltonian        = esdf_boolean('read_hamiltonian', .false.)
       pub_write_overlap           = esdf_boolean('write_overlap', .false.)
       pub_ngwf_cg_rotate          = esdf_boolean('ngwf_cg_rotate', .false.)
       pub_edft_round_evals        = esdf_integer('edft_round_evals', -1)

       ! kkbd: If we're fixing the spin, we want to make sure spin polarisation is enabled.
       if (pub_edft_spin_fix .ge. 0) then
          pub_spin_polarised = .true.
       end if

       ! ab: grand canonical ensemble dft
       pub_edft_grand_canonical    = esdf_boolean('edft_grand_canonical', .false.)
       pub_edft_reference_potential= esdf_physical('edft_reference_potential', mu_ref_SHE,'ha')
       pub_edft_electrode_potential= esdf_physical('edft_electrode_potential', 0.0_DP, 'ha/e')
       if (pub_edft_grand_canonical) then
          pub_edft = .true.
       end if

!       pub_pbc_smeared_ion_rep     = esdf_boolean('pbc_smeared_ion_rep',.false.) !ja531-> these are devel codes for now
!       pub_chemical_softness       = esdf_boolean('chemical_softness',.false.)   !ja531-> these are devel codes for now

       pub_energy_components_interval = esdf_integer(&
            'energy_components_interval',5)
       ! ab: no ngwf optimization in dftb-gfn0, do not calculate non-scf forces
       if (pub_dftb) then
          if (pub_dftb_method == DFTB_GFN0) then
             pub_maxit_ngwf_cg    = 0
             pub_nonsc_forces     = .false.
          else
             call utils_abort("Unsupported DFTB method. Only GFN0 is supported.")
          end if
       else
          pub_maxit_ngwf_cg    = esdf_integer('maxit_ngwf_cg',60)
       end if
       pub_freeze_switch_steps = esdf_integer('freeze_switch_steps',-1)
       pub_do_fandt            = esdf_boolean('do_fandt', .false.)
       pub_freeze_envir_ngwfs  = esdf_boolean('freeze_envir_ngwfs', .false.)
       pub_embed_debug         = esdf_boolean('embed_debug', .false.)
       pub_project_embed       = esdf_boolean('project_embed', .false.)
       pub_ngwf_cg_type        = esdf_string('ngwf_cg_type','NGWF_FLETCHER')
       pub_ngwf_cg_max_step    = esdf_double('ngwf_cg_max_step',-8.00_DP)
#ifdef OLD_DEFAULTS
       pub_elec_cg_max         = esdf_integer('elec_cg_max', 5)
#else
       pub_elec_cg_max         = esdf_integer('elec_cg_max', 3)
#endif
       pub_precond_scheme      = esdf_string('precond_scheme','TETER')
       pub_k_zero              = esdf_physical('k_zero', 3.0_dp,'1/bohr')
       pub_r_precond           = esdf_physical('r_precond', 2.0_dp,'bohr')
       pub_precond_recip       = esdf_boolean('precond_recip',.true.)
       pub_precond_real        = esdf_boolean('precond_real',.false.)
       pub_smooth_scheme       = esdf_string('smooth_scheme','NONE')
       pub_r_smooth            = esdf_physical('r_smooth',1.5_dp,'bohr')
       pub_k_smooth            = esdf_physical('k_smooth',5.0_dp,'1/bohr')
       pub_occ_mix                 = esdf_double('occ_mix', 0.25_DP)
       pub_kernel_update       = esdf_boolean('kernel_update', .false.)
       pub_kernel_christoffel_update = &
            &esdf_boolean('kernel_christoffel_update',.false.)
       pub_maxit_hotelling     = esdf_integer('maxit_hotelling', 50)
       pub_max_resid_hotelling     = esdf_double('max_resid_hotelling', 1.0e-12_DP)
       ! ddor: Use root mean squared measure of commutator and delta K
       pub_rms_kernel_measure  = esdf_boolean('rms_kernel_measure',.false.)

       ngwf_threshold_orig     = esdf_double('ngwf_threshold_orig',2.0e-6_dp)
       pub_elec_energy_tol     = esdf_physical('elec_energy_tol',-0.001_DP, &
            'hartree')
       pub_elec_force_tol      = esdf_physical('elec_force_tol',-0.001_DP, &
            'ha/bohr')
       pub_ngwf_max_grad       = esdf_double('ngwf_max_grad',-2.0e-5_dp)
       pub_delta_e_conv            = esdf_boolean('delta_e_conv', .true.)
       pub_kernel_force_conv      = esdf_boolean('kernel_force_conv', .false.)

       pub_fftbox_batch_size   = esdf_integer('fftbox_batch_size', 16)
       pub_dense_threshold     = esdf_double('dense_threshold',0.40_DP)
       pub_comms_group_size    = esdf_integer('comms_group_size',-1)
       pub_ovlp_for_nonlocal   = esdf_boolean('ovlp_for_nonlocal',.false.)
       pub_use_space_filling_curve = esdf_boolean('use_space_filling_curve', .true.)
       pub_coreham_denskern_guess  = esdf_boolean('coreham_denskern_guess', .true.)
       pub_check_atoms         = esdf_boolean('check_atoms', .true.)
       pub_locpot_scheme       = esdf_string('locpot_scheme','FULL')
       pub_smooth_projectors   = esdf_double('smooth_projectors', -0.4_DP)
       pub_smooth_loc_pspot    = esdf_double('smooth_loc_pspot', -0.4_DP)
       pub_odd_psinc_grid      = esdf_boolean('odd_psinc_grid',.false.)
       pub_even_psinc_grid      = esdf_boolean('even_psinc_grid',.false.)
       pub_realspace_projectors= esdf_boolean('realspace_projectors',.false.)
       pub_projectors_precalculate= esdf_boolean('projectors_precalculate',.true.)
       pub_parallel_scheme     = esdf_string('parallel_scheme', 'NONE')

       pub_num_images  = esdf_integer('num_images',1)
       pub_image_sizes = esdf_string('image_sizes','DEFAULT')

       pub_num_images  = esdf_integer('num_images',1)
       pub_image_sizes = esdf_string('image_sizes','DEFAULT')

       output_detail           = esdf_string('output_detail','NORMAL')
       paw_output_detail       = esdf_string('paw_output_detail','DEFAULT')
       md_output_detail        = esdf_string('md_output_detail','DEFAULT')
       swx_output_detail       = esdf_string('swx_output_detail','DEFAULT')
       pub_timings_level       = esdf_integer('timings_level', 1)
       hfx_output_detail       = esdf_string('hfx_output_detail','DEFAULT')
       ! lk: if geom_output_detail is not specified in the input
       !     file then it is equal to 'DEFAULT'.
       geom_output_detail      = esdf_string('geom_output_detail','DEFAULT')
       pub_timings_order       = esdf_string('timings_order', 'TIME')
       pub_max_runtime         = esdf_double('run_time', -1.0_DP)
       pub_write_params        = esdf_boolean('write_params', .false.)
       pub_esdf_dump           = esdf_boolean('esdf_dump', .false.)
       pub_write_forces        = esdf_boolean('write_forces',.false.)
       pub_write_positions     = esdf_boolean('write_positions',.true.)
       pub_write_velocities    = esdf_boolean('write_velocities',.false.)
       pub_write_xyz           = esdf_boolean('write_xyz', .false.)
       pub_print_qc            = esdf_boolean('print_qc',.false.)
#ifdef ACCELRYS
       pub_cube_format         = esdf_boolean('cube_format', .false.)
#else
       pub_cube_format         = esdf_boolean('cube_format', .true.)
#endif
       pub_dx_format           = esdf_boolean('dx_format', .false.)
       pub_dx_coarse           = esdf_boolean('dx_format_coarse', .false.)
       pub_dx_sig_digits       = esdf_integer('dx_format_digits', 7)
#ifdef ACCELRYS
       pub_grd_format          = esdf_boolean('grd_format', .true.)
#else
       pub_grd_format          = esdf_boolean('grd_format', .false.)
#endif
       pub_write_density_plot  = esdf_boolean('write_density_plot', .true.)
       pub_write_polarisation_plot  = esdf_boolean('write_polarisation_plot', .false.)
       pub_write_ngwf_plot     = esdf_boolean('write_ngwf_plot',.false.)
       pub_write_ngwf_grad_plot = esdf_boolean('write_ngwf_grad_plot',.false.)
       pub_write_ngwf_radial    = esdf_integer('write_ngwf_radial',0)
       pub_write_ngwf_grad_radial = esdf_integer('write_ngwf_grad_radial',0)
       pub_write_radial_step    = esdf_physical('write_radial_step',0.005_dp,'bohr')
       pub_write_radial_smear   = esdf_physical('write_radial_smear',0.01_dp,'bohr')
       pub_write_initial_radial_ngwfs = esdf_boolean('write_initial_radial_ngwfs',.false.)

       pub_hub_proj_read_only          = esdf_boolean('hubbard_proj_read_only',.false.)
       pub_dmft_fully_sc_h             = esdf_boolean('dmft_fully_sc_h',.false.)
       pub_dmft_fully_sc               = esdf_boolean('dmft_fully_sc',.false.)
       pub_dmft_kernel                 = esdf_integer('dmft_kernel',0)
       pub_dmft_nbo                    = esdf_boolean('dmft_nbo',.false.)
       pub_dmft_nkpoints               = esdf_integer('dmft_nkpoints',1)
       pub_dmft_optics                 = esdf_boolean('dmft_optics',.false.)
       pub_dmft_order_proj             = esdf_double('dmft_order_proj',0.00_dp)
       pub_dmft_plot_real_space        = esdf_boolean('dmft_plot_real_space',.false.)
       pub_dmft_points                 = esdf_integer('dmft_points',0)
       pub_dmft_purify_sc              = esdf_boolean('dmft_purify_sc',.false.)
       pub_dmft_sc                     = esdf_boolean('dmft_sc',.false.)
       pub_dmft_spoil_kernel           = esdf_boolean('dmft_spoil_kernel',.false.)
       pub_dmft_switch_off_proj_order  = esdf_boolean('dmft_switch_off_proj_order',.false.)

       ! New dmft words
       pub_dmft_complex_freq           = esdf_boolean('dmft_complex_freq', .true.)
       pub_dmft_cutoff_small           = esdf_physical('dmft_cutoff_small',0.0_dp,'hartree')
       pub_dmft_dos_max                = esdf_physical('dmft_dos_max', 10.0_dp,'hartree')
       pub_dmft_dos_min                = esdf_physical('dmft_dos_min',-10.0_dp,'hartree')
       pub_dmft_emax                   = esdf_physical('dmft_emax', 1.0_dp,'hartree')
       pub_dmft_emin                   = esdf_physical('dmft_emin',-1.0_dp,'hartree')
       pub_dmft_kernel_mix             = esdf_double('dmft_kernel_mix',0.1d0)
       pub_dmft_kpoints_sym            = esdf_boolean('dmft_kpoints_sym',.false.)
       pub_dmft_ks_shift               = esdf_boolean('dmft_ks_shift',.true.)
       pub_dmft_mu_diff_max            = esdf_double('dmft_mu_diff_max',0.00d0)
       pub_dmft_mu_order               = esdf_integer('dmft_mu_order',2)
       pub_dmft_nmu_loop               = esdf_integer('dmft_nmu_loop',1)
       pub_dmft_nval                   = esdf_integer('dmft_nval',40)
       pub_dmft_optics_i1              = esdf_integer('dmft_optics_i1',1)
       pub_dmft_optics_i2              = esdf_integer('dmft_optics_i2',1)
       pub_dmft_optics_window          = esdf_physical('dmft_optics_window',0.10_dp,'hartree')
       pub_dmft_paramagnetic           = esdf_boolean('dmft_paramagnetic',.false.)
       pub_dmft_rotate_green           = esdf_boolean('dmft_rotate_green',.false.)
       pub_dmft_scaling_cutoff         = esdf_double('dmft_scaling_cutoff',1.d-8)
       pub_dmft_scaling_meth           = esdf_integer('dmft_scaling_meth',4)
       pub_dmft_scaling_nmpi           = esdf_integer('dmft_scaling_nmpi',1)
       pub_dmft_scaling_tail           = esdf_double('dmft_scaling_tail',10.d0)
       pub_dmft_skip_energy            = esdf_boolean('dmft_skip_energy',.false.)
       pub_dmft_smear                  = esdf_physical('dmft_smear',0.00018_dp,'hartree')
       pub_dmft_smear_eta              = esdf_physical('dmft_smear_eta',0.01_dp,'hartree')
       pub_dmft_smear_shift            = esdf_physical('dmft_smear_shift',0.000_dp,'hartree')
       pub_dmft_smear_T                = esdf_physical('dmft_smear_T',0.008_dp,'hartree')
       pub_dmft_smear_w                = esdf_physical('dmft_smear_w',0.035_dp,'hartree')
       pub_dmft_temp                   = esdf_physical('dmft_temp',-0.01_dp,'hartree')
       pub_dmft_win                    = esdf_double('dmft_win',0.3d0)
       pub_dmft_write                  = esdf_boolean('dmft_write',.true.)
       pub_dmft_read                   = esdf_boolean('dmft_read',.true.)

       ! DMFT keywords currently only included at the devel_code level
       ! pub_dmft_chem_shift             = esdf_physical('dmft_chem_shift',0.0_dp,'hartree')
       ! pub_dmft_cutoff_tail            = esdf_physical('dmft_cutoff_tail',10.0_dp,'hartree')
       ! pub_dmft_doping                 = esdf_double('dmft_doping',0.00d0)
       ! pub_dmft_embed_iter             = esdf_integer('dmft_embed_iter',10)
       ! pub_dmft_embed_mix              = esdf_double('dmft_embed_mix',0.7d0)
       ! pub_dmft_free_green_frequ       = esdf_physical('dmft_free_green_frequ',200.0_dp,'hartree')
       ! pub_dmft_gpu_num                = esdf_integer('dmft_gpu_num',0)
       ! pub_dmft_impose_chem_spin       = esdf_boolean('dmft_impose_chem_spin',.false.)
       ! pub_dmft_impose_same_coeffs     = esdf_boolean('dmft_impose_same_coeffs',.false.)
       ! pub_dmft_in_bohr                = esdf_boolean('dmft_in_bohr',.false.)
       ! pub_dmft_integrate_green        = esdf_boolean('dmft_integrate_green',.false.)
       ! pub_dmft_invert_overlap         = esdf_boolean('dmft_invert_overlap',.false.)
       ! pub_dmft_kpoints_kernel_gamma   = esdf_boolean('dmft_kpoints_kernel_gamma',.false.)
       ! pub_dmft_lin_scaling            = esdf_boolean('dmft_lin_scaling',.false.)
       ! pub_dmft_local_scratch          = esdf_boolean('dmft_local_scratch',.false.)
       ! pub_dmft_norm_proj              = esdf_integer('dmft_norm_proj',0)
       ! pub_dmft_optics_x1              = esdf_double('dmft_optics_x1',0.00_dp)
       ! pub_dmft_optics_y1              = esdf_double('dmft_optics_y1',0.00_dp)
       ! pub_dmft_optics_z1              = esdf_double('dmft_optics_z1',0.00_dp)
       ! pub_dmft_plot_all_proj          = esdf_boolean('dmft_plot_all_proj',.false.)
       ! pub_dmft_plot_real_space_sigma  = esdf_boolean('dmft_plot_real_space_sigma',.false.)
       ! pub_dmft_scaling_cutoff_h       = esdf_double('dmft_scaling_cutoff_h',1.d-5)
       ! pub_dmft_scaling_iter           = esdf_integer('dmft_scaling_iter',2000)
       ! pub_dmft_scaling_maxspace       = esdf_integer('dmft_scaling_maxspace',20)
       ! pub_dmft_scaling_tol            = esdf_double('dmft_scaling_tol',1.d-8)
       ! pub_dmft_split                  = esdf_boolean('dmft_split',.false.)
       ! pub_dmft_splitk                 = esdf_boolean('dmft_splitk',.false.)

       pub_read_denskern           = esdf_boolean('read_denskern', .false.)
       if (pub_dftb) then
          pub_write_denskern       = esdf_boolean('write_denskern', .false.)
       else
          pub_write_denskern       = esdf_boolean('write_denskern', .true.)
       end if
       pub_read_tightbox_ngwfs     = esdf_boolean('read_tightbox_ngwfs', .false.)
       pub_read_sub_denskern       = esdf_boolean('read_sub_denskern', .false.)
       if (pub_maxit_ngwf_cg==0) then
          pub_write_tightbox_ngwfs = esdf_boolean('write_tightbox_ngwfs', .false.)
       else
          pub_write_tightbox_ngwfs = esdf_boolean('write_tightbox_ngwfs', .true.)
       end if
       pub_write_converged_dk_ngwfs = esdf_boolean('write_converged_dk_ngwfs', .false.)
       pub_read_sw_ngwfs       = esdf_boolean('read_sw_ngwfs', .false.)
       pub_write_sw_ngwfs      = esdf_boolean('write_sw_ngwfs', .false.)
       pub_write_max_l         = esdf_integer('write_max_l', 3)
       pub_read_max_l          = esdf_integer('read_max_l', 3)
       pub_extra_n_sw          = esdf_integer('extra_n_sw', 0)

       pub_paw                 = esdf_boolean('paw', .false.)
       pub_fine_grid_scale     = esdf_double('fine_grid_scale',2.0_DP)
       pub_dbl_grid_scale      = esdf_double('dbl_grid_scale',2.0_DP)
       pub_aug_funcs_recip     = esdf_boolean('aug_funcs_recip',.true.)

#ifdef ACCELRYS
       ! Check if there is a defined pseudopotentials directory in the
       ! environment variable PSPOT_DIR
       call getenv('PSPOT_DIR', pub_pseudo_path)
#else
       default_pseudo_path = ''
       pub_pseudo_path             = esdf_string('pseudo_path',default_pseudo_path, &
            .true.)
#endif

       ! ab: use D2 vdw correction in DFTB
       if (pub_dftb) then
          pub_dispersion       = esdf_string('dispersion','4')
       else
          pub_dispersion       = esdf_string('dispersion','0')
       end if
       pub_vdw_dcoeff          = esdf_double('vdw_dcoeff',-1.0_DP)
       pub_vdw_radial_cutoff   = esdf_physical('vdw_radial_cutoff', 100.0_DP, 'bohr')

       pub_hub_max_iter        = esdf_integer('hubbard_max_iter', 10)
       pub_hub_energy_tol      = esdf_physical('hubbard_energy_tol', 1.0e-8_dp,&
            'hartree')
       pub_hub_conv_win        = esdf_integer('hubbard_conv_win', 2)
       pub_hub_proj_mixing     = esdf_double('hubbard_proj_mixing',0.0_DP)
       pub_hub_functional      = esdf_integer('hubbard_functional',1)
       pub_hub_tensor_corr     = esdf_integer('hubbard_tensor_corr',1)
       pub_hub_ngwf_spin_thr   = esdf_double('hubbard_ngwf_spin_threshold', &
            2.0e-5_dp)
       pub_hub_on_the_fly      = esdf_boolean('hubbardscf_on_the_fly',.false.)
       pub_hub_read_projectors = esdf_boolean('hubbard_read_projectors',.false.)
       pub_hubbard_compute_u_or_j  = esdf_boolean('hubbard_compute_u_or_j',.false.)
       !gibo: for U-(J-)-calculate runs, maintain DFT+U spin-splitting
       if (pub_hubbard_compute_u_or_j) pub_hub_ngwf_spin_thr = 1.E-50_DP
       pub_hubbard_j_minority_term  = esdf_boolean('hubbard_j_minority_term',.false.)

       !gibo: cDFT stuff===== START
       pub_cdft_atom_charge       = esdf_boolean('cdft_atom_charge',.false.)
       pub_cdft_atom_spin         = esdf_boolean('cdft_atom_spin',.false.)
       pub_cdft_group_charge_diff = esdf_boolean('cdft_group_charge_diff', &
            .false.)
       pub_cdft_group_spin_diff   = esdf_boolean('cdft_group_spin_diff', &
            .false.)
       pub_cdft_group_charge_acceptor=esdf_boolean('cdft_group_charge_acceptor',&
            .false.)
       pub_cdft_group_charge_donor   =esdf_boolean('cdft_group_charge_donor',   &
            .false.)
       pub_cdft_group_spin_acceptor  =esdf_boolean('cdft_group_spin_acceptor',  &
            .false.)
       pub_cdft_group_spin_donor     =esdf_boolean('cdft_group_spin_donor',     &
            .false.)
       pub_cdft_hubbard        = esdf_boolean('cdft_hubbard',.false.)
       pub_ci_cdft             = esdf_boolean('ci_cdft',.false.)
       pub_cdft_print_all_occ  = esdf_boolean('cdft_print_all_occ',.false.)
       pub_cdft_read_projectors= esdf_boolean('cdft_read_projectors',.false.)
       pub_cdft_group_charge_up_only =   &
                               esdf_boolean('cdft_group_charge_up_only',.false.)
       pub_cdft_group_charge_down_only = &
                               esdf_boolean('cdft_group_charge_down_only',.false.)
       pub_cdft_multi_proj =   esdf_boolean('cdft_multi_proj',.false.)
       pub_ci_cdft_num_conf    = esdf_integer('ci_cdft_num_conf',0)

       ! count number of active cdft_modes
       ! [1: constrained-charge OR constrained-spin only]
       pub_cdft_modes = 1
       !  [2: constrained-charge AND constrained-spin only]
       ! new with group_charge_up/down_only
       if (pub_cdft_group_charge_up_only .OR. pub_cdft_group_charge_down_only) then
         pub_cdft_modes = 1
       elseif ( ((pub_cdft_group_charge_acceptor.OR.pub_cdft_group_charge_donor).AND.&
           (pub_cdft_group_spin_acceptor .OR. pub_cdft_group_spin_donor)) .OR.   &
            (pub_cdft_group_charge_diff.AND.pub_cdft_group_spin_diff ) ) then
          pub_cdft_modes = 2
       endif
       ! new with group_charge_up/down_only

       ! gibo: check for the 'impossible'
       if ((pub_cdft_modes .NE. 1) .AND. (pub_cdft_modes .NE. 2)) then
          call utils_abort(&
               'Iinvalid number of (charge-spin) active mode (neither 1 or 2), &
               &check cdft_atom_charge/spin and &
               &cdft_group_charge/spin_acceptor/donor commands in input file')
       endif

       pub_cdft_group_charge_diff_target = esdf_double( &
            'cdft_group_charge_diff_target',0._DP)
       pub_cdft_group_spin_diff_target   = esdf_double( &
            'cdft_group_spin_diff_target',0._DP)
       pub_cdft_group_charge_acceptor_target = esdf_double( &
            'cdft_charge_acceptor_target',0._DP)
       pub_cdft_group_charge_donor_target    = esdf_double( &
            'cdft_charge_donor_target',0._DP)
       pub_cdft_group_spin_acceptor_target   = esdf_double( &
            'cdft_spin_acceptor_target',0._DP)
       pub_cdft_group_spin_donor_target      = esdf_double( &
            'cdft_spin_donor_target',0._DP)
       pub_maxit_cdft_u_cg         = esdf_integer('maxit_cdft_u_cg',60)
       pub_cdft_cg_type           = esdf_string('cdft_cg_type','NGWF_FLETCHER')
       pub_cdft_cg_threshold  = esdf_double('cdft_cg_threshold',1.0e-3_dp)
       pub_cdft_trial_length      = esdf_double('cdft_trial_length',0.1_dp)
       pub_cdft_cg_max         = esdf_integer('cdft_cg_max', 5)
       pub_cdft_max_grad       = esdf_double('cdft_max_grad', 1.0e-3_dp)
       pub_cdft_elec_energy_tol= esdf_physical('cdft_elec_energy_tol', -0.0001_DP, &
            'hartree')
       !pub_cdft_cg_max_step        = esdf_double('cdft_cg_max_step',8.00_DP)
       pub_cdft_cg_max_step        = esdf_double('cdft_cg_max_step',50.00_DP)
       pub_cdft_guru           = esdf_boolean('cdft_guru', .false.)
       pub_cdft_continuation   = esdf_boolean('cdft_continuation', .false.)
       pub_cdft_tight         = esdf_boolean('cdft_tight', .false.)
       pub_cdft_write_potentials = esdf_boolean('cdft_write_potentials', .true.)
       !gibo: cDFT stuff===== END

       ! ep: mermin section
       pub_mermin                    = esdf_boolean('mermin', .false.)
       pub_check_mermin              = esdf_boolean('check_mermin', .false.)
       pub_mermin_temp               = esdf_boolean('mermin_temp', .false.)
       pub_mermin_smearing_width     = esdf_physical('mermin_smearing_width', &
           3.166811429e-3_DP, 'hartree')
       pub_mermin_cheb               = esdf_integer('mermin_cheb', 11)
       maxit_mermin                  = esdf_integer('maxit_mermin',10)
       mermin_threshold_orig         = esdf_double('mermin_threshold_orig', &
           1.0e-9_dp)
       pub_mermin_cg_type            = esdf_string('mermin_cg_type', &
           'MERMIN_FLETCHER')
       pub_mermin_cg_max_step        = esdf_double('mermin_cg_max_step', &
           3.00_DP)
       !pub_mermin_check_trial_steps = esdf_boolean('mermin_check_trial_steps',&
       !     .false.)
       pub_mermin_check_trial_steps = esdf_boolean('mermin_check_trial_steps', &
            .true.)
       pub_mermin_round_evals        = esdf_integer('mermin_round_evals', -1)
       pub_mermin_free_energy_thres = esdf_physical('mermin_free_energy_thres', &
            1.0e-6_DP,'hartree')
       pub_mermin_mu_sq              = esdf_boolean('mermin_mu_sq',.false.)
       !if (pub_mermin_spin_fix .ge. 0) then
       !   pub_spin_polarised = .true.
       !end if
       !to be added in the future
       ! ep: mermin section

       pub_dft_nu_opt_u1_only   = esdf_boolean('dft_nu_opt_u1_only', .false.)
       pub_dft_nu_opt_u2_only   = esdf_boolean('dft_nu_opt_u2_only', .false.)
       pub_dft_nu_continuation  = esdf_boolean('dft_nu_continuation', .false.)

       ! jd: Changing defaults to 'unspecified' to help clueless users
       pub_coulomb_radius      = esdf_physical('coulomb_cutoff_radius', &
            -1.0_DP,'bohr')
       pub_coulomb_length      = esdf_physical('coulomb_cutoff_length', &
            -1.0_DP,'bohr')
       pub_coulomb_cutoff_type = esdf_string('coulomb_cutoff_type', 'NONE')
       pub_coulomb_cutoff_write_int = esdf_boolean('coulomb_cutoff_write_int',&
            .false.)
       pub_mt_cutoff           = esdf_physical('pbc_correction_cutoff', 0.0_DP,'bohr')

       pub_do_properties       = esdf_boolean('do_properties', .false.)
       pub_num_eigenvalues     = esdf_integer('num_eigenvalues', 10)
       pub_homo_dens_plot      = esdf_integer('homo_dens_plot', -1)
       pub_lumo_dens_plot      = esdf_integer('lumo_dens_plot', -1)
       pub_homo_plot           = esdf_integer('homo_plot', 5)
       pub_lumo_plot           = esdf_integer('lumo_plot', 5)
       pub_dos_smear           = esdf_physical('dos_smear',0.1_DP / &
            HARTREE_IN_EVS, 'hartree')
       pub_popn_calculate      = esdf_boolean('popn_calculate', &
            pub_do_properties)
       pub_lowdin_popn_calculate = esdf_boolean('lowdin_popn_calculate', &
            .false.)
       pub_popn_bond_cutoff    = esdf_physical('popn_bond_cutoff',&
            3.0_DP*ANGSTROM,'bohr')
       pub_popn_mulliken_partial = esdf_boolean('popn_mulliken_partial', &
            .false.)
       pub_ngwf_analysis       = esdf_boolean('ngwf_analysis', .false.)
       pub_polarisation_calculate = esdf_boolean('polarisation_calculate',&
            .false.)
       pub_polarisation_berry = esdf_boolean('polarisation_berry',&
            .false.)
       pub_polarisation_local = esdf_boolean('polarisation_local',&
            .false.)
       pub_polarisation_simcell_calculate = esdf_boolean(&
            'polarisation_simcell_calculate', .false.)
       buf = esdf_string('polarisation_simcell_refpt','0.0 0.0 0.0')
       ! jd: append "bohr" to polarisation_simcell_refpt as default units
       buf = trim(adjustl(buf))//' bohr'
       read(buf,*) pub_polarisation_simcell_refpt(:), tmp_unit
       pub_polarisation_simcell_refpt(:) = pub_polarisation_simcell_refpt(:) * &
            esdf_convfac(tmp_unit,'bohr')
       pub_efield_calculate = esdf_boolean('efield_calculate', .false.)
       pub_spread_calculate    = esdf_boolean('spread_calculate',.false.)
       pub_anharmonic_calculate= esdf_boolean('anharmonic_calculate',.false.)
       pub_print_potential_noxc= esdf_boolean('print_potential_noxc',.false.)

       pub_pdos_optados_output= esdf_boolean('pdos_optados_output',.true.)

       pub_pdos_max_l= esdf_integer('pdos_max_l',-1)

       pub_pdos_max_n = esdf_integer('pdos_max_n',-1)
       pub_pdos_stride_n = 1 ! esdf_integer('pdos_stride_n',1) ! ja531-> Set as devel code for now

       pub_pdos_lowdin = esdf_boolean('pdos_lowdin',.true.)

       pub_pdos_lcao_optimize = esdf_boolean('pdos_lcao_optimize',.false.)
       pub_pdos_orth_atom_blocks = esdf_boolean('pdos_orth_atom_blocks',.false.)
       pub_pdos_output_basis = esdf_boolean('pdos_output_basis',.false.)
       pub_pdos_output_swopt_kernham = esdf_boolean('pdos_output_swopt_kernham',.false.)

       ! ja531-> pseudo truth table of how the following two options should be set
       ! r       pseudoatomic
       ! e      T  | F  |  A |
       ! d  |T| UU | II | IO |
       ! u  |F| II | II | IO |
       ! c  |A| OI | OI | DD |
       ! e
       !    T:true F:false A:absent
       !    U:undefined I:input O:opposite D:default
       !    Eg. IO: reduce=input, pseudoatomic=.not.reduce
       !
       if(esdf_defined('pdos_pseudoatomic')) then
          if(esdf_defined('pdos_reduce_sws')) then
             ! in this case we always take the input
             ! the undefined TT case will be dealt with later
             pub_pdos_pseudoatomic = esdf_boolean('pdos_pseudoatomic',.true.)
             pub_pdos_reduce_sws = esdf_boolean('pdos_reduce_sws',.false.)
          else
             pub_pdos_pseudoatomic = esdf_boolean('pdos_pseudoatomic',.true.)
             pub_pdos_reduce_sws = esdf_boolean('pdos_reduce_sws',.not.pub_pdos_pseudoatomic)
          end if
       else
          if(esdf_defined('pdos_reduce_sws')) then
             pub_pdos_reduce_sws = esdf_boolean('pdos_reduce_sws',.true.)
             pub_pdos_pseudoatomic = esdf_boolean('pdos_pseudoatomic',.not.pub_pdos_reduce_sws)
          else
             pub_pdos_pseudoatomic = esdf_boolean('pdos_pseudoatomic',.true.)
             pub_pdos_reduce_sws = esdf_boolean('pdos_reduce_sws',.false.)
          end if
       end if


       ! ja531-> by default we should construct the basis in PDOS if we are reducing SWs
       pub_pdos_construct_basis = esdf_boolean('pdos_construct_basis',pub_pdos_reduce_sws)
       pub_pdos_sum_mag = esdf_boolean('pdos_sum_mag', .true.)


       pub_geom_modulus_est        = esdf_physical('geom_modulus_est', 0.017_dp, 'ha/bohr**3')
       pub_geom_frequency_est      = esdf_physical('geom_frequency_est', 0.0076_dp, 'hartree')
       pub_geom_energy_tol         = esdf_physical('geom_energy_tol', 1.0e-6_dp, 'hartree')
       pub_geom_force_tol          = esdf_physical('geom_force_tol', 0.002_dp, 'ha/bohr')
       pub_geom_disp_tol           = esdf_physical('geom_disp_tol', 0.005_dp, 'bohr')
       pub_geom_max_iter           = esdf_integer('geom_max_iter', 50)
       pub_geom_convergence_win    = esdf_integer('geom_convergence_win', 2)
       pub_geom_continuation       = esdf_boolean('geom_continuation', .false.)
       pub_geom_backup_iter        = esdf_integer('geom_backup_iter', 1)
       pub_geom_reset_dk_ngwfs_iter= esdf_integer('geom_reset_dk_ngwfs_iter',6)
       pub_geom_lbfgs_max_updates  = esdf_integer('geom_lbfgs_max_updates',30)
       pub_geom_lbfgs_block_length = esdf_integer('geom_lbfgs_block_length',30)
       if (pub_dftb) then
          ! ab: GFN0 is non-self-consistent method, do not reuse dk or ngwfs
          if (pub_dftb_method == DFTB_GFN0) then
             pub_geom_reuse_dk_ngwfs = .false.
          else
             call utils_abort("Unsupported DFTB method. Only GFN0 is supported.")
          end if
       else
          pub_geom_reuse_dk_ngwfs = esdf_boolean('geom_reuse_dk_ngwfs', .true.)
       end if
       pub_geom_print_inv_hessian  = esdf_boolean('geom_print_inv_hessian',.false.)
       pub_geom_method             = esdf_string('geom_method', 'CARTESIAN')
       pub_geom_lbfgs              = esdf_boolean('geom_lbfgs',.false.)

       ! lam81
       pub_geom_precond_type       = esdf_string('geom_precond_type', 'NONE')
       pub_geom_precond_scale_cell = esdf_boolean('geom_precond_scale_cell', .false.)
       pub_geom_precond_exp_c_stab = esdf_double('geom_precond_exp_c_stab', 0.1_dp)
       pub_geom_precond_exp_A      = esdf_double('geom_precond_exp_A', 3.0_dp)
       pub_geom_precond_exp_r_NN   = esdf_physical('geom_precond_exp_r_NN', 0.0_dp, 'bohr')
       pub_geom_precond_exp_r_cut  = esdf_physical('geom_precond_exp_r_cut', 0.0_dp, 'bohr')
       pub_geom_precond_exp_mu     = esdf_physical('geom_precond_exp_mu', 0.0_dp, 'ha/bohr**2')
       pub_geom_precond_ff_c_stab  = esdf_physical('geom_precond_ff_c_stab', 0.1_dp, 'ha/bohr**2')
       pub_geom_precond_ff_r_cut   = esdf_physical('geom_precond_ff_r_cut', 3.8_dp, 'bohr')


       ! Set up LBFGS block allocation size...
       if(pub_geom_lbfgs_max_updates.ne.0) then
          pub_geom_lbfgs_block_length = pub_geom_lbfgs_max_updates
       end if

       buf                     = esdf_string('tssearch_method', 'LSTQST')
       tssearch_method         = buf(1:8)
       buf                     = esdf_string('tssearch_lstqst_protocol', &
            'LSTMAXIMUM')
       tssearch_lstqst_protocol= buf(1:20)
       tssearch_qst_max_iter   = esdf_integer('tssearch_qst_max_iter', 5)
       tssearch_cg_max_iter    = esdf_integer('tssearch_cg_max_iter', 20)
       tssearch_force_tol      = esdf_physical('tssearch_force_tol', 0.005_DP,&
            'ha/bohr')
       tssearch_disp_tol       = esdf_physical('tssearch_disp_tol', 0.01_DP, &
            'bohr')
       tssearch_energy_tol     = esdf_physical('tssearch_energy_tol', 1.0e-5_dp, 'hartree')

       ! kkbd: If we run NEB though geomopt the neb_ci_delay should be one higher
       pub_neb_ci_delay            = esdf_integer('neb_ci_delay', -2)
       pub_neb_print_summary       = esdf_boolean('neb_print_summary',.true.)
       pub_neb_spring_constant     = esdf_physical('neb_spring_constant',0.02_dp,'ha/bohr**2')
       pub_neb_read_xyz            = esdf_boolean('neb_read_xyz',.false.)
       pub_neb_converge_all        = esdf_boolean('neb_converge_all',.false.)
       pub_neb_update_method       = esdf_string('neb_update_method','GLBFGS')
       pub_neb_max_iter            = esdf_integer('neb_max_iter',-1)
       pub_neb_glbfgs_history_size = esdf_integer('neb_glbfgs_history_size',10)
       pub_neb_continuation        = esdf_boolean('neb_continuation',.false.)

       pub_reactant_rootname = esdf_string('reactant_rootname','NONE',force_no_case_change=.true.)
       pub_product_rootname  = esdf_string('product_rootname', 'NONE',force_no_case_change=.true.)
       pub_reactant_energy   = esdf_physical('reactant_energy',100.0_DP,'hartree')
       pub_product_energy    = esdf_physical('product_energy',100.0_DP,'hartree')

       pub_md_properties       = esdf_boolean('md_properties',.false.)
       md_restart              = esdf_boolean('md_restart',.false.)
       md_restart_thermo       = esdf_boolean('md_restart_thermo',.true.)
       md_global_restart       = esdf_boolean('md_global_restart',.false.)
       md_write_out            = esdf_boolean('md_write_out',.true.)
       md_delta_t              = esdf_physical('md_delta_t',40.0_DP,'aut') !~1fs
       md_num_iter             = esdf_integer('md_num_iter',100)
       md_reset_history        = esdf_integer('md_reset_history',100)
       md_write_history        = esdf_integer('md_write_history',-1)
       buf                     = esdf_string('md_aux_rep','ASYM')
       md_aux_rep              = buf(1:5)
       md_aux_dkn_t            = esdf_physical('md_aux_dkn_t',100.0_dp, &
            '1/aut**2')
       md_aux_beren_tc         = esdf_physical('md_aux_beren_tc',413.41105_dp, &
            'aut') !~0.01 ps
       md_lnv_threshold        = esdf_double('md_lnv_threshold',1.0e-9_DP)
       md_ngwf_threshold       = esdf_double('md_ngwf_threshold',2.0e-6_DP)
       md_autocorr             = esdf_boolean('md_autocorr',.false.)

       mts_xi                  = esdf_boolean('mts_xi',.false.)
       mts_nstep               = esdf_integer('mts_nstep',1)
       mts_ngwf_threshold      = esdf_double('mts_ngwf_threshold',5.0e-4_DP)
       mts_lnv_threshold       = esdf_double('mts_lnv_threshold',5.0e-6_DP)
       mts_elec_energy_tol     = esdf_physical('mts_elec_energy_tol',-0.001_DP, &
            'hartree')
       mts_elec_force_tol      = esdf_physical('mts_elec_force_tol',-0.001_DP, &
            'ha/bohr')
       mts_ngwf_max_grad       = esdf_double('mts_ngwf_max_grad',-2.0e-5_dp)
       mts_maxit_ngwf_cg       = esdf_integer('mts_maxit_ngwf_cg',50)
       mts_minit_lnv           = esdf_integer('mts_minit_lnv',5)
       mts_maxit_lnv           = esdf_integer('mts_maxit_lnv',5)
       mts_maxit_pen           = esdf_integer('mts_maxit_pen',3)

       mts_mix_inc             = esdf_boolean('mts_mix_inc',.false.)

       buf                     = esdf_string('mix_dkn_type','NONE')
       mix_dkn_type            = buf(1:20)
       mix_dkn_num             = esdf_integer('mix_dkn_num',0)
       buf                     = esdf_string('mix_dkn_init_type','NONE')
       mix_dkn_init_type       = buf(1:20)
       mix_dkn_init_num        = esdf_integer('mix_dkn_init_num',0)
       mix_dkn_reset           = esdf_integer('mix_dkn_reset',50)
        buf                     = esdf_string('mix_ngwfs_type','NONE')
       mix_ngwfs_type          = buf(1:20)
       mix_ngwfs_num           = esdf_integer('mix_ngwfs_num',0)
       buf                     = esdf_string('mix_ngwfs_init_type','NONE')
       mix_ngwfs_init_type     = buf(1:20)
       mix_ngwfs_init_num      = esdf_integer('mix_ngwfs_init_num',0)
       mix_ngwfs_reset         = esdf_integer('mix_ngwfs_reset',50)
       mix_ngwfs_coeff         = esdf_double('mix_ngwfs_coeff',0.1_dp)
       mix_local_length        = esdf_physical('mix_local_length',10.0_DP,'bohr')
       mix_local_smear         = esdf_physical('mix_local_smear',5.0_DP,'bohr')

       cond_read_denskern      = esdf_boolean('cond_read_denskern',.false.)
       cond_read_tightbox_ngwfs= esdf_boolean('cond_read_tightbox_ngwfs', &
            .false.)
       cond_fixed_shift        = esdf_boolean('cond_fixed_shift',.false.)
       cond_calc_max_eigen     = esdf_boolean('cond_calc_max_eigen',.true.)
       cond_num_states         = esdf_integer('cond_num_states', 0)
       cond_maxit_lnv          = esdf_integer('cond_maxit_lnv', 10)
       cond_minit_lnv          = esdf_integer('cond_minit_lnv', 10)
       cond_kernel_cutoff      = esdf_physical('cond_kernel_cutoff',1000.0_DP,'bohr')
       cond_init_shift         = esdf_physical('cond_init_shift',0.0_DP,'hartree')
       cond_shift_buffer       = esdf_physical('cond_shift_buffer',0.1_DP,'hartree')
       cond_num_extra_states   = esdf_integer('cond_num_extra_states', 0)
       cond_num_extra_its      = esdf_integer('cond_num_extra_its', 0)
       cond_plot_joint_orbitals = esdf_boolean('cond_plot_joint_orbitals',.true.)
       cond_plot_vc_orbitals    = esdf_boolean('cond_plot_vc_orbitals',.false.)
       cond_energy_range        = esdf_physical('cond_energy_range',-1.0_DP,'Ha')
       cond_energy_gap          = esdf_physical('cond_energy_gap',0.001_DP,'Ha')

       pub_spectra_calculate    = esdf_boolean('cond_calc_optical_spectra',.false.)
       pub_calc_mom_mat_els     = esdf_boolean('cond_spec_calc_mom_mat_els',.true.)
       pub_calc_nonloc_comm     = esdf_boolean('cond_spec_calc_nonloc_comm',.true.)
       pub_spec_cont_deriv      = esdf_boolean('cond_spec_cont_deriv',.true.)
       pub_spec_nonloc_fin_diff = esdf_double('cond_spec_nonloc_comm_shift',0.0001_dp)
       pub_spectra_print_mat_els= esdf_boolean('cond_spec_print_mat_els',.false.)
       pub_opt_smear            = esdf_physical('cond_spec_opt_smear',0.1_DP / &
            HARTREE_IN_EVS, 'hartree')
       pub_spec_scissor_op      = esdf_physical('cond_spec_scissor_op',0.0_DP / &
            HARTREE_IN_EVS, 'hartree')

       pub_eels_calculate       = esdf_boolean('cond_calc_eels',.false.)
       pub_eels_fine_projectors = &
            esdf_boolean('cond_eels_fine_projectors',.false.)
       pub_eels_realspace       = &
            esdf_boolean('cond_eels_realspace',.false.)


       ! gcc32:
       pub_lr_phonons_calculate  = esdf_boolean('lr_phonons_calculate', .false.)
       pub_lr_phonons_zero_dim   = esdf_boolean('lr_phonons_zero_dim',.false.)
       pub_lr_phonons_kernel_cutoff=esdf_physical('lr_phonons_kernel_cutoff',&
            1000.0_DP,'bohr')
       pub_lr_phonons_restart = esdf_boolean('lr_phonons_restart',.false.)

       ! tjz07: pub_tddfttest included:
       pub_lr_tddft_calculate    = esdf_boolean('lr_tddft_calculate', .false.)
       pub_lr_tddft_num_states   =esdf_integer('lr_tddft_num_states', 1)
       pub_lr_tddft_cg_threshold =esdf_double('lr_tddft_cg_threshold',&
            1.0e-6_DP)
       pub_lr_tddft_penalty_tol  =esdf_double('lr_tddft_penalty_tol',1.0e-8_DP)
       pub_lr_tddft_maxit_pen    =esdf_integer('lr_tddft_maxit_pen', 20)
       pub_lr_tddft_maxit_cg     =esdf_integer('lr_tddft_maxit_cg', 60)
       pub_lr_tddft_reset_cg     =esdf_integer('lr_tddft_reset_cg', 100)
       pub_lr_tddft_kernel_cutoff=esdf_physical('lr_tddft_kernel_cutoff',&
            1000.0_DP,'bohr')
       pub_lr_tddft_write_densities=esdf_boolean('lr_tddft_write_densities'&
            ,.true.)
       pub_lr_tddft_write_kernels=esdf_boolean('lr_tddft_write_kernels',&
            .false.)
       pub_lr_tddft_restart      =esdf_boolean('lr_tddft_restart',.false.)
       pub_lr_tddft_triplet      =esdf_boolean('lr_tddft_triplet',.false.)
       pub_lr_tddft_projector    =esdf_boolean('lr_tddft_projector',.true.)
       pub_lr_tddft_joint_set    =esdf_boolean('lr_tddft_joint_set',.true.)
       pub_lr_tddft_preopt       =esdf_boolean('lr_tddft_preopt', .false.)
       pub_lr_tddft_preopt_iter  =esdf_integer('lr_tddft_preopt_iter',20)
       pub_lr_tddft_analysis     =esdf_boolean('lr_tddft_analysis',.false.)
       pub_lr_tddft_num_conv_states =esdf_integer('lr_tddft_num_conv_states',0)
       pub_lr_tddft_check_conv_iter =esdf_integer('lr_tddft_check_conv_iter',5)
       pub_lr_tddft_precond      =esdf_boolean('lr_tddft_precond',.true.)
       pub_lr_tddft_precond_tol  =esdf_double('lr_tddft_precond_tol', 1.0e-6_DP)
       pub_lr_tddft_precond_iter =esdf_integer('lr_tddft_precond_iter',20)
       pub_lr_tddft_HOMO_num     =esdf_integer('lr_tddft_HOMO_num',20)
       pub_lr_tddft_LUMO_num     =esdf_integer('lr_tddft_LUMO_num',20)
       pub_lr_tddft_xc_finite_diff =esdf_boolean('lr_tddft_xc_finite_diff',.true.)
       pub_lr_tddft_subsystem_coupling=esdf_boolean('lr_tddft_subsystem_coupling',.false.)
       pub_lr_tddft_restart_from_TDA=esdf_boolean('lr_tddft_restart_from_TDA', .false.)
       pub_lr_tddft_ct_length=esdf_physical('lr_tddft_ct_length',20.0_DP,'bohr')
       pub_lr_tddft_mlwf_analysis=esdf_boolean('lr_tddft_mlwf_analysis',.false.)
       pub_lr_tddft_mom_mat_els=esdf_boolean('lr_tddft_mom_mat_els',.false.)
       pub_lr_tddft_penalty_func=esdf_boolean('lr_tddft_penalty_func',.true.)
       if(pub_is_implicit_solvent) then
          ! do not choose a default optical permittivity even if implict solvent is
          ! set to T. The user has to switch on dynamic screening of the excitation
          ! himself by setting it to a different value
          lr_tddft_default_permittivity = 1.0_DP
       else
          lr_tddft_default_permittivity = 1.0_DP
       end if
       pub_lr_optical_permittivity=esdf_double('lr_optical_permittivity',&
          lr_tddft_default_permittivity)
       pub_lr_tddft_RPA=esdf_boolean('lr_tddft_RPA',.false.)
       pub_lr_tddft_init_random=esdf_boolean('lr_tddft_init_random',.true.)
       pub_lr_tddft_init_max_overlap=esdf_boolean('lr_tddft_init_max_overlap',.false.)
       pub_lr_tddft_spectrum_smear=esdf_physical('lr_tddft_spectrum_smear',0.1_DP / &
            HARTREE_IN_EVS, 'hartree')

       ! jhl52: qnto
       pub_qnto_analysis          =esdf_boolean('qnto_analysis',.false.)
       pub_qnto_num_transition    =esdf_integer('qnto_num_transition',2)
       pub_qnto_write_orbitals    =esdf_boolean('qnto_write_orbitals',.false.)
       pub_qnto_svd_method        =esdf_integer('qnto_svd_method',2)
       pub_qnto_nbo_proj          =esdf_boolean('qnto_nbo_proj',.false.)
       pub_qnto_ref_dir           =esdf_string('qnto_ref_dir','')
       pub_qnto_num_ref_states    =esdf_integer('qnto_num_ref_states',1)
       pub_qnto_num_core_atoms    =esdf_integer('qnto_num_core_atoms',0)

       ! smmd: Transmission coefficients
       pub_etrans_bulk          = esdf_boolean('etrans_bulk',.false.)
       pub_etrans_lcr           = esdf_boolean('etrans_lcr',.false.)
       pub_etrans_same_leads    = esdf_boolean('etrans_same_leads',.false.)
       pub_etrans_write_xyz     = esdf_boolean('etrans_write_xyz',.true.)
       pub_etrans_write_hs      = esdf_boolean('etrans_write_hs',.false.)
       pub_etrans_ecmplx        = esdf_physical('etrans_ecmplx',1.0e-6_DP, 'hartree')
       pub_etrans_emax          = esdf_physical('etrans_emax',0.2_DP, 'hartree')
       pub_etrans_emin          = esdf_physical('etrans_emin',-0.2_DP, 'hartree')
       pub_etrans_enum          = esdf_integer('etrans_enum',50)
       pub_etrans_eref          = esdf_physical('etrans_eref',0.0_DP,'hartree')
       pub_etrans_eref_method   = esdf_string('etrans_eref_method','LEADS')
       pub_etrans_calc_lead_pot = esdf_boolean('etrans_calculate_lead_mu',.false.)
       pub_etrans_lead_nkpoints = esdf_integer('etrans_lead_nkpoints',32)
       pub_etrans_lead_disp_tol = esdf_physical('etrans_lead_disp_tol',1.0_DP,'bohr')
       pub_etrans_num_eigchan   = esdf_integer('etrans_num_eigchan',0)
       pub_etrans_plot_eigchan  = esdf_boolean('etrans_plot_eigchan',.false.)
       pub_etrans_lead_lcheck   = esdf_boolean('etrans_lead_size_check',.true.)
       pub_etrans_seed_lead     = esdf_integer('etrans_seed_lead',1)
       call utils_assert(pub_etrans_emax > pub_etrans_emin, &
         & 'etrans_emax must be greater than etrans_emin')
       call utils_assert(pub_etrans_enum >= 0, 'etrans_enum must be positive')
       ! rab207, 08/07/14: noticed a bug using etrans_same_leads
       ! leading to incorrect transmission. This functionality shouldn't
       ! really be used as it doesn't offer any real computational speed up
       ! and can introduce errors due to misalignement, or different atomic
       ! ordering ==> not worth fixing at the minute
       call utils_assert(.not.pub_etrans_same_leads, &
                         'etrans_same_leads = T does not currently work')

       ! vv: anharmonic parameters(default)
       pub_rand_seed                = esdf_integer('rand_seed',-1)
       pub_anh_qc_factor            = esdf_string('anh_qc_factor','HARMONIC')
       pub_anh_acf_factor           = esdf_string('anh_acf_factor','NONE')
       pub_anh_type                 = esdf_string('anh_type','IR_CALCULATION')
       pub_anh_first_iter           = esdf_integer('anh_first_iter',0)
       pub_anh_last_iter            = esdf_integer('anh_last_iter',0)
       pub_anh_md_temp              = esdf_physical('anh_md_temp',0.0_DP,'kelvin')
       pub_anh_plot_firstfreq       = esdf_physical('anh_plot_firstfreq',0.0_DP,'1/cm')
       pub_anh_plot_lastfreq        = esdf_physical('anh_plot_lastfreq',0.0_DP,'1/cm')
       pub_anh_plot_all             = esdf_boolean('anh_plot_all',.false.)
       pub_anh_apply_filter         = esdf_boolean('anh_apply_filter',.false.)

       ! fc: parameters for phonon pub_task
       pub_phonon_farming_task = esdf_integer('phonon_farming_task',0)
       pub_phonon_sampling = esdf_integer('phonon_sampling',1)
       pub_phonon_vib_free = esdf_integer('phonon_vib_free',7)
       pub_phonon_disp = esdf_physical('phonon_finite_disp',0.1_DP,'bohr')
       pub_phonon_fmax = esdf_physical('phonon_fmax',0.005_DP,'ha/bohr')
       pub_phonon_energy_check = esdf_boolean('phonon_energy_check',.false.)
       pub_phonon_tmin = esdf_physical('phonon_tmin',0.0_DP,'hartree')
       pub_phonon_tmax = esdf_physical('phonon_tmax',0.002_DP,'hartree')
       pub_phonon_deltat = esdf_physical('phonon_deltat',1.5e-5_DP,'hartree')
       pub_phonon_min_freq = esdf_physical('phonon_min_freq',3.6e-06_DP,'hartree')
       pub_phonon_write_eigenvecs = esdf_boolean('phonon_write_eigenvecs',.false.)
       pub_phonon_animate_scale = esdf_double('phonon_animate_scale',1.0_DP)
       pub_phonon_DOS_min = esdf_double('phonon_DOS_min',0.0_DP)
       pub_phonon_DOS_max = esdf_double('phonon_DOS_max',1000.0_DP)
       pub_phonon_DOS_delta = esdf_double('phonon_DOS_delta',10.0_DP)
       pub_phonon_DOS = esdf_boolean('phonon_DOS',.true.)
       pub_phonon_SK = esdf_boolean('phonon_SK',.false.)

       pub_bs_kpoint_path_spacing = esdf_physical('bs_kpoint_path_spacing', &
            0.1889727_DP,'1/bohr')
       pub_bs_num_eigenvalues  = esdf_integer('bs_num_eigenvalues',-1)
       buf                     = esdf_string('bs_unfold','0 0 0')
       read(buf,*) pub_bs_unfold(:)
       pub_bs_method           = esdf_string('bs_method', 'TB')

       ! gcc32: band-structure unfolding
       pub_bsunfld_num_kpts_path  = esdf_integer('bsunfld_num_kpts_path',2)
       pub_bsunfld_num_atoms_prim = esdf_integer('bsunfld_num_atoms_prim',0)
       pub_bsunfld_restart        = esdf_boolean('bsunfld_restart',.false.)
       pub_bsunfld_num_eigenvalues  = esdf_integer('bsunfld_num_eigenvalues',-1)
       buf                     = esdf_string('bsunfld_transformation', &
            '1 0 0 0 1 0 0 0 1') ! identity transform as default
       read(buf,*) pub_bsunfld_transfo(:)

       ! jme: Perturbative Spin-Orbit Couplings in bandstructure.
       pub_perturbative_soc = esdf_boolean('bs_perturbative_soc',.false.)

       pub_do_tddft                   = esdf_boolean('do_tddft', .false.)
       pub_tddft_maximum_energy       = &
            esdf_physical('tddft_maximum_energy',1.0_DP,'hartree')
       pub_tddft_resolution           = &
            esdf_physical('tddft_resolution',0.001_DP,'hartree')
       buf = esdf_string('tddft_propagation_method', 'CRANKNICHOLSON')
       pub_tddft_propagation_method   = buf(1:20)
       pub_tddft_sparsity_level       = esdf_integer('tddft_sparsity_level', 0)
       pub_tddft_tammdancoff          = &
            esdf_boolean('tddft_tammdancoff', .false.)
       tdbuf                          = &
            esdf_string('tddft_dipole_kick_strength','0.0 0.0 0.0')
       read(tdbuf,*) pub_tddft_dipole_kick_strength(:)
       pub_tddft_xc_functional        = esdf_string('tddft_xc_functional','LDA')
       pub_tddft_hamiltonian_mixing   = &
            esdf_integer('tddft_hamiltonian_mixing', 0)
       pub_tddft_damping              = &
            esdf_physical('tddft_damping',0.0_DP,'hartree')
       pub_tddft_enforced_idempotency = &
            esdf_boolean('tddft_enforced_idempotency', .false.)
       pub_tddft_maxit_hotelling      = &
            esdf_integer('tddft_maxit_hotelling', 50)
       pub_tddft_max_resid_hotelling  = &
            esdf_double('tddft_max_resid_hotelling', 1.0e-18_DP)
       pub_tddft_inv_overlap_exact    = &
            esdf_boolean('tddft_inv_overlap_exact', .true.)

       ! jd: --- QM/MM ---
       pub_pol_emb_pot_filename = esdf_string('pol_emb_pot_filename', &
            '*undefined*', force_no_case_change = .true.)
       ! jd: Turn on polarisable embedding potential automatically if a filename
       !     has been specified.
       pub_pol_emb_pot = (pub_pol_emb_pot_filename(1:1) /= '*')
       pub_pol_emb_qmstar = esdf_boolean('pol_emb_qmstar', pub_pol_emb_pot)
       pub_pol_emb_write_vacuum_restart = &
            esdf_boolean('pol_emb_write_vacuum_restart', .false.)
       pub_pol_emb_dma_min_l = esdf_integer('pol_emb_dma_min_l',0)
       pub_pol_emb_dma_max_l = esdf_integer('pol_emb_dma_max_l',-1)
       pub_pol_emb_vacuum_dma_min_l = esdf_integer('pol_emb_vacuum_dma_min_l',0)
       pub_pol_emb_vacuum_dma_max_l = esdf_integer('pol_emb_vacuum_dma_max_l',-1)
       pub_pol_emb_vacuum_qmstar = esdf_boolean('pol_emb_vacuum_qmstar', &
             pub_pol_emb_qmstar .and. (pub_pol_emb_vacuum_dma_max_l >=0) .and. &
             .not. pub_pol_emb_write_vacuum_restart)
       pub_pol_emb_mpole_exclusion_radius = esdf_physical(&
            'pol_emb_mpole_exclusion_radius', 0.25_DP, 'bohr')
       ! jd: Rationale for default. Average polarisability in amoeba09 is
       !     1.059A^3. (\alpha_i \alpha_j) ^ (1/6) yields 1.01929 A for this
       !     value, which is 1.92618 a0.
       pub_pol_emb_pairwise_polarisability = &
            esdf_physical('pol_emb_pairwise_polarisability', 1.92618_DP, 'bohr')
       ! jd: Thole damping as per AMOEBA
       pub_pol_emb_thole_a = esdf_double('pol_emb_thole_a', 0.39_DP)
       pub_pol_emb_smearing_a = &
            esdf_physical('pol_emb_smearing_a', 0.2_DP, 'bohr')
       pub_pol_emb_perm_scaling = &
            esdf_double('pol_emb_perm_scaling', 1.0_DP)
       pub_pol_emb_fixed_charge = esdf_boolean('pol_emb_fixed_charge', .false.)
       pub_pol_emb_dbl_grid = esdf_boolean('pol_emb_dbl_grid', .false.)
       pub_pol_emb_repulsive_mm_pot_write = esdf_boolean(&
            'pol_emb_repulsive_mm_pot_write', .false.)
       pub_pol_emb_repulsive_mm_pot_verbose = esdf_boolean(&
            'pol_emb_repulsive_mm_pot_verbose', .false.)
       pub_pol_emb_repulsive_mm_pot_cutoff = esdf_physical(&
            'pol_emb_repulsive_mm_pot_cutoff', 10.0_DP, 'bohr')
       pub_pol_emb_repulsive_mm_pot_a = esdf_double(&
            'pol_emb_repulsive_mm_pot_a', 6.97_DP)
       pub_pol_emb_repulsive_mm_pot_b = esdf_double(&
            'pol_emb_repulsive_mm_pot_b', -11.87_DP)
       pub_pol_emb_repulsive_mm_pot_c = esdf_double(&
            'pol_emb_repulsive_mm_pot_c', 5.64_DP)
       pub_pol_emb_repulsive_mm_pot_alpha = esdf_double(&
            'pol_emb_repulsive_mm_pot_alpha', 146869.0_DP)
       pub_pol_emb_repulsive_mm_pot_beta = esdf_double(&
            'pol_emb_repulsive_mm_pot_beta', 8.897_DP)
       pub_pol_emb_repulsive_mm_pot_r0 = esdf_physical(&
            'pol_emb_repulsive_mm_pot_r0', 7.804568_DP, 'bohr')
       pub_pol_emb_polscal = esdf_double('pol_emb_polscal', 1.0_DP)

       if(pub_pol_emb_pot) then
          call utils_assert(.not. pub_write_xyz, &
               '"write_xyz T" cannot be used with pol_emb_pot (your &
               &ONETEP .xyz file would clash with the externally generated &
               &.xyz file).')
          call utils_assert(pub_pol_emb_dma_max_l /= -1, "'pol_emb_dma_max_l' &
               &must be specified when using polarisable embedding")
       else
          call utils_assert(.not. pub_pol_emb_qmstar, &
               '"pol_emb_qmstar" requires polarisable embedding to be &
               &turned on (must specify a "pol_emb_pot_filename")')
          call utils_assert(.not. pub_pol_emb_vacuum_qmstar, &
               '"pol_emb_vacuum_qmstar" requires polarisable embedding to be &
               &turned on (must specify a "pol_emb_pot_filename")')
       end if

       ! jd: -------------

       ! jd: Implicit solvation
       pub_is_implicit_solvent        = esdf_boolean('is_implicit_solvent', &
            .false.)
       pub_is_smeared_ion_rep         = esdf_boolean('is_smeared_ion_rep', &
             .false.)
       pub_is_separate_restart_files  = esdf_boolean('is_separate_restart_files', &
            .false.)
       pub_is_auto_solvation          = esdf_boolean('is_auto_solvation', &
            .false.)
       if(pub_is_auto_solvation) then
          ! jd: Important defaults change depending on whether is_implicit_solvent
          !     is T or F: is_include_apolar, is_bulk_permittivity. Follow
          !     the principle of least surprise when auto solvation is on, but
          !     is_implicit_solvent is omitted, and turn it on. Otherwise we'll
          !     get no apolar term and a permittivity of 1 in the solvent calc.
          pub_is_implicit_solvent = .true.
       end if
       pub_is_include_apolar      = esdf_boolean('is_include_apolar', &
            pub_is_implicit_solvent)
       pub_is_apolar_sasa_definition= esdf_string('is_apolar_sasa_definition', &
            'DENSITY')
       pub_is_apolar_method       = esdf_string('is_apolar_method', 'SASA')
       pub_is_steric_write            = esdf_boolean('is_steric_write', .false.)
       pub_is_pbe                     = esdf_string('is_pbe', 'NONE')
       if(pub_is_pbe == 'LINEARIZED') then
          pub_is_pbe = 'LINEARISED'
       end if
       if(pub_is_pbe == 'NONE') then
          pub_is_pbe_bc_debye_screening  = esdf_boolean(&
               'is_pbe_bc_debye_screening', .false.)
       else
          pub_is_pbe_bc_debye_screening  = esdf_boolean(&
               'is_pbe_bc_debye_screening', .true.)
       end if
                                      ! jd: 0.0 will be converted to max_expcap
       pub_is_pbe_exp_cap             = esdf_double('is_pbe_exp_cap', 0.0_DP)

       ! jd: Tolerances for differences between two formally equivalent energy
       !     expressions. In FULL they should match to numerical precision.
       !     In LINEARISED, at least in PBCs, some things only match to 1st order
       pub_is_pbe_energy_tolerance = -1.0_DP ! jd: Handles "NONE"
       if(pub_is_pbe == "FULL") then
          pub_is_pbe_energy_tolerance = esdf_physical(&
               'is_pbe_energy_tolerance', 1.5936D-6, 'hartree') ! 0.001 kcal/mol
       end if
       if(pub_is_pbe == "LINEARISED") then
          pub_is_pbe_energy_tolerance = esdf_physical(&
               'is_pbe_energy_tolerance', 0.00007968_DP, 'hartree') ! 0.05 kcal/mol
       end if

       ! jd: We need a BC-dependent default, but BCs will only be known later.
       !     Defer this and let check_inputs() adjust it.
       pub_is_pbe_neutralisation_scheme = &
            esdf_string('is_pbe_neutralisation_scheme', '*undefined*')

       ! jd: Default bulk permittivity is 78.54 if using implicit solvent,
       !     1.0 otherwise
       if(pub_is_implicit_solvent) then
          is_default_bulk_permittivity = 78.54_DP
       else
          is_default_bulk_permittivity = 1.0_DP
       end if
       pub_is_bulk_permittivity = esdf_double('is_bulk_permittivity', &
            is_default_bulk_permittivity)
       pub_is_smeared_ion_width = esdf_physical('is_smeared_ion_width', 0.8_DP,&
            'bohr')
       pub_is_core_width = esdf_physical('is_core_width', 1.2_DP,'bohr')
       pub_is_surface_thickness = esdf_double('is_surface_thickness', 0.0002_DP)
       ! jd: Default surface tension is an estimate for H2O from
       !     J.Chem.Phys.124 074103 (2006).
       pub_is_solvent_surf_tension = &
            esdf_physical('is_solvent_surf_tension', 4.7624D-5,'ha/bohr**2')
       pub_is_solvent_pressure = &
            esdf_physical('is_solvent_pressure', -1.1896D-5,'ha/bohr**3')
       ! jd: Default scaling factor converts the physical surface tension for
       !     H2O to the effective value, taking dis-rep into account.
       pub_is_apolar_scaling_factor = esdf_double(&
            'is_apolar_scaling_factor', 0.281075_DP)
       ! jd: Defaults for beta and rho0 now correspond to ONETEP's model.
       pub_is_density_threshold = esdf_double('is_density_threshold',0.00035_DP)
       pub_is_solvation_beta = esdf_double('is_solvation_beta',1.3_DP)
       ! jd: Defaults for rho_min and rho_max follow Andreussi's paper
       pub_is_density_min_threshold = esdf_double('is_density_min_threshold',0.0001_DP)
       pub_is_density_max_threshold = esdf_double('is_density_max_threshold',0.0050_DP)
       ! jd: Defaults for soft-sphere model
       pub_is_soft_sphere_delta = esdf_double('is_soft_sphere_delta',0.5_DP)
       pub_is_soft_sphere_scale = esdf_double('is_soft_sphere_scale',1.33_DP)

       if(pub_is_pbe == 'NONE') then
          pub_is_pbe_temperature = esdf_double(&
               'is_pbe_temperature', -1.0_DP)
       else
          pub_is_pbe_temperature = esdf_double(&
               'is_pbe_temperature', 300.0_DP)
       end if
       pub_is_sc_steric_magnitude = esdf_physical('is_sc_steric_magnitude', &
            -1.0_DP, 'ha*bohr**12')
       pub_is_sc_steric_cutoff = esdf_physical('is_sc_steric_cutoff', -1.0_DP, &
            'bohr')
       pub_is_sc_steric_smoothing_alpha = esdf_physical(&
            'is_sc_steric_smoothing_alpha', 1.5_DP, '1/bohr')
       pub_is_steric_pot_type = esdf_string('is_steric_pot_type', 'X')
       pub_is_hc_steric_smearing = esdf_physical('is_hc_steric_smearing', &
            0.4_DP, 'bohr')
       pub_is_hc_steric_dens_isovalue = &
            esdf_double('is_hc_steric_dens_isovalue', 0.0030_DP)

       if(esdf_block('species_solvent_radius',nlines) .or. &
            esdf_defined('is_hc_steric_smearing') .or. &
            esdf_defined('is_hc_steric_dens_isovalue')) then
          call utils_assert(pub_is_steric_pot_type /= 'X', &
               myself//': You cannot have no steric potential ("is_steric_pot_&
               &type" absent or set to "X") while specifying the details of the&
               & steric potential ("is_hc_steric_smearing" or "is_hc_steric_&
               &dens_isovalue" or "species_solvent_radius" at the same time. &
               &You probably forgot "is_steric_pot_type M".')
       end if


       pub_is_multigrid_verbose = esdf_boolean('is_multigrid_verbose', .false.)
       pub_is_multigrid_verbose_y = esdf_physical('is_multigrid_verbose_y', &
            0.0_DP, 'bohr')
       pub_is_multigrid_verbose_z = esdf_physical('is_multigrid_verbose_z', &
            0.0_DP, 'bohr')
       ! jd: New default for nlevels. This is no longer directly used by the
       !     MG solver, but still affects FD grid granularity. With PBCs,
       !     MG grids must be odd-sized or else there is just 1 MG level.
       !     Moreover, when granularity is insufficient, DL_MG copes by reducing
       !     the number of levels. With 2 levels performance is abysmal, with 3
       !     it's much better. With 5 levels and beyond the gain is minimal.
       !     By setting a default value for nlevels to 2, we get 3 (!) DL_MG
       !     levels, while only carving out at most 7 fine grid slices, which
       !     is an acceptable compromise (previous setting carved up to 31).
       pub_is_bc_coarseness = esdf_integer('is_bc_coarseness', 5)
       pub_is_bc_surface_coarseness = esdf_integer('is_bc_surface_coarseness', &
            1)
       pub_is_bc_threshold = esdf_double('is_bc_threshold', 1D-9)
       pub_is_bc_allow_frac_charge = esdf_boolean('is_bc_allow_frac_charge', &
            .false.)
       pub_is_check_solv_energy_grad = &
            esdf_boolean('is_check_solv_energy_grad', .false.)
       pub_is_solvation_method = esdf_string('is_solvation_method', 'direct')
       pub_is_dielectric_model = esdf_string('is_dielectric_model', &
            'fix_initial')
       pub_is_dielectric_function = esdf_string('is_dielectric_function', 'fgf')
       pub_is_solvation_output_detail = &
            esdf_string('is_solvation_output_detail', 'none')

       pub_is_dielectric_exclusions_smearing = esdf_physical('is_dielectric_exclusions_smear',0.0_DP,'bohr')
       pub_is_solvation_properties = esdf_boolean('is_solvation_properties', &
            .false.)

       ! jd: Turn on smeared ions automatically if using implicit solvent
       if (pub_is_implicit_solvent .and. .not. pub_is_smeared_ion_rep) then
          write(stdout,'(/a)') 'WARNING: is_smeared_ion_rep was automatically &
               &set to T because you specified'//CRLF//'         is_implicit_&
               &solvent T (or is_auto_solvation T), and a solvation'//CRLF//&
               '         calculation cannot be performed without the smeared-&
               &ion representation.'&
               //CRLF//'         Add ''is_smeared_ion_rep T'' to get rid&
               & of this warning.'
          pub_is_smeared_ion_rep = .true.
       end if

       ! jcap: When doing an auto-solvation IS calculation, toggle for
       ! whether the vacuum calculation should be restarted from
       ! vacuum_* files or not
       pub_is_restart_vac_from_vac = esdf_boolean('is_restart_vac_from_vac', .false.)
       ! jcap: Select whether the EMFT kernel or the normal kernel is
       ! used to define the cavity in an IS calculation
       pub_is_emft_cavity = esdf_boolean('is_emft_cavity', .false.)

       ! gab: Warn the user that the newer versions of implicit solvent use a
       !      different defition of the SASA used to calculate the nonpolar
       !      term in the solvation free energy.
       if ((.not.esdf_defined('is_apolar_sasa_definition')) &
            .and. (.not. pub_is_dielectric_function=='SOFT_SPHERE') &
            .and. pub_is_include_apolar) then
          write(stdout,'(/a)') 'WARNING: &
               &the method used to calculate the nonpolar term in implicit &
               &'//CRLF//'         solvent was changed in version 6.1.3. &
               &Therefore the total &
               &'//CRLF//'         energy values for implicit solvent may not &
               &be compatible with &
               &'//CRLF//'         earlier versions. &
               &'//CRLF//CRLF//'         For backwards compatibility,&
               & set is_apolar_sasa_definition to &
               &'//CRLF//'         ''isodensity''. &
               &'//CRLF//CRLF//'         To use the new definition and suppress&
               & this warning, '//CRLF//'&
               &         set is_apolar_sasa_definition to ''density''.'
       end if

       ! jd: Trap use of old-style convention for input parameters, that is
       !     the use of 'is_solvent_surface_tension' or 'is_include_cavitation'.
       obsolete_solvent_surface_tension = esdf_physical(&
            'is_solvent_surface_tension', garbage_real, &
            'ha/bohr**2')
       if(obsolete_solvent_surface_tension /= garbage_real .or. &
            esdf_boolean('is_include_cavitation',.false.) .or. &
            .not. esdf_boolean('is_include_cavitation',.true.)) then
          errmsg = "Your input file is using obsolete syntax for implicit &
               &solvation parameters,"//CRLF//"where the solvent surface &
               &tension used to be pre-scaled to account for"//CRLF//"&
               &solute-solvent dispersion-repulsion. This pre-scaled value &
               &used to be specified"//CRLF//"with the keyword &
               &'is_solvent_surface_tension'."//CRLF//CRLF//"In this version &
               &of ONETEP, the solvent surface tension is specified as &
               &the"//CRLF//"actual physical value that can be obtained from &
               &physical tables, and the "//CRLF//"scaling factor for including&
               & dispersion-repulsion is specified separately."//CRLF//"The &
               &surface tension is specified via 'is_solvent_surf_tension' &
               &(and not"//CRLF//"'is_solvent_surface_tension'), and the &
               &scaling factor -- via"//CRLF//"'is_apolar_scaling_factor'&
               &."//CRLF//CRLF//"Also, 'is_include_cavitation' has been renamed&
               & to 'is_include_apolar'."//CRLF//CRLF//"In the most common &
               &scenario of water solvent, &
               &simply adjust your input file"//CRLF//"from"//CRLF//CRLF//"is_&
               &solvent_surface_tension 0.0000133859 ha/bohr**2"//CRLF//"is_&
               &include_cavitation T"//CRLF//CRLF//"to"//CRLF//CRLF//"is_&
               &solvent_surf_tension 0.07415 N/m"//CRLF//"is_&
               &apolar_scaling_factor 0.281075"//CRLF//"is_include_apolar &
               &T"//CRLF//CRLF//"or simply drop 'is_solvent_surface_tension' &
               &and 'is_include_cavitation'"//CRLF//"altogether, since the abov&
               &e are"//CRLF//"now defaults."//CRLF//CRLF//"If you do not want &
               &to include the solute-solvent dispersion-repulsion &
               &energy"//CRLF//"term, be sure to specify 'is_apolar_scaling_&
               &factor 1.0'."//CRLF//CRLF//"For more details see:"//CRLF//"www.&
               &onetep.org/onetep/pmwiki/uploads/Main/Documentation/implicit_&
               &solvation_changes.pdf"//CRLF//CRLF//"Please update your input&
               & file now."//CRLF
          write(stdout,'(a)') CRLF//errmsg
          call utils_abort(errmsg)
       end if

       ! JCW: is_* keywords made obsolete in update to DL_MG v2.0
       obsolete_parameter_found_DLMGv2 = .false.
       obsolete_parameter_found_DLMGv2 = obsolete_parameter_found_DLMGv2 .or. &
           ( garbage_int /= esdf_integer('is_multigrid_nlevels', garbage_int) )
       obsolete_parameter_found_DLMGv2 = obsolete_parameter_found_DLMGv2 .or. &
           ( garbage_int /= esdf_integer('is_discretization_order', garbage_int) )
       obsolete_parameter_found_DLMGv2 = obsolete_parameter_found_DLMGv2 .or. &
           ( garbage_int /= esdf_integer('is_multigrid_max_iters', garbage_int) )
       obsolete_parameter_found_DLMGv2 = obsolete_parameter_found_DLMGv2 .or. &
           ( esdf_boolean('is_multigrid_error_damping', .false.) .or. &
           .not. esdf_boolean('is_multigrid_error_damping', .true.) )
           ! if either of these evaluate to .true., then keyword must have been
           ! set by the user, but if they both evaluate to false, then the
           ! keyword must not have been set
       obsolete_parameter_found_DLMGv2 = obsolete_parameter_found_DLMGv2 .or. &
           ( garbage_real /= esdf_double('is_multigrid_defect_error_tol', garbage_real) )
       obsolete_parameter_found_DLMGv2 = obsolete_parameter_found_DLMGv2 .or. &
           ( garbage_real /= esdf_double('is_multigrid_error_tol', garbage_real) )
       obsolete_parameter_found_DLMGv2 = obsolete_parameter_found_DLMGv2 .or. &
           ( esdf_boolean('is_pbe_use_fas',.false.) .or. &
            .not. esdf_boolean('is_pbe_use_fas',.true.) )
           ! if either of these evaluate to .true., then keyword must have been
           ! set by the user, but if they both evaluate to false, then the
           ! keyword must not have been set
       if (obsolete_parameter_found_DLMGv2) then
          ! JCW: One of the parameters provided by the user was made obsolete
          ! JCW: when updating the ONETEP / DL_MG interface for the DL_MG v2.0 API
          errmsg = &  ! 1255 characters
              &"Your input file is using obsolete syntax for keywords which &
              &control the "//CRLF//&
              "multigrid solver. The obsolete keywords are as follows:"//CRLF//&
              "    is_multigrid_nlevels,"//CRLF//&
              "    is_discretization_order,"//CRLF//&
              "    is_multigrid_max_iters,"//CRLF//&
              "    is_multigrid_error_damping,"//CRLF//&
              "    is_multigrid_defect_error_tol,"//CRLF//&
              "    is_multigrid_error_tol,"//CRLF//&
              "    is_pbe_use_fas."//CRLF//&
              "The following new keywords have been introduced:"//CRLF//&
              "    mg_granularity_power,"//CRLF//&
              "    mg_defco_fd_order,"//CRLF//&
              "    mg_tol_res_rel,"//CRLF//&
              "    mg_tol_res_abs,"//CRLF//&
              "    mg_tol_cg_res_rel,"//CRLF//&
              "    mg_tol_cg_res_abs,"//CRLF//&
              "    mg_tol_pot_rel,"//CRLF//&
              "    mg_tol_pot_abs,"//CRLF//&
              "    mg_tol_vcyc_rel,"//CRLF//&
              "    mg_tol_vcyc_abs,"//CRLF//&
              "    mg_tol_newton_rel,"//CRLF//&
              "    mg_tol_newton_abs,"//CRLF//&
              "    mg_vcyc_smoother_iter_pre,"//CRLF//&
              "    mg_vcyc_smoother_iter_post,"//CRLF//&
              "    mg_tol_mu_rel,"//CRLF//&
              "    mg_tol_mu_abs,"//CRLF//&
              "    mg_max_res_ratio,"//CRLF//&
              "    mg_max_iters_defco,"//CRLF//&
              "    mg_max_iters_vcycle,"//CRLF//&
              "    mg_max_iters_newton,"//CRLF//&
              "    mg_max_iters_cg,"//CRLF//&
              "    mg_use_error_damping,"//CRLF//&
              "    mg_pbe_use_fas."//CRLF//&
              "These are expert keywords for finely controlling the behaviour &
              &of the multigrid"//CRLF//&
              "solver. If you need to finely control the multigrid solver, you &
              &should remove"//CRLF//&
              "obsolete keywords from your input file and replace with the &
              &appropriate new "//CRLF//&
              "keywords. See the documentation for the implicit solvent &
              &model (online at "//CRLF//&
              "www.onetep.org or packaged with the ONETEP source code) for &
              &details of these"//CRLF//&
              "new keywords. Finally, if you want to control the order of &
              &finite differences"//CRLF//&
              "used elsewhere in ONETEP (not in the defect correction), you &
              &can use the new"//CRLF//&
              "keyword finite_difference_order."//CRLF
          call utils_abort(errmsg)
       end if
       ! jd: --------- End of implicit solvent stuff ---------


       ! JCW:--------- Multigrid solver keywords ---------
       ! JCW: Updated for DL_MG v2.0 API
       ! JCW: pub_mg_granularity_power replaces pub_is_multigrid_nlevels with
       ! JCW: pub_mg_granularity_power = pub_is_multigrid_nlevels + 1
       ! JCW: Setting this = 3 as default leads to usual default granulariy of
       ! JCW: 2**3 = 8.
       pub_mg_granularity_power = esdf_integer('mg_granularity_power', 3)
       ! JCW: pub_mg_defco_fd_order replaces pub_is_discretization_order
       pub_mg_defco_fd_order = esdf_integer('mg_defco_fd_order', 8)
       ! JCW: mg_tol_* keyword defaults are derived from the defaults in DL_MG's
       ! JCW: dl_mg_convergence_params module
       ! JCW: These replace the old keywords
       ! JCW:   pub_is_multigrid_error_tol
       ! JCW:   pub_is_multigrid_defect_error_t
       pub_mg_tol_res_rel = esdf_double('mg_tol_res_rel', 1.0e-2_DP)
       pub_mg_tol_res_abs = esdf_double('mg_tol_res_abs', 5.0e-2_DP)
          ! => DL_MG default tol_res_abs = 1.0e-2, but this was found to be too tight
          !    for the existing ONETEP QC tests, causing some of the solvation tests
          !    to fail due to non-convergence defect correction. Raising this to
          !    5.0e-2 allows the tests to pass, as with earlier versions of DL_MG.
          ! => The convergence test for the defect correction in DL_MG v2.0 is
          !       error_norm < max( tol_pot_abs, tol_pot_rel * pot_norm ) . and.
          !       defect_norm < max( tol_res_abs, tol_res_rel* defect_0_norm )
          !    where pot_norm is the norm of the potential for the previous
          !    iteration of the defect loop and defect_0_norm is the norm of the
          !    defect computed from the second order solver (no defect correction).
          ! => pot_norm is typically ~1.0e+1 in solvent and ~1.0e+2 or
          !    ~1.0e+3 in vacuum for the existing QC tests. This means that with
          !    the default tol_pot_abs and tol_pot_rel, tol_pot_rel will determine
          !    the convergence of the error_norm in all cases and will represent
          !    about the same, or a slightly looser tolerance than in the previous
          !    version of the ONETEP/DL_MG interface where convergence of the
          !    error_norm was determined only by an absolute tolerance, with
          !    default value 1.0e-5.
       pub_mg_tol_pot_rel = esdf_double('mg_tol_pot_rel', 1.0e-6_DP)
       pub_mg_tol_pot_abs = esdf_double('mg_tol_pot_abs', 1.0e-6_DP)
       pub_mg_tol_vcyc_rel = esdf_double('mg_tol_vcyc_rel',1.0e-8_DP)
       pub_mg_tol_vcyc_abs = esdf_double('mg_tol_vcyc_abs',1.0e-5_DP)
       pub_mg_tol_newton_rel = esdf_double('mg_tol_newton_rel',1.0e-8_DP)
       pub_mg_tol_newton_abs = esdf_double('mg_tol_newton_abs',1.0e-5_DP)
       pub_mg_max_res_ratio = esdf_double('mg_max_res_ratio',0.999_DP)
       pub_mg_vcyc_smoother_iter_pre = esdf_integer('mg_vcyc_smoother_iter_pre',2)
       pub_mg_vcyc_smoother_iter_post = esdf_integer('mg_vcyc_smoother_iter_post',1)
       pub_mg_tol_mu_rel = esdf_double('mg_tol_mu_rel',1.0e-3_DP)
       pub_mg_tol_mu_abs = esdf_double('mg_tol_mu_abs',1.0e-3_DP)
       ! JCW: When using DL_MG only as a second order solver (using the old
       ! JCW: defect correction code in ONETEP), via the devel_code
       ! JCW: MG:USE_ONETEP_DEFCO=T:MG then the tolerances used are:
       ! JCW:   pub_mg_tol_pot_abs: error norm convergence threshold
       ! JCW:     --> previously pub_is_multigrid_error_tol
       ! JCW:   pub_mg_tol_res_abs: tolerance parameter used in initial second
       ! JCW:                       solution to Poisson equation
       ! JCW:     --> previously pub_is_multigrid_defect_error_t
       ! JCW: ----------------------------------------------------------------
       ! JCW: mg_max_iters_* keyword defaults are derived from the defaults in
       ! JCW: DL_MG's dl_mg_convergence_params module
       ! JCW: These replace the old keyword
       ! JCW:   pub_is_multigrid_max_iters
       pub_mg_max_iters_defco  = esdf_integer('mg_max_iters_defco',30)
       pub_mg_max_iters_vcycle = esdf_integer('mg_max_iters_vcycle',50)
       pub_mg_max_iters_newton = esdf_integer('mg_max_iters_newton',30)
       pub_mg_max_iters_cg     = esdf_integer('mg_max_iters_cg',50)
       ! JCW: Error damping is automatically activated for PB equation, otherwise
       ! JCW: it is not used unless specifically requested
       ! JCW: This replaces the old keyword
       ! JCW:   pub_is_multigrid_error_damping
       pub_mg_use_error_damping = esdf_boolean('mg_use_error_damping', &
            merge(.true., .false., pub_is_pbe == 'FULL'))
       ! JCW: Full Approximation Scheme is not used (Newton method used by
       ! JCW: default for PB equation) unless specifically requested
       pub_mg_pbe_use_fas = esdf_boolean('mg_use_fas',.false.)

       pub_mg_continue_on_error = esdf_boolean('mg_continue_on_error', &
            merge(.true., .false., pub_is_pbe /= 'NONE'))

       ! jd: Introduced in v6.1.3.6
       pub_mg_tol_cg_res_rel = esdf_double('mg_tol_cg_res_rel', 1.0e-2_DP)
       pub_mg_tol_cg_res_abs = esdf_double('mg_tol_cg_res_abs', 5.0e-2_DP)

       pub_mg_use_cg  = esdf_boolean('mg_use_cg', .false.)

       ! JCW:----- End of multigrid solver keywords ------

       ! JCW:--------- Finite differences keywords ---------
       ! JCW: Control the internal finite_differences module, but do not affect
       ! JCW: the finite differences used in the defect correction (controlled
       ! JCW: by pub_mg_defco_fd_order.
       ! JCW:
       ! JCW: To avoid changing previous default behaviour for implicit solvent
       ! JCW: calculations, pub_finite_difference_order is set equal to
       ! JCW: pub_mg_defco_fd_order if the user does not explicitly set a value.
       ! JCW: This means that if the user does not set any FD order parameters
       ! JCW: then both will have the default value of pub_mg_defco_fd_order (8).
       ! JCW: If the user only sets pub_mg_defco_fd_order, then both will have the
       ! JCW: value of pub_mg_defco_fd_order, but if the user sets
       ! JCW: finite_difference_order, then this can be varied independently.
       ! JCW: This is necessary since the internal finite_differences module is
       ! JCW: used in other components of the solvent model (e.g. non-polar
       ! JCW: contribution).
       ! JCW: **This must be assigned AFTER pub_mg_defco_fd_order**
       pub_finite_difference_order = esdf_integer('finite_difference_order',&
           pub_mg_defco_fd_order)
       ! JCW: Note that the order of FD used in the defect correction when the
       ! JCW: devel_code MG:USE_ONETEP_DEFCO=T:MG is used (to enable the old
       ! JCW: ONETEP implementation of the defect correction) is determined by
       ! JCW: pub_mg_defco_fd_order, NOT pub_finite_difference_order.
       ! JCW:--------- End of finite differences keywords ---------

       ! jd: --------- Open BC keywords ---------
       pub_openbc_ion_ion = esdf_boolean('openbc_ion_ion', .false.)
       pub_openbc_hartree = esdf_boolean('openbc_hartree', .false.)
       pub_openbc_pspot = esdf_boolean('openbc_pspot', .false.)
       pub_openbc_pspot_finetune_nptsx = &
            esdf_integer('openbc_pspot_finetune_nptsx', 100000)
       pub_openbc_pspot_finetune_alpha = &
            esdf_double('openbc_pspot_finetune_alpha', 0.3_DP)
       pub_openbc_pspot_finetune_f = &
            esdf_integer('openbc_pspot_finetune_f', -1)
       ! jd: --------- End of open BC stuff ---------

       ! JCW:--------- Periodic/mixed BC keywords ---------
       pub_multigrid_bc   = esdf_string('multigrid_bc','')
       pub_ion_ion_bc     = esdf_string('ion_ion_bc','')
       pub_pspot_bc       = esdf_string('pspot_bc','')
       pub_smeared_ion_bc = esdf_string('smeared_ion_bc','')
       pub_external_bc_from_cube = esdf_boolean('external_bc_from_cube',.false.)
       pub_vdw_bc         = esdf_string('vdw_bc','')
       pub_dftb_bc        = esdf_string('dftb_bc','')
       pub_hfx_bc         = esdf_string('hfx_bc','')
       ! JCW:--------- End of periodic/mixed BC keywords ---------

       pub_permit_unusual_ngwf_count = &
            esdf_boolean('permit_unusual_ngwf_count', .false.)

       ! ars: Select between tightbox and spherical waves representation when
       ! ars: doing restart
       if (pub_read_sw_ngwfs.and.pub_read_tightbox_ngwfs) then
          ! ars : the program uses tightbox_ngwfs unless we specify in the
          ! ars : input read_tightbox_ngwfs = FALSE pub_read_sw_ngwfs = TRUE
          write(stdout,'(/a)') 'WARNING: pub_read_sw_ngwfs set to FALSE. &
               &Restart will be done using .tightbox_ngwfs file'
          pub_read_sw_ngwfs=.false.
       endif

       ! jd: --- Spherical-wave expansion and spherical-wave RI keywords ---
       ! JCW: swri_cheb_batchsize == 0 indicates that per-integration scheme
       !      default should be selected on initialisation of
       !      swri_resolution_of_identity
       pub_swri_cheb_batchsize = esdf_integer('swri_cheb_batchsize', 0)
       pub_swri_swop_smoothing = esdf_boolean('swri_swop_smoothing', .false.)
       pub_swri_overlap_indirect = esdf_boolean('swri_overlap_indirect',.false.)
       pub_swri_improve_inverse = esdf_boolean('swri_improve_inverse', .false.)
       pub_swri_print_eigenvalues = esdf_boolean('swri_print_eigenvalues', &
            .false.)
       pub_swri_assembly_prefix = esdf_string('swri_assembly_prefix', &
            pub_rootname, force_no_case_change = .true.)
       buf = esdf_string('swri_proximity_sort_point','0.0 0.0 0.0')
       ! jd: append "bohr" as default units
       buf = trim(adjustl(buf))//' bohr'
       read(buf,*) pub_swri_proximity_sort_point(:), tmp_unit
       pub_swri_proximity_sort_point(:) = pub_swri_proximity_sort_point(:) * &
            esdf_convfac(tmp_unit,'bohr')
       pub_swri_verbose = esdf_boolean('swri_verbose',.false.)
       pub_swx_c_threshold = esdf_double('swx_c_threshold',0.0_DP)
       pub_swx_dbl_grid = esdf_boolean('swx_dbl_grid',.false.)

       ! JCW: --- Spherical harmonic rotation ---
       pub_use_sph_harm_rot = esdf_boolean('use_sph_harm_rot',.false.)

       ! jd: --- HF exchange keywords ---
       pub_hfx_max_l = esdf_integer('hfx_max_l',-1)
       pub_hfx_max_q = esdf_integer('hfx_max_q',-1)
       pub_hfx_use_ri = esdf_string('hfx_use_ri', '<unset>',.true.)
       pub_hfx_cutoff = esdf_physical('hfx_cutoff',1000.0_DP,'bohr')
       pub_hfx_nlpp_for_exchange = esdf_boolean('hfx_nlpp_for_exchange',.false.)
       pub_hfx_debug = esdf_boolean('hfx_debug',.false.)
       pub_hfx_read_xmatrix = esdf_boolean('hfx_read_xmatrix', .false.)
       pub_hfx_write_xmatrix = esdf_boolean('hfx_write_xmatrix', .false.)
       hfx_metric_string = esdf_string('hfx_metric', 'electrostatic')
       pub_hfx_bessel_rad_nptsx = esdf_integer('hfx_bessel_rad_nptsx',100000)
       pub_hfx_memory_limit = esdf_integer('hfx_memory_limit',4096)
       buf                  = esdf_string('hfx_memory_weights','-1.0 -1.0 -1.0')
       read(buf,*) pub_hfx_memory_weights(:)
       pub_cache_limit_for_swops = esdf_integer(&
            'cache_limit_for_swops', 1024)                ! MiB/rank
       pub_cache_limit_for_swops2 = esdf_integer(&
            'cache_limit_for_swops2', 1024)               ! MiB/rank
       pub_cache_limit_for_expansions = esdf_integer(&
            'cache_limit_for_expansions', 1024)           ! MiB/rank
       pub_cache_limit_for_prods = esdf_integer(&
            'cache_limit_for_prods', 1024)                ! MiB/rank
       pub_cache_limit_for_dknblks = esdf_integer(&
            'cache_limit_for_dknblks', -1)                ! = use default

       pub_ht_stash_size = esdf_integer('ht_stash_size', 256) ! MiB/thread

       select case(trim(hfx_metric_string))
       case ('OVERLAP')
          pub_hfx_metric = METRIC_OVERLAP
       case ('ELECTROSTATIC')
          pub_hfx_metric = METRIC_ELECTROSTATIC
       case default
          call utils_abort('hfx_metric must be OVERLAP or ELECTROSTATIC.')
       end select
       ! jd: ------ End of SWX, SWRI and Hartree-Fock exchange keywords ------

       ! jd: Zero total force?
       pub_zero_total_force = esdf_boolean('zero_total_force',.true.)
       pub_mw_total_force = esdf_boolean('mw_total_force',.false.)

       ! jd: Are the ions going to move at any time?
       if(index(pub_all_tasks, 'GEOMETRYOPTIMIZATION') > 0 .or. &
            index(pub_all_tasks, 'TRANSITIONSTATESEARCH') > 0 .or. &
            index(pub_all_tasks, 'MOLECULARDYNAMICS') > 0 .or. &
            index(pub_all_tasks, 'PHONON') > 0 .or. &
            index(pub_all_tasks, 'FORCETEST') > 0) then
          pub_ions_will_move = .true.
       else
          pub_ions_will_move = .false.
       end if

       ! jd: Do we need to calculate forces at any time?
       if (pub_write_forces .or. pub_ions_will_move) then
          pub_forces_needed = .true.
       else
          pub_forces_needed = .false.
       end if

       ! jd: pub_psinc_spacing cannot be zero in any dimension, unless it has not
       !     been specified altogether. Negative values are not OK either.
       call utils_assert(all(pub_psinc_spacing == 0D0) .or. &
            .not. any(pub_psinc_spacing == 0D0), &
            'None of the components in psinc_spacing is allowed to be zero')
       call utils_assert(.not. any(pub_psinc_spacing < 0D0), &
            'None of the components in psinc_spacing is allowed to be negative')

       ! jd: Handle the case where KE cutoff energy has not been specified
       if (pub_cutoff_energy < 0D0) then
          if(any(pub_psinc_spacing > 0D0)) then
             ! jd: If psinc spacing has been specified, determine pub_cutoff_energy
             !     for the user. It will be used in the calculation of
             !     pub_ngwf_cg_max_step in a moment.
             pub_cutoff_energy = (((6.0_DP*(PI**2))**(2.0_DP/3.0_DP)) / 2.0_DP) / &
                  (maxval(pub_psinc_spacing)**2)
             write(stdout,'(/a,f10.5,a)') CRLF//&
                  'WARNING: "cutoff_energy" not specified, but "psinc_spacing" &
                  &specified.'//CRLF//'         Using the latter to set &
                  &the KE cutoff energy to ', pub_cutoff_energy*HARTREE_IN_eVs, &
                  'eV,'//CRLF//'         thus *overriding* the default.'//CRLF
          else
             ! jd: If psinc spacing hasn't been specified either, go on with the
             !     default KE cutoff.
             pub_cutoff_energy = -pub_cutoff_energy
             write(stdout,'(/a,f10.5,a)') CRLF//&
                  'WARNING: "cutoff_energy" not specified, using a &
                  &default of ', pub_cutoff_energy*HARTREE_IN_eVs, &
                  'eV.'//CRLF
          end if
       end if

       ! ndmh: do we need conduction NGWFs at any time?
       if ( (index(pub_all_tasks,'COND') > 0) .or. &
            (index(pub_all_tasks,'PROPERTIES_COND') > 0) .or. &
            (index(pub_all_tasks,'LR_TDDFT') > 0)) then

          pub_cond_calculate_any_task = .true.
       else
          pub_cond_calculate_any_task = .false.
       end if

       pub_do_bandstructure = esdf_block('bs_kpoint_path',nlines)
       if (pub_do_bandstructure) then
          bs_start = .true.
          pub_bs_kpoint_path_length = 0
          do iline=1,nlines
             if (index(block_data(iline),'BREAK') > 0) then
                bs_start = .true.
             else
                if (.not. bs_start) pub_bs_kpoint_path_length = &
                     pub_bs_kpoint_path_length + 1
                bs_start = .false.
             end if
          end do
          if (pub_bs_kpoint_path_length < 1) then
             write(stdout,'(/a/23x,a)') 'WARNING in get_rundat: &
                  &invalid path in bs_kpoint_path block', &
                  'bandstructure calculation will be skipped'
             pub_do_bandstructure = .false.
          end if
          if (pub_bs_method /= 'TB' .and. pub_bs_method /= 'KP') then
             write(stdout,'(/2a/23x,a)') 'WARNING in get_rundat: &
                  &unknown bs_method: ',pub_bs_method, &
                  'bandstructure calculation will be skipped'
             pub_do_bandstructure = .false.
          end if
       end if

       ! gcc32: band-structure unfolding
       pub_bsunfld_calculate = esdf_block('bsunfld_kpoint_path',nlines)
       if (pub_bsunfld_calculate) then
          bsunfld_start = .true.
          pub_bsunfld_kpoint_path_length = 0
          do iline=1,nlines
             if (index(block_data(iline),'BREAK') > 0) then
                bsunfld_start = .true.
             else
                if (.not. bsunfld_start) pub_bsunfld_kpoint_path_length = &
                     pub_bsunfld_kpoint_path_length + 1
                bsunfld_start = .false.
             end if
          end do
          if (pub_bsunfld_kpoint_path_length < 1) then
             write(stdout,'(/a/23x,a)') 'WARNING in get_rundat: &
                  &invalid path in bs_kpoint_path block', &
                  'bandstructure calculation will be skipped'
             pub_bsunfld_calculate = .false.
          end if
       end if



       ! agreco: read number of k-points for BZ sampling from kpoint_list block
       pub_kpoints_specified_list = esdf_block('kpoint_list',nlines)
       if (pub_kpoints_specified_list) then
          ! agreco: abort if k-points are specified both as a list and
          ! as a mesh
          if (any(pub_kp_grid_size /= 1).or.any(pub_kp_grid_shift /= 0)) then
             call utils_abort('K-point mesh and list bost specified in the &
                &input file: please choose one of the two')
          end if

          pub_num_kpoints_temp = nlines
          if (pub_num_kpoints_temp < 1) then
             write(stdout,'(/a/23x,a)') 'WARNING in get_rundat: &
                  &invalid list in kpoint_list block', &
                  'Gamma point calculation will be performed'
             pub_kpoints_specified_list = .false.
             pub_num_kpoints_temp = 1
             pub_gamma_point_only = .true.
          else
             pub_gamma_point_only = .false.
          end if
       else
          pub_num_kpoints_temp = 1
          pub_gamma_point_only = .true.
       end if

       ! lpl: DDEC parameters
       pub_ddec = esdf_boolean('ddec_calculate',.false.)
       pub_ddec_write = esdf_boolean('ddec_write_rad',.false.)
       pub_ddec_hirshfeld = esdf_boolean('ddec_classical_hirshfeld',.false.)

       pub_ddec_rad_npts = esdf_integer('ddec_rad_npts',100)
       pub_ddec_rad_rcut = esdf_physical('ddec_rad_rcut',5.0_DP*ANGSTROM,'bohr')

       pub_ddec_maxit = esdf_integer('ddec_maxit',2000)
       pub_ddec_core_maxit = esdf_integer('ddec_core_maxit',2000)
       pub_ddec_conv_thresh = esdf_double('ddec_conv_threshold',1e-5_DP)
       pub_ddec_IH_frac = esdf_double('ddec_IH_fraction',(3.0_DP/14.0_DP))
       pub_ddec_zero_thresh = esdf_double('ddec_zero_threshold',1e-10_DP)
       pub_ddec_use_coredens = esdf_boolean('ddec_use_coredens',.true.) ! Changed from F

       pub_ddec_core_correction = esdf_boolean('ddec_core_correction',.true.) ! Changed from F
       pub_ddec_reshape_dens = esdf_boolean('ddec_reshape_dens',.true.) ! Changed from F
       pub_ddec_core_corr_maxit = esdf_integer('ddec_core_corr_maxit',40)

       pub_ddec_ionic_range = esdf_integer('ddec_IH_ionic_range',2)
! lpl: Refer to COMMENT#01 in ddec_mod
!       pub_ddec_atomsolve_maxit = esdf_integer('ddec_atomsolve_maxit',480)
!       pub_ddec_rcomp_maxit = esdf_integer('ddec_rcomp_maxit',50)
!       pub_ddec_target_radius = esdf_physical('ddec_target_radius',11.5_DP,'bohr')
!       pub_ddec_rcomp_econv = esdf_physical('ddec_rcomp_econv',5e-8_DP,'hartree')
!       pub_ddec_write_rcomp = esdf_boolean('ddec_write_rcomp',.false.)

       pub_ddec_refdens_path = esdf_string('ddec_refdens_path','',.true.)
       pub_ddec_refdens_init = esdf_boolean('ddec_refdens_init',.true.)
       pub_ddec_renormalize_refdens = esdf_boolean('ddec_renormalize_refdens',.false.)
       pub_ddec_c3_refdens = esdf_boolean('ddec_c3_refdens',.true.) ! Changed from F
       pub_ddec_rad_shell_mode = esdf_string('ddec_rad_shell_mode','MIDDLE')

       pub_ddec_multipole = esdf_boolean('ddec_multipole',.false.)
! lpl: Refer to COMMENT#07 in ddec_mod
       pub_ddec_moment = esdf_integer('ddec_moment',-1)
!       pub_ddec_moment_order = esdf_integer('ddec_moment_order',1)

       pub_ddec_interp_rad_dens = esdf_boolean('ddec_interp_rad_dens',.false.)
       pub_ddec_avg_rad = esdf_boolean('ddec_avg_rad',.false.)
       pub_ddec_min_shell_dens = esdf_double('ddec_min_shell_dens',100.0_DP)
       pub_ddec_ref_shell_mode = esdf_integer('ddec_ref_shell_mode',0)

       pub_ddec_format_dens = esdf_boolean('ddec_format_dens',.false.)

       pub_ddec_eff_decay_exp = esdf_boolean('ddec_eff_decay_exp',.true.) ! Changed from F
       pub_ddec_eff_decay_rmin = esdf_physical('ddec_eff_decay_rmin',2.0_DP*ANGSTROM,'bohr')

! lpl: Refer to COMMENT#06 in ddec_mod
!       ! lpl: This is related ro the still-dodgy subroutine 'ddec_rmse'
!       pub_ddec_rmse = esdf_boolean('ddec_rmse',.false.)
       ! lpl: END DDEC parameters

       ! lpl: NBO stuff
       ! lpl: 26092011 - nbo_init_lclowdin now defaults to TRUE
       pub_write_nbo = esdf_boolean('write_nbo',.false.)
       pub_nbo_init_lclowdin  = esdf_boolean('nbo_init_lclowdin',.true.)
       pub_nbo_write_lclowdin = esdf_boolean('nbo_write_lclowdin',.false.)
       pub_nbo_write_npacomp  = esdf_boolean('nbo_write_npacomp',.false.)
       pub_nbo_scale_dm       = esdf_boolean('nbo_scale_dm',.true.)
       pub_nbo_write_dipole   = esdf_boolean('nbo_write_dipole',.false.)
       pub_nbo_scale_spin     = esdf_boolean('nbo_scale_spin',.true.)

       pub_nbo_pnao_analysis  = esdf_boolean('nbo_pnao_analysis',.false.)

       pub_nbo_aopnao_scheme = esdf_string('nbo_aopnao_scheme','ORIGINAL')
       if(pub_nbo_aopnao_scheme /= 'ORIGINAL' .and. &
            pub_nbo_aopnao_scheme /= 'DIAGONALIZATION' .and. &
            pub_nbo_aopnao_scheme /= 'NONE') then
          call utils_abort('Invalid nbo_aopnao_scheme specified ('//&
               trim(pub_nbo_aopnao_scheme)//')')
       end if

       pub_plot_nbo   = esdf_boolean('plot_nbo',.false.)
       if(pub_write_nbo .and. pub_plot_nbo) write(stdout,'(a)') 'ERROR: &
            &pub_write_nbo and pub_plot_nbo cannot be simulataneously T'
       buf                 = esdf_string('nbo_plot_orbtype','')
       pub_nbo_plotorbtype = buf(1:8)
       ! lpl: NBO stuff

       ! aeaa: DDEC anisotropy off center charges
       pub_ddec_aniso = esdf_boolean('ddec_aniso',.false.)
       pub_ddec_aniso_max_dis = esdf_double('ddec_aniso_max_dis', &
            0.8_DP)
       pub_ddec_aniso_max_dis_halogen = esdf_double( &
            'ddec_aniso_max_dis_halogen', 1.5_DP)
       pub_ddec_aniso_error_thres = esdf_double( &
            'ddec_aniso_error_thres', 0.9025_DP)
       pub_ddec_aniso_error_reduce = esdf_double( &
            'ddec_aniso_error_reduce', 0.0625_DP)

       ! mjsp: EDA
       pub_eda_continuation = esdf_boolean('eda_continuation',.false.)
       ! mjsp: correct pub_eda_continuation contradictions
       if (pub_eda_continuation .and. .not.(pub_write_tightbox_ngwfs)) then
          write(stdout,'(/a)') 'WARNING in get_rundat: &
               &eda_continuation is true, and so write_tightbox_ngwfs will be'//CRLF//'&
               &overridden to true.'
          pub_write_tightbox_ngwfs = .true.
       end if
       if (pub_eda_continuation .and. .not.(pub_write_denskern)) then
          write(stdout,'(/a)') 'WARNING in get_rundat: &
               &eda_continuation is true, and so write_denskern will be'//CRLF//'&
               &overridden to true.'
          pub_write_denskern = .true.
       end if
       pub_eda_continuation_loaded_dkn = .false.
       pub_eda_continuation_loaded_ngwfs = .false.
       pub_eda_write = esdf_boolean('eda_write',.false.)
       pub_eda_read_frags = esdf_boolean('eda_read_frags',.false.)
       pub_eda_read_super = esdf_boolean('eda_read_super',.false.)
       pub_eda_mode = -1
       pub_frag_counter = 0
       pub_frag_counter2 = 0
       pub_eda_have_sp_frags = .false.
       pub_eda_deltadens = esdf_boolean('eda_deltadens',.false.)
!       pub_eda_scfmi_restrict_ngwfs = esdf_boolean('scfmi_restrict_ngwfs',.false.)
       pub_eda_frag_isol_pol   = esdf_boolean('eda_frag_isol_pol',.false.)
       pub_eda_frag_isol_ct    = esdf_boolean('eda_frag_isol_ct',.false.)
       pub_eda_reset_ngwfs_pol = esdf_boolean('eda_reset_ngwfs_pol',.false.)
       pub_eda_reset_ngwfs_ct  = esdf_boolean('eda_reset_ngwfs_ct',.true.)
       pub_eda_split_atoms     = esdf_boolean('eda_split_atoms',.false.)
       pub_eda_nodiag = esdf_integer('eda_nodiag',3)
       pub_eda_scfmi     = .false.
       pub_eda_scfmi_any = .false.
       pub_eda_frz_x = 0_DP  ! Intra+interfragmental exchange (FRZ)
       pub_eda_frz_c = 0_DP  ! Intra+interfragmental correlation (FRZ)
       pub_eda_rep_x = 0_DP  ! Intra+interfragmental exchange (REP)
       pub_eda_rep_c = 0_DP  ! Intra+interfragmental correlation (REP)
       pub_eda_intrafrag_x_frz = 0_DP  ! Intrafragmental dE exchange contrib (isol-frz)
       pub_eda_intrafrag_c_frz = 0_DP  ! Intrafragmental dE correlation contrib (isol-frz)
       pub_eda_dE_x_rep = 0_DP  ! Intra+interfragmental dE exchange contrib (rep-frz)
       pub_eda_dE_c_rep = 0_DP  ! Intra+interfragmental dE correlation contrib (rep-frz)

       ! az, jd: --- DMA related variables ---
       pub_dma_calculate = esdf_boolean('dma_calculate',.false.)
       pub_dma_output_potential = esdf_boolean('dma_output_potential',.false.)
       pub_dma_output_potential_reference = esdf_boolean(&
            'dma_output_potential_reference',.false.)
       pub_dma_max_l = esdf_integer('dma_max_l',-1)
       pub_dma_max_q = esdf_integer('dma_max_q',-1)
       pub_dma_scale_charge = esdf_boolean('dma_scale_charge',.true.)
       pub_dma_precise_gdma_output = esdf_boolean('dma_precise_gdma_output', &
            .true.)
       pub_dma_target_num_val_elec = esdf_double('dma_target_num_val_elec', &
            -999999.0_DP)
       pub_dma_bessel_averaging = esdf_boolean('dma_bessel_averaging',.true.)
       pub_dma_multipole_scaling = esdf_double('dma_multipole_scaling', 1.0_DP)
       pub_dma_dipole_scaling = esdf_double('dma_dipole_scaling', 1.0_DP)
       pub_dma_quadrupole_scaling = esdf_double('dma_quadrupole_scaling', 1.0_DP)
       dma_metric_string = esdf_string('dma_metric', 'electrostatic')
       pub_dma_use_ri = esdf_string('dma_use_ri', '<unset>',.true.)

       ! jcap: Parameters for embedding mean field theory
       pub_active_xc_functional = esdf_string('active_xc_functional',pub_xc_functional)
       pub_emft = esdf_boolean('use_emft',.false.)
       pub_emft_follow = esdf_boolean('use_emft_follow',.false.)
       pub_emft_lnv_only = esdf_boolean('use_emft_lnv_only',.false.)
       pub_active_region = esdf_integer('active_region',1)
       pub_block_orthogonalise = esdf_boolean('block_orthogonalise',.false.)
       pub_emft_lnv_steps = esdf_integer('emft_lnv_steps',10)

       select case(trim(dma_metric_string))
       case ('OVERLAP')
          pub_dma_metric = METRIC_OVERLAP
       case ('ELECTROSTATIC')
          pub_dma_metric = METRIC_ELECTROSTATIC
       case default
          call utils_abort('dma_metric must be OVERLAP or ELECTROSTATIC.')
       end select

       if(pub_dma_calculate) then
          if(trim(pub_dma_use_ri) == trim('<unset>')) then
             call utils_abort(&
               '''dma_use_ri'' is mandatory when using DMA.')
          end if
          if(pub_dma_max_l == -1) then
             call utils_abort(&
               '''dma_max_l'' is mandatory when using DMA.')
          end if
          if(pub_dma_max_q == -1) then
             call utils_abort(&
               '''dma_max_q'' is mandatory when using DMA.')
          end if

          call utils_assert(pub_dma_max_q > 0,'''dma_max_q'' must be positive')
          call utils_assert(pub_dma_max_l >= 0, &
               '''dma_max_l'' must be non-negative')

          if(pub_dma_bessel_averaging) then
             call utils_assert(pub_dma_max_q > 1, &
                  '''dma_max_q 1'' is incompatible with Bessel averaging.')
          end if

          call utils_assert(pub_dma_multipole_scaling > SAFE_DIV_EPS, &
               '''dma_multipole_scaling'' cannot be &
               &negative, zero or very close to zero')
               ! because it is used in divisions
       end if

       ! jd: --- end of DMA related variables ---

       pub_turn_off_ewald  = esdf_boolean('turn_off_ewald',.false.)

       ! gom: Omit Hartree terms in energy and potential
       pub_turn_off_hartree  = esdf_boolean('turn_off_hartree',.false.)

       ! gom: By default we do not perform DFT+nu
       pub_dft_nu = .false.

       ! gom: Combine all projectors into one Hubbard site
       pub_hubbard_unify_sites = esdf_boolean('hubbard_unify_sites',.false.)

       ! gom: Reset tensor_corr to 4 if unifying sites
       if (pub_hubbard_unify_sites .eqv. .true.) then
           pub_hub_tensor_corr = 4
       endif

       ! ddor: Decide whether to compute Hubbard metric tensor forces
       if ((pub_hub_tensor_corr == 4) .or. (pub_hub_tensor_corr == 5)) then
          pub_hub_tensor_forces = .true.
       else
          pub_hub_tensor_forces = .false.
       endif
       ! ddor: The input flag overrides this decision if present
       pub_hub_tensor_forces = esdf_boolean('hubbard_tensor_forces', &
            & pub_hub_tensor_forces)

       ! cks: Hyper Hartree-Fock parameters
       pub_hhf_nstates = esdf_integer('hhf_nstates', 0)

       ! rjc: Electron Localisation Descriptor and Kinetic Energy Density output:
       pub_eld_calculate = esdf_boolean('eld_calculate',.false.)
       pub_eld_function = esdf_string('eld_function','ELF')
       pub_ke_density_calculate = esdf_boolean('ke_density_calculate',.false.)

       ! jme: set output detail integers in root proc (avoids broadcasting the
       ! strings later, and pub_output_detail is used in the call below to
       ! util_default_threads.

       ! cks: set output detail integer parameter
       if (output_detail == 'BRIEF') then
          pub_output_detail = BRIEF
       else if (output_detail == 'NORMAL') then
          pub_output_detail = NORMAL
       else if (output_detail == 'VERBOSE') then
          pub_output_detail = VERBOSE
       end if

       ! ndmh: set PAW output detail integer parameter
       if (paw_output_detail == 'DEFAULT') then
          if (pub_output_detail==VERBOSE) then
             pub_paw_output_detail = NORMAL
          else if (pub_output_detail==NORMAL) then
             pub_paw_output_detail = BRIEF
          else if (pub_output_detail==BRIEF) then
             pub_paw_output_detail = BRIEF
          end if
       else if (paw_output_detail == 'BRIEF') then
          pub_paw_output_detail = BRIEF
       else if (paw_output_detail == 'NORMAL') then
          pub_paw_output_detail = NORMAL
       else if (paw_output_detail == 'VERBOSE') then
          pub_paw_output_detail = VERBOSE
       end if

       ! smmd: set MD output detail integer parameter
       if (md_output_detail == 'DEFAULT') then
          if (pub_output_detail==VERBOSE) then
             pub_md_output_detail = NORMAL
          else if (pub_output_detail==NORMAL) then
             pub_md_output_detail = BRIEF
          else if (pub_output_detail==BRIEF) then
             pub_md_output_detail = BRIEF
          end if
       else if (md_output_detail == 'BRIEF') then
          pub_md_output_detail = BRIEF
       else if (md_output_detail == 'NORMAL') then
          pub_md_output_detail = NORMAL
       else if (md_output_detail == 'VERBOSE') then
          pub_md_output_detail = VERBOSE
       end if

       ! jd: set SWX output detail integer parameter
       if (swx_output_detail == 'BRIEF') then
          pub_swx_output_detail = BRIEF
       else if (swx_output_detail == 'NORMAL') then
          pub_swx_output_detail = NORMAL
       else if (swx_output_detail == 'VERBOSE') then
          pub_swx_output_detail = VERBOSE
       else if (hfx_output_detail == 'DEFAULT') then
          pub_swx_output_detail = pub_output_detail
       end if

       ! jd: set HFx output detail integer parameter
       if (hfx_output_detail == 'BRIEF') then
          pub_hfx_output_detail = BRIEF
       else if (hfx_output_detail == 'NORMAL') then
          pub_hfx_output_detail = NORMAL
       else if (hfx_output_detail == 'VERBOSE') then
          pub_hfx_output_detail = VERBOSE
       else if (hfx_output_detail == 'PROLIX') then
          pub_hfx_output_detail = PROLIX
       else if (hfx_output_detail == 'MAXIMUM') then
          pub_hfx_output_detail = MAXIMUM
       else if (hfx_output_detail == 'DEFAULT') then
          pub_hfx_output_detail = pub_output_detail
       end if

       ! lk: set GEOM output detail integer parameter
       if (geom_output_detail == 'BRIEF') then
          pub_geom_output_detail = BRIEF
       else if (geom_output_detail == 'NORMAL') then
          pub_geom_output_detail = NORMAL
       else if (geom_output_detail == 'VERBOSE') then
          pub_geom_output_detail = VERBOSE
       else if (geom_output_detail == 'PROLIX') then
          pub_geom_output_detail = PROLIX
       else if (geom_output_detail == 'MAXIMUM') then
          pub_geom_output_detail = MAXIMUM
       else if (geom_output_detail == 'DEFAULT') then
          ! lk: the default behaviour is to follow
          !     the overall output detail.
          pub_geom_output_detail = pub_output_detail
       end if

       ! ndmh: remove quotes from paths
       call utils_strip_char_from_string(pub_pseudo_path,'"')
       call utils_strip_char_from_string(pub_pseudo_path,'''')
       call utils_strip_char_from_string(pub_ddec_refdens_path,'"')
       call utils_strip_char_from_string(pub_ddec_refdens_path,'''')

       ! rab207: enable kernel checks
       if (pub_kernel_check_all) then
          write(stdout,'(/a)') &
               'KERNEL_CHECK_ALL enabled. Modifying input parameters:'
          if (pub_maxit_kernel_occ_check == 0) then
             write(stdout,'(3x,a)')  'maxit_kernel_occ_check = 1'
             pub_maxit_kernel_occ_check = 1
          endif
          if (.not.pub_kernel_update) then
             write(stdout,'(3x,a)')  'kernel_update = T'
             pub_kernel_update = .true.
          endif
          if (.not.pub_kernel_track_mid_occ) then
             write(stdout,'(3x,a)')  'kernel_track_mid_occ = T'
             pub_kernel_track_mid_occ = .true.
          endif
          if (pub_kerfix < 2) then
             write(stdout,'(3x,a)')  'kerfix = 2'
             pub_kerfix = 2
          endif
       endif

       ! JCW: Additional exchange-correlation functional options
       ! JCW: * Minimum threshold value for tau (KE density)
       pub_xc_min_tau = esdf_double("xc_min_tau",1e-12_DP)
       ! JCW: * Alternative functional for initial guess (meta-GGAs only)
       ! JCW: 31/05/18:
       ! JCW: Default changed from "NULL" to "PBE", since using PBE for in
       ! JCW: xc_radial demonstrably improves NGWF convergence when used to
       ! JCW: produce the initial NGWFs, as compared to omitting the XC
       ! JCW: contribution.
       ! JCW:
       ! JCW: Setting the value manually to "NULL" will no longer work ---
       ! JCW: to omit the XC contribution from the initial guess, use the
       ! JCW: string "NONE".
       pub_xc_initial_functional= esdf_string('xc_initial_functional','PBE')


       ! ab: analytical overlap integrals in DFTB
       if (pub_dftb) then
          pub_dftb_overlap_analytical= esdf_boolean('dftb_overlap_analytical', &
          .true.)
       else
          pub_dftb_overlap_analytical= .false.
       end if

       ! ab: use cartesian NGWFs for analytical overlap integrals
       pub_dftb_cartesian_ngwfs = pub_dftb_overlap_analytical

       ! jd: DFTB param filename. Defaults defined in constants_mod.
       !     Third param prevents case conversion.
       pub_dftb_method_param_file = esdf_string('dftb_method_param_file', &
            DFTB_DEFAULT_METHOD_PARAM_FILES(pub_dftb_method),.true.)
       pub_dftb_common_param_file = esdf_string('dftb_common_param_file', &
            DFTB_DEFAULT_COMMON_PARAM_FILE,.true.)
       pub_dftb_coord_cutoff = esdf_physical('dftb_coord_cutoff', 40.0_DP, 'bohr')
       pub_dftb_rep_cutoff = esdf_physical('dftb_rep_cutoff', 40.0_DP, 'bohr')
       pub_dftb_srb_cutoff = esdf_physical('dftb_srb_cutoff', sqrt(200.0_DP), 'bohr')
       pub_dftb_overlap_cutoff = esdf_physical('dftb_overlap_cutoff', sqrt(2000.0_DP), 'bohr')
       pub_dftb_ewald_replicate_xtb = esdf_boolean('dftb_ewald_replicate_xtb', .false.)
       pub_dftb_ewald_parameter = esdf_physical('dftb_ewald_parameter', -1.0_DP, '1/ang')
       pub_dftb_calculate_force = pub_dftb .and. pub_forces_needed

    end if ! pub_on_root


    ! gcc32: ROOT finished reading the input bandstructure path
    ! gcc32: now, each proc must allocate pub_bs_kpoint_path_start
    ! gcc32: and pub_bs_kpoint_path_end
    call comms_bcast(pub_root_proc_id, pub_do_bandstructure)
    call comms_bcast(pub_root_proc_id, pub_bs_kpoint_path_length)

    ! gcc32: do the same for bandstructure-unfolding
    call comms_bcast(pub_root_proc_id, pub_bsunfld_calculate)
    call comms_bcast(pub_root_proc_id, pub_bsunfld_kpoint_path_length)

    if (pub_do_bandstructure) then
       allocate(pub_bs_kpoint_path_start(3,pub_bs_kpoint_path_length), &
            stat=ierr)
       call utils_alloc_check('get_rundat','pub_bs_kpoint_path_start',ierr)
       allocate(pub_bs_kpoint_path_end(3,pub_bs_kpoint_path_length), &
            stat=ierr)
       call utils_alloc_check('get_rundat','pub_bs_kpoint_path_end',ierr)

       if (pub_on_root) then
          pub_do_bandstructure = esdf_block('bs_kpoint_path',nlines)
          bs_start = .true.
          ipath = 0
          do iline=1,nlines
             if (index(block_data(iline),'BREAK') > 0) then
                bs_start = .true.
             else
                read(block_data(iline),*) dummy_real
                if (.not. bs_start) then
                   ipath = ipath + 1
                   pub_bs_kpoint_path_end(:,ipath) = dummy_real
                else
                   bs_start = .false.
                end if
                if (ipath < pub_bs_kpoint_path_length) &
                     pub_bs_kpoint_path_start(:,ipath+1) = dummy_real
             end if
          end do
       end if

       call comms_bcast(pub_root_proc_id,pub_bs_kpoint_path_start)
       call comms_bcast(pub_root_proc_id,pub_bs_kpoint_path_end)
    end if

    ! gcc32: similarly for band-structure unfolding
    if (pub_bsunfld_calculate) then
       allocate(pub_bsunfld_kpoint_path_start(3, &
            pub_bsunfld_kpoint_path_length), stat=ierr)
       call utils_alloc_check('get_rundat','pub_bsunfld_kpoint_path_start',ierr)
       allocate(pub_bsunfld_kpoint_path_end(3,pub_bsunfld_kpoint_path_length), &
            stat=ierr)
       call utils_alloc_check('get_rundat','pub_bsunfld_kpoint_path_end',ierr)

       if (pub_on_root) then
          pub_bsunfld_calculate = esdf_block('bsunfld_kpoint_path',nlines)
          bsunfld_start = .true.
          ipath = 0
          do iline=1,nlines
             if (index(block_data(iline),'BREAK') > 0) then
                bsunfld_start = .true.
             else
                read(block_data(iline),*) dummy_real
                if (.not. bsunfld_start) then
                   ipath = ipath + 1
                   pub_bsunfld_kpoint_path_end(:,ipath) = dummy_real
                else
                   bsunfld_start = .false.
                end if
                if (ipath < pub_bsunfld_kpoint_path_length) &
                     pub_bsunfld_kpoint_path_start(:,ipath+1) = dummy_real
             end if
          end do
       end if

       call comms_bcast(pub_root_proc_id,pub_bsunfld_kpoint_path_start)
       call comms_bcast(pub_root_proc_id,pub_bsunfld_kpoint_path_end)
    end if


    ! agreco: let the other procs know about k-point list
    call comms_bcast(pub_root_proc_id, pub_kpoints_specified_list)
    call comms_bcast(pub_root_proc_id, pub_num_kpoints_temp)
    call comms_bcast(pub_root_proc_id, pub_gamma_point_only)

    if (pub_kpoints_specified_list) then
       allocate(pub_kpoint_list(4,pub_num_kpoints_temp), &
                stat=ierr)
       call utils_alloc_check('get_rundat','pub_kpoint_list',ierr)

       if (pub_on_root) then
          pub_kpoints_specified_list = esdf_block('kpoint_list',nlines)
          do iline=1,nlines
             ! agreco: append default weight for this k-point
             buf = trim(adjustl(block_data(iline)))//' 1.0'
             read(buf,*) dummy_real, kpoint_weight
             pub_kpoint_list(1:3, iline) = dummy_real
             pub_kpoint_list(4, iline) = kpoint_weight
          end do
       end if

       call comms_bcast(pub_root_proc_id,pub_kpoint_list)

    end if

    !...and broadcasts to all procs
    call comms_bcast(pub_root_proc_id, pub_check_stack_size)
    call comms_bcast(pub_root_proc_id, pub_all_tasks)
    call comms_bcast(pub_root_proc_id, pub_cutoff_energy)
    call comms_bcast(pub_root_proc_id, pub_kernel_cutoff)
    call comms_bcast(pub_root_proc_id, pub_xc_functional)
    call comms_bcast(pub_root_proc_id, pub_libxc_x_func_id)
    call comms_bcast(pub_root_proc_id, pub_libxc_c_func_id)
    call comms_bcast(pub_root_proc_id, pub_charge)
    call comms_bcast(pub_root_proc_id, pub_real_spin)
    call comms_bcast(pub_root_proc_id, pub_spin)
    call comms_bcast(pub_root_proc_id, pub_spin_polarised)
    call comms_bcast(pub_root_proc_id, pub_confined_ngwfs)
    call comms_bcast(pub_root_proc_id, pub_confined_ngwfs_barrier)
    call comms_bcast(pub_root_proc_id, pub_maxit_ngwf_cg_confined)
    call comms_bcast(pub_root_proc_id, pub_constant_efield)
    call comms_bcast(pub_root_proc_id, pub_fftbox_pref)
    call comms_bcast(pub_root_proc_id, pub_augbox_pref)
    call comms_bcast(pub_root_proc_id, pub_ppd_npoints)
    call comms_bcast(pub_root_proc_id, pub_extend_ngwf)
    call comms_bcast(pub_root_proc_id, pub_full_rand_ngwf)
    call comms_bcast(pub_root_proc_id, pub_rand_seed_ngwf)
    call comms_bcast(pub_root_proc_id, pub_rand_sigma)
    call comms_bcast(pub_root_proc_id, pub_rand_conv)
    call comms_bcast(pub_root_proc_id, pub_conv_func)
    call comms_bcast(pub_root_proc_id, pub_conv_region_width)
    call comms_bcast(pub_root_proc_id, pub_ngwf_plot_every_it)
    call comms_bcast(pub_root_proc_id, pub_cmplx_ngwfs)
    call comms_bcast(pub_root_proc_id, pub_read_real_ngwfs)
    call comms_bcast(pub_root_proc_id, pub_show_overlap)
    call comms_bcast(pub_root_proc_id, pub_efield_origin)
    call comms_bcast(pub_root_proc_id, pub_single_kpt)
    call comms_bcast(pub_root_proc_id, pub_kpoint_method)
    call comms_bcast(pub_root_proc_id, pub_imag_thr)
    call comms_bcast(pub_root_proc_id, pub_ngwfs_phase)
    call comms_bcast(pub_root_proc_id, pub_ngwfs_rand_phase)
    call comms_bcast(pub_root_proc_id, pub_check_hermitian)
    call comms_bcast(pub_root_proc_id, pub_check_density)
    call comms_bcast(pub_root_proc_id, pub_kp_grid_size)
    call comms_bcast(pub_root_proc_id, pub_kp_grid_shift)
    call comms_bcast(pub_root_proc_id, pub_contracoham_radmult)
    call comms_bcast(pub_root_proc_id, pub_psinc_spacing)
    call comms_bcast(pub_root_proc_id, pub_devel_code)
    call comms_bcast(pub_root_proc_id, pub_input_xyz_file)
    call comms_bcast(pub_root_proc_id, pub_nnho)
    call comms_bcast(pub_root_proc_id, pub_ngwf_halo)
    call comms_bcast(pub_root_proc_id, pub_nonsc_forces)
    call comms_bcast(pub_root_proc_id, pub_external_pressure)
    call comms_bcast(pub_root_proc_id, pub_isosurface_cutoff)
    call comms_bcast(pub_root_proc_id, pub_smoothing_factor)
    call comms_bcast(pub_root_proc_id, pub_maxit_palser_mano)
    call comms_bcast(pub_root_proc_id, pub_maxit_kernel_occ_check)
    call comms_bcast(pub_root_proc_id, pub_maxit_pen)
    call comms_bcast(pub_root_proc_id, pub_pen_param)
    call comms_bcast(pub_root_proc_id, minit_lnv)
    call comms_bcast(pub_root_proc_id, maxit_lnv)
    call comms_bcast(pub_root_proc_id, lnv_threshold_orig)
    call comms_bcast(pub_root_proc_id, pub_lnv_cg_type)
    call comms_bcast(pub_root_proc_id, pub_lnv_cg_max_step)
    call comms_bcast(pub_root_proc_id, pub_exact_lnv)
    call comms_bcast(pub_root_proc_id, pub_old_lnv)
    call comms_bcast(pub_root_proc_id, pub_lnv_check_trial_steps)
    call comms_bcast(pub_root_proc_id, pub_kerfix)
    call comms_bcast(pub_root_proc_id, pub_maxit_kernel_fix)
    call comms_bcast(pub_root_proc_id, pub_initial_dens_realspace)
    call comms_bcast(pub_root_proc_id, pub_kernel_diis_size)
    call comms_bcast(pub_root_proc_id, pub_kernel_diis_maxit)
    call comms_bcast(pub_root_proc_id, pub_kernel_diis_threshold)
    call comms_bcast(pub_root_proc_id, pub_kernel_diis_linear_iter)
    call comms_bcast(pub_root_proc_id, pub_kernel_diis_coeff)
    call comms_bcast(pub_root_proc_id, pub_kernel_diis_conv_criteria)
    call comms_bcast(pub_root_proc_id, pub_kernel_diis_lshift)
    call comms_bcast(pub_root_proc_id, pub_kernel_diis_ls_iter)
    call comms_bcast(pub_root_proc_id, pub_kernel_diis_scheme)
    call comms_bcast(pub_root_proc_id, pub_foe)
    call comms_bcast(pub_root_proc_id, pub_dense_foe)
    call comms_bcast(pub_root_proc_id, pub_H2denskern_sparsity)
    call comms_bcast(pub_root_proc_id, pub_foe_mu_tol)
    call comms_bcast(pub_root_proc_id, pub_foe_test_sparsity)
    call comms_bcast(pub_root_proc_id, pub_foe_cheby_thres)
    call comms_bcast(pub_root_proc_id, pub_foe_avoid_inversions)
    call comms_bcast(pub_root_proc_id, pub_foe_check_entropy)
    call comms_bcast(pub_root_proc_id, pub_edft)
    call comms_bcast(pub_root_proc_id, pub_edft_maxit)
    call comms_bcast(pub_root_proc_id, pub_edft_max_step)
    call comms_bcast(pub_root_proc_id, pub_edft_smearing_width)
    call comms_bcast(pub_root_proc_id, pub_edft_free_energy_thres)
    call comms_bcast(pub_root_proc_id, pub_edft_energy_thres)
    call comms_bcast(pub_root_proc_id, pub_edft_entropy_thres)
    call comms_bcast(pub_root_proc_id, pub_edft_rms_gradient_thres)
    call comms_bcast(pub_root_proc_id, pub_edft_commutator_thres)
    call comms_bcast(pub_root_proc_id, pub_edft_fermi_thres)
    call comms_bcast(pub_root_proc_id, pub_edft_trial_step)
    call comms_bcast(pub_root_proc_id, pub_edft_nelec_thres)
    call comms_bcast(pub_root_proc_id, pub_edft_write_occ)
    call comms_bcast(pub_root_proc_id, pub_edft_extra_bands)
    call comms_bcast(pub_root_proc_id, pub_edft_init_maxit)
    call comms_bcast(pub_root_proc_id, pub_edft_spin_fix)
    call comms_bcast(pub_root_proc_id, pub_edft_spin_fix_orig)
    call comms_bcast(pub_root_proc_id, pub_edft_update_scheme)
    call comms_bcast(pub_root_proc_id, pub_edft_ham_diis_size)
    call comms_bcast(pub_root_proc_id, pub_eigensolver_orfac)
    call comms_bcast(pub_root_proc_id, pub_eigensolver_abstol)
    call comms_bcast(pub_root_proc_id, pub_write_hamiltonian)
    call comms_bcast(pub_root_proc_id, pub_read_hamiltonian)
    call comms_bcast(pub_root_proc_id, pub_write_overlap)
    call comms_bcast(pub_root_proc_id, pub_ngwf_cg_rotate)
    call comms_bcast(pub_root_proc_id, pub_edft_round_evals)
    call comms_bcast(pub_root_proc_id, pub_edft_grand_canonical)
    call comms_bcast(pub_root_proc_id, pub_edft_electrode_potential)
    call comms_bcast(pub_root_proc_id, pub_edft_reference_potential)
    call comms_bcast(pub_root_proc_id, pub_pbc_smeared_ion_rep)
    call comms_bcast(pub_root_proc_id, pub_chemical_softness)

    call comms_bcast(pub_root_proc_id, pub_energy_components_interval)
    call comms_bcast(pub_root_proc_id, pub_maxit_ngwf_cg)
    call comms_bcast(pub_root_proc_id, pub_freeze_switch_steps)
    call comms_bcast(pub_root_proc_id, pub_do_fandt)
    call comms_bcast(pub_root_proc_id, pub_freeze_envir_ngwfs)
    call comms_bcast(pub_root_proc_id, pub_embed_debug)
    call comms_bcast(pub_root_proc_id, pub_quit_region)
    call comms_bcast(pub_root_proc_id, pub_project_embed)
    call comms_bcast(pub_root_proc_id, pub_ngwf_cg_type)
    call comms_bcast(pub_root_proc_id, pub_ngwf_cg_max_step)
    call comms_bcast(pub_root_proc_id, pub_elec_cg_max)
    call comms_bcast(pub_root_proc_id, pub_precond_scheme)
    call comms_bcast(pub_root_proc_id, pub_k_zero)
    call comms_bcast(pub_root_proc_id, pub_r_precond)
    call comms_bcast(pub_root_proc_id, pub_precond_recip)
    call comms_bcast(pub_root_proc_id, pub_precond_real)
    call comms_bcast(pub_root_proc_id, pub_smooth_scheme)
    call comms_bcast(pub_root_proc_id, pub_r_smooth)
    call comms_bcast(pub_root_proc_id, pub_k_smooth)
    call comms_bcast(pub_root_proc_id, pub_occ_mix)
    call comms_bcast(pub_root_proc_id, pub_kernel_update)
    call comms_bcast(pub_root_proc_id, pub_kernel_christoffel_update)
    call comms_bcast(pub_root_proc_id, pub_kernel_track_mid_occ)
    call comms_bcast(pub_root_proc_id, pub_maxit_hotelling)
    call comms_bcast(pub_root_proc_id, pub_max_resid_hotelling)
    call comms_bcast(pub_root_proc_id, pub_rms_kernel_measure)
    call comms_bcast(pub_root_proc_id, ngwf_threshold_orig)
    call comms_bcast(pub_root_proc_id, pub_elec_energy_tol)
    call comms_bcast(pub_root_proc_id, pub_elec_force_tol)
    call comms_bcast(pub_root_proc_id, pub_ngwf_max_grad)
    call comms_bcast(pub_root_proc_id, pub_delta_e_conv)
    call comms_bcast(pub_root_proc_id, pub_kernel_force_conv)
    call comms_bcast(pub_root_proc_id, pub_dense_threshold)
    call comms_bcast(pub_root_proc_id, pub_comms_group_size)
    call comms_bcast(pub_root_proc_id, pub_ovlp_for_nonlocal)
    call comms_bcast(pub_root_proc_id, pub_use_space_filling_curve)
    call comms_bcast(pub_root_proc_id, pub_parallel_scheme)
    call comms_bcast(pub_root_proc_id, pub_coreham_denskern_guess)
    call comms_bcast(pub_root_proc_id, pub_check_atoms)
    call comms_bcast(pub_root_proc_id, pub_locpot_scheme)
    call comms_bcast(pub_root_proc_id, pub_smooth_projectors)
    call comms_bcast(pub_root_proc_id, pub_smooth_loc_pspot)
    call comms_bcast(pub_root_proc_id, pub_odd_psinc_grid)
    call comms_bcast(pub_root_proc_id, pub_even_psinc_grid)
    call comms_bcast(pub_root_proc_id, pub_realspace_projectors)
    call comms_bcast(pub_root_proc_id, pub_projectors_precalculate)
    call comms_bcast(pub_root_proc_id, pub_num_images)
    call comms_bcast(pub_root_proc_id, pub_image_sizes)
    call comms_bcast(pub_root_proc_id, pub_output_detail)
    call comms_bcast(pub_root_proc_id, pub_paw_output_detail)
    call comms_bcast(pub_root_proc_id, pub_md_output_detail)
    call comms_bcast(pub_root_proc_id, pub_timings_level)
    call comms_bcast(pub_root_proc_id, pub_timings_order)
    call comms_bcast(pub_root_proc_id, pub_max_runtime)
    call comms_bcast(pub_root_proc_id, pub_write_params)
    call comms_bcast(pub_root_proc_id, pub_esdf_dump)
    call comms_bcast(pub_root_proc_id, pub_write_forces)
    call comms_bcast(pub_root_proc_id, pub_write_positions)
    call comms_bcast(pub_root_proc_id, pub_write_velocities)
    call comms_bcast(pub_root_proc_id, pub_write_xyz)
    call comms_bcast(pub_root_proc_id, pub_print_qc)
    call comms_bcast(pub_root_proc_id, pub_cube_format)
    call comms_bcast(pub_root_proc_id, pub_dx_format)
    call comms_bcast(pub_root_proc_id, pub_dx_coarse)
    call comms_bcast(pub_root_proc_id, pub_dx_sig_digits)
    call comms_bcast(pub_root_proc_id, pub_grd_format)
    call comms_bcast(pub_root_proc_id, pub_write_density_plot)
    call comms_bcast(pub_root_proc_id, pub_write_polarisation_plot)
    call comms_bcast(pub_root_proc_id, pub_write_ngwf_plot)
    call comms_bcast(pub_root_proc_id, pub_write_ngwf_grad_plot)
    call comms_bcast(pub_root_proc_id, pub_write_ngwf_radial)
    call comms_bcast(pub_root_proc_id, pub_write_ngwf_grad_radial)
    call comms_bcast(pub_root_proc_id, pub_write_radial_step)
    call comms_bcast(pub_root_proc_id, pub_write_radial_smear)
    call comms_bcast(pub_root_proc_id, pub_write_initial_radial_ngwfs)
    call comms_bcast(pub_root_proc_id, pub_read_denskern)
    call comms_bcast(pub_root_proc_id, pub_write_denskern)
    call comms_bcast(pub_root_proc_id, pub_read_sub_denskern)
    call comms_bcast(pub_root_proc_id, pub_read_tightbox_ngwfs)
    call comms_bcast(pub_root_proc_id, pub_write_tightbox_ngwfs)
    call comms_bcast(pub_root_proc_id, pub_write_converged_dk_ngwfs)
    call comms_bcast(pub_root_proc_id, pub_read_sw_ngwfs)
    call comms_bcast(pub_root_proc_id, pub_write_sw_ngwfs)
    call comms_bcast(pub_root_proc_id, pub_write_max_l)
    call comms_bcast(pub_root_proc_id, pub_read_max_l)
    call comms_bcast(pub_root_proc_id, pub_extra_n_sw)
    call comms_bcast(pub_root_proc_id, pub_paw)
    call comms_bcast(pub_root_proc_id, pub_fine_grid_scale)
    call comms_bcast(pub_root_proc_id, pub_dbl_grid_scale)
    call comms_bcast(pub_root_proc_id, pub_aug_funcs_recip)
    call comms_bcast(pub_root_proc_id, pub_pseudo_path)
    call comms_bcast(pub_root_proc_id, pub_dispersion)
    call comms_bcast(pub_root_proc_id, pub_vdw_dcoeff)
    call comms_bcast(pub_root_proc_id, pub_vdw_radial_cutoff)
    call comms_bcast(pub_root_proc_id, pub_hub_max_iter)
    call comms_bcast(pub_root_proc_id, pub_hub_energy_tol)
    call comms_bcast(pub_root_proc_id, pub_hub_conv_win)
    call comms_bcast(pub_root_proc_id, pub_hub_proj_mixing)
    call comms_bcast(pub_root_proc_id, pub_hub_functional)
    call comms_bcast(pub_root_proc_id, pub_hub_tensor_corr)
    call comms_bcast(pub_root_proc_id, pub_hub_ngwf_spin_thr)
    call comms_bcast(pub_root_proc_id, pub_hub_on_the_fly)
    call comms_bcast(pub_root_proc_id, pub_hub_read_projectors)
    call comms_bcast(pub_root_proc_id, pub_hubbard_compute_u_or_j)
    call comms_bcast(pub_root_proc_id, pub_hubbard_j_minority_term)
    call comms_bcast(pub_root_proc_id, pub_hub_tensor_forces)

    call comms_bcast(pub_root_proc_id, pub_fftbox_batch_size)

    ! gibo: cDFT stuff ===== START
    call comms_bcast(pub_root_proc_id, pub_cdft_atom_charge)
    call comms_bcast(pub_root_proc_id, pub_cdft_atom_spin)
    call comms_bcast(pub_root_proc_id, pub_cdft_group_charge_diff)
    call comms_bcast(pub_root_proc_id, pub_cdft_group_spin_diff)
    call comms_bcast(pub_root_proc_id, pub_cdft_group_charge_acceptor)
    call comms_bcast(pub_root_proc_id, pub_cdft_group_charge_donor)
    call comms_bcast(pub_root_proc_id, pub_cdft_group_spin_acceptor)
    call comms_bcast(pub_root_proc_id, pub_cdft_group_spin_donor)
    call comms_bcast(pub_root_proc_id, pub_cdft_hubbard)
    call comms_bcast(pub_root_proc_id, pub_ci_cdft)
    call comms_bcast(pub_root_proc_id, pub_ci_cdft_num_conf)
    call comms_bcast(pub_root_proc_id, pub_cdft_print_all_occ)
    call comms_bcast(pub_root_proc_id, pub_cdft_read_projectors)
    call comms_bcast(pub_root_proc_id, pub_cdft_group_charge_up_only)
    call comms_bcast(pub_root_proc_id, pub_cdft_group_charge_down_only)
    call comms_bcast(pub_root_proc_id, pub_cdft_multi_proj)
    call comms_bcast(pub_root_proc_id, pub_cdft_modes)
    call comms_bcast(pub_root_proc_id, pub_cdft_group_charge_diff_target)
    call comms_bcast(pub_root_proc_id, pub_cdft_group_spin_diff_target)
    call comms_bcast(pub_root_proc_id, pub_cdft_group_charge_acceptor_target)
    call comms_bcast(pub_root_proc_id, pub_cdft_group_charge_donor_target)
    call comms_bcast(pub_root_proc_id, pub_cdft_group_spin_acceptor_target)
    call comms_bcast(pub_root_proc_id, pub_cdft_group_spin_donor_target)
    call comms_bcast(pub_root_proc_id, pub_maxit_cdft_u_cg)
    call comms_bcast(pub_root_proc_id, pub_cdft_cg_type)
    call comms_bcast(pub_root_proc_id, pub_cdft_cg_threshold)
    call comms_bcast(pub_root_proc_id, pub_cdft_trial_length)
    call comms_bcast(pub_root_proc_id, pub_cdft_cg_max)
    call comms_bcast(pub_root_proc_id, pub_cdft_max_grad)
    call comms_bcast(pub_root_proc_id, pub_cdft_elec_energy_tol)
    call comms_bcast(pub_root_proc_id, pub_cdft_cg_max_step)
    call comms_bcast(pub_root_proc_id, pub_cdft_guru)
    call comms_bcast(pub_root_proc_id, pub_cdft_continuation)
    call comms_bcast(pub_root_proc_id, pub_cdft_tight)
    call comms_bcast(pub_root_proc_id, pub_cdft_write_potentials)
    ! gibo: cDFT stuff ===== END

    ! ep: mermin section
    call comms_bcast(pub_root_proc_id, pub_mermin)
    call comms_bcast(pub_root_proc_id, pub_check_mermin)
    call comms_bcast(pub_root_proc_id, pub_mermin_temp)
    call comms_bcast(pub_root_proc_id, pub_mermin_smearing_width)
    call comms_bcast(pub_root_proc_id, pub_mermin_cheb)
    call comms_bcast(pub_root_proc_id, pub_mermin_mu_sq)
    call comms_bcast(pub_root_proc_id, maxit_mermin)
    call comms_bcast(pub_root_proc_id, pub_firstit_mermin)
    call comms_bcast(pub_root_proc_id, mermin_threshold_orig)
    call comms_bcast(pub_root_proc_id, pub_mermin_cg_type)
    call comms_bcast(pub_root_proc_id, pub_mermin_cg_max_step)
    call comms_bcast(pub_root_proc_id, pub_mermin_check_trial_steps)
    call comms_bcast(pub_root_proc_id, pub_mermin_round_evals)
    call comms_bcast(pub_root_proc_id, pub_mermin_free_energy_thres)
    ! ep: mermin section

    ! ebl DMFT start
    call comms_bcast(pub_root_proc_id, pub_hub_proj_read_only)
    call comms_bcast(pub_root_proc_id, pub_dmft_fully_sc)
    call comms_bcast(pub_root_proc_id, pub_dmft_fully_sc_h)
    call comms_bcast(pub_root_proc_id, pub_dmft_kernel)
    call comms_bcast(pub_root_proc_id, pub_dmft_nbo)
    call comms_bcast(pub_root_proc_id, pub_dmft_nkpoints)
    call comms_bcast(pub_root_proc_id, pub_dmft_optics)
    call comms_bcast(pub_root_proc_id, pub_dmft_order_proj)
    call comms_bcast(pub_root_proc_id, pub_dmft_plot_real_space)
    call comms_bcast(pub_root_proc_id, pub_dmft_points)
    call comms_bcast(pub_root_proc_id, pub_dmft_purify_sc)
    call comms_bcast(pub_root_proc_id, pub_dmft_sc)
    call comms_bcast(pub_root_proc_id, pub_dmft_spoil_kernel)
    call comms_bcast(pub_root_proc_id, pub_dmft_switch_off_proj_order)

    call comms_bcast(pub_root_proc_id, pub_dmft_complex_freq)
    call comms_bcast(pub_root_proc_id, pub_dmft_cutoff_small)
    call comms_bcast(pub_root_proc_id, pub_dmft_dos_max)
    call comms_bcast(pub_root_proc_id, pub_dmft_dos_min)
    call comms_bcast(pub_root_proc_id, pub_dmft_emax)
    call comms_bcast(pub_root_proc_id, pub_dmft_emin)
    call comms_bcast(pub_root_proc_id, pub_dmft_kernel_mix)
    call comms_bcast(pub_root_proc_id, pub_dmft_kpoints_sym)
    call comms_bcast(pub_root_proc_id, pub_dmft_ks_shift)
    call comms_bcast(pub_root_proc_id, pub_dmft_mu_diff_max)
    call comms_bcast(pub_root_proc_id, pub_dmft_mu_order)
    call comms_bcast(pub_root_proc_id, pub_dmft_nmu_loop)
    call comms_bcast(pub_root_proc_id, pub_dmft_nval)
    call comms_bcast(pub_root_proc_id, pub_dmft_optics_i1)
    call comms_bcast(pub_root_proc_id, pub_dmft_optics_i2)
    call comms_bcast(pub_root_proc_id, pub_dmft_optics_window)
    call comms_bcast(pub_root_proc_id, pub_dmft_paramagnetic)
    call comms_bcast(pub_root_proc_id, pub_dmft_rotate_green)
    call comms_bcast(pub_root_proc_id, pub_dmft_scaling_cutoff)
    call comms_bcast(pub_root_proc_id, pub_dmft_scaling_meth)
    call comms_bcast(pub_root_proc_id, pub_dmft_scaling_nmpi)
    call comms_bcast(pub_root_proc_id, pub_dmft_scaling_tail)
    call comms_bcast(pub_root_proc_id, pub_dmft_skip_energy)
    call comms_bcast(pub_root_proc_id, pub_dmft_smear)
    call comms_bcast(pub_root_proc_id, pub_dmft_smear_eta)
    call comms_bcast(pub_root_proc_id, pub_dmft_smear_shift)
    call comms_bcast(pub_root_proc_id, pub_dmft_smear_T)
    call comms_bcast(pub_root_proc_id, pub_dmft_smear_w)
    call comms_bcast(pub_root_proc_id, pub_dmft_temp)
    call comms_bcast(pub_root_proc_id, pub_dmft_win)
    call comms_bcast(pub_root_proc_id, pub_dmft_write)
    call comms_bcast(pub_root_proc_id, pub_dmft_read)

    ! DMFT keywords presently set as devel_code; will eventually be declared
    ! here
    ! call comms_bcast(pub_root_proc_id, pub_dmft_doping)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_scaling_cutoff_h)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_scaling_tol)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_scaling_maxspace)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_scaling_iter)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_impose_same_coeffs)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_impose_chem_spin)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_split)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_splitk)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_kpoints_kernel_gamma)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_lin_scaling)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_invert_overlap)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_local_scratch)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_optics_x1)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_optics_y1)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_optics_z1)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_embed_iter)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_embed_mix)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_in_bohr)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_plot_all_proj)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_integrate_green)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_plot_real_space_sigma)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_norm_proj)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_chem_shift)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_cutoff_tail)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_free_green_frequ)
    ! call comms_bcast(pub_root_proc_id, pub_dmft_gpu_num)
    ! ebl: DMFT end

    ! gcc32: LR_PHONONS start
    call comms_bcast(pub_root_proc_id, pub_lr_phonons_calculate)
    call comms_bcast(pub_root_proc_id, pub_lr_phonons_zero_dim)
    call comms_bcast(pub_root_proc_id, pub_lr_phonons_kernel_cutoff)
    call comms_bcast(pub_root_proc_id, pub_lr_phonons_restart)
    ! gcc32 LR_PHONONS end

    ! tjz07: LR_TDDFT start
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_calculate)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_num_states)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_cg_threshold)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_maxit_cg)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_maxit_pen)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_reset_cg)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_penalty_tol)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_kernel_cutoff)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_write_densities)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_write_kernels)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_restart)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_triplet)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_projector)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_joint_set)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_preopt)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_preopt_iter)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_analysis)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_num_conv_states)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_check_conv_iter)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_precond)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_precond_tol)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_precond_iter)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_HOMO_num)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_LUMO_num)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_xc_finite_diff)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_subsystem_coupling)
    call comms_bcast(pub_root_proc_id, pub_lr_optical_permittivity)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_RPA)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_init_random)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_init_max_overlap)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_spectrum_smear)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_restart_from_TDA)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_ct_length)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_mlwf_analysis)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_mom_mat_els)
    call comms_bcast(pub_root_proc_id, pub_lr_tddft_penalty_func)
    ! tjz07: LR_TDDFT end

    ! jhl52: QNTO start
    call comms_bcast(pub_root_proc_id, pub_qnto_analysis)
    call comms_bcast(pub_root_proc_id, pub_qnto_num_transition)
    call comms_bcast(pub_root_proc_id, pub_qnto_write_orbitals)
    call comms_bcast(pub_root_proc_id, pub_qnto_svd_method)
    call comms_bcast(pub_root_proc_id, pub_qnto_nbo_proj)
    call comms_bcast(pub_root_proc_id, pub_qnto_ref_dir)
    call comms_bcast(pub_root_proc_id, pub_qnto_num_ref_states)
    call comms_bcast(pub_root_proc_id, pub_qnto_num_core_atoms)
    ! jhl52: QNTO end

    call comms_bcast(pub_root_proc_id, pub_permit_unusual_ngwf_count)

    call comms_bcast(pub_root_proc_id, pub_coulomb_radius)
    call comms_bcast(pub_root_proc_id, pub_coulomb_length)
    call comms_bcast(pub_root_proc_id, pub_coulomb_cutoff_type)
    call comms_bcast(pub_root_proc_id, pub_coulomb_cutoff_write_int)
    call comms_bcast(pub_root_proc_id, pub_mt_cutoff)
    call comms_bcast(pub_root_proc_id, pub_do_properties)
    call comms_bcast(pub_root_proc_id, pub_num_eigenvalues)
    call comms_bcast(pub_root_proc_id, pub_homo_dens_plot)
    call comms_bcast(pub_root_proc_id, pub_lumo_dens_plot)
    call comms_bcast(pub_root_proc_id, pub_homo_plot)
    call comms_bcast(pub_root_proc_id, pub_lumo_plot)
    call comms_bcast(pub_root_proc_id, pub_dos_smear)
    call comms_bcast(pub_root_proc_id, pub_popn_calculate)
    call comms_bcast(pub_root_proc_id, pub_lowdin_popn_calculate)
    call comms_bcast(pub_root_proc_id, pub_popn_bond_cutoff)
    call comms_bcast(pub_root_proc_id, pub_popn_mulliken_partial)
    call comms_bcast(pub_root_proc_id, pub_ngwf_analysis)
    call comms_bcast(pub_root_proc_id, pub_polarisation_calculate)
    call comms_bcast(pub_root_proc_id, pub_polarisation_berry)
    call comms_bcast(pub_root_proc_id, pub_polarisation_local)
    call comms_bcast(pub_root_proc_id, pub_polarisation_simcell_calculate)
    call comms_bcast(pub_root_proc_id, pub_polarisation_simcell_refpt)
    call comms_bcast(pub_root_proc_id, pub_efield_calculate)
    call comms_bcast(pub_root_proc_id, pub_spread_calculate)
    call comms_bcast(pub_root_proc_id, pub_anharmonic_calculate)
    call comms_bcast(pub_root_proc_id, pub_print_potential_noxc)
    call comms_bcast(pub_root_proc_id, pub_pdos_optados_output)
    call comms_bcast(pub_root_proc_id, pub_pdos_max_l)
    call comms_bcast(pub_root_proc_id, pub_pdos_reduce_sws)
    call comms_bcast(pub_root_proc_id, pub_pdos_max_n)
    call comms_bcast(pub_root_proc_id, pub_pdos_stride_n)
    call comms_bcast(pub_root_proc_id, pub_pdos_construct_basis)
    call comms_bcast(pub_root_proc_id, pub_pdos_lowdin)
    call comms_bcast(pub_root_proc_id, pub_pdos_lcao_optimize)
    call comms_bcast(pub_root_proc_id, pub_pdos_orth_atom_blocks)
    call comms_bcast(pub_root_proc_id, pub_pdos_output_basis)
    call comms_bcast(pub_root_proc_id, pub_pdos_output_swopt_kernham)
    call comms_bcast(pub_root_proc_id, pub_pdos_pseudoatomic)
    call comms_bcast(pub_root_proc_id, pub_pdos_sum_mag)

    call comms_bcast(pub_root_proc_id, pub_geom_modulus_est)
    call comms_bcast(pub_root_proc_id, pub_geom_frequency_est)
    call comms_bcast(pub_root_proc_id, pub_geom_energy_tol)
    call comms_bcast(pub_root_proc_id, pub_geom_force_tol)
    call comms_bcast(pub_root_proc_id, pub_geom_disp_tol)
    call comms_bcast(pub_root_proc_id, pub_geom_max_iter)
    call comms_bcast(pub_root_proc_id, pub_geom_convergence_win)
    call comms_bcast(pub_root_proc_id, pub_geom_continuation)
    call comms_bcast(pub_root_proc_id, pub_geom_backup_iter)
    call comms_bcast(pub_root_proc_id, pub_geom_reset_dk_ngwfs_iter)
    call comms_bcast(pub_root_proc_id, pub_geom_lbfgs_max_updates)
    call comms_bcast(pub_root_proc_id, pub_geom_lbfgs_block_length)
    call comms_bcast(pub_root_proc_id, pub_geom_reuse_dk_ngwfs)
    call comms_bcast(pub_root_proc_id, pub_geom_print_inv_hessian)
    call comms_bcast(pub_root_proc_id, pub_geom_method)
    call comms_bcast(pub_root_proc_id, pub_geom_output_detail)
    call comms_bcast(pub_root_proc_id, pub_geom_lbfgs)
    ! lam81
    call comms_bcast(pub_root_proc_id, pub_geom_precond_type)
    call comms_bcast(pub_root_proc_id, pub_geom_precond_exp_c_stab)
    call comms_bcast(pub_root_proc_id, pub_geom_precond_scale_cell)
    call comms_bcast(pub_root_proc_id, pub_geom_precond_exp_A)
    call comms_bcast(pub_root_proc_id, pub_geom_precond_exp_r_NN)
    call comms_bcast(pub_root_proc_id, pub_geom_precond_exp_r_cut)
    call comms_bcast(pub_root_proc_id, pub_geom_precond_exp_mu)
    call comms_bcast(pub_root_proc_id, pub_geom_precond_ff_c_stab)
    call comms_bcast(pub_root_proc_id, pub_geom_precond_ff_r_cut)

    call comms_bcast(pub_root_proc_id, tssearch_method)
    call comms_bcast(pub_root_proc_id, tssearch_lstqst_protocol)
    call comms_bcast(pub_root_proc_id, tssearch_qst_max_iter)
    call comms_bcast(pub_root_proc_id, tssearch_cg_max_iter)
    call comms_bcast(pub_root_proc_id, tssearch_force_tol)
    call comms_bcast(pub_root_proc_id, tssearch_energy_tol)
    call comms_bcast(pub_root_proc_id, tssearch_disp_tol)
    call comms_bcast(pub_root_proc_id, pub_neb_ci_delay)
    call comms_bcast(pub_root_proc_id, pub_neb_print_summary)
    call comms_bcast(pub_root_proc_id, pub_neb_spring_constant)
    call comms_bcast(pub_root_proc_id, pub_neb_read_xyz)
    call comms_bcast(pub_root_proc_id, pub_neb_converge_all)
    call comms_bcast(pub_root_proc_id, pub_neb_update_method)
    call comms_bcast(pub_root_proc_id, pub_neb_max_iter)
    call comms_bcast(pub_root_proc_id, pub_neb_glbfgs_history_size)
    call comms_bcast(pub_root_proc_id, pub_neb_continuation)
    call comms_bcast(pub_root_proc_id, pub_reactant_rootname)
    call comms_bcast(pub_root_proc_id, pub_product_rootname)
    call comms_bcast(pub_root_proc_id, pub_reactant_energy)
    call comms_bcast(pub_root_proc_id, pub_product_energy)
    call comms_bcast(pub_root_proc_id, md_restart)
    call comms_bcast(pub_root_proc_id, md_restart_thermo)
    call comms_bcast(pub_root_proc_id, md_global_restart)
    call comms_bcast(pub_root_proc_id, md_write_out)
    call comms_bcast(pub_root_proc_id, pub_md_properties)
    call comms_bcast(pub_root_proc_id, pub_isthermo)
    call comms_bcast(pub_root_proc_id, md_delta_t)
    call comms_bcast(pub_root_proc_id, md_num_iter)
    call comms_bcast(pub_root_proc_id, md_iter_global)
    call comms_bcast(pub_root_proc_id, md_reset_history)
    call comms_bcast(pub_root_proc_id, md_write_history)
    call comms_bcast(pub_root_proc_id, md_aux_rep)
    call comms_bcast(pub_root_proc_id, md_aux_dkn_t)
    call comms_bcast(pub_root_proc_id, md_aux_beren_tc)
    call comms_bcast(pub_root_proc_id, md_lnv_threshold)
    call comms_bcast(pub_root_proc_id, md_ngwf_threshold)
    call comms_bcast(pub_root_proc_id, md_autocorr)
    call comms_bcast(pub_root_proc_id, mts_xi)
    call comms_bcast(pub_root_proc_id, mts_nstep)
    call comms_bcast(pub_root_proc_id, mts_ngwf_threshold)
    call comms_bcast(pub_root_proc_id, mts_lnv_threshold)
    call comms_bcast(pub_root_proc_id, mts_elec_energy_tol)
    call comms_bcast(pub_root_proc_id, mts_elec_force_tol)
    call comms_bcast(pub_root_proc_id, mts_ngwf_max_grad)
    call comms_bcast(pub_root_proc_id, mts_minit_lnv)
    call comms_bcast(pub_root_proc_id, mts_maxit_lnv)
    call comms_bcast(pub_root_proc_id, mts_maxit_pen)
    call comms_bcast(pub_root_proc_id, mts_maxit_ngwf_cg)
    call comms_bcast(pub_root_proc_id, mts_mix_inc)
    call comms_bcast(pub_root_proc_id, mix_dkn_type)
    call comms_bcast(pub_root_proc_id, mix_dkn_num)
    call comms_bcast(pub_root_proc_id, mix_dkn_init_type)
    call comms_bcast(pub_root_proc_id, mix_dkn_init_num)
    call comms_bcast(pub_root_proc_id, mix_dkn_reset)
    call comms_bcast(pub_root_proc_id, mix_ngwfs_type)
    call comms_bcast(pub_root_proc_id, mix_ngwfs_num)
    call comms_bcast(pub_root_proc_id, mix_ngwfs_init_type)
    call comms_bcast(pub_root_proc_id, mix_ngwfs_init_num)
    call comms_bcast(pub_root_proc_id, mix_ngwfs_reset)
    call comms_bcast(pub_root_proc_id, mix_ngwfs_coeff)
    call comms_bcast(pub_root_proc_id, mix_local_length)
    call comms_bcast(pub_root_proc_id, mix_local_smear)
    call comms_bcast(pub_root_proc_id, cond_energy_range)
    call comms_bcast(pub_root_proc_id, cond_energy_gap)
    call comms_bcast(pub_root_proc_id, cond_read_denskern)
    call comms_bcast(pub_root_proc_id, cond_read_tightbox_ngwfs)
    call comms_bcast(pub_root_proc_id, cond_fixed_shift)
    call comms_bcast(pub_root_proc_id, cond_calc_max_eigen)
    call comms_bcast(pub_root_proc_id, cond_num_states)
    call comms_bcast(pub_root_proc_id, cond_num_extra_states)
    call comms_bcast(pub_root_proc_id, cond_num_extra_its)
    call comms_bcast(pub_root_proc_id, cond_maxit_lnv)
    call comms_bcast(pub_root_proc_id, cond_minit_lnv)
    call comms_bcast(pub_root_proc_id, cond_kernel_cutoff)
    call comms_bcast(pub_root_proc_id, cond_init_shift)
    call comms_bcast(pub_root_proc_id, cond_shift_buffer)
    call comms_bcast(pub_root_proc_id, cond_plot_joint_orbitals)
    call comms_bcast(pub_root_proc_id, cond_plot_vc_orbitals)
    call comms_bcast(pub_root_proc_id, pub_spectra_calculate)
    call comms_bcast(pub_root_proc_id, pub_calc_mom_mat_els)
    call comms_bcast(pub_root_proc_id, pub_calc_nonloc_comm)
    call comms_bcast(pub_root_proc_id, pub_spec_cont_deriv)
    call comms_bcast(pub_root_proc_id, pub_spec_nonloc_fin_diff)
    call comms_bcast(pub_root_proc_id, pub_spectra_print_mat_els)
    call comms_bcast(pub_root_proc_id, pub_spec_scissor_op)
    call comms_bcast(pub_root_proc_id, pub_eels_calculate)
    call comms_bcast(pub_root_proc_id, pub_eels_fine_projectors)
    call comms_bcast(pub_root_proc_id, pub_eels_realspace)
    call comms_bcast(pub_root_proc_id, pub_opt_smear)
    call comms_bcast(pub_root_proc_id, pub_bs_kpoint_path_spacing)
    call comms_bcast(pub_root_proc_id, pub_bs_num_eigenvalues)
    call comms_bcast(pub_root_proc_id, pub_bsunfld_num_eigenvalues)
    call comms_bcast(pub_root_proc_id, pub_bs_unfold)
    call comms_bcast(pub_root_proc_id, pub_bsunfld_transfo)
    call comms_bcast(pub_root_proc_id, pub_bsunfld_num_kpts_path)
    call comms_bcast(pub_root_proc_id, pub_bsunfld_num_atoms_prim)
    call comms_bcast(pub_root_proc_id, pub_bsunfld_restart)
    call comms_bcast(pub_root_proc_id, pub_bs_method)
    call comms_bcast(pub_root_proc_id, pub_perturbative_soc)
    call comms_bcast(pub_root_proc_id, pub_do_tddft)
    call comms_bcast(pub_root_proc_id, pub_tddft_maximum_energy)
    call comms_bcast(pub_root_proc_id, pub_tddft_resolution)
    call comms_bcast(pub_root_proc_id, pub_tddft_propagation_method)
    call comms_bcast(pub_root_proc_id, pub_tddft_sparsity_level)
    call comms_bcast(pub_root_proc_id, pub_tddft_tammdancoff)
    call comms_bcast(pub_root_proc_id, pub_tddft_dipole_kick_strength)
    call comms_bcast(pub_root_proc_id, pub_tddft_xc_functional)
    call comms_bcast(pub_root_proc_id, pub_tddft_hamiltonian_mixing)
    call comms_bcast(pub_root_proc_id, pub_tddft_damping)
    call comms_bcast(pub_root_proc_id, pub_tddft_enforced_idempotency)
    call comms_bcast(pub_root_proc_id, pub_tddft_maxit_hotelling)
    call comms_bcast(pub_root_proc_id, pub_tddft_max_resid_hotelling)
    call comms_bcast(pub_root_proc_id, pub_tddft_inv_overlap_exact)
    call comms_bcast(pub_root_proc_id, pub_swx_dbl_grid)
    call comms_bcast(pub_root_proc_id, pub_swx_c_threshold)
    call comms_bcast(pub_root_proc_id, pub_swx_output_detail)
    call comms_bcast(pub_root_proc_id, pub_use_sph_harm_rot)
    call comms_bcast(pub_root_proc_id, pub_swri_cheb_batchsize)
    call comms_bcast(pub_root_proc_id, pub_swri_swop_smoothing)
    call comms_bcast(pub_root_proc_id, pub_swri_overlap_indirect)
    call comms_bcast(pub_root_proc_id, pub_swri_verbose)
    call comms_bcast(pub_root_proc_id, pub_swri_proximity_sort_point)
    call comms_bcast(pub_root_proc_id, pub_swri_improve_inverse)
    call comms_bcast(pub_root_proc_id, pub_swri_print_eigenvalues)
    call comms_bcast(pub_root_proc_id, pub_swri_assembly_prefix)
    call comms_bcast(pub_root_proc_id, pub_hfx_cutoff)
    call comms_bcast(pub_root_proc_id, pub_hfx_nlpp_for_exchange)
    call comms_bcast(pub_root_proc_id, pub_hfx_debug)
    call comms_bcast(pub_root_proc_id, pub_hfx_read_xmatrix)
    call comms_bcast(pub_root_proc_id, pub_hfx_write_xmatrix)
    call comms_bcast(pub_root_proc_id, pub_hfx_metric)
    call comms_bcast(pub_root_proc_id, pub_hfx_max_l)
    call comms_bcast(pub_root_proc_id, pub_hfx_max_q)
    call comms_bcast(pub_root_proc_id, pub_hfx_bessel_rad_nptsx)
    call comms_bcast(pub_root_proc_id, pub_hfx_memory_limit)
    call comms_bcast(pub_root_proc_id, pub_hfx_memory_weights)
    call comms_bcast(pub_root_proc_id, pub_hfx_use_ri)
    call comms_bcast(pub_root_proc_id, pub_hfx_output_detail)
    call comms_bcast(pub_root_proc_id, pub_cache_limit_for_swops)
    call comms_bcast(pub_root_proc_id, pub_cache_limit_for_swops2)
    call comms_bcast(pub_root_proc_id, pub_cache_limit_for_expansions)
    call comms_bcast(pub_root_proc_id, pub_cache_limit_for_prods)
    call comms_bcast(pub_root_proc_id, pub_cache_limit_for_dknblks)
    call comms_bcast(pub_root_proc_id, pub_ht_stash_size)
    call comms_bcast(pub_root_proc_id, pub_zero_total_force)
    call comms_bcast(pub_root_proc_id, pub_mw_total_force)
    call comms_bcast(pub_root_proc_id, pub_ions_will_move)
    call comms_bcast(pub_root_proc_id, pub_forces_needed)
    call comms_bcast(pub_root_proc_id, pub_cond_calculate_any_task)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_pot)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_qmstar)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_write_vacuum_restart)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_vacuum_qmstar)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_pot_filename)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_dma_min_l)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_dma_max_l)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_vacuum_dma_min_l)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_vacuum_dma_max_l)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_mpole_exclusion_radius)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_pairwise_polarisability)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_thole_a)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_smearing_a)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_perm_scaling)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_fixed_charge)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_dbl_grid)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_polscal)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_repulsive_mm_pot_write)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_repulsive_mm_pot_verbose)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_repulsive_mm_pot_cutoff)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_repulsive_mm_pot_a)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_repulsive_mm_pot_b)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_repulsive_mm_pot_c)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_repulsive_mm_pot_alpha)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_repulsive_mm_pot_beta)
    call comms_bcast(pub_root_proc_id, pub_pol_emb_repulsive_mm_pot_r0)
    call comms_bcast(pub_root_proc_id, pub_is_implicit_solvent)
    call comms_bcast(pub_root_proc_id, pub_is_smeared_ion_rep)
    call comms_bcast(pub_root_proc_id, pub_is_separate_restart_files)
    call comms_bcast(pub_root_proc_id, pub_is_auto_solvation)
    call comms_bcast(pub_root_proc_id, pub_is_include_apolar)
    call comms_bcast(pub_root_proc_id, pub_is_apolar_sasa_definition)
    call comms_bcast(pub_root_proc_id, pub_is_apolar_method)
    call comms_bcast(pub_root_proc_id, pub_is_density_threshold)
    call comms_bcast(pub_root_proc_id, pub_is_solvation_beta)
    call comms_bcast(pub_root_proc_id, pub_is_density_min_threshold)
    call comms_bcast(pub_root_proc_id, pub_is_density_max_threshold)
    call comms_bcast(pub_root_proc_id, pub_is_soft_sphere_delta)
    call comms_bcast(pub_root_proc_id, pub_is_soft_sphere_scale)
    call comms_bcast(pub_root_proc_id, pub_is_bulk_permittivity)
    call comms_bcast(pub_root_proc_id, pub_is_smeared_ion_width)
    call comms_bcast(pub_root_proc_id, pub_is_core_width)
    call comms_bcast(pub_root_proc_id, pub_is_surface_thickness)
    call comms_bcast(pub_root_proc_id, pub_is_solvent_surf_tension)
    call comms_bcast(pub_root_proc_id, pub_is_solvent_pressure)
    call comms_bcast(pub_root_proc_id, pub_is_apolar_scaling_factor)
    call comms_bcast(pub_root_proc_id, pub_is_multigrid_verbose)
    call comms_bcast(pub_root_proc_id, pub_is_multigrid_verbose_y)
    call comms_bcast(pub_root_proc_id, pub_is_multigrid_verbose_z)
    call comms_bcast(pub_root_proc_id, pub_is_bc_coarseness)
    call comms_bcast(pub_root_proc_id, pub_is_bc_surface_coarseness)
    call comms_bcast(pub_root_proc_id, pub_is_bc_threshold)
    call comms_bcast(pub_root_proc_id, pub_is_bc_allow_frac_charge)
    call comms_bcast(pub_root_proc_id, pub_is_check_solv_energy_grad)
    call comms_bcast(pub_root_proc_id, pub_is_solvation_method)
    call comms_bcast(pub_root_proc_id, pub_is_dielectric_model)
    call comms_bcast(pub_root_proc_id, pub_is_dielectric_function)
    call comms_bcast(pub_root_proc_id, pub_is_solvation_output_detail)
    call comms_bcast(pub_root_proc_id, pub_is_pbe)
    call comms_bcast(pub_root_proc_id, pub_is_pbe_temperature)
    call comms_bcast(pub_root_proc_id, pub_is_sc_steric_magnitude)
    call comms_bcast(pub_root_proc_id, pub_is_sc_steric_cutoff)
    call comms_bcast(pub_root_proc_id, pub_is_sc_steric_smoothing_alpha)
    call comms_bcast(pub_root_proc_id, pub_is_steric_pot_type)
    call comms_bcast(pub_root_proc_id, pub_is_hc_steric_smearing)
    call comms_bcast(pub_root_proc_id, pub_is_hc_steric_dens_isovalue)
    call comms_bcast(pub_root_proc_id, pub_is_steric_write)
    call comms_bcast(pub_root_proc_id, pub_is_dielectric_exclusions_smearing)
    call comms_bcast(pub_root_proc_id, pub_is_solvation_properties)
    call comms_bcast(pub_root_proc_id, pub_is_pbe_bc_debye_screening)
    call comms_bcast(pub_root_proc_id, pub_is_pbe_exp_cap)
    call comms_bcast(pub_root_proc_id, pub_is_pbe_energy_tolerance)
    call comms_bcast(pub_root_proc_id, pub_is_pbe_neutralisation_scheme)
    call comms_bcast(pub_root_proc_id, pub_is_restart_vac_from_vac)
    call comms_bcast(pub_root_proc_id, pub_is_emft_cavity)
    call comms_bcast(pub_root_proc_id, pub_openbc_ion_ion)
    call comms_bcast(pub_root_proc_id, pub_openbc_hartree)
    call comms_bcast(pub_root_proc_id, pub_openbc_pspot)
    call comms_bcast(pub_root_proc_id, pub_openbc_pspot_finetune_nptsx)
    call comms_bcast(pub_root_proc_id, pub_openbc_pspot_finetune_alpha)
    call comms_bcast(pub_root_proc_id, pub_openbc_pspot_finetune_f)

    ! JCW: Finite differences (ONETEP finite_difference module)
    call comms_bcast(pub_root_proc_id, pub_finite_difference_order)

    ! JCW: Multigrid solver (updated for DL_MG v2.0)
    call comms_bcast(pub_root_proc_id, pub_mg_granularity_power)
    call comms_bcast(pub_root_proc_id, pub_mg_defco_fd_order)
    call comms_bcast(pub_root_proc_id, pub_mg_tol_res_rel)
    call comms_bcast(pub_root_proc_id, pub_mg_tol_res_abs)
    call comms_bcast(pub_root_proc_id, pub_mg_tol_pot_rel)
    call comms_bcast(pub_root_proc_id, pub_mg_tol_pot_abs)
    call comms_bcast(pub_root_proc_id, pub_mg_tol_vcyc_rel)
    call comms_bcast(pub_root_proc_id, pub_mg_tol_vcyc_abs)
    call comms_bcast(pub_root_proc_id, pub_mg_tol_newton_rel)
    call comms_bcast(pub_root_proc_id, pub_mg_tol_newton_abs)
    call comms_bcast(pub_root_proc_id, pub_mg_max_res_ratio)
    call comms_bcast(pub_root_proc_id, pub_mg_vcyc_smoother_iter_pre)
    call comms_bcast(pub_root_proc_id, pub_mg_vcyc_smoother_iter_post)
    call comms_bcast(pub_root_proc_id, pub_mg_max_iters_defco)
    call comms_bcast(pub_root_proc_id, pub_mg_max_iters_vcycle)
    call comms_bcast(pub_root_proc_id, pub_mg_max_iters_newton)
    call comms_bcast(pub_root_proc_id, pub_mg_use_error_damping)
    call comms_bcast(pub_root_proc_id, pub_mg_pbe_use_fas)
    call comms_bcast(pub_root_proc_id, pub_mg_tol_mu_rel)
    call comms_bcast(pub_root_proc_id, pub_mg_tol_mu_abs)
    call comms_bcast(pub_root_proc_id, pub_mg_continue_on_error)
    call comms_bcast(pub_root_proc_id, pub_mg_tol_cg_res_rel)
    call comms_bcast(pub_root_proc_id, pub_mg_tol_cg_res_abs)
    call comms_bcast(pub_root_proc_id, pub_mg_max_iters_cg)
    call comms_bcast(pub_root_proc_id, pub_mg_use_cg)

    ! JCW: Periodic/mixed BC
    call comms_bcast(pub_root_proc_id, pub_multigrid_bc)
    call comms_bcast(pub_root_proc_id, pub_ion_ion_bc)
    call comms_bcast(pub_root_proc_id, pub_pspot_bc)
    call comms_bcast(pub_root_proc_id, pub_smeared_ion_bc)
    call comms_bcast(pub_root_proc_id, pub_external_bc_from_cube)
    call comms_bcast(pub_root_proc_id, pub_vdw_bc)
    call comms_bcast(pub_root_proc_id, pub_dftb_bc)
    call comms_bcast(pub_root_proc_id, pub_hfx_bc)

    ! lpl: DDEC parameters
    call comms_bcast(pub_root_proc_id, pub_ddec)
    call comms_bcast(pub_root_proc_id, pub_ddec_write)
    call comms_bcast(pub_root_proc_id, pub_ddec_hirshfeld)

    call comms_bcast(pub_root_proc_id, pub_ddec_rad_npts)
    call comms_bcast(pub_root_proc_id, pub_ddec_rad_rcut)

    call comms_bcast(pub_root_proc_id, pub_ddec_maxit)
    call comms_bcast(pub_root_proc_id, pub_ddec_core_maxit)
    call comms_bcast(pub_root_proc_id, pub_ddec_conv_thresh)
    call comms_bcast(pub_root_proc_id, pub_ddec_IH_frac)
    call comms_bcast(pub_root_proc_id, pub_ddec_zero_thresh)
    call comms_bcast(pub_root_proc_id, pub_ddec_use_coredens)

    call comms_bcast(pub_root_proc_id, pub_ddec_core_correction)
    call comms_bcast(pub_root_proc_id, pub_ddec_reshape_dens)
    call comms_bcast(pub_root_proc_id, pub_ddec_core_corr_maxit)

    call comms_bcast(pub_root_proc_id, pub_ddec_ionic_range)
! lpl: Refer to COMMENT#01 in ddec_mod
!    call comms_bcast(pub_root_proc_id, pub_ddec_atomsolve_maxit)
!    call comms_bcast(pub_root_proc_id, pub_ddec_target_radius)
!    call comms_bcast(pub_root_proc_id, pub_ddec_rcomp_econv)
!    call comms_bcast(pub_root_proc_id, pub_ddec_rcomp_maxit)
!    call comms_bcast(pub_root_proc_id, pub_ddec_write_rcomp)

    call comms_bcast(pub_root_proc_id, pub_ddec_refdens_path)
    call comms_bcast(pub_root_proc_id, pub_ddec_refdens_init)
    call comms_bcast(pub_root_proc_id, pub_ddec_renormalize_refdens)
    call comms_bcast(pub_root_proc_id, pub_ddec_c3_refdens)
    call comms_bcast(pub_root_proc_id, pub_ddec_rad_shell_mode)

    call comms_bcast(pub_root_proc_id, pub_ddec_multipole)
    call comms_bcast(pub_root_proc_id, pub_ddec_moment)
! lpl: Refer to COMMENT#07 in ddec_mod
!    call comms_bcast(pub_root_proc_id, pub_ddec_moment_order)

    call comms_bcast(pub_root_proc_id, pub_ddec_interp_rad_dens)
    call comms_bcast(pub_root_proc_id, pub_ddec_avg_rad)
    call comms_bcast(pub_root_proc_id, pub_ddec_min_shell_dens)
    call comms_bcast(pub_root_proc_id, pub_ddec_ref_shell_mode)

    call comms_bcast(pub_root_proc_id, pub_ddec_format_dens)

    call comms_bcast(pub_root_proc_id, pub_ddec_eff_decay_exp)
    call comms_bcast(pub_root_proc_id, pub_ddec_eff_decay_rmin)

    ! aeaa: DDEC aniso
    call comms_bcast(pub_root_proc_id, pub_ddec_aniso)
    call comms_bcast(pub_root_proc_id, pub_ddec_aniso_max_dis)
    call comms_bcast(pub_root_proc_id, pub_ddec_aniso_max_dis_halogen)
    call comms_bcast(pub_root_proc_id, pub_ddec_aniso_error_thres)
    call comms_bcast(pub_root_proc_id, pub_ddec_aniso_error_reduce)

! lpl: Refer to COMMENT#06 in ddec_mod
!    ! lpl: This is related ro the still-dodgy subroutine 'ddec_rmse'
!    call comms_bcast(pub_root_proc_id, pub_ddec_rmse)
    ! lpl: END DDEC parameters

    ! lpl: NBO stuff
    call comms_bcast(pub_root_proc_id, pub_write_nbo)
    call comms_bcast(pub_root_proc_id, pub_nbo_init_lclowdin)
    call comms_bcast(pub_root_proc_id, pub_nbo_write_lclowdin)
    call comms_bcast(pub_root_proc_id, pub_nbo_write_npacomp)
    call comms_bcast(pub_root_proc_id, pub_nbo_scale_dm)
    call comms_bcast(pub_root_proc_id, pub_nbo_write_dipole)
    call comms_bcast(pub_root_proc_id, pub_nbo_scale_spin)

    call comms_bcast(pub_root_proc_id, pub_nbo_pnao_analysis)

    call comms_bcast(pub_root_proc_id, pub_nbo_aopnao_scheme)

    call comms_bcast(pub_root_proc_id, pub_plot_nbo)
    call comms_bcast(pub_root_proc_id, pub_nbo_plotorbtype)
    ! lpl: NBO stuff

    ! mjsp: EDA
    call comms_bcast(pub_root_proc_id, pub_eda_continuation)
    call comms_bcast(pub_root_proc_id, pub_eda_continuation_loaded_dkn)
    call comms_bcast(pub_root_proc_id, pub_eda_continuation_loaded_ngwfs)
    call comms_bcast(pub_root_proc_id, pub_eda_write)
    call comms_bcast(pub_root_proc_id, pub_eda_preptool)
    call comms_bcast(pub_root_proc_id, pub_eda_read_frags)
    call comms_bcast(pub_root_proc_id, pub_eda_read_super)
    call comms_bcast(pub_root_proc_id, pub_eda_mode)
    call comms_bcast(pub_root_proc_id, pub_frag_counter)
    call comms_bcast(pub_root_proc_id, pub_frag_counter2)
    call comms_bcast(pub_root_proc_id, pub_eda_have_sp_frags)
    call comms_bcast(pub_root_proc_id, pub_eda_deltadens)
    call comms_bcast(pub_root_proc_id, pub_eda_frag_isol_pol)
    call comms_bcast(pub_root_proc_id, pub_eda_frag_isol_ct)
    call comms_bcast(pub_root_proc_id, pub_eda_reset_ngwfs_pol)
    call comms_bcast(pub_root_proc_id, pub_eda_reset_ngwfs_ct)
    call comms_bcast(pub_root_proc_id, pub_eda_split_atoms)
!    call comms_bcast(pub_root_proc_id, pub_eda_scfmi_restrict_ngwfs)
    call comms_bcast(pub_root_proc_id, pub_eda_nodiag)
    call comms_bcast(pub_root_proc_id, pub_eda_scfmi)
    call comms_bcast(pub_root_proc_id, pub_eda_scfmi_any)
    call comms_bcast(pub_root_proc_id, pub_eda_frz_x)
    call comms_bcast(pub_root_proc_id, pub_eda_frz_c)
    call comms_bcast(pub_root_proc_id, pub_eda_rep_x)
    call comms_bcast(pub_root_proc_id, pub_eda_rep_c)
    call comms_bcast(pub_root_proc_id, pub_eda_intrafrag_x_frz)
    call comms_bcast(pub_root_proc_id, pub_eda_intrafrag_c_frz)
    call comms_bcast(pub_root_proc_id, pub_eda_dE_x_rep)
    call comms_bcast(pub_root_proc_id, pub_eda_dE_c_rep)

    ! smmd: Transmission coefficients
    call comms_bcast(pub_root_proc_id, pub_etrans_lcr)
    call comms_bcast(pub_root_proc_id, pub_etrans_bulk)
    call comms_bcast(pub_root_proc_id, pub_etrans_same_leads)
    call comms_bcast(pub_root_proc_id, pub_etrans_write_xyz)
    call comms_bcast(pub_root_proc_id, pub_etrans_write_hs)
    call comms_bcast(pub_root_proc_id, pub_etrans_ecmplx)
    call comms_bcast(pub_root_proc_id, pub_etrans_emax)
    call comms_bcast(pub_root_proc_id, pub_etrans_emin)
    call comms_bcast(pub_root_proc_id, pub_etrans_enum)
    call comms_bcast(pub_root_proc_id, pub_etrans_eref)
    call comms_bcast(pub_root_proc_id, pub_etrans_eref_method)
    call comms_bcast(pub_root_proc_id, pub_etrans_calc_lead_pot)
    call comms_bcast(pub_root_proc_id, pub_etrans_lead_nkpoints)
    call comms_bcast(pub_root_proc_id, pub_etrans_lead_disp_tol)
    call comms_bcast(pub_root_proc_id, pub_etrans_num_eigchan)
    call comms_bcast(pub_root_proc_id, pub_etrans_plot_eigchan)
    call comms_bcast(pub_root_proc_id, pub_etrans_lead_lcheck)
    call comms_bcast(pub_root_proc_id, pub_etrans_seed_lead)

    ! vv: anharmonic parameters
    call comms_bcast(pub_root_proc_id, pub_anh_qc_factor)
    call comms_bcast(pub_root_proc_id, pub_anh_acf_factor)
    call comms_bcast(pub_root_proc_id, pub_anh_first_iter)
    call comms_bcast(pub_root_proc_id, pub_anh_last_iter)
    call comms_bcast(pub_root_proc_id, pub_anh_plot_firstfreq)
    call comms_bcast(pub_root_proc_id, pub_anh_plot_lastfreq)
    call comms_bcast(pub_root_proc_id, pub_anh_md_temp)
    call comms_bcast(pub_root_proc_id, pub_anh_plot_all)
    call comms_bcast(pub_root_proc_id, pub_anh_apply_filter)
    call comms_bcast(pub_root_proc_id, pub_anh_type)

    ! jd: DFTB
    call comms_bcast(pub_root_proc_id, pub_dftb)
    call comms_bcast(pub_root_proc_id, pub_dftb_method)
    call comms_bcast(pub_root_proc_id, pub_dftb_method_param_file)
    call comms_bcast(pub_root_proc_id, pub_dftb_common_param_file)
    call comms_bcast(pub_root_proc_id, pub_dftb_coord_cutoff)
    call comms_bcast(pub_root_proc_id, pub_dftb_rep_cutoff)
    call comms_bcast(pub_root_proc_id, pub_dftb_srb_cutoff)
    call comms_bcast(pub_root_proc_id, pub_dftb_overlap_cutoff)
    call comms_bcast(pub_root_proc_id, pub_dftb_overlap_analytical)
    call comms_bcast(pub_root_proc_id, pub_dftb_cartesian_ngwfs)
    call comms_bcast(pub_root_proc_id, pub_dftb_ewald_replicate_xtb)
    call comms_bcast(pub_root_proc_id, pub_dftb_ewald_parameter)
    call comms_bcast(pub_root_proc_id, pub_dftb_calculate_force)

    ! fc: deal with phonon parameters
    call comms_bcast(pub_root_proc_id, pub_phonon_farming_task)
    call comms_bcast(pub_root_proc_id, pub_phonon_sampling)
    call comms_bcast(pub_root_proc_id, pub_phonon_vib_free)
    call comms_bcast(pub_root_proc_id, pub_phonon_disp)
    call comms_bcast(pub_root_proc_id, pub_phonon_fmax)
    call comms_bcast(pub_root_proc_id, pub_phonon_energy_check)
    call comms_bcast(pub_root_proc_id, pub_phonon_tmin)
    call comms_bcast(pub_root_proc_id, pub_phonon_tmax)
    call comms_bcast(pub_root_proc_id, pub_phonon_deltat)
    call comms_bcast(pub_root_proc_id, pub_phonon_min_freq)
    call comms_bcast(pub_root_proc_id, pub_phonon_write_eigenvecs)
    call comms_bcast(pub_root_proc_id, pub_phonon_animate_scale)
    call comms_bcast(pub_root_proc_id, pub_phonon_DOS_min)
    call comms_bcast(pub_root_proc_id, pub_phonon_DOS_max)
    call comms_bcast(pub_root_proc_id, pub_phonon_DOS_delta)
    call comms_bcast(pub_root_proc_id, pub_phonon_DOS)
    call comms_bcast(pub_root_proc_id, pub_phonon_SK)
    if (pub_on_root) pub_have_phonon_disp_list = esdf_block('phonon_disp_list',pub_num_disp)
    call comms_bcast(pub_root_proc_id, pub_have_phonon_disp_list)
    if (pub_have_phonon_disp_list) then
       call comms_bcast(pub_root_proc_id, pub_num_disp)
       allocate(pub_phonon_disp_list(1:pub_num_disp),stat=ierr)
       call utils_alloc_check('get_rundat','pub_phonon_disp_list',ierr)
       if (pub_on_root) then
          do iline=1,pub_num_disp
             read(block_data(iline),*) pub_phonon_disp_list(iline)
          end do
       end if
       call comms_bcast(pub_root_proc_id, pub_phonon_disp_list)
    end if
    if (pub_on_root) pub_have_phonon_except_list = esdf_block('phonon_exception_list',pub_num_except)
    call comms_bcast(pub_root_proc_id, pub_have_phonon_except_list)
    if (pub_have_phonon_except_list) then
       call comms_bcast(pub_root_proc_id, pub_num_except)
       allocate(pub_phonon_iexcept_list(1:4,1:pub_num_except),stat=ierr)
       call utils_alloc_check('get_rundat','pub_phonon_iexcept_list',ierr)
       allocate(pub_phonon_dexcept_list(1:pub_num_except),stat=ierr)
       call utils_alloc_check('get_rundat','pub_phonon_dexcept_list',ierr)
       if (pub_on_root) then
          do iline=1,pub_num_except
             read(block_data(iline),*) pub_phonon_iexcept_list(1:4,iline), &
                  pub_phonon_dexcept_list(iline)
          end do
       end if
       call comms_bcast(pub_root_proc_id, pub_phonon_iexcept_list)
       call comms_bcast(pub_root_proc_id, pub_phonon_dexcept_list)
    end if
    if (pub_on_root) pub_have_phonon_animate_list = esdf_block('phonon_animate_list',pub_num_ani)
    call comms_bcast(pub_root_proc_id, pub_have_phonon_animate_list)
    if (pub_have_phonon_animate_list) then
       call comms_bcast(pub_root_proc_id, pub_num_ani)
       allocate(pub_phonon_animate_list(1:pub_num_ani),stat=ierr)
       call utils_alloc_check('get_rundat','pub_phonon_animate_list',ierr)
       if (pub_on_root) then
          do iline=1,pub_num_ani
             read(block_data(iline),*) pub_phonon_animate_list(iline)
          end do
       end if
       call comms_bcast(pub_root_proc_id, pub_phonon_animate_list)
    end if
    if (pub_on_root) pub_have_supercell = esdf_block('supercell',pub_nat_unit)
    call comms_bcast(pub_root_proc_id, pub_have_supercell)
    if (pub_have_supercell) then
       call comms_bcast(pub_root_proc_id, pub_nat_unit)
       allocate(pub_supercell_unit_list(1:pub_nat_unit-1),stat=ierr)
       call utils_alloc_check('get_rundat','pub_supercell_unit_list',ierr)
       if (pub_on_root) then
          do iline=1,pub_nat_unit
             if (iline==1) then
                read(block_data(iline),*) pub_supercell(1:3)
             else
                read(block_data(iline),*) pub_supercell_unit_list(iline-1)
             end if
          end do
       end if
       call comms_bcast(pub_root_proc_id, pub_supercell)
       call comms_bcast(pub_root_proc_id, pub_supercell_unit_list)
       pub_nat_unit=pub_nat_unit-1
    else
       pub_supercell=1
    end if
    if (pub_on_root) pub_have_phonon_grid = esdf_block('phonon_grid',pub_num_grid)
    call comms_bcast(pub_root_proc_id, pub_have_phonon_grid)
    if (pub_have_phonon_grid) then
       if (pub_on_root) then
          read(block_data(1),*) pub_phonon_grid(1:3)
       end if
       call comms_bcast(pub_root_proc_id, pub_phonon_grid)
    else
       pub_phonon_grid=(/1,1,1/)
    end if
    if (pub_on_root) pub_have_phonon_qpoints = esdf_block('phonon_qpoints',pub_num_qpoints)
    call comms_bcast(pub_root_proc_id, pub_have_phonon_qpoints)
    if (pub_have_phonon_qpoints) then
       call comms_bcast(pub_root_proc_id, pub_num_qpoints)
       allocate(pub_phonon_qpoints(1:3,1:pub_num_qpoints),stat=ierr)
       call utils_alloc_check('get_rundat','pub_phonon_qpoints',ierr)
       if (pub_on_root) then
          do iline=1,pub_num_qpoints
             read(block_data(iline),*) pub_phonon_qpoints(1:3,iline)
          end do
       end if
       call comms_bcast(pub_root_proc_id, pub_phonon_qpoints)
    end if

    pub_hubbard = .false. ! ddor: by default we do not carry out DFT+U
    pub_hub_calculating_u = .false. ! ddor: by default we renew the Hamiltonian
    if ( pub_hub_max_iter .gt. 1 ) then
       pub_task = 'HUBBARDSCF'
    elseif (pub_hub_max_iter .eq. -1) then
       pub_hub_calculating_u = .true.
    endif

    pub_cdft = .false.              ! do not carry out cDFT by default

    ! az, jd: dma
    call comms_bcast(pub_root_proc_id,pub_dma_metric)
    call comms_bcast(pub_root_proc_id,pub_dma_calculate)
    call comms_bcast(pub_root_proc_id,pub_dma_output_potential)
    call comms_bcast(pub_root_proc_id,pub_dma_output_potential_reference)
    call comms_bcast(pub_root_proc_id,pub_dma_scale_charge)
    call comms_bcast(pub_root_proc_id,pub_dma_target_num_val_elec)
    call comms_bcast(pub_root_proc_id,pub_dma_bessel_averaging)
    call comms_bcast(pub_root_proc_id,pub_dma_multipole_scaling)
    call comms_bcast(pub_root_proc_id,pub_dma_dipole_scaling)
    call comms_bcast(pub_root_proc_id,pub_dma_quadrupole_scaling)
    call comms_bcast(pub_root_proc_id,pub_dma_precise_gdma_output)
    call comms_bcast(pub_root_proc_id,pub_dma_max_l)
    call comms_bcast(pub_root_proc_id,pub_dma_max_q)
    call comms_bcast(pub_root_proc_id,pub_dma_use_ri)

    call comms_bcast(pub_root_proc_id, pub_turn_off_ewald)
    call comms_bcast(pub_root_proc_id, pub_turn_off_hartree)
    call comms_bcast(pub_root_proc_id, pub_hubbard_unify_sites)
    call comms_bcast(pub_root_proc_id, pub_dft_nu_opt_u1_only)
    call comms_bcast(pub_root_proc_id, pub_dft_nu_opt_u2_only)
    call comms_bcast(pub_root_proc_id, pub_dft_nu_continuation)

    ! JCW: Additional exchange-correlation functional options
    ! JCW: * Minimum threshold value for tau (KE density)
    call comms_bcast(pub_root_proc_id, pub_xc_min_tau)
    ! JCW: * Alternative functional for initial guess (meta-GGAs only)
    call comms_bcast(pub_root_proc_id, pub_xc_initial_functional)

    call comms_bcast(pub_root_proc_id, pub_rand_seed)

    ! cks: Hyper Hartree-Fock parameters
    call comms_bcast(pub_root_proc_id, pub_hhf_nstates)

    ! jcap: Embedding mean field theory parameters
    call comms_bcast(pub_root_proc_id, pub_active_xc_functional)
    call comms_bcast(pub_root_proc_id, pub_emft)
    call comms_bcast(pub_root_proc_id, pub_emft_follow)
    call comms_bcast(pub_root_proc_id, pub_emft_lnv_only)
    call comms_bcast(pub_root_proc_id, pub_active_region)
    call comms_bcast(pub_root_proc_id, pub_block_orthogonalise)
    call comms_bcast(pub_root_proc_id, pub_singlet_triplet_split)
    call comms_bcast(pub_root_proc_id, pub_build_bo)
    call comms_bcast(pub_root_proc_id, pub_emft_lnv_steps)

    !rjc: Electron Localisation Descriptor and Kinetic Energy Density output
    call comms_bcast(pub_root_proc_id, pub_eld_calculate)
    call comms_bcast(pub_root_proc_id, pub_eld_function)
    call comms_bcast(pub_root_proc_id, pub_ke_density_calculate)

  end subroutine get_rundat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rundat_threads_init

    use comms, only: pub_on_root, pub_total_num_procs, pub_root_process_count, &
                   pub_root_proc_id, comms_bcast

    use constants, only: stdout
    use esdf, only: esdf_integer
    implicit none

    integer          :: default_threads
    character(len=3) :: mpi_string
    character(len=:), allocatable :: process_string

!$  integer, external :: omp_get_num_procs
!$  external :: omp_set_dynamic, omp_set_num_threads

    if (pub_on_root) then

      ! set the default number of threads
      default_threads = internal_default_threads()

      ! read input file settings
      pub_threads_max =  esdf_integer('threads_max',default_threads)
      pub_threads_num_fftboxes = esdf_integer('threads_num_fftboxes', &
           default_threads)
      pub_threads_per_fftbox = esdf_integer('threads_per_fftbox',1)
      pub_threads_per_cellfft = esdf_integer('threads_per_cellfft',1)
      pub_threads_max_possible = max(pub_threads_max,pub_threads_per_cellfft, &
           pub_threads_per_fftbox * pub_threads_num_fftboxes)

      ! smmd: OpenMP threads dedicated to internal mkl threading
      pub_threads_num_mkl = esdf_integer('threads_num_mkl',1)

    end if

    ! communicate or compute all threads-related public variables
    call comms_bcast(pub_root_proc_id, pub_threads_max)
    call comms_bcast(pub_root_proc_id, pub_threads_per_fftbox)
    call comms_bcast(pub_root_proc_id, pub_threads_num_fftboxes)
    call comms_bcast(pub_root_proc_id, pub_threads_per_cellfft)
    call comms_bcast(pub_root_proc_id, pub_threads_max_possible)

    call comms_bcast(pub_root_proc_id, pub_threads_num_mkl)

!$  if (pub_threads_per_fftbox.gt.1) then
!$     pub_threads_fftbox = .true.
!$  else
       pub_threads_fftbox = .false.
!$  end if

!$  if (pub_threads_per_cellfft.gt.1) then
!$     pub_threads_cellfft = .true.
!$  else
       pub_threads_cellfft = .false.
!$  end if

    ! print info on MPI processes
    if (pub_on_root) then

      mpi_string = 'MPI'
#ifdef MPI
      process_string = mpi_string // ' process'
      ! Write out how many processes we are using.
      if (pub_total_num_procs == 1) then
        write(stdout,'(3x,a)') 'Running with 1 MPI process.'
      else
        write(stdout,'(3x,a,i0,a)') 'Running with ', pub_total_num_procs, &
                                    ' MPI processes.'
      end if

      write(stdout,'(3x,a,i0,2a)') 'There are ', pub_root_process_count, &
            ' ' // mpi_string , &
            ' processes running on the same node as the root process.'
#else
      process_string = 'process'
      write(stdout,'(3x,3a)') 'This binary does not support ', mpi_string, &
                            ' parallelization.'
#endif

!! jme: Test disabled.  We don't seem to be able to control the scope of
!!      omp_get_num_procs() (total # cores, # cores in NUMA region, or even
!!      the integer part of # cores / # MPI processes).  The warning pops up
!!      so frequently that will not be taken seriously.
!!
!!      ! kaw: Sanity check to make sure that the number of cores requested does not
!!      !      exceed the number available.
!!!$    if (pub_threads_max_possible.gt.omp_get_num_procs()) then
!!!$      write(stdout,'(3x,a)') 'WARNING: The number of threads required for this &
!!!$          &calculation may be larger'
!!!$      write(stdout,'(3x,a,i0,a)') 'than the number of cores available (', &
!!!$          omp_get_num_procs(),')'
!!!$    end if
!!

    end if

    ! kaw: Set up OpenMP threads
!$  call omp_set_dynamic(.false.)
!$  call omp_set_num_threads(pub_threads_max)

#ifndef FFTW3_NO_OMP
    if (pub_threads_fftbox) then
!!$    call omp_set_nested(1)
!!$    call omp_set_max_active_levels(2)
    else
!!$    call omp_set_nested(0)
!!$    call omp_set_max_active_levels(1)
    end if

!$  if (pub_on_root) then
#ifdef MPI
!$    write(stdout,'(3x,a)') 'Each ' // mpi_string // ' process is using: '
#else
!$    write(stdout,'(3x,a)') 'The process is using: '
#endif
!$    write(stdout,'(12x,i4,a)') pub_threads_per_cellfft, &
!$         ' threads for simulation cell FFTs.'
!$    write(stdout,'(12x,i4,a)') pub_threads_per_fftbox, &
!$         ' threads for parallel FFT box operations.'
!$    write(stdout,'(12x,i4,a)') pub_threads_num_fftboxes, &
!$         ' threads for loops over batched FFT box operations.'
!$    write(stdout,'(12x,i4,a)') pub_threads_max, &
!$         ' threads in other parallel regions.'
!$  end if
#else
!!$  call omp_set_nested(0)
!!$  call omp_set_max_active_levels(1)
    if (pub_threads_fftbox) then
      if (pub_on_root) then
        write(stdout,'(3x,a)') 'Cannot perform parallel FFT box operations without'
        write(stdout,'(3x,a)') 'the OpenMP enabled version of the FFTW library.'
        write(stdout,'(3x,a)') '*Defaulting to standard threaded version of ONETEP*'
      end if
      pub_threads_fftbox = .false.
      pub_threads_per_fftbox = 1
    end if
!$  if (pub_on_root) write(stdout,'(3x,2a,i2,a)') 'Each ' // process_string , &
!$            ' is using ', pub_threads_max,' threads.'
#endif

    ! jme: Final warnings displayed from the root proc
    ! # of threads > 1 but no OpenMP is present.
    if (pub_on_root) then
!$    if (.False.) then
        write(stdout,'(3x,a)') 'This binary does not support OpenMP threading.'
        if (pub_threads_max_possible > 1) then
          write(stdout,'(3x,a)') 'WARNING: One or more threading variables have been'
          write(stdout,'(3x,a)') 'set to values > 1, but your ONETEP compilation did'
          write(stdout,'(3x,a)') 'not enable OpenMP threading.  Therefore, these'
          write(stdout,'(3x,a)') 'variables will have no effect, and the code will run'
          write(stdout,'(3x,a)') 'with just one thread per (MPI) process.'
          write(stdout,'(3x,a)') 'ONETEP may be recompiled (with a different or'
          write(stdout,'(3x,a)') 'modified config file) in order to enable threading.'
        end if
!$    endif
!$    if (pub_threads_max_possible.gt.128) then
!$      write(stdout,'(3x,a)') 'WARNING: The number of threads you are using is > 128.'
!$      write(stdout,'(3x,a)') 'This may cause issues with timer_mod.F90 due to the'
!$      write(stdout,'(3x,a)') 'wired value for max_threads.'
!$    end if
    end if

  contains

    ! returns the default number of threads to be used.
    integer function internal_default_threads() result(default_threads)

!$    use constants, only: NORMAL, stdout
!$    use comms, only: pub_total_num_procs
      implicit none

!$    integer, external :: omp_get_max_threads
!$    integer :: env_threads, omp_num_threads_int
!$    character(len=80) :: omp_num_threads
!$    character(len=3) :: mpi_string
!$    integer :: ierr

      ! jme: Determine a sensible default
!$    if (pub_output_detail > NORMAL) then
!$      write(stdout,'(3x,a)') 'Determining the default number of threads...'
!$    end if
#ifdef DEFAULT_THREADS
      default_threads = DEFAULT_THREADS
!$    if (pub_output_detail > NORMAL) then
!$      write(stdout,'(6x,a,i0)') 'Default threads set at compilation time: ', &
!$           default_threads
!$    end if
#else
      default_threads = 1
!$    if (pub_output_detail > NORMAL) then
!$      write(stdout,'(6x,a)') 'No default threads set at compilation time:' // &
!$           ' defaulting to 1.'
!$    end if
#endif
      mpi_string = 'MPI'
!$    env_threads = omp_get_max_threads()
!$    if (pub_output_detail > NORMAL) then
!$      write(stdout,'(6x,a,i0,a)') 'OpenMP library call suggests ', env_threads, &
!$           ' threads.'
#ifdef MPI
!$      write(stdout,'(6X,a)') 'Hybrid ' // mpi_string // '-OpenMP compilation detected.'
!$    end if
!$      if (pub_total_num_procs > 1) then
!$        if (pub_output_detail > NORMAL) then
!$          write(stdout,'(6X,a)') 'Number of ' // mpi_string // ' processes > 1.'
!$        end if
!$        call get_environment_variable('OMP_NUM_THREADS', value=omp_num_threads, &
!$             status=ierr)
!$        if (ierr == 0) then
!$          read(omp_num_threads, *, iostat=ierr) omp_num_threads_int
!$          if ((ierr == 0) .AND. (omp_num_threads_int == env_threads)) then
!$            default_threads = env_threads
!$            if (pub_output_detail > NORMAL) then
!$              write(stdout,'(6x,a)') 'OMP_NUM_THREADS environment variable present and'
!$              write(stdout,'(6x,a)') 'its value matches the one obtained from &
!$                   &the OpenMP library call.'
!$              write(stdout,'(6x,a)') 'Using this value as default.'
!$            end if
!$          else
!$            if (pub_output_detail >= NORMAL) then
!$              write(stdout,'(6x,a)') 'WARNING: Error or mismatch comparing the value &
!$                   &obtained via the library call'
!$              write(stdout,'(6x,a)') 'with the OMP_NUM_THREADS environment variable.'
!$              write(stdout,'(6x,a)') 'Not changing the default number of threads &
!$                   &(risk of CPU oversubscription).'
!$            end if
!$          end if
!$        else
!$          if (pub_output_detail >= NORMAL) then
!$            write(stdout,'(6x,a)') 'WARNING: Error getting the value of the &
!$                 &OMP_NUM_THREADS environment variable.'
!$            write(stdout,'(6x,a)') 'Not changing the default number of threads &
!$                 &(risk of CPU oversubscription).'
!$          end if
!$        end if
!$      else
!$        default_threads = env_threads
!$        if (pub_output_detail > NORMAL) then
!$          write(stdout,'(6x,a)') 'Only 1 process detected: using this value as default.'
!$        end if
!$      end if
#else
!$      write(stdout,'(6x,a)') 'Pure OpenMP compilation: using this value as default.'
!$    end if
!$    default_threads = env_threads
#endif

!$    if (pub_output_detail >= NORMAL) then
!$      write(stdout,'(6x,a,i0)') 'Default threads: ', default_threads
!$    end if
!$    if (pub_output_detail > NORMAL) then
!$      write(stdout,'(6x,a)') 'Default applies to threads_max' //  &
!$             ' and threads_num_fftboxes.'
!$      write(stdout,'(6x,a)') 'Settings in the input file override this default.'
!$      write(stdout,'(3x,a/)') '... done'
!$    end if

    end function internal_default_threads

  end subroutine rundat_threads_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rundat_check_inputs

    use constants, only: DP, stdout
    use comms, only: comms_bcast, pub_on_root, pub_comms_group_size, &
         pub_total_num_procs, pub_root_proc_id, pub_root_process_count
    use utils, only: utils_assert, utils_abort, utils_parse_bc

    implicit none

    ! Local Variables
    logical :: inconsistent_bcs
    logical :: mg_full_obc_hartree, mg_full_pbc_hartree, mg_mixed_bc_hartree
    integer :: isize

    ! -------------------------------------------------------------------------

    ! rc2013: if we're reading in the kernel from subsystem calculations,
    ! check that the arguments are consistent
    if(.not. pub_read_denskern .and. pub_read_sub_denskern) then
         write(stdout,'(/a,a,a/)') 'WARNING: read_denskern overridden to T&
              & because read_sub_denskern = T in input file.'
         pub_read_denskern = .true.
    end if

    ! ndmh: Logic for setting up kernel initialisation iteration
    ! counts, to avoid spoiling denskern in single-point energy restarts
    if (pub_task == 'SINGLEPOINT') then
       if(pub_read_sub_denskern) then
          !ep: mermin method switch
          if (pub_mermin) then
             pub_firstit_mermin = maxit_mermin
          else
             pub_firstit_lnv = minit_lnv
          end if
       elseif (pub_read_denskern .and. pub_read_tightbox_ngwfs) then
          !ep: mermin method switch
          if (pub_mermin) then
             pub_firstit_mermin = maxit_mermin
          else
             pub_firstit_lnv = 0
          end if
          pub_maxit_pen = 0
       elseif(pub_read_denskern .and. pub_read_sw_ngwfs) then ! ars
          !ep: mermin method switch
          if (pub_mermin) then
             pub_firstit_mermin = maxit_mermin
          else
             pub_firstit_lnv = 0
          end if
          pub_maxit_pen = 0                                   ! ars
          ! ddor: Bare response calculation for DFT+U parameters
       elseif (pub_hub_max_iter .eq. -1) then
          !ep: mermin method switch
          if (pub_mermin) then
             pub_firstit_mermin = maxit_mermin
          else
             pub_firstit_lnv = 0
          end if
          pub_maxit_pen = 0
       else
          !ep: mermin method switch
          if (pub_mermin) then
             pub_firstit_mermin = maxit_mermin
          else
             pub_firstit_lnv = minit_lnv
          end if
       endif
    else
       !ep: mermin method switch
       if (pub_mermin) then
          pub_firstit_mermin = maxit_mermin
       else
          pub_firstit_lnv = minit_lnv
       end if
    endif

    ! ndmh: Logic for setting up kernel initialisation iteration
    ! counts, to avoid spoiling denskern in COND restarts
    if (pub_task == 'COND') then
       if (cond_read_denskern .and. cond_read_tightbox_ngwfs) then
          cond_firstit_lnv = 0
       elseif(pub_read_denskern .and. pub_read_sw_ngwfs) then ! ars
          cond_firstit_lnv = 0
          ! ddor: Bare response calculation for DFT+U parameters
       elseif (pub_hub_max_iter .eq. -1) then
          cond_firstit_lnv = 0
       else
          cond_firstit_lnv = cond_minit_lnv
       endif
       ! Cannot use EDFT initialisation for COND pub_task
       if (pub_edft_init_maxit>0) then
          if (pub_on_root) write(stdout,'(a)') 'WARNING in rundat_check_inputs: &
               &edft_init_maxit cannot be used for COND calculations'
          if (pub_on_root) write(stdout,'(a)') 'WARNING in rundat_check_inputs: &
               &Setting edft_init_maxit = 0'
          pub_edft_init_maxit = 0
       end if
    else
       cond_firstit_lnv = cond_minit_lnv
    endif

    ! smmd : impose consistency between mix_ngwfs_type,
    ! smmd : mix_ngwfs_num and mix_ngwfs_init_num
    call internal_check_mix(mix_ngwfs_type, mix_ngwfs_num, &
         mix_ngwfs_init_type, mix_ngwfs_init_num, mix_ngwfs_reset)
    ! vv: Check for Extended-Lagrangian with dissipation for the NGWFs
    if (mix_ngwfs_type=='XLD') then
       if (pub_on_root) write(stdout,'(a)') 'WARNING in rundat_check_inputs: &
            &mix_ngwfs_type cannot be set to XLD. This functionality is no longer available'
       if (pub_on_root) write(stdout,'(a)') 'WARNING in rundat_check_inputs: &
            &Setting mix_ngwfs_type to NONE '
       mix_ngwfs_type = 'NONE'
       mix_ngwfs_init_num = 0
       call internal_check_mix(mix_ngwfs_type, mix_ngwfs_num, &
          mix_ngwfs_init_type, mix_ngwfs_init_num, mix_ngwfs_reset)
    end if

    ! smmd : impose consistency between mix_dkn_type,
    ! smmd : mix_dkn_num and mix_dkn_init_num
    call internal_check_mix(mix_dkn_type, mix_dkn_num, &
         mix_dkn_init_type, mix_dkn_init_num, mix_dkn_reset, &
         mix_ngwfs_num)

    ! vv: Check md_restart is false when md_global_restart is true.
    if(pub_task == 'MOLECULARDYNAMICS' .and. md_global_restart) then
       pub_read_denskern = .false.
       pub_read_tightbox_ngwfs = .false.
       md_restart = .false.
       md_write_out = .false.
    end if

    if (pub_task == 'MOLECULARDYNAMICS' .and. md_write_history .gt. 0 .and. &
        (.not. md_global_restart) ) then
       if (md_write_history < max(mix_ngwfs_init_num,mix_dkn_init_num)) then
           md_write_history = max(md_write_history,&
                max(mix_ngwfs_init_num,mix_dkn_init_num))
       end if
    end if

    ! vv: Check the type of the auxiliary dofs
    if (pub_task == 'MOLECULARDYNAMICS' .and. &
        (mix_dkn_type=='XLI' .or. mix_dkn_type=='XLD' .or. &
         mix_dkn_type=='NAIVE') .or. mix_dkn_type=='XLIS') then
         call utils_assert(md_aux_rep=='ASYM' .or. md_aux_rep=='ORTHO', &
              'Error in get_rundat: &
             &for a molecular dynamics calculation with density kernel propagation, &
             &md_aux_rep must be set to either ASYM or ORTHO')
    end if

    ! lr408: Set value of pub_cond_calculate according to pub_task
    ! lr408: and check parameters are sensible
    if (pub_task == 'COND' .or. pub_task == 'PROPERTIES_COND') then
       pub_cond_calculate = .true.
       if (.not. cond_fixed_shift) then
          if (.not. cond_calc_max_eigen) then
             if (pub_on_root) write(stdout,*) 'Error in get_rundat: &
                  &for a conduction calculation with updating shift'
             if (pub_on_root) write(stdout,*) 'cond_calc_max_eigen must be set to true'
             cond_calc_max_eigen = .true.
          end if
       end if
    else
       ! dhpt: Coupling code breaks otherwise if reading
       ! dhpt: TDDFT states
       if (pub_task /= 'COUPLINGS') then
          pub_cond_calculate = .false.
       endif
    end if

    ! ndmh: setup value of comms_group_size
    ! kaw: modified slightly to get default behaviour if an incompatible value
    ! is present in the input file.
    if (modulo(pub_total_num_procs,pub_comms_group_size)/=0) then
       if (pub_on_root) write(stdout,*) 'WARNING in get_rundat: &
            &Comms group size is incompatible with number of procs'
       if (pub_on_root) write(stdout,*) 'WARNING in get_rundat: &
            &Overriding to comms_group_size : -1'
       pub_comms_group_size = -1
    end if
    if (pub_comms_group_size<0) then
       do isize=int(sqrt(real(pub_total_num_procs,kind=DP))),1,-1
          if (int(pub_total_num_procs/isize)*isize==pub_total_num_procs) then
             pub_comms_group_size = isize
             exit
          end if
       end do
       ! Better method, which should be used if possible: pick number
       ! of processes on proc holding root process
       if (modulo(pub_total_num_procs,pub_root_process_count)==0) then
          pub_comms_group_size = pub_root_process_count
       end if
    end if
    if (pub_total_num_procs/pub_comms_group_size>100) then
       if (pub_on_root) then
          write(stdout,'(a)') 'WARNING in get_rundat: &
               &number of comms groups greater than 100'
          write(stdout,'(a)') '                   : &
               &may result in serious comms inefficiency'
       end if
    end if

    ! aam: checks
    if ( pub_task=='MOLECULARDYNAMICS' ) then
       if (md_delta_t.le.0.0_dp) then
          if (pub_on_root) write(stdout,*) 'Error in get_rundat: &
               &md_delta_t must be positive'
       endif
       if (md_num_iter.lt.0) then
          if (pub_on_root) write(stdout,*) 'Error in get_rundat: &
               &md_num_iter must be positive'
       endif
    endif

    ! pdh: sanity checks
    call utils_assert(.not. pub_old_lnv .or. .not. pub_exact_lnv, &
       'Error in get_rundat: setting both old_lnv and exact_lnv is forbidden')

    ! aam: more sanity checking
    call utils_assert(pub_geom_convergence_win >= 2, &
         'Error in get_rundat: geom_convergence_win < 2')

    if (pub_task=='GEOMETRYOPTIMIZATION') then
       if (pub_geom_reuse_dk_ngwfs.and.any(pub_extend_ngwf)) then
          if (pub_on_root) then
             write(stdout,'(/a,a,a/)') 'WARNING: extend_ngwf is true, so &
                  & geom_reuse_dk_ngwfs cannot be true.'
             write(stdout,'(/a,a,a/)') 'WARNING: setting geom_reuse_dk_ngwfs &
                  & to false'
          end if
          pub_geom_reuse_dk_ngwfs = .false.
       end if
    end if

    if (pub_task=='GEOMETRYOPTIMIZATION' .and. &
         (pub_geom_max_iter < pub_geom_convergence_win)) then
       call utils_abort('Error in get_rundat: &
            &geom_max_iter < geom_convergence_win')
    end if

    ! lam81
    if (pub_task=='GEOMETRYOPTIMIZATION') then
       select case(pub_geom_precond_type)
       case ('NONE')
       case ('ID')
          if (.not. trim(adjustl(pub_geom_method))=='LBFGS') then
             call utils_abort('ID preconditioner is currently implemented only &
                              &for LBFGS optimizer')
          end if
       case ('EXP')
          if (.not. trim(adjustl(pub_geom_method))=='LBFGS') then
             call utils_abort('EXP preconditioner is currently implemented only &
                              &for LBFGS optimizer')
          end if
       case ('FF')
          if (.not. trim(adjustl(pub_geom_method))=='LBFGS') then
             call utils_abort('FF preconditioner is currently implemented only &
                              &for LBFGS optimizer')
          end if
       end select
    end if

    ! cks: automatic setup for properties task
    ! ddor: included TDDFT task
    ! lr408: Included conduction task
    ! ndmh: included properties_cond task
    if (.not. pub_anharmonic_calculate .and. pub_task == 'PROPERTIES'&
         .or. pub_task == 'TDDFT' .or. &
         pub_task == 'COND' .or. pub_task == 'PROPERTIES_COND') then
       ! cks: force ngwf and denskern read from file
       ! jd: Warn the user that this is happening
       if(.not. pub_read_denskern .and. pub_on_root) then
          write(stdout,'(/a,a,a/)') 'WARNING: read_denskern overridden to T&
               & because of task ',trim(pub_task),'.'
       end if
       pub_read_denskern =.true.

       ! ars: ONETEP uses tightbox_ngwfs unless we specify in the input
       ! ars: read_tightbox_ngwfs = FALSE;  pub_read_sw_ngwfs = TRUE
       if (pub_read_tightbox_ngwfs) then
          pub_read_sw_ngwfs =.false.
       elseif (.not. pub_read_sw_ngwfs) then
          ! jd: Warn the user that this is happening
          if(pub_on_root) then
             write(stdout,'(/a,a,a/)') 'WARNING: read_tightbox_ngwfs overridden&
                  & to T because of task ',trim(pub_task),'.'
          end if
          pub_read_tightbox_ngwfs =.true.       ! ars
       endif                                ! ars
       pub_maxit_pen =0
    endif

    ! vv: reset read_denskern and read_tightbox to false for anharmonic
    ! vv: calculations
    if(pub_task == 'PROPERTIES' .and. pub_anharmonic_calculate) then
      if( pub_read_denskern .and. pub_read_tightbox_ngwfs) then
          if(pub_on_root) then
          write(stdout,'(/a,a,a/)') 'WARNING: read_denskern overridden to F&
               & because of anharmonic_calculate = .true.'
             write(stdout,'(/a,a,a/)') 'WARNING: read_tightbox_ngwfs overridden&
                  & to F because of anharmonic_calculate = .true.'
          end if
        pub_read_denskern = .false.
        pub_read_tightbox_ngwfs = .false.
      end if
    end if

    ! ja531-> sanity check PDOS input
    call utils_assert(.not.(pub_pdos_pseudoatomic.and.pub_pdos_reduce_sws), &
         'Both pdos_pseudoatomic and pdos_reduce_sws enabled - this is incompatible.')


    ! ars: set flag for kernel DIIS
    ! ebl: exception for DMFT
    if (pub_dmft_points /= 0) then
       pub_kernel_diis_scheme = 'DKN_PULAY'
       pub_kernel_diis = .false.
    else if (pub_kernel_diis_scheme.eq.'NONE') then
       pub_kernel_diis = .false.
    else
       pub_kernel_diis = .true.
       ! ars: note: if the scheme is nonsense the code will stop when entering
       !            kernel_diis_calculate
    end if

    ! ars: check that EDFT does not interfere with kernel DIIS
    call utils_assert(.not. pub_kernel_diis .or. .not. pub_edft, &
         'Both EDFT and kernel DIIS are enabled -- this is incompatible.')

    ! ars: set ngwf_cg_rotate to true for EDFT calculations
    if (pub_edft.and.(.not.pub_ngwf_cg_rotate)) then
        pub_ngwf_cg_rotate = .true.
        if (pub_on_root) write(stdout,'(/a)') &
           "Resetting ngwf_cg_rotate to TRUE for EDFT calculations"
    end if

    if (pub_task == 'PROPERTIES_COND') then
       ! ndmh: force cond ngwf and denskern read from file
       if(.not. cond_read_denskern .and. pub_on_root) then
          write(stdout,'(/a,a,a/)') 'WARNING: cond_read_denskern overridden &
               &to T because of task ',trim(pub_task),'.'
       end if
       cond_read_denskern =.true.

       if(pub_on_root) then
          write(stdout,'(/a,a,a/)') 'WARNING: cond_read_tightbox_ngwfs &
               &overridden to T because of task ',trim(pub_task),'.'
       end if
       cond_read_tightbox_ngwfs = .true.
    endif

    ! tjz07: Setup for LR_TDDFT calculation
    ! tjz07: Force reading of denskern + NGWFs
    if (pub_task == 'LR_TDDFT') then
       pub_lr_tddft_calculate = .true.

       ! need to set this to true in order to read and define the conduction
       ! NGWFs correctly
       pub_cond_calculate = .true.

       if (.not. cond_fixed_shift) then
          if (.not. cond_calc_max_eigen) then
             if (pub_on_root) write(stdout,*) 'Error in get_rundat: &
                  &for a conduction calculation with updating shift'
             if (pub_on_root) write(stdout,*) 'cond_calc_max_eigen must be set &
                  &to true'
             cond_calc_max_eigen = .true.
          end if
       end if

       ! overwrite read_denskern and write a warning
       if(.not. pub_read_denskern .and. pub_on_root) then
          write(stdout,'(/a,a,a/)') 'WARNING: read_denskern overridden to T &
               &because of task ',trim(pub_task),'.'
       end if
       pub_read_denskern = .true.

       !if(.not. cond_read_denskern .and. pub_on_root) then
       !   write(stdout,'(/a,a,a/)') 'WARNING: cond_read_denskern overridden &
       !        &to because of task ',trim(pub_task),'.'
       !end if
       !cond_read_denskern = .true.

       ! tjz07: Force reading of both conduction + valence NGWFs
       if (.not. pub_read_tightbox_ngwfs .and. pub_on_root) then
          write(stdout,'(/a,a,a/)') 'WARNING: read_tightbox_ngwfs overridden &
               &to T because of task ',trim(pub_task),'.'
       end if
       pub_read_tightbox_ngwfs=.true.

       if (.not. cond_read_tightbox_ngwfs .and. pub_on_root) then
          write(stdout,'(/a,a,a/)') 'WARNING: cond_read_tightbox_ngwfs &
               &overridden to T because of task ',trim(pub_task),'.'
       end if
       cond_read_tightbox_ngwfs = .true.
       pub_maxit_pen = 0

       ! sanity check regarding projector use
       if(.not. pub_lr_tddft_projector .and. pub_lr_tddft_joint_set) then
          pub_lr_tddft_projector=.true.
          write(stdout,'(/a,a,a/)') 'WARNING: lr_tddft_projector &
               &overridden to T. Cannot use joint basis set without projector.'
       endif

       ! sanity check regarding projector use
       if(pub_lr_tddft_projector .and. .not. pub_lr_tddft_joint_set) then
          pub_lr_tddft_projector =.false.
          write(stdout,'(/a,a,a/)') 'WARNING: lr_tddft_projector overridden &
               &to F. Cannot use projector without the joint basis set.'
       endif

       ! sanity check regarding restarts:
       if(pub_lr_tddft_restart .and. pub_lr_tddft_preopt) then
          pub_lr_tddft_preopt=.false.
          write(stdout,'(/a,a,a/)') 'WARNING: lr_tddft_preopt overridden &
               &to F. When restarting from a previous calculation, &
               &a preoptimisation step is not appropriate.'
       endif

       ! sanity check regarding triplets. If this is an open-shell
       ! calculation the triplet keyword should be overwritten to false
       ! as the separation between singlet and triplet subspaces is only
       ! possible for spin-degenerate systems
       if(pub_num_spins==2 .and. pub_lr_tddft_triplet) then
          pub_lr_tddft_triplet=.false.
          write(stdout,'(/a,a,a/)') 'WARNING: lr_tddft_triplet overridden &
               &to F. This is an open-shell TDDFT calculation and the &
               &singlet and triplet subspaces cannot be separated.'
       endif

       ! sanity check regarding number of NTO transitions: !jhl52
       if(pub_qnto_num_transition > 5) then
          pub_qnto_num_transition = 5
          if(pub_on_root) &
               write(stdout,*) 'Number of NTO transitions is reset to 5.'
       endif

    endif

    ! ndmh: when doing auto-solvation, when doing COND or LR-TDDFT pub_task,
    ! ndmh: the vacuum NGWFs and kernel already exist, so there is no
    ! ndmh: need to run the vacuum calculation. However, the vacuum results
    ! ndmh: must be loaded to generate the correct cavity, so set
    ! ndmh: pub_is_separate_restart_files to true.
    if (pub_is_auto_solvation .and. &
         (pub_cond_calculate.or.pub_lr_tddft_calculate)) then
       pub_is_separate_restart_files = .true.
       if (pub_on_root) write(stdout,'(/a,a,a/)') &
            'WARNING: pub_is_separate_restart_files has &
            &been set to true due to is_auto_solvation = true.'
    end if

    ! jcap: check that if pub_is_emft_cavity is set to true, pub_emft
    ! or pub_emft_follow is also set to true -- otherwise the required
    ! dkn file may not be read/exist
    if (pub_is_emft_cavity .and. .not.(pub_emft .or. pub_emft_follow)) &
         call utils_abort("is_emft_cavity set to T without EMFT being requested")

    ! mjsp: EDA input check
    if ( pub_task=='EDA' .or. pub_task=='EDA_PREP') then
       ! mjsp: sanity checks:
       if (pub_hubbard) call utils_abort("DFT+U not implemented &
            &in energy decomposition analysis or derived tasks.")

       call utils_assert((pub_kernel_diis .eqv. .false.),'Energy decomposition &
             &analysis is not available with DIIS.')
       if (pub_task == 'EDA' .or. pub_task == 'EDA_PREP') &
         call utils_assert((pub_dispersion == '0' .or. pub_dispersion == 'NONE'), &
           'Dispersion is not yet implemented with EDA.')

       pub_eda = .true.
       if (pub_task ==  'EDA_PREP') then
          pub_eda_preptool = .true.
       else
          pub_eda_preptool = .false.
       end if

    else
       ! mjsp: by default do not carry out EDA
       pub_eda = .false.
       pub_eda_preptool = .false.
    endif


    ! gcc32: Setup for LR_PHONONS calculation
    ! gcc32: Force reading of denskern + NGWFs
    if (pub_task == 'LR_PHONONS') then
       pub_lr_phonons_calculate = .true.


       ! overwrite denskernel-read and write a warning
       if (.not. pub_read_denskern .and. pub_on_root) then
          write(stdout,'(/a,a,a/)') 'WARNING: read_denskern overridden to T &
               &because of task ',trim(pub_task),'.'
       end if

       pub_read_denskern = .true.
       if (.not. pub_read_tightbox_ngwfs .and. pub_on_root) then
          write(stdout,'(/a,a,a/)') 'WARNING: read_tightbox_ngwfs overridden &
               &to T because of task ',trim(pub_task),'.'
       end if
       pub_read_tightbox_ngwfs=.true.

       pub_maxit_pen = 0

    endif

    ! ndmh: set boolean flag for cutoff coulomb
    if ((index(pub_coulomb_cutoff_type,'NONE')>0).or. &
         (index(pub_coulomb_cutoff_type,'none')>0))then
       pub_coulomb_cutoff=.false.
    else
       pub_coulomb_cutoff=.true.
    end if

    ! ndmh: sanity check for Cutoff Coulomb
    ! jd: More detailed checks introduced so that users do not forget or
    !     confuse radius and length parameters.
    if (pub_coulomb_cutoff) then
       if(trim(pub_coulomb_cutoff_type) == 'SPHERE' .or. &
            trim(pub_coulomb_cutoff_type) == 'sphere') then
          call utils_assert(pub_coulomb_radius /= -1.0_DP, "'coulomb_cutoff_rad&
               &ius' must be specified when using 'coulomb_cutoff_type SPHERE'")
          call utils_assert(pub_coulomb_radius >= 0.0_DP, &
               'coulomb_cutoff_radius must be positive')
       elseif(trim(pub_coulomb_cutoff_type) == 'WIRE' .or. &
            trim(pub_coulomb_cutoff_type) == 'wire') then
          call utils_assert(pub_coulomb_radius /= -1.0_DP, "'coulomb_cutoff_rad&
               &ius' must be specified when using 'coulomb_cutoff_type WIRE'")
          call utils_assert(pub_coulomb_radius >= 0.0_DP, &
               'coulomb_cutoff_radius must be positive')
       elseif(trim(pub_coulomb_cutoff_type) == 'SLAB' .or. &
            trim(pub_coulomb_cutoff_type) == 'slab') then
          call utils_assert(pub_coulomb_length /= -1.0_DP, "'coulomb_cutoff_len&
               &gth' must be specified when using 'coulomb_cutoff_type SLAB'")
          call utils_assert(pub_coulomb_length >= 0.0_DP, &
               'coulomb_cutoff_length must be positive')
       elseif(trim(pub_coulomb_cutoff_type) == 'CYLINDER' .or. &
            trim(pub_coulomb_cutoff_type) == 'cylinder') then
          call utils_assert(pub_coulomb_radius /= -1.0_DP .and. &
               pub_coulomb_length /= -1.0_DP, "'coulomb_cutoff_len&
               &gth' *and* 'coulomb_cutoff_radius' must be specified when using&
               & 'coulomb_cutoff_type CYLINDER'")
          call utils_assert(pub_coulomb_radius >= 0.0_DP, &
               'coulomb_cutoff_radius must be positive')
          call utils_assert(pub_coulomb_length >= 0.0_DP, &
               'coulomb_cutoff_length must be positive')
       else
          call utils_abort("Unrecognized string for coulomb_cutoff_type. Allowe&
              &d values are: 'SPHERE', 'WIRE', 'SLAB', 'CYLINDER' and 'NONE'")
       endif
    end if

    ! qoh: Temporary assignment of tightbox FFT variables
    pub_tightbox_fft_coarse = ( pub_write_sw_ngwfs .or. pub_read_sw_ngwfs )

    ! ddor: Hubbard DFT+U input check
    if ( pub_task == 'HUBBARDSCF' ) then
       pub_hubbard = .true.
       pub_delta_e_conv = .false.
       call utils_assert(.not. pub_nnho, 'Error in get_rundat: pub_nnho==.true.&
            & is incompatible with projector-self-consistent DFT+U.')
       call utils_assert(pub_hub_proj_mixing < 1.0_DP, &
          'Error in get_rundat: pub_hub_proj_mixing out of bounds.')
       if (pub_hub_on_the_fly) call utils_assert(pub_hub_on_the_fly, &
          'Error in get_rundat: HUBBARDSCF + on-the-fly are incompatible.')
    endif
    pub_hubbard_restart = .false.
    !gibo: modified for hubbard/cdft_read_projectors:T restarts
    if ((pub_hub_proj_mixing .lt. 0.0_DP) .OR. pub_hub_read_projectors .OR. &
         pub_cdft_read_projectors)  pub_hubbard_restart = .true.

    !gibo: if hub/cdft_read_projectors:T, make sure (pub_hub_proj_mixing <0)
    !      to read proj. from file withouth modifying "hubbard_projector_update"
    if (pub_hub_read_projectors .OR. pub_cdft_read_projectors) then
        if (pub_hub_proj_mixing .eq. 0.0_DP) then
           pub_hub_proj_mixing = -1.0_DP
        elseif (pub_hub_proj_mixing .gt. 0.0_DP) then
           pub_hub_proj_mixing = -pub_hub_proj_mixing
        endif
    endif
    !gibo: make sure forgotten HUBBARD/cDFT-restart keywords do not mess with
    !      each other
    if (pub_hub_read_projectors .AND. pub_cdft_read_projectors) then
       call utils_abort('Error in get_rundat: &
            &HUBBARD_READ_PROJECTORS .AND. CDFT_READ_PROJECTORS are NOT &
            &compatible with each other. Comment either of them.')
    endif

    ! gom : Complain if tensor_corr=4 but we have not unified sites in
    !       Hubbard, cdft+U or DFT+nu cases
    if((pub_hub_tensor_corr .eq. 4) .and. (.not. pub_hubbard_unify_sites) &
        & .and. (pub_hubbard_unify_sites .or. pub_cdft_hubbard .or. &
        & pub_dft_nu)) then
       write(stdout,'(a)') 'WARNING: hubbard_tensor_corr = 4 but &
       &hubbard_unify_sites : F. Hamiltonian is non-Hermitian.'
   endif

    ! qoh: pub_hfx_nlpp_for_exchange isn't compatible with Hubbard
    pub_hfx_nlpp_for_exchange = ( pub_hfx_nlpp_for_exchange .and. &
         (.not. pub_hubbard) )
    if (pub_hfx_nlpp_for_exchange .and. pub_hubbard .and. pub_on_root) then
       write(stdout,'(/a)') 'WARNING: pub_hfx_nlpp_for_exchange set to FALSE, &
            &because it''s not compatible with HUBBARDSCF.'
    end if

    ! jd: Sanity checks for HFx and SWX
    call utils_assert(pub_hfx_bessel_rad_nptsx > 1, &
         'hfx_bessel_rad_nptsx must be be >1')
    call utils_assert(pub_swri_cheb_batchsize >= 0, &
         'swri_cheb_batchsize must be non-negative.')
    call utils_assert(pub_cache_limit_for_swops >= 0, &
         'cache_limit_for_swops must be non-negative.')
    call utils_assert(pub_cache_limit_for_swops2 >= 0, &
         'cache_limit_for_swops2 must be non-negative.')
    call utils_assert(pub_cache_limit_for_expansions >= 0, &
         'cache_limit_for_expansions must be non-negative.')
    call utils_assert(pub_cache_limit_for_prods >= 0, &
         'cache_limit_for_prods must be non-negative.')
    call utils_assert(pub_hfx_memory_limit >= -1, &
         'hfx_memory_limit must be non-negative or -1.')

    ! gibo: CONSTRAINED_DFT checks
    ! gibo:******************* CONSTRAINED_DFT checks ****************** START

    ! gibo: avoid activation of conflicting atom DFT-modes ==== START
    if ( pub_cdft_atom_charge .AND. pub_cdft_atom_spin) then
       call utils_abort('Error in get_rundat: &
            &cdft_atom_charge .AND. cdft_atom_spin are NOT compatible with each &
            &other.')
    endif

    if ( (pub_cdft_atom_charge .OR. pub_cdft_atom_spin) .AND. &
          (pub_cdft_group_charge_acceptor .OR. pub_cdft_group_charge_donor .OR.&
           pub_cdft_group_spin_acceptor .OR. pub_cdft_group_spin_donor.OR.     &
           pub_cdft_group_charge_diff.OR.pub_cdft_group_spin_diff) ) then
       call utils_abort('Error in get_rundat: &
            &cdft_atom_charge .AND. cdft_atom_spin are NOT compatible with other &
            &cdft_GROUP_CHARGE/SPIN[_DIFF] constrained-DFT modes.')
    end if

    if ((pub_cdft_group_charge_acceptor .OR. pub_cdft_group_charge_donor .OR.&
          pub_cdft_group_spin_acceptor .OR. pub_cdft_group_spin_donor).AND.&
           (pub_cdft_group_charge_diff.OR.pub_cdft_group_spin_diff) )  then
       call utils_abort('Error in get_rundat: &
            &cdft_GROUP_CHARGE/SPIN modes are NOT compatible with &
            &cdft_GROUP_CHARGE/SPIN_DIFF modes.')
    end if

    if ( (pub_cdft_group_charge_up_only .OR. pub_cdft_group_charge_down_only) .AND. &
          (pub_cdft_group_spin_diff.OR.pub_cdft_group_spin_acceptor .OR. &
           pub_cdft_group_spin_donor) ) then
       call utils_abort( 'Error in get_rundat: &
            &cdft_GROUP_SPIN_ACCEPTOR/DONOR/DIFF modes are NOT compatible with &
            &cdft_GROUP_CHARGE_UP/DOWN_ONLY.')
    end if
    if ( pub_cdft_group_charge_up_only .AND.pub_cdft_group_charge_down_only) then
       call utils_abort('Error in get_rundat: &
            &cdft_GROUP_CHARGE_UP_ONLY and cdft_GROUP_CHARGE_DOWN_ONLY are NOT &
            &compatible.')
    endif
    ! gibo: avoid simultanous activation of conflicting cDFT-modes ==== END

    ! gibo: force the user to set cdft_group_charge_acceptor=.T. if
    ! gibo: cdft_group_charge_acceptor_target is found
    if ((ABS(pub_cdft_group_charge_acceptor_target) > 0._DP) .AND. &
         (.not.pub_cdft_group_charge_acceptor))  then
       call utils_abort('Error in get_rundat: &
            &set pub_cdft_group_charge_acceptor=.T. for a group-charge-&
            &acceptor-constrained cDFT simulation')
    end if

    ! gibo: force the user to set cdft_group_charge_donor=.T. if
    ! gibo: cdft_group_charge_donor_target is found
    if ((ABS(pub_cdft_group_charge_donor_target) > 0._DP) .AND. &
         (.not.pub_cdft_group_charge_donor))  then
       call utils_abort('Error in get_rundat: &
            &set pub_cdft_group_charge_donor=.T. for a group-charge-&
            &donor-constrained cDFT simulation.')
    end if

    ! gibo: force the user to set cdft_group_spin_acceptor=.T. if
    ! gibo: cdft_group_spin_acceptor_target is found
    if ((ABS(pub_cdft_group_spin_acceptor_target) > 0._DP) .AND. &
         (.not.pub_cdft_group_spin_acceptor))  then
       call utils_abort('Error in get_rundat: &
            &set pub_cdft_group_spin_acceptor=.T. for a group-spin-&
            &acceptor-constrained cDFT simulation.')
    end if

    ! gibo: force the user to set cdft_group_spin_donor=.T. if
    ! gibo: cdft_group_spin_donor_target is found
    if ((ABS(pub_cdft_group_spin_donor_target) > 0._DP) .AND. &
         (.not.pub_cdft_group_spin_donor))  then
       call utils_abort('Error in get_rundat: &
            &set pub_cdft_group_spin_donor=.T. for a group-spin-&
            &donor-constrained cDFT simulation.')
    end if

    ! gibo: force the user to set cdft_group_charge_diff=.T. if
    ! gibo: cdft_group_charge_diff_target is found
    if ((ABS(pub_cdft_group_charge_diff_target) > 0._DP) .AND. &
         (.not.pub_cdft_group_charge_diff))  then
       call utils_abort('Error in get_rundat: &
            &set pub_cdft_group_charge_diff=.T. for a group-charge-&
            &difference-constrained cDFT simulation.')
    end if

    ! gibo: force the user to set cdft_group_spin_diff=.T. if
    ! gibo: cdft_group_spin_diff_target is found
    if ((ABS(pub_cdft_group_spin_diff_target) > 0._DP) .AND. &
         (.not.pub_cdft_group_spin_diff))  then
       call utils_abort('Error in get_rundat: &
            & set pub_cdft_group_spin_diff=.T. for a group-spin-&
            &difference-constrained cDFT simulation.')
    end if

    !force experienced user [cdft_guru=T[ to think of the best CG-frequency reset
    if (pub_cdft_guru.AND.(pub_cdft_cg_max==0)) then
       call utils_abort('Error in get_rundat: &
            &As cDFT_guru [cdft_guru=T], you need to specify cdft_cg_max in &
            &your input file.')
    endif
    ! gibo:******************* CONSTRAINED_DFT checks ****************** END

    !gom: DFT+nu sanity checks
    if ((pub_dft_nu_opt_u1_only) .AND. (pub_dft_nu_opt_u2_only)) then
       call utils_abort('Error in get_rundat: &
            & dft_nu_opt_u1_only and dft_nu_opt_u2_only cannot be &
            &simultanously .TRUE. in input file.')
    endif


    ! jd: Sanity check for pbc_correction_cutoff
    if (pub_mt_cutoff /= 0.0_DP) then

       call utils_assert(pub_mt_cutoff >= 0.0_DP, 'Error in get_rundat: &
               &pbc_correction_cutoff must be non-negative.')

       call utils_assert(.not. pub_coulomb_cutoff, 'Error in get_rundat: &
               &Cannot use cutoff Coulomb and pbc_correction_cutoff &
               &simultaneously.')

    end if

    ! jd: Sanity check for dx_format_digits
    call utils_assert(pub_dx_sig_digits >= 1, 'Error in get_rundat: &
            &number of significant digits specified with dx_format_digits &
            &must be positive.')

    ! jd: Ion-ion energy is calculated by direct summation iff
    !     (MT correction in use or cutoff Coulomb in use or smeared ions in use
    !      or implicit solvent in use or user overrode by pub_openbc_ion_ion.)
    pub_ii_energy_direct = &
         (pub_mt_cutoff /= 0.0_DP) .or. pub_coulomb_cutoff &
         .or. pub_is_smeared_ion_rep .or. pub_is_implicit_solvent &
         .or. pub_openbc_ion_ion

    ! jd: Hartree is calculated with multigrid iff
    !     (smeared ions in use or implicit solvent in use or user overrode by
    !     pub_openbc_hartree)
    pub_multigrid_hartree = &
         pub_is_smeared_ion_rep .or. pub_is_implicit_solvent &
         .or. pub_openbc_hartree

    ! jd: Local pseudo is calculated with open boundary conditions iff
    !     (smeared ions in use or implicit solvent in use or user overrode by
    !     pub_openbc_pspot)
    pub_open_localpseudo = &
         pub_is_smeared_ion_rep .or. pub_is_implicit_solvent &
         .or. pub_openbc_pspot

    ! JCW: Change to meaning of pub_multigrid_hartree, pub_openbc_hartree
    ! JCW: * pub_multigrid_hartree is .true. when multigrid is used to compute
    ! JCW:   Hartree / electrostatic potential, regardless of BCs (previously
    ! JCW:   this could be used as an indicator for open BCs)
    ! JCW: * pub_openbc_hartree is .true. based on user input, and must be
    ! JCW:   in agreement with any user-provided pub_multigrid_bc values
    ! JCW: * pub_openbc_hartree sets pub_multigrid_bc, if this was not explicitly
    ! JCW:   set by the user
    ! JCW: => pub_multigrid_bc should be used to check consistency of BCs, since
    ! JCW:    this will always be set: either the user explicitly sets this or
    ! JCW:    it is set based on pub_openbc_hartree

    ! JCW: Defining pub_multigrid_bc activates the multigrid solver
    if ( len_trim(pub_multigrid_bc) > 0) pub_multigrid_hartree = .true.
    ! => pub_multigrid_hartree is now .true. for use of the multigrid solver
    !    regardless of the BCs

    ! JCW: If pub_multigrid_bc is empty but pub_multigrid_hartree is .true.,
    ! JCW: then we assume fully open BCs (this does not change the behaviour
    ! JCW: with respect to earlier code versions)
    if (pub_multigrid_hartree.and.len_trim(pub_multigrid_bc) <= 0) then
       pub_multigrid_bc_is_periodic(1:3) = .false.
       pub_multigrid_bc_is_zero(1:3)     = .false.
    end if

    ! jd: Check consistency of boundary conditions
    call internal_check_bc

    ! jd: This flag tells cell_grid_distribute to ensure that the slab
    !     distribution is compatible with the multigrid solver. Set this flag
    !     to .true. anytime multigrid is used. Currently this is equivalent
    !     to using multigrid_hartree, but if there are other uses for the
    !     multigrid in the future (Simon?), this can be extended.
    pub_multigrid_in_use = pub_multigrid_hartree

    ! jd: Sanity check for implicit solvent keywords
    call utils_assert(pub_is_density_threshold > 0.0_DP, &
         'is_density_threshold must be positive')
    call utils_assert(pub_is_solvation_beta > 0.0_DP, &
         'is_solvation_beta must be positive')
    call utils_assert(pub_is_bulk_permittivity >= 1.0_DP, &
         'is_bulk_permittivity must not be below 1.')
    call utils_assert(pub_is_smeared_ion_width > 0.0_DP, &
         'is_smeared_ion_width must be positive')
    call utils_assert(pub_is_solvation_method == 'DIRECT' .or. &
         pub_is_solvation_method == 'CORRECTIVE', &
         'is_solvation_method must be either ''direct'' or ''corrective''.')
    call utils_assert(pub_is_dielectric_model == 'FIX_INITIAL' .or. &
         pub_is_dielectric_model == 'SELF_CONSISTENT', &
         'is_dielectric_model must be one of ''fix_initial'', &
         &'' or ''self_consistent''.')
    call utils_assert(pub_is_pbe == 'NONE' .or. pub_is_pbe == 'FULL' .or. &
         pub_is_pbe == 'LINEARISED', &
         'is_pbe must be one of ''none'', ''full'', or ''linearised''.')
    call utils_assert(.not. pub_is_pbe_bc_debye_screening .or. &
         (pub_is_pbe == 'FULL' .or. pub_is_pbe == 'LINEARISED'), &
         'It does not make sense to have is_bc_debye_screening with is_pbe ''none''.')
    call utils_assert(pub_is_apolar_sasa_definition == 'ISODENSITY' .or. &
         pub_is_apolar_sasa_definition == 'DENSITY', &
         'is_apolar_sasa_definition must be either ''density'', &
         &''isodensity''.')
    call utils_assert(pub_is_dielectric_function == 'FGF' .or. &
         pub_is_dielectric_function == 'GAUSSIAN' .or. &
         pub_is_dielectric_function == 'ANDREUSSI' .or. &
         pub_is_dielectric_function == 'SOFT_SPHERE', &
         'is_dielectric_function must be either ''fgf'', ''andreussi'', &
         &''gaussian'' or ''soft_sphere''.')
    call utils_assert(pub_is_surface_thickness > 0.0_DP, &
         'is_surface_thickness must be positive')
    call utils_assert(pub_is_solvent_surf_tension >= 0.0_DP, &
         'is_solvent_surf_tension must be non-negative')
    call utils_assert(pub_is_bc_coarseness >= 0, &
         'is_bc_coarseness must be non-negative')
    call utils_assert(pub_is_bc_surface_coarseness > 0, &
         'is_bc_surface_coarseness must be positive')
    call utils_assert(pub_is_bc_threshold >= 0, &
         'is_bc_threshold must be non-negative')
    call utils_assert(pub_is_soft_sphere_delta > 0, &
         'is_soft_sphere_delta must be greater than 0')
    call utils_assert(pub_is_soft_sphere_scale > 0, &
         'is_soft_sphere_scale must be greater than 0')
    call utils_assert(.not. (pub_is_smeared_ion_rep .and. pub_coulomb_cutoff), &
         'is_smeared_ion_rep cannot be used simultaneously with &
         &coulomb_cutoff_type other than "NONE"')
    call utils_assert(.not. (pub_is_smeared_ion_rep .and. pub_mt_cutoff>0.0D0),&
         'is_smeared_ion_rep cannot be used simultaneously with &
         &pbc_correction_cutoff')
    call utils_assert(.not. (pub_is_include_apolar .and. .not. &
         pub_is_implicit_solvent),'Cannot have is_include_apolar without &
         &is_implicit_solvent')
    call utils_assert(pub_is_apolar_method == 'SASA' .or. &
         pub_is_apolar_method == 'SAV', &
         'is_apolar_method must be one of ''SASA'', &
         &'' or ''SAV'' .')
    call utils_assert(pub_is_pbe == 'NONE' .or. &
         pub_is_pbe_temperature >= 0D0, &
         'is_pbe_temperature must be specified and be non-negative if &
         &is_pbe is not ''none''.')
    call utils_assert(pub_is_hc_steric_smearing >= 0D0 .and. &
         pub_is_hc_steric_smearing <= 5D-1, &
         'is_hc_steric_smearing must be non-negative and not bigger than 0.5')
    call utils_assert(pub_is_sc_steric_smoothing_alpha >= 0D0, &
         'is_sc_steric_smoothing_alpha must be non-negative')
    call utils_assert(pub_is_dielectric_exclusions_smearing >= 0.0_DP, &
         'is_dielectric_exclusions_smearing must be non-negative')

    ! jd: Sanity check for open BC keywords
    call utils_assert(pub_openbc_pspot_finetune_nptsx > 1, &
         'openbc_pspot_finetune_nptsx must be >1')
    call utils_assert(pub_openbc_pspot_finetune_alpha >= 0.0_DP, &
         'openbc_pspot_finetune_alpha must be non-negative')
    call utils_assert(pub_openbc_pspot_finetune_f > 0 .or. &
         pub_openbc_pspot_finetune_f == -1, &
         'openbc_pspot_finetune_f must be positive or -1')

    ! JCW: Sanity check for finite differences keywords
    call utils_assert(pub_finite_difference_order > 0 .and. &
         mod(pub_finite_difference_order,2) == 0, &
         "finite_difference_order must be positive and even")

    ! JCW: Sanity check for multigrid solver keywords (updated for DL_MG v2.0)
    call utils_assert(pub_mg_granularity_power > 0, &
         "mg_granularity_power must be > 0")
    call utils_assert(pub_mg_defco_fd_order > 0 .and. &
         mod(pub_mg_defco_fd_order,2) == 0, &
         "mg_defco_fd_order must be positive and even")
    call utils_assert(pub_mg_tol_res_rel >= 0.0_DP,&
         "mg_tol_res_rel must be >= 0.0")
    call utils_assert(pub_mg_tol_res_abs >= 0.0_DP,&
         "mg_tol_res_abs must be >= 0.0")
    call utils_assert(pub_mg_tol_pot_rel >= 0.0_DP,&
         "mg_tol_pot_rel must be >= 0.0")
    call utils_assert(pub_mg_tol_pot_abs >= 0.0_DP,&
         "mg_tol_pot_abs must be >= 0.0")
    call utils_assert(pub_mg_tol_vcyc_rel >= 0.0_DP,&
         "mg_tol_vcyc_rel must be >= 0.0")
    call utils_assert(pub_mg_tol_vcyc_abs >= 0.0_DP,&
         "mg_tol_vcyc_abs must be >= 0.0")
    call utils_assert(pub_mg_tol_newton_rel >= 0.0_DP,&
         "mg_tol_newton_rel must be >= 0.0")
    call utils_assert(pub_mg_tol_newton_abs >= 0.0_DP,&
         "mg_tol_newton_abs must be >= 0.0")
    call utils_assert(pub_mg_max_iters_defco > 0,&
         "mg_max_iters_defco must be > 0")
    call utils_assert(pub_mg_max_iters_vcycle > 0,&
         "mg_max_iters_vcycle must be > 0")
    call utils_assert(pub_mg_max_iters_newton > 0,&
         "mg_max_iters_newton must be > 0")
    call utils_assert(pub_mg_max_iters_cg > 0,&
         "mg_max_iters_cg must be > 0")
    call utils_assert(pub_mg_max_res_ratio > 0.0_DP,&
         "mg_max_res_ratio must be positive")
    call utils_assert(pub_mg_tol_cg_res_rel >= 0.0_DP,&
         "mg_tol_cg_res_rel must be >= 0.0")
    call utils_assert(pub_mg_tol_cg_res_abs >= 0.0_DP,&
         "mg_tol_cg_res_abs must be >= 0.0")

    ! jd: If is_pbe_conserve_ion_number is specified, leave it alone.
    !     But if not, adjust the default depending on MG BCs and
    !     whether we do Boltzmann solvation.
    if(pub_is_pbe_neutralisation_scheme == '*UNDEFINED*') then
       if(all(pub_multigrid_bc_is_periodic(:)) .and. pub_is_pbe /= 'NONE') then
          ! jd: For PBCs with electrolyte
          pub_is_pbe_neutralisation_scheme = 'COUNTERIONS_AUTO'
       else
          if(all(pub_multigrid_bc_is_periodic(:))) then
             ! jd: For PBCs with no electrolyte
             pub_is_pbe_neutralisation_scheme = 'JELLIUM'
          else
             ! jd: for OBCs
             pub_is_pbe_neutralisation_scheme = 'NONE'
          end if
       end if
    end if

    call utils_assert(pub_is_pbe_neutralisation_scheme == 'NONE' .or. &
         pub_is_pbe_neutralisation_scheme == 'JELLIUM' .or. &
         pub_is_pbe_neutralisation_scheme == 'ACCESSIBLE_JELLIUM' .or. &
         pub_is_pbe_neutralisation_scheme == 'COUNTERIONS_FIXED' .or. &
         pub_is_pbe_neutralisation_scheme == 'COUNTERIONS_AUTO' .or. &
         pub_is_pbe_neutralisation_scheme == 'COUNTERIONS_AUTO_LINEAR', &
         'is_pbe_neutralisation_scheme must be one of ''none'', ''jellium'', &
         &''accessible_jellium'', ''counterions_fixed'', ''counterions_auto'', &
         &or ''counterions_auto_linear''. You asked for '''//&
         trim(pub_is_pbe_neutralisation_scheme)//'''.')

    if(all(pub_multigrid_bc_is_periodic(:))) then
       call utils_assert(pub_is_pbe_neutralisation_scheme /= 'NONE', &
            'Solvation in PBCs requires is_pbe_neutralisation_scheme different &
            &from ''none'' or else electrostatic energy will diverge.')
    end if

    if(pub_is_pbe == 'NONE') then
       call utils_assert(pub_is_pbe_neutralisation_scheme == 'NONE' .or. &
            pub_is_pbe_neutralisation_scheme == 'JELLIUM', &
            'Non-Boltzmann solvation (is_pbe ''none'') does not need or support&
            &neutralisation schemes other than ''jellium'' (for PBCs) or &
            &''none'' (for OBCs).')
    end if

    if(.not. any(pub_multigrid_bc_is_periodic(:))) then
       call utils_assert(pub_is_pbe_neutralisation_scheme == 'NONE', &
            'Boltzmann solvation in OBCs does not need or support neutralisatio&
            &n schemes. Please set is_pbe_neutralisation_scheme to ''none''.')
    end if

    ! jd: Advise against using is_bc_surface_coarseness for charged systems
    if (pub_is_bc_surface_coarseness > 1 .and. pub_charge /= 0.0_DP &
         .and. pub_on_root) then
       write(stdout,'(//a)') 'WARNING: is_bc_surface_coarseness is larger than&
            & 1. Values larger than 1'
       write(stdout,'(a)') '         should only be used for charge-neutral&
            & systems -- for charged systems'
       write(stdout,'(a)') '         the accuracy of the energy&
            & calculation will likely be negatively'
       write(stdout,'(a)') '         impacted. Seriously&
            & consider setting is_bc_surface_coarseness to 1'
       write(stdout,'(a)') '         and if you find&
            & too much time is spent in multigrid_prepare_bound_cond,'
       write(stdout,'(a)') '         it''s&
            & better to (further) increase is_bc_coarseness rather than'
       write(stdout,'(a/)') '         is_bc_surface_coarseness.'
    end if

    ! jd: Check pub_timings_order is fine, support some aliases.
    if(trim(pub_timings_order)==trim('WALLTIME')) pub_timings_order = 'TIME'
    if(trim(pub_timings_order)==trim('CPUTIME')) pub_timings_order = 'TIME'
    if(trim(pub_timings_order)==trim('ALPHABETIC')) pub_timings_order = 'NAME'
    call utils_assert(trim(pub_timings_order) == trim('TIME') .or. &
         trim(pub_timings_order) == trim('NCALLS') .or. &
         trim(pub_timings_order) == trim('NONE') .or. &
         trim(pub_timings_order) == 'NAME', &
         'timings_order must be either ''TIME'' or ''NCALLS'' or ''NAME''&
         & or ''NONE''.')

    ! jd: Disallow EDFT in combination with conduction. It doesn't work yet.
    ! cf. email exchange of 2022.05.09 between jd and ndmh , and also
    ! https://bitbucket.org/onetep/onetep/issues/1893/edft-conduction-is-fundamentally-broken
    if(pub_edft .and. (pub_task == 'COND' .or. pub_task == 'PROPERTIES_COND')) then
       call utils_abort('Ensemble DFT ("edft T") cannot be used with conduction&
            & ("task COND" and/or "task PROPERTIES_COND").')
    end if

    ! ndmh: if pub_ngwf_cg_max_step is left at its default, rescale for varying ecut
    ! ndmh: 2.0 is about right at 600eV - higher cutoff requires higher max step
    if (pub_ngwf_cg_max_step < 0.0_DP) then
       pub_ngwf_cg_max_step = -pub_ngwf_cg_max_step * (pub_cutoff_energy / 22.04959837_DP)
    end if

    ! kaw: Sanity check for relationship between threading and batch size keywords
    if ((modulo(pub_fftbox_batch_size,pub_threads_num_fftboxes).ne.0).and.pub_on_root) then
       write(stdout,'(/a)') 'WARNING: The number of FFT boxes in a batch is &
            &not a multiple of the number of'
       write(stdout,'(a)') 'OpenMP threads. Consider changing fftbox_batch_size.'
       if (pub_threads_fftbox) then
          write(stdout,'(a)') 'However, you are also using threaded FFTs so this may be &
               &deliberate.'
       end if
    end if

    ! ebl: disabling DMFT capabilities while code is under development
    call utils_assert(.not.pub_dmft_fully_sc_h, 'This DMFT functionality is under development')
    call utils_assert(.not.pub_dmft_fully_sc, 'This DMFT functionality is under development')
    call utils_assert(pub_dmft_nkpoints==1, 'This DMFT functionality is under development')

    ! vv:
    md_iter_global = 0 ! This gets modified in md_mod
    call utils_assert(.not. (pub_mw_total_force .and. pub_zero_total_force), &
         'The keywords "zero_total_force" and "mw_total_force" are &
         &mutually exclusive. They cannot both be set to T.')

    ! vv:
    call utils_assert(pub_rand_seed == -1 .or.pub_rand_seed >= 0, &
         &'rand_seed must be either -1 or a positive integer.')

    ! jme
    call utils_assert(pub_paw .or. (.not. pub_perturbative_soc), 'ERROR in &
         &rundat_check_inputs: Perturbative Spin-Orbit Couplings require PAW.')

    ! jme
    call utils_assert( .not. (pub_write_nbo .and. pub_cmplx_ngwfs), &
         'WRITE_NBO=T is not yet compatible with complex NGWFs.')
    call utils_assert( .not. (pub_plot_nbo .and. pub_cmplx_ngwfs), &
         'PLOT_NBO=T is not yet compatible with complex NGWFs.')

    ! rjc: Check input for Electron Localisation Descriptor and Kinetic Energy
    ! rjc: Density output
    call utils_assert( .not. (pub_paw .and. pub_eld_calculate), &
         'Electron Localisation Descriptors not currently implemented with PAW. &
         & Calculation will abort.')
    if (pub_eld_calculate) then
         call utils_assert(pub_eld_function=='ELF' .or. pub_eld_function=='LOL', &
         & 'Unknown ELD_FUNCTION: please use "ELF" or "LOL". Calculation will &
         &abort.')
    end if
    call utils_assert( .not. (pub_paw .and. pub_ke_density_calculate), &
         'Kinetic Energy Density not currently implemented with PAW. &
         & Calculation will abort.')

  contains

     ! ************************************************************************
     ! ************************************************************************

     subroutine internal_check_bc
       !====================================================================!
       ! Ensures that boundary conditions are consistent between the diffe- !
       ! rent components. Prints out a warning, if not. Outputs a summary   !
       ! of the BCs in VERBOSE or higher.                                   !
       !--------------------------------------------------------------------!
       ! Extracted from rundate_check_inputs by Jacek Dziedzic in July 2022.!
       !====================================================================!

       use constants, only: VERBOSE, CRLF

       implicit none

       ! JCW: If pub_multigrid_bc is not empty, then check its validity
       if (pub_multigrid_hartree.and.len_trim(pub_multigrid_bc) > 0) then
          ! User has explicitly set BCs for multigrid solver
          ! - check each character is an allowed character
          ! - check that the number of non-whitespace characters is 3
          if (.not.internal_validate_string(pub_multigrid_bc,"OoPpZz",3,.true.)) then
             call utils_abort("multigrid_bc is not a valid BC specification: "//&
                  trim(pub_multigrid_bc))
          end if
          ! - check pub_openbc_hartree is .false. (the default) - if not, then
          !   check pub_multigrid_bc is consistent with this (all "O"s)
          if (pub_openbc_hartree.and.verify(pub_multigrid_bc," Oo")/=0) then
               call utils_abort("multigrid_bc == "//trim(pub_multigrid_bc)//" is not &
                    &compatible with openbc_hartree = .true.")
          end if
          ! - populate pub_multigrid_bc_is_{periodic,zero}
          call utils_parse_bc(pub_multigrid_bc,&
               bc_is_periodic = pub_multigrid_bc_is_periodic, &
               bc_is_zero     = pub_multigrid_bc_is_zero)
       end if

       ! - abort if cutoff Coulomb is in use, since use of multigrid with cutoff
       !   Coulomb is unsupported.
       if (pub_multigrid_hartree.and.pub_coulomb_cutoff) then
          call utils_abort("cannot use multigrid solver when cutoff Coulomb is in use")
       end if

       ! JCW: To allow flexible periodic and mixed BCs, but not break existing
       ! JCW: behaviour, pub_ion_ion_bc and pub_pspot_bc allow specification of
       ! JCW: periodic/open BCs as a character string.
       ! JCW: * If explicitly set by the user, these override pub_ii_energy_direct
       ! JCW:   and pub_open_localpseudo
       ! JCW: * If not explicitly set by user, these are set based on the values
       ! JCW:   of pub_ii_energy_direct and pub_open_localpseudo

       ! JCW: If pub_ion_ion_bc is not empty, then check its validity
       if (len_trim(pub_ion_ion_bc) > 0) then
          ! - check each character is an allowed character
          ! - check that the number of non-whitespace characters is 3
          if (.not.internal_validate_string(pub_ion_ion_bc,"OoPp",3,.true.)) then
             call utils_abort("ion_ion_bc is not a valid BC specification: "//&
                  trim(pub_ion_ion_bc))
          end if
          ! - populate pub_ion_ion_bc_is_periodic
          call utils_parse_bc(pub_ion_ion_bc,&
               bc_is_periodic = pub_ion_ion_bc_is_periodic)
       end if

       ! JCW: If pub_ion_ion_bc is empty but pub_ii_energy_direct is .true.,
       ! JCW: then we determine the boundary conditions based on the following
       ! JCW: logic:
       ! JCW: * If pub_coulomb_cutoff is .false., then we must be doing a solvation
       ! JCW:   calculation or Martyna-Tuckerman, which defaults to fully open BCs
       ! JCW: * If pub_coulomb_cutoff is .true., then we determine the BCs based
       ! JCW:   on the value of pub_coulomb_cutoff_type
       ! JCW: * Otherwise, use the usual fully periodic BCs
       if (len_trim(pub_ion_ion_bc) <= 0) then
          if (pub_ii_energy_direct) then
             if (.not.pub_coulomb_cutoff) then
                pub_ion_ion_bc_is_periodic(1:3) = .false.
             else
                select case (pub_coulomb_cutoff_type)

                case ("SPHERE","sphere","CYLINDER","cylinder")
                     ! Fully open BC
                     pub_ion_ion_bc_is_periodic(1:3) = .false.
                case ("WIRE","wire")
                     ! 1-D periodicity (z-direction hard-coded)
                     pub_ion_ion_bc_is_periodic = [ .false., .false., .true. ]
                case ("SLAB","slab")
                     ! 2-D periodicity (xy-plane hard-coded)
                     pub_ion_ion_bc_is_periodic = [ .true., .true., .false. ]
                case default
                     call utils_abort("unexpected value of coulomb_cutoff_type: "//&
                          pub_coulomb_cutoff_type)
                end select
             end if
          else
                pub_ion_ion_bc_is_periodic(1:3) = .true.
          end if
       end if


       ! JCW: If pub_pspot_bc is not empty, then check its validity
       if (len_trim(pub_pspot_bc) > 0) then
          ! - check each character is an allowed character
          ! - check that the number of non-whitespace characters is 3
          if (.not.internal_validate_string(pub_pspot_bc,"OoPp",3,.true.)) then
             call utils_abort("pspot_bc is not a valid BC specification: "//&
                  trim(pub_pspot_bc))
          end if
          ! - populate pub_pspot_bc_is_periodic
          call utils_parse_bc(pub_pspot_bc,&
               bc_is_periodic = pub_pspot_bc_is_periodic)
       end if

       ! JCW: If pub_pspot_bc is empty, then we determine the boundary conditions
       ! JCW: based on the following logic:
       ! JCW: * If pub_open_localpseudo, then we are doing a fully open BC
       ! JCW:   calculation
       ! JCW: * If pub_coulomb_cutoff is .true., then we the BCs are
       ! JCW:   periodic, but in the padded cell --> use the same BC designators
       ! JCW:   as pub_ion_ion_bc as the padding approximates open BCs
       ! JCW: * If pub_coulomb_cutoff is .false., and pub_open_localpseudo is .false.
       ! JCW:   then calculation uses the usual fully periodic BCs.
       if (len_trim(pub_pspot_bc) <= 0) then
          if (pub_open_localpseudo) then
                pub_pspot_bc_is_periodic(1:3) = .false.
          else if (pub_coulomb_cutoff) then
             select case (pub_coulomb_cutoff_type)

             case ("SPHERE","sphere","CYLINDER","cylinder")
                  ! Fully open BC
                  pub_pspot_bc_is_periodic(1:3) = .false.
             case ("WIRE","wire")
                  ! 1-D periodicity (z-direction hard-coded)
                  pub_pspot_bc_is_periodic = [.false.,.false.,.true.]
             case ("SLAB","slab")
                  ! 2-D periodicity (xy-plane hard-coded)
                  pub_pspot_bc_is_periodic = [.true.,.true.,.false.]
             case default
                  call utils_abort("unexpected value of coulomb_cutoff_type: "//&
                       pub_coulomb_cutoff_type)
             end select
          else
             pub_pspot_bc_is_periodic(1:3) = .true.
          end if
       end if

       ! JCW: If pub_smeared_ion_bc is not empty, then check its validity
       if (len_trim(pub_smeared_ion_bc) > 0) then
          ! - check each character is an allowed character
          ! - check that the number of non-whitespace characters is 3
          if (.not. internal_validate_string(pub_smeared_ion_bc, &
               "OoPp",3,.true.)) then
             call utils_abort(&
                  "smeared_ion_bc is not a valid BC specification: "//&
                  trim(pub_smeared_ion_bc))
          end if
          ! - populate pub_smeared_ion_bc_is_periodic
          call utils_parse_bc(pub_smeared_ion_bc,&
               bc_is_periodic = pub_smeared_ion_bc_is_periodic)
       end if

       ! JCW: If pub_smeared_ion_bc is empty then we determine the boundary
       ! JCW: conditions based on the following logic:
       ! JCW: * Since smeared ion representation is intended to be used with
       ! JCW:   multigrid solver, make pub_smeared_ion_bc_is_periodic equal to
       ! JCW:   pub_multigrid_bc_is_periodic (if pub_multigrid_hartree is .true.)
       ! JCW: * Otherwise, leave undefined (smeared ions are not used without
       ! JCW:   multigrid solver).
       if (len_trim(pub_smeared_ion_bc) <= 0) then
          if (pub_multigrid_hartree) then
             pub_smeared_ion_bc_is_periodic(1:3) = &
                  pub_multigrid_bc_is_periodic(1:3)
          end if
       end if

       ! JCW: Now, check for actually supported BCs for multigrid, ion-ion and
       ! JCW: local pseudopotential
       if (pub_multigrid_hartree) then
          ! JCW: multigrid BCs and smeared ion BCs are only defined when using the multigrid solver
          ! JCW: i.e. pub_multigrid_hartree is .true.
          if (.not. all(pub_multigrid_bc_is_periodic(1:3)) .and. &
               any(pub_multigrid_bc_is_periodic(1:3))) then
             call utils_abort("Mixed BC operation of the multigrid solver &
                  &is not currently supported")
          end if
          if (.not. all(pub_smeared_ion_bc_is_periodic(1:3)) .and. &
               any(pub_smeared_ion_bc_is_periodic(1:3))) then
             call utils_abort("Mixed BC smeared ion representation is not &
                  &currently supported")
          end if
       end if
       if (.not. all(pub_ion_ion_bc_is_periodic(1:3)) .and. &
            any(pub_ion_ion_bc_is_periodic(1:3)) .and. &
            .not.pub_coulomb_cutoff) then
          call utils_abort("Mixed BC ion-ion interaction is not &
               &currently supported when cutoff Coulomb is not in use")
       end if
       if (.not. all(pub_pspot_bc_is_periodic(1:3)) .and. &
            any(pub_pspot_bc_is_periodic(1:3)) .and. &
           .not.pub_coulomb_cutoff) then
          call utils_abort("Mixed BC local pseudopotential is not &
               &currently supported when cutoff Coulomb is not in use")
       end if

       ! ab: If pub_vdw_bc is not empty, then check its validity
       if (len_trim(pub_vdw_bc) > 0) then
          ! - check each character is an allowed character
          ! - check that the number of non-whitespace characters is 3
          if (.not.internal_validate_string(pub_vdw_bc,"OoPp",3,.true.)) then
             call utils_abort("vdw_bc is not a valid BC specification: "//&
                  trim(pub_vdw_bc))
          end if
          ! - populate pub_vdw_bc_is_periodic
          call utils_parse_bc(pub_vdw_bc,&
               bc_is_periodic = pub_vdw_bc_is_periodic)
       end if

       ! ab: If pub_vdw_bc is empty, set it same as ion_ion_bc
       if (len_trim(pub_vdw_bc) <= 0) then
          pub_vdw_bc_is_periodic(1:3) = pub_ion_ion_bc_is_periodic(1:3)
       end if

       ! ab: If pub_dftb_bc is not empty, then check its validity
       if (len_trim(pub_dftb_bc) > 0) then
          ! - check each character is an allowed character
          ! - check that the number of non-whitespace characters is 3
          if (.not.internal_validate_string(pub_dftb_bc,"OoPp",3,.true.)) then
             call utils_abort("dftb_bc is not a valid BC specification: "//&
                  trim(pub_dftb_bc))
          end if
          ! - populate pub_dftb_bc_is_periodic
          call utils_parse_bc(pub_dftb_bc,&
               bc_is_periodic = pub_dftb_bc_is_periodic)
       end if

       ! ab: If pub_dftb_bc is empty, set it same as ion_ion_bc
       if (len_trim(pub_dftb_bc) <= 0) then
          pub_dftb_bc_is_periodic(1:3) = pub_ion_ion_bc_is_periodic(1:3)
       end if

       ! jd: If pub_hfx_bc is not empty, then check its validity
       if (len_trim(pub_hfx_bc) > 0) then
          ! - check each character is an allowed character
          ! - check that the number of non-whitespace characters is 3
          if (.not.internal_validate_string(pub_hfx_bc,"OoPp",3,.true.)) then
             call utils_abort("hfx_bc is not a valid BC specification: "//&
                  trim(pub_hfx_bc))
          end if
          ! - populate pub_hfx_bc_is_periodic
          call utils_parse_bc(pub_hfx_bc,&
               bc_is_periodic = pub_hfx_bc_is_periodic)
       end if

       ! jd: If pub_hfx_bc is empty, set it same as ion_ion_bc
       if (len_trim(pub_hfx_bc) <= 0) then
          pub_hfx_bc_is_periodic(1:3) = pub_ion_ion_bc_is_periodic(1:3)
       end if

       ! jd: Warn the user that either PBC or OBC should be used, but not a mix
       ! ndmh: Warning was catching standard CC runs, so needed modified criteria
       ! JCW: Update 08/05/17
       ! JCW: With addition of mixed/periodic BCs for multigrid solver,
       ! JCW: pub_multigrid_hartree is no longer a good indicator of use of fully
       ! JCW: open BCs for the multigrid solver.
       ! JCW: => pub_multigrid_hartree is .true. regardless of BCs used by multigrid
       ! JCW:    solver
       ! JCW: pub_openbc_hartree is not a reliable indicator either, since this is
       ! JCW: only set when the user explicitly sets its value
       ! JCW: => pub_multigrid_bc is unambiguous, since it is always set to the
       ! JCW:    the BCs used by the multigrid solver, either by explicit user input
       ! JCW:    or by default, based on pub_openbc_hartree and pub_multigrid_hartree
       if (pub_multigrid_hartree) then
          mg_full_obc_hartree = .not. any(pub_multigrid_bc_is_periodic(1:3))
          mg_full_pbc_hartree = all(pub_multigrid_bc_is_periodic(1:3))
          mg_mixed_bc_hartree = .not. (mg_full_obc_hartree) .and. &
               (.not. mg_full_pbc_hartree)
       else
          mg_full_obc_hartree = .false.
          mg_full_pbc_hartree = .false.
          mg_mixed_bc_hartree = .false.
       end if

       ! JCW: Check for inconsistent BCs in ion-ion, local pseudo and multigrid solver
       ! JCW: (if in use)...
       inconsistent_bcs = .false.

       ! JCW: Check for inconsistent BCs in ion-ion and local pseudo
       ! JCW: Unlike pub_multigrid_bc, these are always set (pub_multigrid_bc is
       ! JCW: empty when the multigrid solver is not in use)
       ! ab: also check consistency with van der Waals BCs, DFTB BCs
       ! jd: ... and HFx BCs
       if (.not.all(pub_ion_ion_bc_is_periodic(1:3) .eqv. pub_pspot_bc_is_periodic(1:3)) .or. &
            .not.all(pub_ion_ion_bc_is_periodic(1:3) .eqv. pub_vdw_bc_is_periodic(1:3)) .or. &
            .not.all(pub_ion_ion_bc_is_periodic(1:3) .eqv. pub_dftb_bc_is_periodic(1:3)) .or. &
            .not.all(pub_ion_ion_bc_is_periodic(1:3) .eqv. pub_hfx_bc_is_periodic(1:3))) then
          inconsistent_bcs = .true.
       end if

       if (pub_multigrid_hartree) then
          ! JCW: Compare periodicity of multigrid to ion-ion, local pseudo and smeared ion terms
          ! jd: ... and vdW, DFTB, HFx
          if (.not. all(pub_multigrid_bc_is_periodic(1:3) .eqv. &
               pub_ion_ion_bc_is_periodic(1:3)) ) then
             inconsistent_bcs = .true.
          end if
          if (.not. all(pub_multigrid_bc_is_periodic(1:3) .eqv. &
               pub_pspot_bc_is_periodic(1:3)) ) then
             inconsistent_bcs = .true.
          end if
          if (.not. all(pub_multigrid_bc_is_periodic(1:3) .eqv. &
               pub_smeared_ion_bc_is_periodic(1:3)) ) then
             inconsistent_bcs = .true.
          end if
          if (.not. all(pub_multigrid_bc_is_periodic(1:3) .eqv. &
               pub_vdw_bc_is_periodic(1:3)) ) then
             inconsistent_bcs = .true.
          end if
          if (.not. all(pub_multigrid_bc_is_periodic(1:3) .eqv. &
               pub_dftb_bc_is_periodic(1:3)) ) then
             inconsistent_bcs = .true.
          end if
          if (.not. all(pub_multigrid_bc_is_periodic(1:3) .eqv. &
               pub_hfx_bc_is_periodic(1:3)) ) then
             inconsistent_bcs = .true.
          end if
       end if

       if(pub_debug_on_root .or. &
            (pub_on_root .and. &
            (pub_output_detail >= VERBOSE .or. inconsistent_bcs))) then
          call internal_print_bc
       end if

       if (inconsistent_bcs .and. pub_on_root) then
          pub_mix_bcs = .true.
          write(stdout,'(//a/)') 'WARNING: You are attempting a calculation &
               &with inconsistent boundary conditions'//CRLF//&
               '         (different energy components are calculated with &
               &different BCs).'//CRLF//'         Make sure you know what &
               &you are doing!'
       end if
       call comms_bcast(pub_root_proc_id, pub_mix_bcs)

       ! JCW: Add a further warning in the case that mixed BCs are selected for
       ! JCW: multigrid solver until the complementary mixed BC terms for
       ! JCW: ion-ion and local pseudopotential are implemented
       if (mg_mixed_bc_hartree .and. pub_on_root) then
          write(stdout,'(/a/)') 'WARNING: You have requested mixed boundary &
            &conditions for the multigrid solver. The corresponding mixed BC &
            &terms for other electrostatic interactions are not yet implemented &
            &so using the multigrid solver with mixed BCs will always give &
            &a total energy expression with inconsistent BCs.'
       end if

     end subroutine internal_check_bc

     ! ************************************************************************
     ! ************************************************************************

     ! ************************************************************************
     ! ************************************************************************
     subroutine internal_print_bc
       !====================================================================!
       ! Prints a summary of the boundary conditions in effect.             !
       ! To be called on the root proc only.                                !
       !--------------------------------------------------------------------!
       ! Written by Jacek Dziedzic in July 2022.                            !
       !====================================================================!

       implicit none

       ! ----------------------------------------------------------------------

       write(stdout,'(/a)') 'BC: +++++++++++++++ Summary of boundary conditions ++++++++++++++'
       write(stdout,'(a)') 'BC: |       Component       |     Keyword     |   Periodicity   |'
       write(stdout,'(a)') 'BC: |                       |                 |     X  Y  Z [1] |'
       write(stdout,'(a)') 'BC: |-----------------------------------------------------------|'
       write(stdout,'(a,a,a,a,a)') &
            'BC: | Ion-ion               | ion_ion_bc [2]  |     ', &
            merge('P  ','O  ',pub_ion_ion_bc_is_periodic(1)), &
            merge('P  ','O  ',pub_ion_ion_bc_is_periodic(2)), &
            merge('P  ','O  ',pub_ion_ion_bc_is_periodic(3)), '   |'
       write(stdout,'(a,a,a,a,a)') &
            'BC: | Local pseudopotential | pspot_bc [3]    |     ', &
            merge('P  ','O  ',pub_pspot_bc_is_periodic(1)), &
            merge('P  ','O  ',pub_pspot_bc_is_periodic(2)), &
            merge('P  ','O  ',pub_pspot_bc_is_periodic(3)), '   |'

      ! jd: We cannot rely on pub_hfx or pub_activehfx here because they will
      !     only be set in energy_and_force_init_cell().
      if(trim(pub_xc_functional) == 'HF' .or. &
            trim(pub_xc_functional) == 'B1LYP' .or. &
            trim(pub_xc_functional) == 'B1PW91' .or. &
            trim(pub_xc_functional) == 'PBE0' .or. &
            trim(pub_xc_functional) == 'B3LYP' .or. &
            trim(pub_xc_functional) == 'B3PW91' .or. &
            trim(pub_xc_functional) == 'X3LYP' .or. &
            trim(pub_active_xc_functional) == 'HF' .or. &
            trim(pub_active_xc_functional) == 'B1LYP' .or. &
            trim(pub_active_xc_functional) == 'B1PW91' .or. &
            trim(pub_active_xc_functional) == 'PBE0' .or. &
            trim(pub_active_xc_functional) == 'B3LYP' .or. &
            trim(pub_active_xc_functional) == 'B3PW91' .or. &
            trim(pub_active_xc_functional) == 'X3LYP') then
          write(stdout,'(a,a,a,a,a)') &
               'BC: | Hartree-Fock exchange | hfx_bc          |     ', &
               merge('P  ','O  ',pub_hfx_bc_is_periodic(1)), &
               merge('P  ','O  ',pub_hfx_bc_is_periodic(2)), &
               merge('P  ','O  ',pub_hfx_bc_is_periodic(3)), '   |'
       end if
       ! jd: DMA uses SWx, and so follows HFx. Only print this if it is in use
       if(pub_dma_max_l > -1) then
          write(stdout,'(a,a,a,a,a)') &
               'BC: | DMA                   | hfx_bc          |     ', &
               merge('P  ','O  ',pub_hfx_bc_is_periodic(1)), &
               merge('P  ','O  ',pub_hfx_bc_is_periodic(2)), &
               merge('P  ','O  ',pub_hfx_bc_is_periodic(3)), '   |'
       end if
       if(pub_multigrid_hartree) then
          write(stdout,'(a,a,a,a,a)') &
               'BC: | Multigrid solver      | multigrid_bc    |     ', &
               merge('P  ','O  ',pub_multigrid_bc_is_periodic(1)), &
               merge('P  ','O  ',pub_multigrid_bc_is_periodic(2)), &
               merge('P  ','O  ',pub_multigrid_bc_is_periodic(3)), '   |'
       end if
       if(pub_is_smeared_ion_rep) then
          write(stdout,'(a,a,a,a,a)') &
               'BC: | Smeared ions          | smeared_ion_bc  |     ', &
               merge('P  ','O  ',pub_multigrid_bc_is_periodic(1)), &
               merge('P  ','O  ',pub_multigrid_bc_is_periodic(2)), &
               merge('P  ','O  ',pub_multigrid_bc_is_periodic(3)), '   |'
       end if
       write(stdout,'(a,a,a,a,a)') &
            'BC: | Van der Waals         | vdw_bc          |     ', &
            merge('P  ','O  ',pub_vdw_bc_is_periodic(1)), &
            merge('P  ','O  ',pub_vdw_bc_is_periodic(2)), &
            merge('P  ','O  ',pub_vdw_bc_is_periodic(3)), '   |'
       if(pub_dftb) then
          write(stdout,'(a,a,a,a,a)') &
               'BC: | DFTB                  | dftb_bc         |     ', &
               merge('P  ','O  ',pub_dftb_bc_is_periodic(1)), &
               merge('P  ','O  ',pub_dftb_bc_is_periodic(2)), &
               merge('P  ','O  ',pub_dftb_bc_is_periodic(3)), '   |'
       end if
       write(stdout,'(a)') 'BC: +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(stdout,'(a)') 'BC: P: Periodic, O: Open.'
       write(stdout,'(a)') 'BC: [1] Not strictly along X, Y, Z but rather along lattice vectors A1, A2, A3.'
       write(stdout,'(a)') 'BC: [2] Also affected by: pbc_correction_cutoff, cutoff_coulomb_type,'
       write(stdout,'(a)') 'BC:     is_smeared_ion_rep, is_implicit_solvent, openbc_ion_ion.'
       write(stdout,'(a)') 'BC: [3] Also affected by: is_smeared_ion_rep, is_implicit_solvent, openbc_pspot.'

     end subroutine internal_print_bc

     subroutine internal_check_mix(mix_type, mix_num, mix_init_type, &
                         mix_init_num, mix_reset, mix_num_aux)

       !====================================================================!
       ! This subroutine imposes consistency between mixing parameters      !
       !--------------------------------------------------------------------!
       ! Written by Simon M.-M. Dubois on 25/03/2013.                       !
       !====================================================================!

       use utils, only: utils_use_var

       implicit none

       ! Arguments
       character(len=20), intent(inout) :: mix_type
       character(len=20), intent(inout) :: mix_init_type
       integer, intent(inout) :: mix_num
       integer, intent(inout) :: mix_init_num
       integer, intent(inout) :: mix_reset
       integer, optional, intent(in) :: mix_num_aux

       call utils_use_var(mix_init_type) ! jd: Kill 'arg unused' warning, consider removing the arg

       ! JA: avoid checking optional argument if not passed.
       if(present(mix_num_aux)) then
          if ((mix_type == 'PROJ' &
               .OR. mix_type == 'TENSOR')&
               .AND. mix_num_aux == 0) then
             mix_type = 'REUSE'
          endif
       else
          if (mix_type == 'PROJ' &
               .OR. mix_type == 'TENSOR') then
             mix_type = 'REUSE'
          endif
       end if

       if (mix_type == 'NONE') then
          if (mix_init_num > 0) then
             mix_num = 1
             mix_init_num = max(mix_init_num,0)
          else
             mix_num = 0
             mix_init_num = 0
          end if
       else if (mix_type == 'REUSE'&
           .OR. mix_type == 'PROJ' &
           .OR. mix_type == 'TENSOR') then
          mix_num = 1
          mix_init_num = max(mix_init_num,0)
       else if (mix_type == 'LINEAR') then
          mix_num = 2
          mix_init_num = max(mix_init_num,1)
       else if (mix_type == 'XLI' .or. &
          mix_type == 'NAIVE' .or. &
          mix_type == 'XLIS' ) then
          mix_num = 2
          mix_init_num = max(mix_num,mix_init_num)
       else if (mix_type == 'MULTID' &
           .OR. mix_type == 'POLY' &
           .OR. mix_type == 'LOCAL') then
          mix_num = max(2,mix_num)
          mix_init_num = max(mix_init_num,mix_num-1)
       else if (mix_type == 'TRPROP') then
          mix_num = 2
          mix_init_num = max(mix_init_num,2)
       else if (mix_type == 'XLD') then
          mix_num = min(max(mix_num,4),8)
          mix_init_num = max(mix_init_num,mix_num)
       else
          mix_type = 'NONE'
          mix_num = 0
          mix_init_num = 0
       endif
       mix_reset = max(mix_reset,mix_init_num+1)

     end subroutine internal_check_mix

     ! ************************************************************************
     ! ************************************************************************

     function internal_validate_string(string,allowed_chars,allowed_nchars,&
                ignore_whitespace) result(is_valid)
       !====================================================================!
       ! Check that string only contains characters in allowed_chars and is !
       ! that the number of allowed chararacters is equal to allowed_nchars.!
       !                                                                    !
       ! If the string fails either of these tests, is_valid is returned    !
       ! .false., otherwise .true.                                          !
       !                                                                    !
       ! If ignore_whitespace = .true., then whitespace characters are not  !
       ! checked or counted against allowed_nchars. In this case,           !
       ! allowed_chars should not include a whitespace character.           !
       !--------------------------------------------------------------------!
       ! Written by James C. Womack on 08/05/2017.                          !
       !====================================================================!

       implicit none

       ! Result
       logical                      :: is_valid

       ! Arguments
       character(len=*), intent(in) :: string
       character(len=*), intent(in) :: allowed_chars
       integer         , intent(in) :: allowed_nchars
       logical         , intent(in) :: ignore_whitespace

       ! Local variables
       character(len=4)  :: charbuffer ! 4 integer digits should be enough for anyone
       character(len=1)  :: ws_char
       integer           :: i, nchar

       is_valid = .true.

       ! Write allowed_nchars to a string so it can be output in error messages
       write(charbuffer,'(i4)') allowed_nchars

       ! If ignoring whitespace, set ws_char to " ", otherwise ""
       ! --> this will be concatenated with allowed_chars for verifying the
       !     characters in the string
       if (ignore_whitespace) then
          ws_char = " "
       else
          ws_char = ""


       if (verify(string,trim(ws_char)//allowed_chars)/=0) then
          ! One of the characters is a non-allowed character
          is_valid = .false.
          return
       end if

       end if
       ! - check that the number of non-whitespace characters is 3
       nchar = 0
       do i = 1, len(string)
          if (ignore_whitespace.and.string(i:i) == " ") cycle
          if (scan(string(i:i),allowed_chars)/=0) nchar = nchar + 1
       end do
       if (nchar /= allowed_nchars) then
          is_valid = .false.
          return
       end if

       ! If you make it to here, then the string is valid
       return

     end function internal_validate_string

  end subroutine rundat_check_inputs

  subroutine rundat_exit()

    !=======================================================================!
    ! This subroutine deallocates those the allocatable module variables in !
    ! the rundat module (if required).                                      !
    !-----------------------------------------------------------------------!
    ! Written by Quintin Hill on 27/05/2009.                                !
    !=======================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    integer :: ierr

    if (allocated(pub_ldos_groups)) then
       deallocate(pub_ldos_groups,stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_ldos_groups',ierr)
    end if
    if (allocated(pub_ldos_group_nsp)) then
       deallocate(pub_ldos_group_nsp,stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_ldos_group_nsp',ierr)
    end if
    if (allocated(pub_bsunfld_groups)) then
       deallocate(pub_bsunfld_groups,stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_bsunfld_groups',ierr)
    end if
    if (allocated(pub_bsunfld_projatoms)) then
       deallocate(pub_bsunfld_projatoms,stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_bsunfld_projatoms',ierr)
    end if
    if (allocated(pub_bsunfld_group_nsp)) then
       deallocate(pub_bsunfld_group_nsp, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_bsunfld_group_nsp',ierr)
    end if
    if (allocated(pub_ldos_group_nsp)) then
       deallocate(pub_ldos_group_nsp,stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_ldos_group_nsp',ierr)
    end if
    if (allocated(pub_bs_kpoint_path_start)) then
       deallocate(pub_bs_kpoint_path_start, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_bs_kpoint_path_start',ierr)
    end if
    if (allocated(pub_bs_kpoint_path_end)) then
       deallocate(pub_bs_kpoint_path_end, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_bs_kpoint_path_end',ierr)
    end if
    ! agreco: deallocate k-point array if necessary
    if (allocated(pub_kpoint_list)) then
       deallocate(pub_kpoint_list, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_kpoint_list',ierr)
    end if

    ! ja531 -> local pdos
    if (allocated(pub_pdos_groups)) then
       deallocate(pub_pdos_groups,stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_pdos_groups',ierr)
    end if
    if (allocated(pub_pdos_group_nsp)) then
       deallocate(pub_pdos_group_nsp,stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_pdos_group_nsp',ierr)
    end if

    if (allocated(pub_scissor_groups)) then
       deallocate(pub_scissor_groups,stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_scissor_groups',ierr)
    end if
    if (allocated(pub_scissor_group_nsp)) then
       deallocate(pub_scissor_group_nsp,stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_scissor_group_nsp',ierr)
    end if
    if (allocated(pub_scissor_group_shifts)) then
       deallocate(pub_scissor_group_shifts,stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_scissor_group_shifts',ierr)
    end if

    ! lpl: DDEC stuff
    if(allocated(pub_ddec_rcomp)) then
       deallocate(pub_ddec_rcomp, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_ddec_rcomp',ierr)
    end if

    ! lpl: Related to the still-dodgy subroutine 'ddec_rmse'
    !if (allocated(pub_ddec_rmse_vdW)) then
    !   deallocate(pub_ddec_rmse_vdW, stat=ierr)
    !   call utils_dealloc_check('rundat_exit','pub_ddec_rmse_vdW',ierr)
    !end if
    ! lpl: DDEC stuff

    ! lpl: NBO stuff
    if (allocated(pub_nbo_write_species)) then
       deallocate(pub_nbo_write_species, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_nbo_write_species',ierr)
    end if
    if (allocated(pub_nbo_ngwf_label)) then
       deallocate(pub_nbo_ngwf_label, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_nbo_ngwf_label',ierr)
    end if
    if (allocated(pub_nbo_list_plotnbo)) then
       deallocate(pub_nbo_list_plotnbo, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_nbo_list_plotnbo',ierr)
    end if
    ! lpl: NBO stuff

    ! ndmh: TDDFT kernel groups
    if (allocated(pub_tddft_kernel_group_nsp)) then
       deallocate(pub_tddft_kernel_group_nsp,stat=ierr)
       call utils_dealloc_check('rundat_exit', &
            'pub_tddft_kernel_group_nsp',ierr)
    end if
    if (allocated(pub_tddft_kernel_groups)) then
       deallocate(pub_tddft_kernel_groups,stat=ierr)
       call utils_dealloc_check('rundat_exit', &
            'pub_tddft_kernel_groups',ierr)
    end if

    if (allocated(pub_tddft_ct_group_nsp)) then
       deallocate(pub_tddft_ct_group_nsp,stat=ierr)
       call utils_dealloc_check('rundat_exit', &
            'pub_tddft_ct_group_nsp',ierr)
    end if
    if (allocated(pub_tddft_ct_groups)) then
       deallocate(pub_tddft_ct_groups,stat=ierr)
       call utils_dealloc_check('rundat_exit', &
            'pub_tddft_ct_groups',ierr)
    end if
    ! ndmh: TDDFT kernel groups

    ! mjsp: EDA metadata
    if (allocated(pub_frag_iatm)) then
       deallocate(pub_frag_iatm, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_frag_iatm',ierr)
    end if
    if (allocated(pub_frag_atoms)) then
       deallocate(pub_frag_atoms, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_frag_atoms',ierr)
    end if
    if (allocated(pub_frag_charge)) then
       deallocate(pub_frag_charge, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_frag_charge',ierr)
    end if
    if (allocated(pub_frag_file_prefix)) then
       deallocate(pub_frag_file_prefix, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_frag_file_prefix',ierr)
    end if
    ! mjsp: EDA metadata

    ! rc2013: embedding stuff
    if (allocated(pub_ngwf_regions_groups)) then
       deallocate(pub_ngwf_regions_groups, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_ngwf_regions_groups',ierr)
    end if
    if (allocated(pub_ngwf_regions_group_nsp)) then
       deallocate(pub_ngwf_regions_group_nsp, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_ngwf_regions_group_nsp',ierr)
    end if
    ! rc2013: embedding stuff

  end subroutine rundat_exit

  ! --------------------------------------------------------------------------
  subroutine rundat_set_pub_use_hfx(value)
    !=======================================================================!
    ! Sets the protected variable pub_use_hfx. Needed in xc.                !
    !-----------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 2014/07/22.                                !
    !=======================================================================!

    use utils, only: utils_abort

    implicit none

    logical, intent(in) :: value

    pub_use_hfx = value

    if(value) then
       if(trim(pub_hfx_use_ri) == trim('<unset>')) then
          call utils_abort(&
               '''hfx_use_ri'' is mandatory when using Hartree-Fock exchange.')
       end if

       if(pub_hfx_max_l == -1) then
          call utils_abort(&
               '''hfx_max_l'' is mandatory when using Hartree-Fock exchange.')
       end if

       if(pub_hfx_max_q == -1) then
          call utils_abort(&
               '''hfx_max_q'' is mandatory when using Hartree-Fock exchange.')
       end if

       if(pub_paw) then
          call utils_abort('Hartree-Fock exchange cannot be used &
               &simultaneously with PAW, sorry.')
       end if
    else
       if(pub_hhf_nstates > 0) then
          call utils_abort('Hartree-Fock exchange is mandatory when using &
               &hyper Hartree-Fock (hhf_nstates)')
       end if
    end if

  end subroutine rundat_set_pub_use_hfx

  ! --------------------------------------------------------------------------
  subroutine rundat_set_pub_use_swx(value)
    !=======================================================================!
    ! Sets the protected variable pub_use_swx. Needed in energy_and_force.  !
    !-----------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 2014/07/22.                                !
    !=======================================================================!

    implicit none

    logical, intent(in) :: value

    pub_use_swx = value

  end subroutine rundat_set_pub_use_swx

  ! --------------------------------------------------------------------------
  subroutine rundat_set_pub_use_activehfx(value)
    !=======================================================================!
    ! Sets the protected variable pub_use_activehfx. Needed in xc.          !
    !-----------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 2014/07/22.                                !
    ! Copied and modified for active subregion by Joseph Prentice, 23/03/18 !
    !=======================================================================!

    use utils, only: utils_abort

    implicit none

    logical, intent(in) :: value

    ! rc2013: if we're not doing EMFT then active HFx should always be false
    if(pub_emft .or. pub_emft_follow) then
       pub_use_activehfx = value
    else
       pub_use_activehfx = .false.
    end if

    if(value) then
       if(trim(pub_hfx_use_ri) == trim('<unset>')) then
          call utils_abort(&
               '''hfx_use_ri'' is mandatory when using Hartree-Fock exchange.')
       end if

       if(pub_hfx_max_l == -1) then
          call utils_abort(&
               '''hfx_max_l'' is mandatory when using Hartree-Fock exchange.')
       end if

       if(pub_hfx_max_q == -1) then
          call utils_abort(&
               '''hfx_max_q'' is mandatory when using Hartree-Fock exchange.')
       end if

       if(pub_paw) then
          call utils_abort('Hartree-Fock exchange cannot be used &
               &simultaneously with PAW, sorry.')
       end if
    else
       if(pub_hhf_nstates > 0) then
          call utils_abort('Hartree-Fock exchange is mandatory when using &
               &hyper Hartree-Fock (hhf_nstates)')
       end if
    end if

  end subroutine rundat_set_pub_use_activehfx

  ! --------------------------------------------------------------------------
  subroutine rundat_set_pub_use_activeswx(value)
    !=======================================================================!
    ! Sets the protected variable pub_use_activeswx.                        !
    !-----------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 2014/07/22.                                !
    ! Copied and modified for active subregion by Joseph Prentice, 23/03/18 !
    !=======================================================================!

    implicit none

    logical, intent(in) :: value

    ! rc2013: if we're not doing EMFT then use_activeswx should always be false
    if(pub_emft .or. pub_emft_follow) then
       pub_use_activeswx = value
    else
       pub_use_activeswx = .false.
    end if


  end subroutine rundat_set_pub_use_activeswx

  ! --------------------------------------------------------------------------
  subroutine rundat_set_pub_use_sph_harm_rot_true
    !=======================================================================!
    ! Sets the protected variable pub_use_sph_harm_rot to .true.            !
    ! This allows functionality to activate the sph_harm_rotation module,   !
    ! during a calculation, e.g. sw_resolution_of_identity can set the      !
    ! variable to .true. when a metric matrix evaluation scheme requiring   !
    ! spherical harmonic rotation is selected.                              !
    !                                                                       !
    ! This is a ratchet. If the variable is already .true. (because some    !
    ! other functionality set it), then it cannot be returned to .false.    !
    !-----------------------------------------------------------------------!
    ! Written by James Womack, 26/07/2018 (based on rundat_set_pub_use_swx  !
    ! by J. Dziedzic)                                                       !
    !=======================================================================!

    implicit none

    pub_use_sph_harm_rot = .true.

  end subroutine rundat_set_pub_use_sph_harm_rot_true

end module rundat
