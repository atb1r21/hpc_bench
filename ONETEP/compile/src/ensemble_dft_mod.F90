! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
module ensemble_dft

  !==========================================================================!
  ! This module contains the subroutines to perform a direct energy          !
  ! minimisation of the occupancies based on the Ensemble-DFT method.        !
  !                                                                          !
  ! See Marzari, Vanderbilt and Payne, Phys. Rev. Lett. 79, 1337-1340 (1997).!
  !--------------------------------------------------------------------------!
  ! Written by Alvaro Ruiz Serrano in April 2012.                            !
  !==========================================================================!

  use constants, only: DP
  use ensemble_dft_type, only: EDFT_MODEL

  implicit none

  private

  ! ars: public subroutines
  public :: edft_create
  public :: edft_destroy
  public :: edft_calculate
  public :: edft_lagrangian
  public :: edft_matrix_from_eigs
  public :: edft_ngwf_gradient
  public :: edft_reorthogonalise_mo
  public :: edft_diag_ngwf_ham
  public :: edft_fermi_level
  public :: edft_zerok_fermi
  public :: edft_fermidirac_entropy
  public :: edft_pulay_history_flush
  public :: edft_rms_grad
  public :: edft_commutator
  public :: edft_check_convergence
  public :: edft_check_params
  public :: edft_print_iter
  public :: edft_print_line_search
  public :: edft_trial_step
  public :: edft_apply_smearing_operator
  public :: edft_fit_second_order_poly_3points
  public :: edft_fit_third_order_poly_4points
  public :: edft_fit_fourth_order_poly_5points
  public :: edft_diis_shift
  public :: edft_diis_pulay_mix
  public :: edft_diis_find_ientry
  public :: edft_print_qc


  ! ars: edft_matrix_from_eigs
  interface edft_matrix_from_eigs
     module procedure edft_spam3_matrix_from_eigs
     module procedure edft_dem_matrix_from_eigs
  end interface edft_matrix_from_eigs

  ! rab207: edft_fermi_level
  interface edft_fermi_level
     module procedure edft_fermi_level_gamma
     module procedure edft_fermi_level_kpoints
  end interface edft_fermi_level


  ! ja531-> FOE Chemical potential search history variables.
  real(kind=dp), parameter :: delta_mu_fac=10.0_dp ! Factor to multiply chemical potential search window by.
  ! This is set to be 10.0. I don't know if this is good or bad...
  ! Could do with some testing, if it turns out to be important then is could become an input parameter.
  integer, parameter :: mu_lookup_trace=1, mu_lookup_norm=2, mu_lookup_mu=3
  real(kind=DP), save, allocatable :: mu_lookup(:,:,:) ! A table to correspond chemical potentials with Hamiltonians
  integer, parameter :: pub_mu_lookup_maxlen=10 ! Number of mu_lookup entries to store (should be moved to rundat)
  integer, save, allocatable :: mu_lookup_len(:) ! Current number of entries in the mu_lookup array (num spin channels).
  integer, save, allocatable :: mu_lookup_idx(:) ! Current index in the mu_lookup array (num spin channels).

contains




  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  !~~~~~~~~~~~~~~~~~~~~~ ENSEMBLE DFT OCCUPANCIES OPTIMISATION ~~~~~~~~~~~~~~~~!
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!


  subroutine edft_calculate(edft, denskern, ham, rep, lhxc_fine, &
       ngwf_basis, hub_proj_basis, hub, mdl, hfxstate, current_edft_maxit, &
       restart_rootname, force_no_IO)

    !==========================================================================!
    ! Self-consistent optimisation of the occupancies and the unitary rotations!
    ! that transform from NGWF to Kohn-Sham representation, based on the       !
    ! Ensemble-DFT method.                                                     !
    ! @docargs                                                                 !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2012.                            !
    ! Modified for embedding by Joseph Prentice, May 2018.                     !
    ! Minor fix to avoid premature convergence by Jacek Dziedzic, January 2021.!
    !==========================================================================!

    use comms, only: pub_on_root, pub_total_num_procs, comms_reduce
    use constants, only: DP, stdout, paw_en_size, VERBOSE, NORMAL, BRIEF, &
         very_big_double
    use dense, only: DEM, dense_create, dense_destroy, dense_convert, dense_copy, &
         dense_axpy, dense_scale, dense_write
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_dens_dep_matrices, &
         hamiltonian_build_matrix
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: kernel_rescale_spam3
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use restart, only: restart_kernel_write, restart_kernel_read, &
         restart_hamiltonian_write, restart_hamiltonian_read
    use rundat, only: pub_debug_on_root, pub_edft_smearing_width, &
         pub_output_detail, pub_edft_max_step, &
         pub_write_converged_dk_ngwfs, pub_write_denskern, &
         pub_read_hamiltonian, pub_write_hamiltonian, &
         pub_edft_round_evals, pub_print_qc, pub_num_spins, pub_read_denskern, &
         pub_spin_fac, pub_foe, pub_num_kpoints, PUB_1K, &
         pub_xc_ke_density_required, pub_edft_rms_gradient_thres, &
         pub_contracoham_radmult, &
         pub_dense_foe, pub_inner_loop_iteration, pub_edft_spin_fix, &
         pub_real_spin, pub_edft_trial_step, pub_edft_update_scheme, &
         pub_edft_ham_diis_size, pub_devel_code, pub_edft_grand_canonical
    use smearing_operator, only: apply_fermi_function, calculate_fermi_entropy
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
         sparse_embed_create, sparse_embed_destroy, sparse_embed_trace,&
         sparse_embed_copy, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_product, sparse_embed_transpose
    use timer, only: timer_clock, wrappers_etime
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, utils_abort, &
                     utils_devel_code

    implicit none

    ! ars: arguments
    type(EDFT_MODEL), intent(inout) :: edft
    type(SPAM3_EMBED_ARRAY),intent(inout) :: denskern
    type(NGWF_HAM),   intent(inout) :: ham
    type(NGWF_REP),   intent(in)    :: rep
    type(MODEL),      intent(in   ) :: mdl
    type(HFX_STATE),  intent(inout), target :: hfxstate
    real(kind=DP),    intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub)
    type(FUNC_BASIS), intent(in   ) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in   ) :: hub_proj_basis(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    integer, intent(inout) :: current_edft_maxit
    character(len=*),optional,intent(in) :: restart_rootname
    logical, optional, intent(in) :: force_no_IO


    ! ars: iteration counters, spin, orbital and error flag
    integer :: edft_iter, is, iorb, ierr, icol, nrows, ientry, diis_size

    ! ars: isolated components of the energy
    real(kind=DP) :: lhxc_energy,hubbard_energy,paw_sphere_energies(paw_en_size)

    ! ars: energy
    real(kind=DP) :: e0, e1, e2, e3, e4, e5, e6, eaccepted

    ! ars: entropy
    real(kind=DP) :: s0, s1, s2, s3, s4, s5, s6, saccepted

    ! ars: lagrangian
    real(kind=DP) :: l0, l1, l2, l3, l4, l5, l6, laccepted

    ! ars: steps
    real(kind=DP) :: t0, t1, t2, t3, t4, t5, t6, taccepted

    ! ars: predicted lagrangian after line search
    real(kind=DP) :: predicted_l3, predicted_l4, predicted_l5

    ! ars: change in the accepted lagrangian between consecutive iterations
    real(kind=DP) :: deltal, deltae, deltas
    real(kind=DP), save :: prev_laccepted = 0.0_DP
    real(kind=DP), save :: prev_eaccepted = 0.0_DP
    real(kind=DP), save :: prev_saccepted = 0.0_DP
    real(kind=DP), save :: prev_fermi1 = 0.0_DP
    real(kind=DP), save :: prev_fermi2 = 0.0_DP
    real(kind=DP)       :: prev_integrated_ne(1:pub_num_spins)
    real(kind=DP)       :: delta_integrated_ne(1:pub_num_spins)

    ! kkbd temporary total occupancy
    real(kind=DP),dimension(1:pub_num_spins,pub_num_kpoints) :: t_occ

    ! ars: gradients at every step
    real(kind=DP) :: rms_gradient

    ! ars: [h,f] commutator
    real(kind=DP) :: commutator

    ! ars: Fermi level
    real(kind=DP) :: delta_fermi(1:pub_num_spins)

    ! ars: converged?
    logical :: edft_converged

    ! ars: flags for accepted steps
    logical :: quadratic_success, cubic_success, quartic_success
    logical :: quadratic_accepted,cubic_accepted,quartic_accepted,tiny_accepted

    ! ars: NGWF Hamiltonian
    type(SPAM3_EMBED), allocatable :: initial_ham(:), trial_ham(:)

    ! ars: variables used during the besy step
    logical :: L_increase(3)
    integer :: cc

    ! ars: rescale density kernel

    ! ars: golden section search
    real(kind=DP), parameter :: golden_section = (3.0_DP-sqrt(5.0_DP))*0.5_DP
    real(kind=DP) :: xx2(2), xx3(3), xx4(4)
    real(kind=DP) :: yy2(2), yy3(3), yy4(4)
    integer :: npts

    type(DEM) :: denswork(4)

    type(SPAM3_EMBED),allocatable :: tmp_spam(:)

    type(SPAM3_EMBED),allocatable :: hopping(:), testmat(:),testmat2(:)

    real(kind=DP) :: tmp_fermi(2)

    real(kind=DP) :: tmp_s
    integer       :: foe_entropy_approx=2   ! Whether to approximate the entropy (0: exact, 1: approx, 2:refined approx)
    logical       :: foe_line_search_entropy=.true. ! Whether to calculate entropy during the line-search

    real(kind=DP) :: boundl
    real(kind=DP) :: boundu, bounds(2)
    real(kind=DP) :: fermi_tol
    ! agrecocmplx
    logical :: loc_cmplx
    character(len=5) :: foe_mode
    logical :: loc_force_no_IO

    real(kind=dp) :: norm
    real(kind=dp) :: inv_time, start_time
    type(SPAM3_EMBED) :: tmpmat, tmpmat2

    ! gab: Shifted switch for Pulay.
    logical :: shifted

    ! gab: Parameters for safety
    integer :: n_safety
    logical :: safety_ls


    !-------------------------------------------------------------------------
    ! ars: initialise

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering edft_calculate'

    ! jme: KPOINTS_DANGER
    ! First k-point enforced in all uses of n_occ via PUB_1K
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine edft_calculate not yet compatible with more&
         & than one k-point.')
    ! JCW: Abort if KE density required (meta-GGA + ensemble DFT untested).
    call utils_assert(.not.pub_xc_ke_density_required, "Error in &
         &edft_calculate: pub_xc_ke_density_required is true, but &
         &combination of KE-density-dependent XC functionals with ensemble &
         &DFT has not been implemented/tested.")

    ! ars: start timer
    call timer_clock('edft_calculate',1)

    ! ab: for grand canonical ensemble, integrated number of electrons
    ! per spin channel will change so store them for later use
    if (pub_edft_grand_canonical) then
          prev_integrated_ne(1:pub_num_spins) = &
          edft%integrated_ne(1:pub_num_spins)
    end if

    boundl = -5.0_dp
    boundu = 5.0_dp
    bounds(1)=boundl
    bounds(2)=boundu
    fermi_tol = 1e-9_dp

    if(pub_dense_foe) then
       foe_mode="DENSE"
    else
       foe_mode="SPAM3"
    end if

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(1)%iscmplx

    ! ars: print header
    if(pub_on_root.and.pub_output_detail.ge.NORMAL) then
       write(stdout,'(/,a80)') "--------------------------- &
            &Ensemble-DFT optimisation --------------------------"


       if(pub_foe) then
          if(pub_dense_foe) then
#ifdef SCALAPACK
             write(stdout,'(a41, i4, a34)') 'Running the operator expansion method on ', &
                  pub_total_num_procs, ' processor(s), with parallel BLAS.'
#else
             write(stdout,'(a41, i4, a25)') 'Running the operator expansion method on ', &
                  pub_total_num_procs, ' processor(s), with BLAS.'
             write(stdout,'(a76)') 'Beware of serious memory limitations &
                  &due to replicated large dense matrices.'
             write(stdout,'(a)') 'Compile with -DSCALAPACK to use parallel BLAS for large systems.'
#endif
          else
             write(stdout,'(a41, i4, a34)') 'Running the operator expansion method on ', &
                  pub_total_num_procs, ' processor(s), using SPAM3.'
          end if
       else
#ifdef SCALAPACK
          write(stdout,'(a21, i4, a14)') 'Running ScaLAPACK on ', &
               pub_total_num_procs, ' processor(s).'
#else
          write(stdout,'(a18, i4, a14)') 'Running LAPACK on ', &
               pub_total_num_procs, ' processor(s).'
          write(stdout,'(a76)') 'Beware of serious memory limitations &
               &due to replicated large dense matrices.'
          write(stdout,'(a)') 'Use ScaLAPACK for large systems.'
#endif
       end if

    end if

    ! ars: check that edft has been initialised before entering this routine
    call utils_assert(edft%allocd,"Error in edft_calculate:&
         & please allocate 'edft' before calling this routine")

    ! ars: check input parameters
    call edft_check_params(current_edft_maxit)

    ! ars: allocate structures
    allocate(initial_ham(1:pub_num_spins), stat=ierr)
    call utils_alloc_check('edft_calculate', 'initial_ham', ierr)
    allocate(trial_ham(1:pub_num_spins), stat=ierr)
    call utils_alloc_check('edft_calculate', 'trial_ham', ierr)


    ! ars: create DEM structures
    do is = 1, pub_num_spins
       call sparse_embed_create(initial_ham(is), ham%ham(is),arrname='initial_ham')
       call sparse_embed_create(trial_ham(is), ham%ham(is),arrname='trial_ham')
    end do

    ! gab: Flush the hamiltonian history for Pulay mixing

    call edft_pulay_history_flush(edft,ham,cmplx=loc_cmplx)

    if(pub_foe) then
       ! agrecocmplx
       if(foe_mode=="DENSE") then
          call dense_create(denswork(1),edft%num,edft%num,iscmplx=loc_cmplx)
          call dense_create(denswork(2),edft%num,edft%num,iscmplx=loc_cmplx)
          call dense_create(denswork(3),edft%num,edft%num,iscmplx=loc_cmplx)
          call dense_create(denswork(4),edft%num,edft%num,iscmplx=loc_cmplx)
       end if
    end if

    ! ars: initialise variables
    e0=0.0_DP;e1=0.0_DP;e2=0.0_DP;e3=0.0_DP;e4=0.0_DP;e5=0.0_DP;e6=0.0_DP
    s0=0.0_DP;s1=0.0_DP;s2=0.0_DP;s3=0.0_DP;s4=0.0_DP;s5=0.0_DP;s6=0.0_DP
    l0=0.0_DP;e1=0.0_DP;l2=0.0_DP;l3=0.0_DP;l4=0.0_DP;l5=0.0_DP;l6=0.0_DP
    t0=0.0_DP;t1=0.0_DP;t2=0.0_DP;t3=0.0_DP;t4=0.0_DP;t5=0.0_DP;t6=0.0_DP
    eaccepted=0.0_DP;saccepted=0.0_DP;laccepted=0.0_DP
    taccepted=golden_section*pub_edft_max_step
    predicted_l3=0.0_DP;predicted_l4=0.0_DP;predicted_l5=0.0_DP
    deltal = very_big_double         ! jd: We don't initialise to zero, but
    deltae = very_big_double         !     rather to a big value to prevent
    deltas = very_big_double         !     premature convergence after step #0
    delta_fermi(:) = very_big_double !     when other criteria are satisfied, but
                                     !     a delta is not yet available
    rms_gradient=0.0_DP;commutator=0.0_DP
    edft_converged=.false.;L_increase(:)=.false.;cc=0
    quadratic_success=.false.;cubic_success=.false.;quartic_success=.false.
    quadratic_accepted=.false.;cubic_accepted=.false.
    quartic_accepted=.false.;tiny_accepted=.false.


    !=======================================================================!
    !                          OCCUPANCIES LINE SEARCH                      !
    !=======================================================================!
    ! ars: At every step we perform the following operations:               !
    !      1) Calculate non-FD H_\alpha\beta.                               !
    !      2) Diagonalise H_\alpha\beta. -> store e_i and M^\alpha_i        !
    !         (eigenvalues and rotation from MO to NGWF).                   !
    !      3) Build FD occupancies from e_i and calculate Fermi level.      !
    !      4) Build FD K^{\alpha\beta}.                                     !
    !      5) Build density and calculate energy.                           !
    !      6) Calculate FD entropy.                                         !
    !      7) Evaluate Lagrangian.                                          !
    !      8) Check for success in the step.                                !
    !         * If successful, start next iteration:                        !
    !          H_initial=H_trial (both non-FD).                             !
    !         * If not successful, move to next step or, if it was the last !
    !          step, choose the one that returns the lowest Lagrangian.     !
    !=======================================================================!
    ! ars: Note: the search direction in the Hamiltonian space is:          !
    !      D_\alpha\beta = H_\alpha\beta (FD) - H_\alpha\beta(initial)      !
    !      were H(FD) is the Hamiltonian constructed from a density that    !
    !      obeys the Fermi-Dirac distribution.                              !
    !=======================================================================!



    !============================ INITIALISATION ==============================!

    loc_force_no_IO = .false.
    if(present(force_no_IO)) then
       loc_force_no_IO=force_no_IO
    end if

    if(.not.loc_force_no_IO) then
       if (pub_read_hamiltonian.and.(.not.edft%initialised)) then
          ! ars: read initial Hamiltonian from a file
          if (present(restart_rootname)) then
             call restart_hamiltonian_read(ham%ham, &
                  restart_rootname=restart_rootname)
          else
             call restart_hamiltonian_read(ham%ham)
          end if

       else if (pub_read_denskern.and.(.not.edft%initialised)) then
          ! rab207: read denskern and build the hamiltonian
          call restart_kernel_read(denskern,read_cond=.false.)

          call hamiltonian_dens_dep_matrices(ham, lhxc_fine, eaccepted, &
               lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
               ngwf_basis, hub_proj_basis, hub, &
               denskern, mdl, hfxstate, &
               ham_update=.true., lhxc_fixed=.false.)
          call hamiltonian_build_matrix(ham, rep)

       end if
    end if


    ! ars: copy incoming Hamiltonian
    do is = 1, pub_num_spins
       call sparse_embed_copy(trial_ham(is), ham%ham(is))
    end do


    tmp_fermi(1:pub_num_spins)=edft%fermi(1:pub_num_spins)

    if (edft%initialised) then
       if(pub_foe) then
          allocate(tmp_spam(2),stat=ierr)
          call utils_alloc_check('edft_calculate','tmp_spam',ierr)
          call sparse_embed_create(tmp_spam(1),ham%ham(1), edft%rot,arrname='tmp_spam1')
          call sparse_embed_create(tmp_spam(2),edft%rot,arrname='tmp_spam2')
       end if

       ! ars: this is not the first NGWF iteration
       ! ars: assumes ngwf_cg_rotate=T
       ! ars: rescale the density kernel
       ! ars: rebuild the Hamiltonian with the rotated MOs

       ! kkbd: rescale by number of electons
       t_occ(:,PUB_1K) = edft%integrated_ne(:)

       if (pub_debug_on_root)write(stdout,*) &
               &"DEBUG: EDFT: Rescaling denskern to ",t_occ(:,PUB_1K)
       call kernel_rescale_spam3(denskern, rep%overlap, t_occ, &
            silent=.true.)
       do is = 1, pub_num_spins
          !     sparse_trace(denskern%m(is,PUB_1K),rep%overlap)
          if(.not.pub_foe) then
             call edft_matrix_from_eigs(trial_ham(is), &
                  edft%mo(is), edft%h_evals(:,is), edft%num, rep%overlap)

          else
             !            edft%rot^T * ham * edft%rot
             call sparse_embed_product(tmp_spam(1), edft%old_ham(is), edft%rot)
             call sparse_embed_transpose(tmp_spam(2), edft%rot)
             call sparse_embed_product(trial_ham(is), tmp_spam(2), tmp_spam(1))
          end if
       end do

       if(pub_foe) then
          call sparse_embed_destroy(tmp_spam(1),arrname='tmp_spam1')
          call sparse_embed_destroy(tmp_spam(2),arrname='tmp_spam2')
          deallocate(tmp_spam,stat=ierr)
          call utils_dealloc_check('edft_calculate','tmp_spam',ierr)
          saccepted=edft%entropy
       end if

    else ! not edft%initialised

       if(pub_foe) then
          saccepted=0.0_dp
       end if
       ! ars: we haven't done any NGWF iteration yet
       ! ars: perform the initial step
       ! ja531-> Build the Fermi-Dirac Density Kernel

       ! JA: For the next step we want all Hams to be the same.
       do is = 1, pub_num_spins
          call sparse_embed_copy(initial_ham(is), trial_ham(is))
       end do

       ientry = 0

       ! JA: Make density kernel from Ham. t0 was already initialised to 0.
       call edft_evaluate_step(t0, trial_ham, initial_ham, ham, rep, edft, &
            mdl, hfxstate, denskern, hub, hub_proj_basis, lhxc_fine, ngwf_basis,  &
            paw_sphere_energies, lhxc_energy, hubbard_energy,      &
            e0, s0, l0, foe_line_search_entropy, foe_entropy_approx,    &
            bounds, fermi_tol, tmp_fermi, denswork, tmp_s, foe_mode, ientry)

       do is = 1, pub_num_spins
          call sparse_embed_trace(edft%integrated_ne(is),denskern%m(is,PUB_1K),&
               rep%overlap)
       end do

       ! ars: rescale density kernel
       ! JA: Changed to s_occ on restart because this is definitely set up.
       t_occ(:,PUB_1K) = edft%s_occ(:)
       if (pub_debug_on_root)write(stdout,*) &
               &"DEBUG: EDFT: Rescaling denskern to ",t_occ(:,PUB_1K)
       call kernel_rescale_spam3(denskern, rep%overlap, t_occ, &
            silent=.true.)
       do is = 1, pub_num_spins
          call sparse_embed_trace(edft%integrated_ne(is),denskern%m(is,PUB_1K),&
               rep%overlap)
       end do

       if (pub_foe) then
          saccepted = s0
       end if
    end if ! edft%initialised

    ! ars: Update density-dependent terms with current NGWFs and K.
    ! ars: Calculate energy.
    ! ars: The Hamiltonian is not updated.

    call hamiltonian_dens_dep_matrices(ham, lhxc_fine, eaccepted, &
         lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
         ngwf_basis, hub_proj_basis, hub, denskern, mdl, hfxstate, &
         ham_update=.false., lhxc_fixed=.false.)

    ! ars: update entropy from eigenvalues
    if(.not.pub_foe) then
       saccepted = edft_fermidirac_entropy(edft%occ, edft%num, edft%nbands)
    end if

    ! ars: validate lagrangian
    laccepted = edft_lagrangian(edft, eaccepted, saccepted)

    if(pub_foe) then
       if(laccepted<edft%free_energy) then
          edft%fermi(1:pub_num_spins)=tmp_fermi(1:pub_num_spins)
       else
          tmp_fermi(1:pub_num_spins)=edft%fermi(1:pub_num_spins)
       end if
    end if

    ! gab: Set ientry used for the Pulay Mixing history
    ientry = 1

    diis_size = pub_edft_ham_diis_size

    ! gab: Set number of safety line searches used.
    n_safety = 0
    shifted = .false.

    ! ars: set flag for initialisation
    edft%initialised = .true.
    !========================== END INITIALISATION ============================!

    start_time=wrappers_etime()

    ! ars: loop over occupation numbers
    edft_loop: do edft_iter = 0, current_edft_maxit
       ! jd: Globally accessible copy of iteration number, only used for contro-
       !     lling devel_code dumps of grid quantities from various modules
       pub_inner_loop_iteration = edft_iter

       if(pub_foe.and..not.foe_line_search_entropy) then
          saccepted=0.0_dp
          s1=0.0_dp
          s2=0.0_dp
          s3=0.0_dp
          s4=0.0_dp
          s5=0.0_dp
          s6=0.0_dp

          if(foe_mode=="DENSE") then
             do is=1, pub_num_spins
                s0=0.0_dp
                call dense_convert(denswork(2),rep%overlap)
                call dense_convert(denswork(3),denskern%m(is,PUB_1K))
                call calculate_fermi_entropy(denswork(3),denswork(2),s0)
                saccepted=saccepted+s0
             end do
          else
             do is=1, pub_num_spins
                s0=0.0_dp
                call calculate_fermi_entropy(denskern%m(is,PUB_1K),rep%overlap,s0)
                saccepted=saccepted+s0
             end do
          end if
          saccepted = pub_spin_fac*pub_edft_smearing_width*saccepted

          laccepted = edft_lagrangian(edft, eaccepted, saccepted)

          s1=saccepted
          s2=saccepted
          s3=saccepted
          s4=saccepted
          s5=saccepted
          s6=saccepted
       end if

       ! ars: save the last line search (or init) Hamiltonian as initial_ham
       do is = 1, pub_num_spins
          call sparse_embed_copy(initial_ham(is), trial_ham(is))
          if(pub_foe) then
             call sparse_embed_copy(edft%old_ham(is),trial_ham(is))
          end if
       end do

       ! ars: Avoid updating density-dependent terms with current NGWFs and K.
       !      The density was updated using the last, accepted density kernel.
       ! ars: Calculate accepted_energy - identical to the last accepted_energy.
       ! ars: Update the Hamiltonian.
       call hamiltonian_dens_dep_matrices(ham, lhxc_fine, eaccepted, &
            lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
            ngwf_basis, hub_proj_basis, hub, denskern, &
            mdl, hfxstate, .true., .true.)
       call hamiltonian_build_matrix(ham, rep)

       ! ars: save last accepted values as initial
       e0 = eaccepted
       s0 = saccepted
       l0 = laccepted

       !====================================================================!
       !                     CHECK CONVERGENCE AND PRINT                    !
       !====================================================================!


       ! ja531 -> Remove gradient criterion by default, since it's wrong and can diverge.
       if(pub_edft_rms_gradient_thres.gt.0.0_dp) then
          ! ars: calculate gradient wrt occupancies
          rms_gradient = edft_rms_grad(initial_ham, ham%ham)
       end if

       ! ars: calculate commutator
       commutator = edft_commutator(ham%ham, denskern%m(:,PUB_1K), &
            rep%overlap, rep%inv_overlap)

       ! ars: calculate change in variables
       deltal = laccepted - prev_laccepted
       prev_laccepted = laccepted
       deltae = eaccepted - prev_eaccepted
       prev_eaccepted = eaccepted
       deltas = saccepted - prev_saccepted
       prev_saccepted = saccepted

       ! jd: Take care to avoid premature convergence after step #0.
       !     This happens when some deltas are zero by construction at step #0
       !     and all other criteria are satisfied. This only happens for
       !     deltas and delta_fermi; deltal and deltae are not zero at step #0.
       if(deltas == 0.0_DP .and. edft_iter == 0) then
          deltas = very_big_double
       end if

       if (pub_edft_grand_canonical) then
          delta_integrated_ne(1:pub_num_spins) = &
               edft%integrated_ne(1:pub_num_spins) - &
               prev_integrated_ne(1:pub_num_spins)
          prev_integrated_ne(1:pub_num_spins) = &
               edft%integrated_ne(1:pub_num_spins)
       else
          delta_fermi(1) = edft%fermi(1) - prev_fermi1
          prev_fermi1 = edft%fermi(1)
          if(delta_fermi(1) == 0.0_DP .and. edft_iter == 0) then
             delta_fermi(1) = very_big_double
          end if
          if (pub_num_spins.eq.2) then
             delta_fermi(2) = edft%fermi(2) - prev_fermi2
             prev_fermi2 = edft%fermi(2)
             if(delta_fermi(2) == 0.0_DP .and. edft_iter == 0) then
                delta_fermi(2) = very_big_double
             end if
          end if
       end if

       ! ars: print step on the screen
       if (pub_output_detail>BRIEF) then
          if (pub_edft_grand_canonical) then
             call edft_print_iter(edft_iter, edft,  &
                  eaccepted, saccepted, laccepted, deltal, &
                  rms_gradient, commutator, delta_integrated_ne, taccepted)
          else
             call edft_print_iter(edft_iter, edft,  &
                  eaccepted, saccepted, laccepted, deltal, &
                  rms_gradient, commutator, delta_fermi, taccepted)
          end if
       end if

       ! -----> ars: check convergence and exit if appropriate
       if (pub_edft_grand_canonical) then
          edft_converged = edft_check_convergence(deltae, &
               deltas, deltal, delta_integrated_ne, &
               rms_gradient, commutator, mdl%nat)
       else
          edft_converged = edft_check_convergence(deltae, &
               deltas, deltal, delta_fermi, rms_gradient, &
               commutator, mdl%nat)
       end if

       ! ars: check for increase in the value of the Lagrangian
       if (deltal.gt.0.0_DP) then
          cc = cc + 1
          L_increase(cc) = .true.
       else
          cc = 0
          L_increase(:) = .false.
       end if

       safety_ls = .false.
       if (L_increase(1).and.L_increase(2).and.(pub_edft_trial_step.gt.0)) then
          ! gab: Acts as a safety valve if input lambda value leads to energy increases
          ! for two iterations. Performs a line search to give an optimised lambda
          ! value.

          safety_ls = .true.
          n_safety = n_safety+1

          if (.not.(n_safety.eq.3)) then
             if (pub_on_root) then
                write(stdout,'(a,i3,a)') &
                     'WARNING: Consecutive energy increases have occured ', n_safety, ' time(s).'
                write(stdout,'(a)') &
                     "WARNING: If this keeps occuring, your input lambda value is likely too high."
                write(stdout,'(a)') &
                     "WARNING: Performing line search to improve convergence."
             end if
          end if

          if (n_safety.eq.3) then
             if (pub_on_root) then
                write(stdout,'(a,i3,a)') &
                     'WARNING: Consecutive energy increases have occured ', n_safety, ' times.'
                write(stdout,'(a)') &
                     "WARNING: Your input lambda value is likely too high for reliable convergence."
                write(stdout,'(a)') &
                     "WARNING: Exiting inner loop."
             end if
             exit edft_loop
          end if

       end if

       if (L_increase(1).and.L_increase(2).and.L_increase(3)) then
          ! ars: This is the third time in a row that the Lagrangian
          !      increases. We stop the calculation honorably.
          if (pub_on_root) &
               write(stdout,'(a)') &
               "WARNING: EDFT could not optimise the occupancies any further"
          exit edft_loop
       end if

       ! ars: print convergence message if finished
       if (edft_converged) then
          if (pub_on_root.and.pub_output_detail.ge.NORMAL) &
               write(stdout,'(a38,i4,a12,/)') &
               "EDFT occupancies loop converged after ",&
               edft_iter, " iterations."
          exit edft_loop
       end if

       ! ars: exit if we exceeded the maximum number of iterations
       if (edft_iter.ge.current_edft_maxit) then
          if (pub_on_root.and.pub_output_detail.ge.NORMAL) &
               write(stdout,'(a,i4,a/)') &
               'Finished EDFT density kernel iterations (',edft_iter, ')'
          exit edft_loop
       endif

       ! ars: print line-search step if verbose
       if (pub_output_detail.ge.VERBOSE) then
          if (pub_on_root) then
             write(stdout,'(/,a)') &
                  "  Line search of the free energy minimum:"
             write(stdout,'(a)') &
                  "  Step_type        Step             &
                  &Free_energy   Predicted_free_energy"
          end if
          call edft_print_line_search('  Initial', t0, l0)
       end if

       !------------------------------------------------------------------------
       !--------------------------- EDFT LINE SEARCH ---------------------------
       !------------------------------------------------------------------------

       ! gab: Selects method to construct the new residual for H_{n+1} = H_n + \alpha R[H_n]
       select case (pub_edft_update_scheme)

       case ('PULAY_MIX')

          ! gab: Performs Pulay mixing (H_in[n+1] = \sum_{n} sum_j^i c_j*H_in[j] + c_j*R(H_in[n](j))

          ! gab: Generates residual history and history of input hamiltonians.
          do is = 1, pub_num_spins

             call sparse_embed_copy(edft%ham_in_his(ientry,is), initial_ham(is))

             call sparse_embed_copy(edft%ham_err_his(ientry,is), ham%ham(is))
             call sparse_embed_axpy(edft%ham_err_his(ientry,is), initial_ham(is), -1.0_DP)
          end do

       case ('DAMP_FIXPOINT')

          ! gab: Standard line search routine - continue as normal

       case default

          call utils_abort("EDFT mixing scheme not recogised.")

       end select

       ! gab: if optimum mixing parameter already specified by user - ignore line search.
       ! The first steps tend to be larger, so we allow the line search to use the line
       ! search optimised values for the first few steps. Performs a line search if
       ! energy increases twice for specified lambda value
       if ((pub_edft_trial_step.gt.0).and.(edft_iter.gt.1) &
            .and.(.not.(safety_ls))) then
          t1 = pub_edft_trial_step

          call edft_evaluate_step(t1, trial_ham, initial_ham, ham, rep, edft, &
               mdl, hfxstate, denskern, hub, hub_proj_basis, lhxc_fine, ngwf_basis, &
               paw_sphere_energies, lhxc_energy, hubbard_energy,      &
               e1, s1, l1, foe_line_search_entropy, foe_entropy_approx,    &
               bounds, fermi_tol, tmp_fermi, denswork, tmp_s, foe_mode, ientry)

          if (pub_foe) then
             if(l1 <= minval([l0,l1]))then
                edft%fermi(1:pub_num_spins)=tmp_fermi(1:pub_num_spins)
             else
                tmp_fermi(1:pub_num_spins)=edft%fermi(1:pub_num_spins)
             end if
          end if

          if (pub_output_detail.ge.VERBOSE) &
               call edft_print_line_search('   Input Trial', t1, l1)


          eaccepted = e1; saccepted=s1; laccepted=l1; taccepted=t1

       else ! gab: Use default line search for mixing parameter.


             !! FIRST TRIAL STEP !!
             t1 = taccepted
             call edft_evaluate_step(t1, trial_ham, initial_ham, ham, rep, edft, &
                  mdl, hfxstate, denskern, hub, hub_proj_basis, lhxc_fine, ngwf_basis, &
                  paw_sphere_energies, lhxc_energy, hubbard_energy,      &
                  e1, s1, l1, foe_line_search_entropy, foe_entropy_approx,    &
                  bounds, fermi_tol, tmp_fermi, denswork, tmp_s, foe_mode, ientry)

             if (pub_foe) then
                if(l1<l0) then
                   edft%fermi(1:pub_num_spins)=tmp_fermi(1:pub_num_spins)
                else
                   tmp_fermi(1:pub_num_spins)=edft%fermi(1:pub_num_spins)
                end if
             end if

             if (pub_output_detail.ge.VERBOSE) &
                  call edft_print_line_search('   Trial1', t1, l1)

             !! SECOND TRIAL STEP !!
             npts = 2
             xx2 = (/ t0, t1 /); yy2 = (/ l0, l1 /)
             t2 = edft_trial_step(npts, xx2, yy2)

             call edft_evaluate_step(t2, trial_ham, initial_ham, ham, rep, edft, &
                  mdl, hfxstate, denskern, hub, hub_proj_basis, lhxc_fine, ngwf_basis,  &
                  paw_sphere_energies, lhxc_energy, hubbard_energy,      &
                  e2, s2, l2, foe_line_search_entropy, foe_entropy_approx,    &
                  bounds, fermi_tol, tmp_fermi, denswork, tmp_s, foe_mode, ientry)

             if (pub_foe) then
                if(l2 <= minval([l0,l1]))then
                   edft%fermi(1:pub_num_spins)=tmp_fermi(1:pub_num_spins)
                else
                   tmp_fermi(1:pub_num_spins)=edft%fermi(1:pub_num_spins)
                end if
             end if

             if (pub_output_detail.ge.VERBOSE) &
                  call edft_print_line_search('   Trial2', t2, l2)

             !! QUADTRATIC STEP !!

             ! set quadtratic step
             call edft_fit_second_order_poly_3points(quadratic_success, &
                  t3, predicted_l3, t0, t1, t2, l0, l1, l2)

             ! validate quadratic step
             if ( (.not.quadratic_success).or.(t3.lt.0.0_DP).or.&
                  (t3.gt.pub_edft_max_step) ) then
                quadratic_success = .false.
                ! ars: take trial step
                npts=3; xx3=(/ t0, t1, t2 /); yy3=(/ l0, l1, l2 /)
                t3=edft_trial_step(npts, xx3, yy3)
             end if

             call edft_evaluate_step(t3, trial_ham, initial_ham, ham, rep, edft, &
                  mdl, hfxstate, denskern, hub, hub_proj_basis, lhxc_fine, ngwf_basis, &
                  paw_sphere_energies, lhxc_energy, hubbard_energy,      &
                  e3, s3, l3, foe_line_search_entropy, foe_entropy_approx,    &
                  bounds, fermi_tol, tmp_fermi, denswork, tmp_s, foe_mode, ientry)

             if(pub_foe) then
                if(l3<=minval([l0,l1,l2])) then
                   edft%fermi(1:pub_num_spins)=tmp_fermi(1:pub_num_spins)
                else
                   tmp_fermi(1:pub_num_spins)=edft%fermi(1:pub_num_spins)
                end if
             end if

             if ( (quadratic_success).and.(l3.lt.l0).and.&
                  (l3.lt.l1).and.(l3.lt.l2) ) then
                quadratic_accepted = .true.
                eaccepted = e3; saccepted=s3; laccepted=l3; taccepted=t3
             else
                quadratic_accepted = .false.
             end if

             ! ars: print line search if verbose
             if (pub_output_detail.ge.VERBOSE) then
                if (quadratic_success) then
                   call edft_print_line_search('Quadratic', t3, l3, predicted_l3)
                else
                   call edft_print_line_search('   Trial3', t3, l3)
                end if
             end if

             if (.not. quadratic_accepted) then
                !! CUBIC STEP !!

                ! ars: set cubic step
                call edft_fit_third_order_poly_4points(cubic_success, &
                     t4, predicted_l4, t0, t1, t2, t3, l0, l1, l2, l3)
                ! ars: validate cubic step
                if ( (.not.cubic_success).or.(t4.lt.0.0_DP).or.&
                     (t4.gt.pub_edft_max_step) ) then
                   cubic_success = .false.
                   ! ars: take trial step
                   npts=4; xx4=(/ t0, t1, t2, t3 /); yy4=(/ l0, l1, l2, l3 /)
                   t4=edft_trial_step(npts, xx4, yy4)
                end if

                call edft_evaluate_step(t4, trial_ham, initial_ham, ham, rep, edft, &
                     mdl, hfxstate, denskern, hub, hub_proj_basis, lhxc_fine, ngwf_basis, &
                     paw_sphere_energies, lhxc_energy, hubbard_energy,      &
                     e4, s4, l4, foe_line_search_entropy, foe_entropy_approx,    &
                     bounds, fermi_tol, tmp_fermi, denswork, tmp_s, foe_mode, ientry)

                if(pub_foe) then
                   if(l4<=minval([l0,l1,l2,l3])) then
                      edft%fermi(1:pub_num_spins)=tmp_fermi(1:pub_num_spins)
                   else
                      tmp_fermi(1:pub_num_spins)=edft%fermi(1:pub_num_spins)
                   end if
                end if

                ! ars: validate cubic_step
                if ( (cubic_success).and.(l4.lt.l0).and.(l4.lt.l1).and.&
                     (l4.lt.l2).and.(l4.lt.l3) ) then
                   cubic_accepted = .true.
                   eaccepted = e4; saccepted=s4; laccepted=l4; taccepted=t4
                else
                   cubic_accepted = .false.
                end if

                ! ars: print line search if verbose
                if (pub_output_detail.ge.VERBOSE) then
                   if (cubic_success) then
                      call edft_print_line_search('    Cubic', t4, l4, predicted_l4)
                   else
                      call edft_print_line_search('   Trial4', t4, l4)
                   end if
                end if

                if (.not. cubic_accepted) then
                   !! QUARTIC STEP !!

                   ! ars: set quartic step
                   call edft_fit_fourth_order_poly_5points(quartic_success, &
                        t5, predicted_l5, t0, t1, t2, t3, t4, l0, l1, l2, l3, l4)
                   if(t5/=t5) then
                      if(pub_on_root) then
                         !                   write(*,*) "ja531 --> fit fourth order NaN : ",quartic_success, &
                         !                        & t5, predicted_l5, t0, t1, t2, t3, t4, l0, l1, l2, l3, l4
                      end if
                   end if
                   ! ars: validate quartic step
                   if ( (.not.quartic_success).or.(t5.lt.0.0_DP).or.&
                        (t5.gt.pub_edft_max_step) ) then
                      quartic_success = .false.
                      ! ars: set step to tiny value
                      t5 = 0.01_DP*min(t1,t2,t3,t4)
                   end if

                   call edft_evaluate_step(t5, trial_ham, initial_ham, ham, rep, edft, &
                        mdl, hfxstate, denskern, hub, hub_proj_basis, lhxc_fine, ngwf_basis, &
                        paw_sphere_energies, lhxc_energy, hubbard_energy,      &
                        e5, s5, l5, foe_line_search_entropy, foe_entropy_approx,    &
                        bounds, fermi_tol, tmp_fermi, denswork, tmp_s, foe_mode, ientry)

                   if(pub_foe) then
                      if(l5<=minval([l0,l1,l2,l3,l4])) then
                         edft%fermi(1:pub_num_spins)=tmp_fermi(1:pub_num_spins)
                      else
                         tmp_fermi(1:pub_num_spins)=edft%fermi(1:pub_num_spins)
                      end if
                   end if

                   ! ars: validate quartic_step
                   if ( (quartic_success).and.(l5.lt.l0).and.(l5.lt.l1).and.&
                        (l5.lt.l2).and.(l5.lt.l3).and.(l5.lt.l4) ) then
                      quartic_accepted = .true.
                      eaccepted = e5; saccepted=s5; laccepted=l5; taccepted=t5
                   else
                      quartic_accepted = .false.
                   end if

                   ! ars: print line search if verbose
                   if (pub_output_detail.ge.VERBOSE) then
                      if (quartic_success) then
                         call edft_print_line_search('  Quartic',t5,l5,predicted_l5)
                      else
                         call edft_print_line_search('   Trial5',t5,l5)
                      end if
                   end if

                   tiny_accepted = .false.
                   if ( (.not.quartic_accepted).and.(quartic_success).and.&
                        (l0.lt.l1).and.(l0.lt.l2).and.(l0.lt.l3).and.&
                        (l0.lt.l4).and.(l0.lt.l5) ) then
                      !! TINY STEP !!

                      ! ars: set step to a tiny value
                      t6 = 0.01_DP*min(t1,t2,t3,t4,t5)

                      call edft_evaluate_step(t6, trial_ham, initial_ham, ham, rep, edft, &
                           mdl, hfxstate, denskern, hub, hub_proj_basis, lhxc_fine, ngwf_basis, &
                           paw_sphere_energies, lhxc_energy, hubbard_energy,      &
                           e6, s6, l6, foe_line_search_entropy, foe_entropy_approx,    &
                           bounds, fermi_tol, tmp_fermi, denswork, tmp_s, foe_mode, ientry)

                      if(pub_foe) then
                         if(l6<=minval([l0,l1,l2,l3,l4,l5])) then
                            edft%fermi(1:pub_num_spins)=tmp_fermi(1:pub_num_spins)
                         else
                            tmp_fermi(1:pub_num_spins)=edft%fermi(1:pub_num_spins)
                         end if
                      end if

                      ! ars: validate quartic_step
                      if ( (l6.lt.l0).and.(l6.lt.l1).and.(l6.lt.l2).and.&
                           (l6.lt.l3).and.(l6.lt.l4).and.(l6.lt.l5) ) then
                         tiny_accepted = .true.
                         eaccepted = e6; saccepted=s6; laccepted=l6; taccepted=t6
                      else
                         tiny_accepted = .false.
                      end if

                      ! ars: print line search if verbose
                      if (pub_output_detail.ge.VERBOSE) then
                         call edft_print_line_search('   Trial6', t6, l6)
                      end if

                   end if ! .not.quartic_accepted etc
                end if ! .not.cubic_accepted
             end if ! .not.quadratic accepted

             ! kkbd: Choose best of trial steps if none accepted
             if ( (.not.tiny_accepted).and.(.not.quartic_accepted).and.&
                  (.not.cubic_accepted).and.(.not.quadratic_accepted) ) then

                ! ars: set best step
                if      ( (l6.lt.l1).and.(l6.lt.l2).and.&
                     (l6.lt.l3).and.(l6.lt.l4).and.(l6.lt.l5) ) then
                   taccepted = t6
                else if ( (l5.lt.l1).and.(l5.lt.l2).and.&
                     (l5.lt.l3).and.(l5.lt.l4).and.(l5.lt.l6) ) then
                   taccepted = t5
                else if ( (l4.lt.l1).and.(l4.lt.l2).and.&
                     (l4.lt.l3).and.(l4.lt.l5).and.(l4.lt.l6) ) then
                   taccepted = t4
                else if ( (l3.lt.l1).and.(l3.lt.l2).and.&
                     (l3.lt.l4).and.(l3.lt.l5).and.(l3.lt.l6) ) then
                   taccepted = t3
                else if ( (l2.lt.l1).and.(l2.lt.l3).and.&
                     (l2.lt.l4).and.(l2.lt.l5).and.(l2.lt.l6) ) then
                   taccepted = t2
                else if ( (l1.lt.l2).and.(l1.lt.l3).and.&
                     (l1.lt.l4).and.(l1.lt.l5).and.(l1.lt.l6) ) then
                   taccepted = t1
                end if

                call edft_evaluate_step(taccepted, trial_ham, initial_ham, ham, rep, &
                     edft, mdl, hfxstate, denskern, hub, hub_proj_basis, lhxc_fine, ngwf_basis, &
                     paw_sphere_energies, lhxc_energy, hubbard_energy, eaccepted, &
                     saccepted, laccepted, foe_line_search_entropy, foe_entropy_approx,&
                     bounds, fermi_tol, tmp_fermi, denswork, tmp_s, foe_mode, ientry)

                if(pub_foe) then
                   !       if(laccepted<=minval([l0,l1,l2,l3,l4,l5,l6])) then
                   edft%fermi(1:pub_num_spins)=tmp_fermi(1:pub_num_spins)
                   !       else
                   !          tmp_fermi(1:pub_num_spins)=edft%fermi(1:pub_num_spins)
                   !       end if
                end if

                ! ars: print line search if verbose
                if (pub_output_detail.ge.VERBOSE) &
                     call edft_print_line_search(' NoMinFnd', taccepted, laccepted)

             end if ! no steps accepted

          end if ! input trial step specified

          ! gab: shift arrays if we have reached pub_edft_ham_diis_size iterations
          if ((edft_iter).ge.diis_size) then
             shifted = .true.
             call edft_diis_shift(edft%ham_in_his, diis_size, pub_num_spins)
             call edft_diis_shift(edft%ham_err_his, diis_size, pub_num_spins)
          end if

          ! gab: find next position in history to work with
          call edft_diis_find_ientry(ientry, edft_iter, shifted)

    end do edft_loop


    !-------------------------------------------------------------------------
    ! ars: terminate

    ! ars: set final values
    edft%free_energy = laccepted
    edft%energy = eaccepted
    edft%entropy = saccepted

    ! ars: print QC test output
    if (pub_print_qc) call edft_print_qc(edft, commutator, rms_gradient, &
         deltal, deltae, deltas, edft_iter)

    ! ars: save the density kernel and Hamiltonian for future restarts
    if(.not.loc_force_no_IO) then
       if (.not.pub_write_converged_dk_ngwfs) then
          if (pub_write_denskern) call restart_kernel_write(denskern)
          ! JA: Changed to initial_ham - otherwise restarts aren't right!
          if (pub_write_hamiltonian) call restart_hamiltonian_write(initial_ham)
       end if
    end if

    ! ars: destroy structures
    do is = 1, pub_num_spins
       call sparse_embed_destroy(initial_ham(is),arrname='initial_ham')
       call sparse_embed_destroy(trial_ham(is),arrname='trial_ham')
    end do

    if (pub_foe) then
       if(foe_mode=="DENSE") then
          call dense_destroy(denswork(1))
          call dense_destroy(denswork(2))
          call dense_destroy(denswork(3))
          call dense_destroy(denswork(4))
       end if
    end if

    ! ars: deallocate structures
    deallocate(initial_ham, stat=ierr)
    call utils_dealloc_check('edft_calculate', 'initial_ham', ierr)
    deallocate(trial_ham, stat=ierr)
    call utils_dealloc_check('edft_calculate', 'trial_ham', ierr)

    ! kkbd: Reduce spin_fix timer by 1 if it's > 0
    if (pub_edft_spin_fix > 0) then
       if (pub_on_root .and. (pub_edft_spin_fix == 1)) write(stdout,*) &
               &"EDFT: Fixed-spin phase complete. Freeing spin."
       pub_edft_spin_fix = pub_edft_spin_fix - 1
    end if

    if (pub_edft_spin_fix == 0) then
       if (pub_debug_on_root) write(stdout,*) "DEBUG: EDFT: Setting s_occ = integrated_ne"
       edft%s_occ = edft%integrated_ne

       !kkbd: Default behavior is to retain the current net spin when we reset NGWFs etc.
       !      An alternative would be to reset the net spin to the value in the input file for a
       !      'true' reset. The issue there is if the net spin is quite wrong this can affect the
       !      first few steps of the calculation considerably.
       if (utils_devel_code(.true.,"EDFT","RETAIN_SPIN",pub_devel_code)) then
         if (pub_debug_on_root) write(stdout,*) "kkbd: setting pub_real_spin from ",pub_real_spin
         pub_real_spin = edft%integrated_ne(1) - edft%integrated_ne(2)
         if (pub_debug_on_root) write(stdout,*) "kkbd: setting pub_real_spin to ",pub_real_spin
       end if
    end if

    ! ars: stop timer
    call timer_clock('edft_calculate',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving edft_calculate'

  end subroutine edft_calculate

  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine edft_apply_smearing_operator(mode,  &
       rep, ham, edft, chempot, denskern, &
       search_entropy, entropy, ignore_history, entropy_approx, chemtol, bounds, &
       densework, contracoham, cleanup)
    !==========================================================================!
    ! Call a Fermi Operator Expansion (FOE) method from smearing mod to        !
    ! calculate a Fermi-Dirac occupancy distribution density kernel.           !
    ! The routine tries to make a sensible guess for an initial chemical       !
    ! potential because the chemical potential search (to conserve electron    !
    ! number) is expensive. This is done by assuming that of a history of      !
    ! FOE applications over previous iterations, the best chemical potential   !
    ! to initialise with will be the one associated with the Hamiltonian which !
    ! most closely matches the Hamiltonian matrix of the current iteration.    !
    ! This is achieved by creating an array of Hamiltonian traces and          !
    ! Frobenius norms, together with the optimised chemical potential which    !
    ! conserves electron number for these Hamiltonians. The chemical potential !
    ! associated with the Hamiltonian trace and norm which is closest to those !
    ! of the current Hamiltonian is used to initialise.                        !
    ! !
    ! The routine has a lot of arguments which may become rundat parameters or !
    ! global EDFT mod variables.                                               !
    !--------------------------------------------------------------------------!
    ! Written by Jolyon Aarons in April 2016.                                  !
    !==========================================================================!
    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use dense, only: DEM, dense_convert
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_debug_on_root, pub_edft_smearing_width, &
         pub_num_spins,  pub_num_kpoints, PUB_1K, pub_edft_round_evals, &
         pub_contracoham_radmult, pub_H2denskern_sparsity, pub_edft_spin_fix, &
         pub_foe_test_sparsity, pub_debug, pub_FOE_avoid_inversions
    use smearing_operator, only: apply_fermi_function
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_trace,&
         sparse_embed_max_abs_element, sparse_embed_create, &
         sparse_embed_axpy, sparse_embed_product, sparse_embed_transpose, &
         sparse_embed_transpose_structure, sparse_embed_destroy
    use sparse, only: sparse_solve2
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, utils_assert
   implicit none
    character(len=5),        intent(in)    :: mode ! Apply the operator in SPAM3 or convert to DENSE first and back after.
    type(NGWF_REP),          intent(in)    :: rep ! NGWF representation.
    type(SPAM3_EMBED)   ,          intent(in)    :: ham(pub_num_spins)             ! Hamiltonian matrix.
    type(EDFT_MODEL),        intent(inout) :: edft   ! EDFT eig stuff
    real(kind=dp),           intent(inout) :: chempot(pub_num_spins) ! Chemical potential. Initial (and final) value
    type(SPAM3_EMBED_ARRAY),       intent(inout) :: denskern ! Output density kernel.

    !newopts
    logical,       optional, intent(in)    :: search_entropy ! Whether to calculate entropy.
    real(kind=dp), optional, intent(inout) :: entropy ! Entropy (unchanged if search_entropy is false).
    logical,       optional, intent(in)    :: ignore_history ! Force use of the initial chemical potential in chempot and bounds.
    integer,       optional, intent(in)    :: entropy_approx ! Approximate the entropy? (0: exact, 1: approx, 2:refined approx)
    real(kind=dp), optional, intent(in)    :: chemtol ! Tolerance on Chemical potential search (e.g. 1e-9).
    real(kind=dp), optional, intent(inout) :: bounds(2) ! Lower(1) and Upper(2) deltas giving search window in chempot.

    type(DEM),     optional, intent(inout) :: densework(4) ! Three dense work arrays if in DENSE mode. Must be set up already!
    ! NB. These must be the right size for ham/overlap/denskern!
    logical,       optional, intent(in)    :: contracoham ! Use contra-covariant Ham instead of Lowdin orthog.
    logical,       optional, intent(in)    :: cleanup ! Tells the routine to free all arrays and reset everything.


    real(kind=DP),allocatable :: tmp_fermi(:) ! Temporary storage of chemical potential.

    logical, save :: first_pass=.true. ! Allocate everything we need if this is the first pass.
    real(kind=dp) :: boundl, boundu
    logical :: loc_cleanup
    integer :: ierr
    integer :: iorb
    integer :: spin_channel

    logical       :: loc_search_entropy
    logical       :: loc_ignore_history
    integer       :: loc_entropy_approx
    real(kind=dp) :: loc_chemtol
    real(kind=dp) :: loc_bounds(2)
    logical :: loc_contracoham

    integer :: is
    character(len=30) :: contracodkern_trans_struc, dkern_trans_struc, overlap_trans_struc
    real(kind=dp) :: norm
    type(SPAM3_EMBED), dimension(:), allocatable :: hopping
    type(SPAM3_EMBED), dimension(:), allocatable :: testmat, testmat2
    type(SPAM3_EMBED), dimension(:), allocatable :: contracodkern, contracodkern_trans
    type(SPAM3_EMBED), dimension(:), allocatable :: dkern_trans
    type(SPAM3_EMBED) :: overlap_trans


    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering edft_apply_smearing_operator'
    call timer_clock('edft_apply_smearing_operator',1)

    loc_contracoham=.false.
    if(present(contracoham)) then
       loc_contracoham=contracoham
    end if

    if(mode/="SPAM3".and.mode/="DENSE".and.mode/="EIGEN") then
       call utils_abort('Unknown mode in subroutine edft_apply_smearing_operator')
    end if
    if(mode=="DENSE") then
       if(.not.present(densework)) then
          call utils_abort('Dense mode in edft_apply_smearing_operator requires a densework work array.')
       end if
    end if

    ! jme: KPOINTS_DANGER
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine edft_apply_smearing_operator not yet compatible with more&
         & than one k-point.')

    loc_search_entropy=.false.
    if(present(search_entropy)) then
       loc_search_entropy=search_entropy
    end if
    loc_ignore_history=.true.
    if(present(ignore_history)) then
       loc_ignore_history=ignore_history
    end if
    loc_entropy_approx=0
    if(present(entropy_approx)) then
       loc_entropy_approx=entropy_approx
    end if
    loc_chemtol=1e-9_dp
    if(present(chemtol)) then
       loc_chemtol=chemtol
    end if
    loc_bounds=[-5.0_dp,5.0_dp]
    if(present(bounds)) then
       loc_bounds=bounds
    end if

    if(loc_search_entropy) then
       call utils_assert(present(entropy), &
            'Error in edft_apply_smearing_operator: if optimising free-energy, then &
            &entropy argument must be passed.')
    end if

    loc_cleanup=.false.
    if(present(cleanup)) then
       loc_cleanup=cleanup
    end if

    if(.not.allocated(mu_lookup).or..not.allocated(mu_lookup_idx).or..not.allocated(mu_lookup_len)) then
       first_pass=.true.
    end if

    if(loc_cleanup.or.first_pass) then
       if(allocated(mu_lookup)) then
          deallocate(mu_lookup,stat=ierr)
          call utils_dealloc_check('edft_apply_smearing_operator','mu_lookup',ierr)
       end if
       if(allocated(mu_lookup_idx)) then
          deallocate(mu_lookup_idx,stat=ierr)
          call utils_dealloc_check('edft_apply_smearing_operator','mu_lookup_idx',ierr)
       end if
       if(allocated(mu_lookup_len)) then
          deallocate(mu_lookup_len,stat=ierr)
          call utils_dealloc_check('edft_apply_smearing_operator','mu_lookup_len',ierr)
       end if
    end if
    if(present(cleanup)) then
       if(cleanup) then
          first_pass=.true.
          call timer_clock('edft_apply_smearing_operator',2)
          if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving edft_apply_smearing_operator'
          return
       end if
    end if

    if(first_pass) then
       if(allocated(mu_lookup)) then
          deallocate(mu_lookup,stat=ierr)
          call utils_dealloc_check('edft_apply_smearing_operator','mu_lookup',ierr)
       end if
       if(allocated(mu_lookup_idx)) then
          deallocate(mu_lookup_idx,stat=ierr)
          call utils_dealloc_check('edft_apply_smearing_operator','mu_lookup_idx',ierr)
       end if
       if(allocated(mu_lookup_len)) then
          deallocate(mu_lookup_len,stat=ierr)
          call utils_dealloc_check('edft_apply_smearing_operator','mu_lookup_len',ierr)
       end if

       allocate(mu_lookup(pub_mu_lookup_maxlen,3,pub_num_spins),stat=ierr)
       call utils_alloc_check('edft_apply_smearing_operator','mu_lookup',ierr)
       allocate(mu_lookup_idx(pub_num_spins),stat=ierr)
       call utils_alloc_check('edft_apply_smearing_operator','mu_lookup_idx',ierr)
       allocate(mu_lookup_len(pub_num_spins),stat=ierr)
       call utils_alloc_check('edft_apply_smearing_operator','mu_lookup_len',ierr)

       ! Initialise so that we start on the first entry after the cyclic index increment.
       mu_lookup_idx=pub_mu_lookup_maxlen
       mu_lookup_len=0

    end if


    allocate(tmp_fermi(pub_num_spins),stat=ierr)
    call utils_alloc_check('edft_apply_smearing_operator','tmp_fermi',ierr)

    do spin_channel = 1, pub_num_spins
       ! Update the current index and the total number of entries in the mu_lookup array.
       mu_lookup_idx(spin_channel)=mod((mu_lookup_idx(spin_channel)),pub_mu_lookup_maxlen)+1
       mu_lookup_len(spin_channel)=min(mu_lookup_len(spin_channel)+1,pub_mu_lookup_maxlen)
       ! Add the trace and F.norm of this Hamiltonian to the mu_lookup array.
       call sparse_embed_trace( mu_lookup(mu_lookup_idx(spin_channel),mu_lookup_trace,spin_channel), &
            ham(spin_channel) )
       call sparse_embed_trace( mu_lookup(mu_lookup_idx(spin_channel),mu_lookup_norm,spin_channel), &
            ham(spin_channel), ham(spin_channel) )
       mu_lookup(mu_lookup_idx(spin_channel),mu_lookup_norm,spin_channel) = &
            sqrt(mu_lookup(mu_lookup_idx(spin_channel),mu_lookup_norm,spin_channel))
    end do ! spin_channel

    ! kkbd TODO are we treating these bounds correctly?
    boundl=loc_bounds(1)
    boundu=loc_bounds(2)

    ! ja531 -> Fix for the case when we have more than one spin channel.
    ! kkbd - no longer required
    !if(mu_lookup_len(spin_channel)<2) then
    !   first_pass=.true.
    !end if

    ! Find the closest (least-squares) Hamiltonian in terms of trace and norm that
    ! has already been calculated to this Hamiltonian.
    ! Use the chemical potential associated with the old Hamiltonian as the starting guess
    ! for the chemical potential associated with this Hamiltonian.
    if(.not.first_pass) then
       do spin_channel = 1, pub_num_spins
          if(mu_lookup_len(spin_channel)<pub_mu_lookup_maxlen) then
             tmp_fermi(spin_channel) = mu_lookup(minloc(sqrt(&
                  (mu_lookup(1:mu_lookup_idx(spin_channel)-1,mu_lookup_trace,spin_channel)-&
                  mu_lookup(mu_lookup_idx(spin_channel), &
                  mu_lookup_trace,spin_channel))**2 +(mu_lookup(1:mu_lookup_idx(spin_channel)-1,&
                  mu_lookup_norm,spin_channel)-mu_lookup(mu_lookup_idx(spin_channel), &
                  mu_lookup_norm,spin_channel))**2),1), mu_lookup_mu, spin_channel)

             if(mu_lookup_idx(spin_channel)-1>1) then
                if(mu_lookup_idx(spin_channel)-1==2) then
                   ! This is just a bit of a contingency in case we end up with degenerate mu's when we've
                   ! only operated twice... might not be needed, but it shouldn't hurt.
                   if(abs(maxval([mu_lookup(1:mu_lookup_idx(spin_channel)-1,&
                        mu_lookup_mu,spin_channel)],1)-&
                        minval(mu_lookup(1:mu_lookup_idx(spin_channel)-1,mu_lookup_mu,spin_channel),1))*&
                        delta_mu_fac>1e-4) then
                      boundu=abs(maxval([mu_lookup(1:mu_lookup_idx(spin_channel)-1,&
                           mu_lookup_mu,spin_channel)],1)-&
                           minval(mu_lookup(1:mu_lookup_idx(spin_channel)-1,mu_lookup_mu,spin_channel),1))*&
                           delta_mu_fac
                   end if
                else
                   boundu=abs(maxval(mu_lookup(1:mu_lookup_idx(spin_channel)-1,&
                        mu_lookup_mu,spin_channel),1)-&
                        minval(mu_lookup(1:mu_lookup_idx(spin_channel)-1,mu_lookup_mu,spin_channel),1))*&
                        delta_mu_fac
                end if
             end if

          else
             tmp_fermi(spin_channel) = mu_lookup(minloc(&
                  [sqrt((mu_lookup(1:mu_lookup_idx(spin_channel)-1,&
                  mu_lookup_trace,spin_channel)-mu_lookup(mu_lookup_idx(spin_channel), &
                  mu_lookup_trace,spin_channel))**2 +(mu_lookup(1:mu_lookup_idx(spin_channel)-1,&
                  mu_lookup_norm,spin_channel)-mu_lookup(mu_lookup_idx(spin_channel), &
                  mu_lookup_norm,spin_channel))**2),huge(tmp_fermi(spin_channel)),&
                  sqrt((mu_lookup(mu_lookup_idx(spin_channel)+1:pub_mu_lookup_maxlen,&
                  mu_lookup_trace,spin_channel)-mu_lookup(mu_lookup_idx(spin_channel), &
                  mu_lookup_trace,spin_channel))**2 + &
                  (mu_lookup(mu_lookup_idx(spin_channel)+1:pub_mu_lookup_maxlen, &
                  mu_lookup_norm,spin_channel)-mu_lookup(mu_lookup_idx(spin_channel), &
                  mu_lookup_norm,spin_channel))**2)],1),mu_lookup_mu, spin_channel)

             boundu=abs(maxval([mu_lookup(1:mu_lookup_idx(spin_channel)-1,&
                  mu_lookup_mu,spin_channel),&
                  mu_lookup(mu_lookup_idx(spin_channel)+1:pub_mu_lookup_maxlen,&
                  mu_lookup_mu,spin_channel)],1)-&
                  minval([mu_lookup(1:mu_lookup_idx(spin_channel)-1,mu_lookup_mu,spin_channel),&
                  mu_lookup(mu_lookup_idx(spin_channel)+1:pub_mu_lookup_maxlen,&
                  mu_lookup_mu,spin_channel)],1))*delta_mu_fac

          end if
          boundl=-boundu
          boundu=boundu+tmp_fermi(spin_channel)
          boundl=boundl+tmp_fermi(spin_channel)
       end do ! pub_num_spins
    end if

    first_pass=.false.

    ! If requested, just ignore the value we just found and use the one supplied in chempot.
    if(loc_ignore_history) then
       do spin_channel = 1, pub_num_spins
          tmp_fermi(spin_channel)=chempot(spin_channel)
          boundl=loc_bounds(1)
          boundu=loc_bounds(2)
       end do ! spin_channel
    end if

    ! Perform the FOE calculation.
    if(mode=="DENSE") then
       do spin_channel = 1, pub_num_spins
          ! Assuming densework has been set up correctly!!!
          ! Putting ham in densework(1), overlap in densework(2) and denskern in densework(3).
          call dense_convert(densework(1),ham(spin_channel))
          call dense_convert(densework(2),rep%overlap)
          call dense_convert(densework(3),denskern%m(spin_channel,PUB_1K))
          call dense_convert(densework(4),rep%inv_overlap)
          if(loc_search_entropy) then
             call apply_fermi_function(densework(1),densework(2),densework(4),&
                  tmp_fermi(spin_channel), &
                  pub_edft_smearing_width,densework(3), &
                  correct_chemical_potential=.true., &
                  Ne=edft%s_occ(spin_channel), &
                  boundl=boundl, boundu=boundu, &
                  fermi_tol=loc_chemtol, entropy=entropy, &
                  entropy_approx=loc_entropy_approx)


          else
             call apply_fermi_function(densework(1),densework(2),densework(4),&
                  tmp_fermi(spin_channel), &
                  pub_edft_smearing_width,densework(3), &
                  correct_chemical_potential=.true., &
                  Ne=edft%s_occ(spin_channel), &
                  boundl=boundl, boundu=boundu, &
                  fermi_tol=loc_chemtol)
          end if

          call dense_convert(denskern%m(spin_channel,PUB_1K),densework(3))
       end do ! spin_channel

    else if(mode=="SPAM3") then

       if(pub_H2denskern_sparsity) then
          ! ja531-> set up hopping matrix
          allocate(hopping(pub_num_spins),stat=ierr)
          call utils_alloc_check('edft_apply_smearing_operator','hopping',ierr)
          allocate(testmat(pub_num_spins),stat=ierr)
          call utils_alloc_check('edft_apply_smearing_operator','testmat',ierr)
          allocate(testmat2(pub_num_spins),stat=ierr)
          call utils_alloc_check('edft_apply_smearing_operator','testmat2',ierr)

          ! ja531-> set up contracodenskern
          allocate(contracodkern(pub_num_spins),stat=ierr)
          call utils_alloc_check('edft_apply_smearing_operator','contracodkern',ierr)


          spin_loop: do spin_channel = 1, pub_num_spins
             if(abs(pub_contracoham_radmult-1.0_dp)>epsilon(1.0_dp)) then
                hopping(spin_channel)%structure = 'ccH'
                call sparse_embed_create(hopping(spin_channel),arrname='hopping')
             else
                call sparse_embed_create(hopping(spin_channel),ham(spin_channel),arrname='hopping')
             end if

             ! ja531-> The "hopping matrix" is the contra-covariant Hamiltonian matrix.
             ! ja531-> we form this by solving S_{\gamma\alpha}P^{\alpha}_{\beta} = H_{\gamma\beta}
             ! ja531-> for the hopping matrix, P. We could also do this by inverting S and multiplying.
             ! rc2013: EMBED_FIX -- this enforces 1 subsystem
             call sparse_solve2(hopping(spin_channel)%p, &
                  ham(spin_channel)%p, rep%overlap%p, &
                  aplusb_sparsity=.false.)

             if(pub_debug.or.pub_foe_test_sparsity) then
                call sparse_embed_create(testmat(spin_channel),ham(spin_channel),arrname='testmat')
                testmat2(spin_channel)%structure = 'K'
                call sparse_embed_create(testmat2(spin_channel),arrname='testmat2')

                call sparse_embed_product(testmat(spin_channel),rep%overlap,hopping(spin_channel))
                call sparse_embed_axpy(testmat(spin_channel),ham(spin_channel),-1.0_dp)
                !       norm=sparse_entrywspin_channele_norm(testmat(spin_channel),2)
                norm=sparse_embed_max_abs_element(testmat(spin_channel))
                if(pub_debug_on_root) then
                   write(stdout,*) "Norm of delta matrix for approximate contra-covariant Ham: ",norm,&
                        &" spin channel:",spin_channel
                end if

                call sparse_embed_product(testmat2(spin_channel),rep%overlap,hopping(spin_channel))
                call sparse_embed_axpy(testmat2(spin_channel),ham(spin_channel),-1.0_dp)
                !       norm=sparse_entrywise_norm(testmat2(spin_channel),2)
                norm=sparse_embed_max_abs_element(testmat2(spin_channel))
                if(pub_debug_on_root) then
                   write(stdout,*) "Norm of delta matrix for contra-covariant Ham: ",norm,&
                        &" spin channel:",spin_channel
                end if

                !ja531-> Clean up test mats.
                call sparse_embed_destroy(testmat(spin_channel),arrname='testmat')
                call sparse_embed_destroy(testmat2(spin_channel),arrname='testmat2')
             end if

             call sparse_embed_create(contracodkern(spin_channel), &
                  denskern%m(spin_channel,PUB_1K),arrname='contracodkern')

             ! Use the contracovariant FOE interface
             if(pub_FOE_avoid_inversions) then
                call apply_fermi_function(ham(spin_channel),rep%overlap, &
                     hopping(spin_channel),tmp_fermi(spin_channel), &
                     pub_edft_smearing_width, rep%n_occ(spin_channel,PUB_1K), &
                     entropy,denskern%m(spin_channel,PUB_1K))
             else
                call apply_fermi_function(ham(spin_channel),rep%overlap, &
                     hopping(spin_channel),tmp_fermi(spin_channel), &
                     pub_edft_smearing_width, rep%n_occ(spin_channel,PUB_1K), &
                     entropy,denskern%m(spin_channel,PUB_1K), &
                     inverse_overlap=rep%inv_overlap)
             end if


             !ja531-> clean up hopping matrix
             call sparse_embed_destroy(hopping(spin_channel),arrname='hopping')

          end do spin_loop

          deallocate(hopping,stat=ierr)
          call utils_dealloc_check('edft_apply_smearing_operator','hopping',ierr)
          deallocate(testmat,stat=ierr)
          call utils_dealloc_check('edft_apply_smearing_operator','testmat',ierr)
          deallocate(testmat2,stat=ierr)
          call utils_dealloc_check('edft_apply_smearing_operator','testmat2',ierr)


       else

          do spin_channel = 1, pub_num_spins
             if(loc_search_entropy) then
                call apply_fermi_function(ham(spin_channel),rep%overlap,rep%inv_overlap,&
                     tmp_fermi(spin_channel), &
                     pub_edft_smearing_width,denskern%m(spin_channel,PUB_1K), &
                     correct_chemical_potential=.true., &
                     Ne=edft%s_occ(spin_channel), &
                     boundl=boundl, boundu=boundu, &
                     fermi_tol=loc_chemtol, entropy=entropy, &
                     entropy_approx=loc_entropy_approx,lowdin_orthog=.false.)
             else
                call apply_fermi_function(ham(spin_channel),rep%overlap,rep%inv_overlap,&
                     tmp_fermi(spin_channel), &
                     pub_edft_smearing_width,denskern%m(spin_channel,PUB_1K), &
                     correct_chemical_potential=.true., &
                     Ne=edft%s_occ(spin_channel), &
                     boundl=boundl, boundu=boundu, &
                     fermi_tol=loc_chemtol,lowdin_orthog=.false.)
             end if
          end do ! spin_channel
       end if

    else if(mode=="EIGEN") then
       do spin_channel = 1, pub_num_spins
          ! ars: diagonalise trial_ham and store e_i and M^\alpha_i
          call edft_diag_ngwf_ham(ham(spin_channel), rep%overlap, edft%num, &
               edft%mo(spin_channel), edft%h_evals(:,spin_channel))

          ! ars: round h_evals if required
          if (pub_edft_round_evals.ge.0) then
             do iorb = 1, edft%num
                call edft_round_real(edft%h_evals(iorb,spin_channel),pub_edft_round_evals)
             end do
          end if

          ! ars: check orthogonality of MOs and re-orthogonalise if necessary
          call edft_reorthogonalise_mo(edft%mo(spin_channel), denskern%m(spin_channel,PUB_1K), &
               ham(spin_channel), rep%overlap, edft%num, edft%orth_residue(spin_channel))
       end do ! pub_num_spins

       ! ars: build smeared occupancies (diagonal)
       call edft_fermi_level(edft%occ, edft%fermi, &
            edft%integrated_ne, edft%h_evals, &
            edft%s_occ(:), edft%num, edft%nbands, &
            edft%nelec, pub_edft_smearing_width, pub_edft_spin_fix)

       do spin_channel = 1, pub_num_spins
          ! ars: build density kernel
          call edft_matrix_from_eigs(denskern%m(spin_channel,PUB_1K), &
               edft%mo(spin_channel), edft%occ(:,spin_channel), edft%num)

          tmp_fermi(spin_channel)=edft%fermi(spin_channel)

!          entropy = edft_fermidirac_entropy(edft%occ, edft%num, edft%nbands, spin_channel)
       end do ! pub_num_spins

    end if


    do spin_channel = 1, pub_num_spins
       mu_lookup(mu_lookup_idx(spin_channel),mu_lookup_mu,spin_channel)=tmp_fermi(spin_channel)

       chempot(spin_channel) = tmp_fermi(spin_channel)

       loc_bounds(1)=boundl
       loc_bounds(2)=boundu

    end do ! pub_num_spins

    deallocate(tmp_fermi,stat=ierr)
    call utils_dealloc_check('edft_apply_smearing_operator','tmp_fermi',ierr)

    call timer_clock('edft_apply_smearing_operator',2)
    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving edft_apply_smearing_operator'

  end subroutine edft_apply_smearing_operator



  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  !~~~~~~~~~~~~~~~~~~~ EDFT_MODEL CREATE / DESTROY ROUTINES ~~~~~~~~~~~~~~~~~~~!
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine edft_create(edft, ham, num, nelec)

    !========================================================================!
    ! Allocates workspace for EDFT_MODEL container.                          !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                           !
    ! Modified for embedding by Joseph Prentice, May 2018
    !========================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use dense, only: dense_create
    use ngwf_representation, only: NGWF_HAM
    use rundat, only: pub_edft_extra_bands, pub_num_spins, pub_num_kpoints, &
         PUB_1K, pub_foe, pub_edft_spin_fix, pub_real_spin, pub_debug_on_root, &
         pub_ngwf_regions_ngroups, pub_edft_ham_diis_size
    use sparse, only: sparse_create
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy
    use utils, only: utils_abort, utils_alloc_check

    implicit none

    type(EDFT_MODEL), intent(inout) :: edft
    type(NGWF_HAM),   intent(in   ) :: ham
    integer,          intent(in   ) :: num
    real(kind=dp),    intent(in   ) :: nelec

    integer :: is, ierr, nsub, icntr, diis_size
    ! agrecocmplx
    logical :: loc_cmplx
    ! jcap: SPAM3_EMBED matrices for obtaining correct structure for rot
    type(SPAM3_EMBED) :: k_tmp,s_tmp

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering edft_create'

    ! agrecocmplx: consider first spin channel, since either both real
    ! or both complex
    loc_cmplx = ham%ham(1)%p%iscmplx

    ! jme: KPOINTS_DANGER
    ! This subroutine has not been checked for pub_num_kpoints > 1
    if (pub_num_kpoints /= PUB_1K) then
       call utils_abort('Subroutine edft_create not yet compatible with more&
            & than one k-point.')
    end if

    ! ars: set number of molecular orbitals
    edft%num = num
    ! jcap: set the number of regions
    nsub = pub_ngwf_regions_ngroups

    ! kkbd: Set total number of electrons
    edft%nelec = nelec

    ! kkbd: Following n_occ convention for non-spin-polarized systems.
    if (pub_num_spins == 1) then
        edft%nelec = edft%nelec/2.0_DP
    end if

    if( .not.allocated(edft%integrated_ne) ) then
       allocate(edft%integrated_ne(1:pub_num_spins), stat=ierr)
       call utils_alloc_check('edft_create', 'edft%integrated_ne', ierr)
    end if
    edft%integrated_ne(:) = 0.0_DP

    if( .not.allocated(edft%s_occ) ) then
       allocate(edft%s_occ(1:pub_num_spins), stat=ierr)
       call utils_alloc_check('edft_create', 'edft%s_occ', ierr)
    end if
    edft%s_occ(:) = 0.0_DP

    if (pub_num_spins == 2) then
       edft%integrated_ne(1) = (edft%nelec + pub_real_spin) / 2.0_dp
       edft%integrated_ne(2) = (edft%nelec - pub_real_spin) / 2.0_dp
    else
       edft%integrated_ne(:) = edft%nelec &
               &/ real(pub_num_spins,kind=dp)
    end if
    if (pub_debug_on_root) write(stdout,*) &
            "DEBUG: EDFT: Initialized spin occupancies: ", edft%integrated_ne

    edft%s_occ = edft%integrated_ne

    ! ars: set number of MOs with non-zero occupancy
    if (pub_edft_extra_bands.lt.0) then
       edft%nbands = edft%num
    else if ((int(maxval(edft%integrated_ne)) + 1 + pub_edft_extra_bands)&
            .gt. edft%num) then
       pub_edft_extra_bands = edft%num-int(maxval(edft%integrated_ne)) + 1
       if (pub_on_root) then
          write(stdout,'(a)')     "WARNING: the total number of &
               &EDFT bands exceeds the total number of NGWFs"
          write(stdout,'(a,i10)') &
               "         resetting edft_extra_bands to maximum value ", &
               pub_edft_extra_bands
       end if
       edft%nbands = edft%num
    else
       edft%nbands = int(maxval(edft%integrated_ne)) + pub_edft_extra_bands
    end if

    ! ars: initialise free energy, energy and entropy
    edft%free_energy = 0.0_DP
    edft%energy = 0.0_DP
    edft%entropy = 0.0_DP

    ! ars: allocate arrays
    if( .not.allocated(edft%occ) ) then
       allocate(edft%occ(1:edft%num,1:pub_num_spins), stat=ierr)
       call utils_alloc_check('edft_create', 'edft%occ', ierr)
    end if
    edft%occ(:,:) = 0.0_DP

    if( .not.allocated(edft%h_evals) ) then
       allocate(edft%h_evals(1:edft%num,1:pub_num_spins), stat=ierr)
       call utils_alloc_check('edft_create', 'edft%h_evals', ierr)
    end if
    edft%h_evals(:,:) = 0.0_DP

    if( .not.allocated(edft%fermi) ) then
       allocate(edft%fermi(1:pub_num_spins), stat=ierr)
       call utils_alloc_check('edft_create', 'edft%fermi', ierr)
    end if
    edft%fermi(:) = 0.0_DP

    if( .not.allocated(edft%orth_residue) ) then
       allocate(edft%orth_residue(1:pub_num_spins), stat=ierr)
       call utils_alloc_check('edft_create', 'edft%orth_residue', ierr)
    end if
    edft%orth_residue(:) = 0.0_DP


    !ja531 --> hack
    !    if(.not.pub_foe) then
    if( .not.allocated(edft%mo) ) then
       allocate(edft%mo(1:pub_num_spins), stat=ierr)
       call utils_alloc_check('edft_create', 'edft%mo', ierr)
       do is = 1, pub_num_spins
          ! agrecocmplx
          call dense_create(edft%mo(is), num, edft%num, iscmplx=loc_cmplx)
        end do
    end if
    !    else
    k_tmp%structure='K'
    s_tmp%structure='S'
    ! agrecocmplx: this should be complex for complex NGWFs (CHECK THIS!!!)
    ! jcap: need to get correct KS structure codes - this is the easiest way
    call sparse_embed_create(k_tmp, iscmplx=loc_cmplx,arrname='k_tmp')
    call sparse_embed_create(s_tmp, iscmplx=loc_cmplx,arrname='s_tmp')
!    edft%rot%structure = 'KS'
    call sparse_embed_create(edft%rot, k_tmp, s_tmp, arrname='edft%rot')
    call sparse_embed_destroy(k_tmp,arrname='k_tmp')
    call sparse_embed_destroy(s_tmp,arrname='s_tmp')

    if(.not.allocated(edft%old_ham)) then
       allocate(edft%old_ham(1:pub_num_spins), stat=ierr)
       call utils_alloc_check('edft_create', 'edft%old_ham', ierr)
       do is = 1, pub_num_spins
          call sparse_embed_create(edft%old_ham(is),ham%ham(is),arrname='edft%old_ham')
       end do
    end if

    diis_size = pub_edft_ham_diis_size

    if(.not.allocated(edft%ham_in_his)) then
       allocate(edft%ham_in_his(1:diis_size, 1:pub_num_spins), stat=ierr)
       call utils_alloc_check('edft_create', 'edft%ham_in_his', ierr)
       do is = 1, pub_num_spins
          do icntr = 1, diis_size
             ! agrecocmplx
             call sparse_embed_create(edft%ham_in_his(icntr, is), ham%ham(is),arrname='edft%ham_in_his')
          end do
       end do
    end if

    ! gab: Pulay history arrays
    if(.not.allocated(edft%ham_err_his)) then
       allocate(edft%ham_err_his(1:diis_size, 1:pub_num_spins), stat=ierr)
       call utils_alloc_check('edft_create', 'edft%ham_err_his', ierr)
       do is = 1, pub_num_spins
          do icntr = 1, diis_size
             ! agrecocmplx
             call sparse_embed_create(edft%ham_err_his(icntr, is), ham%ham(is),arrname='edft%ham_err_his')
          end do
       end do
    end if

    ! kkbd: For now I'm not able to test spin relaxing + FOE. Just exit.
    if (pub_foe .and. (pub_edft_spin_fix >= 0))then
       call utils_abort("Spin relaxation + FOE not currently supported.")
    end if

    ! ars: set flags
    edft%allocd = .true.
    edft%initialised = .false.

    if (pub_debug_on_root) then
       if(.not.pub_foe) then
          do is = 1, pub_num_spins
             write(stdout,'(a,i6)') "DEBUG: EDFT mo%nrows = ", edft%mo(is)%nrows
             write(stdout,'(a,i6)') "DEBUG: EDFT mo%mcols = ", edft%mo(is)%mcols
          end do
       end if
       write(stdout,'(a,i6)') "DEBUG: EDFT nbands = ", edft%nbands
       write(stdout,'(a,i6)') "DEBUG: EDFT num = ", edft%num
       write(stdout,'(a,i6)') "DEBUG: EDFT spin_fix = ", pub_edft_spin_fix
       write(stdout,'(a)') 'DEBUG: Leaving edft_create'
    end if

  end subroutine edft_create

  !.............................................................................

  subroutine edft_destroy(edft)

    !========================================================================!
    ! Deallocates workspace for EDFT_MODEL container.                        !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                           !
    !========================================================================!

    use constants, only: stdout
    use dense, only: dense_destroy
    use rundat, only: pub_num_spins, pub_foe, pub_edft_ham_diis_size, &
         pub_debug_on_root
    use sparse_embed, only: sparse_embed_destroy
    use utils, only: utils_dealloc_check


    implicit none

    type(EDFT_MODEL), intent(inout) :: edft

    integer :: ierr, is, icntr, diis_size

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering edft_destroy'


    !ja531-> Clean up chemical potential search stuff
    if(allocated(mu_lookup)) then
       deallocate(mu_lookup,stat=ierr)
       call utils_dealloc_check('edft_apply_smearing_operator','mu_lookup',ierr)
    end if
    if(allocated(mu_lookup_idx)) then
       deallocate(mu_lookup_idx,stat=ierr)
       call utils_dealloc_check('edft_apply_smearing_operator','mu_lookup_idx',ierr)
    end if
    if(allocated(mu_lookup_len)) then
       deallocate(mu_lookup_len,stat=ierr)
       call utils_dealloc_check('edft_apply_smearing_operator','mu_lookup_len',ierr)
    end if

    ! ars: deallocate arrays
    if( allocated(edft%occ) ) then
       deallocate(edft%occ, stat=ierr)
       call utils_dealloc_check('edft_destroy', 'edft%occ', ierr)
    end if

    if( allocated(edft%h_evals) ) then
       deallocate(edft%h_evals, stat=ierr)
       call utils_dealloc_check('edft_destroy', 'edft%h_evals', ierr)
    end if

    if( allocated(edft%fermi) ) then
       deallocate(edft%fermi, stat=ierr)
       call utils_dealloc_check('edft_destroy', 'edft%fermi', ierr)
    end if

    if( allocated(edft%integrated_ne) ) then
       deallocate(edft%integrated_ne, stat=ierr)
       call utils_dealloc_check('edft_destroy', 'edft%integrated_ne', ierr)
    end if

    if( allocated(edft%s_occ) ) then
       deallocate(edft%s_occ, stat=ierr)
       call utils_dealloc_check('edft_destroy', 'edft%s_occ', ierr)
    end if

    if( allocated(edft%orth_residue) ) then
       deallocate(edft%orth_residue, stat=ierr)
       call utils_dealloc_check('edft_destroy', 'edft%orth_residue', ierr)
    end if

    call sparse_embed_destroy(edft%rot,arrname='edft%rot')

    if(allocated(edft%old_ham)) then
       do is=1, pub_num_spins
          call sparse_embed_destroy(edft%old_ham(is),'edft%old_ham')
       end do
       deallocate(edft%old_ham, stat=ierr)
       call utils_dealloc_check('edft_destroy', 'edft%old_ham', ierr)
    end if

    if(.not. pub_foe) then
       if( allocated(edft%mo) ) then
          do is = 1, pub_num_spins
             call dense_destroy(edft%mo(is))
          end do
          deallocate(edft%mo, stat=ierr)
          call utils_dealloc_check('edft_destroy', 'edft%mo', ierr)
       end if
    end if

    ! gab: deallocate history of ham for Pulay mixing

    diis_size = pub_edft_ham_diis_size

    if (allocated(edft%ham_in_his)) then
       do is = 1, pub_num_spins
          do icntr = 1, diis_size
             call sparse_embed_destroy(edft%ham_in_his(icntr, is),arrname='edft%ham_in_his')
          end do
       end do
       deallocate(edft%ham_in_his, stat=ierr)
       call utils_dealloc_check('edft_destroy', 'edft%ham_in_his', ierr)
    end if

    if (allocated(edft%ham_err_his)) then
       do is = 1, pub_num_spins
          do icntr = 1, diis_size
             call sparse_embed_destroy(edft%ham_err_his(icntr, is),arrname='edft%ham_err_his')
          end do
       end do
       deallocate(edft%ham_err_his, stat=ierr)
       call utils_dealloc_check('edft_destroy', 'edft%ham_err_his', ierr)
    end if

    ! ars: set flags
    edft%allocd = .false.
    edft%initialised = .false.

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving edft_destroy'

  end subroutine edft_destroy

  subroutine edft_evaluate_step(t, trial_ham, initial_ham, ham, rep, edft, mdl,&
       hfxstate, denskern, hub, hub_proj_basis, lhxc_fine, ngwf_basis,         &
       paw_sphere_energies, lhxc_energy, hubbard_energy, e, s, l,     &
       foe_line_search_entropy, foe_entropy_approx, bounds, fermi_tol,&
       tmp_fermi, denswork, tmp_s, foe_mode, ientry)

    !========================================================================!
    ! Subroutine for a Hamiltonian trial step.                               !
    !------------------------------------------------------------------------!
    ! Created by Kevin Duff based on work by Alvaro Ruiz Serrano             !
    ! Adapted for embedding by Robert Charlton, 10th August 2018.            !
    !========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout
    use dense,  only: DEM, dense_convert
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_dens_dep_matrices
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: kernel_rescale_spam3
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use sparse_embed, only: sparse_embed_copy, sparse_embed_scale, &
         sparse_embed_axpy, sparse_embed_trace, SPAM3_EMBED, SPAM3_EMBED_ARRAY
    use rundat, only: pub_edft_round_evals, pub_edft_smearing_width, &
         pub_spin_fac, PUB_1K, pub_num_spins, pub_foe, pub_debug_on_root, &
         pub_edft_spin_fix, pub_edft_update_scheme
    use smearing_operator, only: apply_fermi_function
    use utils, only: utils_abort

    implicit none

    real(kind=DP),intent(in   ) :: t

    ! kkbd TODO fix the intent of these variables!!
    real(kind=DP),      intent(inout) :: e, s, l
    type(NGWF_HAM),     intent(inout) :: ham
    type(SPAM3_EMBED),  intent(inout) :: trial_ham(:)
    type(SPAM3_EMBED),  intent(inout) :: initial_ham(:)
    type(NGWF_REP),     intent(in   ) :: rep
    type(EDFT_MODEL),   intent(inout) :: edft
    type(MODEL),        intent(in   ) :: mdl
    type(HFX_STATE),    intent(inout), target :: hfxstate
    type(SPAM3_EMBED_ARRAY),  intent(inout) :: denskern
    type(HUBBARD_MODEL),intent(inout) :: hub
    type(FUNC_BASIS),   intent(in   ) :: hub_proj_basis(:), ngwf_basis(:)
    real(kind=DP),      intent(inout) :: lhxc_fine(:,:,:,:,:),&
                                         lhxc_energy,&
                                         hubbard_energy,&
                                         paw_sphere_energies(:),&
                                         tmp_fermi(:), &
                                         fermi_tol, &
                                         tmp_s, &
                                         bounds(2)
    logical,            intent(in   ) :: foe_line_search_entropy
    integer,            intent(in   ) :: foe_entropy_approx
    character(len=5),   intent(in   ) :: foe_mode
    type(DEM),          intent(inout) :: denswork(:)
    integer,            intent(in   ) :: ientry

    ! Local variables
    integer :: iorb, is
    ! kpoints issue - change this from 1 to pub_num_kpoints
    real(kind=dp),dimension(pub_num_spins, 1) :: t_occ

    if (pub_foe.and.foe_line_search_entropy) s = 0.0_dp

    select case (pub_edft_update_scheme)

    case ('PULAY_MIX')

       if (ientry.lt.1) then

          ! gab: Don't perform step in the initialisation step
          do is = 1, pub_num_spins
             ! ars: take step along the search direction in the H space
             call sparse_embed_copy(trial_ham(is),ham%ham(is))
             call sparse_embed_scale(trial_ham(is), t)
             call sparse_embed_axpy(trial_ham(is), initial_ham(is),1.0_DP-t)
          end do ! pub_num_spins

       else

          do is = 1, pub_num_spins
             call sparse_embed_copy(trial_ham(is),initial_ham(is))
          end do

          ! gab: Pulay mixing: H_in (i+1) = sum_j^i[H_in(j) t*C_j*R[H_in(j)]]
          call edft_diis_pulay_mix(trial_ham, edft%ham_in_his, &
               edft%ham_err_his, rep%inv_overlap, ientry, t, ham)

       end if

    case ('DAMP_FIXPOINT')

       do is = 1, pub_num_spins
          ! ars: take step along the search direction in the H space
          call sparse_embed_copy(trial_ham(is),ham%ham(is))
          call sparse_embed_scale(trial_ham(is), t)
          call sparse_embed_axpy(trial_ham(is), initial_ham(is),1.0_DP-t)
       end do ! pub_num_spins

    case default

       call utils_abort("EDFT mixing scheme not recogised.")

    end select

    ! ja531-> Build a new smeared density kernel based on trial_ham.
    if(.not.pub_foe) then
       ! Using Diagonalization
       call edft_apply_smearing_operator('EIGEN',  &
            rep, trial_ham, edft, tmp_fermi, denskern)
    else
       ! Or FOE
       call edft_apply_smearing_operator(foe_mode,  &
            rep, trial_ham, edft, tmp_fermi, denskern, &
            search_entropy=foe_line_search_entropy, entropy=tmp_s, ignore_history=.false., &
            entropy_approx=foe_entropy_approx, chemtol=fermi_tol, bounds=bounds, &
            densework=denswork, cleanup=.false.)
       if(foe_line_search_entropy) then
          s=s+tmp_s
       end if
    end if

    ! ars: rescale density kernel
    t_occ(:,PUB_1K) = edft%integrated_ne(:)
    if(pub_debug_on_root) write(stdout,*) &
         "DEBUG: EDFT: Rescaling denskern to ",t_occ(:,PUB_1K)
    call kernel_rescale_spam3(denskern, rep%overlap, t_occ, &
         silent=.true.)

    ! ars: Update density-dependent terms with current NGWFs and K.
    ! ars: Calculate energy.
    ! ars: Avoid updating the Hamiltonian.
    call hamiltonian_dens_dep_matrices(ham, lhxc_fine, e, &
         lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
         ngwf_basis, hub_proj_basis, hub, denskern, mdl, hfxstate, &
         .false., .false.)

    ! ars: calculate entropy from occupancies
    if(.not.pub_foe) then
       s = edft_fermidirac_entropy(edft%occ, edft%num, edft%nbands)
    else
       if(foe_line_search_entropy) then
          s = pub_spin_fac*pub_edft_smearing_width*s
       end if
    end if

    ! ars: evaluate lagrangian
    l = edft_lagrangian(edft, e, s)

  end subroutine edft_evaluate_step


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  !~~~~~~~~~~~~~~~~~~~~~ FERMI-DIRAC SMEARING AND ENTROPY ~~~~~~~~~~~~~~~~~~~~~!
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!


  subroutine edft_fermi_level_gamma(occ, fermi_energy, integrated_Ne, &
       ham_eigenvals, s_occ, num, nbands, nelec, temperature, spin_fix)

    !========================================================================!
    ! Wrapper to edft_fermi_level_kpoints in the case that the Hamiltonian   !
    ! is sampled at a single k-point only.                                   !
    !------------------------------------------------------------------------!
    ! Written by Robert Bell, November 2013.                                 !
    !========================================================================!

    use constants, only: DP
    use rundat, only: pub_num_spins

    implicit none

    ! ars: Arguments
    real(kind=dp), intent(inout) :: nelec
    integer,       intent(in   ) :: num, nbands
    real(kind=DP), intent(inout) :: occ(1:num, 1:pub_num_spins)
    real(kind=DP), intent(inout) :: fermi_energy(1:pub_num_spins)
    real(kind=DP), intent(  out) :: integrated_Ne(1:pub_num_spins)
    real(kind=DP), intent(in   ) :: ham_eigenvals(1:num, 1:pub_num_spins)
    real(kind=DP), intent(in   ) :: s_occ(1:pub_num_spins)
    real(kind=DP), intent(in   ) :: temperature
    integer,       intent(in   ) :: spin_fix

    ! local arguments
    integer, parameter :: nkpoints = 1
    real(kind=DP) :: occ_loc(1:num,nkpoints,1:pub_num_spins)
    real(kind=DP) :: ham_eigenvals_loc(1:num,nkpoints,1:pub_num_spins)
    real(kind=DP) :: weights(nkpoints)

    ! copy data
    !            num, kpoints, spins
    ham_eigenvals_loc(:,1,:) = ham_eigenvals(:,:)
    weights(1) = 1.0_DP

    call edft_fermi_level_kpoints(occ_loc, fermi_energy, integrated_Ne,   &
         ham_eigenvals_loc, s_occ, num, nbands, nkpoints, nelec, weights, &
         temperature, spin_fix)

    ! copy data
    occ(:,:) = occ_loc(:,1,:)


  end subroutine edft_fermi_level_gamma

  !.............................................................................

  subroutine edft_fermi_level_kpoints(occ, fermi_energy, integrated_Ne, &
       ham_eigenvals, s_occ, num, nbands, nkpoints, nelec, weights, temperature, &
       spin_fix)

    !========================================================================!
    ! ars: This subroutine calculates the Fermi level and the occupancies    !
    !      based on the constraint \sum_i f_i = N_e using a bisection method.!
    !      It works for any smearing scheme.                                 !
    ! rab207: modified to sample multiple k-points                           !
    ! ab: modified to set Fermi level from electrode potential in the case of!
    !     grand canonical ensemble dft                                       !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                           !
    !========================================================================!


    use comms, only: pub_on_root, pub_root_proc_id, comms_bcast
    use constants, only: DP, stdout, VERBOSE, HARTREE_IN_EVS
    use rundat, only: pub_debug_on_root, pub_output_detail, pub_num_spins, &
         pub_edft_grand_canonical, pub_edft_electrode_potential, &
         pub_edft_reference_potential, pub_charge
    use utils, only: utils_assert

    implicit none

    ! ars: Arguments
    real(kind=dp), intent(inout) :: nelec
    integer,       intent(in   ) :: num, nbands, nkpoints
    real(kind=DP), intent(inout) :: occ(1:num,1:nkpoints,1:pub_num_spins)
    real(kind=DP), intent(inout) :: fermi_energy(1:pub_num_spins)
    real(kind=DP), intent(  out) :: integrated_Ne(1:pub_num_spins)
    real(kind=DP), intent(in   ) :: ham_eigenvals(1:num,1:nkpoints,1:pub_num_spins)
    real(kind=DP), intent(in   ) :: s_occ(1:pub_num_spins)
    real(kind=DP), intent(in   ) :: weights(1:nkpoints)
    real(kind=DP), intent(in   ) :: temperature
    integer,       intent(in   ) :: spin_fix

    ! ars: Local
    integer :: iorb, iter, ik, is
    real(kind=DP) :: res, step, rel_res, sum_occ, sum_ne, delta_ne

    ! kkbd: This will indicate which way we're searching - up or down in energy.
    !       That way we can halve the search step only when we've changed direction.
    integer :: search_direction = 0
    ! ars: Parameters
    real(kind=DP), parameter :: threshold = 1.0E-12_DP
    integer, parameter :: max_iter = 1000000
    real(kind=DP), parameter :: occ_threshold = 1.0e-4
    logical, save :: said_before = .false.
    integer :: lev1, lev2
    character(len=*), parameter :: myself = 'edft_fermi_level_kpoints'

    ! --------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering edft_fermi_level_kpoints'

    ! ab: for grand canonical ensemble set fermi level from electrode potential
    ! correspondingly update the number of electrons and the quantum charge
    if (pub_edft_grand_canonical) then
       sum_ne = 0.0_DP
       do is = 1, pub_num_spins
          sum_occ = 0.0_DP
          ! ab: set Fermi level from the electrode potential
          fermi_energy(is) = pub_edft_reference_potential - &
               pub_edft_electrode_potential
          ! ab: update fractional occupancies
          do ik = 1, nkpoints
             do iorb = nbands, 1, -1
                occ(iorb,ik,is) = edft_fermidirac_distr(ham_eigenvals(iorb,ik,is),&
                     fermi_energy(is),temperature)
                sum_occ = sum_occ + weights(ik)*occ(iorb,ik,is)
             end do
          end do
          ! ab: update integrated number of electrons per spin channel
          integrated_Ne(is) = sum_occ
          sum_ne = sum_ne + integrated_Ne(is)
       end do
       ! ab: find the change in number of electrons
       delta_ne = sum_ne - nelec
       ! ab: update the quantum charge
       pub_charge = pub_charge - delta_ne
       ! ab: update the total number of electrons
       nelec = sum_ne
       ! ab: root node broadcasts
       call comms_bcast(pub_root_proc_id, occ, num*nkpoints*pub_num_spins)
       call comms_bcast(pub_root_proc_id, fermi_energy, pub_num_spins)
       call comms_bcast(pub_root_proc_id, integrated_ne, pub_num_spins)
       call comms_bcast(pub_root_proc_id, nelec)
       call comms_bcast(pub_root_proc_id, pub_charge)

       if (pub_debug_on_root) write(stdout,'(a)') &
            'DEBUG: Leaving edft_fermi_level_kpoints'
       return

    end if

    if (spin_fix == 0) then
       ! kkbd TODO this could possibly be moved into this routine
       call edft_fermi_level_kpoints_free(occ, s_occ, fermi_energy, integrated_Ne, &
            ham_eigenvals, num, nbands, nkpoints, nelec, weights, temperature)
       return
    end if

    ! ars: init trial Fermi energy
    ! rab207: changed initial energy guess to middle of band gap
    !         prevents E_f being a neat fraction of 1 Hartree
    !         when smearing temperature is small (as used in transport)
    ! kkbd: Changed to use s_occ, need to make sure it's a sensible intialization
    do is = 1, pub_num_spins
       search_direction = 0
       ! ab: check occupancy of this spin channel
       if (s_occ(is) == 0.0_DP) then ! ab: no electrons in this spin channel
          fermi_energy(is) = 0.0_DP
          integrated_ne(is) = 0.0_DP
          occ(:,:,is) = 0.0_DP
          cycle
       end if
       if (nint(s_occ(is))+1 .le. num) then
          ! jd: Doing this more defensively to ensure we don't overrun the array
          lev1 = nint(s_occ(is))+1
          lev2 = nint(s_occ(is))
          call utils_assert(lev1 <= num .and. lev1 >= 1 .and. &
               lev2 <= num .and. lev2 >=1, myself//': Level leads to &
               &out of bound access in ham_eigenvals', lev1, lev2)
          fermi_energy(is) = 0.5_DP * (minval(ham_eigenvals(lev1,:,is)) + &
               maxval(ham_eigenvals(lev2,:,is)))
       else
          fermi_energy(is) = -5.0_DP
       endif
       !fermi_energy = -5.0_DP
       step = 1.0_DP
       occ(:,:,is) = 0.0_DP


       ! ars: solve by bisection
       bisection: do iter = 1, max_iter

          ! ars: update fermi and occ
          sum_occ = 0.0_DP
          res = 0.0_DP
          do ik = 1, nkpoints
             do iorb = nbands, 1, -1
                occ(iorb,ik,is) = edft_fermidirac_distr(ham_eigenvals(iorb,ik,is),&
                     fermi_energy(is),temperature)
                sum_occ = sum_occ + weights(ik)*occ(iorb,ik,is)
             end do
          end do
          integrated_Ne(is) = sum_occ

          ! ars: calculate residue
          res = sum_occ - real(s_occ(is),kind=DP)
          rel_res = res/real(s_occ(is),kind=DP)

          ! ars: exit if converged
          if (abs(rel_res).le.threshold) then ! ars: converged ==> exit
             exit bisection
          elseif (res.gt.0.0_DP) then ! ars: we passed it - return
             fermi_energy(is) = fermi_energy(is) - step
             ! kkbd only halve search direction if we were searching up in energy last time
             if (search_direction >= 0) then
                step = 0.5_DP*step
                search_direction = -1
             end if
          elseif (res.lt.0.0_DP) then ! ars: we didn't get there, continue
             fermi_energy(is) = fermi_energy(is) + step
             ! kkbd only halve search direction if we were searching down in energy last time
             if (search_direction <= 0) then
                step = 0.5_DP*step
                search_direction = 1
             end if
          end if

       end do bisection


       ! ars: print banner
       if (pub_on_root.and.pub_output_detail.ge.VERBOSE) then
          if (abs(rel_res).ge.threshold) then
             write(stdout,'(a37,i7,a33,es10.3e2,a1)')&
                  "Fermi level failed to converge after ", iter, &
                  " iterations with relative error =", rel_res , "."
          end if
       end if


       ! ars: check that the occupancies are physical
       !      ie, the FD distr is not suddenly truncated
       if (any(occ(nbands,:,is).gt.occ_threshold).and.(.not.said_before)) then
          if (pub_on_root) then
             write(stdout,'(a)') "***** EDFT WARNING:"
             write(stdout,'(a)') "* The Fermi-Dirac distribution is truncated"
             write(stdout,'(a)') "  and the occupancies might be unphysical,"
             write(stdout,'(a)') "  leading to wrong results."
             write(stdout,'(a)') &
                  "* It is recommended to increase edft_extra_bands."
             write(stdout,'(a)') &
                  "* If the problem persists, try setting edft_extra_bands = -1"
             write(stdout,'(a)') &
                  "  to use all possible molecular orbitals."
             write(stdout,'(a)') &
                  "* If still it does not work, try adding more NGWFs per atom."
          end if
          said_before = .true.
       end if

    end do ! pub_num_spins

    ! ars: root proc broadcasts
    call comms_bcast(pub_root_proc_id, occ, num*nkpoints*pub_num_spins)
    call comms_bcast(pub_root_proc_id, fermi_energy, pub_num_spins)
    call comms_bcast(pub_root_proc_id, integrated_ne, pub_num_spins)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving edft_fermi_level_kpoints'

  end subroutine edft_fermi_level_kpoints

  subroutine edft_fermi_level_kpoints_free(occ, s_occ, fermi_energy_is, &
       integrated_Ne, ham_eigenvals, num, nbands, nkpoints, nelec, weights, &
       temperature)
    !============================================================================!
    ! This subroutine calculates the shared fermi level across spin channels for !
    ! free-spin runs. In principle it could just be combined with                !
    ! edft_fermi_level_kpoints but it is separate for development purposes for   !
    ! now.                                                                       !
    !----------------------------------------------------------------------------!
    ! Written by Kevin Duff.                                                     !
    !============================================================================!

    use comms, only: pub_on_root, pub_root_proc_id, comms_bcast
    use constants, only: DP, stdout, VERBOSE
    use rundat, only: pub_output_detail, pub_num_spins, pub_debug_on_root
    use utils, only: utils_assert

    ! kkbd: arguments
    real(kind=dp), intent(in   ) :: nelec
    integer,       intent(in   ) :: num, nbands, nkpoints

    ! kkbd: 2 fermi levels sometimes but only 1 for free spin
    real(kind=dp), intent(inout), dimension(1:pub_num_spins) :: fermi_energy_is
    real(kind=dp), intent(in   ) :: temperature

    ! kkbd: spinchannel-dep kpoint weights?
    real(kind=dp),dimension(1:nkpoints), intent(in) :: weights
    real(kind=dp),dimension(1:pub_num_spins), intent(out) :: integrated_Ne
    real(kind=dp),dimension(1:num,1:nkpoints,1:pub_num_spins),intent(in   ) :: ham_eigenvals
    real(kind=dp),dimension(1:num,1:nkpoints,1:pub_num_spins),intent(inout) :: occ
    real(kind=DP), intent(in   ) :: s_occ(1:pub_num_spins)

    !----
    ! kkbd: convergence parameters, taken from ars
    real(kind=dp), parameter :: threshold = 1.0E-12_dp
    integer,       parameter :: max_iter = 1000000
    real(kind=dp), parameter :: occ_threshold = 1.0E-4_dp
    logical,       save      :: said_before = .false.

    !----
    ! kkbd: local variables
    integer       :: iorb, iter, ik, is
    real(kind=dp) :: step, sum_occ, fermi_energy, rel_res, res

    ! kkbd: This will indicate which way we're searching - up or down in energy.
    !       That way we can halve the search step only when we've changed direction.
    integer :: search_direction = 0

    integer :: lev1, lev2
    character(len=*), parameter :: myself = 'edft_fermi_level_kpoints_free'

    ! --------------------------------------------------------------------------

    !----
    if(pub_debug_on_root) write(stdout,*) "DEBUG: Entering edft_fermi_level_kpoints_free"

    ! kkbd: Initialize fermi energy, based on ars and rab207
    ! kkbd TODO check this initialization is at all sensible??
    if (maxval((s_occ(:)))+1 .le. num) then
       ! jd: Doing this more defensively to ensure we don't overrun the array
       lev1 = nint(s_occ(1))+1
       lev2 = nint(s_occ(1))
       call utils_assert(lev1 <= num .and. lev1 >= 1 .and. &
            lev2 <= num .and. lev2 >=1, myself//': Level leads to &
            &out of bound access in ham_eigenvals', lev1, lev2)
       fermi_energy = 0.5_DP * (minval(ham_eigenvals(lev1,:,:)) + &
               maxval(ham_eigenvals(lev2,:,:)))
    else
       fermi_energy = -5.0_DP
    end if

    ! Initialize occupancies and step size for bisection
    step = 1.0_dp
    occ(:,:,:) = 0.0_dp

    ! kkbd: bisection loop, adapted from ars
    do iter = 1, max_iter
       sum_occ = 0.0_dp
       res = 0.0_dp
       integrated_ne = 0.0_dp

       do is = 1, pub_num_spins
          do ik = 1, nkpoints
             do iorb = nbands, 1, -1
                occ(iorb,ik,is) = edft_fermidirac_distr(ham_eigenvals(iorb,ik,is),&
                        fermi_energy, temperature)
                sum_occ = sum_occ + weights(ik) * occ(iorb,ik,is)
                integrated_ne(is) = integrated_ne(is) + &
                                    weights(ik) * occ(iorb,ik,is)
             end do !nbands
          end do !nkpoints
       end do !pub_num_spins

       ! kkbd: calculate residue
       res = sum_occ - nelec
       rel_res = res/nelec

       ! kkbd: exit if converged
       if (abs(rel_res).le.threshold) then
          exit
       else if (res.gt.0.0_dp) then
          fermi_energy = fermi_energy - step
          ! kkbd only halve search direction if we were searching up in energy last time
          if (search_direction >= 0) then
             step = 0.5_DP*step
             search_direction = -1
          end if
       else if (res.lt.0.0_dp) then
          fermi_energy = fermi_energy + step
          ! kkbd only halve search direction if we were searching down in energy last time
          if (search_direction <= 0) then
                step = 0.5_DP*step
                search_direction = 1
             end if
       end if

    end do !bisection

    ! print banner
    if (pub_on_root.and.pub_output_detail.ge.VERBOSE)then
       if (abs(rel_res).ge.threshold) then
          write(stdout,'(a37,i7,a33,es10.3e2,a1)')&
                "Fermi level failed to converge after ",iter, &
                " iterations with relative error = ",rel_res, "."
       end if
    end if

    ! check that the occupancies are physical
    ! i.e. the FD distribution is not suddenly truncated
    if (any(occ(nbands,:,:).gt.occ_threshold).and.(.not.said_before)) then
       if (pub_on_root) then
          write(stdout,'(a)') "***** EDFT WARNING:"
          write(stdout,'(a)') "* The Fermi-Dirac distribution is truncated"
          write(stdout,'(a)') "  and the occupancies might be unphysical,"
          write(stdout,'(a)') "  leading to wrong results."
          write(stdout,'(a)') &
               "* It is recommended to increase edft_extra_bands."
          write(stdout,'(a)') &
               "* If the problem persists, try setting edft_extra_bands = -1"
          write(stdout,'(a)') &
               "  to use all possible molecular orbitals."
          write(stdout,'(a)') &
               "* If still it does not work, try adding more NGWFs per atom."
       end if
       said_before = .true.
    end if

    ! kkbd: write this fermi energy into both spin channels
    fermi_energy_is = fermi_energy

    call comms_bcast(pub_root_proc_id, occ, num*nkpoints*pub_num_spins)
    call comms_bcast(pub_root_proc_id, fermi_energy_is, pub_num_spins)
    call comms_bcast(pub_root_proc_id, integrated_ne, pub_num_spins)

  end subroutine edft_fermi_level_kpoints_free


  !.............................................................................


  real(kind=DP) function edft_fermidirac_distr(energy,fermi,temperature)

    !========================================================================!
    ! ars: This is the Fermi-Dirac distribution of occupancies.              !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2012.                          !
    !========================================================================!

    use constants, only: DP
    use rundat, only: pub_debug_on_root

    implicit none

    real(kind=DP), intent(in) :: energy, fermi, temperature

    real(kind=DP) ::  exp_arg, exp_val

    ! ars: check for DP overflow
    exp_arg = (energy-fermi) / temperature
    if(exp_arg < 700.0_DP) then ! exp(700) is just below max_DP
       exp_val = exp(exp_arg)
       edft_fermidirac_distr = 1.0_DP / (1.0_DP + exp_val)
    else
       exp_val = huge(1.0_DP) * 0.99_DP
       edft_fermidirac_distr = 0.0_DP
    end if

    !edft_fermidirac_distr = 1.0_DP / (1.0_DP + exp_val)

  end function edft_fermidirac_distr

  !.............................................................................

  real(kind=DP) function edft_fermidirac_entropy(occ, num, nbands)

    !========================================================================!
    ! ars: This is the entropy term due to the Fermi-Dirac distribution of   !
    !      occupancies.                                                      !
    ! ars: -TS = KT* \sum_l f_l* ln(f_l) + (1-f_l)*ln(1-f_l)                 !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2012.                          !
    !========================================================================!

    use constants, only: DP, stdout
    use rundat, only: pub_edft_smearing_width, pub_debug_on_root, &
         pub_num_spins, pub_spin_fac

    implicit none

    integer, intent(in)        :: num, nbands
    real(kind=DP), intent(in)  :: occ(num, pub_num_spins)

    integer :: is, iorb

    real(kind=DP), parameter :: lthreshold = 1.0e-12_DP
    real(kind=DP), parameter :: uthreshold = 1.0_DP-lthreshold

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering edft_fermidirac_entropy'

    edft_fermidirac_entropy = 0.0_DP

    ! ars: loop over spins (Note: multiply by 2 if no spin polarisation)
    do is = 1, pub_num_spins
       do iorb = 1, nbands
          if ((occ(iorb,is) .gt. lthreshold).and.&
               (occ(iorb,is) .lt. uthreshold)) then
             edft_fermidirac_entropy = edft_fermidirac_entropy +&
                  occ(iorb,is)*log(occ(iorb,is)) +&
                  (1.0_DP-occ(iorb,is))*log(1.0_DP-occ(iorb,is))
          end if
       end do
    end do

    ! ars: fix sign and magnitude
    edft_fermidirac_entropy = &
         pub_spin_fac*pub_edft_smearing_width*edft_fermidirac_entropy

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving edft_fermidirac_entropy'

  end function edft_fermidirac_entropy

  !.............................................................................

  subroutine edft_round_real(input, n_figures)

    !=====================================================!
    ! Rounds a real number to a given number of figures.  !
    !-----------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in January 2013.     !
    !=====================================================!

    use constants, only: DP

    implicit none

    real(kind=DP), intent(inout) :: input
    integer, intent(in) :: n_figures

    real(kind=DP) :: long_input, short_input, sign, rfig


    long_input = input * 10**(real(n_figures, kind=DP))

    ! ars: get rounding figure
    rfig = mod(input * 10**(real(n_figures+1, kind=DP)),10.0_DP)
    sign = 1.0_DP; if (input.lt.0.0_DP) sign = -1.0_DP
    if (abs(rfig).lt.5.0) then
       short_input = aint(long_input)
    else
       short_input = aint(long_input)+sign
    end if
    ! ars: set input
    input = short_input * 10**(-real(n_figures, kind=DP))

  end subroutine edft_round_real

  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine edft_diag_ngwf_ham(ham, overlap, num, evecs, evals)

    !========================================================================!
    ! ars: This subroutine diagonalises one spin channel of the Hamiltonian, !
    !      ideally expressed in the NGWFs basis set.                         !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                           !
    ! Modified for embedding by Joseph Prentice, May 2018                    !
    !========================================================================!

    use constants, only: DP, stdout
    use dense, only: DEM, dense_create, dense_destroy, &
         dense_eigensolve, dense_convert
    use rundat, only: pub_eigensolver_orfac, pub_eigensolver_abstol, &
         pub_debug_on_root
    use sparse_embed, only: SPAM3_EMBED
    use timer, only: timer_clock


    implicit none

    integer,        intent(in   ) :: num
    type(SPAM3_EMBED), intent(in) :: ham
    type(SPAM3_EMBED), intent(in) :: overlap
    type(DEM),      intent(inout) :: evecs
    real(kind=DP),  intent(  out) :: evals(1:num)
    ! agrecocmplx
    logical :: loc_cmplx

    type(DEM) :: ham_buffer, overlap_buffer

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering edft_diag_ngwf_ham'

    ! ars: start timer
    call timer_clock('edft_diag_ngwf_ham',1)

    ! agrecocmplx
    loc_cmplx = overlap%p%iscmplx

    ! ars: create buffers
    ! agrecocmplx
    call dense_create(ham_buffer, num, num, iscmplx=loc_cmplx)
    call dense_create(overlap_buffer, num, num, iscmplx=loc_cmplx)

    ! ars: diagonalise Hamiltonian
    call dense_convert(ham_buffer, ham)
    call dense_convert(overlap_buffer, overlap)
    call dense_eigensolve(num, evals, ham_buffer, &
         overlap_buffer,1,evecs,pub_eigensolver_orfac, pub_eigensolver_abstol)

    ! ars: destroy buffers
    call dense_destroy(overlap_buffer)
    call dense_destroy(ham_buffer)

    ! ars: stop timer
    call timer_clock('edft_diag_ngwf_ham',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving edft_diag_ngwf_ham'

  end subroutine edft_diag_ngwf_ham


  !_____________________________________________________________________________
  !_____________________________________________________________________________


  real(kind=DP) function edft_lagrangian(edft, total_energy, entropy, &
                 & orth_residue, integrated_ne, spin_fix)

    !========================================================================!
    ! This subroutine calculates the lagrangian subject to the constraints   !
    ! of orthonormality of MO states and conservation of number of electrons.!
    ! It works in the NGWF representation.                                   !
    !------------------------------------------------------------------------!
    ! See Freysoldt, Boeck and Neugebauer, Phys. Rev. B 79 (2009).           !
    ! L = E - TS - sum_{ij} f_i * h_ji * [<\phi_i|\phi_j> - \delta_{ij}] -   !
    !                 - fermi* [\sum_l f_l - N_e]                            !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in June 2012.                           !
    ! Modified for grand canonical ensemble by Arihant Bhandari in July 2020 !
    ! by relaxing the constraint on number of electrons                      !
    !========================================================================!

    use constants, only: DP, stdout
    use rundat, only: pub_num_spins, pub_spin_fac, pub_num_kpoints, &
         PUB_1K, pub_edft_spin_fix, pub_edft_grand_canonical, pub_debug_on_root
    use utils, only: utils_abort

    implicit none

    ! ars: Arguments
    type(EDFT_MODEL), intent(in) :: edft
    real(kind=DP), optional, intent(in) :: total_energy
    real(kind=DP), optional, intent(in) :: entropy
    real(kind=DP), optional, intent(in) :: orth_residue(pub_num_spins)
    real(kind=DP), optional, intent(in) :: integrated_ne(pub_num_spins)
    integer,       optional, intent(in) :: spin_fix

    ! ars: Local variables
    real(kind=DP) :: orth_constraint, ne_constraint
    real(kind=DP) :: loc_total_energy, loc_entropy
    real(kind=DP) :: loc_orth_residue(pub_num_spins)
    real(kind=DP) :: loc_integrated_ne(pub_num_spins)
    integer :: is, loc_spin_fix

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering edft_lagrangian'

    ! jme: KPOINTS_DANGER
    ! First k-point enforced in all uses of n_occ via PUB_1K
    if (pub_num_kpoints /= PUB_1K) then
       call utils_abort('Function edft_lagrangian not yet compatible with more&
            & than one k-point.')
    end if

    orth_constraint = 0.0_DP
    ne_constraint = 0.0_DP

    ! ars: check optional variables
    if (present(total_energy)) then
       loc_total_energy = total_energy
    else
       loc_total_energy = edft%energy
    end if
    if (present(entropy)) then
       loc_entropy = entropy
    else
       loc_entropy = edft%entropy
    end if
    if (present(orth_residue)) then
       loc_orth_residue(:) = orth_residue(:)
    else
       loc_orth_residue(:) = edft%orth_residue(:)
    end if
    if (present(integrated_ne)) then
       loc_integrated_ne(:) = integrated_ne(:)
    else
       loc_integrated_ne(:) = edft%integrated_ne(:)
    end if
    if (present(spin_fix)) then
       loc_spin_fix = spin_fix
    else
       loc_spin_fix = pub_edft_spin_fix
    end if

    ! ars: loop over spins
    do is = 1, pub_num_spins

       ! ars: first the orthogonality constraint (already has spin_fac)
       orth_constraint = orth_constraint + loc_orth_residue(is)

    end do

    if (loc_spin_fix /= 0) then
       ! kkbd: note that for fixed spins we still have an integer net spin requirement...
       ! ab: check which type of ensemble
       if (pub_edft_grand_canonical) then
          do is = 1, pub_num_spins
             ne_constraint = ne_constraint + pub_spin_fac * edft%fermi(is) * &
               ( loc_integrated_ne(is) )
          end do
       else
          ! ars: and now the constraint on the conservation of number of electrons
          ! kkbd: ...which is now a real
          ! ab: while using a canonical ensemble
          do is = 1, pub_num_spins
             ne_constraint = ne_constraint + pub_spin_fac * edft%fermi(is) * &
               ( loc_integrated_ne(is) - edft%s_occ(is) )
          end do
       end if
    else ! free spin
       ! Only 1 fermi level
       if (pub_edft_grand_canonical) then
          ne_constraint = pub_spin_fac * edft%fermi(1) * &
               (sum(loc_integrated_ne(:)))
       else
       ! Electron constraint is total electrons in edft model - correct # electrons
          ne_constraint = pub_spin_fac * edft%fermi(1) * &
               (sum(loc_integrated_ne(:)) - edft%nelec)
       end if
    end if

    ! ars: lagrangian = energy + entropy - orth_constraint - ne_constraint
    edft_lagrangian = loc_total_energy + &
         loc_entropy - orth_constraint - ne_constraint

    if (pub_debug_on_root) then
       write(stdout,'(a, f22.14)') 'DEBUG: loc_total_energy = ',loc_total_energy
       write(stdout,'(a, f22.14)') 'DEBUG: loc_entropy = ', loc_entropy
       write(stdout,'(a, f22.14)') 'DEBUG: orth_constraint = ', orth_constraint
       write(stdout,'(a, f22.14)') 'DEBUG: ne_constraint = ', ne_constraint
       write(stdout,'(a)') 'DEBUG: Leaving edft_lagrangian'
    end if

  end function edft_lagrangian



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  real(kind=DP) function edft_rms_grad(initial_ham, fd_ham)

    !========================================================================!
    ! This subroutine calculates the RMS gradient of the free energy with    !
    ! respect to the occupation numbers.                                     !
    !------------------------------------------------------------------------!
    ! dL/df = h_\alpha\beta(initial) - h_\alpha\beta(FD)                     !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                           !
    ! Modified for embedding by Joseph Prentice, May 2018                    !
    !========================================================================!

    use constants, only: DP, stdout
    use rundat, only: pub_num_spins, pub_spin_fac, pub_debug_on_root
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_copy, &
         sparse_embed_axpy, sparse_embed_rms_element

    implicit none

    type(SPAM3_EMBED),   intent(in   ) :: initial_ham(1:pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: fd_ham(1:pub_num_spins)

    type(SPAM3_EMBED) :: delta_ham
    integer :: is
    real(kind=DP) :: rms

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering edft_rms_grad'

    ! ars:
    edft_rms_grad = 0.0_DP
    rms = 0.0_DP

    ! ars: create temporary
    call sparse_embed_create(delta_ham, fd_ham(1),arrname='delta_ham')

    ! ars: accumulate RMS
    do is = 1, pub_num_spins
       call sparse_embed_copy(delta_ham, fd_ham(is))
       call sparse_embed_axpy(delta_ham, initial_ham(is), -1.0_DP)
       rms = sparse_embed_rms_element(delta_ham)
       edft_rms_grad = edft_rms_grad + rms*rms
    end do

    ! ars: destroy temporaries
    call sparse_embed_destroy(delta_ham,arrname='delta_ham')

    ! ars: calculate rms_gradient
    edft_rms_grad = pub_spin_fac * &
         abs(sqrt(edft_rms_grad * 0.5_DP * pub_spin_fac))

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving edft_rms_grad'

  end function edft_rms_grad


  !_____________________________________________________________________________
  !_____________________________________________________________________________


  real(kind=DP) function edft_commutator(ham, denskern, overlap, inv_overlap)

    !========================================================================!
    ! This subroutine calculates the HKS-SKH commutator in NGWF rep.         !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                           !
    ! Modified by D.D. O'Regan for commutator invariance in July 2013.       !
    ! Modified for embedding by Joseph Prentice, May 2018                    !
    !========================================================================!

    use constants, only: DP, stdout
    use rundat, only: pub_rms_kernel_measure, pub_num_spins, pub_spin_fac, &
         pub_debug_on_root
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product, &
         sparse_embed_transpose, sparse_embed_axpy, &
         sparse_embed_rms_element, sparse_embed_trace

    implicit none

    ! ars: Arguments
    type(SPAM3_EMBED),   intent(in   ) :: ham(1:pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: denskern(1:pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: overlap
    type(SPAM3_EMBED),   intent(in   ) :: inv_overlap !ddor

    ! ars: local variables
    integer :: is
    type(SPAM3_EMBED) :: sk, skh, hks, mixed_commutator
    real(kind=DP) :: rms, trks, invariant

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering edft_commutator'

    ! ars: init
    edft_commutator= 0.0_DP

    ! ars: create matrices
    ! ddor: Changed to HKH structure
    call sparse_embed_create(sk, ham(1), denskern(1),arrname='sk')
    call sparse_embed_create(skh, sk, ham(1),arrname='skh')
    call sparse_embed_create(hks, skh, arrname='hks')
    if (.not. pub_rms_kernel_measure) then ! ddor
       call sparse_embed_create(mixed_commutator,hks,inv_overlap,arrname='mixed_commutator')
    endif

    do is = 1 , pub_num_spins

       ! ars: calculate matrices
       call sparse_embed_product(sk, overlap, denskern(is))
       call sparse_embed_product(skh, sk, ham(is))
       call sparse_embed_transpose(hks, skh)

       ! ars: calculate rms(HKS - SKH)
       call sparse_embed_axpy(hks, skh, -1.0_DP)

       ! ddor: Changed for invariant commutator
       if (pub_rms_kernel_measure) then

          rms = sparse_embed_rms_element(hks) * pub_spin_fac
          edft_commutator = edft_commutator + rms*rms

       else

          ! ddor: Calculate the mixed-index commutator C = (SKH-HKS)
          call sparse_embed_product(mixed_commutator,hks,inv_overlap)

          ! ddor: Compute occupancy of this spin channel
          call sparse_embed_trace(trks,denskern(is),overlap)

          ! ddor: Trace the squared mixed-index commutator and renormalise
          !     : to give -Tr[C^2] / N^2
          ! ddor: Allow for case of empty density kernel for one spin
          if (trks .ne. 0.0_DP) then
             call sparse_embed_trace(invariant, mixed_commutator, mixed_commutator)
             invariant = -1.0_DP * invariant / (trks)**2.0_DP
          else
             invariant = 0.0_DP
          endif

          ! ddor: Add the square of the commutator measure for this spin
          edft_commutator = edft_commutator + invariant

       endif

    end do

    ! ars: calculate commutator
    edft_commutator = abs(sqrt(edft_commutator / &
         real(pub_num_spins,kind=DP)))

    ! ars: destroy matrices
    if (.not. pub_rms_kernel_measure) then ! ddor
       call sparse_embed_destroy(mixed_commutator,arrname='mixed_commutator')
    endif
    call sparse_embed_destroy(hks,arrname='hks')
    call sparse_embed_destroy(skh,arrname='skh')
    call sparse_embed_destroy(sk,arrname='sk')

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving edft_commutator'

  end function edft_commutator


  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine edft_print_iter(iter, edft, energy, entropy, &
       lagrangian, delta_lagrangian, &
       rms_gradient, commutator, delta_conv, step)

    !========================================================================!
    ! This subroutine prints the result of an iteration on the screen.       !
    ! It is adapted to display values related to finite-temperature DFT      !
    ! calculations.                                                          !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                           !
    ! Modified for grand canonical ensemble by Arihant Bhandari in July 2020 !
    !========================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, VERBOSE, stdout
    use rundat, only: pub_output_detail, pub_edft_write_occ, pub_num_spins, &
         pub_spin_fac, pub_num_kpoints, PUB_1K, pub_foe, &
         pub_edft_rms_gradient_thres, pub_edft_spin_fix, &
         pub_edft_grand_canonical, pub_charge
    use utils, only: utils_abort, utils_unit

    implicit none

    ! ars: Arguments
    integer,          intent(in) :: iter
    type(EDFT_MODEL), intent(in) :: edft
    real(kind=DP),    intent(in) :: energy
    real(kind=DP),    intent(in) :: entropy
    real(kind=DP),    intent(in) :: lagrangian, delta_lagrangian
    real(kind=DP),    intent(in) :: commutator
    real(kind=DP),    intent(in) :: rms_gradient
    real(kind=DP),    intent(in) :: delta_conv(1:pub_num_spins)
    real(kind=DP),    intent(in) :: step


    integer :: is, occ_unit
    character(len=80) :: fmt
    character(len=2) :: symbol
    real(kind=DP) :: ne_constraint, orth_constraint

    ! jme: KPOINTS_DANGER
    ! First k-point enforced in all uses of n_occ via PUB_1K
    if (pub_num_kpoints /= PUB_1K) then
       call utils_abort('Subroutine edft_print_iter not yet compatible with more&
            & than one k-point.')
    end if

    ! ars: calculate residual contributions to the Lagrangian, due to
    !      1) Orthogonality contraint
    !      2) Conservation of the number of electrons in canonical ensemble
    ! kkbd: updated ne_constraint for free-spin runs. Note that we still have
    !       an integer net spin constraint for fixed spin runs...

    ne_constraint = 0.0_DP

    if (pub_edft_spin_fix /= 0) then
       if (pub_edft_grand_canonical) then
          do is = 1, pub_num_spins
             ne_constraint = ne_constraint + &
                  edft%fermi(is) * (edft%integrated_ne(is))
          end do
       else
          do is = 1, pub_num_spins
             ne_constraint = ne_constraint - &
                  edft%fermi(is) * (edft%integrated_ne(is) - edft%s_occ(is))
          end do
       end if
    else ! free spin
       if (pub_edft_grand_canonical) then
          ne_constraint = edft%fermi(1) * &
               (sum(edft%integrated_ne(:)))
       else
          ne_constraint = edft%fermi(1) * &
               (sum(edft%integrated_ne(:)) - edft%nelec)
       end if
    end if ! spin_fix

    ne_constraint = pub_spin_fac * ne_constraint
    orth_constraint = -sum(edft%orth_residue)

    if (pub_on_root) then

       ! ars: if detail=NORMAL -> print one line only
       if ((iter.eq.0).or.(pub_output_detail.ge.VERBOSE)) then
          write(stdout,'(/,a80)') "-----------------------------------&
               &---------------------------------------------"
          if (pub_edft_grand_canonical) then
             if(pub_edft_rms_gradient_thres.gt.0.0_dp) then
                write(stdout,'(a80)')&
                     "  Iter        RMS_gradient          Commutator&
                     &    Grand pot.(L=E-TS-muN)  DeltaL"
             else
                write(stdout,'(a80)')&
                     "  Iter                              Commutator&
                     &    Grand pot.(L=E-TS-muN)  DeltaL"
             end if
          else
             if(pub_edft_rms_gradient_thres.gt.0.0_dp) then
                write(stdout,'(a80)')&
                     "  Iter        RMS_gradient          Commutator&
                     &    Free_energy(A=E-TS)     DeltaA"
             else
                write(stdout,'(a80)')&
                     "  Iter                              Commutator&
                     &    Free_energy(A=E-TS)     DeltaA"
             end if
          end if
       end if

       ! ars: write line with the correct format
       if(pub_edft_rms_gradient_thres.gt.0.0_dp) then
          write(fmt,'(a)') '(a2, i4, 2(1x, f19.12)'
       else
          write(fmt,'(a)') '(a2, i4, (20x), (1x, f19.12)'
       end if
       if (abs(lagrangian)<100000.0_DP) then
          write(fmt,'(a,a)') trim(fmt), ', 1x, f22.14, 1x, es10.2e2)'
       else
          write(fmt,'(a,a)') trim(fmt), ', 1x, f22.12, 1x, es10.2e2)'
       end if

       symbol = "# "
       if(pub_edft_rms_gradient_thres.gt.0.0_dp) then
          write(stdout,fmt) symbol, iter, rms_gradient, commutator, &
               lagrangian, delta_lagrangian
       else
          write(stdout,fmt) symbol, iter, commutator, &
               lagrangian, delta_lagrangian
       end if


       ! ars: write line with the correct format
       write(fmt,'(a)') '(a,f20.12, a5'
       if (abs(lagrangian)<100000.0_DP) then
          write(fmt,'(a,a)') trim(fmt), ', f22.14)'
       else
          write(fmt,'(a,a)') trim(fmt), ', f22.12)'
       end if

       ! ars: if detail=VERBOSE -> print extra info
       if (pub_output_detail.ge.VERBOSE) then

          ! ars: print spin-independent info
          write(stdout,'(a)') ""
          write(stdout,fmt) &
               "Step                       = ", step
          write(stdout,fmt)   &
               "Energy (E)                 = ", energy
          write(stdout,fmt)   &
               "Entropy (-TS)              = ", entropy

          if (pub_edft_grand_canonical) then
             write(stdout,fmt)   &
                  "Chemical potential (-muN)  = ", -ne_constraint
             write(stdout,fmt)   &
                  "Grand potential(L=E-TS-muN)= ", lagrangian
             write(stdout,fmt)   &
                  "Est. 0K Energy 0.5*(E+L)   = ", 0.5_DP*(lagrangian+energy)
             write(stdout,fmt)   &
                  "Charge on quantum system   = ", pub_charge
          else
             write(stdout,fmt)   &
                  "Free energy (A=E-TS)       = ", lagrangian
             write(stdout,fmt)   &
                  "Est. 0K Energy 0.5*(E+A)   = ", 0.5_DP*(lagrangian+energy)
          end if

          write(stdout,fmt)   &
               "Residual Non-orthogonality = ", orth_constraint

          if (.not.pub_edft_grand_canonical) then
             write(stdout,fmt)   &
                  "Residual N_electrons       = ", ne_constraint
          end if

          write(stdout,'(a)') ""

          ! ars: write spin-dependent info
          if (pub_edft_grand_canonical) then
             write(stdout,'(a72)') &
                  "Spin          Integrated_Ne&
                  &            Fermi_level   Delta_Integrated_Ne"
          else
             write(stdout,'(a72)') &
                  "Spin          Integrated_Ne&
                  &            Fermi_level     Delta_Fermi_level"
          end if

          do is = 1, pub_num_spins
             write(stdout,'(i4, 2(1x, f22.12), 12x, es10.2e2)') &
                  is, pub_spin_fac * edft%integrated_ne(is), &
                  edft%fermi(is), delta_conv(is)
          end do

          write(stdout,'(a)') ""

       end if

       ! ars: write occupancies and eigenvalues of the Hamiltonian
       if ((.not.pub_foe).and.(pub_output_detail.ge.VERBOSE)) call edft_print_evals(edft)

       ! ars: print file with the occupancies if required
       if ((.not.pub_foe).and.pub_edft_write_occ) then
          occ_unit = utils_unit()
          call edft_print_evals(edft, iter=iter, print_all = .true., &
               newfile=.false., outunit=occ_unit)
       end if


    end if

  end subroutine edft_print_iter


  !_____________________________________________________________________________
  !_____________________________________________________________________________


  logical function edft_check_convergence(delta_energy, delta_entropy, &
       delta_lagrangian, delta_conv6, rms_gradient, commutator, nat)

    !========================================================================!
    ! This subroutine checks convergence of SCF free energy calculations.    !
    ! It assumes that the occupation numbers are smeared.                    !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                           !
    ! Modified by Gabriel Bramley - per atom energy tolerances, Nov 2020.    !
    !========================================================================!

    use constants, only: DP, stdout, VERBOSE, very_big_double
    use comms, only: pub_on_root
    use rundat, only: pub_output_detail, pub_edft_grand_canonical, &
         pub_edft_free_energy_thres, pub_edft_energy_thres, &
         pub_edft_entropy_thres, pub_edft_rms_gradient_thres, &
         pub_edft_commutator_thres, pub_edft_fermi_thres, pub_num_spins, &
         pub_edft_nelec_thres

    implicit none

    real(kind=DP), intent(in) :: delta_energy, delta_entropy, delta_lagrangian
    real(kind=DP), intent(in) :: delta_conv6(1:pub_num_spins)
    real(kind=DP), intent(in) :: rms_gradient, commutator
    integer, intent(in) :: nat

    integer :: is, per
    integer, parameter :: ncriteria = 6
    logical :: conv_criteria(1:ncriteria)
    real(kind=DP) :: thres

    character(len=4), parameter :: yes=" Yes"
    character(len=4), parameter :: no= "  No"
    character(len=4), parameter :: na= "n.a."

    character(len=4) :: conv_text(1:ncriteria)

    edft_check_convergence = .false.
    conv_criteria(1:ncriteria) = .false.
    conv_text(1:ncriteria)=no

    ! ars: validate convergence criteria:

    ! change in the Lagrangian
    if (pub_edft_free_energy_thres.gt.0.0_DP) then
       if(abs(delta_lagrangian/nat).lt.pub_edft_free_energy_thres) then
          conv_criteria(1) = .true.
          conv_text(1)=yes
       end if
    else
       conv_criteria(1) = .true.
       conv_text(1)=na
    end if

    ! change in the energy
    if (pub_edft_energy_thres.gt.0.0_DP) then
       if(abs(delta_energy/nat).lt.pub_edft_energy_thres) then
          conv_criteria(2) = .true.
          conv_text(2)=yes
       end if
    else
       conv_criteria(2) = .true.
       conv_text(2)=na
    end if

    ! change in the entropy
    if (pub_edft_entropy_thres.gt.0.0_DP) then
       if(abs(delta_entropy/nat).lt.pub_edft_entropy_thres) then
          conv_criteria(3) = .true.
          conv_text(3)=yes
       end if
    else
       conv_criteria(3) = .true.
       conv_text(3)=na
    end if

    ! change in the RMS gradient
    if (pub_edft_rms_gradient_thres.gt.0.0_DP) then
       if(abs(rms_gradient).lt.pub_edft_rms_gradient_thres) then
          conv_criteria(4) = .true.
          conv_text(4)=yes
       end if
    else
       conv_criteria(4) = .true.
       conv_text(4)=na
    end if

    ! change in the RMS commutator hf-fh
    if (pub_edft_commutator_thres.gt.0.0_DP) then
       if(abs(commutator).lt.pub_edft_commutator_thres) then
          conv_criteria(5) = .true.
          conv_text(5)=yes
       end if
    else
       conv_criteria(5) = .true.
       conv_text(5)=na
    end if

    ! ab: final convergence criteria conv6
    ! convergence of the number of electrons in grand canonical ensemble
    ! or the convergence of the fermi_level in canonical ensemble
    if (pub_edft_grand_canonical) then
       thres = pub_edft_nelec_thres
       per = nat
    else
       thres = pub_edft_fermi_thres
       per = 1
    end if

    if (thres.gt.0.0_DP) then
       fermi_loop: do is = 1, pub_num_spins
          if (abs(delta_conv6(is)/per).lt.thres) then
             conv_criteria(6) = .true.
             conv_text(6)=yes
          else
             conv_criteria(6) = .false.
             conv_text(6)=no
             exit fermi_loop
          end if
       end do fermi_loop
    else
       conv_criteria(6) = .true.
       conv_text(6)=na
    end if

    if ((pub_on_root).and.(pub_output_detail.ge.VERBOSE)) then
       write(stdout,'(/,5x,a)') &
            "                       |  Abs_value |  Threshold | Converged"
       write(stdout,'(5x,a,2(es10.2e2," | "),5x,a4)') &
            "Delta_Free_energy/atom | ", abs(delta_lagrangian/nat), &
            pub_edft_free_energy_thres, conv_text(1)
       write(stdout,'(5x,a,2(es10.2e2," | "),5x,a4)') &
            "Delta_Energy/atom      | ", abs(delta_energy/nat), &
            pub_edft_energy_thres, conv_text(2)
       if(delta_entropy /= very_big_double) then
          write(stdout,'(5x,a,2(es10.2e2," | "),5x,a4)') &
               "Delta_Entropy/atom     | ", abs(delta_entropy/nat), &
               pub_edft_entropy_thres, conv_text(3)
       else
          write(stdout,'(5x,a,es10.2e2,a3,5x,a4)') &
               "Delta_Entropy/atom     |     n.a.   | ", &
               pub_edft_fermi_thres, " | ", na
       end if

       if (pub_edft_rms_gradient_thres.gt.0.0_DP) then
          write(stdout,'(5x,a,2(es10.2e2," | "),5x,a4)') &
               "RMS_gradient           | ", rms_gradient, &
               pub_edft_rms_gradient_thres, conv_text(4)
       end if

       write(stdout,'(5x,a,2(es10.2e2," | "),5x,a4)') &
            "Commutator             | ", commutator, &
            pub_edft_commutator_thres, conv_text(5)

       if (pub_edft_grand_canonical) then
          write(stdout,'(5x,a,2(es10.2e2," | "),5x,a4)') &
               "Delta_Int_Ne/atom   | ", maxval(abs(delta_conv6/nat)), &
               pub_edft_nelec_thres, conv_text(6)
       else
          if(maxval(abs(delta_conv6)) /= very_big_double) then
             write(stdout,'(5x,a,2(es10.2e2," | "),5x,a4)') &
                  "Delta_Fermi_level      | ", maxval(abs(delta_conv6)), &
                  pub_edft_fermi_thres, conv_text(6)
          else
             write(stdout,'(5x,a,es10.2e2,a3,5x,a4)') &
                  "Delta_Fermi_level      |     n.a.   | ", &
                  pub_edft_fermi_thres, " | ", na
          end if
       end if
    end if

    edft_check_convergence = conv_criteria(1).and.conv_criteria(2).and.&
         conv_criteria(3).and.conv_criteria(4).and.conv_criteria(5).and.&
         conv_criteria(6)


  end function edft_check_convergence



  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine edft_polynomial_fit(success,poly_coeffs,x_axis,y_axis,poly_order)

    !========================================================================!
    ! Performs a polynomial fit of order 'poly_order' given (poly_order+1)   !
    ! fitting points in X-Y coordinates.                                     !
    ! P(x) = \sum_0^(p+1) k_p x^p                                            !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in October 2012.                        !
    !========================================================================!

    use constants, only: DP, stdout
    use rundat, only: pub_debug_on_root
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    integer,       intent(in   ) :: poly_order
    logical,       intent(  out) :: success
    real(kind=DP), intent(inout) :: poly_coeffs(poly_order+1)
    real(kind=DP), intent(in   ) :: x_axis(poly_order+1)
    real(kind=DP), intent(in   ) :: y_axis(poly_order+1)

    ! LAPACK subroutine
    external :: dgesv

    ! ars: local
    real(kind=DP), allocatable :: matrix(:,:), vector(:)
    integer, allocatable :: pivots(:)
    integer :: col, row, ierr, info

    if (pub_debug_on_root) write(stdout,*) &
         "DEBUG: entering edft_polynomial_fit"

    ! ars: allocate
    allocate(matrix(poly_order+1, poly_order+1), stat=ierr)
    call utils_alloc_check('edft_polynomial_fit','matrix',ierr)
    allocate(vector(poly_order+1), stat=ierr)
    call utils_alloc_check('edft_polynomial_fit','vector',ierr)
    allocate(pivots(poly_order+1), stat=ierr)
    call utils_alloc_check('edft_polynomial_fit','pivots',ierr)

    ! ars: init
    success = .false.
    poly_coeffs = 0.0_DP

    ! ars: construct Vandervonde matrix and vector
    do col= 1, poly_order+1
       do row = 1, poly_order+1
          matrix(row,col) = x_axis(row)**(col-1)
       end do
       vector(col) = y_axis(col)
    end do

    if (pub_debug_on_root) then
       write(stdout,*) "DEBUG: entering edft_polynomial_fit"
       write(stdout,*) "matrix = "
       do row = 1, poly_order+1
          do col= 1, poly_order+1
             write(stdout,'(e16.8)',advance='no') matrix(row,col)
          end do
          write(stdout,*)
       end do
       write(stdout,*) "vector = "
       do row = 1, poly_order+1
          write(stdout,*) vector(row)
       end do
       write(stdout,*) "x_axis = "
       do row = 1, poly_order+1
          write(stdout,*) x_axis(row)
       end do
       write(stdout,*) "y_axis = "
       do row = 1, poly_order+1
          write(stdout,*) y_axis(row)
       end do
    end if

    ! ars: solve with LAPACK
    call dgesv(poly_order+1,1,matrix,poly_order+1,pivots,vector,&
         poly_order+1,info)

    if (info.eq.0) then
       poly_coeffs(:) = vector(:)
       success = .true.
    end if

    ! ars: deallocate
    deallocate(pivots, stat=ierr)
    call utils_dealloc_check('edft_polynomial_fit','pivots',ierr)
    deallocate(vector, stat=ierr)
    call utils_dealloc_check('edft_polynomial_fit','vector',ierr)
    deallocate(matrix, stat=ierr)
    call utils_dealloc_check('edft_polynomial_fit','matrix',ierr)

    if (pub_debug_on_root) then
       write(stdout,*) "poly_order = ", poly_order
       write(stdout,*) "success = ", success
       write(stdout,*) "DEBUG: leaving edft_polynomial_fit"
    end if

  end subroutine edft_polynomial_fit

  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine edft_fit_second_order_poly_3points(success, step, &
       predicted_functional, x1, x2, x3, y1, y2, y3)

    !==========================================================================!
    ! Fit a second-order polynomial given three trial points in X-Y coordinates!
    ! ars: P(x) = a + bx + cx^2                                                !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in October 2012.                          !
    !==========================================================================!

    use constants, only: DP, stdout
    use rundat, only: pub_debug_on_root

    implicit none

    ! ars: arguments
    logical,       intent(  out) :: success
    real(kind=DP), intent(  out) :: step, predicted_functional
    real(kind=DP), intent(in   ) :: x1, x2, x3, y1, y2, y3

    ! ars: local
    integer, parameter :: poly_order = 2
    real(kind=DP), dimension(poly_order+1) :: x_axis, y_axis, poly_coeffs
    real(kind=DP) :: aa, bb, cc
    logical :: step_is_a_min

    if (pub_debug_on_root) write(stdout,*) &
         "DEBUG: entering edft_fit_second_order_poly_3points"

    ! ars: init
    success = .false.; step_is_a_min = .false.
    aa=0.0_DP; bb=0.0_DP; cc=0.0_DP
    step=1.0e4_DP; predicted_functional=1.0e6_DP

    ! ars: fit second-order polynomial to data
    x_axis = (/ x1, x2, x3 /); y_axis = (/ y1, y2, y3 /)
    call edft_polynomial_fit(success, poly_coeffs, x_axis, y_axis, poly_order)

    ! ars: find minimum and predicted functional
    if (success) then
       aa = poly_coeffs(1)
       bb = poly_coeffs(2)
       cc = poly_coeffs(3)

       ! ars: check that this is a minimum (second derivative criterion)
       if (cc.gt.0.0_DP) then
          ! ars: it's a minimum :)
          step_is_a_min = .true.
          success = .true.
          step = -0.5_DP*bb/cc
          predicted_functional = aa - 0.25_DP*bb*bb/cc
       else
          ! ars: it's a maximum :(
          step_is_a_min = .false.
          success = .false.
          step = 1.0e6_DP
          predicted_functional = 1.0e6_DP
       end if
    else
       step = 1.0e6_DP
       predicted_functional = 1.0e6_DP
    end if

    if (pub_debug_on_root) then
       write(stdout,*) "aa = ", aa
       write(stdout,*) "bb = ", bb
       write(stdout,*) "cc = ", cc
       write(stdout,*) "step = ", step
       write(stdout,*) &
            "predicted_functional = ", predicted_functional
       write(stdout,*) "step_is_a_min = ", step_is_a_min
       write(stdout,*) "success = ", success
       write(stdout,*) &
            "DEBUG: leaving edft_fit_second_order_poly_3points"
    end if

  end subroutine edft_fit_second_order_poly_3points

  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine edft_fit_third_order_poly_4points(success, step, &
       predicted_functional, x1, x2, x3, x4, y1, y2, y3, y4)

    !==========================================================================!
    ! Fit a third-order polynomial given three trial points in X-Y coordinates !
    ! ars: P(x) = a + bx + cx^2 + dx^3                                         !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in October 2012.                          !
    !==========================================================================!

    ! ars:
    use constants, only: DP, stdout
    use rundat, only: pub_debug_on_root

    implicit none

    ! ars: arguments
    logical,       intent(  out) :: success
    real(kind=DP), intent(  out) :: step, predicted_functional
    real(kind=DP), intent(in   ) :: x1, x2, x3, x4, y1, y2, y3, y4

    ! ars: local
    integer, parameter :: poly_order = 3
    real(kind=DP), dimension(poly_order+1) :: x_axis, y_axis, poly_coeffs
    real(kind=DP) :: aa, bb, cc, dd
    real(kind=DP) :: qq, sqrt_arg, root1, root2, func1, func2, cc_sign
    logical :: step_is_a_min

    if (pub_debug_on_root) write(stdout,*) &
         "DEBUG: entering edft_fit_third_order_poly_4points"

    ! ars: init
    success = .false.; step_is_a_min = .false.
    aa=0.0_DP; bb=0.0_DP; cc=0.0_DP; dd=0.0_DP; sqrt_arg=0.0_DP; qq=0.0_DP
    cc_sign=0.0_DP; root1=0.0_DP; root2=0.0_DP; func1=0.0_DP; func2=0.0_DP

    ! ars: fit third-order polynomial to data
    x_axis = (/ x1, x2, x3, x4 /); y_axis = (/ y1, y2, y3, y4 /)
    call edft_polynomial_fit(success, poly_coeffs, x_axis, y_axis, poly_order)

    ! ars: calculate step and predicted functional
    if (success) then

       ! ars: set factors
       aa = poly_coeffs(1); bb = poly_coeffs(2)
       cc = poly_coeffs(3); dd = poly_coeffs(4)

       ! ars: determine whether the sqrt has a positive argument
       sqrt_arg = 4.0_DP*cc*cc - 12.0_DP*bb*dd

       if (sqrt_arg.gt.0.0_DP) then

          ! ars: calculate the two candidates for a minimum
          cc_sign = 1.0_DP
          if (cc.lt.0.0_DP) cc_sign = -1.0_DP
          qq = -0.5_DP*( 2.0_DP*cc + cc_sign *sqrt(sqrt_arg))
          root1 = qq/(3.0_DP*dd); root2 = bb/qq

          ! ars: evaluate root1
          func1 = dd*root1**3 + cc*root1**2 + bb*root1 + aa
          func2 = dd*root2**3 + cc*root2**2 + bb*root2 + aa

          ! ars: select minimum
          if (func1.lt.func2) then
             predicted_functional = func1
             step = root1
          else
             predicted_functional = func2
             step = root2
          end if

          ! ars: check that the chosen root is an actual minimum
          !      (second derivative criterion)
          if ((2.0_DP*cc + 6.0_DP*dd*step).gt.0.0_DP) then
             ! ars: it's a minimum :)
             success = .true.
             step_is_a_min = .true.
          else
             ! ars: it's a maximum :(
             success = .false.
             step_is_a_min = .false.
             step = 1.0e6_DP
             predicted_functional = 1.0e6_DP
          end if

       else
          ! ars: the argument for the sqrt is negative - complex roots - reject
          success = .false.
       end if

    else
       ! ars: fitting failed
       step = 1.0e6_DP
       predicted_functional = 1.0e6_DP
    end if

    ! ars: re-check that the fitting was successful
    if (.not.success) then
       step = 1.0e6_DP
       predicted_functional = 1.0e6_DP
    end if

    if(pub_debug_on_root) then
       write(stdout,*) "aa = ", aa
       write(stdout,*) "bb = ", bb
       write(stdout,*) "cc = ", cc
       write(stdout,*) "dd = ", dd
       write(stdout,*) "sqrt_arg = ", sqrt_arg
       write(stdout,*) "cc_sign = ", cc_sign
       write(stdout,*) "qq = ", qq
       write(stdout,*) "root1 = ", root1
       write(stdout,*) "func1 = ", func1
       write(stdout,*) "root2 = ", root2
       write(stdout,*) "func2 = ", func2
       write(stdout,*) "step = ", step
       write(stdout,*) &
            "predicted_functional = ", predicted_functional
       write(stdout,*) "step_is_a_min = ", step_is_a_min
       write(stdout,*) "success = ", success
       write(stdout,*) &
            "DEBUG: leaving edft_fit_third_order_poly_4points"
    end if

  end subroutine edft_fit_third_order_poly_4points


  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine edft_fit_fourth_order_poly_5points(success, step,&
       predicted_functional, x1, x2, x3, x4, x5, y1, y2, y3, y4, y5)

    !==========================================================================!
    ! Fit a fourth-order polynomial given three trial points in X-Y coordinates!
    ! ars: P(x) = a + bx + cx^2 + dx^3 + ex^4                                  !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in October 2012.                          !
    !==========================================================================!

    use constants, only: DP, stdout, TWO_PI, pi
    use rundat, only: pub_debug_on_root

    implicit none

    ! ars: arguments
    logical,       intent(  out) :: success
    real(kind=DP), intent(  out) :: step, predicted_functional
    real(kind=DP), intent(in   ) :: x1, x2, x3, x4, x5, y1, y2, y3, y4, y5

    ! ars: local
    integer, parameter :: poly_order = 4
    real(kind=DP), dimension(poly_order+1) :: x_axis, y_axis, poly_coeffs
    real(kind=DP) :: aa, bb, cc, dd, ee
    real(kind=DP) :: alpha, beta, gamma
    real(kind=DP) :: qq, qq_cubed, rr, rr_squared, cap_A, cap_B, theta, sgn_rr
    real(kind=DP) :: root1, root2, root3, func1, func2, func3
    logical :: all_roots_are_real, step_is_a_min

    ! ja531 --> tmp buffer for principal cube root
    complex(kind=dp) :: pcr

    if (pub_debug_on_root) write(stdout,*) &
         "DEBUG: entering edft_fit_fourth_order_poly_5points"

    ! ars: init
    success = .false.; step_is_a_min = .false.; all_roots_are_real = .false.
    step = 1.0e4_DP; predicted_functional = 1.0e6_DP
    aa=0.0_DP; bb=0.0_DP; cc=0.0_DP; dd=0.0_DP; ee=0.0_DP
    qq=0.0_DP; rr=0.0_DP; cap_A=0.0_DP; cap_B=0.0_DP; theta = 0.0_DP
    qq_cubed = 0.0_DP; rr_squared=0.0_DP; sgn_rr = 1.0_DP
    root1=1.0e4_DP; root2=1.0e4_DP; root3=1.0e4_DP
    func1=1.0e6_DP; func2=1.0e6_DP; func3=1.0e6_DP

    ! ars: fit fourth-order polynomial to data
    x_axis = (/ x1, x2, x3, x4, x5 /); y_axis = (/ y1, y2, y3, y4, y5 /)
    call edft_polynomial_fit(success, poly_coeffs, x_axis, y_axis, poly_order)

    ! ars: calculate step and predicted functional
    if (success) then

       ! ars: set factors
       aa = poly_coeffs(1); bb = poly_coeffs(2)
       cc = poly_coeffs(3); dd = poly_coeffs(4); ee = poly_coeffs(5)

       !if(pub_on_root) write(*,*) "ja531 --> coeffs : ", aa,bb,cc,dd,ee

       ! ars: calculate alpha, beta, gamma
       alpha = 0.75_DP*dd/ee ; beta = 0.5_DP*cc/ee ; gamma = 0.25_DP*bb/ee

       ! ars: calculate Q and R coefficients
       qq = (alpha*alpha -3.0_DP*beta)/9.0_DP
       qq_cubed = qq*qq*qq
       rr = (2.0_DP*alpha*alpha*alpha-9.0_DP*alpha*beta+27.0_DP*gamma)/54.0_DP
       rr_squared = rr*rr

       ! ars: determine if all the roots of the polynomial are real
       if (rr_squared.lt.qq_cubed) then
          all_roots_are_real = .true.
       else
          all_roots_are_real = .false.
       end if

       ! ars: proceed different if there are complex roots
       if (all_roots_are_real) then

          if (qq_cubed.gt.0.0_DP) then
             !if(pub_on_root) write(*,*) "ja531 --> qq_cubed : ", qq_cubed
             ! ars: calculate theta
             theta = acos(rr/sqrt(qq_cubed))

             ! ars: calculate roots
             root1 = -2.0_DP*sqrt(qq) *cos(theta/3.0_DP) - alpha/3.0_DP
             root2 = -2.0_DP*sqrt(qq) *cos((theta+TWO_PI)/3.0_DP) - alpha/3.0_DP
             root3 = -2.0_DP*sqrt(qq) *cos((theta-TWO_PI)/3.0_DP) - alpha/3.0_DP

             ! ars: evaluate lagrangian
             func1 = aa + bb*root1 + cc*root1*root1 + dd*root1*root1*root1 +&
                  ee*root1*root1*root1*root1
             func2 = aa + bb*root2 + cc*root2*root2 + dd*root2*root2*root2 +&
                  ee*root2*root2*root2*root2
             func3 = aa + bb*root3 + cc*root3*root3 + dd*root3*root3*root3 +&
                  ee*root3*root3*root3*root3

             ! ars: check which root gives the minimum
             if ( (func1.lt.func2).and.(func1.lt.func3) ) then
                step = root1; predicted_functional = func1
             else if ( (func2.lt.func1).and.(func2.lt.func3) ) then
                step = root2; predicted_functional = func2
             else if ( (func3.lt.func1).and.(func3.lt.func2) ) then
                step = root3; predicted_functional = func3
             end if

             ! ars: check that the chosen root is an actual minimum
             !      (second derivative criterion)
             if ((2.0_DP*cc+6.0_DP*dd*step+12.0_DP*ee*step*step).gt.0.0_DP) then
                ! ars: it's a minimum :)
                success = .true.
                step_is_a_min = .true.
             else
                ! ars: it's a maximum :(
                step_is_a_min = .false.
                success = .false.
                step = 1.0e4_DP; predicted_functional = 1.0e6_DP
             end if

          else
             success = .false.
          end if

       else

          ! ars: determine sign of R
          if (rr.gt.0.0_DP) then
             sgn_rr = 1.0_DP
          else
             sgn_rr = -1.0_DP
          end if

          ! ars: determine cap_A and cap_B

          ! ja531 --> convert principal cube root to real cube root... otherwise NaN for -ve arg!
          pcr=cmplx(rr + sqrt(rr_squared - qq_cubed),0.0_dp,dp)
          pcr=pcr**(1.0_DP/3.0_DP)
          if(abs(aimag(pcr)+aimag(pcr*exp(cmplx(0.0_dp,(2.0_dp/3.0_dp)*pi,dp))))<epsilon(1.0_dp)) then
             pcr=pcr*exp(cmplx(0.0_dp,-(2.0_dp/3.0_dp)*pi,dp))
          else if(abs(aimag(pcr)+aimag(pcr*exp(cmplx(0.0_dp,-(2.0_dp/3.0_dp)*pi,dp))))<epsilon(1.0_dp)) then
             pcr=pcr*exp(cmplx(0.0_dp,(2.0_dp/3.0_dp)*pi,dp))
          end if

          cap_A = -1.0_DP * real(pcr,dp)


          if (cap_A.lt.tiny(1.0_DP)) then
             cap_B = 0.0_DP
          else
             cap_B = qq/cap_A
          end if

          ! ars: determine the real root, and set absurd values for the others
          root1 = cap_A + cap_B - alpha/3.0_DP
          root2 = 1.0e4_DP
          root3 = 1.0e4_DP

          ! ars: evaluate root1
          func1 = aa + bb*root1 + cc*root1*root1 + dd*root1*root1*root1 +&
               ee*root1*root1*root1*root1

          ! ars: check that the chosen root is an actual minimum
          !      (second derivative criterion)
          if ((2.0_DP*cc+6.0_DP*dd*step+12.0_DP*ee*step*step).gt.0.0_DP) then
             ! ars: it's a minimum :)
             success = .true.
             step_is_a_min = .true.
             step = root1; predicted_functional = func1
          else
             ! ars: it's a maximum :(
             step_is_a_min = .false.
             success = .false.
             step = 1.0e4_DP; predicted_functional = 1.0e6_DP
          end if


       end if

    else
       success = .false.
       step = 1.0e4_DP; predicted_functional = 1.0e6_DP
    end if

    if(pub_debug_on_root) then
       write(stdout,*) "aa = ", aa
       write(stdout,*) "bb = ", bb
       write(stdout,*) "cc = ", cc
       write(stdout,*) "dd = ", dd
       write(stdout,*) "ee = ", ee
       write(stdout,*) "alpha = ", alpha
       write(stdout,*) "beta = ", beta
       write(stdout,*) "gamma = ", gamma
       write(stdout,*) "rr = ", rr
       write(stdout,*) "rr_squared = ", rr_squared
       write(stdout,*) "rr_sign = ", sgn_rr
       write(stdout,*) "qq = ", qq
       write(stdout,*) "qq_cubed = ", qq_cubed
       write(stdout,*) "theta = ", theta
       write(stdout,*) "cap_A = ", cap_A
       write(stdout,*) "cap_B = ", cap_B
       write(stdout,*) "all_roots_are_real = ", all_roots_are_real
       write(stdout,*) "root1 = ", root1
       write(stdout,*) "func1 = ", func1
       write(stdout,*) "root2 = ", root2
       write(stdout,*) "func2 = ", func2
       write(stdout,*) "root3 = ", root3
       write(stdout,*) "func3 = ", func3
       write(stdout,*) "step = ", step
       write(stdout,*) &
            "predicted_functional = ", predicted_functional
       write(stdout,*) "step_is_a_min = ", step_is_a_min
       write(stdout,*) "success = ", success
       write(stdout,*) &
            "DEBUG: leaving edft_fit_fourth_order_poly_5points"
    end if

  end subroutine edft_fit_fourth_order_poly_5points

  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine edft_check_params(current_edft_maxit)

    !========================================================================!
    ! Validates input parameters.                                            !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                           !
    !========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout
    use rundat, only: pub_edft_max_step

    implicit none

    integer, intent(inout) :: current_edft_maxit

    logical, save :: success = .false.

    if (.not.success) then

       ! ars: check correct number of inner loop iterations
       if (current_edft_maxit.lt.0) then

          if (pub_on_root) then
             write(stdout,'(a)') &
                  "WARNING: edft_maxit must be an integer greater than zero"
             write(stdout,'(a)') &
                  "         setting to safe value edft_maxit = 10"
          end if
          current_edft_maxit = 10

       end if

       ! ars: check retrial iterations
       if (pub_edft_max_step.lt.0.0_DP) then

          if (pub_on_root) then
             write(stdout,'(a)') &
                  "WARNING: edft_max_steps must be positive"
             write(stdout,'(a)') &
                  "         setting to safe value edft_max_step = 1.0"
          end if
          pub_edft_max_step = 1.0_DP

       end if

    end if

  end subroutine edft_check_params


  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine edft_print_evals(edft, iter, print_all, newfile, outunit)

    !========================================================================!
    ! Prints the eigenvalues of the density kernel and Hamiltonian.          !
    ! For fractional occupancies, it prints only the smeared ones.           !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                           !
    ! Revisited by Alvaro Ruiz Serrano in October 2012 for improvements.     !
    !========================================================================!

    use constants, only: stdout, VERBOSE
    use rundat, only: pub_output_detail, pub_rootname, pub_num_spins,  &
         pub_num_kpoints, PUB_1K, pub_devel_code, pub_edft_spin_fix, &
         pub_edft_grand_canonical, pub_debug_on_root
    use utils, only: utils_open_unit_check, utils_close_unit_check,    &
         utils_unit, utils_abort, utils_devel_code

    implicit none

    ! ars: Arguments
    type(EDFT_MODEL), intent(in) :: edft
    integer, optional, intent(in) :: iter
    logical, optional, intent(in) :: print_all, newfile
    integer, optional, intent(in) :: outunit

    ! Local variables
    integer :: iorb, myunit, ierr, is
    integer :: first, last, top, bottom, lower_homo, higher_lumo, diff
    integer, save :: ncalls = 0
    character(len=10) :: ncalls_char
    character(len=100) :: outfile, fmt
    character(len=27) :: gap, fl
    character(len=6) :: sep
    logical :: do_print_all, do_newfile
    integer :: end_of_range
    integer :: levels_around_fermi
    integer,dimension(pub_num_spins) :: homo_orb, homo_orb_finite
    logical,dimension(2) :: found_gap
    logical :: printed_interm, exist
    real(kind=dp), dimension(pub_num_spins) :: zerok_fermi

    if (pub_debug_on_root) write(stdout,'(a)') "DEBUG: Entering edft_print_evals"

    found_gap = .false.

    ! jme: KPOINTS_DANGER
    ! First k-point enforced in all uses of n_occ via PUB_1K
    if (pub_num_kpoints /= PUB_1K) then
       call utils_abort('Function edft_print_evals not yet compatible with more&
            & than one k-point.')
    end if

    ! ars: bump up number of calls
    ncalls = ncalls + 1
    write(ncalls_char,'(i3)') ncalls

    ! ab: text to be written in the fermi level
    if (pub_edft_grand_canonical) then
       fl = "        Fermi level       "
    else
       fl = " Finite temp. Fermi level "
    end if

    ! ars: text to be written in the gap
    if (.not. pub_edft_grand_canonical) then
       gap = "- Gap at zero temperature -"
    end if

    sep = "------"

    ! ars: determine optionals
    if (present(print_all)) then
       do_print_all = print_all
    else
       do_print_all = .false.
    end if
    if (present(newfile)) then
       do_newfile = newfile
    else
       do_newfile = .false.
    end if
    if (present(outunit)) then
       ! ars: get unit and open it
       myunit = outunit
    else
       ! ars: write on the screen
       myunit = stdout
    end if

    ! ars: set filename
    if (do_newfile) then
       outfile = trim(adjustl(pub_rootname))//"_iter"//&
            trim(adjustl(ncalls_char))//".occ"
    else
       outfile = trim(adjustl(pub_rootname))//".occ"
    end if

    ! ars: open unit
    if (myunit.ne.stdout) then

       inquire(file=trim(outfile), exist=exist)
       if (exist) then

          open(unit=myunit, status='old' , form="formatted" , &
               position='append' ,file=trim(outfile) ,action="write", &
               iostat=ierr)

       else

          open(unit=myunit, form="formatted" ,file=trim(outfile), &
               action="write", iostat=ierr)

       end if

       call utils_open_unit_check('edft_print_evals',outfile,ierr)
    end if

    ! ars: print banner on screen
    if ((pub_output_detail.ge.VERBOSE).and.(myunit.ne.stdout)) &
         write(stdout,'(/,a29,a42)',advance='no') &
         "Writing occupancies to file: ", outfile

    ! kkbd: This doesn't really need to be an input variable
    !       but should be devel code at least
    levels_around_fermi = utils_devel_code(5,'EDFT','print_levels',pub_devel_code,no_bcast=.true.)
    !levels_around_fermi = 5

    ! kkbd: Increase range to encompass all MOs if we're printing everything
    if (do_print_all) levels_around_fermi = edft%num

    call edft_zerok_fermi(zerok_fermi, edft%s_occ, edft%h_evals)

    do is=1,pub_num_spins
       do iorb=1,edft%num
          if (edft%h_evals(iorb,is) > zerok_fermi(is)) then
             ! kkbd: previous orbital was HOMO of this spin channel
             homo_orb(is) = iorb - 1
             exit
          end if
       end do
    end do

    ! ab: finite temperature homo based on finite temperature Fermi level
    do is=1,pub_num_spins
       do iorb=1,edft%num
          if (edft%h_evals(iorb,is) > edft%fermi(is)) then
             ! kkbd: previous orbital was HOMO of this spin channel
             homo_orb_finite(is) = iorb - 1
             exit
          end if
       end do
    end do

    ! ars: print header
    if (present(iter)) then
       write(myunit,'(a, i3)') "Iteration ", iter
    end if

    if (pub_num_spins.eq.1) then ! ars: only one spin channel
       fmt = '(i8, " | ", 2(1x,f15.10), " | ")'

       write(myunit,'(a)') "                           Spin 1          "
       write(myunit,'(a)') "     Orb      H-eigenvalues     Occupancies"

       ! kkbd: Print first eigenvalue and the gap if we need it.
       write(myunit, fmt) 1, edft%h_evals(1,1), edft%occ(1,1)

       if(.not.pub_edft_grand_canonical) then
          if (homo_orb(1) == 1) write(myunit,'(17x,a27)') gap
       end if

       ! ab: similarly print finite temperature fermi level
       if (homo_orb_finite(1) == 1) write(myunit,'(17x,a27)') fl

       if (MINVAL(homo_orb_finite)-levels_around_fermi > 2) write(myunit,'(26x,a6)') sep

    else if (pub_num_spins.eq.2) then ! ars: two spin channels
       fmt = '(i8, " | ", 2(1x,f15.10), " | ", 2(1x,f15.10), " | ")'

       write(myunit,'(a)') "                           Spin 1           |&
            &                 Spin 2           |"
       write(myunit,'(a)') "     Orb |    H-eigenvalues     Occupancies |&
            &    H-eigenvalues     Occupancies |"

        write(myunit,fmt) 1, edft%h_evals(1,1), edft%occ(1,1), &
            edft%h_evals(1,2), edft%occ(1,2)

        ! kkbd: If we're unfortunate enough to have our fermi level after the first orb
        if (.not.pub_edft_grand_canonical) then
           if ((homo_orb(1) == 1) .and. (homo_orb(2) == 1)) then
              write(myunit,'(17x,a27," | ",4x,a27)')  gap, gap
           else if (homo_orb(1) == 1) then
              write(myunit,'(16x, a27)') gap
           else if (homo_orb(2) == 1) then
              write(myunit,'(51x,a27)') gap
           end if
        end if
        ! ab: print finite temperature fermi level
        if ((homo_orb_finite(1) == 1) .and. (homo_orb_finite(2) == 1)) then
           write(myunit,'(17x,a27," | ",4x,a27)')  fl, fl
        else if (homo_orb_finite(1) == 1) then
           write(myunit,'(16x, a27)') fl
        else if (homo_orb_finite(2) == 1) then
           write(myunit,'(51x,a27)') fl
        end if


        ! kkbd: Unless the fist orbital is within 5 levels of a fermi level,
        !       print the separator
        if (MINVAL(homo_orb_finite)-levels_around_fermi > 2) then
           write(myunit,'(26x,a6,29x,a6)') sep, sep
        end if
    end if

    printed_interm = .false.

    ! kkbd: Print the intermediate orbitals
    do iorb = MAX(MINVAL(homo_orb_finite)-levels_around_fermi+1,2), &
              MIN(MAXVAL(homo_orb_finite)+levels_around_fermi,edft%num-1)

       ! kkbd: Determine if we're in the range of orbitals we want to print.
       if ((iorb >= MINVAL(homo_orb_finite)-levels_around_fermi+1 .and. &
            iorb <= MINVAL(homo_orb_finite)+levels_around_fermi) &
         .or. (iorb >= MAXVAL(homo_orb_finite)-levels_around_fermi+1 .and. &
               iorb <= MAXVAL(homo_orb_finite)+levels_around_fermi)) then
          ! kkbd: iorb is an orbital we want to print
          if (pub_num_spins == 1) then
             write(myunit, fmt) iorb, edft%h_evals(iorb,1), edft%occ(iorb,1)
             ! If we're at the 0 K fermi level, print
             if (.not.pub_edft_grand_canonical) then
                if (iorb == homo_orb(1)) write(myunit,'(17x,a27)') gap
             end if
             ! If we are at the finite temperature fermi level, print
             if (iorb == homo_orb_finite(1)) write(myunit,'(17x,a27)') fl
          else
             write(myunit, fmt) iorb, edft%h_evals(iorb,1), edft%occ(iorb,1), &
                 edft%h_evals(iorb,2), edft%occ(iorb,2)
             ! If we're at a 0 K fermi level, print
             if (.not.pub_edft_grand_canonical) then
                if ((homo_orb(1) == iorb) .and. (homo_orb(2)) == iorb) then
                   write(myunit,'(17x,a27," | ",4x,a27)')  gap, gap
                else if (homo_orb(1) == iorb) then
                   write(myunit,'(16x, a27)') gap
                else if (homo_orb(2) == iorb) then
                   write(myunit,'(51x,a27)') gap
                end if
             end if
             ! ab: If we are at a finite temperature fermi level, print
             if ((homo_orb_finite(1) == iorb) .and. (homo_orb_finite(2)) == iorb) then
                write(myunit,'(17x,a27," | ",4x,a27)')  fl, fl
             else if (homo_orb_finite(1) == iorb) then
                write(myunit,'(16x, a27)') fl
             else if (homo_orb_finite(2) == iorb) then
                write(myunit,'(51x,a27)') fl
             end if
          end if
       else if ((.not.printed_interm) .and. &
               (iorb > MINVAL(homo_orb_finite)+levels_around_fermi) .and. &
               (iorb < MAXVAL(homo_orb_finite)-levels_around_fermi)) then
          ! kkbd: iorb is between levels and we haven't printed the separator yet
          printed_interm = .true.
          if (pub_num_spins == 1) then
             write(myunit,'(26x,a6)') sep
          else
             write(myunit,'(26x,a6,29x,a6)') sep, sep
          end if
       end if
    end do

    if (pub_num_spins.eq.1) then ! ars: only one spin channel
       ! kkbd: Print final eigenvalue and the separator if we need it.
       if (MAXVAL(homo_orb_finite)+levels_around_fermi < edft%num-1) write(myunit,'(26x,a6)') sep

       write(myunit, fmt) edft%num, edft%h_evals(edft%num,1), edft%occ(edft%num,1)

    else if (pub_num_spins.eq.2) then ! ars: two spin channels
       ! kkbd: Unless the last orbital is within 5 levels of a fermi level,
       !       print the separator
       if (MAXVAL(homo_orb_finite)+levels_around_fermi < edft%num-1) write(myunit,'(26x,a6,29x,a6)') sep, sep

       write(myunit,fmt) edft%num, edft%h_evals(edft%num,1), edft%occ(edft%num,1), &
           edft%h_evals(edft%num,2), edft%occ(edft%num,2)
    end if

    ! gab: Print footer
    if (present(iter)) then
       write(myunit,'(a, i3)') "End Occupancy", iter
    end if

    ! ars: close unit happily
    if (myunit.ne.stdout) then
       close(unit=myunit, iostat=ierr)
       call utils_close_unit_check('edft_print_evals','outfile',ierr)
    end if

    ! ars: done!
    if ((pub_output_detail.ge.VERBOSE).and.(myunit.ne.stdout)) &
         write(stdout,'(a9,/)') " ... done"

    if (pub_debug_on_root) write(stdout,'(a)') "DEBUG: Leaving edft_print_evals"

  end subroutine edft_print_evals

  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine edft_zerok_fermi(zerok_fermi, tot_occ, evals)
    !========================================================================!
    ! Idenities the 0K Fermi level based on evals and occs and some context. !
    !------------------------------------------------------------------------!
    ! Intent(in   )                                                          !
    !  - tot_occ:   Total occupancy of each spin channel                     !
    !  - evals: Eigenvalues of each spin channel                             !
    ! Intent(  out)                                                          !
    !  - zerok_fermi: 0K Fermi level of each spin channel                    !
    !------------------------------------------------------------------------!
    ! Written by Kevin Duff, Apr 2018                                        !
    !========================================================================!

    use rundat, only: pub_num_spins, pub_edft_spin_fix
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_heapsort

    ! Arguments
    real(kind=dp), dimension(pub_num_spins),   intent(in   ) :: tot_occ
    real(kind=dp), dimension(:,:), intent(in   ) :: evals
    real(kind=dp), dimension(pub_num_spins),   intent(  out) :: zerok_fermi

    ! Local variables
    real(kind=dp), dimension(:), allocatable :: all_evals
    integer, dimension(:), allocatable :: temp_idx
    integer :: is, ierr, nint_nelec, num

    ! Parameters
    real(kind=DP), parameter :: integer_tol = 1E-10

    ! kkbd: Determine the 0K Fermi Energy
    if (pub_edft_spin_fix /= 0) then
        ! If we have fixed spin, the 0K fermi level is obtained from filling
        ! the energy levels of each spin channel in order.
        do is=1,pub_num_spins
           if (modulo(abs(tot_occ(is) - nint(tot_occ(is))),1.0_DP) < integer_tol) then
              ! kkbd: We have about integer occupancy
              nint_nelec = nint(tot_occ(is))
              zerok_fermi(is) = (evals(nint_nelec + 1,is) + &
                                 evals(nint_nelec,is))/2.0_dp
           else
              ! kkbd: Non-integer occupancy. Fermi level is equal to the partially
              !       filled energy level
              zerok_fermi(is) = evals(int(tot_occ(is)) + 1, is)
           end if
        end do
     else
        ! Free spin - 0K fermi level is obtained from filling the lowest num_electrons
        ! energy levels across both spin channels.
        num = size(evals(:,1))

        allocate(all_evals(pub_num_spins*num),&
                 temp_idx(pub_num_spins*num),stat=ierr)
        call utils_alloc_check("edft_zerok_fermi","all_evals",ierr)

        ! Combine all evals and sort them
        do is=1,pub_num_spins
           all_evals((is-1)*num+1:is*num) = evals(:,is)
        end do

        call utils_heapsort(pub_num_spins*num,all_evals,temp_idx)

        if (modulo(abs(sum(tot_occ) - nint(sum(tot_occ))),1.0_DP) < integer_tol) then
           ! kkbd: Integer number of electrons in the system.
           nint_nelec = nint(sum(tot_occ))
           zerok_fermi = (all_evals(nint_nelec+1) + all_evals(nint_nelec))/2.0_dp
        else
           ! kkbd: Non-integer number of electrons in the system, Fermi level is equal
           !       to the partially filled energy level
           zerok_fermi = all_evals(int(sum(tot_occ))+1)
        end if

        deallocate(all_evals,temp_idx,stat=ierr)
        call utils_dealloc_check("edft_zerok_fermi","all_evals",ierr)
     end if

  end subroutine edft_zerok_fermi

  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine edft_spam3_matrix_from_eigs(mat, eigenvecs, eigenvals,&
       num, overlap)

    !==========================================================================!
    ! ars: this subroutine builds any matrix in SPAM3 format, provided its     !
    !      eigenvalues and eigenvectors.                                       !
    ! jcap: modified so this now builds a matrix in SPAM3_EMBED format
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   mat  (output)     : The SPAM3 matrix to construct                      !
    !   eigenvecs (input) : The eigenvectors (rotation MO to NGWF)             !
    !   eigenvals (input) : The eigenvalues                                    !
    !   num (input)       : Total number of eigenstates                        !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                             !
    ! Modified for embedding by Joseph Prentice, May 2018                      !
    !==========================================================================!

    use constants, only: DP, stdout
    use dense, only: DEM, dense_create, dense_destroy, dense_put_element, &
         dense_product, dense_convert
    use rundat, only: pub_debug_on_root
    use sparse_embed, only: SPAM3_EMBED

    implicit none

    ! ars: Arguments
    integer,       intent(in   ) :: num
    type(SPAM3_EMBED), intent(inout) :: mat
    type(DEM),     intent(in   ) :: eigenvecs
    real(kind=DP), intent(in   ) :: eigenvals(num)
    type(SPAM3_EMBED), optional, intent(in) :: overlap


    ! ars: local variables
    integer :: icol
    type(DEM) :: eval_dens, evec_eval, mat_dens
    type(DEM) :: tc_buffer, overlap_dense
    logical :: do_tc
    ! agrecocmplx
    logical :: loc_cmplx
    character :: opB_loc

    if (pub_debug_on_root) write(stdout,'(a)') &
         "DEBUG: Entering edft_spam3_matrix_from_eigs"

    ! agrecocmplx
    loc_cmplx = mat%p%iscmplx

    ! agrecocmplx
    if (loc_cmplx) then
       opB_loc = 'C'
    else
       opB_loc = 'T'
    end if

    ! ars: check optionals
    if (present(overlap)) then
       do_tc = .true.
    else
       do_tc = .false.
    end if

    ! ars: create DEM buffers
    ! agrecocmplx
    call dense_create(evec_eval, num, num, iscmplx=loc_cmplx)
    call dense_create(eval_dens, num, num, iscmplx=loc_cmplx)
    call dense_create(mat_dens, num, num, iscmplx=loc_cmplx)

    ! ars: put eigenvals in diagonal matrix
    do icol=1,num
       ! agrecocmplx
       if (loc_cmplx) then
          call dense_put_element(cmplx(eigenvals(icol),kind=DP), eval_dens, icol, icol)
       else
          call dense_put_element(eigenvals(icol), eval_dens, icol, icol)
       end if
    end do

    ! ars: build mat = evec * eval * evec^t
    call dense_product(evec_eval, eigenvecs, eval_dens, &
         opA='N',opB='N',first_k=1,last_k=num)
    ! agrecocmplx: need to take the hermitian instead of the transpose
    ! for complex matrices? I think so.....
    call dense_product(mat_dens, evec_eval, eigenvecs, &
         opA='N',opB=opB_loc,first_k=1,last_k=num)

    ! ars: apply tensorial correction
    if (do_tc) then
       ! agrecocmplx
       call dense_create(tc_buffer, num, num, iscmplx=loc_cmplx)
       call dense_create(overlap_dense, num, num, iscmplx=loc_cmplx)
       call dense_convert(overlap_dense, overlap)
       call dense_product(tc_buffer, mat_dens, overlap_dense,&
            opA='N',opB='N',first_k=1,last_k=num)
       call dense_product(mat_dens, overlap_dense, tc_buffer, &
            opA='N',opB='N',first_k=1,last_k=num)
       call dense_destroy(overlap_dense)
       call dense_destroy(tc_buffer)
    end if

    ! ars: convert to SPAM3
    ! jcap: now SPAM3_EMBED
    call dense_convert(mat, mat_dens)


    ! ars: destroy DEM buffers
    call dense_destroy(evec_eval)
    call dense_destroy(eval_dens)
    call dense_destroy(mat_dens)

    if (pub_debug_on_root) write(stdout,'(a)') &
         "DEBUG: Leaving edft_spam3_matrix_from_eigs"

  end subroutine edft_spam3_matrix_from_eigs

  !.............................................................................

  subroutine edft_dem_matrix_from_eigs(mat, eigenvecs, eigenvals,&
       num, overlap)


    !==========================================================================!
    ! ars: this subroutine builds any matrix in DEM format, provided its       !
    !      eigenvalues and eigenvectors.                                       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   mat  (output)     : The DEM matrix to construct                        !
    !   eigenvecs (input) : The eigenvectors (rotation MO to NGWF)             !
    !   eigenvals (input) : The eigenvalues                                    !
    !   num (input)       : Total number of eigenstates                        !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2012.                             !
    ! Modified for embedding by Joseph Prentice, May 2018                      !
    !==========================================================================!

    use constants, only: DP, stdout
    use dense, only: DEM, dense_create, dense_destroy, dense_put_element, &
         dense_product, dense_convert
    use rundat, only: pub_debug_on_root
    use sparse_embed, only: SPAM3_EMBED

    implicit none

    ! ars: Arguments
    integer,       intent(in   ) :: num
    type(DEM),     intent(inout) :: mat
    type(DEM),     intent(in   ) :: eigenvecs
    real(kind=DP), intent(in   ) :: eigenvals(num)
    type(SPAM3_EMBED), optional, intent(in) :: overlap


    ! ars: local variables
    integer :: icol
    type(DEM) :: eval_dens, evec_eval, tc_buffer, overlap_dense
    logical :: do_tc
    ! agrecocmplx
    logical :: loc_cmplx
    character :: opB_loc

    if (pub_debug_on_root) write(stdout,'(a)') &
         "DEBUG: Entering edft_dem_matrix_from_eigs"

    ! agrecocmplx
    loc_cmplx = mat%iscmplx

    ! agrecocmplx
    if (loc_cmplx) then
       opB_loc = 'C'
    else
       opB_loc = 'T'
    end if

    ! ars: check optionals
    if (present(overlap)) then
       do_tc = .true.
    else
       do_tc = .false.
    end if

    ! ars: create DEM buffers
    ! agrecocmplx
    call dense_create(evec_eval, num, num, iscmplx=loc_cmplx)
    call dense_create(eval_dens, num, num, iscmplx=loc_cmplx)

    ! ars: put eigenvals in diagonal matrix
    do icol=1,num
       ! agrecocmplx
       if (loc_cmplx) then
          call dense_put_element(cmplx(eigenvals(icol),kind=DP), eval_dens, icol, icol)
       else
          call dense_put_element(eigenvals(icol), eval_dens, icol, icol)
       end if
    end do

    ! ars: build mat = evec * eval * evev^t
    call dense_product(evec_eval, eigenvecs, eval_dens, &
         opA='N',opB='N',first_k=1,last_k=num)
    ! agrecocmplx: need to take hermitian instead of transpose
    ! for complex matrices? I think so.....
    call dense_product(mat, evec_eval, eigenvecs, &
         opA='N',opB=opB_loc,first_k=1,last_k=num)

    ! ars: apply tensorial correction
    if (do_tc) then
       ! agrecocmplx
       call dense_create(tc_buffer, num, num, iscmplx=loc_cmplx)
       call dense_create(overlap_dense, num, num, iscmplx=loc_cmplx)
       call dense_convert(overlap_dense, overlap)
       call dense_product(tc_buffer, mat, overlap_dense,&
            opA='N',opB='N',first_k=1,last_k=num)
       call dense_product(mat, overlap_dense, tc_buffer, &
            opA='N',opB='N',first_k=1,last_k=num)
       call dense_destroy(overlap_dense)
       call dense_destroy(tc_buffer)
    end if


    ! ars: destroy DEM buffers
    call dense_destroy(evec_eval)
    call dense_destroy(eval_dens)

    if (pub_debug_on_root) write(stdout,'(a)') &
         "DEBUG: Leaving edft_dem_matrix_from_eigs"

  end subroutine edft_dem_matrix_from_eigs



  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine edft_ngwf_gradient(pmat, qmat, denskern, ham, inv_overlap)

    !==========================================================================!
    ! Calculates the NGWF gradient during EDFT calculations.                   !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in October 2012.                          !
    ! Modified for embedding by Joseph Prentice, May 2018                      !
    ! KHS_ expression for Q matrix changed to S_HK by Arihant Bhandari, Sep 2020
    ! based on derivation by C.-K. Skylaris and several tests.                 !
    !==========================================================================!

    use constants, only: DP, stdout
    use rundat, only: pub_num_spins, pub_debug_on_root
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product, &
         sparse_embed_copy, sparse_embed_scale, sparse_embed_axpy

    implicit none

    type(SPAM3_EMBED), intent(inout) :: pmat(1:pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: qmat(1:pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: denskern(1:pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: ham(1:pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: inv_overlap

    integer :: is
    type(SPAM3_EMBED) :: hk

    if (pub_debug_on_root) write(stdout,'(a)') &
         "DEBUG: Entering edft_ngwf_gradient"

    ! ars: create temporary structures
    call sparse_embed_create(hk, ham(1), denskern(1),arrname='hk')


    do is = 1, pub_num_spins

       ! ars: copy denskern to pmat
       call sparse_embed_copy(pmat(is), denskern(is))

       ! ars: calculate qmat
       call sparse_embed_product(hk, ham(is), denskern(is))
       call sparse_embed_product(qmat(is), inv_overlap, hk)
       call sparse_embed_scale(qmat(is), -1.0_DP)

    end do

    ! ars: destroy structures
    call sparse_embed_destroy(hk,arrname='hk')

    if (pub_debug_on_root) write(stdout,'(a)') &
         "DEBUG: Leaving edft_ngwf_gradient"

  end subroutine edft_ngwf_gradient


  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine edft_reorthogonalise_mo(mo, denskern, ham, overlap, num, residue)

    !==========================================================================!
    ! Re-orthogonalises the Kohn-Sham molecular orbitals.                      !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in October 2012.                          !
    ! Modified for embedding by Joseph Prentice, May 2018                      !
    !==========================================================================!

    use constants, only: stdout, DP
    use dense, only: DEM, dense_create, dense_destroy, dense_copy
    use rundat, only: pub_debug_on_root
    use sparse_embed, only: SPAM3_EMBED

    implicit none

    ! ars: arguments
    type(DEM), intent(inout) :: mo
    type(SPAM3_EMBED), intent(in) :: denskern, ham, overlap
    integer, intent(in) :: num
    real(kind=DP), intent(out) :: residue

    ! ars: local
    logical :: orthogonal, mo_reorthogonalised
    real(kind=DP) :: residue_before, residue_after
    type(DEM) :: mo_buffer
    ! agrecocmplx
    logical :: loc_cmplx

    if (pub_debug_on_root) write(stdout,'(a)') &
         "DEBUG: EDFT: Entering edft_reorthogonalise_mo"

    ! agrecocmplx
    loc_cmplx = overlap%p%iscmplx

    ! ars: only calculate the residue of non-orthogonality
    residue = 0.0_DP
    mo_reorthogonalised = .false.
    orthogonal = edft_check_orthogonality(residue_before, mo, denskern, &
         ham, overlap, num)

    if (pub_debug_on_root) write(stdout,*) "Orthogonality before = ", &
         residue_before

    ! ars: note that the non-orthonality residue is measured based on the
    !      contribution to the free energy Lagrangian.
    ! ars: we assume that if this residue is less than the free energy threshold
    !      then we shouldn't worry about it.
    if (.not.orthogonal) then

       ! ars: save current MOs in buffer for safety
       ! agrecocmplx
       call dense_create(mo_buffer, num, num, iscmplx=loc_cmplx)
       call dense_copy(mo_buffer, mo)

       ! ars: use Lowdin method for orthogonalisation
       call edft_lowdin_orthogonalise(mo, overlap, num)

       ! ars: and calculate residue again for tracking
       orthogonal = edft_check_orthogonality(residue_after, mo, denskern, &
            ham, overlap, num)

       if (abs(residue_before).gt.abs(residue_after)) then

          ! ars: Lowdin orthogonalisation succeeded
          mo_reorthogonalised = .true.
          residue = residue_after

       else

          ! ars: Lowdin orthogonalisation failed
          mo_reorthogonalised = .false.
          call dense_copy(mo, mo_buffer)
          residue = residue_before

       end if

       if (pub_debug_on_root) then
          write(stdout,*) "Orthogonality after = ", residue
          write(stdout,*) "mo_reorthogonalised = ", mo_reorthogonalised
       end if

       call dense_destroy(mo_buffer)

    end if

    if (pub_debug_on_root) write(stdout,'(a)') &
         "DEBUG: EDFT: Leaving edft_reorthogonalise_mo"

  end subroutine edft_reorthogonalise_mo

  !-----------------------------------------------------------------------------

  subroutine edft_lowdin_orthogonalise(mo, overlap, num)

    !==========================================================================!
    ! Performs Lowdin orthogonalisation of the Kohn-Sham molecular orbitals.   !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in October 2012.                          !
    ! Modified for embedding by Joseph Prentice, May 2018                      !
    !==========================================================================!

    use constants, only: DP, stdout
    use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
         dense_normal_eigensolve, dense_product, dense_scale, dense_copy
    use rundat, only: pub_eigensolver_orfac, pub_eigensolver_abstol, &
         pub_debug_on_root
    use sparse_embed, only: SPAM3_EMBED
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! ars: arguments
    type(DEM), intent(inout) :: mo
    type(SPAM3_EMBED), intent(in) :: overlap
    integer, intent(in) :: num

    ! ars: local
    type(DEM) :: overlap_dense, sm, msm, evecs
    real(kind=DP), allocatable :: evals(:)
    integer :: iii, ierr
    ! agrecocmplx
    logical :: loc_cmplx
    character :: opA_loc

    if (pub_debug_on_root) write(stdout,'(a)') &
         "DEBUG: EDFT: Entering edft_lowdin_orthogonalise"

    ! agrecocmplx
    loc_cmplx = overlap%p%iscmplx

    ! agrecocmplx
    if (loc_cmplx) then
       opA_loc = 'C'
    else
       opA_loc = 'T'
    end if

    ! ars: allocate evals
    allocate(evals(1:num), stat=ierr)
    call utils_alloc_check('edft_lowdin_orthogonalise', 'evals', ierr)
    evals(:) = 0.0_DP

    ! ars: create buffers
    ! agrecocmplx
    call dense_create(overlap_dense, num, num, iscmplx=loc_cmplx)
    call dense_convert(overlap_dense, overlap)
    call dense_create(sm, num, num, iscmplx=loc_cmplx)
    call dense_create(msm, num, num, iscmplx=loc_cmplx)
    call dense_create(evecs, num, num, iscmplx=loc_cmplx)

    ! ars: construct MO overlap
    call dense_product(sm, overlap_dense, mo, 'N', 'N', 1, num)
    ! agrecocmplx: need to take hermitian instead of transpose for
    ! complex matrices? I think so....
    call dense_product(msm, mo, sm, opA_loc, 'N', 1, num)

    ! ars: diagonalise
    call dense_normal_eigensolve(num,evals,msm,evecs,&
         arg_orfac=pub_eigensolver_orfac, arg_abstol=pub_eigensolver_abstol)

    ! ars: construct s^{-1/2} and store in msm
    call dense_scale(msm, 0.0_DP)
    do iii = 1, num
       evals(iii) = 1.0_DP/sqrt(evals(iii))
    end do
    call edft_matrix_from_eigs(msm, evecs, evals, num)

    ! ars: copy old, non-orthogonal MOs in sm
    call dense_copy(sm, mo)

    ! ars: build new, orthogonalised MOs
    call dense_product(mo, sm, msm, 'N', 'N', 1, num)

    ! ars: destroy buffers
    call dense_destroy(evecs)
    call dense_destroy(msm)
    call dense_destroy(sm)
    call dense_destroy(overlap_dense)

    ! ars: deallocate evals
    deallocate(evals, stat=ierr)
    call utils_dealloc_check('edft_lowdin_orthogonalise', 'evals', ierr)

    if (pub_debug_on_root) write(stdout,'(a)')&
         "DEBUG: EDFT: Leaving edft_lowdin_orthogonalise"

  end subroutine edft_lowdin_orthogonalise

  !-----------------------------------------------------------------------------

  logical function edft_check_orthogonality(residue, mo, denskern, ham,&
       overlap, num)

    !==========================================================================!
    ! Checks that the MOs are orthogonal by evaluating the constraint in the   !
    ! Lagrangian.                                                              !
    ! L_o = tr[ H (XSK - K) ]                                                  !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in October 2012.                          !
    ! Modified for embedding by Joseph Prentice, May 2018                      !
    !==========================================================================!

    use constants, only: DP, stdout
    use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
         dense_product
    use rundat, only: pub_edft_free_energy_thres, pub_spin_fac, &
         pub_debug_on_root
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product,&
         sparse_embed_axpy, sparse_embed_trace

    implicit none

    ! ars: arguments
    real(kind=DP), intent(out) :: residue
    type(DEM), intent(in) :: mo
    type(SPAM3_EMBED), intent(in) :: denskern, ham, overlap
    integer, intent(in) :: num

    ! ars: local variables
    real(kind=DP) :: threshold
    type(DEM) :: xprod_dense
    type(SPAM3_EMBED) :: xprod, xs, xsk
    real(kind=DP), parameter :: tolerance = 1.0e-7_DP
    ! agrecocmplx
    logical :: loc_cmplx
    character :: opB_loc

    if (pub_debug_on_root) write(stdout,'(a)') &
         "DEBUG: EDFT: Entering edft_check_orthogonality"

    ! agrecocmplx
    loc_cmplx = overlap%p%iscmplx

    ! agrecocmplx
    if (loc_cmplx) then
       opB_loc = 'C'
    else
       opB_loc = 'T'
    end if

    ! ars: init
    edft_check_orthogonality = .false.
    residue = 0.0_DP
    if (pub_edft_free_energy_thres.lt.0.0_DP) then
       threshold = tolerance
    else
       threshold = pub_edft_free_energy_thres
    end if

    ! ars: create temporaries
    call sparse_embed_create(xprod, denskern, arrname='xprod')
    call sparse_embed_create(xs, xprod, overlap, arrname='xs')
    call sparse_embed_create(xsk, xs, denskern, arrname='xsk')
    ! agrecocmplx
    call dense_create(xprod_dense, num, num, iscmplx=loc_cmplx)

    ! ars: calculate xprod and matrices
    ! agrecocmplx: need to take hermitian instead of transpose
    ! for complex matrices? I think so.....
    call dense_product(xprod_dense, mo, mo, 'N', opB_loc, 1, num)
    call dense_convert(xprod, xprod_dense)
    call sparse_embed_product(xs, xprod, overlap)
    call sparse_embed_product(xsk, xs, denskern)
    call sparse_embed_axpy(xsk, denskern, -1.0_DP)

    ! ars: calculate residue - this will affect the Lagrangian
    call sparse_embed_trace(residue, ham, xsk)
    residue = residue * pub_spin_fac


    ! ars: note that the non-orthonality residue is measured based on the
    !      contribution to the free energy Lagrangian.
    ! ars: we assume that if this residue is less than the free energy threshold
    !      then we shouldn't worry about it.
    if (abs(residue).lt.threshold) then
       edft_check_orthogonality = .true.
    else
       edft_check_orthogonality = .false.
    end if


    ! ars: destroy temporaries
    call dense_destroy(xprod_dense)
    call sparse_embed_destroy(xprod,arrname='xprod')
    call sparse_embed_destroy(xsk,arrname='xsk')
    call sparse_embed_destroy(xs,arrname='xs')

    if (pub_debug_on_root) then
       write(stdout,'(a,f22.14)')&
            "DEBUG: EDFT: Orthogonality residue = ", residue
       write(stdout,'(a)')&
            "DEBUG: EDFT: Leaving edft_check_orthogonality"
    end if

  end function edft_check_orthogonality


  !-----------------------------------------------------------------------------

  subroutine edft_print_line_search(name, step, lagrangian, predicted)

    !==========================================================================!
    ! Prints the result of the current line search step on the screen.         !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in October 2012.                          !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout

    implicit none

    ! ars: arguments
    character(len=9), intent(in) :: name
    real(kind=DP), intent(in) :: step, lagrangian
    real(kind=DP), intent(in), optional :: predicted

    if (present(predicted)) then
       if (pub_on_root) then
          write(stdout,'(". ", a9, 2x, f10.6, 2x, f22.14, 2x, f22.14)')&
               name, step, lagrangian, predicted
       end if
    else
       if (pub_on_root) then
          write(stdout,'(". ", a9, 2x, f10.6, 2x, f22.14)')&
               name, step, lagrangian
       end if
    end if

  end subroutine edft_print_line_search


  !-----------------------------------------------------------------------------


  subroutine edft_print_qc(edft, commutator, rms_gradient, &
       delta_lagrangian, delta_energy, delta_entropy, iter)

    !==========================================================================!
    ! Prints the QC test results for the Ensemble-DFT module.                  !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in October 2012.                          !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: DP
    use rundat, only: pub_num_spins, pub_edft_rms_gradient_thres, &
         pub_edft_grand_canonical, pub_charge
    use utils, only: utils_qc_print

    implicit none

    type(EDFT_MODEL), intent(in) :: edft
    real(kind=DP), intent(in) :: commutator, rms_gradient
    real(kind=DP), intent(in) :: delta_lagrangian, delta_energy, delta_entropy
    integer :: iter
    character(len=3) :: spin_string

    integer :: is

    if(pub_on_root) then

       call utils_qc_print('commutator',commutator)
       if (pub_edft_rms_gradient_thres.gt.0.0_DP) then
          call utils_qc_print('rms_gradient',rms_gradient)
       end if
       call utils_qc_print('delta_lagrangian',delta_lagrangian)
       call utils_qc_print('delta_energy',delta_energy)
       call utils_qc_print('delta_entropy',delta_entropy)
       call utils_qc_print('iter',iter)
       call utils_qc_print('edft%num',edft%num)
       call utils_qc_print('edft%nbands',edft%nbands)
       call utils_qc_print('edft%free_energy',edft%free_energy)
       call utils_qc_print('edft%energy',edft%energy)
       call utils_qc_print('edft%entropy',edft%entropy)

       if (pub_edft_grand_canonical) then
          call utils_qc_print('edft%nelec',edft%nelec)
          call utils_qc_print('pub_charge',pub_charge)
       end if

       do is = 1, pub_num_spins

          write(spin_string,'(a1,i1,a1)') '(', is, ')'
          call utils_qc_print('orth_residue'//spin_string,edft%orth_residue(is))
          call utils_qc_print('fermi'//spin_string,edft%fermi(is))
          call utils_qc_print('integrated_ne'//spin_string,edft%integrated_ne(is))

       end do

    end if


  end subroutine edft_print_qc


  !-----------------------------------------------------------------------------


  real(kind=DP) function edft_trial_step(npts, xx, yy)

    !==========================================================================!
    ! Calculates the next trial step based on the history of previous steps,   !
    ! using the golden section polynomial search method.                       !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in January 2013.                          !
    !==========================================================================!

    use rundat, only: pub_edft_max_step, pub_debug_on_root
    use constants, only: DP, stdout

    implicit none

    ! ars: arguments
    integer, intent(in) :: npts
    real(kind=DP), intent(in) :: xx(npts), yy(npts)

    ! ars: local variables
    real(kind=DP) :: min_xx, left_xx, right_xx
    integer :: min_pos, nn
    logical :: left_found, right_found
    ! ars: golden section
    real(kind=DP), parameter :: golden_section = (3.0_DP-sqrt(5.0_DP))*0.5_DP

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') &
         "DEBUG: EDFT: Entering edft_trial_step"

    ! ars: find the position of the minimum in the yy array
    min_pos = minloc(yy(:),dim=1)

    ! ars: find the value of xx at the minimum
    min_xx=xx(min_pos)

    ! ars: find left neighbour of the minimum (if any)
    left_xx=0.0_DP; left_found=.false.
    do nn = 1, npts
       if (xx(nn).lt.min_xx) then
          if (xx(nn).gt.left_xx) then
             left_xx = xx(nn)
             left_found = .true.
          end if
       else
          left_xx = 0.0_DP
          left_found = .false.
       end if
    end do

    ! ars: find right neighbour of the minimum (if any)
    right_xx=pub_edft_max_step; right_found=.false.
    do nn = 1, npts
       if (xx(nn).gt.min_xx) then
          if (xx(nn).lt.right_xx) then
             right_xx = xx(nn)
             right_found = .true.
          end if
       else
          right_xx = pub_edft_max_step
          right_found = .false.
       end if
    end do

    ! ars: calculate step - move from the minimum into the largest subinterval
    if ( (left_found).and.(right_found) ) then
       if ( (abs(right_xx-min_xx)).gt.(abs(left_xx-min_xx)) ) then
          edft_trial_step = min_xx + golden_section*abs(right_xx-min_xx)
       else
          edft_trial_step = min_xx - golden_section*abs(left_xx-min_xx)
       end if
    else
       edft_trial_step = min_xx + golden_section*abs(right_xx-min_xx)
    end if

    if(pub_debug_on_root) then
       write(stdout,'(a,f22.14)') "DEBUG: min_xx = ", min_xx
       write(stdout,'(a,f22.14)') "DEBUG: left_xx = ", left_xx
       write(stdout,'(a,l1)') "DEBUG: left_found = ", left_found
       write(stdout,'(a,f22.14)') "DEBUG: right_xx = ", right_xx
       write(stdout,'(a,l1)') "DEBUG: right_found = ", right_found
       write(stdout,'(a,f22.14)')&
            "DEBUG: edft_trial_step = ", edft_trial_step
       write(stdout,'(a)') "DEBUG: EDFT: Leaving edft_trial_step"
    end if

  end function edft_trial_step

!|------------------------------------------------------------------------|!
!|                   DIIS ROUTINES FROM KERNEL_DIIS_MOD                   |!
!|------------------------------------------------------------------------|!

  subroutine edft_diis_pulay_mix(mat_in, mat_his, err_his, &
       inv_overlap, ientry, trial, ham)

    !==============================================================!
    ! This subroutine performs Pulay mixing of matrices.           !
    !--------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010.                !
    ! Modified for embedding by Joseph Prentice, May 2018          !
    ! Tranferred to EDFT by Gabriel Bramley, Nov 2020              !
    !==============================================================!


    use constants, only: DP
    use rundat, only: pub_num_spins, pub_edft_ham_diis_size
    use sparse_embed, only: SPAM3_EMBED
    use utils, only: utils_alloc_check, utils_dealloc_check
    use ngwf_representation, only: NGWF_HAM

    implicit none

    ! ars: Arguments
    type(SPAM3_EMBED), intent(inout) :: mat_in(pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: mat_his(pub_edft_ham_diis_size,pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: err_his(pub_edft_ham_diis_size,pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: inv_overlap
    integer,     intent(in   ) :: ientry
    real(kind=DP),     intent(in   ) :: trial
    type(NGWF_HAM),    intent(in   ) :: ham

    ! ars: local variables
    real(kind=DP), allocatable :: diis_mat(:,:,:)
    real(kind=DP), allocatable :: diis_coeffs(:,:)
    integer :: ierr, diis_size

    diis_size = pub_edft_ham_diis_size

    ! ars: allocate diis_mat and coeffs
    allocate(diis_mat(1:ientry+1, 1:ientry+1, 1:pub_num_spins),stat=ierr)
    call utils_alloc_check('edft_diis_pulay_mix', 'diis_mat', ierr)
    allocate(diis_coeffs(1:ientry+1, 1:pub_num_spins), stat=ierr)
    call utils_alloc_check('edft_diis_pulay_mix', 'diis_coeffs', ierr)


    ! ars: calculate DIIS matrix based on the error vectors.
    call edft_diis_pulay_matrix(diis_mat, err_his, inv_overlap, ientry)

    ! ars: get DIIS coefficients
    call edft_diis_coeffs(diis_coeffs, diis_mat, ientry)
    call edft_diis_coeffs_summ(diis_coeffs, ientry)

    ! ars: mix kernels
    call edft_mix_ham(mat_in, mat_his, err_his, diis_coeffs, ientry, trial, ham)

    ! ars: deallocate Bmat and coeffs
    deallocate(diis_coeffs, stat=ierr)
    call utils_dealloc_check('edft_diis_pulay_mix', 'diis_coeffs', ierr)
    deallocate(diis_mat, stat=ierr)
    call utils_dealloc_check('edft_diis_pulay_mix', 'diis_mat', ierr)

  end subroutine edft_diis_pulay_mix


  subroutine edft_diis_pulay_matrix(Bmat,residues,overlap,ientry)

    !==============================================================!
    ! This subroutine calculates matrix B of the system of linear  !
    ! equations to be solved in order to mix the Hamiltonian       !
    ! according to the Pulay method.                               !
    !--------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in April 2010.                 !
    ! Modified for embedding by Joseph Prentice, May 2018          !
    ! Tranferred for EDFT by Gabriel Bramley, Nov. 2020            !
    !==============================================================!

    use constants, only: DP
    use rundat, only: pub_num_spins, pub_edft_ham_diis_size
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_transpose, &
         sparse_embed_product, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_trace


    implicit none

    ! ars: arguments
    real(kind=DP), intent(  out) :: Bmat(:,:,:)
    type(SPAM3_EMBED),   intent(in   ) :: residues(pub_edft_ham_diis_size, pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: overlap
    integer,       intent(in   ) :: ientry

    ! ars: local variables
    integer     :: is, ii, jj
    type(SPAM3_EMBED) :: res_ov, res_t, res_ov_res_t


    ! ars: init spam3buffer blocking scheme
    call sparse_embed_create(res_ov, residues(1,1), overlap, arrname='res_ov')
    call sparse_embed_create(res_t, residues(1,1), arrname='res_t')
    call sparse_embed_create(res_ov_res_t, res_ov, res_t, arrname='res_ov_res_t')


    do is = 1, pub_num_spins
       do jj = 1, ientry
          do ii = 1, ientry

             ! ars: res_ov = R_i x S
             call sparse_embed_product(res_ov, residues(ii, is), overlap)

             ! ars: res_t = (R_j)^t
             call sparse_embed_transpose(res_t,residues(jj, is))

             ! ars: res_ov_res_t = (R_i) x S x (R_j)^t
             call sparse_embed_product(res_ov_res_t, res_ov, res_t)

             ! ars: Bmat(is,ii,jj) = tr[(R_i) x S x (R_j)^t x S]
             call sparse_embed_trace(Bmat(ii,jj,is),res_ov_res_t,overlap)

          end do
       end do
    end do

    ! ars: destroy SPAM3 buffers
    call sparse_embed_destroy(res_ov,arrname='res_ov')
    call sparse_embed_destroy(res_t,arrname='res_t')
    call sparse_embed_destroy(res_ov_res_t,arrname='res_ov_res_t')

    ! ars: add space for Langrange multiplier
    Bmat(1:ientry,ientry+1,:) = -1.0_DP
    Bmat(ientry+1,1:ientry,:) = -1.0_DP
    Bmat(ientry+1,ientry+1,:) = 0.0_DP

  end subroutine edft_diis_pulay_matrix

  subroutine edft_diis_coeffs(coeffs,Bmat,ientry)

    !==============================================================!
    ! This subroutine finds the coefficients for the Pulay mixing  !
    ! after solving a system of linear equations Bd=0.             !
    !--------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in May 2010.                  !
    ! Transferred to EDFT by Gabriel Bramley, Nov 2020             !
    !==============================================================!

    use comms, only: pub_on_root, comms_barrier
    use constants, only: DP, stdout
    use rundat, only: pub_num_spins


    implicit none

    ! ars: arguments
    real(kind=DP), intent(inout) :: coeffs(:,:)
    real(kind=DP), intent(in   ) :: Bmat(:,:,:)
    integer      , intent(in   ) :: ientry

    ! LAPACK subroutines
    external :: dgetrf, dgetrs

    ! ars: library wrapper variables
    integer :: INFO, LDA, LDB, M, N, NRHS
    integer :: IPIV(1:ientry+1)
    real(kind=DP) :: A(ientry+1,ientry+1)
    real(kind=DP) :: B(ientry+1)
    character(LEN=1) :: TRANS

    ! ars: local variables
    integer :: is

    ! ars: set up parameters for LAPACK DGETRF
    M = ientry+1
    N = ientry+1
    LDA = ientry+1
    ! ars: set up parameters for LAPACK DGETRS
    TRANS ='N'
    NRHS = 1
    LDB = ientry+1

    ! ars: solve linear system
    do is = 1, pub_num_spins

       ! ars: call LAPACK DGETRF
       A = Bmat(:,:,is)
       call DGETRF(M, N, A, LDA, IPIV, INFO)
       if((INFO.ne.0).and.pub_on_root) write(stdout,*) &
            "Error calling DGETRF. INFO|node = ", INFO

       ! ars: call LAPACK DGETRS
       B(1:ientry) = 0.0_DP
       B(ientry+1) = -1.0_DP
       call DGETRS(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       if((INFO.ne.0).and.pub_on_root) write(stdout,*) &
            "Error calling DGETRS. INFO|node = ", INFO

       ! ars: set coeffs(is) before exit
       coeffs(:, is) = B(:)

    end do

    call comms_barrier

  end subroutine edft_diis_coeffs

  subroutine edft_diis_coeffs_summ(coeffs, ientry)

    !==================================================================!
    ! This subroutine prints a summary of the coefficients that are    !
    ! used during the density kernel mixing.                           !
    !------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2010.                     !
    ! Transferred to EDFT by Gabriel Bramley, Nov 2020                 !
    !==================================================================!

    use comms, only : pub_on_root
    use constants, only: DP, stdout, VERBOSE
    use rundat, only: pub_output_detail, pub_num_spins

    implicit none

    real(kind=DP), intent(in   ) :: coeffs(:,:)
    integer,       intent(in   ) :: ientry

    integer :: icntr
    character(len=15) :: mat_name
    character(len=80) :: fmt


    ! ars: prepare format
    write(fmt,'(a)') '(/, a28, a30)'

    ! ars: banner
    if (pub_on_root.and.pub_output_detail.ge.VERBOSE) then
       write(stdout,fmt) &
            "Performing  Pulay mixing of ", &
            "Hamiltonian with coefficients:"
       if (pub_num_spins.eq.1) then
          write(stdout,'(10x,a20)') "              Spin 1"
          do icntr = 1, ientry
             write(stdout,'(10x,f20.12)') coeffs(icntr,1)
          end do
          write(stdout,'(a10,f20.12)',advance='no') "Lambda--> ", coeffs(ientry+1,1)
       else if (pub_num_spins.eq.2) then
          write(stdout,'(10x,a20,1x,a20)') "              Spin 1", &
               "              Spin 2"
          do icntr = 1, ientry
             write(stdout,'(10x,f20.12,1x,f20.12)') coeffs(icntr,1), coeffs(icntr,2)
          end do
          write(stdout,'(a10,f20.12,1x,f20.12)',advance='no') &
               "Lambda--> ", coeffs(ientry+1,1), coeffs(ientry+1,2)
       end if
       write(stdout,'(a,/)') "     ... done"
    end if

  end subroutine edft_diis_coeffs_summ


  subroutine edft_mix_ham(next_ham_in, ham_in, ham_er, coeffs,&
       ientry, trial, ham)

    !==================================================================!
    ! This subroutine mixes the Hamiltonian according to the           !
    ! Pulay method.                                                    !
    !------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in May 2010.                      !
    ! Modified for embedding by Joseph Prentice, May 2018              !
    ! Adapted for EDFT by Gabriel Bramley, Nov 2020                    !
    !==================================================================!

    use constants, only: DP
    use rundat, only: pub_num_spins, pub_edft_ham_diis_size
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_create, sparse_embed_destroy
    use ngwf_representation, only: NGWF_HAM

    implicit none

    type(SPAM3_EMBED),   intent(inout) :: next_ham_in(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: ham_in(pub_edft_ham_diis_size, pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: ham_er(pub_edft_ham_diis_size, pub_num_spins)
    real(kind=DP), intent(in   ) :: coeffs(:,:), trial
    integer,       intent(in   ) :: ientry
    type(NGWF_HAM), intent(in   ) :: ham

    real(kind=DP) :: resid_scale
    integer :: is, counter
    type(SPAM3_EMBED), allocatable :: ham_err_sum(:)

    ! gab: Assign from global variables
    resid_scale = trial

    ! gab: Buffer for error sum test
    allocate(ham_err_sum(1:pub_num_spins))
    do is = 1, pub_num_spins
       ! agrecocmplx
       call sparse_embed_create(ham_err_sum(is),ham%ham(is),arrname='ham_err_sum')
    end do

    ! ars: H_in(n+1) = 0
    do is = 1, pub_num_spins
       call sparse_embed_scale(next_ham_in(is),0.0_DP)
       call sparse_embed_scale(ham_err_sum(is),0.0_DP)
    end do

    ! gab: H_in(n+1) = sum_i d_i ( H_in(i) + (\alpha*precond) R[H_in(i)] )
    do is = 1, pub_num_spins
       do counter = 1, ientry
          call sparse_embed_axpy(next_ham_in(is),ham_in(counter,is),&
               coeffs(counter,is))
       end do
    end do

    do is = 1, pub_num_spins
       do counter = 1, ientry
          call sparse_embed_axpy(ham_err_sum(is),ham_er(counter,is),&
               coeffs(counter,is))
       end do
    end do

    do is = 1, pub_num_spins
       call sparse_embed_axpy(next_ham_in(is),ham_err_sum(is),resid_scale)
    end do

    do is = 1, pub_num_spins
       call sparse_embed_destroy(ham_err_sum(is),arrname='ham_err_sum')
    end do
    deallocate(ham_err_sum)

  end subroutine edft_mix_ham

  subroutine edft_diis_shift(mat_his, his_size, nspins)

    !================================================================!
    ! This subroutine shifts the matrices in arrays in case iter has !
    ! reached pub_edft_ham_diis_size value.                            !
    !----------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in May 2010.                    !
    ! Modified for embedding by Joseph Prentice, May 2018            !
    ! Transferred to EDFT by Gabriel Bramley, Nov. 2020              !
    !================================================================!

    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy

    implicit none

    ! ars: arguments
    integer,     intent(in   ) :: his_size, nspins
    type(SPAM3_EMBED), intent(inout) :: mat_his(his_size, nspins)

    ! ars: local variables
    integer :: icntr, is


    ! ars: move old entries one position back
    do is = 1, nspins
       do icntr = 1, his_size-1
          call sparse_embed_copy(mat_his(icntr, is), mat_his(icntr+1, is))
       end do
    end do


  end subroutine edft_diis_shift

  subroutine edft_diis_find_ientry(ientry, iter, shifted)

    !==============================================================!
    ! This subroutine finds the next ientry to work with.          !
    !--------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in May 2010.                  !
    ! Tranferred to EDFT, Gabriel Bramley Nov. 2020                !
    !==============================================================!

    use comms, only: comms_barrier
    use rundat, only: pub_edft_ham_diis_size

    implicit none

    ! ars: arguments
    integer, intent(  out) :: ientry
    integer, intent(in   ) :: iter
    logical, intent(in   ) :: shifted


    ! ars: find next position
    if (shifted) then
       ientry = pub_edft_ham_diis_size
    else
       ientry = iter+1
    end if
    call comms_barrier


  end subroutine edft_diis_find_ientry

  subroutine edft_pulay_history_flush(edft, ham, cmplx)

    !==============================================================!
    ! Initialises or flushes hamiltonian history for Pulay mixing. !
    !--------------------------------------------------------------!
    ! Written by Gabriel Bramley Nov. 2020                         !
    !==============================================================!

    use constants, only: stdout
    use rundat, only: pub_num_spins, pub_edft_ham_diis_size
    use sparse_embed, only: sparse_embed_destroy, sparse_embed_create, &
         SPAM3_EMBED
    use ngwf_representation, only: NGWF_HAM
    use utils, only: utils_dealloc_check, utils_abort, utils_alloc_check

    implicit none

    type(EDFT_MODEL),  intent(inout) :: edft
    type(NGWF_HAM), intent(in   ) :: ham
    logical, intent(in   ) :: cmplx

    logical :: loc_cmplx
    integer :: ierr, is, icntr, diis_size

    diis_size = pub_edft_ham_diis_size
    loc_cmplx = cmplx

    do is = 1, pub_num_spins
       do icntr = 1, diis_size
          call sparse_embed_destroy(edft%ham_in_his(icntr, is),arrname='edft%ham_in_his')
       end do
    end do

    do is = 1, pub_num_spins
       do icntr = 1, diis_size
          call sparse_embed_destroy(edft%ham_err_his(icntr, is),arrname='edft%ham_err_his')
       end do
    end do

    do is = 1, pub_num_spins
       do icntr = 1, diis_size
          ! agrecocmplx
          call sparse_embed_create(edft%ham_in_his(icntr, is), ham%ham(is),arrname='edft%ham_in_his')
       end do
    end do

    ! gab: Pulay history arrays
    do is = 1, pub_num_spins
       do icntr = 1, diis_size
          ! agrecocmplx
          call sparse_embed_create(edft%ham_err_his(icntr, is), ham%ham(is),arrname='edft%ham_err_his')
       end do
    end do

  end subroutine edft_pulay_history_flush

end module ensemble_dft
