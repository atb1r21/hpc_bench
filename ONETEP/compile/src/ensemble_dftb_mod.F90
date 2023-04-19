
module ensemble_dftb

  !============================================================================!
  ! Following subroutines are extracted from ensemble_DFT and adapted for DFTB !
  ! by Arihant Bhandari in Feb 2022.                                           !
  ! Changes in ensemble_dft_mod should also be reflected here.                 !
  !============================================================================!

  use ensemble_dft

  implicit none

  private

  public :: edftb_calculate
  public :: edftb_evaluate_step

contains

  subroutine edftb_calculate(edft, denskern, ham, rep, &
       ngwf_basis, mdl, current_edft_maxit, restart_rootname, force_no_IO)

    !==========================================================================!
    ! Self-consistent optimisation of the occupancies and the unitary rotations!
    ! that transform from NGWF to Kohn-Sham representation, based on the       !
    ! Ensemble-DFTB.                                                           !
    !==========================================================================!

    use comms, only: pub_on_root, pub_total_num_procs
    use constants, only: DP, stdout, paw_en_size, VERBOSE, NORMAL, BRIEF, &
         very_big_double, DFTB_GFN0
    use dense, only: DEM, dense_create, dense_destroy, dense_convert
    use dftb, only: dftb_calculate_ham, dftb_eht_energy
    use ensemble_dft_type, only: EDFT_MODEL
    use function_basis, only: FUNC_BASIS
    use kernel, only: kernel_rescale_spam3
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use restart, only: restart_kernel_write, restart_kernel_read, &
         restart_hamiltonian_write, restart_hamiltonian_read
    use rundat, only: pub_debug_on_root, pub_edft_smearing_width, &
         pub_output_detail, pub_edft_max_step, &
         pub_write_converged_dk_ngwfs, pub_write_denskern, &
         pub_read_hamiltonian, pub_write_hamiltonian, &
         pub_print_qc, pub_num_spins, pub_read_denskern, &
         pub_spin_fac, pub_foe, pub_num_kpoints, PUB_1K, &
         pub_edft_rms_gradient_thres, pub_dftb_method, &
         pub_dense_foe, pub_inner_loop_iteration, pub_edft_spin_fix, &
         pub_real_spin, pub_edft_trial_step, pub_edft_update_scheme, &
         pub_edft_ham_diis_size, pub_devel_code, pub_edft_grand_canonical
    use smearing_operator, only: calculate_fermi_entropy
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
         sparse_embed_create, sparse_embed_destroy, sparse_embed_trace,&
         sparse_embed_copy, sparse_embed_axpy, &
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
    type(FUNC_BASIS), intent(in   ) :: ngwf_basis(:)
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

    real(kind=DP), allocatable :: tmp_vec(:,:)

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

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering edftb_calculate'

    ! jme: KPOINTS_DANGER
    ! First k-point enforced in all uses of n_occ via PUB_1K
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine edftb_calculate not yet compatible with more&
         & than one k-point.')

    ! ars: start timer
    call timer_clock('edftb_calculate',1)

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
            &Ensemble-DFTB optimisation -------------------------"


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
    call utils_assert(edft%allocd,"Error in edftb_calculate:&
         & please allocate 'edft' before calling this routine")

    ! ars: check input parameters
    call edft_check_params(current_edft_maxit)

    ! ars: allocate structures
    allocate(initial_ham(1:pub_num_spins), stat=ierr)
    call utils_alloc_check('edftb_calculate', 'initial_ham', ierr)
    allocate(trial_ham(1:pub_num_spins), stat=ierr)
    call utils_alloc_check('edftb_calculate', 'trial_ham', ierr)


    ! ars: create DEM structures
    do is = 1, pub_num_spins
       call sparse_embed_create(initial_ham(is), ham%ham(is))
       call sparse_embed_create(trial_ham(is), ham%ham(is))
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
    edft_iter = 0
    if (pub_dftb_method == DFTB_GFN0) then ! ab: no deltas or convergence
       deltal = 0.0_DP
       deltae = 0.0_DP
       deltas = 0.0_DP
       delta_fermi(:) = 0.0_DP
    else
       deltal = very_big_double         ! jd: We don't initialise to zero, but
       deltae = very_big_double         !     rather to a big value to prevent
       deltas = very_big_double         !     premature convergence after step #0
       delta_fermi(:) = very_big_double !     when other criteria are satisfied, but
    end if                              !     a delta is not yet available

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
          ! rab207: read denskern
          call restart_kernel_read(denskern,read_cond=.false.)

          ! ab: build the hamiltonian
          call dftb_calculate_ham(mdl, ham, rep, ngwf_basis)
          call dftb_eht_energy(denskern, ham, eaccepted)

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
          call utils_alloc_check('edftb_calculate','tmp_spam',ierr)
          call sparse_embed_create(tmp_spam(1),ham%ham(1), edft%rot)
          call sparse_embed_create(tmp_spam(2),edft%rot)
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
          call sparse_embed_destroy(tmp_spam(1))
          call sparse_embed_destroy(tmp_spam(2))
          deallocate(tmp_spam,stat=ierr)
          call utils_dealloc_check('edftb_calculate','tmp_spam',ierr)
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
       call edftb_evaluate_step(t0, trial_ham, initial_ham, ham, rep, edft, &
            mdl, denskern, ngwf_basis,  &
            e0, s0, l0, foe_line_search_entropy, foe_entropy_approx,    &
            bounds, fermi_tol, tmp_fermi, denswork, tmp_s, foe_mode, ientry)

       do is = 1, pub_num_spins
          call sparse_embed_trace(edft%integrated_ne(is),&
               denskern%m(is,PUB_1K),rep%overlap)
       end do

       ! ars: rescale density kernel
       ! JA: Changed to s_occ on restart because this is definitely set up.
       t_occ(:,PUB_1K) = edft%s_occ(:)
       if (pub_debug_on_root)write(stdout,*) &
               &"DEBUG: EDFT: Rescaling denskern to ",t_occ(:,PUB_1K)
       call kernel_rescale_spam3(denskern, rep%overlap, t_occ, &
            silent=.true.)
       do is = 1, pub_num_spins
          call sparse_embed_trace(edft%integrated_ne(is),&
               denskern%m(is,PUB_1K),rep%overlap)
       end do

       if (pub_foe) then
          saccepted = s0
       end if
    end if ! edft%initialised

    ! ars: Update density-dependent terms with current NGWFs and K.
    ! ars: Calculate energy.
    ! ars: The Hamiltonian is not updated.
    ! ab: check again
    call dftb_calculate_ham(mdl, ham, rep, ngwf_basis)
    call dftb_eht_energy(denskern, ham, eaccepted)

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

    if (pub_dftb_method == DFTB_GFN0) then ! ab: print iteration and terminate

       call edft_print_iter(edft_iter, edft,  &
            eaccepted, saccepted, laccepted, deltal, &
            rms_gradient, commutator, delta_fermi, 0.0_DP)

    else ! ab: perform edft minimization

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
          ! ab: calculate using DFTB subroutines
          call dftb_calculate_ham(mdl, ham, rep, ngwf_basis)
          call dftb_eht_energy(denskern, ham, eaccepted)

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

             call edftb_evaluate_step(t1, trial_ham, initial_ham, ham, rep, edft, &
                  mdl, denskern, ngwf_basis, &
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
             call edftb_evaluate_step(t1, trial_ham, initial_ham, ham, rep, edft, &
                  mdl, denskern, ngwf_basis, &
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

             call edftb_evaluate_step(t2, trial_ham, initial_ham, ham, rep, edft, &
                  mdl, denskern, ngwf_basis,  &
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

             call edftb_evaluate_step(t3, trial_ham, initial_ham, ham, rep, edft, &
                  mdl, denskern, ngwf_basis, &
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

                call edftb_evaluate_step(t4, trial_ham, initial_ham, ham, rep, edft, &
                     mdl, denskern, ngwf_basis, &
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

                   call edftb_evaluate_step(t5, trial_ham, initial_ham, ham, rep, edft, &
                        mdl, denskern, ngwf_basis, &
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

                      call edftb_evaluate_step(t6, trial_ham, initial_ham, ham, rep, edft, &
                           mdl, denskern, ngwf_basis, &
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

                call edftb_evaluate_step(taccepted, trial_ham, initial_ham, ham, rep, &
                     edft, mdl, denskern, ngwf_basis, &
                     eaccepted, saccepted, laccepted, foe_line_search_entropy, &
                     foe_entropy_approx, bounds, fermi_tol, tmp_fermi, denswork, &
                     tmp_s, foe_mode, ientry)

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

    end if

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
       call sparse_embed_destroy(initial_ham(is))
       call sparse_embed_destroy(trial_ham(is))
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
    call utils_dealloc_check('edftb_calculate', 'initial_ham', ierr)
    deallocate(trial_ham, stat=ierr)
    call utils_dealloc_check('edftb_calculate', 'trial_ham', ierr)

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
    call timer_clock('edftb_calculate',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving edftb_calculate'

  end subroutine edftb_calculate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine edftb_evaluate_step(t, trial_ham, initial_ham, ham, rep, edft, &
       mdl, denskern, ngwf_basis, e, s, l, foe_line_search_entropy, &
       foe_entropy_approx, bounds, fermi_tol, tmp_fermi, denswork, tmp_s, &
       foe_mode, ientry)

    !========================================================================!
    ! Subroutine for a Hamiltonian trial step.                               !
    !========================================================================!

    use constants, only: stdout, DP
    use dense,  only: DEM
    use dftb, only: dftb_eht_energy
    use ensemble_dft_type, only: EDFT_MODEL
    use function_basis, only: FUNC_BASIS
    use kernel, only: kernel_rescale_spam3
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use sparse_embed, only: sparse_embed_copy, sparse_embed_scale, &
         sparse_embed_axpy, SPAM3_EMBED, SPAM3_EMBED_ARRAY
    use rundat, only: pub_edft_smearing_width, &
         pub_spin_fac, PUB_1K, pub_num_spins, pub_foe, &
         pub_edft_update_scheme, pub_debug_on_root
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
    type(SPAM3_EMBED_ARRAY),  intent(inout) :: denskern
    type(FUNC_BASIS),   intent(in   ) :: ngwf_basis(:)
    real(kind=DP),      intent(inout) :: tmp_fermi(:), &
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
    if(pub_debug_on_root)write(stdout,*) &
            &"DEBUG: EDFT: Rescaling denskern to ",t_occ(:,PUB_1K)
    call kernel_rescale_spam3(denskern, rep%overlap, t_occ, &
         silent=.true.)

    ! ars: Update density-dependent terms with current NGWFs and K.
    ! ars: Calculate energy.
    ! ars: Avoid updating the Hamiltonian.
    call dftb_eht_energy(denskern, ham, e)

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

  end subroutine edftb_evaluate_step

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ensemble_dftb
