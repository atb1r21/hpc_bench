! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!          Mermin Density kernel optimisation module             !
!                                                                !
! This module implements a Mermin (Free Energy) approach to      !
! optimisei  the density kernel K for a fixed set of NGWFs.      !
!----------------------------------------------------------------!
! Written by Poli Emiliano,                                      !
! Based on the lnv module by Chris-Kriton Skylaris, Arash Mostofi!
! and Peter Haynes.                                              !
!================================================================!

module mermin

  implicit none

  private

  public :: mermin_denskernel_optimise_cg
  public :: mermin_calculate_mu
  public :: mermin_reset_trial_step_to_initial
  public :: mermin_ngwf_gradient
  public :: mermin_lagrangian
  public :: mermin_switch_mu

  logical, save :: mermin_reset_trial_step_to_initial = .true.

contains


  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////


  subroutine mermin_denskernel_optimise_cg(aux, pur_denskern, ham, &
       lhxc_fine, mu, total_energy, rep, ngwf_basis, hub_proj_basis, hub, &
       mdl, hfxstate, mermin_threshold, current_maxit_mermin, conv_status, &
       dfdtau_fine, kpt, force_no_IO)

    !========================================================================!
    ! This subroutine implements the optimisation of  the density kernel K   !
    ! for a fixed set of NGWFs using the "Mermin" scheme.                    !
    ! In this module the Free Energy is defined as:                          !
    ! A(phi_ab,K^ab,T)=tr(K^ab <phi_a|'T'+V_ext|phi_b> + E_h[n]+ E_xc[n] - & !
    ! -T*s[{(K^**S_**)}]+ mu*(tr(K^ab S_ab)-N_e)                             !
    ! where:                                                                 !
    ! A(phi_ab,K^ab,T) = Free energy                                         !
    ! phi_a, phi_b = NGWFs                                                   !
    ! K^ab = Density Kernel                                                  !
    ! 'T'= Kinetic operator                                                  !
    ! V_ext = External potential operator                                    !
    ! E_h = Hartree Energy                                                   !
    ! n = Electronic Density                                                 !
    ! E_xc = Exchange-Correlation energy                                     !
    ! T = Temperature                                                        !
    ! s[{(K^**S_**)}] = Chebyshev expansion of the Entropy                   !
    ! mu = Chemical Potential (Defined by Millam Scuseria method)            !
    ! S_ab =overlap matrix                                                   !
    ! N_e = Number of electron                                               !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    ! aux           (inout): Auxiliary density kernel in DKERN format        !
    ! pur_denskern  (inout): Purified density kernel in DKERN format         !
    ! ham           (inout): Hamiltonian wrapper (contains SPAM3 matrices)   !
    ! lhxc_fine       (out): Local-Hartree-Exhange-Correlation potential     !
    ! mu              (out): Chemical potential type parameter               !
    ! total_energy    (out): Total energy                                    !
    ! step rep         (in): NGWF Representation (functions and matrices)    !
    ! ngwf_basis       (in): Function basis type describing the NGWFs.       !
    ! hub_proj_basis   (in): Function basis type describing Hubbard projs    !
    ! grid             (in): The GRID_INFO for all the whole-cell arrays     !
    ! mermin_threshold (in): Convergence threshold for mermin gradient       !
    ! current_maxit_mermin(in): Maximum mermin iterations set                !
    ! dfdtau_fine   (inout): Gradient of XC energy per unit volume wrt       !
    !                        KE density energy (optional)                    !
    !------------------------------------------------------------------------!
    ! This version written by Emiliano Poli.                                 !
    ! Based on the lnv module by Chris-Kriton Skylaris                       !
    !========================================================================!

    use comms, only: pub_on_root, pub_root_proc_id, comms_bcast
    use constants, only: DP, stdout, VERBOSE, NORMAL
    use ensemble_dft, only: edft_print_line_search, edft_fit_second_order_poly_3points, &
         edft_fit_third_order_poly_4points, edft_fit_fourth_order_poly_5points
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_build_matrix, hamiltonian_energy_components
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: DKERN, kernel_workspace_invalidate, &
         kernel_workspace_create, kernel_workspace_destroy, kernel_normalise, &
         kernel_create, kernel_destroy, kernel_rescale
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use restart, only : restart_kernel_write
    use rundat, only: pub_output_detail, pub_mermin_cg_type, pub_write_denskern, &
         pub_mermin_check_trial_steps,pub_mermin_cg_max_step, &
         pub_hub_calculating_u,pub_write_converged_dk_ngwfs,pub_devel_code, &
         pub_cond_calculate, pub_debug_on_root, &
         pub_num_spins, pub_num_kpoints, PUB_1K, &
         pub_dmft_fully_sc, pub_dmft_fully_sc_h, pub_inner_loop_iteration, &
         pub_emft, pub_emft_follow, pub_emft_lnv_only, &
         pub_mermin_smearing_width, pub_check_mermin, pub_mermin_cheb, &
         pub_mermin_temp, pub_mermin_cg_max_step, pub_print_qc
    use services, only: services_flush
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_copy, &
         sparse_embed_destroy, sparse_embed_axpy, sparse_embed_trace, &
         sparse_embed_array_create, sparse_embed_array_copy, &
         sparse_embed_array_axpy, sparse_embed_array_trace, &
         sparse_embed_array_destroy, SPAM3_EMBED_ARRAY, &
         sparse_embed_array_scale, sparse_embed_array_product
    use smearing_operator, only: calculate_fermi_entropy_mermin
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! ** Arguments
    ! Auxiliary density kernel:
    type(DKERN), intent(inout)    :: aux
    ! Purified density kernel:
    type(DKERN), intent(inout)    :: pur_denskern
    ! Hamiltonian
    type(NGWF_HAM), intent(inout) :: ham
    type(MODEL), intent(in)       :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    ! Local PSP, Hartree and XC potential on grid:
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub)
    real(kind=DP), intent(out) :: mu(pub_num_spins, pub_num_kpoints) ! Chemical potential
    real(kind=DP), intent(out) :: total_energy     ! Total energy

    type(NGWF_REP), intent(in) :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(mdl%nsub)  ! rc2013: pass all bases
    type(FUNC_BASIS), intent(in) :: hub_proj_basis(mdl%nsub)
    type(HUBBARD_MODEL), intent(inout) :: hub
    integer, intent(in) :: current_maxit_mermin       ! Max number of iterations
    real(kind=DP), intent(in) :: mermin_threshold     ! mermin threshold
    integer, intent(out)      :: conv_status
    ! JCW: Optional gradient of XC energy per unit volume wrt kinetic energy
    ! JCW: density. Required to evaluate lhxc matrix elements for meta-GGAs.
    real(kind=DP), optional, intent(inout) :: dfdtau_fine(:,:,:,:)
    ! JCW: In general, we have
    ! JCW:    dfdtau_fine(mdl%fine_grid%ld1,&
    ! JCW:        mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    ! JCW: or
    ! JCW:    dfdtau_fine(1,1,1,1)
    ! JCW: when it is unneeded.
    ! agrecokpt: needed in TB method in internal_energy
    real(kind=DP), optional, intent(in) :: kpt(3)
    logical, optional, intent(in) :: force_no_IO

    ! ** Local variables
    ! Trial density kernel
    type(DKERN) :: trial_aux

    ! Conjugate gradient variables
    type(SPAM3_EMBED_ARRAY) :: co_gradient       ! Covariant gradient
    type(SPAM3_EMBED_ARRAY) :: con_gradient      ! Contravariant gradient
    type(SPAM3_EMBED_ARRAY) :: con_direction     ! Contra search direction
    type(SPAM3_EMBED_ARRAY) :: old_con_direction ! Previous contra direction
    type(SPAM3_EMBED_ARRAY) :: old_co_gradient   ! Previous co gradient
    !QC ENTROPY SLOPE
    type(SPAM3_EMBED_ARRAY) :: entr_grad_co      ! Contravariant gradient
    type(SPAM3_EMBED_ARRAY) :: entr_grad_comp    ! Contravariant gradient
    type(SPAM3_EMBED_ARRAY) :: entr_dir_comp     ! Contra search direction
    type(SPAM3_EMBED_ARRAY) :: old_entr_dir_comp     ! Contra search direction
    !QC DFT SLOPE
    type(SPAM3_EMBED_ARRAY) :: dft_grad_co      ! Contravariant gradient
    type(SPAM3_EMBED_ARRAY) :: dft_grad_comp    ! Contravariant gradient
    type(SPAM3_EMBED_ARRAY) :: dft_dir_comp     ! Contra search direction
    type(SPAM3_EMBED_ARRAY) :: old_dft_dir_comp     ! Contra search direction
    !QC MU SLOPE
    type(SPAM3_EMBED_ARRAY) :: mu_grad_co      ! Contravariant gradient
    type(SPAM3_EMBED_ARRAY) :: mu_grad_comp    ! Contravariant gradient
    type(SPAM3_EMBED_ARRAY) :: mu_dir_comp     ! Contra search direction
    type(SPAM3_EMBED_ARRAY) :: old_mu_dir_comp     ! Contra search direction

    real(kind=DP) :: ne_pure          ! Electron number from purified kernel
    real(kind=DP) :: old_total_energy ! Total energy at previous mermin step
    real(kind=DP) :: commutator       ! RMS value of [L,H] commutator
    real(kind=DP) :: rms_gradient     ! RMS value of gradient of E wrt L
    real(kind=DP) :: line_search_step ! Line search optimal step length
    real(kind=DP) :: cg_coeff         ! Conjugate gradients coefficient
    real(kind=DP) :: lhxc_energy      ! Local PSP, Hartree and XC energy
    real(kind=DP) :: mu_term          ! Term associated with chemical potential
    real(kind=DP) :: rms_dkern        ! RMS value of change in density kernel
    real(kind=DP) :: old_g_dot_g      ! Previous scalar product of gradients
    real(kind=DP) :: new_g_dot_g      ! Current scalar product of gradients
    real(kind=DP) :: old_d_dot_g      ! Current scalar product of gradients
    real(kind=DP) :: new_d_dot_g      ! Current scalar product of gradients
    integer :: iteration              ! Iteration number
    integer :: cg_count               ! Number of CG steps since last reset
    integer :: mult
!==== ep: TO BE EVALUATED IN THE FUTURE
    integer, parameter :: cg_max=5   ! Maximum number of CG steps before reset
!==== ep: TO BE EVALUATED IN THE FUTURE
    logical :: converged              ! Flag to indicate convergence
    integer :: ierr
    character :: opB_loc

    ! Line search variables
    real(kind=DP) :: Qinitial         ! Initial value of function
    real(kind=DP) :: QFD_op         ! Initial value of function
    real(kind=DP) :: Qslope           ! Initial slope along search direction
    real(kind=DP) :: Eslope           ! Initial slope along search direction
    real(kind=DP) :: DFTslope           ! Initial slope along search direction
    real(kind=DP) :: Muslope           ! Initial slope along search direction
    real(kind=DP), save :: trial_step ! Trial step length
    real(kind=DP) :: Qtrial           ! Function value at trial step
    real(kind=DP) :: optimal_step     ! Optimal step length
    real(kind=DP) :: Qoptimal         ! Function value at optimal step
    real(kind=DP) :: Qpredict         ! Predicted value at optimal step
    real(kind=DP) :: Qminimum         ! Predicted value of function at minimum
    real(kind=DP) :: total_energy_at_length
    real(kind=DP) :: total_energy_at_trial_length
    real(kind=DP) :: total_energy_at_trial_length2
    real(kind=DP) :: total_energy_at_trial_length3
    real(kind=DP) :: total_energy_at_trial_length4
    real(kind=DP) :: total_energy_at_trial_length5
    real(kind=DP) :: total_energy_at_trial_length6
    real(kind=DP), parameter :: golden_section = (3.0_DP-sqrt(5.0_DP))*0.5_DP !1.88
    real(kind=DP) :: xx2(2), xx3(3), xx4(4)
    real(kind=DP) :: yy2(2), yy3(3), yy4(4)
    real(kind=DP) :: trial_step2, trial_step3, trial_step4, trial_step5, trial_step6
    real(kind=DP) :: predicted_e3, predicted_e4, predicted_e5
    real(kind=DP) :: accepted_e, accepted_step
    logical :: quadratic_success, cubic_success, quartic_success
    logical :: quadratic_accepted,cubic_accepted,quartic_accepted,tiny_accepted
    integer :: npts
    real(kind=DP) :: hubbard_energy   ! Hubbard correction energy

    ! General variables
    real(kind=DP), parameter :: eps=epsilon(1.0_DP)
    real(kind=DP) :: norm_fac(pub_num_spins, pub_num_kpoints) ! Normalisation factor
    real(kind=DP), allocatable :: trace_array(:,:)
    integer :: is                     ! Spin counter
    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecokpt
    real(kind=DP) :: loc_kpt(3)

    ! ndmh: for finite differencing
    integer, parameter :: nfd=10
    real(kind=DP) :: fd_trial_step(nfd)=(/0.000001_DP, 0.00001_DP,0.0001_DP, &
         0.001_DP,0.01_DP,0.1_DP,0.2_DP,0.3_DP,0.4_DP,0.5_DP/)
    integer :: ifd
    real(kind=DP) :: QFD(nfd)
    real(kind=DP) :: fd_slope_op
    real(kind=DP) :: fd_slope(nfd)

    !ep QC components of the gradient variables
    real(kind=DP)   :: steppo
    !QC FD ENTROPY SLOPE
    real(kind=DP) :: ent_energy
    real(kind=DP) :: ent_energy_at_length
    real(kind=DP) :: QFD_ent(nfd)
    real(kind=DP) :: fd_slope_ent
    !QC FD DFT SLOPE
    real(kind=DP) :: norm_energy
    real(kind=DP) :: norm_energy_at_length
    real(kind=DP) :: QFD_dft(nfd)
    real(kind=DP) :: fd_slope_op_dft
    real(kind=DP) :: fd_slope_dft
    !QC FD DFT SLOPE
    real(kind=DP) :: mu_energy
    real(kind=DP) :: mu_energy_at_length
    real(kind=DP) :: QFD_mu(nfd)
    real(kind=DP) :: fd_slope_op_mu
    real(kind=DP) :: fd_slope_mu
    !ep QC components of the gradient variables

    !ep : to be checked if needed in future/cleaner version
    character(len=80) :: filename
    integer           :: io_unit, io_stat
    logical           :: loc_force_no_IO
    !ep : to be checked if needed in future/cleaner version

    !ep: initialise minimisation alg. variables
    trial_step2=0.0_DP; trial_step3=0.0_DP; trial_step4=0.0_DP;
    trial_step5=0.0_DP; trial_step6=0.0_DP
    total_energy_at_length=0.0_DP
    total_energy_at_trial_length=0.0_DP
    total_energy_at_trial_length2=0.0_DP
    total_energy_at_trial_length3=0.0_DP
    total_energy_at_trial_length4=0.0_DP
    total_energy_at_trial_length5=0.0_DP
    total_energy_at_trial_length6=0.0_DP
    accepted_e=0.0_DP;accepted_step=0.0_DP
    predicted_e3=0.0_DP;predicted_e4=0.0_DP;predicted_e5=0.0_DP
    quadratic_success=.false.;cubic_success=.false.;quartic_success=.false.
    quadratic_accepted=.false.;cubic_accepted=.false.
    quartic_accepted=.false.;tiny_accepted=.false.
    steppo=0.000001_DP ! QC variable
    !ep: initialise minimisation alg. variables


    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
         &Entering mermin_denskernel_optimise_cg'

    ! Flush output
    call services_flush

    ! Start the timer
    call timer_clock('mermin_denskernel_optimise_cg', 1)

    ! jme: KPOINTS_DANGER
    ! The algorithms in this module have not been checked for more than 1 k-point
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine mermin_denskernel_optimise_cg not checked yet for more&
         & than one k-point.')

    ! agrecocmplx
    loc_cmplx = rep%overlap%iscmplx

    ! agrecokpt: default is Gamma point
    if (present(kpt)) then
       loc_kpt(:) = kpt(:)
    else
       loc_kpt(:) = 0.0_DP
    end if

    ! jd: Reset trial step for each new energy and force calculation
    if(mermin_reset_trial_step_to_initial) then
       !ep : adaptable step with mermin method due to high chebyshev expansions
       trial_step=1.0_DP/(pub_mermin_cheb*0.5_DP)
       mermin_reset_trial_step_to_initial = .false.
    end if

    if (.not. aux%workspace_created) then
       call kernel_workspace_create(aux, rep%overlap)
    end if


    ! Allocate matrices with the sparsity pattern of the density kernel
    ! agrecocmplx
    call kernel_create(trial_aux, aux%kern%m(1,1)%structure//rep%postfix, &
         is_cmplx=loc_cmplx)

    call sparse_embed_array_create(co_gradient, aux%kern)
    call sparse_embed_array_create(con_gradient, aux%kern)
    call sparse_embed_array_create(con_direction, aux%kern)
    call sparse_embed_array_create(old_con_direction, aux%kern)
    call sparse_embed_array_create(old_co_gradient, aux%kern)

    !ep: INITIALISE COMPONENTS TO THE GRADIENT
    if( pub_print_qc) then
        call sparse_embed_array_create(entr_grad_co, aux%kern)
        call sparse_embed_array_create(entr_grad_comp, aux%kern)
        call sparse_embed_array_create(entr_dir_comp, aux%kern)
        call sparse_embed_array_create(old_entr_dir_comp, aux%kern)
        call sparse_embed_array_create(dft_grad_co, aux%kern)
        call sparse_embed_array_create(dft_grad_comp, aux%kern)
        call sparse_embed_array_create(dft_dir_comp, aux%kern)
        call sparse_embed_array_create(old_dft_dir_comp, aux%kern)
        call sparse_embed_array_create(mu_grad_co, aux%kern)
        call sparse_embed_array_create(mu_grad_comp, aux%kern)
        call sparse_embed_array_create(mu_dir_comp, aux%kern)
        call sparse_embed_array_create(old_mu_dir_comp, aux%kern)
    endif
    !ep: INITIALISE COMPONENTS TO THE GRADIENT

    allocate(trace_array(co_gradient%num_spins,co_gradient%num_kpoints),&
         stat=ierr)
    call utils_alloc_check('mermin_denskernel_optimise_cg','trace_array',ierr)

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    ! Ensure the kernel will be stable before optimisation starts
    ! ebl : the dmft dens kernel spoil flag may be needed here
    if (index(pub_devel_code,'RESCALE_AUX')>0) then
       call kernel_rescale(aux, rep%overlap, rep%n_occ)
    else
       call kernel_normalise(aux, rep%overlap, rep%inv_overlap, rep%n_occ)
    end if

    ! agrecocmplx
    loc_cmplx = aux%kern%m(1,1)%iscmplx

    ! agrecocmplx
    if (loc_cmplx) then
       opB_loc = 'C'
    else
       opB_loc = 'T'
    end if

    ! Initialise variables for conjugate gradients
    old_g_dot_g = 1.0e100_DP           ! Set previous product to dummy value
    line_search_step = 0.0_DP          ! Set line search step length to zero
    cg_count = 0                       ! Reset conjugate gradients counter
    converged = .false.                ! Initially no convergence

    if (pub_on_root .and. (pub_output_detail >= NORMAL)) then

       write(stdout,'(a)')'   iter  |      energy (Eh)       | rms gradient |&
            &  commutator  |   dE (Eh)'
    end if

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! ** Start of main loop
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    do iteration=1,current_maxit_mermin
       ! jd: Globally accessible copy of iteration number, only used for contro-
       !     lling devel_code dumps of grid quantities from various modules
       pub_inner_loop_iteration = iteration

       ! Flush output
       call services_flush

       ! cks: ======================= ENERGY AT POINT 0 =======================
       total_energy = internal_energy(aux,.true.)
       ! cks: ==================== END ENERGY AT POINT 0 ======================

       ! cks: ******** CALCULATE KOHN-SHAM MATRIX AND MU **********************
       ! ddor: Hamiltonian renewal is not required when calculating DFT+U parameter
       if ((.not. pub_hub_calculating_u) .or. ((pub_hub_calculating_u).and. &
            (iteration .eq. 1))) then

          ! lr408: no need to recalculate ham for a conduction calculation
          ! ebl: or for certain dmft calculations
          if (.not.pub_cond_calculate) then
             if(.not.pub_dmft_fully_sc.or.(pub_dmft_fully_sc.and.pub_dmft_fully_sc_h)) then
                ! agrecokpt: each component of ham and overlap contains the
                ! appropriate k-dependent terms at this stage
                call hamiltonian_build_matrix(ham, rep)
             end if
          end if

          !ep: calculate chemical potential
          call mermin_calculate_mu(mu, aux, ham, rep, ngwf_basis)

          !ep: calculate QC contributions to the free energy
          if( pub_print_qc) then
            ent_energy = internal_ent(aux,.true.)
            norm_energy = internal_dft(aux,.true.)
            mu_energy = internal_muen(aux,.true.)
          end if

          do is=1,pub_num_spins
             if (pub_on_root .and. (pub_output_detail >= VERBOSE)) then
                write(stdout,'(a,i1,a,f20.12)') '      ideal mu',is,': ', &
                     mu(is,PUB_1K)
             end if
          end do
       endif
       ! cks: **** END CALCULATE KOHN-SHAM MATRIX AND MU **********************


       ! cks: %%%%%%%%%%%%%%%%%%%%%%%% GRADIENT AT POINT 0 %%%%%%%%%%%%%%%%%%%%
       rms_gradient = internal_gradient_norm()
       ! pdh: hack to make consistent with spin-unpolarised case
       rms_gradient = rms_gradient * pub_num_spins
       ! cks: %%%%%%%%%%%%%%%%%%%%%%% END GRADIENT AT POINT 0 %%%%%%%%%%%%%%%%%


       ! cks: ########################## TEST CONVERGENCE #####################
       converged = internal_test_convergence()
       call comms_bcast(pub_root_proc_id,converged)

       if (converged) then
          exit
       end if
       ! cks: ######################## END TEST CONVERGENCE ###################

       ! Line search is skipped in the last iteration
       if (iteration /= current_maxit_mermin) then

          ! cks: ===========================LINE SEARCH========================
          ! Initial function value including constraint term
          Qinitial = total_energy

          if ( (pub_output_detail >= VERBOSE) .and. pub_on_root) &
             write(stdout,'(a,f20.12)') '             E0: ', Qinitial

          call internal_search_direction ! <-- sets Qslope (JCW)

          if ( (pub_output_detail >= VERBOSE) .and. pub_on_root) &
             write(stdout,'(a,f20.12)') '   RMS gradient: ', rms_gradient


          !ep QC components to the energy and FD gradient at steppo
          if (pub_print_qc) then
             if (pub_mermin_check_trial_steps) then
                call internal_kernel_protected_step(trial_aux, &
                     steppo,aux,con_direction)
             else
                call sparse_embed_array_copy(trial_aux%kern, aux%kern)
                call sparse_embed_array_axpy(trial_aux%kern, &
                     con_direction, steppo)
                call kernel_workspace_invalidate(trial_aux)
             end if
             total_energy_at_length = internal_energy(trial_aux, &
                   .false.)
             ent_energy_at_length = internal_ent(trial_aux,.false.)
             norm_energy_at_length = internal_dft(trial_aux,.false.)
             mu_energy_at_length = internal_muen(trial_aux,.false.)
             fd_slope_op = (total_energy_at_length - total_energy)/steppo
             fd_slope_ent = (ent_energy_at_length - ent_energy)/steppo
             fd_slope_dft = (norm_energy_at_length - norm_energy)/steppo
             fd_slope_mu =  (mu_energy_at_length - mu_energy)/steppo
          end if
          !ep QC components to the energy and FD gradient at steppo

          !ep: FD check for the slope derived by the gradient
          if (index(pub_devel_code,'MERMIN_FD')>0) then

             if ( (pub_output_detail < VERBOSE) .and. pub_on_root) &
                  write(stdout,'(a,f20.12)') '             G0: ',Qslope

             ! ndmh: evaluate and print energy at several FD trial steps
             do ifd=1,nfd
                if (pub_mermin_check_trial_steps) then
                   call internal_kernel_protected_step(trial_aux, &
                        fd_trial_step(ifd),aux,con_direction)
                else
                   ! Take a trial step along the search direction
                   call sparse_embed_array_copy(trial_aux%kern, aux%kern)
                   call sparse_embed_array_axpy(trial_aux%kern, &
                        con_direction, fd_trial_step(ifd))
                   call kernel_workspace_invalidate(trial_aux)
                end if

                ! ndmh: %%%%%%%%%%%%%%%%%% ENERGY AT FD LENGTH %%%%%%%%%%%%%%
                total_energy_at_trial_length = internal_energy(trial_aux, &
                     .false.)
                QFD(ifd) = total_energy_at_trial_length
                ! rc2013: check that the slopes are sane
                fd_slope(ifd) = (QFD(ifd) - Qinitial)/fd_trial_step(ifd)
                !ep Print components to the energy
                if ((pub_output_detail>=VERBOSE).and. &
                     (.not.pub_cond_calculate)) then
                   call hamiltonian_energy_components( &
                        trial_aux%kern%m(:,PUB_1K), rep, mdl, &
                        ngwf_basis, hub_proj_basis, hub, &
                        ham%hfexchange, mu(:,PUB_1K))
                end if
                !ep STILL SOME PROBLEM WITH HAM ENTROPY COMP DOUBLE CHECK

                if (pub_on_root) then
                   write(stdout,'(a,3f20.12)') ' FD step, QFD, FD slope  : ', &
                       fd_trial_step(ifd),QFD(ifd),fd_slope(ifd)
                end if
                ! ndmh: %%%%%%%%%%%%%%%%%% END ENERGY AT FD LENGTH %%%%%%%%%%

             end do
          end if
          !ep: FD check for the slope derived by the gradient

          if (pub_mermin_check_trial_steps) then
             call internal_kernel_protected_step(trial_aux, trial_step, &
                  aux,con_direction)
          else
             ! Take a trial step along the search direction
             call sparse_embed_array_copy(trial_aux%kern, aux%kern)
             call sparse_embed_array_axpy(trial_aux%kern, &
                  con_direction, trial_step)
             call kernel_workspace_invalidate(trial_aux)
          end if

          ! ep: %%%%%%%%%%%%%%%%%%%%% NEW LINE SEARCH (ADAPTED FROM EDFT) %%%%%%%%%%%%%%

          ! ep: STEP 1
          total_energy_at_trial_length = internal_energy(trial_aux,.false.)
          if (pub_output_detail.ge.VERBOSE) &
               call edft_print_line_search('   Trial1', trial_step, total_energy_at_trial_length)
          npts = 2
          xx2 = (/ 0.0_DP, trial_step /); yy2 = (/ total_energy, total_energy_at_trial_length /)
          trial_step2 = mermin_trial_step(npts, xx2, yy2)
          !trial_step2 =trial_step*2.0_DP

          ! ep: STEP 2
          call mermin_evaluate_step(aux, trial_aux, con_direction,trial_step2)
          total_energy_at_trial_length2 = internal_energy(trial_aux,.false.)
          if (pub_output_detail.ge.VERBOSE) then
             call edft_print_line_search('   Trial2', trial_step2, total_energy_at_trial_length2)
          end if
          call edft_fit_second_order_poly_3points(quadratic_success, trial_step3, predicted_e3, 0.0_DP, trial_step, &
               trial_step2, total_energy, total_energy_at_trial_length, total_energy_at_trial_length2)
          ! validate quadratic step
          if ( (.not.quadratic_success).or.(trial_step3.lt.0.0_DP).or.&
                (trial_step3.gt.pub_mermin_cg_max_step) ) then
             quadratic_success = .false.
             ! ars: take trial step
             npts=3; xx3=(/ 0.0_DP, trial_step, trial_step2 /); yy3=(/ total_energy, total_energy_at_trial_length, &
             total_energy_at_trial_length2/)
             trial_step3=mermin_trial_step(npts, xx3, yy3)
             !trial_step3 =trial_step*3.0_DP
          end if

          ! ep: STEP 3
          call mermin_evaluate_step(aux, trial_aux, con_direction,trial_step3)
          total_energy_at_trial_length3 = internal_energy(trial_aux,.false.)

          if ( (quadratic_success).and.(total_energy_at_trial_length3.lt.total_energy).and.&
               (total_energy_at_trial_length3.lt.total_energy_at_trial_length).and.&
               (total_energy_at_trial_length3.lt.total_energy_at_trial_length2) ) then
             quadratic_accepted = .true.
             accepted_e = total_energy_at_trial_length3; accepted_step=trial_step3
          else
             quadratic_accepted = .false.
          end if

          if (pub_output_detail.ge.VERBOSE) then
             if (quadratic_success) then
                call edft_print_line_search('Quadratic', trial_step3, total_energy_at_trial_length3, &
                     predicted_e3)
             else
                call edft_print_line_search('   Trial3', trial_step3, total_energy_at_trial_length3)
             end if
          end if


          if (.not. quadratic_accepted) then
             !! CUBIC STEP !!

             call edft_fit_third_order_poly_4points(cubic_success, &
                  trial_step4, predicted_e4, 0.0_DP, trial_step, trial_step2, trial_step3, &
                  total_energy, total_energy_at_trial_length, total_energy_at_trial_length2, &
                  total_energy_at_trial_length3)
             ! validate cubic step
             if ( (.not.cubic_success).or.(trial_step4.lt.0.0_DP).or.&
                  (trial_step4.gt.pub_mermin_cg_max_step) ) then
                cubic_success = .false.
                ! take trial step
                npts=4;
                xx4=(/ 0.0_DP, trial_step, trial_step2, trial_step3 /);
                yy4=(/ total_energy, total_energy_at_trial_length, total_energy_at_trial_length2, &
                    total_energy_at_trial_length3 /)
                trial_step4=mermin_trial_step(npts, xx4, yy4)
                !trial_step4 =trial_step*4.0_DP
             end if

             !ep: STEP 4
             call mermin_evaluate_step(aux, trial_aux, con_direction,trial_step4)
             total_energy_at_trial_length4 = internal_energy(trial_aux,.false.)

             ! Validate cubic_step
             if ( (cubic_success).and.(total_energy_at_trial_length4.lt.total_energy).and.&
                (total_energy_at_trial_length4.lt.total_energy_at_trial_length).and.&
                (total_energy_at_trial_length4.lt.total_energy_at_trial_length2).and.&
                (total_energy_at_trial_length4.lt.total_energy_at_trial_length3) ) then
                cubic_accepted = .true.
                accepted_e = total_energy_at_trial_length4; accepted_step=trial_step4
             else
                cubic_accepted = .false.
             end if

             ! print line search if verbose
             if (pub_output_detail.ge.VERBOSE) then
                if (cubic_success) then
                   call edft_print_line_search('    Cubic', trial_step4, total_energy_at_trial_length4, &
                        predicted_e4)
                else
                   call edft_print_line_search('   Trial4', trial_step4, total_energy_at_trial_length4)
                end if
             end if


             if (.not. cubic_accepted) then
                !! QUARTIC STEP !!

                ! set quartic step
                call edft_fit_fourth_order_poly_5points(quartic_success, &
                     trial_step5, predicted_e5, 0.0_DP, trial_step, trial_step2, &
                     trial_step3, trial_step4, total_energy, total_energy_at_trial_length, &
                     total_energy_at_trial_length2, total_energy_at_trial_length3, &
                     total_energy_at_trial_length4)
                ! validate quartic step
                if ( (.not.quartic_success).or.(trial_step5.lt.0.0_DP).or.&
                     (trial_step5.gt.pub_mermin_cg_max_step) ) then
                   quartic_success = .false.
                   trial_step5 = 0.01_DP*min(trial_step,trial_step2,trial_step3,trial_step4)
                end if

                ! ep: STEP 5
                call mermin_evaluate_step(aux, trial_aux, con_direction,trial_step5)
                total_energy_at_trial_length5= internal_energy(trial_aux,.false.)

                ! validate quartic_step
                if ( (quartic_success).and.(total_energy_at_trial_length5.lt.total_energy).and.&
                   (total_energy_at_trial_length5.lt.total_energy_at_trial_length).and.&
                   (total_energy_at_trial_length5.lt.total_energy_at_trial_length2).and.&
                   (total_energy_at_trial_length5.lt.total_energy_at_trial_length3).and.&
                   (total_energy_at_trial_length5.lt.total_energy_at_trial_length4) ) then
                   quartic_accepted = .true.
                   accepted_e = total_energy_at_trial_length5; accepted_step=trial_step5
                else
                   quartic_accepted = .false.
                end if

                ! print line search if verbose
                if (pub_output_detail.ge.VERBOSE) then
                   if (quartic_success) then
                      call edft_print_line_search('  Quartic',trial_step5,total_energy_at_trial_length5,&
                           predicted_e5)
                   else
                      call edft_print_line_search('   Trial5',trial_step5,total_energy_at_trial_length5)
                   end if
                end if

                tiny_accepted = .false.
                if ( (.not.quartic_accepted).and.(quartic_success).and.&
                     (total_energy.lt.total_energy_at_trial_length).and.&
                     (total_energy.lt.total_energy_at_trial_length2).and.&
                     (total_energy.lt.total_energy_at_trial_length3).and.&
                     (total_energy.lt.total_energy_at_trial_length4).and.&
                     (total_energy.lt.total_energy_at_trial_length5) ) then
                   !! TINY STEP !!

                   ! ep : STEP 6
                   ! set step to a tiny value
                   trial_step6 = 0.01_DP*min(trial_step,trial_step2,trial_step3,trial_step4,&
                                 trial_step5)

                   call mermin_evaluate_step(aux, trial_aux, con_direction,trial_step6)
                   total_energy_at_trial_length6= internal_energy(trial_aux,.false.)

                   ! validate quartic_step
                   if ( (total_energy_at_trial_length6.lt.total_energy).and.&
                        (total_energy_at_trial_length6.lt.total_energy_at_trial_length).and.&
                        (total_energy_at_trial_length6.lt.total_energy_at_trial_length2).and.&
                        (total_energy_at_trial_length6.lt.total_energy_at_trial_length3).and.&
                        (total_energy_at_trial_length6.lt.total_energy_at_trial_length4).and.&
                        (total_energy_at_trial_length6.lt.total_energy_at_trial_length5) ) then
                      tiny_accepted = .true.
                      accepted_e =total_energy_at_trial_length6; accepted_step=trial_step6
                   else
                      tiny_accepted = .false.
                   end if

                   ! print line search if verbose
                   if (pub_output_detail.ge.VERBOSE) then
                      call edft_print_line_search('   Trial6', trial_step6, total_energy_at_trial_length6)
                   end if

                end if ! .not.quartic_accepted etc
             end if ! .not.cubic_accepted
          end if ! .not.quadratic accepted

          ! Choose best of trial steps if none accepted
          if ( (.not.tiny_accepted).and.(.not.quartic_accepted).and.&
               (.not.cubic_accepted).and.(.not.quadratic_accepted) ) then

             ! set best step
             if ((total_energy_at_trial_length6.lt.total_energy_at_trial_length).and.&
                (total_energy_at_trial_length6.lt.total_energy_at_trial_length2).and.&
                (total_energy_at_trial_length6.lt.total_energy_at_trial_length3).and.&
                (total_energy_at_trial_length6.lt.total_energy_at_trial_length4).and.&
                (total_energy_at_trial_length6.lt.total_energy_at_trial_length5)) then
                accepted_step = trial_step6
             else if ( (total_energy_at_trial_length5.lt.total_energy_at_trial_length).and.&
                (total_energy_at_trial_length5.lt.total_energy_at_trial_length2).and.&
                (total_energy_at_trial_length5.lt.total_energy_at_trial_length3).and.&
                (total_energy_at_trial_length5.lt.total_energy_at_trial_length4).and.&
                (total_energy_at_trial_length5.lt.total_energy_at_trial_length6) ) then
                accepted_step = trial_step5
             else if ( (total_energy_at_trial_length4.lt.total_energy_at_trial_length).and.&
                (total_energy_at_trial_length4.lt.total_energy_at_trial_length2).and.&
                (total_energy_at_trial_length4.lt.total_energy_at_trial_length3).and.&
                (total_energy_at_trial_length4.lt.total_energy_at_trial_length5).and.&
                (total_energy_at_trial_length4.lt.total_energy_at_trial_length6) ) then
                accepted_step = trial_step4
             else if ( (total_energy_at_trial_length3.lt.total_energy_at_trial_length).and.&
                (total_energy_at_trial_length3.lt.total_energy_at_trial_length2).and.&
                (total_energy_at_trial_length3.lt.total_energy_at_trial_length4).and.&
                (total_energy_at_trial_length3.lt.total_energy_at_trial_length5).and.&
                (total_energy_at_trial_length3.lt.total_energy_at_trial_length6) ) then
                accepted_step = trial_step3
             else if ( (total_energy_at_trial_length2.lt.total_energy_at_trial_length).and.&
                (total_energy_at_trial_length2.lt.total_energy_at_trial_length3).and.&
                (total_energy_at_trial_length2.lt.total_energy_at_trial_length4).and.&
                (total_energy_at_trial_length2.lt.total_energy_at_trial_length5).and.&
                (total_energy_at_trial_length2.lt.total_energy_at_trial_length6) ) then
                accepted_step = trial_step2
             else if ( (total_energy_at_trial_length.lt.total_energy_at_trial_length2).and.&
                (total_energy_at_trial_length.lt.total_energy_at_trial_length3).and.&
                (total_energy_at_trial_length.lt.total_energy_at_trial_length4).and.&
                (total_energy_at_trial_length.lt.total_energy_at_trial_length5).and.&
                (total_energy_at_trial_length.lt.total_energy_at_trial_length6) ) then
                accepted_step = trial_step
             end if

             if (pub_output_detail.ge.VERBOSE) &
                  call edft_print_line_search(' NoMinFnd', accepted_step, accepted_e)

          end if ! no steps accepted

          ! ep: guess for next step to be thoroughly tested before using
          !if (iteration .gt. 1) then
          !
          !   call sparse_embed_array_trace(trace_array, old_con_direction, old_co_gradient, &
          !           opA=opB_loc)
          !   old_d_dot_g = sum(trace_array)
          !
          !   call sparse_embed_array_trace(trace_array, con_direction, co_gradient, &
          !           opA=opB_loc)
          !   new_d_dot_g = sum(trace_array)
          !
          !    trial_step= line_search_step * (old_d_dot_g/new_d_dot_g)
          !else
          !   trial_step=accepted_step
          !end if
          ! ep: guess for next step to be thoroughly tested before using

          ! ep: QC_TEST
          if (pub_print_qc) call mermin_print_qc(mu, iteration, rms_gradient, &
              aux, con_gradient, line_search_step, Qslope, Eslope, DFTslope, &
              Muslope, total_energy, rep, fd_slope_op, fd_slope_ent, fd_slope_dft, &
              fd_slope_mu )
          ! ep: QC_TEST

          ! Update the auxiliary matrix
          if (pub_mermin_check_trial_steps) then
             line_search_step = accepted_step
             call internal_kernel_protected_step(trial_aux, line_search_step, &
                  aux,con_direction)
             call sparse_embed_array_copy(aux%kern, trial_aux%kern)
             trial_step = accepted_step
          else
             line_search_step = accepted_step
             call sparse_embed_array_axpy(aux%kern, con_direction, line_search_step)
             trial_step = accepted_step
          end if

          call kernel_workspace_invalidate(aux)

       end if

       call internal_print

    end do

1  format(4x,I4.1,T25,' <--   MERMIN_STEPS')

    ! Report convergence failure
    if (pub_on_root .and. (.not. converged) .and. pub_output_detail >= NORMAL) then

       write(stdout,'(a,i4,a)') &
           'Finished density kernel iterations (',current_maxit_mermin, ')'
    endif

    ! Write density kernel to file if requested
    loc_force_no_IO=.false.
    if(present(force_no_IO)) then
       loc_force_no_IO=force_no_IO
    end if

    if(.not.loc_force_no_IO) then
       ! jcap: if we are doing an EMFT follow calculation AND we are
       ! doing the final mermin part of the calculation, change the name of
       ! the dkn file
       if (pub_write_denskern.and.(.not.pub_write_converged_dk_ngwfs)) &
            call restart_kernel_write(aux%kern, write_cond=pub_cond_calculate, &
            write_emft=(pub_emft_follow.and.pub_emft.and.pub_emft_lnv_only))
    end if
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! ** De-Allocate structures for sparse matrices

    ! Destroy CG workspace
    if(pub_mermin_cg_type/='MERMIN_FLETCHER')then
       call sparse_embed_array_destroy(old_co_gradient)
    endif
    call sparse_embed_array_destroy(con_direction)
    call sparse_embed_array_destroy(con_gradient)
    call sparse_embed_array_destroy(co_gradient)
    call sparse_embed_array_destroy(old_con_direction)

    ! ep : destroy QC variable
    if( pub_print_qc) then
        call sparse_embed_array_destroy(entr_grad_co)
        call sparse_embed_array_destroy(entr_grad_comp)
        call sparse_embed_array_destroy(entr_dir_comp)
        call sparse_embed_array_destroy(old_entr_dir_comp)
        call sparse_embed_array_destroy(dft_grad_co)
        call sparse_embed_array_destroy(dft_grad_comp)
        call sparse_embed_array_destroy(dft_dir_comp)
        call sparse_embed_array_destroy(old_dft_dir_comp)
        call sparse_embed_array_destroy(mu_grad_co)
        call sparse_embed_array_destroy(mu_grad_comp)
        call sparse_embed_array_destroy(mu_dir_comp)
        call sparse_embed_array_destroy(old_mu_dir_comp)
    endif
    ! ep : destroy QC variable

    call kernel_destroy(trial_aux)

    ! Deallocate ks and ksk workspace
    call kernel_workspace_destroy(aux)

    deallocate(trace_array,stat=ierr)
    call utils_dealloc_check('mermin_denskernel_optimise_cg','trace_array',ierr)

    ! Stop the timer
    call timer_clock('mermin_denskernel_optimise_cg', 2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
         &Leaving mermin_denskernel_optimise_cg'

    ! Flush output
    call services_flush

    return

  contains

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    real(kind=DP) function internal_energy(kernel_var,update_ham,spoil_force)

      use constants, only: paw_en_size, k_B
      use hamiltonian, only: hamiltonian_dens_dep_matrices
      use kernel, only: DKERN, kernel_workspace_invalidate
      use rundat, only: pub_spin_fac, pub_mermin_cheb
      use sparse_embed, only: sparse_embed_trace, &
         sparse_embed_array_create, sparse_embed_array_scale, &
         sparse_embed_array_product, sparse_embed_array_destroy, &
         SPAM3_EMBED_ARRAY
      use smearing_operator, only: calculate_fermi_entropy_mermin

      implicit none

      ! Argument
      type(DKERN), intent(inout) :: kernel_var
      logical, intent(in) :: update_ham
      logical,optional, intent(in) :: spoil_force

      ! Local variables
      integer :: is, ik
      type(SPAM3_EMBED_ARRAY) :: app
      real(kind=DP) :: smermin(pub_num_spins)
      real(kind=DP), &
           dimension(aux%kern%num_spins, aux%kern%num_kpoints) :: trace
      real(kind=DP) :: paw_sphere_energies(paw_en_size)

      ! Stored versions of ks and ksk no longer valid (except on first call)
      if (.not.update_ham) call kernel_workspace_invalidate(kernel_var)

      norm_fac(:,:) = 1.0_DP
      smermin(:) = 0.0_DP

      call hamiltonian_dens_dep_matrices(ham, lhxc_fine, internal_energy, &
           lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
           ngwf_basis, hub_proj_basis, hub, kernel_var%kern, mdl,&
           hfxstate, update_ham, lhxc_fixed=.false.,spoil_force=spoil_force,&
           dfdtau_fine = dfdtau_fine, kpt=loc_kpt)
      call sparse_embed_array_create(app, kernel_var%kern)
      call sparse_embed_array_product(app, kernel_var%kern, rep%overlap)

      do ik = 1, aux%kern%num_kpoints
         do is = 1, pub_num_spins
            !trace(is, ik) = sparse_embed_trace(app%m(is,ik))
            call sparse_embed_trace(trace(is, ik), app%m(is,ik))
            call calculate_fermi_entropy_mermin(kernel_var%kern%m(is,ik), rep%overlap,1,smermin(is), pub_mermin_cheb)
         enddo
      enddo

      do ik = 1, aux%kern%num_kpoints
         do is = 1, pub_num_spins
            if (pub_num_spins .eq. 2) then
                internal_energy = internal_energy + (smermin(is)*(pub_mermin_smearing_width/k_B)) &
                                  + mu(is,ik) * (trace(is,PUB_1K)-rep%n_occ(is,ik))
            else
                internal_energy = internal_energy + (smermin(is)*pub_spin_fac*(pub_mermin_smearing_width/k_B)) &
                                  + mu(is,ik) * ((sum(trace(:,ik))*pub_spin_fac) &
                                  -(sum(rep%n_occ(:,ik))*pub_spin_fac))
            end if
         enddo
      enddo

      call sparse_embed_array_destroy(app)

    end function internal_energy

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    real(kind=DP) function internal_ent(kernel_var,update_ham,spoil_force)

      use constants, only: paw_en_size, k_B
      use kernel, only: DKERN, kernel_workspace_invalidate
      use rundat, only: pub_spin_fac
      use sparse_embed, only: sparse_embed_trace, &
         sparse_embed_array_create, sparse_embed_array_product, &
         sparse_embed_array_destroy, SPAM3_EMBED_ARRAY
      use smearing_operator, only: calculate_fermi_entropy_mermin

      implicit none

      ! Argument
      type(DKERN), intent(inout) :: kernel_var
      logical, intent(in) :: update_ham
      logical,optional, intent(in) :: spoil_force

      ! Local variables
      integer :: is, ik
      type(SPAM3_EMBED_ARRAY) :: app
      real(kind=DP) :: smermin(pub_num_spins)
      real(kind=DP), &
           dimension(aux%kern%num_spins, aux%kern%num_kpoints) :: trace
      real(kind=DP) :: paw_sphere_energies(paw_en_size)

      ! Stored versions of ks and ksk no longer valid (except on first call)
      if (.not.update_ham) call kernel_workspace_invalidate(kernel_var)

      norm_fac(:,:) = 1.0_DP
      smermin(:) = 0.0_DP
      internal_ent=0.0_DP

      call sparse_embed_array_create(app, kernel_var%kern)
      call sparse_embed_array_product(app, kernel_var%kern, rep%overlap)

      do ik = 1, aux%kern%num_kpoints
         do is = 1, pub_num_spins
            !trace(is, ik) = sparse_embed_trace(app%m(is,ik))
            call sparse_embed_trace(trace(is, ik),app%m(is,ik))
            call calculate_fermi_entropy_mermin(kernel_var%kern%m(is,ik), rep%overlap,1,smermin(is), pub_mermin_cheb)
         enddo
      enddo

      do ik = 1, aux%kern%num_kpoints
         do is = 1, pub_num_spins
            if (pub_num_spins .eq. 2) then
                internal_ent = internal_ent + (smermin(is)*(pub_mermin_smearing_width/k_B))
            else
                internal_ent = internal_ent + (smermin(is)*pub_spin_fac*(pub_mermin_smearing_width/k_B))
            end if
         enddo
      enddo

      call sparse_embed_array_destroy(app)

    end function internal_ent

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    ! ep : DFT QC internal energy contribution
    real(kind=DP) function internal_dft(kernel_var,update_ham,spoil_force)

      use constants, only: paw_en_size
      use hamiltonian, only: hamiltonian_dens_dep_matrices
      use kernel, only: DKERN, kernel_workspace_invalidate
      use sparse_embed, only: SPAM3_EMBED_ARRAY

      implicit none

      ! Argument
      type(DKERN), intent(inout) :: kernel_var
      logical, intent(in) :: update_ham
      logical,optional, intent(in) :: spoil_force

      ! Local variables
      integer :: is, ik
      type(SPAM3_EMBED_ARRAY) :: app
      real(kind=DP) :: smermin(pub_num_spins)
      real(kind=DP) :: app_energy
      real(kind=DP), &
           dimension(aux%kern%num_spins, aux%kern%num_kpoints) :: trace
      real(kind=DP) :: paw_sphere_energies(paw_en_size)

      ! Stored versions of ks and ksk no longer valid (except on first call)
      if (.not.update_ham) call kernel_workspace_invalidate(kernel_var)

      norm_fac(:,:) = 1.0_DP
      smermin(:) = 0.0_DP

      call hamiltonian_dens_dep_matrices(ham, lhxc_fine, app_energy, &
           lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
           ngwf_basis, hub_proj_basis, hub, kernel_var%kern, mdl,&
           hfxstate, update_ham, lhxc_fixed=.false.,spoil_force=spoil_force,&
           dfdtau_fine = dfdtau_fine, kpt=loc_kpt)
      internal_dft = app_energy

    end function internal_dft

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    real(kind=DP) function internal_muen(kernel_var,update_ham,spoil_force)

      use constants, only: paw_en_size
      use kernel, only: DKERN, kernel_workspace_invalidate
      use rundat, only: pub_spin_fac, pub_xc_ke_density_required
      use sparse_embed, only: sparse_embed_trace, &
         sparse_embed_array_create, sparse_embed_array_product, &
         sparse_embed_array_destroy, SPAM3_EMBED_ARRAY, &
         sparse_embed_array_trace
      implicit none

      ! Argument
      type(DKERN), intent(inout) :: kernel_var
      logical, intent(in) :: update_ham
      logical,optional, intent(in) :: spoil_force

      ! Local variables
      integer :: is, ik
      type(SPAM3_EMBED_ARRAY) :: app
      real(kind=DP), &
           dimension(aux%kern%num_spins, aux%kern%num_kpoints) :: trace
      real(kind=DP) :: paw_sphere_energies(paw_en_size)

      ! Stored versions of ks and ksk no longer valid (except on first call)
      if (.not.update_ham) call kernel_workspace_invalidate(kernel_var)

      call sparse_embed_array_create(app, kernel_var%kern)
      call sparse_embed_array_product(app, kernel_var%kern, rep%overlap)

      do ik = 1, aux%kern%num_kpoints
         do is = 1, pub_num_spins
             call sparse_embed_trace(trace(is, ik), app%m(is,ik))
         enddo
      enddo

      do ik = 1, aux%kern%num_kpoints
         do is = 1, pub_num_spins
            if (pub_num_spins .eq. 2) then
                internal_muen =  mu(is,ik) * (trace(is,PUB_1K)-rep%n_occ(is,ik))
            else
                internal_muen =  mu(is,ik) * ((sum(trace(:,ik))*pub_spin_fac) &
                                  -(sum(rep%n_occ(:,ik))*pub_spin_fac))
            end if
         enddo
      enddo

      call sparse_embed_array_destroy(app)

    end function internal_muen

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    real(kind=DP) function internal_gradient_norm()

      use constants, only: SAFE_DIV_EPS
      use kernel, only: kernel_commutator
      use rundat, only: pub_rms_kernel_measure
      use sparse_embed, only: sparse_embed_product, sparse_embed_rms_element, &
          sparse_embed_create, sparse_embed_num_element, SPAM3_EMBED

      implicit none

      ! Local variables
      integer :: is
      real(kind=DP) :: rms_el
      real(kind=DP) :: mu_extra(pub_num_spins, pub_num_kpoints) ! Additional mu
                                                    ! resulting from truncation
      type(SPAM3_EMBED) :: mixed_old_con_dir
      ! agrecocmplx
      character(len=1) :: loc_opA

      ! Start Timer
      call timer_clock('mermin_gradient_norm',1)

      ! Calculate covariant gradient (no tensor correction)
      call internal_co_gradient
      ! Apply tensor correction to obtain contravariant gradient
      call internal_con_gradient

      ! Calculate inner product between gradients to obtain RMS gradient
      ! KPOINTS_DANGER: the following sum also sums over k-points
      ! agrecocmplx: need to take conjugate of con_gradient here?
      ! I think so since this is essentially a dot product
      if (loc_cmplx) then
         loc_opA = 'C'
      else
         loc_opA = 'N'
      end if

      call sparse_embed_array_trace(trace_array, con_gradient, co_gradient, &
              opA=loc_opA)
      new_g_dot_g = sum(trace_array)

      if (new_g_dot_g < 0.0_DP .and. pub_on_root) write(stdout,'(a)') &
              'WARNING in internal_gradient_norm (mermin_mod.F90): &
              &inverse overlap is not positive definite!'

      internal_gradient_norm = sqrt(abs(new_g_dot_g) / &
           sparse_embed_num_element(con_gradient%m(1,1)))

      ! Calculate RMS change in the density kernel from previous iteration
      ! and RMS commutator
      rms_dkern = 0.0_DP
      do is=1,pub_num_spins
         if (pub_rms_kernel_measure) then
            rms_el = sparse_embed_rms_element(old_con_direction%m(is,PUB_1K))
            rms_dkern = rms_dkern + rms_el * rms_el
         else
            ! ddor: Change to invariant, per-electron version
            call sparse_embed_create(mixed_old_con_dir,old_con_direction%m(is,PUB_1K),&
                 rep%overlap)
            call sparse_embed_product(mixed_old_con_dir,old_con_direction%m(is,PUB_1K),&
                 rep%overlap)
            ! ddor: Allow for case of empty density kernel for one spin
            if (abs(real(rep%n_occ(is,PUB_1K),kind=DP)) > SAFE_DIV_EPS) then
               ! agrecocmplx: need to take conjugate of mixed_old_con_dir here?
               ! I think so since this is essentially a dot product
               call sparse_embed_trace(rms_el, mixed_old_con_dir, mixed_old_con_dir, &
                    opA=loc_opA)
               rms_el = rms_el / (real(rep%n_occ(is,PUB_1K),kind=DP))**2.0_DP
            else
               rms_el = 0.0_DP
            endif
            call sparse_embed_destroy(mixed_old_con_dir)
            rms_dkern = rms_dkern + rms_el
         endif
      end do

      rms_dkern = sqrt(rms_dkern / pub_num_spins) * line_search_step

      ! rab207: moved kernel choice into kernel_mod.F90
      commutator = kernel_commutator(aux, ham%ham, rep%overlap, &
                                     rep%inv_overlap)

      ! Stop Timer
      call timer_clock('mermin_gradient_norm',2)

    end function internal_gradient_norm

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_co_gradient

      !========================================================================!
      ! This subroutine calculates the gradient of the Free Energy             !
      ! wrt to the density kernel.                                             !
      ! In this module the Free Energy gradient is defined as:                 !
      ! dA/dK^cd= <phi_d|'T'+V_ext|phi_c> + V_dc + T*{(S_**D(K^**S_**)} +      !
      !           + mu*S_dc                                                    !
      ! where:                                                                 !
      ! dA/dK^ab = Derivative of the Free energy wrt the Density Kernel        !
      ! phi_c, phi_d = NGWFs                                                   !
      ! K^cd = Density Kernel                                                  !
      ! 'T'= Kinetic operator                                                  !
      ! V_ext = External potential operator                                    !
      ! V_LHXC = LHXC potential in NGWF basis                                  !
      ! T = Temperature                                                        !
      ! D(K^**S_**) = Chebyshev expansion of the derivative of the Entropy     !
      ! mu = Chemical Potential (Defined by Millam Scuseria method)            !
      ! S_dc, S_** = overlap matrix                                            !
      !------------------------------------------------------------------------!
      ! Written by Poli Emiliaino                                              !
      ! Based on an original version by Chris-Kriton Skylaris                  !
      !========================================================================!

      use kernel, only: kernel_validate_ks
      !ep
      use constants, only: k_B
      use smearing_operator, only: calculate_fermi_entropy_mermin
      use rundat, only: pub_spin_fac, pub_mermin_smearing_width, pub_mermin_cheb, &
          PUB_1K, pub_print_qc
      use sparse_embed, only: sparse_embed_array_scale,  sparse_embed_scale, &
          sparse_embed_array_copy, SPAM3_EMBED
      !ep

      implicit none

      ! Local variables
      integer :: is
      integer       :: ierr
      real(kind=DP) :: smermin
      type(SPAM3_EMBED), allocatable ::  hamtmp(:)
      type(SPAM3_EMBED_ARRAY) :: deriv, derivok
      real(kind=DP), dimension(aux%kern%num_spins) :: trace


      call kernel_validate_ks(aux, rep%overlap)
      call sparse_embed_array_create(deriv, aux%kern)
      call sparse_embed_array_create(derivok, aux%kern)

      !ep: copy Hamiltonian
      allocate(hamtmp(1:pub_num_spins), stat=ierr)
      do is = 1, pub_num_spins
         call sparse_embed_create(hamtmp(is), aux%kern%m(is,PUB_1K))
         call sparse_embed_copy(hamtmp(is), ham%ham(is))
      end do

      !ep: copy Hamiltonian for QC components of the gradient
      if (pub_print_qc) then
         do is = 1, pub_num_spins
            call sparse_embed_copy(dft_grad_co%m(is,PUB_1K), ham%ham(is))
         end do
      end if

      !ep: CALCULATE dS
      do is = 1, pub_num_spins
         call calculate_fermi_entropy_mermin(aux%kern%m(is,PUB_1K), &
         & rep%overlap,2,smermin,pub_mermin_cheb,deriv%m(is,PUB_1K))
      enddo
      !ep: CALCULATE S * dS
      call sparse_embed_array_product(derivok, rep%overlap, deriv)
      !ep : CALCULATE T * S * dS
      call sparse_embed_array_scale(derivok, (pub_mermin_smearing_width/k_B), 0.0_DP)
      !ep: copy T * S * dS for QC components of the gradient
      if (pub_print_qc) then
         call sparse_embed_array_copy(entr_grad_co, derivok)
      end if
      do is = 1, pub_num_spins
         !ep : CALCULATE H + T * S * dS
         call sparse_embed_axpy(hamtmp(is), derivok%m(is,PUB_1K), 1.0_DP)
         !ep: CALCULATE (H + T * dS * S) + mu * S
         call sparse_embed_axpy(hamtmp(is), rep%overlap, mu(is,PUB_1K))
      end do
      !ep: copy mu * S for QC components of the gradient
      if (pub_print_qc) then
         do is = 1, pub_num_spins
         call sparse_embed_scale(mu_grad_co%m(is,PUB_1K), 0.0_DP)
         call sparse_embed_axpy(mu_grad_co%m(is,PUB_1K), rep%overlap, mu(is,PUB_1K))
         enddo
      end if
      do is = 1, pub_num_spins
         call sparse_embed_copy(co_gradient%m(is,PUB_1K), hamtmp(is))
      end do
      call sparse_embed_array_scale(co_gradient, pub_spin_fac)
      if (pub_print_qc) then
         call sparse_embed_array_scale(mu_grad_co, pub_spin_fac)
         call sparse_embed_array_scale(entr_grad_co, pub_spin_fac)
         call sparse_embed_array_scale(dft_grad_co, pub_spin_fac)

      end if

      !print (spin-resolved) trace of covariant gradient for testing
      !do is = 1, pub_num_spins
      !   trace(is) = sparse_embed_trace(co_gradient%m(is,PUB_1K))
      !end do

      call sparse_embed_array_destroy(deriv)
      call sparse_embed_array_destroy(derivok)
      do is = 1, pub_num_spins
         call sparse_embed_destroy(hamtmp(is))
      end do


    end subroutine internal_co_gradient

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_con_gradient

      !===================================================================!
      ! This subroutine converts a covariant gradient to a contravariant  !
      ! gradient by pre- and post-multiplying by the inverse overlap      !
      ! matrix. Note that truncation will introduce errors at this stage. !
      !-------------------------------------------------------------------!
      ! Written by Emiliano Poli                                          !
      ! Based on an original version by Chris-Kriton Skylaris             !
      !===================================================================!

      implicit none

      ! Local variables
      type(SPAM3_EMBED_ARRAY) :: temp, temp2
      type(SPAM3_EMBED_ARRAY) :: tempS, tempDFT, tempMU
      !ep: debug trace
      !real(kind=DP), &
      !     dimension(aux%kern%num_spins) :: trace


      !ep: Allocate temporary sparse matrix
      call sparse_embed_array_create(temp, aux%kern )
      !call sparse_embed_array_create(temp2, aux%kern )
      !ep: Allocate temporary sparse matrix for QC control
      if (pub_print_qc) then
         call sparse_embed_array_create(tempS, aux%kern )
         call sparse_embed_array_create(tempDFT, aux%kern )
         call sparse_embed_array_create(tempMU, aux%kern )
      end if

      !ep: Copy co_grad in con_grad
      call sparse_embed_array_copy(con_gradient, co_gradient)

      !ep: Copy co_grad in con_grad  components for QC
      if (pub_print_qc) then
         call sparse_embed_array_copy(entr_grad_comp,entr_grad_co)
         call sparse_embed_array_copy(mu_grad_comp,mu_grad_co)
         call sparse_embed_array_copy(dft_grad_comp, dft_grad_co)
      end if
      !ep: Calculate temp := con_grad . S^-1
      call sparse_embed_array_product(temp, con_gradient, rep%inv_overlap)

      !ep: Calculate S^-1 . (co_grad + mu_extra S) . S^-1 =
      !                        = S^-1 . temp -> con_grad
      call sparse_embed_array_product(con_gradient, rep%inv_overlap, temp)

      !ep: Calculate temp := con_grad . S^-1
      !ep: Calculate S^-1 . (co_grad + mu_extra S) . S^-1 =
      != S^-1 . temp -> con_grad for QC components
      if (pub_print_qc) then
         call sparse_embed_array_product(tempS, entr_grad_comp, rep%inv_overlap)
         call sparse_embed_array_product(entr_grad_comp, rep%inv_overlap, tempS)
         call sparse_embed_array_product(tempDFT, dft_grad_comp, rep%inv_overlap)
         call sparse_embed_array_product(dft_grad_comp, rep%inv_overlap,tempDFT)
         call sparse_embed_array_product(tempMU, mu_grad_comp, rep%inv_overlap)
         call sparse_embed_array_product(mu_grad_comp, rep%inv_overlap, tempMU)
      endif

      !ep: DEBUG TRACE
      !ep: Calculate G^ S_
      !call sparse_embed_array_product(temp2, con_gradient, rep%overlap)

      !ep: DEBUG PRINT
      !ep: Calculate Tr( G^ S_ )
      !do is = 1, pub_num_spins
      !   trace(is) = sparse_embed_trace(temp2%m(is,PUB_1K))
      !end do

      !ep: DEBUG PRINT
      !if (pub_debug_on_root) then
         !if (pub_on_root) then
         !   write(stdout,*) "Trace of G^ S_:", trace(1), trace(2)
         !endif
      !endif

      ! Deallocate workspace
      call sparse_embed_array_destroy(temp)
      if (pub_print_qc) then
         call sparse_embed_array_destroy(tempS)
         call sparse_embed_array_destroy(tempDFT)
         call sparse_embed_array_destroy(tempMU)
      end if
      !ep:  DEBUG PRINT RELATED
      !call sparse_embed_array_destroy(temp2)
      !ep:  DEBUG PRINT RELATED

    end subroutine internal_con_gradient

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    logical function internal_test_convergence()

      implicit none

      ! Convergence is assumed when either the RMS gradient or RMS commutator
      ! falls below the threshold

      if ((rms_gradient < mermin_threshold .and. iteration > 0) .or. &
           (commutator < mermin_threshold .and. iteration > 1)) then

         if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then

            if (iteration == 1) then
               write(stdout,'(i5,es27.12,2es15.3)') iteration, &
                    total_energy, commutator, rms_gradient
            else
               write(stdout,'(i5,es27.12,3es15.3)') iteration, total_energy, &
                    rms_gradient, commutator, total_energy - old_total_energy
            endif

            write(stdout,'(/2a)') repeat(' ',17), repeat('.',46)
            write(stdout,'(a,f19.14,a)')'                 | RMS MERMIN GRADIENT&
                 &= ', rms_gradient, '      |'
            write(stdout,'(a)') '                 | MERMIN density kernel &
                 &optimisation converged! |'
            write(stdout,'(2a)') repeat(' ',17), repeat('~',46)

         end if

         internal_test_convergence = .true.
         conv_status = 0

      else

         internal_test_convergence = .false.
         conv_status = 1

      end if

    end function internal_test_convergence

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_search_direction

      use rundat, only: pub_spin_fac

      implicit none

      ! Local variables
      integer :: is
      real(kind=DP) :: polak
      real(kind=DP) :: dltne

      ! Check previous scalar product of gradients is not vanishing

      if (abs(old_g_dot_g) < tiny(1.0_DP)) then

         if (pub_on_root) write(stdout,'(a/a)') &
              'WARNING in mermin_denskernel_optimise_cg: &
              &vanishing previous gradient norm -', &
              '  conjugate gradients coefficient reset'

         cg_coeff = 0.0_DP

      else

         ! Calculate correction factor for Polak-Ribiere formula if necessary
         if (pub_mermin_cg_type == 'MERMIN_FLETCHER') then
            polak = 0.0_DP
         else if (pub_mermin_cg_type == 'MERMIN_POLAK') then
            ! KPOINTS_DANGER
            call sparse_embed_array_trace(trace_array, con_gradient, old_co_gradient)
            polak = sum(trace_array)
         else
            if (pub_on_root) write(stdout,'(a/3a)') &
              'WARNING in mermin_denskernel_optimise_cg: &
              &unknown conjugate gradients method','  "',trim(pub_mermin_cg_type), &
              '" - switching to steepest descents'
            polak = new_g_dot_g
         end if

         ! Calculate conjugate gradients coefficient
         cg_coeff = (new_g_dot_g - polak) / old_g_dot_g

         ! Reset CG coefficient if too large
         if (abs(cg_coeff) > 5.0_DP) then
            if (pub_on_root .and. (pub_output_detail >= VERBOSE)) &
                 write(stdout,'(a/a,e10.3,a)') &
                 'WARNING in mermin_denskernel_optimise_cg: ', &
                 '  conjugate gradient coefficient (', &
                 cg_coeff, ') too large - resetting'
            cg_coeff = 0.0_DP
            cg_count = 0
         end if

      end if

      ! Increment CG counter or reset as appropriate
      if (abs(cg_coeff) > eps) then
         cg_count = cg_count + 1
      else
         cg_count = 0
      end if

      ! Reset conjugate gradients after cg_max steps
      if (cg_count > cg_max) then
         cg_coeff = 0.0_DP
         cg_count = 0
      end if

      ! Calculate the line search direction
      call sparse_embed_array_copy(con_direction, con_gradient)
      !ep: Calculate diffent directions for the different components
      if (pub_print_qc) then
         call sparse_embed_array_copy(entr_dir_comp, entr_grad_comp)
         call sparse_embed_array_scale(entr_dir_comp, -real(pub_num_spins,kind=DP))
         call sparse_embed_array_axpy(entr_dir_comp, old_entr_dir_comp, cg_coeff)
         call sparse_embed_array_copy(dft_dir_comp, dft_grad_comp)
         call sparse_embed_array_scale(dft_dir_comp, -real(pub_num_spins,kind=DP))
         call sparse_embed_array_axpy(dft_dir_comp, old_dft_dir_comp, cg_coeff)
         call sparse_embed_array_copy(mu_dir_comp, mu_grad_comp)
         call sparse_embed_array_scale(mu_dir_comp, -real(pub_num_spins,kind=DP))
         call sparse_embed_array_axpy(mu_dir_comp, old_mu_dir_comp, cg_coeff)
      end if
      ! pdh: hack to keep step lengths the same for spin polarised
      ! pdh: systems: multiply by a factor of two
      call sparse_embed_array_scale(con_direction, -real(pub_num_spins,kind=DP))
      call sparse_embed_array_axpy(con_direction, old_con_direction, cg_coeff)

      ! Ensure that search direction preserves constraint to first order
      ! pdh: (not necessary for current algorithms)
      if  (pub_output_detail >= VERBOSE)  then
         ! KPOINTS_DANGER
         call sparse_embed_array_trace(trace_array, con_direction, rep%overlap)
         dltne = sum(trace_array)
         dltne = dltne * pub_spin_fac
         if (pub_on_root) write(stdout,'(a,f20.12)') 'DltNe tr(Dir*S): ',dltne
      end if

      ! Calculate the function derivative along the search direction
      ! KPOINTS_DANGER
      call sparse_embed_array_trace(trace_array, co_gradient, con_direction)
      Qslope = sum(trace_array)
      !ep: Calculate diffent directions for the different components
      if (pub_print_qc) then
          call sparse_embed_array_trace(trace_array,entr_grad_co,con_direction)
          Eslope = sum(trace_array)
          call sparse_embed_array_trace(trace_array, dft_grad_co,con_direction)
          DFTslope = sum(trace_array)
          call sparse_embed_array_trace(trace_array, mu_grad_co,con_direction)
          Muslope = sum(trace_array)
      end if

      ! Ensure search direction points 'downhill'
      if (Qslope > 0.0_DP) then

         ! First try steepest descent direction i.e. reset conjugate gradients
         if ( (pub_output_detail >= VERBOSE) .and. pub_on_root) &
              write(stdout,'(a/a)') 'WARNING in mermin_denskernel_optimise_cg: &
              &positive line search gradient -', &
              '  resetting conjugate gradients'

         ! KPOINTS_DANGER
         call sparse_embed_array_copy(con_direction, con_gradient)
         call sparse_embed_array_scale(con_direction, -1.0_DP)
         !ep: Calculate diffent directions for the different components
         if (pub_print_qc) then
            call sparse_embed_array_copy(entr_dir_comp, entr_grad_comp)
            call sparse_embed_array_scale(entr_grad_comp, -1.0_DP)
            call sparse_embed_array_copy(dft_dir_comp, dft_grad_comp)
            call sparse_embed_array_scale(dft_grad_comp, -1.0_DP)
            call sparse_embed_array_copy(mu_dir_comp, mu_grad_comp)
            call sparse_embed_array_scale(mu_grad_comp, -1.0_DP)
         end if
         ! KPOINTS_DANGER
         call sparse_embed_array_trace(trace_array, co_gradient, con_direction)
         Qslope = sum(trace_array)
         !ep: Calculate diffent directions for the different components
         if (pub_print_qc) then
            call sparse_embed_array_trace(trace_array,entr_grad_co,con_direction)
            Eslope = sum(trace_array)
            call sparse_embed_array_trace(trace_array, dft_grad_co,con_direction)
            DFTslope = sum(trace_array)
            call sparse_embed_array_trace(trace_array, mu_grad_co,con_direction)
            Muslope = sum(trace_array)
         end if

         ! If this is still positive then things are really pear-shaped i.e.
         ! overlap matrix isn't positive definite or something...
         ! Reverse search direction!
        if (Qslope > 0.0_DP) then
           if ( (pub_output_detail >= VERBOSE ) .and. pub_on_root) &
                write(stdout,'(a)') 'WARNING in mermin_denskernel_optimise_cg: &
                &reversing search gradient'
           call sparse_embed_array_scale(con_direction, -1.0_DP)
           Qslope = -Qslope
        end if

     end if

      if ((pub_output_detail >= VERBOSE) .and. pub_on_root) &
           write(stdout,'(a,f20.12)') '             G0: ',Qslope

       ! Store quantities for next iteration
      call sparse_embed_array_copy(old_con_direction, con_direction)
      !ep: Calculate diffent directions for the different components
      if (pub_print_qc) then
         call sparse_embed_array_copy(old_entr_dir_comp, entr_dir_comp)
         call sparse_embed_array_copy(old_dft_dir_comp, dft_dir_comp)
         call sparse_embed_array_copy(old_mu_dir_comp, mu_dir_comp)
      end if
      if (pub_mermin_cg_type == 'MERMIN_POLAK') then
         call sparse_embed_array_copy(old_co_gradient, co_gradient)
      end if
      old_g_dot_g = new_g_dot_g

    end subroutine internal_search_direction

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_kernel_protected_step(trial_kernel,trial_step, &
         kernel_var, con_dirn)

      !===================================================================!
      ! This subroutine ensures that a chosen line step does not cause    !
      ! kernel occupation numbers to become unstable, by checking the     !
      ! occupancy bounds and RMS occupancy error and reducing the line    !
      ! step until they are within an acceptable stability range.         !
      !-------------------------------------------------------------------!
      ! Written by Nicholas Hine taken by the lnv module.                 !
      !===================================================================!

      use kernel, only: kernel_rms_err, kernel_occupancy_bounds

      implicit none

      ! Arguments
      type(DKERN), intent(inout) :: trial_kernel
      real(kind=DP), intent(inout)     :: trial_step
      type(DKERN), intent(inout) :: kernel_var
      type(SPAM3_EMBED_ARRAY), intent(in) :: con_dirn

      ! Local Variables
      !poli: here we use LNV parameter but can/should be changed
      real(kind=DP), parameter :: minbound = -0.36602540378443864676_DP
      real(kind=DP), parameter :: maxbound = 1.36602540378443864676_DP
      !real(kind=DP), parameter :: minbound = -1.0_DP
      !real(kind=DP), parameter :: maxbound = 2.0_DP
      !poli
      real(kind=DP) :: rms_err, trial_rms_err
      real(kind=DP), &
           dimension(kernel_var%kern%num_spins, kernel_var%kern%num_kpoints) :: &
           min_occ, max_occ
      real(kind=DP), &
           dimension(trial_kernel%kern%num_spins, trial_kernel%kern%num_kpoints) :: &
           trial_min_occ, trial_max_occ
      integer, parameter :: max_trials=10
      integer :: step

      ! ndmh: find rms occupancy error of starting kernel
      rms_err = kernel_rms_err(kernel_var, rep%overlap)

      ! ndmh: find occupancy bounds of starting kernel
      call kernel_occupancy_bounds(max_occ, min_occ, kernel_var, rep%overlap)

      ! ndmh: starting at input value of trial_step, keep reducing trial_step
      ! ndmh: until an acceptable value is found
      do step=1,max_trials

         ! Take trial step along the search direction
         call sparse_embed_array_copy(trial_kernel%kern, kernel_var%kern)
         call sparse_embed_array_axpy(trial_kernel%kern, con_dirn, trial_step)
         call kernel_workspace_invalidate(trial_kernel)

         ! ndmh: find rms occupancy error of trial kernel
         trial_rms_err = kernel_rms_err(trial_kernel, rep%overlap)

         ! ndmh: find occupancy bounds of trial kernel
         call kernel_occupancy_bounds(trial_max_occ, trial_min_occ, &
              trial_kernel, rep%overlap)

         ! ndmh: reduce step length by half if rms err increased by >20% or if
         ! ndmh: the new kernel has occupancies outside stable bounds
         ! KPOINTS_DANGER
        if ((trial_rms_err > rms_err+0.2_DP).or. &
              (any(trial_min_occ(:,:)<minbound) .or. &
              any(trial_max_occ(:,:)>maxbound))) then
            if (pub_on_root .and. (pub_output_detail>=VERBOSE)) then
               if (step==1) then
                  write(stdout,'(a)') '========================= Kernel &
                       &Line Step Protection =========================='
                  write(stdout,'(a)') '|  Iter |  Trial Step |   RMS &
                       &Occupancy Error |   Occupancy Bounds | Accepted? |'
                  write(stdout,'(i5,f16.6,f24.6,4x,"[",f7.4,":",f8.4,"]")') &
                       0,0.0_DP,rms_err,maxval(max_occ),minval(min_occ)
               end if
               write(stdout,'(i5,f16.6,f24.6,4x,"[",f7.4,":",f8.4,"]",a)') &
                    step,trial_step,trial_rms_err,maxval(trial_max_occ), &
                    minval(trial_min_occ),'        No'
            end if
            trial_step = trial_step * 0.5_DP

         else
            ! ndmh: this trial step was OK, so exit
            if ((step > 1) .and. (pub_output_detail>=VERBOSE) .and. &
                 pub_on_root) then
               write(stdout,'(i5,f16.6,f24.6,4x,"[",f7.4,":",f8.4,"]",a)') &
                    step,trial_step,trial_rms_err,maxval(trial_max_occ), &
                    minval(trial_min_occ),'       Yes'
               write(stdout,'(a)') '================================&
                    &================================================'
            end if
            return
         end if

      end do

      ! ndmh: finished loop without finding a stable kernel, so warn about this
      if (pub_on_root) then
         write(stdout,'(a,i2)') 'WARNING in mermin_denskernel_optimise_cg: Trial &
              &step reduced ',max_trials
         write(stdout,'(a)') 'times but kernel may still be unstable'
      end if

    end subroutine internal_kernel_protected_step

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

      !===================================================================!
      ! Taken by the lnv module.                                          !
      !===================================================================!

    subroutine internal_print

      implicit none

      ! Local variables
      integer :: is
      integer :: iter,numiter


      if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then

         if (iteration == 1) then
            write(stdout,'(i5,es27.12,2es15.4)') iteration, &
                 total_energy, rms_gradient, commutator
         else
            write(stdout,'(i5,es27.12,3es15.4)') iteration,total_energy, &
                 rms_gradient, commutator, total_energy - old_total_energy
         endif
      end if

      old_total_energy = total_energy

    end subroutine internal_print

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

   subroutine mermin_print_qc(mu, iter, rms_gradient, &
              denskern, gradient, step, slope, ent_slope, &
              dft_slope, mu_slope, free_energy, rep, &
              fd_slope, fd_ent_slope, fd_dft_slope, fd_mu_slope)

   !===================================================================!
   ! Prints the QC test results for the Mermin module.                 !
   !-------------------------------------------------------------------!
   ! Written by Emiliano Poli                                          !
   !===================================================================!


      use comms, only: pub_on_root
      use constants, only: k_B, DP, paw_en_size
      use hamiltonian, only: hamiltonian_dens_dep_matrices
      use rundat, only: pub_spin_fac, pub_num_spins, &
          pub_mermin_smearing_width, &
          pub_mermin_cheb
      use sparse_embed, only: sparse_embed_array_create, &
          sparse_embed_array_product, sparse_embed_array_destroy, &
          SPAM3_EMBED_ARRAY
      use smearing_operator, only: calculate_fermi_entropy_mermin

      use utils, only: utils_qc_print

      implicit none

      real(kind=DP), intent(in) :: mu(pub_num_spins, pub_num_kpoints)
      integer, intent(in)       :: iter
      real(kind=DP), intent(in) :: rms_gradient
      type(DKERN), intent(inout)   :: denskern
      type(SPAM3_EMBED_ARRAY), intent(in) :: gradient
      real(kind=DP), intent(in) :: step
      real(kind=DP), intent(in) :: slope
      real(kind=DP), intent(in) :: ent_slope
      real(kind=DP), intent(in) :: dft_slope
      real(kind=DP), intent(in) :: mu_slope
      real(kind=DP), intent(in) :: free_energy
      type(NGWF_REP), intent(in) :: rep
      real(kind=DP), intent(in) :: fd_slope
      real(kind=DP), intent(in) :: fd_ent_slope
      real(kind=DP), intent(in) :: fd_dft_slope
      real(kind=DP), intent(in) :: fd_mu_slope
      character(len=3) :: spin_string

      integer :: is
      real(kind=DP)  :: smermin(pub_num_spins)
      real(kind=DP), &
           dimension(denskern%kern%num_spins) :: trace
      real(kind=DP), &
           dimension(denskern%kern%num_spins) :: trace_kern
      real(kind=DP), &
           dimension(denskern%kern%num_spins) :: entropy_energy
      real(kind=DP), &
           dimension(denskern%kern%num_spins) :: mu_energy
      type(SPAM3_EMBED_ARRAY) :: temporary
      real(kind=DP) :: dft_energy
      real(kind=DP) :: paw_sphere_energies(paw_en_size)
      logical :: spoil_force=.false.
      logical :: update_ham=.false.

      call hamiltonian_dens_dep_matrices(ham, lhxc_fine, dft_energy, &
           lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
           ngwf_basis, hub_proj_basis, hub, denskern%kern, mdl,&
           hfxstate, update_ham, lhxc_fixed=.false.,spoil_force=spoil_force,&
           dfdtau_fine = dfdtau_fine, kpt=loc_kpt)

      ! Calculate ( G^ S_ )
      call sparse_embed_array_create(temporary, denskern%kern )
      call sparse_embed_array_product(temporary, con_gradient, rep%overlap)

      do is = 1, pub_num_spins
         call calculate_fermi_entropy_mermin(denskern%kern%m(is,PUB_1K), rep%overlap,1,smermin(is), pub_mermin_cheb)
         call sparse_embed_trace(trace(is), temporary%m(is,PUB_1K))
         call sparse_embed_trace(trace_kern(is),denskern%kern%m(is,PUB_1K),rep%overlap)
         entropy_energy(is) =smermin(is)*pub_spin_fac*(pub_mermin_smearing_width/k_B)
         mu_energy(is) =  mu(is,PUB_1K) * ((trace_kern(is)*pub_spin_fac) &
                         -(rep%n_occ(is,PUB_1K)*pub_spin_fac))
      end do

      if (pub_on_root) then
         call utils_qc_print('iter',iter)
         call utils_qc_print('Chebyshev expansion',pub_mermin_cheb)
         call utils_qc_print('rms_gradient',rms_gradient)
         call utils_qc_print('final_slope', slope)
         call utils_qc_print('final_fd_slope', fd_slope)
         call utils_qc_print('entropy_slope', ent_slope)
         call utils_qc_print('entropy_fd_slope', fd_ent_slope)
         call utils_qc_print('dft_slope', dft_slope)
         call utils_qc_print('dft_fd_slope', fd_dft_slope)
         call utils_qc_print('free_energy (a.u.)', free_energy )
         call utils_qc_print('dft_contrib. (a.u.)',dft_energy)
         do is = 1, pub_num_spins
            write(spin_string,'(a1,i1,a1)') '(', is, ')'
            call utils_qc_print('entr_contr. (a.u.)'//spin_string,entropy_energy(is))
            call utils_qc_print('mu_contrib. (a.u.)'//spin_string,mu_energy(is))
            call utils_qc_print('mu'//spin_string,mu(is,PUB_1K))
         end do
      end if

      call sparse_embed_array_destroy(temporary)

    end subroutine mermin_print_qc

    !_____________________________________________________________________________
    !_____________________________________________________________________________


    real(kind=DP) function mermin_trial_step(npts, xx, yy)

      !==========================================================================!
      ! Calculates the next trial step based on the history of previous steps,   !
      ! using the golden section polynomial search method.                       !
      !--------------------------------------------------------------------------!
      ! Adapted by EDFT module written by Alvaro Ruiz Serrano                    !
      !==========================================================================!

      use rundat, only: pub_mermin_cg_max_step
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
      real(kind=DP), parameter :: golden_section = 2.0_DP*(3.0_DP-sqrt(5.0_DP))*0.5_DP

      ! -------------------------------------------------------------------------

      if (pub_debug_on_root) write(stdout,'(a)') &
           "DEBUG: MERMIN: Entering mermin_trial_step"

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
      right_xx=pub_mermin_cg_max_step; right_found=.false.
      do nn = 1, npts
         if (xx(nn).gt.min_xx) then
            if (xx(nn).lt.right_xx) then
               right_xx = xx(nn)
               right_found = .true.
            end if
         else
            !right_xx = pub_mermin_cg_max_step
            right_xx = trial_step * 5.0_DP
            right_found = .false.
         end if
      end do

      ! ars: calculate step - move from the minimum into the largest subinterval
      if ( (left_found).and.(right_found) ) then
         if ( (abs(right_xx-min_xx)).gt.(abs(left_xx-min_xx)) ) then
            mermin_trial_step = min_xx + golden_section*abs(right_xx-min_xx)
         else
            mermin_trial_step = min_xx - golden_section*abs(left_xx-min_xx)
         end if
      else
         mermin_trial_step = min_xx + golden_section*abs(right_xx-min_xx)
      end if

      if(pub_debug_on_root) then
         write(stdout,'(a,f22.14)') "DEBUG: min_xx = ", min_xx
         write(stdout,'(a,f22.14)') "DEBUG: left_xx = ", left_xx
         write(stdout,'(a,l1)') "DEBUG: left_found = ", left_found
         write(stdout,'(a,f22.14)') "DEBUG: right_xx = ", right_xx
         write(stdout,'(a,l1)') "DEBUG: right_found = ", right_found
         write(stdout,'(a,f22.14)')&
              "DEBUG: mermin_trial_step = ", mermin_trial_step
         write(stdout,'(a)') "DEBUG: MERMIN: Leaving mermin_trial_step"
      end if

    end function mermin_trial_step

  !-----------------------------------------------------------------------------

    subroutine mermin_evaluate_step(aux, trial_aux, con_direction, trial_step)

      !==========================================================================!
      ! Prints the result of the current line search step on the screen.         !
      !--------------------------------------------------------------------------!
      ! Adapted by EDFT module written by Alvaro Ruiz Serrano                    !
      !==========================================================================!

      use rundat, only: pub_mermin_check_trial_steps

      implicit none

      ! ars: arguments
      type(DKERN), intent(inout)    :: aux
      type(DKERN), intent(inout)    :: trial_aux
      type(SPAM3_EMBED_ARRAY), intent(in) :: con_direction     ! Contra search direction
      real(kind=DP), intent(inout) :: trial_step


      if (pub_mermin_check_trial_steps) then
         call internal_kernel_protected_step(trial_aux, trial_step, &
              aux,con_direction)
      else

         ! Take a trial step along the search direction
         call sparse_embed_array_copy(trial_aux%kern, aux%kern)
         call sparse_embed_array_axpy(trial_aux%kern, &
              con_direction, trial_step)
         call kernel_workspace_invalidate(trial_aux)

      end if

    end subroutine mermin_evaluate_step


  !-----------------------------------------------------------------------------

  end subroutine mermin_denskernel_optimise_cg

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////


  subroutine mermin_calculate_mu(mu, aux, ham, rep, ngwf_basis)

    !========================================================================!
    ! This subroutine calculates the Chemical Potential using the Millam- &  !
    ! Scuseria method to obtain a traceless gradient.                        !
    ! J. Chem. Phys. 106, 5569 (1997)                                        !
    ! In this module the Chemical Potential is defined as:                   !
    ! mu = -(tr(H^ab S_ba) +  T * tr({(S_**D(K^**S_**)}^ab S_ba ))) / N      !
    ! where:                                                                 !
    ! mu = Chemical Potential (Defined by Millam Scuseria method)            !
    ! H^ab = Hamiltonian                                                     !
    ! T = Temperature                                                        !
    ! D(K^**S_**) = Chebyshev expansion of the derivative of the Entropy     !
    ! S_ba, S_** = overlap matrix                                            !
    ! N = num NGWFs                                                          !
    !------------------------------------------------------------------------!
    ! Written by Poli Emiliaino                                              !
    ! Based on an original version by Chris-Kriton Skylaris                  !
    !========================================================================!

    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use kernel, only: DKERN, &
         kernel_workspace_create, &
         kernel_workspace_destroy
    use ngwf_representation, only: NGWF_HAM, NGWF_REP
    use rundat, only: pub_debug_on_root, &
         pub_num_kpoints, PUB_1K
    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: mu(:,:)    ! Chemical potential
    type(DKERN), intent(inout) :: aux
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(NGWF_HAM), intent(inout) :: ham
    type(NGWF_REP), intent(in) :: rep

    ! Local variables
    logical :: deallocate_workspace


    if (pub_debug_on_root) then
       write(stdout,'(a)') 'DEBUG: Entering mermin_calculate_mu'
    end if

    ! jme: KPOINTS_DANGER
    ! The algorithms in this module have not been checked for more than 1 k-point
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine mermin_calculate_mu not checked yet for more&
         & than one k-point.')

    ! Check for workspace
    deallocate_workspace = .false.
    if (.not.aux%workspace_created) then
       call kernel_workspace_create(aux, rep%overlap)
       deallocate_workspace = .true.
    end if

    call internal_mu()

    if (deallocate_workspace) call kernel_workspace_destroy(aux)

    if (pub_debug_on_root) then
       write(stdout,'(a)') 'DEBUG: Leaving mermin_calculate_mu'
    end if

contains

    !ep
    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_mu()

      use constants, only: k_B
      use smearing_operator, only: calculate_fermi_entropy_mermin
      use rundat, only: pub_num_spins, pub_mermin_smearing_width, &
          PUB_1K, pub_mermin_cheb
      use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_copy, &
         sparse_embed_destroy, sparse_embed_trace, &
         sparse_embed_array_create, sparse_embed_array_destroy, &
         SPAM3_EMBED_ARRAY, sparse_embed_product
      use kernel, only: kernel_validate_ks
      use sparse_embed, only: sparse_embed_trace
      implicit none

      ! Local variables
      integer :: is, ik, ierr
      real(kind=DP) :: smermin
      type(SPAM3_EMBED_ARRAY) :: deriv, derivok, temp, tempb
      type(SPAM3_EMBED), allocatable :: hamtmp(:), hamtmp2(:)
      real(kind=DP), &
           dimension(aux%kern%num_spins) :: trace
      real(kind=DP), &
           dimension(aux%kern%num_spins) :: trace2

      call kernel_validate_ks(aux, rep%overlap)

      allocate(hamtmp(1:pub_num_spins), stat=ierr)
      allocate(hamtmp2(1:pub_num_spins), stat=ierr)

      do is = 1, pub_num_spins
          call sparse_embed_create(hamtmp(is), aux%kern%m(is,PUB_1K))
          call sparse_embed_create(hamtmp2(is), aux%kern%m(is,PUB_1K))
          call sparse_embed_copy(hamtmp(is), ham%ham(is))
          call sparse_embed_copy(hamtmp2(is), ham%ham(is))
      end do

      call sparse_embed_array_create(deriv, aux%kern)
      call sparse_embed_array_create(derivok, aux%kern)
      call sparse_embed_array_create(temp, aux%kern)
      call sparse_embed_array_create(tempb, aux%kern)

      !ep:  CALCULATE tr(H^ S_) = tr(S^ H_ S^ S_)
      do is = 1, pub_num_spins
       !ep: hamtmp = S^ H_
         call sparse_embed_product(hamtmp(is), rep%inv_overlap, ham%ham(is))
       !ep: hamtmp2 = hamtmp S^ = S^ H_ S^
         call sparse_embed_product(hamtmp2(is), hamtmp(is), rep%inv_overlap)
       !ep: tr(hamtmp2 S_) = S^ H_ S^ S_
         call sparse_embed_trace(trace(is),hamtmp2(is), rep%overlap)
      end do

      do ik = 1, aux%kern%num_kpoints
         do is = 1, pub_num_spins
            call calculate_fermi_entropy_mermin(aux%kern%m(is,ik), &
            & rep%overlap,2,smermin, pub_mermin_cheb,deriv%m(is,ik))
            !ep: CALCULATE S dS
            call sparse_embed_product(temp%m(is,ik), rep%overlap, deriv%m(is,ik))
            !ep: CALCULATE S^-1 * (S dS)
            call sparse_embed_product(tempb%m(is,ik), rep%inv_overlap, temp%m(is,ik))
            !ep: CALCULATE S^-1 * (S dS) * S^-1
            call sparse_embed_product(temp%m(is,ik), tempb%m(is,ik), rep%inv_overlap)
            !ep: CALCULATE [S^-1 * (S dS) * S^-1] * S
            call sparse_embed_product(derivok%m(is,ik), temp%m(is,ik), rep%overlap)
            !ep: CALCULATE T*tr{ [S^-1 * (S dS) * S^-1 ] * S) and sum to tr(H^ S_) = tr(S^ H_ S^ S_)
            call sparse_embed_trace(trace2(is), derivok%m(is,ik))
            mu(is,ik) = (-1.0 * ((trace(is)) + ((pub_mermin_smearing_width/k_B)*trace2(is)))) / &
                        (sum(ngwf_basis(:)%num))
         end do
      end do

      ! Deallocate workspace
      do is = 1, pub_num_spins
         call sparse_embed_destroy(hamtmp(is))
         call sparse_embed_destroy(hamtmp2(is))
      enddo
      call sparse_embed_array_destroy(deriv)
      call sparse_embed_array_destroy(derivok)
      call sparse_embed_array_destroy(temp)
      call sparse_embed_array_destroy(tempb)

    end subroutine internal_mu
    !ep

  end subroutine mermin_calculate_mu

  !====ep
  subroutine mermin_ngwf_gradient(pmat, qmat, denskern, inv_overlap, overlap)
  !====ep

    !========================================================================!
    ! This subroutine calculates the gradient of the Free Energy             !
    ! wrt to the NGWFs                                                       !
    ! In this module the Free Energy gradient is defined as:                 !
    ! dA/dphi*_c =  H|phi_a>K^ac - T * dS/dphi*_c + mu*|phi_a>K^ac           !
    ! where:                                                                 !
    ! dA/dphi*_c = Derivative of the Free energy wrt NGWFs                   !
    ! phi_a, phi*_c = NGWFs                                                  !
    ! K^ac = Density Kernel                                                  !
    ! H = Hamiltonian                                                        !
    ! T = Temperature                                                        !
    ! dS/dphi*_c = Derivative of the Entropy wrt to the NGWFs, for the  &    !
    ! equation see the smearing module calculate_fermi_entropy_mermin sub    !
    ! mu = Chemical Potential (Defined by Millam Scuseria method)            !
    ! Note: the mu contribution to the gradient is caclculated in the NGWF & !
    ! gradient module                                                        !
    !------------------------------------------------------------------------!
    ! Written by Poli Emiliano                                               !
    !========================================================================!

    use constants, only: DP, stdout, k_B
    use kernel, only: DKERN
    use rundat, only: pub_num_spins, pub_mermin_smearing_width, &
        PUB_1K, pub_debug_on_root, pub_mermin_cheb
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy, &
         SPAM3_EMBED_ARRAY, sparse_embed_array_create, &
         sparse_embed_array_scale, &
         sparse_embed_array_destroy
    use smearing_operator, only: calculate_fermi_entropy_mermin
    use utils, only:  utils_unit

    implicit none

    type(SPAM3_EMBED), intent(inout) :: pmat(1:pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: qmat(1:pub_num_spins)
    type(DKERN), intent(in   ) :: denskern
    type(SPAM3_EMBED), intent(in   ) :: inv_overlap
    type(SPAM3_EMBED), intent(in   ) :: overlap

    integer :: is
    real(kind=DP) :: smermin
    type(SPAM3_EMBED_ARRAY) :: deriv

    if (pub_debug_on_root) write(stdout,'(a)') &
         "DEBUG: Entering mermin_ngwf_gradient"


    call sparse_embed_array_create(deriv, denskern%kern)

    ! ep: calculate dA/dphi*_c
    do is = 1, pub_num_spins
       call calculate_fermi_entropy_mermin(denskern%kern%m(is,PUB_1K), overlap, 3 ,smermin, pub_mermin_cheb, &
            deriv%m(is,PUB_1K), inv_overlap)
    end do

    ! ep: calculate T * dA/dphi*_c
    call sparse_embed_array_scale(deriv, (pub_mermin_smearing_width/k_B), 0.0_DP)

    do is = 1, pub_num_spins

       call sparse_embed_copy(pmat(is), denskern%kern%m(is,PUB_1K))
       ! ep: copy T * dA/dphi*_c to QMATRIX
       call sparse_embed_copy(qmat(is), deriv%m(is,PUB_1K))

    end do

    call sparse_embed_array_destroy(deriv)

    if (pub_debug_on_root) write(stdout,'(a)') &
         "DEBUG: Leaving mermin_ngwf_gradient"

  end subroutine mermin_ngwf_gradient


  real(kind=DP) function mermin_lagrangian(total_energy, denskern, &
                 & overlap, mu, n_occ)

    !========================================================================!
    ! This subroutine calculates the mermin lagrangian                       !
    !------------------------------------------------------------------------!
    ! Written by Emiliano Poli                                               !
    !========================================================================!

    use constants, only: DP, stdout, k_B, max_spins
    use rundat, only: pub_num_spins, pub_spin_fac, pub_num_kpoints, &
        PUB_1K, pub_debug_on_root, pub_mermin_smearing_width, &
        pub_mermin_cheb
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_trace
    use smearing_operator, only: calculate_fermi_entropy_mermin
    use utils, only: utils_abort

    implicit none

    ! ars: Arguments

    real(kind=DP), intent(in) :: total_energy
    type(SPAM3_EMBED),   intent(in) :: denskern(:)
    type(SPAM3_EMBED),   intent(in) :: overlap
    real(kind=DP), intent(in) :: mu(max_spins)
    integer,       intent(in) :: n_occ(max_spins)

    ! Local variable
    integer :: is
    real(kind=DP) :: smermin(2)
    real(kind=DP) :: tr_tmp(2)
    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering mermin_lagrangian'

    ! KPOINTS_DANGER
    ! First k-point enforced in all uses of n_occ via PUB_1K
    if (pub_num_kpoints /= PUB_1K) then
       call utils_abort('Function mermin_lagrangian not yet compatible with more&
            & than one k-point.')
    end if



    mermin_lagrangian = total_energy
    !ep: CHECK COMPONENTS DEBUG
    !mermin_lagrangian = 0.0_DP
    !ep: CHECK COMPONENTS DEBUG
    do is = 1, pub_num_spins
      call calculate_fermi_entropy_mermin(denskern(is), overlap, 1, smermin(is), pub_mermin_cheb)
      call sparse_embed_trace(tr_tmp(is),denskern(is),overlap)
      if (pub_num_spins .eq. 2) then
          mermin_lagrangian = mermin_lagrangian  + (smermin(is) *(pub_mermin_smearing_width/k_B)) + &
                           mu(is) * ( tr_tmp(is) - real(n_occ(is),kind=DP))
          !ep CHECK COMPONENTS ENTRO
          !mermin_lagrangian  = mermin_lagrangian  +(smermin(is) * (pub_mermin_smearing_width/k_B))
          !ep CHECK COMPONENTS MU
          !mermin_lagrangian  = mermin_lagrangian  + mu(is) *  (sparse_embed_trace(denskern(is),overlap) - &
          !                     real(n_occ(is),kind=DP))
      else
          mermin_lagrangian = mermin_lagrangian  + (smermin(is) *(pub_spin_fac*(pub_mermin_smearing_width/k_B))) + &
                           pub_spin_fac*mu(is) *  (tr_tmp(is) - real(n_occ(is),kind=DP))
          !ep: CHECK COMPONENTS ENTRO
          !mermin_lagrangian = mermin_lagrangian  + (smermin(is) * (pub_spin_fac*(pub_mermin_smearing_width/k_B)))
          !ep: CHECK COMPONENTS MU
          !mermin_lagrangian = mermin_lagrangian  + pub_spin_fac*mu(is) *  (sparse_embed_trace(denskern(is),overlap) - &
          !                    real(n_occ(is),kind=DP))
      end if
    end do


  end function mermin_lagrangian

  subroutine mermin_switch_mu(total_energy, overlap, denskern, &
                 &  mu, n_occ)

    !========================================================================!
    ! This subroutine calculates the mermin lagrangian                       !
    !------------------------------------------------------------------------!
    ! Written by Emiliano Poli                                               !
    !========================================================================!

    use constants, only: DP, stdout, k_B, max_spins
    use kernel, only: DKERN
    use rundat, only: pub_num_spins, pub_spin_fac, pub_num_kpoints, &
        PUB_1K, pub_debug_on_root, pub_mermin_smearing_width, &
        pub_mermin_cheb
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_trace
    use smearing_operator, only: calculate_fermi_entropy_mermin
    use utils, only: utils_abort

    implicit none

    ! ars: Arguments

    real(kind=DP), intent(inout) :: total_energy
    type(SPAM3_EMBED),   intent(in) :: denskern(:)
    type(SPAM3_EMBED),   intent(in) :: overlap
    real(kind=DP), intent(in) :: mu(max_spins)
    integer,       intent(in) :: n_occ(max_spins)

    ! Local variable
    integer :: is
    real(kind=DP) :: smermin(2)
    real(kind=DP) :: tr_tmp(2)
    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering mermin_lagrangian'

    ! KPOINTS_DANGER
    ! First k-point enforced in all uses of n_occ via PUB_1K
    if (pub_num_kpoints /= PUB_1K) then
       call utils_abort('Function mermin_lagrangian not yet compatible with more&
            & than one k-point.')
    end if

    do is = 1, pub_num_spins
      call calculate_fermi_entropy_mermin(denskern(is), &
           overlap, 1, smermin(is), pub_mermin_cheb)
      call sparse_embed_trace(tr_tmp(is),denskern(is),overlap)
      if (pub_num_spins .eq. 2) then
         total_energy = total_energy - (smermin(is) * (pub_mermin_smearing_width/k_B)) - &
                        mu(is) * (tr_tmp(is) - real(n_occ(is),kind=DP))
      else
         total_energy = total_energy - (smermin(is) * (pub_spin_fac*(pub_mermin_smearing_width/k_B))) - &
                        pub_spin_fac*mu(is) * (tr_tmp(is) - real(n_occ(is),kind=DP))
        end if
    enddo

  end subroutine mermin_switch_mu

end module mermin
