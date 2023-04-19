! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!               Density kernel optimisation module               !
!                                                                !
! This module implements several variants of the LNV method for  !
! optimising the density kernel K for a fixed set of NGWFs. The  !
! density kernel K is expressed in terms of an auxiliary matrix  !
! L by a generalised purifying transformation:                   !
!   K = f(L) [ 3 L.S.L - 2 L.S.L.S.L ]                           !
! where f(L) is a scalar function of L. The variants differ only !
! in the method used to impose the normalisation constraint.     !
!----------------------------------------------------------------!
! Written by Peter Haynes, 17/11/04                              !
! Based on an original version by Chris-Kriton Skylaris in 2000  !
! with subsequent modifications by Chris-Kriton Skylaris,        !
! Arash Mostofi and Peter Haynes.                                !
! Modified by Nicholas Hine to use the NGWF_REP container for    !
! NGWF-related data and matrices, October 2010.                  !
! Modified by Robert Charlton to use SPAM3_EMBED structures and  !
! multiple NGWF bases, July 2017.                                !
!================================================================!

module lnv

  implicit none

  private

  public :: lnv_denskernel_optimise_cg
  public :: lnv_calculate_mu
  public :: lnv_reset_trial_step_to_initial

  logical, save :: lnv_reset_trial_step_to_initial = .false.

contains


  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////


  subroutine lnv_denskernel_optimise_cg(aux, pur_denskern, ham, lhxc_fine, mu, &
       total_energy, rep, ngwf_basis, hub_proj_basis, hub, &
       mdl, hfxstate, lnv_threshold, current_maxit_lnv, conv_status, &
       dfdtau_fine, kpt, force_no_IO)

    !======================================================================!
    ! This subroutine implements several variants of the LNV method for    !
    ! optimising the density kernel K for a fixed set of NGWFs. The        !
    ! density kernel K is expressed in terms of an auxiliary matrix L by a !
    ! generalised purifying transformation:                                !
    !   K = f(L) [ 3 L.S.L - 2 L.S.L.S.L ]                                 !
    ! where f(L) is a scalar function of L. The variants differ only in    !
    ! the method used to impose the normalisation constraint.              !
    !                                                                      !
    ! Methods implemented so far with references:                          !
    ! * Original Li-Nunes-Vanderbilt scheme: Phys. Rev. B 47, 10891 (1993) !
    ! * Millam-Scuseria variant: J. Chem. Phys. 106, 5569 (1997)           !
    !----------------------------------------------------------------------!
    ! Arguments:                                                           !
    ! aux           (inout): Auxiliary density kernel in DKERN format      !
    ! pur_denskern  (inout): Purified density kernel in DKERN format       !
    ! ham           (inout): Hamiltonian wrapper (contains SPAM3 matrices) !
    ! lhxc_fine       (out): Local-Hartree-Exhange-Correlation potential   !
    ! mu              (out): Chemical potential type parameter             !
    ! total_energy    (out): Total energy                                  !
    ! rep              (in): NGWF Representation (functions and matrices)  !
    ! ngwf_basis       (in): Function basis type describing the NGWFs.     !
    ! hub_proj_basis   (in): Function basis type describing Hubbard projs  !
    ! grid             (in): The GRID_INFO for all the whole-cell arrays   !
    ! lnv_threshold    (in): Convergence threshold for LNV gradient        !
    ! current_maxit_lnv(in): Maximum LNV iterations set in electronic mod  !
    ! dfdtau_fine   (inout): Gradient of XC energy per unit volume wrt     !
    !                        KE density energy (optional)                  !
    !----------------------------------------------------------------------!
    ! This version written by Peter Haynes, November 2004.                 !
    ! Based on an original version by Chris-Kriton Skylaris in 2000 with   !
    !   subsequent modifications by Chris-Kriton Skylaris, Arash Mostofi   !
    !   and Peter Haynes.                                                  !
    ! Modified by Peter Haynes for parallel SPAM 2, July 2006.             !
    ! Modified to include Nonlinear Core Corrections by Nicholas Hine in   !
    ! January 2009                                                         !
    ! DFT+U added by David O'Regan, April 2009                             !
    ! Adapted to use SPAM3 matrices by Nicholas Hine, June 2009            !
    ! Modified by Laura Ratcliff for conduction calculations, Oct 2010     !
    ! Tiny fix by Jacek Dziedzic to reset trial step for every new calc    !
    ! and not once per ONETEP run, Oct 2013.                               !
    ! Modified to use embedding matrices by Robert Charlton, July 2017.    !
    !======================================================================!

    use comms, only: pub_on_root, pub_root_proc_id, comms_bcast
    use constants, only: DP, stdout, VERBOSE, NORMAL
#ifdef GPU_PGI
#ifdef GPU_SP_TEST
    use fourier_gpu_wrapper_mod, only: GPU_SP, GPU_SP_switched, GPU_SP_switch, &
                                       deltaE_gpu
#endif
#endif
    use fragment_data, only: eda_dfdtau_fine
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_build_matrix, hamiltonian_energy_components
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: DKERN, kernel_workspace_invalidate, &
         kernel_workspace_create, kernel_workspace_destroy, kernel_normalise, &
         kernel_fix, kernel_create, kernel_destroy, kernel_rescale
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use restart, only : restart_kernel_write
    use rundat, only: pub_output_detail, pub_lnv_cg_type, pub_write_denskern, &
         pub_exact_lnv, pub_lnv_check_trial_steps,pub_lnv_cg_max_step, &
         pub_hub_calculating_u,pub_write_converged_dk_ngwfs,pub_devel_code, &
         pub_cond_calculate, pub_debug_on_root, &
         pub_num_spins, pub_num_kpoints, PUB_1K, &
         pub_rootname,mix_dkn_init_num,md_iter_global, &
         mix_dkn_type, pub_maxit_ngwf_cg, pub_dmft_spoil_kernel,pub_dmft_fully_sc, &
         pub_dmft_fully_sc_h, pub_inner_loop_iteration, &
         pub_xc_ke_density_required, &
         pub_emft, pub_emft_follow, pub_emft_lnv_only
    use services, only: services_flush
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_copy, &
         sparse_embed_destroy, sparse_embed_axpy, sparse_embed_trace, &
         sparse_embed_array_create, sparse_embed_array_copy, &
         sparse_embed_array_axpy, sparse_embed_array_trace, &
         sparse_embed_array_scale, sparse_embed_array_product, &
         sparse_embed_array_destroy, SPAM3_EMBED_ARRAY
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_banner, &
         utils_unit, utils_open_unit_check, utils_close_unit_check, utils_assert

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
    integer, intent(in) :: current_maxit_lnv       ! Max number of iterations
    real(kind=DP), intent(in) :: lnv_threshold     ! LNV threshold
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
    real(kind=DP) :: ne_pure          ! Electron number from purified kernel
    real(kind=DP) :: old_total_energy ! Total energy at previous LNV step
    real(kind=DP) :: commutator       ! RMS value of [L,H] commutator
    real(kind=DP) :: rms_gradient     ! RMS value of gradient of E wrt L
    real(kind=DP) :: line_search_step ! Line search optimal step length
    real(kind=DP) :: cg_coeff         ! Conjugate gradients coefficient
    real(kind=DP) :: lhxc_energy      ! Local PSP, Hartree and XC energy
    real(kind=DP) :: mu_term          ! Term associated with chemical potential
    real(kind=DP) :: rms_dkern        ! RMS value of change in density kernel
    real(kind=DP) :: old_g_dot_g      ! Previous scalar product of gradients
    real(kind=DP) :: new_g_dot_g      ! Current scalar product of gradients
    integer :: iteration              ! Iteration number
    integer :: cg_count               ! Number of CG steps since last reset
    integer, parameter :: cg_max=5    ! Maximum number of CG steps before reset
    logical :: converged              ! Flag to indicate convergence

    ! Line search variables
    real(kind=DP) :: Qinitial         ! Initial value of function
    real(kind=DP) :: Qslope           ! Initial slope along search direction
    real(kind=DP), save :: trial_step ! Trial step length
    real(kind=DP) :: Qtrial           ! Function value at trial step
    real(kind=DP) :: optimal_step     ! Optimal step length
    real(kind=DP) :: Qoptimal         ! Function value at optimal step
    real(kind=DP) :: Qpredict         ! Predicted value at optimal step
    real(kind=DP) :: Qminimum         ! Predicted value of function at minimum
    real(kind=DP) :: total_energy_at_trial_length
    real(kind=DP) :: hubbard_energy   ! Hubbard correction energy

    ! General variables
    real(kind=DP), parameter :: eps=epsilon(1.0_DP)
    real(kind=DP) :: norm_fac(pub_num_spins, pub_num_kpoints) ! Normalisation factor
    real(kind=DP), allocatable :: trace_array(:,:)
    integer :: is                     ! Spin counter
    integer :: ierr
    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecokpt
    real(kind=DP) :: loc_kpt(3)

    ! ndmh: for finite differencing
    integer, parameter :: nfd=10
    real(kind=DP) :: fd_trial_step(nfd)=(/0.000001_DP, 0.00001_DP,0.0001_DP, &
         0.001_DP,0.01_DP,0.1_DP,0.3_DP,0.4_DP,0.5_DP,1.0_DP/)
    ! jd: More detailed LNV_FD, with focus on regions where it matters:
!    integer, parameter :: nfd=19
!    real(kind=DP), parameter :: logd=1.58489_DP
!    real(kind=DP), parameter :: start = 0.00001_DP
!    real(kind=DP) :: fd_trial_step(nfd)=(/&
!         start*logd**0, start*logd**1, start*logd**2, start*logd**3, &
!         start*logd**4, start*logd**5, start*logd**6, start*logd**7, &
!         start*logd**8, start*logd**9, start*logd**10, start*logd**11, &
!         start*logd**12, start*logd**13, start*logd**14, start*logd**15, &
!         start*logd**17, start*logd**19, start*logd**21/)
    ! JCW: More closely-spaced steps between 0.0001 and 0.1 for use in debugging
    !integer, parameter :: nfd=17
    !real(kind=DP) :: fd_trial_step(nfd)=(/0.000001_DP, 0.00001_DP,&
    !     0.0001_DP, 0.00019952623149688788_DP, 0.00039810717055349735_DP, &
    !     0.0007943282347242813_DP, 0.001584893192461114_DP, &
    !     0.0031622776601683794_DP, 0.00630957344480193_DP, &
    !     0.012589254117941675_DP, 0.025118864315095794_DP, &
    !     0.0501187233627272_DP,&
    !     0.1_DP,0.3_DP,0.4_DP,0.5_DP,1.0_DP/)
    integer :: ifd
    real(kind=DP) :: QFD(nfd)
    real(kind=DP) :: fd_slope(nfd)

    ! vv : to output the number of LNV steps
    character(len=80) :: filename
    integer           :: io_unit, io_stat
    logical           :: loc_force_no_IO


    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
         &Entering lnv_denskernel_optimise_cg'

    ! Flush output
    call services_flush

    ! Start the timer
    call timer_clock('lnv_denskernel_optimise_cg', 1)

    if (.not.pub_xc_ke_density_required) then
       ! JCW: Non-tau-dependent XC functional
       ! JCW: dfdtau_fine is either not present, or has a length of 1
       if ( present(dfdtau_fine ) ) then
          call utils_assert( &
               size(dfdtau_fine) == 1,&
               "Error in lnv_denskernel_optimise_cg: &
               &pub_xc_ke_density_required is false, but dfdtau_fine is present &
               &and has a size > 1.")
       end if
    else
       ! JCW: tau-dependent meta-GGA XC functional
       ! JCW: Check dfdtau_fine has been passed in.
       call utils_assert(present(dfdtau_fine),"Error in &
            &lnv_denskernel_optimise_cg: pub_xc_ke_density_required is true &
            &but no dfdtau_fine array has been passed in.")
    end if

    ! jme: KPOINTS_DANGER
    ! The algorithms in this module have not been checked for more than 1 k-point
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine lnv_denskernel_optimise_cg not checked yet for more&
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
    if(lnv_reset_trial_step_to_initial) then
       trial_step=0.5_DP
       lnv_reset_trial_step_to_initial = .false.
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
    if(pub_lnv_cg_type/='LNV_FLETCHER') then
       call sparse_embed_array_create(old_co_gradient, aux%kern)
    endif

    allocate(trace_array(co_gradient%num_spins,co_gradient%num_kpoints),&
         stat=ierr)
    call utils_alloc_check('lnv_denskernel_optimise_cg','trace_array',ierr)

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    ! Ensure the kernel will be stable before optimisation starts
    ! ebl : the dmft dens kernel spoil flag may be needed here
    if (index(pub_devel_code,'RESCALE_AUX')>0) then
       call kernel_rescale(pur_denskern, rep%overlap, rep%n_occ)
    else
       call kernel_normalise(aux, rep%overlap, rep%inv_overlap, rep%n_occ)
       call kernel_fix(aux, rep%overlap, rep%inv_overlap, miniter=3)
    end if

    ! Initialise variables for conjugate gradients
    old_g_dot_g = 1.0e100_DP           ! Set previous product to dummy value
    line_search_step = 0.0_DP          ! Set line search step length to zero
    cg_count = 0                       ! Reset conjugate gradients counter
    converged = .false.                ! Initially no convergence

    if (pub_on_root .and. (pub_output_detail >= NORMAL)) then

       if (pub_output_detail >= VERBOSE) then
          write(stdout,'(a)') repeat('.',80)
          if (pub_exact_lnv) then
             write(stdout,'(a)') utils_banner('<', &
                  'LNV (Original version) density kernel optimisation')
          else
             write(stdout,'(a)') utils_banner('<', &
                  'LNV (Millam-Scuseria variant) density kernel optimisation')
          end if
          write(stdout,'(a)') repeat('~',80)
       end if

       write(stdout,'(a)')'   iter  |      energy (Eh)       | rms gradient |&
            &  commutator  |   dE (Eh)'
    end if

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! ** Start of main loop
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    do iteration=1,current_maxit_lnv
       ! jd: Globally accessible copy of iteration number, only used for contro-
       !     lling devel_code dumps of grid quantities from various modules
       pub_inner_loop_iteration = iteration

#ifdef GPU_SP_TEST
    ! kaw: Accuracy switching toggles
    if ((.not.GPU_SP_switched).and.(abs(deltaE_gpu).gt.(GPU_SP_switch))) then
       GPU_SP=.true.                                                    ! Turn on SP code
    else
       GPU_SP=.false.                                                   ! kaw: Switch to DP
       GPU_SP_switched=.true.
    end if
    if (iteration.eq.1) GPU_SP=.true.
    if (iteration.eq.1) GPU_SP_switched=.false.
#endif

       ! Flush output
       call services_flush

       ! cks: ======================= ENERGY AT POINT 0 =======================
       ! CW
       if(pub_dmft_spoil_kernel.and.iteration==1)then
          total_energy = internal_energy(aux,.true.,spoil_force=.true.)
       else
          total_energy = internal_energy(aux,.true.)
       endif
       ! END CW

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

          call lnv_calculate_mu(mu, aux, ham, rep, ngwf_basis)

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
       if (iteration /= current_maxit_lnv) then

          ! cks: ===========================LINE SEARCH========================
          ! Initial function value including constraint term
          mu_term = internal_mu_term(aux%kern)
          Qinitial = total_energy + mu_term

          if ( (pub_output_detail >= VERBOSE) .and. pub_on_root) &
             write(stdout,'(a,f24.12)') '             Q0: ', Qinitial

          call internal_search_direction ! <-- sets Qslope (JCW)

          if ( (pub_output_detail >= VERBOSE) .and. pub_on_root) &
             write(stdout,'(a,f20.12)') '   RMS gradient: ', rms_gradient

          ! Sophisticated line search is not suitable initially
          if (iteration > 4 .or. rms_gradient <= 0.05_DP) then

             if (index(pub_devel_code,'LNV_FD')>0) then

                if ( (pub_output_detail < VERBOSE) .and. pub_on_root) &
                     write(stdout,'(a,f20.12)') '             G0: ',Qslope

                ! ndmh: evaluate and print energy at several FD trial steps
                do ifd=1,nfd

                   if (pub_lnv_check_trial_steps) then
                      call internal_kernel_protected_step(trial_aux, &
                           fd_trial_step(ifd),aux,con_direction)
                   else
!
                      ! Take a trial step along the search direction
                      call sparse_embed_array_copy(trial_aux%kern, aux%kern)
                      call sparse_embed_array_axpy(trial_aux%kern, &
                           con_direction, fd_trial_step(ifd))
                      call kernel_workspace_invalidate(trial_aux)

                   end if

                   ! ndmh: %%%%%%%%%%%%%%%%%% ENERGY AT FD LENGTH %%%%%%%%%%%%%%
                   total_energy_at_trial_length = internal_energy(trial_aux, &
                        .false.)

                   mu_term = internal_mu_term(trial_aux%kern)

                   QFD(ifd) = total_energy_at_trial_length + mu_term
                   ! rc2013: check that the slopes are sane
                   fd_slope(ifd) = (QFD(ifd) - Qinitial)/fd_trial_step(ifd)
                   if ((pub_output_detail>=VERBOSE).and. &
                        (.not.pub_cond_calculate)) then
                      call hamiltonian_energy_components( &
                           pur_denskern%kern%m(:,PUB_1K), rep, mdl, &
                           ngwf_basis, hub_proj_basis, hub, &
                           ham%hfexchange)
                   end if

                   if (pub_on_root) then
                      write(stdout,'(a,3f20.12)') ' FD step, QFD, FD slope  : ', &
                          fd_trial_step(ifd),QFD(ifd),fd_slope(ifd)
                   end if
                   ! ndmh: %%%%%%%%%%%%%%%%%% END ENERGY AT FD LENGTH %%%%%%%%%%

                end do
             end if

             if (pub_lnv_check_trial_steps) then
                call internal_kernel_protected_step(trial_aux, trial_step, &
                     aux,con_direction)
             else

                ! Take a trial step along the search direction
                call sparse_embed_array_copy(trial_aux%kern, aux%kern)
                call sparse_embed_array_axpy(trial_aux%kern, &
                     con_direction, trial_step)
                call kernel_workspace_invalidate(trial_aux)

             end if

             ! cks: %%%%%%%%%%%%%%%%%%%%% ENERGY AT TRIAL LENGTH %%%%%%%%%%%%%%
             total_energy_at_trial_length = internal_energy(trial_aux,.false.)

             mu_term = internal_mu_term(trial_aux%kern)

             Qtrial = total_energy_at_trial_length + mu_term

             if ( (pub_output_detail >= VERBOSE) .and. pub_on_root) then
                  write(stdout,'(a,f20.12)') '     Trial step: ', trial_step
                  write(stdout,'(a,f24.12)') '             Q1: ', Qtrial
             end if
             ! cks: %%%%%%%%%%%%%%%%%%% END ENERGY AT TRIAL LENGTH %%%%%%%%%%%%


             ! Calculate the optimal step by fitting a parabola to the two
             ! function evaluations and initial slope
             optimal_step = internal_quadratic_step()


             ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
             ! If the predicted step is negative, the parabolic fit is poor,
             !     so take a second trial step and fit a cubic instead
             ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

             if (optimal_step < 0.0_DP) then

                if ((pub_output_detail >= VERBOSE) .and. pub_on_root) &
                     write(stdout,'(a)')'         Quadratic fit unsuccessful: &
                     &proceeding to cubic fit'

                optimal_step = 2.0_DP * trial_step


                ! cks: %%%%%%%%%%%%%%%%%%% ENERGY AT TRIAL LENGTH 2 %%%%%%%%%%%
                if (pub_lnv_check_trial_steps) then
                   call internal_kernel_protected_step(trial_aux,optimal_step, &
                        aux,con_direction)
                else
                   ! Take a second trial step along the search direction
                   call sparse_embed_array_copy(trial_aux%kern, aux%kern)
                   call sparse_embed_array_axpy(trial_aux%kern, con_direction, &
                        optimal_step)
                   call kernel_workspace_invalidate(trial_aux)
                end if

                ! Calculate the energy
                total_energy_at_trial_length = internal_energy(trial_aux,.false.)

                mu_term = internal_mu_term(trial_aux%kern)

                Qoptimal = total_energy_at_trial_length + mu_term

                ! cks: %%%%%%%%%%%%%%%% END ENERGY AT TRIAL LENGTH 2 %%%%%%%%%%

                if ( (pub_output_detail >= VERBOSE) .and. pub_on_root) &
                     write(stdout,'(a,f24.12)') '             Q2: ',Qoptimal


                ! Recalculate the optimal step by fitting a cubic to the three
                ! function evaluations and initial slope

                optimal_step = internal_cubic_step()

             else

                ! ndmh: re-ordered next section a bit for improved readability
                ! ndmh: of warnings and supression of warnings at BRIEF level
                if ((pub_output_detail >= VERBOSE) .and. pub_on_root) &
                     write(stdout,'(a,f24.12)') '           Qmin: ',Qpredict
                if (optimal_step > pub_lnv_cg_max_step) then

                   ! cks: prevent exceedingly long quadratic steps
                   if (pub_on_root .and. (pub_output_detail >= VERBOSE)) then
                      write(stdout,'(a)') 'WARNING in &
                           &lnv_denskernel_optimise_cg: setting quadratic &
                           &optimal_step to safe'
                      write(stdout,'(a,f16.12,a)') 'value (Calculated optimal &
                           &quadratic step=',optimal_step,')'
                   end if
                   cg_coeff = 0.0_DP
                   cg_count = 0
                   optimal_step = 0.15_DP
                elseif ( (pub_output_detail >= VERBOSE) .and. pub_on_root) then
                   ! cks: If quadratic step accepted, print verbose info
                   write(stdout,'(a,f16.12,a)') 'Quadratic minimum step length &
                        &(',optimal_step,') accepted'
                end if

             end if


             ! Update trial step for next line minimisation
             if (optimal_step > 0.0_DP) &
                  trial_step = max(sqrt(trial_step*optimal_step),eps)

             line_search_step = optimal_step

          else

             line_search_step = 0.1_DP

          end if

          ! Update the auxiliary matrix
          if (pub_lnv_check_trial_steps) then
              optimal_step = line_search_step
              call internal_kernel_protected_step(trial_aux, line_search_step, &
                   aux,con_direction)
              call sparse_embed_array_copy(aux%kern, trial_aux%kern)
              ! ndmh: reset CG if line step protection truncated the step, since
              ! ndmh: we have not fully minimised in that search direction
              if (line_search_step < optimal_step) then
                 cg_coeff = 0.0_DP
                 cg_count = 0
              end if
          else
             call sparse_embed_array_axpy(aux%kern, con_direction, line_search_step)
          end if

          call kernel_workspace_invalidate(aux)

          ! cks: ========================== END OF LINE SEARCH ================

       end if

       call internal_purify_and_print

    end do

    ! vv: Print out the number of LNV steps in the log file for the auxiliary dofs
    if ((mix_dkn_type == 'XLD' .or. mix_dkn_type == 'XLI' .or. &
         mix_dkn_type == 'NAIVE' .or. mix_dkn_type == 'XLIS') &
         .and. (md_iter_global > mix_dkn_init_num + 1)&
         .and. (pub_maxit_ngwf_cg == 0)) then

       if (pub_on_root) then
          write(filename,'(a,a)') trim(pub_rootname),'.aux.log'
          io_unit = utils_unit()
          open(unit=io_unit,iostat=io_stat,file=filename,status='UNKNOWN',&
               access='SEQUENTIAL',form='FORMATTED',position='APPEND',&
               action='WRITE')
          call utils_open_unit_check('lnv_mod',filename,&
               io_stat)

          write(io_unit,1) iteration - 1

          close(io_unit,iostat=io_stat)
          call utils_close_unit_check('lnv_mod',filename,io_stat)
       end if
    end if

1  format(4x,I4.1,T25,' <--   LNV_STEPS')

#ifdef GPU_SP_TEST
    GPU_SP=.false.                                                    ! Turn off SP code
#endif

    ! Report convergence failure
    if (pub_on_root .and. (.not. converged) .and. pub_output_detail >= NORMAL) then

       write(stdout,'(a,i4,a)') &
           'Finished density kernel iterations (',current_maxit_lnv, ')'
    endif

    ! Write density kernel to file if requested
    loc_force_no_IO=.false.
    if(present(force_no_IO)) then
       loc_force_no_IO=force_no_IO
    end if

    if(.not.loc_force_no_IO) then
       ! jcap: if we are doing an EMFT follow calculation AND we are
       ! doing the final LNV part of the calculation, change the name of
       ! the dkn file
       if (pub_write_denskern.and.(.not.pub_write_converged_dk_ngwfs)) &
            call restart_kernel_write(aux%kern, write_cond=pub_cond_calculate, &
            write_emft=(pub_emft_follow.and.pub_emft.and.pub_emft_lnv_only))
    end if
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! ** De-Allocate structures for sparse matrices

    ! Destroy CG workspace
    if(pub_lnv_cg_type/='LNV_FLETCHER')then
       call sparse_embed_array_destroy(old_co_gradient)
    endif
    call sparse_embed_array_destroy(con_direction)
    call sparse_embed_array_destroy(con_gradient)
    call sparse_embed_array_destroy(co_gradient)
    call sparse_embed_array_destroy(old_con_direction)

    call kernel_destroy(trial_aux)

    ! Deallocate ks and ksk workspace
    call kernel_workspace_destroy(aux)

    deallocate(trace_array,stat=ierr)
    call utils_dealloc_check('lnv_denskernel_optimise_cg','trace_array',ierr)


    ! Stop the timer
    call timer_clock('lnv_denskernel_optimise_cg', 2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
         &Leaving lnv_denskernel_optimise_cg'

    ! Flush output
    call services_flush

    return

  contains

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    real(kind=DP) function internal_energy(kernel_var,update_ham,spoil_force)

      use constants, only: paw_en_size
      use hamiltonian, only: hamiltonian_dens_dep_matrices
      use kernel, only: kernel_purify
      use rundat, only: pub_spin_fac, pub_eda_scfmi, &
           pub_xc_ke_density_required
      use fragment_scfmi, only: denskern_R, scfmi_construct_nonorth_kernel

      implicit none

      ! Argument
      type(DKERN), intent(inout) :: kernel_var
      logical, intent(in) :: update_ham
      logical,optional, intent(in) :: spoil_force

      ! Local variables
      integer :: is, ik
      real(kind=DP) :: paw_sphere_energies(paw_en_size), trace

      ! Calculate the purified density kernel
      call kernel_purify(pur_denskern%kern, kernel_var, rep%overlap, &
           rep%inv_overlap, rep%n_occ)

      ! Stored versions of ks and ksk no longer valid (except on first call)
      if (.not.update_ham) call kernel_workspace_invalidate(kernel_var)


      ! Scale kernel for revised LNV (store norm_fac for later use in gradient)
      if (pub_exact_lnv) then
         call kernel_rescale(pur_denskern, rep%overlap, rep%n_occ, &
              silent=.true., norm_fac=norm_fac)
      else
         norm_fac(:,:) = 1.0_DP
      end if

      ! ndmh: calculate density dependent energies and matrices
      if (.not.pub_cond_calculate .and. .not.pub_eda_scfmi) then

         ! agrecokpt: call with kpoint to include terms in TB method
         call hamiltonian_dens_dep_matrices(ham, lhxc_fine, internal_energy, &
              lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
              ngwf_basis, hub_proj_basis, hub, pur_denskern%kern, mdl,&
              hfxstate, update_ham, lhxc_fixed=.false.,spoil_force=spoil_force,&
              dfdtau_fine = dfdtau_fine, kpt=loc_kpt)

      ! mjsp: for SCF MI calculation, calculate the energy using the
      ! antisymmetry corrected representation of the density kernel
      else if (pub_eda_scfmi) then

         ! mjsp: update denskern_R: proper representation of the density kernel in the
         ! full overlap
         call scfmi_construct_nonorth_kernel(pur_denskern,rep)

         ! mjsp: calculate the density dependent energies and matrices
         ! using the denskern_R kernel
         ! (NOTE: rep%overlap is block diagonal, but this quantity is not used
         ! within hamiltonian_dens_dep_matrices so we do not need to
         ! update it to the full matrix here.)
         ! agrecokpt
         call hamiltonian_dens_dep_matrices(ham, lhxc_fine, internal_energy, &
              lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
              ngwf_basis, hub_proj_basis, hub, denskern_R%kern, mdl, hfxstate, &
              update_ham, lhxc_fixed=.false.,kpt=loc_kpt, &
              dfdtau_fine = eda_dfdtau_fine)

      else if (pub_cond_calculate) then

         ! lr408: In a conduction calculation, there is no need to recalculate
         ! lr408: any matrices, we just want the trace of the kernel with the
         ! lr408: projected Hamiltonian
         ! KPOINTS_DANGER
         internal_energy = 0.0_DP
         do ik = 1, pur_denskern%kern%num_kpoints
            do is = 1, pur_denskern%kern%num_spins
               call sparse_embed_trace(trace, ham%ham(is), pur_denskern%kern%m(is,ik))
               internal_energy = internal_energy + pub_spin_fac * trace
            end do
         end do

      end if


    end function internal_energy

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    real(kind=DP) function internal_gradient_norm()

      use constants, only: SAFE_DIV_EPS
      use kernel, only: kernel_commutator
      use rundat, only: pub_old_lnv, pub_rms_kernel_measure
      use sparse_embed, only: sparse_embed_product, sparse_embed_rms_element, &
           sparse_embed_array2mat, sparse_embed_create, sparse_embed_num_element

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
      call timer_clock('lnv_gradient_norm',1)

      ! Calculate covariant gradient (no tensor correction)
      call internal_co_gradient

      ! Calculate change in mu resulting from truncation for the
      ! covariant gradient
      if (.not. (pub_old_lnv .or. pub_exact_lnv)) then
         call internal_fix_co_direction(co_gradient, mu_extra)
      else
         mu_extra = 0.0_DP
      end if

      ! Apply tensor correction to obtain contravariant gradient
      call internal_con_gradient(mu_extra)

      ! Ensure contravariant gradient preserves electron number constraint
      ! even when matrix multiplication truncation is effective
      if (.not. (pub_old_lnv .or. pub_exact_lnv)) &
           call internal_fix_con_direction(con_gradient)

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
              'WARNING in internal_gradient_norm (lnv_mod.F90): &
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
      call timer_clock('lnv_gradient_norm',2)

    end function internal_gradient_norm

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_co_gradient

      !======================================================================!
      ! This subroutine calculates the gradient of the LNV functional with   !
      ! respect to the auxiliary matrix.                                     !
      !----------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                               !
      ! Based on an original version by Chris-Kriton Skylaris modified by    !
      ! Arash Mostofi and Peter Haynes                                       !
      ! Spin polarised by Peter Haynes, July 2006                            !
      ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009     !
      ! Improvements to reduce number of matrix products required, by        !
      ! Nicholas Hine in November 2010.                                      !
      ! Modified for embedding by Robert Charlton, July 2017.                !
      !======================================================================!

      use kernel, only: kernel_validate_ks
      use rundat, only: pub_spin_fac
      use sparse_embed, only: sparse_embed_array_transpose, sparse_embed_product
      implicit none

      ! Local variables
      integer :: is
      type(SPAM3_EMBED_ARRAY) :: sl, slh, hls, lsc

      ! KPOINTS_DANGER

      ! Note: con_gradient is used as workspace in this routine since it is
      !       not used until this routine has completed

      ! Create sparse matrix structures
      call sparse_embed_array_create(sl, rep%overlap, aux%kern)
      call sparse_embed_array_create(lsc, aux%kern, rep%overlap)
      call sparse_embed_array_create(slh, sl, ham%ham(1))
      call sparse_embed_array_create(hls, ham%ham(1), aux%ks)

      call kernel_validate_ks(aux, rep%overlap)

      ! In revised LNV, replace H by H - mu * S
      if (pub_exact_lnv) then
         do is = 1, pub_num_spins
            call sparse_embed_axpy(ham%ham(is), rep%overlap, -mu(is,PUB_1K))
         end do
      end if

      ! Calculate product S.L
      call sparse_embed_array_transpose(sl, aux%ks)

      ! Calculate products S.L.H and H.L.S
      do is = 1, slh%num_spins
         call sparse_embed_product(slh%m(is,PUB_1K), &
              sl%m(is,PUB_1K), ham%ham(is))
      end do
      call sparse_embed_array_transpose(hls, slh)

      ! Calculate H.L.S.(3I - 2L.S) in con_gradient
      call sparse_embed_array_copy(lsc, aux%ks)
      call sparse_embed_array_scale(lsc, -2.0_DP, 3.0_DP)
      call sparse_embed_array_product(con_gradient, hls, lsc)

      ! Calculate (3I - 2S.L).S.L.H in co_gradient
      call sparse_embed_array_transpose(co_gradient, con_gradient)

      ! Add up H.L.S.(3I - 2L.S) and (3I - 2S.L).S.L.H in co_gradient
      call sparse_embed_array_axpy(co_gradient, con_gradient, 1.0_DP)

      ! Calculate product S.L.H.L.S in con_gradient and add
      ! -2 S.L.H.L.S to co_gradient
      call sparse_embed_array_product(con_gradient, slh, aux%ks)
      call sparse_embed_array_axpy(co_gradient, con_gradient,-2.0_DP)

      ! Add -mu*S (For revised LNV, mu = 0 here).
      ! Total grad: (3I-2S.L).S.L.H + H.L.S.(3I-2L.S) - 2S.L.H.L.S - mu*S
      if (.not. pub_exact_lnv) &
           call sparse_embed_array_axpy(co_gradient, rep%overlap, -mu)

      ! Multiply by spin and normalisation factor
      call sparse_embed_array_scale(co_gradient, pub_spin_fac*norm_fac)

      ! In revised LNV, return H - mu * S to just H
      if (pub_exact_lnv) then
         do is = 1, pub_num_spins
            call sparse_embed_axpy(ham%ham(is), rep%overlap, mu(is,PUB_1K))
         end do
      end if

      ! Deallocate workspace
      call sparse_embed_array_destroy(hls)
      call sparse_embed_array_destroy(slh)
      call sparse_embed_array_destroy(lsc)
      call sparse_embed_array_destroy(sl)

    end subroutine internal_co_gradient

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_con_gradient(mu_extra)

      !===================================================================!
      ! This subroutine converts a covariant gradient to a contravariant  !
      ! gradient by pre- and post-multiplying by the inverse overlap      !
      ! matrix. Note that truncation will introduce errors at this stage. !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                            !
      ! Based on an original version by Chris-Kriton Skylaris and         !
      ! modified by Peter Haynes                                          !
      ! Spin polarised by Peter Haynes, July 2006                         !
      ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009  !
      !===================================================================!

      implicit none

      ! Argument: additional mu from truncation
      real(kind=DP), intent(in) :: mu_extra(pub_num_spins, pub_num_kpoints)

      ! Local variables
      type(SPAM3_EMBED_ARRAY) :: temp

      ! Allocate temporary sparse matrix
      call sparse_embed_array_create(temp, co_gradient, rep%inv_overlap)

      ! Calculate con_grad := co_grad + mu_extra S
      call sparse_embed_array_copy(con_gradient, co_gradient)
      call sparse_embed_array_axpy(con_gradient, rep%overlap, mu_extra)

      ! Calculate temp := con_grad . S^-1
      call sparse_embed_array_product(temp, con_gradient, rep%inv_overlap)

      ! Calculate S^-1 . (co_grad + mu_extra S) . S^-1 =
      !                        = S^-1 . temp -> con_grad
      call sparse_embed_array_product(con_gradient, rep%inv_overlap, temp)

      ! Deallocate workspace
      call sparse_embed_array_destroy(temp)

    end subroutine internal_con_gradient

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_fix_co_direction(co_dir, mu_extra)

      !===================================================================!
      ! This subroutine ensures that a covariant search direction         !
      ! preserves the electron number constraint (to first order).        !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                            !
      ! Spin polarised by Peter Haynes, July 2006                         !
      ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009  !
      !===================================================================!

      implicit none

      ! Argument
      type(SPAM3_EMBED_ARRAY), intent(inout) :: co_dir  ! co direction
      real(kind=DP), intent(out) :: mu_extra(:,:)  ! additional mu required

      ! Local variables
      integer :: is
      real(kind=DP) :: trsinvs     ! tr[S.S^-1]
      real(kind=DP) :: trdinvs(co_dir%num_spins, co_dir%num_kpoints) ! tr[D.S^-1]

      ! Calculate tr[D.S^-1] and tr[S.S^-1]
      call sparse_embed_trace(trsinvs, rep%overlap, rep%inv_overlap)

      call sparse_embed_array_trace(trdinvs, co_dir, rep%inv_overlap)

      ! Calculate and report correction to mu
      mu_extra = - trdinvs / trsinvs

      do is = 1, pub_num_spins
         if (pub_on_root .and. (pub_output_detail >= VERBOSE)) &
              write(stdout,'(a,i1,a,f20.12)') '     fix cov mu', &
              is,': ', mu_extra(is, PUB_1K)
      end do


    end subroutine internal_fix_co_direction

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_fix_con_direction(con_dir)

      !===================================================================!
      ! This subroutine ensures that a contravariant search direction     !
      ! preserves the electron number constraint (to first order).        !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                            !
      ! Spin polarised by Peter Haynes, July 2006                         !
      ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009  !
      !===================================================================!

      implicit none

      ! Argument
      type(SPAM3_EMBED_ARRAY), intent(inout) :: con_dir

      ! Local variables
      integer :: is
      real(kind=DP), dimension(con_dir%num_spins, con_dir%num_kpoints) :: &
           mu_extra, &           ! additional mu required
           trds                  ! tr[D.S]
      real(kind=DP) :: trsinvs   ! tr[S^-1.S]

      ! Calculate tr[D.S] and tr[S.S^-1]
      trds = 0.0_DP
      call sparse_embed_trace(trsinvs, rep%inv_overlap, rep%overlap)

      call sparse_embed_array_trace(trds, con_dir, rep%overlap)

      ! Calculate and report correction to mu
      mu_extra = - trds / trsinvs

      do is = 1, pub_num_spins
         if (pub_on_root .and. (pub_output_detail >= VERBOSE)) &
              write(stdout,'(a,i1,a,f20.12)') '     fix con mu', &
              is,': ', mu_extra(is,PUB_1K)

      end do

      ! Correct gradient
      call sparse_embed_array_axpy(con_dir, rep%inv_overlap, mu_extra)

    end subroutine internal_fix_con_direction

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    logical function internal_test_convergence()

      implicit none

      ! Convergence is assumed when either the RMS gradient or RMS commutator
      ! falls below the threshold

      if ((rms_gradient < lnv_threshold .and. iteration > 0) .or. &
           (commutator < lnv_threshold .and. iteration > 1)) then

         if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then

            if (iteration == 1) then
               write(stdout,'(i5,es27.12,2es15.3)') iteration, &
                    total_energy, commutator, rms_gradient
            else
               write(stdout,'(i5,es27.12,3es15.3)') iteration, total_energy, &
                    rms_gradient, commutator, total_energy - old_total_energy
            endif

            write(stdout,'(/2a)') repeat(' ',17), repeat('.',46)
            write(stdout,'(a,f19.14,a)')'                 | RMS LNV GRADIENT&
                 &= ', rms_gradient, '      |'
            write(stdout,'(a)') '                 | LNV density kernel &
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

    real(kind=DP) function internal_mu_term(trial_kern)
      !===================================================================!
      ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009  !
      !===================================================================!

      use rundat, only: pub_spin_fac, pub_old_lnv

      implicit none

      ! Argument
      type(SPAM3_EMBED_ARRAY), intent(in) :: trial_kern

      ! cks: local variable
      integer :: is
      real(kind=DP) :: trks(trial_kern%num_spins, trial_kern%num_kpoints)

      internal_mu_term = 0.0_DP

      if (pub_exact_lnv) return

      ! jme: KPOINTS_DANGER
      ! jme: this algorithm has not been checked for more than 1 k-point

      ! Calculate 2 tr(L.S)

      call sparse_embed_array_trace(trks, trial_kern, rep%overlap)
      trks = trks * pub_spin_fac

      ! Calculate [2 tr(L.S) - Ne]
      trks = trks - pub_spin_fac * real(rep%n_occ, kind=DP)

      if ((.not. pub_old_lnv) .and. pub_on_root) then
         do is = 1, pub_num_spins
            if ( abs(trks(is, PUB_1K)) > 1.0e-6_DP ) then
               write(stdout,'(a,f20.12)') &
                    'WARNING: electron number wrong by ', trks(is, PUB_1K)
            end if
         end do
      end if

      internal_mu_term = - sum (mu(:,PUB_1K) * trks(:,PUB_1K))

      if (pub_on_root .and. (pub_output_detail >= VERBOSE)) &
           write(stdout,'(a,f20.12)') '-mu(2tr(KS)-Ne): ', internal_mu_term

    end function internal_mu_term

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
              'WARNING in lnv_denskernel_optimise_cg: &
              &vanishing previous gradient norm -', &
              '  conjugate gradients coefficient reset'

         cg_coeff = 0.0_DP

      else

         ! Calculate correction factor for Polak-Ribiere formula if necessary
         if (pub_lnv_cg_type == 'LNV_FLETCHER') then
            polak = 0.0_DP
         else if (pub_lnv_cg_type == 'LNV_POLAK') then
            ! KPOINTS_DANGER
            call sparse_embed_array_trace(trace_array, con_gradient, old_co_gradient)
            polak = sum(trace_array)
         else
            if (pub_on_root) write(stdout,'(a/3a)') &
              'WARNING in lnv_denskernel_optimise_cg: &
              &unknown conjugate gradients method','  "',trim(pub_lnv_cg_type), &
              '" - switching to steepest descents'
            polak = new_g_dot_g
         end if

         ! Calculate conjugate gradients coefficient
         cg_coeff = (new_g_dot_g - polak) / old_g_dot_g

         ! Reset CG coefficient if too large
         if (abs(cg_coeff) > 5.0_DP) then
            if (pub_on_root .and. (pub_output_detail >= VERBOSE)) &
                 write(stdout,'(a/a,e10.3,a)') &
                 'WARNING in lnv_denskernel_optimise_cg: ', &
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
      ! pdh: hack to keep step lengths the same for spin polarised
      ! pdh: systems: multiply by a factor of two
      call sparse_embed_array_scale(con_direction, -real(pub_num_spins,kind=DP))
      call sparse_embed_array_axpy(con_direction, old_con_direction, cg_coeff)

      ! Ensure that search direction preserves constraint to first order
      ! pdh: (not necessary for current algorithms)
      !      if (.not. (pub_old_lnv .or. pub_exact_lnv)) &
      !           call internal_fix_con_direction(con_direction)
      if ( (pub_output_detail >= VERBOSE) .and.(.not. pub_exact_lnv)) then
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

      ! Ensure search direction points 'downhill'
      if (Qslope > 0.0_DP) then

         ! First try steepest descent direction i.e. reset conjugate gradients
         if ( (pub_output_detail >= VERBOSE) .and. pub_on_root) &
              write(stdout,'(a/a)') 'WARNING in lnv_denskernel_optimise_cg: &
              &positive line search gradient -', &
              '  resetting conjugate gradients'

         ! KPOINTS_DANGER
         call sparse_embed_array_copy(con_direction, con_gradient)
         call sparse_embed_array_scale(con_direction, -1.0_DP)
         ! KPOINTS_DANGER
         call sparse_embed_array_trace(trace_array, co_gradient, con_direction)
         Qslope = sum(trace_array)

         ! If this is still positive then things are really pear-shaped i.e.
         ! overlap matrix isn't positive definite or something...
         ! Reverse search direction!
        if (Qslope > 0.0_DP) then
           if ( (pub_output_detail >= VERBOSE ) .and. pub_on_root) &
                write(stdout,'(a)') 'WARNING in lnv_denskernel_optimise_cg: &
                &reversing search gradient'
           call sparse_embed_array_scale(con_direction, -1.0_DP)
           Qslope = -Qslope
        end if

     end if

      if ((pub_output_detail >= VERBOSE) .and. pub_on_root) &
           write(stdout,'(a,f20.12)') '             G0: ',Qslope

       ! Store quantities for next iteration
      call sparse_embed_array_copy(old_con_direction, con_direction)
      if (pub_lnv_cg_type == 'LNV_POLAK') then
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
      ! step until they are within the LNV stability range.               !
      !-------------------------------------------------------------------!
      ! Written by Nicholas Hine, November 2009.                          !
      !===================================================================!

      use kernel, only: kernel_rms_err, kernel_occupancy_bounds

      implicit none

      ! Arguments
      type(DKERN), intent(inout) :: trial_kernel
      real(kind=DP), intent(inout)     :: trial_step
      type(DKERN), intent(inout) :: kernel_var
      type(SPAM3_EMBED_ARRAY), intent(in) :: con_dirn

      ! Local Variables
      real(kind=DP), parameter :: minbound = -0.36602540378443864676_DP
      real(kind=DP), parameter :: maxbound = 1.36602540378443864676_DP
      real(kind=DP) :: rms_err, trial_rms_err
      real(kind=DP), &
           dimension(kernel_var%kern%num_spins, kernel_var%kern%num_kpoints) :: &
           min_occ, max_occ
      real(kind=DP), &
           dimension(trial_kernel%kern%num_spins, trial_kernel%kern%num_kpoints) :: &
           trial_min_occ, trial_max_occ
      integer, parameter :: max_trials=8
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
         ! ndmh: the new kernel has occupancies outside stable LNV bounds
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
            trial_step = trial_step * 0.50_DP
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
         write(stdout,'(a,i2)') 'WARNING in lnv_denskernel_optimise_cg: Trial &
              &step reduced ',max_trials
         write(stdout,'(a)') 'times but kernel may still be unstable'
      end if

    end subroutine internal_kernel_protected_step

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    real(kind=DP) function internal_quadratic_step()

      !===================================================================!
      ! This function fits the quadratic form                             !
      !       2                                                           !
      !    a x  + b x + c = Q(x)                                          !
      !                                                                   !
      ! given Q(0) = Qinitial, Q(trial_step) = Qtrial and Q'(0) = Qslope  !
      !                                                                   !
      ! and then finds the minimum position.                              !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                            !
      !===================================================================!

      implicit none

      ! Local variables
      real(kind=DP) :: a,b            ! Coefficients of quadratic equation
      !real(kind=DP) :: qdiff          ! Qtrial - Qinitial qoh: Unused
      real(kind=DP) :: linear_term    ! b * trial_step
      real(kind=DP) :: quadratic_term ! a * trial_step**2

      ! Determine coefficients of quadratic equation

      linear_term = Qslope * trial_step
      quadratic_term = Qtrial - Qinitial - linear_term
      a = quadratic_term / (trial_step * trial_step)
      b = Qslope

      if ((pub_output_detail >= VERBOSE) .and. pub_on_root) then
         write(stdout,'(a,f20.12)') ' Quadratic coef: ',quadratic_term
         write(stdout,'(a,f20.12)') ' Linear coef   : ',linear_term
      end if


      ! Find minimum position x = -b / (2 a)
      if (abs(a) > abs(Qinitial)*eps*eps) then

         internal_quadratic_step = -0.5_DP * b / a

         Qpredict = Qinitial + 0.5_DP * b * internal_quadratic_step ! c-b^2/(4a)

      else

         internal_quadratic_step = -1.0_DP
         if ((pub_output_detail >= VERBOSE) .and. pub_on_root) &
              write(stdout,'(a)') '  Quadratic term too small: &
              &proceeding to cubic fit'

      end if

    end function internal_quadratic_step

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    real(kind=DP) function internal_cubic_step()

      !===================================================================!
      ! This function fits the cubic form                                 !
      !       3      2                                                    !
      !    a x  + b x  + c x + d = Q(x)                                   !
      !                                                                   !
      ! given Q(0) = Qinitial, Q(trial_step) = Qtrial, Q'(0) = Qslope     !
      ! and now Q(optimal_step) = Qoptimal                                !
      !                                                                   !
      ! and then finds the minimum position.                              !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                            !
      !===================================================================!

      implicit none

      ! Local variables
      real(kind=DP) :: x1,x2   ! The trial step lengths
      real(kind=DP) :: a,b,c   ! Coefficients of the cubic equation
      real(kind=DP) :: bon3a   ! b / (3 a)
      real(kind=DP) :: con3a   ! c / (3 a)
      real(kind=DP) :: disc    ! Discriminant
      real(kind=DP) :: q

      ! Introduce shorthand working variables:

      !   x1 and x2 are the points at which the energy has been evaluated

      x1 = trial_step
      x2 = optimal_step

      !                                          2
      ! The derivative of the cubic fit is  3 a x  + b x + c

      a = ((x2*x2*Qtrial - x1*x1*Qoptimal)/(x1-x2) &
           + (x1+x2)*Qinitial + x1*x2*Qslope) / (x1*x1*x2*x2)
      b = ((x2*x2*x2*Qtrial - x1*x1*x1*Qoptimal)/(x2-x1) &
           - (x1*x1+x1*x2+x2*x2)*Qinitial - x1*x2*(x1+x2)*Qslope) &
           / (x1*x1*x2*x2)
      c = Qslope

      ! Solve this quadratic equation

      if (abs(a*x2/b) > eps) then                      ! Avoid division by zero

         bon3a = b / (3.0_DP * a)
         con3a = c / (3.0_DP * a)
         disc = bon3a * bon3a - con3a                  ! Discriminant

         if (disc >= 0.0_DP) then

            q = -(bon3a + sign(sqrt(disc),bon3a))
            if (b < 0.0_DP) then
               internal_cubic_step = q
            else
               internal_cubic_step = con3a / q
            end if

            if (pub_on_root .and. (pub_output_detail >= VERBOSE)) &
                 write(stdout,'(a,f20.12)') '  Cubic LS coeff: ', &
                 internal_cubic_step

            Qminimum = Qinitial + internal_cubic_step * &
                 (c + internal_cubic_step * (b + internal_cubic_step * a))

            if (pub_on_root .and. (pub_output_detail >= VERBOSE)) &
                 write(stdout,'(a,f20.12)') 'Predicted energy: ',Qminimum

            ! Do not allow the step length to become too long

            if ( abs(internal_cubic_step) > pub_lnv_cg_max_step/2.0_DP) then
               if (pub_on_root .and. (pub_output_detail>=VERBOSE)) &
                    write(stdout,'(a)') &
                    'WARNING in lnv_denskernel_optimise_cg: &
                    &setting cubic step to safe value'
               internal_cubic_step = 0.15_DP
            else
               if (pub_on_root .and. (pub_output_detail >= VERBOSE)) &
                    write(stdout,'(a)') 'Cubic minimum step length accepted'
            end if

         else

            if (pub_on_root .and. (pub_output_detail >= VERBOSE)) &
                 write(stdout,'(a/a)') 'WARNING in lnv_denskernel_optimise_cg: &
                 &line search unsuccessful -','  optimal step length set to &
                 &safe value'

            internal_cubic_step = 0.15_DP

         end if

      else

         if (pub_on_root .and. (pub_output_detail >= VERBOSE)) &
              write(stdout,'(a/a)') 'WARNING in lnv_denskernel_optimise_cg: &
              &line search unsuccessful -', '  optimal step length set to &
              &safe value'

         internal_cubic_step = 0.15_DP

      end if

    end function internal_cubic_step

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_purify_and_print

      use kernel, only: kernel_purify
      use rundat, only: pub_spin_fac, pub_old_lnv
      implicit none

      ! Local variables
      integer :: is
      integer :: iter,numiter

      if (iteration < current_maxit_lnv) then

         if (pub_exact_lnv) then

               ! Perform some "adaptive" purifications
               if (rms_gradient > 0.00001_DP) then
                  numiter = 3
               else
                  numiter = 1
               end if

               ! cks: First fix Ne of L
               call kernel_rescale(aux, rep%overlap, rep%n_occ, &
                    can_rescale_ks=.true.)
               ! cks: Now fix the idempotency of L
               call kernel_fix(aux, rep%overlap, rep%inv_overlap, numiter)

         else

            if (pub_old_lnv) then

               ! Perform some ad hoc purifications
               if (rms_gradient > 0.00001_DP) then
                  numiter = 3
               else
                  numiter = 1
               end if
               do iter=1,numiter
                  call kernel_purify(aux%kern, aux, rep%overlap, rep%inv_overlap, &
                       rep%n_occ)
                  call kernel_workspace_invalidate(aux)
               end do

            else

               ! Perform some "adaptive" purifications
               if (rms_gradient > 0.00001_DP) then
                  numiter = 3
               else
                  numiter = 1
               end if

               call kernel_fix(aux, rep%overlap, rep%inv_overlap, numiter)

               ! cks: to ensure the kernel is purified adequately even when
               ! cks: the kernel_fix has exited without doing any iterations
               call kernel_purify(aux%kern, aux, rep%overlap, rep%inv_overlap, &
                    rep%n_occ)
               ! aux has changed => invalidate ks and ksk
               call kernel_workspace_invalidate(aux)

               ! cks: The "unpurified" denskern must have the correct Ne.
               ! cks: Enforce this in the face of truncation.
               call kernel_rescale(aux, rep%overlap, rep%n_occ)

            end if

         end if

         ! Calculate purified density kernel
         call kernel_purify(pur_denskern%kern, aux, rep%overlap, rep%inv_overlap, &
                            rep%n_occ)
         ! ks and ksk are still valid for aux after this

         ! Calculate electron number from purified kernel
         if (pub_exact_lnv) then
            ne_pure = sum(rep%n_occ)
         else
            ! KPOINTS_DANGER (this sums over k-points)
            call sparse_embed_array_trace(trace_array, pur_denskern%kern, rep%overlap)
            ne_pure = sum(trace_array)
            ne_pure = ne_pure * pub_spin_fac
         end if

      else

         ! Last iteration
         ! KPOINTS_DANGER (this sums over k-points)
         call sparse_embed_array_trace(trace_array, pur_denskern%kern, rep%overlap)
         ne_pure = sum(trace_array)
         ne_pure = ne_pure * pub_spin_fac

         cg_coeff = 0.0_DP
         line_search_step = 0.0_DP

      end if

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

    end subroutine internal_purify_and_print

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

  end subroutine lnv_denskernel_optimise_cg

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////


  subroutine lnv_calculate_mu(mu, aux, ham, rep, ngwf_basis)
      !===================================================================!
      ! Author unknown.                                                   !
      ! @docme                                                            !
      !-------------------------------------------------------------------!
      ! Rudimentary notes on operation by Jacek Dziedzic, June 2017:      !
      ! Calculates mu, a chemical-potential-like auxiliary quantity.      !
      ! In the default mode (pub_exact_lnv), this is                          !
      ! mu = tr[LHLS(3I-2LS)] / tr[LSLS(3I-2LS)] = tr[KH]/tr[KS].         !
      ! cf. Alvaro's thesis, eq. (2.125) and around.                      !
      ! 'rep' is only needed for its %overlap. 'aux' is inout because of  !
      ! its workspace (storing LS). 'ngwf_basis' is *not* needed in the   !
      ! default mode.                                                     !
      ! Fixed intent for 'ham' argument, Jacek Dziedzic, June 2017.       !
      !===================================================================!

    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use kernel, only: DKERN, kernel_workspace_invalidate, &
         kernel_workspace_create, &
         kernel_workspace_destroy
    use ngwf_representation, only: NGWF_HAM, NGWF_REP
    use rundat, only: pub_spin_fac, pub_exact_lnv, pub_debug_on_root, &
         pub_num_kpoints, PUB_1K
    use sparse_embed, only: sparse_embed_product, sparse_embed_array_create, &
         sparse_embed_array_destroy, sparse_embed_array_axpy, &
         sparse_embed_array_copy, sparse_embed_array_product, &
         sparse_embed_array_trace, sparse_embed_array_scale, SPAM3_EMBED_ARRAY
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: mu(:,:)    ! Chemical potential
    type(DKERN), intent(inout) :: aux
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(NGWF_HAM), intent(inout) :: ham
    type(NGWF_REP), intent(in) :: rep

    ! Local variables
    logical :: deallocate_workspace
    real(kind=DP), allocatable :: trace_array(:,:)
    integer :: ierr

    if (pub_debug_on_root) then
       write(stdout,'(a)') 'DEBUG: Entering lnv_calculate_mu'
    end if

    ! jme: KPOINTS_DANGER
    ! The algorithms in this module have not been checked for more than 1 k-point
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine lnv_calculate_mu not checked yet for more&
         & than one k-point.')

    ! Check for workspace
    deallocate_workspace = .false.
    if (.not.aux%workspace_created) then
       call kernel_workspace_create(aux, rep%overlap)
       deallocate_workspace = .true.
    end if

    allocate(trace_array(aux%ks%num_spins,aux%ks%num_kpoints),stat=ierr)
    call utils_alloc_check('lnv_calculate_mu','trace_array',ierr)

    ! Calculate mu
    if (pub_exact_lnv) then
       call internal_lnv_mu_ideal()
    else
       call internal_ms_mu_ideal()
    end if

    if (deallocate_workspace) call kernel_workspace_destroy(aux)

    deallocate(trace_array,stat=ierr)
    call utils_dealloc_check('lnv_calculate_mu','trace_array',ierr)

    if (pub_debug_on_root) then
       write(stdout,'(a)') 'DEBUG: Leaving lnv_calculate_mu'
    end if

contains

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    real(kind=DP) function internal_lnv_mu_old()

      !===================================================================!
      ! This function calculates mu, the chemical potential-like Lagrange !
      ! multiplier which preserves the electron number to first order in  !
      ! the original version of the LNV.  If there is no                  !
      ! truncation applied to the gradient, then this would be:           !
      !            mu = 6 * tr[LSLS(I-LS)(I-LS)]                          !
      ! which turns out to be proportional to the penalty functional.     !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                            !
      !===================================================================!

      use kernel, only: kernel_validate_ks
      implicit none

      ! Local variables
      integer :: is
      real(kind=DP) :: trace
      type(SPAM3_EMBED_ARRAY) :: lslsc, lsc

      ! Create sparse matrix structures
      call sparse_embed_array_create(lsc, aux%ks)
      call sparse_embed_array_create(lslsc, aux%ks, lsc)

      call kernel_validate_ks(aux, rep%overlap)

      ! Calculate (I - L.S)
      call sparse_embed_array_copy(lsc, aux%ks)
      call sparse_embed_array_scale(lsc, -1.0_DP, 1.0_DP)

      ! Calculate L.S.(I - L.S)
      call sparse_embed_array_product(lslsc, aux%ks, lsc)

      ! Calculate tr[L.S.L.S.(I-L.S).(I-L.S)]
      ! KPOINTS_DANGER: the following sum also sums over k-points!!
      call sparse_embed_array_trace(trace_array, lslsc, lslsc)
      trace = sum(trace_array)

      ! Function result
      internal_lnv_mu_old = -3.0_DP * pub_spin_fac * trace / sum(ngwf_basis(:)%num)

      ! Deallocate workspace
      call sparse_embed_array_destroy(lslsc)
      call sparse_embed_array_destroy(lsc)

    end function internal_lnv_mu_old

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_lnv_mu_ideal()

      !===================================================================!
      ! This function calculates mu, the chemical potential-like Lagrange !
      ! multiplier which is used in the revised LNV and is simply given   !
      ! by:                                                               !
      !        mu = tr[LHLS(3I-2LS)] / tr[LSLS(3I-2LS)]                   !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, December 2004                            !
      ! Modified for spin polarisation by Peter Haynes, July 2006         !
      ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009  !
      ! Modified for embedding by Robert Charlton, August 2017.           !
      !===================================================================!

      use kernel, only: kernel_validate_ks
      implicit none

      ! Local variables
      integer :: is, ik
      real(kind=DP), &
           dimension(aux%kern%num_spins, aux%kern%num_kpoints) :: trh, trs
      type(SPAM3_EMBED_ARRAY) :: lslsc, lh, lsc

      ! Create sparse matrix structures
      call sparse_embed_array_create(lsc, aux%ks)
      call sparse_embed_array_create(lh, aux%kern, ham%ham(1))
      call sparse_embed_array_create(lslsc, aux%ks, lsc)

      call kernel_validate_ks(aux, rep%overlap)

      ! Calculate products L.H, L.S and (3I - 2L.S)
      do ik = 1, aux%kern%num_kpoints
         do is = 1, aux%kern%num_spins
            call sparse_embed_product(lh%m(is,ik), aux%kern%m(is,ik), ham%ham(is))
         end do
      end do
      call sparse_embed_array_copy(lsc, aux%ks)
      call sparse_embed_array_scale(lsc, -2.0_DP, 3.0_DP)

      ! Calculate L.S.(3I - 2L.S)
      call sparse_embed_array_product(lslsc, aux%ks, lsc)

      ! Calculate traces
      call sparse_embed_array_trace(trs, lslsc, aux%ks)
      call sparse_embed_array_trace(trh, lslsc, lh)

      where (abs(trs) > tiny(1.0_DP))
         mu = trh / trs
      else where
         mu = 0.0_DP
      end where

      ! Deallocate workspace
      call sparse_embed_array_destroy(lslsc)
      call sparse_embed_array_destroy(lh)
      call sparse_embed_array_destroy(lsc)

    end subroutine internal_lnv_mu_ideal

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_ms_mu_ideal()

      !===================================================================!
      ! This function calculates mu, the chemical potential-like Lagrange !
      ! multiplier which makes the contravariant gradient traceless in    !
      ! the Milliam-Scuseria variant of the LNV.  If there is no          !
      ! truncation applied to the gradient, then this would be:           !
      !            mu = 6 * tr[H(L-L.S.L)]/num which                      !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                            !
      ! Based on a version by Chris-Kriton Skylaris, 2001                 !
      ! Modified for spin polarisation by Peter Haynes, July 2006         !
      ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009  !
      !===================================================================!

      use kernel, only: kernel_validate_ks, kernel_validate_ksk
      use sparse_embed, only: sparse_embed_trace
      implicit none

      ! Local variables
      integer :: is, ik
      real(kind=DP), &
           dimension(aux%kern%num_spins, aux%kern%num_kpoints) :: trace
      type(SPAM3_EMBED_ARRAY) :: lslc

      ! Create sparse matrix structures
      call sparse_embed_array_create(lslc, aux%ks, aux%kern)

      call kernel_validate_ks(aux, rep%overlap)
      call kernel_validate_ksk(aux, rep%overlap)

      ! Calculate L.S.L - L
      call sparse_embed_array_copy(lslc, aux%ksk)
      call sparse_embed_array_axpy(lslc, aux%kern, -1.0_DP)

      ! Calculate tr[(L.S.L - L)H]
      ! KPOINTS_DANGER
      do ik = 1, lslc%num_kpoints
         do is = 1, lslc%num_spins
            call sparse_embed_trace(trace(is, ik), lslc%m(is,ik), ham%ham(is))
         end do
      end do
      mu = -3.0_DP * pub_spin_fac * trace / sum(ngwf_basis(:)%num)

      ! Deallocate workspace
      call sparse_embed_array_destroy(lslc)

    end subroutine internal_ms_mu_ideal

  end subroutine lnv_calculate_mu

end module lnv
