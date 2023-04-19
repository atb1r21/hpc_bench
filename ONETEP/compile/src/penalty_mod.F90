! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                     Penalty functional module                  !
!                                                                !
! This module implements an alternative to the LNV method for    !
! optimising the density kernel K for a fixed set of NGWFs,      !
! observing the idempotency constraint approximately by means of !
! a penalty functional related to McWeeny's purifying            !
! transformation.                                                !
!----------------------------------------------------------------!
! Written by Peter Haynes, 14/02/06                              !
! Based on an earlier version by Peter Haynes and Chris-Kriton   !
!   Skylaris.                                                    !
! Adapted to use SPAM3 matrices by Nicholas Hine, July 2009      !
!================================================================!

module penalty

  use rundat, only: pub_debug_on_root

  implicit none

  private

  public :: penalty_denskernel_optimise_cg
  public :: pen_reset_trial_step_to_initial

  logical, save :: pen_reset_trial_step_to_initial = .false.

contains

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine penalty_denskernel_optimise_cg(denskern, rep, ham, &
       lhxc_fine, total_energy, ngwf_basis, hub_proj_basis, hub, mdl, hfxstate,&
       pen_threshold)

    !==========================================================================
    ! This subroutine optimises the density kernel for a fixed set of NGWFs
    ! using the penalty functional to impose idempotency approximately.
    !
    ! See P.D. Haynes and M.C. Payne, Phys. Rev. B 59, 12173 (1999).
    !
    !==========================================================================
    ! Arguments:
    ! denskern (inout)        : Density kernel of DKERN type.
    ! rep (input)             : NGWF representation (functions and matrices)
    ! lhxc_fine (output)      : Local-Hartree-Exhange-Correlation potential
    ! total_energy (output)   : Total energy
    ! ngwf_basis (input)      : Function basis of NGWFs
    ! pen_threshold (input)   : Convergence threshold for LNV gradient
    !==========================================================================
    ! Based on lnv_mod.F90 written by Chris-Kriton Skylaris.
    ! Original version written by Peter Haynes.
    ! Parallelised by Peter Haynes and Chris-Kriton Skylaris in December 2003.
    ! Rewritten to use new sparse matrix format by Peter Haynes, February 2006.
    ! Modified to use parallel SPAM 2 throughout by Peter Haynes, July 2006
    ! Modified to include Nonlinear Core Corrections by Nicholas Hine in
    ! January 2009
    ! DFT+U added by David O'Regan, April 2009
    ! Adapted to use SPAM3 matrices by Nicholas Hine, July 2009
    ! Tiny fix by Jacek Dziedzic to reset trial step for every new calc
    ! and not once per ONETEP run, Oct 2013.
    ! Adapted for SPAM3_EMBED by Robert Charlton, August 2018.
    !==========================================================================

    use comms, only: pub_on_root, pub_root_proc_id, comms_bcast
    use constants, only: DP, VERBOSE, stdout, NORMAL, max_spins
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_build_matrix
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: DKERN, kernel_fix, kernel_copy, &
         kernel_workspace_invalidate, kernel_workspace_create, &
         kernel_workspace_destroy, kernel_normalise, &
         kernel_occupancy_bounds, kernel_create, kernel_destroy
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use restart, only: restart_kernel_write
    use rundat, only: pub_maxit_pen, pub_pen_param, pub_output_detail, &
         pub_write_denskern, pub_write_converged_dk_ngwfs, &
         pub_num_spins, pub_num_kpoints, PUB_1K, &
         pub_xc_ke_density_required
    use services, only: services_flush
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_axpy, sparse_embed_copy, &
         sparse_embed_destroy, SPAM3_EMBED_ARRAY, sparse_embed_array_create, &
         sparse_embed_array_axpy, sparse_embed_array_copy, sparse_embed_array_product, &
         sparse_embed_array_scale, sparse_embed_array_destroy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_banner

    implicit none

    ! Arguments: input/output
    type(DKERN), intent(inout) :: denskern
    type(NGWF_HAM), intent(inout) :: ham
    type(MODEL), intent(in) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    real(kind=DP), intent(out) :: lhxc_fine(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins,mdl%nsub)
    real(kind=DP), intent(out) :: total_energy
    type(HUBBARD_MODEL), intent(inout) :: hub

    ! Arguments: input only
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in) :: hub_proj_basis(:)
    type(NGWF_REP), intent(in) :: rep
    real(kind=DP), intent(in) :: pen_threshold

    ! Local variables
    real(kind=DP) :: commutator
    real(kind=DP) :: penalty_energy
    real(kind=DP) :: rms_gradient,line_search_coeff,cg_coeff,eps
    real(kind=DP) :: total_energy_at_trial_step,lhxc_energy
    real(kind=DP), dimension(pub_num_spins, pub_num_kpoints) :: &
         max_occ, min_occ
    real(kind=DP), allocatable :: trace_array(:,:), trace_array2(:,:)
    integer       :: is, ierr
    integer       :: iteration,cg_count,cg_max
    logical       :: converged
    ! agrecocmplx
    logical :: loc_cmplx

    ! Conjugate gradient variables
    type(SPAM3_EMBED_ARRAY) :: co_gradient       ! Covariant gradient
    type(SPAM3_EMBED_ARRAY) :: con_gradient      ! Contravariant gradient
    type(SPAM3_EMBED_ARRAY) :: con_direction     ! Contra search direction
    type(SPAM3_EMBED_ARRAY) :: old_con_direction ! Prev contra search dirn
    type(SPAM3_EMBED_ARRAY) :: old_co_gradient   ! Prev covariant gradient
    type(DKERN) :: trial_denskern    ! Trial density kernel

    !--------------------------------------------------------------------------
    ! pdh: variables introduced purely for line minimisation
    ! pdh: the function to be minimised is generally called 'Q'
    !--------------------------------------------------------------------------
    real(kind=DP) :: Qinitial          ! initial value of function
    real(kind=DP) :: Qslope            ! derivative along search direction
    real(kind=DP),save :: trial_step   ! trial step length
    real(kind=DP) :: Qtrial            ! value of function at trial step
    real(kind=DP) :: optimal_step      ! optimal step length
    real(kind=DP) :: Qoptimal          ! value of function at optimal step
    real(kind=DP) :: Qpredict          ! predicted value at optimal step
    real(kind=DP) :: minimum_step      ! minimum step predicted by cubic
    real(kind=DP) :: Qminimum          ! predicted value of function at minimum
    real(kind=DP) :: linear_term
    real(kind=DP) :: quadratic_term
    real(kind=DP) :: xx, yy, aa, bb, aon3b, disc
    logical       :: trial2
    !--------------------------------------------------------------------------
    ! pdh: end variables for line minimisation
    !--------------------------------------------------------------------------

    ! cks : CONVENTIONS:
    ! cks : search directions are ALWAYS CONTRAVARIANT.
    ! cks : gradients can be either covariant or contravariant.
    ! cks : their name will indicate this in each case.

    ! pdh : All density kernels are for spin alpha only, i.e. all spin
    !       degeneracy factors are explicitly included in this module

    ! JCW: Abort if KE density required (meta-GGA + penalty functional
    ! JCW: not implemented).
    call utils_assert(.not.pub_xc_ke_density_required, "Error in &
         &penalty_denskernel_optimise_cg: pub_xc_ke_density_required is true, &
         &but combination of KE-density-dependent XC functionals with penalty &
         &functional has not been implemented/tested.")

    ! cks: find time spent in this subroutine
    call timer_clock('penalty_denskernel_optimise_cg',1)

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine penalty_denskernel_optimise_cg not ready yet for more&
         & than one k-point.')

    ! agrecocmplx
    loc_cmplx = rep%overlap%p%iscmplx

    if(pen_reset_trial_step_to_initial) then
       trial_step = 1.0_DP
       pen_reset_trial_step_to_initial = .false.
    end if

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &penalty_denskernel_optimise_cg'

    call services_flush

    ! ndmh: allocate kernel work arrays
    call kernel_workspace_create(denskern, rep%overlap)

    ! agrecocmplx
    call kernel_create(trial_denskern, denskern%kern%m(1,1)%structure//rep%postfix, &
         is_cmplx=loc_cmplx)
!    call kernel_workspace_create(trial_denskern, rep%overlap)

    ! Allocate matrices with the sparsity pattern of the density kernel
    call sparse_embed_array_create(co_gradient, denskern%kern)
    call sparse_embed_array_create(con_gradient, denskern%kern)
    call sparse_embed_array_create(con_direction, denskern%kern)
    call sparse_embed_array_create(old_con_direction, denskern%kern)
    call sparse_embed_array_create(old_co_gradient, denskern%kern)

    allocate(trace_array(denskern%kern%num_spins, denskern%kern%num_kpoints), stat=ierr)
    call utils_alloc_check('penalty_denskernel_optimise_cg','trace_array',ierr)
    allocate(trace_array2(denskern%kern%num_spins, denskern%kern%num_kpoints), stat=ierr)
    call utils_alloc_check('penalty_denskernel_optimise_cg','trace_array2',ierr)

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    ! Make sure the density-matrix corresponds to exactly the right
    ! number of electrons
    call kernel_normalise(denskern, rep%overlap, rep%inv_overlap, rep%n_occ)

    ! Initialise conjugate gradients
    rms_gradient = 0.0_DP
    cg_count = 0
    cg_max = 5
    eps = epsilon(1.0_DP)
    converged = .false.

    if (pub_on_root) then
       write(stdout,'(a)') repeat('.',80)
       write(stdout,'(a)') utils_banner('<', &
            'Penalty functional density kernel optimisation')
       write(stdout,'(a)') repeat('~',80)

       write(stdout,'(a)')'iter|       energy      |   rms grad  |&
            &  commutator |LScoef|CGcoef|      Ne    |'
    endif

    do iteration=1,pub_maxit_pen

       if (pub_debug_on_root) write(stdout,'(2(a,i3))') 'DEBUG: In &
            &penalty_denskernel_optimise_cg: iteration ',iteration,' of ', &
            pub_maxit_pen

       call services_flush

       call kernel_copy(trial_denskern, denskern)

       !=========================== ENERGY AT POINT 0 =========================
       total_energy = internal_energy(trial_denskern, .true.)
       penalty_energy = internal_penalty_functional(trial_denskern)
       Qinitial = total_energy + penalty_energy ! for line minimisation

       if (pub_output_detail >= VERBOSE .and. pub_on_root) &
            write(stdout,'(3(a,f20.12))') '      Q0=E0+P0: ',Qinitial,'=', &
            total_energy,'+',penalty_energy

       if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: In &
            &penalty_denskernel_optimise_cg: &
            &initial functional evaluation complete'

       !======================== END ENERGY AT POINT 0 ========================

       !*************** CALCULATE KOHN-SHAM MATRIX ****************************
       call hamiltonian_build_matrix(ham,rep)

       if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: In &
            &penalty_denskernel_optimise_cg: &
            &Kohn-Sham matrix constructed'

       !*********** END CALCULATE KOHN-SHAM MATRIX ****************************


       ! %%%%%%%%%%%%%%%%%%%%%%%% GRADIENT AT POINT 0 %%%%%%%%%%%%%%%%%%%%%%%%%

       rms_gradient = internal_gradient_norm()
       ! pdh: hack to make consistent with spin-unpolarised case
       rms_gradient = rms_gradient * pub_num_spins


       if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: In &
            &penalty_denskernel_optimise_cg: &
            &initial gradient calculated'

       ! %%%%%%%%%%%%%%%%%%%%%%% END GRADIENT AT POINT 0 %%%%%%%%%%%%%%%%%%%%%%


       ! ########################## TEST CONVERGENCE ##########################

       converged = internal_test_convergence()
       call comms_bcast(pub_root_proc_id,converged)
       if (converged) exit

       ! ######################## END TEST CONVERGENCE ########################


       !============================= LINE SEARCH =============================

       if (pub_output_detail >= VERBOSE .and. pub_on_root) &
            write(stdout,'(a,f20.12)') '    Trial step: ', trial_step

       call internal_search_direction

       if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: In &
            &penalty_denskernel_optimise_cg: &
            &search direction constructed'

       ! Take a trial step along the search direction
       call sparse_embed_array_axpy(trial_denskern%kern, con_direction, trial_step)
       ! ndmh: mark kernel workspace variables as invalid
       call kernel_workspace_invalidate(trial_denskern)

       ! %%%%%%%%%%%%%%%%%%%%% ENERGY AT TRIAL LENGTH %%%%%%%%%%%%%%%%%%%%%%

       total_energy_at_trial_step = internal_energy(trial_denskern,.false.)
       penalty_energy = internal_penalty_functional(trial_denskern)
       Qtrial = total_energy_at_trial_step + penalty_energy ! line minimisation

       if (pub_output_detail >= VERBOSE .and. pub_on_root) &
            write(stdout,'(3(a,f20.12))') '      Q1=E1+P1: ',Qtrial,'=', &
            total_energy_at_trial_step,'+',penalty_energy

       if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: In &
            &penalty_denskernel_optimise_cg: &
            &trial functional evaluation complete'

       ! %%%%%%%%%%%%%%%%%%% END ENERGY AT TRIAL LENGTH %%%%%%%%%%%%%%%%%%%%

       ! Fit a quadratic to find the optimal step

       optimal_step = internal_quadratic_step()

       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !pdh: If the predicted step is negative, the parabolic fit is poor,
       !     so take a second trial step and fit a cubic instead
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

       if (optimal_step < 0.0_DP) then

          if (pub_output_detail >= VERBOSE .and. pub_on_root) &
               write(stdout,'(a)')'    !!! Quadratic fit unsuccessful; &
               &proceeding to cubic fit !!!'

          ! Take a second trial step

          trial2 = .true.
          optimal_step = 2.0_DP * trial_step
          call sparse_embed_array_copy(trial_denskern%kern, denskern%kern)
          call sparse_embed_array_axpy(trial_denskern%kern, con_direction, optimal_step)
          ! ndmh: mark kernel workspace variables as invalid
          call kernel_workspace_invalidate(trial_denskern)

       ! ndmh: better protection against dodgy optimal steps
       else if ((optimal_step < trial_step*0.01_DP) .or. &
            (optimal_step > trial_step*20.0_DP)) then

          if (pub_output_detail >= VERBOSE .and. pub_on_root) &
               write(stdout,'(a)')'    !!! Quadratic fit may be unreliable; &
               &taking second trial step !!!'

          ! Take a second trial step

          trial2 = .true.
          call sparse_embed_array_copy(trial_denskern%kern, denskern%kern)
          call sparse_embed_array_axpy(trial_denskern%kern, con_direction, optimal_step)
          ! ndmh: mark kernel workspace variables as invalid
          call kernel_workspace_invalidate(trial_denskern)


       else

          trial2 = .false.
          if (pub_output_detail >= VERBOSE .and. pub_on_root) &
               write(stdout,'(a,f20.12)') ' Predicted min: ',Qpredict

       end if

       if (trial2) then ! if a second trial step is necessary

          ! %%%%%%%%%%%%%%%%%%% ENERGY AT TRIAL LENGTH 2 %%%%%%%%%%%%%%%%%%%

          total_energy_at_trial_step = internal_energy(trial_denskern,.false.)
          penalty_energy = internal_penalty_functional(trial_denskern)
          Qoptimal = total_energy_at_trial_step + penalty_energy

          if (pub_output_detail >= VERBOSE .and. pub_on_root) &
               write(stdout,'(3(a,f20.12))') '      Q2=E2+P2: ',Qoptimal,'=', &
               total_energy_at_trial_step,'+',penalty_energy

          if (pub_debug_on_root) then
             write(stdout,'(a)') 'DEBUG: In penalty_denskernel_optimise_cg: &
                  &second trial functional evaluation complete'
          end if

          ! %%%%%%%%%%%%%%%% END ENERGY AT TRIAL LENGTH 2 %%%%%%%%%%%%%%%%%%

          ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          ! pdh: Now have sufficient information to fit a cubic
          ! ndmh: (as long as we have not already decided to re-do trial step
          ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

          if (optimal_step == trial_step*2.0_DP) then
             optimal_step = internal_cubic_step()
          else
             ! ndmh: find new optimal step, ignoring first trial step
             Qtrial = Qoptimal
             trial_step = optimal_step
             optimal_step = internal_quadratic_step()
          end if

       else

          ! ndmh: prevent exceedingly long quadratic steps
          if (optimal_step > 3.0_DP) then
             if (pub_on_root .and. (pub_output_detail >= NORMAL)) then
                write(stdout,'(a,f16.12)') &
                    'Calculated optimal quadratic step=',optimal_step
                write(stdout,'(a/a)') &
                    'WARNING in penalty_denskernel_optimise_cg: ', &
                    'setting quadratic optimal_step to safe value'
             end if
             optimal_step = 0.15_DP
             cg_coeff = 0.0_DP
             cg_count = 0
          else if ( (pub_output_detail >= VERBOSE) .and. pub_on_root) then
             ! cks: If quadratic step accepted, print verbose info
             write(stdout,'(a)') '    !!! Quadratic minimum step &
                 &length accepted !!!'

          end if

       end if

       ! Update trial step for next line minimisation

       if (optimal_step > 0.0_DP) &
            trial_step = max(sqrt(trial_step*optimal_step),epsilon(1.0_DP))

       line_search_coeff = optimal_step

       ! Set the new density kernel
       call sparse_embed_array_axpy(denskern%kern, con_direction, line_search_coeff)
       ! ndmh: mark kernel workspace variables as invalid
       call kernel_workspace_invalidate(denskern)

       if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: In &
            &penalty_denskernel_optimise_cg: &
            &optimal step taken'

       ! Store gradients for determination of next CG coefficient
       call sparse_embed_array_copy(old_con_direction, con_direction)
       call sparse_embed_array_copy(old_co_gradient, co_gradient)

       ! ========================== END OF LINE SEARCH ========================

       call internal_print()

    end do

    ! Print out information concerning occupation number errors etc

    total_energy = internal_energy(denskern,.false.)
    penalty_energy = internal_penalty_functional(denskern)
    call kernel_occupancy_bounds(max_occ(:,:), min_occ(:,:), denskern, rep%overlap)
    if (pub_on_root) then
       write(stdout,'(/a)') '======================== Penalty functional &
            &minimisation ======================='
       write(stdout,'(a,f20.12)') '   Final total energy     : ',total_energy
       write(stdout,'(a,f20.12)') '   Final penalty energy   : ',penalty_energy
       write(stdout,'(a,f14.6)') '   RMS occupancy error    : ',&
            sqrt(penalty_energy / (pub_pen_param*ngwf_basis%num))

       if (pub_num_spins == 1) then
          write(stdout,'(a,2(f6.3,a))') '   Occupancy bounds       : [', &
               min_occ(1,PUB_1K),',',max_occ(1,PUB_1K),']'
       else
          do is=1,pub_num_spins
             write(stdout,'(a,i1,a,2(f6.3,a))') '   Occupancy bounds spin ', &
                  is,': [',min_occ(is,PUB_1K),',',max_occ(is,PUB_1K),']'
          end do
       end if
       write(stdout,'(a/)') '=============================================&
            &==================================='

       if (minval(min_occ(:,PUB_1K)) < 0.5_DP*(1.0_DP-sqrt(5.0_DP)) .or. &
            maxval(max_occ(:,PUB_1K)) > 0.5_DP*(1.0_DP+sqrt(5.0_DP))) then
          write(stdout,'(/a)') &
               'WARNING: one or more occupancies may be outside the stable range!'
          write(stdout,'(a,i5,a,f10.4)') &
               'WARNING: Recommend maxit_pen >= ',2*pub_maxit_pen, &
               ' and penparam >= ',2.0_DP*pub_pen_param
       end if
    end if

    ! Fix occupancies if necessary
    if (minval(min_occ(:,PUB_1K)) < 0.5_DP*(1.0_DP-sqrt(3.0_DP)) .or. &
         maxval(max_occ(:,PUB_1K)) > 0.5_DP*(1.0_DP+sqrt(3.0_DP))) then
       call kernel_fix(denskern, rep%overlap, rep%inv_overlap)
       call kernel_occupancy_bounds(max_occ(:,:), min_occ(:,:), &
            denskern, rep%overlap)
    end if

    if (pub_on_root) then
       if (minval(min_occ(:,PUB_1K)) < 0.5_DP*(1.0_DP-sqrt(3.0_DP)) .or. &
            maxval(max_occ(:,PUB_1K)) > 0.5_DP*(1.0_DP+sqrt(3.0_DP))) then
          write(stdout,'(a)') &
               'WARNING: one or more occupancies remains outside the stable range!'
       end if
    end if

    ! Stop if failed to converge
    if (.not. converged .and. pub_on_root) then
        write(stdout,'(a,i4,a)') 'Finished density kernel iterations (', &
        pub_maxit_pen, ')'
    end if

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! ** Deallocate structures for sparse matrices

    call sparse_embed_array_destroy(old_co_gradient)
    call sparse_embed_array_destroy(con_direction)
    call sparse_embed_array_destroy(con_gradient)
    call sparse_embed_array_destroy(co_gradient)
    call sparse_embed_array_destroy(old_con_direction)

    call kernel_destroy(trial_denskern)

    ! ndmh: deallocate kernel work arrays
    call kernel_workspace_destroy(denskern)

    deallocate(trace_array,stat=ierr)
    call utils_dealloc_check('penalty_denskernel_optimise_cg','trace_array',ierr)
    deallocate(trace_array2,stat=ierr)
    call utils_dealloc_check('penalty_denskernel_optimise_cg','trace_array2',ierr)

    ! cks: output density kernel to file if this is requested
    if (pub_write_denskern.and.(.not.pub_write_converged_dk_ngwfs)) &
         call restart_kernel_write(denskern%kern)

    ! Stop the timer
    call timer_clock('penalty_denskernel_optimise_cg',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &penalty_denskernel_optimise_cg'

    call services_flush

    return

  contains

    !******************************************************
    !******************************************************

    logical function internal_test_convergence()

      ! ########################## TEST CONVERGENCE ##########################

      internal_test_convergence = (rms_gradient < pen_threshold) .and. &
           (iteration > 0)

      if (internal_test_convergence) then
         if (pub_on_root) then
            write(stdout,'(i2,f20.14,2f15.12)') iteration,total_energy,&
                 rms_gradient,commutator
            write(stdout,'(/a)') repeat(' ',17) // repeat('.',46)
            write(stdout,'(a,f19.14,a)')'                 | RMS PF GRADIENT   = ',&
                 rms_gradient,'    |'
            write(stdout,'(a)') '                 | PF density kernel &
                 &optimisation converged!  |'
            write(stdout,'(a)') repeat(' ',17) // repeat('~',46)
         end if
      endif

      ! ######################## END TEST CONVERGENCE ########################

    end function internal_test_convergence


    !******************************************************
    !******************************************************


    real(kind=DP) function internal_energy(current_denskern,update_ham)

      use constants, only: paw_en_size
      use fragment_scfmi, only: denskern_R, scfmi_construct_nonorth_kernel
      use hamiltonian, only: hamiltonian_dens_dep_matrices
      use rundat, only: pub_num_kpoints, PUB_1K, pub_eda_scfmi
      use utils, only: utils_assert

      implicit none

      ! Argument
      type(DKERN), intent(inout) :: current_denskern
      logical, intent(in) :: update_ham

      ! Local variables
      real(kind=DP) :: hubbard_energy
      real(kind=DP) :: paw_sphere_energies(paw_en_size)

      ! jme: KPOINTS_DANGER
      ! Some parts of this subroutine enforce a single k-point.
      call utils_assert(pub_num_kpoints == PUB_1K, &
           'Function internal_energy (penalty_mod) not ready yet for more&
           & than one k-point.')

      ! mjsp: for SCF MI calculation, calculate the internal energy using the
      ! orthogonalised representation
      if (pub_eda_scfmi) then

         ! mjsp: update denskern_R: proper representation of the density kernel
         ! in the full overlap
         call scfmi_construct_nonorth_kernel(current_denskern,rep)

         ! mjsp: calculate the density using the denskern_R kernel
         call hamiltonian_dens_dep_matrices(ham, lhxc_fine, internal_energy, &
              lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
              ngwf_basis, hub_proj_basis, hub, denskern_R%kern, mdl, hfxstate, &
              update_ham, lhxc_fixed=.false.)

      else

         ! ndmh: calculate density dependent energies and matrices
         call hamiltonian_dens_dep_matrices(ham, lhxc_fine, internal_energy, &
              lhxc_energy, hubbard_energy, paw_sphere_energies, rep, ngwf_basis, &
              hub_proj_basis, hub, current_denskern%kern, mdl, hfxstate, &
              update_ham, lhxc_fixed=.false.)

      end if

    end function internal_energy


    !******************************************************
    !******************************************************


    real(kind=DP) function internal_penalty_functional(current_denskern)

      use kernel, only: kernel_validate_ks
      use rundat, only: pub_spin_fac
      use sparse_embed, only: sparse_embed_trace, sparse_embed_product, &
          sparse_embed_scale, sparse_embed_array_trace
      implicit none

      ! Arguments
      type(DKERN), intent(inout) :: current_denskern

      ! Local variables
      integer :: is
      !integer :: ierr,i,j,num
      real(kind=DP) :: pen
      type(SPAM3_EMBED_ARRAY) :: ksksc, ksc

      ! Allocate local arrays
      call sparse_embed_array_create(ksc, current_denskern%ks)
      call sparse_embed_array_create(ksksc, current_denskern%ks, current_denskern%ks)

      pen = 0.0_DP

      call kernel_validate_ks(current_denskern, rep%overlap)

      ! Calculate penalty funtional
      ! ndmh: denskern%ks is still valid
      call sparse_embed_array_copy(ksc, current_denskern%ks)
      call sparse_embed_array_scale(ksc, -1.0_DP, 1.0_DP)              ! I - K.S
      call sparse_embed_array_product(ksksc, current_denskern%ks, ksc) ! K.S.(I - K.S)
      ! KPOINTS_DANGER: this sum is also over k-point indices!
      call sparse_embed_array_trace(trace_array, ksksc, ksksc)
      pen = sum(trace_array)

      ! Deallocate local arrays
      call sparse_embed_array_destroy(ksksc)
      call sparse_embed_array_destroy(ksc)

      internal_penalty_functional = 0.5_DP * pub_spin_fac * pub_pen_param * pen

    end function internal_penalty_functional



    !******************************************************
    !******************************************************


    real(kind=DP) function internal_gradient_norm()

      use kernel, only: kernel_commutator
      use sparse_embed, only: sparse_embed_trace, &
          sparse_embed_num_element, sparse_embed_array_trace
      implicit none

      ! Local variable
      integer :: is, ik
      real(kind=DP) :: num_els

      ! COVARIANT gradient (not tensor-corrected)
      call internal_co_gradient()

      ! CONTRAVARIANT (tensor-corrected) gradient
      call sparse_embed_array_copy(con_gradient, co_gradient)
      call internal_con_gradient(con_gradient)

      ! Project contravariant gradient to conserve electron number
      call internal_correct_ne_gradients(co_gradient, con_gradient, &
           rep%inv_overlap, rep%overlap)

      ! RMS value of gradient is a measure of the error in the current
      ! density kernel
      ! jme: KPOINTS_DANGER: these two sums are also over k-point indices!
      call sparse_embed_array_trace(trace_array, con_gradient, co_gradient)
      internal_gradient_norm = sum(trace_array)
      num_els = 0.0_DP
      do ik = 1, denskern%kern%num_kpoints
         do is = 1, denskern%kern%num_spins
            num_els = num_els + sparse_embed_num_element(denskern%kern%m(is,ik))
         end do
      end do
      internal_gradient_norm = sqrt(internal_gradient_norm / num_els)

      ! AAM: diagnostics
      ! rab207: moved kernel choice into kernel_mod.F90
      commutator = kernel_commutator(denskern, ham%ham, &
           rep%overlap, rep%inv_overlap)

    end function internal_gradient_norm


    !******************************************************
    !******************************************************


    subroutine internal_search_direction

      use sparse_embed, only: sparse_embed_trace, sparse_embed_scale, &
           sparse_embed_array_trace
      implicit none

      ! Local variable
      integer :: is

      if (iteration > 1 .and. rms_gradient <= 1.01_DP) then
         cg_coeff = internal_polak_cg_coeff(old_con_direction, &
              old_co_gradient,co_gradient,con_gradient)
      else
         cg_coeff = 0.0_DP
      end if

      if (abs(cg_coeff) > eps) then
         cg_count = cg_count + 1
      else
         cg_count = 0
      end if

      ! reset cg_coeff every 'cg_max' steps
      if (cg_count > cg_max) then
         cg_coeff = 0.0_DP
         cg_count = 0
      end if

      ! Determine the line search direction
      call sparse_embed_array_copy(con_direction, con_gradient)
      ! jme: KPOINTS_DANGER: is this hack affected by the number of k-points??
      ! pdh: hack to keep step lengths the same for spin polarised
      ! pdh: systems: multiply by a factor of two
      call sparse_embed_array_scale(con_direction, -real(pub_num_spins,kind=DP))
      call sparse_embed_array_axpy(con_direction, old_con_direction, cg_coeff)

      ! Find the derivative of the function along the search direction
      ! KPOINTS_DANGER: check this for more than 1 k-point
      call sparse_embed_array_trace(trace_array, co_gradient, con_direction)
      Qslope = sum(trace_array)

      ! Check search direction points downhill!!!

      if (Qslope > 0.0_DP) then
         ! First try steepest descent direction i.e. reset conjugate gradients
         if (pub_on_root .and. (pub_output_detail >= NORMAL)) &
              write(stdout,'(a)') '    !!! WARNING: positive line search &
              &gradient: resetting CG !!!'
         call sparse_embed_array_copy(con_direction, con_gradient)
         call sparse_embed_array_scale(con_direction, -1.0_DP)
         ! KPOINTS_DANGER: check this for more than 1 k-point
         call sparse_embed_array_trace(trace_array, co_gradient, con_direction)
         Qslope = sum(trace_array)

         ! If this is still positive then things are really pear-shaped i.e.
         ! overlap matrix isn't positive definite or something...
         ! Reverse search direction!
         if (Qslope > 0.0_DP) then
            ! KPOINTS_DANGER: check this for more than 1 k-point
            call sparse_embed_array_scale(con_direction,-1.0_DP)
            Qslope = -Qslope
         end if
      end if

      if (pub_output_detail >= VERBOSE .and. pub_on_root) &
           write(stdout,'(a,f20.12)') '            G0: ',Qslope

    end subroutine internal_search_direction


    !******************************************************
    !******************************************************


    subroutine internal_print

      use rundat, only: pub_spin_fac
      use sparse_embed, only: sparse_embed_trace, sparse_embed_array_trace
      implicit none

      ! Local variable
      character(55) :: fmt, tmp
      real(kind=DP) :: ne_new

      ! KPOINTS_DANGER: this sum is also over k-point indices!
      call sparse_embed_array_trace(trace_array, denskern%kern, rep%overlap)
      ne_new = sum(trace_array)
      ne_new = ne_new * pub_spin_fac

      if (pub_on_root) then

         ! ndmh: adapt format string so that it always displays sensible
         ! ndmh: numbers of digits for total energy and Ne
         write(tmp,'(a)')'(i2,'
         if(abs(total_energy)<100000.0_DP) then
            write(fmt,'(a,a)') trim(tmp),'f22.14'
         else
            write(fmt,'(a,a)') trim(tmp),'f22.12'
         end if
         write(tmp,'(a,a)') trim(fmt),',2f14.11,f7.4,f7.3,'
         if(ne_new<10000_DP) then
            write(fmt,'(a,a)') trim(tmp),'f14.8)'
         else if (ne_new<100000.0_DP) then
            write(fmt,'(a,a)') trim(tmp),'f14.7)'
         else
            write(fmt,'(a,a)') trim(tmp),'f14.5)'
         end if

         write(stdout,fmt) iteration, total_energy, rms_gradient, commutator, &
             line_search_coeff, cg_coeff, ne_new

      end if

    end subroutine internal_print


    !******************************************************
    !******************************************************


    real(kind=DP) function internal_quadratic_step()

      ! pdh: Three pieces of information: sufficient to fit a quadratic
      ! cks: Q=ax^2+bx+c
      ! cks: c=Qinitial
      linear_term =  Qslope * trial_step  ! cks: bx
      quadratic_term = Qinitial + linear_term - Qtrial  ! cks: ax^2

      if (pub_output_detail >= VERBOSE .and. pub_on_root) then
         write(stdout,'(a,f20.12)') 'Quadratic term: ',quadratic_term
         write(stdout,'(a,f20.12)') '   Linear term: ',linear_term
      end if

      if (abs(quadratic_term / Qinitial) > epsilon(1.0_DP)) then
         internal_quadratic_step = &
              0.5_DP * trial_step * linear_term / quadratic_term ! b/2a
         if (pub_output_detail >= VERBOSE .and. pub_on_root) &
              write(stdout,'(a,f20.12)') '  Optimal step: ',internal_quadratic_step
         Qpredict = Qinitial + &
              0.25_DP*linear_term**2/quadratic_term ! c+b^2/2a
      else
         internal_quadratic_step = -1.0_DP
         if (pub_output_detail >= VERBOSE .and. pub_on_root) &
              write(stdout,'(a)') '     !!! Quadratic term too small;' &
              // 'proceeding to cubic fit !!!'
      end if

    end function internal_quadratic_step


    !******************************************************
    !******************************************************


    real(kind=DP) function internal_cubic_step()

      ! Introduce shorthand working variables:
      !   xx & yy are the points at which the energy has been evaluated
      !   aa and bb are the quadratic and cubic polynomial coefficients
      ! cks: here's what PDH's code means: f= a*x^3+b*x^2+c*x+d
      ! cks: d=Qinitial=f1, c=Qslope=f1', b=aa, a=bb,
      !      f2=Qtrial, f3=Qoptimal

      internal_cubic_step = 0.15_DP
      xx = trial_step ; yy = optimal_step
      aa = ((yy*yy*yy*Qtrial - xx*xx*xx*Qoptimal)/(yy-xx) &
           - (xx*xx+xx*yy+yy*yy)*Qinitial - xx*yy*(xx+yy)*Qslope) &
           / (xx*xx*yy*yy)
      bb = ((yy*yy*Qtrial - xx*xx*Qoptimal)/(xx-yy) &
           + (xx+yy)*Qinitial + xx*yy*Qslope) / (xx*xx*yy*yy)
      disc=-1.0_DP
      if (abs(bb*yy/aa) > epsilon(1.0_DP)) then    ! avoid div by zero
         aon3b = aa / (3.0_DP * bb)
         disc = aon3b * aon3b - Qslope / (3.0_DP * bb)  ! discriminant
         if (disc >= 0.0_DP) then
            minimum_step = -aon3b + sign(sqrt(disc), bb)
            if (pub_output_detail >= VERBOSE .and. pub_on_root) &
                 write(stdout,'(a,f20.12)') 'Cubic LS coeff: ',minimum_step
            ! cks: energy according to cubic minimum
            Qminimum = Qinitial + minimum_step * &
                 (Qslope + minimum_step * (aa + minimum_step * bb))
            if (pub_output_detail >= VERBOSE .and. pub_on_root) &
                 write(stdout,'(a,f20.12)') ' Predicted min: ',Qminimum
            ! cks: set optimal_step to cubic minimum
            internal_cubic_step = minimum_step
            if (pub_output_detail >= VERBOSE .and. pub_on_root) &
                 write(stdout,'(a)') '    !!! Cubic minimum step &
                 &length accepted !!!'
         end if
      end if

      ! Choose safe length if unsuccessful

      if (disc < 0.0_DP ) then
         if (pub_on_root .and. (pub_output_detail >= NORMAL)) then
            write(stdout,'(a)') '    !!! WARNING: Line search was &
                 &unsuccessful: coefficient set manually !!!'
         end if
         !internal_cubic_step = 0.15_DP
      end if

    end function internal_cubic_step


    !******************************************************
    !******************************************************

    ! Written by Arash Mostofi, March 2003
    ! Adapted for penalty functional by Peter Haynes, July 2003
    ! Rearranged to use pre-calculated ks and thus save 1 sparse_product
    ! by Nick Hine, Nov 2007
    subroutine internal_co_gradient

      use kernel, only: kernel_validate_ks
      use rundat, only: pub_spin_fac, pub_num_kpoints, PUB_1K
      use sparse, only: sparse_product, sparse_scale
      implicit none

      ! Local variables
      integer :: is
      type(SPAM3_EMBED_ARRAY) :: sks, ksks

      ! jme: KPOINTS_DANGER
      ! Some parts of this subroutine enforce a single k-point.
      call utils_assert(pub_num_kpoints == PUB_1K, &
           'Subroutine internal_co_gradient (penalty_mod) not ready yet for more&
           & than one k-point.')

      ! Allocate workspace
      call sparse_embed_array_create(sks, rep%overlap, denskern%ks)
      call sparse_embed_array_create(ksks, denskern%ks, denskern%ks)

      call kernel_validate_ks(denskern, rep%overlap)


      ! Calculate gradient := H + (alpha/2) S.K.(I - 2 S.K).(I - S.K).S
      call sparse_embed_array_product(ksks, denskern%ks, denskern%ks) ! KSKS := KS.KS
      call sparse_embed_array_scale(ksks, 2.0_DP, 1.0_DP)         ! KSKS := 1+2KS.KS
      call sparse_embed_array_axpy(ksks, denskern%ks, -3.0_DP)     ! KSKS := 1-3KS+2KS.KS
      call sparse_embed_array_product(sks, rep%overlap, denskern%ks)   ! SKS := S.KS
      call sparse_embed_array_product(co_gradient, sks, ksks) ! G := SKS.(1-3KS+2KSKS)
      call sparse_embed_array_scale(co_gradient, pub_spin_fac * pub_pen_param)
      do is=1,co_gradient%num_spins
         call sparse_embed_axpy(co_gradient%m(is,PUB_1K), ham%ham(is), pub_spin_fac)
      end do

      ! Deallocate workspace
      call sparse_embed_array_destroy(ksks)
      call sparse_embed_array_destroy(sks)

    end subroutine internal_co_gradient


    !******************************************************
    !******************************************************


    subroutine internal_con_gradient(gradient)

      use sparse, only: sparse_product

      implicit none

      type(SPAM3_EMBED_ARRAY), intent(inout) :: gradient

      ! Local variables
      integer :: is
      type(SPAM3_EMBED_ARRAY) :: sinvg

      ! Allocate workspace
      call sparse_embed_array_create(sinvg, rep%inv_overlap, gradient)

      ! Calculate contravariant gradient
      call sparse_embed_array_product(sinvg,rep%inv_overlap, gradient)
      call sparse_embed_array_product(gradient, sinvg, rep%inv_overlap)

      ! Deallocate workspace
      call sparse_embed_array_destroy(sinvg)

    end subroutine internal_con_gradient


    !******************************************************
    !******************************************************


    subroutine internal_correct_ne_gradients(co_gradient, &
         con_gradient, inv_overlap, overlap)

      use rundat, only: pub_num_kpoints, PUB_1K
      use sparse_embed, only: sparse_embed_trace, sparse_embed_array_trace, &
           sparse_embed_array_axpy

      implicit none

      ! Arguments
      type(SPAM3_EMBED_ARRAY), intent(inout) :: co_gradient
      type(SPAM3_EMBED_ARRAY), intent(inout) :: con_gradient
      type(SPAM3_EMBED), intent(in)    :: inv_overlap
      type(SPAM3_EMBED), intent(in)    :: overlap

      ! Local variables
      integer :: is
      real(kind=DP) :: step
      real(kind=DP) :: rate_of_change_ne ! Rate of change of electron number
      real(kind=DP) :: trunc_rank

      ! jme: KPOINTS_DANGER
      ! Some parts of this subroutine are untested for more than 1 k-point.
      call utils_assert(pub_num_kpoints == PUB_1K, &
           'Subroutine internal_correct_ne_gradients (penalty_mod) not checked yet for more&
           & than one k-point.')

      ! Calculate amount of S^-1 to add to correct contravariant gradient
      call sparse_embed_trace(trunc_rank,overlap,inv_overlap)
      trunc_rank = trunc_rank * pub_num_spins
      ! KPOINTS_DANGER: this sum is also over k-point indices!
      call sparse_embed_array_trace(trace_array, con_gradient, overlap)
      rate_of_change_ne = sum(trace_array)
      step = -rate_of_change_ne / trunc_rank

      ! Correct contravariant gradient
      call sparse_embed_array_axpy(con_gradient, inv_overlap, step)

      ! Correct covariant gradient
      call sparse_embed_array_axpy(co_gradient, overlap, step)

    end subroutine internal_correct_ne_gradients



    !******************************************************
    !******************************************************

    ! cks: written by cks on 29/6/2001
    real(kind=DP) function internal_polak_cg_coeff(old_con_direction, &
         old_co_gradient,co_gradient,con_gradient)

      use rundat, only: PUB_1K
      use sparse_embed, only: sparse_embed_array_trace

      implicit none

      ! cks: this subroutine returns the coefficient according to Polak for the
      ! cks: combination of the previous search direction with the negative of
      ! cks: the gradient in the conjugate gradient method.
      ! cks: b_{r+1}=g_{r+1}*(g_{r+1}-g_r)/( p_r* (g_{r+1}-g_r) )
      ! cks: It is written so that covariant vectors are contracted with
      ! cks: contravariant vectors.

      type(SPAM3_EMBED_ARRAY), intent(in) :: old_con_direction
      type(SPAM3_EMBED_ARRAY), intent(in) :: old_co_gradient
      type(SPAM3_EMBED_ARRAY), intent(in) :: co_gradient
      type(SPAM3_EMBED_ARRAY), intent(in) :: con_gradient

      ! cks: internal declarations
      integer :: is
      real(kind=DP) :: denominator

      ! jme: KPOINTS_DANGER
      ! Some parts of this subroutine are untested for more than 1 k-point.
      call utils_assert(pub_num_kpoints == PUB_1K, &
           'Function internal_polak_cg_coeff (penalty_mod) not checked yet for more&
           & than one k-point.')

      ! These sums are also over k-point indices.
      call sparse_embed_array_trace(trace_array, old_con_direction, co_gradient)
      call sparse_embed_array_trace(trace_array2, old_con_direction, old_co_gradient)
      denominator = sum(trace_array) - sum(trace_array2)
      ! pdh: hack to maintain consistency with spin polarisation
      denominator = denominator / pub_num_spins

      if (abs(denominator) > 1.0e-6_DP) then ! avoid div by zero

         call sparse_embed_array_trace(trace_array, con_gradient, co_gradient)
         call sparse_embed_array_trace(trace_array2, con_gradient, old_co_gradient)
         internal_polak_cg_coeff = sum(trace_array) - sum(trace_array2)
         internal_polak_cg_coeff = internal_polak_cg_coeff / denominator

      else

         if (pub_on_root .and. pub_output_detail >= NORMAL) then
            write(stdout,'(a)') &
                 'WARNING: zero denominator in internal_polak_cg_coef'
            write(stdout,'(a)') &
                 'WARNING: setting internal_polak_cg_coeff to zero'
         end if
         internal_polak_cg_coeff = 0.0_DP

      end if

      if (abs(internal_polak_cg_coeff) > 5.0_DP) then
         if (pub_on_root .and. pub_output_detail >= NORMAL) then
            write(stdout,'(a,f8.5)') 'WARNING: internal_polak_cg_coeff = ', &
                 internal_polak_cg_coeff
            write(stdout,'(a)') &
                 'WARNING: setting internal_polak_cg_coeff to zero'
         end if
         internal_polak_cg_coeff = 0.0_DP
      end if

    end function internal_polak_cg_coeff

  end subroutine penalty_denskernel_optimise_cg

end module penalty
