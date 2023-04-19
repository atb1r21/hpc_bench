! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!            Electronic energy initialisation module             !
!                                                                !
! This module contains routines used in the initialisation of the!
! electronic energy with respect to the NGWFs and the density    !
! kernel.                                                        !
!----------------------------------------------------------------!
! Originally written by Chris-Kriton Skylaris in 2000.           !
! Subsequent modifications by Chris-Kriton Skylaris,             !
! Arash A. Mostofi, Peter D. Haynes and Nicholas D.M. Hine.      !
! Modified by Nicholas Hine to use the NGWF_REP container for    !
! NGWF-related data and matrices, October 2010.                  !
! Converted into a separate module by Chris-Kriton Skylaris,     !
! 23/07/2014.                                                    !
!================================================================!


module electronic_init

  use constants, only: DP
  use rundat, only: pub_debug_on_root

  implicit none

  private

  public :: electronic_init_denskern

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine electronic_init_denskern(rep, ham, denskern, &
       lhxc_fine, ngwf_basis, proj_basis, nl_projectors, &
       hub_proj_basis, hub, mdl, hfxstate, val_rep, val_ngwf_basis, &
       val_dkn, val_ham, lhxc_fixed_in, edft_for_init, kpt, &
       dfdtau_fine, restart_rootname, force_no_IO, force_init_pot, &
       allow_pseuinvS)

    !=================================================================!
    ! This subroutine initialises the denskern DKERN, using one of    !
    ! several techniques: diagonalisation of the initial Hamiltonian, !
    ! Palser-Manolopoulos Canonical purification, or minimisation of  !
    ! the penalty functional. It also calls electronic_init_pot which !
    ! initialises the whole-cell arrays.                              !
    !-----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 21/2/2004.                  !
    ! Modified by Peter D. Haynes on 5/7/2004 for Fourier             !
    ! parallelisation.                                                !
    ! Spin polarised by Peter Haynes, July 2006                       !
    ! Modified for NLCC by Nicholas Hine, January 2009                !
    ! DFT+U added by D. D. O'Regan, April 2009                        !
    ! Adapted for SPAM3 and function basis by Nicholas Hine in        !
    ! May-July 2009.                                                  !
    ! Modified by Jacek Dziedzic on 13/05/2010 to include correction  !
    ! due to smeared ions and implicit solvent.                       !
    ! Reorganised, made to use NGWF_REP and NGWF_HAM, split into two  !
    ! parts and generally tidied-up by Nicholas Hine in October 2010. !
    ! Modified by Laura Ratcliff in October 2010 to allow for         !
    ! use in conduction calculations.                                 !
    ! Modified by Simon Dubois in May 2010 to include various         !
    ! mixing schemes                                                  !
    ! Modified by Chris-Kriton Skylaris to fully optimise the density !
    ! kernel when starting a new calculation, 22/07/2014.             !
    ! Modified for embedding structures by Robert Charlton, Aug 2017. !
    ! Modified further for embedding by Joseph Prentice, Sep 2020     !
    !=================================================================!

    use augmentation, only: augmentation_screen_dij, aug_projector_denskern, &
         aug_nonlocal_mat
    use comms, only: pub_on_root
    use constants, only: DP, VERBOSE, stdout, DN, UP, paw_en_size, &
         EDA_ISOLATED, EDA_CTFRAGLOC_DEVEL, EDA_INIT
    use datatypes, only: data_functions_dot, COEF
    use dense, only: DEM
    use dftb, only: dftb_calculate_ham
    use ensemble_dft, only: edft_calculate
    use ensemble_dft_type, only: EDFT_MODEL
    use ensemble_dftb, only: edftb_calculate
    use electronic_history, only: elec_history_compose_dkn
    use fragment_matrix_ops, only: fmo_load_internal_blkdiag_denskern, &
         fmo_freeze_mat_spam_nxn
    use function_basis, only: FUNC_BASIS
    ! agrecokpt: needed in hamiltonian_apply_phases_tb_method
    use geometry, only: POINT
    use hamiltonian, only: hamiltonian_dens_indep_matrices, &
         hamiltonian_proj_cond_matrix, hamiltonian_dens_dep_matrices, &
         hamiltonian_build_matrix, hamiltonian_apply_phases_tb_method, &
         hamiltonian_proj_embed_matrix
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL, hubbard_projector_update, &
         hubbard_ham_matrix
    use integrals, only: integrals_locpot, integrals_div_locpot_grad_functions
    use kernel, only: DKERN, kernel_init_core_ham, kernel_workspace_invalidate
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use palser_mano, only: palser_mano_kernel_optimise
    use parallel_strategy, only: PARAL_INFO
    use paw, only: paw_projector_denskern_init
    use penalty, only: penalty_denskernel_optimise_cg
    use projectors, only: PROJECTOR_SET
    use pseudopotentials, only: PSEUDO_SPECIES
    use restart, only: restart_kernel_read, restart_kernel_write
    use rundat, only: pub_output_detail, pub_maxit_pen, pub_dftb, &
         lnv_threshold_orig, pub_write_denskern, pub_maxit_palser_mano, &
         pub_read_denskern, cond_read_denskern, pub_any_nl_proj, &
         pub_coreham_denskern_guess, pub_hubbard_restart, pub_aug, pub_paw, &
         pub_write_converged_dk_ngwfs, pub_cond_calculate, &
         pub_devel_code, pub_hubbard_atomsolve, pub_hubbard, &
         pub_hub_on_the_fly, pub_task, pub_edft, pub_num_spins, pub_spin_fac, &
         pub_num_kpoints, PUB_1K, pub_kpoint_method, pub_edft_init_maxit, &
         pub_eda_scfmi, pub_eda, pub_eda_mode, pub_eda_read_frags, pub_eda_read_super, &
         ! agrecokpt: kpoint dependence included in definition of projs_on_grid
         pub_eda_split_atoms, pub_frag_counter, pub_frag_counter2, &
         pub_xc_ke_density_required, pub_check_hermitian, pub_imag_thr, &
         pub_active_region, pub_emft, pub_debug, pub_inner_loop_iteration
    use services, only: services_flush
    use sparse, only: SPAM3, sparse_scale, sparse_copy, &
         sparse_axpy, sparse_create, sparse_destroy, sparse_trace, &
         sparse_check_hermitian
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy, &
         sparse_embed_axpy, sparse_embed_copy, sparse_embed_scale, &
         sparse_embed_create, sparse_embed_mat2array, sparse_embed_trace, &
         SPAM3_EMBED_ARRAY, sparse_embed_destroy, sparse_embed_check_hermitian, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, &
         utils_assert
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, &
         utils_assert

    implicit none

    ! Arguments (inputs)
    type(EDFT_MODEL), intent(inout), optional  :: edft_for_init
    type(MODEL), intent(inout) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    type(FUNC_BASIS), intent(inout) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(:)
    type(FUNC_BASIS), intent(inout) :: hub_proj_basis(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(in), optional :: kpt(3)
    logical, intent(in), optional :: force_no_IO
    logical, intent(in), optional :: force_init_pot
    logical, intent(in), optional :: allow_pseuinvS ! ja531-> allow inverse S to be a pseudoinverse

    ! Arguments (outputs)
    type(DKERN), intent(inout) :: denskern
    type(NGWF_REP), intent(inout) :: rep
    type(NGWF_HAM), intent(inout) :: ham
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, rep%nsub)
    ! JCW: Optional gradient of XC energy per unit volume wrt kinetic energy
    ! JCW: density. Required to evaluate lhxc matrix elements for meta-GGAs.
    real(kind=DP), optional, intent(inout) :: dfdtau_fine(:,:,:,:)
    ! JCW: In general, we have
    ! JCW:    dfdtau_fine(mdl%fine_grid%ld1,&
    ! JCW:    mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    ! JCW: or
    ! JCW:    dfdtau_fine(1,1,1,1)
    ! JCW: when it is unneeded.

    ! lr408: Optional arguments needed for conduction calculation
    type(NGWF_REP), optional, intent(in) :: val_rep
    type(NGWF_HAM), optional, intent(in) :: val_ham
    type(FUNC_BASIS), optional, intent(in) :: val_ngwf_basis(:)
    type(SPAM3_EMBED_ARRAY), optional, intent(in) :: val_dkn
    logical, optional, intent(in) :: lhxc_fixed_in
    character(len=*),optional,intent(in) :: restart_rootname

    ! Local Variables

    logical :: loc_allow_pseuinvS
    logical :: loc_force_no_IO
    type(SPAM3_EMBED), allocatable :: rho_ij(:)
    type(SPAM3), allocatable :: rho_ij_tmp(:)
    type(SPAM3) :: ham_array(pub_num_spins), kern_array(pub_num_spins)

    real(kind=DP) :: filling
    real(kind=DP) :: Ecurr
    real(kind=DP) :: num_elec
    integer :: is, ik
    integer :: ierr
    logical :: init_pot, read_dkn, cond_call, frag_call, super_call, lhxc_fixed
    logical :: initialise_new
    logical :: emft_call   ! jcap
    type(DEM), allocatable :: eigs_dens(:)
    real(kind=DP), allocatable :: eigen_en(:,:)
    real(kind=DP), allocatable :: eigen_occ(:,:)
    real(kind=DP) :: lhxc_energy, hubbard_energy
    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    real(kind=DP) :: total_energy ! cks: dummy for kernel optimisation
    ! agrecokpt
    real(kind=DP) :: loc_kpt(3)
    type(POINT) :: loc_kpt_point
    ! agrecocmplx
    logical :: loc_cmplx
    real(kind=DP) :: rms_hermitian_check
    logical :: loc_force_init_pot
    type(PARAL_INFO), pointer :: par
    integer :: isub, jsub, mrows, ncols, lhxc_index
    real(kind=DP) :: debug_trace, debug_val
    type(COEF)    :: debug_dot

    ! ------------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &electronic_init_denskern'

    pub_inner_loop_iteration = 0

    loc_force_no_IO=.false.
    if(present(force_no_IO)) then
       loc_force_no_IO=force_no_IO
    end if

    if (.not.pub_xc_ke_density_required) then
       ! JCW: Non-tau-dependent XC functional
       ! JCW: dfdtau_fine is either not present, or has a length of 1
       if ( present(dfdtau_fine) ) then
          call utils_assert( &
               size(dfdtau_fine) == 1,&
               "Error in electronic_init_denskern: &
               &pub_xc_ke_density_required is false, but dfdtau_fine is present &
               &and has a size > 1.")
       end if
    else
       ! JCW: tau-dependent XC functional
       ! JCW: Check dfdtau_fine is present
       call utils_assert( present(dfdtau_fine),&
            "Error in electronic_init_denskern: &
            &pub_xc_ke_density_required is true, but dfdtau_fine is not &
            &present.")
    end if

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine electronic_init_denskern not ready yet for more&
         & than one k-point.')

    ! agrecokpt: currently only KP method is being implemented
    ! agrecokpt: now also TB method implemented for the basic cases
    ! (i.e. no Hubbard, ecc...). Need to check everything is correct!
    !call utils_assert(pub_kpoint_method == 'KP', &
    !     'Subroutine electronic_init_denskern currently supports&
    !     & only KP method for BZ sampling')

    ! agrecokpt: default is Gamma
    if (present(kpt)) then
       loc_kpt = kpt
       loc_kpt_point%x = kpt(1)
       loc_kpt_point%y = kpt(2)
       loc_kpt_point%z = kpt(3)
    else
       loc_kpt(:) = 0.0_DP
       loc_kpt_point%x = 0.0_DP
       loc_kpt_point%y = 0.0_DP
       loc_kpt_point%z = 0.0_DP
    end if

    ! rc2013: extract size of embedding system
    mrows = denskern%kern%mrows
    ncols = denskern%kern%ncols

    ! agrecocmplx
    loc_cmplx = rep%overlap%iscmplx

    ! lr408: If optional argument is present, we are initialising a conduction
    ! lr480: kernel, and so the potential does not need initialising.
    ! lr408: Also, read_dkn is set to match read_cond_denskern, rather
    ! lr408: than read_denskern
    if (present(val_rep)) then
       init_pot = .false.
       read_dkn = cond_read_denskern
       cond_call = .true.
       frag_call = .false.
       super_call = .false.
       lhxc_fixed = .true.
    ! mjsp: If this is a frozen density calculation, then
    ! mjsp: then kernel is constructed either:
    ! mjsp: from fragment density kernels stored on disk, or
    ! mjsp: from internally stored fragment density kernels.
    else if (pub_eda) then
       ! Isolated fragments:
       if (pub_eda_mode == EDA_ISOLATED) then
          read_dkn = pub_eda_read_frags
          frag_call = pub_eda_read_frags
          super_call = .false.
       ! Reinitialisation without read-in:
       else if (pub_eda_mode == EDA_INIT) then
          read_dkn = .false.
       ! Supermolecule:
       else
          read_dkn = pub_eda_read_super
          super_call = pub_eda_read_super
          frag_call = .false.
       end if
       init_pot = .true.
       cond_call = .false.
       lhxc_fixed = .false.
    else
       init_pot = .true.
       read_dkn = pub_read_denskern
       cond_call = .false.
       frag_call = .false.
       super_call = .false.
       lhxc_fixed = .false.
    end if
    if (present(lhxc_fixed_in)) then
       lhxc_fixed = lhxc_fixed_in
       if (lhxc_fixed) init_pot = .false.
    end if

    if (pub_dftb) init_pot = .false.

    loc_allow_pseuinvS=.false.
    if(present(allow_pseuinvS)) then
       loc_allow_pseuinvS=allow_pseuinvS
    end if

    ! ndmh: prevent initialisation of more cond states than the number of
    ! ndmh: (total states) - (valence states)
    ! tjz07: Commented out as it prevents generating Kernels in the case
    ! tjz07: that the number of COND NGWFs is smaller than the number of
    ! tjz07: val ngwfs
    !if (cond_call) then
    !   do is=1,pub_num_spins
    !      if (rep%n_occ(is)>(ngwf_basis%num - val_rep%n_occ(is))) then
    !         call utils_abort('Error in electronic_init_denskern: number &
    !                 &of ''occupied'' orbitals requested in conduction density &
    !                 &kernel = n, but only m states are available above Fermi &
    !                 &level in conduction NGWF set. Values for n and m were: ',&
    !                 rep%n_occ(is), ngwf_basis%num - val_rep%n_occ(is))
    !      end if
    !   end do
    !end if

    ! ndmh: Initialise parts of Hamiltonian that do not depend on the density
    ! lr408: If this is a conduction call, we need to pass on the valence rep
    ! lr408: and basis
    ! agrecokpt: added optional argument for single kpoint
    ! in KP method, uses k-point dependent projectors and adds extra terms
    ! to kinetic matrix; in TB method, includes k-point terms in overlap and
    ! hamiltonian matrices, projectors are kept kpt NON-dependent
    if (pub_dftb) then
       call dftb_calculate_ham(mdl, ham, rep, ngwf_basis)
    else
       call hamiltonian_dens_indep_matrices(rep, ngwf_basis, proj_basis, &
            nl_projectors, hub_proj_basis, hub, mdl, val_rep, val_ngwf_basis, &
            kpt=loc_kpt,allow_pseuinvS=loc_allow_pseuinvS)
    end if

    ! ddor: Build DFT+U projectors in the case where a restart
    !       file is used with a task such as PROPERTIES
    if ((pub_hubbard_restart .or. pub_hubbard_atomsolve .or. &
         & pub_hub_on_the_fly) .and. (pub_task .ne. 'HUBBARDSCF') .and. &
         & (.not. cond_call)) then

       if (any(loc_kpt/=0.0_DP)) then
          if (pub_debug_on_root) write(stdout,'(a)') &
             'WARNING: computation of hubbard projectors in &
             &electronic_init_denskern not yet checked for &
             &non-Gamma k-point'
       end if

       ! agrecokpt: not sure if kpt dependence is needed here...
       call hubbard_projector_update( &
            rep%ngwfs_on_grid(1), ngwf_basis(1), &
            nl_projectors(1), proj_basis(1), hub_proj_basis(1), hub, &
            rep%inv_overlap%p, rep%hub_overlap%p, &
            rep%hub_overlap_t%p, &
            rep%sp_overlap%p, rep%hub_proj_paw_overlap%p, &
            rep%hub_ngwf_paw_overlap%p, mdl)
    endif

    ! ndmh: Initialise the local pseudopotential, core density, and
    ! ndmh: if required, create the lhxc matrix in SPAM3 format
    if (init_pot) then
       loc_force_init_pot=.false.
       if(present(force_init_pot)) then
          loc_force_init_pot=force_init_pot
       end if
       call electronic_init_pot(lhxc_fine, ham%dijhat, mdl, loc_force_init_pot,&
            dfdtau_fine = dfdtau_fine)
    end if

    ! load density kernel from file if required
    initialise_new = loc_force_no_IO
    if(.not.loc_force_no_IO) then
       if (read_dkn .or. present(restart_rootname)) then

          ! ------- READ A DENSITY KERNEL FROM FILE -----------------------
          ! lr408: Pass extra flag to determine which kernel file to read from
          call restart_kernel_read(denskern%kern, read_cond=cond_call, &
               read_frag=frag_call,read_super=super_call,read_emft=pub_emft, &
               restart_rootname=restart_rootname)
          ! ------- END READ A DENSITY KERNEL FROM FILE--------------------

          ! lr408: Adding safeguard to ensure that the kernel contains the
          ! lr408: correct number of of 'occupied' states
          if (cond_call.and.(pub_task/='PROPERTIES_COND')) then
             ! Calculate current electron number
             do is=1,pub_num_spins
                ! KPOINTS_DANGER: only 1 k-point considered!
                ! agrecokpt: implementation prototype
                !num_elec = 0.0_DP
                !do ik=1,pub_num_kpoints
                !   this_kpt = mdl%unique_kpoints(ik)
                !   num_elec = num_elec + sparse_trace(denskern%kern%m(is,ik), &
                !              rep%overlap%m(1,ik)) * this_kpt%weight
                !end do
                call sparse_embed_trace(num_elec, denskern%kern%m(is,PUB_1K), &
                     rep%overlap)
                ! agrecokpt: here we need to do weight sum over k-points before
                ! dividing for rep%n_occ for current spin; here we assume that
                ! rep%n_occ has not yet been rescaled by k-point weight; if that's
                ! not the case, simply do sum(rep%n_occ(is,:))
                !reduced_occ = 0.0_DP
                !do ik=1,pub_num_kpoints
                !   reduced_occ = reduced_occ + rep%n_occ(is,ik)*mdl%unique_kpoints(ik)%weight
                !end do
                if (abs(num_elec/real(rep%n_occ(is,PUB_1K),DP)-1.0_DP)>0.05_DP) then
                   call utils_abort('Error in electronic_init_denskern: &
                           &number of ''occupied'' orbitals in input kernel is &
                           &''x'', which does not match the expected value of &
                           &''y'' for spin ''s''. Verify that cond_num_extra_&
                           &states and cond_num_extra_its are set correctly. &
                           &Values for y, s and x, follow.', rep%n_occ(is,PUB_1K), is, &
                           opt_real_to_print1=num_elec)
                end if
             end do
          end if

          ! ndmh: build Hamiltonian with this kernel, either because we are
          ! ndmh: re-diagonalising it, analysing it
          if ((index(pub_devel_code,'DIAGONALISE_ON_RESTART')>0).or. &
               (index(pub_devel_code,'INITIAL_EIGENSTATES')>0)) then

             ! KPOINTS_DANGER: only 1-kpoint considered!
             ! agrecokpt: included kpt argument to add extra terms
             ! when using the TB method
             call hamiltonian_dens_dep_matrices(ham, lhxc_fine, &
                  total_energy, lhxc_energy, hubbard_energy, &
                  paw_sphere_energies, rep, ngwf_basis, &
                  hub_proj_basis, hub, denskern%kern, mdl, hfxstate, &
                  .true., lhxc_fixed=.false., &
                  dfdtau_fine = dfdtau_fine, kpt=loc_kpt)

             ! agrecokpt: currently at this stage, all k-point
             ! dependence is included in each ham term separately
             ! as required by the method specified
             call hamiltonian_build_matrix(ham, rep)

          end if

          if (index(pub_devel_code,'DIAGONALISE_ON_RESTART')>0) then
             ! KPOINTS_DANGER: this must be checked for more than one k-point
             ! agrecokpt: hamiltonian and overlap should be already k-point dependent
             ! in TB method at this point? I suppose we simply use k-point
             ! dependent overlap and hamiltonian in TB to generate corresponding
             ! KP-dependent kernel
             do ik = 1, denskern%kern%num_kpoints
                do is = 1, pub_num_spins
                   call kernel_init_core_ham(denskern%kern%m(is,ik), &
                        ham%ham(is), rep%overlap, rep%n_occ(is,ik))


                   ! agrecocmplx: check denskern is hermitian if requested
                   if (pub_check_hermitian) then
                      rms_hermitian_check = &
                         sparse_embed_check_hermitian(denskern%kern%m(is,ik))

                      if (pub_on_root) write(stdout,'(a,i3,a,f8.5)') 'RMS denskern spin',&
                         is,' hermitian check: ', rms_hermitian_check

                      call utils_assert(rms_hermitian_check < pub_imag_thr, &
                           'WARNING: hermitian check for denskern matrix above threshold')

                   end if
                end do
             end do

          end if

       else if ((pub_eda) .and. .not.((pub_eda_mode .eq. EDA_ISOLATED) .or. &
                (pub_eda_mode .eq. EDA_CTFRAGLOC_DEVEL).or. &
                (pub_eda_mode .eq. EDA_INIT)) .and. .not.(pub_eda_split_atoms)) then

          ! mjsp: If this is a frozen density calculation and the density is
          ! stored in memory, then initialise
          ! to the frozen density kernel constructed from the monomer calculation:
          ! mjsp: Copy the elements of the internally stored fragment density
          ! kernels into the appropriate elements of the supermolecule density kernel
          call fmo_load_internal_blkdiag_denskern(denskern%kern, ngwf_basis)
          call kernel_workspace_invalidate(denskern)

       else
          ! agrecokpt: use kpoint dependence in case projection method is used
          ! this is mainly for testing, in the future simply take kpoints from
          ! denskern array
          call elec_history_compose_dkn(denskern, rep, ngwf_basis, &
               nl_projectors, proj_basis, mdl, ierr, kpt=loc_kpt)
          if (ierr .ne. 0) initialise_new = .true.

       endif


    else if ((pub_eda) .and. .not.((pub_eda_mode .eq. EDA_ISOLATED) .or. &
         (pub_eda_mode .eq. EDA_CTFRAGLOC_DEVEL).or. &
         (pub_eda_mode .eq. EDA_INIT)) .and. .not.(pub_eda_split_atoms)) then

       ! mjsp: If this is a frozen density calculation and the density is
       ! stored in memory, then initialise
       ! to the frozen density kernel constructed from the monomer calculation:
       ! mjsp: Copy the elements of the internally stored fragment density
       ! kernels into the appropriate elements of the supermolecule density kernel
       call fmo_load_internal_blkdiag_denskern(denskern%kern, ngwf_basis)
       call kernel_workspace_invalidate(denskern)

    else
       ! agrecokpt: use kpoint dependence in case projection method is used
       ! this is mainly for testing, in the future simply take kpoints from
       ! denskern array
       call elec_history_compose_dkn(denskern, rep, ngwf_basis, &
            nl_projectors, proj_basis, mdl, ierr, kpt=loc_kpt)
       if (ierr .ne. 0) initialise_new = .true.

    endif

    if (pub_coreham_denskern_guess .and. initialise_new) then

       ! ------- GUESS A DENSITY KERNEL --------------------------------

       if (.not. pub_dftb) then
          ! ndmh: Initialise the PAW nonlocal matrix, if required
          if (pub_aug) then
             ! JCW: Exit with error if KE density required, since KE density +
             ! augmentation is untested (temporary)
             if(pub_xc_ke_density_required) then
              call utils_abort('Error in electronic_init_denskern: &
                   &KE density evaluation and augmentation not tested.')
             end if

             ! Allocate storage for projector density kernel
             ! agrecokpt: should this depend on k-point as well?
             allocate(rho_ij(pub_num_spins),stat=ierr)
             call utils_alloc_check('electronic_init_denskern','rho_ij',ierr)
             do is=1,pub_num_spins
                rho_ij(is)%structure='E'
                call sparse_embed_create(rho_ij(is))
             end do

             ! For conduction calculation, copy in nonlocal energies from val ham
             if (cond_call) then
                do is=1,pub_num_spins
                   call sparse_embed_copy(ham%dijhat(is),val_ham%dijhat(is))
                end do
                ! jcap: need temporary arrays to pass in
                allocate(rho_ij_tmp(pub_num_spins),stat=ierr)
                call utils_alloc_check('electronic_init_denskern',&
                     'rho_ij_tmp',ierr)
                call sparse_embed_extract_from_array(rho_ij_tmp,rho_ij,1,1)
                call sparse_embed_extract_from_array(kern_array,val_dkn%m(:,PUB_1K))

                call aug_projector_denskern(rho_ij_tmp, kern_array, &
                     val_rep%sp_overlap%p)

                call sparse_embed_destroy_extracted_array(rho_ij_tmp,rho_ij,&
                     .true.,1,1)
                call sparse_embed_destroy_extracted_array(kern_array)
                deallocate(rho_ij_tmp,stat=ierr)
                call utils_dealloc_check('electronic_init_denskern',&
                     'rho_ij_tmp',ierr)

                do is=1,pub_num_spins
                   call sparse_embed_scale(rho_ij(is),pub_spin_fac)
                end do
             end if

             ! If this is a brand-new calculation, guess a kernel
             if (init_pot) then
                if (pub_paw) then
                   allocate(rho_ij_tmp(pub_num_spins),stat=ierr)
                   call utils_alloc_check('electronic_init_denskern',&
                        'rho_ij_tmp',ierr)
                   do isub=1,mdl%nsub
                      call sparse_embed_extract_from_array(rho_ij_tmp,rho_ij,&
                           isub,isub)
                      call paw_projector_denskern_init(rho_ij_tmp,&
                           mdl%regions(isub)%paw_sp,&
                           mdl%regions(isub)%radial_densities)
                      call sparse_embed_destroy_extracted_array(rho_ij_tmp,&
                           rho_ij,.true.,isub,isub)
                   end do
                   deallocate(rho_ij_tmp,stat=ierr)
                   call utils_dealloc_check('electronic_init_denskern',&
                        'rho_ij_tmp',ierr)
                end if
             end if

             ! Now calculate nonlocal matrix
             ! jcap: need temporary arrays to pass in
             allocate(rho_ij_tmp(pub_num_spins),stat=ierr)
             call utils_alloc_check('electronic_init_denskern',&
                  'rho_ij_tmp',ierr)
             call sparse_embed_extract_from_array(rho_ij_tmp,rho_ij,1,1)
             ! jcap: note that here kern_array is used because it is a
             ! spare sparse array, NOT because ham%dijhat is a kernel
             call sparse_embed_extract_from_array(kern_array,ham%dijhat,1,1)
             call sparse_embed_extract_from_array(ham_array,ham%nonlocpot,1,1)

             call aug_nonlocal_mat(ham_array,kern_array,rho_ij_tmp, &
                  rep%sp_overlap%p,mdl%regions(1)%pseudo_sp, &
                  mdl%regions(1)%paw_sp,show_matrices=.true.)

             call sparse_embed_destroy_extracted_array(ham_array,ham%nonlocpot,&
                  .true.,1,1)
             call sparse_embed_destroy_extracted_array(kern_array)
             call sparse_embed_destroy_extracted_array(rho_ij_tmp)
             deallocate(rho_ij_tmp,stat=ierr)
             call utils_dealloc_check('electronic_init_denskern',&
                  'rho_ij_tmp',ierr)

             ! agrecokpt: add TB method extra terms if requested
             ! check this is correct in PAW case
             do is=1,pub_num_spins
                if ((pub_kpoint_method == 'TB') .and. (loc_kpt_point%x/=0.0_DP.or.&
                     loc_kpt_point%y/=0.0_DP.or.loc_kpt_point%z/=0.0_DP)) then
                   if (pub_debug_on_root) write(stdout,'(a,i3)') 'DEBUG: &
                      &Adding TB terms in nonlocpot matrix spin ', is

                   ! rc2013: abort if using multiple subsystems
                   call utils_assert(mdl%nsub == 1, &
                        'Error in electronic_init_denskern: non-Gamma &
                        &k-points not allowed with more than 1 subsystem.')
                   call hamiltonian_apply_phases_tb_method(ham%nonlocpot(is)%p, &
                        ngwf_basis(1), mdl, loc_kpt_point)

                end if

                ! agrecocmplx: check nonlocpot matrix is hermitian if requested
                if (pub_check_hermitian) then
                   rms_hermitian_check = sparse_embed_check_hermitian(ham%nonlocpot(is))

                   if (pub_on_root) write(stdout,'(a,i3,a,f8.5)') 'RMS nonlocpot spin ',&
                      is,' hermitian check: ', rms_hermitian_check

                   call utils_assert(rms_hermitian_check < pub_imag_thr, &
                        'WARNING: hermitian check for nonlocpot matrix above threshold')

                end if
             end do

             ! Deallocate storage for projector density kernel
             do is=pub_num_spins,1,-1
                rho_ij(is)%structure='E'
                call sparse_embed_destroy(rho_ij(is))
             end do
             deallocate(rho_ij,stat=ierr)
             call utils_dealloc_check('electronic_init_denskern','rho_ij',ierr)

          end if

          ! ndmh: Calculate and include spin splitting part of Hubbard
          ! ndmh: Hamiltonian if required
          if (pub_hubbard) then
             if (any(loc_kpt/=0.0_DP)) then
                if (pub_debug_on_root) write(stdout,'(a)') &
                   'WARNING: computation of hubbard matrix in &
                   &electronic_init_denskern not yet checked for &
                   &non-Gamma k-point'
             end if

             if (pub_num_spins==2) then
                call sparse_copy(hub%projector_ham(1),hub%up_matrix)
                call sparse_copy(hub%projector_ham(2),hub%down_matrix)
             end if

             ! jcap: make temporary hamiltonian array
             call sparse_embed_extract_from_array(ham_array,ham%hubbard_ham)
             call hubbard_ham_matrix(hub,ham_array, &
                  rep%hub_overlap%p,rep%hub_overlap_t%p)
             call sparse_embed_destroy_extracted_array(ham_array,ham%hubbard_ham,.true.)
          end if

          ! ndmh: construct the rest of the initial guess hamiltonian
          do is=1,pub_num_spins

             ! ndmh: calculate lhxc matrix
             ! rc2013: build each matrix block separately
             do jsub=1,ncols
                do isub=1,mrows
                   ! rc2013: the LHXC potential required will depend on the region we're in
                   lhxc_index = isub
                   if((isub .ne. jsub) .and. (isub == pub_active_region)) &
                        lhxc_index = jsub
                   call integrals_locpot(ham%lhxc(is)%m(isub,jsub), &
                        rep%ngwfs_on_grid(isub),ngwf_basis(isub), &
                        rep%ngwfs_on_grid(jsub),ngwf_basis(jsub), &
                        mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, &
                        lhxc_fine(:,:,:,is,lhxc_index))
                enddo
             enddo

             ! agrecokpt: add TB method extra terms if requested
             if ((pub_kpoint_method == 'TB') .and. (loc_kpt_point%x/=0.0_DP.or.&
                  loc_kpt_point%y/=0.0_DP.or.loc_kpt_point%z/=0.0_DP)) then
                if (pub_debug_on_root) write(stdout,'(a,i3)') 'DEBUG: &
                   &Adding TB terms in lhxc matrix spin ', is

                ! rc2013: abort if using multiple subsystems
                call utils_assert(mdl%nsub == 1, &
                     'Error in electronic_init_denskern: non-Gamma &
                     &k-points not allowed with more than 1 subsystem.')
                call hamiltonian_apply_phases_tb_method(ham%lhxc(is)%p, &
                     ngwf_basis(1), mdl, loc_kpt_point)

             end if

             ! agrecocmplx: check lhxc matrix is hermitian if requested
             if (pub_check_hermitian) then
                rms_hermitian_check = sparse_embed_check_hermitian(ham%lhxc(is))

                if (pub_on_root) write(stdout,'(a,i3,a,f8.5)') 'RMS lhxc spin',&
                   is, ' hermitian check: ', rms_hermitian_check

                call utils_assert(rms_hermitian_check < pub_imag_thr, &
                     'WARNING: hermitian check for lhxc matrix above threshold')

             end if

             ! ham = Hartree + locpot
             ! agrecokpt: TB method extra terms already included in each term
             ! at this point
             call sparse_embed_copy(ham%ham(is),ham%lhxc(is))

             !!! JCW: At present (01/2016), the initial ke_density_fine used
             !!! JCW: by xc_energy_potential in electronic_init_pot, is
             !!! JCW: ke_density_fine = 0.0_DP.
             !!! JCW: Consequently, dfdtau_fine should be 0.0_DP at all grid points
             !!! JCW: and the dfdtau integrals will consequently be all zero.
             !!! JCW: Until a non-zero initial guess for ke_density_fine is
             !!! JCW: implemented, the dfdtau integrals do not need to be
             !!! JCW: evaluated in electronic_init_denskern.
             !!! JCW: Uncomment the following if/when this is implemented:
             !if (pub_xc_ke_density_required) then
             !   call utils_assert(present(dfdtau_fine),"Error in &
             !        &electronic_init_denskern: KE-density dependent XC &
             !        &functional requested, but dfdtau_fine not present.")
             !   ! JCW: Calculate dfdtau integrals if KE-density-dependent XC
             !   ! JCW: functional is used.
             !   call integrals_div_locpot_grad_functions(ham%dfdtau(is),&
             !        rep%ngwfs_on_grid,ngwf_basis,rep%ngwfs_on_grid,ngwf_basis, &
             !        mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, &
             !        dfdtau_fine(:,:,:,is))
             !!    JCW: Add dfdtau matrix to Hamiltonian (correction for
             !!    JCW: KE-density-dependent part of XC functional
             !    if (pub_debug_on_root) then
             !      write(stdout,'(a)') "DEBUG: Adding dfdtau matrix to Hamiltonian"
             !   end if
             !   call sparse_axpy(ham%ham(is),ham%dfdtau(is),1.0_DP)
             !end if
             !!! NOTE ON INITIAL GUESS XC CONTRIBUTION
             !!! JCW: Unless XC_INITIAL_FUNCTIONAL keyword is set, the initial
             !!! JCW: guess will not include contributions from XC for meta-GGAs.

             ! ham = (Hartree + locpot) + kinet
             call sparse_embed_axpy(ham%ham(is),rep%kinet,1.0_DP)

             ! ham = (Hartree + locpot + kinet) + nonlocpot
             if (pub_any_nl_proj.and.(.not.pub_aug)) then
                call sparse_embed_axpy(ham%ham(is),rep%nonlocpot,1.0_DP)
             end if
             if (pub_aug) then
                call sparse_embed_axpy(ham%ham(is),ham%nonlocpot(is),1.0_DP)
             end if

             ! ndmh: add the hubbard hamiltonian if necessary (spin-splitting only)
             if (pub_hubbard) call sparse_embed_axpy(ham%ham(is),ham%hubbard_ham(is), &
                  1.0_DP)

          end do

          ! lr408: If this is a conduction call, we need to calculate the
          ! lr408: projected conduction Hamiltonian
          if (cond_call) then
             ! KPOINTS_DANGER: only 1 k-point considered
             call hamiltonian_proj_cond_matrix(rep, ham, &
                  val_rep, val_ham, val_dkn%m(:,PUB_1K))
          end if
       end if

       do is=1,pub_num_spins

          ! KPOINTS_DANGER: only one k-point considered!
          if (cond_call) then
             if (pub_on_root) then
                write(stdout,'(a,I5,a)') 'Generating conduction density kernel &
                     &for ',rep%n_occ(is,PUB_1K),' ''occupied'' orbitals'
             end if
          end if

          if (pub_on_root) then
             if (is == UP) then
                write(stdout,'(a)',advance='no') &
                     'Up spin density kernel initialisation ...'
             elseif (is == DN) then
                write(stdout,'(a)',advance='no') &
                     'Down spin density kernel initialisation ...'
             endif
          endif

          ! mjsp: If SCF-MI, zero interfragmental Hamiltonian blocks to
          ! ensure initialised fragment density kernel is fragment localised.
          if (pub_eda_scfmi) then
             if (pub_eda_mode == EDA_CTFRAGLOC_DEVEL) then
               ! include certain interfragmental blocks if a delocalisation
               ! calculation
               call fmo_freeze_mat_spam_nxn(ham%ham(is),pub_frag_counter, &
                     pub_frag_counter2)
             else
               call fmo_freeze_mat_spam_nxn(ham%ham(is),0)
             end if
          end if

          ! KPOINTS_DANGER: only first k-point considered!
          if ((.not.pub_edft).and.(.not.(pub_edft_init_maxit>0))) then
             if (pub_maxit_palser_mano < 1) then
                ! Diagonalise the initial Hamiltonian
                ! agrecokpt: k-point dependence already included in hamiltonian
                ! and overlap, no need to explicitly include it in density kernel
                call kernel_init_core_ham(denskern%kern%m(is,PUB_1K), &
                     ham%ham(is), rep%overlap, rep%n_occ(is,PUB_1K))
             else
                ! Palser-Manolopoulos Canonical Purification
                call palser_mano_kernel_optimise(denskern%kern%m(is,PUB_1K), &
                     ham%ham(is), rep%overlap, rep%inv_overlap, &
                     rep%n_occ(is,PUB_1K), &
                     num_iter=pub_maxit_palser_mano)
             endif
          else if (pub_edft_init_maxit>0) then
             if (pub_dftb) then
                call edftb_calculate(edft_for_init, denskern%kern, ham, rep, &
                     ngwf_basis, mdl, pub_edft_init_maxit)
             else
                call edft_calculate(edft_for_init, denskern%kern, ham, rep, &
                     lhxc_fine, ngwf_basis, hub_proj_basis, hub, mdl, hfxstate, &
                     pub_edft_init_maxit)
             end if
          end if

          if (pub_on_root) write(stdout,'(a)') '... done'

       end do

       ! cks: output density kernel to file if this is requested
       if(.not.loc_force_no_IO) then
          if (pub_write_denskern.and.(.not.pub_write_converged_dk_ngwfs)) &
               call restart_kernel_write(denskern%kern,write_cond=cond_call)
       end if

       ! --- END GUESS A DENSITY KERNEL ----------------------------

    else if (.not.pub_coreham_denskern_guess .and. initialise_new) then

       ! --- SCALE INVERSE OVERLAP MATRIX --------------------------

       ! KPOINTS_DANGER: this only deals with 1 k-point
       do is=1,pub_num_spins
          filling = real(rep%n_occ(is,PUB_1K),DP) / &
                      real(sum(ngwf_basis(:)%num),DP)
          ! agrecokpt: k-point dependence already included in overlap matrix,
          ! hence no need to explicitly include it in density kernel
          call sparse_embed_copy(denskern%kern%m(is,PUB_1K),rep%inv_overlap)
          call sparse_embed_scale(denskern%kern%m(is,PUB_1K),filling)
          if (pub_output_detail >= VERBOSE .and. pub_on_root) &
               write(stdout,'(a,i1,a,f6.3)') &
               ' Initialising density kernel for spin ',is, &
               ' to inverse overlap scaled by filling ',filling
       end do

       ! --- END SCALE INVERSE OVERLAP MATRIX ----------------------

    end if

    ! agrecocmplx: check density kernel is hermitian if requested
    ! (before kernel optimisation)
    if (pub_check_hermitian) then
       do is=1,pub_num_spins
          rms_hermitian_check = &
             sparse_embed_check_hermitian(denskern%kern%m(is,PUB_1K))

          if (pub_on_root) write(stdout,'(a,i3,a,f8.5)') 'RMS denskern spin ',&
             is,' hermitian check (before optimisation): ', rms_hermitian_check

          call utils_assert(rms_hermitian_check < pub_imag_thr, &
             'WARNING: hermitian check for denskern matrix above threshold')

       end do
    end if


    ! &&&&&&&& USE PENALTY FUNCTIONAL TO IMPROVE DENSITY KERNEL &&&&
    ! ddor: This is unnecessary when restarting with pub_task HUBBARDSCF
    ! lr408: or when restarting with pub_task COND
    if ((pub_maxit_pen > 0) .and. (.not. pub_hubbard_restart ) .and. &
         (.not. pub_cond_calculate) .and. (.not.pub_edft)) then

       call penalty_denskernel_optimise_cg(denskern, rep, ham, &
            lhxc_fine, Ecurr, &
            ngwf_basis, hub_proj_basis, hub, mdl, hfxstate, &
            lnv_threshold_orig)

       ! agrecocmplx: check density kernel is hermitian if requested
       ! (after kernel optimisation)
       if (pub_check_hermitian) then
          do is=1,pub_num_spins
             rms_hermitian_check = &
                sparse_embed_check_hermitian(denskern%kern%m(is,PUB_1K))

             if (pub_on_root) write(stdout,'(a,i3,a,f8.5)') 'RMS denskern spin ',&
                is,' hermitian check (after optimisation): ', rms_hermitian_check

             call utils_assert(rms_hermitian_check < pub_imag_thr, &
                'WARNING: hermitian check for denskern matrix above threshold')

          end do
       end if

    end if
    ! -&&&&& END USE PENALTY FUNCTIONAL TO IMPROVE DENSITY KERNEL &&&&-

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &electronic_init_denskern'

    call services_flush

  end subroutine electronic_init_denskern


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine electronic_init_pot(lhxc_fine, dijhat, mdl, force_init_pot, dfdtau_fine)

    !===============================================================!
    ! This subroutine initialises the local pseudopotential, the    !
    ! core density, the initial electron density, and uses them to  !
    ! initialise the initial-guess Hamiltonian if required.         !
    !---------------------------------------------------------------!
    ! This subroutine was created by Nicholas Hine in October 2010. !
    ! It was based around parts of the previous routine             !
    ! electronic_init_denskern_pot, which was written by the ODG in !
    ! 2000-2010.                                                    !
    ! Support for initial kinetic energy density added by           !
    ! James C. Womack                                               !
    ! Edited for multiple parallel strategies and subsystem calcs   !
    ! by Robert Charlton, January 2017.                             !
    !===============================================================!

    use augmentation, only: augmentation_screen_dij
    use classical_pot, only: classical_pot_struct_fac
    use comms, only: pub_on_root
    use constants, only: DP, VERBOSE, stdout, ANGSTROM
    use cutoff_coulomb, only: cutoff_coulomb_struct_fac, &
         cutoff_coulomb_localpseudo, cutoff_coulomb_hartree, &
         cutoff_coulomb_initial_guess, cutoff_coulomb_core_density, &
         cutoff_coulomb_classical_sfac
    use density_init, only: density_init_guess_recip, density_init_guess_real
    use function_basis, only: FUNC_BASIS
    use hartree, only: hartree_on_grid, hartree_via_multigrid
    use is_smeared_ions, only: smeared_ion_apply_vloc_corr
    !use ke_density, only: weizsacker_ke_density_on_grid ! JCW: experimental
    use model_type, only: MODEL
    use openbc_locps, only: openbc_locps_on_grid
    use paw, only: paw_tcore_hartree_on_grid, paw_tcore_density, &
         paw_tcore_density_recip
    use potential, only: potential_sawtooth_efield
    use pseudopotentials, only: pseudo_make_structure_factor,&
         pseudopotentials_local_on_grid, pseudopotentials_core_density
    use rundat, only: pub_coulomb_cutoff, pub_output_detail,  &
         pub_coreham_denskern_guess, pub_constant_efield, &
         pub_nlcc, pub_paw, pub_read_denskern, pub_multigrid_hartree, &
         pub_nhat_in_xc, pub_initial_dens_realspace, pub_num_spins, &
         pub_turn_off_hartree, pub_is_smeared_ion_rep, pub_spin_fac, &
         pub_pspot_bc_is_periodic, &
         pub_xc_initial_functional, pub_devel_code, &
         pub_xc_ke_density_required, pub_emft, pub_emft_follow, &
         pub_active_region, pub_inner_loop_iteration
    use services, only: services_flush
    use sparse, only: SPAM3, sparse_copy, sparse_create, sparse_destroy
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy, &
        sparse_embed_create, sparse_embed_destroy
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_abort, utils_devel_code
#ifdef TESTCOREDEN
    ! ndmh 28/08/15: to be removed once there is no further risk of any
    ! ndmh 28/08/15: bugs being found in the new core density code.
    use visual, only: visual_scalarfield
#endif
    use xc, only: xc_energy_potential, xc_emft_calculate

    implicit none

    ! Arguments (inputs)
    type(MODEL), intent(inout) :: mdl

    ! Arguments (outputs)
    real(kind=DP), intent(out) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub)
    type(SPAM3_EMBED), intent(inout) :: dijhat(pub_num_spins)

    logical, intent(in) :: force_init_pot
    real(kind=DP), intent(out), optional :: dfdtau_fine(:,:,:,:)
    ! JCW: In general, we have
    ! JCW:    dfdtau_fine(mdl%fine_grid%ld1,&
    ! JCW:    mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    ! JCW: or
    ! JCW:    dfdtau_fine(1,1,1,1)
    ! JCW: when it is unneeded.

    ! Local Variables
    integer :: is
    integer :: ierr
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: struct_fac
    complex(kind=DP), allocatable, dimension(:,:,:) :: struct_fac_classical !ars
    real(kind=DP), allocatable, dimension(:,:,:,:) :: density_fine
    real(kind=DP), allocatable, dimension(:,:,:,:) :: sub_density_fine
    real(kind=DP), allocatable, dimension(:,:,:,:) :: active_density_fine
    real(kind=DP), allocatable, dimension(:,:,:,:) :: pot_fine
    real(kind=DP), allocatable, dimension(:,:,:,:) :: pot_fine_emft
    ! JCW: ke_density_fine contains kinetic energy density
    real(kind=DP), allocatable, dimension(:,:,:,:) :: ke_density_fine
    ! jcap: regional localpseudo_fine
    real(kind=DP), allocatable, dimension(:,:,:) :: reg_localpseudo_fine
    real(kind=DP) :: xc_energy
    real(kind=DP) :: xc_energy_emft
    integer :: isub, num_pspecies
    real(kind=DP) :: debug_val
    ! rc2013: temporary matrix to hold dijhat info; avoids array temporaries
    type(SPAM3), allocatable :: dij_array(:)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &electronic_init_pot'

    ! Allocate storage for structure factor(s) and calculate it
    if (.not.pub_coulomb_cutoff) then

       allocate(struct_fac(mdl%num_pspecies, mdl%fine_grid%ld3, &
            mdl%fine_grid%ld2, mdl%fine_grid%max_slabs23, mdl%nsub),stat=ierr)
    else
       ! allocate padded structure factor
       allocate(struct_fac(mdl%num_pspecies, mdl%padded_grid%ld3, &
            mdl%padded_grid%ld2, mdl%padded_grid%max_slabs23, mdl%nsub),stat=ierr)
    end if
    call utils_alloc_check('electronic_init_pot','struct_fac',ierr)

    ! rc2013: allocate the array to hold the dijhat components
    allocate(dij_array(pub_num_spins),stat=ierr)
    call utils_alloc_check('electronic_init_pot','dij_array',ierr)

    ! rc2013: initialise to 0 before adding regional contributions
    mdl%localpseudo_fine = 0.0_DP
    ! jcap: allocate and deallocate temporary regional
    ! localpseudo_fine here to save memory - not used elsewhere
    allocate(reg_localpseudo_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12),stat=ierr)
    call utils_alloc_check('electronic_init_pot','reg_localpseudo_fine', ierr)
    do isub=1,mdl%nsub
       num_pspecies=mdl%regions(isub)%par%num_pspecies
       ! Allocate storage for structure factor(s) and calculate it
       if (.not.pub_coulomb_cutoff) then

          ! calculate structure factor for all species
          ! jcap: only send in the appropriate amount of the
          ! struct_fac array
          call pseudo_make_structure_factor(struct_fac(1:num_pspecies,:,:,:,isub), &    ! output
               mdl%regions(isub)%elements,mdl%fine_grid, &
               num_pspecies)    ! input

          ! cks: allocate structure factor for classical atoms if it is needed
          if (mdl%nat_classical.gt.0) then

             allocate(struct_fac_classical(mdl%fine_grid%ld3,mdl%fine_grid%ld2, &
                  mdl%fine_grid%max_slabs23),stat=ierr)
             call utils_alloc_check('electronic_init_pot', &
                  'struct_fac_classical',ierr)

             ! ars: calculate structure factor for the classical atoms
             call classical_pot_struct_fac(struct_fac_classical, &
                  mdl%classical_elements,mdl%nat_classical,mdl%fine_grid)
          else
             allocate(struct_fac_classical(0,0,0),stat=ierr)
             call utils_alloc_check('electronic_init_pot', &
                  'struct_fac_classical',ierr)
          end if

          ! calculate the local pseudopotential on the fine grid
          mdl%locps_gzero_term = 0.0_DP
          if (.not.any(pub_pspot_bc_is_periodic(1:3))) then
             ! JCW: If fully open, use open BC local pseudopotential routine
             ! rc2013: pass info specific to this region
             call openbc_locps_on_grid(reg_localpseudo_fine, &
                  mdl%fine_grid,mdl%cell,mdl%regions(isub)%elements, &
                  mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp, &
                  mdl%regions(isub)%par)
             ! jd: Include the correction due to smeared ions, if any, in v_local_fine
             if (pub_is_smeared_ion_rep) then
                call smeared_ion_apply_vloc_corr(mdl%regions(isub)%elements, &
                     reg_localpseudo_fine, mdl%fine_grid, &
                     pub_pspot_bc_is_periodic, mdl%regions(isub)%par)
             end if
          else if (all(pub_pspot_bc_is_periodic(1:3))) then
             ! JCW: If fully periodic, use periodic local pseudopotential routine
             if (.not.pub_paw) then
                call pseudopotentials_local_on_grid(reg_localpseudo_fine, &
                     mdl%locps_gzero_term, struct_fac(1:num_pspecies,:,:,:,isub),&
                     struct_fac_classical, &
                     mdl%fine_grid,mdl%cell,mdl%regions(isub)%pseudo_sp, &
                     mdl%regions(isub)%par, mdl%nat_classical)
             else
                call paw_tcore_hartree_on_grid(reg_localpseudo_fine, &
                     mdl%locps_gzero_term, struct_fac(1:num_pspecies,:,:,:,isub), &
                     struct_fac_classical, &
                     mdl%regions(isub)%paw_sp,mdl%fine_grid,mdl%regions(isub)%par)
             end if

             ! Include the correction due to smeared ions, if any, in v_local_fine
             if (pub_is_smeared_ion_rep) then
                call smeared_ion_apply_vloc_corr(mdl%regions(isub)%elements, &
                     reg_localpseudo_fine, mdl%fine_grid, &
                     pub_pspot_bc_is_periodic, mdl%regions(isub)%par)
             end if

          else
             ! JCW: If mixed periodic/open, then abort (not supported outside of
             ! JCW: cutoff Coulomb currently)
             call utils_abort("Error in electronic_init_pot: &
                  &Mixed BCs for non-cutoff-Coulomb calculations are not yet supported.")
          end if

       else ! Cutoff Coulomb - use padding wrappers

          ! calculate structure factor for all species
          call cutoff_coulomb_struct_fac(struct_fac(1:num_pspecies,:,:,:,isub), &   ! output
               mdl, isub)                                      ! input

          ! ndmh: allocate structure factor for classical atoms if it is needed
          mdl%locps_gzero_term = 0.0_DP
          if (mdl%nat_classical.gt.0) then
             allocate(struct_fac_classical(mdl%padded_grid%ld3, &
                  mdl%padded_grid%ld2, mdl%padded_grid%max_slabs23),stat=ierr)
             call utils_alloc_check('electronic_init_pot', &
                  'struct_fac_classical',ierr)

             ! ndmh: wrapper to calculate structure factor for the classical atoms
             call cutoff_coulomb_classical_sfac(struct_fac_classical,mdl) ! output
          else
             allocate(struct_fac_classical(0,0,0),stat=ierr)
             call utils_alloc_check('electronic_init_pot', &
                  'struct_fac_classical',ierr)
          endif

          ! calculate the local pseudopotential on the fine grid
          call cutoff_coulomb_localpseudo(reg_localpseudo_fine, &
               mdl%locps_gzero_term,struct_fac(1:num_pspecies,:,:,:,isub),&
               struct_fac_classical,mdl,isub)

       end if
       ! rc2013: build up total localpseudo_fine from subsystems
       mdl%localpseudo_fine = mdl%localpseudo_fine + reg_localpseudo_fine
       deallocate(struct_fac_classical,stat=ierr)
       call utils_dealloc_check('electronic_init_pot','struct_fac_classical',ierr)
    end do
    deallocate(reg_localpseudo_fine,stat=ierr)
    call utils_dealloc_check('electronic_init_pot','reg_localpseudo_fine',ierr)

    ! ndmh: calculate the core charge density
    if (pub_nlcc) then

       if (pub_on_root .and. pub_output_detail>=VERBOSE) then
          write(stdout,'(a)',advance='no') 'Calculating core density ...'
       end if

       ! jcap: for each subsystem
       mdl%core_density_fine = 0.0_DP
       do isub=1,mdl%nsub
          num_pspecies=mdl%regions(isub)%par%num_pspecies
          ! ndmh: use padding wrapper if using coulomb cutoff
          if (.not.pub_coulomb_cutoff) then
             if (.not.pub_paw) then
                call pseudopotentials_core_density(mdl%regions(isub)%core_density_fine,&
                     struct_fac(1:num_pspecies,:,:,:,isub), &
                     mdl%regions(isub)%pseudo_sp,mdl%fine_grid, &
                     mdl%regions(isub)%par)
             else
#ifdef TESTCOREDEN
                call paw_tcore_density(mdl%regions(isub)%core_density_fine, mdl%fine_grid, &
                     mdl%cell,mdl%regions(isub)%paw_sp,mdl%regions(isub)%par)
                call visual_scalarfield(mdl%regions(isub)%core_density_fine, mdl%fine_grid, &
                     mdl%cell, 'Core density (in e/ang^3) for:', '_coredensitynew', &
                     mdl%regions(isub)%elements, ANGSTROM**3)
                call paw_tcore_density_recip(mdl%regions(isub)%core_density_fine, &
                     struct_fac(1:num_pspecies,:,:,:,isub), &
                     mdl%regions(isub)%paw_sp,mdl%fine_grid,mdl%regions(isub)%par)
                call visual_scalarfield(mdl%regions(isub)%core_density_fine, mdl%fine_grid, &
                     mdl%cell, 'Core density (in e/ang^3) for:', '_coredensityold', &
                     mdl%regions(isub)%elements, ANGSTROM**3)
#else
                call paw_tcore_density(mdl%regions(isub)%core_density_fine, mdl%fine_grid, &
                     mdl%cell,mdl%regions(isub)%paw_sp,mdl%regions(isub)%par)
#endif
             end if
          else
             ! jcap: will need to edit this to make it more consistent - replace last argument
             call cutoff_coulomb_core_density(mdl%regions(isub)%core_density_fine, &
                  struct_fac(1:num_pspecies,:,:,:,isub),mdl,mdl%regions(isub))
          end if

          ! jcap: sum up to get full system core density
          mdl%core_density_fine = mdl%core_density_fine + mdl%regions(isub)%core_density_fine

          ! jcap: with the exception of the active region in an
          ! embedding calculation, regions%core_density_fine is not
          ! used again, so deallocate
          if ( .not. ((pub_emft.or.pub_emft_follow) .and. (isub.eq.pub_active_region)) ) then
             deallocate(mdl%regions(isub)%core_density_fine,stat=ierr)
             call utils_dealloc_check('electronic_init_pot', &
                  'mdl%regions%core_density_fine', ierr, &
                  allocated_in = 'energy_and_force_calculate')
          end if

       end do ! isub

       if (pub_on_root .and. pub_output_detail>=VERBOSE) then
          write(stdout,'(a)') ' ... done'
       end if

    end if

    ! cks: apply constant electric field as sawtooth contribution to local
    ! cks: potential
    if (any(pub_constant_efield(:) /= 0.0_DP)) then
       call potential_sawtooth_efield(mdl%localpseudo_fine,mdl%fine_grid, &
            mdl%cell)
    end if

    if ((pub_coreham_denskern_guess).and.((.not.pub_read_denskern).or.force_init_pot)) then

       ! Allocate temporary arrays for potential and density
       allocate(density_fine(mdl%fine_grid%ld1, &
            mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins), stat=ierr)
       call utils_alloc_check('electronic_init_pot','density_fine',ierr)
       if (pub_emft) then
          allocate(active_density_fine(mdl%fine_grid%ld1, &
               mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins), &
               stat=ierr)
          call utils_alloc_check('electronic_init_pot','active_density_fine',&
               ierr)
       end if

       ! JCW: Allocate temporary arrays for kinetic energy density, if needed
       if (pub_xc_ke_density_required) then
          allocate(ke_density_fine(mdl%fine_grid%ld1, &
               mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, &
               pub_num_spins),stat=ierr)
          call utils_alloc_check('electronic_init_pot','ke_density_fine',ierr)
       else
          ! JCW: If no tau-dependent XC functional is requested, then allocate a
          ! JCW: 1 element dummy array (this avoids the need for conditional
          ! JCW: subroutine calls and optional arguments)
          allocate(ke_density_fine(1,1,1,1),stat=ierr)
          call utils_alloc_check('electronic_init_pot','ke_density_fine',ierr)
       end if

       ! Calculate initial density guess if required
       if (pub_initial_dens_realspace) then
          if (pub_emft) then
             call density_init_guess_real(density_fine, mdl, &
                  mdl%fine_grid, add_aug_den=.true.,&
                  active_density_on_grid=active_density_fine)  ! jcap
             allocate(pot_fine_emft(mdl%fine_grid%ld1,&
                  mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins),&
                  stat=ierr)
             call utils_alloc_check('electronic_init_pot','pot_fine_emft',ierr)
          else
             call density_init_guess_real(density_fine, mdl, &
               mdl%fine_grid, add_aug_den=.true.)
          end if
          allocate(pot_fine(mdl%fine_grid%ld1,&
               mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins), &
               stat=ierr)
          call utils_alloc_check('electronic_init_pot','pot_fine',ierr)
       else

          density_fine=0.d0
          do isub=1,mdl%nsub
             num_pspecies=mdl%regions(isub)%par%num_pspecies
             ! First, allocate subregion density
             allocate(sub_density_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
                  mdl%fine_grid%max_slabs12,pub_num_spins), stat=ierr)
             call utils_alloc_check('electronic_init_pot','sub_density_fine',&
                  ierr)

             if (pub_coulomb_cutoff) then
                call cutoff_coulomb_initial_guess(sub_density_fine, mdl, &
                     struct_fac(1:num_pspecies,:,:,:,isub), isub)
             else
                call density_init_guess_recip(sub_density_fine, &
                     mdl%regions(isub), struct_fac(1:num_pspecies,:,:,:,isub), &
                     mdl%fine_grid)
             end if

             density_fine=density_fine+sub_density_fine
             if (pub_emft.and.(isub==pub_active_region)) active_density_fine=&
                  sub_density_fine
          end do
       end if

       ! JCW: Initialize kinetic energy density, if required.
       if (pub_xc_ke_density_required) then
          ! JCW: Initial guess is zero
          ke_density_fine = 0.0_DP
          ! JCW: Weizsacker KE density initial guess [Experimental]
          ! JCW: Use von Weizacker functional for initial guess
          ! JCW:   \tau^{W} = (1/8) * |grad(rho)|^{2} / rho
          ! JCW: Warning: using Weizscaker KE density for initial guess is
          ! JCW:          untested, as is the routine
          ! JCW:          weizsacker_ke_density_on_grid
          !call weizsacker_ke_density_on_grid(ke_density_fine, &
          !     density_fine, mdl%fine_grid )
       end if

       ! cks: density_fine_init is converted to Hartree potential
       ! jd: Calculate either:
       !     (1) the usual Hartree potential
       !     (2) Hartree potential with cutoff Coulomb
       !     (3) Hartree potential with multigrid, which could be
       !         a) electronic Hartree potential, with open BCs,
       !         b) molecular Hartree potential, with smeared-ions, in vacuo,
       !         c) molecular Hartree potential, in implicit solvent.
       !     All of the above are mutually exclusive.
       !     Iff (3c) and the dielectric is not fixed, but instead depends on
       !     the electronic density, an extra term caused by this dependence
       !     (V_eps(r)+V_apolar(r)) is included in the potential to get the
       !     correct functional derivative. The Hartree energy then needs to be
       !     calculated differently, and hartree_via_multigrid() takes care of
       !     that. However, here the actual energy value is not needed.

       if (.not. pub_turn_off_hartree) then
          if (pub_coulomb_cutoff) then
             call cutoff_coulomb_hartree(lhxc_fine(:,:,:,:,1), density_fine,mdl)      ! (2)
          else if (pub_multigrid_hartree) then
             call hartree_via_multigrid(lhxc_fine(:,:,:,:,1), density_fine, &
                  mdl%fine_grid, mdl%cell, elements=mdl%elements)           ! (3)
          else
             call hartree_on_grid(lhxc_fine(:,:,:,:,1), density_fine, mdl%fine_grid, &
                  mdl%cell)  ! (1)
             debug_val =  maxval(density_fine)
          end if
          ! rc2013: the Hartree potential will be the same for all regions
          do isub=1,mdl%nsub
             lhxc_fine(:,:,:,:,isub) = lhxc_fine(:,:,:,:,1)
          enddo
       else
          lhxc_fine = 0.0_DP
       endif

       if (pub_xc_ke_density_required.and.pub_xc_initial_functional=="NONE") then
          ! JCW: If a KE-density-dependent meta-GGA XC functional is used AND
          ! JCW: an alternative functional for the initial guess is not specified:
          ! JCW: Do not include the XC potential, since we currently do not have
          ! JCW: a good initial KE density guess. Using ke_density_fine = 0.0_DP
          ! JCW: causes discrepancies between different functional implementations
          ! JCW: due to different treatments of the condition where tau = 0.0 and
          ! JCW: den > 0.0.
          if (pub_on_root) then
             write(stdout,*) "WARNING: XC potential not used to generate initial &
                  &density guess when a meta-GGA exchange-correlation functional &
                  &is used and no alternative initial functional is specified. &
                  &The Hamiltonian is constructed without the XC contribution."
          end if
       else
          ! ndmh: Add up local pseudo + Hartree
          do is=1,pub_num_spins
             do isub=1,mdl%nsub
                lhxc_fine(:,:,:,is,isub) = lhxc_fine(:,:,:,is,isub) &
                     + mdl%localpseudo_fine
             enddo
          end do

          ! ndmh: screen nonlocal projector energies now if the augmentation density
          ! ndmh: is not going to be included in the exchange energy calculation
          if (pub_paw.and.(.not.pub_nhat_in_xc)) then
             ! rc2013: EMBED_WARNING!
             do isub=1,mdl%nsub
                do is=1,pub_num_spins
                   call sparse_create(dij_array(is), dijhat(is)%m(isub,isub))
                   call sparse_copy(dij_array(is), dijhat(is)%m(isub,isub))
                end do
                call augmentation_screen_dij(dij_array, &
                     lhxc_fine(:,:,:,:,isub),mdl%aug_box,mdl%cell,mdl%fine_grid, &
                     mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp)
                do is=1,pub_num_spins
                   call sparse_copy(dijhat(is)%m(isub,isub), dij_array(is))
                   call sparse_destroy(dij_array(is))
                end do
             end do
          end if

          ! If we have a realistic density guess, we can include the XC potential
          if (pub_initial_dens_realspace) then
             if (pub_paw.and.(.not.pub_nhat_in_xc)) then
                ! Re-generate density but leave out augmentation density this time
                if (pub_emft) then
                   call density_init_guess_real(density_fine, mdl, &
                        mdl%fine_grid, add_aug_den=.false., &
                        active_density_on_grid=active_density_fine)  ! jcap
                else
                   call density_init_guess_real(density_fine, mdl, &
                        mdl%fine_grid, add_aug_den=.false.)
                end if
                ! Otherwise keep previous density
             end if

             ! ndmh: add on core density before calculating xc potential for NLCC
             if (pub_nlcc) then
                do is=1,pub_num_spins
                   density_fine(:,:,:,is) = density_fine(:,:,:,is) + &
                        mdl%core_density_fine * 0.5_DP * pub_spin_fac
                   if (pub_emft) active_density_fine(:,:,:,is) = &
                        active_density_fine(:,:,:,is) + &
                        mdl%core_density_fine * 0.5_DP * pub_spin_fac
                end do
             end if

             ! ndmh: no augmentation density to pass in
             if (index(pub_devel_code,'NO_XC_INIT')>0) then
                ! If devel_code: NO_XC_INIT then, skip XC potential
                ! evaluation in initial guess
                if (pub_on_root) write(stdout,*) "WARNING: NO_XC_INIT, &
                   &electronic_init_pot setting pot_fine = 0.0_DP"
                xc_energy = 0.0_DP
                pot_fine = 0.0_DP
             else if (pub_xc_ke_density_required) then
                ! JCW: Uncomment when a non-zero KE density initial guess is
                ! JCW: available
                !call xc_energy_potential(density_fine, xc_energy, pot_fine, &
                !  mdl%fine_grid,mdl%cell,0, &
                !  ke_density_fine = ke_density_fine, &
                !  dfdtau_fine = dfdtau_fine )

                ! JCW: WORKAROUND while KE density initial guess not available:
                ! JCW: Use the alternative functional set in
                ! JCW: pub_xc_initial_functional
                ! jcap: This will need looking at for EMFT embedding
                ! once the workaround is dealt with
                call xc_energy_potential(density_fine, xc_energy, pot_fine, &
                     mdl%fine_grid,mdl%cell,0,initial=.true.)
                dfdtau_fine = 0.0_DP
             else
                call xc_energy_potential(density_fine, xc_energy, pot_fine, &
                     mdl%fine_grid,mdl%cell,0)

                ! jcap: If we are using EMFT, we need to
                ! call this twice more, in order to get the xc energy and
                ! potential for the active region
                if (pub_emft) call xc_emft_calculate(active_density_fine, &
                     xc_energy_emft, pot_fine_emft, mdl%fine_grid, mdl%cell)
             end if

             ! rc2013: with EMFT each region will have different potentials
             do isub=1,mdl%nsub
                lhxc_fine(:,:,:,:,isub) = lhxc_fine(:,:,:,:,isub) + pot_fine
                if(isub == pub_active_region .and. pub_emft) &
                     lhxc_fine(:,:,:,:,isub) = lhxc_fine(:,:,:,:,isub) + pot_fine_emft
             enddo

             ! Deallocate temporary density array
             deallocate(pot_fine,stat=ierr)
             call utils_dealloc_check('electronic_init_pot','pot_fine',ierr)
             if (pub_emft) then
                deallocate(pot_fine_emft,stat=ierr)
                call utils_dealloc_check('electronic_init_pot','pot_fine_emft',&
                     ierr)
             end if

          end if
       end if

       ! ndmh: screen nonlocal projector energies now if the augmentation density
       ! ndmh: is going to be included in the exchange energy calculation
       if (pub_paw.and.pub_nhat_in_xc) then
          ! rc2013: EMBED_WARNING!
          do isub=1,mdl%nsub
             do is=1,pub_num_spins
                call sparse_create(dij_array(is), dijhat(is)%m(isub,isub))
                call sparse_copy(dij_array(is), dijhat(is)%m(isub,isub))
             end do
             call augmentation_screen_dij(dij_array, &
                  lhxc_fine(:,:,:,:,isub),mdl%aug_box,mdl%cell,mdl%fine_grid, &
                  mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp)
             do is=1,pub_num_spins
                call sparse_copy(dijhat(is)%m(isub,isub), dij_array(is))
                call sparse_destroy(dij_array(is))
             end do
          end do
       end if

       ! JCW: Deallocate kinetic energy density on grid
       deallocate(ke_density_fine,stat=ierr)
       call utils_dealloc_check('electronic_init_pot','ke_density_fine',ierr)
       deallocate(density_fine,stat=ierr)
       call utils_dealloc_check('electronic_init_pot','density_fine',ierr)

    end if

    ! Deallocate storage for structure factor(s)
    deallocate(struct_fac,stat=ierr)
    call utils_dealloc_check('electronic_init_pot','struct_fac',ierr)

    ! rc2013: deallocate the array to hold the dijhat components
    deallocate(dij_array,stat=ierr)
    call utils_dealloc_check('electronic_init_pot','dij_array',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &electronic_init_pot'

  end subroutine electronic_init_pot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module electronic_init
