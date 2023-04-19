! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=================================================================!
!                                                                 !
!         Conduction States Calculation Module                    !
!                                                                 !
! This module optimises a set of NGWFs to describe conduction     !
! states and calculates the properties of those states.           !
!-----------------------------------------------------------------!
! Written by Laura Ratcliff on 14/10/2010                         !
! Minor improvements and cleanup by Nicholas Hine, October 2010.  !
! Modified for embedding by Joseph Prentice, June 2018            !
!=================================================================!

module conduction_properties

  use constants, only: DP

  implicit none

  private

  public :: conduction_properties_calculate

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine conduction_properties_calculate(val_denskern, cond_denskern, &
        val_rep, cond_rep, val_ham, cond_ham, val_ngwf_basis, &
        proj_basis, nl_projectors, hub_proj_basis, hub, cond_ngwf_basis, &
        joint_ngwf_basis, core_basis, core_wvfns, mdl, hfxstate, lhxc_fine, &
        kpt) ! agrecokpt

    !========================================================================!
    ! This subroutine calculates properties for a conduction calculation.    !
    !------------------------------------------------------------------------!
    ! Written by Laura Ratcliff in September 2011.                           !
    ! Modified for embedding by Joseph Prentice, June 2018                   !
    ! Extended by Jacek Dziedzic in July 2018 to support hybrid conduction.  !
    !========================================================================!

    use augmentation, only: augmentation_overlap
    use datatypes, only: data_functions_alloc, data_functions_dealloc
    use comms, only: comms_barrier, pub_on_root
    use constants, only: stdout, NORMAL, REP_SWEX_HFX_OTHER
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_brappd_ketppd
    use hamiltonian, only: hamiltonian_dens_indep_matrices, &
         hamiltonian_dens_dep_nonsc
    use hf_exchange, only: HFX_STATE, hf_exchange_init, &
         hf_exchange_dkn_indep_stage
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: DKERN, kernel_create, kernel_destroy
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM, ngwf_rep_create, &
         ngwf_rep_destroy, ngwf_ham_create, ngwf_ham_destroy
    use ngwfs, only: ngwfs_merge_sets
    use projectors, only: PROJECTOR_SET
    use properties, only: properties_calculate
    use rundat, only: pub_homo_plot, pub_lumo_plot, &
         pub_homo_dens_plot, pub_lumo_dens_plot, cond_plot_joint_orbitals, &
         cond_plot_vc_orbitals, pub_output_detail, pub_aug, pub_debug_on_root, &
         pub_num_spins, pub_num_kpoints, PUB_1K, pub_kpoint_method, &
         pub_use_hfx, pub_active_region, pub_use_activehfx, &
         pub_inner_loop_iteration
    use services, only: services_flush
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_copy, sparse_embed_product, &
         sparse_embed_transpose, sparse_embed_transpose_structure, &
         sparse_embed_axpy
    use utils, only: utils_assert, utils_banner

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(FUNC_BASIS), intent(in) :: core_basis(:)
    type(FUNC_BASIS), intent(inout) :: hub_proj_basis(:)
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis(:)
    type(FUNC_BASIS), intent(inout) :: joint_ngwf_basis(:)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(:)
    type(PROJECTOR_SET), intent(inout) :: core_wvfns(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(MODEL), intent(inout)    :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    type(DKERN), intent(inout)   :: val_denskern
    type(DKERN), intent(inout)   :: cond_denskern
    type(NGWF_REP), intent(in)   :: val_rep
    type(NGWF_REP), intent(inout)   :: cond_rep
    type(NGWF_HAM), intent(inout) :: val_ham
    type(NGWF_HAM), intent(inout) :: cond_ham
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, cond_rep%nsub)
    real(kind=DP), optional, intent(in) :: kpt(3)

    ! Local Variables
    type(NGWF_REP) :: joint_rep
    type(NGWF_HAM) :: joint_ham
    type(DKERN) :: joint_denskern
    type(SPAM3_EMBED) :: cross_overlap_cj
    integer :: n_occ_cond(pub_num_spins, pub_num_kpoints)
    integer :: is,isub,jsub
    integer :: tmp_homo_plot, tmp_lumo_plot
    integer :: tmp_homo_dens_plot, tmp_lumo_dens_plot
    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecokpt
    real(kind=DP) :: loc_kpt(3)

#if 0
    ! jd: For the calculation of joint denskern. Currently not needed.
    type (SPAM3_EMBED) :: eff_cond_denskern, S_crossS_inv
    type (SPAM3_EMBED) :: S_crossKv, S_invS_crossKv, S_cross_trans
    type (SPAM3_EMBED) :: eff_cond_denskern2
#endif

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &conduction_properties'

    pub_inner_loop_iteration = -3

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine conduction_properties_calculate not ready yet for more&
         & than one k-point.')

    ! agrecokpt: currently only KP method is being implemented
    call utils_assert(pub_kpoint_method == 'KP', &
         'Subroutine conduction_ngwf_optimise currently supports&
         & only KP method for BZ sampling')

    ! agrecocmplx
    loc_cmplx = val_rep%ngwfs_on_grid(1)%iscmplx .or. &
         cond_rep%ngwfs_on_grid(1)%iscmplx

    ! agrecokpt: default is Gamma
    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt(:) = 0.0_DP
    end if

    if (pub_on_root) then
       write(stdout,'(/a)') utils_banner('-', delimiter='+')
       write(stdout,'(a)') utils_banner(' ', &
            'Starting conduction properties calculation', '|')
       write(stdout,'(a)') utils_banner('-', delimiter='+')
    end if

    ! lr408: Calculate properties for different basis sets/Hamiltonians

    ! lr408: Store original values so they can be reset after temporary overrides
    tmp_homo_plot = pub_homo_plot
    tmp_lumo_plot = pub_lumo_plot
    tmp_homo_dens_plot = pub_homo_dens_plot
    tmp_lumo_dens_plot = pub_lumo_dens_plot

    if (.not. cond_plot_vc_orbitals) then
       pub_homo_plot = -1
       pub_lumo_plot = -1
       pub_homo_dens_plot = -1
       pub_lumo_dens_plot = -1
    end if

    call properties_calculate(val_denskern, val_ham, &
         val_rep, val_ngwf_basis, proj_basis, nl_projectors, &
         hub_proj_basis, hub, core_basis, core_wvfns, mdl, hfxstate, &
         lhxc_fine, 'valence', cond_rep%n_occ)

    if (cond_plot_vc_orbitals) then

       ! never print orbitals for unprojected cond hamiltonian
       pub_homo_plot = -1
       pub_lumo_plot = -1
       pub_homo_dens_plot = -1
       pub_lumo_dens_plot = -1

       ! Unprojected Cond Properties calculation
       n_occ_cond = cond_rep%n_occ
       cond_rep%n_occ = val_rep%n_occ
       call properties_calculate(cond_denskern, cond_ham, &
            cond_rep, cond_ngwf_basis, proj_basis, nl_projectors, &
            hub_proj_basis, hub, core_basis, core_wvfns, mdl, hfxstate, &
            lhxc_fine, 'cond')
       cond_rep%n_occ = n_occ_cond

       ! print orbitals for projected cond hamiltonian if requested
       pub_lumo_plot = tmp_lumo_plot
       pub_lumo_dens_plot = tmp_lumo_dens_plot

       ! Projected Cond Properties calculation
       call properties_calculate(cond_denskern, cond_ham, &
            cond_rep, cond_ngwf_basis, proj_basis, nl_projectors, &
            hub_proj_basis, hub, core_basis, core_wvfns, mdl, hfxstate, &
            lhxc_fine, 'proj')
    end if

    ! plot joint orbitals if requested
    if (cond_plot_joint_orbitals) then
       pub_homo_plot = tmp_homo_plot
       pub_lumo_plot = tmp_lumo_plot
    else
       pub_homo_plot = -1
       pub_lumo_plot = -1
    end if

    ! lr408: Plotting orbital densities not yet implemented for joint basis
    pub_homo_dens_plot = -1
    pub_lumo_dens_plot = -1

    ! ndmh: Allocate SPAM3 structures for overlap matrix, kinetic energy,
    ! ndmh: density kernel, inverse overlap
    if (pub_on_root) write(stdout,*)
    ! agrecocmplx
    call ngwf_rep_create(joint_rep,'j',mdl,is_cmplx=loc_cmplx)

    ! agrecocmplx
    call kernel_create(joint_denskern, 'K'//joint_rep%postfix, &
            is_cmplx=loc_cmplx)

    if ((pub_on_root).and.(pub_output_detail>=NORMAL)) &
         write(stdout,'(a)') ' done'

    ! ndmh: Allocate joint rep NGWFs
    do isub=1,mdl%nsub
       call data_functions_alloc(joint_rep%ngwfs_on_grid(isub), &
            joint_ngwf_basis(isub)%n_ppds*mdl%cell%n_pts, &
            iscmplx=loc_cmplx)
    end do
    call comms_barrier

    ! lr408: Initialise conduction NGWFs to fireballs or read from file
    if (pub_on_root) write(stdout,'(a)',advance='no') 'Joint Valence + &
         &Conduction NGWF initialisation ...'

    ! ndmh: create a merged set of NGWFs containing the valence and conduction
    ! ndmh: NGWFs previously calculated
    do isub=1,mdl%nsub
       call ngwfs_merge_sets(joint_rep%ngwfs_on_grid(isub),joint_ngwf_basis(isub), &
            mdl%cell,val_rep%ngwfs_on_grid(isub),val_ngwf_basis(isub), &
            cond_rep%ngwfs_on_grid(isub), cond_ngwf_basis(isub), &
            mdl%regions(isub)%par)
    end do
    call comms_barrier
    if (pub_on_root) write(stdout,'(a/)') '... done'
    call services_flush

    ! ndmh: copy in occupation numbers from the valence rep
    joint_rep%n_occ(:,:) = val_rep%n_occ(:,:)

    ! ndmh: allocate matrices in the joint hamiltonian
    call ngwf_ham_create(joint_ham,joint_rep)

    ! ndmh: copy in dijhat from valence hamiltonian
    if (pub_aug) then
       do is=1,pub_num_spins
          call sparse_embed_copy(joint_ham%dijhat(is),val_ham%dijhat(is))
       end do
    end if

    ! jd: Initialise HFx for joint rep, if needed. Expand necessary NGWF pairs.
    do isub=1,joint_rep%nsub
       ! rc2013: only set up HFx in the active region if needed
       if (pub_use_hfx.or.(pub_use_activehfx.and.(isub==pub_active_region))) then
          call hf_exchange_init(hfxstate, joint_rep, &
               joint_ham%hfexchange(1)%m(isub,isub), mdl%cell, &
               size(mdl%regions(isub)%elements), init_rep_only = .true.) ! [**]
          ! jd: The first time around conduction NGWFs have to be expanded manually
          !     because by the time they are initialised at [*] cond-HFx has not
          !     been initialised yet (this happens at [**]).
          !     We only want the the mixed (cond-val) expansion here, so HFX_OTHER
          !     for the cond rep.
          call hf_exchange_dkn_indep_stage(hfxstate, mdl, isub, &
               mdl%regions(isub)%par, joint_rep, joint_ngwf_basis(isub), &
               rep2 = val_rep, ngwf_basis2 = val_ngwf_basis(isub), &
               basis_selector = (/1,1,2,2,REP_SWEX_HFX_OTHER/))
               !                  ^ qoh's alpha-beta-gamma-delta conv'n
               ! basis for alpha, beta: 1 (joint),
               ! basis for gamma, delta: 2 (val).
               ! swex selector: HFX_OTHER (of joint_rep)
       end if
    end do

    ! ndmh: calculate density-kernel-independent parts of the joint hamiltonian
    ! agrecokpt: at specified k-point?
    call hamiltonian_dens_indep_matrices(joint_rep, joint_ngwf_basis, &
         proj_basis, nl_projectors, hub_proj_basis, hub, mdl, &
         val_rep, val_ngwf_basis, kpt=loc_kpt)

#if 0
    ! jd: The code to build joint denskern. Seems this is not needed.
    !     Code below is based on lr_tddft_initialise() by tjz.
    do is = 1, pub_num_spins
       call sparse_embed_create(eff_cond_denskern, joint_denskern%kern%m(1,1))
       call sparse_embed_create(S_crossS_inv,joint_rep%cross_overlap, &
            joint_rep%inv_overlap)
       call sparse_embed_transpose_structure(S_cross_trans%structure, &
            joint_rep%cross_overlap)
       call sparse_embed_create(S_cross_trans, iscmplx=joint_rep%cross_overlap%iscmplx)
       call sparse_embed_create(S_crossKv, S_cross_trans, val_denskern%kern%m(1,1))
       call sparse_embed_create(S_invS_crossKv, joint_rep%inv_overlap, S_crossKv)
       call sparse_embed_create(eff_cond_denskern2, S_invS_crossKv, S_crossS_inv)

       call sparse_embed_product(S_crossS_inv, joint_rep%cross_overlap, &
            joint_rep%inv_overlap)
       call sparse_embed_transpose(S_cross_trans, joint_rep%cross_overlap)
       call sparse_embed_product(S_crossKv, S_cross_trans, val_denskern%kern%m(is,PUB_1K))
       call sparse_embed_product(S_invS_crossKv, joint_rep%inv_overlap, S_crossKv)
       call sparse_embed_product(eff_cond_denskern2, S_invS_crossKv, S_crossS_inv)
       call sparse_embed_copy(eff_cond_denskern, joint_rep%inv_overlap)
       call sparse_embed_axpy(eff_cond_denskern,eff_cond_denskern2, -1.0_DP)
       call sparse_embed_copy(joint_denskern%kern%m(is,PUB_1K), eff_cond_denskern)

       call sparse_embed_destroy(eff_cond_denskern)
       call sparse_embed_destroy(S_crossS_inv)
       call sparse_embed_destroy(S_cross_trans)
       call sparse_embed_destroy(S_crossKv)
       call sparse_embed_destroy(S_invS_crossKv)
       call sparse_embed_destroy(eff_cond_denskern2)
    enddo
#endif

    ! ndmh: calculate the rest of the joint hamiltonian
    call hamiltonian_dens_dep_nonsc(joint_ham, joint_rep, joint_ngwf_basis, &
         mdl, hfxstate, lhxc_fine, hub, val_rep, val_ham, &
         val_denskern%kern%m(:,PUB_1K), &
         val_ngwf_basis = val_ngwf_basis, &              ! jd: needed for HFx
         cond_dkn = joint_denskern%kern%m(:,PUB_1K))     ! jd: needed for HFx
         ! jd: ^ joint rep kernel passed here, but is all zero currently

    ! ndmh: calculate overlap matrix between joint NGWFs and cond NGWFs
    cross_overlap_cj%structure = 'T'//trim(adjustl(cond_rep%postfix)) // &
         trim(adjustl(joint_rep%postfix))
    ! agrecocmplx
    call sparse_embed_create(cross_overlap_cj, iscmplx=loc_cmplx)

    ! jcap: loop over regions
    do isub=1,mdl%nsub
       do jsub=1,mdl%nsub
          call function_ops_brappd_ketppd(cross_overlap_cj%m(isub,jsub), &
               cond_rep%ngwfs_on_grid(isub), cond_ngwf_basis(isub), &
               joint_rep%ngwfs_on_grid(jsub), joint_ngwf_basis(jsub), mdl%cell)
       end do
    end do

    if (pub_aug) then
       ! ndmh: Calculate the augmentation of the cross overlap matrix due to
       ! ndmh: the augmentation region part of the overlap operator
       ! jcap: with embedding, this needs to be done on a system-wide level
       call augmentation_overlap(cross_overlap_cj%p,mdl%pseudo_sp, &
            mdl%paw_sp,cond_rep%sp_overlap%p,joint_rep%sp_overlap%p)
    end if

    ! ndmh: perform a properties calculation using the joint hamiltonian
    call properties_calculate(val_denskern, joint_ham, &
         joint_rep, joint_ngwf_basis, proj_basis, nl_projectors, &
         hub_proj_basis, hub, core_basis, core_wvfns, mdl, hfxstate, &
         lhxc_fine, 'joint', cond_rep%n_occ, cross_overlap_cj, cond_denskern)

    call sparse_embed_destroy(cross_overlap_cj)

    ! ndmh: memory deallocation for joint rep/ham
    call ngwf_ham_destroy(joint_ham)
    do isub=1,mdl%nsub
       call data_functions_dealloc(joint_rep%ngwfs_on_grid(isub))
    end do
    call ngwf_rep_destroy(joint_rep)
    call kernel_destroy(joint_denskern)

    if (pub_on_root) then
       write(stdout,'(a)') utils_banner('-', delimiter='+')
       write(stdout,'(a)') utils_banner(' ', &
            'Conduction properties calculation complete', '|')
       write(stdout,'(a)') utils_banner('-', delimiter='+')
    end if

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &conduction_properties'

    return

  end subroutine conduction_properties_calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module conduction_properties
