! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                              Forces Module                                  !
!=============================================================================!
!                                                                             !
! This module calculates the total forces on each ion, after a converged      !
! total energy calculation.                                                   !
!                                                                             !
!-----------------------------------------------------------------------------!
! Module created by Nicholas Hine, 28th June 2010, based largely on existing  !
! routines from energy_and_force_mod, written by Chris-Kriton Skylaris, Arash !
! Mostofi, Peter Haynes and Nicholas Hine over the period 2000-2010.          !
!-----------------------------------------------------------------------------!

module forces

  implicit none

  private

  public :: forces_calculate
  public :: forces_apply_constraints
  public :: forces_calculate_correction

contains

  subroutine forces_calculate(total_forces,denskern,ham,lhxc_fine, &
       rep,ngwf_basis,projector_basis,nl_projectors,hub_proj_basis, &
       hub, mdl, hfxstate, ngwf_nonsc_forces, kpt, dfdtau_fine) ! agrecokpt

    !================================================================!
    ! This subroutine calculates the ionic forces                    !
    !----------------------------------------------------------------!
    ! Written by Nicholas Hine, June 2010, based on the previous     !
    ! routine internal_forces in energy_and_force_calculate, written !
    ! by Arash A. Mostofi, July 2004.                                !
    ! Modified to deal with embedding by Joseph Prentice, May 2018   !
    !================================================================!

    use augmentation, only: augmentation_density_on_grid, &
         aug_nl_calculate_forces, augmentation_density_forces
    use datatypes, only: FUNCTIONS, data_functions_alloc, data_functions_dealloc
    use bibliography, only: bibliography_cite
    use comms, only: comms_barrier, pub_on_root, pub_total_num_procs, &
         pub_my_proc_id
    use constants, only: DP, stdout, VERBOSE, max_spins,periodic_table_mass, &
         CRLF
    use classical_pot, only: classical_pot_ii_forces
    use cutoff_coulomb, only: cutoff_coulomb_ii_forces,  &
         cutoff_coulomb_locps_forces, cutoff_coulomb_nlcc_forces
    use density, only: density_on_grid
    use ewald, only: ewald_calculate_forces
    use function_basis, only: FUNC_BASIS
    ! agrecokpt
    use geometry, only: POINT
    use hamiltonian, only: hamiltonian_lhxc_calculate
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL, hubbard_calculate_forces
    use is_smeared_ions, only: smeared_ion_ii_forces, &
         smeared_ion_force_correction
    use kernel, only: DKERN
    use kernel_diis, only: kernel_diis_mu
    use lnv, only: lnv_calculate_mu
    use model_type, only: MODEL
    use multigrid_methods, only: phi_vac_for_forces
    use ngwf_gradient, only: ngwf_gradient_lnv, &
         ngwf_gradient_paw_precond_init, ngwf_gradient_paw_precond_exit
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use nonsc_forces, only: nonsc_forces_ngwfs_calc
    use openbc_locps, only: openbc_locps_calc_forces
    use parallel_strategy, only: PARAL_INFO
    use paw, only: paw_nlcc_calculate_forces, paw_nlcc_calculate_forces_recip, &
         paw_tcore_hartree_calc_forces
    use projectors, only: PROJECTOR_SET
    use pseudopotentials, only: pseudo_local_calculate_forces, &
         pseudo_nlcc_calculate_forces, pseudo_nl_calculate_forces, pseudo_get_dij
    use rundat, only: pub_any_nl_proj, pub_coulomb_cutoff, pub_dispersion, &
         pub_hubbard, pub_nlcc, pub_paw, pub_aug, pub_output_detail, &
         pub_write_forces, pub_print_qc, &
         pub_nonsc_forces, pub_aug_den_dim, pub_nhat_in_xc, pub_debug_on_root, &
         pub_zero_total_force, pub_is_smeared_ion_rep, pub_is_implicit_solvent,&
         pub_num_spins, pub_is_pbe, pub_edft, pub_kernel_diis, &
         pub_pspot_bc_is_periodic, &
         pub_num_kpoints, PUB_1K, pub_mw_total_force, &
         pub_turn_off_ewald, pub_emft, pub_active_region, &
         pub_ion_ion_bc_is_periodic, pub_is_dielectric_function
    use sparse_array, only: sparse_array_scale
    use sparse, only : SPAM3, sparse_create, sparse_copy, sparse_destroy, &
         sparse_scale
    use sparse_embed, only : SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_array_scale, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, &
         utils_assert, utils_flushed_string_output
    use vdwcorrection, only: vdwcorrection_calculate_forces

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    real(kind=DP), intent(out) :: total_forces(3,mdl%nat)
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in) :: projector_basis(:)
    type(FUNC_BASIS), intent(in) :: hub_proj_basis(:)
    ! jd: 2016.11.
    !     Temporarily relaxing rep to 'inout' to be able to update its NGWF
    !     cache in pol_emb's gradient. Nothing else in rep is modified.
    !     Will attempt to revert to intent(in) once I code up a way to defer
    !     the cache updates to a later point.
    type(NGWF_REP), intent(in) :: rep
    type(NGWF_HAM), intent(inout) :: ham
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, rep%nsub)
    type(DKERN), intent(inout) :: denskern
    ! ars: ngwf_nonsc_forces
    real(kind=DP), target, intent(inout)  :: ngwf_nonsc_forces(1:3,mdl%nat)
    ! agrecokpt
    type(POINT), optional, intent(in) :: kpt
    ! JCW: Optional gradient of XC energy per unit volume wrt kinetic energy
    ! JCW: density.
    real(kind=DP), optional, intent(in)     :: dfdtau_fine(:,:,:,:)
    ! JCW: In general, we have
    ! JCW:    dfdtau_fine(mdl%fine_grid%ld1,&
    ! JCW:    mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    ! JCW: or
    ! JCW:    dfdtau_fine(1,1,1,1)
    ! JCW: when it is unneeded.


    ! Local Variables
    integer :: ierr,atom,is,isub,jsub,iat
    real(kind=DP) :: average_force(1:3)
    ! vv
    real(kind=DP)   :: ion_mass(mdl%nat - mdl%nat_classical)
    real(kind=DP)   :: mtot
    integer         :: atom_Z
    real(kind=DP), allocatable, target, dimension(:,:) :: locps_forces
    real(kind=DP), allocatable, dimension(:,:) :: nlps_forces
    real(kind=DP), allocatable, dimension(:,:) :: nhat_forces
    real(kind=DP), allocatable, dimension(:,:) :: ewald_forces
    real(kind=DP), allocatable, dimension(:,:) :: vdw_forces
    real(kind=DP), allocatable, target, dimension(:,:) :: nlcc_forces
    real(kind=DP), allocatable, dimension(:,:) :: hub_forces
    real(kind=DP), allocatable, dimension(:,:) :: classical_forces
    real(kind=DP), allocatable, dimension(:,:) :: smeared_ion_forces
    ! ndmh: density on fine grid
    real(kind=DP), dimension(:,:,:,:), allocatable, target :: density_fine
    real(kind=DP), dimension(:,:,:,:,:), allocatable, target :: nhat_den_grad
    real(kind=DP) :: lhxc_energy   ! dummy for lhxc_calculate
    type(FUNCTIONS) :: cov_grad_on_grid(mdl%nsub)
    type(FUNCTIONS) :: contra_grad_on_grid(mdl%nsub)
    real(kind=DP) :: mu(pub_num_spins, pub_num_kpoints)
    real(kind=DP) :: f(3)
    integer :: fine_max_slabs12, fine_ld1, fine_ld2 ! jd: shortcut notation
    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecokpt
    type(POINT) :: loc_kpt
    ! jcap: temporary variables
    type(SPAM3_EMBED) :: kb_denom
    integer, allocatable :: sub_first_atom_on_proc(:,:), sub_orig_atom(:,:), &
         sub_num_atoms_on_proc(:,:), sub_glob_atom_no(:,:)
    ! jcap: subsystem forces and densities
    real(kind=DP), pointer, dimension(:,:) :: sub_locps_forces
    real(kind=DP), pointer, dimension(:,:) :: sub_nlcc_forces
    real(kind=DP), pointer, dimension(:,:) :: sub_nonsc_forces
    real(kind=DP), dimension(:,:,:,:), pointer :: sub_density_fine
    real(kind=DP), dimension(:,:,:,:), allocatable :: active_density_fine
    real(kind=DP), dimension(:,:,:,:,:), pointer :: sub_nhat_den_grad
    real(kind=DP), dimension(:,:,:,:,:), allocatable :: active_nhat_den_grad
    ! jcap: subsystem suffix
    character(len=3) :: suffix
    ! jcap: temporary arrays for use with augmentation forces
    type(SPAM3), dimension(pub_num_spins) :: aug_array, ham_array, kern_array

    ! Module parameters
    real(kind=dp), parameter :: electron_mass_si = 9.1093826e-31_dp
    real(kind=dp), parameter :: avogadro_si = 6.0221415e23_dp

    ! ------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering forces_calculate'

    call timer_clock('forces_calculate',1)
    call bibliography_cite('FORCES')

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(1)%iscmplx

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine forces_calculate not ready yet for more&
         & than one k-point.')

    ! agrecokpt
    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    ! agrecokpt: specification of kpt allowed only in complex case
    if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.&
        loc_kpt%z/=0.0_DP) then
       if (.not.loc_cmplx) then
          call utils_abort('Error in forces_calculate: &
               &non-Gamma k-point requires complex NGWFs.')
       end if
    end if


    ! Allocate arrays for each force component
    allocate(ewald_forces(3,mdl%nat),stat=ierr)
    call utils_alloc_check('forces_calculate','ewald_forces',ierr)
    allocate(locps_forces(3,mdl%nat),stat=ierr)
    call utils_alloc_check('forces_calculate','locps_forces',ierr)
    if (pub_any_nl_proj.or.pub_paw) then
       allocate(nlps_forces(3,mdl%nat),stat=ierr)
       call utils_alloc_check('forces_calculate','nlps_forces',ierr)
    end if
    if (pub_aug) then
       allocate(nhat_forces(3,mdl%nat),stat=ierr)
       call utils_alloc_check('forces_calculate','nhat_forces',ierr)
    end if
    if (pub_dispersion/='0') then
       allocate(vdw_forces(3,mdl%nat),stat=ierr)
       call utils_alloc_check('forces_calculate','vdw_forces',ierr)
    end if
    if (pub_nlcc) then
       allocate(nlcc_forces(3,mdl%nat),stat=ierr)
       call utils_alloc_check('forces_calculate','nlcc_forces',ierr)
    end if
    if (pub_hubbard) then
       allocate(hub_forces(3,mdl%nat),stat=ierr)
       call utils_alloc_check('forces_calculate','hub_forces',ierr)
    end if

    if(pub_is_smeared_ion_rep) then
       allocate(smeared_ion_forces(3,mdl%nat),stat=ierr)
       call utils_alloc_check('forces_calculate','smeared_ion_forces',ierr)
    end if

    ! kaw: Allocate array for classical forces
    if (mdl%nat_classical > 0) then
       allocate(classical_forces(3,mdl%nat),stat=ierr)
       call utils_alloc_check('forces_calculate','classical_forces',ierr)
       classical_forces = 0.0_DP
    end if

    ! jd: ion-ion forces
    ! JCW: Smeared ion representation forces
    if (pub_is_smeared_ion_rep) then
       if (.not.all(pub_ion_ion_bc_is_periodic(1:3)) ) then
          ! JCW: Forces in non-fully periodic BCs
          ! TODO Add support for mixed BC ion-ion forces (currently only
          !      fully open BCs are supported)
          call utils_assert(.not.any(pub_ion_ion_bc_is_periodic(1:3)),&
               "Error in forces_calculate: mixed BC ion-ion forces are not &
               &currently supported.")
          ! Calculate forces due to direct or Ewald contribution
          call smeared_ion_ii_forces(ewald_forces,mdl)

          ! kaw: Calculate classical forces
          if (mdl%nat_classical > 0) then
             call classical_pot_ii_forces(classical_forces,mdl%nat_classical, &
                  mdl%nat, mdl%classical_elements, mdl%elements)
          end if

       else
          ! JCW: Forces in fully periodic BCs
          if (.not.pub_turn_off_ewald) then
             ! TODO Add support for PBC smeared ion forces (currently only
             !      fully open BCs are supported)
             call utils_assert(.not.any(pub_ion_ion_bc_is_periodic(1:3)),&
                  "Error in forces_calculate: PBC smeared ion forces are not &
                  &currently supported.")
          else
             ewald_forces = 0D0
          end if
       end if
    else
       ! JCW: Forces without smeared ion representation
       if (.not.all(pub_ion_ion_bc_is_periodic(1:3)) ) then
          ! JCW: If not fully periodic...
          ! Calculate forces due to direct or Ewald contribution
          call cutoff_coulomb_ii_forces(ewald_forces,mdl)

          ! kaw: Calculate classical forces
          if (mdl%nat_classical > 0) then
             call classical_pot_ii_forces(classical_forces,mdl%nat_classical, &
                  mdl%nat, mdl%classical_elements, mdl%elements)
          end if
       else
          ! JCW: If fully periodic...
          if (.not.pub_turn_off_ewald) then
             ! Calculate forces due to Ewald contribution
             call ewald_calculate_forces(mdl%elements,mdl%classical_elements, &
                  mdl%nat,mdl%nat_classical,mdl%cell,ewald_forces,.false.)
             ! Clean up Ewald module allocatables (use locps_forces as temp dummy)
             call ewald_calculate_forces(mdl%elements,mdl%classical_elements, &
                  mdl%nat,mdl%nat_classical,mdl%cell,locps_forces,.true.)
          else
             ewald_forces = 0D0
          end if
       end if
    end if

    ! jd: Smeared ion force correction
    if(pub_is_smeared_ion_rep) then
       call smeared_ion_force_correction(phi_vac_for_forces, mdl%fine_grid, &
            mdl%elements, mdl%nat, smeared_ion_forces)
    end if

    ! jd: Smeared ion force correction
    if(pub_is_pbe /= 'NONE') then
       call utils_abort('Forces are not yet supported for Poisson-Boltzmann &
            &solvation (only for Poisson solvation)')
    end if

    if(pub_is_dielectric_function=='SOFT_SPHERE') then
       if(pub_on_root) then
          write(stdout,'(a)') 'WARNING! Forces with the soft-sphere cavity&
               &are an experimental feature. Caveat emptor!'
       end if
    end if

    ! aam: Allocate charge density slabs
    allocate(density_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
    call utils_alloc_check('forces_calculate','density_fine',ierr)
    if (mdl%nsub.gt.1) then
       allocate(sub_density_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('forces_calculate','sub_density_fine',ierr)
    else
       sub_density_fine => density_fine
    end if

    if (pub_emft) then
       allocate(active_density_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('forces_calculate','active_density_fine',ierr)
       if (pub_aug) then
          allocate(active_nhat_den_grad(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
               mdl%fine_grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim), &
               stat=ierr)
          call utils_alloc_check('forces_calculate','active_nhat_den_grad',ierr)
       end if
    end if
    if (pub_aug) then
       fine_ld1 = mdl%fine_grid%ld1
       fine_ld2 = mdl%fine_grid%ld2
       fine_max_slabs12 = mdl%fine_grid%max_slabs12
       allocate(nhat_den_grad(fine_ld1,fine_ld2,fine_max_slabs12,pub_num_spins,&
            0:pub_aug_den_dim), stat=ierr)
       call utils_alloc_check('forces_calculate','nhat_den_grad',ierr)
       if (mdl%nsub.gt.1) then
          allocate(sub_nhat_den_grad(fine_ld1,fine_ld2,fine_max_slabs12,pub_num_spins,&
               0:pub_aug_den_dim), stat=ierr)
          call utils_alloc_check('forces_calculate','sub_nhat_den_grad',ierr)
       else
          sub_nhat_den_grad => nhat_den_grad
       end if
       nhat_den_grad=0.0_DP
    end if

    ! aam: Scale density kernel to the correct number of electrons
    if (denskern%kern%num_spins == 1) then
       call sparse_embed_array_scale(denskern%kern,2.0_DP)
    end if

    ! jcap: calculate density_fine by summing over subsystems
    if (pub_aug) nhat_den_grad=0.0_DP
    density_fine=0.0_DP
    do jsub=1,mdl%nsub
       do isub=1,mdl%nsub

          call sparse_embed_extract_from_array(kern_array,&
               denskern%kern%m(:,PUB_1K),isub,jsub)

          ! aam: Calculate data-parallelised charge density
          call density_on_grid(sub_density_fine,mdl%fine_grid,mdl%dbl_grid, &
               mdl%cell,mdl%fftbox,kern_array, &
               rep%ngwf_overlap%m(isub,jsub),rep%ngwfs_on_grid(isub), &
               ngwf_basis(isub),rep%ngwfs_on_grid(jsub), ngwf_basis(jsub))

          ! ndmh: Calculate compensation density
          if (pub_aug) then
             sub_nhat_den_grad = 0.0_DP
             call augmentation_density_on_grid(sub_nhat_den_grad, &
                  mdl%fine_grid,mdl%cell,mdl%regions(isub)%pseudo_sp,&
                  mdl%regions(isub)%paw_sp,mdl%aug_box, &
                  kern_array,rep%sp_overlap%m(isub,isub))
             sub_density_fine = sub_density_fine + sub_nhat_den_grad(:,:,:,:,0)
          end if

          call sparse_embed_destroy_extracted_array(kern_array)

          if (mdl%nsub.gt.1) density_fine = density_fine + sub_density_fine
          if (pub_emft.and.(isub==jsub).and.(isub==pub_active_region)) &
               active_density_fine = sub_density_fine
          if (pub_aug.and.(mdl%nsub.gt.1)) nhat_den_grad = nhat_den_grad + &
               sub_nhat_den_grad
          if (pub_aug.and.pub_emft.and.(isub==jsub).and.(isub==pub_active_region)) &
               active_nhat_den_grad = sub_nhat_den_grad
       end do ! isub
    end do ! jsub

    if (mdl%nsub.gt.1) then
       deallocate(sub_density_fine,stat=ierr)
       call utils_dealloc_check('forces_calculate','sub_density_fine',ierr)
       if (pub_aug) then
          deallocate(sub_nhat_den_grad,stat=ierr)
          call utils_dealloc_check('forces_calculate','sub_nhat_den_grad',ierr)
       end if
    end if

    do isub=1,mdl%nsub
       ! jcap: allocate appropriate arrays
       if (mdl%nsub.gt.1) then
          allocate(sub_locps_forces(3,mdl%regions(isub)%par%nat),stat=ierr)
          call utils_alloc_check('forces_calculate','sub_locps_forces',ierr)
          if (pub_nlcc) then
             allocate(sub_nlcc_forces(3,mdl%regions(isub)%par%nat),stat=ierr)
             call utils_alloc_check('forces_calculate','sub_nlcc_forces',ierr)
          end if
       else
          ! jcap: just point to full arrays if mdl%nsub=1
          sub_locps_forces => locps_forces
          if (pub_nlcc) sub_nlcc_forces => nlcc_forces
       end if


       if (pub_coulomb_cutoff) then

          ! Calculate forces due to local part of ionic pseudopotential
          call cutoff_coulomb_locps_forces(density_fine,mdl,sub_locps_forces,isub,&
               mdl%regions(isub)%par)

          ! Calculate forces due to the NLCC core charge
          if (pub_nlcc) then
             if (pub_aug) then
                if (pub_emft) then
                   call cutoff_coulomb_nlcc_forces(density_fine, &
                        mdl%core_density_fine,mdl,sub_nlcc_forces,nhat_den_grad,&
                        ireg=isub,par=mdl%regions(isub)%par,&
                        active_density_fine=active_density_fine,&
                        active_nhat_den_grad_fine=active_nhat_den_grad)
                else
                   call cutoff_coulomb_nlcc_forces(density_fine, &
                        mdl%core_density_fine,mdl,sub_nlcc_forces,nhat_den_grad,&
                        ireg=isub,par=mdl%regions(isub)%par)
                end if
             else
                if (pub_emft) then
                   call cutoff_coulomb_nlcc_forces(density_fine,&
                        mdl%core_density_fine,mdl,sub_nlcc_forces,ireg=isub,&
                        par=mdl%regions(isub)%par,&
                        active_density_fine=active_density_fine)
                else
                   call cutoff_coulomb_nlcc_forces(density_fine,&
                        mdl%core_density_fine,mdl,sub_nlcc_forces,ireg=isub,&
                        par=mdl%regions(isub)%par)
                end if
             end if
          end if

       else

          ! Calculate forces due to local part of ionic pseudopotential
          ! jd: Handle open BC local pseudo forces separately
          if (.not.any(pub_pspot_bc_is_periodic)) then
             ! JCW: If fully open, use open BC local pseudopotential routine
             call openbc_locps_calc_forces(density_fine, &
                  mdl%fine_grid, mdl%cell,mdl%regions(isub)%elements,&
                  mdl%regions(isub)%pseudo_sp, mdl%regions(isub)%paw_sp,&
                  sub_locps_forces,mdl%regions(isub)%par)
          else
             ! JCW: Periodic or mixed periodic/open BCs
             if (.not.pub_paw) then
                call pseudo_local_calculate_forces(density_fine,mdl%fine_grid, &
                     mdl%cell,mdl%regions(isub)%elements,mdl%regions(isub)%pseudo_sp,&
                     sub_locps_forces,mdl%regions(isub)%par)
             else
                call paw_tcore_hartree_calc_forces(density_fine,mdl%fine_grid, &
                     mdl%regions(isub)%elements,mdl%regions(isub)%paw_sp,&
                     sub_locps_forces,mdl%regions(isub)%par)
             end if
          end if

          if (pub_nlcc) then
             ! Calculate forces due to the NLCC core charge
             if (.not.pub_paw) then
                if (pub_emft) then
                   call pseudo_nlcc_calculate_forces(density_fine,&
                        mdl%core_density_fine,mdl%fine_grid,mdl%cell, &
                        mdl%regions(isub)%elements,mdl%regions(isub)%pseudo_sp,&
                        sub_nlcc_forces,mdl%regions(isub)%par,&
                        active_density_fine=active_density_fine)
                else
                   call pseudo_nlcc_calculate_forces(density_fine,&
                        mdl%core_density_fine,mdl%fine_grid,mdl%cell, &
                        mdl%regions(isub)%elements,mdl%regions(isub)%pseudo_sp,&
                        sub_nlcc_forces,mdl%regions(isub)%par)
                end if
             else
                if (pub_emft) then
                   call paw_nlcc_calculate_forces(density_fine,&
                        mdl%core_density_fine,nhat_den_grad,mdl%fine_grid, &
                        mdl%cell,mdl%regions(isub)%elements,mdl%regions(isub)%paw_sp,&
                        sub_nlcc_forces,mdl%regions(isub)%par,&
                        active_density=active_density_fine,&
                        active_nhat_den_grad=active_nhat_den_grad)
                else
                   call paw_nlcc_calculate_forces(density_fine,&
                        mdl%core_density_fine,nhat_den_grad,mdl%fine_grid, &
                        mdl%cell,mdl%regions(isub)%elements,mdl%regions(isub)%paw_sp,&
                        sub_nlcc_forces,mdl%regions(isub)%par)
                end if
#ifdef TESTCOREDEN
                call internal_print_forces(nlcc_forces, &
                     'Non-Recip NLCC forces','Unconstrained')
                if (pub_emft) then
                   call paw_nlcc_calculate_forces_recip(density_fine,&
                        mdl%core_density_fine,nhat_den_grad,mdl%fine_grid, &
                        mdl%cell,mdl%regions(isub)%elements,mdl%regions(isub)%paw_sp,&
                        sub_nlcc_forces,mdl%regions(isub)%par,&
                        active_density=active_density_fine,&
                        active_nhat_den_grad=active_nhat_den_grad)
                else
                   call paw_nlcc_calculate_forces_recip(density_fine,&
                        mdl%core_density_fine,nhat_den_grad,mdl%fine_grid, &
                        mdl%cell,mdl%regions(isub)%elements,mdl%regions(isub)%paw_sp,&
                        sub_nlcc_forces,mdl%regions(isub)%par)
                end if
#endif
             end if
          end if

       end if
       ! jcap: if mdl%nsub=1, the sub_ arrays already point to the
       ! total arrays, so everything is already in the correct place
       if (mdl%nsub.gt.1) then
          do iat=1,mdl%regions(isub)%par%nat
             locps_forces(:,mdl%regions(isub)%elements(iat)%global_atom_number)=&
                  sub_locps_forces(:,iat)
             if (pub_nlcc) nlcc_forces(:,mdl%regions(isub)%elements(iat)%global_atom_number)=&
                  sub_nlcc_forces(:,iat)
          end do
          ! jcap: deallocate appropriate arrays
          deallocate(sub_locps_forces,stat=ierr)
          call utils_dealloc_check('forces_calculate','sub_locps_forces',ierr)
          if (pub_nlcc) then
             deallocate(sub_nlcc_forces,stat=ierr)
             call utils_dealloc_check('forces_calculate','sub_nlcc_forces',ierr)
          end if
       end if

    end do ! isub

    ! ndmh: Calculate nonlocal pseudopotential forces if necessary
    if ((.not.pub_aug).and.pub_any_nl_proj) then

       ! agrecokpt: at specified k-point
       ! jcap: nasty workaround here to get the quantities we need
       ! without having to pass mdl or regions
       allocate(sub_first_atom_on_proc(mdl%nsub,0:pub_total_num_procs-1),stat=ierr)
       call utils_alloc_check('forces_calculate','sub_first_atom_on_proc',ierr)
       allocate(sub_num_atoms_on_proc(mdl%nsub,0:pub_total_num_procs-1),stat=ierr)
       call utils_alloc_check('forces_calculate','sub_num_atoms_on_proc',ierr)
       allocate(sub_orig_atom(mdl%nsub,mdl%nat),stat=ierr)
       call utils_alloc_check('forces_calculate','sub_orig_atom',ierr)
       allocate(sub_glob_atom_no(mdl%nsub,mdl%nat),stat=ierr)
       call utils_alloc_check('forces_calculate','sub_glob_atom_no',ierr)
       sub_first_atom_on_proc=0
       sub_num_atoms_on_proc=0
       sub_orig_atom=0
       sub_glob_atom_no=0
       ! jcap: precalculate kb_denom too
       kb_denom%structure = 'E'
       call sparse_embed_create(kb_denom,iscmplx=rep%ngwfs_on_grid(1)%iscmplx)
       do isub=1,mdl%nsub
          call pseudo_get_dij(kb_denom%m(isub,isub),mdl%regions(isub)%pseudo_sp)
          ! rc2013: avoid array issues by doing this across each proc separately
          sub_first_atom_on_proc(isub,pub_my_proc_id) = &
               mdl%regions(isub)%par%first_atom_on_proc(pub_my_proc_id)
          sub_num_atoms_on_proc(isub,pub_my_proc_id) = &
               mdl%regions(isub)%par%num_atoms_on_proc(pub_my_proc_id)
          sub_orig_atom(isub,1:mdl%regions(isub)%par%nat)=&
               mdl%regions(isub)%par%orig_atom
          do iat=1,mdl%regions(isub)%par%nat
             sub_glob_atom_no(isub,iat)=&
                  mdl%regions(isub)%elements(iat)%global_atom_number
          end do
       end do
       call pseudo_nl_calculate_forces(nlps_forces, &
            rep%sp_overlap,mdl%pseudo_sp,rep%ngwfs_on_grid,ngwf_basis, &
            projector_basis,nl_projectors,mdl%cell,mdl%fftbox, &
            denskern%kern%m(:,PUB_1K),mdl%nsub,mdl%nat,&
            sub_first_atom_on_proc,sub_num_atoms_on_proc,sub_orig_atom,&
            sub_glob_atom_no,kb_denom, kpt=loc_kpt)
       ! kb_denom is destroyed within pseudo_nl_calculate_forces
       deallocate(sub_first_atom_on_proc,stat=ierr)
       call utils_dealloc_check('forces_calculate','sub_first_atom_on_proc',ierr)
       deallocate(sub_num_atoms_on_proc,stat=ierr)
       call utils_dealloc_check('forces_calculate','sub_num_atoms_on_proc',ierr)
       deallocate(sub_orig_atom,stat=ierr)
       call utils_dealloc_check('forces_calculate','sub_orig_atom',ierr)
       deallocate(sub_glob_atom_no,stat=ierr)
       call utils_dealloc_check('forces_calculate','sub_glob_atom_no',ierr)

    ! ndmh: Calculate PAW force terms
    else if (pub_aug) then

       ! ndmh: recalculate LHXC and ham%dijhat
       call hamiltonian_lhxc_calculate(lhxc_fine,lhxc_energy,ham%dijhat, &
            mdl,rep%ngwfs_on_grid,ngwf_basis,denskern%kern,rep%ngwf_overlap, &
            rep%sp_overlap,add_xc_pot=pub_nhat_in_xc)

       ! Calculate PAW nonlocal force terms
       ! agrecokpt: at specified k-point
       ! jcap: need temporary arrays to pass in
       call sparse_embed_extract_from_array(ham_array,ham%ham)
       call sparse_embed_extract_from_array(aug_array,ham%dijhat)
       call sparse_embed_extract_from_array(kern_array,&
            denskern%kern%m(:,PUB_1K))

       call aug_nl_calculate_forces(nlps_forces, &
            rep%ngwfs_on_grid(1),ngwf_basis(1),projector_basis(1),&
            nl_projectors(1),rep%sp_overlap%p,mdl,&
            rep%inv_overlap%p,kern_array, &
            ham_array,aug_array,kpt=loc_kpt)

       ! Calculate PAW compensation density forces
       ! rc2013: EMBED_FIX: this enforces 1 lhxc potential
       call augmentation_density_forces(nhat_forces, kern_array,&
            rep%sp_overlap%p,lhxc_fine(:,:,:,:,1),mdl)

       call sparse_embed_destroy_extracted_array(ham_array)
       call sparse_embed_destroy_extracted_array(aug_array)
       call sparse_embed_destroy_extracted_array(kern_array)

    end if

    ! qoh: Calculate vdw_forces if necessary
    if (pub_dispersion /= '0') then
       call vdwcorrection_calculate_forces(vdw_forces,mdl%elements,mdl%cell,&
            mdl%par)
    end if

    ! ddor: Calculate DFT+U forces if necessary
    if (pub_hubbard) then
       call sparse_embed_extract_from_array(kern_array,denskern%kern%m(:,PUB_1K))
       call hubbard_calculate_forces(hub_forces, &
            rep%ngwfs_on_grid(1),ngwf_basis(1),hub_proj_basis(1),hub, &
            mdl%cell,mdl%fftbox,kern_array, &
            rep%hub_overlap%p, rep%hub_overlap_t%p)
       call sparse_embed_destroy_extracted_array(kern_array)
    endif

    ! aam: Scale density kernel back to half Ne
    if (denskern%kern%num_spins == 1) then
       call sparse_embed_array_scale(denskern%kern, 0.5_DP)
    end if

    !==============================================================!
    !*** Uncomment in order to convert forces from Eh/a to eV/A ***!
    !     ewald_forces =      ewald_forces*HARTREE_IN_EVS*ANGSTROM !
    !     locps_forces =      locps_forces*HARTREE_IN_EVS*ANGSTROM !
    !      nlps_forces =       nlps_forces*HARTREE_IN_EVS*ANGSTROM !
    !      nhat_forces =       nhat_forces*HARTREE_IN_EVS*ANGSTROM !
    !       vdw_forces =        vdw_forces*HARTREE_IN_EVS*ANGSTROM !
    !      nlcc_forces =       nlcc_forces*HARTREE_IN_EVS*ANGSTROM !
    !      nhat_forces =       nhat_forces*HARTREE_IN_EVS*ANGSTROM !
    !ngwf_nonsc_forces = ngwf_nonsc_forces*HARTREE_IN_EVS*ANGSTROM !
    !==============================================================!

    ! Sum all of the force contributions
    total_forces = 0.0_DP

    !kaw: Add classical forces when relevant
    if ((mdl%nat_classical > 0).and.(pub_coulomb_cutoff)) then
       ewald_forces = ewald_forces + classical_forces
    end if
    total_forces = total_forces + ewald_forces + locps_forces

    ! jd: Smeared ion force correction
    if(pub_is_smeared_ion_rep) then
       total_forces = total_forces + smeared_ion_forces
    end if

    if (pub_any_nl_proj.or.pub_paw) then
       total_forces = total_forces + nlps_forces
    end if
    if (pub_aug) then
       total_forces = total_forces + nhat_forces
    end if
    if (pub_dispersion/='0') then
       total_forces = total_forces + vdw_forces
    end if
    if (pub_nlcc) then
       total_forces = total_forces + nlcc_forces
    end if
    if (pub_nonsc_forces) then
       if (all(ngwf_nonsc_forces==-999_DP)) then
          if(pub_output_detail > VERBOSE) then
             call utils_flushed_string_output('NGWF non self-consistent forces &
                  &have not been calculated yet and were requested.'//CRLF//&
                  &'They will be calculated now.'//CRLF)
          end if
          ! Allocate workspace
          do isub=1,mdl%nsub
             call data_functions_alloc(contra_grad_on_grid(isub), &
                  ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
             call data_functions_alloc(cov_grad_on_grid(isub), &
                  ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
          end do

          ! Calculate mu term for gradient
          if ((.not.pub_kernel_diis).and.(.not.pub_edft)) then
             call lnv_calculate_mu(mu, denskern, ham, rep, ngwf_basis)
          else if (pub_kernel_diis) then
             call kernel_diis_mu(mu(:,PUB_1K), denskern%kern, &
                  ham%ham, rep%overlap)
          else
             mu(:,:) = 0.0_DP
          end if

          do isub=1,mdl%nsub
             ! jcap: allocate temporary array
             if (mdl%nsub.gt.1) then
                allocate(sub_nonsc_forces(1:3,mdl%regions(isub)%par%nat),stat=ierr)
                call utils_alloc_check('forces_calculate','sub_nonsc_forces', &
                     ierr)
             else
                sub_nonsc_forces => ngwf_nonsc_forces
             end if
             sub_nonsc_forces=0.0_DP

             ! gcc32: initialize preconditioning for PAW
             call ngwf_gradient_paw_precond_init(nl_projectors(isub), &
                  mdl%regions(isub)%paw_sp,mdl%cell,mdl%fftbox)

             ! Calculate Gradient
             call ngwf_gradient_lnv(contra_grad_on_grid, cov_grad_on_grid, &
                  denskern, rep, ngwf_basis, projector_basis, hub_proj_basis, &
                  nl_projectors, lhxc_fine, ham, hub, mu(:,PUB_1K), mdl, &
                  hfxstate, isub, dfdtau_fine = dfdtau_fine)

             ! Calculate NonSC Forces
             if((pub_output_detail > VERBOSE).and.(mdl%nsub.gt.1)) then
                write(suffix,'(i0)')isub
                call utils_flushed_string_output('Region '//trim(suffix)//':')
             end if
             ! jcap: subregion suffix
             ! rc2013: only if we actually need it
             suffix = trim(rep%postfix)
             if(mdl%nsub .gt. 1) &
                  write(suffix,'(a,i0,i0)') trim(rep%postfix),isub,isub
             call nonsc_forces_ngwfs_calc(sub_nonsc_forces,&
                  rep%ngwfs_on_grid(isub), contra_grad_on_grid(isub), ngwf_basis(isub), &
                  mdl%cell, mdl%fftbox, mdl%regions(isub)%par, suffix=trim(suffix))

             ! jcap: copy data into full array
             ! jcap: unless we only have one region, in which case the
             ! pointer means it is already in the full array
             if (mdl%nsub.gt.1) then
                do iat=1,mdl%regions(isub)%par%nat
                   ngwf_nonsc_forces(:,mdl%regions(isub)%elements(iat)%global_atom_number)&
                        = sub_nonsc_forces(:,iat)
                end do
                deallocate(sub_nonsc_forces,stat=ierr)
                call utils_dealloc_check('forces_calculate','sub_nonsc_forces', &
                     ierr)
             end if

             ! gcc32: exit preconditioning for PAW
             call ngwf_gradient_paw_precond_exit(nl_projectors(isub))

          end do ! isub
          ! Deallocate workspace
          do isub=1,mdl%nsub
             call data_functions_dealloc(cov_grad_on_grid(isub))
             call data_functions_dealloc(contra_grad_on_grid(isub))
          end do

       end if

       total_forces = total_forces + ngwf_nonsc_forces

    end if
    if (pub_hubbard) then
       total_forces = total_forces + hub_forces
    end if

    ! Write forces (root proc only)
    if (pub_on_root) then
       if (pub_write_forces) then
          if (pub_output_detail > VERBOSE) then
             if (.not.all(pub_ion_ion_bc_is_periodic(1:3)) ) then
                call internal_print_forces(ewald_forces, &
                     'Ion-Ion forces','Unconstrained')
             else
                call internal_print_forces(ewald_forces, &
                     'Ewald forces','Unconstrained')
             end if
             call internal_print_forces(locps_forces, &
                  'Local potential forces','Unconstrained')
             if (pub_is_smeared_ion_rep) then
                if (pub_is_implicit_solvent) then
                   call internal_print_forces(smeared_ion_forces, &
                        'Smeared ion correction and solvation forces', &
                        'Unconstrained')
                else
                   call internal_print_forces(smeared_ion_forces, &
                        'Smeared ion correction forces','Unconstrained')
                end if
             end if
             if (pub_aug) then
                call internal_print_forces(nhat_forces, &
                     'Compensation density forces','Unconstrained')
             end if
             if (pub_paw.or.pub_any_nl_proj) then
                call internal_print_forces(nlps_forces, &
                     'Non-local potential forces','Unconstrained')
             end if
             if(pub_nonsc_forces) then
                call internal_print_forces(ngwf_nonsc_forces, &
                     'NGWF non self-consistent forces','Unconstrained')
             end if
             if (pub_dispersion /='0') then
                call internal_print_forces(vdw_forces, &
                     'Dispersion forces','Unconstrained')
             end if
             if (pub_nlcc) then
                call internal_print_forces(nlcc_forces,'NLCC forces', &
                     'Unconstrained')
             end if
             if (pub_hubbard) then
                call internal_print_forces(hub_forces,'DFT+U forces', &
                     'Unconstrained')
             end if
           end if
       end if
    end if

    ! rc2013: get the correction to the total force
    call forces_calculate_correction(total_forces,mdl)

    ! Deallocate
    if(pub_is_smeared_ion_rep) then
       deallocate(smeared_ion_forces,stat=ierr)
       call utils_dealloc_check('forces_calculate','smeared_ion_forces',ierr)
    end if
    if (pub_aug) then
       deallocate(nhat_den_grad,stat=ierr)
       call utils_dealloc_check('forces_calculate','nhat_den_grad',ierr)
       if (pub_emft) then
          deallocate(active_nhat_den_grad,stat=ierr)
          call utils_dealloc_check('forces_calculate','active_nhat_den_grad',ierr)
       end if
    end if
    deallocate(density_fine,stat=ierr)
    call utils_dealloc_check('forces_calculate','density_fine',ierr)
    if (pub_emft) then
       deallocate(active_density_fine,stat=ierr)
       call utils_dealloc_check('forces_calculate','active_density_fine',ierr)
    end if
    if (pub_hubbard) then
       deallocate(hub_forces,stat=ierr)
       call utils_dealloc_check('forces_calculate','hub_forces',ierr)
    end if
    if (pub_nlcc) then
       deallocate(nlcc_forces,stat=ierr)
       call utils_dealloc_check('forces_calculate','nlcc_forces',ierr)
    end if
    if (pub_dispersion/='0') then
       deallocate(vdw_forces,stat=ierr)
       call utils_dealloc_check('forces_calculate','vdw_forces',ierr)
    end if
    if (pub_aug) then
       deallocate(nhat_forces,stat=ierr)
       call utils_dealloc_check('forces_calculate','nhat_forces',ierr)
    end if
    if (pub_paw.or.pub_any_nl_proj) then
       deallocate(nlps_forces,stat=ierr)
       call utils_dealloc_check('forces_calculate','nlps_forces',ierr)
    end if
    deallocate(locps_forces,stat=ierr)
    call utils_dealloc_check('forces_calculate','locps_forces',ierr)
    deallocate(ewald_forces,stat=ierr)
    call utils_dealloc_check('forces_calculate','ewald_forces',ierr)

    if (pub_print_qc) call forces_qc_output(total_forces,mdl%nat)

    call timer_clock('forces_calculate',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving forces_calculate'

    return

contains

    subroutine internal_print_forces(forces,component_name,title)

      use utils, only: utils_banner

      ! Arguments
      real(kind=DP), intent(in) :: forces(:,:)
      character(*), intent(in) :: component_name, title

      ! Local variables
      integer :: bef, aft, w, l1, l2
      w = 58
      l1 = len(title)
      l2 = len(component_name)

      if (pub_on_root) then
         write(stdout,'(a)') ' '
         bef = (w-l1-2)/2
         aft = (w-l1-2)/2 + modulo(w-l1-2,2)
         write(stdout,'(a)') repeat('*',bef)//' '//title//' '//repeat('*',aft)
         bef = (w-l2-2)/2
         aft = (w-l2-2)/2 + modulo(w-l2-2,2)
         write(stdout,'(a)') repeat('*',bef)//' '//component_name// &
              ' '//repeat('*',aft)
         write(stdout,'(a)') '*'//repeat(' ',w-2)//'*'
         write(stdout,'(a)') '* Element  Atom         &
              &Cartesian components (Eh/a)      *'
         write(stdout,'(a)') '* '//repeat('-',w-4)//' *'
         write(stdout,'(a)') '*                       x            &
              &y            z      *'
         write(stdout,'(a)') '*'//repeat(' ',w-2)//'*'
         do atom=1,mdl%nat
            write(stdout,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
                 mdl%elements(atom)%symbol,' ',&
                 atom,'   ',forces(:,atom),' *'
         end do
         write(stdout,'(a)') '*'//repeat(' ',w-2)//'*'
         write(stdout,'(a)') repeat('*',w)
         write(stdout,'(a,3f13.8)') 'TOTAL:           ',&
              sum(forces(1,1:mdl%nat)),  &
              sum(forces(2,1:mdl%nat)),  &
              sum(forces(3,1:mdl%nat))
      end if

    end subroutine internal_print_forces

  end subroutine forces_calculate

  subroutine forces_apply_constraints(forces,mdl)
    !=========================================================================!
    ! Apply ionic constraints to the ionic forces                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   forces, intent=inout, the ionic forces                                !
    !   elements, intent=in, the element data array                           !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi, 02/08/2005                                  !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use model_type, only: MODEL

    implicit none

    ! Arguments
    type(MODEL), intent(in)      :: mdl
    real(kind=dp), intent(inout) :: forces(3,1,mdl%nat)

    ! Local Variables
    integer :: iat
    real(kind=dp) :: proj,norm
    real(kind=dp) :: f(1:3),v(1:3)

    do iat=1,mdl%nat

       f(1:3) = forces(1:3,1,iat)
       v(1:3) = mdl%elements(iat)%ion_constraint(1:3)

       select case (mdl%elements(iat)%ion_constraint_type)
       case ('NONE') ; continue                    ! aam: no constraint
       case ('PLANE')                              ! aam: constrained perp. to v
          proj = f(1)*v(1) + f(2)*v(2) + f(3)*v(3)           ! aam: f.v
          norm = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)           ! aam: v.v
          forces(1:3,1,iat) = f(1:3) - (proj/norm)*v(1:3)    ! aam: f-(f.v/v.v)v
       case ('LINE')                               ! aam: constrained para. to v
          proj = f(1)*v(1) + f(2)*v(2) + f(3)*v(3)           ! aam: f.v
          norm = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)           ! aam: v.v
          forces(1:3,1,iat) = (proj/norm)*v(1:3)             ! aam: (f.v/v.v)v
       case ('FIXED') ; forces(1:3,1,iat) = 0.0_dp ! aam: fixed
       case default
          if (pub_on_root) then
             write(stdout,'(a,i6,a)') 'Error in forces_apply_constraints: &
                  &illegal value for mdl%elements(',iat,')%ion_constraint_type'
          endif
       end select

    enddo

  end subroutine forces_apply_constraints

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine forces_qc_output(total_forces,nat)

    !=========================================================================!
    ! Prints out quality-control lines <QC> for qc-testing the forces.        !
    !-------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano, 19/11/2010.                             !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    use comms, only: pub_on_root, comms_barrier
    use constants, only: DP
    use model_type, only: MODEL
    use utils, only: utils_qc_print

    implicit none

    ! Arguments
    integer, intent(in) :: nat  ! jcap
    real(kind=DP), intent(in) ::  total_forces(1:3, 1:nat)

    ! Local Variables
    integer :: icomp, iatom
    character(len=64) :: icomp_iatom_string

    if(pub_on_root) then

       do icomp = 1, 3
          do iatom = 1, nat
             write(icomp_iatom_string,'(a1,i0,a1,i0,a1)') '(', icomp , ',', &
                  iatom, ')'
             call utils_qc_print('total_force'//trim(icomp_iatom_string), &
                  total_forces(icomp,iatom))
          end do
       end do

    end if

    call comms_barrier

  end subroutine forces_qc_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine forces_calculate_correction(total_forces,mdl)

    !=========================================================================!
    ! Calculates, applies and prints correction to give zero total force      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   forces, intent=inout, the ionic forces                                !
    !   mdl, intent=in, the model                                             !
    !-------------------------------------------------------------------------!
    ! Adapted from existing code for embedding by Joseph Prentice, May 2018   !
    !=========================================================================!

    use constants, only : DP,stdout,periodic_table_mass,VERBOSE
    use rundat, only : pub_zero_total_force,pub_mw_total_force,&
         pub_write_forces,pub_output_detail, &
         pub_task, pub_geom_output_detail
    use model_type, only : MODEL
    use comms, only : comms_barrier, pub_on_root

    implicit none

    real(DP), intent(inout) :: total_forces(:,:)
    type(MODEL), intent(in) :: mdl

    ! Local variables
    real(DP) :: average_force(3)
    real(DP) :: ion_mass(mdl%nat-mdl%nat_classical)
    real(DP) :: mtot,f(3)
    integer :: atom,atom_Z

    ! Module parameters
    real(kind=dp), parameter :: electron_mass_si = 9.1093826e-31_dp
    real(kind=dp), parameter :: avogadro_si = 6.0221415e23_dp

    ! Calculate the average force
    ! kaw: Account for classical atoms
    average_force = 0.0_dp
    if(pub_zero_total_force) then ! jd: this is on by default
       do atom=1,mdl%nat - mdl%nat_classical
          average_force(:) = average_force(:) + total_forces(:,atom)
       end do
       average_force(:) = average_force(:) / (mdl%nat - mdl%nat_classical)

       ! Remove the average force
       do atom=1,mdl%nat - mdl%nat_classical
          total_forces(:,atom) = total_forces(:,atom) - average_force(:)
       end do
    end if

    ! Calculate weighted-correction
    mtot = 0.0_DP
    if (pub_mw_total_force) then
       do atom=1,mdl%nat - mdl%nat_classical
          average_force(:) = average_force(:) + total_forces(:,atom)
          atom_Z = mdl%elements(atom)%atomic_number
          ion_mass(atom) = periodic_table_mass(atom_Z)*1e-3_dp/avogadro_si/electron_mass_si
          mtot = mtot + ion_mass(atom)
       end do

       ! apply correction
       do atom=1,mdl%nat - mdl%nat_classical
          total_forces(:,atom) = total_forces(:,atom) - (ion_mass(atom)/mtot)*average_force(:)
       end do
    end if

    ! lk: write forces during a geometry optimization only if the
    !     user has opted for an ouput detail higher than VERBOSE.
    if (pub_on_root .and. pub_write_forces &
       .and. pub_geom_output_detail > VERBOSE) then
       write(stdout,'(a)') ' '
       write(stdout,'(a)') '********************* Unconstrained &
            &**********************'
       write(stdout,'(a)') '****** Correction to ensure the&
            & total force is zero ******'
       write(stdout,'(a)') '*                                  &
            &                      *'
       write(stdout,'(a)') '* Element  Atom         &
            &Cartesian components (Eh/a)      *'
       write(stdout,'(a)') '* ----------------------------------&
            &-------------------- *'
       write(stdout,'(a)') '*                       x            &
            &y            z      *'
       write(stdout,'(a)') '*                                 &
            &                       *'
       do atom=1,mdl%nat
          if(pub_zero_total_force) then
             f = -average_force(:)
          elseif(pub_mw_total_force) then
             f = - (ion_mass(atom)/mtot)*average_force(:)
          else
             f(:) = 0D0
          end if
          write(stdout,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
               mdl%elements(atom)%symbol,' ',atom,&
               '   ',f,' *'
       end do
       write(stdout,'(a)') '*                                  &
            &                      *'
       write(stdout,'(a)') '***********************************&
            &***********************'
       write(stdout,'(a,3f13.8)') 'TOTAL:           ',&
            -average_force(1) * (mdl%nat - mdl%nat_classical),  &
            -average_force(2) * (mdl%nat - mdl%nat_classical),  &
            -average_force(3) * (mdl%nat - mdl%nat_classical)
       ! jd: NB we output the negative, as this term has been subtracted
       !     rather than added to total force. In this way all components
       !     listed in the output will *add* up to the final force.
       write(stdout,'(/a)') '********************* Unconstrained &
            &**********************'
       write(stdout,'(a)') '************************* Forces &
            &*************************'
       write(stdout,'(a)') '*                                &
            &                        *'
       write(stdout,'(a)') '* Element  Atom         Cartesian &
            &components (Eh/a)      *'
       write(stdout,'(a)') '* -------------------------------&
            &----------------------- *'
       write(stdout,'(a)') '*                       x            &
            &y            z      *'
       write(stdout,'(a)') '*                                 &
            &                       *'
       do atom=1,mdl%nat
          write(stdout,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
               mdl%elements(atom)%symbol,' ',atom,'   ', &
               total_forces(:,atom),' *'
       end do
       write(stdout,'(a)') '*                                 &
            &                       *'
       write(stdout,'(a)') '**********************************&
            &************************'
       write(stdout,'(a,3f13.8,a)') '* TOTAL:         ',&
            sum(total_forces(1,1:mdl%nat)), &
            sum(total_forces(2,1:mdl%nat)), &
            sum(total_forces(3,1:mdl%nat)),' *'
       write(stdout,'(a)') '**********************************&
            &************************'
    end if
    call comms_barrier

  end subroutine forces_calculate_correction

end module forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
