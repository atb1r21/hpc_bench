! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!==============================================================================!
!                                                                              !
!         Linear-Response Phonons module                                       !
!                                                                              !
!------------------------------------------------------------------------------!
! This module was created by Gabriel Constantinescu in 2014-2017               !
!==============================================================================!

module lr_phonons

  use constants, only: DP

  implicit none


  private

  ! Conversion constants
  real(kind=DP), parameter :: electron_mass_u =5.4857990943E-4_DP
  real(kind=DP), parameter :: au2THz = 1.0_DP/2.41888468E-5_DP
  real(kind=DP), parameter :: au2inv_cm = au2THz/2.99792458E-2_DP

  real(kind=DP), allocatable, dimension(:,:,:) :: precond_func_recip

  public :: lr_phonons_calculate

contains


  subroutine lr_phonons_calculate(denskern, val_rep, val_basis, &
       hub_proj_basis, proj_basis, nl_projectors, mdl, hfxstate, lhxc_fine, hub)

    ! ========================================================================!
    ! Main subroutine, calculates the linear-response of the system under     !
    ! atomic and electric field perturbations                                 !
    ! ------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in 2014-2017                          !
    ! Modified for embedding by Joseph Prentice, September 2018               !
    ! ========================================================================!


    use augmentation, only: aug_projector_denskern, &
         augmentation_density_on_grid, augmentation_overlap
    use datatypes, only: data_functions_alloc, FUNCTIONS, COEF, &
         data_functions_dealloc, data_set_to_zero
    use comms, only: pub_on_root, comms_barrier, pub_my_proc_id, &
         pub_total_num_procs, pub_root_proc_id, comms_bcast
    use constants, only: DP, stdout, NORMAL, cmplx_1, cmplx_0, PI, &
         periodic_table_mass, TWO_PI
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS, function_basis_deallocate
    use function_ops, only: function_ops_brappd_ketppd
    use geometry, only: POINT, operator(*), operator(+), operator(.DOT.), &
         operator(.CROSS.)
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use integrals, only: integrals_kinetic, integrals_locpot
    use kernel, only: DKERN, kernel_create, kernel_destroy
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, ngwf_rep_create, NGWF_HAM, &
        ngwf_rep_destroy, ngwf_ham_create, ngwf_ham_destroy
    use ngwfs, only: ngwfs_initialise
    use parallel_strategy, only: PARAL_INFO
    use paw, only: paw_nonlocal_energies, paw_species_calc_proj_prec_mat, &
         paw_tcore_hartree_on_grid, paw_projector_overlap, &
         paw_tcore_density
    use projectors, only: PROJECTOR_SET
    use pseudopotentials, only: pseudo_get_dij, pseudo_make_structure_factor, &
         pseudopotentials_local_on_grid, pseudopotentials_core_density
    use rundat, only: pub_output_detail, pub_spin_fac, pub_aug_den_dim, &
         pub_debug_on_root, pub_num_spins, pub_num_kpoints, PUB_1K, pub_aug, &
         pub_nlcc, pub_nhat_in_xc, pub_paw, pub_usp, &
         pub_any_nl_proj, pub_maxit_ngwf_cg, pub_lr_phonons_restart, pub_rootname
    use services, only: services_flush
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_scale, sparse_embed_axpy, &
         sparse_embed_trace, sparse_embed_product, sparse_embed_scale, &
         sparse_embed_transpose_structure, sparse_embed_transpose, &
         sparse_embed_copy, sparse_embed_conjugate, sparse_embed_array_scale, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use sparse, only: SPAM3, sparse_element_exists, sparse_get_element, &
         sparse_put_element, sparse_create, sparse_copy, sparse_destroy
    use spherical_wave, only: sw_init
    use xc, only: xc_energy_potential
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_abort, utils_int_to_str, utils_unit
    use visual, only: visual_ngwfs

    ! Arguments
    type(SPAM3_EMBED_ARRAY), intent(inout) :: denskern ! denskern matrix
    type(NGWF_REP), intent(in) :: val_rep
    type(FUNC_BASIS), intent(in) :: val_basis(1)
    type(FUNC_BASIS), intent(inout) :: hub_proj_basis(1)
    type(FUNC_BASIS), intent(in) :: proj_basis(1)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(1)
    type(MODEL), intent(inout) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    type(HUBBARD_MODEL), intent(inout) :: hub


    ! Local Variables
    type(NGWF_HAM) :: ham, joint_ham
    type(FUNC_BASIS), dimension(1) :: FOngwf_basis, joint_basis
    type(NGWF_REP) :: FOngwf_rep, joint_rep

    type(DKERN) :: condkern_new

    type(SPAM3_EMBED) :: temp_mat_1, temp_mat_2, temp_mat_3

    ! r indicates Response (first-order) NGWFs
    ! v indicates Valence NGWFs
    ! j indicates Joint (r+v, in that order on a per-atom basis) NGWFs

    type(SPAM3_EMBED) :: iGpv_overlap, viGp_overlap, pv_overlap, pj_overlap, &
         iGpj_overlap, jiGp_overlap

    type(SPAM3_EMBED), allocatable :: FO_dij(:), dij(:), rho_ij_dproj(:), rho_ij(:)

    type(SPAM3_EMBED) :: gradient_matrices(7), gradient_matrices_cov(7)

    type(SPAM3), allocatable :: kern_array(:), dij_tmp(:), rhoij_tmp(:)

    real(kind=DP) :: dummy_gzero

    real(kind=DP), allocatable :: FO_lhxc_fine(:,:,:,:)
    real(kind=DP), allocatable :: fxc_buffer(:,:,:,:)
    real(kind=DP), allocatable :: FO_locpot(:,:,:), SO_locpot(:,:,:)
    real(kind=DP), allocatable :: FO_tcore_density(:,:,:)
    real(kind=DP), allocatable :: vxc_buffer(:,:,:,:)
    real(kind=DP), allocatable, dimension(:,:,:,:) :: dens_fine
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: nhat_den_grad
    real(kind=DP), allocatable, dimension(:,:,:,:) :: tcore_den_grad

    complex(kind=DP), allocatable, dimension(:,:) :: force_constant_matrix
    real(kind=DP), allocatable, dimension(:,:) :: Born_charges, &
         dielectric_tensor

    real(kind=DP), allocatable, dimension(:) :: frequencies
    complex(kind=DP), allocatable, dimension(:,:) :: eigenvecs

    logical :: present_file

    type(SPAM3_EMBED) :: A_mat, B_mat, C_mat, D_mat

    ! CT indicates complex transpose (or just transpose for the real case)

    type(SPAM3_EMBED) :: respkern_jv(pub_num_spins), respkern_jv_CT(pub_num_spins), &
         auxrespkern_jv(pub_num_spins), auxrespkern_jv_CT(pub_num_spins)
    real(kind=DP) :: lhxc_energy, xc_energy
    integer :: is,ierr, atom, cart, num_ngwfs, perturbation
    logical :: loc_cmplx

    integer :: iter, iter1, iter2, iter1_eff, iter2_eff, iter_ngwfs, iatom, &
        iter3, iter4, proc, perturbation_start

    ! Variables for CG optimisation of first-order NGWFs
    integer :: cg_count
    logical :: converged_cg, line_search_success
    real(kind=DP) :: rms_gradient, trial_length, line_search_coeff
    type(COEF) :: previous_g_dot_g
    type(FUNCTIONS) :: prev_direction_on_grid

    real(kind=DP) :: var_energy, element, temp_energy, volume

    type(SPAM3_EMBED) :: FO_overlap, FO_overlap_jv, FO_ham, FO_ham_jv, SO_locpot_mat,&
         dipole_mat_jv(3)

    complex(kind=DP) :: element_cmplx

    complex(kind=DP), allocatable, dimension(:,:,:,:) :: struct_fac
    complex(kind=DP), allocatable, dimension(:,:,:) :: struct_fac_classical

    ! perturbation directions and weights
    integer, allocatable :: directions(:,:)
    real(kind=DP), allocatable :: weights(:,:)

    real(kind=DP) :: current_qpt(3)

    character(len=256) :: output_file, disp_number
    integer :: output_unit

    logical :: at_gamma

    num_ngwfs = val_basis(1)%num

    if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: Entering LR_PHONONS'

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine lr_phonons_calculate not ready yet for more&
         & than one k-point.')

    loc_cmplx = denskern%m(1,1)%p%iscmplx

    allocate(frequencies(1:3*mdl%nat),stat=ierr)
    call utils_alloc_check('lr_phonons_calculate','frequencies',ierr)
    allocate(eigenvecs(3*mdl%nat,3*mdl%nat),stat=ierr)
    call utils_alloc_check('lr_phonons_calculate','eigenvecs',ierr)

    allocate(directions(3*mdl%nat+3, mdl%nat), stat=ierr) ! + 3 for E-field
    call utils_alloc_check('lr_phonons_calculate', 'directions', ierr)
    allocate(weights(3*mdl%nat+3, mdl%nat), stat=ierr) ! + 3 for E-field
    call utils_alloc_check('lr_phonons_calculate', 'weights', ierr)

    allocate(dens_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_phonons_calculate', 'dens_fine', ierr)

    allocate(FO_lhxc_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_phonons_calculate','FO_lhxc_fine',ierr)

    allocate(FO_locpot(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12), stat=ierr)
    call utils_alloc_check('lr_phonons_calculate','FO_locpot',ierr)

    allocate(SO_locpot(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12), stat=ierr)
    call utils_alloc_check('lr_phonons_calculate','SO_locpot',ierr)

    allocate(FO_tcore_density(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12), stat=ierr)
    call utils_alloc_check('lr_phonons_calculate','FO_tcore_density',ierr)

    allocate(fxc_buffer(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_phonons_calculate','fxc_buffer',ierr)

    allocate(vxc_buffer(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_phonons_calculate','vxc_buffer',ierr)

    ! allocate dij
    if (pub_any_nl_proj.or.pub_paw) then
       allocate(dij(pub_num_spins), stat=ierr)
       call utils_alloc_check('lr_phonons_calculate', 'dij', ierr)
    end if

    if (pub_paw) then
       allocate(nhat_den_grad(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12, pub_num_spins,0:pub_aug_den_dim), &
            stat=ierr)
       call utils_alloc_check('lr_phonons_calculate', 'nhat_den_grad', ierr)

    end if

    ! allocate first-order dij, rho_ij
    ! jcap: this needs to be done always, even if PAW isn't being
    ! used, because otherwise we are passing uninitialised arrays to
    ! various routines
    allocate(rho_ij(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_phonons_calculate', 'rho_ij', ierr)
    allocate(rho_ij_dproj(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_phonons_calculate', 'rho_ij_dproj', ierr)
    allocate(FO_dij(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_phonons_calculate', 'FO_dij', ierr)

    allocate(force_constant_matrix(3*mdl%nat+3,3*mdl%nat+3), stat=ierr)
    call utils_alloc_check('lr_phonons_calculate', 'force_constant_matrix', &
         ierr)
    allocate(Born_charges(3,3*mdl%nat), stat=ierr)
    call utils_alloc_check('lr_phonons_calculate', 'Born_charges', ierr)
    allocate(dielectric_tensor(3,3), stat=ierr)
    call utils_alloc_check('lr_phonons_calculate', 'dielectrix_matrix', ierr)


    ! initialise conduction kernel, first-order NGWFs and joint NGWfs
    call lr_phonons_initialise(denskern, val_rep, val_basis, proj_basis, &
         hub_proj_basis, mdl, hfxstate, ham, hub, nl_projectors, lhxc_fine, &
         FOngwf_basis, FOngwf_rep, joint_basis, joint_rep, joint_ham, &
         condkern_new)

    call comms_barrier

    ! #########################################################################
    ! transform from a 1-atom 1-perturbation basis to all-atoms 1-perturbation
    ! basis, such that the first-order NGWFs can be initialised efficiently
    call lr_phonons_transform('F', directions, weights, mdl%nat)
    ! #########################################################################


    ! Create sparse matrices

    call sparse_embed_transpose_structure(iGpv_overlap%structure, val_rep%sp_overlap)
    call sparse_embed_create(iGpv_overlap)
    call sparse_embed_transpose_structure(iGpj_overlap%structure, &
         joint_rep%sp_overlap)
    call sparse_embed_create(iGpj_overlap)

    call sparse_embed_create(viGp_overlap, val_rep%sp_overlap)
    call sparse_embed_create(jiGp_overlap, joint_rep%sp_overlap)

    call sparse_embed_create(FO_ham, ham%ham(1))
    call sparse_embed_transpose_structure(FO_ham_jv%structure, &
         joint_rep%cross_overlap)
    call sparse_embed_create(FO_ham_jv)
    call sparse_embed_create(FO_overlap, val_rep%overlap)
    call sparse_embed_transpose_structure(FO_overlap_jv%structure, &
         joint_rep%cross_overlap)
    call sparse_embed_create(FO_overlap_jv)

    do iter = 1,3
       call sparse_embed_create(dipole_mat_jv(iter), FO_overlap_jv)
       call sparse_embed_scale(dipole_mat_jv(iter), 0.0_DP)
    end do

    do iter = 1,7
       call sparse_embed_create(gradient_matrices(iter), val_rep%inv_overlap)
       call sparse_embed_create(gradient_matrices_cov(iter), val_rep%inv_overlap, &
            val_rep%overlap)
    end do

    call sparse_embed_create(SO_locpot_mat, ham%ham(1))

    ! calculate the following:
    if (pub_any_nl_proj.or.pub_paw) then
       do is = 1,pub_num_spins
          dij(is)%structure = 'E'
          call sparse_embed_create(dij(is))
       end do
    end if

    if (pub_paw) then
       do is=1,pub_num_spins
          call sparse_embed_create(FO_dij(is), dij(is))
          call sparse_embed_create(rho_ij(is), dij(is))
          call sparse_embed_create(rho_ij_dproj(is), dij(is))
       end do
    end if

    ! gcc32 : scale denskern to correct number of electrons
    if(pub_num_spins==1) call sparse_embed_array_scale(denskern, 2.0_DP)

    if (pub_paw) then
       allocate(dij_tmp(pub_num_spins),stat=ierr)
       call utils_alloc_check('lr_phonons_calculate','dij_tmp',ierr)
       allocate(rhoij_tmp(pub_num_spins),stat=ierr)
       call utils_alloc_check('lr_phonons_calculate','rhoij_tmp',ierr)

       ! the GS rho_ij
       call sparse_embed_extract_from_array(rhoij_tmp,rho_ij)
       call sparse_embed_extract_from_array(kern_array,denskern%m(:,1))
       call aug_projector_denskern(rhoij_tmp, kern_array, val_rep%sp_overlap%p)
       call sparse_embed_destroy_extracted_array(kern_array)
       ! the Dij matrix
       call sparse_embed_extract_from_array(dij_tmp,dij)
       call paw_nonlocal_energies(dij_tmp,rhoij_tmp,mdl%regions(1)%paw_sp, &
            mdl%par)

       call sparse_embed_destroy_extracted_array(dij_tmp,dij,.true.)
       call sparse_embed_destroy_extracted_array(rhoij_tmp,rho_ij,.true.)
       ! add dijhat to dij
       do is=1,pub_num_spins
          call sparse_embed_axpy(dij(is), ham%dijhat(is), 1.0_DP)
       end do
    end if

    if (pub_any_nl_proj.and.(.not.pub_paw)) then
       do is = 1,pub_num_spins
          call pseudo_get_dij(dij(is)%p, mdl%regions(1)%pseudo_sp)
       end do
    end if


    ! obtain ground state density
    dens_fine = 0.0_DP

    ! jcap: use hack for spins
    allocate(kern_array(pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_phonons_calculate','kern_array',ierr)

    call sparse_embed_extract_from_array(kern_array,denskern%m(:,1))
    call density_on_grid(dens_fine, mdl%fine_grid, mdl%dbl_grid, &
         mdl%cell, mdl%fftbox, kern_array, val_rep%ngwf_overlap%p, &
         val_rep%ngwfs_on_grid(1), val_basis(1), &
         val_rep%ngwfs_on_grid(1), val_basis(1))

    if (pub_paw) then
       nhat_den_grad = 0.0_DP
       call augmentation_density_on_grid(nhat_den_grad, mdl%fine_grid, &
            mdl%cell, mdl%regions(1)%pseudo_sp, mdl%regions(1)%paw_sp, &
            mdl%aug_box, kern_array, val_rep%sp_overlap%p)

       dens_fine = dens_fine + nhat_den_grad(:,:,:,:,0)
    end if

    call sparse_embed_destroy_extracted_array(kern_array)



    if (pub_nlcc) then
       do is = 1,pub_num_spins
          dens_fine(:,:,:,is) = dens_fine(:,:,:,is) + &
               mdl%core_density_fine * 0.5_DP * pub_spin_fac
       end do
    end if

    ! get XC potential
    vxc_buffer(:,:,:,:) = 0.0_DP
    if(pub_aug) then
       ! if \hat{\rho} not in XC, remove it from the GS density
       if (.not.pub_nhat_in_xc) then
          dens_fine = dens_fine - nhat_den_grad(:,:,:,:,0)
       end if

       call xc_energy_potential(dens_fine, xc_energy, vxc_buffer, &
            mdl%fine_grid, mdl%cell, pub_aug_den_dim, nhat_den_grad)

       ! if \hat{\rho} previously removed, add it back
       if (.not.pub_nhat_in_xc) then
          dens_fine = dens_fine + nhat_den_grad(:,:,:,:,0)
       end if


    else
       call xc_energy_potential(dens_fine, xc_energy, vxc_buffer, &
            mdl%fine_grid, mdl%cell, 0)
    end if

    ! gcc32 : scale back denskern to half number of electrons
    if(pub_num_spins==1) call sparse_embed_array_scale(denskern, 0.5_DP)

    ! qpt logic        ! WORK AT GAMMA-POINT ONLY AT THE MOMENT
    current_qpt = 0.0_DP
    ! NOTE THAT THIS Q-POINT WILL BE IN CARTESIAN COORDINATES, NOT IN FRACTIONAL
    ! ALWAYS COMPUTE THE Q=0 CASE FIRST, SINCE WE NEED THE DYN-MATRIX AND BORN
    ! CHARGES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( (abs(current_qpt(1)).le.1.0E-3_DP) .and. &
         (abs(current_qpt(2)).le.1.0E-3_DP) .and. &
         (abs(current_qpt(3)).le.1.0E-3_DP) ) then
       at_gamma = .true.
    else
       at_gamma = .false.
    end if

    ! now for each atom and cartesian direction obtain the response kernel and
    ! second order energies

    do is = 1,pub_num_spins
       respkern_jv(is)%structure = 'TDRA'//trim(joint_rep%postfix)
       call sparse_embed_create(respkern_jv(is), iscmplx = loc_cmplx)
       call sparse_embed_transpose_structure(respkern_jv_CT(is)%structure, &
            respkern_jv(is))
       call sparse_embed_create(respkern_jv_CT(is), iscmplx = loc_cmplx)

       auxrespkern_jv(is)%structure = 'TDRA'//trim(joint_rep%postfix)
       call sparse_embed_create(auxrespkern_jv(is), iscmplx = loc_cmplx)
       call sparse_embed_transpose_structure(auxrespkern_jv_CT(is)%structure, &
            auxrespkern_jv(is))
       call sparse_embed_create(auxrespkern_jv_CT(is), iscmplx = loc_cmplx)
    end do


    ! =========================================================================
    ! RESTART: before calculation begins, see if we are restarting from
    ! previous force constants or Born Charges
    force_constant_matrix(:,:) = cmplx_0
    Born_charges(:,:) = 0.0_DP

    if (pub_lr_phonons_restart) then
       ! read on root proc force constant and Born charges
       if (pub_on_root) then
          do perturbation = 1, 3*mdl%nat+3

             ! First the Born charges
             if (perturbation.le.(3*mdl%nat)) then
                present_file = .FALSE.
                write(disp_number,*) perturbation
                write(output_file,'(a)') trim(adjustl(pub_rootname))//&
                     '.Born_charges_'//trim(adjustl(disp_number))

                ! see if file exists
                inquire(file = output_file, EXIST = present_file)
                if (present_file) then
                   output_unit=utils_unit()
                   open(unit=output_unit,form='unformatted',&
                       file=trim(adjustl(output_file)),action='read',&
                       status='old')
                   read(output_unit) Born_charges(:,perturbation)
                   close(output_unit)

                else
                   exit ! exists loop
                end if
             end if


             ! Now the force constant
             present_file = .FALSE.
             write(disp_number,*) perturbation
             write(output_file,'(a)') trim(adjustl(pub_rootname))//&
                  '.force_consts_'//trim(adjustl(disp_number))

             ! see if file exists
             inquire(file = output_file, EXIST = present_file)
             if (present_file) then
                perturbation_start = perturbation + 1
                output_unit=utils_unit()
                open(unit=output_unit,form='unformatted',&
                    file=trim(adjustl(output_file)),action='read',&
                    status='old')
                read(output_unit) force_constant_matrix(perturbation,:)
                close(output_unit)

             else
                exit ! exists loop
             end if

          end do
       end if
       ! communicate to other procs
       call comms_bcast(pub_root_proc_id, force_constant_matrix)
       call comms_bcast(pub_root_proc_id, Born_charges)
       call comms_bcast(pub_root_proc_id, perturbation_start)
    else
       perturbation_start = 1
    end if
    ! End of restart
    ! ========================================================================

    call comms_barrier

    ! ======================== LOOP OVER PERTURBATIONS =======================

    do perturbation = perturbation_start,3*mdl%nat+3

       ! E-field perturbation (perturbation > 3*nat) relevant only at Gamma
       if ((perturbation.gt.(3*mdl%nat)).and.(.not.at_gamma)) cycle

       call data_set_to_zero(FOngwf_rep%ngwfs_on_grid(1))
       call data_set_to_zero(joint_rep%ngwfs_on_grid(1))

       FO_locpot = 0.0_DP
       SO_locpot = 0.0_DP
       FO_tcore_density = 0.0_DP

       if (perturbation.le.(3*mdl%nat)) then ! For atomic displacements only

          allocate(struct_fac(mdl%num_pspecies, mdl%fine_grid%ld3, &
               mdl%fine_grid%ld2, mdl%fine_grid%max_slabs23),stat=ierr)
          call utils_alloc_check('lr_phonons_calculate','struct_fac',ierr)

          allocate(struct_fac_classical(0,0,0), stat=ierr)
          call utils_alloc_check('lr_phonons_calculate','struct_fac_classical',&
               ierr)

          ! Calculate structure factor for first-order pseudopotential
          ! and core density
          call pseudo_make_structure_factor(struct_fac, mdl%elements, &
               mdl%fine_grid, mdl%num_pspecies, directions(perturbation,:), &
               weights(perturbation,:))

          if (pub_paw) then
             call paw_tcore_hartree_on_grid(FO_locpot, dummy_gzero, struct_fac,&
                  struct_fac_classical, mdl%regions(1)%paw_sp, mdl%fine_grid, &
                  mdl%par)

             call paw_tcore_density(FO_tcore_density, mdl%fine_grid, mdl%cell, &
                  mdl%regions(1)%paw_sp, mdl%par, &
                  directions(perturbation,:), weights(perturbation,:))
          else
             call pseudopotentials_local_on_grid(FO_locpot, dummy_gzero, &
                  struct_fac, struct_fac_classical, mdl%fine_grid, mdl%cell, &
                  mdl%regions(1)%pseudo_sp, mdl%par, mdl%par%nat_classical)

             call pseudopotentials_core_density(FO_tcore_density, struct_fac, &
                  mdl%regions(1)%pseudo_sp, mdl%fine_grid, par=mdl%par)
          end if


          struct_fac = cmplx(0.0_DP, 0.0_DP, kind = DP)
          struct_fac_classical = cmplx(0.0_DP, 0.0_DP, kind = DP)

          ! Calculate structure factor for second-order pseudopotential
          call pseudo_make_structure_factor(struct_fac,mdl%elements, &
               mdl%fine_grid, mdl%par%num_pspecies, directions(perturbation,:), &
               weights(perturbation,:), directions(perturbation,:), &
               weights(perturbation,:))

          if (pub_paw) then
             call paw_tcore_hartree_on_grid(SO_locpot, dummy_gzero, struct_fac, &
                  struct_fac_classical, mdl%regions(1)%paw_sp, mdl%fine_grid, &
                  mdl%par)
          else
             call pseudopotentials_local_on_grid(SO_locpot, dummy_gzero, &
                  struct_fac, struct_fac_classical, mdl%fine_grid, mdl%cell, &
                  mdl%regions(1)%pseudo_sp, mdl%par, mdl%par%nat_classical)
          end if

          deallocate(struct_fac,stat=ierr)
          call utils_alloc_check('lr_phonons_calculate','struct_fac',ierr)

          deallocate(struct_fac_classical,stat=ierr)
          call utils_dealloc_check('lr_phonons_calculate', &
               'struct_fac_classical', ierr)

          ! Calculate the second-order pseudopotential in val-val rep
          call integrals_locpot(SO_locpot_mat%p, val_rep%ngwfs_on_grid(1), &
               val_basis(1), val_rep%ngwfs_on_grid(1), val_basis(1), &
               mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, SO_locpot)
       end if

       ! Build overlap between 'v' and 'j' NGWFs and first-order projectors
       ! Non-zero only for atomic perturbations, since projectors are not
       ! moved by infinitesimal E-fields
       if (pub_any_nl_proj.or.pub_paw) then
          call sparse_embed_scale(viGp_overlap, 0.0_DP)
          call sparse_embed_scale(jiGp_overlap, 0.0_DP)

          if (perturbation.le.(3*mdl%nat)) then ! for atomic displacements only
             call grad_proj_ovlp(mdl, viGp_overlap%p, val_rep%ngwfs_on_grid(1), &
                  val_basis(1), proj_basis(1), nl_projectors(1), &
                  directions(perturbation,:), weights(perturbation,:))

             call grad_proj_ovlp(mdl, jiGp_overlap%p, joint_rep%ngwfs_on_grid(1), &
                  joint_basis(1), proj_basis(1), nl_projectors(1), &
                  directions(perturbation,:), weights(perturbation,:))
          end if

          call sparse_embed_transpose(iGpv_overlap, viGp_overlap)
          call sparse_embed_transpose(iGpj_overlap, jiGp_overlap)
       end if


       ! Initialise CG variables
       converged_cg = .FALSE.
       previous_g_dot_g%iscmplx = val_rep%ngwfs_on_grid(1)%iscmplx
       call data_set_to_zero(previous_g_dot_g)
       trial_length = 10.0_DP
       line_search_success = .TRUE.
       line_search_coeff = 0.15_DP
       cg_count = 0

       call data_functions_alloc(prev_direction_on_grid, &
            FOngwf_basis(1)%size_on_grid, iscmplx=val_rep%ngwfs_on_grid(1)%iscmplx)
       call data_set_to_zero(prev_direction_on_grid)


       ! ==================== OPTIMISE FIRST-ORDER NGWFS =====================
       do iter_ngwfs = 1, pub_maxit_ngwf_cg

          if (pub_on_root) then
             write(stdout,'(a,i4.3,a,i5.4,a)')'####################### FO-NGWF &
                  &CG ',iter_ngwfs,' | PERTURBATION ',perturbation, &
                  ' ###################'
          end if


          ! E-field: set weights to 1/sqrt(Natoms) only for the
          ! initialisation of the first-order NGWFs
          if ((iter_ngwfs == 1).and.(perturbation.gt.(3*mdl%nat))) then
             do iter1 = (3*mdl%nat + 1), (3*mdl%nat +3)
                weights(iter1, :) = 1.0_DP / sqrt(real(mdl%nat,kind=DP))
                directions(iter1, :) = mod(iter1,(3*mdl%nat))
             end do
          end if

          call calc_FO_NGWFs(mdl, hfxstate, val_rep, val_basis(1), directions, &
               weights, perturbation, proj_basis(1), nl_projectors(1), ham, &
               dij, rho_ij, dens_fine, denskern, hub_proj_basis(1), hub, &
               lhxc_fine, FOngwf_basis(1), FOngwf_rep, FO_lhxc_fine, FO_dij, &
               viGp_overlap, respkern_jv, auxrespkern_jv, vxc_buffer, &
               fxc_buffer, FO_ham, FO_ham_jv, iter_ngwfs, gradient_matrices, &
               gradient_matrices_cov, converged_cg, previous_g_dot_g, &
               rms_gradient, line_search_success, trial_length, &
               prev_direction_on_grid, line_search_coeff, cg_count, &
               joint_basis(1), joint_rep, joint_ham, condkern_new%kern, FO_locpot,&
               FO_tcore_density, SO_locpot_mat, dipole_mat_jv, at_gamma)


          ! E-field: set weights to 0 for the remainder of the CG
          if ((iter_ngwfs == 1).and.(perturbation.gt.(3*mdl%nat))) then
             do iter1 = (3*mdl%nat + 1), (3*mdl%nat +3)
                weights(iter1, :) = 0.0_DP
             end do
          end if

          call comms_barrier

          ! Calculate dipole matrix in j-v representation, only valid
          ! for E-field perturbations
          if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then
             ! calculate dipole matrix in jv rep
             call sparse_embed_transpose_structure(temp_mat_1%structure, &
                  joint_rep%cross_overlap)
             call sparse_embed_create(temp_mat_1)
             call sparse_embed_transpose(temp_mat_1, joint_rep%cross_overlap)

             call calculate_dipole_mat(dipole_mat_jv, mdl, joint_rep, &
                  joint_basis, val_rep, val_basis, temp_mat_1, proj_basis, &
                  nl_projectors, .false.)

             call sparse_embed_destroy(temp_mat_1)
          end if

          ! CREATE MATRICES FOR RESPONSE KERNEL OPTIMISATION (STERNHEIMER EQ.)
          ! D_mat
          call sparse_embed_create(temp_mat_1, val_rep%overlap, denskern%m(1,1))
          call sparse_embed_product(temp_mat_1, val_rep%overlap, denskern%m(1,1))
          call sparse_embed_create(temp_mat_2, temp_mat_1, ham%ham(1))
          call sparse_embed_product(temp_mat_2, temp_mat_1, ham%ham(1))
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_create(D_mat, temp_mat_2, val_rep%inv_overlap)
          call sparse_embed_product(D_mat, temp_mat_2, val_rep%inv_overlap)
          call sparse_embed_destroy(temp_mat_2)

          ! B_mat
          call sparse_embed_create(B_mat, D_mat)
          call sparse_embed_product(B_mat, val_rep%overlap, denskern%m(1,1))

          ! A_mat
          call sparse_embed_create(temp_mat_1, condkern_new%kern%m(1,1), &
               joint_ham%ham(1))
          call sparse_embed_product(temp_mat_1, condkern_new%kern%m(1,1), &
               joint_ham%ham(1))
          call sparse_embed_create(temp_mat_2, temp_mat_1, condkern_new%kern%m(1,1))
          call sparse_embed_product(temp_mat_2, temp_mat_1, condkern_new%kern%m(1,1))
          call sparse_embed_create(A_mat, temp_mat_2, joint_rep%overlap)
          call sparse_embed_product(A_mat, temp_mat_2, joint_rep%overlap)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)

          ! C_mat
          call sparse_embed_create(C_mat, A_mat)
          call sparse_embed_product(C_mat, condkern_new%kern%m(1,1), &
               joint_rep%overlap)
          call sparse_embed_scale(C_mat, -1.0_DP)

          ! Reset response kernels, first-order hamiltonian, first-order overlap
          ! and first-order Dij-matrix
          if (iter_ngwfs == 1) then
             do is = 1, pub_num_spins
                call sparse_embed_scale(respkern_jv(is), 0.0_DP)
                call sparse_embed_scale(auxrespkern_jv(is), 0.0_DP)
             end do
          end if
          call sparse_embed_scale(FO_ham, 0.0_DP)
          call sparse_embed_scale(FO_ham_jv, 0.0_DP)
          call sparse_embed_scale(FO_overlap, 0.0_DP)
          call sparse_embed_scale(FO_overlap_jv, 0.0_DP)

          do is = 1,pub_num_spins
             if (pub_paw) then
                call sparse_embed_scale(FO_dij(is),0.0_DP)
             end if
             fxc_buffer(:,:,:,is) = 0.0_DP
             FO_lhxc_fine(:,:,:,is) = 0.0_DP
          end do

          if (pub_on_root) then
             write(stdout,'(/a)') '============== RESPONSE KERNEL OPTIMISATION &
                  & ==================================='
          end if

          ! Optimise response kernel for the current set of first-order NGWFs
          call lr_phonons_solve_rkn(mdl, val_basis, val_rep, ham, denskern, &
               A_mat, B_mat, C_mat, D_mat, respkern_jv, auxrespkern_jv, &
               joint_basis, joint_rep, joint_ham, condkern_new%kern, &
               directions(perturbation,:), weights(perturbation,:), &
               perturbation, proj_basis, nl_projectors, dij, dens_fine, &
               FO_tcore_density, FO_locpot, rho_ij, lhxc_fine, FO_lhxc_fine, &
               fxc_buffer, vxc_buffer, FO_ham, FO_ham_jv, FO_dij, &
               iGpv_overlap, FO_overlap, FO_overlap_jv, dipole_mat_jv, &
               at_gamma)

          if (pub_on_root) then
             write(stdout,'(/a)') '============================================&
                  &===================================='
          end if

          ! Destroy A_mat, B_mat, C_mat, D_mat
          call sparse_embed_destroy(A_mat)
          call sparse_embed_destroy(B_mat)
          call sparse_embed_destroy(C_mat)
          call sparse_embed_destroy(D_mat)

          ! Calculate adjoint response/auxresponse kernels
          do is = 1, pub_num_spins
             call sparse_embed_transpose(respkern_jv_CT(is), respkern_jv(is))
             call sparse_embed_transpose(auxrespkern_jv_CT(is), auxrespkern_jv(is))
             if (respkern_jv(is)%iscmplx) then
                call sparse_embed_conjugate(respkern_jv_CT(is))
                call sparse_embed_conjugate(auxrespkern_jv_CT(is))
             end if
          end do

          ! Exit if first-order NGWFs are converged
          if (converged_cg) then
             if (pub_on_root) then
                write(stdout,*) 'CG CONVERGED FOR PERTURBATION ', perturbation
             end if
             exit
          else
             call internal_create_gradient_mat
             call comms_barrier
          end if


       end do ! iterations over NGWF opt

       ! ================= END OPTIMISATION FOR FIRST-ORDER NGWFS =============

       ! Calculate Born charges only at Gamma and only for atomic displacements
       if ((perturbation.le.(3*mdl%nat)).and.at_gamma) then

          Born_charges(:,perturbation) = 0.0_DP
          call dynamical_charges_BERRY(Born_charges(:, perturbation), mdl, &
               val_rep, val_basis, directions(perturbation,:), &
               weights(perturbation,:), perturbation, denskern%m(1,1), &
               respkern_jv(1), proj_basis, nl_projectors, viGp_overlap, &
               joint_basis, joint_rep, FO_overlap)
       end if

       ! Get electronic contribution of non-variational second-order energy
       force_constant_matrix(perturbation,:) = cmplx(0.0_DP,0.0_DP,kind=DP)
       call second_order_energy(mdl, val_rep, val_basis, ham, joint_basis,&
            joint_rep, joint_ham, directions, weights, perturbation, &
            dij, FO_dij, FO_ham, FO_ham_jv, lhxc_fine, FO_lhxc_fine, &
            vxc_buffer, fxc_buffer, denskern%m(1,1), respkern_jv(1), &
            proj_basis, nl_projectors, iGpv_overlap, dipole_mat_jv, at_gamma, &
            force_constant_matrix)

       ! Add ion (Ewald) contribution to non-variational second-order energy
       if (perturbation.le.(3*mdl%nat)) then ! only for atomic displacements
          call second_order_Ewald(force_constant_matrix, mdl, directions, &
                weights, perturbation, current_qpt)
       end if

       ! if (pub_on_root) then
       !    do iter1 = 1,3*par%nat
       !       write(stdout,*) 'FORCE CONSTANT TERM ', perturbation, iter1, &
       !            force_constant_matrix(perturbation,iter1)
       !    end do
       ! end if


       call data_functions_dealloc(prev_direction_on_grid)

       ! ======================================================================
       ! RESTART: Write force constant and Born charge to file
       if (pub_on_root) then
          write(disp_number,*) perturbation
          ! force constant
          write(output_file,'(a)') trim(adjustl(pub_rootname))//&
               '.force_consts_'//trim(adjustl(disp_number))
          output_unit=utils_unit()
          open(unit=output_unit,form='unformatted',&
               file=trim(adjustl(output_file)),action='write')
          write(output_unit) force_constant_matrix(perturbation,:)
          close(output_unit)


          ! Born charge if atomic displacement
          if (perturbation.le.(3*mdl%nat)) then
             write(output_file,'(a)') trim(adjustl(pub_rootname))//&
                  '.Born_charges_'//trim(adjustl(disp_number))
             output_unit=utils_unit()
             open(unit=output_unit,form='unformatted',&
                  file=trim(adjustl(output_file)),action='write')
             write(output_unit) Born_charges(:,perturbation)
             close(output_unit)
          end if

       end if
       ! End RESTART
       ! ======================================================================

    end do ! end loop over perturbations

    ! ======================== LOOP OVER PERTURBATIONS =======================

    ! Get dielectric tensor out if q = 0
    ! The "-" appears since the polarisation contains the dipole matrix (D)
    ! with a minus in front

    volume = (mdl%cell%a1 .CROSS. mdl%cell%a2).DOT.mdl%cell%a3
    if (at_gamma) then
       dielectric_tensor(:,:) = -4.0_DP * (PI / volume) * &
            force_constant_matrix(3*mdl%nat+1:3*mdl%nat+3, &
                 3*mdl%nat+1:3*mdl%nat+3)
       do iter1 = 1,3
          dielectric_tensor(iter1,iter1) = 1.0_DP + dielectric_tensor(iter1, &
               iter1)
       end do
    end if


    ! #########################################################################
    ! transform the force constant matrix and Born charges from the
    ! compound perturbation basis to the initial 1-atom 1 direction basis
    call lr_phonons_transform('B', directions, weights, mdl%nat, at_gamma, &
         Born_charges, force_constant_matrix)
    ! #########################################################################

   ! Deal with Born charges / Dielectric matrix only at Gamma
   if (at_gamma) then

       ! Enforce acoustic sum rule on Born Charges
       call  lr_phonons_born_asr(Born_charges, mdl)

       ! print Born Charges
       if (pub_on_root) then
          write(stdout,'()')
          write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
               &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
          write(stdout,'(a)') '            Born charges                     &
               &                 '
          write(stdout,'()')
          do iter1 = 1,3*mdl%nat
             write(stdout,'(3x,i5,a,3(6x,f16.11))') iter1, ":", &
                  Born_charges(:,iter1)
          end do
          write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
               &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
       end if

       ! Print Dielectric tensor
       if (pub_on_root) then
          write(stdout,'()')
          write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
               &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
          write(stdout,'(a)') '            Dielectric permitivity tensor    &
               &                 '
          write(stdout,'()')
          do iter1 = 1,3
             write(stdout,'(3x,i5,a,3(6x,f16.11))') iter1, ":", &
                  dielectric_tensor(iter1,:)
          end do
          write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
               &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
       end if

    end if


    ! Creates, diagonalises and corrects (for sum-rules, if at Gamma)
    ! the dynamical matrix. Also adds the Non-Analytical term if q -> 0
    call lr_phonons_dynmat(force_constant_matrix(1:3*mdl%nat,1:3*mdl%nat), &
         eigenvecs, frequencies, current_qpt, at_gamma, dielectric_tensor, &
         Born_charges, volume, mdl)

    ! =======================================================================



    ! Destroy and deallocate sparse matrices
    do is = 1, pub_num_spins
       call sparse_embed_destroy(respkern_jv(is))
       call sparse_embed_destroy(auxrespkern_jv(is))
       call sparse_embed_destroy(respkern_jv_CT(is))
       call sparse_embed_destroy(auxrespkern_jv_CT(is))
    end do

    call sparse_embed_destroy(FO_overlap)
    call sparse_embed_destroy(FO_overlap_jv)
    call sparse_embed_destroy(FO_ham)
    call sparse_embed_destroy(FO_ham_jv)

    call kernel_destroy(condkern_new)

    do iter = 1,7
       call sparse_embed_destroy(gradient_matrices(iter))
       call sparse_embed_destroy(gradient_matrices_cov(iter))
    end do

    call sparse_embed_destroy(SO_locpot_mat)

    do iter = 1,3
       call sparse_embed_destroy(dipole_mat_jv(iter))
    end do

    call data_functions_dealloc(FOngwf_rep%ngwfs_on_grid(1))
    call function_basis_deallocate(FOngwf_basis(1))
    call ngwf_rep_destroy(FOngwf_rep)

    call data_functions_dealloc(joint_rep%ngwfs_on_grid(1))
    call function_basis_deallocate(joint_basis(1))
    call ngwf_rep_destroy(joint_rep)

    if (pub_paw.or.pub_any_nl_proj) then
       call sparse_embed_destroy(viGp_overlap)
       do is = 1,pub_num_spins
          call sparse_embed_destroy(dij(is))
       end do
    end if

    if (pub_paw) then
       do is = 1,pub_num_spins
          call sparse_embed_destroy(rho_ij(is))
          call sparse_embed_destroy(rho_ij_dproj(is))
          call sparse_embed_destroy(FO_dij(is))
       end do
    end if

    if(pub_paw.or.pub_any_nl_proj) then
       call sparse_embed_destroy(iGpv_overlap)
    end if

    if (pub_paw.or.pub_any_nl_proj) then
       deallocate(dij, stat=ierr)
       call utils_dealloc_check('lr_phonons_calculate', 'dij', ierr)
    end if


    deallocate(directions, stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate', 'directions', ierr)
    deallocate(weights, stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate', 'weights', ierr)

    deallocate(frequencies,stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate','frequencies',ierr)
    deallocate(eigenvecs,stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate','eigenvecs',ierr)


    ! deallocate first-order dij(pub_num_spins)
    deallocate(FO_dij, stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate', 'FO_dij', ierr)
    deallocate(rho_ij, stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate', 'rho_ij', ierr)
    deallocate(rho_ij_dproj, stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate', 'rho_ij_dproj', ierr)
    if (pub_paw) then
       deallocate(nhat_den_grad, stat=ierr)
       call utils_dealloc_check('lr_phonons_calculate', 'nhat_den_grad', ierr)
    end if

    deallocate(force_constant_matrix, stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate', 'force_constant_matrix', &
         ierr)
    deallocate(Born_charges, stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate','Born_charges',ierr)
    deallocate(FO_lhxc_fine, stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate','FO_lhxc_fine',ierr)
    deallocate(FO_locpot, stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate','FO_locpot',ierr)
    deallocate(SO_locpot, stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate','SO_locpot',ierr)
    deallocate(FO_tcore_density, stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate','FO_tcore_density',ierr)
    deallocate(fxc_buffer, stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate','fxc_buffer',ierr)
    deallocate(vxc_buffer, stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate','vxc_buffer',ierr)
    deallocate(dens_fine, stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate', 'dens_fine', ierr)
    deallocate(dielectric_tensor, stat=ierr)
    call utils_dealloc_check('lr_phonons_calculate', 'dielectric_tensor', ierr)

  contains

    subroutine internal_create_gradient_mat

      ! This internal subroutine calculate the matrix coefficients that go
      ! into the optimisation of the First-order NGWFs

      implicit none

          type(SPAM3_EMBED) :: gradient_matrices_vj(4), gradient_matrices_vj_cov(4), &
               gradient_matrices_jj(2), gradient_matrices_jj_cov(2) , &
               X_mat_jv
          integer :: iat, loc_iat, jat, loc_jat, ingwf, jngwf, loc_ingwf, &
               loc_jngwf, ingwf_eff, jngwf_eff

          do iter = 1,4
             call sparse_embed_create(temp_mat_1, val_rep%inv_overlap, &
                  joint_rep%cross_overlap)
             call sparse_embed_create(gradient_matrices_vj(iter), temp_mat_1, &
                  joint_rep%inv_overlap, iscmplx = .false.)
             call sparse_embed_create(gradient_matrices_vj_cov(iter), temp_mat_1, &
                  iscmplx = .false.)
             call sparse_embed_destroy(temp_mat_1)
          end do

          do iter = 1,2
             call sparse_embed_create(gradient_matrices_jj(iter), &
                  joint_rep%inv_overlap, iscmplx = .false.)
             call sparse_embed_create(gradient_matrices_jj_cov(iter), &
                  joint_rep%inv_overlap, joint_rep%overlap, &
                  iscmplx = .false.)
          end do

          do iter = 1,7
             call sparse_embed_scale(gradient_matrices(iter), 0.0_DP)
             call sparse_embed_scale(gradient_matrices_cov(iter), 0.0_DP)
          end do

          !=================================================================
          !============ CALCULATE GRADIENT MATRICES ========================

          ! calculate X_mat_jv first
          call sparse_embed_create(temp_mat_1, joint_ham%ham(1), respkern_jv(1))
          call sparse_embed_create(X_mat_jv, temp_mat_1, val_rep%overlap)
          call sparse_embed_destroy(temp_mat_1)

          call sparse_embed_create(temp_mat_1, joint_ham%ham(1), respkern_jv(1))
          call sparse_embed_product(temp_mat_1, joint_ham%ham(1), respkern_jv(1))
          call sparse_embed_product(X_mat_jv, temp_mat_1, val_rep%overlap)
          call sparse_embed_destroy(temp_mat_1)

          call sparse_embed_create(temp_mat_1, joint_rep%overlap, respkern_jv(1))
          call sparse_embed_product(temp_mat_1, joint_rep%overlap, respkern_jv(1))
          call sparse_embed_create(temp_mat_2, temp_mat_1, ham%ham(1))
          call sparse_embed_product(temp_mat_2, temp_mat_1, ham%ham(1))
          call sparse_embed_axpy(X_mat_jv, temp_mat_2, -1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)

          call sparse_embed_axpy(X_mat_jv, FO_ham_jv, 1.0_DP)

          if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then
             call sparse_embed_axpy(X_mat_jv, &
                  dipole_mat_jv(mod(perturbation,3*mdl%nat)), 1.0_DP)
          end if

          call sparse_embed_create(temp_mat_1, FO_overlap_jv, val_rep%inv_overlap)
          call sparse_embed_product(temp_mat_1, FO_overlap_jv, val_rep%inv_overlap)
          call sparse_embed_create(temp_mat_2, temp_mat_1, ham%ham(1))
          call sparse_embed_product(temp_mat_2, temp_mat_1, ham%ham(1))
          call sparse_embed_axpy(X_mat_jv, temp_mat_2, -1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
          ! end of X_mat

          call comms_barrier

          ! grad_jj(1) -> <r|S|theta>
          call sparse_embed_create(temp_mat_1, respkern_jv(1), ham%ham(1))
          call sparse_embed_product(temp_mat_1, respkern_jv(1), ham%ham(1))
          call sparse_embed_product(gradient_matrices_jj(1), temp_mat_1, &
               respkern_jv_CT(1))
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_scale(gradient_matrices_jj(1), -2.0_DP)
          call sparse_embed_product(gradient_matrices_jj_cov(1), &
               gradient_matrices_jj(1), joint_rep%overlap)

          ! grad_jj(2) -> <r|H|theta>
          call sparse_embed_create(temp_mat_1, respkern_jv(1), val_rep%overlap)
          call sparse_embed_product(temp_mat_1, respkern_jv(1), val_rep%overlap)
          call sparse_embed_product(gradient_matrices_jj(2), temp_mat_1, &
               respkern_jv_CT(1))
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_scale(gradient_matrices_jj(2), 2.0_DP)
          call sparse_embed_product(gradient_matrices_jj_cov(2), &
               gradient_matrices_jj(2), joint_rep%overlap)

          ! grad_vj(1) -> <r|S|phi>
          call sparse_embed_create(temp_mat_1, denskern%m(1,1), &
               joint_rep%cross_overlap)
          call sparse_embed_product(temp_mat_1, denskern%m(1,1), &
               joint_rep%cross_overlap)
          call sparse_embed_create(temp_mat_2, temp_mat_1, joint_rep%inv_overlap)
          call sparse_embed_product(temp_mat_2, temp_mat_1, joint_rep%inv_overlap)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_create(temp_mat_1, temp_mat_2, X_mat_jv)
          call sparse_embed_product(temp_mat_1, temp_mat_2, X_mat_jv)
          call sparse_embed_destroy(temp_mat_2)
          call sparse_embed_create(temp_mat_2, temp_mat_1, denskern%m(1,1))
          call sparse_embed_product(temp_mat_2, temp_mat_1, denskern%m(1,1))
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_create(temp_mat_1, temp_mat_2, val_rep%overlap)
          call sparse_embed_product(temp_mat_1, temp_mat_2, val_rep%overlap)
          call sparse_embed_destroy(temp_mat_2)
          call sparse_embed_product(gradient_matrices_vj(1), temp_mat_1, &
               auxrespkern_jv_CT(1))
          call sparse_embed_scale(gradient_matrices_vj(1), -2.0_DP)
          call sparse_embed_destroy(temp_mat_1)

          call sparse_embed_product(gradient_matrices_vj_cov(1), &
               gradient_matrices_vj(1), joint_rep%overlap)

          ! grad_vj(2) -> <r|FO_overlap|phi>
          call sparse_embed_create(temp_mat_1, val_rep%inv_overlap, ham%ham(1))
          call sparse_embed_product(temp_mat_1, val_rep%inv_overlap, ham%ham(1))
          call sparse_embed_product(gradient_matrices_vj(2), temp_mat_1, &
               respkern_jv_CT(1))
          call sparse_embed_scale(gradient_matrices_vj(2), -2.0_DP)
          call sparse_embed_destroy(temp_mat_1)

          call sparse_embed_product(gradient_matrices_vj_cov(2), &
               gradient_matrices_vj(2), joint_rep%overlap)

          ! grad_vj(3) -> <r|FO_ham|phi>
          call sparse_embed_copy(gradient_matrices_vj(3), respkern_jv_CT(1))
          call sparse_embed_scale(gradient_matrices_vj(3), 2.0_DP)

          call sparse_embed_product(gradient_matrices_vj_cov(3), &
               gradient_matrices_vj(3), joint_rep%overlap)

          ! grad_vj(4) -> <r|FO_E\cdot D|phi>
          if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then
             call sparse_embed_copy(gradient_matrices_vj(4), respkern_jv_CT(1))
             call sparse_embed_scale(gradient_matrices_vj(4), 2.0_DP)

             call sparse_embed_product(gradient_matrices_vj_cov(4), &
                  gradient_matrices_vj(4), joint_rep%overlap)
          else
             call sparse_embed_scale(gradient_matrices_vj(4), 0.0_DP)
             call sparse_embed_scale(gradient_matrices_vj_cov(4), 0.0_DP)
          end if
          ! ---------------------------------------------------------------

          call comms_barrier

          ! EXTRACT ONLY THE BLOCKS WE NEED, SINCE THE MATRIX COEFFICIENTS
          ! ARE ALWAYS REPRESENTED IN THE FIRST-ORDER BASIS ON THE RIGHT

          call lr_phonons_get_basis_block(gradient_matrices_jj(1)%p,&
               gradient_matrices(1)%p, joint_basis(1), joint_basis(1), &
               FOngwf_basis(1), FOngwf_basis(1), 1, 1)

          call lr_phonons_get_basis_block(gradient_matrices_jj_cov(1)%p, &
               gradient_matrices_cov(1)%p, joint_basis(1), joint_basis(1), &
               FOngwf_basis(1), FOngwf_basis(1), 1, 1)

          call lr_phonons_get_basis_block(gradient_matrices_jj(2)%p,&
               gradient_matrices(2)%p, joint_basis(1), joint_basis(1), &
               FOngwf_basis(1), FOngwf_basis(1), 1, 1)

          call lr_phonons_get_basis_block(gradient_matrices_jj_cov(2)%p,&
               gradient_matrices_cov(2)%p, joint_basis(1), joint_basis(1), &
               FOngwf_basis(1), FOngwf_basis(1), 1, 1)

          call lr_phonons_get_basis_block(gradient_matrices_jj(1)%p,&
               gradient_matrices(3)%p, joint_basis(1), joint_basis(1), &
               val_basis(1), FOngwf_basis(1), 2, 1)

          call lr_phonons_get_basis_block(gradient_matrices_jj_cov(1)%p,&
               gradient_matrices_cov(3)%p, joint_basis(1), joint_basis(1), &
               val_basis(1), FOngwf_basis(1), 2, 1)

          call lr_phonons_get_basis_block(gradient_matrices_jj(2)%p,&
               gradient_matrices(4)%p, joint_basis(1), joint_basis(1), &
               val_basis(1), FOngwf_basis(1), 2, 1)

          call lr_phonons_get_basis_block(gradient_matrices_jj_cov(2)%p,&
               gradient_matrices_cov(4)%p, joint_basis(1), joint_basis(1), &
               val_basis(1), FOngwf_basis(1), 2, 1)

          call sparse_embed_create(temp_mat_1, gradient_matrices(3))
          call sparse_embed_create(temp_mat_2, gradient_matrices_cov(3))

          call lr_phonons_get_basis_block(gradient_matrices_vj(1)%p,&
               temp_mat_1%p, val_basis(1), joint_basis(1), &
               val_basis(1), FOngwf_basis(1), 1, 1)

          call lr_phonons_get_basis_block(gradient_matrices_vj_cov(1)%p,&
               temp_mat_2%p, val_basis(1), joint_basis(1), &
               val_basis(1), FOngwf_basis(1), 1, 1)

          call lr_phonons_get_basis_block(gradient_matrices_vj(2)%p,&
               gradient_matrices(5)%p, val_basis(1), joint_basis(1), &
               val_basis(1), FOngwf_basis(1), 1, 1)

          call lr_phonons_get_basis_block(gradient_matrices_vj_cov(2)%p,&
               gradient_matrices_cov(5)%p, val_basis(1), joint_basis(1), &
               val_basis(1), FOngwf_basis(1), 1, 1)

          call lr_phonons_get_basis_block(gradient_matrices_vj(3)%p,&
               gradient_matrices(6)%p, val_basis(1), joint_basis(1), &
               val_basis(1), FOngwf_basis(1), 1, 1)

          call lr_phonons_get_basis_block(gradient_matrices_vj_cov(3)%p,&
               gradient_matrices_cov(6)%p, val_basis(1), joint_basis(1), &
               val_basis(1), FOngwf_basis(1), 1, 1)

          call lr_phonons_get_basis_block(gradient_matrices_vj(4)%p,&
               gradient_matrices(7)%p, val_basis(1), joint_basis(1), &
               val_basis(1), FOngwf_basis(1), 1, 1)

          call lr_phonons_get_basis_block(gradient_matrices_vj_cov(4)%p,&
               gradient_matrices_cov(7)%p, val_basis(1), joint_basis(1), &
               val_basis(1), FOngwf_basis(1), 1, 1)

          call sparse_embed_axpy(gradient_matrices(3), temp_mat_1, 1.0_DP)
          call sparse_embed_axpy(gradient_matrices_cov(3), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)

          call sparse_embed_destroy(X_mat_jv)
          do iter = 1,4
             call sparse_embed_destroy(gradient_matrices_vj(iter))
             call sparse_embed_destroy(gradient_matrices_vj_cov(iter))
          end do

          do iter = 1,2
             call sparse_embed_destroy(gradient_matrices_jj(iter))
             call sparse_embed_destroy(gradient_matrices_jj_cov(iter))
          end do

          call comms_barrier

          ! multiply all contravar and covar gradient matrices by 2*s
          do iter = 1,7
             call sparse_embed_scale(gradient_matrices(iter), 2.0_DP*pub_spin_fac)
             call sparse_embed_scale(gradient_matrices_cov(iter), &
                  2.0_DP*pub_spin_fac)
          end do

          call comms_barrier

          !============= END GRADIENT MATRICES =============================
          !=================================================================

    end subroutine internal_create_gradient_mat


  end subroutine lr_phonons_calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_phonons_initialise(denskern, val_rep, val_basis, &
       proj_basis, hub_proj_basis, mdl, hfxstate, ham, hub, nl_projectors, &
       lhxc_fine, FOngwf_basis, FOngwf_rep, joint_basis, joint_rep, joint_ham, &
       condkern_new)

    ! ========================================================================!
    ! Initialises the spheres, tightboxes and representations for the first-  !
    ! order NGWFs and joint NGWFs                                             !
    ! ------------------------------------------------------------------------!
    ! Modified for embedding and to remove pub_par by Joseph Prentice,        !
    ! September 2018                                                          !
    ! ========================================================================!

    use datatypes, only: data_functions_alloc, FUNCTIONS, &
         data_functions_dealloc, data_set_to_zero
    use comms, only: comms_barrier, pub_on_root
    use density, only: density_on_grid
    use electronic_init, only: electronic_init_denskern
    use function_basis, only: FUNC_BASIS, function_basis_allocate, &
         function_basis_deallocate, function_basis_gath_all_tbs, &
         function_basis_estimate_size, function_basis_init_tight_boxes, &
         function_basis_copy_spheres, function_basis_distribute, &
         function_basis_init_spheres, function_basis_init_uni_tb
    use constants, only: stdout, NORMAL, paw_en_size, UP, DN, VERBOSE, LONG
    use hamiltonian, only: hamiltonian_dens_dep_matrices, &
         hamiltonian_build_matrix
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: DKERN, kernel_create
    use model_type, only: MODEL
    use ngwfs, only: ngwfs_initialise, ngwfs_merge_sets
    use ngwf_representation, only: NGWF_REP, ngwf_rep_create, NGWF_HAM, &
         ngwf_ham_create, ngwf_rep_destroy
    use parallel_strategy, only: PARAL_INFO
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_output_detail, cond_num_states, pub_spin, pub_aug, &
         pub_debug_on_root, cond_init_shift, pub_num_spins, PUB_1K, pub_paw
    use services, only: services_flush
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_scale, sparse_embed_copy, &
         sparse_embed_axpy, sparse_embed_trace, sparse_embed_product, &
         sparse_embed_transpose_structure, sparse_embed_transpose
    use sparse, only: sparse_put_element, BLKS_JOINTPH, sparse_init_blocking_scheme
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3_EMBED_ARRAY), intent(inout) :: denskern
    type(NGWF_REP), intent(in) :: val_rep
    type(FUNC_BASIS), intent(in) :: val_basis(1)
    type(FUNC_BASIS), intent(in) :: proj_basis(1)
    type(FUNC_BASIS), intent(inout) :: hub_proj_basis(1)
    type(MODEL), target, intent(inout) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    type(NGWF_HAM), intent(inout) :: ham
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(1)
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    type(FUNC_BASIS), intent(inout) :: FOngwf_basis(1)
    type(NGWF_REP), intent(inout)   :: FOngwf_rep
    type(FUNC_BASIS), intent(inout) :: joint_basis(1)
    type(NGWF_REP), intent(inout) :: joint_rep
    type(NGWF_HAM), intent(inout) :: joint_ham
    type(DKERN), intent(inout) :: condkern_new

    ! Local Variables
    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    real(kind=DP) :: lhxc_energy
    real(kind=DP) :: hubbard_energy
    integer :: ierr
    integer :: is
    logical :: shift_changed
    real(kind=DP) :: total_energy
    type(SPAM3_EMBED) :: temp_mat_1, temp_mat_2, temp_mat_3
    integer :: local_size
    integer(kind=LONG) :: global_size

    if (pub_debug_on_root) write(stdout, '(a)') &
         'DEBUG: Entering lr_phonons_initialise routine.'

    ! tjz07: Try to create the valence hamiltonian
    call ngwf_ham_create(ham, val_rep)

    !ndmh: calculate density dependent energies and matrices for valence NGWFS
    call hamiltonian_dens_dep_matrices(ham, lhxc_fine, total_energy, &
         lhxc_energy, hubbard_energy, paw_sphere_energies, val_rep, &
         val_basis, hub_proj_basis, hub, denskern, mdl, hfxstate, &
         ham_update=.true., lhxc_fixed=.false.)

    call hamiltonian_build_matrix(ham, val_rep)

    if (pub_on_root .and. pub_output_detail>=VERBOSE) then
       write(stdout, '(a)') 'Successfully built hamiltonian matrix &
            &for valence hamiltonian.'
    endif

    ! ====================  CREATE FOngwf BASIS =============================
    call function_basis_allocate(FOngwf_basis(1),mdl%par%num_ngwfs,'ngwfs_FOngwf', &
         par=mdl%par)
    call function_basis_distribute(FOngwf_basis(1),mdl%elements, par=mdl%par)

    call comms_barrier

    call function_basis_copy_spheres(FOngwf_basis(1), val_basis(1), par=mdl%par)

    call function_basis_init_tight_boxes(FOngwf_basis(1),mdl%fftbox,mdl%cell,&
            par=mdl%par)
    if (pub_on_root) &
         write(stdout,'(a,i0)') '... FOngwf NGWF tight boxes initialised '

    call function_basis_gath_all_tbs(FOngwf_basis(1), par=mdl%par)

    call function_basis_estimate_size(FOngwf_basis(1),local_size,global_size)
    if (pub_on_root) write(stdout,'(a,2i13)') &
         '... FOngwf NGWF basis size (local,global): ', local_size, global_size

    call comms_barrier

    call ngwf_rep_create(FOngwf_rep,'',mdl,is_cmplx=val_rep%overlap%iscmplx)

    call comms_barrier

    call data_functions_alloc(FOngwf_rep%ngwfs_on_grid(1), &
         FOngwf_basis(1)%size_on_grid, iscmplx=val_rep%overlap%iscmplx)

    call comms_barrier


    if (pub_on_root) write(stdout,'(a/)') '... done'

    ! tjz07: Set the correct occupation numbers for the joint set.
    ! Using the joint set means that ALL conduction states representable
    ! by the joint representation are included in the calculation
    FOngwf_rep%n_occ(:,:) = val_rep%n_occ(:,:)

    call comms_barrier
    ! ================== END OF FOngwf BASIS CREATION =======================

    ! ====================  CREATE joint BASIS =============================
    call function_basis_allocate(joint_basis(1), 2*mdl%par%num_ngwfs,'ngwfs_jointph',&
         par=mdl%par)
    call function_basis_distribute(joint_basis(1),mdl%elements, par=mdl%par)

    call comms_barrier

    call function_basis_copy_spheres(joint_basis(1), FOngwf_basis(1), val_basis(1),&
         par=mdl%par)

    call function_basis_init_tight_boxes(joint_basis(1),mdl%fftbox,mdl%cell, &
         par=mdl%par)
    if (pub_on_root) &
         write(stdout,'(a)') '... joint NGWF tight boxes initialised'

    call function_basis_gath_all_tbs(joint_basis(1), par=mdl%par)
    call function_basis_estimate_size(joint_basis(1),local_size,global_size)
    if (pub_on_root) write(stdout,'(a,2i13)') &
         '... joint NGWF basis size (local,global): ', local_size, global_size

    call comms_barrier

    call sparse_init_blocking_scheme(BLKS_JOINTPH,joint_basis(1)%num, &
         joint_basis(1)%num_on_proc, joint_basis(1)%num_on_atom, &
         joint_basis(1)%first_on_proc, joint_basis(1)%first_on_atom, &
         joint_basis(1)%atom_of_func, joint_basis(1)%proc_of_func, &
         mdl%par)


    call ngwf_rep_create(joint_rep,'z',mdl,is_cmplx=val_rep%overlap%iscmplx)

    call comms_barrier

    call data_functions_alloc(joint_rep%ngwfs_on_grid(1), &
         joint_basis(1)%size_on_grid, iscmplx=val_rep%overlap%iscmplx)

    call comms_barrier


    if (pub_on_root) write(stdout,'(a/)') '... done'

    ! tjz07: Set the correct occupation numbers for the joint set.
    ! Using the joint set means that ALL conduction states representable
    ! by the joint representation are included in the calculation
    joint_rep%n_occ(:,:) = val_rep%n_occ(:,:)

    call comms_barrier

    call ngwf_ham_create(joint_ham, joint_rep)
    if (pub_aug) then
       do is = 1,pub_num_spins
          call sparse_embed_copy(joint_ham%dijhat(is), ham%dijhat(is))
       end do
    end if
    ! ================== END OF joint BASIS CREATION =======================


    call kernel_create(condkern_new,'K'//joint_rep%postfix, &
         is_cmplx = denskern%m(1,1)%iscmplx)

  end subroutine lr_phonons_initialise


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_phonons_solve_rkn(mdl, val_basis, val_rep, val_ham, denskern, &
        A_mat, B_mat, C_mat, D_mat, respkern_jv, auxrespkern_jv, &
        joint_basis, joint_rep, joint_ham, condkern_new, direction, &
        weight, perturbation, proj_basis, nl_projectors, dij, dens_fine, &
        FO_tcore_density, FO_lpseudo_fine, rho_ij, lhxc_fine, FO_lhxc_fine, &
        fxc_buffer, vxc_buffer, FO_ham, FO_ham_jv, FO_dij, iGpv_overlap, &
        FO_overlap, FO_overlap_jv, dipole_mat_jv, at_gamma)

    ! ========================================================================!
    ! This subroutine optimises the response kernel by solving the Sternheimer!
    ! Equation (in the form of a generalised Sylvester equation) for an       !
    ! auxiliary response kernel. The response kernel is obtained from its     !
    ! auxiliary counterpart by applying the correct gauge constraint          !
    ! ------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in 2014-2017                          !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,        !
    ! September 2018                                                          !
    ! ========================================================================!

    use datatypes, only: FUNCTIONS
    use comms, only: pub_on_root, pub_my_proc_id, comms_barrier
    use function_basis, only: FUNC_BASIS
    use constants, only: stdout
    use dense,only: DEM, dense_create, dense_destroy, dense_product, &
         dense_axpy,dense_scale,dense_get_element, dense_convert,&
         dense_copy, dense_transpose, dense_put_element, dense_norm, &
         dense_write, dense_trace, dense_invert
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_num_spins, pub_debug_on_root,pub_aug,pub_paw, &
         pub_any_nl_proj, maxit_lnv
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_axpy, sparse_embed_product, &
         sparse_embed_transpose_structure, sparse_embed_transpose, &
         sparse_embed_scale, sparse_embed_copy, sparse_embed_trace, &
         sparse_embed_conjugate, sparse_embed_array_scale
    use utils, only: utils_abort, utils_alloc_check,utils_dealloc_check, &
         utils_unit, utils_isnan
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    type(FUNC_BASIS), intent(in) :: val_basis(1)
    type(NGWF_REP), intent(in) :: val_rep
    type(NGWF_HAM), intent(in) :: val_ham
    type(SPAM3_EMBED_ARRAY), intent(inout) :: denskern
    type(SPAM3_EMBED), intent(in) :: A_mat, B_mat, C_mat, D_mat
    type(SPAM3_EMBED), intent(inout) :: respkern_jv(pub_num_spins), &
         auxrespkern_jv(pub_num_spins)
    type(FUNC_BASIS), intent(in) :: joint_basis(1)
    type(NGWF_REP), intent(in) :: joint_rep
    type(NGWF_HAM), intent(in) :: joint_ham
    type(SPAM3_EMBED_ARRAY), intent(in) :: condkern_new
    integer, intent(in) :: direction(mdl%nat)
    real(kind=DP), intent(in) :: weight(mdl%nat)
    integer, intent(in) :: perturbation
    type(FUNC_BASIS), intent(in) :: proj_basis(1)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(1)
    type(SPAM3_EMBED), intent(in) :: dij(pub_num_spins)
    real(kind=DP), intent(inout) :: dens_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in) :: FO_tcore_density(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12)
    real(kind=DP), intent(in) :: FO_lpseudo_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12)
    type(SPAM3_EMBED), intent(in) :: rho_ij(pub_num_spins)
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(inout) :: FO_lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(inout) :: fxc_buffer(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(inout) :: vxc_buffer(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: FO_ham, FO_ham_jv
    type(SPAM3_EMBED), intent(inout) :: FO_dij(pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: iGpv_overlap
    type(SPAM3_EMBED), intent(inout) :: FO_overlap, FO_overlap_jv
    type(SPAM3_EMBED), intent(in) :: dipole_mat_jv(3)
    logical, intent(in) :: at_gamma

    ! Local Variables
    real(kind=DP) :: delta_t,a_val,b_val,c_val,d_val, &
                     e_val,f_val,beta_1,beta_2,beta_3
    real(kind=DP) :: g_val,mu_val,nu_val,mixing
    real(kind=DP) :: lambda,lambda_old,lambda_0
    real(kind=DP) :: sigma_min,sigma_max, gamma_val, element
    real(kind=DP) :: alpha_plus,alpha_val,alpha_p_new, norm_z, norm_gFZ
    real(kind=DP) :: residual,norm_sqr_zz, norm_sqr_zy, residual_temp
    real(kind=DP) :: f_k,f_new,norm_temp_plus, grad_f_new, f_old

    real(kind=DP) :: dummy_gzero

    integer :: num_ngwfs, iter, iter1, iter2, iter1_eff, iter2_eff
    integer :: is,ierr,n_max_iter,n_prec_iter,n_prec_iter_0,k, piter,step_ok

    type(SPAM3_EMBED) :: XMINZ_mat,XplusK_1,XplusK_2

    real(kind=DP) :: residual_test, f_new_test

    type(SPAM3_EMBED) :: A_mat_CT, C_mat_CT, D_mat_CT, B_mat_CT
    type(SPAM3_EMBED) :: A_lin_op_temp, F_mat_temp, F_mat_temp_CT
    type(SPAM3_EMBED) :: Z_mat_CT, G_mat, F_mat_CT
    type(SPAM3_EMBED) :: E_mat_fixed, E_mat_var

    type(SPAM3_EMBED) :: temp_mat_1, temp_mat_2, temp_mat_3, temp_mat_4, temp_mat_5, &
         temp_mat_6, temp_mat_7, temp_mat_8, temp_mat_9

    type(SPAM3_EMBED) :: respkern_jv_old(pub_num_spins), &
         auxrespkern_jv_old(pub_num_spins), respkern_jv_CT(pub_num_spins), &
         auxrespkern_jv_CT(pub_num_spins)
    type(SPAM3_EMBED) :: F_mat, A_lin_op, Z_mat, Y_mat_diff, Y_mat, F_mat_old
    type(SPAM3_EMBED) :: K_mat(3)
    type(SPAM3_EMBED) :: Af_lin_op, Af_lin_op_sq, Af_lin_op_cb
    type(SPAM3_EMBED) :: Af_lin_op_CT, Af_lin_op_sq_CT, Af_lin_op_cb_CT
    type(SPAM3_EMBED) :: grad_F_mat, grad_F_mat_old, grad_F_mat_CT

    type(SPAM3_EMBED) :: rho_ij_dproj(pub_num_spins)
    type(SPAM3_EMBED) :: viGp_overlap
    type(SPAM3_EMBED) :: FO_ham_jv_fixed, FO_ham_jv_fixed_temp
    type(SPAM3_EMBED) :: dijps, pv_overlap
    type(SPAM3_EMBED) :: pj_overlap, iGpj_overlap, jiGp_overlap

    if (pub_debug_on_root) write(stdout, '(a)') &
         'DEBUG: Entering lr_phonons_solve_rkn routine.'

    ! Start timer
    call timer_clock('lr_phonons_solve_rkn',1)

    num_ngwfs = val_basis(1)%num

    if (pub_paw.or.pub_any_nl_proj) then
       call sparse_embed_create(viGp_overlap, val_rep%sp_overlap)
       call sparse_embed_transpose_structure(pj_overlap%structure, &
            joint_rep%sp_overlap)
       call sparse_embed_create(pj_overlap)
       call sparse_embed_transpose(pj_overlap, joint_rep%sp_overlap)
       call sparse_embed_create(iGpj_overlap, pj_overlap)
       call sparse_embed_create(jiGp_overlap, joint_rep%sp_overlap)
    end if

    if (pub_paw) then
       do is=1,pub_num_spins
          call sparse_embed_create(rho_ij_dproj(is), rho_ij(is))
          call sparse_embed_scale(rho_ij_dproj(is), 0.0_DP)
       end do
    end if

    call sparse_embed_create(FO_ham_jv_fixed, FO_ham_jv)

    if (pub_paw.or.pub_any_nl_proj) then
       call sparse_embed_create(FO_ham_jv_fixed_temp, FO_ham_jv)
    end if

    ! create and initialise sparse matrices

    call sparse_embed_create(temp_mat_1, joint_rep%overlap, condkern_new%m(1,1))
    call sparse_embed_create(temp_mat_2, temp_mat_1, joint_ham%ham(1))
    call sparse_embed_create(A_mat_CT, temp_mat_2, condkern_new%m(1,1))
    call sparse_embed_create(C_mat_CT, A_mat_CT)
    call sparse_embed_destroy(temp_mat_1)
    call sparse_embed_destroy(temp_mat_2)

    call sparse_embed_create(temp_mat_1, val_rep%inv_overlap, val_ham%ham(1))
    call sparse_embed_create(temp_mat_2, temp_mat_1, denskern%m(1,1))
    call sparse_embed_create(D_mat_CT, temp_mat_2, val_rep%overlap)
    call sparse_embed_create(B_mat_CT, D_mat_CT)
    call sparse_embed_destroy(temp_mat_1)
    call sparse_embed_destroy(temp_mat_2)

    call sparse_embed_transpose(A_mat_CT, A_mat)
    if (respkern_jv(1)%iscmplx) then
       call sparse_embed_conjugate(A_mat_CT)
    end if

    call sparse_embed_transpose(B_mat_CT, B_mat)
    if (respkern_jv(1)%iscmplx) then
       call sparse_embed_conjugate(B_mat_CT)
    end if

    call sparse_embed_transpose(C_mat_CT, C_mat)
    if (respkern_jv(1)%iscmplx) then
       call sparse_embed_conjugate(C_mat_CT)
    end if

    call sparse_embed_transpose(D_mat_CT, D_mat)
    if (respkern_jv(1)%iscmplx) then
       call sparse_embed_conjugate(D_mat_CT)
    end if

    do is = 1, pub_num_spins

       call sparse_embed_transpose_structure(respkern_jv_CT(is)%structure, &
            respkern_jv(is))
       call sparse_embed_create(respkern_jv_CT(is))

       call sparse_embed_transpose_structure(auxrespkern_jv_CT(is)%structure, &
            auxrespkern_jv(is))
       call sparse_embed_create(auxrespkern_jv_CT(is))

       call sparse_embed_create(respkern_jv_old(is), respkern_jv(is))
       call sparse_embed_create(auxrespkern_jv_old(is), auxrespkern_jv(is))

    end do

    call sparse_embed_create(temp_mat_1, A_mat, respkern_jv(1))
    call sparse_embed_create(F_mat, temp_mat_1, B_mat)
    call sparse_embed_destroy(temp_mat_1)

    call sparse_embed_create(temp_mat_1, B_mat_CT, respkern_jv_CT(1))
    call sparse_embed_create(F_mat_CT, temp_mat_1, A_mat_CT)
    call sparse_embed_destroy(temp_mat_1)

    call sparse_embed_create(F_mat_old, F_mat)
    call sparse_embed_create(F_mat_temp, F_mat)
    call sparse_embed_create(F_mat_temp_CT, F_mat_CT)

    call sparse_embed_create(E_mat_var, F_mat)
    call sparse_embed_create(E_mat_fixed, F_mat)

    call sparse_embed_create(Z_mat, auxrespkern_jv(1))
    call sparse_embed_create(Z_mat_CT, auxrespkern_jv_CT(1))

    call sparse_embed_create(XMINZ_mat, auxrespkern_jv(1))

    call sparse_embed_create(G_mat, F_mat)
    call sparse_embed_create(Y_mat, F_mat)
    call sparse_embed_create(Y_mat_diff, F_mat)

    call sparse_embed_create(K_mat(1), F_mat)
    call sparse_embed_create(K_mat(2), F_mat)
    call sparse_embed_create(K_mat(3), F_mat)

    call sparse_embed_create(XplusK_1, F_mat)
    call sparse_embed_create(XplusK_2, F_mat)

    call sparse_embed_create(A_lin_op, F_mat)
    call sparse_embed_create(A_lin_op_temp, F_mat)

    call sparse_embed_create(temp_mat_1, A_mat, F_mat)
    call sparse_embed_create(Af_lin_op, temp_mat_1, B_mat)
    call sparse_embed_destroy(temp_mat_1)

    call sparse_embed_create(temp_mat_1, B_mat_CT, F_mat_CT)
    call sparse_embed_create(Af_lin_op_CT, temp_mat_1, A_mat_CT)
    call sparse_embed_destroy(temp_mat_1)

    call sparse_embed_create(Af_lin_op_sq, Af_lin_op)
    call sparse_embed_create(Af_lin_op_cb, Af_lin_op)
    call sparse_embed_create(Af_lin_op_sq_CT, Af_lin_op_CT)
    call sparse_embed_create(Af_lin_op_cb_CT, Af_lin_op_CT)

    call sparse_embed_create(temp_mat_1, A_mat_CT, joint_rep%overlap)
    call sparse_embed_create(temp_mat_2, temp_mat_1, F_mat)
    call sparse_embed_destroy(temp_mat_1)
    call sparse_embed_create(temp_mat_1, temp_mat_2, val_rep%overlap)
    call sparse_embed_create(grad_F_mat, temp_mat_1, B_mat_CT)
    call sparse_embed_destroy(temp_mat_1)
    call sparse_embed_destroy(temp_mat_2)

    call sparse_embed_create(grad_F_mat_old, grad_F_mat)

    call sparse_embed_create(temp_mat_1, B_mat, val_rep%overlap)
    call sparse_embed_create(temp_mat_2, temp_mat_1, F_mat_CT)
    call sparse_embed_destroy(temp_mat_1)
    call sparse_embed_create(temp_mat_1, temp_mat_2, joint_rep%overlap)
    call sparse_embed_create(grad_F_mat_CT, temp_mat_1, A_mat)
    call sparse_embed_destroy(temp_mat_1)
    call sparse_embed_destroy(temp_mat_2)

    do is = 1, pub_num_spins
       call sparse_embed_destroy(respkern_jv_CT(is))
       call sparse_embed_destroy(auxrespkern_jv_CT(is))
    end do

    n_max_iter = maxit_lnv
    n_prec_iter_0 = 20 ! number of preconditioning time-marching steps

    mixing = 1.0_DP

    gamma_val = 0.01_DP
    lambda_0=1.0_DP  ! default step length

    sigma_min = 1.0D-10 ! min cap for step length
    sigma_max = 1.0D+10 ! max cap for step length

    !! BEGIN BUILD OF E_mat_fixed
    !===========================================================
    if (pub_paw.or.pub_any_nl_proj) then
       call sparse_embed_scale(iGpv_overlap, 0.0_DP)
       call sparse_embed_scale(jiGp_overlap, 0.0_DP)
    end if

    call sparse_embed_scale(FO_overlap, 0.0_DP)
    call sparse_embed_scale(FO_overlap_jv, 0.0_DP)

    if (pub_paw.or.pub_any_nl_proj) then
       if (perturbation.le.(3*mdl%nat)) then ! atomic displacements only
          call calculate_first_overlap(mdl, val_rep%ngwfs_on_grid, &
               val_rep%sp_overlap, val_basis, joint_rep%ngwfs_on_grid, &
               joint_rep%sp_overlap, joint_basis, proj_basis, &
               nl_projectors, FO_overlap_jv, FO_overlap, iGpv_overlap, &
               jiGp_overlap, direction, weight, rho_ij_dproj, denskern%m(:,1))
       end if

       call sparse_embed_transpose(iGpj_overlap, jiGp_overlap)
       call sparse_embed_transpose(viGp_overlap, iGpv_overlap)
    end if

    call sparse_embed_create(temp_mat_1, FO_overlap_jv, denskern%m(1,1))
    call sparse_embed_product(temp_mat_1, FO_overlap_jv, denskern%m(1,1))
    call sparse_embed_create(temp_mat_2, temp_mat_1, val_ham%ham(1))
    call sparse_embed_product(temp_mat_2, temp_mat_1, val_ham%ham(1))
    call sparse_embed_destroy(temp_mat_1)
    call sparse_embed_create(temp_mat_1, temp_mat_2, val_rep%inv_overlap)
    call sparse_embed_product(temp_mat_1, temp_mat_2, val_rep%inv_overlap)
    call sparse_embed_product(E_mat_fixed, condkern_new%m(1,1), temp_mat_1)
    call sparse_embed_destroy(temp_mat_1)
    call sparse_embed_destroy(temp_mat_2)

    ! for E-field perturbations at gamma only
    if (at_gamma.and.(perturbation.gt.(3*mdl%nat))) then
       call sparse_embed_create(temp_mat_1, condkern_new%m(1,1), &
            dipole_mat_jv(mod(perturbation,(3*mdl%nat))) )
       call sparse_embed_product(temp_mat_1, condkern_new%m(1,1), &
            dipole_mat_jv(mod(perturbation,(3*mdl%nat))) )
       call sparse_embed_create(temp_mat_2, temp_mat_1, denskern%m(1,1))
       call sparse_embed_product(temp_mat_2, temp_mat_1, denskern%m(1,1))
       call sparse_embed_axpy(E_mat_fixed, temp_mat_2, -1.0_DP)

       call sparse_embed_destroy(temp_mat_1)
       call sparse_embed_destroy(temp_mat_2)
    end if
    !===========================================================
    !! END BUILD OF E_mat_fixed

    !! BUILD FIXED PART OF FO HAMILTONIAN ================================
    ! calculate this local component for the FO hamiltonian
    call sparse_embed_scale(FO_ham_jv_fixed, 0.0_DP)

    if (pub_paw.or.pub_any_nl_proj) then
       call sparse_embed_create(dijps, dij(1), iGpv_overlap)
       call sparse_embed_product(dijps, dij(1), iGpv_overlap)
       call sparse_embed_product(FO_ham_jv_fixed, joint_rep%sp_overlap, dijps)

       call sparse_embed_transpose_structure(pv_overlap%structure, val_rep%sp_overlap)
       call sparse_embed_create(pv_overlap)
       call sparse_embed_transpose(pv_overlap, val_rep%sp_overlap)

       call sparse_embed_scale(dijps, 0.0_DP)
       call sparse_embed_product(dijps, dij(1), pv_overlap)
       call sparse_embed_product(FO_ham_jv_fixed_temp, jiGp_overlap, dijps)
       call sparse_embed_destroy(dijps)

       call sparse_embed_axpy(FO_ham_jv_fixed, FO_ham_jv_fixed_temp, 1.0_DP)
       call sparse_embed_destroy(FO_ham_jv_fixed_temp)
    end if

    !! END BUILD FIXED PART OF FO HAMILTONIAN===========================

    !  call comms_barrier

    do k = 1,n_max_iter  ! optimisation loop

       call sparse_embed_scale(FO_ham_jv, 0.0_DP)
       ! call comms_barrier

       !! BUILD VARIABLE FO HAMILTONIAN
       call calc_FO_ham(mdl, val_rep, val_basis, direction, weight, dij, &
            FO_dij, dens_fine, denskern, rho_ij_dproj, rho_ij, respkern_jv, &
            FO_lpseudo_fine, lhxc_fine, FO_lhxc_fine, fxc_buffer, FO_ham_jv, &
            joint_rep, joint_basis, proj_basis, nl_projectors, &
            FO_tcore_density)

       ! add fixed part of FO_ham_jv
       call sparse_embed_axpy(FO_ham_jv, FO_ham_jv_fixed, 1.0_DP)

       ! get FO_ham from FO_ham_jv
       call sparse_embed_scale(FO_ham, 0.0_DP)
       call lr_phonons_get_basis_block(FO_ham_jv%p, FO_ham%p, joint_basis(1), &
            val_basis(1), val_basis(1), val_basis(1), 2, 1)
       !! END OF BUILD VARIABLE FO HAMILTONIAN

       ! call comms_barrier

       !! BEGIN BUILD OF E_mat_var
       call sparse_embed_create(temp_mat_1, condkern_new%m(1,1), FO_ham_jv)
       call sparse_embed_product(temp_mat_1, condkern_new%m(1,1), FO_ham_jv)
       call sparse_embed_product(E_mat_var, temp_mat_1, denskern%m(1,1))
       call sparse_embed_scale(E_mat_var, -1.0_DP)
       call sparse_embed_destroy(temp_mat_1)
       !! END BUILD OF E_mat_var

       ! call comms_barrier

       !A_lin_op
       call sparse_embed_create(temp_mat_1, A_mat, auxrespkern_jv(1))
       call sparse_embed_product(temp_mat_1, A_mat, auxrespkern_jv(1))
       call sparse_embed_create(temp_mat_2, temp_mat_1, B_mat)
       call sparse_embed_product(temp_mat_2, temp_mat_1, B_mat)
       call sparse_embed_copy(A_lin_op, temp_mat_2)
       call sparse_embed_destroy(temp_mat_1)
       call sparse_embed_destroy(temp_mat_2)

       call sparse_embed_create(temp_mat_1, C_mat, auxrespkern_jv(1))
       call sparse_embed_product(temp_mat_1, C_mat, auxrespkern_jv(1))
       call sparse_embed_create(temp_mat_2, temp_mat_1, D_mat)
       call sparse_embed_product(temp_mat_2, temp_mat_1, D_mat)
       call sparse_embed_axpy(A_lin_op, temp_mat_2, 1.0_DP)
       call sparse_embed_destroy(temp_mat_1)
       call sparse_embed_destroy(temp_mat_2)
       ! end of A_lin_op

       call sparse_embed_copy(F_mat, A_lin_op)
       call sparse_embed_axpy(F_mat, E_mat_fixed, -1.0_DP)
       call sparse_embed_axpy(F_mat, E_mat_var, -1.0_DP)

       call sparse_embed_transpose(F_mat_CT, F_mat)
       if (respkern_jv(1)%iscmplx) then
          call sparse_embed_conjugate(F_mat_CT)
       end if

       ! call comms_barrier

       ! build grad F_mat
       call sparse_embed_create(temp_mat_1, A_mat_CT, joint_rep%overlap)
       call sparse_embed_product(temp_mat_1, A_mat_CT, joint_rep%overlap)
       call sparse_embed_create(temp_mat_2, temp_mat_1, F_mat)
       call sparse_embed_product(temp_mat_2, temp_mat_1, F_mat)
       call sparse_embed_destroy(temp_mat_1)
       call sparse_embed_create(temp_mat_1, temp_mat_2, val_rep%overlap)
       call sparse_embed_product(temp_mat_1, temp_mat_2, val_rep%overlap)
       call sparse_embed_destroy(temp_mat_2)
       call sparse_embed_create(temp_mat_2, temp_mat_1, B_mat_CT)
       call sparse_embed_product(temp_mat_2, temp_mat_1, B_mat_CT)
       call sparse_embed_copy(grad_F_mat, temp_mat_2)
       call sparse_embed_destroy(temp_mat_1)
       call sparse_embed_destroy(temp_mat_2)

       call sparse_embed_create(temp_mat_1, C_mat_CT, joint_rep%overlap)
       call sparse_embed_product(temp_mat_1, C_mat_CT, joint_rep%overlap)
       call sparse_embed_create(temp_mat_2, temp_mat_1, F_mat)
       call sparse_embed_product(temp_mat_2, temp_mat_1, F_mat)
       call sparse_embed_destroy(temp_mat_1)
       call sparse_embed_create(temp_mat_1, temp_mat_2, val_rep%overlap)
       call sparse_embed_product(temp_mat_1, temp_mat_2, val_rep%overlap)
       call sparse_embed_destroy(temp_mat_2)
       call sparse_embed_create(temp_mat_2, temp_mat_1, D_mat_CT)
       call sparse_embed_product(temp_mat_2, temp_mat_1, D_mat_CT)
       call sparse_embed_axpy(grad_F_mat, temp_mat_2, 1.0_DP)
       call sparse_embed_destroy(temp_mat_1)
       call sparse_embed_destroy(temp_mat_2)

       call sparse_embed_create(temp_mat_1, grad_F_mat)
       call sparse_embed_copy(temp_mat_1, grad_F_mat)
       if (respkern_jv(1)%iscmplx) then
          call sparse_embed_conjugate(temp_mat_1)
       end if

       call sparse_embed_axpy(grad_F_mat, temp_mat_1, 1.0_DP)
       call sparse_embed_destroy(temp_mat_1)

       call sparse_embed_transpose(grad_F_mat_CT, grad_F_mat)
       if (respkern_jv(1)%iscmplx) then
          call sparse_embed_conjugate(grad_F_mat_CT)
       end if
       ! end build of grad F mat

       call sparse_embed_create(temp_mat_1, joint_rep%overlap, F_mat)
       call sparse_embed_product(temp_mat_1, joint_rep%overlap, F_mat)
       call sparse_embed_trace(f_new, F_mat_CT, temp_mat_1)
       call sparse_embed_destroy(temp_mat_1)
       residual = sqrt(f_new)

       if(pub_on_root) then
          write(stdout,'(a,i4.3,a,f20.10)') '   ITERATION ',k,'  | RESIDUAL ', &
               residual
       end if

       if (k > 1) then ! calculate lambda, if not the first iteration

          norm_sqr_zz = 0.0_DP
          norm_sqr_zy = 0.0_DP

          call sparse_embed_transpose(Z_mat_CT, Z_mat)
          if (respkern_jv(1)%iscmplx) then
             call sparse_embed_conjugate(Z_mat_CT)
          end if
          call sparse_embed_create(temp_mat_1, joint_rep%overlap, F_mat_old)
          call sparse_embed_product(temp_mat_1, joint_rep%overlap, F_mat_old)
          call sparse_embed_trace(norm_sqr_zz, Z_mat_CT, temp_mat_1)
          call sparse_embed_destroy(temp_mat_1)

          call sparse_embed_copy(Y_mat_diff, F_mat)
          call sparse_embed_axpy(Y_mat_diff, F_mat_old, -1.0_DP)

          call sparse_embed_create(temp_mat_1, joint_rep%overlap, Y_mat_diff)
          call sparse_embed_product(temp_mat_1, joint_rep%overlap, Y_mat_diff)
          call sparse_embed_trace(norm_sqr_zy, Z_mat_CT, temp_mat_1)
          call sparse_embed_destroy(temp_mat_1)

          lambda_old = lambda
          lambda = lambda_old * abs(norm_sqr_zz / norm_sqr_zy)

          if (lambda.lt.sigma_min) then
             if (pub_on_root) then
                write(stdout,*) '|lambda| too small. fixing'
             end if
                lambda=sigma_min
          elseif (lambda.gt.sigma_max) then
             if (pub_on_root) then
                write(stdout,*) '|lambda| too large. fixing'
             end if
                lambda = sigma_max
          endif
       end if
       ! call comms_barrier


       if (mod(k, 10) == 1) then
           lambda = lambda_0
       end if

       if (utils_isnan(lambda)) exit

       ! secondly, employ the Enhanced Cauchy 2 scheme to get Z_mat
       if (residual < 0.1_DP) then
          n_prec_iter = n_prec_iter_0 * ceiling(-log10(residual) + 1.0_DP)
       else
          n_prec_iter = n_prec_iter_0
       end if

       ! call comms_barrier

       if(mod(k,10) == 1) then
          call sparse_embed_scale(Y_mat, 0.0_DP)
       else
          call sparse_embed_copy(Y_mat, Z_mat)
       end if


       do piter = 1, n_prec_iter

          call sparse_embed_copy(XMINZ_mat, auxrespkern_jv(1))
          call sparse_embed_axpy(XMINZ_mat, Y_mat, -1.0_DP)

          call sparse_embed_create(temp_mat_1, A_mat, XMINZ_mat)
          call sparse_embed_product(temp_mat_1, A_mat, XMINZ_mat)
          call sparse_embed_create(temp_mat_2, temp_mat_1, B_mat)
          call sparse_embed_product(temp_mat_2, temp_mat_1, B_mat)
          call sparse_embed_copy(G_mat, temp_mat_2)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)

          call sparse_embed_create(temp_mat_1, C_mat, XMINZ_mat)
          call sparse_embed_product(temp_mat_1, C_mat, XMINZ_mat)
          call sparse_embed_create(temp_mat_2, temp_mat_1, D_mat)
          call sparse_embed_product(temp_mat_2, temp_mat_1, D_mat)
          call sparse_embed_axpy(G_mat, temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)

          call sparse_embed_axpy(G_mat, E_mat_fixed, -1.0_DP)
          call sparse_embed_axpy(G_mat, E_mat_var, -1.0_DP)

          ! Af_lin operators and powers
          call sparse_embed_create(temp_mat_1, A_mat, G_mat)
          call sparse_embed_product(temp_mat_1, A_mat, G_mat)
          call sparse_embed_create(temp_mat_2, temp_mat_1, B_mat)
          call sparse_embed_product(temp_mat_2, temp_mat_1, B_mat)
          call sparse_embed_copy(Af_lin_op, temp_mat_2)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)

          call sparse_embed_create(temp_mat_1, C_mat, G_mat)
          call sparse_embed_product(temp_mat_1, C_mat, G_mat)
          call sparse_embed_create(temp_mat_2, temp_mat_1, D_mat)
          call sparse_embed_product(temp_mat_2, temp_mat_1, D_mat)
          call sparse_embed_axpy(Af_lin_op, temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)


          call sparse_embed_create(temp_mat_1, A_mat, Af_lin_op)
          call sparse_embed_product(temp_mat_1, A_mat, Af_lin_op)
          call sparse_embed_create(temp_mat_2, temp_mat_1, B_mat)
          call sparse_embed_product(temp_mat_2, temp_mat_1, B_mat)
          call sparse_embed_copy(Af_lin_op_sq, temp_mat_2)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)

          call sparse_embed_create(temp_mat_1, C_mat, Af_lin_op)
          call sparse_embed_product(temp_mat_1, C_mat, Af_lin_op)
          call sparse_embed_create(temp_mat_2, temp_mat_1, D_mat)
          call sparse_embed_product(temp_mat_2, temp_mat_1, D_mat)
          call sparse_embed_axpy(Af_lin_op_sq, temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)


          call sparse_embed_create(temp_mat_1, A_mat, Af_lin_op_sq)
          call sparse_embed_product(temp_mat_1, A_mat, Af_lin_op_sq)
          call sparse_embed_create(temp_mat_2, temp_mat_1, B_mat)
          call sparse_embed_product(temp_mat_2, temp_mat_1, B_mat)
          call sparse_embed_copy(Af_lin_op_cb, temp_mat_2)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)

          call sparse_embed_create(temp_mat_1, C_mat, Af_lin_op_sq)
          call sparse_embed_product(temp_mat_1, C_mat, Af_lin_op_sq)
          call sparse_embed_create(temp_mat_2, temp_mat_1, D_mat)
          call sparse_embed_product(temp_mat_2, temp_mat_1, D_mat)
          call sparse_embed_axpy(Af_lin_op_cb, temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
          ! end of Af_lin operators and powers

          call sparse_embed_transpose(Af_lin_op_CT, Af_lin_op)
          call sparse_embed_transpose(Af_lin_op_sq_CT, Af_lin_op_sq)
          call sparse_embed_transpose(Af_lin_op_cb_CT, Af_lin_op_cb)
          if (respkern_jv(1)%iscmplx) then
             call sparse_embed_conjugate(Af_lin_op_CT)
             call sparse_embed_conjugate(Af_lin_op_sq_CT)
             call sparse_embed_conjugate(Af_lin_op_cb_CT)
          end if


          a_val = 0.0_DP
          call sparse_embed_create(temp_mat_1, joint_rep%overlap, Af_lin_op)
          call sparse_embed_product(temp_mat_1, joint_rep%overlap, Af_lin_op)
          call sparse_embed_trace(a_val, Af_lin_op_CT, temp_mat_1)
          call sparse_embed_destroy(temp_mat_1)

          b_val = 0.0_DP
          call sparse_embed_create(temp_mat_1, joint_rep%overlap, Af_lin_op_sq)
          call sparse_embed_product(temp_mat_1, joint_rep%overlap, Af_lin_op_sq)
          call sparse_embed_trace(b_val, Af_lin_op_sq_CT, temp_mat_1)
          call sparse_embed_destroy(temp_mat_1)

          c_val = 0.0_DP
          call sparse_embed_create(temp_mat_1, joint_rep%overlap, Af_lin_op_cb)
          call sparse_embed_product(temp_mat_1, joint_rep%overlap, Af_lin_op_cb)
          call sparse_embed_trace(c_val, Af_lin_op_cb_CT, temp_mat_1)
          call sparse_embed_destroy(temp_mat_1)

          d_val = 0.0_DP
          call sparse_embed_create(temp_mat_1, joint_rep%overlap, G_mat)
          call sparse_embed_product(temp_mat_1, joint_rep%overlap, G_mat)
          call sparse_embed_trace(d_val, Af_lin_op_CT, temp_mat_1)

          e_val = 0.0_DP
          call sparse_embed_trace(e_val, Af_lin_op_sq_CT, temp_mat_1)

          f_val = 0.0_DP
          call sparse_embed_trace(f_val, Af_lin_op_cb_CT, temp_mat_1)

          call sparse_embed_destroy(temp_mat_1)

          g_val = 0.0_DP
          call sparse_embed_create(temp_mat_1, joint_rep%overlap, Af_lin_op)
          call sparse_embed_product(temp_mat_1, joint_rep%overlap, Af_lin_op)
          call sparse_embed_trace(g_val, Af_lin_op_sq_CT, temp_mat_1)

          mu_val = 0.0_DP
          call sparse_embed_trace(mu_val, Af_lin_op_cb_CT, temp_mat_1)

          call sparse_embed_destroy(temp_mat_1)

          nu_val = 0.0_DP
          call sparse_embed_create(temp_mat_1, joint_rep%overlap, Af_lin_op_sq)
          call sparse_embed_product(temp_mat_1, joint_rep%overlap, Af_lin_op_sq)
          call sparse_embed_trace(nu_val, Af_lin_op_cb_CT, temp_mat_1)
          call sparse_embed_destroy(temp_mat_1)

        delta_t = (-mu_val*nu_val*e_val-g_val*nu_val*f_val+mu_val*f_val*b_val+ &
          d_val*nu_val*nu_val-d_val*c_val*b_val+g_val*c_val*e_val ) / &
          (g_val*g_val*c_val+mu_val*mu_val*b_val-a_val*c_val*b_val+a_val* &
          nu_val*nu_val-2.0_DP*mu_val*nu_val*g_val)

        beta_1 = (nu_val*f_val-nu_val*delta_t*mu_val+delta_t*delta_t*nu_val* &
          nu_val-delta_t*delta_t*b_val*c_val-e_val*c_val+delta_t*g_val*c_val)/ &
          (-delta_t*b_val*c_val+delta_t*nu_val*nu_val)

        beta_2 = -(delta_t*nu_val*f_val-mu_val*delta_t*delta_t*nu_val-delta_t* &
          e_val*c_val+delta_t*delta_t*g_val*c_val-f_val*b_val+delta_t*mu_val* &
          b_val+nu_val*e_val-nu_val*delta_t*g_val) / (-delta_t* &
          delta_t*b_val*c_val+delta_t*delta_t*nu_val*nu_val)

          beta_3 = delta_t-beta_1-beta_2

          if (utils_isnan(delta_t).or.utils_isnan(beta_1).or. &
               utils_isnan(beta_2).or.utils_isnan(beta_3)) exit


          ! part 1
           call sparse_embed_copy(K_mat(1), G_mat)
          ! end of part 1

          ! call comms_barrier

          ! part 2
          call sparse_embed_copy(XplusK_1, XMINZ_mat)
          call sparse_embed_axpy(XplusK_1, K_mat(1), -delta_t)
          call sparse_embed_create(temp_mat_1, A_mat, XplusK_1)
          call sparse_embed_product(temp_mat_1, A_mat, XplusK_1)
          call sparse_embed_create(temp_mat_2, temp_mat_1, B_mat)
          call sparse_embed_product(temp_mat_2, temp_mat_1, B_mat)
          call sparse_embed_copy(K_mat(2), temp_mat_2)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)

          call sparse_embed_create(temp_mat_1, C_mat, XplusK_1)
          call sparse_embed_product(temp_mat_1, C_mat, XplusK_1)
          call sparse_embed_create(temp_mat_2, temp_mat_1, D_mat)
          call sparse_embed_product(temp_mat_2, temp_mat_1, D_mat)
          call sparse_embed_axpy(K_mat(2), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)

          call sparse_embed_axpy(K_mat(2), E_mat_var, -1.0_DP)
          call sparse_embed_axpy(K_mat(2), E_mat_fixed, -1.0_DP)
          ! end of part 2

          ! call comms_barrier

          ! part 3
          call sparse_embed_copy(XplusK_2, XMINZ_mat)
          call sparse_embed_axpy(XplusK_2, K_mat(2), -delta_t)
          call sparse_embed_create(temp_mat_1, A_mat, XplusK_2)
          call sparse_embed_product(temp_mat_1, A_mat, XplusK_2)
          call sparse_embed_create(temp_mat_2, temp_mat_1, B_mat)
          call sparse_embed_product(temp_mat_2, temp_mat_1, B_mat)
          call sparse_embed_copy(K_mat(3), temp_mat_2)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)

          call sparse_embed_create(temp_mat_1, C_mat, XplusK_2)
          call sparse_embed_product(temp_mat_1, C_mat, XplusK_2)
          call sparse_embed_create(temp_mat_2, temp_mat_1, D_mat)
          call sparse_embed_product(temp_mat_2, temp_mat_1, D_mat)
          call sparse_embed_axpy(K_mat(3), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)

          call sparse_embed_axpy(K_mat(3), E_mat_var, -1.0_DP)
          call sparse_embed_axpy(K_mat(3), E_mat_fixed, -1.0_DP)
          ! end of part 3

          ! call comms_barrier

          call sparse_embed_axpy(Y_mat, K_mat(1), beta_1)
          call sparse_embed_axpy(Y_mat, K_mat(2), beta_2)
          call sparse_embed_axpy(Y_mat, K_mat(3), beta_3)

       end do

      call sparse_embed_copy(Z_mat, Y_mat)

      ! NOW FOR THE GLOBALISATION SCHEME
      f_k = f_new

      call sparse_embed_scale(Z_mat, -1.0_DP * lambda)

      norm_gFZ = 0.0_DP
      call sparse_embed_create(temp_mat_1, joint_rep%overlap, Z_mat)
      call sparse_embed_product(temp_mat_1, joint_rep%overlap, Z_mat)
      call sparse_embed_trace(norm_gFZ, grad_F_mat_CT, temp_mat_1)
      call sparse_embed_destroy(temp_mat_1)

      alpha_plus = mixing

      step_ok = 0

      do while (step_ok .eq. 0)

         !F_mat_temp
         call sparse_embed_create(temp_mat_3, auxrespkern_jv(1))
         call sparse_embed_copy(temp_mat_3, auxrespkern_jv(1))
         call sparse_embed_axpy(temp_mat_3, Z_mat, alpha_plus)

         ! call comms_barrier

         call sparse_embed_create(temp_mat_1, A_mat, temp_mat_3)
         call sparse_embed_product(temp_mat_1, A_mat, temp_mat_3)
         call sparse_embed_create(temp_mat_2, temp_mat_1, B_mat)
         call sparse_embed_product(temp_mat_2, temp_mat_1, B_mat)
         call sparse_embed_copy(A_lin_op_temp, temp_mat_2)
         call sparse_embed_destroy(temp_mat_1)
         call sparse_embed_destroy(temp_mat_2)

         call sparse_embed_create(temp_mat_1, C_mat, temp_mat_3)
         call sparse_embed_product(temp_mat_1, C_mat, temp_mat_3)
         call sparse_embed_create(temp_mat_2, temp_mat_1, D_mat)
         call sparse_embed_product(temp_mat_2, temp_mat_1, D_mat)
         call sparse_embed_axpy(A_lin_op_temp, temp_mat_2, 1.0_DP)
         call sparse_embed_destroy(temp_mat_1)
         call sparse_embed_destroy(temp_mat_2)

         call sparse_embed_copy(F_mat_temp, A_lin_op_temp)
         call sparse_embed_axpy(F_mat_temp, E_mat_fixed, -1.0_DP)
         call sparse_embed_axpy(F_mat_temp, E_mat_var, -1.0_DP)

         call sparse_embed_destroy(temp_mat_3)
         ! end of F_mat_temp

         call sparse_embed_transpose(F_mat_temp_CT, F_mat_temp)
         if (respkern_jv(1)%iscmplx) then
            call sparse_embed_conjugate(F_mat_temp_CT)
         end if

         norm_temp_plus = 0.0_DP
         call sparse_embed_create(temp_mat_1, joint_rep%overlap, F_mat_temp)
         call sparse_embed_product(temp_mat_1, joint_rep%overlap, F_mat_temp)
         call sparse_embed_trace(norm_temp_plus, F_mat_temp_CT, temp_mat_1)
         call sparse_embed_destroy(temp_mat_1)

         if (utils_isnan(norm_temp_plus)) exit

         if (norm_temp_plus .le. (f_k + gamma_val * alpha_plus * norm_gFZ)) then
            call sparse_embed_copy(auxrespkern_jv_old(1), auxrespkern_jv(1))
            call sparse_embed_copy(F_mat_old, F_mat)
            call sparse_embed_copy(grad_F_mat_old, grad_F_mat)
            alpha_val = max(alpha_plus, 0.20_DP) ! gcc32 TEST

            call sparse_embed_axpy(auxrespkern_jv(1), Z_mat, alpha_val)
            call sparse_embed_scale(Z_mat, -1.0_DP / lambda)
            lambda = lambda * alpha_val
            step_ok = 1
            if (pub_on_root) then
               write(stdout,'(a,f20.10)') '      Step size: ', lambda
            end if
         else

           if (pub_on_root) then
              ! write(stdout,*) norm_temp_plus, " > ", f_k, " + ", &
              !     gamma_val * alpha_plus * norm_gFZ
           end if

           alpha_p_new = -0.5_DP * alpha_plus * alpha_plus * norm_gFZ / &
                 (norm_temp_plus - f_k - alpha_plus * norm_gFZ)

           if((alpha_p_new < 0.7_DP * alpha_plus).and. &
              (alpha_p_new > 0.1_DP * alpha_plus)) then
              alpha_plus = alpha_p_new
              if (pub_on_root) then
                 !  write(stdout,*) 'QUADRATIC INTERPOLATION'
              end if
           else
              alpha_p_new = alpha_plus / 2.0_DP
              alpha_plus = alpha_p_new
              if (pub_on_root) then
                 ! write(stdout,*) 'NOT OK'
              end if
           end if

         end if

      end do ! end do while


      !=======================================================
      ! Obtain respkern from auxrespkern by applying gauge constraint
      do is = 1, pub_num_spins
         call sparse_embed_create(temp_mat_1, condkern_new%m(is,1), &
              joint_rep%overlap)
         call sparse_embed_product(temp_mat_1, condkern_new%m(is,1), &
              joint_rep%overlap)
         call sparse_embed_create(temp_mat_2, temp_mat_1, auxrespkern_jv(is))
         call sparse_embed_product(temp_mat_2, temp_mat_1, auxrespkern_jv(is))
         call sparse_embed_destroy(temp_mat_1)
         call sparse_embed_create(temp_mat_1, temp_mat_2, val_rep%overlap)
         call sparse_embed_product(temp_mat_1, temp_mat_2, val_rep%overlap)
         call sparse_embed_destroy(temp_mat_2)
         call sparse_embed_create(temp_mat_2, temp_mat_1, denskern%m(is,1))
         call sparse_embed_product(temp_mat_2, temp_mat_1, denskern%m(is,1))
         call sparse_embed_destroy(temp_mat_1)

         call sparse_embed_transpose_structure(temp_mat_3%structure, &
              joint_rep%cross_overlap)
         call sparse_embed_create(temp_mat_3)
         call sparse_embed_transpose(temp_mat_3, joint_rep%cross_overlap)
         call sparse_embed_create(temp_mat_1, joint_rep%inv_overlap, temp_mat_3)
         call sparse_embed_product(temp_mat_1, joint_rep%inv_overlap, temp_mat_3)
         call sparse_embed_destroy(temp_mat_3)
         call sparse_embed_create(temp_mat_3, temp_mat_1, denskern%m(is,1))
         call sparse_embed_product(temp_mat_3, temp_mat_1, denskern%m(is,1))
         call sparse_embed_destroy(temp_mat_1)
         call sparse_embed_create(temp_mat_1, temp_mat_3, FO_overlap)
         call sparse_embed_product(temp_mat_1, temp_mat_3, FO_overlap)
         call sparse_embed_destroy(temp_mat_3)
         call sparse_embed_create(temp_mat_3, temp_mat_1, denskern%m(is,1))
         call sparse_embed_product(temp_mat_3, temp_mat_1, denskern%m(is,1))
         call sparse_embed_destroy(temp_mat_1)

         call sparse_embed_axpy(temp_mat_2, temp_mat_3, -0.5_DP)
         call sparse_embed_copy(respkern_jv(is), temp_mat_2)

         call sparse_embed_destroy(temp_mat_2)
         call sparse_embed_destroy(temp_mat_3)
      end do
      !! end of restriction
      !======================================================

    end do  ! end loop for optimisation

    !=======================================================

    ! recreate FO_ham_jv, using the previously converged values
    call sparse_embed_scale(FO_ham_jv, 0.0_DP)

    if (pub_paw) then
       do is = 1,pub_num_spins
          call sparse_embed_scale(FO_dij(is), 0.0_DP)
       end do
    end if

    FO_lhxc_fine = 0.0_DP

    call calc_FO_ham(mdl, val_rep, val_basis, direction, weight, dij, &
         FO_dij, dens_fine, denskern, rho_ij_dproj, rho_ij, respkern_jv, &
         FO_lpseudo_fine, lhxc_fine, FO_lhxc_fine, fxc_buffer, FO_ham_jv, &
         joint_rep, joint_basis, proj_basis, nl_projectors,FO_tcore_density)

    ! add previously-calculated fixed value
    call sparse_embed_axpy(FO_ham_jv, FO_ham_jv_fixed, 1.0_DP)

    ! get FO_ham from FO_ham_jv
    call sparse_embed_scale(FO_ham, 0.0_DP)
    call lr_phonons_get_basis_block(FO_ham_jv%p, FO_ham%p, joint_basis(1), &
         val_basis(1), val_basis(1), val_basis(1), 2, 1)
    !=====================================================

    call sparse_embed_destroy(A_mat_CT)
    call sparse_embed_destroy(B_mat_CT)
    call sparse_embed_destroy(C_mat_CT)
    call sparse_embed_destroy(D_mat_CT)

    do is = 1, pub_num_spins
       call sparse_embed_destroy(auxrespkern_jv_old(is))
       call sparse_embed_destroy(respkern_jv_old(is))
    end do

    call sparse_embed_destroy(F_mat)
    call sparse_embed_destroy(F_mat_CT)
    call sparse_embed_destroy(F_mat_old)
    call sparse_embed_destroy(F_mat_temp)
    call sparse_embed_destroy(F_mat_temp_CT)
    call sparse_embed_destroy(E_mat_var)
    call sparse_embed_destroy(E_mat_fixed)
    call sparse_embed_destroy(Z_mat)
    call sparse_embed_destroy(Z_mat_CT)
    call sparse_embed_destroy(XMINZ_mat)
    call sparse_embed_destroy(G_mat)
    call sparse_embed_destroy(Y_mat)
    call sparse_embed_destroy(Y_mat_diff)
    call sparse_embed_destroy(K_mat(1))
    call sparse_embed_destroy(K_mat(2))
    call sparse_embed_destroy(K_mat(3))
    call sparse_embed_destroy(XplusK_1)
    call sparse_embed_destroy(XplusK_2)

    call sparse_embed_destroy(A_lin_op)
    call sparse_embed_destroy(A_lin_op_temp)
    call sparse_embed_destroy(Af_lin_op)
    call sparse_embed_destroy(Af_lin_op_CT)
    call sparse_embed_destroy(Af_lin_op_sq)
    call sparse_embed_destroy(Af_lin_op_cb)
    call sparse_embed_destroy(Af_lin_op_sq_CT)
    call sparse_embed_destroy(Af_lin_op_cb_CT)
    call sparse_embed_destroy(grad_F_mat)
    call sparse_embed_destroy(grad_F_mat_old)
    call sparse_embed_destroy(grad_F_mat_CT)


    ! Start timer
    call timer_clock('lr_phonons_solve_rkn',2)

    if (pub_debug_on_root) write(stdout, '(a)') &
         'DEBUG: Leaving lr_phonons_solve_rkn routine.'

    call comms_barrier

  end subroutine lr_phonons_solve_rkn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calculate_first_overlap(mdl, val_ngwfs_on_grid, val_sp_overlap, &
       ngwf_basis, joint_ngwfs_on_grid, joint_sp_overlap, joint_basis, &
       proj_basis, nl_projectors, first_overlap_jv, first_overlap_vv, &
       iGps_overlap, jsiGp_overlap, direction, weight, rho_ij_dproj, denskern)

    !=========================================================================!
    ! This subroutine calculates the first-order overlap operator in v-v and  !
    ! j-v representations. This quantity exists due to the movement of non-   !
    ! local projectors in the PAW and USPP formalisms                         !
    !-------------------------------------------------------------------------!
    ! It can also compute the part of the first-order projector density-kernel!
    ! that derives from the rigid movement of the projector centers           !
    !-------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in 2015                               !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,        !
    ! September 2018                                                          !
    !=========================================================================!

    use datatypes, only: FUNCTIONS
    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use paw, only: paw_projector_overlap
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_num_spins, pub_paw, pub_usp, pub_debug_on_root, &
         pub_any_nl_proj
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_axpy, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product, sparse_embed_transpose, &
         sparse_embed_transpose_structure, sparse_embed_scale
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    integer, intent(in) :: direction(mdl%nat)
    real(kind=DP), intent(in) :: weight(mdl%nat)
    type(FUNCTIONS), intent(in) :: val_ngwfs_on_grid(1)
    type(SPAM3_EMBED), intent(in) :: val_sp_overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    type(FUNCTIONS), intent(in) :: joint_ngwfs_on_grid(1)
    type(SPAM3_EMBED), intent(in) :: joint_sp_overlap
    type(FUNC_BASIS), intent(in) :: joint_basis(1)
    type(FUNC_BASIS), intent(in) :: proj_basis(1)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(1)
    type(SPAM3_EMBED), intent(inout) :: first_overlap_jv
    type(SPAM3_EMBED), intent(inout) :: first_overlap_vv
    type(SPAM3_EMBED), intent(inout) :: iGps_overlap, jsiGp_overlap

    type(SPAM3_EMBED), intent(inout), optional :: rho_ij_dproj(pub_num_spins)
    type(SPAM3_EMBED), intent(inout), optional :: denskern(pub_num_spins)

    ! Local Variables
    type(SPAM3_EMBED) :: siGp_overlap
    type(SPAM3_EMBED) :: ps_overlap
    type(SPAM3_EMBED) :: oij
    type(SPAM3_EMBED) :: rho_ij_dproj_temp(pub_num_spins)
    type(SPAM3_EMBED) :: Ksp, KsiGp
    type(SPAM3_EMBED) :: prod_1, prod_2, prod_3, prod_4

    integer :: is

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering calculate_first_overlap'

    ! Start timer
    call timer_clock('calculate_first_overlap',1)

    if (pub_paw) then
       ! Create matrix arrays and structures
       oij%structure = 'E'
       call sparse_embed_create(oij)
       ! Get projector overlap matrix
       call paw_projector_overlap(oij%p,mdl%regions(1)%paw_sp)
    end if

    ! Create matrices to hold <phi|iG*proj> overlap matrix and transpose
    call sparse_embed_transpose_structure(ps_overlap%structure,val_sp_overlap)
    call sparse_embed_create(ps_overlap)
    call sparse_embed_transpose(ps_overlap,val_sp_overlap)
    call sparse_embed_create(siGp_overlap,val_sp_overlap)

    call sparse_embed_scale(siGp_overlap, 0.0_DP)

    ! Calculate <phi|iG*proj> overlap matrix
    call grad_proj_ovlp(mdl, siGp_overlap%p, val_ngwfs_on_grid(1), ngwf_basis(1), &
         proj_basis(1), nl_projectors(1), direction, weight)

    call grad_proj_ovlp(mdl, jsiGp_overlap%p, joint_ngwfs_on_grid(1), joint_basis(1), &
         proj_basis(1), nl_projectors(1), direction, weight)

    ! Transpose it to get <iG*proj|phi> overlap matrix
    call sparse_embed_transpose(iGps_overlap, siGp_overlap)

    if (pub_paw) then
       call sparse_embed_create(prod_1, oij, iGps_overlap)
       call sparse_embed_create(prod_2, prod_1)
       call sparse_embed_create(prod_3, jsiGp_overlap, prod_2)
       call sparse_embed_create(prod_4, siGp_overlap, prod_2)

       ! Calculate <proj_i|phi_a> O_ij <phi_b|-iG.proj_j>
       call sparse_embed_product(prod_1, oij, iGps_overlap)
       call sparse_embed_product(first_overlap_jv, joint_sp_overlap, prod_1)
       call sparse_embed_product(first_overlap_vv, val_sp_overlap, prod_1)

       ! Calculate <-iGproj_i|phi_a> O_ij <phi_b|proj_j>
       call sparse_embed_product(prod_2, oij, ps_overlap)
       call sparse_embed_product(prod_3, jsiGp_overlap, prod_2)
       call sparse_embed_axpy(first_overlap_jv, prod_3, 1.0_DP)

       call sparse_embed_product(prod_4, siGp_overlap, prod_2)
       call sparse_embed_axpy(first_overlap_vv, prod_4, 1.0_DP)

       call sparse_embed_destroy(prod_1)
       call sparse_embed_destroy(prod_2)
       call sparse_embed_destroy(prod_3)
       call sparse_embed_destroy(prod_4)
       call sparse_embed_destroy(oij)
    end if

    if (pub_paw.and.present(rho_ij_dproj).and.present(denskern)) then

       ! gcc32: scale denskern to correct number of electrons
       if(pub_num_spins == 1) call sparse_embed_scale(denskern(1), 2.0_DP)

       ! now create rho_ij_dproj
       do is = 1,pub_num_spins
          call sparse_embed_create(rho_ij_dproj_temp(is), rho_ij_dproj(is))

          call sparse_embed_create(KsiGp, denskern(is), siGp_overlap)
          call sparse_embed_product(KsiGp, denskern(is), siGp_overlap)
          call sparse_embed_product(rho_ij_dproj(is), ps_overlap, KsiGp)

          call sparse_embed_create(Ksp, denskern(is), val_sp_overlap)
          call sparse_embed_product(Ksp, denskern(is), val_sp_overlap)
          call sparse_embed_product(rho_ij_dproj_temp(is), iGps_overlap, Ksp)

          call sparse_embed_axpy(rho_ij_dproj(is), rho_ij_dproj_temp(is), 1.0_DP)

          call sparse_embed_destroy(KsiGp)
          call sparse_embed_destroy(Ksp)

          call sparse_embed_destroy(rho_ij_dproj_temp(is))

       end do

       ! gcc32: scale denskern to half number of electrons
       if(pub_num_spins == 1) call sparse_embed_scale(denskern(1), 0.5_DP)

    end if

    ! Destroy temporary matrices and deallocate arrays
    call sparse_embed_destroy(ps_overlap)
    call sparse_embed_destroy(siGp_overlap)

    ! Stop timer
    call timer_clock('calculate_first_overlap',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving calculate_first_overlap'

  end subroutine calculate_first_overlap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calculate_second_overlap(mdl, val_rep, ngwf_basis, proj_basis, &
       nl_projectors, second_overlap, iG1ps_overlap, iG2ps_overlap, &
       iG1iG2ps_overlap, direction1, weight1, direction2, weight2)

    !=========================================================================!
    ! This subroutine calculates the second-order overlap operator in v-v     !
    ! representation. This quantity exists due to the movement of non-local   !
    ! projectors in the PAW and USPP formalisms                               !
    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in 2015                               !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,        !
    ! October 2018                                                            !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use paw, only: paw_projector_overlap
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_num_spins, pub_paw, pub_usp, pub_debug_on_root,&
         pub_any_nl_proj
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_axpy, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product, sparse_embed_transpose, &
         sparse_embed_transpose_structure, sparse_embed_scale
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    type(NGWF_REP), intent(in) :: val_rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    type(FUNC_BASIS), intent(in) :: proj_basis(1)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(1)
    type(SPAM3_EMBED), intent(inout) :: second_overlap
    type(SPAM3_EMBED), intent(in) :: iG1ps_overlap
    type(SPAM3_EMBED), intent(inout) :: iG2ps_overlap
    type(SPAM3_EMBED), intent(inout) :: iG1iG2ps_overlap
    integer, intent(in) :: direction1(mdl%nat)
    real(kind=DP), intent(in) :: weight1(mdl%nat)
    integer, intent(in) :: direction2(mdl%nat)
    real(kind=DP), intent(in) :: weight2(mdl%nat)

    ! Local Variables
    type(SPAM3_EMBED) :: siG1p_overlap, siG2p_overlap, siG1iG2p_overlap
    type(SPAM3_EMBED) :: ps_overlap, oij
    type(SPAM3_EMBED) :: prod_1, prod_2

    integer :: is

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering calculate_second_overlap'

    ! Start timer
    call timer_clock('calculate_second_overlap',1)

    if (pub_paw) then
       ! Create matrix arrays and structures
       oij%structure = 'E'
       call sparse_embed_create(oij)
       ! Get projector overlap matrix
       call paw_projector_overlap(oij%p,mdl%regions(1)%paw_sp)
    end if


    ! Create matrices to hold <phi|iG*proj> overlap matrix and transpose
    call sparse_embed_transpose_structure(ps_overlap%structure, val_rep%sp_overlap)
    call sparse_embed_create(ps_overlap)
    call sparse_embed_transpose(ps_overlap, val_rep%sp_overlap)

    call sparse_embed_create(siG1p_overlap, val_rep%sp_overlap)
    call sparse_embed_create(siG2p_overlap, val_rep%sp_overlap)
    call sparse_embed_create(siG1iG2p_overlap, val_rep%sp_overlap)


    call grad_proj_ovlp(mdl, siG2p_overlap%p, val_rep%ngwfs_on_grid(1), &
         ngwf_basis(1), proj_basis(1), nl_projectors(1), direction2, weight2, &
         direction1, weight1, siG1iG2p_overlap%p)

    ! Transpose it to get <iG1/2*proj|phi> overlap matrix
    call sparse_embed_transpose(siG1p_overlap, iG1ps_overlap)
    call sparse_embed_transpose(iG2ps_overlap, siG2p_overlap)
    call sparse_embed_transpose(iG1iG2ps_overlap, siG1iG2p_overlap)


    if (pub_paw) then
       ! multiply with O_ij
       call sparse_embed_create(prod_1, oij, ps_overlap)
       call sparse_embed_create(prod_2, val_rep%sp_overlap, prod_1)
       call sparse_embed_scale(prod_2,0.0_DP)
       call sparse_embed_scale(prod_1,0.0_DP)

       ! Calculate <iG1proj_j|phi_b> O_ij <phi/theta_a|iG2proj_i>
       call sparse_embed_product(prod_1, oij, iG1ps_overlap)
       call sparse_embed_product(second_overlap, siG2p_overlap, prod_1)

       ! Calculate <iG2proj_j|phi_b> O_ij <phi/theta_a|iG1proj_i>
       call sparse_embed_product(prod_1, oij, iG2ps_overlap)
       call sparse_embed_product(prod_2, siG1p_overlap, prod_1)

       call sparse_embed_axpy(second_overlap, prod_2, 1.0_DP)
       call sparse_embed_scale(prod_2,0.0_DP)
       call sparse_embed_scale(prod_1,0.0_DP)

       ! Calculate <proj_j|phi_b> O_ij <phi/theta_a|iG1iG2proj_i>
       call sparse_embed_product(prod_1, oij, iG1iG2ps_overlap)
       call sparse_embed_product(prod_2, val_rep%sp_overlap, prod_1)

       call sparse_embed_axpy(second_overlap, prod_2, 1.0_DP)
       call sparse_embed_scale(prod_2,0.0_DP)
       call sparse_embed_scale(prod_1,0.0_DP)

       ! Calculate <proj_j|phi_b> O_ij <phi/theta_a|iG1iG2proj_i>
       call sparse_embed_product(prod_1, oij, ps_overlap)
       call sparse_embed_product(prod_2, siG1iG2p_overlap, prod_1)

       call sparse_embed_axpy(second_overlap, prod_2, 1.0_DP)

       call sparse_embed_destroy(prod_1)
       call sparse_embed_destroy(prod_2)
       call sparse_embed_destroy(oij)
    end if


    ! Destroy temporary matrices and deallocate arrays
    call sparse_embed_destroy(ps_overlap)
    call sparse_embed_destroy(siG1p_overlap)
    call sparse_embed_destroy(siG2p_overlap)
    call sparse_embed_destroy(siG1iG2p_overlap)

    ! Stop timer
    call timer_clock('calculate_second_overlap',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving calculate_second_overlap'

  end subroutine calculate_second_overlap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine grad_proj_ovlp(mdl, siG1p_overlap, ngwfs_on_grid, ngwf_basis, &
       proj_basis, proj_set, direction1, weight1, direction2, weight2, &
       siG1iG2p_overlap, factor_type)

    !========================================================================!
    ! This subroutine calculates the overlap matrix between the NGWFs and    !
    ! the first-order projectors (perturbed in a gived direction), and stores!
    ! it in a SPAM3 matrix.                                                  !
    !------------------------------------------------------------------------!
    ! Writted by Gabriel Constantinescu in 2015                              !
    ! Modified to remove pub_par by Joseph Prentice, October 2018            !
    !========================================================================!

    use basis, only: basis_phase_on_fftbox_recip, basis_point_wrt_box, &
         basis_start_of_box_wrt_cell
    use comms, only: pub_my_proc_id
    use constants, only: cmplx_i,DP,stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc, data_set_to_zero
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_brappd_ketfftbox
    use geometry, only: POINT, operator(*), operator(+)
    use model_type, only: MODEL
    use parallel_strategy, only: PARAL_INFO
    use projectors, only: PROJECTOR_SET, projectors_init_fftbox_recip, &
         projector_in_box_recip,projectors_exit_fftbox_recip
    use rundat, only: pub_fftbox_batch_size, pub_projectors_precalculate, &
         pub_threads_num_fftboxes, pub_threads_fftbox, pub_debug_on_root
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index, &
         sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_use_var, &
         utils_abort

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    type(SPAM3),intent(inout) :: siG1p_overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(PROJECTOR_SET), intent(inout) :: proj_set
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
    integer, intent(in) :: direction1(mdl%nat)
    real(kind=DP), intent(in) :: weight1(mdl%nat)
    integer, intent(in), optional :: direction2(mdl%nat)
    real(kind=DP), intent(in), optional :: weight2(mdl%nat)
    type(SPAM3), intent(inout), optional :: siG1iG2p_overlap
    integer, intent(in), optional :: factor_type

    ! Local Variables
    integer :: n_batches
    integer :: batch_size
    integer :: batch_count
    integer :: local_start
    integer :: local_end
    integer :: local_proj
    integer :: global_proj
    integer :: batch_proj
    integer :: atom_of_proj
    integer :: loc_atom_of_proj
    integer :: isp, atom
    integer :: iproj, shell, am, azim
    integer :: proj_count, proj_count_atom
    integer :: ierr             ! Error flag
    integer :: idx_len          ! Length of sparse index
    integer :: n1,n2,n3,ld1,ld2
    logical :: found
    ! {P}rojector {C}entre wrt fft {B}ox in terms of number of {G}ridpoints
    ! in each lattice vector direction. Not necessarily integer values.
    real(kind=DP) :: pcbg1,pcbg2,pcbg3
    type(FFTBOX_DATA), allocatable :: projector_fftbox(:)
    type(FFTBOX_DATA), allocatable :: projector_fftbox_so(:)
    integer, allocatable :: fftbox_start(:,:)              ! Ket start positions
    integer, allocatable :: sp_overlap_idx(:)              ! Sparse index
    complex(kind=DP), allocatable :: proj_complex(:,:,:)   ! Workspace
    complex(kind=DP), allocatable :: proj_complex_so(:,:,:)   ! Workspace

    ! jcap: parallel strategy for projector region
    type(PARAL_INFO), pointer :: proj_par

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &grad_proj_ovlp'

    ! Start timer
    call timer_clock('grad_proj_ovlp',1)

    ! ATTENTION: since they are projectors, it is assumed that they are real

    ! Initialise shorthand variables
    ld1 = mdl%fftbox%total_ld1
    ld2 = mdl%fftbox%total_ld2
    n1 = mdl%fftbox%total_pt1
    n2 = mdl%fftbox%total_pt2
    n3 = mdl%fftbox%total_pt3

    ! jcap: get correct par
    call sparse_get_par(proj_par, siG1p_overlap, 'C')

    ! Obtain index of siG1p_overlap
    idx_len = sparse_index_length(siG1p_overlap)
    allocate(sp_overlap_idx(idx_len), stat=ierr)
    call utils_alloc_check('grad_proj_ovlp','sp_overlap_idx',ierr)
    call sparse_generate_index(sp_overlap_idx, siG1p_overlap)

    ! Set batch size
    batch_size = pub_fftbox_batch_size

    ! Allocate workspace
    allocate(projector_fftbox(batch_size), stat=ierr)
    call utils_alloc_check('grad_proj_ovlp', 'projector_fftbox', ierr)
    do batch_count = 1,batch_size
       call data_fftbox_alloc(projector_fftbox(batch_count), &
            mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, mdl%fftbox%total_pt3, &
            iscmplx=.false.)
    end do

    if ( present(direction2) .and. (.not.present(siG1iG2p_overlap)) ) then
       call utils_abort('Must provide BOTH direction2 and siG1iG2p_overlap')
    end if

    if ( present(direction2) ) then
       allocate(projector_fftbox_so(batch_size), stat=ierr)
       call utils_alloc_check('grad_proj_ovlp','projector_fftbox_so',ierr)
       do batch_count = 1,batch_size
          call data_fftbox_alloc(projector_fftbox_so(batch_count), &
               mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, &
               mdl%fftbox%total_pt3, iscmplx=.false.)
       end do
    end if

    allocate(fftbox_start(3,batch_size), stat=ierr)
    call utils_alloc_check('grad_proj_ovlp', 'fftbox_start', ierr)

    if (pub_projectors_precalculate) then
       call projectors_init_fftbox_recip(proj_set, mdl%cell, mdl%fftbox)
    end if

    ! Number of row-steps per row-block
    n_batches = proj_basis%max_on_proc / batch_size
    if (mod(proj_basis%max_on_proc,batch_size) .gt. 0) n_batches = n_batches + 1

    local_start = 1
    do batch_count = 1,n_batches
       ! initialise:
       call data_set_to_zero(projector_fftbox(1:batch_size))
       if (present(direction2)) then
          call data_set_to_zero(projector_fftbox_so(1:batch_size))
       end if

       local_end = min(local_start + batch_size - 1, proj_basis%proc_num)

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(local_proj,batch_proj,global_proj,atom_of_proj, &
!$OMP      loc_atom_of_proj,isp,proj_count_atom,proj_count,pcbg1,pcbg2,pcbg3, &
!$OMP      proj_complex,proj_complex_so,ierr,iproj,found,am,azim,shell,atom) &
!$OMP SHARED(local_start,local_end,proj_set,ld1,ld2,n1,n2,n3,batch_size, &
!$OMP      proj_basis,pub_my_proc_id,proj_par, &
!$OMP      mdl,pub_projectors_precalculate,direction1,weight1,direction2, &
!$OMP      weight2, pub_threads_num_fftboxes, fftbox_start,projector_fftbox, &
!$OMP      projector_fftbox_so,pub_threads_fftbox)

       allocate(proj_complex(ld1, ld2, n3), stat=ierr)
       call utils_alloc_check('grad_proj_ovlp', 'proj_complex', ierr)
       proj_complex = cmplx(0.0_DP,0.0_DP,kind=DP)

       if (present(direction2)) then
          allocate(proj_complex_so(ld1, ld2, n3), stat=ierr)
          call utils_alloc_check('grad_proj_ovlp', 'proj_complex_so', ierr)
          proj_complex_so = cmplx(0.0_DP,0.0_DP,kind=DP)
       end if


!$OMP DO
       do local_proj = local_start,local_end

          batch_proj = local_proj - local_start + 1


          proj_complex = cmplx(0.0_DP,0.0_DP,kind=DP)
          if (present(direction2)) then
             proj_complex_so = cmplx(0.0_DP,0.0_DP,kind=DP)
          end if

          ! Find information about this projector
          global_proj = local_proj + proj_basis%first_on_proc(pub_my_proc_id)- 1
          atom_of_proj = proj_basis%atom_of_func(global_proj)
          loc_atom_of_proj = atom_of_proj - &
               proj_par%first_atom_on_proc(pub_my_proc_id) + 1
          isp = proj_set%proj_species(atom_of_proj)
          proj_count_atom = global_proj - &
               proj_basis%first_on_atom(atom_of_proj) + 1
          proj_count = proj_set%species_first_proj(isp) + proj_count_atom - 1

          ! Centre of projector wrt fftbox in terms of grid spacings
          call basis_point_wrt_box(pcbg1, pcbg2, pcbg3, &  ! output
               proj_set%proj_centre(atom_of_proj), n1, n2, n3, mdl%cell, &
               mdl%fftbox)

          ! Start of fftbox wrt cell in terms of grid-point number
          call basis_start_of_box_wrt_cell(fftbox_start(1,batch_proj), &
               fftbox_start(2,batch_proj), fftbox_start(3,batch_proj), &
               proj_set%proj_centre(atom_of_proj), pcbg1, pcbg2, &
               pcbg3, mdl%cell)

          if (pub_projectors_precalculate) then

             ! Apply phase factor e^-(ik.proj_centre_wrt_box) to projector
             ! stored in pre-calculated FFTbox
             call basis_phase_on_fftbox_recip(proj_complex, n1, n2, n3, ld1, &
                  ld2, -pcbg1, -pcbg2, -pcbg3, &
                  proj_set%fftbox_proj_recip(:,:,:,proj_count))
          else
             ! Find angular momentum of this projector
             iproj = 1
             found = .false.
             do shell = 1,proj_set%num_shells(isp)
                am = proj_set%ang_mom(shell,isp)
                do azim=-am,am
                   if (iproj==proj_count_atom) found=.true.
                   if (found) exit
                   iproj = iproj + 1
                end do
                if (found) exit
             end do

             ! Generate projector in FFTbox
             call projector_in_box_recip(proj_complex, &
                  proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                  shell,isp), am, azim, mdl%fftbox%total_pt1, &
                  mdl%fftbox%total_pt2, mdl%fftbox%total_pt3, &
                  proj_set%n_rad_pts(isp), proj_set%gmax(isp), mdl%cell)

             ! Apply phase factor e^-(ik.proj_centre_wrt_box) to projector
             call basis_phase_on_fftbox_recip(proj_complex, n1, n2, n3, ld1, &
                  ld2, -pcbg1, -pcbg2, -pcbg3, proj_complex)

          end if

          if ( present(direction2) ) then
             ! (optional) second order perturbation
             ! Multiply projector workspace by (-iG_dir1)*(-iG_dir2)
             proj_complex_so(:,:,:) = cmplx(-1.0_DP * mdl%fftbox%recip_grid( &
                  direction1(proj_par%orig_atom(atom_of_proj)),:,:,:) * &
                  mdl%fftbox%recip_grid( &
                  direction2(proj_par%orig_atom(atom_of_proj)),:,:,:) * &
                  weight1(proj_par%orig_atom(atom_of_proj)) * &
                  weight2(proj_par%orig_atom(atom_of_proj)), &
                  0.0_DP, kind=DP) * proj_complex(:,:,:)
          end if

          ! mandatory first order perturbation:
          ! Multiply projector workspace by -iG_dir1 and store as proj_complex
          ! scale by -1.0_DP since d/dR_I = -d/dr
          proj_complex(:,:,:) = cmplx(0.0_DP, -1.0_DP * &
               mdl%fftbox%recip_grid( &
               direction1(proj_par%orig_atom(atom_of_proj)),:,:,:) * &
               weight1(proj_par%orig_atom(atom_of_proj)), kind=DP)*&
               proj_complex(:,:,:)


          ! g=0 element must be real
          proj_complex(1,1,1) = cmplx( &
               real(proj_complex(1,1,1),kind=DP), 0.0_DP, kind=DP)

          if (present(direction2)) then
             ! g=0 element must be real
             proj_complex_so(1,1,1) = cmplx( &
                  real(proj_complex_so(1,1,1),kind=DP), 0.0_DP, kind=DP)
          end if

          ! Fourier transform to real space:
          ! ndmh: explicitly in-place transform
          call fourier_apply_box('Coarse','Backward',proj_complex, &
               omp=pub_threads_fftbox)

          if(present(direction2)) then
             call fourier_apply_box('Coarse','Backward',proj_complex_so, &
                  omp=pub_threads_fftbox)
          end if

          ! Put real part into projector_box
          ! assume from the start that they are real: use only d component
          projector_fftbox(batch_proj)%d(1:n1,1:n2,1:n3) &
               = real(proj_complex(1:n1,1:n2,1:n3),kind=DP)/mdl%cell%weight

          if (present(direction2)) then
             projector_fftbox_so(batch_proj)%d(1:n1,1:n2,1:n3) = &
                  real(proj_complex_so(1:n1,1:n2,1:n3),kind=DP)/mdl%cell%weight
          end if


       end do
!$OMP END DO

       deallocate(proj_complex,stat=ierr)
       call utils_dealloc_check('grad_proj_ovlp','proj_complex',ierr)

       if (present(direction2)) then
          deallocate(proj_complex_so,stat=ierr)
          call utils_dealloc_check('grad_proj_ovlp','proj_complex_so',ierr)
       end if

!$OMP END PARALLEL

       ! Calculate overlap integrals
       call function_ops_brappd_ketfftbox(siG1p_overlap, &           ! inout
            ngwfs_on_grid, ngwf_basis, mdl%cell, mdl%fftbox, &      ! input
            projector_fftbox(:), fftbox_start, &              ! input
            batch_size, local_start, local_end, idx_len, &          ! input
            sp_overlap_idx, 'FULL',.false., factor_type)

       if (present(direction2).and.present(siG1iG2p_overlap)) then
          ! Calculate overlap integrals
          call function_ops_brappd_ketfftbox(siG1iG2p_overlap, &       ! inout
               ngwfs_on_grid, ngwf_basis, mdl%cell, mdl%fftbox, &      ! input
               projector_fftbox_so(:), fftbox_start, &                 ! input
               batch_size, local_start, local_end, idx_len, &          ! input
               sp_overlap_idx, 'FULL', .false., factor_type)
       end if

       local_start = local_start + batch_size

    end do

    if (pub_projectors_precalculate) then
       call projectors_exit_fftbox_recip(proj_set)
    end if

    ! Deallocate workspace
    deallocate(fftbox_start,stat=ierr)
    call utils_dealloc_check('grad_proj_ovlp', &
         'fftbox_start',ierr)

    do batch_count = 1,batch_size
        call data_fftbox_dealloc(projector_fftbox(batch_count))
    end do
    deallocate(projector_fftbox,stat=ierr)
    call utils_dealloc_check('grad_proj_ovlp','projector_fftbox',ierr)

    if (present(direction2).and.present(siG1iG2p_overlap)) then
       do batch_count = 1,batch_size
           call data_fftbox_dealloc(projector_fftbox_so(batch_count))
       end do
       deallocate(projector_fftbox_so,stat=ierr)
       call utils_dealloc_check('grad_proj_ovlp','projector_fftbox_so',ierr)
    end if

    deallocate(sp_overlap_idx,stat=ierr)
    call utils_dealloc_check('grad_proj_ovlp','sp_overlap_idx',ierr)

    ! Stop timer
    call timer_clock('grad_proj_ovlp',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &grad_proj_ovlp'

  end subroutine grad_proj_ovlp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_FO_NGWFs(mdl, hfxstate, val_rep, ngwf_basis, directions, &
       weights, perturbation, proj_basis, proj_set, val_ham, dij, rho_ij, &
       dens_fine, denskern, hub_proj_basis, hub, lhxc_fine, FOngwf_basis, &
       FOngwf_rep, FO_lhxc_fine, FO_dij, viGp_overlap, &
       respkern_jv, auxrespkern_jv, vxc_buffer, fxc_buffer, FO_ham, &
       FO_ham_jv, iter_ngwfs, gradient_matrices, gradient_matrices_cov, &
       converged_cg, previous_g_dot_g, rms_gradient, line_search_success, &
       trial_length, prev_direction_on_grid, line_search_coeff, cg_count, &
       joint_basis, joint_rep, joint_ham, condkern_new, FO_locpot, &
       FO_tcore_density, SO_locpot_mat, dipole_mat_jv, at_gamma)


    ! ========================================================================!
    ! This subroutine performes the conjugate gradient (CG) optimisation of   !
    ! the first-order NGWFs. The algorithm is preconditioned (including the   !
    ! necessary corrections for the PAW projectors)                           !
    ! ------------------------------------------------------------------------!
    ! In the first CG iteration, the first-order NGWFs are simply intialised  !
    ! to the negative of their gradient (- d phi(r) / dr), resembling a Pulay !
    ! term.                                                                   !
    ! ------------------------------------------------------------------------!
    ! WARNING: while this initialisation is perfect for atomic displacements, !
    ! it seems to be less than ideal for E-field perturbations. More work is  !
    ! needed in that direction                                                !
    ! ------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in 2016-2017                          !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,        !
    ! October 2018                                                            !
    ! ========================================================================!

    use augmentation, only: augmentation_overlap
    use basis, only: basis_extract_function_from_box, basis_clean_function
    use comms, only: comms_barrier,pub_my_proc_id, comms_reduce, pub_on_root
    use constants, only: DP, LONG, stdout, paw_en_size, cmplx_0
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc, data_set_to_zero, data_functions_alloc, &
         data_functions_dealloc, data_functions_dot, data_coef_abs, COEF, &
         data_functions_copy, data_functions_axpy, data_functions_scale
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS, function_basis_ppds_to_sph_waves, &
         function_basis_est_num_psincs
    use function_ops, only: function_ops_batch_col_start, &
         function_ops_brappd_ketppd
    use hf_exchange, only: HFX_STATE
    use kinetic, only: kinetic_grad_to_func_batch
    use hamiltonian, only: hamiltonian_dens_indep_matrices, &
        hamiltonian_dens_dep_nonsc, hamiltonian_build_matrix, &
        hamiltonian_dens_dep_matrices
    use hubbard_build, only: HUBBARD_MODEL
    use model_type, only: MODEL
    use ngwfs, only: ngwfs_merge_sets
    use ngwf_representation, only: NGWF_REP, NGWF_HAM, ngwf_rep_create, &
        ngwf_rep_destroy, ngwf_ham_create
    use projectors, only: PROJECTOR_SET
    use restart, only: restart_ngwfs_tightbox_input
    use rundat, only: pub_num_spins,pub_debug_on_root,pub_fftbox_batch_size, &
         pub_paw, pub_aug, pub_elec_cg_max, pub_output_detail, &
         pub_ngwf_cg_max_step, pub_precond_scheme, pub_smooth_scheme, pub_k_zero, &
         pub_any_nl_proj
    use services, only: services_line_search_parabola, services_flush, &
         services_cubic_fit_minimum
    use sparse, only: sparse_generate_index, sparse_index_length
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_axpy, sparse_embed_copy, &
         sparse_embed_transpose, sparse_embed_transpose_structure, &
         sparse_embed_product, sparse_embed_scale, sparse_embed_array_scale
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit
    use visual, only: visual_ngwfs

    implicit none

    ! Arguments
    type(MODEL), intent(inout) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    type(NGWF_REP), intent(in) :: val_rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    integer, intent(in) :: directions(3*mdl%nat+3, mdl%nat)
    real(kind=DP), intent(in) :: weights(3*mdl%nat+3, mdl%nat)
    integer, intent(in) :: perturbation
    type(FUNC_BASIS), intent(in) :: proj_basis(1)
    type(PROJECTOR_SET), intent(inout) :: proj_set(1)
    type(NGWF_HAM), intent(inout) :: val_ham
    type(SPAM3_EMBED), intent(inout) :: dij(pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: rho_ij(pub_num_spins)
    real(kind=DP), intent(inout) :: dens_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins)
    type(SPAM3_EMBED_ARRAY), intent(inout) :: denskern
    type(FUNC_BASIS), intent(inout) :: hub_proj_basis(1)
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    type(NGWF_REP), intent(inout) :: FOngwf_rep
    type(FUNC_BASIS), intent(inout) :: FOngwf_basis(1)
    real(kind=DP), intent(in) :: FO_lhxc_fine(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: FO_dij(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: viGp_overlap
    type(SPAM3_EMBED), intent(inout) :: respkern_jv(pub_num_spins), &
         auxrespkern_jv(pub_num_spins)
    real(kind=DP), intent(in) :: vxc_buffer(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in) :: fxc_buffer(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: FO_ham, FO_ham_jv
    integer, intent(in) :: iter_ngwfs
    type(SPAM3_EMBED), intent(in) :: gradient_matrices(7), gradient_matrices_cov(7)
    integer, intent(inout) :: cg_count
    real(kind=DP), intent(inout) :: line_search_coeff
    logical, intent(inout) :: converged_cg, line_search_success
    real(kind=DP), intent(inout) :: rms_gradient, trial_length
    type(COEF), intent(inout) :: previous_g_dot_g
    type(FUNCTIONS), intent(inout) :: prev_direction_on_grid
    type(FUNC_BASIS), intent(inout) :: joint_basis(1)
    type(NGWF_REP), intent(inout) :: joint_rep
    type(NGWF_HAM), intent(inout) :: joint_ham
    type(SPAM3_EMBED_ARRAY), intent(inout) :: condkern_new
    real(kind=DP), intent(in) :: FO_locpot(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12)
    real(kind=DP), intent(in) :: FO_tcore_density(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12)
    type(SPAM3_EMBED), intent(in) :: SO_locpot_mat, dipole_mat_jv(3)
    logical, intent(in) :: at_gamma

    ! Local Variables
    type(FFTBOX_DATA), allocatable, dimension(:,:) :: grad_fftbox_batch
    integer :: batch_size
    integer :: batch_count, local_fa, batch_fa
    integer :: n_batches, atom_idx, global_idx, list_idx
    type(SPAM3_EMBED) :: pr_overlap, pv_overlap
    integer :: atom_of_func, global_func
    integer :: local_start, local_end
    integer :: ierr         ! pdh: error flag
    integer :: idx_len      ! pdh: length of sparse index
    integer, allocatable :: ngwf_box_start(:,:)
    integer, allocatable :: ngwf_start_in_box(:,:)
    character(20) :: pattern
    integer :: is,i1,i2,i3
    logical :: loc_cmplx
    type(SPAM3_EMBED) :: FO_dijps
    integer :: iter, iter_1, iter_2, output_unit_test, io_status
    type(SPAM3_EMBED) :: temp_mat_1, temp_mat_2, temp_mat_3, iGpv_overlap

    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    real(kind=DP) :: total_energy, lhxc_energy, hubbard_energy

    type(FUNCTIONS) :: contra_grad, cov_grad, start_ngwfs_on_grid
    type(COEF) :: current_g_dot_g
    integer(kind=LONG) :: est_num_psincs
    type(COEF) :: G_init, cg_coeff
    logical :: trial2, retrial1, reversing
    real(kind=DP) :: quadratic_coeff, predicted_functional, cubic_coeff, &
         rejected_quadratic_coeff, F0, F1, F2

    type(FUNCTIONS) :: direction_on_grid

    ! FD variables, for dev use only
    integer, parameter :: nfd=1
    integer :: ifd
    real(kind=DP) :: fd_trial_step(25)=(/0.000001_DP,0.00001_DP,0.0001_DP, &
         0.001_DP,0.01_DP,0.02_DP,0.04_DP,0.06_DP,0.08_DP,0.1_DP, &
         0.15_DP,0.2_DP,0.3_DP,0.4_DP,0.5_DP,0.6_DP,0.7_DP,0.8_DP,0.9_DP, &
         1.0_DP,1.2_DP,1.4_DP,1.6_DP,1.8_DP,2.0_DP/)
    real(kind=DP) :: FFD(nfd)

    est_num_psincs = 0_LONG
    est_num_psincs = function_basis_est_num_psincs(FOngwf_basis(1),mdl%cell)
    call data_functions_alloc(direction_on_grid, FOngwf_basis(1)%size_on_grid, &
         iscmplx=FOngwf_rep%ngwfs_on_grid(1)%iscmplx)


    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering calc_FO_NGWFs'

    ! Start timer
    call timer_clock('calc_FO_NGWFs',1)

    loc_cmplx = FOngwf_rep%ngwfs_on_grid(1)%iscmplx

    ! Define shorthand variables
    batch_size = pub_fftbox_batch_size

    ! ddor: FULL integral is calculated for non-square matrices
    ! lr408: changed to different basis names
    pattern = 'FULL'

    call sparse_embed_transpose_structure(iGpv_overlap%structure, viGp_overlap)
    call sparse_embed_create(iGpv_overlap)
    call sparse_embed_transpose(iGpv_overlap, viGp_overlap)

    if (iter_ngwfs == 1) then

       if (loc_cmplx) then
          FOngwf_rep%ngwfs_on_grid(1)%z(:) = cmplx_0
       else
          FOngwf_rep%ngwfs_on_grid(1)%d(:) = 0.0_DP
       end if

       ! Allocate workspace
       allocate(grad_fftbox_batch(batch_size,3), stat=ierr)
       call utils_alloc_check('calc_FO_NGWFs', 'grad_fftbox_batch', ierr)
       do batch_count = 1,batch_size
          call data_fftbox_alloc(grad_fftbox_batch(batch_count,1), &
               mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, &
               mdl%fftbox%total_pt3, iscmplx=loc_cmplx)
          call data_fftbox_alloc(grad_fftbox_batch(batch_count,2), &
               mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, &
               mdl%fftbox%total_pt3, iscmplx=loc_cmplx)
          call data_fftbox_alloc(grad_fftbox_batch(batch_count,3), &
               mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, &
               mdl%fftbox%total_pt3, iscmplx=loc_cmplx)
       end do

       allocate(ngwf_box_start(3,batch_size), stat=ierr)
       call utils_alloc_check('calc_FO_NGWFs', 'ngwf_box_start', ierr)
       allocate(ngwf_start_in_box(3,batch_size), stat=ierr)
       call utils_alloc_check('calc_FO_NGWFs', 'ngwf_start_in_box', ierr)

       ! cks: number of row-steps per row-block
       n_batches = FOngwf_basis(1)%proc_num / batch_size
       if (mod(FOngwf_basis(1)%max_on_proc,batch_size) .gt. 0) &
            n_batches = n_batches + 1

       call comms_barrier

       local_start = 1
       do batch_count = 1,n_batches ! loop over batches

          local_end = min(local_start + batch_size - 1, FOngwf_basis(1)%proc_num)

          ! gcc32: initialise:
          call data_set_to_zero(grad_fftbox_batch(1:batch_size,1))
          call data_set_to_zero(grad_fftbox_batch(1:batch_size,2))
          call data_set_to_zero(grad_fftbox_batch(1:batch_size,3))


          call function_ops_batch_col_start(ngwf_box_start, ngwf_start_in_box, &
               batch_size, local_start, local_end,mdl%fftbox, mdl%cell, &
               FOngwf_basis(1))

          ! cks: apply grad operator on to kets
          call kinetic_grad_to_func_batch(grad_fftbox_batch, &
               val_rep%ngwfs_on_grid(1), FOngwf_basis(1), mdl%cell,mdl%fftbox, &
               ngwf_start_in_box, batch_size, local_start, local_end)

          do local_fa = local_start,local_end

             batch_fa = local_fa - local_start + 1
             global_func = local_fa+FOngwf_basis(1)%first_on_proc(pub_my_proc_id)-1
             atom_of_func = FOngwf_basis(1)%atom_of_func(global_func)

             if (.not.loc_cmplx) then
                grad_fftbox_batch(batch_fa, &
                     directions(perturbation,mdl%par%orig_atom(atom_of_func)))%d(:,:,:) = &
                     grad_fftbox_batch(batch_fa, &
                     directions(perturbation,mdl%par%orig_atom(atom_of_func)))%d * &
                     (-weights(perturbation,mdl%par%orig_atom(atom_of_func)))
             else
                grad_fftbox_batch(batch_fa, &
                     directions(perturbation,mdl%par%orig_atom(atom_of_func)))%z(:,:,:) = &
                     grad_fftbox_batch(batch_fa, &
                     directions(perturbation,mdl%par%orig_atom(atom_of_func)))%z * &
                     cmplx(-weights(perturbation, mdl%par%orig_atom(atom_of_func)),&
                     0.0_DP,kind=DP)
             end if

             ! the '-' is due to the fact that we calculated the GRAD
             ! of the ngwfs, but we need d/dR = -GRAD

             ! extract the \nabla NGWF_dir from FFTboxes into PPDs
             call basis_extract_function_from_box(FOngwf_rep%ngwfs_on_grid(1), &
                  grad_fftbox_batch(batch_fa, &
                  directions(perturbation, mdl%par%orig_atom(atom_of_func))), &
                  FOngwf_basis(1)%spheres(local_fa), &
                  FOngwf_basis(1)%tight_boxes(local_fa), &
                  ngwf_start_in_box(1,batch_fa),ngwf_start_in_box(2,batch_fa), &
                  ngwf_start_in_box(3,batch_fa), &
                  FOngwf_basis(1)%spheres(local_fa)%offset, mdl%cell, mdl%fftbox)

              ! shave
              call basis_clean_function(FOngwf_rep%ngwfs_on_grid(1), &
                   FOngwf_basis(1)%spheres(local_fa), mdl%cell,mdl%fftbox)

          end do

          local_start = local_start + batch_size

       end do ! end loop over batches

       ! pdh: deallocate workspace
       deallocate(ngwf_start_in_box,stat=ierr)
       call utils_dealloc_check('calc_FO_NGWFs', 'ngwf_start_in_box', ierr)
       deallocate(ngwf_box_start,stat=ierr)
       call utils_dealloc_check('calc_FO_NGWFs', 'ngwf_box_start', ierr)

       do batch_count = 1,batch_size
          call data_fftbox_dealloc(grad_fftbox_batch(batch_count,1))
          call data_fftbox_dealloc(grad_fftbox_batch(batch_count,2))
          call data_fftbox_dealloc(grad_fftbox_batch(batch_count,3))
       end do
       deallocate(grad_fftbox_batch, stat=ierr)
       call utils_dealloc_check('calc_FO_NGWFs', 'grad_fftbox_batch', ierr)

       call comms_barrier

    else

       trial2 = .false.
       retrial1 = .false.
       reversing = .false.

       quadratic_coeff = 0.0_DP
       predicted_functional = 0.0_DP
       rejected_quadratic_coeff = 0.0_DP
       cubic_coeff = 0.0_DP

       cg_coeff%iscmplx = val_rep%ngwfs_on_grid(1)%iscmplx
       G_init%iscmplx = val_rep%ngwfs_on_grid(1)%iscmplx
       call data_set_to_zero(G_init)
       call data_set_to_zero(cg_coeff)

       call data_functions_alloc(contra_grad, FOngwf_basis(1)%size_on_grid, &
            loc_cmplx)
       call data_functions_alloc(cov_grad, FOngwf_basis(1)%size_on_grid, &
            loc_cmplx)

       call data_set_to_zero(contra_grad)
       call data_set_to_zero(cov_grad)

       call FO_ngwf_gradient_lnv(contra_grad, cov_grad, proj_basis, proj_set, &
            lhxc_fine, FO_lhxc_fine, dij, FO_dij, mdl, val_rep, ngwf_basis, &
            viGp_overlap, gradient_matrices, gradient_matrices_cov, &
            FOngwf_basis, FOngwf_rep, directions(perturbation,:), &
            weights(perturbation,:), perturbation, at_gamma)

       call comms_barrier


       current_g_dot_g%iscmplx = FOngwf_rep%ngwfs_on_grid(1)%iscmplx
       call data_set_to_zero(current_g_dot_g)

       call comms_barrier

       current_g_dot_g = data_functions_dot(contra_grad, cov_grad)

       if (current_g_dot_g%iscmplx) then
          call comms_reduce('SUM', current_g_dot_g%z)
       else
          call comms_reduce('SUM', current_g_dot_g%d)
       end if

       ! if(pub_on_root) write(stdout,*) 'VIEW current_g_dot_g: ', &
       !    current_g_dot_g%d

       rms_gradient = 0.0_DP
       rms_gradient = sqrt(data_coef_abs(current_g_dot_g) / &
            real(est_num_psincs, kind=DP))

       if (pub_on_root) write(stdout,*) 'RMS: ', rms_gradient

       if (rms_gradient < 0.000001_DP) converged_cg = .TRUE.

       call second_order_energy(mdl, val_rep, ngwf_basis, val_ham, &
            joint_basis, joint_rep, joint_ham, directions, weights, &
            perturbation, dij, FO_dij, FO_ham, FO_ham_jv, lhxc_fine, &
            FO_lhxc_fine, vxc_buffer, fxc_buffer, denskern%m(1,1), &
            respkern_jv(1), proj_basis, proj_set, iGpv_overlap, &
            dipole_mat_jv, at_gamma, only_SO_energy = F0, &
            locpot_SO_mat = SO_locpot_mat)

       call FO_find_direction

       if(pub_on_root) write(stdout,*) 'G_init: ', G_init%d
       ! if(pub_on_root) write(stdout,*) 'VIEW cg_coeff: ', cg_coeff%d

       ! cks: store direction if doing conjugate gradients
       if (pub_elec_cg_max > 0) then
          call data_functions_copy(prev_direction_on_grid, direction_on_grid)
       endif

       previous_g_dot_g = current_g_dot_g

       call data_functions_alloc(start_ngwfs_on_grid,FOngwf_basis(1)%size_on_grid,&
            iscmplx=FOngwf_rep%ngwfs_on_grid(1)%iscmplx)

       call data_functions_copy(start_ngwfs_on_grid, FOngwf_rep%ngwfs_on_grid(1))
!      if(pub_on_root) write(stdout,*) 'START OF FD TEST ================='
!      !!!!! =================== FD TEST ====================================
!      do ifd = 1,nfd
!
!         if (pub_on_root) write(stdout,*) 'FD STEP: ',  fd_trial_step(ifd)
!
!         FFD(ifd) = 0.0_DP
!
!         call data_functions_copy(FOngwf_rep%ngwfs_on_grid, &
!              start_ngwfs_on_grid)
!         call data_functions_axpy(FOngwf_rep%ngwfs_on_grid, &
!              direction_on_grid,fd_trial_step(ifd))
!
!         call comms_barrier
!
!         call ngwfs_merge_sets(joint_rep%ngwfs_on_grid, joint_basis, &
!              mdl%cell, FOngwf_rep%ngwfs_on_grid, FOngwf_basis, &
!              val_rep%ngwfs_on_grid, ngwf_basis)
!
!         call hamiltonian_dens_indep_matrices(FOngwf_rep, FOngwf_basis, &
!              proj_basis, proj_set, hub_proj_basis, hub, mdl)
!
!         call hamiltonian_dens_indep_matrices(joint_rep, joint_basis, &
!              proj_basis, proj_set, hub_proj_basis, hub, mdl, val_rep, &
!              ngwf_basis)
!         call hamiltonian_dens_dep_nonsc(joint_ham, joint_rep, &
!              joint_basis, mdl, lhxc_fine, hub, val_rep, val_ham, &
!              denskern%m(:,1))
!
!         call comms_barrier
!
!         call lr_phonons_regenerate_SO_energy(mdl, ngwf_basis, val_rep, &
!              val_ham, joint_basis, joint_rep, joint_ham,auxrespkern_jv,&
!              directions, weights, perturbation, proj_basis, proj_set, &
!              dij, denskern, dens_fine, lhxc_fine, vxc_buffer, &
!              FO_locpot, FO_tcore_density, SO_locpot_mat, at_gamma, FFD(ifd))
!      end do
!      !!!!! =================== END FD TEST ================================
!      if(pub_on_root) write(stdout,*) 'end OF FD TEST ================='

       call data_functions_copy(FOngwf_rep%ngwfs_on_grid(1), &
            start_ngwfs_on_grid)
       call data_functions_axpy(FOngwf_rep%ngwfs_on_grid(1), &
            direction_on_grid, trial_length)

       call comms_barrier

       call ngwfs_merge_sets(joint_rep%ngwfs_on_grid(1), joint_basis(1), &
            mdl%cell, FOngwf_rep%ngwfs_on_grid(1), FOngwf_basis(1), &
            val_rep%ngwfs_on_grid(1), ngwf_basis(1), par=mdl%par)

       call comms_barrier

       call hamiltonian_dens_indep_matrices(FOngwf_rep, FOngwf_basis, &
            proj_basis, proj_set, hub_proj_basis, hub, mdl)

       call hamiltonian_dens_indep_matrices(joint_rep, joint_basis, &
            proj_basis, proj_set, hub_proj_basis, hub, mdl, val_rep, ngwf_basis)
       call hamiltonian_dens_dep_nonsc(joint_ham, joint_rep, &
            joint_basis, mdl, hfxstate, lhxc_fine, hub, val_rep, val_ham, &
            denskern%m(:,1))

       call comms_barrier

       call lr_phonons_regenerate_SO_energy(mdl, ngwf_basis, val_rep, val_ham, &
            joint_basis, joint_rep, joint_ham, auxrespkern_jv, &
            directions, weights, perturbation, proj_basis, proj_set, dij, &
            denskern, dens_fine, lhxc_fine, vxc_buffer, FO_locpot, &
            FO_tcore_density, SO_locpot_mat, at_gamma, F1)

       call services_line_search_parabola(quadratic_coeff,  &
            predicted_functional,line_search_success, G_init%d, F0, F1, &
            trial_length, pub_ngwf_cg_max_step)

       line_search_coeff = quadratic_coeff

       if (pub_on_root) write(stdout,*) 'QUADRATIC LINE_SEARCH_COEFF: ', &
            line_search_coeff

       ! cks: CUBIC CUBIC CUBIC ----- SEARCH BY FITTING CUBIC ---- CUBIC CUBIC
       if (quadratic_coeff * G_init%d > 0.0_DP) then !jmecmplx

          trial2 = .true.

          ! ndmh_pointerfun
          call data_functions_copy(FOngwf_rep%ngwfs_on_grid(1), &
               start_ngwfs_on_grid)
          call data_functions_axpy(FOngwf_rep%ngwfs_on_grid(1), &
               direction_on_grid, 2.0_DP*trial_length)
       end if

       if ( ((line_search_coeff < 0.05_DP*trial_length) .or. &
            (line_search_coeff > 10.0_DP*trial_length) .or. &
            (.not. line_search_success) .or. &
            (reversing .and. (line_search_coeff > trial_length))) &
            .and. .not. trial2 ) then

          retrial1 = .true.
       end if

       if (retrial1) then
          rejected_quadratic_coeff = quadratic_coeff
          ! ndmh_pointerfun
          call data_functions_copy(FOngwf_rep%ngwfs_on_grid(1), &
               start_ngwfs_on_grid)
          call data_functions_axpy(FOngwf_rep%ngwfs_on_grid(1), &
               direction_on_grid, quadratic_coeff)
       end if
       ! ndmh: SECOND TRIAL STEP ----- SECOND TRIAL STEP ----- SECOND TRIAL STEP

       if (trial2 .or. retrial1) then

          call comms_barrier

          call hamiltonian_dens_indep_matrices(FOngwf_rep, FOngwf_basis, &
               proj_basis, proj_set, hub_proj_basis, hub, mdl)

          call ngwfs_merge_sets(joint_rep%ngwfs_on_grid(1), joint_basis(1), &
               mdl%cell, FOngwf_rep%ngwfs_on_grid(1), FOngwf_basis(1), &
               val_rep%ngwfs_on_grid(1), ngwf_basis(1), par=mdl%par)

          call comms_barrier

          call hamiltonian_dens_indep_matrices(joint_rep, joint_basis, &
               proj_basis, proj_set, hub_proj_basis, hub, mdl, val_rep, &
               ngwf_basis)
          call hamiltonian_dens_dep_nonsc(joint_ham, joint_rep, &
               joint_basis, mdl, hfxstate, lhxc_fine, hub, val_rep, val_ham, &
               denskern%m(:,1))

          call comms_barrier

          call lr_phonons_regenerate_SO_energy(mdl, ngwf_basis, val_rep, &
               val_ham, joint_basis, joint_rep, joint_ham,auxrespkern_jv,&
               directions, weights, perturbation, proj_basis, proj_set, dij, &
               denskern, dens_fine, lhxc_fine, vxc_buffer, FO_locpot, &
               FO_tcore_density, SO_locpot_mat, at_gamma, F2)

          if (trial2) then
             ! Cubic Fit
             call services_cubic_fit_minimum( &
                  cubic_coeff, predicted_functional, line_search_success, &
                  F0, F1, F2, G_init%d, trial_length, 2.0_DP*trial_length, &
                  pub_ngwf_cg_max_step)

             line_search_coeff = cubic_coeff
             if (pub_on_root) write(stdout,*) 'CUBIC LINE_SEARCH_COEFF: ', &
                line_search_coeff
          end if

          if (retrial1) then
             ! ndmh: quadratic fit at new trial length
             call services_line_search_parabola(quadratic_coeff, &
                  predicted_functional, line_search_success, G_init%d, F0, F2,&
                  line_search_coeff, pub_ngwf_cg_max_step)

             line_search_coeff = quadratic_coeff
          end if

       else
          ! ndmh: no second trial step
          F2 = 0.0_DP
       end if
       ! cks: CUBIC CUBIC --------- END SEARCH BY FITTING CUBIC -------- CUBIC

       call data_functions_copy(FOngwf_rep%ngwfs_on_grid(1), &
            start_ngwfs_on_grid)
       call data_functions_axpy(FOngwf_rep%ngwfs_on_grid(1), &
            direction_on_grid, line_search_coeff)

       ! ndmh: reset flags
       trial2 = .false.
       retrial1 = .false.
       reversing = .false.
       F2 = 0.0_DP

       call data_functions_dealloc(start_ngwfs_on_grid)
       call data_functions_dealloc(contra_grad)
       call data_functions_dealloc(cov_grad)

    end if

    call comms_barrier

    call ngwfs_merge_sets(joint_rep%ngwfs_on_grid(1), joint_basis(1), &
         mdl%cell, FOngwf_rep%ngwfs_on_grid(1), FOngwf_basis(1), &
         val_rep%ngwfs_on_grid(1), ngwf_basis(1), par=mdl%par)

    call hamiltonian_dens_indep_matrices(FOngwf_rep, FOngwf_basis, &
         proj_basis, proj_set, hub_proj_basis, hub, mdl)

    call comms_barrier

    call hamiltonian_dens_indep_matrices(joint_rep, joint_basis, &
         proj_basis, proj_set, hub_proj_basis, hub, mdl, val_rep, ngwf_basis)
    call hamiltonian_dens_dep_nonsc(joint_ham, joint_rep, joint_basis, &
         mdl, hfxstate, lhxc_fine, hub, val_rep, val_ham, denskern%m(:,1))

    ! pdh: sync procs
    call comms_barrier

    ! create condkern_new
    call sparse_embed_copy(condkern_new%m(1,1), joint_rep%inv_overlap)
    call sparse_embed_create(temp_mat_1, joint_rep%cross_overlap, &
         joint_rep%inv_overlap)
    call sparse_embed_product(temp_mat_1, joint_rep%cross_overlap, &
         joint_rep%inv_overlap)
    call sparse_embed_transpose_structure(temp_mat_3%structure, &
         joint_rep%cross_overlap)
    call sparse_embed_create(temp_mat_3)
    call sparse_embed_create(temp_mat_2, joint_rep%inv_overlap, temp_mat_3)
    call sparse_embed_destroy(temp_mat_3)
    call sparse_embed_transpose(temp_mat_2, temp_mat_1)
    call sparse_embed_create(temp_mat_3, temp_mat_2, denskern%m(1,1))
    call sparse_embed_product(temp_mat_3, temp_mat_2, denskern%m(1,1))
    call sparse_embed_destroy(temp_mat_2)
    call sparse_embed_create(temp_mat_2, temp_mat_3, temp_mat_1)
    call sparse_embed_product(temp_mat_2, temp_mat_3, temp_mat_1)
    call sparse_embed_axpy(condkern_new%m(1,1), temp_mat_2, -1.0_DP)
    call sparse_embed_destroy(temp_mat_1)
    call sparse_embed_destroy(temp_mat_2)
    call sparse_embed_destroy(temp_mat_3)

    call comms_barrier

    call sparse_embed_destroy(iGpv_overlap)

    ! Stop timer
    call timer_clock('calc_FO_NGWFs',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving calc_FO_NGWFs'


  contains

    subroutine FO_find_direction

      use datatypes, only: data_coef_abs, data_functions_axpy, &
           data_set_to_zero, data_functions_scale, data_coef_scale
      use constants, only: cmplx_0

      implicit none

       if (.not.line_search_success) then
          trial_length = trial_length * 0.5_DP
       else
          trial_length = max(sqrt(trial_length*line_search_coeff), &
               epsilon(1.0_DP))
       end if
      trial_length = max(trial_length,0.0001_DP)

      if (pub_on_root) then
         write(stdout,'(a)') '>>> Improving FO NGWFs using line search:'
      endif

      if(pub_on_root) write(stdout,*) 'TRIAL LENGTH: ', trial_length

      ! calculate CG coefficient
      if ( iter_ngwfs > 2 ) then
         if ((cg_count >= pub_elec_cg_max) .or. (.not.line_search_success)) then
            ! cks: reset CG after "cg_max" steps or after a fitting failure
            call data_set_to_zero(cg_coeff)
            cg_count = 0
            if (pub_on_root  .and. (.not.line_search_success) ) &
                 write(stdout,'(a)') 'Resetting FO_NGWF CG'
         else
            ! cks: original Fletcher-Reeves formula (cheaper in memory)
            if (data_coef_abs(previous_g_dot_g) > epsilon(1.0_DP)) then
               ! ndmh_pointerfun
               cg_coeff%d = current_g_dot_g%d / previous_g_dot_g%d
               ! cks: protection from crazy coefficients
               if (data_coef_abs(cg_coeff) > 10.0_DP ) then
                  if (pub_on_root) then
                     if (cg_coeff%iscmplx) then
                        write(stdout,'(a,f8.4,1x,f8.4,a)') 'WARNING: FO_NGWF&
                             & Fletcher-Reeves CG coeff too large (',  &
                             cg_coeff%z, ') - setting to zero'
                     else
                        write(stdout,'(a,f8.4,a)') 'WARNING: FO_NGWF&
                             & Fletcher-Reeves CG coeff too large (',  &
                             cg_coeff%d, ') - setting to zero'
                     end if
                  end if
                  call data_set_to_zero(cg_coeff)
                  cg_count = 0
               end if
            else
               if (pub_on_root)  write(stdout,*) ' previous_g_dot_g=', &
                    previous_g_dot_g%d, &
                    'CG coeffient set to zero in ngwf_cg_optimise'
               call data_set_to_zero(cg_coeff)
               cg_count = 0
            end if

            ! cks: re-initialise the periodic reset process if cg_coeff was zero
            ! cks: otherwise increase cg_count
            if ( (cg_coeff%iscmplx .and. (cg_coeff%z == cmplx_0)) .or. &
                 ((.not. cg_coeff%iscmplx) .and. (cg_coeff%d == 0.0_DP)) ) then
               cg_count = 0
            else
               cg_count = cg_count + 1
            end if

         end if
      else
         call data_set_to_zero(cg_coeff)
      endif

      ! Find search direction
      if (pub_elec_cg_max > 0) then
         ! ndmh_pointerfun
         call data_functions_copy(direction_on_grid, cov_grad)
         call data_functions_scale(direction_on_grid, -1.0_DP)
         call data_functions_axpy(direction_on_grid, prev_direction_on_grid, &
              cg_coeff)
         !direction_on_grid = cg_coeff * prev_direction_on_grid - cov_grad
      else
         ! cks:steepest descents
         call data_set_to_zero(direction_on_grid)
         call data_functions_axpy(direction_on_grid, cov_grad, -1.0_DP)
      endif

      ! Slope of energy in search direction
      G_init = data_functions_dot(contra_grad, direction_on_grid)

      ! cks: collect the work of each proc
      if (G_init%iscmplx) then
         call comms_reduce('SUM', G_init%z)
      else
         call comms_reduce('SUM', G_init%d)
      end if

      ! take action in case of positive slope along search direction
      if (G_init%d > 0.0_DP) then                  !jmecmplx

         if (pub_on_root) then
            write(stdout,'(a,e16.6)') &
                 'WARNING: slope along search direction is positive:', G_init%d
            write(stdout,'(a)') '         Resetting conjugate gradients!'
         end if
         call data_set_to_zero(direction_on_grid)
         call data_functions_axpy(direction_on_grid, cov_grad, -1.0_DP)

         cg_count = 0
         G_init = data_functions_dot(contra_grad, direction_on_grid)

         ! cks: collect the work of each proc
         if (G_init%iscmplx) then
            call comms_reduce('SUM', G_init%z)
         else
            call comms_reduce('SUM', G_init%d)
         end if

         if (G_init%d > 0.0_DP) then
            if (pub_on_root) then
               write(stdout,'(a)') &
                    'WARNING: slope along search direction is still positive.'
               write(stdout,'(a)') '         Reversing search direction!!'
            end if
            call data_functions_scale(direction_on_grid, -1.0_DP)
            call data_coef_scale(G_init, -1.0_DP)
            ! ndmh: if searching 'uphill', always re-check final step
            ! ndmh: before accepting, to avoid accepting very bad steps.
            reversing = .true.
         end if

      end if

      call services_flush

    end subroutine FO_find_direction


  end subroutine calc_FO_NGWFs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine calc_FO_ham(mdl, val_rep, ngwf_basis, direction, weight, dij, &
       FO_dij, full_GS_density, denskern, rho_ij_dproj, rho_ij, &
       response_kernel_jv, FO_lpseudo_fine, lhxc_fine, FO_lhxc_fine, &
       fxc_buffer, FO_ham_jv, joint_rep, joint_basis, proj_basis, &
       nl_projectors, tcore_FO_density)

    ! ========================================================================!
    ! This subroutine calculates the first-order hamiltonian (excluding the   !
    ! fixed \sum_ij d\tilde{p}^i/d\lambda D_ij \tilde{p}^j + \tilde{p}^i D_ij !
    ! d\tilde{p}^j term, which is calculated beforehand for efficiency        !
    ! ------------------------------------------------------------------------!
    ! As output, one has the first-order hamiltonian in j+v and v+v reps, the !
    ! first-order lhxc potential, the first-order xc potential, and the first-!
    ! order D_ij matrix                                                       !
    ! ------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in 2016                               !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,        !
    ! October 2018                                                            !
    ! ========================================================================!

    use augmentation, only: aug_FO_density_on_grid,augmentation_FO_screen_dij, &
         augmentation_density_on_grid, augmentation_overlap
    use datatypes, only: FUNCTIONS
    use comms, only: comms_barrier,pub_on_root,pub_my_proc_id
    use constants, only: DP, stdout, paw_en_size
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_brappd_ketppd
    use integrals, only: integrals_locpot
    use hartree, only: hartree_on_grid
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use paw, only:   paw_nonlocal_energies, &
         paw_tcore_density, paw_projector_overlap
    use rundat, only: pub_num_spins, pub_debug_on_root,pub_spin_fac, pub_aug, &
         pub_aug_den_dim, pub_paw, pub_nlcc, pub_any_nl_proj, pub_nhat_in_xc
    use projectors, only: PROJECTOR_SET
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_scale, &
         sparse_embed_create, sparse_embed_destroy, sparse_embed_axpy, &
         sparse_embed_copy, sparse_embed_transpose, sparse_embed_transpose_structure,&
         sparse_embed_product, sparse_embed_trace, sparse_embed_conjugate, &
         sparse_embed_array_scale, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use sparse, only: SPAM3, sparse_get_element, sparse_create, sparse_copy, &
         sparse_destroy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_abort

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    type(NGWF_REP), intent(in) :: val_rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    integer, intent(in) :: direction(mdl%nat)
    real(kind=DP), intent(in) :: weight(mdl%nat)
    type(SPAM3_EMBED), intent(in) :: dij(pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: FO_dij(pub_num_spins)
    real(kind=DP), intent(inout) :: full_GS_density(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    type(SPAM3_EMBED_ARRAY), intent(inout) :: denskern
    type(SPAM3_EMBED), intent(in) :: rho_ij_dproj(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: rho_ij(pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: response_kernel_jv(pub_num_spins)
    real(kind=DP), intent(in) ::FO_lpseudo_fine(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12)
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(inout) :: FO_lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(inout) :: fxc_buffer(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: FO_ham_jv
    type(NGWF_REP), intent(in) :: joint_rep
    type(FUNC_BASIS), intent(in) :: joint_basis(1)
    type(FUNC_BASIS), intent(in) :: proj_basis(1)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(1)
    real(kind=DP), intent(in) :: tcore_FO_density(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12)

    ! Local Variables

    real(kind=DP), allocatable :: nhat_GS_den_grad(:,:,:,:,:)
    real(kind=DP), allocatable :: nhat_FO_den_grad(:,:,:,:,:)

    real(kind=DP), dimension(:,:,:,:), allocatable :: dens_fine_FO
    real(kind=DP), dimension(:,:,:,:), allocatable :: FO_hartree_fine

    type(SPAM3), allocatable :: kern_array(:)

    type(SPAM3_EMBED), allocatable :: FO_dijhat(:)
    type(SPAM3_EMBED), allocatable :: FO_dijps(:)

    integer :: is, i1, i2, i3, ierr, iter1, iter2, iter1_eff

    real(kind=DP), allocatable :: dummy_sphere_energies(:)
    type(SPAM3), allocatable :: dij_tmp(:), rhoij_tmp(:), rhoij_FO_tmp(:)

    type(SPAM3_EMBED) :: temp_mat_1
    type(SPAM3_EMBED) :: rho_ij_FO(pub_num_spins)
    type(SPAM3_EMBED) :: pv_overlap, pj_overlap

    real(kind=DP) :: element

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering calc_FO_ham'

    ! Start timer
    call timer_clock('calc_FO_ham',1)

    ! gcc32: scale denskern to correct number of electrons
    if(pub_num_spins == 1) call sparse_embed_array_scale(denskern, 2.0_DP)

    if (pub_paw) then
       allocate(dummy_sphere_energies(paw_en_size), stat=ierr)
       call utils_alloc_check('calc_FO_ham', 'dummy_sphere_energies', ierr)
       allocate(FO_dijhat(pub_num_spins), stat=ierr)
       call utils_alloc_check('calc_FO_ham', 'FO_dijhat', ierr)
       allocate(FO_dijps(pub_num_spins), stat=ierr)
       call utils_alloc_check('calc_FO_ham', 'FO_dijps', ierr)
    end if

    allocate(dens_fine_FO(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
    call utils_alloc_check('calc_FO_ham', 'dens_fine_FO', ierr)

    allocate(FO_hartree_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
    call utils_alloc_check('calc_FO_ham', 'FO_hartree_fine', ierr)

    ! create and reset temporary matrix for  hamiltonian components
    call sparse_embed_scale(FO_ham_jv, 0.0_DP)

    ! gcc32: scale response kernel according to spin degeneracy
    if (pub_num_spins == 1) then
       call sparse_embed_scale(response_kernel_jv(1), 2.0_DP)
    end if

    if (pub_paw.or.pub_any_nl_proj) then
       call sparse_embed_transpose_structure(pv_overlap%structure, val_rep%sp_overlap)
       call sparse_embed_create(pv_overlap)
       call sparse_embed_transpose(pv_overlap, val_rep%sp_overlap)

       call sparse_embed_transpose_structure(pj_overlap%structure, &
            joint_rep%sp_overlap)
       call sparse_embed_create(pj_overlap)
       call sparse_embed_transpose(pj_overlap, joint_rep%sp_overlap)
    end if

    if (pub_paw) then
       do is = 1,pub_num_spins
          call sparse_embed_create(rho_ij_FO(is), rho_ij(is))
       end do

       do is = 1,pub_num_spins
          call sparse_embed_create(FO_dijhat(is), dij(is))
          call sparse_embed_create(FO_dijps(is), FO_dij(is), pv_overlap)
       end do
    end if

    ! calculate FO density
    call sparse_embed_transpose_structure(temp_mat_1%structure, &
         joint_rep%cross_overlap)
    call sparse_embed_create(temp_mat_1)
    call sparse_embed_transpose(temp_mat_1, joint_rep%cross_overlap)

    dens_fine_FO = 0.0_DP
    allocate(kern_array(pub_num_spins),stat=ierr)
    call utils_alloc_check('calc_FO_ham','kern_array',ierr)
    call sparse_embed_extract_from_array(kern_array,response_kernel_jv)
    call density_on_grid(dens_fine_FO, mdl%fine_grid, mdl%dbl_grid, mdl%cell, &
         mdl%fftbox, kern_array, temp_mat_1%p, joint_rep%ngwfs_on_grid(1), &
         joint_basis(1), val_rep%ngwfs_on_grid(1), ngwf_basis(1))
    call sparse_embed_destroy_extracted_array(kern_array)

    call sparse_embed_destroy(temp_mat_1)

    ! multiply by 2 due to symmetry
    do is = 1,pub_num_spins
       do i1 = 1, mdl%fine_grid%max_slabs12
          do i2 = 1, mdl%fine_grid%ld2
             do i3 = 1,mdl%fine_grid%ld1
               dens_fine_FO(i3,i2,i1,is) = 2.0_DP * dens_fine_FO(i3,i2,i1,is)
             end do
          end do
       end do
    end do

    if (pub_paw) then
       ! now, obtain the first-order rho_ij
       do is=1,pub_num_spins

          ! from respkern rv
          call sparse_embed_create(temp_mat_1, response_kernel_jv(is), &
               val_rep%sp_overlap)
          call sparse_embed_product(temp_mat_1, response_kernel_jv(is), &
               val_rep%sp_overlap)
          call sparse_embed_product(rho_ij_FO(is), pj_overlap, temp_mat_1)

          call sparse_embed_scale(rho_ij_FO(is), 2.0_DP)

          call sparse_embed_destroy(temp_mat_1)

          ! also, add the dproj contribution
          call sparse_embed_axpy(rho_ij_FO(is), rho_ij_dproj(is), 1.0_DP)

       end do

    end if

    if (pub_aug) then
       ! allocate and create data structures for nhat_GS_den_grad
       allocate(nhat_GS_den_grad(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12, pub_num_spins,0:pub_aug_den_dim), &
            stat=ierr)
       call utils_alloc_check('calc_FO_ham', 'nhat_GS_den_grad', ierr)

       allocate(nhat_FO_den_grad(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12, pub_num_spins, 0:pub_aug_den_dim), &
            stat=ierr)
       call utils_alloc_check('calc_FO_ham', 'nhat_FO_den_grad', ierr)

       allocate(rhoij_FO_tmp(pub_num_spins),stat=ierr)
       call utils_alloc_check('calc_FO_ham','rhoij_FO_tmp',ierr)

       nhat_GS_den_grad = 0.0_DP
       nhat_FO_den_grad = 0.0_DP

       call sparse_embed_extract_from_array(kern_array,denskern%m(:,1))
       call sparse_embed_extract_from_array(rhoij_FO_tmp,rho_ij_FO)

       ! now, obtain the GROUND-STATE augmentation(compensation) density on grid
       call augmentation_density_on_grid(nhat_GS_den_grad, mdl%fine_grid, &
            mdl%cell, mdl%regions(1)%pseudo_sp, mdl%regions(1)%paw_sp, &
            mdl%aug_box, kern_array, val_rep%sp_overlap%p)

       ! now, obtain the FIRST-ORDER augmentation(compensation) density on grid
       call aug_FO_density_on_grid(nhat_FO_den_grad, mdl%aug_box, &
            mdl%fine_grid, mdl%cell, kern_array, val_rep%sp_overlap%p, &
            rhoij_FO_tmp, direction, weight, mdl%regions(1)%paw_sp)

       call sparse_embed_destroy_extracted_array(kern_array)
       call sparse_embed_destroy_extracted_array(rhoij_FO_tmp)

       deallocate(rhoij_FO_tmp,stat=ierr)
       call utils_dealloc_check('calc_FO_ham','rhoij_FO_tmp',ierr)

       ! add nhat FOdens to dens FOfine: IMPLICITLY ASSUME NHAT included in fxc
       if (.not. pub_nhat_in_xc) call utils_abort('Implicit assumption FALSE')

       do is = 1,pub_num_spins
          do i1 = 1, mdl%fine_grid%max_slabs12
             do i2 = 1, mdl%fine_grid%ld2
                do i3 = 1,mdl%fine_grid%ld1
                   dens_fine_FO(i3,i2,i1,is) = dens_fine_FO(i3,i2,i1,is) + &
                        nhat_FO_den_grad(i3,i2,i1,is,0)
                end do
             end do
          end do
       end do


    end if

    FO_hartree_fine = 0.0_DP
    fxc_buffer = 0.0_DP

    ! calculate the FO hartree pot
    call hartree_on_grid(FO_hartree_fine, dens_fine_FO, mdl%fine_grid, mdl%cell)

    ! add tcore FO dens to dens FO fine
    if (pub_nlcc) then
       do is = 1,pub_num_spins
          do i1 = 1, mdl%fine_grid%max_slabs12
             do i2 = 1, mdl%fine_grid%ld2
                do i3 = 1,mdl%fine_grid%ld1
                   dens_fine_FO(i3,i2,i1,is) = dens_fine_FO(i3,i2,i1,is) + &
                        tcore_FO_density(i3,i2,i1) * 0.5_DP * pub_spin_fac
                end do
             end do
          end do
       end do
    end if


    if(pub_aug) then
       ! jcap: no emft potential here, as EMFT doesn't work yet with PAW
       call phonons_xc_finite_diff(fxc_buffer, dens_fine_FO, &
            full_GS_density, mdl%fine_grid, mdl%cell, nhat_FO_den_grad, &
            nhat_GS_den_grad)
    else
       call phonons_xc_finite_diff(fxc_buffer, dens_fine_FO, &
            full_GS_density, mdl%fine_grid, mdl%cell)
    endif

    ! add the previously calculated hartree component
    do is = 1,pub_num_spins
       do i1 = 1, mdl%fine_grid%max_slabs12
          do i2 = 1, mdl%fine_grid%ld2
             do i3 = 1, mdl%fine_grid%ld1
                FO_lhxc_fine(i3,i2,i1,is) = fxc_buffer(i3,i2,i1,is) + &
                     FO_hartree_fine(i3,i2,i1,is)
             end do
          end do
       end do
    end do

    ! add the local FO pseudo potential
    do is = 1,pub_num_spins
       do i1 = 1, mdl%fine_grid%max_slabs12
          do i2 = 1, mdl%fine_grid%ld2
             do i3 = 1, mdl%fine_grid%ld1
                FO_lhxc_fine(i3,i2,i1,is) = FO_lhxc_fine(i3,i2,i1,is) + &
                   FO_lpseudo_fine(i3,i2,i1)
             end do
          end do
       end do
    end do

    ! calculate the corresponding matrix term
    call sparse_embed_create(temp_mat_1, FO_ham_jv)
    call integrals_locpot(temp_mat_1%p, joint_rep%ngwfs_on_grid(1), &
         joint_basis(1), val_rep%ngwfs_on_grid(1), ngwf_basis(1), mdl%fine_grid, &
         mdl%dbl_grid, mdl%cell, mdl%fftbox, FO_lhxc_fine)

    ! add temporary hamiltonian to full FO hamiltonian
    call sparse_embed_axpy(FO_ham_jv, temp_mat_1, 1.0_DP)

    call sparse_embed_scale(temp_mat_1, 0.0_DP)

    if (pub_paw) then

       ! reset dij matrices
       do is = 1,pub_num_spins
          call sparse_embed_scale(FO_dijhat(is), 0.0_DP)
          call sparse_embed_scale(FO_dij(is), 0.0_DP)
          call sparse_embed_scale(FO_dijps(is), 0.0_DP)
       end do

       ! calculate FO dijhat
       call sparse_embed_extract_from_array(kern_array,FO_dijhat)
       call augmentation_FO_screen_dij(kern_array, lhxc_fine, FO_lhxc_fine, &
            mdl%aug_box, mdl%cell, mdl%fine_grid, direction, weight, &
            mdl%regions(1)%paw_sp, only_FO_locpot = .false., only_FO_Qs = .false., &
            only_SO_Qs =.false.)
       call sparse_embed_destroy_extracted_array(kern_array,FO_dijhat,.true.)

       ! calculate FO dij
       dummy_sphere_energies = 0.0_DP

       allocate(dij_tmp(pub_num_spins),stat=ierr)
       call utils_alloc_check('calc_FO_ham','dij_tmp',ierr)
       allocate(rhoij_tmp(pub_num_spins),stat=ierr)
       call utils_alloc_check('calc_FO_ham','rhoij_tmp',ierr)
       allocate(rhoij_FO_tmp(pub_num_spins),stat=ierr)
       call utils_alloc_check('calc_FO_ham','rhoij_FO_tmp',ierr)
       call sparse_embed_extract_from_array(dij_tmp,FO_dij)
       call sparse_embed_extract_from_array(rhoij_tmp,rho_ij)
       call sparse_embed_extract_from_array(rhoij_FO_tmp,rho_ij_FO)

       call paw_nonlocal_energies(dij_tmp, rhoij_tmp, mdl%regions(1)%paw_sp, &
            mdl%par, dummy_sphere_energies,.false., .true., rhoij_FO_tmp)

       call sparse_embed_destroy_extracted_array(dij_tmp,FO_dij,.true.)
       call sparse_embed_destroy_extracted_array(rhoij_tmp)
       call sparse_embed_destroy_extracted_array(rhoij_FO_tmp)

       ! add FO dijhat to FO dij
       do is = 1,pub_num_spins
          call sparse_embed_axpy(FO_dij(is), FO_dijhat(is), 1.0_DP)
          call sparse_embed_product(FO_dijps(is), FO_dij(is), pv_overlap)
       end do

       call sparse_embed_product(temp_mat_1, joint_rep%sp_overlap, FO_dijps(1))

       do is = 1,pub_num_spins
          call sparse_embed_destroy(FO_dijps(is))
          call sparse_embed_destroy(FO_dijhat(is))
       end do

       call sparse_embed_axpy(FO_ham_jv, temp_mat_1, 1.0_DP)

    end if ! if pub_paw

    ! destroy matrices
    call sparse_embed_destroy(temp_mat_1)

    deallocate(kern_array,stat=ierr)
    call utils_dealloc_check('calc_FO_ham','kern_array',ierr)

    ! gcc32: scale response kernel back to half-number of electrons
    if (pub_num_spins == 1) then
       call sparse_embed_scale(response_kernel_jv(1), 0.5_DP)
    end if

    if (pub_paw.or.pub_any_nl_proj) then
       call sparse_embed_destroy(pv_overlap)
       call sparse_embed_destroy(pj_overlap)
    end if

    if (pub_paw) then
       deallocate(dummy_sphere_energies, stat=ierr)
       call utils_dealloc_check('calc_FO_ham', 'dummy_sphere_energies', ierr)
       deallocate(FO_dijhat, stat=ierr)
       call utils_dealloc_check('calc_FO_ham', 'FO_dijhat', ierr)
       deallocate(FO_dijps, stat=ierr)
       call utils_dealloc_check('calc_FO_ham', 'FO_dijps', ierr)
    end if

    deallocate(dens_fine_FO, stat=ierr)
    call utils_dealloc_check('calc_FO_ham', 'dens_fine_FO', ierr)
    deallocate(FO_hartree_fine, stat=ierr)
    call utils_dealloc_check('calc_FO_ham', 'FO_hartree_fine', ierr)

    if (pub_aug) then
       deallocate(nhat_GS_den_grad, stat=ierr)
       call utils_dealloc_check('calc_FO_ham', 'nhat_GS_den_grad', ierr)
       deallocate(nhat_FO_den_grad, stat=ierr)
       call utils_dealloc_check('calc_FO_ham', 'nhat_FO_den_grad', ierr)
    end if

    ! gcc32: scale back denskern to half number of electrons
    if(pub_num_spins == 1) call sparse_embed_array_scale(denskern, 0.5_DP)

    ! pdh: sync procs
    call comms_barrier

    ! Stop timer
    call timer_clock('calc_FO_ham',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving calc_FO_ham'

  end subroutine calc_FO_ham


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! ========================================================================!
    ! This subroutine calculates the 'B' operator in the j+v rep              !
    ! ------------------------------------------------------------------------!
    ! For details regarding this operator, read the PhD thesis of Gabriel     !
    ! Constantinescu                                                          !
    ! ------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in 2016                               !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,        !
    ! October 2018                                                            !
    ! ========================================================================!

  subroutine calc_B_op(mdl, joint_rep, joint_basis, val_rep, ngwf_basis, &
       direction, weight, dij, B_op_nv, iGps_overlap, newsiGp_overlap, &
       lhxc_fine, locpot_FO_tchart)

    use augmentation, only: augmentation_FO_screen_dij
    use datatypes, only: FUNCTIONS
    use comms, only: comms_barrier, pub_on_root, pub_my_proc_id
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_locpot
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use paw, only:  paw_tcore_hartree_on_grid
    use rundat, only: pub_num_spins, pub_debug_on_root, pub_spin_fac, pub_aug, &
         pub_aug_den_dim, pub_paw, pub_any_nl_proj
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_scale, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_axpy, sparse_embed_copy, &
         sparse_embed_transpose, sparse_embed_transpose_structure, &
         sparse_embed_product, sparse_embed_trace, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use sparse, only: SPAM3, sparse_create, sparse_copy, sparse_destroy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    type(NGWF_REP), intent(in) :: joint_rep
    type(FUNC_BASIS), intent(in) :: joint_basis(1)
    type(NGWF_REP), intent(in) :: val_rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    integer, intent(in) :: direction(mdl%nat)
    real(kind=DP), intent(in) :: weight(mdl%nat)
    type(SPAM3_EMBED), intent(in) :: dij(pub_num_spins)
    type(SPAM3_EMBED),intent(inout) :: B_op_nv
    type(SPAM3_EMBED), intent(in) :: iGps_overlap     ! iG-epsilon
    type(SPAM3_EMBED), intent(in) :: newsiGp_overlap
    real(kind=DP), intent(in) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in) :: locpot_FO_tchart(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12) ! d/d_eps locpot

    ! Local Variables
    integer :: is, ierr
    type(SPAM3_EMBED), allocatable :: dij_FO(:)

    type(SPAM3_EMBED) :: newsp_dij_FO, newspdij, newsiGpdij, temp_mat_1, temp_mat_2, &
         ps_overlap, B_op_nv_temp

    real(kind=DP) :: dummy_gzero
    type(SPAM3), allocatable :: dij_tmp(:)


    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering calc_B_op'

    ! Start timer
    call timer_clock('calc_B_op',1)

    if (pub_paw) then
       allocate(dij_FO(pub_num_spins), stat=ierr)
       call utils_alloc_check('calc_B_op', 'dij_FO', ierr)

       ! create dij_FO matrices
       do is = 1,pub_num_spins
          call sparse_embed_create(dij_FO(is), dij(is))
          call sparse_embed_scale(dij_FO(is), 0.0_DP)
       end do
    end if

    ! create and reset temporary B_op matrices
    call sparse_embed_create(B_op_nv_temp, B_op_nv)

    call sparse_embed_scale(B_op_nv_temp, 0.0_DP)
    call sparse_embed_scale(B_op_nv, 0.0_DP)

    if (pub_paw.or.pub_any_nl_proj) then
       ! create ps overlap in val NGWF rep
       call sparse_embed_transpose_structure(ps_overlap%structure, val_rep%sp_overlap)
       call sparse_embed_create(ps_overlap)
       call sparse_embed_transpose(ps_overlap, val_rep%sp_overlap)

       ! ========== calculate d/d_eps(|pa><pB|) Dab terms in both NGWF reps ===
       call sparse_embed_create(newsiGpdij, newsiGp_overlap, dij(1))
       call sparse_embed_product(newsiGpdij, newsiGp_overlap, dij(1))

       call sparse_embed_create(newspdij, joint_rep%sp_overlap, dij(1))
       call sparse_embed_product(newspdij, joint_rep%sp_overlap, dij(1))

       call sparse_embed_product(B_op_nv_temp, newsiGpdij, ps_overlap)
       call sparse_embed_destroy(newsiGpdij)

       call sparse_embed_create(temp_mat_1, newspdij, iGps_overlap)
       call sparse_embed_product(temp_mat_1, newspdij, iGps_overlap)
       call sparse_embed_axpy(B_op_nv_temp, temp_mat_1, 1.0_DP)
       call sparse_embed_destroy(temp_mat_1)
       call sparse_embed_destroy(newspdij)

       ! copy to B_op main matrices
       call sparse_embed_copy(B_op_nv, B_op_nv_temp)
       call sparse_embed_scale(B_op_nv_temp,0.0_DP)
       ! ======================================================================

    end if ! if pub_paw

    ! =========== calculate the VhtrhoZc matrix terms in both NGWF reps =======
    ! calculate the FO hartree potential of the tcore density on grid
    call integrals_locpot(B_op_nv_temp%p, joint_rep%ngwfs_on_grid(1), &
         joint_basis(1), val_rep%ngwfs_on_grid(1), ngwf_basis(1), mdl%fine_grid, &
         mdl%dbl_grid, mdl%cell, mdl%fftbox, locpot_FO_tchart)

    ! add to B_op main matrices
    call sparse_embed_axpy(B_op_nv, B_op_nv_temp,1.0_DP)
    call sparse_embed_scale(B_op_nv_temp,0.0_DP)
    ! =========================================================================

    ! ========== calculate the proj_Qs term in both NGWF reps =================
    if (pub_aug) then

       allocate(dij_tmp(pub_num_spins),stat=ierr)
       call utils_alloc_check('calc_B_op','dij_tmp',ierr)
       call sparse_embed_extract_from_array(dij_tmp,dij_FO)
       call augmentation_FO_screen_dij(dij_tmp, lhxc_fine, locpot_FO_tchart, &
            mdl%aug_box, mdl%cell, mdl%fine_grid, direction, weight, &
            mdl%regions(1)%paw_sp, only_FO_locpot = .false., only_FO_Qs = .false., &
            only_SO_Qs = .false.)
       call sparse_embed_destroy_extracted_array(dij_tmp,dij_FO,.true.)

       deallocate(dij_tmp,stat=ierr)
       call utils_dealloc_check('calc_B_op','dij_tmp',ierr)

       call sparse_embed_create(newsp_dij_FO, joint_rep%sp_overlap, dij_FO(1))
       call sparse_embed_product(newsp_dij_FO, joint_rep%sp_overlap, dij_FO(1))
       call sparse_embed_product(B_op_nv_temp, newsp_dij_FO, ps_overlap)
       call sparse_embed_destroy(newsp_dij_FO)

       ! add to B_op main matrices
       call sparse_embed_axpy(B_op_nv, B_op_nv_temp,1.0_DP)

    end if
    ! =========================================================================

    call sparse_embed_destroy(B_op_nv_temp)

    if (pub_paw.or.pub_any_nl_proj) then
       call sparse_embed_destroy(ps_overlap)
    end if

    if (pub_paw) then
       ! destroy dij_FO matrices
       do is = 1,pub_num_spins
          call sparse_embed_destroy(dij_FO(is))
       end do

       deallocate(dij_FO, stat=ierr)
       call utils_alloc_check('calc_B_op', 'dij_FO', ierr)
    end if

    ! pdh: sync procs
    call comms_barrier

    ! Stop timer
    call timer_clock('calc_B_op',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving calc_B_op'

  end subroutine calc_B_op

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_A_op(mdl, val_rep, ngwf_basis, direction1, weight1, &
       direction2, weight2, dij, FO_dij, A_op, lhxc_fine, FO_lhxc_fine, &
       locpot_FO_tchart, iG1ps_overlap, iG2ps_overlap, iG1iG2ps_overlap, &
       locpot_SO_mat)

    ! ========================================================================!
    ! This subroutine calculates the 'A' operator in the v+v rep              !
    ! ------------------------------------------------------------------------!
    ! For details regarding this operator, read the PhD thesis of Gabriel     !
    ! Constantinescu                                                          !
    ! ------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in 2016                               !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,        !
    ! October 2018                                                            !
    ! ========================================================================!

    use augmentation, only: augmentation_FO_screen_dij
    use datatypes, only: FUNCTIONS
    use comms, only: comms_barrier, pub_on_root, pub_my_proc_id
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_locpot
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use paw, only: paw_tcore_hartree_on_grid
    use pseudopotentials, only: pseudopotentials_local_on_grid, &
         pseudo_make_structure_factor
    use rundat, only: pub_num_spins, pub_debug_on_root, pub_spin_fac, pub_aug, &
         pub_aug_den_dim, pub_paw, pub_any_nl_proj
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_scale, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_axpy, sparse_embed_transpose, &
         sparse_embed_transpose_structure, sparse_embed_product, sparse_embed_trace, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use sparse, only: SPAM3, sparse_create, sparse_copy, sparse_destroy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    type(NGWF_REP), intent(in) :: val_rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    integer, intent(in) :: direction1(mdl%nat)
    real(kind=DP), intent(in) :: weight1(mdl%nat)
    integer, intent(in) :: direction2(mdl%nat)
    real(kind=DP), intent(in) :: weight2(mdl%nat)
    type(SPAM3_EMBED), intent(in) :: dij(pub_num_spins) ! D^1+\tilde{D}^1+\hat{D}
    type(SPAM3_EMBED), intent(in) :: FO_dij(pub_num_spins) ! FO w.r.t. lambda
    type(SPAM3_EMBED),intent(inout) :: A_op
    real(kind=DP), intent(in) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in) :: FO_lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in) :: locpot_FO_tchart(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12) ! (d/d_eps) locpot
    type(SPAM3_EMBED), intent(in) :: iG1ps_overlap     ! iG_lambda
    type(SPAM3_EMBED), intent(in) :: iG2ps_overlap     ! iG_epsilon
    type(SPAM3_EMBED), intent(in) :: iG1iG2ps_overlap     ! -G_lambda*G_epsilon
    type(SPAM3_EMBED), intent(in), optional :: locpot_SO_mat

    ! Local Variables
    integer :: is, ierr
    type(SPAM3_EMBED), allocatable :: dij_FO_temp_1(:), dij_FO_temp_2(:)

    ! \frac{d^2}{d\epsilon~d\lambda} locpot_tcore_hartree
    real(kind=DP), allocatable :: locpot_SO_tchart(:,:,:)

    type(SPAM3_EMBED) :: siG1p_overlap, siG2p_overlap, siG1iG2p_overlap

    complex(kind=DP), allocatable, dimension(:,:,:,:) :: struct_fac
    complex(kind=DP), allocatable, dimension(:,:,:) :: struct_fac_classical

    real(kind=DP) :: dummy_gzero

    type(SPAM3_EMBED) :: temp_mat_1
    type(SPAM3_EMBED) :: ps_overlap
    type(SPAM3_EMBED) :: A_op_temp
    type(SPAM3), allocatable :: dij_tmp(:)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering calc_A_op'

    ! Start timer
    call timer_clock('calc_A_op',1)

    allocate(locpot_SO_tchart(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12), stat=ierr)
    call utils_alloc_check('calc_A_op','locpot_SO_tchart', ierr)

    if (pub_paw) then

       allocate(dij_FO_temp_1(pub_num_spins), stat=ierr)
       call utils_alloc_check('calc_A_op', 'dij_FO_temp_1', ierr)
       allocate(dij_FO_temp_2(pub_num_spins), stat=ierr)
       call utils_alloc_check('calc_A_op', 'dij_FO_temp_2', ierr)

       ! create dij_FO_temp matrices
       do is = 1,pub_num_spins
          call sparse_embed_create(dij_FO_temp_1(is), dij(is))
          call sparse_embed_create(dij_FO_temp_2(is), dij(is))
       end do

    end if

    ! create and reset a temporary A_op
    call sparse_embed_create(A_op_temp, A_op)
    call sparse_embed_scale(A_op_temp, 0.0_DP)
    call sparse_embed_scale(A_op,0.0_DP)

    if (pub_paw.or.pub_any_nl_proj) then
       ! create temp_mat_1
       call sparse_embed_create(temp_mat_1,val_rep%sp_overlap,dij(1))

       ! create ps overlap in val NGWF rep
       call sparse_embed_transpose_structure(ps_overlap%structure, val_rep%sp_overlap)
       call sparse_embed_create(ps_overlap)
       call sparse_embed_transpose(ps_overlap, val_rep%sp_overlap)

       ! create transpose for iG1ps_overlap, iG2ps_overlap and iG1iG2ps_overlap
       call sparse_embed_transpose_structure(siG1p_overlap%structure, iG1ps_overlap)
       call sparse_embed_create(siG1p_overlap)
       call sparse_embed_transpose(siG1p_overlap, iG1ps_overlap)

       call sparse_embed_transpose_structure(siG2p_overlap%structure, iG2ps_overlap)
       call sparse_embed_create(siG2p_overlap)
       call sparse_embed_transpose(siG2p_overlap, iG2ps_overlap)

       call sparse_embed_transpose_structure(siG1iG2p_overlap%structure, &
            iG1iG2ps_overlap)
       call sparse_embed_create(siG1iG2p_overlap)
       call sparse_embed_transpose(siG1iG2p_overlap, iG1iG2ps_overlap)

    end if

    ! =========== calculate the second-order tcore_hartree matrix terms ====
    ! reset locpot_SO
    locpot_SO_tchart = 0.0_DP

    allocate(struct_fac(mdl%num_pspecies, mdl%fine_grid%ld3, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs23),stat=ierr)
    call utils_alloc_check('calc_A_op','struct_fac',ierr)
    ! ars: calculate structure factor for the classical atoms
    allocate(struct_fac_classical(0,0,0),stat=ierr)
    call utils_alloc_check('calc_A_op', 'struct_fac_classical', ierr)

    ! calculate structure factor for all species
    call pseudo_make_structure_factor(struct_fac, mdl%elements, &
         mdl%fine_grid, mdl%par%num_pspecies, direction1, weight1, &
         direction2, weight2)

    ! calculate the SO hartree potential of the tcore density on grid
    if (pub_paw) then
       call paw_tcore_hartree_on_grid(locpot_SO_tchart, dummy_gzero, &
            struct_fac, struct_fac_classical, mdl%regions(1)%paw_sp, &
            mdl%fine_grid, mdl%par)
    else
       call pseudopotentials_local_on_grid(locpot_SO_tchart, dummy_gzero, &
            struct_fac, struct_fac_classical, mdl%fine_grid, mdl%cell, &
            mdl%regions(1)%pseudo_sp, mdl%par, mdl%par%nat_classical)
    end if


    deallocate(struct_fac,stat=ierr)
    call utils_alloc_check('calc_A_op','struct_fac',ierr)

    deallocate(struct_fac_classical,stat=ierr)
    call utils_dealloc_check('calc_A_op','struct_fac_classical', &
         ierr)

    if (.not.present(locpot_SO_mat)) then
       call integrals_locpot(A_op_temp%p, val_rep%ngwfs_on_grid(1), &
            ngwf_basis(1), val_rep%ngwfs_on_grid(1), ngwf_basis(1), mdl%fine_grid, &
            mdl%dbl_grid, mdl%cell, mdl%fftbox, locpot_SO_tchart)

       call sparse_embed_axpy(A_op,A_op_temp, 1.0_DP) ! copy to full A_op
       call sparse_embed_scale(A_op_temp, 0.0_DP)
    else
       call sparse_embed_axpy(A_op, locpot_SO_mat, 1.0_DP)
    end if
    ! =========================================================================


    if (pub_paw.or.pub_any_nl_proj) then

       ! ========== calculate d2/d_eps d_lam (|pa><pB|) Dab terms  ============

       ! |pa> D_ab [ \frac{d^2}{d\epsilon~d\lambda} <pb| ]
       call sparse_embed_product(temp_mat_1, val_rep%sp_overlap, dij(1))
       call sparse_embed_product(A_op_temp, temp_mat_1, iG1iG2ps_overlap)

       call sparse_embed_axpy(A_op, A_op_temp, 1.0_DP)
       call sparse_embed_scale(A_op_temp, 0.0_DP)

       ! [\frac{d^2}{d\epsilon~d\lambda} |pa> ] D_ab <pb|
       call sparse_embed_product(temp_mat_1, siG1iG2p_overlap, dij(1))
       call sparse_embed_product(A_op_temp, temp_mat_1, ps_overlap)

       call sparse_embed_axpy(A_op, A_op_temp, 1.0_DP)
       call sparse_embed_scale(A_op_temp, 0.0_DP)

       ! [\frac{d}{d\epsilon} |pa> ] D_ab [\frac{d}{d\lambda} <pb| ]
       call sparse_embed_product(temp_mat_1, siG2p_overlap, dij(1))
       call sparse_embed_product(A_op_temp, temp_mat_1, iG1ps_overlap)

       call sparse_embed_axpy(A_op, A_op_temp, 1.0_DP)
       call sparse_embed_scale(A_op_temp, 0.0_DP)

       ! [\frac{d}{d\lambda} |pa> ] D_ab [\frac{d}{d\epsilon} <pb| ]
       call sparse_embed_product(temp_mat_1, siG1p_overlap, dij(1))
       call sparse_embed_product(A_op_temp, temp_mat_1, iG2ps_overlap)

       call sparse_embed_axpy(A_op, A_op_temp, 1.0_DP)
       call sparse_embed_scale(A_op_temp, 0.0_DP)

       ! ======================================================================
    end if

    if (pub_paw) then

       ! ========== calculate d/d_lam (|pa><pb|) [\int lhxc d/d_eps
       !              \hat{Q}_ab dr + \int FO_locpot \hat{Q}_ab dr] terms  ====

       if (pub_aug) then

          do is = 1,pub_num_spins
             call sparse_embed_scale(dij_FO_temp_1(is), 0.0_DP)
          end do

          ! for \int lhxc(r) \frac{d}{d\epsilon} \hat{Q}_ab (r) dr +
          !      \int [\frac{d}{d\epsilon} locpot] \hat{Q}_ab (r) dr
          allocate(dij_tmp(pub_num_spins),stat=ierr)
          call utils_alloc_check('calc_A_op','dij_tmp',ierr)
          call sparse_embed_extract_from_array(dij_tmp,dij_FO_temp_1)
          call augmentation_FO_screen_dij(dij_tmp, lhxc_fine, &
               locpot_FO_tchart, mdl%aug_box, mdl%cell, mdl%fine_grid, &
               direction2, weight2, mdl%regions(1)%paw_sp, only_FO_locpot = .false., &
               only_FO_Qs = .false., only_SO_Qs = .false.)
          call sparse_embed_destroy_extracted_array(dij_tmp,dij_FO_temp_1,.true.)

          ! [\frac{d}{d\lambda} |pa> ] D_ab  <pb|
          call sparse_embed_product(temp_mat_1, siG1p_overlap, dij_FO_temp_1(1))
          call sparse_embed_product(A_op_temp, temp_mat_1, ps_overlap)

          call sparse_embed_axpy(A_op, A_op_temp, 1.0_DP)
          call sparse_embed_scale(A_op_temp, 0.0_DP)

          ! |pa> D_ab [\frac{d}{d\lambda} <pb| ]
          call sparse_embed_product(temp_mat_1, val_rep%sp_overlap, dij_FO_temp_1(1))
          call sparse_embed_product(A_op_temp, temp_mat_1, iG1ps_overlap)

          call sparse_embed_axpy(A_op, A_op_temp, 1.0_DP)
          call sparse_embed_scale(A_op_temp, 0.0_DP)

       end if
       ! ======================================================================

       ! ========== calculate d/d_eps (|pa><pb|) [FO_dij] terms  ==============

       if (pub_aug) then

          ! [\frac{d}{d\epsilon} |pa> ] FO_D_ab  <pb|
          call sparse_embed_product(temp_mat_1, siG2p_overlap, FO_dij(1))
          call sparse_embed_product(A_op_temp, temp_mat_1, ps_overlap)

          call sparse_embed_axpy(A_op, A_op_temp, 1.0_DP)
          call sparse_embed_scale(A_op_temp, 0.0_DP)

          ! |pa> FO_D_ab [\frac{d}{d\epsilon} <pb| ]
          call sparse_embed_product(temp_mat_1, val_rep%sp_overlap, FO_dij(1))
          call sparse_embed_product(A_op_temp, temp_mat_1, iG2ps_overlap)

          call sparse_embed_axpy(A_op, A_op_temp, 1.0_DP)
          call sparse_embed_scale(A_op_temp,0.0_DP)

       end if
       ! ======================================================================


       ! ========== calculate |pa><pb [...] terms  =================

       if (pub_aug) then

          do is = 1,pub_num_spins
             call sparse_embed_scale(dij_FO_temp_1(is), 0.0_DP)
             call sparse_embed_scale(dij_FO_temp_2(is), 0.0_DP)
          end do

          call sparse_embed_extract_from_array(dij_tmp,dij_FO_temp_1)
          ! for \int FO_lhxc(r) \frac{d}{d\eps} \hat{Q}_ab (r) dr
          call augmentation_FO_screen_dij(dij_tmp, FO_lhxc_fine, &
               locpot_FO_tchart, mdl%aug_box, mdl%cell, mdl%fine_grid, &
               direction2, weight2, mdl%regions(1)%paw_sp, only_FO_locpot = .false., &
               only_FO_Qs = .true., only_SO_Qs = .false.)
          call sparse_embed_destroy_extracted_array(dij_tmp,dij_FO_temp_1,.true.)
          call sparse_embed_extract_from_array(dij_tmp,dij_FO_temp_2)
          ! for \int lhxc(r) \frac{d^2}{d\eps~d\lam} \hat{Q}_ab (r) dr
          call augmentation_FO_screen_dij(dij_tmp, lhxc_fine, &
               locpot_FO_tchart, mdl%aug_box, mdl%cell, mdl%fine_grid, &
               direction1, weight1, mdl%regions(1)%paw_sp, only_FO_locpot = .false., &
               only_FO_Qs = .false., only_SO_Qs = .true., direction2 = &
               direction2, weight2 = weight2)
          call sparse_embed_destroy_extracted_array(dij_tmp,dij_FO_temp_2,.true.)

          do is = 1,pub_num_spins
             call sparse_embed_axpy(dij_FO_temp_1(is), dij_FO_temp_2(is), 1.0_DP)
             call sparse_embed_scale(dij_FO_temp_2(is), 0.0_DP)
          end do

          call sparse_embed_extract_from_array(dij_tmp,dij_FO_temp_2)
          ! for \int \frac{d}{d\eps} locpot(r) \frac{d}{d\lam} \hat{Q}_ab (r) dr
          call augmentation_FO_screen_dij(dij_tmp, locpot_FO_tchart, &
               locpot_FO_tchart, mdl%aug_box, mdl%cell, mdl%fine_grid, &
               direction1, weight1, mdl%regions(1)%paw_sp, only_FO_locpot = .false., &
               only_FO_Qs = .true., only_SO_Qs = .false.)
          call sparse_embed_destroy_extracted_array(dij_tmp,dij_FO_temp_2,.true.)

          do is = 1,pub_num_spins
             call sparse_embed_axpy(dij_FO_temp_1(is), dij_FO_temp_2(is), 1.0_DP)
             call sparse_embed_scale(dij_FO_temp_2(is), 0.0_DP)
          end do

          call sparse_embed_extract_from_array(dij_tmp,dij_FO_temp_2)
          ! for \int \frac{d^2}{d\eps~d\lam} locpot(r)  \hat{Q}_ab (r) dr
          call augmentation_FO_screen_dij(dij_tmp, lhxc_fine, &
               locpot_SO_tchart, mdl%aug_box, mdl%cell, mdl%fine_grid, &
               direction1, weight1, mdl%regions(1)%paw_sp, only_FO_locpot = .true., &
               only_FO_Qs = .false.,only_SO_Qs = .false.)
          call sparse_embed_destroy_extracted_array(dij_tmp,dij_FO_temp_2,.true.)

          do is = 1,pub_num_spins
             call sparse_embed_axpy(dij_FO_temp_1(is), dij_FO_temp_2(is), 1.0_DP)
             call sparse_embed_scale(dij_FO_temp_2(is), 0.0_DP)
          end do

          ! |pa> [...]  <pb|
          call sparse_embed_product(temp_mat_1, val_rep%sp_overlap, dij_FO_temp_1(1))
          call sparse_embed_product(A_op_temp, temp_mat_1, ps_overlap)
          call sparse_embed_destroy(temp_mat_1)

          call sparse_embed_axpy(A_op, A_op_temp, 1.0_DP)
          call sparse_embed_scale(A_op_temp, 0.0_DP)
       end if
       ! ======================================================================

    end if ! if pub_paw

    if (pub_paw.or.pub_any_nl_proj) then
       ! destroy transpose for iG1ps_overlap, iG2ps_overlap and iG1iG2ps_overlap
       call sparse_embed_destroy(siG1p_overlap)
       call sparse_embed_destroy(siG2p_overlap)
       call sparse_embed_destroy(siG1iG2p_overlap)
    end if

    call sparse_embed_destroy(A_op_temp)

    if (pub_paw.or.pub_any_nl_proj) then
        call sparse_embed_destroy(ps_overlap)
    end if

    if (pub_paw) then
       ! destroy dij_FO matrices
       do is = 1,pub_num_spins
          call sparse_embed_destroy(dij_FO_temp_1(is))
          call sparse_embed_destroy(dij_FO_temp_2(is))
       end do

       deallocate(dij_FO_temp_1, stat=ierr)
       call utils_alloc_check('calc_A_op', 'dij_FO_temp_1', ierr)
       deallocate(dij_FO_temp_2, stat=ierr)
       call utils_alloc_check('calc_A_op', 'dij_FO_temp_2', ierr)

    end if

    deallocate(locpot_SO_tchart, stat=ierr)
    call utils_dealloc_check('calc_A_op','locpot_SO_tchart', ierr)

    ! pdh: sync procs
    call comms_barrier

    ! Stop timer
    call timer_clock('calc_A_op',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving calc_A_op'

  end subroutine calc_A_op

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine second_order_energy(mdl, val_rep, ngwf_basis, val_ham, &
       joint_basis, joint_rep, joint_ham, directions, weights, &
       perturbation1, dij, FO_dij, FO_ham, FO_ham_jv, lhxc_fine, &
       FO_lhxc_fine, vxc_buffer, fxc_buffer, val_denskern, &
       response_denskern_jv, proj_basis, nl_projectors, iG1pv_overlap, &
       dipole_mat_jv, at_gamma, force_constant_matrix, only_SO_energy, &
       locpot_SO_mat)

    ! ========================================================================!
    ! This subroutine calculates the second-order energy with respect to both !
    ! atomic displacements and E-field perturbations. It can calculate the    !
    ! variational form (if only_SO_energy is present) or the non-variational  !
    ! form (if force_constant_matrix) is given                                !
    ! ------------------------------------------------------------------------!
    ! When the variational form is required, the second-order term with       !
    ! to a single perturbation is calculated, i.e. no cross-terms are computed!
    ! ------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in 2016                               !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,        !
    ! October 2018                                                            !
    ! ========================================================================!

    use augmentation, only: augmentation_FO_screen_dij
    use datatypes, only: FUNCTIONS
    use comms, only: comms_barrier, pub_on_root, pub_my_proc_id
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_kinetic, integrals_locpot
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use paw, only: paw_FO_nlcc_energy, paw_tcore_hartree_on_grid
    use pseudopotentials, only: pseudopotentials_FO_nlcc_energy, &
         pseudo_make_structure_factor, pseudopotentials_local_on_grid
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_num_spins, pub_debug_on_root, pub_spin_fac, pub_aug, &
         pub_aug_den_dim, pub_nlcc, pub_paw, pub_any_nl_proj
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_scale, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_axpy, sparse_embed_transpose, &
         sparse_embed_transpose_structure, sparse_embed_product, sparse_embed_trace, &
         sparse_embed_conjugate
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    use visual, only: visual_ngwfs, visual_ngwfs_radial

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    type(NGWF_REP), intent(in) :: val_rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    type(NGWF_HAM), intent(in) :: val_ham
    type(FUNC_BASIS), intent(in) :: joint_basis(1)
    type(NGWF_REP), intent(in) :: joint_rep
    type(NGWF_HAM), intent(in) :: joint_ham
    integer, intent(in) :: directions(3*mdl%nat+3, mdl%nat)
    real(kind=DP), intent(in) :: weights(3*mdl%nat+3, mdl%nat)
    integer, intent(in) :: perturbation1
    type(SPAM3_EMBED), intent(in) :: dij(pub_num_spins) ! D^1+\tilde{D}^1+\hat{D}
    type(SPAM3_EMBED), intent(in) :: FO_dij(pub_num_spins) ! FO w.r.t. lambda
    type(SPAM3_EMBED), intent(in) :: FO_ham, FO_ham_jv   ! FO w.r.t. lambda
    real(kind=DP), intent(in) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    ! FO lhxc potential w.r.t. lambda:
    real(kind=DP), intent(in) :: FO_lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in) :: vxc_buffer(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    ! FO vxc potential w.r.t. lambda:
    real(kind=DP), intent(in) :: fxc_buffer(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: val_denskern
    type(SPAM3_EMBED), intent(inout) :: response_denskern_jv
    type(FUNC_BASIS), intent(in) :: proj_basis(1)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(1)
    type(SPAM3_EMBED), intent(in) :: iG1pv_overlap
    type(SPAM3_EMBED), intent(in) :: dipole_mat_jv(3)
    logical, intent(in) :: at_gamma
    complex(kind=DP), intent(inout), optional :: force_constant_matrix( &
         mdl%nat*3+3,mdl%nat*3+3)
    real(kind=DP), intent(inout), optional :: only_SO_energy
    type(SPAM3_EMBED), intent(in), optional :: locpot_SO_mat

    ! Local Variables
    integer :: is, ierr
    type(SPAM3_EMBED) :: temp_mat_1, temp_mat_2, temp_mat_3
    type(SPAM3_EMBED) :: A_op, B_op_jv

    type(SPAM3_EMBED) :: pv_overlap
    type(SPAM3_EMBED) :: iG2pv_overlap, iG1iG2pv_overlap, dummy_iGpv_overlap
    type(SPAM3_EMBED) :: second_overlap ! lambda epsilon

    ! \frac{d}{d\epsilon} locpot_tcore_hartree
    real(kind=DP), allocatable :: FO_locpot_eps(:,:,:)
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: struct_fac
    complex(kind=DP), allocatable, dimension(:,:,:) :: struct_fac_classical

    ! epsilon only:
    type(SPAM3_EMBED) :: first_overlap_vv, first_overlap_jv

    type(SPAM3_EMBED) :: pj_overlap ! with respngwfs NGWFs w.r.t. lambda
    type(SPAM3_EMBED) :: jiG2p_overlap ! respngwfs NGWFs w.r.t. lambda and
    ! perturbed projectors w.r.t. epsilon
    type(SPAM3_EMBED) :: viG2p_overlap, ham_jv, response_denskern_jv_CT

    real(kind=DP) :: temp_energy, element, dummy_gzero
    integer :: perturbation2, perturbation2_min,perturbation2_max  ! epsilon

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &second_order_energy'

    ! Start timer
    call timer_clock('second_order_energy',1)

    ! pdh: sync procs
    call comms_barrier

    ! scale denskern and respkern to twice the number of electrons if
    ! degenerate spins
    if (pub_num_spins == 1) call sparse_embed_scale(val_denskern, 2.0_DP)
    if (pub_num_spins == 1) call sparse_embed_scale(response_denskern_jv, 2.0_DP)

    allocate(FO_locpot_eps(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12), stat=ierr)
    call utils_alloc_check('second_order_energy','FO_locpot_eps', ierr)

    allocate(struct_fac(mdl%num_pspecies, mdl%fine_grid%ld3, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs23),stat=ierr)
    call utils_alloc_check('second_order_energy','struct_fac',ierr)
    allocate(struct_fac_classical(0,0,0),stat=ierr)
    call utils_alloc_check('second_order_energy', 'struct_fac_classical', &
         ierr)


    call sparse_embed_transpose_structure(B_op_jv%structure,joint_rep%cross_overlap)
    call sparse_embed_create(B_op_jv)

    call sparse_embed_transpose_structure(response_denskern_jv_CT%structure, &
         response_denskern_jv)
    call sparse_embed_create(response_denskern_jv_CT)
    call sparse_embed_transpose(response_denskern_jv_CT,response_denskern_jv)
    if (response_denskern_jv%iscmplx) then
       call sparse_embed_conjugate(response_denskern_jv_CT)
    end if

    if (pub_paw.or.pub_any_nl_proj) then
       call sparse_embed_transpose_structure(pv_overlap%structure, val_rep%sp_overlap)
       call sparse_embed_create(pv_overlap)
       call sparse_embed_transpose(pv_overlap, val_rep%sp_overlap)

       call sparse_embed_create(iG2pv_overlap, pv_overlap)
       call sparse_embed_create(viG2p_overlap, val_rep%sp_overlap)
       call sparse_embed_create(dummy_iGpv_overlap, pv_overlap)
       call sparse_embed_create(iG1iG2pv_overlap, pv_overlap)

       call sparse_embed_create(jiG2p_overlap, joint_rep%sp_overlap)
       call sparse_embed_transpose_structure(pj_overlap%structure, &
            joint_rep%sp_overlap)
       call sparse_embed_create(pj_overlap)
       call sparse_embed_transpose(pj_overlap, joint_rep%sp_overlap)

    end if

    call sparse_embed_create(second_overlap, val_rep%overlap)
    call sparse_embed_create(first_overlap_vv, val_rep%overlap)
    call sparse_embed_transpose_structure(first_overlap_jv%structure, &
         joint_rep%cross_overlap)
    call sparse_embed_create(first_overlap_jv)

    call sparse_embed_create(A_op, val_ham%ham(1))

    ! prepare ham_jv only if the non-variational second-order energy is required
    if (.not.present(only_SO_energy)) then
       call sparse_embed_create(ham_jv, first_overlap_jv)
       call sparse_embed_scale(ham_jv, 0.0_DP)

       call integrals_kinetic(ham_jv%p, joint_rep%ngwfs_on_grid(1), joint_basis(1), &
            val_rep%ngwfs_on_grid(1), ngwf_basis(1), mdl%cell, mdl%fftbox)

       call sparse_embed_create(temp_mat_1, ham_jv)
       call integrals_locpot(temp_mat_1%p, joint_rep%ngwfs_on_grid(1), &
            joint_basis(1), val_rep%ngwfs_on_grid(1), ngwf_basis(1), mdl%fine_grid, &
            mdl%dbl_grid, mdl%cell, mdl%fftbox, lhxc_fine)
       call sparse_embed_axpy(ham_jv, temp_mat_1, 1.0_DP)
       call sparse_embed_destroy(temp_mat_1)

       if (pub_paw.or.pub_any_nl_proj) then
          call sparse_embed_create(temp_mat_1, ham_jv)
          call sparse_embed_create(temp_mat_2, joint_rep%sp_overlap, dij(1))
          call sparse_embed_product(temp_mat_2, joint_rep%sp_overlap, dij(1))
          call sparse_embed_product(temp_mat_1, temp_mat_2, pv_overlap)
          call sparse_embed_destroy(temp_mat_2)
          call sparse_embed_axpy(ham_jv, temp_mat_1, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
       end if
    end if

    if (present(only_SO_energy)) then
       perturbation2_min = perturbation1
       perturbation2_max = perturbation1
    else
       if (perturbation1.le.(3*mdl%nat)) then
          perturbation2_min = 1
          perturbation2_max = 3*mdl%nat
       else
          perturbation2_min = 3*mdl%nat + 1
          perturbation2_max = 3*mdl%nat + 3
       end if
    end if

    ! loop over epsilon perturbation
    do perturbation2 = perturbation2_min, perturbation2_max

       if (perturbation2.le.(3*mdl%nat)) then
          ! calculate FO_locpot_eps only if atomic displacement
          FO_locpot_eps = 0.0_DP

          ! calculate structure factor for all species
          struct_fac = cmplx(0.0_DP, 0.0_DP, kind=DP)
          call pseudo_make_structure_factor(struct_fac, mdl%elements, &
               mdl%fine_grid, mdl%par%num_pspecies, directions(perturbation2,:), &
               weights(perturbation2,:))

          if (pub_paw) then
             call paw_tcore_hartree_on_grid(FO_locpot_eps, dummy_gzero, &
                  struct_fac, struct_fac_classical, mdl%regions(1)%paw_sp, mdl%fine_grid, &
                  mdl%par)
          else
             call pseudopotentials_local_on_grid(FO_locpot_eps, dummy_gzero, &
                  struct_fac, struct_fac_classical, mdl%fine_grid, mdl%cell, &
                  mdl%regions(1)%pseudo_sp, mdl%par, mdl%par%nat_classical)
          end if

       end if
       ! end of FO_locpot_eps

       temp_energy = 0.0_DP

       call sparse_embed_scale(second_overlap,0.0_DP)

       if (pub_paw.or.pub_any_nl_proj) then
          call sparse_embed_scale(iG2pv_overlap, 0.0_DP)
          call sparse_embed_scale(iG1iG2pv_overlap, 0.0_DP)

          if ( (perturbation1.le.(3*mdl%nat)) .and. &
               (perturbation2.le.(3*mdl%nat))) then
             call calculate_second_overlap(mdl, val_rep, ngwf_basis, &
                  proj_basis, nl_projectors, second_overlap, iG1pv_overlap, &
                  iG2pv_overlap, iG1iG2pv_overlap, directions(perturbation1,:),&
                  weights(perturbation1,:), directions(perturbation2,:) , &
                  weights(perturbation2,:))
          end if
          call sparse_embed_transpose(viG2p_overlap, iG2pv_overlap)
       end if

       call sparse_embed_scale(first_overlap_vv, 0.0_DP)
       call sparse_embed_scale(first_overlap_jv,0.0_DP)

       if (pub_paw.or.pub_any_nl_proj) then
          call sparse_embed_scale(jiG2p_overlap, 0.0_DP)
          call sparse_embed_scale(dummy_iGpv_overlap, 0.0_DP)

          ! calculate first overlap w.r.t epsilon for rv
          if (perturbation2.le.(3*mdl%nat)) then
             call calculate_first_overlap(mdl, val_rep%ngwfs_on_grid, &
                  val_rep%sp_overlap, ngwf_basis, joint_rep%ngwfs_on_grid, &
                  joint_rep%sp_overlap, joint_basis, proj_basis, &
                  nl_projectors, first_overlap_jv, first_overlap_vv, &
                  dummy_iGpv_overlap, jiG2p_overlap, &
                  directions(perturbation2,:), weights(perturbation2,:))
          end if
       end if

       call sparse_embed_scale(A_op,0.0_DP)
       ! calculate A_op (vv rep) only if atomic displacement
       if ( (perturbation1.le.(3*mdl%nat)) .and. &
            (perturbation2.le.(3*mdl%nat))) then
          call calc_A_op(mdl, val_rep, ngwf_basis, directions(perturbation1,:),&
               weights(perturbation1,:), directions(perturbation2,:), &
               weights(perturbation2,:), dij, FO_dij, A_op, lhxc_fine, &
               FO_lhxc_fine, FO_locpot_eps, iG1pv_overlap, iG2pv_overlap, &
               iG1iG2pv_overlap, locpot_SO_mat)
       end if

       call sparse_embed_scale(B_op_jv, 0.0_DP)

       ! calculate B_op (v rep) w.r.t. epsilon only if atomic displacement
       if (perturbation2.le.(3*mdl%nat)) then
          call calc_B_op(mdl, joint_rep, joint_basis, val_rep, ngwf_basis, &
               directions(perturbation2,:), weights(perturbation2,:), dij, &
               B_op_jv, iG2pv_overlap, jiG2p_overlap, lhxc_fine, FO_locpot_eps)
       end if
       ! START BUILDING ENERGY TERMS

       if (.not.present(only_SO_energy)) then
          force_constant_matrix(perturbation1, perturbation2) = 0.0_DP
       else
          only_SO_energy = 0.0_DP
       end if

       call sparse_embed_create(temp_mat_1, second_overlap, val_rep%inv_overlap)
       call sparse_embed_create(temp_mat_2, temp_mat_1, val_ham%ham(1))

       call sparse_embed_product(temp_mat_1, second_overlap, val_rep%inv_overlap)
       call sparse_embed_product(temp_mat_2, temp_mat_1, val_ham%ham(1))
       call sparse_embed_scale(temp_mat_2, -1.0_DP)
       call sparse_embed_axpy(temp_mat_2, A_op, 1.0_DP)

       call sparse_embed_trace(temp_energy, temp_mat_2, val_denskern)

       call sparse_embed_destroy(temp_mat_1)
       call sparse_embed_destroy(temp_mat_2)

       if (.not.present(only_SO_energy)) then
          force_constant_matrix(perturbation1, perturbation2) = &
               force_constant_matrix(perturbation1, perturbation2) + &
               temp_energy
       else
          only_SO_energy = only_SO_energy + temp_energy
       end if

       call sparse_embed_create(temp_mat_1, first_overlap_jv, val_rep%inv_overlap)
       call sparse_embed_create(temp_mat_2, temp_mat_1, val_ham%ham(1))

       call sparse_embed_product(temp_mat_1, first_overlap_jv, val_rep%inv_overlap)
       call sparse_embed_product(temp_mat_2, temp_mat_1, val_ham%ham(1))
       call sparse_embed_scale(temp_mat_2, -1.0_DP)
       call sparse_embed_axpy(temp_mat_2, B_op_jv, 1.0_DP)

       call sparse_embed_trace(temp_energy, temp_mat_2, response_denskern_jv_CT)
       temp_energy = 2.0_DP * temp_energy

       call sparse_embed_destroy(temp_mat_1)
       call sparse_embed_destroy(temp_mat_2)

       if (.not.present(only_SO_energy)) then
          force_constant_matrix(perturbation1, perturbation2) = &
               force_constant_matrix(perturbation1, perturbation2) + &
               temp_energy
       else
          only_SO_energy = only_SO_energy + temp_energy
       end if

       ! only in the E-field case:
       if ( (perturbation1.gt.(3*mdl%nat)) .and. &
            (perturbation2.gt.(3*mdl%nat)).and. at_gamma) then

          call sparse_embed_trace(temp_energy, &
               dipole_mat_jv(mod(perturbation2,3*mdl%nat)), response_denskern_jv_CT)
          temp_energy = 2.0_DP * temp_energy

          if (.not.present(only_SO_energy)) then
             force_constant_matrix(perturbation1, perturbation2) = &
                  force_constant_matrix(perturbation1, perturbation2) + &
                  temp_energy
          else
             only_SO_energy = only_SO_energy + temp_energy
          end if
       end if
       ! end E-field case


       ! NLCC contribution present only for atomic displacements
       if ( (perturbation1.le.(3*mdl%nat)) .and. &
            (perturbation2.le.(3*mdl%nat)) ) then
          if (pub_paw) then

             if (pub_nlcc) then

                call paw_FO_nlcc_energy(fxc_buffer, mdl%fine_grid, mdl%cell, &
                     mdl%elements, mdl%regions(1)%paw_sp, mdl%par, temp_energy, &
                     directions(perturbation2,:), weights(perturbation2,:))

                if (.not.present(only_SO_energy)) then
                   force_constant_matrix(perturbation1, perturbation2) = &
                        force_constant_matrix(perturbation1, perturbation2) + &
                        temp_energy
                else
                   only_SO_energy = only_SO_energy + temp_energy
                end if

                call paw_FO_nlcc_energy(vxc_buffer,mdl%fine_grid, mdl%cell, &
                     mdl%elements, mdl%regions(1)%paw_sp, mdl%par, temp_energy, &
                     directions(perturbation2,:), weights(perturbation2,:), &
                     directions(perturbation1,:), weights(perturbation1,:))

                if (.not.present(only_SO_energy)) then
                   force_constant_matrix(perturbation1, perturbation2) = &
                        force_constant_matrix(perturbation1, perturbation2) + &
                        temp_energy
                else
                   only_SO_energy = only_SO_energy + temp_energy
                end if

             end if
          else
             if (pub_nlcc) then

                call pseudopotentials_FO_nlcc_energy(fxc_buffer, mdl%fine_grid,&
                     mdl%cell, mdl%elements, mdl%regions(1)%pseudo_sp, temp_energy, &
                     directions(perturbation2,:), weights(perturbation2,:))

                if (.not.present(only_SO_energy)) then
                   force_constant_matrix(perturbation1, perturbation2) = &
                        force_constant_matrix(perturbation1, perturbation2) + &
                        temp_energy
                else
                   only_SO_energy = only_SO_energy + temp_energy
                end if

                call pseudopotentials_FO_nlcc_energy(vxc_buffer, &
                     mdl%fine_grid, mdl%cell, mdl%elements, mdl%regions(1)%pseudo_sp, &
                     temp_energy, directions(perturbation2,:), &
                     weights(perturbation2,:), directions(perturbation1,:), &
                     weights(perturbation1,:))

                if (.not.present(only_SO_energy)) then
                   force_constant_matrix(perturbation1, perturbation2) = &
                        force_constant_matrix(perturbation1, perturbation2) + &
                        temp_energy
                else
                   only_SO_energy = only_SO_energy + temp_energy
                end if

             end if
          end if ! atomic displacements only

          ! first-order Lagrange multiplier terms
          if (.not.present(only_SO_energy)) then
             call sparse_embed_create(temp_mat_1, first_overlap_vv, val_denskern)
             call sparse_embed_create(temp_mat_2, val_denskern, temp_mat_1)
             call sparse_embed_create(temp_mat_3, FO_ham, temp_mat_2)

             call sparse_embed_product(temp_mat_1, first_overlap_vv, val_denskern)
             call sparse_embed_product(temp_mat_2, val_denskern, temp_mat_1)
             call sparse_embed_product(temp_mat_3, FO_ham, temp_mat_2)

             call sparse_embed_trace(temp_energy, temp_mat_3)

             call sparse_embed_destroy(temp_mat_1)
             call sparse_embed_destroy(temp_mat_2)
             call sparse_embed_destroy(temp_mat_3)

             if (.not.present(only_SO_energy)) then
                force_constant_matrix(perturbation1, perturbation2) = &
                     force_constant_matrix(perturbation1, perturbation2) - &
                     0.5_DP * temp_energy
             end if

             call sparse_embed_create(temp_mat_1, first_overlap_vv, &
                  response_denskern_jv_CT)
             call sparse_embed_create(temp_mat_2, val_denskern, temp_mat_1)
             call sparse_embed_create(temp_mat_3, ham_jv, temp_mat_2)

             call sparse_embed_product(temp_mat_1, first_overlap_vv, &
                  response_denskern_jv_CT)
             call sparse_embed_product(temp_mat_2, val_denskern, temp_mat_1)
             call sparse_embed_product(temp_mat_3, ham_jv, temp_mat_2)

             call sparse_embed_trace(temp_energy, temp_mat_3)

             call sparse_embed_destroy(temp_mat_1)
             call sparse_embed_destroy(temp_mat_2)
             call sparse_embed_destroy(temp_mat_3)

             if (.not.present(only_SO_energy)) then
                force_constant_matrix(perturbation1, perturbation2) = &
                     force_constant_matrix(perturbation1, perturbation2) - &
                     2.0_DP * 0.5_DP * temp_energy
             end if
          end if

          ! quadratic and FO_ham terms only if variational SO_energy required
          if (present(only_SO_energy)) then

             call sparse_embed_create(temp_mat_1,first_overlap_jv, &
                  val_rep%inv_overlap)
             call sparse_embed_create(temp_mat_2, temp_mat_1, val_ham%ham(1))

             call sparse_embed_product(temp_mat_1, first_overlap_jv, &
                  val_rep%inv_overlap)
             call sparse_embed_product(temp_mat_2, temp_mat_1, val_ham%ham(1))
             call sparse_embed_scale(temp_mat_2, -1.0_DP)
             call sparse_embed_axpy(temp_mat_2, FO_ham_jv, 1.0_DP)

             call sparse_embed_trace(temp_energy, temp_mat_2, response_denskern_jv_CT)
             temp_energy = 2.0_DP * temp_energy

             call sparse_embed_destroy(temp_mat_1)
             call sparse_embed_destroy(temp_mat_2)

             only_SO_energy = only_SO_energy + temp_energy

             ! only in the E-field case at Gamma:
             if ( (perturbation1.gt.(3*mdl%nat)) .and. &
                  (perturbation2.gt.(3*mdl%nat)) .and. at_gamma) then

                call sparse_embed_trace(temp_energy, &
                     dipole_mat_jv(mod(perturbation1,3*mdl%nat)), response_denskern_jv_CT)
                temp_energy = 2.0_DP * temp_energy

                only_SO_energy = only_SO_energy + temp_energy
             end if
             ! end E-field case

             call sparse_embed_create(temp_mat_1, response_denskern_jv, val_rep%overlap)
             call sparse_embed_product(temp_mat_1, response_denskern_jv, val_rep%overlap)
             call sparse_embed_create(temp_mat_2, temp_mat_1, response_denskern_jv_CT)
             call sparse_embed_product(temp_mat_2, temp_mat_1,response_denskern_jv_CT)

             call sparse_embed_trace(temp_energy, temp_mat_2, joint_ham%ham(1))
             call sparse_embed_destroy(temp_mat_1)
             call sparse_embed_destroy(temp_mat_2)

             only_SO_energy = only_SO_energy + temp_energy

             call sparse_embed_create(temp_mat_1, response_denskern_jv, val_ham%ham(1))
             call sparse_embed_product(temp_mat_1, response_denskern_jv, val_ham%ham(1))
             call sparse_embed_create(temp_mat_2, temp_mat_1, response_denskern_jv_CT)
             call sparse_embed_product(temp_mat_2, temp_mat_1,response_denskern_jv_CT)

             call sparse_embed_trace(temp_energy, temp_mat_2, joint_rep%overlap)
             call sparse_embed_destroy(temp_mat_1)
             call sparse_embed_destroy(temp_mat_2)

             only_SO_energy = only_SO_energy - temp_energy
          end if

          if (pub_on_root.and.(present(only_SO_energy))) then
             write(stdout,'(a,f20.10)') 'Second-order energy ', only_SO_energy
          end if

       end if
    end do

    if (pub_paw.or.pub_any_nl_proj) then
       call sparse_embed_destroy(iG2pv_overlap)
       call sparse_embed_destroy(viG2p_overlap)
       call sparse_embed_destroy(dummy_iGpv_overlap)
       call sparse_embed_destroy(iG1iG2pv_overlap)
       call sparse_embed_destroy(jiG2p_overlap)
       call sparse_embed_destroy(pj_overlap)
       call sparse_embed_destroy(pv_overlap)
    end if

    call sparse_embed_destroy(first_overlap_vv)
    call sparse_embed_destroy(first_overlap_jv)
    call sparse_embed_destroy(second_overlap)

    call sparse_embed_destroy(A_op)
    call sparse_embed_destroy(B_op_jv)

    if (.not.present(only_SO_energy)) then
       call sparse_embed_destroy(ham_jv)
    end if
    call sparse_embed_destroy(response_denskern_jv_CT)
    ! we shall not use the ion-ion interaction AT THE MOMENT


    deallocate(FO_locpot_eps, stat=ierr)
    call utils_dealloc_check('second_order_energy','FO_locpot_eps', ierr)
    deallocate(struct_fac,stat=ierr)
    call utils_alloc_check('second_order_energy','struct_fac',ierr)

    deallocate(struct_fac_classical,stat=ierr)
    call utils_dealloc_check('second_order_energy', 'struct_fac_classical', &
         ierr)

    ! scale back denskern and respkern to half the number of electrons if
    ! degenerate spins
    if(pub_num_spins==1) call sparse_embed_scale(val_denskern, 0.5_DP)
    if(pub_num_spins==1) call sparse_embed_scale(response_denskern_jv, 0.5_DP)

    ! pdh: sync procs
    call comms_barrier

    ! Stop timer
    call timer_clock('second_order_energy',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving second_order_energy'

  end subroutine second_order_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine phonons_xc_finite_diff(fxc_buffer,rho_perturbed,&
       ground_state_dens,grid,cell,nhat_den_grad,nhat_den_grad_gs)

    !==============================================================!
    ! ATTENTION: copied from linear_response, such that it is kept !
    ! safe from external modifications related to the TDDFT module !
    ! -------------------------------------------------------------!
    ! Subroutine calculates rho_perturbed*fxc in a finite diff     !
    ! approximation. This is to ensure that no higher order de-    !
    ! rivatives of the exchange-correlation functional are needed  !
    ! which can be troublesome to compute on the standard grid for !
    ! gradient-corrected functionals.                              !
    !==============================================================!
    ! Modified for embedding by Joseph Prentice, October 2018      !
    !==============================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_barrier
    use rundat, only: pub_lr_tddft_triplet, pub_num_spins,&
        pub_aug_den_dim
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    use xc, only: xc_fxc_finite_diff

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    real(kind=DP), intent(inout) :: fxc_buffer(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: rho_perturbed(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: ground_state_dens(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(in) :: nhat_den_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)
    real(kind=DP), optional, intent(inout) :: nhat_den_grad_gs(grid%ld1, &
         grid%ld2, grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)

    ! local variables
    real(kind=DP) :: epsilon_val
    integer :: ierr
    real(kind=DP), allocatable, dimension(:,:,:,:) :: temp_fxc
    real(kind=DP), allocatable, dimension(:,:,:,:) :: eff_dens

    ! Start timer
    call timer_clock('phonons_xc_finite_diff',1)

    ! finite difference parameter
    epsilon_val = 0.001_DP
    allocate(temp_fxc(grid%ld1,grid%ld2,grid%max_slabs12,2), stat=ierr)
    call utils_alloc_check('phonons_xc_finite_diff','temp_fxc',ierr)
    allocate(eff_dens(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins), &
         stat=ierr)
    call utils_alloc_check('phonons_xc_finite_diff','eff_dens',ierr)

    ! use a central difference scheme: rho0+eps*rho1
    eff_dens=ground_state_dens+epsilon_val*rho_perturbed

    temp_fxc = 0.0_DP

    ! if nhat_den_grad is supplied, the finite difference approximation
    ! has to be applied to it as well. It is stored temporarily in
    ! nhat_den_grad_gs
    if(present(nhat_den_grad) .and. present(nhat_den_grad_gs)) then
      nhat_den_grad_gs=nhat_den_grad_gs+epsilon_val*nhat_den_grad

      call xc_fxc_finite_diff(eff_dens, temp_fxc, grid, cell, &
           ground_state_dens,pub_aug_den_dim,nhat_den_grad_gs,nhat_den_grad_gs)
    else
       call xc_fxc_finite_diff(eff_dens,temp_fxc,grid,cell,ground_state_dens,0)
    endif

    fxc_buffer(:,:,:,1)=temp_fxc(:,:,:,1)

    temp_fxc = 0.0_DP

    ! second part of finite difference scheme: rho0-eps*rho1
    eff_dens=ground_state_dens-epsilon_val*rho_perturbed

    if(present(nhat_den_grad) .and. present(nhat_den_grad_gs)) then
      nhat_den_grad_gs=nhat_den_grad_gs-2.0_DP*epsilon_val*nhat_den_grad
      call xc_fxc_finite_diff(eff_dens,temp_fxc,grid,cell, &
           ground_state_dens,pub_aug_den_dim,nhat_den_grad_gs,nhat_den_grad_gs)
       ! restore nhat_dens_grad_gs
       nhat_den_grad_gs=nhat_den_grad_gs+epsilon_val*nhat_den_grad
    else
       call xc_fxc_finite_diff(eff_dens,temp_fxc,grid,cell,ground_state_dens,0)
    endif

    fxc_buffer(:,:,:,1)=fxc_buffer(:,:,:,1)-temp_fxc(:,:,:,1)

    ! scale by epsilon to compute the finite difference approximation
    fxc_buffer = 0.5_DP * fxc_buffer / epsilon_val

    deallocate(temp_fxc,stat=ierr)
    call utils_dealloc_check('phonons_xc_finite_diff','temp_fxc',ierr)
    deallocate(eff_dens,stat=ierr)
    call utils_dealloc_check('phonons_xc_finite_diff','eff_dens',ierr)

    call comms_barrier

    ! Stop timer
    call timer_clock('phonons_xc_finite_diff',2)


  end subroutine phonons_xc_finite_diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine FO_ngwf_gradient_lnv(contra_grad, cov_grad, proj_basis, &
       nl_projectors, lhxc_fine, FO_lhxc_fine, dij, FO_dij, mdl, val_rep, &
       val_basis, viGp_overlap, gradient_coeff_mat, gradient_coeff_mat_cov, &
       FOngwf_basis, FOngwf_rep, directions, weights, perturbation, at_gamma)

    !=========================================================================!
    ! This subroutine calculates the gradient of the LNV total energy         !
    ! function with respect to the expansion coefficients of the NGWFs.       !
    !=========================================================================!
    ! WARNING (gcc32): Many elements copied from ngwf_gradient_lnv. Efforts   !
    ! are being made to integrate it back into ngwf_gradient_lnv              !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,        !
    ! October 2018                                                            !
    !=========================================================================!

    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero, &
         data_fftbox_alloc, data_fftbox_dealloc
    use comms, only: comms_reduce, pub_on_root, comms_barrier, pub_my_proc_id
    use constants, only: DP, stdout, max_spins
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_batch_col_start
    use hf_exchange, only: hf_exchange_calculate
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: DKERN
    use model_type, only: MODEL
    use paw, only: paw_projector_overlap, paw_position_operator
    use potential, only: potential_input_to_workspace
    use projectors, only: PROJECTOR_SET
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_fftbox_batch_size, pub_any_nl_proj, pub_aug, pub_paw, &
         pub_debug_on_root, pub_num_spins, pub_num_kpoints, &
         PUB_1K
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_transpose, &
         sparse_embed_transpose_structure, sparse_embed_copy, sparse_embed_axpy, &
         sparse_embed_product, sparse_embed_scale, sparse_embed_array_create, &
         sparse_embed_array_destroy, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use sparse, only: SPAM3, sparse_create, sparse_copy, sparse_destroy, &
         sparse_get_element, sparse_put_element
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    type(FUNC_BASIS), intent(in) :: proj_basis(1)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(1)
    real(kind=DP), intent(in) :: lhxc_fine(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in) :: FO_lhxc_fine(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: dij(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: FO_dij(pub_num_spins)
    type(FUNCTIONS), intent(inout) :: contra_grad
    type(FUNCTIONS), intent(inout) :: cov_grad
    type(FUNC_BASIS), intent(in) :: val_basis(1)
    type(NGWF_REP), intent(in) :: val_rep
    type(SPAM3_EMBED), intent(in) :: viGp_overlap
    type(SPAM3_EMBED), intent(in) :: gradient_coeff_mat(7), gradient_coeff_mat_cov(7)
    type(NGWF_REP), intent(inout) :: FOngwf_rep
    type(FUNC_BASIS), intent(inout) :: FOngwf_basis(1)
    integer, intent(in) :: directions(mdl%nat)
    real(kind=DP), intent(in) :: weights(mdl%nat)
    integer, intent(in) :: perturbation
    logical, intent(in) :: at_gamma

    ! Local Variables
    real(kind=DP), allocatable :: extra_FO_lhxc_fine(:,:,:,:)
    type(SPAM3_EMBED), allocatable :: local_gradient_coeff_mat_cov(:,:)
    type(SPAM3_EMBED), allocatable :: local_gradient_coeff_mat(:,:)
    type(SPAM3_EMBED), allocatable :: local_gradient_FOval_coeff_mat_cov(:,:)
    type(SPAM3_EMBED), allocatable :: local_gradient_FOval_coeff_mat(:,:)
    type(SPAM3_EMBED), allocatable :: local_proj_coeff_mat(:,:)
    type(SPAM3_EMBED), allocatable :: local_proj_coeff_mat_cov(:,:)
    type(SPAM3_EMBED) :: pv_overlap, iGpv_overlap, pr_overlap
    type(FFTBOX_DATA), allocatable :: fftbox_batch(:,:,:), &
         fftbox_batch_proj(:,:,:), fftbox_batch_FOval(:,:,:)
    type(FFTBOX_DATA), allocatable :: fftbox_batch_cov(:,:,:), &
         fftbox_batch_proj_cov(:,:,:), fftbox_batch_FOval_cov(:,:,:)
    real(kind=DP), allocatable :: lhxc_dbl(:,:,:,:)
    real(kind=DP), allocatable :: FO_lhxc_dbl(:,:,:,:)
    integer :: nmat_FOval
    integer :: nmat
    integer :: nmat_proj
    integer :: is, iter
    integer :: ierr
    integer :: imat
    integer :: batch_size
    integer :: batch_count, ib
    integer :: n_batches
    integer :: local_start, local_end, local_len
    integer :: max_current_size ! maximum batch size over all procs
    integer, allocatable :: fa_box_start(:,:)
    integer, allocatable :: fa_start_in_box(:,:)
    real(kind=DP) :: R_atom(3), R_FFT(3), coc(3), r_el, o_el, fscale
    integer :: loc_iproj, iproj, iat, loc_iat, jproj, xyz

    logical :: loc_cmplx
    type(SPAM3_EMBED) :: temp_mat_1, temp_mat_2, Omat, Omat_pos(3)
    type(SPAM3) :: Omat_tmp(3)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering FO_ngwf_gradient_lnv'

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine FO_ngwf_gradient_lnv not ready yet for more&
         & than one k-point.')


    ! Start timer
    call timer_clock('FO_ngwf_gradient_lnv',1)

    ! ######## INITIALISATIONS #################################
    loc_cmplx = val_rep%ngwfs_on_grid(1)%iscmplx
    call data_set_to_zero(contra_grad)
    call data_set_to_zero(cov_grad)
    batch_size = pub_fftbox_batch_size

    nmat_FOval = 3
    nmat_proj = 2

    if (perturbation.le.(3*mdl%nat)) then
       nmat = 4
    else
       nmat = 5
    end if

    ! ndmh: projector-ngwf overlap matrix
    if (pub_any_nl_proj.or.pub_paw) then
       call sparse_embed_transpose_structure(pv_overlap%structure,val_rep%sp_overlap)
       call sparse_embed_create(pv_overlap)
       call sparse_embed_transpose(pv_overlap,val_rep%sp_overlap)

       call sparse_embed_transpose_structure(pr_overlap%structure, &
            FOngwf_rep%sp_overlap)
       call sparse_embed_create(pr_overlap)
       call sparse_embed_transpose(pr_overlap, FOngwf_rep%sp_overlap)

       call sparse_embed_create(iGpv_overlap, pv_overlap)
       call sparse_embed_transpose(iGpv_overlap, viGp_overlap)
    end if

    ! ndmh: allocate storage for coefficient matrices
    allocate(local_proj_coeff_mat(pub_num_spins,nmat_proj),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv','local_proj_coeff_mat',ierr)
    allocate(local_proj_coeff_mat_cov(pub_num_spins,nmat_proj),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv','local_proj_coeff_mat_cov', &
         ierr)
    allocate(local_gradient_coeff_mat(pub_num_spins,nmat),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv','local_gradient_coeff_mat', &
         ierr)
    allocate(local_gradient_coeff_mat_cov(pub_num_spins,nmat),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv', &
         'local_gradient_coeff_mat_cov', ierr)
    allocate(local_gradient_FOval_coeff_mat(pub_num_spins,nmat_FOval), &
         stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv', &
         'local_gradient_FOval_coeff_mat',ierr)
    allocate(local_gradient_FOval_coeff_mat_cov(pub_num_spins,nmat_FOval), &
         stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv', &
         'local_gradient_NEWval_coeff_mat_cov', ierr)

    ! explanation for components of coeff_mats and corresponding fftboxes
    ! local_gradient_coeff_mat(:,1) corresponds to phi(r)
    ! local_gradient_coeff_mat(:,2) corresponds to (T phi)(r)
    ! local_gradient_coeff_mat(:,3) corresponds to (V_lhxc phi)(r)
    ! local_gradient_coeff_mat(:,4) corresponds to (FO_lhxc)phi(r)
    ! local_gradient_coeff_mat(:,5) corresponds to (dipole_op)phi(r)
    ! local_gradient_FOval_coeff_mat(:,1) corresponds to FO_phi(r)
    ! local_gradient_FOval_coeff_mat(:,2) corresponds to (T FO_phi) (r)
    ! local_gradient_FOval_coeff_mat(:,3) corresponds to (V_lhxc FO_phi)(r)


    ! local_proj_coeff_mat(:,1) corresponds to projector(r)
    ! local_proj_coeff_mat(:,2) corresponds to FO_projector(r)
    ! the projector coeff mats correspond to either the NL part of the
    ! hamiltonian or the augmentation part of the overlap

    do is = 1,pub_num_spins
       call sparse_embed_create(local_gradient_coeff_mat(is,1), &
            gradient_coeff_mat(3))
       call sparse_embed_create(local_gradient_coeff_mat(is,2),&
            gradient_coeff_mat(4))
       call sparse_embed_create(local_gradient_coeff_mat(is,3),&
            gradient_coeff_mat(4))
       call sparse_embed_create(local_gradient_coeff_mat(is,4),&
            gradient_coeff_mat(6))
       if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then
          call sparse_embed_create(local_gradient_coeff_mat(is,5),&
               gradient_coeff_mat(7))
       end if
       call sparse_embed_create(local_gradient_FOval_coeff_mat(is,1), &
            gradient_coeff_mat(1))
       call sparse_embed_create(local_gradient_FOval_coeff_mat(is,2), &
            gradient_coeff_mat(2))
       call sparse_embed_create(local_gradient_FOval_coeff_mat(is,3), &
            gradient_coeff_mat(2))

       call sparse_embed_copy(local_gradient_coeff_mat(is,1), &
            gradient_coeff_mat(3))
       call sparse_embed_copy(local_gradient_coeff_mat(is,2),&
            gradient_coeff_mat(4))
       call sparse_embed_copy(local_gradient_coeff_mat(is,3),&
            gradient_coeff_mat(4))
       call sparse_embed_copy(local_gradient_coeff_mat(is,4),&
            gradient_coeff_mat(6))
       if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then
          call sparse_embed_copy(local_gradient_coeff_mat(is,5),&
               gradient_coeff_mat(7))
       end if
       call sparse_embed_copy(local_gradient_FOval_coeff_mat(is,1), &
            gradient_coeff_mat(1))
       call sparse_embed_copy(local_gradient_FOval_coeff_mat(is,2), &
            gradient_coeff_mat(2))
       call sparse_embed_copy(local_gradient_FOval_coeff_mat(is,3), &
            gradient_coeff_mat(2))

       call sparse_embed_create(local_gradient_coeff_mat_cov(is,1), &
            gradient_coeff_mat_cov(3))
       call sparse_embed_create(local_gradient_coeff_mat_cov(is,2),&
            gradient_coeff_mat_cov(4))
       call sparse_embed_create(local_gradient_coeff_mat_cov(is,3),&
            gradient_coeff_mat_cov(4))
       call sparse_embed_create(local_gradient_coeff_mat_cov(is,4),&
            gradient_coeff_mat_cov(6))
       if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then
          call sparse_embed_create(local_gradient_coeff_mat_cov(is,5),&
               gradient_coeff_mat_cov(7))
       end if
       call sparse_embed_create(local_gradient_FOval_coeff_mat_cov(is,1), &
            gradient_coeff_mat_cov(1))
       call sparse_embed_create(local_gradient_FOval_coeff_mat_cov(is,2), &
            gradient_coeff_mat_cov(2))
       call sparse_embed_create(local_gradient_FOval_coeff_mat_cov(is,3), &
            gradient_coeff_mat_cov(2))

       call sparse_embed_copy(local_gradient_coeff_mat_cov(is,1), &
            gradient_coeff_mat_cov(3))
       call sparse_embed_copy(local_gradient_coeff_mat_cov(is,2),&
            gradient_coeff_mat_cov(4))
       call sparse_embed_copy(local_gradient_coeff_mat_cov(is,3),&
            gradient_coeff_mat_cov(4))
       call sparse_embed_copy(local_gradient_coeff_mat_cov(is,4),&
            gradient_coeff_mat_cov(6))
       if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then
          call sparse_embed_copy(local_gradient_coeff_mat_cov(is,5),&
               gradient_coeff_mat_cov(7))
       end if
       call sparse_embed_copy(local_gradient_FOval_coeff_mat_cov(is,1), &
            gradient_coeff_mat_cov(1))
       call sparse_embed_copy(local_gradient_FOval_coeff_mat_cov(is,2), &
            gradient_coeff_mat_cov(2))
       call sparse_embed_copy(local_gradient_FOval_coeff_mat_cov(is,3), &
            gradient_coeff_mat_cov(2))


    end do

    fscale = 1.0_DP/sum(mdl%elements(:)%ion_charge)
    coc(1) =  sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%X)*fscale
    coc(2) =  sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%Y)*fscale
    coc(3) =  sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%Z)*fscale

    if (pub_paw) then
       ! Create matrix arrays and structures
       Omat%structure = 'E'
       call sparse_embed_create(Omat)
       ! Get projector overlap matrix
       call paw_projector_overlap(Omat%p,mdl%regions(1)%paw_sp)

       if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then ! if E-field only
          do xyz = 1,3
             call sparse_embed_create(Omat_pos(xyz),Omat)
          end do
          call sparse_embed_extract_from_array(Omat_tmp,Omat_pos)

          call paw_position_operator(Omat_tmp,mdl%regions(1)%paw_sp)

          call sparse_embed_destroy_extracted_array(Omat_tmp,Omat_pos,.true.)


          ! ndmh: cycle over projectors on this proc, applying correction to
          ! ndmh: r_sphere to move it to the atom centre
          do loc_iproj=1,proj_basis(1)%num_on_proc(pub_my_proc_id)
             iproj = loc_iproj + proj_basis(1)%first_on_proc(pub_my_proc_id) - 1
             iat = proj_basis(1)%atom_of_func(iproj)
             loc_iat = iat - mdl%par%first_atom_on_proc(pub_my_proc_id) + 1
             do jproj=proj_basis(1)%first_on_atom(iat), &
                  proj_basis(1)%first_on_atom(iat)+proj_basis(1)%num_on_atom(iat)-1

                ! Extract overlap element
                call sparse_get_element(o_el,Omat%p,jproj,iproj)

                R_atom(1) = mdl%par%elements_on_proc(loc_iat)%centre%x
                R_atom(2) = mdl%par%elements_on_proc(loc_iat)%centre%y
                R_atom(3) = mdl%par%elements_on_proc(loc_iat)%centre%z

                ! Extract element from r_sphere and shift by (R_atom-coc)*o_el
                ! ddor: get elements for one direction only if axis is specified
                do xyz = 1,3
                   call sparse_get_element(r_el,Omat_pos(xyz)%p,jproj,iproj)
                   r_el = (R_atom(xyz)-coc(xyz)) * o_el + r_el
                   call sparse_put_element(r_el,Omat_pos(xyz)%p,jproj,iproj)
                end do

             end do
          end do
       end if ! if E-field

    end if ! if paw

    do is = 1, pub_num_spins

       if(pub_paw) then
          call sparse_embed_create(temp_mat_1, pr_overlap, &
               local_gradient_FOval_coeff_mat(is,1))
          call sparse_embed_product(temp_mat_1, pr_overlap, &
               local_gradient_FOval_coeff_mat(is,1))
          call sparse_embed_create(local_proj_coeff_mat(is,1), Omat, temp_mat_1)
          call sparse_embed_product(local_proj_coeff_mat(is,1), Omat, temp_mat_1)
          call sparse_embed_destroy(temp_mat_1)
       end if

       if(pub_paw) then
          call sparse_embed_create(temp_mat_1, pr_overlap, &
               local_gradient_FOval_coeff_mat_cov(is,1))
          call sparse_embed_product(temp_mat_1, pr_overlap, &
               local_gradient_FOval_coeff_mat_cov(is,1))
          call sparse_embed_create(local_proj_coeff_mat_cov(is,1), Omat, temp_mat_1)
          call sparse_embed_product(local_proj_coeff_mat_cov(is,1), Omat, temp_mat_1)
          call sparse_embed_destroy(temp_mat_1)
       end if

       if(pub_any_nl_proj.or.pub_paw) then
          call sparse_embed_create(temp_mat_1, pr_overlap, &
               local_gradient_FOval_coeff_mat(is,2))
          call sparse_embed_product(temp_mat_1, pr_overlap, &
               local_gradient_FOval_coeff_mat(is,2))
          call sparse_embed_create(temp_mat_2, dij(is), temp_mat_1)
          call sparse_embed_product(temp_mat_2, dij(is), temp_mat_1)
          if (.not.pub_paw) then
             call sparse_embed_create(local_proj_coeff_mat(is,1), temp_mat_2)
             call sparse_embed_scale(local_proj_coeff_mat(is,1), 0.0_DP)
          end if
          call sparse_embed_axpy(local_proj_coeff_mat(is,1), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
       end if

       if(pub_any_nl_proj.or.pub_paw) then
          call sparse_embed_create(temp_mat_1, pr_overlap, &
               local_gradient_FOval_coeff_mat_cov(is,2))
          call sparse_embed_product(temp_mat_1, pr_overlap, &
               local_gradient_FOval_coeff_mat_cov(is,2))
          call sparse_embed_create(temp_mat_2, dij(is), temp_mat_1)
          call sparse_embed_product(temp_mat_2, dij(is), temp_mat_1)
          if (.not.pub_paw) then
             call sparse_embed_create(local_proj_coeff_mat_cov(is,1), temp_mat_2)
             call sparse_embed_scale(local_proj_coeff_mat_cov(is,1), 0.0_DP)
          end if
          call sparse_embed_axpy(local_proj_coeff_mat_cov(is,1), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
       end if


       if(pub_paw) then
          call sparse_embed_create(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat(is,1))
          call sparse_embed_product(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat(is,1))
          call sparse_embed_create(temp_mat_2, Omat, temp_mat_1)
          call sparse_embed_product(temp_mat_2, Omat, temp_mat_1)
          call sparse_embed_axpy(local_proj_coeff_mat(is,1), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
       end if

       if(pub_paw) then
          call sparse_embed_create(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat_cov(is,1))
          call sparse_embed_product(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat_cov(is,1))
          call sparse_embed_create(temp_mat_2, Omat, temp_mat_1)
          call sparse_embed_product(temp_mat_2, Omat, temp_mat_1)
          call sparse_embed_axpy(local_proj_coeff_mat_cov(is,1), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
       end if


       if(pub_any_nl_proj.or.pub_paw) then
          call sparse_embed_create(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat(is,2))
          call sparse_embed_product(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat(is,2))
          call sparse_embed_create(temp_mat_2, dij(is), temp_mat_1)
          call sparse_embed_product(temp_mat_2, dij(is), temp_mat_1)
          call sparse_embed_axpy(local_proj_coeff_mat(is,1), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
       end if


       if(pub_any_nl_proj.or.pub_paw) then
          call sparse_embed_create(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat_cov(is,2))
          call sparse_embed_product(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat_cov(is,2))
          call sparse_embed_create(temp_mat_2, dij(is), temp_mat_1)
          call sparse_embed_product(temp_mat_2, dij(is), temp_mat_1)
          call sparse_embed_axpy(local_proj_coeff_mat_cov(is,1), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
       end if


       if(pub_paw) then
          call sparse_embed_create(temp_mat_1, pv_overlap, gradient_coeff_mat(5))
          call sparse_embed_product(temp_mat_1, pv_overlap, gradient_coeff_mat(5))
          call sparse_embed_create(local_proj_coeff_mat(is,2), Omat, temp_mat_1)
          call sparse_embed_product(local_proj_coeff_mat(is,2), Omat, temp_mat_1)
          call sparse_embed_destroy(temp_mat_1)
       end if

       if(pub_paw) then
          call sparse_embed_create(temp_mat_1, pv_overlap, gradient_coeff_mat_cov(5))
          call sparse_embed_product(temp_mat_1, pv_overlap, gradient_coeff_mat_cov(5))
          call sparse_embed_create(local_proj_coeff_mat_cov(is,2), Omat, temp_mat_1)
          call sparse_embed_product(local_proj_coeff_mat_cov(is,2), Omat, temp_mat_1)
          call sparse_embed_destroy(temp_mat_1)
       end if

       if(pub_paw) then
          call sparse_embed_create(temp_mat_1, iGpv_overlap, gradient_coeff_mat(5))
          call sparse_embed_product(temp_mat_1, iGpv_overlap, gradient_coeff_mat(5))
          call sparse_embed_create(temp_mat_2, Omat, temp_mat_1)
          call sparse_embed_product(temp_mat_2, Omat, temp_mat_1)
          call sparse_embed_axpy(local_proj_coeff_mat(is,1), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
       end if

       if(pub_paw) then
          call sparse_embed_create(temp_mat_1, iGpv_overlap, &
               gradient_coeff_mat_cov(5))
          call sparse_embed_product(temp_mat_1, iGpv_overlap, &
               gradient_coeff_mat_cov(5))
          call sparse_embed_create(temp_mat_2, Omat, temp_mat_1)
          call sparse_embed_product(temp_mat_2, Omat, temp_mat_1)
          call sparse_embed_axpy(local_proj_coeff_mat_cov(is,1), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
       end if

       if(pub_any_nl_proj.or.pub_paw) then
          call sparse_embed_create(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat(is,4))
          call sparse_embed_product(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat(is,4))
          call sparse_embed_create(temp_mat_2, dij(is), temp_mat_1)
          call sparse_embed_product(temp_mat_2, dij(is), temp_mat_1)
          if (.not.pub_paw) then
             call sparse_embed_create(local_proj_coeff_mat(is,2), temp_mat_2)
             call sparse_embed_scale(local_proj_coeff_mat(is,2), 0.0_DP)
          end if
          call sparse_embed_axpy(local_proj_coeff_mat(is,2), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
       end if

       if(pub_any_nl_proj.or.pub_paw) then
          call sparse_embed_create(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat_cov(is,4))
          call sparse_embed_product(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat_cov(is,4))
          call sparse_embed_create(temp_mat_2, dij(is), temp_mat_1)
          call sparse_embed_product(temp_mat_2, dij(is), temp_mat_1)
          if (.not.pub_paw) then
             call sparse_embed_create(local_proj_coeff_mat_cov(is,2), temp_mat_2)
             call sparse_embed_scale(local_proj_coeff_mat_cov(is,2), 0.0_DP)
          end if
          call sparse_embed_axpy(local_proj_coeff_mat_cov(is,2), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
       end if


       if(pub_any_nl_proj.or.pub_paw) then
          call sparse_embed_create(temp_mat_1, iGpv_overlap, &
               local_gradient_coeff_mat(is,4))
          call sparse_embed_product(temp_mat_1, iGpv_overlap, &
               local_gradient_coeff_mat(is,4))
          call sparse_embed_create(temp_mat_2, dij(is), temp_mat_1)
          call sparse_embed_product(temp_mat_2, dij(is), temp_mat_1)
          call sparse_embed_axpy(local_proj_coeff_mat(is,1), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
       end if

       if(pub_any_nl_proj.or.pub_paw) then
          call sparse_embed_create(temp_mat_1, iGpv_overlap, &
               local_gradient_coeff_mat_cov(is,4))
          call sparse_embed_product(temp_mat_1, iGpv_overlap, &
               local_gradient_coeff_mat_cov(is,4))
          call sparse_embed_create(temp_mat_2, dij(is), temp_mat_1)
          call sparse_embed_product(temp_mat_2, dij(is), temp_mat_1)
          call sparse_embed_axpy(local_proj_coeff_mat_cov(is,1), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
       end if

       ! note that FO_dij exists only for PAW
       if(pub_paw) then
          call sparse_embed_create(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat(is,4))
          call sparse_embed_product(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat(is,4))
          call sparse_embed_create(temp_mat_2, FO_dij(is), temp_mat_1)
          call sparse_embed_product(temp_mat_2, FO_dij(is), temp_mat_1)
          call sparse_embed_axpy(local_proj_coeff_mat(is,1), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
       end if

       if(pub_paw) then
          call sparse_embed_create(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat_cov(is,4))
          call sparse_embed_product(temp_mat_1, pv_overlap, &
               local_gradient_coeff_mat_cov(is,4))
          call sparse_embed_create(temp_mat_2, FO_dij(is), temp_mat_1)
          call sparse_embed_product(temp_mat_2, FO_dij(is), temp_mat_1)
          call sparse_embed_axpy(local_proj_coeff_mat_cov(is,1), temp_mat_2, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
       end if


       if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then ! if E-field only
          if(pub_paw) then
             call sparse_embed_create(temp_mat_1, pv_overlap, &
                  local_gradient_coeff_mat(is,5))
             call sparse_embed_product(temp_mat_1, pv_overlap, &
                  local_gradient_coeff_mat(is,5))
             call sparse_embed_create(temp_mat_2, &
                  Omat_pos(mod(perturbation,3*mdl%nat)), temp_mat_1)
             call sparse_embed_product(temp_mat_2, &
                  Omat_pos(mod(perturbation,3*mdl%nat)), temp_mat_1)
             call sparse_embed_axpy(local_proj_coeff_mat(is,1), temp_mat_2, 1.0_DP)
             call sparse_embed_destroy(temp_mat_1)
             call sparse_embed_destroy(temp_mat_2)
          end if

          if(pub_paw) then
             call sparse_embed_create(temp_mat_1, pv_overlap, &
                  local_gradient_coeff_mat_cov(is,5))
             call sparse_embed_product(temp_mat_1, pv_overlap, &
                  local_gradient_coeff_mat_cov(is,5))
             call sparse_embed_create(temp_mat_2, &
                  Omat_pos(mod(perturbation,3*mdl%nat)), temp_mat_1)
             call sparse_embed_product(temp_mat_2, &
                  Omat_pos(mod(perturbation,3*mdl%nat)), temp_mat_1)
             call sparse_embed_axpy(local_proj_coeff_mat_cov(is,1), temp_mat_2,1.0_DP)
             call sparse_embed_destroy(temp_mat_1)
             call sparse_embed_destroy(temp_mat_2)
          end if
       end if

    end do ! is

    if (pub_paw) then
       call sparse_embed_destroy(Omat)
       if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then ! if E-field only
          do xyz = 1,3
             call sparse_embed_destroy(Omat_pos(xyz))
          end do
       end if
    end if

    ! ### END INITIALISATIONS #################################

    ! ndmh: in case the fine grid is denser than the double grid, create a
    ! ndmh: temporary double grid
    allocate(lhxc_dbl(mdl%dbl_grid%ld1,mdl%dbl_grid%ld2,&
         mdl%dbl_grid%max_group_slabs12,pub_num_spins), stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv','lhxc_dbl',ierr)

    allocate(FO_lhxc_dbl(mdl%dbl_grid%ld1,mdl%dbl_grid%ld2,&
         mdl%dbl_grid%max_group_slabs12,pub_num_spins), stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv','FO_lhxc_dbl',ierr)

    ! ndmh: fb start positions in box and start positions of FFT boxes
    allocate(fa_box_start(3, batch_size), stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv','fa_box_start',ierr)
    allocate(fa_start_in_box(3, batch_size), stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv','fa_start_in_box',ierr)

    ! ndmh: allocate storage for fftboxes for this batch
    allocate(fftbox_batch(pub_num_spins, nmat, batch_size), stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv','fftbox_batch',ierr)
    allocate(fftbox_batch_FOval(pub_num_spins, nmat_FOval,batch_size),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv','fftbox_batch_FOval',ierr)
    allocate(fftbox_batch_proj(pub_num_spins, nmat_proj, batch_size), stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv','fftbox_batch_proj',ierr)

    allocate(fftbox_batch_cov(pub_num_spins, nmat, batch_size), stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv','fftbox_batch_cov',ierr)
    allocate(fftbox_batch_FOval_cov(pub_num_spins, nmat_FOval,batch_size), &
         stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv','fftbox_batch_FOval_cov', &
         ierr)
    allocate(fftbox_batch_proj_cov(pub_num_spins, nmat_proj, batch_size), &
         stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_lnv','fftbox_batch_proj_cov',ierr)

    do batch_count = 1, batch_size
       do is = 1,pub_num_spins
          do imat = 1, nmat
             call data_fftbox_alloc(fftbox_batch(is, imat, batch_count), &
                  mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, &
                  mdl%fftbox%total_pt3, loc_cmplx)
          end do

          do imat = 1, nmat_FOval
             call data_fftbox_alloc(fftbox_batch_FOval(is, imat, batch_count),&
                  mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, &
                  mdl%fftbox%total_pt3, loc_cmplx)
          end do
       end do

       do imat = 1, nmat_proj
          do is = 1, pub_num_spins
             call data_fftbox_alloc(fftbox_batch_proj(is, imat, batch_count), &
                  mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, &
                  mdl%fftbox%total_pt3, loc_cmplx)
          end do
       end do
    end do

    do batch_count = 1, batch_size
       do is = 1,pub_num_spins
          do imat = 1, nmat
             call data_fftbox_alloc(fftbox_batch_cov(is, imat, batch_count), &
                  mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, &
                  mdl%fftbox%total_pt3, loc_cmplx)
          end do

          do imat = 1, nmat_FOval
             call data_fftbox_alloc(fftbox_batch_FOval_cov(is, imat, &
                  batch_count), mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, &
                  mdl%fftbox%total_pt3, loc_cmplx)
          end do
       end do

       do imat = 1, nmat_proj
          do is = 1, pub_num_spins
             call data_fftbox_alloc(fftbox_batch_proj_cov(is, imat, &
                  batch_count), mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, &
                  mdl%fftbox%total_pt3, loc_cmplx)
          end do
       end do
    end do


    ! ndmh: filter (or just copy) the lhxc potential to the double grid
    do is=1,pub_num_spins
       call potential_input_to_workspace(lhxc_dbl(:,:,:,is), &
            lhxc_fine(:,:,:,is),mdl%dbl_grid,mdl%fine_grid)
       call potential_input_to_workspace(FO_lhxc_dbl(:,:,:,is), &
            FO_lhxc_fine(:,:,:,is),mdl%dbl_grid,mdl%fine_grid)
    end do

    ! cks: number of row-steps per row-block
    n_batches = FOngwf_basis(1)%max_on_proc / batch_size
    if (mod(FOngwf_basis(1)%max_on_proc,batch_size) > 0) &
         n_batches = n_batches + 1

    ! ndmh: loop over batches of NGWFs
    local_start = 1
    do batch_count=1,n_batches

       if (pub_debug_on_root) write(stdout,'(a,2(i4,a))') 'DEBUG: Batch loop ',&
            batch_count, ' of ',n_batches,' in FO_ngwf_gradient_lnv'

       local_end = min(local_start + batch_size - 1, FOngwf_basis(1)%proc_num)
       local_len = local_end - local_start + 1

       ! cks: maximum size of current batch over all procs
       max_current_size = local_len
       call comms_reduce('MAX', max_current_size)

       call function_ops_batch_col_start(fa_box_start,fa_start_in_box, &
            batch_size,local_start,local_end,mdl%fftbox,mdl%cell, FOngwf_basis(1))

       ! cks: zero before accumulation
       do ib = 1, local_len
          do is = 1, pub_num_spins
             do imat = 1, nmat
                call data_set_to_zero(fftbox_batch(is,imat,ib))
             end do
             do imat = 1, nmat_FOval
                call data_set_to_zero(fftbox_batch_FOval(is,imat,ib))
             end do
             do imat = 1, nmat_proj
                call data_set_to_zero(fftbox_batch_proj(is,imat,ib))
             end do
          end do
       end do

       do ib = 1, local_len
          do is = 1, pub_num_spins
             do imat = 1, nmat
                call data_set_to_zero(fftbox_batch_cov(is,imat,ib))
             end do
             do imat = 1, nmat_FOval
                call data_set_to_zero(fftbox_batch_FOval_cov(is,imat,ib))
             end do
             do imat = 1, nmat_proj
                call data_set_to_zero(fftbox_batch_proj_cov(is,imat,ib))
             end do
          end do
       end do

       ! ndmh: deposit sums of various functions (the NGWFs, the Hamiltonian
       ! ndmh: acting on the NGWFS, nonlocal projectors, Hubbard projectors,
       ! ndmh: valence NGWFs in a conduction optimisation) to the accumulating
       ! ndmh: FFTboxes, then precondition the gradient, extract it from the
       ! ndmh: FFTBoxes to PPD storage, and then shave them according to the
       ! ndmh: NGWFs radii.

       call comms_barrier

       call FO_ngwf_gradient_batch(contra_grad, cov_grad, &  ! output
            nmat, nmat_FOval, nmat_proj, fftbox_batch, fftbox_batch_FOval, &
            fftbox_batch_proj, &
            fftbox_batch_cov, fftbox_batch_FOval_cov, fftbox_batch_proj_cov, &
            lhxc_dbl, FO_lhxc_dbl, mdl%dbl_grid, val_rep, &
            val_basis, proj_basis, nl_projectors, fa_box_start, &
            fa_start_in_box, batch_size, local_start, local_end, &
            local_proj_coeff_mat, local_proj_coeff_mat_cov, &
            local_gradient_coeff_mat, local_gradient_coeff_mat_cov, &
            local_gradient_FOval_coeff_mat, local_gradient_FOval_coeff_mat_cov,&
            mdl%fftbox, mdl%cell, pr_overlap, max_current_size, mdl, &
            FOngwf_basis, FOngwf_rep, directions, weights, coc, perturbation, &
            at_gamma)

       call comms_barrier

       local_start = local_start + batch_size
    end do

    call comms_barrier

    ! pdh: deallocate workspace
    do batch_count = batch_size, 1, -1
       do is = pub_num_spins, 1, -1
          do imat = nmat, 1, -1
             call data_fftbox_dealloc(fftbox_batch(is, imat, batch_count))
          end do
          do imat = nmat_FOval, 1, -1
             call data_fftbox_dealloc(fftbox_batch_FOval(is,imat,batch_count))
          end do
          do imat = nmat_proj, 1, -1
             call data_fftbox_dealloc(fftbox_batch_proj(is, imat, batch_count))
          end do
       end do
    end do

    do batch_count = batch_size, 1, -1
       do is = pub_num_spins, 1, -1
          do imat = nmat, 1, -1
             call data_fftbox_dealloc(fftbox_batch_cov(is, imat, batch_count))
          end do
          do imat = nmat_FOval, 1, -1
             call data_fftbox_dealloc(fftbox_batch_FOval_cov(is, imat, &
                  batch_count))
          end do
          do imat = nmat_proj, 1, -1
             call data_fftbox_dealloc(fftbox_batch_proj_cov(is, imat, &
                  batch_count))
          end do
       end do
    end do

    deallocate(fftbox_batch,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv','fftbox_batch',ierr)
    deallocate(fftbox_batch_FOval,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv','fftbox_batch_FOval',ierr)
    deallocate(fftbox_batch_proj,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv','fftbox_batch_proj',ierr)

    deallocate(fftbox_batch_cov,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv','fftbox_batch_cov',ierr)
    deallocate(fftbox_batch_FOval_cov,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv','fftbox_batch_FOval_cov', &
         ierr)
    deallocate(fftbox_batch_proj_cov,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv','fftbox_batch_proj_cov', &
         ierr)

    deallocate(fa_start_in_box,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv','fa_start_in_box',ierr)
    deallocate(fa_box_start,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv','fa_box_start',ierr)
    deallocate(lhxc_dbl, stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv','lhxc_dbl',ierr)
    deallocate(FO_lhxc_dbl, stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv','FO_lhxc_dbl',ierr)

    ! ndmh: deallocate coefficient matrices
    do is = 1,pub_num_spins
       do iter = 1,nmat
          call sparse_embed_destroy(local_gradient_coeff_mat(is,iter))
          call sparse_embed_destroy(local_gradient_coeff_mat_cov(is,iter))
       end do
       do iter = 1,nmat_FOval
          call sparse_embed_destroy(local_gradient_FOval_coeff_mat(is,iter))
          call sparse_embed_destroy(local_gradient_FOval_coeff_mat_cov(is,iter))
       end do
    end do

    deallocate(local_gradient_coeff_mat, stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv','local_gradient_coeff_mat',&
          ierr)
    deallocate(local_gradient_coeff_mat_cov, stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv', &
         'local_gradient_coeff_mat_cov', ierr)
    deallocate(local_gradient_FOval_coeff_mat, stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv', &
         'local_gradient_FOval_coeff_mat', ierr)
    deallocate(local_gradient_FOval_coeff_mat_cov, stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv', &
         'local_gradient_FOval_coeff_mat_cov', ierr)

    ! ndmh: deallocate coefficient matrices
    do is = 1,pub_num_spins
       do iter = 1,nmat_proj
          call sparse_embed_destroy(local_proj_coeff_mat(is,iter))
          call sparse_embed_destroy(local_proj_coeff_mat_cov(is,iter))
       end do
    end do

    deallocate(local_proj_coeff_mat, stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv','local_proj_coeff_mat',ierr)
    deallocate(local_proj_coeff_mat_cov, stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_lnv','local_proj_coeff_mat_cov',&
         ierr)


    if (pub_any_nl_proj.or.pub_aug) then
       call sparse_embed_destroy(pv_overlap)
       call sparse_embed_destroy(pr_overlap)
       call sparse_embed_destroy(iGpv_overlap)
    end if

    call timer_clock('FO_ngwf_gradient_lnv', 2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
       &FO_ngwf_gradient_lnv'

  end subroutine FO_ngwf_gradient_lnv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine FO_ngwf_gradient_batch(contra_grad, cov_grad, &  ! output
       nmat, nmat_FOval, nmat_proj, fftbox_batch, fftbox_batch_FOval, &
       fftbox_batch_proj, &
       fftbox_batch_cov, fftbox_batch_FOval_cov, fftbox_batch_proj_cov, &
       lhxc_dbl, FO_lhxc_dbl, dbl_grid, val_rep, val_basis, &
       proj_basis, nl_projectors, fa_box_start, fa_start_in_box, batch_size, &
       local_start, local_end, local_proj_coeff_mat, local_proj_coeff_mat_cov, &
       local_gradient_coeff_mat, local_gradient_coeff_mat_cov, &
       local_gradient_FOval_coeff_mat, local_gradient_FOval_coeff_mat_cov, &
       fftbox, cell, pr_overlap, max_current_size, mdl, FOngwf_basis, &
       FOngwf_rep, directions, weights, coc, perturbation, at_gamma)

    !=========================================================================!
    ! This subroutine returns the contravariant and covariant NGWF            !
    ! gradients for the NGWFs of the current batch. It does this              !
    ! by applying the Hamiltonian operator to the functions accumulated       !
    ! in the fftbox of each batch, then applying kinetic energy               !
    ! preconditioning if required and then by extracting the                  !
    ! relevant ppds from the fftboxes and shaving their values                !
    ! so that they are non-zero only within their spheres.                    !
    !-------------------------------------------------------------------------!
    ! WARNING (gcc32): Many elements copied from ngwf_gradient_batch. Efforts !
    ! are being made to integrate it back into ngwf_gradient_batch            !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,        !
    ! October 2018                                                            !
    !=========================================================================!

    use augmentation, only: augmentation_overlap
    use basis, only: basis_extract_function_from_box, basis_clean_function, &
         basis_location_func_wrt_cell, basis_ket_start_wrt_fftbox, &
         basis_point_wrt_box, basis_start_of_box_wrt_cell
    use cell_grid, only: GRID_INFO, cell_grid_extract_box
    use comms, only: comms_barrier, pub_my_proc_id
    use constants, only: DP, stdout, EDA_POLFRAGLOC_DEVEL
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero, &
         data_fftbox_alloc, data_fftbox_axpy, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_sum_fftbox_batch, &
        function_ops_brappd_ketppd
    use geometry, only: POINT, operator(*), operator(+)
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use paw, only: paw_species_calc_proj_prec_mat
    use projectors, only: PROJECTOR_SET, projectors_gradient_batch, &
         projectors_grad_precond_batch
    use rundat, only: pub_any_nl_proj, pub_aug, pub_paw, &
         pub_debug_on_root, pub_dbl_grid_scale, pub_num_spins, &
         pub_threads_fftbox, pub_spin_fac
!$  use rundat, only: pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_transpose, &
         sparse_embed_transpose_structure, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index, &
         sparse_create, sparse_copy, sparse_destroy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    integer, intent(in) :: batch_size
    type(NGWF_REP), intent(in) :: val_rep
    type(FUNC_BASIS), intent(in) :: val_basis(1)
    type(FUNC_BASIS), intent(in) :: proj_basis(1)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(1)
    integer, intent(in) :: fa_box_start(3,batch_size)
    integer, intent(in) :: fa_start_in_box(3,batch_size)
    integer, intent(in) :: local_start, local_end
    integer, intent(in) :: max_current_size
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FUNCTIONS), intent(inout) :: contra_grad
    type(FUNCTIONS), intent(inout) :: cov_grad
    integer, intent(in) :: nmat, nmat_FOval, nmat_proj
    type(FFTBOX_DATA), intent(inout) :: fftbox_batch(pub_num_spins, &
         nmat,batch_size)
    type(FFTBOX_DATA), intent(inout) :: fftbox_batch_FOval(pub_num_spins, &
         nmat_FOval, batch_size)
    type(FFTBOX_DATA), intent(inout) :: fftbox_batch_proj(pub_num_spins, &
         nmat_proj, batch_size)
    type(FFTBOX_DATA), intent(inout) :: fftbox_batch_cov(pub_num_spins, nmat, &
         batch_size)
   type(FFTBOX_DATA), intent(inout) :: fftbox_batch_FOval_cov(pub_num_spins, &
         nmat_FOval, batch_size)
    type(FFTBOX_DATA), intent(inout) :: fftbox_batch_proj_cov(pub_num_spins, &
         nmat_proj, batch_size)
    type(SPAM3_EMBED), intent(in) :: pr_overlap
    type(GRID_INFO), intent(in) :: dbl_grid
    real(kind=DP), intent(in) :: lhxc_dbl(dbl_grid%ld1, &
         dbl_grid%ld2, dbl_grid%max_group_slabs12, pub_num_spins)
    real(kind=DP), intent(in) :: FO_lhxc_dbl(dbl_grid%ld1, &
         dbl_grid%ld2, dbl_grid%max_group_slabs12, pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: local_gradient_coeff_mat(pub_num_spins,nmat)
    type(SPAM3_EMBED), intent(in) :: local_gradient_coeff_mat_cov(pub_num_spins,nmat)
    type(SPAM3_EMBED), intent(in) :: local_gradient_FOval_coeff_mat(pub_num_spins, &
         nmat_FOval)
    type(SPAM3_EMBED), intent(in) :: local_gradient_FOval_coeff_mat_cov( &
         pub_num_spins,nmat_FOval)

    type(SPAM3_EMBED), intent(in) :: local_proj_coeff_mat(pub_num_spins,nmat_proj)
    type(SPAM3_EMBED), intent(in) :: local_proj_coeff_mat_cov(pub_num_spins,nmat_proj)
    type(FUNC_BASIS), intent(in) :: FOngwf_basis(1)
    type(NGWF_REP), intent(in) :: FOngwf_rep
    integer, intent(in) :: directions(mdl%nat)
    real(kind=DP), intent(in) :: weights(mdl%nat)
    real(kind=DP), intent(in) :: coc(3)
    integer, intent(in) :: perturbation
    logical, intent(in) :: at_gamma

    ! Local Variables
    logical :: i_need_potential
    integer :: prev_start1, prev_start2, prev_start3
    integer :: fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl
    integer :: local_fa, global_fa
    integer :: batch_count, global_func
    integer :: is, i1, i2, i3, xyz, bs1, bs2, bs3
    type(POINT) :: R_FFT,r1,r2,r3
    type(POINT) :: a1,a2,a3, a1dbl, a2dbl, a3dbl
    integer :: ierr, iter
    integer :: idx_len
    real(kind=DP) :: common_fac, cbg1, cbg2, cbg3
    integer, allocatable, dimension(:) :: overlap_idx
    real(kind=DP), allocatable, dimension(:,:,:) :: lhxc_fftbox_dbl
    real(kind=DP), allocatable, dimension(:,:,:) :: FO_lhxc_fftbox_dbl
    real(kind=DP), allocatable, dimension(:,:,:) :: kin_buffer_1
    real(kind=DP), allocatable, dimension(:,:,:) :: kin_buffer_2
    real(kind=DP), allocatable, dimension(:,:,:) :: buffer_dbl
    real(kind=DP), allocatable, dimension(:,:,:) :: FO_buffer_dbl
    real(kind=DP), allocatable, dimension(:,:,:,:) :: dwork_box_dbl
    complex(kind=DP), allocatable, dimension(:,:,:) :: coarse_work,fine_work
    real(kind=DP), allocatable  :: r_op(:,:,:,:)

    ! jcap: temp coefficient matrices
    type(SPAM3), allocatable, dimension(:,:) :: coeff_mat, coeff_mat_cov, &
         FOval_coeff_mat, FOval_coeff_mat_cov, proj_coeff_mat, &
         proj_coeff_mat_cov
    integer :: imat

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering FO_ngwf_gradient_batch'

    call timer_clock('FO_ngwf_gradient_batch',1)

    !jmecmplx
    call utils_assert(.not. any(fftbox_batch%iscmplx), 'Error in&
         & FO_ngwf_gradient_batch: not ready yet for complex NGWFs.')

    call gcc32_ngwf_grad_init_precond_recip(mdl%fftbox)

    if (pub_paw) then
       call paw_species_calc_proj_prec_mat(nl_projectors(1),mdl%cell,mdl%fftbox, &
            precond_func_recip,mdl%regions(1)%paw_sp)
    end if

    common_fac = cell%weight

    ! jcap: allocate temporary matrices
    allocate(coeff_mat(pub_num_spins,nmat),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','coeff_mat',ierr)
    allocate(coeff_mat_cov(pub_num_spins,nmat),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','coeff_mat_cov',ierr)
    allocate(FOval_coeff_mat(pub_num_spins,nmat_FOval),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','FOval_coeff_mat',ierr)
    allocate(FOval_coeff_mat_cov(pub_num_spins,nmat_FOval),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','FOval_coeff_mat_cov',ierr)
    allocate(proj_coeff_mat(pub_num_spins,nmat),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','proj_coeff_mat',ierr)
    allocate(proj_coeff_mat_cov(pub_num_spins,nmat),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','proj_coeff_mat_cov',ierr)

    ! jcap: copy coeff_mat into temporaries
    ! rc2013: extract the relevant coefficient matrix for this region
    do is=1,pub_num_spins
       do imat=1,nmat
          call sparse_create(coeff_mat(is,imat),&
               local_gradient_coeff_mat(is,imat)%p)
          call sparse_copy(coeff_mat(is,imat),&
               local_gradient_coeff_mat(is,imat)%p)
          call sparse_create(coeff_mat_cov(is,imat),&
               local_gradient_coeff_mat_cov(is,imat)%p)
          call sparse_copy(coeff_mat_cov(is,imat),&
               local_gradient_coeff_mat_cov(is,imat)%p)
       end do
       do imat=1,nmat_FOval
          call sparse_create(FOval_coeff_mat(is,imat),&
               local_gradient_FOval_coeff_mat(is,imat)%p)
          call sparse_copy(FOval_coeff_mat(is,imat),&
               local_gradient_coeff_mat(is,imat)%p)
          call sparse_create(FOval_coeff_mat_cov(is,imat),&
               local_gradient_FOval_coeff_mat_cov(is,imat)%p)
          call sparse_copy(FOval_coeff_mat_cov(is,imat),&
               local_gradient_FOval_coeff_mat_cov(is,imat)%p)
       end do
       do imat=1,nmat_proj
          call sparse_create(proj_coeff_mat(is,imat),&
               local_proj_coeff_mat(is,imat)%p)
          call sparse_copy(proj_coeff_mat(is,imat),&
               local_proj_coeff_mat(is,imat)%p)
          call sparse_create(proj_coeff_mat_cov(is,imat),&
               local_proj_coeff_mat_cov(is,imat)%p)
          call sparse_copy(proj_coeff_mat_cov(is,imat),&
               local_proj_coeff_mat_cov(is,imat)%p)
       end do
    end do

    ! pdh: overlap matrix index
    idx_len = sparse_index_length(val_rep%overlap%p)
    allocate(overlap_idx(idx_len),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','overlap_idx',ierr)
    call sparse_generate_index(overlap_idx,val_rep%overlap%p)

    ! for val_basis functions
    call function_ops_sum_fftbox_batch(fftbox_batch, fftbox, cell, &
         val_rep%ngwfs_on_grid(1), val_basis(1), fa_box_start, batch_size, &
         local_start, local_end, overlap_idx, idx_len, &
         coeff_mat, 1, nmat, common_fac)
    call function_ops_sum_fftbox_batch(fftbox_batch_cov, fftbox, cell, &
         val_rep%ngwfs_on_grid(1), val_basis(1), fa_box_start, batch_size, &
         local_start, local_end, overlap_idx, idx_len, &
         coeff_mat_cov, 1, nmat, common_fac)

    ! for FOval_basis functions
    call function_ops_sum_fftbox_batch(fftbox_batch_FOval, fftbox, cell, &
         FOngwf_rep%ngwfs_on_grid(1), FOngwf_basis(1), fa_box_start, &
         batch_size, local_start, local_end, overlap_idx, idx_len, &
         FOval_coeff_mat, 1, nmat_FOval, common_fac)
    call function_ops_sum_fftbox_batch(fftbox_batch_FOval_cov, fftbox, cell, &
         FOngwf_rep%ngwfs_on_grid(1), FOngwf_basis(1), fa_box_start, &
         batch_size, local_start, local_end, overlap_idx, idx_len, &
         FOval_coeff_mat_cov, 1, nmat_FOval, common_fac)

    call comms_barrier

    deallocate(overlap_idx,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_batch','overlap_idx',ierr)

    call comms_barrier

    ! ndmh: allocate generic complex FFTbox workspace
    allocate(coarse_work(fftbox%total_ld1,fftbox%total_ld2,&
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','coarse_work',ierr)


    ! KKKKKKKKKKK KINETIC KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

    allocate(kin_buffer_1(fftbox%total_ld1,fftbox%total_ld2,&
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','kin_buffer_1',ierr)
    allocate(kin_buffer_2(fftbox%total_ld1,fftbox%total_ld2, &
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','kin_buffer_2',ierr)

    do is=1,pub_num_spins
       do local_fa=local_start,local_start+max_current_size+1
          batch_count = local_fa - local_start + 1

          if (local_fa <= local_end) then
             ! Calculate T on sum of phi_beta
             kin_buffer_1 = 0.0_DP
             kin_buffer_2 = 0.0_DP
             coarse_work = cmplx(0.0_DP,0.0_DP,kind=DP)

             call FO_ngwf_gradient_kinetic(kin_buffer_1, kin_buffer_2, &
                  fftbox_batch_FOval(is,2,batch_count)%d,  &
                  fftbox_batch(is,2,batch_count)%d, &
                  coarse_work,fftbox)

             fftbox_batch_FOval(is,2,batch_count)%d(:,:,:) = kin_buffer_1(:,:,:)
             fftbox_batch(is,2,batch_count)%d(:,:,:) = kin_buffer_2(:,:,:)

             kin_buffer_1 = 0.0_DP
             kin_buffer_2 = 0.0_DP
             coarse_work = cmplx(0.0_DP,0.0_DP,kind=DP)

             call FO_ngwf_gradient_kinetic(kin_buffer_1, kin_buffer_2, &
                  fftbox_batch_FOval_cov(is,2,batch_count)%d,  &
                  fftbox_batch_cov(is,2,batch_count)%d, &
                  coarse_work,fftbox)

             fftbox_batch_FOval_cov(is,2,batch_count)%d(:,:,:) = kin_buffer_1(:,:,:)
             fftbox_batch_cov(is,2,batch_count)%d(:,:,:) = kin_buffer_2(:,:,:)


          end if

       end do
    end do

    deallocate(kin_buffer_1,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_batch','kin_buffer_1',ierr)
    deallocate(kin_buffer_2,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_batch','kin_buffer_2',ierr)

    ! KKKKKKKK END KINETIC KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

    ! LLLLLL LHXC POTENTIAL LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    allocate(r_op(3,fftbox%total_ld1_dbl,&
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','r_op',ierr)
    allocate(lhxc_fftbox_dbl(fftbox%total_ld1_dbl,&
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','lhxc_fftbox_dbl',ierr)
    allocate(FO_lhxc_fftbox_dbl(fftbox%total_ld1_dbl,&
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','FO_lhxc_fftbox_dbl',ierr)
    allocate(buffer_dbl(fftbox%total_ld1_dbl,fftbox%total_ld2_dbl,&
         dbl_grid%max_group_slabs12),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','buffer_dbl',ierr)
    allocate(FO_buffer_dbl(fftbox%total_ld1_dbl,fftbox%total_ld2_dbl,&
         dbl_grid%max_group_slabs12),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','FO_buffer_dbl',ierr)
    allocate(dwork_box_dbl(fftbox%total_ld1_dbl,fftbox%total_ld2_dbl,&
         fftbox%total_pt3_dbl,2),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','dwork_box_dbl',ierr)
    allocate(fine_work(fftbox%total_ld1_dbl,fftbox%total_ld2_dbl,&
         fftbox%total_pt3_dbl),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','fine_work',ierr)

    ! cks: initialisation
    buffer_dbl = 0.0_DP
    FO_buffer_dbl = 0.0_DP

    a1dbl = (1.0_DP / dbl_grid%n1) * cell%a1
    a2dbl = (1.0_DP / dbl_grid%n2) * cell%a2
    a3dbl = (1.0_DP / dbl_grid%n3) * cell%a3

    a1 = (1.0_DP / cell%total_pt1) * cell%a1
    a2 = (1.0_DP / cell%total_pt2) * cell%a2
    a3 = (1.0_DP / cell%total_pt3) * cell%a3

    ! pdh: loop over spins
    do is=1,pub_num_spins

       prev_start1= -1111111
       prev_start2= -2222222
       prev_start3= -3333333

       do local_fa=local_start,local_start+max_current_size+1
          batch_count = local_fa - local_start + 1

          if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then ! for E-field

             global_fa = val_basis(1)%first_on_proc(pub_my_proc_id) + local_fa - 1

             ! Centre of function wrt fftbox in terms of grid spacings
             call basis_point_wrt_box(cbg1, cbg2, cbg3, &  ! output
                  val_basis(1)%spheres(local_fa)%centre, fftbox%total_pt1, &
                  fftbox%total_pt2, fftbox%total_pt3, cell, fftbox)

             ! Start of fftbox wrt cell in terms of grid-point number
             call basis_start_of_box_wrt_cell(bs1,bs2,bs3, &
                  val_basis(1)%spheres(local_fa)%centre,cbg1,cbg2,cbg3,cell)

             ! Find vector to origin of FFTbox
             R_FFT = real(bs1-1,kind=DP) * a1 + real(bs2-1,kind=DP) * a2 &
                  + real(bs3-1,kind=DP) * a3
          end if


          ! ndmh: find FFTbox start position on double grid
          if (local_fa <= local_end) then
             if (pub_dbl_grid_scale > 1.0_DP) then
                fftbox_start1_dbl = 2*fa_box_start(1,batch_count) - 1
                fftbox_start2_dbl = 2*fa_box_start(2,batch_count) - 1
                fftbox_start3_dbl = 2*fa_box_start(3,batch_count) - 1
             else ! ndmh: dbl_grid is scale 1.0
                fftbox_start1_dbl = fa_box_start(1,batch_count)
                fftbox_start2_dbl = fa_box_start(2,batch_count)
                fftbox_start3_dbl = fa_box_start(3,batch_count)
             end if
          else
             fftbox_start1_dbl = -1234
             fftbox_start2_dbl = -1234
             fftbox_start3_dbl = -1234
          end if

          ! Extract potential from cell to FFT box if position of FFT box
          ! is different from last time.
          i_need_potential = .false.
          if ((fftbox_start1_dbl /= prev_start1 .or. &
               fftbox_start2_dbl /= prev_start2 .or. &
               fftbox_start3_dbl /= prev_start3) .and. local_fa <= local_end) &
               i_need_potential = .true.

          call cell_grid_extract_box(lhxc_fftbox_dbl,&
               buffer_dbl, lhxc_dbl(:,:,:,is), dbl_grid, &
               fftbox%total_pt1_dbl, fftbox%total_pt2_dbl, &
               fftbox%total_pt3_dbl, fftbox%total_ld1_dbl, &
               fftbox%total_ld2_dbl, fftbox_start1_dbl, &
               fftbox_start2_dbl, fftbox_start3_dbl, i_need_potential, .true.)

          call cell_grid_extract_box(FO_lhxc_fftbox_dbl,&
               FO_buffer_dbl, FO_lhxc_dbl(:,:,:,is), dbl_grid, &
               fftbox%total_pt1_dbl, fftbox%total_pt2_dbl, &
               fftbox%total_pt3_dbl, fftbox%total_ld1_dbl, &
               fftbox%total_ld2_dbl, fftbox_start1_dbl, &
               fftbox_start2_dbl, fftbox_start3_dbl, i_need_potential, .true.)

          r_op = 0.0_DP

          if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then
             do i3=1,fftbox%total_pt3_dbl
                r3 = real(i3-1,kind=DP) * a3dbl
                do i2=1,fftbox%total_ld2_dbl
                   r2 = r3 + real(i2-1,kind=DP) * a2dbl
                   do i1=1,fftbox%total_ld1_dbl
                      r1 = r2 + real(i1-1,kind=DP) * a1dbl
                      r_op(1,i1,i2,i3) = r1%X + R_FFT%X - coc(1)
                      r_op(2,i1,i2,i3) = r1%Y + R_FFT%Y - coc(2)
                      r_op(3,i1,i2,i3) = r1%Z + R_FFT%Z - coc(3)
                   enddo
                enddo
             enddo
          end if

          prev_start1 = fftbox_start1_dbl
          prev_start2 = fftbox_start2_dbl
          prev_start3 = fftbox_start3_dbl

          if (local_fa <= local_end) then

             ! Calculate V_loc on sum of phi_beta
             if (pub_dbl_grid_scale>1.0_DP) then
                ! ndmh: interpolate function sums to fine grid, multiply
                ! ndmh: by fine grid potential, then filter back to coarse grid

                coarse_work = cmplx(0.0_DP,0.0_DP,kind=DP)
                fine_work = cmplx(0.0_DP, 0.0_DP, kind=DP)

                call FO_ngwf_gradient_local( &
                     fftbox_batch(is,3,batch_count)%d, &
                     fftbox_batch_cov(is,3,batch_count)%d,&
                     lhxc_fftbox_dbl, lhxc_fftbox_dbl, &
                     dwork_box_dbl,coarse_work,fine_work,fftbox)

                coarse_work = cmplx(0.0_DP,0.0_DP,kind=DP)
                fine_work = cmplx(0.0_DP, 0.0_DP, kind=DP)

                call FO_ngwf_gradient_local( &
                     fftbox_batch(is,4,batch_count)%d, &
                     fftbox_batch_cov(is,4,batch_count)%d,&
                     FO_lhxc_fftbox_dbl, FO_lhxc_fftbox_dbl, &
                     dwork_box_dbl,coarse_work,fine_work,fftbox)

                coarse_work = cmplx(0.0_DP,0.0_DP,kind=DP)
                fine_work = cmplx(0.0_DP, 0.0_DP, kind=DP)

                call FO_ngwf_gradient_local( &
                     fftbox_batch_FOval(is,3,batch_count)%d, &
                     fftbox_batch_FOval_cov(is,3,batch_count)%d,&
                     lhxc_fftbox_dbl, lhxc_fftbox_dbl, &
                     dwork_box_dbl,coarse_work,fine_work,fftbox)

                if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then ! E-field
                   coarse_work = cmplx(0.0_DP,0.0_DP,kind=DP)
                   fine_work = cmplx(0.0_DP, 0.0_DP, kind=DP)

                   call FO_ngwf_gradient_local( &
                        fftbox_batch(is,5,batch_count)%d, &
                        fftbox_batch_cov(is,5,batch_count)%d,&
                        r_op(mod(perturbation,3*mdl%nat),:,:,:), &
                        r_op(mod(perturbation,3*mdl%nat),:,:,:), &
                        dwork_box_dbl,coarse_work,fine_work,fftbox)
                end if

             else
                ! ndmh: dbl_grid_scale:1, so no need to interpolate
                fftbox_batch(is,3,batch_count)%d(:,:,:) = &
                     fftbox_batch(is,3,batch_count)%d * &
                     lhxc_fftbox_dbl(:,:,:)
                fftbox_batch(is,4,batch_count)%d(:,:,:) = &
                     fftbox_batch(is,4,batch_count)%d * &
                     FO_lhxc_fftbox_dbl(:,:,:)
                fftbox_batch_FOval(is,3,batch_count)%d(:,:,:) = &
                     fftbox_batch_FOval(is,3,batch_count)%d * &
                     lhxc_fftbox_dbl(:,:,:)
                if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then ! E-field
                   fftbox_batch(is,5,batch_count)%d(:,:,:) = &
                        fftbox_batch(is,5,batch_count)%d * &
                        r_op(mod(perturbation,3*mdl%nat),:,:,:)
                end if

                fftbox_batch_cov(is,3,batch_count)%d(:,:,:) = &
                     fftbox_batch_cov(is,3,batch_count)%d * &
                     lhxc_fftbox_dbl(:,:,:)
                fftbox_batch_cov(is,4,batch_count)%d(:,:,:) = &
                     fftbox_batch_cov(is,4,batch_count)%d * &
                     FO_lhxc_fftbox_dbl(:,:,:)
                fftbox_batch_FOval_cov(is,3,batch_count)%d(:,:,:) = &
                     fftbox_batch_FOval_cov(is,3,batch_count)%d * &
                     lhxc_fftbox_dbl(:,:,:)
                if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then ! E-field
                   fftbox_batch_cov(is,5,batch_count)%d(:,:,:) = &
                        fftbox_batch_cov(is,5,batch_count)%d * &
                        r_op(mod(perturbation,3*mdl%nat),:,:,:)
                end if
             end if

          end if

       end do
    end do

    deallocate(fine_work,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_batch','fine_work',ierr)
    deallocate(dwork_box_dbl,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_batch','dwork_box_dbl',ierr)
    deallocate(buffer_dbl,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_batch','buffer_dbl',ierr)
    deallocate(FO_buffer_dbl,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_batch','FO_buffer_dbl',ierr)
    deallocate(lhxc_fftbox_dbl,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_batch','lhxc_fftbox_dbl',ierr)
    deallocate(FO_lhxc_fftbox_dbl,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_batch','FO_lhxc_fftbox_dbl',ierr)
    deallocate(coarse_work,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_batch','coarse_work',ierr)
    deallocate(r_op,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_batch','r_op',ierr)
    ! LLLLLL END LHXC POTENTIAL LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    ! DDDDDD DIPOLE MATRIX PART DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
    ! gcc32 WARNING: NOT YET IMPLEMENTED, BUT IMPACT SHOULD BE NEGLIGIBLE !!!
    ! DDDDDD END DIPOLE MATRIX PART DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD


    ! NNNNNN NON-LOCAL POTENTIAL NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    if (pub_any_nl_proj.or.pub_aug) then
       call projectors_gradient_batch(fftbox_batch_proj, 1,nmat_proj, &
            local_start, local_end, batch_size, fa_box_start, &
            FOngwf_basis(1), proj_basis(1), cell, fftbox, proj_coeff_mat, &
            pr_overlap%p, nl_projectors(1), directions=directions, weights=weights)
       call projectors_gradient_batch(fftbox_batch_proj_cov, 1,nmat_proj, &
            local_start, local_end, batch_size, fa_box_start, &
            FOngwf_basis(1), proj_basis(1), cell, fftbox, &
            proj_coeff_mat_cov, pr_overlap%p, nl_projectors(1), &
            directions=directions, weights=weights)
    end if
    ! NNNNNN END NON-LOCAL POTENTIAL NNNNNNNNNNNNNNNNNNNNNNNNNNNN

    call comms_barrier

    allocate(coarse_work(mdl%fftbox%total_ld1,mdl%fftbox%total_ld2,&
         mdl%fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('FO_ngwf_gradient_batch','coarse_work',ierr)

    ! ndmh: allocate generic complex FFTbox workspace

    do local_fa=local_start,local_end
       batch_count = local_fa - local_start + 1

       ! gcc32: add all fftbox contributions in fftbox_batch(is,1,batch_count)
       do is = 1,pub_num_spins
          fftbox_batch(is,1,batch_count)%d(:,:,:) = &
               fftbox_batch(is,1,batch_count)%d +&
               fftbox_batch(is,2,batch_count)%d + &
               fftbox_batch(is,3,batch_count)%d + &
               fftbox_batch(is,4,batch_count)%d + &
               fftbox_batch_proj(is,1,batch_count)%d + &
               fftbox_batch_proj(is,2,batch_count)%d + &
               fftbox_batch_FOval(is,1,batch_count)%d + &
               fftbox_batch_FOval(is,2,batch_count)%d + &
               fftbox_batch_FOval(is,3,batch_count)%d

          if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then ! E-field
             fftbox_batch(is,1,batch_count)%d(:,:,:) = &
                  fftbox_batch(is,1,batch_count)%d +&
                  fftbox_batch(is,5,batch_count)%d
          end if

          fftbox_batch_cov(is,1,batch_count)%d(:,:,:) = &
               fftbox_batch_cov(is,1,batch_count)%d + &
               fftbox_batch_cov(is,2,batch_count)%d + &
               fftbox_batch_cov(is,3,batch_count)%d + &
               fftbox_batch_cov(is,4,batch_count)%d + &
               fftbox_batch_proj_cov(is,1,batch_count)%d + &
               fftbox_batch_proj_cov(is,2,batch_count)%d + &
               fftbox_batch_FOval_cov(is,1,batch_count)%d + &
               fftbox_batch_FOval_cov(is,2,batch_count)%d + &
               fftbox_batch_FOval_cov(is,3,batch_count)%d

          if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then ! E-field
             fftbox_batch_cov(is,1,batch_count)%d(:,:,:) = &
                  fftbox_batch_cov(is,1,batch_count)%d +&
                  fftbox_batch_cov(is,5,batch_count)%d
          end if

       end do

    end do

    call comms_barrier

    do local_fa=local_start,local_end
       batch_count = local_fa - local_start + 1

       ! ndmh: Average gradients over spins for spin-polarised calculations.
       if (pub_num_spins == 2) then
          fftbox_batch(1,1,batch_count)%d(:,:,:) = fftbox_batch(1,1,batch_count)%d + &
               fftbox_batch(2,1,batch_count)%d
          fftbox_batch_cov(1,1,batch_count)%d(:,:,:) = &
               fftbox_batch_cov(1,1,batch_count)%d + &
               fftbox_batch_cov(2,1,batch_count)%d
       end if

    end do

    call comms_barrier

!  gcc32 TEST
   ! PPPPPPP PAW PRECONDITIONER PROJ-OVERLAP TERM PPPPPPPPPPPP
  if (pub_paw) then

     ! Calculate overlap of this batch of NGWF covariant gradients
     ! with preconditioned projectors
     ! Multiply <g_a|q_i> matrix by O^prec_ij to get coefficients
     ! of preconditioned projectors to add back into gradient
     call projectors_grad_precond_batch(fftbox_batch_cov, 1, &
          local_start, local_end, batch_size, fa_box_start, &
          mdl%cell, mdl%fftbox, FOngwf_basis(1), pr_overlap%p, &
          nl_projectors(1), precond_func_recip)
  end if
  ! PPPPPPP END PAW PRECONDITIONER PROJ-OVERLAP TERM PPPPPPPP

  call comms_barrier

  ! rc2013: deallocate temporary matrices
  do is=1,pub_num_spins
     do imat=1,nmat
        call sparse_destroy(coeff_mat(is,imat))
        call sparse_destroy(coeff_mat_cov(is,imat))
     end do
     do imat=1,nmat_FOval
        call sparse_destroy(FOval_coeff_mat(is,imat))
        call sparse_destroy(FOval_coeff_mat_cov(is,imat))
     end do
     do imat=1,nmat_proj
        call sparse_destroy(proj_coeff_mat(is,imat))
        call sparse_destroy(proj_coeff_mat_cov(is,imat))
     end do
  end do

  deallocate(coeff_mat,stat=ierr)
  call utils_dealloc_check('FO_ngwf_gradient_batch','coeff_mat',ierr)
  deallocate(coeff_mat_cov,stat=ierr)
  call utils_dealloc_check('FO_ngwf_gradient_batch','coeff_mat_cov',ierr)
  deallocate(FOval_coeff_mat,stat=ierr)
  call utils_dealloc_check('FO_ngwf_gradient_batch','FOval_coeff_mat',ierr)
  deallocate(FOval_coeff_mat_cov,stat=ierr)
  call utils_dealloc_check('FO_ngwf_gradient_batch','FOval_coeff_mat_cov',ierr)
  deallocate(proj_coeff_mat,stat=ierr)
  call utils_dealloc_check('FO_ngwf_gradient_batch','proj_coeff_mat',ierr)
  deallocate(proj_coeff_mat_cov,stat=ierr)
  call utils_dealloc_check('FO_ngwf_gradient_batch','proj_coeff_mat_cov',ierr)

  do local_fa = local_start,local_end, 2
     batch_count = local_fa - local_start + 1

     ! ======= KIN. ENERGY PRECONDITIONING

     ! Copy the covariant gradient into complex workspace array
     if (local_fa < local_end) then
        coarse_work = &
             cmplx(fftbox_batch_cov(1,1,batch_count)%d, &
             fftbox_batch_cov(1,1,batch_count+1)%d,kind=DP)
     else
        coarse_work = cmplx(fftbox_batch_cov(1,1,batch_count)%d, &
             0.0_DP,kind=DP)
     end if

     ! Forward FFT the covariant gradient to reciprocal space
     call fourier_apply_box('Coarse','Forward',coarse_work, &
          omp=pub_threads_fftbox)

     ! Apply kinetic energy preconditioning to covariant gradient
     coarse_work = coarse_work * precond_func_recip

     ! Backward FFT the covariant gradient to real space
     call fourier_apply_box('Coarse','Backward',coarse_work, &
          omp=pub_threads_fftbox)

     ! Copy the preconditioned covariant gradient out of the complex
     ! workspace array
     fftbox_batch_cov(1,1,batch_count)%d(:,:,:) = real(coarse_work,kind=DP)

     if (local_fa < local_end) then
        fftbox_batch_cov(1,1,batch_count+1)%d(:,:,:) = aimag(coarse_work)
     end if

     ! ======= END OF KIN. ENERGY PRECONDITIONING
  end do

     call comms_barrier

    do local_fa=local_start,local_end
       batch_count = local_fa - local_start + 1

       ! cks: shaving - stage 1 / extract ppds from the FFT box
       ! extract the contravariant gradient from ham_box to PPDs
       call basis_extract_function_from_box(contra_grad, &
            fftbox_batch(1,1,batch_count), FOngwf_basis(1)%spheres(local_fa),&
            FOngwf_basis(1)%tight_boxes(local_fa), fa_start_in_box(1,batch_count),&
            fa_start_in_box(2,batch_count), fa_start_in_box(3,batch_count),&
            FOngwf_basis(1)%spheres(local_fa)%offset,cell,fftbox)

       ! extract the covariant gradient from ham_box to PPDs
       call basis_extract_function_from_box(cov_grad, &
            fftbox_batch_cov(1,1,batch_count), FOngwf_basis(1)%spheres(local_fa), &
            FOngwf_basis(1)%tight_boxes(local_fa), fa_start_in_box(1,batch_count),&
            fa_start_in_box(2,batch_count),fa_start_in_box(3,batch_count), &
            FOngwf_basis(1)%spheres(local_fa)%offset, cell, fftbox)

       ! cks: shaving - stage 2 / zero points outside NGWF sphere in PPD rep.
       ! smmd: if smoothing scheme applied previously, no needs to clean the
       ! smmd: covariant gradient (nor the contravariant gradient which is
       ! smmd: always used within contra.DOT.cov expressions)
       if (cell%n_pts > 1) then
          call basis_clean_function(contra_grad, &
               FOngwf_basis(1)%spheres(local_fa), cell,fftbox)
          call basis_clean_function(cov_grad, &
               FOngwf_basis(1)%spheres(local_fa), cell,fftbox)
       endif


    end do

    ! ndmh: END APPLY SPIN-AVERAGING, THEN EXTRACT AND SHAVE GRADIENT

    deallocate(coarse_work,stat=ierr)
    call utils_dealloc_check('FO_ngwf_gradient_batch','coarse_work',ierr)

    call comms_barrier

    if (allocated(precond_func_recip)) then
       deallocate(precond_func_recip,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_exit','precond_func_recip',ierr)
    end if

    call comms_barrier

    call timer_clock('FO_ngwf_gradient_batch',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving FO_ngwf_gradient_batch'

  end subroutine FO_ngwf_gradient_batch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine FO_ngwf_gradient_kinetic(out1,out2,in1,in2,zwork_box,fftbox)

    !=========================================================================!
    ! This subroutine applies the kinetic energy operator to an NGWF in a     !
    ! FFTbox.                                                                 !
    !-------------------------------------------------------------------------!
    ! WARNING (gcc32): Many elements copied from ngwf_gradient_kinetic.       !
    ! Efforts are being made to integrate it back into ngwf_gradient_kinetic  !
    !=========================================================================!

    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box_pair
    use kinetic, only: kinetic_apply_on_box
    use rundat, only: pub_threads_fftbox

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(out)   :: out1(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)
    real(kind=DP), intent(out)   :: out2(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)
    real(kind=DP), intent(inout) :: in1(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)
    real(kind=DP), intent(inout) :: in2(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)
    complex(kind=DP), intent(out) :: zwork_box(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)

    ! Forward Fourier transform to reciprocal space
    call fourier_apply_box_pair('Coarse','Forward',in1,in2,zwork_box, &
         omp=pub_threads_fftbox)

    ! Apply kinetic operator in reciprocal space
    call kinetic_apply_on_box(zwork_box,fftbox)

    ! Backward Fourier transform to real space
    call fourier_apply_box_pair('Coarse','Backward',out1,out2,zwork_box, &
         omp=pub_threads_fftbox)

    return

  end subroutine FO_ngwf_gradient_kinetic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine FO_ngwf_gradient_local(data_box1, data_box2, pot_dbl1, pot_dbl2, &
       dwork_box_dbl,coarse_work,fine_work,fftbox)

    !=========================================================================!
    ! This subroutine applies the local potential to a pair of NGWFs in       !
    ! FFTboxes.                                                               !
    !-------------------------------------------------------------------------!
    ! WARNING (gcc32): Many elements copied from ngwf_gradient_local. Efforts !
    ! are being made to integrate it back into ngwf_gradient_local            !
    !=========================================================================!

    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_interpolate, fourier_filter
    use rundat, only: pub_threads_per_fftbox
    use utils, only: utils_use_var

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(inout) :: data_box1(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)
    real(kind=DP), intent(inout) :: data_box2(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)
    real(kind=DP), intent(in)    :: pot_dbl1(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl)
    real(kind=DP), intent(in)    :: pot_dbl2(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl)
    real(kind=DP), intent(out)   :: dwork_box_dbl(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl,2)
    complex(kind=DP), intent(out) :: fine_work(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl)
    complex(kind=DP), intent(out) :: coarse_work(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)

    call fourier_interpolate(coarse_work,fine_work,data_box1,data_box2, &
         dwork_box_dbl(:,:,:,1),dwork_box_dbl(:,:,:,2))

!!$OMP PARALLEL WORKSHARE NUM_THREADS(pub_threads_per_fftbox)
    dwork_box_dbl(:,:,:,1) = dwork_box_dbl(:,:,:,1) * pot_dbl1
    dwork_box_dbl(:,:,:,2) = dwork_box_dbl(:,:,:,2) * pot_dbl2
!!$OMP END PARALLEL WORKSHARE

    call fourier_filter(coarse_work,fine_work,dwork_box_dbl(:,:,:,1), &
         dwork_box_dbl(:,:,:,2),data_box1, data_box2)

    ! jd: Silence var unused warning when in non-OMP compiles
    call utils_use_var(pub_threads_per_fftbox)

  end subroutine FO_ngwf_gradient_local


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_phonons_regenerate_SO_energy(mdl, val_basis, val_rep, val_ham, &
        joint_basis, joint_rep, joint_ham, auxresponse_kernel_jv, &
        directions, weights, perturbation, proj_basis, nl_projectors, dij, &
        denskern, dens_fine, lhxc_fine, vxc_buffer, FO_lpseudo_fine, &
        FO_tcore_density, SO_locpot_mat, at_gamma, energy)

    !=========================================================================!
    ! This subroutine regenerates the variational second-order energy for a   !
    ! set of new first-order NGWFs. This is used in the CG optimisation of the!
    ! first-order NGWFs                                                       !
    !-------------------------------------------------------------------------!
    ! One also regenerates the response kernel for a set of new first-order   !
    ! NGWFs, by applying the gauge constraint on the auxiliary response kernel!
    !-------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in 2016-2017                          !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,        !
    ! October 2018                                                            !
    !=========================================================================!

    use augmentation, only: aug_projector_denskern, augmentation_overlap
    use comms, only: pub_on_root, pub_my_proc_id, comms_barrier
    use function_basis, only: FUNC_BASIS
    use constants, only: stdout
    use density, only:  density_on_grid
    use function_ops, only: function_ops_brappd_ketppd
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use paw, only: paw_tcore_hartree_on_grid
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_num_spins, pub_debug_on_root,pub_aug,pub_paw, &
         pub_any_nl_proj
    use sparse, only: SPAM3
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_axpy, sparse_embed_product, &
         sparse_embed_transpose_structure, sparse_embed_transpose, &
         sparse_embed_scale, sparse_embed_copy, sparse_embed_conjugate, &
         sparse_embed_array_scale, sparse_embed_extract_from_array, &
         sparse_embed_destroy_extracted_array
    use utils, only: utils_abort, utils_alloc_check,utils_dealloc_check, &
         utils_unit
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    type(FUNC_BASIS), intent(in) :: val_basis(1)
    type(NGWF_REP), intent(in) :: val_rep
    type(NGWF_HAM), intent(in) :: val_ham
    type(FUNC_BASIS), intent(in) :: joint_basis(1)
    type(NGWF_REP), intent(in) :: joint_rep
    type(NGWF_HAM), intent(in) :: joint_ham
    type(SPAM3_EMBED), intent(in) :: auxresponse_kernel_jv(pub_num_spins)
    integer, intent(in) :: directions(3*mdl%nat+3, mdl%nat)
    real(kind=DP), intent(in) :: weights(3*mdl%nat+3, mdl%nat)
    integer, intent(in) :: perturbation
    type(FUNC_BASIS), intent(in) :: proj_basis(1)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(1)
    type(SPAM3_EMBED), intent(in) :: dij(pub_num_spins)
    type(SPAM3_EMBED_ARRAY), intent(inout) :: denskern
    real(kind=DP), intent(inout) :: dens_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in) :: vxc_buffer(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in) :: FO_lpseudo_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12)
    real(kind=DP), intent(in) :: FO_tcore_density(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12)
    type(SPAM3_EMBED), intent(in) :: SO_locpot_mat
    logical, intent(in) :: at_gamma
    real(kind=DP), intent(inout) :: energy

    ! Local Variables
    integer :: is, i1, i2, i3, ierr, iter, val_num_ngwfs

    type(SPAM3_EMBED) :: FO_ham, FO_ham_jv
    type(SPAM3_EMBED) :: FO_dij(pub_num_spins)
    type(SPAM3_EMBED) :: response_kernel_jv(pub_num_spins), &
         response_kernel_jv_CT(pub_num_spins)
    type(SPAM3_EMBED) :: FO_overlap_jv, FO_overlap
    type(SPAM3_EMBED) :: rho_ij(pub_num_spins), rho_ij_dproj(pub_num_spins)

    real(kind=DP), allocatable :: FO_lhxc_fine(:,:,:,:)
    real(kind=DP), allocatable :: fxc_buffer(:,:,:,:)

    real(kind=DP) :: dummy_gzero

    type(SPAM3_EMBED) :: temp_mat_1, temp_mat_2, temp_mat_3, temp_mat_4, temp_mat_5
    type(SPAM3_EMBED) :: viGp_overlap, iGpv_overlap
    type(SPAM3_EMBED) :: FO_ham_jv_fixed, FO_ham_jv_fixed_temp
    type(SPAM3_EMBED) :: dijps, pv_overlap
    type(SPAM3_EMBED) :: jiGp_overlap
    type(SPAM3_EMBED) :: condkern_new, dipole_mat_jv(3)

    integer :: iter1, iter2, iter1_eff, iter2_eff, num_ngwfs
    real(kind=DP) :: element

    ! jcap: temporary arrays for input to aug_projector_denskern
    type(SPAM3), dimension(pub_num_spins) :: kern_array, aug_array

    if (pub_debug_on_root) write(stdout, '(a)') &
         'DEBUG: Entering lr_phonons_regenerate_SO_energy routine.'

    ! Start timer
    call timer_clock('lr_phonons_regenerate_SO_energy',1)

    num_ngwfs = val_basis(1)%num

    energy = 0.0_DP

    allocate(FO_lhxc_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
             mdl%fine_grid%max_slabs12, pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_phonons_regenerate_SO_energy','FO_lhxc_fine',&
         ierr)
    allocate(fxc_buffer(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
             mdl%fine_grid%max_slabs12, pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_phonons_regenerate_SO_energy','fxc_buffer',&
         ierr)

    val_num_ngwfs = val_basis(1)%num

    if (pub_paw) then
       do is=1,pub_num_spins
          call sparse_embed_create(FO_dij(is),dij(is))
       end do
    end if

    if (pub_paw) then
       do is=1,pub_num_spins
          call sparse_embed_create(rho_ij(is),FO_dij(is))
          call sparse_embed_create(rho_ij_dproj(is),FO_dij(is))
       end do
    end if

    if (pub_paw.or.pub_any_nl_proj) then
       call sparse_embed_create(viGp_overlap, val_rep%sp_overlap)

       call sparse_embed_create(jiGp_overlap, joint_rep%sp_overlap)

       call sparse_embed_transpose_structure(pv_overlap%structure, val_rep%sp_overlap)
       call sparse_embed_create(pv_overlap)
       call sparse_embed_transpose(pv_overlap, val_rep%sp_overlap)

       call sparse_embed_create(iGpv_overlap, pv_overlap)
    end if

    call sparse_embed_create(FO_ham, val_rep%overlap)
    call sparse_embed_transpose_structure(FO_ham_jv%structure, &
         joint_rep%cross_overlap)
    call sparse_embed_create(FO_ham_jv)

    call sparse_embed_create(FO_ham_jv_fixed, FO_ham_jv)
    call sparse_embed_create(FO_ham_jv_fixed_temp, FO_ham_jv)

    call sparse_embed_create(FO_overlap, val_rep%overlap)
    call sparse_embed_create(FO_overlap_jv, FO_ham_jv)

    do iter1 = 1, 3
       call sparse_embed_create(dipole_mat_jv(iter1), FO_overlap_jv)
    end do

    if ((perturbation.gt.(3*mdl%nat)).and.at_gamma) then
       ! calculate dipole matrix in jv rep
       call sparse_embed_transpose_structure(temp_mat_1%structure, &
            joint_rep%cross_overlap)
       call sparse_embed_create(temp_mat_1)
       call sparse_embed_transpose(temp_mat_1, joint_rep%cross_overlap)

       call calculate_dipole_mat(dipole_mat_jv, mdl, joint_rep, &
            joint_basis, val_rep, val_basis, temp_mat_1, proj_basis, &
            nl_projectors, .false.)

       call sparse_embed_destroy(temp_mat_1)
    end if



    if (pub_paw) then
       if (pub_num_spins == 1) call sparse_embed_array_scale(denskern,2.0_DP)
       call sparse_embed_extract_from_array(aug_array,rho_ij)
       call sparse_embed_extract_from_array(kern_array,denskern%m(:,1))
       call aug_projector_denskern(aug_array, kern_array, val_rep%sp_overlap%p)
       call sparse_embed_destroy_extracted_array(aug_array,rho_ij,.true.)
       call sparse_embed_destroy_extracted_array(kern_array)
       if (pub_num_spins == 1) call sparse_embed_array_scale(denskern,0.5_DP)
    end if

    !===========================================================
    if (pub_paw.or.pub_any_nl_proj) then
       call sparse_embed_scale(iGpv_overlap,0.0_DP)
       call sparse_embed_scale(jiGp_overlap, 0.0_DP)
    end if

    call sparse_embed_scale(FO_overlap_jv,0.0_DP)
    call sparse_embed_scale(FO_overlap,0.0_DP)

    if (pub_paw.or.pub_any_nl_proj) then
       if (perturbation.le.(3*mdl%nat)) then
          call calculate_first_overlap(mdl, val_rep%ngwfs_on_grid, &
               val_rep%sp_overlap, val_basis, joint_rep%ngwfs_on_grid, &
               joint_rep%sp_overlap, joint_basis, proj_basis, &
               nl_projectors, FO_overlap_jv, FO_overlap, iGpv_overlap, &
               jiGp_overlap, directions(perturbation,:), &
               weights(perturbation,:), rho_ij_dproj, denskern%m(:,1))
       end if

       call sparse_embed_transpose(viGp_overlap, iGpv_overlap)
    end if
    !===========================================================

    !! ============== BUILD FIXED PART OF FO HAMILTONIAN ======================
    call sparse_embed_scale(FO_ham_jv_fixed, 0.0_DP)

    if (pub_paw.or.pub_any_nl_proj) then
       call sparse_embed_create(dijps, dij(1), iGpv_overlap)
       call sparse_embed_product(dijps, dij(1), iGpv_overlap)
       call sparse_embed_product(FO_ham_jv_fixed, joint_rep%sp_overlap, &
            dijps)

       call sparse_embed_scale(dijps, 0.0_DP)
       call sparse_embed_product(dijps, dij(1), pv_overlap)
       call sparse_embed_product(FO_ham_jv_fixed_temp, jiGp_overlap, dijps)

       call sparse_embed_axpy(FO_ham_jv_fixed, FO_ham_jv_fixed_temp, 1.0_DP)
    end if

    !! ============== END BUILD FIXED PART OF FO HAMILTONIAN =================
    call comms_barrier

    ! now the actual conduction projector #####################################
    call sparse_embed_create(condkern_new, joint_rep%inv_overlap)
    call sparse_embed_copy(condkern_new, joint_rep%inv_overlap)

    call sparse_embed_transpose_structure(temp_mat_1%structure, &
         joint_rep%cross_overlap)
    call sparse_embed_create(temp_mat_1)
    call sparse_embed_transpose(temp_mat_1, joint_rep%cross_overlap)


    call sparse_embed_create(temp_mat_2, joint_rep%inv_overlap, temp_mat_1)
    call sparse_embed_product(temp_mat_2, joint_rep%inv_overlap, temp_mat_1)
    call sparse_embed_destroy(temp_mat_1)

    call sparse_embed_create(temp_mat_1, temp_mat_2, denskern%m(1,1))
    call sparse_embed_product(temp_mat_1, temp_mat_2, denskern%m(1,1))
    call sparse_embed_destroy(temp_mat_2)
    call sparse_embed_create(temp_mat_2, temp_mat_1, joint_rep%cross_overlap)
    call sparse_embed_product(temp_mat_2, temp_mat_1, joint_rep%cross_overlap)
    call sparse_embed_destroy(temp_mat_1)
    call sparse_embed_create(temp_mat_1, temp_mat_2, joint_rep%inv_overlap)
    call sparse_embed_product(temp_mat_1, temp_mat_2, joint_rep%inv_overlap)
    call sparse_embed_axpy(condkern_new, temp_mat_1, -1.0_DP)
    call sparse_embed_destroy(temp_mat_1)
    call sparse_embed_destroy(temp_mat_2)
    ! end of condkern creation ###############################################


    ! ========= restrict response kernel to new gauge ====================
    do is = 1, pub_num_spins
       call sparse_embed_create(response_kernel_jv(is), auxresponse_kernel_jv(is))
       call sparse_embed_transpose_structure(response_kernel_jv_CT(is)%structure, &
            response_kernel_jv(is))
       call sparse_embed_create(response_kernel_jv_CT(is))

       call sparse_embed_create(temp_mat_1, condkern_new, joint_rep%overlap)
       call sparse_embed_product(temp_mat_1, condkern_new, joint_rep%overlap)
       call sparse_embed_create(temp_mat_2, temp_mat_1, auxresponse_kernel_jv(is))
       call sparse_embed_product(temp_mat_2, temp_mat_1, auxresponse_kernel_jv(is))
       call sparse_embed_destroy(temp_mat_1)
       call sparse_embed_create(temp_mat_1, temp_mat_2, val_rep%overlap)
       call sparse_embed_product(temp_mat_1, temp_mat_2, val_rep%overlap)
       call sparse_embed_destroy(temp_mat_2)
       call sparse_embed_create(temp_mat_2, temp_mat_1, denskern%m(1,1))
       call sparse_embed_product(temp_mat_2, temp_mat_1, denskern%m(1,1))
       call sparse_embed_destroy(temp_mat_1)

       call sparse_embed_transpose_structure(temp_mat_3%structure, &
            joint_rep%cross_overlap)
       call sparse_embed_create(temp_mat_3)
       call sparse_embed_transpose(temp_mat_3, joint_rep%cross_overlap)
       call sparse_embed_create(temp_mat_1, joint_rep%inv_overlap, temp_mat_3)
       call sparse_embed_product(temp_mat_1, joint_rep%inv_overlap, temp_mat_3)
       call sparse_embed_destroy(temp_mat_3)
       call sparse_embed_create(temp_mat_3, temp_mat_1, denskern%m(1,1))
       call sparse_embed_product(temp_mat_3, temp_mat_1, denskern%m(1,1))
       call sparse_embed_destroy(temp_mat_1)
       call sparse_embed_create(temp_mat_1, temp_mat_3, FO_overlap)
       call sparse_embed_product(temp_mat_1, temp_mat_3, FO_overlap)
       call sparse_embed_destroy(temp_mat_3)
       call sparse_embed_create(temp_mat_3, temp_mat_1, denskern%m(1,1))
       call sparse_embed_product(temp_mat_3, temp_mat_1, denskern%m(1,1))
       call sparse_embed_destroy(temp_mat_1)

       call sparse_embed_axpy(temp_mat_2, temp_mat_3, -0.5_DP)
       call sparse_embed_copy(response_kernel_jv(is), temp_mat_2)

       call sparse_embed_destroy(temp_mat_2)
       call sparse_embed_destroy(temp_mat_3)

       call sparse_embed_transpose(response_kernel_jv_CT(is), response_kernel_jv(is))
       if (response_kernel_jv(is)%iscmplx) then
          call sparse_embed_conjugate(response_kernel_jv_CT(is))
       end if
    end do
    ! end of gauge restriction
    !=======================================================

    ! =================== create remaining FO_ham terms =======================
    if (pub_paw) then
       do is = 1,pub_num_spins
          call sparse_embed_scale(FO_dij(is), 0.0_DP)
       end do
    end if
    FO_lhxc_fine = 0.0_DP
    fxc_buffer = 0.0_DP

    call sparse_embed_scale(FO_ham_jv, 0.0_DP)

    call calc_FO_ham(mdl, val_rep, val_basis, directions(perturbation,:), &
         weights(perturbation,:), dij, FO_dij, dens_fine, denskern, &
         rho_ij_dproj, rho_ij, response_kernel_jv, FO_lpseudo_fine, &
         lhxc_fine, FO_lhxc_fine, fxc_buffer, FO_ham_jv, joint_rep, &
         joint_basis, proj_basis, nl_projectors, FO_tcore_density)

    ! add previously-calculated fixed value
    call sparse_embed_axpy(FO_ham_jv, FO_ham_jv_fixed, 1.0_DP)
    ! get FO_ham from FO_ham_jv
    call sparse_embed_scale(FO_ham, 0.0_DP)
    call lr_phonons_get_basis_block(FO_ham_jv%p, FO_ham%p, joint_basis(1), &
         val_basis(1), val_basis(1), val_basis(1), 2, 1)
    !==========================================================================

    ! ==== get new Second-Order energy
    call second_order_energy(mdl, val_rep, val_basis, val_ham, joint_basis, &
         joint_rep, joint_ham, directions, weights, perturbation, dij, &
         FO_dij, FO_ham, FO_ham_jv, lhxc_fine, FO_lhxc_fine, vxc_buffer, &
         fxc_buffer, denskern%m(1,1), response_kernel_jv(1), proj_basis, &
         nl_projectors, iGpv_overlap, dipole_mat_jv, at_gamma, &
         only_SO_energy = energy, locpot_SO_mat = SO_locpot_mat)

    deallocate(FO_lhxc_fine,stat=ierr)
    call utils_dealloc_check('lr_phonons_regenerate_SO_energy','FO_lhxc_fine',&
         ierr)
    deallocate(fxc_buffer,stat=ierr)
    call utils_dealloc_check('lr_phonons_regenerate_SO_energy','fxc_buffer',&
         ierr)

    if (pub_paw.or.pub_any_nl_proj) then
       call sparse_embed_destroy(jiGp_overlap)
       call sparse_embed_destroy(pv_overlap)
       call sparse_embed_destroy(iGpv_overlap)
       call sparse_embed_destroy(viGp_overlap)
    end if

    if (pub_paw) then
       do is=1,pub_num_spins
          call sparse_embed_destroy(rho_ij_dproj(is))
          call sparse_embed_destroy(rho_ij(is))
       end do
    end if

    if (pub_paw) then
       do is=1,pub_num_spins
          call sparse_embed_destroy(FO_dij(is))
       end do
    end if


    call sparse_embed_destroy(dijps)

    call sparse_embed_destroy(FO_ham)
    call sparse_embed_destroy(FO_ham_jv)
    call sparse_embed_destroy(FO_ham_jv_fixed_temp)
    call sparse_embed_destroy(FO_ham_jv_fixed)
    call sparse_embed_destroy(FO_overlap)
    call sparse_embed_destroy(FO_overlap_jv)


    do is = 1,pub_num_spins
       call sparse_embed_destroy(response_kernel_jv(is))
       call sparse_embed_destroy(response_kernel_jv_CT(is))
    end do

    do iter1 = 1, 3
       call sparse_embed_destroy(dipole_mat_jv(iter1))
    end do

    ! Start timer
    call timer_clock('lr_phonons_regenerate_SO_energy',2)

    if (pub_debug_on_root) write(stdout, '(a)') &
         'DEBUG: Leaving lr_phonons_regenerate_SO_energy routine.'

    call comms_barrier

  end subroutine lr_phonons_regenerate_SO_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gcc32_ngwf_grad_init_precond_recip(fftbox)

    !=========================================================================!
    ! Initialises the preconditioning function in reciprocal space.           !
    !-------------------------------------------------------------------------!
    ! WARNING(gcc32): interm subroutine, will be erased as soon as I integrate!
    ! some functions back into the ngwf_gradient module                       !
    !=========================================================================!

    use constants, only: DP
    use fft_box, only: FFTBOX_INFO
    use rundat, only: pub_precond_scheme, pub_smooth_scheme, pub_k_zero
    use utils, only: utils_alloc_check, utils_assert, utils_abort

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local variables
    integer :: ierr
    integer :: i1,i2,i3,i1start
    real(kind=DP) :: scale
    real(kind=DP) :: x,temp

    ! cks: test to avoid divisions by zero and other problems
    call utils_assert(pub_k_zero > 0.0_DP, &
         'Error in ngwf_gradient_precond_recip: invalid k_zero =', pub_k_zero)

    allocate(precond_func_recip(fftbox%total_ld1,fftbox%total_ld2, &
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('gcc32_ngwf_grad_init_precond_recip', &
         'precond_func_recip',ierr)

    scale = 1.0_DP / (0.5_DP * pub_k_zero * pub_k_zero)
    precond_func_recip = 0.0_DP

    select case (pub_precond_scheme)
    case ('BG')    ! B-G method
       do i3=1,fftbox%total_pt3
          do i2=1,fftbox%total_pt2
             do i1=1,fftbox%total_pt1
                x = scale * fftbox%recip_grid(5,i1,i2,i3)
                if (pub_smooth_scheme .eq. 'NONE' .or. pub_smooth_scheme .eq. 'GAUSSIAN') then
                   precond_func_recip(i1,i2,i3) = 1.0_DP / (1.0_DP + x)
                else
                   precond_func_recip(i1,i2,i3) = sqrt(1.0_DP / (1.0_DP + x))
                endif
             end do
          end do
       end do
    case ('MAURI') ! Mauri method
       precond_func_recip(1,1,1) = 1.0_DP   ! G=0
       i1start = 2                          ! miss G=0
       do i3=1,fftbox%total_pt3
          do i2=1,fftbox%total_pt2
             do i1=i1start,fftbox%total_pt1
                x = scale * fftbox%recip_grid(5,i1,i2,i3)
                if (pub_smooth_scheme .eq. 'NONE' .or. pub_smooth_scheme .eq. 'GAUSSIAN') then
                   precond_func_recip(i1,i2,i3) = min(1.0_DP,1.0_DP/x)
                else
                   precond_func_recip(i1,i2,i3) = sqrt(min(1.0_DP,1.0_DP/x))
                endif
             end do
             i1start = 1
          end do
       end do
    case ('TETER') ! Teter-Allen-Payne method
       do i3=1,fftbox%total_pt3
          do i2=1,fftbox%total_pt2
             do i1=1,fftbox%total_pt1
                x = scale * fftbox%recip_grid(5,i1,i2,i3)
                temp = 27.0_DP + x * (18.0_DP + x * (12.0_DP + 8.0_DP * x))
                if (pub_smooth_scheme .eq. 'NONE' .or. pub_smooth_scheme .eq. 'GAUSSIAN') then
                   precond_func_recip(i1,i2,i3) = temp / (temp + 16.0_DP * x**4)
                else
                   precond_func_recip(i1,i2,i3) = sqrt(temp / (temp + 16.0_DP * x**4))
                endif
             end do
          end do
       end do
    case ('NONE')
       precond_func_recip = 1.0_DP
    case default
       call utils_abort('Error in gcc32_ngwf_grad_init_precond_recip: preconditioning&
            & scheme "'//trim(pub_precond_scheme)//'" not recognised')
    end select

  end subroutine gcc32_ngwf_grad_init_precond_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dynamical_charges_BERRY(charge, mdl, val_rep, ngwf_basis, &
       direction, weight, perturbation, val_denskern, response_denskern_jv, &
       proj_basis, nl_projectors, viGp_overlap, joint_basis, joint_rep, &
       FO_overlap)

    ! ========================================================================!
    ! This subroutine calculates the Born effective charges (also called      !
    ! dynamical charges) for both periodic and non-periodic systems           !
    ! ------------------------------------------------------------------------!
    ! For detailed explanations please read the PhD thesis of Gabriel         !
    ! Constantinescu                                                          !
    ! ------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in 2016                               !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,        !
    ! October 2018                                                            !
    ! ========================================================================!

    use datatypes, only: FUNCTIONS
    use comms, only: comms_barrier, pub_on_root, pub_my_proc_id
    use constants, only: DP, stdout, cmplx_1, cmplx_0, cmplx_i
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_brappd_ketppd
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use optics, only: optics_pos_mat_els
    use paw, only: paw_projector_overlap
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_num_spins, pub_debug_on_root, pub_spin_fac, pub_paw, &
         pub_any_nl_proj
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_scale, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_axpy, sparse_embed_transpose, &
         sparse_embed_transpose_structure, sparse_embed_product, sparse_embed_trace, &
         sparse_embed_copy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: charge(3)
    type(MODEL), intent(in) :: mdl
    type(NGWF_REP), intent(in) :: val_rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    integer, intent(in) :: direction(mdl%nat)
    real(kind=DP), intent(in) :: weight(mdl%nat)
    integer, intent(in) :: perturbation
    type(SPAM3_EMBED), intent(inout) :: val_denskern
    type(SPAM3_EMBED), intent(inout) :: response_denskern_jv
    type(FUNC_BASIS), intent(in) :: proj_basis(1)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(1)
    type(SPAM3_EMBED), intent(in) :: viGp_overlap
    type(FUNC_BASIS), intent(in) :: joint_basis(1)
    type(NGWF_REP), intent(in) :: joint_rep
    type(SPAM3_EMBED), intent(in) :: FO_overlap

    ! Local Variables

    type(SPAM3_EMBED)  :: FO_rpos(3), FO_rpos_jv(3), FO_rpos_vj(3)
    integer      :: iter, xyz, ierr, atom
    type(SPAM3_EMBED)  :: response_denskern_jv_TR, temp_mat_1
    real(kind=DP)  :: temp_charge

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &dynamical_charges_BERRY'

    ! Start timer
    call timer_clock('dynamical_charges_BERRY',1)

    ! pdh: sync procs
    call comms_barrier

    ! scale denskern and respkern to twice the number of electrons if
    ! degenerate spins
    if (pub_num_spins == 1) call sparse_embed_scale(val_denskern, 2.0_DP)
    if (pub_num_spins == 1) call sparse_embed_scale(response_denskern_jv, 2.0_DP)

    do iter = 1,3
       call sparse_embed_create(FO_rpos(iter),val_rep%overlap)
       call sparse_embed_create(FO_rpos_vj(iter), joint_rep%cross_overlap)
       call sparse_embed_transpose_structure(FO_rpos_jv(iter)%structure, &
            joint_rep%cross_overlap)
       call sparse_embed_create(FO_rpos_jv(iter))
    end do

    call calculate_dipole_mat(FO_rpos_vj, mdl, val_rep, ngwf_basis, &
         joint_rep, joint_basis, joint_rep%cross_overlap, proj_basis, &
         nl_projectors, .false.)

    call sparse_embed_transpose_structure(temp_mat_1%structure, &
         joint_rep%cross_overlap)
    call sparse_embed_create(temp_mat_1)
    call sparse_embed_transpose(temp_mat_1, joint_rep%cross_overlap)

    call calculate_dipole_mat(FO_rpos_jv, mdl, joint_rep, joint_basis, &
         val_rep, ngwf_basis, temp_mat_1, proj_basis, &
         nl_projectors, .false.)

    call sparse_embed_destroy(temp_mat_1)

    if (pub_paw) then
       call calculate_dipole_mat(FO_rpos, mdl, val_rep, ngwf_basis, &
            val_rep, ngwf_basis, val_rep%overlap, proj_basis, nl_projectors, &
            .true., viGp_overlap, viGp_overlap, FO_overlap, direction, weight)
    end if

    call sparse_embed_transpose_structure(response_denskern_jv_TR%structure, &
         response_denskern_jv)
    call sparse_embed_create(response_denskern_jv_TR)
    call sparse_embed_transpose(response_denskern_jv_TR, response_denskern_jv)

    charge = 0.0_DP

    do xyz = 1,3

       temp_charge = 0.0_DP
       call sparse_embed_trace(temp_charge, response_denskern_jv, FO_rpos_vj(xyz))
       charge(xyz) = charge(xyz) - temp_charge

       temp_charge = 0.0_DP
       call sparse_embed_trace(temp_charge, FO_rpos_jv(xyz),response_denskern_jv_TR)
       charge(xyz) = charge(xyz) - temp_charge

       temp_charge = 0.0_DP
       call sparse_embed_trace(temp_charge, FO_rpos(xyz), val_denskern)
       charge(xyz) = charge(xyz) - temp_charge

    end do


    ! THE ION PART IS ZERO SINCE WE USED R_C

    do iter = 1,3
       call sparse_embed_destroy(FO_rpos(iter))
       call sparse_embed_destroy(FO_rpos_jv(iter))
       call sparse_embed_destroy(FO_rpos_vj(iter))
    end do

    call sparse_embed_destroy(response_denskern_jv_TR)

    ! scale back denskern and respkern to half the number of electrons if
    ! degenerate spins
    if(pub_num_spins==1) call sparse_embed_scale(val_denskern, 0.5_DP)
    if(pub_num_spins==1) call sparse_embed_scale(response_denskern_jv, 0.5_DP)

    ! pdh: sync procs
    call comms_barrier

    ! Stop timer
    call timer_clock('dynamical_charges_BERRY',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &dynamical_charges_BERRY'

  end subroutine dynamical_charges_BERRY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine second_order_Ewald(force_constant_matrix, mdl, directions, &
       weights, perturbation1, qpt)

    ! ========================================================================!
    ! This subroutine calculates the second-order ion-ion (Ewald) interaction !
    ! due to atomic displacements. The equations and reasoning behind the     !
    ! implementation can be found in the PhD thesis of Gabriel Constantinescu !
    ! ------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in Dec 2016                           !
    ! Modified to remove pub_par by Joseph Prentice, October 2018             !
    ! ========================================================================!

    use comms, only: comms_barrier, pub_on_root, pub_my_proc_id, &
         comms_reduce, pub_total_num_procs
    use constants, only: DP, stdout, PI, cmplx_i, cmplx_1, cmplx_0
    use geometry, only: POINT,  operator(.CROSS.), operator(.DOT.)
    use model_type, only: MODEL
    use rundat, only: pub_debug_on_root
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_erfc

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    complex(kind=DP), intent(inout) :: force_constant_matrix( &
         mdl%nat*3+3,mdl%nat*3+3)
    integer, intent(in) :: directions(3*mdl%nat+3, mdl%nat)
    real(kind=DP), intent(in) :: weights(3*mdl%nat+3, mdl%nat)
    integer, intent(in) :: perturbation1 ! lambda
    real(kind=DP), intent(in) :: qpt(3)

    ! Local Variables
    integer :: is, ierr, iter1,iter2,iter3, ion_i, ion_j, perturbation2, &
         max_realcell_1, max_realcell_2, max_realcell_3, max_recipcell_1, &
         max_recipcell_2, max_recipcell_3, n1, n2, n3
    real(kind=DP) :: lattice_R(3), lattice_G(3), eta, a1, a2, a3, b1, b2, b3, &
         mean_a, skewness, cutoff_real, cutoff_recip, min_real, min_recip, &
         sum_precision, volume, a_val, b_val, mod_G, g_dot_rij, &
         ri(3), rj(3), r_ijR(3), mod_r_ijR, eta_r_ijR, q_dot_R

    complex(kind=DP) :: real_energy, recip_energy, temp_energy, temp_func_1, &
         temp_func_2
    real(kind=DP) :: weight_factor_recip, weight_factor_real1, &
         weight_factor_real2

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &second_order_Ewald'

    ! Start timer
    call timer_clock('second_order_Ewald',1)

    ! pdh: sync procs
    call comms_barrier

    sum_precision = 36.0_DP
    min_real = 1.0e-10_DP
    min_recip = 1.0e-10_DP

    a1 = sqrt(mdl%cell%a1%X * mdl%cell%a1%X + mdl%cell%a1%Y * mdl%cell%a1%Y + &
         mdl%cell%a1%Z * mdl%cell%a1%Z)
    a2 = sqrt(mdl%cell%a2%X * mdl%cell%a2%X + mdl%cell%a2%Y * mdl%cell%a2%Y + &
         mdl%cell%a2%Z * mdl%cell%a2%Z)
    a3 = sqrt(mdl%cell%a3%X * mdl%cell%a3%X + mdl%cell%a3%Y * mdl%cell%a3%Y + &
         mdl%cell%a3%Z * mdl%cell%a3%Z)

    b1 = sqrt(mdl%cell%b1%X * mdl%cell%b1%X + mdl%cell%b1%Y * mdl%cell%b1%Y + &
         mdl%cell%b1%Z * mdl%cell%b1%Z)
    b2 = sqrt(mdl%cell%b2%X * mdl%cell%b2%X + mdl%cell%b2%Y * mdl%cell%b2%Y + &
         mdl%cell%b2%Z * mdl%cell%b2%Z)
    b3 = sqrt(mdl%cell%b3%X * mdl%cell%b3%X + mdl%cell%b3%Y * mdl%cell%b3%Y + &
         mdl%cell%b3%Z * mdl%cell%b3%Z)

    volume = (mdl%cell%a1 .CROSS. mdl%cell%a2).DOT.mdl%cell%a3

    mean_a = (a1*a2*a3)**(1.0_DP/3.0_DP)

    eta = sqrt(PI) * (real(mdl%nat,kind=DP)**(1.0_DP/6.0_DP)) / mean_a

    cutoff_real = sqrt(sum_precision) * max(a1,a2,a3) / mean_a
    cutoff_recip = 2.0_DP * cutoff_real

    a_val = sqrt( (cutoff_real/eta)**2.0_DP + &
         (abs(mdl%cell%a1%X)+abs(mdl%cell%a2%X)+abs(mdl%cell%a3%X))**2.0_DP + &
         (abs(mdl%cell%a1%Y)+abs(mdl%cell%a2%Y)+abs(mdl%cell%a3%Y))**2.0_DP + &
         (abs(mdl%cell%a1%Z)+abs(mdl%cell%a2%Z)+abs(mdl%cell%a3%Z))**2.0_DP )

    b_val = sqrt( (cutoff_recip * 2.0_DP * eta)**2.0_DP + &
         (abs(mdl%cell%b1%X)+abs(mdl%cell%b2%X)+abs(mdl%cell%b3%X))**2.0_DP + &
         (abs(mdl%cell%b1%Y)+abs(mdl%cell%b2%Y)+abs(mdl%cell%b3%Y))**2.0_DP + &
         (abs(mdl%cell%b1%Z)+abs(mdl%cell%b2%Z)+abs(mdl%cell%b3%Z))**2.0_DP )

    max_realcell_1 = int(a_val / a1) + 1
    max_realcell_2 = int(a_val / a2) + 1
    max_realcell_3 = int(a_val / a3) + 1

    max_recipcell_1 = int(b_val / b1) + 1
    max_recipcell_2 = int(b_val / b2) + 1
    max_recipcell_3 = int(b_val / b3) + 1

    ! #########################################################################
    ! THIS SUBROUTINE COULD USE A LOT OF OPTIMISATION,
    ! BUT IT WILL SUFFICE FOR NOW
    ! ########################################################################

    ! loop over epsilon perturbation
    do perturbation2 = 1, 3*mdl%nat

       real_energy = cmplx_0
       recip_energy = cmplx_0

       ! MPI-parallelise sum over ion_i
       do ion_i = pub_my_proc_id + 1, mdl%nat, pub_total_num_procs

          do ion_j = 1, mdl%nat

             ri(1) = mdl%elements(ion_i)%centre%X
             ri(2) = mdl%elements(ion_i)%centre%Y
             ri(3) = mdl%elements(ion_i)%centre%Z

             rj(1) = mdl%elements(ion_j)%centre%X
             rj(2) = mdl%elements(ion_j)%centre%Y
             rj(3) = mdl%elements(ion_j)%centre%Z

             ! RECIP-sum term
             do n1 = -max_recipcell_1, max_recipcell_1
                do n2 = -max_recipcell_2, max_recipcell_2
                   do n3 = -max_recipcell_3, max_recipcell_3

                      ! first part
                      lattice_G(1)  = qpt(1) + &
                           real(n1,kind=DP) * mdl%cell%b1%X + &
                           real(n2,kind=DP) * mdl%cell%b2%X + &
                           real(n3,kind=DP) * mdl%cell%b3%X
                      lattice_G(2)  = qpt(2) + &
                           real(n1,kind=DP) * mdl%cell%b1%Y + &
                           real(n2,kind=DP) * mdl%cell%b2%Y + &
                           real(n3,kind=DP) * mdl%cell%b3%Y
                      lattice_G(3)  = qpt(3) + &
                           real(n1,kind=DP) * mdl%cell%b1%Z + &
                           real(n2,kind=DP) * mdl%cell%b2%Z + &
                           real(n3,kind=DP) * mdl%cell%b3%Z
                      mod_G = sqrt( lattice_G(1) * lattice_G(1) + &
                           lattice_G(2) * lattice_G(2) + lattice_G(3) * &
                           lattice_G(3) )
                      ! safety measure:
                      ! if (mod_G .le. 1.0E-10_DP) mod_G = 0.0_DP

                      mod_G = mod_G / (2.0_DP * eta)

                      g_dot_rij = lattice_G(1) * (ri(1) - rj(1)) + &
                           lattice_G(2) * (ri(2) - rj(2)) + &
                           lattice_G(3) * (ri(3) - rj(3))


                      if ((mod_G.le.cutoff_recip).and.(mod_G.ge.min_recip)) then

                         weight_factor_recip = &
                              lattice_G(directions(perturbation2, ion_j)) * &
                              lattice_G(directions(perturbation1, ion_i)) * &
                              weights(perturbation2, ion_j) * &
                              weights(perturbation1, ion_i) - &
                              lattice_G(directions(perturbation2, ion_i)) * &
                              lattice_G(directions(perturbation1, ion_i)) * &
                              weights(perturbation2, ion_i) * &
                              weights(perturbation1, ion_i)


                         temp_energy = exp(-mod_G*mod_G) * &
                              weight_factor_recip * &
                              exp(-cmplx_i*g_dot_rij) / &
                              (4.0_DP*eta*eta*mod_G*mod_G)

                         recip_energy = recip_energy + temp_energy * 2.0_DP * &
                              PI * mdl%elements(ion_i)%ion_charge * &
                              mdl%elements(ion_j)%ion_charge / volume

                      end if

                      ! second part
                      lattice_G(1)  = real(n1,kind=DP) * mdl%cell%b1%X + &
                           real(n2,kind=DP) * mdl%cell%b2%X + &
                           real(n3,kind=DP) * mdl%cell%b3%X
                      lattice_G(2)  = real(n1,kind=DP) * mdl%cell%b1%Y + &
                           real(n2,kind=DP) * mdl%cell%b2%Y + &
                           real(n3,kind=DP) * mdl%cell%b3%Y
                      lattice_G(3)  = real(n1,kind=DP) * mdl%cell%b1%Z + &
                           real(n2,kind=DP) * mdl%cell%b2%Z + &
                           real(n3,kind=DP) * mdl%cell%b3%Z
                      mod_G = sqrt( lattice_G(1) * lattice_G(1) + &
                           lattice_G(2) * lattice_G(2) + lattice_G(3) * &
                           lattice_G(3) )
                      mod_G = mod_G / (2.0_DP * eta)
                      ! safety measure:
                      ! if (mod_G .le. 1.0E-10_DP) mod_G = 0.0_DP


                      g_dot_rij = lattice_G(1) * (ri(1) - rj(1)) + &
                           lattice_G(2) * (ri(2) - rj(2)) + &
                           lattice_G(3) * (ri(3) - rj(3))

                      if ((mod_G.le.cutoff_recip).and.(mod_G.ge.min_recip)) then

                         weight_factor_recip = &
                              lattice_G(directions(perturbation2, ion_j)) * &
                              lattice_G(directions(perturbation1, ion_i)) * &
                              weights(perturbation2, ion_j) * &
                              weights(perturbation1, ion_i) - &
                              lattice_G(directions(perturbation2, ion_i)) * &
                              lattice_G(directions(perturbation1, ion_i)) * &
                              weights(perturbation2, ion_i) * &
                              weights(perturbation1, ion_i)


                         temp_energy = exp(-mod_G*mod_G) * &
                              weight_factor_recip * &
                              exp(-cmplx_i*g_dot_rij) / &
                              (4.0_DP*eta*eta*mod_G*mod_G)

                         recip_energy = recip_energy + temp_energy * 2.0_DP * &
                              PI * mdl%elements(ion_i)%ion_charge * &
                              mdl%elements(ion_j)%ion_charge / volume

                      end if

                   end do ! n3
                end do ! n2
             end do ! n1

             !================================================================

             ! REAL-sum term
             do n1 = -max_realcell_1, max_realcell_1
                do n2 = -max_realcell_2, max_realcell_2
                   do n3 = -max_realcell_3, max_realcell_3

                      lattice_R(1)  = real(n1,kind=DP) * mdl%cell%a1%X + &
                           real(n2,kind=DP) * mdl%cell%a2%X + &
                           real(n3,kind=DP) * mdl%cell%a3%X
                      lattice_R(2)  = real(n1,kind=DP) * mdl%cell%a1%Y + &
                           real(n2,kind=DP) * mdl%cell%a2%Y + &
                           real(n3,kind=DP) * mdl%cell%a3%Y
                      lattice_R(3)  = real(n1,kind=DP) * mdl%cell%a1%Z + &
                           real(n2,kind=DP) * mdl%cell%a2%Z + &
                           real(n3,kind=DP) * mdl%cell%a3%Z

                      r_ijR(1) = ri(1) - rj(1) + lattice_R(1)
                      r_ijR(2) = ri(2) - rj(2) + lattice_R(2)
                      r_ijR(3) = ri(3) - rj(3) + lattice_R(3)

                      mod_r_ijR = sqrt( r_ijR(1)**2.0_DP + r_ijR(2)**2.0_DP + &
                           r_ijR(3)**2.0_DP )

                      ! safety measure:
                      ! if (mod_r_ijR .le. 1.0E-10_DP) mod_r_ijR = 0.0_DP

                      eta_r_ijR = eta * mod_r_ijR

                      if ( (eta_r_ijR .le. cutoff_real) .and. &
                           (eta_r_ijR .ge. min_real) ) then


                         weight_factor_real1 = 0.0_DP
                         if (directions(perturbation2,ion_j) == &
                             directions(perturbation1,ion_i)) then
                             weight_factor_real1 = weight_factor_real1 + &
                                  weights(perturbation2,ion_j) * &
                                  weights(perturbation1,ion_i)
                         end if

                         if (directions(perturbation2,ion_i) == &
                             directions(perturbation1,ion_i)) then
                             weight_factor_real1 = weight_factor_real1 - &
                                  weights(perturbation2,ion_i) * &
                                  weights(perturbation1,ion_i)
                         end if


                         weight_factor_real2 = &
                              r_ijR(directions(perturbation2, ion_j)) * &
                              r_ijR(directions(perturbation1, ion_i)) * &
                              weights(perturbation2, ion_j) * &
                              weights(perturbation1, ion_i) - &
                              r_ijR(directions(perturbation2, ion_i)) * &
                              r_ijR(directions(perturbation1, ion_i)) * &
                              weights(perturbation2, ion_i) * &
                              weights(perturbation1, ion_i)



                         temp_func_1 = utils_erfc(eta_r_ijR) + 2.0_DP * &
                              (eta_r_ijR / sqrt(PI)) * exp(-eta_r_ijR*eta_r_ijR)

                         temp_func_2 = 3.0_DP*utils_erfc(eta_r_ijR) + 2.0_DP * &
                              (eta_r_ijR / sqrt(PI))*exp(-eta_r_ijR*eta_r_ijR)*&
                              (3.0_DP + 2.0_DP*eta_r_ijR*eta_r_ijR)

                         temp_func_1 = temp_func_1 / (2.0_DP*mod_r_ijR**3.0_DP)
                         temp_func_2 = temp_func_2 / (2.0_DP*mod_r_ijR**5.0_DP)

                         q_dot_R = qpt(1) * lattice_R(1) + &
                              qpt(2) * lattice_R(2) + &
                              qpt(3) * lattice_R(3)


                         ! first_part
                         temp_energy = weight_factor_real1 * temp_func_1 * &
                              (1.0_DP + exp(cmplx_i*q_dot_R))

                         temp_energy = temp_energy - weight_factor_real2 * &
                              temp_func_2 * (1.0_DP + exp(cmplx_i*q_dot_R))

                         real_energy = real_energy + temp_energy * &
                              mdl%elements(ion_i)%ion_charge * &
                              mdl%elements(ion_j)%ion_charge

                      end if

                   end do ! n3
                end do ! n2
             end do ! n1

          end do ! over ion_j

       end do ! over ion_i

       call comms_barrier
       call comms_reduce('SUM',real_energy)
       call comms_reduce('SUM',recip_energy)

       ! add energies to force constant
       force_constant_matrix(perturbation1, perturbation2) = &
            force_constant_matrix(perturbation1, perturbation2) + &
            real_energy + recip_energy

    end do ! over perturbation2


    ! pdh: sync procs
    call comms_barrier

    ! Stop timer
    call timer_clock('second_order_Ewald',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving second_order_Ewald'

  end subroutine second_order_Ewald

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_phonons_born_asr(Born_charges, mdl)

    ! ========================================================================!
    ! This subroutine enforces onto the Born charges the acoustic and         !
    ! rotational (if applicable) sum rules, according to the MSc thesis of    !
    ! Nicolas Mounet ("Structural, vibrational and thermodynamic properties of!
    ! carbon allotropes from first-principles: diamond, graphite, and         !
    ! nanotubes"), and as implemented in Quantum Espresso.                    !
    ! ------------------------------------------------------------------------!
    ! ------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in Dec 2016/ Jan 2017                 !
    ! Modified to remove pub_par by Joseph Prentice, October 2018             !
    ! ========================================================================!

    use constants, only: DP, TWO_PI
    use geometry, only: operator(.DOT.), POINT
    use model_type, only : MODEL
    use rundat, only: pub_lr_phonons_zero_dim
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! arguments
    real(kind=DP), dimension(:,:), intent(inout)  :: Born_charges
    type(MODEL) :: mdl

    ! Local variables

    real(kind=DP), allocatable :: local_Born_charges(:,:,:), local_U(:,:,:,:), &
         local_X(:,:,:), local_W(:,:,:), atom_pos(:,:)
    integer :: local_less(18)
    real(kind=DP) :: val
    integer :: i, j, p, k, iless, nless, r, iter, ierr, q

    call timer_clock('lr_phonons_born_asr',1)

    allocate(atom_pos(3,1:mdl%nat),stat=ierr)
    call utils_alloc_check('lr_phonons_born_asr','atom_pos',ierr)
    allocate(local_Born_charges(3,3,1:mdl%nat),stat=ierr)
    call utils_alloc_check('lr_phonons_born_asr','local_Born_charges',ierr)
    allocate(local_U(18,3,3,1:mdl%nat),stat=ierr)
    call utils_alloc_check('lr_phonons_born_asr','local_U',ierr)
    allocate(local_W(3,3,1:mdl%nat),stat=ierr)
    call utils_alloc_check('lr_phonons_born_asr','local_W',ierr)
    allocate(local_X(3,3,1:mdl%nat),stat=ierr)
    call utils_alloc_check('lr_phonons_born_asr','local_X',ierr)


    local_Born_charges = 0.0_DP
    local_U = 0.0_DP
    local_X = 0.0_DP
    local_W = 0.0_DP
    local_less = 0

    ! tranfer Born charges to local copy
    do iter = 1, 3*mdl%nat
       i = ((iter - 1) / 3) + 1
       j = mod((iter - 1), 3) + 1
       local_Born_charges(:, j, i) = Born_charges(:,iter)
    end do

    do iter = 1,mdl%nat
       atom_pos(1,iter) = mdl%elements(iter)%centre%X
       atom_pos(2,iter) = mdl%elements(iter)%centre%Y
       atom_pos(3,iter) = mdl%elements(iter)%centre%Z
    end do

    p = 0
    ! acoustical branch first
    do i = 1, 3
       do j = 1,3
          p = p +1
          local_U(p,i,j,:) = 1.0_DP
       end do
    end do

    ! if zero-dimensional case, deal with rotational modes
    if (pub_lr_phonons_zero_dim) then
       do i = 1,3
          do j = 1,3
             p = p + 1
             local_U(p,i,mod(j,3)+1,:) = -atom_pos(mod(j+1,3)+1,:)
             local_U(p,i,mod(j+1,3)+1,:) = atom_pos(mod(j,3)+1,:)

          end do
       end do
    end if

    ! perform Gram-Schmidt orthogonalisation

    nless = 0
    do k = 1, p
       local_W = local_U(k,:,:,:)
       local_X = local_U(k,:,:,:)
       do q = 1, k-1
          r = 1
          do iless = 1, nless
             if (local_less(iless) == q) r = 0
          end do
          if (r.ne.0) then
             local_W = local_W - local_U(q,:,:,:) * sum( local_X(:,:,:) * &
                  local_U(q,:,:,:) )
          end if
       end do

       if (sum(local_W*local_W).ge.1.0E-16_DP) then
          local_U(k,:,:,:) = local_W(:,:,:) / sqrt(sum(local_W*local_W))
       else
          nless = nless + 1
          local_less(nless) = k
       end if
    end do


    ! perform projection
    local_W = 0.0_DP
    do k = 1, p
       r = 1
       do iless = 1, nless
          if (local_less(iless) == k) r = 0
       end do

       if (r.ne.0) then
          local_X(:,:,:) = local_U(k,:,:,:)
          local_W = local_W + local_U(k,:,:,:) * sum(local_X*local_Born_charges)
       end if
    end do

    local_Born_charges = local_Born_charges - local_W

    ! plug the modified values back in the original Born charges routine
    Born_charges = 0.0_DP
    do iter = 1, 3*mdl%nat
       i = ((iter - 1) / 3) + 1
       j = mod((iter - 1), 3) + 1
       Born_charges(:,iter) = local_Born_charges(:, j, i)
    end do

    deallocate(atom_pos,stat=ierr)
    call utils_dealloc_check('lr_phonons_born_asr','atom_pos',ierr)
    deallocate(local_Born_charges,stat=ierr)
    call utils_dealloc_check('lr_phonons_born_asr','local_Born_charges',ierr)
    deallocate(local_U,stat=ierr)
    call utils_dealloc_check('lr_phonons_born_asr','local_U',ierr)
    deallocate(local_W,stat=ierr)
    call utils_dealloc_check('lr_phonons_born_asr','local_W',ierr)
    deallocate(local_X,stat=ierr)
    call utils_dealloc_check('lr_phonons_born_asr','local_X',ierr)

    call timer_clock('lr_phonons_born_asr',2)

  end subroutine lr_phonons_born_asr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_phonons_dynmat(force_constant_matrix, eigenvecs, eigenvals,&
       current_qpt, at_gamma, dielectric_tensor, Born_charges, volume, mdl)

    ! ========================================================================!
    ! 1) If the calculation is at q = 0, one enforces the acoustic and        !
    ! rotational (if applicable) sum rules, according to the MSc thesis of    !
    ! Nicolas Mounet ("Structural, vibrational and thermodynamic properties of!
    ! carbon allotropes from first-principles: diamond, graphite, and         !
    ! nanotubes"), and as implemented in Quantum Espresso.                    !
    ! ------------------------------------------------------------------------!
    ! 2) If the calculation is NOT at q = 0, one adds the non-analytical part !
    ! of the force constant matrix, allowing the correct representation of the!
    ! LO-TO mode splitting, according to Eq. 60 of X. Gonze and C. Lee        !
    ! PRB (1997) 55, 10355-10368.                                             !
    ! ------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in Dec 2016/ Jan 2017                 !
    ! Modified to remove pub_par by Joseph Prentice, October 2018             !
    ! ========================================================================!

    use comms, only: comms_barrier, pub_on_root
    use constants, only: DP, cmplx_1, cmplx_0, TWO_PI, periodic_table_mass, &
         stdout, TWO_PI, PI
    use dense, only: DEM, dense_create, dense_destroy, dense_eigensolve, &
         dense_scale, dense_put_element, dense_get_element
    use geometry, only: operator(.DOT.), POINT
    use model_type, only: MODEL
    use rundat, only: pub_eigensolver_orfac, pub_eigensolver_abstol, &
         pub_lr_phonons_zero_dim
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    complex(kind=DP), intent(inout) :: force_constant_matrix(3*mdl%nat, &
         3*mdl%nat)
    complex(kind=DP), intent(inout) :: eigenvecs(3*mdl%nat,3*mdl%nat)
    real(kind=DP), intent(inout) :: eigenvals(3*mdl%nat)
    real(kind=DP), intent(in) :: current_qpt(3)
    logical, intent(in) :: at_gamma
    real(kind=DP), intent(in) :: dielectric_tensor(3,3)
    real(kind=DP), intent(in) :: Born_charges(3,3*mdl%nat)
    real(kind=DP), intent(in) :: volume

    ! Local variables
    integer :: acoustics(3), ierr, iter1, iter2, atom_1, atom_2, &
         dir_1, dir_2
    real(kind=DP) :: element, mass_1, mass_2
    complex(kind=DP) :: element_cmplx
    complex(kind=DP), allocatable :: force_constant_matrix_temp(:,:), &
         dynamical_matrix(:,:)
    type(DEM) :: dyn_matrix_dem, unit_mat_dem, eigenvecs_dem

    real(kind=DP), allocatable :: local_fcm(:,:,:,:), local_w(:,:,:,:), &
         local_x(:,:,:,:), local_u(:,:,:,:,:), local_v(:,:), atom_pos(:,:), &
         local_born(:,:,:)
    integer, allocatable :: local_ind_v(:,:,:), local_u_less(:)
    integer :: i,j,k,p,q,r,m,l,nless,iless
    real(kind=DP) :: scal, norm2

    call timer_clock('lr_phonons_dynmat', 1)

    call dense_create(dyn_matrix_dem, 3*mdl%nat, 3*mdl%nat, iscmplx = .true.)
    call dense_create(unit_mat_dem, 3*mdl%nat, 3*mdl%nat, iscmplx = .true.)
    call dense_create(eigenvecs_dem, 3*mdl%nat, 3*mdl%nat, iscmplx = .true.)

    allocate(atom_pos(3,1:mdl%nat),stat=ierr)
    call utils_alloc_check('lr_phonons_born_asr','atom_pos',ierr)

    allocate(dynamical_matrix(3*mdl%nat,3*mdl%nat), stat=ierr)
    call utils_alloc_check('lr_phonons_dynmat', 'dynamical_matrix', ierr)
    allocate(force_constant_matrix_temp(3*mdl%nat,3*mdl%nat), stat=ierr)
    call utils_alloc_check('lr_phonons_dynmat', &
         'force_constant_matrix_temp', ierr)
    allocate(local_fcm(3,3,mdl%nat,mdl%nat), stat=ierr)
    call utils_alloc_check('lr_phonons_dynmat', 'local_fcm', ierr)
    allocate(local_born(3,3,mdl%nat), stat=ierr)
    call utils_alloc_check('lr_phonons_dynmat', 'local_born', ierr)
    allocate(local_x(3,3,mdl%nat,mdl%nat), stat=ierr)
    call utils_alloc_check('lr_phonons_dynmat', 'local_x', ierr)
    allocate(local_w(3,3,mdl%nat,mdl%nat), stat=ierr)
    call utils_alloc_check('lr_phonons_dynmat', 'local_w', ierr)
    allocate(local_u(18*mdl%nat,3,3,mdl%nat,mdl%nat), stat=ierr)
    call utils_alloc_check('lr_phonons_dynmat', 'local_u', ierr)
    allocate(local_v(9*mdl%nat*mdl%nat,2), stat=ierr)
    call utils_alloc_check('lr_phonons_dynmat', 'local_v', ierr)
    allocate(local_ind_v(9*mdl%nat*mdl%nat,2,4), stat=ierr)
    call utils_alloc_check('lr_phonons_dynmat', 'local_ind_v', ierr)
    allocate(local_u_less(18*mdl%nat), stat=ierr)
    call utils_alloc_check('lr_phonons_dynmat', 'local_u_less', ierr)

    call comms_barrier

    force_constant_matrix_temp(:,:) = cmplx_0
    ! ENSURE FORCE CONSTANT MATRIX IS HERMITIAN
    do iter1 = 1,3*mdl%nat
       do iter2 = iter1, 3*mdl%nat ! start from iter1 to include diagonal terms

          element_cmplx = 0.5_DP * force_constant_matrix(iter1, iter2) + &
               0.5_DP * conjg(force_constant_matrix(iter2,iter1))
          force_constant_matrix_temp(iter1, iter2) = element_cmplx
          force_constant_matrix_temp(iter2, iter1) = conjg(element_cmplx)
       end do
    end do
    force_constant_matrix(:,:) = force_constant_matrix_temp(:,:)

    ! SCALE FORCE CONSTANT MATRICES ACCORDING TO ATOM MASSES -> DYNAMICAL MATRIX
    do iter1 = 1, 3*mdl%nat
       do iter2 = 1, 3*mdl%nat
          atom_1 = ((iter1-1) / 3) + 1
          atom_2 = ((iter2-1) / 3) + 1
          mass_1 = periodic_table_mass(mdl%elements(atom_1)%atomic_number) / &
               electron_mass_u
          mass_2 = periodic_table_mass(mdl%elements(atom_2)%atomic_number) / &
               electron_mass_u
          dynamical_matrix(iter1,iter2) = &
               force_constant_matrix(iter1,iter2) / sqrt(mass_1*mass_2)

       end do
    end do

    call comms_barrier

    ! ########################################################################
    ! CONVERT FORCE CONSTANT / DYNAMICAL MATRIX TO DENSE FORMAT
    ! AFTERWARDS, DIAGONALISE AND GET FREQUENCIES

    do iter1 = 1,3*mdl%nat
       do iter2 = 1,3*mdl%nat
          call dense_put_element(dynamical_matrix(iter1,iter2), &
               dyn_matrix_dem, iter1, iter2)
          if (iter1 == iter2) then
             call dense_put_element(cmplx_1, unit_mat_dem, iter1, iter2)
          else
             call dense_put_element(cmplx_0, unit_mat_dem, iter1, iter2)
          end if
       end do
    end do

    call comms_barrier

    eigenvals = 0.0_DP
    call dense_scale(eigenvecs_dem, cmplx_0)

    call dense_eigensolve(3*mdl%nat, eigenvals, dyn_matrix_dem, &
         unit_mat_dem, 1, eigenvecs_dem, pub_eigensolver_orfac, &
         pub_eigensolver_abstol)

    ! copy eigenvectors to complex arrays from DEM format
    do iter2 = 1, 3*mdl%nat
       do iter1 = 1, 3*mdl%nat

          element_cmplx = cmplx_0
          call dense_get_element(element_cmplx, eigenvecs_dem, iter1, iter2)

          eigenvecs(iter1, iter2) = element_cmplx

       end do
    end do
    ! ########################################################################


    ! print out phonon frequencies at this qpt
    if (pub_on_root) then
       write(stdout,'()')
       write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
            &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
       write(stdout,'(a)') 'Phonon frequencies BEFORE ASR    (cm^-1)     &
            &            (THz)'
       write(stdout,'()')
       do iter1 = 1,3*mdl%nat
          element = sign(sqrt(abs(eigenvals(iter1))) / TWO_PI, &
               eigenvals(iter1))
          write(stdout,'(12x,i5,a,2(6x,f16.11))') iter1, ":", &
               element*au2inv_cm, element*au2THz
       end do
       write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
            &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
    end if

    ! ########################################################################


    if (at_gamma) then
    !-------------------------------------------------------------------------!
    !------------------- APPLY ASR CORRECTION IF AT GAMMA --------------------!
    !-------------------------------------------------------------------------!

       ! move real part of dynamical matrix to local_fcm
       do iter1 = 1,3*mdl%nat
          do iter2 = 1,3*mdl%nat
             atom_1 = ((iter1 - 1) / 3) + 1
             dir_1 = mod((iter1 - 1), 3) + 1
             atom_2 = ((iter2 - 1) / 3) + 1
             dir_2 = mod((iter2 - 1), 3) + 1

             local_fcm(dir_1,dir_2,atom_1,atom_2) = &
                  real(force_constant_matrix(iter1,iter2), kind=DP)
          end do
       end do


       do iter1 = 1,mdl%nat
          atom_pos(1,iter1) = mdl%elements(iter1)%centre%X
          atom_pos(2,iter1) = mdl%elements(iter1)%centre%Y
          atom_pos(3,iter1) = mdl%elements(iter1)%centre%Z
       end do

       local_u = 0.0_DP

       p = 0
       ! acoustical branch
       do dir_1 = 1,3
          do dir_2 = 1,3
             do atom_1 = 1,mdl%nat
                p = p + 1
                do atom_2 = 1,mdl%nat
                   local_u(p,dir_1,dir_2,atom_1,atom_2) = 1.0_DP
                end do
             end do
          end do
       end do

       ! if we are dealing with a molecule, consider rotational modes also
       if (pub_lr_phonons_zero_dim) then
          do dir_1 = 1,3
             do dir_2 = 1,3
                do atom_1 = 1,mdl%nat
                   p=p+1
                   do atom_2 = 1,mdl%nat
                      local_u(p,dir_1,dir_2,atom_1,atom_2) = 0.0_DP
                      local_u(p,dir_1,MOD(dir_2,3)+1,atom_1,atom_2) = &
                           -atom_pos(MOD(dir_2+1,3)+1,atom_2)
                      local_u(p,dir_1,MOD(dir_2+1,3)+1,atom_1,atom_2) = &
                           atom_pos(MOD(dir_2,3)+1,atom_2)
                   end do
                end do
             end do
          end do
       end if

       ! vectors associated to symmetry constraints
       m = 0
       do dir_1 = 1,3
          do dir_2 = 1,3
             do atom_1 = 1,mdl%nat
                do atom_2 = 1,mdl%nat
                   ! These are the vectors associated with the symmetry constraints
                   q = 1
                   l = 1
                   do while((l.le.m).and.(q.ne.0))
                      if ((local_ind_v(l,1,1).eq.dir_1) .and. &
                          (local_ind_v(l,1,2).eq.dir_2) .and. &
                          (local_ind_v(l,1,3).eq.atom_1) .and. &
                          (local_ind_v(l,1,4).eq.atom_2)) q=0

                      if ((local_ind_v(l,2,1).eq.dir_1) .and. &
                          (local_ind_v(l,2,2).eq.dir_2) .and. &
                          (local_ind_v(l,2,3).eq.atom_1) .and. &
                          (local_ind_v(l,2,4).eq.atom_2)) q=0

                      l=l+1
                   end do
                   if ((dir_1.eq.dir_2).and.(atom_1.eq.atom_2)) q=0
                   if (q.ne.0) then
                      m=m+1
                      local_ind_v(m,1,1) = dir_1
                      local_ind_v(m,1,2) = dir_2
                      local_ind_v(m,1,3) = atom_1
                      local_ind_v(m,1,4) = atom_2
                      local_v(m,1) = 1.0_DP / sqrt(2.0_DP)
                      local_ind_v(m,2,1) = dir_2
                      local_ind_v(m,2,2) = dir_1
                      local_ind_v(m,2,3) = atom_2
                      local_ind_v(m,2,4) = atom_1
                      local_v(m,2)= -1.0_DP / sqrt(2.0_DP)
                   end if
                end do
             end do
          end do
       end do

       ! Gram-Schmidt orthogonalisation
       nless = 0
       do k = 1,p
          local_w(:,:,:,:) = local_u(k,:,:,:,:)
          local_x(:,:,:,:) = local_u(k,:,:,:,:)
          do l = 1,m
             call scal_prod2(local_x,local_v(l,:),local_ind_v(l,:,:),scal)
             do r=1,2
                dir_1 = local_ind_v(l,r,1)
                dir_2 = local_ind_v(l,r,2)
                atom_1 = local_ind_v(l,r,3)
                atom_2 = local_ind_v(l,r,4)
                local_w(dir_1,dir_2,atom_1,atom_2) = &
                     local_w(dir_1,dir_2,atom_1,atom_2) - scal * local_v(l,r)
             end do
          end do
          if (k.le.(9*mdl%nat)) then
             atom_1=MOD(k,mdl%nat)
             if (atom_1.eq.0) atom_1=mdl%nat
             dir_2=MOD((k-atom_1)/mdl%nat,3)+1
             dir_1=MOD((((k-atom_1)/mdl%nat)-dir_2+1)/3,3)+1
          else
             q=k-9*mdl%nat
             atom_1=MOD(q,mdl%nat)
             if (atom_1.eq.0) atom_1=mdl%nat
             dir_2=MOD((q-atom_1)/mdl%nat,3)+1
             dir_1=MOD((((q-atom_1)/mdl%nat)-dir_2+1)/3,3)+1
          endif
          do q=1,k-1
             r=1
             do iless=1,nless
                if (local_u_less(iless).eq.q) r=0
             end do
             if (r.ne.0) then
                scal = sum( local_x(dir_1,:,atom_1,:) * &
                     local_u(q,dir_1,:,atom_1,:))
                local_w(:,:,:,:) = local_w(:,:,:,:) - scal * local_u(q,:,:,:,:)
             end if
          end do
          norm2 = sum(local_w * local_w)
          if (norm2.gt.1.0E-16_DP) then
             local_u(k,:,:,:,:) = local_w(:,:,:,:) / sqrt(norm2)
          else
             nless = nless+1
             local_u_less(nless) = k
          end if
       end do


       ! project on orthogonal subspace
       local_w(:,:,:,:) = 0.0_DP
       do l = 1,m
          call scal_prod2(local_fcm(:,:,:,:),local_v(l,:),local_ind_v(l,:,:), &
               scal)
          do r=1,2
             dir_1 = local_ind_v(l,r,1)
             dir_2 = local_ind_v(l,r,2)
             atom_1 = local_ind_v(l,r,3)
             atom_2 = local_ind_v(l,r,4)
             local_w(dir_1,dir_2,atom_1,atom_2) = &
                  local_w(dir_1,dir_2,atom_1,atom_2) + scal*local_v(l,r)
          end do
       end do
       do k=1,p
          r=1
          do iless=1,nless
             if (local_u_less(iless).eq.k) r=0
          end do
          if (r.ne.0) then
             local_x(:,:,:,:) = local_u(k,:,:,:,:)
             scal = sum(local_x * local_fcm)
             local_w(:,:,:,:) = local_w(:,:,:,:) + scal * local_u(k,:,:,:,:)
          end if
       end do

      local_fcm = local_fcm - local_w


       ! replace real part of force constant matrix with local_fcm
       do iter1 = 1,3*mdl%nat
          do iter2 = 1,3*mdl%nat
             atom_1 = ((iter1 - 1) / 3) + 1
             dir_1 = mod((iter1 - 1), 3) + 1
             atom_2 = ((iter2 - 1) / 3) + 1
             dir_2 = mod((iter2 - 1), 3) + 1

             force_constant_matrix(iter1,iter2) = &
                  cmplx(local_fcm(dir_1,dir_2,atom_1,atom_2), &
                       aimag(force_constant_matrix(iter1,iter2)), kind=DP)
          end do
       end do

    else ! not at GAMMA
    !-------------------------------------------------------------------------!
    !------------------- APPLY LO-TO CORRECTION IF NOT AT GAMMA --------------!
    !-------------------------------------------------------------------------!
       local_fcm = 0.0_DP

       ! tranfer Born charges to local copy
       do iter1 = 1, 3*mdl%nat
          atom_1 = ((iter1 - 1) / 3) + 1
          dir_1 = mod((iter1 - 1), 3) + 1
          local_born(:, dir_1, atom_1) = Born_charges(:,iter1)
       end do

       ! calculate denominator
       scal = sum(matmul(dielectric_tensor, current_qpt) * current_qpt)

       do atom_1 = 1, mdl%nat
          do dir_1 = 1,3
             do atom_2 = 1, mdl%nat
                do dir_2 = 1,3
                   local_fcm(dir_1,dir_2,atom_1,atom_2) = (4.0_DP*PI/volume) * &
                        sum(current_qpt(:) * local_born(:, dir_1, atom_1)) * &
                        sum(current_qpt(:) * local_born(:, dir_2, atom_2)) / &
                        scal
                end do
             end do
          end do
       end do

       ! add local_fcm to real part of force constant matrix
       do iter1 = 1,3*mdl%nat
          do iter2 = 1,3*mdl%nat
             atom_1 = ((iter1 - 1) / 3) + 1
             dir_1 = mod((iter1 - 1), 3) + 1
             atom_2 = ((iter2 - 1) / 3) + 1
             dir_2 = mod((iter2 - 1), 3) + 1

             force_constant_matrix(iter1,iter2) = &
                  force_constant_matrix(iter1,iter2) + &
                  cmplx(local_fcm(dir_1,dir_2,atom_1,atom_2), 0.0_DP, kind=DP)
          end do
       end do


    end if ! if at q=0 or not
    ! ########################################################################

    force_constant_matrix_temp(:,:) = cmplx_0
    ! ENSURE FORCE CONSTANT MATRIX IS HERMITIAN
    do iter1 = 1,3*mdl%nat
       do iter2 = iter1, 3*mdl%nat ! start from iter1 to include diagonal terms

          element_cmplx = 0.5_DP * force_constant_matrix(iter1, iter2) + &
               0.5_DP * conjg(force_constant_matrix(iter2,iter1))
          force_constant_matrix_temp(iter1, iter2) = element_cmplx
          force_constant_matrix_temp(iter2, iter1) = conjg(element_cmplx)
       end do
    end do
    force_constant_matrix(:,:) = force_constant_matrix_temp(:,:)

    ! SCALE FORCE CONSTANT MATRICES ACCORDING TO ATOM MASSES -> DYNAMICAL MATRIX
    do iter1 = 1, 3*mdl%nat
       do iter2 = 1, 3*mdl%nat
          atom_1 = ((iter1-1) / 3) + 1
          atom_2 = ((iter2-1) / 3) + 1
          mass_1 = periodic_table_mass(mdl%elements(atom_1)%atomic_number) / &
               electron_mass_u
          mass_2 = periodic_table_mass(mdl%elements(atom_2)%atomic_number) / &
               electron_mass_u
          dynamical_matrix(iter1,iter2) = &
               force_constant_matrix(iter1,iter2) / sqrt(mass_1*mass_2)

       end do
    end do

    ! CONVERT DYNAMICAL MATRIX TO DENSE FORM, DIAGONALISE AND GET FREQUENCIES

    do iter1 = 1,3*mdl%nat
       do iter2 = 1,3*mdl%nat
          call dense_put_element(dynamical_matrix(iter1,iter2), &
               dyn_matrix_dem, iter1, iter2)
          if (iter1 == iter2) then
             call dense_put_element(cmplx_1, unit_mat_dem, iter1, iter2)
          else
             call dense_put_element(cmplx_0, unit_mat_dem, iter1, iter2)
          end if
       end do
    end do

    call comms_barrier

    eigenvals = 0.0_DP
    call dense_scale(eigenvecs_dem, cmplx_0)

    call dense_eigensolve(3*mdl%nat, eigenvals, dyn_matrix_dem, &
         unit_mat_dem, 1, eigenvecs_dem, pub_eigensolver_orfac, &
         pub_eigensolver_abstol)

    ! copy eigenvectors to complex arrays from DEM format
    do iter2 = 1, 3*mdl%nat
       do iter1 = 1, 3*mdl%nat

          element_cmplx = cmplx_0
          call dense_get_element(element_cmplx, eigenvecs_dem, iter1, iter2)

          eigenvecs(iter1, iter2) = element_cmplx

       end do
    end do
    ! ########################################################################


    ! print out phonon frequencies at this qpt
    if (pub_on_root) then
       write(stdout,'()')
       write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
            &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
       write(stdout,'(a)') 'Phonon frequencies AFTER ASR / LO-TO (cm^-1)   &
            &          (THz)'
       write(stdout,'()')
       do iter1 = 1,3*mdl%nat
          element = sign(sqrt(abs(eigenvals(iter1))) / TWO_PI, &
               eigenvals(iter1))
          write(stdout,'(12x,i5,a,2(6x,f16.11))') iter1, ":", &
               element*au2inv_cm, element*au2THz
       end do
       write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
            &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
    end if

    ! ########################################################################

    deallocate(atom_pos,stat=ierr)
    call utils_dealloc_check('lr_phonons_born_asr','atom_pos',ierr)

    deallocate(dynamical_matrix, stat=ierr)
    call utils_dealloc_check('lr_phonons_dynmat', 'dynamical_matrix', ierr)
    deallocate(force_constant_matrix_temp, stat=ierr)
    call utils_dealloc_check('lr_phonons_dynmat', &
         'force_constant_matrix_temp', ierr)
    deallocate(local_fcm, stat=ierr)
    call utils_dealloc_check('lr_phonons_dynmat', 'local_fcm', ierr)
    deallocate(local_born, stat=ierr)
    call utils_dealloc_check('lr_phonons_dynmat', 'local_born', ierr)
    deallocate(local_x, stat=ierr)
    call utils_dealloc_check('lr_phonons_dynmat', 'local_x', ierr)
    deallocate(local_w, stat=ierr)
    call utils_dealloc_check('lr_phonons_dynmat', 'local_w', ierr)
    deallocate(local_u, stat=ierr)
    call utils_dealloc_check('lr_phonons_dynmat', 'local_u', ierr)
    deallocate(local_v, stat=ierr)
    call utils_dealloc_check('lr_phonons_dynmat', 'local_v', ierr)
    deallocate(local_ind_v, stat=ierr)
    call utils_dealloc_check('lr_phonons_dynmat', 'local_ind_v', ierr)
    deallocate(local_u_less, stat=ierr)
    call utils_dealloc_check('lr_phonons_dynmat', 'local_u_less', ierr)

    call dense_destroy(dyn_matrix_dem)
    call dense_destroy(unit_mat_dem)
    call dense_destroy(eigenvecs_dem)


    call timer_clock('lr_phonons_dynmat', 2)

  contains

     subroutine scal_prod2(u,v,ind_v,scal)

        implicit none

        ! arguments
        real(kind=DP), intent(in) :: u(3,3,mdl%nat,mdl%nat)
        integer, intent(in) :: ind_v(2,4)
        real(kind=DP), intent(in) :: v(2)
        real(kind=DP), intent(inout) :: scal

        ! local
        integer :: i

        scal = 0.0_DP
        do i=1,2
          scal = scal + u(ind_v(i,1),ind_v(i,2),ind_v(i,3),ind_v(i,4))*v(i)
        end do

        return

     end subroutine scal_prod2

  end subroutine lr_phonons_dynmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_phonons_get_basis_block(in_mat, out_mat, in_basis_left, &
      in_basis_right, out_basis_left, out_basis_right, order_left, &
      order_right)

    ! ========================================================================!
    ! This subroutine retrieves a block corresponding to one basis set from a !
    ! matrix represented by a larger set that contains the smaller basis      !
    ! ------------------------------------------------------------------------!
    ! The 'order-left' and 'order-right' parameters are used to obtain the    !
    ! position of the smaller set inside the larger set. For instance, if the !
    ! joint basis is formed out of the response+valence set, the order of the !
    ! response set would be 1, while the order of the valence one would be 2  !
    ! ------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in Dec 2016                           !
    ! Modified to remove pub_par by Joseph Prentice, October 2018             !
    ! ========================================================================!

    use comms, only: comms_barrier, pub_on_root, pub_my_proc_id, &
         pub_total_num_procs
    use constants, only: DP, stdout, cmplx_0
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: PARAL_INFO
    use sparse, only: SPAM3, sparse_scale, sparse_get_element, &
         sparse_put_element, sparse_element_exists, sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: in_mat
    type(SPAM3), intent(inout) :: out_mat
    type(FUNC_BASIS), intent(in) :: in_basis_left, in_basis_right, &
         out_basis_left, out_basis_right
    integer, intent(in) :: order_left, order_right

    ! Local variables
    integer :: ingwf, iat, loc_ingwf, ingwf_eff, jngwf, jat, loc_jngwf, &
         jngwf_eff, proc, ratio_left, ratio_right, loc_iat, loc_jat
    real(kind=DP) :: element
    complex(kind=DP) :: element_cmplx
    logical :: loc_cmplx

    type(PARAL_INFO), pointer :: row_par, col_par

    call timer_clock('lr_phonons_get_basis_block', 1)

    if (out_mat%iscmplx) then
       call sparse_scale(out_mat, cmplx_0)
    else
       call sparse_scale(out_mat, 0.0_DP)
    end if

    call comms_barrier

    if (out_mat%iscmplx.neqv.in_mat%iscmplx) then
       call utils_abort('Both matrices must be either real or complex')
    end if

    loc_cmplx = in_mat%iscmplx

    call sparse_get_par(row_par,in_mat,'R')
    call sparse_get_par(col_par,in_mat,'C')

    ! out assumption is that the any of the number of NGWFs in
    ! out_basis_left/right is either the same or half of the number of NGWFs
    ! in in_basis_left/right

    ratio_left = in_basis_left%num / out_basis_left%num
    ratio_right = in_basis_right%num / out_basis_right%num

    if ((ratio_right * out_basis_right%num) .ne. in_basis_right%num) then
       call utils_abort('Right in-basis not a perfect multiple of right &
            &out-basis')
    end if

    if ((ratio_left * out_basis_left%num) .ne. in_basis_left%num) then
       call utils_abort('Left in-basis not a perfect multiple of left &
            &out-basis')
    end if

    if ((ratio_left.ne.1).and.(ratio_left.ne.2)) then
       call utils_abort('Not one of the accepted left basis cases')
    end if

    if ((ratio_right.ne.1).and.(ratio_right.ne.2)) then
       call utils_abort('Not one of the accepted right basis cases')
    end if


    if ((order_left < 1) .or. (order_left > ratio_left)) then
       call utils_abort('Specified left basis order is out of bounds')
    end if

    if ((order_right < 1) .or. (order_right > ratio_right)) then
       call utils_abort('Specified right basis order is out of bounds')
    end if


    do proc = 0, pub_total_num_procs - 1
       do loc_iat = 1, row_par%num_atoms_on_proc(proc)
          iat = row_par%first_atom_on_proc(proc) + loc_iat - 1
          do loc_ingwf = 1, in_basis_left%num_on_atom(iat)
             ingwf = in_basis_left%first_on_atom(iat) + loc_ingwf - 1

             do loc_jat = 1, col_par%num_atoms_on_proc(pub_my_proc_id)
                jat = col_par%first_atom_on_proc(pub_my_proc_id) + loc_jat - 1
                do loc_jngwf = 1, in_basis_right%num_on_atom(jat)
                   jngwf = in_basis_right%first_on_atom(jat) + loc_jngwf - 1

                   if (.not.sparse_element_exists(in_mat,ingwf,jngwf)) cycle

                   if ((order_left.eq.1).and.(order_right.eq.1)) then

                      if ((loc_ingwf.le.out_basis_left%num_on_atom(iat)).and.&
                          (loc_jngwf.le.out_basis_right%num_on_atom(jat))) then

                         ingwf_eff = out_basis_left%first_on_atom(iat) + &
                               loc_ingwf - 1
                         jngwf_eff = out_basis_right%first_on_atom(jat) + &
                               loc_jngwf - 1

                         if (sparse_element_exists(out_mat, &
                              ingwf_eff, jngwf_eff)) then

                            if (.not. loc_cmplx) then
                               call sparse_get_element(element, in_mat, ingwf, &
                                    jngwf)
                               call sparse_put_element(element, out_mat, &
                                    ingwf_eff, jngwf_eff)
                            else
                               call sparse_get_element(element_cmplx, in_mat, &
                                    ingwf, jngwf)
                               call sparse_put_element(element_cmplx, out_mat, &
                                    ingwf_eff, jngwf_eff)
                            end if

                         end if
                      end if

                   else if ((order_left.eq.2).and.(order_right.eq.1)) then

                      if ((loc_ingwf.gt.out_basis_left%num_on_atom(iat)).and.&
                          (loc_jngwf.le.out_basis_right%num_on_atom(jat))) then

                         ingwf_eff = out_basis_left%first_on_atom(iat) + &
                               loc_ingwf - 1 - out_basis_left%num_on_atom(iat)
                         jngwf_eff = out_basis_right%first_on_atom(jat) + &
                               loc_jngwf - 1

                         if (sparse_element_exists(out_mat, &
                              ingwf_eff, jngwf_eff)) then

                            if (.not. loc_cmplx) then
                               call sparse_get_element(element, in_mat, ingwf, &
                                    jngwf)
                               call sparse_put_element(element, out_mat, &
                                    ingwf_eff, jngwf_eff)
                            else
                               call sparse_get_element(element_cmplx, in_mat, &
                                    ingwf, jngwf)
                               call sparse_put_element(element_cmplx, out_mat, &
                                    ingwf_eff, jngwf_eff)
                            end if

                         end if
                      end if

                   else if ((order_left.eq.1).and.(order_right.eq.2)) then

                      if ((loc_ingwf.le.out_basis_left%num_on_atom(iat)).and.&
                          (loc_jngwf.gt.out_basis_right%num_on_atom(jat))) then

                         ingwf_eff = out_basis_left%first_on_atom(iat) + &
                               loc_ingwf - 1
                         jngwf_eff = out_basis_right%first_on_atom(jat) + &
                               loc_jngwf - 1 - out_basis_right%num_on_atom(jat)

                         if (sparse_element_exists(out_mat, &
                              ingwf_eff, jngwf_eff)) then

                            if (.not. loc_cmplx) then
                               call sparse_get_element(element, in_mat, ingwf, &
                                    jngwf)
                               call sparse_put_element(element, out_mat, &
                                    ingwf_eff, jngwf_eff)
                            else
                               call sparse_get_element(element_cmplx, in_mat, &
                                    ingwf, jngwf)
                               call sparse_put_element(element_cmplx, out_mat, &
                                    ingwf_eff, jngwf_eff)
                            end if

                         end if
                      end if

                   else if ((order_left.eq.2).and.(order_right.eq.2)) then

                      if ((loc_ingwf.gt.out_basis_left%num_on_atom(iat)).and.&
                          (loc_jngwf.gt.out_basis_right%num_on_atom(jat))) then

                         ingwf_eff = out_basis_left%first_on_atom(iat) + &
                               loc_ingwf - 1 - out_basis_left%num_on_atom(iat)
                         jngwf_eff = out_basis_right%first_on_atom(jat) + &
                               loc_jngwf - 1 - out_basis_right%num_on_atom(jat)

                         if (sparse_element_exists(out_mat, &
                              ingwf_eff, jngwf_eff)) then

                            if (.not. loc_cmplx) then
                               call sparse_get_element(element, in_mat, ingwf, &
                                    jngwf)
                               call sparse_put_element(element, out_mat, &
                                    ingwf_eff, jngwf_eff)
                            else
                               call sparse_get_element(element_cmplx, in_mat, &
                                    ingwf, jngwf)
                               call sparse_put_element(element_cmplx, out_mat, &
                                    ingwf_eff, jngwf_eff)
                            end if

                         end if
                      end if

                   end if


                end do
             end do
          end do
       end do
    end do


    call comms_barrier

    call timer_clock('lr_phonons_get_basis_block', 2)

  end subroutine lr_phonons_get_basis_block

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calculate_dipole_mat(dipole_mat, mdl, bra_rep, bra_basis, &
       ket_rep, ket_basis, braket_overlap, proj_basis, nl_projectors, &
       first_order, braiGp_overlap, ketiGp_overlap, FO_overlap_braket, &
       direction, weight)

    ! ========================================================================!
    ! This subroutine calculates the dipole matrix (including corrections for !
    ! periodic systems), used in the computation of Born effective charges and!
    ! in the first-order NGWF gradient terms for E-field perturbations        !
    ! ------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in Dec 2016/ Jan 2017                 !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,        !
    ! October 2018                                                            !
    ! ========================================================================!

    use datatypes, only: FUNCTIONS
    use comms, only: comms_barrier, pub_on_root, pub_my_proc_id
    use constants, only: DP, stdout, cmplx_1, cmplx_0, cmplx_i
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_brappd_ketppd
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use optics, only: optics_pos_mat_els
    use paw, only: paw_projector_overlap
    use projectors, only: PROJECTOR_SET, projectors_func_ovlp_box
    use rundat, only: pub_num_spins, pub_debug_on_root, pub_spin_fac, pub_paw, &
         pub_any_nl_proj
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_scale, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_axpy, sparse_embed_transpose, &
         sparse_embed_transpose_structure, sparse_embed_product, sparse_embed_trace, &
         sparse_embed_copy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: dipole_mat(3)
    type(MODEL), intent(in) :: mdl
    type(NGWF_REP), intent(in) :: bra_rep
    type(FUNC_BASIS), intent(in) :: bra_basis(1)
    type(NGWF_REP), intent(in) :: ket_rep
    type(FUNC_BASIS), intent(in) :: ket_basis(1)
    type(SPAM3_EMBED), intent(in) :: braket_overlap
    type(FUNC_BASIS), intent(in) :: proj_basis(1)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(1)
    logical, intent(in) :: first_order
    type(SPAM3_EMBED), intent(in), optional :: braiGp_overlap
    type(SPAM3_EMBED), intent(in), optional :: ketiGp_overlap
    type(SPAM3_EMBED), intent(in), optional :: FO_overlap_braket
    integer, intent(in), optional :: direction(mdl%nat)
    real(kind=DP), intent(in), optional :: weight(mdl%nat)

    ! Local Variables

    type(SPAM3_EMBED)  :: fact_braiGp_ovlp(3), fact_brap_ovlp(3), fact_ovlp_braket(3)
    type(SPAM3_EMBED)  :: O_mat, pket_overlap, iGpket_overlap
    integer      :: iter, xyz, ierr, atom
    type(SPAM3_EMBED)  :: temp_mat_1, temp_mat_2
    real(kind=DP)  :: fscale, coc(3), factor

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &calculate_dipole_mat'

    ! Start timer
    call timer_clock('calculate_dipole_mat',1)

    ! pdh: sync procs
    call comms_barrier

    ! check input integrity
    if (first_order) then
       if (.not. (present(braiGp_overlap).and.present(ketiGp_overlap).and. &
            present(FO_overlap_braket).and.present(direction).and. &
            present(weight)) ) &
            call utils_abort('Incomplete first-order dipole matrix')
    end if


    if (pub_paw) then
       call sparse_embed_transpose_structure(pket_overlap%structure, &
            ket_rep%sp_overlap)
       call sparse_embed_create(pket_overlap)
       call sparse_embed_transpose(pket_overlap, ket_rep%sp_overlap)

       if (first_order) then
          call sparse_embed_transpose_structure(iGpket_overlap%structure, &
               ketiGp_overlap)
          call sparse_embed_create(iGpket_overlap)
          call sparse_embed_transpose(iGpket_overlap, ketiGp_overlap)

          do iter = 1,3
             call sparse_embed_create(fact_braiGp_ovlp(iter), bra_rep%sp_overlap)
          end do
       end if

       do iter = 1,3
          call sparse_embed_create(fact_brap_ovlp(iter), bra_rep%sp_overlap)
       end do

       O_mat%structure = 'E'
       call sparse_embed_create(O_mat)
       call paw_projector_overlap(O_mat%p, mdl%regions(1)%paw_sp)
    end if

    do iter = 1,3
       call sparse_embed_create(fact_ovlp_braket(iter), braket_overlap)

       call sparse_embed_scale(dipole_mat(iter), 0.0_DP)
    end do

    if (.not.first_order) then
       call optics_pos_mat_els(dipole_mat, bra_rep, bra_basis, ket_rep, &
            ket_basis, braket_overlap, braket_overlap, proj_basis, mdl)
    else
       if (pub_paw) then ! this term would be 0 outside the PAW/USP framework
          call optics_pos_mat_els(dipole_mat, bra_rep, bra_basis, ket_rep, &
               ket_basis, braket_overlap, braket_overlap, proj_basis, mdl, .true., &
               braiGp_overlap, ketiGp_overlap, direction, weight)
       end if
    end if

    ! calculate center of ionic charge
    fscale = 1.0_DP/sum(mdl%elements(:)%ion_charge)
    coc(1) =  sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%X)*fscale
    coc(2) =  sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%Y)*fscale
    coc(3) =  sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%Z)*fscale

    do xyz = 1, 3
       if (.not. first_order) then
          call function_ops_brappd_ketppd(fact_ovlp_braket(xyz)%p, &
               bra_rep%ngwfs_on_grid(1), bra_basis(1), ket_rep%ngwfs_on_grid(1), &
               ket_basis(1), mdl%cell, factor_type = xyz)
       end if

       if (pub_paw) then
          call projectors_func_ovlp_box(fact_brap_ovlp(xyz)%p, &
               bra_rep%ngwfs_on_grid(1), bra_basis(1), proj_basis(1), nl_projectors(1), &
               mdl%fftbox, mdl%cell, factor_type = xyz)

          if (first_order) then
             call grad_proj_ovlp(mdl, fact_braiGp_ovlp(xyz)%p, &
                  bra_rep%ngwfs_on_grid(1), bra_basis(1), proj_basis(1), nl_projectors(1), &
                  direction, weight, factor_type = xyz)
          end if

          if (.not. first_order) then
             call sparse_embed_create(temp_mat_1, fact_brap_ovlp(xyz), O_mat)
             call sparse_embed_product(temp_mat_1, fact_brap_ovlp(xyz), O_mat)
             call sparse_embed_create(temp_mat_2, temp_mat_1, pket_overlap)
             call sparse_embed_product(temp_mat_2, temp_mat_1, pket_overlap)
             call sparse_embed_axpy(fact_ovlp_braket(xyz), temp_mat_2, 1.0_DP)
             call sparse_embed_destroy(temp_mat_1)
             call sparse_embed_destroy(temp_mat_2)
          else
             call sparse_embed_create(temp_mat_1, fact_brap_ovlp(xyz), O_mat)
             call sparse_embed_product(temp_mat_1, fact_brap_ovlp(xyz), O_mat)
             call sparse_embed_create(temp_mat_2, temp_mat_1, iGpket_overlap)
             call sparse_embed_product(temp_mat_2, temp_mat_1, iGpket_overlap)
             call sparse_embed_axpy(fact_ovlp_braket(xyz), temp_mat_2, 1.0_DP)
             call sparse_embed_destroy(temp_mat_1)
             call sparse_embed_destroy(temp_mat_2)

             call sparse_embed_create(temp_mat_1, fact_braiGp_ovlp(xyz), O_mat)
             call sparse_embed_product(temp_mat_1, fact_braiGp_ovlp(xyz), O_mat)
             call sparse_embed_create(temp_mat_2, temp_mat_1, pket_overlap)
             call sparse_embed_product(temp_mat_2, temp_mat_1, pket_overlap)
             call sparse_embed_axpy(fact_ovlp_braket(xyz), temp_mat_2, 1.0_DP)
             call sparse_embed_destroy(temp_mat_1)
             call sparse_embed_destroy(temp_mat_2)
          end if

       end if

       if (first_order) then
          factor = 0.0_DP
          do iter = 1, mdl%nat
             if (direction(iter) == xyz) then
                factor = factor + mdl%elements(iter)%ion_charge * &
                     weight(iter) * fscale
             end if
          end do

          call sparse_embed_axpy(dipole_mat(xyz), FO_overlap_braket, &
               -1.0_DP * coc(xyz))
          call sparse_embed_axpy(dipole_mat(xyz), braket_overlap, -1.0_DP * factor)

          call sparse_embed_axpy(dipole_mat(xyz), fact_ovlp_braket(xyz), -1.0_DP)
       else
          call sparse_embed_axpy(dipole_mat(xyz), braket_overlap, -1.0_DP * coc(xyz))
          call sparse_embed_axpy(dipole_mat(xyz), fact_ovlp_braket(xyz), -1.0_DP)
       end if

    end do

    call comms_barrier

    ! destroy matrices
    if (pub_paw) then
       call sparse_embed_destroy(pket_overlap)

       if (first_order) then
          call sparse_embed_destroy(iGpket_overlap)
          do iter = 1,3
             call sparse_embed_destroy(fact_braiGp_ovlp(iter))
          end do
       end if

       do iter = 1,3
          call sparse_embed_destroy(fact_brap_ovlp(iter))
       end do

       call sparse_embed_destroy(O_mat)
    end if

    do iter = 1,3
       call sparse_embed_destroy(fact_ovlp_braket(iter))
    end do

    ! pdh: sync procs
    call comms_barrier

    ! Stop timer
    call timer_clock('calculate_dipole_mat',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &calculate_dipole_mat'

  end subroutine calculate_dipole_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_phonons_transform(transtype, directions, weights, nat, &
       at_gamma, Born_charges, force_constant_matrix)

    ! ========================================================================!
    ! This subroutine is used to transform to/from the 1-atom 1-direction     !
    ! perturbation-basis from/to the all-atoms 1-direction perturbation-basis !
    ! ------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in Jan 2017                           !
    ! Modified to remove pub_par by Joseph Prentice, October 2018             !
    ! ========================================================================!

    use comms, only: comms_barrier
    use constants, only: stdout, PI, DP, cmplx_0, cmplx_1
    use dense, only: DEM,dense_create,dense_convert,dense_destroy,dense_scale, &
         dense_product,dense_get_element,dense_put_element,dense_transpose, &
         dense_copy, dense_axpy
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    ! Arguments
    character, intent(in) :: transtype
    integer, intent(in) :: nat
    integer, intent(inout) :: directions(3*nat+3, nat)
    real(kind=DP), intent(inout) :: weights(3*nat+3, nat)
    logical, intent(in) , optional:: at_gamma
    real(kind=DP), intent(inout), optional :: Born_charges(3,3*nat)
    complex(kind=DP), intent(inout), optional :: force_constant_matrix(3* &
         nat+3,nat)

    ! Local Variables
    type(DEM) :: den_mat_pp_1, den_mat_pp_2, transfo_mat, transfo_mat_TR, &
         den_mat_pp_cmplx_1, den_mat_pp_cmplx_2, transfo_mat_cmplx, &
         transfo_mat_cmplx_TR

    integer :: iter1, iter2, iter3, iatom, iter1_eff, iter2_eff, ierr
    real(kind=DP) :: element

    real(kind=DP), allocatable :: Born_charges_temp(:,:)

    call dense_create(den_mat_pp_1, 3*nat, 3*nat)
    call dense_create(den_mat_pp_2, 3*nat, 3*nat)
    call dense_create(transfo_mat, 3*nat, 3*nat)
    call dense_create(transfo_mat_TR, 3*nat, 3*nat)

    call dense_create(den_mat_pp_cmplx_1, 3*nat, 3*nat, iscmplx=.true.)
    call dense_create(den_mat_pp_cmplx_2, 3*nat, 3*nat, iscmplx=.true.)
    call dense_create(transfo_mat_cmplx, 3*nat, 3*nat, iscmplx=.true.)
    call dense_create(transfo_mat_cmplx_TR, 3*nat, 3*nat,iscmplx=.true.)

    if (transtype == 'B') then
       if ((.not.present(at_gamma)).or.(.not.present(Born_charges)).or. &
            (.not.present(force_constant_matrix))) then
          call utils_abort('Insufficient variables for backward transformation')
       end if

    end if

    if (transtype == 'B') then
       allocate(Born_charges_temp(3,3*nat), stat=ierr)
       call utils_alloc_check('lr_phonons_transform', 'Born_charges_temp', ierr)
    end if


    call comms_barrier

    ! initialise transformation mat as identity matrix
    call dense_scale(transfo_mat, 0.0_DP)
    do iter1 = 1, 3*nat
       do iter2 = 1, 3*nat
          if (iter1 == iter2) then
             call dense_put_element(1.0_DP, transfo_mat, iter1, iter2)
          end if
       end do
    end do

    call comms_barrier

    ! compute rotation matrices and multiply them one at a time
    do iter1 = 1, 3*(nat-1)
       do iter2 = iter1+3, 3*nat, 3

          call dense_scale(den_mat_pp_1, 0.0_DP)

          element = PI / real(nat+1,kind=DP) ! angle of multidim rotation
          ! note that element must be < PI / nat, condition which is
          ! satisfied
          call dense_put_element(cos(element), den_mat_pp_1, iter1, iter1)
          call dense_put_element(cos(element), den_mat_pp_1, iter2, iter2)
          call dense_put_element(sin(element), den_mat_pp_1, iter2, iter1)
          call dense_put_element(-sin(element), den_mat_pp_1, iter1, iter2)

          do iter1_eff = 1,3*nat
             if ( (iter1_eff.ne.iter1) .and. (iter1_eff.ne.iter2) ) then
                call dense_put_element(1.0_DP, den_mat_pp_1, iter1_eff, &
                     iter1_eff)
             end if
          end do

          call comms_barrier

          call dense_product(den_mat_pp_2, transfo_mat, den_mat_pp_1)
          call dense_copy(transfo_mat, den_mat_pp_2)

          call comms_barrier

       end do
    end do

    call dense_transpose(transfo_mat_TR, transfo_mat)

    ! FORWARD TRANSFORMATION
    if (transtype == 'F') then
       ! in the end, the transformation matrix is transfo_mat
       ! now obtain the weights and directions
       directions = 0
       weights = 0.0_DP

       do iter1 = 1, 3*nat
          do iter2 = 1, 3*nat
             iatom = ((iter2 - 1) / 3) + 1
             call dense_get_element(element, transfo_mat, iter1, iter2)

             if (abs(element) > 1.0d-9) then
                if (directions(iter1, iatom).ne.0) then
                   call utils_abort('Multiple directions in a single atom')
                else
                   directions(iter1, iatom) = 1 + mod((iter2 - 1), 3)
                   weights(iter1, iatom) = element
                end if
             end if

          end do
       end do

       call comms_barrier

       ! check integrity of directions and weights
       do iter1 = 1,3*nat
          do iter2 = 1, nat
             if (directions(iter1,iter2) == 0) then
                call utils_abort('MULTIDIMENSIONAL ROTATION FAILED')
             end if
          end do
       end do


       do iter1 = 1,3*nat
          do iter2 = 1, nat
             if (abs(weights(iter1,iter2)) < 0.000000001_DP) then
                call utils_abort('MULTIDIMENSIONAL ROTATION FAILED')
             end if
          end do
       end do

       ! gcc32 TEST
       !do iter1 = 1,3*nat
       !   if (pub_on_root) then
       !      element = 0.0_DP
       !      do iter2 = 1,nat
       !         element = element + weights(iter1,iter2)*weights(iter1,iter2)
       !      end do
       !      write(stdout,*) 'PERT ', iter1, ' NORM: ', sqrt(element)
       !      write(stdout,*) 'PERT ', iter1, ' DIR: ', directions(iter1,:)
       !      write(stdout,*) 'PERT ', iter1, ' WEIGHT: ', weights(iter1,:)
       !   end if
       !end do
       ! gcc32 END TEST


    else ! BACKWARD TRANSFORMATION

       call dense_scale(transfo_mat_cmplx, cmplx_0)
       call dense_axpy(transfo_mat_cmplx, transfo_mat, cmplx_1)
       call dense_scale(transfo_mat_cmplx_TR, cmplx_0)
       call dense_axpy(transfo_mat_cmplx_TR, transfo_mat_TR, cmplx_1)


       ! TRANSFORM FORCE_CONSTANT MATRIX INTO THE ORIGINAL BASIS OF 1-ATOM 1-DiR
       ! PERTURBATIONS
       ! =======================================================================

       call dense_scale(den_mat_pp_cmplx_1, 0.0_DP)
       do iter1 = 1,3*nat
          do iter2 = 1,3*nat
             call dense_put_element(force_constant_matrix(iter1,iter2), &
                  den_mat_pp_cmplx_1, iter1, iter2)
          end do
       end do

       call dense_product(den_mat_pp_cmplx_2, den_mat_pp_cmplx_1, &
            transfo_mat_cmplx)
       call dense_product(den_mat_pp_cmplx_1, transfo_mat_cmplx_TR, &
            den_mat_pp_cmplx_2)

       force_constant_matrix(1:3*nat,1:3*nat) = cmplx_0

       do iter1 = 1,3*nat
          do iter2 = 1,3*nat
             call dense_get_element(force_constant_matrix(iter1,iter2), &
                  den_mat_pp_cmplx_1, iter1, iter2)
             ! if (pub_on_root) write(stdout,*) 'PERTURBATION PAIR ', iter1, &
             !       iter2, ' : ', force_constant_matrix(iter1,iter2)
          end do
       end do

       ! If at Gamma, transform the Born charges back to the original
       ! perturbation basis
       if (at_gamma) then

           Born_charges_temp = 0.0_DP
           do iter1 = 1,3*nat
              do iter2 = 1,3*nat
                 call dense_get_element(element, transfo_mat, iter2, iter1)
                 do iter3 = 1,3
                    Born_charges_temp(iter3, iter1) = &
                         Born_charges_temp(iter3, iter1) + &
                         Born_charges(iter3, iter2) * element
                 end do
              end do
           end do

           ! copy temporary matrix to initial Born-charges array
           Born_charges(:,:) = Born_charges_temp(:,:)
       end if

    end if

    call dense_destroy(transfo_mat)
    call dense_destroy(transfo_mat_TR)
    call dense_destroy(den_mat_pp_1)
    call dense_destroy(den_mat_pp_2)

    call dense_destroy(transfo_mat_cmplx)
    call dense_destroy(transfo_mat_cmplx_TR)
    call dense_destroy(den_mat_pp_cmplx_1)
    call dense_destroy(den_mat_pp_cmplx_2)

    if (transtype == 'B') then
       deallocate(Born_charges_temp, stat=ierr)
       call utils_dealloc_check('lr_phonons_transform', 'Born_charges_temp', &
            ierr)
    end if


  end subroutine lr_phonons_transform

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module lr_phonons

