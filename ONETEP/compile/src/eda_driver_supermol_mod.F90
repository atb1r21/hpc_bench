! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                                                                             !
!                   E N E R G Y    D E C O M P O S I T I O N                  !
!                   S U P E R M O L E C U L E    D R I V E R                  !
!                                                                             !
!=============================================================================!
!                                                                             !
! This module handles the fragment and supermolecule density data methods,    !
! and the supermolecule EDA driver routine (polarisation, charge transfer).   !
!                                                                             !
!-----------------------------------------------------------------------------!
! Originally written by Max Phipps August 2014                                !
!-----------------------------------------------------------------------------!

module eda_driver_supermol

  use constants, only: dp, stdout
  use model_type, only: MODEL

  implicit none

  private

  ! subroutines
  ! main drivers
  public :: eda_driver_supermol_run
  public :: eda_driver_prep_supermol

  ! delta densities construction and storage routines:
  public :: eda_allocate_densities
  public :: eda_store_density
  public :: eda_destroy_frz_density
  public :: eda_destroy_pol_density
  public :: eda_destroy_poliso_density
  public :: eda_destroy_full_density

  ! supermolecule restart IO routines:
  public :: eda_super_restart_write_metadata ! dkn, NGWFs, metadata
  !public :: eda_super_restart_read_dkn_ngwfs ! dkn, NGWFs
  public :: eda_super_restart_read_metadata  ! metadata

  integer :: eda_cont_mode = 0 ! continution file pub_eda_mode
  integer :: eda_cont_fc1 = 0  ! continution file pub_frag_counter
  integer :: eda_cont_fc2 = 0  ! continution file pub_frag_counter2
  integer :: eda_cont_ctdeloc_it = 1 ! continuation file delocalisation
                                     ! calculation iterator

contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine eda_driver_supermol_run(mdl, hfxstate, ngwf_basis, aux_ngwf_basis,&
       cond_ngwf_basis, joint_ngwf_basis, proj_basis, hub_proj_basis, hub, &
       core_basis, core_wvfns, rep, ham, denskern, edft, nl_projectors, &
       ngwf_nonsc_forces, lhxc_fine)

    !======================================================================!
    ! This suboutine is the driver for the supermolecule stage of the EDA  !
    ! calculation.  This driver handles the fragment-specific and/or       !
    ! full-supermolecule-only polarisation calculation, the charge         !
    ! transfer calculation, and optional delta density calculations.       !
    ! TODO: Implement SCF-MI with DIIS optimisation.                       !
    ! TODO: Spin polarised cases.                                          !
    !----------------------------------------------------------------------!
    ! Written by Max Phipps July 2015.                                     !
    ! Modified for embedding by Joseph Prentice, September 2018            !
    !======================================================================!

    use comms, only: pub_on_root
    use constants, only: EDA_FROZENIDEM, EDA_POLSIMUL, &
         EDA_POLFRAGLOC_DEVEL, EDA_CTFRAGLOC_DEVEL, &
         EDA_POLSIMUL_OFF, EDA_POLFRAGLOC_DEVEL_OFF, &
         EDA_FULL, paw_en_size, ANGSTROM
    use electronic_init, only: electronic_init_denskern
    use ensemble_dft_type, only: EDFT_MODEL
    use fragment_data, only: pub_frag_data, &
         eda_frzIdem_density_fine, eda_pol_density_fine, &
         eda_poliso_density_fine, eda_full_density_fine, &
         complex_energy_frz, complex_energy_frzIdem, &
         complex_energy_pol_simul, complex_energy_pol, &
         complex_energy_ct, complex_energy_full, &
         eda_dfdtau_fine, &
         fragment_data_alloc_paw_core_densities
    use fragment_matrix_ops, only: fmo_construct_masks, &
         fmo_load_internal_blkdiag_denskern
    use fragment_scfmi, only: scfmi_construct_nonorth_kernel, denskern_R
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use hamiltonian, only: hamiltonian_dens_indep_matrices, &
         hamiltonian_dens_dep_matrices, hamiltonian_energy_components, &
         hamiltonian_build_matrix
    use kernel, only: DKERN, kernel_workspace_invalidate
    use model_type, only: MODEL
    use ngwf_cg, only: ngwf_cg_optimise
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use ngwfs, only: ngwfs_initialise, ngwfs_initialise_from_fragments
    use projectors, only: PROJECTOR_SET
    use properties, only: properties_plot_delta_density, &
         properties_calculate
    use restart, only: restart_kernel_write, restart_ngwfs_tightbox_output
    use rundat, only: pub_eda_isol, pub_eda_isol_x, pub_eda_isol_c, &
         pub_eda_intrafrag_x_frz, pub_eda_intrafrag_c_frz, &
         pub_eda_frz_x, pub_eda_rep_x, &
         pub_eda_frz_c, pub_eda_rep_c, &
         pub_eda_dE_x_rep, pub_eda_dE_c_rep, &
         pub_frag_deloc, &
         pub_eda_preptool, pub_eda_write, pub_eda_read_super, &
         pub_eda_continuation, &
         pub_eda_scfmi, pub_eda_scfmi_any, &
         pub_eda_frag_isol_pol, pub_eda_reset_ngwfs_pol, &
         pub_eda_frag_isol_ct, pub_eda_reset_ngwfs_ct, &
         pub_frag_iatm, pub_eda_mode, pub_frag_counter, pub_frag_counter2, &
         pub_eda_deltadens, &
         pub_num_spins, pub_nonsc_forces, &
         pub_maxit_ngwf_cg, &
         pub_write_density_plot, pub_use_aux_ngwfs, &
         pub_write_ngwf_plot, pub_cube_format, pub_grd_format, &
         pub_dx_format, pub_do_properties, pub_kerfix, &
         pub_exact_lnv, PUB_1K
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_copy, sparse_embed_array_scale, sparse_embed_array_copy
    use utils, only: utils_banner
    use visual, only: visual_ngwfs, visual_scalarfield

    ! Arguments
    type(MODEL), intent(inout)         :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    type(FUNC_BASIS), intent(inout)    :: ngwf_basis(1)
    type(FUNC_BASIS), intent(inout)    :: aux_ngwf_basis(1)
    type(FUNC_BASIS), intent(inout)    :: cond_ngwf_basis(1)
    type(FUNC_BASIS), intent(inout)    :: joint_ngwf_basis(1)
    type(FUNC_BASIS), intent(inout)    :: proj_basis(1)
    type(FUNC_BASIS), intent(inout)    :: hub_proj_basis(1)
    type(HUBBARD_MODEL)                :: hub
    type(FUNC_BASIS), intent(inout)    :: core_basis(1)
    type(NGWF_REP), intent(inout)      :: rep
    type(NGWF_HAM), intent(inout)      :: ham
    type(DKERN), intent(inout)         :: denskern
    type(EDFT_MODEL), intent(inout)    :: edft  ! ensemble DFT container (dummy)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(1)
    type(PROJECTOR_SET), intent(inout) :: core_wvfns(1)
    real(kind=DP), dimension(:,:), intent(inout) :: ngwf_nonsc_forces
    real(kind=DP), dimension(:,:,:,:), intent(inout) :: lhxc_fine

    ! Local Variables
    real(kind=DP)     :: dummy_energy, dummy_energy2
    real(kind=DP)     :: dummy_paw_sphere_energies(paw_en_size)
    character(len=10) :: fnstr
    logical           :: exact_lnv_temp_off = .false.
    integer           :: pub_kerfix_temp
    logical           :: ngwfs_converged
    integer           :: it
    logical           :: out_of_runtime
    logical           :: eda_write_density_plot
    logical           :: ngwfs_changed  ! to avoid unneccesary hamiltonian routine calls
    integer           :: deloc_iterator_start
    type(SPAM3_EMBED) :: olap_backup, inv_olap_backup
    integer           :: ierr
    character(255)    :: errmsg = ""

    eda_dfdtau_fine = eda_get_dfdtau_fine( &
         ham, &
         mdl, &
         rep, &
         ngwf_basis, &
         denskern%kern)

    ! mjsp: Initialisations:
    ! mjsp: Initialise ngwfs_changed flag (used to avoid
    ! unneccesary calls to hamiltonian_dens_indep_matrices)
    ngwfs_changed = .false.

    ! mjsp: Setup write cube file logical
    eda_write_density_plot = (pub_write_density_plot .and. &
         .not.(pub_use_aux_ngwfs) .and. &
         (pub_cube_format .or. pub_grd_format .or. pub_dx_format))

    ! mjsp: Calculation involves SCF-MI at some point:
    pub_eda_scfmi_any = .true.

    ! initialise start index of iterator for the
    ! delocalisation calculation mode (EDA_CTFRAGLOC_DEVEL)
    if (pub_eda_continuation) then
       deloc_iterator_start = eda_cont_ctdeloc_it
    else
       deloc_iterator_start = 1
    end if

    ! mjsp: If electron density and/or electron density difference
    ! mjsp: to be plotted, then allocate intermediate density storage
    if (pub_eda_deltadens .or. eda_write_density_plot) &
       call eda_allocate_densities(mdl)


    ! mjsp: If we are reading in the supermolecule data then we can skip
    ! the frozen calculation, as the '.eda' file includes this data.
    ! otherwise...
    if (.not.(pub_eda_read_super)) then

       ! mjsp: (1) non-variational calculation of frozen density interaction
       ! energy

       ! mjsp: Calculate the energy using the frozen density kernel and store in
       ! complex_energy_frz:
       call hamiltonian_dens_dep_matrices(ham, lhxc_fine, complex_energy_frz, &
            dummy_energy, dummy_energy2, dummy_paw_sphere_energies, rep, &
            ngwf_basis, hub_proj_basis, hub, denskern%kern, mdl, hfxstate, &
            .true., .false., dfdtau_fine = eda_dfdtau_fine)

       ! mjsp: Plot NGWFs if requested
       if (pub_write_ngwf_plot .and. &
            (pub_cube_format .or. pub_grd_format .or. pub_dx_format)) &
          call visual_ngwfs(rep%ngwfs_on_grid(1), ngwf_basis(1), &
                  'frz_', mdl%regions(1)%elements, mdl%cell, &
                  mdl%fftbox, mdl%par)

       if(pub_on_root) write(stdout,'(/a)') "EDA: Printing energy components for &
            &pre-orthogonalised state..."

       ! Printout energy components:
       call hamiltonian_energy_components( &
            denskern%kern%m(:,PUB_1K), rep, mdl, &
            ngwf_basis, hub_proj_basis, hub, &
            ham%hfexchange)

       ! interfragment X,C frozen interaction calculation:
       ! accumulate corrections for fragment X,C energies (i.e. -interfragmental):
       do it=1,pub_frag_iatm(0)
          pub_eda_intrafrag_x_frz = pub_eda_intrafrag_x_frz - pub_eda_isol_x(it)
          pub_eda_intrafrag_c_frz = pub_eda_intrafrag_c_frz - pub_eda_isol_c(it)
       end do
       ! +interfragmental +intrafragmental (to leave just intrafragmental
       ! interaction)
       pub_eda_intrafrag_x_frz = pub_eda_intrafrag_x_frz + pub_eda_frz_x
       pub_eda_intrafrag_c_frz = pub_eda_intrafrag_c_frz + pub_eda_frz_c

    end if  ! not(pub_eda_read_super)

    ! mjsp: Check not the case that we have loaded in continuation data
    ! mjsp: and calculation stopped beyond this stage
    if (.not.(pub_eda_continuation) .or. &
         (pub_eda_continuation .and. (eda_cont_mode .lt. EDA_FROZENIDEM))) then

       ! mjsp: Write out the supermolecule metadata
       ! mjsp: for supermolecule calculation continuation.
       if (pub_eda_write) then
          call eda_super_restart_write_metadata(is_preptool=.false.)
       else if (pub_eda_preptool) then
          call eda_super_restart_write_metadata(is_preptool=.true.)
       end if

       ! mjsp: Write out the kernel and NGWFs
       if (pub_eda_write .or. pub_eda_preptool) then
          call restart_kernel_write(denskern%kern)
          call restart_ngwfs_tightbox_output(rep%ngwfs_on_grid(1), &
               ngwf_basis(1), mdl%cell, mdl%fftbox, mdl%regions(1)%elements, &
                'tightbox_'//trim(ngwf_basis(1)%name),mdl%regions(1))
       end if

    end if

    ! If this is a supermolecule preparation calculation then
    ! do not do the supermolecule analysis.
    ! otherwise...
    if (.not.(pub_eda_preptool)) then

       ! mjsp: Check not the case that we have loaded in continuation data
       ! mjsp: and calculation stopped at or beyond this stage
       ! mjsp: (.lt. EDA_FROZENIDEM to avoid loading the continuation
       ! mjsp: NGWFs/dkn instead of the EDA_PREP stage NGWFs/dkn)
       if (.not.(pub_eda_continuation) .or. &
            (pub_eda_continuation .and. (eda_cont_mode .lt. EDA_FROZENIDEM))) then

          ! mjsp: (2) Properly antisymmetrized frozen density interaction energy
          ! mjsp: This energy is stored in the public variable complex_energy_frzIdem
          pub_eda_mode = EDA_FROZENIDEM

          ! Effectively a single Stoll locally projected SCF pass:
          pub_eda_scfmi = .true.

          ! Rebuild the density-indep part of the Hamiltonian for rebuilding
          ! the antisymmetrized density kernel
          call hamiltonian_dens_indep_matrices(rep, ngwf_basis,proj_basis, &
               nl_projectors, hub_proj_basis, hub, mdl)! val_rep,val_ngwf_basis)


          ! Construct the kernel via denskern_R: proper representation of
          ! the density kernel in the full overlap
          call scfmi_construct_nonorth_kernel(denskern,rep)

          ! mjsp: Stored versions of ks and ksk no longer valid
          call kernel_workspace_invalidate(denskern)


          ! mjsp: Calculate the energy using the new frozen density kernel and
          ! store in complex_energy_frzIdem:
          call hamiltonian_dens_dep_matrices(ham, lhxc_fine, complex_energy_frzIdem, &
               dummy_energy, dummy_energy2, dummy_paw_sphere_energies, rep, &
               ngwf_basis, hub_proj_basis, hub, denskern_R%kern, mdl, hfxstate,&
               .true., .false., dfdtau_fine = eda_dfdtau_fine)

          if(pub_on_root) write(stdout,'(/a)') "EDA: Printing energy components for &
               &orthogonalised state..."

          ! mjsp: Printout energy components:
          call hamiltonian_energy_components( &
               denskern_R%kern%m(:,PUB_1K), rep, mdl, &
               ngwf_basis, hub_proj_basis, hub, &
               ham%hfexchange)


          ! interfragment X,C MO occupancy fixing (Pauli repulsion) interaction
          ! calculation:
          pub_eda_dE_x_rep = pub_eda_rep_x - pub_eda_frz_x
          pub_eda_dE_c_rep = pub_eda_rep_c - pub_eda_frz_c

          ! mjsp: If continuation write is true, then write data
          if (pub_eda_write) &
             call eda_super_restart_write_metadata(is_preptool=.false.)

          ! mjsp: If electron density and/or electron density difference
          ! mjsp: to be plotted
          if (pub_eda_deltadens .or. eda_write_density_plot) then

             if(pub_on_root) write(stdout,'(a)') &
                  "EDA Density: Calculating density on grid for frozen &
                  &idempotent state..."

             ! mjsp: construct and store the frozen state density
             ! on the grid for delta polarisation density calculation
             call eda_store_density(eda_frzIdem_density_fine, rep, mdl, &
                  ngwf_basis, denskern_R%kern%m(:,PUB_1K))

             if (eda_write_density_plot) then
                ! mjsp: output delta density in plot format file
                ! vm: output density in Angstrom rather than in Bohr
                ! jd: ... but leave the unit conversion to visual_scalarfield
                call visual_scalarfield( &
                     eda_frzIdem_density_fine(:,:,:,1), mdl%fine_grid, mdl%cell, &
                     'Electronic density (in e/ang^3) for:', &
                     '_frzidem_ed_density', &
                     mdl%elements, ANGSTROM**3)
                if (pub_num_spins == 2) call visual_scalarfield( &
                     eda_frzIdem_density_fine(:,:,:,2), mdl%fine_grid, mdl%cell, &
                     'Electronic spin density (in e/ang^3) for:', &
                     '_frzidem_ed_spindensity', &
                     mdl%elements, ANGSTROM**3)
             end if

             if(pub_on_root) write(stdout,'(a)') "EDA Density: ...done"

          end if

       end if ! continuation check


       ! ================ INFORMATION PRINTOUT ============
       if (pub_on_root) then
          write(stdout,'(//a)') utils_banner('=')
          write(stdout,'(a)') utils_banner('=', &
               'Energy Decomposition Analysis (EDA)')
          write(stdout,'(a/)') utils_banner('=')
          write(stdout,'(2x,a,f10.3,a)') &
               '(1) (Nonorthogonal) Frozen Interaction Energy  : ', &
               (complex_energy_frz - pub_eda_isol(0))*627.509438736, ' kcal/mol'
          write(stdout,'(15x,a,f10.3,a)') &
               'Frozen Electrostatics E_ES        : ', &
               (complex_energy_frz - pub_eda_isol(0) - pub_eda_intrafrag_x_frz &
               - pub_eda_intrafrag_c_frz)*627.509438736, ' kcal/mol'
          write(stdout,'(15x,a,f10.3,a)') &
               'Frozen Exchange E_EX              : ', &
               (pub_eda_intrafrag_x_frz) *627.509438736, ' kcal/mol'
          write(stdout,'(15x,a,f10.3,a/)') &
               'Frozen Correlation E_FRZ-CORR     : ', &
               (pub_eda_intrafrag_c_frz) *627.509438736, ' kcal/mol'

          write(stdout,'(2x,a,f10.3,a)') &
               '(2) Repulsion Interaction Energy               : ', &
               (complex_energy_frzIdem - complex_energy_frz)*627.509438736, &
               ' kcal/mol'
          write(stdout,'(15x,a,f10.3,a)') &
               'Pauli-Repulsion E_REP             : ', &
               (complex_energy_frzIdem - complex_energy_frz - pub_eda_dE_c_rep) &
               *627.509438736, ' kcal/mol'
          write(stdout,'(15x,a,f10.3,a/)') &
               'Repulsion Correlation E_REP-CORR  : ', &
               (pub_eda_dE_c_rep)*627.509438736, ' kcal/mol'

          write(stdout,'(2x,a,f10.3,a/)') &
               '(1+2) Frozen Density Component E_FRZ           : ', &
               (complex_energy_frzIdem - pub_eda_isol(0)) *627.509438736, &
               ' kcal/mol'
          write(stdout,'(a/)') utils_banner('=')
       end if
       ! ================ END INFORMATION PRINTOUT ============

       ! mjsp: (3) variational calculation of polarisation energy

       ! mjsp: exact LNV involves normalisation rescales using norm_fac.
       ! The code currently is able to handle non-exact LNV by using a
       ! modified fragment rescaling routine.  Exact LNV requires more
       ! extensive modifications to handle this.
       ! --> Therefore force exact LNV off:
       if (pub_exact_lnv) then
          if (pub_on_root) write(stdout,'(/a)') 'WARNING: Turning off exact &
               &LNV for fragment polarisation calculation(s).'
          exact_lnv_temp_off = .true.
          pub_exact_lnv = .false.
       end if

       ! mjsp: pub_kerfix must be set to 2 to prevent interfragment
       ! delocalisations (in kernel_mod:internal_correct_gradient())
       ! if we have any fragments with no occupied orbitals.
       ! loop the fragments:
       pub_kerfix_temp = pub_kerfix
       if (pub_kerfix .ne. 2) then
          if (pub_on_root) write(stdout,'(a/)') &
               "WARNING: Parameter 'kerfix' overridden to 2 for &
               &SCF-MI calculation."
          pub_kerfix = 2
       end if


       ! mjsp: SCF MI approach: individual polarisation
       ! contributions from each of the fragments.
       ! Each polarisation contribution is
       ! calculated by SCF-MI of the individual fragments in the field of the
       ! potential resulting from all other fragments in their 'frozen' state.
       if (pub_eda_frag_isol_pol) then

          ! mjsp: Check not the case that we have loaded in continuation data
          ! mjsp: and calculation stopped at or beyond this stage
          if (.not.(pub_eda_continuation) .or. &
               (pub_eda_continuation .and. &
               (eda_cont_mode .le. EDA_POLFRAGLOC_DEVEL))) then

             ! Stoll locally projected approach:
             pub_eda_scfmi = .true.

             complex_energy_pol(0) = 0.0_DP

             pub_eda_mode = EDA_POLFRAGLOC_DEVEL

             ! loop the fragments:
             do pub_frag_counter=1,pub_frag_iatm(0)

                ! mjsp: check variational freedoms available (is not polarisable)
                if (sum(pub_frag_data(pub_frag_counter)%rep%n_occ(:,PUB_1K)) &
                     .eq. 0) then

                   if (pub_on_root) then
                      write(stdout,'(//a)') utils_banner('=')
                      write(stdout,'(a)') utils_banner('=', &
                           'Energy Decomposition Analysis (EDA)')
                      write(stdout,'(a/)') utils_banner('=')
                      write(stdout,'(a,I3/)') '      Skipping polarisation &
                           &calculation for fragment ID: ', pub_frag_counter
                      write(stdout,'(3x,a)') '=> Fragment has too few &
                           &variational freedoms (no occupied NGWFs).'
                      write(stdout,'(a/)') utils_banner('=')
                   end if

                   ! mjsp: No optimisation to perform, so energy of this state equals
                   ! mjsp: the frozen idempotent state's energy, to give
                   ! mjsp: zero polarisation overall
                   complex_energy_pol(pub_frag_counter) = complex_energy_frzIdem

                else ! is polarisable

                   if (pub_on_root) then
                      write(stdout,'(//a)') utils_banner('=')
                      write(stdout,'(a)') utils_banner('=', &
                           'Energy Decomposition Analysis (EDA)')
                      write(stdout,'(a/)') utils_banner('=')
                      write(stdout,'(a,I3/)') '      Beginning polarisation &
                           &calculation for fragment ID: ', pub_frag_counter
                      write(stdout,'(a/)') utils_banner('=')
                   end if

                   ! mjsp: avoid reloading continuation NGWFs and density kernel if
                   ! mjsp: has already been done in energy_and_force_calculate()
                   ! mjsp: by checking that we are ahead of
                   ! mjsp: the stage of the continuation file
                   if (.not. pub_eda_continuation .or. &
                        (pub_eda_continuation .and. &
                        ((pub_eda_mode .ne. eda_cont_mode) .and. &
                         (pub_frag_counter .ne. eda_cont_fc1)) ) ) &
                      ! mjsp: Reload the block-diagonal kernel as the initial guess
                      ! mjsp: kernel:
                      call internal_reload_frzidem_state()

                   ! mjsp: Construct the SCF-MI fragment masks with the appropriate
                   ! mjsp: pattern
                   call fmo_construct_masks()

!                   ! mjsp: freeze the off-diagonal blocks of the overlap matrix
!                   ! and recompute the inverse overlap.
!                   ! Rebuild the density-indep part of the Hamiltonian with the
!                   ! frozen overlap matrix
!                   ! (TODO: hacky --could be made more efficient)
!                   rep%inv_overlap_init = .false.
!                   rep%inv_overlap_b_init = .false.
!                   call hamiltonian_dens_indep_matrices(rep, ngwf_basis,proj_basis, &
!                        nl_projectors, hub_proj_basis, hub, mdl)! val_rep,val_ngwf_basis)

                   ! mjsp: optimise the kernel and NGWFs to find the polarised state.
                   call ngwf_cg_optimise(complex_energy_pol(pub_frag_counter), &
                        ngwfs_converged, out_of_runtime, ham, &
                        denskern, edft, rep, ngwf_nonsc_forces, lhxc_fine, ngwf_basis, &
                        proj_basis,nl_projectors, hub_proj_basis, hub, mdl, &
                        hfxstate, dfdtau_fine = eda_dfdtau_fine)

                   ! mjsp: NGWFs have changed, so update ngwfs_changed flag
                   ngwfs_changed = .true.

                   ! mjsp: Gather the total polarisation energy of all the fragments
                   complex_energy_pol(0) = complex_energy_pol(0) + &
                        complex_energy_pol(pub_frag_counter) - complex_energy_frzIdem

                   ! mjsp: If continuation write is true, then write data
                   if (pub_eda_write) &
                        call eda_super_restart_write_metadata(is_preptool=.false.)

                   ! mjsp: If electron density and/or electron density difference
                   ! mjsp: to be plotted
                   if (pub_eda_deltadens .or. eda_write_density_plot) then

                      if(pub_on_root) write(stdout,'(a,I3,a)') &
                           "EDA Density: Calculating density for fragment ID ", &
                           pub_frag_counter ," polarised idempotent state..."

                      ! mjsp: construct and store the polarisation state density
                      ! on the grid
                      call eda_store_density(eda_pol_density_fine(:,:,:,:, &
                           pub_frag_counter), rep, mdl, &
                           ngwf_basis, denskern_R%kern%m(:,PUB_1K))

                      ! construct the filename string
                      fnstr = ''
                      write( fnstr, '(a,i3.3,a)' ) '_pol_', pub_frag_counter, '_'

                      if (eda_write_density_plot) then
                         ! mjsp: output delta density in plot format file
                         ! vm: output density in Angstrom rather than in Bohr
                         ! jd: ... but leave the unit conversion to visual_scalarfield
                         call visual_scalarfield( &
                              eda_pol_density_fine(:,:,:,1,pub_frag_counter), &
                              mdl%fine_grid, mdl%cell, &
                              'Electronic density (in e/ang^3) for:', &
                              trim(fnstr)//'ed_density', &
                              mdl%elements, ANGSTROM**3)
                         if (pub_num_spins == 2) call visual_scalarfield( &
                              eda_pol_density_fine(:,:,:,2,pub_frag_counter), &
                              mdl%fine_grid, mdl%cell, &
                              'Electronic spin density (in e/ang^3) for:', &
                              trim(fnstr)//'ed_spindensity', &
                              mdl%elements, ANGSTROM**3)
                      end if

                      if (pub_eda_deltadens) then
                         ! output the frozen-to-polarisation delta density
                         call properties_plot_delta_density(trim(fnstr)//'edd_', &
                              eda_pol_density_fine(:,:,:,:,pub_frag_counter), &
                              eda_frzIdem_density_fine, mdl)
                      end if

                      if(pub_on_root) write(stdout,'(a)') "EDA Density: ...done"

                   end if  ! deltadens


                   ! mjsp: properties:
                   if (pub_do_properties) then

                      ! mjsp: Copy the proper representation of the full
                      ! mjsp: overlap matrix to the rep container:
                      call sparse_embed_create(olap_backup,rep%overlap)
                      call sparse_embed_create(inv_olap_backup,rep%inv_overlap)
                      call sparse_embed_copy(olap_backup,rep%overlap)
                      call sparse_embed_copy(inv_olap_backup,rep%inv_overlap)
                      call sparse_embed_copy(rep%overlap,rep%overlap_scfmi_full)
                      call sparse_embed_copy(rep%inv_overlap,rep%inv_overlap_scfmi_full)

                      ! mjsp: Compute the unprojected Hamiltonian
                      pub_eda_mode = EDA_POLFRAGLOC_DEVEL_OFF
                      call hamiltonian_build_matrix(ham,rep)

                      ! mjsp: Do the properties calculation(s):
                      call properties_calculate(denskern_R, ham, &
                           rep, ngwf_basis, proj_basis, nl_projectors, &
                           hub_proj_basis, hub, core_basis, core_wvfns, &
                           mdl, hfxstate, lhxc_fine, 'valence')

                      pub_eda_mode = EDA_POLFRAGLOC_DEVEL

                      ! mjsp: Revert overlap
                      call sparse_embed_copy(rep%overlap, olap_backup)
                      call sparse_embed_copy(rep%inv_overlap, inv_olap_backup)
                      call sparse_embed_destroy(olap_backup)
                      call sparse_embed_destroy(inv_olap_backup)

                   end if

                end if ! is polarisable check

             end do

          end if ! continuation check

       end if ! pub_eda_frag_isol_pol

       ! mjsp: Check not the case that we have loaded in continuation data
       ! mjsp: and calculation stopped at or beyond this stage
       if (.not.(pub_eda_continuation) .or. &
            (pub_eda_continuation .and. (eda_cont_mode .le. EDA_POLSIMUL))) then

          if (pub_on_root) then
             write(stdout,'(//a)') utils_banner('=')
             write(stdout,'(a)') utils_banner('=', &
                  'Energy Decomposition Analysis (EDA)')
             write(stdout,'(a/)') utils_banner('=')
             write(stdout,'(a/)')  '           Beginning polarisation &
                  &calculation for supermolecule...'
             write(stdout,'(a/)') utils_banner('=')
          end if

          ! mjsp: avoid reloading continuation NGWFs and density kernel if has
          ! mjsp: already has been done in energy_and_force_calculate() by checking
          ! mjsp: that we are ahead of the stage of the continuation file
          if (.not. pub_eda_continuation .or. &
               (pub_eda_continuation .and. &
               (pub_eda_mode .ne. eda_cont_mode))) then

             ! mjsp: if using the NGWFs calculated in the fragment calculations
             ! mjsp: as the initial guess of the polarised state's NGWFs:
             ! mjsp: (i.e. fully reinitialising the NGWFs)
             if (pub_eda_reset_ngwfs_pol) then

                ! mjsp: Initialise to initial-guess NGWFs:

                ! mjsp: Initialise the NGWF values (either from PAOs, STO-3G,
                ! mjsp: or by loading from disk (or memory)
                call ngwfs_initialise(rep,ngwf_basis,mdl,hfxstate)

             else

                ! mjsp: Initialise to optimised fragments' NGWFs:

                if (pub_eda_read_super) then

                   ! (note, code below is same as above, however
                   ! the state of pub_eda_reset_ngwfs_pol is accounted for in
                   ! ngwfs_initialise().  I've separated the logic out just
                   ! for clarity.)

                   ! mjsp: Initialise to the frozen density NGWFs from disk
                   call ngwfs_initialise(rep,ngwf_basis,mdl,hfxstate)

                else

                   ! mjsp: Reload the frozen density NGWFs from
                   ! mjsp: the optimized fragment spheres
                   ngwfs_changed = .true. ! jd: Ensures hamiltonian_dens_indep_matrices
                                          !     gets called below, regenerating
                                          !     overlap and inv_overlap.
                   call ngwfs_initialise_from_fragments(ngwf_basis, rep, mdl, &
                        hfxstate)

                end if

             end if

             ! mjsp: Reload the block-diagonal kernel as the initial guess
             ! mjsp: kernel:
             call sparse_embed_array_scale(denskern%kern, 0.0_DP)

             rep%inv_overlap_init = .false.
             rep%inv_overlap_b_init = .false.

             if (pub_eda_read_super) then
                call fragment_data_alloc_paw_core_densities(mdl)

                ! mjsp: Initialise the frozen density kernel from disk
                call electronic_init_denskern(rep, ham, denskern, &
                     lhxc_fine, ngwf_basis, proj_basis, nl_projectors, &
                     hub_proj_basis, hub, mdl, hfxstate, &
                     dfdtau_fine = eda_dfdtau_fine)
             else

                ! mjsp: Insert the elements of the internally stored fragment
                ! mjsp: density kernel into the appropriate elements of the
                ! mjsp: supermolecule density kernel.
                call fmo_load_internal_blkdiag_denskern(denskern%kern, ngwf_basis)
                call kernel_workspace_invalidate(denskern)

                ! Rebuild the density-indep part Hamiltonian with the
                ! full overlap matrix (TODO: hacky --could be made more efficient)
                if (ngwfs_changed) then
                   call hamiltonian_dens_indep_matrices(rep, ngwf_basis,proj_basis, &
                        nl_projectors, hub_proj_basis, hub, mdl)! val_rep,val_ngwf_basis)
                end if

             end if

          end if ! NGWFs/dkn require loading check

          ! mjsp: Continue with standard SCF MI approach (simultaneous
          ! mjsp: polarisation of all fragments)
          ! - If we performed fragment localised polarisation calculations, then
          ! perform standard SCF MI to calculate higher order polarisation
          ! remainder.
          ! - Else, perform standard SCF MI to calculate total polarisation.

          ! Simultaneous polarisation of all fragments:
          pub_frag_counter = 0
          pub_eda_mode = EDA_POLSIMUL

          ! Stoll locally projected approach:
          pub_eda_scfmi = .true.

          ! mjsp: Construct the SCF-MI fragment masks with the appropriate
          ! pattern
          call fmo_construct_masks()

!          ! mjsp: freeze the off-diagonal blocks of the overlap matrix
!          ! and recompute the inverse overlap.
!          ! Rebuild the density-indep part of the Hamiltonian with the frozen
!          ! overlap matrix (TODO: hacky --could be made more efficient)
!          rep%inv_overlap_init = .false.
!          rep%inv_overlap_b_init = .false.
!          call hamiltonian_dens_indep_matrices(rep, ngwf_basis,proj_basis, &
!               nl_projectors, hub_proj_basis, hub, mdl)! val_rep,val_ngwf_basis)

          ! mjsp: optimise the kernel and NGWFs to find the polarised state.
          call ngwf_cg_optimise(complex_energy_pol_simul, ngwfs_converged, &
               out_of_runtime, ham, &
               denskern, edft, rep, ngwf_nonsc_forces, lhxc_fine, ngwf_basis, &
               proj_basis,nl_projectors, hub_proj_basis, hub, mdl, hfxstate, &
               dfdtau_fine = eda_dfdtau_fine)

          ! mjsp: NGWFs have changed, so update ngwfs_changed flag
          ngwfs_changed = .true.

          ! mjsp: If continuation write is true, then write data
          if (pub_eda_write) &
             call eda_super_restart_write_metadata(is_preptool=.false.)

          ! mjsp: Plot NGWFs if requested
          if (pub_write_ngwf_plot .and. &
               (pub_cube_format .or. pub_grd_format .or. pub_dx_format)) &
             call visual_ngwfs(rep%ngwfs_on_grid(1), ngwf_basis(1), &
                  'pol_', mdl%regions(1)%elements, mdl%cell, &
                  mdl%fftbox, mdl%par)

          ! mjsp: If electron density and/or electron density difference
          ! mjsp: to be plotted
          if (pub_eda_deltadens .or. eda_write_density_plot) then

             if(pub_on_root) write(stdout,'(a)') "EDA Density: Calculating &
                  &density for the fully polarised idempotent state..."

             ! mjsp: construct and store the polarisation state density
             ! on the grid
             call eda_store_density(eda_pol_density_fine(:,:,:,:,0), &
                  rep, mdl, ngwf_basis, denskern_R%kern%m(:,PUB_1K))

             if (eda_write_density_plot) then
                ! mjsp: output delta density in plot format file
                ! vm: output density in Angstrom rather than in Bohr
                ! jd: ... but leave the unit conversion to visual_scalarfield
                call visual_scalarfield( &
                     eda_pol_density_fine(:,:,:,1,0), mdl%fine_grid, mdl%cell, &
                     'Electronic density (in e/ang^3) for:', &
                     '_pol_ed_density', &
                     mdl%elements, ANGSTROM**3)
                if (pub_num_spins == 2) call visual_scalarfield( &
                     eda_pol_density_fine(:,:,:,2,0), mdl%fine_grid, mdl%cell, &
                     'Electronic spin density (in e/ang^3) for:', &
                     '_pol_ed_spindensity', &
                     mdl%elements, ANGSTROM**3)
             end if

             if(pub_on_root) write(stdout,'(a)') "EDA Density: ...done"

             if (pub_eda_deltadens) then
                ! output the frozen-to-polarised delta density
                call properties_plot_delta_density('_pol_edd_', &
                   eda_pol_density_fine(:,:,:,:,0), eda_frzIdem_density_fine, &
                   mdl)
             end if

             ! if polarisation in a fragment-wise manner then also
             ! compute higher order polarisation delta density
             if (pub_eda_frag_isol_pol .and. pub_eda_deltadens) then

                ! - Effectively, the following accumulates the delta densities of
                ! all the fragments' polarisations into index 0 and
                ! adds 1x the frozen density to this to obtain the net fully
                ! polarised density of all the individual fragment polarisations.
                ! (i.e. frozen density + polarisation perturbation)
                if(pub_on_root) write(stdout,'(a)') "EDA Density: &
                     &Calculating delta density for the fully polarised &
                     &idempotent state (higher order polarisation)..."

                ! initialise
                eda_poliso_density_fine(:,:,:,:) = 0.0_DP

                ! accumulate the densities for the polarisations of all
                ! the fragments
                do it=1,pub_frag_iatm(0)
                   eda_poliso_density_fine(:,:,:,:) = &
                      eda_poliso_density_fine(:,:,:,:) + &
                      eda_pol_density_fine(:,:,:,:,it)
                end do

                ! subtract (number of fragments - 1)*frozen density to obtain the
                ! polarised density in index 0
                eda_poliso_density_fine(:,:,:,:) = &
                   eda_poliso_density_fine(:,:,:,:) - &
                   real((pub_frag_iatm(0)-1),kind=DP) * &
                   eda_frzIdem_density_fine(:,:,:,:)

                ! output the isolated-fragment-polarised-to-fully-polarised
                ! delta density
                call properties_plot_delta_density('_pol_higher_order_edd_', &
                   eda_pol_density_fine(:,:,:,:,0), &
                   eda_poliso_density_fine(:,:,:,:), &
                   mdl)

                if(pub_on_root) write(stdout,'(a)') "EDA Density: ...done"

                ! mjsp: higher order polarized density on grid is
                ! no longer required so deallocate
                call eda_destroy_poliso_density()

             end if

             ! mjsp: frozen density on grid is no longer required
             ! so deallocate
             call eda_destroy_frz_density()

          end if  ! deltadens

          ! mjsp: properties:
          if (pub_do_properties) then

             ! mjsp: Copy the proper representation of the density kernel and
             ! the full overlap matrix to the denskern and rep containers:
             call sparse_embed_create(olap_backup,rep%overlap)
             call sparse_embed_create(inv_olap_backup,rep%inv_overlap)
             call sparse_embed_copy(olap_backup,rep%overlap)
             call sparse_embed_copy(inv_olap_backup,rep%inv_overlap)
             call sparse_embed_copy(rep%overlap,rep%overlap_scfmi_full)
             call sparse_embed_copy(rep%inv_overlap,rep%inv_overlap_scfmi_full)

             ! mjsp: Compute the unprojected Hamiltonian
             pub_eda_mode = EDA_POLSIMUL_OFF
             call hamiltonian_build_matrix(ham,rep)

             ! mjsp: Do the properties calculation(s):
             call properties_calculate(denskern_R, ham, &
                  rep, ngwf_basis, proj_basis, nl_projectors, &
                  hub_proj_basis, hub, &
                  core_basis, core_wvfns, mdl, hfxstate, lhxc_fine, 'valence')

             pub_eda_mode = EDA_POLSIMUL

             ! mjsp: Revert overlap
             call sparse_embed_copy(rep%overlap, olap_backup)
             call sparse_embed_copy(rep%inv_overlap, inv_olap_backup)
             call sparse_embed_destroy(olap_backup)
             call sparse_embed_destroy(inv_olap_backup)

          end if

          ! mjsp: printout
          if (pub_on_root) then
             write(stdout,'(//a)') utils_banner('=')
             write(stdout,'(a)') utils_banner('=', &
                  'Energy Decomposition Analysis (EDA)')
             write(stdout,'(a)') utils_banner('=')
             write(stdout,'(a)') '==> Polarisation finished '

             ! individual fragment polarisations
             if (pub_eda_frag_isol_pol) then

                write(stdout,'(6x,a70)') '--------------&
                     &P-O-L-A-R-I-S-A-T-I-O-N----A-N-A-L-Y-S-I-S&
                     &----------------'

                write(stdout,'(51x,a11)') 'Fragment ID'
                do it=1,pub_frag_iatm(0)
                   write(stdout,'(19x,a31,i5,f10.3,a)') &
                        'E_POL energy :', it, &
                        (complex_energy_pol(it) - &
                        complex_energy_frzIdem)*627.509438736, ' kcal/mol'
                end do

                write(stdout,'(6x,a49,f10.3,a)') &
                     'Higher order E_POL energy (Stoll SCF-MI) :', &
                     (complex_energy_pol_simul - complex_energy_pol(0)&
                     -complex_energy_frzIdem)*627.509438736, ' kcal/mol'

                write(stdout,'(6x,a70)') utils_banner('-')

             end if

             write(stdout,'(6x,a49,f10.3,a)') &
                  'E_POL energy (Stoll SCF-MI) total :', &
                  (complex_energy_pol_simul - &
                  complex_energy_frzIdem)*627.509438736, ' kcal/mol'

             write(stdout,'(/11x,a/)') &
                  'Beginning charge transfer calculation for supermolecule...'

             write(stdout,'(a/)') utils_banner('=')

          end if

       end if ! continuation check



       ! mjsp: Check not the case that we have loaded in continuation data
       ! mjsp: and calculation stopped at or beyond this stage
       if (.not.(pub_eda_continuation) .or. &
            (pub_eda_continuation .and. &
            (eda_cont_mode .le. EDA_CTFRAGLOC_DEVEL))) then

          ! mjsp: if calculating (charge transfer) delocalisations
          ! between fragment pairs
          if (pub_eda_frag_isol_ct) then

             pub_eda_mode = EDA_CTFRAGLOC_DEVEL

             ! Stoll locally projected approach:
             pub_eda_scfmi = .true.

             ! loop the fragment pairs given in the input file:
             do it=deloc_iterator_start,SIZE(pub_frag_deloc(1,:))

                pub_frag_counter = pub_frag_deloc(1,it)
                pub_frag_counter2 = pub_frag_deloc(2,it)

                if (pub_on_root) then
                   write(stdout,'(//a)') utils_banner('=')
                   write(stdout,'(a)') utils_banner('=', &
                        'Energy Decomposition Analysis (EDA)')
                   write(stdout,'(a/)') utils_banner('=')
                   write(stdout,'(6x,a,I3,a,I3/)') &
                        'Beginning delocalisation calculation for fragment pair: ', &
                        pub_frag_counter, ' and ', pub_frag_counter2
                   write(stdout,'(a/)') utils_banner('=')
                end if

                ! mjsp: avoid reloading continuation NGWFs and density kernel if
                ! mjsp: already done in energy_and_force_calculate() by checking
                ! mjsp: that we are ahead of the stage of the continuation file
                if (.not. pub_eda_continuation .or. &
                     (pub_eda_continuation .and. &
                     ((pub_eda_mode .ne. eda_cont_mode) .and. &
                      (pub_frag_counter .ne. eda_cont_fc1) .and. &
                      (pub_frag_counter2 .ne. eda_cont_fc2)) ) ) then


                   if (pub_eda_read_super) then

                      ! mjsp: Initialise the NGWFs
                      ! mjsp: read_super=t: NGWFs handled by ngwfs_initialise
                      ! pub_eda_reset_ngwfs_ct == .true. => initial guess NGWFs
                      ! pub_eda_reset_ngwfs_ct == .false. => frozen state NGWFs
                      call ngwfs_initialise(rep,ngwf_basis,mdl,hfxstate)

                   else

                      if (pub_eda_reset_ngwfs_ct) then

                         ! mjsp: initial guess NGWFs
                         call ngwfs_initialise(rep,ngwf_basis,mdl,hfxstate)

                      else

                         ! mjsp: frozen state NGWFs
                         call ngwfs_initialise_from_fragments(ngwf_basis,rep,&
                              mdl,hfxstate)

                      end if

                   end if

                   call sparse_embed_array_scale(denskern%kern, 0.0_DP)

                   rep%inv_overlap_init = .false.
                   rep%inv_overlap_b_init = .false.

                   if (pub_eda_read_super) then
                      call fragment_data_alloc_paw_core_densities(mdl)

                      ! mjsp: Initialise to the frozen density kernel from disk
                      call electronic_init_denskern(rep, ham, denskern, &
                           lhxc_fine, ngwf_basis, proj_basis, nl_projectors, &
                           hub_proj_basis, hub, mdl, hfxstate, &
                           dfdtau_fine = eda_dfdtau_fine)

                   else

                      ! mjsp: Insert the elements of the internally stored
                      ! mjsp: fragment density kernel into the appropriate
                      ! mjsp: elements of the supermolecule density kernel.
                      call fmo_load_internal_blkdiag_denskern(denskern%kern, &
                           ngwf_basis)
                      call kernel_workspace_invalidate(denskern)

                      ! mjsp: Freeze the off-diagonal blocks of the overlap
                      ! mjsp: matrix and recompute the inverse overlap:
                      ! mjsp: rebuild the density-indep part of the Hamiltonian
                      ! mjsp: with the frozen overlap matrix.
                      ! (TODO: hacky --could be made more efficient)
                      if (ngwfs_changed) then
                         call hamiltonian_dens_indep_matrices(rep, ngwf_basis,&
                              proj_basis, nl_projectors, hub_proj_basis, &
                              hub, mdl)! val_rep,val_ngwf_basis)
                      end if

                   end if


                   ! mjsp: Construct the SCF-MI fragment masks with the appropriate
                   ! pattern
                   call fmo_construct_masks()

!                   ! mjsp: freeze the off-diagonal blocks of the overlap matrix
!                   ! and recompute the inverse overlap.
!                   ! --> hack: Rebuild the density-indep part of the Hamiltonian with the frozen overlap matrix
!                   ! (TODO: hacky --could be made more efficient)
!                   rep%inv_overlap_init = .false.
!                   rep%inv_overlap_b_init = .false.
!                   call hamiltonian_dens_indep_matrices(rep, ngwf_basis,proj_basis, &
!                        nl_projectors, hub_proj_basis, hub, mdl)! val_rep,val_ngwf_basis)

                   ! mjsp: optimise the kernel and NGWFs to find the
                   ! pair-delocalised state.
                   call ngwf_cg_optimise( &
                        complex_energy_ct(pub_frag_counter,pub_frag_counter2), &
                        ngwfs_converged, out_of_runtime, ham, denskern, edft, &
                        rep, ngwf_nonsc_forces, lhxc_fine, ngwf_basis, &
                        proj_basis,nl_projectors, hub_proj_basis, hub, mdl, &
                        hfxstate, dfdtau_fine = eda_dfdtau_fine)

                   ! mjsp: If continuation write is true, then write data
                   if (pub_eda_write) call eda_super_restart_write_metadata( &
                        is_preptool=.false.,ct_deloc_it=it)


                end if ! inner continuation check

             end do ! fragment loop

          end if ! fragment loop

       end if ! continuation check

       ! mjsp: turn off global fragment flags
       pub_frag_counter = 0
       pub_frag_counter2 = 0


       ! mjsp: Revert to exact LNV if temporarily turned off
       if (pub_on_root) write(stdout,*)
       if (exact_lnv_temp_off .eqv. .true.) then
          if (pub_on_root) write(stdout,'(a)') 'WARNING: Turning exact LNV back &
               &on for charge transfer calculation.'
          exact_lnv_temp_off = .false.
          pub_exact_lnv = .true.
       end if

       ! mjsp: Revert pub_kerfix if temporarily modified
       if (pub_kerfix_temp .ne. 2) then
          pub_kerfix = pub_kerfix_temp
          if(pub_on_root) write(stdout,'(a,i2)') "WARNING: Parameter 'kerfix' &
              &reverted to ", pub_kerfix
       end if
       if (pub_on_root) write(stdout,*)



       ! mjsp: (4) variational calculation of the fully relaxed state.
       pub_eda_mode = EDA_FULL
       pub_eda_scfmi = .false.

       if (pub_on_root) then
          write(stdout,'(//a)') utils_banner('=')
          write(stdout,'(a)') utils_banner('=', &
               'Energy Decomposition Analysis (EDA)')
          write(stdout,'(a)') utils_banner('=')
          write(stdout,'(a)') &
               '==> Beginning fully relaxed supermolecule calculation '
          write(stdout,'(a/)') utils_banner('=')
       end if


       ! mjsp: if using the NGWFs calculated in the fragment calculations as
       ! the initial guess of the fully optimised state NGWFs
       ! (i.e. fully reinitialising the NGWFs)
       if (pub_eda_reset_ngwfs_ct) then

          ! mjsp: Initialise to initial-guess NGWFs:
          call ngwfs_initialise(rep,ngwf_basis,mdl,hfxstate)

          ! mjsp: Reload the block-diagonal kernel as the initial guess kernel
          call sparse_embed_array_scale(denskern%kern, 0.0_DP)

          rep%inv_overlap_init = .false.
          rep%inv_overlap_b_init = .false.

          if (pub_eda_read_super) then
             call fragment_data_alloc_paw_core_densities(mdl)

             ! mjsp: Initialise the frozen density kernel from disk
             call electronic_init_denskern(rep, ham, denskern, &
                  lhxc_fine, ngwf_basis, proj_basis, nl_projectors,&
                  hub_proj_basis, hub, mdl, hfxstate, &
                  dfdtau_fine = eda_dfdtau_fine)
          else

             ! mjsp:Insert the elements of the internally stored fragment density
             ! kernel into the appropriate elements of the supermolecule density kernel.
             call fmo_load_internal_blkdiag_denskern(denskern%kern, ngwf_basis)
             call kernel_workspace_invalidate(denskern)

             ! Rebuild the density-indep part Hamiltonian with the full overlap matrix
             ! (TODO: hacky --could be made more efficient)
             call hamiltonian_dens_indep_matrices(rep, ngwf_basis,proj_basis, &
                  nl_projectors, hub_proj_basis, hub, mdl)! val_rep,val_ngwf_basis)

          end if

       ! mjsp: if using the NGWFs calculated in the polarisation calculation as
       ! the initial guess of the fully optimised state NGWFs:
       else

          ! mjsp: The previous Stoll equations used to calculate the polarisation
          ! contributions used a block diagonal form of the overlap matrix.
          ! mjsp: Copy the proper representation of the density kernel in the full
          ! overlap matrix to the density kernel:
          call sparse_embed_array_copy(denskern%kern, denskern_R%kern)

          ! Rebuild the density-indep part Hamiltonian with the full overlap matrix
          ! (TODO: hacky --could be made more efficient)
          rep%inv_overlap_init = .false.
          rep%inv_overlap_b_init = .false.
          call hamiltonian_dens_indep_matrices(rep, ngwf_basis,proj_basis, &
               nl_projectors, hub_proj_basis, hub, mdl)! val_rep,val_ngwf_basis)

       end if


       ! optimise density kernel and NGWF coefficients
       if ((pub_maxit_ngwf_cg .gt. 0).or.pub_nonsc_forces) then
          call ngwf_cg_optimise(complex_energy_full, ngwfs_converged, &
               out_of_runtime, ham, &
               denskern, edft, rep, ngwf_nonsc_forces, lhxc_fine, ngwf_basis, &
               proj_basis,nl_projectors, hub_proj_basis, hub, mdl, hfxstate, &
               dfdtau_fine = eda_dfdtau_fine)
       end if

       ! mjsp: If continuation write is true, then write data
       if (pub_eda_write) call eda_super_restart_write_metadata(is_preptool=.false.)


       ! mjsp: Plot NGWFs if requested
       if (pub_write_ngwf_plot .and. &
            (pub_cube_format .or. pub_grd_format .or. pub_dx_format)) &
          call visual_ngwfs(rep%ngwfs_on_grid(1), ngwf_basis(1), &
               'ct_', mdl%regions(1)%elements, mdl%cell, &
               mdl%fftbox, mdl%par)


       ! mjsp: If electron density and/or electron density difference
       ! mjsp: to be plotted
       if (pub_eda_deltadens .or. eda_write_density_plot) then

          if(pub_on_root) write(stdout,'(a)') "EDA Density: Calculating density &
               &for fully variationally optimized state..."

          ! mjsp: construct and store the charge transfer state density
          ! on the grid for delta charge transfer density calculation
          call eda_store_density(eda_full_density_fine, rep, mdl, ngwf_basis, &
               denskern%kern%m(:,PUB_1K))

          if (eda_write_density_plot) then
             ! mjsp: output delta density in plot format file
             ! vm: output density in Angstrom rather than in Bohr
             ! jd: ... but leave the unit conversion to visual_scalarfield
             call visual_scalarfield( &
                  eda_full_density_fine(:,:,:,1), mdl%fine_grid, mdl%cell, &
                  'Electronic density (in e/ang^3) for:', &
                  '_relaxed_ed_density', &
                  mdl%elements, ANGSTROM**3)
             if (pub_num_spins == 2) call visual_scalarfield( &
                  eda_full_density_fine(:,:,:,2), mdl%fine_grid, mdl%cell, &
                  'Electronic spin density (in e/ang^3) for:', &
                  '_relaxed_ed_spindensity', &
                  mdl%elements, ANGSTROM**3)
          end if

          if (pub_eda_deltadens) then
             ! mjsp: Output the polarised-to-fully-delocalised delta density:
             call properties_plot_delta_density('_ct_edd_', &
                  eda_full_density_fine, eda_pol_density_fine(:,:,:,:,0), mdl)
          end if

          if(pub_on_root) write(stdout,'(a)') "EDA Density: ...done"

          ! mjsp: polarized and fully optimized densities on grid are
          ! no longer required so deallocate
          call eda_destroy_pol_density()
          call eda_destroy_full_density()

       end if  ! deltadens

    end if  ! preptool

!    ! mjsp: Initialisation for NGWF reloading:
!    if ((pub_eda_reset_ngwfs_pol) .or. (pub_eda_read_super)) then
!        ! mjsp: deallocate storage for radial densities
!        if (pub_initial_dens_realspace) call density_init_radial_exit(mdl)
!    end if

  contains

    subroutine internal_reload_frzidem_state()

      ! Reloads the kernel and NGWFs for initialising the polarisation EDA
      ! stage (to the frozen idempotent state)

      if (pub_eda_read_super) then
         ! mjsp: Initialise the NGWFs
         ! mjsp: read_super=t: NGWFs handled by ngwfs_initialise
         ! mjsp: pub_eda_reset_ngwfs_pol=T => Initialise to guess NGWFs
         ! mjsp: pub_eda_reset_ngwfs_pol=F => Initialise to the frozen density NGWFs
         call ngwfs_initialise(rep,ngwf_basis,mdl,hfxstate)
      else
         ! mjsp: read_super=f: NGWFs handled here:
         if (pub_eda_reset_ngwfs_pol) then
            ! mjsp: pub_eda_reset_ngwfs_pol=T => Initialise to guess NGWFs
            call ngwfs_initialise(rep,ngwf_basis,mdl,hfxstate)
         else
            ! mjsp: pub_eda_reset_ngwfs_pol=F => Initialise to the frozen density NGWFs:
            ! mjsp: Reload the frozen density NGWFs from the optimized fragment spheres
            call ngwfs_initialise_from_fragments(ngwf_basis, rep, mdl, hfxstate)
         end if
      end if

      call sparse_embed_array_scale(denskern%kern, 0.0_DP)

      rep%inv_overlap_init = .false.
      rep%inv_overlap_b_init = .false.

      if (pub_eda_read_super) then
         call fragment_data_alloc_paw_core_densities(mdl)

         ! mjsp: Initialise the frozen density kernel from disk
         call electronic_init_denskern(rep, ham, denskern, &
              lhxc_fine, ngwf_basis, proj_basis, nl_projectors,&
              hub_proj_basis, hub, mdl, hfxstate, &
              dfdtau_fine = eda_dfdtau_fine)
      else

         ! mjsp:Insert the elements of the internally stored fragment density
         ! kernel into the appropriate elements of the supermolecule density kernel.
         call fmo_load_internal_blkdiag_denskern(denskern%kern, ngwf_basis)
         call kernel_workspace_invalidate(denskern)

         ! mjsp: freeze the off-diagonal blocks of the overlap matrix
         ! and recompute the inverse overlap.
         ! Rebuild the density-indep part of the Hamiltonian with the frozen
         ! overlap matrix (TODO: hacky --could be made more efficient)
         if (ngwfs_changed) then
            call hamiltonian_dens_indep_matrices(rep, ngwf_basis,proj_basis, &
                 nl_projectors, hub_proj_basis, hub, mdl)! val_rep,val_ngwf_basis)
         end if

      end if

    end subroutine internal_reload_frzidem_state

  end subroutine eda_driver_supermol_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eda_driver_prep_supermol(mdl, hfxstate, rep, ngwf_basis)

    !======================================================================!
    ! Prepares the supermolecule-fragment relationship arrays, and the     !
    ! supermolecule NGWFs.                                                 !
    !----------------------------------------------------------------------!
    ! Written by Max Phipps January 2016 (from energy_and_force_calculate) !
    ! Modified for embedding by Joseph Prentice, September 2018            !
    !======================================================================!

    use comms, only: pub_my_proc_id
    use dense, only: dense_create
#ifdef SCALAPACK
    use dense, only: dense_redist
#endif
    use fragment_data, only: NGWF_INDEX_MAP_TYPE, pub_frag_data, &
         pub_super2fragid_on_proc, pub_super2frag_ngwf_log, &
         fragment_data_get_ngwf_index_map
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use ngwfs, only: ngwfs_initialise, ngwfs_initialise_from_fragments
    use rundat, only: pub_num_spins, pub_num_kpoints, &
         pub_eda_continuation, pub_eda_read_super, pub_frag_iatm
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(MODEL), intent(inout)      :: mdl ! fragment model
    type(HFX_STATE), intent(inout), target :: hfxstate
    type(FUNC_BASIS), intent(inout) :: ngwf_basis(1) ! fragment NGWFs
    type(NGWF_REP), intent(inout)   :: rep ! fragment representation

    ! Local Variables
    integer :: ingwf, it, ierr, is
    integer :: ingwf_frag, ingwf_super, super_proc, ingwf_super_onproc

    type(NGWF_INDEX_MAP_TYPE) :: ngwf_index_map

    ! Construct the supermolecule atom to fragment id referencing array
    ! (proc-specific)
    if (.not. allocated(pub_super2fragid_on_proc)) then
       allocate(pub_super2fragid_on_proc(ngwf_basis(1)%proc_num),stat=ierr)
       call utils_alloc_check('eda_driver_prep_supermol','pub_super2fragid_on_proc', &
            ierr)
    end if

    ! Allocate n_occ according to public numbers of spins and k-points.
    allocate(pub_frag_data(0)%rep%n_occ(pub_num_spins, pub_num_kpoints), &
         stat=ierr)
    call utils_alloc_check('eda_driver_prep_supermol', 'pfd(0)%rep%n_occ', ierr)

    ! initialise number of occupied MOs
    pub_frag_data(0)%rep%n_occ = rep%n_occ

    if (pub_eda_read_super) then

       ! Allocate n_occ according to public numbers of spins and k-points.
       do it=1,pub_frag_iatm(0) ! loop fragments
          allocate(pub_frag_data(it)%rep%n_occ(pub_num_spins, pub_num_kpoints), &
               stat=ierr)
          call utils_alloc_check('eda_driver_prep_supermol', &
               'pfd(it)%rep%n_occ', ierr)
       end do

       ! mjsp: Initialise the frozen density NGWFs from disk
       call ngwfs_initialise(rep,ngwf_basis,mdl,hfxstate)

       ! mjsp: Read in the required EDA calculation metadata from the '.eda' file
       if (pub_eda_continuation) then
          ! mjsp: Continuation: read from continuation file
          call eda_super_restart_read_metadata( &
               mdl, &
               ngwf_basis, &
               is_preptool = .false.)
       else
          ! mjsp: Else read from preptool file
          call eda_super_restart_read_metadata( &
               mdl, &
               ngwf_basis, &
               is_preptool = .true.)
       end if

       ngwf_index_map = fragment_data_get_ngwf_index_map()

       ! mjsp: Populate the pub_super2fragid_on_proc array
       ! using the pub_frag2super_ngwf_idx array we have just loaded:
       ingwf = 1
       do it=1,pub_frag_iatm(0) ! loop fragments

          do ingwf_frag=1,pub_frag_data(it)%ngwf_basis(1)%num  ! loop NGWFs on fragment

             ! mjsp: obtain ingwf in the current supermolecular
             ! NGWF distribution:
             ingwf_super = ngwf_index_map%supers_by_cum_frag(ingwf)

             ! mjsp: identify the proc this fragment NGWF's data is on
             super_proc = ngwf_basis(1)%proc_of_func(ingwf_super)

             ! mjsp: get the array index of the NGWF after having been
             ! distributed to its proc:
             ingwf_super_onproc = ingwf_super - &
                ngwf_basis(1)%first_on_proc(super_proc) + 1

             if ( pub_my_proc_id == super_proc ) &
                pub_super2fragid_on_proc(ingwf_super_onproc) = it

             ingwf = ingwf + 1

          end do

       end do

    else

       ! mjsp: Initialise the frozen density NGWFs from the internally
       ! stored fragment NGWFs (frags->super)
       call ngwfs_initialise_from_fragments(ngwf_basis, rep, mdl, hfxstate)

       ngwf_index_map = fragment_data_get_ngwf_index_map()

    end if

    ! =========

    ! mjsp: Construct the fragment NGWF logical table, pub_super2frag_ngwf_log
    ! [NGWF,frag]
    allocate(pub_super2frag_ngwf_log(ngwf_basis(1)%num,pub_frag_iatm(0)),stat=ierr)
    call utils_alloc_check('eda_driver_prep_supermol','pub_super2frag_ngwf_log', &
         ierr)

    ! mjsp: Initialise array as no fragment NGWFs on all tables,
    ! mjsp: and modify to identify fragment NGWFs on each table.
    pub_super2frag_ngwf_log(:,:) = .false.

    ingwf = 1
    do it=1,pub_frag_iatm(0) ! loop fragments
       do ingwf_frag=1,pub_frag_data(it)%ngwf_basis(1)%num  ! loop NGWFs on fragment
          ! mjsp: Get supermolecule NGWF index
          ingwf_super = ngwf_index_map%supers_by_cum_frag(ingwf)
          ! mjsp: Update pub_super2frag_ngwf_log to identify the fragment NGWF
          pub_super2frag_ngwf_log(ingwf_super,it) = .true.
          ingwf = ingwf + 1
       end do
    end do


#ifdef SCALAPACK
    ! mjsp: If internally storing density kernels:
    if (.not.(pub_eda_read_super)) then

       ! mjsp: Distribute the fragment density kernels over the procs
       ! mjsp: (for later construction of the supermolecule
       ! mjsp: density kernel from the fragment density kernels)
       do it=1,pub_frag_iatm(0) ! loop fragments

          ! mjsp: Allocate ScaLAPACK distributed fragment density kernels
          allocate(pub_frag_data(it)%denskern_dens(pub_num_spins), stat=ierr)
          call utils_alloc_check('energy_and_force_calculate','pub_frag_data(it)%denskern_dens',ierr)

          ! mjsp: redistribute fragment DEM density kernels (from root)
          do is=1,pub_num_spins

             call dense_create(pub_frag_data(it)%denskern_dens(is), &
                               pub_frag_data(it)%mdl%par%num_ngwfs, &
                               pub_frag_data(it)%mdl%par%num_ngwfs)
             call dense_redist(pub_frag_data(it)%denskern_dens_loc(is),pub_frag_data(it)%denskern_dens(is))

          end do

       end do
    end if
#endif

  end subroutine eda_driver_prep_supermol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eda_allocate_densities(mdl)

    ! Allocates the intermediate (frozen, polarized and fully optimized)
    ! densities for later use in constructing the delta densities on grid.
    ! Written by Max Phipps Feb 2015

    use rundat, only: pub_num_spins, pub_frag_iatm, pub_eda_frag_isol_pol
    use utils, only: utils_alloc_check
    use fragment_data, only: eda_frzIdem_density_fine, eda_pol_density_fine, &
       eda_poliso_density_fine, eda_full_density_fine

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl

    ! Local Variables
    integer :: ierr

    ! mjsp: Allocate frozen density storage
    allocate(eda_frzIdem_density_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
    call utils_alloc_check('eda_allocate_densities', &
         'eda_frzIdem_density_fine',ierr)

    ! mjsp: If isolated fragment polarisation analysis included within EDA
    if (pub_eda_frag_isol_pol) then

       ! mjsp: Allocate polarized density storage: one density for each fragment
       allocate(eda_pol_density_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12,pub_num_spins,0:pub_frag_iatm(0)),stat=ierr)
       call utils_alloc_check('eda_allocate_densities', &
            'eda_pol_density_fine',ierr)

       ! mjsp: Allocate net polarized fragments density storage
       allocate(eda_poliso_density_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('eda_allocate_densities', &
            'eda_poliso_density_fine',ierr)

    else

       ! mjsp: Allocate polarized density storage (only one element required for
       ! final index due to not performing isolated fragment polarisation
       ! calculations)
       allocate(eda_pol_density_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12,pub_num_spins,0:1),stat=ierr)
       call utils_alloc_check('eda_allocate_densities', &
            'eda_pol_density_fine',ierr)

    end if

    ! mjsp: Allocate fully optimized density storage
    allocate(eda_full_density_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
    call utils_alloc_check('eda_allocate_densities', &
         'eda_full_density_fine',ierr)

  end subroutine eda_allocate_densities

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eda_store_density(density_fine, rep, mdl, ngwf_basis, denskern)

    ! Constructs and stores densities on the grid for later
    ! use in constructing the delta densities on grid.
    ! Written by Max Phipps Feb 2015
    ! Modified for embedding by Joseph Prentice, September 2018

    use sparse, only: SPAM3
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_scale, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use density, only: density_on_grid
    use rundat, only: pub_num_spins, pub_spin_fac
    use ngwf_representation, only: NGWF_REP
    use function_basis, only: FUNC_BASIS
    use constants, only: UP, DN

    implicit none

    ! Arguments
    ! the density to store:
    real(kind=DP), intent(inout) :: density_fine(:,:,:,:)
    type(NGWF_REP), intent(in) :: rep
    type(MODEL), intent(in) :: mdl
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    type(SPAM3_EMBED), intent(inout) :: denskern(pub_num_spins)

    ! Local variables
    type(SPAM3) :: kern_array(pub_num_spins)

    ! mjsp: Scale denskern for spin degeneracy
    if (pub_num_spins == 1) then
       call sparse_embed_scale(denskern(1),pub_spin_fac)
    end if

    call sparse_embed_extract_from_array(kern_array,denskern)

    ! mjsp: Calculate data-parallelised charge density
    call density_on_grid(density_fine,mdl%fine_grid,mdl%dbl_grid, &
         mdl%cell, mdl%fftbox, kern_array, rep%ngwf_overlap%p, &
         rep%ngwfs_on_grid(1), ngwf_basis(1), &
         rep%ngwfs_on_grid(1), ngwf_basis(1))

    call sparse_embed_destroy_extracted_array(kern_array)

    ! mjsp: spin polarisation: calculate total density and spin density
    if (pub_num_spins == 2) then
       density_fine(:,:,:,1) = density_fine(:,:,:,UP) + &
            density_fine(:,:,:,DN)
       density_fine(:,:,:,2) = density_fine(:,:,:,UP) - &
            2.0_DP * density_fine(:,:,:,DN)
    end if


    ! mjsp: restore normalisation of denskern
    if (pub_num_spins == 1) then
       call sparse_embed_scale(denskern(1),1.0_DP/pub_spin_fac)
    end if

  end subroutine eda_store_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eda_destroy_frz_density()

    ! Destroys the stored frozen density on the grid
    ! used to construct the delta polarisation density on grid.
    ! Written by Max Phipps Feb 2015

    use utils, only: utils_dealloc_check
    use fragment_data, only: eda_frzIdem_density_fine

    implicit none

    ! Local Variables
    integer :: ierr

    ! mjsp: Deallocate frozen density storage
    deallocate(eda_frzIdem_density_fine, stat=ierr)
    call utils_dealloc_check('eda_destroy_frz_density', &
         'eda_frzIdem_density_fine',ierr)

  end subroutine eda_destroy_frz_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eda_destroy_pol_density()

    ! Destroys the stored polarized density on the grid
    ! used to construct the delta polarisation and charge transfer
    ! densities on grid.
    ! Written by Max Phipps Feb 2015

    use utils, only: utils_dealloc_check
    use fragment_data, only: eda_pol_density_fine

    implicit none

    ! Local Variables
    integer :: ierr

    ! mjsp: Deallocate frozen density storage
    deallocate(eda_pol_density_fine, stat=ierr)
    call utils_dealloc_check('eda_destroy_pol_density', &
         'eda_pol_density_fine',ierr)

  end subroutine eda_destroy_pol_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eda_destroy_poliso_density()

    ! Destroys the stored net polarized density constructed from the individual
    ! polarised fragment densities used to construct the higher order
    ! polarisation delta density.
    ! Written by Max Phipps July 2015

    use utils, only: utils_dealloc_check
    use fragment_data, only: eda_poliso_density_fine

    implicit none

    ! Local Variables
    integer :: ierr

    ! mjsp: Deallocate frozen density storage
    deallocate(eda_poliso_density_fine, stat=ierr)
    call utils_dealloc_check('eda_destroy_poliso_density', &
         'eda_poliso_density_fine',ierr)

  end subroutine eda_destroy_poliso_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eda_destroy_full_density()

    ! Destroys the stored fully optimized density on the grid
    ! used to construct the delta charge transfer density on grid.
    ! Written by Max Phipps Feb 2015

    use utils, only: utils_dealloc_check
    use fragment_data, only: eda_full_density_fine

    implicit none

    ! Local Variables
    integer :: ierr

    ! mjsp: Deallocate frozen density storage
    deallocate(eda_full_density_fine, stat=ierr)
    call utils_dealloc_check('eda_destroy_full_density', &
         'eda_full_density_fine',ierr)

  end subroutine eda_destroy_full_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eda_super_restart_write_metadata(ct_deloc_it,is_preptool)

    !======================================================================!
    ! This subroutine writes out the continuation file data. This is data  !
    ! for the supermolecule, i.e. metadata relating to the fragment/       !
    ! supermolecule/intermediate state NGWFs, kernel, energies and other   !
    ! metadata.                                                            !
    !----------------------------------------------------------------------!
    ! Written by Max Phipps September 2015.                                !
    ! Modified by Max Phipps December 2015 to improve usability.           !
    ! Modified for embedding by Joseph Prentice, September 2018            !
    !======================================================================!

    use comms, only: pub_on_root
    use rundat, only: pub_eda_isol, pub_eda_isol_x, pub_eda_isol_c, &
         pub_eda_frz_x, pub_eda_frz_c, pub_eda_rep_x, pub_eda_rep_c, &
         pub_eda_intrafrag_x_frz, pub_eda_intrafrag_c_frz, pub_eda_dE_x_rep, &
         pub_eda_dE_c_rep, pub_frag_iatm, &
         pub_frag_counter, pub_frag_counter2, &
         pub_frag_deloc, pub_eda_preptool, &
         pub_eda_frag_isol_pol, pub_eda_frag_isol_ct, &
         pub_eda_mode, pub_rootname
    use fragment_data, only: NGWF_INDEX_MAP_TYPE, pub_frag_data, &
         complex_energy_frz, complex_energy_frzIdem, &
         complex_energy_pol_simul, complex_energy_pol, complex_energy_ct, &
         fragment_data_get_ngwf_index_map
    use utils, only: utils_unit, utils_open_unit_check, utils_close_unit_check

    implicit none

    ! Arguments (inputs)
    logical           :: is_preptool ! whether we are writing preptool data
    integer, optional :: ct_deloc_it ! delocalisation calculation iteration of
                                     ! pub_frag_deloc array (from input file)

    ! Local Variables
    integer           :: it
    integer           :: io_unit,io_stat
    integer           :: frag1, frag2
    character(len=80) :: filename

    type(NGWF_INDEX_MAP_TYPE) :: ngwf_index_map

    ngwf_index_map = fragment_data_get_ngwf_index_map()

    ! mjsp: Write out the EDA data
    if (pub_on_root) then

       ! Find available unit specifier
       io_unit = utils_unit()

       if (is_preptool) then
          ! mjsp: If preptool
          write(filename,'(a,a)') trim(pub_rootname),'.eda'

          write(stdout,'(a,a)') 'Writing EDA preparation data &
              &for supermolecule to file ',filename
       else
          ! mjsp: Else, continuation file
          write(filename,'(a,a)') trim(pub_rootname),'.eda_continuation'

          write(stdout,'(a,a)') 'Writing EDA continuation data &
              &for supermolecule to file ',filename
       end if

       open(unit=io_unit,iostat=io_stat,file=filename,status='REPLACE',&
            form='UNFORMATTED',position='REWIND',action='WRITE')
       call utils_open_unit_check('eda_super_restart_write_metadata',&
            filename,io_stat)

       ! Write the data
       ! Unformatted:

       ! mjsp: Write stage we are up to
       write(io_unit) pub_eda_preptool
       write(io_unit) pub_eda_mode, pub_frag_counter, pub_frag_counter2

       ! mjsp: EDA_FROZENNONIDEM or EDA_PREP data
       write(io_unit) pub_eda_isol(0)
       do it=1,pub_frag_iatm(0)
          write(io_unit) pub_eda_isol(it), pub_eda_isol_x(it), pub_eda_isol_c(it)
          write(io_unit) pub_frag_data(it)%ngwf_basis(1)%num
          write(io_unit) pub_frag_data(it)%rep%n_occ
       end do
       write(io_unit) complex_energy_frz
       write(io_unit) pub_eda_frz_x, pub_eda_frz_c
       !write(io_unit) pub_eda_rep_x, pub_eda_rep_c
       write(io_unit) pub_eda_intrafrag_x_frz,  pub_eda_intrafrag_c_frz
       !write(io_unit) pub_eda_dE_x_rep, pub_eda_dE_c_rep

       write(io_unit) size(ngwf_index_map%atom_ids_by_cum_frag)
       write(io_unit) ngwf_index_map%atom_ids_by_cum_frag

       ! If not writing to an EDA_PREP calculation file
       if (.not.(pub_eda_preptool)) then

          ! Frozen, antisymmetrized
          write(io_unit) complex_energy_frzIdem
          write(io_unit) pub_eda_dE_x_rep, pub_eda_rep_x, pub_eda_frz_x
          write(io_unit) pub_eda_dE_c_rep, pub_eda_rep_c, pub_eda_frz_c

          ! Polarisation:
          write(io_unit) complex_energy_pol_simul
          ! individual fragment polarisations
          if (pub_eda_frag_isol_pol) then
             do it=1,pub_frag_iatm(0)
                write(io_unit) complex_energy_pol(it)
             end do
          end if

          ! Charge transfer:
          ! individual fragment-pair delocalisations
          if (pub_eda_frag_isol_ct) then
             ! loop the fragment pairs given in the input file:
             do it=1,SIZE(pub_frag_deloc(1,:))
                frag1 = pub_frag_deloc(1,it)
                frag2 = pub_frag_deloc(2,it)
                write(io_unit) complex_energy_ct(frag1,frag2)
             end do
             ! write the charge transfer delocalisation stage
             ! iteration we are on (if we are
             ! currently at the CT delocalisations stage)
             if (present(ct_deloc_it)) then
                write(io_unit) ct_deloc_it
             else
                write(io_unit) 0
             end if
          end if

       end if ! preptool

       ! Close the file
       close(io_unit,iostat=io_stat)
       call utils_close_unit_check('eda_super_restart_write_metadata',&
                                    filename,io_stat)

    end if ! pub_on_root

  end subroutine eda_super_restart_write_metadata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eda_super_restart_read_metadata(mdl, ngwf_basis, is_preptool)

    !======================================================================!
    ! This subroutine reads in data from the EDA continuation file.        !
    !----------------------------------------------------------------------!
    ! Written by Max Phipps September 2015.                                !
    ! Modified by Max Phipps December 2015 to improve EDA usability.       !
    ! Modified for embedding by Joseph Prentice, September 2018            !
    !======================================================================!

    use comms, only: pub_on_root, comms_bcast, pub_root_proc_id
    use constants, only: EDA_FROZENIDEM, EDA_POLSIMUL, &
         EDA_POLFRAGLOC_DEVEL, EDA_CTFRAGLOC_DEVEL, EDA_CTFRAGLOC_DEVEL_OFF, &
         EDA_POLSIMUL_OFF, EDA_POLFRAGLOC_DEVEL_OFF
    use fragment_data, only: pub_frag_data, complex_energy_frz, &
         complex_energy_frzIdem, complex_energy_pol_simul, &
         complex_energy_pol, complex_energy_ct, &
         fragment_data_set_ngwf_index_map
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_eda_isol, pub_eda_isol_x, pub_eda_isol_c, &
         pub_eda_frz_x, pub_eda_frz_c, pub_eda_rep_x, pub_eda_rep_c, &
         pub_eda_intrafrag_x_frz, pub_eda_intrafrag_c_frz, pub_eda_dE_x_rep, &
         pub_eda_dE_c_rep, pub_frag_iatm, &
         pub_frag_counter, pub_frag_counter2, &
         pub_frag_deloc, &
         pub_eda_frag_isol_pol, pub_eda_frag_isol_ct, &
         pub_eda_mode, pub_eda_scfmi, &
         pub_super_file_prefix, pub_rootname
    use services, only: services_flush
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_close_unit_check, utils_open_unit_check, utils_unit

    implicit none

    ! Arguments (inputs)
    type(MODEL), intent(in)      :: mdl
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)  ! Supermolecular NGWF basis.
    logical                      :: is_preptool    ! whether we are reading
                                                   ! preptool data

    ! Local Variables
    integer           :: it
    integer           :: io_unit,io_stat
    integer           :: frag1, frag2
    character(len=80) :: filename
    logical           :: preptool_value ! value of preptool logical in the EDA file

    integer :: num_cum_frag_ngwfs
    integer, allocatable :: atom_ids_by_cum_frag_ngwf(:)

    integer :: ierr
    character(255) :: errmsg = ""

    ! Allocate the array that stores NGWF bases for each fragment if the array
    ! is not already allocated.
    do it = 1, pub_frag_iatm(0)
       if (.not. allocated(pub_frag_data(it)%ngwf_basis)) then
          allocate( &
               pub_frag_data(it)%ngwf_basis(1), &
               stat = ierr, &
               errmsg = errmsg)

          call utils_alloc_check( &
               "eda_super_restart_read_metadata", &
               trim(errmsg), &
               ierr)
       end if
    end do

    ! mjsp: Write out the EDA data
    if (pub_on_root) then

       ! mjsp: Find available unit specifier
       io_unit = utils_unit()

       if (is_preptool) then
          ! mjsp: If preptool, then filename is given by eda_super block
          write(filename,'(a,a)') trim(pub_super_file_prefix),'.eda'

          write(stdout,'(a,a)',advance='no') 'Reading EDA preparation data &
              &for supermolecule from file ',trim(filename),' ...'
       else
          ! mjsp: Else, continuation file within this directory.
          write(filename,'(a,a)') trim(pub_rootname),'.eda_continuation'

          write(stdout,'(a,a)',advance='no') 'Reading EDA continuation data &
              &for supermolecule from file ',trim(filename),' ...'
       end if


       open(unit=io_unit,iostat=io_stat,file=filename,&
            access='SEQUENTIAL',form='UNFORMATTED',position='REWIND', &
            action='READ')
       call utils_open_unit_check('eda_super_restart_read_metadata',&
            filename,io_stat)

       ! Read the data
       ! Unformatted:

       ! mjsp: Read the stage we are up to
       read(io_unit) preptool_value
       read(io_unit) eda_cont_mode, eda_cont_fc1, eda_cont_fc2
       pub_eda_mode = eda_cont_mode
       pub_frag_counter = eda_cont_fc1
       pub_frag_counter2 = eda_cont_fc2

       ! mjsp: EDA_FROZENNONIDEM or EDA_PREP data
       read(io_unit) pub_eda_isol(0)
       do it=1,pub_frag_iatm(0)
          read(io_unit) pub_eda_isol(it), pub_eda_isol_x(it), pub_eda_isol_c(it)
          read(io_unit) pub_frag_data(it)%ngwf_basis(1)%num
          read(io_unit) pub_frag_data(it)%rep%n_occ
       end do
       read(io_unit) complex_energy_frz
       read(io_unit) pub_eda_frz_x, pub_eda_frz_c
       !read(io_unit) pub_eda_rep_x, pub_eda_rep_c
       read(io_unit) pub_eda_intrafrag_x_frz, &
         pub_eda_intrafrag_c_frz
       !read(io_unit) pub_eda_dE_x_rep, pub_eda_dE_c_rep

       ! Construct the map of cumulative fragment NGWF index to atom identifier
       ! for this calculation.
       read(io_unit) num_cum_frag_ngwfs

       allocate( &
            atom_ids_by_cum_frag_ngwf(num_cum_frag_ngwfs), &
            stat = ierr, &
            errmsg = errmsg)

       call utils_alloc_check( &
            "eda_super_restart_read_metadata", &
            trim(errmsg), &
            ierr)

       read(io_unit) atom_ids_by_cum_frag_ngwf

       ! If not reading from an EDA_PREP calculation
       if (.not.(preptool_value)) then

          ! Frozen, antisymmetrized
          read(io_unit) complex_energy_frzIdem
          read(io_unit) pub_eda_dE_x_rep, pub_eda_rep_x, pub_eda_frz_x
          read(io_unit) pub_eda_dE_c_rep, pub_eda_rep_c, pub_eda_frz_c

          ! Polarisation:
          read(io_unit) complex_energy_pol_simul
          ! individual fragment polarisations
          if (pub_eda_frag_isol_pol) then
             do it=1,pub_frag_iatm(0)
                read(io_unit) complex_energy_pol(it)
             end do
          end if

          ! Charge transfer:
          ! individual fragment-pair delocalisations
          if (pub_eda_frag_isol_ct) then
             ! loop the fragment pairs given in the input file:
             do it=1,SIZE(pub_frag_deloc(1,:))
                frag1 = pub_frag_deloc(1,it)
                frag2 = pub_frag_deloc(2,it)
                read(io_unit) complex_energy_ct(frag1,frag2)
             end do
             read(io_unit) eda_cont_ctdeloc_it
          end if

       end if ! preptool

       ! Close the file
       close(io_unit,iostat=io_stat)
       call utils_close_unit_check('eda_super_restart_read_metadata',&
            filename,io_stat)

       write(stdout,'(a)')' done'

    end if ! pub_on_root
    call services_flush


    ! Broadcast data from root

    call comms_bcast(pub_root_proc_id, eda_cont_mode)
    call comms_bcast(pub_root_proc_id, eda_cont_fc1)
    call comms_bcast(pub_root_proc_id, eda_cont_fc2)
    call comms_bcast(pub_root_proc_id, pub_eda_mode)
    call comms_bcast(pub_root_proc_id, pub_frag_counter)
    call comms_bcast(pub_root_proc_id, pub_frag_counter2)

    do it=1,pub_frag_iatm(0)
       call comms_bcast(pub_root_proc_id, pub_eda_isol(it))
       call comms_bcast(pub_root_proc_id, pub_eda_isol_x(it))
       call comms_bcast(pub_root_proc_id, pub_eda_isol_c(it))
       call comms_bcast(pub_root_proc_id, pub_frag_data(it)%ngwf_basis(1)%num)
       call comms_bcast(pub_root_proc_id, pub_frag_data(it)%rep%n_occ)
    end do
    call comms_bcast(pub_root_proc_id, complex_energy_frz)
    call comms_bcast(pub_root_proc_id, pub_eda_frz_x)
    call comms_bcast(pub_root_proc_id, pub_eda_frz_c)
    !call comms_bcast(pub_root_proc_id, pub_eda_rep_x)
    !call comms_bcast(pub_root_proc_id, pub_eda_rep_c)
    call comms_bcast(pub_root_proc_id, pub_eda_intrafrag_x_frz)
    call comms_bcast(pub_root_proc_id, pub_eda_intrafrag_c_frz)
    !call comms_bcast(pub_root_proc_id, pub_eda_dE_x_rep)
    !call comms_bcast(pub_root_proc_id, pub_eda_dE_c_rep)

    call comms_bcast(pub_root_proc_id, num_cum_frag_ngwfs)

    if (.not. pub_on_root) then
       allocate( &
            atom_ids_by_cum_frag_ngwf(num_cum_frag_ngwfs), &
            stat = ierr, &
            errmsg = errmsg)

       call utils_alloc_check( &
            "eda_super_restart_read_metadata", &
            trim(errmsg), &
            ierr)
    end if

    call comms_bcast(pub_root_proc_id, atom_ids_by_cum_frag_ngwf)

    call comms_bcast(pub_root_proc_id, preptool_value)
    if (.not.(preptool_value)) then
       call comms_bcast(pub_root_proc_id, complex_energy_pol_simul)
       if (pub_eda_frag_isol_pol) then
          do it=1,pub_frag_iatm(0)
             call comms_bcast(pub_root_proc_id, complex_energy_pol(it))
          end do
       end if

       ! Charge transfer:
       ! individual fragment-pair delocalisations
       if (pub_eda_frag_isol_ct) then
          ! loop the fragment pairs given in the input file:
          do it=1,SIZE(pub_frag_deloc(1,:))
             frag1 = pub_frag_deloc(1,it)
             frag2 = pub_frag_deloc(2,it)
             call comms_bcast(pub_root_proc_id, complex_energy_ct(frag1,frag2))
          end do
          call comms_bcast(pub_root_proc_id, eda_cont_ctdeloc_it)
       end if
    end if

    call fragment_data_set_ngwf_index_map( &
         atom_ids_by_cum_frag_ngwf, &
         mdl, &
         ngwf_basis(1))

    deallocate( &
         atom_ids_by_cum_frag_ngwf, &
         stat = ierr, &
         errmsg = errmsg)

    call utils_dealloc_check( &
         "eda_super_restart_read_metadata", &
         trim(errmsg), &
         ierr)

    ! mjsp: Derive pub_eda_scfmi from pub_eda
    ! mjsp: (not included in the eda file for backwards-compatibility
    ! mjsp: reasons so must be calculated logically)
    pub_eda_scfmi = (pub_eda_mode .eq. EDA_FROZENIDEM .or. &
         pub_eda_mode .eq. EDA_POLFRAGLOC_DEVEL .or. &
         pub_eda_mode .eq. EDA_POLFRAGLOC_DEVEL_OFF .or. &
         pub_eda_mode .eq. EDA_POLSIMUL .or. &
         pub_eda_mode .eq. EDA_POLSIMUL_OFF .or. &
         pub_eda_mode .eq. EDA_CTFRAGLOC_DEVEL .or. &
         pub_eda_mode .eq. EDA_CTFRAGLOC_DEVEL_OFF)

  end subroutine eda_super_restart_read_metadata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine eda_super_restart_read_dkn_ngwfs(rep,denskern,ngwf_basis,mdl)
!
!    !======================================================================!
!    ! This subroutine reads in data from the 'EDA_PREP' task.              !
!    ! Reads in the supermolecule NGWFs and density kernel                  !
!    !----------------------------------------------------------------------!
!    ! Written by Max Phipps September 2015.                                !
!    !======================================================================!
!
!    use comms, only: pub_on_root
!    use restart, only: restart_kernel_read, restart_ngwfs_tightbox_input
!    use ngwf_representation, only: NGWF_REP
!    use kernel, only: DKERN
!    use function_basis, only: FUNC_BASIS
!    use model_type, only: MODEL
!
!    implicit none
!
!    ! Arguments (inputs)
!    type(NGWF_REP), intent(inout) :: rep
!    type(DKERN), intent(inout) :: denskern
!    type(FUNC_BASIS), intent(inout) :: ngwf_basis
!    type(MODEL), intent(inout) :: mdl
!
!    ! mjsp: Read in the kernel and NGWFs
!    call restart_kernel_read(denskern%kern)
!    call restart_ngwfs_tightbox_input(rep%ngwfs_on_grid, ngwf_basis, &
!         mdl%cell, mdl%fftbox, mdl%elements, &
!         'tightbox_'//trim(ngwf_basis%name)//'_supermolecule')
!
!  end subroutine eda_super_restart_read_dkn_ngwfs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !>
  !!  @brief Gets `dfdtau_fine` from @link
  !!  hamiltonian.hamiltonian_lhxc_calculate @endlink.
  !!
  !!  @param [inout] ham
  !!    Hamiltonian.
  !!
  !!  @param [in] mdl
  !!    Model.
  !!
  !!  @param [in] rep
  !!    NGWF representation.
  !!
  !!  @param [in] ngwf_basis
  !!    NGWF basis.
  !!
  !!  @param [in] kern
  !!    Density kernel.
  !!
  !!  @return
  !!    `dfdtau_fine`. If @link rundat.pub_xc_ke_density_required @endlink is
  !!    `false`, it will have a size of `1` in each dimension.
  !!
  !!  @internal
  !!    @author hc5n17
  !!    @since 2019
  !!  @endinternal
  !!
  function eda_get_dfdtau_fine(ham, mdl, rep, ngwf_basis, kern)
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_lhxc_calculate
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_num_spins, pub_xc_ke_density_required
    use sparse_embed, only: SPAM3_EMBED_ARRAY
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    type(NGWF_HAM), intent(inout) :: ham
    type(MODEL), intent(in) :: mdl
    type(NGWF_REP), intent(in) :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(mdl%nsub)
    type(SPAM3_EMBED_ARRAY), intent(in) :: kern
    real(kind = DP), allocatable :: eda_get_dfdtau_fine(:, :, :, :)

    real(kind = DP) :: energy
    real(kind = DP), allocatable :: lhxc_fine(:, :, :, :, :)

    integer :: ierr
    character(255) :: errmsg = ""

    if (allocated(eda_get_dfdtau_fine)) then
       deallocate( &
            eda_get_dfdtau_fine, &
            stat = ierr, &
            errmsg = errmsg)

       call utils_dealloc_check( &
            "eda_get_dfdtau_fine", &
            trim(errmsg), &
            ierr)
    end if

    allocate( &
         lhxc_fine( &
         mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, &
         pub_num_spins, &
         mdl%nsub), &
         stat = ierr, &
         errmsg = errmsg)

    call utils_alloc_check( &
         "eda_get_dfdtau_fine", &
         trim(errmsg), &
         ierr)

    if (pub_xc_ke_density_required) then
       allocate( &
            eda_get_dfdtau_fine( &
            mdl%fine_grid%ld1, &
            mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12, &
            pub_num_spins), &
            stat = ierr, &
            errmsg = errmsg)
    else
       allocate( &
            eda_get_dfdtau_fine(1, 1, 1, 1), &
            stat = ierr, &
            errmsg = errmsg)
    end if

    call utils_alloc_check( &
         "eda_get_dfdtau_fine", &
         trim(errmsg), &
         ierr)

    call hamiltonian_lhxc_calculate( &
         lhxc_fine, &
         energy, &
         ham%dijhat, &
         mdl, &
         rep%ngwfs_on_grid, &
         ngwf_basis, &
         kern, &
         rep%ngwf_overlap, &
         rep%sp_overlap, &
         dfdtau_fine = eda_get_dfdtau_fine)

    deallocate( &
         lhxc_fine, &
         stat = ierr, &
         errmsg = errmsg)

    call utils_dealloc_check( &
         "eda_get_dfdtau_fine", &
         trim(errmsg), &
         ierr)
  end function eda_get_dfdtau_fine

end module eda_driver_supermol

