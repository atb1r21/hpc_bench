! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                    Property analysis module                    !
!                                                                !
! This module performs property analysis using canonical         !
! molecular orbitals obtained from converged ONETEP calculations.!
!----------------------------------------------------------------!
! Written by Chris-Kriton Skylaris on 10/05/2006                 !
! Population analysis added by Peter Haynes, July 2006           !
! NGWF analysis added by Mark Robinson, January 2008             !
! General re-writing and updating by Nicholas Hine, 2009-2011.   !
! Further modifications by Laura Ratcliff, David O'Regan,        !
! Louis Lee, Robert Bell, Jacek Dziedzic and Simon Dubois.       !
!================================================================!

module properties

  use constants, only: DP
  use rundat, only: pub_debug_on_root

  implicit none

  private

  public :: properties_calculate
  public :: properties_polarisation
  public :: properties_plot_delta_density

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine properties_calculate(denskern, ham, &
       rep, ngwf_basis, proj_basis, nl_projectors, hub_proj_basis, hub, &
       core_basis, core_wvfns, mdl, hfxstate, &
       lhxc_fine, ham_type, cond_n_occ, cross_overlap_cj, cond_denskern, &
       dfdtau_fine)

    !======================================================================!
    ! This subroutine generates and diagonalises a Hamiltonian and         !
    ! performs various kinds of analysis with the resulting eigenvectors   !
    ! and eigenenergies.                                                   !
    !----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 10/05/2006.                      !
    ! Modified by Nicholas Hine in July 2009 to use function_basis and     !
    ! SPAM3 routines.                                                      !
    ! Modified by Nicholas Hine in October 2009 to make use of ScaLAPACK   !
    ! for diagonalisation of large systems.                                !
    ! Modified by Nicholas Hine in February 2010 to remove #ifdef's by     !
    ! moving code relating to eigensolving to dense_mod.                   !
    ! Modified by Nicholas Hine in October 2010 to use NGWF_REP.           !
    ! Modified by Laura Ratcliff in October 2010 for conduction            !
    ! calculations.                                                        !
    ! Simplified by Nicholas Hine in April 2011 to reduce branching of     !
    ! code in spectra calculations.                                        !
    ! Modified for embedding structures by Robert Charlton, 24/06/2018.    !
    !======================================================================!

    use augmentation, only: aug_projector_denskern
    use bandstructure, only: bandstructure_calculate, bsunfold_calculate
    use comms, only: pub_on_root, comms_barrier
    use constants, only: stdout, paw_en_size, max_spins, &
         REP_SWEX_PROPERTIES_DMA_1, REP_SWEX_PROPERTIES_DMA_2
    use datatypes, only: FUNCTIONS, data_set_to_zero, &
         data_functions_alloc, data_functions_dealloc
    use dense, only: DEM
    use dma, only: dma_calculate, dma_output_potential, dma_free_multipoles
    use dmft, only: hubbard_dmft_interface
    use eigenstates, only: eigenstates_calculate, eigenstates_create_storage, &
         eigenstates_destroy_storage
    use etrans, only: etrans_calculate
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_sum_ppd_funcs ! lpl: Used NBO NGWF analysis
    use hamiltonian, only: hamiltonian_dens_dep_matrices, &
         hamiltonian_build_matrix
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: DKERN, kernel_purify
    use model_type, only: MODEL
    use multipole_ops, only: SPHERICAL_MULTIPOLE_SET
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use npa, only: npa_main
    use population, only: population_analysis_lowdin, &
         population_analysis_mulliken, population_analysis_ddec
    use projectors, only: PROJECTOR_SET
    use restart, only: restart_kernel_write,restart_kernel_read
    use rundat, only: pub_homo_dens_plot, pub_lumo_dens_plot, pub_homo_plot, &
         pub_lumo_plot, pub_dos_smear, pub_edft, &
         pub_num_eigenvalues, pub_write_density_plot, pub_grd_format, &
         pub_cube_format, pub_dx_format, pub_popn_calculate, &
         pub_ngwf_analysis, pub_nnho, pub_do_bandstructure, &
         pub_polarisation_calculate, pub_popn_mulliken_partial, &
         pub_spread_calculate, pub_cond_calculate, pub_etrans_lcr,&
         pub_etrans_bulk, pub_paw,pub_plot_nbo, pub_write_nbo,&
         pub_nbo_pnao_analysis, pub_nbo_write_dipole, pub_ddec, pub_task, &
         pub_devel_code, pub_dma_calculate, pub_dma_output_potential_reference,&
         pub_dma_output_potential, pub_use_aux_ngwfs, &
         pub_lowdin_popn_calculate, pub_num_spins, PUB_1K, &
         pub_polarisation_simcell_calculate, pub_spin_fac, pub_num_kpoints, &
         pub_dmft_points, pub_dmft_sc, pub_dmft_nkpoints, pub_dmft_kernel, &
         pub_dmft_purify_sc, pub_efield_calculate, pub_pol_emb_pot, &
         pub_xc_ke_density_required, pub_polarisation_berry, &
         pub_bsunfld_calculate, pub_rootname, pub_dmft_read, &
         pub_inner_loop_iteration
    use sparse, only: SPAM3, sparse_create, sparse_destroy, &
         sparse_read, sparse_copy
    use sparse_array, only: SPAM3_ARRAY, sparse_array_destroy
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_copy, &
         SPAM3_EMBED_ARRAY, sparse_embed_array_create, &
         sparse_embed_array_copy, sparse_embed_array_destroy, &
         sparse_embed_array_scale, sparse_embed_array_to_sparse_array, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, &
         utils_devel_code, utils_assert, utils_unit, utils_banner

    implicit none

    ! Arguments
    type(NGWF_REP), intent(in) :: rep
    type(NGWF_HAM), intent(inout) :: ham   ! NGWF Hamiltonian type - changed to argument by lr408
    type(DKERN), intent(inout) :: denskern                ! denskern matrix
    type(MODEL), intent(inout)   :: mdl ! ja531-> made inout for PDOS in eigenstates
    type(HFX_STATE), intent(inout), target :: hfxstate
    type(FUNC_BASIS), intent(in) :: ngwf_basis(mdl%nsub)
    type(FUNC_BASIS), intent(in) :: proj_basis(mdl%nsub)
    type(FUNC_BASIS), intent(in) :: core_basis(mdl%nsub)
    type(FUNC_BASIS), intent(in) :: hub_proj_basis(mdl%nsub)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(mdl%nsub)
    type(PROJECTOR_SET), intent(inout) :: core_wvfns(mdl%nsub)
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub)
    ! lr408: Optional string defining which type of Hamiltonian this is for
    character(len=*), intent(in) :: ham_type
    ! lr408: Optional integer defining the total number of states to include for optical spectra
    integer, optional, intent(in) :: cond_n_occ(:,:)
    ! ndmh: Optional matrix arg for joint basis calculations, overlap of joint and cross NGWFs
    type(SPAM3_EMBED), optional, intent(in) :: cross_overlap_cj
    ! ndmh: Optional DKERN arg for joint basis calculations, giving joint denskern
    type(DKERN), optional, intent(in) :: cond_denskern
    ! JCW: Optional gradient of XC energy per unit volume wrt kinetic energy
    ! JCW: density. Required to evaluate lhxc matrix elements for meta-GGAs.
    real(kind=DP), optional, intent(inout) :: dfdtau_fine(:,:,:,:)
    ! JCW: In general, we have
    ! JCW:    dfdtau_fine(mdl%fine_grid%ld1,&
    ! JCW:    mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    ! JCW: or
    ! JCW:    dfdtau_fine(1,1,1,1)
    ! JCW: when it is unneeded.

    ! Local Variables
    type(FUNCTIONS) :: orth_on_grid(1:1) ! on-site orthogonal NGWFs
    real(kind=DP) :: total_energy      ! total energy
    real(kind=DP) :: lhxc_energy       ! pseudo+hartree+xc energy
    real(kind=DP) :: hubbard_energy    ! Hubbard energy
    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    type(SPAM3_EMBED), allocatable :: rhoij(:)
    type(SPAM3_EMBED_ARRAY) :: purkern         ! Purified kernel K
    integer :: ierr                    ! memory allocation error flag

    ! dense matrix for eigenvectors
    type(DEM), allocatable :: eigs_dens(:)
    ! hamiltonian eigenvalues
    real(kind=DP), allocatable :: eigen_en(:,:)
    ! denskern occupancies
    real(kind=DP), allocatable :: eigen_occ(:,:)

    ! lpl: NBO local variables
    type(SPAM3_EMBED) :: pnao_tr(1)          ! NBO NGWF analysis
    type(SPAM3_EMBED) :: dipole_mat(3)       ! NBO NGWF dipole matrix input
    integer :: ixyz, is

    logical :: read_ham1, dmft_sc_back, store_exists
    integer :: store_unit
    type(SPAM3_EMBED) :: store_overlap_mat

    ! jd: DMA multipoles
    type(SPHERICAL_MULTIPOLE_SET) :: dma_multipoles(3)

    ! agrecocmplx
    logical :: loc_cmplx
    ! jcap: region counters
    integer :: isub, jsub
    ! jcap: temporary kernel
    type(SPAM3_ARRAY) :: temp_kern
    type(SPAM3) :: ham_array(pub_num_spins), kern_array(pub_num_spins), &
         aug_array(pub_num_spins)

    ! ja531-> did we create a kernel workspace?
    logical :: made_kerws

    if(pub_debug_on_root) write(stdout,'(a)') 'Entering properties_calculate.'

    ! cks: start timer
    call timer_clock("properties_calculate",1)

    pub_inner_loop_iteration = -2

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine properties_calculate not ready yet for more&
         & than one k-point.')

    if (.not.pub_xc_ke_density_required) then
       ! JCW: Non-tau-dependent XC functional
       ! JCW: dfdtau_fine is either not present, or has a length of 1
       if ( present(dfdtau_fine) ) then
          call utils_assert( &
               size(dfdtau_fine) == 1,&
               "Error in properties_calculate: &
               &pub_xc_ke_density_required is false, but dfdtau_fine is present &
               &and has a size > 1.")
       end if
    else
       ! JCW: tau-dependent XC functional
       ! JCW: Check dfdtau_fine is present
       call utils_assert( present(dfdtau_fine),&
            "Error in properties_calculate: &
            &pub_xc_ke_density_required is true, but dfdtau_fine is not &
            &present.")
    end if

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(1)%iscmplx

    ! ebl: For DMFT, skipping hamiltonian calculation if store file exists
    !      or pub_dmft_read is not true
    read_ham1 = .false.
    if(pub_dmft_points>0 .and. pub_dmft_read)then
       inquire(file=trim(pub_rootname)//'.ham1',exist=read_ham1)
    endif

    if (.not.read_ham1) then

       ! JCW: Abort if tau-dependent XC functional requested, since this has not
       ! JCW: been implemented in properties_calculate (yet)
       if (pub_xc_ke_density_required) then
          if (pub_on_root) write(stdout,*) "WARNING: &
              &Support for tau-dependent XC functionals in &
              &properties_calculate is experimental."
       end if

       ! cks: === BUILD HAMILTONIAN & WRITE DENSITY/POTENTIAL ===================

       ! ndmh: calculate density dependent energies and matrices
       dmft_sc_back = pub_dmft_sc
       if ((.not.pub_cond_calculate) .and. (.not.(ham_type=='aux'))) then
          ! ebl: added writing of store files for DFT+DMFT
          pub_dmft_sc=.false.
          if(pub_dmft_purify_sc.and.pub_dmft_kernel/=0)then
             call sparse_embed_array_create(purkern,denskern%kern)
             call kernel_purify(purkern,denskern,rep%overlap,rep%inv_overlap,rep%n_occ)
             if (pub_num_spins == 1) call sparse_embed_array_scale(purkern,2.0_DP)
             ! JCW: Abort if a tau-dependent XC functional is requested
             ! JCW: (untested functionality)
             call utils_assert(.not.pub_xc_ke_density_required,"Error in &
                  &properties_calculate: the combination of tau-dependent &
                  &XC functionals and DMFT has not been tested.")
             call hamiltonian_dens_dep_matrices(ham, lhxc_fine, total_energy, &
                  lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
                  ngwf_basis, hub_proj_basis, hub, purkern, mdl, hfxstate, &
                  ham_update=.true., lhxc_fixed=.false.)
             if (pub_num_spins == 1) call sparse_embed_array_scale(purkern,0.5_DP)
             if (pub_dmft_points>0) then
                call utils_assert(.not. pub_pol_emb_pot, "Error in &
                     &properties_calculate: DMFT store functionality is &
                     &oblivious to polarisable embedding.")
                store_unit = utils_unit()
                open(unit=store_unit,file=trim(pub_rootname)//'.tot_energy')
                write(store_unit,*) total_energy,lhxc_energy,hubbard_energy
                close(store_unit)
             endif
             call hamiltonian_build_matrix(ham, rep)
             call sparse_embed_array_destroy(purkern)
          else
             call hamiltonian_dens_dep_matrices(ham, lhxc_fine, total_energy, &
                  lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
                  ngwf_basis, hub_proj_basis, hub, denskern%kern, mdl, &
                  hfxstate, ham_update=.true., lhxc_fixed=.false., &
                  dfdtau_fine = dfdtau_fine )
             if (pub_dmft_points>0) then
                call utils_assert(.not. pub_pol_emb_pot, "Error in &
                     &properties_calculate: DMFT store functionality is &
                     &oblivious to polarisable embedding.")
                store_unit = utils_unit()
                open(unit=store_unit,file=trim(pub_rootname)//'.tot_energy')
                write(store_unit,*) total_energy,lhxc_energy,hubbard_energy
                close(store_unit)
             endif

             ! ndmh: build Hamiltonian matrix from its component matrices
             call hamiltonian_build_matrix(ham, rep)

          end if
       end if
       pub_dmft_sc=dmft_sc_back

       ! ndmh: in PAW, rebuild rhoij as we need it again (very cheap)
       allocate(rhoij(pub_num_spins),stat=ierr)
       call utils_alloc_check('properties_calculate','rhoij',ierr)
       if (pub_paw) then
          do is=1,pub_num_spins
             rhoij(is)%structure = 'E'
             call sparse_embed_create(rhoij(is))
          end do
          call sparse_embed_extract_from_array(aug_array,rhoij)
          call sparse_embed_extract_from_array(kern_array,denskern%kern%m(:,PUB_1K))
          if (pub_cond_calculate) then
             !call aug_projector_denskern(rhoij,denskern,val_rep%sp_overlap)
          else
             call aug_projector_denskern(aug_array,kern_array,rep%sp_overlap%p)
          end if
          call sparse_embed_destroy_extracted_array(aug_array,rhoij,.true.)
          call sparse_embed_destroy_extracted_array(kern_array)
       end if

       if (pub_write_density_plot .and. (ham_type=='valence') .and. &
            (pub_cube_format .or. pub_grd_format .or. pub_dx_format)) then
          ! cks: write potential and density plot files
          call internal_plot_dens_and_pot
       endif

       ! cks: END BUILD HAMILTONIAN & WRITE DENSITY/POTENTIAL ===================

       ! ndmh: population analysis, NGWF analysis and polarisation analysis
       ! ndmh: all only work for valence Hamiltonian
       if (ham_type=='valence') then

          ! Allocate space for purified kernel (according to overlap sparsity)
          call sparse_embed_array_create(purkern, rep%overlap, n_spins=pub_num_spins, &
               n_kpoints=pub_num_kpoints)

          if (.not.pub_edft) then
             ! Calculate purified kernel
             call kernel_purify(purkern, denskern, rep%overlap, rep%inv_overlap, &
                                rep%n_occ)
          else
             ! ars: occupancies are fractional in EDFT
             call sparse_embed_array_copy(purkern, denskern%kern)
          end if

          ! pdh: POPULATION ANALYSIS
          if (pub_popn_calculate) call &
               population_analysis_mulliken(rep%overlap, purkern%m(:,PUB_1K), &
               mdl, ngwf_basis, partial_charge_orbital=-1)
          if (pub_lowdin_popn_calculate) call &
               population_analysis_lowdin(rep%overlap, purkern%m(:,PUB_1K), &
               mdl, ngwf_basis)
          if(pub_ddec) call &
               population_analysis_ddec(denskern%kern%m(:,PUB_1K), &
                    rep,ngwf_basis,mdl)
          ! pdh: END POPULATION ANALYSIS

          call sparse_embed_array_destroy(purkern)

          ! mr: NGWF ANALYSIS
          if ((pub_nnho) .and. (pub_ngwf_analysis .or. pub_spread_calculate)) then
             ! orthogonalise NNHO before characterising
             ! agrecocmplx
             ! rc2013: for the moment do "subsystem" orthogonalisations
             do isub=1,mdl%nsub
                call data_functions_alloc(orth_on_grid(1), &
                     ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
                if(mdl%nsub .gt. 1) then
                   call ngwfs_orthogonalise(orth_on_grid,rep%ngwfs_on_grid(isub), &
                        ngwf_basis(isub), rep%overlap%m(isub,isub), isub)
                else
                   call ngwfs_orthogonalise(orth_on_grid,rep%ngwfs_on_grid(isub), &
                        ngwf_basis(isub), rep%overlap%m(isub,isub))
                end if
                ! mr: print s/p/d character of NGWFs
                if (pub_ngwf_analysis) then
                   call properties_ngwfs_char(orth_on_grid(1), &
                        ngwf_basis(isub),mdl,3,ireg=isub)
                endif
                ! ddor: Print NGWF spreads
                ! rc2013: subsystem spread (whatever that means...)
                if (pub_spread_calculate) then
                   call properties_spread(orth_on_grid(1),ngwf_basis(isub),mdl, &
                        rep%inv_overlap%m(isub,isub),rep%overlap%m(isub,isub),.false.)
                endif
                call data_functions_dealloc(orth_on_grid(1))
             end do
          else
             ! mr: print s/p/d character of NGWFs
             if (pub_ngwf_analysis) then
                do isub=1,mdl%nsub
                   call properties_ngwfs_char(rep%ngwfs_on_grid(isub), &
                        ngwf_basis(isub),mdl,3,ireg=isub)
                end do
             endif
             ! ddor: Print NGWF spreads
             if (pub_spread_calculate) then
                do isub=1,mdl%nsub
                   call properties_spread(rep%ngwfs_on_grid(isub),ngwf_basis(isub),mdl, &
                        rep%inv_overlap%m(isub,isub),rep%overlap%m(isub,isub),.false.)
                end do
             endif
          end if
          ! mr: END NGWF ANALYSIS


          ! mr: POLARISATION
          ! ndmh: flag to activate polarisation calculations
          if (pub_polarisation_calculate) then

             ! smmd: get devel_code parameters
             if (pub_on_root.and.pub_polarisation_berry) &
                write(stdout,'(a,l1)') 'Properties devel_code : BERRY = ', &
                     pub_polarisation_berry

             ! lpl: Don't repeat calculation if NBO is going to call it later
             if(.not.(pub_write_nbo .and. pub_nbo_write_dipole)) then
                if (.not.pub_polarisation_berry) then
                   call properties_polarisation(rep,ngwf_basis,proj_basis, &
                        mdl, denskern%kern)
                endif
             endif

          end if
          ! mr: END POLARISATION

          if(pub_polarisation_simcell_calculate) then ! jd
             call properties_polarisation_simcell(rep%ngwfs_on_grid, ngwf_basis, &
                  mdl, denskern%kern, rep%ngwf_overlap)
          end if

       end if

       ! pdh: Only diagonalise Hamiltonian if required
       if ((pub_num_eigenvalues > 0) .or. (pub_dos_smear > 0.0_DP) .or. &
            (pub_homo_dens_plot >= 0) .or. (pub_lumo_dens_plot >= 0) .or. &
            (pub_homo_plot >= 0) .or. (pub_lumo_plot >= 0) .or. &
            pub_cond_calculate .or. pub_popn_mulliken_partial .or. &
            pub_dmft_points>0 ) then
          call eigenstates_create_storage(ngwf_basis, eigen_en, eigen_occ, &
               eigs_dens, loc_cmplx)

          if (ham_type/='joint') then

             call eigenstates_calculate(denskern%kern, ham, &
                  rep, ngwf_basis, mdl, hfxstate, ham_type, eigs_dens, &
                  eigen_en, eigen_occ, nl_projectors, proj_basis, core_wvfns, &
                  core_basis, cond_n_occ)

          else

             call eigenstates_calculate(denskern%kern, ham, &
                  rep, ngwf_basis, mdl, hfxstate, ham_type, eigs_dens, &
                  eigen_en, eigen_occ, nl_projectors, proj_basis, core_wvfns, &
                  core_basis, cond_n_occ, &
                  cross_overlap_cj, cond_denskern%kern%m(:,PUB_1K))

          end if

          if (pub_polarisation_calculate.and.pub_polarisation_berry) then
              call properties_polarisation_berry(rep,ngwf_basis,proj_basis, &
                   mdl, nl_projectors, denskern%kern)

          end if

          ! keeping eigen_en if required for DMFT interface
          if (pub_dmft_points == 0) then
             call eigenstates_destroy_storage(eigen_en, eigen_occ, &
                  eigs_dens)
          end if

       end if

       ! rc2013: bandstructure calculations not compatible with embedding
       ! pdh: BANDSTRUCTURE
       if (( (pub_do_bandstructure.and.(ham_type/='cond')) .or.  abs(pub_dmft_nkpoints)>1 ) &
            .and. (mdl%nsub == 1)) then
          call bandstructure_calculate(ham,rep,ngwf_basis(1),proj_basis(1), &
               nl_projectors(1),rhoij,mdl,ham_type)
       end if
       ! pdh: BANDSTRUCTURE

       ! gcc32: BANDSTRUCTURE UNFOLDING
       if (pub_bsunfld_calculate.and.(ham_type/='cond').and.(mdl%nsub == 1)) then
          call bsunfold_calculate(ham,rep,ngwf_basis(1),proj_basis(1), &
               nl_projectors(1),rhoij,mdl)
       end if
       ! gcc32: BANDSTRUCTURE UNFOLDING
    end if

    ! ebl: setting up files for TOSCAM DMFT calculation
    if (pub_dmft_points > 0 .and. mdl%nsub == 1) then

       call sparse_embed_extract_from_array(ham_array,ham%ham)
       call hubbard_dmft_interface(eigen_en,rep%n_occ(:,PUB_1K), &
            ham_array,rep%overlap,rep%inv_overlap,ngwf_basis, &
            hub_proj_basis(1),hub,rep,mdl%regions(1)%elements,denskern,mdl, &
            nl_projectors,proj_basis=proj_basis)
       call sparse_embed_destroy_extracted_array(ham_array,ham%ham,.true.)

       ! ebl: reading/writing store_overlap and store_inv_overlap
       call comms_barrier
       call sparse_embed_create(store_overlap_mat, rep%overlap, iscmplx=.false.)
       inquire(file=trim(pub_rootname)//'.overlap', exist=store_exists)
       if (store_exists .and. pub_dmft_read) then
          if (pub_on_root) write(stdout,'(3a)') "Loading ", trim(pub_rootname),&
               ".overlap"
          call sparse_read(store_overlap_mat%m(1,1), trim(pub_rootname)//'.overlap')
       else
          if (pub_on_root) write(stdout,'(3a)') "Not loading ", &
               trim(pub_rootname), ".overlap"
          call sparse_embed_copy(store_overlap_mat, rep%overlap)
       end if
       call comms_barrier

       ! ebl: doing away with the nopurify argument by handing denskern to
       !      population_analysis_mulliken in place of purkern
       if (.not. read_ham1) then
          call population_analysis_mulliken(store_overlap_mat, &
               denskern%kern%m(:,PUB_1K), mdl,ngwf_basis,partial_charge_orbital=-1)
       end if
       call comms_barrier
    else if (pub_dmft_points > 0 .and. mdl%nsub .gt. 1) then
       ! rc2013: not compatible with subsystem calculations yet
       write(stdout,'(/a)') utils_banner('=', &
            'WARNING: Dynamical mean field theory requested with embedding, &
            &but this combination is not implemented/tested. Skipping DMFT.')
    end if

    ! lpl: NBO
    if (pub_write_nbo .and. mdl%nsub == 1) then
       ! lpl: NGWF dipole matrix for NBO
       if (pub_nbo_write_dipole) then
          do ixyz=1,3
             call sparse_embed_create(dipole_mat(ixyz),rep%overlap)
          end do

          call properties_polarisation(rep,ngwf_basis,proj_basis,mdl, &
               denskern%kern, dipole_mat=dipole_mat)
       end if

       ! lpl: PNAO analysis + NBO
       if(pub_nbo_pnao_analysis) then
          pnao_tr(1)%structure = 'D'
          if(pub_use_aux_ngwfs) pnao_tr(1)%structure = 'Da' ! jhl52
          call sparse_embed_create(pnao_tr(1))
          if(pub_nbo_write_dipole) then
             call npa_main(denskern, rep, ham, mdl%elements, &
                  ngwf_basis(1), dipole_mat, pnao_tr(1))
          else
             call npa_main(denskern, rep, ham, mdl%elements, &
                  ngwf_basis(1), pnao_tr_output=pnao_tr(1))
          end if

          do isub=1,mdl%nsub
             call data_functions_alloc(orth_on_grid(1), &
                  ngwf_basis(isub)%size_on_grid,rep%ngwfs_on_grid(isub)%iscmplx)

             call data_set_to_zero(orth_on_grid(1))
             call function_ops_sum_ppd_funcs(orth_on_grid,ngwf_basis(isub), &
                  pnao_tr(1)%m(isub,isub),1,1,pnao_tr(1)%m(isub,isub), &
                  rep%ngwfs_on_grid(isub),ngwf_basis(isub))

             ! lpl: NGWF analysis
             if(pub_on_root) write(stdout,'(a)') &
                  '============ PNAO s/p/d/f Character ============'
             call properties_ngwfs_char(orth_on_grid(1),ngwf_basis(isub),mdl,3)
             if(pub_on_root) write(stdout,'(a)') &
                  '================================================'
             call data_functions_dealloc(orth_on_grid(1))
          end do

          call sparse_embed_destroy(pnao_tr(1))
       else ! No PNAO analysis
          if(pub_nbo_write_dipole) then
             call npa_main(denskern, rep, ham, mdl%elements, &
                  ngwf_basis(1), dipole_mat)
          else
             call npa_main(denskern, rep, ham, mdl%elements, ngwf_basis(1))
          end if
       end if

       ! lpl: Deallocate NGWF dipole matrix
       if(pub_nbo_write_dipole) then
          do ixyz=3,1,-1
             call sparse_embed_destroy(dipole_mat(ixyz))
          end do
       end if
    else if (pub_write_nbo .and. mdl%nsub .gt. 1) then
       ! rc2013: not compatible with subsystem calculations yet
       write(stdout,'(/a)') utils_banner('=', &
            'WARNING: NBO analysis requested with embedding, but this &
            &combination is not implemented/tested. Skipping NBO analysis.')
    end if
    if(pub_plot_nbo) call internal_plot_nbo
    ! lpl: NBO

    ! smmd: TRANSMISSION FUNCTION
    if ((pub_etrans_lcr .or. pub_etrans_bulk) .and. (mdl%nsub == 1)) then
       if (.not.pub_cond_calculate .and. ham_type == 'valence') then
          call etrans_calculate(ham, ham_type, rep, denskern%kern, ngwf_basis, mdl)

       elseif ((pub_task == 'COND' .or. pub_task == 'PROPERTIES_COND') &
            .and. ham_type == 'joint') then
          ! rab207: conduction transport calculations currently not working
          ! can override the removal of this functionality using the devel code
          if (utils_devel_code(.false.,'ETRANS','USE_COND',pub_devel_code)) then
             call etrans_calculate(ham, ham_type, rep, denskern%kern, ngwf_basis, mdl)
          else
             call utils_abort("Transmission using conduction NGWFs is still under development.&
                  & To use this functionality, use the devel code 'devel_code ETRANS:USE_COND=T:ETRANS'")
          endif
       endif
    else if ((pub_etrans_lcr .or. pub_etrans_bulk) .and. (mdl%nsub .gt. 1)) then
       ! rc2013: not compatible with subsystem calculations yet
       write(stdout,'(/a)') utils_banner('=', &
            'WARNING: electronic transport requested with embedding, but this &
            &combination is not implemented/tested. &
            &Skipping Electronic Transport analysis.')
    end if
    ! smmd: TRANSMISSION FUNCTION

    ! jd: DISTRIBUTED MULTIPOLE ANALYSIS
    if(pub_dma_calculate) then
       ! jcap: extract the denskern into a SPAM3_ARRAY
       call sparse_embed_array_to_sparse_array(temp_kern,denskern%kern,1)
       call dma_calculate(dma_multipoles, 'properties_dma', &
            rep%swexes(REP_SWEX_PROPERTIES_DMA_1:REP_SWEX_PROPERTIES_DMA_2), &
            rep%n_occ, rep%overlap, ngwf_basis, mdl, REP_SWEX_PROPERTIES_DMA_1,&
            'P', denskern%kern, dkn_scale_factor = pub_spin_fac*0.5_DP)
       if (pub_dma_output_potential) then
          call dma_output_potential(dma_multipoles, mdl, &
               pub_dma_output_potential_reference, ngwf_basis(1), &
               temp_kern, rep%ngwf_overlap%m(1,1), rep%ngwfs_on_grid(1))
       end if
       call sparse_array_destroy(temp_kern)
       call dma_free_multipoles(dma_multipoles)
    end if
    ! jd: DISTRIBUTED MULTIPOLE ANALYSIS

    if (.not. read_ham1) then
       if (pub_paw) then
          do is=pub_num_spins,1,-1
             call sparse_embed_destroy(rhoij(is))
          end do
       end if
       deallocate(rhoij,stat=ierr)
       call utils_dealloc_check('properties_calculate','rhoij',ierr)
    end if

    ! cks: stop timer
    call timer_clock("properties_calculate",2)

    if(pub_debug_on_root) write(stdout,'(a)') 'Leaving properties_calculate.'

    contains

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine internal_plot_dens_and_pot
      ! Trivially modified by Jacek Dziedzic on 14/05/2010 to delegate   !
      ! the unit conversion to visual_scalarfield.                       !
      ! ndmh: added PAW compensation density, and changed variable names !
      ! to match other parts of the code 29/09/2010.                     !
      ! cks: added output of electrostatic potential from electron       !
      ! density + Gaussian smeared ions under open BC. 16/7/2014/        !
      ! Modified for embedding by Joseph Prentice, August 2018           !
      ! rjc: added kinetic energy density output and electron            !
      ! localisation descriptor 07/2019

        use augmentation, only: augmentation_density_on_grid
        use comms, only: comms_bcast, pub_my_proc_id
        use constants, only: DP, ANGSTROM, HARTREE_IN_EVS, UP, DN, CRLF, &
             HARTREE_PER_BOHR_TO_MV_PER_CM
        use cutoff_coulomb, only: cutoff_coulomb_hartree
        use density, only: density_on_grid
        use eld, only: eld_on_grid
        use hartree, only: hartree_on_grid, hartree_via_multigrid
        use is_smeared_ions, only: smeared_ion_hartree, smeared_ion_density
        use kernel, only: kernel_validate_ksk, kernel_workspace_create,&
             kernel_workspace_destroy
        use ke_density, only: ke_density_on_grid
        use paw, only: paw_sphere_density_on_grid
        use rundat, only: pub_aug, pub_aug_den_dim, pub_is_smeared_ion_rep, &
             pub_spin_fac, pub_is_implicit_solvent, pub_finite_difference_order, &
             pub_edft_smearing_width, pub_pbc_smeared_ion_rep, &
             pub_chemical_softness, pub_print_qc, pub_coulomb_cutoff, &
             pub_active_region, pub_eld_calculate, &
             pub_eld_function, pub_ke_density_calculate
        use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
             sparse_embed_create, sparse_embed_destroy, &
             sparse_embed_array_scale, sparse_embed_array_create, &
             sparse_embed_array_copy, sparse_embed_array_axpy, &
             sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
        use sparse, only: SPAM3, sparse_create, sparse_destroy
        use utils, only: utils_banner, utils_qc_print, utils_assert
        use visual, only: visual_scalarfield, visual_efield_on_grid_obc, &
             visual_efield_on_grid_pbc

        implicit none

        ! Local Variables
        real(kind=DP), allocatable :: density_fine(:,:,:,:)! density
        real(kind=DP), allocatable :: density_plot(:,:,:,:)! density
        real(kind=DP), allocatable :: ke_density_fine(:,:,:,:) !ke density
        real(kind=DP), allocatable :: eld_fine(:,:,:,:)! eld
        real(kind=DP), allocatable :: chemsoft_fine(:,:,:,:)! chemical softness
        real(kind=DP), allocatable :: chemsoft_plot(:,:,:,:)! chemical softness
        real(kind=DP), allocatable :: efield_fine(:,:,:,:) ! efield
        real(kind=DP), allocatable :: nhat_den_grad(:,:,:,:,:)! compensation den
        real(kind=DP), allocatable :: nhat_chemsoft_grad(:,:,:,:,:)! compensation chemsoft
        real(kind=DP), allocatable :: lhxc_si_fine(:,:,:,:)! for elec pot output

        integer :: ierr          ! memory allocation error flag
        integer :: fine_max_slabs12, fine_ld1, fine_ld2 ! jd: shortcut notation
        character(len=3) :: spin_string
        integer :: qc_npoints
        real(kind=dp), dimension(:), allocatable :: qc_point
        character(len=4) :: qc_point_char
        integer :: qc_owner_proc, qc_first_on_grid, ii
        integer, dimension(:,:), allocatable :: qc_coords
        integer       :: isub

        type(SPAM3_EMBED_ARRAY) :: KmKSK

        real(kind=DP), allocatable :: rho_tot(:,:,:,:)
        type(SPAM3) :: kern_array(pub_num_spins)  ! temporary kernel array
        real(kind=DP), allocatable :: reg_density_fine(:,:,:,:) ! density
        real(kind=DP), allocatable :: sub_nhat_den_grad(:,:,:,:,:) ! compensation den
        real(kind=DP), allocatable :: sub_nhat_chemsoft_grad(:,:,:,:,:) ! compensation den

        ! ---------------------------------------------------------------------

        if (pub_on_root) then
           if(pub_efield_calculate) then
              write(stdout,'(/a)') utils_banner('=', &
                   'Writing density, potential and electric field')
           else
              write(stdout,'(/a)') utils_banner('=', &
                   'Writing density and potential')
           end if
        end if

        ! cks: Allocate charge density slabs
        allocate(density_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
             mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
        call utils_alloc_check('internal_plot_dens_and_pot &
             &(properties_calculate)','density_fine',ierr)
        allocate(reg_density_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
             mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
        call utils_alloc_check('internal_plot_dens_and_pot &
             &(properties_calculate)','reg_density_fine',ierr)

        ! cks:########## DENSITY OUTPUT ##########

        ! ndmh: Scale denskern for spin degeneracy
        ! jme: KPOINTS_DANGER: check this factor for more than 1 k-point
        if (pub_num_spins == 1) then
           call sparse_embed_array_scale(denskern%kern, pub_spin_fac)
        end if


        ! jcap: loop over regions
        density_fine=0.0_DP
        do isub=1,mdl%nsub
           do jsub=1,mdl%nsub
              reg_density_fine=0.0_DP
              ! jcap: create temporary density kernel
              call sparse_embed_extract_from_array(kern_array,&
                   denskern%kern%m(:,PUB_1K),isub,jsub)

              ! cks: Calculate data-parallelised charge density
              call density_on_grid(reg_density_fine, mdl%fine_grid, mdl%dbl_grid, &
                   mdl%cell, mdl%fftbox, kern_array, rep%ngwf_overlap%m(isub,jsub), &
                   rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                   rep%ngwfs_on_grid(jsub), ngwf_basis(jsub))

              density_fine=density_fine+reg_density_fine

              ! jcap: destroy temporary density kernel
              call sparse_embed_destroy_extracted_array(kern_array)
           end do
        end do

        ! ndmh: Calculate compensation density on fine grid
        if (pub_aug) then ! PAW version

           fine_ld1 = mdl%fine_grid%ld1
           fine_ld2 = mdl%fine_grid%ld2
           fine_max_slabs12 = mdl%fine_grid%max_slabs12
           allocate(nhat_den_grad(fine_ld1,fine_ld2,fine_max_slabs12,&
                pub_num_spins,0:pub_aug_den_dim),stat=ierr)
           call utils_alloc_check('internal_plot_dens_and_pot &
                &(properties_calculate)','nhat_den_grad',ierr)
           allocate(sub_nhat_den_grad(fine_ld1,fine_ld2,fine_max_slabs12,&
                pub_num_spins,0:pub_aug_den_dim),stat=ierr)
           call utils_alloc_check('internal_plot_dens_and_pot &
                &(properties_calculate)','sub_nhat_den_grad',ierr)
           allocate(density_plot(fine_ld1,fine_ld2,fine_max_slabs12,&
                pub_num_spins),stat=ierr)
           call utils_alloc_check('internal_plot_dens_and_pot &
                &(properties_calculate)','density_plot',ierr)

           ! CREATE PSEUDO DENSITY (INCLUDING NHAT) AND OUTPUT

           nhat_den_grad = 0.0_DP
           ! jcap: loop over regions
           do isub=1,mdl%nsub
              do jsub=1,mdl%nsub
                 sub_nhat_den_grad = 0.0_DP
                 call sparse_embed_extract_from_array(kern_array,&
                      denskern%kern%m(:,PUB_1K),isub,jsub)
                 call augmentation_density_on_grid(sub_nhat_den_grad, mdl%fine_grid, &
                      mdl%cell,mdl%regions(isub)%pseudo_sp,mdl%regions(jsub)%paw_sp,&
                      mdl%aug_box,kern_array, rep%sp_overlap%m(isub,jsub))
                 nhat_den_grad=nhat_den_grad+sub_nhat_den_grad
                 call sparse_embed_destroy_extracted_array(kern_array)
              end do
           end do

           density_plot = density_fine + nhat_den_grad(:,:,:,:,0)

           ! ndmh: spin polarisation: calculate total density and spin density
           if (pub_num_spins == 2) then
              density_plot(:,:,:,1) = density_plot(:,:,:,UP) + &
                   density_plot(:,:,:,DN)
              density_plot(:,:,:,2) = density_plot(:,:,:,UP) - &
                   2.0_DP * density_plot(:,:,:,DN)
           end if

           ! cks: output density in plot format file
           ! vm: output density in Angstrom rather than in Bohr
           ! jd: ... but leave the unit conversion to visual_scalarfield
           call visual_scalarfield( &
                density_plot(:,:,:,1), mdl%fine_grid, mdl%cell, &
                'Electronic pseudo density (in e/ang^3) for:', '_ps_density', &
                mdl%elements, ANGSTROM**3)
           if (pub_num_spins == 2) call visual_scalarfield( &
                density_plot(:,:,:,2), mdl%fine_grid, mdl%cell, &
                'Electronic pseudo spin density (in e/ang^3) for:', &
                '_ps_spindensity', mdl%elements, ANGSTROM**3)


           ! CREATE AE DENSITY AND OUTPUT

           nhat_den_grad(:,:,:,:,0) = 0.0_DP
           ! rc2013: copy rhoij across to avoid temporary arrays
           call sparse_embed_extract_from_array(kern_array,rhoij)
           call paw_sphere_density_on_grid(nhat_den_grad(:,:,:,:,0), &
                mdl%fine_grid,mdl%cell,mdl%aug_box, &
                kern_array,1.0_DP,-1.0_DP,mdl%regions(1)%paw_sp)
           call sparse_embed_destroy_extracted_array(kern_array)
           density_plot = density_fine + nhat_den_grad(:,:,:,:,0)

           ! ndmh: spin polarisation: calculate total density and spin density
           if (pub_num_spins == 2) then
              density_plot(:,:,:,1) = density_plot(:,:,:,UP) + &
                   density_plot(:,:,:,DN)
              density_plot(:,:,:,2) = density_plot(:,:,:,UP) - &
                   2.0_DP * density_plot(:,:,:,DN)
           end if

           ! cks: output density in plot format file
           ! vm: output density in Angstrom rather than in Bohr
           ! jd: ... but leave the unit conversion to visual_scalarfield
           call visual_scalarfield( &
                density_plot(:,:,:,1), mdl%fine_grid, mdl%cell, &
                'All-electron density (in e/ang^3) for:', '_ae_density', &
                mdl%elements, ANGSTROM**3)
           if (pub_num_spins == 2) call visual_scalarfield( &
                density_plot(:,:,:,2), mdl%fine_grid, mdl%cell, &
                'All-electron spin density (in e/ang^3) for:', &
                '_ae_spindensity', mdl%elements, ANGSTROM**3)

        else ! non-PAW version

           ! ndmh: spin polarisation: calculate total density and spin density
           if (pub_num_spins == 2) then
              density_fine(:,:,:,1) = density_fine(:,:,:,UP) + &
                   density_fine(:,:,:,DN)
              density_fine(:,:,:,2) = density_fine(:,:,:,UP) - &
                   2.0_DP * density_fine(:,:,:,DN)
           end if

           ! cks: output density in plot format file
           ! vm: output density in Angstrom rather than in Bohr
           ! jd: ... but leave the unit conversion to visual_scalarfield
           call visual_scalarfield( &
                density_fine(:,:,:,1), mdl%fine_grid, mdl%cell, &
                'Electronic density (in e/ang^3) for:', '_density', mdl%elements, &
                ANGSTROM**3)
           if (pub_num_spins == 2) call visual_scalarfield( &
                density_fine(:,:,:,2), mdl%fine_grid, mdl%cell, &
                'Electronic spin density (in e/ang^3) for:', '_spindensity', &
                mdl%elements, ANGSTROM**3)

        end if

        ! cks:####### END DENSITY OUTPUT ##########


        ! rjc:####### KE DENSITY & ELD OUTPUT ##########

        ! rjc: Nested conditionals.
        ! rjc: Outer if: calculates and outputs kinetic energy density to plot
        ! rjc: format file, if keywords KE_DENSITY_CALCULATE or ELD_CALCULATE
        ! rjc: chosen.
        ! rjc: Inner if: calculates electron localisation descriptors using
        ! rjc: kinetic energy density if only ELD_CALCULATE chosen, and
        ! rjc: outputs to plot format file

        if (pub_ke_density_calculate .or. pub_eld_calculate) then

           allocate(ke_density_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
                mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
           call utils_alloc_check('internal_plot_dens_and_pot &
                &(properties_calculate)','ke_density_fine',ierr)

           ! rjc: abort if more than one subsystem
           call utils_assert(mdl%nsub.eq.1, 'Error in properties_calculate: &
                &kinetic energy density and electron localisation descriptors &
                &not yet ready for more than one subsystem')

           ! rjc: Set isub = 1 as we are not doing embedding yet
           isub = 1

           ! jcap: Need to use a temporary array to store denskern in
           call sparse_embed_extract_from_array(kern_array,&
                denskern%kern%m(:,PUB_1K),isub,isub)

           call ke_density_on_grid(ke_density_fine, &
                mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, &
                kern_array, rep%ngwf_overlap%p, &
                rep%ngwfs_on_grid(isub), ngwf_basis(isub))

           ! jcap: Destroy temporary kernel
           call sparse_embed_destroy_extracted_array(kern_array)

           ! rjc: Calculate total and spin density from individual spins
           if (pub_num_spins == 2) then
              ke_density_fine(:,:,:,1) = ke_density_fine(:,:,:,UP) + &
                   ke_density_fine(:,:,:,DN)
              ke_density_fine(:,:,:,2) = ke_density_fine(:,:,:,UP) - &
                   2.0_DP * ke_density_fine(:,:,:,DN)
           end if

           ! rjc: Output kinetic energy density in plot format file
           ! rjc: output kinetic energy density in Angstrom rather than in Bohr
           ! rjc: but leave the unit conversion to visual_scalarfield
           if (pub_ke_density_calculate) then
           call visual_scalarfield( &
                ke_density_fine(:,:,:,1), mdl%fine_grid, mdl%cell, &
                'Kinetic energy density for:', '_ke_density', mdl%elements, &
                ANGSTROM**3)
           if (pub_num_spins == 2) call visual_scalarfield( &
                ke_density_fine(:,:,:,2), mdl%fine_grid, mdl%cell, &
                'Kinetic energy spin density for:', '_ke_spindensity', &
                mdl%elements, ANGSTROM**3)
           end if

           if (pub_eld_calculate) then

              allocate(eld_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
                   mdl%fine_grid%max_slabs12, pub_num_spins),stat=ierr)
              call utils_alloc_check('internal_plot_dens_and_pot &
                   &(properties_calculate)','eld_fine',ierr)

              ! rjc: Use kinetic energy density and charge density called above

              ! rjc: Apply spin polarisation to ke density and charge density
              ! rjc: to obtain individual spins from total and spin densities
              ! rjc: for ELD
              if (pub_num_spins == 2) then
                 ke_density_fine(:,:,:,1) = 0.5_DP*(ke_density_fine(:,:,:,1) + &
                      ke_density_fine(:,:,:,2))
                 ke_density_fine(:,:,:,2) = ke_density_fine(:,:,:,1) - &
                      ke_density_fine(:,:,:,2)

                 density_fine(:,:,:,1) = 0.5_DP*(density_fine(:,:,:,1) + &
                      density_fine(:,:,:,2))
                 density_fine(:,:,:,2) = density_fine(:,:,:,1) - &
                      density_fine(:,:,:,2)
              end if

              ! rjc: For non spin-polarised case, halve the total density
              ! rjc: such that the output is of the same magnitude as the
              ! rjc: individual up and down spins.
              ! rjc: Do the same for ke density to be consistent
              if (pub_num_spins == 1) then
                 ke_density_fine(:,:,:,1) = 0.5_DP*(ke_density_fine(:,:,:,1))
                 density_fine(:,:,:,1) = 0.5_DP*(density_fine(:,:,:,1))
              end if

              ! rjc: Calculate electron localisation descriptor on fine grid

              call eld_on_grid(eld_fine, mdl%fine_grid, density_fine, &
                   ke_density_fine, pub_eld_function)

              ! rjc: QC tests [test73]. Applied to spin and non-spin case of H2O
              if (pub_print_qc) then
                 qc_npoints=7
                 allocate(qc_coords(3,qc_npoints),stat=ierr)
                 call utils_alloc_check('internal_plot_dens_and_pot &
                      &(properties_calculate)','qc_point',ierr)
                 allocate(qc_point(qc_npoints),stat=ierr)
                 call utils_alloc_check('internal_plot_dens_and_pot &
                      &(properties_calculate)','qc_point',ierr)

                 ! rjc: Selection of points across grid  to sample in ELF,
                 ! rjc: chosen for a range of values between 0 and 1
                 qc_coords(:,1) = [63,63,63]
                 qc_coords(:,2) = [63,63,60]
                 qc_coords(:,3) = [64,64,60]
                 qc_coords(:,4) = [65,65,59]
                 qc_coords(:,5) = [66,65,50]
                 qc_coords(:,6) = [63,63,70]
                 qc_coords(:,7) = [70,91,91]

                 do is=1,pub_num_spins
                    do ii=1,qc_npoints

                       qc_owner_proc = mdl%fine_grid%proc_slab12(&
                            qc_coords(3,ii))

                       qc_first_on_grid = mdl%fine_grid%first_slab12(&
                            qc_owner_proc)

                       qc_point(ii)=0.0_dp

                       if(qc_owner_proc==pub_my_proc_id) then
                          qc_point(ii)=eld_fine(qc_coords(1,ii), &
                               qc_coords(2,ii), qc_coords(3,ii) - &
                               qc_first_on_grid+1,is)
                       end if

                       call comms_bcast(qc_owner_proc,qc_point(ii))

                       write(qc_point_char,'(i0)') ii
                       write(spin_string,'(i1)') is

                       if(pub_on_root) then
                          call utils_qc_print('props_ELD_spin'//&
                               trim(adjustl(spin_string))//'_n'//&
                               trim(adjustl(qc_point_char)),qc_point(ii))
                       end if
                    end do
                 end do

                 deallocate(qc_coords,stat=ierr)
                 call utils_dealloc_check('internal_plot_dens_and_pot &
                      &(properties_calculate)','qc_coords',ierr)
                 deallocate(qc_point,stat=ierr)
                 call utils_dealloc_check('internal_plot_dens_and_pot &
                      &(properties_calculate)','qc_point',ierr)

              end if

              ! rjc: output ELD in plot format file
              ! rjc: output ELD in unit volume, already consistent with
              ! rjc: calculation
              if (pub_num_spins == 1) then
                 call visual_scalarfield( &
                      eld_fine(:,:,:,1), mdl%fine_grid, mdl%cell, &
                      trim(pub_eld_function)//' for:', '_'//&
                      trim(pub_eld_function), mdl%elements, 1.0_DP )

              else if (pub_num_spins == 2) then
                 call visual_scalarfield(eld_fine(:,:,:,1), mdl%fine_grid,&
                      mdl%cell, trim(pub_eld_function)//' for spin UP:', '_'//&
                      trim(pub_eld_function)//'_UP', mdl%elements,&
                      1.0_DP)
                 call visual_scalarfield(eld_fine(:,:,:,2), mdl%fine_grid,&
                      mdl%cell, trim(pub_eld_function)//' for spin DOWN:','_'//&
                      trim(pub_eld_function)//'_DN', mdl%elements,&
                      1.0_DP )
              end if

              ! rjc: Deallocate memory for electron localisation descriptor
              deallocate(eld_fine,stat=ierr)
              call utils_dealloc_check('internal_plot_dens_and_pot &
                   &(properties_calculate)', 'eld_fine',ierr)

              ! rjc: Revert spin polarisation. Return to total charge density
              ! rjc: and spin density for later use in the subroutine
              if (pub_num_spins == 1) then
                 density_fine(:,:,:,1) = 2.0_DP*(density_fine(:,:,:,1))
              end if

              if (pub_num_spins == 2) then
                 density_fine(:,:,:,1) = density_fine(:,:,:,UP) + &
                      density_fine(:,:,:,DN)
                 density_fine(:,:,:,2) = density_fine(:,:,:,UP) - &
                      2.0_DP * density_fine(:,:,:,DN)
              end if

           end if !eld_calculate loop

           ! rjc: Deallocate memory for electron localisation descriptor
           deallocate(ke_density_fine,stat=ierr)
           call utils_dealloc_check('internal_plot_dens_and_pot &
                &(properties_calculate)','ke_density_fine',ierr)

        end if !eld_calculate or ke_density loop

        ! rjc:############## END KE DENSITY & ELD OUTPUT #############


        ! cks:########## POTENTIAL OUTPUT ##########

        ! cks: output potential in plot format file
        ! vm: output potential in eV rather than in Hartree
        !     also use "-" to do ESP: potential seen by a positive probe
        ! jd: ... but leave the unit conversion to visual_scalarfield

        ! pdh: only one potential even if spin-polarised
        ! rc2013: if this is an embedding run only print out the active potential
        call visual_scalarfield( &
             lhxc_fine(:,:,:,1,pub_active_region), mdl%fine_grid, mdl%cell, &
             'Local potential (eV) (Ion+Hartree+XC) for:', &
             '_potential', mdl%elements, -HARTREE_IN_EVS)


        ! cks: compute and output "proper" electrostatic potential, i.e. the
        !      local potential without the XC term.
        !      electron density

        call properties_electrostatic_potential(mdl, denskern, rep, ngwf_basis,&
             reg_density_fine, sub_nhat_den_grad)

        ! cks:###### END POTENTIAL OUTPUT ##########


        ! ndmh: restore normalisation of denskern
        ! jme: KPOINTS_DANGER: check this factor for more than 1 k-point
        if (pub_num_spins == 1) then
           call sparse_embed_array_scale(denskern%kern, 1.0_DP/pub_spin_fac)
        end if



        ! ja531-> Chemical softness
        ! ja531-> This is an experiment to see if we can get (local) Chemical softness from the
        ! derivative wrt chemical potiental... seems to work... might be useful as a descriptor
        ! needs testing.

        pub_chemical_softness = utils_devel_code(.false.,'PROPERTIES','CHEMICAL_SOFTNESS',pub_devel_code)

        if(pub_chemical_softness) then
           if (pub_on_root) then
              if(pub_efield_calculate) then
                 write(stdout,'(/a)') utils_banner('=', &
                      'Chemical softness')
              else
                 write(stdout,'(/a)') utils_banner('=', &
                      'Chemical softness')
              end if
           end if


           ! Construct K-KSK

           made_kerws=.false.
           if(.not.denskern%workspace_created) then
              made_kerws=.true.
              call kernel_workspace_create(denskern, rep%overlap)
           end if
           call kernel_validate_ksk(denskern, rep%overlap)
           call sparse_embed_array_create(KmKSK, denskern%ksk)
           call sparse_embed_array_copy(KmKSK, denskern%kern)
           call sparse_embed_array_axpy(KmKSK, denskern%ksk, -1.0_dp)

           if(made_kerws) then
              call kernel_workspace_destroy(denskern)
           end if

           ! scale by beta(?)
           call sparse_embed_array_scale(KmKSK,1.0_dp/pub_edft_smearing_width)

           ! construct dn/dmu

           allocate(chemsoft_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
                mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
           call utils_alloc_check('internal_plot_dens_and_pot &
                &(properties_calculate)','chemsoft_fine',ierr)


           ! jme: KPOINTS_DANGER: check this factor for more than 1 k-point
           if (pub_num_spins == 1) then
              call sparse_embed_array_scale(KmKSK, pub_spin_fac)
           end if

           ! jcap: loop over regions
           chemsoft_fine=0.0_DP
           do isub=1,mdl%nsub
              do jsub=1,mdl%nsub
                 reg_density_fine=0.0_DP
                 call sparse_embed_extract_from_array(kern_array,&
                      KmKSK%m(:,PUB_1K),isub,jsub)
                 call density_on_grid(reg_density_fine, mdl%fine_grid, mdl%dbl_grid, &
                      mdl%cell, mdl%fftbox, kern_array, rep%ngwf_overlap%m(isub,jsub), &
                      rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                      rep%ngwfs_on_grid(jsub), ngwf_basis(jsub))
                 chemsoft_fine = chemsoft_fine + reg_density_fine
                 call sparse_embed_destroy_extracted_array(kern_array)
              end do
           end do

           ! output

           if (pub_aug) then ! PAW version

              fine_ld1 = mdl%fine_grid%ld1
              fine_ld2 = mdl%fine_grid%ld2
              fine_max_slabs12 = mdl%fine_grid%max_slabs12
              allocate(nhat_chemsoft_grad(fine_ld1,fine_ld2,fine_max_slabs12,&
                   pub_num_spins,0:pub_aug_den_dim),stat=ierr)
              call utils_alloc_check('internal_plot_dens_and_pot &
                   &(properties_calculate)','nhat_chemsoft_grad',ierr)
              allocate(sub_nhat_chemsoft_grad(fine_ld1,fine_ld2,fine_max_slabs12,&
                   pub_num_spins,0:pub_aug_den_dim),stat=ierr)
              call utils_alloc_check('internal_plot_dens_and_pot &
                   &(properties_calculate)','sub_nhat_chemsoft_grad',ierr)
              allocate(chemsoft_plot(fine_ld1,fine_ld2,fine_max_slabs12,&
                   pub_num_spins),stat=ierr)
              call utils_alloc_check('internal_plot_dens_and_pot &
                   &(properties_calculate)','chemsoft_plot',ierr)

              ! CREATE PSEUDO DENSITY (INCLUDING NHAT) AND OUTPUT

              nhat_chemsoft_grad = 0.0_DP
              do isub=1,mdl%nsub
                 do jsub=1,mdl%nsub
                    sub_nhat_chemsoft_grad=0.0_DP
                    call sparse_embed_extract_from_array(kern_array,&
                         KmKSK%m(:,PUB_1K),isub,jsub)

                    call augmentation_density_on_grid(sub_nhat_chemsoft_grad, mdl%fine_grid, &
                         mdl%cell,mdl%regions(isub)%pseudo_sp,mdl%regions(jsub)%paw_sp,mdl%aug_box, &
                         kern_array, rep%sp_overlap%m(isub,jsub))

                    nhat_chemsoft_grad = nhat_chemsoft_grad + sub_nhat_chemsoft_grad
                    call sparse_embed_destroy_extracted_array(kern_array)
                 end do
              end do

              chemsoft_plot = chemsoft_fine + nhat_chemsoft_grad(:,:,:,:,0)

              ! ndmh: spin polarisation: calculate total density and spin density
              if (pub_num_spins == 2) then
                 chemsoft_plot(:,:,:,1) = chemsoft_plot(:,:,:,UP) + &
                      chemsoft_plot(:,:,:,DN)
                 chemsoft_plot(:,:,:,2) = chemsoft_plot(:,:,:,UP) - &
                      2.0_DP * chemsoft_plot(:,:,:,DN)
              end if

              ! cks: output chemsoft in plot format file
              ! vm: output chemsoft in Angstrom rather than in Bohr
              ! jd: ... but leave the unit conversion to visual_scalarfield
              call visual_scalarfield( &
                   chemsoft_plot(:,:,:,1), mdl%fine_grid, mdl%cell, &
                   'Electronic pseudo chemical softness (in e/ang^3) for:', '_ps_chemsoft', &
                   mdl%elements, ANGSTROM**3)
              if (pub_num_spins == 2) call visual_scalarfield( &
                   chemsoft_plot(:,:,:,2), mdl%fine_grid, mdl%cell, &
                   'Electronic pseudo spin chemical softness (in e/ang^3) for:', &
                   '_ps_spinchemsoft', mdl%elements, ANGSTROM**3)


              ! CREATE AE Chemical Softness AND OUTPUT

              nhat_chemsoft_grad(:,:,:,:,0) = 0.0_DP
              ! rc2013: copy rhoij across to avoid temporary arrays
              call sparse_embed_extract_from_array(kern_array,rhoij)
              call paw_sphere_density_on_grid(nhat_chemsoft_grad(:,:,:,:,0), &
                   mdl%fine_grid,mdl%cell,mdl%aug_box, &
                   kern_array,1.0_DP,-1.0_DP,mdl%regions(1)%paw_sp)
              call sparse_embed_destroy_extracted_array(kern_array)
              chemsoft_plot = chemsoft_fine + nhat_chemsoft_grad(:,:,:,:,0)

              ! ndmh: spin polarisation: calculate total chemsoft and spin chemsoft
              if (pub_num_spins == 2) then
                 chemsoft_plot(:,:,:,1) = chemsoft_plot(:,:,:,UP) + &
                      chemsoft_plot(:,:,:,DN)
                 chemsoft_plot(:,:,:,2) = chemsoft_plot(:,:,:,UP) - &
                      2.0_DP * chemsoft_plot(:,:,:,DN)
              end if

              ! cks: output chemsoft in plot format file
              ! vm: output chemsoft in Angstrom rather than in Bohr
              ! jd: ... but leave the unit conversion to visual_scalarfield
              call visual_scalarfield( &
                   chemsoft_plot(:,:,:,1), mdl%fine_grid, mdl%cell, &
                   'All-electron chemical softness (in e/ang^3) for:', '_ae_chemsoft', &
                   mdl%elements, ANGSTROM**3)
              if (pub_num_spins == 2) call visual_scalarfield( &
                   density_plot(:,:,:,2), mdl%fine_grid, mdl%cell, &
                   'All-electron spin chemical softness (in e/ang^3) for:', &
                   '_ae_spinchemsoft', mdl%elements, ANGSTROM**3)

              deallocate(chemsoft_plot,stat=ierr)
              call utils_dealloc_check('internal_plot_dens_and_pot &
                   &(properties_calculate)','chemsoft_plot',ierr)
              deallocate(nhat_chemsoft_grad,stat=ierr)
              call utils_dealloc_check('internal_plot_dens_and_pot &
                   &(properties_calculate)','nhat_chemsoft_grad',ierr)
              deallocate(sub_nhat_chemsoft_grad,stat=ierr)
              call utils_dealloc_check('internal_plot_dens_and_pot &
                   &(properties_calculate)','sub_nhat_chemsoft_grad',ierr)

           else ! non-PAW version

              ! ndmh: spin polarisation: calculate total chemsoft and spin chemsoft
              if (pub_num_spins == 2) then
                 chemsoft_fine(:,:,:,1) = chemsoft_fine(:,:,:,UP) + &
                      chemsoft_fine(:,:,:,DN)
                 chemsoft_fine(:,:,:,2) = chemsoft_fine(:,:,:,UP) - &
                      2.0_DP * chemsoft_fine(:,:,:,DN)
              end if

              ! cks: output chemsoft in plot format file
              ! vm: output chemsoft in Angstrom rather than in Bohr
              ! jd: ... but leave the unit conversion to visual_scalarfield
              call visual_scalarfield( &
                   chemsoft_fine(:,:,:,1), mdl%fine_grid, mdl%cell, &
                   'Electronic chemical softness (in e/ang^3) for:', '_chemsoft', mdl%elements, &
                   ANGSTROM**3)
              if (pub_num_spins == 2) call visual_scalarfield( &
                   chemsoft_fine(:,:,:,2), mdl%fine_grid, mdl%cell, &
                   'Electronic spin chemical softness (in e/ang^3) for:', '_spinchemsoft', &
                   mdl%elements, ANGSTROM**3)

           end if


           call sparse_embed_array_destroy(KmKSK)

           deallocate(chemsoft_fine,stat=ierr)
           call utils_dealloc_check('internal_plot_dens_and_pot &
                &(properties_calculate)','chemsoft_fine',ierr)

        end if

        ! ja531-> end of chemical softness


        ! jd: ###### EFIELD OUTPUT ##########
        if(pub_efield_calculate) then
           ! rc2013: if this is an embedding run only print out active potential
           isub = pub_active_region
           allocate(efield_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
                mdl%fine_grid%max_slabs12, 3), stat=ierr)
           call utils_alloc_check('internal_plot_dens_and_pot &
                &(properties_calculate)','efield_fine_obc', ierr)
           call visual_efield_on_grid_obc(efield_fine, lhxc_fine(:,:,:,1,isub), &
                mdl%fine_grid, mdl%cell, pub_finite_difference_order)
           call visual_scalarfield(efield_fine(:,:,:,1), &
                mdl%fine_grid, mdl%cell, &
                'Local electric field (Ion+Hartree+XC) [MV/cm] for:', &
                '_efield_obc_x', mdl%elements, -HARTREE_PER_BOHR_TO_MV_PER_CM)
           call visual_scalarfield(efield_fine(:,:,:,2), &
                mdl%fine_grid, mdl%cell, &
                'Local electric field (Ion+Hartree+XC) [MV/cm] for:', &
                '_efield_obc_y', mdl%elements, -HARTREE_PER_BOHR_TO_MV_PER_CM)
           call visual_scalarfield(efield_fine(:,:,:,3), &
                mdl%fine_grid, mdl%cell, &
                'Local electric field (Ion+Hartree+XC) [MV/cm] for:', &
                '_efield_obc_z', mdl%elements, -HARTREE_PER_BOHR_TO_MV_PER_CM)

           call visual_efield_on_grid_pbc(efield_fine, lhxc_fine(:,:,:,1,isub), &
                mdl%fine_grid, mdl%cell)
           call visual_scalarfield(efield_fine(:,:,:,1), &
                mdl%fine_grid, mdl%cell, &
                'Local electric field (Ion+Hartree+XC) [MV/cm] for:', &
                '_efield_pbc_x', mdl%elements, -HARTREE_PER_BOHR_TO_MV_PER_CM)
           call visual_scalarfield(efield_fine(:,:,:,2), &
                mdl%fine_grid, mdl%cell, &
                'Local electric field (Ion+Hartree+XC) [MV/cm] for:', &
                '_efield_pbc_y', mdl%elements, -HARTREE_PER_BOHR_TO_MV_PER_CM)
           call visual_scalarfield(efield_fine(:,:,:,3), &
                mdl%fine_grid, mdl%cell, &
                'Local electric field (Ion+Hartree+XC) [MV/cm] for:', &
                '_efield_pbc_z', mdl%elements, -HARTREE_PER_BOHR_TO_MV_PER_CM)
           deallocate(efield_fine, stat=ierr)
           call utils_dealloc_check('internal_plot_dens_and_pot &
                &(properties_calculate)','efield_fine', ierr)
        end if
        ! jd: ###### END EFIELD OUTPUT ##########

        if(pub_aug) then
           if(allocated(density_plot)) then
              deallocate(density_plot,stat=ierr)
              call utils_dealloc_check('internal_plot_dens_and_pot &
                   &(properties_calculate)','density_plot',ierr)
           end if
           if(allocated(nhat_den_grad)) then
              deallocate(nhat_den_grad,stat=ierr)
              call utils_dealloc_check('internal_plot_dens_and_pot &
                   &(properties_calculate)','nhat_den_grad',ierr)
           end if
           if(allocated(sub_nhat_den_grad)) then
              deallocate(sub_nhat_den_grad,stat=ierr)
              call utils_dealloc_check('internal_plot_dens_and_pot &
                   &(properties_calculate)','sub_nhat_den_grad',ierr)
           end if
        end if

        ! cks: Deallocate charge density slabs
        deallocate(density_fine,stat=ierr)
        call utils_dealloc_check('internal_plot_dens_and_pot &
             &(properties_calculate)','density_fine',ierr)
        deallocate(reg_density_fine,stat=ierr)
        call utils_dealloc_check('internal_plot_dens_and_pot &
             &(properties_calculate)','buffer_fine',ierr)

        if (pub_on_root) write(stdout,'(a)') repeat('=',80)


      end subroutine internal_plot_dens_and_pot


      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine internal_plot_nbo
      ! lpl: Reads GENNBO output files to plot relevant orbitals
      !      Uses properties_plot_mo
      ! lpl: 01052012 - Updated to read spin-polarized orbitals
      !                 nbongwf_tr initialized after getting nspin
      !                 from appropriate NBO output file
      ! jcap: modified for embedding 28/09/18

        use augmentation, only: augmentation_overlap
        use comms, only: comms_barrier, comms_bcast, &
             pub_root_proc_id
        use constants, only: max_spins
        use dense, only: DEM, dense_create, dense_destroy, dense_get_col
        use eigenstates, only: eigenstates_plot_mo
        use npa, only: npa_read_nbodata
        use rundat, only: pub_aug, pub_nbo_list_plotnbo, pub_nbo_plotorbtype, &
             pub_rootname, pub_num_kpoints, PUB_1K
        use utils, only: utils_abort, utils_assert

        implicit none

        ! Local Variables
        type(SPAM3_EMBED) :: mo_kern(max_spins), aug_overlap
        type(DEM), allocatable :: nbongwf_tr(:)
        ! agrecocmplx: use FUNCTIONS type instead
        !real(kind=DP), allocatable :: col_buffer(:)
        type(FUNCTIONS) :: col_buffer
        integer :: it, num_ao, num_ngwfs
        character(len=64) :: output_tag, output_file, title, nbo_trfile, &
             spin_flag
        logical :: pre_ao

        integer :: is, nspin ! # spins of orbitals in NBO's output files

        if(pub_debug_on_root) write(stdout,'(a)') 'Entering internal_plot_nbo'

        ! agrecocmplx
        call utils_assert(.not.any(rep%ngwfs_on_grid(:)%iscmplx), 'Error in &
             & properties_calculate: subroutine not yet checked for complex NGWFs.')

        ! jme: KPOINTS_DANGER
        ! Some parts of this subroutine enforce a single k-point.
        call utils_assert(pub_num_kpoints == PUB_1K, &
             'Subroutine internal_plot_nbo not ready yet for more&
             & than one k-point.')

        ! jcap: abort if more than one subsystem
        call utils_assert(mdl%nsub.eq.1, 'Error in properties_calculate: &
             &NBO not yet ready for more than one subsystem')

        ! jcap: calculate total number of ngwfs
        num_ngwfs=ngwf_basis(1)%num

        ! lpl: Initialize matrices and arrays
        ! agrecocmplx: allocate using appropriate routine
        call data_functions_alloc(col_buffer,num_ngwfs,iscmplx=.false.)

        ! ndmh: NBO 'kernel' and aug overlap part for dealing with augmentation
        if (pub_aug) then
           call sparse_embed_create(mo_kern(1), denskern%kern%m(1,PUB_1K))
           call sparse_embed_create(aug_overlap,rep%overlap)
           call augmentation_overlap(aug_overlap%p,mdl%pseudo_sp,mdl%paw_sp, &
                rep%sp_overlap%p)
        end if

        ! lpl: Determine type of orbital and appropriate NBO input filename
        pub_nbo_plotorbtype = adjustl(pub_nbo_plotorbtype)
        if(pub_on_root) then
           select case(trim(pub_nbo_plotorbtype))
              case ('PNAO')
                 pre_ao = .true.
                 nbo_trfile = 'NONE'
                 output_tag = '_pnao'
              case ('NAO')
                 pre_ao = .false.
                 nbo_trfile = 'NONE'
                 output_tag = '_nao'
              case ('PNHO')
                 pre_ao = .true.
                 nbo_trfile = '.35'
                 output_tag = '_pnho'
              case ('NHO')
                 pre_ao = .false.
                 nbo_trfile = '.35'
                 output_tag = '_nho'
              case ('PNBO')
                 pre_ao = .true.
                 nbo_trfile = '.37'
                 output_tag = '_pnbo'
              case ('NBO')
                 pre_ao = .false.
                 nbo_trfile = '.37'
                 output_tag = '_nbo'
              case ('NLMO')
                 pre_ao = .false.
                 nbo_trfile = '.39'
                 output_tag = '_nlmo'
              case default
                 call utils_abort('Invalid nbo_plot_orbtype')
           end select

           nbo_trfile = adjustl(nbo_trfile)
           if(nbo_trfile /= 'NONE') then
              write(nbo_trfile,'(a)') &
                   trim(pub_rootname)//'_nao'//trim(nbo_trfile)
              nbo_trfile = adjustl(nbo_trfile)
           end if

           output_tag = adjustl(output_tag)

        end if

        call comms_bcast(pub_root_proc_id,pre_ao)
        call comms_bcast(pub_root_proc_id,nbo_trfile)
        call comms_bcast(pub_root_proc_id,output_tag)

        ! lpl: Inquire if input file from nbo_trfile is spin-polarized
        if(nbo_trfile /= 'NONE') then
           call npa_read_nbodata(nbongwf_tr,ngwf_basis(1),trim(nbo_trfile), &
                pre_ao,num_ao,spin_inquire=.true.,nspin=nspin,&
                par=mdl%par)
        else
           nspin = 1
        end if
        ! lpl: Allocate nbongwf_tr based on result of nspin (only 1 & 2 valid)
        allocate(nbongwf_tr(nspin),stat=ierr)
        call utils_alloc_check('properties_plot_nbo','nbongwf_tr',ierr)
        do is=1,nspin
           call dense_create(nbongwf_tr(is),num_ngwfs,num_ngwfs)
        end do

        ! lpl: Obtain full NGWF to NBO's orbital transformation
        ! lpl: 01052012 - Re-read if spin-polarized flag
        call npa_read_nbodata(nbongwf_tr,ngwf_basis(1),trim(nbo_trfile), &
             pre_ao,num_ao,spin_inquire=.false.,nspin=nspin,&
             par=mdl%par)

        if(pub_debug_on_root) then
           write(stdout,'(a)') 'DEBUG info: NBOs read for plotting'
           write(stdout,'(a)') '----------------------------------'
           do it=1,pub_nbo_list_plotnbo(0)
              write(stdout,'(1x,I5,1x,I5)') it,pub_nbo_list_plotnbo(it)
           end do
           write(stdout,'(a)') '----------------------------------'
        end if

        ! lpl: Plot orbitals
        do is=1,nspin

           ! lpl: Spin label
           if(nspin > 1) then
              if(is == 1) then
                 spin_flag = '_alpha_'
              else if(is == 2) then
                 spin_flag = '_beta_'
              end if
           else
              spin_flag = ''
           end if

           do it=1,pub_nbo_list_plotnbo(0)
              write(title,'(I5)') pub_nbo_list_plotnbo(it)

              output_file &
                   = trim(output_tag)//trim(spin_flag)//trim(adjustl(title))
              title = trim(pub_nbo_plotorbtype)//' # '//trim(adjustl(title))

              if(pub_nbo_list_plotnbo(it) <= num_ao) then
                 ! lpl: Get NBO in run-time NGWF basis
                 ! agrecocmplx: currently only real case implemented
                 call dense_get_col(col_buffer%d,nbongwf_tr(is), &
                      pub_nbo_list_plotnbo(it))
                 call comms_barrier

                 ! agrecocmplx
                 call eigenstates_plot_mo(col_buffer,rep%ngwfs_on_grid(1), &
                      ngwf_basis(1),mdl,trim(adjustl(title)), &
                      trim(adjustl(output_file)),aug_overlap%m(1,1),mo_kern(1)%m(1,1))

              else
                 write(stdout,'(a,I5)') &
                      'ERROR: # orbs > num_ao. Skipping...',it
              end if

              call comms_barrier
           end do

        end do ! END do is=1,nspin

        do is=nspin,1,-1
           call dense_destroy(nbongwf_tr(is))
        end do
        deallocate(nbongwf_tr,stat=ierr)
        call utils_dealloc_check('properties_plot_nbo','nbongwf_tr',ierr)

        if (pub_aug) then
           call sparse_embed_destroy(aug_overlap)
           call sparse_embed_destroy(mo_kern(1))
        end if

        ! agrecocmplx: deallocate using appropriate routine
        call data_functions_dealloc(col_buffer)

        if(pub_debug_on_root) write(stdout,'(a)') 'Exiting internal_plot_nbo'

      end subroutine internal_plot_nbo


  end subroutine properties_calculate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ngwfs_orthogonalise(orth_on_grid,ngwfs_on_grid,&
       ngwf_basis,overlap,ireg)

    !==========================================================================!
    ! Symmetrically orthogonalise NGWFs on each atom using the Lowdin          !
    ! transformation.                                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! ngwfs_on_grid (input)    : NGWFs in ppd representation on this proc      !
    ! ngwf_basis (input)       : The function basis for the NGWFs              !
    ! orth_on_grid (output)    : Orthogonalised NGWFs in ppd representation    !
    !--------------------------------------------------------------------------!
    ! Written by Mark Robinson, January 2008                                   !
    ! Cleaned up and made to use new code by Nicholas Hine, 08/06/2010         !
    ! Moved to ngwfs_mod by Nicholas Hine.                                     !
    !==========================================================================!

    use datatypes, only: FUNCTIONS, data_set_to_zero
    use comms, only: pub_my_proc_id
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_sum_ppd_funcs
    use linalg, only: linalg_dsygv_lt, linalg_zhegv_lt
    use parallel_strategy, only: par=>pub_par
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_get_block, &
         sparse_put_block
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNCTIONS), intent(inout)    :: orth_on_grid(1:1)
    type(FUNCTIONS), intent(in)     :: ngwfs_on_grid
    type(SPAM3), intent(in) :: overlap
    integer, intent(in), optional :: ireg

    ! BLAS subroutines
    external :: zgemm, dgemm

    ! Local variables
    type(SPAM3) :: dmat(1)
    integer :: num_ngwfs
    integer :: loc_iat, iat
    integer :: i
    integer :: ierr
    ! agrecocmplx
    logical :: loc_cmplx

    ! matrices and arrays
    real(kind=DP), allocatable, dimension(:,:) :: overlap_block_real, &
        transform_real,d_real,temp_real
    ! agrecocmplx: needed for complex case
    complex(kind=DP), allocatable, dimension(:,:) :: overlap_block_cmplx, &
        transform_cmplx,d_cmplx,temp_cmplx
    real(kind=DP), allocatable, dimension(:)   :: evals
    character(len=1) :: ireg_str

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering ngwfs_orthogonalise'

    ! Start timer
    call timer_clock('ngwfs_orthogonalise',1)

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid%iscmplx

    ! Allocate temporary arrays
    num_ngwfs = ngwf_basis%max_on_atom

    ! agrecocmplx
    ! complex case
    if (loc_cmplx) then
       allocate(overlap_block_cmplx(num_ngwfs,num_ngwfs),stat=ierr)
       call utils_alloc_check('ngwfs_orthogonalise','overlap_block_cmplx',ierr)
       allocate(transform_cmplx(num_ngwfs,num_ngwfs),stat=ierr)
       call utils_alloc_check('ngwfs_orthogonalise','transform_cmplx',ierr)
       allocate(d_cmplx(num_ngwfs,num_ngwfs),stat=ierr)
       call utils_alloc_check('ngwfs_orthogonalise','d_cmplx',ierr)
       allocate(temp_cmplx(num_ngwfs,num_ngwfs),stat=ierr)
       call utils_alloc_check('ngwfs_orthogonalise','temp_cmplx',ierr)
    ! real case
    else
       allocate(overlap_block_real(num_ngwfs,num_ngwfs),stat=ierr)
       call utils_alloc_check('ngwfs_orthogonalise','overlap_block_real',ierr)
       allocate(transform_real(num_ngwfs,num_ngwfs),stat=ierr)
       call utils_alloc_check('ngwfs_orthogonalise','transform_real',ierr)
       allocate(d_real(num_ngwfs,num_ngwfs),stat=ierr)
       call utils_alloc_check('ngwfs_orthogonalise','d_real',ierr)
       allocate(temp_real(num_ngwfs,num_ngwfs),stat=ierr)
       call utils_alloc_check('ngwfs_orthogonalise','temp_real',ierr)
    end if


    allocate(evals(num_ngwfs),stat=ierr)
    call utils_alloc_check('ngwfs_orthogonalise','evals',ierr)

    ! agrecocmplx
    if(present(ireg)) then
       write(ireg_str,'(i1)') ireg
       dmat(1)%structure = 'D'//ireg_str//ireg_str
    else
       dmat(1)%structure = 'D'
    end if
    call sparse_create(dmat(1),iscmplx=loc_cmplx)

    ! Loop over atoms on this proc
    ! agrecocmplx
    ! complex case
    if (loc_cmplx) then
       do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)

          iat = loc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1

          ! Find number of NGWFs on atom
          num_ngwfs = ngwf_basis%num_on_atom(iat)

          ! Only need to orthogonalise if more than one NGWF
          if (num_ngwfs > 1) then

             ! Get lower triangle of overlap matrix of NGWFs on atom
             overlap_block_cmplx = (0.0_DP,0.0_DP)
             call sparse_get_block(overlap_block_cmplx,overlap,iat,iat)

             ! Solve the eigenvalue problem
             transform_cmplx = overlap_block_cmplx
             temp_cmplx = (0.0_DP,0.0_DP)
             do i=1,num_ngwfs
                temp_cmplx(i,i) = (1.0_DP,0.0_DP)
             enddo

             ! agrecocmplx: call complex version of routine
             call linalg_zhegv_lt(transform_cmplx,evals,temp_cmplx,num_ngwfs)

             ! Construct d(-1/2) from evals
             d_cmplx = (0.0_DP,0.0_DP)
             do i=1,num_ngwfs
                d_cmplx(i,i) = cmplx(1 / sqrt(evals(i)),0.0,kind=DP)
             enddo

             ! Transform d^(-1/2) to S^(-1/2)
             ! X=d*transform'
             ! agrecocmplx: call complex version and take conjugate
             call zgemm('N','C',num_ngwfs,num_ngwfs,num_ngwfs,(1.0_DP,0.0_DP), &
                        d_cmplx,num_ngwfs,transform_cmplx,num_ngwfs,&
                        (0.0_DP,0.0_DP),temp_cmplx,num_ngwfs)
             ! S^-1/2=transform*X
             call zgemm('N','N',num_ngwfs,num_ngwfs,num_ngwfs,(1.0_DP,0.0_DP), &
                        transform_cmplx,num_ngwfs,temp_cmplx,num_ngwfs,&
                        (0.0_DP,0.0_DP),d_cmplx,num_ngwfs)

             ! Put resulting matrix block into dmat
             call sparse_put_block(d_cmplx,dmat(1),iat,iat)

          else
             ! Copy block straight from overlap matrix
             call sparse_get_block(d_cmplx,overlap,iat,iat)
             call sparse_put_block(d_cmplx,dmat(1),iat,iat)
          endif

       enddo
    ! real case
    else
       do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)

          iat = loc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1

          ! Find number of NGWFs on atom
          num_ngwfs = ngwf_basis%num_on_atom(iat)

          ! Only need to orthogonalise if more than one NGWF
          if (num_ngwfs > 1) then

             ! Get lower triangle of overlap matrix of NGWFs on atom
             overlap_block_real = 0.0_DP
             call sparse_get_block(overlap_block_real,overlap,iat,iat)

             ! Solve the eigenvalue problem
             transform_real = overlap_block_real
             temp_real = 0.0_DP
             do i=1,num_ngwfs
                temp_real(i,i) = 1.0_DP
             enddo
             call linalg_dsygv_lt(transform_real,evals,temp_real,num_ngwfs)

             ! Construct d(-1/2) from evals
             d_real = 0.0_DP
             do i=1,num_ngwfs
                d_real(i,i) = 1 / sqrt(evals(i))
             enddo

             ! Transform d^(-1/2) to S^(-1/2)
             ! X=d*transform'
             call dgemm('N','T',num_ngwfs,num_ngwfs,num_ngwfs,1.0_DP, &
                        d_real,num_ngwfs,transform_real,num_ngwfs,&
                        0.0_DP,temp_real,num_ngwfs)
             ! S^-1/2=transform*X
             call dgemm('N','N',num_ngwfs,num_ngwfs,num_ngwfs,1.0_DP, &
                        transform_real,num_ngwfs,temp_real,num_ngwfs,&
                        0.0_DP,d_real,num_ngwfs)

             ! Put resulting matrix block into dmat
             call sparse_put_block(d_real,dmat(1),iat,iat)

          else
             ! Copy block straight from overlap matrix
             call sparse_get_block(d_real,overlap,iat,iat)
             call sparse_put_block(d_real,dmat(1),iat,iat)
          endif

       enddo
    end if

    call data_set_to_zero(orth_on_grid(1))
    call function_ops_sum_ppd_funcs(orth_on_grid,ngwf_basis, &
         dmat,1,1,dmat(1),ngwfs_on_grid,ngwf_basis)

    ! Deallocate temporary arrays
    call sparse_destroy(dmat(1))
    deallocate(evals,stat=ierr)
    call utils_dealloc_check('ngwfs_orthogonalise','evals',ierr)

    ! agrecocmplx
    if (loc_cmplx) then
       deallocate(temp_cmplx,stat=ierr)
       call utils_dealloc_check('ngwfs_orthogonalise','temp_cmplx',ierr)
       deallocate(d_cmplx,stat=ierr)
       call utils_dealloc_check('ngwfs_orthogonalise','d_cmplx',ierr)
       deallocate(transform_cmplx,stat=ierr)
       call utils_dealloc_check('ngwfs_orthogonalise','transform_cmplx',ierr)
       deallocate(overlap_block_cmplx,stat=ierr)
       call utils_dealloc_check('ngwfs_orthogonalise','overlap_block_cmplx',ierr)
    else
       deallocate(temp_real,stat=ierr)
       call utils_dealloc_check('ngwfs_orthogonalise','temp_real',ierr)
       deallocate(d_real,stat=ierr)
       call utils_dealloc_check('ngwfs_orthogonalise','d_real',ierr)
       deallocate(transform_real,stat=ierr)
       call utils_dealloc_check('ngwfs_orthogonalise','transform_real',ierr)
       deallocate(overlap_block_real,stat=ierr)
       call utils_dealloc_check('ngwfs_orthogonalise','overlap_block_real',ierr)
    end if

    ! Stop timer
    call timer_clock('ngwfs_orthogonalise',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving ngwfs_orthogonalise'

  end subroutine ngwfs_orthogonalise


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine properties_ngwfs_char(ngwfs_on_grid,ngwf_basis,mdl,max_l, &
       nmax,output,ireg)

    !==========================================================================!
    ! Prints out an analysis of the s/p/d/f character of the NGWFs             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! ngwfs_on_grid (input)    : NGWFs in ppd representation on this proc      !
    ! ngwf_basis (input)       : The function basis for the NGWFs              !
    ! lmax (input)             : Largest angular momentum SW to include        !
    ! nmax (input)             : Size of spherical wave basis set              !
    ! output (input)           : Should we write out SW coefficients to file   !
    ! ireg (input)             : NGWF region that we are analysing             !
    !--------------------------------------------------------------------------!
    ! Written by Mark Robinson, November 2007                                  !
    ! Modified to remove pub_par by Robert Charlton, 23/08/2018.               !
    !==========================================================================!

    use basis, only: basis_extract_function_from_box, &
         basis_location_func_wrt_cell
    use comms, only: pub_on_root, pub_root_proc_id, pub_total_num_procs, &
         pub_my_proc_id, comms_recv, comms_send, comms_barrier, comms_reduce
    use constants, only: DP, stdout, two_pi
    use datatypes, only: COEF, data_coef_abs, data_functions_dot, &
         FUNCTIONS, data_functions_alloc, data_functions_dealloc, &
         FFTBOX_DATA, data_fftbox_alloc, data_fftbox_dealloc
    use function_basis, only: FUNC_BASIS
    use geometry, only: point, operator(.DOT.), operator(*), operator(+)
    use model_type, only: MODEL
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_open_unit_check, utils_close_unit_check, utils_assert
    use spherical_wave, only: sw_bessel_zeros,sw_init, sw_recp_generate
    use timer, only: timer_clock

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(MODEL), intent(in) :: mdl
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
    integer, intent(in)       :: max_l
    integer, optional, intent(in) :: nmax
    logical, optional, intent(in) :: output
    integer, optional, intent(in) :: ireg

    ! default to not printing out SW coefficients
    logical :: sw2file                           ! if true write SW coeffs to file
    integer :: max_n

    ! counters
    integer       :: atom,orig_atom,loc_ngwf,ngwf
    integer       :: l,m,n,proc        ! loop counters
    integer       :: ngwf_on_atom      ! more counters

    ! variables
    type(POINT)   :: fbcentre,centre,frac,disp
    real(kind=DP) :: radius,norm
    type(COEF)    :: tot
    integer       :: ierr,offset,npoints !,lastp
    integer       :: sx,sy,sz
    integer       :: loc_ireg

    ! spherical bessel function zeros array
    real(kind=DP), allocatable :: qnl(:)

    ! SW FFTbox and tightbox
    type(FFTBOX_DATA) :: sw_fftbox
    type(FUNCTIONS) :: sw_on_grid
    complex(kind=DP), allocatable :: sw_work(:,:,:)

    ! summation arrays
    real(kind=DP), allocatable :: my_spd_sum(:,:)
    real(kind=DP), allocatable :: snd_spd_sum(:,:)
    real(kind=DP), allocatable :: spd_sum(:,:)

    ! file data (when printing out SW coefficients)
    integer           :: fileunit
    character(len=16) :: filename

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering properties_ngwfs_char'

    ! start timer
    call timer_clock('properties_ngwfs_char',1)

    !jmecmplx
    call utils_assert(.not. ngwfs_on_grid%iscmplx, 'Error:&
         & properties_ngwfs_char not ready yet for complex NGWFs.')

    ! rc2013: set the regions indicator if allocated
    loc_ireg = 1
    if(present(ireg)) loc_ireg = ireg

    ! check optional arguments
    if (present (nmax)) then
       max_n = nmax
    else
       ! default is number of psincs across largest NGWF
       max_n = ceiling(maxval(ngwf_basis%spheres(:)%radius) / &
                       min(mdl%cell%d1,mdl%cell%d2,mdl%cell%d3))
       call comms_reduce('MAX',max_n)
    endif

    sw2file = .false.
    if (present (output)) sw2file = output

    ! allocate qnl array
    allocate(qnl(max_n),stat=ierr)
    call utils_alloc_check('properties_ngwfs_char','qnl',ierr)

    ! init SW module if not already - with max kr
    radius = maxval(ngwf_basis%spheres(:)%radius)**2 * 2 &
         * maxval(mdl%fftbox%recip_grid(5,:,:,:))
    call comms_reduce('MAX',radius)
    call sw_init(max_l,max_n,radius)

    ! calculate shift required to centre SW in FFTbox on a grid point
    fbcentre = real(mdl%fftbox%total_pt1-1,kind=DP)/2/mdl%fftbox%total_pt1 * mdl%fftbox%a1 + &
               real(mdl%fftbox%total_pt2-1,kind=DP)/2/mdl%fftbox%total_pt2 * mdl%fftbox%a2 + &
               real(mdl%fftbox%total_pt3-1,kind=DP)/2/mdl%fftbox%total_pt3 * mdl%fftbox%a3

    ! obtain free unit and number of ngwfs prior to this proc
    ! (if writing SW coefficients to file)
    fileunit = utils_unit ()

    ! allocate workspace
    tot%iscmplx = ngwfs_on_grid%iscmplx
    allocate(my_spd_sum(0:max_l+2,ngwf_basis%proc_num),stat=ierr)
    call utils_alloc_check('properties_ngwfs_char','my_spd_sum',ierr)
    call data_fftbox_alloc(sw_fftbox, mdl%fftbox%total_ld1, mdl%fftbox%total_ld2,&
         mdl%fftbox%total_pt3, iscmplx=ngwfs_on_grid%iscmplx) !jmecmplx
    allocate(sw_work(mdl%fftbox%total_ld1,mdl%fftbox%total_ld2,&
         mdl%fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('properties_ngwfs_char','sw_work',ierr)
    ! allocate sw_on_grid array
    call data_functions_alloc(sw_on_grid, &
         ngwf_basis%max_n_ppds_sphere*mdl%cell%n_pts, ngwfs_on_grid%iscmplx)

    ! zero summation array
    my_spd_sum = 0.0_DP

    ! loop over ngwfs on this proc
    loc_ngwf = 0
    do ngwf=ngwf_basis%first_on_proc(pub_my_proc_id), &
         ngwf_basis%first_on_proc(pub_my_proc_id+1)-1
       loc_ngwf = loc_ngwf + 1

       orig_atom = mdl%regions(loc_ireg)%par%orig_atom(ngwf_basis%atom_of_func(ngwf))
       ngwf_on_atom = ngwf - &
            ngwf_basis%first_on_atom(ngwf_basis%atom_of_func(ngwf)) + 1

       ! open file to write SW coefficients to if requested
       if (sw2file) then
          write(filename,'(a5,i5.5,a6)') 'atom_',orig_atom, &
               'ngwf_',ngwf_on_atom, '.swdat'
          open(unit=fileunit,file=trim(filename),action='write',iostat=ierr)
          call utils_open_unit_check('properties_ngwfs_char','fileunit',ierr)
       endif

       ! NGWF center and location in large array
       centre = ngwf_basis%spheres(loc_ngwf)%centre
       radius = ngwf_basis%spheres(loc_ngwf)%radius
       offset = ngwf_basis%spheres(loc_ngwf)%offset
       npoints = ngwf_basis%spheres(loc_ngwf)%n_ppds_sphere * mdl%cell%n_pts
!       lastp = offset + npoints -1

       ! calculate product of NGWF with itself
       my_spd_sum(max_l+1,loc_ngwf) = data_coef_abs(data_functions_dot( &
            ngwfs_on_grid, ngwfs_on_grid, offset, offset, npoints))  !jmecmplx

       ! calculate fractional coordinates of ngwf centre
       frac%X = (centre .DOT. mdl%cell%b1) / two_pi * real(mdl%cell%total_pt1,kind=DP)
       frac%Y = (centre .DOT. mdl%cell%b2) / two_pi * real(mdl%cell%total_pt2,kind=DP)
       frac%Z = (centre .DOT. mdl%cell%b3) / two_pi * real(mdl%cell%total_pt3,kind=DP)

       ! calculate displacement of centre from nearest grid point
       disp = mod(frac%X,1.0_DP) / mdl%fftbox%total_pt1 * mdl%fftbox%a1 + &
              mod(frac%Y,1.0_DP) / mdl%fftbox%total_pt2 * mdl%fftbox%a2 + &
              mod(frac%Z,1.0_DP) / mdl%fftbox%total_pt3 * mdl%fftbox%a3

       ! calculate first grid point of tightbox wrt cell
       call basis_location_func_wrt_cell(sx,sy,sz, &
            ngwf_basis%tight_boxes(loc_ngwf),mdl%cell)

       ! calculate location of tightbox wrt FFTbox and displacement of centre wrt
       ! FFTbox depending whether FFTbox side coincides with sim cell
       if(mdl%fftbox%coin1) then
          sx = 1
          disp%X = centre%X
       else
          sx = (mdl%fftbox%total_pt1-1)/2 - floor(frac%X) + sx
          disp%X = fbcentre%X + disp%X
       endif
       if(mdl%fftbox%coin2) then
          sy = 1
          disp%Y = centre%Y
       else
          sy = (mdl%fftbox%total_pt2-1)/2 - floor(frac%Y) + sy
          disp%Y = fbcentre%Y + disp%Y
       endif
       if(mdl%fftbox%coin3) then
          sz = 1
          disp%Z = centre%Z
       else
          sz = (mdl%fftbox%total_pt3-1)/2 - floor(frac%Z) + sz
          disp%Z = fbcentre%Z + disp%Z
       endif

       ! loop over angular momentum l
       do l=0,max_l

          ! generate bessel function roots for current NGWF radius and l
          qnl(:) = sw_bessel_zeros(:,l) / radius

          ! loop over azimuthal angular momentum m
          do m=-l,+l

             ! write l,m to file if requested
             if (sw2file) write(fileunit,'(/2i4,1x)',advance='no') l, m

             ! loop over n
             do n=1,max_n

                ! calculate SW in an FFTbox
                call sw_recp_generate(sw_fftbox%d,sw_work,l,m,qnl(n),radius, & !jmecmplx
                     disp,mdl%fftbox)

                ! extract SW from FFTbox into ppds
                call basis_extract_function_from_box(sw_on_grid, sw_fftbox, &
                     ngwf_basis%spheres(loc_ngwf), &
                     ngwf_basis%tight_boxes(loc_ngwf),sx,sy,sz,1, &
                     mdl%cell,mdl%fftbox)

                ! perform integrals !jmecmplx
                tot = data_functions_dot(sw_on_grid, ngwfs_on_grid, 1, offset, npoints) !jmecmplx
                norm = data_coef_abs(data_functions_dot(sw_on_grid, sw_on_grid, 1, 1, npoints))

                ! sum the square of the integral after normalising spherical wave
                my_spd_sum(l,loc_ngwf) = my_spd_sum(l,loc_ngwf) + data_coef_abs(tot)**2/norm

                ! write SW coefficient to file if requested
                if (sw2file) then
                   if (tot%iscmplx) then
                      write(fileunit,'(2(f10.5))',advance='no') tot%z*sqrt(mdl%cell%weight/norm)
                   else
                      write(fileunit,'(f10.5)',advance='no') tot%d*sqrt(mdl%cell%weight/norm)
                   end if
                end if

             enddo  ! end n loop
          enddo     ! end m loop
       enddo        ! end l loop

       ! calculate total and spilling
       my_spd_sum(max_l+2,loc_ngwf) = sum(my_spd_sum(:max_l,loc_ngwf))
       my_spd_sum(max_l+1,loc_ngwf) = my_spd_sum(max_l+1,loc_ngwf) - &
            my_spd_sum(max_l+2,loc_ngwf)

       if (sw2file) then
          close(unit=fileunit,iostat=ierr)
          call utils_close_unit_check('properties_ngwfs_char','fileunit',ierr)
       end if

    enddo  ! end ngwf loop

    ! weight integrals
    my_spd_sum(:,:) = my_spd_sum(:,:) * mdl%cell%weight

    ! deallocate sw_on_grid array
    call data_functions_dealloc(sw_on_grid)
    ! deallocate workspace
    deallocate(sw_work,stat=ierr)
    call utils_dealloc_check('properties_ngwfs_char','sw_work',ierr)
    call data_fftbox_dealloc(sw_fftbox)

    ! deallocate qnl array
    deallocate(qnl,stat=ierr)
    call utils_dealloc_check('properties_ngwfs_char','qnl',ierr)

    ! allocate workspace
    allocate(snd_spd_sum(0:max_l+2,ngwf_basis%max_on_proc),stat=ierr)
    call utils_alloc_check('properties_ngwfs_char','snd_spd_sum',ierr)
    allocate(spd_sum(0:max_l+2,ngwf_basis%num),stat=ierr)
    call utils_alloc_check('properties_ngwfs_char','spd_sum',ierr)

    ! gather my_spd_sums to root proc
    ! this is effectively an MPI_GATHERV
    if (pub_on_root) then
       offset = 1
       do proc=0,pub_total_num_procs-1
          if (proc .eq. pub_root_proc_id) then
             spd_sum(:,offset:offset+ngwf_basis%proc_num-1) = &
                my_spd_sum(:,:ngwf_basis%proc_num)
          else
             call comms_recv(proc,snd_spd_sum)
             spd_sum(:,offset:offset+ngwf_basis%num_on_proc(proc)-1) = &
                snd_spd_sum(:,:ngwf_basis%num_on_proc(proc))
          endif
          offset = offset + ngwf_basis%num_on_proc(proc)
       enddo
    else
       call comms_send(pub_root_proc_id,my_spd_sum)
    endif
    call comms_barrier

    ! write out table
    if (pub_on_root) then
       write(stdout,'(/a,i2.2,a)') '************************ NGWF s/p/d/f Character Analysis &
                             &*********** (n=',max_n,') ****'
       write(stdout,'(a1,3x,2x,1x,a6,1x,a6,4x,a1,2x,4a9,4x,a9,3x,a1)') &
          '|','Atom','NGWF',':','%s','%p','%d','%f','spilling','|'
       do atom=1,mdl%regions(loc_ireg)%par%nat
          ngwf = ngwf_basis%first_on_atom(mdl%regions(loc_ireg)%par%distr_atom(atom)) - 1
          do ngwf_on_atom=1,ngwf_basis%num_on_atom( &
               mdl%regions(loc_ireg)%par%distr_atom(atom))
             ngwf = ngwf + 1           ! original NGWF counting
             if (ngwf_on_atom .eq. 1) then
                write(stdout,'(a1,3x,a2,1x,i6,1x,i6,4x,a1,2x,4f9.3,4x,e9.2,3x,a1)') &
                   '|',mdl%regions(loc_ireg)%elements(atom)%symbol,atom,ngwf_on_atom,':',&
                   100.0_DP*spd_sum(:3,ngwf)/spd_sum(max_l+2,ngwf),&
                   spd_sum(max_l+1,ngwf),'|'
             else
                write(stdout,'(a1,3x,2x,1x,6x,1x,i6,4x,a1,2x,4f9.3,4x,e9.2,3x,a1)') &
                   '|',ngwf_on_atom,':',&
                   100.0_DP*spd_sum(:3,ngwf)/spd_sum(max_l+2,ngwf),&
                   spd_sum(max_l+1,ngwf),'|'
             endif
          enddo
       enddo
       write(stdout,'(a/)') '****************************************&
                             &****************************************'
    endif

    ! deallocate workspace
    deallocate(snd_spd_sum,stat=ierr)
    call utils_dealloc_check('properties_ngwfs_char','snd_spd_sum',ierr)
    deallocate(spd_sum,stat=ierr)
    call utils_dealloc_check('properties_ngwfs_char','spd_sum',ierr)
    deallocate(my_spd_sum,stat=ierr)
    call utils_dealloc_check('properties_ngwfs_char','my_spd_sum',ierr)

    ! stop timer
    call timer_clock('properties_ngwfs_char',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving properties_ngwfs_char'

  end subroutine properties_ngwfs_char


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine properties_polarisation_simcell(ngwfs_on_grid,ngwf_basis,&
       mdl,denskern,overlap)

    !==========================================================================!
    ! Calculates dipole moment of a molecule by integrating over the           !
    ! simulation cell (cannot handle atoms moving across cell boundary)        !
    !                                                                          !
    ! For large systems this will be slower than properties_polarisation and   !
    ! will also need a larger stack size                                       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! ngwfs_on_grid (input)    : NGWFs in ppd representation on this proc      !
    ! ngwf_basis (input)       : The function basis of NGWFs                   !
    ! denskern (input)         : Density kernel                                !
    ! overlap (input)          : Overlap matrix                                !
    !--------------------------------------------------------------------------!
    ! Written by Mark Robinson, July 2008                                      !
    ! Modified 22/01/2011 by Nicholas Hine to use cell_grid_real_pt routine to !
    ! save on storage.                                                         !
    ! Modified 12/08/2014 by Jacek Dziedzic to calculate the quadrupole, based !
    ! on CKS's notes.                                                          !
    ! Modified for embedding by Joseph Prentice, August 2018                   !
    !==========================================================================!

    use datatypes, only: FUNCTIONS
    use cell_grid, only: cell_grid_real_pt
    use comms, only: comms_reduce, pub_on_root
    use constants, only: DP, stdout, CRLF
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use density, only: density_on_grid
    use multipole_ops, only: multipole_2_cart_traceless_to_sph_real
    use rundat, only: pub_num_spins, pub_spin_fac, pub_polarisation_simcell_refpt, &
         pub_num_kpoints, PUB_1K
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use sparse, only: SPAM3
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_banner, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(MODEL), intent(in) :: mdl
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid(:)
    type(SPAM3_EMBED_ARRAY), intent(in)   :: denskern
    type(SPAM3_EMBED), intent(in)   :: overlap

    ! Local Variables
    integer       :: islab12,i2,i1,is,isub,jsub,ierr
    real(kind=DP), allocatable :: density_fine(:,:,:,:)
    real(kind=DP), allocatable :: buffer_fine(:,:,:,:)
    type(SPAM3), allocatable :: kern_array(:)
    real(kind=DP) :: d(3),q(3),t(3),rpt(3)

    real(kind=DP) :: refpt(3) ! = pub_polarisation_simcell_refpt
    real(kind=DP) :: rpt0(3) ! = rpt-refpt

    ! jd: Cartesian primitive quadrupole tensor, electronic
    real(kind=DP) :: qpole_elec_p(3,3)
    ! jd: Cartesian primitive quadrupole tensor, ionic
    real(kind=DP) :: qpole_ionic_p(3,3)
    ! jd: Cartesian primitive quadrupole tensor, total
    real(kind=DP) :: qpole_total_p(3,3)
    ! jd: Cartesian traceless quadrupole tensor, electronic, CKS's convention
    real(kind=DP) :: qpole_elec_t(3,3)
    ! jd: Cartesian traceless quadrupole tensor, ionic, CKS's convention
    real(kind=DP) :: qpole_ionic_t(3,3)
    ! jd: Cartesian traceless quadrupole tensor, total, CKS's convention
    real(kind=DP) :: qpole_total_t(3,3)
    ! jd: Cartesian traceless quadrupole tensor, electronic, Gaussian convention
    real(kind=DP) :: qpole_elec_tg(3,3)
    ! jd: Cartesian traceless quadrupole tensor, ionic, Gaussian convention
    real(kind=DP) :: qpole_ionic_tg(3,3)
    ! jd: Cartesian traceless quadrupole tensor, total, Gaussian convention
    real(kind=DP) :: qpole_total_tg(3,3)
    ! jd: Components of the spherical quadrupole tensor, electronic
    !     c: real, s: imaginary
    real(kind=DP) :: qpole_elec_20, qpole_elec_21c, qpole_elec_21s, &
         qpole_elec_22c, qpole_elec_22s
    ! jd: Components of the spherical quadrupole tensor, ionic
    !     c: real, s: imaginary
    real(kind=DP) :: qpole_ionic_20, qpole_ionic_21c, qpole_ionic_21s, &
         qpole_ionic_22c, qpole_ionic_22s
    ! jd: Components of the spherical quadrupole tensor, total
    !     c: real, s: imaginary
    real(kind=DP) :: qpole_total_20, qpole_total_21c, qpole_total_21s, &
         qpole_total_22c, qpole_total_22s

    integer       :: qpole_i, qpole_j ! jd: placeholders for x, y, z
    integer       :: out_unit
    character(len=47) :: unit_name

    real(kind=DP), parameter :: THIRD = 1.0_DP/3.0_DP
    real(kind=DP), parameter :: au_to_DAsq = 1.345034536_DP

    ! -------------------------------------------------------------------------

    ! Start timer
    call timer_clock('properties_polarisation_simcell',1)

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine properties_polarisation_simcell not ready yet for more&
         & than one k-point.')

    refpt = pub_polarisation_simcell_refpt

    ! calculate dipole moment for each coordinate
    d = 0.0_DP ! ddor: Set dipole moment initially to zero
    qpole_elec_p = 0.0_DP   ! jd: Same for quadrupole moment (primitive)
    qpole_elec_t = 0.0_DP   !     (and traceless)
    qpole_elec_tg = 0.0_DP

    ! calculate density on fine grid

    ! rc2013: allocate arrays
    allocate(density_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,pub_num_spins), stat=ierr)
    call utils_alloc_check('properties_polarisation_simcell',&
         'density_fine',ierr)
    allocate(buffer_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,pub_num_spins), stat=ierr)
    call utils_alloc_check('properties_polarisation_simcell',&
         'buffer_fine',ierr)

    ! jcap: allocate temporary kernel
    allocate(kern_array(pub_num_spins),stat=ierr)
    call utils_alloc_check('properties_polarisation_simcell',&
         'kern_array',ierr)

    ! jcap: loop over regions
    density_fine=0.0_DP
    do isub=1,mdl%nsub
       do jsub=1,mdl%nsub
          ! jcap: create temporary kernel
          call sparse_embed_extract_from_array(kern_array,denskern%m(:,PUB_1K),&
               isub,jsub)

          call density_on_grid(buffer_fine,mdl%fine_grid,mdl%dbl_grid, &
               mdl%cell, mdl%fftbox, kern_array,overlap%m(isub,jsub),&
               ngwfs_on_grid(isub),ngwf_basis(jsub), &
               ngwfs_on_grid(jsub),ngwf_basis(jsub))

          density_fine=density_fine+buffer_fine

          ! jcap: destroy temporary kernel
          call sparse_embed_destroy_extracted_array(kern_array)
       end do
    end do

    deallocate(kern_array,stat=ierr)
    call utils_dealloc_check('properties_polarisation_simcell',&
         'kern_array',ierr)

    ! loop over fine grid points and over spins
    do islab12=1,mdl%fine_grid%num_my_slabs12
       do i2=1,mdl%fine_grid%n2
          do i1=1,mdl%fine_grid%n1
             call cell_grid_real_pt(rpt,i1,i2,islab12,mdl%fine_grid)
             rpt0 = rpt - refpt
             do is=1,pub_num_spins
                ! integrate density*position over simulation cell
                 d(1:3) = d(1:3) + density_fine(i1,i2,islab12,is) * rpt0(1:3)
                 ! jd: Calculate primitive and traceless quadrupole moment
                 !     tensors, cf. (8) in cks's notes "DMA expansionn of the
                 !     charge density in ONETEP" for the latter.
                 do qpole_i = 1, 3
                    do qpole_j = 1, 3
                       qpole_elec_p(qpole_i,qpole_j) = &
                            qpole_elec_p(qpole_i,qpole_j) + &
                            density_fine(i1,i2,islab12,is) * &
                            rpt0(qpole_i)*rpt0(qpole_j)
                       qpole_elec_t(qpole_i,qpole_j) = &
                            qpole_elec_t(qpole_i,qpole_j) + &
                            density_fine(i1,i2,islab12,is) * ( &
                            1.5_DP*rpt0(qpole_i)*rpt0(qpole_j) - &
                            merge(0.5_DP, 0.0_DP, qpole_i == qpole_j) * &
                            sum(rpt0(:)*rpt0(:)))
                       qpole_elec_tg(qpole_i,qpole_j) = &
                            qpole_elec_tg(qpole_i,qpole_j) + &
                            density_fine(i1,i2,islab12,is) * ( &
                            1.0_DP*rpt0(qpole_i)*rpt0(qpole_j) - &
                            merge(THIRD, 0.0_DP, qpole_i == qpole_j) * &
                            sum(rpt0(:)*rpt0(:)))
                    end do
                 end do
             end do
          end do
       enddo
    enddo

    ! rc2013: free up memory
    deallocate(density_fine, stat=ierr)
    call utils_dealloc_check('properties_polarisation_simcell',&
         'density_fine',ierr)
    deallocate(buffer_fine, stat=ierr)
    call utils_dealloc_check('properties_polarisation_simcell',&
         'buffer_fine',ierr)

    ! sum contributions from all procs
    call comms_reduce('SUM',d)
    call comms_reduce('SUM',qpole_elec_p)
    call comms_reduce('SUM',qpole_elec_t)
    call comms_reduce('SUM',qpole_elec_tg)

    ! weight integrals
    d = d * mdl%fine_grid%weight
    qpole_elec_p = qpole_elec_p * mdl%fine_grid%weight
    qpole_elec_t = qpole_elec_t * mdl%fine_grid%weight
    qpole_elec_tg = qpole_elec_tg * mdl%fine_grid%weight

    ! multiplication factor depending on number of spins
    if (pub_num_spins /= 2) then
       d(:) = pub_spin_fac * d(:)
       qpole_elec_p(:,:) = pub_spin_fac * qpole_elec_p(:,:)
       qpole_elec_t(:,:) = pub_spin_fac * qpole_elec_t(:,:)
       qpole_elec_tg(:,:) = pub_spin_fac * qpole_elec_tg(:,:)
    end if

    ! jd: Quadrupole convention conversion, all these are electronic
    call multipole_2_cart_traceless_to_sph_real(qpole_elec_22s, &      ! output
         qpole_elec_21s,qpole_elec_20,qpole_elec_21c,qpole_elec_22c, & ! output
         qpole_elec_t(1,1), qpole_elec_t(1,2), qpole_elec_t(1,3), &    ! input
         qpole_elec_t(2,2), qpole_elec_t(2,3), qpole_elec_t(3,3), &    ! input
         .true.)

    ! calculate dipole moment from ions
    q(1) = sum(mdl%elements(:)%ion_charge * (mdl%elements(:)%centre%X-refpt(1)))
    q(2) = sum(mdl%elements(:)%ion_charge * (mdl%elements(:)%centre%Y-refpt(2)))
    q(3) = sum(mdl%elements(:)%ion_charge * (mdl%elements(:)%centre%Z-refpt(3)))

    ! calculate quadrupole primitive Cartesian moment from ions
    qpole_ionic_p(1,1) = sum(mdl%elements(:)%ion_charge * &
         (mdl%elements(:)%centre%X-refpt(1)) * &
         (mdl%elements(:)%centre%X-refpt(1)))
    qpole_ionic_p(1,2) = sum(mdl%elements(:)%ion_charge * &
         (mdl%elements(:)%centre%X-refpt(1)) * &
         (mdl%elements(:)%centre%Y-refpt(2)))
    qpole_ionic_p(1,3) = sum(mdl%elements(:)%ion_charge * &
         (mdl%elements(:)%centre%X-refpt(1)) * &
         (mdl%elements(:)%centre%Z-refpt(3)))
    qpole_ionic_p(2,1) = sum(mdl%elements(:)%ion_charge * &
         (mdl%elements(:)%centre%Y-refpt(2)) * &
         (mdl%elements(:)%centre%X-refpt(1)))
    qpole_ionic_p(2,2) = sum(mdl%elements(:)%ion_charge * &
         (mdl%elements(:)%centre%Y-refpt(2)) * &
         (mdl%elements(:)%centre%Y-refpt(2)))
    qpole_ionic_p(2,3) = sum(mdl%elements(:)%ion_charge * &
         (mdl%elements(:)%centre%Y-refpt(2)) * &
         (mdl%elements(:)%centre%Z-refpt(3)))
    qpole_ionic_p(3,1) = sum(mdl%elements(:)%ion_charge * &
         (mdl%elements(:)%centre%Z-refpt(3)) * &
         (mdl%elements(:)%centre%X-refpt(1)))
    qpole_ionic_p(3,2) = sum(mdl%elements(:)%ion_charge * &
         (mdl%elements(:)%centre%Z-refpt(3)) * &
         (mdl%elements(:)%centre%Y-refpt(2)))
    qpole_ionic_p(3,3) = sum(mdl%elements(:)%ion_charge * &
         (mdl%elements(:)%centre%Z-refpt(3)) * &
         (mdl%elements(:)%centre%Z-refpt(3)))

    ! calculate quadrupole traceless Cartesian moment from ions, CKS's conv'ion
    qpole_ionic_t(1,1) = sum(mdl%elements(:)%ion_charge * &
         (1.5_DP*(mdl%elements(:)%centre%X-refpt(1)) * &
         (mdl%elements(:)%centre%X-refpt(1)) - &
         0.5_DP * ((mdl%elements(:)%centre%X-refpt(1)) * &
         (mdl%elements(:)%centre%X-refpt(1)) + &
         (mdl%elements(:)%centre%Y-refpt(2)) * &
         (mdl%elements(:)%centre%Y-refpt(2)) + &
         (mdl%elements(:)%centre%Z-refpt(3)) * &
         (mdl%elements(:)%centre%Z-refpt(3)))))
    qpole_ionic_t(1,2) = sum(mdl%elements(:)%ion_charge * &
         1.5_DP*(mdl%elements(:)%centre%X-refpt(1)) * &
         (mdl%elements(:)%centre%Y-refpt(2)))
    qpole_ionic_t(1,3) = sum(mdl%elements(:)%ion_charge * &
         1.5_DP*(mdl%elements(:)%centre%X-refpt(1)) * &
         (mdl%elements(:)%centre%Z-refpt(3)))
    qpole_ionic_t(2,1) = sum(mdl%elements(:)%ion_charge * &
         1.5_DP*(mdl%elements(:)%centre%Y-refpt(2)) * &
         (mdl%elements(:)%centre%X-refpt(1)))
    qpole_ionic_t(2,2) = sum(mdl%elements(:)%ion_charge * &
         (1.5_DP*(mdl%elements(:)%centre%Y-refpt(2)) * &
         (mdl%elements(:)%centre%Y-refpt(2)) - &
         0.5_DP * ((mdl%elements(:)%centre%X-refpt(1)) * &
         (mdl%elements(:)%centre%X-refpt(1)) + &
         (mdl%elements(:)%centre%Y-refpt(2)) * &
         (mdl%elements(:)%centre%Y-refpt(2)) + &
         (mdl%elements(:)%centre%Z-refpt(3)) * &
         (mdl%elements(:)%centre%Z-refpt(3)))))
    qpole_ionic_t(2,3) = sum(mdl%elements(:)%ion_charge * &
         1.5_DP*(mdl%elements(:)%centre%Y-refpt(2)) * &
         (mdl%elements(:)%centre%Z-refpt(3)))
    qpole_ionic_t(3,1) = sum(mdl%elements(:)%ion_charge * &
         1.5_DP*(mdl%elements(:)%centre%Z-refpt(3)) * &
         (mdl%elements(:)%centre%X-refpt(1)))
    qpole_ionic_t(3,2) = sum(mdl%elements(:)%ion_charge * &
         1.5_DP*(mdl%elements(:)%centre%Z-refpt(3)) * &
         (mdl%elements(:)%centre%Y-refpt(2)))
    qpole_ionic_t(3,3) = sum(mdl%elements(:)%ion_charge * &
         (1.5_DP*(mdl%elements(:)%centre%Z-refpt(3)) * &
         (mdl%elements(:)%centre%Z-refpt(3)) - &
         0.5_DP * ((mdl%elements(:)%centre%X-refpt(1)) * &
         (mdl%elements(:)%centre%X-refpt(1)) + &
         (mdl%elements(:)%centre%Y-refpt(2)) * &
         (mdl%elements(:)%centre%Y-refpt(2)) + &
         (mdl%elements(:)%centre%Z-refpt(3)) * &
         (mdl%elements(:)%centre%Z-refpt(3)))))

    ! calculate quadrupole traceless Cartesian moment from ions, Gaussian conv'n
    qpole_ionic_tg(1,1) = sum(mdl%elements(:)%ion_charge * &
         (1.0_DP*(mdl%elements(:)%centre%X-refpt(1)) * &
         (mdl%elements(:)%centre%X-refpt(1)) - &
         THIRD * ((mdl%elements(:)%centre%X-refpt(1)) * &
         (mdl%elements(:)%centre%X-refpt(1)) + &
         (mdl%elements(:)%centre%Y-refpt(2)) * &
         (mdl%elements(:)%centre%Y-refpt(2)) + &
         (mdl%elements(:)%centre%Z-refpt(3)) * &
         (mdl%elements(:)%centre%Z-refpt(3)))))
    qpole_ionic_tg(1,2) = sum(mdl%elements(:)%ion_charge * &
         1.0_DP*(mdl%elements(:)%centre%X-refpt(1)) * &
         (mdl%elements(:)%centre%Y-refpt(2)))
    qpole_ionic_tg(1,3) = sum(mdl%elements(:)%ion_charge * &
         1.0_DP*(mdl%elements(:)%centre%X-refpt(1)) * &
         (mdl%elements(:)%centre%Z-refpt(3)))
    qpole_ionic_tg(2,1) = sum(mdl%elements(:)%ion_charge * &
         1.0_DP*(mdl%elements(:)%centre%Y-refpt(2)) * &
         (mdl%elements(:)%centre%X-refpt(1)))
    qpole_ionic_tg(2,2) = sum(mdl%elements(:)%ion_charge * &
         (1.0_DP*(mdl%elements(:)%centre%Y-refpt(2)) * &
         (mdl%elements(:)%centre%Y-refpt(2)) - &
         THIRD * ((mdl%elements(:)%centre%X-refpt(1)) * &
         (mdl%elements(:)%centre%X-refpt(1)) + &
         (mdl%elements(:)%centre%Y-refpt(2)) * &
         (mdl%elements(:)%centre%Y-refpt(2)) + &
         (mdl%elements(:)%centre%Z-refpt(3)) * &
         (mdl%elements(:)%centre%Z-refpt(3)))))
    qpole_ionic_tg(2,3) = sum(mdl%elements(:)%ion_charge * &
         1.0_DP*(mdl%elements(:)%centre%Y-refpt(2)) * &
         (mdl%elements(:)%centre%Z-refpt(3)))
    qpole_ionic_tg(3,1) = sum(mdl%elements(:)%ion_charge * &
         1.0_DP*(mdl%elements(:)%centre%Z-refpt(3)) * &
         (mdl%elements(:)%centre%X-refpt(1)))
    qpole_ionic_tg(3,2) = sum(mdl%elements(:)%ion_charge * &
         1.0_DP*(mdl%elements(:)%centre%Z-refpt(3)) * &
         (mdl%elements(:)%centre%Y-refpt(2)))
    qpole_ionic_tg(3,3) = sum(mdl%elements(:)%ion_charge * &
         (1.0_DP*(mdl%elements(:)%centre%Z-refpt(3)) * &
         (mdl%elements(:)%centre%Z-refpt(3)) - &
         THIRD * ((mdl%elements(:)%centre%X-refpt(1)) * &
         (mdl%elements(:)%centre%X-refpt(1)) + &
         (mdl%elements(:)%centre%Y-refpt(2)) * &
         (mdl%elements(:)%centre%Y-refpt(2)) + &
         (mdl%elements(:)%centre%Z-refpt(3)) * &
         (mdl%elements(:)%centre%Z-refpt(3)))))

    ! jd: Quadrupole convention conversion, all these are ionic
    call multipole_2_cart_traceless_to_sph_real(qpole_ionic_22s, &      ! output
         qpole_ionic_21s,qpole_ionic_20,qpole_ionic_21c,qpole_ionic_22c, & ! out
         qpole_ionic_t(1,1), qpole_ionic_t(1,2), qpole_ionic_t(1,3), &  ! input
         qpole_ionic_t(2,2), qpole_ionic_t(2,3), qpole_ionic_t(3,3), &  ! input
         .true.)

    ! jd: Flip electronic quantities to the 'cores are positive' convention
    d = -d
    qpole_elec_p = -qpole_elec_p
    qpole_elec_t = -qpole_elec_t
    qpole_elec_tg = -qpole_elec_tg

    ! calculate total dipole moment
    t = q + d

    ! calculate total quadrupole moment
    qpole_total_p = qpole_ionic_p + qpole_elec_p
    qpole_total_t = qpole_ionic_t + qpole_elec_t
    qpole_total_tg = qpole_ionic_tg + qpole_elec_tg

    ! jd: Quadrupole convention conversion, all these are total
    call multipole_2_cart_traceless_to_sph_real(qpole_total_22s, &      ! output
         qpole_total_21s,qpole_total_20,qpole_total_21c,qpole_total_22c, & ! out
         qpole_total_t(1,1), qpole_total_t(1,2), qpole_total_t(1,3), &  ! input
         qpole_total_t(2,2), qpole_total_t(2,3), qpole_total_t(3,3), &  ! input
         .true.)


    ! print pretty table
    if(pub_on_root) then
       write(stdout,'(/a)') utils_banner('=',&
            'Dipole Moment Calculation (SimCell)')
       write(stdout,'(a34,4x,a)') 'Type:','Single Molecule'
       write(stdout,*)
       write(stdout,'(a34,4x,a4,1x,f24.12)') 'Electronic dipole moment &
            &(e.bohr):','dx =',d(1)
       write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',d(2)
       write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',d(3)
       write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
            sqrt(sum(d(:)**2))
       write(stdout,*)
       write(stdout,'(a34,4x,a4,1x,f24.12)') 'Ionic dipole moment (e.bohr):', &
            'dx =',q(1)
       write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',q(2)
       write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',q(3)
       write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
            sqrt(sum(q(:)**2))
       write(stdout,*)
       write(stdout,'(a34,4x,a4,1x,f24.12)') 'Total dipole moment (e.bohr):', &
            'dx =',t(1)
       write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',t(2)
       write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',t(3)
       write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
            sqrt(sum(t(:)**2))
       write(stdout,*)
       write(stdout,'(a)') repeat('=',80)
       write(stdout,'(/a)') CRLF//'======================== Quadrupole Moment &
                            &Calculation (SimCell) ==============='
       write(stdout,'(a34,4x,a)') 'Type:','Single Molecule'
       write(stdout,*)

       do out_unit = 1, 2

          if(out_unit == 1) then
             unit_name = '(e.bohr^2)'
          else
             unit_name = '(Debye.A)  '
          end if
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '*** ELECTRONIC ***',           'Qxx =',qpole_elec_p(1,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Cartesian primitive quadrupole moments', &
                                               'Qxy =',qpole_elec_p(1,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               adjustr(unit_name),             'Qxz =',qpole_elec_p(1,3)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qyx =',qpole_elec_p(2,1)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qyy =',qpole_elec_p(2,2)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qyz =',qpole_elec_p(2,3)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qzx =',qpole_elec_p(3,1)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qzy =',qpole_elec_p(3,2)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qzz =',qpole_elec_p(3,3)
          write(stdout,'(41x,4x,a11,1x,f24.9)') 'magnitude =', &
               sqrt(sum(qpole_elec_p(:,:)*qpole_elec_p(:,:)))
          write(stdout,'(45x,4x,a7,1x,f24.9)') 'trace =', &
               qpole_elec_p(1,1)+qpole_elec_p(2,2)+qpole_elec_p(3,3)
          write(stdout,*)
          ! -------------------------------------------------------------------
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '*** ELECTRONIC ***',           'Qxx =',qpole_elec_t(1,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Cartesian traceless quadrupole moments', &
                                               'Qxy =',qpole_elec_t(1,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               adjustr(unit_name),             'Qxz =',qpole_elec_t(1,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'following the convention of',  'Qyx =',qpole_elec_t(2,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Gray-Gubbins, whereby', &
                                               'Qyy =',qpole_elec_t(2,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'traceless Qij = 3/2 * Qij - (Qxx+Qyy+Qzz)/2', &
                                               'Qyz =',qpole_elec_t(2,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i=j,                   ', &
                                               'Qzx =',qpole_elec_t(3,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '         and  = 3/2 * Qij                  ', &
                                               'Qzy =',qpole_elec_t(3,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i/=j                   ', &
                                               'Qzz =',qpole_elec_t(3,3)
          write(stdout,'(41x,4x,a11,1x,f24.9)') 'magnitude =', &
               sqrt(sum(qpole_elec_t(:,:)*qpole_elec_t(:,:)))
          write(stdout,'(28x,4x,a24,1x,f24.9)') 'trace (should be zero) =', &
               qpole_elec_t(1,1)+qpole_elec_t(2,2)+qpole_elec_t(3,3)
          write(stdout,*)
          ! -------------------------------------------------------------------
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '*** ELECTRONIC ***',           'Qxx =',2.0_DP*qpole_elec_t(1,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Cartesian traceless quadrupole moments', &
                                               'Qxy =',2.0_DP*qpole_elec_t(1,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               adjustr(unit_name),             'Qxz =',2.0_DP*qpole_elec_t(1,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'following the convention of',  'Qyx =',2.0_DP*qpole_elec_t(2,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Jackson, whereby', &
                                               'Qyy =',2.0_DP*qpole_elec_t(2,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'traceless Qij = 3 * Qij - (Qxx+Qyy+Qzz)', &
                                               'Qyz =',2.0_DP*qpole_elec_t(2,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i=j,               ', &
                                               'Qzx =',2.0_DP*qpole_elec_t(3,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '         and  = 3 * Qij                ', &
                                               'Qzy =',2.0_DP*qpole_elec_t(3,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i/=j               ', &
                                               'Qzz =',2.0_DP*qpole_elec_t(3,3)
          write(stdout,'(41x,4x,a11,1x,f24.9)') 'magnitude =', &
               2.0_DP*sqrt(sum(qpole_elec_t(:,:)*qpole_elec_t(:,:)))
          write(stdout,'(28x,4x,a24,1x,f24.9)') 'trace (should be zero) =', &
               2.0_DP*qpole_elec_t(1,1)+ &
               2.0_DP*qpole_elec_t(2,2)+ &
               2.0_DP*qpole_elec_t(3,3)
          write(stdout,*)
          ! -------------------------------------------------------------------
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '*** ELECTRONIC ***',           'Qxx =',qpole_elec_tg(1,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Cartesian traceless quadrupole moments', &
                                               'Qxy =',qpole_elec_tg(1,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               adjustr(unit_name),             'Qxz =',qpole_elec_tg(1,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'following the convention of',  'Qyx =',qpole_elec_tg(2,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Gaussian, whereby', &
                                               'Qyy =',qpole_elec_tg(2,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'traceless Qij = Qij - (Qxx+Qyy+Qzz)/3', &
                                               'Qyz =',qpole_elec_tg(2,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i=j,             ', &
                                               'Qzx =',qpole_elec_tg(3,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '          and = Qij                  ', &
                                               'Qzy =',qpole_elec_tg(3,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i/=j             ', &
                                               'Qzz =',qpole_elec_tg(3,3)
          write(stdout,'(41x,4x,a11,1x,f24.9)') 'magnitude =', &
               sqrt(sum(qpole_elec_tg(:,:)*qpole_elec_tg(:,:)))
          write(stdout,'(28x,4x,a24,1x,f24.9)') 'trace (should be zero) =', &
               qpole_elec_tg(1,1)+qpole_elec_tg(2,2)+qpole_elec_tg(3,3)
          write(stdout,*)
          ! -------------------------------------------------------------------
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '*** ELECTRONIC ***',           'Q20 =', qpole_elec_20
          write(stdout,'(a46,4x,a6,1x,f24.9)') &
               'Spherical quadrupole moments', 'Q21c =', qpole_elec_21c
          write(stdout,'(46x,4x,a6,1x,f24.9)') 'Q21s =', qpole_elec_21s
          write(stdout,'(46x,4x,a6,1x,f24.9)') 'Q22c =', qpole_elec_22c
          write(stdout,'(46x,4x,a6,1x,f24.9)') 'Q22s =', qpole_elec_22s
          write(stdout,*)
          write(stdout,'(a)') repeat('-',80)
          write(stdout,*)
          ! -------------------------------------------------------------------
          ! -------------------------------------------------------------------
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '*** IONIC ***',                'Qxx =',qpole_ionic_p(1,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Cartesian primitive quadrupole moments', &
                                               'Qxy =',qpole_ionic_p(1,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               adjustr(unit_name),             'Qxz =',qpole_ionic_p(1,3)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qyx =',qpole_ionic_p(2,1)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qyy =',qpole_ionic_p(2,2)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qyz =',qpole_ionic_p(2,3)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qzx =',qpole_ionic_p(3,1)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qzy =',qpole_ionic_p(3,2)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qzz =',qpole_ionic_p(3,3)
          write(stdout,'(41x,4x,a11,1x,f24.9)') 'magnitude =', &
               sqrt(sum(qpole_ionic_p(:,:)*qpole_ionic_p(:,:)))
          write(stdout,'(45x,4x,a7,1x,f24.9)') 'trace =', &
               qpole_ionic_p(1,1)+qpole_ionic_p(2,2)+qpole_ionic_p(3,3)
          write(stdout,*)
          ! ---------------------------------------------------------------------
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '*** IONIC ***',                'Qxx =',qpole_ionic_t(1,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Cartesian traceless quadrupole moments', &
                                               'Qxy =',qpole_ionic_t(1,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               adjustr(unit_name),             'Qxz =',qpole_ionic_t(1,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'following the convention of',  'Qyx =',qpole_ionic_t(2,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Gray-Gubbins, whereby', &
                                               'Qyy =',qpole_ionic_t(2,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'traceless Qij = 3/2 * Qij - (Qxx+Qyy+Qzz)/2', &
                                               'Qyz =',qpole_ionic_t(2,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i=j,                   ', &
                                               'Qzx =',qpole_ionic_t(3,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '         and  = 3/2 * Qij                  ', &
                                               'Qzy =',qpole_ionic_t(3,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i/=j                   ', &
                                               'Qzz =',qpole_ionic_t(3,3)
          write(stdout,'(41x,4x,a11,1x,f24.9)') 'magnitude =', &
               sqrt(sum(qpole_ionic_t(:,:)*qpole_ionic_t(:,:)))
          write(stdout,'(28x,4x,a24,1x,f24.9)') 'trace (should be zero) =', &
               qpole_ionic_t(1,1)+qpole_ionic_t(2,2)+qpole_ionic_t(3,3)
          write(stdout,*)
          ! -------------------------------------------------------------------
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '*** IONIC ***',                'Qxx =',2.0_DP*qpole_ionic_t(1,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Cartesian traceless quadrupole moments', &
                                               'Qxy =',2.0_DP*qpole_ionic_t(1,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               adjustr(unit_name),             'Qxz =',2.0_DP*qpole_ionic_t(1,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'following the convention of',  'Qyx =',2.0_DP*qpole_ionic_t(2,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Jackson, whereby', &
                                               'Qyy =',2.0_DP*qpole_ionic_t(2,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'traceless Qij = 3 * Qij - (Qxx+Qyy+Qzz)', &
                                               'Qyz =',2.0_DP*qpole_ionic_t(2,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i=j,               ', &
                                               'Qzx =',2.0_DP*qpole_ionic_t(3,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '         and  = 3 * Qij                ', &
                                               'Qzy =',2.0_DP*qpole_ionic_t(3,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i/=j               ', &
                                               'Qzz =',2.0_DP*qpole_ionic_t(3,3)
          write(stdout,'(41x,4x,a11,1x,f24.9)') 'magnitude =', &
               2.0_DP*sqrt(sum(qpole_ionic_t(:,:)*qpole_ionic_t(:,:)))
          write(stdout,'(28x,4x,a24,1x,f24.9)') 'trace (should be zero) =', &
               2.0_DP*qpole_ionic_t(1,1)+ &
               2.0_DP*qpole_ionic_t(2,2)+ &
               2.0_DP*qpole_ionic_t(3,3)
          write(stdout,*)
          ! -------------------------------------------------------------------
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '*** IONIC ***',                'Qxx =',qpole_ionic_tg(1,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Cartesian traceless quadrupole moments', &
                                               'Qxy =',qpole_ionic_tg(1,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               adjustr(unit_name),             'Qxz =',qpole_ionic_tg(1,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'following the convention of',  'Qyx =',qpole_ionic_tg(2,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Gaussian, whereby', &
                                               'Qyy =',qpole_ionic_tg(2,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'traceless Qij = Qij - (Qxx+Qyy+Qzz)/3', &
                                               'Qyz =',qpole_ionic_tg(2,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i=j,             ', &
                                               'Qzx =',qpole_ionic_tg(3,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '          and = Qij                  ', &
                                               'Qzy =',qpole_ionic_tg(3,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i/=j             ', &
                                               'Qzz =',qpole_ionic_tg(3,3)
          write(stdout,'(41x,4x,a11,1x,f24.9)') 'magnitude =', &
               sqrt(sum(qpole_ionic_tg(:,:)*qpole_ionic_tg(:,:)))
          write(stdout,'(28x,4x,a24,1x,f24.9)') 'trace (should be zero) =', &
               qpole_ionic_tg(1,1)+qpole_ionic_tg(2,2)+qpole_ionic_tg(3,3)
          write(stdout,*)
          ! -------------------------------------------------------------------
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '*** IONIC ***',                'Q20 =', qpole_ionic_20
          write(stdout,'(a46,4x,a6,1x,f24.9)') &
               'Spherical quadrupole moments', 'Q21c =', qpole_ionic_21c
          write(stdout,'(46x,4x,a6,1x,f24.9)') 'Q21s =', qpole_ionic_21s
          write(stdout,'(46x,4x,a6,1x,f24.9)') 'Q22c =', qpole_ionic_22c
          write(stdout,'(46x,4x,a6,1x,f24.9)') 'Q22s =', qpole_ionic_22s
          write(stdout,*)
          write(stdout,'(a)') repeat('-',80)
          write(stdout,*)
          ! -------------------------------------------------------------------
          ! -------------------------------------------------------------------
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '*** TOTAL ***',                'Qxx =',qpole_total_p(1,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Cartesian primitive quadrupole moments', &
                                               'Qxy =',qpole_total_p(1,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               adjustr(unit_name),             'Qxz =',qpole_total_p(1,3)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qyx =',qpole_total_p(2,1)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qyy =',qpole_total_p(2,2)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qyz =',qpole_total_p(2,3)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qzx =',qpole_total_p(3,1)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qzy =',qpole_total_p(3,2)
          write(stdout,'(47x,4x,a5,1x,f24.9)') 'Qzz =',qpole_total_p(3,3)
          write(stdout,'(41x,4x,a11,1x,f24.9)') 'magnitude =', &
               sqrt(sum(qpole_total_p(:,:)*qpole_total_p(:,:)))
          write(stdout,'(45x,4x,a7,1x,f24.9)') 'trace =', &
               qpole_total_p(1,1)+qpole_total_p(2,2)+qpole_total_p(3,3)
          write(stdout,*)
          ! -------------------------------------------------------------------
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '*** TOTAL ***',                'Qxx =',qpole_total_t(1,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Cartesian traceless quadrupole moments', &
                                               'Qxy =',qpole_total_t(1,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               adjustr(unit_name),             'Qxz =',qpole_total_t(1,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'following the convention of',  'Qyx =',qpole_total_t(2,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Gray-Gubbins, whereby', &
                                               'Qyy =',qpole_total_t(2,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'traceless Qij = 3/2 * Qij - (Qxx+Qyy+Qzz)/2', &
                                               'Qyz =',qpole_total_t(2,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i=j,                   ', &
                                               'Qzx =',qpole_total_t(3,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '         and  = 3/2 * Qij                  ', &
                                               'Qzy =',qpole_total_t(3,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i/=j                   ', &
                                               'Qzz =',qpole_total_t(3,3)
          write(stdout,'(41x,4x,a11,1x,f24.9)') 'magnitude =', &
               sqrt(sum(qpole_total_t(:,:)*qpole_total_t(:,:)))
          write(stdout,'(28x,4x,a24,1x,f24.9)') 'trace (should be zero) =', &
               qpole_total_t(1,1)+qpole_total_t(2,2)+qpole_total_t(3,3)
          write(stdout,*)
          ! -------------------------------------------------------------------
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '*** TOTAL ***',                'Qxx =',2.0_DP*qpole_total_t(1,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Cartesian traceless quadrupole moments', &
                                               'Qxy =',2.0_DP*qpole_total_t(1,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               adjustr(unit_name),             'Qxz =',2.0_DP*qpole_total_t(1,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'following the convention of',  'Qyx =',2.0_DP*qpole_total_t(2,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Jackson, whereby', &
                                               'Qyy =',2.0_DP*qpole_total_t(2,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'traceless Qij = 3 * Qij - (Qxx+Qyy+Qzz)', &
                                               'Qyz =',2.0_DP*qpole_total_t(2,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i=j,               ', &
                                               'Qzx =',2.0_DP*qpole_total_t(3,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '         and  = 3 * Qij                ', &
                                               'Qzy =',2.0_DP*qpole_total_t(3,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i/=j               ', &
                                               'Qzz =',2.0_DP*qpole_total_t(3,3)
          write(stdout,'(41x,4x,a11,1x,f24.9)') 'magnitude =', &
               2.0_DP*sqrt(sum(qpole_total_t(:,:)*qpole_total_t(:,:)))
          write(stdout,'(28x,4x,a24,1x,f24.9)') 'trace (should be zero) =', &
               2.0_DP*qpole_total_t(1,1)+ &
               2.0_DP*qpole_total_t(2,2)+ &
               2.0_DP*qpole_total_t(3,3)
          write(stdout,*)
          ! -------------------------------------------------------------------
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '*** TOTAL ***',                'Qxx =',qpole_total_tg(1,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Cartesian traceless quadrupole moments', &
                                               'Qxy =',qpole_total_tg(1,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               adjustr(unit_name),             'Qxz =',qpole_total_tg(1,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'following the convention of',  'Qyx =',qpole_total_tg(2,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'Gaussian, whereby', &
                                               'Qyy =',qpole_total_tg(2,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               'traceless Qij = Qij - (Qxx+Qyy+Qzz)/3', &
                                               'Qyz =',qpole_total_tg(2,3)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i=j,             ', &
                                               'Qzx =',qpole_total_tg(3,1)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '          and = Qij                  ', &
                                               'Qzy =',qpole_total_tg(3,2)
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '                for i/=j             ', &
                                               'Qzz =',qpole_total_tg(3,3)
          write(stdout,'(41x,4x,a11,1x,f24.9)') 'magnitude =', &
               sqrt(sum(qpole_total_tg(:,:)*qpole_total_tg(:,:)))
          write(stdout,'(28x,4x,a24,1x,f24.9)') 'trace (should be zero) =', &
               qpole_total_tg(1,1)+qpole_total_tg(2,2)+qpole_total_tg(3,3)
          write(stdout,*)
          ! -------------------------------------------------------------------
          write(stdout,'(a47,4x,a5,1x,f24.9)') &
               '*** TOTAL ***',                'Q20 =', qpole_total_20
          write(stdout,'(a46,4x,a6,1x,f24.9)') &
               'Spherical quadrupole moments', 'Q21c =', qpole_total_21c
          write(stdout,'(46x,4x,a6,1x,f24.9)') 'Q21s =', qpole_total_21s
          write(stdout,'(46x,4x,a6,1x,f24.9)') 'Q22c =', qpole_total_22c
          write(stdout,'(46x,4x,a6,1x,f24.9)') 'Q22s =', qpole_total_22s
          write(stdout,*)
          write(stdout,'(a)') '************************************************&
                              &********************************'
          write(stdout,'(a)') '************************************************&
                              &********************************'
          write(stdout,*)

          qpole_elec_p = qpole_elec_p * au_to_DAsq
          qpole_ionic_p = qpole_ionic_p * au_to_DAsq
          qpole_total_p = qpole_total_p * au_to_DAsq
          qpole_elec_t = qpole_elec_t * au_to_DAsq
          qpole_ionic_t = qpole_ionic_t * au_to_DAsq
          qpole_total_t = qpole_total_t * au_to_DAsq
          qpole_elec_tg = qpole_elec_tg * au_to_DAsq
          qpole_ionic_tg = qpole_ionic_tg * au_to_DAsq
          qpole_total_tg = qpole_total_tg * au_to_DAsq
          qpole_elec_20 = qpole_elec_20 * au_to_DAsq
          qpole_elec_21c = qpole_elec_21c * au_to_DAsq
          qpole_elec_21s = qpole_elec_21s * au_to_DAsq
          qpole_elec_22c = qpole_elec_22c * au_to_DAsq
          qpole_elec_22s = qpole_elec_22s * au_to_DAsq
          qpole_ionic_20 = qpole_ionic_20 * au_to_DAsq
          qpole_ionic_21c = qpole_ionic_21c * au_to_DAsq
          qpole_ionic_21s = qpole_ionic_21s * au_to_DAsq
          qpole_ionic_22c = qpole_ionic_22c * au_to_DAsq
          qpole_ionic_22s = qpole_ionic_22s * au_to_DAsq
          qpole_total_20 = qpole_total_20 * au_to_DAsq
          qpole_total_21c = qpole_total_21c * au_to_DAsq
          qpole_total_21s = qpole_total_21s * au_to_DAsq
          qpole_total_22c = qpole_total_22c * au_to_DAsq
          qpole_total_22s = qpole_total_22s * au_to_DAsq

       end do
       ! ---------------------------------------------------------------------
       ! ---------------------------------------------------------------------

       write(stdout,*)
       write(stdout,'(a)') repeat('=',80)

    endif

    ! Stop timer
    call timer_clock('properties_polarisation_simcell',2)

  end subroutine properties_polarisation_simcell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GCC32: THIS SUBROUTINE IS OUT OF DATE AND NOT USED ANYMORE
!        I have written a new subroutine to replace this
! subroutine properties_polarisation_berry(mdl, rep, &
!      ngwf_basis, eigs_dens, axis, dipole)
!
!   !==========================================================================!
!   ! Calculates dipole moment of a molecule by computing hamiltonian          !
!   ! eigenstates and Berry phases in the limit of Gamma-point only.           !
!   !--------------------------------------------------------------------------!
!   ! Arguments:                                                               !
!   ! rep (input)        : NGWF representation for valence                     !
!   ! ngwf_basis (input) : Function basis type describing the NGWFs            !
!   !--------------------------------------------------------------------------!
!   ! Written by Simon M.-M. Dubois, October 2013                              !
!   !==========================================================================!
!
!   use basis, only : basis_location_func_wrt_cell, basis_ket_start_wrt_fftbox, &
!        basis_point_wrt_box, basis_start_of_box_wrt_cell
!   use comms, only: pub_on_root, pub_my_proc_id
!   use constants, only: DP, stdout, cmplx_1
!   use dense, only: DEM, dense_create, dense_destroy, dense_get_col, &
!        dense_put_col, dense_product, dense_determinant, &
!        dense_convert, dense_transpose
!   use eigenstates, only: eigenstates_calculate
!   use function_basis, only: FUNC_BASIS
!   use geometry, only: POINT, operator(*), operator(+), magnitude, &
!        operator(.DOT.)
!   use integrals, only: integrals_exr
!   use md_ions, only: md_ionic_polarisation, md_elec_polarisation, &
!        md_total_polarisation
!   use model_type, only: MODEL
!   use ngwf_representation, only: NGWF_REP
!   use rundat, only: pub_print_qc, pub_task, pub_num_spins, pub_spin_fac, &
!        pub_num_kpoints, PUB_1K
!   use sparse, only: sparse_create, sparse_get_element, SPAM3, &
!        sparse_destroy, sparse_put_element, sparse_element_exists
!   use timer, only: timer_clock
!   use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
!        utils_qc_print, utils_assert, utils_banner
!
!   ! Arguments
!   type(MODEL), intent(in)      :: mdl
!   type(NGWF_REP), intent(in)   :: rep
!   type(FUNC_BASIS), intent(in) :: ngwf_basis
!   type(DEM), intent(in)        :: eigs_dens(pub_num_spins)
!
!   integer, optional, intent(in) :: axis
!   real(kind=DP), optional, intent(out) :: dipole
!
!   ! Internal variables
!   real(kind=DP), allocatable  :: tmpcol_real(:)
!   complex(kind=DP), allocatable  :: tmpcol_cmplx(:)
!   complex(kind=DP)  :: berry_phase(3), sberry_el, sberry_el_shift
!   type(POINT)    :: qpt(3), r
!   type(POINT)    :: a1, a2, a3
!   real(kind=DP)  :: cbg1, cbg2, cbg3
!ifdef FFTBOX_OLD_POS
!   integer        :: cs1,cs2,cs3
!endif
!   integer        :: bs1,bs2,bs3
!   real(kind=DP)  :: qr(3), coc(3)
!   real(kind=DP)  :: d(3),q(3),t(3)
!   real(kind=DP)  :: norm, fscale
!   type(SPAM3)    :: sberry(3)
!   type(DEM)      :: sberry_dense(3)
!   type(DEM)      :: eigtot, eigtott, sbo, osbo
!   integer        :: icol, jcol, numocc
!   integer        :: xyz, axmin, axmax, is
!   integer        :: loc_ingwf, ingwf, jngwf
!   integer        :: ierr
!    character(len=40) :: filename
!
!   if (pub_debug_on_root) write(stdout,'(a)') &
!        'DEBUG: Entering properties_polarisation_berry'
!
!   ! Start timer
!   call timer_clock('properties_polarisation_berry',1)
!
!   ! jme: KPOINTS_DANGER
!   ! Some parts of this subroutine enforce a single k-point.
!   call utils_assert(pub_num_kpoints == PUB_1K, &
!        'Subroutine properties_polarisation_berry not ready yet for more&
!        & than one k-point.')
!
!   ! ddor: Allocate only one matrix if axis is specified
!   if (present(axis)) then
!      axmin = axis
!      axmax = axis
!   else
!      axmin = 1
!      axmax = 3
!   end if
!   if (pub_on_root) then
!      write(stdout,'(a)') "Properties polarisation_berry... "
!      write(stdout,'(a,i2,1x,i2)') "   axis :  ", axmin, axmax
!   endif
!
!   ! Reciprocal cell vectors
!   qpt(1) = mdl%cell%b1
!   qpt(2) = mdl%cell%b2
!   qpt(3) = mdl%cell%b3
!   if (pub_on_root) then
!      write(stdout,'(a,3(f12.8,1x))') "   cell : ", mdl%cell%a1%X, mdl%cell%a1%Y, mdl%cell%a1%Z
!      write(stdout,'(a,3(f12.8,1x))') "   cell : ", mdl%cell%a2%X, mdl%cell%a2%Y, mdl%cell%a2%Z
!      write(stdout,'(a,3(f12.8,1x))') "   cell : ", mdl%cell%a3%X, mdl%cell%a3%Y, mdl%cell%a3%Z
!      write(stdout,'(a,3(f12.8,1x))') "   qp1  : ", qpt(1)%X, qpt(1)%Y, qpt(1)%Z
!      write(stdout,'(a,3(f12.8,1x))') "   qp2  : ", qpt(2)%X, qpt(2)%Y, qpt(2)%Z
!      write(stdout,'(a,3(f12.8,1x))') "   qp3  : ", qpt(3)%X, qpt(3)%Y, qpt(3)%Z
!   endif
!
!   ! Number of occupied eigenstates
!   ! jme: KPOINTS_DANGER: only valid if there's one k-point
!   numocc  = sum( rep%n_occ(:, PUB_1K) )
!   !if (pub_on_root) &
!   !   write(stdout,'(a,i)') "   numocc : ", numocc
!
!   ! Create the exponential position matrix with NGWF-overlap sparsity
!   do xyz=axmin,axmax
!      call sparse_create(sberry(xyz),rep%ngwf_overlap,iscmplx=.true.)
!   end do
!
!   ! calculate matrix elements < phi_a | exp[-iG.r] | phi_b >
!   call integrals_exr(sberry,rep%ngwfs_on_grid,ngwf_basis, &
!        rep%ngwfs_on_grid,ngwf_basis,mdl%dbl_grid,mdl%cell,mdl%fftbox, &
!        qpt,axis)
!
!   !if (pub_on_root) &
!   !   write(stdout,'(a)') "   +++> sberry computed! "
!
!
!   ! Calculate vectors between grid points
!   a1 = (1.0_DP / mdl%cell%total_pt1) * mdl%cell%a1
!   a2 = (1.0_DP / mdl%cell%total_pt2) * mdl%cell%a2
!   a3 = (1.0_DP / mdl%cell%total_pt3) * mdl%cell%a3
!
!   ! calculate center of charge
!   fscale = 1.0_DP/sum(mdl%elements(:)%ion_charge)
!   coc(1) =  sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%X)*fscale
!   coc(2) =  sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%Y)*fscale
!   coc(3) =  sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%Z)*fscale
!   !if (pub_on_root) &
!   !   write(stdout,'(a,3(f12.8,x))') "   coc : ", coc(1), coc(2), coc(3)
!
!   ! Loop over all ngwfs on this proc
!   do loc_ingwf=1,ngwf_basis%num_on_proc(pub_my_proc_id)
!      ingwf = loc_ingwf + ngwf_basis%first_on_proc(pub_my_proc_id) - 1
!
!ifdef FFTBOX_OLD_POS
!      ! Determine where tightbox begins wrt fftbox
!      call basis_ket_start_wrt_fftbox(bs1,bs2,bs3,mdl%fftbox%total_pt1,&
!           mdl%fftbox%total_pt2,mdl%fftbox%total_pt3,mdl%fftbox)
!
!      ! Find start of tightbox of ingwf
!      call basis_location_func_wrt_cell(cs1,cs2,cs3, &
!           ngwf_basis%tight_boxes(loc_ingwf),mdl%cell)
!
!      ! Find vector to origin of FFTbox
!      r = real(cs1-bs1,kind=DP) * a1 &
!           + real(cs2-bs2,kind=DP) * a2 &
!           + real(cs3-bs3,kind=DP) * a3
!
!else
!      ! Centre of function wrt fftbox in terms of grid spacings
!      call basis_point_wrt_box(cbg1, cbg2, cbg3, &  ! output
!           ngwf_basis%spheres(loc_ingwf)%centre, mdl%fftbox%total_pt1, &
!           mdl%fftbox%total_pt2, mdl%fftbox%total_pt3, mdl%cell, mdl%fftbox)
!
!      ! Start of fftbox wrt cell in terms of grid-point number
!      call basis_start_of_box_wrt_cell(bs1,bs2,bs3, &
!           ngwf_basis%spheres(loc_ingwf)%centre,cbg1,cbg2,cbg3,mdl%cell)
!
!      ! Find vector to origin of FFTbox
!      r = real(bs1-1,kind=DP) * a1 &
!           + real(bs2-1,kind=DP) * a2 &
!           + real(bs3-1,kind=DP) * a3
!endif
!
!      do xyz=axmin,axmax
!         qr(xyz) = qpt(xyz).DOT.r
!         qr(xyz) = qr(xyz) - coc(xyz)
!      enddo
!
!      ! Loop over all row ngwfs
!      do jngwf=1,ngwf_basis%num
!
!         ! Test if these NGWFs overlap
!         if (.not.sparse_element_exists(rep%ngwf_overlap,jngwf,ingwf)) cycle
!
!         ! Extract element from sberry and shift it
!         do xyz=axmin,axmax
!            call sparse_get_element(sberry_el,sberry(xyz),jngwf,ingwf)
!            sberry_el_shift = exp(cmplx(0.0_DP,qr(xyz),kind=DP)) * sberry_el
!            call sparse_put_element(sberry_el_shift,sberry(xyz),jngwf,ingwf)
!         end do
!
!      end do
!   end do
!
!   !if (pub_on_root) &
!   !   write(stdout,'(a)') "   +++> sberry shifted! "
!
!   ! convert sparse sberry to dense format
!   do xyz=axmin,axmax
!      call dense_create(sberry_dense(xyz),ngwf_basis%num,ngwf_basis%num,iscmplx=.true.)
!      call dense_convert(sberry_dense(xyz),sberry(xyz))
!      call sparse_destroy(sberry(xyz))
!
!      !if (pub_on_root) then
!      !   write(filename,'(a,i1,a)') 'sberry',xyz,'.dense'
!      !   call dense_write(sberry_dense(xyz),trim(adjustl(filename)))
!      !endif
!
!   end do
!
!   !if (pub_on_root) &
!   !   write(stdout,'(a)') "   +++> hamiltonian eigenvectors! "
!
!   ! Keep only the valence eigenvectors regardless of their spin
!   call dense_create(eigtot,ngwf_basis%num,numocc,iscmplx=.true.)
!   allocate(tmpcol_real(ngwf_basis%num),stat=ierr)
!   call utils_alloc_check('properties_polarisation_berry','tmpcol_real',ierr)
!   allocate(tmpcol_cmplx(ngwf_basis%num),stat=ierr)
!   call utils_alloc_check('properties_polarisation_berry','tmpcol_cmplx',ierr)
!   jcol = 0
!   do is = 1, pub_num_spins
!      do icol = 1, rep%n_occ(is,PUB_1K)
!         jcol = jcol + 1
!         call dense_get_col(tmpcol_real,eigs_dens(is),icol)
!         tmpcol_cmplx(:) = cmplx_1*tmpcol_real(:)
!         call dense_put_col(tmpcol_cmplx,eigtot,jcol)
!      enddo
!   enddo
!   deallocate(tmpcol_cmplx,stat=ierr)
!   call utils_dealloc_check('properties_polarisation_berry','tmpcol_cmplx',ierr)
!   deallocate(tmpcol_real,stat=ierr)
!   call utils_dealloc_check('properties_polarisation_berry','tmpcol_real',ierr)
!
!   call dense_create(eigtott,numocc,ngwf_basis%num,iscmplx=.true.)
!   call dense_transpose(eigtott,eigtot)
!
!   ! Calculate Berry phase
!   call dense_create(sbo,ngwf_basis%num,numocc,iscmplx=.true.)
!   call dense_create(osbo,numocc,numocc,iscmplx=.true.)
!   do xyz=axmin,axmax
!      call dense_product(sbo,sberry_dense(xyz),eigtot)
!      !call dense_product(osbo,eigtot,sbo,transpose_amat=.true.)
!      call dense_product(osbo,eigtott,sbo)
!
!      call dense_determinant(osbo,berry_phase(xyz))
!
!   enddo
!   call dense_destroy(osbo)
!   call dense_destroy(sbo)
!   call dense_destroy(eigtott)
!   call dense_destroy(eigtot)
!   do xyz=axmin,axmax
!      call dense_destroy(sberry_dense(xyz))
!   end do
!
!   ! Calculate dipole moment from electrons by Im{ln(berry_phase)}
!   d = 0.0_DP
!
!   ! ddor: Calculate for one direction only if axis is specified
!   axiscase: if (present(axis)) then
!
!      norm = magnitude(qpt(axis))
!      d(axis) = d(axis) - aimag(log(berry_phase(axis)))/norm
!
!      ! multiplication factor depending on number of spins
!      if (pub_num_spins /= 2) d(:) = pub_spin_fac * d(:)
!
!      ! calculate dipole moment from ions
!      q(1) = sum(mdl%elements(:)%ion_charge * (mdl%elements(:)%centre%X - coc(1)))
!      q(2) = sum(mdl%elements(:)%ion_charge * (mdl%elements(:)%centre%Y - coc(2)))
!      q(3) = sum(mdl%elements(:)%ion_charge * (mdl%elements(:)%centre%Z - coc(3)))
!
!      ! calculate total dipole moment
!      t = q + d
!
!      ! print pretty table
!      if (pub_on_root) then
!         write(stdout,'(/a)') utils_banner('=', 'Dipole Moment Calculation')
!         write(stdout,'(a34,4x,a)') 'Type:','Single Molecule'
!         write(stdout,'(a34,4x,a2,1x,i1,1x,a3,1x,f24.18)') 'Electronic dipole moment (e.bohr):',&
!              &'d(',axis,') =',d(axis)
!         write(stdout,'(a34,4x,a2,1x,i1,1x,a3,1x,f24.18)') 'Ionic dipole moment (e.bohr):',&
!              &'q(',axis,') =',q(axis)
!         write(stdout,'(a34,4x,a2,1x,i1,1x,a3,1x,f24.18)') 'Total dipole moment (e.bohr):',&
!              &'t(',axis,') =',t(axis)
!         write(stdout,'(a)') repeat('=',80)
!      end if
!
!      dipole = d(axis)
!
!   else
!
!      do xyz=1,3
!         norm = magnitude(qpt(xyz))
!         d(xyz) = d(xyz) - aimag(log(berry_phase(xyz)))/norm
!      end do
!
!      ! Multiplication factor depending on number of spins
!      if (pub_num_spins /= 2) d(:) = pub_spin_fac * d(:)
!
!      ! Calculate dipole moment from ions
!      q(1) = sum(mdl%elements(:)%ion_charge * (mdl%elements(:)%centre%X- coc(1)))
!      q(2) = sum(mdl%elements(:)%ion_charge * (mdl%elements(:)%centre%Y- coc(2)))
!      q(3) = sum(mdl%elements(:)%ion_charge * (mdl%elements(:)%centre%Z- coc(3)))
!
!      ! Calculate total dipole moment
!      t = q + d
!
!      ! Print pretty table
!      if(pub_on_root) then
!         write(stdout,'(/a)') utils_banner('=', 'Dipole Moment Calculation')
!         write(stdout,'(a34,4x,a)') 'Type:','Single Molecule'
!         write(stdout,*)
!         write(stdout,'(a34,4x,a4,1x,f24.12)') 'Electronic dipole moment &
!              &(e.bohr):','dx =',d(1)
!         write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',d(2)
!         write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',d(3)
!         write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
!              sqrt(sum(d(:)**2))
!         write(stdout,*)
!         write(stdout,'(a34,4x,a4,1x,f24.12)') 'Ionic dipole moment (e.bohr):', &
!              'dx =',q(1)
!         write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',q(2)
!         write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',q(3)
!         write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
!              sqrt(sum(q(:)**2))
!         write(stdout,*)
!         write(stdout,'(a34,4x,a4,1x,f24.12)') 'Total dipole moment (e.bohr):', &
!              'dx =',t(1)
!         write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',t(2)
!         write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',t(3)
!         write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
!              sqrt(sum(t(:)**2))
!         write(stdout,*)
!         if (pub_print_qc) then
!            call utils_qc_print('elec_polarisation_d1',d(1))
!            call utils_qc_print('elec_polarisation_d2',d(2))
!            call utils_qc_print('elec_polarisation_d3',d(3))
!            call utils_qc_print('ionic_polarisation_d1',q(1))
!            call utils_qc_print('ionic_polarisation_d2',q(2))
!            call utils_qc_print('ionic_polarisation_d3',q(3))
!            call utils_qc_print('total_polarisation_d1',t(1))
!            call utils_qc_print('total_polarisation_d2',t(2))
!            call utils_qc_print('total_polarisation_d3',t(3))
!         end if
!         write(stdout,'(a)') repeat('=',80)
!
!         ! smmd: print data in md file
!         if (pub_task == 'MOLECULARDYNAMICS') then
!            md_elec_polarisation(:) = d(:)
!            md_ionic_polarisation(:) = q(:)
!            md_total_polarisation(:) = t(:)
!         end if
!
!      endif
!
!   end if axiscase
!
!   ! Stop timer
!   call timer_clock('properties_polarisation_berry',2)
!
!   if (pub_debug_on_root) write(stdout,'(a)') &
!        'DEBUG: Leaving properties_polarisation_berry'
!
! end subroutine properties_polarisation_berry
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine properties_polarisation(rep,ngwf_basis,proj_basis,mdl, &
       denskern,axis,dipole,dipole_mat)

    !==========================================================================!
    ! Calculates dipole moment of a molecule by computing matrix elements      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! rep (input)        : NGWF representation for valence                     !
    ! ngwf_basis (input) : Function basis type describing the NGWFs            !
    ! denskern (input)   : Density kernel                                      !
    !--------------------------------------------------------------------------!
    ! Written by Mark Robinson, April 2008                                     !
    ! Modified by David O'Regan in March 2009 so that it can be                !
    ! called for a single direction only.                                      !
    ! Extended for PAW sphere terms by Nicholas Hine in November 2011.         !
    ! lpl: Modified to allow passing of dipole matrices outside subroutine     !
    !      17-04-2012                                                          !
    ! Polarisation density output added by David O'Regan in July 2013.         !
    ! Modified for embedding and to remove pub_par by Joseph Prentice, August  !
    ! 2018                                                                     !
    !==========================================================================!

    use basis, only : basis_location_func_wrt_cell, &
         basis_point_wrt_box, basis_start_of_box_wrt_cell
    use comms, only: pub_on_root, pub_my_proc_id, comms_reduce
    use constants, only: DP, stdout, ANGSTROM, UP, DN, max_spins
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(*), operator(+)
    use optics, only: optics_pos_mat_els
    use md_ions, only: md_ionic_polarisation, md_elec_polarisation, &
         md_total_polarisation
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use paw, only: paw_position_operator
    use rundat, only: pub_print_qc, pub_write_polarisation_plot, pub_task, &
         pub_num_spins, pub_spin_fac, pub_polarisation_local, &
         pub_num_kpoints, PUB_1K
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create, &
         sparse_embed_trace, sparse_embed_destroy, sparse_embed_copy, &
         sparse_embed_product, sparse_embed_transpose, sparse_embed_axpy, &
         sparse_embed_scale, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use sparse, only: SPAM3, sparse_get_element
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_qc_print, utils_assert, utils_banner
    use visual, only: visual_scalarfield

    ! Arguments
    type(NGWF_REP), intent(in)   :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(MODEL), intent(in)      :: mdl
    type(SPAM3_EMBED_ARRAY), intent(in):: denskern
    ! ddor: When axis is present, only the electronic contribution to the dipole
    !       moment in direction axis is calculated, and this is output to dipole.
    integer, optional, intent(in) :: axis
    real(kind=DP), optional, intent(out) :: dipole

    ! lpl: Accepts external dipole_mat buffer, and returns dipole mat
    type(SPAM3_EMBED), optional, intent(inout) :: dipole_mat(3)

    ! Local Variables
    type(SPAM3_EMBED) :: r_elements(3)
    type(SPAM3_EMBED) :: kr_elements, srk_elements, krs_elements(max_spins)
    type(SPAM3), allocatable :: kern_array(:)
    real(kind=DP) :: d(3),q(3),t(3), dtot(3),qtot(3),ttot(3),&
                     dsol(3),qsol(3),tsol(3),doff(3),qoff(3),toff(3)
#ifdef FFTBOX_OLD_POS
    integer       :: cs1,cs2,cs3
#endif
    real(kind=DP) :: xji,kij,sji,trace_el,trsol_el,trace_nu,trsol_nu,&
                     traceoff_el,trace_all,trace
    integer       :: xyz,axmin,axmax,is
    integer       :: jngwf,loc_ingwf,ingwf
    integer       :: ierr ! memory allocation error flag
    integer       :: iat,isub,jsub
    integer :: fine_max_slabs12, fine_ld1, fine_ld2 ! jd: shortcut notation
    character(len=1) :: dir_string

    ! Allocatable local variables
    real(kind=DP), allocatable :: polarisation_density_fine(:,:,:,:)
    real(kind=DP), allocatable :: buffer_fine(:,:,:,:)

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering properties_polarisation'

    ! Start timer
    call timer_clock('properties_polarisation',1)

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine properties_polarisation not ready yet for more&
         & than one k-point.')

    ! ddor: Allocate only one matrix if axis is specified
    if (present(axis)) then
       axmin = axis
       axmax = axis
    else
       axmin = 1
       axmax = 3
    end if

    ! ndmh: create position elements matrix with NGWF-overlap sparsity
    do xyz=axmin,axmax
       call sparse_embed_create(r_elements(xyz),rep%ngwf_overlap)
    end do

    ! Calculate matrix elements of position operator in NGWF representation
    call optics_pos_mat_els(r_elements,rep,ngwf_basis,rep,ngwf_basis,&
         rep%overlap,rep%ngwf_overlap,proj_basis, mdl)

    ! Calculate dipole moment from electrons by Tr(KR)
    d = 0.0_DP

    ! ddor: Calculate for one direction only if axis is specified
    axiscase: if (present(axis)) then

       do is=1,pub_num_spins
          call sparse_embed_trace(trace,denskern%m(is,PUB_1K),r_elements(axis))
          d(axis) = d(axis) - trace
       end do

       ! multiplication factor depending on number of spins
       if (pub_num_spins /= 2) d(:) = pub_spin_fac * d(:)

       ! calculate dipole moment from ions
       q(1) = sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%X)
       q(2) = sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%Y)
       q(3) = sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%Z)

       ! calculate total dipole moment
       t = q + d

       ! print pretty table
       if (pub_on_root) then
          write(stdout,'(/a)') utils_banner('=', 'Dipole Moment Calculation')
          write(stdout,'(a34,4x,a)') 'Type:','Single Molecule'
          write(stdout,'(a34,4x,a2,1x,i1,1x,a3,1x,f24.18)') 'Electronic dipole moment (e.bohr):',&
               &'d(',axis,') =',d(axis)
          write(stdout,'(a34,4x,a2,1x,i1,1x,a3,1x,f24.18)') 'Ionic dipole moment (e.bohr):',&
               &'q(',axis,') =',q(axis)
          write(stdout,'(a34,4x,a2,1x,i1,1x,a3,1x,f24.18)') 'Total dipole moment (e.bohr):',&
               &'t(',axis,') =',t(axis)
          write(stdout,'(a)') repeat('=',80)
       end if

       dipole = d(axis)

    ! vv : Calculate the partial dipole moment by manually multipling
    ! the matrix elements
    elseif (pub_polarisation_local) then
       ! jcap: loop over regions
       do isub=1,mdl%nsub
          do loc_ingwf=1,ngwf_basis(isub)%num_on_proc(pub_my_proc_id)
             ingwf = loc_ingwf + ngwf_basis(isub)%first_on_proc(pub_my_proc_id) - 1
             ! jcap: loop over regions
             do jsub=1,mdl%nsub
                do jngwf=1,ngwf_basis(jsub)%num
                   if(mdl%regions(jsub)%elements(mdl%regions(jsub)%par%orig_atom(&
                        ngwf_basis(jsub)%atom_of_func(jngwf)))%loc_dipole &
                        .AND. mdl%regions(isub)%elements(mdl%regions(isub)%par%orig_atom(&
                        ngwf_basis(isub)%atom_of_func(ingwf)))&
                        %loc_dipole) then
                      do xyz=1,3
                         call sparse_get_element(xji,r_elements(xyz)%m(jsub,isub),jngwf,ingwf)
                         do is=1,pub_num_spins
                            call sparse_get_element(kij,denskern%m(is,PUB_1K)%m(jsub,isub),&
                                 jngwf,ingwf)
                            d(xyz) = d(xyz) - kij * xji
                         end do
                      end do
                   end if
                end do
             end do
          end do
       end do

       call comms_reduce('SUM',d)

       ! Multiplication factor depending on number of spins
       if (pub_num_spins /= 2) d(:) = pub_spin_fac * d(:)
       ! vv: Calculate dipole moment from ions that belong to the
       ! loc_dipole_group
       q(:) = 0.0_DP
       do iat = 1,mdl%nat
          if(mdl%elements(iat)%loc_dipole) then
             q(1) = q(1) + mdl%elements(iat)%ion_charge * mdl%elements(iat)%centre%X
             q(2) = q(2) + mdl%elements(iat)%ion_charge * mdl%elements(iat)%centre%Y
             q(3) = q(3) + mdl%elements(iat)%ion_charge * mdl%elements(iat)%centre%Z
          end if
       end do

       ! vv: off-diagonal contribution
       ! Calculate the dipole moment of the off-diagonal blocks if any
       if (.not. all(mdl%elements%loc_dipole)) then
          doff(:) = 0.0_DP
          toff(:) = 0.0_DP
          qoff(:) = 0.0_DP
          traceoff_el = 0.0_DP
          ! jcap: loop over regions
          do isub=1,mdl%nsub
             do loc_ingwf=1,ngwf_basis(isub)%num_on_proc(pub_my_proc_id)
                ingwf = loc_ingwf + ngwf_basis(isub)%first_on_proc(pub_my_proc_id) - 1
                ! jcap: loop over regions
                do jsub=1,mdl%nsub
                   do jngwf=1,ngwf_basis(jsub)%num
                      if((mdl%regions(jsub)%elements(mdl%regions(jsub)%par%orig_atom(&
                           ngwf_basis(jsub)%atom_of_func(jngwf)))%loc_dipole &
                           .and. .not.(mdl%regions(isub)%elements(mdl%regions(isub)%par%orig_atom(&
                           ngwf_basis(isub)%atom_of_func(ingwf)))&
                           %loc_dipole)) .or.(.not.(mdl%regions(jsub)%elements(&
                           mdl%regions(jsub)%par%orig_atom(ngwf_basis(jsub)%atom_of_func(jngwf)))%loc_dipole)&
                           .and. mdl%regions(isub)%elements(mdl%regions(isub)%par%orig_atom(&
                           ngwf_basis(isub)%atom_of_func(ingwf)))%loc_dipole)) then
                         do xyz=1,3
                            call sparse_get_element(xji,r_elements(xyz)%m(jsub,isub),jngwf,ingwf)
                            do is=1,pub_num_spins
                               call sparse_get_element(kij, &
                                    denskern%m(is,PUB_1K)%m(jsub,isub),jngwf,ingwf)
                               doff(xyz) = doff(xyz) - kij * xji
                            end do
                         end do
                         call sparse_get_element(sji,rep%overlap%m(jsub,isub),jngwf,ingwf)
                         do is=1,pub_num_spins
                            call sparse_get_element(kij, &
                                 denskern%m(is,PUB_1K)%m(jsub,isub),jngwf,ingwf)
                            traceoff_el = traceoff_el + kij*sji
                         end do
                      end if
                   end do
                end do
             end do
          end do

          call comms_reduce('SUM',traceoff_el)
          call comms_reduce('SUM',doff)

          ! Multiplication factor depending on number of spins
          if (pub_num_spins /= 2) then
             doff(:) = pub_spin_fac * doff(:)
             traceoff_el = pub_spin_fac * traceoff_el
          end if

          ! vv: Calculate dipole moment from ions that belong to the
          ! loc_dipole_group
          do iat = 1,mdl%nat
             if(.not. mdl%elements(iat)%loc_dipole) then
                qoff(1) = qoff(1) + mdl%elements(iat)%ion_charge * mdl%elements(iat)%centre%X
                qoff(2) = qoff(2) + mdl%elements(iat)%ion_charge * mdl%elements(iat)%centre%Y
                qoff(3) = qoff(3) + mdl%elements(iat)%ion_charge * mdl%elements(iat)%centre%Z
             end if
          end do

          ! Calculate the off-diagonal total dipole moment
          toff = qoff + doff
       end if

       ! Apply correction
       d(:) = d(:) + 0.5_dp * doff(:)


       ! Calculate the dipole moment of the solvent if any
       if (.not. all(mdl%elements%loc_dipole)) then
          dsol(:) = 0.0_DP
          tsol(:) = 0.0_DP
          qsol(:) = 0.0_DP
          ! jcap: loop over regions
          do isub=1,mdl%nsub
             do loc_ingwf=1,ngwf_basis(isub)%num_on_proc(pub_my_proc_id)
                ingwf = loc_ingwf + ngwf_basis(isub)%first_on_proc(pub_my_proc_id) - 1
                ! jcap: loop over regions
                do jsub=1,mdl%nsub
                   do jngwf=1,ngwf_basis(jsub)%num
                      if(.not.(mdl%regions(jsub)%elements(mdl%regions(jsub)%par%orig_atom(&
                           ngwf_basis(jsub)%atom_of_func(jngwf)))%loc_dipole) &
                           .and. .not.(mdl%regions(isub)%elements(mdl%regions(isub)%par%orig_atom(&
                           ngwf_basis(isub)%atom_of_func(ingwf)))&
                           %loc_dipole)) then
                         do xyz=1,3
                            call sparse_get_element(xji,r_elements(xyz)%m(jsub,isub),jngwf,ingwf)
                            do is=1,pub_num_spins
                               call sparse_get_element(kij, &
                                    denskern%m(is,PUB_1K)%m(jsub,isub),jngwf,ingwf)
                               dsol(xyz) = dsol(xyz) - kij * xji
                            end do
                         end do
                      end if
                   end do
                end do
             end do
          end do

          call comms_reduce('SUM',dsol)

          ! Multiplication factor depending on number of spins
          if (pub_num_spins /= 2) dsol(:) = pub_spin_fac * dsol(:)

          ! Apply correction
          dsol(:) = dsol(:) + 0.5_dp * doff(:)

          ! vv: Calculate dipole moment from ions that do not belong to the
          ! loc_dipole_group
          qsol(:) = 0.0_DP
          do iat = 1,mdl%nat
             if(.not. mdl%elements(iat)%loc_dipole) then
                qsol(1) = qsol(1) + mdl%elements(iat)%ion_charge * mdl%elements(iat)%centre%X
                qsol(2) = qsol(2) + mdl%elements(iat)%ion_charge * mdl%elements(iat)%centre%Y
                qsol(3) = qsol(3) + mdl%elements(iat)%ion_charge * mdl%elements(iat)%centre%Z
             end if
          end do

       end if



       ! Calculate the electronic charge of the solute, or more rigorously, the Tr(KS)
       trace_el = 0.0_DP

       ! Loop over NGWFs on this proc
       ! jcap: loop over regions
       do isub=1,mdl%nsub
          do loc_ingwf=1,ngwf_basis(isub)%num_on_proc(pub_my_proc_id)
             ingwf = loc_ingwf + ngwf_basis(isub)%first_on_proc(pub_my_proc_id) - 1
             ! Loop over all rows
             ! jcap: loop over regions
             do jsub=1,mdl%nsub
                do jngwf=1,ngwf_basis(jsub)%num
                   if(mdl%regions(jsub)%elements(mdl%regions(jsub)%par%orig_atom(&
                        ngwf_basis(jsub)%atom_of_func(jngwf)))%loc_dipole &
                        .AND. mdl%regions(isub)%elements(mdl%regions(isub)%par%orig_atom(&
                        ngwf_basis(isub)%atom_of_func(ingwf)))&
                        %loc_dipole) then
                      call sparse_get_element(sji,rep%overlap%m(jsub,isub),jngwf,ingwf)
                      do is = 1,pub_num_spins
                         call sparse_get_element(kij, &
                              denskern%m(is,PUB_1K)%m(jsub,isub),jngwf,ingwf)
                         trace_el = trace_el + kij*sji
                      end do
                   end if
                end do
             end do
          end do
       end do

       call comms_reduce('SUM',trace_el)
       ! Multiplication factor depending on number of spins
       if (pub_num_spins /= 2) trace_el = pub_spin_fac * trace_el


       ! vv: Calculate the ionic charge of the solute
       trace_nu = 0.0_DP
       do iat = 1,mdl%nat
          if(mdl%elements(iat)%loc_dipole) then
             trace_nu = trace_nu + mdl%elements(iat)%ion_charge
          end if
       end do


       trsol_el = 0.0_DP

       ! Loop over NGWFs on this proc
       ! jcap: loop over regions
       do isub=1,mdl%nsub
          do loc_ingwf=1,ngwf_basis(isub)%num_on_proc(pub_my_proc_id)
             ingwf = loc_ingwf + ngwf_basis(isub)%first_on_proc(pub_my_proc_id) - 1
             ! Loop over all rows
             ! jcap: loop over regions
             do jsub=1,mdl%nsub
                do jngwf=1,ngwf_basis(jsub)%num
                   if(.not.(mdl%regions(jsub)%elements(mdl%regions(jsub)%par%orig_atom(&
                        ngwf_basis(jsub)%atom_of_func(jngwf)))%loc_dipole) &
                        .and. .not.(mdl%regions(isub)%elements(mdl%regions(isub)%par%orig_atom(&
                        ngwf_basis(isub)%atom_of_func(ingwf)))&
                        %loc_dipole)) then
                      call sparse_get_element(sji,rep%overlap%m(jsub,isub),jngwf,ingwf)
                      do is = 1,pub_num_spins
                         call sparse_get_element(kij, &
                              denskern%m(is,PUB_1K)%m(jsub,isub),jngwf,ingwf)
                         trsol_el = trsol_el + kij*sji
                      end do
                   end if
                end do
             end do
          end do
       end do

       call comms_reduce('SUM',trsol_el)
       ! Multiplication factor depending on number of spins
       if (pub_num_spins /= 2) trsol_el = pub_spin_fac * trsol_el

       ! vv: Calculate dipole moment from ions that belong to the
       ! loc_dipole_group
       trsol_nu = 0.0_DP
       do iat = 1,mdl%nat
          if(mdl%elements(iat)%loc_dipole) then
             trsol_nu = trsol_nu + mdl%elements(iat)%ion_charge
          end if
       end do

       ! Apply correction to the electronic dipole moments in order to
       ! obtain a consistent charge, i.e. integer number of electrons

       trace_el = (trace_el + 0.5_dp * traceoff_el)
       trsol_el = (trsol_el + 0.5_dp * traceoff_el)

       d(:) = d(:) * (trace_nu/trace_el)
       dsol(:) = dsol(:) * (trsol_nu / trsol_el)

       ! Calculate the solute total dipole moment with the correction
       t = q + d

       ! Calculate the solvent total dipole moment with correction
       tsol = qsol + dsol

       ! Calculate the electronic dipole moment of the entire system as well
       dtot(:) = 0.0_DP
       do xyz=1,3
          do is=1,pub_num_spins
             call sparse_embed_trace(trace,denskern%m(is,PUB_1K),r_elements(xyz))
             dtot(xyz) = dtot(xyz) - trace
          end do
       end do

       ! Multiplication factor depending on number of spins
       if (pub_num_spins /= 2) dtot(:) = pub_spin_fac * dtot(:)

       ! Calculate dipole moment from ions
       qtot(1) = sum( mdl%elements(:)%ion_charge * mdl%elements(:)%centre%X)
       qtot(2) = sum( mdl%elements(:)%ion_charge * mdl%elements(:)%centre%Y)
       qtot(3) = sum( mdl%elements(:)%ion_charge * mdl%elements(:)%centre%Z)

       ! Calculate total dipole moment
       ttot = qtot + dtot

       ! Calculate the electronic charge of the entire system
       trace_all = 0.0_DP
       do is=1,pub_num_spins
          call sparse_embed_trace(trace,denskern%m(is,PUB_1K),rep%overlap)
          trace_all = trace_all + trace
       end do

       ! Multiplication factor depending on number of spins
       if (pub_num_spins /= 2) trace_all = pub_spin_fac * trace_all

       ! Print pretty table
       if(pub_on_root) then
          write(stdout,'(/a)') utils_banner('=','Dipole Moment Calculation')
          write(stdout,'(a34,4x,a)') 'Type:','Selected subsystem'
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Electronic dipole moment &
               &(e.bohr):','dx =',d(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',d(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',d(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(d(:)**2))
          write(stdout,*)
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Ionic dipole moment (e.bohr):', &
               'qx =',q(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'qy =',q(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'qz =',q(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(q(:)**2))
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Total solute dipole moment (e.bohr):', &
               'tx =',t(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'ty =',t(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'tz =',t(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(t(:)**2))
          write(stdout,*)
          write(stdout,'(a34,4x,a)') 'Type:','Complementary subsystem'
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Electronic dipole moment &
               &(e.bohr):','dx =',dsol(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',dsol(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',dsol(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(dsol(:)**2))
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Ionic dipole moment (e.bohr):', &
               'qx =',qsol(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'qy =',qsol(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'qz =',qsol(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(qsol(:)**2))
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Total solute dipole moment (e.bohr):', &
               'tx =',tsol(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'ty =',tsol(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'tz =',tsol(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(tsol(:)**2))
          write(stdout,*)
          write(stdout,'(a34,4x,a)') 'Type:','Entire system'
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Electronic dipole moment &
               &(e.bohr):','dx =',dtot(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',dtot(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',dtot(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(dtot(:)**2))
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Ionic dipole moment (e.bohr):', &
               'qx =',qtot(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'qy =',qtot(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'qz =',qtot(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(qtot(:)**2))
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Total solute dipole moment (e.bohr):', &
               'tx =',ttot(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'ty =',ttot(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'tz =',ttot(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(ttot(:)**2))
          write(stdout,*)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Ionic charge of selected atoms (e.bohr):', &
               trace_nu
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Electronic charge of selected atoms (e.bohr):', &
               trace_el
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Total Trace(KS) :', &
               trace_all
          write(stdout,*)
          if (pub_print_qc) then
             call utils_qc_print('elec_polarisation_d1',d(1))
             call utils_qc_print('elec_polarisation_d2',d(2))
             call utils_qc_print('elec_polarisation_d3',d(3))
             call utils_qc_print('ionic_polarisation_d1',q(1))
             call utils_qc_print('ionic_polarisation_d2',q(2))
             call utils_qc_print('ionic_polarisation_d3',q(3))
             call utils_qc_print('total_polarisation_d1',t(1))
             call utils_qc_print('total_polarisation_d2',t(2))
             call utils_qc_print('total_polarisation_d3',t(3))
          end if
          write(stdout,'(a)') repeat('=',80)
          ! vv: print data in md file
          if (pub_task == 'MOLECULARDYNAMICS') then
             md_elec_polarisation(:) = d(:)
             md_ionic_polarisation(:) = q(:)
             md_total_polarisation(:) = t(:)
          end if

       end if

    else

       do xyz=1,3

          ! ddor: Write out the electronic "polarisation density" defined as
          !     : -1.0 * phi_a(r) K^ab <r>^i_bc S^cd phi_d(r) to .cube file
          !     : N.B. only the change in this quantity between different
          !     : densities is of any interest, and NOT the quantity itself.
          if (pub_write_polarisation_plot) then

             ! ddor: Allocate polarisation density slabs
             fine_ld1 = mdl%fine_grid%ld1
             fine_ld2 = mdl%fine_grid%ld2
             fine_max_slabs12 = mdl%fine_grid%max_slabs12
             allocate(polarisation_density_fine(fine_ld1,fine_ld2,&
                  fine_max_slabs12,pub_num_spins),stat=ierr)
             call utils_alloc_check('properties_polarisation',&
                  'polarisation_density_fine',ierr)

             ! ddor: Create matrices with sparsity KSK for symmetrisation
             call sparse_embed_create(kr_elements,denskern%m(1,1),r_elements(xyz))
             do is=1,pub_num_spins
                call sparse_embed_create(krs_elements(is),kr_elements,denskern%m(1,1))
             enddo
             call sparse_embed_create(srk_elements,krs_elements(1))

             ! ddor: Loop over spins
             do is=1,pub_num_spins

                ! ddor: Compute the dipole moment
                call sparse_embed_product(kr_elements,denskern%m(is,PUB_1K),r_elements(xyz))
                call sparse_embed_trace(trace,kr_elements)
                d(xyz) = d(xyz) - trace

                ! ddor: Compute the symmetrised matrix elements of K<r>S^^
                call sparse_embed_product(krs_elements(is),kr_elements,&
                     rep%inv_overlap)
                call sparse_embed_transpose(srk_elements,krs_elements(is))
                call sparse_embed_axpy(krs_elements(is),srk_elements,1.0_DP)
                if (pub_num_spins .ne. 1) & ! Scale for spin if necessary
                     call sparse_embed_scale(krs_elements(is),0.5_DP)

             enddo

             ! jcap: temporary arrays for calculating density
             allocate(buffer_fine(fine_ld1,fine_ld2,&
                  fine_max_slabs12,pub_num_spins),stat=ierr)
             call utils_alloc_check('properties_polarisation',&
                  'buffer_fine',ierr)
             allocate(kern_array(pub_num_spins),stat=ierr)
             call utils_alloc_check('properties_polarisation',&
                  'kern_array',ierr)

             polarisation_density_fine=0.0_DP
             ! jcap: loop over regions
             do isub=1,mdl%nsub
                do jsub=1,mdl%nsub
                   ! jcap: create temporary kernel array
                   call sparse_embed_extract_from_array(kern_array,krs_elements,&
                        isub,jsub)
                   ! ddor: Calculate data-parallelised polarisation density
                   call density_on_grid(buffer_fine,mdl%fine_grid, &
                        mdl%dbl_grid,mdl%cell, mdl%fftbox, kern_array,&
                        rep%ngwf_overlap%m(isub,jsub),&
                        rep%ngwfs_on_grid(isub),ngwf_basis(isub), &
                        rep%ngwfs_on_grid(jsub),ngwf_basis(jsub))

                   polarisation_density_fine=polarisation_density_fine+buffer_fine

                   ! jcap: destroy temporary kernel array
                   call sparse_embed_destroy_extracted_array(kern_array)
                end do
             end do


             ! ddor: spin polarisation: calculate total polarisation density
             !     : and polarisation spin density
             if (pub_num_spins == 2) then
                polarisation_density_fine(:,:,:,1) = &
                     polarisation_density_fine(:,:,:,UP) + &
                     polarisation_density_fine(:,:,:,DN)
                polarisation_density_fine(:,:,:,2) = &
                     polarisation_density_fine(:,:,:,UP) - &
                     2.0_DP * polarisation_density_fine(:,:,:,DN)
             end if

             ! ddor: Direction string for file name
             if (xyz==1) dir_string='x'
             if (xyz==2) dir_string='y'
             if (xyz==3) dir_string='z'

             ! cks: output density in plot format file
             ! vm: output density in Angstrom rather than in Bohr
             ! jd: ... but leave the unit conversion to visual_scalarfield
             ! ddor: The minus sign for electrons occurs here
             call visual_scalarfield( &
                  polarisation_density_fine(:,:,:,1), mdl%fine_grid, mdl%cell, &
                  'Electronic polarisation density (in e/ang^2) for:', &
                  '_polarisation_density_'//trim(dir_string), &
                  mdl%elements, -1.0_DP * ANGSTROM**2)
             if (pub_num_spins == 2) call visual_scalarfield( &
                  polarisation_density_fine(:,:,:,2), mdl%fine_grid, mdl%cell, &
                  'Electronic polarisation spin density (in e/ang^2) for:', &
                  '_polarisation_spindensity_'//trim(dir_string), &
                  mdl%elements, -1.0_DP * ANGSTROM**2)

             !ddor: Destroy matrices
             call sparse_embed_destroy(srk_elements)
             do is=pub_num_spins,1,-1
                call sparse_embed_destroy(krs_elements(is))
             enddo
             call sparse_embed_destroy(kr_elements)

             ! ddor: Deallocate polarisation density slabs
             deallocate(polarisation_density_fine,stat=ierr)
             call utils_dealloc_check('properties_polarisation',&
                  'polarisation_density_fine',ierr)
             deallocate(buffer_fine,stat=ierr)
             call utils_dealloc_check('properties_polarisation',&
                  'buffer_fine',ierr)
             deallocate(kern_array,stat=ierr)
             call utils_dealloc_check('properties_polarisation',&
                  'kern_array',ierr)

          else

             do is=1,pub_num_spins
                call sparse_embed_trace(trace,denskern%m(is,PUB_1K),r_elements(xyz))
                d(xyz) = d(xyz) - trace
             end do

          endif

       end do

       ! Multiplication factor depending on number of spins
       if (pub_num_spins /= 2) d(:) = pub_spin_fac * d(:)

       ! Calculate dipole moment from ions
       q(1) = sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%X)
       q(2) = sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%Y)
       q(3) = sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%Z)

       ! Calculate total dipole moment
       t = q + d

       ! Print pretty table
       if(pub_on_root) then
          write(stdout,'(/a)') utils_banner('=', 'Dipole Moment Calculation')
          write(stdout,'(a34,4x,a)') 'Type:','Single Molecule'
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Electronic dipole moment &
               &(e.bohr):','dx =',d(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',d(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',d(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(d(:)**2))
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Ionic dipole moment (e.bohr):', &
               'dx =',q(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',q(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',q(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(q(:)**2))
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Total dipole moment (e.bohr):', &
               'dx =',t(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',t(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',t(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(t(:)**2))
          write(stdout,*)
          if (pub_print_qc) then
             call utils_qc_print('elec_polarisation_d1',d(1))
             call utils_qc_print('elec_polarisation_d2',d(2))
             call utils_qc_print('elec_polarisation_d3',d(3))
             call utils_qc_print('ionic_polarisation_d1',q(1))
             call utils_qc_print('ionic_polarisation_d2',q(2))
             call utils_qc_print('ionic_polarisation_d3',q(3))
             call utils_qc_print('total_polarisation_d1',t(1))
             call utils_qc_print('total_polarisation_d2',t(2))
             call utils_qc_print('total_polarisation_d3',t(3))
          end if
          write(stdout,'(a)') repeat('=',80)

          ! smmd: print data in md file
          if (pub_task == 'MOLECULARDYNAMICS') then
             md_elec_polarisation(:) = d(:)
             md_ionic_polarisation(:) = q(:)
             md_total_polarisation(:) = t(:)
          end if

       endif

    end if axiscase

    ! ndmh: copy to output matrix if required
    if (present(dipole_mat)) then
       do xyz=axmin,axmax
          call sparse_embed_copy(dipole_mat(xyz),r_elements(xyz))
       end do
    end if

    ! Deallocate sparse matrices
    do xyz=axmax,axmin,-1
       call sparse_embed_destroy(r_elements(xyz))
    end do

    ! Stop timer
    call timer_clock('properties_polarisation',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving properties_polarisation'

  end subroutine properties_polarisation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine properties_spread(ngwfs_on_grid,ngwf_basis, &
       mdl,inv_overlap,overlap,nonorthogonal)

    !==========================================================================!
    ! Calculates the spread (second central moment) of the NGWFs.              !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! ngwfs_on_grid (input)    : NGWFs in ppd representation on this proc      !
    ! ngwf_basis (input)       : Function basis type describing the NGWFs      !
    ! inv_overlap (input)      : Inverse overlap matrix                        !
    ! overlap (input)          : Overlap matrix                                !
    !--------------------------------------------------------------------------!
    ! Written by David O'Regan in September 2009,                              !
    ! based on properties_polarisation by M. Robinson.                         !
    ! Automatic arrays removed by Nicholas Hine, December 2009.                !
    !==========================================================================!

    use basis, only : basis_location_func_wrt_cell, basis_ket_start_wrt_fftbox, &
         basis_point_wrt_box, basis_start_of_box_wrt_cell
    use comms, only: pub_on_root, pub_my_proc_id, comms_reduce, comms_barrier
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_brappd_ketppd
    use geometry, only: POINT, operator(*), operator(+)
    use integrals, only: integrals_pos
    use model_type, only: MODEL
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_aug
    use sparse, only: sparse_create, sparse_get_element, SPAM3, &
         sparse_destroy, sparse_put_element, sparse_product, &
         sparse_axpy, sparse_element_exists
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(MODEL), intent(in)              :: mdl
    type(FUNCTIONS), intent(in)            :: ngwfs_on_grid
    type(SPAM3), intent(in)              :: inv_overlap
    type(SPAM3), intent(in)              :: overlap
    ! ddor: Whether or not the input NGWFs are nonorthogonal
    logical, intent(in)                  :: nonorthogonal

    ! Local Variables
    type(SPAM3)   :: r_elements(3)
    type(SPAM3)   :: r2_elements(3)
    type(SPAM3)   :: rs(3), rsr(3), thespread(3)
    type(SPAM3)   :: closetounity
    type(POINT)   :: a1,a2,a3
    type(POINT)   :: r
    real(kind=DP) :: R_fft(3),r_el(3),r2_el(3),o_el
    real(kind=DP) :: cbg1, cbg2, cbg3
#ifdef FFTBOX_OLD_POS
    integer       :: cs1,cs2,cs3
#endif
    integer       :: bs1,bs2,bs3
    integer       :: xyz
    integer       :: loc_ingwf, ingwf, jngwf
    integer       :: atom, ngwf_on_atom
    integer       :: ierr
    real(kind=DP), allocatable :: covariant(:,:)
    real(kind=DP), allocatable :: invariant(:,:)

    ! Start timer
    call timer_clock('properties_spread',1)

    ! ndmh: fail if called for augmented matrices
    if (pub_aug) then
       call utils_abort('ERROR in properties_spread: routine does not &
            &support augmentation')
    end if


    ! ndmh: allocate temporary arrays
    allocate(covariant(ngwf_basis%num,4),stat=ierr)
    call utils_alloc_check('properties_spread','covariant',ierr)
    allocate(invariant(ngwf_basis%num,4),stat=ierr)
    call utils_alloc_check('properties_spread','invariant',ierr)

    do xyz=1,3
       call sparse_create(r_elements(xyz),overlap)
       call sparse_create(r2_elements(xyz),overlap)
       if (nonorthogonal) then
          call sparse_create(rs(xyz),r_elements(xyz),inv_overlap)
          call sparse_create(rsr(xyz),rs(xyz),r_elements(xyz))
          call sparse_create(thespread(xyz),rsr(xyz),inv_overlap)
       else
          call sparse_create(thespread(xyz),overlap,overlap)
       endif
    enddo

    if (.not. nonorthogonal) then
       call sparse_create(closetounity,overlap)
       call function_ops_brappd_ketppd(closetounity, &  !input-output
            ngwfs_on_grid, ngwf_basis, &
            ngwfs_on_grid, ngwf_basis, mdl%cell)
    endif

    covariant = 0.0_DP
    invariant = 0.0_DP

    ! Calculate vectors between grid points
    a1 = (1.0_DP / mdl%cell%total_pt1) * mdl%cell%a1
    a2 = (1.0_DP / mdl%cell%total_pt2) * mdl%cell%a2
    a3 = (1.0_DP / mdl%cell%total_pt3) * mdl%cell%a3

    ! calculate matrix elements < phi_a | r_op - R_fft | phi_b >
    ! ndmh: use integrals_pos routine
    call integrals_pos(r_elements,&
         ngwfs_on_grid,ngwf_basis,ngwfs_on_grid,ngwf_basis,mdl%dbl_grid, &
         mdl%cell,mdl%fftbox,1)

    ! calculate matrix elements < phi_a | ( r_op - R_fft )^2 | phi_b >
    call integrals_pos(r2_elements,&
         ngwfs_on_grid,ngwf_basis,ngwfs_on_grid,ngwf_basis,mdl%dbl_grid, &
         mdl%cell,mdl%fftbox,2)

    ! Add R_fft * S_ab to matrix elements, taking into account atoms moving
    ! across cell boundary (NOT IMPLEMENTED YET)
    do loc_ingwf=1,ngwf_basis%num_on_proc(pub_my_proc_id)

#ifdef FFTBOX_OLD_POS
       ! Determine where tightbox begins wrt fftbox
       call basis_ket_start_wrt_fftbox(bs1,bs2,bs3,mdl%fftbox%total_pt1,&
            mdl%fftbox%total_pt2,mdl%fftbox%total_pt3,mdl%fftbox)

       ! Find start of tightbox of ingwf
       call basis_location_func_wrt_cell(cs1,cs2,cs3, &
            ngwf_basis%tight_boxes(loc_ingwf),mdl%cell)

       ! Find vector to origin of FFTbox
       r = real(cs1-bs1,kind=DP) * a1 &
            + real(cs2-bs2,kind=DP) * a2 &
            + real(cs3-bs3,kind=DP) * a3

#else
       ! Centre of function wrt fftbox in terms of grid spacings
       call basis_point_wrt_box(cbg1, cbg2, cbg3, &  ! output
            ngwf_basis%spheres(loc_ingwf)%centre, mdl%fftbox%total_pt1, &
            mdl%fftbox%total_pt2, mdl%fftbox%total_pt3, mdl%cell, mdl%fftbox)

       ! Start of fftbox wrt cell in terms of grid-point number
       call basis_start_of_box_wrt_cell(bs1,bs2,bs3, &
            ngwf_basis%spheres(loc_ingwf)%centre,cbg1,cbg2,cbg3,mdl%cell)

       ! Find vector to origin of FFTbox
       r = real(bs1-1,kind=DP) * a1 &
            + real(bs2-1,kind=DP) * a2 &
            + real(bs3-1,kind=DP) * a3
#endif

       R_fft(1) = r%X
       R_fft(2) = r%Y
       R_fft(3) = r%Z

       ! ddor: add shifts to matrix elements of position
       !       operator and position operator squared
       ingwf = loc_ingwf + ngwf_basis%first_on_proc(pub_my_proc_id) - 1
       do jngwf=ingwf,ngwf_basis%num

          ! Test if these NGWFs overlap
          if (.not.sparse_element_exists(overlap,jngwf,ingwf)) cycle

          if (nonorthogonal) then
             call sparse_get_element(o_el,overlap,jngwf,ingwf)
          else
             call sparse_get_element(o_el,closetounity,jngwf,ingwf)
          endif
          do xyz=1,3
             call sparse_get_element(r_el(xyz),r_elements(xyz),jngwf,ingwf)
             call sparse_get_element(r2_el(xyz),r2_elements(xyz),jngwf,ingwf)
             !ddor: <a| r | b>  =  <a| (r - R) | b> +  R <a | b>
             r_el(xyz) = (R_fft(xyz) * o_el) + r_el(xyz)
             !ddor: <a| r** 2 |b>  =  <a| (r - R)**2 | b> +  2 <a| R.r |b> -
             !                        <a| R**2 |b>
             !                     =  <a| (r - R)**2 | b> +  2 R.<a| r |b> -
             !                        R**2 <a | b>
             r2_el(xyz) = r2_el(xyz) + &
                  (2.0_DP * R_fft(xyz) * r_el(xyz)) - &
                  (R_fft(xyz) * R_fft(xyz) * o_el )
             call sparse_put_element(r2_el(xyz),r2_elements(xyz),jngwf,ingwf)
             call sparse_put_element(r_el(xyz),r_elements(xyz),jngwf,ingwf)
          enddo
       enddo
    enddo

    call comms_barrier

    do xyz=1,3

       if (nonorthogonal) then

          ! ddor: <a| r |b> S^bc
          call sparse_product(rs(xyz),r_elements(xyz),inv_overlap)
          ! ddor: <a| r |b> S^bc <c| r | d>
          call sparse_product(rsr(xyz),rs(xyz),r_elements(xyz))
          ! ddor: The negative of the covariant NGWF spread matrix
          !       <a| r |b> S^bc <c| r | d> - <a| r^2 |d>
          call sparse_axpy(rsr(xyz),r2_elements(xyz),-1.0_DP)

          ! ddor: The negative of the invariant NGWF spread matrix
          !       ( <a| r |b> S^bc <c| r | d> - <a| r^2 |d> ) S^de
          call sparse_product(thespread(xyz),rsr(xyz),inv_overlap)

       else
          ! ddor: Case of orthogonalised NGWFs

          ! ddor: <a| r |c> <c| r |b>
          call sparse_product(thespread(xyz),r_elements(xyz),r_elements(xyz))
          ! ddor: The negative of the NGWF spread matrix
          !       <a| r |c> <c| r | b> - <a| r^2 |b>
          call sparse_axpy(thespread(xyz),r2_elements(xyz),-1.0_DP)

       endif

       ! ddor: Collect diagonal elements of Spread matrices
       do loc_ingwf=1,ngwf_basis%num_on_proc(pub_my_proc_id)
          ingwf = loc_ingwf + ngwf_basis%first_on_proc(pub_my_proc_id) - 1
          if (nonorthogonal) then
             call sparse_get_element(covariant(ingwf,xyz),&
                  rsr(xyz),ingwf,ingwf)
             covariant(ingwf,xyz) = -1.0_DP * covariant(ingwf,xyz)
          endif
          call sparse_get_element(invariant(ingwf,xyz),&
               thespread(xyz),ingwf,ingwf)
          invariant(ingwf,xyz) = -1.0_DP * invariant(ingwf,xyz)
       enddo

    enddo

    ! ddor: Sum of spreads in each direction
    !       Could also calculate spread anisotropy and output Wannier centres
    do loc_ingwf=1,ngwf_basis%num_on_proc(pub_my_proc_id)
       ingwf = loc_ingwf + ngwf_basis%first_on_proc(pub_my_proc_id) - 1
       if (nonorthogonal) then
          covariant(ingwf,4) = SUM(covariant(ingwf,1:3))
       endif
       invariant(ingwf,4) = SUM(invariant(ingwf,1:3))
    enddo

    call comms_reduce('SUM',invariant)

    ! ddor: If we need to print out both covariant and invariant spreads
    if (nonorthogonal) then

       call comms_reduce('SUM',covariant)

       ! Print table of spreads, adapted from properties_ngwfs_char by mr
       if (pub_on_root) then
          write(stdout,'(/a)') '****************************************&
               &****************************************'
          write(stdout,'(a1,2x,a4,2x,a4,1x,8x,a1,1x,a25,1x,a25,4x,a1)') &
               &'|','Atom','NGWF',':','Covariant Spread (bohr^2)',&
               &'Invariant Spread (bohr^2)','|'
          do atom=1,par%nat
             ingwf = ngwf_basis%first_on_atom(par%distr_atom(atom))
             do ngwf_on_atom=1,ngwf_basis%num_on_atom(par%distr_atom(atom))
                if (ngwf_on_atom==1) then
                   write(stdout,'(a1,2x,a2,2x,i6,1x,i6,2x,a1,2x,f24.10,2x,&
                        &f24.10,4x,a1)') &
                        &'|',mdl%elements(atom)%symbol,atom,ngwf_on_atom,':',&
                        covariant(ingwf,4),invariant(ingwf,4),'|'
                else
                   write(stdout,'(a1,2x,2x,2x,6x,1x,i6,2x,a1,2x,f24.10,2x,&
                        &f24.10,4x,a1)') &
                        &'|',ngwf_on_atom,':',&
                        covariant(ingwf,4),invariant(ingwf,4),'|'
                endif
                ingwf = ingwf + 1           ! original NGWF counting
             enddo
          enddo

          write(stdout,'(a/)') '****************************************&
               &****************************************'
       endif

    else

       ! Print table of spreads, adapted from properties_ngwfs_char by mr
       if (pub_on_root) then
          write(stdout,'(/a)') '****************************************&
               &****************************************'
          write(stdout,'(a1,2x,a4,2x,a4,1x,8x,a1,1x,a25,4x,a1)') &
               &'|','Atom','NGWF',':','Spread (bohr^2)','|'
          do atom=1,par%nat
             ingwf = ngwf_basis%first_on_atom(par%distr_atom(atom))
             do ngwf_on_atom=1,ngwf_basis%num_on_atom(par%distr_atom(atom))
                if (ngwf_on_atom==1) then
                   write(stdout,'(a1,2x,a2,2x,i6,1x,i6,2x,a1,2x,&
                        &f24.10,4x,a1)') &
                        &'|',mdl%elements(atom)%symbol,atom,ngwf_on_atom,':',&
                        invariant(ingwf,4),'|'
                else
                   write(stdout,'(a1,2x,2x,2x,6x,1x,i6,2x,a1,2x,&
                        &f24.10,4x,a1)') &
                        &'|',ngwf_on_atom,':',&
                        invariant(ingwf,4),'|'
                endif
                ingwf = ingwf + 1           ! original NGWF counting
             enddo
          enddo

          write(stdout,'(a/)') '****************************************&
               &****************************************'
       endif

    endif

    ! Deallocate sparse matrices
    if (.not. nonorthogonal) then
       call sparse_destroy(closetounity)
    endif

    do xyz=1,3
       call sparse_destroy(thespread(xyz))
       if (nonorthogonal) then
          call sparse_destroy(rsr(xyz))
          call sparse_destroy(rs(xyz))
       endif
       call sparse_destroy(r2_elements(xyz))
       call sparse_destroy(r_elements(xyz))
    enddo

    ! ndmh: deallocate temporary arrays
    deallocate(covariant,stat=ierr)
    call utils_dealloc_check('properties_spread','covariant',ierr)
    deallocate(invariant,stat=ierr)
    call utils_dealloc_check('properties_spread','invariant',ierr)

    ! Stop timer
    call timer_clock('properties_spread',2)

  end subroutine properties_spread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine properties_plot_delta_density(filename_core, density_fine1, &
        density_fine2, mdl)

    !==========================================================================!
    ! Calculates the difference in density between two input densities.        !
    ! (= density_fine1 - density_fine2).                                       !
    ! Used by the EDA to calculate delta densities between the EDA states in   !
    ! order to visualise the EDA processes (e.g. charge transfer,              !
    ! polarisation...).                                                        !
    !--------------------------------------------------------------------------!
    ! Written 09/02/2014 by Max Phipps                                         !
    ! based on properties_calculate:internal_plot_dens_and_pot                 !
    !==========================================================================!

    use constants, only: ANGSTROM
    use model_type, only: MODEL
    use rundat, only: pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check
    use visual, only: visual_scalarfield

    implicit none

    ! Arguments
    ! the file name core string:
    character(len=*), intent(in) :: filename_core
    ! supermolecule system variables:
    type(MODEL), intent(in) :: mdl
    ! densities used to calculate delta density:
    real(kind=DP), intent(in) :: density_fine1(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in) :: density_fine2(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins)

    ! Local Variables
    integer :: is

    real(kind=DP), allocatable :: d_density_fine(:,:,:,:)! delta density
    integer :: ierr          ! memory allocation error flag
!   real(kind=DP), allocatable :: nhat_den_grad(:,:,:,:,:)! compensation den
!   type(SPAM3), allocatable :: rhoij(:)


    ! mjsp: allocate delta density
    allocate(d_density_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
    call utils_alloc_check('properties_plot_delta_density', &
         'd_density_fine', ierr)

    ! mjsp:########## DENSITY OUTPUT ##########

    ! mjsp: calculate the delta density from the passed densities:
    do is=1,pub_num_spins
       d_density_fine(:,:,:,is) = density_fine1(:,:,:,is) - &
             density_fine2(:,:,:,is)
    end do

    ! mjsp: output delta density in plot format file
    ! vm: output density in Angstrom rather than in Bohr
    ! jd: ... but leave the unit conversion to visual_scalarfield
    call visual_scalarfield( &
         d_density_fine(:,:,:,1), mdl%fine_grid, mdl%cell, &
         'Electronic delta density (in e/ang^3) for:', filename_core//'density', mdl%elements, &
         ANGSTROM**3)
    if (pub_num_spins == 2) call visual_scalarfield( &
         d_density_fine(:,:,:,2), mdl%fine_grid, mdl%cell, &
         'Electronic delta spin density (in e/ang^3) for:', filename_core//'spindensity', &
         mdl%elements, ANGSTROM**3)

    ! mjsp:####### END DENSITY OUTPUT ##########

    ! mjsp: deallocate delta density
    deallocate(d_density_fine, stat=ierr)
    call utils_dealloc_check('properties_plot_delta_density', &
         'd_density_fine', ierr)

  end subroutine properties_plot_delta_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine properties_polarisation_berry(rep,ngwf_basis,proj_basis,mdl, &
       proj_set,denskern,axis,dipole,dipole_mat)

    !==========================================================================!
    ! Calculates dipole moment of a molecule by computing matrix elements      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! rep (input)        : NGWF representation for valence                     !
    ! ngwf_basis (input) : Function basis type describing the NGWFs            !
    ! denskern (input)   : Density kernel                                      !
    !--------------------------------------------------------------------------!
    ! Built upon the old properties_polarisation routine, but contains         !
    ! the necessary corrections for periodic systems at the Gamma point        !
    ! If one needs to extend this to multiple k-points, one additional term    !
    ! needs to be implemented                                                  !
    ! -------------------------------------------------------------------------!
    ! Written by G.C. Constantinescu in Dec 2016                               !
    ! Modified for embedding and to remove pub_par by Joseph Prentice, August  !
    ! 2018                                                                     !
    !==========================================================================!

    use basis, only : basis_location_func_wrt_cell, &
         basis_point_wrt_box, basis_start_of_box_wrt_cell
    use comms, only: pub_on_root, comms_reduce
    use constants, only: DP, stdout, ANGSTROM, UP, DN, max_spins
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_brappd_ketppd
    use geometry, only: POINT, operator(*), operator(+)
    use optics, only: optics_pos_mat_els
    use md_ions, only: md_ionic_polarisation, md_elec_polarisation, &
         md_total_polarisation
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use paw, only: paw_position_operator, paw_projector_overlap
    use projectors, only: PROJECTOR_SET, projectors_func_ovlp_box
    use rundat, only: pub_print_qc, pub_task, &
         pub_num_spins, pub_spin_fac, pub_num_kpoints, PUB_1K, pub_paw
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create, &
         sparse_embed_trace, sparse_embed_destroy, sparse_embed_copy, &
         sparse_embed_product, sparse_embed_transpose, sparse_embed_axpy, &
         sparse_embed_scale, sparse_embed_transpose_structure
    use timer, only: timer_clock
    use utils, only: utils_qc_print, utils_assert, utils_banner

    ! Arguments
    type(NGWF_REP), intent(in)   :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(PROJECTOR_SET), intent(inout) :: proj_set(:)
    type(MODEL), intent(in)      :: mdl
    type(SPAM3_EMBED_ARRAY), intent(in):: denskern
    ! ddor: When axis is present, only the electronic contribution to the dipole
    !       moment in direction axis is calculated, and this is output to dipole.
    integer, optional, intent(in) :: axis
    real(kind=DP), optional, intent(out) :: dipole

    ! lpl: Accepts external dipole_mat buffer, and returns dipole mat
    type(SPAM3_EMBED), optional, intent(inout) :: dipole_mat(3)

    ! Local Variables
    type(SPAM3_EMBED) :: r_elements(3), factored_overlap(3), factored_sp_overlap(3)
    type(SPAM3_EMBED) :: r_elements_corr(3), O_mat
    type(SPAM3_EMBED) :: temp_mat_1, temp_mat_2, temp_mat_3
    real(kind=DP) :: d(3),q(3),t(3), trace
#ifdef FFTBOX_OLD_POS
    integer       :: cs1,cs2,cs3
#endif
    integer       :: xyz,axmin,axmax,is
    integer       :: isub, jsub
    real(kind=DP) :: coc(3), fscale

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering properties_polarisation_berry'

    ! Start timer
    call timer_clock('properties_polarisation_berry',1)

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine properties_polarisation_berry not ready yet for more&
         & than one k-point.')

    ! ddor: Allocate only one matrix if axis is specified
    if (present(axis)) then
       axmin = axis
       axmax = axis
    else
       axmin = 1
       axmax = 3
    end if

    ! ndmh: create position elements matrix with NGWF-overlap sparsity
    do xyz=axmin,axmax
       call sparse_embed_create(r_elements(xyz),rep%ngwf_overlap)
       call sparse_embed_create(r_elements_corr(xyz),rep%ngwf_overlap)
       call sparse_embed_create(factored_overlap(xyz), rep%overlap)
       if (pub_paw) then
          call sparse_embed_create(factored_sp_overlap(xyz), rep%sp_overlap)
       end if
    end do

    if (pub_paw) then
       O_mat%structure='E'
       call sparse_embed_create(O_mat)
       ! rc2013; temporary fix
       ! jcap: not ready for embedding yet
       call paw_projector_overlap(O_mat%p, mdl%paw_sp)
!       call paw_projector_overlap(O_mat, mdl%paw_sp, mdl%regions(1)%par)
    end if

    ! Calculate matrix elements of position operator in NGWF representation
    call optics_pos_mat_els(r_elements,rep,ngwf_basis,rep,ngwf_basis,&
         rep%overlap,rep%ngwf_overlap,proj_basis, mdl)

    ! calculate center of ionic charge
    fscale = 1.0_DP/sum(mdl%elements(:)%ion_charge)
    coc(1) =  sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%X)*fscale
    coc(2) =  sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%Y)*fscale
    coc(3) =  sum(mdl%elements(:)%ion_charge * mdl%elements(:)%centre%Z)*fscale

    ! PERIODIC SYSTEM CORRECTIONS ============================================
    do xyz = axmin, axmax
       call sparse_embed_copy(r_elements_corr(xyz), rep%overlap)
       call sparse_embed_scale(r_elements_corr(xyz), -1.0_DP * coc(xyz))

       ! jcap: loop over regions
       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub
             call function_ops_brappd_ketppd(factored_overlap(xyz)%m(isub,jsub), &
                  rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                  rep%ngwfs_on_grid(jsub), ngwf_basis(jsub), &
                  mdl%cell, factor_type = xyz)
          end do
       end do

       if (pub_paw) then
          ! jcap: loop over regions
          ! jcap: THIS MIGHT NOT BE CORRECT - temporary fix for now
          do isub=1,mdl%nsub
             do jsub=1,mdl%nsub
                call projectors_func_ovlp_box(factored_sp_overlap(xyz)%m(isub,jsub), &
                     rep%ngwfs_on_grid(isub), ngwf_basis(isub),&
                     proj_basis(jsub),proj_set(jsub),mdl%fftbox, &
                     mdl%cell,factor_type = xyz)
             end do
          end do
          call sparse_embed_create(temp_mat_1, factored_sp_overlap(xyz), O_mat)
          call sparse_embed_product(temp_mat_1, factored_sp_overlap(xyz), O_mat)
          call sparse_embed_transpose_structure(temp_mat_2%structure, rep%sp_overlap)
          call sparse_embed_create(temp_mat_2)
          call sparse_embed_transpose(temp_mat_2, rep%sp_overlap)
          call sparse_embed_create(temp_mat_3, temp_mat_1, temp_mat_2)
          call sparse_embed_product(temp_mat_3, temp_mat_1, temp_mat_2)
          call sparse_embed_axpy(factored_overlap(xyz), temp_mat_3, 1.0_DP)
          call sparse_embed_destroy(temp_mat_1)
          call sparse_embed_destroy(temp_mat_2)
          call sparse_embed_destroy(temp_mat_3)
       end if

       call sparse_embed_axpy(r_elements_corr(xyz), factored_overlap(xyz), -1.0_DP)

       call sparse_embed_axpy(r_elements(xyz), r_elements_corr(xyz), 1.0_DP)
    end do

    ! END OF PERIODIC SYSTEM CORRECTIONS =====================================

    ! Calculate dipole moment from electrons by Tr(KR)
    d = 0.0_DP

    ! ddor: Calculate for one direction only if axis is specified
    axiscase: if (present(axis)) then

       do is=1,pub_num_spins
          call sparse_embed_trace(trace,denskern%m(is,PUB_1K),r_elements(axis))
          d(axis) = d(axis) - trace
       end do

       ! multiplication factor depending on number of spins
       if (pub_num_spins /= 2) d(:) = pub_spin_fac * d(:)

       ! calculate dipole moment from ions
       q(1) = sum(mdl%elements(:)%ion_charge * (mdl%elements(:)%centre%X - &
            coc(1)))
       q(2) = sum(mdl%elements(:)%ion_charge * (mdl%elements(:)%centre%Y - &
            coc(2)))
       q(3) = sum(mdl%elements(:)%ion_charge * (mdl%elements(:)%centre%Z - &
            coc(3)))

       ! calculate total dipole moment
       t = q + d

       ! print pretty table
       if (pub_on_root) then
          write(stdout,'(/a)') utils_banner('=', 'Dipole Moment Calculation')
          write(stdout,'(a34,4x,a)') 'Type:','Single Molecule'
          write(stdout,'(a34,4x,a2,1x,i1,1x,a3,1x,f24.18)') 'Electronic dipole moment (e.bohr):',&
               &'d(',axis,') =',d(axis)
          write(stdout,'(a34,4x,a2,1x,i1,1x,a3,1x,f24.18)') 'Ionic dipole moment (e.bohr):',&
               &'q(',axis,') =',q(axis)
          write(stdout,'(a34,4x,a2,1x,i1,1x,a3,1x,f24.18)') 'Total dipole moment (e.bohr):',&
               &'t(',axis,') =',t(axis)
          write(stdout,'(a)') repeat('=',80)
       end if

       dipole = d(axis)

    else

       do xyz=1,3

          do is=1,pub_num_spins
             call sparse_embed_trace(trace,denskern%m(is,PUB_1K),r_elements(xyz))
             d(xyz) = d(xyz) - trace
          end do

       end do

       ! Multiplication factor depending on number of spins
       if (pub_num_spins /= 2) d(:) = pub_spin_fac * d(:)

       ! Calculate dipole moment from ions
       q(1) = sum(mdl%elements(:)%ion_charge * (mdl%elements(:)%centre%X - &
              coc(1)))
       q(2) = sum(mdl%elements(:)%ion_charge * (mdl%elements(:)%centre%Y - &
              coc(2)))
       q(3) = sum(mdl%elements(:)%ion_charge * (mdl%elements(:)%centre%Z - &
              coc(3)))

       ! Calculate total dipole moment
       t = q + d

       ! Print pretty table
       if(pub_on_root) then
          write(stdout,'(/a)') utils_banner('=', 'Dipole Moment Calculation')
          write(stdout,'(a34,4x,a)') 'Type:','Single Molecule'
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Electronic dipole moment &
               &(e.bohr):','dx =',d(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',d(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',d(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(d(:)**2))
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Ionic dipole moment (e.bohr):', &
               'dx =',q(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',q(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',q(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(q(:)**2))
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Total dipole moment (e.bohr):', &
               'dx =',t(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',t(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',t(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(t(:)**2))
          write(stdout,*)
          if (pub_print_qc) then
             call utils_qc_print('elec_polarisation_d1',d(1))
             call utils_qc_print('elec_polarisation_d2',d(2))
             call utils_qc_print('elec_polarisation_d3',d(3))
             call utils_qc_print('ionic_polarisation_d1',q(1))
             call utils_qc_print('ionic_polarisation_d2',q(2))
             call utils_qc_print('ionic_polarisation_d3',q(3))
             call utils_qc_print('total_polarisation_d1',t(1))
             call utils_qc_print('total_polarisation_d2',t(2))
             call utils_qc_print('total_polarisation_d3',t(3))
          end if
          write(stdout,'(a)') repeat('=',80)

          ! smmd: print data in md file
          if (pub_task == 'MOLECULARDYNAMICS') then
             md_elec_polarisation(:) = d(:)
             md_ionic_polarisation(:) = q(:)
             md_total_polarisation(:) = t(:)
          end if

       endif

    end if axiscase

    ! ndmh: copy to output matrix if required
    if (present(dipole_mat)) then
       do xyz=axmin,axmax
          call sparse_embed_copy(dipole_mat(xyz),r_elements(xyz))
       end do
    end if

    ! Deallocate sparse matrices
    do xyz=axmax,axmin,-1
       call sparse_embed_destroy(r_elements(xyz))
       call sparse_embed_destroy(r_elements_corr(xyz))
       call sparse_embed_destroy(factored_overlap(xyz))
       if (pub_paw) then
          call sparse_embed_destroy(factored_sp_overlap(xyz))
          call sparse_embed_destroy(O_mat)
       end if
    end do

    ! Stop timer
    call timer_clock('properties_polarisation_berry',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving properties_polarisation_berry'

  end subroutine properties_polarisation_berry


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine properties_electrostatic_potential(mdl, denskern, rep, ngwf_basis,&
       work_reg_density_on_fine, work_reg_nhat_den_grad)
    !==========================================================================!
    ! Calculates and outputs the 'electrostatic potential', that is:           !
    ! (SI: smeared ions, CC: cutoff Coulomb, MT: Martyna-Tuckerman)            !
    ! - if SI are in use: the potential seen by the multigrid solver -- due to !
    !      electrons, smeared ions, and -- if solvent present -- due to the    !
    !      reaction field of the solvent and electrolyte.                      !
    ! - if SI are not in use: the potential due to electrons and pseudoions.   !
    ! In either case, PAW contributions are added, if present.                 !
    ! When SI are not in use, it works fine with CC and MT too.                !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! mdl (input):       : the model                                           !
    ! denskern (input)   : the density kernel
    ! rep (input)        : NGWF representation for valence                     !
    ! ngwf_basis (input) : Function basis type describing the NGWFs            !
    ! denskern (input)   : Density kernel                                      !
    ! work_reg_density_on_fine: } work arrays allocated in the caller, holding !
    ! work_reg_nhat_den_grad:   } per-region density and PAW density on fine.  !
    ! -------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !  - work_reg_* arrays are assumed to have been allocated to the correct   !
    !    size in the caller.                                                   !
    !  - trying qc_print with very small grids will cause an array overrun --  !
    !    check the hardcoded gridpoint indices below in                        !
    !    internal_qc_check_of_grid_data.                                       !
    ! -------------------------------------------------------------------------!
    ! Originally written by C.-K. Skylaris and J. Aarons. Extracted to here    !
    ! generalised, bugfixed, and refactored by J. Dziedzic in June 2021.       !
    !==========================================================================!

    use augmentation, only: augmentation_density_on_grid
    use constants, only: UP, DN, HARTREE_IN_EVS
    use cutoff_coulomb, only: cutoff_coulomb_hartree
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use hartree, only: hartree_via_multigrid, hartree_on_grid
    use kernel, only: DKERN
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: PUB_1K, pub_multigrid_bc_is_periodic, pub_num_spins, &
         pub_aug, pub_is_smeared_ion_rep, pub_coulomb_cutoff
    use sparse, only: SPAM3
    use sparse_embed, only: sparse_embed_extract_from_array, &
         sparse_embed_destroy_extracted_array
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check
    use visual, only: visual_scalarfield

    implicit none

    ! Arguments
    type(MODEL), intent(in)      :: mdl
    type(DKERN), intent(in)      :: denskern
    type(NGWF_REP), intent(in)   :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(mdl%nsub)
    real(kind=DP), allocatable, intent(inout) :: work_reg_density_on_fine(:,:,:,:)
    real(kind=DP), allocatable, intent(inout) :: work_reg_nhat_den_grad(:,:,:,:,:)

    ! Local variables
    real(kind=DP), allocatable :: pot_on_fine(:,:,:,:)
    real(kind=DP), allocatable :: density_on_fine(:,:,:,:)
    type(SPAM3) :: kern_array(pub_num_spins)
    logical :: pbc
    character(len=128) :: pot_detail_string
    integer :: isub, jsub
    integer :: ierr
    character(len=*), parameter :: myself = "properties_electrostatic_potential"

    ! -------------------------------------------------------------------------

    if(all(pub_multigrid_bc_is_periodic)) then
       pbc = .true.
    else if(.not. any(pub_multigrid_bc_is_periodic)) then
       pbc = .false.
    else
       call utils_abort(myself//": Mixed boundary conditions are not supported")
    end if

    allocate(pot_on_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
    call utils_alloc_check(myself,'pot_on_fine',ierr)
    allocate(density_on_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
    call utils_alloc_check(myself,'pot_on_fine',ierr)

    pot_on_fine(:,:,:,:) = 0.0_DP
    density_on_fine(:,:,:,:) = 0.0_DP

    ! jd: Accumulate density on the fine grid from all regions.
    do isub = 1, mdl%nsub
       do jsub = 1, mdl%nsub
          work_reg_density_on_fine(:,:,:,:) = 0.0_DP
          ! jcap: create temporary density kernel
          call sparse_embed_extract_from_array(kern_array, &
               denskern%kern%m(:,PUB_1K),isub,jsub)

          call density_on_grid(work_reg_density_on_fine, mdl%fine_grid, &
               mdl%dbl_grid, mdl%cell, mdl%fftbox, kern_array, &
               rep%ngwf_overlap%m(isub,jsub), rep%ngwfs_on_grid(isub), &
               ngwf_basis(isub), rep%ngwfs_on_grid(jsub), ngwf_basis(jsub))

          density_on_fine = density_on_fine + work_reg_density_on_fine

          call sparse_embed_destroy_extracted_array(kern_array)
       end do
    end do

    ! jd: If PAW, include additional density nhat_den_grad
    if(pub_aug) then

       do isub = 1, mdl%nsub
          do jsub = 1, mdl%nsub

             work_reg_nhat_den_grad(:,:,:,:,:) = 0.0_DP

             ! jcap: create temporary density kernel
             call sparse_embed_extract_from_array(kern_array, &
                  denskern%kern%m(:,PUB_1K),isub,jsub)

             call augmentation_density_on_grid(work_reg_nhat_den_grad, &
                  mdl%fine_grid, mdl%cell, mdl%regions(isub)%pseudo_sp, &
                  mdl%regions(jsub)%paw_sp, mdl%aug_box, kern_array, &
                  rep%sp_overlap%m(isub,jsub))

             density_on_fine = density_on_fine + work_reg_nhat_den_grad(:,:,:,:,0)

             call sparse_embed_destroy_extracted_array(kern_array)
          end do
       end do

    end if ! PAW

    ! jd: Case with smeared ions
    if(pub_is_smeared_ion_rep) then
       call hartree_via_multigrid(pot_on_fine, density_on_fine, &
            mdl%fine_grid, mdl%cell, elements = mdl%elements, &
            no_dump = .true.)
       pot_detail_string = "electrons + smeared ions"
    else
       ! jd: Case without smeared ions. Cannot call hartree_via_multigrid()
       !     because pub_multigrid_hartree is F. We thus call hartree_on_grid()
       !     for the electronic part and include the local PS term manually.
       if(pub_coulomb_cutoff) then
          call cutoff_coulomb_hartree(pot_on_fine, density_on_fine, mdl)
          pot_detail_string = "electrons + pseudoions (CC)"
       else
          call hartree_on_grid(pot_on_fine, density_on_fine, &
               mdl%fine_grid, mdl%cell)
          pot_detail_string = "electrons + pseudoions"
       end if

       pot_on_fine(:,:,:,1) = pot_on_fine(:,:,:,1) + mdl%localpseudo_fine(:,:,:)

    end if

    call internal_qc_check_of_grid_data()

    if(pub_aug) pot_detail_string = trim(pot_detail_string)//" + PAW"

    call visual_scalarfield(pot_on_fine(:,:,:,1), mdl%fine_grid, mdl%cell, &
         'Electrostatic pot. (eV) ('//trim(pot_detail_string)//') for:', &
         '_electrostatic_potential', mdl%elements, -HARTREE_IN_EVS)

    ! jd: Clean up
    deallocate(density_on_fine, stat=ierr)
    call utils_dealloc_check(myself,'density_on_fine',ierr)
    deallocate(pot_on_fine, stat=ierr)
    call utils_dealloc_check(myself,'pot_on_fine',ierr)

    return

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine internal_qc_check_of_grid_data

      use comms, only: pub_on_root, pub_my_proc_id, comms_bcast
      use rundat, only: pub_print_qc, pub_coulomb_cutoff
      use utils, only: utils_qc_print

      implicit none

      ! jd: Local variables
      integer :: third_index
      integer :: ii
      integer :: qc_owner_proc, qc_first_on_grid
      integer, parameter :: qc_npoints = 10
      character(len=4) :: qc_point_char
      real(kind=DP) :: qc_point(qc_npoints)
      ! jd: Coordinates of the gridpoints for printing in QC tests.
      !     1st set for cutoff Coulomb, 2nd set otherwise.
      !     NB the reshape is needed ("The result of an array constructor
      !     expression is always a 1D array") and gfortran will catch this.
      integer, parameter :: qc_coords(3,10,2) = reshape([ &
           [ [39, 39, 10], [39, 39, 20], [39, 39, 30], [39, 39, 39], &
           [39, 39, 50], [20, 39, 39], [50, 39, 39], &
           [39, 20, 39], [39, 50, 39], [50, 50, 50] ], &
           [ [62, 53, 41], [38, 59, 49], [38, 47, 41], [50, 40, 39], &
           [53, 40, 39], [42, 51, 63], [38, 46, 59], &
           [62, 58, 53], [49, 41, 62], [45, 39, 40] ] ], (/3,10,2/) )

      ! -----------------------------------------------------------------------

      if(pub_print_qc) then

         if(pub_on_root .and. pub_print_qc .and. pub_aug) then
            call utils_qc_print('props_ES_ld1',mdl%fine_grid%ld1)
            call utils_qc_print('props_ES_ld2',mdl%fine_grid%ld2)
            call utils_qc_print('props_ES_n3',mdl%fine_grid%n3)
            call utils_qc_print('props_ES_num_spins',pub_num_spins)
         end if

         do ii=1,qc_npoints

            third_index = merge(1,2,pub_coulomb_cutoff)
            qc_owner_proc = mdl%fine_grid%proc_slab12(qc_coords(3,ii,third_index))
            qc_first_on_grid = mdl%fine_grid%first_slab12(qc_owner_proc)

            if(qc_owner_proc == pub_my_proc_id) then
               qc_point(ii) = pot_on_fine(qc_coords(1,ii,third_index), &
                    qc_coords(2,ii,third_index), qc_coords(3,ii,third_index) &
                    - qc_first_on_grid+1,1)
            end if

            call comms_bcast(qc_owner_proc,qc_point(ii))

            write(qc_point_char,'(i0)') ii

            if(pub_on_root) then
               call utils_qc_print('props_ES_n'//trim(adjustl(qc_point_char)),&
                    qc_point(ii))
            end if
         end do

      end if

    end subroutine internal_qc_check_of_grid_data

  end subroutine properties_electrostatic_potential

end module properties
