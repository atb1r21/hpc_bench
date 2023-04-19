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
!=================================================================!

module conduction

  use constants, only: DP

  implicit none

  private

  public :: conduction_ngwf_optimise

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine conduction_ngwf_optimise(total_energy, val_denskern, val_rep, &
       val_ham, val_ngwf_basis, proj_basis, nl_projectors, hub_proj_basis, &
       hub, edft, cond_ngwf_basis, joint_ngwf_basis, core_basis, core_wvfns, &
       mdl, hfxstate, lhxc_fine, cond_properties_only, kpt) ! agrecokpt

    !==========================================================================!
    ! This subroutine optimises a set of conduction NGWFs by forming a         !
    ! Hamiltonian which projects out valence contributions to the energy.      !
    ! It first initialises the appropriate valence variables and also calls    !
    ! properties_calculate at the end in which the various Hamiltonian forms   !
    ! are diagonalised to give the individual eigenvalues.                     !
    !------------------------------------------------------------------------  !
    ! total_energy    (out) : The final projected conduction energy            !
    ! val_denskern  (inout) : The ground state density kernel                  !
    ! val_rep          (in) : Valence NGWF representation                      !
    ! val_ham          (in) : Valence NGWF Hamiltonian                         !
    ! val_ngwf_basis   (in) : Array of function basis type for valence NGWFs   !
    ! proj_basis       (in) : Array of function basis type for nonlocal        !
    !                              pseudopotential projectors                  !
    ! nl_projectors    (in) : Array of type describing nonlocal projectors     !
    ! hub_proj_basis   (in) : Array of function basis type for Hubbard projs   !
    ! hub              (in) : Hubbard Model type for Hubbard projectors        !
    ! core_basis       (in) : Array of function basis type for core wavefuncs  !
    ! core_wvfns       (in) : Projector Set array for the core wavefunctions   !
    ! cond_ngwf_basis  (in) : Function basis type array for conduction NGWF    !
    ! joint_ngwf_basis (in) : Function basis type array for the joint NGWFs    !
    ! lhxc_fine     (inout) : Local-Hartree-Exhange-Correlation potential      !
    ! @docme: cond_properties_only, kpt
    !--------------------------------------------------------------------------!
    ! Written by Laura Ratcliff in October 2010.                               !
    ! Modified calls to properties to enable the calculation of matrix         !
    ! elements for optical absorption spectra by Laura Ratcliff March 2011.    !
    ! Minor modification for hybrid conduction by Jacek Dziedzic, July 2018.   !
    ! Modified for embedding by Joseph Prentice, June 2018                     !
    !==========================================================================!

    use datatypes, only: data_functions_alloc, data_functions_dealloc
    use bibliography, only: bibliography_cite
    use comms, only: comms_barrier, pub_on_root
    use conduction_properties, only: conduction_properties_calculate
    use constants, only: stdout, paw_en_size, NORMAL, REP_SWEX_HFX_OTHER
    use electronic_init, only: electronic_init_denskern
    use ensemble_dft_type, only: EDFT_MODEL
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_dens_dep_matrices, &
         hamiltonian_build_matrix, hamiltonian_dens_dep_nonsc
    use hf_exchange, only: HFX_STATE, hf_exchange_init, &
         hf_exchange_dkn_indep_stage
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: DKERN, kernel_create, kernel_destroy
    use model_type, only: MODEL
    use ngwf_cg, only: ngwf_cg_optimise
    use ngwf_representation, only: NGWF_REP, NGWF_HAM, ngwf_rep_create, &
         ngwf_rep_destroy, ngwf_ham_create, ngwf_ham_destroy
    use ngwfs, only: ngwfs_initialise
    use parallel_strategy, only: PARAL_INFO
    use projectors, only: PROJECTOR_SET
    use rundat, only: cond_num_extra_states, cond_num_extra_its, &
         pub_maxit_ngwf_cg, pub_output_detail, pub_aug, &
         cond_read_denskern, cond_init_shift, pub_debug_on_root, &
         pub_num_spins, pub_num_kpoints, PUB_1K, pub_kpoint_method, &
         pub_xc_ke_density_required, pub_use_hfx, cond_num_states, &
         pub_emft, pub_use_activehfx, pub_active_region
    use services, only: services_flush
    use sparse_embed, only: sparse_embed_copy
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(FUNC_BASIS), intent(in) :: core_basis(:)
    type(FUNC_BASIS), intent(inout) :: hub_proj_basis(:)
    type(FUNC_BASIS), intent(inout) :: cond_ngwf_basis(:)
    type(FUNC_BASIS), intent(inout) :: joint_ngwf_basis(:)
    type(MODEL), intent(inout), target :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    type(DKERN), intent(inout)   :: val_denskern
    type(NGWF_REP), intent(in)   :: val_rep
    type(NGWF_HAM), intent(inout):: val_ham
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(:)
    type(PROJECTOR_SET), intent(inout) :: core_wvfns(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(out) :: total_energy
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, val_rep%nsub)
    logical, intent(in) :: cond_properties_only
    type(EDFT_MODEL), intent(inout) :: edft ! ars: ensemble-DFT container
    ! agrecokpt: optional argument to specify non-Gamma k-point
    real(kind=DP), optional, intent(in) :: kpt(3)

    ! Local Variables
    type(NGWF_REP) :: cond_rep
    type(NGWF_HAM) :: cond_ham
    type(DKERN) :: cond_denskern
    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    real(kind=DP) :: hubbard_energy
    real(kind=DP) :: lhxc_energy
    integer :: is, isub
    integer :: tmp_maxit_ngwf_cg
    logical :: converged
    logical :: out_of_runtime
    logical :: shift_changed
    ! ars: ngwf_nonsc_forces
    real(kind=DP) :: ngwf_nonsc_forces(1:3,mdl%nat)
    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecokpt
    real(kind=DP) :: loc_kpt(3)
    type(PARAL_INFO), pointer :: par

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &conduction_ngwf_optimise'

    ! JCW: Abort if KE density required (meta-GGA + valence NGWFs not
    ! JCW: implemented).
    call utils_assert(.not.pub_xc_ke_density_required, "Error in &
         &confuction_ngwf_optimise: pub_xc_ke_density_required is true, but &
         &combination of KE-density-dependent XC functionals with valence &
         &NGWFs has not been implemented/tested.")

    call bibliography_cite('COND')

    ! agrecocmplx
    loc_cmplx = val_rep%ngwfs_on_grid(1)%iscmplx

    ! jme: SPINS_DANGER: this subroutine mixes spins
    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine internal_cond_occupancies (conduction) not ready yet for more&
         & than one k-point.')

    ! agrecokpt: currently only KP method is being implemented
    call utils_assert(pub_kpoint_method == 'KP', &
         'Subroutine conduction_ngwf_optimise currently supports&
         & only KP method for BZ sampling')

    ! agrecokpt: default is Gamma
    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt(:) = 0.0_DP
    end if

    ! ndmh: Allocate SPAM3 structures for overlap matrix, kinetic energy,
    ! ndmh: density kernel, inverse overlap
    ! agrecocmplx
    call ngwf_rep_create(cond_rep,'c',mdl,is_cmplx=loc_cmplx)

    ! ndmh: ngwf_basis%n_ppds is the number of ppds for the
    ! ndmh: ngwfs belonging to pub_my_proc_id
    ! agrecocmplx
    ! jcap: loop over regions
    do isub=1,mdl%nsub
       call data_functions_alloc(cond_rep%ngwfs_on_grid(isub), &
            cond_ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
    end do
    call comms_barrier

    ! agrecocmplx
    call kernel_create(cond_denskern, 'K'//cond_rep%postfix, &
            is_cmplx=loc_cmplx)

    if ((pub_on_root).and.(pub_output_detail>=NORMAL)) &
         write(stdout,'(a)') ' done'

    ! lr408: Initialise conduction NGWFs to fireballs or read from file
    call ngwfs_initialise(cond_rep, cond_ngwf_basis, mdl, hfxstate, &
         do_not_reexpand = .true.) ! [*]
    call comms_barrier
    call services_flush

    ! ndmh: calculate density dependent energies and matrices for valence NGWFS
    call hamiltonian_dens_dep_matrices(val_ham, lhxc_fine, total_energy, &
         lhxc_energy, hubbard_energy, paw_sphere_energies, val_rep, &
         val_ngwf_basis, hub_proj_basis, hub, val_denskern%kern, &
         mdl, hfxstate, ham_update=.true., lhxc_fixed=.false.)

    call hamiltonian_build_matrix(val_ham, val_rep)

    ! lr408: At some point might be nice to print out the final valence energy
    ! lr408: (components) here, but it's not essential to the calculation so
    ! lr408: for now we'll leave it

    ! ndmh: create conduction hamiltonian
    call ngwf_ham_create(cond_ham,cond_rep)

    ! jd: Initialise HFx for conduction, if needed. Expand necessary NGWF pairs.
    do isub=1,cond_rep%nsub
       ! rc2013: only set up HFx in the active region if needed
       if (pub_use_hfx.or.(pub_use_activehfx.and.(isub==pub_active_region))) then
          call hf_exchange_init(hfxstate, cond_rep, &
               cond_ham%hfexchange(1)%m(isub,isub), mdl%cell, &
               size(mdl%regions(isub)%elements), init_rep_only = .true.) ! [**]
          ! jd: The first time around conduction NGWFs have to be expanded manually
          !     because by the time they are initialised at [*] cond-HFx has not
          !     been initialised yet (this happens at [**]).
          !     We only want the the mixed (cond-val) expansion here, so HFX_OTHER
          !     for the cond rep.
          call hf_exchange_dkn_indep_stage(hfxstate, mdl, isub, &
               mdl%regions(isub)%par, cond_rep, cond_ngwf_basis(isub), &
               rep2 = val_rep, ngwf_basis2 = val_ngwf_basis(isub), &
               basis_selector = (/1,1,2,2,REP_SWEX_HFX_OTHER/))
               !                  ^ qoh's alpha-beta-gamma-delta conv'n
               ! 1:cond, 2: val
               ! swex selector: HFX_OTHER (of cond_rep)
       end if
    end do

    ! ndmh: skip pre-optimisation if only doing a properties calculation
    if (.not. cond_properties_only) then

       if (pub_on_root) then
          write(stdout,'(a)') '+-----------------------------------------------&
               &-------------------------------+'
          write(stdout,'(a)') '|                   Starting conduction NGWF &
               &optimisation                      |'
          write(stdout,'(a)') '+-----------------------------------------------&
               &-------------------------------+'
       end if

       ! lr408: Do a few NGWF iterations optimising for a higher
       ! lr408: number of conduction states than required as a form of
       ! lr408: preconditioning
       ! ndmh:  or, if specifying a window, do the extra iterations first then
       ! ndmh:  re-initialise conduction kernel
       if ((cond_num_extra_its > 0) .and. &
            ((cond_num_extra_states > 0).or.(cond_num_states==0))) then
          if (pub_on_root) then
             write(stdout,'(a)') ''
             write(stdout,'(a)') 'Starting conduction NGWF pre-optimisation &
                  &process'
             write(stdout,'(a)') ''
          end if

          tmp_maxit_ngwf_cg = pub_maxit_ngwf_cg
          pub_maxit_ngwf_cg = cond_num_extra_its

          ! lr408: Set up 'occupation' numbers for conduction states
          call internal_cond_occupancies

          ! lr408: Initialise conduction kernel
          call electronic_init_denskern(cond_rep, cond_ham, cond_denskern, &
               lhxc_fine, cond_ngwf_basis, proj_basis, nl_projectors, &
               hub_proj_basis, hub, mdl, hfxstate, val_rep, val_ngwf_basis, &
               val_denskern%kern, val_ham, kpt=loc_kpt)

          call comms_barrier

          ! ndmh: copy in PAW nonlocal energy screening terms from valence
          ! ndmh: hamiltonian
          if (pub_aug) then
             do is=1,pub_num_spins
                call sparse_embed_copy(cond_ham%dijhat(is),val_ham%dijhat(is))
             end do
          end if

          ! lr408: Now we can optimise the conduction NGWFs
          ! agrecokpt: at specified k-point?
          call ngwf_cg_optimise(total_energy, converged, out_of_runtime, &
               cond_ham, cond_denskern, edft, cond_rep, ngwf_nonsc_forces, &
               lhxc_fine,cond_ngwf_basis, proj_basis, nl_projectors, &
               hub_proj_basis, hub, mdl, hfxstate, &
               val_rep, val_ngwf_basis, val_denskern%kern, val_ham, &
               kpt=loc_kpt)

          pub_maxit_ngwf_cg = tmp_maxit_ngwf_cg
          cond_num_extra_states = 0
          cond_read_denskern = .false.

          if (pub_on_root) then
             write(stdout,'(a)') ''
             write(stdout,'(a)') 'Conduction NGWF pre-optimisation process &
                  &complete.'
             write(stdout,'(a)') 'Proceeding with the main conduction NGWF &
                  &optimisation process.'
             write(stdout,'(a)') ''
          end if

          call comms_barrier

       end if

    end if ! not properties only

    ! lr408: Set up 'occupation' numbers for conduction states
    call internal_cond_occupancies

    ! lr408: Initialise conduction kernel
    ! ndmh: reset cond shift, as we are now re-using same cond_ham as before
    cond_ham%cond_shift = cond_init_shift
    ! agrecokpt: at specified k-point?
    call electronic_init_denskern(cond_rep, cond_ham, cond_denskern, &
         lhxc_fine, cond_ngwf_basis, &
         proj_basis, nl_projectors, hub_proj_basis, hub, mdl, hfxstate, &
         val_rep, val_ngwf_basis, val_denskern%kern, val_ham, &
         kpt=loc_kpt)

    call comms_barrier

    ! ndmh: copy in PAW nonlocal energy screening terms from valence
    ! ndmh: hamiltonian
    if (pub_aug) then
       do is=1,pub_num_spins
          call sparse_embed_copy(cond_ham%dijhat(is),val_ham%dijhat(is))
       end do
    end if

    ! ndmh: skip pre-optimisation if only doing a properties calculation
    if (.not. cond_properties_only) then

       ! lr408: Now we can optimise the conduction NGWFs
       ! agrecokpt: at specified k-point?
       call ngwf_cg_optimise(total_energy, converged, out_of_runtime, cond_ham, &
            cond_denskern, edft, cond_rep, ngwf_nonsc_forces, lhxc_fine, &
            cond_ngwf_basis, proj_basis, nl_projectors, hub_proj_basis, hub, &
            mdl, hfxstate, val_rep, val_ngwf_basis, val_denskern%kern, val_ham,&
            kpt=loc_kpt)

       if (pub_on_root) then
          write(stdout,'(3a)') '+', repeat('-',78), '+'
          write(stdout,'(5a)') '|', repeat(' ',20), &
               'Conduction NGWF optimisation complete', repeat(' ',21), '|'
          write(stdout,'(3a)') '+', repeat('-',78), '+'
       end if
    else

       ! ndmh: just calculate the cond hamiltonian for properties calculation
       call hamiltonian_dens_dep_nonsc(cond_ham, cond_rep, cond_ngwf_basis, &
            mdl, hfxstate, lhxc_fine,hub,val_rep,val_ham, &
            val_denskern%kern%m(:,PUB_1K),shift_changed, &
            val_ngwf_basis = val_ngwf_basis, &         ! jd: needed for HFx
            cond_dkn = cond_denskern%kern%m(:,PUB_1K)) ! jd: needed for HFx

       out_of_runtime = .false.

    end if ! not properties only

    ! call conduction properties module
    if (.not. out_of_runtime) then
       ! agrecokpt: at specified k-point?
       call conduction_properties_calculate(val_denskern, cond_denskern, &
            val_rep, cond_rep, val_ham, cond_ham, val_ngwf_basis, &
            proj_basis, nl_projectors, hub_proj_basis, hub, cond_ngwf_basis, &
            joint_ngwf_basis, core_basis, core_wvfns, mdl, hfxstate, lhxc_fine,&
            kpt=loc_kpt)
    endif

    ! ndmh: memory deallocation for cond_ham
    call ngwf_ham_destroy(cond_ham)

    call kernel_destroy(cond_denskern)

    ! ndmh: memory deallocation for cond_rep
    do isub=1,mdl%nsub
       call data_functions_dealloc(cond_rep%ngwfs_on_grid(isub))
    end do
    call ngwf_rep_destroy(cond_rep)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &conduction_ngwf_optimise'

    return

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_cond_occupancies

      !======================================================================!
      ! This suboutine determines the number of conduction states of each    !
      ! spin.                                                                !
      !----------------------------------------------------------------------!
      ! Written by Laura Ratcliff, June 2010                                 !
      !======================================================================!

      use constants, only: DN, UP
      use kernel, only: kernel_num_bound_states
      use rundat, only: cond_num_states, pub_spin, cond_energy_range,&
          cond_energy_gap, pub_num_spins

      implicit none

      integer :: num_ngwfs

      ! jme: SPINS_DANGER: this subroutine mixes spins
      ! jme: KPOINTS_DANGER
      ! Some parts of this subroutine enforce a single k-point.
      call utils_assert(pub_num_kpoints == PUB_1K, &
           'Subroutine internal_cond_occupancies (conduction) not ready yet for more&
           & than one k-point.')

      ! tjz21: If cond_num_states is not specified, we initialise cond%n_occ to
      ! tjz21: the number of bound states in the system
      ! lr408: Switches which has more if odd number of electrons
      ! lr408: and remains the same if even number of electrons
      ! ndmh: added logic to ensure only cond_num_states used when doing
      ! ndmh: cond properties calculation
      if (cond_num_states == 0) then
         num_ngwfs=sum(val_ngwf_basis(:)%num)
         if(cond_energy_range>0.0_DP) then
            cond_rep%n_occ(:,PUB_1K) = &
                 kernel_num_bound_states(val_denskern,val_ham%ham,&
                 val_rep%overlap,val_rep%inv_overlap,num_ngwfs,&
                 val_rep%n_occ(:,PUB_1K), &
                 energy_range=cond_energy_range,&
                 energy_gap=cond_energy_gap)
         else
            cond_rep%n_occ(:,PUB_1K) = &
                 kernel_num_bound_states(val_denskern,val_ham%ham,&
                 val_rep%overlap,val_rep%inv_overlap,num_ngwfs,&
                 val_rep%n_occ(:,PUB_1K),energy_gap=cond_energy_gap)
         endif
      else if ((cond_num_extra_states /= 0).and.(.not.cond_properties_only)) then
         cond_rep%n_occ(UP,PUB_1K) = (cond_num_states + &
              cond_num_extra_states - pub_spin)/2
         if (pub_num_spins > 1) then
            cond_rep%n_occ(DN,PUB_1K) = (cond_num_states + &
                 cond_num_extra_states + pub_spin)/2
         end if
      else
         cond_rep%n_occ(UP,PUB_1K) = (cond_num_states - pub_spin)/2
         if (pub_num_spins > 1) then
            cond_rep%n_occ(DN,PUB_1K) = (cond_num_states + pub_spin)/2
         end if
      end if

    end subroutine internal_cond_occupancies

  end subroutine conduction_ngwf_optimise


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module conduction
