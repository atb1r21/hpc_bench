! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================== !
!                                                                   !
! Constrain-potential Conjugate Gradients optimisation module       !
!                                                                   !
! This module performs the optimisation (i.e.**maximisation**) of   !
! the electronic energy with respect to the constraining-potentials !
! for the constrained-DFT functionality of ONETEP.                  !
!-------------------------------------------------------------------!
! Written by Gilberto Teobaldi in December 2011 using the           !
! ngwf_cg_mod.F90 module as inspiration...                          !
!===================================================================!

module cdft_intermediate_cg

  implicit none

  private

  public :: cdft_intermediate_u_cg_optimise
  public :: cdft_intermediate_restart_read
  public :: cdft_intermediate_update_matrices
  public :: dft_nu_intermediate_restart_read
  public :: dft_nu_intermediate_update_matrices

!gibo (4 humans)
! Be aware that this module contains **3** public subroutines
! They are listed above.
! In turn, the cdft_u_cg_optimise subroutine (below) contains
! several internal (private) subroutines/functions,
! whose name start with 'internal_'
CONTAINS

    subroutine cdft_intermediate_u_cg_optimise(total_energy, cdft_converged, &
       denskern, pur_denskern, ham, &
       lhxc_fine, mu, edft, rep, ngwf_basis, hub_proj_basis, hub, &
       mdl, hfxstate, lnv_threshold, current_maxit_lnv, conv_stat, &
       dfdtau_fine, kpt)

    !==========================================================================!
    ! This subroutine maximise the total energy with respect to the            !
    ! constraining-potential of the constrained-DFT mode                       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! total_energy    (output): the total energy                               !
    ! cdft_converged  (output): whether the CDFT U optimisation is converged   !
    ! rep              (inout): NGWF Representation (functions and matrices)   !
    ! denskern         (inout): density kernel for only alpha (or beta)        !
    !                           electrons                                      !
    ! ngwf_basis       (input): Function basis type describing the NGWFs       !
    ! proj_basis       (input): Function basis type describing the nonlocal    !
    !                           pseudopotential projectors                     !
    ! hub_proj_basis   (input): Function basis type describing the Hubbard     !
    !                           projectors                                     !
    ! lhxc_fine        (inout): Local-Hartree-Exhange-Correlation potential    !
    ! ngwf_nonsc_forces(inout): Outer loop non self-consistent force correction!
    ! For Conduction NGWF optimisation only:                                   !
    ! val_rep          (input): Valence NGWF representation (optional)         !
    ! val_ngwf_basis   (input): Valence NGWF basis (optional)                  !
    ! val_dkn          (input): Valence density kernel (optional)              !
    ! conv_stat        (input)  Whether the electronic optimisation converged  !
    !--------------------------------------------------------------------------!
    ! Adapted from ngwf_cg_optimise by Gilberto Teobaldi in April 2011         !
    !==========================================================================!


    use comms, only: pub_on_root, pub_root_proc_id, comms_bcast, &
         comms_reduce,comms_barrier
    use constants, only: DP, verbose, normal, stdout, brief, &
         HARTREE_IN_EVS, max_spins
    use electronic, only: electronic_energy
    use ensemble_dft_type, only: EDFT_MODEL
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL, cdft_energy_info,&
          dft_nu_energy_info
    use hubbard_init, only: h_species
    use linalg, only: linalg_ddot
    use kernel, only: DKERN
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_output_detail, pub_print_qc, &
        pub_cdft_atom_charge, pub_cdft_atom_spin, &
        pub_cdft_group_charge_acceptor, pub_cdft_group_charge_donor, &
        pub_cdft_group_spin_acceptor, pub_cdft_group_spin_donor, &
        pub_cdft_group_charge_acceptor_u, pub_cdft_group_charge_donor_u, &
        pub_cdft_group_spin_acceptor_u, pub_cdft_group_spin_donor_u, &
        pub_cdft_group_charge_diff, pub_cdft_group_spin_diff, &
        pub_cdft_group_charge_diff_u, pub_cdft_group_spin_diff_u, &
        pub_cdft_modes, pub_maxit_cdft_u_cg, pub_cdft_cg_type, &
        pub_cdft_cg_threshold, pub_cdft_cg_max, pub_cdft_cg_max_step, &
        pub_cdft_trial_length, &
        pub_cdft_group_charge_up_only, pub_cdft_group_charge_down_only, &
        pub_cdft_write_potentials, pub_debug_on_root, pub_num_spins, &
        pub_dft_nu, pub_hubbard_unify_sites, pub_dft_nu_opt_u1_only, &
        pub_dft_nu_opt_u2_only
    use services, only: services_flush, services_line_search_parabola, &
         services_cubic_fit_maximum
    use sparse, only: SPAM3
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments (inputs)
    type(MODEL), intent(in)      :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    type(NGWF_REP), intent(in)   :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    type(FUNC_BASIS), intent(in) :: hub_proj_basis(1)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(EDFT_MODEL), intent(inout) :: edft  !gibo: new 28.12.12
    real(kind=DP), intent(in) :: lnv_threshold
    real(kind=DP), intent(out)   :: total_energy
    integer, intent(in) :: current_maxit_lnv
    ! JCW: Optional gradient of XC energy per unit volume wrt kinetic energy
    ! JCW: density.
    real(kind=DP), optional, intent(out) :: dfdtau_fine(:,:,:,:)
    ! agrecokpt
    real(kind=DP), optional, intent(in) :: kpt(3)

    ! Arguments (outputs)
    type(DKERN), intent(inout) :: denskern
    type(DKERN), intent(inout) :: pur_denskern
    type(NGWF_HAM), intent(inout) :: ham
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, 1)
    real(kind=DP), intent(inout) :: mu(max_spins)
    logical, intent(inout) :: cdft_converged   !gibo: cDFT convergence flag


    ! Local Variables
    !local-array for cdft_gradient
    real(kind=DP), allocatable, dimension(:) :: local_cdft_gradient
    !local-array for previous step cdft_gradient
    real(kind=DP), allocatable, dimension(:) :: local_previous_cdft_gradient
    !local-array for U-opt direction
    real(kind=DP), allocatable, dimension(:) :: opt_direction
    !local-array for PREVIOUS U-opt direction
    real(kind=DP), allocatable, dimension(:) :: prev_opt_direction
    !local-array for (being-)optimised U-potentials
    real(kind=DP), allocatable, dimension(:) :: local_start_cdft_u
    real(kind=DP), allocatable, dimension(:) :: delta_u
    real(kind=DP) :: total_energy_current
    real(kind=DP) :: last_n_energies(3) ! space to store the 3 most recent energies
    real(kind=DP) :: rms_gradient
    real(kind=DP) :: max_gradient
    real(kind=DP) :: previous_rms_gradient ! RMS NGWF grad of previous iteration
    real(kind=DP) :: line_search_coeff
    real(kind=DP) :: F0,F1,F2,G_init,trial_length
    real(kind=DP) :: trial_length_x2
    real(kind=DP) :: cg_coeff
    real(kind=DP) :: quadratic_coeff,cubic_coeff,rejected_quadratic_coeff
    real(kind=DP) :: predicted_functional
    real(kind=DP) :: previous_g_dot_g
    real(kind=DP) :: current_g_dot_g
    real(kind=DP) :: previous_dir_dot_g
    real(kind=DP) :: cdft_threshold

    character(len=80), allocatable, dimension(:) :: summary_lines

    integer :: iteration         ! current iteration
    integer :: cg_count          ! current number of steps since CG reset
    integer :: is         ! pdh: spin loop counter
    integer :: ierr       ! error flag
    integer :: minit      ! qoh: Minimum number of CG iterations
    integer :: hat, species
    integer :: g_counter ! 1D U-gradient array counter
    integer :: n_mode    ! cDFT-gradient mode-counter
    integer :: cdft_gradient_size

    logical :: trial2     ! ndmh: flag to perform second trial step
    logical :: retrial1   ! ndmh: flag to perform repeat first trial step
    logical :: reversing  ! ndmh: line search is going uphill
    logical :: line_search_success ! ndmh: line search fit success flag
    logical :: check_conjugacy     ! ndmh: flag for doing cg conjugacy check

    integer :: group

    integer, optional, intent(out) :: conv_stat
    ! agrecokpt
    real(kind=DP) :: loc_kpt(3)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering cdft_u_cg_optimise'

    ! Flush output
    call services_flush

    ! Start timer
    call timer_clock('cdft_u_cg_optimise',1)

    ! agrecokpt: default is Gamma point
    if (present(kpt)) then
       loc_kpt(:) = kpt(:)
    else
       loc_kpt(:) = 0.0_DP
    end if

    ! allocate workspace for U-optimisation
    !further modified to allow parsing of U_hub corrections on all atoms
    ! and cDFT-potential only on selected atoms
    if (pub_cdft_group_charge_up_only.OR.pub_cdft_group_charge_down_only) then
     ! number_cDFT_atoms
     cdft_gradient_size = mdl%par%nat_cdft_active
    else
     ! number of spin (2) x number_cDFT_atoms x number of cDFT-modes active
     cdft_gradient_size = pub_num_spins*mdl%par%nat_cdft_active*pub_cdft_modes
    endif

    !gom: set size of gradient for DFT+nu
    if (pub_dft_nu .and. .not. pub_hubbard_unify_sites) then
       if ((pub_dft_nu_opt_u1_only) .OR. (pub_dft_nu_opt_u2_only)) then
         cdft_gradient_size = pub_num_spins*mdl%par%nat_hub
       else
         cdft_gradient_size = pub_num_spins*mdl%par%nat_hub*2
       endif
    else if (pub_dft_nu .and. pub_hubbard_unify_sites) then
       if ((pub_dft_nu_opt_u1_only) .OR. (pub_dft_nu_opt_u2_only)) then
         cdft_gradient_size = pub_num_spins*hub%num_groups
       else
         cdft_gradient_size = pub_num_spins*hub%num_groups*2
       endif
    endif

    allocate(local_cdft_gradient(cdft_gradient_size),stat=ierr)
    call utils_alloc_check('cdft_u_cg_optimise', &
         'local_cdft_gradient',ierr)
    allocate(opt_direction(cdft_gradient_size),stat=ierr)
    call utils_alloc_check('cdft_u_cg_optimise', &
         'opt_direction',ierr)
    allocate(local_start_cdft_u(cdft_gradient_size),stat=ierr)
    call utils_alloc_check('cdft_u_cg_optimise', &
         'local_start_cdft_u',ierr)
    allocate(delta_u(cdft_gradient_size),stat=ierr)
    call utils_alloc_check('cdft_u_cg_optimise', &
         'delta_u',ierr)


    if (pub_cdft_cg_max > 0) then
     allocate(prev_opt_direction(cdft_gradient_size),stat=ierr)
     call utils_alloc_check('cdft_u_cg_optimise', &
          'prev_opt_direction',ierr)
    endif

     if (pub_cdft_cg_type == 'NGWF_POLAK') then
      allocate(local_previous_cdft_gradient(cdft_gradient_size), stat=ierr)
      call utils_alloc_check('cdft_u_cg_optimise', &
           'local_previous_cdft_gradient',ierr)
     endif


    ! cks: <<< parameter initialisations >>>
    cdft_threshold     = pub_cdft_cg_threshold

    ! cks: <<< variable intialisations >>>
    previous_rms_gradient = huge(1.0_DP)
    max_gradient          = huge(1.0_DP)
    line_search_coeff  = 0.15_DP
    !trial_length       = 0.1_DP  !ORIGINAL
    ! control initial trial_length from input (default 0.1)
    trial_length = pub_cdft_trial_length
    F0=0.0_DP ; F1=0.0_DP ; F2=0.0_DP ; G_init=0.0_DP ; cg_coeff=0.0_DP
    rms_gradient=1.0_DP ;  cg_count=0
    total_energy=0.0_DP ; total_energy_current=0.0_DP
    predicted_functional =0.0_DP
    current_g_dot_g  =0.0_DP
    previous_g_dot_g =0.0_DP
    line_search_success = .true.
    trial2 = .false.
    retrial1 = .false.
    reversing = .false.
    if (pub_cdft_cg_max > 0) prev_opt_direction(:) = 0.0_DP
    last_n_energies(:) = -huge(1.0_DP)
    quadratic_coeff = 0.0_DP
    rejected_quadratic_coeff = 0.0_DP
    cubic_coeff = 0.0_DP
    cdft_converged = .false.
    local_cdft_gradient = 0.0_DP

    !allocate(summary_lines(max(pub_maxit_ngwf_cg+1,2)),stat=ierr)
    allocate(summary_lines(max(pub_maxit_cdft_u_cg+1,2)),stat=ierr)
    call utils_alloc_check('cdft_u_cg_optimise','summary_lines',ierr)


    !gibo: prepare to output summary of cDFT-U optimisation
    if (pub_on_root) write(summary_lines(1),'(a80)') '|ITER|    RMS GRADIENT   &
         &|     TOTAL ENERGY    |   step   |     Epredicted  '
    summary_lines(1) = adjustl(summary_lines(1))

    !gibo: First message for brief output level
    if (pub_on_root .and. pub_output_detail == BRIEF) then
       write(stdout,'(/a)')'oooooooooooooooooooooooooooooooooooooooo&
            &oooooooooooooooooooooooooooooooooooooooo'
       write(stdout,'(a)')'ooooooooooooooooooooo &
            & cDFT U-optimisation  oooooooooooooooooooooooooooooooooooo'
       write(stdout,'(a)')'oooooooooooooooooooooooooooooooooooooooo&
            &oooooooooooooooooooooooooooooooooooooooo'
       write(stdout,'(a80)') summary_lines(1)
    end if

    !Allow blank calculations for timings purposes if pub_maxit_cdft_u_cg < 0
    minit = 1
    if (pub_maxit_cdft_u_cg < 0) minit = 0

    iteration = 0 ! ndmh: prevent unitialised variable when pub_maxit_ngwf_cg = 0

!##############################################################################
!##############################################################################

    !gibo: U-optimisation loop: START
    OPT_LOOP: do iteration=1,max(pub_maxit_cdft_u_cg,minit)

!##############################################################################
!##############################################################################


       ! cks: ++++++++++ START ITERATION HEADER +++++++++++++++++++++++++++++
       if (pub_on_root .and. pub_output_detail >= NORMAL .and. &
            pub_maxit_cdft_u_cg > 0) then
          write(stdout,'(/a)')'oooooooooooooooooooooooooooooooooooooooo&
               &oooooooooooooooooooooooooooooooooooooooo'
          write(stdout,'(a,i4.3,a)')'ooooooooooooooooooooooooooo cDFT &
                  &CG iteration ',iteration,' ooooooooooooooooooooooooooooo'
          write(stdout,'(a/)')'oooooooooooooooooooooooooooooooooooooooo&
               &oooooooooooooooooooooooooooooooooooooooo'
       end if
       call services_flush
       ! cks: ++++++ END ITERATION HEADER +++++++++++++++++++++++++++++

       CDFT_vs_DFTnu_0: if (.not.pub_dft_nu) then

           !gibo: make a local copy of the hub%cdft_u... potentials
           !====== cdft_ATOM modes: START
           CDFT_ATOM: if (pub_cdft_atom_charge) then

             g_counter = 0
             do hat = 1, mdl%par%nat_hub
               species = hub%species_number(hat)
               ! cycle loop if present species is not cDFT-active
               if (.not.h_species(species)%cdft_active) CYCLE
               g_counter = g_counter + 2
               local_start_cdft_u(g_counter-1) = hub%cdft_u_charge_up(hat)
               local_start_cdft_u(g_counter)   = hub%cdft_u_charge_down(hat)
             enddo

            elseif (pub_cdft_atom_spin) then

             g_counter = 0
             do hat = 1, mdl%par%nat_hub
               species = hub%species_number(hat)
               ! cycle loop if present species is not cDFT-active
               if (.not.h_species(species)%cdft_active) CYCLE
               g_counter = g_counter + 2
               local_start_cdft_u(g_counter-1) = hub%cdft_u_spin(hat)
               local_start_cdft_u(g_counter)   = hub%cdft_u_spin(hat)
             enddo

           endif CDFT_ATOM
           !====== cdft_ATOM modes: END

           !====== cdft_GROUP modes: START
           ! initialise gradient-component counter to zero
           g_counter = 0

           ! GROUP-CHARGE-ACCEPTOR/DONOR mode
           if (pub_cdft_group_charge_acceptor .OR. pub_cdft_group_charge_donor .OR. &
               pub_cdft_group_charge_diff) then

    !gibo-start-02.10.12
    !ORIGINAL-OK
    !         do hat = 1, mdl%par%nat_hub
    !           g_counter = g_counter + 2
    !           local_start_cdft_u(g_counter-1) = hub%cdft_u_charge_up(hat)
    !           local_start_cdft_u(g_counter)   = hub%cdft_u_charge_down(hat)
    !         enddo
    !ORIGINAL-OK
    ! modified charge_up/down-only
             if (pub_cdft_group_charge_up_only) then
                  do hat = 1, mdl%par%nat_hub
                    species = hub%species_number(hat)
                    ! cycle loop if present species is not cDFT-active
                    if (.not.h_species(species)%cdft_active) CYCLE
                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter)   = hub%cdft_u_charge_up(hat)
                  enddo

             elseif (pub_cdft_group_charge_down_only) then
                  species = hub%species_number(hat)
                  ! cycle loop if present species is not cDFT-active
                  if (.not.h_species(species)%cdft_active) CYCLE
                  do hat = 1, mdl%par%nat_hub
                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter)   = hub%cdft_u_charge_down(hat)
                  enddo

             else
                  do hat = 1, mdl%par%nat_hub
                    species = hub%species_number(hat)
                    ! cycle loop if present species is not cDFT-active
                    if (.not.h_species(species)%cdft_active) CYCLE
                    g_counter = g_counter + 2
                    local_start_cdft_u(g_counter-1) = hub%cdft_u_charge_up(hat)
                    local_start_cdft_u(g_counter)   = hub%cdft_u_charge_down(hat)
                  enddo
             endif
    ! modified charge_up/down-only
    !gibo-end-02.10.12

           endif

           ! GROUP-SPIN-ACCEPTOR/DONOR mode
           if (pub_cdft_group_spin_acceptor .OR.&
                 pub_cdft_group_spin_donor .OR. &
               pub_cdft_group_spin_diff) then

             do hat = 1, mdl%par%nat_hub
               species = hub%species_number(hat)
               ! cycle loop if present species is not cDFT-active
               if (.not.h_species(species)%cdft_active) CYCLE
               g_counter = g_counter + 2
               local_start_cdft_u(g_counter-1) = hub%cdft_u_spin(hat)
               local_start_cdft_u(g_counter)   = hub%cdft_u_spin(hat)
             enddo

           endif
           !====== cdft_GROUP modes: END
       elseif (pub_dft_nu) then

           if (.not. pub_hubbard_unify_sites) then

             g_counter = 0
             do hat = 1, mdl%par%nat_hub

                if (pub_dft_nu_opt_u1_only) then
                    ! U1 for atom "hat" and spin "ispin"
                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter) = hub%nu_u1_up(hat)

                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter) = hub%nu_u1_down(hat)

                elseif (pub_dft_nu_opt_u2_only) then
                    ! U2 for atom "hat" and spin "ispin"
                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter) = hub%nu_u2_up(hat)

                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter) = hub%nu_u2_down(hat)
              else
                    ! U1, U2 for atom hat and spin "ispin"

                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter) = hub%nu_u1_up(hat)

                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter) = hub%nu_u2_up(hat)

                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter) = hub%nu_u1_down(hat)

                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter) = hub%nu_u2_down(hat)

                endif

             enddo ! hat

           else if (pub_hubbard_unify_sites) then

              g_counter = 0
              do group = 1, hub%num_groups

                if (pub_dft_nu_opt_u1_only) then
                    ! U1 for atom "hat" and spin "ispin"
                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter)=hub%nu_u1_up(group)

                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter)=hub%nu_u1_down(group)

                elseif (pub_dft_nu_opt_u2_only) then
                    ! U2 for atom "hat" and spin "ispin"
                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter)=hub%nu_u2_up(group)

                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter)=hub%nu_u2_down(group)

                else
                    ! U1, U2 for atom hat and spin "ispin"

                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter)=hub%nu_u1_up(group)

                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter)=hub%nu_u2_up(group)

                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter)=hub%nu_u1_down(group)

                    g_counter = g_counter + 1
                    local_start_cdft_u(g_counter)=hub%nu_u2_down(group)

                endif

             enddo ! groups
          endif
       endif CDFT_vs_DFTnu_0


       !gibo: write list of the latest cDFT U-potentials into [.cdft] file
        if (pub_cdft_write_potentials) then
            if(.not.pub_dft_nu) call internal_cdft_restart_write()
            if(pub_dft_nu) call internal_dft_nu_restart_write()
        endif


       !gibo: ============ OPTIMISE DENSITY KERNEL ONLY  ===========
       ! agrecokpt: need k-point dependence in TB method?
       total_energy = electronic_energy(denskern, pur_denskern, ham, &
            lhxc_fine, mu, edft, rep, ngwf_basis, hub_proj_basis, hub, &
            mdl, hfxstate, lnv_threshold, current_maxit_lnv, &
            kernel_update=.true., conv_status=conv_stat, &
            dfdtau_fine = dfdtau_fine, kpt=loc_kpt)

       total_energy_current = total_energy
       !gibo: ============ OPTIMISE DENSITY KERNEL ONLY  ===========


       !gibo: reduce hub%cdft_gradient(:,:,:) across the procs
       if (pub_dft_nu) then
        call comms_reduce('SUM',hub%nu_gradient)
       else
        call comms_reduce('SUM',hub%cdft_gradient)
       endif


       !======= CDFT (U) GRADIENT=============== START
       !gibo: store locally the hub%cdft_gradient
       ![previously calculated in sbrtne cdft_energy_total,
       ! see hubbard_build_mod.F90]
!gibo-start-02.10.12
!ORIGINAL-OK
!       g_counter = 0
!       do n_mode = 1, pub_cdft_modes
!        do hat = 1, mdl%par%nat_hub
!         do is = 1, mdl%par%num_spins
!          g_counter = g_counter + 1
!          !local_cdft_gradient(g_counter) = hub%cdft_gradient(is, hat)
!          local_cdft_gradient(g_counter) = hub%cdft_gradient(is, hat, n_mode)
!         enddo
!        enddo
!       enddo
!ORIGINAL-OK
!modified charge_up/down_only
    CDFT_vs_DFTnu: if (.not.pub_dft_nu) then

        if (pub_cdft_group_charge_up_only) then
            n_mode = 1
            g_counter = 0
            do hat = 1, mdl%par%nat_hub
              species = hub%species_number(hat)
              ! cycle loop if present species is not cDFT-active
              if (.not.h_species(species)%cdft_active) CYCLE
              g_counter = g_counter + 1
              local_cdft_gradient(g_counter) = hub%cdft_gradient(1, hat, n_mode)  ! UP-electrons only
            enddo

        elseif (pub_cdft_group_charge_down_only) then
            n_mode = 1
            g_counter = 0
            do hat = 1, mdl%par%nat_hub
              species = hub%species_number(hat)
              ! cycle loop if present species is not cDFT-active
              if (.not.h_species(species)%cdft_active) CYCLE
              g_counter = g_counter + 1
              local_cdft_gradient(g_counter) = hub%cdft_gradient(2, hat, n_mode)  ! DOWN-electrons only
            enddo

        else
           g_counter = 0
           do n_mode = 1, pub_cdft_modes
            LATOMS: do hat = 1, mdl%par%nat_hub
              species = hub%species_number(hat)
              ! cycle loop if present species is not cDFT-active
              if (.not.h_species(species)%cdft_active) CYCLE LATOMS
             do is = 1, pub_num_spins
              g_counter = g_counter + 1
              local_cdft_gradient(g_counter) = hub%cdft_gradient(is, hat, n_mode)
             enddo
            enddo LATOMS
           enddo

        endif

    elseif(pub_dft_nu) then

        if (.not. pub_hubbard_unify_sites) then

            g_counter = 0

            do hat = 1, mdl%par%nat_hub
              do is = 1, pub_num_spins

               if (pub_dft_nu_opt_u1_only) then
                  g_counter = g_counter + 1
                  ! U1 for atom "hat" and spin "ispin"
                  local_cdft_gradient(g_counter) = &
                    hub%nu_gradient(is, hat, 1)

               elseif (pub_dft_nu_opt_u2_only) then
                  g_counter = g_counter + 1
                  ! U2 for atom "hat" and spin "ispin"
                  local_cdft_gradient(g_counter) = &
                    hub%nu_gradient(is, hat, 2)

               else
                  g_counter = g_counter + 1
                  ! U1 for atom "hat" and spin "ispin"
                  local_cdft_gradient(g_counter) = &
                    hub%nu_gradient(is, hat, 1)

                  g_counter = g_counter + 1
                  ! U2 for atom "hat" and spin "ispin"
                  local_cdft_gradient(g_counter) = &
                     hub%nu_gradient(is, hat, 2)

               endif

              enddo ! is
            enddo ! hat

        else if (pub_hubbard_unify_sites) then
            g_counter = 0

            do group = 1, hub%num_groups
              do is = 1, pub_num_spins

               if (pub_dft_nu_opt_u1_only) then
                  g_counter = g_counter + 1
                  ! U1 for atom "hat" and spin "ispin"
                  local_cdft_gradient(g_counter) = &
                    hub%nu_gradient(is, group, 1)

               elseif (pub_dft_nu_opt_u2_only) then
                  g_counter = g_counter + 1
                  ! U2 for atom "hat" and spin "ispin"
                  local_cdft_gradient(g_counter) = &
                    hub%nu_gradient(is, group, 2)

               else
                  g_counter = g_counter + 1
                  ! U1 for atom "hat" and spin "ispin"
                  local_cdft_gradient(g_counter) = &
                    hub%nu_gradient(is, group, 1)

                  g_counter = g_counter + 1
                  ! U2 for atom "hat" and spin "ispin"
                  local_cdft_gradient(g_counter) = &
                    hub%nu_gradient(is, group, 2)

               endif

              enddo ! is
            enddo ! hat

        end if
    endif CDFT_vs_DFTnu


!modified charge_up/down_only
       !======= CDFT (U) GRADIENT=============== END


       ! ndmh: check conjugacy condition on search direction
       check_conjugacy = .false.
       if (check_conjugacy) then
        !gibo(4 humans): ddot is the double precision BLAS
        !                subroutine dedicated to the dot_product between two vectors...
        previous_dir_dot_g = linalg_ddot(cdft_gradient_size,&
             local_cdft_gradient,1,prev_opt_direction,1)
       end if


       !======= TEST CDFT-OPT CONVERGENCE=============== START
       cdft_converged = internal_test_convergence()
       call comms_bcast(pub_root_proc_id,cdft_converged)
       if (cdft_converged) then
         total_energy = total_energy_current
         EXIT OPT_LOOP ! exit U-optimisation loop
       endif
       !======= TEST CDFT-OPT CONVERGENCE=============== END


       ! pdh: exit if no cDFT U-optimisation is required
       if (pub_maxit_cdft_u_cg == 0) then
          cdft_converged = .false.
          total_energy = total_energy_current
          EXIT OPT_LOOP ! exit U-optimisation loop
       end if

       ! cks: **************** FUNCTIONAL AT INITIAL POINT *********************
       !gibo: use total_energy [from last NGWF optimisation] for F0
       F0 = total_energy_current
       ! cks: ************** END FUNCTIONAL AT INITIAL POINT *******************

       ! cks: ************************ START LINE SEARCH ********************
       ! find new CG-(maximisation) direction
       call internal_find_direction

       ! cks: store direction if doing conjugate gradients
       if (pub_cdft_cg_max > 0) then
        prev_opt_direction(:) = opt_direction(:)
       endif

       if (pub_cdft_cg_type == 'NGWF_POLAK') &
            local_previous_cdft_gradient(:) = local_cdft_gradient(:)

       previous_g_dot_g = current_g_dot_g


       !===== UPDATE cDFT POTENTIALS & MATRICES ===== START
       !update cDFT U-potentials on the basis of local_start_cdft_u and delta_u
       !gom
       if (.not.pub_dft_nu) then
         call internal_cdft_update_u(trial_length)
       elseif (pub_dft_nu) then
         call internal_dft_nu_update_u(trial_length)
       endif

       !gibo: write a list of the latest cDFT U-potentials into file
       if (pub_cdft_write_potentials) then
            if(.not.pub_dft_nu) call internal_cdft_restart_write()
            if(pub_dft_nu) call internal_dft_nu_restart_write()
       endif

       ! update hub%up/down_matrices, which contain the
       ! (diagonal) cDFT Hamiltonian
       !call cdft_update_matrices(hub, hub_proj_basis)
       !gom
       if (.not.pub_dft_nu) then
         call cdft_intermediate_update_matrices(hub, hub_proj_basis(1))
       elseif (pub_dft_nu) then
         call dft_nu_intermediate_update_matrices(hub, hub_proj_basis(1))
       endif

       !===== UPDATE cDFT POTENTIALS & MATRICES ===== END


       ! cks: %%%%%%%%%%%%%%%%%% FUNCTIONAL AT TRIAL LENGTH 1 %%%%%%%%%%%%%%%%%%

       !gibo: ============ OPTIMISE DENSITY KERNEL ONLY  ===========
       ! agrecokpt: need k-point dependence in TB method?
       total_energy = electronic_energy(denskern, pur_denskern, ham, &
            lhxc_fine, mu, edft, rep, ngwf_basis, hub_proj_basis, hub, &
            mdl, hfxstate, lnv_threshold, current_maxit_lnv, &
            kernel_update=.true., conv_status=conv_stat, &
            dfdtau_fine = dfdtau_fine, kpt=loc_kpt)

       F1 = total_energy
       !gibo: ============ OPTIMISE DENSITY KERNEL ONLY  ===========

       ! cks: %%%%%%%%%%%%%%%% END FUNCTIONAL AT TRIAL LENGTH 1 %%%%%%%%%%%%%%%%


       ! cks: ===================== SEARCH BY FITTING PARABOLA =================
       !gibo: follow same line-search strategy as for Denskern (internal_quadratic_step)
       call internal_line_search_parabola(&
            quadratic_coeff, predicted_functional,line_search_success, & !output
            G_init, F0, F1, trial_length, pub_cdft_cg_max_step)              !input

       line_search_coeff = quadratic_coeff
       ! cks: ================= END SEARCH BY FITTING PARABOLA =================

       ! cks: CUBIC CUBIC CUBIC ----- SEARCH BY FITTING CUBIC ---- CUBIC CUBIC
       !gibo: inverted because we are (hopefully) dealing with a convex parabola
       !      and the gradient is positive (G_init >0 )along a good search direction
       !      (quadratic_coeff >0). Therefore we need trial2=.true. *ONLY*
       !      if G_init * quadratic_coeff < 0.
       !CUBIC: if (G_init * quadratic_coeff > 0.0_DP) then
       CUBIC: if (G_init * quadratic_coeff < 0.0_DP) then

          trial2 = .true.

          !update cDFT U-potentials from local_start_cdft_u and delta_u
          trial_length_x2 = 2.0_DP*trial_length

          !gom
          if (.not.pub_dft_nu) then
            call internal_cdft_update_u(trial_length_x2)
          elseif (pub_dft_nu) then
            call internal_dft_nu_update_u(trial_length_x2)
          endif

          trial_length = 0.5_DP*trial_length_x2

          ! update hub%up/down_matrices, which contain the
          ! (diagonal) cDFT Hamiltonian
          !gom
          if (.not.pub_dft_nu) then
            call cdft_intermediate_update_matrices(hub, hub_proj_basis(1))
          elseif (pub_dft_nu) then
            call dft_nu_intermediate_update_matrices(hub, hub_proj_basis(1))
          endif

       end if CUBIC


       ! ndmh: SECOND TRIAL STEP ----- SECOND TRIAL STEP ----- SECOND TRIAL STEP
       ! ndmh: protection against bad line search results: redo trial step if
       ! ndmh: line search result is too much bigger or smaller than trial step
       ! gibo: or if we are searching 'downhill' and went past trial step position
       if ( ((line_search_coeff < 0.05_DP*trial_length) .or. &
            (line_search_coeff > 20.0_DP*trial_length) .or. &
            (.not. line_search_success) .or. &
            (reversing .and. (line_search_coeff > trial_length))) &
            .and. .not. trial2 ) then

          retrial1 = .true.
       end if

       RETRIAL_1: if (retrial1) then

        rejected_quadratic_coeff = quadratic_coeff

        !update cDFT U-potentials on the basis of local_start_cdft_u and delta_u
        !gom
        if (.not.pub_dft_nu) then
          call internal_cdft_update_u(quadratic_coeff)
        elseif (pub_dft_nu) then
          call internal_dft_nu_update_u(quadratic_coeff)
        endif


        ! update hub%up/down_matrices, which contain the
        ! (diagonal) cDFT Hamiltonian
        !gom
        if (.not.pub_dft_nu) then
          call cdft_intermediate_update_matrices(hub, hub_proj_basis(1))
        elseif (pub_dft_nu) then
          call dft_nu_intermediate_update_matrices(hub, hub_proj_basis(1))
        endif


       end if RETRIAL_1
       ! ndmh: SECOND TRIAL STEP ----- SECOND TRIAL STEP ----- SECOND TRIAL STEP

       if (trial2 .or. retrial1) then

          ! cks: %%%%%%%%%%%%%%%% FUNCTIONAL AT TRIAL LENGTH 2 %%%%%%%%%%%%%%%%%

       !gibo: write the latest cDFT U-potentials into file
       if (pub_cdft_write_potentials) then
            if(.not.pub_dft_nu) call internal_cdft_restart_write()
            if(pub_dft_nu) call internal_dft_nu_restart_write()
       endif

          !gibo: ============ OPTIMISE DENSITY KERNEL ONLY  ===========
          ! agrecokpt: need k-point dependence in TB method?
          total_energy = electronic_energy(denskern, pur_denskern, ham, &
               lhxc_fine, mu, edft, rep, ngwf_basis, hub_proj_basis, hub, &
               mdl, hfxstate, lnv_threshold, current_maxit_lnv, &
               kernel_update=.true., conv_status=conv_stat, &
               dfdtau_fine = dfdtau_fine, kpt=loc_kpt)

          F2 = total_energy
          !gibo: ============ OPTIMISE DENSITY KERNEL ONLY  ===========

          ! cks: %%%%%%%%%%%%% END FUNCTIONAL AT TRIAL LENGTH 2 %%%%%%%%%%%%%%%%

          if (trial2) then
             call services_cubic_fit_MAXIMUM( &
                  cubic_coeff, predicted_functional, line_search_success, & !output
                  F0, F1, F2, G_init, trial_length, 2.0_DP*trial_length, &  !input
                  pub_cdft_cg_max_step)                                         !input
             line_search_coeff = cubic_coeff
          end if

          if (retrial1) then
             ! ndmh: quadratic fit at new trial length
             !gibo: follow same line-search strategy as for Denskern (internal_quadratic_step)
             call internal_line_search_parabola(&
                  quadratic_coeff, predicted_functional, &               !output
                  line_search_success, G_init, F0, F2, &                 !input
                  line_search_coeff, pub_cdft_cg_max_step)                   !input

             line_search_coeff = quadratic_coeff
          end if

          ! AAM: quadratic fit using two trial steps
          !          quadratic_coeff=services_parabolic_step( &
          !               F0,F1,F2,trial_length,2.0_DP*trial_length)
          !          line_search_coeff=quadratic_coeff
       else
          ! ndmh: no second trial step
          F2 = 0.0_DP
       end if
       ! cks: CUBIC CUBIC --------- END SEARCH BY FITTING CUBIC -------- CUBIC


       !update cDFT U-potentials on the basis of local_start_cdft_u and delta_u
       if (.not.pub_dft_nu) then
         call internal_cdft_update_u(line_search_coeff)
       elseif (pub_dft_nu) then
         call internal_dft_nu_update_u(line_search_coeff)
       endif

       ! update hub%up/down_matrices, which contain the
       ! (diagonal) cDFT Hamiltonian
       !gom
       if (.not.pub_dft_nu) then
         call cdft_intermediate_update_matrices(hub, hub_proj_basis(1))
       elseif (pub_dft_nu) then
         call dft_nu_intermediate_update_matrices(hub, hub_proj_basis(1))
       endif


       !gibo: write intermediate summary of U-optimisation
       call internal_cdft_opt_summary

       !if required, print out information on current line search
       if (pub_on_root .and. pub_output_detail >= NORMAL) &
           call internal_print_search_info

       ! cks: ************************* END LINE SEARCH ***********************

       ! ndmh: reset flags
       trial2 = .false.
       retrial1 = .false.
       reversing = .false.
       F2 = 0.0_DP

       ! cks: write line with summary info for curent iteration
       if (pub_on_root) then
          if(abs(total_energy_current)<100000.0_DP) then
             write(summary_lines(iteration+1),'(i4,f21.14,f22.14,f11.6,f22.14)')&
                iteration,rms_gradient,total_energy_current,line_search_coeff, &
                predicted_functional
          else
             write(summary_lines(iteration+1),'(i4,f21.14,f22.12,f11.6,f22.12)')&
                iteration,rms_gradient,total_energy_current,line_search_coeff, &
                predicted_functional
          end if
          if (pub_output_detail == BRIEF) then
             write(stdout,'(a80)') summary_lines(iteration+1)
          end if
       end if


!##############################################################################
!##############################################################################

    !gibo: U-optimisation loop: END
    enddo OPT_LOOP

!##############################################################################
!##############################################################################


    if (.not. cdft_converged) then
       ! cks: reset iteration number just for storing final line of calculation
       !      summary
       iteration = iteration - 1
       ! cks: print warning that calculation failed to converge
       if (pub_on_root) write(stdout,'(a,i4,a)') &
            'WARNING: maximum number of cDFT CG iterations (',pub_maxit_cdft_u_cg, &
            ') exceeded!'
    end if


    ! cks: print calculation summary
    call internal_calculation_summary

    ! cks: print quality control information
    if (pub_print_qc) call internal_qc_output

    ! gibo: write out cDFT(+U) occupancies if it hasn't already been done
    if (pub_output_detail .ne. VERBOSE) then
       if(.not.pub_dft_nu) then
       call cdft_energy_info(hub,hub_proj_basis(1))
       else if(pub_dft_nu) then
       call dft_nu_energy_info(hub,hub_proj_basis(1))
       end if
    endif

    !gibo: write out atom-resolved  summary of U-optimisation
    call internal_cdft_opt_summary


    ! deallocate workspace for U-optimisation
    deallocate(local_cdft_gradient,stat=ierr)
    call utils_dealloc_check('cdft_u_cg_optimise','local_cdft_gradient',ierr)
    deallocate(opt_direction,stat=ierr)
    call utils_dealloc_check('cdft_u_cg_optimise','opt_direction',ierr)
    deallocate(local_start_cdft_u,stat=ierr)
    call utils_dealloc_check('cdft_u_cg_optimise','local_start_cdft_u',ierr)
    deallocate(delta_u,stat=ierr)
    call utils_dealloc_check('cdft_u_cg_optimise','delta_u',ierr)

    if (pub_cdft_cg_max > 0) then
     deallocate(prev_opt_direction,stat=ierr)
     call utils_dealloc_check('cdft_u_cg_optimise', &
          'previous_opt_direction',ierr)
    endif

    if (pub_cdft_cg_type == 'NGWF_POLAK') then
     deallocate(local_previous_cdft_gradient,stat=ierr)
     call utils_dealloc_check('cdft_u_cg_optimise', &
          'local_previous_cdft_gradient',ierr)
    endif

    deallocate(summary_lines,stat=ierr)
    call utils_dealloc_check('cdft_u_cg_optimise','summary_lines',ierr)

    call timer_clock('cdft_u_cg_optimise',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving cdft_u_cg_optimise'

    ! Flush output
    call services_flush

    return

! gibo (4 humans)
! Below, the (internal) subrtuines/functions contained by the
! 'cdft_u_cg_optimise' subroutine
  contains

!##############################################################################
!##############################################################################

    subroutine internal_cdft_update_u(step)

    !==========================================================!
    ! This subroutine updates the cDFT U potentials.           !
    ! It also protects agains crazy steps which would turn     !
    ! acceptor-atoms into donor-atoms [U_old*U_new <0]         !
    ! and viceversa...                                         !
    !----------------------------------------------------------!
    ! Written for cDFT-OPT module by Gilberto Teobaldi         !
    ! in December 2011                                         !
    !==========================================================!

    use utils, only: utils_assert

    implicit none

    real(kind=DP), intent(inout) :: step

    integer :: hat, g_counter, sp, delta_counter
    real(kind=DP) :: delta_scale
    real(kind=DP), parameter :: tol_u = 1.E-9_DP
    logical :: good_step

    delta_scale = 1.0_DP

    !CHARGE only cDFT run
    !update hub%cdft_u_charge_up/down .AND. the (energy-)suggested cDFT step

    ! For cdft_atom_charge, allow cdft-potential to change sign (needed for
    ! cDFT-EMBEDDING)
    CDFT_TYPE: if ((pub_cdft_atom_charge) .AND. (pub_cdft_modes == 1)  ) then

         !scale line_search step according to delta_scale
         step = step * delta_scale

            g_counter = 0
            atoms_cdft_3: do hat = 1, mdl%par%nat_hub
              species = hub%species_number(hat)
              ! cycle loop if present species is not cDFT-active
              if (.not.h_species(species)%cdft_active) CYCLE atoms_cdft_3
              g_counter = g_counter + 2

              ! work out tentative change for cDFT potentials
              delta_u(g_counter-1) = step*opt_direction(g_counter-1)
              delta_u(g_counter)   = step*opt_direction(g_counter)

              ! update cDFT potentials
              hub%cdft_u_charge_up(hat) = local_start_cdft_u(g_counter-1) &
                       + delta_scale*delta_u(g_counter-1)
              hub%cdft_u_charge_down(hat) = local_start_cdft_u(g_counter) &
                       + delta_scale*delta_u(g_counter)

            enddo atoms_cdft_3

    elseif ((pub_cdft_group_charge_acceptor.OR.&
          pub_cdft_group_charge_donor.OR.pub_cdft_group_charge_diff) .AND.    &
          (pub_cdft_modes == 1)  ) then


      !protect against sign-inversion of U-potentials
      CHECK_DONOR_ACCEPTOR_charge: do

         good_step =  .TRUE.

         !scale line_search step according to delta_scale
         step = step * delta_scale

         ! UP-electrons only
         SPLIT_UP_DOWN_charge: if (pub_cdft_group_charge_up_only) then

            g_counter = 0
            atoms_cdft_4: do hat = 1, mdl%par%nat_hub
              species = hub%species_number(hat)
              ! cycle loop if present species is not cDFT-active
              if (.not.h_species(species)%cdft_active) CYCLE atoms_cdft_4
              !g_counter = g_counter + 2
              g_counter = g_counter + 1

              ! work out tentative change for cDFT potentials
              !delta_u(g_counter-1) = step*opt_direction(g_counter-1)
              delta_u(g_counter)   = step*opt_direction(g_counter)

              ! update cDFT potentials
              !hub%cdft_u_charge_up(hat) = local_start_cdft_u(g_counter-1) &
              !         + delta_scale*delta_u(g_counter-1)
              !hub%cdft_u_charge_down(hat) = local_start_cdft_u(g_counter) &
              !         + delta_scale*delta_u(g_counter)
              hub%cdft_u_charge_up(hat) = local_start_cdft_u(g_counter) &   ! UP-electrons only
                       + delta_scale*delta_u(g_counter)
              !hub%cdft_u_charge_down(hat) = local_start_cdft_u(g_counter) &   ! DOWN-electrons only
              !         + delta_scale*delta_u(g_counter)

              !change of sign for potential_UP of atom hat! Reduce delta_scale
              if (hub%cdft_u_charge_up(hat)*local_start_cdft_u(g_counter) < &
                   epsilon(1.0_DP) ) then
                   delta_scale = 0.1_DP*delta_scale
                   good_step = .FALSE.
                   CYCLE CHECK_DONOR_ACCEPTOR_charge
              endif

              !!change of sign for potential_DOWN of atom hat! Reduce delta_scale
              !if (hub%cdft_u_charge_down(hat)*local_start_cdft_u(g_counter) < &
              !     epsilon(1.0_DP) ) then
              !     delta_scale = 0.1_DP*delta_scale
              !     good_step = .FALSE.
              !     CYCLE CHECK_DONOR_ACCEPTOR_charge
              !endif

            enddo atoms_cdft_4


         ! DOWN-electrons only
         elseif (pub_cdft_group_charge_down_only) then

            g_counter = 0
            atoms_cdft_5: do hat = 1, mdl%par%nat_hub
              species = hub%species_number(hat)
              ! cycle loop if present species is not cDFT-active
              if (.not.h_species(species)%cdft_active) CYCLE atoms_cdft_5
              !g_counter = g_counter + 2
              g_counter = g_counter + 1

              ! work out tentative change for cDFT potentials
              !delta_u(g_counter-1) = step*opt_direction(g_counter-1)
              delta_u(g_counter)   = step*opt_direction(g_counter)

              ! update cDFT potentials
              !hub%cdft_u_charge_up(hat) = local_start_cdft_u(g_counter-1) &
              !         + delta_scale*delta_u(g_counter-1)
              !hub%cdft_u_charge_down(hat) = local_start_cdft_u(g_counter) &
              !         + delta_scale*delta_u(g_counter)
              !hub%cdft_u_charge_up(hat) = local_start_cdft_u(g_counter) &   ! UP-electrons only
              !         + delta_scale*delta_u(g_counter)
              hub%cdft_u_charge_down(hat) = local_start_cdft_u(g_counter) &   ! DOWN-electrons only
                       + delta_scale*delta_u(g_counter)

              !!change of sign for potential_UP of atom hat! Reduce delta_scale
              !if (hub%cdft_u_charge_up(hat)*local_start_cdft_u(g_counter) < &
              !     epsilon(1.0_DP) ) then
              !     delta_scale = 0.1_DP*delta_scale
              !     good_step = .FALSE.
              !     CYCLE CHECK_DONOR_ACCEPTOR_charge
              !endif

              !change of sign for potential_DOWN of atom hat! Reduce delta_scale
              if (hub%cdft_u_charge_down(hat)*local_start_cdft_u(g_counter) < &
                   epsilon(1.0_DP) ) then
                   delta_scale = 0.1_DP*delta_scale
                   good_step = .FALSE.
                   CYCLE CHECK_DONOR_ACCEPTOR_charge
              endif

            enddo atoms_cdft_5


         ! DOWN+UP electrons
         else

            g_counter = 0
            atoms_cdft_6: do hat = 1, mdl%par%nat_hub
              species = hub%species_number(hat)
              ! cycle loop if present species is not cDFT-active
              if (.not.h_species(species)%cdft_active) CYCLE atoms_cdft_6
              g_counter = g_counter + 2

              ! work out tentative change for cDFT potentials
              delta_u(g_counter-1) = step*opt_direction(g_counter-1)
              delta_u(g_counter)   = step*opt_direction(g_counter)

              ! update cDFT potentials
              hub%cdft_u_charge_up(hat) = local_start_cdft_u(g_counter-1) &
                       + delta_scale*delta_u(g_counter-1)
              hub%cdft_u_charge_down(hat) = local_start_cdft_u(g_counter) &
                       + delta_scale*delta_u(g_counter)

              !change of sign for potential_UP of atom hat! Reduce delta_scale
              if (hub%cdft_u_charge_up(hat)*local_start_cdft_u(g_counter-1) < &
                   epsilon(1.0_DP) ) then
                   delta_scale = 0.1_DP*delta_scale
                   good_step = .FALSE.
                   CYCLE CHECK_DONOR_ACCEPTOR_charge
              endif

              !change of sign for potential_DOWN of atom hat! Reduce delta_scale
              if (hub%cdft_u_charge_down(hat)*local_start_cdft_u(g_counter) < &
                   epsilon(1.0_DP) ) then
                   delta_scale = 0.1_DP*delta_scale
                   good_step = .FALSE.
                   CYCLE CHECK_DONOR_ACCEPTOR_charge
              endif

            enddo atoms_cdft_6

         endif SPLIT_UP_DOWN_charge

         if (good_step) EXIT CHECK_DONOR_ACCEPTOR_charge

         call utils_assert(delta_scale >= tol_u, &
              'Error in internal_cdft_update_u: prevented change of sign for &
              &cDFT potentials!!! Disaster may loom large...')

      enddo CHECK_DONOR_ACCEPTOR_charge


    !SPIN only cDFT run
    !update hub%cdft_u_spin .AND. the (energy-)suggested cDFT step

    ! For cdft_atom_spine, allow cdft-potential to change sign
    elseif ( (pub_cdft_atom_spin) .AND. (pub_cdft_modes == 1) ) then

         !scale line_search step according to delta_scale
         step = step * delta_scale

         g_counter = 0
         atoms_cdft_7: do hat = 1, mdl%par%nat_hub
           species = hub%species_number(hat)
           ! cycle loop if present species is not cDFT-active
           if (.not.h_species(species)%cdft_active) CYCLE atoms_cdft_7
           g_counter = g_counter + 2

           ! work-out tentative change for cDFT potentials
           delta_u(g_counter-1) = step*opt_direction(g_counter-1)
           !delta_u(g_counter)   = step*opt_direction(g_counter)

           ! update cDFT potentials
           hub%cdft_u_spin(hat) = local_start_cdft_u(g_counter-1) &
                    + delta_scale*delta_u(g_counter-1)
           ! gibo: no need to update for spin-beta
           !       [implictly accounted for in cdft_energy_info]

         enddo atoms_cdft_7


    elseif ((pub_cdft_group_spin_acceptor .OR.       &
             pub_cdft_group_spin_donor.OR.pub_cdft_group_spin_diff)  .AND. &
             (pub_cdft_modes == 1) ) then

      !protect against sign-inversion of U-potentials
      CHECK_DONOR_ACCEPTOR_spin: do

         good_step =  .TRUE.

         !scale line_search step according to delta_scale
         step = step * delta_scale

         g_counter = 0
         atoms_cdft_8: do hat = 1, mdl%par%nat_hub
           species = hub%species_number(hat)
           ! cycle loop if present species is not cDFT-active
           if (.not.h_species(species)%cdft_active) CYCLE atoms_cdft_8
           g_counter = g_counter + 2

           ! work-out tentative change for cDFT potentials
           delta_u(g_counter-1) = step*opt_direction(g_counter-1)
           !delta_u(g_counter)   = step*opt_direction(g_counter)

           ! update cDFT potentials
           hub%cdft_u_spin(hat) = local_start_cdft_u(g_counter-1) &
                    + delta_scale*delta_u(g_counter-1)
           ! gibo: no need to update for spin-beta
           !       [implictly accounted for in cdft_energy_info]

           !change of sign for spin-potential of atom hat! Reduce delta_scale
           if (hub%cdft_u_spin(hat)*local_start_cdft_u(g_counter-1) < &
                epsilon(1.0_DP) ) then
                delta_scale = 0.1_DP*delta_scale
                good_step = .FALSE.
                CYCLE CHECK_DONOR_ACCEPTOR_spin
           endif

         enddo atoms_cdft_8

         if (good_step) EXIT CHECK_DONOR_ACCEPTOR_spin

         call utils_assert(delta_scale >= tol_u, &
              'Error in internal_cdft_update_u: prevented change of sign for &
              &cDFT potentials!!! Disaster may loom large...')

      enddo CHECK_DONOR_ACCEPTOR_spin

    ! CHARGE+SPIN cDFT run
    elseif (pub_cdft_modes == 2) then

      delta_counter = INT(cdft_gradient_size/2)

      !protect against sign-inversion of U-potentials
      CHECK_DONOR_ACCEPTOR_charge_spin: do

         good_step =  .TRUE.

         !scale line_search step according to delta_scale
         step = step * delta_scale

         g_counter = 0
         atoms_cdft_9: do hat = 1, mdl%par%nat_hub
           species = hub%species_number(hat)
           ! cycle loop if present species is not cDFT-active
           if (.not.h_species(species)%cdft_active) CYCLE atoms_cdft_9
           g_counter = g_counter + 2

           ! work out tentative change for cDFT CHARGE-potentials
           delta_u(g_counter-1) = step*opt_direction(g_counter-1)
           delta_u(g_counter)   = step*opt_direction(g_counter)

           ! update cDFT CHARGE-potentials
           hub%cdft_u_charge_up(hat) = local_start_cdft_u(g_counter-1) &
                    + delta_scale*delta_u(g_counter-1)
           hub%cdft_u_charge_down(hat) = local_start_cdft_u(g_counter) &
                    + delta_scale*delta_u(g_counter)


           ! work-out tentative change for cDFT SPIN-potentials
           ! MIND that for pub_cdft_modes=2, spin-related terms are stored in
           ! [1+pub_num_spins*mdl%par%nat_hub: 2*pub_num_spins*mdl%par%nat_hub]
           delta_u(g_counter-1+delta_counter) = step*opt_direction(g_counter-1+delta_counter)
           !delta_u(g_counter+delta_counter)   = step*opt_direction(g_counter+delta_counter)

           ! update cDFT SPIN-potentials
           hub%cdft_u_spin(hat) = local_start_cdft_u(g_counter-1+delta_counter) &
                    + delta_scale*delta_u(g_counter-1+delta_counter)
           ! gibo: no need to update for spin-beta
           !       [implictly accounted for in cdft_energy_info]


           !change of sign for potential_UP of atom hat! Reduce delta_scale
           if (hub%cdft_u_charge_up(hat)*local_start_cdft_u(g_counter-1) < &
                epsilon(1.0_DP) ) then
                delta_scale = 0.1_DP*delta_scale
                good_step = .FALSE.
                CYCLE CHECK_DONOR_ACCEPTOR_charge_spin
           endif

           !change of sign for potential_DOWN of atom hat! Reduce delta_scale
           if (hub%cdft_u_charge_down(hat)*local_start_cdft_u(g_counter) < &
                epsilon(1.0_DP) ) then
                delta_scale = 0.1_DP*delta_scale
                good_step = .FALSE.
                CYCLE CHECK_DONOR_ACCEPTOR_charge_spin
           endif

           !change of sign for spin-potential of atom hat! Reduce delta_scale
           if (hub%cdft_u_spin(hat)*local_start_cdft_u(g_counter-1+delta_counter) < &
                epsilon(1.0_DP) ) then
                delta_scale = 0.1_DP*delta_scale
                good_step = .FALSE.
                CYCLE CHECK_DONOR_ACCEPTOR_charge_spin
           endif


         enddo atoms_cdft_9

         if (good_step) EXIT CHECK_DONOR_ACCEPTOR_charge_spin

         call utils_assert(delta_scale >= tol_u, &
              'Error in internal_cdft_update_u: prevented change of sign for &
              &cDFT potentials!!! Disaster may loom large...')

      enddo CHECK_DONOR_ACCEPTOR_charge_spin


    endif CDFT_TYPE


    !for group_charge/spin_acceptor/donor runs,
    !update group cDFT U-potential
    GROUP_CHARGE_ACCEPTOR: if (pub_cdft_group_charge_acceptor) then

       if (.not.pub_cdft_group_charge_down_only) then

         LOOP_atom_1_1: do hat = 1, mdl%par%nat_hub
           sp = hub%species_number(hat)
           ! cycle loop if present species is not cDFT-active
           if (.not.h_species(sp)%cdft_active) CYCLE LOOP_atom_1_1
           if (h_species(sp)%cdft_charge_acceptor) then
             pub_cdft_group_charge_acceptor_u = hub%cdft_u_charge_up(hat)
             EXIT LOOP_atom_1_1
           endif
         enddo LOOP_atom_1_1

       elseif (pub_cdft_group_charge_down_only) then

        LOOP_atom_1_2: do hat = 1, mdl%par%nat_hub
           sp = hub%species_number(hat)
           ! cycle loop if present species is not cDFT-active
           if (.not.h_species(sp)%cdft_active) CYCLE LOOP_atom_1_2
           if (h_species(sp)%cdft_charge_acceptor) then
             !pub_cdft_group_charge_acceptor_u = hub%cdft_u_charge_up(hat)
             pub_cdft_group_charge_acceptor_u = hub%cdft_u_charge_down(hat) ! UP-electrons are neglected...
             EXIT LOOP_atom_1_2
           endif
         enddo LOOP_atom_1_2

       endif

    endif GROUP_CHARGE_ACCEPTOR

    GROUP_CHARGE_DONOR: if (pub_cdft_group_charge_donor) then

       if (.not.pub_cdft_group_charge_down_only) then

         LOOP_atom_2_1: do hat = 1, mdl%par%nat_hub
           sp = hub%species_number(hat)
           ! cycle loop if present species is not cDFT-active
           if (.not.h_species(sp)%cdft_active) CYCLE LOOP_atom_2_1
           if (h_species(sp)%cdft_charge_donor) then
             pub_cdft_group_charge_donor_u = hub%cdft_u_charge_up(hat)
             EXIT LOOP_atom_2_1
           endif
         enddo LOOP_atom_2_1

       elseif (pub_cdft_group_charge_down_only) then

         LOOP_atom_2_2: do hat = 1, mdl%par%nat_hub
           sp = hub%species_number(hat)
           ! cycle loop if present species is not cDFT-active
           if (.not.h_species(sp)%cdft_active) CYCLE LOOP_atom_2_2
           if (h_species(sp)%cdft_charge_donor) then
             !pub_cdft_group_charge_donor_u = hub%cdft_u_charge_up(hat)
             pub_cdft_group_charge_donor_u = hub%cdft_u_charge_down(hat) ! UP-electrons are neglected...
             EXIT LOOP_atom_2_2
           endif
         enddo LOOP_atom_2_2

       endif

    endif GROUP_CHARGE_DONOR

    GROUP_SPIN_ACCEPTOR: if (pub_cdft_group_spin_acceptor) then
      LOOP_atom_3: do hat = 1, mdl%par%nat_hub
        sp = hub%species_number(hat)
        ! cycle loop if present species is not cDFT-active
        if (.not.h_species(sp)%cdft_active) CYCLE LOOP_atom_3
        if (h_species(sp)%cdft_spin_acceptor) then
          pub_cdft_group_spin_acceptor_u = hub%cdft_u_spin(hat)
          EXIT LOOP_atom_3
        endif
      enddo LOOP_atom_3
    endif GROUP_SPIN_ACCEPTOR

    GROUP_SPIN_DONOR: if (pub_cdft_group_spin_donor) then
      LOOP_atom_4: do hat = 1, mdl%par%nat_hub
        sp = hub%species_number(hat)
        ! cycle loop if present species is not cDFT-active
        if (.not.h_species(sp)%cdft_active) CYCLE LOOP_atom_4
        if (h_species(sp)%cdft_spin_donor) then
          pub_cdft_group_spin_donor_u = hub%cdft_u_spin(hat)
          EXIT LOOP_atom_4
        endif
      enddo LOOP_atom_4
    endif GROUP_SPIN_DONOR

    GROUP_CHARGE_DIFF: if (pub_cdft_group_charge_diff) then
          ! take just the first U-potential (they are all equal..)
          ! modified for charge_spin_up/down_only

          atoms_cdft_10: do hat = 1, mdl%par%nat_hub

             sp = hub%species_number(hat)
             ! cycle loop if present species is not cDFT-active
             if (h_species(sp)%cdft_active) then

               if (.not.pub_cdft_group_charge_down_only) then
                pub_cdft_group_charge_diff_u = ABS(hub%cdft_u_charge_up(hat))
               elseif (pub_cdft_group_charge_down_only) then
                pub_cdft_group_charge_diff_u = ABS(hub%cdft_u_charge_down(hat)) ! DOWN-electrons are neglected
               endif

               EXIT atoms_cdft_10

             endif

          enddo atoms_cdft_10

    endif GROUP_CHARGE_DIFF

    GROUP_SPIN_DIFF: if (pub_cdft_group_spin_diff) then

          atoms_cdft_11: do hat = 1, mdl%par%nat_hub

             sp = hub%species_number(hat)
             ! cycle loop if present species is not cDFT-active
             if (h_species(sp)%cdft_active) then

               ! take just the first U-potential (they are all equal..)
               pub_cdft_group_spin_diff_u = ABS(hub%cdft_u_spin(hat))

               EXIT atoms_cdft_11

             endif

          enddo atoms_cdft_11

    endif GROUP_SPIN_DIFF


    !we do not really need this since we sum_reduce the cdft_gradient
    !at each cDFT-step. However, I leave it for extra safety...
    if (pub_cdft_group_charge_acceptor) &
           call comms_bcast(pub_root_proc_id,pub_cdft_group_charge_acceptor_u)

    if (pub_cdft_group_charge_donor) &
           call comms_bcast(pub_root_proc_id,pub_cdft_group_charge_donor_u)

    if (pub_cdft_group_spin_acceptor) &
           call comms_bcast(pub_root_proc_id,pub_cdft_group_spin_acceptor_u)

    if (pub_cdft_group_spin_donor) &
           call comms_bcast(pub_root_proc_id,pub_cdft_group_spin_donor_u)

    if (pub_cdft_group_charge_diff) &
           call comms_bcast(pub_root_proc_id,pub_cdft_group_charge_diff_u)

    if (pub_cdft_group_spin_diff) &
           call comms_bcast(pub_root_proc_id,pub_cdft_group_spin_diff_u)

    end subroutine internal_cdft_update_u

!##############################################################################
!##############################################################################

    subroutine internal_dft_nu_update_u(step)

    !==========================================================!
    ! This subroutine updates the cDFT_nu U potentials.        !
    !----------------------------------------------------------!
    ! Written for cDFT-OPT module by Gilberto Teobaldi         !
    ! in December 2011                                         !
    !==========================================================!

    use utils, only: utils_assert

    implicit none

    real(kind=DP), intent(inout) :: step

    integer :: hat, g_counter, sp, delta_counter, group
    real(kind=DP) :: delta_scale
    real(kind=DP), parameter :: tol_u = 1.E-9_DP
    logical :: good_step

    delta_scale = 1.0_DP

    !scale line_search step according to delta_scale
    step = step * delta_scale

    if (.not. pub_hubbard_unify_sites) then
        g_counter = 0
        do hat = 1, mdl%par%nat_hub

           if (pub_dft_nu_opt_u1_only) then

              g_counter = g_counter + 1

              local_start_cdft_u(g_counter) = hub%nu_u1_up(hat)
              delta_u(g_counter) = step*opt_direction(g_counter)
              hub%nu_u1_up(hat) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)

              g_counter = g_counter + 1
              local_start_cdft_u(g_counter) = hub%nu_u1_down(hat)
              delta_u(g_counter) = step*opt_direction(g_counter)
              hub%nu_u1_down(hat) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)


           elseif (pub_dft_nu_opt_u2_only) then
              g_counter = g_counter + 1
              local_start_cdft_u(g_counter) = hub%nu_u2_up(hat)
              delta_u(g_counter) = step*opt_direction(g_counter)
              hub%nu_u2_up(hat) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)

              g_counter = g_counter + 1
              local_start_cdft_u(g_counter) = hub%nu_u2_down(hat)
              delta_u(g_counter) = step*opt_direction(g_counter)
              hub%nu_u2_down(hat) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)
           else
              g_counter = g_counter + 1
              local_start_cdft_u(g_counter) = hub%nu_u1_up(hat)
              delta_u(g_counter) = step*opt_direction(g_counter)
              ! update cDFT potentials
              hub%nu_u1_up(hat) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)

              g_counter = g_counter + 1
              local_start_cdft_u(g_counter) = hub%nu_u2_up(hat)
              delta_u(g_counter) = step*opt_direction(g_counter)
              hub%nu_u2_up(hat) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)

              g_counter = g_counter + 1
              local_start_cdft_u(g_counter) = hub%nu_u1_down(hat)
              delta_u(g_counter) = step*opt_direction(g_counter)
              hub%nu_u1_down(hat) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)

              g_counter = g_counter + 1
              local_start_cdft_u(g_counter) = hub%nu_u2_down(hat)
              delta_u(g_counter) = step*opt_direction(g_counter)
              hub%nu_u2_down(hat) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)

           endif

        enddo !hat
  else if (pub_hubbard_unify_sites) then
        g_counter = 0
        do group = 1, hub%num_groups

           if (pub_dft_nu_opt_u1_only) then

              g_counter = g_counter + 1

              local_start_cdft_u(g_counter) = hub%nu_u1_up(group)
              delta_u(g_counter) = step*opt_direction(g_counter)
              hub%nu_u1_up(group) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)

              g_counter = g_counter + 1
              local_start_cdft_u(g_counter) = hub%nu_u1_down(group)
              delta_u(g_counter) = step*opt_direction(g_counter)
              hub%nu_u1_down(group) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)


           elseif (pub_dft_nu_opt_u2_only) then
              g_counter = g_counter + 1
              local_start_cdft_u(g_counter) = hub%nu_u2_up(group)
              delta_u(g_counter) = step*opt_direction(g_counter)
              hub%nu_u2_up(group) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)

              g_counter = g_counter + 1
              local_start_cdft_u(g_counter) = hub%nu_u2_down(group)
              delta_u(g_counter) = step*opt_direction(g_counter)
              hub%nu_u2_down(group) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)

           else
              g_counter = g_counter + 1
              local_start_cdft_u(g_counter) = hub%nu_u1_up(group)
              delta_u(g_counter) = step*opt_direction(g_counter)
              hub%nu_u1_up(group) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)
              g_counter = g_counter + 1
              local_start_cdft_u(g_counter) = hub%nu_u2_up(group)
              delta_u(g_counter) = step*opt_direction(g_counter)
              hub%nu_u2_up(group) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)

              g_counter = g_counter + 1
              local_start_cdft_u(g_counter) = hub%nu_u1_down(group)
              delta_u(g_counter) = step*opt_direction(g_counter)
              hub%nu_u1_down(group) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)

              g_counter = g_counter + 1
              local_start_cdft_u(g_counter) = hub%nu_u2_down(group)
              delta_u(g_counter) = step*opt_direction(g_counter)
              hub%nu_u2_down(group) = local_start_cdft_u(g_counter) + &
                                  delta_scale*delta_u(g_counter)

           endif

        enddo !group

    endif

    end subroutine internal_dft_nu_update_u


!##############################################################################
!##############################################################################

    logical function internal_test_convergence()


      use linalg, only: linalg_ddot
      use rundat, only: pub_delta_e_conv, pub_cdft_max_grad, pub_cdft_elec_energy_tol, &
           pub_kernel_force_conv, maxit_lnv
      use sparse, only: sparse_is_dense

      implicit none

      logical :: energy_change_converged
      logical :: all_converged
      logical :: test_e_conv
      logical :: test_rms_grad, test_max_grad, test_dkn_conv
      real(kind=DP) :: max_energy_diff
      real(kind=DP) :: diff12, diff23

      logical :: maxit_lnv_reached
      logical :: energy_decreasing
      logical :: cdft_grad_converged
      logical :: cdft_max_grad_converged
      logical :: dkn_converged

      ! cks: Update storage of most recent three energies
      last_n_energies(1) = last_n_energies(2)
      last_n_energies(2) = last_n_energies(3)
      last_n_energies(3) = total_energy_current

      ! ndmh: Check energy change per atom vs tolerance
      !gibo:adapted to check energy change per CDFT-atom vs tolerance
      energy_change_converged = .false.
      !test_e_conv = (pub_elec_energy_tol > 0.0_DP).and.(iteration>1)
      test_e_conv = (pub_cdft_elec_energy_tol > 0.0_DP).and.(iteration>1)
      if (test_e_conv) then
         diff12 = abs(last_n_energies(1)-last_n_energies(2)) &
              / real(mdl%par%nat_hub,kind=DP)
         diff23 = abs(last_n_energies(2)-last_n_energies(3)) &
              / real(mdl%par%nat_hub,kind=DP)
         max_energy_diff = max(diff12,diff23)
         if (max_energy_diff<pub_cdft_elec_energy_tol) then
            energy_change_converged = .true.
         end if
      end if

      ! smmd: Check convergence of density kernel
      test_dkn_conv = pub_kernel_force_conv
      if (test_dkn_conv) then
         maxit_lnv_reached = (current_maxit_lnv==maxit_lnv)
         if (conv_stat.gt.0 .and. (.not.maxit_lnv_reached)) then
            dkn_converged = .false.
         else
            dkn_converged = .true.
         endif
      endif

      ! calculate rms gradient
      ! cks: parallel calculation of NGWF rms_gradient
      !gibo(4 humans): ddot is the double precision BLAS
      !                subroutine dedicated to the dot_product between two vectors...
      rms_gradient = linalg_ddot(cdft_gradient_size,&
           local_cdft_gradient,1,local_cdft_gradient,1)

      !gibo:hub%cdft_gradient has already sum-reduced (cdft_u_cg_optimise) so
      !     no need to repeat the operation here...

      ! cks: store for Fletcher-Reeves calculation of CG coeff
      current_g_dot_g = rms_gradient
      rms_gradient = sqrt(abs(rms_gradient) / &
                     real(cdft_gradient_size, kind=DP))

      !gibo: Calculate max value of cdft_gradient
      max_gradient = 0.0_DP
      test_max_grad = (pub_cdft_max_grad > 0.0_DP)
      if (test_max_grad) then
        do g_counter = 1, cdft_gradient_size
          max_gradient = max(max_gradient, &
          sqrt(local_cdft_gradient(g_counter)*local_cdft_gradient(g_counter) ) )
        enddo
        !gibo:hub%cdft_gradient has already sum-reduced (cdft_u_cg_optimise) so
        !     no need to repeat the operation here
      end if

      ! cks: Test for if allowed to use energy gain as convergence criterion
      !gibo: adapted to CDFT Ecdft **OPTIMISATION**
      energy_decreasing = .false.
      !gibo-start-03.04.12
      !energy_decreasing = pub_delta_e_conv .and. &
      !     ( last_n_energies(3) < last_n_energies(2)) .and. &
      !     ( last_n_energies(2) < last_n_energies(1))
      !gibo-end-03.04.12

      !cDFT RMS gradient convergence criterion
      test_rms_grad = (cdft_threshold > 0.0_DP)
      if (test_rms_grad) then
         cdft_grad_converged = (rms_gradient < cdft_threshold)
      else
         cdft_grad_converged = .false.
      end if

      !cDFT MAX gradient convergence criterion
      cdft_max_grad_converged = .false.
      if (test_max_grad) then
         cdft_max_grad_converged = (max_gradient < pub_cdft_max_grad)
      end if

      ! ndmh: Check all relevant criteria for convergence
      all_converged = .true.
      if ((pub_cdft_elec_energy_tol > 0.0_DP).and.(iteration<2)) &
           all_converged = .false.
      if (test_e_conv) &
           all_converged = all_converged.and.energy_change_converged
      if (test_rms_grad) &
           all_converged = all_converged.and.cdft_grad_converged
      if (test_max_grad) &
           all_converged = all_converged.and.cdft_max_grad_converged
      if (test_dkn_conv) &
           all_converged = all_converged.and.dkn_converged
      if (pub_delta_e_conv) &
           all_converged = all_converged.or.energy_decreasing

      ALLCONVERGED: if (all_converged) then

         ! cks: print details only when output is not set to brief
         if (pub_output_detail >= NORMAL) then

            if (pub_on_root) then
               write(stdout,'(/a)')'           oooooooooooooooooooooooooooo&
                    &oooooooooooooooooooooooooooo'
               write(stdout,'(a)')'           | *** cDFT optimisation &
                    &converged ***                  |'
               if (test_rms_grad) write(stdout,'(a,f19.14,a)') &
                    '           o RMS cDFT-gradient = ',rms_gradient,'              o'
               if (test_max_grad) write(stdout,'(a,f19.14,a)') &
                    '           o MAX cDFT-gradient = ',max_gradient,'              o'
               if (test_e_conv) write(stdout,'(a,f17.14,a)') &
                    '           o Maximum change in energy = ',max_energy_diff,' / atom  o'
               write(stdout,'(a)')'           o Criteria satisfied: &
                    &                                 o'
               if (test_rms_grad.and.cdft_grad_converged) write(stdout,'(a)') &
                    '           o -> RMS cDFT-gradient lower than set threshold.       o'
               if (test_max_grad.and.cdft_max_grad_converged) write(stdout,'(a)') &
                    '           o -> MAX cDFT-gradient lower than set threshold.       o'
               if (test_e_conv.and.energy_change_converged) write(stdout,'(a)') &
                    '           o -> Energy change per atom lower than set threshold.  o'
               if (test_dkn_conv.and.conv_stat.le.0) write(stdout,'(a)') &
                    '           o -> Density kernel threshold reached                  o'
               if (energy_decreasing) then
                  write(stdout,'(a)') &
                       '           o -> WARNING!!!!                                       o'
                  write(stdout,'(a)') &
                       '           o -> Energy DECREASE over the last 3 cDFT CG-steps     o'
                  write(stdout,'(a)') &
                       '           o -> Disaster may be looming large....                 o'
               endif
               write(stdout,'(a)') '           ooooooooooooooooooooooooooo&
                    &ooooooooooooooooooooooooooooo'
            end if
         end if

       internal_test_convergence = .true.

      else

         ! ndmh: print details only when output is set to verbose
         if (pub_output_detail >= VERBOSE) then

            if (pub_on_root) then
               write(stdout,'(/a)')'           oooooooooooooooooooooooooooo&
                    &oooooooooooooooooooooooooooo'
               write(stdout,'(a)')'           o *** cDFT optimisation &
                    &not converged ***              o'
               if (test_rms_grad) write(stdout,'(a,f19.14,a)') &
                    '           o RMS cDFT-gradient = ',rms_gradient,'              o'
               if (test_max_grad) write(stdout,'(a,f19.14,a)') &
                    '           o MAX cDFT-gradient = ',max_gradient,'              o'
               if (test_e_conv) write(stdout,'(a,f17.14,a)') &
                    '           o Maximum change in energy = ',max_energy_diff,' / atom  o'
               write(stdout,'(a)')'           o Criteria not satisfied: &
                    &                             o'
               if (test_rms_grad.and.(.not.cdft_grad_converged)) &
                    write(stdout,'(a)') &
                    '           o -> RMS cDFT-gradient higher than set threshold.      o'
               if (test_max_grad.and.(.not.cdft_max_grad_converged)) &
                    write(stdout,'(a)') &
                    '           o -> MAX cDFT-gradient higher than set threshold.      o'
               if (test_e_conv.and.(.not.energy_change_converged)) &
                    write(stdout,'(a)') &
                    '           o -> Energy change per atom higher than set threshold. o'
               if (test_dkn_conv.and. (conv_stat.gt.0)) &
                    write(stdout,'(a)') &
                    '           o -> Density kernel threshold unsatisfied.             o'
               if (energy_decreasing) then
                  write(stdout,'(a)') &
                       '           o -> WARNING!!!!                                       o'
                  write(stdout,'(a)') &
                       '           o -> Energy DECREASE over the last 3 cDFT CG-steps     o'
                  write(stdout,'(a)') &
                       '           o -> Disaster may be looming large....                 o'
               endif
               write(stdout,'(a)') '           ooooooooooooooooooooooooooo&
                    &ooooooooooooooooooooooooooooo'
            end if
         end if

         internal_test_convergence = .false.

      endif ALLCONVERGED


    end function internal_test_convergence

!##############################################################################
!##############################################################################

    subroutine internal_find_direction

      use linalg, only : linalg_ddot

      implicit none

      if (.not.line_search_success) then
         trial_length = trial_length * 0.5_DP
      else if (line_search_coeff > 0.0_DP) then
         trial_length = &
              max(sqrt(trial_length * line_search_coeff),epsilon(1.0_DP))
      end if
      trial_length = max(trial_length,0.0001_DP)

      if (pub_on_root .and. (pub_output_detail >= NORMAL) ) write(stdout,'(/a)') &
           'ooooooooooooooooooooooooooooooo cDFT line search &
           &ooooooooooooooooooooooooooooooo'


      ! calculate CG coefficient
      if ( iteration > 1 ) then
         if ((cg_count >= pub_cdft_cg_max) .or. (.not.line_search_success)) then
            ! cks: reset CG after "cg_max" steps
            ! ndmh: or after a fitting failure
            ! lr408: or if conduction shift has just been updated
            cg_coeff = 0.0_DP
            cg_count = 0
            if (pub_on_root  .and. (.not.line_search_success) .and. &
                 (pub_output_detail >= NORMAL) ) write(stdout,'(a)') &
                 'Resetting cDFT CG'
         else
            ! cks <<< FLETCHER >>>
            if (pub_cdft_cg_type == 'NGWF_FLETCHER') then

               ! cks: original Fletcher-Reeves formula (cheaper in memory)
               if (abs(previous_g_dot_g) > epsilon(1.0_DP)) then
                  cg_coeff = current_g_dot_g / previous_g_dot_g
                  ! cks: protection from crazy coefficients
                  if (abs(cg_coeff) > 2.0_DP ) then
                     if (pub_on_root  .and.(pub_output_detail >= NORMAL) ) &
                          write(stdout,'(a,f8.4,a)') &
                          'WARNING: cDFT Fletcher-Reeves CG coeff too large (',&
                          cg_coeff, ') - setting to zero'
                     cg_coeff = 0.0_DP
                     cg_count = 0
                  end if
               else
                  if (pub_on_root .and.(pub_output_detail >= NORMAL) ) &
                       write(stdout,*) ' previous_g_dot_g=', &
                       previous_g_dot_g, &
                       'CG coeffient set to zero in cdft_cg_u_optimise'
                  cg_coeff = 0.0_DP
                  cg_count = 0
               end if

               ! cks: <<< POLAK >>>
            else if (pub_cdft_cg_type == 'NGWF_POLAK') then

               ! POLAK FORMULA
               !gibo: calculate the POLAK-RIBIERE CG coefficient
               cg_coeff = internal_cdft_polak_cg_coeff()

            end if

            ! cks: re-initialise the periodic reset process if cg_coeff was zero
            ! cks: otherwise increase cg_count
            if (cg_coeff == 0.0_DP) then
               cg_count = 0
            else
               cg_count = cg_count + 1
            end if

         end if
      else
         cg_coeff = 0.0_DP
      endif

      ! Find search direction
      if (pub_cdft_cg_max > 0) then
         !****WARNING****
         ! Since we want to maximise E_cDFT we need to go uphill
         ! [i.e. follow to gradient i.e. no minus the gradient
         !opt_direction(:) = -local_cdft_gradient(:) + &  ! OK for minimisation
         opt_direction(:) =  local_cdft_gradient(:) + &  ! OK for MAXIMISATION
                                 cg_coeff * prev_opt_direction(:)

      else
         ! cks:steepest descents
         !opt_direction(:) = -local_cdft_gradient(:) ! OK for minimisation
         opt_direction(:) =  local_cdft_gradient(:) ! OK for MAXIMISATION
      endif


      ! Slope of energy in search direction
      G_init = linalg_ddot(cdft_gradient_size, &
              local_cdft_gradient,1,opt_direction,1)

      !gibo: no need to collect anything since we are working on
      !      the sum-reduced cdft_gradient

      ! take action in case of positive slope along search direction
      !gibo: for cDFT, take action against NEGATIVE slope along search direction
      if (G_init < 0.0_DP) then

         if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then
            write(stdout,'(a,e16.6)') &
               'WARNING: slope along cDFT search direction is negative:', G_init
            write(stdout,'(a)') '         Resetting conjugate gradients!'
         end if
         !opt_direction(:) = -local_cdft_gradient(:) ! OK for minimisation
         opt_direction(:) =  local_cdft_gradient(:) ! OK for MAXIMISATION

         cg_count = 0
         G_init = linalg_ddot(cdft_gradient_size, &
                 local_cdft_gradient,1,opt_direction,1)

         !gibo: no need to collect anything since we are working on
         !      the sum-reduced cdft_gradient

         !gibo: for cDFT, act against NEGATIVE slope along search direction
         if (G_init < 0.0_DP) then
            if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then
               write(stdout,'(a)') 'WARNING: slope along cDFT search direction &
                    &is still negative.'
               write(stdout,'(a)') '         Reversing search direction!!'
            end if
            opt_direction(:) = -opt_direction(:)
            G_init = -G_init
            ! ndmh: if searching 'uphill', (or downhill for cDFT)
            ! ndmh: always re-check final step
            ! ndmh: before accepting, to avoid accepting very bad steps.
            reversing = .true.
         end if

      end if

      call services_flush


    end subroutine internal_find_direction

!##############################################################################
!##############################################################################

   function internal_cdft_polak_cg_coeff()

     !===============================================================!
     ! This subroutine calculates the conjugate gradient coefficient !
     ! according to the Polak formula:                               !
     !  b_{i+1}=g_{i}^T * [g_{i}-g{i-1}]/[G_{i-1}^T * G_{i-1}]       !
     !                                                               !
     ! for the cDFT optimisation of the U-opt (1D) array             !
     !---------------------------------------------------------------!
     ! Adapted from services_polak_cg_coeff by Gilberto Teobaldi     !
     ! in December 2011                                              !
     !===============================================================!

     use comms, only: pub_on_root, comms_reduce
     use constants, only: DP, stdout, NORMAL
     use linalg, only : linalg_ddot
     use rundat, only : pub_output_detail

     real(kind=DP)  :: internal_cdft_polak_cg_coeff
     real(kind=DP)  :: eps
     real(kind=DP)  :: denominator

     eps =epsilon(1.0_DP)

     denominator =  linalg_ddot(cdft_gradient_size,&
              local_cdft_gradient,1,local_cdft_gradient,1)

     if ( abs(denominator) .gt. eps ) then
        internal_cdft_polak_cg_coeff = &
           linalg_ddot(cdft_gradient_size,local_cdft_gradient,1, &
           (local_cdft_gradient - local_previous_cdft_gradient),1) / denominator

     else
        if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
             write(stdout,'(a/a)') &
             'WARNING: zero denominator in internal_cdft_polak_cg_coeff', &
             '         Setting to zero'
        internal_cdft_polak_cg_coeff = 0.0_DP
     endif

     if ( abs(internal_cdft_polak_cg_coeff).gt.(2.0_DP) ) then

        if (pub_on_root .and. pub_output_detail >= NORMAL) &
             write(stdout,'(a,f11.5,a)') &
             'WARNING: internal_cdft_polak_cg_coeff too large &
             &(',internal_cdft_polak_cg_coeff,'). Setting to zero'
        internal_cdft_polak_cg_coeff = 0.0_DP
     endif

   end function internal_cdft_polak_cg_coeff

!##############################################################################
!##############################################################################

  subroutine internal_line_search_parabola( &
       min_length, min_value, success, &                       ! output
       grad_at_zero, value_at_zero, value_at_trial_length, &   ! input
       trial_length,max_step)                                  ! input

    !==================================================================!
    ! Do line search by fitting a parabola.                            !
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/6/2001.                   !
    ! Improved by Chris-Kriton Skylaris on 20/6/2003 and on 26/4/2004. !
    ! Return value to indicate success added by Nick Hine on 04/10/2008!
    ! Adapted to cDFT-line search from services_line_search_parabola   !
    ! by Gilberto Teobaldi in April 2011                               !
    !==================================================================!

    use comms, only : pub_on_root
    use constants, only: DP, stdout, NORMAL
    use rundat, only: pub_output_detail
    implicit none

    real(kind=DP), intent(out) :: min_length, min_value
    logical, intent(out) :: success
    real(kind=DP), intent(in) :: trial_length, grad_at_zero, value_at_zero &
         , value_at_trial_length,max_step

    ! cks: internal declarations
    real(kind=DP), parameter :: eps=epsilon(1.0_DP)
    real(kind=DP) :: a,b             ! Coefficients of quadratic equation
    real(kind=DP) :: linear_term    ! b * trial_step
    !real(kind=DP) :: quadratic_term ! a * trial_step**2
    real(kind=DP) :: q_diff         ! value_at_trial_length - value_at_zero
    real(kind=DP) :: trial_length_2 ! square of trial_length

    ! Determine coefficients of quadratic equation
    linear_term =  grad_at_zero*trial_length  ! gibo: b*x1

    q_diff = value_at_trial_length - value_at_zero   ! y1- y0

    !=== numerically unstable if qdiff gets dangerously close to linear_term
    !    quadratic_term = q_diff - linear_term            ! y1 - y0 - b*x1
    !    a = quadratic_term / (trial_length*trial_length) !(y1-y0-b*x1)/(x1*x1)
    !alternative...
    trial_length_2 = trial_length*trial_length
    ! less efficient, yet numerically more stable if q_diff~linear_term
    ! [which will happen close to the maximum...]
    a = q_diff/trial_length_2 - linear_term/trial_length_2
    !=====

    b = grad_at_zero

    if (abs(a) > abs(value_at_zero)*eps*eps) then

       min_length = -0.5_DP * b / a ! -b/2a

       min_value = value_at_zero + 0.5_DP * b * min_length ! c-b^2/(4a)

       success = .true.

    else

      if ((pub_output_detail >= VERBOSE) .and. pub_on_root) &
           write(stdout,'(a)') '  cDFT-quadratic term too small: &
           &increase trial_length/proceed to cubic fitting'

      min_length = 0.1_DP
      success = .false.

    end if

    ! cks: protection from crazy line search coefficients
    ! ndmh: changed to use keyword for max line search step
    if (abs(min_length) > max_step) then
       if (pub_on_root .and. pub_output_detail >= NORMAL) &
            !write(stdout,'(a/a,f11.5,a)') &
            write(stdout,'(a/a,f18.8,a)') &
            'WARNING in internal_line_search_parabola:', &
            '  quadratic step (', min_length, &
            ') too large - setting to safe value'
            write(stdout,'(a,f14.8,a)') &
            'value at zero =', value_at_zero
            write(stdout,'(a,f14.8,a)') &
            'value at trial_length =', value_at_trial_length
            write(stdout,'(a,g14.8,a)') &
            'q_diff =', q_diff
            write(stdout,'(a,g14.8,a)') &
            'linear_term=', linear_term
            !write(stdout,'(a,g14.8,a)') &
            !'quadratic_term=', quadratic_term
            write(stdout,'(a,f14.8,a)') &
            'trial_length=', trial_length
            write(stdout,'(a,f14.8,a)') &
            'max_step=', max_step
            write(stdout,'(a,g14.8,a)') &
            'fitted "a"=', a
            write(stdout,'(a,g14.8,a)') &
            'fitted "b"=', b
       min_length = 0.15_DP
       success = .false.
    end if

  end subroutine internal_line_search_parabola

!##############################################################################
!##############################################################################


    subroutine internal_print_search_info

      use rundat, only: pub_cdft_max_grad

      implicit none

      write(stdout,'(/a,f22.14)')   'cDFT:RMS gradient                = ',rms_gradient
      if (pub_cdft_max_grad>0.0_DP) write(stdout,'(a,f22.14)') &
           'cDFT MAX gradient                = ',max_gradient
      write(stdout,'(a,f22.6)')    'cDFT:Trial step length           = ',trial_length
      write(stdout,'(a,f22.11)')   'cDFT:Gradient along search dir.  = ',G_init
      if (check_conjugacy) then
         write(stdout,'(a,f22.11)')   'cDFT:Conjugacy Test              = ', &
              previous_dir_dot_g
      end if
      write(stdout,'(a,f22.14)')   'cDFT:Functional at step 0        = ',F0
      write(stdout,'(a,f22.14)')   'cDFT:Functional at step 1        = ',F1
      if (trial2) then
         write(stdout,'(a,f22.14)')'cDFT:Functional at step 2        = ',F2
         write(stdout,'(a,f22.14)')'cDFT:Functional predicted        = ', &
              predicted_functional
         write(stdout,'(a,f22.6)') 'cDFT:Rejected quadratic step     = ', &
              quadratic_coeff
         write(stdout,'(a,f22.6)') 'cDFT:Selected cubic step         = ',cubic_coeff
      else if (retrial1) then
         write(stdout,'(a,f22.6)') 'cDFT:Rejected quadratic step     = ', &
              rejected_quadratic_coeff
         write(stdout,'(a,f22.14)')'cDFT:Functional at new step 1    = ',F2
         write(stdout,'(a,f22.14)')'cDFT:Functional predicted        = ', &
              predicted_functional
         write(stdout,'(a,f22.6)') 'cDFT:Selected quadratic step     = ', &
              quadratic_coeff
      else
         write(stdout,'(a,f22.14)')'cDFT:Functional predicted        = ', &
              predicted_functional
         write(stdout,'(a,f22.6)') 'cDFT:Selected quadratic step     = ', &
              quadratic_coeff
      endif
      if ( pub_cdft_cg_max > 0 ) then
         write(stdout,'(a,f22.6)')    'cDFT:Conjugate gradients coeff.  = ',cg_coeff
      endif
      write(stdout,'(a)')'oLooooooooooooooooooooooooo cDFT line search &
           &finished oooooooooooooooooooooooooo'
      write(stdout,'(a)')'                                               '

    end subroutine internal_print_search_info

!##############################################################################
!##############################################################################

   !gibo: print out intermediate summary of cDFT-optimisation
   !      with atom resolved break-down of the gradient
   subroutine internal_cdft_opt_summary

    implicit none

    real(kind=DP), parameter :: zero=0.0_DP

    integer :: delta_counter

    if (pub_on_root .and. pub_output_detail >= NORMAL .and. &
         pub_maxit_cdft_u_cg > 0) then

        write(stdout,'(/a)')'oSoooooooooooooooooooooooooooooooooooooo&
             &oooooooooooooooooooooooooooooooooooooooo'
        write(stdout,'(a,i4.3,a,f24.14/)') &
                 & '    cDFT summary, CG iteration: ', iteration, &
                 & '       Energy (Ha): ', total_energy_current
        write(stdout,'(a,f15.9)')  'TOTAL RMS cDFT-GRADIENT(UP+DOWN):', &
               rms_gradient
        write(stdout,'(a,f15.9)')  'MAX cDFT-GRADIENT(UP+DOWN):      ', &
               max_gradient
    endif

    if (pub_on_root .and. pub_output_detail >= VERBOSE .and. &
         pub_maxit_cdft_u_cg > 0) then

        cDFT_vs_DFTnu_1: if (.not.pub_dft_nu) then

            !CHARGE only cDFT run
            CHARGE_SPIN_CDFT: if ((pub_cdft_atom_charge .OR.            &
                                   pub_cdft_group_charge_acceptor .OR.  &
                                   pub_cdft_group_charge_donor    .OR.  &
                                   pub_cdft_group_charge_diff)   .AND.  &
                                    (pub_cdft_modes == 1)) then

               write(stdout,'(/a)') '  site    label       Uq(UP)    Uq(DOWN) &
               &      Us      dE/dUq(UP)  dE/dUq(DOWN) '


    !gibo-start-02.10.12
               ! UP-electrons only
               SPLIT_UP_DOWN_charge_1: if (pub_cdft_group_charge_up_only) then

                  SPECIES_1: do species = 1, mdl%par%num_hub_species
                     ! cycle loop if present species is not cDFT-active
                     if (.not.h_species(species)%cdft_active) CYCLE SPECIES_1
                     do hat = 1 , mdl%par%nat_hub

                      if (hub%species_number(hat) == species) then
                          !g_counter = 2*hat
                          g_counter = hat
                          write(stdout,'(i5,7x,a,5(f12.6))') hat, &
                            & h_species(species)%hub_species, &
                            & local_start_cdft_u(g_counter) * HARTREE_IN_EVS, &
                            & zero,   &
                            & zero,   &
                            & local_cdft_gradient(g_counter), &
                            & zero
                      endif

                     enddo
                  enddo SPECIES_1


               ! DOWN-electrons only
               elseif (pub_cdft_group_charge_down_only) then


                  SPECIES_2: do species = 1, mdl%par%num_hub_species
                     ! cycle loop if present species is not cDFT-active
                     if (.not.h_species(species)%cdft_active) CYCLE SPECIES_2
                     do hat = 1 , mdl%par%nat_hub

                      if (hub%species_number(hat) == species) then
                          g_counter = 2*hat
                          write(stdout,'(i5,7x,a,5(f12.6))') hat, &
                            & h_species(species)%hub_species, &
                            & zero,   &
                            & local_start_cdft_u(g_counter)   * HARTREE_IN_EVS, &
                            & zero,   &
                            & zero,   &
                            & local_cdft_gradient(g_counter)
                      endif

                     enddo
                  enddo SPECIES_2


               else
               ! UP+ DOWN-electrons
    !gibo-end-02.10.12
                  SPECIES_3: do species = 1, mdl%par%num_hub_species
                     ! cycle loop if present species is not cDFT-active
                     if (.not.h_species(species)%cdft_active) CYCLE SPECIES_3
                     do hat = 1 , mdl%par%nat_hub

                      if (hub%species_number(hat) == species) then
                          g_counter = 2*hat
                          write(stdout,'(i5,7x,a,5(f12.6))') hat, &
                            & h_species(species)%hub_species, &
                            & local_start_cdft_u(g_counter-1) * HARTREE_IN_EVS, &
                            & local_start_cdft_u(g_counter)   * HARTREE_IN_EVS, &
                            & zero                                        , &
                            & local_cdft_gradient(g_counter-1), &
                            & local_cdft_gradient(g_counter)
                      endif

                     enddo
                  enddo SPECIES_3

               endif SPLIT_UP_DOWN_charge_1


            !SPIN only cDFT run
            elseif ((pub_cdft_atom_spin .OR.            &
                     pub_cdft_group_spin_acceptor .OR.  &
                     pub_cdft_group_spin_donor    .OR.  &
                     pub_cdft_group_spin_diff)   .AND.  &
                      (pub_cdft_modes == 1)) then

             write(stdout,'(/a)') '  site    label       Uq(UP)    Uq(DOWN) &
             &      Us      dE/dUs(UP)  dE/dUs(DOWN) '

               SPECIES_4: do species = 1, mdl%par%num_hub_species
                  ! cycle loop if present species is not cDFT-active
                  if (.not.h_species(species)%cdft_active) CYCLE SPECIES_4
                  do hat = 1 , mdl%par%nat_hub

                   if (hub%species_number(hat) == species) then
                       g_counter = 2*hat
                       write(stdout,'(i5,7x,a,5(f12.6))') hat, &
                         & h_species(species)%hub_species, &
                         & zero                                        , &
                         & zero                                        , &
                         & local_start_cdft_u(g_counter-1) * HARTREE_IN_EVS, &
                         & local_cdft_gradient(g_counter-1), &
                         & local_cdft_gradient(g_counter)
                   endif

                  enddo
               enddo SPECIES_4

            ! CHARGE+SPIN cDFT run
            elseif (pub_cdft_modes == 2) then

               ! shift to acess spin-terms within the same loop
               delta_counter = INT(cdft_gradient_size/2)

               write(stdout,'(/a)') '  site    label       Uq(UP)    Uq(DOWN)      &
                  & Us      dE/dUq(UP)  dE/dUq(DOWN)    dE/dUs(UP)  dE/dUs(DOWN) '

               SPECIES_5: do species = 1, mdl%par%num_hub_species
                  ! cycle loop if present species is not cDFT-active
                  if (.not.h_species(species)%cdft_active) CYCLE SPECIES_5
                  do hat = 1 , mdl%par%nat_hub

                   if (hub%species_number(hat) == species) then
                       g_counter = 2*hat
                       write(stdout,'(i5,7x,a,7(f12.6))') hat, &
                         & h_species(species)%hub_species, &
                         & local_start_cdft_u(g_counter-1) * HARTREE_IN_EVS, &
                         & local_start_cdft_u(g_counter)   * HARTREE_IN_EVS, &
                         & local_start_cdft_u(g_counter-1  +delta_counter) * &
                         & HARTREE_IN_EVS, &
                         & local_cdft_gradient(g_counter-1),                 &
                         & local_cdft_gradient(g_counter),                   &
                         & local_cdft_gradient(g_counter-1 +delta_counter),  &
                         & local_cdft_gradient(g_counter   +delta_counter)
                   endif

                  enddo
               enddo SPECIES_5


            endif CHARGE_SPIN_CDFT

        elseif (pub_dft_nu) then

           if (pub_dft_nu_opt_u1_only) then
           write(stdout,'(/a)') '  site    label     &
             & U1(UP)     U1(DOWN)    dE/dU1(UP)  dE/dU1(DOWN) '

           elseif (pub_dft_nu_opt_u2_only) then
           write(stdout,'(/a)') '  site    label     &
             & U2(UP)     U2(DOWN)    dE/dU2(UP)  dE/dU2(DOWN) '

           else
           write(stdout,'(/a)') '  site    label       U1(UP)    U1(DOWN)      &
             & U2(UP)     U2(DOWN)      dE/dU1(UP)  dE/dU1(DOWN)    dE/dU2(UP) &
             & dE/dU2(DOWN) '
           endif

           if (.not. pub_hubbard_unify_sites) then
               SPECIES_6: do species = 1, mdl%par%num_hub_species

                      do hat = 1 , mdl%par%nat_hub

                       if (hub%species_number(hat) == species) then

                        if (pub_dft_nu_opt_u1_only) then
                           g_counter = 2*hat
                           write(stdout,'(i5,7x,a,8(f12.6))') hat, &
                             & h_species(species)%hub_species, &
                             & local_start_cdft_u(g_counter-1)*HARTREE_IN_EVS,&
                             & local_start_cdft_u(g_counter)  *HARTREE_IN_EVS,&
                             & local_cdft_gradient(g_counter-1),&
                             & local_cdft_gradient(g_counter)

                        elseif (pub_dft_nu_opt_u2_only) then
                           g_counter = 2*hat
                           write(stdout,'(i5,7x,a,8(f12.6))') hat, &
                             & h_species(species)%hub_species, &
                             & local_start_cdft_u(g_counter-1)*HARTREE_IN_EVS,&
                             & local_start_cdft_u(g_counter)  *HARTREE_IN_EVS,&
                             & local_cdft_gradient(g_counter-1),&
                             & local_cdft_gradient(g_counter)

                        else
                           g_counter = 4*hat
                           write(stdout,'(i5,7x,a,8(f12.6))') hat, &
                             & h_species(species)%hub_species, &
                             & local_start_cdft_u(g_counter-3)*HARTREE_IN_EVS,&
                             & local_start_cdft_u(g_counter-1)*HARTREE_IN_EVS,&
                             & local_start_cdft_u(g_counter-2)*HARTREE_IN_EVS,&
                             & local_start_cdft_u(g_counter)  *HARTREE_IN_EVS,&
                             & local_cdft_gradient(g_counter-3),&
                             & local_cdft_gradient(g_counter-1),&
                             & local_cdft_gradient(g_counter-2),&
                             & local_cdft_gradient(g_counter)

                        endif
                       endif

                      enddo !hat

               enddo SPECIES_6
          else if (pub_hubbard_unify_sites) then
              SPECIES_7: do species = 1, mdl%par%num_hub_species

                  do group = 1 , hub%num_groups

                      if (h_species(species)%nu_group == group) then

                          if (pub_dft_nu_opt_u1_only) then
                             g_counter = 2*group
                             write(stdout,'(i5,7x,a,8(f12.6))') group, &
                               & h_species(species)%hub_species, &
                               & local_start_cdft_u(g_counter-1) *&
                                 HARTREE_IN_EVS, &
                               & local_start_cdft_u(g_counter)   *&
                                 HARTREE_IN_EVS, &
                               & local_cdft_gradient(g_counter-1),&
                               & local_cdft_gradient(g_counter)

                          elseif (pub_dft_nu_opt_u2_only) then
                             g_counter = 2*group
                             write(stdout,'(i5,7x,a,8(f12.6))') group, &
                               & h_species(species)%hub_species, &
                               & local_start_cdft_u(g_counter-1) *&
                                 HARTREE_IN_EVS, &
                               & local_start_cdft_u(g_counter)   *&
                                 HARTREE_IN_EVS, &
                               & local_cdft_gradient(g_counter-1),&
                               & local_cdft_gradient(g_counter)

                          else
                             g_counter = 4*group
                             write(stdout,'(i5,7x,a,8(f12.6))') group, &
                               & h_species(species)%hub_species, &
                               & local_start_cdft_u(g_counter-3) *&
                                 HARTREE_IN_EVS, &
                               & local_start_cdft_u(g_counter-1) *&
                                 HARTREE_IN_EVS, &
                               & local_start_cdft_u(g_counter-2) *&
                                 HARTREE_IN_EVS, &
                               & local_start_cdft_u(g_counter)   *&
                                 HARTREE_IN_EVS, &
                               & local_cdft_gradient(g_counter-3),&
                               & local_cdft_gradient(g_counter-1),&
                               & local_cdft_gradient(g_counter-2),  &
                               & local_cdft_gradient(g_counter)

                         endif

                      endif

                  enddo !group
              enddo SPECIES_7
           endif
        endif cDFT_vs_DFTnu_1


    endif

    if (pub_on_root .and. pub_output_detail >= NORMAL .and. &
         pub_maxit_cdft_u_cg > 0) then

        write(stdout,'(a/)')'oooooooooooooooooooooooooooooooooooooooo&
             &oooooooooooooooooooooooooooooooooooooooo'

    endif

    call services_flush

   end subroutine internal_cdft_opt_summary

!##############################################################################
!##############################################################################

   ! cks: print a concise summary of the calculation
   subroutine internal_calculation_summary

     implicit none

     integer :: sumrow
     integer :: lastrow

     lastrow = max(iteration+1,2)

     if (pub_on_root) then
        ! ndmh: adapt number of digits in total energy to suit its magnitude
        if(abs(total_energy)<100000.0_DP) then
           write(summary_lines(lastrow),'(i4, f21.14, f22.14, a)') &
                iteration, rms_gradient, total_energy, '  <-- CG'
        else
           write(summary_lines(lastrow),'(i4, f21.14, f22.12, a)') &
                iteration, rms_gradient, total_energy, '  <-- CG'
        end if

        if (pub_output_detail >= NORMAL) then
           write(stdout,'(/20x,a)') 'ooooo CDFT OPTIMISATION SUMMARY ooooo'
           do sumrow=1,lastrow
              write(stdout,'(a80)') summary_lines(sumrow)
           end do
        else
           write(stdout,'(a80)') summary_lines(iteration+1)
        endif

     end if

   end subroutine internal_calculation_summary

!##############################################################################
!##############################################################################

    subroutine internal_qc_output

      !====================================================!
      ! This subroutine prints out quality control info in !
      ! a form that can be easily accessed and compared.   !
      !----------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 12/03/2005.    !
      ! Adapted to cDFT by Gilberto Teobaldi on 7/12/2011  !
      !====================================================!

      use comms, only: comms_barrier
      use constants, only: stdout
      use utils, only: utils_qc_print ! Jacek's revision

      ! Jacek's revision === start
      if (pub_on_root) then
         write(stdout,'(a1)')' '
         call utils_qc_print('cDFT iterations',iteration)
         call utils_qc_print('cDFT total_energy',total_energy_current)
         call utils_qc_print('cDFT rms_gradient',rms_gradient)
      endif
      ! Jacek's revision === end

      !!gibo: previous verbose output
      !if (pub_on_root) then
      !   write(stdout,'(a1)')' '
      !   write(stdout, '(a30,   i9)')'<QC>        [cDFT iterations]:', &
      !        iteration
      !   write(stdout, '(a30,f18.8)')'<QC>      [cDFT total_energy]:', &
      !        total_energy_current
      !   write(stdout, '(a30,f18.8)')'<QC>                [cDFT F0]:', &
      !        F0
      !   write(stdout, '(a30,f18.8)')'<QC>                [cDFT F1]:', &
      !        F1
      !   write(stdout, '(a30,f18.8)')'<QC>                [cDFT F2]:', &
      !        F2
      !   write(stdout, '(a30,f18.8)')'<QC>[cDFT predicted_functinl]:', &
      !        predicted_functional
      !   write(stdout,'(a30,f22.12)')'<QC>      [cDFT rms_gradient]:', &
      !        rms_gradient
      !   write(stdout,'(a30,f22.12)')'<QC>            [cDFT G_init]:', &
      !        G_init
      !   write(stdout, '(a30,f16.6)')'<QC>   [cDFT quadratic_coeff]:', &
      !        quadratic_coeff
      !   write(stdout, '(a30,f16.6)')'<QC>       [cDFT cubic_coeff]:', &
      !        cubic_coeff
      !   write(stdout, '(a30,f16.6)')'<QC>          [cDFT cg_coeff]:', &
      !        cg_coeff
      !   write(stdout, '(a30,f16.6)')'<QC>      [cDFT trial_length]:', &
      !        trial_length
      !endif
      !!gibo: previous verbose output

      call comms_barrier


    end subroutine internal_qc_output

!##############################################################################
!##############################################################################

    subroutine internal_cdft_restart_write()

    !==========================================================================!
    ! This subroutine write the cDFT U-potentials                              !
    ! [i.e. the hub%cdft_u_charge_up/down or hub%cdft_u_spin arrays]           !
    ! to (.cdft) file                                                          !
    !--------------------------------------------------------------------------!
    ! Written by Gilberto Teobaldi in December 2011                            !
    !==========================================================================!

      use constants,       only: stdout
      use comms,           only: pub_on_root
      use rundat,          only: pub_rootname, pub_debug_on_root
      use utils,           only: utils_unit, utils_open_unit_check, &
             utils_close_unit_check

      implicit none

      ! Local variables
      integer           :: io_unit,io_stat
      integer           :: species, hat, g_counter
      character(len=80) :: filename
      real(KIND=DP), parameter :: zero= 0.0_DP

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering internal_cdft_restart_write'

      if (pub_on_root) then

         ! Find available unit specifier
         io_unit = utils_unit()
         write(filename,'(a,a)') trim(pub_rootname),'.cdft'

         open(unit=io_unit,iostat=io_stat,file=filename,status='REPLACE',&
              form='FORMATTED',position='REWIND',action='WRITE')
         call utils_open_unit_check('internal_cdft_restart_write',&
                                     filename,io_stat)

         ! write cDFT type for current run
         CDFT_MODE: if (pub_cdft_modes == 1) then

            if (pub_cdft_atom_charge) then
              write(io_unit,'(a/)') 'cDFT TYPE: CDFT_ATOM_CHARGE'

            elseif (pub_cdft_atom_spin) then
              write(io_unit,'(a/)') 'cDFT TYPE: CDFT_ATOM_SPIN'

            elseif (pub_cdft_group_charge_acceptor .AND. &
                    (.not.pub_cdft_group_charge_donor)) then
              write(io_unit,'(a/)') 'cDFT TYPE: CDFT_GROUP_CHARGE_ACCEPTOR'

            elseif (pub_cdft_group_charge_donor .AND.    &
                    (.not.pub_cdft_group_charge_acceptor)) then
              write(io_unit,'(a/)') 'cDFT TYPE: CDFT_GROUP_CHARGE_DONOR'

            elseif (pub_cdft_group_charge_donor .AND.    &
                    pub_cdft_group_charge_acceptor) then
              write(io_unit,'(a/)') 'cDFT TYPE: CDFT_GROUP_CHARGE_ACCEPTOR/&
                                    &DONOR [i.e. CHARGE-DIFFERENCE]'

            elseif (pub_cdft_group_spin_acceptor .AND.   &
                    (.not.pub_cdft_group_spin_donor)) then
              write(io_unit,'(a/)') 'cDFT TYPE: CDFT_GROUP_SPIN_ACCEPTOR'

            elseif (pub_cdft_group_spin_donor .AND.      &
                    (.not.pub_cdft_group_spin_acceptor)) then
              write(io_unit,'(a/)') 'cDFT TYPE: CDFT_GROUP_SPIN_DONOR'

            elseif (pub_cdft_group_spin_donor .AND.      &
                    pub_cdft_group_spin_acceptor) then
              write(io_unit,'(a/)') 'cDFT TYPE: CDFT_GROUP_SPIN_ACCEPTOR/DONOR &
                                    &[i.e. SPIN-DIFFERENCE]'

            elseif (pub_cdft_group_charge_diff) then
              write(io_unit,'(a/)') 'cDFT TYPE: CDFT_GROUP_CHARGE_DIFF'

            elseif (pub_cdft_group_spin_diff) then
              write(io_unit,'(a/)') 'cDFT TYPE: CDFT_GROUP_SPIN_DIFF'

            endif

         elseif (pub_cdft_modes == 2) then
              ! to be possibly extended with all the combinations...
              write(io_unit,'(a/)') 'cDFT TYPE: CDFT_GROUP_CHARGE+SPIN_[ACCEPTOR&
                                     &/DONOR/DIFFERENCE]'

         endif CDFT_MODE

         write(io_unit,'(a/)') '  site    label       l        Z        Uq(UP) &
                &    Uq(DOWN)       Us         Nup        Ndown      Nup-Ndown'

        !CHARGE only cDFT run
        CHARGE_SPIN_CDFT: if ((pub_cdft_atom_charge .OR.            &
                               pub_cdft_group_charge_acceptor .OR.  &
                               pub_cdft_group_charge_donor    .OR.  &
                               pub_cdft_group_charge_diff)   .AND.  &
                                (pub_cdft_modes == 1)) then

          do species = 1, mdl%par%num_hub_species
            do hat = 1 , mdl%par%nat_hub

             if (hub%species_number(hat) == species) then
                 g_counter = 2*hat
                 write(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') hat,    &
                   & h_species(species)%hub_species,               &
                   & h_species(species)%hub_ang_mom,               &
                   & h_species(species)%hub_charge,                &
                   & hub%cdft_u_charge_up(hat)   * HARTREE_IN_EVS, &
                   & hub%cdft_u_charge_down(hat) * HARTREE_IN_EVS, &
                   & zero,                                         &
                   & h_species(species)%cdft_target_up,            &
                   & h_species(species)%cdft_target_down,          &
                   & h_species(species)%cdft_target_spin
             endif

            enddo
          enddo


        !SPIN only cDFT run
        elseif ((pub_cdft_atom_spin .OR.            &
                 pub_cdft_group_spin_acceptor .OR.  &
                 pub_cdft_group_spin_donor    .OR.  &
                 pub_cdft_group_spin_diff)   .AND.  &
                  (pub_cdft_modes == 1)) then

          do species = 1, mdl%par%num_hub_species
            do hat = 1 , mdl%par%nat_hub

             if (hub%species_number(hat) == species) then
                 g_counter = 2*hat
                 write(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') hat,    &
                   & h_species(species)%hub_species,               &
                   & h_species(species)%hub_ang_mom,               &
                   & h_species(species)%hub_charge,                &
                   & zero,                                         &
                   & zero,                                         &
                   & hub%cdft_u_spin(hat)        * HARTREE_IN_EVS, &
                   & h_species(species)%cdft_target_up,            &
                   & h_species(species)%cdft_target_down,          &
                   & h_species(species)%cdft_target_spin
             endif

            enddo
          enddo

        ! CHARGE+SPIN cDFT run
        elseif (pub_cdft_modes == 2) then

          do species = 1, mdl%par%num_hub_species
            do hat = 1 , mdl%par%nat_hub

             if (hub%species_number(hat) == species) then
                 g_counter = 2*hat
                 write(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') hat,    &
                   & h_species(species)%hub_species,               &
                   & h_species(species)%hub_ang_mom,               &
                   & h_species(species)%hub_charge,                &
                   & hub%cdft_u_charge_up(hat)   * HARTREE_IN_EVS, &
                   & hub%cdft_u_charge_down(hat) * HARTREE_IN_EVS, &
                   & hub%cdft_u_spin(hat)        * HARTREE_IN_EVS, &
                   & h_species(species)%cdft_target_up,            &
                   & h_species(species)%cdft_target_down,          &
                   & h_species(species)%cdft_target_spin
             endif

            enddo
          enddo


        endif CHARGE_SPIN_CDFT

        close(io_unit,iostat=io_stat)
        call utils_open_unit_check('internal_cdft_restart_write',&
                                     filename,io_stat)

      endif

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving internal_cdft_restart_write'

    end subroutine internal_cdft_restart_write

!##############################################################################
!##############################################################################

    subroutine internal_dft_nu_restart_write()

    !=================================================================!
    ! This subroutine write the DFTnu U-potentials                    !
    ! [i.e. the hub%cdft_u_charge_up/down or hub%cdft_u_spin arrays]  !
    ! to (.cdft) file                                                 !
    !-----------------------------------------------------------------!
    ! Written by Glenn Moynihan in 2015 based on                      !
    ! internal_cdft_restart_write() written by Gilberto Teobaldi      !
    !=================================================================!

      use constants,       only: stdout
      use comms,           only: pub_on_root
      use rundat,          only: pub_rootname, pub_debug_on_root, &
                                 pub_hubbard_unify_sites
      use utils,           only: utils_unit, utils_open_unit_check, &
                                 utils_close_unit_check

      implicit none

      ! Local variables
      integer           :: group
      integer           :: io_unit,io_stat
      integer           :: species, hat, g_counter
      character(len=80) :: filename
      real(KIND=DP), parameter :: zero= 0.0_DP


      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering internal_dft_nu_restart_write'
      if (pub_on_root) then

         ! Find available unit specifier
         io_unit = utils_unit()
         write(filename,'(a,a)') trim(pub_rootname),'.dft_nu'

         open(unit=io_unit,iostat=io_stat,file=filename,status='REPLACE',&
              form='FORMATTED',position='REWIND',action='WRITE')
         call utils_open_unit_check('internal_cdft_restart_write',&
                                     filename,io_stat)

         ! write DFT_nu mode for current run
         DFTnu_MODE: if (pub_dft_nu_opt_u1_only) then
            write(io_unit,'(a/)') 'DFT_nu TYPE: dft_nu_opt_u1_only'

         elseif (pub_dft_nu_opt_u2_only) then
            write(io_unit,'(a/)') 'DFT_nu TYPE: dft_nu_opt_u2_only'

         else
            write(io_unit,'(a/)') 'DFT_nu TYPE: optimise U1 and U2'

         endif DFTnu_MODE


         write(io_unit,'(a/)') '   site    label     l      Z      U1(UP)    &
                &U1(DOWN)     U2(UP)    U2(DOWN)    Nup    Ndown'

         if (.not. pub_hubbard_unify_sites) then
              do species = 1, mdl%par%num_hub_species
                do hat = 1 , mdl%par%nat_hub

                 if (hub%species_number(hat) == species) then
                     write(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') hat,    &
                       & h_species(species)%hub_species,               &
                       & h_species(species)%hub_ang_mom,               &
                       & h_species(species)%hub_charge,                &
                       & hub%nu_u1_up(hat)   * HARTREE_IN_EVS,         &
                       & hub%nu_u1_down(hat) * HARTREE_IN_EVS,         &
                       & hub%nu_u2_up(hat)   * HARTREE_IN_EVS,         &
                       & hub%nu_u2_down(hat) * HARTREE_IN_EVS,         &
                       & hub%nu_target_n_up(hat)+hub%nu_nu_up(hat),    &
                       & hub%nu_target_n_down(hat)+hub%nu_nu_down(hat)
                 endif
                enddo
              enddo
         else if (pub_hubbard_unify_sites) then
              do species = 1, mdl%par%num_hub_species
                do group = 1 ,hub%num_groups

                 if (h_species(species)%nu_group == group) then
                     write(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') group,    &
                       & h_species(species)%hub_species,               &
                       & h_species(species)%hub_ang_mom,               &
                       & h_species(species)%hub_charge,                &
                       & hub%nu_u1_up(group)   * HARTREE_IN_EVS,         &
                       & hub%nu_u1_down(group) * HARTREE_IN_EVS,         &
                       & hub%nu_u2_up(group)   * HARTREE_IN_EVS,         &
                       & hub%nu_u2_down(group) * HARTREE_IN_EVS,         &
                       & hub%nu_target_n_up(group)+hub%nu_nu_up(group),  &
                       & hub%nu_target_n_down(group)+hub%nu_nu_down(group)
                 endif

                enddo
              enddo

         endif

        close(io_unit,iostat=io_stat)
        call utils_open_unit_check('internal_cdft_restart_write',&
                                     filename,io_stat)

      endif

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving internal_dft_nu_restart_write'

    end subroutine internal_dft_nu_restart_write

!##############################################################################
!##############################################################################


  end subroutine cdft_intermediate_u_cg_optimise

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

    subroutine cdft_intermediate_update_matrices(hub, hub_proj_basis)

    !==========================================================================!
    ! This subroutine updates the cDFT[+U] SPAM3 matrices                      !
    !--------------------------------------------------------------------------!
    ! Written for cDFT-OPT module,and based on                                 !
    ! hubbard_build_matrices by Gilberto Teobaldi in Dec. 2011                 !
    !==========================================================================!
    ! Potentially useful comment:                                              !
    ! for cDFT run, we neglect hub_spin_splitting and hub_alpha terms.         !
    ! This allows us to safely update directly hub%down_matrix and             !
    ! hub%down_matrix with the udpated cDFT U-potentials, as we do below...    !
    ! Modified to remove pub_par by Robert Charlton, 25/10/2018.               !
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use hubbard_build, only: HUBBARD_MODEL
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: &
        pub_cdft_atom_charge, pub_cdft_atom_spin, &
        pub_cdft_group_charge_acceptor, pub_cdft_group_charge_donor, &
        pub_cdft_group_spin_acceptor, pub_cdft_group_spin_donor, &
        pub_cdft_group_charge_diff, pub_cdft_group_spin_diff, &
        pub_cdft_modes, pub_debug_on_root
    use sparse, only: SPAM3, sparse_put_element, &
         sparse_get_element, sparse_get_par

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    integer :: sp, hat_on_proc, hat, theatom
    integer :: hub_proj
    !integer :: is, ierr
    real(kind=DP) :: u_charge_up, u_charge_down, u_spin_up, u_spin_down
    type(PARAL_INFO), pointer :: par

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_build_matrices'

    ! rc2013: get the parallel strategy for this SPAM3
    call sparse_get_par(par, hub%up_matrix)

    ! gibo: atom_charge .OR. one group_charge mode only
    CHARGE_CDFT: if (pub_cdft_atom_charge .OR. pub_cdft_group_charge_diff .OR. &
        (pub_cdft_group_charge_acceptor.OR.pub_cdft_group_charge_donor) ) then


    ! ddor: Loop over Hubbard atoms on my proc
    cdft_atoms_1 : do hat_on_proc = 1, &
         par%num_hub_atoms_on_proc(pub_my_proc_id)

       hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
       sp = hub%species_number(hat)
       theatom = par%distr_atom(hub%orig(hat))

       ! cycle loop if present species is not cDFT-active
       if (.not.h_species(sp)%cdft_active) CYCLE cdft_atoms_1

       u_charge_up   = hub%cdft_u_charge_up(hat)
       u_charge_down = hub%cdft_u_charge_down(hat)

       ! ddor: Loop over Hubbard projectors on my proc
       hubprojs_cdft_1: do hub_proj = hub_proj_basis%first_on_atom(theatom), &
            hub_proj_basis%first_on_atom(theatom) + &
            hub_proj_basis%num_on_atom(theatom) - 1

          call sparse_put_element(u_charge_up,&
               hub%up_matrix,hub_proj,hub_proj)
          call sparse_put_element(u_charge_down,&
               hub%down_matrix,hub_proj,hub_proj)

       end do hubprojs_cdft_1

    end do cdft_atoms_1

    endif CHARGE_CDFT

    ! gibo: atom_spin .OR. one group_spin mode only
    SPIN_CDFT: if (pub_cdft_atom_spin .OR. pub_cdft_group_spin_diff .OR. &
        (pub_cdft_group_spin_acceptor.OR.pub_cdft_group_spin_donor) ) then

       ! depending on how many cdft-modes are active, either populate or augment
       ! the cdft_potential matrices with spin-terms.
       ! Spin-only single-cDFT mode case
       CDFT_SINGLE_MULTIPLE: if (pub_cdft_modes == 1) then

          ! ddor: Loop over Hubbard atoms on my proc
          cdft_atoms_2 : do hat_on_proc = 1, &
               par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
             sp = hub%species_number(hat)
             theatom = par%distr_atom(hub%orig(hat))

             ! cycle loop if present species is not cDFT-active
             if (.not.h_species(sp)%cdft_active) CYCLE cdft_atoms_2

             ! MIND THE SIGNS: Ecdft_spin = Us*[Nup -Ndown - DN_target]
             ! thus V_up = Us, and V_down = -Us
             u_spin_up   =  hub%cdft_u_spin(hat)
             u_spin_down = -hub%cdft_u_spin(hat)

             ! ddor: Loop over Hubbard projectors on my proc
             hubprojs_cdft_2: do hub_proj = hub_proj_basis%first_on_atom(theatom), &
                  hub_proj_basis%first_on_atom(theatom) + &
                  hub_proj_basis%num_on_atom(theatom) - 1

                call sparse_put_element(u_spin_up,&
                     hub%up_matrix,hub_proj,hub_proj)
                call sparse_put_element(u_spin_down,&
                     hub%down_matrix,hub_proj,hub_proj)

             end do hubprojs_cdft_2

          end do cdft_atoms_2

       ! multiple-cDFt modes case
       elseif (pub_cdft_modes == 2) then

          ! ddor: Loop over Hubbard atoms on my proc
          cdft_atoms_3 : do hat_on_proc = 1, &
               par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
             sp = hub%species_number(hat)
             theatom = par%distr_atom(hub%orig(hat))

             ! cycle loop if present species is not cDFT-active
             if (.not.h_species(sp)%cdft_active) CYCLE cdft_atoms_3

             ! ddor: Loop over Hubbard projectors on my proc
             hubprojs_cdft_3: do hub_proj = hub_proj_basis%first_on_atom(theatom), &
                  hub_proj_basis%first_on_atom(theatom) + &
                  hub_proj_basis%num_on_atom(theatom) - 1

                ! step_1: retrieve (charge-)potential terms
                call sparse_get_element(u_spin_up,&
                     hub%up_matrix,hub_proj,hub_proj)
                call sparse_get_element(u_spin_down,&
                     hub%down_matrix,hub_proj,hub_proj)

                ! step_2: AUGMENT (with spin-contributions) (charge-)potential terms
                ! MIND THE SIGNS: Ecdft_spin = Us*[Nup -Ndown - DN_target]
                ! thus V_up = Us, and V_down = -Us
                u_spin_up   = u_spin_up   + hub%cdft_u_spin(hat)
                u_spin_down = u_spin_down - hub%cdft_u_spin(hat)

                ! step_3: update potential matrices with augmented terms
                call sparse_put_element(u_spin_up,&
                     hub%up_matrix,hub_proj,hub_proj)
                call sparse_put_element(u_spin_down,&
                     hub%down_matrix,hub_proj,hub_proj)

             end do hubprojs_cdft_3

          end do cdft_atoms_3

       endif CDFT_SINGLE_MULTIPLE

    endif SPIN_CDFT

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_build_matrices'

    end subroutine cdft_intermediate_update_matrices

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

    subroutine cdft_intermediate_restart_read(hub,mdl)

    !==========================================================================!
    ! This subroutine reads the cDFT U-potential from (.cdft) file and update  !
    ! the hub%cdft_u_charge_up/down or hub%cdft_u_spin arrays                  !
    !--------------------------------------------------------------------------!
    ! Written by Gilberto Teobaldi in December 2011                            !
    ! Modified to remove pub_par by Robert Charlton, 25/10/2018.               !
    !==========================================================================!

      use constants,       only: DP, stdout, HARTREE_IN_EVS
      use comms,           only: pub_on_root, pub_root_proc_id, comms_bcast
      use hubbard_build,   only: HUBBARD_MODEL
      use hubbard_init, only: h_species
      use model_type, only: MODEL
      use rundat, only: pub_rootname, &
        pub_cdft_atom_charge, pub_cdft_atom_spin, &
        pub_cdft_group_charge_acceptor, pub_cdft_group_charge_donor,    &
        pub_cdft_group_spin_acceptor, pub_cdft_group_spin_donor,        &
        pub_cdft_group_charge_diff, pub_cdft_group_spin_diff,           &
        pub_cdft_group_charge_acceptor_u, pub_cdft_group_charge_donor_u,&
        pub_cdft_group_spin_acceptor_u, pub_cdft_group_spin_donor_u,    &
        pub_cdft_group_charge_diff_u, pub_cdft_group_spin_diff_u,       &
        pub_cdft_modes, pub_cdft_continuation, pub_task,                    &
        pub_cdft_group_charge_up_only, pub_cdft_group_charge_down_only, &
        pub_debug_on_root
      use utils,           only: utils_unit, utils_open_unit_check, &
                                 utils_close_unit_check, utils_abort

      implicit none

      ! Arguments
      type(HUBBARD_MODEL), intent(inout) :: hub
      type(MODEL), intent(in) :: mdl

      ! Local variables
      integer           :: io_unit,io_stat
      integer           :: species, hat, sp, dummy_int
      character(len=80) :: filename
      ! ars: dummies
      character(len=4)  :: dummy_id
      integer           :: dummy_l
      real(kind=dp)     :: dummy_c
      real(kind=dp)     :: dummy_u_s, dummy_s_target
      real(kind=dp)     :: dummy_u_q_up, dummy_u_q_down
      real(kind=dp)     :: dummy_n_up, dummy_n_down

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering cdft_restart_read'

      if (pub_on_root) then

         ! Find available unit specifier
         io_unit = utils_unit()
         write(filename,'(a,a)') trim(pub_rootname),'.cdft'

         open(unit=io_unit,iostat=io_stat,file=filename,&
              access='SEQUENTIAL',form='FORMATTED',position='REWIND', &
              action='READ')


         RESTART_OK: if (io_stat == 0) then

           read(io_unit,'(/,/,/)')  ! skip the 4-line header of .cdft file


        !CHARGE only cDFT run
        CHARGE_SPIN_CDFT: if ((pub_cdft_atom_charge .OR.            &
                               pub_cdft_group_charge_acceptor .OR.  &
                               pub_cdft_group_charge_donor    .OR.  &
                               pub_cdft_group_charge_diff)   .AND.  &
                                (pub_cdft_modes == 1)) then


             do species = 1, mdl%par%num_hub_species
               do hat = 1 , mdl%par%nat_hub

                if (hub%species_number(hat) == species) then
                    read(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') dummy_int, &
                      & dummy_id,                                       &
                      & dummy_l,                                        &
                      & dummy_c,                                        &
                      & dummy_u_q_up,                                   &
                      & dummy_u_q_down,                                 &
                      & dummy_u_s,                                      &
                      & dummy_n_up,                                     &
                      & dummy_n_down,                                   &
                      & dummy_s_target

                    ! internally convert U-potentials in Hartree
!gibo-start-02.12.10
                    if (.not.pub_cdft_group_charge_down_only) &
!gibo-end-02.12.10
                       hub%cdft_u_charge_up(hat)  = dummy_u_q_up/HARTREE_IN_EVS
!gibo-start-02.12.10
                    if (.not.pub_cdft_group_charge_up_only) &
!gibo-end-02.12.10
                       hub%cdft_u_charge_down(hat)= dummy_u_q_down/HARTREE_IN_EVS
                endif

               enddo
             enddo


        !SPIN only cDFT run
        elseif ((pub_cdft_atom_spin .OR.            &
                 pub_cdft_group_spin_acceptor .OR.  &
                 pub_cdft_group_spin_donor    .OR.  &
                 pub_cdft_group_spin_diff)   .AND.  &
                  (pub_cdft_modes == 1)) then

             do species = 1, mdl%par%num_hub_species
               do hat = 1 , mdl%par%nat_hub

                if (hub%species_number(hat) == species) then
                    read(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') dummy_int, &
                      & dummy_id,                                       &
                      & dummy_l,                                        &
                      & dummy_c,                                        &
                      & dummy_u_q_up,                                   &
                      & dummy_u_q_down,                                 &
                      & dummy_u_s,                                      &
                      & dummy_n_up,                                     &
                      & dummy_n_down,                                   &
                      & dummy_s_target

                    ! internally convert U-potentials in Hartree
                    hub%cdft_u_spin(hat)   = dummy_u_s/HARTREE_IN_EVS
                endif

               enddo
             enddo

        ! CHARGE+SPIN cDFT run
        elseif (pub_cdft_modes == 2) then

          do species = 1, mdl%par%num_hub_species
            do hat = 1 , mdl%par%nat_hub

                if (hub%species_number(hat) == species) then
                    read(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') dummy_int, &
                      & dummy_id,                                       &
                      & dummy_l,                                        &
                      & dummy_c,                                        &
                      & dummy_u_q_up,                                   &
                      & dummy_u_q_down,                                 &
                      & dummy_u_s,                                      &
                      & dummy_n_up,                                     &
                      & dummy_n_down,                                   &
                      & dummy_s_target


                    ! internally convert U-potentials in Hartree
                    hub%cdft_u_charge_up(hat)  = dummy_u_q_up/HARTREE_IN_EVS
                    hub%cdft_u_charge_down(hat)= dummy_u_q_down/HARTREE_IN_EVS
                    hub%cdft_u_spin(hat)   = dummy_u_s/HARTREE_IN_EVS
             endif

            enddo
          enddo


           endif CHARGE_SPIN_CDFT

           close(unit=io_unit,iostat=io_stat)
           call utils_close_unit_check('cdft_restart_read',filename,io_stat)

         else
           write(stdout,'(/,3a)') 'ERROR cdft_restart_read could &
             &not be opened "', trim(filename),'"'
           write(stdout,'(/a)') '...aborting...'

           if (pub_cdft_continuation) then
             write(stdout,'(/,a)') 'ERROR: it looks like you have set&
               &"cdft_continuation=.TRUE." in your input file withouth having &
               &the .cdft file in your work-directory'
             write(stdout,'(/,a)') 'SOLUTION: comment out &
               &"cdft_continuation=.TRUE." in your input file'

           elseif ((pub_task=='GEOMETRYOPTIMIZATION')  .OR. &
                   (pub_task=='MOLECULARDYNAMICS')     .OR. &
                   (pub_task=='TRANSITIONSTATESEARCH') .OR. &
                   (pub_task=='PROPERTIES')            .OR. &
                   (pub_task=='COND')                  .OR. &
                   (pub_task=='COUPLINGS')             .OR. &
                   (pub_task=='PROPERTIES_COND') ) then
             call utils_abort('A .cdft file must be present &
               &to perform a constrained-DFT GEOMETRYOPTIMIZATION / &
               &MOLECULARDYNAMICS / TRANSITIONSTATESEARCH / PROPERTIES / &
               &COND / PROPERTIES_COND task. SOLUTION: set &
               &"task : SINGLEPOINT" in your input file to obtain an optimised &
               & .cdft file first.')
           end if
         endif RESTART_OK

      endif


      !broadcast updated U-potentials
      !CHARGE only cDFT run
      if ((pub_cdft_atom_charge .OR.            &
           pub_cdft_group_charge_acceptor .OR.  &
           pub_cdft_group_charge_donor    .OR.  &
           pub_cdft_group_charge_diff)   .AND.  &
            (pub_cdft_modes == 1)) then

!gibo-start-02.12.10
        if (.not.pub_cdft_group_charge_down_only) &
!gibo-end-02.12.10
         call comms_bcast(pub_root_proc_id, hub%cdft_u_charge_up(:))

!gibo-start-02.12.10
        if (.not.pub_cdft_group_charge_up_only) &
!gibo-end-02.12.10
        call comms_bcast(pub_root_proc_id, hub%cdft_u_charge_down(:))


      !SPIN only cDFT run
      elseif ((pub_cdft_atom_spin .OR.            &
               pub_cdft_group_spin_acceptor .OR.  &
               pub_cdft_group_spin_donor    .OR.  &
               pub_cdft_group_spin_diff)   .AND.  &
                (pub_cdft_modes == 1)) then

        call comms_bcast(pub_root_proc_id, hub%cdft_u_spin(:))


      ! CHARGE+SPIN cDFT run
      elseif (pub_cdft_modes == 2) then
        call comms_bcast(pub_root_proc_id, hub%cdft_u_charge_up(:))
        call comms_bcast(pub_root_proc_id, hub%cdft_u_charge_down(:))
        call comms_bcast(pub_root_proc_id, hub%cdft_u_spin(:))

      endif

    !for group_charge/spin_acceptor/donor runs,
    !update group cDFT U-potential
    GROUP_CHARGE_ACCEPTOR: if (pub_cdft_group_charge_acceptor) then
! gibo-start-02.10.12
! ORIGINAL-OK
!      LOOP_atom_1: do hat = 1, mdl%par%nat_hub
!        sp = hub%species_number(hat)
!        if (h_species(sp)%cdft_charge_acceptor) then
!          pub_cdft_group_charge_acceptor_u = hub%cdft_u_charge_up(hat)
!          EXIT LOOP_atom_1
!        endif
!      enddo LOOP_atom_1
! ORIGINAL-OK
! modified group_charge_UP/DOWN_only
       if (.not.pub_cdft_group_charge_down_only) then

         LOOP_atom_1_1: do hat = 1, mdl%par%nat_hub
           sp = hub%species_number(hat)
           if (h_species(sp)%cdft_charge_acceptor) then
             pub_cdft_group_charge_acceptor_u = hub%cdft_u_charge_up(hat)
             EXIT LOOP_atom_1_1
           endif
         enddo LOOP_atom_1_1

       elseif (pub_cdft_group_charge_down_only) then

        LOOP_atom_1_2: do hat = 1, mdl%par%nat_hub
           sp = hub%species_number(hat)
           if (h_species(sp)%cdft_charge_acceptor) then
             !pub_cdft_group_charge_acceptor_u = hub%cdft_u_charge_up(hat)
             pub_cdft_group_charge_acceptor_u = hub%cdft_u_charge_down(hat) ! UP-electrons are neglected...
             EXIT LOOP_atom_1_2
           endif
         enddo LOOP_atom_1_2

       endif
! modified group_charge_UP/DOWN_only
! gibo-end-02.10.12
    endif GROUP_CHARGE_ACCEPTOR

    GROUP_CHARGE_DONOR: if (pub_cdft_group_charge_donor) then
! gibo-start-02.10.12
! ORIGINAL-OK
!      LOOP_atom_2: do hat = 1, mdl%par%nat_hub
!        sp = hub%species_number(hat)
!        if (h_species(sp)%cdft_charge_donor) then
!          pub_cdft_group_charge_donor_u = hub%cdft_u_charge_up(hat)
!          EXIT LOOP_atom_2
!        endif
!      enddo LOOP_atom_2
! ORIGINAL-OK
! modified group_charge_UP/DOWN_only
       if (.not.pub_cdft_group_charge_down_only) then

         LOOP_atom_2_1: do hat = 1, mdl%par%nat_hub
           sp = hub%species_number(hat)
           if (h_species(sp)%cdft_charge_donor) then
             pub_cdft_group_charge_donor_u = hub%cdft_u_charge_up(hat)
             EXIT LOOP_atom_2_1
           endif
         enddo LOOP_atom_2_1

       elseif (pub_cdft_group_charge_down_only) then

         LOOP_atom_2_2: do hat = 1, mdl%par%nat_hub
           sp = hub%species_number(hat)
           if (h_species(sp)%cdft_charge_donor) then
             !pub_cdft_group_charge_donor_u = hub%cdft_u_charge_up(hat)
             pub_cdft_group_charge_donor_u = hub%cdft_u_charge_down(hat) ! UP-electrons are neglected...
             EXIT LOOP_atom_2_2
           endif
         enddo LOOP_atom_2_2

       endif
! modified group_charge_UP/DOWN_only
! gibo-end-02.10.12
    endif GROUP_CHARGE_DONOR

    GROUP_SPIN_ACCEPTOR: if (pub_cdft_group_spin_acceptor) then
      LOOP_atom_3: do hat = 1, mdl%par%nat_hub
        sp = hub%species_number(hat)
        if (h_species(sp)%cdft_spin_acceptor) then
          pub_cdft_group_spin_acceptor_u = hub%cdft_u_spin(hat)
          EXIT LOOP_atom_3
        endif
      enddo LOOP_atom_3
    endif GROUP_SPIN_ACCEPTOR

    GROUP_SPIN_DONOR: if (pub_cdft_group_spin_donor) then
      LOOP_atom_4: do hat = 1, mdl%par%nat_hub
        sp = hub%species_number(hat)
        if (h_species(sp)%cdft_spin_donor) then
          pub_cdft_group_spin_donor_u = hub%cdft_u_spin(hat)
          EXIT LOOP_atom_4
        endif
      enddo LOOP_atom_4
    endif GROUP_SPIN_DONOR

    GROUP_CHARGE_DIFF: if (pub_cdft_group_charge_diff) then
          ! take just the first U-potential (they are all equal..)
          ! modified for charge_spin_up/down_only

          atoms_cdft_10: do hat = 1, mdl%par%nat_hub

             sp = hub%species_number(hat)
             ! cycle loop if present species is not cDFT-active
             if (h_species(sp)%cdft_active) then

               if (.not.pub_cdft_group_charge_down_only) then
                pub_cdft_group_charge_diff_u = ABS(hub%cdft_u_charge_up(hat))
               elseif (pub_cdft_group_charge_down_only) then
                pub_cdft_group_charge_diff_u = ABS(hub%cdft_u_charge_down(hat)) ! DOWN-electrons are neglected
               endif

               EXIT atoms_cdft_10

             endif

          enddo atoms_cdft_10

    endif GROUP_CHARGE_DIFF

    GROUP_SPIN_DIFF: if (pub_cdft_group_spin_diff) then

          atoms_cdft_11: do hat = 1, mdl%par%nat_hub

             sp = hub%species_number(hat)
             ! cycle loop if present species is not cDFT-active
             if (h_species(sp)%cdft_active) then

               ! take just the first U-potential (they are all equal..)
               pub_cdft_group_spin_diff_u = ABS(hub%cdft_u_spin(hat))

               EXIT atoms_cdft_11

             endif

          enddo atoms_cdft_11

    endif GROUP_SPIN_DIFF

    ! broadcat group cDFT U-potential
    if (pub_cdft_group_charge_acceptor) &
           call comms_bcast(pub_root_proc_id,pub_cdft_group_charge_acceptor_u)

    if (pub_cdft_group_charge_donor) &
           call comms_bcast(pub_root_proc_id,pub_cdft_group_charge_donor_u)

    if (pub_cdft_group_spin_acceptor) &
           call comms_bcast(pub_root_proc_id,pub_cdft_group_spin_acceptor_u)

    if (pub_cdft_group_spin_donor) &
           call comms_bcast(pub_root_proc_id,pub_cdft_group_spin_donor_u)

    if (pub_cdft_group_charge_diff) &
           call comms_bcast(pub_root_proc_id,pub_cdft_group_charge_diff_u)

    if (pub_cdft_group_spin_diff) &
           call comms_bcast(pub_root_proc_id,pub_cdft_group_spin_diff_u)


      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving cdft_restart_read'

    end subroutine cdft_intermediate_restart_read

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

    subroutine dft_nu_intermediate_restart_read(hub,mdl)

    !==========================================================================!
    ! This subroutine reads the DFT-nu U1/U2 values from (.dft_nu) file and    !
    ! updates the hub%dft_nu_u1 etc. matrices                                  !
    !--------------------------------------------------------------------------!
    ! Written by Gilberto Teobaldi in December 2011                            !
    ! Ammended for DFT+nu by Glenn Moynihan in 2015                            !
    !==========================================================================!

      use constants,       only: DP, stdout, HARTREE_IN_EVS
      use comms,           only: pub_on_root, pub_root_proc_id, comms_bcast
      use hubbard_build,   only: HUBBARD_MODEL
      use hubbard_init, only: h_species
      use model_type, only: MODEL
      use rundat, only: pub_rootname, &
        pub_cdft_atom_charge, pub_cdft_atom_spin, &
        pub_cdft_group_charge_acceptor, pub_cdft_group_charge_donor, &
        pub_cdft_group_spin_acceptor, pub_cdft_group_spin_donor, &
        pub_cdft_group_charge_diff, pub_cdft_group_spin_diff, &
        pub_cdft_group_charge_acceptor_u, &
        pub_cdft_group_charge_donor_u,&
        pub_cdft_group_spin_acceptor_u, pub_cdft_group_spin_donor_u, &
        pub_cdft_group_charge_diff_u, pub_cdft_group_spin_diff_u,&
        pub_cdft_modes, pub_cdft_continuation, pub_task,&
        pub_cdft_group_charge_up_only,pub_cdft_group_charge_down_only,&
        pub_debug_on_root, pub_dft_nu_continuation, &
        pub_hubbard_unify_sites, pub_cdft_multi_proj
      use utils,           only: utils_unit, utils_open_unit_check, &
                                 utils_close_unit_check, utils_abort

      implicit none

      ! Arguments
      type(HUBBARD_MODEL), intent(inout) :: hub
      type(MODEL), intent(in) :: mdl

      ! Local variables
      integer           :: io_unit,io_stat
      integer           :: species, hat, sp, dummy_int
      character(len=80) :: filename
      ! ars: dummies
      character(len=4)  :: dummy_id
      integer           :: dummy_l
      real(kind=dp)     :: dummy_c
      real(kind=dp)     :: dummy_u_s, dummy_s_target
      real(kind=dp)     :: dummy_u_q_up, dummy_u_q_down
      real(kind=dp)     :: dummy_n_up, dummy_n_down
      real(kind=dp)     :: dummy_nu_up, dummy_nu_down
      real(kind=dp)     :: dummy_u1_up,dummy_u1_down
      real(kind=dp)     :: dummy_u2_up,dummy_u2_down
      integer           :: group

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering dft_nu_intermediate_restart_read'

      if (pub_on_root) then

         ! Find available unit specifier
         io_unit = utils_unit()
         write(filename,'(a,a)') trim(pub_rootname),'.dft_nu'

         open(unit=io_unit,iostat=io_stat,file=filename,&
              access='SEQUENTIAL',form='FORMATTED',position='REWIND', &
              action='READ')


         RESTART: if (io_stat == 0) then

           read(io_unit,'(/,/,/)')  !skip the 4-line header of .dft_nu

             if (.not. pub_hubbard_unify_sites) then

                 do species = 1, mdl%par%num_hub_species

                 read(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') dummy_int, &
                   & dummy_id,         &
                   & dummy_l,          &
                   & dummy_c,          &
                   & dummy_u1_up,      &
                   & dummy_u1_down,    &
                   & dummy_u2_up,      &
                   & dummy_u2_down,    &
                   & dummy_n_up,       &
                   & dummy_n_down

                    ! This is a restart so all figures are fixed
                    dummy_nu_up = dummy_n_up - floor(dummy_n_up)
                    dummy_n_up = floor(dummy_n_up)
                    dummy_nu_down = dummy_n_down - floor(dummy_n_down)
                    dummy_n_down = floor(dummy_n_down)

                   do hat = 1 , mdl%par%nat_hub

                    if (hub%species_number(hat) == species) then

                           ! internally convert U-potentials in Hartree
                           hub%nu_u1_up(hat)    = dummy_u1_up/&
                                                  HARTREE_IN_EVS
                           hub%nu_u1_down(hat)  = dummy_u1_down/&
                                                  HARTREE_IN_EVS
                           hub%nu_u2_up(hat)    = dummy_u2_up/&
                                                  HARTREE_IN_EVS
                           hub%nu_u2_down(hat)  = dummy_u2_down/&
                                                  HARTREE_IN_EVS

                           hub%nu_target_n_up(hat)  = dummy_n_up
                           hub%nu_target_n_down(hat)= dummy_n_down
                           hub%nu_nu_up(hat)        = dummy_nu_up
                           hub%nu_nu_down(hat)      = dummy_nu_down

                           hub%nu_group(hat)        = dummy_int
                    endif

                   enddo
                 enddo

             else if (pub_hubbard_unify_sites) then

                 hub%nu_projectors(:)=0

                 do species = 1, mdl%par%num_hub_species
                   do group = 1 , hub%num_groups

                    if (h_species(species)%nu_group == group) then

                        read(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') dummy_int, &
                          & dummy_id,         &
                          & dummy_l,          &
                          & dummy_c,          &
                          & dummy_u1_up,      &
                          & dummy_u1_down,    &
                          & dummy_u2_up,      &
                          & dummy_u2_down,    &
                          & dummy_n_up,       &
                          & dummy_n_down

                           ! This is a restart so all figures are fixed
                           dummy_nu_up = dummy_n_up - floor(dummy_n_up)
                           dummy_n_up = floor(dummy_n_up)
                           dummy_nu_down = dummy_n_down -&
                                            floor(dummy_n_down)
                           dummy_n_down = floor(dummy_n_down)

                           ! internally convert U-potentials in Hartree
                           hub%nu_u1_up(group)    = dummy_u1_up/&
                                                    HARTREE_IN_EVS
                           hub%nu_u1_down(group)  = dummy_u1_down/&
                                                    HARTREE_IN_EVS
                           hub%nu_u2_up(group)    = dummy_u2_up/&
                                                    HARTREE_IN_EVS
                           hub%nu_u2_down(group)  = dummy_u2_down/&
                                                    HARTREE_IN_EVS

                           hub%nu_target_n_up(group)   = dummy_n_up
                           hub%nu_target_n_down(group) = dummy_n_down
                           hub%nu_nu_up(group)         = dummy_nu_up
                           hub%nu_nu_down(group)       = dummy_nu_down
                           hub%nu_group(group)         = dummy_int

                          if (pub_cdft_multi_proj) then
                              hub%nu_projectors(group)   = &
                                hub%nu_projectors(group) + &
                                h_species(species)%cdft_num_proj
                          else
                              hub%nu_projectors(group)   = &
                                hub%nu_projectors(group) + &
                                h_species(species)%hub_ang_mom*2+ 1
                          endif

                      endif

                   enddo
                 enddo

             end if

          close(unit=io_unit,iostat=io_stat)
          call utils_close_unit_check('dft_nu_restart_read',filename,io_stat)


         else
           write(stdout,'(/,3a)') 'ERROR dft_nu_restart_read could &
             &not be opened "', trim(filename),'"'
           write(stdout,'(/a)') '...aborting...'

! gibo: this will be needed when we talk about geometry optimizations...

!           if (pub_dft_nu_continuation) then
!             write(stdout,'(/,a)') 'ERROR: it looks like you have set&
!               &"cdft_continuation=.TRUE." in your input file withouth having &
!               &the .cdft file in your work-directory'
!             write(stdout,'(/,a)') 'SOLUTION: comment out &
!               &"cdft_continuation=.TRUE." in your input file'
!
!           elseif ((pub_task=='GEOMETRYOPTIMIZATION')  .OR. &
!                   (pub_task=='MOLECULARDYNAMICS')     .OR. &
!                   (pub_task=='TRANSITIONSTATESEARCH') .OR. &
!                   (pub_task=='PROPERTIES')            .OR. &
!                   (pub_task=='COND')                  .OR. &
!                   (pub_task=='PROPERTIES_COND') ) then
!             call utils_abort('A .cdft file must be present &
!               &to perform a constrained-DFT GEOMETRYOPTIMIZATION / &
!               &MOLECULARDYNAMICS / TRANSITIONSTATESEARCH / PROPERTIES / &
!               &COND / PROPERTIES_COND task. SOLUTION: set &
!               &"task : SINGLEPOINT" in your input file to obtain an optimised
!               &
!               & .cdft file first.')
!           end if

         endif RESTART

      endif

      !broadcast updated U-potentials
      call comms_bcast(pub_root_proc_id, hub%nu_u1_up(:))
      call comms_bcast(pub_root_proc_id, hub%nu_u1_down(:))
      call comms_bcast(pub_root_proc_id, hub%nu_u2_up(:))
      call comms_bcast(pub_root_proc_id, hub%nu_u2_down(:))

     if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving dft_nu_intermediate_restart_read'

    end subroutine dft_nu_intermediate_restart_read
!

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

    subroutine dft_nu_intermediate_update_matrices(hub, hub_proj_basis)

    !==========================================================================!
    ! This subroutine updates the DFT_nu SPAM3 matrices                        !
    !--------------------------------------------------------------------------!
    ! Written for DFT+nu module by Glenn Moynihan in 2015 and based on         !
    ! hubbard_build_matrices by Gilberto Teobaldi in Dec. 2011                 !
    !==========================================================================!
    ! Potentially useful comment:                                              !
    ! for cDFT run, we neglect hub_spin_splitting and hub_alpha terms.         !
    ! This allows us to safely update directly hub%down_matrix and             !
    ! hub%down_matrix with the udpated cDFT U-potentials, as we do below...    !
    ! Modified to remove pub_par by Robert Charlton, 25/10/2018.               !
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use hubbard_build, only: HUBBARD_MODEL
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_cdft, &
        pub_cdft_atom_charge, pub_cdft_atom_spin, &
        pub_cdft_group_charge_acceptor, pub_cdft_group_charge_donor, &
        pub_cdft_group_spin_acceptor, pub_cdft_group_spin_donor, &
        pub_cdft_group_charge_diff, pub_cdft_group_spin_diff, &
        pub_cdft_modes, pub_debug_on_root, pub_hubbard_unify_sites
    use sparse, only: SPAM3, sparse_put_element, &
         sparse_get_element, sparse_get_par


    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    integer :: sp, hat_on_proc, hat, theatom
    integer :: hub_proj
    !integer :: is, ierr
    real(kind=DP) :: u_charge_up, u_charge_down, u_spin_up, u_spin_down
    real(kind=DP) :: half_u2_up, half_u2_down
    real(kind=DP) :: potential_up, potential_down

    integer :: group
    type(PARAL_INFO), pointer :: par

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &dft_nu_intermediate_update_matrices'

    ! rc2013: get the parallel strategy for this SPAM3
    call sparse_get_par(par, hub%u_matrix_up)

!----
    if (.not.pub_hubbard_unify_sites) then
       ! ddor: Loop over Hubbard atoms on my proc
        do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

           hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
           sp = hub%species_number(hat)
           theatom = par%distr_atom(hub%orig(hat))


                 half_u2_up      = hub%nu_u2_up(hat) * (-1.0_DP)
                 half_u2_down    = hub%nu_u2_down(hat) * (-1.0_DP)

                 ! gom: V_up = U1/2
                 potential_up    = hub%nu_u1_up(hat) * 0.5_DP
                 potential_down  = hub%nu_u1_down(hat) * 0.5_DP

                 ! gom: Loop over Hubbard projectors on my proc
                 do hub_proj = hub_proj_basis%first_on_atom(theatom), &
                      hub_proj_basis%first_on_atom(theatom) + &
                      hub_proj_basis%num_on_atom(theatom) - 1

                    call sparse_put_element(half_u2_up,&
                         hub%u_matrix_up,hub_proj,hub_proj)
                    call sparse_put_element(half_u2_down,&
                         hub%u_matrix_down,hub_proj,hub_proj)

                    call sparse_put_element(potential_up,&
                         hub%up_matrix,hub_proj,hub_proj)
                    call sparse_put_element(potential_down,&
                         hub%down_matrix,hub_proj,hub_proj)

                 enddo


        enddo !hat_on_proc

    else if (pub_hubbard_unify_sites) then

       ! grom: Loop over groups
        do group = 1, hub%num_groups

            half_u2_up      = hub%nu_u2_up(group) * (-1.0_DP)
            half_u2_down    = hub%nu_u2_down(group) * (-1.0_DP)

            ! gom: V_up = U1/2
            potential_up    = hub%nu_u1_up(group) * 0.5_DP
            potential_down  = hub%nu_u1_down(group) * 0.5_DP

            ! gom: Loop over Hubbard projectors in my group
            do hub_proj = 1,hub%nu_projectors(group)

               call sparse_put_element(half_u2_up,&
                    hub%u_matrix_up,hub_proj,hub_proj)
               call sparse_put_element(half_u2_down,&
                    hub%u_matrix_down,hub_proj,hub_proj)

               call sparse_put_element(potential_up,&
                    hub%up_matrix,hub_proj,hub_proj)
               call sparse_put_element(potential_down,&
                    hub%down_matrix,hub_proj,hub_proj)

            enddo

        enddo !groups

    end if

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &dft_nu_intermediate_update_matrices'

    end subroutine dft_nu_intermediate_update_matrices

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------



end module cdft_intermediate_cg
