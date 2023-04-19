! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                   E N E R G Y      D E C O M P O S I T I O N                !
!                               A N A L Y S I S                               !
!                          __                                                 !
!                         / ()        _  ,_   _,                              !
!                         >-   /|/|  |/ /  | / | |  |                         !
!                         \___/ | |_/|_/   |/\/|/ \/|/                        !
!                                             (|   (|                         !
!             ___                                                             !
!            (|  \  _  _   _               _   ,  o_|_ o  _                   !
!             |   ||/ /   / \_/|/|/|  |/\_/ \_/ \_| |  | / \_/|/|             !
!            (\__/ |_/\__/\_/  | | |_/|_/ \_/  \/ |/|_/|/\_/  | |_/           !
!                                    (|                                       !
!                       __,                                                   !
!                      /  |         _,  |\       ,  o  ,                      !
!                     |   |  /|/|  / |  |/ |  | / \_| / \_                    !
!                      \_/\_/ | |_/\/|_/|_/ \/|/ \/ |/ \/                     !
!                                            (|                               !
!                                                                             !
!=============================================================================!
!                                                                             !
! This module is the driver for energy decomposition analysis (EDA) within    !
! ONETEP.                                                                     !
!                                                                             !
! Notes:                                                                      !
! -The module stores fragment data in a public array of pub_frag_data objects.!
! -Fragment energies are calculated by calling energy_and_force_calculate().  !
! -A second driver, eda_driver_supermol_run(), handles the                    !
!   optimisations of the supermolecule states. This is called within          !
!   energy_and_force_calculate(). In this sense,                              !
!   energy_and_force_calculate() is used as an initialiser for the            !
!   fragment and supermolecule states' NGWFs, density kernels and other       !
!   quantities.                                                               !
! -Polarisation is calculated using the ALMO approach.  This involves using   !
!   the Stoll SCF-MI equations to effectively localise the MOs to specific    !
!   pre-determined fragments.                                                 !
!                                                                             !
! @TODO:                                                                      !
! - A number of XC functionals are not currently compatible with EDA.         !
! - The EDA module has been modified to not work with classical point         !
!   charges. The number of classical atoms in a fragment is overriden to zero.!
! Known bugs:                                                                 !
! - If pub_eda_read_super=.true. and pub_eda_reset_ngwfs_pol=.true. then      !
!   NGWFs appear not to always be reset.                                      !
!                                                                             !
!-----------------------------------------------------------------------------!
! Originally written by Max Phipps July 2013                                  !
!-----------------------------------------------------------------------------!

module eda_main


  implicit none

  private

  ! subroutines
  public :: eda_task         ! Main EDA driver (for both the EDA and
                             ! supermolecule preparation tool)
!  public :: eda_prep_frags  ! For preparing fragment data objects
!  public :: eda_prep_super  ! For preparing supermolecule data object

contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  subroutine eda_task(full_mdl)

    !======================================================================!
    ! This suboutine is the main driver for the energy decomposition       !
    ! analysis in ONETEP.  This calculates the singlepoint energies of     !
    ! the fragments in isolation, the 'frozen' state constructed by        !
    ! combining the densities of the optimised fragments, the polarised    !
    ! state by restricted variational optimisation, and the fully relaxed  !
    ! density.  These energies are used to calculate the electrostatics,   !
    ! exchange, correlation, Pauli-repulsion, polarisation and charge      !
    ! transfer contributions to the interaction energy (as defined in the  !
    ! ALMO EDA [1] and LMO EDA [2] scheme literature).                     !
    !                                                                      !
    ! If the task is set to 'EDA_PREP' then the frozen supermolecule       !
    ! state is written to disk for later calculation. This relieves the    !
    ! maximum processor limit of supermolecule-stage EDA calculations from !
    ! running on up to the minimum of the number of fragment atoms to      !
    ! running on up to the number of supermolecule atoms.                  !
    !                                                                      !
    ! [1] R. Z. Khaliullin, R. Lochan, E. Cobar, A. T. Bell, and M.        !
    ! Head-Gordon, J. Phys. Chem. A 111, 8753 (2007)                       !
    ! [2] P.Su and H.Li, J. Chem. Phys. 131, 014102/1-15(2009)             !
    !                                                                      !
    !======================================================================!
    ! Theory Overview                                                      !
    ! The ONETEP EDA is a hybrid ALMOEDA/LMOEDA scheme that includes an    !
    ! expansion of the ALMO EDA scheme frozen density component, FRZ, in   !
    ! its separated electrostatic, exchange, correlation, and Pauli        !
    ! repulsion parts.                                                     !
    !                                                                      !
    ! Briefly,                                                             !
    !   dE_ALMO = FRZ + POL + CT                                           !
    !   dE_LMO  = ES + EX + REP + POL' + DISP                              !
    ! where FRZ is the frozen density interaction, POL is the ALMOEDA      !
    ! SCF-MI polarisation energy, and CT is the (remainder) charge transfer!
    ! term, and where ES, EX, REP are the electrostatics, exchange, and    !
    ! Pauli repulsions contributions, and POL' the LMOEDA "polarisation"   !
    ! contribution that also includes delocalisation effects between       !
    ! fragments, and DISP is the difference in correlation energy of the   !
    ! isolated fragments and the fully relaxed supermolecule system.       !
    !                                                                      !
    ! Redefining DISP = CORR + CORR_mix, where CORR is the correlation     !
    ! energy change on forming the frozen density state, and CORR_mix is   !
    ! the correlation contribution to the ALMOEDA POL and CT energies then,!
    !   dE_ONETEP = ES + EX + REP + CORR + POL + CT                        !
    !                                                                      !
    !----------------------------------------------------------------------!
    ! Written by Max Phipps, August 2013-2015                              !
    ! Modified by Max Phipps December 2015 to improve code structure.      !
    !======================================================================!

    use constants, only: EDA_ISOLATED, EDA_FROZENNONIDEM, &
         DP, stdout
    use comms, only: pub_on_root, comms_barrier
    use energy_and_force, only: energy_and_force_calculate, &
        energy_and_force_init_cell, energy_and_force_exit_cell
    use fragment_data, only: pub_frag_data, fragment_data_destroy, &
         complex_energy_frz, complex_energy_frzIdem, &
         complex_energy_pol_simul, complex_energy_pol, &
         complex_energy_ct, complex_energy_full
    use model_type, only: MODEL
    use parallel_strategy, only: PARAL_INFO, pub_par
    use rundat, only: pub_debug_on_root, pub_any_nl_proj, &
         pub_eda_mode, pub_frag_counter, pub_frag_charge, &
         pub_eda_split_atoms, pub_frag_atoms, pub_frag_deloc, &
         pub_charge, pub_eda_read_super, &
         pub_eda_intrafrag_x_frz, pub_eda_intrafrag_c_frz, &
         pub_eda_dE_c_rep, pub_frag_iatm, &
         pub_eda_isol_x, pub_eda_isol_c, pub_eda_isol, &
         pub_eda_have_sp_frags, pub_num_spins, pub_spin_polarised, &
         pub_eda_frag_isol_pol, &
         pub_eda_frag_isol_ct, pub_eda_preptool, &
         pub_eda_scfmi, pub_eda_scfmi_any, pub_print_qc
         !, pub_eda_scfmi_restrict_ngwfs
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_banner, utils_qc_print


    implicit none


    ! Arguments
    type(MODEL), target, intent(inout)  :: full_mdl

    ! Local Variables
    real(kind=DP)     :: total_forces(1:3,1:full_mdl%nat)
    real(kind=DP)     :: dummy_energy
    integer           :: it, frag1, frag2, ierr
    integer           :: frag_index, nsubfrags
    character(len=3)  :: qc_intstr

    ! if the user forced spin polarisation, or if ONETEP has determined
    ! that spin polarisation is neccesary:
    logical :: user_forced_spin_polarised

    type(PARAL_INFO), pointer :: par


    ! Warnings/checks:
    call internal_warnings(full_mdl)

    ! jcap: initialise par here
    call internal_initialise_par(full_mdl)

    ! mjsp: Initialisations
    user_forced_spin_polarised = pub_spin_polarised

    ! mjsp: Calculation involves SCF-MI at some point:
    pub_eda_scfmi_any = .true.


    if(pub_debug_on_root) then
       write(stdout,'(a)') 'DEBUG info: fragment atom definitions'
       write(stdout,'(a)') '-------------------------------------'
       write(stdout,'(1x,a13,1x,a10)') 'Fragment ID', 'Atom Count'
       do it=1,pub_frag_iatm(0)
          write(stdout,'(8x,I3,3x,I5)') it,pub_frag_iatm(it)
       end do
       write(stdout,'(a)') '-------------------------------------'
    end if


    ! mjsp: allocate energy storage arrays:
    allocate(complex_energy_pol(0:pub_frag_iatm(0)), stat=ierr)
    call utils_alloc_check('eda_task','complex_energy_pol',ierr)
    allocate(complex_energy_ct(0:pub_frag_iatm(0),0:pub_frag_iatm(0)), stat=ierr)
    call utils_alloc_check('eda_task','complex_energy_ct',ierr)

    ! mjsp: isolated fragment energies.
    ! mjsp: Note: idx 0 stores total of fragment energies.
    allocate(pub_eda_isol(0:pub_frag_iatm(0)), stat=ierr)
    call utils_alloc_check('eda_task','pub_eda_isol',ierr)

    ! mjsp: initialise the storage of the sum of the individual fragment
    ! mjsp: energies in index 0
    pub_eda_isol(0) = 0.0_DP

    ! mjsp: isolated fragment x,c energies.
    allocate(pub_eda_isol_x(1:pub_frag_iatm(0)), stat=ierr)
    call utils_alloc_check('eda_task','pub_eda_isol_x',ierr)
    allocate(pub_eda_isol_c(1:pub_frag_iatm(0)), stat=ierr)
    call utils_alloc_check('eda_task','pub_eda_isol_c',ierr)

    ! mjsp: Supermolecule initialisation:
    call eda_prep_super(full_mdl)

    ! =========================== TEST SUPERMOLECULE SYSTEM ================

    ! mjsp: Initialise and exit the system cells to check padded cell etc.
    ! for issues
    call internal_check_supermolecule(pub_frag_data(0)%mdl, &
         user_forced_spin_polarised)

    ! ========================= END TEST SUPERMOLECULE SYSTEM ==============


    ! mjsp: if all-in-one calculation that includes fragment calculations
    if (.not.(pub_eda_read_super)) then

       ! mjsp: Fragments initialisation:
       call eda_prep_frags(full_mdl)

       call comms_barrier


       ! ========================= TEST MONOMER SYSTEMS =====================
       ! mjsp: Initialise and exit the system cells to check padded cell etc.
       ! for issues (NOTE: duplicated in energy_and_force_calculate())
       do it=1,pub_frag_iatm(0)
          call internal_check_fragment(it, pub_frag_data(it)%mdl, &
               user_forced_spin_polarised)
       end do

       if(pub_on_root) then
          write(stdout,'(/a)') utils_banner('=')
          write(stdout,'(a)') &
               ' EDA: Finished initialisation checks for fragments'
          write(stdout,'(a/)') utils_banner('=')
       end if

       ! ========================== END TEST MONOMER SYSTEMS ================


       ! ========================== MONOMER ENERGY CALCULATIONS =============

       ! mjsp: Update flag to isolated fragment mode - this ensures
       ! pub_frag_data is updated with the fragment data at the end of each
       ! calculation.
       pub_eda_mode = EDA_ISOLATED

       ! mjsp: Calculate monomer energies
       do pub_frag_counter=1,pub_frag_iatm(0) ! idx 0 stores # atoms

          ! unless the user has requested otherwise, or if impossible to
          ! run the supermolecule calculation without spin polarisation, then
          ! assume spin un-polarised calculation is possible for the fragments:
          ! this assumption will be corrected in energy_and_force_init_cell
          ! if this is not possible.
          if (.not. user_forced_spin_polarised) pub_spin_polarised = .false.

          if(pub_on_root) then
             write(stdout,'(//a)') utils_banner('=')
             write(stdout,'(a,i3)') &
                  ' EDA: Monomer calculation for fragment #', pub_frag_counter

             ! Check if any intrafragmental splits
             if (pub_on_root .and. pub_eda_split_atoms) then
                if (ANY(pub_frag_atoms(1,:) .eq. pub_frag_counter)) then
                   do frag_index=1,size(pub_frag_atoms(1,:))
                      if (pub_frag_atoms(1,frag_index) .eq. pub_frag_counter) &
                         write(stdout,'(a,i6)') &
                              ' --> Intrafragment split at atom: ', &
                              pub_frag_atoms(2,frag_index)
                   end do
                end if
             end if

             write(stdout,'(a/)') utils_banner('=')
          end if

          ! direct public parallel strategy to fragment parallel strategy
          pub_par=>pub_frag_data(pub_frag_counter)%mdl%par

          ! mjsp: Update public charge
          pub_charge = pub_frag_charge(pub_frag_counter)

          ! mjsp: initialise the cell with our system
          call energy_and_force_init_cell(pub_frag_data(pub_frag_counter)%mdl)

          ! update fragment data storage with the updated pub_any_nl_proj logical
          pub_frag_data(pub_frag_counter)%any_nl_proj = pub_any_nl_proj

          ! update fragment data storage with the newly calculated spin
          pub_frag_data(pub_frag_counter)%num_spins = pub_num_spins

          pub_eda_scfmi = .false.
          ! mjsp: If intrafragment partitioning involved:
          if (pub_eda_split_atoms) then
             if (ANY(pub_frag_atoms(1,:).eq.pub_frag_counter)) then

                ! mjsp: Number of subfragments:
                nsubfrags = count(pub_frag_atoms(1,:).eq.pub_frag_counter)

!                pub_eda_scfmi = .true.

             end if

          end if


          if(pub_debug_on_root) then
             write(stdout,'(a)') &
                  'DEBUG info: pub_cell for EDA monomer calculation'
             write(stdout,*) "DEBUG info: num_spec:      ", &
                  pub_frag_data(pub_frag_counter)%mdl%par%num_species
             write(stdout,*) "DEBUG info: num_pspec:     ", &
                  pub_frag_data(pub_frag_counter)%mdl%par%num_pspecies
             write(stdout,*) "DEBUG info: num_ngwfs:     ", &
                  pub_frag_data(pub_frag_counter)%mdl%par%num_ngwfs
             write(stdout,*) "DEBUG info: num_spins:     ", &
                  pub_frag_data(pub_frag_counter)%num_spins
             write(stdout,*) "DEBUG info: nat:           ", &
                  pub_frag_data(pub_frag_counter)%mdl%nat
             write(stdout,*) "DEBUG info: nat_classical: ", &
                  pub_frag_data(pub_frag_counter)%mdl%nat_classical
          end if


          ! mjsp: Calculate the energy now that we have finished initialising
          call comms_barrier
          call energy_and_force_calculate(pub_eda_isol(pub_frag_counter), &
               total_forces,pub_frag_data(pub_frag_counter)%mdl)
          call comms_barrier

          ! mjsp: store the sum of the individual fragment energies in index 0
          pub_eda_isol(0) = pub_eda_isol(0) + pub_eda_isol(pub_frag_counter)

          ! cleanup for next monomer
          call energy_and_force_exit_cell(pub_frag_data(pub_frag_counter)%mdl)

       end do


       ! ====================== END MONOMER ENERGY CALCULATIONS ================

       ! ================ INFORMATION PRINTOUT ============
       if (pub_on_root) then
          write(stdout,'(//a)') utils_banner('=')
          write(stdout,'(a)') utils_banner('=', &
               'Energy Decomposition Analysis (EDA)')
          write(stdout,'(a)') utils_banner('=')
          write(stdout,'(a)') '==> Monomer Energy Calculations Finished '
          write(stdout,'(1x,a,9x,a)') &
               'Fragment ID     Energy (Ha)        Charge','Spins'
          do it=1,pub_frag_iatm(0)
             write(stdout,'(4x,I3,3x,f22.14,5x,f4.1,12x,I1)') &
                  it, pub_eda_isol(it), &
                  pub_frag_charge(it), pub_frag_data(it)%num_spins
          end do
          write(stdout,'(a,2x,f22.14)') '   Total', pub_eda_isol(0)  !, pub_frag_charge(0)
          if (pub_spin_polarised) write(stdout,'(a)') '    ==> Supermolecule &
              &calculation will be ran spin-polarised.'
          write(stdout,'(a//)') utils_banner('=')
       end if
       ! ================ END INFORMATION PRINTOUT ============

       ! mjsp: update flag for if we have any spin polarised fragments
       pub_eda_have_sp_frags = .false.
       do it=1,pub_frag_iatm(0) ! loop fragments
          if (pub_frag_data(it)%num_spins == 2) then

            ! update the flag
            pub_eda_have_sp_frags = .true.

            ! force spin polarised calculation
            pub_spin_polarised = .true.

            exit

          end if
       end do

    end if  ! not(pub_eda_read_super)


    ! ============== FULL (COMPLEX) SYSTEM ENERGY CALCULATIONS =================
    ! mjsp: We have finished the isolated fragment calculations.  Continue with
    ! the supermolecule calculations (frozen, frozen idempotent, polarised, fully
    ! optimized).

    ! mjsp: A reconstruction of the density kernel from the fragment submatrices
    ! is performed in energy_and_force_calculate if the flag is in this state
    pub_eda_mode = EDA_FROZENNONIDEM
    pub_frag_counter = 0  ! supermolecule

    ! mjsp: point par to the supermolecule model's parallel strategy
    pub_par=>pub_frag_data(0)%mdl%par

    ! mjsp: Update public charge
    pub_charge = pub_frag_charge(0)

    ! mjsp: Modify cell to the full system
    call energy_and_force_init_cell(pub_frag_data(0)%mdl)

    ! mjsp: Ensure global flag controlling initialisation of density kernel from
    ! the core Hamiltonian is set to false (as we load these from the internally
    ! stored fragment NGWFs)
    !coreham_denskern_guess = .false.

    ! mjsp: energy calculation using pub_cell and frag_mdl%elements
    ! mjsp: parameters we have now set
    call comms_barrier
    call energy_and_force_calculate(dummy_energy,total_forces,full_mdl)
    call comms_barrier

    ! ========== END FULL (COMPLEX) SYSTEM ENERGY CALCULATIONS ==============

    ! ================ INFORMATION PRINTOUT ============
    if (pub_on_root) then

       write(stdout,'(//a)') utils_banner('=')
       write(stdout,'(a)') &
            utils_banner('=', 'Energy Decomposition Analysis (EDA)')
       write(stdout,'(a)') utils_banner('=')

       if (pub_eda_preptool) then

          ! TASK : EDA_PREP
          ! This calculation task is for preparing supermolecule-stage data.

          write(stdout,'(/2x,a/)') ' Supermolecule data written to disk &
               &for supermolecule-stage calculation.'

       else

          ! Full EDA output:

          write(stdout,'(/a)') &
               '==> (Fully) Variationally Optimized Supermolecule Finished '
          write(stdout,'(/a)') utils_banner('=')
          write(stdout,'(/a/)') '==> Final EDA Calculation Results:'

          write(stdout,'(2x,a,f10.3,a/)') &
               '(1+2) Frozen Density Component E_FRZ           : ', &
               (complex_energy_frzIdem - pub_eda_isol(0)) *627.509438736, &
               ' kcal/mol'

          write(stdout,'(6x,a,f10.3,a)') &
               '(1) (Nonortho.) Frozen Interaction Energy  : ', &
               (complex_energy_frz - pub_eda_isol(0))*627.509438736, ' kcal/mol'
          write(stdout,'(15x,a,f10.3,a)') &
               'Frozen Electrostatics E_ES        : ', &
               (complex_energy_frz - pub_eda_isol(0) - pub_eda_intrafrag_x_frz &
               - pub_eda_intrafrag_c_frz ) *627.509438736, ' kcal/mol'
          write(stdout,'(15x,a,f10.3,a)') &
               'Frozen Exchange E_EX              : ', &
               (pub_eda_intrafrag_x_frz) *627.509438736, ' kcal/mol'
          write(stdout,'(15x,a,f10.3,a/)') &
               'Frozen Correlation E_FRZ-CORR     : ', &
               (pub_eda_intrafrag_c_frz) *627.509438736, ' kcal/mol'

          write(stdout,'(6x,a,f10.3,a)') &
               '(2) Repulsion Interaction Energy           : ', &
               (complex_energy_frzIdem - complex_energy_frz)*627.509438736, &
               ' kcal/mol'
          write(stdout,'(15x,a,f10.3,a)') &
               'Pauli-Repulsion E_REP             : ', &
               (complex_energy_frzIdem - complex_energy_frz &
               - pub_eda_dE_c_rep)*627.509438736, ' kcal/mol'
          write(stdout,'(15x,a,f10.3,a/)') &
               'Repulsion Correlation E_REP-CORR  : ', &
               (pub_eda_dE_c_rep)*627.509438736, ' kcal/mol'


          write(stdout,'(/2x,a,f10.3,a)') &
               '(3) Delocalisation Interaction Energy          : ', &
               (complex_energy_full - complex_energy_frzIdem)*627.509438736, &
               ' kcal/mol'


          write(stdout,'(/6x,a,f10.3,a)') &
               '(Stoll SCF-MI) Polarisation E_POL          : ', &
               (complex_energy_pol_simul - complex_energy_frzIdem) &
               *627.509438736, ' kcal/mol'

          ! individual fragment polarisations
          if (pub_eda_frag_isol_pol) then

             write(stdout,'(6x,a70)') '--------------&
                  &P-O-L-A-R-I-S-A-T-I-O-N----A-N-A-L-Y-S-I-S----------------'

             write(stdout,'(51x,a11)') 'Fragment ID'
             do it=1,pub_frag_iatm(0)
               write(stdout,'(19x,a22,i5,f10.3,a)') &
                    'E_POL energy :', it, &
                    (complex_energy_pol(it) &
                    - complex_energy_frzIdem)*627.509438736, ' kcal/mol'
             end do

             write(stdout,'(14x,a,5x,f10.3,a)') &
                  'Higher Order E_POL energy :', &
                  (complex_energy_pol_simul-complex_energy_pol(0)&
                  -complex_energy_frzIdem)*627.509438736, ' kcal/mol'

             write(stdout,'(6x,a70)') '--------------&
                  &-----------------------------------------------------------'

          end if


          write(stdout,'(/6x,a,f10.3,a)') &
               'Charge Transfer E_CT                       : ', &
               (complex_energy_full - complex_energy_pol_simul)*627.509438736, &
               ' kcal/mol'

          ! individual fragment pair delocalisations
          if (pub_eda_frag_isol_ct) then

             write(stdout,'(6x,a70)') &
                  '----------D-E-L-O-C-A-L-I-S-A-T-I-O-N----&
                  &A-N-A-L-Y-S-I-S----------------'

             write(stdout,'(10x,a60)') &
                  '  Fragment ID1  Fragment ID2  Delocn. Energy (kcal/mol)'

! User-specified fragment pairs:
             ! loop the fragment pairs given in the input file:
             do it=1,SIZE(pub_frag_deloc(1,:))
                frag1 = pub_frag_deloc(1,it)
                frag2 = pub_frag_deloc(2,it)
                write(stdout,'(20x,i5,8x,i5,8x,f10.3)') frag1, frag2, &
                     (complex_energy_ct(frag1,frag2) - &
                     complex_energy_pol_simul)*627.509438736
             end do


! All fragment pairs:
!             do it=1,pub_frag_iatm(0)
!               do it2=2,pub_frag_iatm(0)
!                 if (it.ge.it2) cycle
!                 write(stdout,'(20x,i5,8x,i5,8x,f10.3)') it, it2, &
!                      (complex_energy_ct(it,it2) - &
!                      complex_energy_pol_simul)*627.509438736
!               end do
!             end do
!
!             write(stdout,'(/6x,a50,f10.3,a)'), &
!                  'Higher Order Delocalisation (Stoll SCF-MI) :', &
!                  (complex_energy_full-complex_energy_ct(0,0)-&
!                  complex_energy_pol_simul)*627.509438736, ' kcal/mol'

             write(stdout,'(6x,a70)') '--------------&
                  &-----------------------------------------------------------'

          end if


          write(stdout,'(/2x,a,f10.3,a/)') &
               '(1+2+3) Total Interaction Energy               : ', &
                (complex_energy_full - pub_eda_isol(0))*627.509438736, &
               ' kcal/mol'

       end if

       write(stdout,'(a//)') utils_banner('=')

    end if  ! pub_on_root

    ! mjsp: QC tests
    if (pub_print_qc .and. pub_on_root) then
       write(stdout,'(a1)')' '
       do it=1,pub_frag_iatm(0)
          write(qc_intstr,'(I3)') it
          call utils_qc_print('EDA frag x '//qc_intstr, pub_eda_isol_x(it))
          call utils_qc_print('EDA frag c '//qc_intstr, pub_eda_isol_c(it))
          ! mjsp: If fragment polarisation QC test, then include fragment
          ! mjsp: polarisation energies:
          if (pub_eda_frag_isol_pol) then
            call utils_qc_print('EDA cmplx E pol '//qc_intstr, complex_energy_pol(it))
          end if
       end do
       call utils_qc_print('EDA intrfrg x frz',pub_eda_intrafrag_x_frz)
       call utils_qc_print('EDA intrfrg c frz',pub_eda_intrafrag_c_frz)
       call utils_qc_print('EDA dE c rep',pub_eda_dE_c_rep)
       call utils_qc_print('EDA cmplx E frz', complex_energy_frz)
       call utils_qc_print('EDA cmplx E frzIdem', complex_energy_frzIdem)
       call utils_qc_print('EDA cmplx E pol_simul', complex_energy_pol_simul)
       call utils_qc_print('EDA cmplx E full', complex_energy_full)
    endif


    ! ================ END INFORMATION PRINTOUT ============

    ! cleanup
    call energy_and_force_exit_cell(pub_frag_data(0)%mdl)

    ! mjsp: destroy (public) fragment container
    call fragment_data_destroy()

    ! mjsp: deallocate energy storage arrays:
    deallocate(complex_energy_pol, stat=ierr)
    call utils_dealloc_check('eda_task','complex_energy_pol',ierr)
    deallocate(complex_energy_ct, stat=ierr)
    call utils_dealloc_check('eda_task','complex_energy_ct',ierr)

    deallocate(pub_eda_isol, stat=ierr)
    call utils_dealloc_check('eda_task','pub_eda_isol',ierr)
    deallocate(pub_eda_isol_c, stat=ierr)
    call utils_dealloc_check('eda_task','pub_eda_isol_c',ierr)
    deallocate(pub_eda_isol_x, stat=ierr)
    call utils_dealloc_check('eda_task','pub_eda_isol_x',ierr)

  end subroutine eda_task

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eda_prep_super(full_mdl)

    !======================================================================!
    ! This suboutine prepares the supermolecule data object.               !
    !======================================================================!
    ! Written by Max Phipps, August 2013-2015                              !
    ! Modified for embedding by Joseph Prentice, September 2018            !
    !======================================================================!

    use utils, only: utils_alloc_check
    use fragment_data, only: pub_frag_data
    use parallel_strategy, only: pub_par
    use model_type, only: MODEL

    implicit none


    ! Arguments
    type(MODEL), target, intent(inout)  :: full_mdl

    ! Local Variables
    integer :: num_atoms, ierr


    ! mjsp: Supermolecule initialisation:
    pub_par=>pub_frag_data(0)%mdl%par
    ! jcap: need to copy everything from full model, including regions
    ! - allocate in each region

    ! number of atoms in fragment/supermolecule
    num_atoms = full_mdl%nat
    pub_frag_data(0)%mdl%nat = full_mdl%nat
    pub_frag_data(0)%mdl%nat_classical = 0
    pub_frag_data(0)%mdl%nsub = full_mdl%nsub

    ! copy species from the supermolecule model
    allocate(pub_frag_data(0)%mdl%species(1:size(full_mdl%species)),stat=ierr)
    call utils_alloc_check('eda_task','pub_frag_data(0)%mdl%species',ierr)
    pub_frag_data(0)%mdl%species = full_mdl%species
    pub_frag_data(0)%mdl%num_pspecies=full_mdl%num_pspecies

    ! allocate fragment storage for elements
    allocate(pub_frag_data(0)%mdl%elements(1:num_atoms),stat=ierr)
    call utils_alloc_check('eda_task','pub_frag_data(0)%mdl%elements',ierr)

    ! mjsp: Fill the fragment elements
    pub_frag_data(0)%mdl%elements(1:num_atoms) = full_mdl%elements

    ! jcap: allocate regions - only compatible with 1 region
    allocate(pub_frag_data(0)%mdl%regions(1),stat=ierr)
    call utils_alloc_check('eda_task','pub_frag_data(0)%mdl%regions',ierr)

    ! jcap: copy the appropriate information over
    ! copy full model into the supermolecule fragment storage
    call internal_copy_par(pub_frag_data(0)%mdl%regions(1)%par,&
         full_mdl%regions(1)%par)
    ! jcap: zero number of classical atoms
    pub_frag_data(0)%mdl%regions(1)%par%nat_classical = 0
    pub_frag_data(0)%mdl%regions(1)%reg_count = 1

    ! copy species from the supermolecule model
    allocate(pub_frag_data(0)%mdl%regions(1)%species(1:size(full_mdl%regions(1)%species)),&
         stat=ierr)
    call utils_alloc_check('eda_task','pub_frag_data(0)%mdl%regions%species',ierr)
    pub_frag_data(0)%mdl%regions(1)%species = full_mdl%regions(1)%species

    ! allocate fragment storage for elements
    allocate(pub_frag_data(0)%mdl%regions(1)%elements(1:size(full_mdl%regions(1)%elements)),&
         stat=ierr)
    call utils_alloc_check('eda_task','pub_frag_data(0)%mdl%regions%elements',ierr)

    ! mjsp: Fill the fragment elements
    pub_frag_data(0)%mdl%regions(1)%elements = full_mdl%regions(1)%elements

    ! jcap: point to the appropriate parallel strategy
    pub_frag_data(0)%mdl%par => pub_frag_data(0)%mdl%regions(1)%par
    pub_par => pub_frag_data(0)%mdl%regions(1)%par

    ! jcap: allocate ngwf_basis array
    allocate(pub_frag_data(0)%ngwf_basis(1),stat=ierr)
    call utils_alloc_check('eda_task','pub_frag_data(0)%ngwf_basis',ierr)

    ! jcap: distribute the ngwfs
    ! mjsp: for EDA calculations, the parallel strategies of the
    ! supermolecule and all the fragments are calculated in one block to
    ! ensure that incompatibility errors relating to the user's parallelisation
    ! come early in the calculation.
    call internal_distribute_ngwfs(pub_frag_data(0)%mdl%par, &
         pub_frag_data(0)%mdl, &
         pub_frag_data(0)%ngwf_basis(1),0)


  end subroutine eda_prep_super

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eda_prep_frags(full_mdl)

    !======================================================================!
    ! This suboutine prepares the fragment data objects.                   !
    !======================================================================!
    ! Written by Max Phipps, August 2013-2015                              !
    ! Modified for embedding by Joseph Prentice, September 2018            !
    !======================================================================!

    use utils, only: utils_alloc_check
    use rundat, only: pub_frag_iatm
    use fragment_data, only: pub_frag_data
    use parallel_strategy, only: pub_par
    use model_type, only: MODEL

    implicit none


    ! Arguments
    type(MODEL), target, intent(inout)  :: full_mdl

    ! Local Variables
    integer :: at ! atom counter
    integer :: it, at_id_min, at_id_max, ierr

    integer :: row
    logical :: new_species
    integer :: pseudo_counter, iat

    ! number of atoms used for supermolecule/fragment initialisation loop:
    integer :: num_atoms


    ! mjsp: Fragment initialisations:

    at_id_min = 1
    do it=1,pub_frag_iatm(0) ! mjsp: loop over fragments


       ! jcap: need to copy everything from full model, including regions
       ! - allocate in each region

       ! number of atoms in fragment/supermolecule
       num_atoms = pub_frag_iatm(it)
       pub_frag_data(it)%mdl%nat = num_atoms
       pub_frag_data(it)%mdl%nat_classical = 0
       pub_frag_data(it)%mdl%nsub = 1
       ! atom ranges
       at_id_max = at_id_min + num_atoms - 1

       ! allocate fragment storage for elements
       allocate(pub_frag_data(it)%mdl%elements(1:num_atoms),stat=ierr)
       call utils_alloc_check('eda_task','pub_frag_data(it)%mdl%elements',ierr)

       ! mjsp: Fill the fragment elements
       pub_frag_data(it)%mdl%elements(1:num_atoms) = &
            full_mdl%elements(at_id_min:at_id_max)

       ! jcap: need to reduce global_atom_number by the number of
       ! atoms in fragments already treated
       if (it.gt.1) &
            pub_frag_data(it)%mdl%elements(1:num_atoms)%global_atom_number = &
            pub_frag_data(it)%mdl%elements(1:num_atoms)%global_atom_number - &
            sum( pub_frag_iatm(1:(it-1)) )

       ! jcap: allocate regions - only compatible with 1 region
       allocate(pub_frag_data(it)%mdl%regions(1),stat=ierr)
       call utils_alloc_check('eda_task','pub_frag_data(it)%mdl%regions',ierr)

       pub_frag_data(it)%mdl%regions(1)%reg_count = 1

       ! allocate fragment storage for elements
       allocate(pub_frag_data(it)%mdl%regions(1)%elements(1:num_atoms),&
            stat=ierr)
       call utils_alloc_check('eda_task','pub_frag_data(it)%mdl%elements',&
            ierr)

       ! mjsp: Fill the fragment elements
       do iat=1,num_atoms
          pub_frag_data(it)%mdl%regions(1)%elements(iat) = &
               pub_frag_data(it)%mdl%elements(iat)
       end do

       ! jcap: set up pointer for mdl%par
       pub_frag_data(it)%mdl%par => pub_frag_data(it)%mdl%regions(1)%par

       pub_frag_data(it)%mdl%par%nat = num_atoms
       pub_frag_data(it)%mdl%par%nat_classical = 0

       ! mjsp: Modifications must be made to correct for fragment models:

 !!!!!!!!!!!!!!!!!!!!!!!!!!! NUM_NGWFS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! update par%num_ngwfs for this fragment:
       ! loop over atoms on this fragment and count up the total number of ngwfs
       ! on the fragment
       pub_frag_data(it)%mdl%par%num_ngwfs = 0
       do at=1,num_atoms
          pub_frag_data(it)%mdl%par%num_ngwfs = &
               pub_frag_data(it)%mdl%par%num_ngwfs + &
               pub_frag_data(it)%mdl%regions(1)%elements(at)%nfunctions
       end do


 !!!!!!!!!!!!!!!!!!!!!!!!!!! NUM_SPECIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! modified from rundat_blocks_mod.F90:
       ! mjsp: Re-initialise total number of different atomic species
       !       (as listed in % block species).
       pub_frag_data(it)%mdl%par%num_species = 1
       do at=2,num_atoms
          new_species = .true.
          do row=1,at-1
             if (pub_frag_data(it)%mdl%regions(1)%elements(at)%species_number &
                  == pub_frag_data(it)%mdl%regions(1)%elements(row)%species_number) &
                new_species = .false.
          end do

          if (new_species) &
             pub_frag_data(it)%mdl%par%num_species &
                  = pub_frag_data(it)%mdl%par%num_species + 1
       end do

       allocate(pub_frag_data(it)%mdl%regions(1)%species(&
            pub_frag_data(it)%mdl%par%num_species),stat=ierr)
       call utils_alloc_check('eda_task','pub_frag_data(it)%mdl%species',&
            ierr)

 !!!!!!!!!!!!!!!!!!!!!!!!!!! NUM_PSPECIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! modified from rundat_blocks_mod.F90:
          ! cks: initialise total number of different atomic pseudopotential species
          ! cks, 25/1/2004: A distint "species" is determined by the
          ! cks: pseudopotential name rather than the atomic number

       pub_frag_data(it)%mdl%par%num_pspecies = 1
       do at=2,num_atoms

          new_species = .true.
          do row=1,at-1
             if (full_mdl%species(&
                  pub_frag_data(it)%mdl%regions(1)%elements(at)&
                  %species_number)%pseudo_name &
                  == full_mdl%species(&
                  pub_frag_data(it)%mdl%regions(1)%elements(row)&
                  %species_number)%pseudo_name) &
                new_species = .false.
          end do

          if (new_species) &
             pub_frag_data(it)%mdl%par%num_pspecies &
                  = pub_frag_data(it)%mdl%par%num_pspecies + 1
       end do

       ! jcap: assign mdl%num_species
       pub_frag_data(it)%mdl%num_pspecies=pub_frag_data(it)%mdl%par%num_pspecies

 !!!!!!!!!!!!!!!!!!!!!!!!!!! SPECIES_NUMBER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! mjsp: Cull unrequired species numbers from the elements:

       ! Init:
       pseudo_counter = 1

       ! Iterate over the supermolecule species (the species will be culled to
       ! the range of num_species found in the fragment here)
       do row=1,full_mdl%par%num_species

          ! mjsp: Find an example of this species
          do iat=1,pub_frag_data(it)%mdl%par%nat
             if (pub_frag_data(it)%mdl%elements(iat)%species_number==row) exit
          end do

          ! mjsp: If an example has been found
          if (iat .le. pub_frag_data(it)%mdl%par%nat) then

             ! mjsp: Update any element with species_number=row to
             ! species_number=-pseudo_counter  (temporarily negative to avoid
             ! double correcting)
             do iat=1,pub_frag_data(it)%mdl%par%nat
                if (pub_frag_data(it)%mdl%elements(iat)%species_number==row) &
                     pub_frag_data(it)%mdl%elements(iat)%species_number=-pseudo_counter
             end do

             pseudo_counter = pseudo_counter + 1
          end if

       end do

       ! mjsp: Re-correct the fixed species_number values from negative
       ! mjsp: back to positive
       pub_frag_data(it)%mdl%elements(:)%species_number = &
            abs(pub_frag_data(it)%mdl%elements(:)%species_number)

       ! jcap: copy this back to the region element array
       do iat=1,pub_frag_data(it)%mdl%nat
          pub_frag_data(it)%mdl%regions(1)%elements(iat)%species_number = &
               pub_frag_data(it)%mdl%elements(iat)%species_number
       end do

 !!!!!!!!!!!!!!!!!!!!!!!!!!! SPECIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! mjsp: Allocate fragment species
       allocate(pub_frag_data(it)%mdl%species( &
            1:pub_frag_data(it)%mdl%par%num_species),stat=ierr)
       call utils_alloc_check('eda_task','pub_frag_data(it)%mdl%species',ierr)

       ! mjsp: Copy species from the supermolecule model to the fragment, whilst
       ! culling to construct a (correctly ordered) minimal set:
       ! TODO: This loop could be made more efficient...

!       ! Counter of species copied to the fragment set:
!       pseudo_counter = 0

       ! loop fragment atoms
       do iat=1,pub_frag_data(it)%mdl%par%nat

          ! mjsp: Find this species in the supermolecule species set
          do row=1,full_mdl%par%num_species
             if (pub_frag_data(it)%mdl%elements(iat)%species_id== &
                  full_mdl%species(row)%species_id) &
                exit
          end do

          ! mjsp: Retain this species in our fragment model
          pub_frag_data(it)%mdl%species( &
               pub_frag_data(it)%mdl%elements(iat)%species_number ) = &
               full_mdl%species( row )

          ! jcap: also do this for the regions array
          pub_frag_data(it)%mdl%regions(1)%species( &
               pub_frag_data(it)%mdl%regions(1)%elements(iat)%species_number ) = &
               full_mdl%species( row )

          pseudo_counter = pseudo_counter + 1

!          ! mjsp: Check if we have copied over all species to the fragment set:
!          if (pseudo_counter .eq. pub_frag_data(it)%mdl%par%num_species) &
!             exit

       end do

       ! jcap: let's initialise par now - bits have already been
       ! initialised above, but this shouldn't affect those
       call internal_initialise_par(pub_frag_data(it)%mdl)

       ! jcap: allocate ngwf_basis array
       allocate(pub_frag_data(it)%ngwf_basis(1),stat=ierr)
       call utils_alloc_check('eda_task','pub_frag_data(it)%ngwf_basis',ierr)

       ! Distribute the NGWFs. We loop over each system here
       ! before we have invested any calculation time in order to
       ! check that we won't be oversubscribed with cores later on
       ! in the EDA calculation.
       call internal_distribute_ngwfs(pub_frag_data(it)%mdl%par, &
            pub_frag_data(it)%mdl, &
            pub_frag_data(it)%ngwf_basis(1),it)

       ! iterate fragment atom counter
       at_id_min = at_id_max + 1

       ! jcap: point pub_par to the appropriate parallel strategy
       pub_par => pub_frag_data(it)%mdl%regions(1)%par

    end do


  end subroutine eda_prep_frags

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine internal_distribute_ngwfs(par,mdl,ngwf_basis, frag_it)

    ! jcap: modified for embedding by Joseph Prentice, September 2018

    use constants, only: NORMAL, CRLF, stdout
    use rundat, only: pub_output_detail
    use parallel_strategy, only: PARAL_INFO, parallel_strategy_distr_atoms
    use function_basis, only: FUNC_BASIS, function_basis_allocate, &
         function_basis_distribute
    use model_type, only: MODEL
    use comms, only: pub_on_root, comms_barrier
    use utils, only: utils_assert

    ! Arguments
    type(PARAL_INFO), intent(inout), target :: par
    type(MODEL), intent(inout)      :: mdl
    type(FUNC_BASIS), intent(inout) :: ngwf_basis
    integer, intent(in)          :: frag_it

    ! Local Variables
    type(PARAL_INFO), pointer :: point_par

    ! rc2013: hack to avoid silly changes
    point_par=>par

    ! Adopt parallel strategy for atom distribution
    call comms_barrier
    if (pub_on_root .and. (pub_output_detail >= NORMAL)) then
       if (frag_it==0) then
         write(stdout,'(/a)',advance='no') 'Determining parallel &
               &strategy for supermolecule ...'
       else
         write(stdout,'(a,i3,a)',advance='no') 'Determining parallel &
               &strategy for fragment ID: ',frag_it,'...'
       end if
    end if

    call parallel_strategy_distr_atoms(par,mdl%regions(1)%elements,mdl%cell)

    if (pub_on_root .and. (pub_output_detail >= NORMAL)) &
         write(stdout,'(a)') '... done'

    ! Distribute the NGWFs
    call function_basis_allocate(ngwf_basis,par%num_ngwfs,'ngwfs',point_par)
    call function_basis_distribute(ngwf_basis,mdl%regions(1)%elements,point_par)


    ! Test that this system is compatible with the number of cores
    ! the user has chosen (modified from energy_and_force_calculate).
    ! Exit if there is a proc with no NGWFs:
    call utils_assert (minval(ngwf_basis%num_on_proc) >= 1, &
        'Too many processors (MPI ranks) for this system size. Sorry!'//CRLF//&
        'Your options include: '//CRLF//&
        '1) Reducing the number of processors (MPI ranks), '//CRLF//&
        '2) Increasing the number of atoms, '//CRLF//&
        '3) Using OpenMP parallelisation to your advantage, by reducing'//CRLF//&
        '   the number of MPI ranks, and simultaneously using threads_max'//CRLF//&
        '   and other threads_* keywords in your input file to have each'//CRLF//&
        '   MPI rank spawn multiple threads.'//CRLF//CRLF//&
        'Occasionally, turning off the space-filling curve by specifying'//CRLF//&
        '"use_space_filling_curve F" in the input file might help --'//CRLF//&
        'but only if the resultant atom ordering better lends itself to '//CRLF//&
        'distributing across processors. If you do this, be aware that'//CRLF//&
        'if the ions move (geometry optimisation, TS search, MD), the'//CRLF//&
        'ordering may once again become unfavourable and this error might'//CRLF//&
        'come back in the course of the calculation.')


  end subroutine internal_distribute_ngwfs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine internal_warnings(full_mdl)

    !======================================================================!
    ! This suboutine provides the user with input/system errors and/or     !
    ! warnings.                                                            !
    !======================================================================!
    ! Written by Max Phipps, February 2016                                 !
    !======================================================================!

    use comms, only: pub_on_root
    use constants, only: CRLF, stdout
    use rundat, only: pub_exact_lnv, pub_all_tasks, pub_charge, &
         pub_eda_preptool, pub_eda_frag_isol_pol, &
         pub_eda_frag_isol_ct, pub_frag_iatm, pub_frag_charge, &
         maxit_lnv, minit_lnv, pub_maxit_ngwf_cg, pub_xc_functional
    use model_type, only: MODEL
    use utils, only: utils_assert

    ! Arguments
    type(MODEL), intent(inout)      :: full_mdl


    if (pub_eda_preptool) then
       call utils_assert(.not.(minit_lnv .gt. 0 .or. maxit_lnv .gt. 0), &
            'Please ensure that the number of LNV iterations is'//CRLF// &
            'set to zero for an EDA preparation stage calculation.')
       call utils_assert(.not.(pub_maxit_ngwf_cg .gt. 0), &
            ' Please ensure that the number of NGWF iterations is'//CRLF// &
            'set to zero for an EDA preparation stage calculation.')
    end if

    call utils_assert((pub_charge == pub_frag_charge(0)), &
         'Please ensure that the charge keyword &
         &is set to equal the total of the fragment'//CRLF//'charges.')

    call utils_assert((full_mdl%nat == sum(pub_frag_iatm(1:pub_frag_iatm(0)))),&
         'Please ensure that the total number of atoms defined by &
         &the fragment input block'//CRLF//'eda_iatm is in agreement with &
         &the total number of atoms in the system.')

    call utils_assert(.not.(pub_eda_frag_isol_ct.and.(pub_frag_iatm(0).eq.2)),&
         'Charge transfer analysis requires greater than two fragments. &
         & Please'//CRLF//'&
         &increase the number of fragments, or set EDA_FRAG_ISOL_CT to false.')

    call utils_assert(.not.(pub_xc_functional .eq. 'VDWDF' .or. &
                            pub_xc_functional .eq. 'VDWDF2' .or. &
                            pub_xc_functional .eq. 'VDWC09' .or. &
                            pub_xc_functional .eq. 'C09-2' .or. &
                            pub_xc_functional .eq. 'OPTPBE' .or. &
                            pub_xc_functional .eq. 'OPTB88' .or. &
                            pub_xc_functional .eq. 'PBEK1' .or. &
                            pub_xc_functional .eq. 'VV10' .or. &  ! (dispersion)
                            pub_xc_functional .eq. 'AVV10S' ), &
          'The following XC functionals are not available&
          & for use with SCF-MI and/or EDA:'//CRLF//'&
          &VDWDF, VDWDF2, VDWC09, C09-2, OPTPBE, OPTB88, &
          &PBEK1, VV10, AVV10S .'//CRLF//'&
          &Please modify the XC functional.')

    if ((pub_eda_frag_isol_pol) .and. (pub_exact_lnv) .and. (pub_on_root) .and. &
         (index(pub_all_tasks,'EDA_PREP') .eq. 0)) &
       write(stdout,'(a)') 'WARNING: Exact LNV will be turned off during &
            &SCF-MI and/or EDA polarisation'//CRLF//'calculations.'


  end subroutine internal_warnings

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine internal_check_supermolecule(super_mdl, user_forced_spin_polarised)

    !======================================================================!
    ! This suboutine initialises and exits the supermolecule cell to       !
    ! check for issues.                                                    !
    !======================================================================!
    ! Written by Max Phipps, February 2016                                 !
    ! Modified for embedding by Joseph Prentice, September 2018            !
    !======================================================================!

    use constants, only: DP, stdout
    use comms, only: pub_on_root
    use cutoff_coulomb, only: cutoff_coulomb_check_boundaries
    use energy_and_force, only: &
         energy_and_force_init_cell, energy_and_force_exit_cell
    use parallel_strategy, only: pub_par, parallel_strategy_check_atoms
    use pbc_corrections, only: pbc_corr_check_cell_size
    use rundat, only: pub_spin_polarised, &
         pub_coulomb_cutoff, pub_mt_cutoff, &
         pub_charge, pub_frag_charge
    use model_type, only: MODEL
    use utils, only: utils_banner

    ! Arguments
    type(MODEL), target, intent(inout)      :: super_mdl
    logical, intent(in)                     :: user_forced_spin_polarised


    if(pub_on_root) then
       write(stdout,'(/a)') utils_banner('=')
       write(stdout,'(a)') &
            ' EDA: Performing initialisation checks for supermolecule system'
       write(stdout,'(a/)') utils_banner('=')
    end if

    if (.not. user_forced_spin_polarised) pub_spin_polarised = .false.
    pub_par=>super_mdl%par
    pub_charge = pub_frag_charge(0)

    call energy_and_force_init_cell(super_mdl)

    ! mjsp: for cutoff coulomb, check no NGWFs extend beyond cell boundaries
    if (pub_coulomb_cutoff) then
       if (pub_on_root) &
            write(stdout,'(a)',advance='no') 'Checking NGWF boundaries...'
       call cutoff_coulomb_check_boundaries(super_mdl)
       if (pub_on_root) write(stdout,*)' done'
    end if

    ! mjsp: make sure that no atom is outside the simulation cell
    if (pub_on_root) write(stdout,'(a)',advance='no') &
         'Checking atom boundaries...'
    call parallel_strategy_check_atoms(super_mdl%par,&
         super_mdl%regions(1)%elements,&
         super_mdl%cell)
    if (pub_on_root) write(stdout,*)' done'

    ! mjsp: make sure that the cell is large enough if
    ! mjsp: MT correction is in effect
    if (pub_on_root .and. pub_mt_cutoff /= 0.0_DP) then
       if (pub_on_root) write(stdout,'(a)',advance='no') &
            'Checking periodic boundary conditions...'
       call pbc_corr_check_cell_size(super_mdl%regions(1)%elements,&
            super_mdl%cell)
       if (pub_on_root) write(stdout,*)' done'
    end if

    call energy_and_force_exit_cell(super_mdl)

    if(pub_on_root) then
       write(stdout,'(/a)') utils_banner('=')
       write(stdout,'(a)') &
            ' EDA: Finished initialisation checks for supermolecule system'
       write(stdout,'(a/)') utils_banner('=')
    end if

  end subroutine internal_check_supermolecule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine internal_check_fragment(frag_id, frag_mdl, &
       user_forced_spin_polarised)

    !======================================================================!
    ! This suboutine initialises and exits the given fragment cell to      !
    ! check for issues.                                                    !
    !======================================================================!
    ! Written by Max Phipps, February 2016                                 !
    !======================================================================!

    use constants, only: DP, stdout
    use comms, only: pub_on_root
    use cutoff_coulomb, only: cutoff_coulomb_check_boundaries
    use energy_and_force, only: &
         energy_and_force_init_cell, energy_and_force_exit_cell
    use parallel_strategy, only: pub_par, parallel_strategy_check_atoms
    use pbc_corrections, only: pbc_corr_check_cell_size
    use rundat, only: pub_spin_polarised, &
         pub_coulomb_cutoff, pub_mt_cutoff, &
         pub_charge, pub_frag_charge
    use model_type, only: MODEL
    use utils, only: utils_banner

    ! Arguments
    integer, intent(in)                     :: frag_id
    type(MODEL), target, intent(inout)      :: frag_mdl
    logical, intent(in)                     :: user_forced_spin_polarised


    if(pub_on_root) then
       write(stdout,'(//a)') utils_banner('=')
       write(stdout,'(a,i3)') &
            ' EDA: Performing initialisation checks for fragment #', frag_id
       write(stdout,'(a/)') utils_banner('=')
    end if

    if (.not. user_forced_spin_polarised) pub_spin_polarised = .false.
    pub_par=>frag_mdl%par
    pub_charge = pub_frag_charge(frag_id)

    call energy_and_force_init_cell(frag_mdl)

    ! mjsp: for cutoff coulomb, check no NGWFs extend beyond cell boundaries
    if (pub_coulomb_cutoff) then
       if (pub_on_root) &
            write(stdout,'(/a)',advance='no') 'Checking NGWF boundaries...'
       call cutoff_coulomb_check_boundaries(frag_mdl)
       if (pub_on_root) write(stdout,*)' done'
    end if

    ! mjsp: make sure that no atom is outside the simulation cell
    if (pub_on_root) &
         write(stdout,'(a)',advance='no') 'Checking atom boundaries...'
    call parallel_strategy_check_atoms(frag_mdl%par,&
                                       frag_mdl%regions(1)%elements,&
                                       frag_mdl%cell)
    if (pub_on_root) write(stdout,*)' done'

    ! mjsp: make sure that the cell is large enough if
    ! mjsp: MT correction is in effect
    if (pub_on_root .and. pub_mt_cutoff /= 0.0_DP) then
       if (pub_on_root) write(stdout,'(a)',advance='no') &
            'Checking periodic boundary conditions...'
       call pbc_corr_check_cell_size(frag_mdl%regions(1)%elements,&
            frag_mdl%cell)
       if (pub_on_root) write(stdout,*)' done'
    end if

    call energy_and_force_exit_cell(frag_mdl)


  end subroutine internal_check_fragment

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine internal_copy_par(par_out, par_in)

    !======================================================================!
    ! This suboutine allocates and copies data from one par to another     !
    !======================================================================!
    ! Written by Joseph Prentice, September 2018                           !
    !======================================================================!

    use utils, only: utils_alloc_check
    use parallel_strategy, only: PARAL_INFO

    ! Arguments
    type(PARAL_INFO), intent(inout) :: par_out
    type(PARAL_INFO), intent(in)    :: par_in

    ! Local variables
    integer :: ierr

    if (allocated(par_in%num_atoms_on_proc).and.&
         .not.allocated(par_out%num_atoms_on_proc)) then
       allocate(par_out%num_atoms_on_proc(size(par_in%num_atoms_on_proc)),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','num_atoms_on_proc',&
            ierr)
    end if

    if (allocated(par_in%first_atom_on_proc).and.&
         .not.allocated(par_out%first_atom_on_proc)) then
       allocate(par_out%first_atom_on_proc(size(par_in%first_atom_on_proc)),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','first_atom_on_proc',&
            ierr)
    end if

    if (allocated(par_in%proc_of_atom).and.&
         .not.allocated(par_out%proc_of_atom)) then
       allocate(par_out%proc_of_atom(size(par_in%proc_of_atom)),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','proc_of_atom',&
            ierr)
    end if

    if (allocated(par_in%num_hub_atoms_on_proc).and.&
         .not.allocated(par_out%num_hub_atoms_on_proc)) then
       allocate(par_out%num_hub_atoms_on_proc(size(par_in%num_hub_atoms_on_proc)),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','num_hub_atoms_on_proc',&
            ierr)
    end if

    if (allocated(par_in%proc_of_hub_atom).and.&
         .not.allocated(par_out%proc_of_hub_atom)) then
       allocate(par_out%proc_of_hub_atom(size(par_in%proc_of_hub_atom)),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','proc_of_hub_atom',&
            ierr)
    end if

    if (allocated(par_in%hub_atoms_on_proc).and.&
         .not.allocated(par_out%hub_atoms_on_proc)) then
       allocate(par_out%hub_atoms_on_proc(size(par_in%hub_atoms_on_proc(:,1)),&
            size(par_in%hub_atoms_on_proc(1,:))),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','hub_atoms_on_proc',&
            ierr)
    end if

    if (allocated(par_in%hub_atom_orig_atom).and.&
         .not.allocated(par_out%hub_atom_orig_atom)) then
       allocate(par_out%hub_atom_orig_atom(size(par_in%hub_atom_orig_atom)),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','hub_atom_orig_atom',&
            ierr)
    end if

    if (allocated(par_in%atoms_on_proc).and.&
         .not.allocated(par_out%atoms_on_proc)) then
       allocate(par_out%atoms_on_proc(size(par_in%atoms_on_proc(:,1)),&
            size(par_in%atoms_on_proc(1,:))),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','atoms_on_proc',&
            ierr)
    end if

    if (allocated(par_in%orig_atom).and.&
         .not.allocated(par_out%orig_atom)) then
       allocate(par_out%orig_atom(size(par_in%orig_atom)),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','orig_atom',&
            ierr)
    end if

    if (allocated(par_in%distr_atom).and.&
         .not.allocated(par_out%distr_atom)) then
       allocate(par_out%distr_atom(size(par_in%distr_atom)),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','distr_atom',&
            ierr)
    end if

    if (allocated(par_in%elements_on_proc).and.&
         .not.allocated(par_out%elements_on_proc)) then
       allocate(par_out%elements_on_proc(size(par_in%elements_on_proc)),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','elements_on_proc',&
            ierr)
    end if

    if (allocated(par_in%num_kpoints_on_proc).and.&
         .not.allocated(par_out%num_kpoints_on_proc)) then
       allocate(par_out%num_kpoints_on_proc(size(par_in%num_kpoints_on_proc)),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','num_kpoints_on_proc',&
            ierr)
    end if

    if (allocated(par_in%first_kpoint_on_proc).and.&
         .not.allocated(par_out%first_kpoint_on_proc)) then
       allocate(par_out%first_kpoint_on_proc(size(par_in%first_kpoint_on_proc)),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','first_kpoint_on_proc',&
            ierr)
    end if

    if (allocated(par_in%proc_of_kpoint).and.&
         .not.allocated(par_out%proc_of_kpoint)) then
       allocate(par_out%proc_of_kpoint(size(par_in%proc_of_kpoint)),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','proc_of_kpoint',&
            ierr)
    end if

    if (allocated(par_in%kp_indices_on_proc).and.&
         .not.allocated(par_out%kp_indices_on_proc)) then
       allocate(par_out%kp_indices_on_proc(size(par_in%kp_indices_on_proc(:,1)),&
            size(par_in%kp_indices_on_proc(1,:))),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','kp_indices_on_proc',&
            ierr)
    end if

    if (allocated(par_in%kpoints_on_proc).and.&
         .not.allocated(par_out%kpoints_on_proc)) then
       allocate(par_out%kpoints_on_proc(size(par_in%kpoints_on_proc)),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','kpoints_on_proc',&
            ierr)
    end if

    if (allocated(par_in%first_elem_on_atom).and.&
         .not.allocated(par_out%first_elem_on_atom)) then
       allocate(par_out%first_elem_on_atom(size(par_in%first_elem_on_atom(:,1)),&
            size(par_in%first_elem_on_atom(1,:))),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','first_elem_on_atom',&
            ierr)
    end if

    if (allocated(par_in%first_elem_on_proc).and.&
         .not.allocated(par_out%first_elem_on_proc)) then
       allocate(par_out%first_elem_on_proc(size(par_in%first_elem_on_proc(:,1)),&
            size(par_in%first_elem_on_proc(1,:))),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','first_elem_on_proc',&
            ierr)
    end if

    if (allocated(par_in%num_elems_on_atom).and.&
         .not.allocated(par_out%num_elems_on_atom)) then
       allocate(par_out%num_elems_on_atom(size(par_in%num_elems_on_atom(:,1)),&
            size(par_in%num_elems_on_atom(1,:))),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','num_elems_on_atom',&
            ierr)
    end if

    if (allocated(par_in%num_elems_on_proc).and.&
         .not.allocated(par_out%num_elems_on_proc)) then
       allocate(par_out%num_elems_on_proc(size(par_in%num_elems_on_proc(:,1)),&
            size(par_in%num_elems_on_proc(1,:))),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','num_elems_on_proc',&
            ierr)
    end if

    if (allocated(par_in%atom_of_elem).and.&
         .not.allocated(par_out%atom_of_elem)) then
       allocate(par_out%atom_of_elem(size(par_in%atom_of_elem(:,1)),&
            size(par_in%atom_of_elem(1,:))),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','atom_of_elem',&
            ierr)
    end if

    if (allocated(par_in%proc_of_elem).and.&
         .not.allocated(par_out%proc_of_elem)) then
       allocate(par_out%proc_of_elem(size(par_in%proc_of_elem(:,1)),&
            size(par_in%proc_of_elem(1,:))),&
            stat=ierr)
       call utils_alloc_check('internal_copy_par','proc_of_elem',&
            ierr)
    end if

    par_out = par_in

  end subroutine internal_copy_par

!-------------------------------------------------------------------------------

    !==========================================================================!
    ! Allocate the number of procs to be included for each parallel strategy.  !
    !--------------------------------------------------------------------------!
    ! Written by Robert Charlton, 11/05/2018.                                  !
    ! Copied from energy_and_force by Joseph Prentice, September 2018          !
    !==========================================================================!

    subroutine internal_initialise_par(mdl)

      use constants, only: DP, CRLF
      use comms, only: pub_total_num_procs
      use model_type, only: MODEL
      use rundat, only: pub_parallel_scheme
      use utils, only: utils_assert

      implicit none

      type(MODEL), intent(inout) :: mdl

      ! local variables
      integer :: procs_per_par(mdl%nsub)
      integer :: proc_to_give(1), ii, isub
      real(kind=DP) :: remainder(mdl%nsub)

      do isub=1,mdl%nsub
         ! rc2013: set the region index
         mdl%regions(isub)%par%par_index = isub
      end do

      ! rc2013: decide how to distribute the procs
      select case(pub_parallel_scheme)
      ! Split processore equally between regions
      case('SENATE')
         ! Make sure we have enough resources for this
         call utils_assert(pub_total_num_procs .ge. mdl%nsub, &
              ' Error in energy_and_force_mod: pub_parallel_scheme: '//CRLF// &
              '   SENATE requested, but the number of regions in '//CRLF// &
              '   %block species_ngwf_regions exceeds the number of '//CRLF// &
              '   processors. Please either allocate additional '//CRLF// &
              '   processors,  or set pub_parallel_scheme: NONE in input file.')
         procs_per_par = pub_total_num_procs/mdl%nsub

      ! Distribute processors proportionally between regions
      case('HOUSE')
         ! Make sure we have enough resources for this
         call utils_assert(pub_total_num_procs .ge. mdl%nsub, &
              ' Error in energy_and_force_mod: pub_parallel_scheme: '//CRLF// &
              '   HOUSE requested, but the number of regions in '//CRLF// &
              '   %block species_ngwf_regions exceeds the number of '//CRLF// &
              '   processors. Please either allocate additional '//CRLF// &
              '   processors,  or set pub_parallel_scheme: NONE in input file.')
         procs_per_par = 1

         ! Calculate the number of procs per par with Huntingdon-Hill algorithm
         ! Each region already has 1 proc, so subtract that from total
         do ii=1,pub_total_num_procs-mdl%nsub
            do isub=1,mdl%nsub
               remainder(isub) = mdl%regions(isub)%par%nat/ &
                    sqrt(real(procs_per_par(isub)*(procs_per_par(isub)+1)))
            end do
            ! Give the next proc to the region with the highest remainder
            proc_to_give = maxloc(remainder)
            procs_per_par(proc_to_give(1)) = procs_per_par(proc_to_give(1)) + 1
         end do

      ! Default: spread atoms across all processors regardless of region
      case('NONE')
         procs_per_par = pub_total_num_procs
      end select

      mdl%regions(1)%par%first_proc = 0
      mdl%regions(1)%par%num_procs = procs_per_par(1)
      do isub=2,mdl%nsub
         if(pub_parallel_scheme == 'NONE') then
            mdl%regions(isub)%par%first_proc = 0
         else
            mdl%regions(isub)%par%first_proc = &
                 mdl%regions(isub-1)%par%first_proc + procs_per_par(isub-1)
         end if
         mdl%regions(isub)%par%num_procs = procs_per_par(isub)
      end do

      if(mdl%nsub == 1) then
         call utils_assert(mdl%regions(1)%par%first_proc == 0, &
              'Error in internal_initialise_par: failure to allocate &
              &parallel strategies.')
         call utils_assert(procs_per_par(1) == pub_total_num_procs, &
              'Error in internal_initialise_par: there is only one parallel &
              &strategy, but the total number of procs allocated to this par &
              &does not match the number of procs available.')
      endif

    end subroutine internal_initialise_par

end module eda_main

