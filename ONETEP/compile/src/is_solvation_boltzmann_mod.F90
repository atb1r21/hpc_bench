!=============================================================================!
!                     S O L V A T I O N _ B O L T Z M A N N                   !
!=============================================================================!
! This module collects the Boltzmann (ion-concentration-dependent) facet of   !
! implicit solvent. It mostly encompasses variables, datatypes and subroutines!
! related to the Boltzmann ions and the steric potential.                     !
!-----------------------------------------------------------------------------!
! Written by Jacek Dziedzic, with help from Lucian Anton, 10-11/2014.         !
! Major changes by Jacek Dziedzic, 2019.08-11, implementing formulas derived  !
! by Arihant Bhandari, Chris-Kriton Skylaris, Lucian Anton and Jacek Dziedzic.!
!-----------------------------------------------------------------------------!
! Recommended reading:                                                        !
! [1] Nattino, Truscott, Marzari, Andreussi, "Continuum models of the electro-!
!     chemical diffuse layer in electronic-structure calculations",           !
!     J. Chem. Phys. 150, 041722 (2019).                                      !
! [2] C.G. Gray, P.J. Stiles, "Nonlinear electrostatics: the Poisson-Boltzmann!
!     equation", European Journal of Physics 39 (5) (2018).                   !
!-----------------------------------------------------------------------------!
! Algorithms implemented here follow:                                         !
! [3] J. Dziedzic, A. Bhandari, L. Anton, C. Peng, J.C. Womack, M. Famili,    !
!     D. Kramer, and C.-K. Skylaris, "Practical Approach to Large-Scale       !
!     Electronic Structure Calculations in Electrolyte Solutions via          !
!     Continuum-Embedded Linear-Scaling Density Functional Theory",           !
!     J. Phys. Chem. C 124, 14, 7860-7872 (2020).                             !
! [4] A. Bhandari, L. Anton, J. Dziedzic, C. Peng, D. Kramer, and             !
!     C.-K. Skylaris, "Electronic Structure Calculations in Electrolyte       !
!     Solutions: Methods for Neutralization of Extended Charged Interfaces",  !
!     J. Chem. Phys. [accepted (2020)].                                       !
!-----------------------------------------------------------------------------!

module is_solvation_boltzmann

  use constants, only: DP, stdout
#ifdef HAVE_DL_MG
  use dl_mg !! External dependency
  use dl_mg_nonlinear_model, only: max_expcap   !! External dependency
  use dl_mg_params, only: dl_mg_neutralise_none, &
       dl_mg_neutralise_with_jellium_uniform, &
       dl_mg_neutralise_with_jellium_vacc, dl_mg_neutralise_with_ions_fixed, &
       dl_mg_neutralise_with_ions_auto, dl_mg_neutralise_with_ions_auto_linear
#endif

  implicit none

  private

  !---------------------------------------------------------------------------!
  !                       P u b l i c   r o u t i n e s                       !
  !                              a n d   t y p e s                            !
  !---------------------------------------------------------------------------!

  ! The call sequence is as follows:
  !   1) rundat_blocks calls init_ions,
  !   2) energy_and_force_calculate calls init,
  !   3) multigrid_methods calls main, possibly shifted_hartree,
  !   4) hamiltonian, e&f and multigrid_methods call energy_contrib,
  !   5) energy_and_force_calculate calls exit.
  !   2..5 are called multiple times if ions move.
  !   6) energy_and_force_exit_cell calls exit_ions.

  public :: implicit_solvent_boltzmann_init_ions
  public :: implicit_solvent_boltzmann_init
  public :: implicit_solvent_boltzmann_defect
  public :: implicit_solvent_boltzmann_main
  public :: implicit_solvent_boltzmann_energy_contrib
  public :: implicit_solvent_boltzmann_shifted_hartree
  public :: implicit_solvent_boltzmann_ion_conc
  public :: implicit_solvent_boltzmann_ion_stats
  public :: implicit_solvent_boltzmann_exit
  public :: implicit_solvent_boltzmann_exit_ions
  public :: implicit_solvent_boltzmann_debug_dump
  public :: implicit_solvent_boltzmann_mg_neutralisation_method

  !---------------------------------------------------------------------------!
  !                      P u b l i c   v a r i a b l e s                      !
  !---------------------------------------------------------------------------!
  ! jd: The maximum number of (kinds of) mean-field ions in implicit solvent
  !     Needed in rundat_blocks and to dimension IS_PBE_MG_RETURN.
  integer, parameter, public   :: implicit_solvent_boltzmann_max_n_ions = 1000

  integer, parameter, public   :: n_boltzmann_energy_terms = 11

  ! jd: The calculated Debye length.
  !     Needed in multigrid_methods.
  real(kind=DP), save, public  :: implicit_solvent_boltzmann_debye_length

  ! ******************** Boltzmann ions ********************

  ! jd: Describes a mean-field ion in implicit solvent
  type BOLTZMANN_ION
     character(len=10) :: name
     real(kind=DP)     :: charge             ! jd: in |e|
     real(kind=DP)     :: neutralisation_ratio ! jd: generalisation of 'x' frac
     real(kind=DP)     :: concentration_user ! jd: Conc. specified by user (au)
     real(kind=DP)     :: concentration_adj  ! jd: } Conc. after adjustment in
  end type BOLTZMANN_ION                     ! jd: } DL_MG (counterions) (au)

  public :: BOLTZMANN_ION

  ! jd: All Boltzmann ion kinds in solvent. A copy is held on all procs.
  type(BOLTZMANN_ION), allocatable, save, protected, public :: boltzmann_ions(:)

  ! jd: Number of Boltzmann ions
  integer, save, protected, public   :: n_boltzmann_ions

  ! jd: Collects information relevant to PBE solvation returned by DL_MG.
  !     Cannot use n_boltzmann_ions for dimensioning (not a compile-time constant).
  type IS_PBE_MG_RETURN
     real(kind=DP) :: steric_weight_average
     real(kind=DP) :: betamu_elec(implicit_solvent_boltzmann_max_n_ions)
     real(kind=DP) :: used_ion_concentration(implicit_solvent_boltzmann_max_n_ions)
     real(kind=DP) :: used_neutralisation_ratios(implicit_solvent_boltzmann_max_n_ions)
  end type IS_PBE_MG_RETURN

  public :: IS_PBE_MG_RETURN

  ! ******************** Steric potential and weight ********************

  ! jd: Steric potential (slab-distributed on MG grid)
  real(kind=DP), allocatable, save, protected, public :: steric_pot(:,:,:)

  ! jd: Accessibility (slab-distributed on MG grid)
  real(kind=DP), allocatable, save, protected, public :: gamma(:,:,:)

#ifdef HAVE_DL_MG
  ! JCW: max_expcap_dl_mg is derived from the value of max_expcap in
  ! JCW: dl_mg_nonlinear_model and used if pub_is_pbe_exp_cap == 0.0_DP.
  ! JCW: max_expcap has the parameter attribute in dl_mg_nonlinear_model, so
  ! JCW: we can safely set the duplicate of it as a parameter.
  ! JCW: Must be public, since required to call solver routines in
  ! JCW: multigrid_methods.
  real(kind=DP), parameter, public   :: max_expcap_dl_mg = max_expcap
#else
  ! jd: Dummy value, only so the code compiles with no DL_MG
  real(kind=DP), parameter, public   :: max_expcap_dl_mg = -1.0_DP
#endif

  ! ---------------------------------------------------------------------------
  ! -------------------------- P R I V A T E ----------------------------------
  ! ---------------------------------------------------------------------------


  logical, save   :: boltzmann_ions_reported = .false.

  character(len=80), save :: steric_pot_type ! 'H': hard-core, 'S': soft-core, 'X': none
                                     ! 'M': smoothed hard-core

contains

  !---------------------------------------------------------------------------!
  !                      P u b l i c   s u b r o u t i n e s                  !
  !---------------------------------------------------------------------------!

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_boltzmann_init_ions(cbuf, dbuf, &
       runbl_n_sol_ions)
    !=========================================================================!
    ! Initializes the sol_ions structure based on raw data passed from        !
    ! rundat_blocks. The structure is allocated here, to be deallocated in    !
    ! implicit_solvent_exit.
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   cbuf, dbuf (in): See description in rundat_blocks.                    !
    !   runbl_n_sol_ions (in): The number of solvent ions from rundat_blocks. !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic in 07/2014.                                   !
    !=========================================================================!

    use constants, only: CONC_MOL_PER_L_TO_AU
    use rundat, only: pub_is_pbe
    use utils, only: utils_alloc_check, utils_assert

    implicit none

    ! jd: Arguments
    character(len=10), intent(in) :: cbuf(implicit_solvent_boltzmann_max_n_ions)
    real(kind=DP), intent(in)   :: dbuf(implicit_solvent_boltzmann_max_n_ions*3)
    integer                     :: runbl_n_sol_ions

    ! jd: Local variables
    integer :: ierr
    integer :: i
    character(len=*), parameter   :: myself = &
         'implicit_solvent_boltzmann_init_ions'

    ! -------------------------------------------------------------------------

    call utils_assert(pub_is_pbe /= 'NONE', trim(myself)//&
         ' must only be called when pub_is_pbe is not ''none''.')

    n_boltzmann_ions = runbl_n_sol_ions

    allocate(boltzmann_ions(n_boltzmann_ions), stat=ierr)
    call utils_alloc_check(myself,'boltzmann_ions',ierr)

    do i=1, n_boltzmann_ions
       boltzmann_ions(i)%name = cbuf(i)(1:10)
       boltzmann_ions(i)%charge = dbuf(i*3-2)
       boltzmann_ions(i)%concentration_user = CONC_MOL_PER_L_TO_AU * dbuf(i*3-1)
       boltzmann_ions(i)%neutralisation_ratio = dbuf(i*3)

       ! jd: Do not allow Poisson-Boltzmann calculations with zero concentration
       !     unless in the LINEARISED case (which is fine).
       if(pub_is_pbe == "FULL") then
          call utils_assert(boltzmann_ions(i)%concentration_user > 0.0_DP, &
               myself//": Boltzmann ion '"//trim(boltzmann_ions(i)%name)//"' &
               &has zero bulk concentration. All ions present must have &
               &non-zero concentrations. If your system has no Boltzmann ions, &
               &you must use 'is_pbe NONE'.")
       end if

       call utils_assert(boltzmann_ions(i)%neutralisation_ratio >= 0.0_DP, &
               myself//": Boltzmann ion '"//trim(boltzmann_ions(i)%name)//"' &
               &has a negative neutralisation fraction.")
    end do

  end subroutine implicit_solvent_boltzmann_init_ions
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_boltzmann_init(mdl)
    !=========================================================================!
    ! Initialises Boltzmann solvation in ONETEP.                              !
    ! - Reports on Boltzmann ions (these are initialised earlier from         !
    !   rundat_blocks -> sol_ions_report), but it's neater to display the     !
    !   report only now.                                                      !
    ! - Flips the ions' charges to internal charge convention.                !
    ! - Calculates steric potential and accessibility (gamma), allocating     !
    !   them, if necessary.                                                   !
    ! - Calculates the Debye length for a dielectric with the permittivity    !
    !   pub_is_bulk_permittivity                                              !
    ! ----------------------------------------------------------------------- !
    ! Notes:                                                                  !
    !   This can (and needs to) be called each time ions move, as the steric  !
    !   potential must be updated                                             !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   mdl (in/out): Current model. We need elements, cell, fine_grid,       !
    !                 regions and radial_densities (the latter we update).    !
    !                 All this happens in init_steric_pot().                  !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic in 07/2014.                                   !
    ! Updated by Jacek Dziedzic in 11/2014.                                   !
    ! Updated by James C. Womack in 02/2017.                                  !
    !=========================================================================!

    use bibliography, only: bibliography_cite
    use comms, only: pub_on_root
    use constants, only : FOUR_PI, k_B, VERBOSE, CRLF
    use model_type, only: MODEL
    use rundat, only: pub_is_pbe, pub_mg_pbe_use_fas, pub_is_pbe_exp_cap, &
         pub_output_detail, pub_is_pbe_temperature, pub_is_bulk_permittivity
    use utils, only: utils_abort, utils_feature_not_supported

    implicit none

    ! jd: Arguments
    type(MODEL), intent(inout)     :: mdl

    ! jd: Local variables
    real(kind=DP) :: exp_cap
    character(len=64) :: steric_pot_string
    character(len=*), parameter :: myself = 'implicit_solvent_boltzmann_init'

    !------------------------------------------------------------------------

#ifdef HAVE_DL_MG

    if(pub_is_pbe /= 'NONE') then

       call bibliography_cite('Electrolyte')

       if(pub_on_root .and. pub_output_detail >= VERBOSE) then
          write(stdout,'(a)') &
               'IS: Initialising Boltzmann terms in implicit solvation...'
       end if

       ! jd: Report on Boltzmann ions
       call boltzmann_ions_report

       ! jd: Now that we printed the report, switch the ions to internal
       !     convention, where electrons are positive.
       boltzmann_ions(:)%charge = -boltzmann_ions(:)%charge

       ! jd: Initialise steric potential
       call boltzmann_init_steric_pot(mdl)

       ! JCW: In the new DL_MG interface (01/2017), the Debye length is no
       ! JCW: longer provided by dl_mg_init_nonlin so we need to calculate this
       ! JCW: inside ONETEP. The following form is used:
       ! JCW:
       ! JCW:  Debye length =
       ! JCW:     [ ( k_B * T \epsilon_{r} ) / (4\pi \sum_{i} c_{i} q_{i}^{2} ) ]^{1/2}
       ! JCW:
       ! JCW: where vacuum permittivity = 1/(4\pi) in atomic units
       implicit_solvent_boltzmann_debye_length = &
            sqrt( k_B * pub_is_pbe_temperature * pub_is_bulk_permittivity / &
            ( FOUR_PI * sum ( boltzmann_ions(1:n_boltzmann_ions)%concentration_user * &
            boltzmann_ions(1:n_boltzmann_ions)%charge**2 ) ) )

       if(pub_on_root .and. pub_output_detail >= VERBOSE) then
          ! Determine exp_cap based on user input / DL_MG default
          if(pub_is_pbe_exp_cap == 0.0_DP) then
             exp_cap = max_expcap_dl_mg
          else
             exp_cap = pub_is_pbe_exp_cap
          end if
          select case(steric_pot_type)
          case('X')
             steric_pot_string = 'NONE'
          case('H')
             steric_pot_string = 'HARD-CORE SPHERES'
          case('M')
             steric_pot_string = 'SMOOTHED SPHERES'
          case('S')
             steric_pot_string = 'SOFT-CORE'
          case default
               call utils_abort(trim(myself)//': Unrecognized steric potential.')
          end select

          write(stdout,'(a,f7.2,a,a,a,l1,a,f7.2,a,a,a)') &
               'IS: Debye length is ',&
               implicit_solvent_boltzmann_debye_length, &
               ' bohr. PB equation variant: ', trim(pub_is_pbe), &
               '. Use_fas: ', pub_mg_pbe_use_fas, '.'//CRLF//'IS: Exp_cap: ', &
               exp_cap, '. Steric pot: ', trim(steric_pot_string), '.'
       end if

    else
       call utils_abort(trim(myself)//' must only be called &
            &when pub_is_pbe is not ''none''.')
    end if
#else
    call utils_feature_not_supported('DL_MG multigrid solver')
#endif

  end subroutine implicit_solvent_boltzmann_init
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_boltzmann_exit
    !=========================================================================!
    ! Cleans up after Boltzmann implicit solvent.                             !
    !=========================================================================!

    use rundat, only: pub_is_pbe
    use utils, only: utils_dealloc_check, utils_assert

    implicit none

    integer                     :: ierr
    character(len=*), parameter :: myself = 'implicit_solvent_boltzmann_exit'

    ! ------------------------------------------------------------------------

    call utils_assert(pub_is_pbe /= 'NONE',trim(myself)//': Unexpected call')

    call utils_assert(allocated(gamma), &
         'Internal error in '//trim(myself)//': accessibility not allocated')
    call utils_assert(allocated(steric_pot), &
         'Internal error in '//trim(myself)//': steric_pot not allocated')
    deallocate(gamma,stat=ierr)
    call utils_dealloc_check(myself,'gamma',ierr)
    deallocate(steric_pot,stat=ierr)
    call utils_dealloc_check(myself,'steric_pot',ierr)

  end subroutine implicit_solvent_boltzmann_exit
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_boltzmann_exit_ions
    !=========================================================================!
    ! Cleans up after Boltzmann implicit solvent (global state).              !
    !=========================================================================!

    use rundat, only: pub_is_pbe
    use utils, only: utils_dealloc_check, utils_assert

    implicit none

    integer                     :: ierr
    character(len=*), parameter :: myself = &
         'implicit_solvent_boltzmann_exit_ions'

    ! ------------------------------------------------------------------------

    call utils_assert(pub_is_pbe /= 'NONE',trim(myself)//': Unexpected call')

    call utils_assert(allocated(boltzmann_ions), &
         'Internal error in '//trim(myself)//': sol_ions not allocated')
    deallocate(boltzmann_ions, stat=ierr)
    call utils_dealloc_check(myself,'boltzmann_ions',ierr)

  end subroutine implicit_solvent_boltzmann_exit_ions
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_boltzmann_main(solvation_terms, &
       mg, grid, mg_return, phi, hartree_from_mg, pbc)
    !=========================================================================!
    ! Main subroutine for handling the Boltzmann part of implicit solvation.  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   solvation_terms (in/out): Individual solvation terms will be written  !
    !                             (updated) here.                             !
    !   mg (in): The grid on which MG operates.                               !
    !   grid (in): The ONETEP grid.                                           !
    !   mg_return (in): Details of the solution returned by DL_MG.            !
    !      %steric_weight_average (in): DL_MG's idea of accessibility fraction!
    !                                   (\Gamma). Only used to check against  !
    !                                   the value calculated in ONETEP.       !
    !      %betamu_elec (in): Boltzmann \beta * \mu_elec.                     !
    !      %used_ion_concentration (in): Bulk concentrations of Boltzmann ions!
    !                                    after having been potentially        !
    !                                    rescaled inside DL_MG. This happens  !
    !                                    in PBCs with the counterion models.  !
    !   phi (in): The electrostatic potential returned by DL_MG.              !
    !   hartree_from_mg (in): The electrostatic term ("molecular Hartree")    !
    !                         obtained in multigrid_methods.                  !
    !   pbc (in): .true. if working in PBCs, .false. for OBCs.                !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2019.                             !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root
    use constants, only: VERBOSE, DP, stdout
    use finite_differences, only: FD_GRID_INFO
    use is_poisson, only: SOLVATION_ENERGY_TERMS
    use rundat, only: pub_is_pbe, pub_output_detail, pub_is_pbe_energy_tolerance
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(SOLVATION_ENERGY_TERMS), intent(inout) :: solvation_terms
    type(FD_GRID_INFO), intent(in)              :: mg
    type(GRID_INFO), intent(in)                 :: grid
    real(kind=DP), intent(in)                   :: phi(&
         mg%ld1f, mg%ld2f, mg%max_slabs12)
    type(IS_PBE_MG_RETURN), intent(in)          :: mg_return
    real(kind=DP), intent(in)                   :: hartree_from_mg
    logical, intent(in)                         :: pbc

    ! jd: Local variables
    real(kind=DP) :: E_half_rho_f_phi
    real(kind=DP) :: Etot_expr1
    real(kind=DP) :: Etot_expr2
    real(kind=DP) :: big_gamma
    character(len=*), parameter :: myself = 'implicit_solvent_boltzmann_main'

    ! -------------------------------------------------------------------------

    call utils_assert(pub_is_pbe == "FULL" .or. pub_is_pbe == "LINEARISED", &
         myself//': Illegal pub_is_pbe.')

    Etot_expr1 = 0.0_DP
    Etot_expr2 = 0.0_DP

    ! jd: Overwrite bulk concentrations with current, possibly scaled, concent-
    !     rations from DL_MG. Currently they are only scaled in the counterion
    !     model. In other cases DL_MG returns the same concentrations that we
    !     passed to it, making this statement a no-op.
    boltzmann_ions(:)%concentration_adj = &
         mg_return%used_ion_concentration(1:n_boltzmann_ions)

    if(pub_is_pbe == 'FULL' .or. pub_is_pbe == 'LINEARISED') then

       E_half_rho_f_phi = hartree_from_mg

       big_gamma = boltzmann_calc_big_gamma(mg, mg_return%steric_weight_average)

       call implicit_solvent_boltzmann_free_energy_terms(phi, mg, &
            big_gamma, mg_return%betamu_elec(1:n_boltzmann_ions), &
            solvation_terms)

       Etot_expr1 = E_half_rho_f_phi + &
            implicit_solvent_boltzmann_energy_contrib(solvation_terms)

       if(pub_on_root  .and. pub_output_detail >= VERBOSE) then
          write(stdout,'(a,f14.6)') "MG Fixed charge (1/2 \int rho_f phi dr) &
               &energy term:             ", E_half_rho_f_phi
          write(stdout,'(a,f14.6)') "MG Mobile charge (1/2 \int rho_m phi dr)&
               & energy term:            ", solvation_terms%E_elec_mob
          write(stdout,'(a,f14.6)') "MG Osmotic pressure (-kT \int dr \sum_i &
               &c_i(r)) energy term:     ", solvation_terms%E_osmo
          write(stdout,'(a,f14.6)') "MG Mobile ion acc. corr. (\int dr \sum_i&
               & c_i(r) * V_steric(r)):  ", solvation_terms%E_acc
          write(stdout,'(a,f14.6)') "MG Ion atmosph. rearr. (kT \int dr \sum_&
               &i c_i(r)*ln(c_i(r)/c^0): ", solvation_terms%E_atmo
          write(stdout,'(a,f14.6)') "MG Chemical potential (-\sum_i \mu_i \in&
               &t dr c_i(r)) energy term:", solvation_terms%E_chempot
          write(stdout,'(a)') "MG                                              &
               &                 --------------"
          write(stdout,'(a,f14.6)') "MG Mean-field contribution to grand pote&
               &ntial                    ", Etot_expr1
       end if

       call implicit_solvent_boltzmann_free_energy_terms_check(phi, mg, &
            big_gamma, mg_return%betamu_elec(1:n_boltzmann_ions), &
            solvation_terms)

       if(pbc) then
          Etot_expr2 = E_half_rho_f_phi + &
               solvation_terms%E_elec_minus_mob + &
               solvation_terms%E_minus_ktv_gamma_c_i_bulk
       else
          Etot_expr2 = E_half_rho_f_phi + &
               solvation_terms%E_elec_minus_mob + &
               solvation_terms%E_minus_kt_int_c_i
       end if

       if(pub_on_root  .and. pub_output_detail >= VERBOSE) then
          write(stdout,'(/a)') "MG Check from alternate expression:"
          write(stdout,'(a,f14.6)') "MG Fixed charge (1/2 \int rho_f phi dr) &
               &energy term:             ", E_half_rho_f_phi
          write(stdout,'(a,f14.6)') "MG Minus mobile charge (-1/2 \int rho_m &
               &phi dr) energy term:     ", solvation_terms%E_elec_minus_mob
          if(pbc) then
             write(stdout,'(a,f14.6)') "MG -kTV\Gamma \sum_i c^bulk_i(r) (osmot&
                  &ic pressure) energy term: ", &
                  solvation_terms%E_minus_ktv_gamma_c_i_bulk
          else
             write(stdout,'(a,f14.6)') "MG -kT \int \sum_i c_i(r) dr (osmotic&
                  & pressure) energy term:     ", &
                  solvation_terms%E_minus_kt_int_c_i
          end if
          write(stdout,'(a)') "MG                                              &
               &                 --------------"
          write(stdout,'(a,f14.6)') "MG Mean-field contribution to grand &
               &potential:                   ", Etot_expr2
       end if

       call utils_assert(abs(Etot_expr1 - Etot_expr2) < &
            pub_is_pbe_energy_tolerance, myself//': Discrepancy between two &
            &energy expressions that should be identical is too large. This &
            &could mean your system is not well-representable on the fine grid.&
            & Particularly, the steric potential might be to steep or, if you&
            & are using the linearised approximation, it might have broken&
            & down. The values of the two expressions, and the tolerance &
            &(set with "is_pbe_energy_tolerance") (all in Ha) follow.', &
            Etot_expr1, Etot_expr2, pub_is_pbe_energy_tolerance)

    end if

    call implicit_solvent_boltzmann_ion_stats(phi, mg, grid, mg_return, &
         big_gamma)

  end subroutine implicit_solvent_boltzmann_main
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_boltzmann_defect(bt, pot, n1, n2, n3)
    !=========================================================================!
    ! Computes the Boltzmann term to defect correction.                       !
    !                                                                         !
    ! Updated by James C. Womack, 02/2017, to fix mistake in kappa2 for       !
    ! linearized PB equation defect.                                          !
    !=========================================================================!

    use constants, only : DP, PI, k_B
    use rundat, only: pub_is_bulk_permittivity, &
         pub_is_pbe, pub_is_pbe_temperature, pub_is_pbe_exp_cap
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: bt(:,:,:)
    real(kind=DP), intent(in)  :: pot(:,:,:)
    integer, intent(in)        :: n1, n2, n3

    ! jd: Internal variables
    integer                    :: i
    real(kind=DP)              :: q, beta, kappa2
    real(kind=DP), parameter   :: FOURPI = 4.0_DP * PI
    real(kind=DP)              :: exp_cap

    ! -----------------------------------------------------------------------

    beta = 1.0_DP/(k_B*pub_is_pbe_temperature)

    if(pub_is_pbe_exp_cap == 0.0_DP) then
       exp_cap = max_expcap_dl_mg
    else
       exp_cap = pub_is_pbe_exp_cap
    end if

    bt(:, :, :) = 0.0_DP
    if (pub_is_pbe == "FULL") then
       do i = 1, size(boltzmann_ions)
          q = boltzmann_ions(i)%charge
          bt(1:n1, 1:n2, 1:n3) = bt(1:n1, 1:n2, 1:n3) + &
               FOURPI*boltzmann_ions(i)%concentration_adj * q * &
               gamma(1:n1, 1:n2, 1:n3) * &
               exp(min(exp_cap, -q * beta * pot(1:n1, 1:n2, 1:n3)))
       enddo
    elseif (pub_is_pbe == "LINEARISED") then
       kappa2 = pub_is_bulk_permittivity/implicit_solvent_boltzmann_debye_length**2
       bt(1:n1, 1:n2, 1:n3) = - kappa2 * &
            gamma(1:n1, 1:n2, 1:n3) * pot(1:n1, 1:n2, 1:n3)
    else
       call utils_abort("Internal error in implicit_solvent_boltzmann_defect")
    end if

  end subroutine implicit_solvent_boltzmann_defect
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_boltzmann_free_energy_terms(pot, mg_grid, &
       big_gamma, betamu_elec, solvation_terms)
    !=========================================================================!
    ! Calculates all free energy terms relevant in Boltzmann solvation:       !
    ! - mobile ion energy term: 1/2 * \int dr sum_i z_i c_i(r) \phi(r).       !
    ! - osmotic pressure term, -kT \int dr \sum_i c_i(r).                     !
    ! - accessibility correction term, \int dr \sum_i c_i(r) * V_steric(r).   !
    ! - ion atmos. entropic term, kT \int dr \sum_i c_i(r) * ln(c_i(r) / c^0).!
    ! - the chemical potential energy term, -\sum_i \mu_i \int dr c_i(r).     !
    ! - total energy of pure electrolyte, -kTV \sum_i c_i_bulk.               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   pot (in): Electrostatic potential from the solution of the PBE.       !
    !   mg_grid (in): The grid on which this happens.                         !
    !   big_gamma (in): Accessible fraction (\Gamma), calculated in ONETEP.   !

    !   betamu_elec (in): Boltzmann beta * \mu_elec, from DL_MG.              !
    !   solvation_terms (inout): Results are returned here.                   !
    !-------------------------------------------------------------------------!
    ! Note: The remaining terms in solvation_terms are left untouched.        !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2019.                             !
    !=========================================================================!

    use comms, only: comms_reduce
    use constants, only: k_B, CONC_0, CONC_MOL_PER_L_TO_AU, CONC_AU_TO_MOL_PER_L
    use finite_differences, only: FD_GRID_INFO
    use is_poisson, only: SOLVATION_ENERGY_TERMS
    use rundat, only: pub_is_pbe_temperature, pub_is_pbe_exp_cap, pub_is_pbe, &
         pub_multigrid_bc_is_periodic
    use utils, only: utils_assert, utils_abort

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in)      :: pot(:,:,:)
    type(FD_GRID_INFO), intent(in) :: mg_grid
    real(kind=DP), intent(in)      :: big_gamma
    real(kind=DP), intent(in)      :: betamu_elec(1:n_boltzmann_ions)
    type(SOLVATION_ENERGY_TERMS), intent(inout) :: solvation_terms

    ! jd: Local variables
    integer       :: i1, i2, islab12
    integer       :: b_ion
    integer       :: tot_n_points
    real(kind=DP) :: beta
    real(kind=DP) :: pot_here, gamma_here, steric_pot_here
    real(kind=DP) :: z_i
    real(kind=DP) :: c_i, c_i_bulk, c_i_bulk_undiluted
    real(kind=DP) :: exp_cap
    real(kind=DP) :: mu_ideal(1:n_boltzmann_ions)
    real(kind=DP) :: mu_excess(1:n_boltzmann_ions)
    real(kind=DP) :: mu_total(1:n_boltzmann_ions)
    real(kind=DP) :: e_elec_mob, e_osmo, e_acc, e_atmo, e_chempot, e_pure
    logical       :: pbc
    logical       :: is_full
    real(kind=DP) :: negative_threshold = -0.001_DP * CONC_MOL_PER_L_TO_AU
    !                                   ^ -0.001 mol/L
    character(len=*), parameter :: myself = &
         'implicit_solvent_boltzmann_free_energy_terms'

    ! -------------------------------------------------------------------------

    if(any(pub_multigrid_bc_is_periodic)) then
       call utils_assert(all(pub_multigrid_bc_is_periodic), &
            myself//': Mixed boundary conditions are not supported.')
       pbc = .true.
    else
       pbc = .false.
    end if

    if(pub_is_pbe_exp_cap == 0.0_DP) then
       exp_cap = max_expcap_dl_mg
    else
       exp_cap = pub_is_pbe_exp_cap
    end if

    beta = 1.0_DP/(k_B*pub_is_pbe_temperature)

    tot_n_points = mg_grid%pq1f * mg_grid%pq2f * mg_grid%pq3f

    if(pub_is_pbe == "FULL") then
       is_full = .true.
    else if(pub_is_pbe == "LINEARISED") then
       is_full = .false.
    else
       call utils_abort(myself//': Illegal "is_pbe": '//trim(pub_is_pbe))
    end if

    call implicit_solvent_boltzmann_chem_pot(mu_ideal, mu_excess, mu_total, &
         betamu_elec, beta, big_gamma)

    e_elec_mob = 0.0_DP
    e_osmo = 0.0_DP
    e_acc = 0.0_DP
    e_atmo = 0.0_DP
    e_chempot = 0.0_DP
    e_pure = 0.0_DP

    do islab12 = 1, mg_grid%num_slabs12_pq
       do i2 = 1, mg_grid%pq2f
          do i1 = 1, mg_grid%pq1f
             pot_here = pot(i1,i2,islab12)
             gamma_here = gamma(i1,i2,islab12)
             steric_pot_here = steric_pot(i1,i2,islab12)
             do b_ion = 1, n_boltzmann_ions
                z_i = boltzmann_ions(b_ion)%charge
                c_i_bulk = boltzmann_ions(b_ion)%concentration_adj
                c_i = implicit_solvent_boltzmann_c_i(z_i, pot_here, c_i_bulk, &
                     beta, gamma_here, betamu_elec(b_ion), big_gamma, exp_cap, &
                     pbc, is_full)
                e_osmo = e_osmo + c_i
                e_elec_mob = e_elec_mob + c_i * z_i * pot_here
                e_acc = e_acc + c_i * steric_pot_here
                e_chempot = e_chempot + mu_total(b_ion) * c_i
                e_pure = e_pure + c_i_bulk ! Note: no longer contains Gamma

                ! jd: Ignore tiny numerical noise in c_i that can make it zero
                !     or negligibly negative. Bail out when seriously negative
                !     c_i is detected.
                if(c_i > 0.0_DP) then
                   e_atmo = e_atmo + c_i * log(c_i/CONC_0)
                else if(c_i < negative_threshold) then
                   call utils_abort(myself//&
                        ': Unphysical negative concentration (mol/L): ', &
                        opt_real_to_print1 = c_i * CONC_AU_TO_MOL_PER_L)
                end if

             end do
          end do
       end do
    end do

    call comms_reduce('SUM',e_elec_mob)
    call comms_reduce('SUM',e_osmo)
    call comms_reduce('SUM',e_acc)
    call comms_reduce('SUM',e_atmo)
    call comms_reduce('SUM',e_chempot)
    call comms_reduce('SUM',e_pure)

    ! jd: Include constant factors and weights
    solvation_terms%E_elec_mob = 0.5_DP * e_elec_mob * mg_grid%weight
    solvation_terms%E_osmo = -e_osmo / beta * mg_grid%weight
    solvation_terms%E_acc = e_acc * mg_grid%weight
    solvation_terms%E_atmo = e_atmo / beta * mg_grid%weight
    solvation_terms%E_chempot = -e_chempot * mg_grid%weight
    solvation_terms%E_pure = -e_pure / beta * mg_grid%weight

  end subroutine implicit_solvent_boltzmann_free_energy_terms
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_boltzmann_free_energy_terms_check(pot, mg_grid, &
       big_gamma, betamu_elec, solvation_terms)
    !=========================================================================!
    ! Calculates all free energy terms relevant in Boltzmann solvation, from  !
    ! a different expression, to serve as a check.                            !
    ! - minus mobile ion energy term: -1/2 * \int dr sum_i z_i c_i(r) \phi(r).!
    ! - -kTV \Gamma \sum_i c_i^bulk (only in PBC),                            !
    ! - -kT \int dr sum_i c_i(r) (only in OBC).                               !
    !-------------------------------------------------------------------------!
    ! The expressions implemented are Ref. 3, eq. 31 (for PBC), and eq. 35    !
    ! (for OBC), *EXCEPT* for the fixed electrostatic term (first term in     !
    ! brackets in eq. 10), which is omitted here, and is added in             !
    ! implicit_solvent_boltzmann_main() as E_half_rho_f_phi, and is known     !
    ! elsewhere in the code as the "molecular Hartree" term.                  !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2019.                             !
    !=========================================================================!

    use comms, only: comms_reduce
    use constants, only: k_B
    use finite_differences, only: FD_GRID_INFO
    use is_poisson, only: SOLVATION_ENERGY_TERMS
    use rundat, only: pub_is_pbe_temperature, pub_is_pbe_exp_cap, pub_is_pbe, &
         pub_multigrid_bc_is_periodic
    use utils, only: utils_assert, utils_abort

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in)      :: pot(:,:,:)
    type(FD_GRID_INFO), intent(in) :: mg_grid
    real(kind=DP), intent(in)      :: big_gamma
    real(kind=DP), intent(in)      :: betamu_elec(1:n_boltzmann_ions)
    type(SOLVATION_ENERGY_TERMS), intent(inout) :: solvation_terms

    ! jd: Local variables
    integer       :: i1, i2, islab12
    integer       :: b_ion
    integer       :: tot_n_points
    real(kind=DP) :: beta
    real(kind=DP) :: pot_here, gamma_here
    real(kind=DP) :: z_i
    real(kind=DP) :: c_i, c_i_bulk
    real(kind=DP) :: exp_cap
    real(kind=DP) :: box_volume
    real(kind=DP) :: e_elec_minus_mob
    real(kind=DP) :: e_minus_ktv_gamma_c_bulk
    real(kind=DP) :: e_minus_kt_int_c_i
    logical       :: pbc
    logical       :: is_full
    character(len=*), parameter :: myself = &
         'implicit_solvent_boltzmann_free_energy_terms_check'

    ! -------------------------------------------------------------------------

    if(any(pub_multigrid_bc_is_periodic)) then
       call utils_assert(all(pub_multigrid_bc_is_periodic), &
            myself//': Mixed boundary conditions are not supported.')
       pbc = .true.
    else
       pbc = .false.
    end if

    if(pub_is_pbe_exp_cap == 0.0_DP) then
       exp_cap = max_expcap_dl_mg
    else
       exp_cap = pub_is_pbe_exp_cap
    end if

    beta = 1.0_DP/(k_B*pub_is_pbe_temperature)

    tot_n_points = mg_grid%pq1f * mg_grid%pq2f * mg_grid%pq3f
    box_volume = real(tot_n_points,kind=DP) * mg_grid%weight

    if(pub_is_pbe == "FULL") then
       is_full = .true.
    else if(pub_is_pbe == "LINEARISED") then
       is_full = .false.
    else
       call utils_abort(myself//': Illegal "is_pbe": '//trim(pub_is_pbe))
    end if

    e_elec_minus_mob = 0.0_DP
    e_minus_ktv_gamma_c_bulk = 0.0_DP
    e_minus_kt_int_c_i = 0.0_DP

    do islab12 = 1, mg_grid%num_slabs12_pq
       do i2 = 1, mg_grid%pq2f
          do i1 = 1, mg_grid%pq1f
             pot_here = pot(i1,i2,islab12)
             gamma_here = gamma(i1,i2,islab12)
             do b_ion = 1, n_boltzmann_ions
                z_i = boltzmann_ions(b_ion)%charge
                c_i_bulk = boltzmann_ions(b_ion)%concentration_adj
                c_i = implicit_solvent_boltzmann_c_i(z_i, pot_here, c_i_bulk, &
                     beta, gamma_here, betamu_elec(b_ion), big_gamma, exp_cap, &
                     pbc, is_full)
                e_elec_minus_mob = e_elec_minus_mob + c_i * z_i * pot_here
                e_minus_kt_int_c_i = e_minus_kt_int_c_i + c_i
             end do
          end do
       end do
    end do

    call comms_reduce('SUM',e_elec_minus_mob)
    call comms_reduce('SUM',e_minus_kt_int_c_i)

    ! jd: Include constant factors and weights
    solvation_terms%E_elec_minus_mob = -0.5_DP * e_elec_minus_mob * mg_grid%weight
    if(pbc) then
       solvation_terms%E_minus_kt_int_c_i = 0.0_DP
    else
       solvation_terms%E_minus_kt_int_c_i = -e_minus_kt_int_c_i / beta * mg_grid%weight
    end if

    if(pbc) then
       do b_ion = 1, n_boltzmann_ions
          c_i_bulk = boltzmann_ions(b_ion)%concentration_adj
          e_minus_ktv_gamma_c_bulk = e_minus_ktv_gamma_c_bulk - &
               c_i_bulk / beta * big_gamma * box_volume
       end do
       solvation_terms%E_minus_ktv_gamma_c_i_bulk = e_minus_ktv_gamma_c_bulk
    else
       solvation_terms%E_minus_ktv_gamma_c_i_bulk = 0.0_DP
    end if

  end subroutine implicit_solvent_boltzmann_free_energy_terms_check
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  real(kind=DP) function implicit_solvent_boltzmann_energy_contrib(&
       solvation_terms)
    !=========================================================================!
    ! Calculates the total Boltzmann contribution to the energy, that is all  !
    ! the terms except E_elec_f aka molecular Hartree. This function exists   !
    ! for the convenience of callers -- we'd rather not burden them with the  !
    ! details of *which* energy terms need to be summed depending on BCs or   !
    ! whether we do FULL or LINEARISED Boltzmann.                             !
    !-------------------------------------------------------------------------!
    ! The expression implemented is Ref. 3, eq. 5, where individual terms are !
    ! calculated according to eqs. 10, 12, 13, 15, and 16 *EXCEPT* for the    !
    ! fixed electrostatic term (first term in brackets in eq. 10), which is   !
    ! omitted here, and is added in implicit_solvent_boltzmann_main() as      !
    ! E_half_rho_f_phi, and is known elsewhere in the code as the "molecular  !
    ! Hartree" term.                                                          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   solvation_terms (in): The datastructure containing energy terms pre-  !
    !                         viously calculated via isb_free_energy_terms(). !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2019.                             !
    !=========================================================================!

    use is_poisson, only: SOLVATION_ENERGY_TERMS
    use rundat, only: pub_is_pbe

    implicit none

    ! jd: Arguments
    type(SOLVATION_ENERGY_TERMS), intent(in) :: solvation_terms

    ! jd: Local variables
    character(len=*), parameter :: myself = &
         'implicit_solvent_boltzmann_energy_contrib'

    ! -------------------------------------------------------------------------

    ! jd: Ultimately, we found that the energy expression is identical
    !     between FULL and LINEARISED.
    !     BCs only matter in the energy check performed in multigrid_methods.
    if(pub_is_pbe == "FULL") then
       implicit_solvent_boltzmann_energy_contrib = &
            solvation_terms%E_elec_mob + &
            solvation_terms%E_osmo + &
            solvation_terms%E_acc + &
            solvation_terms%E_atmo + &
            solvation_terms%E_chempot
    else if(pub_is_pbe == "LINEARISED") then
       implicit_solvent_boltzmann_energy_contrib = &
            solvation_terms%E_elec_mob + &
            solvation_terms%E_osmo + &
            solvation_terms%E_acc + &
            solvation_terms%E_atmo + &
            solvation_terms%E_chempot
    else if(pub_is_pbe == "NONE") then
       implicit_solvent_boltzmann_energy_contrib = 0.0_DP
    end if

  end function implicit_solvent_boltzmann_energy_contrib
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  real(kind=DP) function implicit_solvent_boltzmann_shifted_hartree(&
       grid, mg, mg_return, rho, phi)
    !=========================================================================!
    ! Calculates the 'molecular Hartree' energy term in PBC solvation. This   !
    ! is also called 'electrostatic contribution' in some places.             !
    ! Normally this term is calculated in multigrid_methods, as               !
    ! 1/2 * \int \rho(r) \phi(r) d3r, but in PBCs and simultaneously in the   !
    ! absence of Boltzmann counterions we need to shift \rho so that <rho>=0. !
    ! This can be done in more than one way, depending on how we neutralise.  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   grid (in): The ONETEP grid.                                           !
    !   mg (in): The grid on which MG operates.                               !
    !   mg_return (in): Details of the solution returned by DL_MG.            !
    !      %steric_weight_average (in): DL_MG's idea of accessibility fraction!
    !                                   (\Gamma). Only used to check against  !
    !                                   the value calculated in ONETEP.       !
    !   rho (in): The original (unshifted) density.                           !
    !   phi (in): Electrostatic potential from DL_MG.                         !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in April 2020.                                !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use finite_differences, only: FD_GRID_INFO
    use integrals, only: integrals_product_on_grid
    use rundat, only: pub_is_pbe_neutralisation_scheme
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! jd: Arguments
    type(GRID_INFO), intent(in)        :: grid
    type(FD_GRID_INFO), intent(in)     :: mg
    type(IS_PBE_MG_RETURN), intent(in) :: mg_return
    real(kind=DP), intent(in)          :: rho(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in)          :: phi(mg%ld1f,mg%ld2f,mg%max_slabs12)

    ! jd: Local variables
    real(kind=DP), allocatable :: rho_shifted(:,:,:)
    real(kind=DP) :: rho_avg
    real(kind=DP) :: big_gamma
    real(kind=DP) :: v_acc
    integer :: n_points
    integer :: ierr
    character(len=*), parameter :: myself = &
         'implicit_solvent_boltzmann_shifted_hartree'

    ! -------------------------------------------------------------------------

    n_points = mg%pq1f * mg%pq2f * mg%pq3f
    rho_avg = integrals_product_on_grid(grid, rho, &
         m1=mg%pq1f, m2=mg%pq2f, m3=mg%num_slabs12_pq) / real(n_points, kind=DP)

    allocate(rho_shifted(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
    call utils_alloc_check(myself,'rho_shifted',ierr)

    select case(pub_is_pbe_neutralisation_scheme)
    case('JELLIUM')
       ! jd: For jellium we shift according to [4, eq. (5)].
       rho_shifted(:,:,:) = rho(:,:,:) - rho_avg / grid%weight
    case('ACCESSIBLE_JELLIUM')
       ! jd: For accessible jellium we follow [4, eq. (6)].
       big_gamma = boltzmann_calc_big_gamma(mg, mg_return%steric_weight_average)
       v_acc = real(n_points,kind=DP) * grid%weight * big_gamma
       rho_avg = rho_avg / big_gamma
       rho_shifted(:,:,:) = rho(:,:,:) - gamma(:,:,:) * rho_avg / grid%weight
    case default
       ! jd: For counterion schemes we should never have been called.
       call utils_abort(myself//": 'is_pbe_neutralisation_scheme' unrecognised &
            &or illegal in this context: '"//&
            trim(pub_is_pbe_neutralisation_scheme)//"'.")
    end select

    implicit_solvent_boltzmann_shifted_hartree = &
         0.5_DP * integrals_product_on_grid(grid, &
         phi, rho_shifted, mg%pq1f, mg%pq2f, mg%num_slabs12_pq)

    deallocate(rho_shifted,stat=ierr)
    call utils_dealloc_check(myself,'rho_shifted',ierr)

  end function implicit_solvent_boltzmann_shifted_hartree
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  real(kind=DP) function implicit_solvent_boltzmann_c_i(z_i, pot_here, &
       c_i_bulk, beta, gamma_here, betamu_elec, big_gamma, exp_cap, pbc, &
       is_full)
    !=========================================================================!
    ! Calculates the concentration of a Boltzmann ion at a point.             !
    ! Formulas differ depending on BCs and whether we are following the FULL  !
    ! or LINEARISED approach.                                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   z_i (in): The charge on Boltzmann ion i.                              !
    !   pot_here (in): The value of the potential from PBE at this point.     !
    !   c_i_bulk (in): Bulk concentration of Boltzmann ion i.                 !
    !   beta (in): 1/kT.                                                      !
    !   gamma_here (in): The value of electrolyte accessibility at this point.!
    !   betamu_elec (in): 1/kT * \mu_elec, obtained from DL_MG.               !
    !   big_gamma (in): Accessible fraction of volume, calculated in ONETEP.  !
    !   exp_cap (in): Capping value for the exponentials.                     !
    !   pbc (in): .true. in PBCs, .false. in OBCs.                            !
    !   is_full (in): .true. if doing the FULL eq'n, .false. if LINEARISED.   !
    !                                                                         !
    !   The different expressions for c_i reference different sets of the     !
    !   above arguments, e.g. exp_cap is not needed in the LINEARISED case.   !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2019.                             !
    !=========================================================================!

    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in) :: z_i
    real(kind=DP), intent(in) :: pot_here
    real(kind=DP), intent(in) :: c_i_bulk
    real(kind=DP), intent(in) :: beta
    real(kind=DP), intent(in) :: gamma_here
    real(kind=DP), intent(in) :: betamu_elec
    real(kind=DP), intent(in) :: big_gamma
    real(kind=DP), intent(in) :: exp_cap
    logical, intent(in)       :: pbc
    logical, intent(in)       :: is_full

    ! jd: Local variables
    real(kind=DP) :: exp_term
    character(len=*), parameter :: myself = 'implicit_solvent_boltzmann_c_i'

    ! -------------------------------------------------------------------------

    if(gamma_here < 0.0_DP .or. gamma_here > 1.000001_DP) then
       call utils_abort(myself//': Illegal accessibility: ', &
            opt_real_to_print1 = gamma_here)
    end if

    if(is_full) then
       if(pbc) then
          exp_term = &
               exp(min(exp_cap,-z_i * pot_here * beta + betamu_elec))
          implicit_solvent_boltzmann_c_i = c_i_bulk * exp_term * gamma_here
       else
          exp_term = exp(min(exp_cap,-z_i * pot_here * beta))
          implicit_solvent_boltzmann_c_i = c_i_bulk * exp_term * gamma_here
       end if
    else
       if(pbc) then
          ! jd: We need an exp version even in the linearised case, since
          !     otherwise we risk negative concentrations, which has problematic
          !     together with the c*ln(c) term.
          implicit_solvent_boltzmann_c_i = c_i_bulk * gamma_here * &
               exp(min(exp_cap,- beta * z_i * pot_here))
       else
          ! jd: We need an exp version even in the linearised case, since
          !     otherwise we risk negative concentrations, which has problematic
          !     together with the c*ln(c) term.
          implicit_solvent_boltzmann_c_i = c_i_bulk * gamma_here * &
               exp(min(exp_cap,- beta * z_i * pot_here))
       end if
    end if

  end function implicit_solvent_boltzmann_c_i
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_boltzmann_ion_conc(bion_conc, &
       big_gamma, phi, betamu_elec, mg_grid, opt_bion_conc_avg)
    !=========================================================================!
    ! Calculates the Boltzmann ion concentrations (cf. [1], eqs. (8), (12))   !
    ! from the electrostatic potential and accessibility (gamma). The         !
    ! potential is specified via an argument, gamma is a module-wide variable.!
    ! Each kind (species) of Boltzmann ion is done separately (last index of  !
    ! bion_conc).                                                             !
    !                                                                         !
    ! Usual slab distribution (over index 3) is used for bion_conc.           !
    ! The optionals: average and accessible fraction -- are automatically     !
    ! averaged across procs and returned on all procs.                        !
    ! This subroutine is MPI-collective (potentially, when opt'nals are used).!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   bion_conc (output): The Boltzmann ion ("bion") concentrations are     !
    !                       returned here. Internal units of concentration    !
    !                       are used.                                         !
    !   big_gamma (input): Accessibility fraction, calculated by ONETEP.      !
    !   phi (input): The electrostatic potential from the solution of the PB  !
    !                equation.                                                !
    !   betamu_elec (input): Boltzmann \beta * \mu_elec, calculated in DL_MG. !
    !   mg_grid (input): The grid on which this happens.                      !
    !   opt_bion_conc_avg (output, opt): If specified, average concentrations !
    !                                    over the entire box (up to pq) are   !
    !                                    returned here for every bion.        !
    !-------------------------------------------------------------------------!
    ! Caveats:                                                                !
    !   It is the responsibility of the caller to supply a suitably allocated !
    !   output buffer for bion_conc.                                          !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2019.                                 !
    !=========================================================================!

    use comms, only: comms_reduce
    use constants, only : DP, k_B
    use finite_differences, only: FD_GRID_INFO
    use rundat, only: pub_is_pbe, pub_is_pbe_temperature, pub_is_pbe_exp_cap, &
         pub_multigrid_bc_is_periodic
    use utils, only: utils_assert, utils_abort, utils_real_to_str

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)           :: bion_conc(:,:,:,:) ! last idx is ion #
    real(kind=DP), intent(in)            :: big_gamma
    real(kind=DP), intent(in)            :: phi(:,:,:)
    real(kind=DP), intent(in)            :: betamu_elec(1:n_boltzmann_ions)
    type(FD_GRID_INFO), intent(in)       :: mg_grid
    real(kind=DP), intent(out), optional :: opt_bion_conc_avg(:)

    ! jd: Local variables
    integer :: i1, i2, islab12
    integer :: b_ion
    integer :: d
    integer :: tot_n_points ! over *all* MPI ranks
    real(kind=DP) :: beta
    real(kind=DP) :: z_i
    real(kind=DP) :: c_i_bulk
    real(kind=DP) :: lin_term
    real(kind=DP) :: exp_cap
    real(kind=DP) :: phi_here, gamma_here
    real(kind=DP) :: my_sum, tot_sum
    real(kind=DP) :: betamu_el(1:n_boltzmann_ions)
    logical       :: pbc
    logical       :: is_full
    character(len=*), parameter :: myself = 'implicit_solvent_boltzmann_ion_conc'

    ! -------------------------------------------------------------------------

    do d = 1, 3
       call utils_assert(lbound(bion_conc,d) == lbound(phi,d), &
            myself//': Incorrectly dimensioned input (lbound).')
       call utils_assert(ubound(bion_conc,d) == ubound(phi,d), &
            myself//': Incorrectly dimensioned input (ubound).')
    end do

    call utils_assert(mg_grid%pq1f > 0, myself//': Badly dimensioned MG grid.')
    call utils_assert(mg_grid%pq2f > 0, myself//': Badly dimensioned MG grid.')
    call utils_assert(mg_grid%num_slabs12_pq > 0, &
         myself//': Badly dimensioned MG grid.')

    ! -------------------------------------------------------------------------

    beta = 1.0_DP/(k_B*pub_is_pbe_temperature)

    if(pub_is_pbe_exp_cap == 0.0_DP) then
       exp_cap = max_expcap_dl_mg
    else
       exp_cap = pub_is_pbe_exp_cap
    end if

    bion_conc(:,:,:,:) = 0.0_DP ! jd: takes care of padding

    if(any(pub_multigrid_bc_is_periodic)) then
       call utils_assert(all(pub_multigrid_bc_is_periodic), &
            myself//': Mixed boundary conditions are not supported.')
       pbc = .true.
    else
       pbc = .false.
    end if

    if(pub_is_pbe == "FULL") then
       is_full = .true.
    else if(pub_is_pbe == "LINEARISED") then
       is_full = .false.
    else
       call utils_abort(myself//': Illegal "is_pbe": '//trim(pub_is_pbe))
    end if

    if(pbc) then
       betamu_el(1:n_boltzmann_ions) = betamu_elec(1:n_boltzmann_ions)
    else
       betamu_el(1:n_boltzmann_ions) = 0.0_DP
    end if

    if(pub_is_pbe == "FULL" .or. pub_is_pbe == "LINEARISED") then
       do islab12 = 1, mg_grid%num_slabs12_pq
          do i2 = 1, mg_grid%pq2f
             do i1 = 1, mg_grid%pq1f
                phi_here = phi(i1,i2,islab12)
                gamma_here = gamma(i1,i2,islab12)
                if(gamma_here < 0.0_DP .or. gamma_here > 1.000001_DP) then
                   call utils_abort(myself//': Illegal accessibility: '//&
                        trim(utils_real_to_str(gamma_here)), i1, i2, islab12)
                end if
                do b_ion = 1, n_boltzmann_ions
                   z_i = boltzmann_ions(b_ion)%charge
                   c_i_bulk = boltzmann_ions(b_ion)%concentration_adj
                   bion_conc(i1,i2,islab12,b_ion) = &
                        implicit_solvent_boltzmann_c_i(z_i, phi_here, &
                        c_i_bulk, beta, gamma_here, betamu_elec(b_ion), &
                        big_gamma, exp_cap, pbc, is_full)
                end do
             end do
          end do
       end do
    else
       call utils_abort(myself//': unrecognised is_pbe: "'//&
            trim(pub_is_pbe)//'"')
    end if

    tot_n_points = mg_grid%pq1f * mg_grid%pq2f * mg_grid%pq3f

    if(present(opt_bion_conc_avg)) then
       do b_ion = 1, n_boltzmann_ions
          my_sum = sum(bion_conc(:,:,:,b_ion))
          call comms_reduce('SUM',my_sum)
          tot_sum = my_sum
          opt_bion_conc_avg(b_ion) = tot_sum / tot_n_points
       end do
    end if

  end subroutine implicit_solvent_boltzmann_ion_conc
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_boltzmann_ion_dens(bion_dens, phi, betamu_elec, &
       big_gamma, mg_grid, opt_bion_charge_tot)
    !=========================================================================!
    ! Calculates the Boltzmann ion charge density, ie. z_i c_i, where c_i is  !
    ! obtained from implicit_solvent_boltzmann_ion_conc().                    !
    ! Inputs are the electrostatic potential and accessibility (gamma). The   !
    ! potential is specified via an argument, gamma is a module-wide variable.!
    ! Each kind (species) of Boltzmann ion is done separately (last index of  !
    ! bion_dens).                                                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   bion_dens (output): The Boltzmann ion ("bion") charge densities are   !
    !                       returned here. Internal (atomic) units of density !
    !                       are used. Internal convention (electrons positive)!
    !                       is used for charge.                               !
    !   phi (input): The electrostatic potential from the solution of the PB  !
    !                equation.                                                !
    !   betamu_elec (input): Boltzmann \beta * \mu_elec, returned by DL_MG.   !
    !   big_gamma (input): Accessible fraction of vol., calculated in ONETEP. !
    !   mg_grid (input): The grid on which this happens.                      !
    !   opt_bion_charge_tot (output, opt): If specified, total charge over    !
    !                                      the entire box (up to pq) is       !
    !                                      returned here for every bion.      !
    !-------------------------------------------------------------------------!
    ! Caveats:                                                                !
    !   It is the responsibility of the caller to supply a suitably allocated !
    !   output buffer for bion_dens.                                          !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2019.                                 !
    !=========================================================================!

    use comms, only: comms_reduce
    use finite_differences, only: FD_GRID_INFO

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)     :: bion_dens(:,:,:,:) ! last index is ion #
    real(kind=DP), intent(in)      :: phi(:,:,:)
    real(kind=DP), intent(in)      :: betamu_elec(1:n_boltzmann_ions)
    real(kind=DP), intent(in)      :: big_gamma
    type(FD_GRID_INFO), intent(in) :: mg_grid
    real(kind=DP), intent(out), optional :: opt_bion_charge_tot(:)

    ! jd: Local variables
    integer :: b_ion
    real(kind=DP) :: z_i
    real(kind=DP) :: my_sum

    ! -------------------------------------------------------------------------

    call implicit_solvent_boltzmann_ion_conc(bion_dens, &
         big_gamma, phi, betamu_elec, mg_grid)

    do b_ion = 1, n_boltzmann_ions
       z_i = boltzmann_ions(b_ion)%charge
       bion_dens(:,:,:,b_ion) = z_i * bion_dens(:,:,:,b_ion)
    end do

    if(present(opt_bion_charge_tot)) then
       do b_ion = 1, n_boltzmann_ions
          my_sum = sum(bion_dens(:,:,:,b_ion))
          call comms_reduce('SUM',my_sum)
          opt_bion_charge_tot(b_ion) = my_sum * mg_grid%weight
       end do
    end if

  end subroutine implicit_solvent_boltzmann_ion_dens
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_boltzmann_debug_dump(phi, mg_return, grid, &
       mg_grid, cell)
    !=========================================================================!
    ! Produces a debug dump of quantities relevant to Boltzmann solvation.    !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2019.                                 !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: ANGSTROM, CONC_AU_TO_MOL_PER_L, k_B
    use finite_differences, only: FD_GRID_INFO
    use rundat, only: pub_is_pbe, pub_is_pbe_temperature
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_int_to_str
    use visual, only: visual_scalarfield

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in)          :: phi(:,:,:)
    type(IS_PBE_MG_RETURN), intent(in) :: mg_return
    type(GRID_INFO), intent(in)        :: grid
    type(FD_GRID_INFO), intent(in)     :: mg_grid
    type(CELL_INFO), intent(in)        :: cell

    ! jd: Local variables
    real(kind=DP), allocatable :: bion_conc(:,:,:,:)
    real(kind=DP) :: bion_conc_avg(1:n_boltzmann_ions)
    real(kind=DP) :: bion_charge_tot(1:n_boltzmann_ions)
    real(kind=DP) :: big_gamma
    real(kind=DP) :: boltzmann_beta
    integer :: b_ion
    integer :: ierr
    character(len=*), parameter :: myself = &
         'implicit_solvent_boltzmann_debug_dump'

    ! -------------------------------------------------------------------------

    allocate(bion_conc(mg_grid%ld1f, mg_grid%ld2f, mg_grid%max_slabs12, &
         n_boltzmann_ions),stat=ierr)
    call utils_alloc_check(myself,'bion_conc',ierr)

    big_gamma = boltzmann_calc_big_gamma(mg_grid)

    ! jd: Overwrite bulk concentrations with current, possibly scaled, concent-
    !     rations from DL_MG. Currently they are only scaled in the counterion
    !     model. In other cases DL_MG returns the same concentrations that we
    !     passed to it, making this statement a no-op.
    !     This is normally done in implicit_solvent_boltzmann_main(), but we
    !     are called before that call.
    boltzmann_ions(1:n_boltzmann_ions)%concentration_adj = &
        mg_return%used_ion_concentration(1:n_boltzmann_ions)

    ! jd: Bion concentrations. Each separately, but also accumulate to print
    !     out total.
    call implicit_solvent_boltzmann_ion_conc(bion_conc, &
         big_gamma, phi, mg_return%betamu_elec(1:n_boltzmann_ions), mg_grid, &
         bion_conc_avg)

    do b_ion = 1, n_boltzmann_ions
       call visual_scalarfield(bion_conc(:,:,:,b_ion), grid, cell, &
            'Bion concentration (mol/L), species '//&
            trim(boltzmann_ions(b_ion)%name), '_out_bion_conc_species_'//&
            trim(utils_int_to_str(b_ion)), &
            conversion_factor = CONC_AU_TO_MOL_PER_L)
       if(b_ion > 1) then
          bion_conc(:,:,:,1) = bion_conc(:,:,:,1) + bion_conc(:,:,:,b_ion)
       end if
    end do
    call visual_scalarfield(bion_conc(:,:,:,1), grid, cell, &
         'Bion concentration (mol/L), total', '_out_bion_conc_total', &
         conversion_factor = CONC_AU_TO_MOL_PER_L)

    ! jd: Bion charge densities. Each separately, but also accumulate to print
    !     out total.
    call implicit_solvent_boltzmann_ion_dens(bion_conc, phi, &
         mg_return%betamu_elec(1:n_boltzmann_ions), big_gamma, mg_grid, &
         bion_charge_tot)
    do b_ion = 1, n_boltzmann_ions
       call visual_scalarfield(bion_conc(:,:,:,b_ion), grid, cell, &
            'Bion charge density (e/ang^3), species '//&
            trim(boltzmann_ions(b_ion)%name), '_out_bion_density_species_'//&
            trim(utils_int_to_str(b_ion)), conversion_factor = ANGSTROM**3)
       if(b_ion > 1) then   ! converts to electrons -ve conv'n ^
          bion_conc(:,:,:,1) = bion_conc(:,:,:,1) + bion_conc(:,:,:,b_ion)
       end if
       ! jd: If LINEARISED, also output z_i*phi(r)/kT so that users can check
       !     how the linearisation condition holds. The absolute value of this
       !     should be << 1.
       if(pub_is_pbe == "LINEARISED") then
          boltzmann_beta = 1.0_DP / (k_B*pub_is_pbe_temperature)
          ! jd: Linearisation condition itself
          call visual_scalarfield(phi(:,:,:) * boltzmann_beta * &
               boltzmann_ions(b_ion)%charge, grid, cell, &
               'Linearisation condition, species '//&
               trim(boltzmann_ions(b_ion)%name), &
               '_out_linearisation_cond_species_'//&
               trim(utils_int_to_str(b_ion)), conversion_factor = -1.0_DP)
                               ! converts to electrons -ve conv'n ^
          ! jd: ... and weighted by accessibility
          call visual_scalarfield(phi(:,:,:) * boltzmann_beta * &
               boltzmann_ions(b_ion)%charge * gamma(:,:,:), grid, cell, &
               'Acc.-weighted lin. cond., species '//&
               trim(boltzmann_ions(b_ion)%name), &
               '_out_lin_cond_gamma_species_'//&
               trim(utils_int_to_str(b_ion)), conversion_factor = -1.0_DP)
                               ! converts to electrons -ve conv'n ^
       end if
    end do
    call visual_scalarfield(bion_conc(:,:,:,1), grid, cell, &
         'Bion charge density (e/ang^3), total', '_out_bion_density_total', &
         conversion_factor = ANGSTROM**3)
                           ! ^ converts to electrons -ve convention
    call boltzmann_ions_print_stats(bion_conc_avg, bion_charge_tot, &
         big_gamma, mg_return, mg_grid)

    deallocate(bion_conc,stat=ierr)
    call utils_dealloc_check(myself,'bion_conc',ierr)

  end subroutine implicit_solvent_boltzmann_debug_dump
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_boltzmann_ion_stats(phi, mg_grid, grid, &
       mg_return, big_gamma)
    !=========================================================================!
    ! Calculates Boltzmann ion concentrations and produces a short report on  !
    ! the root proc.                                                          !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   phi (input): The electrostatic potential from the solution of the PB  !
    !                equation.                                                !
    !   mg_grid (input): The MG grid on which this happens.                   !
    !   grid (input): The ONETEP grid.                                        !
    !   mg_return (in): Details of the solution returned by DL_MG.            !
    !      %betamu_elec (in): Boltzmann \beta * \mu_elec.                     !
    !   big_gamma (in): Accessible fraction (\Gamma), calculated in ONETEP.   !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic in October 2019.                              !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use finite_differences, only: FD_GRID_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in)          :: phi(:,:,:)
    type(FD_GRID_INFO), intent(in)     :: mg_grid
    type(GRID_INFO), intent(in)        :: grid
    type(IS_PBE_MG_RETURN), intent(in) :: mg_return
    real(kind=DP), intent(in)          :: big_gamma

    ! jd: Local variables
    real(kind=DP), allocatable :: bion_conc(:,:,:,:)
    real(kind=DP), allocatable :: bion_dens(:,:,:,:)
    real(kind=DP) :: bion_conc_avg(1:n_boltzmann_ions)
    real(kind=DP) :: bion_charge_tot(1:n_boltzmann_ions)
    integer :: ierr
    character(len=*), parameter :: myself = 'implicit_solvent_boltzmann_ion_stats'

    ! -------------------------------------------------------------------------

    allocate(bion_conc(mg_grid%ld1f, mg_grid%ld2f, mg_grid%max_slabs12, &
         n_boltzmann_ions),stat=ierr)
    call utils_alloc_check(myself,'bion_conc',ierr)
    allocate(bion_dens(mg_grid%ld1f, mg_grid%ld2f, mg_grid%max_slabs12, &
         n_boltzmann_ions),stat=ierr)
    call utils_alloc_check(myself,'bion_dens',ierr)

    call implicit_solvent_boltzmann_ion_conc(bion_conc, &
         big_gamma, phi, mg_return%betamu_elec(1:n_boltzmann_ions), mg_grid, &
         bion_conc_avg)

    call implicit_solvent_boltzmann_ion_dens(bion_dens, phi, &
         mg_return%betamu_elec(1:n_boltzmann_ions), big_gamma, mg_grid, &
         bion_charge_tot)

    call boltzmann_ions_print_stats(bion_conc_avg, bion_charge_tot, &
         big_gamma, mg_return, mg_grid)

    call boltzmann_conc_discrepancy_stats(bion_conc, grid, mg_grid)

    deallocate(bion_dens,stat=ierr)
    call utils_dealloc_check(myself,'bion_dens',ierr)
    deallocate(bion_conc,stat=ierr)
    call utils_dealloc_check(myself,'bion_conc',ierr)

  end subroutine implicit_solvent_boltzmann_ion_stats
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  integer function implicit_solvent_boltzmann_mg_neutralisation_method( &
       mg, grid, rho)
    !=========================================================================!
    ! Returns a selector recognised by DL_MG corresponding to a text ID       !
    ! recognised by ONETEP for is_pbe_neutralisation_scheme.                  !
    ! Also, ensures that if the system is sufficiently close to being charge  !
    ! neutral, counterions are not used and jellium is used instead.          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mg (in): The grid on which MG operates.                               !
    !   grid (in): The ONETEP grid.                                           !
    !   rho (in): The density passed to DL_MG (normally elec + smeared ions). !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in April 2020.                                !
    ! Extended by Jacek Dziedzic in January 2021 to avoid counterions if the  !
    ! system is charge-neutral.                                               !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root
    use constants, only: SAFE_DIV_EPS, CRLF
    use finite_differences, only: FD_GRID_INFO
    use integrals, only: integrals_product_on_grid
    use rundat, only: pub_is_pbe_neutralisation_scheme
    use utils, only: utils_abort, utils_feature_not_supported

    implicit none

    ! jd: Arguments
    type(FD_GRID_INFO), intent(in) :: mg
    type(GRID_INFO), intent(in)    :: grid
    real(kind=DP), intent(in)      :: rho(mg%ld1f,mg%ld2f,mg%max_slabs12)

    ! jd: Local variables
    real(kind=DP) :: tot_quantum_charge
    character(len=*), parameter :: myself = &
         'implicit_solvent_boltzmann_mg_neutralisation_method'

    ! -------------------------------------------------------------------------

#ifdef HAVE_DL_MG
    select case(pub_is_pbe_neutralisation_scheme)
    case('NONE')
       implicit_solvent_boltzmann_mg_neutralisation_method = &
            dl_mg_neutralise_none
    case('JELLIUM')
       implicit_solvent_boltzmann_mg_neutralisation_method = &
            dl_mg_neutralise_with_jellium_uniform
    case('ACCESSIBLE_JELLIUM')
       implicit_solvent_boltzmann_mg_neutralisation_method = &
            dl_mg_neutralise_with_jellium_vacc
    case('COUNTERIONS_FIXED')
       implicit_solvent_boltzmann_mg_neutralisation_method = &
            dl_mg_neutralise_with_ions_fixed
    case('COUNTERIONS_AUTO')
       implicit_solvent_boltzmann_mg_neutralisation_method = &
            dl_mg_neutralise_with_ions_auto
    case('COUNTERIONS_AUTO_LINEAR')
       implicit_solvent_boltzmann_mg_neutralisation_method = &
            dl_mg_neutralise_with_ions_auto_linear
    case default
       call utils_abort(myself//": Unrecognised 'is_pbe_neutralisation_scheme'"&
            //": '"//trim(pub_is_pbe_neutralisation_scheme)//"'.")
    end select

    ! jd: Counterion shifts {x_i} become singular if the total quantum charge
    !     is zero. In practice, for neutral systems the total quantum charge is
    !     very close to zero due to numerical noise. If the noise is sufficient
    !     to avoid the singularity, we keep counterion neutralisation. But if
    !     the noise is dangerously low, we switch to jellium to avoid division
    !     by zero.
    if(pub_is_pbe_neutralisation_scheme(1:11) == 'COUNTERIONS') then
       tot_quantum_charge = integrals_product_on_grid(grid, rho, &
            m1=mg%pq1f, m2=mg%pq2f, m3=mg%num_slabs12_pq)
       if(abs(tot_quantum_charge) < SAFE_DIV_EPS) then
          ! jd: Override to jellium
          implicit_solvent_boltzmann_mg_neutralisation_method = &
               dl_mg_neutralise_with_jellium_uniform
          if(pub_on_root) then
             write(stdout,'(a)') "WARNING: Electrolyte neutralisation scheme &
                  &overridden to 'JELLIUM' because"//CRLF//"         the system&
                  & is charge-neutral. This is because counterion shifts are"//&
                  CRLF//"         singular in charge-neutral systems."
          end if
       end if
    end if

#else
    call utils_feature_not_supported('DL_MG multigrid solver')
     implicit_solvent_boltzmann_mg_neutralisation_method = -1 ! jd: Silence compiler warning
#endif

  end function implicit_solvent_boltzmann_mg_neutralisation_method
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !---------------------------------------------------------------------------!
  !                     P r i v a t e   s u b r o u t i n e s                 !
  !---------------------------------------------------------------------------!

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine boltzmann_conc_discrepancy_stats(bion_conc, grid, mg_grid)
    !=========================================================================!
    ! Calculates and prints statistics on the discrepancy (difference)        !
    ! between the Boltzmann ion concentrations and the user-defined bulk      !
    ! concentratons on the faces of the simulation cell. Only fully electro-  !
    ! lyte-accessible surfaces are considered. A surface is deemed fully      !
    ! accessible if the accessibility is 1.0 on all its points.               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   bion_conc (in): Concentrations of Boltzmann ions on the grid.         !
    !   grid (in): The ONETEP grid.                                           !
    !   mg_grid (in): The DL_MG grid.                                         !
    !-------------------------------------------------------------------------!
    ! Used state:                                                             !
    !   gamma (module-wide): The electrolyte accessibility on the grid.       !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2020.                                !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root, comms_reduce, comms_bcast
    use constants, only: CONC_AU_TO_MOL_PER_L, CRLF, garbage_real
    use finite_differences, only: FD_GRID_INFO
    use rundat, only: pub_print_qc
    use utils, only: utils_qc_print, utils_int_to_str

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in)      :: bion_conc(:,:,:,:)
    type(GRID_INFO), intent(in)    :: grid
    type(FD_GRID_INFO), intent(in) :: mg_grid

    ! jd: Local variables
    logical :: used_faces(6) ! left, right, front, back, bottom, top
    integer :: face_starts(3,6)
    integer :: face_ends(3,6)
    integer :: face
    integer :: n_points_on_face
    integer :: b_ion
    real(kind=DP) :: conc_wanted
    real(kind=DP) :: conc_avg
    real(kind=DP) :: discrep_sum(8)   !} 1..6: six faces, 7: sum over faces for
    real(kind=DP) :: discrep_sumsq(8) !} current bion, 8: sum over faces & bions
    real(kind=DP) :: face_ds_weights(6)
    real(kind=DP) :: total_surf_area
    real(kind=DP) :: gamma_prod
    character(len=4), parameter :: face_names(7) = &
         ['Xmin','Xmax','Ymin','Ymax','Zmin','Zmax',' AVG']
    character(len=*), parameter :: myself = 'boltzmann_conc_discrepancy_stats'

    ! -------------------------------------------------------------------------

    used_faces(1:6) = .true.

    face_starts = reshape([&
         [1,1,1],                       &  ! X==1
         [mg_grid%pq1f,1,1],            &  ! X==max
         [1,1,1],                       &  ! Y==1
         [1,mg_grid%pq2f,1],            &  ! Y==max
         [1,1,1],                       &  ! Z==1
         [1,1,mg_grid%num_slabs12_pq]], &  ! Z==max
         (/3,6/))

    face_ends = reshape([&
         [1,mg_grid%pq2f,mg_grid%num_slabs12_pq],             &  ! X==1
         [mg_grid%pq1f,mg_grid%pq2f,mg_grid%num_slabs12_pq],  &  ! X==max
         [mg_grid%pq1f,1,mg_grid%num_slabs12_pq],             &  ! Y==1
         [mg_grid%pq1f,mg_grid%pq2f,mg_grid%num_slabs12_pq],  &  ! Y==max
         [mg_grid%pq1f,mg_grid%pq2f,1],                       &  ! Z==1
         [mg_grid%pq1f,mg_grid%pq2f,mg_grid%num_slabs12_pq]], &  ! Z==max
         (/3,6/))

    ! jd: Weights for dS -- the grid is orthorhombic, but not necessarily cubic.
    !     These would normally be used for integration, but we are looking at
    !     averages of intensive quantities, so we don't actually ever integrate
    !     over dS -- just sum points. However, to calculate averages over faces
    !     we need to take the surface area into account in the average weights.
    face_ds_weights = [&
         grid%da2%y*grid%da3%z, & ! X==1
         grid%da2%y*grid%da3%z, & ! X==max
         grid%da1%x*grid%da3%z, & ! Y==1
         grid%da1%x*grid%da3%z, & ! Y==max
         grid%da1%x*grid%da2%y, & ! Z==1
         grid%da1%x*grid%da2%y]   ! Z==max

    if(pub_on_root) then
       write(stdout,'(/a)') '+---- Concentration difference from bulk value on&
            & simcell faces (mol/L) ---+'//CRLF//'|   # | Name | Face | Bulk c&
            &onc. | Actual conc. |  Difference |  RMS diff. |'//CRLF//&
            '|------------------------------------------------&
            &--------------------------|'
    end if

    ! jd: Look at the accessibilities over the six faces. If a face is not
    !     fully electrolyte-accessible, exclude it from further consideration.
    do face = 1, 6
       gamma_prod = product(gamma(&
            face_starts(1,face):face_ends(1,face), &
            face_starts(2,face):face_ends(2,face), &
            face_starts(3,face):face_ends(3,face)))

       ! jd: For all faces except Z faces we need to reduce over MPI procs
       !     due to the slab distribution.
       if(face<=4) call comms_reduce('PROD',gamma_prod)

       ! jd: For Z faces we need to get the result from the procs that own
       !     them. Results on remaining procs are meaningless and are
       !     ignored.
       if(face==5) call comms_bcast(grid%proc_slab12(1),gamma_prod)
       if(face==6) call comms_bcast(grid%proc_slab12(mg_grid%pq3f),gamma_prod)

       if(gamma_prod < 1.0_DP) then
          used_faces(face) = .false.
       end if
    end do

    if(any(used_faces(:))) then

       discrep_sum(8) = 0.0_DP
       discrep_sumsq(8) = 0.0_DP
       conc_avg = 0.0_DP

       ! -----------------------------------------------------------------------
       ! jd: Compute and print he difference on the surfaces of the box.
       ! -----------------------------------------------------------------------
       loop_bions:                                                             &
       ! -----------------------------------------------------------------------
       do b_ion = 1, n_boltzmann_ions
          conc_wanted = boltzmann_ions(b_ion)%concentration_user

          discrep_sum(1:6) = garbage_real
          discrep_sumsq(1:6) = garbage_real
          discrep_sum(7) = 0.0_DP
          discrep_sumsq(7) = 0.0_DP
          total_surf_area = 0.0_DP

          ! --------------------------------------------------------------------
          loop_faces:                                                          &
          ! --------------------------------------------------------------------
          do face = 1, 6
             if(used_faces(face)) then
                discrep_sum(face) = sum(bion_conc(&
                     face_starts(1,face):face_ends(1,face), &
                     face_starts(2,face):face_ends(2,face), &
                     face_starts(3,face):face_ends(3,face),b_ion) &
                     - conc_wanted)
                discrep_sumsq(face) = sum((bion_conc(&
                     face_starts(1,face):face_ends(1,face), &
                     face_starts(2,face):face_ends(2,face), &
                     face_starts(3,face):face_ends(3,face),b_ion) &
                     - conc_wanted)**2)

                ! jd: Number of points on face (but mind slab distribution).
                n_points_on_face = &
                     (face_ends(1,face)-face_starts(1,face)+1) * &
                     (face_ends(2,face)-face_starts(2,face)+1) * &
                     (face_ends(3,face)-face_starts(3,face)+1)

                ! jd: For all faces except Z faces we need to reduce over MPI
                !     procs due to the slab distribution.
                if(face<=4) then
                   call comms_reduce('SUM', discrep_sum(face))
                   call comms_reduce('SUM', discrep_sumsq(face))
                   call comms_reduce('SUM', n_points_on_face)
                end if

                ! jd: For Z faces we need to get the result from the procs that
                !     own them. Results on remaining procs are meaningless and
                !     are ignored.
                if(face==5) then
                   call comms_bcast(grid%proc_slab12(1), discrep_sum(5))
                   call comms_bcast(grid%proc_slab12(1), discrep_sumsq(5))
                   call comms_bcast(grid%proc_slab12(1), n_points_on_face)
                end if
                if(face==6) then
                   call comms_bcast(grid%proc_slab12(mg_grid%pq3f), &
                        discrep_sum(6))
                   call comms_bcast(grid%proc_slab12(mg_grid%pq3f), &
                        discrep_sumsq(6))
                   call comms_bcast(grid%proc_slab12(1), n_points_on_face)
                end if

                ! jd: Divide by global number of points to get averages
                discrep_sum(face) = discrep_sum(face) / n_points_on_face
                discrep_sumsq(face) = discrep_sumsq(face) / n_points_on_face

                ! jd: Accumulate average, weighted by the surface area of each
                !     face
                discrep_sum(7) = discrep_sum(7) + &
                     discrep_sum(face) * n_points_on_face * &
                     face_ds_weights(face)
                discrep_sumsq(7) = discrep_sumsq(7) + &
                     discrep_sumsq(face) * n_points_on_face * &
                     face_ds_weights(face)

                total_surf_area = total_surf_area + &
                     n_points_on_face * face_ds_weights(face)

                if(pub_on_root) then
                   write(stdout, &
                        '(a,i3,a,a4,a,a4,a,f10.7,a,f12.7,a,f11.7,a,f10.7,a)') &
                        '| ', b_ion, ' | ', trim(boltzmann_ions(b_ion)%name), &
                        ' | ', face_names(face), ' | ', &
                        conc_wanted * CONC_AU_TO_MOL_PER_L, ' | ', &
                        (conc_wanted+discrep_sum(face)) * CONC_AU_TO_MOL_PER_L,&
                        ' | ', discrep_sum(face) * CONC_AU_TO_MOL_PER_L, ' | ',&
                        sqrt(discrep_sumsq(face)) * CONC_AU_TO_MOL_PER_L, ' | '
                   if(pub_print_qc) then
                      call utils_qc_print('conc_d_face_'//&
                           trim(utils_int_to_str(face))//'_bion_'//&
                           trim(utils_int_to_str(b_ion)),discrep_sum(face))
                   end if
                end if

             end if

          end do loop_faces

          ! jd: Summary for average over faces.
          if(pub_on_root) then

             ! jd: Divide by total surface area to get the average from the
             !     weighted sum.
             discrep_sum(7) = discrep_sum(7) / total_surf_area
             discrep_sumsq(7) = discrep_sumsq(7) / total_surf_area

             write(stdout,'(a)') '+.........................................&
                  &.................................+'
             write(stdout,'(a,i3,a,a4,a,a4,a,f10.7,a,f12.7,a,f11.7,a,f10.7,a)')&
                  '| ', b_ion, ' | ', trim(boltzmann_ions(b_ion)%name), ' | ', &
                  face_names(7), ' | ', &
                  conc_wanted * CONC_AU_TO_MOL_PER_L, ' | ', &
                  (conc_wanted + discrep_sum(7)) * CONC_AU_TO_MOL_PER_L, &
                  ' | ', discrep_sum(7) * CONC_AU_TO_MOL_PER_L, ' | ', &
                  sqrt(discrep_sumsq(7)) * CONC_AU_TO_MOL_PER_L, ' | '
             write(stdout,'(a)') '+-----------------------------------------&
                  &---------------------------------+'
          end if

          ! jd: Sums over bions
          discrep_sum(8) = discrep_sum(8) + discrep_sum(7)
          discrep_sumsq(8) = discrep_sumsq(8) + discrep_sumsq(7)
          conc_avg = conc_avg + conc_wanted

       end do loop_bions

       discrep_sum(8) = discrep_sum(8) / real(n_boltzmann_ions,kind=DP)
       discrep_sumsq(8) = discrep_sumsq(8) / real(n_boltzmann_ions,kind=DP)
       conc_avg = conc_avg / real(n_boltzmann_ions,kind=DP)

       ! jd: Summary for average over Boltzmann ions.
       if(pub_on_root) then
          write(stdout,'(a)') '+-----------------------------------------------&
               &---------------------------+'
          write(stdout,'(a,f10.7,a,f12.7,a,f11.7,a,f10.7,a)')&
               '|      OVERALL      | ', &
               conc_avg * CONC_AU_TO_MOL_PER_L, ' | ', &
               (conc_avg + discrep_sum(8)) * CONC_AU_TO_MOL_PER_L, &
               ' | ', discrep_sum(8) * CONC_AU_TO_MOL_PER_L, ' | ', &
               sqrt(discrep_sumsq(8)) * CONC_AU_TO_MOL_PER_L, ' | '
          write(stdout,'(a)') '+-----------------------------------------&
               &---------------------------------+'
       end if
    else
       if(pub_on_root) then
          write(stdout,'(a)') '+         There are no fully electrolyte-accessi&
               &ble simcell faces.         +'
          write(stdout,'(a)') '+-----------------------------------------------&
               &---------------------------+'
       end if
    end if

  end subroutine boltzmann_conc_discrepancy_stats
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_boltzmann_chem_pot(mu_ideal, mu_excess, &
         mu_total, betamu_elec, beta, big_gamma)
    !=========================================================================!
    ! Calculates the components to the chemical potential and the total       !
    ! chemical potential -- relevant in Boltzmann solvation.                  !
    !=========================================================================!

    use constants, only: CONC_0, SAFE_DIV_EPS
    use rundat, only: pub_multigrid_bc_is_periodic
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: mu_ideal(1:n_boltzmann_ions)
    real(kind=DP), intent(out) :: mu_excess(1:n_boltzmann_ions)
    real(kind=DP), intent(out) :: mu_total(1:n_boltzmann_ions)
    real(kind=DP), intent(in)  :: betamu_elec(1:n_boltzmann_ions)
    real(kind=DP), intent(in)  :: beta
    real(kind=DP), intent(in)  :: big_gamma

    ! jd: Local variables
    integer :: b_ion
    logical :: pbc
    real(kind=DP) :: c_i_bulk
    character(len=*), parameter :: myself = 'implicit_solvent_boltzmann_chem_pot'

    ! -------------------------------------------------------------------------

    if(any(pub_multigrid_bc_is_periodic)) then
       call utils_assert(all(pub_multigrid_bc_is_periodic), &
            myself//': Mixed boundary conditions are not supported.')
       pbc = .true.
    else
       pbc = .false.
    end if

    do b_ion = 1, n_boltzmann_ions
       c_i_bulk = boltzmann_ions(b_ion)%concentration_adj
       if(c_i_bulk < SAFE_DIV_EPS) then
          ! jd: Corner case of zero or near-zero concentration. In that case
          !     the number of particles is zero, and \mu_ideal is also zero.
          mu_ideal(b_ion) = 0.0_DP
       else
          mu_ideal(b_ion) = log(c_i_bulk/CONC_0) / beta
       end if
       if(pbc) then
          mu_excess(b_ion) = betamu_elec(b_ion) / beta
          mu_total(b_ion) = mu_ideal(b_ion) + mu_excess(b_ion)
       else
          mu_excess(b_ion) = 0.0_DP
          mu_total(b_ion) = mu_ideal(b_ion)
       end if
    end do

  end subroutine implicit_solvent_boltzmann_chem_pot
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine boltzmann_init_steric_pot(mdl)
    !=========================================================================!
    ! Generates the steric potential, stores it into the module-wide variable !
    ! 'steric_pot', then calculates the accessibility (gamma), which in this  !
    ! model is e^(-beta*steric_pot). It is also stored in a module-wide var.  !
    !                                                                         !
    ! Currently the steric potential can be:                                  !
    ! - soft-core: the repulsive part of the Lennard-Jones potential, i.e.    !
    !   A/r^12 centred around each ion, truncated and shifted to zero at the  !
    !   truncation radius.                                                    !
    ! - hard-core: \inf within a radius around each atom, 0 elsewhere.        !
    !              In practice 1 is stored in steric_pot, not \inf, to avoid  !
    !              IEEE734 traps.                                             !
    ! - smeared hard-core: Practically \inf within a radius around each atom, !
    !                      0 elsewehere, with a smooth transition.            !
    !                                                                         !
    ! The steric potential and weight are generated on the mg grid.           !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   mdl (in/out): The model in use. Writable, because we potentially      !
    !                 update its hc_steric_cutoffs.                           !
    ! ----------------------------------------------------------------------- !
    ! Caveats:                                                                !
    !  - 'steric_pot' and 'gamma' are allocated here. Deallocation happens in !
    !    happens in implicit_solvent_boltzmann_exit().                        !
    !  - a copy of the mg grid is constructed from 'grid', as the original is !
    !    unavailable (cannot use multigrid_methods here). Cannot pass 'mg'    !
    !    either, because 'grid' is needed in the visual_scalar call.          !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 07/2014.                                     !
    ! Modified to remove pub_par, Joseph Prentice September 2018              !
    ! Extended for smeared hc steric potential by Jacek Dziedzic, July 2019.  !
    !=========================================================================!

    use comms, only: pub_on_root, comms_bcast, pub_root_proc_id
!$  use comms, only: pub_my_proc_id
    use constants, only: stdout, VERBOSE, PI, k_B
    use finite_differences, only: FD_GRID_INFO, finite_difference_set_geometry
    use geometry, only: POINT, magnitude, OPERATOR(-), OPERATOR(+), OPERATOR(*)
    use ion, only: RADIAL_DENSITY_TYPE
    use model_type, only: MODEL
    use rundat, only: pub_num_spins, pub_is_pbe, pub_is_sc_steric_magnitude, &
         pub_is_steric_write, pub_is_sc_steric_cutoff, &
         pub_is_steric_pot_type, pub_is_hc_steric_smearing, &
         pub_is_sc_steric_smoothing_alpha, pub_output_detail, &
         pub_multigrid_bc_is_periodic, pub_mg_granularity_power, &
         pub_is_pbe_temperature, pub_initial_dens_realspace, &
         pub_is_hc_steric_dens_isovalue
!$  use rundat, only: pub_threads_max
    use simulation_cell, only: minimum_image_distance_ortho
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_abort, utils_assert, utils_erf, &
         utils_trace_in, utils_trace_out, utils_real_to_str, utils_int_to_str
    use visual, only: visual_scalarfield

    implicit none

    ! jd: Arguments
    type(MODEL), intent(inout), target :: mdl

    ! jd: Local variables
    type(FD_GRID_INFO) :: mg_grid
    type(RADIAL_DENSITY_TYPE), pointer :: rad_dens
    integer       :: islab12, i1, i2, i3, ipt, I
    integer       :: ierr, nat
    integer       :: multigrid_granularity
    integer       :: isp, jsp
    integer       :: isub
    integer       :: iradpt
    type(POINT)   :: r, R_I
    real(kind=DP) :: boltzmann_beta
    real(kind=DP) :: r_here
    real(kind=DP) :: dens_here, max_dens
    real(kind=DP) :: gamma_here
    real(kind=DP) :: d, accum
    real(kind=DP) :: cutoff
    real(kind=DP) :: shift ! jd: Shift of the potential to get zero at cutoff
    real(kind=DP) :: val_at_zero ! jd: Limit of the potential at d->0
    logical       :: pbc
    character(len=*), parameter :: myself = 'boltzmann_init_steric_pot'

    !------------------------------------------------------------------------

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    if(pub_is_pbe == 'NONE') then
       call utils_abort(trim(myself)//' called with ''is_pbe none''.')
    end if

    ! jcap: find number of atoms
    nat=size(mdl%elements)

    ! jd: Initialise a copy of the MG grid. Ugly, I know.
    multigrid_granularity = 2**pub_mg_granularity_power
    mg_grid = finite_difference_set_geometry(mdl%fine_grid, &
         multigrid_granularity, is_periodic = pub_multigrid_bc_is_periodic)

    ! ab: read the type of steric potential
    steric_pot_type = pub_is_steric_pot_type

    select case(steric_pot_type)
    case('X')
       if(pub_on_root) write(stdout,'(/a)') &
            'IS: No steric potential will be applied!'
    case('H')
       if(pub_on_root) write(stdout,'(/a)') &
            'IS: Generating hard-core steric potential ...'
    case('M')
       if(pub_on_root) write(stdout,'(/a)') &
            'IS: Generating smoothed hard-core steric potential ...'
    case('S')
       if(pub_on_root) write(stdout,'(/a)') &
            'IS: Generating soft-core steric potential ...'
    case default
       call utils_abort(trim(myself)//': Unrecognized steric potential')
    end select

    allocate(steric_pot(mg_grid%ld1f,mg_grid%ld2f,mg_grid%max_slabs12), &
         stat=ierr)
    call utils_alloc_check(myself,'steric_pot',ierr)
    allocate(gamma(mg_grid%ld1f,mg_grid%ld2f,mg_grid%max_slabs12), &
         stat=ierr)
    call utils_alloc_check(myself,'gamma',ierr)

    ! jd: Take care of padding between pt and ld
    steric_pot = 0.0_DP
    gamma = 1.0_DP

    ! jd: If no steric potential desired, we just exit here with 0 everywhere,
    !     but mind to output it, if desired, and stop the timer
    if(steric_pot_type == 'X') goto 999

    if(steric_pot_type == 'S') then
       shift = pub_is_sc_steric_cutoff**(-12.0_DP) * &
            utils_erf(pub_is_sc_steric_smoothing_alpha * &
            pub_is_sc_steric_cutoff)**12.0_DP
       val_at_zero = 4096.0_DP * pub_is_sc_steric_magnitude * &
            pub_is_sc_steric_smoothing_alpha**12.0_DP / PI**6.0_DP

       if(pub_on_root .and. pub_output_detail >= VERBOSE) then
          write(stdout,'(a,e12.5,a)') &
               'IS: - maximum value (at core centres):        ', val_at_zero,' Ha.'
          write(stdout,'(a,e12.5,a)') &
               'IS: - constant shift (to get zero at cutoff): ', shift,' Ha.'
       end if
    end if


    if(steric_pot_type == 'M') then
       if(pub_on_root .and. pub_output_detail >= VERBOSE) then
          write(stdout,'(a,e12.5,a)') &
               'IS: - smearing:       ', pub_is_hc_steric_smearing,' bohr.'
       end if
    end if

    ! ab: For both hard-core steric potentials, the cutoff has two components:
    !     one specified by species_solvent_radius and another calculated from
    !     isovalue of electronic density which is established from the
    !     atomsolver radial densities. The component from radial density
    !     is not calculated for a negative value of is_hc_steric_dens_isovalue

    if((steric_pot_type == 'H' .or. steric_pot_type == 'M') .and. &
         pub_is_hc_steric_dens_isovalue >= 0.0_DP) then
       call utils_assert(pub_initial_dens_realspace, myself//&
            ": When implicit solvent electrolyte accessibility is determined &
            &through a hard-core potential, with automatically determined radii&
            &('is_hc_steric_dens_isovalue' is present and non-negative), 'initi&
            &al_dens_realspace' must be T, because the radii for all species ar&
            &e determined from atomsolver radial densities")

       if(pub_on_root) then
          ! jd: This code courtesy of Joe Prentice scans region-wide radial
          !     densities for the species we need.
          ! jcap: loop over all species
          do isp=1,size(mdl%species)
             ! jcap: loop over regions
             do isub=1,mdl%nsub
                ! jcap: loop over species within this region
                do jsp=1,size(mdl%regions(isub)%species)
                   ! jcap: Check if species match
                   if(mdl%species(isp)%species_id == &
                        mdl%regions(isub)%species(jsp)%species_id) then

                      call utils_assert(mdl%regions(isub)%&
                           radial_densities(jsp)%present, myself//&
                           ': Radial density for species '//&
                           trim(utils_int_to_str(jsp))//' in region '//&
                           trim(utils_int_to_str(isub))//' absent on root proc.')

                      ! jcap: If they match, use this instance of
                      !       regions%radial_densities
                      rad_dens => mdl%regions(isub)%radial_densities(jsp)

                      r_here = -1D0
                      max_dens = -1D0
                      do iradpt = rad_dens%npts, 1, -1
                         if(pub_num_spins == 1) then
                            dens_here = rad_dens%den(iradpt,1)
                         else
                            dens_here = 0.5_DP * (rad_dens%den(iradpt,1) + &
                                 rad_dens%den(iradpt,2))
                         end if
                         max_dens = max(max_dens,dens_here)
                         if(dens_here >= pub_is_hc_steric_dens_isovalue) then
                            r_here = rad_dens%rad(iradpt)
                            mdl%species(jsp)%hc_steric_cutoff = r_here
                            exit
                         end if
                      end do
                      call utils_assert(r_here /= -1D0, myself//&
                           ': Could not determine a suitable steric potential &
                           &cut-off radius from the radial density. Desired &
                           &density isovalue: '//trim(utils_real_to_str(&
                           pub_is_hc_steric_dens_isovalue))//'. Number of &
                           &radial points: '//trim(utils_int_to_str(&
                           rad_dens%npts))//'. Maximum found density: '//&
                           trim(utils_real_to_str(max_dens))//'.')
                   end if ! regions match
                end do
             end do
          end do

       end if ! on root

       do isp=1,size(mdl%species)
          call comms_bcast(pub_root_proc_id, mdl%species(isp)%hc_steric_cutoff)
       end do

    end if ! steric pot is hc

    !ab: print radial cutoff for hard core steric potential
    if (steric_pot_type == 'H' .or. steric_pot_type == 'M') then
       if (pub_on_root) then
          write(stdout,'(a)') &
               '---------------------------------------------------------------'
          write(stdout,'(a)') &
               '|      Atomic hard-core steric potential radial cut-offs      |'
          write(stdout,'(a)') &
               '+-------------------------------------------------------------+'
          write(stdout,'(a)') &
               '|   Species |   Symbol |            Radial cut-off (bohr)     |'
          write(stdout,'(a)') &
               '|           |          |     solute      solvent        total |'
          do isp = 1, size(mdl%species)
             write(stdout,'(a,i9,a,a4,a,f10.6,a,f10.6,a,f10.6,a)') &
                  '| ', isp, ' |     ', adjustr(mdl%species(isp)%symbol), &
                  ' | ', mdl%species(isp)%hc_steric_cutoff, ' + ', &
                  mdl%species(isp)%solvent_radius,' = ', &
                  mdl%species(isp)%solvent_radius + &
                  mdl%species(isp)%hc_steric_cutoff, ' |'
          end do
          write(stdout,'(a)') &
               '---------------------------------------------------------------'
       end if ! on root
    end if !steric pot is hc

    if(any(pub_multigrid_bc_is_periodic)) then
       call utils_assert(all(pub_multigrid_bc_is_periodic), &
            myself//': Mixed boundary conditions are not supported.')
       pbc = .true.
    else
       pbc = .false.
    end if

    ! jd: Generate steric potential on every point of fine grid
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,i3,r,accum,I,R_I,d,gamma_here,isp,cutoff) &
!$OMP SHARED(nat,mdl,mg_grid,pub_is_sc_steric_cutoff, pbc, &
!$OMP      pub_is_pbe_temperature, pub_is_hc_steric_smearing, &
!$OMP      pub_is_sc_steric_magnitude, &
!$OMP      steric_pot,pub_my_proc_id, val_at_zero, steric_pot_type, &
!$OMP      pub_is_sc_steric_smoothing_alpha, shift, pub_threads_max)
    do ipt=1,mg_grid%num_slabs12_pq*mg_grid%pq2f*mg_grid%pq1f
       i1 = modulo(ipt-1,mg_grid%pq1f) + 1
       i2 = modulo((ipt-i1)/mg_grid%pq1f,mg_grid%pq2f) + 1
       islab12 = (ipt-(i2-1)*mg_grid%pq1f-i1) / (mg_grid%pq1f*mg_grid%pq2f) + 1
       i3 = mg_grid%my_first_slab12 + islab12 - 1

       r%X = real((i1-1),kind=DP) * mg_grid%d1f
       r%Y = real((i2-1),kind=DP) * mg_grid%d2f
       r%Z = real((i3-1),kind=DP) * mg_grid%d3f

       accum = 0.0_DP
       select case(steric_pot_type)
       case('S')
          ! jd: Soft-core potential
          do I = 1, nat

             R_I = mdl%elements(I)%centre

             if(pbc) then
                d = minimum_image_distance_ortho(r,R_I,mdl%cell)
             else
                d = magnitude(r-R_I)
             end if

             ! jd: Only calculate steric potential within a cutoff around each
             !     atom, to retain O(N) scaling
             if(d <= pub_is_sc_steric_cutoff) then
                if(d < 1D-10) then ! jd: Very close to gridpoint, avoid singularity
                   accum = val_at_zero / pub_is_sc_steric_magnitude
                else ! jd: Usual case
                   accum = accum + d**(-12.0_DP) * &
                        utils_erf(pub_is_sc_steric_smoothing_alpha*d)**12.0_DP
                end if
             end if

          end do ! over I

          ! jd: Multiply by prefactor (A), apply constant shift (from each atom)
          !     Values that are negative (rounding errors in erf) are zeroed.
          steric_pot(i1,i2,islab12) = &
               max(pub_is_sc_steric_magnitude * accum - nat * shift,0.0_DP)

        case('H')
          ! jd: Hard-core potential
          do I = 1, nat

             R_I = mdl%elements(I)%centre
             if(pbc) then
                d = minimum_image_distance_ortho(r,R_I,mdl%cell)
             else
                d = magnitude(r-R_I)
             end if

             isp = mdl%elements(I)%global_species_number
             cutoff = mdl%elements(I)%solvent_radius + &
                  mdl%species(isp)%hc_steric_cutoff

             if(d <= cutoff) then
                accum = 1.0_DP
                exit ! jd: Enough if within one cutoff
             end if

          end do ! over I

          steric_pot(i1,i2,islab12) = accum

        case('M')
          ! jd: Smoothed hard-core potential
          do I = 1, nat

             R_I = mdl%elements(I)%centre
             if(pbc) then
                d = minimum_image_distance_ortho(r,R_I,mdl%cell)
             else
                d = magnitude(r-R_I)
             end if

             isp = mdl%elements(I)%global_species_number
             cutoff = mdl%elements(I)%solvent_radius + &
                  mdl%species(isp)%hc_steric_cutoff

             ! jd: Points very far from ions do not contribute. Elide them to
             !     improve scaling.
             if(d > 2D0 * cutoff) then
                cycle
             end if

             ! jd: Usual points -- use smoothed accessibility
             gamma_here = 0.5_DP + utils_erf(&
                  (d-cutoff)/pub_is_hc_steric_smearing) * 0.5_DP
             if(gamma_here < 1D-7) gamma_here = 1D-90
             gamma_here = min(1.0_DP, gamma_here)   ! jd: Ensure we stay within
             gamma_here = max(1D-90, gamma_here)    !     [0+tiny,1]

             ! jd: Accessibility, gamma, is 1 - 1/(1+exp((d-d0)/s))
             !     Convert to *potential* (=-kT ln gamma) so that if we are
             !     within multiple spheres, we accumulate correctly
             accum = accum -k_B * pub_is_pbe_temperature * log(gamma_here)

          end do ! over I

          steric_pot(i1,i2,islab12) = accum

        end select

    end do ! ipt
!$OMP END PARALLEL DO

    ! jd: Transform steric potential to accessibility (gamma), favoured
    !     by DL_MG (called 'steric_weight' there).
    boltzmann_beta = 1.0_DP / (k_B*pub_is_pbe_temperature)
    select case(steric_pot_type)
    case('S')
       gamma = exp(-boltzmann_beta*steric_pot)
    case('M')
       gamma = exp(-boltzmann_beta*steric_pot)
    case('H')
       ! jd: Already generated as zeros and ones.
       !     Translate pot=0 -> weight of 1
       !               pot=1 -> weight of 0
       gamma = 1.0_DP - steric_pot
    case('X')
       ! jd: Already one everywhere. Just a sanity check.
       call utils_assert(maxval(steric_pot)==1.0_DP .and. &
            minval(steric_pot)==1.0_DP,trim(myself)//': Zero steric &
            &potential (1 steric weight) sanity check fail.', &
            minval(gamma),maxval(gamma))
    end select

    if(pub_on_root) then
       write(stdout,'(a)') ''
    end if

999 if(pub_is_steric_write) then
       call visual_scalarfield(steric_pot, mdl%fine_grid, mdl%cell, &
            'Steric potential for PB equation, atomic units', '_steric_pot')
       call visual_scalarfield(gamma, mdl%fine_grid, mdl%cell, &
            'Steric weight for PB equation', '_steric_weight')
    end if

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine boltzmann_init_steric_pot
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  real(kind=DP) function boltzmann_calc_pot_average(mg_grid, phi, &
       big_gamma)
    !=========================================================================!
    ! Calculates the average of the electrostatic potential, taking into      !
    ! account the accessibility. Currently this is never used.                !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   mg_grid (in): The MG grid over which we calculate.                    !
    !   phi (in): The electrostatic potential on the grid.                    !
    !   big_gamma (in): Accessible fraction \Gamma, as calculated by ONETEP.  !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic in November 2019.                             !
    !=========================================================================!

    use comms, only: comms_reduce
    use finite_differences, only: FD_GRID_INFO

    implicit none

    ! jd: Arguments
    type(FD_GRID_INFO), intent(in) :: mg_grid
    real(kind=DP), intent(in)      :: phi(:,:,:)
    real(kind=DP), intent(in)      :: big_gamma

    ! jd: Local variables
    integer       :: i1, i2, islab12
    integer       :: tot_n_points
    real(kind=DP) :: box_volume
    real(kind=DP) :: phi_here
    real(kind=DP) :: gamma_here
    real(kind=DP) :: pot_sum
    real(kind=DP) :: pot_int
    character(len=*), parameter :: myself = 'boltzmann_calc_pot_average'

    ! -------------------------------------------------------------------------

    tot_n_points = mg_grid%pq1f * mg_grid%pq2f * mg_grid%pq3f

    box_volume = real(tot_n_points,kind=DP) * mg_grid%weight

    pot_sum = 0.0_DP
    do islab12 = 1, mg_grid%num_slabs12_pq
       do i2 = 1, mg_grid%pq2f
          do i1 = 1, mg_grid%pq1f
             phi_here = phi(i1,i2,islab12)
             gamma_here = gamma(i1,i2,islab12)
             pot_sum = pot_sum + gamma_here * phi_here
          end do
       end do
    end do

    call comms_reduce('SUM',pot_sum)

    pot_int = pot_sum * mg_grid%weight

    boltzmann_calc_pot_average = 1.0_DP / (box_volume * big_gamma) * pot_int

  end function boltzmann_calc_pot_average

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  real(kind=DP) function boltzmann_calc_big_gamma(mg_grid, steric_weight_average)
    !=========================================================================!
    ! Calculates ONETEP's equivalent of steric_weight_average (fraction of    !
    ! accessible volume), Gamma = \int dr gamma(r)/V.                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mg_grid (in): The grid on which this happens.                         !
    !   steric_weight_average(in): DL_MG's idea of the same quantity, to      !
    !                              check if they match (in PBCs).             !
    ! Return value:                                                           !
    !   ONETEP's version of \Gamma.                                           !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in December 2019.                             !
    !=========================================================================!

    use comms, only: comms_reduce
    use finite_differences, only: FD_GRID_INFO
    use rundat, only: pub_multigrid_bc_is_periodic
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(FD_GRID_INFO), intent(in)      :: mg_grid
    real(kind=DP), intent(in), optional :: steric_weight_average

    ! jd: Local variables
    integer :: tot_n_points
    logical :: pbc
    real(kind=DP) :: gamma_sum
    real(kind=DP) :: big_gamma
    character(len=*), parameter :: myself = 'boltzmann_calc_big_gamma'

    ! -------------------------------------------------------------------------

    if(any(pub_multigrid_bc_is_periodic)) then
       call utils_assert(all(pub_multigrid_bc_is_periodic), &
            myself//': Mixed boundary conditions are not supported.')
       pbc = .true.
    else
       pbc = .false.
    end if

    tot_n_points = mg_grid%pq1f * mg_grid%pq2f * mg_grid%pq3f

    gamma_sum = &
         sum(gamma(1:mg_grid%pq1f,1:mg_grid%pq2f,1:mg_grid%num_slabs12_pq))

    call comms_reduce('SUM',gamma_sum)

    big_gamma = gamma_sum / real(tot_n_points,kind=DP)

    ! jd: Trap scenario where we excluded practically everything
    call utils_assert(big_gamma > 1D-6, myself//&
         ': Unexpectedly low solvent-accessible fraction of volume (\Gamma): ',&
         big_gamma)

    ! jd: Trap scenario where we somehow have negative exclusions
    call utils_assert(big_gamma < 1.000001_DP, myself//&
         ': Solvent-accessible fraction of volume (\Gamma) in excess of 1.0: ',&
         big_gamma)

    ! jd: Trap scenario where our calculated \Gamma does not match DL_MG's
    !     In OBCs we do *NOT* do this comparison, as DL_MG calculates \Gamma
    !     excluding the first and last points of the grid. This is because in
    !     OBCs DL_MG does not treat them as part of the domain, only as a
    !     boundary condition. Becoming consistent with this interpretation would
    !     be complicated, because it would involve changing the integration
    !     domain in ONETEP not only for \Gamma, but also for energies and it
    !     would redefine volume. The change in energy terms could potentially
    !     affect non-Boltzmann OBC solvation results (in cases where there is
    !     non-negligible charge density on the boundary). We agreed to leave
    !     this be and simply ignore the different \Gamma reported by DL_MG in
    !     OBCs. We do not worry about this "discrepancy" since \Gamma does not
    !     enter any energy expressions in OBCs.
    if(present(steric_weight_average)) then
       if(pbc) then
          call utils_assert(abs(steric_weight_average-big_gamma) < 0.00001_DP , &
               myself//': Solvent-accessible fraction of volume (\Gamma) does not &
               &match DL_MG''s idea: ', big_gamma, steric_weight_average)
       end if
    end if

    boltzmann_calc_big_gamma = big_gamma

  end function boltzmann_calc_big_gamma
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine boltzmann_ions_report
    !=========================================================================!
    ! Writes out a report on what Boltzmann ions are present in the system.   !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   None.                                                                 !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic in 07/2014.                                   !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, CRLF, CONC_MOL_PER_L_TO_AU
    use rundat, only: pub_is_pbe
    use utils, only: utils_assert

    implicit none

    ! jd: Local variables
    integer :: i

    ! -------------------------------------------------------------------------

    call utils_assert(pub_is_pbe /= 'NONE', &
         'boltzmann_ions_report() must only be called when pub_is_pbe is not &
         &''none''.')

    if(boltzmann_ions_reported) return

    if(pub_on_root) then
       write(stdout,'(a)') &
            '----------------------------------------------------------'&
            //CRLF//&
            '|           Boltzmann ions in implicit solvent:          |'&
            //CRLF//&
            '+--------------------------------------------------------+'&
            //CRLF//&
            '|       Name | Charge |  Conc. (mol/L)  |  Conc. (a.u.)  |'
       do i=1, n_boltzmann_ions
          write(stdout,'(a2,a10,a2,f7.3,a6,e11.5,a3,e14.5,a3)') '| ', &
               adjustr(boltzmann_ions(i)%name), ' |',boltzmann_ions(i)%charge,&
               ' |    ',&
               boltzmann_ions(i)%concentration_user / CONC_MOL_PER_L_TO_AU,'  |', &
               boltzmann_ions(i)%concentration_user,'  |'
       end do
       write(stdout,'(a)') '---------------------------------------------------&
            &-------'
    end if

    boltzmann_ions_reported = .true.

  end subroutine boltzmann_ions_report
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine boltzmann_ions_print_stats(bion_conc_avg, bion_charge_tot, &
       big_gamma, mg_return, mg_grid)
    !=========================================================================!
    ! Writes out a report on Boltzmann ion concentrations on root proc.       !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   bion_conc_avg(in): Average concentrations for all Boltzmann ions,     !
    !                      obtained from implicit_solvent_boltzman_ion_conc().!
    !   bion_charge_tot(in): Total charges for all Boltzmann ions, obtained   !
    !                        from implicit_solvent_boltzman_ion_dens().       !
    !   big gamma (in): Accessible volume fraction from ONETEP.               !
    !
    !   mg_grid (input): The grid on which this happens.                      !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic in October 2019.                              !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: CONC_AU_TO_MOL_PER_L, k_B, CRLF, garbage_real
    use finite_differences, only: FD_GRID_INFO
    use rundat, only: pub_multigrid_bc_is_periodic, pub_is_pbe_temperature, &
         pub_is_pbe_neutralisation_scheme
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in)          :: bion_conc_avg(1:n_boltzmann_ions)
    real(kind=DP), intent(in)          :: bion_charge_tot(1:n_boltzmann_ions)
    real(kind=DP), intent(in)          :: big_gamma
    type(IS_PBE_MG_RETURN), intent(in) :: mg_return
    type(FD_GRID_INFO), intent(in)     :: mg_grid

    ! jd: Local variables
    real(kind=DP) :: mu_ideal(1:n_boltzmann_ions)
    real(kind=DP) :: mu_excess(1:n_boltzmann_ions)
    real(kind=DP) :: mu_total(1:n_boltzmann_ions)
    real(kind=DP) :: box_volume
    real(kind=DP) :: beta
    character(len=3) :: bc_string
    logical :: pbc
    integer :: b_ion
    character(len=*), parameter :: myself = 'boltzmann_ions_print_stats'

    ! -------------------------------------------------------------------------

    if(pub_on_root) then

       if(any(pub_multigrid_bc_is_periodic)) then
          call utils_assert(all(pub_multigrid_bc_is_periodic), &
               myself//': Mixed boundary conditions are not supported.')
          bc_string = 'PBC'
          pbc = .true.
       else
          bc_string = 'OBC'
          pbc = .false.
       end if

       beta = 1.0_DP/(k_B*pub_is_pbe_temperature)

       call implicit_solvent_boltzmann_chem_pot(mu_ideal, mu_excess, &
            mu_total, mg_return%betamu_elec(1:n_boltzmann_ions), beta, &
            big_gamma)

       box_volume = mg_grid%pq1f * mg_grid%pq2f * mg_grid%pq3f * mg_grid%weight

       write(stdout,'(/a,f14.3,a,f8.6,a,a,a)') 'IS: MG box volume: ', &
            box_volume, ' a0^3. Accessible fraction: ', big_gamma, &
            '. ', bc_string,'.'
       if(pbc) then
          write(stdout,'(a)') 'IS: Charge neutralisation scheme: '//&
               trim(pub_is_pbe_neutralisation_scheme)//'.'
       end if
       write(stdout,'(/a)') '+--------------------------- Chemical potential --&
            &------------------------+'//CRLF//'|   # | Name | Bulk conc. |   &
            & mu_ideal   |   mu_excess   |    mu_total   |'//CRLF//'|---------&
            &----------------------------------------------------------------|'
       do b_ion = 1, n_boltzmann_ions
          write(stdout,'(a,i3,a,a4,a,f10.7,a,f13.9,a,f13.9,a,f13.9,a)')&
               '| ', b_ion, ' | ', trim(boltzmann_ions(b_ion)%name), ' | ', &
               boltzmann_ions(b_ion)%concentration_adj * CONC_AU_TO_MOL_PER_L, &
               ' | ', mu_ideal(b_ion), ' | ', mu_excess(b_ion), ' | ', &
               mu_total(b_ion), ' |'
       end do

       write(stdout,'(a)') '|-------------------------------------------------&
            &------------------------|'

       write(stdout,'(/a)') '+----------------- Boltzmann ion concentration and&
            & charge -----------------+'//CRLF//'|   # | Name | Bulk &
            &conc. | Average conc. | Neutr. ratio |    Total charge |'//CRLF//&
            '|---------------------------------------------------------&
            &-----------------|'
       do b_ion = 1, n_boltzmann_ions
          ! jd: Show information on neutralisation ratios only if counterions
          !     were actually used. When neutr. ratios are 'garbage_real', it
          !     means we actually used jellium after having determined the
          !     system is too close to charge-neutral for using counterions.
          if(pub_is_pbe_neutralisation_scheme(1:11) == 'COUNTERIONS' .and. &
               .not. all(mg_return%used_neutralisation_ratios(1:n_boltzmann_ions) &
               == garbage_real)) then
             if(abs(mg_return%used_neutralisation_ratios(b_ion)) < 1D4) then
                write(stdout,'(a,i3,a,a4,a,f10.7,a,f13.9,a,f12.6,a,f15.9,a)')&
                     '| ', b_ion, ' | ', trim(boltzmann_ions(b_ion)%name), ' | ', &
                     boltzmann_ions(b_ion)%concentration_adj * CONC_AU_TO_MOL_PER_L, &
                     ' | ', bion_conc_avg(b_ion)/big_gamma * CONC_AU_TO_MOL_PER_L, &
                     ' | ', mg_return%used_neutralisation_ratios(b_ion), ' | ', &
                     -bion_charge_tot(b_ion), ' | '
             else
                write(stdout,'(a,i3,a,a4,a,f10.7,a,f13.9,a,e12.6,a,f15.9,a)')&
                     '| ', b_ion, ' | ', trim(boltzmann_ions(b_ion)%name), ' | ', &
                     boltzmann_ions(b_ion)%concentration_adj * CONC_AU_TO_MOL_PER_L, &
                     ' | ', bion_conc_avg(b_ion)/big_gamma * CONC_AU_TO_MOL_PER_L, &
                     ' | ', mg_return%used_neutralisation_ratios(b_ion), ' | ', &
                     -bion_charge_tot(b_ion), ' | '

             end if
          else
             write(stdout,'(a,i3,a,a4,a,f10.7,a,f13.9,a,a,a,f15.9,a)')&
                  '| ', b_ion, ' | ', trim(boltzmann_ions(b_ion)%name), ' | ', &
                  boltzmann_ions(b_ion)%concentration_adj * CONC_AU_TO_MOL_PER_L, &
                  ' | ', bion_conc_avg(b_ion)/big_gamma * CONC_AU_TO_MOL_PER_L, &
                  ' | ', '------------', ' | ', &
                  -bion_charge_tot(b_ion), ' | '
          end if
       end do
       write(stdout,'(a)') '+-------------------------------------------------&
            &-------------------------+'
       write(stdout,'(a,f15.9,a)') &
            '                                                    total: ', &
            sum(-bion_charge_tot(:)), merge('*',' ',pbc)

    end if

  end subroutine boltzmann_ions_print_stats

end module is_solvation_boltzmann
